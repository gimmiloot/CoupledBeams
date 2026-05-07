from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Sequence

import numpy as np


# ============================================================
# USER PARAMETERS
# Change these values for ordinary runs.
# CLI arguments, if provided, override these defaults.
# ============================================================

DEFAULT_BRANCH_ID = "bending_desc_05"
DEFAULT_BRANCH_NUMBER = 5
DEFAULT_BETA = 15.0
DEFAULT_MU = 0.2
DEFAULT_EPSILON = 0.0025
DEFAULT_L_TOTAL = 2.0
DEFAULT_NUM_ANALYTIC_ROOTS = 20
DEFAULT_OUTPUT_PREFIX = None
DEFAULT_SAVE_DEBUG_CSV = True
DEFAULT_SCAN_THETA_CONVENTIONS = False
DEFAULT_LOCALIZE_RESIDUAL = False
DEFAULT_SCAN_BENDING_BASIS_PAIRINGS = False
DEFAULT_SCAN_AXIAL_CONVENTIONS = False
DEFAULT_SCAN_AXIAL_SCALES = False
DEFAULT_SCAN_RIGHT_TRANSFORM_CONVENTIONS = False
DEFAULT_AUDIT_FREQUENCY_SCALING = False


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import (  # noqa: E402
    BeamParams,
    assemble_clamped_coupled_matrix,
)
from my_project.fem import python_fem as fem  # noqa: E402
from scripts.analyze_flat_mu_bending_energy import (  # noqa: E402
    N_SOLVE as FEM_N_SOLVE,
    N_TRACK as FEM_N_TRACK,
    assign_by_mac_and_frequency,
    build_branch_metadata,
)
from scripts.compare_beta0_analytic_vs_fem import (  # noqa: E402
    BASE_PARAMS,
    fem_axial_fractions,
    fem_parameter_override,
)
from scripts.lib.analytic_branch_tracking import (  # noqa: E402
    DEFAULT_BETA_STEPS,
    DEFAULT_MAC_WARNING_THRESHOLD,
    DEFAULT_MU_STEPS,
    DEFAULT_N_TRACK,
    base_sorted_index_from_branch_id,
    branch_id_from_base_sorted_index,
    dense_mu_values_for_targets,
    track_mu_sweep,
)
from scripts.lib.analytic_coupled_rods_shapes import (  # noqa: E402
    NEAR_ZERO_NORM,
    ROW_LABELS,
    analytic_null_vector,
    reconstruct_analytic_components,
)
from scripts.lib.tracked_bending_descendant_shapes import compute_local_components  # noqa: E402
from scripts.sweep_grid_policy import analysis_beta_grid, analysis_mu_grid  # noqa: E402


RESULTS_DIR = REPO_ROOT / "results"
RIGHT_THETA_SIGNS = (1.0, -1.0)
RESIDUAL_GOOD_THRESHOLD = 1e-4
RESIDUAL_BAD_THRESHOLD = 1e-3
FREQUENCY_GOOD_THRESHOLD = 1e-2
GLOBAL_MAC_GOOD_THRESHOLD = 0.98
GLOBAL_MAC_BAD_THRESHOLD = 0.95
GLOBAL_REL_BAD_THRESHOLD = 0.25
DEFAULT_MODAL_PROJECTION_MODES = 30

SUMMARY_FIELDNAMES = [
    "branch_id",
    "beta",
    "mu",
    "epsilon",
    "base_sorted_index",
    "current_sorted_index",
    "analytic_lambda",
    "fem_lambda",
    "rel_lambda_diff",
    "branch_tracking_mac_min",
    "branch_tracking_warning_flag",
    "fem_sorted_index",
    "residual_l2",
    "relative_residual",
    "residual_inf",
    "relative_residual_inf",
    "global_dof_relative_l2",
    "global_dof_mac_euclidean",
    "global_dof_mac_mass_weighted",
    "translational_relative_l2",
    "rotational_relative_l2",
    "left_arm_relative_l2",
    "right_arm_relative_l2",
    "joint_relative_l2",
    "conclusion_flag",
    "l_total",
    "radius",
    "analytic_omega_nd",
    "analytic_eigenvalue_omega_sq",
    "fem_omega_nd",
    "fem_eigenvalue_omega_sq",
    "fem_selected_self_residual_l2",
    "fem_selected_self_relative_residual",
    "fem_frequency_hz",
    "right_theta_sign_used",
    "right_theta_sign_best_by_fem_residual",
    "residual_l2_theta_pos_fem_omega",
    "relative_residual_theta_pos_fem_omega",
    "residual_inf_theta_pos_fem_omega",
    "relative_residual_inf_theta_pos_fem_omega",
    "residual_l2_theta_neg_fem_omega",
    "relative_residual_theta_neg_fem_omega",
    "residual_inf_theta_neg_fem_omega",
    "relative_residual_inf_theta_neg_fem_omega",
    "residual_l2_theta_pos_analytic_omega",
    "relative_residual_theta_pos_analytic_omega",
    "residual_l2_theta_neg_analytic_omega",
    "relative_residual_theta_neg_analytic_omega",
    "selected_analytic_omega_relative_residual",
    "nearest_fem_sorted_index_by_lambda",
    "nearest_fem_lambda_by_lambda",
    "nearest_fem_rel_lambda_diff_by_lambda",
    "analytic_fem_sorted_index_delta",
    "branch_tracking_mac_target",
    "branch_tracking_mac_mean",
    "branch_tracking_path_points",
    "smallest_singular_value",
    "singular_value_ratio",
    "analytic_matrix_residual_inf",
    "analytic_embedding_local_relative_l2",
    "analytic_embedding_local_max_abs",
    "analytic_vs_fem_local_relative_l2",
    "analytic_vs_fem_local_mac_euclidean",
    "rel_l2_theta_vs_dwds_left",
    "rel_l2_theta_vs_minus_dwds_left",
    "rel_l2_theta_vs_dwds_right",
    "rel_l2_theta_vs_minus_dwds_right",
    "rel_l2_theta_vs_dwdxi_left",
    "rel_l2_theta_vs_minus_dwdxi_left",
    "rel_l2_theta_vs_dwdxi_right",
    "rel_l2_theta_vs_minus_dwdxi_right",
    "rel_l2_theta_vs_dwdz_left",
    "rel_l2_theta_vs_minus_dwdz_left",
    "rel_l2_theta_vs_dwdz_right",
    "rel_l2_theta_vs_minus_dwdz_right",
    "best_theta_convention_left",
    "best_theta_convention_right",
    "best_theta_scale_left",
    "best_theta_scale_right",
    "theta_convention_scan_csv",
    "theta_convention_scan_best",
    "theta_convention_scan_best_relative_residual",
    "analytic_rayleigh_omega_sq",
    "fem_selected_omega_sq",
    "relative_rayleigh_omega_sq_diff",
    "analytic_rayleigh_Lambda",
    "relative_rayleigh_Lambda_diff",
    "top_modal_projection_mode",
    "top_modal_projection_share",
    "cumulative_projection_share_first_3",
    "residual_largest_group",
    "residual_largest_group_share",
    "residual_largest_element_arm",
    "residual_largest_element_index",
    "residual_largest_element_share",
    "q_joint_translation_mismatch",
    "q_joint_rotation_mismatch",
    "q_external_clamp_max_abs",
    "modal_projection_csv",
    "residual_groups_csv",
    "element_residual_energy_csv",
    "bending_basis_pairing_scan_csv",
    "full_basis_column_audit_csv",
    "bending_basis_canonical_relative_residual",
    "bending_basis_best_variant",
    "bending_basis_best_relative_residual",
    "bending_basis_best_preserves_matrix_consistency",
    "bending_basis_pairing_mismatch_flag",
    "axial_convention_scan_csv",
    "axial_convention_baseline_relative_residual",
    "axial_convention_baseline_rayleigh_lambda",
    "axial_convention_best_variant",
    "axial_convention_best_relative_residual",
    "axial_convention_best_rayleigh_lambda",
    "axial_convention_mismatch_flag",
    "axial_scale_scan_csv",
    "axial_scale_baseline_relative_residual",
    "axial_scale_baseline_rayleigh_lambda",
    "axial_scale_best_residual_variant",
    "axial_scale_best_residual_alpha_left",
    "axial_scale_best_residual_alpha_right",
    "axial_scale_best_residual_relative_residual",
    "axial_scale_best_residual_rayleigh_lambda",
    "axial_scale_best_rayleigh_variant",
    "axial_scale_best_rayleigh_alpha_left",
    "axial_scale_best_rayleigh_alpha_right",
    "axial_scale_best_rayleigh_relative_residual",
    "axial_scale_best_rayleigh_lambda",
    "axial_scale_best_mass_mac_variant",
    "axial_scale_best_mass_mac",
    "axial_scale_mismatch_flag",
    "right_transform_convention_scan_csv",
    "right_transform_sanity_csv",
    "right_transform_baseline_relative_residual",
    "right_transform_baseline_rayleigh_lambda",
    "right_transform_best_residual_variant",
    "right_transform_best_residual_relative_residual",
    "right_transform_best_residual_rayleigh_lambda",
    "right_transform_best_rayleigh_variant",
    "right_transform_best_rayleigh_relative_residual",
    "right_transform_best_rayleigh_lambda",
    "right_transform_best_mass_mac_variant",
    "right_transform_best_mass_mac",
    "right_transform_mismatch_flag",
    "fem_transform_inference",
    "frequency_scaling_audit_csv",
    "frequency_scaling_best_analytic_factor",
    "frequency_scaling_best_analytic_relative_residual",
    "frequency_scaling_best_fem_factor",
    "frequency_scaling_best_fem_relative_residual",
    "joint_translation_embedding_mismatch",
    "joint_theta_embedding_mismatch",
    "global_alignment_scalar",
    "summary_csv",
    "debug_nodes_csv",
]

NODE_FIELDNAMES = [
    "node_index",
    "arm",
    "s_norm",
    "x0",
    "y0",
    "analytic_ux",
    "analytic_uy",
    "analytic_theta",
    "fem_ux",
    "fem_uy",
    "fem_theta",
    "diff_ux",
    "diff_uy",
    "diff_theta",
]

THETA_CONVENTION_SCAN_FIELDNAMES = [
    "theta_convention",
    "theta_left_sign",
    "theta_right_sign",
    "theta_left_scale_description",
    "theta_right_scale_description",
    "relative_residual",
    "global_dof_mac_euclidean",
    "mass_weighted_mac",
    "translational_relative_l2",
    "rotational_relative_l2",
    "rotational_mac",
    "left_arm_relative_l2",
    "right_arm_relative_l2",
    "translational_mac",
    "residual_l2",
    "residual_inf",
    "relative_residual_inf",
    "translational_residual_norm",
    "rotational_residual_norm",
    "left_arm_residual_norm",
    "right_arm_residual_norm",
    "joint_residual_norm",
    "global_dof_relative_l2",
    "joint_relative_l2",
    "left_theta_convention",
    "right_theta_convention",
    "global_alignment_scalar",
]

MODAL_PROJECTION_FIELDNAMES = [
    "rank_by_mass_share",
    "mode_index",
    "omega_nd",
    "lambda_value",
    "coefficient_euclidean",
    "coefficient_abs_euclidean",
    "euclidean_projection_share",
    "coefficient_mass_weighted",
    "coefficient_abs_mass_weighted",
    "coefficient_abs",
    "coefficient_share",
    "mass_weighted_projection_share",
    "cumulative_mass_weighted_share_by_mode_order",
    "cumulative_mass_weighted_share_by_rank",
    "modal_mass",
    "modal_euclidean_norm_sq",
]

RESIDUAL_GROUP_FIELDNAMES = [
    "group_name",
    "residual_scope",
    "dof_count",
    "group_residual_norm",
    "group_residual_share",
    "group_relative_residual",
    "group_residual_inf",
    "group_kq_norm",
    "group_omega_mq_norm",
]

ELEMENT_RESIDUAL_ENERGY_FIELDNAMES = [
    "arm",
    "element_index",
    "global_node0",
    "global_node1",
    "midpoint_s_norm",
    "analytic_local_u1",
    "analytic_local_v1",
    "analytic_local_theta1",
    "analytic_local_u2",
    "analytic_local_v2",
    "analytic_local_theta2",
    "fem_local_u1",
    "fem_local_v1",
    "fem_local_theta1",
    "fem_local_u2",
    "fem_local_v2",
    "fem_local_theta2",
    "analytic_local_axial_energy",
    "analytic_local_bending_energy",
    "fem_local_axial_energy",
    "fem_local_bending_energy",
    "local_residual_norm",
    "local_residual_norm_share",
    "local_residual_norm_sq_share",
    "local_residual_inf",
    "local_residual_axial_norm",
    "local_residual_bending_norm",
]

BENDING_BASIS_PAIRING_SCAN_FIELDNAMES = [
    "variant_name",
    "left_mapping",
    "right_mapping",
    "right_z_sign",
    "external_clamp_residual",
    "matrix_residual_norm",
    "independent_endpoint_residual_norm",
    "matrix_vs_independent_endpoint_mismatch",
    "analytic_in_fem_relative_residual",
    "global_dof_mac_euclidean",
    "mass_weighted_mac",
    "translational_relative_l2",
    "rotational_relative_l2",
    "conclusion_flag",
]

FULL_BASIS_COLUMN_AUDIT_FIELDNAMES = [
    "basis_column",
    "row_label",
    "matrix_column_value",
    "independent_full_basis_value",
    "abs_difference",
    "relative_difference",
    "variant_name",
    "right_z_sign",
]

BASIS_COLUMN_LABELS = ("A1", "B1", "A2", "B2", "P1", "P2")


@dataclass
class BendingBasisMapping:
    name: str
    swap: bool
    sign_A: float
    sign_B: float


AXIAL_CONVENTION_SCAN_FIELDNAMES = [
    "variant_name",
    "left_u_sign",
    "right_u_mode",
    "relative_residual",
    "rayleigh_omega_sq",
    "rayleigh_lambda",
    "relative_rayleigh_lambda_diff",
    "global_dof_mac_euclidean",
    "mass_weighted_mac",
    "translational_relative_l2",
    "rotational_relative_l2",
    "right_arm_translational_residual_share",
    "axial_to_bending_residual_ratio",
    "analytic_right_axial_energy",
    "fem_right_axial_energy",
    "right_axial_energy_ratio",
    "analytic_total_axial_fraction",
    "fem_total_axial_fraction",
    "left_axial_energy",
    "right_axial_share",
    "rel_l2_du_ds_left",
    "rel_l2_du_ds_right",
    "max_abs_du_ds_left",
    "max_abs_du_ds_right",
    "best_du_ds_right_convention",
    "best_du_ds_right_rel_l2",
    "conclusion_flag",
]

AXIAL_SCALE_SCAN_FIELDNAMES = [
    "variant_name",
    "alpha_left",
    "alpha_right",
    "scale_source",
    "relative_residual",
    "rayleigh_omega_sq",
    "rayleigh_lambda",
    "rayleigh_lambda_relative_diff",
    "global_dof_mac_euclidean",
    "mass_weighted_mac",
    "translational_relative_l2",
    "rotational_relative_l2",
    "right_arm_translational_residual_share",
    "axial_to_bending_residual_ratio",
    "analytic_total_axial_fraction",
    "analytic_right_axial_share",
    "right_axial_energy_ratio_to_fem",
    "conclusion_flag",
]

RIGHT_TRANSFORM_CONVENTION_SCAN_FIELDNAMES = [
    "variant_name",
    "transform_description",
    "relative_residual",
    "rayleigh_lambda",
    "rayleigh_lambda_relative_diff",
    "global_dof_mac_euclidean",
    "mass_weighted_mac",
    "translational_relative_l2",
    "rotational_relative_l2",
    "right_arm_translational_residual_share",
    "axial_to_bending_residual_ratio",
    "joint_translation_mismatch",
    "joint_rotation_mismatch",
    "external_clamp_max_abs",
    "conclusion_flag",
]

RIGHT_TRANSFORM_SANITY_FIELDNAMES = [
    "variant_name",
    "transform_description",
    "fem_assembly_inference",
    "local_unit_axial_global_ux",
    "local_unit_axial_global_uy",
    "local_unit_transverse_global_ux",
    "local_unit_transverse_global_uy",
]

FREQUENCY_SCALING_AUDIT_FIELDNAMES = [
    "vector_type",
    "spectral_factor_name",
    "spectral_factor_value",
    "relative_residual",
    "residual_l2",
    "residual_inf",
    "analytic_lambda",
    "fem_selected_omega",
    "fem_selected_omega_sq",
    "fem_selected_lambda_from_omega",
    "fem_selected_lambda_from_omega_sq",
    "analytic_lambda_squared",
    "analytic_lambda_fourth",
    "rayleigh_omega_sq",
    "rayleigh_lambda",
    "lambda_fourth_to_fem_omega_sq_ratio",
    "lambda_squared_to_fem_omega_ratio",
    "fem_omega_minus_lambda_squared",
    "fem_omega_sq_minus_lambda_fourth",
    "note",
]


@dataclass
class FemBranchState:
    branch_id: str
    current_sorted_index: int
    omega_nd: float
    lambda_value: float
    frequency_hz: float
    vec_free: np.ndarray
    full_vec: np.ndarray
    tracking_mac_min: float
    tracking_mac_mean: float
    tracking_path_points: int


def option_was_provided(argv: Sequence[str], option_name: str) -> bool:
    return any(arg == option_name or arg.startswith(f"{option_name}=") for arg in argv)


def filename_number_token(value: float) -> str:
    return f"{float(value):g}".replace("-", "m").replace("+", "").replace(".", "p")


def default_output_prefix(branch_id: str, beta_deg: float, mu_value: float, epsilon: float) -> Path:
    return RESULTS_DIR / (
        f"analytic_shape_in_fem_residual_beta{filename_number_token(beta_deg)}_{branch_id}"
        f"_mu{filename_number_token(mu_value)}_eps{filename_number_token(epsilon)}"
    )


def resolve_output_prefix(value: str | None, branch_id: str, beta_deg: float, mu_value: float, epsilon: float) -> Path:
    if value is None:
        return default_output_prefix(branch_id=branch_id, beta_deg=beta_deg, mu_value=mu_value, epsilon=epsilon)
    path = Path(value)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def suffixed_path(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}{suffix}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    raw_args = list(sys.argv[1:] if argv is None else argv)
    branch_id_explicit = option_was_provided(raw_args, "--branch-id")
    branch_number_explicit = option_was_provided(raw_args, "--branch-number")

    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Embed a canonically tracked analytic coupled-rods shape into the baseline FEM DOF vector "
            "and compute the direct FEM residual Kq - omega^2 Mq."
        ),
    )
    parser.add_argument("--branch-id", default=DEFAULT_BRANCH_ID)
    parser.add_argument("--branch-number", type=int, default=DEFAULT_BRANCH_NUMBER)
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA)
    parser.add_argument("--mu", type=float, default=DEFAULT_MU)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--l-total", type=float, default=DEFAULT_L_TOTAL)
    parser.add_argument("--num-analytic-roots", type=int, default=DEFAULT_NUM_ANALYTIC_ROOTS)
    parser.add_argument("--output-prefix", default=DEFAULT_OUTPUT_PREFIX)
    parser.add_argument(
        "--scan-theta-conventions",
        action="store_true",
        default=DEFAULT_SCAN_THETA_CONVENTIONS,
        help="Write a diagnostic scan of analytic rotational DOF theta conventions.",
    )
    parser.add_argument(
        "--localize-residual",
        action="store_true",
        default=DEFAULT_LOCALIZE_RESIDUAL,
        help="Write Rayleigh, modal projection, residual-group, and element localization diagnostics.",
    )
    parser.add_argument(
        "--scan-bending-basis-pairings",
        action="store_true",
        default=DEFAULT_SCAN_BENDING_BASIS_PAIRINGS,
        help="Audit full cos/sin/cosh/sinh bending-basis coefficient pairings against determinant columns.",
    )
    parser.add_argument(
        "--scan-axial-conventions",
        action="store_true",
        default=DEFAULT_SCAN_AXIAL_CONVENTIONS,
        help="Scan analytic axial-component sign/reversal conventions in the FEM embedding.",
    )
    parser.add_argument(
        "--scan-axial-scales",
        action="store_true",
        default=DEFAULT_SCAN_AXIAL_SCALES,
        help="Scan analytic axial-component amplitude scale factors in the FEM embedding.",
    )
    parser.add_argument(
        "--scan-right-transform-conventions",
        action="store_true",
        default=DEFAULT_SCAN_RIGHT_TRANSFORM_CONVENTIONS,
        help="Scan right-arm local/global transform conventions for analytic FEM embedding.",
    )
    parser.add_argument(
        "--audit-frequency-scaling",
        action="store_true",
        default=DEFAULT_AUDIT_FREQUENCY_SCALING,
        help="Audit Lambda/omega/omega^2 spectral-factor choices in the FEM residual.",
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--save-debug-csv", dest="save_debug_csv", action="store_true")
    group.add_argument("--no-save-debug-csv", dest="save_debug_csv", action="store_false")
    parser.set_defaults(save_debug_csv=DEFAULT_SAVE_DEBUG_CSV)
    args = parser.parse_args(raw_args)

    if branch_id_explicit and branch_number_explicit:
        try:
            branch_number_from_id = base_sorted_index_from_branch_id(args.branch_id)
        except ValueError as exc:
            parser.error(str(exc))
        if int(branch_number_from_id) != int(args.branch_number):
            parser.error("--branch-id suffix and --branch-number must describe the same base branch.")
    elif branch_number_explicit and not branch_id_explicit:
        args.branch_id = branch_id_from_base_sorted_index(int(args.branch_number))
    elif branch_id_explicit and not branch_number_explicit:
        try:
            args.branch_number = base_sorted_index_from_branch_id(args.branch_id)
        except ValueError as exc:
            parser.error(str(exc))
    else:
        try:
            args.branch_number = base_sorted_index_from_branch_id(args.branch_id)
        except ValueError as exc:
            parser.error(str(exc))

    if args.branch_number <= 0:
        parser.error("--branch-number must be positive.")
    if args.beta < 0.0:
        parser.error("--beta must be non-negative; canonical continuation starts at beta = 0.")
    if args.mu < -1e-12:
        parser.error("--mu must be non-negative for this diagnostic path.")
    if args.epsilon <= 0.0:
        parser.error("--epsilon must be positive.")
    if args.l_total <= 0.0:
        parser.error("--l-total must be positive.")
    if args.num_analytic_roots < args.branch_number:
        parser.error("--num-analytic-roots must be at least the requested branch number.")
    return args


def radius_from_epsilon(epsilon: float, l_total: float) -> float:
    return float(l_total) * float(epsilon)


def build_params_for_case(epsilon: float, l_total: float) -> BeamParams:
    return BeamParams(
        E=float(BASE_PARAMS.E),
        rho=float(BASE_PARAMS.rho),
        r=radius_from_epsilon(epsilon, l_total),
        L_total=float(l_total),
    )


def free_dof_indices() -> np.ndarray:
    n = fem.N_ELEM
    ndof = 3 * (2 * n + 1)
    bc = {0, 1, 2, 3 * (2 * n), 3 * (2 * n) + 1, 3 * (2 * n) + 2}
    return np.array(sorted(set(range(ndof)) - bc), dtype=int)


def expand_full_vector(vec_free: np.ndarray) -> np.ndarray:
    full = np.zeros(3 * (2 * fem.N_ELEM + 1), dtype=float)
    full[free_dof_indices()] = np.asarray(vec_free, dtype=float)
    return full


def assemble_fem_matrices(
    params: BeamParams,
    *,
    mu: float,
    beta_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with fem_parameter_override(params):
        n = fem.N_ELEM
        ndof = 3 * (2 * n + 1)
        length_left = max(float(fem.ell) * (1.0 - float(mu)), 1e-3)
        length_right = max(float(fem.ell) * (1.0 + float(mu)), 1e-3)
        le_left = length_left / n
        le_right = length_right / n
        t_right = fem.rotation_matrix_6x6(np.deg2rad(float(beta_deg)))
        k_left = fem.elem_K(le_left)
        m_left = fem.elem_M(le_left)
        k_right = t_right.T @ fem.elem_K(le_right) @ t_right
        m_right = t_right.T @ fem.elem_M(le_right) @ t_right
        k_full = np.zeros((ndof, ndof), dtype=float)
        m_full = np.zeros((ndof, ndof), dtype=float)

        def assemble(dof_map: list[int], ke: np.ndarray, me: np.ndarray) -> None:
            for row, dof_row in enumerate(dof_map):
                for col, dof_col in enumerate(dof_map):
                    k_full[dof_row, dof_col] += ke[row, col]
                    m_full[dof_row, dof_col] += me[row, col]

        for elem in range(n):
            dofs = [
                3 * elem,
                3 * elem + 1,
                3 * elem + 2,
                3 * (elem + 1),
                3 * (elem + 1) + 1,
                3 * (elem + 1) + 2,
            ]
            assemble(dofs, k_left, m_left)

        for elem in range(n):
            node0 = n + elem
            node1 = n + elem + 1
            dofs = [
                3 * node0,
                3 * node0 + 1,
                3 * node0 + 2,
                3 * node1,
                3 * node1 + 1,
                3 * node1 + 2,
            ]
            assemble(dofs, k_right, m_right)

    free = free_dof_indices()
    return (
        k_full[np.ix_(free, free)],
        m_full[np.ix_(free, free)],
        k_full,
        m_full,
        free,
    )


def solve_fem_modes_nd(
    params: BeamParams,
    *,
    mu: float,
    beta_deg: float,
    n_modes: int,
) -> tuple[np.ndarray, np.ndarray, float]:
    with fem_parameter_override(params):
        omega_nd, vecs = fem.fem_solve(mu=float(mu), beta_deg=float(beta_deg), n_modes=int(n_modes))
        return np.asarray(omega_nd, dtype=float), np.asarray(vecs, dtype=float), float(fem.scale)


def track_fem_branch(
    params: BeamParams,
    *,
    branch_id: str,
    beta_deg: float,
    mu: float,
) -> FemBranchState:
    seed_omega, seed_vecs, scale_hz = solve_fem_modes_nd(
        params,
        mu=0.0,
        beta_deg=0.0,
        n_modes=FEM_N_SOLVE,
    )
    seed_freqs = seed_omega * scale_hz
    seed_fracs = fem_axial_fractions(seed_vecs)
    metadata = build_branch_metadata(seed_freqs=seed_freqs[:FEM_N_TRACK], seed_fracs=seed_fracs[:FEM_N_TRACK])
    branch_lookup = {str(meta["branch_id"]): idx for idx, meta in enumerate(metadata)}
    if branch_id not in branch_lookup:
        available = ", ".join(str(meta["branch_id"]) for meta in metadata)
        raise KeyError(f"FEM branch metadata does not contain {branch_id!r}. Available: {available}")
    branch_pos = int(branch_lookup[branch_id])

    prev_vecs = seed_vecs[:, :FEM_N_TRACK].copy()
    prev_omega = seed_omega[:FEM_N_TRACK].copy()
    prev_freqs = seed_freqs[:FEM_N_TRACK].copy()
    prev_indices = np.arange(1, FEM_N_TRACK + 1, dtype=int)
    mac_history: list[float] = []

    beta_values = analysis_beta_grid(start=0.0, stop=float(beta_deg), anchor_points=(float(beta_deg),))
    for beta_value in beta_values[1:]:
        cur_omega, cur_vecs, _ = solve_fem_modes_nd(
            params,
            mu=0.0,
            beta_deg=float(beta_value),
            n_modes=FEM_N_SOLVE,
        )
        cur_freqs = cur_omega * scale_hz
        assignment, assigned_mac = assign_by_mac_and_frequency(
            prev_vecs=prev_vecs,
            prev_freqs=prev_freqs,
            cur_vecs=cur_vecs,
            cur_freqs=cur_freqs,
        )
        prev_vecs = cur_vecs[:, assignment]
        prev_omega = cur_omega[assignment]
        prev_freqs = cur_freqs[assignment]
        prev_indices = assignment + 1
        mac_history.append(float(assigned_mac[branch_pos]))

    mu_values = analysis_mu_grid(start=0.0, stop=float(mu), anchor_points=(float(mu),))
    for mu_value in mu_values[1:]:
        cur_omega, cur_vecs, _ = solve_fem_modes_nd(
            params,
            mu=float(mu_value),
            beta_deg=float(beta_deg),
            n_modes=FEM_N_SOLVE,
        )
        cur_freqs = cur_omega * scale_hz
        assignment, assigned_mac = assign_by_mac_and_frequency(
            prev_vecs=prev_vecs,
            prev_freqs=prev_freqs,
            cur_vecs=cur_vecs,
            cur_freqs=cur_freqs,
        )
        prev_vecs = cur_vecs[:, assignment]
        prev_omega = cur_omega[assignment]
        prev_freqs = cur_freqs[assignment]
        prev_indices = assignment + 1
        mac_history.append(float(assigned_mac[branch_pos]))

    omega_nd = float(prev_omega[branch_pos])
    finite_mac = [value for value in mac_history if np.isfinite(value)]
    return FemBranchState(
        branch_id=str(branch_id),
        current_sorted_index=int(prev_indices[branch_pos]),
        omega_nd=omega_nd,
        lambda_value=float(np.sqrt(max(omega_nd, 0.0))),
        frequency_hz=float(prev_freqs[branch_pos]),
        vec_free=prev_vecs[:, branch_pos].copy(),
        full_vec=expand_full_vector(prev_vecs[:, branch_pos]),
        tracking_mac_min=float(min(finite_mac)) if finite_mac else np.nan,
        tracking_mac_mean=float(np.mean(finite_mac)) if finite_mac else np.nan,
        tracking_path_points=int(len(mac_history) + 1),
    )


def branch_points_up_to_target(points: Sequence[object], *, beta: float, target_mu: float) -> list[object]:
    out: list[object] = []
    for point in points:
        step_type = str(getattr(point, "step_type"))
        point_beta = float(getattr(point, "beta"))
        point_mu = float(getattr(point, "mu"))
        if step_type in {"base", "beta"} and point_beta <= float(beta) + 1e-10:
            out.append(point)
        elif abs(point_beta - float(beta)) <= 1e-10 and point_mu <= float(target_mu) + 1e-10:
            out.append(point)
    return out


def resolve_analytic_branch(
    *,
    branch_id: str,
    beta_deg: float,
    mu: float,
    epsilon: float,
    num_analytic_roots: int,
):
    base_index = base_sorted_index_from_branch_id(branch_id)
    n_track = max(DEFAULT_N_TRACK, int(base_index))
    n_solve = max(int(num_analytic_roots), n_track)
    mu_values = dense_mu_values_for_targets([float(mu)], mu_steps=DEFAULT_MU_STEPS)
    result = track_mu_sweep(
        epsilon=float(epsilon),
        beta=float(beta_deg),
        mu_values=mu_values,
        n_track=n_track,
        n_solve=n_solve,
        mac_warning_threshold=DEFAULT_MAC_WARNING_THRESHOLD,
        shape_metric="full",
        beta_steps=DEFAULT_BETA_STEPS,
        allow_low_mac=False,
        required_branch_ids=[branch_id],
    )
    target_point = result.point_at(branch_id, beta=float(beta_deg), mu=float(mu))
    branch_points = branch_points_up_to_target(
        result.points_for_branch(branch_id),
        beta=float(beta_deg),
        target_mu=float(mu),
    )
    mac_values = [
        float(getattr(point, "mac_to_previous"))
        for point in branch_points
        if str(getattr(point, "step_type")) != "base" and np.isfinite(float(getattr(point, "mac_to_previous")))
    ]
    warning_flag = "needs_review" if any(str(getattr(point, "warning_flag")) != "ok" for point in branch_points) else "ok"
    return {
        "point": target_point,
        "branch_points": branch_points,
        "mac_min": float(min(mac_values)) if mac_values else np.nan,
        "mac_mean": float(np.mean(mac_values)) if mac_values else np.nan,
        "warning_flag": warning_flag,
    }


def slope_values(a_value: float, b_value: float, z_values: np.ndarray) -> np.ndarray:
    z = np.asarray(z_values, dtype=float)
    return float(a_value) * (-np.sin(z) - np.sinh(z)) + float(b_value) * (np.cos(z) - np.cosh(z))


def analytic_local_fields(
    Lambda: float,
    *,
    mu: float,
    epsilon: float,
    coeff: np.ndarray,
    s_norm: np.ndarray,
    right_theta_sign: float,
) -> dict[str, np.ndarray]:
    components = reconstruct_analytic_components(
        float(Lambda),
        mu_value=float(mu),
        epsilon=float(epsilon),
        coeff=np.asarray(coeff, dtype=float),
        s_norm=np.asarray(s_norm, dtype=float),
        right_coordinate="joint-to-external",
    )
    a1, b1, a2, b2, _, _ = [float(value) for value in coeff]
    xi = np.asarray(s_norm, dtype=float)
    z_left = float(Lambda) * (1.0 - float(mu)) * xi
    xi_external_to_joint = 1.0 - xi
    z_right = -float(Lambda) * (1.0 + float(mu)) * xi_external_to_joint
    theta_left = float(Lambda) * slope_values(a1, b1, z_left)
    theta_right = float(right_theta_sign) * float(Lambda) * slope_values(a2, b2, z_right)
    return {
        "u_left": np.asarray(components["u_left"], dtype=float),
        "w_left": np.asarray(components["w_left"], dtype=float),
        "theta_left": theta_left,
        "u_right": np.asarray(components["u_right"], dtype=float),
        "w_right": np.asarray(components["w_right"], dtype=float),
        "theta_right": theta_right,
    }


def build_analytic_global_q(
    fields: dict[str, np.ndarray],
    *,
    beta_deg: float,
) -> tuple[np.ndarray, dict[str, float]]:
    n = fem.N_ELEM
    n_nodes = 2 * n + 1
    beta_rad = float(np.deg2rad(beta_deg))
    cb = float(np.cos(beta_rad))
    sb = float(np.sin(beta_rad))
    ux = np.zeros(n_nodes, dtype=float)
    uy = np.zeros(n_nodes, dtype=float)
    theta = np.zeros(n_nodes, dtype=float)

    ux[: n + 1] = np.asarray(fields["u_left"], dtype=float)
    uy[: n + 1] = np.asarray(fields["w_left"], dtype=float)
    theta[: n + 1] = np.asarray(fields["theta_left"], dtype=float)

    u_right = np.asarray(fields["u_right"], dtype=float)
    w_right = np.asarray(fields["w_right"], dtype=float)
    theta_right = np.asarray(fields["theta_right"], dtype=float)
    right_ux = cb * u_right - sb * w_right
    right_uy = sb * u_right + cb * w_right
    ux[n + 1 :] = right_ux[1:]
    uy[n + 1 :] = right_uy[1:]
    theta[n + 1 :] = theta_right[1:]

    joint_translation_mismatch = max(
        abs(float(ux[n]) - float(right_ux[0])),
        abs(float(uy[n]) - float(right_uy[0])),
    )
    joint_theta_mismatch = abs(float(theta[n]) - float(theta_right[0]))
    q = np.empty(3 * n_nodes, dtype=float)
    q[0::3] = ux
    q[1::3] = uy
    q[2::3] = theta
    return q, {
        "joint_translation_embedding_mismatch": float(joint_translation_mismatch),
        "joint_theta_embedding_mismatch": float(joint_theta_mismatch),
    }


def local_fields_from_q(q_full: np.ndarray, *, beta_deg: float) -> dict[str, np.ndarray]:
    q = np.asarray(q_full, dtype=float)
    n = fem.N_ELEM
    ux = q[0::3]
    uy = q[1::3]
    theta = q[2::3]
    local = compute_local_components(u_global=ux, v_global=uy, beta_deg=float(beta_deg))
    return {
        "u_left": np.asarray(local["u_left_local"], dtype=float),
        "w_left": np.asarray(local["w_left_local"], dtype=float),
        "theta_left": np.asarray(theta[: n + 1], dtype=float),
        "u_right": np.asarray(local["u_right_local"], dtype=float),
        "w_right": np.asarray(local["w_right_local"], dtype=float),
        "theta_right": np.asarray(theta[n:], dtype=float),
    }


def local_component_vector(fields: dict[str, np.ndarray]) -> np.ndarray:
    keys = ("u_left", "w_left", "theta_left", "u_right", "w_right", "theta_right")
    return np.concatenate([np.asarray(fields[key], dtype=float) for key in keys])


def numerical_derivative(values: np.ndarray, coordinates: np.ndarray) -> np.ndarray:
    y = np.asarray(values, dtype=float)
    x = np.asarray(coordinates, dtype=float)
    if y.size != x.size:
        raise ValueError("values and coordinates must have the same length.")
    if y.size < 2:
        return np.zeros_like(y, dtype=float)
    if y.size == 2:
        slope = (y[1] - y[0]) / max(x[1] - x[0], NEAR_ZERO_NORM)
        return np.full_like(y, float(slope), dtype=float)
    return np.gradient(y, x, edge_order=2)


def convention_sign_and_scale(label: str) -> tuple[float, str]:
    sign = -1.0 if label.startswith("-") else 1.0
    return sign, label[1:]


def theta_candidate_labels() -> tuple[str, ...]:
    return ("+dw/ds", "-dw/ds", "+dw/dxi", "-dw/dxi", "+dw/dz", "-dw/dz")


def theta_candidate_scale_description(label: str) -> str:
    _, base = convention_sign_and_scale(label)
    return base


def theta_convention_candidates(
    *,
    slope_z: np.ndarray,
    Lambda: float,
    length_factor: float,
    l_base: float,
) -> dict[str, np.ndarray]:
    slope = np.asarray(slope_z, dtype=float)
    Lambda_f = float(Lambda)
    length_factor_f = float(length_factor)
    l_base_f = max(float(l_base), NEAR_ZERO_NORM)
    positive = {
        "+dw/ds": (Lambda_f / l_base_f) * slope,
        "+dw/dxi": Lambda_f * length_factor_f * slope,
        "+dw/dz": slope,
    }
    out: dict[str, np.ndarray] = {}
    for label, values in positive.items():
        out[label] = np.asarray(values, dtype=float)
        out[label.replace("+", "-", 1)] = -np.asarray(values, dtype=float)
    return out


def analytic_theta_candidate_sets(
    Lambda: float,
    *,
    mu: float,
    coeff: np.ndarray,
    s_norm: np.ndarray,
    l_base: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    a1, b1, a2, b2, _, _ = [float(value) for value in coeff]
    xi = np.asarray(s_norm, dtype=float)
    left_factor = 1.0 - float(mu)
    right_factor = 1.0 + float(mu)
    z_left = float(Lambda) * left_factor * xi
    z_right = -float(Lambda) * right_factor * (1.0 - xi)
    left = theta_convention_candidates(
        slope_z=slope_values(a1, b1, z_left),
        Lambda=float(Lambda),
        length_factor=left_factor,
        l_base=float(l_base),
    )
    right = theta_convention_candidates(
        slope_z=slope_values(a2, b2, z_right),
        Lambda=float(Lambda),
        length_factor=right_factor,
        l_base=float(l_base),
    )
    return left, right


def bending_basis_mappings() -> tuple[BendingBasisMapping, ...]:
    out: list[BendingBasisMapping] = []
    for swap in (False, True):
        prefix = "swap" if swap else "identity"
        for sign_a in (1.0, -1.0):
            for sign_b in (1.0, -1.0):
                out.append(
                    BendingBasisMapping(
                        name=f"{prefix}_A{'+' if sign_a > 0 else '-'}_B{'+' if sign_b > 0 else '-'}",
                        swap=swap,
                        sign_A=float(sign_a),
                        sign_B=float(sign_b),
                    )
                )
    return tuple(out)


def canonical_bending_basis_mapping() -> BendingBasisMapping:
    return BendingBasisMapping(name="identity_A+_B+", swap=False, sign_A=1.0, sign_B=1.0)


def effective_full_basis_ab(
    A_value: float,
    B_value: float,
    mapping: BendingBasisMapping,
) -> tuple[float, float]:
    if bool(mapping.swap):
        return float(mapping.sign_B) * float(B_value), float(mapping.sign_A) * float(A_value)
    return float(mapping.sign_A) * float(A_value), float(mapping.sign_B) * float(B_value)


def full_basis_values(
    A_value: float,
    B_value: float,
    z_values: np.ndarray | float,
    mapping: BendingBasisMapping,
) -> dict[str, np.ndarray]:
    z = np.asarray(z_values, dtype=float)
    a_value, b_value = effective_full_basis_ab(float(A_value), float(B_value), mapping)
    c_value = -a_value
    d_value = -b_value
    cos_z = np.cos(z)
    sin_z = np.sin(z)
    cosh_z = np.cosh(z)
    sinh_z = np.sinh(z)
    return {
        "w": a_value * cos_z + b_value * sin_z + c_value * cosh_z + d_value * sinh_z,
        "slope": -a_value * sin_z + b_value * cos_z + c_value * sinh_z + d_value * cosh_z,
        "moment": -a_value * cos_z - b_value * sin_z + c_value * cosh_z + d_value * sinh_z,
        "third": a_value * sin_z - b_value * cos_z + c_value * sinh_z + d_value * cosh_z,
    }


def right_z_sign_value(right_z_sign: str) -> float:
    if right_z_sign == "-":
        return -1.0
    if right_z_sign == "+":
        return 1.0
    raise ValueError(f"Unknown right_z_sign: {right_z_sign!r}")


def full_basis_external_clamp_residual(
    coeff: np.ndarray,
    *,
    left_mapping: BendingBasisMapping,
    right_mapping: BendingBasisMapping,
) -> float:
    A1, B1, A2, B2, _, _ = [float(value) for value in coeff]
    left = full_basis_values(A1, B1, np.array([0.0], dtype=float), left_mapping)
    right = full_basis_values(A2, B2, np.array([0.0], dtype=float), right_mapping)
    values = [
        float(left["w"][0]),
        float(left["slope"][0]),
        0.0,
        float(right["w"][0]),
        float(right["slope"][0]),
        0.0,
    ]
    return float(max(abs(value) for value in values))


def full_basis_endpoint_residual(
    Lambda: float,
    *,
    beta_rad: float,
    mu: float,
    epsilon: float,
    coeff: np.ndarray,
    left_mapping: BendingBasisMapping,
    right_mapping: BendingBasisMapping,
    right_z_sign: str,
) -> np.ndarray:
    A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]
    left_factor = 1.0 - float(mu)
    right_factor = 1.0 + float(mu)
    x1 = float(Lambda) * left_factor
    x2 = float(Lambda) * right_factor
    th1 = float(epsilon) * float(Lambda) ** 2 * left_factor
    th2 = float(epsilon) * float(Lambda) ** 2 * right_factor
    z2 = right_z_sign_value(right_z_sign) * x2
    left = full_basis_values(A1, B1, np.array([x1], dtype=float), left_mapping)
    right = full_basis_values(A2, B2, np.array([z2], dtype=float), right_mapping)
    w1 = float(left["w"][0])
    slope1 = float(left["slope"][0])
    moment1 = float(left["moment"][0])
    shear1 = -float(epsilon) * float(Lambda) * float(left["third"][0])
    u1 = float(P1) * float(np.sin(th1))
    axial1 = float(P1) * float(np.cos(th1))

    w2 = float(right["w"][0])
    slope2 = float(right["slope"][0])
    moment2 = float(right["moment"][0])
    shear2 = -float(epsilon) * float(Lambda) * float(right["third"][0])
    u2 = float(P2) * float(np.sin(-th2))
    axial2 = float(P2) * float(np.cos(-th2))

    cb = float(np.cos(beta_rad))
    sb = float(np.sin(beta_rad))
    return np.array(
        [
            w1 - w2 * cb - u2 * sb,
            u1 + w2 * sb - u2 * cb,
            slope1 - slope2,
            moment1 - moment2,
            shear1 - shear2 * cb - axial2 * sb,
            axial1 + shear2 * sb - axial2 * cb,
        ],
        dtype=float,
    )


def full_basis_fields(
    Lambda: float,
    *,
    mu: float,
    epsilon: float,
    coeff: np.ndarray,
    s_norm: np.ndarray,
    l_base: float,
    left_mapping: BendingBasisMapping,
    right_mapping: BendingBasisMapping,
    right_z_sign: str,
) -> dict[str, np.ndarray]:
    A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]
    xi_external_to_joint = np.asarray(s_norm, dtype=float)
    left_factor = 1.0 - float(mu)
    right_factor = 1.0 + float(mu)
    left_z = float(Lambda) * left_factor * xi_external_to_joint
    right_z = right_z_sign_value(right_z_sign) * float(Lambda) * right_factor * xi_external_to_joint
    left_values = full_basis_values(A1, B1, left_z, left_mapping)
    right_values_external = full_basis_values(A2, B2, right_z, right_mapping)
    u_left = float(P1) * np.sin(float(epsilon) * float(Lambda) ** 2 * left_factor * xi_external_to_joint)
    u_right_external = float(P2) * np.sin(
        -float(epsilon) * float(Lambda) ** 2 * right_factor * xi_external_to_joint
    )
    theta_left = (float(Lambda) / max(float(l_base), NEAR_ZERO_NORM)) * np.asarray(
        left_values["slope"], dtype=float
    )
    theta_right_external_to_joint = (
        -right_z_sign_value(right_z_sign)
        * float(Lambda)
        / max(float(l_base), NEAR_ZERO_NORM)
        * np.asarray(right_values_external["slope"], dtype=float)
    )
    return {
        "u_left": np.asarray(u_left, dtype=float),
        "w_left": np.asarray(left_values["w"], dtype=float),
        "theta_left": np.asarray(theta_left, dtype=float),
        "u_right": np.asarray(u_right_external[::-1], dtype=float),
        "w_right": np.asarray(right_values_external["w"][::-1], dtype=float),
        "theta_right": np.asarray(theta_right_external_to_joint[::-1], dtype=float),
    }


def best_scaled_theta_convention(theta: np.ndarray, candidates: dict[str, np.ndarray]) -> tuple[str, float, float]:
    theta_ref = np.asarray(theta, dtype=float)
    best_label = ""
    best_scale = np.nan
    best_error = np.inf
    for label, candidate_values in candidates.items():
        candidate = np.asarray(candidate_values, dtype=float)
        scale = optimal_scalar(candidate, theta_ref)
        error = safe_relative_l2(scale * candidate - theta_ref, theta_ref)
        if np.isfinite(error) and error < best_error:
            best_label = str(label)
            best_scale = float(scale)
            best_error = float(error)
    return best_label, best_scale, best_error


def direct_theta_relative_errors(theta: np.ndarray, candidates: dict[str, np.ndarray]) -> dict[str, float]:
    theta_ref = np.asarray(theta, dtype=float)
    return {
        "dwds": safe_relative_l2(theta_ref - candidates["+dw/ds"], theta_ref),
        "minus_dwds": safe_relative_l2(theta_ref - candidates["-dw/ds"], theta_ref),
        "dwdxi": safe_relative_l2(theta_ref - candidates["+dw/dxi"], theta_ref),
        "minus_dwdxi": safe_relative_l2(theta_ref - candidates["-dw/dxi"], theta_ref),
        "dwdz": safe_relative_l2(theta_ref - candidates["+dw/dz"], theta_ref),
        "minus_dwdz": safe_relative_l2(theta_ref - candidates["-dw/dz"], theta_ref),
    }


def fem_theta_derivative_diagnostics(
    fem_local: dict[str, np.ndarray],
    *,
    Lambda: float,
    mu: float,
    l_base: float,
) -> dict[str, float | str]:
    n = fem.N_ELEM
    xi = np.linspace(0.0, 1.0, n + 1)
    left_factor = 1.0 - float(mu)
    right_factor = 1.0 + float(mu)
    s_left = max(float(l_base) * left_factor, NEAR_ZERO_NORM) * xi
    s_right = max(float(l_base) * right_factor, NEAR_ZERO_NORM) * xi
    dwdxi_left = numerical_derivative(np.asarray(fem_local["w_left"], dtype=float), xi)
    dwdxi_right = numerical_derivative(np.asarray(fem_local["w_right"], dtype=float), xi)
    dwds_left = numerical_derivative(np.asarray(fem_local["w_left"], dtype=float), s_left)
    dwds_right = numerical_derivative(np.asarray(fem_local["w_right"], dtype=float), s_right)
    left_candidates = {
        "+dw/ds": dwds_left,
        "-dw/ds": -dwds_left,
        "+dw/dxi": dwdxi_left,
        "-dw/dxi": -dwdxi_left,
        "+dw/dz": dwdxi_left / max(float(Lambda) * left_factor, NEAR_ZERO_NORM),
        "-dw/dz": -dwdxi_left / max(float(Lambda) * left_factor, NEAR_ZERO_NORM),
    }
    right_candidates = {
        "+dw/ds": dwds_right,
        "-dw/ds": -dwds_right,
        "+dw/dxi": dwdxi_right,
        "-dw/dxi": -dwdxi_right,
        "+dw/dz": dwdxi_right / max(float(Lambda) * right_factor, NEAR_ZERO_NORM),
        "-dw/dz": -dwdxi_right / max(float(Lambda) * right_factor, NEAR_ZERO_NORM),
    }
    theta_left = np.asarray(fem_local["theta_left"], dtype=float)
    theta_right = np.asarray(fem_local["theta_right"], dtype=float)
    best_left, scale_left, _ = best_scaled_theta_convention(theta_left, left_candidates)
    best_right, scale_right, _ = best_scaled_theta_convention(theta_right, right_candidates)
    left_errors = direct_theta_relative_errors(theta_left, left_candidates)
    right_errors = direct_theta_relative_errors(theta_right, right_candidates)
    return {
        "rel_l2_theta_vs_dwds_left": float(left_errors["dwds"]),
        "rel_l2_theta_vs_minus_dwds_left": float(left_errors["minus_dwds"]),
        "rel_l2_theta_vs_dwds_right": float(right_errors["dwds"]),
        "rel_l2_theta_vs_minus_dwds_right": float(right_errors["minus_dwds"]),
        "rel_l2_theta_vs_dwdxi_left": float(left_errors["dwdxi"]),
        "rel_l2_theta_vs_minus_dwdxi_left": float(left_errors["minus_dwdxi"]),
        "rel_l2_theta_vs_dwdxi_right": float(right_errors["dwdxi"]),
        "rel_l2_theta_vs_minus_dwdxi_right": float(right_errors["minus_dwdxi"]),
        "rel_l2_theta_vs_dwdz_left": float(left_errors["dwdz"]),
        "rel_l2_theta_vs_minus_dwdz_left": float(left_errors["minus_dwdz"]),
        "rel_l2_theta_vs_dwdz_right": float(right_errors["dwdz"]),
        "rel_l2_theta_vs_minus_dwdz_right": float(right_errors["minus_dwdz"]),
        "best_theta_convention_left": best_left,
        "best_theta_convention_right": best_right,
        "best_theta_scale_left": float(scale_left),
        "best_theta_scale_right": float(scale_right),
    }


def safe_relative_l2(diff: np.ndarray, reference: np.ndarray) -> float:
    denominator = float(np.linalg.norm(reference))
    return float(np.linalg.norm(diff) / denominator) if denominator > NEAR_ZERO_NORM else np.nan


def euclidean_mac(a_values: np.ndarray, b_values: np.ndarray) -> float:
    a = np.asarray(a_values, dtype=float)
    b = np.asarray(b_values, dtype=float)
    denominator = float(np.dot(a, a) * np.dot(b, b))
    return float((np.dot(a, b) ** 2) / denominator) if denominator > NEAR_ZERO_NORM else np.nan


def mass_weighted_mac(a_values: np.ndarray, b_values: np.ndarray, mass_matrix: np.ndarray) -> float:
    a = np.asarray(a_values, dtype=float)
    b = np.asarray(b_values, dtype=float)
    m = np.asarray(mass_matrix, dtype=float)
    aa = float(a @ m @ a)
    bb = float(b @ m @ b)
    ab = float(a @ m @ b)
    denominator = aa * bb
    return float((ab * ab) / denominator) if denominator > NEAR_ZERO_NORM else np.nan


def optimal_scalar(a_values: np.ndarray, b_values: np.ndarray) -> float:
    a = np.asarray(a_values, dtype=float)
    b = np.asarray(b_values, dtype=float)
    denominator = float(np.dot(a, a))
    return float(np.dot(a, b) / denominator) if denominator > NEAR_ZERO_NORM else 1.0


def residual_metrics(
    k_free: np.ndarray,
    m_free: np.ndarray,
    q_full: np.ndarray,
    free: np.ndarray,
    *,
    omega_sq: float,
) -> dict[str, float]:
    q_free = np.asarray(q_full, dtype=float)[np.asarray(free, dtype=int)]
    kq = np.asarray(k_free, dtype=float) @ q_free
    mq = np.asarray(m_free, dtype=float) @ q_free
    residual = kq - float(omega_sq) * mq
    denominator_l2 = float(np.linalg.norm(kq) + abs(float(omega_sq)) * np.linalg.norm(mq))
    denominator_inf = float(np.max(np.abs(kq)) + abs(float(omega_sq)) * np.max(np.abs(mq)))
    return {
        "residual_l2": float(np.linalg.norm(residual)),
        "relative_residual": float(np.linalg.norm(residual) / denominator_l2) if denominator_l2 > NEAR_ZERO_NORM else np.nan,
        "residual_inf": float(np.max(np.abs(residual))),
        "relative_residual_inf": float(np.max(np.abs(residual)) / denominator_inf)
        if denominator_inf > NEAR_ZERO_NORM
        else np.nan,
    }


def residual_vector_free(
    k_free: np.ndarray,
    m_free: np.ndarray,
    q_full: np.ndarray,
    free: np.ndarray,
    *,
    omega_sq: float,
) -> np.ndarray:
    q_free = np.asarray(q_full, dtype=float)[np.asarray(free, dtype=int)]
    return np.asarray(k_free, dtype=float) @ q_free - float(omega_sq) * (np.asarray(m_free, dtype=float) @ q_free)


def free_group_residual_norm(residual_free: np.ndarray, free: np.ndarray, full_mask: np.ndarray) -> float:
    selected = np.asarray(full_mask, dtype=bool)[np.asarray(free, dtype=int)]
    values = np.asarray(residual_free, dtype=float)[selected]
    return float(np.linalg.norm(values))


def safe_relative_difference(value: float, reference: float) -> float:
    return abs(float(value) - float(reference)) / abs(float(reference)) if abs(float(reference)) > NEAR_ZERO_NORM else np.nan


def rayleigh_diagnostics(
    k_free: np.ndarray,
    m_free: np.ndarray,
    q_full: np.ndarray,
    free: np.ndarray,
    *,
    fem_omega_sq: float,
    fem_lambda: float,
) -> dict[str, float]:
    q_free = np.asarray(q_full, dtype=float)[np.asarray(free, dtype=int)]
    k = np.asarray(k_free, dtype=float)
    m = np.asarray(m_free, dtype=float)
    kq = k @ q_free
    mq = m @ q_free
    denominator = float(q_free @ mq)
    omega_rayleigh_sq = float((q_free @ kq) / denominator) if abs(denominator) > NEAR_ZERO_NORM else np.nan
    rayleigh_lambda = float(omega_rayleigh_sq ** 0.25) if omega_rayleigh_sq >= 0.0 else np.nan
    return {
        "analytic_rayleigh_omega_sq": omega_rayleigh_sq,
        "fem_selected_omega_sq": float(fem_omega_sq),
        "relative_rayleigh_omega_sq_diff": safe_relative_difference(omega_rayleigh_sq, float(fem_omega_sq)),
        "analytic_rayleigh_Lambda": rayleigh_lambda,
        "relative_rayleigh_Lambda_diff": safe_relative_difference(rayleigh_lambda, float(fem_lambda)),
    }


def safe_sqrt(value: float) -> float:
    return float(np.sqrt(float(value))) if float(value) >= 0.0 else np.nan


def safe_quarter_power(value: float) -> float:
    return float(float(value) ** 0.25) if float(value) >= 0.0 else np.nan


def frequency_scaling_audit_rows(
    k_free: np.ndarray,
    m_free: np.ndarray,
    q_analytic: np.ndarray,
    q_fem: np.ndarray,
    free: np.ndarray,
    *,
    analytic_lambda: float,
    fem_selected_omega: float,
    fem_selected_omega_sq: float,
    rayleigh_omega_sq: float,
) -> list[dict[str, float | str]]:
    analytic_lambda_squared = float(analytic_lambda) ** 2
    analytic_lambda_fourth = float(analytic_lambda) ** 4
    fem_lambda_from_omega = safe_sqrt(float(fem_selected_omega))
    fem_lambda_from_omega_sq = safe_quarter_power(float(fem_selected_omega_sq))
    rayleigh_lambda = safe_quarter_power(float(rayleigh_omega_sq))
    lambda_fourth_ratio = (
        float(analytic_lambda_fourth) / float(fem_selected_omega_sq)
        if abs(float(fem_selected_omega_sq)) > NEAR_ZERO_NORM
        else np.nan
    )
    lambda_squared_ratio = (
        float(analytic_lambda_squared) / float(fem_selected_omega)
        if abs(float(fem_selected_omega)) > NEAR_ZERO_NORM
        else np.nan
    )

    common = {
        "analytic_lambda": float(analytic_lambda),
        "fem_selected_omega": float(fem_selected_omega),
        "fem_selected_omega_sq": float(fem_selected_omega_sq),
        "fem_selected_lambda_from_omega": float(fem_lambda_from_omega),
        "fem_selected_lambda_from_omega_sq": float(fem_lambda_from_omega_sq),
        "analytic_lambda_squared": float(analytic_lambda_squared),
        "analytic_lambda_fourth": float(analytic_lambda_fourth),
        "rayleigh_omega_sq": float(rayleigh_omega_sq),
        "rayleigh_lambda": float(rayleigh_lambda),
        "lambda_fourth_to_fem_omega_sq_ratio": float(lambda_fourth_ratio),
        "lambda_squared_to_fem_omega_ratio": float(lambda_squared_ratio),
        "fem_omega_minus_lambda_squared": float(fem_selected_omega - analytic_lambda_squared),
        "fem_omega_sq_minus_lambda_fourth": float(fem_selected_omega_sq - analytic_lambda_fourth),
    }

    analytic_factors = [
        ("analytic_lambda", float(analytic_lambda), "diagnostic only: factor inserted directly as Kq - factor*Mq"),
        ("analytic_lambda_squared", float(analytic_lambda_squared), "equals analytic omega if omega = Lambda^2"),
        ("analytic_lambda_fourth", float(analytic_lambda_fourth), "equals analytic omega^2 if omega = Lambda^2"),
        ("fem_selected_omega", float(fem_selected_omega), "FEM omega, not the eigenvalue in K phi = omega^2 M phi"),
        ("fem_selected_omega_squared", float(fem_selected_omega) ** 2, "omega squared from returned FEM omega"),
        ("fem_selected_omega_sq", float(fem_selected_omega_sq), "selected FEM eigenvalue for K phi = omega^2 M phi"),
        ("rayleigh_omega_sq", float(rayleigh_omega_sq), "Rayleigh-minimizing scalar for this analytic q"),
    ]
    fem_factors = [
        ("fem_selected_omega_sq", float(fem_selected_omega_sq), "selected FEM eigenvalue for K phi = omega^2 M phi"),
        ("fem_selected_omega", float(fem_selected_omega), "intentionally wrong if used as eigenvalue factor"),
        ("analytic_lambda_fourth", float(analytic_lambda_fourth), "analytic Lambda^4 mapped to omega^2"),
        ("analytic_lambda_squared", float(analytic_lambda_squared), "analytic Lambda^2 mapped to omega"),
    ]

    rows: list[dict[str, float | str]] = []

    def append_rows(vector_type: str, q_full: np.ndarray, factors: list[tuple[str, float, str]]) -> None:
        for factor_name, factor_value, note in factors:
            metrics = residual_metrics(k_free, m_free, q_full, free, omega_sq=float(factor_value))
            rows.append(
                {
                    "vector_type": vector_type,
                    "spectral_factor_name": factor_name,
                    "spectral_factor_value": float(factor_value),
                    "relative_residual": float(metrics["relative_residual"]),
                    "residual_l2": float(metrics["residual_l2"]),
                    "residual_inf": float(metrics["residual_inf"]),
                    "note": note,
                    **common,
                }
            )

    append_rows("analytic_q", np.asarray(q_analytic, dtype=float), analytic_factors)
    append_rows("fem_selected_q", np.asarray(q_fem, dtype=float), fem_factors)
    return rows


def frequency_scaling_best_row(
    rows: Sequence[dict[str, float | str]],
    *,
    vector_type: str,
) -> dict[str, float | str] | None:
    candidates = [
        row
        for row in rows
        if str(row["vector_type"]) == vector_type and np.isfinite(float(row["relative_residual"]))
    ]
    return min(candidates, key=lambda row: float(row["relative_residual"])) if candidates else None


def modal_projection_rows(
    params: BeamParams,
    *,
    mu: float,
    beta_deg: float,
    q_full: np.ndarray,
    free: np.ndarray,
    mass_matrix: np.ndarray,
    n_modes: int = DEFAULT_MODAL_PROJECTION_MODES,
) -> list[dict[str, float | int]]:
    omega_values, vecs, _ = solve_fem_modes_nd(
        params,
        mu=float(mu),
        beta_deg=float(beta_deg),
        n_modes=int(n_modes),
    )
    q_free = np.asarray(q_full, dtype=float)[np.asarray(free, dtype=int)]
    mass = np.asarray(mass_matrix, dtype=float)
    q_mass_norm_sq = float(q_free @ mass @ q_free)
    q_euclidean_norm_sq = float(q_free @ q_free)
    rows_by_mode: list[dict[str, float | int]] = []
    cumulative_by_mode = 0.0
    for idx in range(vecs.shape[1]):
        phi = np.asarray(vecs[:, idx], dtype=float)
        modal_mass = float(phi @ mass @ phi)
        modal_euclidean = float(phi @ phi)
        mass_coeff = float((phi @ mass @ q_free) / modal_mass) if abs(modal_mass) > NEAR_ZERO_NORM else np.nan
        euclidean_coeff = (
            float((phi @ q_free) / modal_euclidean) if abs(modal_euclidean) > NEAR_ZERO_NORM else np.nan
        )
        mass_share = (
            float((mass_coeff**2) * modal_mass / q_mass_norm_sq)
            if np.isfinite(mass_coeff) and abs(q_mass_norm_sq) > NEAR_ZERO_NORM
            else np.nan
        )
        euclidean_share = (
            float((euclidean_coeff**2) * modal_euclidean / q_euclidean_norm_sq)
            if np.isfinite(euclidean_coeff) and abs(q_euclidean_norm_sq) > NEAR_ZERO_NORM
            else np.nan
        )
        if np.isfinite(mass_share):
            cumulative_by_mode += float(mass_share)
        rows_by_mode.append(
            {
                "rank_by_mass_share": 0,
                "mode_index": int(idx + 1),
                "omega_nd": float(omega_values[idx]),
                "lambda_value": float(np.sqrt(max(float(omega_values[idx]), 0.0))),
                "coefficient_euclidean": euclidean_coeff,
                "coefficient_abs_euclidean": abs(euclidean_coeff) if np.isfinite(euclidean_coeff) else np.nan,
                "euclidean_projection_share": euclidean_share,
                "coefficient_mass_weighted": mass_coeff,
                "coefficient_abs_mass_weighted": abs(mass_coeff) if np.isfinite(mass_coeff) else np.nan,
                "coefficient_abs": abs(mass_coeff) if np.isfinite(mass_coeff) else np.nan,
                "coefficient_share": mass_share,
                "mass_weighted_projection_share": mass_share,
                "cumulative_mass_weighted_share_by_mode_order": float(cumulative_by_mode),
                "cumulative_mass_weighted_share_by_rank": 0.0,
                "modal_mass": modal_mass,
                "modal_euclidean_norm_sq": modal_euclidean,
            }
        )

    rows = sorted(
        rows_by_mode,
        key=lambda row: float(row["mass_weighted_projection_share"])
        if np.isfinite(float(row["mass_weighted_projection_share"]))
        else -np.inf,
        reverse=True,
    )
    cumulative_by_rank = 0.0
    for rank, row in enumerate(rows, start=1):
        share = float(row["mass_weighted_projection_share"])
        if np.isfinite(share):
            cumulative_by_rank += share
        row["rank_by_mass_share"] = int(rank)
        row["cumulative_mass_weighted_share_by_rank"] = float(cumulative_by_rank)
    return rows


def dof_mask_for_nodes(
    *,
    dof_count: int,
    nodes: Sequence[int],
    offsets: Sequence[int],
) -> np.ndarray:
    mask = np.zeros(int(dof_count), dtype=bool)
    for node in nodes:
        for offset in offsets:
            idx = 3 * int(node) + int(offset)
            if 0 <= idx < dof_count:
                mask[idx] = True
    return mask


def residual_group_rows(
    k_full: np.ndarray,
    m_full: np.ndarray,
    q_full: np.ndarray,
    free: np.ndarray,
    *,
    omega_sq: float,
) -> list[dict[str, float | int | str]]:
    n = fem.N_ELEM
    dof_count = len(q_full)
    q = np.asarray(q_full, dtype=float)
    kq = np.asarray(k_full, dtype=float) @ q
    mq = np.asarray(m_full, dtype=float) @ q
    residual = kq - float(omega_sq) * mq
    free_mask = np.zeros(dof_count, dtype=bool)
    free_mask[np.asarray(free, dtype=int)] = True
    total_free_norm = float(np.linalg.norm(residual[free_mask]))
    total_full_norm = float(np.linalg.norm(residual))
    left_internal_nodes = list(range(1, n))
    right_internal_nodes = list(range(n + 1, 2 * n))
    group_specs = [
        (
            "left_arm_translational",
            "free_equilibrium",
            dof_mask_for_nodes(dof_count=dof_count, nodes=left_internal_nodes, offsets=(0, 1)),
        ),
        (
            "left_arm_rotational",
            "free_equilibrium",
            dof_mask_for_nodes(dof_count=dof_count, nodes=left_internal_nodes, offsets=(2,)),
        ),
        (
            "right_arm_translational",
            "free_equilibrium",
            dof_mask_for_nodes(dof_count=dof_count, nodes=right_internal_nodes, offsets=(0, 1)),
        ),
        (
            "right_arm_rotational",
            "free_equilibrium",
            dof_mask_for_nodes(dof_count=dof_count, nodes=right_internal_nodes, offsets=(2,)),
        ),
        ("joint_dofs", "free_equilibrium", dof_mask_for_nodes(dof_count=dof_count, nodes=(n,), offsets=(0, 1, 2))),
        (
            "left_external_clamp_dofs",
            "constraint_reaction",
            dof_mask_for_nodes(dof_count=dof_count, nodes=(0,), offsets=(0, 1, 2)),
        ),
        (
            "right_external_clamp_dofs",
            "constraint_reaction",
            dof_mask_for_nodes(dof_count=dof_count, nodes=(2 * n,), offsets=(0, 1, 2)),
        ),
        (
            "external_clamp_dofs",
            "constraint_reaction",
            dof_mask_for_nodes(dof_count=dof_count, nodes=(0, 2 * n), offsets=(0, 1, 2)),
        ),
    ]
    rows: list[dict[str, float | int | str]] = []
    for name, scope, mask in group_specs:
        selected = np.asarray(mask, dtype=bool)
        if scope == "free_equilibrium":
            selected = selected & free_mask
            share_denominator = total_free_norm
        else:
            share_denominator = total_full_norm
        group_residual = residual[selected]
        group_kq = kq[selected]
        group_mq = mq[selected]
        residual_norm = float(np.linalg.norm(group_residual))
        denominator = float(np.linalg.norm(group_kq) + abs(float(omega_sq)) * np.linalg.norm(group_mq))
        rows.append(
            {
                "group_name": name,
                "residual_scope": scope,
                "dof_count": int(np.count_nonzero(selected)),
                "group_residual_norm": residual_norm,
                "group_residual_share": float(residual_norm / share_denominator)
                if share_denominator > NEAR_ZERO_NORM
                else np.nan,
                "group_relative_residual": float(residual_norm / denominator)
                if denominator > NEAR_ZERO_NORM
                else np.nan,
                "group_residual_inf": float(np.max(np.abs(group_residual))) if group_residual.size else 0.0,
                "group_kq_norm": float(np.linalg.norm(group_kq)),
                "group_omega_mq_norm": float(abs(float(omega_sq)) * np.linalg.norm(group_mq)),
            }
        )
    return rows


def local_stiffness_energy(q_local: np.ndarray, k_local: np.ndarray, indices: np.ndarray) -> float:
    q_block = np.asarray(q_local, dtype=float)[indices]
    k_block = np.asarray(k_local, dtype=float)[np.ix_(indices, indices)]
    return float(0.5 * q_block @ k_block @ q_block)


def element_residual_energy_rows(
    params: BeamParams,
    *,
    mu: float,
    beta_deg: float,
    q_analytic: np.ndarray,
    q_fem: np.ndarray,
    omega_sq: float,
) -> list[dict[str, float | int | str]]:
    n = fem.N_ELEM
    axial_idx = np.array([0, 3], dtype=int)
    bending_idx = np.array([1, 2, 4, 5], dtype=int)
    rows: list[dict[str, float | int | str]] = []
    with fem_parameter_override(params):
        length_left = max(float(fem.ell) * (1.0 - float(mu)), 1e-3)
        length_right = max(float(fem.ell) * (1.0 + float(mu)), 1e-3)
        le_left = length_left / n
        le_right = length_right / n
        k_left = np.asarray(fem.elem_K(le_left), dtype=float)
        m_left = np.asarray(fem.elem_M(le_left), dtype=float)
        k_right = np.asarray(fem.elem_K(le_right), dtype=float)
        m_right = np.asarray(fem.elem_M(le_right), dtype=float)
        t_left = fem.rotation_matrix_6x6(0.0)
        t_right = fem.rotation_matrix_6x6(np.deg2rad(float(beta_deg)))

    def append_element(arm: str, elem_idx: int, node0: int, node1: int) -> None:
        dofs = np.array(
            [
                3 * node0,
                3 * node0 + 1,
                3 * node0 + 2,
                3 * node1,
                3 * node1 + 1,
                3 * node1 + 2,
            ],
            dtype=int,
        )
        if arm == "left":
            transform = t_left
            k_local = k_left
            m_local = m_left
        else:
            transform = t_right
            k_local = k_right
            m_local = m_right
        q_analytic_local = transform @ np.asarray(q_analytic, dtype=float)[dofs]
        q_fem_local = transform @ np.asarray(q_fem, dtype=float)[dofs]
        local_residual = k_local @ q_analytic_local - float(omega_sq) * (m_local @ q_analytic_local)
        rows.append(
            {
                "arm": arm,
                "element_index": int(elem_idx),
                "global_node0": int(node0),
                "global_node1": int(node1),
                "midpoint_s_norm": float((elem_idx + 0.5) / n),
                "analytic_local_u1": float(q_analytic_local[0]),
                "analytic_local_v1": float(q_analytic_local[1]),
                "analytic_local_theta1": float(q_analytic_local[2]),
                "analytic_local_u2": float(q_analytic_local[3]),
                "analytic_local_v2": float(q_analytic_local[4]),
                "analytic_local_theta2": float(q_analytic_local[5]),
                "fem_local_u1": float(q_fem_local[0]),
                "fem_local_v1": float(q_fem_local[1]),
                "fem_local_theta1": float(q_fem_local[2]),
                "fem_local_u2": float(q_fem_local[3]),
                "fem_local_v2": float(q_fem_local[4]),
                "fem_local_theta2": float(q_fem_local[5]),
                "analytic_local_axial_energy": local_stiffness_energy(q_analytic_local, k_local, axial_idx),
                "analytic_local_bending_energy": local_stiffness_energy(q_analytic_local, k_local, bending_idx),
                "fem_local_axial_energy": local_stiffness_energy(q_fem_local, k_local, axial_idx),
                "fem_local_bending_energy": local_stiffness_energy(q_fem_local, k_local, bending_idx),
                "local_residual_norm": float(np.linalg.norm(local_residual)),
                "local_residual_norm_share": 0.0,
                "local_residual_norm_sq_share": 0.0,
                "local_residual_inf": float(np.max(np.abs(local_residual))),
                "local_residual_axial_norm": float(np.linalg.norm(local_residual[axial_idx])),
                "local_residual_bending_norm": float(np.linalg.norm(local_residual[bending_idx])),
            }
        )

    for elem in range(n):
        append_element("left", elem, elem, elem + 1)
    for elem in range(n):
        append_element("right", elem, n + elem, n + elem + 1)

    total_norm = sum(float(row["local_residual_norm"]) for row in rows)
    total_norm_sq = sum(float(row["local_residual_norm"]) ** 2 for row in rows)
    for row in rows:
        norm = float(row["local_residual_norm"])
        row["local_residual_norm_share"] = float(norm / total_norm) if total_norm > NEAR_ZERO_NORM else np.nan
        row["local_residual_norm_sq_share"] = (
            float(norm**2 / total_norm_sq) if total_norm_sq > NEAR_ZERO_NORM else np.nan
        )
    return rows


def q_kinematic_checks(q_full: np.ndarray, embedding: dict[str, float]) -> dict[str, float]:
    n = fem.N_ELEM
    q = np.asarray(q_full, dtype=float)
    left_clamp = q[0:3]
    right_clamp = q[3 * (2 * n) : 3 * (2 * n) + 3]
    return {
        "q_joint_translation_mismatch": float(embedding["joint_translation_embedding_mismatch"]),
        "q_joint_rotation_mismatch": float(embedding["joint_theta_embedding_mismatch"]),
        "q_external_clamp_max_abs": float(max(np.max(np.abs(left_clamp)), np.max(np.abs(right_clamp)))),
    }


def axial_right_u_variant(u_right: np.ndarray, mode: str) -> np.ndarray:
    values = np.asarray(u_right, dtype=float)
    if mode == "as_is":
        return values.copy()
    if mode == "negative":
        return -values
    if mode == "reversed":
        return values[::-1].copy()
    if mode == "reversed_negative":
        return -values[::-1]
    raise ValueError(f"Unknown right_u_mode: {mode!r}")


def axial_variant_fields(
    fields_base: dict[str, np.ndarray],
    *,
    left_u_sign: float,
    right_u_mode: str,
) -> dict[str, np.ndarray]:
    return {
        "u_left": float(left_u_sign) * np.asarray(fields_base["u_left"], dtype=float),
        "w_left": np.asarray(fields_base["w_left"], dtype=float).copy(),
        "theta_left": np.asarray(fields_base["theta_left"], dtype=float).copy(),
        "u_right": axial_right_u_variant(np.asarray(fields_base["u_right"], dtype=float), right_u_mode),
        "w_right": np.asarray(fields_base["w_right"], dtype=float).copy(),
        "theta_right": np.asarray(fields_base["theta_right"], dtype=float).copy(),
    }


def element_energy_sums(rows: Sequence[dict[str, float | int | str]]) -> dict[str, float]:
    left_axial = sum(float(row["analytic_local_axial_energy"]) for row in rows if row["arm"] == "left")
    right_axial = sum(float(row["analytic_local_axial_energy"]) for row in rows if row["arm"] == "right")
    left_bending = sum(float(row["analytic_local_bending_energy"]) for row in rows if row["arm"] == "left")
    right_bending = sum(float(row["analytic_local_bending_energy"]) for row in rows if row["arm"] == "right")
    fem_left_axial = sum(float(row["fem_local_axial_energy"]) for row in rows if row["arm"] == "left")
    fem_right_axial = sum(float(row["fem_local_axial_energy"]) for row in rows if row["arm"] == "right")
    fem_left_bending = sum(float(row["fem_local_bending_energy"]) for row in rows if row["arm"] == "left")
    fem_right_bending = sum(float(row["fem_local_bending_energy"]) for row in rows if row["arm"] == "right")
    total_axial = left_axial + right_axial
    total_bending = left_bending + right_bending
    fem_total_axial = fem_left_axial + fem_right_axial
    fem_total_bending = fem_left_bending + fem_right_bending
    axial_residual = sum(float(row["local_residual_axial_norm"]) for row in rows)
    bending_residual = sum(float(row["local_residual_bending_norm"]) for row in rows)
    return {
        "left_axial_energy": float(left_axial),
        "right_axial_energy": float(right_axial),
        "right_axial_share": float(right_axial / total_axial) if total_axial > NEAR_ZERO_NORM else np.nan,
        "total_axial_fraction": float(total_axial / (total_axial + total_bending))
        if total_axial + total_bending > NEAR_ZERO_NORM
        else np.nan,
        "fem_right_axial_energy": float(fem_right_axial),
        "fem_total_axial_fraction": float(fem_total_axial / (fem_total_axial + fem_total_bending))
        if fem_total_axial + fem_total_bending > NEAR_ZERO_NORM
        else np.nan,
        "axial_to_bending_residual_ratio": float(axial_residual / bending_residual)
        if bending_residual > NEAR_ZERO_NORM
        else np.nan,
    }


def axial_strain_diagnostics(
    q_analytic: np.ndarray,
    q_fem: np.ndarray,
    *,
    beta_deg: float,
    mu: float,
    l_base: float,
) -> dict[str, float | str]:
    analytic_local = local_fields_from_q(q_analytic, beta_deg=float(beta_deg))
    fem_local = local_fields_from_q(q_fem, beta_deg=float(beta_deg))
    n = fem.N_ELEM
    xi = np.linspace(0.0, 1.0, n + 1)
    s_left = max(float(l_base) * (1.0 - float(mu)), NEAR_ZERO_NORM) * xi
    s_right = max(float(l_base) * (1.0 + float(mu)), NEAR_ZERO_NORM) * xi
    du_left = numerical_derivative(np.asarray(analytic_local["u_left"], dtype=float), s_left)
    du_right = numerical_derivative(np.asarray(analytic_local["u_right"], dtype=float), s_right)
    fem_du_left = numerical_derivative(np.asarray(fem_local["u_left"], dtype=float), s_left)
    fem_du_right = numerical_derivative(np.asarray(fem_local["u_right"], dtype=float), s_right)
    right_candidates = {
        "+du/ds": du_right,
        "-du/ds": -du_right,
        "reversed +du/ds": du_right[::-1],
        "reversed -du/ds": -du_right[::-1],
    }
    best_label = ""
    best_error = np.inf
    for label, candidate in right_candidates.items():
        error = safe_relative_l2(np.asarray(candidate, dtype=float) - fem_du_right, fem_du_right)
        if np.isfinite(error) and error < best_error:
            best_label = label
            best_error = float(error)
    return {
        "rel_l2_du_ds_left": safe_relative_l2(du_left - fem_du_left, fem_du_left),
        "rel_l2_du_ds_right": safe_relative_l2(du_right - fem_du_right, fem_du_right),
        "max_abs_du_ds_left": float(np.max(np.abs(du_left - fem_du_left))),
        "max_abs_du_ds_right": float(np.max(np.abs(du_right - fem_du_right))),
        "best_du_ds_right_convention": best_label,
        "best_du_ds_right_rel_l2": float(best_error),
    }


def residual_group_share_by_name(
    rows: Sequence[dict[str, float | int | str]],
    group_name: str,
) -> float:
    for row in rows:
        if str(row["group_name"]) == group_name:
            return float(row["group_residual_share"])
    return np.nan


def axial_convention_conclusion(
    *,
    variant_name: str,
    residual: float,
    baseline_residual: float,
    rayleigh_diff: float,
    baseline_rayleigh_diff: float,
) -> str:
    if variant_name == "left+_right_as_is":
        return "baseline"
    residual_dramatic = np.isfinite(residual) and float(residual) < 0.1 * float(baseline_residual)
    rayleigh_dramatic = (
        np.isfinite(rayleigh_diff)
        and np.isfinite(baseline_rayleigh_diff)
        and float(rayleigh_diff) < 0.1 * float(baseline_rayleigh_diff)
    )
    if residual_dramatic or rayleigh_dramatic:
        return "possible_axial_convention_mismatch"
    if np.isfinite(residual) and float(residual) < float(baseline_residual) * (1.0 - 1e-6):
        return "needs_review"
    if (
        np.isfinite(rayleigh_diff)
        and np.isfinite(baseline_rayleigh_diff)
        and float(rayleigh_diff) < float(baseline_rayleigh_diff) * (1.0 - 1e-6)
    ):
        return "needs_review"
    return "no_improvement"


def axial_convention_scan_rows(
    *,
    fields_base: dict[str, np.ndarray],
    beta_deg: float,
    mu: float,
    l_base: float,
    params: BeamParams,
    k_free: np.ndarray,
    m_free: np.ndarray,
    k_full: np.ndarray,
    m_full: np.ndarray,
    free: np.ndarray,
    omega_sq: float,
    fem_lambda: float,
    q_fem: np.ndarray,
    translational_mask: np.ndarray,
    rotational_mask: np.ndarray,
) -> list[dict[str, float | str]]:
    right_modes = ("as_is", "negative", "reversed", "reversed_negative")
    rows: list[dict[str, float | str]] = []
    q_fem_full = np.asarray(q_fem, dtype=float)
    q_fem_free = q_fem_full[np.asarray(free, dtype=int)]
    for left_u_sign in (1.0, -1.0):
        for right_u_mode in right_modes:
            variant_name = f"left{'+' if left_u_sign > 0 else '-'}_right_{right_u_mode}"
            fields = axial_variant_fields(
                fields_base,
                left_u_sign=float(left_u_sign),
                right_u_mode=right_u_mode,
            )
            q_analytic, _ = build_analytic_global_q(fields, beta_deg=float(beta_deg))
            q_analytic_free = q_analytic[np.asarray(free, dtype=int)]
            alignment_scalar = optimal_scalar(q_analytic_free, q_fem_free)
            q_aligned = alignment_scalar * q_analytic
            q_aligned_free = q_aligned[np.asarray(free, dtype=int)]
            residual = residual_metrics(k_free, m_free, q_analytic, free, omega_sq=float(omega_sq))
            rayleigh = rayleigh_diagnostics(
                k_free,
                m_free,
                q_analytic,
                free,
                fem_omega_sq=float(omega_sq),
                fem_lambda=float(fem_lambda),
            )
            group_rows = residual_group_rows(k_full, m_full, q_aligned, free, omega_sq=float(omega_sq))
            element_rows = element_residual_energy_rows(
                params,
                mu=float(mu),
                beta_deg=float(beta_deg),
                q_analytic=q_aligned,
                q_fem=q_fem_full,
                omega_sq=float(omega_sq),
            )
            energy = element_energy_sums(element_rows)
            strain = axial_strain_diagnostics(
                q_aligned,
                q_fem_full,
                beta_deg=float(beta_deg),
                mu=float(mu),
                l_base=float(l_base),
            )
            right_axial_energy = float(energy["right_axial_energy"])
            fem_right_axial_energy = float(energy["fem_right_axial_energy"])
            rows.append(
                {
                    "variant_name": variant_name,
                    "left_u_sign": float(left_u_sign),
                    "right_u_mode": right_u_mode,
                    "relative_residual": float(residual["relative_residual"]),
                    "rayleigh_omega_sq": float(rayleigh["analytic_rayleigh_omega_sq"]),
                    "rayleigh_lambda": float(rayleigh["analytic_rayleigh_Lambda"]),
                    "relative_rayleigh_lambda_diff": float(rayleigh["relative_rayleigh_Lambda_diff"]),
                    "global_dof_mac_euclidean": euclidean_mac(q_analytic_free, q_fem_free),
                    "mass_weighted_mac": mass_weighted_mac(q_analytic_free, q_fem_free, m_free),
                    "translational_relative_l2": subset_relative_l2(q_aligned, q_fem_full, translational_mask),
                    "rotational_relative_l2": subset_relative_l2(q_aligned, q_fem_full, rotational_mask),
                    "right_arm_translational_residual_share": residual_group_share_by_name(
                        group_rows, "right_arm_translational"
                    ),
                    "axial_to_bending_residual_ratio": float(energy["axial_to_bending_residual_ratio"]),
                    "analytic_right_axial_energy": right_axial_energy,
                    "fem_right_axial_energy": fem_right_axial_energy,
                    "right_axial_energy_ratio": float(right_axial_energy / fem_right_axial_energy)
                    if fem_right_axial_energy > NEAR_ZERO_NORM
                    else np.nan,
                    "analytic_total_axial_fraction": float(energy["total_axial_fraction"]),
                    "fem_total_axial_fraction": float(energy["fem_total_axial_fraction"]),
                    "left_axial_energy": float(energy["left_axial_energy"]),
                    "right_axial_share": float(energy["right_axial_share"]),
                    "rel_l2_du_ds_left": float(strain["rel_l2_du_ds_left"]),
                    "rel_l2_du_ds_right": float(strain["rel_l2_du_ds_right"]),
                    "max_abs_du_ds_left": float(strain["max_abs_du_ds_left"]),
                    "max_abs_du_ds_right": float(strain["max_abs_du_ds_right"]),
                    "best_du_ds_right_convention": str(strain["best_du_ds_right_convention"]),
                    "best_du_ds_right_rel_l2": float(strain["best_du_ds_right_rel_l2"]),
                    "conclusion_flag": "",
                }
            )

    baseline = next(row for row in rows if row["variant_name"] == "left+_right_as_is")
    baseline_residual = float(baseline["relative_residual"])
    baseline_rayleigh_diff = float(baseline["relative_rayleigh_lambda_diff"])
    for row in rows:
        row["conclusion_flag"] = axial_convention_conclusion(
            variant_name=str(row["variant_name"]),
            residual=float(row["relative_residual"]),
            baseline_residual=baseline_residual,
            rayleigh_diff=float(row["relative_rayleigh_lambda_diff"]),
            baseline_rayleigh_diff=baseline_rayleigh_diff,
        )
    rows.sort(key=lambda row: float(row["relative_residual"]))
    return rows


def write_axial_convention_scan_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=AXIAL_CONVENTION_SCAN_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def print_axial_convention_scan_top(rows: Sequence[dict[str, float | str]], *, limit: int = 6) -> None:
    if not rows:
        print("axial convention scan: no rows")
        return
    baseline = next((row for row in rows if row["variant_name"] == "left+_right_as_is"), None)
    if baseline is not None:
        print(
            "baseline axial convention: "
            f"residual={float(baseline['relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(baseline['rayleigh_lambda']):.6e}"
        )
    print("Top axial convention variants by analytic-in-FEM relative residual:")
    for row in list(rows)[:limit]:
        print(
            f"  residual={float(row['relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(row['rayleigh_lambda']):.6e}, "
            f"right_share={float(row['right_arm_translational_residual_share']):.6e}, "
            f"{row['variant_name']} ({row['conclusion_flag']})"
        )
    rayleigh_sorted = sorted(
        rows,
        key=lambda row: float(row["relative_rayleigh_lambda_diff"])
        if np.isfinite(float(row["relative_rayleigh_lambda_diff"]))
        else np.inf,
    )
    print("Top axial convention variants by Rayleigh Lambda closeness:")
    for row in rayleigh_sorted[:limit]:
        print(
            f"  rel_Lambda_diff={float(row['relative_rayleigh_lambda_diff']):.6e}, "
            f"Rayleigh Lambda={float(row['rayleigh_lambda']):.6e}, "
            f"residual={float(row['relative_residual']):.6e}, {row['variant_name']} "
            f"({row['conclusion_flag']})"
        )
    if any(str(row["conclusion_flag"]) == "possible_axial_convention_mismatch" for row in rows):
        print("possible axial convention mismatch")
    else:
        print("simple sign/reversal axial convention does not explain mismatch")


def axial_scale_variant_fields(
    fields_base: dict[str, np.ndarray],
    *,
    alpha_left: float,
    alpha_right: float,
) -> dict[str, np.ndarray]:
    return {
        "u_left": float(alpha_left) * np.asarray(fields_base["u_left"], dtype=float),
        "w_left": np.asarray(fields_base["w_left"], dtype=float).copy(),
        "theta_left": np.asarray(fields_base["theta_left"], dtype=float).copy(),
        "u_right": float(alpha_right) * np.asarray(fields_base["u_right"], dtype=float),
        "w_right": np.asarray(fields_base["w_right"], dtype=float).copy(),
        "theta_right": np.asarray(fields_base["theta_right"], dtype=float).copy(),
    }


def combine_scale_label(existing: str, label: str) -> str:
    if not existing:
        return str(label)
    labels = [part.strip() for part in existing.split("|")]
    return existing if str(label) in labels else f"{existing}|{label}"


def axial_scale_candidate_values(
    *,
    Lambda: float,
    epsilon: float,
    mu: float,
) -> list[tuple[str, float]]:
    Lambda_f = max(abs(float(Lambda)), NEAR_ZERO_NORM)
    epsilon_f = float(epsilon)
    left_factor = max(1.0 - float(mu), NEAR_ZERO_NORM)
    right_factor = max(1.0 + float(mu), NEAR_ZERO_NORM)
    raw = [
        ("0", 0.0),
        ("epsilon", epsilon_f),
        ("sqrt(epsilon)", float(np.sqrt(epsilon_f))),
        ("epsilon*Lambda", epsilon_f * Lambda_f),
        ("epsilon*Lambda^2", epsilon_f * Lambda_f**2),
        ("1/Lambda", 1.0 / Lambda_f),
        ("1/Lambda^2", 1.0 / Lambda_f**2),
        ("1/(epsilon*Lambda)", 1.0 / (epsilon_f * Lambda_f)),
        ("1/(epsilon*Lambda^2)", 1.0 / (epsilon_f * Lambda_f**2)),
        ("1", 1.0),
        ("Lambda", Lambda_f),
        ("Lambda^2", Lambda_f**2),
        ("1/(1-mu)", 1.0 / left_factor),
        ("1/(1+mu)", 1.0 / right_factor),
        ("epsilon/(1-mu)", epsilon_f / left_factor),
        ("epsilon/(1+mu)", epsilon_f / right_factor),
    ]
    candidates: list[tuple[str, float]] = []
    for label, value in raw:
        value_f = float(value)
        if not np.isfinite(value_f) or value_f < 0.0:
            continue
        match_idx = None
        for idx, (_, existing_value) in enumerate(candidates):
            if np.isclose(value_f, existing_value, rtol=1e-12, atol=1e-15):
                match_idx = idx
                break
        if match_idx is None:
            candidates.append((str(label), value_f))
        else:
            old_label, old_value = candidates[match_idx]
            candidates[match_idx] = (combine_scale_label(old_label, str(label)), old_value)
    return candidates


def optimal_l2_relative_alpha(
    analytic_values: np.ndarray,
    fem_values: np.ndarray,
) -> float:
    a = np.asarray(analytic_values, dtype=float)
    b = np.asarray(fem_values, dtype=float)
    denominator = float(np.dot(a, a))
    return float(np.dot(a, b) / denominator) if denominator > NEAR_ZERO_NORM else np.nan


def l2_optimal_axial_alphas(
    fields_base: dict[str, np.ndarray],
    *,
    beta_deg: float,
    q_fem: np.ndarray,
    free: np.ndarray,
) -> tuple[float, float]:
    q_base, _ = build_analytic_global_q(fields_base, beta_deg=float(beta_deg))
    base_scalar = optimal_scalar(q_base[np.asarray(free, dtype=int)], np.asarray(q_fem, dtype=float)[np.asarray(free, dtype=int)])
    base_local = local_fields_from_q(base_scalar * q_base, beta_deg=float(beta_deg))
    fem_local = local_fields_from_q(q_fem, beta_deg=float(beta_deg))
    alpha_left = optimal_l2_relative_alpha(base_local["u_left"], fem_local["u_left"])
    alpha_right = optimal_l2_relative_alpha(base_local["u_right"], fem_local["u_right"])
    return float(alpha_left), float(alpha_right)


def element_diagnostic_context(
    params: BeamParams,
    *,
    mu: float,
    beta_deg: float,
) -> dict[str, np.ndarray]:
    n = fem.N_ELEM
    with fem_parameter_override(params):
        length_left = max(float(fem.ell) * (1.0 - float(mu)), 1e-3)
        length_right = max(float(fem.ell) * (1.0 + float(mu)), 1e-3)
        le_left = length_left / n
        le_right = length_right / n
        return {
            "k_left": np.asarray(fem.elem_K(le_left), dtype=float),
            "m_left": np.asarray(fem.elem_M(le_left), dtype=float),
            "k_right": np.asarray(fem.elem_K(le_right), dtype=float),
            "m_right": np.asarray(fem.elem_M(le_right), dtype=float),
            "t_left": fem.rotation_matrix_6x6(0.0),
            "t_right": fem.rotation_matrix_6x6(np.deg2rad(float(beta_deg))),
            "axial_idx": np.array([0, 3], dtype=int),
            "bending_idx": np.array([1, 2, 4, 5], dtype=int),
        }


def element_energy_sums_from_context(
    context: dict[str, np.ndarray],
    *,
    q_analytic: np.ndarray,
    q_fem: np.ndarray,
    omega_sq: float,
) -> dict[str, float]:
    n = fem.N_ELEM
    axial_idx = np.asarray(context["axial_idx"], dtype=int)
    bending_idx = np.asarray(context["bending_idx"], dtype=int)
    q_analytic_full = np.asarray(q_analytic, dtype=float)
    q_fem_full = np.asarray(q_fem, dtype=float)
    sums = {
        "left_axial": 0.0,
        "right_axial": 0.0,
        "left_bending": 0.0,
        "right_bending": 0.0,
        "fem_left_axial": 0.0,
        "fem_right_axial": 0.0,
        "fem_left_bending": 0.0,
        "fem_right_bending": 0.0,
        "axial_residual": 0.0,
        "bending_residual": 0.0,
    }

    def add_element(arm: str, node0: int, node1: int) -> None:
        dofs = np.array(
            [
                3 * node0,
                3 * node0 + 1,
                3 * node0 + 2,
                3 * node1,
                3 * node1 + 1,
                3 * node1 + 2,
            ],
            dtype=int,
        )
        if arm == "left":
            transform = np.asarray(context["t_left"], dtype=float)
            k_local = np.asarray(context["k_left"], dtype=float)
            m_local = np.asarray(context["m_left"], dtype=float)
        else:
            transform = np.asarray(context["t_right"], dtype=float)
            k_local = np.asarray(context["k_right"], dtype=float)
            m_local = np.asarray(context["m_right"], dtype=float)
        q_analytic_local = transform @ q_analytic_full[dofs]
        q_fem_local = transform @ q_fem_full[dofs]
        local_residual = k_local @ q_analytic_local - float(omega_sq) * (m_local @ q_analytic_local)
        sums[f"{arm}_axial"] += local_stiffness_energy(q_analytic_local, k_local, axial_idx)
        sums[f"{arm}_bending"] += local_stiffness_energy(q_analytic_local, k_local, bending_idx)
        sums[f"fem_{arm}_axial"] += local_stiffness_energy(q_fem_local, k_local, axial_idx)
        sums[f"fem_{arm}_bending"] += local_stiffness_energy(q_fem_local, k_local, bending_idx)
        sums["axial_residual"] += float(np.linalg.norm(local_residual[axial_idx]))
        sums["bending_residual"] += float(np.linalg.norm(local_residual[bending_idx]))

    for elem in range(n):
        add_element("left", elem, elem + 1)
    for elem in range(n):
        add_element("right", n + elem, n + elem + 1)

    total_axial = sums["left_axial"] + sums["right_axial"]
    total_bending = sums["left_bending"] + sums["right_bending"]
    fem_total_axial = sums["fem_left_axial"] + sums["fem_right_axial"]
    fem_total_bending = sums["fem_left_bending"] + sums["fem_right_bending"]
    return {
        "left_axial_energy": float(sums["left_axial"]),
        "right_axial_energy": float(sums["right_axial"]),
        "right_axial_share": float(sums["right_axial"] / total_axial) if total_axial > NEAR_ZERO_NORM else np.nan,
        "total_axial_fraction": float(total_axial / (total_axial + total_bending))
        if total_axial + total_bending > NEAR_ZERO_NORM
        else np.nan,
        "fem_right_axial_energy": float(sums["fem_right_axial"]),
        "fem_total_axial_fraction": float(fem_total_axial / (fem_total_axial + fem_total_bending))
        if fem_total_axial + fem_total_bending > NEAR_ZERO_NORM
        else np.nan,
        "axial_to_bending_residual_ratio": float(sums["axial_residual"] / sums["bending_residual"])
        if sums["bending_residual"] > NEAR_ZERO_NORM
        else np.nan,
    }


def axial_scale_conclusion(
    *,
    variant_name: str,
    scale_source: str,
    alpha_left: float,
    alpha_right: float,
    residual: float,
    baseline_residual: float,
    rayleigh_diff: float,
    baseline_rayleigh_diff: float,
) -> str:
    if variant_name == "baseline_alpha_1_1":
        return "baseline"
    residual_dramatic = np.isfinite(residual) and float(residual) < 0.1 * float(baseline_residual)
    rayleigh_dramatic = (
        np.isfinite(rayleigh_diff)
        and np.isfinite(baseline_rayleigh_diff)
        and float(rayleigh_diff) < 0.1 * float(baseline_rayleigh_diff)
    )
    both_zero = abs(float(alpha_left)) <= NEAR_ZERO_NORM and abs(float(alpha_right)) <= NEAR_ZERO_NORM
    if (residual_dramatic or rayleigh_dramatic) and both_zero:
        return "zero_axial_improves_rayleigh"
    if (residual_dramatic or rayleigh_dramatic) and str(scale_source).startswith("candidate_grid"):
        return "possible_axial_amplitude_scaling_mismatch"
    if residual_dramatic or rayleigh_dramatic:
        return "needs_review"
    if np.isfinite(residual) and float(residual) < float(baseline_residual) * 0.9:
        return "needs_review"
    if (
        np.isfinite(rayleigh_diff)
        and np.isfinite(baseline_rayleigh_diff)
        and float(rayleigh_diff) < float(baseline_rayleigh_diff) * 0.9
    ):
        return "needs_review"
    return "no_improvement"


def evaluate_axial_scale_variant(
    *,
    fields_base: dict[str, np.ndarray],
    beta_deg: float,
    params: BeamParams,
    k_free: np.ndarray,
    m_free: np.ndarray,
    k_full: np.ndarray,
    m_full: np.ndarray,
    free: np.ndarray,
    omega_sq: float,
    fem_lambda: float,
    q_fem: np.ndarray,
    translational_mask: np.ndarray,
    rotational_mask: np.ndarray,
    element_context: dict[str, np.ndarray],
    variant_name: str,
    scale_source: str,
    alpha_left: float,
    alpha_right: float,
) -> dict[str, float | str]:
    del params
    fields = axial_scale_variant_fields(
        fields_base,
        alpha_left=float(alpha_left),
        alpha_right=float(alpha_right),
    )
    q_analytic, _ = build_analytic_global_q(fields, beta_deg=float(beta_deg))
    q_fem_full = np.asarray(q_fem, dtype=float)
    free_idx = np.asarray(free, dtype=int)
    q_analytic_free = q_analytic[free_idx]
    q_fem_free = q_fem_full[free_idx]
    alignment_scalar = optimal_scalar(q_analytic_free, q_fem_free)
    q_aligned = alignment_scalar * q_analytic
    residual = residual_metrics(k_free, m_free, q_analytic, free_idx, omega_sq=float(omega_sq))
    rayleigh = rayleigh_diagnostics(
        k_free,
        m_free,
        q_analytic,
        free_idx,
        fem_omega_sq=float(omega_sq),
        fem_lambda=float(fem_lambda),
    )
    group_rows = residual_group_rows(k_full, m_full, q_aligned, free_idx, omega_sq=float(omega_sq))
    energy = element_energy_sums_from_context(
        element_context,
        q_analytic=q_aligned,
        q_fem=q_fem_full,
        omega_sq=float(omega_sq),
    )
    right_axial_energy = float(energy["right_axial_energy"])
    fem_right_axial_energy = float(energy["fem_right_axial_energy"])
    return {
        "variant_name": str(variant_name),
        "alpha_left": float(alpha_left),
        "alpha_right": float(alpha_right),
        "scale_source": str(scale_source),
        "relative_residual": float(residual["relative_residual"]),
        "rayleigh_omega_sq": float(rayleigh["analytic_rayleigh_omega_sq"]),
        "rayleigh_lambda": float(rayleigh["analytic_rayleigh_Lambda"]),
        "rayleigh_lambda_relative_diff": float(rayleigh["relative_rayleigh_Lambda_diff"]),
        "global_dof_mac_euclidean": euclidean_mac(q_analytic_free, q_fem_free),
        "mass_weighted_mac": mass_weighted_mac(q_analytic_free, q_fem_free, m_free),
        "translational_relative_l2": subset_relative_l2(q_aligned, q_fem_full, translational_mask),
        "rotational_relative_l2": subset_relative_l2(q_aligned, q_fem_full, rotational_mask),
        "right_arm_translational_residual_share": residual_group_share_by_name(group_rows, "right_arm_translational"),
        "axial_to_bending_residual_ratio": float(energy["axial_to_bending_residual_ratio"]),
        "analytic_total_axial_fraction": float(energy["total_axial_fraction"]),
        "analytic_right_axial_share": float(energy["right_axial_share"]),
        "right_axial_energy_ratio_to_fem": float(right_axial_energy / fem_right_axial_energy)
        if fem_right_axial_energy > NEAR_ZERO_NORM
        else np.nan,
        "conclusion_flag": "",
    }


def axial_scale_scan_rows(
    *,
    fields_base: dict[str, np.ndarray],
    Lambda: float,
    beta_deg: float,
    mu: float,
    epsilon: float,
    params: BeamParams,
    k_free: np.ndarray,
    m_free: np.ndarray,
    k_full: np.ndarray,
    m_full: np.ndarray,
    free: np.ndarray,
    omega_sq: float,
    fem_lambda: float,
    q_fem: np.ndarray,
    translational_mask: np.ndarray,
    rotational_mask: np.ndarray,
) -> list[dict[str, float | str]]:
    element_context = element_diagnostic_context(params, mu=float(mu), beta_deg=float(beta_deg))
    rows: list[dict[str, float | str]] = []

    def append_variant(variant_name: str, scale_source: str, alpha_left: float, alpha_right: float) -> None:
        if not (np.isfinite(float(alpha_left)) and np.isfinite(float(alpha_right))):
            return
        rows.append(
            evaluate_axial_scale_variant(
                fields_base=fields_base,
                beta_deg=float(beta_deg),
                params=params,
                k_free=k_free,
                m_free=m_free,
                k_full=k_full,
                m_full=m_full,
                free=free,
                omega_sq=float(omega_sq),
                fem_lambda=float(fem_lambda),
                q_fem=q_fem,
                translational_mask=translational_mask,
                rotational_mask=rotational_mask,
                element_context=element_context,
                variant_name=variant_name,
                scale_source=scale_source,
                alpha_left=float(alpha_left),
                alpha_right=float(alpha_right),
            )
        )

    candidates = axial_scale_candidate_values(Lambda=float(Lambda), epsilon=float(epsilon), mu=float(mu))
    for left_idx, (left_label, alpha_left) in enumerate(candidates):
        for right_idx, (right_label, alpha_right) in enumerate(candidates):
            variant_name = (
                "baseline_alpha_1_1"
                if np.isclose(alpha_left, 1.0, rtol=1e-12, atol=1e-15)
                and np.isclose(alpha_right, 1.0, rtol=1e-12, atol=1e-15)
                else f"candidate_L{left_idx:02d}_R{right_idx:02d}"
            )
            append_variant(
                variant_name,
                f"candidate_grid:left={left_label};right={right_label}",
                alpha_left,
                alpha_right,
            )

    alpha_left_l2, alpha_right_l2 = l2_optimal_axial_alphas(
        fields_base,
        beta_deg=float(beta_deg),
        q_fem=np.asarray(q_fem, dtype=float),
        free=np.asarray(free, dtype=int),
    )
    append_variant(
        "l2_optimal_local_axial",
        "l2_optimal_after_baseline_global_alignment",
        alpha_left_l2,
        alpha_right_l2,
    )

    log_values = np.logspace(-6.0, 6.0, 49)
    for left_idx, alpha_left in enumerate(log_values):
        for right_idx, alpha_right in enumerate(log_values):
            append_variant(
                f"log_grid_L{left_idx:02d}_R{right_idx:02d}",
                "coarse_log_grid_1e-6_to_1e6_49x49",
                float(alpha_left),
                float(alpha_right),
            )

    baseline = next(
        row
        for row in rows
        if np.isclose(float(row["alpha_left"]), 1.0, rtol=1e-12, atol=1e-15)
        and np.isclose(float(row["alpha_right"]), 1.0, rtol=1e-12, atol=1e-15)
    )
    baseline_residual = float(baseline["relative_residual"])
    baseline_rayleigh_diff = float(baseline["rayleigh_lambda_relative_diff"])
    for row in rows:
        row["conclusion_flag"] = axial_scale_conclusion(
            variant_name=str(row["variant_name"]),
            scale_source=str(row["scale_source"]),
            alpha_left=float(row["alpha_left"]),
            alpha_right=float(row["alpha_right"]),
            residual=float(row["relative_residual"]),
            baseline_residual=baseline_residual,
            rayleigh_diff=float(row["rayleigh_lambda_relative_diff"]),
            baseline_rayleigh_diff=baseline_rayleigh_diff,
        )
    rows.sort(key=lambda row: float(row["relative_residual"]))
    return rows


def write_axial_scale_scan_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=AXIAL_SCALE_SCAN_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def best_axial_scale_by_rayleigh(rows: Sequence[dict[str, float | str]]) -> dict[str, float | str] | None:
    finite_rows = [
        row
        for row in rows
        if np.isfinite(float(row["rayleigh_lambda_relative_diff"]))
    ]
    if not finite_rows:
        return None
    return min(finite_rows, key=lambda row: float(row["rayleigh_lambda_relative_diff"]))


def best_axial_scale_by_mass_mac(rows: Sequence[dict[str, float | str]]) -> dict[str, float | str] | None:
    finite_rows = [row for row in rows if np.isfinite(float(row["mass_weighted_mac"]))]
    if not finite_rows:
        return None
    return max(finite_rows, key=lambda row: float(row["mass_weighted_mac"]))


def axial_scale_overall_flag(rows: Sequence[dict[str, float | str]]) -> str:
    if not rows:
        return ""
    baseline = next((row for row in rows if str(row["variant_name"]) == "baseline_alpha_1_1"), None)
    if baseline is None:
        return "needs_review"
    simple_possible = any(
        str(row["conclusion_flag"]) == "possible_axial_amplitude_scaling_mismatch"
        for row in rows
    )
    if simple_possible:
        return "possible_axial_amplitude_scaling_mismatch_detected"
    best_rayleigh = best_axial_scale_by_rayleigh(rows)
    if best_rayleigh is not None:
        both_zero = (
            abs(float(best_rayleigh["alpha_left"])) <= NEAR_ZERO_NORM
            and abs(float(best_rayleigh["alpha_right"])) <= NEAR_ZERO_NORM
        )
        if both_zero and str(best_rayleigh["conclusion_flag"]) == "zero_axial_improves_rayleigh":
            return "analytic_axial_components_likely_over_penalized_no_physical_scale_found_yet"
    if any(str(row["conclusion_flag"]) == "needs_review" for row in rows):
        return "scale_variant_needs_review"
    return "axial_scale_mismatch_not_supported_likely_operator_model_mismatch"


def print_axial_scale_scan_top(rows: Sequence[dict[str, float | str]], *, limit: int = 6) -> None:
    if not rows:
        print("axial scale scan: no rows")
        return
    baseline = next((row for row in rows if str(row["variant_name"]) == "baseline_alpha_1_1"), None)
    if baseline is not None:
        print(
            "baseline axial scale: "
            f"residual={float(baseline['relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(baseline['rayleigh_lambda']):.6e}"
        )
    print("Top axial scale variants by analytic-in-FEM relative residual:")
    for row in list(rows)[:limit]:
        print(
            f"  residual={float(row['relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(row['rayleigh_lambda']):.6e}, "
            f"alpha=({float(row['alpha_left']):.6g}, {float(row['alpha_right']):.6g}), "
            f"{row['variant_name']} ({row['conclusion_flag']})"
        )
    rayleigh_sorted = sorted(
        rows,
        key=lambda row: float(row["rayleigh_lambda_relative_diff"])
        if np.isfinite(float(row["rayleigh_lambda_relative_diff"]))
        else np.inf,
    )
    print("Top axial scale variants by Rayleigh Lambda closeness:")
    for row in rayleigh_sorted[:limit]:
        print(
            f"  rel_Lambda_diff={float(row['rayleigh_lambda_relative_diff']):.6e}, "
            f"Rayleigh Lambda={float(row['rayleigh_lambda']):.6e}, "
            f"residual={float(row['relative_residual']):.6e}, "
            f"alpha=({float(row['alpha_left']):.6g}, {float(row['alpha_right']):.6g}), "
            f"{row['variant_name']} ({row['conclusion_flag']})"
        )
    mass_sorted = sorted(
        rows,
        key=lambda row: float(row["mass_weighted_mac"])
        if np.isfinite(float(row["mass_weighted_mac"]))
        else -np.inf,
        reverse=True,
    )
    print("Top axial scale variants by mass-weighted MAC:")
    for row in mass_sorted[:limit]:
        print(
            f"  mass_MAC={float(row['mass_weighted_mac']):.6e}, "
            f"residual={float(row['relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(row['rayleigh_lambda']):.6e}, "
            f"alpha=({float(row['alpha_left']):.6g}, {float(row['alpha_right']):.6g}), "
            f"{row['variant_name']} ({row['conclusion_flag']})"
        )
    flag = axial_scale_overall_flag(rows)
    if flag == "possible_axial_amplitude_scaling_mismatch_detected":
        print("possible axial amplitude scaling mismatch detected")
    elif flag == "analytic_axial_components_likely_over_penalized_no_physical_scale_found_yet":
        print("analytic axial components likely over-penalized in FEM embedding, but no physically justified scale found yet")
    elif flag == "axial_scale_mismatch_not_supported_likely_operator_model_mismatch":
        print("axial scale mismatch not supported; likely operator/model mismatch")
    else:
        print(f"axial scale scan conclusion: {flag}")


def fem_transform_inference_text() -> str:
    return "K_global=T.T K_local T implies q_local=T q_global; therefore FEM-consistent global_from_local is T.T."


def right_transform_matrix(beta_deg: float, transform_kind: str) -> np.ndarray:
    beta_rad = float(np.deg2rad(beta_deg))
    cb = float(np.cos(beta_rad))
    sb = float(np.sin(beta_rad))
    t = np.array([[cb, -sb], [sb, cb]], dtype=float)
    if transform_kind == "T":
        return t
    if transform_kind == "T_transpose":
        return t.T
    raise ValueError(f"Unknown right-arm transform kind: {transform_kind!r}")


def right_transform_variants() -> tuple[tuple[str, str, str], ...]:
    return (
        (
            "baseline_current",
            "current embedding: q_global = R(+beta) q_local = T q_local",
            "T",
        ),
        (
            "right_global_from_local_T",
            "q_global = T q_local with T = rotation_matrix_6x6(+beta)",
            "T",
        ),
        (
            "right_global_from_local_T_transpose",
            "q_global = T.T q_local",
            "T_transpose",
        ),
        (
            "right_global_from_local_R_beta",
            "q_global = R(+beta) q_local",
            "T",
        ),
        (
            "right_global_from_local_R_minus_beta",
            "q_global = R(-beta) q_local = T.T q_local",
            "T_transpose",
        ),
        (
            "right_global_from_local_inverse_of_fem_T",
            "q_global = inv(T) q_local = T.T q_local because T is orthogonal",
            "T_transpose",
        ),
        (
            "right_global_from_local_same_as_fem_energy_convention",
            "FEM assembly Kglob=T.T Kloc T gives q_local=T q_global, so q_global=T.T q_local",
            "T_transpose",
        ),
    )


def build_analytic_global_q_with_right_transform(
    fields: dict[str, np.ndarray],
    *,
    beta_deg: float,
    transform_kind: str,
) -> tuple[np.ndarray, dict[str, float]]:
    n = fem.N_ELEM
    n_nodes = 2 * n + 1
    ux = np.zeros(n_nodes, dtype=float)
    uy = np.zeros(n_nodes, dtype=float)
    theta = np.zeros(n_nodes, dtype=float)

    ux[: n + 1] = np.asarray(fields["u_left"], dtype=float)
    uy[: n + 1] = np.asarray(fields["w_left"], dtype=float)
    theta[: n + 1] = np.asarray(fields["theta_left"], dtype=float)

    u_right = np.asarray(fields["u_right"], dtype=float)
    w_right = np.asarray(fields["w_right"], dtype=float)
    theta_right = np.asarray(fields["theta_right"], dtype=float)
    transform = right_transform_matrix(float(beta_deg), transform_kind)
    right_global = transform @ np.vstack([u_right, w_right])
    right_ux = np.asarray(right_global[0], dtype=float)
    right_uy = np.asarray(right_global[1], dtype=float)
    ux[n + 1 :] = right_ux[1:]
    uy[n + 1 :] = right_uy[1:]
    theta[n + 1 :] = theta_right[1:]

    joint_translation_mismatch = max(
        abs(float(ux[n]) - float(right_ux[0])),
        abs(float(uy[n]) - float(right_uy[0])),
    )
    joint_theta_mismatch = abs(float(theta[n]) - float(theta_right[0]))
    q = np.empty(3 * n_nodes, dtype=float)
    q[0::3] = ux
    q[1::3] = uy
    q[2::3] = theta
    return q, {
        "joint_translation_embedding_mismatch": float(joint_translation_mismatch),
        "joint_theta_embedding_mismatch": float(joint_theta_mismatch),
    }


def external_clamp_max_abs_from_q(q_full: np.ndarray) -> float:
    n = fem.N_ELEM
    q = np.asarray(q_full, dtype=float)
    left_clamp = q[0:3]
    right_clamp = q[3 * (2 * n) : 3 * (2 * n) + 3]
    return float(max(np.max(np.abs(left_clamp)), np.max(np.abs(right_clamp))))


def right_transform_conclusion(
    *,
    variant_name: str,
    residual: float,
    baseline_residual: float,
    rayleigh_diff: float,
    baseline_rayleigh_diff: float,
    joint_translation_mismatch: float,
) -> str:
    if variant_name == "baseline_current":
        return "baseline"
    residual_dramatic = np.isfinite(residual) and float(residual) < 0.1 * float(baseline_residual)
    rayleigh_dramatic = (
        np.isfinite(rayleigh_diff)
        and np.isfinite(baseline_rayleigh_diff)
        and float(rayleigh_diff) < 0.1 * float(baseline_rayleigh_diff)
    )
    joint_ok = np.isfinite(joint_translation_mismatch) and float(joint_translation_mismatch) <= 1e-7
    if (residual_dramatic or rayleigh_dramatic) and joint_ok:
        return "possible_right_arm_transform_convention_mismatch"
    if residual_dramatic or rayleigh_dramatic:
        return "improves_but_breaks_joint_continuity"
    residual_review = np.isfinite(residual) and float(residual) < float(baseline_residual) * 0.9
    rayleigh_review = (
        np.isfinite(rayleigh_diff)
        and np.isfinite(baseline_rayleigh_diff)
        and float(rayleigh_diff) < float(baseline_rayleigh_diff) * 0.9
    )
    if (residual_review or rayleigh_review) and joint_ok:
        return "needs_review"
    if residual_review or rayleigh_review:
        return "rayleigh_or_residual_closer_but_breaks_joint_continuity"
    return "no_improvement"


def right_transform_convention_scan_rows(
    *,
    fields_base: dict[str, np.ndarray],
    beta_deg: float,
    params: BeamParams,
    k_free: np.ndarray,
    m_free: np.ndarray,
    k_full: np.ndarray,
    m_full: np.ndarray,
    free: np.ndarray,
    omega_sq: float,
    fem_lambda: float,
    q_fem: np.ndarray,
    translational_mask: np.ndarray,
    rotational_mask: np.ndarray,
    element_context: dict[str, np.ndarray],
) -> list[dict[str, float | str]]:
    del params
    rows: list[dict[str, float | str]] = []
    q_fem_full = np.asarray(q_fem, dtype=float)
    free_idx = np.asarray(free, dtype=int)
    q_fem_free = q_fem_full[free_idx]
    for variant_name, description, transform_kind in right_transform_variants():
        q_analytic, embedding = build_analytic_global_q_with_right_transform(
            fields_base,
            beta_deg=float(beta_deg),
            transform_kind=transform_kind,
        )
        q_analytic_free = q_analytic[free_idx]
        alignment_scalar = optimal_scalar(q_analytic_free, q_fem_free)
        q_aligned = alignment_scalar * q_analytic
        residual = residual_metrics(k_free, m_free, q_analytic, free_idx, omega_sq=float(omega_sq))
        rayleigh = rayleigh_diagnostics(
            k_free,
            m_free,
            q_analytic,
            free_idx,
            fem_omega_sq=float(omega_sq),
            fem_lambda=float(fem_lambda),
        )
        group_rows = residual_group_rows(k_full, m_full, q_aligned, free_idx, omega_sq=float(omega_sq))
        energy = element_energy_sums_from_context(
            element_context,
            q_analytic=q_aligned,
            q_fem=q_fem_full,
            omega_sq=float(omega_sq),
        )
        rows.append(
            {
                "variant_name": str(variant_name),
                "transform_description": str(description),
                "relative_residual": float(residual["relative_residual"]),
                "rayleigh_lambda": float(rayleigh["analytic_rayleigh_Lambda"]),
                "rayleigh_lambda_relative_diff": float(rayleigh["relative_rayleigh_Lambda_diff"]),
                "global_dof_mac_euclidean": euclidean_mac(q_analytic_free, q_fem_free),
                "mass_weighted_mac": mass_weighted_mac(q_analytic_free, q_fem_free, m_free),
                "translational_relative_l2": subset_relative_l2(q_aligned, q_fem_full, translational_mask),
                "rotational_relative_l2": subset_relative_l2(q_aligned, q_fem_full, rotational_mask),
                "right_arm_translational_residual_share": residual_group_share_by_name(
                    group_rows,
                    "right_arm_translational",
                ),
                "axial_to_bending_residual_ratio": float(energy["axial_to_bending_residual_ratio"]),
                "joint_translation_mismatch": float(embedding["joint_translation_embedding_mismatch"]),
                "joint_rotation_mismatch": float(embedding["joint_theta_embedding_mismatch"]),
                "external_clamp_max_abs": external_clamp_max_abs_from_q(q_analytic),
                "conclusion_flag": "",
            }
        )
    baseline = next(row for row in rows if row["variant_name"] == "baseline_current")
    baseline_residual = float(baseline["relative_residual"])
    baseline_rayleigh_diff = float(baseline["rayleigh_lambda_relative_diff"])
    for row in rows:
        row["conclusion_flag"] = right_transform_conclusion(
            variant_name=str(row["variant_name"]),
            residual=float(row["relative_residual"]),
            baseline_residual=baseline_residual,
            rayleigh_diff=float(row["rayleigh_lambda_relative_diff"]),
            baseline_rayleigh_diff=baseline_rayleigh_diff,
            joint_translation_mismatch=float(row["joint_translation_mismatch"]),
        )
    rows.sort(key=lambda row: float(row["relative_residual"]))
    return rows


def right_transform_sanity_rows(beta_deg: float) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for variant_name, description, transform_kind in right_transform_variants():
        transform = right_transform_matrix(float(beta_deg), transform_kind)
        axial = transform @ np.array([1.0, 0.0], dtype=float)
        transverse = transform @ np.array([0.0, 1.0], dtype=float)
        rows.append(
            {
                "variant_name": str(variant_name),
                "transform_description": str(description),
                "fem_assembly_inference": fem_transform_inference_text(),
                "local_unit_axial_global_ux": float(axial[0]),
                "local_unit_axial_global_uy": float(axial[1]),
                "local_unit_transverse_global_ux": float(transverse[0]),
                "local_unit_transverse_global_uy": float(transverse[1]),
            }
        )
    return rows


def write_right_transform_convention_scan_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RIGHT_TRANSFORM_CONVENTION_SCAN_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_right_transform_sanity_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RIGHT_TRANSFORM_SANITY_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def best_right_transform_by_rayleigh(rows: Sequence[dict[str, float | str]]) -> dict[str, float | str] | None:
    finite_rows = [row for row in rows if np.isfinite(float(row["rayleigh_lambda_relative_diff"]))]
    if not finite_rows:
        return None
    return min(finite_rows, key=lambda row: float(row["rayleigh_lambda_relative_diff"]))


def best_right_transform_by_mass_mac(rows: Sequence[dict[str, float | str]]) -> dict[str, float | str] | None:
    finite_rows = [row for row in rows if np.isfinite(float(row["mass_weighted_mac"]))]
    if not finite_rows:
        return None
    return max(finite_rows, key=lambda row: float(row["mass_weighted_mac"]))


def right_transform_overall_flag(rows: Sequence[dict[str, float | str]]) -> str:
    if not rows:
        return ""
    if any(
        str(row["conclusion_flag"]) == "possible_right_arm_transform_convention_mismatch"
        for row in rows
    ):
        return "possible_right_arm_transform_convention_mismatch_detected"
    if any(str(row["conclusion_flag"]) == "needs_review" for row in rows):
        return "right_transform_variant_needs_review"
    return "simple_right_arm_transform_convention_mismatch_not_supported"


def print_right_transform_convention_scan_top(
    rows: Sequence[dict[str, float | str]],
    sanity_rows: Sequence[dict[str, float | str]],
    *,
    limit: int = 6,
) -> None:
    if not rows:
        print("right transform convention scan: no rows")
        return
    print(f"FEM transform inference: {fem_transform_inference_text()}")
    print("Right transform sanity, local unit vectors -> global translations:")
    for row in sanity_rows:
        print(
            f"  {row['variant_name']}: "
            f"axial=({float(row['local_unit_axial_global_ux']):+.6f}, "
            f"{float(row['local_unit_axial_global_uy']):+.6f}), "
            f"transverse=({float(row['local_unit_transverse_global_ux']):+.6f}, "
            f"{float(row['local_unit_transverse_global_uy']):+.6f})"
        )
    baseline = next((row for row in rows if str(row["variant_name"]) == "baseline_current"), None)
    if baseline is not None:
        print(
            "baseline right transform: "
            f"residual={float(baseline['relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(baseline['rayleigh_lambda']):.6e}"
        )
    print("Top right transform variants by analytic-in-FEM relative residual:")
    for row in list(rows)[:limit]:
        print(
            f"  residual={float(row['relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(row['rayleigh_lambda']):.6e}, "
            f"mass_MAC={float(row['mass_weighted_mac']):.6e}, "
            f"joint_mismatch={float(row['joint_translation_mismatch']):.6e}, "
            f"{row['variant_name']} ({row['conclusion_flag']})"
        )
    rayleigh_sorted = sorted(
        rows,
        key=lambda row: float(row["rayleigh_lambda_relative_diff"])
        if np.isfinite(float(row["rayleigh_lambda_relative_diff"]))
        else np.inf,
    )
    print("Top right transform variants by Rayleigh Lambda closeness:")
    for row in rayleigh_sorted[:limit]:
        print(
            f"  rel_Lambda_diff={float(row['rayleigh_lambda_relative_diff']):.6e}, "
            f"Rayleigh Lambda={float(row['rayleigh_lambda']):.6e}, "
            f"residual={float(row['relative_residual']):.6e}, "
            f"{row['variant_name']} ({row['conclusion_flag']})"
        )
    mass_sorted = sorted(
        rows,
        key=lambda row: float(row["mass_weighted_mac"])
        if np.isfinite(float(row["mass_weighted_mac"]))
        else -np.inf,
        reverse=True,
    )
    print("Top right transform variants by mass-weighted MAC:")
    for row in mass_sorted[:limit]:
        print(
            f"  mass_MAC={float(row['mass_weighted_mac']):.6e}, "
            f"residual={float(row['relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(row['rayleigh_lambda']):.6e}, "
            f"{row['variant_name']} ({row['conclusion_flag']})"
        )
    flag = right_transform_overall_flag(rows)
    if flag == "possible_right_arm_transform_convention_mismatch_detected":
        print("possible right-arm transform convention mismatch detected")
    elif flag == "simple_right_arm_transform_convention_mismatch_not_supported":
        print("simple right-arm transform convention mismatch not supported")
    else:
        print(f"right transform scan conclusion: {flag}")


def write_modal_projection_csv(path: Path, rows: list[dict[str, float | int]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MODAL_PROJECTION_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_residual_groups_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=RESIDUAL_GROUP_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_element_residual_energy_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ELEMENT_RESIDUAL_ENERGY_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def full_basis_column_audit_rows(
    matrix: np.ndarray,
    *,
    Lambda: float,
    beta_rad: float,
    mu: float,
    epsilon: float,
) -> list[dict[str, float | str]]:
    left_mapping = canonical_bending_basis_mapping()
    right_mapping = canonical_bending_basis_mapping()
    rows: list[dict[str, float | str]] = []
    for right_z_sign in ("-", "+"):
        for col_idx, basis_column in enumerate(BASIS_COLUMN_LABELS):
            coeff = np.zeros(6, dtype=float)
            coeff[col_idx] = 1.0
            matrix_values = np.asarray(matrix, dtype=float) @ coeff
            endpoint_values = full_basis_endpoint_residual(
                float(Lambda),
                beta_rad=float(beta_rad),
                mu=float(mu),
                epsilon=float(epsilon),
                coeff=coeff,
                left_mapping=left_mapping,
                right_mapping=right_mapping,
                right_z_sign=right_z_sign,
            )
            for row_idx, row_label in enumerate(ROW_LABELS):
                matrix_value = float(matrix_values[row_idx])
                endpoint_value = float(endpoint_values[row_idx])
                diff = abs(matrix_value - endpoint_value)
                denominator = max(abs(matrix_value), abs(endpoint_value), NEAR_ZERO_NORM)
                rows.append(
                    {
                        "basis_column": str(basis_column),
                        "row_label": str(row_label),
                        "matrix_column_value": matrix_value,
                        "independent_full_basis_value": endpoint_value,
                        "abs_difference": float(diff),
                        "relative_difference": float(diff / denominator),
                        "variant_name": "canonical_full_basis",
                        "right_z_sign": right_z_sign,
                    }
                )
    return rows


def pairing_scan_conclusion(
    *,
    variant_name: str,
    residual: float,
    canonical_residual: float,
    matrix_mismatch: float,
    external_clamp_residual: float,
) -> str:
    consistency_ok = matrix_mismatch <= 1e-7 and external_clamp_residual <= 1e-9
    improves = float(residual) < float(canonical_residual) * (1.0 - 1e-6)
    if variant_name == "canonical_full_basis_right_z-" and consistency_ok:
        return "canonical_consistent"
    if improves and consistency_ok:
        return "improves_residual_and_consistent"
    if improves and not consistency_ok:
        return "improves_residual_but_breaks_matrix_consistency"
    if not improves:
        return "no_improvement"
    return "needs_review"


def bending_basis_pairing_scan_rows(
    *,
    Lambda: float,
    beta_deg: float,
    beta_rad: float,
    mu: float,
    epsilon: float,
    coeff: np.ndarray,
    s_norm: np.ndarray,
    l_base: float,
    matrix: np.ndarray,
    k_free: np.ndarray,
    m_free: np.ndarray,
    free: np.ndarray,
    omega_sq: float,
    q_fem: np.ndarray,
    translational_mask: np.ndarray,
    rotational_mask: np.ndarray,
) -> list[dict[str, float | str]]:
    matrix_residual = np.asarray(matrix, dtype=float) @ np.asarray(coeff, dtype=float)
    matrix_residual_norm = float(np.linalg.norm(matrix_residual))
    mappings = bending_basis_mappings()
    q_fem_full = np.asarray(q_fem, dtype=float)
    q_fem_free = q_fem_full[np.asarray(free, dtype=int)]
    rows: list[dict[str, float | str]] = []
    for left_mapping in mappings:
        for right_mapping in mappings:
            for right_z_sign in ("-", "+"):
                variant_name = f"left_{left_mapping.name}__right_{right_mapping.name}__right_z{right_z_sign}"
                if (
                    left_mapping.name == "identity_A+_B+"
                    and right_mapping.name == "identity_A+_B+"
                    and right_z_sign == "-"
                ):
                    variant_name = "canonical_full_basis_right_z-"
                fields = full_basis_fields(
                    float(Lambda),
                    mu=float(mu),
                    epsilon=float(epsilon),
                    coeff=np.asarray(coeff, dtype=float),
                    s_norm=np.asarray(s_norm, dtype=float),
                    l_base=float(l_base),
                    left_mapping=left_mapping,
                    right_mapping=right_mapping,
                    right_z_sign=right_z_sign,
                )
                q_analytic, _ = build_analytic_global_q(fields, beta_deg=float(beta_deg))
                q_analytic_free = q_analytic[np.asarray(free, dtype=int)]
                alignment_scalar = optimal_scalar(q_analytic_free, q_fem_free)
                q_aligned = alignment_scalar * q_analytic
                endpoint_residual = full_basis_endpoint_residual(
                    float(Lambda),
                    beta_rad=float(beta_rad),
                    mu=float(mu),
                    epsilon=float(epsilon),
                    coeff=np.asarray(coeff, dtype=float),
                    left_mapping=left_mapping,
                    right_mapping=right_mapping,
                    right_z_sign=right_z_sign,
                )
                residual = residual_metrics(k_free, m_free, q_analytic, free, omega_sq=float(omega_sq))
                external_clamp = full_basis_external_clamp_residual(
                    np.asarray(coeff, dtype=float),
                    left_mapping=left_mapping,
                    right_mapping=right_mapping,
                )
                rows.append(
                    {
                        "variant_name": variant_name,
                        "left_mapping": left_mapping.name,
                        "right_mapping": right_mapping.name,
                        "right_z_sign": right_z_sign,
                        "external_clamp_residual": float(external_clamp),
                        "matrix_residual_norm": matrix_residual_norm,
                        "independent_endpoint_residual_norm": float(np.linalg.norm(endpoint_residual)),
                        "matrix_vs_independent_endpoint_mismatch": float(
                            np.linalg.norm(matrix_residual - endpoint_residual)
                        ),
                        "analytic_in_fem_relative_residual": float(residual["relative_residual"]),
                        "global_dof_mac_euclidean": euclidean_mac(q_analytic_free, q_fem_free),
                        "mass_weighted_mac": mass_weighted_mac(q_analytic_free, q_fem_free, m_free),
                        "translational_relative_l2": subset_relative_l2(q_aligned, q_fem_full, translational_mask),
                        "rotational_relative_l2": subset_relative_l2(q_aligned, q_fem_full, rotational_mask),
                        "conclusion_flag": "",
                    }
                )
    canonical_rows = [
        row for row in rows if str(row["variant_name"]) == "canonical_full_basis_right_z-"
    ]
    if canonical_rows:
        canonical_residual = float(canonical_rows[0]["analytic_in_fem_relative_residual"])
    else:
        canonical_residual = float(
            next(
                row["analytic_in_fem_relative_residual"]
                for row in rows
                if row["left_mapping"] == "identity_A+_B+"
                and row["right_mapping"] == "identity_A+_B+"
                and row["right_z_sign"] == "-"
            )
        )
    for row in rows:
        row["conclusion_flag"] = pairing_scan_conclusion(
            variant_name=str(row["variant_name"]),
            residual=float(row["analytic_in_fem_relative_residual"]),
            canonical_residual=canonical_residual,
            matrix_mismatch=float(row["matrix_vs_independent_endpoint_mismatch"]),
            external_clamp_residual=float(row["external_clamp_residual"]),
        )
    rows.sort(key=lambda row: float(row["analytic_in_fem_relative_residual"]))
    return rows


def canonical_pairing_row(rows: Sequence[dict[str, float | str]]) -> dict[str, float | str] | None:
    for row in rows:
        if (
            row["left_mapping"] == "identity_A+_B+"
            and row["right_mapping"] == "identity_A+_B+"
            and row["right_z_sign"] == "-"
        ):
            return row
    return None


def write_bending_basis_pairing_scan_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=BENDING_BASIS_PAIRING_SCAN_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_full_basis_column_audit_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FULL_BASIS_COLUMN_AUDIT_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_frequency_scaling_audit_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FREQUENCY_SCALING_AUDIT_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def print_bending_basis_pairing_scan_top(rows: Sequence[dict[str, float | str]], *, limit: int = 6) -> None:
    if not rows:
        print("bending basis pairing scan: no rows")
        return
    canonical = canonical_pairing_row(rows)
    best = rows[0]
    if canonical is not None:
        print(
            "canonical full-basis residual: "
            f"{float(canonical['analytic_in_fem_relative_residual']):.6e}, "
            f"matrix mismatch={float(canonical['matrix_vs_independent_endpoint_mismatch']):.6e}"
        )
    print(
        "best bending basis variant: "
        f"{best['variant_name']} residual={float(best['analytic_in_fem_relative_residual']):.6e}, "
        f"matrix mismatch={float(best['matrix_vs_independent_endpoint_mismatch']):.6e}, "
        f"flag={best['conclusion_flag']}"
    )
    print("Top bending basis variants by analytic-in-FEM relative residual:")
    for row in list(rows)[:limit]:
        print(
            f"  residual={float(row['analytic_in_fem_relative_residual']):.6e}, "
            f"matrix_mismatch={float(row['matrix_vs_independent_endpoint_mismatch']):.6e}, "
            f"{row['variant_name']} ({row['conclusion_flag']})"
        )
    best_flag = str(best["conclusion_flag"])
    if best_flag in {"improves_residual_and_consistent", "improves_residual_but_breaks_matrix_consistency"}:
        print("possible coefficient-basis pairing mismatch detected")
    else:
        print("bending basis pairing scan: no supported A/B pairing mismatch detected")


def print_frequency_scaling_audit(rows: Sequence[dict[str, float | str]]) -> None:
    if not rows:
        print("frequency scaling audit: no rows")
        return
    first = rows[0]
    print("Frequency scaling audit")
    print(f"analytic Lambda = {float(first['analytic_lambda']):.12g}")
    print(f"analytic Lambda^2 = {float(first['analytic_lambda_squared']):.12g}")
    print(f"analytic Lambda^4 = {float(first['analytic_lambda_fourth']):.12g}")
    print(f"FEM omega = {float(first['fem_selected_omega']):.12g}")
    print(f"FEM omega^2 = {float(first['fem_selected_omega_sq']):.12g}")
    print(f"FEM Lambda from omega = {float(first['fem_selected_lambda_from_omega']):.12g}")
    print(f"FEM Lambda from omega_sq = {float(first['fem_selected_lambda_from_omega_sq']):.12g}")
    print(
        "ratios: "
        f"Lambda^4/FEM omega^2={float(first['lambda_fourth_to_fem_omega_sq_ratio']):.12g}, "
        f"Lambda^2/FEM omega={float(first['lambda_squared_to_fem_omega_ratio']):.12g}"
    )
    print("Residuals with spectral factor inserted as Kq - factor*Mq:")
    print("  vector_type       factor                         factor_value      relative_residual")
    for row in rows:
        print(
            f"  {str(row['vector_type']):16s} "
            f"{str(row['spectral_factor_name']):29s} "
            f"{float(row['spectral_factor_value']):.8e} "
            f"{float(row['relative_residual']):.8e}"
        )


def node_coordinates(params: BeamParams, *, mu: float, beta_deg: float) -> tuple[np.ndarray, np.ndarray]:
    n = fem.N_ELEM
    ell = float(params.L_base)
    length_left = ell * (1.0 - float(mu))
    length_right = ell * (1.0 + float(mu))
    beta_rad = float(np.deg2rad(beta_deg))
    x = np.zeros(2 * n + 1, dtype=float)
    y = np.zeros(2 * n + 1, dtype=float)
    x[: n + 1] = np.linspace(0.0, length_left, n + 1)
    s_right = np.linspace(0.0, length_right, n + 1)
    x[n:] = length_left + s_right * np.cos(beta_rad)
    y[n:] = s_right * np.sin(beta_rad)
    return x, y


def subset_relative_l2(q_analytic_aligned: np.ndarray, q_fem: np.ndarray, mask: np.ndarray) -> float:
    return safe_relative_l2(
        np.asarray(q_analytic_aligned, dtype=float)[mask] - np.asarray(q_fem, dtype=float)[mask],
        np.asarray(q_fem, dtype=float)[mask],
    )


def conclusion_flag(row: dict[str, float | int | str]) -> str:
    if str(row["branch_tracking_warning_flag"]) != "ok":
        return "branch_warning"
    relative_residual = float(row["relative_residual"])
    global_mac = float(row["global_dof_mac_euclidean"])
    global_relative = float(row["global_dof_relative_l2"])
    rel_lambda = float(row["rel_lambda_diff"])
    if relative_residual > RESIDUAL_BAD_THRESHOLD:
        return "residual_bad"
    if relative_residual <= RESIDUAL_GOOD_THRESHOLD and (
        global_mac < GLOBAL_MAC_BAD_THRESHOLD or global_relative > GLOBAL_REL_BAD_THRESHOLD
    ):
        return "residual_good_component_bad"
    if (
        rel_lambda <= FREQUENCY_GOOD_THRESHOLD
        and relative_residual <= RESIDUAL_GOOD_THRESHOLD
        and global_mac >= GLOBAL_MAC_GOOD_THRESHOLD
    ):
        return "good_match"
    return "needs_review"


def write_summary_csv(path: Path, row: dict[str, float | int | str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerow(row)


def build_node_debug_rows(
    *,
    q_analytic_aligned: np.ndarray,
    q_fem: np.ndarray,
    params: BeamParams,
    mu: float,
    beta_deg: float,
) -> list[dict[str, float | int | str]]:
    n = fem.N_ELEM
    x0, y0 = node_coordinates(params, mu=mu, beta_deg=beta_deg)
    rows: list[dict[str, float | int | str]] = []

    def append_row(node_index: int, arm: str, s_norm: float) -> None:
        base = 3 * int(node_index)
        analytic = np.asarray(q_analytic_aligned, dtype=float)[base : base + 3]
        fem_values = np.asarray(q_fem, dtype=float)[base : base + 3]
        rows.append(
            {
                "node_index": int(node_index),
                "arm": arm,
                "s_norm": float(s_norm),
                "x0": float(x0[node_index]),
                "y0": float(y0[node_index]),
                "analytic_ux": float(analytic[0]),
                "analytic_uy": float(analytic[1]),
                "analytic_theta": float(analytic[2]),
                "fem_ux": float(fem_values[0]),
                "fem_uy": float(fem_values[1]),
                "fem_theta": float(fem_values[2]),
                "diff_ux": float(analytic[0] - fem_values[0]),
                "diff_uy": float(analytic[1] - fem_values[1]),
                "diff_theta": float(analytic[2] - fem_values[2]),
            }
        )

    for idx in range(n + 1):
        append_row(idx, "left", idx / n)
    for local_idx in range(n + 1):
        append_row(n + local_idx, "right", local_idx / n)
    return rows


def write_node_debug_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=NODE_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def evaluate_theta_convention_scan(
    *,
    analytic_lambda: float,
    mu: float,
    epsilon: float,
    coeff: np.ndarray,
    s_norm: np.ndarray,
    l_base: float,
    beta_deg: float,
    k_free: np.ndarray,
    m_free: np.ndarray,
    free: np.ndarray,
    q_fem: np.ndarray,
    omega_sq: float,
    translational_mask: np.ndarray,
    rotational_mask: np.ndarray,
    left_mask: np.ndarray,
    right_mask: np.ndarray,
    joint_mask: np.ndarray,
) -> list[dict[str, float | int | str]]:
    components = reconstruct_analytic_components(
        float(analytic_lambda),
        mu_value=float(mu),
        epsilon=float(epsilon),
        coeff=np.asarray(coeff, dtype=float),
        s_norm=np.asarray(s_norm, dtype=float),
        right_coordinate="joint-to-external",
    )
    left_candidates, right_candidates = analytic_theta_candidate_sets(
        float(analytic_lambda),
        mu=float(mu),
        coeff=np.asarray(coeff, dtype=float),
        s_norm=np.asarray(s_norm, dtype=float),
        l_base=float(l_base),
    )
    q_fem_full = np.asarray(q_fem, dtype=float)
    q_fem_free = q_fem_full[np.asarray(free, dtype=int)]
    rows: list[dict[str, float | int | str]] = []
    for left_label in theta_candidate_labels():
        for right_label in theta_candidate_labels():
            fields = {
                "u_left": np.asarray(components["u_left"], dtype=float),
                "w_left": np.asarray(components["w_left"], dtype=float),
                "theta_left": np.asarray(left_candidates[left_label], dtype=float),
                "u_right": np.asarray(components["u_right"], dtype=float),
                "w_right": np.asarray(components["w_right"], dtype=float),
                "theta_right": np.asarray(right_candidates[right_label], dtype=float),
            }
            q_analytic, _ = build_analytic_global_q(fields, beta_deg=float(beta_deg))
            q_analytic_free = q_analytic[np.asarray(free, dtype=int)]
            alignment_scalar = optimal_scalar(q_analytic_free, q_fem_free)
            q_aligned = alignment_scalar * q_analytic
            q_aligned_free = q_aligned[np.asarray(free, dtype=int)]
            residual = residual_metrics(k_free, m_free, q_analytic, free, omega_sq=float(omega_sq))
            residual_free = residual_vector_free(k_free, m_free, q_analytic, free, omega_sq=float(omega_sq))
            left_sign, _ = convention_sign_and_scale(left_label)
            right_sign, _ = convention_sign_and_scale(right_label)
            rows.append(
                {
                    "theta_convention": f"left={left_label}; right={right_label}",
                    "theta_left_sign": float(left_sign),
                    "theta_right_sign": float(right_sign),
                    "theta_left_scale_description": theta_candidate_scale_description(left_label),
                    "theta_right_scale_description": theta_candidate_scale_description(right_label),
                    "relative_residual": float(residual["relative_residual"]),
                    "global_dof_mac_euclidean": euclidean_mac(q_analytic_free, q_fem_free),
                    "mass_weighted_mac": mass_weighted_mac(q_analytic_free, q_fem_free, m_free),
                    "translational_relative_l2": subset_relative_l2(q_aligned, q_fem_full, translational_mask),
                    "rotational_relative_l2": subset_relative_l2(q_aligned, q_fem_full, rotational_mask),
                    "rotational_mac": euclidean_mac(q_analytic[rotational_mask], q_fem_full[rotational_mask]),
                    "left_arm_relative_l2": subset_relative_l2(q_aligned, q_fem_full, left_mask),
                    "right_arm_relative_l2": subset_relative_l2(q_aligned, q_fem_full, right_mask),
                    "translational_mac": euclidean_mac(q_analytic[translational_mask], q_fem_full[translational_mask]),
                    "residual_l2": float(residual["residual_l2"]),
                    "residual_inf": float(residual["residual_inf"]),
                    "relative_residual_inf": float(residual["relative_residual_inf"]),
                    "translational_residual_norm": free_group_residual_norm(
                        residual_free, free, translational_mask
                    ),
                    "rotational_residual_norm": free_group_residual_norm(residual_free, free, rotational_mask),
                    "left_arm_residual_norm": free_group_residual_norm(residual_free, free, left_mask),
                    "right_arm_residual_norm": free_group_residual_norm(residual_free, free, right_mask),
                    "joint_residual_norm": free_group_residual_norm(residual_free, free, joint_mask),
                    "global_dof_relative_l2": safe_relative_l2(q_aligned_free - q_fem_free, q_fem_free),
                    "joint_relative_l2": subset_relative_l2(q_aligned, q_fem_full, joint_mask),
                    "left_theta_convention": left_label,
                    "right_theta_convention": right_label,
                    "global_alignment_scalar": float(alignment_scalar),
                }
            )
    rows.sort(key=lambda row: float(row["relative_residual"]))
    return rows


def write_theta_convention_scan_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=THETA_CONVENTION_SCAN_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def print_theta_convention_scan_top(
    rows: Sequence[dict[str, float | int | str]],
    *,
    baseline_relative_residual: float,
    limit: int = 6,
) -> None:
    if not rows:
        print("theta convention scan: no rows")
        return
    print("Top theta convention variants by analytic-in-FEM relative residual:")
    for row in list(rows)[:limit]:
        print(
            f"  residual={float(row['relative_residual']):.6e}, "
            f"global_MAC={float(row['global_dof_mac_euclidean']):.6e}, "
            f"rot_MAC={float(row['rotational_mac']):.6e}, {row['theta_convention']}"
        )
    best = rows[0]
    best_residual = float(best["relative_residual"])
    if best_residual < 0.1 * float(baseline_relative_residual):
        print("theta convention candidate found")
    else:
        print("theta convention scan: no dramatic residual reduction found")


def print_summary(row: dict[str, float | int | str]) -> None:
    print("Analytic shape in FEM residual diagnostic")
    print(f"branch_id: {row['branch_id']}")
    print(f"beta: {float(row['beta']):g} deg")
    print(f"mu: {float(row['mu']):g}")
    print(f"epsilon: {float(row['epsilon']):g}")
    print(f"analytic Lambda: {float(row['analytic_lambda']):.10g}")
    print(f"FEM tracked Lambda: {float(row['fem_lambda']):.10g}")
    print(f"relative Lambda difference: {float(row['rel_lambda_diff']):.6e}")
    print(f"analytic current_sorted_index: {int(row['current_sorted_index'])}")
    print(f"FEM current_sorted_index: {int(row['fem_sorted_index'])}")
    print(f"right theta sign used: {float(row['right_theta_sign_used']):+.0f}")
    print(f"relative residual with FEM omega: {float(row['relative_residual']):.6e}")
    print(f"selected FEM self relative residual: {float(row['fem_selected_self_relative_residual']):.6e}")
    print(f"relative residual with analytic omega: {float(row['selected_analytic_omega_relative_residual']):.6e}")
    print(f"global DOF MAC: {float(row['global_dof_mac_euclidean']):.6e}")
    print(f"mass-weighted DOF MAC: {float(row['global_dof_mac_mass_weighted']):.6e}")
    print(f"global DOF relative L2: {float(row['global_dof_relative_l2']):.6e}")
    print(f"analytic q->local relative L2: {float(row['analytic_embedding_local_relative_l2']):.6e}")
    print(f"analytic q->local vs FEM local relative L2: {float(row['analytic_vs_fem_local_relative_l2']):.6e}")
    print(
        "FEM theta derivative best: "
        f"left={row['best_theta_convention_left']} ({row['best_theta_scale_left']}), "
        f"right={row['best_theta_convention_right']} ({row['best_theta_scale_right']})"
    )
    if str(row.get("theta_convention_scan_csv", "")):
        print(f"theta convention scan CSV: {row['theta_convention_scan_csv']}")
        print(
            "theta convention scan best: "
            f"{row['theta_convention_scan_best']} "
            f"relative_residual={float(row['theta_convention_scan_best_relative_residual']):.6e}"
        )
    if str(row.get("modal_projection_csv", "")):
        print(
            "Rayleigh omega_sq: "
            f"{float(row['analytic_rayleigh_omega_sq']):.6e} "
            f"(relative diff {float(row['relative_rayleigh_omega_sq_diff']):.6e})"
        )
        print(
            "top mass projection: "
            f"mode {int(row['top_modal_projection_mode'])}, "
            f"share={float(row['top_modal_projection_share']):.6e}, "
            f"top3={float(row['cumulative_projection_share_first_3']):.6e}"
        )
        print(
            "largest residual group: "
            f"{row['residual_largest_group']} "
            f"share={float(row['residual_largest_group_share']):.6e}"
        )
        print(
            "largest element residual: "
            f"{row['residual_largest_element_arm']}[{int(row['residual_largest_element_index'])}] "
            f"share={float(row['residual_largest_element_share']):.6e}"
        )
        print(f"modal projection CSV: {row['modal_projection_csv']}")
        print(f"residual groups CSV: {row['residual_groups_csv']}")
        print(f"element residual/energy CSV: {row['element_residual_energy_csv']}")
    if str(row.get("bending_basis_pairing_scan_csv", "")):
        print(
            "bending basis scan best: "
            f"{row['bending_basis_best_variant']} "
            f"relative_residual={float(row['bending_basis_best_relative_residual']):.6e}, "
            f"matrix_consistent={row['bending_basis_best_preserves_matrix_consistency']}"
        )
        print(f"bending basis pairing scan CSV: {row['bending_basis_pairing_scan_csv']}")
        print(f"full basis column audit CSV: {row['full_basis_column_audit_csv']}")
    if str(row.get("axial_convention_scan_csv", "")):
        print(
            "axial convention scan best: "
            f"{row['axial_convention_best_variant']} "
            f"relative_residual={float(row['axial_convention_best_relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(row['axial_convention_best_rayleigh_lambda']):.6e}"
        )
        print(f"axial convention scan CSV: {row['axial_convention_scan_csv']}")
    if str(row.get("axial_scale_scan_csv", "")):
        print(
            "axial scale scan best by residual: "
            f"{row['axial_scale_best_residual_variant']} "
            f"alpha=({float(row['axial_scale_best_residual_alpha_left']):.6g}, "
            f"{float(row['axial_scale_best_residual_alpha_right']):.6g}), "
            f"relative_residual={float(row['axial_scale_best_residual_relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(row['axial_scale_best_residual_rayleigh_lambda']):.6e}"
        )
        print(
            "axial scale scan best by Rayleigh: "
            f"{row['axial_scale_best_rayleigh_variant']} "
            f"alpha=({float(row['axial_scale_best_rayleigh_alpha_left']):.6g}, "
            f"{float(row['axial_scale_best_rayleigh_alpha_right']):.6g}), "
            f"relative_residual={float(row['axial_scale_best_rayleigh_relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(row['axial_scale_best_rayleigh_lambda']):.6e}"
        )
        print(f"axial scale scan CSV: {row['axial_scale_scan_csv']}")
    if str(row.get("right_transform_convention_scan_csv", "")):
        print(f"FEM transform inference: {row['fem_transform_inference']}")
        print(
            "right transform scan best by residual: "
            f"{row['right_transform_best_residual_variant']} "
            f"relative_residual={float(row['right_transform_best_residual_relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(row['right_transform_best_residual_rayleigh_lambda']):.6e}"
        )
        print(
            "right transform scan best by Rayleigh: "
            f"{row['right_transform_best_rayleigh_variant']} "
            f"relative_residual={float(row['right_transform_best_rayleigh_relative_residual']):.6e}, "
            f"Rayleigh Lambda={float(row['right_transform_best_rayleigh_lambda']):.6e}"
        )
        print(f"right transform convention scan CSV: {row['right_transform_convention_scan_csv']}")
        print(f"right transform sanity CSV: {row['right_transform_sanity_csv']}")
    if str(row.get("frequency_scaling_audit_csv", "")):
        print(
            "frequency scaling best for analytic q: "
            f"{row['frequency_scaling_best_analytic_factor']} "
            f"relative_residual={float(row['frequency_scaling_best_analytic_relative_residual']):.6e}"
        )
        print(
            "frequency scaling best for FEM self residual: "
            f"{row['frequency_scaling_best_fem_factor']} "
            f"relative_residual={float(row['frequency_scaling_best_fem_relative_residual']):.6e}"
        )
        print(f"frequency scaling audit CSV: {row['frequency_scaling_audit_csv']}")
    print(f"conclusion_flag: {row['conclusion_flag']}")
    print(f"summary CSV: {row['summary_csv']}")
    print(f"debug nodes CSV: {row['debug_nodes_csv']}")


def main(argv: Sequence[str] | None = None) -> dict[str, float | int | str]:
    args = parse_args(argv)
    branch_id = str(args.branch_id)
    params = build_params_for_case(epsilon=float(args.epsilon), l_total=float(args.l_total))
    analytic_tracking = resolve_analytic_branch(
        branch_id=branch_id,
        beta_deg=float(args.beta),
        mu=float(args.mu),
        epsilon=float(args.epsilon),
        num_analytic_roots=int(args.num_analytic_roots),
    )
    point = analytic_tracking["point"]
    analytic_lambda = float(point.Lambda)
    beta_rad = float(np.deg2rad(args.beta))
    matrix = assemble_clamped_coupled_matrix(analytic_lambda, beta_rad, float(args.mu), float(args.epsilon))
    coeff, smallest_singular_value, singular_value_ratio = analytic_null_vector(matrix)
    matrix_residual = matrix @ coeff

    fem_branch = track_fem_branch(
        params,
        branch_id=branch_id,
        beta_deg=float(args.beta),
        mu=float(args.mu),
    )
    k_free, m_free, k_full, m_full, free = assemble_fem_matrices(params, mu=float(args.mu), beta_deg=float(args.beta))
    target_omega, _, _ = solve_fem_modes_nd(
        params,
        mu=float(args.mu),
        beta_deg=float(args.beta),
        n_modes=FEM_N_SOLVE,
    )
    target_lambdas = np.sqrt(np.maximum(target_omega, 0.0))
    nearest_idx0 = int(np.argmin(np.abs(target_lambdas - analytic_lambda)))
    nearest_lambda = float(target_lambdas[nearest_idx0])
    nearest_rel_diff = (
        abs(nearest_lambda - analytic_lambda) / abs(analytic_lambda)
        if abs(analytic_lambda) > NEAR_ZERO_NORM
        else np.nan
    )

    s_norm = np.linspace(0.0, 1.0, fem.N_ELEM + 1)
    analytic_omega_nd = analytic_lambda**2
    analytic_omega_sq = analytic_lambda**4
    fem_omega_sq = fem_branch.omega_nd**2

    sign_results: dict[float, dict[str, object]] = {}
    for theta_sign in RIGHT_THETA_SIGNS:
        fields = analytic_local_fields(
            analytic_lambda,
            mu=float(args.mu),
            epsilon=float(args.epsilon),
            coeff=coeff,
            s_norm=s_norm,
            right_theta_sign=theta_sign,
        )
        q_analytic, embedding = build_analytic_global_q(fields, beta_deg=float(args.beta))
        fem_residual = residual_metrics(k_free, m_free, q_analytic, free, omega_sq=fem_omega_sq)
        analytic_residual = residual_metrics(k_free, m_free, q_analytic, free, omega_sq=analytic_omega_sq)
        sign_results[theta_sign] = {
            "fields": fields,
            "q": q_analytic,
            "embedding": embedding,
            "fem_residual": fem_residual,
            "analytic_residual": analytic_residual,
        }

    selected_sign = min(
        RIGHT_THETA_SIGNS,
        key=lambda sign: float(sign_results[sign]["fem_residual"]["relative_residual"]),  # type: ignore[index]
    )
    selected = sign_results[selected_sign]
    q_analytic = np.asarray(selected["q"], dtype=float)
    fields_selected = selected["fields"]  # type: ignore[assignment]
    q_fem = np.asarray(fem_branch.full_vec, dtype=float)
    fem_self_residual = residual_metrics(k_free, m_free, q_fem, free, omega_sq=fem_omega_sq)
    q_analytic_free = q_analytic[free]
    q_fem_free = q_fem[free]
    alignment_scalar = optimal_scalar(q_analytic_free, q_fem_free)
    q_analytic_aligned = alignment_scalar * q_analytic
    q_analytic_free_aligned = q_analytic_aligned[free]

    n = fem.N_ELEM
    full_dof_count = len(q_fem)
    all_indices = np.arange(full_dof_count, dtype=int)
    translational_mask = (all_indices % 3) != 2
    rotational_mask = (all_indices % 3) == 2
    left_mask = all_indices <= (3 * n + 2)
    right_mask = all_indices >= 3 * n
    joint_mask = (all_indices >= 3 * n) & (all_indices <= 3 * n + 2)

    analytic_q_local = local_fields_from_q(q_analytic, beta_deg=float(args.beta))
    direct_local_vector = local_component_vector(fields_selected)  # type: ignore[arg-type]
    q_local_vector = local_component_vector(analytic_q_local)
    local_embedding_diff = q_local_vector - direct_local_vector
    fem_local = local_fields_from_q(q_fem, beta_deg=float(args.beta))
    analytic_aligned_local = local_fields_from_q(q_analytic_aligned, beta_deg=float(args.beta))
    fem_local_vector = local_component_vector(fem_local)
    analytic_aligned_local_vector = local_component_vector(analytic_aligned_local)
    theta_diagnostics = fem_theta_derivative_diagnostics(
        fem_local,
        Lambda=float(fem_branch.lambda_value),
        mu=float(args.mu),
        l_base=float(params.L_base),
    )

    rel_lambda_diff = (
        abs(analytic_lambda - fem_branch.lambda_value) / abs(fem_branch.lambda_value)
        if abs(fem_branch.lambda_value) > NEAR_ZERO_NORM
        else np.nan
    )
    selected_fem_residual = selected["fem_residual"]  # type: ignore[assignment]
    selected_analytic_residual = selected["analytic_residual"]  # type: ignore[assignment]

    pos_fem = sign_results[1.0]["fem_residual"]  # type: ignore[index]
    neg_fem = sign_results[-1.0]["fem_residual"]  # type: ignore[index]
    pos_analytic = sign_results[1.0]["analytic_residual"]  # type: ignore[index]
    neg_analytic = sign_results[-1.0]["analytic_residual"]  # type: ignore[index]
    selected_embedding = selected["embedding"]  # type: ignore[assignment]

    output_prefix = resolve_output_prefix(
        args.output_prefix,
        branch_id=branch_id,
        beta_deg=float(args.beta),
        mu_value=float(args.mu),
        epsilon=float(args.epsilon),
    )
    summary_csv = suffixed_path(output_prefix, ".csv")
    nodes_csv = suffixed_path(output_prefix, "_nodes.csv") if bool(args.save_debug_csv) else Path("")
    theta_scan_csv = suffixed_path(output_prefix, "_theta_convention_scan.csv")
    modal_projection_csv = suffixed_path(output_prefix, "_modal_projection.csv")
    residual_groups_csv = suffixed_path(output_prefix, "_residual_groups.csv")
    element_residual_energy_csv = suffixed_path(output_prefix, "_element_residual_energy.csv")
    bending_basis_scan_csv = suffixed_path(output_prefix, "_bending_basis_pairing_scan.csv")
    full_basis_column_audit_csv = suffixed_path(output_prefix, "_full_basis_column_audit.csv")
    axial_convention_scan_csv = suffixed_path(output_prefix, "_axial_convention_scan.csv")
    axial_scale_scan_csv = suffixed_path(output_prefix, "_axial_scale_scan.csv")
    right_transform_scan_csv = suffixed_path(output_prefix, "_right_transform_convention_scan.csv")
    right_transform_sanity_csv = suffixed_path(output_prefix, "_right_transform_sanity.csv")
    frequency_scaling_audit_csv = suffixed_path(output_prefix, "_frequency_scaling_audit.csv")
    theta_scan_rows: list[dict[str, float | str]] = []
    theta_scan_best = ""
    theta_scan_best_relative_residual = np.nan
    if bool(args.scan_theta_conventions):
        theta_scan_rows = evaluate_theta_convention_scan(
            analytic_lambda=analytic_lambda,
            mu=float(args.mu),
            epsilon=float(args.epsilon),
            coeff=coeff,
            s_norm=s_norm,
            l_base=float(params.L_base),
            beta_deg=float(args.beta),
            k_free=k_free,
            m_free=m_free,
            free=free,
            omega_sq=fem_omega_sq,
            q_fem=q_fem,
            translational_mask=translational_mask,
            rotational_mask=rotational_mask,
            left_mask=left_mask,
            right_mask=right_mask,
            joint_mask=joint_mask,
        )
        write_theta_convention_scan_csv(theta_scan_csv, theta_scan_rows)
        if theta_scan_rows:
            theta_scan_best = str(theta_scan_rows[0]["theta_convention"])
            theta_scan_best_relative_residual = float(theta_scan_rows[0]["relative_residual"])

    rayleigh = rayleigh_diagnostics(
        k_free,
        m_free,
        q_analytic,
        free,
        fem_omega_sq=fem_omega_sq,
        fem_lambda=float(fem_branch.lambda_value),
    )
    frequency_scaling_rows: list[dict[str, float | str]] = []
    frequency_scaling_best_analytic_factor = ""
    frequency_scaling_best_analytic_relative_residual = np.nan
    frequency_scaling_best_fem_factor = ""
    frequency_scaling_best_fem_relative_residual = np.nan
    if bool(args.audit_frequency_scaling):
        frequency_scaling_rows = frequency_scaling_audit_rows(
            k_free,
            m_free,
            q_analytic,
            q_fem,
            free,
            analytic_lambda=analytic_lambda,
            fem_selected_omega=float(fem_branch.omega_nd),
            fem_selected_omega_sq=float(fem_omega_sq),
            rayleigh_omega_sq=float(rayleigh["analytic_rayleigh_omega_sq"]),
        )
        write_frequency_scaling_audit_csv(frequency_scaling_audit_csv, frequency_scaling_rows)
        best_analytic_frequency = frequency_scaling_best_row(frequency_scaling_rows, vector_type="analytic_q")
        if best_analytic_frequency is not None:
            frequency_scaling_best_analytic_factor = str(best_analytic_frequency["spectral_factor_name"])
            frequency_scaling_best_analytic_relative_residual = float(best_analytic_frequency["relative_residual"])
        best_fem_frequency = frequency_scaling_best_row(frequency_scaling_rows, vector_type="fem_selected_q")
        if best_fem_frequency is not None:
            frequency_scaling_best_fem_factor = str(best_fem_frequency["spectral_factor_name"])
            frequency_scaling_best_fem_relative_residual = float(best_fem_frequency["relative_residual"])
    kinematic = q_kinematic_checks(q_analytic, selected_embedding)  # type: ignore[arg-type]
    modal_rows: list[dict[str, float | int]] = []
    residual_group_data: list[dict[str, float | int | str]] = []
    element_rows: list[dict[str, float | int | str]] = []
    top_modal_projection_mode = 0
    top_modal_projection_share = np.nan
    cumulative_projection_share_first_3 = np.nan
    residual_largest_group = ""
    residual_largest_group_share = np.nan
    residual_largest_element_arm = ""
    residual_largest_element_index = -1
    residual_largest_element_share = np.nan
    if bool(args.localize_residual):
        modal_rows = modal_projection_rows(
            params,
            mu=float(args.mu),
            beta_deg=float(args.beta),
            q_full=q_analytic_aligned,
            free=free,
            mass_matrix=m_free,
            n_modes=DEFAULT_MODAL_PROJECTION_MODES,
        )
        write_modal_projection_csv(modal_projection_csv, modal_rows)
        if modal_rows:
            top_modal_projection_mode = int(modal_rows[0]["mode_index"])
            top_modal_projection_share = float(modal_rows[0]["mass_weighted_projection_share"])
            cumulative_projection_share_first_3 = float(
                sum(
                    float(row["mass_weighted_projection_share"])
                    for row in modal_rows[:3]
                    if np.isfinite(float(row["mass_weighted_projection_share"]))
                )
            )

        residual_group_data = residual_group_rows(
            k_full,
            m_full,
            q_analytic_aligned,
            free,
            omega_sq=fem_omega_sq,
        )
        write_residual_groups_csv(residual_groups_csv, residual_group_data)
        free_group_rows = [
            row for row in residual_group_data if str(row["residual_scope"]) == "free_equilibrium"
        ]
        if free_group_rows:
            largest_group = max(free_group_rows, key=lambda row: float(row["group_residual_share"]))
            residual_largest_group = str(largest_group["group_name"])
            residual_largest_group_share = float(largest_group["group_residual_share"])

        element_rows = element_residual_energy_rows(
            params,
            mu=float(args.mu),
            beta_deg=float(args.beta),
            q_analytic=q_analytic_aligned,
            q_fem=q_fem,
            omega_sq=fem_omega_sq,
        )
        write_element_residual_energy_csv(element_residual_energy_csv, element_rows)
        if element_rows:
            largest_element = max(element_rows, key=lambda row: float(row["local_residual_norm_share"]))
            residual_largest_element_arm = str(largest_element["arm"])
            residual_largest_element_index = int(largest_element["element_index"])
            residual_largest_element_share = float(largest_element["local_residual_norm_share"])

    bending_basis_rows: list[dict[str, float | str]] = []
    full_basis_column_rows: list[dict[str, float | str]] = []
    bending_basis_canonical_relative_residual = np.nan
    bending_basis_best_variant = ""
    bending_basis_best_relative_residual = np.nan
    bending_basis_best_preserves_matrix_consistency = ""
    bending_basis_pairing_mismatch_flag = ""
    if bool(args.scan_bending_basis_pairings):
        full_basis_column_rows = full_basis_column_audit_rows(
            matrix,
            Lambda=analytic_lambda,
            beta_rad=beta_rad,
            mu=float(args.mu),
            epsilon=float(args.epsilon),
        )
        write_full_basis_column_audit_csv(full_basis_column_audit_csv, full_basis_column_rows)
        bending_basis_rows = bending_basis_pairing_scan_rows(
            Lambda=analytic_lambda,
            beta_deg=float(args.beta),
            beta_rad=beta_rad,
            mu=float(args.mu),
            epsilon=float(args.epsilon),
            coeff=coeff,
            s_norm=s_norm,
            l_base=float(params.L_base),
            matrix=matrix,
            k_free=k_free,
            m_free=m_free,
            free=free,
            omega_sq=fem_omega_sq,
            q_fem=q_fem,
            translational_mask=translational_mask,
            rotational_mask=rotational_mask,
        )
        write_bending_basis_pairing_scan_csv(bending_basis_scan_csv, bending_basis_rows)
        canonical_basis = canonical_pairing_row(bending_basis_rows)
        if canonical_basis is not None:
            bending_basis_canonical_relative_residual = float(canonical_basis["analytic_in_fem_relative_residual"])
        if bending_basis_rows:
            best_basis = bending_basis_rows[0]
            bending_basis_best_variant = str(best_basis["variant_name"])
            bending_basis_best_relative_residual = float(best_basis["analytic_in_fem_relative_residual"])
            best_consistent = (
                float(best_basis["matrix_vs_independent_endpoint_mismatch"]) <= 1e-7
                and float(best_basis["external_clamp_residual"]) <= 1e-9
            )
            bending_basis_best_preserves_matrix_consistency = "yes" if best_consistent else "no"
            bending_basis_pairing_mismatch_flag = (
                "possible_coefficient_basis_pairing_mismatch_detected"
                if str(best_basis["conclusion_flag"])
                in {"improves_residual_and_consistent", "improves_residual_but_breaks_matrix_consistency"}
                else "no_supported_pairing_mismatch"
            )

    axial_convention_rows: list[dict[str, float | str]] = []
    axial_convention_baseline_relative_residual = np.nan
    axial_convention_baseline_rayleigh_lambda = np.nan
    axial_convention_best_variant = ""
    axial_convention_best_relative_residual = np.nan
    axial_convention_best_rayleigh_lambda = np.nan
    axial_convention_mismatch_flag = ""
    if bool(args.scan_axial_conventions):
        axial_convention_rows = axial_convention_scan_rows(
            fields_base=fields_selected,  # type: ignore[arg-type]
            beta_deg=float(args.beta),
            mu=float(args.mu),
            l_base=float(params.L_base),
            params=params,
            k_free=k_free,
            m_free=m_free,
            k_full=k_full,
            m_full=m_full,
            free=free,
            omega_sq=fem_omega_sq,
            fem_lambda=float(fem_branch.lambda_value),
            q_fem=q_fem,
            translational_mask=translational_mask,
            rotational_mask=rotational_mask,
        )
        write_axial_convention_scan_csv(axial_convention_scan_csv, axial_convention_rows)
        baseline_axial = next(
            (row for row in axial_convention_rows if row["variant_name"] == "left+_right_as_is"),
            None,
        )
        if baseline_axial is not None:
            axial_convention_baseline_relative_residual = float(baseline_axial["relative_residual"])
            axial_convention_baseline_rayleigh_lambda = float(baseline_axial["rayleigh_lambda"])
        if axial_convention_rows:
            best_axial = axial_convention_rows[0]
            axial_convention_best_variant = str(best_axial["variant_name"])
            axial_convention_best_relative_residual = float(best_axial["relative_residual"])
            axial_convention_best_rayleigh_lambda = float(best_axial["rayleigh_lambda"])
            axial_convention_mismatch_flag = (
                "possible_axial_convention_mismatch"
                if any(
                    str(row["conclusion_flag"]) == "possible_axial_convention_mismatch"
                    for row in axial_convention_rows
                )
                else "simple_sign_reversal_axial_convention_does_not_explain_mismatch"
            )

    axial_scale_rows: list[dict[str, float | str]] = []
    axial_scale_baseline_relative_residual = np.nan
    axial_scale_baseline_rayleigh_lambda = np.nan
    axial_scale_best_residual_variant = ""
    axial_scale_best_residual_alpha_left = np.nan
    axial_scale_best_residual_alpha_right = np.nan
    axial_scale_best_residual_relative_residual = np.nan
    axial_scale_best_residual_rayleigh_lambda = np.nan
    axial_scale_best_rayleigh_variant = ""
    axial_scale_best_rayleigh_alpha_left = np.nan
    axial_scale_best_rayleigh_alpha_right = np.nan
    axial_scale_best_rayleigh_relative_residual = np.nan
    axial_scale_best_rayleigh_lambda = np.nan
    axial_scale_best_mass_mac_variant = ""
    axial_scale_best_mass_mac = np.nan
    axial_scale_mismatch_flag = ""
    if bool(args.scan_axial_scales):
        axial_scale_rows = axial_scale_scan_rows(
            fields_base=fields_selected,  # type: ignore[arg-type]
            Lambda=analytic_lambda,
            beta_deg=float(args.beta),
            mu=float(args.mu),
            epsilon=float(args.epsilon),
            params=params,
            k_free=k_free,
            m_free=m_free,
            k_full=k_full,
            m_full=m_full,
            free=free,
            omega_sq=fem_omega_sq,
            fem_lambda=float(fem_branch.lambda_value),
            q_fem=q_fem,
            translational_mask=translational_mask,
            rotational_mask=rotational_mask,
        )
        write_axial_scale_scan_csv(axial_scale_scan_csv, axial_scale_rows)
        baseline_scale = next(
            (row for row in axial_scale_rows if str(row["variant_name"]) == "baseline_alpha_1_1"),
            None,
        )
        if baseline_scale is not None:
            axial_scale_baseline_relative_residual = float(baseline_scale["relative_residual"])
            axial_scale_baseline_rayleigh_lambda = float(baseline_scale["rayleigh_lambda"])
        if axial_scale_rows:
            best_residual_scale = axial_scale_rows[0]
            axial_scale_best_residual_variant = str(best_residual_scale["variant_name"])
            axial_scale_best_residual_alpha_left = float(best_residual_scale["alpha_left"])
            axial_scale_best_residual_alpha_right = float(best_residual_scale["alpha_right"])
            axial_scale_best_residual_relative_residual = float(best_residual_scale["relative_residual"])
            axial_scale_best_residual_rayleigh_lambda = float(best_residual_scale["rayleigh_lambda"])
            best_rayleigh_scale = best_axial_scale_by_rayleigh(axial_scale_rows)
            if best_rayleigh_scale is not None:
                axial_scale_best_rayleigh_variant = str(best_rayleigh_scale["variant_name"])
                axial_scale_best_rayleigh_alpha_left = float(best_rayleigh_scale["alpha_left"])
                axial_scale_best_rayleigh_alpha_right = float(best_rayleigh_scale["alpha_right"])
                axial_scale_best_rayleigh_relative_residual = float(best_rayleigh_scale["relative_residual"])
                axial_scale_best_rayleigh_lambda = float(best_rayleigh_scale["rayleigh_lambda"])
            best_mass_scale = best_axial_scale_by_mass_mac(axial_scale_rows)
            if best_mass_scale is not None:
                axial_scale_best_mass_mac_variant = str(best_mass_scale["variant_name"])
                axial_scale_best_mass_mac = float(best_mass_scale["mass_weighted_mac"])
            axial_scale_mismatch_flag = axial_scale_overall_flag(axial_scale_rows)

    right_transform_rows: list[dict[str, float | str]] = []
    right_transform_sanity_data: list[dict[str, float | str]] = []
    right_transform_baseline_relative_residual = np.nan
    right_transform_baseline_rayleigh_lambda = np.nan
    right_transform_best_residual_variant = ""
    right_transform_best_residual_relative_residual = np.nan
    right_transform_best_residual_rayleigh_lambda = np.nan
    right_transform_best_rayleigh_variant = ""
    right_transform_best_rayleigh_relative_residual = np.nan
    right_transform_best_rayleigh_lambda = np.nan
    right_transform_best_mass_mac_variant = ""
    right_transform_best_mass_mac = np.nan
    right_transform_mismatch_flag = ""
    if bool(args.scan_right_transform_conventions):
        transform_context = element_diagnostic_context(params, mu=float(args.mu), beta_deg=float(args.beta))
        right_transform_rows = right_transform_convention_scan_rows(
            fields_base=fields_selected,  # type: ignore[arg-type]
            beta_deg=float(args.beta),
            params=params,
            k_free=k_free,
            m_free=m_free,
            k_full=k_full,
            m_full=m_full,
            free=free,
            omega_sq=fem_omega_sq,
            fem_lambda=float(fem_branch.lambda_value),
            q_fem=q_fem,
            translational_mask=translational_mask,
            rotational_mask=rotational_mask,
            element_context=transform_context,
        )
        right_transform_sanity_data = right_transform_sanity_rows(float(args.beta))
        write_right_transform_convention_scan_csv(right_transform_scan_csv, right_transform_rows)
        write_right_transform_sanity_csv(right_transform_sanity_csv, right_transform_sanity_data)
        baseline_transform = next(
            (row for row in right_transform_rows if str(row["variant_name"]) == "baseline_current"),
            None,
        )
        if baseline_transform is not None:
            right_transform_baseline_relative_residual = float(baseline_transform["relative_residual"])
            right_transform_baseline_rayleigh_lambda = float(baseline_transform["rayleigh_lambda"])
        if right_transform_rows:
            best_residual_transform = right_transform_rows[0]
            right_transform_best_residual_variant = str(best_residual_transform["variant_name"])
            right_transform_best_residual_relative_residual = float(best_residual_transform["relative_residual"])
            right_transform_best_residual_rayleigh_lambda = float(best_residual_transform["rayleigh_lambda"])
            best_rayleigh_transform = best_right_transform_by_rayleigh(right_transform_rows)
            if best_rayleigh_transform is not None:
                right_transform_best_rayleigh_variant = str(best_rayleigh_transform["variant_name"])
                right_transform_best_rayleigh_relative_residual = float(best_rayleigh_transform["relative_residual"])
                right_transform_best_rayleigh_lambda = float(best_rayleigh_transform["rayleigh_lambda"])
            best_mass_transform = best_right_transform_by_mass_mac(right_transform_rows)
            if best_mass_transform is not None:
                right_transform_best_mass_mac_variant = str(best_mass_transform["variant_name"])
                right_transform_best_mass_mac = float(best_mass_transform["mass_weighted_mac"])
            right_transform_mismatch_flag = right_transform_overall_flag(right_transform_rows)

    row: dict[str, float | int | str] = {
        "branch_id": branch_id,
        "beta": float(args.beta),
        "mu": float(args.mu),
        "epsilon": float(args.epsilon),
        "base_sorted_index": int(point.base_sorted_index),
        "current_sorted_index": int(point.current_sorted_index),
        "analytic_lambda": analytic_lambda,
        "fem_lambda": float(fem_branch.lambda_value),
        "rel_lambda_diff": float(rel_lambda_diff),
        "branch_tracking_mac_min": float(analytic_tracking["mac_min"]),
        "branch_tracking_warning_flag": str(analytic_tracking["warning_flag"]),
        "fem_sorted_index": int(fem_branch.current_sorted_index),
        "residual_l2": float(selected_fem_residual["residual_l2"]),  # type: ignore[index]
        "relative_residual": float(selected_fem_residual["relative_residual"]),  # type: ignore[index]
        "residual_inf": float(selected_fem_residual["residual_inf"]),  # type: ignore[index]
        "relative_residual_inf": float(selected_fem_residual["relative_residual_inf"]),  # type: ignore[index]
        "global_dof_relative_l2": safe_relative_l2(q_analytic_free_aligned - q_fem_free, q_fem_free),
        "global_dof_mac_euclidean": euclidean_mac(q_analytic_free, q_fem_free),
        "global_dof_mac_mass_weighted": mass_weighted_mac(q_analytic_free, q_fem_free, m_free),
        "translational_relative_l2": subset_relative_l2(q_analytic_aligned, q_fem, translational_mask),
        "rotational_relative_l2": subset_relative_l2(q_analytic_aligned, q_fem, rotational_mask),
        "left_arm_relative_l2": subset_relative_l2(q_analytic_aligned, q_fem, left_mask),
        "right_arm_relative_l2": subset_relative_l2(q_analytic_aligned, q_fem, right_mask),
        "joint_relative_l2": subset_relative_l2(q_analytic_aligned, q_fem, joint_mask),
        "conclusion_flag": "",
        "l_total": float(args.l_total),
        "radius": float(params.r),
        "analytic_omega_nd": float(analytic_omega_nd),
        "analytic_eigenvalue_omega_sq": float(analytic_omega_sq),
        "fem_omega_nd": float(fem_branch.omega_nd),
        "fem_eigenvalue_omega_sq": float(fem_omega_sq),
        "fem_selected_self_residual_l2": float(fem_self_residual["residual_l2"]),
        "fem_selected_self_relative_residual": float(fem_self_residual["relative_residual"]),
        "fem_frequency_hz": float(fem_branch.frequency_hz),
        "right_theta_sign_used": float(selected_sign),
        "right_theta_sign_best_by_fem_residual": float(selected_sign),
        "residual_l2_theta_pos_fem_omega": float(pos_fem["residual_l2"]),  # type: ignore[index]
        "relative_residual_theta_pos_fem_omega": float(pos_fem["relative_residual"]),  # type: ignore[index]
        "residual_inf_theta_pos_fem_omega": float(pos_fem["residual_inf"]),  # type: ignore[index]
        "relative_residual_inf_theta_pos_fem_omega": float(pos_fem["relative_residual_inf"]),  # type: ignore[index]
        "residual_l2_theta_neg_fem_omega": float(neg_fem["residual_l2"]),  # type: ignore[index]
        "relative_residual_theta_neg_fem_omega": float(neg_fem["relative_residual"]),  # type: ignore[index]
        "residual_inf_theta_neg_fem_omega": float(neg_fem["residual_inf"]),  # type: ignore[index]
        "relative_residual_inf_theta_neg_fem_omega": float(neg_fem["relative_residual_inf"]),  # type: ignore[index]
        "residual_l2_theta_pos_analytic_omega": float(pos_analytic["residual_l2"]),  # type: ignore[index]
        "relative_residual_theta_pos_analytic_omega": float(pos_analytic["relative_residual"]),  # type: ignore[index]
        "residual_l2_theta_neg_analytic_omega": float(neg_analytic["residual_l2"]),  # type: ignore[index]
        "relative_residual_theta_neg_analytic_omega": float(neg_analytic["relative_residual"]),  # type: ignore[index]
        "selected_analytic_omega_relative_residual": float(selected_analytic_residual["relative_residual"]),  # type: ignore[index]
        "nearest_fem_sorted_index_by_lambda": int(nearest_idx0 + 1),
        "nearest_fem_lambda_by_lambda": float(nearest_lambda),
        "nearest_fem_rel_lambda_diff_by_lambda": float(nearest_rel_diff),
        "analytic_fem_sorted_index_delta": int(point.current_sorted_index) - int(fem_branch.current_sorted_index),
        "branch_tracking_mac_target": float(point.mac_to_previous),
        "branch_tracking_mac_mean": float(analytic_tracking["mac_mean"]),
        "branch_tracking_path_points": int(len(analytic_tracking["branch_points"])),
        "smallest_singular_value": float(smallest_singular_value),
        "singular_value_ratio": float(singular_value_ratio),
        "analytic_matrix_residual_inf": float(np.max(np.abs(matrix_residual))),
        "analytic_embedding_local_relative_l2": safe_relative_l2(local_embedding_diff, direct_local_vector),
        "analytic_embedding_local_max_abs": float(np.max(np.abs(local_embedding_diff))),
        "analytic_vs_fem_local_relative_l2": safe_relative_l2(
            analytic_aligned_local_vector - fem_local_vector,
            fem_local_vector,
        ),
        "analytic_vs_fem_local_mac_euclidean": euclidean_mac(analytic_aligned_local_vector, fem_local_vector),
        "joint_translation_embedding_mismatch": float(selected_embedding["joint_translation_embedding_mismatch"]),  # type: ignore[index]
        "joint_theta_embedding_mismatch": float(selected_embedding["joint_theta_embedding_mismatch"]),  # type: ignore[index]
        "global_alignment_scalar": float(alignment_scalar),
        "rel_l2_theta_vs_dwds_left": float(theta_diagnostics["rel_l2_theta_vs_dwds_left"]),
        "rel_l2_theta_vs_minus_dwds_left": float(theta_diagnostics["rel_l2_theta_vs_minus_dwds_left"]),
        "rel_l2_theta_vs_dwds_right": float(theta_diagnostics["rel_l2_theta_vs_dwds_right"]),
        "rel_l2_theta_vs_minus_dwds_right": float(theta_diagnostics["rel_l2_theta_vs_minus_dwds_right"]),
        "rel_l2_theta_vs_dwdxi_left": float(theta_diagnostics["rel_l2_theta_vs_dwdxi_left"]),
        "rel_l2_theta_vs_minus_dwdxi_left": float(theta_diagnostics["rel_l2_theta_vs_minus_dwdxi_left"]),
        "rel_l2_theta_vs_dwdxi_right": float(theta_diagnostics["rel_l2_theta_vs_dwdxi_right"]),
        "rel_l2_theta_vs_minus_dwdxi_right": float(theta_diagnostics["rel_l2_theta_vs_minus_dwdxi_right"]),
        "rel_l2_theta_vs_dwdz_left": float(theta_diagnostics["rel_l2_theta_vs_dwdz_left"]),
        "rel_l2_theta_vs_minus_dwdz_left": float(theta_diagnostics["rel_l2_theta_vs_minus_dwdz_left"]),
        "rel_l2_theta_vs_dwdz_right": float(theta_diagnostics["rel_l2_theta_vs_dwdz_right"]),
        "rel_l2_theta_vs_minus_dwdz_right": float(theta_diagnostics["rel_l2_theta_vs_minus_dwdz_right"]),
        "best_theta_convention_left": str(theta_diagnostics["best_theta_convention_left"]),
        "best_theta_convention_right": str(theta_diagnostics["best_theta_convention_right"]),
        "best_theta_scale_left": str(theta_diagnostics["best_theta_scale_left"]),
        "best_theta_scale_right": str(theta_diagnostics["best_theta_scale_right"]),
        "theta_convention_scan_csv": str(theta_scan_csv) if bool(args.scan_theta_conventions) else "",
        "theta_convention_scan_best": theta_scan_best,
        "theta_convention_scan_best_relative_residual": float(theta_scan_best_relative_residual),
        "analytic_rayleigh_omega_sq": float(rayleigh["analytic_rayleigh_omega_sq"]),
        "fem_selected_omega_sq": float(rayleigh["fem_selected_omega_sq"]),
        "relative_rayleigh_omega_sq_diff": float(rayleigh["relative_rayleigh_omega_sq_diff"]),
        "analytic_rayleigh_Lambda": float(rayleigh["analytic_rayleigh_Lambda"]),
        "relative_rayleigh_Lambda_diff": float(rayleigh["relative_rayleigh_Lambda_diff"]),
        "top_modal_projection_mode": int(top_modal_projection_mode),
        "top_modal_projection_share": float(top_modal_projection_share),
        "cumulative_projection_share_first_3": float(cumulative_projection_share_first_3),
        "residual_largest_group": residual_largest_group,
        "residual_largest_group_share": float(residual_largest_group_share),
        "residual_largest_element_arm": residual_largest_element_arm,
        "residual_largest_element_index": int(residual_largest_element_index),
        "residual_largest_element_share": float(residual_largest_element_share),
        "q_joint_translation_mismatch": float(kinematic["q_joint_translation_mismatch"]),
        "q_joint_rotation_mismatch": float(kinematic["q_joint_rotation_mismatch"]),
        "q_external_clamp_max_abs": float(kinematic["q_external_clamp_max_abs"]),
        "modal_projection_csv": str(modal_projection_csv) if bool(args.localize_residual) else "",
        "residual_groups_csv": str(residual_groups_csv) if bool(args.localize_residual) else "",
        "element_residual_energy_csv": str(element_residual_energy_csv) if bool(args.localize_residual) else "",
        "bending_basis_pairing_scan_csv": str(bending_basis_scan_csv)
        if bool(args.scan_bending_basis_pairings)
        else "",
        "full_basis_column_audit_csv": str(full_basis_column_audit_csv)
        if bool(args.scan_bending_basis_pairings)
        else "",
        "bending_basis_canonical_relative_residual": float(bending_basis_canonical_relative_residual),
        "bending_basis_best_variant": bending_basis_best_variant,
        "bending_basis_best_relative_residual": float(bending_basis_best_relative_residual),
        "bending_basis_best_preserves_matrix_consistency": bending_basis_best_preserves_matrix_consistency,
        "bending_basis_pairing_mismatch_flag": bending_basis_pairing_mismatch_flag,
        "axial_convention_scan_csv": str(axial_convention_scan_csv) if bool(args.scan_axial_conventions) else "",
        "axial_convention_baseline_relative_residual": float(axial_convention_baseline_relative_residual),
        "axial_convention_baseline_rayleigh_lambda": float(axial_convention_baseline_rayleigh_lambda),
        "axial_convention_best_variant": axial_convention_best_variant,
        "axial_convention_best_relative_residual": float(axial_convention_best_relative_residual),
        "axial_convention_best_rayleigh_lambda": float(axial_convention_best_rayleigh_lambda),
        "axial_convention_mismatch_flag": axial_convention_mismatch_flag,
        "axial_scale_scan_csv": str(axial_scale_scan_csv) if bool(args.scan_axial_scales) else "",
        "axial_scale_baseline_relative_residual": float(axial_scale_baseline_relative_residual),
        "axial_scale_baseline_rayleigh_lambda": float(axial_scale_baseline_rayleigh_lambda),
        "axial_scale_best_residual_variant": axial_scale_best_residual_variant,
        "axial_scale_best_residual_alpha_left": float(axial_scale_best_residual_alpha_left),
        "axial_scale_best_residual_alpha_right": float(axial_scale_best_residual_alpha_right),
        "axial_scale_best_residual_relative_residual": float(axial_scale_best_residual_relative_residual),
        "axial_scale_best_residual_rayleigh_lambda": float(axial_scale_best_residual_rayleigh_lambda),
        "axial_scale_best_rayleigh_variant": axial_scale_best_rayleigh_variant,
        "axial_scale_best_rayleigh_alpha_left": float(axial_scale_best_rayleigh_alpha_left),
        "axial_scale_best_rayleigh_alpha_right": float(axial_scale_best_rayleigh_alpha_right),
        "axial_scale_best_rayleigh_relative_residual": float(axial_scale_best_rayleigh_relative_residual),
        "axial_scale_best_rayleigh_lambda": float(axial_scale_best_rayleigh_lambda),
        "axial_scale_best_mass_mac_variant": axial_scale_best_mass_mac_variant,
        "axial_scale_best_mass_mac": float(axial_scale_best_mass_mac),
        "axial_scale_mismatch_flag": axial_scale_mismatch_flag,
        "right_transform_convention_scan_csv": str(right_transform_scan_csv)
        if bool(args.scan_right_transform_conventions)
        else "",
        "right_transform_sanity_csv": str(right_transform_sanity_csv)
        if bool(args.scan_right_transform_conventions)
        else "",
        "right_transform_baseline_relative_residual": float(right_transform_baseline_relative_residual),
        "right_transform_baseline_rayleigh_lambda": float(right_transform_baseline_rayleigh_lambda),
        "right_transform_best_residual_variant": right_transform_best_residual_variant,
        "right_transform_best_residual_relative_residual": float(right_transform_best_residual_relative_residual),
        "right_transform_best_residual_rayleigh_lambda": float(right_transform_best_residual_rayleigh_lambda),
        "right_transform_best_rayleigh_variant": right_transform_best_rayleigh_variant,
        "right_transform_best_rayleigh_relative_residual": float(right_transform_best_rayleigh_relative_residual),
        "right_transform_best_rayleigh_lambda": float(right_transform_best_rayleigh_lambda),
        "right_transform_best_mass_mac_variant": right_transform_best_mass_mac_variant,
        "right_transform_best_mass_mac": float(right_transform_best_mass_mac),
        "right_transform_mismatch_flag": right_transform_mismatch_flag,
        "fem_transform_inference": fem_transform_inference_text(),
        "frequency_scaling_audit_csv": str(frequency_scaling_audit_csv) if bool(args.audit_frequency_scaling) else "",
        "frequency_scaling_best_analytic_factor": frequency_scaling_best_analytic_factor,
        "frequency_scaling_best_analytic_relative_residual": float(
            frequency_scaling_best_analytic_relative_residual
        ),
        "frequency_scaling_best_fem_factor": frequency_scaling_best_fem_factor,
        "frequency_scaling_best_fem_relative_residual": float(frequency_scaling_best_fem_relative_residual),
        "summary_csv": str(summary_csv),
        "debug_nodes_csv": str(nodes_csv) if bool(args.save_debug_csv) else "",
    }
    row["conclusion_flag"] = conclusion_flag(row)

    write_summary_csv(summary_csv, row)
    if bool(args.save_debug_csv):
        debug_rows = build_node_debug_rows(
            q_analytic_aligned=q_analytic_aligned,
            q_fem=q_fem,
            params=params,
            mu=float(args.mu),
            beta_deg=float(args.beta),
        )
        write_node_debug_csv(nodes_csv, debug_rows)
    if theta_scan_rows:
        print_theta_convention_scan_top(
            theta_scan_rows,
            baseline_relative_residual=float(row["relative_residual"]),
        )
    if bending_basis_rows:
        print_bending_basis_pairing_scan_top(bending_basis_rows)
    if axial_convention_rows:
        print_axial_convention_scan_top(axial_convention_rows)
    if axial_scale_rows:
        print_axial_scale_scan_top(axial_scale_rows)
    if right_transform_rows:
        print_right_transform_convention_scan_top(right_transform_rows, right_transform_sanity_data)
    if frequency_scaling_rows:
        print_frequency_scaling_audit(frequency_scaling_rows)
    print_summary(row)
    return row


if __name__ == "__main__":
    main()
