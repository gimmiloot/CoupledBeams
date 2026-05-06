from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import (  # noqa: E402
    assemble_clamped_coupled_matrix,
    det_clamped_coupled,
)
from my_project.analytic.solvers import (  # noqa: E402
    find_first_n_roots,
    root_by_min_abs_det,
)
from my_project.fem import python_fem as fem  # noqa: E402
from scripts.lib.tracked_bending_descendant_shapes import (  # noqa: E402
    arm_energy_diagnostics,
    collect_single_branch_shape,
    compute_local_components,
    radius_from_epsilon,
)


DEFAULT_BRANCH_NUMBER = 5
DEFAULT_BETA = 30.0
DEFAULT_MU = 0.0
DEFAULT_EPSILON = 0.0025
DEFAULT_L_TOTAL = 2.0
DEFAULT_NUM_ROOTS = 20
DEFAULT_ROOT_WINDOW = 1.0
RIGHT_COORDINATE_CHOICES = ("external-to-joint", "joint-to-external")
DEFAULT_RIGHT_COORDINATE = "external-to-joint"
NEAR_ZERO_NORM = 1e-12
RESULTS_DIR = REPO_ROOT / "results"
COMPONENT_KEYS = ("u_left", "w_left", "u_right", "w_right")
SIGN_CHOICES = (-1.0, 1.0)
ORIENTATION_WARNING_RELATIVE_FACTOR = 0.75
ORIENTATION_WARNING_ABSOLUTE_GAIN = 0.20
ENDPOINT_CONSISTENCY_TOL = 1e-9
ROW_LABELS = (
    "kinematic_transverse_compatibility",
    "kinematic_longitudinal_compatibility",
    "slope_compatibility",
    "moment_equilibrium",
    "transverse_global_force_equilibrium",
    "axial_global_force_equilibrium",
)
ENDPOINT_CONSISTENCY_FIELDNAMES = [
    "check_type",
    "basis_index",
    "row_index",
    "row_label",
    "quantity",
    "convention",
    "matrix_residual",
    "field_residual",
    "matrix_column_value",
    "field_column_value",
    "abs_difference",
    "relative_difference",
    "left_external_w",
    "left_external_slope",
    "left_external_u",
    "right_external_w",
    "right_external_slope",
    "right_external_u",
    "fem_row1_kinematic_residual",
    "fem_row2_kinematic_residual",
    "fem_slope_residual_if_available",
    "fem_joint_kinematic_max_abs_residual",
]


def branch_id_from_number(branch_number: int) -> str:
    return f"bending_desc_{branch_number:02d}"


def branch_number_from_id(branch_id: str) -> int | None:
    prefix = "bending_desc_"
    if branch_id.startswith(prefix):
        suffix = branch_id.removeprefix(prefix)
        if suffix.isdecimal():
            return int(suffix)
    return None


def option_was_provided(argv: Sequence[str], option_name: str) -> bool:
    return any(arg == option_name or arg.startswith(f"{option_name}=") for arg in argv)


def filename_number_token(value: float) -> str:
    return f"{float(value):g}".replace("-", "m").replace("+", "").replace(".", "p")


def default_output_prefix(branch_id: str, beta_deg: float, mu_value: float, epsilon: float) -> Path:
    return RESULTS_DIR / (
        f"analytic_fem_shape_compare_beta{filename_number_token(beta_deg)}_{branch_id}"
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
    branch_number_explicit = option_was_provided(raw_args, "--branch-number")
    branch_id_explicit = option_was_provided(raw_args, "--branch-id")

    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Compare determinant-nullspace analytic local components against one tracked FEM "
            "bending-descendant shape. This is diagnostic postprocessing only."
        ),
    )
    parser.add_argument(
        "--branch-number",
        type=int,
        default=DEFAULT_BRANCH_NUMBER,
        help="Descendant branch number mapped to bending_desc_NN. Default: 5.",
    )
    parser.add_argument(
        "--branch-id",
        default=None,
        help="Tracked branch id alternative to --branch-number, e.g. bending_desc_05.",
    )
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA, help="Coupling angle in degrees.")
    parser.add_argument("--mu", type=float, default=DEFAULT_MU, help="Length-asymmetry parameter.")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON, help="Slenderness/coupling epsilon.")
    parser.add_argument(
        "--l-total",
        type=float,
        default=DEFAULT_L_TOTAL,
        help="Total beam length used only for epsilon -> FEM radius conversion.",
    )
    parser.add_argument(
        "--num-roots",
        type=int,
        default=DEFAULT_NUM_ROOTS,
        help="Number of analytic roots to acquire before nearest-root matching.",
    )
    parser.add_argument(
        "--root-window",
        type=float,
        default=DEFAULT_ROOT_WINDOW,
        help="Half-width of the local determinant scan around the FEM Lambda.",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help="Output path prefix. Relative paths are resolved from the repo root.",
    )
    parser.add_argument(
        "--right-coordinate",
        choices=RIGHT_COORDINATE_CHOICES,
        default=DEFAULT_RIGHT_COORDINATE,
        help=(
            "Comparison convention for right-arm arrays. Default external-to-joint gives "
            "s/L = 0 at the right clamp and s/L = 1 at the joint."
        ),
    )
    parser.add_argument(
        "--use-best-orientation",
        action="store_true",
        help=(
            "Use the best orientation/sign variant for the main overlay and sample CSV. "
            "By default the best variant is reported only."
        ),
    )
    parser.add_argument(
        "--check-endpoint-consistency",
        action="store_true",
        help=(
            "Write endpoint/matrix consistency diagnostics for the analytic null vector, "
            "basis columns, external clamps, and FEM joint kinematics."
        ),
    )
    parser.add_argument(
        "--compare-energies",
        action="store_true",
        help="Write analytic-vs-FEM arm-wise axial/bending energy diagnostics.",
    )
    parser.add_argument(
        "--check-fem-direct-energy",
        action="store_true",
        help="Write a direct FEM element-stiffness energy check against arm_energy_diagnostics.",
    )
    args = parser.parse_args(raw_args)

    if branch_number_explicit and branch_id_explicit:
        parser.error("--branch-number and --branch-id are alternatives; provide only one.")
    if branch_number_explicit and not branch_id_explicit:
        args.branch_id = None
    if args.branch_number is not None and args.branch_number <= 0:
        parser.error("--branch-number must be positive.")
    if args.beta < 0.0:
        parser.error("--beta must be non-negative; tracked beta continuation starts at beta = 0.")
    if args.epsilon <= 0.0:
        parser.error("--epsilon must be positive.")
    if args.l_total <= 0.0:
        parser.error("--l-total must be positive.")
    if args.num_roots <= 0:
        parser.error("--num-roots must be positive.")
    if args.root_window <= 0.0:
        parser.error("--root-window must be positive.")
    return args


def determinant_value(Lambda: float, beta_rad: float, mu_value: float, epsilon: float) -> float:
    return det_clamped_coupled(float(Lambda), float(beta_rad), float(mu_value), float(epsilon))


def root_bisect(
    beta_rad: float,
    mu_value: float,
    epsilon: float,
    left: float,
    right: float,
    *,
    iterations: int = 80,
) -> float:
    a = float(left)
    b = float(right)
    fa = determinant_value(a, beta_rad, mu_value, epsilon)
    fb = determinant_value(b, beta_rad, mu_value, epsilon)
    if not (np.isfinite(fa) and np.isfinite(fb)):
        raise ValueError("non-finite determinant value in local bracket")
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if np.sign(fa) == np.sign(fb):
        raise ValueError("local bracket does not straddle a determinant sign change")

    for _ in range(iterations):
        mid = 0.5 * (a + b)
        fm = determinant_value(mid, beta_rad, mu_value, epsilon)
        if fm == 0.0:
            return mid
        if np.sign(fa) == np.sign(fm):
            a = mid
            fa = fm
        else:
            b = mid
    return 0.5 * (a + b)


def unique_sorted_roots(values: Sequence[float], *, tolerance: float = 1e-7) -> list[float]:
    roots = sorted(float(value) for value in values if np.isfinite(value) and float(value) > 0.0)
    unique: list[float] = []
    for root in roots:
        if not unique or abs(root - unique[-1]) > tolerance:
            unique.append(root)
    return unique


def local_scan_roots(
    beta_rad: float,
    mu_value: float,
    epsilon: float,
    center: float,
    half_width: float,
) -> list[float]:
    left = max(0.05, float(center) - float(half_width))
    right = max(left + 1e-6, float(center) + float(half_width))
    width = right - left
    n_points = int(min(4001, max(401, np.ceil(width / 0.002) + 1)))
    grid = np.linspace(left, right, n_points)
    vals = np.array([determinant_value(x, beta_rad, mu_value, epsilon) for x in grid], dtype=float)

    roots: list[float] = []
    for idx in range(len(grid) - 1):
        fa = vals[idx]
        fb = vals[idx + 1]
        if not (np.isfinite(fa) and np.isfinite(fb)):
            continue
        if fa == 0.0:
            roots.append(float(grid[idx]))
        elif np.sign(fa) != np.sign(fb):
            roots.append(root_bisect(beta_rad, mu_value, epsilon, float(grid[idx]), float(grid[idx + 1])))
    if vals[-1] == 0.0:
        roots.append(float(grid[-1]))
    return unique_sorted_roots(roots)


def global_root_index(root: float, global_roots: Sequence[float], *, tolerance: float = 1e-6) -> int:
    for idx, candidate in enumerate(global_roots, start=1):
        if abs(float(root) - float(candidate)) <= tolerance:
            return idx
    return 0


def find_analytic_root_near_fem(
    fem_lambda: float,
    *,
    beta_rad: float,
    mu_value: float,
    epsilon: float,
    num_roots: int,
    root_window: float,
) -> dict[str, float | int | str | list[float]]:
    raw_global_roots = find_first_n_roots(beta_rad, mu_value, epsilon, num_roots)
    global_roots = unique_sorted_roots(raw_global_roots)
    local_roots = local_scan_roots(
        beta_rad=beta_rad,
        mu_value=mu_value,
        epsilon=epsilon,
        center=fem_lambda,
        half_width=root_window,
    )

    candidates: list[dict[str, float | int | str]] = []
    for idx, root in enumerate(global_roots, start=1):
        candidates.append({"root": root, "source": "global_first_n", "analytic_root_index": idx})
    for root in local_roots:
        candidates.append(
            {
                "root": root,
                "source": "local_sign_scan",
                "analytic_root_index": global_root_index(root, global_roots),
            }
        )

    if not candidates:
        left = max(0.05, fem_lambda - root_window)
        right = max(left + 1e-6, fem_lambda + root_window)
        det_here = lambda x: determinant_value(x, beta_rad, mu_value, epsilon)
        root, det_abs, err = root_by_min_abs_det(det_here, left, right, tol_lambda=1e-12)
        candidates.append(
            {
                "root": float(root),
                "source": "local_min_abs_det",
                "analytic_root_index": 0,
                "local_min_abs_det": float(det_abs),
                "local_min_abs_det_lambda_error_bound": float(err),
            }
        )

    selected = min(candidates, key=lambda item: abs(float(item["root"]) - float(fem_lambda)))
    analytic_lambda = float(selected["root"])
    return {
        "analytic_lambda": analytic_lambda,
        "analytic_root_index": int(selected["analytic_root_index"]),
        "analytic_root_source": str(selected["source"]),
        "analytic_root_abs_det": abs(determinant_value(analytic_lambda, beta_rad, mu_value, epsilon)),
        "global_roots_found": len(global_roots),
        "local_roots_found": len(local_roots),
        "global_roots": global_roots,
        "local_roots": local_roots,
    }


def analytic_null_vector(
    Lambda: float,
    *,
    beta_rad: float,
    mu_value: float,
    epsilon: float,
) -> tuple[np.ndarray, float, float]:
    matrix = assemble_clamped_coupled_matrix(Lambda, beta_rad, mu_value, epsilon)
    _, singular_values, vh = np.linalg.svd(matrix)
    coeff = vh[-1, :].astype(float)
    smallest = float(singular_values[-1])
    ratio = (
        float(singular_values[-1] / singular_values[-2])
        if len(singular_values) >= 2 and abs(float(singular_values[-2])) > NEAR_ZERO_NORM
        else np.nan
    )
    return coeff, smallest, ratio


def reconstruct_analytic_components(
    Lambda: float,
    *,
    mu_value: float,
    epsilon: float,
    coeff: np.ndarray,
    s_norm: np.ndarray,
    right_coordinate: str,
) -> dict[str, np.ndarray]:
    A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]
    xi = np.asarray(s_norm, dtype=float)

    z1 = Lambda * (1.0 - mu_value) * xi
    theta1_arg = epsilon * Lambda**2 * (1.0 - mu_value) * xi
    w_left = A1 * (np.cos(z1) - np.cosh(z1)) + B1 * (np.sin(z1) - np.sinh(z1))
    u_left = P1 * np.sin(theta1_arg)

    # Determinant unknown order is A1, B1, A2, B2, P1, P2. For the right arm,
    # this diagnostic first evaluates the determinant-side external-to-joint
    # convention with a negative local argument, then reverses only if the user
    # explicitly asks for joint-to-external output.
    z2 = -Lambda * (1.0 + mu_value) * xi
    theta2_arg = -epsilon * Lambda**2 * (1.0 + mu_value) * xi
    w_right = A2 * (np.cos(z2) - np.cosh(z2)) + B2 * (np.sin(z2) - np.sinh(z2))
    u_right = P2 * np.sin(theta2_arg)

    if right_coordinate == "joint-to-external":
        u_right = u_right[::-1]
        w_right = w_right[::-1]

    return {
        "u_left": u_left,
        "w_left": w_left,
        "u_right": u_right,
        "w_right": w_right,
    }


def fem_components_for_comparison(
    shape_case: dict[str, object],
    *,
    right_coordinate: str,
) -> tuple[dict[str, np.ndarray], bool]:
    u_left = np.asarray(shape_case["u_left_local"], dtype=float)
    w_left = np.asarray(shape_case["w_left_local"], dtype=float)
    u_right = np.asarray(shape_case["u_right_local"], dtype=float)
    w_right = np.asarray(shape_case["w_right_local"], dtype=float)

    # Current FEM node ids on the right arm run from the joint to the right
    # clamp. The default comparison convention is s/L = 0 at each external
    # clamp and s/L = 1 at the joint, so right-arm local components are reversed.
    right_reversed = right_coordinate == "external-to-joint"
    if right_reversed:
        u_right = u_right[::-1]
        w_right = w_right[::-1]

    return (
        {
            "u_left": u_left,
            "w_left": w_left,
            "u_right": u_right,
            "w_right": w_right,
        },
        right_reversed,
    )


def concatenate_components(components: dict[str, np.ndarray]) -> np.ndarray:
    return np.concatenate([np.asarray(components[key], dtype=float) for key in COMPONENT_KEYS])


def normalize_components(
    fem_components: dict[str, np.ndarray],
    analytic_components: dict[str, np.ndarray],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, float]]:
    fem_vector = concatenate_components(fem_components)
    analytic_vector = concatenate_components(analytic_components)
    fem_scale = max(float(np.max(np.abs(fem_vector))), NEAR_ZERO_NORM)
    analytic_scale = max(float(np.max(np.abs(analytic_vector))), NEAR_ZERO_NORM)

    fem_normed = {key: np.asarray(value, dtype=float) / fem_scale for key, value in fem_components.items()}
    analytic_normed = {key: np.asarray(value, dtype=float) / analytic_scale for key, value in analytic_components.items()}

    fem_vector_normed = concatenate_components(fem_normed)
    analytic_vector_normed = concatenate_components(analytic_normed)
    dot_before = float(np.dot(analytic_vector_normed, fem_vector_normed))
    analytic_sign = -1.0 if dot_before < 0.0 else 1.0
    if analytic_sign < 0.0:
        analytic_normed = {key: -value for key, value in analytic_normed.items()}

    return (
        fem_normed,
        analytic_normed,
        {
            "fem_normalization_max_abs": fem_scale,
            "analytic_normalization_max_abs": analytic_scale,
            "shape_dot_before_sign_alignment": dot_before,
            "analytic_sign_applied": analytic_sign,
        },
    )


def safe_relative_l2(diff: np.ndarray, reference: np.ndarray) -> float:
    denominator = float(np.linalg.norm(reference))
    return float(np.linalg.norm(diff) / denominator) if denominator > NEAR_ZERO_NORM else np.nan


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if abs(float(denominator)) > NEAR_ZERO_NORM else np.nan


def comparison_errors(
    fem_components: dict[str, np.ndarray],
    analytic_components: dict[str, np.ndarray],
) -> dict[str, float]:
    fem_vector = concatenate_components(fem_components)
    analytic_vector = concatenate_components(analytic_components)
    diff_vector = analytic_vector - fem_vector
    dot_denominator = float(np.linalg.norm(analytic_vector) * np.linalg.norm(fem_vector))
    dot_after = float(np.dot(analytic_vector, fem_vector) / dot_denominator) if dot_denominator > NEAR_ZERO_NORM else np.nan

    errors: dict[str, float] = {
        "shape_dot_after_sign_alignment": dot_after,
        "relative_L2_error_all_components": safe_relative_l2(diff_vector, fem_vector),
        "max_abs_error_all_components": float(np.max(np.abs(diff_vector))),
    }
    for component in COMPONENT_KEYS:
        diff = analytic_components[component] - fem_components[component]
        errors[f"rel_l2_{component}"] = safe_relative_l2(diff, fem_components[component])
        errors[f"max_abs_error_{component}"] = float(np.max(np.abs(diff)))

    fem_u_right = float(np.max(np.abs(fem_components["u_right"])))
    fem_w_right = float(np.max(np.abs(fem_components["w_right"])))
    analytic_u_right = float(np.max(np.abs(analytic_components["u_right"])))
    analytic_w_right = float(np.max(np.abs(analytic_components["w_right"])))
    errors.update(
        {
            "fem_max_abs_u_right_normalized": fem_u_right,
            "fem_max_abs_w_right_normalized": fem_w_right,
            "analytic_max_abs_u_right_normalized": analytic_u_right,
            "analytic_max_abs_w_right_normalized": analytic_w_right,
            "fem_right_axial_to_right_transverse": safe_ratio(fem_u_right, fem_w_right),
            "analytic_right_axial_to_right_transverse": safe_ratio(analytic_u_right, analytic_w_right),
        }
    )
    return errors


def symmetric_relative_difference(left: float, right: float) -> float:
    difference = abs(float(left) - float(right))
    scale = max(abs(float(left)), abs(float(right)))
    return difference / scale if scale > NEAR_ZERO_NORM else 0.0


def analytic_endpoint_quantities(
    Lambda: float,
    *,
    beta_rad: float,
    mu_value: float,
    epsilon: float,
    coeff: np.ndarray,
) -> dict[str, float]:
    """Endpoint algebra in the same row-level scaling as the determinant matrix."""
    A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]
    x1 = float(Lambda) * (1.0 - float(mu_value))
    x2 = float(Lambda) * (1.0 + float(mu_value))
    th1 = float(epsilon) * float(Lambda) ** 2 * (1.0 - float(mu_value))
    th2 = float(epsilon) * float(Lambda) ** 2 * (1.0 + float(mu_value))

    components = reconstruct_analytic_components(
        Lambda,
        mu_value=mu_value,
        epsilon=epsilon,
        coeff=coeff,
        s_norm=np.array([0.0, 1.0], dtype=float),
        right_coordinate=DEFAULT_RIGHT_COORDINATE,
    )

    def slope_value(A: float, B: float, z: float) -> float:
        return A * (-np.sin(z) - np.sinh(z)) + B * (np.cos(z) - np.cosh(z))

    def moment_value(A: float, B: float, z: float) -> float:
        return A * (-np.cos(z) - np.cosh(z)) + B * (-np.sin(z) - np.sinh(z))

    def shear_row_value(A: float, B: float, z: float) -> float:
        d3w_dz3 = A * (np.sin(z) - np.sinh(z)) + B * (-np.cos(z) - np.cosh(z))
        return -float(epsilon) * float(Lambda) * d3w_dz3

    return {
        "u1_joint": float(components["u_left"][-1]),
        "w1_joint": float(components["w_left"][-1]),
        "slope1_joint": float(slope_value(A1, B1, x1)),
        "moment1_joint": float(moment_value(A1, B1, x1)),
        "shear1_joint": float(shear_row_value(A1, B1, x1)),
        "axial_force1_joint": float(P1 * np.cos(th1)),
        "u2_joint": float(components["u_right"][-1]),
        "w2_joint": float(components["w_right"][-1]),
        "slope2_joint": float(slope_value(A2, B2, -x2)),
        "moment2_joint": float(moment_value(A2, B2, -x2)),
        "shear2_joint": float(shear_row_value(A2, B2, -x2)),
        "axial_force2_joint": float(P2 * np.cos(-th2)),
        "left_external_w": float(components["w_left"][0]),
        "left_external_slope": float(slope_value(A1, B1, 0.0)),
        "left_external_u": float(components["u_left"][0]),
        "right_external_w": float(components["w_right"][0]),
        "right_external_slope": float(slope_value(A2, B2, 0.0)),
        "right_external_u": float(components["u_right"][0]),
    }


def field_residuals_from_endpoint_quantities(
    endpoint: dict[str, float],
    *,
    beta_rad: float,
) -> np.ndarray:
    cb = float(np.cos(beta_rad))
    sb = float(np.sin(beta_rad))
    return np.array(
        [
            endpoint["w1_joint"] - endpoint["w2_joint"] * cb - endpoint["u2_joint"] * sb,
            endpoint["u1_joint"] + endpoint["w2_joint"] * sb - endpoint["u2_joint"] * cb,
            endpoint["slope1_joint"] - endpoint["slope2_joint"],
            endpoint["moment1_joint"] - endpoint["moment2_joint"],
            endpoint["shear1_joint"] - endpoint["shear2_joint"] * cb - endpoint["axial_force2_joint"] * sb,
            endpoint["axial_force1_joint"] + endpoint["shear2_joint"] * sb - endpoint["axial_force2_joint"] * cb,
        ],
        dtype=float,
    )


def analytic_field_residuals(
    Lambda: float,
    *,
    beta_rad: float,
    mu_value: float,
    epsilon: float,
    coeff: np.ndarray,
) -> np.ndarray:
    endpoint = analytic_endpoint_quantities(
        Lambda,
        beta_rad=beta_rad,
        mu_value=mu_value,
        epsilon=epsilon,
        coeff=coeff,
    )
    return field_residuals_from_endpoint_quantities(endpoint, beta_rad=beta_rad)


def analytic_external_clamp_quantities(
    Lambda: float,
    *,
    beta_rad: float,
    mu_value: float,
    epsilon: float,
    coeff: np.ndarray,
) -> dict[str, float]:
    endpoint = analytic_endpoint_quantities(
        Lambda,
        beta_rad=beta_rad,
        mu_value=mu_value,
        epsilon=epsilon,
        coeff=coeff,
    )
    return {
        "left_external_w": endpoint["left_external_w"],
        "left_external_slope": endpoint["left_external_slope"],
        "left_external_u": endpoint["left_external_u"],
        "right_external_w": endpoint["right_external_w"],
        "right_external_slope": endpoint["right_external_slope"],
        "right_external_u": endpoint["right_external_u"],
    }


def fem_joint_kinematic_rows(
    shape_case: dict[str, object],
    *,
    beta_deg: float,
) -> list[dict[str, float | str]]:
    u_raw = np.asarray(shape_case["u_raw"], dtype=float)
    v_raw = np.asarray(shape_case["v_raw"], dtype=float)
    theta_raw = np.asarray(shape_case["theta_raw"], dtype=float)
    local = compute_local_components(u_global=u_raw, v_global=v_raw, beta_deg=beta_deg)
    joint = fem.N_ELEM
    cb = float(np.cos(np.deg2rad(beta_deg)))
    sb = float(np.sin(np.deg2rad(beta_deg)))
    u_left_joint = float(local["u_left_local"][joint])
    w_left_joint = float(local["w_left_local"][joint])
    theta_left_joint = float(theta_raw[joint])

    rows: list[dict[str, float | str]] = []
    for convention in RIGHT_COORDINATE_CHOICES:
        u_right = np.asarray(local["u_right_local"], dtype=float)
        w_right = np.asarray(local["w_right_local"], dtype=float)
        theta_right = np.asarray(theta_raw[joint:], dtype=float)
        slope_multiplier = 1.0
        joint_index = 0
        if convention == "external-to-joint":
            u_right = u_right[::-1]
            w_right = w_right[::-1]
            theta_right = theta_right[::-1]
            slope_multiplier = -1.0
            joint_index = -1

        u_right_joint = float(u_right[joint_index])
        w_right_joint = float(w_right[joint_index])
        theta_right_joint = float(theta_right[joint_index])
        row1 = w_left_joint - (w_right_joint * cb + u_right_joint * sb)
        row2 = u_left_joint - (-w_right_joint * sb + u_right_joint * cb)
        slope_residual = theta_left_joint - slope_multiplier * theta_right_joint
        rows.append(
            {
                "check_type": "fem_joint_kinematic",
                "convention": convention,
                "fem_row1_kinematic_residual": float(row1),
                "fem_row2_kinematic_residual": float(row2),
                "fem_slope_residual_if_available": float(slope_residual),
                "fem_joint_kinematic_max_abs_residual": float(
                    max(abs(row1), abs(row2), abs(slope_residual))
                ),
            }
        )
    return rows


def endpoint_consistency_rows(
    *,
    Lambda: float,
    beta_rad: float,
    mu_value: float,
    epsilon: float,
    coeff: np.ndarray,
    shape_case: dict[str, object],
) -> tuple[list[dict[str, float | int | str]], dict[str, object]]:
    matrix = assemble_clamped_coupled_matrix(Lambda, beta_rad, mu_value, epsilon)
    matrix_residual = matrix @ coeff
    field_residual = analytic_field_residuals(
        Lambda,
        beta_rad=beta_rad,
        mu_value=mu_value,
        epsilon=epsilon,
        coeff=coeff,
    )
    rows: list[dict[str, float | int | str]] = []
    null_max_difference = 0.0
    for row_index, row_label in enumerate(ROW_LABELS):
        matrix_value = float(matrix_residual[row_index])
        field_value = float(field_residual[row_index])
        abs_difference = abs(matrix_value - field_value)
        null_max_difference = max(null_max_difference, abs_difference)
        rows.append(
            {
                "check_type": "analytic_null_row",
                "row_index": row_index,
                "row_label": row_label,
                "matrix_residual": matrix_value,
                "field_residual": field_value,
                "abs_difference": float(abs_difference),
                "relative_difference": symmetric_relative_difference(matrix_value, field_value),
            }
        )

    basis_max_difference = 0.0
    worst_basis_index = 0
    worst_basis_row_index = 0
    worst_basis_row_label = ROW_LABELS[0]
    for basis_index in range(matrix.shape[1]):
        basis_coeff = np.zeros(matrix.shape[1], dtype=float)
        basis_coeff[basis_index] = 1.0
        field_column = analytic_field_residuals(
            Lambda,
            beta_rad=beta_rad,
            mu_value=mu_value,
            epsilon=epsilon,
            coeff=basis_coeff,
        )
        matrix_column = matrix[:, basis_index]
        for row_index, row_label in enumerate(ROW_LABELS):
            matrix_value = float(matrix_column[row_index])
            field_value = float(field_column[row_index])
            abs_difference = abs(matrix_value - field_value)
            if abs_difference > basis_max_difference:
                basis_max_difference = float(abs_difference)
                worst_basis_index = basis_index
                worst_basis_row_index = row_index
                worst_basis_row_label = row_label
            rows.append(
                {
                    "check_type": "analytic_basis_column",
                    "basis_index": basis_index,
                    "row_index": row_index,
                    "row_label": row_label,
                    "matrix_column_value": matrix_value,
                    "field_column_value": field_value,
                    "abs_difference": float(abs_difference),
                    "relative_difference": symmetric_relative_difference(matrix_value, field_value),
                }
            )

    external = analytic_external_clamp_quantities(
        Lambda,
        beta_rad=beta_rad,
        mu_value=mu_value,
        epsilon=epsilon,
        coeff=coeff,
    )
    external_max_residual = max(abs(float(value)) for value in external.values())
    rows.append(
        {
            "check_type": "analytic_external_clamp",
            **{key: float(value) for key, value in external.items()},
        }
    )

    fem_rows = fem_joint_kinematic_rows(shape_case, beta_deg=float(np.rad2deg(beta_rad)))
    rows.extend(fem_rows)

    fem_consistent_conventions = [
        str(row["convention"])
        for row in fem_rows
        if float(row["fem_joint_kinematic_max_abs_residual"]) <= ENDPOINT_CONSISTENCY_TOL
    ]
    summary: dict[str, object] = {
        "null_max_abs_difference": float(null_max_difference),
        "basis_max_abs_difference": float(basis_max_difference),
        "worst_basis_index": int(worst_basis_index),
        "worst_basis_row_index": int(worst_basis_row_index),
        "worst_basis_row_label": worst_basis_row_label,
        "external_max_residual": float(external_max_residual),
        "fem_rows": fem_rows,
        "fem_consistent_conventions": fem_consistent_conventions,
        "analytic_matrix_consistent": bool(
            max(null_max_difference, basis_max_difference) <= ENDPOINT_CONSISTENCY_TOL
        ),
        "analytic_external_clamp_consistent": bool(external_max_residual <= ENDPOINT_CONSISTENCY_TOL),
    }
    return rows, summary


def write_endpoint_consistency_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ENDPOINT_CONSISTENCY_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def endpoint_consistency_summary_lines(
    *,
    csv_path: Path,
    summary_path: Path,
    summary: dict[str, object],
) -> list[str]:
    null_max = float(summary["null_max_abs_difference"])
    basis_max = float(summary["basis_max_abs_difference"])
    external_max = float(summary["external_max_residual"])
    worst_basis_index = int(summary["worst_basis_index"])
    worst_basis_row_index = int(summary["worst_basis_row_index"])
    worst_basis_row_label = str(summary["worst_basis_row_label"])
    fem_rows = list(summary["fem_rows"])  # type: ignore[arg-type]
    fem_consistent = [str(value) for value in summary["fem_consistent_conventions"]]  # type: ignore[index]
    analytic_text = (
        "analytic reconstruction is consistent with matrix"
        if bool(summary["analytic_matrix_consistent"])
        else "analytic reconstruction is not consistent with matrix"
    )
    fem_text = (
        "FEM local components are consistent with analytic kinematics under convention(s): "
        + ", ".join(fem_consistent)
        if fem_consistent
        else "FEM local components are not consistent with analytic kinematics under the checked conventions"
    )

    lines = [
        "Endpoint consistency diagnostics",
        f"endpoint_consistency_csv: {csv_path}",
        f"endpoint_consistency_summary_txt: {summary_path}",
        f"max_abs_difference_field_residual_vs_matrix_residual: {null_max:.6e}",
        f"max_basis_column_mismatch: {basis_max:.6e}",
        (
            "worst_basis_row_mismatch: "
            f"basis_index={worst_basis_index}, row_index={worst_basis_row_index}, "
            f"row_label={worst_basis_row_label}"
        ),
        f"analytic_external_clamp_max_residual: {external_max:.6e}",
        "FEM joint kinematic residuals:",
    ]
    for row in fem_rows:
        lines.append(
            f"  {row['convention']}: "
            f"row1={float(row['fem_row1_kinematic_residual']):.6e}, "
            f"row2={float(row['fem_row2_kinematic_residual']):.6e}, "
            f"slope={float(row['fem_slope_residual_if_available']):.6e}, "
            f"max={float(row['fem_joint_kinematic_max_abs_residual']):.6e}"
        )
    lines.append(analytic_text)
    lines.append(fem_text)
    return lines


def run_endpoint_consistency_diagnostics(
    *,
    output_prefix: Path,
    analytic_lambda: float,
    beta_rad: float,
    mu_value: float,
    epsilon: float,
    coeff: np.ndarray,
    shape_case: dict[str, object],
) -> dict[str, object]:
    csv_path = suffixed_path(output_prefix, "_endpoint_consistency.csv")
    summary_path = suffixed_path(output_prefix, "_endpoint_consistency_summary.txt")
    rows, summary = endpoint_consistency_rows(
        Lambda=analytic_lambda,
        beta_rad=beta_rad,
        mu_value=mu_value,
        epsilon=epsilon,
        coeff=coeff,
        shape_case=shape_case,
    )
    write_endpoint_consistency_csv(csv_path, rows)
    lines = endpoint_consistency_summary_lines(csv_path=csv_path, summary_path=summary_path, summary=summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print()
    for line in lines:
        print(line)
    summary["endpoint_consistency_csv"] = str(csv_path)
    summary["endpoint_consistency_summary_txt"] = str(summary_path)
    return summary


def apply_analytic_component_signs(
    analytic_components: dict[str, np.ndarray],
    *,
    sign_u_left: float,
    sign_w_left: float,
    sign_u_right: float,
    sign_w_right: float,
) -> dict[str, np.ndarray]:
    signs = {
        "u_left": float(sign_u_left),
        "w_left": float(sign_w_left),
        "u_right": float(sign_u_right),
        "w_right": float(sign_w_right),
    }
    return {key: signs[key] * np.asarray(value, dtype=float) for key, value in analytic_components.items()}


def orientation_variant_label(row: dict[str, float | int | str]) -> str:
    return (
        f"right={row['right_coordinate']}, "
        f"suL={int(float(row['sign_u_left'])):+d}, "
        f"swL={int(float(row['sign_w_left'])):+d}, "
        f"suR={int(float(row['sign_u_right'])):+d}, "
        f"swR={int(float(row['sign_w_right'])):+d}"
    )


def evaluate_orientation_variant(
    *,
    variant_id: int,
    shape_case: dict[str, object],
    analytic_lambda: float,
    coeff: np.ndarray,
    s_norm: np.ndarray,
    mu_value: float,
    epsilon: float,
    right_coordinate: str,
    sign_u_left: float,
    sign_w_left: float,
    sign_u_right: float,
    sign_w_right: float,
) -> dict[str, object]:
    fem_components, right_fem_arrays_reversed = fem_components_for_comparison(
        shape_case,
        right_coordinate=right_coordinate,
    )
    analytic_base = reconstruct_analytic_components(
        analytic_lambda,
        mu_value=mu_value,
        epsilon=epsilon,
        coeff=coeff,
        s_norm=s_norm,
        right_coordinate=right_coordinate,
    )
    analytic_signed = apply_analytic_component_signs(
        analytic_base,
        sign_u_left=sign_u_left,
        sign_w_left=sign_w_left,
        sign_u_right=sign_u_right,
        sign_w_right=sign_w_right,
    )
    fem_normed, analytic_normed, normalization_info = normalize_components(fem_components, analytic_signed)
    errors = comparison_errors(fem_normed, analytic_normed)
    row: dict[str, float | int | str] = {
        "variant_id": int(variant_id),
        "right_coordinate": right_coordinate,
        "sign_u_left": float(sign_u_left),
        "sign_w_left": float(sign_w_left),
        "sign_u_right": float(sign_u_right),
        "sign_w_right": float(sign_w_right),
        "right_fem_arrays_reversed": str(bool(right_fem_arrays_reversed)),
        "fem_normalization_max_abs": float(normalization_info["fem_normalization_max_abs"]),
        "analytic_normalization_max_abs": float(normalization_info["analytic_normalization_max_abs"]),
        "shape_dot_before_sign_alignment": float(normalization_info["shape_dot_before_sign_alignment"]),
        "analytic_sign_applied": float(normalization_info["analytic_sign_applied"]),
    }
    row.update(errors)
    row["variant_label"] = orientation_variant_label(row)
    return {
        "row": row,
        "fem_components": fem_normed,
        "analytic_components": analytic_normed,
        "normalization_info": normalization_info,
        "errors": errors,
        "right_fem_arrays_reversed": right_fem_arrays_reversed,
    }


def scan_orientation_variants(
    *,
    shape_case: dict[str, object],
    analytic_lambda: float,
    coeff: np.ndarray,
    s_norm: np.ndarray,
    mu_value: float,
    epsilon: float,
) -> tuple[list[dict[str, float | int | str]], dict[str, object]]:
    results: list[dict[str, object]] = []
    variant_id = 0
    for right_coordinate in RIGHT_COORDINATE_CHOICES:
        for sign_u_left in SIGN_CHOICES:
            for sign_w_left in SIGN_CHOICES:
                for sign_u_right in SIGN_CHOICES:
                    for sign_w_right in SIGN_CHOICES:
                        variant_id += 1
                        results.append(
                            evaluate_orientation_variant(
                                variant_id=variant_id,
                                shape_case=shape_case,
                                analytic_lambda=analytic_lambda,
                                coeff=coeff,
                                s_norm=s_norm,
                                mu_value=mu_value,
                                epsilon=epsilon,
                                right_coordinate=right_coordinate,
                                sign_u_left=sign_u_left,
                                sign_w_left=sign_w_left,
                                sign_u_right=sign_u_right,
                                sign_w_right=sign_w_right,
                            )
                        )

    results.sort(key=lambda item: float(item["errors"]["relative_L2_error_all_components"]))  # type: ignore[index]
    return [result["row"] for result in results], results[0]


def write_orientation_scan_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("orientation scan produced no rows")
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def trapezoid_integral(values: np.ndarray, coordinates: np.ndarray) -> float:
    y_values = np.asarray(values, dtype=float)
    x_values = np.asarray(coordinates, dtype=float)
    integrate = getattr(np, "trapezoid", None)
    if integrate is not None:
        return float(integrate(y_values, x_values))
    return float(np.sum(0.5 * (y_values[1:] + y_values[:-1]) * np.diff(x_values)))


def analytic_arm_energy_diagnostics(
    components: dict[str, np.ndarray],
    *,
    mu_value: float,
    epsilon: float,
    s_norm: np.ndarray,
) -> dict[str, float]:
    """Compute diagnostic continuous arm energies from normalized analytic samples."""
    xi = np.asarray(s_norm, dtype=float)
    ea_nd = 1.0 / float(epsilon) ** 2
    length_factors = {
        "left": max(1.0 - float(mu_value), NEAR_ZERO_NORM),
        "right": max(1.0 + float(mu_value), NEAR_ZERO_NORM),
    }
    energies: dict[str, float] = {}
    for arm, length_factor in length_factors.items():
        u_values = np.asarray(components[f"u_{arm}"], dtype=float)
        w_values = np.asarray(components[f"w_{arm}"], dtype=float)
        du_dxi = np.gradient(u_values, xi, edge_order=2)
        dw_dxi = np.gradient(w_values, xi, edge_order=2)
        d2w_dxi2 = np.gradient(dw_dxi, xi, edge_order=2)
        axial_energy = 0.5 * ea_nd / length_factor * trapezoid_integral(du_dxi * du_dxi, xi)
        bending_energy = 0.5 / length_factor**3 * trapezoid_integral(d2w_dxi2 * d2w_dxi2, xi)
        energies[f"{arm}_axial_energy"] = max(float(axial_energy), 0.0)
        energies[f"{arm}_bending_energy"] = max(float(bending_energy), 0.0)

    total_axial = energies["left_axial_energy"] + energies["right_axial_energy"]
    total_bending = energies["left_bending_energy"] + energies["right_bending_energy"]
    total = total_axial + total_bending
    energies.update(
        {
            "total_axial_energy": float(total_axial),
            "total_bending_energy": float(total_bending),
            "total_energy": float(total),
            "axial_energy_fraction": safe_ratio(total_axial, total),
            "right_axial_share": safe_ratio(energies["right_axial_energy"], total_axial),
            "right_bending_share": safe_ratio(energies["right_bending_energy"], total_bending),
        }
    )
    return energies


def scaled_shape_case_for_energy(
    shape_case: dict[str, object],
    *,
    normalization_scale: float,
) -> dict[str, object]:
    scale = max(float(normalization_scale), NEAR_ZERO_NORM)
    scaled = dict(shape_case)
    for key in (
        "u_raw",
        "v_raw",
        "theta_raw",
        "u_left_local",
        "w_left_local",
        "u_right_local",
        "w_right_local",
    ):
        scaled[key] = np.asarray(shape_case[key], dtype=float) / scale
    return scaled


def build_energy_comparison_row(
    *,
    branch_id: str,
    branch_number: int | None,
    beta_deg: float,
    mu_value: float,
    epsilon: float,
    l_total: float,
    radius: float,
    fem_lambda: float,
    analytic_lambda: float,
    shape_case: dict[str, object],
    analytic_components: dict[str, np.ndarray],
    s_norm: np.ndarray,
    normalization_info: dict[str, float],
    selected_variant_row: dict[str, float | int | str],
    output_path: Path,
) -> dict[str, float | int | str]:
    analytic_energy = analytic_arm_energy_diagnostics(
        analytic_components,
        mu_value=mu_value,
        epsilon=epsilon,
        s_norm=s_norm,
    )
    fem_energy = arm_energy_diagnostics(
        scaled_shape_case_for_energy(
            shape_case,
            normalization_scale=float(normalization_info["fem_normalization_max_abs"]),
        ),
        epsilon=epsilon,
    )

    analytic_bending_fraction = safe_ratio(
        analytic_energy["total_bending_energy"],
        analytic_energy["total_energy"],
    )
    fem_bending_fraction = safe_ratio(
        float(fem_energy["total_bending_energy"]),
        float(fem_energy["total_energy"]),
    )
    analytic_right_axial_dominant = bool(float(analytic_energy["right_axial_share"]) > 0.5)
    fem_right_axial_dominant = bool(float(fem_energy["right_axial_share"]) > 0.5)
    right_axial_dominance_agreement = (
        "agree" if analytic_right_axial_dominant == fem_right_axial_dominant else "disagree"
    )

    row: dict[str, float | int | str] = {
        "branch_id": branch_id,
        "branch_number": "" if branch_number is None else int(branch_number),
        "beta": float(beta_deg),
        "mu": float(mu_value),
        "epsilon": float(epsilon),
        "l_total": float(l_total),
        "radius": float(radius),
        "fem_lambda": float(fem_lambda),
        "analytic_lambda": float(analytic_lambda),
        "selected_orientation_variant_id": int(selected_variant_row["variant_id"]),
        "selected_orientation_variant": str(selected_variant_row["variant_label"]),
        "analytic_energy_right_coordinate": str(selected_variant_row["right_coordinate"]),
        "fem_energy_right_coordinate": "joint-to-external/native FEM element order",
        "energy_normalization": (
            "component-comparison max-abs normalization; FEM u/v/theta and local components "
            "scaled by fem_normalization_max_abs; analytic samples scaled and sign-aligned "
            "by normalize_components"
        ),
        "fem_normalization_max_abs": float(normalization_info["fem_normalization_max_abs"]),
        "analytic_normalization_max_abs": float(normalization_info["analytic_normalization_max_abs"]),
        "energy_comparison_csv": str(output_path),
        "analytic_left_axial_energy": float(analytic_energy["left_axial_energy"]),
        "analytic_right_axial_energy": float(analytic_energy["right_axial_energy"]),
        "analytic_left_bending_energy": float(analytic_energy["left_bending_energy"]),
        "analytic_right_bending_energy": float(analytic_energy["right_bending_energy"]),
        "analytic_total_axial_energy": float(analytic_energy["total_axial_energy"]),
        "analytic_total_bending_energy": float(analytic_energy["total_bending_energy"]),
        "analytic_total_energy": float(analytic_energy["total_energy"]),
        "analytic_axial_energy_fraction": float(analytic_energy["axial_energy_fraction"]),
        "analytic_bending_energy_fraction": float(analytic_bending_fraction),
        "analytic_right_axial_share": float(analytic_energy["right_axial_share"]),
        "analytic_right_bending_share": float(analytic_energy["right_bending_share"]),
        "fem_left_axial_energy": float(fem_energy["left_axial_energy"]),
        "fem_right_axial_energy": float(fem_energy["right_axial_energy"]),
        "fem_left_bending_energy": float(fem_energy["left_bending_energy"]),
        "fem_right_bending_energy": float(fem_energy["right_bending_energy"]),
        "fem_total_axial_energy": float(fem_energy["total_axial_energy"]),
        "fem_total_bending_energy": float(fem_energy["total_bending_energy"]),
        "fem_total_energy": float(fem_energy["total_energy"]),
        "fem_axial_energy_fraction": float(fem_energy["axial_energy_fraction_from_arm_energies"]),
        "fem_bending_energy_fraction": float(fem_bending_fraction),
        "fem_right_axial_share": float(fem_energy["right_axial_share"]),
        "fem_right_bending_share": float(fem_energy["right_bending_share"]),
    }
    row.update(
        {
            "axial_energy_fraction_difference": float(
                row["analytic_axial_energy_fraction"] - row["fem_axial_energy_fraction"]
            ),
            "right_axial_share_difference": float(
                row["analytic_right_axial_share"] - row["fem_right_axial_share"]
            ),
            "bending_energy_fraction_difference": float(
                row["analytic_bending_energy_fraction"] - row["fem_bending_energy_fraction"]
            ),
            "right_bending_share_difference": float(
                row["analytic_right_bending_share"] - row["fem_right_bending_share"]
            ),
            "analytic_right_axial_dominant": str(analytic_right_axial_dominant),
            "fem_right_axial_dominant": str(fem_right_axial_dominant),
            "right_axial_dominance_agreement": right_axial_dominance_agreement,
        }
    )
    return row


def write_energy_comparison_csv(path: Path, row: dict[str, float | int | str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def fem_element_stiffness_for_epsilon(le: float, epsilon: float) -> np.ndarray:
    saved_ea_nd = float(fem.EA_nd)
    saved_ei_nd = float(fem.EI_nd)
    try:
        fem.EA_nd = 1.0 / float(epsilon) ** 2
        fem.EI_nd = 1.0
        return np.asarray(fem.elem_K(float(le)), dtype=float)
    finally:
        fem.EA_nd = saved_ea_nd
        fem.EI_nd = saved_ei_nd


def fem_direct_element_energy_diagnostics(
    shape_case: dict[str, object],
    *,
    epsilon: float,
) -> dict[str, float]:
    """Compute FEM energies directly from current element stiffness blocks."""
    mu_value = float(shape_case["mu"])
    ell = DEFAULT_L_TOTAL / 2.0
    left_length = max(ell * (1.0 - mu_value), NEAR_ZERO_NORM)
    right_length = max(ell * (1.0 + mu_value), NEAR_ZERO_NORM)
    le_left = left_length / fem.N_ELEM
    le_right = right_length / fem.N_ELEM
    k_left = fem_element_stiffness_for_epsilon(le_left, epsilon)
    k_right = fem_element_stiffness_for_epsilon(le_right, epsilon)
    t_left = fem.rotation_matrix_6x6(0.0)
    t_right = fem.rotation_matrix_6x6(np.deg2rad(float(shape_case["beta"])))
    axial_idx = np.array([0, 3], dtype=int)
    bending_idx = np.array([1, 2, 4, 5], dtype=int)
    k_axial_left = k_left[np.ix_(axial_idx, axial_idx)]
    k_bending_left = k_left[np.ix_(bending_idx, bending_idx)]
    k_axial_right = k_right[np.ix_(axial_idx, axial_idx)]
    k_bending_right = k_right[np.ix_(bending_idx, bending_idx)]
    u_raw = np.asarray(shape_case["u_raw"], dtype=float)
    v_raw = np.asarray(shape_case["v_raw"], dtype=float)
    theta_raw = np.asarray(shape_case["theta_raw"], dtype=float)

    left_axial_energy = 0.0
    right_axial_energy = 0.0
    left_bending_energy = 0.0
    right_bending_energy = 0.0

    def q_global_for_nodes(node0: int, node1: int) -> np.ndarray:
        return np.array(
            [
                u_raw[node0],
                v_raw[node0],
                theta_raw[node0],
                u_raw[node1],
                v_raw[node1],
                theta_raw[node1],
            ],
            dtype=float,
        )

    for elem in range(fem.N_ELEM):
        q_left = t_left @ q_global_for_nodes(elem, elem + 1)
        q_left_axial = q_left[axial_idx]
        q_left_bending = q_left[bending_idx]
        left_axial_energy += 0.5 * float(q_left_axial @ k_axial_left @ q_left_axial)
        left_bending_energy += 0.5 * float(q_left_bending @ k_bending_left @ q_left_bending)

        node0 = fem.N_ELEM + elem
        node1 = fem.N_ELEM + elem + 1
        q_right = t_right @ q_global_for_nodes(node0, node1)
        q_right_axial = q_right[axial_idx]
        q_right_bending = q_right[bending_idx]
        right_axial_energy += 0.5 * float(q_right_axial @ k_axial_right @ q_right_axial)
        right_bending_energy += 0.5 * float(q_right_bending @ k_bending_right @ q_right_bending)

    left_axial_energy = max(float(left_axial_energy), 0.0)
    right_axial_energy = max(float(right_axial_energy), 0.0)
    left_bending_energy = max(float(left_bending_energy), 0.0)
    right_bending_energy = max(float(right_bending_energy), 0.0)
    total_axial = left_axial_energy + right_axial_energy
    total_bending = left_bending_energy + right_bending_energy
    total = total_axial + total_bending
    return {
        "left_axial_energy": left_axial_energy,
        "right_axial_energy": right_axial_energy,
        "left_bending_energy": left_bending_energy,
        "right_bending_energy": right_bending_energy,
        "total_axial_energy": float(total_axial),
        "total_bending_energy": float(total_bending),
        "total_energy": float(total),
        "axial_energy_fraction": safe_ratio(total_axial, total),
        "right_axial_share": safe_ratio(right_axial_energy, total_axial),
        "right_bending_share": safe_ratio(right_bending_energy, total_bending),
    }


def build_fem_direct_energy_check_row(
    *,
    branch_id: str,
    branch_number: int | None,
    beta_deg: float,
    mu_value: float,
    epsilon: float,
    l_total: float,
    radius: float,
    fem_lambda: float,
    shape_case: dict[str, object],
    output_path: Path,
) -> dict[str, float | int | str]:
    existing = arm_energy_diagnostics(shape_case, epsilon=epsilon)
    direct = fem_direct_element_energy_diagnostics(shape_case, epsilon=epsilon)
    energy_keys = (
        "left_axial_energy",
        "right_axial_energy",
        "left_bending_energy",
        "right_bending_energy",
        "total_axial_energy",
        "total_bending_energy",
        "total_energy",
    )
    comparison: dict[str, float] = {}
    max_abs_difference = 0.0
    max_relative_difference = 0.0
    for key in energy_keys:
        existing_value = float(existing[key])
        direct_value = float(direct[key])
        difference = direct_value - existing_value
        relative_difference = abs(difference) / abs(existing_value) if abs(existing_value) > NEAR_ZERO_NORM else np.nan
        comparison[f"fem_direct_minus_existing_{key}"] = float(difference)
        comparison[f"fem_direct_to_existing_{key}_ratio"] = safe_ratio(direct_value, existing_value)
        comparison[f"fem_direct_existing_{key}_relative_difference"] = float(relative_difference)
        max_abs_difference = max(max_abs_difference, abs(difference))
        if np.isfinite(relative_difference):
            max_relative_difference = max(max_relative_difference, float(relative_difference))

    fraction_difference = float(direct["axial_energy_fraction"] - existing["axial_energy_fraction_from_arm_energies"])
    right_share_difference = float(direct["right_axial_share"] - existing["right_axial_share"])
    max_fraction_difference = max(abs(fraction_difference), abs(right_share_difference))
    energy_agrees = bool(max_relative_difference <= 1e-9 or max_abs_difference <= 1e-12)
    fraction_agrees = bool(max_fraction_difference <= 1e-9)
    agrees = energy_agrees and fraction_agrees
    row: dict[str, float | int | str] = {
        "branch_id": branch_id,
        "branch_number": "" if branch_number is None else int(branch_number),
        "beta": float(beta_deg),
        "mu": float(mu_value),
        "epsilon": float(epsilon),
        "l_total": float(l_total),
        "radius": float(radius),
        "fem_lambda": float(fem_lambda),
        "fem_direct_energy_check_csv": str(output_path),
        "fem_direct_energy_source": "fem.elem_K local stiffness blocks with q_local = rotation_matrix_6x6(beta) @ q_global",
        "fem_existing_left_axial_energy": float(existing["left_axial_energy"]),
        "fem_existing_right_axial_energy": float(existing["right_axial_energy"]),
        "fem_existing_left_bending_energy": float(existing["left_bending_energy"]),
        "fem_existing_right_bending_energy": float(existing["right_bending_energy"]),
        "fem_existing_total_axial_energy": float(existing["total_axial_energy"]),
        "fem_existing_total_bending_energy": float(existing["total_bending_energy"]),
        "fem_existing_total_energy": float(existing["total_energy"]),
        "fem_existing_axial_energy_fraction": float(existing["axial_energy_fraction_from_arm_energies"]),
        "fem_existing_right_axial_share": float(existing["right_axial_share"]),
        "fem_existing_right_bending_share": float(existing["right_bending_share"]),
        "fem_direct_left_axial_energy": float(direct["left_axial_energy"]),
        "fem_direct_right_axial_energy": float(direct["right_axial_energy"]),
        "fem_direct_left_bending_energy": float(direct["left_bending_energy"]),
        "fem_direct_right_bending_energy": float(direct["right_bending_energy"]),
        "fem_direct_total_axial_energy": float(direct["total_axial_energy"]),
        "fem_direct_total_bending_energy": float(direct["total_bending_energy"]),
        "fem_direct_total_energy": float(direct["total_energy"]),
        "fem_direct_axial_energy_fraction": float(direct["axial_energy_fraction"]),
        "fem_direct_right_axial_share": float(direct["right_axial_share"]),
        "fem_direct_right_bending_share": float(direct["right_bending_share"]),
        "fem_direct_minus_existing_axial_energy_fraction": fraction_difference,
        "fem_direct_minus_existing_right_axial_share": right_share_difference,
        "fem_direct_to_existing_axial_energy_fraction_ratio": safe_ratio(
            float(direct["axial_energy_fraction"]),
            float(existing["axial_energy_fraction_from_arm_energies"]),
        ),
        "fem_direct_to_existing_right_axial_share_ratio": safe_ratio(
            float(direct["right_axial_share"]),
            float(existing["right_axial_share"]),
        ),
        "fem_direct_energy_max_abs_difference": float(max_abs_difference),
        "fem_direct_energy_max_relative_difference": float(max_relative_difference),
        "fem_direct_energy_agreement": "agree" if agrees else "disagree",
        "arm_energy_diagnostics_reliability": (
            "reliable for this element-stiffness check" if agrees else "not reliable for this element-stiffness check"
        ),
    }
    row.update(comparison)
    return row


def write_fem_direct_energy_check_csv(path: Path, row: dict[str, float | int | str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def write_samples_csv(
    path: Path,
    *,
    s_norm: np.ndarray,
    fem_components: dict[str, np.ndarray],
    analytic_components: dict[str, np.ndarray],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "s_norm",
        "fem_u_left",
        "analytic_u_left",
        "fem_w_left",
        "analytic_w_left",
        "fem_u_right",
        "analytic_u_right",
        "fem_w_right",
        "analytic_w_right",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for idx, s_value in enumerate(s_norm):
            writer.writerow(
                {
                    "s_norm": float(s_value),
                    "fem_u_left": float(fem_components["u_left"][idx]),
                    "analytic_u_left": float(analytic_components["u_left"][idx]),
                    "fem_w_left": float(fem_components["w_left"][idx]),
                    "analytic_w_left": float(analytic_components["w_left"][idx]),
                    "fem_u_right": float(fem_components["u_right"][idx]),
                    "analytic_u_right": float(analytic_components["u_right"][idx]),
                    "fem_w_right": float(fem_components["w_right"][idx]),
                    "analytic_w_right": float(analytic_components["w_right"][idx]),
                }
            )


def write_diagnostics_csv(path: Path, row: dict[str, float | int | str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def plot_components_overlay(
    path: Path,
    *,
    s_norm: np.ndarray,
    fem_components: dict[str, np.ndarray],
    analytic_components: dict[str, np.ndarray],
    title: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), sharex=True)
    panels = [
        ("u_left", "Left axial"),
        ("w_left", "Left transverse"),
        ("u_right", "Right axial"),
        ("w_right", "Right transverse"),
    ]
    for ax, (component, label) in zip(axes.ravel(), panels):
        ax.plot(s_norm, fem_components[component], color="#1f77b4", linewidth=2.0, label="FEM")
        ax.plot(
            s_norm,
            analytic_components[component],
            color="#d62728",
            linestyle="--",
            linewidth=2.0,
            label="analytic",
        )
        ax.set_title(label)
        ax.set_xlabel("s / L")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=9, loc="best")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(pad=0.9)
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def build_diagnostics_row(
    *,
    branch_id: str,
    branch_number: int | None,
    beta_deg: float,
    mu_value: float,
    epsilon: float,
    l_total: float,
    radius: float,
    fem_lambda: float,
    root_info: dict[str, float | int | str | list[float]],
    current_sorted_index: int,
    smallest_singular_value: float,
    singular_value_ratio: float,
    normalization_info: dict[str, float],
    errors: dict[str, float],
    coordinate_convention_used: str,
    right_fem_arrays_reversed: bool,
    output_paths: dict[str, Path],
    selected_variant_row: dict[str, float | int | str],
    current_variant_row: dict[str, float | int | str],
    best_variant_row: dict[str, float | int | str],
    used_best_orientation: bool,
) -> dict[str, float | int | str]:
    analytic_lambda = float(root_info["analytic_lambda"])
    relative_lambda_difference = (
        abs(analytic_lambda - fem_lambda) / abs(fem_lambda) if abs(fem_lambda) > NEAR_ZERO_NORM else np.nan
    )
    row: dict[str, float | int | str] = {
        "branch_id": branch_id,
        "branch_number": "" if branch_number is None else int(branch_number),
        "beta": float(beta_deg),
        "mu": float(mu_value),
        "epsilon": float(epsilon),
        "l_total": float(l_total),
        "radius": float(radius),
        "fem_lambda": float(fem_lambda),
        "analytic_lambda": analytic_lambda,
        "relative_lambda_difference": float(relative_lambda_difference),
        "current_sorted_index": int(current_sorted_index),
        "analytic_root_index": int(root_info["analytic_root_index"]),
        "analytic_root_source": str(root_info["analytic_root_source"]),
        "analytic_root_abs_det": float(root_info["analytic_root_abs_det"]),
        "global_roots_found": int(root_info["global_roots_found"]),
        "local_roots_found": int(root_info["local_roots_found"]),
        "smallest_singular_value": float(smallest_singular_value),
        "singular_value_ratio": float(singular_value_ratio),
        "shape_dot_before_sign_alignment": float(normalization_info["shape_dot_before_sign_alignment"]),
        "analytic_sign_applied": float(normalization_info["analytic_sign_applied"]),
        "fem_normalization_max_abs": float(normalization_info["fem_normalization_max_abs"]),
        "analytic_normalization_max_abs": float(normalization_info["analytic_normalization_max_abs"]),
        "coordinate_convention_used": coordinate_convention_used,
        "right_fem_arrays_reversed": str(bool(right_fem_arrays_reversed)),
        "used_best_orientation": str(bool(used_best_orientation)),
        "selected_orientation_variant_id": int(selected_variant_row["variant_id"]),
        "selected_orientation_variant": str(selected_variant_row["variant_label"]),
        "current_orientation_variant_id": int(current_variant_row["variant_id"]),
        "current_orientation_variant": str(current_variant_row["variant_label"]),
        "current_orientation_relative_L2_error_all_components": float(
            current_variant_row["relative_L2_error_all_components"]
        ),
        "best_orientation_variant_id": int(best_variant_row["variant_id"]),
        "best_orientation_variant": str(best_variant_row["variant_label"]),
        "best_orientation_relative_L2_error_all_components": float(
            best_variant_row["relative_L2_error_all_components"]
        ),
        "best_orientation_analytic_max_abs_u_right_normalized": float(
            best_variant_row["analytic_max_abs_u_right_normalized"]
        ),
        "components_overlay_png": str(output_paths["png"]),
        "components_samples_csv": str(output_paths["samples_csv"]),
        "diagnostics_csv": str(output_paths["diagnostics_csv"]),
        "orientation_scan_csv": str(output_paths["orientation_scan_csv"]),
    }
    row.update(errors)
    return row


def print_summary(row: dict[str, float | int | str]) -> None:
    print("Analytic vs FEM tracked descendant shape comparison")
    print(f"branch_id: {row['branch_id']}")
    print(f"branch_number: {row['branch_number'] if row['branch_number'] != '' else 'not available'}")
    print(f"beta: {float(row['beta']):g} deg")
    print(f"mu: {float(row['mu']):g}")
    print(f"epsilon: {float(row['epsilon']):g}")
    print(f"fem_lambda: {float(row['fem_lambda']):.10g}")
    print(f"analytic_lambda: {float(row['analytic_lambda']):.10g}")
    print(f"relative_lambda_difference: {float(row['relative_lambda_difference']):.6e}")
    print(f"current_sorted_index: {int(row['current_sorted_index'])}")
    print(f"analytic_root_index: {int(row['analytic_root_index'])}")
    print(f"analytic_root_source: {row['analytic_root_source']}")
    print(f"analytic_root_abs_det: {float(row['analytic_root_abs_det']):.6e}")
    print(f"smallest_singular_value: {float(row['smallest_singular_value']):.6e}")
    print(f"singular_value_ratio: {float(row['singular_value_ratio']):.6e}")
    print(f"coordinate_convention_used: {row['coordinate_convention_used']}")
    print(f"right_fem_arrays_reversed: {row['right_fem_arrays_reversed']}")
    print(f"used_best_orientation: {row['used_best_orientation']}")
    print(f"selected_orientation_variant: {row['selected_orientation_variant']}")
    print(f"shape_dot_after_sign_alignment: {float(row['shape_dot_after_sign_alignment']):.6g}")
    print(f"relative_L2_error_all_components: {float(row['relative_L2_error_all_components']):.6e}")
    print(f"max_abs_error_all_components: {float(row['max_abs_error_all_components']):.6e}")
    for component in COMPONENT_KEYS:
        print(f"rel_l2_{component}: {float(row[f'rel_l2_{component}']):.6e}")
        print(f"max_abs_error_{component}: {float(row[f'max_abs_error_{component}']):.6e}")
    print(f"fem_max_abs_u_right_normalized: {float(row['fem_max_abs_u_right_normalized']):.6e}")
    print(f"analytic_max_abs_u_right_normalized: {float(row['analytic_max_abs_u_right_normalized']):.6e}")
    print(f"fem_right_axial_to_right_transverse: {float(row['fem_right_axial_to_right_transverse']):.6e}")
    print(f"analytic_right_axial_to_right_transverse: {float(row['analytic_right_axial_to_right_transverse']):.6e}")
    print(f"current_orientation_relative_L2_error_all_components: {float(row['current_orientation_relative_L2_error_all_components']):.6e}")
    print(f"best_orientation_variant: {row['best_orientation_variant']}")
    print(f"best_orientation_relative_L2_error_all_components: {float(row['best_orientation_relative_L2_error_all_components']):.6e}")
    print(
        "best_orientation_analytic_max_abs_u_right_normalized: "
        f"{float(row['best_orientation_analytic_max_abs_u_right_normalized']):.6e}"
    )
    print(f"saved PNG: {row['components_overlay_png']}")
    print(f"saved samples CSV: {row['components_samples_csv']}")
    print(f"saved diagnostics CSV: {row['diagnostics_csv']}")
    print(f"saved orientation scan CSV: {row['orientation_scan_csv']}")


def print_orientation_scan_top(rows: list[dict[str, float | int | str]], *, limit: int = 5) -> None:
    print(f"top {min(limit, len(rows))} orientation/sign variants by relative_L2_error_all_components:")
    for row in rows[:limit]:
        print(
            f"  #{int(row['variant_id']):02d}: rel_l2_all={float(row['relative_L2_error_all_components']):.6e}, "
            f"dot={float(row['shape_dot_after_sign_alignment']):.6g}, "
            f"uR_analytic={float(row['analytic_max_abs_u_right_normalized']):.6e}, "
            f"{row['variant_label']}"
        )


def print_energy_summary(row: dict[str, float | int | str]) -> None:
    print("Energy comparison diagnostics")
    print(f"analytic_axial_energy_fraction: {float(row['analytic_axial_energy_fraction']):.6e}")
    print(f"fem_axial_energy_fraction: {float(row['fem_axial_energy_fraction']):.6e}")
    print(f"analytic_right_axial_share: {float(row['analytic_right_axial_share']):.6e}")
    print(f"fem_right_axial_share: {float(row['fem_right_axial_share']):.6e}")
    print(f"Energy diagnostics {row['right_axial_dominance_agreement']} on right axial dominance.")
    print(f"saved energy comparison CSV: {row['energy_comparison_csv']}")


def print_fem_direct_energy_check_summary(row: dict[str, float | int | str]) -> None:
    print("Direct FEM element-stiffness energy check")
    print(f"existing_fem_axial_energy_fraction: {float(row['fem_existing_axial_energy_fraction']):.6e}")
    print(f"direct_fem_axial_energy_fraction: {float(row['fem_direct_axial_energy_fraction']):.6e}")
    print(f"existing_fem_right_axial_share: {float(row['fem_existing_right_axial_share']):.6e}")
    print(f"direct_fem_right_axial_share: {float(row['fem_direct_right_axial_share']):.6e}")
    print(f"fem_direct_energy_agreement: {row['fem_direct_energy_agreement']}")
    print(f"arm_energy_diagnostics_reliability: {row['arm_energy_diagnostics_reliability']}")
    print(f"saved FEM direct energy check CSV: {row['fem_direct_energy_check_csv']}")


def main(argv: Sequence[str] | None = None) -> dict[str, float | int | str]:
    args = parse_args(argv)
    if args.branch_id is None:
        branch_number = args.branch_number
        branch_id = branch_id_from_number(branch_number)
    else:
        branch_id = str(args.branch_id)
        branch_number = branch_number_from_id(branch_id)

    radius = radius_from_epsilon(args.epsilon, l_total=args.l_total)
    beta_rad = float(np.deg2rad(args.beta))
    shape_case = collect_single_branch_shape(
        branch_id=branch_id,
        beta_deg=args.beta,
        mu_value=args.mu,
        radius=radius,
    )
    fem_lambda = float(shape_case["lambda"])

    root_info = find_analytic_root_near_fem(
        fem_lambda,
        beta_rad=beta_rad,
        mu_value=args.mu,
        epsilon=args.epsilon,
        num_roots=args.num_roots,
        root_window=args.root_window,
    )
    analytic_lambda = float(root_info["analytic_lambda"])
    coeff, smallest_singular_value, singular_value_ratio = analytic_null_vector(
        analytic_lambda,
        beta_rad=beta_rad,
        mu_value=args.mu,
        epsilon=args.epsilon,
    )

    s_norm = np.linspace(0.0, 1.0, fem.N_ELEM + 1)
    current_result = evaluate_orientation_variant(
        variant_id=0,
        shape_case=shape_case,
        analytic_lambda=analytic_lambda,
        coeff=coeff,
        s_norm=s_norm,
        mu_value=args.mu,
        epsilon=args.epsilon,
        right_coordinate=args.right_coordinate,
        sign_u_left=1.0,
        sign_w_left=1.0,
        sign_u_right=1.0,
        sign_w_right=1.0,
    )
    orientation_scan_rows, best_result = scan_orientation_variants(
        shape_case=shape_case,
        analytic_lambda=analytic_lambda,
        coeff=coeff,
        s_norm=s_norm,
        mu_value=args.mu,
        epsilon=args.epsilon,
    )
    selected_result = best_result if args.use_best_orientation else current_result
    selected_variant_row = selected_result["row"]
    current_variant_row = current_result["row"]
    best_variant_row = best_result["row"]
    fem_normed = selected_result["fem_components"]
    analytic_normed = selected_result["analytic_components"]
    normalization_info = selected_result["normalization_info"]
    errors = selected_result["errors"]
    right_fem_arrays_reversed = bool(selected_result["right_fem_arrays_reversed"])

    output_prefix = resolve_output_prefix(
        args.output_prefix,
        branch_id=branch_id,
        beta_deg=args.beta,
        mu_value=args.mu,
        epsilon=args.epsilon,
    )
    output_paths = {
        "png": suffixed_path(output_prefix, "_components_overlay.png"),
        "samples_csv": suffixed_path(output_prefix, "_components_samples.csv"),
        "diagnostics_csv": suffixed_path(output_prefix, "_diagnostics.csv"),
        "orientation_scan_csv": suffixed_path(output_prefix, "_orientation_scan.csv"),
        "energy_comparison_csv": suffixed_path(output_prefix, "_energy_comparison.csv"),
        "fem_direct_energy_check_csv": suffixed_path(output_prefix, "_fem_direct_energy_check.csv"),
    }

    title = (
        f"{branch_id}, beta={args.beta:g} deg, mu={args.mu:g}, eps={args.epsilon:g}, "
        f"FEM Lambda={fem_lambda:.4g}, analytic Lambda={analytic_lambda:.4g}, "
        f"variant #{int(selected_variant_row['variant_id'])}"
    )
    plot_components_overlay(
        output_paths["png"],
        s_norm=s_norm,
        fem_components=fem_normed,
        analytic_components=analytic_normed,
        title=title,
    )
    write_samples_csv(
        output_paths["samples_csv"],
        s_norm=s_norm,
        fem_components=fem_normed,
        analytic_components=analytic_normed,
    )
    write_orientation_scan_csv(output_paths["orientation_scan_csv"], orientation_scan_rows)
    row = build_diagnostics_row(
        branch_id=branch_id,
        branch_number=branch_number,
        beta_deg=args.beta,
        mu_value=args.mu,
        epsilon=args.epsilon,
        l_total=args.l_total,
        radius=radius,
        fem_lambda=fem_lambda,
        root_info=root_info,
        current_sorted_index=int(shape_case["current_sorted_index"]),
        smallest_singular_value=smallest_singular_value,
        singular_value_ratio=singular_value_ratio,
        normalization_info=normalization_info,
        errors=errors,
        coordinate_convention_used=(
            f"left=external-to-joint; right={selected_variant_row['right_coordinate']}; "
            "analytic component signs are diagnostic orientation factors"
        ),
        right_fem_arrays_reversed=right_fem_arrays_reversed,
        output_paths=output_paths,
        selected_variant_row=selected_variant_row,
        current_variant_row=current_variant_row,
        best_variant_row=best_variant_row,
        used_best_orientation=bool(args.use_best_orientation),
    )
    energy_row: dict[str, float | int | str] | None = None
    if args.compare_energies:
        energy_row = build_energy_comparison_row(
            branch_id=branch_id,
            branch_number=branch_number,
            beta_deg=args.beta,
            mu_value=args.mu,
            epsilon=args.epsilon,
            l_total=args.l_total,
            radius=radius,
            fem_lambda=fem_lambda,
            analytic_lambda=analytic_lambda,
            shape_case=shape_case,
            analytic_components=analytic_normed,
            s_norm=s_norm,
            normalization_info=normalization_info,
            selected_variant_row=selected_variant_row,
            output_path=output_paths["energy_comparison_csv"],
        )
        write_energy_comparison_csv(output_paths["energy_comparison_csv"], energy_row)
        row.update(energy_row)
    fem_direct_energy_row: dict[str, float | int | str] | None = None
    if args.check_fem_direct_energy:
        fem_direct_energy_row = build_fem_direct_energy_check_row(
            branch_id=branch_id,
            branch_number=branch_number,
            beta_deg=args.beta,
            mu_value=args.mu,
            epsilon=args.epsilon,
            l_total=args.l_total,
            radius=radius,
            fem_lambda=fem_lambda,
            shape_case=shape_case,
            output_path=output_paths["fem_direct_energy_check_csv"],
        )
        write_fem_direct_energy_check_csv(output_paths["fem_direct_energy_check_csv"], fem_direct_energy_row)
        row.update(fem_direct_energy_row)
    write_diagnostics_csv(output_paths["diagnostics_csv"], row)
    print_summary(row)
    if energy_row is not None:
        print_energy_summary(energy_row)
    if fem_direct_energy_row is not None:
        print_fem_direct_energy_check_summary(fem_direct_energy_row)
    print_orientation_scan_top(orientation_scan_rows, limit=5)
    current_rel_l2 = float(current_variant_row["relative_L2_error_all_components"])
    best_rel_l2 = float(best_variant_row["relative_L2_error_all_components"])
    best_is_much_better = (
        best_rel_l2 < ORIENTATION_WARNING_RELATIVE_FACTOR * current_rel_l2
        or current_rel_l2 - best_rel_l2 > ORIENTATION_WARNING_ABSOLUTE_GAIN
    )
    if best_is_much_better and not args.use_best_orientation:
        print(
            "warning: A different orientation/sign convention gives much better agreement; "
            "inspect orientation_scan.csv."
        )
    if args.check_endpoint_consistency:
        run_endpoint_consistency_diagnostics(
            output_prefix=output_prefix,
            analytic_lambda=analytic_lambda,
            beta_rad=beta_rad,
            mu_value=args.mu,
            epsilon=args.epsilon,
            coeff=coeff,
            shape_case=shape_case,
        )
    return row


if __name__ == "__main__":
    main()
