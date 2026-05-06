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
    collect_single_branch_shape,
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
        "components_overlay_png": str(output_paths["png"]),
        "components_samples_csv": str(output_paths["samples_csv"]),
        "diagnostics_csv": str(output_paths["diagnostics_csv"]),
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
    print(f"saved PNG: {row['components_overlay_png']}")
    print(f"saved samples CSV: {row['components_samples_csv']}")
    print(f"saved diagnostics CSV: {row['diagnostics_csv']}")


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
    fem_components, right_fem_arrays_reversed = fem_components_for_comparison(
        shape_case,
        right_coordinate=args.right_coordinate,
    )
    analytic_components = reconstruct_analytic_components(
        analytic_lambda,
        mu_value=args.mu,
        epsilon=args.epsilon,
        coeff=coeff,
        s_norm=s_norm,
        right_coordinate=args.right_coordinate,
    )
    fem_normed, analytic_normed, normalization_info = normalize_components(fem_components, analytic_components)
    errors = comparison_errors(fem_normed, analytic_normed)

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
    }

    title = (
        f"{branch_id}, beta={args.beta:g} deg, mu={args.mu:g}, eps={args.epsilon:g}, "
        f"FEM Lambda={fem_lambda:.4g}, analytic Lambda={analytic_lambda:.4g}"
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
        coordinate_convention_used=f"left=external-to-joint; right={args.right_coordinate}",
        right_fem_arrays_reversed=right_fem_arrays_reversed,
        output_paths=output_paths,
    )
    write_diagnostics_csv(output_paths["diagnostics_csv"], row)
    print_summary(row)
    return row


if __name__ == "__main__":
    main()
