"""Check the beta=0 single-rod limit of the current determinant approach.

At beta=0 the two clamped rods are collinear and should reproduce one
clamped-clamped Euler-Bernoulli beam of total nondimensional length 2, with an
artificial internal joint at x = 1 - mu.  This diagnostic intentionally uses the
current determinant matrix and its SVD null vector without changing formulas.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.optimize import brentq, minimize_scalar

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
from my_project.fem import python_fem as fem  # noqa: E402
from scripts.analysis.check_analytic_shape_in_fem_residual import (  # noqa: E402
    assemble_fem_matrices,
    build_params_for_case,
    residual_metrics,
)
from scripts.compare_beta0_analytic_vs_fem import fem_parameter_override  # noqa: E402

DEFAULT_EPSILON = 0.0025
DEFAULT_BETA = 0.0
DEFAULT_MU_VALUES = (0.0, 0.2, 0.6)
DEFAULT_NUM_BENDING_MODES = 5
DEFAULT_N_FEM_MODES = 20
DEFAULT_OUTPUT_PREFIX = Path("results/beta0_single_rod_limit")
DEFAULT_N_ELEM = fem.N_ELEM

KNOWN_CLAMPED_CLAMPED_ALPHA = (
    4.730040744862704,
    7.853204624095838,
    10.99560783800167,
    14.137165491257,
    17.278759657399,
)


@dataclass(frozen=True)
class RootRecord:
    mu: float
    mode_number: int
    expected_alpha: float
    expected_lambda: float
    determinant_lambda: float
    determinant_rel_error: float
    fem_lambda: float
    fem_rel_error: float
    determinant_sorted_index: int
    determinant_abs_det: float


@dataclass(frozen=True)
class ShapeRecord:
    mu: float
    mode_number: int
    shape_mac_det_vs_exact: float
    shape_rel_l2_det_vs_exact: float
    shape_max_abs_det_vs_exact: float
    fem_residual_relative: float
    rayleigh_lambda: float
    rayleigh_lambda_rel_error: float
    joint_w_mismatch: float
    joint_slope_mismatch: float
    external_clamp_max_abs: float


def parse_mu_values(raw: str) -> tuple[float, ...]:
    values = tuple(float(part.strip()) for part in raw.split(",") if part.strip())
    if not values:
        raise argparse.ArgumentTypeError("expected at least one comma-separated mu value")
    return values


def output_paths(prefix: Path) -> dict[str, Path]:
    if not prefix.is_absolute():
        prefix = REPO_ROOT / prefix
    return {
        "roots": prefix.with_name(prefix.name + "_roots.csv"),
        "shapes": prefix.with_name(prefix.name + "_shapes.csv"),
        "samples": prefix.with_name(prefix.name + "_shape_samples.csv"),
        "summary": prefix.with_name(prefix.name + "_summary.csv"),
    }


def determinant_value(lam: float, beta: float, mu: float, eps: float) -> float:
    value = det_clamped_coupled(lam, beta, mu, eps)
    return float(np.real_if_close(value))


def find_determinant_root_near(expected_lambda: float, beta: float, mu: float, eps: float) -> tuple[float, float]:
    width = 0.35
    lower = max(1.0e-8, expected_lambda - width)
    upper = expected_lambda + width
    grid = np.linspace(lower, upper, 801)
    values = np.array([determinant_value(float(x), beta, mu, eps) for x in grid])

    sign_change_candidates: list[tuple[float, float, float]] = []
    for a, b, fa, fb in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if not (np.isfinite(fa) and np.isfinite(fb)):
            continue
        if fa == 0.0:
            return float(a), 0.0
        if fa * fb < 0.0:
            midpoint = 0.5 * (a + b)
            sign_change_candidates.append((abs(midpoint - expected_lambda), float(a), float(b)))

    if sign_change_candidates:
        _, a, b = min(sign_change_candidates, key=lambda item: item[0])
        root = brentq(lambda x: determinant_value(x, beta, mu, eps), a, b, xtol=1.0e-13, rtol=1.0e-13)
        return float(root), abs(determinant_value(root, beta, mu, eps))

    min_index = int(np.nanargmin(np.abs(values)))
    local_lower = float(grid[max(0, min_index - 2)])
    local_upper = float(grid[min(len(grid) - 1, min_index + 2)])
    if local_lower == local_upper:
        local_lower, local_upper = lower, upper
    result = minimize_scalar(
        lambda x: abs(determinant_value(float(x), beta, mu, eps)),
        bounds=(local_lower, local_upper),
        method="bounded",
        options={"xatol": 1.0e-13},
    )
    lam = float(result.x)
    return lam, abs(determinant_value(lam, beta, mu, eps))


def determinant_null_vector(lam: float, beta: float, mu: float, eps: float) -> np.ndarray:
    matrix = assemble_clamped_coupled_matrix(lam, beta, mu, eps)
    _, _, vh = np.linalg.svd(matrix)
    coeff = np.asarray(vh[-1, :], dtype=float)
    pivot = int(np.argmax(np.abs(coeff)))
    if coeff[pivot] < 0.0:
        coeff = -coeff
    norm = np.linalg.norm(coeff)
    if norm == 0.0:
        raise RuntimeError("determinant SVD returned a zero null vector")
    return coeff / norm


def reduced_bending_w(a: float, b: float, z: np.ndarray) -> np.ndarray:
    return a * (np.cos(z) - np.cosh(z)) + b * (np.sin(z) - np.sinh(z))


def reduced_bending_slope_z(a: float, b: float, z: np.ndarray) -> np.ndarray:
    return a * (-np.sin(z) - np.sinh(z)) + b * (np.cos(z) - np.cosh(z))


def determinant_fields(
    lam: float,
    mu: float,
    eps: float,
    coeff: np.ndarray,
    n_elem: int = DEFAULT_N_ELEM,
) -> dict[str, np.ndarray | float]:
    a1, b1, a2, b2, p1, p2 = coeff
    l1 = 1.0 - mu
    l2 = 1.0 + mu

    s_left = np.linspace(0.0, l1, n_elem + 1)
    x_left = s_left
    z_left = lam * s_left
    w_left = reduced_bending_w(a1, b1, z_left)
    theta_left = lam * reduced_bending_slope_z(a1, b1, z_left)
    ux_left = p1 * np.sin(eps * lam**2 * s_left)

    # Right determinant coefficients use an external-to-joint local argument
    # z = -Lambda * xi.  For global x we reverse it to joint-to-right-clamp:
    # x = L1 + s_right, xi = L2 - s_right, z = -Lambda * (L2 - s_right).
    s_right = np.linspace(0.0, l2, n_elem + 1)
    x_right = l1 + s_right
    z_right = -lam * (l2 - s_right)
    w_right = reduced_bending_w(a2, b2, z_right)
    theta_right = lam * reduced_bending_slope_z(a2, b2, z_right)
    ux_right = p2 * np.sin(-eps * lam**2 * (l2 - s_right))

    x_global = np.concatenate([x_left, x_right[1:]])
    w_global = np.concatenate([w_left, w_right[1:]])
    theta_global = np.concatenate([theta_left, theta_right[1:]])
    ux_global = np.concatenate([ux_left, ux_right[1:]])

    return {
        "x": x_global,
        "ux": ux_global,
        "w": w_global,
        "theta": theta_global,
        "left_joint_w": float(w_left[-1]),
        "right_joint_w": float(w_right[0]),
        "left_joint_theta": float(theta_left[-1]),
        "right_joint_theta": float(theta_right[0]),
        "left_clamp": np.array([ux_left[0], w_left[0], theta_left[0]], dtype=float),
        "right_clamp": np.array([ux_right[-1], w_right[-1], theta_right[-1]], dtype=float),
    }


def exact_clamped_clamped_shape(alpha: float, x: np.ndarray, total_length: float = 2.0) -> tuple[np.ndarray, np.ndarray]:
    k = alpha / total_length
    sigma = (math.cosh(alpha) - math.cos(alpha)) / (math.sinh(alpha) - math.sin(alpha))
    kx = k * x
    w = np.cosh(kx) - np.cos(kx) - sigma * (np.sinh(kx) - np.sin(kx))
    theta = k * (np.sinh(kx) + np.sin(kx) - sigma * (np.cosh(kx) - np.cos(kx)))
    return w, theta


def sign_aligned_metrics(det: np.ndarray, exact: np.ndarray) -> tuple[float, float, float, np.ndarray, np.ndarray]:
    det = np.asarray(det, dtype=float)
    exact = np.asarray(exact, dtype=float)
    det_norm = np.linalg.norm(det)
    exact_norm = np.linalg.norm(exact)
    if det_norm == 0.0 or exact_norm == 0.0:
        return float("nan"), float("nan"), float("nan"), det, exact
    det_unit = det / det_norm
    exact_unit = exact / exact_norm
    if float(np.dot(det_unit, exact_unit)) < 0.0:
        det_unit = -det_unit
    dot = float(np.dot(det_unit, exact_unit))
    mac = dot * dot
    diff = det_unit - exact_unit
    rel_l2 = float(np.linalg.norm(diff) / np.linalg.norm(exact_unit))
    max_abs = float(np.max(np.abs(diff)))
    return mac, rel_l2, max_abs, det_unit, exact_unit


def build_full_fem_vector(fields: dict[str, np.ndarray | float]) -> np.ndarray:
    ux = np.asarray(fields["ux"], dtype=float)
    w = np.asarray(fields["w"], dtype=float)
    theta = np.asarray(fields["theta"], dtype=float)
    q_full = np.zeros(3 * len(w), dtype=float)
    q_full[0::3] = ux
    q_full[1::3] = w
    q_full[2::3] = theta
    return q_full


def fem_lambdas_for_case(mu: float, beta_deg: float, epsilon: float, n_modes: int) -> np.ndarray:
    params = build_params_for_case(epsilon, l_total=2.0)
    with fem_parameter_override(params):
        omega, _ = fem.fem_solve(mu=mu, beta_deg=beta_deg, n_modes=n_modes)
    return np.sqrt(np.asarray(omega, dtype=float))


def nearest_fem_lambda(expected_lambda: float, fem_lambdas: np.ndarray) -> float:
    index = int(np.argmin(np.abs(fem_lambdas - expected_lambda)))
    return float(fem_lambdas[index])


def fem_residual_and_rayleigh(
    mu: float,
    epsilon: float,
    q_full: np.ndarray,
    lam: float,
) -> tuple[float, float, float]:
    params = build_params_for_case(epsilon, l_total=2.0)
    k_free, m_free, _, _, free = assemble_fem_matrices(params, mu=mu, beta_deg=0.0)
    metrics = residual_metrics(k_free, m_free, q_full, free, omega_sq=lam**4)
    q_free = q_full[free]
    numerator = float(q_free @ (k_free @ q_free))
    denominator = float(q_free @ (m_free @ q_free))
    if denominator <= 0.0 or numerator < 0.0:
        rayleigh_lambda = float("nan")
        rayleigh_rel = float("nan")
    else:
        rayleigh_lambda = float((numerator / denominator) ** 0.25)
        rayleigh_rel = abs(rayleigh_lambda - lam) / abs(lam)
    return float(metrics["relative_residual"]), rayleigh_lambda, rayleigh_rel


def external_clamp_max_abs(fields: dict[str, np.ndarray | float]) -> float:
    left = np.asarray(fields["left_clamp"], dtype=float)
    right = np.asarray(fields["right_clamp"], dtype=float)
    return float(np.max(np.abs(np.concatenate([left, right]))))


def write_csv(path: Path, rows: Iterable[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def run_diagnostic(
    epsilon: float,
    mu_values: tuple[float, ...],
    num_bending_modes: int,
    n_fem_modes: int,
    output_prefix: Path,
) -> None:
    if num_bending_modes > len(KNOWN_CLAMPED_CLAMPED_ALPHA):
        raise ValueError(f"only {len(KNOWN_CLAMPED_CLAMPED_ALPHA)} exact roots are built in")

    beta = DEFAULT_BETA
    paths = output_paths(output_prefix)
    root_records: list[RootRecord] = []
    shape_records: list[ShapeRecord] = []
    sample_rows: list[dict[str, object]] = []

    for mu in mu_values:
        fem_lambdas = fem_lambdas_for_case(mu, beta_deg=0.0, epsilon=epsilon, n_modes=n_fem_modes)
        for mode_number, alpha in enumerate(KNOWN_CLAMPED_CLAMPED_ALPHA[:num_bending_modes], start=1):
            expected_lambda = alpha / 2.0
            det_lambda, abs_det = find_determinant_root_near(expected_lambda, beta, mu, epsilon)
            fem_lambda = nearest_fem_lambda(expected_lambda, fem_lambdas)
            coeff = determinant_null_vector(det_lambda, beta, mu, epsilon)
            fields = determinant_fields(det_lambda, mu, epsilon, coeff)
            q_full = build_full_fem_vector(fields)
            fem_residual, rayleigh_lambda, rayleigh_rel = fem_residual_and_rayleigh(
                mu, epsilon, q_full, det_lambda
            )

            x = np.asarray(fields["x"], dtype=float)
            w_det = np.asarray(fields["w"], dtype=float)
            w_exact, _ = exact_clamped_clamped_shape(alpha, x)
            mac, rel_l2, max_abs, w_det_unit, w_exact_unit = sign_aligned_metrics(w_det, w_exact)

            joint_w = abs(float(fields["left_joint_w"]) - float(fields["right_joint_w"]))
            joint_slope = abs(float(fields["left_joint_theta"]) - float(fields["right_joint_theta"]))
            clamp_abs = external_clamp_max_abs(fields)

            root_records.append(
                RootRecord(
                    mu=mu,
                    mode_number=mode_number,
                    expected_alpha=alpha,
                    expected_lambda=expected_lambda,
                    determinant_lambda=det_lambda,
                    determinant_rel_error=abs(det_lambda - expected_lambda) / abs(expected_lambda),
                    fem_lambda=fem_lambda,
                    fem_rel_error=abs(fem_lambda - expected_lambda) / abs(expected_lambda),
                    determinant_sorted_index=mode_number,
                    determinant_abs_det=abs_det,
                )
            )
            shape_records.append(
                ShapeRecord(
                    mu=mu,
                    mode_number=mode_number,
                    shape_mac_det_vs_exact=mac,
                    shape_rel_l2_det_vs_exact=rel_l2,
                    shape_max_abs_det_vs_exact=max_abs,
                    fem_residual_relative=fem_residual,
                    rayleigh_lambda=rayleigh_lambda,
                    rayleigh_lambda_rel_error=rayleigh_rel,
                    joint_w_mismatch=joint_w,
                    joint_slope_mismatch=joint_slope,
                    external_clamp_max_abs=clamp_abs,
                )
            )
            for x_value, det_value, exact_value in zip(x, w_det_unit, w_exact_unit):
                sample_rows.append(
                    {
                        "mu": mu,
                        "mode_number": mode_number,
                        "x": float(x_value),
                        "w_det": float(det_value),
                        "w_exact": float(exact_value),
                        "diff": float(det_value - exact_value),
                    }
                )

    write_csv(
        paths["roots"],
        (record.__dict__ for record in root_records),
        [
            "mu",
            "mode_number",
            "expected_alpha",
            "expected_lambda",
            "determinant_lambda",
            "determinant_rel_error",
            "fem_lambda",
            "fem_rel_error",
            "determinant_sorted_index",
            "determinant_abs_det",
        ],
    )
    write_csv(
        paths["shapes"],
        (record.__dict__ for record in shape_records),
        [
            "mu",
            "mode_number",
            "shape_mac_det_vs_exact",
            "shape_rel_l2_det_vs_exact",
            "shape_max_abs_det_vs_exact",
            "fem_residual_relative",
            "rayleigh_lambda",
            "rayleigh_lambda_rel_error",
            "joint_w_mismatch",
            "joint_slope_mismatch",
            "external_clamp_max_abs",
        ],
    )
    write_csv(paths["samples"], sample_rows, ["mu", "mode_number", "x", "w_det", "w_exact", "diff"])

    summary_rows: list[dict[str, object]] = []
    for mode_number in range(1, num_bending_modes + 1):
        det_values = np.array(
            [record.determinant_lambda for record in root_records if record.mode_number == mode_number],
            dtype=float,
        )
        fem_values = np.array(
            [record.fem_lambda for record in root_records if record.mode_number == mode_number],
            dtype=float,
        )
        summary_rows.append(
            {
                "mode_number": mode_number,
                "determinant_lambda_min": float(np.min(det_values)),
                "determinant_lambda_max": float(np.max(det_values)),
                "determinant_lambda_spread": float(np.max(det_values) - np.min(det_values)),
                "fem_lambda_min": float(np.min(fem_values)),
                "fem_lambda_max": float(np.max(fem_values)),
                "fem_lambda_spread": float(np.max(fem_values) - np.min(fem_values)),
                "notes": "beta=0 should be mu-invariant",
            }
        )
    write_csv(
        paths["summary"],
        summary_rows,
        [
            "mode_number",
            "determinant_lambda_min",
            "determinant_lambda_max",
            "determinant_lambda_spread",
            "fem_lambda_min",
            "fem_lambda_max",
            "fem_lambda_spread",
            "notes",
        ],
    )

    print("beta=0 single-rod determinant diagnostic")
    print(f"epsilon={epsilon:g}, mu_values={mu_values}, modes={num_bending_modes}")
    print("Root comparison:")
    for record in root_records:
        print(
            "  "
            f"mu={record.mu:g}, mode={record.mode_number}: "
            f"det rel err={record.determinant_rel_error:.3e}, "
            f"FEM rel err={record.fem_rel_error:.3e}, "
            f"|det|={record.determinant_abs_det:.3e}"
        )
    print("Shape/FEM residual comparison:")
    for record in shape_records:
        print(
            "  "
            f"mu={record.mu:g}, mode={record.mode_number}: "
            f"MAC={record.shape_mac_det_vs_exact:.12f}, "
            f"shape rel L2={record.shape_rel_l2_det_vs_exact:.3e}, "
            f"FEM residual={record.fem_residual_relative:.3e}, "
            f"Rayleigh rel err={record.rayleigh_lambda_rel_error:.3e}"
        )
    print("Mu-invariance spread:")
    for row in summary_rows:
        print(
            "  "
            f"mode={row['mode_number']}: "
            f"det spread={row['determinant_lambda_spread']:.3e}, "
            f"FEM spread={row['fem_lambda_spread']:.3e}"
        )
    print("CSV outputs:")
    for path in paths.values():
        print(f"  {path.relative_to(REPO_ROOT)}")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Check determinant roots/null vectors in the beta=0 single-rod limit."
    )
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--mu-values", type=parse_mu_values, default=DEFAULT_MU_VALUES)
    parser.add_argument("--num-bending-modes", type=int, default=DEFAULT_NUM_BENDING_MODES)
    parser.add_argument("--output-prefix", type=Path, default=DEFAULT_OUTPUT_PREFIX)
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    run_diagnostic(
        epsilon=args.epsilon,
        mu_values=args.mu_values,
        num_bending_modes=args.num_bending_modes,
        n_fem_modes=DEFAULT_N_FEM_MODES,
        output_prefix=args.output_prefix,
    )


if __name__ == "__main__":
    main()
