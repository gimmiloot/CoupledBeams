from __future__ import annotations

from typing import Sequence

import numpy as np


RIGHT_COORDINATE_CHOICES = ("external-to-joint", "joint-to-external")
NEAR_ZERO_NORM = 1e-12
COMPONENT_KEYS = ("u_left", "w_left", "u_right", "w_right")
ROW_LABELS = (
    "kinematic_transverse_compatibility",
    "kinematic_longitudinal_compatibility",
    "slope_compatibility",
    "moment_equilibrium",
    "transverse_global_force_equilibrium",
    "axial_global_force_equilibrium",
)


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if abs(float(denominator)) > NEAR_ZERO_NORM else np.nan


def unique_sorted_roots(values: Sequence[float], *, tolerance: float = 1e-7) -> list[float]:
    roots = sorted(float(value) for value in values if np.isfinite(value) and float(value) > 0.0)
    unique: list[float] = []
    for root in roots:
        if not unique or abs(root - unique[-1]) > tolerance:
            unique.append(root)
    return unique


def analytic_null_vector(matrix: np.ndarray) -> tuple[np.ndarray, float, float]:
    _, singular_values, vh = np.linalg.svd(np.asarray(matrix, dtype=float))
    coeff = vh[-1, :].astype(float)
    if coeff.size:
        pivot = int(np.argmax(np.abs(coeff)))
        if coeff[pivot] < 0.0:
            coeff = -coeff
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
    right_coordinate: str = "joint-to-external",
) -> dict[str, np.ndarray]:
    if right_coordinate not in RIGHT_COORDINATE_CHOICES:
        raise ValueError(f"Unknown right coordinate convention: {right_coordinate}")

    A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]
    xi = np.asarray(s_norm, dtype=float)

    z1 = float(Lambda) * (1.0 - float(mu_value)) * xi
    theta1_arg = float(epsilon) * float(Lambda) ** 2 * (1.0 - float(mu_value)) * xi
    w_left = A1 * (np.cos(z1) - np.cosh(z1)) + B1 * (np.sin(z1) - np.sinh(z1))
    u_left = P1 * np.sin(theta1_arg)

    # This matches the determinant-side reconstruction used by the comparison
    # diagnostic: right arm first in external-to-joint order with negative local
    # argument, then reversed for joint-to-external reporting.
    z2 = -float(Lambda) * (1.0 + float(mu_value)) * xi
    theta2_arg = -float(epsilon) * float(Lambda) ** 2 * (1.0 + float(mu_value)) * xi
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


def concatenate_components(components: dict[str, np.ndarray]) -> np.ndarray:
    return np.concatenate([np.asarray(components[key], dtype=float) for key in COMPONENT_KEYS])


def component_scale(
    components: dict[str, np.ndarray],
    *,
    plot_kind: str,
    normalize: str,
) -> float:
    if normalize == "none":
        return 1.0
    if normalize == "max-transverse":
        value = max(
            float(np.max(np.abs(components["w_left"]))),
            float(np.max(np.abs(components["w_right"]))),
        )
    elif normalize == "max-full":
        left_full = np.sqrt(np.asarray(components["u_left"], dtype=float) ** 2 + np.asarray(components["w_left"], dtype=float) ** 2)
        right_full = np.sqrt(np.asarray(components["u_right"], dtype=float) ** 2 + np.asarray(components["w_right"], dtype=float) ** 2)
        value = max(float(np.max(left_full)), float(np.max(right_full)))
    else:
        raise ValueError(f"Unknown normalization mode for {plot_kind}: {normalize}")
    return max(value, NEAR_ZERO_NORM)


def normalize_components(
    components: dict[str, np.ndarray],
    *,
    plot_kind: str,
    normalize: str,
) -> tuple[dict[str, np.ndarray], float]:
    scale = component_scale(components, plot_kind=plot_kind, normalize=normalize)
    return {key: np.asarray(value, dtype=float) / scale for key, value in components.items()}, scale


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


def slope_value(A: float, B: float, z: float) -> float:
    return float(A * (-np.sin(z) - np.sinh(z)) + B * (np.cos(z) - np.cosh(z)))


def moment_value(A: float, B: float, z: float) -> float:
    return float(A * (-np.cos(z) - np.cosh(z)) + B * (-np.sin(z) - np.sinh(z)))


def shear_row_value(A: float, B: float, z: float, *, Lambda: float, epsilon: float) -> float:
    d3w_dz3 = A * (np.sin(z) - np.sinh(z)) + B * (-np.cos(z) - np.cosh(z))
    return -float(epsilon) * float(Lambda) * float(d3w_dz3)


def analytic_endpoint_quantities(
    Lambda: float,
    *,
    mu_value: float,
    epsilon: float,
    coeff: np.ndarray,
) -> dict[str, float]:
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
        right_coordinate="external-to-joint",
    )
    return {
        "u1_joint": float(components["u_left"][-1]),
        "w1_joint": float(components["w_left"][-1]),
        "slope1_joint": slope_value(A1, B1, x1),
        "moment1_joint": moment_value(A1, B1, x1),
        "shear1_joint": shear_row_value(A1, B1, x1, Lambda=Lambda, epsilon=epsilon),
        "axial_force1_joint": float(P1 * np.cos(th1)),
        "u2_joint": float(components["u_right"][-1]),
        "w2_joint": float(components["w_right"][-1]),
        "slope2_joint": slope_value(A2, B2, -x2),
        "moment2_joint": moment_value(A2, B2, -x2),
        "shear2_joint": shear_row_value(A2, B2, -x2, Lambda=Lambda, epsilon=epsilon),
        "axial_force2_joint": float(P2 * np.cos(-th2)),
        "left_external_w": float(components["w_left"][0]),
        "left_external_slope": slope_value(A1, B1, 0.0),
        "left_external_u": float(components["u_left"][0]),
        "right_external_w": float(components["w_right"][0]),
        "right_external_slope": slope_value(A2, B2, 0.0),
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


def endpoint_consistency_diagnostics(
    matrix: np.ndarray,
    coeff: np.ndarray,
    *,
    Lambda: float,
    beta_rad: float,
    mu_value: float,
    epsilon: float,
) -> dict[str, float]:
    matrix_residual = np.asarray(matrix, dtype=float) @ np.asarray(coeff, dtype=float)
    endpoint = analytic_endpoint_quantities(
        Lambda,
        mu_value=mu_value,
        epsilon=epsilon,
        coeff=coeff,
    )
    field_residual = field_residuals_from_endpoint_quantities(endpoint, beta_rad=beta_rad)
    external_values = [
        endpoint["left_external_w"],
        endpoint["left_external_slope"],
        endpoint["left_external_u"],
        endpoint["right_external_w"],
        endpoint["right_external_slope"],
        endpoint["right_external_u"],
    ]
    return {
        "max_matrix_residual": float(np.max(np.abs(matrix_residual))),
        "max_field_residual": float(np.max(np.abs(field_residual))),
        "matrix_field_residual_max_abs_difference": float(np.max(np.abs(matrix_residual - field_residual))),
        "external_clamp_residual": float(max(abs(float(value)) for value in external_values)),
    }


__all__ = [
    "COMPONENT_KEYS",
    "NEAR_ZERO_NORM",
    "RIGHT_COORDINATE_CHOICES",
    "ROW_LABELS",
    "analytic_arm_energy_diagnostics",
    "analytic_endpoint_quantities",
    "analytic_null_vector",
    "component_scale",
    "concatenate_components",
    "endpoint_consistency_diagnostics",
    "field_residuals_from_endpoint_quantities",
    "normalize_components",
    "reconstruct_analytic_components",
    "safe_ratio",
    "trapezoid_integral",
    "unique_sorted_roots",
]
