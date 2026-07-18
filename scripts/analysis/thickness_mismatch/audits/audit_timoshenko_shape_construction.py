from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    assemble_clamped_coupled_matrix_eta,
)
from scripts.analysis.thickness_mismatch.maps import (  # noqa: E402
    plot_eb_vs_timoshenko_lambda_beta_cases as beta_workflow,
)
from scripts.lib import in_plane_shape_geometry as DISPLAY  # noqa: E402
from scripts.lib import variable_length_timoshenko as TIMO  # noqa: E402
from scripts.lib.analytic_coupled_rods_shapes import analytic_null_vector  # noqa: E402


MODEL_EB = beta_workflow.MODEL_EB
MODEL_TIMO = beta_workflow.MODEL_TIMO
MODELS = (MODEL_EB, MODEL_TIMO)

DEFAULT_BETA_DEG = 45.0
DEFAULT_ETA = 0.0
DEFAULT_OUTPUT_DIR = Path("results") / "timoshenko_shape_construction_audit"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "timoshenko_shape_construction_audit"
DEFAULT_N_ROOTS = 6
DEFAULT_N_POINTS = 801
SMOKE_N_POINTS = 201
FULL_DEFORMATION_FRACTION = 0.16
TRANSVERSE_DEFORMATION_FRACTION = 0.22
DEFAULT_DEFORMATION_SCALE_FRACTION = 0.08
FALLBACK_DEFORMATION_SCALE_FRACTION = 0.05
FIXED_SCALE_FRACTIONS = (0.02, 0.05, 0.10, 0.20)
JOINT_TOL = 1.0e-6
CURVE_SPEED_TOL = 5.0e-2

CONTROL_EPSILON = 0.03
CONTROL_MU = 0.35
CONTROL_SORTED = (1, 2)

SUMMARY_FIELDS = [
    "epsilon",
    "mu",
    "beta_deg",
    "eta",
    "sorted_index",
    "model",
    "Lambda",
    "case_role",
    "full_shape_file",
    "transverse_shape_file",
    "component_file",
    "final_scale_fraction",
    "scale_check_passed",
    "fallback_scale_used",
    "final_figure_file",
    "coefficient_order",
    "rod1_coefficients",
    "rod2_coefficients",
    "rod1_coordinate",
    "rod2_coordinate",
    "rod2_global_transform",
    "plot_object_policy",
    "point_order",
    "rod2_first_point_role",
    "rod2_last_point_role",
    "full_deformation_scale",
    "transverse_deformation_scale",
    "max_abs_u",
    "max_abs_w",
    "axial_energy_fraction",
    "bending_energy_fraction",
    "shear_energy_fraction",
    "classification",
    "bad_x2_plus_l2_max_kinematic_gap",
    "bad_standard_rotation_plot_joint_gap",
    "bad_local_as_global_plot_joint_gap",
    "bad_external_to_joint_connector_length",
    "notes",
]

JOINT_FIELDS = [
    "epsilon",
    "mu",
    "beta_deg",
    "eta",
    "sorted_index",
    "Lambda",
    "gap_w",
    "gap_u",
    "gap_psi",
    "gap_M",
    "gap_Q",
    "gap_N",
    "max_abs_kinematic_gap",
    "max_abs_force_gap",
    "pass_kinematic",
    "pass_force",
    "undeformed_plot_joint_gap",
    "deformed_plot_joint_gap",
    "notes",
]

POINT_ORDER_FIELDS = [
    "epsilon",
    "mu",
    "beta_deg",
    "eta",
    "sorted_index",
    "Lambda",
    "rod_id",
    "local_coordinate_start",
    "local_coordinate_end",
    "is_monotone_in_material_coordinate",
    "plotted_in_material_order",
    "warning",
    "notes",
]

REGULARITY_FIELDS = [
    "epsilon",
    "mu",
    "beta_deg",
    "eta",
    "sorted_index",
    "rod_id",
    "scale_kind",
    "deformation_scale_fraction",
    "deformation_scale",
    "min_abs_dr_plot_ds",
    "max_abs_dr_plot_ds",
    "projected_x_sign_change_count",
    "projected_x_is_monotone",
    "near_overlap_pair_count",
    "min_nonneighbor_distance",
    "visually_folded_at_scale",
    "warning",
    "notes",
]

SCALE_CHECK_FIELDS = [
    "epsilon",
    "mu",
    "beta_deg",
    "eta",
    "sorted_index",
    "Lambda",
    "scale_fraction",
    "min_drds",
    "max_drds",
    "fold_flag",
    "near_overlap_flag",
    "plotted_joint_gap",
    "max_abs_kinematic_gap",
    "max_abs_force_gap",
    "passed_scale_check",
    "notes",
]


@dataclass(frozen=True)
class SuspectSpec:
    epsilon: float
    sorted_index: int
    mu_values: tuple[float, ...]


@dataclass(frozen=True)
class RootPair:
    eb: float
    timo: float
    source: str
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class RodFields:
    x_local: np.ndarray
    u: np.ndarray
    w: np.ndarray
    psi: np.ndarray | None
    gamma: np.ndarray | None
    u_prime: np.ndarray | None
    psi_prime: np.ndarray | None
    w_second: np.ndarray | None


@dataclass(frozen=True)
class ModeResult:
    epsilon: float
    beta_deg: float
    eta: float
    mu: float
    sorted_index: int
    model: str
    Lambda: float
    coeff: np.ndarray
    rod1: RodFields
    rod2: RodFields
    case_role: str
    root_source: str
    warnings: tuple[str, ...]
    sign_factor: float
    energy: dict[str, float]
    joint_row: dict[str, object] | None


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def fmt(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        value_f = float(value)
        if not isfinite(value_f):
            return "nan"
        return f"{value_f:.16e}"
    if isinstance(value, (bool, np.bool_)):
        return "true" if bool(value) else "false"
    return value


def write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})
    return path


def token(value: float, *, min_decimals: int = 2, max_decimals: int = 5) -> str:
    text = f"{float(value):.{max_decimals}f}".rstrip("0").rstrip(".")
    if "." not in text and min_decimals > 0:
        text += "." + "0" * min_decimals
    if "." in text:
        decimals = len(text.split(".", 1)[1])
        if decimals < min_decimals:
            text += "0" * (min_decimals - decimals)
    return text.replace("-", "m").replace(".", "p")


def eps_token(value: float) -> str:
    return token(float(value), min_decimals=2, max_decimals=5)


def mu_token(value: float) -> str:
    return token(float(value), min_decimals=2, max_decimals=5)


def scale_token(value: float) -> str:
    return token(float(value), min_decimals=2, max_decimals=5)


def fixed_scale_fractions(candidate_fraction: float) -> tuple[float, ...]:
    values = list(FIXED_SCALE_FRACTIONS)
    candidate = float(candidate_fraction)
    if not any(abs(candidate - value) <= 1.0e-12 for value in values):
        values.append(candidate)
    return tuple(sorted(values))


def result_scale_key_from_values(epsilon: float, mu: float, sorted_index: int) -> tuple[float, float, int]:
    return (round(float(epsilon), 12), round(float(mu), 12), int(sorted_index))


def result_scale_key(result: ModeResult) -> tuple[float, float, int]:
    return result_scale_key_from_values(result.epsilon, result.mu, result.sorted_index)


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if abs(float(denominator)) > 1.0e-30 else float("nan")


def trapz(values: np.ndarray, coordinates: np.ndarray) -> float:
    integrate = getattr(np, "trapezoid", None)
    if integrate is not None:
        return float(integrate(np.asarray(values, dtype=float), np.asarray(coordinates, dtype=float)))
    y = np.asarray(values, dtype=float)
    x = np.asarray(coordinates, dtype=float)
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def sign_factor_from_components(rod1: RodFields, rod2: RodFields) -> float:
    arrays = [rod1.u, rod1.w, rod2.u, rod2.w]
    if rod1.psi is not None:
        arrays.append(rod1.psi)
    if rod2.psi is not None:
        arrays.append(rod2.psi)
    flat = np.concatenate([np.asarray(array, dtype=float) for array in arrays])
    if flat.size == 0 or not np.any(np.isfinite(flat)):
        return 1.0
    index = int(np.nanargmax(np.abs(flat)))
    value = float(flat[index])
    return -1.0 if value < 0.0 else 1.0


def classify_energy(model: str, axial: float, bending: float, shear: float) -> str:
    bending_shear = float(bending) + float(shear)
    if isfinite(float(axial)) and float(axial) >= 0.70:
        return "longitudinal_dominated"
    if model == MODEL_TIMO and isfinite(bending_shear) and bending_shear >= 0.70:
        return "bending_dominated"
    if model == MODEL_EB and isfinite(float(bending)) and float(bending) >= 0.70:
        return "bending_dominated"
    return "mixed"


def energy_dict(
    *,
    model: str,
    U_axial: float,
    U_bending: float,
    U_shear: float,
    disp_u: float,
    disp_w: float,
) -> dict[str, float]:
    total = float(U_axial) + float(U_bending) + float(U_shear)
    axial = safe_ratio(U_axial, total)
    bending = safe_ratio(U_bending, total)
    shear = safe_ratio(U_shear, total)
    return {
        "U_axial": float(U_axial),
        "U_bending": float(U_bending),
        "U_shear": float(U_shear),
        "U_total": float(total),
        "axial_energy_fraction": float(axial),
        "bending_energy_fraction": float(bending),
        "shear_energy_fraction": float(shear),
        "axial_displacement_fraction": safe_ratio(disp_u, disp_u + disp_w),
        "transverse_displacement_fraction": safe_ratio(disp_w, disp_u + disp_w),
        "classification": classify_energy(model, axial, bending, shear),
    }


class RootProvider:
    def __init__(self, beta_deg: float, eta: float, n_roots: int) -> None:
        self.beta_deg = float(beta_deg)
        self.eta = float(eta)
        self.n_roots = int(n_roots)
        self.cache: dict[tuple[float, float], tuple[np.ndarray, np.ndarray, tuple[str, ...]]] = {}

    def roots(self, epsilon: float, mu: float) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
        key = (round(float(epsilon), 12), round(float(mu), 12))
        if key in self.cache:
            eb, timo, warnings = self.cache[key]
            return eb.copy(), timo.copy(), warnings

        case = beta_workflow.CaseSpec(mu=float(mu), eta=self.eta, epsilon=float(epsilon))
        warnings: list[str] = []

        eb_result = beta_workflow.solve_model(case, self.beta_deg, self.n_roots, MODEL_EB)
        if eb_result.root_count_found < self.n_roots:
            eb_result = beta_workflow.retry_missing_roots(eb_result, case, self.beta_deg, self.n_roots, MODEL_EB)
        warnings.extend(f"EB: {item}" for item in eb_result.warnings)

        finite_eb = [float(root) for root in eb_result.roots if isfinite(float(root))]
        upper_hint = max(finite_eb) if finite_eb else None
        timo_result = beta_workflow.solve_model(
            case,
            self.beta_deg,
            self.n_roots,
            MODEL_TIMO,
            upper_hint=upper_hint,
        )
        if timo_result.root_count_found < self.n_roots:
            timo_result = beta_workflow.retry_missing_roots(
                timo_result,
                case,
                self.beta_deg,
                self.n_roots,
                MODEL_TIMO,
                upper_hint=upper_hint,
            )
        warnings.extend(f"Timoshenko: {item}" for item in timo_result.warnings)

        eb_roots = np.asarray(eb_result.roots, dtype=float)
        timo_roots = np.asarray(timo_result.roots, dtype=float)
        self.cache[key] = (eb_roots.copy(), timo_roots.copy(), tuple(warnings))
        return eb_roots, timo_roots, tuple(warnings)

    def pair(self, epsilon: float, mu: float, sorted_index: int) -> RootPair:
        eb_roots, timo_roots, warnings = self.roots(float(epsilon), float(mu))
        index = int(sorted_index) - 1
        if index < 0 or index >= len(eb_roots) or index >= len(timo_roots):
            raise RuntimeError(f"sorted index {sorted_index} is outside solved root arrays")
        eb = float(eb_roots[index])
        timo = float(timo_roots[index])
        if not (isfinite(eb) and isfinite(timo)):
            raise RuntimeError(
                f"non-finite sorted root for epsilon={epsilon:g}, mu={mu:g}, sorted={sorted_index}"
            )
        return RootPair(eb=eb, timo=timo, source="sorted_global_solver", warnings=warnings)


def eb_mode_result(
    *,
    epsilon: float,
    beta_deg: float,
    eta: float,
    mu: float,
    sorted_index: int,
    Lambda: float,
    n_points: int,
    case_role: str,
    root_source: str,
    root_warnings: Sequence[str],
) -> ModeResult:
    if abs(float(eta)) > 1.0e-14:
        raise ValueError("this focused EB comparison path is restricted to eta=0")
    factors = TIMO.tau_factors(float(mu), float(eta))
    section1 = TIMO.section_from_epsilon_tau(float(epsilon), factors.tau1)
    section2 = TIMO.section_from_epsilon_tau(float(epsilon), factors.tau2)
    l1, l2 = TIMO.segment_lengths(float(mu))
    matrix = assemble_clamped_coupled_matrix_eta(
        float(Lambda),
        float(np.deg2rad(float(beta_deg))),
        float(mu),
        float(epsilon),
        float(eta),
    )
    coeff, _smallest, _ratio = analytic_null_vector(matrix)
    A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]

    x1 = np.linspace(0.0, l1, int(n_points), dtype=float)
    x2 = np.linspace(0.0, -l2, int(n_points), dtype=float)
    z1 = float(Lambda) * x1
    z2 = float(Lambda) * x2
    th1 = float(epsilon) * float(Lambda) ** 2 * x1
    th2 = float(epsilon) * float(Lambda) ** 2 * x2

    w1 = A1 * (np.cos(z1) - np.cosh(z1)) + B1 * (np.sin(z1) - np.sinh(z1))
    u1 = P1 * np.sin(th1)
    u1_prime = P1 * float(epsilon) * float(Lambda) ** 2 * np.cos(th1)
    w1_second = float(Lambda) ** 2 * (
        A1 * (-np.cos(z1) - np.cosh(z1)) + B1 * (-np.sin(z1) - np.sinh(z1))
    )

    w2 = A2 * (np.cos(z2) - np.cosh(z2)) + B2 * (np.sin(z2) - np.sinh(z2))
    u2 = P2 * np.sin(th2)
    u2_prime = P2 * float(epsilon) * float(Lambda) ** 2 * np.cos(th2)
    w2_second = float(Lambda) ** 2 * (
        A2 * (-np.cos(z2) - np.cosh(z2)) + B2 * (-np.sin(z2) - np.sinh(z2))
    )

    rod1 = RodFields(x_local=x1, u=u1, w=w1, psi=None, gamma=None, u_prime=u1_prime, psi_prime=None, w_second=w1_second)
    rod2 = RodFields(x_local=x2, u=u2, w=w2, psi=None, gamma=None, u_prime=u2_prime, psi_prime=None, w_second=w2_second)
    sign_factor = sign_factor_from_components(rod1, rod2)

    x1_abs = np.abs(x1)
    x2_abs = np.abs(x2)
    U_axial = 0.5 * TIMO.E * section1.area * trapz(u1_prime**2, x1_abs)
    U_axial += 0.5 * TIMO.E * section2.area * trapz(u2_prime**2, x2_abs)
    U_bending = 0.5 * section1.bending_stiffness * trapz(w1_second**2, x1_abs)
    U_bending += 0.5 * section2.bending_stiffness * trapz(w2_second**2, x2_abs)
    disp_u = trapz(u1**2, x1_abs) + trapz(u2**2, x2_abs)
    disp_w = trapz(w1**2, x1_abs) + trapz(w2**2, x2_abs)

    return ModeResult(
        epsilon=float(epsilon),
        beta_deg=float(beta_deg),
        eta=float(eta),
        mu=float(mu),
        sorted_index=int(sorted_index),
        model=MODEL_EB,
        Lambda=float(Lambda),
        coeff=np.asarray(coeff, dtype=float),
        rod1=rod1,
        rod2=rod2,
        case_role=case_role,
        root_source=root_source,
        warnings=tuple(root_warnings),
        sign_factor=sign_factor,
        energy=energy_dict(
            model=MODEL_EB,
            U_axial=max(float(U_axial), 0.0),
            U_bending=max(float(U_bending), 0.0),
            U_shear=0.0,
            disp_u=max(float(disp_u), 0.0),
            disp_w=max(float(disp_w), 0.0),
        ),
        joint_row=None,
    )


def endpoint_value(rod: dict[str, object], name: str, index: int) -> float:
    return float(np.asarray(rod[name], dtype=float)[index])


def timo_joint_row(
    *,
    epsilon: float,
    beta_deg: float,
    eta: float,
    mu: float,
    sorted_index: int,
    Lambda: float,
    fields: dict[str, object],
    coeff: np.ndarray,
    sign_factor: float,
) -> dict[str, object]:
    rod1 = fields["rod1"]
    rod2 = fields["rod2"]
    beta_rad = np.deg2rad(float(beta_deg))
    cb = float(np.cos(beta_rad))
    sb = float(np.sin(beta_rad))
    l1, l2 = TIMO.segment_lengths(float(mu))

    i1 = -1
    i2 = -1
    w1 = endpoint_value(rod1, "w", i1)
    u1 = endpoint_value(rod1, "u", i1)
    psi1 = endpoint_value(rod1, "psi", i1)
    m1 = endpoint_value(rod1, "M", i1)
    q1 = endpoint_value(rod1, "Q", i1)
    n1 = endpoint_value(rod1, "N", i1)

    w2 = endpoint_value(rod2, "w", i2)
    u2 = endpoint_value(rod2, "u", i2)
    psi2 = endpoint_value(rod2, "psi", i2)
    m2 = endpoint_value(rod2, "M", i2)
    q2 = endpoint_value(rod2, "Q", i2)
    n2 = endpoint_value(rod2, "N", i2)

    gap_w = w1 - cb * w2 + sb * u2
    gap_u = u1 - sb * w2 - cb * u2
    gap_psi = psi1 - psi2
    gap_m = m1 - m2
    gap_q = -q1 + cb * q2 - sb * n2
    gap_n = n1 - sb * q2 - cb * n2
    max_abs_kinematic_gap = max(abs(gap_w), abs(gap_u), abs(gap_psi))
    max_abs_force_gap = max(abs(gap_m), abs(gap_q), abs(gap_n))

    base_joint_1 = np.array([l1, 0.0], dtype=float)
    base_joint_2 = np.array([l1, 0.0], dtype=float)
    undeformed_plot_joint_gap = float(np.linalg.norm(base_joint_1 - base_joint_2))

    full_scale = deformation_scale_for_arrays(
        mu=float(mu),
        rod1_u=np.asarray(rod1["u"], dtype=float),
        rod1_w=np.asarray(rod1["w"], dtype=float),
        rod2_u=np.asarray(rod2["u"], dtype=float),
        rod2_w=np.asarray(rod2["w"], dtype=float),
        transverse_only=False,
    )
    sign = float(sign_factor)
    dx1, dy1 = DISPLAY.rod1_local_displacement_to_display(
        sign * np.array([u1], dtype=float),
        sign * np.array([w1], dtype=float),
    )
    dx2, dy2 = DISPLAY.rod2_local_displacement_to_display(
        sign * np.array([u2], dtype=float),
        sign * np.array([w2], dtype=float),
        beta_deg=float(beta_deg),
    )
    disp1 = np.array([dx1[0], dy1[0]], dtype=float)
    disp2 = np.array([dx2[0], dy2[0]], dtype=float)
    deformed_plot_joint_gap = float(np.linalg.norm((base_joint_1 + full_scale * disp1) - (base_joint_2 + full_scale * disp2)))

    return {
        "epsilon": float(epsilon),
        "mu": float(mu),
        "beta_deg": float(beta_deg),
        "eta": float(eta),
        "sorted_index": int(sorted_index),
        "Lambda": float(Lambda),
        "gap_w": float(gap_w),
        "gap_u": float(gap_u),
        "gap_psi": float(gap_psi),
        "gap_M": float(gap_m),
        "gap_Q": float(gap_q),
        "gap_N": float(gap_n),
        "max_abs_kinematic_gap": float(max_abs_kinematic_gap),
        "max_abs_force_gap": float(max_abs_force_gap),
        "pass_kinematic": bool(max_abs_kinematic_gap <= JOINT_TOL),
        "pass_force": bool(max_abs_force_gap <= JOINT_TOL),
        "undeformed_plot_joint_gap": undeformed_plot_joint_gap,
        "deformed_plot_joint_gap": deformed_plot_joint_gap,
        "notes": (
            "rod1 endpoint x1=+l1; rod2 endpoint x2=-l2; "
            "display joint uses d1=(u1,-w1), d2=(c*u2+s*w2,s*u2-c*w2); "
            f"coeff_norm={float(np.linalg.norm(np.asarray(coeff, dtype=float))):.6g}"
        ),
    }


def timo_mode_result(
    *,
    epsilon: float,
    beta_deg: float,
    eta: float,
    mu: float,
    sorted_index: int,
    Lambda: float,
    n_points: int,
    case_role: str,
    root_source: str,
    root_warnings: Sequence[str],
) -> ModeResult:
    mode = TIMO.timo_mode_coefficients(
        float(Lambda),
        float(beta_deg),
        float(mu),
        float(epsilon),
        float(eta),
    )
    fields = TIMO.timo_mode_fields(
        float(Lambda),
        float(beta_deg),
        float(mu),
        float(epsilon),
        float(eta),
        coeff=mode.coeff,
        n_points=int(n_points),
    )
    rod1_raw = fields["rod1"]
    rod2_raw = fields["rod2"]

    gamma1 = np.asarray(rod1_raw["w_prime"], dtype=float) - np.asarray(rod1_raw["psi"], dtype=float)
    gamma2 = np.asarray(rod2_raw["w_prime"], dtype=float) - np.asarray(rod2_raw["psi"], dtype=float)
    rod1 = RodFields(
        x_local=np.asarray(rod1_raw["x"], dtype=float),
        u=np.asarray(rod1_raw["u"], dtype=float),
        w=np.asarray(rod1_raw["w"], dtype=float),
        psi=np.asarray(rod1_raw["psi"], dtype=float),
        gamma=gamma1,
        u_prime=np.asarray(rod1_raw["u_prime"], dtype=float),
        psi_prime=np.asarray(rod1_raw["psi_prime"], dtype=float),
        w_second=None,
    )
    rod2 = RodFields(
        x_local=np.asarray(rod2_raw["x"], dtype=float),
        u=np.asarray(rod2_raw["u"], dtype=float),
        w=np.asarray(rod2_raw["w"], dtype=float),
        psi=np.asarray(rod2_raw["psi"], dtype=float),
        gamma=gamma2,
        u_prime=np.asarray(rod2_raw["u_prime"], dtype=float),
        psi_prime=np.asarray(rod2_raw["psi_prime"], dtype=float),
        w_second=None,
    )
    sign_factor = sign_factor_from_components(rod1, rod2)
    joint = timo_joint_row(
        epsilon=epsilon,
        beta_deg=beta_deg,
        eta=eta,
        mu=mu,
        sorted_index=sorted_index,
        Lambda=Lambda,
        fields=fields,
        coeff=mode.coeff,
        sign_factor=sign_factor,
    )

    partition = TIMO.timo_energy_partition(
        float(Lambda),
        float(beta_deg),
        float(mu),
        float(epsilon),
        float(eta),
        coeff=mode.coeff,
        n_points=int(n_points),
    )
    x1_abs = np.abs(rod1.x_local)
    x2_abs = np.abs(rod2.x_local)
    disp_u = trapz(rod1.u**2, x1_abs) + trapz(rod2.u**2, x2_abs)
    disp_w = trapz(rod1.w**2, x1_abs) + trapz(rod2.w**2, x2_abs)
    warnings = list(root_warnings)
    warnings.extend(mode.warnings)
    warnings.extend(fields.get("warnings", ()))
    warnings = list(dict.fromkeys(item for item in warnings if item))

    return ModeResult(
        epsilon=float(epsilon),
        beta_deg=float(beta_deg),
        eta=float(eta),
        mu=float(mu),
        sorted_index=int(sorted_index),
        model=MODEL_TIMO,
        Lambda=float(Lambda),
        coeff=np.asarray(mode.coeff, dtype=float),
        rod1=rod1,
        rod2=rod2,
        case_role=case_role,
        root_source=root_source,
        warnings=tuple(warnings),
        sign_factor=sign_factor,
        energy=energy_dict(
            model=MODEL_TIMO,
            U_axial=max(float(partition["U_a_total"]), 0.0),
            U_bending=max(float(partition["U_b_total"]), 0.0),
            U_shear=max(float(partition["U_s_total"]), 0.0),
            disp_u=max(float(disp_u), 0.0),
            disp_w=max(float(disp_w), 0.0),
        ),
        joint_row=joint,
    )


def max_abs_u_w(result: ModeResult) -> tuple[float, float]:
    max_u = max(float(np.max(np.abs(result.rod1.u))), float(np.max(np.abs(result.rod2.u))))
    max_w = max(float(np.max(np.abs(result.rod1.w))), float(np.max(np.abs(result.rod2.w))))
    return max_u, max_w


def max_full_displacement(result: ModeResult) -> float:
    full1 = np.sqrt(result.rod1.u**2 + result.rod1.w**2)
    full2 = np.sqrt(result.rod2.u**2 + result.rod2.w**2)
    return max(float(np.max(full1)), float(np.max(full2)), 1.0e-14)


def deformation_scale_for_arrays(
    *,
    mu: float,
    rod1_u: np.ndarray,
    rod1_w: np.ndarray,
    rod2_u: np.ndarray,
    rod2_w: np.ndarray,
    transverse_only: bool,
) -> float:
    l1, l2 = TIMO.segment_lengths(float(mu))
    length_scale = max(l1, l2, 1.0)
    if transverse_only:
        denom = max(float(np.max(np.abs(rod1_w))), float(np.max(np.abs(rod2_w))), 1.0e-14)
        return TRANSVERSE_DEFORMATION_FRACTION * length_scale / denom
    full1 = np.sqrt(np.asarray(rod1_u, dtype=float) ** 2 + np.asarray(rod1_w, dtype=float) ** 2)
    full2 = np.sqrt(np.asarray(rod2_u, dtype=float) ** 2 + np.asarray(rod2_w, dtype=float) ** 2)
    denom = max(float(np.max(full1)), float(np.max(full2)), 1.0e-14)
    return FULL_DEFORMATION_FRACTION * length_scale / denom


def deformation_scale(result: ModeResult, *, transverse_only: bool) -> float:
    return deformation_scale_for_arrays(
        mu=result.mu,
        rod1_u=result.rod1.u,
        rod1_w=result.rod1.w,
        rod2_u=result.rod2.u,
        rod2_w=result.rod2.w,
        transverse_only=transverse_only,
    )


def base_coordinates(mu: float, beta_deg: float, n_points: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    l1, l2 = TIMO.segment_lengths(float(mu))
    x1_grid = np.linspace(0.0, l1, int(n_points), dtype=float)
    x2_grid = np.linspace(-l2, 0.0, int(n_points), dtype=float)
    zeros1 = np.zeros_like(x1_grid)
    zeros2 = np.zeros_like(x2_grid)
    rod1 = DISPLAY.rod1_local_fields_to_display(x1_grid, zeros1, zeros1, scale=0.0)
    rod2 = DISPLAY.rod2_local_fields_to_display(
        x2_grid,
        zeros2,
        zeros2,
        l2=l2,
        x_joint=l1,
        beta_deg=float(beta_deg),
        scale=0.0,
    )
    return rod1.x_base, rod1.y_base, rod2.x_base, rod2.y_base


def plot_arrays(result: ModeResult) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    sign = float(result.sign_factor)
    return (
        sign * np.asarray(result.rod1.u, dtype=float),
        sign * np.asarray(result.rod1.w, dtype=float),
        sign * np.asarray(result.rod2.u, dtype=float)[::-1],
        sign * np.asarray(result.rod2.w, dtype=float)[::-1],
    )


def fixed_scale_from_fraction(mu: float, fraction: float) -> float:
    l1, l2 = TIMO.segment_lengths(float(mu))
    return float(fraction) * (l1 + l2)


def deformed_coordinates_with_scale(
    result: ModeResult,
    *,
    transverse_only: bool,
    scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    u1, w1, u2, w2 = plot_arrays(result)
    u1_plot = np.zeros_like(u1) if transverse_only else u1
    u2_plot = np.zeros_like(u2) if transverse_only else u2
    _l1, l2 = TIMO.segment_lengths(result.mu)
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
        u1_plot,
        w1,
        scale=float(scale),
    )
    rod2 = rod2_mapper(
        result.rod2.x_local[::-1],
        u2_plot,
        w2,
        l2=l2,
        x_joint=float(result.rod1.x_local[-1]),
        beta_deg=result.beta_deg,
        scale=float(scale),
    )
    return rod1.x_deformed, rod1.y_deformed, rod2.x_deformed, rod2.y_deformed


def deformed_coordinates(
    result: ModeResult,
    *,
    transverse_only: bool,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    return deformed_coordinates_with_scale(
        result,
        transverse_only=transverse_only,
        scale=deformation_scale(result, transverse_only=transverse_only),
    )


def axis_limits(results: Sequence[ModeResult]) -> tuple[tuple[float, float], tuple[float, float]]:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for result in results:
        x1, y1, x2, y2 = base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
        xs.extend([x1, x2])
        ys.extend([y1, y2])
        for transverse_only in (False, True):
            xd1, yd1, xd2, yd2 = deformed_coordinates(result, transverse_only=transverse_only)
            xs.extend([xd1, xd2])
            ys.extend([yd1, yd2])
    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    x_span = max(float(np.max(x_all) - np.min(x_all)), 1.0)
    y_span = max(float(np.max(y_all) - np.min(y_all)), 0.5)
    return (
        (float(np.min(x_all) - 0.07 * x_span), float(np.max(x_all) + 0.07 * x_span)),
        (float(np.min(y_all) - 0.15 * y_span), float(np.max(y_all) + 0.15 * y_span)),
    )


def draw_centerline(
    ax: plt.Axes,
    result: ModeResult,
    *,
    transverse_only: bool,
    limits: tuple[tuple[float, float], tuple[float, float]],
    labels: bool,
) -> None:
    x1, y1, x2, y2 = base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
    xd1, yd1, xd2, yd2 = deformed_coordinates(result, transverse_only=transverse_only)
    color = "#1f77b4" if result.model == MODEL_EB else "#d55e00"
    ax.plot(x1, y1, color="0.72", linestyle="--", linewidth=1.0, label="undeformed" if labels else None)
    ax.plot(x2, y2, color="0.72", linestyle="--", linewidth=1.0)
    ax.plot(xd1, yd1, color=color, linewidth=2.0, label=result.model if labels else None)
    ax.plot(xd2, yd2, color=color, linewidth=2.0)
    ax.scatter([x1[-1]], [y1[-1]], color="black", s=14, zorder=5, label="joint" if labels else None)
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.88", linewidth=0.6)


def plot_single_centerline(result: ModeResult, path: Path, *, transverse_only: bool) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    limits = axis_limits([result])
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    draw_centerline(ax, result, transverse_only=transverse_only, limits=limits, labels=True)
    scale = deformation_scale(result, transverse_only=transverse_only)
    mode_text = "transverse-only" if transverse_only else "full"
    energy = result.energy
    ax.set_title(
        (
            f"{result.model} {mode_text}, eps={result.epsilon:g}, mu={result.mu:g}, "
            f"sorted {result.sorted_index}\n"
            f"Lambda={result.Lambda:.7g}, axial={energy['axial_energy_fraction']:.3f}, "
            f"bending={energy['bending_energy_fraction']:.3f}, shear={energy['shear_energy_fraction']:.3f}, "
            f"scale={scale:.3g}"
        ),
        fontsize=9.2,
    )
    ax.set_xlabel("global x")
    ax.set_ylabel("global y")
    ax.legend(loc="best", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def timo_component_arrays(result: ModeResult) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    if result.model != MODEL_TIMO:
        raise ValueError("component plots are only defined for Timoshenko modes")
    sign = float(result.sign_factor)
    return {
        "u": (sign * result.rod1.u, sign * result.rod2.u),
        "w": (sign * result.rod1.w, sign * result.rod2.w),
        "psi": (sign * np.asarray(result.rod1.psi, dtype=float), sign * np.asarray(result.rod2.psi, dtype=float)),
        "gamma": (sign * np.asarray(result.rod1.gamma, dtype=float), sign * np.asarray(result.rod2.gamma, dtype=float)),
    }


def plot_timo_components(result: ModeResult, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = timo_component_arrays(result)
    fig, axes = plt.subplots(4, 1, figsize=(7.0, 7.0), sharex=True)
    for ax, name in zip(axes, ("u", "w", "psi", "gamma")):
        values1, values2 = arrays[name]
        ax.plot(result.rod1.x_local, values1, color="#0072b2", linewidth=1.6, label="rod 1")
        ax.plot(result.rod2.x_local, values2, color="#d55e00", linewidth=1.6, label="rod 2")
        ax.scatter(
            [result.rod1.x_local[-1], result.rod2.x_local[-1]],
            [values1[-1], values2[-1]],
            color="black",
            s=14,
            zorder=5,
            label="joint" if name == "u" else None,
        )
        ax.axhline(0.0, color="0.8", linewidth=0.6)
        ax.grid(True, color="0.90", linewidth=0.55)
        ax.set_ylabel(name)
    axes[0].legend(loc="best", fontsize=8, frameon=False)
    axes[-1].set_xlabel("signed local coordinate x_i")
    fig.suptitle(
        (
            f"Timoshenko components, eps={result.epsilon:g}, mu={result.mu:g}, "
            f"sorted {result.sorted_index}, Lambda={result.Lambda:.7g}"
        ),
        fontsize=10,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.965))
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_suspect_shape_grid(output_dir: Path, spec: SuspectSpec, results: Sequence[ModeResult]) -> Path:
    selected = [
        result
        for result in results
        if result.case_role == "suspect"
        and abs(result.epsilon - spec.epsilon) <= 1.0e-12
        and result.sorted_index == spec.sorted_index
    ]
    by_key = {(round(result.mu, 12), result.model): result for result in selected}
    mu_values = list(spec.mu_values)
    output = output_dir / f"suspect_shapes_eps{eps_token(spec.epsilon)}_sorted{spec.sorted_index}_grid_corrected.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    limits = axis_limits(selected)
    columns = [
        (MODEL_EB, False, "EB full"),
        (MODEL_EB, True, "EB w-only"),
        (MODEL_TIMO, False, "Timo full"),
        (MODEL_TIMO, True, "Timo w-only"),
    ]
    fig, axes = plt.subplots(len(mu_values), len(columns), figsize=(13.0, 2.35 * len(mu_values) + 0.7), squeeze=False)
    for row_index, mu in enumerate(mu_values):
        for col_index, (model, transverse_only, label) in enumerate(columns):
            ax = axes[row_index, col_index]
            result = by_key[(round(float(mu), 12), model)]
            draw_centerline(ax, result, transverse_only=transverse_only, limits=limits, labels=False)
            energy = result.energy
            ax.set_title(
                (
                    f"{label}, mu={float(mu):g}, Lambda={result.Lambda:.5g}\n"
                    f"axial={energy['axial_energy_fraction']:.3f}, class={energy['classification']}"
                ),
                fontsize=8.2,
            )
            ax.tick_params(labelsize=7)
            if row_index == len(mu_values) - 1:
                ax.set_xlabel("global x", fontsize=8)
            if col_index == 0:
                ax.set_ylabel("global y", fontsize=8)
    fig.suptitle(
        f"Corrected sorted shapes: eps={spec.epsilon:g}, beta=45 deg, eta=0, sorted {spec.sorted_index}",
        fontsize=11,
        y=0.992,
    )
    handles = [
        Line2D([0], [0], color="0.72", linestyle="--", linewidth=1.0, label="undeformed rods"),
        Line2D([0], [0], color="#1f77b4", linewidth=2.0, label="Euler-Bernoulli"),
        Line2D([0], [0], color="#d55e00", linewidth=2.0, label="Timoshenko"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.tight_layout(rect=(0.0, 0.035, 1.0, 0.955), h_pad=0.8, w_pad=0.55)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_suspect_component_grid(output_dir: Path, spec: SuspectSpec, results: Sequence[ModeResult]) -> Path:
    selected = [
        result
        for result in results
        if result.case_role == "suspect"
        and result.model == MODEL_TIMO
        and abs(result.epsilon - spec.epsilon) <= 1.0e-12
        and result.sorted_index == spec.sorted_index
    ]
    by_mu = {round(result.mu, 12): result for result in selected}
    output = output_dir / f"suspect_components_eps{eps_token(spec.epsilon)}_sorted{spec.sorted_index}_grid.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    component_names = ("u", "w", "psi", "gamma")
    fig, axes = plt.subplots(len(spec.mu_values), len(component_names), figsize=(12.5, 2.0 * len(spec.mu_values) + 0.75), squeeze=False)
    for row_index, mu in enumerate(spec.mu_values):
        result = by_mu[round(float(mu), 12)]
        arrays = timo_component_arrays(result)
        for col_index, name in enumerate(component_names):
            ax = axes[row_index, col_index]
            values1, values2 = arrays[name]
            ax.plot(result.rod1.x_local, values1, color="#0072b2", linewidth=1.25)
            ax.plot(result.rod2.x_local, values2, color="#d55e00", linewidth=1.25)
            ax.scatter(
                [result.rod1.x_local[-1], result.rod2.x_local[-1]],
                [values1[-1], values2[-1]],
                color="black",
                s=10,
                zorder=5,
            )
            ax.axhline(0.0, color="0.82", linewidth=0.55)
            ax.grid(True, color="0.91", linewidth=0.5)
            ax.tick_params(labelsize=7)
            if row_index == 0:
                ax.set_title(name, fontsize=9)
            if col_index == 0:
                ax.set_ylabel(f"mu={float(mu):g}", fontsize=8)
            if row_index == len(spec.mu_values) - 1:
                ax.set_xlabel("signed x_i", fontsize=8)
    handles = [
        Line2D([0], [0], color="#0072b2", linewidth=1.4, label="rod 1"),
        Line2D([0], [0], color="#d55e00", linewidth=1.4, label="rod 2"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False, fontsize=8)
    fig.suptitle(
        f"Timoshenko components: eps={spec.epsilon:g}, beta=45 deg, eta=0, sorted {spec.sorted_index}",
        fontsize=11,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.95), h_pad=0.55, w_pad=0.55)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def suspect_timo_results(results: Sequence[ModeResult], spec: SuspectSpec) -> list[ModeResult]:
    selected = [
        result
        for result in results
        if result.case_role == "suspect"
        and result.model == MODEL_TIMO
        and abs(result.epsilon - spec.epsilon) <= 1.0e-12
        and result.sorted_index == spec.sorted_index
    ]
    by_mu = {round(result.mu, 12): result for result in selected}
    return [by_mu[round(float(mu), 12)] for mu in spec.mu_values]


def fixed_scale_axis_limits(results: Sequence[ModeResult], fractions: Sequence[float]) -> tuple[tuple[float, float], tuple[float, float]]:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for result in results:
        x1, y1, x2, y2 = base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
        xs.extend([x1, x2])
        ys.extend([y1, y2])
        for fraction in fractions:
            scale = fixed_scale_from_fraction(result.mu, float(fraction))
            xd1, yd1, xd2, yd2 = deformed_coordinates_with_scale(result, transverse_only=False, scale=scale)
            xs.extend([xd1, xd2])
            ys.extend([yd1, yd2])
    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    x_span = max(float(np.max(x_all) - np.min(x_all)), 1.0)
    y_span = max(float(np.max(y_all) - np.min(y_all)), 0.5)
    return (
        (float(np.min(x_all) - 0.07 * x_span), float(np.max(x_all) + 0.07 * x_span)),
        (float(np.min(y_all) - 0.15 * y_span), float(np.max(y_all) + 0.15 * y_span)),
    )


def draw_fixed_scale_timo_centerline(
    ax: plt.Axes,
    result: ModeResult,
    *,
    fraction: float,
    limits: tuple[tuple[float, float], tuple[float, float]],
) -> None:
    x1, y1, x2, y2 = base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
    scale = fixed_scale_from_fraction(result.mu, float(fraction))
    xd1, yd1, xd2, yd2 = deformed_coordinates_with_scale(result, transverse_only=False, scale=scale)
    ax.plot(x1, y1, color="0.72", linestyle="--", linewidth=1.0)
    ax.plot(x2, y2, color="0.72", linestyle="--", linewidth=1.0)
    ax.plot(xd1, yd1, color="#d55e00", linewidth=1.8)
    ax.plot(xd2, yd2, color="#d55e00", linewidth=1.8)
    ax.scatter([x1[-1]], [y1[-1]], color="black", s=12, zorder=5)
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.89", linewidth=0.55)


def plot_suspect_fixed_scale_shape_grid(
    output_dir: Path,
    spec: SuspectSpec,
    results: Sequence[ModeResult],
    *,
    fraction: float,
    mu_values: Sequence[float] | None = None,
    suffix: str = "",
) -> Path:
    selected_all = suspect_timo_results(results, spec)
    wanted_mu = (
        {round(float(mu), 12) for mu in mu_values}
        if mu_values is not None
        else {round(float(result.mu), 12) for result in selected_all}
    )
    selected = [result for result in selected_all if round(float(result.mu), 12) in wanted_mu]
    output = output_dir / (
        f"suspect_shapes_eps{eps_token(spec.epsilon)}_sorted{spec.sorted_index}"
        f"_grid_scale{scale_token(fraction)}{suffix}.png"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    limits = fixed_scale_axis_limits(selected, (float(fraction),))
    fig, axes = plt.subplots(len(selected), 1, figsize=(7.3, 2.25 * len(selected) + 0.95), squeeze=False)
    for row_index, result in enumerate(selected):
        ax = axes[row_index, 0]
        draw_fixed_scale_timo_centerline(ax, result, fraction=float(fraction), limits=limits)
        energy = result.energy
        ax.set_title(
            (
                f"mu={result.mu:g}, Lambda={result.Lambda:.6g}, deformation scale = {float(fraction):.2f}\n"
                f"axial={energy['axial_energy_fraction']:.3f}, class={energy['classification']}"
            ),
            fontsize=8.4,
        )
        ax.tick_params(labelsize=7)
        ax.set_xlabel("global x", fontsize=8)
        ax.set_ylabel("global y", fontsize=8)
    handles = [
        Line2D([0], [0], color="0.72", linestyle="--", linewidth=1.0, label="undeformed rods"),
        Line2D([0], [0], color="#d55e00", linewidth=1.8, label="Timoshenko fixed-scale"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.suptitle(
        (
            f"Timoshenko full centerline: eps={spec.epsilon:g}, beta=45 deg, eta=0, "
            f"sorted {spec.sorted_index}, deformation scale = {float(fraction):.2f}"
        ),
        fontsize=11,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.045, 1.0, 0.94), h_pad=0.65, w_pad=0.55)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_scale_comparison_grid(
    output_dir: Path,
    spec: SuspectSpec,
    results: Sequence[ModeResult],
    *,
    fractions: Sequence[float],
) -> Path:
    selected = suspect_timo_results(results, spec)
    output = output_dir / (
        f"suspect_timoshenko_scale_{scale_token(fractions[0])}_vs_{scale_token(fractions[1])}"
        f"_eps{eps_token(spec.epsilon)}_sorted{spec.sorted_index}.png"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    limits = fixed_scale_axis_limits(selected, fractions)
    fig, axes = plt.subplots(
        len(selected),
        len(fractions),
        figsize=(8.6, 2.2 * len(selected) + 0.9),
        squeeze=False,
    )
    for row_index, result in enumerate(selected):
        for col_index, fraction in enumerate(fractions):
            ax = axes[row_index, col_index]
            draw_fixed_scale_timo_centerline(ax, result, fraction=float(fraction), limits=limits)
            ax.set_title(
                f"mu={result.mu:g}, deformation scale = {float(fraction):.2f}",
                fontsize=8.2,
            )
            ax.tick_params(labelsize=7)
            if row_index == len(selected) - 1:
                ax.set_xlabel("global x", fontsize=8)
            if col_index == 0:
                ax.set_ylabel("global y", fontsize=8)
    handles = [
        Line2D([0], [0], color="0.72", linestyle="--", linewidth=1.0, label="undeformed rods"),
        Line2D([0], [0], color="#d55e00", linewidth=1.8, label="Timoshenko fixed-scale"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.suptitle(
        (
            f"Timoshenko scale comparison: eps={spec.epsilon:g}, beta=45 deg, eta=0, "
            f"sorted {spec.sorted_index}"
        ),
        fontsize=11,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.05, 1.0, 0.945), h_pad=0.55, w_pad=0.55)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_scale_sensitivity_grid(
    output_dir: Path,
    spec: SuspectSpec,
    results: Sequence[ModeResult],
    fixed_fractions: Sequence[float],
) -> Path:
    selected = suspect_timo_results(results, spec)
    output = output_dir / f"suspect_timoshenko_scale_sensitivity_eps{eps_token(spec.epsilon)}_sorted{spec.sorted_index}.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    limits = fixed_scale_axis_limits(selected, fixed_fractions)
    fig, axes = plt.subplots(
        len(selected),
        len(fixed_fractions),
        figsize=(3.15 * len(fixed_fractions), 2.25 * len(selected) + 0.75),
        squeeze=False,
    )
    for row_index, result in enumerate(selected):
        for col_index, fraction in enumerate(fixed_fractions):
            ax = axes[row_index, col_index]
            draw_fixed_scale_timo_centerline(ax, result, fraction=float(fraction), limits=limits)
            energy = result.energy
            ax.set_title(
                (
                    f"mu={result.mu:g}, scale={float(fraction):.2f}\n"
                    f"axial={energy['axial_energy_fraction']:.3f}, class={energy['classification']}"
                ),
                fontsize=8.2,
            )
            ax.tick_params(labelsize=7)
            if row_index == len(selected) - 1:
                ax.set_xlabel("global x", fontsize=8)
            if col_index == 0:
                ax.set_ylabel("global y", fontsize=8)
    handles = [
        Line2D([0], [0], color="0.72", linestyle="--", linewidth=1.0, label="undeformed rods"),
        Line2D([0], [0], color="#d55e00", linewidth=1.8, label="Timoshenko fixed-scale"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.suptitle(
        (
            f"Timoshenko full centerline scale sensitivity: eps={spec.epsilon:g}, "
            f"beta=45 deg, eta=0, sorted {spec.sorted_index}"
        ),
        fontsize=11,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.045, 1.0, 0.95), h_pad=0.55, w_pad=0.55)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def plotted_local_coordinate(result: ModeResult, rod_id: int) -> np.ndarray:
    if int(rod_id) == 1:
        return np.asarray(result.rod1.x_local, dtype=float)
    if int(rod_id) == 2:
        return np.asarray(result.rod2.x_local, dtype=float)[::-1]
    raise ValueError("rod_id must be 1 or 2")


def point_order_rows(results: Sequence[ModeResult]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        if result.model != MODEL_TIMO:
            continue
        for rod_id in (1, 2):
            local = plotted_local_coordinate(result, rod_id)
            diffs = np.diff(local)
            monotone = bool(np.all(diffs >= -1.0e-12) or np.all(diffs <= 1.0e-12))
            warning = "" if monotone else "local coordinate is not monotone in plotted order"
            rows.append(
                {
                    "epsilon": result.epsilon,
                    "mu": result.mu,
                    "beta_deg": result.beta_deg,
                    "eta": result.eta,
                    "sorted_index": result.sorted_index,
                    "Lambda": result.Lambda,
                    "rod_id": int(rod_id),
                    "local_coordinate_start": float(local[0]),
                    "local_coordinate_end": float(local[-1]),
                    "is_monotone_in_material_coordinate": monotone,
                    "plotted_in_material_order": monotone,
                    "warning": warning,
                    "notes": (
                        "rod plotted as separate line object; no global-x sort; "
                        "rod2 uses reversed helper order so geometry starts at the joint"
                        if int(rod_id) == 2
                        else "rod plotted as separate line object; no global-x sort"
                    ),
                }
            )
    return rows


def plotted_rod_curve(
    result: ModeResult,
    *,
    rod_id: int,
    scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xd1, yd1, xd2, yd2 = deformed_coordinates_with_scale(result, transverse_only=False, scale=float(scale))
    local = plotted_local_coordinate(result, rod_id)
    if int(rod_id) == 1:
        return local, xd1, yd1
    if int(rod_id) == 2:
        return local, xd2, yd2
    raise ValueError("rod_id must be 1 or 2")


def projected_x_sign_change_count(x_values: np.ndarray) -> tuple[int, bool]:
    dx = np.diff(np.asarray(x_values, dtype=float))
    tol = 1.0e-10 * max(1.0, float(np.max(np.abs(x_values))))
    nonzero = dx[np.abs(dx) > tol]
    if nonzero.size <= 1:
        return 0, True
    signs = np.sign(nonzero)
    count = int(np.sum(signs[1:] * signs[:-1] < 0.0))
    monotone = bool(np.all(nonzero >= 0.0) or np.all(nonzero <= 0.0))
    return count, monotone


def near_overlap_stats(x_values: np.ndarray, y_values: np.ndarray, total_length: float) -> tuple[int, float]:
    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    step = max(1, int(np.ceil(len(x) / 260)))
    points = np.column_stack([x[::step], y[::step]])
    n_points = points.shape[0]
    if n_points < 6:
        return 0, float("nan")
    delta = points[:, None, :] - points[None, :, :]
    distances = np.sqrt(np.sum(delta * delta, axis=2))
    index = np.arange(n_points)
    mask = np.abs(index[:, None] - index[None, :]) > 10
    upper = np.triu(np.ones((n_points, n_points), dtype=bool), 1)
    valid = mask & upper
    if not np.any(valid):
        return 0, float("nan")
    selected = distances[valid]
    min_distance = float(np.min(selected))
    threshold = 0.005 * float(total_length)
    count = int(np.sum(selected < threshold))
    return count, min_distance


def regularity_row(
    result: ModeResult,
    *,
    rod_id: int,
    scale_kind: str,
    deformation_scale_fraction: float,
    deformation_scale_value: float,
) -> dict[str, object]:
    local, x_plot, y_plot = plotted_rod_curve(result, rod_id=int(rod_id), scale=float(deformation_scale_value))
    ds = np.abs(np.diff(local))
    dr = np.sqrt(np.diff(x_plot) ** 2 + np.diff(y_plot) ** 2)
    valid = ds > 1.0e-14
    speeds = dr[valid] / ds[valid] if np.any(valid) else np.array([], dtype=float)
    min_speed = float(np.min(speeds)) if speeds.size else float("nan")
    max_speed = float(np.max(speeds)) if speeds.size else float("nan")
    sign_changes, x_monotone = projected_x_sign_change_count(x_plot)
    l1, l2 = TIMO.segment_lengths(result.mu)
    near_count, min_non = near_overlap_stats(x_plot, y_plot, l1 + l2)
    folded = bool(
        (isfinite(min_speed) and min_speed < CURVE_SPEED_TOL)
        or not x_monotone
        or near_count > 0
    )
    warnings: list[str] = []
    if isfinite(min_speed) and min_speed < CURVE_SPEED_TOL:
        warnings.append("min_curve_speed_close_to_zero")
    if not x_monotone:
        warnings.append("projected_x_nonmonotone")
    if near_count > 0:
        warnings.append("near_overlap_indicator")
    return {
        "epsilon": result.epsilon,
        "mu": result.mu,
        "beta_deg": result.beta_deg,
        "eta": result.eta,
        "sorted_index": result.sorted_index,
        "rod_id": int(rod_id),
        "scale_kind": scale_kind,
        "deformation_scale_fraction": deformation_scale_fraction,
        "deformation_scale": float(deformation_scale_value),
        "min_abs_dr_plot_ds": min_speed,
        "max_abs_dr_plot_ds": max_speed,
        "projected_x_sign_change_count": sign_changes,
        "projected_x_is_monotone": x_monotone,
        "near_overlap_pair_count": near_count,
        "min_nonneighbor_distance": min_non,
        "visually_folded_at_scale": folded,
        "warning": "; ".join(warnings) if warnings else "",
        "notes": "material order, separate rod line object; full centerline geometry",
    }


def regularity_rows(results: Sequence[ModeResult], fixed_fractions: Sequence[float]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        if result.model != MODEL_TIMO or result.case_role != "suspect":
            continue
        auto_scale = deformation_scale(result, transverse_only=False)
        for rod_id in (1, 2):
            rows.append(
                regularity_row(
                    result,
                    rod_id=rod_id,
                    scale_kind="autoscaled_full",
                    deformation_scale_fraction=float("nan"),
                    deformation_scale_value=auto_scale,
                )
            )
        for fraction in fixed_fractions:
            fixed_scale = fixed_scale_from_fraction(result.mu, float(fraction))
            for rod_id in (1, 2):
                rows.append(
                    regularity_row(
                        result,
                        rod_id=rod_id,
                        scale_kind="fixed_fraction",
                        deformation_scale_fraction=float(fraction),
                        deformation_scale_value=fixed_scale,
                    )
                )
    return rows


def regularity_fold_count(rows: Sequence[dict[str, object]], *, scale_kind: str, fraction: float | None = None) -> int:
    count = 0
    for row in rows:
        if str(row["scale_kind"]) != scale_kind:
            continue
        if fraction is not None:
            try:
                if abs(float(row["deformation_scale_fraction"]) - float(fraction)) > 1.0e-12:
                    continue
            except (TypeError, ValueError):
                continue
        if bool(row["visually_folded_at_scale"]):
            count += 1
    return count


def fixed_scale_plot_joint_gap(result: ModeResult, fraction: float) -> float:
    scale = fixed_scale_from_fraction(result.mu, float(fraction))
    xd1, yd1, xd2, yd2 = deformed_coordinates_with_scale(result, transverse_only=False, scale=scale)
    return float(np.linalg.norm(np.array([xd1[-1] - xd2[0], yd1[-1] - yd2[0]], dtype=float)))


def scale_check_rows(
    results: Sequence[ModeResult],
    order_rows: Sequence[dict[str, object]],
    *,
    fraction: float,
) -> list[dict[str, object]]:
    point_warnings: dict[tuple[float, float, int], list[str]] = {}
    for row in order_rows:
        warning = str(row.get("warning", ""))
        if not warning:
            continue
        key = result_scale_key_from_values(
            float(row["epsilon"]),
            float(row.get("mu", 0.0)),
            int(row["sorted_index"]),
        )
        point_warnings.setdefault(key, []).append(warning)

    rows: list[dict[str, object]] = []
    for result in results:
        if result.model != MODEL_TIMO or result.case_role != "suspect":
            continue
        fixed_scale = fixed_scale_from_fraction(result.mu, float(fraction))
        rod_rows = [
            regularity_row(
                result,
                rod_id=rod_id,
                scale_kind="fixed_fraction",
                deformation_scale_fraction=float(fraction),
                deformation_scale_value=fixed_scale,
            )
            for rod_id in (1, 2)
        ]
        min_drds = min((float(row["min_abs_dr_plot_ds"]) for row in rod_rows), default=float("nan"))
        max_drds = max((float(row["max_abs_dr_plot_ds"]) for row in rod_rows), default=float("nan"))
        fold_flag = any(bool(row["visually_folded_at_scale"]) for row in rod_rows)
        near_overlap_flag = any(int(row["near_overlap_pair_count"]) > 0 for row in rod_rows)
        plotted_gap = fixed_scale_plot_joint_gap(result, float(fraction))
        joint = result.joint_row or {}
        max_kin = float(joint.get("max_abs_kinematic_gap", float("nan")))
        max_force = float(joint.get("max_abs_force_gap", float("nan")))
        warnings = [str(row["warning"]) for row in rod_rows if str(row["warning"])]
        point_warning_values = point_warnings.get(result_scale_key(result), [])
        passed = bool(
            not fold_flag
            and not near_overlap_flag
            and plotted_gap <= JOINT_TOL
            and max_kin <= JOINT_TOL
            and max_force <= JOINT_TOL
            and not point_warning_values
        )
        note_parts: list[str] = []
        if warnings:
            note_parts.append("; ".join(warnings))
        if point_warning_values:
            note_parts.append("point_order_warning=" + "; ".join(point_warning_values))
        if passed:
            note_parts.append("passed; material point order; separate rod line objects; no artificial connector")
        rows.append(
            {
                "epsilon": result.epsilon,
                "mu": result.mu,
                "beta_deg": result.beta_deg,
                "eta": result.eta,
                "sorted_index": result.sorted_index,
                "Lambda": result.Lambda,
                "scale_fraction": float(fraction),
                "min_drds": min_drds,
                "max_drds": max_drds,
                "fold_flag": fold_flag,
                "near_overlap_flag": near_overlap_flag,
                "plotted_joint_gap": plotted_gap,
                "max_abs_kinematic_gap": max_kin,
                "max_abs_force_gap": max_force,
                "passed_scale_check": passed,
                "notes": "; ".join(note_parts) if note_parts else "scale check failed without classified warning",
            }
        )
    return rows


def recommended_fixed_fraction(rows: Sequence[dict[str, object]], fixed_fractions: Sequence[float]) -> float:
    for fraction in reversed(tuple(float(value) for value in fixed_fractions)):
        selected = [
            row
            for row in rows
            if str(row["scale_kind"]) == "fixed_fraction"
            and abs(float(row["deformation_scale_fraction"]) - float(fraction)) <= 1.0e-12
        ]
        if selected and all(
            int(row["near_overlap_pair_count"]) == 0 and not bool(row["visually_folded_at_scale"])
            for row in selected
        ):
            return float(fraction)
    return float(tuple(float(value) for value in fixed_fractions)[0])


def bad_x2_plus_l2_gap(result: ModeResult) -> float:
    if result.model != MODEL_TIMO:
        return float("nan")
    l1, l2 = TIMO.segment_lengths(result.mu)
    factors = TIMO.tau_factors(result.mu, result.eta)
    section2 = TIMO.section_from_epsilon_tau(result.epsilon, factors.tau2)
    basis2 = TIMO.timo_basis(result.Lambda, result.epsilon, section2)
    wrong = TIMO.evaluate_timo_rod_fields(
        result.Lambda,
        result.epsilon,
        section2,
        basis2,
        result.coeff[3:6],
        np.array([l2], dtype=float),
    )
    w1 = float(result.rod1.w[-1])
    u1 = float(result.rod1.u[-1])
    psi1 = float(np.asarray(result.rod1.psi, dtype=float)[-1])
    w2 = float(np.asarray(wrong["w"], dtype=float)[0])
    u2 = float(np.asarray(wrong["u"], dtype=float)[0])
    psi2 = float(np.asarray(wrong["psi"], dtype=float)[0])
    beta_rad = np.deg2rad(result.beta_deg)
    cb = float(np.cos(beta_rad))
    sb = float(np.sin(beta_rad))
    gaps = (w1 - cb * w2 + sb * u2, u1 - sb * w2 - cb * u2, psi1 - psi2)
    return max(abs(float(value)) for value in gaps)


def bad_transform_joint_gaps(result: ModeResult) -> tuple[float, float, float]:
    if result.model != MODEL_TIMO:
        return float("nan"), float("nan"), float("nan")
    beta_rad = np.deg2rad(result.beta_deg)
    cb = float(np.cos(beta_rad))
    sb = float(np.sin(beta_rad))
    sign = float(result.sign_factor)
    u1 = sign * float(result.rod1.u[-1])
    w1 = sign * float(result.rod1.w[-1])
    u2 = sign * float(result.rod2.u[-1])
    w2 = sign * float(result.rod2.w[-1])
    scale = deformation_scale(result, transverse_only=False)
    good_dx, good_dy = DISPLAY.rod2_local_displacement_to_display(
        np.array([u2], dtype=float),
        np.array([w2], dtype=float),
        beta_deg=result.beta_deg,
    )
    rod1_dx, rod1_dy = DISPLAY.rod1_local_displacement_to_display(
        np.array([u1], dtype=float),
        np.array([w1], dtype=float),
    )
    good = np.array([good_dx[0], good_dy[0]], dtype=float)
    standard = np.array([cb * u2 + sb * w2, -sb * u2 + cb * w2], dtype=float)
    local = np.array([u2, w2], dtype=float)
    rod1 = np.array([rod1_dx[0], rod1_dy[0]], dtype=float)
    good_gap = float(np.linalg.norm(scale * (rod1 - good)))
    standard_gap = float(np.linalg.norm(scale * (rod1 - standard)))
    local_gap = float(np.linalg.norm(scale * (rod1 - local)))
    return good_gap, standard_gap, local_gap


def bad_external_to_joint_connector_length(result: ModeResult) -> float:
    if result.model != MODEL_TIMO:
        return float("nan")
    _l1, l2 = TIMO.segment_lengths(result.mu)
    return float(l2)


def summary_row(
    result: ModeResult,
    figure_files: dict[tuple[float, int], dict[str, Path]],
    scale_check_by_key: dict[tuple[float, float, int], dict[str, object]],
    *,
    candidate_fraction: float,
    fallback_fraction: float,
) -> dict[str, object]:
    max_u, max_w = max_abs_u_w(result)
    full_scale = deformation_scale(result, transverse_only=False)
    trans_scale = deformation_scale(result, transverse_only=True)
    if result.model == MODEL_TIMO:
        coeff_order = "(A1,B1,P1,A2,B2,P2)"
        rod1_coefficients = "(A1,B1,P1)"
        rod2_coefficients = "(A2,B2,P2)"
        transform = "dX2=c*u2+s*w2; dY2=s*u2-c*w2; rod1 dY1=-w1"
        point_order = "rod1 clamp-to-joint; rod2 joint-to-clamp; rods plotted separately"
        rod2_first = "joint"
        rod2_last = "clamp"
        bad_plus = bad_x2_plus_l2_gap(result)
        _good_gap, bad_standard, bad_local = bad_transform_joint_gaps(result)
        bad_connector = bad_external_to_joint_connector_length(result)
    else:
        coeff_order = "(A1,B1,A2,B2,P1,P2)"
        rod1_coefficients = "(A1,B1,P1) from positions 1,2,5"
        rod2_coefficients = "(A2,B2,P2) from positions 3,4,6"
        transform = "shared display basis: dX2=c*u2+s*w2; dY2=s*u2-c*w2"
        point_order = "rod1 clamp-to-joint; rod2 joint-to-clamp; rods plotted separately"
        rod2_first = "joint"
        rod2_last = "clamp"
        bad_plus = float("nan")
        bad_standard = float("nan")
        bad_local = float("nan")
        bad_connector = float("nan")
    key = (round(result.epsilon, 12), result.sorted_index)
    files = figure_files.get(key, {})
    full_shape = files.get("full_shape") or files.get("shape_grid")
    transverse_shape = files.get("transverse_shape") or files.get("shape_grid")
    component_file = files.get("component_grid")
    scale_row = scale_check_by_key.get(result_scale_key(result)) if result.model == MODEL_TIMO else None
    if scale_row is None:
        final_scale_fraction: object = ""
        scale_check_passed: object = ""
        fallback_scale_used: object = ""
        final_figure_file: object = ""
    else:
        scale_check_passed = bool(scale_row["passed_scale_check"])
        fallback_scale_used = not scale_check_passed
        final_scale_fraction = float(candidate_fraction) if scale_check_passed else float(fallback_fraction)
        final_file = files.get("fixed_scale_shape") if scale_check_passed else files.get("fallback_shape")
        final_figure_file = rel(final_file) if final_file is not None else ""
    return {
        "epsilon": result.epsilon,
        "mu": result.mu,
        "beta_deg": result.beta_deg,
        "eta": result.eta,
        "sorted_index": result.sorted_index,
        "model": result.model,
        "Lambda": result.Lambda,
        "case_role": result.case_role,
        "full_shape_file": rel(full_shape) if full_shape is not None else "",
        "transverse_shape_file": rel(transverse_shape) if transverse_shape is not None else "",
        "component_file": rel(component_file) if result.model == MODEL_TIMO and component_file is not None else "",
        "final_scale_fraction": final_scale_fraction,
        "scale_check_passed": scale_check_passed,
        "fallback_scale_used": fallback_scale_used,
        "final_figure_file": final_figure_file,
        "coefficient_order": coeff_order,
        "rod1_coefficients": rod1_coefficients,
        "rod2_coefficients": rod2_coefficients,
        "rod1_coordinate": "x1 from 0 to +l1",
        "rod2_coordinate": "x2 from 0 to -l2; reversed for geometry plots",
        "rod2_global_transform": transform,
        "plot_object_policy": "separate line objects for rod1 and rod2",
        "point_order": point_order,
        "rod2_first_point_role": rod2_first,
        "rod2_last_point_role": rod2_last,
        "full_deformation_scale": full_scale,
        "transverse_deformation_scale": trans_scale,
        "max_abs_u": max_u,
        "max_abs_w": max_w,
        "axial_energy_fraction": result.energy["axial_energy_fraction"],
        "bending_energy_fraction": result.energy["bending_energy_fraction"],
        "shear_energy_fraction": result.energy["shear_energy_fraction"],
        "classification": result.energy["classification"],
        "bad_x2_plus_l2_max_kinematic_gap": bad_plus,
        "bad_standard_rotation_plot_joint_gap": bad_standard,
        "bad_local_as_global_plot_joint_gap": bad_local,
        "bad_external_to_joint_connector_length": bad_connector,
        "notes": (
            f"root_source={result.root_source}; sign_factor={result.sign_factor:+.0f}; "
            f"warnings={'; '.join(result.warnings) if result.warnings else 'none'}"
        ),
    }


def read_old_energy_rows() -> list[dict[str, str]]:
    old_path = REPO_ROOT / "results" / "eb_vs_timoshenko_longitudinal_suspect_modes" / "suspect_mode_energy_fractions.csv"
    if not old_path.exists():
        return []
    with old_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def old_energy_comparison(results: Sequence[ModeResult]) -> dict[str, object]:
    old_rows = read_old_energy_rows()
    if not old_rows:
        return {"status": "old_energy_csv_not_found", "max_abs_delta": float("nan"), "matched_rows": 0}
    old_by_key: dict[tuple[float, float, str, int], dict[str, str]] = {}
    for row in old_rows:
        try:
            key = (
                round(float(row["epsilon"]), 12),
                round(float(row["mu"]), 12),
                str(row["model"]),
                int(float(row["sorted_index"])),
            )
        except (KeyError, ValueError):
            continue
        old_by_key[key] = row
    deltas: list[float] = []
    for result in results:
        if result.case_role != "suspect":
            continue
        key = (round(result.epsilon, 12), round(result.mu, 12), result.model, result.sorted_index)
        old = old_by_key.get(key)
        if old is None:
            continue
        for field in ("axial_energy_fraction", "bending_energy_fraction", "shear_energy_fraction"):
            try:
                deltas.append(abs(float(old[field]) - float(result.energy[field])))
            except (KeyError, ValueError):
                continue
    if not deltas:
        return {"status": "old_energy_rows_not_matched", "max_abs_delta": float("nan"), "matched_rows": 0}
    return {
        "status": "compared_to_existing_longitudinal_suspect_energy_csv",
        "max_abs_delta": max(deltas),
        "matched_rows": len(deltas) // 3,
    }


def write_report(
    path: Path,
    *,
    results: Sequence[ModeResult],
    joint_rows: Sequence[dict[str, object]],
    point_order_rows: Sequence[dict[str, object]],
    regularity_rows: Sequence[dict[str, object]],
    scale_check_rows: Sequence[dict[str, object]],
    summary_csv: Path,
    joint_csv: Path,
    point_order_csv: Path,
    regularity_csv: Path,
    scale_check_csv: Path,
    figure_paths: Sequence[Path],
    energy_compare: dict[str, object],
    deformation_scale_fraction: float,
    fixed_fractions: Sequence[float],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    timo_results = [result for result in results if result.model == MODEL_TIMO]
    max_kin = max((float(row["max_abs_kinematic_gap"]) for row in joint_rows), default=float("nan"))
    max_force = max((float(row["max_abs_force_gap"]) for row in joint_rows), default=float("nan"))
    max_plot = max(
        (
            max(float(row["undeformed_plot_joint_gap"]), float(row["deformed_plot_joint_gap"]))
            for row in joint_rows
        ),
        default=float("nan"),
    )
    max_bad_standard = max(
        (bad_transform_joint_gaps(result)[1] for result in timo_results),
        default=float("nan"),
    )
    max_bad_local = max(
        (bad_transform_joint_gaps(result)[2] for result in timo_results),
        default=float("nan"),
    )
    max_bad_x2 = max(
        (bad_x2_plus_l2_gap(result) for result in timo_results),
        default=float("nan"),
    )
    max_bad_connector = max(
        (bad_external_to_joint_connector_length(result) for result in timo_results),
        default=float("nan"),
    )
    suspect_timo = [result for result in timo_results if result.case_role == "suspect"]
    longitudinal_count = sum(1 for result in suspect_timo if result.energy["classification"] == "longitudinal_dominated")
    mixed_count = sum(1 for result in suspect_timo if result.energy["classification"] == "mixed")
    bending_count = sum(1 for result in suspect_timo if result.energy["classification"] == "bending_dominated")
    point_order_warning_count = sum(1 for row in point_order_rows if str(row.get("warning", "")))
    autoscale_fold_count = regularity_fold_count(regularity_rows, scale_kind="autoscaled_full")
    fixed_fold_counts = {
        float(fraction): regularity_fold_count(regularity_rows, scale_kind="fixed_fraction", fraction=float(fraction))
        for fraction in fixed_fractions
    }
    fixed_near_overlap_counts = {
        float(fraction): sum(
            int(row["near_overlap_pair_count"])
            for row in regularity_rows
            if str(row["scale_kind"]) == "fixed_fraction"
            and abs(float(row["deformation_scale_fraction"]) - float(fraction)) <= 1.0e-12
        )
        for fraction in fixed_fractions
    }
    autoscale_near_overlap_count = sum(
        int(row["near_overlap_pair_count"])
        for row in regularity_rows
        if str(row["scale_kind"]) == "autoscaled_full"
    )
    recommended_fraction = recommended_fixed_fraction(regularity_rows, fixed_fractions)
    min_fixed_speed = min(
        (
            float(row["min_abs_dr_plot_ds"])
            for row in regularity_rows
            if str(row["scale_kind"]) == "fixed_fraction" and isfinite(float(row["min_abs_dr_plot_ds"]))
        ),
        default=float("nan"),
    )
    min_auto_speed = min(
        (
            float(row["min_abs_dr_plot_ds"])
            for row in regularity_rows
            if str(row["scale_kind"]) == "autoscaled_full" and isfinite(float(row["min_abs_dr_plot_ds"]))
        ),
        default=float("nan"),
    )
    fixed_fold_text = ", ".join(f"{fraction:.2f}: {count}" for fraction, count in fixed_fold_counts.items())
    fixed_overlap_text = ", ".join(f"{fraction:.2f}: {count}" for fraction, count in fixed_near_overlap_counts.items())
    scale_total = len(scale_check_rows)
    scale_pass_count = sum(1 for row in scale_check_rows if bool(row["passed_scale_check"]))
    scale_fail_count = scale_total - scale_pass_count
    scale_fold_flags = sum(1 for row in scale_check_rows if bool(row["fold_flag"]))
    scale_near_flags = sum(1 for row in scale_check_rows if bool(row["near_overlap_flag"]))
    scale_max_plot_gap = max((float(row["plotted_joint_gap"]) for row in scale_check_rows), default=float("nan"))
    failed_scale_modes = [
        (
            f"eps={float(row['epsilon']):g}, mu={float(row['mu']):g}, "
            f"sorted={int(row['sorted_index'])}"
        )
        for row in scale_check_rows
        if not bool(row["passed_scale_check"])
    ]
    if scale_fail_count == 0:
        scale_recommendation = f"{float(deformation_scale_fraction):.2f}"
    else:
        scale_recommendation = (
            f"{FALLBACK_DEFORMATION_SCALE_FRACTION:.2f} for failed shapes "
            f"({'; '.join(failed_scale_modes)}) and {float(deformation_scale_fraction):.2f} for the rest"
        )

    lines = [
        "# Timoshenko Shape Construction Audit",
        "",
        "## Scope",
        "",
        "This diagnostic uses sorted analytic frequencies only. It does not use descendant tracking, FEM, 3D FEM, Gmsh, CalculiX, article workspaces, main.tex, article figures, old determinant code, old solvers, or baseline-result edits.",
        "",
        "No analytic formulas, determinant entries, root solvers, or Timoshenko shear coefficient k' were changed. The script only reconstructs fields from existing helpers and audits plotting geometry.",
        "",
        "## Shape Construction Procedure",
        "",
        "1. Timoshenko coefficient vector order is `(A1,B1,P1,A2,B2,P2)`. The split is `rod1=(A1,B1,P1)` and `rod2=(A2,B2,P2)`.",
        "2. Rod 1 is sampled on `x1` from `0` to `+l1`. Rod 2 is sampled on `x2` from `0` to `-l2`, matching the determinant endpoint `x2=-l2` at the joint.",
        "3. For geometry plotting, rod 2 samples are reversed so the plotted order is joint-to-clamp. This avoids a connector from the joint to the far clamp.",
        "4. Determinant-frame components are not used directly as Cartesian coordinates. The display bases are `t1=(1,0)`, `n1=(0,-1)`, `t2=(c,s)`, and `n2=(s,-c)`, so `d1=(u1,-w1)` and `d2=(c*u2+s*w2,s*u2-c*w2)`.",
        "5. Rods are plotted as separate line objects. No single concatenated polyline is used for the two rods.",
        "6. Deformation scale is global per mode and per visualization type. Full centerline scale uses max `sqrt(u^2+w^2)` over both rods. Transverse-only scale uses max `|w|` over both rods. There is no separate rod normalization.",
        "7. The fixed-scale sensitivity plots do not use per-mode amplitude normalization. They use `scale = deformation_scale_fraction * (l1+l2)` on the reconstructed displacement fields.",
        "",
        "## Joint Continuity",
        "",
        f"- Joint continuity CSV: `{rel(joint_csv)}`",
        f"- Point-order CSV: `{rel(point_order_csv)}`",
        f"- Curve-regularity CSV: `{rel(regularity_csv)}`",
        f"- Scale `{float(deformation_scale_fraction):.2f}` regularity CSV: `{rel(scale_check_csv)}`",
        f"- Shape summary CSV: `{rel(summary_csv)}`",
        f"- Maximum mathematical kinematic gap: `{max_kin:.6e}`",
        f"- Maximum mathematical force/moment gap: `{max_force:.6e}`",
        f"- Maximum plotting joint gap: `{max_plot:.6e}`",
        "",
        "All corrected Timoshenko plots use the same signed rod-2 endpoint as the determinant (`x2=-l2`). If the mathematical gaps are small and the plotting gaps are small, the reconstructed shape is compatible at the joint.",
        "",
        "## Plotting Bug Diagnosis",
        "",
        f"- Corrected reflected-display plotting gap stays at `{max_plot:.6e}`.",
        f"- Using the old determinant transform directly as Cartesian display components gives a maximum artificial plot joint gap of `{max_bad_standard:.6e}` after the required rod-1 reflection.",
        f"- Plotting local `u2,w2` directly as global coordinates gives a maximum artificial plot joint gap of `{max_bad_local:.6e}`.",
        f"- Evaluating the rod-2 joint at `+l2` instead of `-l2` gives a maximum kinematic gap of `{max_bad_x2:.6e}`.",
        f"- Concatenating rod 1 with rod 2 in external-to-joint order would add an artificial connector up to `{max_bad_connector:.6e}` long.",
        "",
        "Conclusion: the old global plots mixed determinant-frame components with Cartesian display coordinates. The corrected plots use `x2=-l2`, the reflected display bases, one global scale, and separate line objects.",
        "",
        "## Invalidated Global Figures",
        "",
        "Timoshenko global centerline and displacement-vector figures generated before this correction used the inconsistent rod-2 display frame and must not be used for physical interpretation. Their local fields, frequencies, null-vector checks, joint equations, and energy fractions remain valid. Corrected global figures use `t2=(cos(beta),sin(beta))` and `n2=(sin(beta),-cos(beta))`; positive determinant `w1` maps to negative display Y.",
        "",
        "Strong axial displacement can still make the full deformed centerline visually unintuitive in longitudinal or mixed modes. That is a visualization limitation, not a joint compatibility failure.",
        "",
        "## Point Order and Scale Sensitivity",
        "",
        f"- Point-order warnings in the corrected plots: `{point_order_warning_count}`.",
        "- Corrected plots keep material-coordinate order for each rod and do not sort by global x.",
        f"- Autoscaled full centerlines: folded/near-overlap regularity rows `{autoscale_fold_count}`; near-overlap pair indicators `{autoscale_near_overlap_count}`; min `|dr/ds|` `{min_auto_speed:.6e}`.",
        f"- Fixed-scale folded/near-overlap rows by deformation fraction: {fixed_fold_text}.",
        f"- Fixed-scale near-overlap pair indicators by deformation fraction: {fixed_overlap_text}.",
        f"- Minimum fixed-scale `|dr/ds|`: `{min_fixed_speed:.6e}`.",
        f"- Common all-shape safe full-centerline diagnostic scale: `{recommended_fraction:.2f}` of total length.",
        "",
        "The fixed-scale figures are the preferred way to inspect full centerlines. The component figures are the primary diagnostic for modal character because they display `u`, `w`, `psi`, and `gamma=w'-psi` directly on signed local coordinates and mark the joint.",
        "",
        "## Choice of deformation scale",
        "",
        "- Scale conclusions from the old inconsistent display frame are invalid. Scale `0.05` and the requested candidate are recomputed with the reflected display geometry in this run.",
        f"- Scale `{float(deformation_scale_fraction):.2f}` was tested as the intermediate candidate.",
        f"- Shapes checked at scale `{float(deformation_scale_fraction):.2f}`: `{scale_total}`.",
        f"- Passed/failed at scale `{float(deformation_scale_fraction):.2f}`: `{scale_pass_count}` / `{scale_fail_count}`.",
        f"- Fold flags at scale `{float(deformation_scale_fraction):.2f}`: `{scale_fold_flags}`.",
        f"- Near-overlap flags at scale `{float(deformation_scale_fraction):.2f}`: `{scale_near_flags}`.",
        f"- Maximum plotted joint gap at scale `{float(deformation_scale_fraction):.2f}`: `{scale_max_plot_gap:.6e}`.",
        f"- Recommended deformation scale: `{scale_recommendation}`.",
        "",
        "The deformation scale affects only visualization of the full centerline. It does not affect sorted frequencies, energy fractions, mode classification, analytic formulas, determinant entries, root solvers, or the Timoshenko shear coefficient k'.",
        "",
        "## Energy Fractions",
        "",
        f"- Energy comparison status: `{energy_compare['status']}`",
        f"- Matched old rows: `{energy_compare['matched_rows']}`",
        f"- Max absolute old/new energy-fraction delta: `{float(energy_compare['max_abs_delta']):.6e}`",
        f"- Suspect Timoshenko classification counts: longitudinal={longitudinal_count}, mixed={mixed_count}, bending={bending_count}.",
        "",
        "The energy-based conclusion is separate from visual plotting. The energy fractions are computed from the same existing stiffness-energy terms as before; only the plotting geometry is audited here.",
        "",
        "## Answers",
        "",
        "1. Are the Timoshenko fields continuous/compatible at the joint? Yes, within the mathematical and plotting gaps reported above.",
        f"2. Are apparent jumps caused by plotting order? In the corrected plots, no point-order bug was found (`{point_order_warning_count}` warnings). A wrong concatenation order would still create an artificial connector, so rods are plotted separately.",
        "3. Are apparent jumps caused by large deformation scale? The regularity audit shows that autoscaled full centerlines can become visually folded/near-overlapping, especially when axial displacement is amplified.",
        f"4. Does a small fixed scale remove the visual jumps? Use the final per-shape scale from the Choice of deformation scale section; the common all-shape safe scale is `{recommended_fraction:.2f}` of total length.",
        "5. For longitudinal/mixed modes, should component plots replace full deformed centerlines? Yes. Use component plots and energy fractions as primary evidence; use fixed-scale full centerlines only as secondary geometry checks.",
        "6. Does the longitudinal/mixed conclusion change? No. The conclusion follows the energy fractions, which are unchanged except for roundoff-level recomputation differences if an older CSV is present.",
        "",
        "## Generated Figures",
        "",
    ]
    for figure in figure_paths:
        lines.append(f"- `{rel(figure)}`")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def default_suspects() -> tuple[SuspectSpec, ...]:
    return (
        SuspectSpec(epsilon=0.03, sorted_index=5, mu_values=(0.0, 0.35, 0.54, 0.7)),
        SuspectSpec(epsilon=0.05, sorted_index=4, mu_values=(0.0, 0.35, 0.7)),
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Focused diagnostic audit of Timoshenko shape construction and plotting geometry.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--beta-deg", type=float, default=DEFAULT_BETA_DEG)
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--n-roots", type=int, default=DEFAULT_N_ROOTS)
    parser.add_argument("--n-points", type=int, default=DEFAULT_N_POINTS)
    parser.add_argument("--deformation-scale-fraction", type=float, default=DEFAULT_DEFORMATION_SCALE_FRACTION)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    if bool(args.smoke):
        args.output_dir = SMOKE_OUTPUT_DIR
        args.n_points = SMOKE_N_POINTS
        args.n_roots = max(args.n_roots, 4)
    args.output_dir = repo_path(Path(args.output_dir))
    if abs(float(args.eta)) > 1.0e-14:
        raise ValueError("this focused audit is restricted to eta=0, as requested")
    if int(args.n_points) < 51:
        raise ValueError("--n-points must be at least 51")
    if int(args.n_roots) < 6 and not bool(args.smoke):
        raise ValueError("--n-roots must be at least 6 for the default suspect set")
    if float(args.deformation_scale_fraction) <= 0.0:
        raise ValueError("--deformation-scale-fraction must be positive")
    return args


def build_results(args: argparse.Namespace) -> tuple[list[ModeResult], tuple[SuspectSpec, ...]]:
    if bool(args.smoke):
        suspects = (SuspectSpec(epsilon=0.05, sorted_index=4, mu_values=(0.35,)),)
        control_sorted = (1,)
    else:
        suspects = default_suspects()
        control_sorted = CONTROL_SORTED

    provider = RootProvider(float(args.beta_deg), float(args.eta), int(args.n_roots))
    results: list[ModeResult] = []

    control_pair = provider.pair(CONTROL_EPSILON, CONTROL_MU, max(control_sorted))
    eb_roots, timo_roots, root_warnings = provider.roots(CONTROL_EPSILON, CONTROL_MU)
    for sorted_index in control_sorted:
        Lambda = float(timo_roots[int(sorted_index) - 1])
        results.append(
            timo_mode_result(
                epsilon=CONTROL_EPSILON,
                beta_deg=float(args.beta_deg),
                eta=float(args.eta),
                mu=CONTROL_MU,
                sorted_index=int(sorted_index),
                Lambda=Lambda,
                n_points=int(args.n_points),
                case_role="control",
                root_source=control_pair.source,
                root_warnings=root_warnings,
            )
        )

    for spec in suspects:
        for mu in spec.mu_values:
            pair = provider.pair(spec.epsilon, float(mu), spec.sorted_index)
            results.append(
                eb_mode_result(
                    epsilon=spec.epsilon,
                    beta_deg=float(args.beta_deg),
                    eta=float(args.eta),
                    mu=float(mu),
                    sorted_index=spec.sorted_index,
                    Lambda=pair.eb,
                    n_points=int(args.n_points),
                    case_role="suspect",
                    root_source=pair.source,
                    root_warnings=pair.warnings,
                )
            )
            results.append(
                timo_mode_result(
                    epsilon=spec.epsilon,
                    beta_deg=float(args.beta_deg),
                    eta=float(args.eta),
                    mu=float(mu),
                    sorted_index=spec.sorted_index,
                    Lambda=pair.timo,
                    n_points=int(args.n_points),
                    case_role="suspect",
                    root_source=pair.source,
                    root_warnings=pair.warnings,
                )
            )
    return results, suspects


def write_figures(
    output_dir: Path,
    results: Sequence[ModeResult],
    suspects: Sequence[SuspectSpec],
    *,
    deformation_scale_fraction: float,
    fixed_fractions: Sequence[float],
    scale_rows: Sequence[dict[str, object]],
) -> tuple[list[Path], dict[tuple[float, int], dict[str, Path]]]:
    figure_paths: list[Path] = []
    figure_files: dict[tuple[float, int], dict[str, Path]] = {}
    control_dir = output_dir / "control_modes"
    for result in results:
        if result.case_role != "control" or result.model != MODEL_TIMO:
            continue
        stem = f"timoshenko_control_sorted{result.sorted_index}_mu{mu_token(result.mu)}_eps{eps_token(result.epsilon)}"
        full = plot_single_centerline(result, control_dir / f"{stem}.png", transverse_only=False)
        transverse = plot_single_centerline(result, control_dir / f"{stem}_transverse_only.png", transverse_only=True)
        components = plot_timo_components(result, control_dir / f"{stem}_components.png")
        figure_paths.extend([full, transverse, components])
        figure_files[(round(result.epsilon, 12), result.sorted_index)] = {
            "full_shape": full,
            "transverse_shape": transverse,
            "component_grid": components,
        }

    for spec in suspects:
        component_grid = plot_suspect_component_grid(output_dir, spec, results)
        shape_grid = plot_suspect_shape_grid(output_dir, spec, results)
        fixed_shape_grid = plot_suspect_fixed_scale_shape_grid(
            output_dir,
            spec,
            results,
            fraction=float(deformation_scale_fraction),
        )
        scale_sensitivity = plot_scale_sensitivity_grid(output_dir, spec, results, fixed_fractions)
        figure_paths.extend([component_grid, shape_grid, fixed_shape_grid, scale_sensitivity])
        files = {
            "shape_grid": shape_grid,
            "component_grid": component_grid,
            "fixed_scale_shape": fixed_shape_grid,
            "scale_sensitivity": scale_sensitivity,
        }
        failed_mus = [
            float(row["mu"])
            for row in scale_rows
            if abs(float(row["epsilon"]) - float(spec.epsilon)) <= 1.0e-12
            and int(row["sorted_index"]) == int(spec.sorted_index)
            and not bool(row["passed_scale_check"])
        ]
        if failed_mus:
            fallback_shape = plot_suspect_fixed_scale_shape_grid(
                output_dir,
                spec,
                results,
                fraction=FALLBACK_DEFORMATION_SCALE_FRACTION,
                mu_values=failed_mus,
                suffix="_fallback",
            )
            figure_paths.append(fallback_shape)
            files["fallback_shape"] = fallback_shape
        if abs(float(spec.epsilon) - 0.05) <= 1.0e-12 and int(spec.sorted_index) == 4:
            comparison = plot_scale_comparison_grid(
                output_dir,
                spec,
                results,
                fractions=(FALLBACK_DEFORMATION_SCALE_FRACTION, float(deformation_scale_fraction)),
            )
            figure_paths.append(comparison)
            files["scale_comparison"] = comparison
        figure_files[(round(spec.epsilon, 12), spec.sorted_index)] = files
    return figure_paths, figure_files


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    fixed_fractions = fixed_scale_fractions(float(args.deformation_scale_fraction))
    results, suspects = build_results(args)
    joint_rows = [dict(result.joint_row) for result in results if result.model == MODEL_TIMO and result.joint_row is not None]
    order_rows = point_order_rows(results)
    curve_rows = regularity_rows(results, fixed_fractions)
    scale_rows = scale_check_rows(results, order_rows, fraction=float(args.deformation_scale_fraction))
    figure_paths, figure_files = write_figures(
        args.output_dir,
        results,
        suspects,
        deformation_scale_fraction=float(args.deformation_scale_fraction),
        fixed_fractions=fixed_fractions,
        scale_rows=scale_rows,
    )
    scale_check_by_key = {
        result_scale_key_from_values(float(row["epsilon"]), float(row["mu"]), int(row["sorted_index"])): row
        for row in scale_rows
    }
    summary_rows = [
        summary_row(
            result,
            figure_files,
            scale_check_by_key,
            candidate_fraction=float(args.deformation_scale_fraction),
            fallback_fraction=FALLBACK_DEFORMATION_SCALE_FRACTION,
        )
        for result in results
    ]
    summary_csv = write_csv(args.output_dir / "timoshenko_shape_construction_summary.csv", summary_rows, SUMMARY_FIELDS)
    joint_csv = write_csv(args.output_dir / "timoshenko_joint_continuity_audit_corrected.csv", joint_rows, JOINT_FIELDS)
    point_order_csv = write_csv(args.output_dir / "timoshenko_shape_point_order_audit.csv", order_rows, POINT_ORDER_FIELDS)
    regularity_csv = write_csv(args.output_dir / "timoshenko_plotted_curve_regularity_audit.csv", curve_rows, REGULARITY_FIELDS)
    scale_check_csv = write_csv(
        args.output_dir / f"timoshenko_plotted_curve_regularity_scale{scale_token(float(args.deformation_scale_fraction))}.csv",
        scale_rows,
        SCALE_CHECK_FIELDS,
    )
    energy_compare = old_energy_comparison(results)
    report = write_report(
        args.output_dir / "timoshenko_shape_construction_audit_report.md",
        results=results,
        joint_rows=joint_rows,
        point_order_rows=order_rows,
        regularity_rows=curve_rows,
        scale_check_rows=scale_rows,
        summary_csv=summary_csv,
        joint_csv=joint_csv,
        point_order_csv=point_order_csv,
        regularity_csv=regularity_csv,
        scale_check_csv=scale_check_csv,
        figure_paths=figure_paths,
        energy_compare=energy_compare,
        deformation_scale_fraction=float(args.deformation_scale_fraction),
        fixed_fractions=fixed_fractions,
    )

    max_kin = max((float(row["max_abs_kinematic_gap"]) for row in joint_rows), default=float("nan"))
    max_force = max((float(row["max_abs_force_gap"]) for row in joint_rows), default=float("nan"))
    max_plot = max(
        (
            max(float(row["undeformed_plot_joint_gap"]), float(row["deformed_plot_joint_gap"]))
            for row in joint_rows
        ),
        default=float("nan"),
    )
    print("generated Timoshenko shape construction audit outputs:")
    for path in [report, joint_csv, point_order_csv, regularity_csv, scale_check_csv, summary_csv, *figure_paths]:
        print(f"  {rel(path)}")
    print(f"max mathematical kinematic gap: {max_kin:.6e}")
    print(f"max mathematical force gap: {max_force:.6e}")
    print(f"max plotting joint gap: {max_plot:.6e}")
    print(f"point-order warnings: {sum(1 for row in order_rows if str(row.get('warning', '')))}")
    print(f"autoscaled folded regularity rows: {regularity_fold_count(curve_rows, scale_kind='autoscaled_full')}")
    for fraction in fixed_fractions:
        print(
            "fixed scale "
            f"{float(fraction):.2f} folded regularity rows: "
            f"{regularity_fold_count(curve_rows, scale_kind='fixed_fraction', fraction=float(fraction))}"
        )
    scale_pass_count = sum(1 for row in scale_rows if bool(row["passed_scale_check"]))
    scale_fold_flags = sum(1 for row in scale_rows if bool(row["fold_flag"]))
    scale_near_flags = sum(1 for row in scale_rows if bool(row["near_overlap_flag"]))
    scale_max_plot_gap = max((float(row["plotted_joint_gap"]) for row in scale_rows), default=float("nan"))
    print(
        f"scale {float(args.deformation_scale_fraction):.2f} shape checks: "
        f"{scale_pass_count}/{len(scale_rows)} passed"
    )
    print(
        f"scale {float(args.deformation_scale_fraction):.2f} fold flags: "
        f"{scale_fold_flags}, near-overlap flags: {scale_near_flags}"
    )
    print(f"scale {float(args.deformation_scale_fraction):.2f} max plotted joint gap: {scale_max_plot_gap:.6e}")
    print(f"energy comparison: {energy_compare['status']}, max delta={float(energy_compare['max_abs_delta']):.6e}")
    return {
        "report": report,
        "joint_csv": joint_csv,
        "point_order_csv": point_order_csv,
        "regularity_csv": regularity_csv,
        "scale_check_csv": scale_check_csv,
        "summary_csv": summary_csv,
        "figure_paths": figure_paths,
        "joint_rows": joint_rows,
        "point_order_rows": order_rows,
        "regularity_rows": curve_rows,
        "scale_check_rows": scale_rows,
        "summary_rows": summary_rows,
        "energy_compare": energy_compare,
    }


if __name__ == "__main__":
    main()
