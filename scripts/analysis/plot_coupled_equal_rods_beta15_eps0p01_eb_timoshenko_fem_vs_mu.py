from __future__ import annotations

import csv
from dataclasses import dataclass
import importlib.util
import math
from pathlib import Path
import sys
from types import ModuleType
from typing import Callable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.optimize import brentq, linear_sum_assignment


# =========================
# Diagnostic run parameters
# =========================
BETA_DEG = 15.0
EPSILON = 0.01
ETA = 0.0
MU_MIN = 0.0
MU_MAX = 0.9
ANALYTIC_MU_VALUES = np.round(np.linspace(MU_MIN, MU_MAX, 91), 10)
FEM_MU_VALUES = [0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9]
N_BRANCHES = 6
N_SOLVE_ANALYTIC = 14
N_SOLID_MODES = 30
FEM_MESH_SIZE_FACTOR = 1.0
USE_PLANAR_POINT_JOINT_FEM = True


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import assemble_clamped_coupled_matrix  # noqa: E402
from scripts.lib.analytic_branch_tracking import (  # noqa: E402
    branch_id_from_base_sorted_index,
    track_mu_sweep,
)


ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"
RESULTS_DIR = REPO_ROOT / "results"
OUTPUT_STEM = "coupled_equal_thickness_beta15_eps0p01_eb_timoshenko_3d_fem_vs_mu_corrected"
OUTPUT_PNG = RESULTS_DIR / f"{OUTPUT_STEM}.png"
OUTPUT_CSV = RESULTS_DIR / f"{OUTPUT_STEM}.csv"
OUTPUT_REPORT = RESULTS_DIR / f"{OUTPUT_STEM}_report.md"
FEM_OUTPUT_ROOT = RESULTS_DIR / "solid_fem_coupled_equal_rods_beta15_eps0p01_planar_mu"


E = 1.0
RHO = 1.0
NU = 0.3
L_SEGMENT = 1.0
THIN_ROD_LIMIT = 0.1
ROOT_SCAN_START = 0.2
ROOT_SCAN_STEP = 0.01
ROOT_DEDUP_TOL = 1.0e-7
N_SLICES_PER_ROD = 80
MAC_STRONG_THRESHOLD = 0.8
MAC_MODERATE_THRESHOLD = 0.5


@dataclass(frozen=True)
class Section:
    radius: float
    area: float
    inertia: float
    shear_modulus: float
    kappa: float

    @property
    def shear_stiffness(self) -> float:
        return self.kappa * self.shear_modulus * self.area

    @property
    def bending_stiffness(self) -> float:
        return E * self.inertia

    @property
    def mass_per_length(self) -> float:
        return RHO * self.area

    @property
    def rotary_inertia_per_length(self) -> float:
        return RHO * self.inertia


@dataclass(frozen=True)
class TimoshenkoBasis:
    a: float
    b: float
    h: float
    q: float
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class RodFrame:
    outer: tuple[float, float, float]
    joint: tuple[float, float, float]
    tangent: tuple[float, float, float]
    in_plane_normal: tuple[float, float, float]
    length: float


@dataclass(frozen=True)
class Geometry:
    joint: tuple[float, float, float]
    rod1: RodFrame
    rod2: RodFrame
    beta_rad: float
    mu: float


@dataclass(frozen=True)
class TimoTrackedPoint:
    branch_id: str
    base_sorted_index: int
    beta_deg: float
    mu: float
    current_sorted_index: int
    Lambda: float
    mac_to_previous: float
    vector: np.ndarray


@dataclass(frozen=True)
class SolidShape:
    mu: float
    solid_mode: int
    omega: float
    Lambda: float
    vector: np.ndarray
    filled_slice_count: int
    empty_slice_count: int
    out_of_plane_fraction: float


@dataclass(frozen=True)
class FemCaseResult:
    mu: float
    case_dir: Path
    success: bool
    message: str
    gmsh_messages: tuple[str, ...]
    mesh_quality_status: str
    negative_jacobian_warning: bool
    nodes: int
    solid_elements: int
    l1: float
    l2: float
    radius: float
    diameter_to_rod1_length: float
    diameter_to_rod2_length: float
    mesh_rod1_length_estimate: float
    mesh_rod2_length_estimate: float
    mesh_rod1_radius_estimate: float
    mesh_rod2_radius_estimate: float
    rod1_fixed_nodes: int
    rod2_fixed_nodes: int
    rod1_inner_nodes: int
    rod2_inner_nodes: int
    joint_coupled_nodes: int
    parsed_modes: int
    solid_shapes: tuple[SolidShape, ...]


def load_point_joint_module() -> ModuleType:
    module_path = ANALYSIS_DIR / "solid_fem_coupled_equal_rods_point_joint.py"
    spec = importlib.util.spec_from_file_location("point_joint_fem_helpers_for_mu_plot", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load point-joint FEM helpers from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def circular_shear_coefficient(nu: float) -> float:
    return (6.0 + 12.0 * nu + 6.0 * nu**2) / (7.0 + 12.0 * nu + 4.0 * nu**2)


def section_from_epsilon(epsilon: float, kappa: float | None = None) -> Section:
    radius = 2.0 * float(epsilon) * L_SEGMENT
    area = math.pi * radius**2
    inertia = math.pi * radius**4 / 4.0
    shear_modulus = E / (2.0 * (1.0 + NU))
    return Section(
        radius=radius,
        area=area,
        inertia=inertia,
        shear_modulus=shear_modulus,
        kappa=circular_shear_coefficient(NU) if kappa is None else float(kappa),
    )


def project_omega(lambda_value: float, epsilon: float) -> float:
    return float(epsilon) * float(lambda_value) ** 2


def lambda_from_omega(omega: float, epsilon: float) -> float:
    if omega <= 0.0 or epsilon <= 0.0:
        return float("nan")
    return math.sqrt(float(omega) / float(epsilon))


def segment_lengths(mu: float) -> tuple[float, float]:
    return 1.0 - float(mu), 1.0 + float(mu)


def diameter_to_segment_ratios(epsilon: float, mu: float) -> tuple[float, float, float]:
    l1, l2 = segment_lengths(float(mu))
    diameter = 4.0 * float(epsilon)
    ratio1 = diameter / l1 if l1 > 0.0 else float("inf")
    ratio2 = diameter / l2 if l2 > 0.0 else float("inf")
    return ratio1, ratio2, max(ratio1, ratio2)


def thin_rod_valid(epsilon: float, mu: float) -> bool:
    return diameter_to_segment_ratios(float(epsilon), float(mu))[2] <= THIN_ROD_LIMIT + 1.0e-15


def omega_cutoff(section: Section) -> float:
    return math.sqrt(section.shear_stiffness / (RHO * section.inertia))


def timoshenko_cutoff_ratio(lambda_value: float, epsilon: float) -> float:
    section = section_from_epsilon(float(epsilon))
    return project_omega(float(lambda_value), float(epsilon)) / omega_cutoff(section)


def safe_float_token(value: float) -> str:
    text = f"{float(value):.6g}"
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text.replace("-", "m").replace(".", "p")


def mu_key(value: float) -> float:
    return round(float(value), 10)


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def mul(value: float, vector: tuple[float, float, float]) -> tuple[float, float, float]:
    return (float(value) * vector[0], float(value) * vector[1], float(value) * vector[2])


def norm(vector: tuple[float, float, float]) -> float:
    return math.sqrt(dot(vector, vector))


def unit(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    length = norm(vector)
    if length <= 0.0:
        raise ValueError("zero vector cannot be normalized")
    return (vector[0] / length, vector[1] / length, vector[2] / length)


def geometry_for_beta_mu(beta_deg: float, mu: float) -> Geometry:
    beta = math.radians(float(beta_deg))
    cos_b = math.cos(beta)
    sin_b = math.sin(beta)
    l1, l2 = segment_lengths(float(mu))
    joint = (0.0, 0.0, 0.0)
    rod1_outer = (-l1, 0.0, 0.0)
    rod2_outer = (-l2 * cos_b, -l2 * sin_b, 0.0)
    rod1_tangent = unit(sub(joint, rod1_outer))
    rod2_tangent = unit(sub(joint, rod2_outer))
    rod1_normal = (-rod1_tangent[1], rod1_tangent[0], 0.0)
    rod2_normal = (-rod2_tangent[1], rod2_tangent[0], 0.0)
    return Geometry(
        joint=joint,
        rod1=RodFrame(outer=rod1_outer, joint=joint, tangent=rod1_tangent, in_plane_normal=rod1_normal, length=l1),
        rod2=RodFrame(outer=rod2_outer, joint=joint, tangent=rod2_tangent, in_plane_normal=rod2_normal, length=l2),
        beta_rad=beta,
        mu=float(mu),
    )


def projection_on_rod(point: tuple[float, float, float], rod: RodFrame) -> tuple[float, float]:
    rel = sub(point, rod.outer)
    s = dot(rel, rod.tangent)
    axis_point = add(rod.outer, mul(s, rod.tangent))
    radial = norm(sub(point, axis_point))
    return s, radial


def mesh_size_for_epsilon(epsilon: float, mesh_size_factor: float = 1.0) -> float:
    radius = 2.0 * float(epsilon)
    baseline = min(L_SEGMENT / 40.0, radius / 2.0)
    return baseline * float(mesh_size_factor)


def row_normalized(matrix: np.ndarray) -> np.ndarray:
    out = np.array(matrix, dtype=float, copy=True)
    norms = np.linalg.norm(out, axis=1)
    for index, row_norm in enumerate(norms):
        if row_norm > 0.0 and math.isfinite(float(row_norm)):
            out[index, :] /= row_norm
    return out


def normalized_det(matrix: np.ndarray) -> float:
    return float(np.linalg.det(row_normalized(matrix)))


def analytic_null_vector(matrix: np.ndarray) -> tuple[np.ndarray, float, float]:
    scaled = row_normalized(matrix)
    _u, singular_values, vh = np.linalg.svd(scaled)
    coeff = vh[-1, :].astype(float)
    if coeff.size:
        pivot = int(np.argmax(np.abs(coeff)))
        if coeff[pivot] < 0.0:
            coeff = -coeff
    smallest = float(singular_values[-1])
    ratio = (
        float(singular_values[-1] / singular_values[-2])
        if len(singular_values) >= 2 and abs(float(singular_values[-2])) > 1.0e-14
        else float("nan")
    )
    return coeff, smallest, ratio


def timo_basis(Lambda: float, epsilon: float, section: Section) -> TimoshenkoBasis:
    omega = project_omega(float(Lambda), float(epsilon))
    k_stiff = section.shear_stiffness
    b_stiff = section.bending_stiffness
    m = section.mass_per_length
    j = section.rotary_inertia_per_length
    c2 = omega**2 * (b_stiff * m / k_stiff + j)
    c0 = j * m * omega**4 / k_stiff - m * omega**2
    roots = np.roots([b_stiff, c2, c0])
    warnings: list[str] = []
    scale = max(1.0, float(np.max(np.abs(roots))))
    if np.max(np.abs(np.imag(roots))) > 1.0e-8 * scale:
        warnings.append(f"Timoshenko lambda^2 roots have non-negligible imaginary parts at Lambda={Lambda:.8g}.")
    z_values = np.real(roots)
    positives = [float(z) for z in z_values if z > 1.0e-12]
    negatives = [float(z) for z in z_values if z < -1.0e-12]
    if len(positives) != 1 or len(negatives) != 1:
        warnings.append(f"Expected one positive and one negative lambda^2 root at Lambda={Lambda:.8g}.")
        z_positive = float(np.max(z_values))
        z_negative = float(np.min(z_values))
        if z_positive <= 0.0 or z_negative >= 0.0:
            raise ValueError(warnings[-1])
    else:
        z_positive = positives[0]
        z_negative = negatives[0]
    a = math.sqrt(z_positive)
    b = math.sqrt(-z_negative)
    h = a + m * omega**2 / (k_stiff * a)
    q = b - m * omega**2 / (k_stiff * b)
    if abs(q) <= 1.0e-12 * max(1.0, abs(b)):
        warnings.append(f"Timoshenko trigonometric coupling q is close to zero at Lambda={Lambda:.8g}.")
    return TimoshenkoBasis(a=float(a), b=float(b), h=float(h), q=float(q), warnings=tuple(warnings))


def bending_endpoint_columns(x: float, basis: TimoshenkoBasis) -> dict[str, np.ndarray]:
    a = basis.a
    b = basis.b
    h = basis.h
    q = basis.q
    ax = a * float(x)
    bx = b * float(x)
    cosh_ax = np.cosh(ax)
    sinh_ax = np.sinh(ax)
    cos_bx = np.cos(bx)
    sin_bx = np.sin(bx)
    return {
        "w": np.array([cosh_ax - cos_bx, sinh_ax - (h / q) * sin_bx], dtype=float),
        "psi": np.array([h * sinh_ax + q * sin_bx, h * (cosh_ax - cos_bx)], dtype=float),
        "w_prime": np.array([a * sinh_ax + b * sin_bx, a * cosh_ax - (h / q) * b * cos_bx], dtype=float),
        "psi_prime": np.array([h * a * cosh_ax + q * b * cos_bx, h * (a * sinh_ax + b * sin_bx)], dtype=float),
    }


def endpoint_columns(x: float, theta: float, basis: TimoshenkoBasis, section: Section) -> dict[str, np.ndarray]:
    bend = bending_endpoint_columns(float(x), basis)
    w = np.zeros(3, dtype=float)
    psi = np.zeros(3, dtype=float)
    w_prime = np.zeros(3, dtype=float)
    psi_prime = np.zeros(3, dtype=float)
    u = np.zeros(3, dtype=float)
    u_prime = np.zeros(3, dtype=float)
    w[:2] = bend["w"]
    psi[:2] = bend["psi"]
    w_prime[:2] = bend["w_prime"]
    psi_prime[:2] = bend["psi_prime"]
    u[2] = np.sin(theta * float(x))
    u_prime[2] = theta * np.cos(theta * float(x))
    q_force = section.shear_stiffness * (w_prime - psi)
    moment = section.bending_stiffness * psi_prime
    normal = E * section.area * u_prime
    return {"w": w, "u": u, "psi": psi, "Q": q_force, "M": moment, "N": normal}


def timo_coupling_matrix(
    Lambda: float,
    beta_deg: float,
    mu: float,
    epsilon: float,
    *,
    kappa: float | None = None,
) -> tuple[np.ndarray, tuple[str, ...]]:
    section = section_from_epsilon(float(epsilon), kappa=kappa)
    basis = timo_basis(float(Lambda), float(epsilon), section)
    theta = project_omega(float(Lambda), float(epsilon))
    l1, l2 = segment_lengths(float(mu))
    rod1 = endpoint_columns(l1, theta, basis, section)
    rod2 = endpoint_columns(-l2, theta, basis, section)
    cb = math.cos(math.radians(float(beta_deg)))
    sb = math.sin(math.radians(float(beta_deg)))
    matrix = np.zeros((6, 6), dtype=float)
    matrix[0, 0:3] = rod1["w"]
    matrix[0, 3:6] = -cb * rod2["w"] + sb * rod2["u"]
    matrix[1, 0:3] = rod1["u"]
    matrix[1, 3:6] = -sb * rod2["w"] - cb * rod2["u"]
    matrix[2, 0:3] = rod1["psi"]
    matrix[2, 3:6] = -rod2["psi"]
    matrix[3, 0:3] = rod1["M"]
    matrix[3, 3:6] = -rod2["M"]
    matrix[4, 0:3] = -rod1["Q"]
    matrix[4, 3:6] = cb * rod2["Q"] - sb * rod2["N"]
    matrix[5, 0:3] = rod1["N"]
    matrix[5, 3:6] = -sb * rod2["Q"] - cb * rod2["N"]
    return matrix, basis.warnings


def timo_det_for_scan(Lambda: float, beta_deg: float, mu: float, epsilon: float) -> float:
    try:
        matrix, _warnings = timo_coupling_matrix(float(Lambda), float(beta_deg), float(mu), float(epsilon))
    except (FloatingPointError, ValueError, OverflowError):
        return float("nan")
    if not np.all(np.isfinite(matrix)):
        return float("nan")
    return normalized_det(matrix)


def find_roots_by_sign_scan(
    func: Callable[[float], float],
    n_roots: int,
    *,
    start: float,
    upper: float,
    scan_step: float,
    grow_factor: float = 1.35,
    max_tries: int = 8,
) -> tuple[np.ndarray, list[str]]:
    warnings: list[str] = []
    roots: list[float] = []
    current_upper = float(upper)
    for _ in range(max_tries):
        roots.clear()
        grid = np.arange(float(start), current_upper + scan_step, scan_step)
        values = np.array([func(float(x)) for x in grid], dtype=float)
        for left, right, f_left, f_right in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
            if not (np.isfinite(f_left) and np.isfinite(f_right)):
                continue
            candidate: float | None = None
            if f_left == 0.0:
                candidate = float(left)
            elif f_left * f_right < 0.0:
                try:
                    candidate = float(brentq(func, float(left), float(right), xtol=1.0e-12, rtol=1.0e-12, maxiter=120))
                except ValueError as exc:
                    warnings.append(f"Skipped bracket [{left:.8g}, {right:.8g}]: {exc}")
                    continue
            if candidate is None:
                continue
            if roots and abs(candidate - roots[-1]) <= ROOT_DEDUP_TOL:
                continue
            roots.append(candidate)
            if len(roots) >= n_roots:
                return np.array(roots[:n_roots], dtype=float), warnings
        current_upper *= grow_factor
    warnings.append(f"Found only {len(roots)} Timoshenko roots below Lambda={current_upper:.8g}.")
    out = np.full(int(n_roots), np.nan, dtype=float)
    out[: len(roots)] = roots[:n_roots]
    return out, warnings


def timo_sorted_roots(beta_deg: float, mu: float, epsilon: float, n_roots: int) -> tuple[np.ndarray, list[str]]:
    l1, _l2 = segment_lengths(float(mu))
    upper = max(18.0, 1.35 * math.pi * (float(n_roots) + 4) / max(l1, 0.08))
    return find_roots_by_sign_scan(
        lambda value: timo_det_for_scan(value, float(beta_deg), float(mu), float(epsilon)),
        n_roots,
        start=ROOT_SCAN_START,
        upper=upper,
        scan_step=ROOT_SCAN_STEP,
    )


def slice_coordinates() -> np.ndarray:
    return (np.arange(N_SLICES_PER_ROD, dtype=float) + 0.5) / float(N_SLICES_PER_ROD)


def global_centerline_vector(
    beta_deg: float,
    mu: float,
    u_left: np.ndarray,
    w_left: np.ndarray,
    u_right: np.ndarray,
    w_right: np.ndarray,
) -> np.ndarray:
    geom = geometry_for_beta_mu(float(beta_deg), float(mu))
    parts: list[np.ndarray] = []
    for u_values, w_values, rod in (
        (np.asarray(u_left, dtype=float), np.asarray(w_left, dtype=float), geom.rod1),
        (np.asarray(u_right, dtype=float), np.asarray(w_right, dtype=float), geom.rod2),
    ):
        ux = u_values * rod.tangent[0] + w_values * rod.in_plane_normal[0]
        uy = u_values * rod.tangent[1] + w_values * rod.in_plane_normal[1]
        interleaved = np.empty(2 * len(ux), dtype=float)
        interleaved[0::2] = ux
        interleaved[1::2] = uy
        parts.append(interleaved)
    return np.concatenate(parts)


def normalize_vector(vector: np.ndarray) -> np.ndarray:
    out = np.asarray(vector, dtype=float)
    scale = float(np.linalg.norm(out))
    if scale <= 1.0e-28 or not math.isfinite(scale):
        return np.full_like(out, np.nan, dtype=float)
    return out / scale


def eb_centerline_vector(beta_deg: float, mu: float, epsilon: float, Lambda: float) -> np.ndarray:
    matrix = assemble_clamped_coupled_matrix(float(Lambda), math.radians(float(beta_deg)), float(mu), float(epsilon))
    coeff, _smallest, _ratio = analytic_null_vector(matrix)
    a1, b1, a2, b2, p1, p2 = [float(value) for value in coeff]
    xi = slice_coordinates()
    l1, l2 = segment_lengths(float(mu))
    z1 = float(Lambda) * l1 * xi
    theta1 = float(epsilon) * float(Lambda) ** 2 * l1 * xi
    z2 = -float(Lambda) * l2 * xi
    theta2 = -float(epsilon) * float(Lambda) ** 2 * l2 * xi
    w_left = a1 * (np.cos(z1) - np.cosh(z1)) + b1 * (np.sin(z1) - np.sinh(z1))
    u_left = p1 * np.sin(theta1)
    w_right = a2 * (np.cos(z2) - np.cosh(z2)) + b2 * (np.sin(z2) - np.sinh(z2))
    u_right = p2 * np.sin(theta2)
    return normalize_vector(global_centerline_vector(float(beta_deg), float(mu), u_left, w_left, u_right, w_right))


def timo_centerline_vector(beta_deg: float, mu: float, epsilon: float, Lambda: float) -> tuple[np.ndarray, tuple[str, ...]]:
    matrix, matrix_warnings = timo_coupling_matrix(float(Lambda), float(beta_deg), float(mu), float(epsilon))
    coeff, _smallest, _ratio = analytic_null_vector(matrix)
    section = section_from_epsilon(float(epsilon))
    basis = timo_basis(float(Lambda), float(epsilon), section)
    a1, b1, p1, a2, b2, p2 = [float(value) for value in coeff]
    xi = slice_coordinates()
    l1, l2 = segment_lengths(float(mu))
    theta = project_omega(float(Lambda), float(epsilon))
    w_left = np.zeros_like(xi)
    w_right = np.zeros_like(xi)
    for index, xi_value in enumerate(xi):
        left_cols = bending_endpoint_columns(l1 * float(xi_value), basis)["w"]
        right_cols = bending_endpoint_columns(-l2 * float(xi_value), basis)["w"]
        w_left[index] = a1 * float(left_cols[0]) + b1 * float(left_cols[1])
        w_right[index] = a2 * float(right_cols[0]) + b2 * float(right_cols[1])
    u_left = p1 * np.sin(theta * l1 * xi)
    u_right = p2 * np.sin(-theta * l2 * xi)
    warnings = tuple(dict.fromkeys([*matrix_warnings, *basis.warnings]))
    return normalize_vector(global_centerline_vector(float(beta_deg), float(mu), u_left, w_left, u_right, w_right)), warnings


def centerline_mac(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    denominator = float(np.dot(a, a) * np.dot(b, b))
    if denominator <= 1.0e-28 or not math.isfinite(denominator):
        return float("nan")
    dot_value = float(np.dot(a, b))
    return float((dot_value * dot_value) / denominator)


def mac_strength(value: object) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "missing"
    if not math.isfinite(number):
        return "missing"
    if number >= MAC_STRONG_THRESHOLD:
        return "strong"
    if number >= MAC_MODERATE_THRESHOLD:
        return "moderate"
    return "weak"


def solve_timo_states(beta_deg: float, mu: float, epsilon: float, n_roots: int) -> tuple[list[dict[str, object]], list[str]]:
    roots, warnings = timo_sorted_roots(float(beta_deg), float(mu), float(epsilon), int(n_roots))
    states: list[dict[str, object]] = []
    for index, root in enumerate(roots, start=1):
        if not math.isfinite(float(root)):
            continue
        try:
            vector, shape_warnings = timo_centerline_vector(float(beta_deg), float(mu), float(epsilon), float(root))
        except Exception as exc:
            warnings.append(f"beta={beta_deg:g}, mu={mu:g}, sorted={index}: Timoshenko vector failed: {exc}")
            continue
        warnings.extend(f"beta={beta_deg:g}, mu={mu:g}, sorted={index}: {item}" for item in shape_warnings)
        states.append(
            {
                "sorted_index": int(index),
                "Lambda": float(root),
                "vector": vector,
            }
        )
    return states, warnings


def assign_previous_to_current(
    branch_ids: Sequence[str],
    previous_vectors: dict[str, np.ndarray],
    current_states: Sequence[dict[str, object]],
) -> dict[str, tuple[dict[str, object], float]]:
    if not branch_ids or not current_states:
        return {}
    cost = np.ones((len(branch_ids), len(current_states)), dtype=float)
    mac = np.zeros_like(cost)
    for row, branch_id in enumerate(branch_ids):
        previous = previous_vectors[branch_id]
        for col, state in enumerate(current_states):
            value = centerline_mac(previous, np.asarray(state["vector"], dtype=float))
            mac[row, col] = value if math.isfinite(value) else 0.0
            cost[row, col] = 1.0 - mac[row, col]
    row_ind, col_ind = linear_sum_assignment(cost)
    assigned: dict[str, tuple[dict[str, object], float]] = {}
    for row, col in zip(row_ind, col_ind):
        assigned[branch_ids[int(row)]] = (current_states[int(col)], float(mac[int(row), int(col)]))
    return assigned


def track_timoshenko_descendants(mu_values: Sequence[float]) -> tuple[dict[str, dict[float, TimoTrackedPoint]], list[str]]:
    warnings: list[str] = []
    branch_ids = [branch_id_from_base_sorted_index(index) for index in range(1, N_BRANCHES + 1)]
    beta_path = list(np.linspace(0.0, BETA_DEG, 16))
    path: list[tuple[float, float]] = [(float(beta), 0.0) for beta in beta_path]
    for mu in mu_values:
        if abs(float(mu)) <= 1.0e-12:
            continue
        path.append((BETA_DEG, float(mu)))

    first_states, first_warnings = solve_timo_states(0.0, 0.0, EPSILON, N_SOLVE_ANALYTIC)
    warnings.extend(first_warnings)
    if len(first_states) < N_BRANCHES:
        raise RuntimeError(f"Only {len(first_states)} Timoshenko roots found at beta=0, mu=0.")

    previous_vectors: dict[str, np.ndarray] = {}
    tracked: dict[str, dict[float, TimoTrackedPoint]] = {branch_id: {} for branch_id in branch_ids}
    for branch_id, state in zip(branch_ids, first_states[:N_BRANCHES]):
        previous_vectors[branch_id] = np.asarray(state["vector"], dtype=float)
        if BETA_DEG == 0.0 and 0.0 in mu_values:
            tracked[branch_id][0.0] = TimoTrackedPoint(
                branch_id=branch_id,
                base_sorted_index=int(branch_id.rsplit("_", 1)[1]),
                beta_deg=0.0,
                mu=0.0,
                current_sorted_index=int(state["sorted_index"]),
                Lambda=float(state["Lambda"]),
                mac_to_previous=1.0,
                vector=np.asarray(state["vector"], dtype=float),
            )

    for beta_deg, mu in path[1:]:
        states, root_warnings = solve_timo_states(beta_deg, mu, EPSILON, N_SOLVE_ANALYTIC)
        warnings.extend(root_warnings)
        assignment = assign_previous_to_current(branch_ids, previous_vectors, states)
        if len(assignment) < len(branch_ids):
            warnings.append(f"Timoshenko assignment incomplete at beta={beta_deg:g}, mu={mu:g}.")
        for branch_id in branch_ids:
            if branch_id not in assignment:
                continue
            state, mac_value = assignment[branch_id]
            previous_vectors[branch_id] = np.asarray(state["vector"], dtype=float)
            if abs(beta_deg - BETA_DEG) <= 1.0e-10 and any(abs(float(mu) - float(item)) <= 1.0e-10 for item in mu_values):
                tracked[branch_id][mu_key(mu)] = TimoTrackedPoint(
                    branch_id=branch_id,
                    base_sorted_index=int(branch_id.rsplit("_", 1)[1]),
                    beta_deg=float(beta_deg),
                    mu=float(mu),
                    current_sorted_index=int(state["sorted_index"]),
                    Lambda=float(state["Lambda"]),
                    mac_to_previous=float(mac_value),
                    vector=np.asarray(state["vector"], dtype=float),
                )

    return tracked, sorted(set(warnings))


def fem_case_paths(pj: ModuleType, mu: float) -> object:
    beta = safe_float_token(BETA_DEG)
    eps = safe_float_token(EPSILON)
    token_mu = safe_float_token(float(mu))
    case_dir = FEM_OUTPUT_ROOT / f"beta_{beta}" / f"eps_{eps}" / f"mu_{token_mu}"
    stem = f"coupled_equal_rods_beta{beta}_eps{eps}_mu{token_mu}_planar"
    return pj.CasePaths(
        case_dir=case_dir,
        geo=case_dir / f"{stem}.geo",
        msh=case_dir / f"{stem}_mesh.msh",
        gmsh_inp=case_dir / f"{stem}_mesh.inp",
        ccx_mesh_inp=case_dir / f"{stem}_ccx_mesh.inp",
        ccx_modal_inp=case_dir / f"{stem}_ccx_modal.inp",
        ccx_dat=case_dir / f"{stem}_ccx_modal.dat",
        ccx_frd=case_dir / f"{stem}_ccx_modal.frd",
        ccx_stdout=case_dir / f"{stem}_ccx_modal.stdout.txt",
        ccx_stderr=case_dir / f"{stem}_ccx_modal.stderr.txt",
    )


def write_gmsh_geo(paths: object, mu: float) -> str:
    radius = 2.0 * EPSILON
    mesh_size = mesh_size_for_epsilon(EPSILON, FEM_MESH_SIZE_FACTOR)
    geom = geometry_for_beta_mu(BETA_DEG, float(mu))
    cos_b = math.cos(geom.beta_rad)
    sin_b = math.sin(geom.beta_rad)
    l1, l2 = segment_lengths(float(mu))
    paths.case_dir.mkdir(parents=True, exist_ok=True)
    paths.geo.write_text(
        f"""// Diagnostic-only planar point-joint 3D solid FEM mesh for Lambda(mu).
// beta = {BETA_DEG:.17g} deg, epsilon = {EPSILON:.17g}, eta = {ETA:.17g}, mu = {float(mu):.17g}
SetFactory("OpenCASCADE");

L1 = {l1:.17g};
L2 = {l2:.17g};
R = {radius:.17g};
h = {mesh_size:.17g};
CosB = {cos_b:.17g};
SinB = {sin_b:.17g};

Cylinder(1) = {{-L1, 0, 0, L1, 0, 0, R, 2*Pi}};
Cylinder(2) = {{-L2*CosB, -L2*SinB, 0, L2*CosB, L2*SinB, 0, R, 2*Pi}};
Physical Volume("ROD1_SOLID", 1) = {{1}};
Physical Volume("ROD2_SOLID", 2) = {{2}};

Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Mesh.ElementOrder = 2;
Mesh.Optimize = 1;
Mesh.HighOrderOptimize = 2;
Mesh.MshFileVersion = 4.1;
""",
        encoding="utf-8",
    )
    return "planar point-joint variable-length separate-cylinder diagnostic"


def rod_assignment_cost(point: tuple[float, float, float], rod: RodFrame) -> float:
    s, radial = projection_on_rod(point, rod)
    outside = max(0.0, -s, s - rod.length)
    return radial + outside


def rod_node_memberships(mesh: object, mu: float) -> tuple[dict[int, set[int]], dict[int, int]]:
    geom = geometry_for_beta_mu(BETA_DEG, float(mu))
    memberships: dict[int, set[int]] = {1: set(), 2: set()}
    element_counts: dict[int, int] = {1: 0, 2: 0}
    for connectivity in mesh.solid_elements.values():
        points = [mesh.nodes[node_id] for node_id in connectivity if node_id in mesh.nodes]
        if not points:
            continue
        count = float(len(points))
        centroid = (
            sum(point[0] for point in points) / count,
            sum(point[1] for point in points) / count,
            sum(point[2] for point in points) / count,
        )
        cost1 = rod_assignment_cost(centroid, geom.rod1)
        cost2 = rod_assignment_cost(centroid, geom.rod2)
        rod_index = 1 if cost1 <= cost2 else 2
        element_counts[rod_index] += 1
        memberships[rod_index].update(connectivity)
    return memberships, element_counts


def coordinate_case_node_sets(mesh: object, mu: float) -> tuple[list[int], list[int], list[int], list[int]]:
    geom = geometry_for_beta_mu(BETA_DEG, float(mu))
    memberships, _element_counts = rod_node_memberships(mesh, float(mu))
    radius = 2.0 * EPSILON
    plane_tol = max(mesh_size_for_epsilon(EPSILON, FEM_MESH_SIZE_FACTOR) / 50.0, 1.0e-6)
    radial_tol = 1.02 * radius + plane_tol
    rod1_outer: list[int] = []
    rod2_outer: list[int] = []
    rod1_inner: list[int] = []
    rod2_inner: list[int] = []
    for rod_index, rod, outer, inner in (
        (1, geom.rod1, rod1_outer, rod1_inner),
        (2, geom.rod2, rod2_outer, rod2_inner),
    ):
        for node_id in memberships[rod_index]:
            point = mesh.nodes.get(node_id)
            if point is None:
                continue
            s, radial = projection_on_rod(point, rod)
            if abs(s) <= plane_tol and radial <= radial_tol:
                outer.append(node_id)
            if abs(s - rod.length) <= plane_tol and radial <= radial_tol:
                inner.append(node_id)
    return sorted(rod1_outer), sorted(rod2_outer), sorted(rod1_inner), sorted(rod2_inner)


def has_negative_jacobian_warning(messages: Sequence[str]) -> bool:
    for message in messages:
        lower = str(message).lower()
        if "jac. < 0" in lower or "negative jac" in lower:
            return True
        if "worst distortion = -" in lower:
            return True
    return False


def estimate_mesh_geometry(mesh: object, mu: float) -> tuple[float, float, float, float]:
    geom = geometry_for_beta_mu(BETA_DEG, float(mu))
    memberships, _element_counts = rod_node_memberships(mesh, float(mu))
    estimates: list[tuple[float, float]] = []
    for rod_index, rod in ((1, geom.rod1), (2, geom.rod2)):
        s_values: list[float] = []
        radial_values: list[float] = []
        for node_id in memberships[rod_index]:
            point = mesh.nodes.get(node_id)
            if point is None:
                continue
            s, radial = projection_on_rod(point, rod)
            s_values.append(float(s))
            if -1.0e-7 <= s <= rod.length + 1.0e-7:
                radial_values.append(float(radial))
        if not s_values:
            estimates.append((float("nan"), float("nan")))
            continue
        length_estimate = max(s_values) - min(s_values)
        radius_estimate = max(radial_values) if radial_values else float("nan")
        estimates.append((float(length_estimate), float(radius_estimate)))
    return estimates[0][0], estimates[1][0], estimates[0][1], estimates[1][1]


def write_calculix_inputs(pj: ModuleType, paths: object, mesh: object, mu: float) -> tuple[int, int, int, int, int]:
    rod1_fixed, rod2_fixed, rod1_inner, rod2_inner = coordinate_case_node_sets(mesh, float(mu))
    joint_coupled = sorted(set(rod1_inner) | set(rod2_inner))
    all_nodes = sorted(mesh.nodes)
    planar_independent_nodes = sorted(set(all_nodes) - set(joint_coupled))
    joint_ref_node_id = max(mesh.nodes.keys(), default=0) + 1
    pj.write_calculix_mesh_include(paths, mesh)
    lines: list[str] = [
        "** Diagnostic-only planar-constrained point-joint 3D solid FEM input.",
        "** Separate variable-length cylinders connected by one shared rigid reference node.",
        "** Planar constraint uses MPC-compatible independent-node set plus JOINT_REF,3,5,0.",
        f"*INCLUDE, INPUT={paths.ccx_mesh_inp.name}",
        "*NODE",
        f"{joint_ref_node_id}, 0.0, 0.0, 0.0",
        "*NSET, NSET=JOINT_REF",
        f"{joint_ref_node_id}",
        "*NSET, NSET=ROD1_OUTER_FIXED",
        *pj.format_id_lines(rod1_fixed),
        "*NSET, NSET=ROD2_OUTER_FIXED",
        *pj.format_id_lines(rod2_fixed),
        "*NSET, NSET=ROD1_INNER_COUPLED",
        *pj.format_id_lines(rod1_inner),
        "*NSET, NSET=ROD2_INNER_COUPLED",
        *pj.format_id_lines(rod2_inner),
        "*NSET, NSET=JOINT_COUPLED_ALL",
        *pj.format_id_lines(joint_coupled),
        "*NSET, NSET=PLANAR_INDEPENDENT_SOLID_NODES",
        *pj.format_id_lines(planar_independent_nodes),
        "*MATERIAL, NAME=MAT",
        "*ELASTIC",
        f"{E:.17g}, {NU:.17g}",
        "*DENSITY",
        f"{RHO:.17g}",
        "*SOLID SECTION, ELSET=SOLID, MATERIAL=MAT",
        f"*RIGID BODY, NSET=JOINT_COUPLED_ALL, REF NODE={joint_ref_node_id}",
        "*BOUNDARY",
        "PLANAR_INDEPENDENT_SOLID_NODES, 3, 3, 0",
        "JOINT_REF, 3, 5, 0",
        "ROD1_OUTER_FIXED, 1, 3, 0",
        "ROD2_OUTER_FIXED, 1, 3, 0",
        "*STEP",
        "*FREQUENCY",
        f"{N_SOLID_MODES}",
        "*NODE FILE",
        "U",
        "*END STEP",
    ]
    paths.ccx_modal_inp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(rod1_fixed), len(rod2_fixed), len(rod1_inner), len(rod2_inner), len(joint_coupled)


def fill_missing_slice_means(values: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, int, int]:
    out = np.array(values, dtype=float, copy=True)
    count_values = np.asarray(counts, dtype=float)
    filled = int(np.count_nonzero(count_values > 0.0))
    empty = int(count_values.size - filled)
    known = np.where(count_values > 0.0)[0]
    if known.size == 0:
        return np.zeros_like(out), 0, int(count_values.size)
    for index in range(count_values.size):
        if count_values[index] > 0.0:
            out[index, :] /= count_values[index]
            continue
        nearest = int(known[np.argmin(np.abs(known - index))])
        out[index, :] = values[nearest, :] / count_values[nearest]
    return out, filled, empty


def solid_centerline_vector(mesh: object, mode_shape: dict[int, tuple[float, float, float]], mu: float) -> tuple[np.ndarray, int, int, float]:
    geom = geometry_for_beta_mu(BETA_DEG, float(mu))
    memberships, _element_counts = rod_node_memberships(mesh, float(mu))
    rod_vectors: list[np.ndarray] = []
    total_filled = 0
    total_empty = 0
    total_energy = 0.0
    out_energy = 0.0
    for rod_index, rod in ((1, geom.rod1), (2, geom.rod2)):
        values = np.zeros((N_SLICES_PER_ROD, 2), dtype=float)
        counts = np.zeros(N_SLICES_PER_ROD, dtype=float)
        for node_id in memberships[rod_index]:
            displacement = mode_shape.get(node_id)
            point = mesh.nodes.get(node_id)
            if displacement is None or point is None:
                continue
            ux, uy, uz = displacement
            total_energy += ux * ux + uy * uy + uz * uz
            out_energy += uz * uz
            s, _radial = projection_on_rod(point, rod)
            if s < -1.0e-8 or s > rod.length + 1.0e-8:
                continue
            slice_index = int(math.floor(max(0.0, min(0.999999999999, s / rod.length)) * N_SLICES_PER_ROD))
            values[slice_index, 0] += float(ux)
            values[slice_index, 1] += float(uy)
            counts[slice_index] += 1.0
        means, filled, empty = fill_missing_slice_means(values, counts)
        rod_vectors.append(means.reshape(-1))
        total_filled += filled
        total_empty += empty
    out_fraction = out_energy / total_energy if total_energy > 0.0 else float("nan")
    return normalize_vector(np.concatenate(rod_vectors)), total_filled, total_empty, float(out_fraction)


def run_fem_cases() -> tuple[list[FemCaseResult], list[str]]:
    warnings: list[str] = []
    pj = load_point_joint_module()
    gmsh_exe, gmsh_note = pj.resolve_gmsh_executable()
    ccx_exe, ccx_note = pj.resolve_ccx_executable()
    if gmsh_exe is None or ccx_exe is None:
        warnings.append(f"FEM skipped: gmsh={gmsh_note}; ccx={ccx_note}")
        return [], warnings

    results: list[FemCaseResult] = []
    for mu in FEM_MU_VALUES:
        l1, l2 = segment_lengths(float(mu))
        radius = 2.0 * EPSILON
        ratio1, ratio2, _ratio_max = diameter_to_segment_ratios(EPSILON, float(mu))
        paths = fem_case_paths(pj, float(mu))
        write_gmsh_geo(paths, float(mu))
        mesh_ok, mesh_message, gmsh_messages = pj.generate_mesh_with_gmsh_cli(paths, gmsh_exe)
        warnings.extend(f"mu={float(mu):g}: {item}" for item in gmsh_messages)
        negative_jacobian = has_negative_jacobian_warning(gmsh_messages)
        mesh_quality_status = "mesh_quality_warning_negative_jacobian" if negative_jacobian else "ok"
        if not mesh_ok:
            results.append(
                FemCaseResult(
                    mu=float(mu),
                    case_dir=paths.case_dir,
                    success=False,
                    message=mesh_message,
                    gmsh_messages=tuple(gmsh_messages),
                    mesh_quality_status=mesh_quality_status,
                    negative_jacobian_warning=negative_jacobian,
                    nodes=0,
                    solid_elements=0,
                    l1=l1,
                    l2=l2,
                    radius=radius,
                    diameter_to_rod1_length=ratio1,
                    diameter_to_rod2_length=ratio2,
                    mesh_rod1_length_estimate=float("nan"),
                    mesh_rod2_length_estimate=float("nan"),
                    mesh_rod1_radius_estimate=float("nan"),
                    mesh_rod2_radius_estimate=float("nan"),
                    rod1_fixed_nodes=0,
                    rod2_fixed_nodes=0,
                    rod1_inner_nodes=0,
                    rod2_inner_nodes=0,
                    joint_coupled_nodes=0,
                    parsed_modes=0,
                    solid_shapes=(),
                )
            )
            continue
        mesh = pj.read_gmsh_inp_mesh_data(paths.gmsh_inp)
        mesh_l1, mesh_l2, mesh_r1, mesh_r2 = estimate_mesh_geometry(mesh, float(mu))
        r1_fixed, r2_fixed, r1_inner, r2_inner, joint_coupled = write_calculix_inputs(pj, paths, mesh, float(mu))
        ccx_result = pj.run_calculix(paths, ccx_exe, BETA_DEG, EPSILON)
        if not ccx_result.success:
            results.append(
                FemCaseResult(
                    mu=float(mu),
                    case_dir=paths.case_dir,
                    success=False,
                    message=ccx_result.message,
                    gmsh_messages=tuple(gmsh_messages),
                    mesh_quality_status=mesh_quality_status,
                    negative_jacobian_warning=negative_jacobian,
                    nodes=len(mesh.nodes),
                    solid_elements=len(mesh.solid_elements),
                    l1=l1,
                    l2=l2,
                    radius=radius,
                    diameter_to_rod1_length=ratio1,
                    diameter_to_rod2_length=ratio2,
                    mesh_rod1_length_estimate=mesh_l1,
                    mesh_rod2_length_estimate=mesh_l2,
                    mesh_rod1_radius_estimate=mesh_r1,
                    mesh_rod2_radius_estimate=mesh_r2,
                    rod1_fixed_nodes=r1_fixed,
                    rod2_fixed_nodes=r2_fixed,
                    rod1_inner_nodes=r1_inner,
                    rod2_inner_nodes=r2_inner,
                    joint_coupled_nodes=joint_coupled,
                    parsed_modes=0,
                    solid_shapes=(),
                )
            )
            continue
        parse_result = pj.parse_calculix_frd_mode_shapes(paths, BETA_DEG, EPSILON)
        if not parse_result.success:
            results.append(
                FemCaseResult(
                    mu=float(mu),
                    case_dir=paths.case_dir,
                    success=False,
                    message=parse_result.message,
                    gmsh_messages=tuple(gmsh_messages),
                    mesh_quality_status=mesh_quality_status,
                    negative_jacobian_warning=negative_jacobian,
                    nodes=len(mesh.nodes),
                    solid_elements=len(mesh.solid_elements),
                    l1=l1,
                    l2=l2,
                    radius=radius,
                    diameter_to_rod1_length=ratio1,
                    diameter_to_rod2_length=ratio2,
                    mesh_rod1_length_estimate=mesh_l1,
                    mesh_rod2_length_estimate=mesh_l2,
                    mesh_rod1_radius_estimate=mesh_r1,
                    mesh_rod2_radius_estimate=mesh_r2,
                    rod1_fixed_nodes=r1_fixed,
                    rod2_fixed_nodes=r2_fixed,
                    rod1_inner_nodes=r1_inner,
                    rod2_inner_nodes=r2_inner,
                    joint_coupled_nodes=joint_coupled,
                    parsed_modes=0,
                    solid_shapes=(),
                )
            )
            continue
        shapes: list[SolidShape] = []
        for mode_index, omega in enumerate(ccx_result.parsed_omegas[:N_SOLID_MODES], start=1):
            mode_shape = parse_result.mode_shapes.get(mode_index)
            if mode_shape is None:
                continue
            vector, filled, empty, out_fraction = solid_centerline_vector(mesh, mode_shape, float(mu))
            shapes.append(
                SolidShape(
                    mu=float(mu),
                    solid_mode=int(mode_index),
                    omega=float(omega),
                    Lambda=lambda_from_omega(float(omega), EPSILON),
                    vector=vector,
                    filled_slice_count=filled,
                    empty_slice_count=empty,
                    out_of_plane_fraction=out_fraction,
                )
            )
        results.append(
            FemCaseResult(
                mu=float(mu),
                case_dir=paths.case_dir,
                success=True,
                message=ccx_result.message,
                gmsh_messages=tuple(gmsh_messages),
                mesh_quality_status=mesh_quality_status,
                negative_jacobian_warning=negative_jacobian,
                nodes=len(mesh.nodes),
                solid_elements=len(mesh.solid_elements),
                l1=l1,
                l2=l2,
                radius=radius,
                diameter_to_rod1_length=ratio1,
                diameter_to_rod2_length=ratio2,
                mesh_rod1_length_estimate=mesh_l1,
                mesh_rod2_length_estimate=mesh_l2,
                mesh_rod1_radius_estimate=mesh_r1,
                mesh_rod2_radius_estimate=mesh_r2,
                rod1_fixed_nodes=r1_fixed,
                rod2_fixed_nodes=r2_fixed,
                rod1_inner_nodes=r1_inner,
                rod2_inner_nodes=r2_inner,
                joint_coupled_nodes=joint_coupled,
                parsed_modes=len(shapes),
                solid_shapes=tuple(shapes),
            )
        )
    return results, sorted(set(warnings))


def build_eb_tracking() -> tuple[dict[str, dict[float, float]], dict[str, dict[float, np.ndarray]], list[str]]:
    branch_ids = [branch_id_from_base_sorted_index(index) for index in range(1, N_BRANCHES + 1)]
    tracking = track_mu_sweep(
        epsilon=EPSILON,
        beta=BETA_DEG,
        mu_values=[float(value) for value in ANALYTIC_MU_VALUES],
        n_track=N_BRANCHES,
        n_solve=N_SOLVE_ANALYTIC,
        allow_low_mac=True,
        required_branch_ids=branch_ids,
    )
    lambdas: dict[str, dict[float, float]] = {branch_id: {} for branch_id in branch_ids}
    vectors: dict[str, dict[float, np.ndarray]] = {branch_id: {} for branch_id in branch_ids}
    for branch_id in branch_ids:
        for point in tracking.points_for_branch(branch_id):
            if abs(float(point.beta) - BETA_DEG) > 1.0e-10:
                continue
            key = mu_key(point.mu)
            lambdas[branch_id][key] = float(point.Lambda)
            vectors[branch_id][key] = eb_centerline_vector(BETA_DEG, float(point.mu), EPSILON, float(point.Lambda))
    return lambdas, vectors, list(tracking.warnings)


def build_analytic_rows(
    eb_lambdas: dict[str, dict[float, float]],
    timo_points: dict[str, dict[float, TimoTrackedPoint]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    branch_ids = [branch_id_from_base_sorted_index(index) for index in range(1, N_BRANCHES + 1)]
    for model in ("Euler-Bernoulli", "Timoshenko"):
        for branch_index, branch_id in enumerate(branch_ids, start=1):
            for mu in ANALYTIC_MU_VALUES:
                key = mu_key(mu)
                if model == "Euler-Bernoulli":
                    lambda_value = eb_lambdas.get(branch_id, {}).get(key, float("nan"))
                    current_sorted = ""
                    mac_to_previous = ""
                    omega_over_cutoff = ""
                    below_cutoff = ""
                else:
                    point = timo_points.get(branch_id, {}).get(key)
                    lambda_value = point.Lambda if point is not None else float("nan")
                    current_sorted = point.current_sorted_index if point is not None else ""
                    mac_to_previous = point.mac_to_previous if point is not None else ""
                    omega_over_cutoff = (
                        timoshenko_cutoff_ratio(float(lambda_value), EPSILON)
                        if math.isfinite(float(lambda_value))
                        else float("nan")
                    )
                    below_cutoff = bool(math.isfinite(float(omega_over_cutoff)) and float(omega_over_cutoff) < 1.0)
                ratio1, ratio2, ratio_max = diameter_to_segment_ratios(EPSILON, float(mu))
                rows.append(
                    {
                        "row_kind": "analytic",
                        "model": model,
                        "branch_index": branch_index,
                        "branch_id": branch_id,
                        "mu": float(mu),
                        "Lambda": lambda_value,
                        "plotted": bool(math.isfinite(float(lambda_value))),
                        "eb_thin_rod_valid": thin_rod_valid(EPSILON, float(mu)),
                        "timoshenko_omega_over_cutoff": omega_over_cutoff,
                        "timoshenko_below_cutoff": below_cutoff,
                        "diameter_to_rod1_length": ratio1,
                        "diameter_to_rod2_length": ratio2,
                        "max_diameter_to_segment_length": ratio_max,
                        "current_sorted_index": current_sorted,
                        "mac_to_previous": mac_to_previous,
                    }
                )
    return rows


def build_fem_match_rows(
    fem_results: Sequence[FemCaseResult],
    eb_lambdas: dict[str, dict[float, float]],
    eb_vectors: dict[str, dict[float, np.ndarray]],
    timo_points: dict[str, dict[float, TimoTrackedPoint]],
) -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    warnings: list[str] = []
    branch_ids = [branch_id_from_base_sorted_index(index) for index in range(1, N_BRANCHES + 1)]
    for case in fem_results:
        if not case.success:
            warnings.append(f"mu={case.mu:g}: FEM case failed: {case.message}")
            continue
        selected_solid_modes: dict[int, list[int]] = {}
        case_rows: list[dict[str, object]] = []
        for branch_index, branch_id in enumerate(branch_ids, start=1):
            key = mu_key(case.mu)
            eb_lambda = eb_lambdas.get(branch_id, {}).get(key, float("nan"))
            eb_vector = eb_vectors.get(branch_id, {}).get(key)
            timo_point = timo_points.get(branch_id, {}).get(key)
            timo_lambda = timo_point.Lambda if timo_point is not None else float("nan")
            timo_vector = timo_point.vector if timo_point is not None else None
            model_candidates: list[dict[str, object]] = []
            for model, analytic_lambda, analytic_vector in (
                ("Euler-Bernoulli", eb_lambda, eb_vector),
                ("Timoshenko", timo_lambda, timo_vector),
            ):
                if analytic_vector is None or not math.isfinite(float(analytic_lambda)):
                    continue
                best_shape: SolidShape | None = None
                best_mac = -1.0
                for shape in case.solid_shapes:
                    value = centerline_mac(shape.vector, np.asarray(analytic_vector, dtype=float))
                    if math.isfinite(value) and value > best_mac:
                        best_mac = float(value)
                        best_shape = shape
                if best_shape is not None:
                    model_candidates.append(
                        {
                            "shape_model": model,
                            "solid_shape": best_shape,
                            "mac": best_mac,
                            "analytic_lambda": analytic_lambda,
                        }
                    )
            if not model_candidates:
                warnings.append(f"mu={case.mu:g}, branch={branch_index}: no FEM/analytic shape match available.")
                continue
            best = max(model_candidates, key=lambda item: float(item["mac"]))
            shape = best["solid_shape"]
            assert isinstance(shape, SolidShape)
            selected_solid_modes.setdefault(shape.solid_mode, []).append(branch_index)
            best_mac = float(best["mac"])
            strength = mac_strength(best_mac)
            abs_diff_eb = abs(shape.Lambda - eb_lambda) if math.isfinite(float(eb_lambda)) else float("nan")
            abs_diff_timo = abs(shape.Lambda - timo_lambda) if math.isfinite(float(timo_lambda)) else float("nan")
            if math.isfinite(abs_diff_eb) and math.isfinite(abs_diff_timo):
                closer = "Timoshenko" if abs_diff_timo < abs_diff_eb else "Euler-Bernoulli"
            else:
                closer = ""
            exclusion_reasons: list[str] = []
            if strength not in {"strong", "moderate"}:
                exclusion_reasons.append("weak_mac_excluded")
            if case.negative_jacobian_warning:
                exclusion_reasons.append("mesh_quality_warning_negative_jacobian")
            plotted = not exclusion_reasons
            if strength not in {"strong", "moderate"}:
                warnings.append(
                    f"mu={case.mu:g}, branch={branch_index}: weak FEM MAC={best_mac:.3f}; FEM point excluded from plot."
                )
            ratio1, ratio2, ratio_max = diameter_to_segment_ratios(EPSILON, case.mu)
            row = (
                {
                    "row_kind": "fem",
                    "model": "planar point-joint 3D solid FEM",
                    "branch_index": branch_index,
                    "branch_id": branch_id,
                    "mu": case.mu,
                    "Lambda": shape.Lambda,
                    "plotted": plotted,
                    "exclusion_reason": ";".join(exclusion_reasons),
                    "audit_flags": ";".join(exclusion_reasons),
                    "eb_thin_rod_valid": thin_rod_valid(EPSILON, case.mu),
                    "diameter_to_rod1_length": ratio1,
                    "diameter_to_rod2_length": ratio2,
                    "max_diameter_to_segment_length": ratio_max,
                    "l1": case.l1,
                    "l2": case.l2,
                    "radius": case.radius,
                    "mesh_quality_status": case.mesh_quality_status,
                    "negative_jacobian_warning": case.negative_jacobian_warning,
                    "mesh_rod1_length_estimate": case.mesh_rod1_length_estimate,
                    "mesh_rod2_length_estimate": case.mesh_rod2_length_estimate,
                    "mesh_rod1_radius_estimate": case.mesh_rod1_radius_estimate,
                    "mesh_rod2_radius_estimate": case.mesh_rod2_radius_estimate,
                    "solid_mode": shape.solid_mode,
                    "omega": shape.omega,
                    "best_mac": best_mac,
                    "mac_strength": strength,
                    "matching_shape_model": best["shape_model"],
                    "Lambda_EB_branch": eb_lambda,
                    "Lambda_Timoshenko_branch": timo_lambda,
                    "fem_match_timoshenko_omega_over_cutoff": (
                        timoshenko_cutoff_ratio(float(timo_lambda), EPSILON)
                        if math.isfinite(float(timo_lambda))
                        else float("nan")
                    ),
                    "abs_diff_FEM_EB": abs_diff_eb,
                    "abs_diff_FEM_Timoshenko": abs_diff_timo,
                    "closer_to": closer,
                    "out_of_plane_fraction": shape.out_of_plane_fraction,
                    "filled_slice_count": shape.filled_slice_count,
                    "empty_slice_count": shape.empty_slice_count,
                    "case_dir": str(case.case_dir),
                }
            )
            case_rows.append(row)
        for solid_mode, branches in selected_solid_modes.items():
            if len(branches) > 1:
                warnings.append(f"mu={case.mu:g}: solid mode {solid_mode} matched multiple branches {branches}.")
                for row in case_rows:
                    if int(row.get("solid_mode", -1)) == int(solid_mode):
                        existing_flags = [item for item in str(row.get("audit_flags", "")).split(";") if item]
                        existing_flags.append("duplicate_match_warning")
                        row["audit_flags"] = ";".join(dict.fromkeys(existing_flags))
                        row["duplicate_match_warning"] = True
        rows.extend(case_rows)
    return rows, sorted(set(warnings))


def plot_results(rows: Sequence[dict[str, object]]) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    branch_ids = [branch_id_from_base_sorted_index(index) for index in range(1, N_BRANCHES + 1)]
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, N_BRANCHES))
    fig, ax = plt.subplots(figsize=(9.5, 5.8), constrained_layout=True)

    for branch_index, branch_id in enumerate(branch_ids, start=1):
        color = colors[branch_index - 1]
        for model, linewidth, alpha in (
            ("Euler-Bernoulli", 1.25, 0.72),
            ("Timoshenko", 2.0, 0.96),
        ):
            model_rows = [
                row
                for row in rows
                if row.get("row_kind") == "analytic"
                and row.get("model") == model
                and row.get("branch_id") == branch_id
                and bool(row.get("plotted"))
            ]
            if model == "Euler-Bernoulli":
                valid_rows = [row for row in model_rows if bool(row.get("eb_thin_rod_valid"))]
                invalid_rows = [row for row in model_rows if not bool(row.get("eb_thin_rod_valid"))]
            else:
                valid_rows = [
                    row
                    for row in model_rows
                    if str(row.get("timoshenko_below_cutoff", "")).lower() not in {"false", "0"}
                ]
                invalid_rows = [
                    row
                    for row in model_rows
                    if str(row.get("timoshenko_below_cutoff", "")).lower() in {"false", "0"}
                ]
            if valid_rows:
                ax.plot(
                    [float(row["mu"]) for row in valid_rows],
                    [float(row["Lambda"]) for row in valid_rows],
                    color=color,
                    linestyle="-",
                    linewidth=linewidth,
                    alpha=alpha,
                )
            if invalid_rows:
                ax.plot(
                    [float(row["mu"]) for row in invalid_rows],
                    [float(row["Lambda"]) for row in invalid_rows],
                    color=color,
                    linestyle="--",
                    linewidth=linewidth,
                    alpha=min(alpha, 0.55),
                )

        fem_rows = [
            row
            for row in rows
            if row.get("row_kind") == "fem"
            and row.get("branch_id") == branch_id
            and bool(row.get("plotted"))
        ]
        if fem_rows:
            ax.scatter(
                [float(row["mu"]) for row in fem_rows],
                [float(row["Lambda"]) for row in fem_rows],
                s=48,
                marker="o",
                color=color,
                edgecolor="black",
                linewidth=0.6,
                zorder=5,
            )

    ax.set_xlabel("mu")
    ax.set_ylabel("Lambda")
    ax.set_title("Coupled rods with equal thickness, beta=15 deg, epsilon=0.01, eta=0")
    ax.grid(True, which="major", color="0.88", linewidth=0.8)
    ax.set_xlim(MU_MIN - 0.02, MU_MAX + 0.02)
    handles = [
        Line2D([0], [0], color="0.15", linestyle="-", linewidth=1.25, alpha=0.72, label="Euler-Bernoulli descendants"),
        Line2D([0], [0], color="0.15", linestyle="-", linewidth=2.0, alpha=0.96, label="Timoshenko descendants"),
        Line2D([0], [0], color="0.4", linestyle="--", linewidth=1.7, label="EB dashed: thin-rod criterion violated"),
        Line2D(
            [0],
            [0],
            color="white",
            marker="o",
            markerfacecolor="0.55",
            markeredgecolor="black",
            linestyle="None",
            markersize=7,
            label="planar 3D FEM points",
        ),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=True)
    fig.savefig(OUTPUT_PNG, dpi=220)
    plt.close(fig)


def write_csv(rows: Sequence[dict[str, object]]) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def format_float(value: object, digits: int = 6) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(number):
        return "nan"
    return f"{number:.{digits}g}"


def eb_timoshenko_difference_summary(rows: Sequence[dict[str, object]]) -> dict[str, object]:
    lookup: dict[tuple[str, float, str], float] = {}
    for row in rows:
        if row.get("row_kind") != "analytic":
            continue
        try:
            lookup[(str(row["branch_id"]), mu_key(float(row["mu"])), str(row["model"]))] = float(row["Lambda"])
        except (KeyError, TypeError, ValueError):
            continue
    max_abs = -1.0
    max_rel = -1.0
    max_abs_info: tuple[str, float, float, float] | None = None
    max_rel_info: tuple[str, float, float, float] | None = None
    for branch_id, mu, model in list(lookup):
        if model != "Euler-Bernoulli":
            continue
        eb_value = lookup.get((branch_id, mu, "Euler-Bernoulli"), float("nan"))
        timo_value = lookup.get((branch_id, mu, "Timoshenko"), float("nan"))
        if not (math.isfinite(eb_value) and math.isfinite(timo_value)):
            continue
        abs_diff = abs(timo_value - eb_value)
        rel_diff = abs_diff / abs(eb_value) if abs(eb_value) > 1.0e-14 else float("nan")
        if abs_diff > max_abs:
            max_abs = abs_diff
            max_abs_info = (branch_id, mu, eb_value, timo_value)
        if math.isfinite(rel_diff) and rel_diff > max_rel:
            max_rel = rel_diff
            max_rel_info = (branch_id, mu, eb_value, timo_value)
    return {
        "max_abs_diff": max_abs,
        "max_abs_info": max_abs_info,
        "max_rel_diff": max_rel,
        "max_rel_info": max_rel_info,
    }


def timoshenko_cutoff_summary(rows: Sequence[dict[str, object]]) -> dict[str, object]:
    max_ratio = -1.0
    max_info: tuple[str, float, float] | None = None
    violations: list[tuple[str, float, float]] = []
    for row in rows:
        if row.get("row_kind") != "analytic" or row.get("model") != "Timoshenko":
            continue
        try:
            ratio = float(row.get("timoshenko_omega_over_cutoff", "nan"))
            branch_id = str(row["branch_id"])
            mu = float(row["mu"])
        except (KeyError, TypeError, ValueError):
            continue
        if not math.isfinite(ratio):
            continue
        if ratio > max_ratio:
            max_ratio = ratio
            max_info = (branch_id, mu, ratio)
        if ratio >= 1.0:
            violations.append((branch_id, mu, ratio))
    return {"max_ratio": max_ratio, "max_info": max_info, "violations": violations}


def large_discrepancy_rows(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in rows:
        if row.get("row_kind") != "fem":
            continue
        try:
            diff_eb = float(row.get("abs_diff_FEM_EB", "nan"))
            diff_timo = float(row.get("abs_diff_FEM_Timoshenko", "nan"))
        except (TypeError, ValueError):
            continue
        if diff_eb > 0.25 and diff_timo > 0.25:
            out.append(row)
    return out


def write_report(
    rows: Sequence[dict[str, object]],
    fem_results: Sequence[FemCaseResult],
    warnings: Sequence[str],
    gmsh_ccx_note: str,
) -> None:
    plotted_fem = [row for row in rows if row.get("row_kind") == "fem" and bool(row.get("plotted"))]
    fem_rows = [row for row in rows if row.get("row_kind") == "fem"]
    excluded_fem = [row for row in fem_rows if not bool(row.get("plotted"))]
    strength_counts: dict[str, int] = {"strong": 0, "moderate": 0, "weak": 0, "missing": 0}
    closer_counts: dict[str, int] = {"Euler-Bernoulli": 0, "Timoshenko": 0, "": 0}
    for row in fem_rows:
        strength_counts[str(row.get("mac_strength", "missing"))] = strength_counts.get(str(row.get("mac_strength", "missing")), 0) + 1
        if bool(row.get("plotted")):
            closer_counts[str(row.get("closer_to", ""))] = closer_counts.get(str(row.get("closer_to", "")), 0) + 1
    difference_summary = eb_timoshenko_difference_summary(rows)
    cutoff_summary = timoshenko_cutoff_summary(rows)
    large_rows = large_discrepancy_rows(rows)
    max_abs_info = difference_summary.get("max_abs_info")
    max_rel_info = difference_summary.get("max_rel_info")
    cutoff_info = cutoff_summary.get("max_info")

    lines: list[str] = [
        "# Coupled equal-thickness rods Lambda(mu): EB, Timoshenko, planar 3D FEM",
        "",
        "## Purpose",
        "",
        "Corrected diagnostic-only plot for coupled circular rods with equal thickness at beta=15 deg, epsilon=0.01, eta=0.",
        "Here `mu` changes the rod lengths, so this is not an equal-length-rods plot except at `mu=0`.",
        "The 3D points use the planar-constrained point-joint solid FEM variant because this figure compares in-plane flexural descendant branches. This is not a replacement for full 3D point-joint validation.",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g} (identical thickness)",
        f"- analytic mu grid: {len(ANALYTIC_MU_VALUES)} points from {MU_MIN:g} to {MU_MAX:g}, step {ANALYTIC_MU_VALUES[1] - ANALYTIC_MU_VALUES[0]:g}",
        f"- FEM mu grid: {', '.join(format_float(value, 3) for value in FEM_MU_VALUES)}",
        f"- descendant branches shown: {N_BRANCHES}",
        f"- FEM formulation: planar point-joint 3D solid FEM, U_z=0 on independent solid nodes plus JOINT_REF,3,5,0",
        f"- tools: {gmsh_ccx_note}",
        "",
        "## Branch Identity",
        "",
        "Analytic curves are descendant branches, not sorted roots. EB branches use `scripts/lib/analytic_branch_tracking.py` with continuation from beta=0, mu=0 to beta=15 and then along mu. Timoshenko branches use the same descendant idea in this diagnostic script: shape-MAC tracking from beta=0, mu=0 through beta=15 and then along mu.",
        "",
        "FEM points are assigned branch-by-branch using centerline MAC-like matching against the EB and Timoshenko descendant shapes at the same mu; weak matches are not plotted.",
        "",
        "## Applicability Rules",
        "",
        f"Euler--Bernoulli diameter criterion: `2*r_i/l_i <= {THIN_ROD_LIMIT:g}`.",
        "For eta=0 and epsilon=0.01, `2*r1/l1 = 4*epsilon/(1 - mu)` and `2*r2/l2 = 4*epsilon/(1 + mu)`.",
        "Rod 1 violates this EB thin-rod criterion for `mu > 0.6`; only EB curves are dashed there.",
        "Timoshenko applicability is checked separately with `Omega_c = sqrt(kappa*G*A/(rho*I))` and `Omega/Omega_c`.",
        "Timoshenko curves are not dashed by the EB diameter criterion. They stay solid while below cut-off.",
        "3D FEM points are not dashed by the EB criterion; FEM plotting is controlled by solver, mesh-quality, and MAC audit flags.",
        "",
        "## Why EB And Timoshenko Are Close For epsilon=0.01",
        "",
        "For slender rods the Timoshenko shear/rotary-inertia correction is expected to be small in the first flexural branches. At epsilon=0.01 the first six descendant branches remain well below the Timoshenko cut-off, so near-overlap of EB and Timoshenko curves is expected and is not by itself an error.",
        f"- max absolute EB-vs-Timoshenko Lambda difference: {format_float(difference_summary.get('max_abs_diff'), 6)}"
        + (
            f" at {max_abs_info[0]}, mu={format_float(max_abs_info[1], 3)}"
            if isinstance(max_abs_info, tuple)
            else ""
        ),
        f"- max relative EB-vs-Timoshenko Lambda difference: {format_float(difference_summary.get('max_rel_diff'), 6)}"
        + (
            f" at {max_rel_info[0]}, mu={format_float(max_rel_info[1], 3)}"
            if isinstance(max_rel_info, tuple)
            else ""
        ),
        f"- max Timoshenko Omega/Omega_c over plotted analytic grid: {format_float(cutoff_summary.get('max_ratio'), 6)}"
        + (
            f" at {cutoff_info[0]}, mu={format_float(cutoff_info[1], 3)}"
            if isinstance(cutoff_info, tuple)
            else ""
        ),
        f"- Timoshenko cut-off violations: {len(cutoff_summary.get('violations', []))}",
        "",
        "## FEM Case Summary",
        "",
        "| mu | status | mesh quality | nodes | elems | l1 | l2 | radius | 2r/l1 | 2r/l2 | mesh l1 | mesh l2 | mesh r1 | mesh r2 | parsed modes | max out-of-plane |",
        "|---:|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for case in fem_results:
        max_out = max((shape.out_of_plane_fraction for shape in case.solid_shapes), default=float("nan"))
        lines.append(
            "| "
            f"{case.mu:g} | {'ok' if case.success else 'failed'} | {case.mesh_quality_status} | "
            f"{case.nodes} | {case.solid_elements} | {format_float(case.l1, 4)} | {format_float(case.l2, 4)} | "
            f"{format_float(case.radius, 4)} | {format_float(case.diameter_to_rod1_length, 4)} | "
            f"{format_float(case.diameter_to_rod2_length, 4)} | {format_float(case.mesh_rod1_length_estimate, 4)} | "
            f"{format_float(case.mesh_rod2_length_estimate, 4)} | {format_float(case.mesh_rod1_radius_estimate, 4)} | "
            f"{format_float(case.mesh_rod2_radius_estimate, 4)} | {case.parsed_modes} | "
            f"{format_float(max_out, 3)} |"
        )
    if not fem_results:
        lines.append("| | FEM skipped | | | | | | | | | | | | | | |")

    lines.extend(
        [
            "",
            "## FEM Matching Summary",
            "",
            f"- strong: {strength_counts.get('strong', 0)}",
            f"- moderate: {strength_counts.get('moderate', 0)}",
            f"- weak: {strength_counts.get('weak', 0)}",
            f"- plotted FEM points: {len(plotted_fem)}",
            f"- excluded FEM branch points: {len(excluded_fem)}",
            f"- plotted points closer to EB by Lambda: {closer_counts.get('Euler-Bernoulli', 0)}",
            f"- plotted points closer to Timoshenko by Lambda: {closer_counts.get('Timoshenko', 0)}",
            "",
            "## FEM Branch Audit",
            "",
            "| mu | branch | plotted | flags | Lambda_FEM | Lambda_EB | Lambda_Timo | solid mode | MAC | strength | out-of-plane | closer |",
            "|---:|---:|:---|:---|---:|---:|---:|---:|---:|:---|---:|:---|",
        ]
    )
    for row in fem_rows:
        lines.append(
            "| "
            f"{format_float(row.get('mu'), 3)} | {row.get('branch_index')} | {bool(row.get('plotted'))} | "
            f"{row.get('audit_flags', '')} | {format_float(row.get('Lambda'), 6)} | "
            f"{format_float(row.get('Lambda_EB_branch'), 6)} | {format_float(row.get('Lambda_Timoshenko_branch'), 6)} | "
            f"{row.get('solid_mode')} | {format_float(row.get('best_mac'), 4)} | {row.get('mac_strength')} | "
            f"{format_float(row.get('out_of_plane_fraction'), 3)} | {row.get('closer_to')} |"
        )

    lines.extend(
        [
            "",
            "## Large FEM Discrepancy Audit",
            "",
        ]
    )
    if large_rows:
        lines.extend(
            [
                "| mu | branch | plotted | flags | solid mode | MAC | abs FEM-EB | abs FEM-Timo | EB invalid? | Timo Omega/Omega_c | interpretation guardrail |",
                "|---:|---:|:---|:---|---:|---:|---:|---:|:---|---:|:---|",
            ]
        )
        for row in large_rows:
            flags = str(row.get("audit_flags", ""))
            if "mesh_quality_warning" in flags:
                interpretation = "mesh issue must be resolved before physical interpretation"
            elif str(row.get("mac_strength")) == "weak":
                interpretation = "mode matching issue likely"
            elif "duplicate_match_warning" in flags:
                interpretation = "duplicate assignment may contaminate branch identity"
            else:
                interpretation = "clean enough for follow-up; could be point-joint or 3D-vs-1D deviation"
            lines.append(
                "| "
                f"{format_float(row.get('mu'), 3)} | {row.get('branch_index')} | {bool(row.get('plotted'))} | "
                f"{flags} | {row.get('solid_mode')} | {format_float(row.get('best_mac'), 4)} | "
                f"{format_float(row.get('abs_diff_FEM_EB'), 5)} | {format_float(row.get('abs_diff_FEM_Timoshenko'), 5)} | "
                f"{not bool(row.get('eb_thin_rod_valid'))} | {format_float(row.get('fem_match_timoshenko_omega_over_cutoff'), 4)} | "
                f"{interpretation} |"
            )
    else:
        lines.append("No FEM branch point exceeded both absolute-difference thresholds of 0.25.")

    lines.extend(
        [
            "",
            "## Warnings",
            "",
        ]
    )
    if warnings:
        lines.extend(f"- {item}" for item in sorted(set(warnings)))
    else:
        lines.append("- none")

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This plot is a diagnostic comparison of descendant analytic branches and planar point-joint 3D FEM points. It should be read as a check of mode identity and FEM/analytic trend consistency, not as final validation.",
            "Large FEM discrepancies are not interpreted as physical unless the mesh quality is clean and the MAC assignment is at least moderate.",
        ]
    )
    if plotted_fem:
        if closer_counts.get("Timoshenko", 0) > closer_counts.get("Euler-Bernoulli", 0):
            lines.append("Among plotted FEM points, more Lambda values are closer to Timoshenko than to Euler-Bernoulli.")
        elif closer_counts.get("Euler-Bernoulli", 0) > closer_counts.get("Timoshenko", 0):
            lines.append("Among plotted FEM points, more Lambda values are closer to Euler-Bernoulli than to Timoshenko.")
        else:
            lines.append("Among plotted FEM points, EB and Timoshenko closeness counts are tied or inconclusive.")
    else:
        lines.append("No FEM points passed the moderate/strong MAC plotting threshold.")

    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
            f"- CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
            f"- FEM generated files: `{FEM_OUTPUT_ROOT.relative_to(REPO_ROOT)}`",
        ]
    )
    OUTPUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    eb_lambdas, eb_vectors, eb_warnings = build_eb_tracking()
    timo_points, timo_warnings = track_timoshenko_descendants([float(value) for value in ANALYTIC_MU_VALUES])
    analytic_rows = build_analytic_rows(eb_lambdas, timo_points)
    fem_results, fem_warnings = run_fem_cases()
    fem_rows, match_warnings = build_fem_match_rows(fem_results, eb_lambdas, eb_vectors, timo_points)
    all_rows = [*analytic_rows, *fem_rows]
    write_csv(all_rows)
    plot_results(all_rows)
    pj = load_point_joint_module()
    gmsh_exe, gmsh_note = pj.resolve_gmsh_executable()
    ccx_exe, ccx_note = pj.resolve_ccx_executable()
    tool_note = f"Gmsh={gmsh_exe or gmsh_note}; CalculiX={ccx_exe or ccx_note}"
    write_report(all_rows, fem_results, [*eb_warnings, *timo_warnings, *fem_warnings, *match_warnings], tool_note)


if __name__ == "__main__":
    main()
