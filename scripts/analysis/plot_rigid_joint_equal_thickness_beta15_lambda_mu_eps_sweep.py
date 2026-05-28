from __future__ import annotations

import csv
from dataclasses import dataclass
import importlib.util
import math
from pathlib import Path
import sys
from types import ModuleType
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.optimize import linear_sum_assignment


BETA_DEG = 15.0
ETA = 0.0
EPSILON_VALUES = [0.01, 0.025, 0.05]
MU_MIN = 0.0
MU_MAX = 0.9
ANALYTIC_MU_VALUES = np.round(np.linspace(MU_MIN, MU_MAX, 91), 10)
FEM_MU_VALUES = np.round(np.linspace(MU_MIN, MU_MAX, 16), 10)
N_BRANCHES = 6
N_SOLVE_ANALYTIC = 14
N_SOLID_MODES = 30
FEM_MESH_SIZE_FACTOR = 1.0
N_SLICES_PER_ROD = 80
MAC_STRONG_THRESHOLD = 0.8
MAC_MODERATE_THRESHOLD = 0.5
THIN_ROD_LIMIT = 0.1
REUSE_EXISTING_FEM = True
LARGE_REL_ERROR_THRESHOLD = 0.15

E = 1.0
RHO = 1.0
NU = 0.3
L_SEGMENT = 1.0

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
from scripts.lib import variable_length_timoshenko as vlt  # noqa: E402


RESULTS_DIR = REPO_ROOT / "results"
FEM_OUTPUT_ROOT = RESULTS_DIR / "solid_fem_rigid_joint_equal_thickness_beta15_lambda_mu"
OUTPUT_CSV = RESULTS_DIR / "rigid_joint_equal_thickness_beta15_lambda_mu_eb_timoshenko_fem.csv"
OUTPUT_REPORT = RESULTS_DIR / "rigid_joint_equal_thickness_beta15_lambda_mu_eb_timoshenko_fem_report.md"


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
    epsilon: float
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
    epsilon: float
    mu: float
    case_dir: Path
    success: bool
    message: str
    parse_source: str
    gmsh_messages: tuple[str, ...]
    solver_warnings: tuple[str, ...]
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
    joint_ref_node_id: int
    rigid_body_line: str
    parsed_modes: int
    reused_existing_outputs: bool
    solid_shapes: tuple[SolidShape, ...]


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_point_joint_module() -> ModuleType:
    return load_module(
        "point_joint_helpers_for_rigid_lambda_mu_eps_sweep",
        REPO_ROOT / "scripts" / "analysis" / "solid_fem_coupled_equal_rods_point_joint.py",
    )


def safe_float_token(value: float) -> str:
    text = f"{float(value):.8g}"
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text.replace("-", "m").replace(".", "p")


def epsilon_output_token(epsilon: float) -> str:
    return "eps" + safe_float_token(float(epsilon))


def output_png(epsilon: float) -> Path:
    return RESULTS_DIR / (
        f"rigid_joint_equal_thickness_beta15_{epsilon_output_token(float(epsilon))}"
        "_lambda_mu_eb_timoshenko_fem.png"
    )


def mu_key(value: float) -> float:
    return round(float(value), 10)


def finite(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def fmt(value: object, digits: int = 6) -> str:
    if value is None or value == "":
        return "-"
    number = finite(value)
    if not math.isfinite(number):
        return "nan"
    return f"{number:.{digits}g}"


def rel_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def mul(value: float, vector: tuple[float, float, float]) -> tuple[float, float, float]:
    return (float(value) * vector[0], float(value) * vector[1], float(value) * vector[2])


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm(vector: tuple[float, float, float]) -> float:
    return math.sqrt(dot(vector, vector))


def unit(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    length = norm(vector)
    if length <= 0.0:
        raise ValueError("zero vector cannot be normalized")
    return (vector[0] / length, vector[1] / length, vector[2] / length)


def segment_lengths(mu: float) -> tuple[float, float]:
    return 1.0 - float(mu), 1.0 + float(mu)


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


def radius_from_epsilon(epsilon: float) -> float:
    return 2.0 * float(epsilon) * L_SEGMENT


def mesh_size_for_epsilon(epsilon: float, mesh_size_factor: float = 1.0) -> float:
    radius = radius_from_epsilon(float(epsilon))
    return min(L_SEGMENT / 40.0, radius / 2.0) * float(mesh_size_factor)


def diameter_to_segment_ratios(epsilon: float, mu: float) -> tuple[float, float, float]:
    l1, l2 = segment_lengths(float(mu))
    diameter = 2.0 * radius_from_epsilon(float(epsilon))
    ratio1 = diameter / l1 if l1 > 0.0 else float("inf")
    ratio2 = diameter / l2 if l2 > 0.0 else float("inf")
    return ratio1, ratio2, max(ratio1, ratio2)


def thin_rod_valid(epsilon: float, mu: float) -> bool:
    return diameter_to_segment_ratios(float(epsilon), float(mu))[2] <= THIN_ROD_LIMIT + 1.0e-12


def lambda_from_omega(omega: float, epsilon: float) -> float:
    if float(omega) <= 0.0 or float(epsilon) <= 0.0:
        return float("nan")
    return math.sqrt(float(omega) / float(epsilon))


def row_normalized(matrix: np.ndarray) -> np.ndarray:
    out = np.asarray(matrix, dtype=float).copy()
    norms = np.linalg.norm(out, axis=1)
    for index, row_norm in enumerate(norms):
        if row_norm > 0.0 and math.isfinite(float(row_norm)):
            out[index, :] /= row_norm
    return out


def analytic_null_vector(matrix: np.ndarray) -> tuple[np.ndarray, float, float]:
    scaled = row_normalized(np.asarray(matrix, dtype=float))
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


def normalize_vector(vector: np.ndarray) -> np.ndarray:
    out = np.asarray(vector, dtype=float)
    scale = float(np.linalg.norm(out))
    if scale <= 1.0e-28 or not math.isfinite(scale):
        return np.full_like(out, np.nan, dtype=float)
    return out / scale


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


def eb_centerline_vector(beta_deg: float, mu: float, epsilon: float, lambda_value: float) -> np.ndarray:
    matrix = assemble_clamped_coupled_matrix(
        float(lambda_value),
        math.radians(float(beta_deg)),
        float(mu),
        float(epsilon),
    )
    coeff, _smallest, _ratio = analytic_null_vector(matrix)
    a1, b1, a2, b2, p1, p2 = [float(value) for value in coeff]
    xi = slice_coordinates()
    l1, l2 = segment_lengths(float(mu))
    z1 = float(lambda_value) * l1 * xi
    theta1 = float(epsilon) * float(lambda_value) ** 2 * l1 * xi
    z2 = -float(lambda_value) * l2 * xi
    theta2 = -float(epsilon) * float(lambda_value) ** 2 * l2 * xi
    w_left = a1 * (np.cos(z1) - np.cosh(z1)) + b1 * (np.sin(z1) - np.sinh(z1))
    u_left = p1 * np.sin(theta1)
    w_right = a2 * (np.cos(z2) - np.cosh(z2)) + b2 * (np.sin(z2) - np.sinh(z2))
    u_right = p2 * np.sin(theta2)
    return normalize_vector(global_centerline_vector(float(beta_deg), float(mu), u_left, w_left, u_right, w_right))


def timo_centerline_vector(beta_deg: float, mu: float, epsilon: float, lambda_value: float) -> tuple[np.ndarray, tuple[str, ...]]:
    factors = vlt.tau_factors(float(mu), ETA)
    mode = vlt.timo_mode_coefficients(
        float(lambda_value),
        float(beta_deg),
        factors.mu,
        float(epsilon),
        factors.eta,
    )
    section1 = vlt.section_from_epsilon_tau(float(epsilon), factors.tau1)
    section2 = vlt.section_from_epsilon_tau(float(epsilon), factors.tau2)
    basis1 = vlt.timo_basis(float(lambda_value), float(epsilon), section1)
    basis2 = vlt.timo_basis(float(lambda_value), float(epsilon), section2)
    l1, l2 = vlt.segment_lengths(factors.mu)
    xi = slice_coordinates()
    fields1 = vlt.evaluate_timo_rod_fields(
        float(lambda_value),
        float(epsilon),
        section1,
        basis1,
        np.asarray(mode.coeff[0:3], dtype=float),
        l1 * xi,
    )
    fields2 = vlt.evaluate_timo_rod_fields(
        float(lambda_value),
        float(epsilon),
        section2,
        basis2,
        np.asarray(mode.coeff[3:6], dtype=float),
        -l2 * xi,
    )
    warnings = tuple(dict.fromkeys([*mode.warnings, *basis1.warnings, *basis2.warnings]))
    return (
        normalize_vector(
            global_centerline_vector(
                float(beta_deg),
                factors.mu,
                np.asarray(fields1["u"], dtype=float),
                np.asarray(fields1["w"], dtype=float),
                np.asarray(fields2["u"], dtype=float),
                np.asarray(fields2["w"], dtype=float),
            )
        ),
        warnings,
    )


def centerline_mac(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    denominator = float(np.dot(a, a) * np.dot(b, b))
    if denominator <= 1.0e-28 or not math.isfinite(denominator):
        return float("nan")
    dot_value = float(np.dot(a, b))
    return float((dot_value * dot_value) / denominator)


def mac_strength(value: object) -> str:
    number = finite(value)
    if not math.isfinite(number):
        return "missing"
    if number >= MAC_STRONG_THRESHOLD:
        return "strong"
    if number >= MAC_MODERATE_THRESHOLD:
        return "moderate"
    return "weak"


def solve_timo_states(beta_deg: float, mu: float, epsilon: float, n_roots: int) -> tuple[list[dict[str, object]], list[str]]:
    roots, warnings = vlt.timo_sorted_roots(
        float(beta_deg),
        float(mu),
        float(epsilon),
        int(n_roots),
        eta=ETA,
    )
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
        states.append({"sorted_index": int(index), "Lambda": float(root), "vector": vector})
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


def track_timoshenko_descendants(
    epsilon: float,
    mu_values: Sequence[float],
) -> tuple[dict[str, dict[float, TimoTrackedPoint]], list[str]]:
    warnings: list[str] = []
    branch_ids = [branch_id_from_base_sorted_index(index) for index in range(1, N_BRANCHES + 1)]
    beta_path = list(np.linspace(0.0, BETA_DEG, 16))
    path: list[tuple[float, float]] = [(float(beta), 0.0) for beta in beta_path]
    for mu in mu_values:
        if abs(float(mu)) <= 1.0e-12:
            continue
        path.append((BETA_DEG, float(mu)))

    first_states, first_warnings = solve_timo_states(0.0, 0.0, float(epsilon), N_SOLVE_ANALYTIC)
    warnings.extend(first_warnings)
    if len(first_states) < N_BRANCHES:
        raise RuntimeError(f"Only {len(first_states)} Timoshenko roots found at beta=0, mu=0, epsilon={epsilon:g}.")

    previous_vectors: dict[str, np.ndarray] = {}
    tracked: dict[str, dict[float, TimoTrackedPoint]] = {branch_id: {} for branch_id in branch_ids}
    for branch_id, state in zip(branch_ids, first_states[:N_BRANCHES]):
        previous_vectors[branch_id] = np.asarray(state["vector"], dtype=float)

    target_mu = {mu_key(value) for value in mu_values}
    for beta_deg, mu in path[1:]:
        states, root_warnings = solve_timo_states(beta_deg, mu, float(epsilon), N_SOLVE_ANALYTIC)
        warnings.extend(root_warnings)
        assignment = assign_previous_to_current(branch_ids, previous_vectors, states)
        if len(assignment) < len(branch_ids):
            warnings.append(f"Timoshenko assignment incomplete at beta={beta_deg:g}, mu={mu:g}, epsilon={epsilon:g}.")
        for branch_id in branch_ids:
            if branch_id not in assignment:
                continue
            state, mac_value = assignment[branch_id]
            previous_vectors[branch_id] = np.asarray(state["vector"], dtype=float)
            key = mu_key(mu)
            if abs(beta_deg - BETA_DEG) <= 1.0e-10 and key in target_mu:
                tracked[branch_id][key] = TimoTrackedPoint(
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


def timoshenko_cutoff_ratio(lambda_value: float, epsilon: float, mu: float) -> float:
    if not math.isfinite(float(lambda_value)):
        return float("nan")
    factors = vlt.tau_factors(float(mu), ETA)
    omega = vlt.project_omega(float(lambda_value), float(epsilon))
    section1 = vlt.section_from_epsilon_tau(float(epsilon), factors.tau1)
    section2 = vlt.section_from_epsilon_tau(float(epsilon), factors.tau2)
    ratios = [omega / vlt.omega_cutoff(section1), omega / vlt.omega_cutoff(section2)]
    return float(max(ratios))


def fem_case_paths(pj: ModuleType, epsilon: float, mu: float) -> object:
    eps = safe_float_token(float(epsilon))
    token_mu = safe_float_token(float(mu))
    case_dir = FEM_OUTPUT_ROOT / f"eps_{eps}" / f"mu_{token_mu}"
    stem = f"rigid_joint_beta{safe_float_token(BETA_DEG)}_eps{eps}_mu{token_mu}"
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


def write_gmsh_geo(paths: object, epsilon: float, mu: float) -> str:
    radius = radius_from_epsilon(float(epsilon))
    mesh_size = mesh_size_for_epsilon(float(epsilon), FEM_MESH_SIZE_FACTOR)
    geom = geometry_for_beta_mu(BETA_DEG, float(mu))
    cos_b = math.cos(geom.beta_rad)
    sin_b = math.sin(geom.beta_rad)
    l1, l2 = segment_lengths(float(mu))
    paths.case_dir.mkdir(parents=True, exist_ok=True)
    paths.geo.write_text(
        f"""// Diagnostic-only full rigid end-face 3D solid FEM mesh for Lambda(mu).
// beta = {BETA_DEG:.17g} deg, epsilon = {float(epsilon):.17g}, eta = {ETA:.17g}, mu = {float(mu):.17g}
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
    return "full rigid end-face variable-length separate-cylinder diagnostic"


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


def coordinate_case_node_sets(mesh: object, epsilon: float, mu: float) -> tuple[list[int], list[int], list[int], list[int]]:
    geom = geometry_for_beta_mu(BETA_DEG, float(mu))
    memberships, _element_counts = rod_node_memberships(mesh, float(mu))
    radius = radius_from_epsilon(float(epsilon))
    plane_tol = max(mesh_size_for_epsilon(float(epsilon), FEM_MESH_SIZE_FACTOR) / 50.0, 1.0e-6)
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
        estimates.append((float(max(s_values) - min(s_values)), float(max(radial_values) if radial_values else float("nan"))))
    return estimates[0][0], estimates[1][0], estimates[0][1], estimates[1][1]


def write_calculix_inputs(
    pj: ModuleType,
    paths: object,
    mesh: object,
    epsilon: float,
    mu: float,
) -> tuple[int, int, int, int, int, int, str]:
    rod1_fixed, rod2_fixed, rod1_inner, rod2_inner = coordinate_case_node_sets(mesh, float(epsilon), float(mu))
    joint_coupled = sorted(set(rod1_inner) | set(rod2_inner))
    joint_ref_node_id = max(mesh.nodes.keys(), default=0) + 1
    rigid_body_line = f"*RIGID BODY, NSET=JOINT_COUPLED_ALL, REF NODE={joint_ref_node_id}"
    pj.write_calculix_mesh_include(paths, mesh)
    lines: list[str] = [
        "** Diagnostic-only full rigid end-face point-joint 3D solid FEM input.",
        "** Separate variable-length cylinders connected by one shared rigid reference node.",
        "** No planar constraint, no patch tuning, no fused volume.",
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
        "*MATERIAL, NAME=MAT",
        "*ELASTIC",
        f"{E:.17g}, {NU:.17g}",
        "*DENSITY",
        f"{RHO:.17g}",
        "*SOLID SECTION, ELSET=SOLID, MATERIAL=MAT",
        rigid_body_line,
        "*BOUNDARY",
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
    return (
        len(rod1_fixed),
        len(rod2_fixed),
        len(rod1_inner),
        len(rod2_inner),
        len(joint_coupled),
        joint_ref_node_id,
        rigid_body_line,
    )


def fill_missing_slice_means(values: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, int, int]:
    out = np.asarray(values, dtype=float).copy()
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


def solid_centerline_vector(
    mesh: object,
    mode_shape: dict[int, tuple[float, float, float]],
    mu: float,
) -> tuple[np.ndarray, int, int, float]:
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


def collect_solver_warnings(paths: object) -> tuple[str, ...]:
    warnings: list[str] = []
    for path in (paths.ccx_stdout, paths.ccx_stderr, paths.ccx_dat):
        if not path.exists():
            continue
        for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            stripped = raw.strip()
            parts = stripped.split()
            if len(parts) == 2 and parts[0].lower() == "error":
                try:
                    float(parts[1].replace("D", "E").replace("d", "E"))
                    continue
                except ValueError:
                    pass
            lower = raw.lower()
            if any(token in lower for token in ("warning", "error", "negative", "distortion", "jacobian")):
                if stripped:
                    warnings.append(f"{path.name}: {stripped}")
    return tuple(dict.fromkeys(warnings))


def cached_outputs_available(paths: object) -> bool:
    return bool(paths.gmsh_inp.exists() and paths.ccx_dat.exists() and paths.ccx_frd.exists())


def failed_case(
    epsilon: float,
    mu: float,
    paths: object,
    message: str,
    gmsh_messages: Sequence[str],
    negative_jacobian: bool,
    mesh: object | None = None,
    mesh_estimates: tuple[float, float, float, float] | None = None,
    node_counts: tuple[int, int, int, int, int, int, str] | None = None,
    parse_source: str = "",
    reused: bool = False,
) -> FemCaseResult:
    l1, l2 = segment_lengths(float(mu))
    ratio1, ratio2, _ratio_max = diameter_to_segment_ratios(float(epsilon), float(mu))
    mesh_l1, mesh_l2, mesh_r1, mesh_r2 = mesh_estimates or (float("nan"), float("nan"), float("nan"), float("nan"))
    r1_fixed, r2_fixed, r1_inner, r2_inner, joint_coupled, joint_ref, rigid_line = node_counts or (0, 0, 0, 0, 0, 0, "")
    return FemCaseResult(
        epsilon=float(epsilon),
        mu=float(mu),
        case_dir=paths.case_dir,
        success=False,
        message=message,
        parse_source=parse_source,
        gmsh_messages=tuple(gmsh_messages),
        solver_warnings=collect_solver_warnings(paths),
        mesh_quality_status="mesh_quality_warning_negative_jacobian" if negative_jacobian else "ok",
        negative_jacobian_warning=bool(negative_jacobian),
        nodes=len(mesh.nodes) if mesh is not None else 0,
        solid_elements=len(mesh.solid_elements) if mesh is not None else 0,
        l1=l1,
        l2=l2,
        radius=radius_from_epsilon(float(epsilon)),
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
        joint_ref_node_id=joint_ref,
        rigid_body_line=rigid_line,
        parsed_modes=0,
        reused_existing_outputs=bool(reused),
        solid_shapes=(),
    )


def run_fem_case(pj: ModuleType, gmsh_exe: str, ccx_exe: str, epsilon: float, mu: float) -> FemCaseResult:
    paths = fem_case_paths(pj, float(epsilon), float(mu))
    gmsh_messages: tuple[str, ...] = ()
    negative_jacobian = False
    reused = False

    if REUSE_EXISTING_FEM and cached_outputs_available(paths):
        reused = True
        mesh = pj.read_gmsh_inp_mesh_data(paths.gmsh_inp)
        mesh_message = "reused existing mesh and CalculiX outputs"
    else:
        write_gmsh_geo(paths, float(epsilon), float(mu))
        mesh_ok, mesh_message, gmsh_messages = pj.generate_mesh_with_gmsh_cli(paths, gmsh_exe)
        negative_jacobian = has_negative_jacobian_warning(gmsh_messages)
        if not mesh_ok:
            return failed_case(
                float(epsilon),
                float(mu),
                paths,
                mesh_message,
                gmsh_messages,
                negative_jacobian,
            )
        mesh = pj.read_gmsh_inp_mesh_data(paths.gmsh_inp)

    mesh_estimates = estimate_mesh_geometry(mesh, float(mu))
    node_counts = write_calculix_inputs(pj, paths, mesh, float(epsilon), float(mu))

    if reused:
        omegas, parse_source = pj.parse_calculix_eigen_omegas(paths.ccx_dat, paths.ccx_stdout)
        ccx_message = f"Reused existing CalculiX outputs; parsed {len(omegas)} frequencies from {parse_source}."
    else:
        ccx_result = pj.run_calculix(paths, ccx_exe, BETA_DEG, float(epsilon))
        omegas = list(ccx_result.parsed_omegas)
        parse_source = ccx_result.parse_source
        ccx_message = ccx_result.message
        if not ccx_result.success:
            return failed_case(
                float(epsilon),
                float(mu),
                paths,
                ccx_message,
                gmsh_messages,
                negative_jacobian,
                mesh=mesh,
                mesh_estimates=mesh_estimates,
                node_counts=node_counts,
                parse_source=parse_source,
                reused=False,
            )

    parse_result = pj.parse_calculix_frd_mode_shapes(paths, BETA_DEG, float(epsilon))
    if not parse_result.success:
        return failed_case(
            float(epsilon),
            float(mu),
            paths,
            parse_result.message,
            gmsh_messages,
            negative_jacobian,
            mesh=mesh,
            mesh_estimates=mesh_estimates,
            node_counts=node_counts,
            parse_source=parse_source,
            reused=reused,
        )

    shapes: list[SolidShape] = []
    for mode_index, omega in enumerate(omegas[:N_SOLID_MODES], start=1):
        mode_shape = parse_result.mode_shapes.get(mode_index)
        if mode_shape is None:
            continue
        vector, filled, empty, out_fraction = solid_centerline_vector(mesh, mode_shape, float(mu))
        shapes.append(
            SolidShape(
                epsilon=float(epsilon),
                mu=float(mu),
                solid_mode=int(mode_index),
                omega=float(omega),
                Lambda=lambda_from_omega(float(omega), float(epsilon)),
                vector=vector,
                filled_slice_count=filled,
                empty_slice_count=empty,
                out_of_plane_fraction=out_fraction,
            )
        )

    l1, l2 = segment_lengths(float(mu))
    ratio1, ratio2, _ratio_max = diameter_to_segment_ratios(float(epsilon), float(mu))
    r1_fixed, r2_fixed, r1_inner, r2_inner, joint_coupled, joint_ref, rigid_line = node_counts
    return FemCaseResult(
        epsilon=float(epsilon),
        mu=float(mu),
        case_dir=paths.case_dir,
        success=True,
        message=ccx_message,
        parse_source=parse_source,
        gmsh_messages=tuple(gmsh_messages),
        solver_warnings=collect_solver_warnings(paths),
        mesh_quality_status="mesh_quality_warning_negative_jacobian" if negative_jacobian else "ok",
        negative_jacobian_warning=bool(negative_jacobian),
        nodes=len(mesh.nodes),
        solid_elements=len(mesh.solid_elements),
        l1=l1,
        l2=l2,
        radius=radius_from_epsilon(float(epsilon)),
        diameter_to_rod1_length=ratio1,
        diameter_to_rod2_length=ratio2,
        mesh_rod1_length_estimate=mesh_estimates[0],
        mesh_rod2_length_estimate=mesh_estimates[1],
        mesh_rod1_radius_estimate=mesh_estimates[2],
        mesh_rod2_radius_estimate=mesh_estimates[3],
        rod1_fixed_nodes=r1_fixed,
        rod2_fixed_nodes=r2_fixed,
        rod1_inner_nodes=r1_inner,
        rod2_inner_nodes=r2_inner,
        joint_coupled_nodes=joint_coupled,
        joint_ref_node_id=joint_ref,
        rigid_body_line=rigid_line,
        parsed_modes=len(shapes),
        reused_existing_outputs=bool(reused),
        solid_shapes=tuple(shapes),
    )


def run_fem_cases() -> tuple[list[FemCaseResult], list[str]]:
    warnings: list[str] = []
    pj = load_point_joint_module()
    pj.N_SOLID_MODES = N_SOLID_MODES
    gmsh_exe, gmsh_note = pj.resolve_gmsh_executable()
    ccx_exe, ccx_note = pj.resolve_ccx_executable()
    if gmsh_exe is None or ccx_exe is None:
        warnings.append(f"FEM skipped: gmsh={gmsh_note}; ccx={ccx_note}")
        return [], warnings

    results: list[FemCaseResult] = []
    for epsilon in EPSILON_VALUES:
        for mu in FEM_MU_VALUES:
            print(f"Running/reusing FEM epsilon={float(epsilon):g}, mu={float(mu):g}")
            case = run_fem_case(pj, gmsh_exe, ccx_exe, float(epsilon), float(mu))
            results.append(case)
            warnings.extend(f"epsilon={float(epsilon):g}, mu={float(mu):g}: {item}" for item in case.gmsh_messages)
            warnings.extend(f"epsilon={float(epsilon):g}, mu={float(mu):g}: {item}" for item in case.solver_warnings)
            if not case.success:
                warnings.append(f"epsilon={float(epsilon):g}, mu={float(mu):g}: FEM failed: {case.message}")
    return results, sorted(set(warnings))


def build_eb_tracking(
    epsilon: float,
) -> tuple[dict[str, dict[float, float]], dict[str, dict[float, np.ndarray]], list[str]]:
    branch_ids = [branch_id_from_base_sorted_index(index) for index in range(1, N_BRANCHES + 1)]
    tracking = track_mu_sweep(
        epsilon=float(epsilon),
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
            vectors[branch_id][key] = eb_centerline_vector(BETA_DEG, float(point.mu), float(epsilon), float(point.Lambda))
    return lambdas, vectors, list(tracking.warnings)


def build_analytic_rows(
    epsilon: float,
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
                    cutoff_status = ""
                else:
                    point = timo_points.get(branch_id, {}).get(key)
                    lambda_value = point.Lambda if point is not None else float("nan")
                    current_sorted = point.current_sorted_index if point is not None else ""
                    mac_to_previous = point.mac_to_previous if point is not None else ""
                    omega_over_cutoff = timoshenko_cutoff_ratio(float(lambda_value), float(epsilon), float(mu))
                    cutoff_status = cutoff_region(omega_over_cutoff)
                ratio1, ratio2, ratio_max = diameter_to_segment_ratios(float(epsilon), float(mu))
                rows.append(
                    {
                        "row_kind": "analytic",
                        "epsilon": float(epsilon),
                        "model": model,
                        "branch_index": branch_index,
                        "branch_id": branch_id,
                        "mu": float(mu),
                        "Lambda": lambda_value,
                        "plotted": bool(math.isfinite(float(lambda_value))),
                        "eb_thin_rod_valid": thin_rod_valid(float(epsilon), float(mu)),
                        "timoshenko_omega_over_cutoff": omega_over_cutoff,
                        "timoshenko_cutoff_status": cutoff_status,
                        "diameter_to_rod1_length": ratio1,
                        "diameter_to_rod2_length": ratio2,
                        "max_diameter_to_segment_length": ratio_max,
                        "current_sorted_index": current_sorted,
                        "mac_to_previous": mac_to_previous,
                    }
                )
    return rows


def best_shape_by_mac(solid_shapes: Sequence[SolidShape], analytic_vector: np.ndarray | None) -> tuple[SolidShape | None, float]:
    if analytic_vector is None:
        return None, float("nan")
    best_shape: SolidShape | None = None
    best_mac = float("nan")
    for shape in solid_shapes:
        value = centerline_mac(shape.vector, np.asarray(analytic_vector, dtype=float))
        if best_shape is None or (math.isfinite(value) and value > best_mac):
            best_shape = shape
            best_mac = float(value)
    return best_shape, best_mac


def build_fem_match_rows(
    fem_results: Sequence[FemCaseResult],
    analytic_by_epsilon: dict[float, dict[str, object]],
) -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    warnings: list[str] = []
    branch_ids = [branch_id_from_base_sorted_index(index) for index in range(1, N_BRANCHES + 1)]
    for case in fem_results:
        eps_data = analytic_by_epsilon.get(case.epsilon)
        if eps_data is None:
            continue
        if not case.success:
            warnings.append(f"epsilon={case.epsilon:g}, mu={case.mu:g}: FEM case failed: {case.message}")
            continue
        eb_lambdas = eps_data["eb_lambdas"]
        eb_vectors = eps_data["eb_vectors"]
        timo_points = eps_data["timo_points"]
        assert isinstance(eb_lambdas, dict)
        assert isinstance(eb_vectors, dict)
        assert isinstance(timo_points, dict)
        case_rows: list[dict[str, object]] = []
        selected_solid_modes: dict[int, list[int]] = {}
        for branch_index, branch_id in enumerate(branch_ids, start=1):
            key = mu_key(case.mu)
            eb_lambda = eb_lambdas.get(branch_id, {}).get(key, float("nan"))
            eb_vector = eb_vectors.get(branch_id, {}).get(key)
            timo_point = timo_points.get(branch_id, {}).get(key)
            timo_lambda = timo_point.Lambda if timo_point is not None else float("nan")
            timo_vector = timo_point.vector if timo_point is not None else None

            eb_shape, eb_mac = best_shape_by_mac(case.solid_shapes, eb_vector)
            timo_shape, timo_mac = best_shape_by_mac(case.solid_shapes, timo_vector)
            if eb_shape is None and timo_shape is None:
                warnings.append(f"epsilon={case.epsilon:g}, mu={case.mu:g}, branch={branch_index}: no MAC candidate.")
                continue

            selected_model = "Timoshenko" if finite(timo_mac) > finite(eb_mac) else "Euler-Bernoulli"
            selected_shape = timo_shape if selected_model == "Timoshenko" else eb_shape
            selected_mac = timo_mac if selected_model == "Timoshenko" else eb_mac
            if selected_shape is None:
                continue

            eb_mode = eb_shape.solid_mode if eb_shape is not None else ""
            timo_mode = timo_shape.solid_mode if timo_shape is not None else ""
            different_best_solid_modes = bool(eb_shape is not None and timo_shape is not None and eb_shape.solid_mode != timo_shape.solid_mode)
            selected_solid_modes.setdefault(selected_shape.solid_mode, []).append(branch_index)
            strength = mac_strength(selected_mac)
            rel_error_eb = abs(selected_shape.Lambda - eb_lambda) / abs(eb_lambda) if math.isfinite(finite(eb_lambda)) else float("nan")
            rel_error_timo = (
                abs(selected_shape.Lambda - timo_lambda) / abs(timo_lambda)
                if math.isfinite(finite(timo_lambda))
                else float("nan")
            )
            closer = (
                "Timoshenko"
                if math.isfinite(rel_error_timo) and math.isfinite(rel_error_eb) and rel_error_timo < rel_error_eb
                else "Euler-Bernoulli"
                if math.isfinite(rel_error_timo) and math.isfinite(rel_error_eb)
                else ""
            )

            exclusion_reasons: list[str] = []
            if strength not in {"strong", "moderate"}:
                exclusion_reasons.append("weak_mac")
            if different_best_solid_modes:
                exclusion_reasons.append("ambiguous_eb_timo_best_modes")
            if case.negative_jacobian_warning:
                exclusion_reasons.append("mesh_quality_warning_negative_jacobian")
            plotted = not exclusion_reasons
            if not plotted:
                warnings.append(
                    f"epsilon={case.epsilon:g}, mu={case.mu:g}, branch={branch_index}: "
                    f"excluded from main plot ({';'.join(exclusion_reasons)}), MAC={fmt(selected_mac, 4)}."
                )
            ratio1, ratio2, ratio_max = diameter_to_segment_ratios(case.epsilon, case.mu)
            case_rows.append(
                {
                    "row_kind": "fem_match",
                    "epsilon": case.epsilon,
                    "model": "3D FEM full rigid end-face joint",
                    "branch_index": branch_index,
                    "branch_id": branch_id,
                    "mu": case.mu,
                    "Lambda": selected_shape.Lambda,
                    "plotted": plotted,
                    "exclusion_reason": ";".join(exclusion_reasons),
                    "eb_thin_rod_valid": thin_rod_valid(case.epsilon, case.mu),
                    "diameter_to_rod1_length": ratio1,
                    "diameter_to_rod2_length": ratio2,
                    "max_diameter_to_segment_length": ratio_max,
                    "solid_mode": selected_shape.solid_mode,
                    "omega": selected_shape.omega,
                    "selected_MAC": selected_mac,
                    "mac_strength": strength,
                    "selected_shape_model_by_MAC": selected_model,
                    "best_EB_solid_mode": eb_mode,
                    "best_EB_MAC": eb_mac,
                    "best_Timoshenko_solid_mode": timo_mode,
                    "best_Timoshenko_MAC": timo_mac,
                    "different_best_solid_modes": different_best_solid_modes,
                    "Lambda_EB_branch": eb_lambda,
                    "Lambda_Timoshenko_branch": timo_lambda,
                    "rel_error_EB": rel_error_eb,
                    "rel_error_Timoshenko": rel_error_timo,
                    "closer_to": closer,
                    "timoshenko_omega_over_cutoff": timoshenko_cutoff_ratio(timo_lambda, case.epsilon, case.mu),
                    "out_of_plane_fraction": selected_shape.out_of_plane_fraction,
                    "filled_slice_count": selected_shape.filled_slice_count,
                    "empty_slice_count": selected_shape.empty_slice_count,
                    "case_dir": rel_path(case.case_dir),
                }
            )
        for solid_mode, branches in selected_solid_modes.items():
            if len(branches) <= 1:
                continue
            warnings.append(
                f"epsilon={case.epsilon:g}, mu={case.mu:g}: solid mode {solid_mode} matched multiple branches {branches}."
            )
            for row in case_rows:
                if int(row.get("solid_mode", -1)) == int(solid_mode):
                    reasons = [item for item in str(row.get("exclusion_reason", "")).split(";") if item]
                    reasons.append("duplicate_match")
                    row["exclusion_reason"] = ";".join(dict.fromkeys(reasons))
                    row["plotted"] = False
                    row["duplicate_match_warning"] = True
        rows.extend(case_rows)
    return rows, sorted(set(warnings))


def build_case_summary_rows(fem_results: Sequence[FemCaseResult]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in fem_results:
        rows.append(
            {
                "row_kind": "fem_case_summary",
                "epsilon": case.epsilon,
                "mu": case.mu,
                "success": case.success,
                "message": case.message,
                "parse_source": case.parse_source,
                "reused_existing_outputs": case.reused_existing_outputs,
                "nodes": case.nodes,
                "solid_elements": case.solid_elements,
                "parsed_modes": case.parsed_modes,
                "l1": case.l1,
                "l2": case.l2,
                "radius": case.radius,
                "diameter_to_rod1_length": case.diameter_to_rod1_length,
                "diameter_to_rod2_length": case.diameter_to_rod2_length,
                "max_diameter_to_segment_length": max(case.diameter_to_rod1_length, case.diameter_to_rod2_length),
                "mesh_rod1_length_estimate": case.mesh_rod1_length_estimate,
                "mesh_rod2_length_estimate": case.mesh_rod2_length_estimate,
                "mesh_rod1_radius_estimate": case.mesh_rod1_radius_estimate,
                "mesh_rod2_radius_estimate": case.mesh_rod2_radius_estimate,
                "rod1_fixed_nodes": case.rod1_fixed_nodes,
                "rod2_fixed_nodes": case.rod2_fixed_nodes,
                "rod1_inner_nodes": case.rod1_inner_nodes,
                "rod2_inner_nodes": case.rod2_inner_nodes,
                "joint_coupled_nodes": case.joint_coupled_nodes,
                "joint_ref_node_id": case.joint_ref_node_id,
                "rigid_body_line": case.rigid_body_line,
                "mesh_quality_status": case.mesh_quality_status,
                "negative_jacobian_warning": case.negative_jacobian_warning,
                "gmsh_warning_count": len(case.gmsh_messages),
                "solver_warning_count": len(case.solver_warnings),
                "warnings": " | ".join([*case.gmsh_messages, *case.solver_warnings]),
                "case_dir": rel_path(case.case_dir),
            }
        )
    return rows


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
        writer.writerows(rows)


def cutoff_region(ratio: object) -> str:
    number = finite(ratio)
    if not math.isfinite(number):
        return "unknown"
    if number >= 1.0:
        return "violation"
    if number >= 0.8:
        return "warning"
    return "ok"


def plot_segmented_eb(ax: plt.Axes, rows: list[dict[str, object]], color: object) -> None:
    if not rows:
        return
    current_valid: bool | None = None
    segment: list[dict[str, object]] = []
    segments: list[tuple[bool, list[dict[str, object]]]] = []
    for row in rows:
        valid = bool(row.get("eb_thin_rod_valid"))
        if current_valid is None or valid == current_valid:
            segment.append(row)
            current_valid = valid
            continue
        segments.append((bool(current_valid), segment))
        segment = [row]
        current_valid = valid
    if segment:
        segments.append((bool(current_valid), segment))
    for valid, segment_rows in segments:
        ax.plot(
            [float(row["mu"]) for row in segment_rows],
            [float(row["Lambda"]) for row in segment_rows],
            color=color,
            linestyle="-" if valid else "--",
            linewidth=1.25,
            alpha=0.72 if valid else 0.5,
        )


def plot_results(rows: Sequence[dict[str, object]], epsilon: float) -> None:
    branch_ids = [branch_id_from_base_sorted_index(index) for index in range(1, N_BRANCHES + 1)]
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, N_BRANCHES))
    fig, ax = plt.subplots(figsize=(9.6, 5.8), constrained_layout=True)
    eps_rows = [row for row in rows if abs(finite(row.get("epsilon")) - float(epsilon)) <= 1.0e-12]
    for branch_index, branch_id in enumerate(branch_ids, start=1):
        color = colors[branch_index - 1]
        eb_rows = [
            row
            for row in eps_rows
            if row.get("row_kind") == "analytic"
            and row.get("model") == "Euler-Bernoulli"
            and row.get("branch_id") == branch_id
            and bool(row.get("plotted"))
        ]
        timo_rows = [
            row
            for row in eps_rows
            if row.get("row_kind") == "analytic"
            and row.get("model") == "Timoshenko"
            and row.get("branch_id") == branch_id
            and bool(row.get("plotted"))
        ]
        plot_segmented_eb(ax, eb_rows, color)
        if timo_rows:
            ax.plot(
                [float(row["mu"]) for row in timo_rows],
                [float(row["Lambda"]) for row in timo_rows],
                color=color,
                linestyle="-",
                linewidth=2.0,
                alpha=0.96,
            )
        fem_rows = [
            row
            for row in eps_rows
            if row.get("row_kind") == "fem_match"
            and row.get("branch_id") == branch_id
            and str(row.get("plotted")).lower() in {"true", "1"}
        ]
        if fem_rows:
            ax.scatter(
                [float(row["mu"]) for row in fem_rows],
                [float(row["Lambda"]) for row in fem_rows],
                s=42,
                marker="o",
                color=color,
                edgecolor="black",
                linewidth=0.55,
                zorder=5,
            )

    ax.set_xlabel("mu")
    ax.set_ylabel("Lambda")
    ax.set_title(f"Full rigid end-face joint, beta=15 deg, epsilon={float(epsilon):g}, eta=0")
    ax.grid(True, which="major", color="0.88", linewidth=0.8)
    ax.set_xlim(MU_MIN - 0.02, MU_MAX + 0.02)
    handles = [
        Line2D([0], [0], color="0.15", linestyle="-", linewidth=1.25, alpha=0.72, label="Euler--Bernoulli descendants"),
        Line2D([0], [0], color="0.15", linestyle="-", linewidth=2.0, alpha=0.96, label="Timoshenko descendants"),
        Line2D([0], [0], color="0.35", linestyle="--", linewidth=1.25, alpha=0.55, label="EB dashed: thin-rod criterion violated"),
        Line2D(
            [0],
            [0],
            color="white",
            marker="o",
            markerfacecolor="0.55",
            markeredgecolor="black",
            linestyle="None",
            markersize=7,
            label="3D FEM full rigid end-face points",
        ),
    ]
    ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=True)
    fig.savefig(output_png(float(epsilon)), dpi=220)
    plt.close(fig)


def table(headers: Sequence[str], body: Sequence[Sequence[str]]) -> list[str]:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in body:
        lines.append("| " + " | ".join(row) + " |")
    return lines


def fem_match_rows(rows: Sequence[dict[str, object]], epsilon: float | None = None) -> list[dict[str, object]]:
    out = [row for row in rows if row.get("row_kind") == "fem_match"]
    if epsilon is not None:
        out = [row for row in out if abs(finite(row.get("epsilon")) - float(epsilon)) <= 1.0e-12]
    return out


def eligible_fem_rows(rows: Sequence[dict[str, object]], epsilon: float | None = None) -> list[dict[str, object]]:
    return [row for row in fem_match_rows(rows, epsilon) if str(row.get("plotted")).lower() in {"true", "1"}]


def mean(values: Sequence[float]) -> float:
    numbers = [float(value) for value in values if math.isfinite(float(value))]
    return float(np.mean(numbers)) if numbers else float("nan")


def write_report(rows: Sequence[dict[str, object]], fem_warnings: Sequence[str], analytic_warnings: Sequence[str]) -> None:
    lines: list[str] = [
        "# Rigid End-Face Equal-Thickness Lambda(mu) Diagnostic",
        "",
        "Diagnostic only. This workflow builds three separate `Lambda(mu)` plots for equal-thickness rods at "
        "`beta=15 deg`, using the full rigid end-face 3D FEM engineering joint. It is not patch-radius "
        "calibration and does not tune the FEM model to match either analytic reference.",
        "",
        "No article file, `paper_dorofeev_style` file, article figure, old determinant, "
        "`src/my_project/analytic/formulas.py`, old solver, existing FEM physical model, baseline result, "
        "or analytic Timoshenko helper formula is modified.",
        "",
        "## Parameters",
        "",
        f"- beta: `{BETA_DEG:g} deg`",
        f"- eta: `{ETA:g}`",
        f"- epsilon values: `{', '.join(fmt(value, 6) for value in EPSILON_VALUES)}`",
        f"- analytic mu grid: `{fmt(float(ANALYTIC_MU_VALUES[0]), 3)}:{fmt(float(ANALYTIC_MU_VALUES[1] - ANALYTIC_MU_VALUES[0]), 3)}:{fmt(float(ANALYTIC_MU_VALUES[-1]), 3)}` ({len(ANALYTIC_MU_VALUES)} values)",
        f"- FEM mu grid: `{', '.join(fmt(float(value), 3) for value in FEM_MU_VALUES)}` ({len(FEM_MU_VALUES)} values)",
        f"- descendant branches shown: `{N_BRANCHES}`",
        f"- solid modes requested per FEM case: `{N_SOLID_MODES}`",
        f"- mesh factor: `{FEM_MESH_SIZE_FACTOR:g}`",
        "",
        "## FEM Joint Model",
        "",
        "The 3D benchmark uses two separate solid cylinders with clamped outer end faces. All inner end-face "
        "nodes are coupled by CalculiX `*RIGID BODY` to a shared `JOINT_REF`. There is no planar constraint, "
        "no fused volume, and no central patch or radius tuning.",
        "",
        "## Outputs",
        "",
        f"- combined CSV: `{rel_path(OUTPUT_CSV)}`",
        f"- report: `{rel_path(OUTPUT_REPORT)}`",
        f"- generated FEM files: `{rel_path(FEM_OUTPUT_ROOT)}`",
    ]
    for epsilon in EPSILON_VALUES:
        lines.append(f"- epsilon `{epsilon:g}` plot: `{rel_path(output_png(float(epsilon)))}`")
    lines.extend(["", "## FEM Point Counts", ""])
    count_rows = []
    for epsilon in EPSILON_VALUES:
        eps_rows = fem_match_rows(rows, float(epsilon))
        eligible = eligible_fem_rows(rows, float(epsilon))
        count_rows.append(
            [
                fmt(epsilon, 6),
                str(len(eps_rows)),
                str(len(eligible)),
                str(len(eps_rows) - len(eligible)),
            ]
        )
    lines.extend(table(["epsilon", "FEM match rows", "main plotted", "excluded"], count_rows))

    lines.extend(["", "## MAC Summary", ""])
    mac_rows = []
    for epsilon in EPSILON_VALUES:
        eps_rows = fem_match_rows(rows, float(epsilon))
        mac_rows.append(
            [
                fmt(epsilon, 6),
                str(sum(1 for row in eps_rows if row.get("mac_strength") == "strong")),
                str(sum(1 for row in eps_rows if row.get("mac_strength") == "moderate")),
                str(sum(1 for row in eps_rows if row.get("mac_strength") == "weak")),
                str(sum(1 for row in eps_rows if "ambiguous" in str(row.get("exclusion_reason", "")))),
                str(sum(1 for row in eps_rows if "duplicate_match" in str(row.get("exclusion_reason", "")))),
            ]
        )
    lines.extend(table(["epsilon", "strong", "moderate", "weak", "ambiguous", "duplicate"], mac_rows))

    excluded = [row for row in fem_match_rows(rows) if str(row.get("plotted")).lower() not in {"true", "1"}]
    lines.extend(["", "## Excluded FEM Rows", ""])
    if excluded:
        lines.extend(
            table(
                ["epsilon", "mu", "branch", "solid", "MAC", "strength", "reason", "closer"],
                [
                    [
                        fmt(row.get("epsilon"), 6),
                        fmt(row.get("mu"), 3),
                        fmt(row.get("branch_index"), 0),
                        fmt(row.get("solid_mode"), 0),
                        fmt(row.get("selected_MAC"), 4),
                        str(row.get("mac_strength", "")),
                        str(row.get("exclusion_reason", "")),
                        str(row.get("closer_to", "")),
                    ]
                    for row in excluded
                ],
            )
        )
    else:
        lines.append("No FEM match rows were excluded from the main plots.")

    lines.extend(["", "## Mesh And Solver Warnings", ""])
    case_rows = [row for row in rows if row.get("row_kind") == "fem_case_summary"]
    warning_cases = [
        row
        for row in case_rows
        if not bool(row.get("success"))
        or int(finite(row.get("gmsh_warning_count"))) > 0
        or int(finite(row.get("solver_warning_count"))) > 0
        or bool(row.get("negative_jacobian_warning"))
    ]
    if warning_cases:
        lines.extend(
            table(
                ["epsilon", "mu", "success", "nodes", "elems", "gmsh warn", "solver warn", "message"],
                [
                    [
                        fmt(row.get("epsilon"), 6),
                        fmt(row.get("mu"), 3),
                        str(row.get("success", "")),
                        fmt(row.get("nodes"), 0),
                        fmt(row.get("solid_elements"), 0),
                        fmt(row.get("gmsh_warning_count"), 0),
                        fmt(row.get("solver_warning_count"), 0),
                        str(row.get("message", ""))[:120],
                    ]
                    for row in warning_cases
                ],
            )
        )
    else:
        lines.append("No failed FEM cases, negative-Jacobian markers, or captured Gmsh/CalculiX warning lines were recorded.")

    lines.extend(["", "## EB Applicability", ""])
    eb_rows = []
    for epsilon in EPSILON_VALUES:
        valid_mu = [
            float(mu)
            for mu in ANALYTIC_MU_VALUES
            if thin_rod_valid(float(epsilon), float(mu))
        ]
        ratio0 = diameter_to_segment_ratios(float(epsilon), 0.0)[2]
        ratio_max = diameter_to_segment_ratios(float(epsilon), float(ANALYTIC_MU_VALUES[-1]))[2]
        eb_rows.append(
            [
                fmt(epsilon, 6),
                fmt(ratio0, 4),
                fmt(ratio_max, 4),
                f"{fmt(min(valid_mu), 3)}..{fmt(max(valid_mu), 3)}" if valid_mu else "none",
                str(len(valid_mu)),
            ]
        )
    lines.extend(table(["epsilon", "max 2r/l at mu=0", "max 2r/l at mu=0.9", "valid mu range", "valid analytic grid pts"], eb_rows))
    lines.append("")
    lines.append("Only EB curves are dashed where this criterion fails. Timoshenko curves are not dashed by the EB thin-rod rule.")

    lines.extend(["", "## Timoshenko Cut-Off", ""])
    cutoff_rows = []
    for epsilon in EPSILON_VALUES:
        timo_rows = [
            row
            for row in rows
            if row.get("row_kind") == "analytic"
            and row.get("model") == "Timoshenko"
            and abs(finite(row.get("epsilon")) - float(epsilon)) <= 1.0e-12
        ]
        max_row = max(timo_rows, key=lambda row: finite(row.get("timoshenko_omega_over_cutoff")), default={})
        max_ratio = finite(max_row.get("timoshenko_omega_over_cutoff"))
        warning_count = sum(1 for row in timo_rows if cutoff_region(row.get("timoshenko_omega_over_cutoff")) == "warning")
        violation_count = sum(1 for row in timo_rows if cutoff_region(row.get("timoshenko_omega_over_cutoff")) == "violation")
        cutoff_rows.append(
            [
                fmt(epsilon, 6),
                fmt(max_ratio, 5),
                fmt(max_row.get("mu"), 3),
                fmt(max_row.get("branch_index"), 0),
                str(warning_count),
                str(violation_count),
                cutoff_region(max_ratio),
            ]
        )
    lines.extend(table(["epsilon", "max Omega/Omega_c", "mu", "branch", "warning rows", "violation rows", "status"], cutoff_rows))

    lines.extend(["", "## EB vs Timoshenko Closeness", ""])
    trend_rows = []
    for epsilon in EPSILON_VALUES:
        eligible = eligible_fem_rows(rows, float(epsilon))
        mean_eb = mean([finite(row.get("rel_error_EB")) for row in eligible])
        mean_timo = mean([finite(row.get("rel_error_Timoshenko")) for row in eligible])
        ratio = mean_timo / mean_eb if mean_eb > 0.0 and math.isfinite(mean_eb) else float("nan")
        closer_eb = sum(1 for row in eligible if row.get("closer_to") == "Euler-Bernoulli")
        closer_timo = sum(1 for row in eligible if row.get("closer_to") == "Timoshenko")
        trend_rows.append(
            [
                fmt(epsilon, 6),
                str(len(eligible)),
                fmt(mean_eb, 5),
                fmt(mean_timo, 5),
                fmt(ratio, 5),
                str(closer_eb),
                str(closer_timo),
            ]
        )
    lines.extend(table(["epsilon", "eligible", "mean rel EB", "mean rel Timo", "Timo/EB", "closer EB", "closer Timo"], trend_rows))

    lines.extend(["", "## Suspicious Large-Gap Rows", ""])
    suspicious = [
        row
        for row in eligible_fem_rows(rows)
        if min(finite(row.get("rel_error_EB")), finite(row.get("rel_error_Timoshenko"))) > LARGE_REL_ERROR_THRESHOLD
    ]
    if suspicious:
        lines.extend(
            table(
                ["epsilon", "mu", "branch", "solid", "MAC", "strength", "rel EB", "rel Timo", "closer", "mesh warnings"],
                [
                    [
                        fmt(row.get("epsilon"), 6),
                        fmt(row.get("mu"), 3),
                        fmt(row.get("branch_index"), 0),
                        fmt(row.get("solid_mode"), 0),
                        fmt(row.get("selected_MAC"), 4),
                        str(row.get("mac_strength", "")),
                        fmt(row.get("rel_error_EB"), 4),
                        fmt(row.get("rel_error_Timoshenko"), 4),
                        str(row.get("closer_to", "")),
                        str(row.get("exclusion_reason", "")) or "none",
                    ]
                    for row in suspicious
                ],
            )
        )
        lines.append("")
        lines.append("These rows are not hidden; they should be treated as follow-up audit candidates if they affect the interpretation.")
    else:
        lines.append("No eligible FEM rows exceeded the large-gap screen.")

    lines.extend(["", "## Interpretation", ""])
    for epsilon in EPSILON_VALUES:
        eligible = eligible_fem_rows(rows, float(epsilon))
        closer_eb = sum(1 for row in eligible if row.get("closer_to") == "Euler-Bernoulli")
        closer_timo = sum(1 for row in eligible if row.get("closer_to") == "Timoshenko")
        mean_eb = mean([finite(row.get("rel_error_EB")) for row in eligible])
        mean_timo = mean([finite(row.get("rel_error_Timoshenko")) for row in eligible])
        lines.append(
            f"- epsilon `{epsilon:g}`: eligible FEM rows `{len(eligible)}`; closer counts EB/Timoshenko "
            f"`{closer_eb}/{closer_timo}`; mean relative errors EB/Timoshenko "
            f"`{fmt(mean_eb, 5)}/{fmt(mean_timo, 5)}`."
        )
    lines.append("")
    lines.append(
        "The intended conclusion is comparative, not exact 1D/3D equivalence. If thicker cases show lower "
        "Timoshenko/EB error ratios and more Timoshenko-closer FEM points, this supports the diagnostic trend; "
        "mesh convergence and visual mode review remain separate article-grade requirements."
    )

    if analytic_warnings:
        lines.extend(["", "## Analytic Warnings", ""])
        for warning in analytic_warnings[:80]:
            lines.append(f"- {warning}")
        if len(analytic_warnings) > 80:
            lines.append(f"- ... {len(analytic_warnings) - 80} more analytic warnings")
    if fem_warnings:
        lines.extend(["", "## FEM Warnings", ""])
        for warning in fem_warnings[:120]:
            lines.append(f"- {warning}")
        if len(fem_warnings) > 120:
            lines.append(f"- ... {len(fem_warnings) - 120} more FEM warnings")

    OUTPUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, object]] = []
    analytic_by_epsilon: dict[float, dict[str, object]] = {}
    analytic_warnings: list[str] = []

    for epsilon in EPSILON_VALUES:
        print(f"Tracking analytic descendants epsilon={float(epsilon):g}")
        eb_lambdas, eb_vectors, eb_warnings = build_eb_tracking(float(epsilon))
        timo_points, timo_warnings = track_timoshenko_descendants(float(epsilon), ANALYTIC_MU_VALUES)
        analytic_warnings.extend(f"epsilon={float(epsilon):g}: {item}" for item in eb_warnings)
        analytic_warnings.extend(f"epsilon={float(epsilon):g}: {item}" for item in timo_warnings)
        analytic_by_epsilon[float(epsilon)] = {
            "eb_lambdas": eb_lambdas,
            "eb_vectors": eb_vectors,
            "timo_points": timo_points,
        }
        all_rows.extend(build_analytic_rows(float(epsilon), eb_lambdas, timo_points))

    fem_results, fem_warnings = run_fem_cases()
    all_rows.extend(build_case_summary_rows(fem_results))
    fem_match, match_warnings = build_fem_match_rows(fem_results, analytic_by_epsilon)
    fem_warnings.extend(match_warnings)
    all_rows.extend(fem_match)

    write_csv(all_rows)
    for epsilon in EPSILON_VALUES:
        plot_results(all_rows, float(epsilon))
    write_report(all_rows, sorted(set(fem_warnings)), sorted(set(analytic_warnings)))

    print(f"Wrote {rel_path(OUTPUT_CSV)}")
    for epsilon in EPSILON_VALUES:
        print(f"Wrote {rel_path(output_png(float(epsilon)))}")
    print(f"Wrote {rel_path(OUTPUT_REPORT)}")
    print(f"Wrote generated FEM files under {rel_path(FEM_OUTPUT_ROOT)}")


if __name__ == "__main__":
    main()
