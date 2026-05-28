from __future__ import annotations

import csv
from dataclasses import dataclass
import importlib.util
import math
from pathlib import Path
import sys
from types import ModuleType
from typing import Sequence

import numpy as np
from scipy.linalg import expm


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))


BETA_DEG = 15.0
MU = 0.0
ETA = 0.0
EPSILON = 0.01
MESH_FACTORS = [1.0]
N_ANALYTIC_MODES = 10
N_SOLID_MODES = 24
L_SEGMENT = 1.0
E = 1.0
RHO = 1.0
NU = 0.3
N_SLICES_PER_ROD = 80
MAC_STRONG_THRESHOLD = 0.8
MAC_MODERATE_THRESHOLD = 0.5

OUTPUT_DIR = REPO_ROOT / "results"
AUDIT_ROOT = OUTPUT_DIR / "solid_fem_point_joint_mu0_eps0p01_audit"
OUTPUT_CSV = OUTPUT_DIR / "3d_fem_point_joint_mu0_eps0p01_audit.csv"
OUTPUT_REPORT = OUTPUT_DIR / "3d_fem_point_joint_mu0_eps0p01_audit.md"


@dataclass(frozen=True)
class AuditCase:
    case_id: str
    label: str
    beta_deg: float
    mesh_factor: float
    straight: bool
    planar_current: bool


@dataclass(frozen=True)
class StraightFrame:
    outer: tuple[float, float, float]
    joint: tuple[float, float, float]
    tangent_outer_to_joint: tuple[float, float, float]
    tangent_report: tuple[float, float, float]


@dataclass(frozen=True)
class StraightSolidShape:
    solid_mode: int
    omega: float
    lambda_fem: float
    classification: str
    vector: np.ndarray
    filled_slice_count: int
    empty_slice_count: int


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
        "point_joint_helpers_for_mu0_eps001_audit",
        REPO_ROOT / "scripts" / "analysis" / "solid_fem_coupled_equal_rods_point_joint.py",
    )


def load_single_rod_module() -> ModuleType:
    return load_module(
        "single_rod_reference_for_mu0_eps001_audit",
        REPO_ROOT / "scripts" / "analysis" / "compare_single_rod_eb_timoshenko.py",
    )


def safe_float_token(value: float) -> str:
    text = f"{float(value):.6g}"
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text.replace("-", "m").replace(".", "p")


def case_paths(pj: ModuleType, case: AuditCase) -> object:
    case_dir = AUDIT_ROOT / case.case_id / f"mesh_{safe_float_token(case.mesh_factor)}"
    beta = safe_float_token(case.beta_deg)
    eps = safe_float_token(EPSILON)
    stem = f"{case.case_id}_beta{beta}_eps{eps}"
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


def configure_point_joint_module(pj: ModuleType) -> None:
    pj.N_ANALYTIC_MODES = N_ANALYTIC_MODES
    pj.N_SOLID_MODES = N_SOLID_MODES
    pj.RUN_BETA_DEG_VALUES = [BETA_DEG]
    pj.RUN_EPSILON_VALUES = [EPSILON]
    pj.MESH_SIZE_FACTORS = MESH_FACTORS


def norm(vector: tuple[float, float, float]) -> float:
    return math.sqrt(sum(float(item) ** 2 for item in vector))


def sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def unit(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    length = norm(vector)
    if length <= 0.0:
        raise ValueError("zero vector cannot be normalized")
    return (vector[0] / length, vector[1] / length, vector[2] / length)


def straight_frames() -> dict[int, StraightFrame]:
    return {
        1: StraightFrame(
            outer=(-L_SEGMENT, 0.0, 0.0),
            joint=(0.0, 0.0, 0.0),
            tangent_outer_to_joint=(1.0, 0.0, 0.0),
            tangent_report=(1.0, 0.0, 0.0),
        ),
        2: StraightFrame(
            outer=(L_SEGMENT, 0.0, 0.0),
            joint=(0.0, 0.0, 0.0),
            tangent_outer_to_joint=(-1.0, 0.0, 0.0),
            tangent_report=(1.0, 0.0, 0.0),
        ),
    }


def write_straight_gmsh_geo(pj: ModuleType, paths: object, mesh_factor: float) -> tuple[bool, str, tuple[str, ...]]:
    radius = 2.0 * EPSILON
    mesh_size = pj.mesh_size_for_epsilon(EPSILON, float(mesh_factor))
    paths.case_dir.mkdir(parents=True, exist_ok=True)
    paths.geo.write_text(
        f"""// Diagnostic-only straight point-joint sanity mesh.
// Two collinear half-rods are separate solid volumes and are connected only by CalculiX constraints.
SetFactory("OpenCASCADE");

L = {L_SEGMENT:.17g};
R = {radius:.17g};
h = {mesh_size:.17g};

Cylinder(1) = {{-L, 0, 0, L, 0, 0, R, 2*Pi}};
Cylinder(2) = {{L, 0, 0, -L, 0, 0, R, 2*Pi}};
Physical Volume("ROD1_SOLID", 1) = {{1}};
Physical Volume("ROD2_SOLID", 2) = {{2}};

Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Mesh.ElementOrder = {pj.GMSH_ELEMENT_ORDER};
Mesh.MshFileVersion = 4.1;
""",
        encoding="utf-8",
    )
    return pj.generate_mesh_with_gmsh_cli(paths, pj.resolve_gmsh_executable()[0])


def projection_on_straight_rod(point: tuple[float, float, float], frame: StraightFrame) -> tuple[float, float]:
    rel = sub(point, frame.outer)
    s_value = dot(rel, frame.tangent_outer_to_joint)
    axis = (
        frame.outer[0] + s_value * frame.tangent_outer_to_joint[0],
        frame.outer[1] + s_value * frame.tangent_outer_to_joint[1],
        frame.outer[2] + s_value * frame.tangent_outer_to_joint[2],
    )
    radial = norm(sub(point, axis))
    return s_value, radial


def straight_rod_node_memberships(mesh: object) -> dict[int, set[int]]:
    frames = straight_frames()
    memberships: dict[int, set[int]] = {1: set(), 2: set()}
    for connectivity in mesh.solid_elements.values():
        points = [mesh.nodes[node_id] for node_id in connectivity if node_id in mesh.nodes]
        if not points:
            continue
        centroid = (
            sum(point[0] for point in points) / len(points),
            sum(point[1] for point in points) / len(points),
            sum(point[2] for point in points) / len(points),
        )
        costs = {}
        for rod_index, frame in frames.items():
            s_value, radial = projection_on_straight_rod(centroid, frame)
            outside = max(0.0, -s_value, s_value - L_SEGMENT)
            costs[rod_index] = radial + outside
        rod_index = 1 if costs[1] <= costs[2] else 2
        memberships[rod_index].update(connectivity)
    return memberships


def straight_node_sets(pj: ModuleType, mesh: object, mesh_factor: float) -> tuple[list[int], list[int], list[int], list[int]]:
    frames = straight_frames()
    memberships = straight_rod_node_memberships(mesh)
    radius = 2.0 * EPSILON
    plane_tol = max(pj.mesh_size_for_epsilon(EPSILON, float(mesh_factor)) / 50.0, 1.0e-6)
    radial_tol = 1.02 * radius + plane_tol
    out: dict[str, list[int]] = {"r1_outer": [], "r2_outer": [], "r1_inner": [], "r2_inner": []}
    for rod_index, frame in frames.items():
        for node_id in memberships[rod_index]:
            point = mesh.nodes.get(node_id)
            if point is None:
                continue
            s_value, radial = projection_on_straight_rod(point, frame)
            if radial > radial_tol:
                continue
            if abs(s_value) <= plane_tol:
                out[f"r{rod_index}_outer"].append(node_id)
            if abs(s_value - L_SEGMENT) <= plane_tol:
                out[f"r{rod_index}_inner"].append(node_id)
    return (
        sorted(out["r1_outer"]),
        sorted(out["r2_outer"]),
        sorted(out["r1_inner"]),
        sorted(out["r2_inner"]),
    )


def write_calculix_inputs(
    pj: ModuleType,
    paths: object,
    mesh: object,
    *,
    case: AuditCase,
) -> object:
    if case.straight:
        rod1_fixed, rod2_fixed, rod1_inner, rod2_inner = straight_node_sets(pj, mesh, case.mesh_factor)
    elif case.planar_current:
        rod1_fixed, rod2_fixed, rod1_inner, rod2_inner = pj.coordinate_case_node_sets(
            case.beta_deg,
            mesh,
            EPSILON,
            case.mesh_factor,
        )
    else:
        return pj.write_calculix_inputs(
            paths,
            case.beta_deg,
            EPSILON,
            mesh,
            case.mesh_factor,
            planar_constraint=False,
            planar_reference_fallback=False,
            planar_mpc_compatible_fallback=False,
        )

    joint_coupled = sorted(set(rod1_inner) | set(rod2_inner))
    all_nodes = sorted(mesh.nodes)
    planar_independent_nodes = sorted(set(all_nodes) - set(joint_coupled))
    joint_ref_node_id = max(mesh.nodes.keys(), default=0) + 1
    pj.write_calculix_mesh_include(paths, mesh)
    lines: list[str] = [
        f"** Diagnostic-only {case.label} CalculiX modal input.",
        "** Generated by audit_3d_fem_point_joint_mu0_eps0p01.py.",
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
    ]
    if case.planar_current:
        lines.extend(["*NSET, NSET=PLANAR_INDEPENDENT_SOLID_NODES", *pj.format_id_lines(planar_independent_nodes)])
    lines.extend(
        [
            "*MATERIAL, NAME=MAT",
            "*ELASTIC",
            f"{E:.17g}, {NU:.17g}",
            "*DENSITY",
            f"{RHO:.17g}",
            "*SOLID SECTION, ELSET=SOLID, MATERIAL=MAT",
            f"*RIGID BODY, NSET=JOINT_COUPLED_ALL, REF NODE={joint_ref_node_id}",
            "*BOUNDARY",
        ]
    )
    if case.planar_current:
        lines.extend(["PLANAR_INDEPENDENT_SOLID_NODES, 3, 3, 0", "JOINT_REF, 3, 5, 0"])
    lines.extend(
        [
            "ROD1_OUTER_FIXED, 1, 3, 0",
            "ROD2_OUTER_FIXED, 1, 3, 0",
            "*STEP",
            "*FREQUENCY",
            f"{N_SOLID_MODES}",
            "*NODE FILE",
            "U",
            "*END STEP",
        ]
    )
    paths.ccx_modal_inp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return pj.CcxInputSummary(
        beta_deg=case.beta_deg,
        epsilon=EPSILON,
        modal_input_path=paths.ccx_modal_inp,
        mesh_include_path=paths.ccx_mesh_inp,
        rod1_fixed_node_count=len(rod1_fixed),
        rod2_fixed_node_count=len(rod2_fixed),
        rod1_inner_node_count=len(rod1_inner),
        rod2_inner_node_count=len(rod2_inner),
        joint_coupled_node_count=len(joint_coupled),
        joint_ref_node_id=joint_ref_node_id,
        constraint_method=f"*RIGID BODY, NSET=JOINT_COUPLED_ALL, REF NODE={joint_ref_node_id}",
        node_sets_present=bool(rod1_fixed and rod2_fixed and rod1_inner and rod2_inner),
        all_solid_node_count=0,
        planar_independent_node_count=len(planar_independent_nodes) if case.planar_current else 0,
        planar_constraint=case.planar_current,
        planar_reference_fallback=case.planar_current,
        planar_mpc_compatible_fallback=case.planar_current,
    )


def parse_input_sets_and_constraints(path: Path) -> dict[str, object]:
    node_sets: dict[str, set[int]] = {}
    boundaries: list[str] = []
    rigid_lines: list[str] = []
    current_set: str | None = None
    in_boundary = False
    if not path.exists():
        return {
            "node_sets": {},
            "boundaries": [],
            "rigid_lines": [],
            "dependent_node_spc_conflict": "missing_input",
            "reference_restrictions": [],
            "has_planar_constraints": False,
        }
    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        upper = line.upper()
        if upper.startswith("*NSET"):
            current_set = ""
            in_boundary = False
            for part in line.split(","):
                if "=" in part:
                    key, value = part.split("=", 1)
                    if key.strip().upper() == "NSET":
                        current_set = value.strip()
                        node_sets.setdefault(current_set, set())
            continue
        if upper.startswith("*BOUNDARY"):
            current_set = None
            in_boundary = True
            continue
        if upper.startswith("*RIGID BODY"):
            rigid_lines.append(line)
            current_set = None
            in_boundary = False
            continue
        if line.startswith("*"):
            current_set = None
            in_boundary = False
            continue
        if current_set:
            for part in line.split(","):
                token = part.strip()
                if token:
                    try:
                        node_sets[current_set].add(int(token))
                    except ValueError:
                        pass
        if in_boundary:
            boundaries.append(line)
    dependent = node_sets.get("JOINT_COUPLED_ALL", set())
    conflicts: list[str] = []
    reference_restrictions: list[str] = []
    for line in boundaries:
        set_name = line.split(",", 1)[0].strip()
        if set_name.upper() == "JOINT_REF":
            reference_restrictions.append(line)
        nodes = node_sets.get(set_name, set())
        overlap = nodes & dependent
        if overlap:
            conflicts.append(f"{set_name}:{len(overlap)}")
    return {
        "node_sets": {key: len(value) for key, value in sorted(node_sets.items())},
        "boundaries": boundaries,
        "rigid_lines": rigid_lines,
        "dependent_node_spc_conflict": ";".join(conflicts) if conflicts else "none",
        "reference_restrictions": reference_restrictions,
        "has_planar_constraints": any("PLANAR" in line.upper() or line.upper().startswith("JOINT_REF") for line in boundaries),
    }


def warning_lines_from_files(paths: object) -> tuple[str, ...]:
    warnings: list[str] = []
    for label, path in (
        ("stdout", paths.ccx_stdout),
        ("stderr", paths.ccx_stderr),
        ("dat", paths.ccx_dat),
    ):
        if not path.exists():
            continue
        for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            stripped = line.strip()
            lower = stripped.lower()
            is_actual_error = (
                "error:" in lower
                or lower.startswith("*error")
                or lower.startswith("error in")
                or lower.startswith("error message")
            )
            is_residual_table_line = lower.startswith("error") and any(char.isdigit() for char in lower)
            if (
                "warning" in lower
                or is_actual_error
                or "negative jacobian" in lower
                or "distort" in lower
            ) and not is_residual_table_line:
                warnings.append(f"{label}: {line.strip()}")
    return tuple(warnings)


def mesh_inferred_geometry(pj: ModuleType, case: AuditCase, mesh: object) -> dict[str, float]:
    if case.straight:
        memberships = straight_rod_node_memberships(mesh)
        frames = straight_frames()
        out: dict[str, float] = {}
        for rod_index, frame in frames.items():
            s_values: list[float] = []
            radial_values: list[float] = []
            for node_id in memberships[rod_index]:
                point = mesh.nodes.get(node_id)
                if point is None:
                    continue
                s_value, radial = projection_on_straight_rod(point, frame)
                s_values.append(s_value)
                if -1.0e-7 <= s_value <= L_SEGMENT + 1.0e-7:
                    radial_values.append(radial)
            out[f"rod{rod_index}_length_estimate"] = max(s_values) - min(s_values) if s_values else float("nan")
            out[f"rod{rod_index}_radius_estimate"] = max(radial_values) if radial_values else float("nan")
        return out

    geom = pj.geometry_for_beta(case.beta_deg)
    memberships, _counts = pj.rod_node_memberships(case.beta_deg, mesh)
    out = {}
    for rod_index, rod in ((1, geom.rod1), (2, geom.rod2)):
        s_values = []
        radial_values = []
        for node_id in memberships[rod_index]:
            point = mesh.nodes.get(node_id)
            if point is None:
                continue
            s_value, radial = pj.projection_on_rod(point, rod)
            s_values.append(s_value)
            if -1.0e-7 <= s_value <= L_SEGMENT + 1.0e-7:
                radial_values.append(radial)
        out[f"rod{rod_index}_length_estimate"] = max(s_values) - min(s_values) if s_values else float("nan")
        out[f"rod{rod_index}_radius_estimate"] = max(radial_values) if radial_values else float("nan")
    return out


def normalize_vector(vector: np.ndarray) -> np.ndarray:
    scale = float(np.linalg.norm(vector))
    if scale <= 1.0e-28 or not math.isfinite(scale):
        return np.full_like(vector, np.nan, dtype=float)
    return np.asarray(vector, dtype=float) / scale


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


def straight_solid_centerline_vector(mesh: object, mode_shape: dict[int, tuple[float, float, float]]) -> tuple[np.ndarray, int, int]:
    memberships = straight_rod_node_memberships(mesh)
    frames = straight_frames()
    parts: list[np.ndarray] = []
    total_filled = 0
    total_empty = 0
    for rod_index in (1, 2):
        values = np.zeros((N_SLICES_PER_ROD, 2), dtype=float)
        counts = np.zeros(N_SLICES_PER_ROD, dtype=float)
        frame = frames[rod_index]
        for node_id in memberships[rod_index]:
            point = mesh.nodes.get(node_id)
            displacement = mode_shape.get(node_id)
            if point is None or displacement is None:
                continue
            if rod_index == 1:
                s_report = point[0] + L_SEGMENT
            else:
                s_report = point[0]
            if s_report < -1.0e-8 or s_report > L_SEGMENT + 1.0e-8:
                continue
            _s_outer, radial = projection_on_straight_rod(point, frame)
            if radial > 2.5 * EPSILON + 2.0e-3:
                continue
            slice_index = int(math.floor(max(0.0, min(0.999999999999, s_report / L_SEGMENT)) * N_SLICES_PER_ROD))
            values[slice_index, 0] += float(displacement[0])
            values[slice_index, 1] += float(displacement[1])
            counts[slice_index] += 1.0
        means, filled, empty = fill_missing_slice_means(values, counts)
        parts.append(means.reshape(-1))
        total_filled += filled
        total_empty += empty
    return normalize_vector(np.concatenate(parts)), total_filled, total_empty


def classify_straight_mode(mode_shape: dict[int, tuple[float, float, float]]) -> tuple[str, float, float, float]:
    totals = np.zeros(3, dtype=float)
    for displacement in mode_shape.values():
        totals += np.asarray(displacement, dtype=float) ** 2
    total = float(np.sum(totals))
    if total <= 0.0:
        return "missing", float("nan"), float("nan"), float("nan")
    axial = float(totals[0] / total)
    in_plane = float(totals[1] / total)
    out_plane = float(totals[2] / total)
    if in_plane >= 0.6:
        label = "in_plane_bending_like"
    elif out_plane >= 0.6:
        label = "out_of_plane_bending_like"
    elif axial >= 0.6:
        label = "axial_like"
    else:
        label = "mixed_or_unclassified"
    return label, axial, in_plane, out_plane


def straight_solid_shapes(
    pj: ModuleType,
    paths: object,
    mesh: object,
    omegas: Sequence[float],
    parse_result: object,
) -> tuple[list[StraightSolidShape], list[dict[str, object]]]:
    shapes: list[StraightSolidShape] = []
    metric_rows: list[dict[str, object]] = []
    for mode_index, omega in enumerate(omegas, start=1):
        mode_shape = parse_result.mode_shapes.get(mode_index)
        if mode_shape is None:
            continue
        vector, filled, empty = straight_solid_centerline_vector(mesh, mode_shape)
        classification, axial, in_plane, out_plane = classify_straight_mode(mode_shape)
        lambda_fem = pj.lambda_fem_from_omega(float(omega), EPSILON)
        shapes.append(
            StraightSolidShape(
                solid_mode=mode_index,
                omega=float(omega),
                lambda_fem=lambda_fem,
                classification=classification,
                vector=vector,
                filled_slice_count=filled,
                empty_slice_count=empty,
            )
        )
        metric_rows.append(
            {
                "solid_mode": mode_index,
                "Omega_FEM": float(omega),
                "Lambda_FEM": lambda_fem,
                "classification": classification,
                "axial_energy_fraction": axial,
                "in_plane_energy_fraction": in_plane,
                "out_of_plane_energy_fraction": out_plane,
                "filled_slice_count": filled,
                "empty_slice_count": empty,
            }
        )
    return shapes, metric_rows


def centerline_mac(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    denominator = float(np.dot(a, a) * np.dot(b, b))
    if denominator <= 1.0e-28 or not math.isfinite(denominator):
        return float("nan")
    value = float(np.dot(a, b))
    return float(value * value / denominator)


def mac_strength(value: float) -> str:
    if not math.isfinite(float(value)):
        return "missing"
    if value >= MAC_STRONG_THRESHOLD:
        return "strong"
    if value >= MAC_MODERATE_THRESHOLD:
        return "moderate"
    return "weak"


def single_rod_references(single: ModuleType) -> list[dict[str, object]]:
    single.N_MODES = N_ANALYTIC_MODES
    section = single.section_from_epsilon(EPSILON)
    alpha = single.fixed_fixed_eb_roots(N_ANALYTIC_MODES)
    eb_omegas = np.array([single.omega_eb(value, section) for value in alpha], dtype=float)
    timo_omegas = single.find_timoshenko_roots(section, eb_omegas, N_ANALYTIC_MODES)
    rows = []
    for mode_index in range(N_ANALYTIC_MODES):
        rows.append(
            {
                "mode": mode_index + 1,
                "alpha": float(alpha[mode_index]),
                "Lambda_EB": float(single.lambda_from_omega(eb_omegas[mode_index], section)),
                "Lambda_Timoshenko": float(single.lambda_from_omega(timo_omegas[mode_index], section)),
                "Omega_EB": float(eb_omegas[mode_index]),
                "Omega_Timoshenko": float(timo_omegas[mode_index]),
            }
        )
    return rows


def clamped_clamped_shape(alpha: float, s_norm: np.ndarray) -> np.ndarray:
    sigma = (math.cosh(alpha) - math.cos(alpha)) / (math.sinh(alpha) - math.sin(alpha))
    z = alpha * s_norm
    return np.cosh(z) - np.cos(z) - sigma * (np.sinh(z) - np.sin(z))


def single_eb_vector(alpha: float) -> np.ndarray:
    xi = (np.arange(N_SLICES_PER_ROD, dtype=float) + 0.5) / N_SLICES_PER_ROD
    s1 = 0.5 * xi
    s2 = 0.5 + 0.5 * xi
    w = np.concatenate([clamped_clamped_shape(alpha, s1), clamped_clamped_shape(alpha, s2)])
    vector = np.empty(2 * len(w), dtype=float)
    vector[0::2] = 0.0
    vector[1::2] = w
    return normalize_vector(vector)


def single_timo_vector(single: ModuleType, omega: float) -> np.ndarray:
    section = single.section_from_epsilon(EPSILON)
    state = single.timoshenko_state_matrix(float(omega), section)
    transfer = expm(state * single.L)
    _u, _singular_values, vh = np.linalg.svd(transfer[:2, 2:4])
    loads = vh[-1, :]
    y0 = np.array([0.0, 0.0, loads[0], loads[1]], dtype=float)
    xi = (np.arange(N_SLICES_PER_ROD, dtype=float) + 0.5) / N_SLICES_PER_ROD
    s_values = np.concatenate([0.5 * xi, 0.5 + 0.5 * xi])
    w_values = np.array([float((expm(state * float(s_value)) @ y0)[0]) for s_value in s_values], dtype=float)
    vector = np.empty(2 * len(w_values), dtype=float)
    vector[0::2] = 0.0
    vector[1::2] = w_values
    return normalize_vector(vector)


def straight_mac_rows(single: ModuleType, solid_shapes: Sequence[StraightSolidShape]) -> list[dict[str, object]]:
    refs = single_rod_references(single)
    rows: list[dict[str, object]] = []
    for ref in refs:
        mode = int(ref["mode"])
        vectors = {
            "EB": single_eb_vector(float(ref["alpha"])),
            "Timoshenko": single_timo_vector(single, float(ref["Omega_Timoshenko"])),
        }
        for model, vector in vectors.items():
            lambda_key = "Lambda_EB" if model == "EB" else "Lambda_Timoshenko"
            lambda_ref = float(ref[lambda_key])
            best: tuple[StraightSolidShape | None, float] = (None, float("nan"))
            candidates = [
                shape
                for shape in solid_shapes
                if shape.classification in {"in_plane_bending_like", "mixed_or_unclassified"}
            ] or list(solid_shapes)
            for shape in candidates:
                mac_value = centerline_mac(shape.vector, vector)
                if best[0] is None or (math.isfinite(mac_value) and mac_value > best[1]):
                    best = (shape, mac_value)
            shape = best[0]
            if shape is None:
                continue
            rows.append(
                {
                    "analytic_model": model,
                    "analytic_mode": mode,
                    "Lambda_analytic": lambda_ref,
                    "Lambda_EB": ref["Lambda_EB"],
                    "Lambda_Timoshenko": ref["Lambda_Timoshenko"],
                    "best_solid_mode": shape.solid_mode,
                    "best_solid_class": shape.classification,
                    "Lambda_solid": shape.lambda_fem,
                    "best_MAC": best[1],
                    "best_MAC_strength": mac_strength(best[1]),
                    "abs_diff": abs(shape.lambda_fem - lambda_ref),
                    "rel_diff": abs(shape.lambda_fem - lambda_ref) / abs(lambda_ref),
                    "abs_diff_EB": abs(shape.lambda_fem - float(ref["Lambda_EB"])),
                    "rel_diff_EB": abs(shape.lambda_fem - float(ref["Lambda_EB"])) / abs(float(ref["Lambda_EB"])),
                    "abs_diff_Timoshenko": abs(shape.lambda_fem - float(ref["Lambda_Timoshenko"])),
                    "rel_diff_Timoshenko": abs(shape.lambda_fem - float(ref["Lambda_Timoshenko"]))
                    / abs(float(ref["Lambda_Timoshenko"])),
                }
            )
    return rows


def append_case_summary_rows(
    rows: list[dict[str, object]],
    *,
    case: AuditCase,
    paths: object,
    mesh: object | None,
    mesh_summary: object | None,
    ccx_input: object | None,
    calculix_result: object | None,
    parse_result: object | None,
    gmsh_messages: Sequence[str],
) -> None:
    prefix = {"case_id": case.case_id, "case_label": case.label, "mesh_factor": case.mesh_factor}
    rows.append(
        {
            **prefix,
            "row_kind": "case_summary",
            "beta_deg": case.beta_deg,
            "epsilon": EPSILON,
            "straight_case": case.straight,
            "planar_current": case.planar_current,
            "case_dir": str(paths.case_dir.relative_to(REPO_ROOT)),
            "gmsh_messages": " | ".join(gmsh_messages),
            "calculix_success": calculix_result.success if calculix_result is not None else False,
            "calculix_message": calculix_result.message if calculix_result is not None else "not run",
            "parsed_frequency_count": len(calculix_result.parsed_omegas) if calculix_result is not None else 0,
            "parse_source": calculix_result.parse_source if calculix_result is not None else "",
            "frd_parse_success": parse_result.success if parse_result is not None else False,
            "frd_parse_message": parse_result.message if parse_result is not None else "",
        }
    )
    if mesh is not None:
        components, largest = (0, 0)
        try:
            components, largest = load_point_joint_module().connected_component_counts(mesh)
        except Exception:
            pass
        geometry = mesh_inferred_geometry(load_point_joint_module(), case, mesh)
        rows.append(
            {
                **prefix,
                "row_kind": "mesh_audit",
                "nodes": len(mesh.nodes),
                "solid_elements": len(mesh.solid_elements),
                "connected_components_before_constraints": components,
                "largest_component_elements": largest,
                "bbox_min": getattr(mesh, "bbox_min", ""),
                "bbox_max": getattr(mesh, "bbox_max", ""),
                **geometry,
            }
        )
    if mesh_summary is not None:
        rows.append(
            {
                **prefix,
                "row_kind": "mesh_warning_audit",
                "negative_jacobian_or_distortion_warning": any(
                    "negative jacobian" in item.lower() or "distort" in item.lower()
                    for item in getattr(mesh_summary, "gmsh_messages", ())
                ),
                "gmsh_warning_count": len(getattr(mesh_summary, "gmsh_messages", ())),
            }
        )
    if ccx_input is not None:
        constraint_info = parse_input_sets_and_constraints(paths.ccx_modal_inp)
        rows.append(
            {
                **prefix,
                "row_kind": "constraint_audit",
                "rod1_fixed_node_count": ccx_input.rod1_fixed_node_count,
                "rod2_fixed_node_count": ccx_input.rod2_fixed_node_count,
                "rod1_inner_node_count": ccx_input.rod1_inner_node_count,
                "rod2_inner_node_count": ccx_input.rod2_inner_node_count,
                "joint_coupled_node_count": ccx_input.joint_coupled_node_count,
                "joint_ref_node_id": ccx_input.joint_ref_node_id,
                "node_sets_present": ccx_input.node_sets_present,
                "planar_independent_node_count": getattr(ccx_input, "planar_independent_node_count", 0),
                "rigid_body_lines": " | ".join(constraint_info["rigid_lines"]),
                "boundary_lines": " | ".join(constraint_info["boundaries"]),
                "reference_restrictions": " | ".join(constraint_info["reference_restrictions"]),
                "dependent_node_spc_conflict": constraint_info["dependent_node_spc_conflict"],
                "has_planar_constraints": constraint_info["has_planar_constraints"],
            }
        )
        for name, count in constraint_info["node_sets"].items():
            rows.append({**prefix, "row_kind": "node_set_count", "node_set": name, "count": count})
    if calculix_result is not None:
        warnings = warning_lines_from_files(paths)
        rows.append(
            {
                **prefix,
                "row_kind": "calculix_warning_audit",
                "warning_count": len(warnings),
                "warnings": " | ".join(warnings[:20]),
            }
        )


def process_case(pj: ModuleType, single: ModuleType, case: AuditCase) -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    messages: list[str] = []
    paths = case_paths(pj, case)
    gmsh_exe, gmsh_source = pj.resolve_gmsh_executable()
    ccx_exe, ccx_source = pj.resolve_ccx_executable()
    if not gmsh_exe:
        raise RuntimeError("Gmsh executable was not found; set GMSH_EXE.")
    if not ccx_exe:
        raise RuntimeError("CalculiX executable was not found; set CCX_EXE.")

    if case.straight:
        mesh_ok, mesh_message, gmsh_messages = write_straight_gmsh_geo(pj, paths, case.mesh_factor)
    else:
        mesh_ok, mesh_message, _geometry_label, gmsh_messages = pj.generate_mesh(
            case.beta_deg,
            EPSILON,
            paths,
            gmsh_exe,
            case.mesh_factor,
        )
    messages.append(f"{case.case_id}: {mesh_message}; Gmsh={gmsh_source}; CalculiX={ccx_source}")
    if not mesh_ok:
        append_case_summary_rows(
            rows,
            case=case,
            paths=paths,
            mesh=None,
            mesh_summary=None,
            ccx_input=None,
            calculix_result=None,
            parse_result=None,
            gmsh_messages=gmsh_messages,
        )
        return rows, messages

    mesh = pj.read_gmsh_inp_mesh_data(paths.gmsh_inp)
    mesh_summary = pj.mesh_summary(case.beta_deg, EPSILON, mesh, gmsh_messages) if not case.straight else None
    ccx_input = write_calculix_inputs(pj, paths, mesh, case=case)
    calculix_result = pj.run_calculix(paths, ccx_exe, case.beta_deg, EPSILON)
    parse_result = pj.parse_calculix_frd_mode_shapes(paths, case.beta_deg, EPSILON)

    append_case_summary_rows(
        rows,
        case=case,
        paths=paths,
        mesh=mesh,
        mesh_summary=mesh_summary,
        ccx_input=ccx_input,
        calculix_result=calculix_result,
        parse_result=parse_result,
        gmsh_messages=gmsh_messages,
    )
    if not calculix_result.success or not calculix_result.parsed_omegas:
        return rows, messages

    if case.straight:
        solid_shapes, metric_rows = straight_solid_shapes(pj, paths, mesh, calculix_result.parsed_omegas, parse_result)
        for metric in metric_rows:
            rows.append({**case_row_prefix(case), "row_kind": "mode_metric", **metric})
        references = single_rod_references(single)
        for mode_index, omega in enumerate(calculix_result.parsed_omegas[:N_ANALYTIC_MODES], start=1):
            ref = references[mode_index - 1]
            lambda_fem = pj.lambda_fem_from_omega(float(omega), EPSILON)
            rows.append(
                {
                    **case_row_prefix(case),
                    "row_kind": "sorted_comparison",
                    "mode": mode_index,
                    "solid_sorted_mode": mode_index,
                    "Omega_FEM": float(omega),
                    "Lambda_FEM": lambda_fem,
                    "Lambda_EB": ref["Lambda_EB"],
                    "Lambda_Timoshenko": ref["Lambda_Timoshenko"],
                    "abs_diff_EB": abs(lambda_fem - float(ref["Lambda_EB"])),
                    "rel_diff_EB": abs(lambda_fem - float(ref["Lambda_EB"])) / abs(float(ref["Lambda_EB"])),
                    "abs_diff_Timoshenko": abs(lambda_fem - float(ref["Lambda_Timoshenko"])),
                    "rel_diff_Timoshenko": abs(lambda_fem - float(ref["Lambda_Timoshenko"]))
                    / abs(float(ref["Lambda_Timoshenko"])),
                    "notes": "straight sorted comparison only",
                }
            )
        for match in straight_mac_rows(single, solid_shapes):
            rows.append({**case_row_prefix(case), "row_kind": "mac_match_summary", **match})
        return rows, messages

    analytic_rows, analytic_warnings = pj.compute_analytic_rows()
    messages.extend(f"{case.case_id}: {item}" for item in analytic_warnings)
    metrics, parse_result = pj.compute_mode_metrics(case.beta_deg, EPSILON, paths, mesh, calculix_result.parsed_omegas)
    solid_shapes = pj.compute_solid_centerline_shapes(
        case.beta_deg,
        EPSILON,
        mesh,
        calculix_result.parsed_omegas,
        parse_result,
        metrics,
    )
    sorted_rows = pj.build_sorted_comparison_rows(
        analytic_rows,
        {(case.beta_deg, EPSILON): list(calculix_result.parsed_omegas)},
    )
    mac_eb_rows, mac_timo_rows, mac_warnings = pj.build_mac_matrix_rows(analytic_rows, solid_shapes)
    messages.extend(f"{case.case_id}: {item}" for item in mac_warnings)
    match_rows = pj.build_mac_match_rows(mac_eb_rows, mac_timo_rows)
    for metric in metrics:
        rows.append(
            {
                **case_row_prefix(case),
                "row_kind": "mode_metric",
                "solid_mode": metric.solid_mode,
                "Omega_FEM": metric.omega,
                "Lambda_FEM": metric.lambda_fem,
                "classification": metric.classification,
                "in_plane_energy_fraction": metric.in_plane_energy_fraction,
                "out_of_plane_energy_fraction": metric.out_of_plane_energy_fraction,
                "axial_energy_fraction": metric.axial_energy_fraction,
                "joint_motion_fraction": metric.joint_motion_fraction,
                "warping_indicator": metric.warping_indicator,
            }
        )
    for row in sorted_rows:
        rows.append(
            {
                **case_row_prefix(case),
                "row_kind": "sorted_comparison",
                "mode": row["mode"],
                "solid_sorted_mode": row["solid_sorted_mode"],
                "Omega_FEM": row["Omega_FEM"],
                "Lambda_FEM": row["Lambda_FEM"],
                "Lambda_EB": row["Lambda_EB"],
                "Lambda_Timoshenko": row["Lambda_Timoshenko"],
                "abs_diff_EB": abs(float(row["Lambda_FEM"]) - float(row["Lambda_EB"]))
                if row["Lambda_FEM"] != ""
                else "",
                "rel_diff_EB": abs(float(row["Lambda_FEM"]) - float(row["Lambda_EB"])) / abs(float(row["Lambda_EB"]))
                if row["Lambda_FEM"] != ""
                else "",
                "abs_diff_Timoshenko": abs(float(row["Lambda_FEM"]) - float(row["Lambda_Timoshenko"]))
                if row["Lambda_FEM"] != ""
                else "",
                "rel_diff_Timoshenko": abs(float(row["Lambda_FEM"]) - float(row["Lambda_Timoshenko"]))
                / abs(float(row["Lambda_Timoshenko"]))
                if row["Lambda_FEM"] != ""
                else "",
                "notes": row["notes"],
            }
        )
    for row in match_rows:
        if row.get("match_direction") != "analytic_mode_summary":
            continue
        lambda_solid = float(row["Lambda_solid"])
        lambda_eb = float(row["best_EB_Lambda"])
        lambda_timo = float(row["best_Timo_Lambda"])
        rows.append(
            {
                **case_row_prefix(case),
                "row_kind": "mac_match_summary",
                "analytic_model": row["analytic_model"],
                "analytic_mode": row["analytic_mode"],
                "Lambda_solid": lambda_solid,
                "Lambda_EB": lambda_eb,
                "Lambda_Timoshenko": lambda_timo,
                "best_solid_mode": row["best_solid_mode"],
                "best_solid_class": row["best_solid_class"],
                "best_MAC": row["best_MAC"],
                "best_MAC_strength": row["best_MAC_strength"],
                "preferred_shape_model_by_MAC": row["preferred_shape_model_by_MAC"],
                "preferred_frequency_model_after_MAC": row["preferred_frequency_model_after_MAC"],
                "abs_diff_EB": abs(lambda_solid - lambda_eb),
                "rel_diff_EB": abs(lambda_solid - lambda_eb) / abs(lambda_eb),
                "abs_diff_Timoshenko": abs(lambda_solid - lambda_timo),
                "rel_diff_Timoshenko": abs(lambda_solid - lambda_timo) / abs(lambda_timo),
            }
        )
    return rows, messages


def case_row_prefix(case: AuditCase) -> dict[str, object]:
    return {
        "case_id": case.case_id,
        "case_label": case.label,
        "mesh_factor": case.mesh_factor,
        "beta_deg": case.beta_deg,
        "epsilon": EPSILON,
        "straight_case": case.straight,
        "planar_current": case.planar_current,
    }


def write_csv(rows: Sequence[dict[str, object]]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def finite(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def fmt(value: object, digits: int = 6) -> str:
    number = finite(value)
    if not math.isfinite(number):
        return "" if value in {None, ""} else str(value)
    return f"{number:.{digits}g}"


def case_rows(rows: Sequence[dict[str, object]], case_id: str, kind: str) -> list[dict[str, object]]:
    return [row for row in rows if row.get("case_id") == case_id and row.get("row_kind") == kind]


def median_mac_diff(rows: Sequence[dict[str, object]], case_id: str) -> tuple[float, float, int]:
    subset = [
        row
        for row in case_rows(rows, case_id, "mac_match_summary")
        if str(row.get("best_MAC_strength")) in {"strong", "moderate"}
    ]
    diffs = [finite(row.get("rel_diff_EB")) for row in subset if math.isfinite(finite(row.get("rel_diff_EB")))]
    max_diff = max(diffs) if diffs else float("nan")
    med_diff = float(np.median(diffs)) if diffs else float("nan")
    return med_diff, max_diff, len(diffs)


def conclusion_from_rows(rows: Sequence[dict[str, object]]) -> str:
    planar_med, _planar_max, planar_count = median_mac_diff(rows, "beta15_planar_current")
    full_med, _full_max, full_count = median_mac_diff(rows, "beta15_full_point_joint")
    straight_med, _straight_max, straight_count = median_mac_diff(rows, "straight_point_joint")
    threshold = 0.05
    if planar_count and full_count and planar_med > threshold and full_med <= threshold:
        return "discrepancy mostly introduced by planar constraint"
    if full_count and full_med > threshold:
        if straight_count and straight_med > threshold:
            return "discrepancy already present in beta=0 straight point-joint, suggesting rigid joint constraint issue"
        return "discrepancy already present in full point-joint"
    if planar_count and planar_med > threshold:
        return "discrepancy mostly introduced by planar constraint"
    if straight_count and straight_med > threshold:
        return "discrepancy already present in beta=0 straight point-joint, suggesting rigid joint constraint issue"
    return "discrepancy remains unexplained"


def markdown_table(rows: Sequence[dict[str, object]], columns: Sequence[tuple[str, str]], limit: int | None = None) -> list[str]:
    selected = list(rows[:limit] if limit is not None else rows)
    lines = ["| " + " | ".join(label for label, _key in columns) + " |"]
    lines.append("| " + " | ".join("---" for _ in columns) + " |")
    for row in selected:
        values = []
        for _label, key in columns:
            value = row.get(key, "")
            values.append(fmt(value, 5) if isinstance(value, (int, float, np.floating)) else str(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def write_report(rows: Sequence[dict[str, object]], messages: Sequence[str]) -> None:
    conclusion = conclusion_from_rows(rows)
    planar_med, planar_max, planar_count = median_mac_diff(rows, "beta15_planar_current")
    full_med, full_max, full_count = median_mac_diff(rows, "beta15_full_point_joint")
    straight_med, straight_max, straight_count = median_mac_diff(rows, "straight_point_joint")
    straight_timo_diffs = [
        finite(row.get("rel_diff_Timoshenko"))
        for row in case_rows(rows, "straight_point_joint", "mac_match_summary")
        if math.isfinite(finite(row.get("rel_diff_Timoshenko")))
    ]
    straight_timo_max = max(straight_timo_diffs) if straight_timo_diffs else float("nan")
    lines: list[str] = [
        "# 3D FEM Point-Joint mu=0 epsilon=0.01 Audit",
        "",
        "## Purpose",
        "",
        "Diagnostic-only isolation audit for the beta=15, mu=0, eta=0, epsilon=0.01 3D point-joint discrepancy. "
        "It compares the current planar-constrained point-joint, the full point-joint without planar constraints, "
        "and a separate straight two-half-rod point-joint sanity case. No presentation plots are produced.",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg for coupled cases",
        f"- straight sanity: two collinear half-rods joined at one rigid reference node",
        f"- mu: {MU:g}",
        f"- eta: {ETA:g}",
        f"- epsilon: {EPSILON:g}",
        f"- mesh factors: {', '.join(f'{value:g}' for value in MESH_FACTORS)}",
        f"- analytic modes compared: {N_ANALYTIC_MODES}",
        f"- solid modes requested: {N_SOLID_MODES}",
        "",
        "## Required Conclusion",
        "",
        conclusion + ".",
        "",
        "Median relative EB difference over moderate/strong MAC rows:",
        "",
        "| case | rows | median rel diff | max rel diff |",
        "| --- | ---: | ---: | ---: |",
        f"| planar current | {planar_count} | {fmt(planar_med)} | {fmt(planar_max)} |",
        f"| full point-joint | {full_count} | {fmt(full_med)} | {fmt(full_max)} |",
        f"| straight point-joint | {straight_count} | {fmt(straight_med)} | {fmt(straight_max)} |",
        "",
        "Do not interpret these rows as physical validation until the discrepancy source is resolved.",
        "",
        "Interpretation guardrails:",
        "",
        "- The planar and full beta=15 point-joint cases both show moderate/strong-MAC relative differences above the 5% diagnostic threshold, so the planar constraint is not the sole source.",
        "- The straight two-half-rod point-joint sanity case has strong MAC rows and small differences "
        f"(max EB relative difference {fmt(straight_max)}, max Timoshenko relative difference {fmt(straight_timo_max)}), "
        "so a generic rigid-reference joint error is not supported by this straight-limit check.",
        "- The remaining suspect is the angled point-joint idealization or the way the 3D rigid end-face coupling maps to the 1D point-joint assumptions.",
        "",
        "## Constraint Audit",
        "",
        "| case | dependent-node SPC conflict | reference restrictions | planar constraints? | boundary lines |",
        "| --- | --- | --- | --- | --- |",
    ]
    for row in [item for item in rows if item.get("row_kind") == "constraint_audit"]:
        lines.append(
            "| "
            f"{row['case_id']} | {row.get('dependent_node_spc_conflict', '')} | "
            f"{row.get('reference_restrictions', '')} | {row.get('has_planar_constraints', '')} | "
            f"{row.get('boundary_lines', '')} |"
        )
    lines.extend(
        [
            "",
            "## Mesh Audit",
            "",
            *markdown_table(
                [item for item in rows if item.get("row_kind") == "mesh_audit"],
                [
                    ("case", "case_id"),
                    ("nodes", "nodes"),
                    ("elements", "solid_elements"),
                    ("components", "connected_components_before_constraints"),
                    ("rod1 L", "rod1_length_estimate"),
                    ("rod2 L", "rod2_length_estimate"),
                    ("rod1 r", "rod1_radius_estimate"),
                    ("rod2 r", "rod2_radius_estimate"),
                ],
            ),
            "",
            "## MAC-Matched Comparison",
            "",
            *markdown_table(
                [item for item in rows if item.get("row_kind") == "mac_match_summary"],
                [
                    ("case", "case_id"),
                    ("mode", "analytic_mode"),
                    ("solid", "best_solid_mode"),
                    ("class", "best_solid_class"),
                    ("MAC", "best_MAC"),
                    ("strength", "best_MAC_strength"),
                    ("Lambda solid", "Lambda_solid"),
                    ("Lambda EB", "Lambda_EB"),
                    ("rel EB", "rel_diff_EB"),
                    ("Lambda Timo", "Lambda_Timoshenko"),
                    ("rel Timo", "rel_diff_Timoshenko"),
                ],
            ),
            "",
            "## Sorted Comparison",
            "",
            "Sorted rows are retained only as diagnostics and are not used as mode identity.",
            "",
            *markdown_table(
                [item for item in rows if item.get("row_kind") == "sorted_comparison"],
                [
                    ("case", "case_id"),
                    ("mode", "mode"),
                    ("solid sorted", "solid_sorted_mode"),
                    ("Lambda FEM", "Lambda_FEM"),
                    ("Lambda EB", "Lambda_EB"),
                    ("rel EB", "rel_diff_EB"),
                    ("Lambda Timo", "Lambda_Timoshenko"),
                    ("rel Timo", "rel_diff_Timoshenko"),
                ],
            ),
            "",
            "## CalculiX Warnings",
            "",
        ]
    )
    warning_rows = [item for item in rows if item.get("row_kind") == "calculix_warning_audit"]
    if warning_rows:
        for row in warning_rows:
            lines.append(f"- {row['case_id']}: count={row.get('warning_count')}; {row.get('warnings', '')}")
    else:
        lines.append("- none recorded")
    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- This is diagnostic-only and does not change the existing FEM physical model or baseline outputs.",
            "- The straight sanity case uses an explicit two-half-rod geometry because the older beta=0 point-joint geometry helper creates overlapping same-side cylinders.",
            "- No single solid fixed-fixed rod epsilon=0.01 run was included in this pass.",
            "- MAC matching is centerline-displacement based and still requires visual review for article use.",
            "",
            "## Outputs",
            "",
            f"- CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
            f"- report: `{OUTPUT_REPORT.relative_to(REPO_ROOT)}`",
            f"- generated FEM files: `{AUDIT_ROOT.relative_to(REPO_ROOT)}`",
            "",
            "## Run Messages",
            "",
        ]
    )
    lines.extend(f"- {message}" for message in messages)
    OUTPUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    pj = load_point_joint_module()
    single = load_single_rod_module()
    configure_point_joint_module(pj)
    cases = [
        AuditCase(
            case_id="beta15_planar_current",
            label="beta=15 current planar-constrained point-joint",
            beta_deg=BETA_DEG,
            mesh_factor=1.0,
            straight=False,
            planar_current=True,
        ),
        AuditCase(
            case_id="beta15_full_point_joint",
            label="beta=15 full point-joint",
            beta_deg=BETA_DEG,
            mesh_factor=1.0,
            straight=False,
            planar_current=False,
        ),
        AuditCase(
            case_id="straight_point_joint",
            label="straight two-half-rod point-joint",
            beta_deg=0.0,
            mesh_factor=1.0,
            straight=True,
            planar_current=False,
        ),
    ]
    rows: list[dict[str, object]] = []
    messages: list[str] = []
    for case in cases:
        case_output, case_messages = process_case(pj, single, case)
        rows.extend(case_output)
        messages.extend(case_messages)
    write_csv(rows)
    write_report(rows, messages)
    print("3D FEM point-joint mu=0 epsilon=0.01 audit")
    print(f"Conclusion: {conclusion_from_rows(rows)}")
    print(f"Wrote {OUTPUT_CSV.relative_to(REPO_ROOT)}")
    print(f"Wrote {OUTPUT_REPORT.relative_to(REPO_ROOT)}")
    print(f"Wrote generated FEM files under {AUDIT_ROOT.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
