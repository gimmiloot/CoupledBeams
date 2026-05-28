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


EPSILON = 0.01
MU = 0.0
ETA = 0.0
BETA_VALUES = [0.0, 15.0, 45.0, 90.0]
MESH_FACTOR = 1.0
N_ANALYTIC_MODES = 10
N_SUMMARY_MODES = 6
N_SOLID_MODES = 24
L_SEGMENT = 1.0
E = 1.0
RHO = 1.0
NU = 0.3
N_SLICES_PER_ROD = 80
MAC_STRONG_THRESHOLD = 0.8
MAC_MODERATE_THRESHOLD = 0.5
SANITY_TOL = 0.03
TARGET_TOL = 0.03
LARGE_ERROR_TOL = 0.05

OUTPUT_DIR = REPO_ROOT / "results"
AUDIT_ROOT = OUTPUT_DIR / "solid_fem_joint_coupling_models_audit"
OUTPUT_CSV = OUTPUT_DIR / "3d_fem_joint_coupling_models_audit.csv"
OUTPUT_REPORT = OUTPUT_DIR / "3d_fem_joint_coupling_models_audit.md"


@dataclass(frozen=True)
class JointModel:
    model_id: str
    label: str
    patch_radius_fraction: float | None


@dataclass(frozen=True)
class AuditCase:
    beta_deg: float
    model: JointModel

    @property
    def case_id(self) -> str:
        return f"beta_{safe_float_token(self.beta_deg)}_{self.model.model_id}"

    @property
    def straight(self) -> bool:
        return abs(float(self.beta_deg)) <= 1.0e-12


@dataclass(frozen=True)
class StraightFrame:
    outer: tuple[float, float, float]
    tangent_outer_to_joint: tuple[float, float, float]


@dataclass(frozen=True)
class StraightSolidShape:
    solid_mode: int
    omega: float
    lambda_fem: float
    classification: str
    vector: np.ndarray
    filled_slice_count: int
    empty_slice_count: int


JOINT_MODELS = [
    JointModel("rigid_end_faces", "all inner end-face nodes tied to JOINT_REF", None),
    JointModel("center_patch_0p25", "central patch radius 0.25*r tied to JOINT_REF", 0.25),
    JointModel("center_patch_0p5", "central patch radius 0.5*r tied to JOINT_REF", 0.5),
]


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
        "point_joint_helpers_for_joint_coupling_audit",
        REPO_ROOT / "scripts" / "analysis" / "solid_fem_coupled_equal_rods_point_joint.py",
    )


def load_single_rod_module() -> ModuleType:
    return load_module(
        "single_rod_reference_for_joint_coupling_audit",
        REPO_ROOT / "scripts" / "analysis" / "compare_single_rod_eb_timoshenko.py",
    )


def configure_point_joint_module(pj: ModuleType) -> None:
    pj.N_ANALYTIC_MODES = N_ANALYTIC_MODES
    pj.N_SOLID_MODES = N_SOLID_MODES
    pj.MESH_SIZE_FACTORS = [MESH_FACTOR]
    pj.RUN_EPSILON_VALUES = [EPSILON]


def safe_float_token(value: float) -> str:
    text = f"{float(value):.6g}"
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text.replace("-", "m").replace(".", "p")


def case_paths(pj: ModuleType, case: AuditCase) -> object:
    beta_token = safe_float_token(case.beta_deg)
    eps_token = safe_float_token(EPSILON)
    case_dir = AUDIT_ROOT / case.model.model_id / f"beta_{beta_token}" / f"mesh_{safe_float_token(MESH_FACTOR)}"
    stem = f"{case.model.model_id}_beta{beta_token}_eps{eps_token}"
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


def norm(vector: tuple[float, float, float]) -> float:
    return math.sqrt(sum(float(item) ** 2 for item in vector))


def sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def straight_frames() -> dict[int, StraightFrame]:
    return {
        1: StraightFrame(outer=(-L_SEGMENT, 0.0, 0.0), tangent_outer_to_joint=(1.0, 0.0, 0.0)),
        2: StraightFrame(outer=(L_SEGMENT, 0.0, 0.0), tangent_outer_to_joint=(-1.0, 0.0, 0.0)),
    }


def projection_on_straight_rod(point: tuple[float, float, float], frame: StraightFrame) -> tuple[float, float]:
    rel = sub(point, frame.outer)
    s_value = dot(rel, frame.tangent_outer_to_joint)
    axis = (
        frame.outer[0] + s_value * frame.tangent_outer_to_joint[0],
        frame.outer[1] + s_value * frame.tangent_outer_to_joint[1],
        frame.outer[2] + s_value * frame.tangent_outer_to_joint[2],
    )
    return s_value, norm(sub(point, axis))


def write_straight_gmsh_geo(pj: ModuleType, paths: object) -> tuple[bool, str, tuple[str, ...]]:
    radius = 2.0 * EPSILON
    mesh_size = pj.mesh_size_for_epsilon(EPSILON, MESH_FACTOR)
    paths.case_dir.mkdir(parents=True, exist_ok=True)
    paths.geo.write_text(
        f"""// Diagnostic-only straight two-half-rod point-joint mesh.
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
    gmsh_exe = pj.resolve_gmsh_executable()[0]
    if gmsh_exe is None:
        return False, "Gmsh executable not found.", ()
    return pj.generate_mesh_with_gmsh_cli(paths, gmsh_exe)


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
            costs[rod_index] = radial + max(0.0, -s_value, s_value - L_SEGMENT)
        memberships[1 if costs[1] <= costs[2] else 2].update(connectivity)
    return memberships


def straight_node_sets(pj: ModuleType, mesh: object) -> tuple[list[int], list[int], list[int], list[int]]:
    frames = straight_frames()
    memberships = straight_rod_node_memberships(mesh)
    radius = 2.0 * EPSILON
    plane_tol = max(pj.mesh_size_for_epsilon(EPSILON, MESH_FACTOR) / 50.0, 1.0e-6)
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


def coupled_node_sets(pj: ModuleType, case: AuditCase, mesh: object) -> tuple[list[int], list[int], list[int], list[int]]:
    if case.straight:
        return straight_node_sets(pj, mesh)
    return pj.coordinate_case_node_sets(case.beta_deg, mesh, EPSILON, MESH_FACTOR)


def node_radial_distance(pj: ModuleType, case: AuditCase, mesh: object, rod_index: int, node_id: int) -> float:
    point = mesh.nodes[node_id]
    if case.straight:
        return projection_on_straight_rod(point, straight_frames()[rod_index])[1]
    geom = pj.geometry_for_beta(case.beta_deg)
    rod = geom.rod1 if rod_index == 1 else geom.rod2
    return pj.projection_on_rod(point, rod)[1]


def patch_nodes(
    pj: ModuleType,
    case: AuditCase,
    mesh: object,
    rod_index: int,
    full_nodes: Sequence[int],
    patch_fraction: float,
) -> tuple[list[int], bool, float]:
    patch_radius = float(patch_fraction) * 2.0 * EPSILON
    scored = [(node_radial_distance(pj, case, mesh, rod_index, node_id), int(node_id)) for node_id in full_nodes]
    selected = [node_id for radial, node_id in scored if radial <= patch_radius + 1.0e-12]
    if selected:
        return sorted(selected), False, patch_radius
    if not scored:
        return [], True, patch_radius
    scored.sort()
    return [scored[0][1]], True, patch_radius


def write_calculix_inputs(
    pj: ModuleType,
    paths: object,
    mesh: object,
    case: AuditCase,
) -> tuple[object, dict[str, object]]:
    r1_fixed, r2_fixed, r1_inner_full, r2_inner_full = coupled_node_sets(pj, case, mesh)
    fallback = False
    patch_radius = ""
    if case.model.patch_radius_fraction is None:
        r1_inner = r1_inner_full
        r2_inner = r2_inner_full
    else:
        r1_inner, fallback1, radius1 = patch_nodes(
            pj, case, mesh, 1, r1_inner_full, case.model.patch_radius_fraction
        )
        r2_inner, fallback2, radius2 = patch_nodes(
            pj, case, mesh, 2, r2_inner_full, case.model.patch_radius_fraction
        )
        fallback = fallback1 or fallback2
        patch_radius = max(radius1, radius2)
    joint_coupled = sorted(set(r1_inner) | set(r2_inner))
    joint_ref_node_id = max(mesh.nodes.keys(), default=0) + 1
    pj.write_calculix_mesh_include(paths, mesh)
    lines = [
        f"** Diagnostic-only joint coupling model: {case.model.model_id}.",
        f"** patch_radius_fraction={case.model.patch_radius_fraction}",
        f"*INCLUDE, INPUT={paths.ccx_mesh_inp.name}",
        "*NODE",
        f"{joint_ref_node_id}, 0.0, 0.0, 0.0",
        "*NSET, NSET=JOINT_REF",
        f"{joint_ref_node_id}",
        "*NSET, NSET=ROD1_OUTER_FIXED",
        *pj.format_id_lines(r1_fixed),
        "*NSET, NSET=ROD2_OUTER_FIXED",
        *pj.format_id_lines(r2_fixed),
        "*NSET, NSET=ROD1_INNER_FULL",
        *pj.format_id_lines(r1_inner_full),
        "*NSET, NSET=ROD2_INNER_FULL",
        *pj.format_id_lines(r2_inner_full),
        "*NSET, NSET=ROD1_INNER_COUPLED",
        *pj.format_id_lines(r1_inner),
        "*NSET, NSET=ROD2_INNER_COUPLED",
        *pj.format_id_lines(r2_inner),
        "*NSET, NSET=JOINT_COUPLED_ALL",
        *pj.format_id_lines(joint_coupled),
        "*MATERIAL, NAME=MAT",
        "*ELASTIC",
        f"{E:.17g}, {NU:.17g}",
        "*DENSITY",
        f"{RHO:.17g}",
        "*SOLID SECTION, ELSET=SOLID, MATERIAL=MAT",
        f"*RIGID BODY, NSET=JOINT_COUPLED_ALL, REF NODE={joint_ref_node_id}",
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
    summary = pj.CcxInputSummary(
        beta_deg=case.beta_deg,
        epsilon=EPSILON,
        modal_input_path=paths.ccx_modal_inp,
        mesh_include_path=paths.ccx_mesh_inp,
        rod1_fixed_node_count=len(r1_fixed),
        rod2_fixed_node_count=len(r2_fixed),
        rod1_inner_node_count=len(r1_inner),
        rod2_inner_node_count=len(r2_inner),
        joint_coupled_node_count=len(joint_coupled),
        joint_ref_node_id=joint_ref_node_id,
        constraint_method=f"*RIGID BODY, NSET=JOINT_COUPLED_ALL, REF NODE={joint_ref_node_id}",
        node_sets_present=bool(r1_fixed and r2_fixed and r1_inner and r2_inner),
    )
    metadata = {
        "rod1_inner_full_count": len(r1_inner_full),
        "rod2_inner_full_count": len(r2_inner_full),
        "patch_radius": patch_radius,
        "patch_fallback_used": fallback,
    }
    return summary, metadata


def parse_input_constraints(path: Path) -> dict[str, object]:
    node_sets: dict[str, set[int]] = {}
    boundaries: list[str] = []
    rigid_lines: list[str] = []
    equation_lines: list[str] = []
    current_set: str | None = None
    in_boundary = False
    in_equation = False
    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        upper = line.upper()
        if upper.startswith("*NSET"):
            current_set = ""
            in_boundary = False
            in_equation = False
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
            in_equation = False
            continue
        if upper.startswith("*RIGID BODY"):
            rigid_lines.append(line)
            current_set = None
            in_boundary = False
            in_equation = False
            continue
        if upper.startswith("*EQUATION"):
            equation_lines.append(line)
            current_set = None
            in_boundary = False
            in_equation = True
            continue
        if line.startswith("*"):
            current_set = None
            in_boundary = False
            in_equation = False
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
        if in_equation:
            equation_lines.append(line)
    dependent = node_sets.get("JOINT_COUPLED_ALL", set())
    conflicts: list[str] = []
    reference_restrictions: list[str] = []
    for line in boundaries:
        set_name = line.split(",", 1)[0].strip()
        if set_name.upper() == "JOINT_REF":
            reference_restrictions.append(line)
        overlap = node_sets.get(set_name, set()) & dependent
        if overlap:
            conflicts.append(f"{set_name}:{len(overlap)}")
    return {
        "node_sets": {key: len(value) for key, value in sorted(node_sets.items())},
        "boundaries": boundaries,
        "rigid_lines": rigid_lines,
        "equation_lines": equation_lines,
        "reference_restrictions": reference_restrictions,
        "dependent_node_spc_conflict": ";".join(conflicts) if conflicts else "none",
    }


def warning_lines_from_files(paths: object) -> tuple[str, ...]:
    warnings: list[str] = []
    for label, path in (("stdout", paths.ccx_stdout), ("stderr", paths.ccx_stderr), ("dat", paths.ccx_dat)):
        if not path.exists():
            continue
        for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            line = raw_line.strip()
            lower = line.lower()
            residual_table = lower.startswith("error") and any(char.isdigit() for char in lower)
            actual_error = "error:" in lower or lower.startswith("*error") or lower.startswith("error in")
            if ("warning" in lower or actual_error or "negative jacobian" in lower or "distort" in lower) and not residual_table:
                warnings.append(f"{label}: {line}")
    return tuple(warnings)


def mesh_inferred_geometry(pj: ModuleType, case: AuditCase, mesh: object) -> dict[str, float]:
    if case.straight:
        memberships = straight_rod_node_memberships(mesh)
        frames = straight_frames()
        out: dict[str, float] = {}
        for rod_index, frame in frames.items():
            s_values = []
            radial_values = []
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
    parts = []
    total_filled = 0
    total_empty = 0
    for rod_index in (1, 2):
        values = np.zeros((N_SLICES_PER_ROD, 2), dtype=float)
        counts = np.zeros(N_SLICES_PER_ROD, dtype=float)
        for node_id in memberships[rod_index]:
            point = mesh.nodes.get(node_id)
            displacement = mode_shape.get(node_id)
            if point is None or displacement is None:
                continue
            s_report = point[0] + L_SEGMENT if rod_index == 1 else point[0]
            if s_report < -1.0e-8 or s_report > L_SEGMENT + 1.0e-8:
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
    mesh: object,
    omegas: Sequence[float],
    parse_result: object,
) -> tuple[list[StraightSolidShape], list[dict[str, object]]]:
    shapes = []
    rows = []
    for mode_index, omega in enumerate(omegas, start=1):
        shape = parse_result.mode_shapes.get(mode_index)
        if shape is None:
            continue
        vector, filled, empty = straight_solid_centerline_vector(mesh, shape)
        label, axial, in_plane, out_plane = classify_straight_mode(shape)
        lambda_fem = pj.lambda_fem_from_omega(float(omega), EPSILON)
        shapes.append(
            StraightSolidShape(
                solid_mode=mode_index,
                omega=float(omega),
                lambda_fem=lambda_fem,
                classification=label,
                vector=vector,
                filled_slice_count=filled,
                empty_slice_count=empty,
            )
        )
        rows.append(
            {
                "solid_mode": mode_index,
                "Omega_FEM": float(omega),
                "Lambda_FEM": lambda_fem,
                "classification": label,
                "axial_energy_fraction": axial,
                "in_plane_energy_fraction": in_plane,
                "out_of_plane_energy_fraction": out_plane,
                "filled_slice_count": filled,
                "empty_slice_count": empty,
            }
        )
    return shapes, rows


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
    return [
        {
            "mode": index + 1,
            "alpha": float(alpha[index]),
            "Lambda_EB": float(single.lambda_from_omega(eb_omegas[index], section)),
            "Lambda_Timoshenko": float(single.lambda_from_omega(timo_omegas[index], section)),
            "Omega_EB": float(eb_omegas[index]),
            "Omega_Timoshenko": float(timo_omegas[index]),
        }
        for index in range(N_ANALYTIC_MODES)
    ]


def clamped_clamped_shape(alpha: float, s_norm: np.ndarray) -> np.ndarray:
    sigma = (math.cosh(alpha) - math.cos(alpha)) / (math.sinh(alpha) - math.sin(alpha))
    z = alpha * s_norm
    return np.cosh(z) - np.cos(z) - sigma * (np.sinh(z) - np.sin(z))


def single_eb_vector(alpha: float) -> np.ndarray:
    xi = (np.arange(N_SLICES_PER_ROD, dtype=float) + 0.5) / N_SLICES_PER_ROD
    s_values = np.concatenate([0.5 * xi, 0.5 + 0.5 * xi])
    w_values = clamped_clamped_shape(alpha, s_values)
    vector = np.empty(2 * len(w_values), dtype=float)
    vector[0::2] = 0.0
    vector[1::2] = w_values
    return normalize_vector(vector)


def single_timo_vector(single: ModuleType, omega: float) -> np.ndarray:
    section = single.section_from_epsilon(EPSILON)
    state = single.timoshenko_state_matrix(float(omega), section)
    transfer = expm(state * single.L)
    _u, _singular_values, vh = np.linalg.svd(transfer[:2, 2:4])
    y0 = np.array([0.0, 0.0, vh[-1, 0], vh[-1, 1]], dtype=float)
    xi = (np.arange(N_SLICES_PER_ROD, dtype=float) + 0.5) / N_SLICES_PER_ROD
    s_values = np.concatenate([0.5 * xi, 0.5 + 0.5 * xi])
    w_values = np.array([float((expm(state * float(s_value)) @ y0)[0]) for s_value in s_values], dtype=float)
    vector = np.empty(2 * len(w_values), dtype=float)
    vector[0::2] = 0.0
    vector[1::2] = w_values
    return normalize_vector(vector)


def straight_mac_rows(single: ModuleType, solid_shapes: Sequence[StraightSolidShape]) -> list[dict[str, object]]:
    refs = single_rod_references(single)
    rows = []
    for ref in refs:
        vectors = {
            "EB": single_eb_vector(float(ref["alpha"])),
            "Timoshenko": single_timo_vector(single, float(ref["Omega_Timoshenko"])),
        }
        candidates = [
            shape for shape in solid_shapes if shape.classification in {"in_plane_bending_like", "mixed_or_unclassified"}
        ] or list(solid_shapes)
        for model, vector in vectors.items():
            lambda_ref = float(ref["Lambda_EB"] if model == "EB" else ref["Lambda_Timoshenko"])
            best_shape = None
            best_mac = float("nan")
            for shape in candidates:
                mac = centerline_mac(shape.vector, vector)
                if best_shape is None or (math.isfinite(mac) and mac > best_mac):
                    best_shape = shape
                    best_mac = mac
            if best_shape is None:
                continue
            rows.append(
                {
                    "analytic_model": model,
                    "analytic_mode": int(ref["mode"]),
                    "Lambda_analytic": lambda_ref,
                    "Lambda_EB": float(ref["Lambda_EB"]),
                    "Lambda_Timoshenko": float(ref["Lambda_Timoshenko"]),
                    "best_solid_mode": best_shape.solid_mode,
                    "best_solid_class": best_shape.classification,
                    "Lambda_solid": best_shape.lambda_fem,
                    "best_MAC": best_mac,
                    "best_MAC_strength": mac_strength(best_mac),
                    "abs_diff_EB": abs(best_shape.lambda_fem - float(ref["Lambda_EB"])),
                    "rel_diff_EB": abs(best_shape.lambda_fem - float(ref["Lambda_EB"])) / abs(float(ref["Lambda_EB"])),
                    "abs_diff_Timoshenko": abs(best_shape.lambda_fem - float(ref["Lambda_Timoshenko"])),
                    "rel_diff_Timoshenko": abs(best_shape.lambda_fem - float(ref["Lambda_Timoshenko"]))
                    / abs(float(ref["Lambda_Timoshenko"])),
                }
            )
    return rows


def case_prefix(case: AuditCase) -> dict[str, object]:
    return {
        "case_id": case.case_id,
        "model_id": case.model.model_id,
        "model_label": case.model.label,
        "beta_deg": case.beta_deg,
        "epsilon": EPSILON,
        "mu": MU,
        "eta": ETA,
        "mesh_factor": MESH_FACTOR,
        "straight_case": case.straight,
    }


def append_audit_rows(
    rows: list[dict[str, object]],
    *,
    case: AuditCase,
    paths: object,
    mesh: object | None,
    mesh_summary: object | None,
    ccx_input: object | None,
    ccx_metadata: dict[str, object],
    calculix_result: object | None,
    parse_result: object | None,
    gmsh_messages: Sequence[str],
) -> None:
    prefix = case_prefix(case)
    rows.append(
        {
            **prefix,
            "row_kind": "case_summary",
            "case_dir": str(paths.case_dir.relative_to(REPO_ROOT)),
            "calculix_success": calculix_result.success if calculix_result is not None else False,
            "calculix_message": calculix_result.message if calculix_result is not None else "not run",
            "parsed_frequency_count": len(calculix_result.parsed_omegas) if calculix_result is not None else 0,
            "parse_source": calculix_result.parse_source if calculix_result is not None else "",
            "frd_parse_success": parse_result.success if parse_result is not None else False,
            "frd_parse_message": parse_result.message if parse_result is not None else "",
            "gmsh_messages": " | ".join(gmsh_messages),
        }
    )
    if mesh is not None:
        components, largest = load_point_joint_module().connected_component_counts(mesh)
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
                **mesh_inferred_geometry(load_point_joint_module(), case, mesh),
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
        constraint = parse_input_constraints(paths.ccx_modal_inp)
        rows.append(
            {
                **prefix,
                "row_kind": "constraint_audit",
                "rod1_fixed_node_count": ccx_input.rod1_fixed_node_count,
                "rod2_fixed_node_count": ccx_input.rod2_fixed_node_count,
                "rod1_inner_coupled_node_count": ccx_input.rod1_inner_node_count,
                "rod2_inner_coupled_node_count": ccx_input.rod2_inner_node_count,
                "rod1_inner_full_node_count": ccx_metadata.get("rod1_inner_full_count", ""),
                "rod2_inner_full_node_count": ccx_metadata.get("rod2_inner_full_count", ""),
                "joint_coupled_node_count": ccx_input.joint_coupled_node_count,
                "patch_radius": ccx_metadata.get("patch_radius", ""),
                "patch_fallback_used": ccx_metadata.get("patch_fallback_used", ""),
                "joint_ref_node_id": ccx_input.joint_ref_node_id,
                "rigid_body_lines": " | ".join(constraint["rigid_lines"]),
                "equation_lines": " | ".join(constraint["equation_lines"]),
                "boundary_lines": " | ".join(constraint["boundaries"]),
                "reference_restrictions": " | ".join(constraint["reference_restrictions"]),
                "dependent_node_spc_conflict": constraint["dependent_node_spc_conflict"],
            }
        )
        for name, count in constraint["node_sets"].items():
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
    if gmsh_exe is None:
        raise RuntimeError("Gmsh executable not found. Set GMSH_EXE.")
    if ccx_exe is None:
        raise RuntimeError("CalculiX executable not found. Set CCX_EXE.")
    if case.straight:
        mesh_ok, mesh_message, gmsh_messages = write_straight_gmsh_geo(pj, paths)
    else:
        mesh_ok, mesh_message, _geometry_label, gmsh_messages = pj.generate_mesh(
            case.beta_deg, EPSILON, paths, gmsh_exe, MESH_FACTOR
        )
    messages.append(f"{case.case_id}: {mesh_message}; Gmsh={gmsh_source}; CalculiX={ccx_source}")
    if not mesh_ok:
        append_audit_rows(
            rows,
            case=case,
            paths=paths,
            mesh=None,
            mesh_summary=None,
            ccx_input=None,
            ccx_metadata={},
            calculix_result=None,
            parse_result=None,
            gmsh_messages=gmsh_messages,
        )
        return rows, messages
    mesh = pj.read_gmsh_inp_mesh_data(paths.gmsh_inp)
    mesh_summary = None if case.straight else pj.mesh_summary(case.beta_deg, EPSILON, mesh, gmsh_messages)
    ccx_input, ccx_metadata = write_calculix_inputs(pj, paths, mesh, case)
    calculix_result = pj.run_calculix(paths, ccx_exe, case.beta_deg, EPSILON)
    parse_result = pj.parse_calculix_frd_mode_shapes(paths, case.beta_deg, EPSILON)
    append_audit_rows(
        rows,
        case=case,
        paths=paths,
        mesh=mesh,
        mesh_summary=mesh_summary,
        ccx_input=ccx_input,
        ccx_metadata=ccx_metadata,
        calculix_result=calculix_result,
        parse_result=parse_result,
        gmsh_messages=gmsh_messages,
    )
    if not calculix_result.success or not calculix_result.parsed_omegas:
        return rows, messages
    if case.straight:
        solid_shapes, metric_rows = straight_solid_shapes(pj, mesh, calculix_result.parsed_omegas, parse_result)
        for row in metric_rows:
            rows.append({**case_prefix(case), "row_kind": "mode_metric", **row})
        refs = single_rod_references(single)
        for mode_index, omega in enumerate(calculix_result.parsed_omegas[:N_ANALYTIC_MODES], start=1):
            ref = refs[mode_index - 1]
            lambda_fem = pj.lambda_fem_from_omega(float(omega), EPSILON)
            rows.append(
                {
                    **case_prefix(case),
                    "row_kind": "sorted_comparison",
                    "mode": mode_index,
                    "solid_sorted_mode": mode_index,
                    "Omega_FEM": float(omega),
                    "Lambda_FEM": lambda_fem,
                    "Lambda_EB": ref["Lambda_EB"],
                    "Lambda_Timoshenko": ref["Lambda_Timoshenko"],
                    "rel_diff_EB": abs(lambda_fem - float(ref["Lambda_EB"])) / abs(float(ref["Lambda_EB"])),
                    "rel_diff_Timoshenko": abs(lambda_fem - float(ref["Lambda_Timoshenko"]))
                    / abs(float(ref["Lambda_Timoshenko"])),
                    "notes": "straight sorted comparison only",
                }
            )
        for row in straight_mac_rows(single, solid_shapes):
            rows.append({**case_prefix(case), "row_kind": "mac_match_summary", **row})
        return rows, messages

    pj.RUN_BETA_DEG_VALUES = [case.beta_deg]
    pj.RUN_EPSILON_VALUES = [EPSILON]
    analytic_rows, analytic_warnings = pj.compute_analytic_rows()
    messages.extend(f"{case.case_id}: {item}" for item in analytic_warnings)
    metrics, parse_result = pj.compute_mode_metrics(case.beta_deg, EPSILON, paths, mesh, calculix_result.parsed_omegas)
    solid_shapes = pj.compute_solid_centerline_shapes(
        case.beta_deg, EPSILON, mesh, calculix_result.parsed_omegas, parse_result, metrics
    )
    sorted_rows = pj.build_sorted_comparison_rows(analytic_rows, {(case.beta_deg, EPSILON): list(calculix_result.parsed_omegas)})
    mac_eb_rows, mac_timo_rows, mac_warnings = pj.build_mac_matrix_rows(analytic_rows, solid_shapes)
    messages.extend(f"{case.case_id}: {item}" for item in mac_warnings)
    match_rows = pj.build_mac_match_rows(mac_eb_rows, mac_timo_rows)
    for metric in metrics:
        rows.append(
            {
                **case_prefix(case),
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
        lambda_fem = row["Lambda_FEM"]
        rows.append(
            {
                **case_prefix(case),
                "row_kind": "sorted_comparison",
                "mode": row["mode"],
                "solid_sorted_mode": row["solid_sorted_mode"],
                "Omega_FEM": row["Omega_FEM"],
                "Lambda_FEM": lambda_fem,
                "Lambda_EB": row["Lambda_EB"],
                "Lambda_Timoshenko": row["Lambda_Timoshenko"],
                "rel_diff_EB": abs(float(lambda_fem) - float(row["Lambda_EB"])) / abs(float(row["Lambda_EB"]))
                if lambda_fem != ""
                else "",
                "rel_diff_Timoshenko": abs(float(lambda_fem) - float(row["Lambda_Timoshenko"]))
                / abs(float(row["Lambda_Timoshenko"]))
                if lambda_fem != ""
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
                **case_prefix(case),
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
                "rel_diff_EB": abs(lambda_solid - lambda_eb) / abs(lambda_eb),
                "rel_diff_Timoshenko": abs(lambda_solid - lambda_timo) / abs(lambda_timo),
            }
        )
    return rows, messages


def finite(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def summary_source_rows(rows: Sequence[dict[str, object]], model_id: str, beta: float) -> list[dict[str, object]]:
    subset = [
        row
        for row in rows
        if row.get("row_kind") == "mac_match_summary"
        and row.get("model_id") == model_id
        and abs(finite(row.get("beta_deg")) - float(beta)) <= 1.0e-12
        and int(float(row.get("analytic_mode", 999))) <= N_SUMMARY_MODES
    ]
    if abs(float(beta)) <= 1.0e-12:
        subset = [row for row in subset if row.get("analytic_model") == "EB"]
    return subset


def summarize_case(rows: Sequence[dict[str, object]], model_id: str, beta: float) -> dict[str, object]:
    subset = summary_source_rows(rows, model_id, beta)
    eligible = [row for row in subset if row.get("best_MAC_strength") in {"strong", "moderate"}]
    rel_eb = [finite(row.get("rel_diff_EB")) for row in eligible if math.isfinite(finite(row.get("rel_diff_EB")))]
    rel_timo = [
        finite(row.get("rel_diff_Timoshenko"))
        for row in eligible
        if math.isfinite(finite(row.get("rel_diff_Timoshenko")))
    ]
    best_rel = [
        min(finite(row.get("rel_diff_EB")), finite(row.get("rel_diff_Timoshenko")))
        for row in eligible
        if math.isfinite(finite(row.get("rel_diff_EB"))) and math.isfinite(finite(row.get("rel_diff_Timoshenko")))
    ]
    return {
        "model_id": model_id,
        "beta_deg": beta,
        "eligible_mac_rows_first6": len(eligible),
        "mean_rel_error_EB": float(np.mean(rel_eb)) if rel_eb else float("nan"),
        "mean_rel_error_Timoshenko": float(np.mean(rel_timo)) if rel_timo else float("nan"),
        "mean_best_rel_error": float(np.mean(best_rel)) if best_rel else float("nan"),
        "max_rel_error_EB": max(rel_eb) if rel_eb else float("nan"),
        "max_rel_error_Timoshenko": max(rel_timo) if rel_timo else float("nan"),
        "max_best_rel_error": max(best_rel) if best_rel else float("nan"),
        "weak_or_missing_rows_first6": N_SUMMARY_MODES - len(eligible),
    }


def classify_models(rows: Sequence[dict[str, object]], summary_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    by_model_beta = {(row["model_id"], float(row["beta_deg"])): row for row in summary_rows}
    baseline = by_model_beta.get(("rigid_end_faces", 15.0), {})
    baseline_error = finite(baseline.get("mean_best_rel_error"))
    classifications = []
    for model in JOINT_MODELS:
        beta0 = by_model_beta.get((model.model_id, 0.0), {})
        beta15 = by_model_beta.get((model.model_id, 15.0), {})
        solver_issue = any(
            row.get("model_id") == model.model_id
            and row.get("row_kind") == "case_summary"
            and str(row.get("calculix_success")) not in {"True", "true", "1"}
            for row in rows
        )
        beta0_pass = (
            finite(beta0.get("max_best_rel_error")) <= SANITY_TOL
            and int(beta0.get("eligible_mac_rows_first6", 0)) >= N_SUMMARY_MODES
        )
        beta15_error = finite(beta15.get("mean_best_rel_error"))
        beta15_max = finite(beta15.get("max_best_rel_error"))
        reduced = math.isfinite(beta15_error) and math.isfinite(baseline_error) and beta15_error < baseline_error
        if solver_issue:
            status = "rejected: solver/constraint issue"
        elif not beta0_pass:
            status = "rejected: fails beta=0 sanity"
        elif beta15_max <= TARGET_TOL and int(beta15.get("eligible_mac_rows_first6", 0)) >= N_SUMMARY_MODES:
            status = "candidate for article validation"
        elif reduced:
            status = "promising: beta=0 passes and beta=15 errors reduced"
        else:
            status = "rejected: too stiff/large errors at beta=15"
        classifications.append(
            {
                "model_id": model.model_id,
                "model_label": model.label,
                "classification": status,
                "beta0_pass": beta0_pass,
                "beta15_reduced_vs_rigid_end_faces": reduced,
                "beta15_mean_best_rel_error": beta15_error,
                "beta15_max_best_rel_error": beta15_max,
            }
        )
    return classifications


def append_summary_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    summary_rows = []
    for model in JOINT_MODELS:
        for beta in BETA_VALUES:
            summary_rows.append({**summarize_case(rows, model.model_id, beta), "row_kind": "model_beta_summary"})
    classifications = classify_models(rows, summary_rows)
    for row in summary_rows:
        rows.append(row)
    for row in classifications:
        rows.append({**row, "row_kind": "model_classification"})
    rows.append(
        {
            "row_kind": "unsupported_model_note",
            "model_id": "equation_average",
            "classification": "not run",
            "notes": "Translation-average *EQUATION coupling was not promoted because rotations are not represented; this remains a separate constraint-syntax probe.",
        }
    )
    return classifications


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


def fmt(value: object, digits: int = 6) -> str:
    number = finite(value)
    if not math.isfinite(number):
        return "" if value in {"", None} else str(value)
    return f"{number:.{digits}g}"


def table(rows: Sequence[dict[str, object]], columns: Sequence[tuple[str, str]]) -> list[str]:
    lines = ["| " + " | ".join(label for label, _key in columns) + " |"]
    lines.append("| " + " | ".join("---" for _ in columns) + " |")
    for row in rows:
        values = []
        for _label, key in columns:
            value = row.get(key, "")
            values.append(fmt(value, 5) if isinstance(value, (int, float, np.floating)) else str(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def best_candidate(classifications: Sequence[dict[str, object]]) -> str:
    candidates = [
        row
        for row in classifications
        if str(row.get("classification", "")).startswith("candidate")
        or str(row.get("classification", "")).startswith("promising")
    ]
    if not candidates:
        return "none"
    return str(min(candidates, key=lambda row: finite(row.get("beta15_mean_best_rel_error")))["model_id"])


def write_report(rows: Sequence[dict[str, object]], classifications: Sequence[dict[str, object]], messages: Sequence[str]) -> None:
    summary_rows = [row for row in rows if row.get("row_kind") == "model_beta_summary"]
    constraint_rows = [row for row in rows if row.get("row_kind") == "constraint_audit"]
    mesh_rows = [row for row in rows if row.get("row_kind") == "mesh_audit"]
    warning_rows = [row for row in rows if row.get("row_kind") == "calculix_warning_audit"]
    lines = [
        "# 3D FEM Joint Coupling Models Audit",
        "",
        "## Purpose",
        "",
        "Diagnostic-only study of alternative 3D point-joint coupling idealizations for two coupled rods. "
        "The goal is to find a 3D FEM joint model that behaves like the analytic 1D point-joint model, not to create presentation figures or article-ready validation.",
        "",
        "## Parameters",
        "",
        f"- epsilon: {EPSILON:g}",
        f"- mu: {MU:g}",
        f"- eta: {ETA:g}",
        f"- beta values: {', '.join(f'{value:g}' for value in BETA_VALUES)}",
        f"- analytic modes: {N_ANALYTIC_MODES}; summary uses first {N_SUMMARY_MODES}",
        f"- solid modes requested: {N_SOLID_MODES}",
        f"- mesh factor: {MESH_FACTOR:g}",
        "",
        "## Joint Models",
        "",
        "- `rigid_end_faces`: all inner end-face nodes tied to one `JOINT_REF` by `*RIGID BODY`.",
        "- `center_patch_0p25`: only nodes inside radius `0.25*r` on each inner end face are tied.",
        "- `center_patch_0p5`: only nodes inside radius `0.5*r` on each inner end face are tied.",
        "- `equation_average`: not run; translation-average `*EQUATION` coupling would omit rotations and remains a separate syntax probe.",
        "",
        "## Required Classification",
        "",
        *table(
            classifications,
            [
                ("model", "model_id"),
                ("classification", "classification"),
                ("beta0 pass", "beta0_pass"),
                ("beta15 reduced", "beta15_reduced_vs_rigid_end_faces"),
                ("beta15 mean best rel", "beta15_mean_best_rel_error"),
                ("beta15 max best rel", "beta15_max_best_rel_error"),
            ],
        ),
        "",
        f"Best current candidate: `{best_candidate(classifications)}`.",
        "",
        "No model is article-ready unless it is classified as `candidate for article validation`; promising rows still require visual mode-shape review, mesh refinement, and unique-assignment review.",
        "",
        "## Model/Beta Metrics",
        "",
        *table(
            summary_rows,
            [
                ("model", "model_id"),
                ("beta", "beta_deg"),
                ("eligible rows", "eligible_mac_rows_first6"),
                ("mean rel EB", "mean_rel_error_EB"),
                ("mean rel Timo", "mean_rel_error_Timoshenko"),
                ("mean best rel", "mean_best_rel_error"),
                ("max best rel", "max_best_rel_error"),
                ("weak/missing", "weak_or_missing_rows_first6"),
            ],
        ),
        "",
        "## Constraint Audit",
        "",
        *table(
            constraint_rows,
            [
                ("model", "model_id"),
                ("beta", "beta_deg"),
                ("r1 fixed", "rod1_fixed_node_count"),
                ("r2 fixed", "rod2_fixed_node_count"),
                ("r1 coupled", "rod1_inner_coupled_node_count"),
                ("r2 coupled", "rod2_inner_coupled_node_count"),
                ("joint", "joint_coupled_node_count"),
                ("patch fallback", "patch_fallback_used"),
                ("SPC conflict", "dependent_node_spc_conflict"),
                ("ref restrictions", "reference_restrictions"),
            ],
        ),
        "",
        "## Mesh Audit",
        "",
        *table(
            mesh_rows,
            [
                ("model", "model_id"),
                ("beta", "beta_deg"),
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
        "## CalculiX Warnings",
        "",
    ]
    if warning_rows:
        for row in warning_rows:
            lines.append(f"- {row['model_id']} beta={fmt(row['beta_deg'])}: count={row.get('warning_count')}; {row.get('warnings', '')}")
    else:
        lines.append("- none recorded")
    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- This is diagnostic-only and does not change existing FEM physical models or baseline outputs.",
            "- Center-patch coupling can reduce end-disk artificial rigidity but may under-transfer rotations/moments; it is not a validated 1D joint by construction.",
            "- The beta=0 cases use an explicit straight two-half-rod geometry to avoid the overlapping-cylinder beta=0 geometry of the older coupled helper.",
            "- MAC matching is centerline-displacement based and still needs visual shape review.",
            "- Only mesh factor 1.0 was run in this pass.",
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
    rows: list[dict[str, object]] = []
    messages: list[str] = []
    for model in JOINT_MODELS:
        for beta in BETA_VALUES:
            case = AuditCase(beta_deg=float(beta), model=model)
            case_rows, case_messages = process_case(pj, single, case)
            rows.extend(case_rows)
            messages.extend(case_messages)
    classifications = append_summary_rows(rows)
    write_csv(rows)
    write_report(rows, classifications, messages)
    print("3D FEM joint coupling models audit")
    for row in classifications:
        print(f"{row['model_id']}: {row['classification']}")
    print(f"Best current candidate: {best_candidate(classifications)}")
    print(f"Wrote {OUTPUT_CSV.relative_to(REPO_ROOT)}")
    print(f"Wrote {OUTPUT_REPORT.relative_to(REPO_ROOT)}")
    print(f"Wrote generated FEM files under {AUDIT_ROOT.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
