from __future__ import annotations

import csv
import importlib.util
import math
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Sequence


# Diagnostic-only 3D solid FEM workflow for two equal coupled circular rods.
#
# This script is intentionally separate from the one-rod solid FEM workflow
# because the geometry, boundary-condition sets, mode classification basis, and
# comparison contract are different. It does not make Gmsh or CalculiX required
# project dependencies.

BETA_DEG = 15.0
EPSILON_VALUES = [0.025, 0.05]
N_ANALYTIC_MODES = 6
N_SOLID_MODES = 30
L_SEGMENT = 1.0
E = 1.0
RHO = 1.0
NU = 0.3
THICKNESS_RATIO_LIMIT = 0.1

GMSH_ELEMENT_ORDER = 2
GMSH_EXE_ENV = "GMSH_EXE"
GMSH_EXE_DEFAULT = r"D:\PHD\gmsh\gmsh-4.15.2-Windows64\gmsh.exe"
CCX_EXE_ENV = "CCX_EXE"
CCX_EXE_DEFAULT = r"D:\PHD\calculix\calculix_2.22_4win\ccx_static.exe"
CCX_SEARCH_ROOTS = [r"D:\PHD\calculix\calculix_2.22_4win"]
CCX_TIMEOUT_SECONDS = 360
MODE_SHAPE_SLICE_COUNT = 31
CUTOFF_WARNING_RATIO = 0.8

MU = 0.0
ETA = 0.0
TAU1 = 1.0
TAU2 = 1.0

REPO_ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"
OUTPUT_DIR = REPO_ROOT / "results"
SOLID_OUTPUT_DIR = OUTPUT_DIR / "solid_fem_coupled_equal_rods_beta15"
OUTPUT_REPORT = OUTPUT_DIR / "coupled_equal_rods_beta15_3d_solid_fem_report.md"
OUTPUT_SORTED_CSV = OUTPUT_DIR / "coupled_equal_rods_beta15_3d_solid_fem_sorted_comparison.csv"
OUTPUT_MODE_METRICS_CSV = OUTPUT_DIR / "coupled_equal_rods_beta15_3d_solid_fem_mode_metrics.csv"
OUTPUT_IN_PLANE_CSV = OUTPUT_DIR / "coupled_equal_rods_beta15_3d_solid_fem_in_plane_comparison.csv"

FRD_NUMBER_PATTERN = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")


@dataclass(frozen=True)
class ToolAudit:
    gmsh_exe: str | None
    gmsh_source: str | None
    gmsh_version: str | None
    ccx_exe: str | None
    ccx_source: str | None
    ccx_version: str | None


@dataclass(frozen=True)
class RodFrame:
    outer: tuple[float, float, float]
    joint: tuple[float, float, float]
    tangent: tuple[float, float, float]
    in_plane_normal: tuple[float, float, float]


@dataclass(frozen=True)
class Geometry:
    joint: tuple[float, float, float]
    rod1: RodFrame
    rod2: RodFrame
    beta_rad: float


@dataclass(frozen=True)
class CasePaths:
    case_dir: Path
    geo: Path
    msh: Path
    gmsh_inp: Path
    ccx_mesh_inp: Path
    ccx_modal_inp: Path
    ccx_dat: Path
    ccx_frd: Path
    ccx_stdout: Path
    ccx_stderr: Path


@dataclass(frozen=True)
class MeshData:
    nodes: dict[int, tuple[float, float, float]]
    solid_elements: dict[int, list[int]]
    bbox_min: tuple[float, float, float] | None
    bbox_max: tuple[float, float, float] | None


@dataclass(frozen=True)
class MeshSummary:
    epsilon: float
    node_count: int
    solid_element_count: int
    bbox_min: tuple[float, float, float] | None
    bbox_max: tuple[float, float, float] | None
    connected_component_count: int
    largest_component_element_count: int
    geometry_connected: bool
    gmsh_messages: tuple[str, ...]


@dataclass(frozen=True)
class CcxInputSummary:
    epsilon: float
    modal_input_path: Path
    mesh_include_path: Path
    rod1_fixed_node_count: int
    rod2_fixed_node_count: int
    node_sets_present: bool


@dataclass(frozen=True)
class CalculixResult:
    epsilon: float
    success: bool
    message: str
    stdout_path: Path
    stderr_path: Path
    dat_path: Path
    frd_path: Path
    parsed_omegas: tuple[float, ...]
    parse_source: str


@dataclass(frozen=True)
class ModeShapeParseResult:
    epsilon: float
    success: bool
    message: str
    frd_path: Path
    mode_shapes: dict[int, dict[int, tuple[float, float, float]]]
    source: str


@dataclass(frozen=True)
class ModeMetric:
    epsilon: float
    solid_mode: int
    omega: float
    lambda_fem: float
    parsed_node_count: int
    in_plane_energy_fraction: float
    out_of_plane_energy_fraction: float
    axial_energy_fraction: float
    transverse_in_plane_energy_fraction: float
    mean_axial_energy_fraction: float
    mean_in_plane_bending_fraction: float
    mean_out_of_plane_bending_fraction: float
    torsion_indicator: float
    joint_motion_fraction: float
    warping_indicator: float
    classification: str
    classification_note: str


def safe_float_token(value: float) -> str:
    return f"{float(value):.6g}".replace("-", "m").replace(".", "p")


def existing_executable(path_text: str | None) -> str | None:
    if not path_text:
        return None
    path = Path(path_text.strip().strip('"'))
    if path.is_file():
        return str(path)
    return None


def resolve_gmsh_executable() -> tuple[str | None, str | None]:
    env_value = existing_executable(os.environ.get(GMSH_EXE_ENV))
    if env_value is not None:
        return env_value, f"environment variable {GMSH_EXE_ENV}"
    default_value = existing_executable(GMSH_EXE_DEFAULT)
    if default_value is not None:
        return default_value, "script constant GMSH_EXE_DEFAULT"
    path_value = shutil.which("gmsh")
    if path_value is not None:
        return path_value, "PATH lookup"
    return None, None


def resolve_ccx_executable() -> tuple[str | None, str | None]:
    env_value = existing_executable(os.environ.get(CCX_EXE_ENV))
    if env_value is not None:
        return env_value, f"environment variable {CCX_EXE_ENV}"
    default_value = existing_executable(CCX_EXE_DEFAULT)
    if default_value is not None:
        return default_value, "script constant CCX_EXE_DEFAULT"
    for root_text in CCX_SEARCH_ROOTS:
        root = Path(root_text)
        for candidate_name in ("ccx.exe", "ccx_static.exe", "ccx_dynamic.exe"):
            candidate = root / candidate_name
            if candidate.is_file():
                return str(candidate), f"search root {root}"
    path_value = shutil.which("ccx")
    if path_value is not None:
        return path_value, "PATH lookup"
    return None, None


def executable_version(executable: str | None, *args: str) -> str | None:
    if executable is None:
        return None
    try:
        completed = subprocess.run(
            [executable, *args],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=20,
        )
    except Exception:
        return None
    text = (completed.stdout or completed.stderr).strip()
    return text.splitlines()[0].strip() if text else None


def audit_tools() -> ToolAudit:
    gmsh_exe, gmsh_source = resolve_gmsh_executable()
    ccx_exe, ccx_source = resolve_ccx_executable()
    return ToolAudit(
        gmsh_exe=gmsh_exe,
        gmsh_source=gmsh_source,
        gmsh_version=executable_version(gmsh_exe, "--version"),
        ccx_exe=ccx_exe,
        ccx_source=ccx_source,
        ccx_version=executable_version(ccx_exe, "-v"),
    )


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def mul(value: float, vector: tuple[float, float, float]) -> tuple[float, float, float]:
    return (value * vector[0], value * vector[1], value * vector[2])


def norm(vector: tuple[float, float, float]) -> float:
    return math.sqrt(dot(vector, vector))


def unit(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    length = norm(vector)
    if length == 0.0:
        raise ValueError("zero vector cannot be normalized")
    return (vector[0] / length, vector[1] / length, vector[2] / length)


def geometry() -> Geometry:
    beta = math.radians(BETA_DEG)
    cos_b = math.cos(beta)
    sin_b = math.sin(beta)
    joint = (0.0, 0.0, 0.0)
    rod1_outer = (-L_SEGMENT, 0.0, 0.0)
    rod2_outer = (-L_SEGMENT * cos_b, -L_SEGMENT * sin_b, 0.0)
    rod1_tangent = unit(sub(joint, rod1_outer))
    rod2_tangent = unit(sub(joint, rod2_outer))
    rod1_normal = (-rod1_tangent[1], rod1_tangent[0], 0.0)
    rod2_normal = (-rod2_tangent[1], rod2_tangent[0], 0.0)
    return Geometry(
        joint=joint,
        rod1=RodFrame(outer=rod1_outer, joint=joint, tangent=rod1_tangent, in_plane_normal=rod1_normal),
        rod2=RodFrame(outer=rod2_outer, joint=joint, tangent=rod2_tangent, in_plane_normal=rod2_normal),
        beta_rad=beta,
    )


def radius_from_epsilon(epsilon: float) -> float:
    return 2.0 * float(epsilon) * L_SEGMENT


def area_from_radius(radius: float) -> float:
    return math.pi * radius**2


def inertia_from_radius(radius: float) -> float:
    return math.pi * radius**4 / 4.0


def shear_modulus() -> float:
    return E / (2.0 * (1.0 + NU))


def circular_shear_coefficient() -> float:
    return (6.0 + 12.0 * NU + 6.0 * NU**2) / (7.0 + 12.0 * NU + 4.0 * NU**2)


def omega_cutoff(epsilon: float) -> float:
    radius = radius_from_epsilon(float(epsilon))
    area = area_from_radius(radius)
    inertia = inertia_from_radius(radius)
    return math.sqrt(circular_shear_coefficient() * shear_modulus() * area / (RHO * inertia))


def lambda_fem_from_omega(omega: float, epsilon: float) -> float:
    if omega <= 0.0 or epsilon <= 0.0:
        return float("nan")
    return math.sqrt(float(omega) / float(epsilon))


def project_omega(lambda_value: float, epsilon: float) -> float:
    return float(epsilon) * float(lambda_value) ** 2


def diameter_to_segment_length(epsilon: float) -> float:
    return 4.0 * float(epsilon)


def thin_rod_valid(epsilon: float) -> bool:
    return diameter_to_segment_length(float(epsilon)) <= THICKNESS_RATIO_LIMIT + 1.0e-15


def mesh_size_for_epsilon(epsilon: float) -> float:
    radius = radius_from_epsilon(float(epsilon))
    return min(L_SEGMENT / 40.0, radius / 2.0)


def case_paths(epsilon: float) -> CasePaths:
    token = safe_float_token(float(epsilon))
    case_dir = SOLID_OUTPUT_DIR / f"eps_{token}"
    stem = f"coupled_equal_rods_beta15_eps{token}"
    return CasePaths(
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


def write_gmsh_geo(epsilon: float, paths: CasePaths, *, boolean_mode: str) -> str:
    radius = radius_from_epsilon(float(epsilon))
    mesh_size = mesh_size_for_epsilon(float(epsilon))
    geom = geometry()
    cos_b = math.cos(geom.beta_rad)
    sin_b = math.sin(geom.beta_rad)
    paths.case_dir.mkdir(parents=True, exist_ok=True)
    if boolean_mode == "union":
        sphere_line = ""
        boolean_line = "fused[] = BooleanUnion{ Volume{1}; Delete; }{ Volume{2}; Delete; };"
        geometry_label = "fused cylinders (OpenCASCADE BooleanUnion)"
    else:
        sphere_line = "Sphere(3) = {0, 0, 0, R};"
        boolean_line = "fused[] = BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Volume{3}; Delete; };"
        geometry_label = "fused cylinders + spherical joint (OpenCASCADE BooleanFragments fallback)"
    bbox = L_SEGMENT + 2.0 * radius
    paths.geo.write_text(
        f"""// Diagnostic-only two equal coupled rods solid mesh.
// beta = {BETA_DEG:.17g} deg, epsilon = {float(epsilon):.17g}
SetFactory("OpenCASCADE");

L = {L_SEGMENT:.17g};
R = {radius:.17g};
h = {mesh_size:.17g};
CosB = {cos_b:.17g};
SinB = {sin_b:.17g};

Cylinder(1) = {{-L, 0, 0, L, 0, 0, R, 2*Pi}};
Cylinder(2) = {{-L*CosB, -L*SinB, 0, L*CosB, L*SinB, 0, R, 2*Pi}};
{sphere_line}

{boolean_line}

vols[] = Volume In BoundingBox {{-L-{bbox:.17g}, -L-{bbox:.17g}, -R-{bbox:.17g}, R+{bbox:.17g}, R+{bbox:.17g}, R+{bbox:.17g}}};
Physical Volume("SOLID", 1) = vols[];

Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Mesh.ElementOrder = {GMSH_ELEMENT_ORDER};
Mesh.MshFileVersion = 4.1;
""",
        encoding="utf-8",
    )
    return geometry_label


def interesting_external_lines(text: str, *, prefix: str) -> list[str]:
    lines: list[str] = []
    for line in text.splitlines():
        lower = line.lower()
        if "warning" in lower or "error" in lower:
            lines.append(f"{prefix}: {line.strip()}")
    return lines


def generate_mesh_with_gmsh_cli(paths: CasePaths, gmsh_exe: str) -> tuple[bool, str, tuple[str, ...]]:
    commands = [
        ("msh4", [gmsh_exe, str(paths.geo), "-3", "-format", "msh4", "-o", str(paths.msh)]),
        ("inp", [gmsh_exe, str(paths.geo), "-3", "-format", "inp", "-o", str(paths.gmsh_inp)]),
    ]
    messages: list[str] = []
    for label, command in commands:
        try:
            completed = subprocess.run(
                command,
                cwd=str(paths.case_dir),
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=240,
            )
        except Exception as exc:
            return False, f"Gmsh {label} generation failed before completion: {exc}", tuple(messages)
        messages.extend(interesting_external_lines(completed.stdout, prefix=f"gmsh {label} stdout"))
        messages.extend(interesting_external_lines(completed.stderr, prefix=f"gmsh {label} stderr"))
        if completed.returncode != 0:
            tail = (completed.stdout + "\n" + completed.stderr)[-800:]
            return False, f"Gmsh {label} generation returned {completed.returncode}: {tail!r}", tuple(messages)
    return True, "mesh generated successfully with Gmsh executable", tuple(messages)


def generate_mesh(paths: CasePaths, epsilon: float, gmsh_exe: str) -> tuple[bool, str, str, tuple[str, ...]]:
    geometry_label = write_gmsh_geo(float(epsilon), paths, boolean_mode="union")
    ok, message, messages = generate_mesh_with_gmsh_cli(paths, gmsh_exe)
    if ok:
        return ok, message, geometry_label, messages
    fallback_label = write_gmsh_geo(float(epsilon), paths, boolean_mode="fragments")
    fallback_ok, fallback_message, fallback_messages = generate_mesh_with_gmsh_cli(paths, gmsh_exe)
    all_messages = tuple(list(messages) + [f"union attempt failed: {message}"] + list(fallback_messages))
    return fallback_ok, fallback_message, fallback_label, all_messages


def update_bbox(
    current_min: tuple[float, float, float] | None,
    current_max: tuple[float, float, float] | None,
    point: tuple[float, float, float],
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    if current_min is None or current_max is None:
        return point, point
    return (
        tuple(min(a, b) for a, b in zip(current_min, point)),
        tuple(max(a, b) for a, b in zip(current_max, point)),
    )


def parse_keyword_args(line: str) -> dict[str, str]:
    args: dict[str, str] = {}
    for part in line.split(",")[1:]:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        args[key.strip().upper()] = value.strip()
    return args


def read_gmsh_inp_mesh_data(path: Path) -> MeshData:
    nodes: dict[int, tuple[float, float, float]] = {}
    solid_elements: dict[int, list[int]] = {}
    bbox_min: tuple[float, float, float] | None = None
    bbox_max: tuple[float, float, float] | None = None
    mode = ""
    current_element_type = ""
    for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("*"):
            upper = line.upper()
            args = parse_keyword_args(line)
            if upper.startswith("*NODE"):
                mode = "NODE"
            elif upper.startswith("*ELEMENT"):
                mode = "ELEMENT"
                current_element_type = args.get("TYPE", "").upper()
            else:
                mode = ""
            continue
        parts = [part.strip() for part in line.split(",") if part.strip()]
        if mode == "NODE" and len(parts) >= 4:
            try:
                node_id = int(parts[0])
                point = (float(parts[1]), float(parts[2]), float(parts[3]))
            except ValueError:
                continue
            nodes[node_id] = point
            bbox_min, bbox_max = update_bbox(bbox_min, bbox_max, point)
        elif mode == "ELEMENT" and current_element_type.startswith("C3D") and len(parts) >= 2:
            try:
                solid_elements[int(parts[0])] = [int(part) for part in parts[1:]]
            except ValueError:
                continue
    return MeshData(nodes=nodes, solid_elements=solid_elements, bbox_min=bbox_min, bbox_max=bbox_max)


def connected_component_counts(mesh: MeshData) -> tuple[int, int]:
    parent: dict[int, int] = {}

    def find(node: int) -> int:
        parent.setdefault(node, node)
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(a: int, b: int) -> None:
        root_a = find(a)
        root_b = find(b)
        if root_a != root_b:
            parent[root_b] = root_a

    for connectivity in mesh.solid_elements.values():
        if not connectivity:
            continue
        first = connectivity[0]
        for node in connectivity[1:]:
            union(first, node)

    component_elements: dict[int, int] = defaultdict(int)
    for connectivity in mesh.solid_elements.values():
        if connectivity:
            component_elements[find(connectivity[0])] += 1
    if not component_elements:
        return 0, 0
    return len(component_elements), max(component_elements.values())


def mesh_summary(epsilon: float, mesh: MeshData, gmsh_messages: Sequence[str]) -> MeshSummary:
    components, largest = connected_component_counts(mesh)
    return MeshSummary(
        epsilon=float(epsilon),
        node_count=len(mesh.nodes),
        solid_element_count=len(mesh.solid_elements),
        bbox_min=mesh.bbox_min,
        bbox_max=mesh.bbox_max,
        connected_component_count=components,
        largest_component_element_count=largest,
        geometry_connected=components == 1,
        gmsh_messages=tuple(gmsh_messages),
    )


def projection_on_rod(point: tuple[float, float, float], rod: RodFrame) -> tuple[float, float]:
    rel = sub(point, rod.outer)
    s = dot(rel, rod.tangent)
    axis_point = add(rod.outer, mul(s, rod.tangent))
    radial = norm(sub(point, axis_point))
    return s, radial


def coordinate_fixed_node_sets(mesh: MeshData, epsilon: float) -> tuple[list[int], list[int]]:
    geom = geometry()
    radius = radius_from_epsilon(float(epsilon))
    plane_tol = max(1.0e-7 * L_SEGMENT, 1.0e-6)
    radial_tol = 1.02 * radius + plane_tol
    rod1_nodes: list[int] = []
    rod2_nodes: list[int] = []
    for node_id, point in mesh.nodes.items():
        s1, r1 = projection_on_rod(point, geom.rod1)
        s2, r2 = projection_on_rod(point, geom.rod2)
        if abs(s1) <= plane_tol and r1 <= radial_tol:
            rod1_nodes.append(node_id)
        if abs(s2) <= plane_tol and r2 <= radial_tol:
            rod2_nodes.append(node_id)
    return sorted(rod1_nodes), sorted(rod2_nodes)


def ccx_float(value: float) -> str:
    return f"{float(value):.12E}"


def format_id_lines(ids: Sequence[int], *, per_line: int = 16) -> list[str]:
    lines: list[str] = []
    for start in range(0, len(ids), per_line):
        lines.append(", ".join(str(item) for item in ids[start : start + per_line]))
    return lines


def write_calculix_mesh_include(paths: CasePaths, mesh: MeshData) -> None:
    lines = [
        "** Sanitized diagnostic-only mesh include for CalculiX.",
        "** Generated from Gmsh .inp by keeping nodes and 3D solid elements only.",
        "*NODE",
    ]
    for node_id in sorted(mesh.nodes):
        x, y, z = mesh.nodes[node_id]
        lines.append(f"{node_id}, {ccx_float(x)}, {ccx_float(y)}, {ccx_float(z)}")
    lines.append("*ELEMENT, TYPE=C3D10, ELSET=SOLID")
    for element_id in sorted(mesh.solid_elements):
        connectivity = ", ".join(str(node_id) for node_id in mesh.solid_elements[element_id])
        lines.append(f"{element_id}, {connectivity}")
    paths.ccx_mesh_inp.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_calculix_inputs(paths: CasePaths, epsilon: float, mesh: MeshData) -> CcxInputSummary:
    rod1_nodes, rod2_nodes = coordinate_fixed_node_sets(mesh, float(epsilon))
    write_calculix_mesh_include(paths, mesh)
    paths.ccx_modal_inp.write_text(
        "\n".join(
            [
                "** Diagnostic-only CalculiX modal input for two fused equal circular rods.",
                "** Fixed node sets are generated from outer-end coordinate projections.",
                f"*INCLUDE, INPUT={paths.ccx_mesh_inp.name}",
                "*NSET, NSET=ROD1_OUTER_FIXED",
                *format_id_lines(rod1_nodes),
                "*NSET, NSET=ROD2_OUTER_FIXED",
                *format_id_lines(rod2_nodes),
                "*MATERIAL, NAME=MAT",
                "*ELASTIC",
                f"{E:.17g}, {NU:.17g}",
                "*DENSITY",
                f"{RHO:.17g}",
                "*SOLID SECTION, ELSET=SOLID, MATERIAL=MAT",
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
        )
        + "\n",
        encoding="utf-8",
    )
    return CcxInputSummary(
        epsilon=float(epsilon),
        modal_input_path=paths.ccx_modal_inp,
        mesh_include_path=paths.ccx_mesh_inp,
        rod1_fixed_node_count=len(rod1_nodes),
        rod2_fixed_node_count=len(rod2_nodes),
        node_sets_present=bool(rod1_nodes and rod2_nodes),
    )


def parse_calculix_eigen_omegas(dat_path: Path, stdout_path: Path) -> tuple[list[float], str]:
    if not dat_path.exists() and not stdout_path.exists():
        return [], "missing .dat/stdout"
    omegas: list[float] = []
    dat_text = dat_path.read_text(encoding="utf-8", errors="ignore") if dat_path.exists() else ""
    for line in dat_text.splitlines():
        parts = line.split()
        if len(parts) < 2 or not parts[0].isdigit():
            continue
        try:
            eigenvalue = float(parts[1].replace("D", "E").replace("d", "E"))
        except ValueError:
            continue
        if eigenvalue > 0.0:
            omegas.append(math.sqrt(eigenvalue))
        if len(omegas) >= N_SOLID_MODES:
            break
    if omegas:
        return omegas, "CalculiX .dat numeric table second column interpreted as Omega^2"

    frequency_pattern = re.compile(r"frequency\s*[:=]?\s*([-+0-9.DEded]+)", re.IGNORECASE)
    stdout_text = stdout_path.read_text(encoding="utf-8", errors="ignore") if stdout_path.exists() else ""
    cyclic_frequencies: list[float] = []
    for text in (dat_text, stdout_text):
        for match in frequency_pattern.finditer(text):
            try:
                value = float(match.group(1).replace("D", "E").replace("d", "E"))
            except ValueError:
                continue
            if value > 0.0:
                cyclic_frequencies.append(value)
    if cyclic_frequencies:
        return [2.0 * math.pi * value for value in cyclic_frequencies[:N_SOLID_MODES]], (
            "CalculiX frequency entries interpreted as cyclic frequency f; Omega=2*pi*f"
        )
    return [], "no recognizable CalculiX eigenfrequency entries"


def run_calculix(paths: CasePaths, ccx: str, epsilon: float) -> CalculixResult:
    job_name = paths.ccx_modal_inp.stem
    for stale_path in (paths.ccx_dat, paths.ccx_frd, paths.ccx_stdout, paths.ccx_stderr):
        try:
            stale_path.unlink(missing_ok=True)
        except OSError:
            pass
    try:
        completed = subprocess.run(
            [ccx, job_name],
            cwd=str(paths.case_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=CCX_TIMEOUT_SECONDS,
        )
    except Exception as exc:
        paths.ccx_stdout.write_text("", encoding="utf-8")
        paths.ccx_stderr.write_text(str(exc), encoding="utf-8")
        return CalculixResult(
            epsilon=float(epsilon),
            success=False,
            message=f"CalculiX run failed before completion: {exc}",
            stdout_path=paths.ccx_stdout,
            stderr_path=paths.ccx_stderr,
            dat_path=paths.ccx_dat,
            frd_path=paths.ccx_frd,
            parsed_omegas=(),
            parse_source="run_exception",
        )
    paths.ccx_stdout.write_text(completed.stdout, encoding="utf-8")
    paths.ccx_stderr.write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        tail = (completed.stdout + "\n" + completed.stderr)[-800:]
        return CalculixResult(
            epsilon=float(epsilon),
            success=False,
            message=f"CalculiX returned exit code {completed.returncode}: {tail!r}",
            stdout_path=paths.ccx_stdout,
            stderr_path=paths.ccx_stderr,
            dat_path=paths.ccx_dat,
            frd_path=paths.ccx_frd,
            parsed_omegas=(),
            parse_source="solver_failed",
        )
    omegas, parse_source = parse_calculix_eigen_omegas(paths.ccx_dat, paths.ccx_stdout)
    return CalculixResult(
        epsilon=float(epsilon),
        success=bool(omegas),
        message=f"CalculiX completed; parsed {len(omegas)} frequencies from {parse_source}.",
        stdout_path=paths.ccx_stdout,
        stderr_path=paths.ccx_stderr,
        dat_path=paths.ccx_dat,
        frd_path=paths.ccx_frd,
        parsed_omegas=tuple(omegas),
        parse_source=parse_source,
    )


def numeric_tokens(text: str) -> list[float]:
    values: list[float] = []
    for token in FRD_NUMBER_PATTERN.findall(text):
        try:
            values.append(float(token.replace("D", "E").replace("d", "E")))
        except ValueError:
            continue
    return values


def parse_calculix_frd_mode_shapes(paths: CasePaths, epsilon: float) -> ModeShapeParseResult:
    if not paths.ccx_frd.exists():
        return ModeShapeParseResult(
            epsilon=float(epsilon),
            success=False,
            message="CalculiX .frd file is missing; no mode-shape vectors parsed.",
            frd_path=paths.ccx_frd,
            mode_shapes={},
            source="missing_frd",
        )
    mode_shapes: dict[int, dict[int, tuple[float, float, float]]] = {}
    current_mode: int | None = None
    collecting = False
    with paths.ccx_frd.open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            stripped = raw_line.strip()
            if not stripped:
                continue
            if stripped.startswith("1PMODE"):
                parts = stripped.split()
                try:
                    current_mode = int(parts[-1])
                except (ValueError, IndexError):
                    current_mode = None
                collecting = False
                continue
            if stripped.startswith("-4") and "DISP" in stripped.upper():
                if current_mode is not None:
                    mode_shapes.setdefault(current_mode, {})
                    collecting = True
                continue
            if collecting and stripped.startswith("-1"):
                values = numeric_tokens(stripped[2:])
                if len(values) >= 4 and current_mode is not None:
                    node_id = int(values[0])
                    mode_shapes[current_mode][node_id] = (float(values[1]), float(values[2]), float(values[3]))
                continue
            if collecting and (
                stripped.startswith("-3")
                or stripped.startswith("-4")
                or stripped.startswith("100C")
                or stripped.startswith("1P")
            ):
                collecting = False
    parsed_modes = sorted(mode for mode, shape in mode_shapes.items() if shape)
    if not parsed_modes:
        return ModeShapeParseResult(
            epsilon=float(epsilon),
            success=False,
            message="No displacement eigenvector blocks were recognized in the CalculiX .frd file.",
            frd_path=paths.ccx_frd,
            mode_shapes={},
            source="frd_no_displacement_blocks",
        )
    return ModeShapeParseResult(
        epsilon=float(epsilon),
        success=True,
        message=f"Parsed displacement eigenvectors for {len(parsed_modes)} modes from CalculiX .frd DISP blocks.",
        frd_path=paths.ccx_frd,
        mode_shapes={mode: mode_shapes[mode] for mode in parsed_modes},
        source="CalculiX .frd DISP blocks from *NODE FILE U",
    )


def assigned_rod(point: tuple[float, float, float], epsilon: float) -> tuple[int, float, float, float] | None:
    geom = geometry()
    radius = radius_from_epsilon(float(epsilon))
    candidates: list[tuple[float, int, float, float]] = []
    for rod_index, rod in ((1, geom.rod1), (2, geom.rod2)):
        s, radial = projection_on_rod(point, rod)
        distance_outside = max(0.0, -s, s - L_SEGMENT) + radial
        if -1.5 * radius <= s <= L_SEGMENT + 1.5 * radius:
            candidates.append((distance_outside, rod_index, s, radial))
    if not candidates:
        return None
    _, rod_index, s, radial = min(candidates, key=lambda item: item[0])
    clamped_s = max(0.0, min(L_SEGMENT, s))
    return rod_index, clamped_s, radial, radius


def classify_mode_metric(
    *,
    axial_fraction: float,
    mean_axial_fraction: float,
    mean_in_plane_bending_fraction: float,
    mean_out_of_plane_bending_fraction: float,
    out_of_plane_fraction: float,
    torsion_indicator: float,
) -> tuple[str, str]:
    if axial_fraction >= 0.55 and mean_axial_fraction >= 0.12:
        return "axial_like", "dominant displacement along local rod tangents"
    if torsion_indicator >= 0.55 and mean_in_plane_bending_fraction + mean_out_of_plane_bending_fraction <= 0.25:
        return "torsion_like", "dominant circumferential pattern about rod axes"
    if mean_in_plane_bending_fraction >= 0.18 and mean_in_plane_bending_fraction >= 1.25 * mean_out_of_plane_bending_fraction:
        return "in_plane_bending_like", "dominant centroidal in-plane transverse motion"
    if (
        mean_out_of_plane_bending_fraction >= 0.18
        and out_of_plane_fraction >= 0.35
        and mean_out_of_plane_bending_fraction >= mean_in_plane_bending_fraction
    ):
        return "out_of_plane_bending_like", "dominant centroidal out-of-plane transverse motion"
    return "mixed_or_unclassified", "no single diagnostic indicator dominates"


def compute_mode_metric(
    epsilon: float,
    solid_mode: int,
    omega: float,
    mesh: MeshData,
    mode_shape: dict[int, tuple[float, float, float]],
) -> ModeMetric | None:
    geom = geometry()
    total_energy = 0.0
    in_plane_energy = 0.0
    out_of_plane_energy = 0.0
    axial_energy = 0.0
    transverse_in_plane_energy = 0.0
    joint_energy = 0.0
    parsed_count = 0
    slice_data: dict[tuple[int, int], list[float]] = defaultdict(lambda: [0.0] * 11)
    radius = radius_from_epsilon(float(epsilon))

    for node_id, point in mesh.nodes.items():
        displacement = mode_shape.get(node_id)
        if displacement is None:
            continue
        ux, uy, uz = displacement
        local_total = ux * ux + uy * uy + uz * uz
        if local_total <= 0.0:
            parsed_count += 1
            continue
        total_energy += local_total
        in_plane_energy += ux * ux + uy * uy
        out_of_plane_energy += uz * uz
        if norm(sub(point, geom.joint)) <= 1.5 * radius:
            joint_energy += local_total
        assignment = assigned_rod(point, float(epsilon))
        if assignment is None:
            parsed_count += 1
            continue
        rod_index, s, _radial, _radius = assignment
        rod = geom.rod1 if rod_index == 1 else geom.rod2
        rel = sub(point, add(rod.outer, mul(s, rod.tangent)))
        y_local = dot(rel, rod.in_plane_normal)
        z_local = point[2]
        u_t = dot(displacement, rod.tangent)
        u_n = dot(displacement, rod.in_plane_normal)
        u_z = uz
        axial_energy += u_t * u_t
        transverse_in_plane_energy += u_n * u_n
        slice_index = int(math.floor(max(0.0, min(0.999999999999, s / L_SEGMENT)) * MODE_SHAPE_SLICE_COUNT))
        bucket = slice_data[(rod_index, slice_index)]
        bucket[0] += 1.0
        bucket[1] += u_t
        bucket[2] += u_n
        bucket[3] += u_z
        bucket[4] += local_total
        bucket[5] += -z_local * u_n + y_local * u_z
        bucket[6] += y_local * y_local + z_local * z_local
        bucket[7] += u_t * u_t
        bucket[8] += u_n * u_n
        bucket[9] += u_z * u_z
        bucket[10] += ux + uy + uz
        parsed_count += 1

    if total_energy <= 0.0 or parsed_count == 0:
        return None

    mean_total_energy = 0.0
    mean_axial_energy = 0.0
    mean_in_plane_bending_energy = 0.0
    mean_out_of_plane_bending_energy = 0.0
    torsion_explained_energy = 0.0
    transverse_energy_for_torsion = 0.0
    for bucket in slice_data.values():
        count, sum_ut, sum_un, sum_uz, _sum_total, torsion_num, torsion_den, _sum_ut2, sum_un2, sum_uz2, _sum_all = bucket
        if count <= 0.0:
            continue
        mean_axial_energy += (sum_ut * sum_ut) / count
        mean_in_plane_bending_energy += (sum_un * sum_un) / count
        mean_out_of_plane_bending_energy += (sum_uz * sum_uz) / count
        mean_total_energy += (sum_ut * sum_ut + sum_un * sum_un + sum_uz * sum_uz) / count
        transverse_energy_for_torsion += sum_un2 + sum_uz2
        if torsion_den > 0.0:
            torsion_explained_energy += (torsion_num * torsion_num) / torsion_den

    axial_fraction = axial_energy / total_energy
    mean_axial_fraction = mean_axial_energy / total_energy
    mean_in_plane_bending_fraction = mean_in_plane_bending_energy / total_energy
    mean_out_of_plane_bending_fraction = mean_out_of_plane_bending_energy / total_energy
    torsion_indicator = (
        torsion_explained_energy / transverse_energy_for_torsion if transverse_energy_for_torsion > 0.0 else 0.0
    )
    classification, note = classify_mode_metric(
        axial_fraction=axial_fraction,
        mean_axial_fraction=mean_axial_fraction,
        mean_in_plane_bending_fraction=mean_in_plane_bending_fraction,
        mean_out_of_plane_bending_fraction=mean_out_of_plane_bending_fraction,
        out_of_plane_fraction=out_of_plane_energy / total_energy,
        torsion_indicator=torsion_indicator,
    )
    return ModeMetric(
        epsilon=float(epsilon),
        solid_mode=int(solid_mode),
        omega=float(omega),
        lambda_fem=lambda_fem_from_omega(float(omega), float(epsilon)),
        parsed_node_count=parsed_count,
        in_plane_energy_fraction=in_plane_energy / total_energy,
        out_of_plane_energy_fraction=out_of_plane_energy / total_energy,
        axial_energy_fraction=axial_fraction,
        transverse_in_plane_energy_fraction=transverse_in_plane_energy / total_energy,
        mean_axial_energy_fraction=mean_axial_fraction,
        mean_in_plane_bending_fraction=mean_in_plane_bending_fraction,
        mean_out_of_plane_bending_fraction=mean_out_of_plane_bending_fraction,
        torsion_indicator=torsion_indicator,
        joint_motion_fraction=joint_energy / total_energy,
        warping_indicator=max(0.0, min(1.0, (total_energy - mean_total_energy) / total_energy)),
        classification=classification,
        classification_note=note,
    )


def compute_mode_metrics(
    epsilon: float,
    paths: CasePaths,
    mesh: MeshData,
    omegas: Sequence[float],
) -> tuple[list[ModeMetric], ModeShapeParseResult]:
    parse_result = parse_calculix_frd_mode_shapes(paths, float(epsilon))
    if not parse_result.success:
        return [], parse_result
    metrics: list[ModeMetric] = []
    for mode_index, omega in enumerate(omegas, start=1):
        shape = parse_result.mode_shapes.get(mode_index)
        if shape is None:
            continue
        metric = compute_mode_metric(float(epsilon), mode_index, float(omega), mesh, shape)
        if metric is not None:
            metrics.append(metric)
    return metrics, parse_result


def load_analytic_module() -> ModuleType:
    module_path = ANALYSIS_DIR / "compare_coupled_equal_rods_eb_timoshenko.py"
    spec = importlib.util.spec_from_file_location("coupled_equal_rods_reference", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load analytic reference module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.BETA_DEG = BETA_DEG
    module.MU = MU
    module.N_MODES = N_ANALYTIC_MODES
    module.L_SEGMENT = L_SEGMENT
    module.E = E
    module.RHO = RHO
    module.NU = NU
    return module


def compute_analytic_rows() -> tuple[list[dict[str, object]], list[str]]:
    reference = load_analytic_module()
    rows: list[dict[str, object]] = []
    warnings: list[str] = []
    for epsilon in EPSILON_VALUES:
        eb_values = reference.eb_roots(float(epsilon))
        timo_values, root_warnings = reference.timo_roots(float(epsilon), eb_values)
        warnings.extend(f"epsilon={float(epsilon):g}: {warning}" for warning in root_warnings)
        omega_c = omega_cutoff(float(epsilon))
        lambda_c = math.sqrt(omega_c / float(epsilon))
        valid = thin_rod_valid(float(epsilon))
        diameter_ratio = diameter_to_segment_length(float(epsilon))
        for mode_index in range(N_ANALYTIC_MODES):
            lambda_eb = float(eb_values[mode_index])
            lambda_timo = float(timo_values[mode_index])
            omega_eb = project_omega(lambda_eb, float(epsilon))
            omega_timo = project_omega(lambda_timo, float(epsilon)) if math.isfinite(lambda_timo) else float("nan")
            rows.append(
                {
                    "epsilon": float(epsilon),
                    "diameter_to_segment_length": diameter_ratio,
                    "thin_rod_valid": valid,
                    "mode": mode_index + 1,
                    "Lambda_EB": lambda_eb,
                    "Lambda_Timoshenko": lambda_timo,
                    "Omega_EB": omega_eb,
                    "Omega_Timoshenko": omega_timo,
                    "Omega_cutoff": omega_c,
                    "Omega_EB_over_cutoff": omega_eb / omega_c,
                    "Omega_Timoshenko_over_cutoff": omega_timo / omega_c if math.isfinite(omega_timo) else float("nan"),
                    "Lambda_cutoff": lambda_c,
                    "below_cutoff_Timoshenko": bool(math.isfinite(omega_timo) and omega_timo < omega_c),
                }
            )
    return rows, sorted(set(warnings))


def analytic_by_epsilon_mode(rows: Sequence[dict[str, object]]) -> dict[tuple[float, int], dict[str, object]]:
    return {(float(row["epsilon"]), int(row["mode"])): row for row in rows}


def build_sorted_comparison_rows(
    analytic_rows: Sequence[dict[str, object]],
    fem_omegas_by_epsilon: dict[float, list[float]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for analytic in analytic_rows:
        epsilon = float(analytic["epsilon"])
        mode = int(analytic["mode"])
        fem_omegas = fem_omegas_by_epsilon.get(epsilon, [])
        omega_fem = float(fem_omegas[mode - 1]) if mode <= len(fem_omegas) else float("nan")
        lambda_fem = lambda_fem_from_omega(omega_fem, epsilon) if math.isfinite(omega_fem) else float("nan")
        omega_c = float(analytic["Omega_cutoff"])
        rows.append(
            {
                **analytic,
                "solid_sorted_mode": mode if math.isfinite(omega_fem) else "",
                "Omega_FEM": omega_fem if math.isfinite(omega_fem) else "",
                "Lambda_FEM": lambda_fem if math.isfinite(lambda_fem) else "",
                "Omega_FEM_over_cutoff": omega_fem / omega_c if math.isfinite(omega_fem) else "",
                "below_cutoff_FEM": bool(math.isfinite(omega_fem) and omega_fem < omega_c),
                "status": "computed_preliminary_sorted" if math.isfinite(omega_fem) else "not_computed",
                "notes": "sorted frequency only; do not use for mode identity claims",
            }
        )
    return rows


def build_in_plane_comparison_rows(
    analytic_rows: Sequence[dict[str, object]],
    mode_metrics: Sequence[ModeMetric],
) -> list[dict[str, object]]:
    lookup = analytic_by_epsilon_mode(analytic_rows)
    rows: list[dict[str, object]] = []
    metrics_by_epsilon: dict[float, list[ModeMetric]] = defaultdict(list)
    for metric in mode_metrics:
        if metric.classification == "in_plane_bending_like":
            metrics_by_epsilon[metric.epsilon].append(metric)
    for epsilon in EPSILON_VALUES:
        in_plane = sorted(metrics_by_epsilon.get(float(epsilon), []), key=lambda item: item.solid_mode)
        for mode_index in range(1, N_ANALYTIC_MODES + 1):
            analytic = lookup[(float(epsilon), mode_index)]
            omega_c = float(analytic["Omega_cutoff"])
            if mode_index <= len(in_plane):
                metric = in_plane[mode_index - 1]
                diff_eb = metric.lambda_fem - float(analytic["Lambda_EB"])
                diff_timo = metric.lambda_fem - float(analytic["Lambda_Timoshenko"])
                closer = "Timoshenko" if abs(diff_timo) < abs(diff_eb) else "Euler-Bernoulli"
                rows.append(
                    {
                        **analytic,
                        "classified_in_plane_index": mode_index,
                        "solid_mode": metric.solid_mode,
                        "Omega_FEM": metric.omega,
                        "Lambda_FEM": metric.lambda_fem,
                        "Omega_FEM_over_cutoff": metric.omega / omega_c,
                        "below_cutoff_FEM": metric.omega < omega_c,
                        "FEM_minus_EB": diff_eb,
                        "FEM_minus_Timoshenko": diff_timo,
                        "closer_to": closer,
                        "classification": metric.classification,
                        "status": "classified_in_plane",
                        "notes": metric.classification_note,
                    }
                )
            else:
                rows.append(
                    {
                        **analytic,
                        "classified_in_plane_index": mode_index,
                        "solid_mode": "",
                        "Omega_FEM": "",
                        "Lambda_FEM": "",
                        "Omega_FEM_over_cutoff": "",
                        "below_cutoff_FEM": "",
                        "FEM_minus_EB": "",
                        "FEM_minus_Timoshenko": "",
                        "closer_to": "",
                        "classification": "",
                        "status": "not_classified",
                        "notes": "not enough classified in-plane bending-like solid modes in requested solid spectrum",
                    }
                )
    return rows


def write_csv(path: Path, rows: Sequence[dict[str, object]], fieldnames: Sequence[str]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_sorted_csv(rows: Sequence[dict[str, object]]) -> None:
    write_csv(
        OUTPUT_SORTED_CSV,
        rows,
        [
            "epsilon",
            "diameter_to_segment_length",
            "thin_rod_valid",
            "mode",
            "Lambda_EB",
            "Lambda_Timoshenko",
            "Omega_EB",
            "Omega_Timoshenko",
            "Omega_cutoff",
            "Omega_EB_over_cutoff",
            "Omega_Timoshenko_over_cutoff",
            "Lambda_cutoff",
            "below_cutoff_Timoshenko",
            "solid_sorted_mode",
            "Omega_FEM",
            "Lambda_FEM",
            "Omega_FEM_over_cutoff",
            "below_cutoff_FEM",
            "status",
            "notes",
        ],
    )


def write_mode_metrics_csv(rows: Sequence[ModeMetric]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "epsilon",
        "solid_mode",
        "Omega_FEM",
        "Lambda_FEM",
        "parsed_node_count",
        "in_plane_energy_fraction",
        "out_of_plane_energy_fraction",
        "axial_energy_fraction",
        "transverse_in_plane_energy_fraction",
        "mean_axial_energy_fraction",
        "mean_in_plane_bending_fraction",
        "mean_out_of_plane_bending_fraction",
        "torsion_indicator",
        "joint_motion_fraction",
        "warping_indicator",
        "classification",
        "classification_note",
    ]
    with OUTPUT_MODE_METRICS_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "epsilon": row.epsilon,
                    "solid_mode": row.solid_mode,
                    "Omega_FEM": row.omega,
                    "Lambda_FEM": row.lambda_fem,
                    "parsed_node_count": row.parsed_node_count,
                    "in_plane_energy_fraction": row.in_plane_energy_fraction,
                    "out_of_plane_energy_fraction": row.out_of_plane_energy_fraction,
                    "axial_energy_fraction": row.axial_energy_fraction,
                    "transverse_in_plane_energy_fraction": row.transverse_in_plane_energy_fraction,
                    "mean_axial_energy_fraction": row.mean_axial_energy_fraction,
                    "mean_in_plane_bending_fraction": row.mean_in_plane_bending_fraction,
                    "mean_out_of_plane_bending_fraction": row.mean_out_of_plane_bending_fraction,
                    "torsion_indicator": row.torsion_indicator,
                    "joint_motion_fraction": row.joint_motion_fraction,
                    "warping_indicator": row.warping_indicator,
                    "classification": row.classification,
                    "classification_note": row.classification_note,
                }
            )


def write_in_plane_csv(rows: Sequence[dict[str, object]]) -> None:
    write_csv(
        OUTPUT_IN_PLANE_CSV,
        rows,
        [
            "epsilon",
            "diameter_to_segment_length",
            "thin_rod_valid",
            "mode",
            "Lambda_EB",
            "Lambda_Timoshenko",
            "Omega_EB",
            "Omega_Timoshenko",
            "Omega_cutoff",
            "Omega_EB_over_cutoff",
            "Omega_Timoshenko_over_cutoff",
            "Lambda_cutoff",
            "below_cutoff_Timoshenko",
            "classified_in_plane_index",
            "solid_mode",
            "Omega_FEM",
            "Lambda_FEM",
            "Omega_FEM_over_cutoff",
            "below_cutoff_FEM",
            "FEM_minus_EB",
            "FEM_minus_Timoshenko",
            "closer_to",
            "classification",
            "status",
            "notes",
        ],
    )


def process_cases(
    audit: ToolAudit,
) -> tuple[
    dict[float, list[float]],
    list[ModeMetric],
    list[ModeShapeParseResult],
    list[MeshSummary],
    list[CcxInputSummary],
    list[CalculixResult],
    list[str],
    str,
]:
    fem_omegas_by_epsilon: dict[float, list[float]] = {}
    mode_metrics: list[ModeMetric] = []
    parse_results: list[ModeShapeParseResult] = []
    mesh_summaries: list[MeshSummary] = []
    ccx_inputs: list[CcxInputSummary] = []
    calculix_results: list[CalculixResult] = []
    messages: list[str] = []
    geometry_label = "not generated"
    for epsilon in EPSILON_VALUES:
        paths = case_paths(float(epsilon))
        if not audit.gmsh_exe:
            write_gmsh_geo(float(epsilon), paths, boolean_mode="union")
            messages.append(f"epsilon={float(epsilon):g}: no Gmsh executable detected; wrote .geo only.")
            continue
        mesh_ok, mesh_message, geometry_label, gmsh_messages = generate_mesh(paths, float(epsilon), audit.gmsh_exe)
        messages.append(f"epsilon={float(epsilon):g}: {mesh_message}")
        if gmsh_messages:
            messages.extend(f"epsilon={float(epsilon):g}: {item}" for item in gmsh_messages)
        if not mesh_ok or not paths.gmsh_inp.exists():
            continue
        mesh = read_gmsh_inp_mesh_data(paths.gmsh_inp)
        mesh_summaries.append(mesh_summary(float(epsilon), mesh, gmsh_messages))
        ccx_input = write_calculix_inputs(paths, float(epsilon), mesh)
        ccx_inputs.append(ccx_input)
        messages.append(
            f"epsilon={float(epsilon):g}: generated ROD1_OUTER_FIXED nodes={ccx_input.rod1_fixed_node_count}, "
            f"ROD2_OUTER_FIXED nodes={ccx_input.rod2_fixed_node_count}"
        )
        if not audit.ccx_exe:
            messages.append(f"epsilon={float(epsilon):g}: CalculiX ccx not found; modal solve not run.")
            continue
        result = run_calculix(paths, audit.ccx_exe, float(epsilon))
        calculix_results.append(result)
        messages.append(f"epsilon={float(epsilon):g}: {result.message}")
        if result.success and result.parsed_omegas:
            fem_omegas_by_epsilon[float(epsilon)] = list(result.parsed_omegas)
            metrics, parse_result = compute_mode_metrics(float(epsilon), paths, mesh, result.parsed_omegas)
            mode_metrics.extend(metrics)
            parse_results.append(parse_result)
            messages.append(f"epsilon={float(epsilon):g}: {parse_result.message}")
            if metrics:
                messages.append(f"epsilon={float(epsilon):g}: computed {len(metrics)} mode-shape metric rows.")
    return (
        fem_omegas_by_epsilon,
        mode_metrics,
        parse_results,
        mesh_summaries,
        ccx_inputs,
        calculix_results,
        messages,
        geometry_label,
    )


def fmt(value: object, digits: int = 8) -> str:
    if value == "" or value is None:
        return "-"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return "-"
    return f"{number:.{digits}g}"


def rel_path(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def markdown_mesh_table(rows: Sequence[MeshSummary]) -> str:
    if not rows:
        return "No mesh summaries were available."
    lines = [
        "| epsilon | nodes | solid elems | bbox min | bbox max | components | largest component | connected |",
        "| ---: | ---: | ---: | --- | --- | ---: | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            "| "
            f"{row.epsilon:g} | {row.node_count} | {row.solid_element_count} | "
            f"{row.bbox_min} | {row.bbox_max} | {row.connected_component_count} | "
            f"{row.largest_component_element_count} | {row.geometry_connected} |"
        )
    return "\n".join(lines)


def markdown_fixed_set_table(rows: Sequence[CcxInputSummary]) -> str:
    if not rows:
        return "No CalculiX input node-set summaries were available."
    lines = [
        "| epsilon | ROD1_OUTER_FIXED nodes | ROD2_OUTER_FIXED nodes | node sets present | modal input |",
        "| ---: | ---: | ---: | --- | --- |",
    ]
    for row in rows:
        lines.append(
            "| "
            f"{row.epsilon:g} | {row.rod1_fixed_node_count} | {row.rod2_fixed_node_count} | "
            f"{row.node_sets_present} | `{rel_path(row.modal_input_path)}` |"
        )
    return "\n".join(lines)


def markdown_calculix_lines(rows: Sequence[CalculixResult]) -> list[str]:
    if not rows:
        return ["CalculiX was not run or no CalculiX result objects were produced."]
    lines: list[str] = []
    for row in rows:
        status = "succeeded" if row.success else "failed"
        lines.append(f"- epsilon={row.epsilon:g}: {status}; {row.message}")
        lines.append(f"  stdout: `{rel_path(row.stdout_path)}`")
        lines.append(f"  stderr: `{rel_path(row.stderr_path)}`")
        lines.append(f"  dat: `{rel_path(row.dat_path)}`")
        lines.append(f"  frd: `{rel_path(row.frd_path)}`")
        lines.append(f"  parse source: {row.parse_source}")
    return lines


def markdown_mode_shape_lines(rows: Sequence[ModeShapeParseResult]) -> list[str]:
    if not rows:
        return ["Mode-shape parsing was not attempted or no successful CalculiX `.frd` results were available."]
    lines: list[str] = []
    for row in rows:
        status = "succeeded" if row.success else "failed"
        lines.append(f"- epsilon={row.epsilon:g}: {status}; {row.message}")
        lines.append(f"  source: {row.source}")
    return lines


def markdown_classification_counts(rows: Sequence[ModeMetric]) -> list[str]:
    if not rows:
        return ["- no mode-shape metric rows"]
    grouped: dict[float, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in rows:
        grouped[row.epsilon][row.classification] += 1
    lines: list[str] = []
    for epsilon in sorted(grouped):
        parts = ", ".join(f"{key}: {value}" for key, value in sorted(grouped[epsilon].items()))
        lines.append(f"- epsilon={epsilon:g}: {parts}")
    return lines


def markdown_metric_table(rows: Sequence[ModeMetric]) -> str:
    if not rows:
        return "No mode-shape metrics were computed."
    lines = [
        "| epsilon | solid mode | Lambda_FEM | class | in-plane | out-of-plane | axial | mean in-plane bend | mean out bend | torsion | joint | warp |",
        "| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in sorted(rows, key=lambda item: (item.epsilon, item.solid_mode)):
        lines.append(
            "| "
            f"{row.epsilon:g} | {row.solid_mode} | {row.lambda_fem:.8g} | {row.classification} | "
            f"{row.in_plane_energy_fraction:.3f} | {row.out_of_plane_energy_fraction:.3f} | "
            f"{row.axial_energy_fraction:.3f} | {row.mean_in_plane_bending_fraction:.3f} | "
            f"{row.mean_out_of_plane_bending_fraction:.3f} | {row.torsion_indicator:.3f} | "
            f"{row.joint_motion_fraction:.3f} | {row.warping_indicator:.3f} |"
        )
    return "\n".join(lines)


def markdown_in_plane_table(rows: Sequence[dict[str, object]]) -> str:
    if not rows:
        return "No classified in-plane comparison rows were available."
    lines = [
        "| epsilon | analytic sorted mode | solid mode | Lambda_EB | Lambda_Timoshenko | Lambda_FEM | FEM-EB | FEM-Timo | closer | cutoff FEM | status |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            "| "
            f"{fmt(row['epsilon'])} | {int(row['mode'])} | {fmt(row['solid_mode'])} | "
            f"{fmt(row['Lambda_EB'])} | {fmt(row['Lambda_Timoshenko'])} | {fmt(row['Lambda_FEM'])} | "
            f"{fmt(row['FEM_minus_EB'], 4)} | {fmt(row['FEM_minus_Timoshenko'], 4)} | "
            f"{row['closer_to'] if row['closer_to'] else '-'} | {fmt(row['Omega_FEM_over_cutoff'], 4)} | "
            f"{row['status']} |"
        )
    return "\n".join(lines)


def markdown_sorted_table(rows: Sequence[dict[str, object]]) -> str:
    if not rows:
        return "No sorted comparison rows were available."
    lines = [
        "| epsilon | analytic sorted mode | solid sorted mode | Lambda_EB | Lambda_Timoshenko | Lambda_FEM | cutoff FEM | status |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            "| "
            f"{fmt(row['epsilon'])} | {int(row['mode'])} | {fmt(row['solid_sorted_mode'])} | "
            f"{fmt(row['Lambda_EB'])} | {fmt(row['Lambda_Timoshenko'])} | {fmt(row['Lambda_FEM'])} | "
            f"{fmt(row['Omega_FEM_over_cutoff'], 4)} | {row['status']} |"
        )
    return "\n".join(lines)


def cutoff_lines(rows: Sequence[dict[str, object]]) -> list[str]:
    ratios: list[float] = []
    for row in rows:
        for key in ("Omega_EB_over_cutoff", "Omega_Timoshenko_over_cutoff", "Omega_FEM_over_cutoff"):
            value = row.get(key, "")
            if value != "":
                try:
                    number = float(value)
                except (TypeError, ValueError):
                    continue
                if math.isfinite(number):
                    ratios.append(number)
    if not ratios:
        return ["- no finite cut-off ratios were available"]
    max_ratio = max(ratios)
    lines = [f"- maximum finite Omega/Omega_c ratio in written comparison rows: {max_ratio:.6g}"]
    if max_ratio >= 1.0:
        lines.append("- ERROR: at least one row reaches or exceeds the Timoshenko cut-off.")
    elif max_ratio >= CUTOFF_WARNING_RATIO:
        lines.append("- WARNING: at least one row is within the 0.8 cut-off warning band.")
    else:
        lines.append("- all written analytic/FEM rows remain below the 0.8 cut-off warning band.")
    return lines


def closeness_lines(rows: Sequence[dict[str, object]]) -> list[str]:
    classified = [row for row in rows if row["status"] == "classified_in_plane"]
    if not classified:
        return ["- no classified in-plane FEM rows were available for closeness comparison"]
    timo_count = sum(1 for row in classified if row["closer_to"] == "Timoshenko")
    eb_count = sum(1 for row in classified if row["closer_to"] == "Euler-Bernoulli")
    missing = [row for row in rows if row["status"] != "classified_in_plane"]
    lines = [f"- {timo_count}/{len(classified)} classified in-plane rows are closer to Timoshenko; {eb_count} are closer to EB."]
    if timo_count == 0 and eb_count > 0:
        lines.append(
            "- this diagnostic run does not confirm the expected Timoshenko-closer trend; "
            "the finite fused 3D joint, solid-vs-point-joint mismatch, mesh resolution, "
            "normalization, and mode-shape identity should be reviewed before interpretation."
        )
    if missing:
        labels = ", ".join(f"epsilon={row['epsilon']:g} mode={int(row['mode'])}" for row in missing)
        lines.append(f"- no classified in-plane FEM mode was available for: {labels}.")
    return lines


def write_report(
    *,
    audit: ToolAudit,
    geometry_label: str,
    analytic_warnings: Sequence[str],
    mesh_summaries: Sequence[MeshSummary],
    ccx_inputs: Sequence[CcxInputSummary],
    calculix_results: Sequence[CalculixResult],
    parse_results: Sequence[ModeShapeParseResult],
    mode_metrics: Sequence[ModeMetric],
    sorted_rows: Sequence[dict[str, object]],
    in_plane_rows: Sequence[dict[str, object]],
    messages: Sequence[str],
) -> None:
    geom = geometry()
    lines: list[str] = []
    lines.append("# Coupled Equal Rods beta=15 3D Solid FEM Diagnostic")
    lines.append("")
    lines.append("Diagnostic only. This is a sorted-frequency and classified-mode comparison, not descendant branch tracking.")
    lines.append("No article file, paper_dorofeev_style file, article figure, old determinant, formulas.py, old solver, existing FEM physical model, or baseline result was changed.")
    lines.append("")
    lines.append("## Tools")
    lines.append("")
    lines.append(f"- Gmsh executable: `{audit.gmsh_exe}`" if audit.gmsh_exe else "- Gmsh executable: not found")
    lines.append(f"- Gmsh resolution source: {audit.gmsh_source if audit.gmsh_source else 'not found'}")
    lines.append(f"- Gmsh version: {audit.gmsh_version if audit.gmsh_version else 'not available'}")
    lines.append(f"- CalculiX executable: `{audit.ccx_exe}`" if audit.ccx_exe else "- CalculiX executable: not found")
    lines.append(f"- CalculiX resolution source: {audit.ccx_source if audit.ccx_source else 'not found'}")
    lines.append(f"- CalculiX version probe: {audit.ccx_version if audit.ccx_version else 'not available'}")
    lines.append("")
    lines.append("## Geometry And Boundary Conditions")
    lines.append("")
    lines.append(f"- beta: `{BETA_DEG:g} deg`")
    lines.append(f"- mu: `{MU:g}`, eta: `{ETA:g}`, tau1=tau2=`1`")
    lines.append(f"- joint: `{geom.joint}`")
    lines.append(f"- rod 1 outer end: `{geom.rod1.outer}`")
    lines.append(f"- rod 2 outer end: `{geom.rod2.outer}`")
    lines.append(f"- joint geometry used: `{geometry_label}`")
    lines.append("- fixed node sets: `ROD1_OUTER_FIXED`, `ROD2_OUTER_FIXED`")
    lines.append("- the joint is not fixed; it is part of the fused 3D solid geometry")
    lines.append("")
    lines.append("The 3D solid joint has finite volume and is not identical to the ideal point joint used in the 1D analytic model. This is a diagnostic validation path, not an article-level equivalence claim.")
    lines.append("")
    lines.append("## Normalization And Applicability")
    lines.append("")
    lines.append("- equal-rod project normalization: `Omega = epsilon*Lambda^2`")
    lines.append("- FEM conversion: `Lambda_FEM = sqrt(Omega/epsilon)`")
    lines.append("- if CalculiX `.dat` gives eigenvalue `Omega^2`, then `Omega = sqrt(eigenvalue)`")
    lines.append("- diameter-to-segment-length ratio: `2r/L_segment = 4*epsilon`")
    lines.append(f"- thin-rod criterion: `4*epsilon <= {THICKNESS_RATIO_LIMIT:g}`")
    lines.append("")
    invalid_eps = [epsilon for epsilon in EPSILON_VALUES if not thin_rod_valid(float(epsilon))]
    if invalid_eps:
        lines.append(f"Thin-rod warning: epsilon values outside the current criterion: {', '.join(f'{value:g}' for value in invalid_eps)}.")
    else:
        lines.append("All tested epsilon values satisfy the current thin-rod criterion.")
    lines.append("")
    lines.append("## Mesh Summary")
    lines.append("")
    lines.append(markdown_mesh_table(mesh_summaries))
    lines.append("")
    if any(summary.gmsh_messages for summary in mesh_summaries):
        lines.append("Captured Gmsh warning/error lines:")
        for summary in mesh_summaries:
            for message in summary.gmsh_messages:
                lines.append(f"- epsilon={summary.epsilon:g}: {message}")
        lines.append("")
    lines.append("## Fixed Node Sets")
    lines.append("")
    lines.append(markdown_fixed_set_table(ccx_inputs))
    lines.append("")
    lines.append("## CalculiX Runs")
    lines.append("")
    lines.extend(markdown_calculix_lines(calculix_results))
    lines.append("")
    parse_sources = sorted({result.parse_source for result in calculix_results if result.parsed_omegas})
    if parse_sources:
        lines.append("Recognized frequency convention:")
        for source in parse_sources:
            lines.append(f"- {source}")
        lines.append("")
    lines.append("## Analytic References")
    lines.append("")
    lines.append("Analytic EB/Timoshenko references are imported from `scripts/analysis/compare_coupled_equal_rods_eb_timoshenko.py`; the old determinant, formulas.py, and old solvers are used read-only and are not edited.")
    if analytic_warnings:
        lines.append("")
        lines.append("Analytic/root-search warnings:")
        for warning in analytic_warnings[:20]:
            lines.append(f"- {warning}")
        if len(analytic_warnings) > 20:
            lines.append(f"- ... {len(analytic_warnings) - 20} more warnings")
    lines.append("")
    lines.append("## Mode-Shape Classification")
    lines.append("")
    lines.append("Mode-shape parsing uses CalculiX `.frd` `DISP` blocks from `*NODE FILE U`. Classification is based on approximate nodal metrics in local rod frames: in-plane/out-of-plane energy, tangent-axis energy, centroidal slice motion, torsion-like circumferential fit, joint energy fraction, and within-slice warping.")
    lines.append("")
    lines.extend(markdown_mode_shape_lines(parse_results))
    lines.append("")
    lines.append("Classification counts:")
    lines.append("")
    lines.extend(markdown_classification_counts(mode_metrics))
    lines.append("")
    lines.append(markdown_metric_table(mode_metrics))
    lines.append("")
    lines.append("## Classified In-Plane Comparison")
    lines.append("")
    lines.append("This is the preferred comparison for the planar 1D coupled-rods model. It uses the first classified `in_plane_bending_like` solid modes rather than blindly comparing sorted solid modes.")
    lines.append("")
    lines.extend(closeness_lines(in_plane_rows))
    lines.append("")
    lines.append(markdown_in_plane_table(in_plane_rows))
    lines.append("")
    lines.append("## Preliminary Sorted Comparison")
    lines.append("")
    lines.append("The sorted comparison is retained only as a frequency-order diagnostic. It should not be used as a mode-identity claim.")
    lines.append("")
    lines.append(markdown_sorted_table(sorted_rows))
    lines.append("")
    lines.append("## Cut-Off Check")
    lines.append("")
    lines.append("Cut-off formula used:")
    lines.append("")
    lines.append("```text")
    lines.append("Omega_c = sqrt(kappa*G*A/(rho*I))")
    lines.append("```")
    lines.append("")
    lines.extend(cutoff_lines(list(sorted_rows) + list(in_plane_rows)))
    lines.append("")
    lines.append("## Case Messages")
    lines.append("")
    for message in messages:
        lines.append(f"- {message}")
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    lines.append(f"- sorted comparison CSV: `{rel_path(OUTPUT_SORTED_CSV)}`")
    lines.append(f"- mode metrics CSV: `{rel_path(OUTPUT_MODE_METRICS_CSV)}`")
    lines.append(f"- in-plane comparison CSV: `{rel_path(OUTPUT_IN_PLANE_CSV)}`")
    lines.append(f"- report: `{rel_path(OUTPUT_REPORT)}`")
    lines.append(f"- mesh/input/output directory: `{rel_path(SOLID_OUTPUT_DIR)}`")
    lines.append("")
    lines.append("## Next Step")
    lines.append("")
    lines.append("Inspect representative in-plane and out-of-plane mode shapes visually or add MAC-like solid/1D shape checks before making stronger validation claims.")
    lines.append("")
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    audit = audit_tools()
    analytic_rows, analytic_warnings = compute_analytic_rows()
    (
        fem_omegas_by_epsilon,
        mode_metrics,
        parse_results,
        mesh_summaries,
        ccx_inputs,
        calculix_results,
        messages,
        geometry_label,
    ) = process_cases(audit)
    sorted_rows = build_sorted_comparison_rows(analytic_rows, fem_omegas_by_epsilon)
    in_plane_rows = build_in_plane_comparison_rows(analytic_rows, mode_metrics)
    write_sorted_csv(sorted_rows)
    write_mode_metrics_csv(mode_metrics)
    write_in_plane_csv(in_plane_rows)
    write_report(
        audit=audit,
        geometry_label=geometry_label,
        analytic_warnings=analytic_warnings,
        mesh_summaries=mesh_summaries,
        ccx_inputs=ccx_inputs,
        calculix_results=calculix_results,
        parse_results=parse_results,
        mode_metrics=mode_metrics,
        sorted_rows=sorted_rows,
        in_plane_rows=in_plane_rows,
        messages=messages,
    )
    print(f"wrote sorted comparison CSV: {OUTPUT_SORTED_CSV}")
    print(f"wrote mode metrics CSV: {OUTPUT_MODE_METRICS_CSV}")
    print(f"wrote in-plane comparison CSV: {OUTPUT_IN_PLANE_CSV}")
    print(f"wrote report: {OUTPUT_REPORT}")
    print(f"wrote mesh/input/output files under: {SOLID_OUTPUT_DIR}")
    if not any(result.success for result in calculix_results):
        print("no 3D solid modal results computed; see report for missing tool or solver failure details")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
