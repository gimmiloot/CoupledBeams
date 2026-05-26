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

import numpy as np


# Diagnostic-only 3D solid FEM workflow for two equal coupled circular rods
# connected by an idealized CalculiX rigid reference-node point joint.
#
# This script is intentionally separate from the one-rod solid FEM workflow
# because the geometry, boundary-condition sets, mode classification basis, and
# comparison contract are different. It is also intentionally separate from the
# fused-cylinder coupled-rods workflow because the joint idealization and input
# contract are different. It does not make Gmsh or CalculiX required project
# dependencies.

BETA_DEG_VALUES = [15.0, 45.0, 90.0]
BETA_DEG = BETA_DEG_VALUES[0]
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
N_SLICES_PER_ROD = 80
CUTOFF_WARNING_RATIO = 0.8
MAC_STRONG_THRESHOLD = 0.8
MAC_MODERATE_THRESHOLD = 0.5

MU = 0.0
ETA = 0.0
TAU1 = 1.0
TAU2 = 1.0

REPO_ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"
OUTPUT_DIR = REPO_ROOT / "results"
SOLID_OUTPUT_DIR = OUTPUT_DIR / "solid_fem_coupled_equal_rods_point_joint"
OUTPUT_GEOMETRY_AUDIT_CSV = OUTPUT_DIR / "coupled_equal_rods_point_joint_3d_solid_fem_geometry_audit.csv"
OUTPUT_REPORT = OUTPUT_DIR / "coupled_equal_rods_point_joint_3d_solid_fem_report.md"
OUTPUT_IN_PLANE_CSV = OUTPUT_DIR / "coupled_equal_rods_point_joint_3d_solid_fem_in_plane_comparison_by_beta.csv"
OUTPUT_MODE_METRICS_CSV = OUTPUT_DIR / "coupled_equal_rods_point_joint_3d_solid_fem_mode_metrics_by_beta.csv"
OUTPUT_SORTED_CSV = OUTPUT_DIR / "coupled_equal_rods_point_joint_3d_solid_fem_sorted_comparison_by_beta.csv"
OUTPUT_MAC_MATRIX_EB_CSV = OUTPUT_DIR / "coupled_equal_rods_point_joint_3d_solid_fem_mac_matrix_eb.csv"
OUTPUT_MAC_MATRIX_TIMO_CSV = OUTPUT_DIR / "coupled_equal_rods_point_joint_3d_solid_fem_mac_matrix_timoshenko.csv"
OUTPUT_MAC_MATCHES_CSV = OUTPUT_DIR / "coupled_equal_rods_point_joint_3d_solid_fem_mac_matches.csv"
FUSED_IN_PLANE_CSV = OUTPUT_DIR / "coupled_equal_rods_3d_solid_fem_in_plane_comparison_by_beta.csv"

JOINT_REF_NODE_NAME = "JOINT_REF"
JOINT_COUPLED_SET_NAME = "JOINT_COUPLED_ALL"
RIGID_CONSTRAINT_METHOD = "*RIGID BODY, NSET=JOINT_COUPLED_ALL, REF NODE=<joint reference node>"

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
class ConstraintCapabilityAudit:
    supported: bool
    selected_method: str
    message: str
    local_search_root: str | None
    local_documentation_hits: tuple[str, ...]
    probe_input_path: Path | None
    probe_stdout_path: Path | None
    probe_stderr_path: Path | None
    probe_dat_path: Path | None
    probe_returncode: int | None


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
    beta_deg: float
    epsilon: float
    radius: float
    sin_beta: float
    estimated_overlap_length: float
    estimated_overlap_fraction: float
    node_count: int
    solid_element_count: int
    bbox_min: tuple[float, float, float] | None
    bbox_max: tuple[float, float, float] | None
    connected_component_count: int
    largest_component_element_count: int
    geometry_connected: bool
    joint_zone_node_count: int
    joint_zone_element_count: int
    gmsh_messages: tuple[str, ...]


@dataclass(frozen=True)
class CcxInputSummary:
    beta_deg: float
    epsilon: float
    modal_input_path: Path
    mesh_include_path: Path
    rod1_fixed_node_count: int
    rod2_fixed_node_count: int
    rod1_inner_node_count: int
    rod2_inner_node_count: int
    joint_coupled_node_count: int
    joint_ref_node_id: int
    constraint_method: str
    node_sets_present: bool


@dataclass(frozen=True)
class CalculixResult:
    beta_deg: float
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
    beta_deg: float
    epsilon: float
    success: bool
    message: str
    frd_path: Path
    mode_shapes: dict[int, dict[int, tuple[float, float, float]]]
    source: str


@dataclass(frozen=True)
class ModeMetric:
    beta_deg: float
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


@dataclass(frozen=True)
class SolidCenterlineShape:
    beta_deg: float
    epsilon: float
    solid_mode: int
    omega: float
    lambda_fem: float
    classification: str
    vector: np.ndarray
    filled_slice_count: int
    empty_slice_count: int


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


def local_constraint_documentation_hits(ccx_exe: str | None) -> tuple[str | None, tuple[str, ...]]:
    roots: list[Path] = []
    if ccx_exe is not None:
        roots.append(Path(ccx_exe).resolve().parent)
    roots.extend(Path(root) for root in CCX_SEARCH_ROOTS)
    search_root: Path | None = None
    hits: list[str] = []
    patterns = ("*RIGID BODY", "*KINEMATIC COUPLING", "*MPC", "*EQUATION")
    suffixes = {".inp", ".txt", ".htm", ".html", ".fbd", ".dat", ""}
    for root in roots:
        if not root.exists() or search_root is not None:
            continue
        search_root = root
        for path in root.rglob("*"):
            if len(hits) >= 40:
                break
            if not path.is_file() or path.suffix.lower() not in suffixes:
                continue
            try:
                if path.stat().st_size > 5_000_000:
                    continue
                text = path.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            upper_text = text.upper()
            if any(pattern in upper_text for pattern in patterns):
                hits.append(str(path))
    return (str(search_root) if search_root is not None else None), tuple(sorted(set(hits)))


def write_constraint_probe_input(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                "** Diagnostic-only capability probe for one rigid reference joint.",
                "** Two solid top-face sets are tied to one shared reference node.",
                "*NODE",
                "1, -1, 0, 0",
                "2, 0, 0, 0",
                "3, 0, 1, 0",
                "4, -1, 1, 0",
                "5, -1, 0, 1",
                "6, 0, 0, 1",
                "7, 0, 1, 1",
                "8, -1, 1, 1",
                "11, 0, 0, 0",
                "12, 1, 0, 0",
                "13, 1, 1, 0",
                "14, 0, 1, 0",
                "15, 0, 0, 1",
                "16, 1, 0, 1",
                "17, 1, 1, 1",
                "18, 0, 1, 1",
                "1000, 0, 0.5, 1",
                "*ELEMENT, TYPE=C3D8, ELSET=EALL",
                "1, 1,2,3,4,5,6,7,8",
                "2, 11,12,13,14,15,16,17,18",
                "*NSET, NSET=BOTTOM",
                "1,2,3,4,11,12,13,14",
                "*NSET, NSET=JOINT_COUPLED_ALL",
                "5,6,7,8,15,16,17,18",
                "*NSET, NSET=JOINT_REF",
                "1000",
                "*MATERIAL, NAME=MAT",
                "*ELASTIC",
                "1.0, 0.3",
                "*DENSITY",
                "1.0",
                "*SOLID SECTION, ELSET=EALL, MATERIAL=MAT",
                "*RIGID BODY, NSET=JOINT_COUPLED_ALL, REF NODE=1000",
                "*BOUNDARY",
                "BOTTOM, 1, 3, 0",
                "*STEP",
                "*FREQUENCY",
                "3",
                "*NODE FILE",
                "U",
                "*END STEP",
            ]
        )
        + "\n",
        encoding="ascii",
    )


def audit_constraint_capability(audit: ToolAudit) -> ConstraintCapabilityAudit:
    search_root, hits = local_constraint_documentation_hits(audit.ccx_exe)
    method = RIGID_CONSTRAINT_METHOD
    if audit.ccx_exe is None:
        return ConstraintCapabilityAudit(
            supported=False,
            selected_method=method,
            message="CalculiX executable is not available; rigid point-joint capability was not tested.",
            local_search_root=search_root,
            local_documentation_hits=hits,
            probe_input_path=None,
            probe_stdout_path=None,
            probe_stderr_path=None,
            probe_dat_path=None,
            probe_returncode=None,
        )

    probe_dir = SOLID_OUTPUT_DIR / "_constraint_capability_probe"
    probe_input = probe_dir / "rigid_body_point_joint_probe.inp"
    probe_stdout = probe_dir / "rigid_body_point_joint_probe.stdout.txt"
    probe_stderr = probe_dir / "rigid_body_point_joint_probe.stderr.txt"
    probe_dat = probe_dir / "rigid_body_point_joint_probe.dat"
    write_constraint_probe_input(probe_input)
    for stale_path in (probe_stdout, probe_stderr, probe_dat, probe_dir / "rigid_body_point_joint_probe.frd"):
        try:
            stale_path.unlink(missing_ok=True)
        except OSError:
            pass
    try:
        completed = subprocess.run(
            [audit.ccx_exe, probe_input.stem],
            cwd=str(probe_dir),
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=120,
        )
    except Exception as exc:
        probe_stdout.write_text("", encoding="utf-8")
        probe_stderr.write_text(str(exc), encoding="utf-8")
        return ConstraintCapabilityAudit(
            supported=False,
            selected_method=method,
            message=f"Rigid-body capability probe failed before completion: {exc}",
            local_search_root=search_root,
            local_documentation_hits=hits,
            probe_input_path=probe_input,
            probe_stdout_path=probe_stdout,
            probe_stderr_path=probe_stderr,
            probe_dat_path=probe_dat,
            probe_returncode=None,
        )
    probe_stdout.write_text(completed.stdout, encoding="utf-8")
    probe_stderr.write_text(completed.stderr, encoding="utf-8")
    combined = f"{completed.stdout}\n{completed.stderr}"
    supported = completed.returncode == 0 and "Job finished" in combined
    message = (
        "CalculiX accepted a shared-reference-node *RIGID BODY probe and completed the modal step."
        if supported
        else f"CalculiX rejected or failed the shared-reference-node *RIGID BODY probe with return code {completed.returncode}."
    )
    return ConstraintCapabilityAudit(
        supported=supported,
        selected_method=method,
        message=message,
        local_search_root=search_root,
        local_documentation_hits=hits,
        probe_input_path=probe_input,
        probe_stdout_path=probe_stdout,
        probe_stderr_path=probe_stderr,
        probe_dat_path=probe_dat,
        probe_returncode=completed.returncode,
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


def geometry_for_beta(beta_deg: float) -> Geometry:
    beta = math.radians(float(beta_deg))
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


def geometry() -> Geometry:
    return geometry_for_beta(BETA_DEG)


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


def beta_token(beta_deg: float) -> str:
    return safe_float_token(float(beta_deg))


def overlap_estimates(beta_deg: float, epsilon: float) -> tuple[float, float, float, float]:
    radius = radius_from_epsilon(float(epsilon))
    sin_beta = abs(math.sin(math.radians(float(beta_deg))))
    if sin_beta <= 1.0e-14:
        overlap_length = float("inf")
    else:
        overlap_length = radius / sin_beta
    return radius, sin_beta, overlap_length, overlap_length / L_SEGMENT


def case_paths(beta_deg: float, epsilon: float) -> CasePaths:
    beta = beta_token(float(beta_deg))
    token = safe_float_token(float(epsilon))
    case_dir = SOLID_OUTPUT_DIR / f"beta_{beta}" / f"eps_{token}"
    stem = f"coupled_equal_rods_beta{beta}_eps{token}"
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


def write_gmsh_geo(epsilon: float, paths: CasePaths) -> str:
    radius = radius_from_epsilon(float(epsilon))
    mesh_size = mesh_size_for_epsilon(float(epsilon))
    geom = geometry()
    cos_b = math.cos(geom.beta_rad)
    sin_b = math.sin(geom.beta_rad)
    paths.case_dir.mkdir(parents=True, exist_ok=True)
    geometry_label = "point-joint constraint model: separate cylinders plus CalculiX rigid reference-node coupling"
    bbox = L_SEGMENT + 2.0 * radius
    paths.geo.write_text(
        f"""// Diagnostic-only two equal coupled rods solid mesh with separate cylinder volumes.
// beta = {BETA_DEG:.17g} deg, epsilon = {float(epsilon):.17g}
SetFactory("OpenCASCADE");

L = {L_SEGMENT:.17g};
R = {radius:.17g};
h = {mesh_size:.17g};
CosB = {cos_b:.17g};
SinB = {sin_b:.17g};

Cylinder(1) = {{-L, 0, 0, L, 0, 0, R, 2*Pi}};
Cylinder(2) = {{-L*CosB, -L*SinB, 0, L*CosB, L*SinB, 0, R, 2*Pi}};
// Deliberately no BooleanUnion/BooleanFragments: the cylinders remain
// separate solid volumes and are connected only by CalculiX constraints.
Physical Volume("ROD1_SOLID", 1) = {{1}};
Physical Volume("ROD2_SOLID", 2) = {{2}};

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
    geometry_label = write_gmsh_geo(float(epsilon), paths)
    ok, message, messages = generate_mesh_with_gmsh_cli(paths, gmsh_exe)
    return ok, message, geometry_label, messages


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


def mesh_summary(beta_deg: float, epsilon: float, mesh: MeshData, gmsh_messages: Sequence[str]) -> MeshSummary:
    components, largest = connected_component_counts(mesh)
    radius, sin_beta, overlap_length, overlap_fraction = overlap_estimates(float(beta_deg), float(epsilon))
    geom = geometry()
    joint_zone_radius = 2.0 * radius
    joint_zone_nodes = {
        node_id for node_id, point in mesh.nodes.items() if norm(sub(point, geom.joint)) <= joint_zone_radius
    }
    joint_zone_element_count = sum(
        1 for connectivity in mesh.solid_elements.values() if any(node_id in joint_zone_nodes for node_id in connectivity)
    )
    return MeshSummary(
        beta_deg=float(beta_deg),
        epsilon=float(epsilon),
        radius=radius,
        sin_beta=sin_beta,
        estimated_overlap_length=overlap_length,
        estimated_overlap_fraction=overlap_fraction,
        node_count=len(mesh.nodes),
        solid_element_count=len(mesh.solid_elements),
        bbox_min=mesh.bbox_min,
        bbox_max=mesh.bbox_max,
        connected_component_count=components,
        largest_component_element_count=largest,
        geometry_connected=components == 1,
        joint_zone_node_count=len(joint_zone_nodes),
        joint_zone_element_count=joint_zone_element_count,
        gmsh_messages=tuple(gmsh_messages),
    )


def projection_on_rod(point: tuple[float, float, float], rod: RodFrame) -> tuple[float, float]:
    rel = sub(point, rod.outer)
    s = dot(rel, rod.tangent)
    axis_point = add(rod.outer, mul(s, rod.tangent))
    radial = norm(sub(point, axis_point))
    return s, radial


def rod_assignment_cost(point: tuple[float, float, float], rod: RodFrame) -> float:
    s, radial = projection_on_rod(point, rod)
    outside = max(0.0, -s, s - L_SEGMENT)
    return radial + outside


def rod_node_memberships(mesh: MeshData) -> tuple[dict[int, set[int]], dict[int, int]]:
    geom = geometry()
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


def coordinate_case_node_sets(mesh: MeshData, epsilon: float) -> tuple[list[int], list[int], list[int], list[int]]:
    geom = geometry()
    memberships, _element_counts = rod_node_memberships(mesh)
    radius = radius_from_epsilon(float(epsilon))
    plane_tol = max(mesh_size_for_epsilon(float(epsilon)) / 50.0, 1.0e-6)
    radial_tol = 1.02 * radius + plane_tol
    rod1_nodes: list[int] = []
    rod2_nodes: list[int] = []
    rod1_inner_nodes: list[int] = []
    rod2_inner_nodes: list[int] = []

    for node_id in memberships[1]:
        point = mesh.nodes.get(node_id)
        if point is None:
            continue
        s, radial = projection_on_rod(point, geom.rod1)
        if abs(s) <= plane_tol and radial <= radial_tol:
            rod1_nodes.append(node_id)
        if abs(s - L_SEGMENT) <= plane_tol and radial <= radial_tol:
            rod1_inner_nodes.append(node_id)
    for node_id in memberships[2]:
        point = mesh.nodes.get(node_id)
        if point is None:
            continue
        s, radial = projection_on_rod(point, geom.rod2)
        if abs(s) <= plane_tol and radial <= radial_tol:
            rod2_nodes.append(node_id)
        if abs(s - L_SEGMENT) <= plane_tol and radial <= radial_tol:
            rod2_inner_nodes.append(node_id)
    return sorted(rod1_nodes), sorted(rod2_nodes), sorted(rod1_inner_nodes), sorted(rod2_inner_nodes)


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


def write_calculix_inputs(paths: CasePaths, beta_deg: float, epsilon: float, mesh: MeshData) -> CcxInputSummary:
    rod1_nodes, rod2_nodes, rod1_inner_nodes, rod2_inner_nodes = coordinate_case_node_sets(mesh, float(epsilon))
    joint_coupled_nodes = sorted(set(rod1_inner_nodes) | set(rod2_inner_nodes))
    joint_ref_node_id = max(mesh.nodes.keys(), default=0) + 1
    write_calculix_mesh_include(paths, mesh)
    paths.ccx_modal_inp.write_text(
        "\n".join(
            [
                "** Diagnostic-only CalculiX modal input for two separate equal circular rods.",
                "** Outer fixed and inner coupled node sets are generated from local rod coordinate projections.",
                "** The inner end faces are tied to one shared rigid reference node.",
                f"*INCLUDE, INPUT={paths.ccx_mesh_inp.name}",
                "*NODE",
                f"{joint_ref_node_id}, 0.0, 0.0, 0.0",
                f"*NSET, NSET={JOINT_REF_NODE_NAME}",
                f"{joint_ref_node_id}",
                "*NSET, NSET=ROD1_OUTER_FIXED",
                *format_id_lines(rod1_nodes),
                "*NSET, NSET=ROD2_OUTER_FIXED",
                *format_id_lines(rod2_nodes),
                "*NSET, NSET=ROD1_INNER_COUPLED",
                *format_id_lines(rod1_inner_nodes),
                "*NSET, NSET=ROD2_INNER_COUPLED",
                *format_id_lines(rod2_inner_nodes),
                f"*NSET, NSET={JOINT_COUPLED_SET_NAME}",
                *format_id_lines(joint_coupled_nodes),
                "*MATERIAL, NAME=MAT",
                "*ELASTIC",
                f"{E:.17g}, {NU:.17g}",
                "*DENSITY",
                f"{RHO:.17g}",
                "*SOLID SECTION, ELSET=SOLID, MATERIAL=MAT",
                f"*RIGID BODY, NSET={JOINT_COUPLED_SET_NAME}, REF NODE={joint_ref_node_id}",
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
        beta_deg=float(beta_deg),
        epsilon=float(epsilon),
        modal_input_path=paths.ccx_modal_inp,
        mesh_include_path=paths.ccx_mesh_inp,
        rod1_fixed_node_count=len(rod1_nodes),
        rod2_fixed_node_count=len(rod2_nodes),
        rod1_inner_node_count=len(rod1_inner_nodes),
        rod2_inner_node_count=len(rod2_inner_nodes),
        joint_coupled_node_count=len(joint_coupled_nodes),
        joint_ref_node_id=joint_ref_node_id,
        constraint_method=RIGID_CONSTRAINT_METHOD.replace("<joint reference node>", str(joint_ref_node_id)),
        node_sets_present=bool(rod1_nodes and rod2_nodes and rod1_inner_nodes and rod2_inner_nodes),
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


def run_calculix(paths: CasePaths, ccx: str, beta_deg: float, epsilon: float) -> CalculixResult:
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
            beta_deg=float(beta_deg),
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
            beta_deg=float(beta_deg),
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
        beta_deg=float(beta_deg),
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


def parse_calculix_frd_mode_shapes(paths: CasePaths, beta_deg: float, epsilon: float) -> ModeShapeParseResult:
    if not paths.ccx_frd.exists():
        return ModeShapeParseResult(
            beta_deg=float(beta_deg),
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
            beta_deg=float(beta_deg),
            epsilon=float(epsilon),
            success=False,
            message="No displacement eigenvector blocks were recognized in the CalculiX .frd file.",
            frd_path=paths.ccx_frd,
            mode_shapes={},
            source="frd_no_displacement_blocks",
        )
    return ModeShapeParseResult(
        beta_deg=float(beta_deg),
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
    beta_deg: float,
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
        beta_deg=float(beta_deg),
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
    beta_deg: float,
    epsilon: float,
    paths: CasePaths,
    mesh: MeshData,
    omegas: Sequence[float],
) -> tuple[list[ModeMetric], ModeShapeParseResult]:
    parse_result = parse_calculix_frd_mode_shapes(paths, float(beta_deg), float(epsilon))
    if not parse_result.success:
        return [], parse_result
    metrics: list[ModeMetric] = []
    for mode_index, omega in enumerate(omegas, start=1):
        shape = parse_result.mode_shapes.get(mode_index)
        if shape is None:
            continue
        metric = compute_mode_metric(float(beta_deg), float(epsilon), mode_index, float(omega), mesh, shape)
        if metric is not None:
            metrics.append(metric)
    return metrics, parse_result


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
    mesh: MeshData,
    mode_shape: dict[int, tuple[float, float, float]],
) -> tuple[np.ndarray, int, int]:
    geom = geometry()
    memberships, _element_counts = rod_node_memberships(mesh)
    rod_vectors: list[np.ndarray] = []
    total_filled = 0
    total_empty = 0
    for rod_index, rod in ((1, geom.rod1), (2, geom.rod2)):
        values = np.zeros((N_SLICES_PER_ROD, 2), dtype=float)
        counts = np.zeros(N_SLICES_PER_ROD, dtype=float)
        for node_id in memberships[rod_index]:
            displacement = mode_shape.get(node_id)
            point = mesh.nodes.get(node_id)
            if displacement is None or point is None:
                continue
            s, _radial = projection_on_rod(point, rod)
            if s < -1.0e-8 or s > L_SEGMENT + 1.0e-8:
                continue
            slice_index = int(math.floor(max(0.0, min(0.999999999999, s / L_SEGMENT)) * N_SLICES_PER_ROD))
            values[slice_index, 0] += float(displacement[0])
            values[slice_index, 1] += float(displacement[1])
            counts[slice_index] += 1.0
        means, filled, empty = fill_missing_slice_means(values, counts)
        rod_vectors.append(means.reshape(-1))
        total_filled += filled
        total_empty += empty
    return np.concatenate(rod_vectors), total_filled, total_empty


def compute_solid_centerline_shapes(
    beta_deg: float,
    epsilon: float,
    mesh: MeshData,
    omegas: Sequence[float],
    parse_result: ModeShapeParseResult,
    metrics: Sequence[ModeMetric],
) -> list[SolidCenterlineShape]:
    class_by_mode = {metric.solid_mode: metric.classification for metric in metrics}
    rows: list[SolidCenterlineShape] = []
    for mode_index, omega in enumerate(omegas, start=1):
        mode_shape = parse_result.mode_shapes.get(mode_index)
        if mode_shape is None:
            continue
        vector, filled, empty = solid_centerline_vector(mesh, mode_shape)
        rows.append(
            SolidCenterlineShape(
                beta_deg=float(beta_deg),
                epsilon=float(epsilon),
                solid_mode=int(mode_index),
                omega=float(omega),
                lambda_fem=lambda_fem_from_omega(float(omega), float(epsilon)),
                classification=class_by_mode.get(mode_index, "not_classified"),
                vector=vector,
                filled_slice_count=filled,
                empty_slice_count=empty,
            )
        )
    return rows


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
    for beta_deg in BETA_DEG_VALUES:
        reference.BETA_DEG = float(beta_deg)
        for epsilon in EPSILON_VALUES:
            eb_values = reference.eb_roots(float(epsilon))
            timo_values, root_warnings = reference.timo_roots(float(epsilon), eb_values)
            warnings.extend(f"beta={float(beta_deg):g}, epsilon={float(epsilon):g}: {warning}" for warning in root_warnings)
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
                        "beta_deg": float(beta_deg),
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


def analytic_by_case_mode(rows: Sequence[dict[str, object]]) -> dict[tuple[float, float, int], dict[str, object]]:
    return {(float(row["beta_deg"]), float(row["epsilon"]), int(row["mode"])): row for row in rows}


def row_normalized_np(matrix: np.ndarray) -> np.ndarray:
    out = np.asarray(matrix, dtype=float).copy()
    norms = np.linalg.norm(out, axis=1)
    for index, row_norm in enumerate(norms):
        if row_norm > 0.0 and math.isfinite(float(row_norm)):
            out[index, :] /= row_norm
    return out


def analytic_null_vector_from_matrix(matrix: np.ndarray) -> tuple[np.ndarray, float, float]:
    scaled = row_normalized_np(matrix)
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


def slice_coordinates() -> np.ndarray:
    return (np.arange(N_SLICES_PER_ROD, dtype=float) + 0.5) / float(N_SLICES_PER_ROD)


def interleaved_global_vector(
    beta_deg: float,
    u_left: np.ndarray,
    w_left: np.ndarray,
    u_right: np.ndarray,
    w_right: np.ndarray,
) -> np.ndarray:
    geom = geometry_for_beta(float(beta_deg))
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


def eb_centerline_vector(
    reference: ModuleType,
    beta_deg: float,
    epsilon: float,
    lambda_value: float,
) -> tuple[np.ndarray | None, str]:
    if not math.isfinite(float(lambda_value)):
        return None, "missing EB Lambda"
    matrix = reference.assemble_clamped_coupled_matrix(
        float(lambda_value), math.radians(float(beta_deg)), MU, float(epsilon)
    )
    coeff, smallest, ratio = analytic_null_vector_from_matrix(matrix)
    a1, b1, a2, b2, p1, p2 = [float(value) for value in coeff]
    xi = slice_coordinates()
    z1 = float(lambda_value) * xi
    theta1 = float(epsilon) * float(lambda_value) ** 2 * xi
    z2 = -float(lambda_value) * xi
    theta2 = -float(epsilon) * float(lambda_value) ** 2 * xi
    w_left = a1 * (np.cos(z1) - np.cosh(z1)) + b1 * (np.sin(z1) - np.sinh(z1))
    u_left = p1 * np.sin(theta1)
    w_right = a2 * (np.cos(z2) - np.cosh(z2)) + b2 * (np.sin(z2) - np.sinh(z2))
    u_right = p2 * np.sin(theta2)
    vector = interleaved_global_vector(float(beta_deg), u_left, w_left, u_right, w_right)
    return vector, f"EB null vector residual smallest_singular={smallest:.3e}, ratio={ratio:.3e}"


def timo_centerline_vector(
    reference: ModuleType,
    beta_deg: float,
    epsilon: float,
    lambda_value: float,
) -> tuple[np.ndarray | None, str]:
    if not math.isfinite(float(lambda_value)):
        return None, "missing Timoshenko Lambda"
    old_beta = getattr(reference, "BETA_DEG", None)
    reference.BETA_DEG = float(beta_deg)
    try:
        matrix, basis_warnings = reference.timo_coupling_matrix(float(lambda_value), float(epsilon))
        coeff, smallest, ratio = analytic_null_vector_from_matrix(matrix)
        section = reference.section_from_epsilon(float(epsilon))
        basis = reference.timo_basis(float(lambda_value), float(epsilon), section)
    except Exception as exc:
        if old_beta is not None:
            reference.BETA_DEG = old_beta
        return None, f"Timoshenko shape reconstruction failed: {exc}"
    finally:
        if old_beta is not None:
            reference.BETA_DEG = old_beta
    a1, b1, p1, a2, b2, p2 = [float(value) for value in coeff]
    xi = slice_coordinates()
    theta = project_omega(float(lambda_value), float(epsilon))
    w_left = np.zeros_like(xi)
    w_right = np.zeros_like(xi)
    for index, s_value in enumerate(xi):
        left_cols = reference.bending_endpoint_columns(float(s_value), basis)["w"]
        right_cols = reference.bending_endpoint_columns(-float(s_value), basis)["w"]
        w_left[index] = a1 * float(left_cols[0]) + b1 * float(left_cols[1])
        w_right[index] = a2 * float(right_cols[0]) + b2 * float(right_cols[1])
    u_left = p1 * np.sin(theta * xi)
    u_right = p2 * np.sin(-theta * xi)
    vector = interleaved_global_vector(float(beta_deg), u_left, w_left, u_right, w_right)
    warning_text = "; ".join(basis_warnings) if basis_warnings else "no Timoshenko basis warnings"
    return vector, f"Timoshenko null vector residual smallest_singular={smallest:.3e}, ratio={ratio:.3e}; {warning_text}"


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


def build_analytic_shape_vectors(
    analytic_rows: Sequence[dict[str, object]],
) -> tuple[dict[tuple[float, float, str, int], tuple[np.ndarray, str]], list[str]]:
    reference = load_analytic_module()
    vectors: dict[tuple[float, float, str, int], tuple[np.ndarray, str]] = {}
    warnings: list[str] = []
    for row in analytic_rows:
        beta_deg = float(row["beta_deg"])
        epsilon = float(row["epsilon"])
        mode = int(row["mode"])
        for model, lambda_key in (("EB", "Lambda_EB"), ("Timoshenko", "Lambda_Timoshenko")):
            lambda_value = float(row[lambda_key])
            if model == "EB":
                vector, note = eb_centerline_vector(reference, beta_deg, epsilon, lambda_value)
            else:
                vector, note = timo_centerline_vector(reference, beta_deg, epsilon, lambda_value)
            if vector is None:
                warnings.append(f"beta={beta_deg:g}, epsilon={epsilon:g}, {model} mode={mode}: {note}")
                continue
            vectors[(beta_deg, epsilon, model, mode)] = (vector, note)
    return vectors, sorted(set(warnings))


def build_mac_matrix_rows(
    analytic_rows: Sequence[dict[str, object]],
    solid_shapes: Sequence[SolidCenterlineShape],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str]]:
    analytic_vectors, warnings = build_analytic_shape_vectors(analytic_rows)
    analytic_lookup = analytic_by_case_mode(analytic_rows)
    solid_by_case: dict[tuple[float, float], list[SolidCenterlineShape]] = defaultdict(list)
    for solid in solid_shapes:
        solid_by_case[(solid.beta_deg, solid.epsilon)].append(solid)
    rows_by_model: dict[str, list[dict[str, object]]] = {"EB": [], "Timoshenko": []}
    for beta_deg in BETA_DEG_VALUES:
        for epsilon in EPSILON_VALUES:
            solids = sorted(solid_by_case.get((float(beta_deg), float(epsilon)), []), key=lambda item: item.solid_mode)
            for mode in range(1, N_ANALYTIC_MODES + 1):
                analytic = analytic_lookup[(float(beta_deg), float(epsilon), mode)]
                for model, lambda_key in (("EB", "Lambda_EB"), ("Timoshenko", "Lambda_Timoshenko")):
                    vector_note = analytic_vectors.get((float(beta_deg), float(epsilon), model, mode))
                    if vector_note is None:
                        continue
                    analytic_vector, note = vector_note
                    lambda_analytic = float(analytic[lambda_key])
                    omega_analytic = project_omega(lambda_analytic, float(epsilon))
                    for solid in solids:
                        frequency_abs_diff = abs(float(solid.omega) - omega_analytic)
                        frequency_rel_diff = (
                            frequency_abs_diff / abs(omega_analytic) if abs(omega_analytic) > 1.0e-14 else float("nan")
                        )
                        rows_by_model[model].append(
                            {
                                "beta": float(beta_deg),
                                "epsilon": float(epsilon),
                                "solid_mode": solid.solid_mode,
                                "solid_class": solid.classification,
                                "analytic_model": model,
                                "analytic_mode": mode,
                                "Lambda_solid": solid.lambda_fem,
                                "Lambda_analytic": lambda_analytic,
                                "MAC": centerline_mac(solid.vector, analytic_vector),
                                "frequency_abs_diff": frequency_abs_diff,
                                "frequency_rel_diff": frequency_rel_diff,
                                "analytic_shape_note": note,
                            }
                        )
    return rows_by_model["EB"], rows_by_model["Timoshenko"], warnings


def candidate_mac_rows(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    preferred = [
        row
        for row in rows
        if row["solid_class"] in {"in_plane_bending_like", "mixed_or_unclassified"}
        and math.isfinite(float(row["MAC"]))
    ]
    if preferred:
        return preferred
    return [row for row in rows if math.isfinite(float(row["MAC"]))]


def best_row_by_mac(rows: Sequence[dict[str, object]]) -> dict[str, object] | None:
    candidates = candidate_mac_rows(rows)
    if not candidates:
        return None
    return max(candidates, key=lambda row: float(row["MAC"]))


def build_mac_match_rows(
    eb_rows: Sequence[dict[str, object]],
    timo_rows: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    all_rows = list(eb_rows) + list(timo_rows)
    by_solid: dict[tuple[float, float, int], dict[str, list[dict[str, object]]]] = defaultdict(lambda: defaultdict(list))
    by_analytic: dict[tuple[float, float, str, int], list[dict[str, object]]] = defaultdict(list)
    for row in all_rows:
        beta = float(row["beta"])
        epsilon = float(row["epsilon"])
        model = str(row["analytic_model"])
        solid_mode = int(row["solid_mode"])
        analytic_mode = int(row["analytic_mode"])
        by_solid[(beta, epsilon, solid_mode)][model].append(row)
        by_analytic[(beta, epsilon, model, analytic_mode)].append(row)

    output_rows: list[dict[str, object]] = []
    for key, model_rows in sorted(by_solid.items()):
        beta, epsilon, solid_mode = key
        best_eb = best_row_by_mac(model_rows.get("EB", []))
        best_timo = best_row_by_mac(model_rows.get("Timoshenko", []))
        solid_class = str((best_eb or best_timo or {}).get("solid_class", ""))
        lambda_solid = (best_eb or best_timo or {}).get("Lambda_solid", "")
        best_eb_mac = float(best_eb["MAC"]) if best_eb is not None else float("nan")
        best_timo_mac = float(best_timo["MAC"]) if best_timo is not None else float("nan")
        preferred_shape = "Timoshenko" if best_timo_mac > best_eb_mac else "EB"
        eb_freq = float(best_eb["frequency_abs_diff"]) if best_eb is not None else float("inf")
        timo_freq = float(best_timo["frequency_abs_diff"]) if best_timo is not None else float("inf")
        preferred_frequency = "Timoshenko" if timo_freq < eb_freq else "EB"
        notes = []
        if mac_strength(max(best_eb_mac, best_timo_mac)) in {"weak", "missing"}:
            notes.append("weak_best_MAC")
        output_rows.append(
            {
                "match_direction": "solid_to_analytic",
                "beta": beta,
                "epsilon": epsilon,
                "solid_mode": solid_mode,
                "solid_class": solid_class,
                "Lambda_solid": lambda_solid,
                "best_EB_mode": int(best_eb["analytic_mode"]) if best_eb is not None else "",
                "best_EB_MAC": best_eb_mac if math.isfinite(best_eb_mac) else "",
                "best_EB_Lambda": best_eb["Lambda_analytic"] if best_eb is not None else "",
                "best_EB_abs_freq_diff": eb_freq if math.isfinite(eb_freq) else "",
                "best_Timo_mode": int(best_timo["analytic_mode"]) if best_timo is not None else "",
                "best_Timo_MAC": best_timo_mac if math.isfinite(best_timo_mac) else "",
                "best_Timo_Lambda": best_timo["Lambda_analytic"] if best_timo is not None else "",
                "best_Timo_abs_freq_diff": timo_freq if math.isfinite(timo_freq) else "",
                "preferred_shape_model_by_MAC": preferred_shape,
                "preferred_frequency_model_after_MAC": preferred_frequency,
                "notes": ";".join(notes),
                "analytic_model": "",
                "analytic_mode": "",
                "Lambda_analytic": "",
                "best_solid_mode": "",
                "best_solid_class": "",
                "best_MAC": "",
                "best_MAC_strength": "",
                "abs_freq_diff": "",
                "rel_freq_diff": "",
                "duplicate_match_warning": "",
            }
        )

    analytic_best_rows: list[dict[str, object]] = []
    for key, rows in sorted(by_analytic.items()):
        beta, epsilon, model, analytic_mode = key
        best = best_row_by_mac(rows)
        if best is None:
            continue
        analytic_best_rows.append(best)

    duplicate_counts: dict[tuple[float, float, str, int], int] = defaultdict(int)
    for row in analytic_best_rows:
        duplicate_counts[(float(row["beta"]), float(row["epsilon"]), str(row["analytic_model"]), int(row["solid_mode"]))] += 1

    for best in analytic_best_rows:
        duplicate_key = (
            float(best["beta"]),
            float(best["epsilon"]),
            str(best["analytic_model"]),
            int(best["solid_mode"]),
        )
        duplicate = duplicate_counts[duplicate_key] > 1
        output_rows.append(
            {
                "match_direction": "analytic_to_solid",
                "beta": best["beta"],
                "epsilon": best["epsilon"],
                "solid_mode": "",
                "solid_class": "",
                "Lambda_solid": "",
                "best_EB_mode": "",
                "best_EB_MAC": "",
                "best_EB_Lambda": "",
                "best_EB_abs_freq_diff": "",
                "best_Timo_mode": "",
                "best_Timo_MAC": "",
                "best_Timo_Lambda": "",
                "best_Timo_abs_freq_diff": "",
                "preferred_shape_model_by_MAC": "",
                "preferred_frequency_model_after_MAC": "",
                "notes": "duplicate_best_solid_mode" if duplicate else "",
                "analytic_model": best["analytic_model"],
                "analytic_mode": best["analytic_mode"],
                "Lambda_analytic": best["Lambda_analytic"],
                "best_solid_mode": best["solid_mode"],
                "best_solid_class": best["solid_class"],
                "best_MAC": best["MAC"],
                "best_MAC_strength": mac_strength(best["MAC"]),
                "abs_freq_diff": best["frequency_abs_diff"],
                "rel_freq_diff": best["frequency_rel_diff"],
                "duplicate_match_warning": "duplicate" if duplicate else "",
            }
        )

    analytic_summary_lookup: dict[tuple[float, float, int], dict[str, dict[str, object]]] = defaultdict(dict)
    for best in analytic_best_rows:
        analytic_summary_lookup[(float(best["beta"]), float(best["epsilon"]), int(best["analytic_mode"]))][
            str(best["analytic_model"])
        ] = best
    for key, model_best in sorted(analytic_summary_lookup.items()):
        beta, epsilon, analytic_mode = key
        eb = model_best.get("EB")
        timo = model_best.get("Timoshenko")
        if eb is None and timo is None:
            continue
        eb_mac = float(eb["MAC"]) if eb is not None else float("nan")
        timo_mac = float(timo["MAC"]) if timo is not None else float("nan")
        preferred_shape = "Timoshenko" if timo_mac > eb_mac else "EB"
        shape_row = timo if preferred_shape == "Timoshenko" else eb
        if shape_row is None:
            shape_row = eb or timo
        lambda_solid = float(shape_row["Lambda_solid"])
        lambda_eb = float(eb["Lambda_analytic"]) if eb is not None else float("nan")
        lambda_timo = float(timo["Lambda_analytic"]) if timo is not None else float("nan")
        omega_solid = project_omega(lambda_solid, epsilon)
        eb_freq = abs(omega_solid - project_omega(lambda_eb, epsilon)) if math.isfinite(lambda_eb) else float("inf")
        timo_freq = abs(omega_solid - project_omega(lambda_timo, epsilon)) if math.isfinite(lambda_timo) else float("inf")
        preferred_frequency = "Timoshenko" if timo_freq < eb_freq else "EB"
        output_rows.append(
            {
                "match_direction": "analytic_mode_summary",
                "beta": beta,
                "epsilon": epsilon,
                "solid_mode": "",
                "solid_class": "",
                "Lambda_solid": lambda_solid,
                "best_EB_mode": analytic_mode,
                "best_EB_MAC": eb_mac if math.isfinite(eb_mac) else "",
                "best_EB_Lambda": lambda_eb if math.isfinite(lambda_eb) else "",
                "best_EB_abs_freq_diff": eb_freq if math.isfinite(eb_freq) else "",
                "best_Timo_mode": analytic_mode,
                "best_Timo_MAC": timo_mac if math.isfinite(timo_mac) else "",
                "best_Timo_Lambda": lambda_timo if math.isfinite(lambda_timo) else "",
                "best_Timo_abs_freq_diff": timo_freq if math.isfinite(timo_freq) else "",
                "preferred_shape_model_by_MAC": preferred_shape,
                "preferred_frequency_model_after_MAC": preferred_frequency,
                "notes": "preferred over simple first-in-plane ordering",
                "analytic_model": "EB_vs_Timoshenko",
                "analytic_mode": analytic_mode,
                "Lambda_analytic": "",
                "best_solid_mode": shape_row["solid_mode"],
                "best_solid_class": shape_row["solid_class"],
                "best_MAC": shape_row["MAC"],
                "best_MAC_strength": mac_strength(shape_row["MAC"]),
                "abs_freq_diff": min(eb_freq, timo_freq),
                "rel_freq_diff": "",
                "duplicate_match_warning": "",
            }
        )
    return output_rows


def build_sorted_comparison_rows(
    analytic_rows: Sequence[dict[str, object]],
    fem_omegas_by_case: dict[tuple[float, float], list[float]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for analytic in analytic_rows:
        beta_deg = float(analytic["beta_deg"])
        epsilon = float(analytic["epsilon"])
        mode = int(analytic["mode"])
        fem_omegas = fem_omegas_by_case.get((beta_deg, epsilon), [])
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
    lookup = analytic_by_case_mode(analytic_rows)
    rows: list[dict[str, object]] = []
    metrics_by_case: dict[tuple[float, float], list[ModeMetric]] = defaultdict(list)
    for metric in mode_metrics:
        if metric.classification == "in_plane_bending_like":
            metrics_by_case[(metric.beta_deg, metric.epsilon)].append(metric)
    for beta_deg in BETA_DEG_VALUES:
        for epsilon in EPSILON_VALUES:
            key = (float(beta_deg), float(epsilon))
            in_plane = sorted(metrics_by_case.get(key, []), key=lambda item: item.solid_mode)
            for mode_index in range(1, N_ANALYTIC_MODES + 1):
                analytic = lookup[(float(beta_deg), float(epsilon), mode_index)]
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
            "beta_deg",
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
        "beta_deg",
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
                    "beta_deg": row.beta_deg,
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
            "beta_deg",
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


def write_geometry_audit_csv(rows: Sequence[MeshSummary]) -> None:
    write_csv(
        OUTPUT_GEOMETRY_AUDIT_CSV,
        [
            {
                "beta_deg": row.beta_deg,
                "epsilon": row.epsilon,
                "radius": row.radius,
                "sin_beta": row.sin_beta,
                "estimated_overlap_length": row.estimated_overlap_length,
                "estimated_overlap_fraction": row.estimated_overlap_fraction,
                "node_count": row.node_count,
                "solid_element_count": row.solid_element_count,
                "bbox_min": row.bbox_min,
                "bbox_max": row.bbox_max,
                "connected_component_count": row.connected_component_count,
                "largest_component_element_count": row.largest_component_element_count,
                "geometry_connected": row.geometry_connected,
                "joint_zone_node_count": row.joint_zone_node_count,
                "joint_zone_element_count": row.joint_zone_element_count,
                "gmsh_messages": " | ".join(row.gmsh_messages),
            }
            for row in rows
        ],
        [
            "beta_deg",
            "epsilon",
            "radius",
            "sin_beta",
            "estimated_overlap_length",
            "estimated_overlap_fraction",
            "node_count",
            "solid_element_count",
            "bbox_min",
            "bbox_max",
            "connected_component_count",
            "largest_component_element_count",
            "geometry_connected",
            "joint_zone_node_count",
            "joint_zone_element_count",
            "gmsh_messages",
        ],
    )


def write_mac_matrix_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    write_csv(
        path,
        rows,
        [
            "beta",
            "epsilon",
            "solid_mode",
            "solid_class",
            "analytic_model",
            "analytic_mode",
            "Lambda_solid",
            "Lambda_analytic",
            "MAC",
            "frequency_abs_diff",
            "frequency_rel_diff",
            "analytic_shape_note",
        ],
    )


def write_mac_matches_csv(rows: Sequence[dict[str, object]]) -> None:
    write_csv(
        OUTPUT_MAC_MATCHES_CSV,
        rows,
        [
            "match_direction",
            "beta",
            "epsilon",
            "solid_mode",
            "solid_class",
            "Lambda_solid",
            "best_EB_mode",
            "best_EB_MAC",
            "best_EB_Lambda",
            "best_EB_abs_freq_diff",
            "best_Timo_mode",
            "best_Timo_MAC",
            "best_Timo_Lambda",
            "best_Timo_abs_freq_diff",
            "preferred_shape_model_by_MAC",
            "preferred_frequency_model_after_MAC",
            "notes",
            "analytic_model",
            "analytic_mode",
            "Lambda_analytic",
            "best_solid_mode",
            "best_solid_class",
            "best_MAC",
            "best_MAC_strength",
            "abs_freq_diff",
            "rel_freq_diff",
            "duplicate_match_warning",
        ],
    )


def process_cases(
    audit: ToolAudit,
    constraint_audit: ConstraintCapabilityAudit,
) -> tuple[
    dict[tuple[float, float], list[float]],
    list[ModeMetric],
    list[ModeShapeParseResult],
    list[SolidCenterlineShape],
    list[MeshSummary],
    list[CcxInputSummary],
    list[CalculixResult],
    list[str],
    str,
]:
    global BETA_DEG
    fem_omegas_by_case: dict[tuple[float, float], list[float]] = {}
    mode_metrics: list[ModeMetric] = []
    parse_results: list[ModeShapeParseResult] = []
    solid_centerline_shapes: list[SolidCenterlineShape] = []
    mesh_summaries: list[MeshSummary] = []
    ccx_inputs: list[CcxInputSummary] = []
    calculix_results: list[CalculixResult] = []
    messages: list[str] = []
    geometry_label = "not generated"
    for beta_deg in BETA_DEG_VALUES:
        BETA_DEG = float(beta_deg)
        for epsilon in EPSILON_VALUES:
            case_label = f"beta={float(beta_deg):g}, epsilon={float(epsilon):g}"
            paths = case_paths(float(beta_deg), float(epsilon))
            if not audit.gmsh_exe:
                write_gmsh_geo(float(epsilon), paths)
                messages.append(f"{case_label}: no Gmsh executable detected; wrote .geo only.")
                continue
            mesh_ok, mesh_message, geometry_label, gmsh_messages = generate_mesh(paths, float(epsilon), audit.gmsh_exe)
            messages.append(f"{case_label}: {mesh_message}")
            if gmsh_messages:
                messages.extend(f"{case_label}: {item}" for item in gmsh_messages)
            if not mesh_ok or not paths.gmsh_inp.exists():
                continue
            mesh = read_gmsh_inp_mesh_data(paths.gmsh_inp)
            mesh_summaries.append(mesh_summary(float(beta_deg), float(epsilon), mesh, gmsh_messages))
            ccx_input = write_calculix_inputs(paths, float(beta_deg), float(epsilon), mesh)
            ccx_inputs.append(ccx_input)
            messages.append(
                f"{case_label}: generated ROD1_OUTER_FIXED nodes={ccx_input.rod1_fixed_node_count}, "
                f"ROD2_OUTER_FIXED nodes={ccx_input.rod2_fixed_node_count}, "
                f"ROD1_INNER_COUPLED nodes={ccx_input.rod1_inner_node_count}, "
                f"ROD2_INNER_COUPLED nodes={ccx_input.rod2_inner_node_count}"
            )
            if not audit.ccx_exe:
                messages.append(f"{case_label}: CalculiX ccx not found; modal solve not run.")
                continue
            if not constraint_audit.supported:
                messages.append(
                    f"{case_label}: rigid point-joint constraint capability was not confirmed; modal solve not run."
                )
                continue
            result = run_calculix(paths, audit.ccx_exe, float(beta_deg), float(epsilon))
            calculix_results.append(result)
            messages.append(f"{case_label}: {result.message}")
            if result.success and result.parsed_omegas:
                fem_omegas_by_case[(float(beta_deg), float(epsilon))] = list(result.parsed_omegas)
                metrics, parse_result = compute_mode_metrics(
                    float(beta_deg), float(epsilon), paths, mesh, result.parsed_omegas
                )
                mode_metrics.extend(metrics)
                parse_results.append(parse_result)
                solid_centerline_shapes.extend(
                    compute_solid_centerline_shapes(
                        float(beta_deg), float(epsilon), mesh, result.parsed_omegas, parse_result, metrics
                    )
                )
                messages.append(f"{case_label}: {parse_result.message}")
                if metrics:
                    messages.append(f"{case_label}: computed {len(metrics)} mode-shape metric rows.")
    return (
        fem_omegas_by_case,
        mode_metrics,
        parse_results,
        solid_centerline_shapes,
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
        "| beta | epsilon | r | sin(beta) | r/sin(beta) | overlap frac | nodes | solid elems | bbox min | bbox max | components | connected | joint-zone nodes | joint-zone elems |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | ---: | --- | ---: | ---: |",
    ]
    for row in sorted(rows, key=lambda item: (item.beta_deg, item.epsilon)):
        lines.append(
            "| "
            f"{row.beta_deg:g} | {row.epsilon:g} | {row.radius:.6g} | {row.sin_beta:.6g} | "
            f"{row.estimated_overlap_length:.6g} | {row.estimated_overlap_fraction:.6g} | "
            f"{row.node_count} | {row.solid_element_count} | {row.bbox_min} | {row.bbox_max} | "
            f"{row.connected_component_count} | {row.geometry_connected} | "
            f"{row.joint_zone_node_count} | {row.joint_zone_element_count} |"
        )
    return "\n".join(lines)


def markdown_fixed_set_table(rows: Sequence[CcxInputSummary]) -> str:
    if not rows:
        return "No CalculiX input node-set summaries were available."
    lines = [
        "| beta | epsilon | ROD1_OUTER_FIXED | ROD2_OUTER_FIXED | ROD1_INNER_COUPLED | ROD2_INNER_COUPLED | joint coupled | joint ref | node sets present | modal input |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for row in sorted(rows, key=lambda item: (item.beta_deg, item.epsilon)):
        lines.append(
            "| "
            f"{row.beta_deg:g} | {row.epsilon:g} | {row.rod1_fixed_node_count} | {row.rod2_fixed_node_count} | "
            f"{row.rod1_inner_node_count} | {row.rod2_inner_node_count} | {row.joint_coupled_node_count} | "
            f"{row.joint_ref_node_id} | {row.node_sets_present} | `{rel_path(row.modal_input_path)}` |"
        )
    return "\n".join(lines)


def markdown_calculix_lines(rows: Sequence[CalculixResult]) -> list[str]:
    if not rows:
        return ["CalculiX was not run or no CalculiX result objects were produced."]
    lines: list[str] = []
    for row in sorted(rows, key=lambda item: (item.beta_deg, item.epsilon)):
        status = "succeeded" if row.success else "failed"
        lines.append(f"- beta={row.beta_deg:g}, epsilon={row.epsilon:g}: {status}; {row.message}")
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
    for row in sorted(rows, key=lambda item: (item.beta_deg, item.epsilon)):
        status = "succeeded" if row.success else "failed"
        lines.append(f"- beta={row.beta_deg:g}, epsilon={row.epsilon:g}: {status}; {row.message}")
        lines.append(f"  source: {row.source}")
    return lines


def markdown_classification_counts(rows: Sequence[ModeMetric]) -> list[str]:
    if not rows:
        return ["- no mode-shape metric rows"]
    grouped: dict[tuple[float, float], dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in rows:
        grouped[(row.beta_deg, row.epsilon)][row.classification] += 1
    lines: list[str] = []
    for beta_deg, epsilon in sorted(grouped):
        parts = ", ".join(f"{key}: {value}" for key, value in sorted(grouped[(beta_deg, epsilon)].items()))
        lines.append(f"- beta={beta_deg:g}, epsilon={epsilon:g}: {parts}")
    return lines


def markdown_metric_table(rows: Sequence[ModeMetric]) -> str:
    if not rows:
        return "No mode-shape metrics were computed."
    lines = [
        "| beta | epsilon | solid mode | Lambda_FEM | class | in-plane | out-of-plane | axial | mean in-plane bend | mean out bend | torsion | joint | warp |",
        "| ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in sorted(rows, key=lambda item: (item.beta_deg, item.epsilon, item.solid_mode)):
        lines.append(
            "| "
            f"{row.beta_deg:g} | {row.epsilon:g} | {row.solid_mode} | {row.lambda_fem:.8g} | {row.classification} | "
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
        "| beta | epsilon | analytic sorted mode | solid mode | Lambda_EB | Lambda_Timoshenko | Lambda_FEM | FEM-EB | FEM-Timo | closer | cutoff FEM | status |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            "| "
            f"{fmt(row['beta_deg'])} | {fmt(row['epsilon'])} | {int(row['mode'])} | {fmt(row['solid_mode'])} | "
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
        "| beta | epsilon | analytic sorted mode | solid sorted mode | Lambda_EB | Lambda_Timoshenko | Lambda_FEM | cutoff FEM | status |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            "| "
            f"{fmt(row['beta_deg'])} | {fmt(row['epsilon'])} | {int(row['mode'])} | {fmt(row['solid_sorted_mode'])} | "
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
    grouped: dict[float, list[dict[str, object]]] = defaultdict(list)
    for row in classified:
        grouped[float(row["beta_deg"])].append(row)
    mean_timo_by_beta: list[tuple[float, float]] = []
    for beta_deg in sorted(grouped):
        beta_rows = grouped[beta_deg]
        beta_timo_count = sum(1 for row in beta_rows if row["closer_to"] == "Timoshenko")
        beta_eb_count = sum(1 for row in beta_rows if row["closer_to"] == "Euler-Bernoulli")
        mean_abs_eb = sum(abs(float(row["FEM_minus_EB"])) for row in beta_rows) / len(beta_rows)
        mean_abs_timo = sum(abs(float(row["FEM_minus_Timoshenko"])) for row in beta_rows) / len(beta_rows)
        mean_timo_by_beta.append((beta_deg, mean_abs_timo))
        lines.append(
            f"- beta={beta_deg:g}: {beta_timo_count}/{len(beta_rows)} closer to Timoshenko, "
            f"{beta_eb_count} closer to EB; mean |FEM-EB|={mean_abs_eb:.6g}, "
            f"mean |FEM-Timoshenko|={mean_abs_timo:.6g}."
        )
    if len(mean_timo_by_beta) >= 2:
        values = [value for _beta, value in sorted(mean_timo_by_beta)]
        if all(later <= earlier for earlier, later in zip(values, values[1:])):
            lines.append("- mean absolute FEM-Timoshenko mismatch decreases as beta increases in this diagnostic sweep.")
        else:
            lines.append("- mean absolute FEM-Timoshenko mismatch is not monotone decreasing with beta in this diagnostic sweep.")
    if timo_count == 0 and eb_count > 0:
        lines.append(
            "- this diagnostic run does not confirm the expected Timoshenko-closer trend; "
            "the constraint idealization, solid-vs-1D model mismatch, mesh resolution, "
            "normalization, and mode-shape identity should be reviewed before interpretation."
        )
    if missing:
        labels = ", ".join(
            f"beta={row['beta_deg']:g} epsilon={row['epsilon']:g} mode={int(row['mode'])}" for row in missing
        )
        lines.append(f"- no classified in-plane FEM mode was available for: {labels}.")
    return lines


def fused_comparison_lines() -> list[str]:
    if not FUSED_IN_PLANE_CSV.exists():
        return ["- previous fused-joint by-beta comparison CSV was not found; no high-level comparison was made."]
    rows: list[dict[str, str]] = []
    with FUSED_IN_PLANE_CSV.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = [row for row in reader if row.get("status") == "classified_in_plane"]
    if not rows:
        return ["- previous fused-joint CSV contains no classified in-plane rows."]
    grouped: dict[float, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        try:
            grouped[float(row["beta_deg"])].append(row)
        except (KeyError, ValueError):
            continue
    lines = ["Previous fused-cylinder diagnostic, for context only:"]
    for beta_deg in sorted(grouped):
        beta_rows = grouped[beta_deg]
        timo_count = sum(1 for row in beta_rows if row.get("closer_to") == "Timoshenko")
        eb_count = sum(1 for row in beta_rows if row.get("closer_to") == "Euler-Bernoulli")
        try:
            mean_abs_timo = sum(abs(float(row["FEM_minus_Timoshenko"])) for row in beta_rows) / len(beta_rows)
            mean_abs_eb = sum(abs(float(row["FEM_minus_EB"])) for row in beta_rows) / len(beta_rows)
        except (KeyError, ValueError):
            mean_abs_timo = float("nan")
            mean_abs_eb = float("nan")
        lines.append(
            f"- fused beta={beta_deg:g}: {timo_count}/{len(beta_rows)} closer to Timoshenko, "
            f"{eb_count} closer to EB; mean |FEM-EB|={mean_abs_eb:.6g}, "
            f"mean |FEM-Timoshenko|={mean_abs_timo:.6g}."
        )
    return lines


def mac_match_summary_lines(rows: Sequence[dict[str, object]], warnings: Sequence[str]) -> list[str]:
    analytic_rows = [row for row in rows if row.get("match_direction") == "analytic_to_solid"]
    summary_rows = [row for row in rows if row.get("match_direction") == "analytic_mode_summary"]
    if not analytic_rows:
        lines = ["- no MAC-like analytic-to-solid match rows were available"]
    else:
        strength_counts: dict[str, int] = defaultdict(int)
        for row in analytic_rows:
            strength_counts[str(row.get("best_MAC_strength", "missing"))] += 1
        lines = [
            "- analytic-to-solid MAC strength counts: "
            + ", ".join(f"{key}: {strength_counts[key]}" for key in ("strong", "moderate", "weak", "missing"))
        ]
        duplicate_count = sum(1 for row in analytic_rows if row.get("duplicate_match_warning") == "duplicate")
        if duplicate_count:
            lines.append(f"- duplicate best-solid-mode warnings: {duplicate_count}")
        else:
            lines.append("- duplicate best-solid-mode warnings: none")
    if summary_rows:
        timo_shape = sum(1 for row in summary_rows if row.get("preferred_shape_model_by_MAC") == "Timoshenko")
        eb_shape = sum(1 for row in summary_rows if row.get("preferred_shape_model_by_MAC") == "EB")
        timo_freq = sum(1 for row in summary_rows if row.get("preferred_frequency_model_after_MAC") == "Timoshenko")
        eb_freq = sum(1 for row in summary_rows if row.get("preferred_frequency_model_after_MAC") == "EB")
        lines.append(
            f"- analytic-mode summaries by MAC: shape preference Timoshenko {timo_shape}/{len(summary_rows)}, "
            f"EB {eb_shape}/{len(summary_rows)}."
        )
        lines.append(
            f"- after MAC-selected mode identity, frequency is closer to Timoshenko {timo_freq}/{len(summary_rows)} "
            f"and closer to EB {eb_freq}/{len(summary_rows)}."
        )
        grouped: dict[float, list[dict[str, object]]] = defaultdict(list)
        for row in summary_rows:
            grouped[float(row["beta"])].append(row)
        for beta in sorted(grouped):
            beta_rows = grouped[beta]
            beta_timo_freq = sum(1 for row in beta_rows if row.get("preferred_frequency_model_after_MAC") == "Timoshenko")
            beta_eb_freq = sum(1 for row in beta_rows if row.get("preferred_frequency_model_after_MAC") == "EB")
            beta_strong = sum(1 for row in beta_rows if row.get("best_MAC_strength") == "strong")
            beta_moderate = sum(1 for row in beta_rows if row.get("best_MAC_strength") == "moderate")
            beta_weak = sum(1 for row in beta_rows if row.get("best_MAC_strength") == "weak")
            lines.append(
                f"- beta={beta:g}: frequency closer after MAC -> Timoshenko {beta_timo_freq}/{len(beta_rows)}, "
                f"EB {beta_eb_freq}; MAC strengths strong/moderate/weak = {beta_strong}/{beta_moderate}/{beta_weak}."
            )
    if warnings:
        lines.append(f"- analytic shape reconstruction warnings: {len(warnings)}")
    else:
        lines.append("- analytic shape reconstruction warnings: none")
    return lines


def markdown_mac_summary_table(rows: Sequence[dict[str, object]]) -> str:
    summary_rows = [row for row in rows if row.get("match_direction") == "analytic_mode_summary"]
    if not summary_rows:
        return "No MAC-like analytic-mode summary rows were available."
    lines = [
        "| beta | epsilon | analytic mode | best solid | solid class | MAC | strength | Lambda_EB | Lambda_Timo | Lambda_solid | shape pref | freq pref |",
        "| ---: | ---: | ---: | ---: | --- | ---: | --- | ---: | ---: | ---: | --- | --- |",
    ]
    for row in sorted(summary_rows, key=lambda item: (float(item["beta"]), float(item["epsilon"]), int(item["analytic_mode"]))):
        lines.append(
            "| "
            f"{fmt(row['beta'])} | {fmt(row['epsilon'])} | {int(row['analytic_mode'])} | "
            f"{fmt(row['best_solid_mode'])} | {row['best_solid_class']} | {fmt(row['best_MAC'], 4)} | "
            f"{row['best_MAC_strength']} | {fmt(row['best_EB_Lambda'])} | {fmt(row['best_Timo_Lambda'])} | "
            f"{fmt(row['Lambda_solid'])} | {row['preferred_shape_model_by_MAC']} | "
            f"{row['preferred_frequency_model_after_MAC']} |"
        )
    return "\n".join(lines)


def write_report(
    *,
    audit: ToolAudit,
    constraint_audit: ConstraintCapabilityAudit,
    geometry_label: str,
    analytic_warnings: Sequence[str],
    mesh_summaries: Sequence[MeshSummary],
    ccx_inputs: Sequence[CcxInputSummary],
    calculix_results: Sequence[CalculixResult],
    parse_results: Sequence[ModeShapeParseResult],
    mode_metrics: Sequence[ModeMetric],
    sorted_rows: Sequence[dict[str, object]],
    in_plane_rows: Sequence[dict[str, object]],
    mac_match_rows: Sequence[dict[str, object]],
    mac_warnings: Sequence[str],
    messages: Sequence[str],
) -> None:
    lines: list[str] = []
    lines.append("# Coupled Equal Rods Point-Joint 3D Solid FEM Diagnostic Sweep")
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
    lines.append("## Constraint Capability Audit")
    lines.append("")
    lines.append(f"- selected constraint method: `{constraint_audit.selected_method}`")
    lines.append(f"- supported by local probe: `{constraint_audit.supported}`")
    lines.append(f"- message: {constraint_audit.message}")
    lines.append(f"- local search root: `{constraint_audit.local_search_root}`" if constraint_audit.local_search_root else "- local search root: not available")
    if constraint_audit.local_documentation_hits:
        lines.append("- local documentation/example hits:")
        for path in constraint_audit.local_documentation_hits[:20]:
            lines.append(f"  - `{path}`")
        if len(constraint_audit.local_documentation_hits) > 20:
            lines.append(f"  - ... {len(constraint_audit.local_documentation_hits) - 20} more")
    else:
        lines.append("- local documentation/example hits: none found in the installed CalculiX folder")
    if constraint_audit.probe_input_path is not None:
        lines.append(f"- probe input: `{rel_path(constraint_audit.probe_input_path)}`")
        lines.append(f"- probe stdout: `{rel_path(constraint_audit.probe_stdout_path)}`")
        lines.append(f"- probe stderr: `{rel_path(constraint_audit.probe_stderr_path)}`")
        lines.append(f"- probe dat: `{rel_path(constraint_audit.probe_dat_path)}`")
        lines.append(f"- probe return code: `{constraint_audit.probe_returncode}`")
    lines.append("")
    lines.append("## Geometry And Boundary Conditions")
    lines.append("")
    lines.append(f"- beta values: `{', '.join(f'{value:g}' for value in BETA_DEG_VALUES)} deg`")
    lines.append(f"- mu: `{MU:g}`, eta: `{ETA:g}`, tau1=tau2=`1`")
    lines.append("- joint: `(0, 0, 0)`")
    lines.append("- rod 1 outer end: `(-L, 0, 0)`")
    lines.append("- rod 2 outer end: `(-L*cos(beta), -L*sin(beta), 0)`")
    lines.append(f"- joint geometry used: `{geometry_label}`")
    lines.append("- cylinders are separate solid volumes; no BooleanUnion, BooleanFragments, or fused overlap volume is used")
    lines.append("- fixed node sets: `ROD1_OUTER_FIXED`, `ROD2_OUTER_FIXED`")
    lines.append("- coupled inner node sets: `ROD1_INNER_COUPLED`, `ROD2_INNER_COUPLED`, combined as `JOINT_COUPLED_ALL`")
    lines.append("- shared reference set: `JOINT_REF`")
    lines.append("- the joint reference node is not fixed; both inner faces are rigidly coupled to it")
    lines.append("")
    lines.append("This point-joint constraint model is closer to the ideal 1D point joint than the fused-cylinder solid geometry, but it remains a diagnostic 3D idealization and not an article-level equivalence claim.")
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
    lines.append("Geometry audit keeps `estimated_overlap_length = r/sin(beta)` only as a comparison scale for the old fused-cylinder geometry. In this point-joint workflow no fused overlap volume is created; the separate meshes should have two geometric connected components before constraints.")
    lines.append("")
    lines.append(markdown_mesh_table(mesh_summaries))
    lines.append("")
    if any(summary.gmsh_messages for summary in mesh_summaries):
        lines.append("Captured Gmsh warning/error lines:")
        for summary in mesh_summaries:
            for message in summary.gmsh_messages:
                lines.append(f"- beta={summary.beta_deg:g}, epsilon={summary.epsilon:g}: {message}")
        lines.append("")
    lines.append("## Fixed Node Sets")
    lines.append("")
    lines.append(markdown_fixed_set_table(ccx_inputs))
    lines.append("")
    lines.append(f"Rigid coupling method used in production inputs: `{RIGID_CONSTRAINT_METHOD}`.")
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
    lines.append("## MAC-like Centerline Shape Matching")
    lines.append("")
    lines.append(
        "The MAC-like comparison is the preferred mode-identity diagnostic in this report. "
        f"For each solid mode, nodal displacements are averaged on `{N_SLICES_PER_ROD}` axial slices per rod. "
        "The planar centerline vector uses only global in-plane components `[Ux, Uy]` for both rods. "
        "Analytic EB and Timoshenko vectors are reconstructed on the same external-to-joint slice coordinates, "
        "using the determinant-side right-arm convention and global transform `u*t + w*n` with "
        "`t_1=(1,0)`, `n_1=(0,1)`, `t_2=(cos(beta), sin(beta))`, and `n_2=(-sin(beta), cos(beta))`."
    )
    lines.append("")
    lines.append("MAC formula:")
    lines.append("")
    lines.append("```text")
    lines.append("MAC(phi, psi) = |phi dot psi|^2 / ((phi dot phi) * (psi dot psi))")
    lines.append("```")
    lines.append("")
    lines.append(
        "This centerline MAC ignores full cross-section deformation, torsional warping, and rotations; "
        "it is a diagnostic shape-identity filter, not a final proof."
    )
    lines.append("")
    lines.extend(mac_match_summary_lines(mac_match_rows, mac_warnings))
    lines.append("")
    lines.append(markdown_mac_summary_table(mac_match_rows))
    lines.append("")
    if mac_warnings:
        lines.append("MAC reconstruction warnings:")
        for warning in mac_warnings[:20]:
            lines.append(f"- {warning}")
        if len(mac_warnings) > 20:
            lines.append(f"- ... {len(mac_warnings) - 20} more warnings")
        lines.append("")
    lines.append("## High-Level Comparison With Previous Fused-Joint Diagnostic")
    lines.append("")
    lines.extend(fused_comparison_lines())
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
    lines.append(f"- geometry audit CSV: `{rel_path(OUTPUT_GEOMETRY_AUDIT_CSV)}`")
    lines.append(f"- mode metrics CSV: `{rel_path(OUTPUT_MODE_METRICS_CSV)}`")
    lines.append(f"- in-plane comparison CSV: `{rel_path(OUTPUT_IN_PLANE_CSV)}`")
    lines.append(f"- EB centerline MAC matrix CSV: `{rel_path(OUTPUT_MAC_MATRIX_EB_CSV)}`")
    lines.append(f"- Timoshenko centerline MAC matrix CSV: `{rel_path(OUTPUT_MAC_MATRIX_TIMO_CSV)}`")
    lines.append(f"- MAC match summary CSV: `{rel_path(OUTPUT_MAC_MATCHES_CSV)}`")
    lines.append(f"- report: `{rel_path(OUTPUT_REPORT)}`")
    lines.append(f"- mesh/input/output directory: `{rel_path(SOLID_OUTPUT_DIR)}`")
    lines.append("")
    lines.append("## Next Step")
    lines.append("")
    lines.append("Inspect representative in-plane and out-of-plane mode shapes visually, and consider a fuller shape metric including rotations/cross-section deformation before making stronger validation claims.")
    lines.append("")
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    audit = audit_tools()
    constraint_audit = audit_constraint_capability(audit)
    analytic_rows, analytic_warnings = compute_analytic_rows()
    (
        fem_omegas_by_case,
        mode_metrics,
        parse_results,
        solid_centerline_shapes,
        mesh_summaries,
        ccx_inputs,
        calculix_results,
        messages,
        geometry_label,
    ) = process_cases(audit, constraint_audit)
    sorted_rows = build_sorted_comparison_rows(analytic_rows, fem_omegas_by_case)
    in_plane_rows = build_in_plane_comparison_rows(analytic_rows, mode_metrics)
    mac_eb_rows, mac_timo_rows, mac_warnings = build_mac_matrix_rows(analytic_rows, solid_centerline_shapes)
    mac_match_rows = build_mac_match_rows(mac_eb_rows, mac_timo_rows)
    write_sorted_csv(sorted_rows)
    write_geometry_audit_csv(mesh_summaries)
    write_mode_metrics_csv(mode_metrics)
    write_in_plane_csv(in_plane_rows)
    write_mac_matrix_csv(OUTPUT_MAC_MATRIX_EB_CSV, mac_eb_rows)
    write_mac_matrix_csv(OUTPUT_MAC_MATRIX_TIMO_CSV, mac_timo_rows)
    write_mac_matches_csv(mac_match_rows)
    write_report(
        audit=audit,
        constraint_audit=constraint_audit,
        geometry_label=geometry_label,
        analytic_warnings=analytic_warnings,
        mesh_summaries=mesh_summaries,
        ccx_inputs=ccx_inputs,
        calculix_results=calculix_results,
        parse_results=parse_results,
        mode_metrics=mode_metrics,
        sorted_rows=sorted_rows,
        in_plane_rows=in_plane_rows,
        mac_match_rows=mac_match_rows,
        mac_warnings=mac_warnings,
        messages=messages,
    )
    print(f"wrote sorted comparison CSV: {OUTPUT_SORTED_CSV}")
    print(f"wrote geometry audit CSV: {OUTPUT_GEOMETRY_AUDIT_CSV}")
    print(f"wrote mode metrics CSV: {OUTPUT_MODE_METRICS_CSV}")
    print(f"wrote in-plane comparison CSV: {OUTPUT_IN_PLANE_CSV}")
    print(f"wrote EB MAC matrix CSV: {OUTPUT_MAC_MATRIX_EB_CSV}")
    print(f"wrote Timoshenko MAC matrix CSV: {OUTPUT_MAC_MATRIX_TIMO_CSV}")
    print(f"wrote MAC matches CSV: {OUTPUT_MAC_MATCHES_CSV}")
    print(f"wrote report: {OUTPUT_REPORT}")
    print(f"wrote mesh/input/output files under: {SOLID_OUTPUT_DIR}")
    if not any(result.success for result in calculix_results):
        print("no 3D solid modal results computed; see report for missing tool or solver failure details")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
