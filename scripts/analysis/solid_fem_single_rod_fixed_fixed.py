from __future__ import annotations

import csv
import os
import importlib.util
import math
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Sequence


# Diagnostic-only 3D solid FEM workflow seed.
#
# This script deliberately does not make Gmsh, CalculiX, or Code_Aster required
# project dependencies. If they are missing, it still writes reproducible input
# templates and a report explaining the next manual step.

L = 1.0
E = 1.0
RHO = 1.0
NU = 0.3
EPSILON_VALUES = [0.025, 0.05, 0.075]
REFERENCE_1D_MODES = 6
SOLID_MODES_REQUESTED = 16
GMSH_ELEMENT_ORDER = 2
CCX_TIMEOUT_SECONDS = 240
MODE_SHAPE_SLICE_COUNT = 41
BENDING_DOUBLET_REL_TOL = 0.08
GMSH_EXE_ENV = "GMSH_EXE"
GMSH_EXE_DEFAULT = r"D:\PHD\gmsh\gmsh-4.15.2-Windows64\gmsh.exe"
CCX_EXE_ENV = "CCX_EXE"
CCX_EXE_DEFAULT = r"D:\PHD\calculix\calculix_2.22_4win\ccx_static.exe"
CCX_SEARCH_ROOTS = [r"D:\PHD\calculix\calculix_2.22_4win"]
END_FACE_NODE_TOL = 1.0e-8

REPO_ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = REPO_ROOT / "scripts" / "analysis"
OUTPUT_DIR = REPO_ROOT / "results"
SOLID_OUTPUT_DIR = OUTPUT_DIR / "solid_fem_single_rod"
OUTPUT_CSV = OUTPUT_DIR / "single_rod_fixed_fixed_3d_solid_fem_comparison.csv"
OUTPUT_MODE_METRICS_CSV = OUTPUT_DIR / "single_rod_fixed_fixed_3d_solid_fem_mode_metrics.csv"
OUTPUT_BENDING_DOUBLET_CSV = OUTPUT_DIR / "single_rod_fixed_fixed_3d_solid_fem_bending_doublet_comparison.csv"
OUTPUT_REPORT = OUTPUT_DIR / "single_rod_fixed_fixed_3d_solid_fem_report.md"


@dataclass(frozen=True)
class ToolAudit:
    gmsh_exe: str | None
    gmsh_source: str | None
    gmsh_version: str | None
    gmsh_python_api: bool
    meshio_python_api: bool
    calculix_ccx: str | None
    calculix_source: str | None
    calculix_version: str | None
    code_aster_as_run: str | None
    code_aster_run_aster: str | None
    salome: str | None
    salome_meca: str | None

    @property
    def has_mesh_generator(self) -> bool:
        return bool(self.gmsh_python_api or self.gmsh_exe)

    @property
    def has_modal_solver(self) -> bool:
        return bool(self.calculix_ccx or self.code_aster_run_aster or self.code_aster_as_run)


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
class MeshSummary:
    epsilon: float
    file_path: Path
    file_kind: str
    node_count: int
    solid_element_count: int
    left_face_element_count: int | None
    right_face_element_count: int | None
    left_face_node_count: int | None
    right_face_node_count: int | None
    bbox_min: tuple[float, float, float] | None
    bbox_max: tuple[float, float, float] | None
    inferred_length: float | None
    inferred_radius: float | None
    expected_length: float
    expected_radius: float
    dimensions_match_expected: bool
    end_face_groups_present: bool
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class InpMeshData:
    nodes: dict[int, tuple[float, float, float]]
    solid_elements: dict[int, list[int]]
    bbox_min: tuple[float, float, float] | None
    bbox_max: tuple[float, float, float] | None
    radial_max: float


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
class CcxInputSummary:
    epsilon: float
    modal_input_path: Path
    mesh_include_path: Path
    solid_element_count: int
    left_fixed_node_count: int
    right_fixed_node_count: int
    node_sets_present: bool


@dataclass(frozen=True)
class ModeShapeParseResult:
    epsilon: float
    frd_path: Path
    success: bool
    message: str
    mode_shapes: dict[int, dict[int, tuple[float, float, float]]]
    source: str


@dataclass(frozen=True)
class ModeMetric:
    epsilon: float
    solid_mode: int
    omega: float
    lambda_fem: float
    parsed_node_count: int
    axial_energy_fraction: float
    transverse_y_energy_fraction: float
    transverse_z_energy_fraction: float
    transverse_energy_fraction: float
    mean_axial_energy_fraction: float
    mean_transverse_energy_fraction: float
    torsion_indicator: float
    bending_indicator: float
    warping_indicator: float
    dominant_transverse_axis: str
    classification: str
    classification_note: str


@dataclass(frozen=True)
class BendingDoubletComparison:
    epsilon: float
    bending_doublet_index: int
    mode_a: int | None
    mode_b: int | None
    class_a: str
    class_b: str
    lambda_a: float | None
    lambda_b: float | None
    lambda_fem_mean: float | None
    lambda_split_relative: float | None
    lambda_eb: float
    lambda_timoshenko: float
    fem_minus_eb: float | None
    fem_minus_timoshenko: float | None
    status: str
    notes: str


def safe_float_token(value: float) -> str:
    return f"{float(value):.6g}".replace("-", "m").replace(".", "p")


def python_module_available(name: str) -> bool:
    return importlib.util.find_spec(name) is not None


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
    if not text:
        return None
    return text.splitlines()[0].strip()


def audit_tools() -> ToolAudit:
    gmsh_exe, gmsh_source = resolve_gmsh_executable()
    ccx_exe, ccx_source = resolve_ccx_executable()
    return ToolAudit(
        gmsh_exe=gmsh_exe,
        gmsh_source=gmsh_source,
        gmsh_version=executable_version(gmsh_exe, "--version"),
        gmsh_python_api=python_module_available("gmsh"),
        meshio_python_api=python_module_available("meshio"),
        calculix_ccx=ccx_exe,
        calculix_source=ccx_source,
        calculix_version=executable_version(ccx_exe, "-v"),
        code_aster_as_run=shutil.which("as_run"),
        code_aster_run_aster=shutil.which("run_aster"),
        salome=shutil.which("salome"),
        salome_meca=shutil.which("salome_meca"),
    )


def load_single_rod_reference_module() -> ModuleType:
    module_path = ANALYSIS_DIR / "compare_single_rod_eb_timoshenko.py"
    spec = importlib.util.spec_from_file_location("single_rod_eb_timoshenko_reference", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load reference module from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def case_paths(epsilon: float) -> CasePaths:
    token = safe_float_token(float(epsilon))
    case_dir = SOLID_OUTPUT_DIR / f"eps_{token}"
    stem = f"single_rod_fixed_fixed_eps{token}"
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


def mesh_size_for_epsilon(epsilon: float) -> float:
    radius = float(epsilon) * L
    return min(L / 80.0, radius / 2.0)


def write_gmsh_geo(epsilon: float, paths: CasePaths) -> None:
    radius = float(epsilon) * L
    mesh_size = mesh_size_for_epsilon(float(epsilon))
    paths.case_dir.mkdir(parents=True, exist_ok=True)
    paths.geo.write_text(
        f"""// Diagnostic-only fixed-fixed circular rod solid mesh.
// Generate manually, for example:
// gmsh "{paths.geo.name}" -3 -format inp -o "{paths.gmsh_inp.name}"
// gmsh "{paths.geo.name}" -3 -format msh4 -o "{paths.msh.name}"

SetFactory("OpenCASCADE");

L = {L:.17g};
R = {radius:.17g};
h = {mesh_size:.17g};
tol = 1e-8;

Cylinder(1) = {{0, 0, 0, L, 0, 0, R, 2*Pi}};

vols[] = Volume In BoundingBox {{-tol, -R-tol, -R-tol, L+tol, R+tol, R+tol}};
left_faces[] = Surface In BoundingBox {{-tol, -R-tol, -R-tol, tol, R+tol, R+tol}};
right_faces[] = Surface In BoundingBox {{L-tol, -R-tol, -R-tol, L+tol, R+tol, R+tol}};

Physical Volume("SOLID", 1) = vols[];
Physical Surface("FIXED_LEFT", 2) = left_faces[];
Physical Surface("FIXED_RIGHT", 3) = right_faces[];

Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Mesh.ElementOrder = {GMSH_ELEMENT_ORDER};
Mesh.MshFileVersion = 4.1;
""",
        encoding="utf-8",
    )


def read_gmsh_inp_mesh_data(path: Path) -> InpMeshData:
    nodes: dict[int, tuple[float, float, float]] = {}
    solid_elements: dict[int, list[int]] = {}
    bbox_min: tuple[float, float, float] | None = None
    bbox_max: tuple[float, float, float] | None = None
    radial_max = 0.0
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
            radial_max = max(radial_max, math.hypot(point[1], point[2]))
        elif mode == "ELEMENT" and current_element_type.startswith("C3D") and len(parts) >= 2:
            try:
                element_id = int(parts[0])
                solid_elements[element_id] = [int(part) for part in parts[1:]]
            except ValueError:
                continue

    return InpMeshData(
        nodes=nodes,
        solid_elements=solid_elements,
        bbox_min=bbox_min,
        bbox_max=bbox_max,
        radial_max=radial_max,
    )


def coordinate_end_node_sets(mesh: InpMeshData) -> tuple[list[int], list[int]]:
    tol = max(END_FACE_NODE_TOL, 1.0e-7 * L)
    left = sorted(node_id for node_id, point in mesh.nodes.items() if abs(point[0]) <= tol)
    right = sorted(node_id for node_id, point in mesh.nodes.items() if abs(point[0] - L) <= tol)
    return left, right


def format_id_lines(ids: Sequence[int], *, per_line: int = 16) -> list[str]:
    lines: list[str] = []
    for start in range(0, len(ids), per_line):
        lines.append(", ".join(str(item) for item in ids[start : start + per_line]))
    return lines


def ccx_float(value: float) -> str:
    return f"{float(value):.12E}"


def write_calculix_mesh_include(paths: CasePaths, mesh: InpMeshData) -> None:
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


def write_calculix_template(paths: CasePaths, mesh: InpMeshData) -> tuple[int, int]:
    left_nodes, right_nodes = coordinate_end_node_sets(mesh)
    write_calculix_mesh_include(paths, mesh)
    paths.ccx_modal_inp.write_text(
        "\n".join(
            [
                "** Diagnostic-only CalculiX modal input.",
                "** End node sets are generated from coordinates, not trusted to Gmsh Abaqus export.",
                f"*INCLUDE, INPUT={paths.ccx_mesh_inp.name}",
                "*NSET, NSET=LEFT_FIXED",
                *format_id_lines(left_nodes),
                "*NSET, NSET=RIGHT_FIXED",
                *format_id_lines(right_nodes),
                "*MATERIAL, NAME=MAT",
                "*ELASTIC",
                f"{E:.17g}, {NU:.17g}",
                "*DENSITY",
                f"{RHO:.17g}",
                "*SOLID SECTION, ELSET=SOLID, MATERIAL=MAT",
                "*BOUNDARY",
                "LEFT_FIXED, 1, 3, 0",
                "RIGHT_FIXED, 1, 3, 0",
                "*STEP",
                "*FREQUENCY",
                f"{SOLID_MODES_REQUESTED}",
                "*NODE FILE",
                "U",
                "*END STEP",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return len(left_nodes), len(right_nodes)


def write_calculix_inputs(paths: CasePaths, epsilon: float) -> CcxInputSummary:
    mesh = read_gmsh_inp_mesh_data(paths.gmsh_inp)
    left_count, right_count = write_calculix_template(paths, mesh)
    return CcxInputSummary(
        epsilon=float(epsilon),
        modal_input_path=paths.ccx_modal_inp,
        mesh_include_path=paths.ccx_mesh_inp,
        solid_element_count=len(mesh.solid_elements),
        left_fixed_node_count=left_count,
        right_fixed_node_count=right_count,
        node_sets_present=left_count > 0 and right_count > 0,
    )


def interesting_gmsh_lines(text: str, *, prefix: str) -> list[str]:
    out: list[str] = []
    for line in text.splitlines():
        lower = line.lower()
        if "warning" in lower or "error" in lower:
            out.append(f"{prefix}: {line.strip()}")
    return out


def generate_mesh_with_gmsh_python(paths: CasePaths) -> tuple[bool, str, tuple[str, ...]]:
    try:
        import gmsh  # type: ignore[import-not-found]
    except ImportError as exc:
        return False, f"gmsh Python API unavailable: {exc}", (str(exc),)

    try:
        gmsh.initialize()
        gmsh.open(str(paths.geo))
        gmsh.model.mesh.generate(3)
        gmsh.write(str(paths.msh))
        gmsh.write(str(paths.gmsh_inp))
        return True, "generated with gmsh Python API", ()
    except Exception as exc:  # pragma: no cover - depends on external gmsh
        return False, f"gmsh Python API mesh generation failed: {exc}", (str(exc),)
    finally:
        try:
            gmsh.finalize()
        except Exception:
            pass


def generate_mesh_with_gmsh_cli(paths: CasePaths, gmsh_exe: str) -> tuple[bool, str, tuple[str, ...]]:
    commands = [
        ("msh4", [gmsh_exe, str(paths.geo), "-3", "-format", "msh4", "-o", str(paths.msh)]),
        ("inp", [gmsh_exe, str(paths.geo), "-3", "-format", "inp", "-o", str(paths.gmsh_inp)]),
    ]
    gmsh_messages: list[str] = []
    try:
        for label, command in commands:
            completed = subprocess.run(
                command,
                cwd=str(paths.case_dir),
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=CCX_TIMEOUT_SECONDS,
            )
            gmsh_messages.extend(interesting_gmsh_lines(completed.stdout, prefix=label))
            gmsh_messages.extend(interesting_gmsh_lines(completed.stderr, prefix=label))
    except subprocess.CalledProcessError as exc:  # pragma: no cover - depends on external gmsh
        gmsh_messages.extend(interesting_gmsh_lines(exc.stdout or "", prefix="failed"))
        gmsh_messages.extend(interesting_gmsh_lines(exc.stderr or "", prefix="failed"))
        gmsh_messages.append(str(exc))
        return False, f"gmsh CLI mesh generation failed: {exc}", tuple(gmsh_messages)
    except Exception as exc:  # pragma: no cover - depends on external gmsh
        gmsh_messages.append(str(exc))
        return False, f"gmsh CLI mesh generation failed: {exc}", tuple(gmsh_messages)
    return True, "mesh generated successfully with Gmsh executable", tuple(gmsh_messages)


def run_calculix(paths: CasePaths, ccx: str, epsilon: float) -> CalculixResult:
    if not paths.ccx_modal_inp.exists():
        return CalculixResult(
            epsilon=float(epsilon),
            success=False,
            message="CalculiX not run because the modal input was not generated.",
            stdout_path=paths.ccx_stdout,
            stderr_path=paths.ccx_stderr,
            dat_path=paths.ccx_dat,
            frd_path=paths.ccx_frd,
            parsed_omegas=(),
            parse_source="not_run",
        )
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
    except Exception as exc:  # pragma: no cover - depends on external ccx
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
        return CalculixResult(
            epsilon=float(epsilon),
            success=False,
            message=(
                f"CalculiX returned exit code {completed.returncode}. "
                f"stdout={completed.stdout[-500:]!r}; stderr={completed.stderr[-500:]!r}"
            ),
            stdout_path=paths.ccx_stdout,
            stderr_path=paths.ccx_stderr,
            dat_path=paths.ccx_dat,
            frd_path=paths.ccx_frd,
            parsed_omegas=(),
            parse_source="ccx_nonzero_exit",
        )
    omegas, parse_source = parse_calculix_eigen_omegas(paths.ccx_dat, paths.ccx_stdout)
    return CalculixResult(
        epsilon=float(epsilon),
        success=True,
        message=f"CalculiX completed; parsed {len(omegas)} frequencies from {parse_source}.",
        stdout_path=paths.ccx_stdout,
        stderr_path=paths.ccx_stderr,
        dat_path=paths.ccx_dat,
        frd_path=paths.ccx_frd,
        parsed_omegas=tuple(omegas),
        parse_source=parse_source,
    )


def parse_calculix_eigen_omegas(dat_path: Path, stdout_path: Path) -> tuple[list[float], str]:
    if not dat_path.exists() and not stdout_path.exists():
        return [], "missing .dat/stdout"
    omegas: list[float] = []
    eigenvalue_pattern = re.compile(r"eigenvalue\s*[:=]?\s*([-+0-9.DEded]+)", re.IGNORECASE)
    frequency_pattern = re.compile(r"frequency\s*[:=]?\s*([-+0-9.DEded]+)", re.IGNORECASE)

    dat_text = dat_path.read_text(encoding="utf-8", errors="ignore") if dat_path.exists() else ""
    for match in eigenvalue_pattern.finditer(dat_text):
        try:
            eigenvalue = float(match.group(1).replace("D", "E").replace("d", "E"))
        except ValueError:
            continue
        if eigenvalue > 0.0:
            omegas.append(math.sqrt(eigenvalue))
    if omegas:
        return omegas[:SOLID_MODES_REQUESTED], "CalculiX .dat eigenvalue entries interpreted as Omega^2"

    cyclic_frequencies: list[float] = []
    for match in frequency_pattern.finditer(dat_text):
        try:
            value = float(match.group(1).replace("D", "E").replace("d", "E"))
        except ValueError:
            continue
        if value > 0.0:
            cyclic_frequencies.append(value)
    if cyclic_frequencies:
        return [2.0 * math.pi * value for value in cyclic_frequencies[:SOLID_MODES_REQUESTED]], (
            "CalculiX .dat frequency entries interpreted as cyclic frequency f; Omega=2*pi*f"
        )

    for line in dat_text.splitlines():
        parts = line.split()
        if len(parts) < 2 or not parts[0].isdigit():
            continue
        try:
            eigenvalue = float(parts[1].replace("D", "E"))
        except ValueError:
            continue
        if eigenvalue > 0.0:
            omegas.append(math.sqrt(eigenvalue))
        if len(omegas) >= SOLID_MODES_REQUESTED:
            break
    if omegas:
        return omegas, "CalculiX .dat numeric table second column interpreted as Omega^2"

    stdout_text = stdout_path.read_text(encoding="utf-8", errors="ignore") if stdout_path.exists() else ""
    for match in frequency_pattern.finditer(stdout_text):
        try:
            value = float(match.group(1).replace("D", "E").replace("d", "E"))
        except ValueError:
            continue
        if value > 0.0:
            cyclic_frequencies.append(value)
    if cyclic_frequencies:
        return [2.0 * math.pi * value for value in cyclic_frequencies[:SOLID_MODES_REQUESTED]], (
            "CalculiX stdout frequency entries interpreted as cyclic frequency f; Omega=2*pi*f"
        )
    return [], "no recognizable CalculiX eigenfrequency entries"


FRD_NUMBER_PATTERN = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?")


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
            frd_path=paths.ccx_frd,
            success=False,
            message="CalculiX .frd file is missing; no mode-shape vectors parsed.",
            mode_shapes={},
            source="missing_frd",
        )

    mode_shapes: dict[int, dict[int, tuple[float, float, float]]] = {}
    current_mode: int | None = None
    collecting_displacements = False

    for raw_line in paths.ccx_frd.read_text(encoding="utf-8", errors="ignore").splitlines():
        stripped = raw_line.strip()
        if not stripped:
            continue
        if stripped.startswith("1PMODE"):
            parts = stripped.split()
            try:
                current_mode = int(parts[-1])
            except (ValueError, IndexError):
                current_mode = None
            collecting_displacements = False
            continue
        if stripped.startswith("-4") and "DISP" in stripped.upper():
            if current_mode is not None:
                mode_shapes.setdefault(current_mode, {})
                collecting_displacements = True
            continue
        if collecting_displacements and stripped.startswith("-1"):
            values = numeric_tokens(stripped[2:])
            if len(values) < 4 or current_mode is None:
                continue
            node_id = int(values[0])
            mode_shapes[current_mode][node_id] = (float(values[1]), float(values[2]), float(values[3]))
            continue
        if collecting_displacements and (
            stripped.startswith("-3")
            or stripped.startswith("-4")
            or stripped.startswith("100C")
            or stripped.startswith("1P")
        ):
            collecting_displacements = False

    parsed_modes = sorted(mode for mode, shape in mode_shapes.items() if shape)
    if not parsed_modes:
        return ModeShapeParseResult(
            epsilon=float(epsilon),
            frd_path=paths.ccx_frd,
            success=False,
            message="No displacement eigenvector blocks were recognized in the CalculiX .frd file.",
            mode_shapes={},
            source="frd_no_displacement_blocks",
        )
    return ModeShapeParseResult(
        epsilon=float(epsilon),
        frd_path=paths.ccx_frd,
        success=True,
        message=f"Parsed displacement eigenvectors for {len(parsed_modes)} modes from CalculiX .frd DISP blocks.",
        mode_shapes={mode: mode_shapes[mode] for mode in parsed_modes},
        source="CalculiX .frd DISP blocks from *NODE FILE U",
    )


def classify_mode_shape(
    *,
    axial_fraction: float,
    transverse_y_fraction: float,
    transverse_z_fraction: float,
    mean_axial_fraction: float,
    mean_transverse_fraction: float,
    torsion_indicator: float,
    warping_indicator: float,
) -> tuple[str, str]:
    dominant_axis = "y" if transverse_y_fraction >= transverse_z_fraction else "z"
    if axial_fraction >= 0.65 and mean_axial_fraction >= 0.25:
        return "axial_like", "dominant axial displacement energy and axial slice-mean motion"
    if (
        torsion_indicator >= 0.65
        and mean_transverse_fraction <= 0.25
        and axial_fraction <= 0.35
    ):
        return "torsion_like", "dominant circumferential transverse pattern with weak centroid translation"
    if mean_transverse_fraction >= 0.20 and (transverse_y_fraction + transverse_z_fraction) >= 0.45:
        if dominant_axis == "y":
            return "bending_y_like", "dominant centroidal transverse motion in y-like direction"
        return "bending_z_like", "dominant centroidal transverse motion in z-like direction"
    if warping_indicator >= 0.75:
        return "mixed_or_unclassified", "large within-slice residual/warping relative to total displacement energy"
    return "mixed_or_unclassified", "no single diagnostic indicator dominates"


def compute_mode_metric(
    *,
    epsilon: float,
    solid_mode: int,
    omega: float,
    mesh: InpMeshData,
    mode_shape: dict[int, tuple[float, float, float]],
) -> ModeMetric | None:
    total_energy = 0.0
    axial_energy = 0.0
    transverse_y_energy = 0.0
    transverse_z_energy = 0.0
    parsed_node_count = 0
    slice_data: dict[int, list[float]] = defaultdict(lambda: [0.0] * 9)

    for node_id, point in mesh.nodes.items():
        displacement = mode_shape.get(node_id)
        if displacement is None:
            continue
        x, y, z = point
        ux, uy, uz = displacement
        local_total = ux * ux + uy * uy + uz * uz
        total_energy += local_total
        axial_energy += ux * ux
        transverse_y_energy += uy * uy
        transverse_z_energy += uz * uz
        parsed_node_count += 1

        if L > 0.0:
            slice_index = int(math.floor(max(0.0, min(0.999999999999, x / L)) * MODE_SHAPE_SLICE_COUNT))
        else:
            slice_index = 0
        bucket = slice_data[slice_index]
        bucket[0] += 1.0
        bucket[1] += ux
        bucket[2] += uy
        bucket[3] += uz
        bucket[4] += local_total
        bucket[5] += -z * uy + y * uz
        bucket[6] += y * y + z * z
        bucket[7] += uy * uy + uz * uz
        bucket[8] += ux * ux

    if total_energy <= 0.0 or parsed_node_count == 0:
        return None

    mean_total_energy = 0.0
    mean_axial_energy = 0.0
    mean_transverse_energy = 0.0
    torsion_explained_energy = 0.0
    transverse_energy = transverse_y_energy + transverse_z_energy
    for bucket in slice_data.values():
        count, sum_ux, sum_uy, sum_uz, _sum_total, torsion_num, torsion_den, _sum_trans, _sum_axial = bucket
        if count <= 0.0:
            continue
        mean_axial_energy += (sum_ux * sum_ux) / count
        mean_transverse_energy += (sum_uy * sum_uy + sum_uz * sum_uz) / count
        mean_total_energy += (sum_ux * sum_ux + sum_uy * sum_uy + sum_uz * sum_uz) / count
        if torsion_den > 0.0:
            torsion_explained_energy += (torsion_num * torsion_num) / torsion_den

    axial_fraction = axial_energy / total_energy
    transverse_y_fraction = transverse_y_energy / total_energy
    transverse_z_fraction = transverse_z_energy / total_energy
    transverse_fraction = transverse_energy / total_energy
    mean_axial_fraction = mean_axial_energy / total_energy
    mean_transverse_fraction = mean_transverse_energy / total_energy
    torsion_indicator = torsion_explained_energy / transverse_energy if transverse_energy > 0.0 else 0.0
    bending_indicator = mean_transverse_fraction
    warping_indicator = max(0.0, min(1.0, (total_energy - mean_total_energy) / total_energy))
    classification, note = classify_mode_shape(
        axial_fraction=axial_fraction,
        transverse_y_fraction=transverse_y_fraction,
        transverse_z_fraction=transverse_z_fraction,
        mean_axial_fraction=mean_axial_fraction,
        mean_transverse_fraction=mean_transverse_fraction,
        torsion_indicator=torsion_indicator,
        warping_indicator=warping_indicator,
    )
    dominant_axis = "y" if transverse_y_fraction >= transverse_z_fraction else "z"
    return ModeMetric(
        epsilon=float(epsilon),
        solid_mode=int(solid_mode),
        omega=float(omega),
        lambda_fem=lambda_from_omega(float(omega), float(epsilon)),
        parsed_node_count=parsed_node_count,
        axial_energy_fraction=axial_fraction,
        transverse_y_energy_fraction=transverse_y_fraction,
        transverse_z_energy_fraction=transverse_z_fraction,
        transverse_energy_fraction=transverse_fraction,
        mean_axial_energy_fraction=mean_axial_fraction,
        mean_transverse_energy_fraction=mean_transverse_fraction,
        torsion_indicator=torsion_indicator,
        bending_indicator=bending_indicator,
        warping_indicator=warping_indicator,
        dominant_transverse_axis=dominant_axis,
        classification=classification,
        classification_note=note,
    )


def compute_mode_metrics_for_case(
    *,
    epsilon: float,
    paths: CasePaths,
    omegas: Sequence[float],
) -> tuple[list[ModeMetric], ModeShapeParseResult]:
    parse_result = parse_calculix_frd_mode_shapes(paths, float(epsilon))
    if not parse_result.success:
        return [], parse_result
    mesh = read_gmsh_inp_mesh_data(paths.gmsh_inp)
    metrics: list[ModeMetric] = []
    for index, omega in enumerate(omegas, start=1):
        mode_shape = parse_result.mode_shapes.get(index)
        if mode_shape is None:
            continue
        metric = compute_mode_metric(
            epsilon=float(epsilon),
            solid_mode=index,
            omega=float(omega),
            mesh=mesh,
            mode_shape=mode_shape,
        )
        if metric is not None:
            metrics.append(metric)
    return metrics, parse_result


def mode_shape_parse_lines(results: Sequence[ModeShapeParseResult]) -> list[str]:
    if not results:
        return ["Mode-shape parsing was not attempted because no successful CalculiX modal result was available."]
    lines: list[str] = []
    for result in results:
        status = "succeeded" if result.success else "failed"
        lines.append(f"- epsilon={result.epsilon:g}: {status}; {result.message}")
        lines.append(f"  frd: `{result.frd_path.relative_to(REPO_ROOT)}`")
        lines.append(f"  source: {result.source}")
    return lines


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


def mesh_dimensions_match(
    inferred_length: float | None,
    inferred_radius: float | None,
    expected_radius: float,
) -> bool:
    if inferred_length is None or inferred_radius is None:
        return False
    length_tol = max(1.0e-7, 1.0e-6 * L)
    radius_tol = max(1.0e-7, 1.0e-5 * expected_radius)
    return abs(inferred_length - L) <= length_tol and abs(inferred_radius - expected_radius) <= radius_tol


def infer_length_and_radius(
    bbox_min: tuple[float, float, float] | None,
    bbox_max: tuple[float, float, float] | None,
    radial_max: float,
) -> tuple[float | None, float | None]:
    if bbox_min is None or bbox_max is None:
        return None, None
    return bbox_max[0] - bbox_min[0], radial_max


def parse_keyword_args(line: str) -> dict[str, str]:
    args: dict[str, str] = {}
    for part in line.split(",")[1:]:
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        args[key.strip().upper()] = value.strip()
    return args


def parse_inp_mesh_summary(path: Path, epsilon: float, gmsh_messages: Sequence[str]) -> MeshSummary:
    expected_radius = float(epsilon) * L
    node_count = 0
    solid_element_count = 0
    nodes: dict[int, tuple[float, float, float]] = {}
    element_nodes: dict[int, list[int]] = {}
    target_element_ids: dict[str, set[int]] = {"FIXED_LEFT": set(), "FIXED_RIGHT": set()}
    target_nodes: dict[str, set[int]] = {"FIXED_LEFT": set(), "FIXED_RIGHT": set()}
    bbox_min: tuple[float, float, float] | None = None
    bbox_max: tuple[float, float, float] | None = None
    radial_max = 0.0
    mode = ""
    current_element_type = ""
    current_elset = ""
    current_set_name = ""
    warnings = list(gmsh_messages)

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
                current_elset = args.get("ELSET", "").upper()
            elif upper.startswith("*NSET"):
                mode = "NSET"
                current_set_name = args.get("NSET", "").upper()
            elif upper.startswith("*ELSET"):
                mode = "ELSET"
                current_set_name = args.get("ELSET", "").upper()
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
            node_count += 1
            bbox_min, bbox_max = update_bbox(bbox_min, bbox_max, point)
            radial_max = max(radial_max, math.hypot(point[1], point[2]))
        elif mode == "ELEMENT" and len(parts) >= 2:
            try:
                element_id = int(parts[0])
                connectivity = [int(part) for part in parts[1:]]
            except ValueError:
                continue
            element_nodes[element_id] = connectivity
            if current_element_type.startswith("C3D"):
                solid_element_count += 1
            if current_elset in target_element_ids:
                target_element_ids[current_elset].add(element_id)
                target_nodes[current_elset].update(connectivity)
        elif mode == "NSET" and current_set_name in target_nodes:
            for part in parts:
                try:
                    target_nodes[current_set_name].add(int(part))
                except ValueError:
                    continue
        elif mode == "ELSET" and current_set_name in target_element_ids:
            for part in parts:
                try:
                    target_element_ids[current_set_name].add(int(part))
                except ValueError:
                    continue

    for group, ids in target_element_ids.items():
        for element_id in ids:
            target_nodes[group].update(element_nodes.get(element_id, []))

    left_elements = len(target_element_ids["FIXED_LEFT"]) or None
    right_elements = len(target_element_ids["FIXED_RIGHT"]) or None
    left_nodes = len(target_nodes["FIXED_LEFT"]) or None
    right_nodes = len(target_nodes["FIXED_RIGHT"]) or None
    end_groups_present = bool(left_nodes and right_nodes)
    if not end_groups_present:
        warnings.append("end-face physical groups not detected; CalculiX boundary conditions may fail")

    inferred_length, inferred_radius = infer_length_and_radius(bbox_min, bbox_max, radial_max)
    return MeshSummary(
        epsilon=float(epsilon),
        file_path=path,
        file_kind="inp",
        node_count=node_count,
        solid_element_count=solid_element_count,
        left_face_element_count=left_elements,
        right_face_element_count=right_elements,
        left_face_node_count=left_nodes,
        right_face_node_count=right_nodes,
        bbox_min=bbox_min,
        bbox_max=bbox_max,
        inferred_length=inferred_length,
        inferred_radius=inferred_radius,
        expected_length=L,
        expected_radius=expected_radius,
        dimensions_match_expected=mesh_dimensions_match(inferred_length, inferred_radius, expected_radius),
        end_face_groups_present=end_groups_present,
        warnings=tuple(warnings),
    )


def parse_msh_physical_names(lines: Sequence[str]) -> dict[tuple[int, int], str]:
    names: dict[tuple[int, int], str] = {}
    try:
        start = lines.index("$PhysicalNames")
        count = int(lines[start + 1].strip())
    except (ValueError, IndexError):
        return names
    for offset in range(count):
        parts = lines[start + 2 + offset].split(maxsplit=2)
        if len(parts) < 3:
            continue
        try:
            dim = int(parts[0])
            tag = int(parts[1])
        except ValueError:
            continue
        names[(dim, tag)] = parts[2].strip().strip('"')
    return names


def parse_msh_entity_physical_tags(lines: Sequence[str]) -> dict[tuple[int, int], list[int]]:
    entity_tags: dict[tuple[int, int], list[int]] = {}
    try:
        start = lines.index("$Entities")
        counts = [int(part) for part in lines[start + 1].split()]
    except (ValueError, IndexError):
        return entity_tags
    index = start + 2
    for dim, count in enumerate(counts):
        for _ in range(count):
            parts = lines[index].split()
            index += 1
            try:
                tag = int(parts[0])
                phys_index = 4 if dim == 0 else 7
                physical_count = int(parts[phys_index])
                physical_tags = [int(value) for value in parts[phys_index + 1 : phys_index + 1 + physical_count]]
            except (ValueError, IndexError):
                continue
            entity_tags[(dim, tag)] = physical_tags
    return entity_tags


def parse_msh_nodes(lines: Sequence[str]) -> tuple[int, dict[int, tuple[float, float, float]], tuple[float, float, float] | None, tuple[float, float, float] | None, float]:
    nodes: dict[int, tuple[float, float, float]] = {}
    bbox_min: tuple[float, float, float] | None = None
    bbox_max: tuple[float, float, float] | None = None
    radial_max = 0.0
    try:
        index = lines.index("$Nodes") + 1
        header = [int(part) for part in lines[index].split()]
    except (ValueError, IndexError):
        return 0, nodes, bbox_min, bbox_max, radial_max
    index += 1
    block_count = header[0]
    for _ in range(block_count):
        entity_dim, entity_tag, parametric, node_count = [int(part) for part in lines[index].split()]
        _ = (entity_dim, entity_tag)
        index += 1
        node_tags = [int(lines[index + offset].strip()) for offset in range(node_count)]
        index += node_count
        for node_tag in node_tags:
            parts = [float(part) for part in lines[index].split()]
            index += 1
            point = (parts[0], parts[1], parts[2])
            if parametric:
                point = (parts[0], parts[1], parts[2])
            nodes[node_tag] = point
            bbox_min, bbox_max = update_bbox(bbox_min, bbox_max, point)
            radial_max = max(radial_max, math.hypot(point[1], point[2]))
    return len(nodes), nodes, bbox_min, bbox_max, radial_max


def parse_msh_elements(
    lines: Sequence[str],
    physical_names: dict[tuple[int, int], str],
    entity_physical_tags: dict[tuple[int, int], list[int]],
) -> tuple[int, dict[str, int], dict[str, set[int]]]:
    solid_element_count = 0
    face_element_counts = {"FIXED_LEFT": 0, "FIXED_RIGHT": 0}
    face_nodes = {"FIXED_LEFT": set(), "FIXED_RIGHT": set()}
    try:
        index = lines.index("$Elements") + 1
        header = [int(part) for part in lines[index].split()]
    except (ValueError, IndexError):
        return solid_element_count, face_element_counts, face_nodes
    index += 1
    block_count = header[0]
    for _ in range(block_count):
        entity_dim, entity_tag, element_type, element_count = [int(part) for part in lines[index].split()]
        _ = element_type
        index += 1
        physical_group_names = {
            physical_names.get((entity_dim, physical_tag), "").upper()
            for physical_tag in entity_physical_tags.get((entity_dim, entity_tag), [])
        }
        for _ in range(element_count):
            parts = [int(part) for part in lines[index].split()]
            index += 1
            connectivity = parts[1:]
            if entity_dim == 3:
                solid_element_count += 1
            for group in ("FIXED_LEFT", "FIXED_RIGHT"):
                if group in physical_group_names:
                    face_element_counts[group] += 1
                    face_nodes[group].update(connectivity)
    return solid_element_count, face_element_counts, face_nodes


def parse_msh_mesh_summary(path: Path, epsilon: float, gmsh_messages: Sequence[str]) -> MeshSummary:
    expected_radius = float(epsilon) * L
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    physical_names = parse_msh_physical_names(lines)
    entity_physical_tags = parse_msh_entity_physical_tags(lines)
    node_count, _nodes, bbox_min, bbox_max, radial_max = parse_msh_nodes(lines)
    solid_element_count, face_element_counts, face_nodes = parse_msh_elements(
        lines,
        physical_names,
        entity_physical_tags,
    )
    left_elements = face_element_counts["FIXED_LEFT"] or None
    right_elements = face_element_counts["FIXED_RIGHT"] or None
    left_nodes = len(face_nodes["FIXED_LEFT"]) or None
    right_nodes = len(face_nodes["FIXED_RIGHT"]) or None
    end_groups_present = bool(left_nodes and right_nodes)
    warnings = list(gmsh_messages)
    if not end_groups_present:
        warnings.append("end-face physical groups not detected; CalculiX boundary conditions may fail")
    inferred_length, inferred_radius = infer_length_and_radius(bbox_min, bbox_max, radial_max)
    return MeshSummary(
        epsilon=float(epsilon),
        file_path=path,
        file_kind="msh",
        node_count=node_count,
        solid_element_count=solid_element_count,
        left_face_element_count=left_elements,
        right_face_element_count=right_elements,
        left_face_node_count=left_nodes,
        right_face_node_count=right_nodes,
        bbox_min=bbox_min,
        bbox_max=bbox_max,
        inferred_length=inferred_length,
        inferred_radius=inferred_radius,
        expected_length=L,
        expected_radius=expected_radius,
        dimensions_match_expected=mesh_dimensions_match(inferred_length, inferred_radius, expected_radius),
        end_face_groups_present=end_groups_present,
        warnings=tuple(warnings),
    )


def summarize_generated_meshes(
    epsilon: float,
    paths: CasePaths,
    gmsh_messages: Sequence[str],
) -> list[MeshSummary]:
    summaries: list[MeshSummary] = []
    if paths.msh.exists():
        summaries.append(parse_msh_mesh_summary(paths.msh, float(epsilon), gmsh_messages))
    if paths.gmsh_inp.exists():
        summaries.append(parse_inp_mesh_summary(paths.gmsh_inp, float(epsilon), gmsh_messages))
    return summaries


def lambda_from_omega(omega: float, epsilon: float) -> float:
    radius = float(epsilon) * L
    area = math.pi * radius**2
    inertia = math.pi * radius**4 / 4.0
    return math.sqrt(float(omega)) * (L / 2.0) * (RHO * area / (E * inertia)) ** 0.25


def compute_1d_references() -> list[dict[str, float | int]]:
    reference = load_single_rod_reference_module()
    alpha_roots = reference.fixed_fixed_eb_roots(REFERENCE_1D_MODES)
    rows: list[dict[str, float | int]] = []
    for epsilon in EPSILON_VALUES:
        section = reference.section_from_epsilon(float(epsilon))
        eb_omegas = [reference.omega_eb(alpha, section) for alpha in alpha_roots]
        timo_omegas = reference.find_timoshenko_roots(section, eb_omegas, REFERENCE_1D_MODES)
        for index in range(REFERENCE_1D_MODES):
            rows.append(
                {
                    "epsilon": float(epsilon),
                    "mode_1d": index + 1,
                    "Omega_EB": float(eb_omegas[index]),
                    "Lambda_EB": float(reference.lambda_from_omega(float(eb_omegas[index]), section)),
                    "Omega_Timoshenko": float(timo_omegas[index]),
                    "Lambda_Timoshenko": float(reference.lambda_from_omega(float(timo_omegas[index]), section)),
                }
            )
    return rows


def build_comparison_rows(
    reference_rows: Sequence[dict[str, float | int]],
    fem_omegas_by_epsilon: dict[float, list[float]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for reference in reference_rows:
        epsilon = float(reference["epsilon"])
        mode_1d = int(reference["mode_1d"])
        fem_omegas = fem_omegas_by_epsilon.get(epsilon, [])
        doublet_indices = [2 * mode_1d - 1, 2 * mode_1d]
        emitted = False
        for member, solid_index in enumerate(doublet_indices, start=1):
            if 1 <= solid_index <= len(fem_omegas):
                omega_fem = float(fem_omegas[solid_index - 1])
                rows.append(
                    {
                        **reference,
                        "solid_sorted_mode": solid_index,
                        "solid_doublet_member": member,
                        "Omega_FEM": omega_fem,
                        "Lambda_FEM": lambda_from_omega(omega_fem, epsilon),
                        "status": "computed_preliminary_sorted_doublet",
                        "notes": "Circular solid bending modes are expected to appear as near-doublets; identity not overclaimed.",
                    }
                )
                emitted = True
        if not emitted:
            rows.append(
                {
                    **reference,
                    "solid_sorted_mode": "",
                    "solid_doublet_member": "",
                    "Omega_FEM": "",
                    "Lambda_FEM": "",
                    "status": "not_computed",
                    "notes": "No 3D solid modal frequency was computed in this run.",
                }
            )
    return rows


def bending_like(metric: ModeMetric) -> bool:
    return metric.classification in {"bending_y_like", "bending_z_like"}


def reference_lookup(reference_rows: Sequence[dict[str, float | int]]) -> dict[tuple[float, int], dict[str, float | int]]:
    lookup: dict[tuple[float, int], dict[str, float | int]] = {}
    for row in reference_rows:
        lookup[(float(row["epsilon"]), int(row["mode_1d"]))] = row
    return lookup


def pair_bending_modes(metrics: Sequence[ModeMetric]) -> list[tuple[ModeMetric, ModeMetric, float]]:
    bending_modes = sorted((metric for metric in metrics if bending_like(metric)), key=lambda item: item.solid_mode)
    pairs: list[tuple[ModeMetric, ModeMetric, float]] = []
    index = 0
    while index + 1 < len(bending_modes):
        first = bending_modes[index]
        second = bending_modes[index + 1]
        mean_lambda = 0.5 * (first.lambda_fem + second.lambda_fem)
        split = abs(second.lambda_fem - first.lambda_fem) / mean_lambda if mean_lambda > 0.0 else math.inf
        if split <= BENDING_DOUBLET_REL_TOL:
            pairs.append((first, second, split))
            index += 2
        else:
            index += 1
    return pairs


def build_bending_doublet_comparison_rows(
    reference_rows: Sequence[dict[str, float | int]],
    mode_metrics: Sequence[ModeMetric],
) -> list[BendingDoubletComparison]:
    references = reference_lookup(reference_rows)
    metrics_by_epsilon: dict[float, list[ModeMetric]] = defaultdict(list)
    for metric in mode_metrics:
        metrics_by_epsilon[float(metric.epsilon)].append(metric)

    rows: list[BendingDoubletComparison] = []
    for epsilon in EPSILON_VALUES:
        epsilon_key = float(epsilon)
        pairs = pair_bending_modes(metrics_by_epsilon.get(epsilon_key, []))
        for mode_index in range(1, REFERENCE_1D_MODES + 1):
            reference = references[(epsilon_key, mode_index)]
            if mode_index <= len(pairs):
                first, second, split = pairs[mode_index - 1]
                mean_lambda = 0.5 * (first.lambda_fem + second.lambda_fem)
                rows.append(
                    BendingDoubletComparison(
                        epsilon=epsilon_key,
                        bending_doublet_index=mode_index,
                        mode_a=first.solid_mode,
                        mode_b=second.solid_mode,
                        class_a=first.classification,
                        class_b=second.classification,
                        lambda_a=first.lambda_fem,
                        lambda_b=second.lambda_fem,
                        lambda_fem_mean=mean_lambda,
                        lambda_split_relative=split,
                        lambda_eb=float(reference["Lambda_EB"]),
                        lambda_timoshenko=float(reference["Lambda_Timoshenko"]),
                        fem_minus_eb=mean_lambda - float(reference["Lambda_EB"]),
                        fem_minus_timoshenko=mean_lambda - float(reference["Lambda_Timoshenko"]),
                        status="classified_bending_doublet",
                        notes=(
                            "Paired from close bending-like solid modes; mode identity remains diagnostic "
                            "until shape plots or MAC checks are added."
                        ),
                    )
                )
            else:
                rows.append(
                    BendingDoubletComparison(
                        epsilon=epsilon_key,
                        bending_doublet_index=mode_index,
                        mode_a=None,
                        mode_b=None,
                        class_a="",
                        class_b="",
                        lambda_a=None,
                        lambda_b=None,
                        lambda_fem_mean=None,
                        lambda_split_relative=None,
                        lambda_eb=float(reference["Lambda_EB"]),
                        lambda_timoshenko=float(reference["Lambda_Timoshenko"]),
                        fem_minus_eb=None,
                        fem_minus_timoshenko=None,
                        status="not_classified",
                        notes="No close classified bending doublet was available for this 1D mode index.",
                    )
                )
    return rows


def write_comparison_csv(rows: Sequence[dict[str, object]]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "epsilon",
        "mode_1d",
        "Omega_EB",
        "Lambda_EB",
        "Omega_Timoshenko",
        "Lambda_Timoshenko",
        "solid_sorted_mode",
        "solid_doublet_member",
        "Omega_FEM",
        "Lambda_FEM",
        "status",
        "notes",
    ]
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_mode_metrics_csv(rows: Sequence[ModeMetric]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "epsilon",
        "solid_mode",
        "Omega_FEM",
        "Lambda_FEM",
        "parsed_node_count",
        "axial_energy_fraction",
        "transverse_y_energy_fraction",
        "transverse_z_energy_fraction",
        "transverse_energy_fraction",
        "mean_axial_energy_fraction",
        "mean_transverse_energy_fraction",
        "torsion_indicator",
        "bending_indicator",
        "warping_indicator",
        "dominant_transverse_axis",
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
                    "axial_energy_fraction": row.axial_energy_fraction,
                    "transverse_y_energy_fraction": row.transverse_y_energy_fraction,
                    "transverse_z_energy_fraction": row.transverse_z_energy_fraction,
                    "transverse_energy_fraction": row.transverse_energy_fraction,
                    "mean_axial_energy_fraction": row.mean_axial_energy_fraction,
                    "mean_transverse_energy_fraction": row.mean_transverse_energy_fraction,
                    "torsion_indicator": row.torsion_indicator,
                    "bending_indicator": row.bending_indicator,
                    "warping_indicator": row.warping_indicator,
                    "dominant_transverse_axis": row.dominant_transverse_axis,
                    "classification": row.classification,
                    "classification_note": row.classification_note,
                }
            )


def write_bending_doublet_csv(rows: Sequence[BendingDoubletComparison]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "epsilon",
        "bending_doublet_index",
        "solid_mode_a",
        "solid_mode_b",
        "classification_a",
        "classification_b",
        "Lambda_FEM_a",
        "Lambda_FEM_b",
        "Lambda_FEM_mean",
        "Lambda_split_relative",
        "Lambda_EB",
        "Lambda_Timoshenko",
        "FEM_minus_EB",
        "FEM_minus_Timoshenko",
        "status",
        "notes",
    ]
    with OUTPUT_BENDING_DOUBLET_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "epsilon": row.epsilon,
                    "bending_doublet_index": row.bending_doublet_index,
                    "solid_mode_a": row.mode_a if row.mode_a is not None else "",
                    "solid_mode_b": row.mode_b if row.mode_b is not None else "",
                    "classification_a": row.class_a,
                    "classification_b": row.class_b,
                    "Lambda_FEM_a": row.lambda_a if row.lambda_a is not None else "",
                    "Lambda_FEM_b": row.lambda_b if row.lambda_b is not None else "",
                    "Lambda_FEM_mean": row.lambda_fem_mean if row.lambda_fem_mean is not None else "",
                    "Lambda_split_relative": row.lambda_split_relative if row.lambda_split_relative is not None else "",
                    "Lambda_EB": row.lambda_eb,
                    "Lambda_Timoshenko": row.lambda_timoshenko,
                    "FEM_minus_EB": row.fem_minus_eb if row.fem_minus_eb is not None else "",
                    "FEM_minus_Timoshenko": row.fem_minus_timoshenko if row.fem_minus_timoshenko is not None else "",
                    "status": row.status,
                    "notes": row.notes,
                }
            )


def tool_text(value: str | None) -> str:
    return f"`{value}`" if value else "not found"


def markdown_reference_table(rows: Sequence[dict[str, object]]) -> str:
    lines = [
        "| epsilon | 1D mode | Lambda_EB | Lambda_Timoshenko | solid mode | Lambda_FEM | status |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        lambda_fem = row["Lambda_FEM"]
        solid_mode = row["solid_sorted_mode"]
        solid_mode_text = str(solid_mode) if solid_mode != "" else "-"
        lambda_fem_text = f"{float(lambda_fem):.8g}" if lambda_fem != "" else "-"
        lines.append(
            "| "
            f"{float(row['epsilon']):.5g} | "
            f"{int(row['mode_1d'])} | "
            f"{float(row['Lambda_EB']):.8g} | "
            f"{float(row['Lambda_Timoshenko']):.8g} | "
            f"{solid_mode_text} | "
            f"{lambda_fem_text} | "
            f"{row['status']} |"
        )
    return "\n".join(lines)


def markdown_mode_metrics_table(rows: Sequence[ModeMetric]) -> str:
    if not rows:
        return "No 3D solid mode-shape metrics were computed."
    lines = [
        "| epsilon | solid mode | Lambda_FEM | class | axial frac | trans frac | mean trans frac | torsion | warping | axis | nodes |",
        "| ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | --- | ---: |",
    ]
    for row in sorted(rows, key=lambda item: (item.epsilon, item.solid_mode)):
        lines.append(
            "| "
            f"{row.epsilon:.5g} | "
            f"{row.solid_mode} | "
            f"{row.lambda_fem:.8g} | "
            f"{row.classification} | "
            f"{row.axial_energy_fraction:.3f} | "
            f"{row.transverse_energy_fraction:.3f} | "
            f"{row.mean_transverse_energy_fraction:.3f} | "
            f"{row.torsion_indicator:.3f} | "
            f"{row.warping_indicator:.3f} | "
            f"{row.dominant_transverse_axis} | "
            f"{row.parsed_node_count} |"
        )
    return "\n".join(lines)


def markdown_bending_doublet_table(rows: Sequence[BendingDoubletComparison]) -> str:
    if not rows:
        return "No classified bending doublets were available."
    lines = [
        "| epsilon | bending index | solid modes | classes | Lambda_EB | Lambda_Timoshenko | Lambda_FEM avg | split | FEM-EB | FEM-Timo | status |",
        "| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        mode_text = f"{row.mode_a}, {row.mode_b}" if row.mode_a is not None and row.mode_b is not None else "-"
        class_text = f"{row.class_a}, {row.class_b}" if row.class_a and row.class_b else "-"
        lambda_mean = f"{row.lambda_fem_mean:.8g}" if row.lambda_fem_mean is not None else "-"
        split = f"{row.lambda_split_relative:.3g}" if row.lambda_split_relative is not None else "-"
        fem_minus_eb = f"{row.fem_minus_eb:.3g}" if row.fem_minus_eb is not None else "-"
        fem_minus_timo = f"{row.fem_minus_timoshenko:.3g}" if row.fem_minus_timoshenko is not None else "-"
        lines.append(
            "| "
            f"{row.epsilon:.5g} | "
            f"{row.bending_doublet_index} | "
            f"{mode_text} | "
            f"{class_text} | "
            f"{row.lambda_eb:.8g} | "
            f"{row.lambda_timoshenko:.8g} | "
            f"{lambda_mean} | "
            f"{split} | "
            f"{fem_minus_eb} | "
            f"{fem_minus_timo} | "
            f"{row.status} |"
        )
    return "\n".join(lines)


def classification_counts(rows: Sequence[ModeMetric]) -> list[str]:
    if not rows:
        return ["- no parsed mode shapes"]
    grouped: dict[float, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for row in rows:
        grouped[row.epsilon][row.classification] += 1
    lines: list[str] = []
    for epsilon in sorted(grouped):
        parts = ", ".join(f"{name}: {count}" for name, count in sorted(grouped[epsilon].items()))
        lines.append(f"- epsilon={epsilon:g}: {parts}")
    return lines


def doublet_observation_lines(rows: Sequence[BendingDoubletComparison]) -> list[str]:
    classified = [row for row in rows if row.status == "classified_bending_doublet"]
    if not classified:
        return ["- no classified bending doublets were available for EB/Timoshenko comparison"]
    closer_to_timo = sum(
        1
        for row in classified
        if row.fem_minus_timoshenko is not None
        and row.fem_minus_eb is not None
        and abs(row.fem_minus_timoshenko) < abs(row.fem_minus_eb)
    )
    lines = [
        (
            f"- {closer_to_timo}/{len(classified)} classified doublets are closer to "
            "the Timoshenko Lambda than to the Euler--Bernoulli Lambda."
        )
    ]
    missing = [row for row in rows if row.status != "classified_bending_doublet"]
    if missing:
        labels = ", ".join(f"epsilon={row.epsilon:g} mode={row.bending_doublet_index}" for row in missing)
        lines.append(
            "- classified doublet comparison is unavailable for "
            f"{labels}; within the requested 16 solid modes, bending candidates are interleaved "
            "with axial/torsion-like modes or not enough close bending candidates remain."
        )
    return lines


def fmt_optional(value: float | int | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.8g}"


def fmt_bbox(value: tuple[float, float, float] | None) -> str:
    if value is None:
        return "-"
    return f"({value[0]:.6g}, {value[1]:.6g}, {value[2]:.6g})"


def markdown_mesh_summary_table(mesh_summaries: Sequence[MeshSummary]) -> str:
    if not mesh_summaries:
        return "No generated mesh files were available for sanity summary."
    lines = [
        "| epsilon | file | nodes | 3D elems | left elems/nodes | right elems/nodes | bbox min | bbox max | length | radius est | dims ok | end groups |",
        "| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for summary in mesh_summaries:
        left_text = f"{fmt_optional(summary.left_face_element_count)}/{fmt_optional(summary.left_face_node_count)}"
        right_text = f"{fmt_optional(summary.right_face_element_count)}/{fmt_optional(summary.right_face_node_count)}"
        lines.append(
            "| "
            f"{summary.epsilon:.5g} | "
            f"`{summary.file_path.relative_to(REPO_ROOT)}` | "
            f"{summary.node_count} | "
            f"{summary.solid_element_count} | "
            f"{left_text} | "
            f"{right_text} | "
            f"{fmt_bbox(summary.bbox_min)} | "
            f"{fmt_bbox(summary.bbox_max)} | "
            f"{fmt_optional(summary.inferred_length)} | "
            f"{fmt_optional(summary.inferred_radius)} | "
            f"{summary.dimensions_match_expected} | "
            f"{summary.end_face_groups_present} |"
        )
    return "\n".join(lines)


def mesh_warning_lines(mesh_summaries: Sequence[MeshSummary]) -> list[str]:
    lines: list[str] = []
    for summary in mesh_summaries:
        if not summary.warnings:
            continue
        label = summary.file_path.relative_to(REPO_ROOT)
        for warning in summary.warnings:
            lines.append(f"- `{label}`: {warning}")
    if not lines:
        return ["No Gmsh warnings/errors or mesh sanity warnings were recorded."]
    return lines


def markdown_ccx_input_table(ccx_inputs: Sequence[CcxInputSummary]) -> str:
    if not ccx_inputs:
        return "No CalculiX modal input files were generated."
    lines = [
        "| epsilon | modal input | mesh include | solid elems | LEFT_FIXED nodes | RIGHT_FIXED nodes | node sets present |",
        "| ---: | --- | --- | ---: | ---: | ---: | --- |",
    ]
    for item in ccx_inputs:
        lines.append(
            "| "
            f"{item.epsilon:.5g} | "
            f"`{item.modal_input_path.relative_to(REPO_ROOT)}` | "
            f"`{item.mesh_include_path.relative_to(REPO_ROOT)}` | "
            f"{item.solid_element_count} | "
            f"{item.left_fixed_node_count} | "
            f"{item.right_fixed_node_count} | "
            f"{item.node_sets_present} |"
        )
    return "\n".join(lines)


def calculix_result_lines(results: Sequence[CalculixResult]) -> list[str]:
    if not results:
        return ["CalculiX was not run."]
    lines: list[str] = []
    for result in results:
        status = "succeeded" if result.success else "failed"
        lines.append(f"- epsilon={result.epsilon:g}: {status}; {result.message}")
        lines.append(f"  stdout: `{result.stdout_path.relative_to(REPO_ROOT)}`")
        lines.append(f"  stderr: `{result.stderr_path.relative_to(REPO_ROOT)}`")
        lines.append(f"  dat: `{result.dat_path.relative_to(REPO_ROOT)}`")
        lines.append(f"  frd: `{result.frd_path.relative_to(REPO_ROOT)}`")
        lines.append(f"  parse source: {result.parse_source}")
    return lines


def write_report(
    audit: ToolAudit,
    comparison_rows: Sequence[dict[str, object]],
    mode_metrics: Sequence[ModeMetric],
    bending_doublet_rows: Sequence[BendingDoubletComparison],
    case_messages: Sequence[str],
    mesh_summaries: Sequence[MeshSummary],
    ccx_inputs: Sequence[CcxInputSummary],
    calculix_results: Sequence[CalculixResult],
    mode_shape_parse_results: Sequence[ModeShapeParseResult],
) -> None:
    computed_count = sum(1 for row in comparison_rows if row["status"] != "not_computed")
    any_computed = computed_count > 0
    classified_doublet_count = sum(1 for row in bending_doublet_rows if row.status == "classified_bending_doublet")
    lines: list[str] = []
    lines.append("# Single Fixed-Fixed Rod 3D Solid FEM Diagnostic")
    lines.append("")
    lines.append("Diagnostic only. This workflow seeds an independent 3D solid FEM check for the one-rod EB/Timoshenko comparison.")
    lines.append("It does not change the article, paper_dorofeev_style, article figures, the old determinant, formulas.py, old solvers, or the existing FEM physical model.")
    lines.append("")
    lines.append("## Tool Audit")
    lines.append("")
    lines.append(f"- Gmsh executable: {tool_text(audit.gmsh_exe)}")
    lines.append(f"- Gmsh resolution source: {audit.gmsh_source if audit.gmsh_source else 'not found'}")
    lines.append(f"- Gmsh version: {audit.gmsh_version if audit.gmsh_version else 'not available'}")
    lines.append(f"- Gmsh Python API: {'found' if audit.gmsh_python_api else 'not found'}")
    lines.append(f"- meshio Python API: {'found' if audit.meshio_python_api else 'not found'}")
    lines.append(f"- CalculiX ccx: {tool_text(audit.calculix_ccx)}")
    lines.append(f"- CalculiX resolution source: {audit.calculix_source if audit.calculix_source else 'not found'}")
    lines.append(f"- CalculiX version probe: {audit.calculix_version if audit.calculix_version else 'not available'}")
    lines.append(f"- Code_Aster as_run: {tool_text(audit.code_aster_as_run)}")
    lines.append(f"- Code_Aster run_aster: {tool_text(audit.code_aster_run_aster)}")
    lines.append(f"- Salome-Meca salome: {tool_text(audit.salome)}")
    lines.append(f"- Salome-Meca salome_meca: {tool_text(audit.salome_meca)}")
    lines.append("")
    lines.append("## Chosen Path")
    lines.append("")
    lines.append("The preferred lightweight path is Gmsh cylinder mesh generation followed by a CalculiX modal frequency run.")
    lines.append("The script writes Gmsh `.geo` files and CalculiX modal templates under `results/solid_fem_single_rod/` for each epsilon.")
    if audit.has_mesh_generator and audit.calculix_ccx:
        lines.append("Both a mesh generator and CalculiX were detected, so the script attempted actual modal runs.")
    elif audit.has_mesh_generator:
        lines.append("A mesh generator was detected, but no CalculiX solver was found, so only mesh/input generation was attempted.")
    else:
        lines.append("No Gmsh executable or Python API was detected, so this run only wrote input templates and analytic reference rows.")
    lines.append("")
    lines.append("## Benchmark Definition")
    lines.append("")
    lines.append("- geometry: straight circular cylinder")
    lines.append(f"- length: `L = {L:g}`")
    lines.append("- radius: `r = epsilon*L`")
    lines.append(f"- epsilons: {', '.join(f'{value:g}' for value in EPSILON_VALUES)}")
    lines.append(f"- material: `E = {E:g}`, `rho = {RHO:g}`, `nu = {NU:g}`")
    lines.append("- boundary conditions: both end faces fully fixed, `Ux = Uy = Uz = 0`")
    lines.append(f"- requested solid modes: {SOLID_MODES_REQUESTED}")
    lines.append("")
    lines.append("3D solid circular rods should produce bending doublets because the two transverse bending planes are equivalent. Small numerical splitting is expected and should not be overinterpreted.")
    lines.append("")
    lines.append("## Normalization")
    lines.append("")
    lines.append("The 3D FEM angular frequency `Omega` is converted to the same single-rod project normalization used in the EB/Timoshenko diagnostic:")
    lines.append("")
    lines.append("```text")
    lines.append("Lambda_FEM = sqrt(Omega)*(L/2)*(rho*A/(E*I))**0.25")
    lines.append("```")
    lines.append("")
    lines.append("## Result Status")
    lines.append("")
    if any_computed:
        lines.append(f"Computed {computed_count} preliminary solid sorted-mode rows.")
        lines.append(f"Computed {len(mode_metrics)} mode-shape metric rows and {classified_doublet_count} classified bending doublets.")
    else:
        lines.append("No 3D solid modal frequencies were computed in this run.")
        lines.append("No fake FEM frequency rows were created; FEM columns in the CSV are blank with `status = not_computed`.")
        if audit.has_mesh_generator and not audit.calculix_ccx:
            lines.append("Mesh generation can still succeed; modal solve remains pending until CalculiX `ccx` is available.")
    lines.append("")
    lines.append("## Mesh Sanity Summary")
    lines.append("")
    lines.append(markdown_mesh_summary_table(mesh_summaries))
    lines.append("")
    lines.append("Expected dimensions are `L = 1` and `r = epsilon*L`. Radius is estimated as `max(sqrt(y^2 + z^2))` over mesh nodes.")
    lines.append("")
    lines.append("Mesh sanity warnings and captured Gmsh warning/error lines:")
    lines.append("")
    lines.append("Raw Gmsh export warnings are listed separately from the generated CalculiX node-set summary below.")
    lines.append("")
    lines.extend(mesh_warning_lines(mesh_summaries))
    lines.append("")
    lines.append("## CalculiX Input Summary")
    lines.append("")
    lines.append(markdown_ccx_input_table(ccx_inputs))
    lines.append("")
    lines.append("Coordinate-derived node sets are used for boundary conditions:")
    lines.append("")
    lines.append("```text")
    lines.append("*BOUNDARY")
    lines.append("LEFT_FIXED, 1, 3, 0")
    lines.append("RIGHT_FIXED, 1, 3, 0")
    lines.append("```")
    lines.append("")
    lines.append("## CalculiX Run Summary")
    lines.append("")
    lines.extend(calculix_result_lines(calculix_results))
    lines.append("")
    if calculix_results:
        parse_sources = sorted({result.parse_source for result in calculix_results if result.parsed_omegas})
        if parse_sources:
            lines.append("Recognized frequency convention:")
            for source in parse_sources:
                lines.append(f"- {source}")
            lines.append("")
    lines.append("## Mode-Shape Classification")
    lines.append("")
    lines.append("The sorted comparison below is retained as a preliminary frequency-only view. For bending-mode comparisons, the classified doublet table is preferred because higher sorted 3D modes can include torsion-like, axial-like, or mixed solid modes.")
    lines.append("")
    lines.append("Mode-shape metrics are approximate nodal diagnostics. Nodes are grouped into axial slices; each slice is reduced to mean `Ux, Uy, Uz`, a torsion-like circumferential fit about the x-axis, and a residual within-slice warping/non-planarity measure after subtracting the slice mean translation.")
    lines.append("")
    lines.extend(mode_shape_parse_lines(mode_shape_parse_results))
    lines.append("")
    lines.append("Classification counts:")
    lines.append("")
    lines.extend(classification_counts(mode_metrics))
    lines.append("")
    lines.append(markdown_mode_metrics_table(mode_metrics))
    lines.append("")
    lines.append("## Classified Bending Doublets")
    lines.append("")
    lines.append("Close bending-like modes are paired in frequency order, their `Lambda_FEM` values are averaged, and the averaged doublet is compared with the 1D EB/Timoshenko references. The pairing is diagnostic only; circular-solid doublets may rotate within the y/z subspace and require shape plots or MAC checks for stronger identity claims.")
    lines.append("")
    lines.extend(doublet_observation_lines(bending_doublet_rows))
    lines.append("")
    lines.append(markdown_bending_doublet_table(bending_doublet_rows))
    lines.append("")
    lines.append("## Comparison Table")
    lines.append("")
    lines.append("This is the older sorted-mode comparison and is intentionally marked preliminary.")
    lines.append("")
    lines.append(markdown_reference_table(comparison_rows))
    lines.append("")
    lines.append("## Case Messages")
    lines.append("")
    if case_messages:
        for message in case_messages:
            lines.append(f"- {message}")
    else:
        lines.append("- No case-specific messages.")
    lines.append("")
    lines.append("## Next Manual Step")
    lines.append("")
    if not audit.has_mesh_generator:
        lines.append("Install or point `GMSH_EXE` to a portable Gmsh executable, then rerun:")
    elif not audit.calculix_ccx:
        lines.append("Mesh generation is available. Install or expose CalculiX `ccx`, then rerun:")
    else:
        lines.append("Review the generated CalculiX `.dat` files and mode doublets, then extend the parser or shape checks if needed:")
    lines.append("")
    lines.append("```text")
    lines.append("python scripts/analysis/solid_fem_single_rod_fixed_fixed.py")
    lines.append("```")
    lines.append("")
    lines.append("Manual Gmsh generation example for one case:")
    lines.append("")
    lines.append("```text")
    lines.append("gmsh results/solid_fem_single_rod/eps_0p025/single_rod_fixed_fixed_eps0p025.geo -3 -format inp -o results/solid_fem_single_rod/eps_0p025/single_rod_fixed_fixed_eps0p025_mesh.inp")
    lines.append("```")
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    lines.append(f"- comparison CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`")
    lines.append(f"- mode metrics CSV: `{OUTPUT_MODE_METRICS_CSV.relative_to(REPO_ROOT)}`")
    lines.append(f"- bending doublet CSV: `{OUTPUT_BENDING_DOUBLET_CSV.relative_to(REPO_ROOT)}`")
    lines.append(f"- report: `{OUTPUT_REPORT.relative_to(REPO_ROOT)}`")
    lines.append("- input templates: `results/solid_fem_single_rod/`")
    lines.append("")
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def process_cases(
    audit: ToolAudit,
) -> tuple[
    dict[float, list[float]],
    list[ModeMetric],
    list[ModeShapeParseResult],
    list[str],
    list[MeshSummary],
    list[CcxInputSummary],
    list[CalculixResult],
]:
    fem_omegas_by_epsilon: dict[float, list[float]] = {}
    mode_metrics: list[ModeMetric] = []
    mode_shape_parse_results: list[ModeShapeParseResult] = []
    messages: list[str] = []
    mesh_summaries: list[MeshSummary] = []
    ccx_inputs: list[CcxInputSummary] = []
    calculix_results: list[CalculixResult] = []
    for epsilon in EPSILON_VALUES:
        paths = case_paths(float(epsilon))
        write_gmsh_geo(float(epsilon), paths)
        messages.append(f"epsilon={float(epsilon):g}: wrote {paths.geo.relative_to(REPO_ROOT)}")

        mesh_generated = False
        gmsh_messages: tuple[str, ...] = ()
        if audit.gmsh_exe:
            mesh_generated, message, gmsh_messages = generate_mesh_with_gmsh_cli(paths, audit.gmsh_exe)
            messages.append(f"epsilon={float(epsilon):g}: {message}")
        elif audit.gmsh_python_api:
            mesh_generated, message, gmsh_messages = generate_mesh_with_gmsh_python(paths)
            messages.append(f"epsilon={float(epsilon):g}: {message}")
        else:
            messages.append(
                f"epsilon={float(epsilon):g}: no Gmsh executable/API detected; mesh was not generated."
            )

        mesh_summaries.extend(summarize_generated_meshes(float(epsilon), paths, gmsh_messages))

        if not mesh_generated:
            continue
        ccx_input = write_calculix_inputs(paths, float(epsilon))
        ccx_inputs.append(ccx_input)
        messages.append(f"epsilon={float(epsilon):g}: wrote {paths.ccx_mesh_inp.relative_to(REPO_ROOT)}")
        messages.append(f"epsilon={float(epsilon):g}: wrote {paths.ccx_modal_inp.relative_to(REPO_ROOT)}")
        messages.append(
            f"epsilon={float(epsilon):g}: generated LEFT_FIXED nodes={ccx_input.left_fixed_node_count}, "
            f"RIGHT_FIXED nodes={ccx_input.right_fixed_node_count}"
        )
        if not audit.calculix_ccx:
            messages.append(
                f"epsilon={float(epsilon):g}: mesh generated successfully; "
                "modal solve pending until CalculiX ccx is available."
            )
            continue

        result = run_calculix(paths, audit.calculix_ccx, float(epsilon))
        calculix_results.append(result)
        messages.append(f"epsilon={float(epsilon):g}: {result.message}")
        if result.success and result.parsed_omegas:
            fem_omegas_by_epsilon[float(epsilon)] = list(result.parsed_omegas)
            messages.append(
                f"epsilon={float(epsilon):g}: parsed {len(result.parsed_omegas)} eigenfrequency rows."
            )
            metrics, parse_result = compute_mode_metrics_for_case(
                epsilon=float(epsilon),
                paths=paths,
                omegas=result.parsed_omegas,
            )
            mode_shape_parse_results.append(parse_result)
            mode_metrics.extend(metrics)
            messages.append(f"epsilon={float(epsilon):g}: {parse_result.message}")
            if metrics:
                messages.append(f"epsilon={float(epsilon):g}: computed {len(metrics)} mode-shape metric rows.")
    return (
        fem_omegas_by_epsilon,
        mode_metrics,
        mode_shape_parse_results,
        messages,
        mesh_summaries,
        ccx_inputs,
        calculix_results,
    )


def main() -> int:
    audit = audit_tools()
    (
        fem_omegas_by_epsilon,
        mode_metrics,
        mode_shape_parse_results,
        case_messages,
        mesh_summaries,
        ccx_inputs,
        calculix_results,
    ) = process_cases(audit)
    reference_rows = compute_1d_references()
    comparison_rows = build_comparison_rows(reference_rows, fem_omegas_by_epsilon)
    bending_doublet_rows = build_bending_doublet_comparison_rows(reference_rows, mode_metrics)
    write_comparison_csv(comparison_rows)
    write_mode_metrics_csv(mode_metrics)
    write_bending_doublet_csv(bending_doublet_rows)
    write_report(
        audit,
        comparison_rows,
        mode_metrics,
        bending_doublet_rows,
        case_messages,
        mesh_summaries,
        ccx_inputs,
        calculix_results,
        mode_shape_parse_results,
    )

    print(f"wrote comparison CSV: {OUTPUT_CSV}")
    print(f"wrote mode metrics CSV: {OUTPUT_MODE_METRICS_CSV}")
    print(f"wrote bending doublet CSV: {OUTPUT_BENDING_DOUBLET_CSV}")
    print(f"wrote report: {OUTPUT_REPORT}")
    print(f"wrote input templates under: {SOLID_OUTPUT_DIR}")
    if not any(row["status"] != "not_computed" for row in comparison_rows):
        if audit.has_mesh_generator and not audit.calculix_ccx:
            print("mesh generated where possible; no 3D solid modal results computed because CalculiX ccx is missing")
        else:
            print("no 3D solid modal results computed; missing mesh generator and/or solver")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
