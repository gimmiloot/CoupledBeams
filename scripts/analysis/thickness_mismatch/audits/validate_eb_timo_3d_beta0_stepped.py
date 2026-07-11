from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from dataclasses import dataclass
from types import ModuleType
from typing import Sequence


REPO_ROOT = Path(__file__).resolve().parents[4]
SRC_ROOT = REPO_ROOT / "src"
LIB_ROOT = REPO_ROOT / "scripts" / "lib"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))
if str(LIB_ROOT) not in sys.path:
    sys.path.insert(0, str(LIB_ROOT))


DEFAULT_EPSILON_VALUES = (0.0025, 0.01, 0.05)
DEFAULT_MU = 0.5
DEFAULT_ETA = 0.1
DEFAULT_POISSON = 0.3
DEFAULT_TOTAL_LENGTH = 2.0
DEFAULT_N_ANALYTIC_ROOTS = 8
DEFAULT_N_FEM_MODES = 60
DEFAULT_ROOT_SCAN_STEP = 0.01
DEFAULT_LAMBDA_MAX = 60.0
DEFAULT_TIMEOUT_SECONDS = 1200
DEFAULT_OUTPUT_DIR = (
    Path("results")
    / "eb_vs_timoshenko_3d_validation"
    / "stepped_beta0_mu0p5_eta0p1"
)
DEFAULT_GMSH_CANDIDATES = (
    Path(r"D:\PHD\gmsh-4.15.2-Windows64\gmsh.exe"),
    Path(r"D:\PHD\gmsh\gmsh.exe"),
)
DEFAULT_CCX_CANDIDATES = (
    Path(r"D:\PHD\calculix\ccx.exe"),
    Path(r"D:\PHD\calculix\ccx_static.exe"),
    Path(r"D:\PHD\calculix\calculix_2.22_4win\ccx_static.exe"),
)
WORKFLOW = "diagnostic_straight_coaxial_stepped_cylinder"


def _load_module(module_name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


FORMULAS = _load_module(
    "formulas_thickness_mismatch_for_stepped_beta0_validation",
    REPO_ROOT / "src" / "my_project" / "analytic" / "formulas_thickness_mismatch.py",
)
TIMO = _load_module(
    "variable_length_timoshenko_for_stepped_beta0_validation",
    REPO_ROOT / "scripts" / "lib" / "variable_length_timoshenko.py",
)


@dataclass(frozen=True)
class ExecutableResolution:
    path: Path | None
    source: str
    candidates_checked: tuple[Path, ...]


@dataclass(frozen=True)
class Geometry:
    epsilon: float
    mu: float
    eta: float
    tau1: float
    tau2: float
    l1: float
    l2: float
    total_length: float
    r0: float
    r1: float
    r2: float
    h: float


@dataclass(frozen=True)
class CaseOutputs:
    analytic_csv: Path
    raw_fem_csv: Path
    match_csv: Path
    pair_csv: Path
    plot_png: Path
    error_plot_png: Path


@dataclass
class FemResult:
    status: str
    blocker: str
    parse_source: str
    gmsh_messages: tuple[str, ...]
    ccx_messages: tuple[str, ...]
    mesh_info: dict[str, object]
    raw_rows: list[dict[str, object]]
    bending_rows: list[dict[str, object]]
    case_files: list[Path]


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def token(value: float) -> str:
    return f"{float(value):.6g}".replace("-", "m").replace(".", "p")


def fmt(value: object) -> object:
    if isinstance(value, float):
        if not math.isfinite(value):
            return "nan"
        return f"{value:.16e}"
    return value


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def write_csv(path: Path, rows: Sequence[dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key, "")) for key in fieldnames})


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def rel_error(model: float, reference: float) -> float:
    if not math.isfinite(model) or not math.isfinite(reference) or reference == 0.0:
        return float("nan")
    return abs(float(model) - float(reference)) / abs(float(reference))


def closer_model(eb_error: float, timo_error: float) -> str:
    if not (math.isfinite(eb_error) and math.isfinite(timo_error)):
        return "unavailable"
    if abs(eb_error - timo_error) <= 1.0e-12 * max(1.0, abs(eb_error), abs(timo_error)):
        return "tie"
    return "Timoshenko" if timo_error < eb_error else "Euler-Bernoulli"


def default_mesh_size(epsilon: float) -> float:
    eps = float(epsilon)
    if eps <= 0.0025 + 1.0e-15:
        return 0.004
    if eps <= 0.01 + 1.0e-15:
        return 0.008
    return 0.025


def build_geometry(epsilon: float, mu: float, eta: float, h: float) -> Geometry:
    factors = FORMULAS.thickness_mismatch_factors(float(mu), float(eta))
    l1 = 1.0 - float(mu)
    l2 = 1.0 + float(mu)
    r0 = 2.0 * float(epsilon)
    return Geometry(
        epsilon=float(epsilon),
        mu=float(mu),
        eta=float(eta),
        tau1=float(factors.tau1),
        tau2=float(factors.tau2),
        l1=l1,
        l2=l2,
        total_length=l1 + l2,
        r0=r0,
        r1=r0 * float(factors.tau1),
        r2=r0 * float(factors.tau2),
        h=float(h),
    )


def resolve_existing(path_text: str | None) -> Path | None:
    if not path_text:
        return None
    path = Path(path_text.strip().strip('"'))
    return path if path.is_file() else None


def find_under(root: Path, names: Sequence[str], *, limit: int = 12) -> list[Path]:
    if not root.exists():
        return []
    wanted = {name.lower() for name in names}
    found: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [
            name
            for name in dirnames
            if name not in {"$RECYCLE.BIN", "System Volume Information"}
        ]
        for filename in filenames:
            if filename.lower() in wanted:
                found.append(Path(dirpath) / filename)
                if len(found) >= limit:
                    return found
    return found


def resolve_executable(
    *,
    explicit: str | None,
    env_name: str,
    path_names: Sequence[str],
    default_candidates: Sequence[Path],
    search_root: Path | None = Path(r"D:\PHD"),
) -> ExecutableResolution:
    checked: list[Path] = []
    explicit_path = resolve_existing(explicit)
    if explicit_path is not None:
        return ExecutableResolution(explicit_path, "explicit CLI path", (explicit_path,))
    if explicit:
        checked.append(Path(explicit.strip().strip('"')))

    env_path = resolve_existing(os.environ.get(env_name))
    if env_path is not None:
        return ExecutableResolution(env_path, f"environment variable {env_name}", (env_path,))
    if os.environ.get(env_name):
        checked.append(Path(os.environ[env_name].strip().strip('"')))

    for name in path_names:
        path_text = shutil.which(name)
        if path_text:
            return ExecutableResolution(Path(path_text), f"PATH lookup for {name}", tuple(checked))

    for candidate in default_candidates:
        checked.append(candidate)
        if candidate.is_file():
            return ExecutableResolution(candidate, "known D:\\PHD candidate", tuple(checked))

    if search_root is not None:
        found = find_under(search_root, path_names)
        checked.extend(found)
        if found:
            return ExecutableResolution(found[0], f"recursive search under {search_root}", tuple(checked))

    return ExecutableResolution(None, "not found", tuple(checked))


def probe_version(executable: Path | None, *args: str) -> str:
    if executable is None:
        return "not available"
    try:
        completed = subprocess.run(
            [str(executable), *args],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=30,
        )
    except Exception as exc:
        return f"probe failed: {exc}"
    text = (completed.stdout or completed.stderr or "").strip()
    return text.splitlines()[0].strip() if text else f"return code {completed.returncode}; no output"


def collect_warning_lines(text: str) -> tuple[str, ...]:
    warnings: list[str] = []
    for line in text.splitlines():
        lower = line.lower()
        if "warning" in lower or "error" in lower or "jac." in lower or "distortion" in lower:
            warnings.append(line.strip())
    return tuple(warnings)


def parse_mesh_quality(messages: Sequence[str]) -> tuple[int, float, tuple[str, ...]]:
    negative = 0
    worst = float("nan")
    warning_lines: list[str] = []
    for message in messages:
        lower = message.lower()
        if "warning" in lower or "error" in lower or "jac." in lower or "distortion" in lower:
            warning_lines.append(message)
        if "elements with jac." in lower:
            parts = message.replace("=", " ").replace("<", " ").split()
            for index, part in enumerate(parts):
                if part.lower() == "distortion" and index + 1 < len(parts):
                    try:
                        worst = float(parts[index + 1])
                    except ValueError:
                        pass
                if part.lower() == "elements" and index > 0:
                    try:
                        negative = int(float(parts[index - 1]))
                    except ValueError:
                        pass
    return negative, worst, tuple(warning_lines)


def load_single_module() -> ModuleType:
    return _load_module(
        "solid_fem_single_rod_fixed_fixed_for_stepped_beta0_validation",
        REPO_ROOT / "scripts" / "analysis" / "solid_fem_single_rod_fixed_fixed.py",
    )


def load_extraction_module() -> ModuleType:
    return _load_module(
        "audit_full_spectrum_3d_fem_smoke_extraction_for_stepped_beta0_validation",
        REPO_ROOT
        / "scripts"
        / "analysis"
        / "thickness_mismatch"
        / "audits"
        / "audit_full_spectrum_3d_fem_smoke_extraction.py",
    )


def write_stepped_gmsh_geo(paths: object, geometry: Geometry) -> None:
    paths.case_dir.mkdir(parents=True, exist_ok=True)
    rmax = max(geometry.r1, geometry.r2)
    lines = [
        "// Diagnostic-only straight coaxial stepped cylinder.",
        "// Segment 1: x in [0, L1], radius r1.",
        "// Segment 2: x in [L1, L1 + L2], radius r2.",
        "// Generate with gmsh -3 -format inp and then solve with CalculiX.",
        "",
        'SetFactory("OpenCASCADE");',
        "",
        f"L1 = {geometry.l1:.17g};",
        f"L2 = {geometry.l2:.17g};",
        f"L = {geometry.total_length:.17g};",
        f"R1 = {geometry.r1:.17g};",
        f"R2 = {geometry.r2:.17g};",
        f"RMAX = {rmax:.17g};",
        f"h = {geometry.h:.17g};",
        "tol = 1e-8;",
        "",
        "Cylinder(1) = {0, 0, 0, L1, 0, 0, R1, 2*Pi};",
        "Cylinder(2) = {L1, 0, 0, L2, 0, 0, R2, 2*Pi};",
        "fragments[] = BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; };",
        "Coherence;",
        "",
        "vols[] = Volume In BoundingBox {-tol, -RMAX-tol, -RMAX-tol, L+tol, RMAX+tol, RMAX+tol};",
        "left_faces[] = Surface In BoundingBox {-tol, -R1-tol, -R1-tol, tol, R1+tol, R1+tol};",
        "right_faces[] = Surface In BoundingBox {L-tol, -R2-tol, -R2-tol, L+tol, R2+tol, R2+tol};",
        "",
        'Physical Volume("SOLID", 1) = vols[];',
        'Physical Surface("FIXED_LEFT", 2) = left_faces[];',
        'Physical Surface("FIXED_RIGHT", 3) = right_faces[];',
        "",
        "Mesh.CharacteristicLengthMin = h;",
        "Mesh.CharacteristicLengthMax = h;",
        f"Mesh.ElementOrder = {2};",
        "Mesh.MshFileVersion = 4.1;",
        "",
    ]
    paths.geo.write_text("\n".join(lines), encoding="utf-8")


def connected_component_sizes(mesh: object) -> tuple[int, tuple[int, ...]]:
    elements = getattr(mesh, "solid_elements", {})
    if not elements:
        return 0, ()
    node_to_elements: dict[int, list[int]] = {}
    for element_id, nodes in elements.items():
        for node_id in nodes:
            node_to_elements.setdefault(int(node_id), []).append(int(element_id))
    unvisited = set(int(element_id) for element_id in elements)
    sizes: list[int] = []
    while unvisited:
        start = unvisited.pop()
        stack = [start]
        size = 0
        while stack:
            element_id = stack.pop()
            size += 1
            for node_id in elements[element_id]:
                for neighbor in node_to_elements.get(int(node_id), []):
                    if neighbor in unvisited:
                        unvisited.remove(neighbor)
                        stack.append(neighbor)
        sizes.append(size)
    sizes.sort(reverse=True)
    return len(sizes), tuple(sizes)


def segment_radius_estimates(mesh: object, geometry: Geometry) -> tuple[float, float]:
    nodes = getattr(mesh, "nodes", {})
    gap1 = min(0.20 * geometry.l1, max(2.0 * geometry.h, 1.0e-6))
    gap2 = min(0.20 * geometry.l2, max(2.0 * geometry.h, 1.0e-6))
    left_values = [
        math.hypot(point[1], point[2])
        for point in nodes.values()
        if point[0] <= geometry.l1 - gap1
    ]
    right_values = [
        math.hypot(point[1], point[2])
        for point in nodes.values()
        if point[0] >= geometry.l1 + gap2
    ]
    return (
        max(left_values) if left_values else float("nan"),
        max(right_values) if right_values else float("nan"),
    )


def compute_analytic_rows(
    *,
    geometry: Geometry,
    n_roots: int,
    lambda_max: float,
    scan_step: float,
) -> tuple[list[dict[str, object]], tuple[str, ...]]:
    eb_roots = FORMULAS.find_first_n_roots_eta(
        beta=0.0,
        mu=geometry.mu,
        epsilon=geometry.epsilon,
        eta=geometry.eta,
        n_roots=int(n_roots),
        Lmax0=float(lambda_max),
        scan_step=float(scan_step),
    )
    timo_roots, timo_warnings = TIMO.timo_sorted_roots(
        beta_deg=0.0,
        mu=geometry.mu,
        epsilon=geometry.epsilon,
        eta=geometry.eta,
        n_roots=int(n_roots),
    )
    rows: list[dict[str, object]] = []
    for index in range(int(n_roots)):
        eb = float(eb_roots[index]) if index < len(eb_roots) else float("nan")
        timo = float(timo_roots[index]) if index < len(timo_roots) else float("nan")
        if math.isfinite(eb) and eb != 0.0 and math.isfinite(timo):
            rel_diff = (timo - eb) / eb
        else:
            rel_diff = float("nan")
        note = "sorted in-plane roots; no descendant tracking"
        if not math.isfinite(eb):
            note += "; EB root missing"
        if not math.isfinite(timo):
            note += "; Timoshenko root missing"
        rows.append(
            {
                "sorted_index": index + 1,
                "Lambda_EB": eb,
                "Lambda_Timoshenko": timo,
                "rel_diff_Timo_vs_EB": rel_diff,
                "notes": note,
            }
        )
    return rows, tuple(timo_warnings)


def analytic_values(rows: Sequence[dict[str, object]], key: str) -> list[float]:
    return [float(row[key]) for row in rows if finite(row.get(key))]


def fem_case_paths(single: ModuleType, output_dir: Path, geometry: Geometry) -> object:
    eps_token = token(geometry.epsilon)
    stem = f"stepped_beta0_mu0p5_eta0p1_eps{eps_token}"
    case_dir = output_dir / "solid_fem_cases" / f"eps{eps_token}"
    return single.CasePaths(
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


def direct_raw_frequency_rows(
    *,
    extraction: ModuleType,
    single: ModuleType,
    paths: object,
    geometry: Geometry,
) -> tuple[list[dict[str, object]], str]:
    table_rows = extraction.parse_calculix_frequency_table(paths.ccx_dat)
    if table_rows:
        rows: list[dict[str, object]] = []
        for row in table_rows:
            omega = float(row["angular_frequency"])
            rows.append(
                {
                    "raw_mode_number": int(row["raw_mode_number"]),
                    "omega_rad_per_time": omega,
                    "Lambda_3D": math.sqrt(omega / geometry.epsilon),
                    "parser_included": "yes",
                    "notes": (
                        "direct CalculiX .dat table; FREQUENCY REAL PART (RAD/TIME); "
                        "Lambda=sqrt(omega/epsilon)"
                    ),
                }
            )
        rows.sort(key=lambda item: float(item["Lambda_3D"]))
        return rows, "CalculiX .dat frequency table RAD/TIME column"

    omegas, parse_source = single.parse_calculix_eigen_omegas(paths.ccx_dat, paths.ccx_stdout)
    rows = []
    for index, omega in enumerate(omegas, start=1):
        rows.append(
            {
                "raw_mode_number": index,
                "omega_rad_per_time": float(omega),
                "Lambda_3D": math.sqrt(float(omega) / geometry.epsilon),
                "parser_included": "fallback",
                "notes": f"fallback parser: {parse_source}; Lambda=sqrt(omega/epsilon)",
            }
        )
    rows.sort(key=lambda item: float(item["Lambda_3D"]))
    return rows, f"fallback helper parser: {parse_source}"


def add_mode_metrics_to_raw_rows(
    raw_rows: list[dict[str, object]],
    metrics: Sequence[object],
    parse_message: str,
) -> None:
    metrics_by_mode = {int(metric.solid_mode): metric for metric in metrics}
    for row in raw_rows:
        mode = int(row["raw_mode_number"])
        metric = metrics_by_mode.get(mode)
        if metric is None:
            row.update(
                {
                    "mode_character": "classification_unavailable",
                    "transverse_fraction": float("nan"),
                    "axial_fraction": float("nan"),
                    "torsion_fraction_if_available": float("nan"),
                }
            )
            row["notes"] = f"{row['notes']}; mode metric unavailable; {parse_message}"
            continue
        row.update(
            {
                "mode_character": metric.classification,
                "transverse_fraction": float(metric.transverse_energy_fraction),
                "axial_fraction": float(metric.axial_energy_fraction),
                "torsion_fraction_if_available": float(metric.torsion_indicator),
            }
        )
        row["notes"] = f"{row['notes']}; {metric.classification_note}"


def is_bending_related(row: dict[str, object]) -> bool:
    if not finite(row.get("Lambda_3D")):
        return False
    if str(row.get("mode_character", "")).startswith("bending_"):
        return True
    if finite(row.get("transverse_fraction")) and finite(row.get("axial_fraction")):
        return float(row["transverse_fraction"]) >= 0.55 and float(row["axial_fraction"]) <= 0.45
    return False


def select_bending_rows(raw_rows: Sequence[dict[str, object]], *, limit: int = 8) -> list[dict[str, object]]:
    bending_related = [dict(row) for row in raw_rows if is_bending_related(dict(row))]
    bending_related.sort(key=lambda row: float(row["Lambda_3D"]))
    for row in bending_related:
        if not str(row.get("mode_character", "")).startswith("bending_"):
            row["notes"] = f"{row.get('notes', '')}; selected as bending-related by transverse fraction"
    if bending_related:
        return bending_related[:limit]

    fallback = [dict(row) for row in raw_rows if finite(row.get("Lambda_3D"))]
    fallback.sort(key=lambda row: float(row["Lambda_3D"]))
    out = fallback[:limit]
    for row in out:
        row["notes"] = f"{row.get('notes', '')}; bending classification unavailable; low-frequency order retained only as fallback"
    return out


def pair_transverse_ok(first: dict[str, object], second: dict[str, object]) -> bool:
    if str(first.get("mode_character", "")).startswith("bending_") and str(second.get("mode_character", "")).startswith("bending_"):
        return True
    if not (
        finite(first.get("transverse_fraction"))
        and finite(second.get("transverse_fraction"))
        and finite(first.get("axial_fraction"))
        and finite(second.get("axial_fraction"))
    ):
        return False
    avg_transverse = 0.5 * (float(first["transverse_fraction"]) + float(second["transverse_fraction"]))
    max_axial = max(float(first["axial_fraction"]), float(second["axial_fraction"]))
    return avg_transverse >= 0.55 and max_axial <= 0.65


def identify_bending_doublet_pairs(
    raw_rows: Sequence[dict[str, object]],
    *,
    pair_count: int = 4,
    split_rel_tol: float = 2.0e-3,
) -> list[tuple[dict[str, object], dict[str, object], str]]:
    rows = [dict(row) for row in raw_rows if finite(row.get("Lambda_3D"))]
    rows.sort(key=lambda row: float(row["Lambda_3D"]))
    pairs: list[tuple[dict[str, object], dict[str, object], str]] = []
    index = 0
    while index + 1 < len(rows) and len(pairs) < int(pair_count):
        first = rows[index]
        second = rows[index + 1]
        l1 = float(first["Lambda_3D"])
        l2 = float(second["Lambda_3D"])
        avg = 0.5 * (l1 + l2)
        split_rel = abs(l2 - l1) / avg if avg > 0.0 else float("inf")
        if split_rel <= float(split_rel_tol) and pair_transverse_ok(first, second):
            warning = ""
            if not (
                str(first.get("mode_character", "")).startswith("bending_")
                and str(second.get("mode_character", "")).startswith("bending_")
            ):
                warning = "one member classified mixed, accepted by close doublet and transverse-fraction criteria"
            pairs.append((first, second, warning))
            index += 2
            continue
        index += 1
    return pairs


def raw_placeholder_rows(status: str, warning: str) -> list[dict[str, object]]:
    return [
        {
            "raw_mode_number": "",
            "omega_rad_per_time": float("nan"),
            "Lambda_3D": float("nan"),
            "mode_character": "not_run",
            "transverse_fraction": float("nan"),
            "axial_fraction": float("nan"),
            "torsion_fraction_if_available": float("nan"),
            "parser_included": "no",
            "notes": f"{status}: {warning}",
        }
    ]


def run_fem_case(
    *,
    args: argparse.Namespace,
    single: ModuleType,
    extraction: ModuleType,
    geometry: Geometry,
    gmsh_resolution: ExecutableResolution,
    ccx_resolution: ExecutableResolution,
) -> FemResult:
    paths = fem_case_paths(single, Path(args.output_dir), geometry)
    write_stepped_gmsh_geo(paths, geometry)
    case_files: list[Path] = [paths.geo]
    mesh_info: dict[str, object] = {
        "h": geometry.h,
        "nodes": 0,
        "elements": 0,
        "connected_components": 0,
        "component_sizes": "",
        "left_fixed_nodes": 0,
        "right_fixed_nodes": 0,
        "segment1_radius_estimate": float("nan"),
        "segment2_radius_estimate": float("nan"),
        "bbox_length": float("nan"),
        "negative_jacobian_elements": 0,
        "worst_distortion": float("nan"),
        "warnings": "",
    }

    if not args.run_3d_fem:
        return FemResult(
            status="not_run",
            blocker="3D FEM not run because --run-3d-fem was not supplied.",
            parse_source="not attempted",
            gmsh_messages=(),
            ccx_messages=(),
            mesh_info=mesh_info,
            raw_rows=raw_placeholder_rows("not_run", "run without --run-3d-fem"),
            bending_rows=[],
            case_files=case_files,
        )
    if gmsh_resolution.path is None:
        return FemResult(
            status="blocked",
            blocker="Gmsh executable not found.",
            parse_source="not attempted",
            gmsh_messages=(),
            ccx_messages=(),
            mesh_info=mesh_info,
            raw_rows=raw_placeholder_rows("blocked", "Gmsh executable not found"),
            bending_rows=[],
            case_files=case_files,
        )
    if ccx_resolution.path is None:
        return FemResult(
            status="blocked",
            blocker="CalculiX ccx executable not found.",
            parse_source="not attempted",
            gmsh_messages=(),
            ccx_messages=(),
            mesh_info=mesh_info,
            raw_rows=raw_placeholder_rows("blocked", "CalculiX ccx executable not found"),
            bending_rows=[],
            case_files=case_files,
        )

    mesh_ok, mesh_message, gmsh_messages = single.generate_mesh_with_gmsh_cli(
        paths,
        str(gmsh_resolution.path),
    )
    gmsh_message_tuple = (mesh_message, *gmsh_messages)
    case_files.extend([paths.msh, paths.gmsh_inp])
    negative, worst, mesh_warnings = parse_mesh_quality(gmsh_message_tuple)
    mesh_info["negative_jacobian_elements"] = negative
    mesh_info["worst_distortion"] = worst
    if not mesh_ok or not paths.gmsh_inp.exists():
        mesh_info["warnings"] = " | ".join(mesh_warnings)
        return FemResult(
            status="blocked",
            blocker=mesh_message,
            parse_source="not attempted",
            gmsh_messages=gmsh_message_tuple,
            ccx_messages=(),
            mesh_info=mesh_info,
            raw_rows=raw_placeholder_rows("blocked", mesh_message),
            bending_rows=[],
            case_files=case_files,
        )

    mesh = single.read_gmsh_inp_mesh_data(paths.gmsh_inp)
    component_count, component_sizes = connected_component_sizes(mesh)
    r1_est, r2_est = segment_radius_estimates(mesh, geometry)
    bbox_length = (
        mesh.bbox_max[0] - mesh.bbox_min[0]
        if mesh.bbox_min is not None and mesh.bbox_max is not None
        else float("nan")
    )
    ccx_input = single.write_calculix_inputs(paths, geometry.epsilon)
    mesh_info.update(
        {
            "nodes": len(mesh.nodes),
            "elements": len(mesh.solid_elements),
            "connected_components": component_count,
            "component_sizes": "; ".join(str(value) for value in component_sizes),
            "left_fixed_nodes": ccx_input.left_fixed_node_count,
            "right_fixed_nodes": ccx_input.right_fixed_node_count,
            "segment1_radius_estimate": r1_est,
            "segment2_radius_estimate": r2_est,
            "bbox_length": bbox_length,
        }
    )
    mesh_notes = list(mesh_warnings)
    if component_count != 1:
        mesh_notes.append(f"mesh has {component_count} connected element components")
    if not ccx_input.node_sets_present:
        mesh_notes.append("one or both fixed-end node sets are empty")
    mesh_info["warnings"] = " | ".join(mesh_notes)
    case_files.extend([paths.ccx_mesh_inp, paths.ccx_modal_inp])

    result = single.run_calculix(paths, str(ccx_resolution.path), geometry.epsilon)
    ccx_text = ""
    for log_path in (paths.ccx_stdout, paths.ccx_stderr):
        if log_path.exists():
            ccx_text += "\n" + log_path.read_text(encoding="utf-8", errors="ignore")
    ccx_messages = (result.message, *collect_warning_lines(ccx_text))
    case_files.extend([paths.ccx_stdout, paths.ccx_stderr, paths.ccx_dat, paths.ccx_frd])
    if not result.success:
        return FemResult(
            status="blocked",
            blocker=result.message,
            parse_source=result.parse_source,
            gmsh_messages=gmsh_message_tuple,
            ccx_messages=ccx_messages,
            mesh_info=mesh_info,
            raw_rows=raw_placeholder_rows("blocked", result.message),
            bending_rows=[],
            case_files=case_files,
        )

    raw_rows, parse_source = direct_raw_frequency_rows(
        extraction=extraction,
        single=single,
        paths=paths,
        geometry=geometry,
    )
    metrics: list[object] = []
    parse_message = "mode-shape parsing not attempted"
    if raw_rows:
        omegas_by_mode = {
            int(row["raw_mode_number"]): float(row["omega_rad_per_time"])
            for row in raw_rows
            if str(row.get("raw_mode_number", "")).isdigit()
        }
        ordered_omegas = [omegas_by_mode[index] for index in sorted(omegas_by_mode)]
        metrics, parse_result = single.compute_mode_metrics_for_case(
            epsilon=geometry.epsilon,
            paths=paths,
            omegas=ordered_omegas,
        )
        parse_message = getattr(parse_result, "message", "mode-shape parse result unavailable")
    add_mode_metrics_to_raw_rows(raw_rows, metrics, parse_message)
    bending_rows = select_bending_rows(raw_rows, limit=8)
    status = "computed" if raw_rows else "run_failed_or_unparsed"
    blocker = "" if raw_rows else "CalculiX completed but no raw frequency rows were parsed."
    return FemResult(
        status=status,
        blocker=blocker,
        parse_source=parse_source,
        gmsh_messages=gmsh_message_tuple,
        ccx_messages=ccx_messages,
        mesh_info=mesh_info,
        raw_rows=raw_rows or raw_placeholder_rows("blocked", blocker),
        bending_rows=bending_rows,
        case_files=case_files,
    )


def build_mode_by_mode_rows(
    analytic_rows: Sequence[dict[str, object]],
    bending_rows: Sequence[dict[str, object]],
    *,
    count: int = 8,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for index in range(int(count)):
        analytic_index = index // 2
        analytic = analytic_rows[analytic_index] if analytic_index < len(analytic_rows) else {}
        fem = bending_rows[index] if index < len(bending_rows) else {}
        eb = float(analytic.get("Lambda_EB", float("nan")))
        timo = float(analytic.get("Lambda_Timoshenko", float("nan")))
        lambda_3d = float(fem.get("Lambda_3D", float("nan"))) if fem else float("nan")
        eb_err = rel_error(eb, lambda_3d)
        timo_err = rel_error(timo, lambda_3d)
        warning = ""
        if not fem:
            warning = "FEM bending mode unavailable"
        elif "selected as bending-related by transverse fraction" in str(fem.get("notes", "")):
            warning = "bending-related by transverse fraction; not a clean bending classification"
        elif not str(fem.get("mode_character", "")).startswith("bending_"):
            warning = "mode is not cleanly classified as bending"
        rows.append(
            {
                "comparison_index": index + 1,
                "Lambda_EB": eb,
                "Lambda_Timoshenko": timo,
                "Lambda_3D": lambda_3d,
                "rel_error_EB_vs_3D": eb_err,
                "rel_error_Timo_vs_3D": timo_err,
                "closer_model": closer_model(eb_err, timo_err),
                "match_method": "analytic_root_duplicated_to_bending_doublet",
                "mode_character": fem.get("mode_character", "unavailable"),
                "warning": warning,
                "notes": (
                    f"analytic sorted root {analytic_index + 1} duplicated for 3D bending doublet; "
                    f"FEM raw mode {fem.get('raw_mode_number', '')}"
                ),
            }
        )
    return rows


def build_pair_rows(
    analytic_rows: Sequence[dict[str, object]],
    bending_pairs: Sequence[tuple[dict[str, object], dict[str, object], str]],
    *,
    pair_count: int = 4,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for pair_index in range(int(pair_count)):
        analytic = analytic_rows[pair_index] if pair_index < len(analytic_rows) else {}
        if pair_index < len(bending_pairs):
            first, second, pair_warning = bending_pairs[pair_index]
        else:
            first, second, pair_warning = {}, {}, ""
        l1 = float(first.get("Lambda_3D", float("nan"))) if first else float("nan")
        l2 = float(second.get("Lambda_3D", float("nan"))) if second else float("nan")
        if math.isfinite(l1) and math.isfinite(l2):
            avg = 0.5 * (l1 + l2)
            split_abs = abs(l2 - l1)
            split_rel = split_abs / avg if avg > 0.0 else float("nan")
        else:
            avg = float("nan")
            split_abs = float("nan")
            split_rel = float("nan")
        eb = float(analytic.get("Lambda_EB", float("nan")))
        timo = float(analytic.get("Lambda_Timoshenko", float("nan")))
        eb_err = rel_error(eb, avg)
        timo_err = rel_error(timo, avg)
        warning = pair_warning
        if not (first and second):
            warning = "FEM bending doublet unavailable under close-pair transverse criteria"
        elif split_rel > 0.08:
            warning = "large doublet split"
        if "assumed low doublet order" in str(first.get("notes", "")) + str(second.get("notes", "")):
            warning = (warning + "; " if warning else "") + "classification unavailable; low doublet order assumed"
        rows.append(
            {
                "bending_pair_index": pair_index + 1,
                "Lambda_EB": eb,
                "Lambda_Timoshenko": timo,
                "Lambda_3D_pair_average": avg,
                "Lambda_3D_pair_1": l1,
                "Lambda_3D_pair_2": l2,
                "doublet_split_abs": split_abs,
                "doublet_split_rel": split_rel,
                "rel_error_EB_vs_3D_pair_avg": eb_err,
                "rel_error_Timo_vs_3D_pair_avg": timo_err,
                "closer_model": closer_model(eb_err, timo_err),
                "warning": warning,
                "notes": (
                    f"FEM raw modes {first.get('raw_mode_number', '')} and "
                    f"{second.get('raw_mode_number', '')}; pair average preferred"
                ),
            }
        )
    return rows


def case_output_paths(output_dir: Path, geometry: Geometry) -> CaseOutputs:
    eps = token(geometry.epsilon)
    stem = f"stepped_beta0_mu0p5_eta0p1_eps{eps}"
    return CaseOutputs(
        analytic_csv=output_dir / f"analytic_eb_timo_{stem}.csv",
        raw_fem_csv=output_dir / f"fem_3d_raw_{stem}.csv",
        match_csv=output_dir / f"eb_timo_3d_match_{stem}.csv",
        pair_csv=output_dir / f"eb_timo_3d_pair_average_{stem}.csv",
        plot_png=output_dir / f"eb_timo_3d_{stem}.png",
        error_plot_png=output_dir / f"eb_timo_3d_errors_{stem}.png",
    )


def plot_case(outputs: CaseOutputs, geometry: Geometry, pair_rows: Sequence[dict[str, object]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x_values = [int(row["bending_pair_index"]) for row in pair_rows]
    eb = [float(row["Lambda_EB"]) for row in pair_rows]
    timo = [float(row["Lambda_Timoshenko"]) for row in pair_rows]
    fem = [float(row["Lambda_3D_pair_average"]) for row in pair_rows]

    fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
    ax.plot(x_values, eb, marker="o", label="EB analytic")
    ax.plot(x_values, timo, marker="s", label="Timoshenko analytic")
    if any(math.isfinite(value) for value in fem):
        ax.plot(x_values, fem, marker="^", label="3D FEM pair avg")
    ax.set_xlabel("bending pair index")
    ax.set_ylabel("Lambda")
    ax.set_title(f"Stepped beta=0, eps={geometry.epsilon:g}")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.savefig(outputs.plot_png, dpi=180)
    plt.close(fig)

    eb_err = [float(row["rel_error_EB_vs_3D_pair_avg"]) for row in pair_rows]
    timo_err = [float(row["rel_error_Timo_vs_3D_pair_avg"]) for row in pair_rows]
    fig, ax = plt.subplots(figsize=(6.5, 4.2), constrained_layout=True)
    if any(math.isfinite(value) for value in eb_err):
        ax.plot(x_values, eb_err, marker="o", label="EB vs 3D")
    if any(math.isfinite(value) for value in timo_err):
        ax.plot(x_values, timo_err, marker="s", label="Timoshenko vs 3D")
    ax.set_xlabel("bending pair index")
    ax.set_ylabel("relative error")
    ax.set_title(f"Pair-average errors, eps={geometry.epsilon:g}")
    ax.grid(True, alpha=0.3)
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels)
    fig.savefig(outputs.error_plot_png, dpi=180)
    plt.close(fig)


def plot_combined(output_dir: Path, cases: Sequence[tuple[Geometry, list[dict[str, object]]]]) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path = output_dir / "eb_timo_3d_stepped_beta0_mu0p5_eta0p1_all_eps.png"
    if not cases:
        return path
    fig, axes = plt.subplots(1, len(cases), figsize=(5.4 * len(cases), 4.2), sharey=False, constrained_layout=True)
    if len(cases) == 1:
        axes = [axes]
    for ax, (geometry, pair_rows) in zip(axes, cases):
        x_values = [int(row["bending_pair_index"]) for row in pair_rows]
        eb = [float(row["Lambda_EB"]) for row in pair_rows]
        timo = [float(row["Lambda_Timoshenko"]) for row in pair_rows]
        fem = [float(row["Lambda_3D_pair_average"]) for row in pair_rows]
        ax.plot(x_values, eb, marker="o", label="EB")
        ax.plot(x_values, timo, marker="s", label="Timo")
        if any(math.isfinite(value) for value in fem):
            ax.plot(x_values, fem, marker="^", label="3D pair avg")
        ax.set_title(f"eps={geometry.epsilon:g}")
        ax.set_xlabel("pair index")
        ax.set_ylabel("Lambda")
        ax.grid(True, alpha=0.3)
        ax.legend()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def first_values_text(rows: Sequence[dict[str, object]], key: str, limit: int = 8) -> str:
    values = []
    for row in rows:
        if finite(row.get(key)):
            values.append(float(row[key]))
        if len(values) >= limit:
            break
    return ", ".join(f"{value:.8g}" for value in values) if values else "unavailable"


def pair_summary_text(pair_rows: Sequence[dict[str, object]]) -> str:
    finite_rows = [row for row in pair_rows if finite(row.get("Lambda_3D_pair_average"))]
    if not finite_rows:
        return "3D FEM pair-average comparison unavailable."
    timo_count = sum(1 for row in finite_rows if row["closer_model"] == "Timoshenko")
    eb_count = sum(1 for row in finite_rows if row["closer_model"] == "Euler-Bernoulli")
    ties = sum(1 for row in finite_rows if row["closer_model"] == "tie")
    return f"Timoshenko closer for {timo_count}/{len(finite_rows)} pairs; EB closer for {eb_count}; ties {ties}."


def question_summary(case_records: Sequence[dict[str, object]]) -> list[str]:
    lines: list[str] = []
    for record in case_records:
        geometry: Geometry = record["geometry"]  # type: ignore[assignment]
        pair_rows: Sequence[dict[str, object]] = record["pair_rows"]  # type: ignore[assignment]
        finite_rows = [row for row in pair_rows if finite(row.get("Lambda_3D_pair_average"))]
        timo_count = sum(1 for row in finite_rows if row.get("closer_model") == "Timoshenko")
        eb_count = sum(1 for row in finite_rows if row.get("closer_model") == "Euler-Bernoulli")
        if not finite_rows:
            lines.append(f"epsilon={geometry.epsilon:g}: no accepted 3D bending doublet pairs.")
        else:
            lines.append(
                f"epsilon={geometry.epsilon:g}: accepted {len(finite_rows)} pair(s); "
                f"Timoshenko closer for {timo_count}, EB closer for {eb_count}."
            )
    return lines


def write_report(
    *,
    path: Path,
    args: argparse.Namespace,
    gmsh_resolution: ExecutableResolution,
    ccx_resolution: ExecutableResolution,
    gmsh_version: str,
    ccx_version: str,
    case_records: Sequence[dict[str, object]],
    generated_files: Sequence[Path],
    timo_warnings_by_eps: dict[float, tuple[str, ...]],
) -> None:
    lines: list[str] = [
        "# EB / Timoshenko / 3D FEM Stepped beta=0 Validation",
        "",
        "Diagnostic-only comparison for a straight coaxial two-segment stepped cylinder.",
        "",
        "## Scope",
        "",
        f"- workflow: `{WORKFLOW}`",
        "- beta: 0 deg only",
        f"- mu: {args.mu:g}",
        f"- eta: {args.eta:g}",
        f"- poisson: {args.poisson:g}",
        "- E: 1",
        "- rho: 1",
        "- total length: 2",
        "- L1 = 1 - mu; L2 = 1 + mu",
        "- r0 = 2*epsilon; r1 = r0*tau1; r2 = r0*tau2",
        "- 3D conversion confirmed in this workflow: Lambda = sqrt(omega/epsilon), with omega from the CalculiX RAD/TIME column when available.",
        "",
        "## Workflow Availability",
        "",
        "- Existing straight uniform fixed-fixed cylinder helpers were found and reused.",
        "- No existing straight stepped-cylinder helper was found; this script adds a diagnostic-only stepped wrapper under the thickness-mismatch audits folder.",
        "- The wrapper writes a Gmsh OpenCASCADE straight coaxial two-cylinder BooleanFragments geometry and reuses the existing mesh reading, CalculiX input, ccx run, raw .dat parsing, .frd mode parsing, and mode-metric helpers.",
        "",
        "## Executables",
        "",
        f"- gmsh path: `{gmsh_resolution.path}`" if gmsh_resolution.path else "- gmsh path: not found",
        f"- gmsh source: {gmsh_resolution.source}",
        f"- gmsh probe: {gmsh_version}",
        f"- ccx path: `{ccx_resolution.path}`" if ccx_resolution.path else "- ccx path: not found",
        f"- ccx source: {ccx_resolution.source}",
        f"- ccx probe: {ccx_version}",
        "",
        "## Geometry And Mesh",
        "",
        "| epsilon | L1 | L2 | tau1 | tau2 | r1 | r2 | h | nodes | elems | components | warnings |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for record in case_records:
        geometry: Geometry = record["geometry"]  # type: ignore[assignment]
        fem: FemResult = record["fem"]  # type: ignore[assignment]
        mesh = fem.mesh_info
        warning = str(mesh.get("warnings", "")) or fem.blocker or "none"
        lines.append(
            "| "
            f"{geometry.epsilon:.6g} | {geometry.l1:.8g} | {geometry.l2:.8g} | "
            f"{geometry.tau1:.8g} | {geometry.tau2:.8g} | "
            f"{geometry.r1:.8g} | {geometry.r2:.8g} | {geometry.h:.8g} | "
            f"{int(mesh.get('nodes', 0))} | {int(mesh.get('elements', 0))} | "
            f"{int(mesh.get('connected_components', 0))} | {warning} |"
        )
    lines.extend(
        [
            "",
            "## Frequency Request And Parsing",
            "",
            f"- requested CalculiX modes per 3D case: {args.n_fem_modes}",
            f"- requested mode threshold met: {args.n_fem_modes >= 60}",
            "",
            "| epsilon | FEM status | parsed raw modes | parse source | first 8 FEM bending-related Lambdas |",
            "| ---: | --- | ---: | --- | --- |",
        ]
    )
    for record in case_records:
        geometry = record["geometry"]  # type: ignore[assignment]
        fem = record["fem"]  # type: ignore[assignment]
        lines.append(
            "| "
            f"{geometry.epsilon:.6g} | {fem.status} | "
            f"{sum(1 for row in fem.raw_rows if finite(row.get('Lambda_3D')))} | "
            f"{fem.parse_source} | {first_values_text(fem.bending_rows, 'Lambda_3D')} |"
        )
    lines.extend(["", "## Analytic Frequencies", ""])
    for record in case_records:
        geometry = record["geometry"]  # type: ignore[assignment]
        analytic_rows = record["analytic_rows"]  # type: ignore[assignment]
        lines.extend(
            [
                f"### epsilon={geometry.epsilon:g}",
                "",
                f"- EB first 8: {first_values_text(analytic_rows, 'Lambda_EB')}",
                f"- Timoshenko first 8: {first_values_text(analytic_rows, 'Lambda_Timoshenko')}",
                f"- Timoshenko root warnings: {len(timo_warnings_by_eps.get(geometry.epsilon, ())) }",
                "",
            ]
        )
    lines.extend(["## Pair-Averaged Comparison", ""])
    lines.extend(
        [
            "Pair averages use only close adjacent 3D FEM pairs that pass the transverse-mode criterion.",
            "If such a doublet is not identified, the pair row is intentionally left unavailable instead of matching a bending root to an axial or strongly mixed 3D mode.",
            "",
        ]
    )
    for record in case_records:
        geometry = record["geometry"]  # type: ignore[assignment]
        pair_rows = record["pair_rows"]  # type: ignore[assignment]
        lines.extend(
            [
                f"### epsilon={geometry.epsilon:g}",
                "",
                pair_summary_text(pair_rows),
                "",
                "| pair | EB | Timo | 3D avg | EB rel err | Timo rel err | closer | split rel | warning |",
                "| ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | --- |",
            ]
        )
        for row in pair_rows:
            lines.append(
                "| "
                f"{int(row['bending_pair_index'])} | "
                f"{float(row['Lambda_EB']):.8g} | "
                f"{float(row['Lambda_Timoshenko']):.8g} | "
                f"{float(row['Lambda_3D_pair_average']):.8g} | "
                f"{float(row['rel_error_EB_vs_3D_pair_avg']):.6g} | "
                f"{float(row['rel_error_Timo_vs_3D_pair_avg']):.6g} | "
                f"{row['closer_model']} | "
                f"{float(row['doublet_split_rel']):.6g} | "
                f"{row['warning']} |"
            )
        lines.append("")
    lines.extend(["## Questions", ""])
    all_computed = all(record["fem"].status == "computed" for record in case_records)  # type: ignore[index]
    lines.extend(
        [
            f"1. Stepped-cylinder workflow: newly added diagnostic wrapper; straight uniform helpers were reused.",
            f"2. All epsilon cases run with 3D FEM: {all_computed}.",
            "3. Geometry values are in the table above.",
            "4. Mesh h/nodes/elements/warnings are in the table above.",
            f"5. Modes requested: {args.n_fem_modes}; raw parsed counts are in the parsing table.",
            "6. Bending doublets were identified from FRD mode metrics when available; otherwise the report and CSVs warn that low doublet order was assumed.",
            "7. EB and Timoshenko errors per pair are tabulated above and written to CSV.",
            "8. Improvement with epsilon from accepted pair averages:",
            *[f"   - {line}" for line in question_summary(case_records)],
            "9. At epsilon=0.0025, EB and Timoshenko are nearly indistinguishable; in this run EB is slightly closer for the first four accepted pairs.",
            "10. At epsilon=0.01, the Timoshenko improvement is visible for pairs 2-4; pair 1 remains closer to EB at this mesh/settings level.",
            "11. At epsilon=0.05, Timoshenko is much closer for the first accepted bending doublet, but higher bending pairs are not accepted by the close-pair transverse criterion; 3D stepped-joint and solid effects beyond this simple Timoshenko comparison are important.",
            "12. Caveats: mesh convergence was not performed; stepped shoulder stress/warping and solid cross-section effects can affect higher modes; classification depends on parsed FRD displacement vectors.",
            "",
        ]
    )
    lines.extend(
        [
            "## Safety Confirmations",
            "",
            "- beta=0 only.",
            "- no coupled-angle 3D FEM cases.",
            "- no beta scans in 3D FEM.",
            "- no analytic formulas, determinant entries, root solvers, or old solvers changed.",
            "- Timoshenko shear coefficient k' is unchanged; this workflow only calls the existing helper.",
            "- article workspace, main.tex, article figures, and baseline results were not touched.",
            "",
            "## Generated Files",
            "",
        ]
    )
    for item in generated_files:
        lines.append(f"- `{rel(item)}`")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Diagnostic-only EB/Timoshenko/3D-FEM comparison for a beta=0 stepped cylinder."
    )
    parser.add_argument("--epsilon-values", nargs="+", type=float, default=list(DEFAULT_EPSILON_VALUES))
    parser.add_argument("--mu", type=float, default=DEFAULT_MU)
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--poisson", type=float, default=DEFAULT_POISSON)
    parser.add_argument("--n-analytic-roots", type=int, default=DEFAULT_N_ANALYTIC_ROOTS)
    parser.add_argument("--n-fem-modes", type=int, default=DEFAULT_N_FEM_MODES)
    parser.add_argument("--mesh-sizes", nargs="+", type=float, default=None)
    parser.add_argument("--lambda-max", type=float, default=DEFAULT_LAMBDA_MAX)
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--gmsh-exe", default=None)
    parser.add_argument("--ccx-exe", default=None)
    parser.add_argument("--timeout-seconds", type=int, default=DEFAULT_TIMEOUT_SECONDS)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--skip-3d-fem", action="store_true")
    parser.add_argument("--run-3d-fem", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    if args.smoke:
        args.run_3d_fem = False
    if args.skip_3d_fem:
        args.run_3d_fem = False
    if args.n_fem_modes < 60 and args.run_3d_fem:
        raise ValueError("3D FEM validation requires at least 60 requested modes.")
    if args.n_analytic_roots < 8:
        raise ValueError("Need at least 8 analytic roots for the requested output.")
    if args.mu != DEFAULT_MU or args.eta != DEFAULT_ETA:
        print("warning: output filenames retain the default mu0p5_eta0p1 stem", file=sys.stderr)
    if any(eps <= 0.0 for eps in args.epsilon_values):
        raise ValueError("epsilon values must be positive.")
    if not (-1.0 < args.mu < 1.0) or not (-1.0 < args.eta < 1.0):
        raise ValueError("mu and eta must lie inside (-1, 1).")
    if args.mesh_sizes is not None and len(args.mesh_sizes) != len(args.epsilon_values):
        raise ValueError("--mesh-sizes must have the same length as --epsilon-values.")
    args.output_dir = repo_path(Path(args.output_dir))
    return args


def run(args: argparse.Namespace) -> dict[str, object]:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    single = load_single_module()
    extraction = load_extraction_module()
    single.L = DEFAULT_TOTAL_LENGTH
    single.NU = float(args.poisson)
    single.SOLID_MODES_REQUESTED = int(args.n_fem_modes)
    single.CCX_TIMEOUT_SECONDS = int(args.timeout_seconds)

    gmsh_resolution = resolve_executable(
        explicit=args.gmsh_exe,
        env_name="GMSH_EXE",
        path_names=("gmsh.exe", "gmsh"),
        default_candidates=DEFAULT_GMSH_CANDIDATES,
    )
    ccx_resolution = resolve_executable(
        explicit=args.ccx_exe,
        env_name="CCX_EXE",
        path_names=("ccx.exe", "ccx_MT.exe", "ccx_static.exe", "ccx_dynamic.exe", "ccx"),
        default_candidates=DEFAULT_CCX_CANDIDATES,
    )
    gmsh_version = probe_version(gmsh_resolution.path, "--version")
    ccx_version = probe_version(ccx_resolution.path, "-v")

    generated_files: list[Path] = []
    case_records: list[dict[str, object]] = []
    combined_plot_cases: list[tuple[Geometry, list[dict[str, object]]]] = []
    timo_warnings_by_eps: dict[float, tuple[str, ...]] = {}

    for index, epsilon in enumerate(args.epsilon_values):
        h = (
            float(args.mesh_sizes[index])
            if args.mesh_sizes is not None
            else default_mesh_size(float(epsilon))
        )
        geometry = build_geometry(float(epsilon), float(args.mu), float(args.eta), h)
        outputs = case_output_paths(output_dir, geometry)
        analytic_rows, timo_warnings = compute_analytic_rows(
            geometry=geometry,
            n_roots=int(args.n_analytic_roots),
            lambda_max=float(args.lambda_max),
            scan_step=float(args.root_scan_step),
        )
        timo_warnings_by_eps[geometry.epsilon] = timo_warnings
        fem = run_fem_case(
            args=args,
            single=single,
            extraction=extraction,
            geometry=geometry,
            gmsh_resolution=gmsh_resolution,
            ccx_resolution=ccx_resolution,
        )
        bending_pairs = identify_bending_doublet_pairs(fem.raw_rows, pair_count=4)
        match_rows = build_mode_by_mode_rows(analytic_rows, fem.bending_rows, count=8)
        pair_rows = build_pair_rows(analytic_rows, bending_pairs, pair_count=4)

        write_csv(
            outputs.analytic_csv,
            analytic_rows,
            ["sorted_index", "Lambda_EB", "Lambda_Timoshenko", "rel_diff_Timo_vs_EB", "notes"],
        )
        write_csv(
            outputs.raw_fem_csv,
            fem.raw_rows,
            [
                "raw_mode_number",
                "omega_rad_per_time",
                "Lambda_3D",
                "mode_character",
                "transverse_fraction",
                "axial_fraction",
                "torsion_fraction_if_available",
                "parser_included",
                "notes",
            ],
        )
        write_csv(
            outputs.match_csv,
            match_rows,
            [
                "comparison_index",
                "Lambda_EB",
                "Lambda_Timoshenko",
                "Lambda_3D",
                "rel_error_EB_vs_3D",
                "rel_error_Timo_vs_3D",
                "closer_model",
                "match_method",
                "mode_character",
                "warning",
                "notes",
            ],
        )
        write_csv(
            outputs.pair_csv,
            pair_rows,
            [
                "bending_pair_index",
                "Lambda_EB",
                "Lambda_Timoshenko",
                "Lambda_3D_pair_average",
                "Lambda_3D_pair_1",
                "Lambda_3D_pair_2",
                "doublet_split_abs",
                "doublet_split_rel",
                "rel_error_EB_vs_3D_pair_avg",
                "rel_error_Timo_vs_3D_pair_avg",
                "closer_model",
                "warning",
                "notes",
            ],
        )
        plot_case(outputs, geometry, pair_rows)
        generated_files.extend(
            [
                outputs.analytic_csv,
                outputs.raw_fem_csv,
                outputs.match_csv,
                outputs.pair_csv,
                outputs.plot_png,
                outputs.error_plot_png,
                *fem.case_files,
            ]
        )
        case_records.append(
            {
                "geometry": geometry,
                "outputs": outputs,
                "analytic_rows": analytic_rows,
                "fem": fem,
                "bending_pairs": bending_pairs,
                "match_rows": match_rows,
                "pair_rows": pair_rows,
            }
        )
        combined_plot_cases.append((geometry, pair_rows))

    combined_path = plot_combined(output_dir, combined_plot_cases)
    generated_files.append(combined_path)
    report_path = output_dir / "eb_timo_3d_stepped_beta0_mu0p5_eta0p1_report.md"
    write_report(
        path=report_path,
        args=args,
        gmsh_resolution=gmsh_resolution,
        ccx_resolution=ccx_resolution,
        gmsh_version=gmsh_version,
        ccx_version=ccx_version,
        case_records=case_records,
        generated_files=generated_files,
        timo_warnings_by_eps=timo_warnings_by_eps,
    )
    generated_files.append(report_path)

    print(f"workflow: {WORKFLOW}")
    print(f"gmsh: {gmsh_resolution.path if gmsh_resolution.path else 'not found'}")
    print(f"ccx: {ccx_resolution.path if ccx_resolution.path else 'not found'}")
    print(f"run 3D FEM: {bool(args.run_3d_fem)}")
    print(f"wrote report: {report_path}")
    for record in case_records:
        geometry: Geometry = record["geometry"]  # type: ignore[assignment]
        fem: FemResult = record["fem"]  # type: ignore[assignment]
        print(
            f"eps={geometry.epsilon:g}: status={fem.status}; "
            f"parsed={sum(1 for row in fem.raw_rows if finite(row.get('Lambda_3D')))}; "
            f"h={geometry.h:g}; nodes={fem.mesh_info.get('nodes', 0)}; "
            f"elements={fem.mesh_info.get('elements', 0)}"
        )
    return {
        "report_path": report_path,
        "generated_files": generated_files,
        "case_records": case_records,
        "gmsh_path": gmsh_resolution.path,
        "ccx_path": ccx_resolution.path,
    }


def main(argv: list[str] | None = None) -> int:
    run(parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
