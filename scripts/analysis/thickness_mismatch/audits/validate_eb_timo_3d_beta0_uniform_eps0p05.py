from __future__ import annotations

import argparse
import csv
import importlib.util
import math
from dataclasses import dataclass
from pathlib import Path
import sys
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


DEFAULT_EPSILON = 0.05
OPTIONAL_STRESS_EPSILON = 0.1
DEFAULT_POISSON = 0.3
DEFAULT_TOTAL_LENGTH = 2.0
DEFAULT_MESH_SIZE = 0.025
DEFAULT_N_ANALYTIC_ROOTS = 8
DEFAULT_N_FEM_MODES = 60
DEFAULT_BENDING_PAIRS = 4
DEFAULT_LAMBDA_MAX = 60.0
DEFAULT_ROOT_SCAN_STEP = 0.01
DEFAULT_TIMEOUT_SECONDS = 1200
DEFAULT_OUTPUT_DIR = (
    Path("results")
    / "eb_vs_timoshenko_3d_validation"
    / "uniform_beta0_eps0p05"
)
WORKFLOW = "diagnostic_straight_uniform_fixed_fixed_cylinder_beta0"


def _load_module(module_name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


STEPPED = _load_module(
    "validate_eb_timo_3d_beta0_stepped_shared_for_uniform",
    REPO_ROOT
    / "scripts"
    / "analysis"
    / "thickness_mismatch"
    / "audits"
    / "validate_eb_timo_3d_beta0_stepped.py",
)
FORMULAS = STEPPED.FORMULAS
TIMO = STEPPED.TIMO
SINGLE_REF = _load_module(
    "compare_single_rod_eb_timoshenko_for_uniform_beta0_eps0p05",
    REPO_ROOT / "scripts" / "analysis" / "compare_single_rod_eb_timoshenko.py",
)


@dataclass(frozen=True)
class Geometry:
    epsilon: float
    h: float
    total_length: float
    radius: float
    poisson: float
    tau1: float = 1.0
    tau2: float = 1.0
    l1: float = 1.0
    l2: float = 1.0
    mu: float = 0.0
    eta: float = 0.0


@dataclass(frozen=True)
class OutputPaths:
    analytic_csv: Path
    raw_fem_csv: Path
    first8_csv: Path
    pair_csv: Path
    first8_png: Path
    first8_error_png: Path
    pair_png: Path
    pair_error_png: Path


@dataclass
class FemResult:
    status: str
    blocker: str
    parse_source: str
    gmsh_messages: tuple[str, ...]
    ccx_messages: tuple[str, ...]
    mesh_info: dict[str, object]
    raw_rows: list[dict[str, object]]
    case_files: list[Path]


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def token(value: float) -> str:
    return f"{float(value):.6g}".replace("-", "m").replace(".", "p")


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def fmt(value: object) -> object:
    if isinstance(value, float):
        if not math.isfinite(value):
            return "nan"
        return f"{value:.16e}"
    return value


def finite(value: object) -> bool:
    return STEPPED.finite(value)


def write_csv(path: Path, rows: Sequence[dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key, "")) for key in fieldnames})


def case_stem(epsilon: float) -> str:
    return f"uniform_beta0_mu0_eta0_eps{token(epsilon)}"


def output_paths(output_dir: Path, epsilon: float) -> OutputPaths:
    stem = case_stem(epsilon)
    return OutputPaths(
        analytic_csv=output_dir / f"analytic_eb_timo_{stem}.csv",
        raw_fem_csv=output_dir / f"fem_3d_raw_{stem}.csv",
        first8_csv=output_dir / f"eb_timo_3d_match_{stem}_first8.csv",
        pair_csv=output_dir / f"eb_timo_3d_pair_average_{stem}.csv",
        first8_png=output_dir / f"eb_timo_3d_{stem}_first8.png",
        first8_error_png=output_dir / f"eb_timo_3d_{stem}_first8_errors.png",
        pair_png=output_dir / f"eb_timo_3d_{stem}_bending_pair_average.png",
        pair_error_png=output_dir / f"eb_timo_3d_{stem}_bending_pair_average_errors.png",
    )


def fem_case_paths(single: ModuleType, output_dir: Path, geometry: Geometry) -> object:
    stem = case_stem(geometry.epsilon)
    case_dir = output_dir / "solid_fem_cases" / f"eps{token(geometry.epsilon)}"
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


def axial_roots(epsilon: float, count: int) -> list[float]:
    return [math.sqrt(index * math.pi / (2.0 * float(epsilon))) for index in range(1, int(count) + 1)]


def bending_reference_rows(epsilon: float, poisson: float, count: int) -> list[dict[str, object]]:
    SINGLE_REF.L = DEFAULT_TOTAL_LENGTH
    SINGLE_REF.NU = float(poisson)
    section = SINGLE_REF.section_from_epsilon(float(epsilon))
    alphas = SINGLE_REF.fixed_fixed_eb_roots(int(count))
    eb_omegas = [SINGLE_REF.omega_eb(float(alpha), section) for alpha in alphas]
    eb_lambdas = [SINGLE_REF.lambda_from_omega(float(omega), section) for omega in eb_omegas]
    timo_omegas = SINGLE_REF.find_timoshenko_roots(section, eb_omegas, int(count))
    timo_lambdas = [
        SINGLE_REF.lambda_from_omega(float(omega), section)
        for omega in timo_omegas
    ]
    rows: list[dict[str, object]] = []
    for index in range(int(count)):
        rows.append(
            {
                "bending_pair_index": index + 1,
                "alpha_EB": float(alphas[index]),
                "Lambda_EB": float(eb_lambdas[index]),
                "Lambda_Timoshenko": float(timo_lambdas[index]),
                "Omega_Timoshenko": float(timo_omegas[index]),
                "analytic_mode_note": "single straight fixed-fixed bending root, total length 2",
            }
        )
    return rows


def classify_analytic_root(lambda_eb: float, bending_rows: Sequence[dict[str, object]], epsilon: float) -> str:
    candidates: list[tuple[float, str]] = []
    for row in bending_rows:
        value = float(row["Lambda_EB"])
        candidates.append((abs(lambda_eb - value), f"bending-like EB root {int(row['bending_pair_index'])}"))
    for index, value in enumerate(axial_roots(float(epsilon), 8), start=1):
        candidates.append((abs(lambda_eb - value), f"axial-like EB root {index}"))
    if not candidates or not math.isfinite(float(lambda_eb)):
        return "sorted in-plane root; unclassified"
    distance, label = min(candidates, key=lambda item: item[0])
    tolerance = 2.0e-6 * max(1.0, abs(float(lambda_eb)))
    if distance <= tolerance:
        return f"sorted in-plane root; {label}"
    return "sorted in-plane root; no clean bending/axial label"


def analytic_rows(
    *,
    geometry: Geometry,
    n_roots: int,
    lambda_max: float,
    scan_step: float,
    bending_rows: Sequence[dict[str, object]],
) -> tuple[list[dict[str, object]], tuple[str, ...]]:
    eb_roots = FORMULAS.find_first_n_roots_eta(
        beta=0.0,
        mu=0.0,
        epsilon=geometry.epsilon,
        eta=0.0,
        n_roots=int(n_roots),
        Lmax0=float(lambda_max),
        scan_step=float(scan_step),
    )
    timo_roots, timo_warnings = TIMO.timo_sorted_roots(
        beta_deg=0.0,
        mu=0.0,
        epsilon=geometry.epsilon,
        eta=0.0,
        n_roots=int(n_roots),
    )
    rows: list[dict[str, object]] = []
    for index in range(int(n_roots)):
        eb = float(eb_roots[index]) if index < len(eb_roots) else float("nan")
        timo = float(timo_roots[index]) if index < len(timo_roots) else float("nan")
        rel_diff = (timo - eb) / eb if math.isfinite(eb) and eb != 0.0 and math.isfinite(timo) else float("nan")
        note = classify_analytic_root(eb, bending_rows, geometry.epsilon)
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
                "analytic_mode_note": note,
            }
        )
    return rows, tuple(timo_warnings)


def configure_single_module(single: ModuleType, args: argparse.Namespace, geometry: Geometry) -> None:
    single.L = float(geometry.total_length)
    single.NU = float(args.poisson)
    single.SOLID_MODES_REQUESTED = int(args.n_fem_modes)
    single.CCX_TIMEOUT_SECONDS = int(args.timeout_seconds)
    single.mesh_size_for_epsilon = lambda _epsilon: float(geometry.h)


def run_fem_case(
    *,
    args: argparse.Namespace,
    single: ModuleType,
    extraction: ModuleType,
    geometry: Geometry,
    gmsh_resolution: object,
    ccx_resolution: object,
) -> FemResult:
    paths = fem_case_paths(single, Path(args.output_dir), geometry)
    configure_single_module(single, args, geometry)
    single.write_gmsh_geo(float(geometry.epsilon), paths)
    case_files: list[Path] = [paths.geo]
    mesh_info: dict[str, object] = {
        "h": geometry.h,
        "nodes": 0,
        "elements": 0,
        "connected_components": 0,
        "component_sizes": "",
        "left_fixed_nodes": 0,
        "right_fixed_nodes": 0,
        "radius_estimate": float("nan"),
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
            raw_rows=STEPPED.raw_placeholder_rows("not_run", "run without --run-3d-fem"),
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
            raw_rows=STEPPED.raw_placeholder_rows("blocked", "Gmsh executable not found"),
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
            raw_rows=STEPPED.raw_placeholder_rows("blocked", "CalculiX ccx executable not found"),
            case_files=case_files,
        )

    mesh_ok, mesh_message, gmsh_messages = single.generate_mesh_with_gmsh_cli(paths, str(gmsh_resolution.path))
    gmsh_message_tuple = (mesh_message, *gmsh_messages)
    case_files.extend([paths.msh, paths.gmsh_inp])
    negative, worst, mesh_warnings = STEPPED.parse_mesh_quality(gmsh_message_tuple)
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
            raw_rows=STEPPED.raw_placeholder_rows("blocked", mesh_message),
            case_files=case_files,
        )

    mesh = single.read_gmsh_inp_mesh_data(paths.gmsh_inp)
    component_count, component_sizes = STEPPED.connected_component_sizes(mesh)
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
            "radius_estimate": float(mesh.radial_max),
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
    ccx_messages = (result.message, *STEPPED.collect_warning_lines(ccx_text))
    case_files.extend([paths.ccx_stdout, paths.ccx_stderr, paths.ccx_dat, paths.ccx_frd])
    if not result.success:
        return FemResult(
            status="blocked",
            blocker=result.message,
            parse_source=result.parse_source,
            gmsh_messages=gmsh_message_tuple,
            ccx_messages=ccx_messages,
            mesh_info=mesh_info,
            raw_rows=STEPPED.raw_placeholder_rows("blocked", result.message),
            case_files=case_files,
        )

    raw_rows, parse_source = STEPPED.direct_raw_frequency_rows(
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
    STEPPED.add_mode_metrics_to_raw_rows(raw_rows, metrics, parse_message)
    return FemResult(
        status="computed" if raw_rows else "run_failed_or_unparsed",
        blocker="" if raw_rows else "CalculiX completed but no raw frequency rows were parsed.",
        parse_source=parse_source,
        gmsh_messages=gmsh_message_tuple,
        ccx_messages=ccx_messages,
        mesh_info=mesh_info,
        raw_rows=raw_rows or STEPPED.raw_placeholder_rows("blocked", "no parsed FEM frequencies"),
        case_files=case_files,
    )


def first8_rows(
    analytic: Sequence[dict[str, object]],
    raw_fem: Sequence[dict[str, object]],
    count: int = 8,
) -> list[dict[str, object]]:
    fem_rows = [row for row in raw_fem if finite(row.get("Lambda_3D"))]
    fem_rows.sort(key=lambda row: float(row["Lambda_3D"]))
    rows: list[dict[str, object]] = []
    for index in range(int(count)):
        analytic_row = analytic[index] if index < len(analytic) else {}
        fem = fem_rows[index] if index < len(fem_rows) else {}
        eb = float(analytic_row.get("Lambda_EB", float("nan")))
        timo = float(analytic_row.get("Lambda_Timoshenko", float("nan")))
        lambda_3d = float(fem.get("Lambda_3D", float("nan"))) if fem else float("nan")
        eb_err = STEPPED.rel_error(eb, lambda_3d)
        timo_err = STEPPED.rel_error(timo, lambda_3d)
        warning = ""
        if not fem:
            warning = "FEM mode unavailable"
        elif not str(fem.get("mode_character", "")).startswith(("bending_", "axial_", "torsion_")):
            warning = "3D mode is mixed or unclassified"
        rows.append(
            {
                "mode_index": index + 1,
                "Lambda_EB": eb,
                "Lambda_Timoshenko": timo,
                "Lambda_3D": lambda_3d,
                "rel_error_EB_vs_3D": eb_err,
                "rel_error_Timo_vs_3D": timo_err,
                "closer_model": STEPPED.closer_model(eb_err, timo_err),
                "mode_character": fem.get("mode_character", "unavailable"),
                "warning": warning,
                "notes": "sorted first-8 diagnostic; 3D bending doublets and torsion/axial modes are not removed",
            }
        )
    return rows


def pair_average_rows(
    bending_rows: Sequence[dict[str, object]],
    bending_pairs: Sequence[tuple[dict[str, object], dict[str, object], str]],
    pair_count: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for pair_index in range(int(pair_count)):
        analytic = bending_rows[pair_index] if pair_index < len(bending_rows) else {}
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
        eb_err = STEPPED.rel_error(eb, avg)
        timo_err = STEPPED.rel_error(timo, avg)
        warning = pair_warning
        if not (first and second):
            warning = "FEM bending doublet unavailable under close-pair transverse criteria"
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
                "closer_model": STEPPED.closer_model(eb_err, timo_err),
                "warning": warning,
                "notes": (
                    f"FEM raw modes {first.get('raw_mode_number', '')} and "
                    f"{second.get('raw_mode_number', '')}; analytic values are single-rod bending roots"
                ),
            }
        )
    return rows


def is_bending_or_transverse(row: dict[str, object]) -> bool:
    character = str(row.get("mode_character", ""))
    if character.startswith("bending_"):
        return True
    if finite(row.get("transverse_fraction")) and finite(row.get("axial_fraction")):
        return float(row["transverse_fraction"]) >= 0.55 and float(row["axial_fraction"]) <= 0.65
    return False


def identify_uniform_bending_doublet_pairs(
    raw_rows: Sequence[dict[str, object]],
    *,
    pair_count: int,
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
        if split_rel <= float(split_rel_tol):
            first_ok = is_bending_or_transverse(first)
            second_ok = is_bending_or_transverse(second)
            if first_ok and second_ok:
                warning = ""
            elif first_ok or second_ok:
                warning = (
                    "one member not cleanly classified as bending; accepted by close-doublet "
                    "pattern for the circular uniform cylinder"
                )
            else:
                warning = (
                    "classification does not identify either member as bending; accepted by "
                    "close-doublet pattern only"
                )
            pairs.append((first, second, warning))
            index += 2
            continue
        index += 1
    return pairs


def finite_values(rows: Sequence[dict[str, object]], key: str) -> list[float]:
    return [float(row[key]) for row in rows if finite(row.get(key))]


def mean_max(values: Sequence[float]) -> tuple[float, float]:
    finite_only = [float(value) for value in values if math.isfinite(float(value))]
    if not finite_only:
        return float("nan"), float("nan")
    return float(sum(finite_only) / len(finite_only)), float(max(finite_only))


def plot_first8(paths: OutputPaths, geometry: Geometry, rows: Sequence[dict[str, object]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x_values = [int(row["mode_index"]) for row in rows]
    eb = [float(row["Lambda_EB"]) for row in rows]
    timo = [float(row["Lambda_Timoshenko"]) for row in rows]
    fem = [float(row["Lambda_3D"]) for row in rows]
    fig, ax = plt.subplots(figsize=(7.0, 4.3), constrained_layout=True)
    ax.plot(x_values, eb, marker="o", label="EB sorted in-plane")
    ax.plot(x_values, timo, marker="s", label="Timoshenko sorted in-plane")
    if any(math.isfinite(value) for value in fem):
        ax.plot(x_values, fem, marker="^", label="3D FEM sorted")
    ax.set_xlabel("mode index")
    ax.set_ylabel("Lambda")
    ax.set_title(f"Uniform beta=0, epsilon={geometry.epsilon:g}: first 8 sorted")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(paths.first8_png, dpi=180)
    plt.close(fig)

    eb_err = [float(row["rel_error_EB_vs_3D"]) for row in rows]
    timo_err = [float(row["rel_error_Timo_vs_3D"]) for row in rows]
    fig, ax = plt.subplots(figsize=(7.0, 4.3), constrained_layout=True)
    if any(math.isfinite(value) for value in eb_err):
        ax.plot(x_values, eb_err, marker="o", label="EB vs 3D")
    if any(math.isfinite(value) for value in timo_err):
        ax.plot(x_values, timo_err, marker="s", label="Timoshenko vs 3D")
    ax.set_xlabel("mode index")
    ax.set_ylabel("relative error")
    ax.set_title(f"Uniform beta=0, epsilon={geometry.epsilon:g}: first 8 errors")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(paths.first8_error_png, dpi=180)
    plt.close(fig)


def plot_pair_average(paths: OutputPaths, geometry: Geometry, rows: Sequence[dict[str, object]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x_values = [int(row["bending_pair_index"]) for row in rows]
    eb = [float(row["Lambda_EB"]) for row in rows]
    timo = [float(row["Lambda_Timoshenko"]) for row in rows]
    fem = [float(row["Lambda_3D_pair_average"]) for row in rows]
    fig, ax = plt.subplots(figsize=(7.0, 4.3), constrained_layout=True)
    ax.plot(x_values, eb, marker="o", label="EB bending")
    ax.plot(x_values, timo, marker="s", label="Timoshenko bending")
    if any(math.isfinite(value) for value in fem):
        ax.plot(x_values, fem, marker="^", label="3D FEM pair avg")
    ax.set_xlabel("bending pair index")
    ax.set_ylabel("Lambda")
    ax.set_title(f"Uniform beta=0, epsilon={geometry.epsilon:g}: bending pairs")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(paths.pair_png, dpi=180)
    plt.close(fig)

    eb_err = [float(row["rel_error_EB_vs_3D_pair_avg"]) for row in rows]
    timo_err = [float(row["rel_error_Timo_vs_3D_pair_avg"]) for row in rows]
    fig, ax = plt.subplots(figsize=(7.0, 4.3), constrained_layout=True)
    if any(math.isfinite(value) for value in eb_err):
        ax.plot(x_values, eb_err, marker="o", label="EB vs 3D pair avg")
    if any(math.isfinite(value) for value in timo_err):
        ax.plot(x_values, timo_err, marker="s", label="Timoshenko vs 3D pair avg")
    ax.set_xlabel("bending pair index")
    ax.set_ylabel("relative error")
    ax.set_title(f"Uniform beta=0, epsilon={geometry.epsilon:g}: pair-average errors")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(paths.pair_error_png, dpi=180)
    plt.close(fig)


def first_values_text(rows: Sequence[dict[str, object]], key: str, limit: int = 8) -> str:
    values = finite_values(rows, key)[: int(limit)]
    return ", ".join(f"{value:.8g}" for value in values) if values else "none"


def pair_summary(rows: Sequence[dict[str, object]]) -> tuple[str, float, float, float, float]:
    finite_rows = [row for row in rows if finite(row.get("Lambda_3D_pair_average"))]
    eb_mean, eb_max = mean_max([float(row["rel_error_EB_vs_3D_pair_avg"]) for row in finite_rows])
    timo_mean, timo_max = mean_max([float(row["rel_error_Timo_vs_3D_pair_avg"]) for row in finite_rows])
    timo_count = sum(1 for row in finite_rows if row.get("closer_model") == "Timoshenko")
    eb_count = sum(1 for row in finite_rows if row.get("closer_model") == "Euler-Bernoulli")
    ties = sum(1 for row in finite_rows if row.get("closer_model") == "tie")
    text = f"Timoshenko closer for {timo_count}/{len(finite_rows)} pairs; EB closer for {eb_count}; ties {ties}."
    return text, eb_mean, eb_max, timo_mean, timo_max


def write_report(
    *,
    path: Path,
    args: argparse.Namespace,
    gmsh_resolution: object,
    ccx_resolution: object,
    gmsh_version: str,
    ccx_version: str,
    case_records: Sequence[dict[str, object]],
    generated_files: Sequence[Path],
    timo_warnings_by_eps: dict[float, tuple[str, ...]],
) -> None:
    lines: list[str] = [
        "# EB / Timoshenko / 3D FEM Uniform beta=0 eps=0.05 Diagnostic",
        "",
        "Diagnostic-only comparison for a straight uniform fixed-fixed circular cylinder.",
        "",
        "## Scope",
        "",
        f"- workflow: `{WORKFLOW}`",
        "- beta: 0 deg only",
        "- mu: 0",
        "- eta: 0",
        f"- poisson: {args.poisson:g}",
        "- E: 1",
        "- rho: 1",
        "- total length L: 2",
        "- l0: 1",
        "- r0 = 2*epsilon",
        "- L1 = L2 = 1; tau1 = tau2 = 1; r1 = r2 = r0",
        "- 3D conversion confirmed in this workflow: Lambda = sqrt(omega/epsilon), with omega from the CalculiX RAD/TIME column when available.",
        "",
        "## Workflow Availability",
        "",
        "- Existing straight uniform fixed-fixed cylinder Gmsh + CalculiX helpers were reused.",
        "- This wrapper is diagnostic-only and uses the single-cylinder helper, not the stepped-cylinder geometry writer.",
        "- The input/output contract is different from the stepped audit: it reports both first-8 sorted low spectrum and bending doublet pair averages for a single uniform cylinder.",
        "",
        "## Executables",
        "",
        f"- gmsh path: `{gmsh_resolution.path}`",
        f"- gmsh source: {gmsh_resolution.source}",
        f"- gmsh probe: {gmsh_version}",
        f"- ccx path: `{ccx_resolution.path}`",
        f"- ccx source: {ccx_resolution.source}",
        f"- ccx probe: {ccx_version}",
        "",
        "## Geometry And Mesh",
        "",
        "| epsilon | radius | L | h | nodes | elems | components | warnings |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for record in case_records:
        geometry: Geometry = record["geometry"]  # type: ignore[assignment]
        fem: FemResult = record["fem"]  # type: ignore[assignment]
        mesh = fem.mesh_info
        warnings = str(mesh.get("warnings", "")) or "none"
        lines.append(
            f"| {geometry.epsilon:g} | {geometry.radius:.8g} | {geometry.total_length:g} | "
            f"{float(mesh.get('h', float('nan'))):.8g} | {int(mesh.get('nodes', 0))} | "
            f"{int(mesh.get('elements', 0))} | {int(mesh.get('connected_components', 0))} | {warnings} |"
        )
    lines.extend(
        [
            "",
            "## Frequency Request And Parsing",
            "",
            f"- requested CalculiX modes per 3D case: {args.n_fem_modes}",
            "",
            "| epsilon | FEM status | parsed raw modes | parse source | first 8 sorted 3D Lambdas |",
            "| ---: | --- | ---: | --- | --- |",
        ]
    )
    for record in case_records:
        geometry = record["geometry"]
        fem = record["fem"]
        raw_rows = fem.raw_rows
        parsed = sum(1 for row in raw_rows if finite(row.get("Lambda_3D")))
        lines.append(
            f"| {geometry.epsilon:g} | {fem.status} | {parsed} | {fem.parse_source} | "
            f"{first_values_text(raw_rows, 'Lambda_3D')} |"
        )
    lines.extend(["", "## Analytic Frequencies", ""])
    for record in case_records:
        geometry = record["geometry"]
        analytic = record["analytic_rows"]
        bending = record["bending_rows"]
        lines.extend(
            [
                f"### epsilon={geometry.epsilon:g}",
                "",
                f"- EB sorted first 8: {first_values_text(analytic, 'Lambda_EB')}",
                f"- Timoshenko sorted first 8: {first_values_text(analytic, 'Lambda_Timoshenko')}",
                f"- EB bending first 4: {first_values_text(bending, 'Lambda_EB', 4)}",
                f"- Timoshenko bending first 4: {first_values_text(bending, 'Lambda_Timoshenko', 4)}",
                f"- Timoshenko root warnings: {len(timo_warnings_by_eps.get(geometry.epsilon, ())) }",
                "",
            ]
        )
    lines.extend(["## Low First-8 Sorted Comparison", ""])
    for record in case_records:
        geometry = record["geometry"]
        rows = record["first8_rows"]
        lines.extend(
            [
                f"### epsilon={geometry.epsilon:g}",
                "",
                "| mode | EB | Timo | 3D | EB rel err | Timo rel err | closer | 3D character | warning |",
                "| ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |",
            ]
        )
        for row in rows:
            lines.append(
                f"| {int(row['mode_index'])} | {float(row['Lambda_EB']):.8g} | "
                f"{float(row['Lambda_Timoshenko']):.8g} | {float(row['Lambda_3D']):.8g} | "
                f"{float(row['rel_error_EB_vs_3D']):.6g} | {float(row['rel_error_Timo_vs_3D']):.6g} | "
                f"{row['closer_model']} | {row['mode_character']} | {row['warning']} |"
            )
        lines.append("")
    lines.extend(
        [
            "## Bending Pair-Averaged Comparison",
            "",
            "Pair averages use close adjacent 3D FEM pairs in the circular-cylinder doublet pattern.",
            "FRD displacement metrics are retained as warnings when one member is not cleanly classified as bending.",
            "This comparison is the cleaner diagnostic for bending behavior.",
            "",
        ]
    )
    for record in case_records:
        geometry = record["geometry"]
        rows = record["pair_rows"]
        summary, eb_mean, eb_max, timo_mean, timo_max = pair_summary(rows)
        lines.extend(
            [
                f"### epsilon={geometry.epsilon:g}",
                "",
                summary,
                f"Mean/max EB relative error: {eb_mean:.6g} / {eb_max:.6g}.",
                f"Mean/max Timoshenko relative error: {timo_mean:.6g} / {timo_max:.6g}.",
                "",
                "| pair | EB bending | Timo bending | 3D avg | EB rel err | Timo rel err | closer | split rel | warning |",
                "| ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | --- |",
            ]
        )
        for row in rows:
            lines.append(
                f"| {int(row['bending_pair_index'])} | {float(row['Lambda_EB']):.8g} | "
                f"{float(row['Lambda_Timoshenko']):.8g} | {float(row['Lambda_3D_pair_average']):.8g} | "
                f"{float(row['rel_error_EB_vs_3D_pair_avg']):.6g} | "
                f"{float(row['rel_error_Timo_vs_3D_pair_avg']):.6g} | {row['closer_model']} | "
                f"{float(row['doublet_split_rel']):.6g} | {row['warning']} |"
            )
        lines.append("")
    lines.extend(
        [
            "## Questions",
            "",
            "1. Case run: straight uniform fixed-fixed circular cylinder with beta=0, mu=0, eta=0, epsilon=0.05 by default.",
            "2. It is straight uniform beta=0; no stepped geometry is used.",
            "3. Mesh and solver settings are tabulated above.",
            f"4. Modes requested: {args.n_fem_modes}; parsed counts are tabulated above.",
            "5. First 8 EB frequencies are listed in the analytic section.",
            "6. First 8 Timoshenko frequencies are listed in the analytic section.",
            "7. First 8 3D FEM frequencies are listed in the parsing section.",
            "8. The first-8 plot is useful as a low-spectrum diagnostic but mixes 3D bending doublets with any torsion/axial modes that enter sorted order.",
            "9. Bending doublets are identified by close adjacent pair splitting; FRD displacement metrics are used for warnings, and inconsistent member classifications are recorded.",
            "10. Pair-average closer-model decisions and mean/max relative errors are tabulated above.",
            "11. At epsilon=0.05, the EB/Timoshenko difference is visually clear in bending pair averages, while 3D effects beyond Timoshenko can already be visible.",
            "12. This uniform case is clear enough as a diagnostic baseline before returning to the stepped case, but it is not final article validation.",
            "",
            "## Safety Confirmations",
            "",
            "- beta=0 only.",
            "- mu=0 and eta=0 only.",
            "- no stepped-cylinder 3D FEM.",
            "- no coupled-angle 3D FEM cases.",
            "- no beta scans.",
            "- no analytic formulas, determinant entries, root solvers, old solvers, or baseline results changed.",
            "- Timoshenko shear coefficient k' is unchanged; this workflow only calls existing helpers.",
            "- article workspace, main.tex, and article figures were not touched.",
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
        description="Diagnostic-only EB/Timoshenko/3D-FEM comparison for a beta=0 straight uniform cylinder."
    )
    parser.add_argument("--poisson", type=float, default=DEFAULT_POISSON)
    parser.add_argument("--n-analytic-roots", type=int, default=DEFAULT_N_ANALYTIC_ROOTS)
    parser.add_argument("--n-fem-modes", type=int, default=DEFAULT_N_FEM_MODES)
    parser.add_argument("--mesh-size", type=float, default=DEFAULT_MESH_SIZE)
    parser.add_argument("--lambda-max", type=float, default=DEFAULT_LAMBDA_MAX)
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--gmsh-exe", default=None)
    parser.add_argument("--ccx-exe", default=None)
    parser.add_argument("--timeout-seconds", type=int, default=DEFAULT_TIMEOUT_SECONDS)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--run-3d-fem", action="store_true")
    parser.add_argument("--skip-3d-fem", action="store_true")
    parser.add_argument("--also-run-eps0p1", action="store_true")
    args = parser.parse_args(argv)
    if args.smoke:
        args.run_3d_fem = False
        if args.output_dir == DEFAULT_OUTPUT_DIR:
            args.output_dir = Path("results") / "_smoke" / "uniform_beta0_eps0p05"
    if args.skip_3d_fem:
        args.run_3d_fem = False
    if args.run_3d_fem and args.n_fem_modes < 40:
        raise ValueError("3D FEM validation requires at least 40 requested modes.")
    if args.n_analytic_roots < 8:
        raise ValueError("Need at least 8 analytic roots for the requested output.")
    if args.mesh_size <= 0.0:
        raise ValueError("mesh size must be positive.")
    args.output_dir = repo_path(Path(args.output_dir))
    return args


def run(args: argparse.Namespace) -> dict[str, object]:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    single = STEPPED.load_single_module()
    extraction = STEPPED.load_extraction_module()
    gmsh_resolution = STEPPED.resolve_executable(
        explicit=args.gmsh_exe,
        env_name="GMSH_EXE",
        path_names=("gmsh.exe", "gmsh"),
        default_candidates=STEPPED.DEFAULT_GMSH_CANDIDATES,
    )
    ccx_resolution = STEPPED.resolve_executable(
        explicit=args.ccx_exe,
        env_name="CCX_EXE",
        path_names=("ccx.exe", "ccx_MT.exe", "ccx_static.exe", "ccx_dynamic.exe", "ccx"),
        default_candidates=STEPPED.DEFAULT_CCX_CANDIDATES,
    )
    gmsh_version = STEPPED.probe_version(gmsh_resolution.path, "--version")
    ccx_version = STEPPED.probe_version(ccx_resolution.path, "-v")

    epsilons = [DEFAULT_EPSILON]
    if args.also_run_eps0p1:
        epsilons.append(OPTIONAL_STRESS_EPSILON)

    generated_files: list[Path] = []
    case_records: list[dict[str, object]] = []
    timo_warnings_by_eps: dict[float, tuple[str, ...]] = {}

    for epsilon in epsilons:
        geometry = Geometry(
            epsilon=float(epsilon),
            h=float(args.mesh_size),
            total_length=DEFAULT_TOTAL_LENGTH,
            radius=2.0 * float(epsilon),
            poisson=float(args.poisson),
        )
        outputs = output_paths(output_dir, geometry.epsilon)
        bending_rows = bending_reference_rows(
            geometry.epsilon,
            float(args.poisson),
            DEFAULT_BENDING_PAIRS,
        )
        analytic, timo_warnings = analytic_rows(
            geometry=geometry,
            n_roots=int(args.n_analytic_roots),
            lambda_max=float(args.lambda_max),
            scan_step=float(args.root_scan_step),
            bending_rows=bending_rows,
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
        bending_pairs = identify_uniform_bending_doublet_pairs(
            fem.raw_rows,
            pair_count=DEFAULT_BENDING_PAIRS,
        )
        first8 = first8_rows(analytic, fem.raw_rows, count=8)
        pair_rows = pair_average_rows(
            bending_rows,
            bending_pairs,
            DEFAULT_BENDING_PAIRS,
        )

        write_csv(
            outputs.analytic_csv,
            analytic,
            ["sorted_index", "Lambda_EB", "Lambda_Timoshenko", "rel_diff_Timo_vs_EB", "analytic_mode_note"],
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
            outputs.first8_csv,
            first8,
            [
                "mode_index",
                "Lambda_EB",
                "Lambda_Timoshenko",
                "Lambda_3D",
                "rel_error_EB_vs_3D",
                "rel_error_Timo_vs_3D",
                "closer_model",
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
        plot_first8(outputs, geometry, first8)
        plot_pair_average(outputs, geometry, pair_rows)
        generated_files.extend(
            [
                outputs.analytic_csv,
                outputs.raw_fem_csv,
                outputs.first8_csv,
                outputs.pair_csv,
                outputs.first8_png,
                outputs.first8_error_png,
                outputs.pair_png,
                outputs.pair_error_png,
                *fem.case_files,
            ]
        )
        case_records.append(
            {
                "geometry": geometry,
                "outputs": outputs,
                "analytic_rows": analytic,
                "bending_rows": bending_rows,
                "fem": fem,
                "first8_rows": first8,
                "pair_rows": pair_rows,
            }
        )

    report_path = output_dir / "eb_timo_3d_uniform_beta0_mu0_eta0_eps0p05_report.md"
    write_report(
        path=report_path,
        args=args,
        gmsh_resolution=gmsh_resolution,
        ccx_resolution=ccx_resolution,
        gmsh_version=gmsh_version,
        ccx_version=ccx_version,
        case_records=case_records,
        generated_files=[*generated_files, report_path],
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
