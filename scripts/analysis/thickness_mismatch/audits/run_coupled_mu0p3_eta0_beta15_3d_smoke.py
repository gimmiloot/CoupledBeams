from __future__ import annotations

import argparse
import csv
import importlib.util
import math
from pathlib import Path
import sys
from types import ModuleType
from typing import Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import find_first_n_roots_eta  # noqa: E402
from my_project.analytic.solvers_out_of_plane import find_first_n_roots_out_of_plane_with_warnings  # noqa: E402


CASE_ID = "coupled_mu0p3_eta0_beta15"
DEFAULT_BETA_DEG = 15.0
DEFAULT_MU = 0.3
DEFAULT_ETA = 0.0
DEFAULT_EPSILON = 0.0025
DEFAULT_POISSON = 0.3
DEFAULT_MESH_SIZE = 0.004
DEFAULT_REQUESTED_MODES = 60
DEFAULT_COMPARE_MODES = 20
DEFAULT_ANALYTIC_ROOTS_PER_SUBSYSTEM = 24
DEFAULT_LAMBDA_MAX = 45.0
DEFAULT_ROOT_SCAN_STEP = 0.02
DEFAULT_OUTPUT_DIR = (
    Path("results")
    / "full_spectrum_analytic_vs_3d_fem"
    / "coupled_mu0p3_eta0_beta15"
)
WORKFLOW = "variable_length_fused_cylinders_reusing_equal_rod_extraction_helpers"
E = 1.0
RHO = 1.0


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def _fmt(value: object) -> object:
    if isinstance(value, float):
        if not math.isfinite(value):
            return "nan"
        return f"{value:.16e}"
    return value


def write_csv(path: Path, rows: Sequence[dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row.get(key, "")) for key in fieldnames})


def rel(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


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
    if length <= 0.0:
        raise ValueError("zero vector cannot be normalized")
    return (vector[0] / length, vector[1] / length, vector[2] / length)


def case_lengths(mu: float) -> tuple[float, float]:
    return 1.0 - float(mu), 1.0 + float(mu)


def radius_from_epsilon(epsilon: float) -> float:
    return 2.0 * float(epsilon)


def endpoint_geometry(beta_deg: float, mu: float) -> dict[str, tuple[float, float, float] | float]:
    beta = math.radians(float(beta_deg))
    cos_b = math.cos(beta)
    sin_b = math.sin(beta)
    l1, l2 = case_lengths(float(mu))
    joint = (0.0, 0.0, 0.0)
    rod1_outer = (-l1, 0.0, 0.0)
    rod2_outer = (-l2 * cos_b, -l2 * sin_b, 0.0)
    rod1_tangent = unit(sub(joint, rod1_outer))
    rod2_tangent = unit(sub(joint, rod2_outer))
    return {
        "beta_rad": beta,
        "joint": joint,
        "rod1_outer": rod1_outer,
        "rod2_outer": rod2_outer,
        "rod1_tangent": rod1_tangent,
        "rod2_tangent": rod2_tangent,
    }


def projection_on_axis(
    point: tuple[float, float, float],
    outer: tuple[float, float, float],
    tangent: tuple[float, float, float],
) -> tuple[float, float]:
    s = dot(sub(point, outer), tangent)
    axis_point = add(outer, mul(s, tangent))
    return s, norm(sub(point, axis_point))


def format_id_lines(ids: Sequence[int], *, per_line: int = 16) -> list[str]:
    lines: list[str] = []
    for start in range(0, len(ids), per_line):
        lines.append(", ".join(str(item) for item in ids[start : start + per_line]))
    return lines


def ccx_float(value: float) -> str:
    return f"{float(value):.12E}"


def parse_mesh_quality(messages: Sequence[str]) -> tuple[int, float, tuple[str, ...]]:
    negative_count = 0
    worst = float("nan")
    warning_lines: list[str] = []
    for message in messages:
        lower = message.lower()
        if "worst distortion" in lower or "jac." in lower or "distortion" in lower:
            warning_lines.append(message)
        if "worst distortion" not in lower:
            continue
        try:
            after = lower.split("worst distortion", 1)[1].split("=", 1)[1]
            worst = float(after.split()[0].strip("(),;"))
        except (IndexError, ValueError):
            pass
        if "elements with jac." in lower:
            try:
                before = lower.split("elements with jac.", 1)[0]
                negative_count = int(before.split()[-1])
            except (IndexError, ValueError):
                pass
    return negative_count, worst, tuple(warning_lines)


def build_analytic_rows(args: argparse.Namespace, compare: ModuleType) -> tuple[list[dict[str, object]], tuple[str, ...]]:
    beta_rad = math.radians(float(args.beta_deg))
    in_plane = find_first_n_roots_eta(
        beta=beta_rad,
        mu=float(args.mu),
        epsilon=float(args.epsilon),
        eta=float(args.eta),
        n_roots=int(args.n_analytic_roots_per_subsystem),
        Lmax0=float(args.lambda_max),
        scan_step=float(args.root_scan_step),
    )
    out = find_first_n_roots_out_of_plane_with_warnings(
        beta=beta_rad,
        mu=float(args.mu),
        epsilon=float(args.epsilon),
        eta=float(args.eta),
        poisson=float(args.poisson),
        n_roots=int(args.n_analytic_roots_per_subsystem),
        lambda_max=float(args.lambda_max),
        scan_step=float(args.root_scan_step),
    )
    union = compare.build_analytic_full_spectrum_union(in_plane, out.roots)
    rows: list[dict[str, object]] = []
    for row in union:
        rows.append(
            {
                "case_id": CASE_ID,
                "beta_deg": float(args.beta_deg),
                "mu": float(args.mu),
                "eta": float(args.eta),
                "epsilon": float(args.epsilon),
                "poisson": float(args.poisson),
                "analytic_full_index": int(row["analytic_full_index"]),
                "analytic_subsystem": row["analytic_subsystem"],
                "analytic_subsystem_index": int(row["analytic_subsystem_index"]),
                "Lambda": float(row["Lambda"]),
                "notes": "out-of-plane root warning" if out.warnings and row["analytic_subsystem"] == "out_of_plane" else "ok",
            }
        )
    warnings = tuple(f"out_of_plane: {item}" for item in out.warnings)
    return rows, warnings


def write_variable_length_geo(
    *,
    paths: object,
    beta_deg: float,
    mu: float,
    epsilon: float,
    mesh_size: float,
    boolean_mode: str,
) -> str:
    geom = endpoint_geometry(beta_deg, mu)
    l1, l2 = case_lengths(mu)
    radius = radius_from_epsilon(epsilon)
    cos_b = math.cos(float(geom["beta_rad"]))
    sin_b = math.sin(float(geom["beta_rad"]))
    bbox = max(l1, l2) + 4.0 * radius
    paths.case_dir.mkdir(parents=True, exist_ok=True)
    if boolean_mode == "union":
        sphere_line = ""
        boolean_line = "fused[] = BooleanUnion{ Volume{1}; Delete; }{ Volume{2}; Delete; };"
        label = "fused variable-length cylinders (OpenCASCADE BooleanUnion)"
    else:
        sphere_line = "Sphere(3) = {0, 0, 0, R};"
        boolean_line = "fused[] = BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Volume{3}; Delete; };"
        label = "fused variable-length cylinders + spherical joint (OpenCASCADE BooleanFragments fallback)"
    paths.geo.write_text(
        f"""// Diagnostic-only variable-length coupled cylinder smoke case.
// beta = {float(beta_deg):.17g} deg, mu = {float(mu):.17g}, eta = 0, epsilon = {float(epsilon):.17g}
SetFactory("OpenCASCADE");

L1 = {l1:.17g};
L2 = {l2:.17g};
R = {radius:.17g};
h = {float(mesh_size):.17g};
CosB = {cos_b:.17g};
SinB = {sin_b:.17g};

Cylinder(1) = {{-L1, 0, 0, L1, 0, 0, R, 2*Pi}};
Cylinder(2) = {{-L2*CosB, -L2*SinB, 0, L2*CosB, L2*SinB, 0, R, 2*Pi}};
{sphere_line}

{boolean_line}

vols[] = Volume In BoundingBox {{-L2-{bbox:.17g}, -L2-{bbox:.17g}, -R-{bbox:.17g}, R+{bbox:.17g}, R+{bbox:.17g}, R+{bbox:.17g}}};
Physical Volume("SOLID", 1) = vols[];

Mesh.CharacteristicLengthMin = h;
Mesh.CharacteristicLengthMax = h;
Mesh.ElementOrder = 2;
Mesh.MshFileVersion = 4.1;
""",
        encoding="utf-8",
    )
    return label


def generate_mesh(
    *,
    coupled: ModuleType,
    paths: object,
    beta_deg: float,
    mu: float,
    epsilon: float,
    mesh_size: float,
    gmsh_exe: str,
) -> tuple[bool, str, str, tuple[str, ...]]:
    label = write_variable_length_geo(
        paths=paths,
        beta_deg=beta_deg,
        mu=mu,
        epsilon=epsilon,
        mesh_size=mesh_size,
        boolean_mode="union",
    )
    ok, message, messages = coupled.generate_mesh_with_gmsh_cli(paths, gmsh_exe)
    if ok:
        return ok, message, label, messages
    fallback_label = write_variable_length_geo(
        paths=paths,
        beta_deg=beta_deg,
        mu=mu,
        epsilon=epsilon,
        mesh_size=mesh_size,
        boolean_mode="fragments",
    )
    fallback_ok, fallback_message, fallback_messages = coupled.generate_mesh_with_gmsh_cli(paths, gmsh_exe)
    all_messages = tuple(list(messages) + [f"union attempt failed: {message}"] + list(fallback_messages))
    return fallback_ok, fallback_message, fallback_label, all_messages


def fixed_node_sets(mesh: object, *, beta_deg: float, mu: float, epsilon: float) -> tuple[list[int], list[int]]:
    geom = endpoint_geometry(beta_deg, mu)
    l1, l2 = case_lengths(mu)
    radius = radius_from_epsilon(epsilon)
    plane_tol = max(1.0e-7 * max(l1, l2), 1.0e-6)
    radial_tol = 1.02 * radius + plane_tol
    rod1_nodes: list[int] = []
    rod2_nodes: list[int] = []
    rod1_outer = geom["rod1_outer"]
    rod2_outer = geom["rod2_outer"]
    rod1_tangent = geom["rod1_tangent"]
    rod2_tangent = geom["rod2_tangent"]
    assert isinstance(rod1_outer, tuple)
    assert isinstance(rod2_outer, tuple)
    assert isinstance(rod1_tangent, tuple)
    assert isinstance(rod2_tangent, tuple)
    for node_id, point in mesh.nodes.items():
        s1, r1 = projection_on_axis(point, rod1_outer, rod1_tangent)
        s2, r2 = projection_on_axis(point, rod2_outer, rod2_tangent)
        if abs(s1) <= plane_tol and r1 <= radial_tol:
            rod1_nodes.append(node_id)
        if abs(s2) <= plane_tol and r2 <= radial_tol:
            rod2_nodes.append(node_id)
    return sorted(rod1_nodes), sorted(rod2_nodes)


def write_calculix_inputs(
    *,
    coupled: ModuleType,
    paths: object,
    mesh: object,
    beta_deg: float,
    mu: float,
    epsilon: float,
    poisson: float,
    requested_modes: int,
) -> dict[str, object]:
    rod1_nodes, rod2_nodes = fixed_node_sets(mesh, beta_deg=beta_deg, mu=mu, epsilon=epsilon)
    coupled.write_calculix_mesh_include(paths, mesh)
    paths.ccx_modal_inp.write_text(
        "\n".join(
            [
                "** Diagnostic-only CalculiX modal input for variable-length coupled cylinders.",
                "** Fixed node sets are generated from outer-end coordinate projections.",
                f"*INCLUDE, INPUT={paths.ccx_mesh_inp.name}",
                "*NSET, NSET=ROD1_OUTER_FIXED",
                *format_id_lines(rod1_nodes),
                "*NSET, NSET=ROD2_OUTER_FIXED",
                *format_id_lines(rod2_nodes),
                "*MATERIAL, NAME=MAT",
                "*ELASTIC",
                f"{E:.17g}, {float(poisson):.17g}",
                "*DENSITY",
                f"{RHO:.17g}",
                "*SOLID SECTION, ELSET=SOLID, MATERIAL=MAT",
                "*BOUNDARY",
                "ROD1_OUTER_FIXED, 1, 3, 0",
                "ROD2_OUTER_FIXED, 1, 3, 0",
                "*STEP",
                "*FREQUENCY",
                f"{int(requested_modes)}",
                "*NODE FILE",
                "U",
                "*END STEP",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "rod1_fixed_node_count": len(rod1_nodes),
        "rod2_fixed_node_count": len(rod2_nodes),
        "node_sets_present": bool(rod1_nodes and rod2_nodes),
        "mesh_include": paths.ccx_mesh_inp,
        "modal_input": paths.ccx_modal_inp,
    }


def classify_shape(compare: ModuleType, mode_shape: dict[int, tuple[float, float, float]]) -> tuple[str, float, float]:
    if not mode_shape:
        return "classification_unavailable", float("nan"), float("nan")
    values = np.asarray(list(mode_shape.values()), dtype=float)
    return compare.classify_3d_displacement_mode(values)


def build_fem_rows(
    *,
    coupled: ModuleType,
    compare: ModuleType,
    omegas: Sequence[float],
    parse_result: object | None,
    epsilon: float,
    status: str,
    warning: str,
) -> list[dict[str, object]]:
    shape_by_mode = getattr(parse_result, "mode_shapes", {}) if parse_result is not None else {}
    if not omegas:
        return [
            {
                "case_id": CASE_ID,
                "fem_mode_index": "",
                "ccx_mode_index": "",
                "Omega_3d_fem": float("nan"),
                "Lambda_3d_fem": float("nan"),
                "fem_mode_character": "classification_unavailable",
                "fem_in_plane_fraction": float("nan"),
                "fem_out_of_plane_fraction": float("nan"),
                "status": status,
                "warning": warning,
                "notes": "no parsed FEM frequencies",
            }
        ]
    raw: list[tuple[int, float, float]] = []
    for ccx_mode_index, omega in enumerate(omegas, start=1):
        Lambda = coupled.lambda_fem_from_omega(float(omega), float(epsilon))
        raw.append((ccx_mode_index, float(omega), Lambda))
    raw.sort(key=lambda item: item[2])
    rows: list[dict[str, object]] = []
    for sorted_index, (ccx_mode_index, omega, Lambda) in enumerate(raw, start=1):
        shape = shape_by_mode.get(ccx_mode_index, {})
        character, in_fraction, out_fraction = classify_shape(compare, shape)
        rows.append(
            {
                "case_id": CASE_ID,
                "fem_mode_index": sorted_index,
                "ccx_mode_index": ccx_mode_index,
                "Omega_3d_fem": omega,
                "Lambda_3d_fem": Lambda,
                "fem_mode_character": character,
                "fem_in_plane_fraction": in_fraction,
                "fem_out_of_plane_fraction": out_fraction,
                "status": status,
                "warning": warning,
                "notes": (
                    "classification from CalculiX FRD displacement vectors"
                    if shape
                    else "FRD displacement vectors unavailable for this mode"
                ),
            }
        )
    return rows


def sorted_match_rows(
    analytic_rows: Sequence[dict[str, object]],
    fem_rows: Sequence[dict[str, object]],
    compare_modes: int,
) -> list[dict[str, object]]:
    finite_fem = [
        row
        for row in fem_rows
        if str(row.get("fem_mode_index", "")).isdigit()
        and math.isfinite(float(row.get("Lambda_3d_fem", float("nan"))))
    ]
    rows: list[dict[str, object]] = []
    for index, analytic in enumerate(analytic_rows[: int(compare_modes)], start=1):
        fem = finite_fem[index - 1] if index <= len(finite_fem) else None
        Lambda_analytic = float(analytic["Lambda"])
        Lambda_fem = float("nan") if fem is None else float(fem["Lambda_3d_fem"])
        abs_error = abs(Lambda_fem - Lambda_analytic) if math.isfinite(Lambda_fem) else float("nan")
        rel_error = abs_error / abs(Lambda_analytic) if math.isfinite(abs_error) and Lambda_analytic != 0.0 else float("nan")
        rows.append(
            {
                "case_id": CASE_ID,
                "beta_deg": analytic["beta_deg"],
                "mu": analytic["mu"],
                "eta": analytic["eta"],
                "epsilon": analytic["epsilon"],
                "poisson": analytic["poisson"],
                "analytic_full_index": analytic["analytic_full_index"],
                "analytic_subsystem": analytic["analytic_subsystem"],
                "analytic_subsystem_index": analytic["analytic_subsystem_index"],
                "Lambda_analytic": Lambda_analytic,
                "fem_mode_index": "" if fem is None else fem["fem_mode_index"],
                "ccx_mode_index": "" if fem is None else fem["ccx_mode_index"],
                "Lambda_3d_fem": Lambda_fem,
                "abs_error": abs_error,
                "rel_error": rel_error,
                "fem_mode_character": "classification_unavailable" if fem is None else fem["fem_mode_character"],
                "fem_in_plane_fraction": float("nan") if fem is None else fem["fem_in_plane_fraction"],
                "fem_out_of_plane_fraction": float("nan") if fem is None else fem["fem_out_of_plane_fraction"],
                "match_method": "sorted_index" if fem is not None else "fem_mode_missing",
                "warning": "" if fem is not None else "No FEM mode at this sorted index.",
                "notes": "frequency-order smoke comparison only; do not use as branch identity",
            }
        )
    return rows


def nearest_rows(compare: ModuleType, analytic_rows: Sequence[dict[str, object]], fem_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    raw_nearest = compare.nearest_frequency_match_rows(analytic_rows, fem_rows)
    fem_by_index = {int(row["fem_mode_index"]): row for row in fem_rows if str(row.get("fem_mode_index", "")).isdigit()}
    enriched: list[dict[str, object]] = []
    for row in raw_nearest:
        fem = fem_by_index.get(int(row["fem_mode_index"])) if str(row.get("fem_mode_index", "")).isdigit() else None
        enriched.append(
            {
                **row,
                "ccx_mode_index": "" if fem is None else fem["ccx_mode_index"],
                "fem_mode_character": "classification_unavailable" if fem is None else fem["fem_mode_character"],
                "fem_in_plane_fraction": float("nan") if fem is None else fem["fem_in_plane_fraction"],
                "fem_out_of_plane_fraction": float("nan") if fem is None else fem["fem_out_of_plane_fraction"],
            }
        )
    return enriched


def error_summary(rows: Sequence[dict[str, object]], *, limit: int | None = None) -> tuple[float, float]:
    source = rows[:limit] if limit is not None else rows
    values = sorted(
        float(row["rel_error"])
        for row in source
        if math.isfinite(float(row.get("rel_error", float("nan"))))
    )
    if not values:
        return float("nan"), float("nan")
    midpoint = len(values) // 2
    median = values[midpoint] if len(values) % 2 else 0.5 * (values[midpoint - 1] + values[midpoint])
    return max(values), median


def first_values(rows: Sequence[dict[str, object]], key: str, *, limit: int) -> list[float]:
    values: list[float] = []
    for row in rows:
        try:
            value = float(row[key])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(value):
            values.append(value)
        if len(values) >= int(limit):
            break
    return values


def close_pair_summary(values: Sequence[float], *, rel_tol: float = 2.0e-3, limit: int = 20) -> str:
    pairs: list[str] = []
    limited = list(values[:limit])
    for index in range(len(limited) - 1):
        left = limited[index]
        right = limited[index + 1]
        mean = 0.5 * (left + right)
        if mean > 0.0 and abs(right - left) / mean <= rel_tol:
            pairs.append(f"{index + 1}-{index + 2} (rel split {abs(right - left) / mean:.3g})")
    return "; ".join(pairs) if pairs else "none"


def class_counts(fem_rows: Sequence[dict[str, object]]) -> str:
    counts: dict[str, int] = {}
    for row in fem_rows:
        label = str(row.get("fem_mode_character", ""))
        if label:
            counts[label] = counts.get(label, 0) + 1
    return ", ".join(f"{key}={counts[key]}" for key in sorted(counts)) if counts else "none"


def compatibility_summary(rows: Sequence[dict[str, object]], *, nearest: bool) -> str:
    total = 0
    compatible = 0
    mixed = 0
    unavailable = 0
    for row in rows:
        subsystem_key = "nearest_analytic_subsystem" if nearest else "analytic_subsystem"
        subsystem = str(row.get(subsystem_key, ""))
        character = str(row.get("fem_mode_character", ""))
        if character == "classification_unavailable":
            unavailable += 1
            continue
        if character == "mixed_or_unclassified":
            mixed += 1
            continue
        if character not in {"mostly_in_plane", "mostly_out_of_plane"}:
            continue
        total += 1
        if (subsystem == "in_plane" and character == "mostly_in_plane") or (
            subsystem == "out_of_plane" and character == "mostly_out_of_plane"
        ):
            compatible += 1
    if total == 0:
        return f"no definite compatible/incompatible labels; mixed={mixed}, unavailable={unavailable}"
    return f"{compatible}/{total} definite labels compatible; mixed={mixed}, unavailable={unavailable}"


def overlap_note(beta_deg: float, epsilon: float, mu: float) -> tuple[float, float, float]:
    radius = radius_from_epsilon(epsilon)
    sin_beta = abs(math.sin(math.radians(beta_deg)))
    overlap_length = float("inf") if sin_beta <= 1.0e-14 else radius / sin_beta
    l1, l2 = case_lengths(mu)
    return overlap_length, overlap_length / min(l1, l2), sin_beta


def write_report(
    *,
    path: Path,
    args: argparse.Namespace,
    analytic_rows: Sequence[dict[str, object]],
    fem_rows: Sequence[dict[str, object]],
    sorted_rows: Sequence[dict[str, object]],
    nearest_match_rows: Sequence[dict[str, object]],
    messages: Sequence[str],
    analytic_warnings: Sequence[str],
    generated_files: Sequence[Path],
    geometry_label: str,
    mesh: object | None,
    mesh_quality: tuple[int, float, tuple[str, ...]],
    fixed_summary: dict[str, object] | None,
    tool_audit: object,
    coupled: ModuleType,
    parse_message: str,
    actual_3d_run: bool,
) -> None:
    l1, l2 = case_lengths(float(args.mu))
    radius = radius_from_epsilon(float(args.epsilon))
    overlap_length, overlap_fraction, sin_beta = overlap_note(float(args.beta_deg), float(args.epsilon), float(args.mu))
    first_analytic_20 = first_values(analytic_rows, "Lambda", limit=20)
    first_fem_20 = first_values(fem_rows, "Lambda_3d_fem", limit=20)
    sorted_max, sorted_median = error_summary(sorted_rows)
    nearest_first20 = nearest_match_rows[: int(args.compare_modes)]
    nearest_max, nearest_median = error_summary(nearest_first20)
    ambiguity_count = sum(int(row.get("ambiguity_count_within_tolerance", 0)) > 1 for row in nearest_first20)
    negative_count, worst_distortion, mesh_warning_lines = mesh_quality
    components = (0, 0) if mesh is None else coupled.connected_component_counts(mesh)
    node_count = 0 if mesh is None else len(mesh.nodes)
    element_count = 0 if mesh is None else len(mesh.solid_elements)
    bbox_min = None if mesh is None else mesh.bbox_min
    bbox_max = None if mesh is None else mesh.bbox_max
    promising = (
        actual_3d_run
        and len(first_fem_20) >= int(args.compare_modes)
        and math.isfinite(nearest_median)
        and math.isfinite(sorted_median)
        and nearest_median < 2.0e-2
        and sorted_median < 2.0e-2
        and negative_count == 0
    )
    lines = [
        "# Coupled-Angle 3D FEM Smoke Comparison",
        "",
        "Diagnostic-only full-spectrum smoke comparison for one variable-length equal-thickness coupled-cylinder case.",
        "This is a 3D smoke/convergence baseline, not article validation.",
        "",
        "## Case",
        "",
        f"- beta: {float(args.beta_deg):g} deg",
        f"- mu: {float(args.mu):g}",
        f"- eta: {float(args.eta):g}",
        f"- epsilon: {float(args.epsilon):g}",
        f"- poisson: {float(args.poisson):g}",
        f"- lengths: L1={l1:.16g}, L2={l2:.16g}, L1+L2={l1 + l2:.16g}",
        f"- radius: r0={radius:.16g} (`2*epsilon` with base length l=1)",
        f"- material: E={E:g}, rho={RHO:g}",
        f"- mesh size h: {float(args.mesh_size):.16g}",
        f"- requested FEM modes: {int(args.n_fem_modes)}",
        f"- compared sorted/nearest modes: {int(args.compare_modes)}",
        "",
        "## Workflow",
        "",
        f"- selected workflow: `{WORKFLOW}`",
        f"- geometry: {geometry_label or 'not generated'}",
        "- reason: existing equal-rod 3D scripts cannot represent mu=0.3 variable lengths directly,",
        "  so this wrapper reuses their Gmsh/CalculiX parsing and extraction helpers only.",
        f"- Gmsh executable: `{tool_audit.gmsh_exe}`" if getattr(tool_audit, "gmsh_exe", None) else "- Gmsh executable: not found",
        f"- CalculiX executable: `{tool_audit.ccx_exe}`" if getattr(tool_audit, "ccx_exe", None) else "- CalculiX executable: not found",
        f"- actual 3D FEM run: {'yes' if actual_3d_run else 'no'}",
        "",
        "## Geometry Audit",
        "",
        f"- rod 1 outer end: {endpoint_geometry(float(args.beta_deg), float(args.mu))['rod1_outer']}",
        f"- rod 2 outer end: {endpoint_geometry(float(args.beta_deg), float(args.mu))['rod2_outer']}",
        f"- sin(beta): {sin_beta:.16g}",
        f"- estimated local fused-cylinder overlap length: {overlap_length:.16g}",
        f"- estimated overlap / min(length): {overlap_fraction:.16g}",
        "- overlap interpretation: the Boolean fused geometry creates a localized intersection near the shared joint; no deliberate overlap is introduced along the remote rod spans.",
        f"- mesh nodes: {node_count}",
        f"- solid elements: {element_count}",
        f"- connected components: {components[0]} (largest element count {components[1]})",
        f"- bounding box min: {bbox_min}",
        f"- bounding box max: {bbox_max}",
        f"- fixed rod-1 outer nodes: {0 if fixed_summary is None else fixed_summary['rod1_fixed_node_count']}",
        f"- fixed rod-2 outer nodes: {0 if fixed_summary is None else fixed_summary['rod2_fixed_node_count']}",
        f"- Gmsh negative-Jacobian elements: {negative_count}",
        f"- Gmsh worst distortion: {worst_distortion:.8g}" if math.isfinite(worst_distortion) else "- Gmsh worst distortion: not reported",
        "",
        "## Roots And Matching",
        "",
        f"- first 20 analytic roots: {', '.join(f'{value:.8g}' for value in first_analytic_20)}",
        f"- first 20 FEM roots: {', '.join(f'{value:.8g}' for value in first_fem_20) if first_fem_20 else 'unavailable'}",
        f"- sorted-index max relative error first {int(args.compare_modes)}: {sorted_max:.8g}" if math.isfinite(sorted_max) else f"- sorted-index max relative error first {int(args.compare_modes)}: unavailable",
        f"- sorted-index median relative error first {int(args.compare_modes)}: {sorted_median:.8g}" if math.isfinite(sorted_median) else f"- sorted-index median relative error first {int(args.compare_modes)}: unavailable",
        f"- nearest-frequency max relative error first {int(args.compare_modes)} FEM modes: {nearest_max:.8g}" if math.isfinite(nearest_max) else f"- nearest-frequency max relative error first {int(args.compare_modes)} FEM modes: unavailable",
        f"- nearest-frequency median relative error first {int(args.compare_modes)} FEM modes: {nearest_median:.8g}" if math.isfinite(nearest_median) else f"- nearest-frequency median relative error first {int(args.compare_modes)} FEM modes: unavailable",
        f"- nearest-frequency ambiguous duplicate/clustered matches in first {int(args.compare_modes)}: {ambiguity_count}",
        f"- analytic close adjacent pairs in first 20: {close_pair_summary(first_analytic_20)}",
        f"- FEM close adjacent pairs in first 20: {close_pair_summary(first_fem_20)}",
        "",
        "## Mode Classification",
        "",
        f"- classes present: {class_counts(fem_rows)}",
        f"- sorted-index analytic-label/3D-character compatibility: {compatibility_summary(sorted_rows, nearest=False)}",
        f"- nearest-frequency analytic-label/3D-character compatibility: {compatibility_summary(nearest_first20, nearest=True)}",
        f"- parser status: {parse_message}",
        "",
        "## Diagnostic Assessment",
        "",
        f"- requesting {int(args.n_fem_modes)} modes recovered at least {len(first_fem_20)} sorted FEM roots for the first-{int(args.compare_modes)} comparison.",
        f"- eta != 0 readiness: {'promising enough for a later eta smoke only, still diagnostic' if promising else 'not yet strong enough to treat as a clean eta != 0 baseline without further review'}",
        "- sorted-index matching is orientation only; nearest-frequency matching is included because duplicate and clustered roots can move sorted positions.",
        "- No analytic determinant, production FEM workflow, old solver, article workspace, or article figure was modified.",
        "",
        "## Messages",
        "",
    ]
    lines.extend([f"- {message}" for message in messages] or ["- none"])
    lines.extend(["", "## Mesh Warnings", ""])
    lines.extend([f"- {message}" for message in mesh_warning_lines] or ["- none"])
    lines.extend(["", "## Analytic Warnings", ""])
    lines.extend([f"- {warning}" for warning in analytic_warnings] or ["- none"])
    lines.extend(["", "## Generated Files", ""])
    for item in [*generated_files, path]:
        lines.append(f"- `{rel(item)}`")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the beta=15, mu=0.3, eta=0 coupled-angle 3D FEM smoke case.")
    parser.add_argument("--beta-deg", type=float, default=DEFAULT_BETA_DEG)
    parser.add_argument("--mu", type=float, default=DEFAULT_MU)
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--poisson", type=float, default=DEFAULT_POISSON)
    parser.add_argument("--mesh-size", type=float, default=DEFAULT_MESH_SIZE)
    parser.add_argument("--n-fem-modes", type=int, default=DEFAULT_REQUESTED_MODES)
    parser.add_argument("--compare-modes", type=int, default=DEFAULT_COMPARE_MODES)
    parser.add_argument("--n-analytic-roots-per-subsystem", type=int, default=DEFAULT_ANALYTIC_ROOTS_PER_SUBSYSTEM)
    parser.add_argument("--lambda-max", type=float, default=DEFAULT_LAMBDA_MAX)
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    parser.add_argument("--ccx-timeout-seconds", type=int, default=600)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args(argv)
    if args.eta != 0.0:
        raise ValueError("this diagnostic wrapper is intentionally scoped to eta=0 equal thickness.")
    if not (-1.0 < args.mu < 1.0):
        raise ValueError("mu must lie inside (-1, 1).")
    if args.epsilon <= 0.0 or args.mesh_size <= 0.0:
        raise ValueError("epsilon and mesh size must be positive.")
    if args.n_fem_modes <= 0 or args.compare_modes <= 0:
        raise ValueError("mode counts must be positive.")
    args.output_dir = _repo_path(args.output_dir)
    return args


def run(args: argparse.Namespace) -> dict[str, object]:
    output_dir = Path(args.output_dir)
    coupled = _load_module("solid_fem_coupled_equal_rods_for_variable_length_smoke", REPO_ROOT / "scripts" / "analysis" / "solid_fem_coupled_equal_rods.py")
    compare = _load_module("compare_full_spectrum_analytic_vs_3d_fem_for_coupled_smoke", REPO_ROOT / "scripts" / "analysis" / "thickness_mismatch" / "audits" / "compare_full_spectrum_analytic_vs_3d_fem.py")
    coupled.N_SOLID_MODES = int(args.n_fem_modes)
    coupled.NU = float(args.poisson)
    coupled.CCX_TIMEOUT_SECONDS = int(args.ccx_timeout_seconds)

    analytic_rows, analytic_warnings = build_analytic_rows(args, compare)

    case_dir = output_dir / "solid_fem_case"
    stem = "coupled_mu0p3_eta0_beta15_eps0p0025_h0p004"
    paths = coupled.CasePaths(
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
    audit = coupled.audit_tools()
    messages: list[str] = []
    actual_3d_run = False
    mesh = None
    geometry_label = ""
    mesh_quality = (0, float("nan"), tuple())
    fixed_summary: dict[str, object] | None = None
    parse_result = None
    parse_message = "mode-shape parsing not attempted"
    omegas: tuple[float, ...] = ()
    fem_status = "not_run"
    fem_warning = ""

    if not audit.gmsh_exe:
        fem_warning = "Gmsh executable not found"
        messages.append(fem_warning)
    elif not audit.ccx_exe:
        fem_warning = "CalculiX executable not found"
        messages.append(fem_warning)
    else:
        mesh_ok, mesh_message, geometry_label, gmsh_messages = generate_mesh(
            coupled=coupled,
            paths=paths,
            beta_deg=float(args.beta_deg),
            mu=float(args.mu),
            epsilon=float(args.epsilon),
            mesh_size=float(args.mesh_size),
            gmsh_exe=audit.gmsh_exe,
        )
        messages.append(f"wrote Gmsh geometry: {rel(paths.geo)}")
        messages.append(mesh_message)
        messages.extend(gmsh_messages)
        mesh_quality = parse_mesh_quality(gmsh_messages)
        if mesh_ok and paths.gmsh_inp.exists():
            mesh = coupled.read_gmsh_inp_mesh_data(paths.gmsh_inp)
            fixed_summary = write_calculix_inputs(
                coupled=coupled,
                paths=paths,
                mesh=mesh,
                beta_deg=float(args.beta_deg),
                mu=float(args.mu),
                epsilon=float(args.epsilon),
                poisson=float(args.poisson),
                requested_modes=int(args.n_fem_modes),
            )
            messages.append(
                f"wrote CalculiX input: rod1_fixed={fixed_summary['rod1_fixed_node_count']}, "
                f"rod2_fixed={fixed_summary['rod2_fixed_node_count']}"
            )
            if not fixed_summary["node_sets_present"]:
                fem_status = "run_failed_or_unparsed"
                fem_warning = "one or both fixed node sets are empty"
                messages.append(fem_warning)
            else:
                result = coupled.run_calculix(paths, audit.ccx_exe, float(args.beta_deg), float(args.epsilon))
                actual_3d_run = True
                fem_status = "computed" if result.success and result.parsed_omegas else "run_failed_or_unparsed"
                fem_warning = "" if result.success and result.parsed_omegas else result.message
                messages.append(result.message)
                omegas = tuple(result.parsed_omegas)
                if result.success and result.parsed_omegas:
                    parse_result = coupled.parse_calculix_frd_mode_shapes(paths, float(args.beta_deg), float(args.epsilon))
                    parse_message = getattr(parse_result, "message", "mode-shape parse result unavailable")
                    messages.append(parse_message)
        else:
            fem_status = "run_failed_or_unparsed"
            fem_warning = mesh_message

    fem_rows = build_fem_rows(
        coupled=coupled,
        compare=compare,
        omegas=omegas,
        parse_result=parse_result,
        epsilon=float(args.epsilon),
        status=fem_status,
        warning=fem_warning,
    )
    sorted_rows = sorted_match_rows(analytic_rows, fem_rows, int(args.compare_modes))
    nearest_match_rows = nearest_rows(compare, analytic_rows, fem_rows)

    analytic_path = output_dir / "analytic_union.csv"
    fem_path = output_dir / "fem_3d_raw_modes.csv"
    sorted_path = output_dir / "sorted_index_match.csv"
    nearest_path = output_dir / "nearest_frequency_match.csv"
    report_path = output_dir / "report.md"
    write_csv(
        analytic_path,
        analytic_rows,
        [
            "case_id",
            "beta_deg",
            "mu",
            "eta",
            "epsilon",
            "poisson",
            "analytic_full_index",
            "analytic_subsystem",
            "analytic_subsystem_index",
            "Lambda",
            "notes",
        ],
    )
    write_csv(
        fem_path,
        fem_rows,
        [
            "case_id",
            "fem_mode_index",
            "ccx_mode_index",
            "Omega_3d_fem",
            "Lambda_3d_fem",
            "fem_mode_character",
            "fem_in_plane_fraction",
            "fem_out_of_plane_fraction",
            "status",
            "warning",
            "notes",
        ],
    )
    write_csv(
        sorted_path,
        sorted_rows,
        [
            "case_id",
            "beta_deg",
            "mu",
            "eta",
            "epsilon",
            "poisson",
            "analytic_full_index",
            "analytic_subsystem",
            "analytic_subsystem_index",
            "Lambda_analytic",
            "fem_mode_index",
            "ccx_mode_index",
            "Lambda_3d_fem",
            "abs_error",
            "rel_error",
            "fem_mode_character",
            "fem_in_plane_fraction",
            "fem_out_of_plane_fraction",
            "match_method",
            "warning",
            "notes",
        ],
    )
    write_csv(
        nearest_path,
        nearest_match_rows,
        [
            "case_id",
            "fem_mode_index",
            "ccx_mode_index",
            "Lambda_3d_fem",
            "nearest_analytic_full_index",
            "nearest_analytic_subsystem",
            "nearest_analytic_subsystem_index",
            "Lambda_analytic_nearest",
            "abs_error",
            "rel_error",
            "ambiguity_count_within_tolerance",
            "fem_mode_character",
            "fem_in_plane_fraction",
            "fem_out_of_plane_fraction",
            "warning",
            "notes",
        ],
    )
    generated_files = [analytic_path, fem_path, sorted_path, nearest_path]
    write_report(
        path=report_path,
        args=args,
        analytic_rows=analytic_rows,
        fem_rows=fem_rows,
        sorted_rows=sorted_rows,
        nearest_match_rows=nearest_match_rows,
        messages=messages,
        analytic_warnings=analytic_warnings,
        generated_files=generated_files,
        geometry_label=geometry_label,
        mesh=mesh,
        mesh_quality=mesh_quality,
        fixed_summary=fixed_summary,
        tool_audit=audit,
        coupled=coupled,
        parse_message=parse_message,
        actual_3d_run=actual_3d_run,
    )
    print(f"saved analytic union CSV: {analytic_path}")
    print(f"saved raw FEM CSV: {fem_path}")
    print(f"saved sorted-index CSV: {sorted_path}")
    print(f"saved nearest-frequency CSV: {nearest_path}")
    print(f"saved report: {report_path}")
    print(f"modes requested: {int(args.n_fem_modes)}")
    print(f"modes parsed: {len(first_values(fem_rows, 'Lambda_3d_fem', limit=1000))}")
    print(f"actual 3D FEM run: {'yes' if actual_3d_run else 'no'}")
    return {
        "analytic_path": analytic_path,
        "fem_path": fem_path,
        "sorted_path": sorted_path,
        "nearest_path": nearest_path,
        "report_path": report_path,
        "actual_3d_run": actual_3d_run,
    }


def main(argv: list[str] | None = None) -> int:
    run(parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
