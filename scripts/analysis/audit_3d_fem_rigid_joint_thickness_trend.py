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
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import assemble_clamped_coupled_matrix  # noqa: E402
from scripts.lib import variable_length_timoshenko as vlt  # noqa: E402


BETA_DEG = 15.0
MU = 0.0
ETA = 0.0
EPSILON_VALUES = [0.005, 0.01, 0.025, 0.05]
N_ANALYTIC_MODES = 10
N_MATCH_MODES = 8
N_SOLID_MODES = 30
MESH_FACTOR = 1.0
N_SLICES_PER_ROD = 80
MAC_STRONG_THRESHOLD = 0.8
MAC_MODERATE_THRESHOLD = 0.5
L_SEGMENT = 1.0
THIN_ROD_LIMIT = 0.1

OUTPUT_DIR = REPO_ROOT / "results"
AUDIT_ROOT = OUTPUT_DIR / "solid_fem_rigid_joint_thickness_trend"
OUTPUT_CSV = OUTPUT_DIR / "3d_fem_rigid_joint_thickness_trend.csv"
OUTPUT_REPORT = OUTPUT_DIR / "3d_fem_rigid_joint_thickness_trend_report.md"
OUTPUT_PNG = OUTPUT_DIR / "3d_fem_rigid_joint_thickness_trend_error_ratio.png"


@dataclass(frozen=True)
class AnalyticMode:
    epsilon: float
    mode: int
    lambda_eb: float
    lambda_timo: float
    omega_eb: float
    omega_timo: float
    omega_cutoff: float
    omega_eb_over_cutoff: float
    omega_timo_over_cutoff: float
    lambda_cutoff: float
    eb_vector: np.ndarray
    timo_vector: np.ndarray
    eb_shape_note: str
    timo_shape_note: str


@dataclass(frozen=True)
class FemCase:
    epsilon: float
    paths: object
    mesh: object | None
    mesh_summary: object | None
    ccx_input: object | None
    calculix_result: object | None
    parse_result: object | None
    mode_metrics: tuple[object, ...]
    solid_shapes: tuple[object, ...]
    gmsh_messages: tuple[str, ...]
    messages: tuple[str, ...]


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
        "point_joint_helpers_for_rigid_joint_thickness_trend",
        REPO_ROOT / "scripts" / "analysis" / "solid_fem_coupled_equal_rods_point_joint.py",
    )


def load_eb_reference_module() -> ModuleType:
    module = load_module(
        "eb_reference_for_rigid_joint_thickness_trend",
        REPO_ROOT / "scripts" / "analysis" / "compare_coupled_equal_rods_eb_timoshenko.py",
    )
    module.BETA_DEG = BETA_DEG
    module.MU = MU
    module.N_MODES = N_ANALYTIC_MODES
    return module


def configure_point_joint_module(pj: ModuleType) -> None:
    pj.N_SOLID_MODES = N_SOLID_MODES
    pj.N_ANALYTIC_MODES = N_ANALYTIC_MODES
    pj.RUN_BETA_DEG_VALUES = [BETA_DEG]
    pj.RUN_EPSILON_VALUES = list(EPSILON_VALUES)
    pj.MESH_SIZE_FACTORS = [MESH_FACTOR]


def safe_float_token(value: float) -> str:
    text = f"{float(value):.6g}"
    if "." not in text and "e" not in text.lower():
        text += ".0"
    return text.replace("-", "m").replace(".", "p")


def case_paths(pj: ModuleType, epsilon: float) -> object:
    eps_token = safe_float_token(float(epsilon))
    beta_token = safe_float_token(BETA_DEG)
    case_dir = AUDIT_ROOT / f"eps_{eps_token}" / f"mesh_{safe_float_token(MESH_FACTOR)}"
    stem = f"rigid_joint_beta{beta_token}_eps{eps_token}"
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


def finite(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def normalize_vector(vector: np.ndarray) -> np.ndarray:
    out = np.asarray(vector, dtype=float)
    norm_value = float(np.linalg.norm(out))
    if norm_value <= 1.0e-28 or not math.isfinite(norm_value):
        return np.full_like(out, np.nan, dtype=float)
    return out / norm_value


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


def slice_coordinates() -> np.ndarray:
    return (np.arange(N_SLICES_PER_ROD, dtype=float) + 0.5) / float(N_SLICES_PER_ROD)


def global_centerline_vector(
    pj: ModuleType,
    u_left: np.ndarray,
    w_left: np.ndarray,
    u_right: np.ndarray,
    w_right: np.ndarray,
) -> np.ndarray:
    geom = pj.geometry_for_beta(BETA_DEG)
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
    return normalize_vector(np.concatenate(parts))


def eb_centerline_vector(pj: ModuleType, epsilon: float, lambda_value: float) -> tuple[np.ndarray, str]:
    matrix = assemble_clamped_coupled_matrix(
        float(lambda_value),
        math.radians(BETA_DEG),
        MU,
        float(epsilon),
    )
    coeff, smallest, ratio = analytic_null_vector(matrix)
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
    note = f"EB normalized determinant null vector smallest_singular={smallest:.3e}, ratio={ratio:.3e}"
    return global_centerline_vector(pj, u_left, w_left, u_right, w_right), note


def timo_centerline_vector(pj: ModuleType, epsilon: float, lambda_value: float) -> tuple[np.ndarray, str]:
    coeff_info = vlt.timo_mode_coefficients(float(lambda_value), BETA_DEG, MU, float(epsilon), ETA)
    factors = vlt.tau_factors(MU, ETA)
    section1 = vlt.section_from_epsilon_tau(float(epsilon), factors.tau1)
    section2 = vlt.section_from_epsilon_tau(float(epsilon), factors.tau2)
    basis1 = vlt.timo_basis(float(lambda_value), float(epsilon), section1)
    basis2 = vlt.timo_basis(float(lambda_value), float(epsilon), section2)
    xi = slice_coordinates()
    x1 = xi
    x2 = -xi
    fields1 = vlt.evaluate_timo_rod_fields(
        float(lambda_value),
        float(epsilon),
        section1,
        basis1,
        coeff_info.coeff[0:3],
        x1,
    )
    fields2 = vlt.evaluate_timo_rod_fields(
        float(lambda_value),
        float(epsilon),
        section2,
        basis2,
        coeff_info.coeff[3:6],
        x2,
    )
    warnings = tuple(dict.fromkeys([*coeff_info.warnings, *basis1.warnings, *basis2.warnings]))
    warning_text = "; ".join(warnings) if warnings else "no Timoshenko basis warnings"
    note = (
        "Timoshenko variable-length helper eta=0 "
        f"smallest_singular={coeff_info.smallest_singular_value:.3e}, "
        f"ratio={coeff_info.singular_value_ratio:.3e}; {warning_text}"
    )
    return (
        global_centerline_vector(
            pj,
            np.asarray(fields1["u"], dtype=float),
            np.asarray(fields1["w"], dtype=float),
            np.asarray(fields2["u"], dtype=float),
            np.asarray(fields2["w"], dtype=float),
        ),
        note,
    )


def compute_analytic_modes(pj: ModuleType, eb_reference: ModuleType) -> tuple[list[AnalyticMode], list[str]]:
    rows: list[AnalyticMode] = []
    warnings: list[str] = []
    for epsilon in EPSILON_VALUES:
        eb_roots = eb_reference.eb_roots(float(epsilon))
        timo_roots, timo_warnings = vlt.timo_sorted_roots(
            BETA_DEG,
            MU,
            float(epsilon),
            N_ANALYTIC_MODES,
            eta=ETA,
        )
        warnings.extend(f"epsilon={float(epsilon):g}: {item}" for item in timo_warnings)
        section = vlt.section_from_epsilon(float(epsilon))
        omega_cutoff = vlt.omega_cutoff(section)
        lambda_cutoff = vlt.lambda_cutoff(float(epsilon), section)
        for mode_index in range(1, N_ANALYTIC_MODES + 1):
            lambda_eb = float(eb_roots[mode_index - 1])
            lambda_timo = float(timo_roots[mode_index - 1])
            omega_eb = vlt.project_omega(lambda_eb, float(epsilon))
            omega_timo = vlt.project_omega(lambda_timo, float(epsilon))
            eb_vector, eb_note = eb_centerline_vector(pj, float(epsilon), lambda_eb)
            try:
                timo_vector, timo_note = timo_centerline_vector(pj, float(epsilon), lambda_timo)
            except Exception as exc:
                warnings.append(f"epsilon={float(epsilon):g}, mode={mode_index}: Timoshenko shape failed: {exc}")
                timo_vector = np.full_like(eb_vector, np.nan, dtype=float)
                timo_note = f"Timoshenko shape failed: {exc}"
            rows.append(
                AnalyticMode(
                    epsilon=float(epsilon),
                    mode=mode_index,
                    lambda_eb=lambda_eb,
                    lambda_timo=lambda_timo,
                    omega_eb=omega_eb,
                    omega_timo=omega_timo,
                    omega_cutoff=omega_cutoff,
                    omega_eb_over_cutoff=omega_eb / omega_cutoff,
                    omega_timo_over_cutoff=omega_timo / omega_cutoff,
                    lambda_cutoff=lambda_cutoff,
                    eb_vector=eb_vector,
                    timo_vector=timo_vector,
                    eb_shape_note=eb_note,
                    timo_shape_note=timo_note,
                )
            )
    return rows, sorted(set(warnings))


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


def parse_input_constraints(path: Path) -> dict[str, object]:
    node_sets: dict[str, set[int]] = {}
    rigid_lines: list[str] = []
    boundaries: list[str] = []
    current_set: str | None = None
    in_boundary = False
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
        if upper.startswith("*RIGID BODY"):
            rigid_lines.append(line)
            current_set = None
            in_boundary = False
            continue
        if upper.startswith("*BOUNDARY"):
            current_set = None
            in_boundary = True
            continue
        if line.startswith("*"):
            current_set = None
            in_boundary = False
            continue
        if current_set:
            for part in line.split(","):
                token = part.strip()
                if not token:
                    continue
                try:
                    node_sets[current_set].add(int(token))
                except ValueError:
                    pass
        if in_boundary:
            boundaries.append(line)
    dependent = node_sets.get("JOINT_COUPLED_ALL", set())
    conflicts = []
    reference_restrictions = []
    for line in boundaries:
        set_name = line.split(",", 1)[0].strip()
        if set_name.upper() == "JOINT_REF":
            reference_restrictions.append(line)
        overlap = node_sets.get(set_name, set()) & dependent
        if overlap:
            conflicts.append(f"{set_name}:{len(overlap)}")
    return {
        "node_sets": {key: len(value) for key, value in sorted(node_sets.items())},
        "rigid_lines": rigid_lines,
        "boundaries": boundaries,
        "reference_restrictions": reference_restrictions,
        "dependent_node_spc_conflict": ";".join(conflicts) if conflicts else "none",
    }


def mesh_inferred_geometry(pj: ModuleType, epsilon: float, mesh: object) -> dict[str, float]:
    geom = pj.geometry_for_beta(BETA_DEG)
    memberships, _counts = pj.rod_node_memberships(BETA_DEG, mesh)
    out: dict[str, float] = {}
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
    out["target_radius"] = pj.radius_from_epsilon(float(epsilon))
    out["target_diameter_to_segment_length"] = pj.diameter_to_segment_length(float(epsilon))
    return out


def run_fem_case(pj: ModuleType, epsilon: float) -> FemCase:
    messages: list[str] = []
    paths = case_paths(pj, float(epsilon))
    gmsh_exe, gmsh_source = pj.resolve_gmsh_executable()
    ccx_exe, ccx_source = pj.resolve_ccx_executable()
    if gmsh_exe is None:
        raise RuntimeError("Gmsh executable not found. Set GMSH_EXE.")
    if ccx_exe is None:
        raise RuntimeError("CalculiX executable not found. Set CCX_EXE.")
    mesh_ok, mesh_message, _geometry_label, gmsh_messages = pj.generate_mesh(
        BETA_DEG,
        float(epsilon),
        paths,
        gmsh_exe,
        MESH_FACTOR,
    )
    messages.append(
        f"epsilon={float(epsilon):g}: {mesh_message}; Gmsh={gmsh_source}; CalculiX={ccx_source}"
    )
    if not mesh_ok:
        return FemCase(
            epsilon=float(epsilon),
            paths=paths,
            mesh=None,
            mesh_summary=None,
            ccx_input=None,
            calculix_result=None,
            parse_result=None,
            mode_metrics=(),
            solid_shapes=(),
            gmsh_messages=tuple(gmsh_messages),
            messages=tuple(messages),
        )
    mesh = pj.read_gmsh_inp_mesh_data(paths.gmsh_inp)
    mesh_summary = pj.mesh_summary(BETA_DEG, float(epsilon), mesh, gmsh_messages)
    ccx_input = pj.write_calculix_inputs(
        paths,
        BETA_DEG,
        float(epsilon),
        mesh,
        mesh_size_factor=MESH_FACTOR,
        planar_constraint=False,
    )
    calculix_result = pj.run_calculix(paths, ccx_exe, BETA_DEG, float(epsilon))
    mode_metrics, parse_result = pj.compute_mode_metrics(
        BETA_DEG,
        float(epsilon),
        paths,
        mesh,
        calculix_result.parsed_omegas,
    )
    solid_shapes = (
        pj.compute_solid_centerline_shapes(
            BETA_DEG,
            float(epsilon),
            mesh,
            calculix_result.parsed_omegas,
            parse_result,
            mode_metrics,
        )
        if calculix_result.success and parse_result.success
        else []
    )
    return FemCase(
        epsilon=float(epsilon),
        paths=paths,
        mesh=mesh,
        mesh_summary=mesh_summary,
        ccx_input=ccx_input,
        calculix_result=calculix_result,
        parse_result=parse_result,
        mode_metrics=tuple(mode_metrics),
        solid_shapes=tuple(solid_shapes),
        gmsh_messages=tuple(gmsh_messages),
        messages=tuple(messages),
    )


def candidate_solid_shapes(solid_shapes: Sequence[object]) -> list[object]:
    preferred = [
        shape
        for shape in solid_shapes
        if getattr(shape, "classification", "") in {"in_plane_bending_like", "mixed_or_unclassified"}
    ]
    return preferred or list(solid_shapes)


def best_solid_by_mac(solid_shapes: Sequence[object], analytic_vector: np.ndarray) -> tuple[object | None, float]:
    best_shape = None
    best_mac = float("nan")
    for shape in candidate_solid_shapes(solid_shapes):
        mac = centerline_mac(getattr(shape, "vector"), analytic_vector)
        if best_shape is None or (math.isfinite(mac) and mac > best_mac):
            best_shape = shape
            best_mac = mac
    return best_shape, best_mac


def match_analytic_modes(
    analytic_modes: Sequence[AnalyticMode],
    fem_cases: Sequence[FemCase],
) -> list[dict[str, object]]:
    by_epsilon = {case.epsilon: case for case in fem_cases}
    rows: list[dict[str, object]] = []
    for analytic in analytic_modes:
        if analytic.mode > N_MATCH_MODES:
            continue
        case = by_epsilon.get(analytic.epsilon)
        if case is None or not case.solid_shapes:
            rows.append(
                {
                    "row_kind": "mac_match",
                    "epsilon": analytic.epsilon,
                    "analytic_mode": analytic.mode,
                    "status": "missing_fem_shapes",
                    "best_MAC_strength": "missing",
                }
            )
            continue
        eb_shape, eb_mac = best_solid_by_mac(case.solid_shapes, analytic.eb_vector)
        timo_shape, timo_mac = best_solid_by_mac(case.solid_shapes, analytic.timo_vector)
        selected_model = "Timoshenko" if finite(timo_mac) > finite(eb_mac) else "EB"
        selected_shape = timo_shape if selected_model == "Timoshenko" else eb_shape
        selected_mac = timo_mac if selected_model == "Timoshenko" else eb_mac
        if selected_shape is None:
            rows.append(
                {
                    "row_kind": "mac_match",
                    "epsilon": analytic.epsilon,
                    "analytic_mode": analytic.mode,
                    "status": "no_mac_candidate",
                    "best_MAC_strength": "missing",
                }
            )
            continue
        lambda_solid = float(getattr(selected_shape, "lambda_fem"))
        rel_eb = abs(lambda_solid - analytic.lambda_eb) / abs(analytic.lambda_eb)
        rel_timo = abs(lambda_solid - analytic.lambda_timo) / abs(analytic.lambda_timo)
        eb_mode = int(getattr(eb_shape, "solid_mode")) if eb_shape is not None else ""
        timo_mode = int(getattr(timo_shape, "solid_mode")) if timo_shape is not None else ""
        rows.append(
            {
                "row_kind": "mac_match",
                "epsilon": analytic.epsilon,
                "analytic_mode": analytic.mode,
                "Lambda_EB": analytic.lambda_eb,
                "Lambda_Timoshenko": analytic.lambda_timo,
                "Omega_EB_over_cutoff": analytic.omega_eb_over_cutoff,
                "Omega_Timoshenko_over_cutoff": analytic.omega_timo_over_cutoff,
                "Lambda_cutoff": analytic.lambda_cutoff,
                "best_EB_solid_mode": eb_mode,
                "best_EB_MAC": eb_mac,
                "best_EB_class": getattr(eb_shape, "classification", "") if eb_shape is not None else "",
                "best_Timoshenko_solid_mode": timo_mode,
                "best_Timoshenko_MAC": timo_mac,
                "best_Timoshenko_class": getattr(timo_shape, "classification", "") if timo_shape is not None else "",
                "selected_shape_model_by_MAC": selected_model,
                "selected_solid_mode": int(getattr(selected_shape, "solid_mode")),
                "selected_solid_class": getattr(selected_shape, "classification", ""),
                "selected_MAC": selected_mac,
                "best_MAC_strength": mac_strength(selected_mac),
                "Lambda_FEM": lambda_solid,
                "Omega_FEM": float(getattr(selected_shape, "omega")),
                "rel_error_EB": rel_eb,
                "rel_error_Timoshenko": rel_timo,
                "closer_model_after_MAC": "Timoshenko" if rel_timo < rel_eb else "EB",
                "different_best_solid_modes": eb_mode != timo_mode,
                "ambiguity_note": "EB and Timoshenko shapes choose different FEM modes" if eb_mode != timo_mode else "",
                "status": "matched",
            }
        )
    return rows


def summarize_matches(match_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for epsilon in EPSILON_VALUES:
        subset = [row for row in match_rows if abs(finite(row.get("epsilon")) - float(epsilon)) <= 1.0e-12]
        eligible = [row for row in subset if row.get("best_MAC_strength") in {"strong", "moderate"}]
        weak = [row for row in subset if row.get("best_MAC_strength") == "weak"]
        missing = [row for row in subset if row.get("best_MAC_strength") == "missing"]
        rel_eb = [finite(row.get("rel_error_EB")) for row in eligible if math.isfinite(finite(row.get("rel_error_EB")))]
        rel_timo = [
            finite(row.get("rel_error_Timoshenko"))
            for row in eligible
            if math.isfinite(finite(row.get("rel_error_Timoshenko")))
        ]
        closer_eb = sum(1 for row in eligible if row.get("closer_model_after_MAC") == "EB")
        closer_timo = sum(1 for row in eligible if row.get("closer_model_after_MAC") == "Timoshenko")
        mean_eb = float(np.mean(rel_eb)) if rel_eb else float("nan")
        mean_timo = float(np.mean(rel_timo)) if rel_timo else float("nan")
        ratio = mean_timo / mean_eb if math.isfinite(mean_eb) and abs(mean_eb) > 1.0e-14 else float("nan")
        rows.append(
            {
                "row_kind": "epsilon_summary",
                "epsilon": float(epsilon),
                "diameter_to_segment_length": 4.0 * float(epsilon),
                "thin_rod_valid": 4.0 * float(epsilon) <= THIN_ROD_LIMIT + 1.0e-15,
                "strong_count": sum(1 for row in subset if row.get("best_MAC_strength") == "strong"),
                "moderate_count": sum(1 for row in subset if row.get("best_MAC_strength") == "moderate"),
                "weak_count": len(weak),
                "missing_count": len(missing),
                "eligible_count": len(eligible),
                "mean_rel_error_EB": mean_eb,
                "mean_rel_error_Timoshenko": mean_timo,
                "mean_error_Timo_over_EB": ratio,
                "max_rel_error_EB": max(rel_eb) if rel_eb else float("nan"),
                "max_rel_error_Timoshenko": max(rel_timo) if rel_timo else float("nan"),
                "closer_EB_count": closer_eb,
                "closer_Timoshenko_count": closer_timo,
                "different_best_solid_mode_count": sum(
                    1 for row in eligible if str(row.get("different_best_solid_modes")) in {"True", "true", "1"}
                ),
            }
        )
    return rows


def case_rows(pj: ModuleType, cases: Sequence[FemCase]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in cases:
        prefix = {
            "epsilon": case.epsilon,
            "beta_deg": BETA_DEG,
            "mu": MU,
            "eta": ETA,
            "mesh_factor": MESH_FACTOR,
        }
        result = case.calculix_result
        parse_result = case.parse_result
        rows.append(
            {
                **prefix,
                "row_kind": "case_summary",
                "case_dir": str(case.paths.case_dir.relative_to(REPO_ROOT)),
                "calculix_success": bool(result.success) if result is not None else False,
                "calculix_message": result.message if result is not None else "not run",
                "parsed_frequency_count": len(result.parsed_omegas) if result is not None else 0,
                "parse_source": result.parse_source if result is not None else "",
                "frd_parse_success": bool(parse_result.success) if parse_result is not None else False,
                "frd_parse_message": parse_result.message if parse_result is not None else "",
            }
        )
        if case.mesh is not None:
            components, largest = pj.connected_component_counts(case.mesh)
            mesh_warning = any(
                "negative jacobian" in item.lower() or "distort" in item.lower() for item in case.gmsh_messages
            )
            rows.append(
                {
                    **prefix,
                    "row_kind": "mesh_audit",
                    "nodes": len(case.mesh.nodes),
                    "solid_elements": len(case.mesh.solid_elements),
                    "connected_components_before_constraints": components,
                    "largest_component_elements": largest,
                    "negative_jacobian_or_distortion_warning": mesh_warning,
                    "gmsh_warning_count": len(case.gmsh_messages),
                    "bbox_min": getattr(case.mesh, "bbox_min", ""),
                    "bbox_max": getattr(case.mesh, "bbox_max", ""),
                    **mesh_inferred_geometry(pj, case.epsilon, case.mesh),
                }
            )
        if case.ccx_input is not None:
            constraint = parse_input_constraints(case.paths.ccx_modal_inp)
            rows.append(
                {
                    **prefix,
                    "row_kind": "constraint_audit",
                    "rod1_fixed_node_count": case.ccx_input.rod1_fixed_node_count,
                    "rod2_fixed_node_count": case.ccx_input.rod2_fixed_node_count,
                    "rod1_inner_coupled_node_count": case.ccx_input.rod1_inner_node_count,
                    "rod2_inner_coupled_node_count": case.ccx_input.rod2_inner_node_count,
                    "joint_coupled_node_count": case.ccx_input.joint_coupled_node_count,
                    "joint_ref_node_id": case.ccx_input.joint_ref_node_id,
                    "rigid_body_lines": " | ".join(constraint["rigid_lines"]),
                    "boundary_lines": " | ".join(constraint["boundaries"]),
                    "reference_restrictions": " | ".join(constraint["reference_restrictions"]),
                    "dependent_node_spc_conflict": constraint["dependent_node_spc_conflict"],
                }
            )
            for name, count in constraint["node_sets"].items():
                rows.append({**prefix, "row_kind": "node_set_count", "node_set": name, "count": count})
        if result is not None:
            warnings = warning_lines_from_files(case.paths)
            rows.append(
                {
                    **prefix,
                    "row_kind": "calculix_warning_audit",
                    "warning_count": len(warnings),
                    "warnings": " | ".join(warnings[:20]),
                }
            )
        for metric in case.mode_metrics:
            rows.append(
                {
                    **prefix,
                    "row_kind": "mode_metric",
                    "solid_mode": metric.solid_mode,
                    "Omega_FEM": metric.omega,
                    "Lambda_FEM": metric.lambda_fem,
                    "classification": metric.classification,
                    "classification_note": metric.classification_note,
                    "in_plane_energy_fraction": metric.in_plane_energy_fraction,
                    "out_of_plane_energy_fraction": metric.out_of_plane_energy_fraction,
                    "axial_energy_fraction": metric.axial_energy_fraction,
                    "mean_in_plane_bending_fraction": metric.mean_in_plane_bending_fraction,
                    "mean_out_of_plane_bending_fraction": metric.mean_out_of_plane_bending_fraction,
                    "torsion_indicator": metric.torsion_indicator,
                    "joint_motion_fraction": metric.joint_motion_fraction,
                    "warping_indicator": metric.warping_indicator,
                }
            )
    return rows


def sorted_comparison_rows(analytic_modes: Sequence[AnalyticMode], cases: Sequence[FemCase]) -> list[dict[str, object]]:
    cases_by_epsilon = {case.epsilon: case for case in cases}
    rows: list[dict[str, object]] = []
    for analytic in analytic_modes:
        case = cases_by_epsilon.get(analytic.epsilon)
        omega = ""
        lambda_fem = ""
        if case is not None and case.calculix_result is not None:
            parsed = list(case.calculix_result.parsed_omegas)
            if analytic.mode <= len(parsed):
                omega = float(parsed[analytic.mode - 1])
                lambda_fem = math.sqrt(omega / analytic.epsilon)
        rows.append(
            {
                "row_kind": "sorted_comparison",
                "epsilon": analytic.epsilon,
                "analytic_mode": analytic.mode,
                "solid_sorted_mode": analytic.mode if lambda_fem != "" else "",
                "Omega_FEM": omega,
                "Lambda_FEM": lambda_fem,
                "Lambda_EB": analytic.lambda_eb,
                "Lambda_Timoshenko": analytic.lambda_timo,
                "rel_error_EB": abs(float(lambda_fem) - analytic.lambda_eb) / abs(analytic.lambda_eb)
                if lambda_fem != ""
                else "",
                "rel_error_Timoshenko": abs(float(lambda_fem) - analytic.lambda_timo) / abs(analytic.lambda_timo)
                if lambda_fem != ""
                else "",
                "notes": "sorted-frequency comparison only; MAC rows are the main identity evidence",
            }
        )
    return rows


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


def markdown_table(rows: Sequence[dict[str, object]], columns: Sequence[tuple[str, str]]) -> list[str]:
    lines = ["| " + " | ".join(label for label, _key in columns) + " |"]
    lines.append("| " + " | ".join("---" for _ in columns) + " |")
    for row in rows:
        values = []
        for _label, key in columns:
            value = row.get(key, "")
            values.append(fmt(value, 5) if isinstance(value, (int, float, np.floating)) else str(value))
        lines.append("| " + " | ".join(values) + " |")
    return lines


def trend_conclusion(summary_rows: Sequence[dict[str, object]]) -> str:
    eligible = [row for row in summary_rows if int(row.get("eligible_count", 0)) > 0]
    ratios = [(finite(row.get("epsilon")), finite(row.get("mean_error_Timo_over_EB"))) for row in eligible]
    ratios = [(eps, ratio) for eps, ratio in ratios if math.isfinite(eps) and math.isfinite(ratio)]
    if not ratios:
        return "No eligible strong/moderate MAC rows were available for the trend conclusion."
    all_timo_better = all(ratio < 1.0 for _eps, ratio in ratios)
    sorted_ratios = sorted(ratios)
    thick_half = [ratio for eps, ratio in sorted_ratios if eps >= 0.025]
    thin_half = [ratio for eps, ratio in sorted_ratios if eps <= 0.01]
    thick_mean = float(np.mean(thick_half)) if thick_half else float("nan")
    thin_mean = float(np.mean(thin_half)) if thin_half else float("nan")
    monotone = all(right <= left + 1.0e-12 for (_e1, left), (_e2, right) in zip(sorted_ratios, sorted_ratios[1:]))
    last_eps, last_ratio = sorted_ratios[-1]
    first_eps, first_ratio = sorted_ratios[0]
    if monotone and last_ratio < 1.0 and math.isfinite(thick_mean) and math.isfinite(thin_mean) and thick_mean < thin_mean:
        return (
            "The mean Timoshenko/EB error ratio decreases over the completed epsilon grid "
            f"({first_ratio:.3g} at epsilon={first_eps:g} to {last_ratio:.3g} at epsilon={last_eps:g}). "
            "The thinnest case is essentially parity/slightly EB-favored, while the thicker cases become "
            "Timoshenko-favored, most clearly at epsilon=0.05. This supports the comparative thickness trend "
            "for the full rigid end-face engineering joint, without claiming exact 1D/3D joint equivalence."
        )
    if all_timo_better and math.isfinite(thick_mean) and math.isfinite(thin_mean) and thick_mean < thin_mean:
        return (
            "For the eligible MAC-matched rows, Timoshenko is closer than EB at every completed epsilon "
            "and the mean Timoshenko/EB error ratio is lower in the thicker epsilon >= 0.025 cases than "
            "in the thinner epsilon <= 0.01 cases. This supports the comparative thickness trend, while "
            "not proving exact 1D/3D joint equivalence."
        )
    if all_timo_better:
        qualifier = "monotonically decreases" if monotone else "does not monotonically decrease"
        return (
            "For the eligible MAC-matched rows, Timoshenko is closer than EB at every completed epsilon, "
            f"but the mean Timoshenko/EB error ratio {qualifier} over the tested grid."
        )
    return (
        "The completed eligible rows do not show Timoshenko closer than EB for every epsilon; inspect the "
        "per-mode MAC and frequency rows before drawing a thickness-trend conclusion."
    )


def write_plot(summary_rows: Sequence[dict[str, object]]) -> None:
    plot_rows = [
        row
        for row in summary_rows
        if math.isfinite(finite(row.get("epsilon"))) and math.isfinite(finite(row.get("mean_error_Timo_over_EB")))
    ]
    if not plot_rows:
        return
    plot_rows.sort(key=lambda row: finite(row.get("epsilon")))
    eps = [finite(row.get("epsilon")) for row in plot_rows]
    ratios = [finite(row.get("mean_error_Timo_over_EB")) for row in plot_rows]
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    ax.plot(eps, ratios, marker="o", linewidth=1.8, color="#1f77b4")
    ax.axhline(1.0, color="0.35", linewidth=1.0, linestyle="--")
    ax.set_xlabel("epsilon")
    ax.set_ylabel("mean relative error ratio, Timoshenko / EB")
    ax.set_title("Rigid end-face 3D FEM trend")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=180)
    plt.close(fig)


def write_report(rows: Sequence[dict[str, object]], messages: Sequence[str]) -> None:
    summary_rows = [row for row in rows if row.get("row_kind") == "epsilon_summary"]
    match_rows = [row for row in rows if row.get("row_kind") == "mac_match"]
    case_summary_rows = [row for row in rows if row.get("row_kind") == "case_summary"]
    mesh_rows = [row for row in rows if row.get("row_kind") == "mesh_audit"]
    constraint_rows = [row for row in rows if row.get("row_kind") == "constraint_audit"]
    warning_rows = [row for row in rows if row.get("row_kind") == "calculix_warning_audit"]
    warning_rows = [row for row in warning_rows if finite(row.get("warning_count")) > 0]
    lines = [
        "# 3D FEM Rigid Joint Thickness Trend Audit",
        "",
        "## Purpose",
        "",
        "Diagnostic-only comparative validation of a clear physical 3D engineering joint: two separate solid cylinders, outer ends clamped, and all inner end-face nodes coupled to one common CalculiX `JOINT_REF` using `*RIGID BODY`.",
        "",
        "This 3D joint model is not tuned to match the analytic model. It uses no planar constraint, no central patch, no patch-radius calibration, and no fused joint volume. Exact agreement with the 1D analytic point-joint is not expected. The validation question here is comparative: as rod thickness increases, does the Timoshenko analytic reference become closer to this 3D FEM model than Euler-Bernoulli?",
        "",
        "This is diagnostic evidence only until mesh convergence, unique-assignment policy, and visual mode-shape review are complete. No article files, `paper_dorofeev_style` files, article figures, old determinants, `src/my_project/analytic/formulas.py`, old solvers, existing FEM physical-model baselines, baseline results, or analytic Timoshenko helper formulas are modified.",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- mu: {MU:g}",
        f"- eta: {ETA:g}",
        f"- epsilon values: {', '.join(f'{value:g}' for value in EPSILON_VALUES)}",
        f"- analytic modes solved: {N_ANALYTIC_MODES}; MAC trend uses first {N_MATCH_MODES}",
        f"- solid modes requested: {N_SOLID_MODES}",
        f"- mesh factor: {MESH_FACTOR:g}",
        "",
        "## Trend Summary",
        "",
        *markdown_table(
            summary_rows,
            [
                ("epsilon", "epsilon"),
                ("D/L", "diameter_to_segment_length"),
                ("thin EB", "thin_rod_valid"),
                ("strong", "strong_count"),
                ("moderate", "moderate_count"),
                ("weak", "weak_count"),
                ("eligible", "eligible_count"),
                ("mean EB", "mean_rel_error_EB"),
                ("mean Timo", "mean_rel_error_Timoshenko"),
                ("Timo/EB", "mean_error_Timo_over_EB"),
                ("closer EB", "closer_EB_count"),
                ("closer Timo", "closer_Timoshenko_count"),
                ("ambiguous", "different_best_solid_mode_count"),
            ],
        ),
        "",
        trend_conclusion(summary_rows),
        "",
        "## Per-Mode MAC Matches",
        "",
        *markdown_table(
            match_rows,
            [
                ("eps", "epsilon"),
                ("mode", "analytic_mode"),
                ("strength", "best_MAC_strength"),
                ("sel model", "selected_shape_model_by_MAC"),
                ("solid", "selected_solid_mode"),
                ("class", "selected_solid_class"),
                ("MAC", "selected_MAC"),
                ("rel EB", "rel_error_EB"),
                ("rel Timo", "rel_error_Timoshenko"),
                ("closer", "closer_model_after_MAC"),
                ("diff modes", "different_best_solid_modes"),
            ],
        ),
        "",
        "Weak rows are recorded but excluded from the main trend conclusion.",
        "",
        "## Solver Status",
        "",
        *markdown_table(
            case_summary_rows,
            [
                ("epsilon", "epsilon"),
                ("success", "calculix_success"),
                ("freqs", "parsed_frequency_count"),
                ("FRD", "frd_parse_success"),
                ("message", "calculix_message"),
            ],
        ),
        "",
        "## Mesh Audit",
        "",
        *markdown_table(
            mesh_rows,
            [
                ("epsilon", "epsilon"),
                ("nodes", "nodes"),
                ("elems", "solid_elements"),
                ("components", "connected_components_before_constraints"),
                ("target r", "target_radius"),
                ("rod1 L", "rod1_length_estimate"),
                ("rod2 L", "rod2_length_estimate"),
                ("rod1 r", "rod1_radius_estimate"),
                ("rod2 r", "rod2_radius_estimate"),
                ("gmsh warn", "gmsh_warning_count"),
                ("bad jac", "negative_jacobian_or_distortion_warning"),
            ],
        ),
        "",
        "## Constraint Audit",
        "",
        *markdown_table(
            constraint_rows,
            [
                ("epsilon", "epsilon"),
                ("r1 fixed", "rod1_fixed_node_count"),
                ("r2 fixed", "rod2_fixed_node_count"),
                ("r1 inner", "rod1_inner_coupled_node_count"),
                ("r2 inner", "rod2_inner_coupled_node_count"),
                ("joint", "joint_coupled_node_count"),
                ("rigid body", "rigid_body_lines"),
                ("SPC conflict", "dependent_node_spc_conflict"),
                ("ref restrictions", "reference_restrictions"),
            ],
        ),
        "",
        "## CalculiX Warnings",
        "",
    ]
    if warning_rows:
        lines.extend(
            f"- epsilon={fmt(row.get('epsilon'))}: count={row.get('warning_count')}; {row.get('warnings', '')}"
            for row in warning_rows
        )
    else:
        lines.append("- none recorded")
    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
            f"- report: `{OUTPUT_REPORT.relative_to(REPO_ROOT)}`",
            f"- plot: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
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
    configure_point_joint_module(pj)
    eb_reference = load_eb_reference_module()
    analytic_modes, analytic_warnings = compute_analytic_modes(pj, eb_reference)
    cases: list[FemCase] = []
    messages: list[str] = list(analytic_warnings)
    for epsilon in EPSILON_VALUES:
        case = run_fem_case(pj, float(epsilon))
        cases.append(case)
        messages.extend(case.messages)
    rows: list[dict[str, object]] = []
    rows.extend(case_rows(pj, cases))
    rows.extend(sorted_comparison_rows(analytic_modes, cases))
    match_rows = match_analytic_modes(analytic_modes, cases)
    rows.extend(match_rows)
    summary_rows = summarize_matches(match_rows)
    rows.extend(summary_rows)
    write_csv(rows)
    write_plot(summary_rows)
    write_report(rows, messages)
    print("3D FEM rigid joint thickness trend audit")
    for row in summary_rows:
        print(
            "epsilon={epsilon:g}: eligible={eligible_count}, mean_EB={mean_rel_error_EB:.6g}, "
            "mean_Timo={mean_rel_error_Timoshenko:.6g}, ratio={mean_error_Timo_over_EB:.6g}, "
            "closer Timo/EB={closer_Timoshenko_count}/{closer_EB_count}".format(**row)
        )
    print(f"Wrote {OUTPUT_CSV.relative_to(REPO_ROOT)}")
    print(f"Wrote {OUTPUT_REPORT.relative_to(REPO_ROOT)}")
    print(f"Wrote {OUTPUT_PNG.relative_to(REPO_ROOT)}")
    print(f"Wrote generated FEM files under {AUDIT_ROOT.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
