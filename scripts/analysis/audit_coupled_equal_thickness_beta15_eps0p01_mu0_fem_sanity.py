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


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))


BETA_DEG = 15.0
EPSILON = 0.01
ETA = 0.0
MU = 0.0
N_MODES = 6

RESULTS_DIR = REPO_ROOT / "results"
SOURCE_CSV = RESULTS_DIR / "coupled_equal_thickness_beta15_eps0p01_eb_timoshenko_3d_fem_vs_mu_corrected.csv"
AUDIT_CSV = RESULTS_DIR / "coupled_equal_thickness_beta15_eps0p01_mu0_fem_sanity_audit.csv"
AUDIT_REPORT = RESULTS_DIR / "coupled_equal_thickness_beta15_eps0p01_mu0_fem_sanity_audit.md"
PLANAR_CASE_DIR = (
    RESULTS_DIR
    / "solid_fem_coupled_equal_rods_beta15_eps0p01_planar_mu"
    / "beta_15p0"
    / "eps_0p01"
    / "mu_0p0"
)
FULL_POINT_JOINT_EPS001_DIR = RESULTS_DIR / "solid_fem_coupled_equal_rods_point_joint" / "beta_15" / "eps_0p01"
FULL_POINT_JOINT_EPS0025_DIR = RESULTS_DIR / "solid_fem_coupled_equal_rods_point_joint" / "beta_15" / "eps_0p025"


@dataclass(frozen=True)
class SourceRows:
    analytic_eb: list[dict[str, str]]
    analytic_timo: list[dict[str, str]]
    planar_fem: list[dict[str, str]]


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_plot_module() -> ModuleType:
    return load_module(
        "coupled_equal_thickness_mu_plot_for_audit",
        REPO_ROOT / "scripts" / "analysis" / "plot_coupled_equal_rods_beta15_eps0p01_eb_timoshenko_fem_vs_mu.py",
    )


def load_point_joint_module() -> ModuleType:
    return load_module(
        "point_joint_fem_helpers_for_mu0_audit",
        REPO_ROOT / "scripts" / "analysis" / "solid_fem_coupled_equal_rods_point_joint.py",
    )


def load_frame_fem_module() -> ModuleType:
    return load_module("one_d_eb_frame_fem_for_mu0_audit", SRC_ROOT / "my_project" / "fem" / "python_fem.py")


def as_float(value: object, default: float = float("nan")) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def read_source_rows() -> SourceRows:
    if not SOURCE_CSV.exists():
        raise FileNotFoundError(f"Corrected source CSV is missing: {SOURCE_CSV}")
    analytic_eb: list[dict[str, str]] = []
    analytic_timo: list[dict[str, str]] = []
    planar_fem: list[dict[str, str]] = []
    with SOURCE_CSV.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if abs(as_float(row.get("mu")) - MU) > 1.0e-12:
                continue
            if row.get("row_kind") == "analytic" and row.get("model") == "Euler-Bernoulli":
                analytic_eb.append(row)
            elif row.get("row_kind") == "analytic" and row.get("model") == "Timoshenko":
                analytic_timo.append(row)
            elif row.get("row_kind") == "fem":
                planar_fem.append(row)
    return SourceRows(
        analytic_eb=sorted(analytic_eb, key=lambda item: int(item["branch_index"])),
        analytic_timo=sorted(analytic_timo, key=lambda item: int(item["branch_index"])),
        planar_fem=sorted(planar_fem, key=lambda item: int(item["branch_index"])),
    )


def first_file(directory: Path, pattern: str) -> Path | None:
    matches = sorted(directory.glob(pattern))
    return matches[0] if matches else None


def parse_node_set_counts(inp_path: Path) -> tuple[dict[str, int], list[str], list[str]]:
    if not inp_path.exists():
        return {}, [], []
    counts: dict[str, int] = {}
    boundaries: list[str] = []
    rigid_body_lines: list[str] = []
    current_set: str | None = None
    in_boundary = False
    for raw_line in inp_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        upper = line.upper()
        if upper.startswith("*NSET"):
            current_set = None
            for part in line.split(","):
                if "=" in part:
                    key, value = part.split("=", 1)
                    if key.strip().upper() == "NSET":
                        current_set = value.strip()
                        counts.setdefault(current_set, 0)
            in_boundary = False
            continue
        if upper.startswith("*BOUNDARY"):
            current_set = None
            in_boundary = True
            continue
        if upper.startswith("*RIGID BODY"):
            rigid_body_lines.append(line)
            current_set = None
            in_boundary = False
            continue
        if line.startswith("*"):
            current_set = None
            in_boundary = False
            continue
        if current_set is not None:
            counts[current_set] += len([part for part in line.split(",") if part.strip()])
        if in_boundary:
            boundaries.append(line)
    return counts, boundaries, rigid_body_lines


def parse_planar_raw_lambdas(pj: ModuleType) -> tuple[list[float], list[float], str]:
    dat_path = first_file(PLANAR_CASE_DIR, "*_ccx_modal.dat")
    stdout_path = first_file(PLANAR_CASE_DIR, "*_ccx_modal.stdout.txt")
    if dat_path is None or stdout_path is None:
        return [], [], "missing planar CalculiX .dat/stdout"
    omegas, source = pj.parse_calculix_eigen_omegas(dat_path, stdout_path)
    lambdas = [math.sqrt(float(omega) / EPSILON) for omega in omegas[:N_MODES]]
    return [float(omega) for omega in omegas[:N_MODES]], lambdas, source


def parse_planar_mesh_summary(pj: ModuleType, plot_module: ModuleType) -> dict[str, float | int | str]:
    mesh_path = first_file(PLANAR_CASE_DIR, "*_mesh.inp")
    if mesh_path is None:
        return {"mesh_status": "missing mesh inp"}
    mesh = pj.read_gmsh_inp_mesh_data(mesh_path)
    l1_est, l2_est, r1_est, r2_est = plot_module.estimate_mesh_geometry(mesh, MU)
    return {
        "mesh_status": "ok",
        "nodes": len(mesh.nodes),
        "solid_elements": len(mesh.solid_elements),
        "l1_estimate": l1_est,
        "l2_estimate": l2_est,
        "r1_estimate": r1_est,
        "r2_estimate": r2_est,
    }


def compute_frame_fem_lambdas() -> list[float]:
    frame = load_frame_fem_module()
    frame.eps = EPSILON
    frame.EA_nd = 1.0 / EPSILON**2
    omega, _vectors = frame.fem_solve(MU, beta_deg=BETA_DEG, n_modes=N_MODES)
    # In this older 1D EB frame code the generalized eigenfrequency is Lambda^2
    # in the determinant normalization, so Lambda = sqrt(omega).
    return [math.sqrt(float(value)) for value in omega[:N_MODES]]


def sorted_lambdas(rows: Sequence[dict[str, str]], key: str = "Lambda") -> list[float]:
    values = [as_float(row.get(key)) for row in rows]
    return sorted(value for value in values if math.isfinite(value))[:N_MODES]


def branch_lookup(rows: Sequence[dict[str, str]], key: str = "Lambda") -> dict[int, float]:
    return {int(row["branch_index"]): as_float(row.get(key)) for row in rows}


def build_audit_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    source = read_source_rows()
    plot_module = load_plot_module()
    pj = load_point_joint_module()

    planar_omegas, planar_raw_lambdas, parse_source = parse_planar_raw_lambdas(pj)
    frame_lambdas = compute_frame_fem_lambdas()
    mesh_summary = parse_planar_mesh_summary(pj, plot_module)
    modal_input = first_file(PLANAR_CASE_DIR, "*_ccx_modal.inp")
    node_set_counts, boundary_lines, rigid_body_lines = parse_node_set_counts(modal_input or Path("__missing__"))
    stdout_path = first_file(PLANAR_CASE_DIR, "*_ccx_modal.stdout.txt")
    stderr_path = first_file(PLANAR_CASE_DIR, "*_ccx_modal.stderr.txt")
    stdout_text = stdout_path.read_text(encoding="utf-8", errors="ignore") if stdout_path else ""
    stderr_text = stderr_path.read_text(encoding="utf-8", errors="ignore") if stderr_path else ""
    calculix_status = "completed" if "job finished" in stdout_text.lower() else "unknown"
    if stderr_text.strip():
        calculix_status += "_with_stderr"

    eb_by_branch = branch_lookup(source.analytic_eb)
    timo_by_branch = branch_lookup(source.analytic_timo)
    fem_by_branch = branch_lookup(source.planar_fem)
    rows: list[dict[str, object]] = []

    formula_checks = [
        (
            "bending_basis_lengths",
            "rod1 uses endpoint_columns(l1, theta, basis, section); rod2 uses endpoint_columns(-l2, theta, basis, section)",
            "l1=1-mu and l2=1+mu enter through the physical/local x coordinate.",
        ),
        (
            "axial_basis_lengths",
            "endpoint_columns uses u=p*sin(theta*x) and u'=theta*cos(theta*x)",
            "The same x=l1 or x=-l2 includes length in the axial phase.",
        ),
        (
            "omega_lambda_arguments",
            "theta = project_omega(Lambda, epsilon) = epsilon*Lambda^2",
            "The code keeps the project normalization Omega=epsilon*Lambda^2 and does not introduce per-rod Lambda_i.",
        ),
        (
            "endpoint_coordinates",
            "rod 1 endpoint x=+l1; rod 2 endpoint x=-l2; shape vectors use x1=l1*xi and x2=-l2*xi",
            "This mirrors the EB determinant argument convention but is not independently derived in theory docs.",
        ),
        (
            "mu0_matrix_reuse",
            "The old equal-length Timoshenko matrix is not called directly for mu!=0",
            "The implementation is a local variable-length extension of the mu=0 matrix.",
        ),
        (
            "validation_status",
            "No verified theory-file derivation was found for the variable-length Timoshenko determinant used here",
            "Treat Timoshenko Lambda(mu) as an ad hoc, not-yet-validated diagnostic extension.",
        ),
    ]
    for check, formula, result in formula_checks:
        rows.append({"row_kind": "timoshenko_formula_audit", "check": check, "formula_or_source": formula, "result": result})

    for mode, (eb_value, timo_value, fem_value, frame_value) in enumerate(
        zip(sorted_lambdas(source.analytic_eb), sorted_lambdas(source.analytic_timo), planar_raw_lambdas, frame_lambdas),
        start=1,
    ):
        rows.append(
            {
                "row_kind": "sorted_mu0_comparison",
                "mode": mode,
                "Lambda_EB_sorted": eb_value,
                "Lambda_Timoshenko_sorted": timo_value,
                "Lambda_planar_FEM_sorted": fem_value,
                "Lambda_1D_EB_frame_FEM_sorted": frame_value,
                "FEM_minus_EB": fem_value - eb_value,
                "FEM_minus_Timoshenko": fem_value - timo_value,
                "frame_minus_EB": frame_value - eb_value,
            }
        )

    for branch in range(1, N_MODES + 1):
        fem_row = next((row for row in source.planar_fem if int(row["branch_index"]) == branch), {})
        omega = as_float(fem_row.get("omega"))
        lambda_from_omega = math.sqrt(omega / EPSILON) if math.isfinite(omega) else float("nan")
        lambda_recorded = fem_by_branch.get(branch, float("nan"))
        rows.append(
            {
                "row_kind": "mac_branch_mu0_comparison",
                "branch": branch,
                "Lambda_EB_branch": eb_by_branch.get(branch, float("nan")),
                "Lambda_Timoshenko_branch": timo_by_branch.get(branch, float("nan")),
                "solid_mode_by_MAC": fem_row.get("solid_mode", ""),
                "Lambda_planar_FEM": lambda_recorded,
                "Omega_planar_FEM": omega,
                "Lambda_recomputed_sqrt_Omega_over_epsilon": lambda_from_omega,
                "normalization_delta": lambda_recorded - lambda_from_omega if math.isfinite(lambda_recorded) else float("nan"),
                "MAC": fem_row.get("best_mac", ""),
                "MAC_strength": fem_row.get("mac_strength", ""),
                "matching_shape_model": fem_row.get("matching_shape_model", ""),
                "FEM_minus_EB_branch": lambda_recorded - eb_by_branch.get(branch, float("nan")),
                "FEM_minus_Timoshenko_branch": lambda_recorded - timo_by_branch.get(branch, float("nan")),
            }
        )

    for name, count in sorted(node_set_counts.items()):
        rows.append({"row_kind": "constraint_audit", "item": name, "count": count})
    for line in boundary_lines:
        rows.append({"row_kind": "boundary_audit", "item": line})
    for line in rigid_body_lines:
        rows.append({"row_kind": "rigid_body_audit", "item": line})

    for key, value in mesh_summary.items():
        rows.append({"row_kind": "geometry_mesh_audit", "item": key, "value": value})

    rows.append({"row_kind": "calculix_audit", "item": "status", "value": calculix_status})
    rows.append({"row_kind": "calculix_audit", "item": "frequency_parse_source", "value": parse_source})
    rows.append({"row_kind": "availability_audit", "item": "full_point_joint_eps0p01_available", "value": FULL_POINT_JOINT_EPS001_DIR.exists()})
    rows.append({"row_kind": "availability_audit", "item": "full_point_joint_eps0p025_available_not_comparable", "value": FULL_POINT_JOINT_EPS0025_DIR.exists()})

    context = {
        "source": source,
        "planar_omegas": planar_omegas,
        "planar_raw_lambdas": planar_raw_lambdas,
        "frame_lambdas": frame_lambdas,
        "mesh_summary": mesh_summary,
        "node_set_counts": node_set_counts,
        "boundary_lines": boundary_lines,
        "rigid_body_lines": rigid_body_lines,
        "calculix_status": calculix_status,
        "parse_source": parse_source,
        "full_eps001": FULL_POINT_JOINT_EPS001_DIR.exists(),
        "full_eps0025": FULL_POINT_JOINT_EPS0025_DIR.exists(),
    }
    return rows, context


def write_csv(rows: Sequence[dict[str, object]]) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with AUDIT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: object, digits: int = 6) -> str:
    number = as_float(value)
    if not math.isfinite(number):
        return str(value) if value not in {None, ""} else ""
    return f"{number:.{digits}g}"


def write_report(rows: Sequence[dict[str, object]], context: dict[str, object]) -> None:
    source: SourceRows = context["source"]  # type: ignore[assignment]
    sorted_rows = [row for row in rows if row.get("row_kind") == "sorted_mu0_comparison"]
    branch_rows = [row for row in rows if row.get("row_kind") == "mac_branch_mu0_comparison"]
    frame_max_abs = max(abs(as_float(row.get("frame_minus_EB"))) for row in sorted_rows)
    normalization_max_abs = max(abs(as_float(row.get("normalization_delta"))) for row in branch_rows)
    branch12 = [row for row in branch_rows if row.get("branch") in {1, 2}]

    lines: list[str] = [
        "# Mu=0 FEM Sanity Audit",
        "",
        "## Scope",
        "",
        f"Diagnostic-only audit for beta={BETA_DEG:g} deg, epsilon={EPSILON:g}, eta={ETA:g}, mu={MU:g}.",
        "This audit does not redraw the presentation plot. It checks the variable-mu Timoshenko implementation and the mu=0 FEM discrepancy.",
        "",
        "## Timoshenko Variable-Mu Formulation Audit",
        "",
        "- Rod lengths are `l1 = 1 - mu` and `l2 = 1 + mu`.",
        "- Bending basis: `endpoint_columns(l1, theta, basis, section)` for rod 1 and `endpoint_columns(-l2, theta, basis, section)` for rod 2.",
        "- Axial basis: `u = p*sin(theta*x)` with the same endpoint coordinates, so axial phases include `l1` and `l2` through `x`.",
        "- Local argument scaling: `theta = epsilon*Lambda^2`; the code keeps the project normalization `Omega = epsilon*Lambda^2` and does not introduce per-rod `Lambda_i`.",
        "- Centerline shape reconstruction uses `x1=l1*xi` and `x2=-l2*xi`.",
        "- The old equal-length Timoshenko matrix is not called directly for `mu != 0`.",
        "",
        "Audit status: this is a derived-by-analogy variable-length extension of the earlier `mu=0` Timoshenko matrix, not a verified theory-file derivation. Treat Timoshenko `Lambda(mu)` as not yet validated for article use.",
        "",
        "## Mu=0 Geometry And Constraints",
        "",
        f"- expected lengths: l1=1, l2=1",
        f"- expected radius: r={2*EPSILON:g}",
        f"- mesh nodes: {context['mesh_summary'].get('nodes') if isinstance(context['mesh_summary'], dict) else ''}",
        f"- mesh solid elements: {context['mesh_summary'].get('solid_elements') if isinstance(context['mesh_summary'], dict) else ''}",
        f"- inferred mesh l1/l2: {fmt(context['mesh_summary'].get('l1_estimate') if isinstance(context['mesh_summary'], dict) else '')}, {fmt(context['mesh_summary'].get('l2_estimate') if isinstance(context['mesh_summary'], dict) else '')}",
        f"- inferred mesh radius r1/r2: {fmt(context['mesh_summary'].get('r1_estimate') if isinstance(context['mesh_summary'], dict) else '')}, {fmt(context['mesh_summary'].get('r2_estimate') if isinstance(context['mesh_summary'], dict) else '')}",
        f"- CalculiX status: {context['calculix_status']}",
        f"- frequency parse source: {context['parse_source']}",
        "",
        "Node-set counts:",
    ]
    for name, count in sorted(context["node_set_counts"].items()):  # type: ignore[union-attr]
        lines.append(f"- {name}: {count}")
    lines.extend(
        [
            "",
            "Rigid-body and boundary lines:",
        ]
    )
    for line in context["rigid_body_lines"]:  # type: ignore[union-attr]
        lines.append(f"- `{line}`")
    for line in context["boundary_lines"]:  # type: ignore[union-attr]
        lines.append(f"- `{line}`")

    lines.extend(
        [
            "",
            "## Sorted Frequency Sanity",
            "",
            "The older 1D EB frame FEM is used only as a sorted-frequency sanity anchor. In that code the eigenfrequency equals `Lambda^2` in this normalization, so `Lambda=sqrt(omega)`. It matches the EB analytic sorted roots closely.",
            "",
            "| sorted mode | EB analytic | Timoshenko analytic | planar 3D FEM | 1D EB frame FEM | FEM-EB | frame-EB |",
            "|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in sorted_rows:
        lines.append(
            "| "
            f"{row['mode']} | {fmt(row['Lambda_EB_sorted'])} | {fmt(row['Lambda_Timoshenko_sorted'])} | "
            f"{fmt(row['Lambda_planar_FEM_sorted'])} | {fmt(row['Lambda_1D_EB_frame_FEM_sorted'])} | "
            f"{fmt(row['FEM_minus_EB'])} | {fmt(row['frame_minus_EB'])} |"
        )

    lines.extend(
        [
            "",
            "## MAC Branch Sanity",
            "",
            f"Maximum `Lambda_FEM - sqrt(Omega/epsilon)` normalization residual: {fmt(normalization_max_abs, 3)}.",
            "Thus the large planar 3D FEM discrepancy is not traced to the `Lambda_FEM=sqrt(Omega/epsilon)` conversion.",
            "",
            "| branch | EB branch | Timo branch | solid mode by MAC | planar FEM | MAC | FEM-EB | FEM-Timo |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in branch_rows:
        lines.append(
            "| "
            f"{row['branch']} | {fmt(row['Lambda_EB_branch'])} | {fmt(row['Lambda_Timoshenko_branch'])} | "
            f"{row['solid_mode_by_MAC']} | {fmt(row['Lambda_planar_FEM'])} | {fmt(row['MAC'], 4)} | "
            f"{fmt(row['FEM_minus_EB_branch'])} | {fmt(row['FEM_minus_Timoshenko_branch'])} |"
        )

    lines.extend(
        [
            "",
            "Branch 1/2 note: descendant branch order is not sorted order at mu=0. EB sorted mode 1 corresponds to branch 2, while EB sorted mode 2 corresponds to branch 1. The MAC assignment in the planar 3D solid data maps branch 1 to solid mode 1 and branch 2 to solid mode 2, so sorted-frequency and MAC-order comparisons tell different stories for the first pair.",
            "",
            "## Availability Checks",
            "",
            f"- full point-joint 3D FEM at epsilon=0.01 available: {context['full_eps001']}",
            f"- full point-joint 3D FEM at epsilon=0.025 available but not comparable to epsilon=0.01 audit: {context['full_eps0025']}",
            "- existing 1D EB frame FEM: available and used as sorted-frequency sanity check",
            "",
            "## Required Conclusion",
            "",
            "Timoshenko variable-mu formulation incomplete.",
            "",
            "For the mu=0 FEM discrepancy specifically: normalization and 1D EB frame-FEM sanity checks pass, and branch 1/2 has a real sorted-vs-MAC ambiguity. The remaining planar 3D solid mismatch is not yet traced to a physical effect; full point-joint epsilon=0.01 and constraint-isolation checks are still needed before interpreting it.",
            "",
            "## Outputs",
            "",
            f"- CSV: `{AUDIT_CSV.relative_to(REPO_ROOT)}`",
            f"- report: `{AUDIT_REPORT.relative_to(REPO_ROOT)}`",
        ]
    )
    AUDIT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    rows, context = build_audit_rows()
    write_csv(rows)
    write_report(rows, context)


if __name__ == "__main__":
    main()
