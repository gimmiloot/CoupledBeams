"""Audit readiness for the design-only limited Chapter-2 3D FEM anchor.

This entry point inventories material data and locally available solver stacks,
then writes frozen design contracts.  It deliberately contains no mesh builder
and no 3D eigenfrequency calculation.
"""

from __future__ import annotations

import argparse
import csv
import importlib.metadata
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib.yartsev_ch2_monoclinic_rod import hms_dx_209_material


DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_limited_3d_fem_anchor"
    / "design_readiness"
)
SOURCE_FILE = "docs/literature/pdf/Глава 1_compressed.pdf"
SOURCE_TABLE_PAGE = "45"
SOURCE_TABLE = "Table 1.1"
SOURCE_CONVENTION_PAGE = "40"
SOURCE_CONVENTION_EQUATION = "equation (1.50)"
POISSON_CONVENTION = (
    "nu_ij=-epsilon_j/epsilon_i under uniaxial sigma_i; "
    "S_ij=-nu_ij/E_i"
)

REQUIRED_ENGINEERING_CONSTANTS = (
    "E1",
    "E2",
    "E3",
    "G12",
    "G13",
    "G23",
    "nu12",
    "nu13",
    "nu23",
    "rho",
)
ELASTIC_MODULI = ("E1", "E2", "E3", "G12", "G13", "G23")
POISSON_RATIOS = ("nu12", "nu13", "nu23")
PROVENANCE_FIELDS = (
    "material_name",
    "units",
    "source_file",
    "source_printed_page",
    "source_table_or_equation",
    "poisson_convention",
    "source_audit_status",
)

STATUS_READY = "READY_FOR_3D_A1"
STATUS_MATERIAL = "BLOCKED_INCOMPLETE_3D_CONSTITUTIVE_DATA"
STATUS_SOLVER = "BLOCKED_NO_CAPABLE_3D_SOLVER"
STATUS_BOTH = "BLOCKED_MATERIAL_AND_SOLVER"
STATUS_INCONCLUSIVE = "INCONCLUSIVE_READINESS_AUDIT"
TENSOR_NUMERICAL_TOLERANCE = 1.0e-12

ONE_D_STATUSES = {
    "phase": "COMPLETE",
    "overall": "PARTIAL_PASS",
    "UT-0": "PASS",
    "UT-1": "PARTIAL_PASS",
    "UT-1a": "INCONCLUSIVE_NUMERICAL_AUDIT",
    "UT-2": "PASS",
    "UT-3": "PASS",
}

FUTURE_THRESHOLDS = {
    "orientation_orthogonality_residual": 1.0e-12,
    "relative_mass_difference": 1.0e-10,
    "target_parity_maximum": -0.99,
    "target_antisymmetric_residual": 1.0e-3,
    "cluster_relative_gap": 1.0e-3,
    "mesh_convergence_modes_1_3": 5.0e-3,
    "mesh_convergence_modes_1_6": 1.0e-2,
    "three_d_vs_section_clamp_modes_1_3": 5.0e-2,
    "three_d_vs_section_clamp_modes_1_6": 1.0e-1,
    "requested_mode_counts": [24, 32, 40],
    "target_mode_count": 7,
}

CAPABILITY_COLUMNS = (
    "headless_scripting",
    "3D_solid_eigenfrequency",
    "full_orthotropic_elasticity",
    "local_material_orientation",
    "quadratic_hex",
    "quadratic_tet",
    "reference_point_with_6_dof",
    "kinematic_face_coupling",
    "MPC_or_equation_constraints",
    "mass_matrix",
    "mode_shape_export",
    "node_element_export",
    "result_parser",
)

EXECUTABLE_CANDIDATES = (
    ("Abaqus", "abaqus"),
    ("Abaqus", "abq*"),
    ("ANSYS", "ansys"),
    ("ANSYS MAPDL", "mapdl"),
    ("CalculiX", "ccx"),
    ("Gmsh", "gmsh"),
    ("Elmer", "ElmerSolver"),
    ("SALOME", "salome"),
    ("FreeCAD", "FreeCADCmd"),
    ("COMSOL", "comsol"),
)
PACKAGE_CANDIDATES = (
    "gmsh",
    "meshio",
    "pyvista",
    "sfepy",
    "dolfinx",
    "fenics",
    "petsc4py",
    "slepc4py",
    "ansys.mapdl.core",
)
PACKAGE_DISTRIBUTIONS = {"ansys.mapdl.core": "ansys-mapdl-core"}


def _tracked_material_inventory() -> dict[str, float | None]:
    """Return only constants actually represented by the tracked material."""

    material = hms_dx_209_material()
    return {
        "E1": float(material.E1_real),
        "E2": float(material.E2_real),
        "E3": None,
        "G12": float(material.G12_real),
        "G13": float(material.G13_real),
        "G23": float(material.G23_real),
        "nu12": float(material.nu12),
        "nu13": None,
        "nu23": None,
        "rho": float(material.rho),
    }


def _missing_engineering_constants(values: dict[str, Any]) -> list[str]:
    return [name for name in REQUIRED_ENGINEERING_CONSTANTS if values.get(name) is None]


def _poisson_convention_is_unambiguous(value: Any) -> bool:
    if not isinstance(value, str) or not value.strip():
        return False
    compact = "".join(value.lower().split()).replace("ν", "nu")
    return (
        ("s_ij=-nu_ij/e_i" in compact or "sij=-nuij/ei" in compact)
        or ("s12=-nu12/e1" in compact)
    )


def _orthotropic_tensor_audit(values: dict[str, Any]) -> dict[str, Any]:
    """Validate a complete orthotropic tensor in engineering notation."""

    missing = _missing_engineering_constants(values)
    errors: list[str] = []
    if missing:
        errors.append("missing constants: " + ", ".join(missing))
        return {
            "complete": False,
            "ready": False,
            "missing_constants": missing,
            "errors": errors,
        }

    numeric: dict[str, float] = {}
    for name in REQUIRED_ENGINEERING_CONSTANTS:
        try:
            numeric[name] = float(values[name])
        except (TypeError, ValueError):
            errors.append(f"{name} is not numeric")
    if errors:
        return {
            "complete": True,
            "ready": False,
            "missing_constants": [],
            "errors": errors,
        }

    for name, value in numeric.items():
        if not math.isfinite(value):
            errors.append(f"{name} is not finite")
    for name in (*ELASTIC_MODULI, "rho"):
        if numeric[name] <= 0.0:
            errors.append(f"{name} must be positive")
    if errors:
        return {
            "complete": True,
            "ready": False,
            "missing_constants": [],
            "errors": errors,
        }

    E1, E2, E3 = numeric["E1"], numeric["E2"], numeric["E3"]
    G12, G13, G23 = numeric["G12"], numeric["G13"], numeric["G23"]
    nu12, nu13, nu23 = numeric["nu12"], numeric["nu13"], numeric["nu23"]
    compliance = np.zeros((6, 6), dtype=float)
    compliance[0, 0] = 1.0 / E1
    compliance[1, 1] = 1.0 / E2
    compliance[2, 2] = 1.0 / E3
    compliance[3, 3] = 1.0 / G23
    compliance[4, 4] = 1.0 / G13
    compliance[5, 5] = 1.0 / G12
    compliance[0, 1] = compliance[1, 0] = -nu12 / E1
    compliance[0, 2] = compliance[2, 0] = -nu13 / E1
    compliance[1, 2] = compliance[2, 1] = -nu23 / E2

    symmetry_residual = float(np.max(np.abs(compliance - compliance.T)))
    compliance_eigenvalues = np.linalg.eigvalsh(compliance)
    compliance_spd = bool(np.all(compliance_eigenvalues > 0.0))
    stiffness: np.ndarray | None = None
    stiffness_eigenvalues = np.full(6, np.nan)
    stiffness_spd = False
    stiffness_condition = float("inf")
    if compliance_spd:
        stiffness = np.linalg.inv(compliance)
        stiffness_eigenvalues = np.linalg.eigvalsh(stiffness)
        stiffness_spd = bool(np.all(stiffness_eigenvalues > 0.0))
        stiffness_condition = float(np.linalg.cond(stiffness))

    nu21 = nu12 * E2 / E1
    nu31 = nu13 * E3 / E1
    nu32 = nu23 * E3 / E2
    reciprocal_pairs = (
        (nu12 / E1, nu21 / E2),
        (nu13 / E1, nu31 / E3),
        (nu23 / E2, nu32 / E3),
    )
    reciprocity_residuals = [
        abs(left - right) / max(abs(left), abs(right), np.finfo(float).tiny)
        for left, right in reciprocal_pairs
    ]
    theta0_recovery_residuals = {
        "E1": abs(1.0 / compliance[0, 0] - E1) / E1,
        "G12": abs(1.0 / compliance[5, 5] - G12) / G12,
        "G13": abs(1.0 / compliance[4, 4] - G13) / G13,
    }
    maximum_reciprocity_residual = float(max(reciprocity_residuals))
    maximum_theta0_recovery_residual = float(
        max(theta0_recovery_residuals.values())
    )
    ready = (
        compliance_spd
        and stiffness_spd
        and symmetry_residual <= TENSOR_NUMERICAL_TOLERANCE
        and maximum_reciprocity_residual <= TENSOR_NUMERICAL_TOLERANCE
        and maximum_theta0_recovery_residual <= TENSOR_NUMERICAL_TOLERANCE
    )
    if not compliance_spd:
        errors.append("compliance matrix is not positive definite")
    if compliance_spd and not stiffness_spd:
        errors.append("stiffness matrix is not positive definite")
    if maximum_reciprocity_residual > TENSOR_NUMERICAL_TOLERANCE:
        errors.append("Poisson reciprocity residual exceeds the numerical tolerance")
    if maximum_theta0_recovery_residual > TENSOR_NUMERICAL_TOLERANCE:
        errors.append("theta=0 reduced-constant recovery exceeds the numerical tolerance")

    return {
        "complete": True,
        "ready": bool(ready),
        "missing_constants": [],
        "errors": errors,
        "numeric_constants": numeric,
        "compliance_matrix": compliance,
        "stiffness_matrix": stiffness,
        "symmetry_residual": symmetry_residual,
        "compliance_positive_definite": compliance_spd,
        "stiffness_positive_definite": stiffness_spd,
        "compliance_minimum_eigenvalue": float(np.min(compliance_eigenvalues)),
        "stiffness_minimum_eigenvalue": float(np.min(stiffness_eigenvalues)),
        "compliance_condition_number": float(np.linalg.cond(compliance)),
        "stiffness_condition_number": stiffness_condition,
        "derived_poisson_ratios": {"nu21": nu21, "nu31": nu31, "nu32": nu32},
        "reciprocity_residuals": reciprocity_residuals,
        "maximum_reciprocity_residual": maximum_reciprocity_residual,
        "theta0_recovery_residuals": theta0_recovery_residuals,
        "maximum_theta0_recovery_residual": maximum_theta0_recovery_residual,
        "tensor_numerical_tolerance": TENSOR_NUMERICAL_TOLERANCE,
    }


def _audit_material_payload(payload: dict[str, Any]) -> dict[str, Any]:
    missing_provenance = [
        name for name in PROVENANCE_FIELDS if not payload.get(name)
    ]
    tensor = _orthotropic_tensor_audit(payload)
    source_confirmed = payload.get("source_audit_status") == "SOURCE_CONFIRMED"
    convention_ok = _poisson_convention_is_unambiguous(
        payload.get("poisson_convention")
    )
    ready = bool(
        tensor.get("ready")
        and source_confirmed
        and convention_ok
        and not missing_provenance
    )
    errors = list(tensor.get("errors", []))
    if missing_provenance:
        errors.append("missing provenance: " + ", ".join(missing_provenance))
    if not source_confirmed:
        errors.append("source_audit_status is not SOURCE_CONFIRMED")
    if not convention_ok:
        errors.append("Poisson convention is absent or ambiguous")
    return {
        **tensor,
        "ready": ready,
        "errors": errors,
        "source_confirmed": source_confirmed,
        "poisson_convention_unambiguous": convention_ok,
        "missing_provenance": missing_provenance,
        "provenance_complete": not missing_provenance,
        "material_name": payload.get("material_name"),
        "provenance": {name: payload.get(name) for name in PROVENANCE_FIELDS},
    }


def _load_material_audit(material_json: Path | None) -> dict[str, Any]:
    if material_json is None:
        values = _tracked_material_inventory()
        return {
            "input_kind": "tracked_reduced_HMS_DX_209",
            "ready": False,
            "status": STATUS_MATERIAL,
            "values": values,
            "missing_constants": _missing_engineering_constants(values),
            "source_confirmed": True,
            "source_evidence": {
                "source_file": SOURCE_FILE,
                "engineering_constants_page": SOURCE_TABLE_PAGE,
                "engineering_constants_table": SOURCE_TABLE,
                "poisson_convention_page": SOURCE_CONVENTION_PAGE,
                "poisson_convention_equation": SOURCE_CONVENTION_EQUATION,
            },
            "errors": [
                "tracked and source-audited data do not provide E3, nu13, or nu23"
            ],
        }

    path = Path(material_json)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return {
            "input_kind": "material_json",
            "ready": False,
            "status": STATUS_INCONCLUSIVE,
            "values": {},
            "missing_constants": list(REQUIRED_ENGINEERING_CONSTANTS),
            "source_confirmed": False,
            "errors": [f"cannot read material JSON: {exc}"],
        }
    if not isinstance(payload, dict):
        return {
            "input_kind": "material_json",
            "ready": False,
            "status": STATUS_MATERIAL,
            "values": {},
            "missing_constants": list(REQUIRED_ENGINEERING_CONSTANTS),
            "source_confirmed": False,
            "errors": ["material JSON root must be an object"],
        }
    audit = _audit_material_payload(payload)
    return {
        "input_kind": "material_json",
        "status": STATUS_READY if audit["ready"] else STATUS_MATERIAL,
        "values": {name: payload.get(name) for name in REQUIRED_ENGINEERING_CONSTANTS},
        **audit,
    }


def _find_abq_executables() -> list[str]:
    matches: set[str] = set()
    for entry in os.environ.get("PATH", "").split(os.pathsep):
        if not entry:
            continue
        directory = Path(entry)
        if not directory.is_dir():
            continue
        for pattern in ("abq*.exe", "abq*.bat", "abq*.cmd"):
            for path in directory.glob(pattern):
                if path.is_file():
                    matches.add(str(path.resolve()))
    return sorted(matches)


def _safe_package_version(package: str) -> str:
    distribution = PACKAGE_DISTRIBUTIONS.get(package, package)
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return "UNKNOWN_METADATA"


def _capability_values(solver: str, detected: bool) -> dict[str, str]:
    if not detected:
        value = "NOT_AUDITED_NOT_DETECTED"
    elif solver in {"Gmsh", "FreeCAD", "meshio", "pyvista"}:
        value = "NO_MESH_OR_POSTPROCESS_ONLY"
    else:
        value = "UNVERIFIED_LOCAL_RUNTIME"
    return {name: value for name in CAPABILITY_COLUMNS}


def _solver_inventory(
    *,
    which: Callable[[str], str | None] = shutil.which,
    find_spec: Callable[[str], Any] = importlib.util.find_spec,
    abq_finder: Callable[[], list[str]] = _find_abq_executables,
) -> list[dict[str, Any]]:
    """Passively inventory candidates; never launch or import a solver."""

    rows: list[dict[str, Any]] = []
    for solver, command in EXECUTABLE_CANDIDATES:
        detected_path: str | None
        if command == "abq*":
            matches = abq_finder()
            detected_path = matches[0] if matches else None
            extra = f"matches={matches}" if matches else "no PATH match"
        else:
            detected_path = which(command)
            extra = "passive PATH lookup only"
        detected = detected_path is not None
        row: dict[str, Any] = {
            "solver": solver,
            "executable_or_package": command,
            "detected": detected,
            "detected_path_or_origin": detected_path or "",
            "version": "NOT_QUERIED_NO_PROCESS_LAUNCH" if detected else "NOT_DETECTED",
            "license_state": "AVAILABLE_LICENSE_UNVERIFIED" if detected else "NOT_DETECTED",
            "mandatory_capabilities_pass": False,
            "notes": (
                f"{extra}; presence alone does not confirm the mandatory capability set"
            ),
        }
        row.update(_capability_values(solver, detected))
        rows.append(row)

    for package in PACKAGE_CANDIDATES:
        try:
            spec = find_spec(package)
        except (ImportError, ModuleNotFoundError, AttributeError, ValueError):
            spec = None
        detected = spec is not None
        row = {
            "solver": package,
            "executable_or_package": f"python:{package}",
            "detected": detected,
            "detected_path_or_origin": getattr(spec, "origin", "") if detected else "",
            "version": _safe_package_version(package) if detected else "NOT_DETECTED",
            "license_state": "PACKAGE_METADATA_ONLY" if detected else "NOT_DETECTED",
            "mandatory_capabilities_pass": False,
            "notes": (
                "find_spec/importlib.metadata only; package was not imported; "
                "a Python package alone is not classified as a ready 3D solver"
            ),
        }
        row.update(_capability_values(package, detected))
        rows.append(row)

    rows.append(
        {
            "solver": "CoupledBeams python_fem.py",
            "executable_or_package": "src/my_project/fem/python_fem.py",
            "detected": (ROOT / "src/my_project/fem/python_fem.py").is_file(),
            "detected_path_or_origin": "src/my_project/fem/python_fem.py",
            "version": "tracked 1D baseline",
            "license_state": "NOT_APPLICABLE",
            **{name: "NO_1D_FRAME_BEAM_BASELINE" for name in CAPABILITY_COLUMNS},
            "mandatory_capabilities_pass": False,
            "notes": "1D planar frame/beam FEM; explicitly excluded as a 3D solid solver",
        }
    )
    return rows


def _anchor_case_manifest() -> list[dict[str, Any]]:
    common = {
        "material": "HMS/DX-209 full 3D tensor required",
        "theta_deg": 0.0,
        "local_e1": "t_i (outer clamp to joint)",
        "local_e2": "n_i=e_z cross t_i",
        "local_e3": "global e_z",
        "outer_boundary": "all three translations fixed on each outer face",
        "primary_1D_comparator": "state_corrected + timoshenko_section_clamp",
        "source_diagnostic": "state_corrected + book_slope_clamp",
        "EB_diagnostic": "rectangular EB",
    }
    return [
        {
            "case_id": "S0_uniform_5mm",
            "future_stage": "3D-A1",
            "geometry": "one straight prism",
            "arm_count": 1,
            "a1_m": 0.005,
            "a2_m": "",
            "b1_m": 0.020,
            "b2_m": "",
            "L1_m": 0.400,
            "L2_m": "",
            "beta_deg": 0.0,
            "joint_contract": "not applicable; both end faces fully clamped",
            "role": "material/orientation/element/parity/mesh-convergence anchor",
            **common,
        },
        {
            "case_id": "J30_equal_5_5",
            "future_stage": "3D-A2",
            "geometry": "two independent overlapping full prisms",
            "arm_count": 2,
            "a1_m": 0.005,
            "a2_m": 0.005,
            "b1_m": 0.020,
            "b2_m": 0.020,
            "L1_m": 0.200,
            "L2_m": 0.200,
            "beta_deg": 30.0,
            "joint_contract": "one massless 6-DOF RP; exact rigid-face kinematics",
            "role": "generic equal-thickness joint anchor",
            **common,
        },
        {
            "case_id": "J30_unequal_4_6",
            "future_stage": "3D-A3",
            "geometry": "two independent overlapping full prisms",
            "arm_count": 2,
            "a1_m": 0.004,
            "a2_m": 0.006,
            "b1_m": 0.020,
            "b2_m": 0.020,
            "L1_m": 0.200,
            "L2_m": 0.200,
            "beta_deg": 30.0,
            "joint_contract": "one massless 6-DOF RP; exact rigid-face kinematics",
            "role": "primary unequal-thickness anchor",
            **common,
        },
        {
            "case_id": "J90_unequal_4_6",
            "future_stage": "3D-A3",
            "geometry": "two independent overlapping full prisms",
            "arm_count": 2,
            "a1_m": 0.004,
            "a2_m": 0.006,
            "b1_m": 0.020,
            "b2_m": 0.020,
            "L1_m": 0.200,
            "L2_m": 0.200,
            "beta_deg": 90.0,
            "joint_contract": "one massless 6-DOF RP; exact rigid-face kinematics",
            "role": "quarter-turn demanding anchor",
            **common,
        },
    ]


def _mesh_policy() -> list[dict[str, Any]]:
    return [
        {
            "level": "M0",
            "element_policy": "quadratic structured hexahedra preferred",
            "N_L_200mm": 50,
            "N_L_S0_400mm": 100,
            "N_b_20mm": 10,
            "N_a_4mm": 2,
            "N_a_5mm": 3,
            "N_a_6mm": 3,
            "tet_fallback": "quadratic tetrahedra with separate convergence only",
            "linear_tet_primary": False,
        },
        {
            "level": "M1",
            "element_policy": "quadratic structured hexahedra preferred",
            "N_L_200mm": 75,
            "N_L_S0_400mm": 150,
            "N_b_20mm": 15,
            "N_a_4mm": 3,
            "N_a_5mm": 4,
            "N_a_6mm": 5,
            "tet_fallback": "quadratic tetrahedra with separate convergence only",
            "linear_tet_primary": False,
        },
        {
            "level": "M2",
            "element_policy": "quadratic structured hexahedra preferred",
            "N_L_200mm": 100,
            "N_L_S0_400mm": 200,
            "N_b_20mm": 20,
            "N_a_4mm": 4,
            "N_a_5mm": 5,
            "N_a_6mm": 6,
            "tet_fallback": "quadratic tetrahedra with separate convergence only",
            "linear_tet_primary": False,
        },
    ]


def _orientation_matrix(beta_deg: float, arm: int) -> np.ndarray:
    if arm not in (1, 2):
        raise ValueError("arm must be 1 or 2")
    e_z = np.array([0.0, 0.0, 1.0])
    t_1 = np.array([1.0, 0.0, 0.0])
    n_1 = np.cross(e_z, t_1)
    if arm == 1:
        tangent, normal = t_1, n_1
    else:
        beta = math.radians(float(beta_deg))
        tangent = -math.cos(beta) * t_1 + math.sin(beta) * n_1
        normal = np.cross(e_z, tangent)
    return np.column_stack((tangent, normal, e_z))


def _orientation_design_audit() -> dict[str, Any]:
    matrices = []
    for beta in (0.0, 30.0, 90.0):
        for arm in (1, 2):
            matrix = _orientation_matrix(beta, arm)
            matrices.append(
                {
                    "beta_deg": beta,
                    "arm": arm,
                    "matrix": matrix.tolist(),
                    "determinant": float(np.linalg.det(matrix)),
                    "orthogonality_residual": float(
                        np.max(np.abs(matrix.T @ matrix - np.eye(3)))
                    ),
                }
            )
    return {
        "matrices": matrices,
        "maximum_determinant_residual": max(
            abs(row["determinant"] - 1.0) for row in matrices
        ),
        "maximum_orthogonality_residual": max(
            row["orthogonality_residual"] for row in matrices
        ),
        "threshold": FUTURE_THRESHOLDS["orientation_orthogonality_residual"],
    }


def _classify_readiness(
    *, material_ready: bool, solver_ready: bool, inconclusive: bool = False
) -> str:
    if inconclusive:
        return STATUS_INCONCLUSIVE
    if material_ready and solver_ready:
        return STATUS_READY
    if not material_ready and not solver_ready:
        return STATUS_BOTH
    if not material_ready:
        return STATUS_MATERIAL
    return STATUS_SOLVER


def _git_context() -> dict[str, str]:
    def git(*args: str) -> str:
        try:
            return subprocess.check_output(
                ["git", *args], cwd=ROOT, text=True, stderr=subprocess.STDOUT
            ).strip()
        except (OSError, subprocess.CalledProcessError) as exc:
            return f"UNAVAILABLE: {exc}"

    return {
        "cwd": str(Path.cwd()),
        "git_root": git("rev-parse", "--show-toplevel"),
        "branch": git("branch", "--show-current"),
        "HEAD": git("rev-parse", "HEAD"),
        "origin_main": git("rev-parse", "origin/main"),
        "status_short": git("status", "--short"),
        "diff_stat": git("diff", "--stat"),
        "untracked": git("ls-files", "--others", "--exclude-standard"),
    }


def _csv_safe(value: Any) -> Any:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (list, dict, tuple)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return value


def _write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: Sequence[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: _csv_safe(row.get(name, "")) for name in fieldnames})


def _constitutive_rows(material_audit: dict[str, Any]) -> list[dict[str, Any]]:
    values = material_audit.get("values", {})
    is_tracked = material_audit.get("input_kind") == "tracked_reduced_HMS_DX_209"
    provenance = material_audit.get("provenance", {})
    rows: list[dict[str, Any]] = []
    for name in REQUIRED_ENGINEERING_CONSTANTS:
        value = values.get(name)
        if name in ELASTIC_MODULI:
            units = "Pa"
        elif name == "rho":
            units = "kg/m^3"
        else:
            units = "dimensionless"
        available = value is not None
        source_supported = (
            is_tracked and name not in {"E3", "nu13", "nu23"}
        ) or bool(
            material_audit.get("source_confirmed")
            and material_audit.get("provenance_complete")
            and available
        )
        rows.append(
            {
                "constant_or_check": name,
                "value": value if available else "MISSING",
                "units": units,
                "availability": "AVAILABLE" if available else "MISSING",
                "source_supported": source_supported,
                "source_file": (
                    SOURCE_FILE if is_tracked else provenance.get("source_file", "")
                )
                if source_supported
                else "",
                "source_printed_page": (
                    SOURCE_TABLE_PAGE
                    if is_tracked
                    else provenance.get("source_printed_page", "")
                )
                if source_supported
                else "",
                "source_table_or_equation": (
                    SOURCE_TABLE
                    if is_tracked
                    else provenance.get("source_table_or_equation", "")
                )
                if source_supported
                else "",
                "poisson_convention": POISSON_CONVENTION if name.startswith("nu") else "",
                "status": "PASS" if source_supported else "BLOCKED",
                "notes": (
                    "No transverse-isotropy completion was assumed"
                    if not available
                    else "tracked/source value"
                ),
            }
        )
    if material_audit.get("complete"):
        for name in (
            "symmetry_residual",
            "compliance_minimum_eigenvalue",
            "stiffness_minimum_eigenvalue",
            "compliance_condition_number",
            "stiffness_condition_number",
            "maximum_reciprocity_residual",
            "maximum_theta0_recovery_residual",
        ):
            rows.append(
                {
                    "constant_or_check": name,
                    "value": material_audit.get(name),
                    "units": "diagnostic",
                    "availability": "COMPUTED",
                    "source_supported": material_audit.get("source_confirmed", False),
                    "source_file": material_audit.get("provenance", {}).get(
                        "source_file", ""
                    ),
                    "source_printed_page": material_audit.get("provenance", {}).get(
                        "source_printed_page", ""
                    ),
                    "source_table_or_equation": material_audit.get("provenance", {}).get(
                        "source_table_or_equation", ""
                    ),
                    "poisson_convention": material_audit.get("provenance", {}).get(
                        "poisson_convention", ""
                    ),
                    "status": "PASS" if material_audit.get("ready") else "FAIL",
                    "notes": "engineering-notation tensor audit",
                }
            )
    return rows


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def _readiness_report(summary: dict[str, Any]) -> str:
    missing = ", ".join(summary["material_readiness"]["missing_constants"]) or "none"
    solver_detected = [
        row["executable_or_package"]
        for row in summary["solver_inventory"]
        if row["detected"] and row["solver"] != "CoupledBeams python_fem.py"
    ]
    solver_text = ", ".join(solver_detected) if solver_detected else "none"
    material_dimension = summary["readiness_dimensions"]["material"]
    solver_dimension = summary["readiness_dimensions"]["solver"]
    material_evidence = (
        "Source-confirmed full tensor passed symmetry, reciprocity, and SPD checks."
        if material_dimension == "READY"
        else f"Required full tensor is incomplete or unconfirmed; missing: {missing}."
    )
    solver_evidence = (
        "A local candidate passes every mandatory capability and runtime gate."
        if solver_dimension == "READY"
        else f"Detected external solver/package candidates: {solver_text}; no candidate passes all mandatory capabilities."
    )
    return f"""# 3D-A0 limited 3D FEM anchor readiness report

## Status

```text
Unequal-thickness 1D validation phase:
COMPLETE

Overall unequal-thickness 1D validation:
PARTIAL_PASS

Limited 3D FEM anchor phase:
DESIGN_ONLY

3D-A0 design and readiness audit:
{summary['readiness_status']}

3D mesh generation:
NOT_STARTED

3D eigenfrequency calculations:
NOT_STARTED
```

## Readiness dimensions

| Dimension | Status | Evidence |
| --- | --- | --- |
| material | {material_dimension} | {material_evidence} |
| solver | {solver_dimension} | {solver_evidence} |
| geometry | DESIGN_READY | Two-full-prism overlapping-part abstraction and axes are frozen. |
| boundary condition | DESIGN_READY | Fully fixed outer faces and one zero-mass exact-kinematic 6-DOF joint RP are frozen. |
| mode extraction | DESIGN_READY | Reflection parity, clustered-space purification, and 24→32→40 policy are frozen. |
| implementation | NOT_STARTED | No mesh, deck, solver run, result parser, or 3D frequency output exists. |

## Material source audit

The tracked HMS/DX-209 constants and Chapter-1 Table 1.1 (printed p.45) provide
`E1=191 GPa`, `E2=5 GPa`, `G12=3 GPa`, `G13=3 GPa`, `G23=2.5 GPa`,
`nu12=0.279`, and `rho=1580 kg/m^3`. Equation (1.50), printed p.40, gives
`S12=-nu12/E1=-nu21/E2`. The audited pages do not provide `E3`, `nu13`, or
`nu23`, and no transverse-isotropy relations were assumed.

## Frozen design

Coupled cases use two independent full prisms, allowed to overlap geometrically
but never united, merged, contacted, or mass-subtracted. Both joint-face centers
are at the origin. Local axes are `e1=t_i`, `e2=n_i=e_z×t_i`, `e3=e_z`, with
`t1=(1,0,0)` and `t2=-cos(beta)t1+sin(beta)n1`. One massless six-DOF reference
point imposes exact rigid-face kinematics on both inner faces; each outer face
has all three translations fixed.

The primary future 1D comparator is `state_corrected` Timoshenko with
`timoshenko_section_clamp`. The source-facing `book_slope_clamp` and rectangular
EB spectra are diagnostics only; existing UT-0--UT-3 outputs are unchanged.

The target family obeys `phi≈-R_z phi`, where `Q_z=diag(1,1,-1)`: `u_x,u_y`
are odd and `u_z` is even in `z`. Future classification uses `p<=-0.99` or
`r_-<=1e-3`; clusters with relative gap `<1e-3` require parity-matrix
diagonalization before classification.

Exactly four cases are frozen: `S0_uniform_5mm`, `J30_equal_5_5`,
`J30_unequal_4_6`, and `J90_unequal_4_6`. Mesh levels are exactly M0/M1/M2,
with quadratic structured hexahedra preferred and linear tetrahedra excluded
as the primary anchor.

## Explicit exclusions

No 3D mesh, geometry file, solver deck, solver or license-consuming job,
eigenfrequency, mode shape, package installation, material completion,
production API, 3D-A1/A2/A3 calculation, sweep, plot, or PDF was created.
"""


def _run_readiness_audit(
    output_dir: Path, material_json: Path | None = None
) -> dict[str, Any]:
    material_audit = _load_material_audit(material_json)
    solver_rows = _solver_inventory()
    solver_ready_rows = [row for row in solver_rows if row["mandatory_capabilities_pass"]]
    solver_ready = bool(solver_ready_rows)
    material_inconclusive = material_audit.get("status") == STATUS_INCONCLUSIVE
    readiness_status = _classify_readiness(
        material_ready=bool(material_audit.get("ready")),
        solver_ready=solver_ready,
        inconclusive=material_inconclusive,
    )
    cases = _anchor_case_manifest()
    mesh_rows = _mesh_policy()
    orientation_audit = _orientation_design_audit()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary: dict[str, Any] = {
        "audit": "3D-A0 limited 3D FEM anchor design and readiness audit",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "git_context": _git_context(),
        "one_dimensional_regression_status": ONE_D_STATUSES,
        "limited_3d_phase": "DESIGN_ONLY",
        "readiness_status": readiness_status,
        "readiness_dimensions": {
            "material": "READY" if material_audit.get("ready") else "BLOCKED",
            "solver": "READY" if solver_ready else "BLOCKED",
            "geometry": "DESIGN_READY",
            "boundary_condition": "DESIGN_READY",
            "mode_extraction": "DESIGN_READY",
            "implementation": "NOT_STARTED",
        },
        "material_readiness": material_audit,
        "required_engineering_constants": list(REQUIRED_ENGINEERING_CONSTANTS),
        "solver_inventory": solver_rows,
        "primary_solver_candidate": solver_ready_rows[0]["solver"] if solver_ready_rows else None,
        "fallback_solver_candidate": solver_ready_rows[1]["solver"] if len(solver_ready_rows) > 1 else None,
        "mandatory_solver_capabilities": list(CAPABILITY_COLUMNS),
        "geometry_contract": {
            "coupled_parts": "two independent full rectangular prisms",
            "overlap": "allowed geometrically; no union, merge, contact, shared elements, or mass subtraction",
            "joint_center_m": [0.0, 0.0, 0.0],
            "outer_center": "-L_i t_i",
            "axes": "e1=t_i, e2=n_i=e_z cross t_i, e3=e_z",
            "t1": [1.0, 0.0, 0.0],
            "t2": "-cos(beta)t1+sin(beta)n1",
            "orientation_design_audit": orientation_audit,
            "future_mass_gate": "rho*(A1*L1+A2*L2), relative residual <=1e-10; S0=rho*A*L",
        },
        "joint_contract": {
            "reference_point": "one massless 6-DOF RP at the origin",
            "mass": 0.0,
            "rotary_inertia": 0.0,
            "coupling": "exact kinematic rigid-face coupling for both inner faces",
            "future_rigid_motion_smoke": "u(r)=U_J+Omega_J cross (r-r_J) for six unit RP motions",
            "forbidden": ["penalty", "contact", "face merge", "fixed RP", "joint mass"],
        },
        "boundary_contract": {
            "solid_outer_faces": "u_x=u_y=u_z=0 at every outer-face node",
            "primary_1D_comparator": "Timoshenko state_corrected + timoshenko_section_clamp",
            "source_diagnostic": "Timoshenko state_corrected + book_slope_clamp",
            "EB_diagnostic": "rectangular EB (w'=psi, clamp variants coincide)",
        },
        "mode_extraction_contract": {
            "reflection": "R_z u(x,y,z)=diag(1,1,-1)u(x,y,-z)",
            "target": "phi approximately -R_z phi; u_x,u_y odd and u_z even",
            "parity_criterion": "p<=-0.99 or r_-<=1e-3",
            "cluster_rule": "gap<1e-3: diagonalize V_c^T M R_z V_c before classification",
            "mode_count_sequence": [24, 32, 40],
            "failure_after_40": "TARGET_FAMILY_NOT_FOUND_WITHIN_LIMIT",
        },
        "anchor_cases": cases,
        "mesh_policy": mesh_rows,
        "future_stage_plan": {
            "3D-A1": "S0_uniform_5mm at M0/M1/M2",
            "3D-A2": "J30_equal_5_5 at M2 only",
            "3D-A3": "J30_unequal_4_6 at M2; J90_unequal_4_6 at M1/M2",
        },
        "frozen_future_thresholds": FUTURE_THRESHOLDS,
        "mesh_generation": "NOT_STARTED",
        "eigenfrequency_calculations": "NOT_STARTED",
        "explicit_exclusions": [
            "3D mesh generation",
            "3D eigenfrequency calculation",
            "solver or GUI launch",
            "license-consuming job",
            "package installation or download",
            "assumed constitutive completion",
            "3D-A1, 3D-A2, or 3D-A3 execution",
            "swapped, negative-beta, intermediate-beta, or parameter sweep cases",
            "plots or PDF",
            "production API changes",
        ],
    }

    solver_fields = [
        "solver",
        "executable_or_package",
        "detected",
        "detected_path_or_origin",
        "version",
        "license_state",
        *CAPABILITY_COLUMNS,
        "mandatory_capabilities_pass",
        "notes",
    ]
    constitutive_fields = [
        "constant_or_check",
        "value",
        "units",
        "availability",
        "source_supported",
        "source_file",
        "source_printed_page",
        "source_table_or_equation",
        "poisson_convention",
        "status",
        "notes",
    ]
    case_fields = list(cases[0])
    mesh_fields = list(mesh_rows[0])
    _write_csv(output_dir / "solver_inventory.csv", solver_rows, solver_fields)
    _write_csv(
        output_dir / "constitutive_audit.csv",
        _constitutive_rows(material_audit),
        constitutive_fields,
    )
    _write_csv(output_dir / "anchor_case_manifest.csv", cases, case_fields)
    _write_csv(output_dir / "mesh_policy.csv", mesh_rows, mesh_fields)
    (output_dir / "readiness_summary.json").write_text(
        json.dumps(_jsonable(summary), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    (output_dir / "readiness_report.md").write_text(
        _readiness_report(summary), encoding="utf-8"
    )
    return summary


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="directory for the six text-only readiness artifacts",
    )
    parser.add_argument(
        "--material-json",
        type=Path,
        help="optional provenance-bearing full engineering-constant JSON",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    summary = _run_readiness_audit(args.output_dir, args.material_json)
    print(f"3D-A0 design and readiness audit: {summary['readiness_status']}")
    print("3D mesh generation: NOT_STARTED")
    print("3D eigenfrequency calculations: NOT_STARTED")
    print(f"Artifacts: {Path(args.output_dir).resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
