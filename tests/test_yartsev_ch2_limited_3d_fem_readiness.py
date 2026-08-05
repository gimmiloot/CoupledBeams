import ast
import hashlib
import importlib
import json
import subprocess
import sys
from pathlib import Path
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis.anisotropic_rods import (
    audit_yartsev_ch2_limited_3d_fem_readiness as readiness,
)


SCRIPT = (
    ROOT
    / "scripts"
    / "analysis"
    / "anisotropic_rods"
    / "audit_yartsev_ch2_limited_3d_fem_readiness.py"
)
UNEQUAL_SCRIPT = (
    ROOT
    / "scripts"
    / "analysis"
    / "anisotropic_rods"
    / "validate_yartsev_ch2_unequal_thickness.py"
)
UNEQUAL_SCRIPT_PRE_3D_A0_SHA256 = (
    "e947aaf4f9cb6845ff030c11400231ed809f5cfe1a1af44f852b64796e9c7e75"
)


def _synthetic_constants() -> dict[str, float]:
    return {
        "E1": 140.0e9,
        "E2": 10.0e9,
        "E3": 8.0e9,
        "G12": 5.0e9,
        "G13": 4.0e9,
        "G23": 3.0e9,
        "nu12": 0.25,
        "nu13": 0.20,
        "nu23": 0.30,
        "rho": 1600.0,
    }


def _provenance_payload() -> dict[str, object]:
    return {
        **_synthetic_constants(),
        "material_name": "synthetic orthotropic test material",
        "units": "SI: Pa, kg/m^3, dimensionless",
        "source_file": "synthetic_test_source.pdf",
        "source_printed_page": "1",
        "source_table_or_equation": "test table",
        "poisson_convention": "S12=-nu12/E1; S_ij=-nu_ij/E_i",
        "source_audit_status": "SOURCE_CONFIRMED",
    }


def test_readiness_script_import_does_not_run_audit() -> None:
    with mock.patch.object(
        Path, "mkdir", side_effect=AssertionError("import attempted to write artifacts")
    ):
        module = importlib.reload(readiness)
    assert callable(module.main)


def test_tracked_reduced_material_inventory_is_exact() -> None:
    values = readiness._tracked_material_inventory()
    assert values == {
        "E1": 191.0e9,
        "E2": 5.0e9,
        "E3": None,
        "G12": 3.0e9,
        "G13": 3.0e9,
        "G23": 2.5e9,
        "nu12": 0.279,
        "nu13": None,
        "nu23": None,
        "rho": 1580.0,
    }


def test_missing_full_3d_constants_are_reported_not_filled() -> None:
    values = readiness._tracked_material_inventory()
    assert readiness._missing_engineering_constants(values) == ["E3", "nu13", "nu23"]
    assert all(values[name] is None for name in ("E3", "nu13", "nu23"))


def test_incomplete_json_is_blocked(tmp_path: Path) -> None:
    payload = _provenance_payload()
    payload.pop("E3")
    path = tmp_path / "material.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    audit = readiness._load_material_audit(path)
    assert audit["status"] == readiness.STATUS_MATERIAL
    assert audit["ready"] is False
    assert audit["missing_constants"] == ["E3"]


def test_json_without_provenance_is_not_ready() -> None:
    audit = readiness._audit_material_payload(_synthetic_constants())
    assert audit["ready"] is False
    assert set(audit["missing_provenance"]) == set(readiness.PROVENANCE_FIELDS)


def test_complete_synthetic_tensor_is_symmetric() -> None:
    audit = readiness._orthotropic_tensor_audit(_synthetic_constants())
    compliance = audit["compliance_matrix"]
    assert audit["symmetry_residual"] == 0.0
    np.testing.assert_array_equal(compliance, compliance.T)


def test_positive_definite_synthetic_tensor_passes() -> None:
    audit = readiness._orthotropic_tensor_audit(_synthetic_constants())
    assert audit["ready"] is True
    assert audit["compliance_positive_definite"] is True
    assert audit["stiffness_positive_definite"] is True
    assert audit["compliance_minimum_eigenvalue"] > 0.0
    assert audit["stiffness_minimum_eigenvalue"] > 0.0


def test_non_positive_definite_tensor_is_rejected() -> None:
    values = _synthetic_constants()
    values["nu12"] = 20.0
    audit = readiness._orthotropic_tensor_audit(values)
    assert audit["ready"] is False
    assert audit["compliance_positive_definite"] is False
    assert audit["compliance_minimum_eigenvalue"] < 0.0


def test_poisson_reciprocity_is_checked() -> None:
    audit = readiness._orthotropic_tensor_audit(_synthetic_constants())
    derived = audit["derived_poisson_ratios"]
    values = _synthetic_constants()
    assert np.isclose(
        derived["nu21"] / values["E2"], values["nu12"] / values["E1"]
    )
    assert np.isclose(
        derived["nu31"] / values["E3"], values["nu13"] / values["E1"]
    )
    assert np.isclose(
        derived["nu32"] / values["E3"], values["nu23"] / values["E2"]
    )
    assert audit["maximum_reciprocity_residual"] < 1.0e-15


def test_source_confirmed_complete_payload_can_be_material_ready() -> None:
    audit = readiness._audit_material_payload(_provenance_payload())
    assert audit["ready"] is True
    assert audit["source_confirmed"] is True
    assert audit["poisson_convention_unambiguous"] is True


def test_solver_detection_does_not_launch_gui_or_job() -> None:
    def fake_which(command: str) -> str | None:
        return r"C:\fake\abaqus.exe" if command == "abaqus" else None

    with mock.patch.object(subprocess, "Popen", side_effect=AssertionError("process launched")):
        rows = readiness._solver_inventory(
            which=fake_which,
            find_spec=lambda _name: None,
            abq_finder=lambda: [],
        )
    abaqus = next(row for row in rows if row["executable_or_package"] == "abaqus")
    assert abaqus["detected"] is True
    assert abaqus["version"] == "NOT_QUERIED_NO_PROCESS_LAUNCH"
    assert abaqus["license_state"] == "AVAILABLE_LICENSE_UNVERIFIED"
    assert abaqus["mandatory_capabilities_pass"] is False


def test_existing_python_fem_is_not_classified_as_3d_solver() -> None:
    rows = readiness._solver_inventory(
        which=lambda _name: None,
        find_spec=lambda _name: None,
        abq_finder=lambda: [],
    )
    baseline = next(row for row in rows if row["solver"] == "CoupledBeams python_fem.py")
    assert baseline["detected"] is True
    assert baseline["mandatory_capabilities_pass"] is False
    assert baseline["3D_solid_eigenfrequency"] == "NO_1D_FRAME_BEAM_BASELINE"


def test_manifest_contains_exactly_the_four_frozen_cases() -> None:
    cases = readiness._anchor_case_manifest()
    assert [case["case_id"] for case in cases] == [
        "S0_uniform_5mm",
        "J30_equal_5_5",
        "J30_unequal_4_6",
        "J90_unequal_4_6",
    ]


def test_mesh_policy_contains_only_m0_m1_m2() -> None:
    policy = readiness._mesh_policy()
    assert [row["level"] for row in policy] == ["M0", "M1", "M2"]
    assert [row["N_L_S0_400mm"] for row in policy] == [100, 150, 200]
    assert all(row["linear_tet_primary"] is False for row in policy)


def test_no_3d_solve_function_is_defined() -> None:
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    names = {
        node.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    assert not any(name == "solve" or name.startswith("solve_") or name.endswith("_solve") for name in names)


def test_no_mesh_generation_function_is_defined() -> None:
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    names = {
        node.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    forbidden = {"generate_mesh", "build_mesh", "create_mesh", "write_mesh"}
    assert names.isdisjoint(forbidden)


def test_legacy_unequal_thickness_cli_is_byte_for_byte_unchanged() -> None:
    digest = hashlib.sha256(UNEQUAL_SCRIPT.read_bytes()).hexdigest()
    assert digest == UNEQUAL_SCRIPT_PRE_3D_A0_SHA256


def test_design_contract_orientation_is_right_handed() -> None:
    audit = readiness._orientation_design_audit()
    assert audit["maximum_determinant_residual"] <= 1.0e-15
    assert audit["maximum_orthogonality_residual"] <= 1.0e-12


def test_readiness_audit_writes_exactly_six_text_artifacts(tmp_path: Path) -> None:
    with mock.patch.object(
        readiness,
        "_solver_inventory",
        return_value=readiness._solver_inventory(
            which=lambda _name: None,
            find_spec=lambda _name: None,
            abq_finder=lambda: [],
        ),
    ):
        summary = readiness._run_readiness_audit(tmp_path)
    assert summary["readiness_status"] == readiness.STATUS_BOTH
    assert {path.name for path in tmp_path.iterdir()} == {
        "solver_inventory.csv",
        "constitutive_audit.csv",
        "anchor_case_manifest.csv",
        "mesh_policy.csv",
        "readiness_summary.json",
        "readiness_report.md",
    }
    assert summary["mesh_generation"] == "NOT_STARTED"
    assert summary["eigenfrequency_calculations"] == "NOT_STARTED"
