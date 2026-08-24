from __future__ import annotations

import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
import pytest


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis.laminated_beams import (  # noqa: E402
    audit_reddy_table_4_3_3_discrepancies as audit,
)


SOURCE_SHA256 = "3F479D04DBD601C4BBFEB9D32463956A5201D72ECC41A067B8234979D885C388"
LEGACY_CSV_SHA256 = "73D4F0D966A6F502F48D3900375B6FC1F0BACE67AC636FC1AA7D9893063B7D92"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def _relative_difference(value: float, reference: float) -> float:
    return abs(value - reference) / abs(reference)


@pytest.fixture(scope="module")
def full_audit() -> SimpleNamespace:
    contract, cases = audit.load_all_direct_cases()
    inventory, metadata = audit.load_legacy_transfer_inventory(contract)
    audits = audit.audit_direct_cases(cases, inventory)
    rotary_rows, source_violations = audit.build_rotary_inertia_rows(
        audits, inventory
    )
    return SimpleNamespace(
        contract=contract,
        cases=cases,
        inventory=inventory,
        metadata=metadata,
        audits=audits,
        rotary_rows=rotary_rows,
        source_violations=source_violations,
    )


def test_all_72_model_cases_and_54_published_cases_are_frozen(
    full_audit: SimpleNamespace,
) -> None:
    assert len(full_audit.cases) == 72
    assert len(full_audit.audits) == 72
    assert len({case.case_id for case in full_audit.cases}) == 72
    assert sum(case.source_value is not None for case in full_audit.cases) == 54
    assert sum(case.source_value is None for case in full_audit.cases) == 18
    assert {case.boundary_condition for case in full_audit.cases} == {
        "HH",
        "CC",
        "CF",
    }
    assert sum(case.rotary_inertia for case in full_audit.cases) == 36
    assert all(item.model_status == "PASS" for item in full_audit.audits)
    assert all(item.ritz_status == "PASS" for item in full_audit.audits)


def test_ritz_matrices_are_symmetric_and_energy_definite(
    full_audit: SimpleNamespace,
) -> None:
    for item in full_audit.audits:
        mode = item.converged_mode
        assert mode.stiffness_symmetry_residual <= 1.0e-14
        assert mode.mass_symmetry_residual <= 1.0e-14
        assert mode.stiffness_spd
        assert mode.stiffness_minimum_eigenvalue > 0.0
        assert mode.translational_mass_spd
        assert mode.translational_mass_minimum_eigenvalue > 0.0
        assert mode.zero_mode_count == 0
        assert mode.lowest_positive_eigenvalue > 0.0
        if item.case.rotary_inertia:
            assert mode.full_mass_spd is True
            assert mode.full_mass_minimum_eigenvalue is not None
            assert mode.full_mass_minimum_eigenvalue > 0.0
        else:
            assert mode.full_mass_spd is None
            assert mode.full_mass_minimum_eigenvalue is None


def test_shifted_legendre_bases_enforce_only_essential_conditions() -> None:
    endpoints = np.array([0.0, 1.0])
    for boundary_condition in ("HH", "CC", "CF"):
        w_values, _w_derivatives = audit.shifted_legendre_essential_basis(
            8, endpoints, boundary_condition, "w"
        )
        psi_values, _psi_derivatives = audit.shifted_legendre_essential_basis(
            8, endpoints, boundary_condition, "psi"
        )
        if boundary_condition == "HH":
            assert_array_equal(w_values, np.zeros_like(w_values))
            assert np.linalg.norm(psi_values) > 0.0
        elif boundary_condition == "CC":
            assert_array_equal(w_values, np.zeros_like(w_values))
            assert_array_equal(psi_values, np.zeros_like(psi_values))
        else:
            assert_array_equal(w_values[:, 0], np.zeros(8))
            assert_array_equal(psi_values[:, 0], np.zeros(8))
            assert np.linalg.norm(w_values[:, 1]) > 0.0
            assert np.linalg.norm(psi_values[:, 1]) > 0.0


def test_all_modes_satisfy_essential_and_hh_cf_natural_boundaries(
    full_audit: SimpleNamespace,
) -> None:
    natural_cases = 0
    for item in full_audit.audits:
        mode = item.converged_mode
        assert mode.essential_boundary_residual <= audit.ESSENTIAL_BOUNDARY_RESIDUAL_TOLERANCE
        states = audit.evaluate_ritz_mode(mode, np.array([0.0, 1.0]))
        denominator = max(float(np.linalg.norm(states)), sys.float_info.min)
        if item.case.boundary_condition == "HH":
            natural = float(np.linalg.norm(states[:, 3]) / denominator)
        elif item.case.boundary_condition == "CF":
            natural = float(np.linalg.norm(states[1, 2:]) / denominator)
        else:
            continue
        natural_cases += 1
        assert natural == pytest.approx(
            mode.natural_boundary_residual, rel=2.0e-13, abs=2.0e-15
        )
        assert natural <= audit.NATURAL_BOUNDARY_RESIDUAL_TOLERANCE
    assert natural_cases == 48


def test_exact_j0_static_condensation_has_zero_rotary_mass_and_full_residual(
    full_audit: SimpleNamespace,
) -> None:
    j0_items = [item for item in full_audit.audits if not item.case.rotary_inertia]
    assert len(j0_items) == 36
    for item in j0_items:
        mode = item.converged_mode
        assert item.case.properties.J == 0.0
        assert mode.static_condensation_used
        assert mode.mass_condition_number is None
        assert mode.full_mass_spd is None
        assert mode.psi_mass_zero_residual == 0.0
        assert mode.j0_mass_structure_pass
        assert mode.static_condensation_residual is not None
        assert mode.static_condensation_residual <= audit.ALGEBRAIC_RESIDUAL_TOLERANCE
        assert mode.recovered_full_mode_residual <= audit.ALGEBRAIC_RESIDUAL_TOLERANCE

    representative = next(
        item
        for item in j0_items
        if item.case.laminate_id == "0_deg"
        and item.case.boundary_condition == "CF"
        and item.case.a_over_h == 20
    )
    order = representative.converged_mode.order
    matrices = audit.assemble_ritz_matrices(
        representative.case.properties,
        representative.case.length,
        representative.case.boundary_condition,
        order,
    )
    assert matrices.rotary_mass_ratio == 0.0
    assert_array_equal(matrices.M_psipsi, np.zeros((order, order)))
    assert_array_equal(matrices.mass[order:, :], np.zeros((order, 2 * order)))
    assert np.linalg.norm(representative.converged_mode.coefficients[order:]) > 0.0


def test_nested_orders_converge_and_n16_guard_is_explicit(
    full_audit: SimpleNamespace,
) -> None:
    guarded: set[str] = set()
    for item in full_audit.audits:
        assert set(audit.RITZ_ORDERS).issubset(item.modes)
        frequencies = [item.modes[order].omega for order in audit.RITZ_ORDERS]
        assert all(math_value > 0.0 and np.isfinite(math_value) for math_value in frequencies)
        assert all(
            right <= left * (1.0 + 1.0e-10)
            for left, right in zip(frequencies, frequencies[1:])
        )
        assert item.used_convergence_relative <= audit.RITZ_CONVERGENCE_RELATIVE_TOLERANCE
        if audit.RITZ_GUARD_ORDER in item.modes:
            guarded.add(item.case.case_id)
            assert item.converged_order == audit.RITZ_GUARD_ORDER
            assert item.guard_relative is not None
            assert item.guard_relative <= audit.RITZ_CONVERGENCE_RELATIVE_TOLERANCE
        else:
            assert item.converged_order == 12
    assert guarded == {
        "90_deg__CC__a_h_100__with_RI",
        "90_deg__CF__a_h_100__with_RI",
    }


def test_representative_and_all_transfer_ritz_agreement_is_below_hard_limit(
    full_audit: SimpleNamespace,
) -> None:
    relatives = []
    for item in full_audit.audits:
        assert item.transfer.accepted
        relative = _relative_difference(
            item.converged_mode.omega, item.transfer.omega
        )
        relatives.append(relative)
        assert relative <= audit.TRANSFER_RITZ_RELATIVE_TOLERANCE
    assert max(relatives) <= 1.0e-8

    representatives = [
        item
        for item in full_audit.audits
        if item.case.laminate_id == "cross_ply_0_90_s"
        and item.case.a_over_h == 20
    ]
    assert len(representatives) == 6
    assert {
        (item.case.boundary_condition, item.case.rotary_inertia)
        for item in representatives
    } == {
        (boundary_condition, rotary_inertia)
        for boundary_condition in ("HH", "CC", "CF")
        for rotary_inertia in (False, True)
    }
    assert max(
        _relative_difference(item.converged_mode.omega, item.transfer.omega)
        for item in representatives
    ) <= 1.0e-8


def test_computed_frequencies_obey_rotary_inertia_monotonicity(
    full_audit: SimpleNamespace,
) -> None:
    rows = full_audit.rotary_rows
    assert len(rows) == 36
    assert sum(bool(row["source_pair_published"]) for row in rows) == 18
    assert sum(not bool(row["source_pair_published"]) for row in rows) == 18
    assert all(row["status"] == "PASS" for row in rows)
    assert all(row["transfer_monotonicity"] is True for row in rows)
    assert all(row["ritz_monotonicity"] is True for row in rows)
    assert all(float(row["transfer_without_minus_with"]) >= 0.0 for row in rows)
    assert all(float(row["ritz_without_minus_with"]) >= 0.0 for row in rows)
    assert len(full_audit.source_violations) == 11


def test_source_json_and_legacy_csv_values_and_hashes_are_preserved(
    full_audit: SimpleNamespace,
) -> None:
    assert _sha256(audit.SOURCE_CONTRACT_PATH) == SOURCE_SHA256
    assert _sha256(audit.LEGACY_TRANSFER_INVENTORY_PATH) == LEGACY_CSV_SHA256
    assert full_audit.metadata == {
        "path": "results/laminated_beams/reddy_symmetric_single_beam/bending_source_comparison.csv",
        "sha256": LEGACY_CSV_SHA256,
        "direct_row_count": 72,
        "published_direct_row_count": 54,
        "source_null_project_diagnostic_count": 18,
        "source_key_validation": "PASS",
    }

    for case in full_audit.cases:
        key = (
            case.laminate_id,
            case.boundary_condition,
            case.a_over_h,
            case.row_role,
        )
        row = full_audit.inventory[key]
        assert row["source_status"] == case.source_status
        assert float(row["K"]) == case.K
        assert row["K_provenance"] == case.K_provenance
        if case.source_value is None:
            assert row["source_value"] == ""
            assert row["printed_tolerance"] == ""
        else:
            assert float(row["source_value"]) == case.source_value
            assert float(row["printed_tolerance"]) == case.printed_tolerance

    assert full_audit.inventory[
        ("0_deg", "CF", 10, "fsdt_frequency_with_RI")
    ]["computed_value"] == "4.559813410416229"
    assert full_audit.inventory[
        ("cross_ply_0_90_s", "CF", 10, "fsdt_frequency_with_RI")
    ]["computed_value"] == "4.178374014380404"


def test_ast_proves_no_yartsev_or_circular_isotropic_imports() -> None:
    module_path = (
        ROOT
        / "scripts"
        / "analysis"
        / "laminated_beams"
        / "audit_reddy_table_4_3_3_discrepancies.py"
    )
    tree = ast.parse(module_path.read_text(encoding="utf-8"))
    imports: set[str] = set()
    imported_from: set[tuple[str, str]] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module)
            imported_from.update((node.module, alias.name) for alias in node.names)
    forbidden = ("yartsev", "circular", "variable_length_timoshenko")
    assert not {
        imported
        for imported in imports
        if any(fragment in imported.lower() for fragment in forbidden)
    }
    project_imports = {
        (module, name)
        for module, name in imported_from
        if module.startswith("scripts.")
    }
    assert project_imports == {("scripts.lib", "reddy_symmetric_laminated_beam")}


def test_import_is_side_effect_free_and_cli_smoke_writes_only_contract_outputs(
    tmp_path: Path,
) -> None:
    default_output = audit.DEFAULT_OUTPUT_DIR

    def snapshot(path: Path) -> dict[str, tuple[int, int]] | None:
        if not path.exists():
            return None
        return {
            str(item.relative_to(path)): (item.stat().st_size, item.stat().st_mtime_ns)
            for item in path.rglob("*")
            if item.is_file()
        }

    before = snapshot(default_output)
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    imported = subprocess.run(
        [
            sys.executable,
            "-B",
            "-c",
            (
                "import scripts.analysis.laminated_beams."
                "audit_reddy_table_4_3_3_discrepancies; print('IMPORTED')"
            ),
        ],
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert imported.returncode == 0, imported.stderr
    assert imported.stdout.strip() == "IMPORTED"
    assert snapshot(default_output) == before

    output_dir = tmp_path / "audit"
    completed = subprocess.run(
        [
            sys.executable,
            "-B",
            str(
                ROOT
                / "scripts"
                / "analysis"
                / "laminated_beams"
                / "audit_reddy_table_4_3_3_discrepancies.py"
            ),
            "--smoke",
            "--output-dir",
            str(output_dir),
        ],
        cwd=ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert {item.name for item in output_dir.iterdir()} == set(audit.REQUIRED_OUTPUTS)
    manifest = json.loads((output_dir / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["command_mode"] == "smoke"
    assert manifest["source_contract"]["sha256"] == SOURCE_SHA256
    assert manifest["immutable_transfer_inventory"]["sha256"] == LEGACY_CSV_SHA256
    assert manifest["immutable_transfer_inventory"]["sha256_after"] == LEGACY_CSV_SHA256
    assert manifest["immutable_transfer_inventory"]["preserved"] is True
    assert manifest["old_RLB0_outputs_modified"] is False
    assert manifest["summary"]["numerical"]["all_model_case_count"] == 24
    assert manifest["summary"]["numerical"]["published_direct_case_count"] == 18
    assert manifest["summary"]["numerical"]["rotary_pair_count"] == 12
    assert manifest["statuses"][audit.STATUS_MODEL] == "PASS"
    assert manifest["statuses"][audit.STATUS_RITZ] == "PASS"
    assert manifest["statuses"][audit.STATUS_RI] == "PASS"
    assert manifest["statuses"][audit.STATUS_OVERALL] == "PARTIAL_PASS"
    assert snapshot(default_output) == before
