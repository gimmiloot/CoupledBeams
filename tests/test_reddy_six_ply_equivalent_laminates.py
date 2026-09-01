from __future__ import annotations

import ast
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.testing import assert_allclose
import pytest

from scripts.analysis.laminated_beams import (
    analyze_reddy_six_ply_equivalent_laminates as analysis,
)


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = (
    ROOT
    / "scripts"
    / "analysis"
    / "laminated_beams"
    / "analyze_reddy_six_ply_equivalent_laminates.py"
)
RESULT_DIR = (
    ROOT / "results" / "laminated_beams" / "reddy_six_ply_equivalent_laminates"
)


@pytest.fixture(scope="module")
def constitutive_data() -> tuple[dict[str, Any], list[dict[str, Any]], dict[float, Any]]:
    return analysis.compute_constitutive_data()


@pytest.fixture(scope="module")
def response_data(
    constitutive_data: tuple[dict[str, Any], list[dict[str, Any]], dict[float, Any]],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    return analysis.compute_ply_response_data(constitutive_data[2])


def test_01_exact_geometry_material_and_six_ply_contract() -> None:
    assert (analysis.MU, analysis.TAU, analysis.L1, analysis.L2, analysis.L_TOTAL) == (
        0.0,
        0.0,
        1.0,
        1.0,
        2.0,
    )
    assert (analysis.WIDTH, analysis.THICKNESS, analysis.PLY_THICKNESS) == pytest.approx(
        (0.20, 0.05, 0.05 / 6.0)
    )
    assert analysis.K == pytest.approx(5.0 / 6.0)
    assert analysis.base_material_contract() == {
        "E1_0": 1.1,
        "E2_0": 0.9,
        "nu12_0": 0.3,
        "G12_0": 1.0 / 2.6,
        "G13_0": 1.0 / 2.6,
        "G23_0": 1.0 / 2.6,
        "rho_0": 1.0,
    }
    assert analysis.STACK_BOTTOM_TO_TOP == (
        "OUTER",
        "MIDDLE",
        "CENTER",
        "CENTER",
        "MIDDLE",
        "OUTER",
    )


def test_02_exact_zeta_grid_and_positive_multipliers() -> None:
    assert_allclose(
        analysis.zeta_grid(), np.arange(-25, 26, dtype=float) / 100.0, rtol=0.0, atol=0.0
    )
    assert len(analysis.zeta_grid()) == 51
    for zeta in analysis.zeta_grid():
        multipliers = analysis.stiffness_multipliers(float(zeta))
        assert multipliers["CENTER"] == pytest.approx(1.0 + 2.0 * zeta)
        assert multipliers["MIDDLE"] == pytest.approx(1.0 - 3.0 * zeta)
        assert multipliers["OUTER"] == pytest.approx(1.0 + zeta)
        assert min(multipliers.values()) > 0.0


def test_03_exact_neutral_direction_and_pair_weights() -> None:
    assert analysis.PAIR_A_WEIGHTS == {
        "CENTER": 1.0 / 3.0,
        "MIDDLE": 1.0 / 3.0,
        "OUTER": 1.0 / 3.0,
    }
    assert analysis.PAIR_D_WEIGHTS == {
        "CENTER": 1.0 / 27.0,
        "MIDDLE": 7.0 / 27.0,
        "OUTER": 19.0 / 27.0,
    }
    for zeta in analysis.zeta_grid():
        s = analysis.stiffness_multipliers(float(zeta))
        assert s["CENTER"] + s["MIDDLE"] + s["OUTER"] == pytest.approx(3.0)
        assert (
            s["CENTER"] + 7.0 * s["MIDDLE"] + 19.0 * s["OUTER"]
        ) == pytest.approx(27.0)
    assert 2 - 3 + 1 == 0
    assert 2 - 3 * 7 + 19 == 0


@pytest.mark.parametrize("zeta", [-0.25, 0.0, 0.25])
def test_04_exact_stack_interfaces_angles_and_material_scaling(zeta: float) -> None:
    beam, _coupled = analysis._physics_modules()
    section = analysis.build_six_ply_section(zeta)
    assert len(section.laminate.plies) == 6
    assert tuple(ply.label for ply in section.laminate.plies) == analysis.STACK_BOTTOM_TO_TOP
    assert all(ply.angle_deg == 0.0 for ply in section.laminate.plies)
    assert_allclose(
        section.laminate.z_interfaces,
        analysis.THICKNESS * np.asarray([-1 / 2, -1 / 3, -1 / 6, 0, 1 / 6, 1 / 3, 1 / 2]),
        rtol=0.0,
        atol=1.0e-17,
    )
    base = analysis.build_scaled_material(1.0, "base")
    Q0 = beam.transformed_reduced_stiffness(base, 0.0)
    shear0 = beam.transformed_transverse_shear_stiffness(base, 0.0)
    for ply in section.laminate.plies:
        scale = section.multipliers[ply.label]
        assert_allclose(
            beam.transformed_reduced_stiffness(ply.material, 0.0),
            scale * Q0,
            rtol=analysis.MATRIX_RELATIVE_TOLERANCE,
            atol=1.0e-16,
        )
        assert_allclose(
            beam.transformed_transverse_shear_stiffness(ply.material, 0.0),
            shear0,
            rtol=0.0,
            atol=0.0,
        )
        assert ply.material.nu12 == pytest.approx(analysis.NU12_0)
        assert ply.material.G13 == pytest.approx(analysis.G13_0)
        assert ply.material.G23 == pytest.approx(analysis.G23_0)
        assert ply.material.rho == pytest.approx(analysis.RHO_0)


def test_05_constitutive_gate_full_grid_passes(
    constitutive_data: tuple[dict[str, Any], list[dict[str, Any]], dict[float, Any]],
) -> None:
    gate, rows, sections = constitutive_data
    assert gate["status"] == "PASS"
    assert len(rows) == len(sections) == 51
    assert gate["neutral_direction"] == [2, -3, 1]
    assert gate["maxima"]["A_matrix_relative"] <= analysis.MATRIX_RELATIVE_TOLERANCE
    assert gate["maxima"]["D_matrix_relative"] <= analysis.MATRIX_RELATIVE_TOLERANCE
    assert gate["maxima"]["shear_matrix_relative"] <= analysis.MATRIX_RELATIVE_TOLERANCE
    assert gate["maxima"]["B_scaled_residual"] <= analysis.SYMMETRY_RELATIVE_TOLERANCE
    assert gate["maxima"]["I1_scaled_residual"] <= analysis.SYMMETRY_RELATIVE_TOLERANCE


def test_06_all_production_beam_fields_are_invariant(
    constitutive_data: tuple[dict[str, Any], list[dict[str, Any]], dict[float, Any]],
) -> None:
    sections = constitutive_data[2]
    baseline = sections[0.0].properties
    for section in sections.values():
        for field in ("A", "D", "S", "m", "J", "K", "width"):
            assert getattr(section.properties, field) == pytest.approx(
                getattr(baseline, field), rel=analysis.REDUCED_PROPERTY_TOLERANCE
            )
    assert baseline.A == pytest.approx(0.011)
    assert baseline.D == pytest.approx(2.291666666666667e-6)
    assert baseline.S == pytest.approx(0.003205128205128205)
    assert baseline.m == pytest.approx(0.01)
    assert baseline.J == pytest.approx(2.083333333333334e-6)


def test_07_state_matrices_and_boundary_matrices_coincide(
    constitutive_data: tuple[dict[str, Any], list[dict[str, Any]], dict[float, Any]],
) -> None:
    gate, rows = analysis.boundary_matrix_equivalence_rows(constitutive_data[2])
    assert gate["status"] == "PASS"
    assert len(rows) == 18
    assert gate["maximum_state_matrix_relative_difference"] <= analysis.MATRIX_RELATIVE_TOLERANCE
    assert gate["maximum_boundary_matrix_relative_difference"] <= analysis.MATRIX_RELATIVE_TOLERANCE


def test_08_full_width_resultant_and_energy_convention(
    constitutive_data: tuple[dict[str, Any], list[dict[str, Any]], dict[float, Any]],
) -> None:
    section = constitutive_data[2][0.0]
    N, M, epsilon0, kappa = analysis._load_generalized_fields(
        section, "UNIT_AXIAL_RESULTANT"
    )
    assert_allclose(analysis.WIDTH * section.laminate.A @ epsilon0, N, atol=1.0e-13)
    assert_allclose(M, np.zeros(3), atol=0.0)
    assert_allclose(kappa, np.zeros(3), atol=0.0)
    N, M, epsilon0, kappa = analysis._load_generalized_fields(
        section, "UNIT_BENDING_MOMENT"
    )
    assert_allclose(analysis.WIDTH * section.laminate.D @ kappa, M, atol=1.0e-13)
    assert_allclose(N, np.zeros(3), atol=0.0)
    assert_allclose(epsilon0, np.zeros(3), atol=0.0)


def test_09_resultant_and_total_energy_identities_pass_all_102_load_states(
    response_data: tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]],
) -> None:
    gate, rows, summaries = response_data
    assert gate["status"] == "PASS"
    assert len(rows) == 612
    assert len(summaries) == 102
    assert gate["maximum_residuals"]["resultant_relative_residual"] <= analysis.RESULTANT_RECONSTRUCTION_TOLERANCE
    assert gate["maximum_residuals"]["energy_identity_relative_residual"] <= analysis.ENERGY_IDENTITY_TOLERANCE
    for summary in summaries:
        if summary["load_case"] == "UNIT_AXIAL_RESULTANT":
            assert_allclose(summary["M_reconstructed"], np.zeros(3), atol=1.0e-12)
        else:
            assert_allclose(summary["N_reconstructed"], np.zeros(3), atol=1.0e-12)


def test_10_analytical_pair_energy_fractions_and_sums(
    response_data: tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]],
) -> None:
    gate, rows, summaries = response_data
    assert gate["maximum_residuals"]["pair_fraction_residual"] <= analysis.ENERGY_FRACTION_TOLERANCE
    assert gate["maximum_residuals"]["pair_fraction_sum_residual"] <= analysis.ENERGY_FRACTION_TOLERANCE
    for summary in summaries:
        assert math.fsum(summary["pair_energy_fractions"].values()) == pytest.approx(1.0)
    for zeta, expected_axial, expected_bending in (
        (-0.25, (1 / 6, 7 / 12, 1 / 4), (0.5 / 27, 12.25 / 27, 14.25 / 27)),
        (0.0, (1 / 3, 1 / 3, 1 / 3), (1 / 27, 7 / 27, 19 / 27)),
        (0.25, (1 / 2, 1 / 12, 5 / 12), (1.5 / 27, 1.75 / 27, 23.75 / 27)),
    ):
        axial = next(
            item for item in summaries if item["zeta"] == zeta and item["load_case"] == "UNIT_AXIAL_RESULTANT"
        )
        bending = next(
            item for item in summaries if item["zeta"] == zeta and item["load_case"] == "UNIT_BENDING_MOMENT"
        )
        assert_allclose(
            [axial["pair_energy_fractions"][pair] for pair in analysis.PAIR_ORDER],
            expected_axial,
            rtol=0.0,
            atol=1.0e-14,
        )
        assert_allclose(
            [bending["pair_energy_fractions"][pair] for pair in analysis.PAIR_ORDER],
            expected_bending,
            rtol=0.0,
            atol=1.0e-14,
        )


def test_11_axial_and_bending_hotspot_transitions_and_ties(
    response_data: tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]],
) -> None:
    gate, rows, summaries = response_data
    assert gate["axial_hotspot_transition"]["status"] == "PASS"
    assert gate["bending_hotspot_transition"] == {
        "analytical_parameter": -1.0 / 9.0,
        "grid_bracket": [-0.12, -0.11],
        "lower_side_pair": "MIDDLE",
        "upper_side_pair": "OUTER",
        "status": "PASS",
    }
    axial_zero = next(
        item for item in summaries if item["zeta"] == 0.0 and item["load_case"] == "UNIT_AXIAL_RESULTANT"
    )
    assert axial_zero["hotspot_pairs"] == list(analysis.PAIR_ORDER)
    assert axial_zero["hotspot_tie_count"] == 6
    bend_minus = next(
        item for item in summaries if item["zeta"] == -0.12 and item["load_case"] == "UNIT_BENDING_MOMENT"
    )
    bend_plus = next(
        item for item in summaries if item["zeta"] == -0.11 and item["load_case"] == "UNIT_BENDING_MOMENT"
    )
    assert bend_minus["hotspot_pairs"] == ["MIDDLE"]
    assert bend_plus["hotspot_pairs"] == ["OUTER"]
    assert all(sum(bool(row["canonical_hotspot_flag"]) for row in rows if row["zeta"] == z and row["load_case"] == load) == 1 for z in analysis.zeta_grid() for load in ("UNIT_AXIAL_RESULTANT", "UNIT_BENDING_MOMENT"))


def test_12_piecewise_stress_jumps_are_retained_without_smoothing(
    response_data: tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]],
) -> None:
    summary = analysis._stress_jump_summary(response_data[1])
    assert summary["piecewise_segments_retained"] is True
    assert summary["smoothing_applied"] is False
    assert summary["interface_record_count"] == 15
    assert summary["nonzero_jump_count"] > 0
    assert summary["maximum_absolute_jump"] > 0.0


def test_13_plot_only_reads_csv_and_performs_zero_calculations(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    response_data: tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]],
) -> None:
    analysis._atomic_write_csv(tmp_path / analysis.PLY_FILENAME, response_data[1], analysis.PLY_FIELDS)
    before_hash = analysis._sha256(tmp_path / analysis.PLY_FILENAME)

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("plot-only entered a physics or root-calculation path")

    monkeypatch.setattr(analysis, "build_six_ply_section", forbidden)
    monkeypatch.setattr(analysis, "compute_constitutive_data", forbidden)
    monkeypatch.setattr(analysis, "compute_ply_response_data", forbidden)
    monkeypatch.setattr(analysis, "spectral_spot_check", forbidden)
    monkeypatch.setattr(analysis, "boundary_matrix_equivalence_rows", forbidden)
    monkeypatch.setattr(analysis, "_physics_modules", forbidden)
    result = analysis.plot_only_workflow(tmp_path)
    assert result["root_calculation_count"] == 0
    assert result["laminate_recalculation_count"] == 0
    assert result["matrix_assembly_count"] == 0
    assert result["determinant_evaluation_count"] == 0
    assert result["SVD_evaluation_count"] == 0
    assert analysis._sha256(tmp_path / analysis.PLY_FILENAME) == before_hash
    for name in (
        analysis.ENERGY_PLOT_FILENAME,
        analysis.STRESS_PLOT_FILENAME,
        analysis.PROFILE_PLOT_FILENAME,
    ):
        assert (tmp_path / name).stat().st_size > 0


def test_14_manifest_only_is_zero_calculation() -> None:
    payload = analysis.manifest_only()
    assert payload["root_calculation_count"] == 0
    assert payload["laminate_calculation_count"] == 0
    assert payload["matrix_assembly_count"] == 0
    assert payload["contract"]["spectral_spot_check"]["root_calculation_case_count"] == 6
    assert payload["contract"]["spectral_spot_check"]["roots_10_plus"] is False


def test_15_no_forbidden_imports_solver_or_smoothing() -> None:
    source = SCRIPT_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    imports: list[str] = []
    locally_defined: set[str] = set()
    calls: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module or "")
        elif isinstance(node, ast.FunctionDef):
            locally_defined.add(node.name)
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                calls.add(node.func.id)
            elif isinstance(node.func, ast.Attribute):
                calls.add(node.func.attr)
    joined = "\n".join(imports)
    for forbidden in (
        "spectral_sweep_runner",
        "ritz",
        "fem",
        "Yartsev",
        "yartsev",
        "ShellBuckling",
        "circular",
        "analytic_branch_tracking",
    ):
        assert forbidden not in joined
    assert not any("characteristic_determinant" in name for name in locally_defined)
    assert not ({"interp", "interp1d", "UnivariateSpline", "savgol_filter"} & calls)


@pytest.mark.skipif(
    not (RESULT_DIR / analysis.MANIFEST_FILENAME).is_file(),
    reason="RLB-2I final manifest not generated",
)
def test_16_generated_outputs_are_complete_and_consistent() -> None:
    manifest = json.loads((RESULT_DIR / analysis.MANIFEST_FILENAME).read_text(encoding="utf-8"))
    section_rows = analysis._read_csv(RESULT_DIR / analysis.SECTION_FILENAME)
    ply_rows = analysis._read_csv(RESULT_DIR / analysis.PLY_FILENAME)
    spectral_rows = analysis._read_csv(RESULT_DIR / analysis.SPECTRAL_FILENAME)
    boundary_rows = analysis._read_csv(RESULT_DIR / analysis.BOUNDARY_FILENAME)
    assert manifest["stage_id"] == "RLB-2I"
    assert manifest["scientific_status"] == "PASS"
    assert set(manifest["status_gates"].values()) == {"PASS"}
    assert len(section_rows) == 51
    assert analysis._audit_ply_csv_rows(ply_rows)["status"] == "PASS"
    assert len(spectral_rows) == 54
    assert len(boundary_rows) == 18
    keys = {
        (float(row["beta_deg"]), float(row["zeta"]), int(row["sorted_position"]))
        for row in spectral_rows
    }
    assert len(keys) == 54
    assert {int(row["sorted_position"]) for row in spectral_rows} == set(range(1, 10))
    assert all(
        (row["root_role"] == "ROOT_9_GUARD") == (int(row["sorted_position"]) == 9)
        for row in spectral_rows
    )
    assert all(row["roots_above_9_computed"] == "false" for row in spectral_rows)
    assert manifest["spectral_spot_check"]["maximum_relative_frequency_difference"] <= analysis.SPECTRAL_RELATIVE_TOLERANCE
    assert manifest["spectral_spot_check"]["maximum_scaled_sigma_ratio"] <= analysis.ROOT_SINGULAR_RATIO_TOLERANCE
    assert manifest["spectral_spot_check"]["maximum_boundary_null_residual"] <= analysis.BOUNDARY_RESIDUAL_TOLERANCE
    assert manifest["explicit_confirmations"]["wide_51_point_frequency_sweep_run"] is False
    assert manifest["explicit_confirmations"]["roots_above_9_computed"] is False
    assert manifest["production_physics_preserved"] is True
    assert manifest["predecessor_result_trees_preserved"] is True
    required = {
        analysis.SECTION_FILENAME,
        analysis.PLY_FILENAME,
        analysis.SPECTRAL_FILENAME,
        analysis.BOUNDARY_FILENAME,
        analysis.ENERGY_PLOT_FILENAME,
        analysis.STRESS_PLOT_FILENAME,
        analysis.PROFILE_PLOT_FILENAME,
        analysis.REPORT_FILENAME,
        analysis.MANIFEST_FILENAME,
    }
    assert required <= {path.name for path in RESULT_DIR.iterdir() if path.is_file()}
    for name, expected in manifest["output_hashes"].items():
        assert analysis._sha256(RESULT_DIR / name) == expected

