"""Targeted gates for the RLB-2J pairwise stiffness-transfer maps.

The tests deliberately keep spectral work synthetic.  Physical calculations
are limited to the laminate constitutive reduction; no determinant root sweep
is started by this file.
"""

from __future__ import annotations

import ast
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.testing import assert_allclose
import pytest

from scripts.analysis.laminated_beams import (
    sweep_reddy_six_ply_pairwise_stiffness_transfer as sweep,
)


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = (
    ROOT
    / "scripts"
    / "analysis"
    / "laminated_beams"
    / "sweep_reddy_six_ply_pairwise_stiffness_transfer.py"
)
RESULT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_six_ply_pairwise_stiffness_transfer"
)


@pytest.fixture(scope="module")
def constitutive_data() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    return sweep.compute_constitutive_data()


def _synthetic_rows() -> list[dict[str, Any]]:
    """Return a complete spectrum that is an exact function of D/D0."""

    rows: list[dict[str, Any]] = []
    for configuration_id in sweep.CONFIGURATIONS:
        for index in range(-40, 41):
            xi = index / 50.0
            q = sweep.TRANSFER_LEVERS[configuration_id] * index
            multipliers = sweep.stiffness_multipliers(configuration_id, xi)
            physical_solve_id = (
                "SHARED_XI0_PHYSICAL_SOLVE"
                if index == 0
                else f"{configuration_id}__xi_{xi:+.2f}"
            )
            for position in range(1, sweep.K_GUARD + 1):
                # Depend only on q, so every exact matched-D comparison must
                # collapse without interpolation.
                Omega = 20.0 * position + 0.01 * q
                rows.append(
                    {
                        "row_id": (
                            f"{configuration_id}__{index:+03d}__p{position:02d}"
                        ),
                        "configuration_id": configuration_id,
                        "xi": xi,
                        "xi_index": index,
                        "s_C": multipliers[sweep.PAIR_CENTER],
                        "s_M": multipliers[sweep.PAIR_MIDDLE],
                        "s_O": multipliers[sweep.PAIR_OUTER],
                        "D_over_D0": float(Fraction(225 + q, 225)),
                        "D_key_q": q,
                        "grid_kind": "BASE",
                        "physical_solve_id": physical_solve_id,
                        "transaction_id": f"synthetic-{configuration_id}-{index}",
                        "continuation_leg": sweep.continuation_leg(xi),
                        "sorted_position": position,
                        "root_role": (
                            "PLOTTED"
                            if position <= sweep.K_PLOT
                            else "ROOT_9_GUARD"
                        ),
                        "guard_flag": position == sweep.K_GUARD,
                        "omega": Omega / sweep.OMEGA_TO_OMEGA_SCALE,
                        "Omega": Omega,
                        "Lambda": math.sqrt(Omega),
                        "predictor_Omega": "",
                        "predictor_used_as_final": False,
                        "locator_interval_left_Omega": Omega - 0.5,
                        "locator_interval_right_Omega": Omega + 0.5,
                        "root_interval_left_Omega": Omega - 1.0e-6,
                        "root_interval_right_Omega": Omega + 1.0e-6,
                        "detector_refiner_provenance": ["synthetic_refiner"],
                        "raw_determinant": 0.0,
                        "scaled_determinant": 0.0,
                        "raw_sigma_ratio": 0.0,
                        "scaled_sigma_ratio": 0.0,
                        "boundary_null_residual": 0.0,
                        "detected_nullity": 1,
                        "cluster_id": f"synthetic-p{position}",
                        "cluster_multiplicity": 1,
                        "cluster_total_nullity": 1,
                        "unresolved_candidates_below_root9": 0,
                        "search_right_Omega": 182.5 + 0.01 * q,
                        "root9_right_margin_Omega": 2.5,
                        "solve_mode": "SYNTHETIC",
                        "fallback_used": False,
                        "quality_status": "PASS",
                        "point_status": "PASS",
                        "shared_xi0_anchor_reused": (
                            index == 0
                            and configuration_id != sweep.CENTER_OUTER_TRANSFER
                        ),
                        "shared_xi0_source_configuration": (
                            sweep.CENTER_OUTER_TRANSFER if index == 0 else ""
                        ),
                        "is_canonical_plot_source": True,
                        "supersedes_row_id": "",
                        "repair_id": "",
                        "roots_above_9_computed": False,
                    }
                )
    return rows


def _write_required_completed_files(
    directory: Path, rows: list[dict[str, Any]]
) -> dict[str, bytes]:
    sweep._atomic_write_csv(
        directory / sweep.SPECTRUM_FILENAME,
        rows,
        sweep.SPECTRUM_FIELDS,
    )
    payloads: dict[str, bytes] = {}
    for filename in (
        sweep.SECTION_FILENAME,
        sweep.AUDIT_FILENAME,
        sweep.MATCHED_D_FILENAME,
        sweep.SLOPE_FILENAME,
    ):
        (directory / filename).write_text("marker\n", encoding="utf-8")
    for filename in (sweep.BENCHMARK_FILENAME, sweep.CHECKPOINT_FILENAME):
        (directory / filename).write_text("{}\n", encoding="utf-8")
    for filename in (sweep.XI_PLOT_FILENAME, sweep.MASTER_PLOT_FILENAME):
        (directory / filename).write_bytes(b"synthetic png marker")
    (directory / sweep.REPORT_FILENAME).write_text(
        "synthetic report\n", encoding="utf-8"
    )
    (directory / sweep.MANIFEST_FILENAME).write_text(
        json.dumps({"stage_id": sweep.STAGE_ID}) + "\n", encoding="utf-8"
    )
    for path in directory.iterdir():
        if path.is_file():
            payloads[path.name] = path.read_bytes()
    return payloads


def test_01_exact_geometry_material_and_six_ply_contract() -> None:
    assert (
        sweep.MU,
        sweep.TAU,
        sweep.BETA_DEG,
        sweep.L1,
        sweep.L2,
        sweep.L_TOTAL,
    ) == (0.0, 0.0, 30.0, 1.0, 1.0, 2.0)
    assert (sweep.WIDTH, sweep.THICKNESS, sweep.PLY_THICKNESS) == pytest.approx(
        (0.20, 0.05, 0.05 / 6.0)
    )
    assert sweep.K == pytest.approx(5.0 / 6.0)
    assert sweep.base_material_contract() == {
        "E1_0": 1.1,
        "E2_0": 0.9,
        "nu12_0": 0.3,
        "G12_0": 1.0 / 2.6,
        "G13_0": 1.0 / 2.6,
        "G23_0": 1.0 / 2.6,
        "rho_0": 1.0,
    }
    assert sweep.STACK_BOTTOM_TO_TOP == (
        "OUTER",
        "MIDDLE",
        "CENTER",
        "CENTER",
        "MIDDLE",
        "OUTER",
    )


def test_02_exact_xi_grid_and_positive_multiplier_range() -> None:
    assert_allclose(
        sweep.xi_grid(),
        np.arange(-40, 41, dtype=float) / 50.0,
        rtol=0.0,
        atol=0.0,
    )
    assert len(sweep.xi_grid()) == 81
    observed: list[float] = []
    for configuration_id in sweep.CONFIGURATIONS:
        for xi in sweep.xi_grid():
            multipliers = sweep.stiffness_multipliers(
                configuration_id, float(xi)
            )
            observed.extend(multipliers.values())
            assert sum(multipliers.values()) == pytest.approx(3.0)
    assert min(observed) == pytest.approx(0.2)
    assert max(observed) == pytest.approx(1.8)
    with pytest.raises(ValueError, match="exact"):
        sweep.xi_index(0.011)


def test_03_exact_transfer_vectors_metadata_and_D_levers() -> None:
    assert sweep.TRANSFER_VECTORS == {
        sweep.CENTER_MIDDLE_TRANSFER: (-1, 1, 0),
        sweep.MIDDLE_OUTER_TRANSFER: (0, -1, 1),
        sweep.CENTER_OUTER_TRANSFER: (-1, 0, 1),
    }
    assert sweep.TRANSFER_LEVERS == {
        sweep.CENTER_MIDDLE_TRANSFER: 1,
        sweep.MIDDLE_OUTER_TRANSFER: 2,
        sweep.CENTER_OUTER_TRANSFER: 3,
    }
    assert np.array_equal(
        np.asarray(sweep.TRANSFER_VECTORS[sweep.CENTER_OUTER_TRANSFER]),
        np.asarray(sweep.TRANSFER_VECTORS[sweep.CENTER_MIDDLE_TRANSFER])
        + np.asarray(sweep.TRANSFER_VECTORS[sweep.MIDDLE_OUTER_TRANSFER]),
    )
    assert [
        sum(
            weight * component
            for weight, component in zip(
                (1, 7, 19), sweep.TRANSFER_VECTORS[configuration_id], strict=True
            )
        )
        for configuration_id in sweep.CONFIGURATIONS
    ] == [6, 12, 18]
    assert [
        sweep.TRANSFER_METADATA[item]["fixed_pair"]
        for item in sweep.CONFIGURATIONS
    ] == ["OUTER", "CENTER", "MIDDLE"]


@pytest.mark.parametrize(
    ("configuration_id", "minus", "plus"),
    (
        (sweep.CENTER_MIDDLE_TRANSFER, Fraction(37, 45), Fraction(53, 45)),
        (sweep.MIDDLE_OUTER_TRANSFER, Fraction(29, 45), Fraction(61, 45)),
        (sweep.CENTER_OUTER_TRANSFER, Fraction(7, 15), Fraction(23, 15)),
    ),
)
def test_04_exact_D_formulas_and_endpoint_ranges(
    configuration_id: str, minus: Fraction, plus: Fraction
) -> None:
    assert sweep.exact_D_ratio(configuration_id, -0.8) == minus
    assert sweep.exact_D_ratio(configuration_id, 0.8) == plus
    for xi in sweep.xi_grid():
        index = sweep.xi_index(float(xi))
        lever = sweep.TRANSFER_LEVERS[configuration_id]
        assert sweep.D_key_q(configuration_id, float(xi)) == lever * index
        assert sweep.exact_D_ratio(configuration_id, float(xi)) == Fraction(
            225 + lever * index, 225
        )


@pytest.mark.parametrize("configuration_id", sweep.CONFIGURATIONS)
def test_05_exact_stack_angles_and_selective_material_scaling(
    configuration_id: str,
) -> None:
    beam, _coupled = sweep._physics_modules()
    section = sweep.build_six_ply_section(configuration_id, 0.8)
    assert len(section.laminate.plies) == 6
    assert tuple(ply.label for ply in section.laminate.plies) == (
        sweep.STACK_BOTTOM_TO_TOP
    )
    assert [ply.angle_deg for ply in section.laminate.plies] == [0.0] * 6
    assert_allclose(
        section.laminate.z_interfaces,
        sweep.THICKNESS
        * np.asarray([-1 / 2, -1 / 3, -1 / 6, 0, 1 / 6, 1 / 3, 1 / 2]),
        rtol=0.0,
        atol=1.0e-17,
    )
    base = sweep.build_scaled_material(1.0, "BASE")
    Q0 = beam.transformed_reduced_stiffness(base, 0.0)
    shear0 = beam.transformed_transverse_shear_stiffness(base, 0.0)
    for ply in section.laminate.plies:
        scale = section.multipliers[str(ply.label)]
        assert_allclose(
            beam.transformed_reduced_stiffness(ply.material, 0.0),
            scale * Q0,
            rtol=sweep.MATRIX_RELATIVE_TOLERANCE,
            atol=1.0e-16,
        )
        assert_allclose(
            beam.transformed_transverse_shear_stiffness(ply.material, 0.0),
            shear0,
            rtol=0.0,
            atol=0.0,
        )
        assert ply.material.nu12 == pytest.approx(sweep.NU12_0)
        assert ply.material.G13 == pytest.approx(sweep.G13_0)
        assert ply.material.G23 == pytest.approx(sweep.G23_0)
        assert ply.material.rho == pytest.approx(sweep.RHO_0)


def test_06_full_constitutive_gate_and_section_schema(
    constitutive_data: tuple[dict[str, Any], list[dict[str, Any]]],
) -> None:
    gate, rows = constitutive_data
    assert gate["status"] == "PASS"
    assert gate["vector_sum_identity"] is True
    assert gate["D_projections"] == {
        sweep.CENTER_MIDDLE_TRANSFER: 6,
        sweep.MIDDLE_OUTER_TRANSFER: 12,
        sweep.CENTER_OUTER_TRANSFER: 18,
    }
    assert len(rows) == 243
    assert len(
        {
            (str(row["configuration_id"]), int(row["xi_index"]))
            for row in rows
        }
    ) == 243
    maxima = gate["maximum_residuals"]
    assert maxima["sum_multiplier_absolute"] <= sweep.EXACT_IDENTITY_TOLERANCE
    assert maxima["A_matrix_invariance_relative"] <= sweep.MATRIX_RELATIVE_TOLERANCE
    assert maxima["D_matrix_formula_relative"] <= sweep.MATRIX_RELATIVE_TOLERANCE
    assert maxima["shear_matrix_invariance_relative"] <= sweep.MATRIX_RELATIVE_TOLERANCE
    assert maxima["A_beam_invariance_relative"] <= sweep.REDUCED_PROPERTY_TOLERANCE
    assert maxima["S_beam_invariance_relative"] <= sweep.REDUCED_PROPERTY_TOLERANCE
    assert maxima["mass_invariance_relative"] <= sweep.REDUCED_PROPERTY_TOLERANCE
    assert maxima["rotary_inertia_invariance_relative"] <= sweep.REDUCED_PROPERTY_TOLERANCE
    assert maxima["B_scaled_residual"] <= sweep.SYMMETRY_RELATIVE_TOLERANCE
    assert maxima["I1_scaled_residual"] <= sweep.SYMMETRY_RELATIVE_TOLERANCE
    assert all(row["constitutive_status"] == "PASS" for row in rows)


def test_07_matched_D_section_properties_are_identical(
    constitutive_data: tuple[dict[str, Any], list[dict[str, Any]]],
) -> None:
    rows = constitutive_data[1]
    by_q: dict[int, list[dict[str, Any]]] = {}
    for row in rows:
        by_q.setdefault(int(row["D_key_q"]), []).append(row)
    for matched in (items for items in by_q.values() if len(items) >= 2):
        reference = matched[0]
        for row in matched[1:]:
            for field in ("A_beam", "D_beam", "S_beam", "m", "J", "K", "width"):
                assert float(row[field]) == pytest.approx(
                    float(reference[field]), rel=sweep.REDUCED_PROPERTY_TOLERANCE
                )


def test_08_identical_arm_provider_reuses_one_arm_map() -> None:
    _provider, metadata = sweep.make_matrix_provider(
        sweep.CENTER_MIDDLE_TRANSFER, 0.2
    )
    assert metadata["identical_arms"] is True
    assert metadata["identical_arm_map_reused"] is True
    assert metadata["cached_vs_public_builder_max_abs"] <= 16.0 * np.finfo(float).eps


def test_09_frequency_map_policy_is_complete_and_local() -> None:
    contract = sweep.contract_payload()
    assert contract["frequency_map_policy"] == {
        "frequency_map_policy": "frequency-map-v1",
        "calculation_mode": "fast_plot",
        "spectrum_semantics": "sorted_positions",
        "sweep_parameter": "xi",
        "parameter_grid": "-0.80:0.02:0.80",
        "continuation_anchor": 0.0,
        "continuation_paths": ["0.00:-0.02:-0.80", "0.00:0.02:0.80"],
        "K_plot": 8,
        "K_guard": 9,
        "guard_root_role": "completeness_only",
        "neighbour_audit": "enabled",
        "local_repair_policy": "triggered_only",
        "strict_audit_default": False,
        "branch_tracking": False,
        "mac": False,
        "mode_shapes": False,
        "energy_analysis": False,
    }
    assert contract["root_contract"] == {
        "plotted_positions": list(range(1, 9)),
        "guard_position": 9,
        "guard_role": "completeness_only",
        "roots_above_9": "not_computed",
    }


def test_10_complete_synthetic_BASE_inventory_and_shared_anchor() -> None:
    rows = _synthetic_rows()
    audit = sweep.audit_spectrum_rows(rows)
    assert audit["status"] == "PASS"
    assert audit["base_group_count"] == 243
    assert audit["base_row_count"] == 2187
    assert audit["roots_above_guard_count"] == 0
    assert len({str(row["physical_solve_id"]) for row in rows}) == 241
    xi0 = [row for row in rows if int(row["xi_index"]) == 0]
    assert len(xi0) == 27
    assert {str(row["physical_solve_id"]) for row in xi0} == {
        "SHARED_XI0_PHYSICAL_SOLVE"
    }
    guards = [row for row in rows if int(row["sorted_position"]) == 9]
    assert len(guards) == 243
    assert all(row["root_role"] == "ROOT_9_GUARD" for row in guards)
    assert all(not bool(row["predictor_used_as_final"]) for row in rows)
    assert all(not bool(row["roots_above_9_computed"]) for row in rows)
    damaged = rows[:-1]
    assert sweep.audit_spectrum_rows(damaged)["status"] == "FAIL"


def test_11_exact_matched_D_key_sets_and_pairwise_rows(
    constitutive_data: tuple[dict[str, Any], list[dict[str, Any]]],
) -> None:
    sets = {
        configuration_id: {
            sweep.D_key_q(configuration_id, float(xi)) for xi in sweep.xi_grid()
        }
        for configuration_id in sweep.CONFIGURATIONS
    }
    cm = sets[sweep.CENTER_MIDDLE_TRANSFER]
    mo = sets[sweep.MIDDLE_OUTER_TRANSFER]
    co = sets[sweep.CENTER_OUTER_TRANSFER]
    assert len(cm & mo) == 41
    assert len(cm & co) == 27
    assert len(mo & co) == 27
    assert len(cm & mo & co) == 13
    assert len((cm & mo) | (cm & co) | (mo & co)) == 69

    matched = sweep.matched_D_rows(_synthetic_rows(), constitutive_data[1])
    assert len(matched) == 855
    assert len(
        {
            (
                str(row["left_configuration"]),
                str(row["right_configuration"]),
                int(row["D_key_q"]),
                int(row["sorted_position"]),
            )
            for row in matched
        }
    ) == 855
    assert len({int(row["D_key_q"]) for row in matched}) == 69
    assert len(
        {
            int(row["D_key_q"])
            for row in matched
            if bool(row["common_to_all_three"])
        }
    ) == 13
    assert all(not bool(row["interpolation_used"]) for row in matched)
    assert all(float(row["relative_Omega_difference"]) == 0.0 for row in matched)
    assert all(row["status"] == "PASS" for row in matched)


def test_12_initial_slopes_recover_lever_ratio_for_synthetic_master_curve() -> None:
    checks = sweep.initial_slope_rows(_synthetic_rows())
    assert len(checks) == 8
    assert all(row["expected_lever_ratio"] == "1:2:3" for row in checks)
    assert all(row["hard_spectral_gate"] is False for row in checks)
    assert all(
        float(row["common_D_secant_relative_spread"]) <= 2.0e-12
        for row in checks
    )
    assert all(
        row["status"] in {"PASS", "INSENSITIVE_AT_BASELINE"}
        for row in checks
    )


def test_13_neighbour_audit_is_trigger_only_and_never_smoothing() -> None:
    rows = _synthetic_rows()
    target = next(
        row
        for row in rows
        if row["configuration_id"] == sweep.CENTER_MIDDLE_TRANSFER
        and int(row["xi_index"]) == 3
        and int(row["sorted_position"]) == 2
    )
    target["Lambda"] = float(target["Lambda"]) * 1.2
    audit = sweep.neighbour_audit_rows(rows)
    assert len(audit) == 3 * 79 * 8
    assert any(bool(row["flagged"]) for row in audit)
    assert all(not bool(row["smoothing_applied"]) for row in audit)
    assert sweep.flagged_repair_points(audit)


def test_14_root9_is_excluded_and_unresolved_plot_value_remains_a_gap() -> None:
    rows = _synthetic_rows()
    target = next(
        row
        for row in rows
        if row["configuration_id"] == sweep.CENTER_OUTER_TRANSFER
        and int(row["xi_index"]) == 7
        and int(row["sorted_position"]) == 3
    )
    target["Lambda"] = math.nan
    target["point_status"] = "UNRESOLVED_AFTER_LOCAL_REPAIR"
    audit = sweep.audit_plot_rows(rows)
    assert audit["status"] == "PASS"
    assert audit["row_count"] == 1944
    assert audit["root9_plotted"] is False
    assert math.isnan(float(target["Lambda"]))


def test_15_plot_only_reads_spectrum_CSV_and_runs_no_physics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rows = _synthetic_rows()
    sweep._atomic_write_csv(
        tmp_path / sweep.SPECTRUM_FILENAME, rows, sweep.SPECTRUM_FIELDS
    )
    source_hash = sweep._sha256(tmp_path / sweep.SPECTRUM_FILENAME)

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("plot-only entered a physics/root path")

    monkeypatch.setattr(sweep, "_physics_modules", forbidden)
    monkeypatch.setattr(sweep, "_base", forbidden)
    monkeypatch.setattr(sweep, "build_six_ply_section", forbidden)
    monkeypatch.setattr(sweep, "compute_constitutive_data", forbidden)
    monkeypatch.setattr(sweep, "solve_point", forbidden)
    result = sweep.create_plots_from_csv(tmp_path)
    assert result["root_calculation_count"] == 0
    assert result["xi_panel_count"] == 3
    assert result["spectral_lines_per_xi_panel"] == 8
    assert result["master_solid_lines"] == 8
    assert result["root9_plotted"] is False
    assert result["interpolation_used"] is False
    assert result["smoothing_applied"] is False
    assert sweep._sha256(tmp_path / sweep.SPECTRUM_FILENAME) == source_hash
    assert (tmp_path / sweep.XI_PLOT_FILENAME).stat().st_size > 0
    assert (tmp_path / sweep.MASTER_PLOT_FILENAME).stat().st_size > 0


def test_16_completed_missing_only_performs_no_calculation_or_rewrite(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    before = _write_required_completed_files(tmp_path, _synthetic_rows())

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("completed missing-only entered a calculation path")

    monkeypatch.setattr(sweep, "prepare_constitutive_outputs", forbidden)
    monkeypatch.setattr(sweep, "run_benchmarks", forbidden)
    monkeypatch.setattr(sweep, "complete_missing_points", forbidden)
    monkeypatch.setattr(sweep, "finalize_outputs", forbidden)
    monkeypatch.setattr(
        sweep, "_refresh_completed_metadata", lambda _target, manifest: manifest
    )
    before_root_count = sweep.ROOT_CALCULATION_COUNT
    result = sweep.run_workflow(tmp_path, missing_only=True)
    after = {
        path.name: path.read_bytes() for path in tmp_path.iterdir() if path.is_file()
    }
    assert result["invocation_root_calculation_count"] == 0
    assert result["completed_missing_only_numerical_outputs_unchanged"] is True
    assert sweep.ROOT_CALCULATION_COUNT == before_root_count
    assert after == before


def test_17_no_forbidden_imports_local_solver_or_smoothing() -> None:
    source = SCRIPT_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    imports: list[str] = []
    local_functions: set[str] = set()
    calls: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imports.append(node.module or "")
        elif isinstance(node, ast.FunctionDef):
            local_functions.add(node.name)
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                calls.add(node.func.id)
            elif isinstance(node.func, ast.Attribute):
                calls.add(node.func.attr)
    joined = "\n".join(imports).lower()
    for forbidden in (
        "spectral_sweep_runner",
        "ritz",
        "fem",
        "yartsev",
        "shellbuckling",
        "circular",
        "analytic_branch_tracking",
    ):
        assert forbidden not in joined
    assert not any("characteristic_determinant" in name for name in local_functions)
    assert not ({"interp", "interp1d", "UnivariateSpline", "savgol_filter"} & calls)


@pytest.mark.skipif(
    not (RESULT_DIR / sweep.MANIFEST_FILENAME).is_file(),
    reason="RLB-2J final manifest not generated",
)
def test_18_generated_outputs_are_complete_and_consistent() -> None:
    manifest = json.loads(
        (RESULT_DIR / sweep.MANIFEST_FILENAME).read_text(encoding="utf-8")
    )
    spectrum = sweep._read_csv(RESULT_DIR / sweep.SPECTRUM_FILENAME)
    sections = sweep._read_csv(RESULT_DIR / sweep.SECTION_FILENAME)
    matched = sweep._read_csv(RESULT_DIR / sweep.MATCHED_D_FILENAME)
    slopes = sweep._read_csv(RESULT_DIR / sweep.SLOPE_FILENAME)
    assert manifest["stage_id"] == "RLB-2J"
    assert manifest["scientific_status"] == "PASS"
    assert sweep.audit_spectrum_rows(spectrum)["status"] == "PASS"
    assert sweep.audit_plot_rows(spectrum)["status"] == "PASS"
    assert len(sections) == 243
    assert len(matched) == 855
    assert len(slopes) == 8
    assert manifest["counts"]["base_groups"] == 243
    assert manifest["counts"]["base_rows"] == 2187
    assert manifest["counts"]["local_refinement_rows"] == 225
    assert manifest["counts"]["root9_guards"] == 243
    assert manifest["counts"]["unique_physical_base_solves"] == 241
    assert manifest["counts"]["unique_physical_base_solve_contract"] == 241
    assert manifest["neighbour_audit"]["row_count"] == 1896
    assert manifest["neighbour_audit"]["repair_count"] == 25
    assert manifest["neighbour_audit"]["unresolved_point_count"] == 0
    assert manifest["matched_D_summary"]["unique_matched_D_values"] == 69
    assert manifest["matched_D_summary"]["common_all_three_D_values"] == 13
    assert manifest["matched_D_summary"]["interpolation_used"] is False
    assert manifest["four_vs_six_ply_equivalence"]["status"] == "PASS"
    assert manifest["four_vs_six_ply_equivalence"]["row_count"] == 135
    assert manifest["plots"]["root9_plotted"] is False
    assert manifest["exclusions_confirmed"]["roots_above_9_computed"] is False
    assert manifest["exclusions_confirmed"]["branch_tracking"] is False
    runtime = manifest["runtime"]
    assert runtime["measurement_scope_version"] == 2
    assert runtime["campaign_total_wall_time_available"] is False
    assert runtime["campaign_total_wall_time_seconds"] is None
    assert runtime["recorded_base_point_wall_time_sum_seconds"] > 0.0
    assert runtime["determinant_evaluations_recorded_lower_bound"] > 0
    assert runtime["accepted_local_refinement_group_count"] == 25
    assert runtime["local_repair_determinant_evaluations"] is None
    assert manifest["metadata_refresh_root_calculation_count"] == 0
    assert manifest["metadata_refresh_numerical_outputs_changed"] is False
    assert manifest["analysis_script_sha256"] == sweep._sha256(SCRIPT_PATH)
    assert manifest["production_physics_hashes_match_frozen"] is True
    assert manifest["predecessor_result_trees_preserved"] is True
    required = {
        sweep.SPECTRUM_FILENAME,
        sweep.SECTION_FILENAME,
        sweep.AUDIT_FILENAME,
        sweep.MATCHED_D_FILENAME,
        sweep.SLOPE_FILENAME,
        sweep.BENCHMARK_FILENAME,
        sweep.CHECKPOINT_FILENAME,
        sweep.XI_PLOT_FILENAME,
        sweep.MASTER_PLOT_FILENAME,
        sweep.REPORT_FILENAME,
        sweep.MANIFEST_FILENAME,
    }
    assert required <= {path.name for path in RESULT_DIR.iterdir() if path.is_file()}
    for name, expected in manifest["output_hashes"].items():
        assert sweep._sha256(RESULT_DIR / name) == expected
