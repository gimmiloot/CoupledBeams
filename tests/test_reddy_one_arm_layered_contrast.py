"""Targeted gates for the RLB-2F signed one-arm contrast map.

The unit tests use synthetic spectra for orchestration and plotting checks.
They never launch the production root search.  The final test audits generated
artifacts only when the RLB-2F result tree already exists.
"""

from __future__ import annotations

import ast
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
from numpy.testing import assert_allclose
import pytest

from scripts.analysis.laminated_beams import (
    sweep_reddy_one_arm_layered_contrast as sweep,
)


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = (
    ROOT
    / "scripts"
    / "analysis"
    / "laminated_beams"
    / "sweep_reddy_one_arm_layered_contrast.py"
)
RESULT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_one_arm_layered_contrast_sweep"
)


def _synthetic_rows() -> list[dict[str, Any]]:
    """Return a complete, quality-passing 81 by 9 BASE inventory."""

    rows: list[dict[str, Any]] = []
    for xi_value in sweep.xi_grid():
        xi = float(xi_value)
        leg = "ANCHOR" if xi == 0.0 else ("NEGATIVE" if xi < 0.0 else "POSITIVE")
        for position in range(1, sweep.K_GUARD + 1):
            Omega = 10.0 * position + 0.2 * xi
            rows.append(
                {
                    "row_id": f"{sweep.CONFIGURATION_ID}__{xi:+.6f}__p{position:02d}",
                    "configuration_id": sweep.CONFIGURATION_ID,
                    "xi": xi,
                    "grid_kind": "BASE",
                    "continuation_leg": leg,
                    "solve_id": "synthetic",
                    "transaction_id": f"synthetic_{xi:+.6f}",
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
                    "locator_interval_left_Omega": Omega - 1.0,
                    "locator_interval_right_Omega": Omega + 1.0,
                    "root_interval_left_Omega": Omega - 0.1,
                    "root_interval_right_Omega": Omega + 0.1,
                    "detector_refiner_provenance": ["synthetic"],
                    "raw_determinant": 0.0,
                    "scaled_determinant": 0.0,
                    "raw_sigma_ratio": 0.0,
                    "scaled_sigma_ratio": 0.0,
                    "boundary_null_residual": 0.0,
                    "detected_nullity": 1,
                    "unresolved_candidates_below_root9": 0,
                    "search_right_Omega": 95.0 + 0.2 * xi,
                    "root9_right_margin_Omega": 5.0,
                    "solve_mode": "SYNTHETIC",
                    "fallback_used": False,
                    "quality_status": "PASS",
                    "point_status": "PASS",
                    "is_canonical_plot_source": True,
                    "supersedes_row_id": "",
                    "repair_id": "",
                    "roots_above_9_computed": False,
                }
            )
    return rows


def _synthetic_solution(
    rows: list[dict[str, Any]],
    xi: float,
    *,
    Omega_shift: float = 0.0,
    swapped: bool = False,
) -> sweep.PointSolution:
    selected: list[dict[str, Any]] = []
    for source in rows:
        if (
            str(source["grid_kind"]) == "BASE"
            and round(float(source["xi"]), 10) == round(float(xi), 10)
        ):
            row = dict(source)
            position = int(row["sorted_position"])
            Omega = float(row["Omega"]) + Omega_shift
            row.update(
                {
                    "row_id": (
                        f"{sweep.CONFIGURATION_ID}__{float(xi):+.6f}__"
                        f"synthetic_refinement__p{position:02d}"
                    ),
                    "grid_kind": "LOCAL_REFINEMENT",
                    "Omega": Omega,
                    "omega": Omega / sweep.OMEGA_TO_OMEGA_SCALE,
                    "Lambda": math.sqrt(Omega),
                    "repair_id": "synthetic_refinement",
                }
            )
            selected.append(row)
    selected.sort(key=lambda row: int(row["sorted_position"]))
    assert len(selected) == sweep.K_GUARD
    return sweep.PointSolution(
        xi=float(xi),
        rows=tuple(selected),
        wall_time_seconds=0.0,
        peak_rss_bytes=0,
        determinant_evaluations=1,
        sigma_evaluations=1,
        search_left_Omega=0.0,
        search_right_Omega=float(selected[-1]["Omega"]) + 5.0,
        local_refinements=1,
        solve_mode="SYNTHETIC_LOCAL",
        fallback_used=False,
        unresolved_candidates_below_root9=0,
        continuation_leg="ARM_SWAP_DIAGNOSTIC" if swapped else "LOCAL_REPAIR",
        swapped_arms=swapped,
    )


def _write_spectrum(path: Path, rows: list[dict[str, Any]]) -> None:
    sweep.rlb2e._atomic_write_csv(path, rows, sweep.SPECTRUM_FIELDS)


def test_01_exact_base_material_contract() -> None:
    assert sweep.base_material_contract() == {
        "delta": 0.1,
        "E1": 1.1,
        "E2": 0.9,
        "nu12": 0.3,
        "G12": 1.0 / 2.6,
        "G13": 1.0 / 2.6,
        "G23": 1.0 / 2.6,
        "rho": 1.0,
    }


@pytest.mark.parametrize("xi", [-0.8, -0.4, 0.0, 0.4, 0.8])
def test_02_signed_material_scaling_and_fixed_nu_density(xi: float) -> None:
    materials = sweep.signed_materials(xi)
    factors = {
        sweep.OUTER_PLY_MATERIAL: 1.0 + xi,
        sweep.INNER_PLY_MATERIAL: 1.0 - xi,
        sweep.HOMOGENEOUS_M0: 1.0,
    }
    for label, factor in factors.items():
        material = materials[label]
        assert material.E1 == pytest.approx(factor * 1.1)
        assert material.E2 == pytest.approx(factor * 0.9)
        assert material.G12 == pytest.approx(factor / 2.6)
        assert material.G13 == pytest.approx(factor / 2.6)
        assert material.G23 == pytest.approx(factor / 2.6)
        assert material.nu12 == 0.3
        assert material.rho == 1.0
    assert sweep.material_scales(xi) == pytest.approx(factors)


def test_03_all_moduli_are_positive_on_the_complete_signed_grid() -> None:
    assert len(sweep.xi_grid()) == 81
    for xi in sweep.xi_grid():
        for material in sweep.signed_materials(float(xi)).values():
            assert min(
                material.E1,
                material.E2,
                material.G12,
                material.G13,
                material.G23,
            ) > 0.0


def test_04_exact_four_equal_zero_degree_plies_and_stack_labels() -> None:
    layered, homogeneous = sweep.build_arm_sections(0.4)
    assert layered.layout == sweep.LAYERED_LAYOUT == (
        sweep.OUTER_PLY_MATERIAL,
        sweep.INNER_PLY_MATERIAL,
        sweep.INNER_PLY_MATERIAL,
        sweep.OUTER_PLY_MATERIAL,
    )
    assert homogeneous.layout == sweep.HOMOGENEOUS_LAYOUT == (
        sweep.HOMOGENEOUS_M0,
    ) * 4
    for section in (layered, homogeneous):
        assert len(section.laminate.plies) == 4
        assert [ply.angle_deg for ply in section.laminate.plies] == [0.0] * 4
        assert_allclose(
            [ply.thickness for ply in section.laminate.plies],
            [0.0125] * 4,
            rtol=0.0,
            atol=0.0,
        )
        assert_allclose(
            section.laminate.z_interfaces,
            [-0.025, -0.0125, 0.0, 0.0125, 0.025],
            rtol=0.0,
            atol=64.0 * np.finfo(float).eps * sweep.THICKNESS,
        )
    swapped = sweep.build_arm_sections(0.4, swapped=True)
    assert swapped[0].layout == sweep.HOMOGENEOUS_LAYOUT
    assert swapped[1].layout == sweep.LAYERED_LAYOUT


@pytest.mark.parametrize("xi", [-0.8, -0.4, 0.0, 0.4, 0.8])
def test_05_constitutive_invariants_and_analytic_D_formula(xi: float) -> None:
    layered = sweep.build_layered_section(xi)
    homogeneous = sweep.build_homogeneous_section()
    expected = 1.0 + 0.75 * xi
    assert_allclose(
        layered.laminate.D,
        expected * homogeneous.laminate.D,
        rtol=sweep.MATRIX_RELATIVE_TOLERANCE,
        atol=1.0e-16,
    )
    assert_allclose(
        layered.laminate.A,
        homogeneous.laminate.A,
        rtol=sweep.MATRIX_RELATIVE_TOLERANCE,
        atol=1.0e-16,
    )
    assert_allclose(
        layered.laminate.shear,
        homogeneous.laminate.shear,
        rtol=sweep.MATRIX_RELATIVE_TOLERANCE,
        atol=1.0e-16,
    )
    assert layered.laminate.I0 == pytest.approx(homogeneous.laminate.I0, rel=1.0e-12)
    assert layered.laminate.I2 == pytest.approx(homogeneous.laminate.I2, rel=1.0e-12)
    for name in ("A", "S", "m", "J"):
        assert getattr(layered.properties, name) == pytest.approx(
            getattr(homogeneous.properties, name),
            rel=sweep.REDUCED_PROPERTY_TOLERANCE,
        )
    assert sweep.rlb2e._scaled_B(layered.laminate) <= sweep.SYMMETRY_RELATIVE_TOLERANCE
    assert sweep.rlb2e._scaled_I1(layered.laminate) <= sweep.SYMMETRY_RELATIVE_TOLERANCE
    assert layered.properties.D / homogeneous.properties.D == pytest.approx(
        expected, rel=sweep.REDUCED_PROPERTY_TOLERANCE
    )
    assert sweep.rlb2e._reduction_max_residual(layered.properties) <= (
        sweep.REDUCTION_ROUTE_TOLERANCE
    )


def test_06_endpoint_D_ratios_and_total_stiffness_are_exact() -> None:
    homogeneous = sweep.build_homogeneous_section()
    for xi, expected in ((-0.8, 0.4), (0.0, 1.0), (0.8, 1.6)):
        layered = sweep.build_layered_section(xi)
        actual = layered.properties.D / homogeneous.properties.D
        assert actual == pytest.approx(expected, rel=1.0e-11)
        assert 1.0 + actual == pytest.approx(2.0 + 0.75 * xi, rel=1.0e-11)


def test_07_constitutive_gate_passes_all_declared_values() -> None:
    gate = sweep.constitutive_gate()
    assert gate["status"] == "PASS"
    assert gate["full_grid_moduli_positive"] is True
    assert [item["xi"] for item in gate["checks"]] == [-0.8, -0.4, 0.0, 0.4, 0.8]
    assert all(item["status"] == "PASS" for item in gate["checks"])
    assert gate["maximum_residuals"]["D_matrix_formula_relative"] <= 1.0e-12
    assert gate["maximum_residuals"]["Dbeam_formula_relative"] <= 1.0e-11
    assert gate["maximum_residuals"]["symmetry_relative"] <= 1.0e-12
    assert gate["maximum_residuals"]["reduction_route_relative"] <= 1.0e-11


def test_08_section_property_table_has_162_unambiguous_arm_rows() -> None:
    rows = sweep.section_property_rows()
    assert len(rows) == 162
    assert len({round(float(row["xi"]), 10) for row in rows}) == 81
    for xi in sweep.xi_grid():
        group = [row for row in rows if float(row["xi"]) == pytest.approx(float(xi))]
        assert len(group) == 2
        assert {int(row["arm_id"]) for row in group} == {1, 2}
        assert {str(row["arm_role"]) for row in group} == {
            "LAYERED",
            "HOMOGENEOUS_REFERENCE",
        }
        assert all(row["ply_index"] == [1, 2, 3, 4] for row in group)
        assert all(row["ply_region"] == ["OUTER", "INNER", "INNER", "OUTER"] for row in group)
        assert all(row["angle_deg"] == [0.0] * 4 for row in group)
        assert all(row["ply_thickness"] == [0.0125] * 4 for row in group)
        assert all(row["constitutive_status"] == "PASS" for row in group)


def test_09_xi_zero_is_the_fully_homogeneous_characteristic_matrix() -> None:
    layered = sweep.build_layered_section(0.0)
    homogeneous = sweep.build_homogeneous_section()
    assert_allclose(layered.laminate.A, homogeneous.laminate.A, rtol=1.0e-13)
    assert_allclose(layered.laminate.D, homogeneous.laminate.D, rtol=1.0e-13)
    assert_allclose(layered.laminate.shear, homogeneous.laminate.shear, rtol=1.0e-13)
    for name in ("A", "D", "S", "m", "J"):
        assert getattr(layered.properties, name) == pytest.approx(
            getattr(homogeneous.properties, name), rel=1.0e-13
        )

    provider, metadata = sweep.make_matrix_provider(0.0)
    _beam, coupled = sweep.rlb2e._physics_modules()
    for omega in (0.731, 3.217):
        direct = coupled.coupled_boundary_matrix(
            omega,
            math.radians(sweep.BETA_DEG),
            sweep.L1,
            homogeneous.properties,
            sweep.L2,
            homogeneous.properties,
        )
        assert_allclose(provider(omega), direct, rtol=0.0, atol=16.0 * np.finfo(float).eps)
    assert metadata["production_modules_only"] is True


def test_10_arm_swap_diagnostic_uses_only_four_declared_points(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = _synthetic_rows()
    calls: list[tuple[float, bool]] = []

    def fake_solve(xi: float, **kwargs: Any) -> sweep.PointSolution:
        calls.append((float(xi), bool(kwargs.get("swapped", False))))
        return _synthetic_solution(rows, xi, swapped=True)

    monkeypatch.setattr(sweep, "solve_point", fake_solve)
    result = sweep.arm_swap_checks(rows)
    assert result["status"] == "PASS"
    assert [item["xi"] for item in result["checks"]] == [-0.8, -0.4, 0.4, 0.8]
    assert calls == [(-0.8, True), (-0.4, True), (0.4, True), (0.8, True)]
    assert all(item["root_count"] == 9 for item in result["checks"])
    assert all(not item["roots_above_9_computed"] for item in result["checks"])


def test_11_exact_81_point_grid_and_two_continuation_legs() -> None:
    grid = sweep.xi_grid()
    negative, positive = sweep.continuation_paths()
    assert_allclose(grid, np.linspace(-0.8, 0.8, 81), rtol=0.0, atol=1.0e-14)
    assert_allclose(negative, np.linspace(0.0, -0.8, 41), rtol=0.0, atol=1.0e-14)
    assert_allclose(positive, np.linspace(0.0, 0.8, 41), rtol=0.0, atol=1.0e-14)
    assert set(np.round(negative, 10)) & set(np.round(positive, 10)) == {0.0}
    assert len(set(np.round(np.concatenate((negative, positive)), 10))) == 81


def test_12_frozen_search_thresholds_are_reused_and_predictor_is_locator_only() -> None:
    policy = sweep.rlb2e._rlb2e_search_policy()
    predecessor = sweep.rlb2e._root_tools().SearchPolicy()
    assert policy.requested_roots == 8
    assert policy.guard_roots == 1
    assert policy.required_slots == 9
    for name in (
        "sigma_ratio_tolerance",
        "rank_relative_tolerance",
        "boundary_residual_tolerance",
        "root_xtol_bar",
        "root_rtol",
        "dedup_atol_bar",
        "dedup_rtol",
        "cluster_atol_bar",
        "cluster_rtol",
    ):
        assert getattr(policy, name) == getattr(predecessor, name)

    slots: list[Any] = []
    for position in range(1, sweep.K_GUARD + 1):
        Omega = 10.0 * position
        diagnostics = SimpleNamespace(
            raw_determinant=0.0,
            scaled_determinant=0.0,
            raw_sigma_ratio=0.0,
            scaled_sigma_ratio=0.0,
            raw_boundary_null_residual=0.0,
            detected_nullity=1,
        )
        candidate = SimpleNamespace(
            interval_left_bar=Omega - 0.1,
            interval_right_bar=Omega + 0.1,
            detection_sources=("determinant", "sigma"),
            diagnostics=diagnostics,
        )
        slots.append(
            SimpleNamespace(
                event=SimpleNamespace(
                    omega_bar=Omega,
                    omega=Omega / sweep.OMEGA_TO_OMEGA_SCALE,
                    candidate=candidate,
                )
            )
        )
    predicted = np.arange(1.0, 10.0) + 0.25
    rows = sweep._root_rows(
        0.2,
        slots,
        solve_id="synthetic",
        transaction_id="synthetic",
        solve_mode="FAST_LOCAL",
        fallback_used=False,
        predicted=predicted,
        search_right=95.0,
        unresolved=0,
        continuation_leg="POSITIVE",
    )
    assert [float(row["Omega"]) for row in rows] == [10.0 * value for value in range(1, 10)]
    assert [float(row["predictor_Omega"]) for row in rows] == pytest.approx(predicted)
    assert all(row["predictor_used_as_final"] is False for row in rows)


def test_13_base_inventory_audit_requires_729_ordered_rows() -> None:
    rows = _synthetic_rows()
    audit = sweep.audit_spectrum_rows(rows)
    assert audit["status"] == "PASS"
    assert audit["base_group_count"] == 81
    assert audit["base_row_count"] == 729
    assert audit["duplicate_row_id_count"] == 0
    assert audit["canonical_source_failures"] == []
    for group in sweep._base_group_index(rows).values():
        ordered = sorted(group, key=lambda row: int(row["sorted_position"]))
        assert [int(row["sorted_position"]) for row in ordered] == list(range(1, 10))
        assert np.all(np.diff([float(row["Omega"]) for row in ordered]) > 0.0)
    assert sweep.audit_spectrum_rows(rows[:-1])["status"] == "FAIL"


def test_14_root9_is_guard_excluded_from_plot_and_no_tail_exists() -> None:
    rows = _synthetic_rows()
    guards = [row for row in rows if int(row["sorted_position"]) == 9]
    assert len(guards) == 81
    assert all(row["root_role"] == "ROOT_9_GUARD" for row in guards)
    assert all(row["guard_flag"] for row in guards)
    assert all(float(row["root9_right_margin_Omega"]) >= 2.0 for row in guards)
    plotted = sweep.canonical_plot_rows(rows)
    assert len(plotted) == 81 * 8
    assert max(int(row["sorted_position"]) for row in plotted) == 8
    assert max(int(row["sorted_position"]) for row in rows) == 9
    assert all(not row["roots_above_9_computed"] for row in rows)
    assert all(not row["predictor_used_as_final"] for row in rows)


def test_15_atomic_completed_group_is_preserved_byte_for_byte(tmp_path: Path) -> None:
    rows = _synthetic_rows()[: sweep.K_GUARD]
    path = tmp_path / sweep.SPECTRUM_FILENAME
    _write_spectrum(path, rows)
    before = path.read_bytes()
    solution = _synthetic_solution(_synthetic_rows(), float(rows[0]["xi"]))
    preserved = sweep._write_point_transaction(tmp_path, rows, solution)
    assert preserved == rows
    assert path.read_bytes() == before


def test_16_neighbour_spike_audit_flags_without_smoothing() -> None:
    rows = _synthetic_rows()
    target = next(
        row
        for row in rows
        if float(row["xi"]) == pytest.approx(0.4)
        and int(row["sorted_position"]) == 4
    )
    target["Omega"] = 45.0
    target["omega"] = 45.0 / sweep.OMEGA_TO_OMEGA_SCALE
    target["Lambda"] = math.sqrt(45.0)
    audit = sweep.neighbour_audit_rows(rows)
    assert any(float(row["xi"]) == pytest.approx(0.4) and row["flagged"] for row in audit)
    assert all(row["smoothing_applied"] is False for row in audit)


def test_17_only_flagged_parameter_points_are_selected_for_repair() -> None:
    audit = [
        {"xi": -0.2, "flagged": False},
        {"xi": 0.4, "flagged": True},
        {"xi": 0.4, "flagged": True},
    ]
    assert sweep.flagged_repair_points(audit) == [0.4]


def test_18_reproduced_feature_is_retained_not_smoothed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = _synthetic_rows()
    original = sweep._rows_for_roots(rows, 0.4).copy()
    audit = [{"xi": 0.4, "sorted_position": 4, "flagged": True}]
    calls: list[float] = []

    def fake_solve(xi: float, **_kwargs: Any) -> sweep.PointSolution:
        calls.append(float(xi))
        return _synthetic_solution(rows, xi)

    monkeypatch.setattr(sweep, "solve_point", fake_solve)
    repaired, updated_audit, records = sweep.apply_local_repairs(rows, audit)
    assert calls == [0.4]
    assert records[0]["status"] == "REPRODUCED_AFTER_LOCAL_REPAIR"
    assert records[0]["smoothing_applied"] is False
    assert updated_audit[0]["repair_status"] == "REPRODUCED_AFTER_LOCAL_REPAIR"
    assert_allclose(sweep._rows_for_roots(repaired, 0.4), original, rtol=0.0, atol=0.0)
    assert sweep.audit_spectrum_rows(repaired)["status"] == "PASS"


def test_19_unresolved_local_repair_creates_a_visible_nan_gap(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = _synthetic_rows()
    audit = [{"xi": -0.4, "sorted_position": 3, "flagged": True}]

    def fail_solve(*_args: Any, **_kwargs: Any) -> sweep.PointSolution:
        raise RuntimeError("synthetic unresolved neighbourhood")

    monkeypatch.setattr(sweep, "solve_point", fail_solve)
    repaired, updated_audit, records = sweep.apply_local_repairs(rows, audit)
    canonical = sweep._canonical_group(repaired, -0.4)
    assert records[0]["status"] == "UNRESOLVED"
    assert updated_audit[0]["repair_status"] == "UNRESOLVED_AFTER_LOCAL_REPAIR"
    assert math.isnan(float(canonical[2]["Lambda"]))
    assert canonical[2]["point_status"] == "UNRESOLVED_AFTER_LOCAL_REPAIR"
    assert sweep.audit_plot_rows(repaired)["status"] == "PASS"


def test_20_plot_only_has_one_panel_eight_lines_and_no_physics_calls(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    _write_spectrum(tmp_path / sweep.SPECTRUM_FILENAME, _synthetic_rows())

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("plot_only called production physics or root search")

    monkeypatch.setattr(sweep, "solve_point", forbidden)
    monkeypatch.setattr(sweep, "make_matrix_provider", forbidden)
    monkeypatch.setattr(sweep, "build_layered_section", forbidden)
    monkeypatch.setattr(sweep.rlb2e, "_physics_modules", forbidden)
    result = sweep.create_plot_from_csv(tmp_path)
    assert result["root_calculation_count"] == 0
    assert result["panel_count"] == 1
    assert result["spectral_line_count"] == 8
    assert result["reference_line_count"] == 1
    assert result["plotted_positions"] == list(range(1, 9))
    assert result["root9_plotted"] is False
    assert (tmp_path / sweep.PLOT_FILENAME).is_file()


def test_21_completed_missing_only_run_is_idempotent(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    spectrum = tmp_path / sweep.SPECTRUM_FILENAME
    _write_spectrum(spectrum, _synthetic_rows())
    for name in sweep.MANDATORY_COMPLETED_OUTPUTS - {sweep.SPECTRUM_FILENAME}:
        (tmp_path / name).write_bytes(f"synthetic {name}\n".encode("utf-8"))
    output_hashes = {
        name: sweep.rlb2e._sha256(tmp_path / name)
        for name in sweep.MANDATORY_COMPLETED_OUTPUTS
    }
    manifest = {
        "stage_id": sweep.STAGE_ID,
        "algorithm_version": sweep.ALGORITHM_VERSION,
        "contract_sha256": sweep.contract_hash(),
        "analysis_script_sha256": sweep.rlb2e._sha256(SCRIPT_PATH),
        "production_physics_hashes": {
            path: sweep.rlb2e._sha256(ROOT / path)
            for path in sweep.PRODUCTION_PHYSICS_PATHS
        },
        "reused_RLB2E_script_sha256": sweep.rlb2e._sha256(
            ROOT / sweep.RLB2E_SCRIPT_PATH
        ),
        "output_hashes": output_hashes,
    }
    sweep.rlb2e._atomic_write_json(tmp_path / sweep.MANIFEST_FILENAME, manifest)
    before = {path.name: path.read_bytes() for path in tmp_path.iterdir()}

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("completed --missing-only attempted a root calculation")

    monkeypatch.setattr(sweep, "solve_point", forbidden)
    monkeypatch.setattr(sweep, "constitutive_gate", forbidden)
    result = sweep.run_workflow(tmp_path, missing_only=True)
    after = {path.name: path.read_bytes() for path in tmp_path.iterdir()}
    assert after == before
    assert result["invocation"] == {
        "missing_only": True,
        "root_calculation_count": 0,
        "outputs_modified": False,
    }


def test_22_manifest_contains_complete_frequency_map_policy_instance() -> None:
    contract = sweep.contract_payload()
    assert contract["frequency_map_policy"] == {
        "frequency_map_policy": "frequency-map-v1",
        "calculation_mode": "fast_plot",
        "spectrum_semantics": "sorted_positions",
        "sweep_parameter": "xi",
        "parameter_grid": "-0.80:0.02:0.80",
        "continuation_anchor": "0.00",
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
    assert contract["root_search_reuse"] == {
        "source": sweep.RLB2E_SCRIPT_PATH,
        "requested_roots": 8,
        "guard_roots": 1,
        "required_slots": 9,
        "physics_formulas_copied": False,
    }
    assert contract["continuation"]["xi0_calculated_once"] is True
    assert len(contract["xi_grid"]) == 81


def test_23_normalization_and_benchmark_contract_are_unchanged() -> None:
    omega = 2.75
    Omega = sweep.omega_to_Omega(omega)
    assert sweep.OMEGA_TO_OMEGA_SCALE == pytest.approx(math.sqrt(12.0 / 0.05**2))
    assert Omega == pytest.approx(omega * sweep.OMEGA_TO_OMEGA_SCALE)
    assert sweep.Omega_to_Lambda(Omega) ** 2 == pytest.approx(Omega)
    benchmark = sweep._benchmark_payload(
        [
            {"xi": 0.0, "wall_time_seconds": 1.0},
            {"xi": -0.8, "wall_time_seconds": 2.0},
            {"xi": 0.8, "wall_time_seconds": 1.5},
        ]
    )
    assert benchmark["total_unique_base_points"] == 81
    assert benchmark["declared_anchor_points"] == 3
    assert benchmark["remaining_base_points"] == 78
    assert benchmark["conservative_eta_seconds"] == pytest.approx(234.0)


def test_24_no_forbidden_imports_local_physics_solver_or_smoothing() -> None:
    source = SCRIPT_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    imported: list[str] = []
    functions: set[str] = set()
    classes: set[str] = set()
    called_names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            functions.add(node.name)
        elif isinstance(node, ast.ClassDef):
            classes.add(node.name)
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                called_names.add(node.func.id)
            elif isinstance(node.func, ast.Attribute):
                called_names.add(node.func.attr)
    import_text = " ".join(imported).lower()
    for forbidden in (
        "spectral_sweep_runner",
        "ritz",
        "fem",
        "yartsev",
        "shellbuckling",
        "isotropic_rectangular_timoshenko",
        "circular",
    ):
        assert forbidden not in import_text
    assert classes == {"PointSolution"}
    assert not any("characteristic_determinant" in name for name in functions)
    assert "determinant" not in functions
    assert not ({"interp", "interp1d", "savgol_filter", "UnivariateSpline"} & called_names)
    assert "spectral_sweep_runner_used\": False" in source
    assert "smoothing\": False" in source


def test_guard_detector_reconciliation_requires_one_overlapping_neighbourhood() -> None:
    def slot(value: float, left: float, right: float, nullity: int = 1) -> Any:
        diagnostics = SimpleNamespace(detected_nullity=nullity)
        candidate = SimpleNamespace(
            interval_left_bar=left,
            interval_right_bar=right,
            diagnostics=diagnostics,
        )
        return SimpleNamespace(
            event=SimpleNamespace(omega_bar=value, candidate=candidate)
        )

    slots = [slot(float(index), index - 0.1, index + 0.1) for index in range(1, 9)]
    slots.extend(
        (
            slot(136.5667868954, 136.53, 136.61),
            slot(136.5667868969, 136.55, 136.59),
        )
    )
    policy = SimpleNamespace(
        reference_detector_reconciliation_atol_bar=2.0e-10,
        reference_detector_reconciliation_rtol=5.0e-8,
    )
    assert sweep._guard_detector_duplicate(slots, policy)

    slots[-1] = slot(136.8, 136.75, 136.85)
    assert not sweep._guard_detector_duplicate(slots, policy)


def test_25_predecessor_output_hashes_are_preserved_if_available() -> None:
    manifest_path = sweep.RLB2E_RESULT_DIR / sweep.rlb2e.MANIFEST_FILENAME
    if not manifest_path.is_file():
        pytest.skip("Frozen RLB-2E output tree is not present.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest.get("output_hashes")
    for name, expected in manifest["output_hashes"].items():
        path = sweep.RLB2E_RESULT_DIR / name
        assert path.is_file()
        assert sweep.rlb2e._sha256(path) == expected


def test_26_generated_outputs_if_present() -> None:
    if not (RESULT_DIR / sweep.MANIFEST_FILENAME).is_file():
        pytest.skip("RLB-2F production outputs have not been finalized yet.")
    required = {
        sweep.SPECTRUM_FILENAME,
        sweep.SECTION_FILENAME,
        sweep.AUDIT_FILENAME,
        sweep.BENCHMARK_FILENAME,
        sweep.CHECKPOINT_FILENAME,
        sweep.ARM_SWAP_FILENAME,
        sweep.PLOT_FILENAME,
        sweep.REPORT_FILENAME,
        sweep.MANIFEST_FILENAME,
    }
    assert required <= {path.name for path in RESULT_DIR.iterdir() if path.is_file()}
    manifest_path = RESULT_DIR / sweep.MANIFEST_FILENAME
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    spectrum = sweep.rlb2e._read_csv(RESULT_DIR / sweep.SPECTRUM_FILENAME)
    sections = sweep.rlb2e._read_csv(RESULT_DIR / sweep.SECTION_FILENAME)
    arm_swap = json.loads((RESULT_DIR / sweep.ARM_SWAP_FILENAME).read_text(encoding="utf-8"))

    inventory = sweep.audit_spectrum_rows(spectrum)
    assert inventory["status"] == "PASS"
    assert inventory["base_group_count"] == 81
    assert inventory["base_row_count"] == 729
    assert inventory["roots_above_guard_count"] == 0
    assert len(sections) == 162
    assert len({round(float(row["xi"]), 10) for row in sections}) == 81
    assert {row["arm_role"] for row in sections} == {
        "LAYERED",
        "HOMOGENEOUS_REFERENCE",
    }
    assert manifest["counts"]["base_points"] == 81
    assert manifest["counts"]["base_rows"] == 729
    assert manifest["counts"]["root9_guards"] == 81
    assert manifest["root_contract"]["root9_plotted"] is False
    assert manifest["root_contract"]["roots_above_9_computed"] is False
    assert manifest["root_contract"]["accepted_candidates_above_root9"] == 0
    assert manifest["baseline_consistency"]["status"] == "PASS"
    assert arm_swap["status"] == "PASS"
    assert [item["xi"] for item in arm_swap["checks"]] == [-0.8, -0.4, 0.4, 0.8]
    assert all(item["root_count"] == 9 for item in arm_swap["checks"])
    assert all(item["status"] == "PASS" for item in arm_swap["checks"])
    assert manifest["contract"]["frequency_map_policy"] == sweep.FREQUENCY_MAP_POLICY
    assert manifest["analysis_script_sha256"] == sweep.rlb2e._sha256(SCRIPT_PATH)
    for path, expected in manifest["production_physics_hashes"].items():
        assert sweep.rlb2e._sha256(ROOT / path) == expected
    for name, expected in manifest["output_hashes"].items():
        assert sweep.rlb2e._sha256(RESULT_DIR / name) == expected
    plot_path = RESULT_DIR / sweep.PLOT_FILENAME
    assert plot_path.stat().st_size > 0
