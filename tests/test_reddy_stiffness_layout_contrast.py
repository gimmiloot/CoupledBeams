"""Targeted gates for the RLB-2E stiffness-layout frequency map."""

from __future__ import annotations

import ast
import csv
import json
import math
from pathlib import Path
import sys
from types import SimpleNamespace
from typing import Any

import numpy as np
from numpy.testing import assert_allclose
import pytest

from scripts.analysis.laminated_beams import (
    sweep_reddy_stiffness_layout_contrast as sweep,
)


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_stiffness_layout_contrast_sweep"
)


def _synthetic_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for configuration_index, configuration in enumerate(sweep.CONFIGURATIONS):
        for chi in sweep.chi_grid():
            for position in range(1, sweep.K_GUARD + 1):
                Omega = 10.0 * position + configuration_index + 0.2 * float(chi)
                rows.append(
                    {
                        "row_id": f"{configuration}_{chi:.2f}_{position}",
                        "configuration_id": configuration,
                        "chi": float(chi),
                        "grid_kind": "BASE",
                        "solve_id": "synthetic",
                        "transaction_id": "synthetic",
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
                        "search_right_Omega": 95.0,
                        "root9_right_margin_Omega": 5.0,
                        "solve_mode": "SYNTHETIC",
                        "fallback_used": False,
                        "quality_status": "PASS",
                        "point_status": "PASS",
                        "shared_chi0_anchor_reused": chi == 0.0,
                        "shared_chi0_source_configuration": (
                            sweep.CONFIG_BOTH_OUTER if chi == 0.0 else ""
                        ),
                        "is_canonical_plot_source": True,
                        "supersedes_row_id": "",
                        "repair_id": "",
                        "roots_above_9_computed": False,
                    }
                )
    return rows


def test_exact_base_material_contract() -> None:
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


@pytest.mark.parametrize("chi", [0.0, 0.4, 0.8])
def test_H_L_material_scaling_and_fixed_density_nu(chi: float) -> None:
    materials = sweep.contrast_materials(chi)
    for label, factor in (("H", 1.0 + chi), ("L", 1.0 - chi)):
        material = materials[label]
        assert material.E1 == pytest.approx(factor * 1.1)
        assert material.E2 == pytest.approx(factor * 0.9)
        assert material.G12 == pytest.approx(factor / 2.6)
        assert material.G13 == pytest.approx(factor / 2.6)
        assert material.G23 == pytest.approx(factor / 2.6)
        assert material.nu12 == 0.3
        assert material.rho == 1.0


def test_all_moduli_are_positive_on_exact_41_point_grid() -> None:
    assert_allclose(sweep.chi_grid(), np.linspace(0.0, 0.8, 41), atol=0.0)
    for chi in sweep.chi_grid():
        for material in sweep.contrast_materials(float(chi)).values():
            assert min(
                material.E1,
                material.E2,
                material.G12,
                material.G13,
                material.G23,
            ) > 0.0


@pytest.mark.parametrize("layout", [sweep.OUTER_LAYOUT, sweep.INNER_LAYOUT])
def test_exact_four_equal_zero_degree_plies(layout: tuple[str, ...]) -> None:
    section = sweep.build_layout_section(layout, 0.4)
    assert len(section.laminate.plies) == 4
    assert tuple(ply.label for ply in section.laminate.plies) == layout
    assert [ply.angle_deg for ply in section.laminate.plies] == [0.0] * 4
    assert_allclose(
        [ply.thickness for ply in section.laminate.plies],
        np.full(4, 0.0125),
        rtol=0.0,
        atol=0.0,
    )
    assert_allclose(
        section.laminate.z_interfaces,
        [-0.025, -0.0125, 0.0, 0.0125, 0.025],
        rtol=0.0,
        atol=64.0 * np.finfo(float).eps * sweep.THICKNESS,
    )


def test_three_configuration_stack_contract() -> None:
    assert sweep.CONFIGURATION_LAYOUTS == {
        "BOTH_OUTER_STIFF": (sweep.OUTER_LAYOUT, sweep.OUTER_LAYOUT),
        "BOTH_INNER_STIFF": (sweep.INNER_LAYOUT, sweep.INNER_LAYOUT),
        "ANTI_PHASE": (sweep.INNER_LAYOUT, sweep.OUTER_LAYOUT),
    }


def test_constitutive_gate_passes_frozen_checks() -> None:
    gate = sweep.constitutive_gate()
    assert gate["status"] == "PASS"
    assert all(item["status"] == "PASS" for item in gate["checks"])
    assert gate["maximum_residuals"]["D_matrix_formula_relative"] <= 1.0e-12
    assert gate["maximum_residuals"]["Dbeam_formula_relative"] <= 1.0e-11
    assert gate["maximum_residuals"]["symmetry_relative"] <= 1.0e-12
    assert gate["maximum_residuals"]["reduction_route_relative"] <= 1.0e-11


@pytest.mark.parametrize("chi", [0.0, 0.4, 0.8])
def test_A_shear_mass_invariant_and_D_formulas(chi: float) -> None:
    baseline = sweep._baseline_section()
    outer = sweep.build_layout_section(sweep.OUTER_LAYOUT, chi)
    inner = sweep.build_layout_section(sweep.INNER_LAYOUT, chi)
    assert_allclose(outer.laminate.A, inner.laminate.A, rtol=1.0e-12)
    assert_allclose(outer.laminate.shear, inner.laminate.shear, rtol=1.0e-12)
    for name in ("A", "S", "m", "J"):
        assert getattr(outer.properties, name) == pytest.approx(
            getattr(inner.properties, name), rel=1.0e-11
        )
    assert outer.properties.D / baseline.properties.D == pytest.approx(
        1.0 + 0.75 * chi, rel=1.0e-11
    )
    assert inner.properties.D / baseline.properties.D == pytest.approx(
        1.0 - 0.75 * chi, rel=1.0e-11
    )
    assert outer.properties.D + inner.properties.D == pytest.approx(
        2.0 * baseline.properties.D, rel=1.0e-11
    )


def test_section_property_rows_are_complete_and_unambiguous() -> None:
    rows = sweep.section_property_rows()
    assert len(rows) == 3 * 41 * 2
    assert all(row["grid_kind"] == "BASE" for row in rows)
    assert all(row["Dbeam0"] > 0.0 for row in rows)
    assert all(row["analytic_D_beam_ratio_residual"] <= 1.0e-11 for row in rows)


def test_old_frequency_normalization_is_exact() -> None:
    assert sweep.OMEGA_TO_OMEGA_SCALE == pytest.approx(math.sqrt(12.0 / 0.05**2))
    omega = 2.75
    Omega = sweep.omega_to_Omega(omega)
    assert sweep.Omega_to_Lambda(Omega) ** 2 == pytest.approx(Omega)


def test_hold_and_secant_predictor_is_locator_only() -> None:
    old = np.arange(1.0, 10.0)
    current = 2.0 * old
    hold = sweep.hold_secant_predictor(0.02, 0.0, old)
    secant = sweep.hold_secant_predictor(0.04, 0.02, current, 0.0, old)
    assert_allclose(hold, old)
    assert_allclose(secant, 3.0 * old)
    assert "predictor_used_as_final" in sweep.SPECTRUM_FIELDS


def test_local_windows_are_ordered_nonoverlapping_and_retain_close_roots() -> None:
    roots = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 5.0001, 7.0, 8.0, 9.0])
    windows = sweep.local_search_windows(roots)
    assert len(windows) == 9
    assert all(left < right for left, right in windows)
    assert all(windows[index][1] == windows[index + 1][0] for index in range(8))
    assert windows[4] != windows[5]
    assert windows[-1][1] == pytest.approx(roots[-1] + 1.2)


def test_search_policy_retains_exactly_first_eight_plus_root9() -> None:
    policy = sweep._rlb2e_search_policy()
    assert policy.requested_roots == 8
    assert policy.guard_roots == 1
    assert policy.required_slots == 9
    # Numerical quality thresholds remain the frozen predecessor values.
    predecessor = sweep._root_tools().SearchPolicy()
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


def test_benchmark_eta_keeps_the_initial_118_point_campaign() -> None:
    records = [
        {"configuration_id": sweep.CONFIG_BOTH_OUTER, "chi": 0.0, "wall_time_seconds": 1.0},
        {"configuration_id": sweep.CONFIG_BOTH_OUTER, "chi": 0.8, "wall_time_seconds": 2.0},
        {"configuration_id": sweep.CONFIG_ANTI_PHASE, "chi": 0.8, "wall_time_seconds": 1.5},
    ]
    benchmark = sweep._benchmark_payload(records)
    assert benchmark["total_unique_base_solves"] == 121
    assert benchmark["declared_anchor_solves"] == 3
    assert benchmark["remaining_unique_root_points"] == 118
    assert benchmark["conservative_eta_seconds"] == pytest.approx(354.0)


def test_base_group_audit_requires_positions_one_through_guard() -> None:
    rows = _synthetic_rows()
    audit = sweep.audit_spectrum_rows(rows)
    assert audit["status"] == "PASS"
    assert audit["base_group_count"] == 123
    assert audit["base_row_count"] == 1107
    damaged = rows[:-1]
    assert sweep.audit_spectrum_rows(damaged)["status"] == "FAIL"


def test_complete_but_quality_bad_group_is_replaced_not_reused(
    tmp_path: Path,
) -> None:
    good = [dict(row) for row in _synthetic_rows()[: sweep.K_GUARD]]
    bad = [dict(row) for row in good]
    bad[0]["scaled_sigma_ratio"] = 1.0
    key = (sweep.CONFIG_BOTH_OUTER, 0.0)
    assert key not in sweep._complete_base_group_index(bad)
    solution = sweep.PointSolution(
        configuration_id=sweep.CONFIG_BOTH_OUTER,
        chi=0.0,
        rows=tuple(good),
        wall_time_seconds=0.0,
        peak_rss_bytes=0,
        determinant_evaluations=1,
        sigma_evaluations=1,
        search_left_Omega=0.0,
        search_right_Omega=95.0,
        local_refinements=0,
        solve_mode="SYNTHETIC_REPAIR",
        fallback_used=False,
        unresolved_candidates_below_root9=0,
    )
    replaced = sweep._write_point_transaction(tmp_path, bad, solution)
    assert sweep._base_group_is_acceptable(replaced)
    assert all(float(row["scaled_sigma_ratio"]) == 0.0 for row in replaced)


def test_shared_chi0_rows_have_one_common_physical_provenance() -> None:
    rows = _synthetic_rows()
    chi0 = [row for row in rows if row["chi"] == 0.0]
    assert len(chi0) == 27
    assert {row["solve_id"] for row in chi0} == {"synthetic"}
    assert all(row["shared_chi0_anchor_reused"] for row in chi0)


def test_root9_is_guard_and_excluded_from_plot_data() -> None:
    rows = _synthetic_rows()
    guards = [row for row in rows if row["sorted_position"] == 9]
    assert len(guards) == 123
    assert all(row["root_role"] == "ROOT_9_GUARD" for row in guards)
    plotted = sweep.canonical_plot_rows(rows)
    assert len(plotted) == 3 * 41 * 8
    assert max(int(row["sorted_position"]) for row in plotted) == 8


def test_no_output_contract_allows_roots_above_nine() -> None:
    rows = _synthetic_rows()
    assert max(int(row["sorted_position"]) for row in rows) == 9
    assert all(not row["roots_above_9_computed"] for row in rows)


def test_internal_candidate_above_root9_rejects_point(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    diagnostics = object()
    slots = [
        SimpleNamespace(
            event=SimpleNamespace(
                omega_bar=float(position),
                candidate=SimpleNamespace(diagnostics=diagnostics),
            )
        )
        for position in range(1, sweep.K_GUARD + 1)
    ]
    canonical = [SimpleNamespace(omega_bar=float(position)) for position in range(1, 11)]
    policy = SimpleNamespace(dedup_atol_bar=1.0e-12, dedup_rtol=1.0e-12)
    monkeypatch.setattr(
        sweep,
        "_root_tools",
        lambda: SimpleNamespace(_candidate_quality=lambda _diagnostics, _policy: (True, "PASS")),
    )
    passed, evidence = sweep._point_is_acceptable(
        canonical,
        [],
        slots,
        11.0,
        policy,
    )
    assert passed is False
    assert evidence["accepted_candidates_above_root9"] == 1
    assert evidence["roots_above_9_computed"] is True


def test_centered_secant_and_mad_flag_are_not_smoothing() -> None:
    assert sweep.centred_secant_residual(1.0, 2.0, 1.0) == 0.5
    rows = _synthetic_rows()
    target = next(
        row
        for row in rows
        if row["configuration_id"] == sweep.CONFIG_BOTH_OUTER
        and row["chi"] == pytest.approx(0.4)
        and row["sorted_position"] == 1
    )
    target["Lambda"] = float(target["Lambda"]) * 1.25
    audit = sweep.neighbour_audit_rows(rows)
    flagged = [row for row in audit if row["flagged"]]
    assert flagged
    assert all(not row["smoothing_applied"] for row in audit)


def test_only_flagged_points_are_selected_for_local_repair() -> None:
    audit = [
        {"configuration_id": sweep.CONFIG_BOTH_OUTER, "chi": 0.2, "flagged": False},
        {"configuration_id": sweep.CONFIG_ANTI_PHASE, "chi": 0.4, "flagged": True},
        {"configuration_id": sweep.CONFIG_ANTI_PHASE, "chi": 0.4, "flagged": True},
    ]
    assert sweep.flagged_repair_points(audit) == [(sweep.CONFIG_ANTI_PHASE, 0.4)]


def test_dense_repair_spacing_is_applied_only_to_affected_positions(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    spacings: list[float] = []

    def fake_scan(*_args: Any, spacing: float, **_kwargs: Any) -> list[Any]:
        spacings.append(spacing)
        return []

    monkeypatch.setattr(sweep, "_scan_candidate_window", fake_scan)
    monkeypatch.setattr(sweep, "_canonical_slots", lambda _raw, _policy: ([], [], []))
    counted = SimpleNamespace()
    sweep._local_candidates(
        counted,
        SimpleNamespace(),
        np.arange(1.0, 10.0),
        solve_id="repair",
        dense=True,
        dense_positions=[3, 7],
    )
    assert spacings == [0.1, 0.1, 0.05, 0.1, 0.1, 0.1, 0.05, 0.1, 0.1]


def test_unresolved_plot_value_is_a_nan_gap() -> None:
    rows = _synthetic_rows()
    row = rows[0]
    row["Lambda"] = math.nan
    row["point_status"] = "UNRESOLVED_AFTER_LOCAL_REPAIR"
    selected = sweep.canonical_plot_rows(rows)
    assert math.isnan(float(selected[0]["Lambda"]))
    assert selected[0]["point_status"] == "UNRESOLVED_AFTER_LOCAL_REPAIR"


def test_atomic_csv_write_and_complete_group_preservation(tmp_path: Path) -> None:
    rows = _synthetic_rows()[:9]
    path = tmp_path / "roots.csv"
    sweep._atomic_write_csv(path, rows, sweep.SPECTRUM_FIELDS)
    before = path.read_bytes()
    solution = sweep.PointSolution(
        configuration_id=sweep.CONFIG_BOTH_OUTER,
        chi=0.0,
        rows=tuple(rows),
        wall_time_seconds=0.0,
        peak_rss_bytes=0,
        determinant_evaluations=0,
        sigma_evaluations=0,
        search_left_Omega=0.0,
        search_right_Omega=95.0,
        local_refinements=0,
        solve_mode="REUSE",
        fallback_used=False,
        unresolved_candidates_below_root9=0,
    )
    preserved = sweep._write_point_transaction(tmp_path, rows, solution)
    assert preserved == rows
    assert path.read_bytes() == before


def test_checkpoint_contract_has_resume_and_missing_fields() -> None:
    payload = sweep._checkpoint_payload(
        [],
        [],
        constitutive={"status": "PASS"},
        started_at="2026-08-28T00:00:00+00:00",
        benchmark_status="PENDING",
    )
    assert payload["completed_base_groups"] == 0
    assert len(payload["missing_points"]) == 123
    assert payload["failed_points"] == []
    assert payload["terminal_unresolved_points"] == []
    assert payload["contract_sha256"] == sweep.contract_hash()


def test_manifest_only_does_not_mutate_existing_checkpoint(tmp_path: Path) -> None:
    checkpoint = tmp_path / sweep.CHECKPOINT_FILENAME
    checkpoint.write_bytes(b"immutable checkpoint\n")
    result = sweep.manifest_only(tmp_path)
    assert result["root_calculation_count"] == 0
    assert result["resume_artifacts_modified"] is False
    assert checkpoint.read_bytes() == b"immutable checkpoint\n"


def test_plot_only_uses_csv_and_never_calls_physics_or_roots(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rows = _synthetic_rows()
    sweep._atomic_write_csv(tmp_path / sweep.SPECTRUM_FILENAME, rows, sweep.SPECTRUM_FIELDS)

    def forbidden() -> None:
        raise AssertionError("plot_only imported or called the root path")

    monkeypatch.setattr(sweep, "_physics_modules", forbidden)
    monkeypatch.setattr(sweep, "_root_tools", forbidden)
    result = sweep.create_plot_from_csv(tmp_path)
    assert result["root_calculation_count"] == 0
    assert result["panel_count"] == 3
    assert result["lines_per_panel"] == 8
    assert result["root9_plotted"] is False
    assert (tmp_path / sweep.PLOT_FILENAME).is_file()


def test_plot_only_rejects_missing_canonical_value(tmp_path: Path) -> None:
    rows = _synthetic_rows()[1:]
    sweep._atomic_write_csv(tmp_path / sweep.SPECTRUM_FILENAME, rows, sweep.SPECTRUM_FIELDS)
    with pytest.raises(RuntimeError, match="incomplete .*data"):
        sweep.create_plot_from_csv(tmp_path)


def test_benchmark_resume_recovers_anchor_metrics_from_checkpoint_records(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    all_rows = _synthetic_rows()
    anchor_keys = {
        (sweep.CONFIG_BOTH_OUTER, 0.0),
        (sweep.CONFIG_BOTH_INNER, 0.0),
        (sweep.CONFIG_ANTI_PHASE, 0.0),
        (sweep.CONFIG_BOTH_OUTER, 0.8),
        (sweep.CONFIG_ANTI_PHASE, 0.8),
    }
    rows = [
        row
        for row in all_rows
        if (row["configuration_id"], round(float(row["chi"]), 10)) in anchor_keys
    ]
    records = [
        {
            "configuration_id": sweep.CONFIG_BOTH_OUTER,
            "chi": 0.0,
            "benchmark": True,
            "wall_time_seconds": 1.0,
        },
        {
            "configuration_id": sweep.CONFIG_BOTH_OUTER,
            "chi": 0.8,
            "benchmark": True,
            "wall_time_seconds": 2.0,
        },
        {
            "configuration_id": sweep.CONFIG_ANTI_PHASE,
            "chi": 0.8,
            "benchmark": True,
            "wall_time_seconds": 1.5,
        },
    ]
    monkeypatch.setattr(
        sweep,
        "solve_point",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("completed benchmark anchor was recalculated")
        ),
    )
    _rows, benchmark = sweep.run_benchmarks(
        tmp_path,
        rows,
        records,
        {"status": "PASS"},
        "2026-08-28T00:00:00+00:00",
    )
    assert len(benchmark["anchors"]) == 3
    assert benchmark["conservative_eta_seconds"] == pytest.approx(354.0)


def test_refresh_rejects_incompatible_legacy_manifest(tmp_path: Path) -> None:
    rows = _synthetic_rows()
    (tmp_path / sweep.MANIFEST_FILENAME).write_text(
        json.dumps(
            {
                "stage_id": sweep.STAGE_ID,
                "algorithm_version": "legacy_v1",
                "contract_sha256": sweep.contract_hash(),
                "root_contract": {
                    "search_policy_requested_roots": 12,
                    "search_policy_guard_roots": 1,
                    "guard_position": 13,
                    "root9_plotted": False,
                },
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError, match="incompatible legacy provenance"):
        sweep.refresh_completed_outputs(tmp_path, rows)


def test_analysis_import_boundary_and_no_local_characteristic_solver() -> None:
    path = (
        ROOT
        / "scripts"
        / "analysis"
        / "laminated_beams"
        / "sweep_reddy_stiffness_layout_contrast.py"
    )
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imported: list[str] = []
    functions: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.FunctionDef):
            functions.add(node.name)
    text = " ".join(imported).lower()
    for forbidden in (
        "spectral_sweep_runner",
        "ritz",
        "fem",
        "yartsev",
        "shellbuckling",
    ):
        assert forbidden not in text
    assert not any("characteristic_determinant" in name for name in functions)
    assert "determinant" not in functions


def test_manifest_contains_complete_frequency_map_policy_instance() -> None:
    contract = sweep.contract_payload()
    policy = contract["frequency_map_policy"]
    assert policy == {
        "frequency_map_policy": "frequency-map-v1",
        "calculation_mode": "fast_plot",
        "spectrum_semantics": "sorted_positions",
        "sweep_parameter": "chi",
        "parameter_grid": "0.00:0.02:0.80",
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


def test_generated_outputs_if_present() -> None:
    manifest_path = RESULT_DIR / sweep.MANIFEST_FILENAME
    if not manifest_path.is_file():
        pytest.skip("RLB-2E production outputs have not been generated yet.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows = sweep._read_csv(RESULT_DIR / sweep.SPECTRUM_FILENAME)
    assert sweep.audit_spectrum_rows(rows)["status"] == "PASS"
    assert manifest["counts"]["base_rows"] == 1107
    assert manifest["counts"]["base_configuration_points"] == 123
    assert manifest["counts"]["root9_guards"] == 123
    assert manifest["root_contract"]["roots_above_9_computed"] is False
    assert manifest["root_contract"]["search_policy_requested_roots"] == 8
    assert manifest["root_contract"]["search_policy_guard_roots"] == 1
    assert manifest["root_contract"]["accepted_candidates_above_root9"] == 0
    assert manifest["root_contract"]["root9_plotted"] is False
    assert all(item["status"] == "PASS" for item in manifest["arm_swap_checks"])
    for name, expected in manifest["output_hashes"].items():
        assert sweep._sha256(RESULT_DIR / name) == expected
    assert (RESULT_DIR / sweep.PLOT_FILENAME).is_file()
