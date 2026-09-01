from __future__ import annotations

import ast
import json
import math
from pathlib import Path
import shutil
from typing import Any

import numpy as np
from numpy.testing import assert_allclose
import pytest

from scripts.analysis.laminated_beams import (
    sweep_reddy_axial_stiffness_visibility as sweep,
)


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = (
    ROOT
    / "scripts"
    / "analysis"
    / "laminated_beams"
    / "sweep_reddy_axial_stiffness_visibility.py"
)
RESULT_DIR = (
    ROOT / "results" / "laminated_beams" / "reddy_axial_stiffness_visibility"
)


def _synthetic_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for beta_deg in sweep.BETA_VALUES_DEG:
        for alpha_A in sweep.alpha_grid():
            alpha = float(alpha_A)
            for position in range(1, sweep.K_GUARD + 1):
                Omega = 10.0 * position + 0.5 * beta_deg + 0.2 * alpha
                rows.append(
                    {
                        "row_id": f"beta{int(beta_deg)}_{alpha:.2f}_{position}",
                        "beta_deg": beta_deg,
                        "alpha_A": alpha,
                        "grid_kind": "BASE",
                        "continuation_leg": "SYNTHETIC",
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
                        "cluster_id": f"synthetic_{beta_deg:.0f}_{alpha:.2f}_{position}",
                        "cluster_semantics": "ISOLATED",
                        "cluster_multiplicity": 1,
                        "cluster_total_nullity": 1,
                        "cluster_center_Omega": Omega,
                        "cluster_metadata_source": "SYNTHETIC",
                        "unresolved_candidates_below_root9": 0,
                        "search_right_Omega": 95.0,
                        "root9_right_margin_Omega": 5.0,
                        "guard_cluster_multiplicity": (
                            1 if position == sweep.K_GUARD else ""
                        ),
                        "guard_cluster_extends_beyond_export": False,
                        "solve_mode": "SYNTHETIC",
                        "fallback_used": False,
                        "quality_status": "PASS",
                        "point_status": "PASS",
                        "is_canonical_plot_source": True,
                        "supersedes_row_id": "",
                        "repair_id": "",
                        "decoupled_reference_classification": "",
                        "roots_above_9_computed": False,
                    }
                )
    return rows


def _synthetic_solution(
    rows: list[dict[str, Any]], beta_deg: float, alpha_A: float
) -> sweep.PointSolution:
    selected = [
        dict(row)
        for row in rows
        if float(row["beta_deg"]) == pytest.approx(beta_deg)
        and float(row["alpha_A"]) == pytest.approx(alpha_A)
        and str(row["grid_kind"]) == "BASE"
    ]
    for row in selected:
        row["row_id"] = str(row["row_id"]) + "_repair"
        row["grid_kind"] = "LOCAL_REFINEMENT"
        row["is_canonical_plot_source"] = True
    return sweep.PointSolution(
        beta_deg=beta_deg,
        alpha_A=alpha_A,
        rows=tuple(selected),
        wall_time_seconds=0.01,
        peak_rss_bytes=1024,
        determinant_evaluations=9,
        sigma_evaluations=9,
        search_left_Omega=0.0,
        search_right_Omega=95.0,
        local_refinements=9,
        solve_mode="SYNTHETIC_LOCAL_REPAIR",
        fallback_used=False,
        unresolved_candidates_below_root9=0,
        continuation_leg="LOCAL_REPAIR",
    )


def _write_spectrum(path: Path, rows: list[dict[str, Any]]) -> None:
    sweep.rlb2e._atomic_write_csv(path, rows, sweep.SPECTRUM_FIELDS)


def test_01_exact_geometry_contract_and_beta_values() -> None:
    assert sweep.BETA_VALUES_DEG == (0.0, 30.0)
    assert (sweep.MU, sweep.TAU, sweep.L1, sweep.L2, sweep.L_TOTAL) == (
        0.0,
        0.0,
        1.0,
        1.0,
        2.0,
    )
    assert (sweep.WIDTH, sweep.THICKNESS, sweep.PLY_THICKNESS) == (
        0.20,
        0.05,
        0.0125,
    )
    assert sweep.K == pytest.approx(5.0 / 6.0)


def test_02_exact_base_material_contract() -> None:
    assert sweep.base_material_contract() == {
        "E1_0": 1.1,
        "E2_0": 0.9,
        "nu12_0": 0.3,
        "G12_0": 1.0 / 2.6,
        "G13_0": 1.0 / 2.6,
        "G23_0": 1.0 / 2.6,
        "rho_0": 1.0,
    }


def test_03_exact_alpha_grid_and_continuation_paths() -> None:
    assert_allclose(
        sweep.alpha_grid(),
        np.linspace(0.70, 1.30, 31),
        rtol=0.0,
        atol=0.0,
    )
    lower, upper = sweep.continuation_paths()
    assert lower[0] == pytest.approx(1.0)
    assert lower[-1] == pytest.approx(0.70)
    assert upper[0] == pytest.approx(1.0)
    assert upper[-1] == pytest.approx(1.30)


@pytest.mark.parametrize("alpha_A", [0.70, 1.00, 1.30])
def test_04_outer_inner_formulas_and_positive_multipliers(alpha_A: float) -> None:
    s_outer, s_inner = sweep.alpha_scales(alpha_A)
    assert s_outer == pytest.approx((4.0 - alpha_A) / 3.0)
    assert s_inner == pytest.approx((7.0 * alpha_A - 4.0) / 3.0)
    assert s_outer > 0.0
    assert s_inner > 0.0
    assert (s_outer + s_inner) / 2.0 == pytest.approx(alpha_A)
    assert (7.0 * s_outer + s_inner) / 8.0 == pytest.approx(1.0)


@pytest.mark.parametrize("alpha_A", [0.70, 1.00, 1.30])
def test_05_only_in_plane_moduli_scale(alpha_A: float) -> None:
    s_outer, s_inner = sweep.alpha_scales(alpha_A)
    for scale in (s_outer, s_inner):
        material = sweep.build_scaled_material(scale, "test")
        assert material.E1 == pytest.approx(scale * 1.1)
        assert material.E2 == pytest.approx(scale * 0.9)
        assert material.G12 == pytest.approx(scale / 2.6)
        assert material.nu12 == pytest.approx(0.3)
        assert material.G13 == pytest.approx(1.0 / 2.6)
        assert material.G23 == pytest.approx(1.0 / 2.6)
        assert material.rho == pytest.approx(1.0)


def test_06_Q_scales_but_transverse_shear_Q_is_invariant() -> None:
    beam, _coupled = sweep._physics_modules()
    base = sweep.build_scaled_material(1.0, "base")
    scaled = sweep.build_scaled_material(1.3, "scaled")
    assert_allclose(
        beam.lamina_reduced_stiffness(scaled),
        1.3 * beam.lamina_reduced_stiffness(base),
        rtol=1.0e-12,
        atol=1.0e-16,
    )
    assert_allclose(
        beam.lamina_transverse_shear_stiffness(scaled),
        beam.lamina_transverse_shear_stiffness(base),
        rtol=1.0e-12,
        atol=1.0e-16,
    )


def test_07_exact_four_equal_zero_degree_plies() -> None:
    section = sweep.build_layered_section(1.2)
    assert section.layout == sweep.LAYOUT
    assert len(section.laminate.plies) == 4
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


@pytest.mark.parametrize("alpha_A", [0.70, 1.00, 1.30])
def test_08_constitutive_A_D_S_m_J_contract(alpha_A: float) -> None:
    baseline = sweep._baseline_section()
    section = sweep.build_layered_section(alpha_A)
    assert_allclose(
        section.laminate.A,
        alpha_A * baseline.laminate.A,
        rtol=sweep.MATRIX_RELATIVE_TOLERANCE,
        atol=1.0e-16,
    )
    assert_allclose(
        section.laminate.D,
        baseline.laminate.D,
        rtol=sweep.MATRIX_RELATIVE_TOLERANCE,
        atol=1.0e-16,
    )
    assert_allclose(
        section.laminate.shear,
        baseline.laminate.shear,
        rtol=sweep.MATRIX_RELATIVE_TOLERANCE,
        atol=1.0e-16,
    )
    assert section.properties.A / baseline.properties.A == pytest.approx(
        alpha_A, rel=1.0e-11
    )
    assert section.properties.D / baseline.properties.D == pytest.approx(
        1.0, rel=1.0e-11
    )
    assert section.properties.S / baseline.properties.S == pytest.approx(
        1.0, rel=1.0e-11
    )
    assert section.properties.m / baseline.properties.m == pytest.approx(
        1.0, rel=1.0e-11
    )
    assert section.properties.J / baseline.properties.J == pytest.approx(
        1.0, rel=1.0e-11
    )


def test_09_B_I1_and_beamproperties_only_A_changes() -> None:
    baseline = sweep._baseline_section()
    maxima = []
    for alpha_A in sweep.alpha_grid():
        section = sweep.build_layered_section(float(alpha_A))
        assert sweep._scaled_B(section.laminate) <= 1.0e-12
        assert sweep._scaled_I1(section.laminate) <= 1.0e-12
        maxima.append(
            max(
                sweep._relative(section.properties.D, baseline.properties.D),
                sweep._relative(section.properties.S, baseline.properties.S),
                sweep._relative(section.properties.m, baseline.properties.m),
                sweep._relative(section.properties.J, baseline.properties.J),
                sweep._relative(section.properties.K, baseline.properties.K),
                sweep._relative(section.properties.width, baseline.properties.width),
            )
        )
    assert max(maxima) <= 1.0e-11


def test_10_constitutive_gate_passes_full_grid() -> None:
    gate = sweep.constitutive_gate()
    assert gate["status"] == "PASS"
    assert all(item["status"] == "PASS" for item in gate["checks"])


def test_11_section_rows_have_exact_31_parameter_rows() -> None:
    rows = sweep.section_property_rows()
    assert len(rows) == 31
    assert all(row["properties_identical_between_arms"] for row in rows)
    assert all(row["constitutive_status"] == "PASS" for row in rows)


def test_12_frequency_map_policy_and_root_contract_are_complete() -> None:
    contract = sweep.contract_payload()
    policies = contract["frequency_map_policy_instances"]
    assert set(policies) == {"beta_0", "beta_30"}
    common = {
        "frequency_map_policy": "frequency-map-v1",
        "calculation_mode": "fast_plot",
        "spectrum_semantics": "sorted_positions",
        "sweep_parameter": "alpha_A",
        "parameter_grid": "0.70:0.02:1.30",
        "continuation_anchor": "1.00",
        "continuation_paths": ["1.00:0.98:0.70", "1.00:1.02:1.30"],
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
    for policy in policies.values():
        for key, expected in common.items():
            assert policy[key] == expected
    assert policies["beta_0"]["beta_deg"] == 0.0
    assert policies["beta_30"]["beta_deg"] == 30.0
    root_contract = contract["root_contract"]
    assert root_contract["plotted_positions"] == list(range(1, 9))
    assert root_contract["guard_position"] == 9
    assert root_contract["root9_plotted"] is False
    assert root_contract["roots_above_9_computed"] is False


def test_13_old_normalization_is_exact() -> None:
    expected = sweep.L_REFERENCE**2 * math.sqrt(12.0 / sweep.THICKNESS**2)
    assert sweep.OMEGA_TO_OMEGA_SCALE == pytest.approx(expected)
    omega = 2.75
    Omega = sweep.omega_to_Omega(omega)
    Lambda = sweep.Omega_to_Lambda(Omega)
    assert_allclose(Lambda * Lambda, Omega, rtol=2.0e-16, atol=0.0)


def test_14_root_inventory_audit_requires_62_groups_and_558_base_rows() -> None:
    rows = _synthetic_rows()
    audit = sweep.audit_spectrum_rows(rows)
    assert audit["status"] == "PASS"
    assert audit["base_group_count"] == 62
    assert audit["base_row_count"] == 558
    assert audit["roots_above_guard_count"] == 0
    assert len({str(row["row_id"]) for row in rows}) == len(rows)
    for group in sweep._base_group_index(rows).values():
        ordered = sorted(group, key=lambda row: int(row["sorted_position"]))
        assert [int(row["sorted_position"]) for row in ordered] == list(range(1, 10))
        assert np.all(np.diff([float(row["Omega"]) for row in ordered]) > 0.0)
        assert all(int(row["cluster_multiplicity"]) == 1 for row in ordered)
    assert sweep.audit_spectrum_rows(rows[:-1])["status"] == "FAIL"
    position10 = dict(rows[-1])
    position10["row_id"] += "_p10"
    position10["sorted_position"] = 10
    position10["Omega"] = float(position10["Omega"]) + 10.0
    position10["Lambda"] = math.sqrt(float(position10["Omega"]))
    assert sweep.audit_spectrum_rows([*rows, position10])["status"] == "FAIL"


def test_15_root9_is_guard_and_excluded_from_plot_data() -> None:
    rows = _synthetic_rows()
    guards = [row for row in rows if int(row["sorted_position"]) == 9]
    assert len(guards) == 62
    assert all(row["root_role"] == "ROOT_9_GUARD" for row in guards)
    assert all(row["guard_flag"] for row in guards)
    plot_rows = sweep.canonical_plot_rows(rows)
    assert len(plot_rows) == 2 * 31 * 8
    assert all(int(row["sorted_position"]) <= 8 for row in plot_rows)
    assert all(not row["predictor_used_as_final"] for row in rows)
    assert all(not row["roots_above_9_computed"] for row in rows)


def test_15b_resume_schema_backfill_changes_metadata_only() -> None:
    rows = _synthetic_rows()
    before = [float(row["Omega"]) for row in rows]
    for row in rows:
        row["cluster_metadata_source"] = ""
        if int(row["sorted_position"]) == sweep.K_GUARD:
            row["guard_cluster_multiplicity"] = ""
            row["guard_cluster_extends_beyond_export"] = ""
    summary = sweep.backfill_resumed_spectrum_metadata(rows)
    assert [float(row["Omega"]) for row in rows] == before
    assert summary == {
        "cluster_metadata_source_backfilled_rows": 558,
        "guard_metadata_backfilled_groups": 62,
        "root_values_modified": 0,
    }
    guards = [row for row in rows if int(row["sorted_position"]) == 9]
    assert all(int(row["guard_cluster_multiplicity"]) == 1 for row in guards)
    assert all(not row["guard_cluster_extends_beyond_export"] for row in guards)


def test_15c_distinct_root_above_guard_cannot_be_silently_truncated() -> None:
    from types import SimpleNamespace

    slots = [
        SimpleNamespace(
            event=SimpleNamespace(
                omega_bar=float(index),
                cluster_id="",
                event_id=f"root_{index}",
            )
        )
        for index in range(1, 11)
    ]
    with pytest.raises(RuntimeError, match="distinct root above the guard"):
        sweep._truncate_inventory_to_root9([], slots, sweep._policy_root9())


def test_15d_unexpected_canonical_parameter_key_is_rejected() -> None:
    rows = _synthetic_rows()
    unexpected = dict(rows[0])
    unexpected["row_id"] = "unexpected_alpha_grid_row"
    unexpected["alpha_A"] = 0.71
    unexpected["grid_kind"] = "LOCAL_REFINEMENT"
    rows.append(unexpected)
    assert sweep.audit_spectrum_rows(rows)["status"] == "FAIL"
    assert sweep.audit_plot_rows(rows)["status"] == "FAIL"
    unexpected["is_canonical_plot_source"] = False
    assert sweep.audit_spectrum_rows(rows)["status"] == "FAIL"
    assert sweep.audit_plot_rows(rows)["status"] == "FAIL"


def test_15e_completed_schema_requires_provenance_fields() -> None:
    rows = _synthetic_rows()
    for row in rows:
        if float(row["beta_deg"]) == 0.0:
            row["decoupled_reference_classification"] = "AXIAL"
    assert sweep.audit_final_spectrum_schema(rows)["status"] == "PASS"
    rows[0]["cluster_metadata_source"] = ""
    assert sweep.audit_final_spectrum_schema(rows)["status"] == "FAIL"


def test_16_neighbour_audit_marks_no_flag_for_smooth_synthetic_data() -> None:
    audit_rows = sweep.neighbour_audit_rows(_synthetic_rows())
    assert audit_rows
    assert not any(bool(row["flagged"]) for row in audit_rows)


def test_16b_cluster_aware_comparison_preserves_multiplicity_and_nullity() -> None:
    def rows(total_nullity: int = 2) -> list[dict[str, Any]]:
        return [
            {
                "sorted_position": position,
                "Omega": 10.0,
                "cluster_id": "exact_double",
                "cluster_semantics": "EXACT_DEGENERATE_SUBSPACE",
                "cluster_multiplicity": 2,
                "cluster_total_nullity": total_nullity,
                "cluster_center_Omega": 10.0,
            }
            for position in (1, 2)
        ]

    records, summary = sweep._compare_cluster_groups(
        rows(), rows(), alpha_A=1.0, comparison_kind="synthetic"
    )
    assert summary["status"] == "PASS"
    assert records[0]["comparison_unit"] == "CLUSTER"
    _records, bad = sweep._compare_cluster_groups(
        rows(), rows(total_nullity=1), alpha_A=1.0, comparison_kind="synthetic"
    )
    assert bad["status"] == "FAIL"


def test_16c_neighbour_spike_is_flagged_without_smoothing() -> None:
    rows = _synthetic_rows()
    target = next(
        row
        for row in rows
        if float(row["beta_deg"]) == 30.0
        and float(row["alpha_A"]) == pytest.approx(1.0)
        and int(row["sorted_position"]) == 4
    )
    target["Omega"] = 80.0
    target["omega"] = 80.0 / sweep.OMEGA_TO_OMEGA_SCALE
    target["Lambda"] = math.sqrt(80.0)
    audit = sweep.neighbour_audit_rows(rows)
    assert any(
        float(row["beta_deg"]) == 30.0
        and float(row["alpha_A"]) == pytest.approx(1.0)
        and bool(row["flagged"])
        for row in audit
    )
    assert all(row["smoothing_applied"] is False for row in audit)


def test_16d_only_flagged_point_is_recomputed_and_reproduced(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = _synthetic_rows()
    audit = [
        {"beta_deg": 0.0, "alpha_A": 0.8, "sorted_position": 4, "flagged": False},
        {"beta_deg": 30.0, "alpha_A": 1.0, "sorted_position": 4, "flagged": True},
    ]
    calls: list[tuple[float, float]] = []

    def fake_solve(beta_deg: float, alpha_A: float, **_kwargs: Any) -> sweep.PointSolution:
        calls.append((float(beta_deg), float(alpha_A)))
        return _synthetic_solution(rows, beta_deg, alpha_A)

    monkeypatch.setattr(sweep, "solve_point", fake_solve)
    repaired, updated, records = sweep.apply_local_repairs(rows, audit)
    assert calls == [(30.0, 1.0)]
    assert records[0]["status"] == "REPRODUCED_AFTER_LOCAL_REPAIR"
    assert records[0]["smoothing_applied"] is False
    assert updated[1]["repair_status"] == "REPRODUCED_AFTER_LOCAL_REPAIR"
    assert sweep.audit_spectrum_rows(repaired)["status"] == "PASS"


def test_16e_unresolved_local_repair_creates_visible_nan_gap(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows = _synthetic_rows()
    audit = [
        {"beta_deg": 30.0, "alpha_A": 1.0, "sorted_position": 3, "flagged": True}
    ]

    def fail_solve(*_args: Any, **_kwargs: Any) -> sweep.PointSolution:
        raise RuntimeError("synthetic unresolved root neighbourhood")

    monkeypatch.setattr(sweep, "solve_point", fail_solve)
    repaired, updated, records = sweep.apply_local_repairs(rows, audit)
    canonical = sweep._canonical_group(repaired, 30.0, 1.0)
    assert records[0]["status"] == "UNRESOLVED"
    assert updated[0]["repair_status"] == "UNRESOLVED_AFTER_LOCAL_REPAIR"
    assert math.isnan(float(canonical[2]["Lambda"]))
    assert canonical[2]["point_status"] == "UNRESOLVED_AFTER_LOCAL_REPAIR"
    assert sweep.audit_plot_rows(repaired)["status"] == "PASS"
    recovered = sweep._recover_existing_repair_records(repaired, updated)
    assert recovered[0]["status"] == "UNRESOLVED"
    assert recovered[0]["affected_positions"] == [3]


def test_16f_six_benchmark_anchors_and_eta_gate() -> None:
    anchors = [
        (0.0, 1.0),
        (0.0, 0.70),
        (0.0, 1.30),
        (30.0, 1.0),
        (30.0, 0.70),
        (30.0, 1.30),
    ]
    records = [
        {
            "beta_deg": beta,
            "alpha_A": alpha,
            "wall_time_seconds": 1.0,
            "peak_rss_bytes": 1024,
            "determinant_evaluations": 1,
        }
        for beta, alpha in anchors
    ]
    payload = sweep._benchmark_payload(records)
    assert [(row["beta_deg"], row["alpha_A"]) for row in payload["anchors"]] == anchors
    assert payload["anchor_count"] == 6
    assert payload["remaining_unique_root_points"] == 56
    assert payload["eta_limit_seconds"] == 2700.0
    assert payload["production_run_permitted"] is True
    records[0]["wall_time_seconds"] = 100.0
    assert sweep._benchmark_payload(records)["production_run_permitted"] is False


def test_17_plot_only_creates_two_panel_figure_and_never_calls_physics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write_spectrum(tmp_path / sweep.SPECTRUM_FILENAME, _synthetic_rows())

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("plot_only called production physics or root search")

    monkeypatch.setattr(sweep, "solve_point", forbidden)
    monkeypatch.setattr(sweep, "make_matrix_provider", forbidden)
    monkeypatch.setattr(sweep, "build_layered_section", forbidden)
    monkeypatch.setattr(sweep, "_physics_modules", forbidden)
    result = sweep.create_plot_from_csv(tmp_path)
    assert result["panel_count"] == 2
    assert result["lines_per_panel"] == 8
    assert result["root9_plotted"] is False
    assert result["root_calculation_count"] == 0
    assert (tmp_path / sweep.PLOT_FILENAME).is_file()


def test_17b_beta0_diagnostic_plot_uses_csv_only(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reference_rows = []
    for alpha_A in sweep.alpha_grid():
        for position in range(1, sweep.K_GUARD + 1):
            family = "axial" if position == 7 else "bending"
            reference_rows.append(
                {
                    "alpha_A": float(alpha_A),
                    "sorted_position": position,
                    "family": family,
                    "family_index": 1 if family == "axial" else position,
                    "Lambda": math.sqrt(10.0 * position),
                }
            )
    sweep.rlb2e._atomic_write_csv(
        tmp_path / sweep.REFERENCE_FILENAME, reference_rows
    )

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("diagnostic plot called physics or a root solver")

    monkeypatch.setattr(sweep, "solve_point", forbidden)
    monkeypatch.setattr(sweep, "make_matrix_provider", forbidden)
    monkeypatch.setattr(sweep, "build_layered_section", forbidden)
    monkeypatch.setattr(sweep, "_physics_modules", forbidden)
    result = sweep.create_beta0_reference_plot(tmp_path)
    assert result["root_calculation_count"] == 0
    assert (tmp_path / sweep.BETA0_PLOT_FILENAME).stat().st_size > 0


def test_18_manifest_only_uses_zero_root_calculations() -> None:
    payload = sweep.manifest_only()
    assert payload["root_calculation_count"] == 0
    assert payload["section_property_row_count"] == 31


def test_18a_completed_fast_path_rejects_stale_output_hash(
    tmp_path: Path,
) -> None:
    rows = _synthetic_rows()
    for row in rows:
        if float(row["beta_deg"]) == 0.0:
            row["decoupled_reference_classification"] = "AXIAL"
    _write_spectrum(tmp_path / sweep.SPECTRUM_FILENAME, rows)
    required = {
        sweep.SECTION_FILENAME,
        sweep.AUDIT_FILENAME,
        sweep.REFERENCE_FILENAME,
        sweep.BENCHMARK_FILENAME,
        sweep.CHECKPOINT_FILENAME,
        sweep.PLOT_FILENAME,
        sweep.BETA0_PLOT_FILENAME,
        sweep.REPORT_FILENAME,
    }
    for name in required:
        (tmp_path / name).write_bytes(f"synthetic:{name}".encode("utf-8"))
    manifest = {
        "stage_id": sweep.STAGE_ID,
        "contract_sha256": sweep.contract_hash(),
        "analysis_script_sha256": sweep.rlb2e._sha256(SCRIPT_PATH),
        "scientific_status": "PASS",
        "status_gates": {
            name: "PASS" for name in sweep.REQUIRED_STATUS_GATE_NAMES
        },
        "production_physics_preserved": True,
        "predecessor_result_trees_preserved": True,
        "production_physics_hashes": {
            path: sweep.rlb2e._sha256(ROOT / path)
            for path in sweep.PRODUCTION_PHYSICS_PATHS
        },
        "predecessor_result_tree_hashes": {
            name: sweep._sha_tree(path)
            for name, path in sweep.PREDECESSOR_RESULT_DIRS.items()
        },
        "output_hashes": sweep._output_hashes(tmp_path),
    }
    sweep.rlb2e._atomic_write_json(tmp_path / sweep.MANIFEST_FILENAME, manifest)
    assert sweep._results_are_complete(rows, tmp_path)
    (tmp_path / sweep.REPORT_FILENAME).write_text("stale", encoding="utf-8")
    assert not sweep._results_are_complete(rows, tmp_path)


def test_18aa_checkpoint_contract_is_validated_before_output_write(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sweep.rlb2e._atomic_write_json(
        tmp_path / sweep.CHECKPOINT_FILENAME,
        {"contract_sha256": "WRONG", "point_records": []},
    )

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("output preparation ran before checkpoint validation")

    monkeypatch.setattr(sweep, "prepare_constitutive_outputs", forbidden)
    with pytest.raises(RuntimeError, match="Checkpoint contract differs"):
        sweep.run_workflow(tmp_path, missing_only=True)


def test_18ab_checkpoint_record_requires_a_completed_spectrum_group() -> None:
    checkpoint = {
        "schema_version": 1,
        "stage_id": sweep.STAGE_ID,
        "algorithm_version": sweep.ALGORITHM_VERSION,
        "contract_sha256": sweep.contract_hash(),
        "point_records": [
            {
                "beta_deg": 0.0,
                "alpha_A": 1.0,
                "roots": [
                    {"sorted_position": position}
                    for position in range(1, sweep.K_GUARD + 1)
                ],
            }
        ],
    }
    with pytest.raises(RuntimeError, match="completed spectrum groups"):
        sweep._validated_checkpoint_point_records(checkpoint, [])


def test_18ac_checkpoint_root_values_match_stored_base_spectrum() -> None:
    rows = _synthetic_rows()
    group = [
        row
        for row in rows
        if float(row["beta_deg"]) == 0.0
        and float(row["alpha_A"]) == pytest.approx(1.0)
    ]
    group.sort(key=lambda row: int(row["sorted_position"]))
    record = {
        "beta_deg": 0.0,
        "alpha_A": 1.0,
        "benchmark": False,
        "continuation_leg": "SYNTHETIC",
        "wall_time_seconds": 1.0,
        "peak_rss_bytes": 1,
        "determinant_evaluations": 9,
        "sigma_evaluations": 9,
        "solve_mode": "SYNTHETIC",
        "fallback_used": False,
        "unresolved_candidates_below_root9": 0,
        "roots": [
            {
                "sorted_position": int(row["sorted_position"]),
                "Omega": float(row["Omega"]),
                "Lambda": float(row["Lambda"]),
                "singular_ratio": float(row["scaled_sigma_ratio"]),
                "boundary_residual": float(row["boundary_null_residual"]),
                "quality_status": str(row["quality_status"]),
            }
            for row in group
        ],
    }
    checkpoint = {
        "schema_version": 1,
        "stage_id": sweep.STAGE_ID,
        "algorithm_version": sweep.ALGORITHM_VERSION,
        "contract_sha256": sweep.contract_hash(),
        "point_records": [record],
    }
    assert len(sweep._validated_checkpoint_point_records(checkpoint, rows)) == 1
    record["roots"][0]["Omega"] *= 1.01
    with pytest.raises(RuntimeError, match="disagrees with the stored BASE"):
        sweep._validated_checkpoint_point_records(checkpoint, rows)


@pytest.mark.skipif(
    not (RESULT_DIR / sweep.MANIFEST_FILENAME).is_file(),
    reason="RLB-2H final manifest not generated",
)
def test_18b_completed_missing_only_is_bytewise_idempotent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    before = {
        path.name: sweep.rlb2e._sha256(path)
        for path in RESULT_DIR.iterdir()
        if path.is_file()
    }

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("completed missing-only called a calculation path")

    monkeypatch.setattr(sweep, "solve_point", forbidden)
    monkeypatch.setattr(sweep, "make_matrix_provider", forbidden)
    monkeypatch.setattr(sweep, "build_layered_section", forbidden)
    result = sweep.run_workflow(RESULT_DIR, missing_only=True)
    after = {
        path.name: sweep.rlb2e._sha256(path)
        for path in RESULT_DIR.iterdir()
        if path.is_file()
    }
    assert after == before
    assert result["resume_root_calculation_count"] == 0
    assert result["resume_outputs_modified"] is False


@pytest.mark.skipif(
    not (RESULT_DIR / sweep.MANIFEST_FILENAME).is_file(),
    reason="RLB-2H final manifest not generated",
)
def test_18c_plot_only_is_zero_physics_and_rejects_modified_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    target = tmp_path / "copied_result"
    shutil.copytree(RESULT_DIR, target)

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("plot_only called a physics or root-calculation path")

    monkeypatch.setattr(sweep, "solve_point", forbidden)
    monkeypatch.setattr(sweep, "make_matrix_provider", forbidden)
    monkeypatch.setattr(sweep, "build_layered_section", forbidden)
    result = sweep.plot_only_workflow(target)
    assert result["root_calculation_count"] == 0
    assert result["matrix_assembly_count"] == 0
    assert result["determinant_evaluation_count"] == 0
    assert result["SVD_evaluation_count"] == 0
    before = {
        name: sweep.rlb2e._sha256(target / name)
        for name in (sweep.PLOT_FILENAME, sweep.BETA0_PLOT_FILENAME)
    }
    (target / sweep.REPORT_FILENAME).write_text("modified", encoding="utf-8")
    with pytest.raises(RuntimeError, match="modified evidence"):
        sweep.plot_only_workflow(target)
    after = {
        name: sweep.rlb2e._sha256(target / name)
        for name in (sweep.PLOT_FILENAME, sweep.BETA0_PLOT_FILENAME)
    }
    assert after == before


def test_19_no_forbidden_imports() -> None:
    source = SCRIPT_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    imported_modules: list[str] = []
    locally_defined: set[str] = set()
    calls: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported_modules.append(node.module or "")
        elif isinstance(node, ast.FunctionDef):
            locally_defined.add(node.name)
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                calls.add(node.func.id)
            elif isinstance(node.func, ast.Attribute):
                calls.add(node.func.attr)
    text = "\n".join(imported_modules)
    forbidden = (
        "spectral_sweep_runner",
        "Ritz",
        "ritz",
        "FEM",
        "fem",
        "analytic_branch_tracking",
    )
    assert all(item not in text for item in forbidden)
    assert not any("characteristic_determinant" in name for name in locally_defined)
    assert not ({"interp", "interp1d", "savgol_filter", "UnivariateSpline"} & calls)
    assert '"smoothing": False' in source
    assert '"interpolation_based_frequencies": False' in source


@pytest.mark.skipif(
    not (RESULT_DIR / sweep.MANIFEST_FILENAME).is_file(),
    reason="RLB-2H final manifest not generated",
)
def test_20_generated_outputs_are_complete_and_consistent() -> None:
    manifest = json.loads((RESULT_DIR / sweep.MANIFEST_FILENAME).read_text(encoding="utf-8"))
    rows = sweep.rlb2e._read_csv(RESULT_DIR / sweep.SPECTRUM_FILENAME)
    audit = sweep.audit_spectrum_rows(rows)
    required = {
        sweep.SECTION_FILENAME,
        sweep.SPECTRUM_FILENAME,
        sweep.AUDIT_FILENAME,
        sweep.REFERENCE_FILENAME,
        sweep.BENCHMARK_FILENAME,
        sweep.CHECKPOINT_FILENAME,
        sweep.PLOT_FILENAME,
        sweep.BETA0_PLOT_FILENAME,
        sweep.REPORT_FILENAME,
        sweep.MANIFEST_FILENAME,
    }
    assert required <= {path.name for path in RESULT_DIR.iterdir() if path.is_file()}
    assert audit["status"] == "PASS"
    assert manifest["scientific_status"] == "PASS"
    assert set(manifest["status_gates"].values()) == {"PASS"}
    assert manifest["counts"]["beta0_points_complete"] == 31
    assert manifest["counts"]["beta30_points_complete"] == 31
    assert manifest["counts"]["base_rows"] == 558
    assert manifest["counts"]["base_groups"] == 62
    assert manifest["counts"]["root9_guards"] == 62
    assert len(sweep.rlb2e._read_csv(RESULT_DIR / sweep.SECTION_FILENAME)) == 31
    reference_rows = sweep.rlb2e._read_csv(RESULT_DIR / sweep.REFERENCE_FILENAME)
    assert len(reference_rows) == 31 * 9
    reference = manifest["beta0_reference"]
    assert reference["status"] == "PASS"
    assert reference["bending_reference_source_group_count"] == 1
    assert reference["bending_reference_root_count"] == 8
    assert reference["root_solver_invocations"] == 0
    assert reference["bending_invariance_error"] <= sweep.DIRECT_REFERENCE_RELATIVE_TOLERANCE
    assert reference["axial_scaling_error"] <= sweep.DIRECT_REFERENCE_RELATIVE_TOLERANCE
    direct = manifest["beta0_direct_coupled_check"]
    assert direct["status"] == "PASS"
    assert direct["direct_reference_check_count"] == 3
    assert direct["direct_reference_solve_count"] == 0
    assert {record["alpha_A"] for record in direct["direct_reference_records"]} == {
        0.7,
        1.0,
        1.3,
    }
    assert all(
        len(record["detected_nullities"]) == sweep.K_GUARD
        for record in direct["direct_reference_records"]
    )
    assert all(row["status"] == "PASS" for row in direct["comparisons"])
    all_grid = manifest["beta0_all_grid_coupled_subsystem_check"]
    assert all_grid["status"] == "PASS"
    assert all_grid["point_count"] == 31
    assert all_grid["unclassified_canonical_slot_count"] == 0
    neighbour = manifest["neighbour_audit"]
    assert neighbour["unresolved_point_count"] == 0
    assert neighbour["smoothing_applied"] is False
    assert len(manifest["beta30_shift_summary"]) == 8
    assert len(manifest["minimum_adjacent_sorted_gaps"]) == 2
    assert manifest["production_physics_preserved"] is True
    assert manifest["predecessor_result_trees_preserved"] is True
    assert manifest["analysis_script_sha256"] == sweep.rlb2e._sha256(SCRIPT_PATH)
    for name, expected in manifest["output_hashes"].items():
        assert sweep.rlb2e._sha256(RESULT_DIR / name) == expected
    for filename in (sweep.PLOT_FILENAME, sweep.BETA0_PLOT_FILENAME, sweep.REPORT_FILENAME):
        assert (RESULT_DIR / filename).stat().st_size > 0
    report = (RESULT_DIR / sweep.REPORT_FILENAME).read_text(encoding="utf-8")
    for position in range(1, 9):
        assert f"`k={position}`" in report
    assert "RLB-2H-BETA0-DIRECT-SUBSYSTEM-REFERENCE" in report
    assert "**PASS**" in report
