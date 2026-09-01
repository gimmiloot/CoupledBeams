"""Targeted gates for RLB-2G mass-layout duality.

Physical unit tests assemble laminate sections but never search for roots.
Spectrum, resume, repair, and plotting tests use synthetic root inventories.
The production-output audit runs only when the result tree already exists.
"""

from __future__ import annotations

import ast
import builtins
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
from numpy.testing import assert_allclose
import pytest

from scripts.analysis.laminated_beams import sweep_reddy_mass_layout_duality as sweep


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = (
    ROOT
    / "scripts"
    / "analysis"
    / "laminated_beams"
    / "sweep_reddy_mass_layout_duality.py"
)
RESULT_DIR = ROOT / "results" / "laminated_beams" / "reddy_mass_layout_duality"


def _synthetic_row(
    *,
    experiment_id: str,
    configuration_id: str,
    parameter_name: str,
    parameter_value: float,
    position: int,
    Omega: float,
) -> dict[str, Any]:
    shared = abs(parameter_value) < 1.0e-15
    return {
        "row_id": (
            f"{experiment_id}__{configuration_id}__{parameter_name}_"
            f"{parameter_value:+.6f}__p{position:02d}"
        ),
        "experiment_id": experiment_id,
        "configuration_id": configuration_id,
        "parameter_name": parameter_name,
        "parameter_value": float(parameter_value),
        "eta": float(parameter_value) if experiment_id == sweep.EXPERIMENT_A else "",
        "xi_rho": (
            float(parameter_value) if experiment_id == sweep.EXPERIMENT_B else ""
        ),
        "grid_kind": "BASE",
        "continuation_leg": "SYNTHETIC",
        "solve_id": (
            "RLB-2G_SHARED_ZERO_CLONE"
            if shared
            else f"synthetic_{experiment_id}_{configuration_id}_{parameter_value:+.2f}"
        ),
        "physical_solve_id": (
            "RLB-2G_SHARED_ZERO_PHYSICAL_SOLVE"
            if shared
            else f"physical_{experiment_id}_{configuration_id}_{parameter_value:+.2f}"
        ),
        "transaction_id": f"synthetic_{experiment_id}_{configuration_id}_{parameter_value:+.2f}",
        "sorted_position": position,
        "root_role": "PLOTTED" if position <= sweep.K_PLOT else "ROOT_9_GUARD",
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
        "shared_zero_contrast_anchor_reused": shared,
        "shared_zero_contrast_source": "RLB-2G_SHARED_ZERO" if shared else "",
        "is_canonical_plot_source": True,
        "supersedes_row_id": "",
        "repair_id": "",
        "roots_above_9_computed": False,
    }


def _synthetic_A_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for configuration_index, configuration_id in enumerate(sweep.CONFIGURATIONS_A):
        for eta_value in sweep.eta_grid():
            eta = float(eta_value)
            for position in range(1, sweep.K_GUARD + 1):
                # Configuration shift vanishes at the common physical anchor.
                Omega = 10.0 * position + 0.2 * configuration_index * eta
                rows.append(
                    _synthetic_row(
                        experiment_id=sweep.EXPERIMENT_A,
                        configuration_id=configuration_id,
                        parameter_name="eta",
                        parameter_value=eta,
                        position=position,
                        Omega=Omega,
                    )
                )
    return rows


def _synthetic_B_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for xi_value in sweep.xi_rho_grid():
        xi_rho = float(xi_value)
        for position in range(1, sweep.K_GUARD + 1):
            Omega = 10.0 * position + 0.2 * xi_rho
            rows.append(
                _synthetic_row(
                    experiment_id=sweep.EXPERIMENT_B,
                    configuration_id=sweep.CONFIG_ONE_ARM,
                    parameter_name="xi_rho",
                    parameter_value=xi_rho,
                    position=position,
                    Omega=Omega,
                )
            )
    return rows


def _synthetic_solution(
    rows: list[dict[str, Any]],
    experiment_id: str,
    configuration_id: str,
    parameter_value: float,
    *,
    Omega_shift: float = 0.0,
    swapped: bool = False,
) -> sweep.PointSolution:
    selected: list[dict[str, Any]] = []
    for source in rows:
        if (
            str(source["configuration_id"]) == configuration_id
            and round(float(source["parameter_value"]), 10)
            == round(float(parameter_value), 10)
            and str(source["grid_kind"]) == "BASE"
        ):
            row = dict(source)
            position = int(row["sorted_position"])
            Omega = float(row["Omega"]) + Omega_shift
            row.update(
                {
                    "row_id": (
                        f"{experiment_id}__{configuration_id}__"
                        f"synthetic_refinement__{parameter_value:+.6f}__p{position:02d}"
                    ),
                    "grid_kind": (
                        "ARM_SWAP_DIAGNOSTIC" if swapped else "LOCAL_REFINEMENT"
                    ),
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
        experiment_id=experiment_id,
        configuration_id=configuration_id,
        parameter_name=sweep._parameter_name(experiment_id),
        parameter_value=float(parameter_value),
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


def _write_spectrum(output_dir: Path, experiment_id: str, rows: list[dict[str, Any]]) -> None:
    sweep._write_spectrum(output_dir, experiment_id, rows)


def test_01_geometry_and_exact_M0_contract() -> None:
    assert (sweep.MU, sweep.TAU, sweep.BETA_DEG) == (0.0, 0.0, 30.0)
    assert (sweep.L1, sweep.L2, sweep.L_TOTAL) == (1.0, 1.0, 2.0)
    assert (sweep.WIDTH, sweep.THICKNESS, sweep.PLY_THICKNESS) == (
        0.20,
        0.05,
        0.0125,
    )
    assert sweep.K == pytest.approx(5.0 / 6.0)
    assert sweep.base_material_contract() == {
        "delta": 0.1,
        "E1": 1.1,
        "E2": 0.9,
        "nu12": 0.3,
        "G12": 1.0 / 2.6,
        "G13": 1.0 / 2.6,
        "G23": 1.0 / 2.6,
        "rho0": 1.0,
    }
    assert sweep.M0 == pytest.approx(0.01)
    assert sweep.J0 == pytest.approx(1.0 * 0.20 * 0.05**3 / 12.0)


@pytest.mark.parametrize("eta", [0.0, 0.4, 0.8])
def test_02_H_L_density_scaling_changes_no_elastic_field(eta: float) -> None:
    materials = sweep.density_materials(eta)
    assert materials[sweep.HEAVY].rho == pytest.approx(1.0 + eta)
    assert materials[sweep.LIGHT].rho == pytest.approx(1.0 - eta)
    assert materials[sweep.HOMOGENEOUS].rho == 1.0
    elastic = {
        sweep._elastic_tuple(material) for material in materials.values()
    }
    assert elastic == {(1.1, 0.9, 0.3, 1.0 / 2.6, 1.0 / 2.6, 1.0 / 2.6)}


@pytest.mark.parametrize("xi_rho", [-0.8, -0.4, 0.0, 0.4, 0.8])
def test_03_signed_one_arm_density_scaling(xi_rho: float) -> None:
    outer, inner = sweep.signed_density_values(xi_rho)
    assert outer == pytest.approx(1.0 + xi_rho)
    assert inner == pytest.approx(1.0 - xi_rho)
    section = sweep.build_one_arm_layered_section(xi_rho)
    assert section.ply_densities == pytest.approx((outer, inner, inner, outer))
    assert sweep.elastic_properties_identical(section)


def test_04_densities_are_positive_on_both_complete_grids() -> None:
    assert len(sweep.eta_grid()) == 41
    assert len(sweep.xi_rho_grid()) == 81
    for eta in sweep.eta_grid():
        values = sweep.density_materials(float(eta))
        assert min(material.rho for material in values.values()) > 0.0
    for xi_rho in sweep.xi_rho_grid():
        assert min(sweep.signed_density_values(float(xi_rho))) > 0.0


def test_05_four_equal_zero_degree_plies_and_declared_layouts() -> None:
    sections = (
        sweep.build_eta_layout_section(sweep.OUTER_HEAVY_LAYOUT, 0.4),
        sweep.build_eta_layout_section(sweep.INNER_HEAVY_LAYOUT, 0.4),
        sweep.build_homogeneous_section(),
        sweep.build_one_arm_layered_section(-0.4),
    )
    assert sweep.OUTER_HEAVY_LAYOUT == ("H", "L", "L", "H")
    assert sweep.INNER_HEAVY_LAYOUT == ("L", "H", "H", "L")
    assert sweep.HOMOGENEOUS_LAYOUT == ("M0",) * 4
    for section in sections:
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


@pytest.mark.parametrize("eta", [0.0, 0.4, 0.8])
def test_06_elastic_matrices_reduced_stiffness_and_mass_are_invariant(
    eta: float,
) -> None:
    baseline = sweep.build_homogeneous_section()
    for layout in (sweep.OUTER_HEAVY_LAYOUT, sweep.INNER_HEAVY_LAYOUT):
        section = sweep.build_eta_layout_section(layout, eta)
        assert_allclose(section.laminate.A, baseline.laminate.A, rtol=1.0e-12)
        assert_allclose(section.laminate.D, baseline.laminate.D, rtol=1.0e-12)
        assert_allclose(section.laminate.shear, baseline.laminate.shear, rtol=1.0e-12)
        for name in ("A", "D", "S", "m"):
            assert getattr(section.properties, name) == pytest.approx(
                getattr(baseline.properties, name), rel=1.0e-11
            )
        assert section.laminate.I0 == pytest.approx(baseline.laminate.I0, rel=1.0e-12)
        assert sweep.rlb2e._scaled_B(section.laminate) <= 1.0e-12
        assert sweep.rlb2e._scaled_I1(section.laminate) <= 1.0e-12
        assert sweep.rlb2e._reduction_max_residual(section.properties) <= 1.0e-11


@pytest.mark.parametrize(
    ("eta", "outer_ratio", "inner_ratio"),
    [(0.0, 1.0, 1.0), (0.4, 1.3, 0.7), (0.8, 1.6, 0.4)],
)
def test_07_outer_inner_J_formulas_and_anti_phase_sum(
    eta: float, outer_ratio: float, inner_ratio: float
) -> None:
    baseline = sweep.build_homogeneous_section()
    outer = sweep.build_eta_layout_section(sweep.OUTER_HEAVY_LAYOUT, eta)
    inner = sweep.build_eta_layout_section(sweep.INNER_HEAVY_LAYOUT, eta)
    assert outer.properties.J / baseline.properties.J == pytest.approx(
        outer_ratio, rel=1.0e-11
    )
    assert inner.properties.J / baseline.properties.J == pytest.approx(
        inner_ratio, rel=1.0e-11
    )
    assert outer.properties.J + inner.properties.J == pytest.approx(
        2.0 * baseline.properties.J, rel=1.0e-11
    )


@pytest.mark.parametrize(
    ("xi_rho", "expected"),
    [(-0.8, 0.4), (-0.4, 0.7), (0.0, 1.0), (0.4, 1.3), (0.8, 1.6)],
)
def test_08_signed_layered_J_formula(xi_rho: float, expected: float) -> None:
    layered = sweep.build_one_arm_layered_section(xi_rho)
    baseline = sweep.build_homogeneous_section()
    assert layered.properties.J / baseline.properties.J == pytest.approx(
        expected, rel=1.0e-11
    )
    for name in ("A", "D", "S", "m"):
        assert getattr(layered.properties, name) == pytest.approx(
            getattr(baseline.properties, name), rel=1.0e-11
        )


def test_09_outer_and_inner_BeamProperties_differ_only_in_J() -> None:
    outer = sweep.build_eta_layout_section(sweep.OUTER_HEAVY_LAYOUT, 0.8)
    inner = sweep.build_eta_layout_section(sweep.INNER_HEAVY_LAYOUT, 0.8)
    assert sweep._properties_only_J_can_change(outer.properties, inner.properties)
    for name in ("A", "D", "S", "m", "K", "width"):
        assert getattr(outer.properties, name) == pytest.approx(
            getattr(inner.properties, name), rel=1.0e-11
        )
    assert outer.properties.J != inner.properties.J
    for name in (
        "axial_reduction",
        "bending_reduction",
        "shear_reduction_before_K",
    ):
        assert getattr(outer.properties, name).value == pytest.approx(
            getattr(inner.properties, name).value, rel=1.0e-11
        )
    beam, _coupled = sweep.rlb2e._physics_modules()
    omega = 2.3
    difference = beam.combined_state_matrix(omega, outer.properties) - (
        beam.combined_state_matrix(omega, inner.properties)
    )
    nonzero = np.argwhere(np.abs(difference) > 1.0e-14)
    assert nonzero.tolist() == [[5, 2]]
    assert difference[5, 2] == pytest.approx(
        -(outer.properties.J - inner.properties.J) * omega**2
    )


def test_10_constitutive_gate_and_408_section_rows() -> None:
    gate = sweep.constitutive_gate()
    assert gate["status"] == "PASS"
    assert gate["full_grid_density_positive"] is True
    assert all(item["status"] == "PASS" for item in gate["checks"])
    assert gate["maximum_residuals"]["J_formula_residual"] <= 1.0e-11
    assert gate["maximum_residuals"]["A_matrix_invariance_residual"] <= 1.0e-12
    assert gate["maximum_residuals"]["D_matrix_invariance_residual"] <= 1.0e-12
    assert gate["maximum_residuals"]["shear_matrix_invariance_residual"] <= 1.0e-12
    rows = sweep.section_property_rows()
    assert len(rows) == 408
    assert sum(row["experiment_id"] == sweep.EXPERIMENT_A for row in rows) == 246
    assert sum(row["experiment_id"] == sweep.EXPERIMENT_B for row in rows) == 162
    assert all(row["elastic_properties_identical_across_plies"] for row in rows)
    assert all(row["A_invariant"] for row in rows)
    assert all(row["D_invariant"] for row in rows)
    assert all(row["S_invariant"] for row in rows)
    assert all(row["m_invariant"] for row in rows)
    assert all(row["only_J_changes"] for row in rows)
    assert all(row["constitutive_status"] == "PASS" for row in rows)


def test_11_exact_grids_and_signed_continuation_paths() -> None:
    assert_allclose(sweep.eta_grid(), np.linspace(0.0, 0.8, 41), atol=1.0e-14)
    assert_allclose(sweep.xi_rho_grid(), np.linspace(-0.8, 0.8, 81), atol=1.0e-14)
    negative, positive = sweep.xi_rho_continuation_paths()
    assert_allclose(negative, np.linspace(0.0, -0.8, 41), atol=1.0e-14)
    assert_allclose(positive, np.linspace(0.0, 0.8, 41), atol=1.0e-14)
    assert set(np.round(negative, 10)) & set(np.round(positive, 10)) == {0.0}


def test_12_shared_zero_anchor_is_cloned_without_another_root_solve(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    all_A = _synthetic_A_rows()
    rows_A = [
        dict(row)
        for row in all_A
        if row["configuration_id"] == sweep.CONFIG_BOTH_OUTER
        and float(row["parameter_value"]) == 0.0
    ]
    rows_B: list[dict[str, Any]] = []

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("shared zero was recalculated")

    monkeypatch.setattr(sweep, "solve_point", forbidden)
    rows_A, rows_B, record = sweep.ensure_shared_zero_anchor(
        tmp_path, rows_A, rows_B, []
    )
    groups = [
        sweep._roots_for(rows_A, configuration, 0.0)
        for configuration in sweep.CONFIGURATIONS_A
    ]
    groups.append(sweep._roots_for(rows_B, sweep.CONFIG_ONE_ARM, 0.0))
    for values in groups[1:]:
        assert_allclose(values, groups[0], rtol=0.0, atol=0.0)
    zero_rows = [
        row for row in [*rows_A, *rows_B] if float(row["parameter_value"]) == 0.0
    ]
    assert len(zero_rows) == 4 * 9
    assert all(row["shared_zero_contrast_anchor_reused"] for row in zero_rows)
    assert {row["physical_solve_id"] for row in zero_rows} == {
        "RLB-2G_SHARED_ZERO_PHYSICAL_SOLVE"
    }
    assert record["logical_group_count"] == 4
    assert record["physical_solve_count"] == 0
    assert record["reused_existing"] is True


def test_13_synthetic_base_inventories_have_all_groups_roots_and_order() -> None:
    rows_A = _synthetic_A_rows()
    rows_B = _synthetic_B_rows()
    audit_A = sweep.audit_spectrum_rows(rows_A, sweep.EXPERIMENT_A)
    audit_B = sweep.audit_spectrum_rows(rows_B, sweep.EXPERIMENT_B)
    assert audit_A["status"] == "PASS"
    assert audit_A["base_group_count"] == 123
    assert audit_A["base_row_count"] == 1107
    assert audit_B["status"] == "PASS"
    assert audit_B["base_group_count"] == 81
    assert audit_B["base_row_count"] == 729
    for groups in (sweep._base_group_index(rows_A), sweep._base_group_index(rows_B)):
        for group in groups.values():
            ordered = sorted(group, key=lambda row: int(row["sorted_position"]))
            assert [int(row["sorted_position"]) for row in ordered] == list(range(1, 10))
            assert np.all(np.diff([float(row["Omega"]) for row in ordered]) > 0.0)
    assert sweep.audit_spectrum_rows(rows_A[:-1], sweep.EXPERIMENT_A)["status"] == "FAIL"


def test_14_root9_is_guard_not_plot_data_and_no_tail_is_allowed() -> None:
    for rows in (_synthetic_A_rows(), _synthetic_B_rows()):
        guards = [row for row in rows if int(row["sorted_position"]) == 9]
        assert all(row["root_role"] == "ROOT_9_GUARD" for row in guards)
        assert all(row["guard_flag"] for row in guards)
        plotted = sweep.canonical_plot_rows(rows)
        assert max(int(row["sorted_position"]) for row in plotted) == 8
        assert max(int(row["sorted_position"]) for row in rows) == 9
        assert all(not row["roots_above_9_computed"] for row in rows)
        assert all(not row["predictor_used_as_final"] for row in rows)
    assert len(sweep.canonical_plot_rows(_synthetic_A_rows())) == 984
    assert len(sweep.canonical_plot_rows(_synthetic_B_rows())) == 648


def test_15_predictor_is_locator_only_and_frozen_search_thresholds_are_reused() -> None:
    policy = sweep.rlb2e._rlb2e_search_policy()
    predecessor = sweep.rlb2e._root_tools().SearchPolicy()
    assert (policy.requested_roots, policy.guard_roots, policy.required_slots) == (8, 1, 9)
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
    for position in range(1, 10):
        Omega = 10.0 * position
        diagnostic = SimpleNamespace(
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
            diagnostics=diagnostic,
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
        sweep.EXPERIMENT_B,
        sweep.CONFIG_ONE_ARM,
        0.2,
        slots,
        solve_id="synthetic",
        physical_solve_id="synthetic",
        transaction_id="synthetic",
        solve_mode="FAST_LOCAL",
        fallback_used=False,
        predicted=predicted,
        search_right=95.0,
        unresolved=0,
        continuation_leg="XI_RHO_POSITIVE",
    )
    assert [float(row["Omega"]) for row in rows] == [10.0 * value for value in range(1, 10)]
    assert [float(row["predictor_Omega"]) for row in rows] == pytest.approx(predicted)
    assert all(row["predictor_used_as_final"] is False for row in rows)


def test_16_arm_swap_orchestration_uses_exactly_six_declared_checks(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows_A, rows_B = _synthetic_A_rows(), _synthetic_B_rows()
    calls: list[tuple[str, str, float, bool]] = []

    def fake_solve(
        experiment_id: str, configuration_id: str, value: float, **kwargs: Any
    ) -> sweep.PointSolution:
        calls.append((experiment_id, configuration_id, float(value), kwargs["swapped"]))
        source = rows_A if experiment_id == sweep.EXPERIMENT_A else rows_B
        return _synthetic_solution(
            source, experiment_id, configuration_id, value, swapped=True
        )

    monkeypatch.setattr(sweep, "solve_point", fake_solve)
    result = sweep.arm_swap_checks(rows_A, rows_B)
    assert result["status"] == "PASS"
    assert len(result["checks"]) == 6
    assert [call[2] for call in calls] == [0.4, 0.8, -0.8, -0.4, 0.4, 0.8]
    assert all(call[3] is True for call in calls)
    assert all(item["root_count"] == 9 for item in result["checks"])


def test_17_neighbour_audit_keeps_experiments_and_configurations_separate() -> None:
    audit = sweep.neighbour_audit_rows(_synthetic_A_rows(), _synthetic_B_rows())
    assert len(audit) == 3 * 39 * 8 + 79 * 8
    assert {row["experiment_id"] for row in audit} == {
        sweep.EXPERIMENT_A,
        sweep.EXPERIMENT_B,
    }
    assert all(row["smoothing_applied"] is False for row in audit)
    manually_flagged = [
        {
            "experiment_id": sweep.EXPERIMENT_A,
            "configuration_id": sweep.CONFIG_BOTH_OUTER,
            "parameter_value": 0.4,
            "flagged": True,
        },
        {
            "experiment_id": sweep.EXPERIMENT_B,
            "configuration_id": sweep.CONFIG_ONE_ARM,
            "parameter_value": 0.4,
            "flagged": False,
        },
    ]
    assert sweep.flagged_repair_points(manually_flagged) == [
        (sweep.EXPERIMENT_A, sweep.CONFIG_BOTH_OUTER, 0.4)
    ]


def test_18_reproduced_local_feature_is_retained_without_cross_contamination(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows_A, rows_B = _synthetic_A_rows(), _synthetic_B_rows()
    before_B = json.dumps(rows_B, sort_keys=True)
    audit = [
        {
            "experiment_id": sweep.EXPERIMENT_A,
            "configuration_id": sweep.CONFIG_BOTH_OUTER,
            "parameter_value": 0.4,
            "sorted_position": 4,
            "flagged": True,
        }
    ]

    def fake_solve(experiment_id: str, configuration_id: str, value: float, **_kwargs: Any) -> sweep.PointSolution:
        return _synthetic_solution(rows_A, experiment_id, configuration_id, value)

    monkeypatch.setattr(sweep, "solve_point", fake_solve)
    rows_A, rows_B, audit, records = sweep.apply_local_repairs(rows_A, rows_B, audit)
    assert records[0]["status"] == "REPRODUCED_AFTER_LOCAL_REPAIR"
    assert records[0]["smoothing_applied"] is False
    assert audit[0]["repair_status"] == "REPRODUCED_AFTER_LOCAL_REPAIR"
    assert json.dumps(rows_B, sort_keys=True) == before_B
    assert sweep.audit_spectrum_rows(rows_A, sweep.EXPERIMENT_A)["status"] == "PASS"


def test_19_unresolved_repair_creates_visible_gap_only_at_affected_position(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    rows_A, rows_B = _synthetic_A_rows(), _synthetic_B_rows()
    audit = [
        {
            "experiment_id": sweep.EXPERIMENT_B,
            "configuration_id": sweep.CONFIG_ONE_ARM,
            "parameter_value": -0.4,
            "sorted_position": 3,
            "flagged": True,
        }
    ]

    def fail(*_args: Any, **_kwargs: Any) -> sweep.PointSolution:
        raise RuntimeError("synthetic unresolved neighbourhood")

    monkeypatch.setattr(sweep, "solve_point", fail)
    rows_A, rows_B, audit, records = sweep.apply_local_repairs(rows_A, rows_B, audit)
    group = sweep._canonical_group(rows_B, sweep.CONFIG_ONE_ARM, -0.4)
    assert records[0]["status"] == "UNRESOLVED"
    assert audit[0]["repair_status"] == "UNRESOLVED"
    assert math.isnan(float(group[2]["Lambda"]))
    assert group[2]["point_status"] == "UNRESOLVED_AFTER_LOCAL_REPAIR"
    assert sweep.audit_plot_rows(rows_B, sweep.EXPERIMENT_B)["status"] == "PASS"


def test_20_atomic_completed_group_is_preserved(tmp_path: Path) -> None:
    all_rows = _synthetic_A_rows()
    rows = all_rows[:9]
    _write_spectrum(tmp_path, sweep.EXPERIMENT_A, rows)
    path = tmp_path / sweep.SPECTRUM_A_FILENAME
    before = path.read_bytes()
    solution = _synthetic_solution(
        all_rows,
        sweep.EXPERIMENT_A,
        sweep.CONFIG_BOTH_OUTER,
        0.0,
    )
    preserved = sweep._write_point_transaction(tmp_path, rows, solution)
    assert preserved == rows
    assert path.read_bytes() == before


def test_fail_closed_durable_benchmark_without_measurements_is_not_rewritten(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    all_rows = _synthetic_A_rows()
    durable_group = [
        row
        for row in all_rows
        if row["configuration_id"] == sweep.CONFIG_BOTH_OUTER
        and float(row["parameter_value"]) == pytest.approx(0.8)
    ]
    assert len(durable_group) == sweep.K_GUARD
    _write_spectrum(tmp_path, sweep.EXPERIMENT_A, durable_group)
    spectrum_path = tmp_path / sweep.SPECTRUM_A_FILENAME
    before = spectrum_path.read_bytes()

    recovered = sweep._solution_record(
        _synthetic_solution(
            all_rows,
            sweep.EXPERIMENT_A,
            sweep.CONFIG_BOTH_OUTER,
            0.8,
        ),
        benchmark=True,
    )
    recovered["wall_time_seconds"] = 0.0
    recovered["determinant_evaluations"] = 0

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("durable benchmark was recomputed")

    monkeypatch.setattr(sweep, "solve_point", forbidden)
    with pytest.raises(RuntimeError, match="without measured timing"):
        sweep._ensure_benchmark_point(
            tmp_path,
            durable_group,
            [recovered],
            experiment_id=sweep.EXPERIMENT_A,
            configuration_id=sweep.CONFIG_BOTH_OUTER,
            parameter_value=0.8,
        )
    assert spectrum_path.read_bytes() == before
    assert {path.name for path in tmp_path.iterdir()} == {
        sweep.SPECTRUM_A_FILENAME
    }


@pytest.mark.parametrize(
    "damage",
    ("missing_root9", "bad_root9", "missing_section_key"),
)
def test_plot_only_rejects_incomplete_inventory_before_matplotlib(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    damage: str,
) -> None:
    rows_A = _synthetic_A_rows()
    rows_B = _synthetic_B_rows()
    sections = sweep.section_property_rows()
    if damage == "missing_root9":
        rows_A = [
            row
            for row in rows_A
            if not (
                row["configuration_id"] == sweep.CONFIG_BOTH_OUTER
                and float(row["parameter_value"]) == pytest.approx(0.8)
                and int(row["sorted_position"]) == sweep.K_GUARD
            )
        ]
    elif damage == "bad_root9":
        target = next(
            row
            for row in rows_B
            if float(row["parameter_value"]) == pytest.approx(0.8)
            and int(row["sorted_position"]) == sweep.K_GUARD
        )
        target["scaled_sigma_ratio"] = 10.0 * sweep.ROOT_SINGULAR_RATIO_TOLERANCE
    else:
        sections.pop()

    _write_spectrum(tmp_path, sweep.EXPERIMENT_A, rows_A)
    _write_spectrum(tmp_path, sweep.EXPERIMENT_B, rows_B)
    sweep.rlb2e._atomic_write_csv(
        tmp_path / sweep.SECTION_FILENAME,
        sections,
        sweep.SECTION_FIELDS,
    )

    imported_matplotlib: list[str] = []
    original_import = builtins.__import__

    def guarded_import(
        name: str,
        globals: Any = None,
        locals: Any = None,
        fromlist: Any = (),
        level: int = 0,
    ) -> Any:
        if name.startswith("matplotlib"):
            imported_matplotlib.append(name)
            raise AssertionError("matplotlib imported before inventory audit")
        return original_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarded_import)
    with pytest.raises(RuntimeError, match="Plot-only requires complete"):
        sweep.create_plots_from_csv(tmp_path)
    assert imported_matplotlib == []
    assert not any(
        (tmp_path / filename).exists()
        for filename in (
            sweep.PLOT_A_FILENAME,
            sweep.PLOT_B_FILENAME,
            sweep.PLOT_J_FILENAME,
        )
    )


def test_21_plot_only_creates_three_figures_and_never_calls_physics(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write_spectrum(tmp_path, sweep.EXPERIMENT_A, _synthetic_A_rows())
    _write_spectrum(tmp_path, sweep.EXPERIMENT_B, _synthetic_B_rows())
    sweep.rlb2e._atomic_write_csv(
        tmp_path / sweep.SECTION_FILENAME,
        sweep.section_property_rows(),
        sweep.SECTION_FIELDS,
    )

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("plot_only called physics or a root solver")

    monkeypatch.setattr(sweep, "ROOT_CALCULATION_COUNT", 0)
    monkeypatch.setattr(sweep, "solve_point", forbidden)
    monkeypatch.setattr(sweep, "make_matrix_provider", forbidden)
    monkeypatch.setattr(sweep, "_build_section", forbidden)
    monkeypatch.setattr(sweep.rlb2e, "_physics_modules", forbidden)
    result = sweep.create_plots_from_csv(tmp_path)
    assert result["root_calculation_count"] == 0
    assert result["matrix_assembly_count"] == 0
    assert result["determinant_evaluations"] == 0
    assert result["SVD_evaluations"] == 0
    assert result["experiment_A_panels"] == 3
    assert result["experiment_B_panels"] == 1
    assert result["frequency_line_count_A"] == 24
    assert result["frequency_line_count_B"] == 8
    assert result["property_line_count"] == 3
    assert result["root9_plotted"] is False
    for filename in (sweep.PLOT_A_FILENAME, sweep.PLOT_B_FILENAME, sweep.PLOT_J_FILENAME):
        assert (tmp_path / filename).is_file()


def test_22_completed_missing_only_is_bytewise_idempotent(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write_spectrum(tmp_path, sweep.EXPERIMENT_A, _synthetic_A_rows())
    _write_spectrum(tmp_path, sweep.EXPERIMENT_B, _synthetic_B_rows())
    for name in sweep.MANDATORY_COMPLETED_OUTPUTS:
        path = tmp_path / name
        if path.exists():
            continue
        path.write_bytes(b"synthetic completed output\n")
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
        "reused_RLB2E_script_sha256": sweep.rlb2e._sha256(ROOT / sweep.RLB2E_SCRIPT_PATH),
        "reused_RLB2F_script_sha256": sweep.rlb2e._sha256(ROOT / sweep.RLB2F_SCRIPT_PATH),
        "output_hashes": output_hashes,
    }
    sweep.rlb2e._atomic_write_json(tmp_path / sweep.MANIFEST_FILENAME, manifest)
    before = {path.name: path.read_bytes() for path in tmp_path.iterdir()}

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("completed missing-only called root or constitutive path")

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


def test_hash_incompatible_complete_tree_fails_without_mutating_any_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _write_spectrum(tmp_path, sweep.EXPERIMENT_A, _synthetic_A_rows())
    _write_spectrum(tmp_path, sweep.EXPERIMENT_B, _synthetic_B_rows())
    for name in sweep.MANDATORY_COMPLETED_OUTPUTS:
        path = tmp_path / name
        if not path.exists():
            path.write_bytes(f"synthetic frozen {name}\n".encode("utf-8"))
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
        "reused_RLB2F_script_sha256": sweep.rlb2e._sha256(
            ROOT / sweep.RLB2F_SCRIPT_PATH
        ),
        "output_hashes": output_hashes,
    }
    sweep.rlb2e._atomic_write_json(tmp_path / sweep.MANIFEST_FILENAME, manifest)

    # Make the tree hash-incompatible only after recording a materially complete
    # manifest.  The fail-closed path must preserve this evidence byte for byte.
    report_path = tmp_path / sweep.REPORT_FILENAME
    report_path.write_bytes(report_path.read_bytes() + b"hash mismatch\n")
    before = {
        path.name: sweep.rlb2e._sha256(path)
        for path in tmp_path.iterdir()
        if path.is_file()
    }

    def forbidden(*_args: Any, **_kwargs: Any) -> Any:
        raise AssertionError("hash-incompatible tree entered a calculation path")

    monkeypatch.setattr(sweep, "solve_point", forbidden)
    monkeypatch.setattr(sweep, "constitutive_gate", forbidden)
    monkeypatch.setattr(sweep, "section_property_rows", forbidden)
    with pytest.raises(RuntimeError, match="hash-incompatible"):
        sweep.run_workflow(tmp_path, missing_only=True)
    after = {
        path.name: sweep.rlb2e._sha256(path)
        for path in tmp_path.iterdir()
        if path.is_file()
    }
    assert after == before


def test_23_manifest_records_both_complete_frequency_map_instances() -> None:
    contract = sweep.contract_payload()
    policies = contract["frequency_map_policy_instances"]
    assert policies["experiment_A"] == sweep.POLICY_EXPERIMENT_A
    assert policies["experiment_B"] == sweep.POLICY_EXPERIMENT_B
    common = {
        "frequency_map_policy": "frequency-map-v1",
        "calculation_mode": "fast_plot",
        "spectrum_semantics": "sorted_positions",
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
    assert policies["experiment_A"]["sweep_parameter"] == "eta"
    assert policies["experiment_A"]["parameter_grid"] == "0.00:0.02:0.80"
    assert policies["experiment_B"]["sweep_parameter"] == "xi_rho"
    assert policies["experiment_B"]["parameter_grid"] == "-0.80:0.02:0.80"
    assert policies["experiment_B"]["continuation_anchor"] == "0.00"
    assert contract["shared_zero_contrast_anchor"] == {
        "logical_group_count": 4,
        "physical_solve_count": 1,
        "clone_count": 3,
    }
    assert contract["normalization"]["density_or_J_dependent"] is False


def test_24_benchmark_accounting_uses_five_physical_anchors() -> None:
    source = ast.parse(SCRIPT_PATH.read_text(encoding="utf-8"))
    text = SCRIPT_PATH.read_text(encoding="utf-8")
    assert "total_logical_groups\": 204" in text
    assert "total_unique_physical_base_points\": 201" in text
    assert "physical_anchor_count\": 5" in text
    assert "remaining_unique_physical_points\": remaining" in text
    assert any(isinstance(node, ast.FunctionDef) and node.name == "run_benchmarks" for node in ast.walk(source))


def test_25_no_forbidden_imports_local_solver_smoothing_or_interpolation() -> None:
    source = SCRIPT_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    imported: list[str] = []
    functions: set[str] = set()
    calls: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.append(node.module)
            imported.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.FunctionDef):
            functions.add(node.name)
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name):
                calls.add(node.func.id)
            elif isinstance(node.func, ast.Attribute):
                calls.add(node.func.attr)
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
    assert not any("characteristic_determinant" in name for name in functions)
    assert "determinant" not in functions
    assert not ({"interp", "interp1d", "savgol_filter", "UnivariateSpline"} & calls)
    assert '"smoothing": False' in source
    assert '"interpolation_based_frequencies": False' in source


def test_26_predecessor_output_hashes_are_preserved_if_available() -> None:
    checked = 0
    for result_dir in (sweep.RLB2E_RESULT_DIR, sweep.RLB2F_RESULT_DIR):
        manifest_path = result_dir / "run_manifest.json"
        if not manifest_path.is_file():
            continue
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for name, expected in manifest.get("output_hashes", {}).items():
            assert sweep.rlb2e._sha256(result_dir / name) == expected
        checked += 1
    if checked == 0:
        pytest.skip("Neither predecessor result tree is available.")


def test_27_generated_outputs_if_present() -> None:
    if not RESULT_DIR.is_dir():
        pytest.skip("RLB-2G production outputs have not been generated yet.")
    required = set(sweep.MANDATORY_COMPLETED_OUTPUTS) | {sweep.MANIFEST_FILENAME}
    assert required <= {path.name for path in RESULT_DIR.iterdir() if path.is_file()}
    manifest = json.loads((RESULT_DIR / sweep.MANIFEST_FILENAME).read_text(encoding="utf-8"))
    rows_A = sweep.rlb2e._read_csv(RESULT_DIR / sweep.SPECTRUM_A_FILENAME)
    rows_B = sweep.rlb2e._read_csv(RESULT_DIR / sweep.SPECTRUM_B_FILENAME)
    sections = sweep.rlb2e._read_csv(RESULT_DIR / sweep.SECTION_FILENAME)
    arm_swap = json.loads((RESULT_DIR / sweep.ARM_SWAP_FILENAME).read_text(encoding="utf-8"))
    assert sweep.audit_spectrum_rows(rows_A, sweep.EXPERIMENT_A)["status"] == "PASS"
    assert sweep.audit_spectrum_rows(rows_B, sweep.EXPERIMENT_B)["status"] == "PASS"
    assert len([row for row in rows_A if row["grid_kind"] == "BASE"]) == 1107
    assert len([row for row in rows_B if row["grid_kind"] == "BASE"]) == 729
    assert len(sections) == 408
    assert manifest["constitutive_gate"]["status"] == "PASS"
    assert manifest["counts"]["experiment_A_groups"] == 123
    assert manifest["counts"]["experiment_B_groups"] == 81
    assert manifest["counts"]["shared_zero_logical_groups"] == 4
    assert manifest["counts"]["shared_zero_physical_solves"] == 1
    assert manifest["root_contract"]["root9_plotted"] is False
    assert manifest["root_contract"]["roots_above_9_computed"] is False
    assert manifest["baseline_consistency"]["status"] == "PASS"
    assert arm_swap["status"] == "PASS"
    assert len(arm_swap["checks"]) == 6
    zero_rows = [
        row for row in [*rows_A, *rows_B] if float(row["parameter_value"]) == 0.0
    ]
    assert len(zero_rows) == 36
    assert all(sweep._as_bool(row["shared_zero_contrast_anchor_reused"]) for row in zero_rows)
    assert len({row["physical_solve_id"] for row in zero_rows}) == 1
    assert manifest["contract"]["frequency_map_policy_instances"] == {
        "experiment_A": sweep.POLICY_EXPERIMENT_A,
        "experiment_B": sweep.POLICY_EXPERIMENT_B,
    }
    assert manifest["plots"]["root_calculation_count"] == 0
    assert manifest["plots"]["root9_plotted"] is False
    for path, expected in manifest["production_physics_hashes"].items():
        assert sweep.rlb2e._sha256(ROOT / path) == expected
    for name, expected in manifest["output_hashes"].items():
        assert sweep.rlb2e._sha256(RESULT_DIR / name) == expected
    for filename in (sweep.PLOT_A_FILENAME, sweep.PLOT_B_FILENAME, sweep.PLOT_J_FILENAME):
        assert (RESULT_DIR / filename).stat().st_size > 0
