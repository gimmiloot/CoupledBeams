from __future__ import annotations

import ast
import inspect
import json
import subprocess
import tempfile
import unittest
from dataclasses import fields
from pathlib import Path
from unittest import mock

import matplotlib.pyplot as plt
import numpy as np

from scripts.analysis.anisotropic_rods import (
    plot_yartsev_ch2_supervisor_figures as supervisor,
)
from scripts.lib.yartsev_ch2_monoclinic_rod import BookMaterial, hms_dx_209_material


class YartsevChapter2SupervisorFiguresTest(unittest.TestCase):
    @staticmethod
    def _quality() -> dict[str, object]:
        return {
            "quality_status": "PASS",
            "quality_basis": "scaled",
            "root_status": "accepted_brent",
            "accepted_determinant_residual": 1.0e-12,
            "accepted_relative_singular_residual": 2.0e-12,
            "scaled_determinant_residual": 1.0e-12,
            "scaled_relative_singular_residual": 2.0e-12,
            "physical_raw_normalized_determinant_residual": 3.0e-12,
            "physical_raw_relative_singular_residual": 4.0e-12,
        }

    def _synthetic_fast_result(
        self, beta_values: tuple[float, ...], frequency_offset: float = 0.0
    ) -> supervisor.FastSweepResult[supervisor.SpectrumResult]:
        spectra: dict[float, supervisor.SpectrumResult] = {}
        for beta in beta_values:
            roots = [
                supervisor.RootResult(
                    omega=complex(2.0 * np.pi * (100.0 * mode + beta + frequency_offset)),
                    frequency_hz=100.0 * mode + beta + frequency_offset,
                    determinant_residual=1.0e-12,
                    raw_determinant_abs=1.0e-12,
                    sigma_min=1.0e-12,
                    sigma_max=1.0,
                    relative_singular_residual=1.0e-12,
                    min_neighbor_distance_hz=100.0,
                    refinements=1,
                    status="accepted_brent",
                )
                for mode in range(1, supervisor.GUARD_ROOT_COUNT + 1)
            ]
            spectra[beta] = supervisor.SpectrumResult(
                roots=roots,
                quality=[self._quality() for _ in roots],
                scaled_matrix_evaluations=0,
                raw_quality_matrix_evaluations=len(roots),
                runtime_seconds=0.0,
            )
        return supervisor.FastSweepResult(
            spectra=spectra,
            data_origins={beta: "fast_local" for beta in beta_values},
            fallback_reasons={},
            anchor_relative_errors={beta_values[0]: 0.0},
            counters=supervisor.PerformanceCounters(),
            runtime_seconds=0.0,
        )

    def _synthetic_comparison_rows(self, figure: int) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        for beta in (0.0, 90.0):
            for mode in range(1, supervisor.GUARD_ROOT_COUNT + 1):
                common = {
                    "beta_deg": beta,
                    "mode": mode,
                    "root_role": (
                        "plotted"
                        if mode <= supervisor.PLOTTED_ROOT_COUNT
                        else "guard"
                    ),
                }
                quality_columns = {
                    f"{prefix}_{key}": value
                    for prefix in ("book_slope", "section_clamp", "timoshenko", "eb")
                    for key, value in self._quality().items()
                }
                if figure == 2:
                    rows.append(
                        {
                            **common,
                            **quality_columns,
                            "book_slope_frequency_hz": 100.0 * mode + beta,
                            "book_slope_lambda": 0.8 * mode + beta / 900.0,
                            "book_slope_quality_status": "PASS",
                            "book_slope_quality_basis": "scaled",
                            "book_slope_root_status": "accepted_brent",
                            "book_slope_accepted_determinant_residual": 1.0e-12,
                            "book_slope_accepted_relative_singular_residual": 2.0e-12,
                            "section_clamp_frequency_hz": 101.0 * mode + beta,
                            "section_clamp_lambda": 0.81 * mode + beta / 900.0,
                            "section_clamp_quality_status": "PASS",
                            "section_clamp_quality_basis": "scaled",
                            "section_clamp_root_status": "accepted_brent",
                            "section_clamp_accepted_determinant_residual": 1.0e-12,
                            "section_clamp_accepted_relative_singular_residual": 2.0e-12,
                            "relative_clamp_difference": 0.0125,
                        }
                    )
                else:
                    rows.append(
                        {
                            **common,
                            **quality_columns,
                            "figure": figure,
                            "clamp_type": (
                                "book_slope_clamp"
                                if figure == 3
                                else "timoshenko_section_clamp"
                            ),
                            "timoshenko_frequency_hz": 100.0 * mode + beta,
                            "timoshenko_lambda": 0.8 * mode + beta / 900.0,
                            "timoshenko_quality_status": "PASS",
                            "timoshenko_quality_basis": "scaled",
                            "timoshenko_root_status": "accepted_brent",
                            "timoshenko_accepted_determinant_residual": 1.0e-12,
                            "timoshenko_accepted_relative_singular_residual": 2.0e-12,
                            "eb_frequency_hz": 102.0 * mode + beta,
                            "eb_lambda": 0.82 * mode + beta / 900.0,
                            "eb_quality_status": "PASS",
                            "eb_quality_basis": "scaled",
                            "eb_root_status": "accepted_brent",
                            "eb_accepted_determinant_residual": 1.0e-12,
                            "eb_accepted_relative_singular_residual": 2.0e-12,
                            "relative_theory_difference": 0.025,
                        }
                    )
        return rows

    def test_import_has_no_computation_or_plotting_entry_side_effect(self) -> None:
        self.assertTrue(callable(supervisor.main))
        self.assertFalse(hasattr(supervisor, "_IMPORT_TIME_RESULTS"))

    def test_lambda_equivalent_formulas(self) -> None:
        omega = np.array([1.0, 2.0 * np.pi * 250.0, 2.0 * np.pi * 1800.0])
        kwargs = {
            "rho_kg_m3": 1580.0,
            "a_m": 0.005,
            "b_m": 0.020,
            "length_1_m": 0.400,
            "length_2_m": 0.400,
            "elastic_ex_pa": 191.0e9,
        }
        direct = supervisor.lambda_from_omega(omega, **kwargs)
        equivalent = supervisor.lambda_from_omega_equivalent(omega, **kwargs)
        np.testing.assert_allclose(direct, equivalent, rtol=2.0e-15, atol=0.0)

    def test_figure_2_canonical_book_material_and_geometry(self) -> None:
        preset = supervisor.FIGURE_2_PRESET
        self.assertIs(preset.material_factory, BookMaterial)
        self.assertEqual(preset.a_m, 9.76e-3)
        self.assertEqual(preset.b_m, 2.524e-2)
        self.assertEqual(preset.length_1_m, 0.295)
        self.assertEqual(preset.length_2_m, 0.295)
        self.assertEqual(preset.theta_1_deg, 0.0)
        self.assertEqual(preset.theta_2_deg, 0.0)
        self.assertEqual(preset.mu, 0.0)

    def test_figures_3_4_hms_rectangular_parameters(self) -> None:
        preset = supervisor.FIGURES_3_4_PRESET
        self.assertIs(preset.material_factory, hms_dx_209_material)
        self.assertEqual(preset.material_name, "HMS/DX-209")
        self.assertEqual(preset.a_m, 0.005)
        self.assertEqual(preset.b_m, 0.020)
        self.assertEqual(preset.b_m / preset.a_m, 4.0)
        self.assertEqual(preset.length_1_m, 0.400)
        self.assertEqual(preset.length_2_m, 0.400)
        self.assertEqual(preset.theta_1_deg, 0.0)
        self.assertEqual(preset.theta_2_deg, 0.0)
        self.assertEqual(preset.mu, 0.0)

    def test_rectangular_area_and_second_moment_are_used(self) -> None:
        area, inertia = supervisor.rectangular_reference_section(0.005, 0.020)
        self.assertEqual(area, 0.005 * 0.020)
        self.assertEqual(inertia, 0.005**3 * 0.020 / 12.0)

    def test_presets_have_no_radius_or_old_thickness_parameter(self) -> None:
        field_names = {field.name for field in fields(supervisor.FigurePreset)}
        self.assertNotIn("r", field_names)
        self.assertNotIn("eta", field_names)
        parser_destinations = {
            action.dest for action in supervisor.parse_args.__globals__["argparse"].ArgumentParser()._actions
        }
        self.assertNotIn("eta", parser_destinations)

    def test_no_isotropic_steel_defaults(self) -> None:
        for preset in supervisor.FIGURE_PRESETS:
            material = preset.material_factory()
            self.assertNotEqual(material.E1_real, 2.1e11)
            self.assertNotEqual(material.rho, 7800.0)

    def test_beta_grid_is_inclusive_and_default_step_is_half_degree(self) -> None:
        self.assertEqual(supervisor.DEFAULT_BETA_STEP_DEG, 0.5)
        grid = supervisor.beta_grid()
        self.assertEqual(grid[0], 0.0)
        self.assertEqual(grid[-1], 90.0)
        self.assertEqual(len(grid), 181)
        np.testing.assert_array_equal(np.diff(grid), np.full(180, 0.5))

    def test_root_contract(self) -> None:
        self.assertEqual(supervisor.GUARD_ROOT_COUNT, 7)
        self.assertEqual(supervisor.PLOTTED_ROOT_COUNT, 6)

    def test_scientific_imports_are_only_chapter_2_anisotropic_modules(self) -> None:
        source_path = Path(inspect.getsourcefile(supervisor) or "")
        tree = ast.parse(source_path.read_text(encoding="utf-8"))
        script_lib_imports = {
            node.module
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
            and node.module is not None
            and node.module.startswith("scripts.lib.")
        }
        scientific_builders = {
            "scripts.lib.yartsev_ch2_monoclinic_rod",
            "scripts.lib.yartsev_ch2_coupled_rods",
            "scripts.lib.yartsev_ch2_rectangular_eb",
        }
        self.assertEqual(
            script_lib_imports,
            scientific_builders | {"scripts.lib.yartsev_ch2_fast_beta_sweep"},
        )
        forbidden = {
            "scripts.lib.formulas",
            "src.my_project.analytic.FreqFromAngle",
            "src.my_project.analytic.FreqFromMu",
            "src.my_project.analytic.FreqMuNet",
            "src.my_project.analytic.formulas",
            "src.my_project.analytic.solvers",
            "src.my_project.fem.python_fem",
        }
        self.assertTrue(script_lib_imports.isdisjoint(forbidden))

    def test_section_builder_calls_existing_transfer_clamp_and_unchanged_joint(self) -> None:
        original_joint = supervisor.joint_matrix_book
        original_transfer = supervisor.physical_state_transfer_matrix
        calls: list[tuple[str, object]] = []
        clamp = np.zeros((6, 3), dtype=np.complex128)
        clamp[3:, :] = np.eye(3)
        joint = np.arange(72, dtype=float).reshape(6, 12)

        with mock.patch.object(
            supervisor,
            "physical_state_transfer_matrix",
            side_effect=lambda omega, point: (
                calls.append(("transfer", point)) or np.eye(6, dtype=np.complex128)
            ),
        ), mock.patch.object(
            supervisor,
            "cantilever_clamp_matrix",
            side_effect=lambda point, variant, scaled=False: (
                calls.append((variant, scaled)) or clamp
            ),
        ), mock.patch.object(
            supervisor,
            "joint_matrix_book",
            side_effect=lambda beta: calls.append(("joint", beta)) or joint,
        ):
            result = supervisor._section_clamp_coupled_boundary_matrix_raw(
                10.0, 0.3, object(), object()  # type: ignore[arg-type]
            )

        self.assertEqual(result.shape, (6, 6))
        self.assertEqual(
            [call for call in calls if call[0] == "timoshenko_section_clamp"],
            [("timoshenko_section_clamp", False)] * 2,
        )
        self.assertEqual(len([call for call in calls if call[0] == "transfer"]), 2)
        self.assertEqual(len([call for call in calls if call[0] == "joint"]), 1)
        self.assertIs(supervisor.joint_matrix_book, original_joint)
        self.assertIs(supervisor.physical_state_transfer_matrix, original_transfer)

    def test_section_clamp_straight_reference_smoke(self) -> None:
        result = supervisor.section_clamp_straight_reference(
            supervisor.FIGURES_3_4_PRESET
        )
        self.assertEqual(result.status, "PASS")
        self.assertLessEqual(result.maximum_relative_frequency_difference, 1.0e-8)
        self.assertEqual(len(result.coupled.roots), 7)
        self.assertEqual(len(result.straight.roots), 7)
        coupled = np.asarray(
            [root.frequency_hz for root in result.coupled.roots], dtype=float
        )
        self.assertAlmostEqual(coupled[4], 762.0131191222644, places=6)
        self.assertAlmostEqual(coupled[5], 766.8613179011050, places=6)

    def test_close_pair_is_complete_on_small_book_material_beta_slice(self) -> None:
        sweep = supervisor._sweep(
            supervisor.FIGURE_2_PRESET,
            np.asarray([32.5, 33.0, 33.5, 34.0]),
            "book_slope_clamp",
        )
        spectrum = sweep.spectra[33.5]
        roots = spectrum.roots
        frequencies = np.asarray([root.frequency_hz for root in roots], dtype=float)
        self.assertEqual(len(frequencies), supervisor.GUARD_ROOT_COUNT)
        self.assertGreater(frequencies[4] - frequencies[3], 0.0)
        self.assertLess(frequencies[4] - frequencies[3], 0.1)
        self.assertTrue(
            all(item["quality_status"] == "PASS" for item in spectrum.quality)
        )

    def test_eb_data_for_figures_3_4_are_exactly_identical(self) -> None:
        rows_3 = self._synthetic_comparison_rows(3)
        rows_4 = self._synthetic_comparison_rows(4)
        self.assertTrue(supervisor.eb_arrays_are_identical(rows_3, rows_4))
        rows_4[0]["eb_lambda"] = float(rows_4[0]["eb_lambda"]) + 1.0e-12
        self.assertFalse(supervisor.eb_arrays_are_identical(rows_3, rows_4))

    def test_figures_2_4_have_twelve_lines_no_legends_and_same_colors(self) -> None:
        rows_2 = self._synthetic_comparison_rows(2)
        rows_3 = self._synthetic_comparison_rows(3)
        rows_4 = self._synthetic_comparison_rows(4)
        shared = supervisor.comparison_ylim(rows_3, rows_4)
        figures = [
            supervisor.create_comparison_figure(
                rows_2,
                solid_key="book_slope_lambda",
                dashed_key="section_clamp_lambda",
                ylim=(0.5, 5.5),
            ),
            supervisor.create_comparison_figure(
                rows_3,
                solid_key="timoshenko_lambda",
                dashed_key="eb_lambda",
                ylim=shared,
            ),
            supervisor.create_comparison_figure(
                rows_4,
                solid_key="timoshenko_lambda",
                dashed_key="eb_lambda",
                ylim=shared,
            ),
        ]
        try:
            for figure in figures:
                self.assertEqual(len(figure.axes), 1)
                self.assertEqual(len(figure.axes[0].lines), 12)
                self.assertIsNone(figure.axes[0].get_legend())
                self.assertFalse(figure.legends)
                colors = [line.get_color() for line in figure.axes[0].lines]
                for mode, color in enumerate(supervisor.MODE_COLORS):
                    self.assertEqual(colors[2 * mode], color)
                    self.assertEqual(colors[2 * mode + 1], color)
            self.assertEqual(figures[1].axes[0].get_ylim(), figures[2].axes[0].get_ylim())
        finally:
            for figure in figures:
                plt.close(figure)

    def test_figure_1_has_four_panels_and_no_legend(self) -> None:
        rows: list[dict[str, object]] = []
        for theta in (0.0, 90.0):
            for mode in range(1, 8):
                rows.append(
                    {
                        "theta_deg": theta,
                        "book_curve_mode_label": mode,
                        "calculated_frequency_khz": 0.4 * mode + theta / 900.0,
                        "calculated_modal_loss_factor_times_100": 0.2 * mode,
                        "digitized_frequency_khz": 0.4 * mode + theta / 900.0,
                        "digitized_modal_loss_factor_times_100": 0.2 * mode,
                    }
                )
        figure = supervisor.create_figure_1(rows)
        try:
            self.assertEqual(len(figure.axes), 4)
            self.assertFalse(figure.legends)
            self.assertTrue(all(axis.get_legend() is None for axis in figure.axes))
        finally:
            plt.close(figure)

    def test_manifest_contract_lists_twelve_figures(self) -> None:
        self.assertEqual(sorted(supervisor.FIGURE_BASENAMES), list(range(1, 13)))
        self.assertEqual(sorted(supervisor.DATA_FILENAMES), list(range(1, 13)))

    def test_figures_5_and_6_exact_direct_geometry(self) -> None:
        figure_5 = supervisor.FIGURE_5_PRESET
        self.assertIs(figure_5.material_factory, hms_dx_209_material)
        self.assertEqual(
            (
                figure_5.length_1_m,
                figure_5.length_2_m,
                figure_5.a_m,
                figure_5.a_2_m,
                figure_5.b_m,
                figure_5.mu,
            ),
            (0.3, 0.5, 0.005, None, 0.020, 0.25),
        )
        figure_6 = supervisor.FIGURE_6_PRESET
        self.assertIs(figure_6.material_factory, hms_dx_209_material)
        self.assertEqual(
            (
                figure_6.length_1_m,
                figure_6.length_2_m,
                figure_6.a_m,
                figure_6.a_2_m,
                figure_6.b_m,
                figure_6.b_2_m,
                figure_6.mu,
            ),
            (0.3, 0.5, 0.004, 0.006, 0.020, 0.020, 0.25),
        )
        manifest = supervisor._preset_manifest(figure_6)
        self.assertEqual(manifest["a_1_m"], 0.004)
        self.assertEqual(manifest["a_2_m"], 0.006)
        self.assertNotIn("eta", manifest)

    def test_figures_7_and_8_exact_material_angles(self) -> None:
        figure_7 = supervisor.FIGURE_7_PRESET
        figure_8 = supervisor.FIGURE_8_PRESET
        for preset, theta in ((figure_7, 5.0), (figure_8, 15.0)):
            self.assertIs(preset.material_factory, hms_dx_209_material)
            self.assertEqual((preset.theta_1_deg, preset.theta_2_deg), (theta, theta))
            self.assertEqual((preset.a_m, preset.b_m), (0.005, 0.020))
            self.assertEqual((preset.length_1_m, preset.length_2_m), (0.4, 0.4))

    def test_figures_9_12_change_only_figure_number_and_material_angle(self) -> None:
        expected = {9: 1.0, 10: 2.0, 11: 3.0, 12: 4.0}
        excluded = {"figure_numbers", "theta_1_deg", "theta_2_deg"}
        for figure, theta in expected.items():
            preset = supervisor.SMALL_THETA_PRESETS[figure]
            self.assertEqual(preset.figure_numbers, (figure,))
            self.assertEqual((preset.theta_1_deg, preset.theta_2_deg), (theta, theta))
            self.assertIs(preset.material_factory, hms_dx_209_material)
            self.assertEqual(supervisor.SMALL_THETA_CLAMP, "book_slope_clamp")
            for field in fields(supervisor.FigurePreset):
                if field.name not in excluded:
                    self.assertEqual(
                        getattr(preset, field.name),
                        getattr(supervisor.FIGURE_7_PRESET, field.name),
                    )

    def test_fixed_reference_lambda_is_used_for_figures_5_12(self) -> None:
        root = self._synthetic_fast_result((0.0,)).spectra[0.0].roots[0]
        values = [
            supervisor._fixed_reference_lambda_for_root(root, preset)
            for preset in (
                supervisor.FIGURE_5_PRESET,
                supervisor.FIGURE_6_PRESET,
                supervisor.FIGURE_7_PRESET,
                supervisor.FIGURE_8_PRESET,
                *supervisor.SMALL_THETA_PRESETS.values(),
            )
        ]
        np.testing.assert_array_equal(values, np.full(len(values), values[0]))
        normalization = supervisor._normalization_manifest(
            supervisor.FIGURE_6_PRESET, fixed_reference=True
        )
        self.assertEqual(normalization["a0_m"], 0.005)
        self.assertEqual(normalization["b0_m"], 0.020)

    def test_figures_7_8_reuse_exact_figure_3_arrays(self) -> None:
        beta_values = np.asarray([0.0, 90.0])
        rows_3 = self._synthetic_comparison_rows(3)
        figure_7 = supervisor._extended_figure_rows(
            7,
            supervisor.FIGURE_7_PRESET,
            beta_values,
            self._synthetic_fast_result((0.0, 90.0), 5.0),
            comparison_type="diagnostic",
            left_model="Chapter2_monoclinic_Timoshenko_theta5",
            left_theta_deg=5.0,
            right_model="rectangular_orthotropic_EB_theta0_Saint_Venant_approximation",
            right_theta_deg=0.0,
            reused_figure_3_rows=rows_3,
            reused_prefix="eb",
        )
        figure_8 = supervisor._extended_figure_rows(
            8,
            supervisor.FIGURE_8_PRESET,
            beta_values,
            self._synthetic_fast_result((0.0, 90.0), 15.0),
            comparison_type="material_rotation",
            left_model="Chapter2_theta15",
            left_theta_deg=15.0,
            right_model="Chapter2_theta0",
            right_theta_deg=0.0,
            reused_figure_3_rows=rows_3,
            reused_prefix="timoshenko",
        )
        np.testing.assert_array_equal(
            [row["right_frequency_hz"] for row in figure_7],
            [row["eb_frequency_hz"] for row in rows_3],
        )
        np.testing.assert_array_equal(
            [row["right_frequency_hz"] for row in figure_8],
            [row["timoshenko_frequency_hz"] for row in rows_3],
        )
        self.assertTrue(all(row["data_origin"] == "reused_figure_03" for row in figure_7))
        self.assertTrue(all(row["data_origin"] == "reused_figure_03" for row in figure_8))

    def test_figures_9_12_reuse_exact_figure_3_eb_array_and_schema(self) -> None:
        beta_values = np.asarray([0.0, 90.0])
        rows_3 = self._synthetic_comparison_rows(3)
        source = np.asarray(
            [
                (
                    row["beta_deg"],
                    row["mode"],
                    row["eb_frequency_hz"],
                    row["eb_lambda"],
                )
                for row in rows_3
            ],
            dtype=float,
        )
        for figure, preset in supervisor.SMALL_THETA_PRESETS.items():
            rows = supervisor._extended_figure_rows(
                figure,
                preset,
                beta_values,
                self._synthetic_fast_result(
                    (0.0, 90.0), preset.theta_1_deg
                ),
                comparison_type="small_theta",
                left_model="Chapter2_monoclinic_Timoshenko",
                left_theta_deg=preset.theta_1_deg,
                right_model="rectangular_orthotropic_EB_theta0",
                right_theta_deg=0.0,
                reused_figure_3_rows=rows_3,
                reused_prefix="eb",
                relative_difference_key=(
                    "relative_difference_to_theta0_EB_baseline"
                ),
            )
            reused = np.asarray(
                [
                    (
                        row["beta_deg"],
                        row["mode"],
                        row["right_frequency_hz"],
                        row["right_lambda"],
                    )
                    for row in rows
                ],
                dtype=float,
            )
            self.assertTrue(np.array_equal(reused, source))
            self.assertTrue(
                all(
                    "relative_difference_to_theta0_EB_baseline" in row
                    and "relative_theory_difference" not in row
                    for row in rows
                )
            )

    def test_figures_5_8_plot_twelve_lines_without_legends(self) -> None:
        rows = [
            {
                "beta_deg": beta,
                "mode": mode,
                "left_lambda": 0.8 * mode + beta / 900.0,
                "right_lambda": 0.82 * mode + beta / 900.0,
            }
            for beta in (0.0, 90.0)
            for mode in range(1, supervisor.GUARD_ROOT_COUNT + 1)
        ]
        figure = supervisor.create_comparison_figure(
            rows,
            solid_key="left_lambda",
            dashed_key="right_lambda",
            ylim=supervisor._own_ylim(rows, "left_lambda", "right_lambda"),
        )
        try:
            self.assertEqual(figure.get_size_inches().tolist(), [7.2, 4.8])
            self.assertEqual(len(figure.axes[0].lines), 12)
            self.assertIsNone(figure.axes[0].get_legend())
            self.assertFalse(figure.legends)
        finally:
            plt.close(figure)

    def test_figures_9_12_plot_style_and_figure_7_y_limits(self) -> None:
        rows = [
            {
                "beta_deg": beta,
                "mode": mode,
                "left_lambda": 0.8 * mode + beta / 900.0,
                "right_lambda": 0.82 * mode + beta / 900.0,
            }
            for beta in (0.0, 90.0)
            for mode in range(1, supervisor.GUARD_ROOT_COUNT + 1)
        ]
        figure_7_limits = (1.208572405120738, 7.808086968882901)
        figures = [
            supervisor.create_comparison_figure(
                rows,
                solid_key="left_lambda",
                dashed_key="right_lambda",
                ylim=figure_7_limits,
            )
            for _ in supervisor.SMALL_THETA_FIGURE_NUMBERS
        ]
        try:
            for figure in figures:
                self.assertEqual(figure.get_size_inches().tolist(), [7.2, 4.8])
                self.assertEqual(figure.axes[0].get_ylim(), figure_7_limits)
                self.assertEqual(len(figure.axes[0].lines), 12)
                self.assertIsNone(figure.axes[0].get_legend())
                self.assertFalse(figure.legends)
        finally:
            for figure in figures:
                plt.close(figure)

    def test_theta_small_character_table_contract_and_energy_fractions(self) -> None:
        rows = []
        for theta in supervisor.SMALL_THETA_ANGLES_DEG:
            for position in range(1, supervisor.GUARD_ROOT_COUNT + 1):
                rows.append(
                    {
                        "theta_deg": theta,
                        "sorted_position": position,
                        "frequency_hz": 100.0 * position,
                        "lambda_ref": float(position),
                        "bending_fraction": 0.7,
                        "shear_fraction": 0.1,
                        "torsion_fraction": 0.2,
                        "dominant_character": "bending-like",
                        "determinant_residual": 1.0e-12,
                        "singular_residual": 2.0e-12,
                    }
                )
        supervisor._validate_theta_small_character_rows(rows)
        self.assertEqual(len(rows), 4 * 7)
        fractions = np.asarray(
            [
                [
                    row["bending_fraction"],
                    row["shear_fraction"],
                    row["torsion_fraction"],
                ]
                for row in rows
            ]
        )
        self.assertTrue(np.all(np.isfinite(fractions)))
        self.assertTrue(np.all(fractions >= 0.0))
        np.testing.assert_allclose(np.sum(fractions, axis=1), 1.0)

    def test_mode_character_evolution_reports_target_pair_changes(self) -> None:
        rows = []
        labels = {
            0.0: {3: "torsion-like", 4: "bending-like", 5: "torsion-like", 6: "bending-like", 7: "bending-like"},
            1.0: {3: "torsion-like", 4: "bending-like", 5: "bending-like", 6: "torsion-like", 7: "bending-like"},
            4.0: {3: "bending-like", 4: "torsion-like", 5: "bending-like", 6: "torsion-like", 7: "bending-like"},
        }
        for theta, by_position in labels.items():
            for position, character in by_position.items():
                rows.append(
                    {
                        "theta_deg": theta,
                        "sorted_position": position,
                        "frequency_hz": 100.0 * position,
                        "dominant_character": character,
                    }
                )
        evolution = supervisor._mode_character_evolution(rows)
        self.assertEqual(
            evolution["first_pair_order_change_theta_deg"],
            {"positions_3_4": 4.0, "positions_5_6": 1.0},
        )

    def test_small_theta_scientific_path_has_no_interpolation_or_fitting(self) -> None:
        source = inspect.getsource(supervisor._compute_extended_figure_data)
        for forbidden in ("interp", "polyfit", "curve_fit", "minimize"):
            self.assertNotIn(forbidden, source)
        self.assertIn("run_family", source)
        self.assertIn("SMALL_THETA_PRESETS", source)

    def test_figure_7_caption_classifies_diagnostic_approximation(self) -> None:
        source = inspect.getsource(supervisor._report_text)
        self.assertIn("диагностикой применимости ортотропного приближения", source)
        self.assertIn("не чистым сравнением", source)
        self.assertIn("Chapter2_theta15", inspect.getsource(supervisor._compute_extended_figure_data))

    def test_legacy_solver_path_remains_selectable(self) -> None:
        sentinel = mock.sentinel.sweep
        with mock.patch.object(supervisor, "_sweep", return_value=sentinel) as legacy:
            result = supervisor._selected_sweep(
                supervisor.FIGURE_2_PRESET,
                np.asarray([0.0]),
                "book_slope_clamp",
                solver_mode="legacy",
            )
        self.assertIs(result, sentinel)
        legacy.assert_called_once()

    def test_fast_cli_contract_and_default(self) -> None:
        arguments = supervisor.parse_args([])
        self.assertEqual(arguments.solver_mode, "fast")
        self.assertEqual(arguments.jobs, 1)
        self.assertFalse(arguments.resume)
        self.assertFalse(arguments.validate_fast_solver)
        self.assertFalse(arguments.benchmark_fast_solver)
        self.assertEqual(
            supervisor.parse_args(["--figure", "theta-small"]).figure,
            "theta-small",
        )

    def test_saved_full_oracle_validation_evidence(self) -> None:
        output_dir = supervisor.DEFAULT_OUTPUT_DIR
        benchmark_path = output_dir / supervisor.FAST_BENCHMARK_FILENAME
        validation_path = output_dir / supervisor.FAST_VALIDATION_FILENAME
        self.assertTrue(benchmark_path.is_file())
        self.assertTrue(validation_path.is_file())
        benchmark = json.loads(benchmark_path.read_text(encoding="utf-8"))
        self.assertEqual(benchmark["status"], "FAST_SOLVER_PASS")
        self.assertEqual(benchmark["beta_point_count"], 181)
        self.assertEqual(benchmark["root_count"], 7)
        self.assertEqual(benchmark["validated_row_count"], 5 * 181 * 7)
        self.assertLessEqual(benchmark["maximum_relative_frequency_error"], 1.0e-8)
        rows = supervisor._read_csv(validation_path)
        self.assertEqual(len(rows), 5 * 181 * 7)
        self.assertTrue(all(row["status"] == "PASS" for row in rows))
        keys = {(row["family_id"], float(row["beta_deg"]), int(row["mode"])) for row in rows}
        self.assertEqual(len(keys), len(rows))

    def test_sampled_oracle_rows_match_and_eb_remains_shared(self) -> None:
        rows = supervisor._read_csv(
            supervisor.DEFAULT_OUTPUT_DIR / supervisor.FAST_VALIDATION_FILENAME
        )
        sampled = [
            row
            for row in rows
            if float(row["beta_deg"]) in (0.0, 45.0, 90.0)
        ]
        self.assertEqual(len(sampled), 5 * 3 * 7)
        self.assertTrue(
            all(float(row["relative_frequency_error"]) <= 1.0e-8 for row in sampled)
        )
        rows_3 = supervisor._comparison_rows_as_numbers(
            supervisor._read_csv(supervisor.DEFAULT_OUTPUT_DIR / "figure_03_data.csv"), 3
        )
        rows_4 = supervisor._comparison_rows_as_numbers(
            supervisor._read_csv(supervisor.DEFAULT_OUTPUT_DIR / "figure_04_data.csv"), 4
        )
        self.assertTrue(supervisor.eb_arrays_are_identical(rows_3, rows_4))

    def test_reuse_and_force_are_mutually_exclusive(self) -> None:
        with self.assertRaises(SystemExit):
            supervisor.parse_args(["--reuse-data", "--force-recompute"])

    def test_reuse_data_does_not_call_root_solver(self) -> None:
        rows = self._synthetic_comparison_rows(2)
        results_root = supervisor.ROOT / "results"
        results_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=results_root) as temporary:
            output_dir = Path(temporary)
            supervisor._write_csv(output_dir / "figure_02_data.csv", rows)
            (output_dir / "plot_manifest.json").write_text(
                json.dumps(
                    {
                        "section_clamp_reference": {
                            "figure_2_book_material": {
                                "maximum_relative_frequency_difference": 0.0,
                                "status": "PASS",
                            }
                        },
                        "runtimes_seconds": {"scientific_total": 0.0},
                        "matrix_evaluation_counts": {},
                    }
                ),
                encoding="utf-8",
            )
            with mock.patch.object(
                supervisor,
                "find_elastic_roots",
                side_effect=AssertionError("root solver called under --reuse-data"),
            ):
                status = supervisor.main(
                    [
                        "--figure",
                        "2",
                        "--output-dir",
                        str(output_dir),
                        "--reuse-data",
                        "--beta-step-deg",
                        "90",
                    ]
                )
            self.assertEqual(status, 0)
            manifest = json.loads(
                (output_dir / "plot_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["execution_mode"], "reuse_data")

    def test_theta_small_reuse_data_does_not_call_scientific_solver(self) -> None:
        beta_values = np.asarray([0.0, 90.0])
        rows_3 = self._synthetic_comparison_rows(3)
        results_root = supervisor.ROOT / "results"
        results_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=results_root) as temporary:
            output_dir = Path(temporary)
            character_rows: list[dict[str, object]] = []
            reuse_checks: dict[str, bool] = {}
            for figure, preset in supervisor.SMALL_THETA_PRESETS.items():
                rows = supervisor._extended_figure_rows(
                    figure,
                    preset,
                    beta_values,
                    self._synthetic_fast_result(
                        (0.0, 90.0), preset.theta_1_deg
                    ),
                    comparison_type="small_theta",
                    left_model="Chapter2_monoclinic_Timoshenko",
                    left_theta_deg=preset.theta_1_deg,
                    right_model="rectangular_orthotropic_EB_theta0",
                    right_theta_deg=0.0,
                    reused_figure_3_rows=rows_3,
                    reused_prefix="eb",
                    relative_difference_key=(
                        "relative_difference_to_theta0_EB_baseline"
                    ),
                )
                supervisor._write_csv(
                    output_dir / supervisor.DATA_FILENAMES[figure], rows
                )
                reuse_checks[
                    f"figure_{figure}_eb_equals_figure_3_eb"
                ] = True
                for position in range(1, supervisor.GUARD_ROOT_COUNT + 1):
                    character_rows.append(
                        {
                            "theta_deg": preset.theta_1_deg,
                            "sorted_position": position,
                            "frequency_hz": 100.0 * position,
                            "lambda_ref": float(position),
                            "bending_fraction": 0.7,
                            "shear_fraction": 0.1,
                            "torsion_fraction": 0.2,
                            "dominant_character": "bending-like",
                            "determinant_residual": 1.0e-12,
                            "singular_residual": 2.0e-12,
                        }
                    )
            supervisor._write_csv(
                output_dir / supervisor.THETA_SMALL_CHARACTER_FILENAME,
                character_rows,
            )
            (output_dir / "plot_manifest.json").write_text(
                json.dumps(
                    {
                        "y_limits": {
                            supervisor.FIGURE_7_Y_LIMITS_KEY: [1.0, 8.0]
                        },
                        "extended_fast_solver": {
                            "reuse_checks": reuse_checks,
                            "families": {},
                        },
                        "runtimes_seconds": {"scientific_total": 0.0},
                        "matrix_evaluation_counts": {},
                    }
                ),
                encoding="utf-8",
            )
            with mock.patch.object(
                supervisor,
                "_run_fast_family",
                side_effect=AssertionError(
                    "scientific solver called under --reuse-data"
                ),
            ), mock.patch.object(
                supervisor,
                "find_elastic_roots",
                side_effect=AssertionError("root solver called under --reuse-data"),
            ):
                status = supervisor.main(
                    [
                        "--figure",
                        "theta-small",
                        "--output-dir",
                        str(output_dir),
                        "--reuse-data",
                        "--beta-step-deg",
                        "90",
                        "--jobs",
                        "1",
                    ]
                )
            self.assertEqual(status, 0)
            manifest = json.loads(
                (output_dir / "plot_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["execution_mode"], "reuse_data")
            self.assertEqual(manifest["small_theta_figures_status"], "PASS")

    def test_figure_1_uses_canonical_csv_without_fitting_logic(self) -> None:
        source = inspect.getsource(supervisor._build_figure_1_rows)
        self.assertIn("CANONICAL_FIGURE_1_ROOTS", source)
        self.assertIn("CANONICAL_FIGURE_1_DIGITIZED", source)
        self.assertNotIn("curve_fit", source)
        self.assertNotIn("minimize", source)
        self.assertNotIn("BookMaterial(", source)

    def test_default_output_is_outside_article_workspaces(self) -> None:
        resolved = supervisor._safe_output_dir(supervisor.DEFAULT_OUTPUT_DIR)
        for workspace in supervisor.ARTICLE_WORKSPACES:
            with self.assertRaises(ValueError):
                resolved.relative_to(workspace.resolve())

    def test_parallel_article_workspace_has_no_git_changes(self) -> None:
        completed = subprocess.run(
            [
                "git",
                "status",
                "--short",
                "--",
                "paper_dorofeev_style",
                "paper_thickness_mismatch_timoshenko",
            ],
            cwd=supervisor.ROOT,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
        )
        self.assertEqual(completed.stdout.strip(), "")


if __name__ == "__main__":
    unittest.main()
