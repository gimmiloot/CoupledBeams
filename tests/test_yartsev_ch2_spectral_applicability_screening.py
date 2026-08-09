from __future__ import annotations

import ast
import subprocess
import sys
import tempfile
import unittest
from dataclasses import asdict, fields
from pathlib import Path
from unittest import mock

import numpy as np

from scripts.analysis.anisotropic_rods import (
    screen_yartsev_ch2_spectral_applicability as screening,
)


ROOT = Path(__file__).resolve().parents[1]


class YartsevChapter2SpectralApplicabilityScreeningTest(unittest.TestCase):
    def test_import_has_no_calculation_or_plotting_side_effect(self) -> None:
        command = (
            "import matplotlib.pyplot as plt; "
            "import scripts.analysis.anisotropic_rods."
            "screen_yartsev_ch2_spectral_applicability; "
            "print(len(plt.get_fignums()))"
        )
        completed = subprocess.run(
            [sys.executable, "-B", "-c", command],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.stdout.strip(), "0")

    def test_fixed_parameter_grids(self) -> None:
        self.assertEqual(screening.MU_VALUES, (0.0, 0.25, 0.5))
        self.assertEqual(screening.TAU_VALUES, (-0.2, 0.0, 0.2))
        self.assertEqual(len(screening.screening_geometries()), 9)
        np.testing.assert_array_equal(
            screening.beta_grid(), np.arange(0.0, 91.0, 5.0)
        )
        self.assertEqual(len(screening.beta_grid()), 19)

    def test_geometry_definitions_and_volume_gate(self) -> None:
        for case in screening.screening_geometries():
            recovered_tau = (case.a2_m - case.a1_m) / (
                case.a1_m + case.a2_m
            )
            self.assertAlmostEqual(recovered_tau, case.tau, places=15)
            self.assertAlmostEqual(
                case.a2_over_a1,
                (1.0 + case.tau) / (1.0 - case.tau),
                places=14,
            )
            self.assertLessEqual(case.volume_relative_residual, 1.0e-12)
            self.assertAlmostEqual(case.b1_m / case.a1_m, 4.0, places=14)
            self.assertAlmostEqual(case.b2_m / case.a2_m, 4.0, places=14)
            self.assertAlmostEqual(case.L1_m + case.L2_m, 0.8, places=15)
            self.assertAlmostEqual(case.mass_kg, 0.1264, places=14)

    def test_baseline_geometry_is_exact_reference_case(self) -> None:
        case = screening.make_geometry(0.0, 0.0)
        self.assertEqual(case.L1_m, 0.4)
        self.assertEqual(case.L2_m, 0.4)
        self.assertEqual(case.a1_m, 0.005)
        self.assertEqual(case.a2_m, 0.005)
        self.assertEqual(case.b1_m, 0.020)
        self.assertEqual(case.b2_m, 0.020)
        self.assertEqual(case.volume_m3, 8.0e-5)

    def test_models_angles_and_root_counts_are_frozen(self) -> None:
        self.assertEqual(screening.MODEL_IDS, ("T2", "T0", "EB0"))
        self.assertEqual(
            screening.MODEL_THETA_DEG, {"T2": 2.0, "T0": 0.0, "EB0": 0.0}
        )
        self.assertEqual(screening.GUARD_ROOT_COUNT, 7)
        self.assertEqual(screening.COMPARED_ROOT_COUNT, 6)
        self.assertEqual(screening.APPLICABILITY_THRESHOLD, 0.10)

    def test_presets_pass_both_real_arm_geometries_directly(self) -> None:
        case = screening.make_geometry(0.5, -0.2)
        for model_id in screening.MODEL_IDS:
            preset = screening.model_preset(case, model_id)
            self.assertEqual(preset.a_m, case.a1_m)
            self.assertEqual(preset.b_m, case.b1_m)
            self.assertEqual(preset.a_2_m, case.a2_m)
            self.assertEqual(preset.b_2_m, case.b2_m)
            self.assertEqual(preset.length_1_m, case.L1_m)
            self.assertEqual(preset.length_2_m, case.L2_m)
            self.assertEqual(preset.material_name, "HMS/DX-209")
            self.assertEqual(preset.material_mode, "elastic")

    def test_fixed_lambda_normalization_is_geometry_and_model_independent(self) -> None:
        frequency = np.asarray([100.0, 250.0])
        reference = screening.fixed_lambda_from_frequency(frequency)
        for case in screening.screening_geometries():
            for model_id in screening.MODEL_IDS:
                np.testing.assert_array_equal(
                    screening.fixed_lambda_from_frequency(frequency), reference
                )
                contract = screening.normalization_contract(case, model_id)
                self.assertEqual(contract["A0_m2"], 0.005 * 0.020)
                self.assertEqual(
                    contract["I_y0_m4"], 0.005**3 * 0.020 / 12.0
                )
                self.assertEqual(contract["E_x0_pa"], 191.0e9)
                self.assertTrue(contract["fixed_for_all_geometries_and_models"])

    def test_metric_denominators_and_gap_definitions(self) -> None:
        t2 = np.asarray([10, 20, 30, 40, 50, 60, 70], dtype=float)
        t0 = np.asarray([11, 22, 33, 44, 55, 66, 77], dtype=float)
        eb0 = np.asarray([9, 18, 27, 36, 45, 54, 63], dtype=float)
        metrics = screening.point_metrics_from_lambdas(t2, t0, eb0)
        expected_beam = np.abs(t0[:6] - eb0[:6]) / t0[:6]
        expected_orient = np.abs(t2[:6] - t0[:6]) / t2[:6]
        expected_simpl = np.abs(t2[:6] - eb0[:6]) / t2[:6]
        self.assertEqual(metrics["Delta_beam"], float(np.max(expected_beam)))
        self.assertEqual(metrics["Delta_orient"], float(np.max(expected_orient)))
        self.assertEqual(metrics["Delta_simpl"], float(np.max(expected_simpl)))
        gaps_0 = 2.0 * np.diff(t0) / (t0[1:] + t0[:-1])
        gaps_2 = 2.0 * np.diff(t2) / (t2[1:] + t2[:-1])
        self.assertEqual(metrics["g_min_0"], float(np.min(gaps_0)))
        self.assertEqual(metrics["g_min_2"], float(np.min(gaps_2)))
        self.assertEqual(
            metrics["G_nearest"],
            float(gaps_2[int(np.argmin(gaps_0))] - np.min(gaps_0)),
        )
        self.assertEqual(metrics["G_open"], float(np.max(gaps_2 - gaps_0)))
        self.assertEqual(metrics["G_close"], float(np.min(gaps_2 - gaps_0)))

    def test_no_old_thickness_parameter_or_circular_article_import(self) -> None:
        source_path = Path(screening.__file__).resolve()
        tree = ast.parse(source_path.read_text(encoding="utf-8"))
        identifiers = {
            node.id for node in ast.walk(tree) if isinstance(node, ast.Name)
        }
        arguments = {
            argument.arg
            for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            for argument in (*node.args.args, *node.args.kwonlyargs)
        }
        self.assertNotIn("eta", identifiers | arguments)
        self.assertNotIn("eta", {field.name for field in fields(screening.GeometryCase)})
        imported_modules = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        } | {
            node.module or ""
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
        }
        forbidden_fragments = (
            "src.my_project.analytic",
            "src.my_project.fem",
            "thickness_mismatch",
            "paper_dorofeev_style",
            "paper_thickness_mismatch_timoshenko",
        )
        self.assertFalse(
            any(
                fragment in module
                for fragment in forbidden_fragments
                for module in imported_modules
            )
        )

    def test_reuse_baseline_exactly_matches_figures_3_and_10(self) -> None:
        case = screening.make_geometry(0.0, 0.0)
        beta = np.asarray([0.0, 5.0, 90.0])
        for model_id, figure in (("T2", 10), ("T0", 3), ("EB0", 3)):
            rows, audit = screening.load_reused_family(case, model_id, beta)
            self.assertEqual(len(rows), 21)
            self.assertEqual(audit["source_figure"], figure)
            self.assertTrue(audit["beta_frequency_lambda_exactly_equal"])
            self.assertEqual(audit["status"], "PASS")

    def test_reuse_unequal_length_exactly_matches_figure_5(self) -> None:
        case = screening.make_geometry(0.25, 0.0)
        beta = np.asarray([0.0, 45.0, 90.0])
        for model_id in ("T0", "EB0"):
            rows, audit = screening.load_reused_family(case, model_id, beta)
            self.assertEqual(len(rows), 21)
            self.assertEqual(audit["source_figure"], 5)
            self.assertTrue(audit["beta_frequency_lambda_exactly_equal"])

    def test_figure_6_is_never_an_approved_reuse_source(self) -> None:
        self.assertNotIn(6, {spec.figure_number for spec in screening.REUSE_SPECS})
        self.assertEqual(len(screening.REUSE_SPECS), 5)

    def test_cli_data_modes_are_mutually_exclusive_and_smoke_is_separate(self) -> None:
        with self.assertRaises(SystemExit):
            screening.parse_args(["--reuse-data", "--force-recompute"])
        self.assertEqual(len(screening.screening_geometries(smoke=True)), 3)
        np.testing.assert_array_equal(
            screening.beta_grid(smoke=True), np.asarray([0.0, 45.0, 90.0])
        )
        self.assertNotEqual(
            screening.DEFAULT_SMOKE_OUTPUT_DIR, screening.DEFAULT_OUTPUT_DIR
        )

    def test_reuse_data_path_never_calls_scientific_solver(self) -> None:
        results_root = ROOT / "results"
        with tempfile.TemporaryDirectory(
            prefix="ap0_reuse_test_", dir=results_root
        ) as temporary:
            output_dir = Path(temporary)
            cases = screening.screening_geometries()
            screening._write_csv(
                output_dir / screening.GEOMETRY_FILENAME,
                [asdict(case) for case in cases],
            )
            spectra = []
            for case_index, case in enumerate(cases):
                for beta in screening.beta_grid():
                    for model_index, model_id in enumerate(screening.MODEL_IDS):
                        for mode in range(1, 8):
                            point_variation = (
                                1.0e-4
                                * (case_index + beta / 5.0)
                                * mode**2
                            )
                            spectra.append(
                                {
                                    "geometry_id": case.geometry_id,
                                    "mu": case.mu,
                                    "tau": case.tau,
                                    "beta_deg": beta,
                                    "model_id": model_id,
                                    "theta_deg": screening.MODEL_THETA_DEG[model_id],
                                    "mode": mode,
                                    "root_role": "plotted" if mode <= 6 else "guard",
                                    "frequency_hz": 100.0 * mode + model_index,
                                    "lambda_ref": (
                                        1.0 * mode
                                        + point_variation
                                        + 0.01
                                        * model_index
                                        * (1.0 + beta / 90.0)
                                    ),
                                    "quality_status": "PASS",
                                    "quality_basis": "scaled",
                                    "root_status": "accepted",
                                    "accepted_determinant_residual": 1.0e-12,
                                    "accepted_relative_singular_residual": 1.0e-12,
                                    "data_origin": "saved_synthetic_test_data",
                                    "global_fallback_reason": "",
                                }
                            )
            screening._write_csv(output_dir / screening.SPECTRA_FILENAME, spectra)
            screening._write_csv(
                output_dir / screening.REUSE_AUDIT_FILENAME,
                [{"status": "PASS", "source_figure": spec.figure_number} for spec in screening.REUSE_SPECS],
            )

            def fake_figures(_summaries, _points, path):
                paths = []
                for basename in (
                    screening.FIGURE_S1_BASENAME,
                    screening.FIGURE_S2_BASENAME,
                    screening.FIGURE_S3_BASENAME,
                ):
                    for suffix in ("pdf", "png"):
                        figure_path = path / f"{basename}.{suffix}"
                        figure_path.touch()
                        paths.append(str(figure_path.relative_to(ROOT)))
                return paths

            with mock.patch.object(
                screening,
                "solve_new_family",
                side_effect=AssertionError("scientific solver called"),
            ), mock.patch.object(
                screening, "create_diagnostic_figures", side_effect=fake_figures
            ):
                result = screening.main(
                    ["--output-dir", str(output_dir), "--reuse-data"]
                )
            self.assertEqual(result, 0)
            self.assertEqual(
                screening._read_csv(output_dir / screening.SPECTRA_FILENAME),
                screening._read_csv(output_dir / screening.SPECTRA_FILENAME),
            )

    def test_supervisor_hashes_are_unchanged_by_reuse_reads(self) -> None:
        before = screening.supervisor_data_hashes()
        screening.load_reused_family(
            screening.make_geometry(0.0, 0.0), "T2", np.asarray([0.0])
        )
        screening.load_reused_family(
            screening.make_geometry(0.25, 0.0), "EB0", np.asarray([90.0])
        )
        self.assertEqual(screening.supervisor_data_hashes(), before)

    def test_workflow_pass_is_independent_of_scientific_classification(self) -> None:
        cases = screening.screening_geometries()
        beta = screening.beta_grid()
        spectra = [
            {"quality_status": "PASS"}
            for _ in range(9 * 19 * 3 * 7)
        ]
        points = [{} for _ in range(9 * 19)]
        summaries = [
            {"screening_classification": "EXCEEDS_10_PERCENT_ON_SCREENING_GRID"}
            for _ in range(9)
        ]
        with tempfile.TemporaryDirectory(
            prefix="ap0_status_test_", dir=ROOT / "results"
        ) as temporary:
            output_dir = Path(temporary)
            for path in screening._expected_output_paths(output_dir):
                path.touch()
            status = screening.determine_workflow_status(
                smoke=False,
                geometries=cases,
                beta_values=beta,
                spectra=spectra,
                point_rows=points,
                summaries=summaries,
                reuse_audit=[{"status": "PASS"} for _ in range(5)],
                supervisor_hashes_unchanged=True,
                output_dir=output_dir,
            )
        self.assertEqual(status, "PASS")


if __name__ == "__main__":
    unittest.main()
