from __future__ import annotations

import ast
import shutil
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
    @staticmethod
    def _synthetic_ap1_spectra() -> list[dict[str, object]]:
        source = screening.spectra_as_numbers(
            screening._read_csv(screening.AP0_SPECTRA_PATH)
        )
        rows: list[dict[str, object]] = []
        for source_row in source:
            row = dict(source_row)
            if row["model_id"] == "T2":
                row["model_id"] = "T5"
                row["theta_deg"] = 5.0
                row["frequency_hz"] = 1.02 * float(row["frequency_hz"])
                row["lambda_ref"] = float(
                    screening.fixed_lambda_from_frequency(row["frequency_hz"])
                )
            row["source_file"] = "synthetic_unit_test"
            row["source_sha256"] = "synthetic_unit_test"
            rows.append(row)
        rows.sort(
            key=lambda row: (
                str(row["geometry_id"]),
                float(row["beta_deg"]),
                screening.AP1_MODEL_IDS.index(str(row["model_id"])),
                int(row["mode"]),
            )
        )
        return rows

    @staticmethod
    def _synthetic_ap2_spectra(theta: int) -> list[dict[str, object]]:
        source = screening.spectra_as_numbers(
            screening._read_csv(screening.AP0_SPECTRA_PATH)
        )
        model_ids = screening.AP2_MODEL_IDS[theta]
        rows: list[dict[str, object]] = []
        factor = 1.0 + 0.005 * theta
        for source_row in source:
            row = dict(source_row)
            if row["model_id"] == "T2":
                row["model_id"] = f"T{theta}"
                row["theta_deg"] = float(theta)
                row["frequency_hz"] = factor * float(row["frequency_hz"])
                row["lambda_ref"] = float(
                    screening.fixed_lambda_from_frequency(row["frequency_hz"])
                )
            row["source_file"] = "synthetic_ap2_unit_test"
            row["source_sha256"] = "synthetic_ap2_unit_test"
            rows.append(row)
        rows.sort(
            key=lambda row: (
                str(row["geometry_id"]),
                float(row["beta_deg"]),
                model_ids.index(str(row["model_id"])),
                int(row["mode"]),
            )
        )
        return rows

    @classmethod
    def _synthetic_sampled_sources(
        cls,
    ) -> dict[int, dict[str, list[dict[str, object]]]]:
        sources: dict[int, dict[str, list[dict[str, object]]]] = {
            2: {
                "spectra": screening._read_csv(screening.AP0_SPECTRA_PATH),
                "points": screening._read_csv(screening.AP0_POINT_METRICS_PATH),
                "summaries": screening._read_csv(
                    screening.AP0_GEOMETRY_SUMMARY_PATH
                ),
            },
            5: {
                "spectra": screening._read_csv(
                    screening.AP1_OUTPUT_DIR / screening.SPECTRA_FILENAME
                ),
                "points": screening._read_csv(
                    screening.AP1_OUTPUT_DIR / screening.POINT_METRICS_FILENAME
                ),
                "summaries": screening._read_csv(
                    screening.AP1_OUTPUT_DIR / screening.GEOMETRY_SUMMARY_FILENAME
                ),
            },
        }
        geometries = screening.screening_geometries()
        for theta in (3, 4):
            spectra = cls._synthetic_ap2_spectra(theta)
            points, _ = screening.build_ap2_metrics(
                spectra, geometries, screening.beta_grid(), theta
            )
            summaries = screening.build_ap2_geometry_summary(
                points, geometries, theta
            )
            sources[theta] = {
                "spectra": spectra,
                "points": points,
                "summaries": summaries,
            }
        return sources

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

    def test_ap1_configuration_preserves_ap0_default_and_selects_theta5(self) -> None:
        default_args = screening.parse_args([])
        self.assertEqual(default_args.target_theta_deg, 2.0)
        ap0 = screening.screening_configuration(default_args.target_theta_deg)
        ap1 = screening.screening_configuration(5.0)
        self.assertEqual(ap0.primary_model_id, "T2")
        self.assertEqual(ap1.primary_model_id, "T5")
        self.assertEqual(ap0.output_dir.name, "theta2_small_grid")
        self.assertEqual(ap1.output_dir.name, "theta5_small_grid")
        self.assertEqual(ap1.smoke_output_dir.name, "theta5_smoke")

    def test_ap1_geometry_and_beta_grid_are_identical_to_ap0(self) -> None:
        source = screening.load_ap0_geometries(smoke=False)
        canonical = screening.screening_geometries()
        self.assertEqual(source, canonical)
        self.assertEqual(len(source), 9)
        np.testing.assert_array_equal(
            screening.beta_grid(), np.arange(0.0, 91.0, 5.0)
        )

    def test_ap1_models_angles_roots_and_normalization(self) -> None:
        self.assertEqual(screening.AP1_MODEL_IDS, ("T5", "T0", "EB0"))
        self.assertEqual(
            screening.AP1_MODEL_THETA_DEG,
            {"T5": 5.0, "T0": 0.0, "EB0": 0.0},
        )
        self.assertEqual(screening.GUARD_ROOT_COUNT, 7)
        self.assertEqual(screening.COMPARED_ROOT_COUNT, 6)
        self.assertEqual(screening.APPLICABILITY_THRESHOLD, 0.10)
        contract = screening.normalization_contract(
            screening.make_geometry(0.5, 0.2), "T5"
        )
        self.assertEqual(contract["A0_m2"], screening.REFERENCE_AREA_M2)
        self.assertEqual(contract["I_y0_m4"], screening.REFERENCE_IY_M4)
        self.assertEqual(contract["E_x0_pa"], 191.0e9)

    def test_ap1_reuses_all_ap0_t0_eb0_families_exactly(self) -> None:
        beta = screening.beta_grid()
        audits = []
        for case in screening.load_ap0_geometries(smoke=False):
            for model_id in ("T0", "EB0"):
                rows, audit = screening.load_ap1_ap0_family(case, model_id, beta)
                self.assertEqual(len(rows), 19 * 7)
                self.assertTrue(audit["frequency_array_equal"])
                self.assertTrue(audit["lambda_array_equal"])
                self.assertTrue(audit["quality_status_equal"])
                self.assertEqual(audit["reuse_status"], "PASS")
                audits.append(audit)
        self.assertEqual(len(audits), 18)

    def test_ap1_reuses_figure7_t5_subset_exactly(self) -> None:
        rows, audit = screening.load_ap1_figure7_t5_family(
            screening.make_geometry(0.0, 0.0), screening.beta_grid()
        )
        self.assertEqual(len(rows), 19 * 7)
        self.assertEqual(audit["model_id"], "T5")
        self.assertTrue(audit["frequency_array_equal"])
        self.assertTrue(audit["lambda_array_equal"])
        self.assertEqual(audit["reuse_status"], "PASS")
        self.assertEqual(
            audit["source_file"],
            str(screening.FIGURE_07_PATH.relative_to(ROOT)),
        )

    def test_ap1_expected_family_and_row_accounting(self) -> None:
        self.assertEqual(8 + 1 + 9 + 9, 27)
        self.assertEqual(8 * 19 * 7, 1064)
        self.assertEqual(9 * 19 * 3 * 7, 3591)
        self.assertEqual(9 * 19, 171)
        self.assertEqual(9 * 19 * 6, 1026)

    def test_ap1_metric_denominators_are_t5_and_delta_beam_is_t0(self) -> None:
        t5 = np.asarray([10, 20, 30, 40, 50, 60, 70], dtype=float)
        t0 = np.asarray([11, 22, 33, 44, 55, 66, 77], dtype=float)
        eb0 = np.asarray([9, 18, 27, 36, 45, 54, 63], dtype=float)
        metrics = screening.ap1_point_metrics_from_lambdas(t5, t0, eb0)
        self.assertEqual(
            metrics["Delta_beam"],
            float(np.max(np.abs(t0[:6] - eb0[:6]) / t0[:6])),
        )
        self.assertEqual(
            metrics["Delta_orient_5"],
            float(np.max(np.abs(t5[:6] - t0[:6]) / t5[:6])),
        )
        self.assertEqual(
            metrics["Delta_simpl_5"],
            float(np.max(np.abs(t5[:6] - eb0[:6]) / t5[:6])),
        )

    def test_ap1_synthetic_metrics_have_exact_counts_and_same_pair_gaps(self) -> None:
        spectra = self._synthetic_ap1_spectra()
        geometries = screening.load_ap0_geometries(smoke=False)
        points, pairs, comparisons = screening.build_ap1_metrics(
            spectra, geometries, screening.beta_grid()
        )
        self.assertEqual(len(spectra), 3591)
        self.assertEqual(len(points), 171)
        self.assertEqual(len(pairs), 1026)
        self.assertEqual(len(comparisons), 171)
        self.assertTrue(all(row["Delta_beam_array_equal"] for row in comparisons))
        self.assertTrue(
            all(
                int(row["pair_upper_mode"]) == int(row["pair_lower_mode"]) + 1
                for row in pairs
            )
        )
        first = pairs[0]
        self.assertAlmostEqual(
            first["delta_g_theta5"],
            first["g_pair_T5"] - first["g_pair_T0"],
        )

    def test_ap1_transitions_use_saved_ap0_classification(self) -> None:
        spectra = self._synthetic_ap1_spectra()
        geometries = screening.load_ap0_geometries(smoke=False)
        points, _, _ = screening.build_ap1_metrics(
            spectra, geometries, screening.beta_grid()
        )
        summaries, comparisons = screening.build_ap1_geometry_summary(
            points, geometries
        )
        ap0 = {
            row["geometry_id"]: row["screening_classification"]
            for row in screening._read_csv(screening.AP0_GEOMETRY_SUMMARY_PATH)
        }
        self.assertEqual(len(comparisons), 9)
        for row in summaries:
            self.assertEqual(row["classification_theta2"], ap0[row["geometry_id"]])
            self.assertEqual(
                row["classification_transition"],
                screening._classification_transition(
                    ap0[row["geometry_id"]],
                    row["screening_classification_theta5"],
                ),
            )

    def test_ap1_checkpoint_fingerprint_is_theta_specific(self) -> None:
        case = screening.make_geometry(0.5, 0.2)
        beta = screening.beta_grid(smoke=True)
        t2 = screening.supervisor._fast_family_fingerprint(
            screening.model_preset(case, "T2"),
            beta,
            screening.MODEL_SOLVER_PATH["T2"],
            lambda_normalization=screening.normalization_contract(case, "T2"),
        )
        t5 = screening.supervisor._fast_family_fingerprint(
            screening.model_preset(case, "T5"),
            beta,
            screening.MODEL_SOLVER_PATH["T5"],
            lambda_normalization=screening.normalization_contract(
                case, "T5", source_hashes={"source": "theta5"}
            ),
        )
        self.assertNotEqual(t2, t5)
        self.assertNotEqual(
            screening._family_id(case, "T2", smoke=False),
            screening._family_id(
                case, "T5", smoke=False, stage_tag="ap1"
            ),
        )

    def test_ap1_has_no_forbidden_scientific_calls_or_identifiers(self) -> None:
        tree = ast.parse(Path(screening.__file__).read_text(encoding="utf-8"))
        identifiers = {
            node.id for node in ast.walk(tree) if isinstance(node, ast.Name)
        }
        self.assertNotIn("eta", identifiers)
        self.assertNotIn("Pi", identifiers)
        calls = []
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            if isinstance(node.func, ast.Name):
                calls.append(node.func.id.lower())
            elif isinstance(node.func, ast.Attribute):
                calls.append(node.func.attr.lower())
        self.assertFalse(any("mac" in call for call in calls))
        self.assertFalse(any("shape" in call for call in calls))
        self.assertFalse(any("energy" in call for call in calls))

    def test_ap1_reuse_data_does_not_call_solver(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="ap1_reuse_test_", dir=ROOT / "results"
        ) as temporary:
            output_dir = Path(temporary)
            geometries = screening.load_ap0_geometries(smoke=False)
            shutil.copyfile(
                screening.AP0_GEOMETRY_PATH,
                output_dir / screening.GEOMETRY_FILENAME,
            )
            spectra = self._synthetic_ap1_spectra()
            points, pairs, point_comparisons = screening.build_ap1_metrics(
                spectra, geometries, screening.beta_grid()
            )
            summaries, geometry_comparisons = screening.build_ap1_geometry_summary(
                points, geometries
            )
            audits = []
            for case in geometries:
                for model_id in ("T0", "EB0"):
                    _, audit = screening.load_ap1_ap0_family(
                        case, model_id, screening.beta_grid()
                    )
                    audits.append(audit)
            _, baseline_audit = screening.load_ap1_figure7_t5_family(
                screening.make_geometry(0, 0), screening.beta_grid()
            )
            audits.append(baseline_audit)
            candidates = screening.select_ap1_candidates(summaries, pairs)
            screening._write_csv(output_dir / screening.SPECTRA_FILENAME, spectra)
            screening._write_csv(
                output_dir / screening.POINT_METRICS_FILENAME, points
            )
            screening._write_csv(
                output_dir / screening.GEOMETRY_SUMMARY_FILENAME, summaries
            )
            screening._write_csv(output_dir / screening.PAIRWISE_FILENAME, pairs)
            screening._write_csv(
                output_dir / screening.POINT_COMPARISON_FILENAME,
                point_comparisons,
            )
            screening._write_csv(
                output_dir / screening.GEOMETRY_COMPARISON_FILENAME,
                geometry_comparisons,
            )
            screening._write_csv(
                output_dir / screening.REUSE_AUDIT_FILENAME, audits
            )
            screening._json_write(
                output_dir / screening.CANDIDATES_FILENAME, candidates
            )

            def fake_figures(_summaries, _comparisons, _pairs, path):
                paths = []
                for basename in (
                    screening.AP1_FIGURE_S1_BASENAME,
                    screening.AP1_FIGURE_S2_BASENAME,
                    screening.AP1_FIGURE_S3_BASENAME,
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
                screening, "create_ap1_figures", side_effect=fake_figures
            ):
                result = screening.main(
                    [
                        "--target-theta-deg",
                        "5",
                        "--output-dir",
                        str(output_dir),
                        "--reuse-data",
                    ]
                )
            self.assertEqual(result, 0)

    def test_ap1_computational_status_ignores_scientific_exceedance(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="ap1_status_test_", dir=ROOT / "results"
        ) as temporary:
            output_dir = Path(temporary)
            for path in screening._ap1_expected_output_paths(output_dir):
                path.touch()
            status = screening.determine_ap1_status(
                smoke=False,
                geometries=screening.screening_geometries(),
                beta_values=screening.beta_grid(),
                spectra=[{"quality_status": "PASS"} for _ in range(3591)],
                point_rows=[{} for _ in range(171)],
                pair_rows=[{} for _ in range(1026)],
                point_comparisons=[
                    {"Delta_beam_array_equal": True} for _ in range(171)
                ],
                summaries=[
                    {
                        "screening_classification_theta5":
                        "EXCEEDS_10_PERCENT_ON_SCREENING_GRID"
                    }
                    for _ in range(9)
                ],
                geometry_comparisons=[{} for _ in range(9)],
                reuse_audit=[{"reuse_status": "PASS"} for _ in range(19)],
                ap0_unchanged=True,
                supervisor_unchanged=True,
                output_dir=output_dir,
            )
        self.assertEqual(status, "PASS")

    def test_ap1_reads_do_not_change_ap0_or_supervisor_hashes(self) -> None:
        ap0_before = screening.preservation_hashes(
            screening.AP0_OUTPUT_DIR, screening.AP0_PRESERVATION_FILENAMES
        )
        supervisor_before = screening.supervisor_preservation_hashes()
        screening.load_ap1_ap0_family(
            screening.make_geometry(0.5, -0.2), "T0", np.asarray([0.0])
        )
        screening.load_ap1_figure7_t5_family(
            screening.make_geometry(0, 0), np.asarray([90.0])
        )
        self.assertEqual(
            screening.preservation_hashes(
                screening.AP0_OUTPUT_DIR,
                screening.AP0_PRESERVATION_FILENAMES,
            ),
            ap0_before,
        )
        self.assertEqual(
            screening.supervisor_preservation_hashes(), supervisor_before
        )

    def test_ap2_cli_accepts_theta3_theta4_and_preserves_other_defaults(self) -> None:
        for theta in (3.0, 4.0):
            args = screening.parse_args(["--target-theta-deg", str(theta)])
            config = screening.screening_configuration(args.target_theta_deg)
            self.assertTrue(config.is_ap2)
            self.assertEqual(config.primary_model_id, f"T{int(theta)}")
            self.assertEqual(config.output_dir.name, f"theta{int(theta)}_small_grid")
            self.assertEqual(config.smoke_output_dir.name, f"theta{int(theta)}_smoke")
        self.assertFalse(screening.parse_args([]).combine_sampled_theta)
        self.assertTrue(
            screening.parse_args(["--combine-sampled-theta"]).combine_sampled_theta
        )

    def test_ap2_geometry_manifest_copy_is_byte_equal_to_ap0(self) -> None:
        with tempfile.TemporaryDirectory(
            prefix="ap2_geometry_copy_", dir=ROOT / "results"
        ) as temporary:
            target = Path(temporary) / screening.GEOMETRY_FILENAME
            shutil.copyfile(screening.AP0_GEOMETRY_PATH, target)
            self.assertEqual(target.read_bytes(), screening.AP0_GEOMETRY_PATH.read_bytes())
        self.assertEqual(
            screening.load_ap0_geometries(smoke=False),
            screening.screening_geometries(),
        )

    def test_ap2_models_and_fixed_lambda_normalization(self) -> None:
        frequency = np.asarray([75.0, 150.0, 300.0])
        expected = (
            screening.MATERIAL_DENSITY_KG_M3
            * screening.REFERENCE_AREA_M2
            * (2.0 * np.pi * frequency) ** 2
            * screening.REFERENCE_LENGTH_M**4
            / (screening.REFERENCE_EX_PA * screening.REFERENCE_IY_M4)
        ) ** 0.25
        np.testing.assert_allclose(
            screening.fixed_lambda_from_frequency(frequency),
            expected,
            rtol=2.0e-15,
            atol=0.0,
        )
        for theta in (3, 4):
            self.assertEqual(
                screening.AP2_MODEL_IDS[theta], (f"T{theta}", "T0", "EB0")
            )
            contract = screening.normalization_contract(
                screening.make_geometry(0.5, -0.2), f"T{theta}"
            )
            self.assertEqual(contract["E_x0_pa"], 191.0e9)
            self.assertEqual(
                contract["screening_fingerprint_context"]["theta_deg"],
                float(theta),
            )

    def test_ap2_reuses_all_t0_eb0_and_supervisor_baselines_exactly(self) -> None:
        beta = screening.beta_grid()
        for case in screening.screening_geometries():
            for model_id in ("T0", "EB0"):
                rows, audit = screening.load_ap1_ap0_family(case, model_id, beta)
                self.assertEqual(len(rows), 19 * 7)
                self.assertEqual(audit["reuse_status"], "PASS")
        for theta, expected_figure in ((3, 11), (4, 12)):
            source = screening.locate_ap2_supervisor_baseline(theta)
            self.assertEqual(source.figure_number, expected_figure)
            rows, audit = screening.load_ap2_supervisor_family(
                screening.make_geometry(0.0, 0.0), f"T{theta}", beta
            )
            self.assertEqual(len(rows), 19 * 7)
            self.assertEqual(audit["source_figure"], expected_figure)
            self.assertTrue(audit["frequency_array_equal"])
            self.assertTrue(audit["lambda_array_equal"])
            self.assertTrue(audit["quality_status_equal"])

    def test_ap2_supervisor_identity_mismatch_is_rejected_for_recompute(self) -> None:
        with mock.patch.object(
            screening,
            "load_ap2_supervisor_family",
            side_effect=RuntimeError("synthetic identity mismatch"),
        ):
            rows, audit = screening.attempt_ap2_supervisor_reuse(
                screening.make_geometry(0.0, 0.0),
                "T3",
                screening.beta_grid(smoke=True),
            )
        self.assertIsNone(rows)
        self.assertEqual(
            audit["reuse_status"], "SUPERVISOR_BASELINE_REUSE_REJECTED"
        )
        self.assertIn("synthetic identity mismatch", audit["reuse_rejection_reason"])

    def test_ap2_metric_denominators_counts_and_pairwise_identities(self) -> None:
        geometries = screening.screening_geometries()
        for theta in (3, 4):
            spectra = self._synthetic_ap2_spectra(theta)
            screening.validate_spectra(
                spectra,
                geometries,
                screening.beta_grid(),
                model_ids=screening.AP2_MODEL_IDS[theta],
            )
            points, pairs = screening.build_ap2_metrics(
                spectra, geometries, screening.beta_grid(), theta
            )
            summaries = screening.build_ap2_geometry_summary(
                points, geometries, theta
            )
            self.assertEqual(len(spectra), 3591)
            self.assertEqual(len(points), 171)
            self.assertEqual(len(summaries), 9)
            self.assertEqual(len(pairs), 1026)
            first = pairs[0]
            self.assertAlmostEqual(
                first[f"delta_g_theta{theta}"],
                first[f"g_pair_T{theta}"] - first["g_pair_T0"],
            )
            self.assertAlmostEqual(
                first[f"change_g_2_to_{theta}"],
                first[f"g_pair_T{theta}"] - first["g_pair_T2"],
            )
            self.assertTrue(
                all(row["Delta_beam"] == points[index]["Delta_beam"] for index, row in enumerate(points))
            )

    def test_sampled_status_first_exceedance_reentry_and_tolerance(self) -> None:
        statuses, label, first, reentry, nondecreasing = screening._sampled_status(
            [0.09, 0.11, 0.08, 0.12]
        )
        self.assertEqual(statuses, ["W", "E", "W", "E"])
        self.assertEqual(label, "FIRST_EXCEEDS_AT_3_DEG")
        self.assertEqual(first, 3.0)
        self.assertTrue(reentry)
        self.assertFalse(nondecreasing)
        self.assertTrue(
            screening._sampled_status([0.1, 0.1 - 5.0e-13, 0.1, 0.1])[4]
        )

    def test_sampled_combined_tables_and_accepted_endpoint_counts(self) -> None:
        sources = self._synthetic_sampled_sources()
        points = screening.build_sampled_point_comparison(sources)
        geometries = screening.build_sampled_geometry_comparison(sources)
        counts = screening.build_sampled_classification_counts(points, geometries)
        pairs = screening.build_sampled_pairwise_metrics(sources)
        self.assertEqual(len(points), 171)
        self.assertEqual(len(geometries), 9)
        self.assertEqual(len(counts), 4)
        self.assertEqual(len(pairs), 1026)
        self.assertEqual(counts[0]["pointwise_within_count"], 171)
        self.assertEqual(counts[0]["pointwise_exceeding_count"], 0)
        self.assertEqual(counts[0]["uniform_families_within_count"], 9)
        self.assertEqual(counts[3]["pointwise_within_count"], 66)
        self.assertEqual(counts[3]["pointwise_exceeding_count"], 105)
        self.assertEqual(counts[3]["uniform_families_within_count"], 0)

    def test_sampled_pairwise_same_pair_interval_identities(self) -> None:
        pairs = screening.build_sampled_pairwise_metrics(
            self._synthetic_sampled_sources()
        )
        row = pairs[137]
        for theta in (3, 4):
            self.assertAlmostEqual(
                row[f"delta_g_theta{theta}"],
                row[f"g_pair_T{theta}"] - row["g_pair_T0"],
            )
        for left, right in ((2, 3), (3, 4), (4, 5)):
            self.assertAlmostEqual(
                row[f"change_g_{left}_to_{right}"],
                row[f"g_pair_T{right}"] - row[f"g_pair_T{left}"],
            )
        self.assertEqual(row["pair_upper_mode"], row["pair_lower_mode"] + 1)

    def test_ap2_classification_decomposition_is_complete(self) -> None:
        for theta in (3, 4):
            spectra = self._synthetic_ap2_spectra(theta)
            points, _ = screening.build_ap2_metrics(
                spectra,
                screening.screening_geometries(),
                screening.beta_grid(),
                theta,
            )
            result = screening._classification_decomposition(points, theta)
            self.assertEqual(result["A_plus_B"] + result["C_simplification_within"], 171)
            self.assertEqual(result["total"], 171)

    def test_ap2_read_paths_preserve_ap0_ap1_and_supervisor_hashes(self) -> None:
        before = screening.accepted_artifact_hashes()
        screening.load_ap2_supervisor_family(
            screening.make_geometry(0.0, 0.0),
            "T4",
            np.asarray([0.0, 90.0]),
        )
        screening._read_csv(screening.AP0_SPECTRA_PATH)
        screening._read_csv(
            screening.AP1_OUTPUT_DIR / screening.SPECTRA_FILENAME
        )
        self.assertEqual(screening.accepted_artifact_hashes(), before)


if __name__ == "__main__":
    unittest.main()
