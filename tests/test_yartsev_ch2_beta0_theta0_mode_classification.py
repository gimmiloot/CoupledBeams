import contextlib
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


from scripts.analysis.anisotropic_rods import (  # noqa: E402
    classify_yartsev_ch2_beta0_theta0_modes as classification,
)
from scripts.lib.yartsev_ch2_coupled_rods import (  # noqa: E402
    coupled_boundary_matrix_raw,
)
from scripts.lib.yartsev_ch2_monoclinic_rod import RootResult  # noqa: E402


class YartsevChapter2Beta0Theta0ModeClassificationTest(unittest.TestCase):
    @staticmethod
    def _root(frequency_hz: float) -> RootResult:
        return RootResult(
            omega=complex(2.0 * math.pi * frequency_hz),
            frequency_hz=float(frequency_hz),
            determinant_residual=0.0,
            raw_determinant_abs=0.0,
            sigma_min=0.0,
            sigma_max=1.0,
            relative_singular_residual=0.0,
            status="accepted_synthetic",
        )

    @classmethod
    def _spectrum(cls, frequencies_hz: list[float]) -> classification.SpectrumResult:
        roots = tuple(cls._root(value) for value in frequencies_hz)
        quality = tuple(
            {
                "quality_status": "PASS",
                "quality_basis": "scaled",
                "root_status": "accepted_synthetic",
            }
            for _ in roots
        )
        return classification.SpectrumResult(
            roots=roots,
            quality=quality,
            matrix_evaluations=0,
            raw_quality_evaluations=0,
            runtime_seconds=0.0,
        )

    def test_exact_case_presets_and_article_geometries(self) -> None:
        self.assertEqual(
            [(case.case_id, case.mu, case.tau) for case in classification.CASES],
            [("A", 0.0, 0.0), ("B", 0.5, -0.2)],
        )
        self.assertEqual(classification.FORMULATION, "state_corrected")
        self.assertEqual(classification.CLAMP, "book_slope_clamp")
        self.assertEqual(classification.MATERIAL_MODE, "elastic")
        self.assertEqual(classification.BETA_DEG, 0.0)
        self.assertEqual(classification.THETA_DEG, 0.0)
        self.assertEqual(classification.FULL_ROOT_COUNT, 7)
        self.assertEqual(classification.PLOTTED_ROOT_COUNT, 6)
        self.assertEqual(classification.PARTIAL_ROOT_COUNT, 7)

        for case in classification.CASES:
            with self.subTest(case=case.case_id):
                denominator = 1.0 + case.tau**2 + 2.0 * case.mu * case.tau
                scale = classification.REFERENCE_A_M / math.sqrt(denominator)
                expected_a_1 = scale * (1.0 - case.tau)
                expected_a_2 = scale * (1.0 + case.tau)
                expected_l_1 = classification.REFERENCE_LENGTH_M * (1.0 - case.mu)
                expected_l_2 = classification.REFERENCE_LENGTH_M * (1.0 + case.mu)

                geometry_1, geometry_2 = classification.article_geometry(
                    case.mu, case.tau
                )
                np.testing.assert_allclose(
                    [
                        geometry_1.a,
                        geometry_1.b,
                        geometry_1.length,
                        geometry_2.a,
                        geometry_2.b,
                        geometry_2.length,
                    ],
                    [
                        expected_a_1,
                        4.0 * expected_a_1,
                        expected_l_1,
                        expected_a_2,
                        4.0 * expected_a_2,
                        expected_l_2,
                    ],
                    rtol=0.0,
                    atol=2.0e-18,
                )
                reference_volume = (
                    2.0
                    * classification.REFERENCE_A_M
                    * classification.REFERENCE_B_M
                    * classification.REFERENCE_LENGTH_M
                )
                self.assertAlmostEqual(
                    geometry_1.area * geometry_1.length
                    + geometry_2.area * geometry_2.length,
                    reference_volume,
                    delta=2.0e-20,
                )

                point_1, point_2 = classification.case_points(case)
                for point in (point_1, point_2):
                    self.assertEqual(point.material_mode, "elastic")
                    self.assertEqual(point.properties.theta_deg, 0.0)
                    self.assertEqual(point.properties.Sbar16, 0.0)
                    self.assertEqual(point.properties.Ex.real, 191.0e9)
                    self.assertEqual(point.material.rho, 1580.0)
                    self.assertEqual(
                        point.geometry.shear_factor, classification.SHEAR_FACTOR
                    )

    def test_fixed_reference_lambda_formula(self) -> None:
        omega = np.asarray([0.0, 2.0 * math.pi * 100.0, 2.0 * math.pi * 750.0])
        area_0 = classification.REFERENCE_A_M * classification.REFERENCE_B_M
        inertia_0 = (
            classification.REFERENCE_A_M**3
            * classification.REFERENCE_B_M
            / 12.0
        )
        expected = (
            1580.0
            * area_0
            * omega**2
            * classification.REFERENCE_LENGTH_M**4
            / (classification.REFERENCE_EX_PA * inertia_0)
        ) ** 0.25
        np.testing.assert_allclose(
            classification.lambda_reference(omega), expected, rtol=2.0e-15, atol=0.0
        )
        self.assertAlmostEqual(
            classification.lambda_reference(float(omega[1])),
            float(expected[1]),
            delta=2.0e-15,
        )

    def test_actual_beta0_matrices_have_the_declared_blocks(self) -> None:
        self.assertEqual(classification.BENDING_ROW_INDICES, (0, 2, 3, 5))
        self.assertEqual(classification.TORSION_ROW_INDICES, (1, 4))
        self.assertEqual(classification.BENDING_COLUMN_INDICES, (0, 1, 3, 4))
        self.assertEqual(classification.TORSION_COLUMN_INDICES, (2, 5))
        omega = 2.0 * math.pi * 500.0

        for case in classification.CASES:
            with self.subTest(case=case.case_id):
                point_1, point_2 = classification.case_points(case)
                raw = coupled_boundary_matrix_raw(omega, 0.0, point_1, point_2)
                permuted, bending, torsion = classification.split_raw_matrix(raw)
                self.assertEqual(permuted.shape, (6, 6))
                self.assertEqual(bending.shape, (4, 4))
                self.assertEqual(torsion.shape, (2, 2))
                np.testing.assert_array_equal(
                    bending,
                    raw[
                        np.ix_(
                            classification.BENDING_ROW_INDICES,
                            classification.BENDING_COLUMN_INDICES,
                        )
                    ],
                )
                np.testing.assert_array_equal(
                    torsion,
                    raw[
                        np.ix_(
                            classification.TORSION_ROW_INDICES,
                            classification.TORSION_COLUMN_INDICES,
                        )
                    ],
                )

                upper_right = permuted[:4, 4:]
                lower_left = permuted[4:, :4]
                off_block = np.concatenate(
                    (upper_right.ravel(), lower_left.ravel())
                )
                relative = float(np.linalg.norm(off_block)) / max(
                    float(np.linalg.norm(permuted)), np.finfo(float).tiny
                )
                self.assertLessEqual(relative, classification.BLOCK_RELATIVE_TOLERANCE)
                audit = classification.block_separation_row(
                    case, 500.0, raw
                )
                self.assertEqual(audit["status"], "PASS")
                self.assertLessEqual(
                    audit["off_block_relative_frobenius_norm"],
                    classification.BLOCK_RELATIVE_TOLERANCE,
                )

                factories = classification.matrix_factories(point_1, point_2)
                np.testing.assert_array_equal(factories["full"][0](omega), raw)
                np.testing.assert_array_equal(factories["bending"][0](omega), bending)
                np.testing.assert_array_equal(factories["torsion"][0](omega), torsion)

    def test_synthetic_partial_union_preserves_order_and_family_number(self) -> None:
        bending = self._spectrum([10.0, 30.0, 50.0, 70.0])
        torsion = self._spectrum([20.0, 40.0, 60.0, 80.0])
        merged = classification.merge_partial_spectra(bending, torsion)
        self.assertEqual(
            [item.root.frequency_hz for item in merged],
            [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0],
        )
        self.assertEqual(
            [(item.family, item.partial_mode_number) for item in merged],
            [
                ("bending", 1),
                ("torsion", 1),
                ("bending", 2),
                ("torsion", 2),
                ("bending", 3),
                ("torsion", 3),
                ("bending", 4),
                ("torsion", 4),
            ],
        )
        self.assertEqual(classification.cross_family_near_coincidences(merged), [])

        coincident = classification.merge_partial_spectra(
            self._spectrum([10.0, 30.0]), self._spectrum([20.0, 30.0])
        )
        self.assertEqual(
            [(item.family, item.partial_mode_number) for item in coincident],
            [("bending", 1), ("torsion", 1), ("bending", 2), ("torsion", 2)],
        )
        self.assertEqual(
            classification.cross_family_near_coincidences(coincident),
            [(2, 3, 0.0)],
        )

    def test_dimensionless_shape_norm_suppresses_the_opposite_exact_block(self) -> None:
        point_1, point_2 = classification.case_points(classification.CASES[0])
        omega = complex(2.0 * math.pi * 100.0)

        with mock.patch.object(classification, "SHAPE_SAMPLES", 7):
            pure_bending_matrix = np.eye(6, dtype=np.complex128)
            pure_bending_matrix[:, 0] = 0.0
            bending = classification.shape_block_norms(
                omega, point_1, point_2, pure_bending_matrix
            )
            self.assertAlmostEqual(
                bending["bending_block_relative_norm"], 1.0, delta=2.0e-14
            )
            self.assertLessEqual(bending["torsion_block_relative_norm"], 2.0e-14)
            self.assertLessEqual(
                bending["full_null_vector_relative_residual"], 2.0e-14
            )

            pure_torsion_matrix = np.eye(6, dtype=np.complex128)
            pure_torsion_matrix[:, 2] = 0.0
            torsion = classification.shape_block_norms(
                omega, point_1, point_2, pure_torsion_matrix
            )
            self.assertAlmostEqual(
                torsion["torsion_block_relative_norm"], 1.0, delta=2.0e-14
            )
            self.assertLessEqual(torsion["bending_block_relative_norm"], 2.0e-14)
            self.assertLessEqual(
                torsion["full_null_vector_relative_residual"], 2.0e-14
            )

    def test_reuse_and_force_recompute_are_mutually_exclusive(self) -> None:
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            classification.parse_args(["--reuse-data", "--force-recompute"])

    def test_reuse_data_validates_hashes_without_calling_root_solver(self) -> None:
        results_root = ROOT / "results"
        results_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=results_root) as temporary:
            output_dir = Path(temporary)
            classification_rows = [
                {
                    "case_id": case.case_id,
                    "sorted_position": mode,
                    "root_quality_status": "PASS",
                }
                for case in classification.CASES
                for mode in range(1, classification.PLOTTED_ROOT_COUNT + 1)
            ]
            scientific_rows = {
                "mode_classification.csv": classification_rows,
                "full_partial_spectrum_comparison.csv": [{"status": "PASS"}],
                "block_separation_audit.csv": [{"status": "PASS"}],
            }
            hashes: dict[str, str] = {}
            for filename, rows in scientific_rows.items():
                path = output_dir / filename
                classification.write_csv(path, rows)
                hashes[filename] = classification.sha256(path)

            manifest = {
                "status": "PASS",
                "scientific_csv_sha256": hashes,
                "gate_results": {
                    "raw_maxima": {
                        "max_spectral_union_relative_difference": 0.0,
                    }
                },
            }
            (output_dir / "manifest.json").write_text(
                json.dumps(manifest, ensure_ascii=False), encoding="utf-8"
            )

            with (
                mock.patch.object(
                    classification,
                    "find_elastic_roots",
                    side_effect=AssertionError("root solver must not run"),
                ) as root_solver,
                mock.patch.object(
                    classification,
                    "run_computation",
                    side_effect=AssertionError("scientific computation must not run"),
                ) as computation,
                contextlib.redirect_stdout(io.StringIO()),
            ):
                self.assertEqual(
                    classification.main(
                        ["--reuse-data", "--output-dir", str(output_dir)]
                    ),
                    0,
                )
            root_solver.assert_not_called()
            computation.assert_not_called()


if __name__ == "__main__":
    unittest.main()
