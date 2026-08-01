import math
import importlib
import sys
import unittest
from pathlib import Path
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis.anisotropic_rods import (  # noqa: E402
    validate_yartsev_ch2_unequal_thickness as validation,
)
from scripts.lib.yartsev_ch2_coupled_rods import joint_matrix_book  # noqa: E402
from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
    cantilever_geometry,
    physical_state_transfer_matrix,
    state_matrix,
)
from scripts.lib.yartsev_ch2_rectangular_eb import eb_state_matrix  # noqa: E402
from scipy.linalg import expm  # noqa: E402


class YartsevChapter2UnequalThicknessUT0Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.points = validation._make_section_points()
        cls.omega = 2.0 * math.pi * validation.DIAGNOSTIC_FREQUENCY_HZ

    def test_section_points_have_expected_geometry(self) -> None:
        for case, a_m in validation.SECTION_A_M.items():
            point = self.points[case]
            geometry = point.geometry
            with self.subTest(case=case):
                self.assertEqual(geometry.a, a_m)
                self.assertEqual(geometry.b, validation.WIDTH_B_M)
                self.assertEqual(geometry.length, validation.LENGTH_M)
                self.assertEqual(geometry.area, geometry.a * geometry.b)
                self.assertEqual(
                    geometry.I_y, geometry.a**3 * geometry.b / 12.0
                )
                self.assertEqual(
                    geometry.I_p,
                    geometry.a
                    * geometry.b
                    * (geometry.a**2 + geometry.b**2)
                    / 12.0,
                )

    def test_a5_matches_existing_cantilever_geometry_baseline(self) -> None:
        self.assertEqual(
            self.points["a5"].geometry, cantilever_geometry(0.2)
        )
        self.assertLessEqual(
            validation._baseline_a5_residual(self.points["a5"]), 1.0e-12
        )

    def test_theta_zero_torsional_identity(self) -> None:
        for case, point in self.points.items():
            with self.subTest(case=case):
                self.assertEqual(point.properties.Sbar16, 0.0j)
                self.assertLessEqual(
                    validation._relative_error(
                        point.torsion.C_T, point.torsion.Cbar
                    ),
                    1.0e-12,
                )

    def test_timoshenko_state_coefficients_are_arm_specific(self) -> None:
        area_coefficients = []
        bending_coefficients = []
        for case, point in self.points.items():
            matrix = state_matrix(self.omega, point)
            expected = validation._expected_timoshenko_state_matrix(
                self.omega, point
            )
            with self.subTest(case=case):
                np.testing.assert_array_equal(matrix, expected)
                self.assertLessEqual(
                    validation._timoshenko_state_coefficient_residual(
                        self.omega, point
                    ),
                    1.0e-12,
                )
            area_coefficients.append(matrix[3, 0])
            bending_coefficients.append(matrix[1, 4])
        self.assertEqual(len(set(area_coefficients)), 3)
        self.assertEqual(len(set(bending_coefficients)), 3)

    def test_eb_state_coefficients_are_arm_specific(self) -> None:
        area_coefficients = []
        bending_coefficients = []
        for case, point in self.points.items():
            matrix = eb_state_matrix(self.omega, point)
            expected = validation._expected_eb_state_matrix(self.omega, point)
            with self.subTest(case=case):
                np.testing.assert_array_equal(matrix, expected)
                self.assertLessEqual(
                    validation._eb_state_coefficient_residual(
                        self.omega, point
                    ),
                    1.0e-12,
                )
            area_coefficients.append(matrix[3, 0])
            bending_coefficients.append(matrix[1, 4])
        self.assertEqual(len(set(area_coefficients)), 3)
        self.assertEqual(len(set(bending_coefficients)), 3)

    def test_physical_transfer_scaling_is_arm_specific(self) -> None:
        for case, point in self.points.items():
            with self.subTest(case=case):
                direct = expm(state_matrix(self.omega, point) * 0.2)
                recovered = physical_state_transfer_matrix(self.omega, point)
                residual = np.linalg.norm(
                    recovered - direct, ord="fro"
                ) / np.linalg.norm(direct, ord="fro")
                self.assertLessEqual(residual, 1.0e-10)
                self.assertLessEqual(
                    validation._physical_transfer_scaling_residual(
                        self.omega, point
                    ),
                    1.0e-10,
                )

    def test_beta_zero_state_map_satisfies_joint_matrix(self) -> None:
        states = (
            np.arange(1.0, 7.0),
            np.array([0.5, -2.0, 3.5, -4.0, 5.25, -6.5]),
            np.array([7.0, 0.0, -1.0, 2.0, 0.0, 9.0]),
        )
        state_map = validation._beta0_state_map()
        for state_1 in states:
            two_arm_state = np.concatenate((state_1, state_map @ state_1))
            np.testing.assert_array_equal(
                joint_matrix_book(0.0) @ two_arm_state, np.zeros(6)
            )

    def test_stepped_builders_are_finite_3_by_3_and_independent(self) -> None:
        point_1, point_2 = self.points["a4"], self.points["a6"]
        with mock.patch.object(
            validation,
            "coupled_boundary_matrix",
            side_effect=AssertionError("coupled builder called"),
        ), mock.patch.object(
            validation,
            "coupled_boundary_matrix_raw",
            side_effect=AssertionError("coupled raw builder called"),
        ), mock.patch.object(
            validation,
            "eb_coupled_boundary_matrix",
            side_effect=AssertionError("EB coupled builder called"),
        ), mock.patch.object(
            validation,
            "eb_coupled_boundary_matrix_raw",
            side_effect=AssertionError("EB coupled raw builder called"),
        ):
            matrices = (
                validation._timoshenko_stepped_boundary_matrix_raw(
                    self.omega, point_1, point_2
                ),
                validation._timoshenko_stepped_boundary_matrix(
                    self.omega, point_1, point_2
                ),
                validation._eb_stepped_boundary_matrix_raw(
                    self.omega, point_1, point_2
                ),
                validation._eb_stepped_boundary_matrix(
                    self.omega, point_1, point_2
                ),
            )
        for matrix in matrices:
            self.assertEqual(matrix.shape, (3, 3))
            self.assertTrue(np.all(np.isfinite(matrix)))

    def test_analysis_module_import_has_not_run_smoke(self) -> None:
        reloaded = importlib.reload(validation)
        self.assertTrue(callable(reloaded.main))


if __name__ == "__main__":
    unittest.main()
