import math
import importlib
import sys
import tempfile
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
from scripts.lib.yartsev_ch2_coupled_rods import (  # noqa: E402
    straight_boundary_matrix_raw,
)
from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
    cantilever_geometry,
    physical_state_transfer_matrix,
    state_matrix,
)
from scripts.lib.yartsev_ch2_rectangular_eb import (  # noqa: E402
    eb_state_matrix,
    eb_straight_boundary_matrix_raw,
)
from scipy.linalg import expm  # noqa: E402


class YartsevChapter2UnequalThicknessUT0Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.points = validation._make_section_points()
        cls.omega = 2.0 * math.pi * validation.DIAGNOSTIC_FREQUENCY_HZ

    def _small_ut1a_pair(self, elements: int = 2):
        point_4 = validation._make_section_point(0.004)
        point_6 = validation._make_section_point(0.006)
        return (
            validation.assemble_two_arm_eb_fem(
                point_4, point_6, validation.BETA_RAD, elements
            ),
            validation.assemble_two_arm_eb_fem(
                point_6, point_4, validation.BETA_RAD, elements
            ),
        )

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

    def test_cli_modes_are_mutually_exclusive_and_no_mode_is_rejected(self) -> None:
        combinations = (
            ["--smoke", "--ut1-beta0"],
            ["--smoke", "--ut1a-fem-exchange-audit"],
            ["--ut1-beta0", "--ut1a-fem-exchange-audit"],
        )
        for arguments in combinations:
            with self.subTest(arguments=arguments), self.assertRaises(
                SystemExit
            ) as raised:
                validation.parse_args(arguments)
            self.assertEqual(raised.exception.code, 2)
        with mock.patch.object(validation, "_run_smoke") as smoke, mock.patch.object(
            validation, "_run_ut1_beta0"
        ) as ut1, mock.patch.object(
            validation, "_run_ut1a_fem_exchange_audit"
        ) as ut1a:
            self.assertEqual(validation.main([]), 2)
            smoke.assert_not_called()
            ut1.assert_not_called()
            ut1a.assert_not_called()

    def test_ut1_section_cases_use_two_independent_points(self) -> None:
        cases = validation._make_ut1_section_cases()
        self.assertEqual(set(cases), set(validation.UT1_SECTION_CASES))
        for case, (point_1, point_2) in cases.items():
            expected_a_1, expected_a_2 = validation.UT1_SECTION_CASES[case]
            with self.subTest(case=case):
                self.assertIsNot(point_1, point_2)
                self.assertIsNot(point_1.geometry, point_2.geometry)
                self.assertEqual(point_1.geometry.a, expected_a_1)
                self.assertEqual(point_2.geometry.a, expected_a_2)

    def test_seven_root_and_fem_mesh_configuration_is_ut1_only(self) -> None:
        self.assertEqual(validation.NUM_ROOTS, 3)
        self.assertEqual(validation.UT1_NUM_ROOTS, 7)
        self.assertEqual(validation.UT1_FEM_MESHES, (16, 32, 64))
        self.assertEqual(validation.UT1A_FEM_ELEMENTS_PER_ARM, 64)
        self.assertEqual(validation.BETA_RAD, 0.0)

    def test_ut1a_full_permutation_swaps_complete_arm_blocks(self) -> None:
        assembly_46, assembly_64 = self._small_ut1a_pair()
        permutation = validation._ut1a_full_swap_permutation(
            assembly_46, assembly_64
        )
        arm_size = 3 * (assembly_46.element_counts[0] + 1)
        self.assertEqual(permutation.shape, (2 * arm_size, 2 * arm_size))
        np.testing.assert_array_equal(permutation.T @ permutation, np.eye(2 * arm_size))
        np.testing.assert_array_equal(permutation @ permutation, np.eye(2 * arm_size))
        vector = np.arange(2 * arm_size, dtype=float)
        np.testing.assert_array_equal(
            permutation @ vector,
            np.concatenate((vector[arm_size:], vector[:arm_size])),
        )

    def test_ut1a_reduced_permutation_has_prescribed_joint_signs(self) -> None:
        assembly_46, assembly_64 = self._small_ut1a_pair()
        permutation = validation._ut1a_reduced_swap_permutation(
            assembly_46, assembly_64
        )
        internal_size = 3 * (assembly_46.element_counts[0] - 1)
        reduced_size = 2 * internal_size + 3
        self.assertEqual(permutation.shape, (reduced_size, reduced_size))
        np.testing.assert_array_equal(permutation.T @ permutation, np.eye(reduced_size))
        np.testing.assert_array_equal(permutation @ permutation, np.eye(reduced_size))
        vector = np.arange(1.0, reduced_size + 1.0)
        expected = np.concatenate(
            (
                vector[internal_size : 2 * internal_size],
                vector[:internal_size],
                [vector[-3], -vector[-2], -vector[-1]],
            )
        )
        np.testing.assert_array_equal(permutation @ vector, expected)

    def test_ut1a_reduction_intertwining_identity_on_small_mesh(self) -> None:
        assembly_46, assembly_64 = self._small_ut1a_pair()
        full = validation._ut1a_full_swap_permutation(assembly_46, assembly_64)
        reduced = validation._ut1a_reduced_swap_permutation(
            assembly_46, assembly_64
        )
        np.testing.assert_array_equal(
            assembly_64.reduction @ reduced,
            full @ assembly_46.reduction,
        )

    def test_ut1a_full_matrix_exchange_on_small_mesh(self) -> None:
        assembly_46, assembly_64 = self._small_ut1a_pair()
        full = validation._ut1a_full_swap_permutation(assembly_46, assembly_64)
        for left, source in (
            (assembly_64.stiffness_full, assembly_46.stiffness_full),
            (assembly_64.mass_full, assembly_46.mass_full),
        ):
            np.testing.assert_array_equal(left, full @ source @ full.T)

    def test_ut1a_reduced_matrix_congruence_on_small_mesh(self) -> None:
        assembly_46, assembly_64 = self._small_ut1a_pair()
        reduced = validation._ut1a_reduced_swap_permutation(
            assembly_46, assembly_64
        )
        for left, source in (
            (assembly_64.stiffness, assembly_46.stiffness),
            (assembly_64.mass, assembly_46.mass),
        ):
            np.testing.assert_array_equal(left, reduced @ source @ reduced.T)
            np.testing.assert_array_equal(source, reduced.T @ left @ reduced)

    def test_ut1a_transported_rayleigh_identity_for_fixed_vector(self) -> None:
        assembly_46, assembly_64 = self._small_ut1a_pair()
        reduced = validation._ut1a_reduced_swap_permutation(
            assembly_46, assembly_64
        )
        vector_46 = np.linspace(0.25, 2.25, assembly_46.stiffness.shape[0])
        vector_64 = reduced @ vector_46
        source = validation._ut1a_transported_rayleigh(
            assembly_46.stiffness, assembly_46.mass, vector_46
        )
        transported = validation._ut1a_transported_rayleigh(
            assembly_64.stiffness, assembly_64.mass, vector_64
        )
        self.assertLessEqual(abs(source - transported) / abs(source), 1.0e-13)

    def test_ut1a_pair_uses_exactly_two_mesh64_production_assemblies(self) -> None:
        production = validation.assemble_two_arm_eb_fem
        with mock.patch.object(
            validation, "assemble_two_arm_eb_fem", wraps=production
        ) as assembler, mock.patch.object(
            validation,
            "solve_two_arm_eb_fem",
            side_effect=AssertionError("mesh-64 eigensolve entered unit test"),
        ) as solver, mock.patch.object(
            validation,
            "find_elastic_roots",
            side_effect=AssertionError("continuum root finder entered UT-1a"),
        ) as root_finder:
            pair = validation._assemble_ut1a_fem_pair()
        self.assertEqual(len(pair), 2)
        self.assertEqual(assembler.call_count, 2)
        for call in assembler.call_args_list:
            self.assertEqual(call.args[3], 64)
            self.assertEqual(call.args[2], 0.0)
        solver.assert_not_called()
        root_finder.assert_not_called()

    def test_ut1a_driver_does_not_call_continuum_root_finder(self) -> None:
        small_pair = self._small_ut1a_pair()
        prior = {
            "fem_exchange_maximum": 5.762284502890923e-08,
            "ut0_regression_status": "PASS",
            "ut1_status": "PARTIAL_PASS",
        }
        with tempfile.TemporaryDirectory(dir=ROOT) as temporary, mock.patch.object(
            validation, "_assemble_ut1a_fem_pair", return_value=small_pair
        ) as pair_builder, mock.patch.object(
            validation, "_load_ut1_regression_evidence", return_value=prior
        ), mock.patch.object(
            validation,
            "find_elastic_roots",
            side_effect=AssertionError("continuum root finder entered UT-1a"),
        ) as root_finder:
            summary = validation._run_ut1a_fem_exchange_audit(Path(temporary))
            generated = {path.name for path in Path(temporary).iterdir()}
        pair_builder.assert_called_once_with()
        root_finder.assert_not_called()
        self.assertEqual(summary["continuum_root_calculations"], 0)
        self.assertEqual(
            generated,
            {
                "matrix_congruence.csv",
                "eigenpair_transport.csv",
                "ut1a_summary.json",
                "ut1a_report.md",
            },
        )

    def test_ut1a_diagnostic_helpers_do_not_modify_production_module(self) -> None:
        library = ROOT / "scripts" / "lib" / "yartsev_ch2_rectangular_eb.py"
        before = library.read_bytes()
        assembly_46, assembly_64 = self._small_ut1a_pair()
        full = validation._ut1a_full_swap_permutation(assembly_46, assembly_64)
        reduced = validation._ut1a_reduced_swap_permutation(
            assembly_46, assembly_64
        )
        rows = validation._ut1a_matrix_congruence_rows(
            assembly_46, assembly_64, full, reduced
        )
        self.assertTrue(all(row["status"] == "PASS" for row in rows))
        self.assertEqual(library.read_bytes(), before)

    def test_equal_section_stepped_matches_homogeneous_straight_matrices(self) -> None:
        arm = self.points["a5"]
        straight = validation._make_section_point(0.005, length_m=0.4)
        for frequency_hz in (100.0, 500.0, 1000.0):
            omega = 2.0 * math.pi * frequency_hz
            pairs = (
                (
                    validation._timoshenko_stepped_boundary_matrix_raw(
                        omega, arm, arm
                    ),
                    straight_boundary_matrix_raw(omega, straight),
                ),
                (
                    validation._eb_stepped_boundary_matrix_raw(
                        omega, arm, arm
                    ),
                    eb_straight_boundary_matrix_raw(omega, straight),
                ),
            )
            for stepped, homogeneous in pairs:
                determinant_residual = abs(
                    np.linalg.det(stepped) - np.linalg.det(homogeneous)
                ) / max(abs(np.linalg.det(homogeneous)), np.finfo(float).tiny)
                self.assertLessEqual(determinant_residual, 1.0e-12)
                stepped_scaled = validation.equilibrate_matrix(stepped)[0]
                homogeneous_scaled = validation.equilibrate_matrix(homogeneous)[0]
                np.testing.assert_allclose(
                    np.linalg.svd(stepped_scaled, compute_uv=False),
                    np.linalg.svd(homogeneous_scaled, compute_uv=False),
                    rtol=1.0e-12,
                    atol=1.0e-14,
                )

    def test_both_unequal_stepped_orders_are_finite_3_by_3(self) -> None:
        for point_1, point_2 in (
            (self.points["a4"], self.points["a6"]),
            (self.points["a6"], self.points["a4"]),
        ):
            for builder in (
                validation._timoshenko_stepped_boundary_matrix_raw,
                validation._eb_stepped_boundary_matrix_raw,
            ):
                matrix = builder(self.omega, point_1, point_2)
                self.assertEqual(matrix.shape, (3, 3))
                self.assertTrue(np.all(np.isfinite(matrix)))

    def test_cluster_matching_handles_separated_and_connected_clusters(self) -> None:
        separated = validation._match_frequency_spectra(
            [100.0, 200.0, 300.0], [100.1, 200.1, 300.1]
        )
        np.testing.assert_allclose(separated[0], [100.1, 200.1, 300.1])
        self.assertEqual(separated[1], ["sorted", "sorted", "sorted"])

        two = validation._match_frequency_spectra(
            [100.0, 100.05, 200.0], [100.04, 100.01, 200.1]
        )
        np.testing.assert_allclose(two[0], [100.01, 100.04, 200.1])
        self.assertTrue(two[1][0].startswith("local_frequency_cluster_modes_1-2"))
        self.assertEqual(two[1][2], "sorted")

        three = validation._match_frequency_spectra(
            [100.0, 100.04, 100.08], [100.07, 100.01, 100.05]
        )
        np.testing.assert_allclose(three[0], [100.01, 100.05, 100.07])
        self.assertTrue(
            all(
                item.startswith("local_frequency_cluster_modes_1-2-3")
                for item in three[1]
            )
        )

    def test_fem_runner_delegates_to_existing_solver(self) -> None:
        point_1, point_2 = self.points["a4"], self.points["a6"]
        existing_solver = validation.solve_two_arm_eb_fem
        with mock.patch.object(
            validation, "solve_two_arm_eb_fem", wraps=existing_solver
        ) as solver:
            result = validation._run_eb_fem_case(
                section_case="unit",
                point_1=point_1,
                point_2=point_2,
                elements_per_arm=4,
                analytic_frequencies=np.linspace(100.0, 700.0, 7),
            )
        solver.assert_called_once()
        self.assertEqual(len(result["solution"].frequencies_hz), 7)


if __name__ == "__main__":
    unittest.main()
