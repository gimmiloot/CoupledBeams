import math
import importlib
import inspect
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

    def _small_ut2_assembly(
        self, a_1: float, a_2: float, beta_deg: float, elements: int = 2
    ):
        return validation.assemble_two_arm_eb_fem(
            validation._make_section_point(a_1),
            validation._make_section_point(a_2),
            math.radians(beta_deg),
            elements,
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
            ["--smoke", "--ut2-beta30"],
            ["--ut1-beta0", "--ut2-beta30"],
            ["--ut1a-fem-exchange-audit", "--ut2-beta30"],
            ["--smoke", "--ut3-beta90"],
            ["--ut1-beta0", "--ut3-beta90"],
            ["--ut1a-fem-exchange-audit", "--ut3-beta90"],
            ["--ut2-beta30", "--ut3-beta90"],
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
        ) as ut1a, mock.patch.object(
            validation, "_run_ut2_beta30"
        ) as ut2, mock.patch.object(
            validation, "_run_ut3_beta90"
        ) as ut3:
            self.assertEqual(validation.main([]), 2)
            smoke.assert_not_called()
            ut1.assert_not_called()
            ut1a.assert_not_called()
            ut2.assert_not_called()
            ut3.assert_not_called()

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

    def test_ut2_scope_constants_are_exact(self) -> None:
        self.assertEqual(validation.UT2_ANGLES_DEG, (-30.0, 30.0))
        self.assertEqual(
            validation.UT2_SECTION_CASES,
            {
                "baseline_5_5": (0.005, 0.005),
                "asymmetric_4_6": (0.004, 0.006),
                "swapped_6_4": (0.006, 0.004),
            },
        )
        self.assertEqual(validation.UT2_NUM_ROOTS, 7)
        self.assertEqual(validation.UT2_FEM_ELEMENTS_PER_ARM, (64, 64))

    def test_ut3_scope_constants_are_exact(self) -> None:
        self.assertEqual(validation.UT3_ANGLES_DEG, (-90.0, 90.0))
        self.assertEqual(
            validation.UT3_SECTION_CASES,
            {
                "baseline_5_5": (0.005, 0.005),
                "asymmetric_4_6": (0.004, 0.006),
                "swapped_6_4": (0.006, 0.004),
            },
        )
        self.assertEqual(validation.UT3_NUM_ROOTS, 7)
        self.assertEqual(validation.UT3_FEM_ELEMENTS_PER_ARM, (64, 64))

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

    def test_ut2_reflection_full_transform_flips_only_phi(self) -> None:
        assembly = self._small_ut2_assembly(0.004, 0.006, 30.0)
        transform = validation._ut2_reflection_full_transform(assembly)
        full_size = assembly.stiffness_full.shape[0]
        self.assertEqual(transform.shape, (full_size, full_size))
        np.testing.assert_array_equal(transform.T @ transform, np.eye(full_size))
        np.testing.assert_array_equal(transform @ transform, np.eye(full_size))
        vector = np.arange(1.0, full_size + 1.0)
        expected = vector * np.tile([1.0, 1.0, -1.0], full_size // 3)
        np.testing.assert_array_equal(transform @ vector, expected)

    def test_ut2_reflection_reduced_transform_has_declared_signs(self) -> None:
        assembly = self._small_ut2_assembly(0.004, 0.006, 30.0)
        transform = validation._ut2_reflection_reduced_transform(assembly)
        reduced_size = assembly.stiffness.shape[0]
        self.assertEqual(transform.shape, (reduced_size, reduced_size))
        np.testing.assert_array_equal(
            transform.T @ transform, np.eye(reduced_size)
        )
        np.testing.assert_array_equal(transform @ transform, np.eye(reduced_size))
        vector = np.arange(1.0, reduced_size + 1.0)
        internal_nodes = (reduced_size - 3) // 3
        signs = np.concatenate(
            (
                np.tile([1.0, 1.0, -1.0], internal_nodes),
                [1.0, -1.0, 1.0],
            )
        )
        np.testing.assert_array_equal(transform @ vector, signs * vector)

    def test_ut2_reflection_endpoint_map_identities(self) -> None:
        rows = [
            row
            for row in validation._ut2_endpoint_map_rows()
            if row["symmetry_type"] == "reflection"
        ]
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row["status"] == "PASS" for row in rows))
        self.assertLessEqual(
            max(row["relative_max_residual"] for row in rows), 1.0e-13
        )

    def test_ut2_relabel_full_transform_swaps_without_local_signs(self) -> None:
        source = self._small_ut2_assembly(0.004, 0.006, 30.0)
        target = self._small_ut2_assembly(0.006, 0.004, -30.0)
        transform = validation._ut2_relabel_full_transform(source, target)
        arm_size = source.stiffness_full.shape[0] // 2
        vector = np.arange(source.stiffness_full.shape[0], dtype=float)
        np.testing.assert_array_equal(
            transform @ vector,
            np.concatenate((vector[arm_size:], vector[:arm_size])),
        )
        np.testing.assert_array_equal(transform.T @ transform, np.eye(len(vector)))

    def test_ut2_relabel_joint_transform_is_declared_matrix(self) -> None:
        beta = math.radians(30.0)
        expected = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, -math.cos(beta), math.sin(beta)],
                [0.0, -math.sin(beta), -math.cos(beta)],
            ]
        )
        np.testing.assert_array_equal(
            validation._ut2_relabel_joint_transform(beta), expected
        )

    def test_ut2_relabel_endpoint_map_identities(self) -> None:
        rows = [
            row
            for row in validation._ut2_endpoint_map_rows()
            if row["symmetry_type"] == "oriented_angle_relabeling"
        ]
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row["status"] == "PASS" for row in rows))
        self.assertLessEqual(
            max(row["relative_max_residual"] for row in rows), 1.0e-13
        )

    def test_ut2_reflection_matrix_identities_on_small_mesh(self) -> None:
        plus = self._small_ut2_assembly(0.004, 0.006, 30.0)
        minus = self._small_ut2_assembly(0.004, 0.006, -30.0)
        rows = validation._ut2_reflection_matrix_rows(
            "asymmetric_4_6", plus, minus
        )
        self.assertEqual(len(rows), 9)
        self.assertTrue(all(row["status"] == "PASS" for row in rows))
        self.assertTrue(
            any(row["check"] == "reflection reduction intertwining" for row in rows)
        )
        self.assertTrue(
            any(row["check"] == "K_full reflection congruence" for row in rows)
        )
        self.assertTrue(
            any(row["check"] == "M_reduced reflection congruence" for row in rows)
        )

    def test_ut2_relabel_matrix_identities_on_small_mesh(self) -> None:
        source = self._small_ut2_assembly(0.004, 0.006, 30.0)
        target = self._small_ut2_assembly(0.006, 0.004, -30.0)
        rows = validation._ut2_relabeling_matrix_rows(
            "asymmetric_4_6", "swapped_6_4", source, target
        )
        self.assertEqual(len(rows), 7)
        self.assertTrue(all(row["status"] == "PASS" for row in rows))
        self.assertTrue(
            any(row["check"] == "relabeling reduction intertwining" for row in rows)
        )
        self.assertTrue(
            any(row["check"] == "K_full relabeling congruence" for row in rows)
        )
        self.assertTrue(
            any(row["check"] == "M_reduced relabeling congruence" for row in rows)
        )

    def test_ut2_continuum_factories_do_not_call_stepped_references(self) -> None:
        point_1, point_2 = self.points["a4"], self.points["a6"]
        forbidden = (
            "_timoshenko_stepped_boundary_matrix",
            "_timoshenko_stepped_boundary_matrix_raw",
            "_eb_stepped_boundary_matrix",
            "_eb_stepped_boundary_matrix_raw",
        )
        patches = [
            mock.patch.object(
                validation,
                name,
                side_effect=AssertionError("stepped reference entered UT-2"),
            )
            for name in forbidden
        ]
        with patches[0], patches[1], patches[2], patches[3]:
            for model in ("Timoshenko", "EB"):
                factory, raw = validation._ut2_continuum_factories(
                    model, math.radians(30.0), point_1, point_2
                )
                self.assertEqual(factory(self.omega).shape, (6, 6))
                self.assertEqual(raw(self.omega).shape, (6, 6))

    def test_ut2_driver_contains_no_refinement_or_ut1a_eigenpair_audit(self) -> None:
        source = inspect.getsource(validation._run_ut2_beta30)
        for forbidden in (
            "UT1_FEM_MESHES",
            "_run_eb_fem_case",
            "_timoshenko_stepped",
            "_eb_stepped",
            "_ut1a_native_and_transport_audit",
            "_ut1a_transported_rayleigh",
        ):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, source)

    def test_ut2_native_fem_spectral_symmetry_is_not_a_hard_gate(self) -> None:
        parameters = inspect.signature(validation._classify_ut2_status).parameters
        self.assertNotIn("native_fem_spectral_symmetry_ok", parameters)
        self.assertEqual(
            validation._classify_ut2_status(
                continuum_hard_ok=True,
                matrix_hard_ok=True,
                fem_hard_ok=True,
                regressions_ok=True,
            ),
            "PASS",
        )
        self.assertEqual(
            validation._classify_ut2_status(
                continuum_hard_ok=True,
                matrix_hard_ok=True,
                fem_hard_ok=False,
                regressions_ok=True,
            ),
            "PARTIAL_PASS",
        )

    def test_ut3_exact_joint_matrices_have_prescribed_integer_entries(self) -> None:
        plus = np.array(
            [
                [1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1],
            ],
            dtype=float,
        )
        minus = np.array(
            [
                [1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0],
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1],
            ],
            dtype=float,
        )
        np.testing.assert_array_equal(validation._ut3_joint_plus_90_exact(), plus)
        np.testing.assert_array_equal(validation._ut3_joint_minus_90_exact(), minus)

    def test_ut3_joint_book_and_endpoint_exact_limits(self) -> None:
        rows = validation._ut3_exact_joint_limit_rows()
        selected = {
            row["check"]: row
            for row in rows
            if row["check"].startswith("J_book")
            or row["check"]
            in ("F2(+90) vs I", "F2(-90) vs diag(1,-1,-1)")
        }
        self.assertEqual(
            set(selected),
            {
                "J_book(+90) vs J_plus_90_exact",
                "J_book(-90) vs J_minus_90_exact",
                "F2(+90) vs I",
                "F2(-90) vs diag(1,-1,-1)",
            },
        )
        self.assertTrue(all(row["status"] == "PASS" for row in selected.values()))
        self.assertLessEqual(
            max(row["absolute_max_residual"] for row in selected.values()),
            1.0e-14,
        )

    def test_ut3_exact_continuum_builders_are_independent(self) -> None:
        point_1, point_2 = self.points["a4"], self.points["a6"]
        with mock.patch.object(
            validation,
            "joint_matrix_book",
            side_effect=AssertionError("J_book entered exact builder"),
        ), mock.patch.object(
            validation,
            "coupled_boundary_matrix",
            side_effect=AssertionError("coupled builder entered exact builder"),
        ), mock.patch.object(
            validation,
            "eb_coupled_boundary_matrix",
            side_effect=AssertionError("EB coupled builder entered exact builder"),
        ):
            for model in ("Timoshenko", "EB"):
                scaled, raw = validation._ut3_exact_continuum_factories(
                    model, point_1, point_2
                )
                self.assertEqual(scaled(self.omega).shape, (6, 6))
                self.assertEqual(raw(self.omega).shape, (6, 6))

    def test_ut3_deterministic_channel_exchange_relations(self) -> None:
        rows = [
            row
            for row in validation._ut3_exact_joint_limit_rows()
            if row["check"].startswith("deterministic channel-exchange")
        ]
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(row["status"] == "PASS" for row in rows))
        self.assertLessEqual(
            max(row["absolute_max_residual"] for row in rows), 1.0e-14
        )

    def test_ut3_exact_relabel_joint_transform_is_not_involutory(self) -> None:
        exact = np.array(
            [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]]
        )
        computed = validation._ut2_relabel_joint_transform(math.pi / 2.0)
        np.testing.assert_array_equal(validation._ut3_relabel_joint_exact(), exact)
        np.testing.assert_allclose(computed, exact, rtol=0.0, atol=1.0e-14)
        np.testing.assert_array_equal(exact.T @ exact, np.eye(3))
        self.assertFalse(np.array_equal(exact @ exact, np.eye(3)))
        np.testing.assert_array_equal(np.linalg.inv(exact), exact.T)

    def test_ut3_reflection_and_relabel_endpoint_identities(self) -> None:
        rows = validation._ut2_endpoint_map_rows(90.0)
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(row["status"] == "PASS" for row in rows))
        exact_rows = validation._ut3_exact_joint_limit_rows()
        endpoint = [
            row
            for row in exact_rows
            if row["check"].startswith(("F1 H_rel", "F2(-90) H_rel"))
        ]
        self.assertEqual(len(endpoint), 2)
        self.assertTrue(all(row["status"] == "PASS" for row in endpoint))

    def test_ut3_reflection_matrix_identities_on_small_mesh(self) -> None:
        plus = self._small_ut2_assembly(0.004, 0.006, 90.0)
        minus = self._small_ut2_assembly(0.004, 0.006, -90.0)
        rows = validation._ut2_reflection_matrix_rows(
            "asymmetric_4_6", plus, minus, beta_deg=90.0
        )
        self.assertEqual(len(rows), 9)
        self.assertTrue(all(row["status"] == "PASS" for row in rows))

    def test_ut3_relabel_matrix_identities_on_small_mesh(self) -> None:
        source = self._small_ut2_assembly(0.004, 0.006, 90.0)
        target = self._small_ut2_assembly(0.006, 0.004, -90.0)
        rows = validation._ut2_relabeling_matrix_rows(
            "asymmetric_4_6",
            "swapped_6_4",
            source,
            target,
            beta_deg=90.0,
        )
        self.assertEqual(len(rows), 7)
        self.assertTrue(all(row["status"] == "PASS" for row in rows))

    def test_ut3_driver_excludes_refinement_timoshenko_fem_and_ut1a_audit(self) -> None:
        source = inspect.getsource(validation._run_ut3_beta90)
        for forbidden in (
            "UT1_FEM_MESHES",
            "_run_eb_fem_case",
            "_ut1a_native_and_transport_audit",
            "_ut1a_transported_rayleigh",
            "Timoshenko FEM",
        ):
            with self.subTest(forbidden=forbidden):
                self.assertNotIn(forbidden, source)
        self.assertIn("_run_ut2_eb_fem_configuration", source)
        self.assertGreaterEqual(
            source.count("for beta_deg in UT3_ANGLES_DEG"), 2
        )

    def test_ut3_native_fem_symmetry_is_diagnostic_only_and_phase_completes(self) -> None:
        parameters = inspect.signature(validation._classify_ut3_status).parameters
        self.assertNotIn("native_fem_spectral_symmetry_ok", parameters)
        self.assertEqual(
            validation._classify_ut3_status(
                continuum_hard_ok=True,
                exact_limit_ok=True,
                matrix_hard_ok=True,
                fem_hard_ok=True,
                regressions_ok=True,
            ),
            ("PASS", "PARTIAL_PASS", "COMPLETE"),
        )
        self.assertEqual(
            validation._classify_ut3_status(
                continuum_hard_ok=True,
                exact_limit_ok=True,
                matrix_hard_ok=True,
                fem_hard_ok=False,
                regressions_ok=True,
            ),
            ("PARTIAL_PASS", "PARTIAL_PASS", "COMPLETE_WITH_LIMITATIONS"),
        )

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
