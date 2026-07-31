import math
import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib.yartsev_ch2_coupled_rods import (  # noqa: E402
    coupled_boundary_matrix,
    coupled_boundary_matrix_raw,
    equilibrate_matrix,
    joint_matrix_book,
    straight_boundary_matrix,
)
from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
    Geometry,
    cantilever_geometry,
    find_elastic_roots,
    hms_dx_209_material,
    make_rod_point,
)
from scripts.lib.yartsev_ch2_rectangular_eb import (  # noqa: E402
    assemble_two_arm_eb_fem,
    eb_bending_coupled_boundary_matrix,
    eb_clamp_matrix,
    eb_coupled_boundary_matrix,
    eb_coupled_boundary_matrix_raw,
    eb_element_matrices,
    eb_joint_dof_maps,
    eb_joint_mapping_residual,
    eb_physical_end_map,
    eb_state_matrix,
    eb_state_transfer_matrix,
    eb_straight_boundary_matrix,
    eb_straight_boundary_matrix_raw,
    eb_straight_right_clamp_matrix,
    fixed_fixed_bending_dimensionless_roots,
    fixed_fixed_bending_frequencies_hz,
    fixed_fixed_torsion_frequencies_hz,
    saint_venant_stiffness,
    solve_two_arm_eb_fem,
)


class YartsevChapter2RectangularEBTest(unittest.TestCase):
    _targeted_cache = {}

    @staticmethod
    def _point(length: float = 0.2, scale: float = 1.0):
        base = cantilever_geometry(length)
        geometry = Geometry(
            a=scale * base.a,
            b=scale * base.b,
            length=length,
            shear_factor=base.shear_factor,
        )
        return make_rod_point(
            0.0,
            geometry=geometry,
            material=hms_dx_209_material(),
        )

    @staticmethod
    def _roots(point, builder, count=7):
        def boundary(omega, _point, _formulation):
            return builder(omega)

        return find_elastic_roots(
            point,
            "state_corrected",
            num_roots=count,
            scan_step_hz=10.0,
            initial_max_hz=5000.0,
            max_hz=100_000.0,
            boundary_matrix_builder=boundary,
        )

    @classmethod
    def _targeted_solution(cls, element_counts):
        cached = cls._targeted_cache.get(element_counts)
        if cached is not None:
            return cached
        point_1, point_2 = cls._point(0.1), cls._point(0.3)
        assembly = assemble_two_arm_eb_fem(
            point_1, point_2, 0.0, element_counts
        )
        solution = solve_two_arm_eb_fem(assembly, num_roots=7)
        straight = cls._point(0.4)
        exact = np.sort(
            np.concatenate(
                [
                    fixed_fixed_bending_frequencies_hz(straight, 7),
                    fixed_fixed_torsion_frequencies_hz(straight, 7),
                ]
            )
        )[:7]
        errors = np.abs(solution.frequencies_hz - exact) / exact
        cached = (assembly, solution, exact, errors)
        cls._targeted_cache[element_counts] = cached
        return cached

    def test_a_theta_zero_is_orthotropic(self) -> None:
        point = self._point()
        self.assertEqual(point.properties.Sbar16, 0.0j)

    def test_b_theta_zero_torsional_stiffness_identity(self) -> None:
        point = self._point()
        self.assertAlmostEqual(point.torsion.C_T.real, point.torsion.Cbar.real)
        self.assertAlmostEqual(saint_venant_stiffness(point), point.torsion.Cbar.real)

    def test_c_eb_state_entries_and_dimensions(self) -> None:
        point = self._point()
        omega = 2.0 * math.pi * 100.0
        matrix = eb_state_matrix(omega, point)
        self.assertEqual(matrix.shape, (6, 6))
        expected_nonzero = {(0, 1), (1, 4), (2, 5), (3, 0), (4, 3), (5, 2)}
        actual_nonzero = set(zip(*np.nonzero(matrix)))
        self.assertEqual(actual_nonzero, expected_nonzero)
        self.assertAlmostEqual(matrix[1, 4].real, 1.0 / (point.properties.Ex.real * point.geometry.I_y))
        self.assertAlmostEqual(matrix[2, 5].real, 1.0 / point.torsion.Cbar.real)
        self.assertAlmostEqual(matrix[3, 0].real, -point.material.rho * point.geometry.area * omega**2)
        self.assertAlmostEqual(matrix[5, 2].real, -point.material.rho * point.geometry.I_p * omega**2)

    def test_d_eb_clamps_and_end_maps(self) -> None:
        point = self._point()
        clamp = eb_clamp_matrix(point)
        np.testing.assert_array_equal(clamp[:3], 0.0)
        np.testing.assert_array_equal(clamp[3:], np.eye(3))
        omega = 2.0 * math.pi * 500.0
        np.testing.assert_allclose(
            eb_physical_end_map(omega, point),
            eb_state_transfer_matrix(omega, point) @ clamp,
            rtol=2e-14,
            atol=1e-13,
        )

    def test_e_coupled_and_straight_matrix_shapes(self) -> None:
        point = self._point()
        omega = 2.0 * math.pi * 500.0
        self.assertEqual(eb_coupled_boundary_matrix_raw(omega, math.pi / 6, point, point).shape, (6, 6))
        self.assertEqual(eb_coupled_boundary_matrix(omega, math.pi / 6, point, point).shape, (6, 6))
        straight = self._point(0.4)
        self.assertEqual(eb_straight_boundary_matrix_raw(omega, straight).shape, (3, 3))
        self.assertEqual(eb_straight_boundary_matrix(omega, straight).shape, (3, 3))
        np.testing.assert_array_equal(eb_straight_right_clamp_matrix(straight)[:, :3], np.eye(3))

    def test_f_fixed_fixed_bending_characteristic_roots(self) -> None:
        expected = np.array([4.730040744862704, 7.853204624095838, 10.99560783800167])
        roots = fixed_fixed_bending_dimensionless_roots(3)
        np.testing.assert_allclose(roots, expected, rtol=2e-14, atol=2e-14)
        np.testing.assert_allclose(np.cosh(roots) * np.cos(roots) - 1.0, 0.0, atol=2e-10)

    def test_g_straight_transfer_matches_exact_bending_family(self) -> None:
        point = self._point(0.4)
        expected = fixed_fixed_bending_frequencies_hz(point, 3)

        def builder(omega):
            raw = eb_straight_boundary_matrix_raw(omega, point)
            return equilibrate_matrix(raw[np.ix_([0, 1], [0, 1])])[0]

        roots = self._roots(point, builder, count=3)
        np.testing.assert_allclose([root.frequency_hz for root in roots], expected, rtol=2e-10, atol=2e-8)

    def test_h_straight_transfer_matches_exact_torsion_family(self) -> None:
        point = self._point(0.4)
        expected = fixed_fixed_torsion_frequencies_hz(point, 3)

        def builder(omega):
            raw = eb_straight_boundary_matrix_raw(omega, point)
            return raw[np.ix_([2], [2])]

        roots = self._roots(point, builder, count=3)
        np.testing.assert_allclose([root.frequency_hz for root in roots], expected, rtol=2e-10, atol=2e-8)

    def test_i_unequal_length_eb_split_invariance(self) -> None:
        straight = self._point(0.4)
        reference = self._roots(straight, lambda omega: eb_straight_boundary_matrix(omega, straight))
        reference_frequencies = np.array([root.frequency_hz for root in reference])
        for length_1, length_2 in ((0.2, 0.2), (0.15, 0.25), (0.1, 0.3), (0.05, 0.35)):
            point_1, point_2 = self._point(length_1), self._point(length_2)
            roots = self._roots(
                point_1,
                lambda omega, p1=point_1, p2=point_2: eb_coupled_boundary_matrix(omega, 0.0, p1, p2),
            )
            np.testing.assert_allclose(
                [root.frequency_hz for root in roots], reference_frequencies, rtol=1e-8, atol=3e-7
            )

    def test_j_unequal_length_timoshenko_split_invariance(self) -> None:
        straight = self._point(0.4)
        reference = self._roots(straight, lambda omega: straight_boundary_matrix(omega, straight))
        reference_frequencies = np.array([root.frequency_hz for root in reference])
        for length_1, length_2 in ((0.2, 0.2), (0.15, 0.25), (0.1, 0.3), (0.05, 0.35)):
            point_1, point_2 = self._point(length_1), self._point(length_2)
            roots = self._roots(
                point_1,
                lambda omega, p1=point_1, p2=point_2: coupled_boundary_matrix(omega, 0.0, p1, p2),
            )
            np.testing.assert_allclose(
                [root.frequency_hz for root in roots], reference_frequencies, rtol=1e-8, atol=3e-7
            )

    def test_k_beta_zero_orthotropic_blocks(self) -> None:
        point = self._point()
        matrix = eb_coupled_boundary_matrix_raw(2.0 * math.pi * 500.0, 0.0, point, point)
        bending_rows, torsion_rows = [0, 2, 3, 5], [1, 4]
        bending_columns, torsion_columns = [0, 1, 3, 4], [2, 5]
        off_block = np.concatenate(
            [
                matrix[np.ix_(bending_rows, torsion_columns)].ravel(),
                matrix[np.ix_(torsion_rows, bending_columns)].ravel(),
            ]
        )
        self.assertLess(np.linalg.norm(off_block), 2e-14 * max(np.linalg.norm(matrix), 1.0))
        self.assertEqual(eb_bending_coupled_boundary_matrix(2.0 * math.pi * 500.0, point, point).shape, (4, 4))

    def test_l_element_matrices_are_symmetric_with_positive_mass(self) -> None:
        stiffness, mass = eb_element_matrices(self._point(), 0.025)
        np.testing.assert_allclose(stiffness, stiffness.T, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(mass, mass.T, rtol=0.0, atol=0.0)
        self.assertGreater(np.linalg.eigvalsh(mass)[0], 0.0)
        self.assertEqual(np.linalg.matrix_rank(stiffness), 3)

    def test_m_fem_joint_maps_are_independent_and_compatible(self) -> None:
        for beta in (0.0, math.pi / 6.0, math.pi / 2.0):
            self.assertLessEqual(eb_joint_mapping_residual(beta), 2e-15)
            maps = eb_joint_dof_maps(beta)
            transform = np.zeros((12, 3))
            for arm, mapping in enumerate(maps):
                transform[6 * arm : 6 * arm + 3] = mapping
            np.testing.assert_allclose(joint_matrix_book(beta)[:3] @ transform, 0.0, atol=2e-15)

    def test_n_fem_global_symmetry_mass_and_no_zero_modes(self) -> None:
        point = self._point()
        assembly = assemble_two_arm_eb_fem(point, point, math.pi / 6.0, 8)
        solution = solve_two_arm_eb_fem(assembly, num_roots=7)
        self.assertLess(solution.stiffness_symmetry_residual, 2e-14)
        self.assertLess(solution.mass_symmetry_residual, 2e-14)
        self.assertGreater(solution.minimum_mass_eigenvalue, 0.0)
        self.assertEqual(solution.zero_mode_count, 0)
        self.assertEqual(len(solution.frequencies_hz), 7)
        self.assertTrue(np.all(np.diff(solution.frequencies_hz) > 0.0))

    def test_o_fem_joint_reactions_satisfy_congruence_equilibrium(self) -> None:
        point = self._point()
        for beta in (0.0, math.pi / 6.0, math.pi / 2.0):
            solution = solve_two_arm_eb_fem(
                assemble_two_arm_eb_fem(point, point, beta, 8), num_roots=7
            )
            self.assertLess(np.max(solution.joint_equilibrium_residuals), 2e-9)

    def test_p_beta_zero_fem_converges_to_exact_straight_spectrum(self) -> None:
        point = self._point()
        straight = self._point(0.4)
        exact = np.sort(
            np.concatenate(
                [fixed_fixed_bending_frequencies_hz(straight, 7), fixed_fixed_torsion_frequencies_hz(straight, 7)]
            )
        )[:7]
        errors = []
        for mesh in (8, 16, 32):
            frequencies = solve_two_arm_eb_fem(
                assemble_two_arm_eb_fem(point, point, 0.0, mesh), num_roots=7
            ).frequencies_hz
            errors.append(np.abs(frequencies - exact) / exact)
        self.assertTrue(np.all(errors[2][:6] <= 1.01 * errors[1][:6] + 2e-12))
        self.assertLess(np.max(errors[2][:3]), 1.1e-4)
        self.assertLess(np.max(errors[2][:6]), 1.0e-3)

    def test_q_fem_converges_for_angled_joint(self) -> None:
        point = self._point()
        reference = solve_two_arm_eb_fem(
            assemble_two_arm_eb_fem(point, point, math.pi / 6.0, 64), num_roots=7
        ).frequencies_hz
        errors = []
        for mesh in (8, 16, 32):
            frequencies = solve_two_arm_eb_fem(
                assemble_two_arm_eb_fem(point, point, math.pi / 6.0, mesh), num_roots=7
            ).frequencies_hz
            errors.append(np.abs(frequencies - reference) / reference)
        self.assertTrue(np.all(errors[2][:6] <= 1.01 * errors[1][:6] + 2e-12))

    def test_r_slender_beta_zero_bending_difference_decreases(self) -> None:
        differences = []
        for scale in (1.0, 0.5, 0.25):
            arm = self._point(scale=scale)
            eb_roots = self._roots(
                arm,
                lambda omega, p=arm: eb_bending_coupled_boundary_matrix(omega, p, p),
                count=3,
            )

            def timo_builder(omega, p=arm):
                raw = coupled_boundary_matrix_raw(omega, 0.0, p, p)
                return equilibrate_matrix(raw[np.ix_([0, 2, 3, 5], [0, 1, 3, 4])])[0]

            timo_roots = self._roots(arm, timo_builder, count=3)
            eb = np.array([root.frequency_hz for root in eb_roots])
            timo = np.array([root.frequency_hz for root in timo_roots])
            differences.append(np.abs(timo - eb) / timo)
        differences = np.asarray(differences)
        self.assertTrue(np.all(differences[1] <= differences[0] * 1.001 + 1e-10))
        self.assertTrue(np.all(differences[2] <= differences[1] * 1.001 + 1e-10))

    def test_s_targeted_a_distinct_arm_element_counts_are_accepted(self) -> None:
        point_1, point_2 = self._point(0.1), self._point(0.3)
        assembly = assemble_two_arm_eb_fem(point_1, point_2, 0.0, (1, 2))
        self.assertEqual(assembly.element_counts, (1, 2))
        self.assertEqual(assembly.stiffness.shape, (6, 6))

    def test_t_targeted_b_equal_count_api_is_backward_compatible(self) -> None:
        point = self._point()
        old = assemble_two_arm_eb_fem(point, point, math.pi / 6.0, 4)
        explicit = assemble_two_arm_eb_fem(point, point, math.pi / 6.0, (4, 4))
        self.assertEqual(old.element_counts, explicit.element_counts)
        np.testing.assert_array_equal(old.stiffness_full, explicit.stiffness_full)
        np.testing.assert_array_equal(old.mass_full, explicit.mass_full)
        np.testing.assert_array_equal(old.reduction, explicit.reduction)
        np.testing.assert_array_equal(old.stiffness, explicit.stiffness)
        np.testing.assert_array_equal(old.mass, explicit.mass)

    def test_u_targeted_c_dof_count_for_8_24(self) -> None:
        point_1, point_2 = self._point(0.1), self._point(0.3)
        assembly = assemble_two_arm_eb_fem(point_1, point_2, 0.0, (8, 24))
        self.assertEqual(assembly.stiffness_full.shape, (102, 102))
        self.assertEqual(assembly.stiffness.shape, (93, 93))
        self.assertEqual(assembly.reduction.shape, (102, 93))

    def test_v_targeted_d_proportional_meshes_have_equal_physical_steps(self) -> None:
        for n_1, n_2 in ((8, 24), (16, 48), (32, 96), (64, 192)):
            with self.subTest(n_1=n_1, n_2=n_2):
                h_1 = 0.1 / n_1
                h_2 = 0.3 / n_2
                self.assertLessEqual(abs(h_1 - h_2), 1e-15)

    def test_w_targeted_e_joint_mapping_survives_unequal_counts(self) -> None:
        point_1, point_2 = self._point(0.1), self._point(0.3)
        assembly = assemble_two_arm_eb_fem(
            point_1, point_2, math.pi / 6.0, (3, 7)
        )
        expected = eb_joint_dof_maps(math.pi / 6.0)
        np.testing.assert_array_equal(assembly.joint_maps[0], expected[0])
        np.testing.assert_array_equal(assembly.joint_maps[1], expected[1])
        self.assertLessEqual(eb_joint_mapping_residual(math.pi / 6.0), 2e-15)

    def test_x_targeted_f_structural_checks_on_proportional_mesh(self) -> None:
        assembly, solution, _, _ = self._targeted_solution((16, 48))
        self.assertEqual(assembly.element_counts, (16, 48))
        self.assertLess(solution.stiffness_symmetry_residual, 1e-12)
        self.assertLess(solution.mass_symmetry_residual, 1e-12)
        self.assertGreater(solution.minimum_mass_eigenvalue, 0.0)
        self.assertEqual(solution.zero_mode_count, 0)
        self.assertLess(np.max(solution.joint_equilibrium_residuals), 1e-7)

    def test_y_targeted_g_first_seven_positive_roots_are_returned(self) -> None:
        _, solution, _, _ = self._targeted_solution((16, 48))
        self.assertEqual(len(solution.frequencies_hz), 7)
        self.assertTrue(np.all(solution.frequencies_hz > 0.0))
        self.assertTrue(np.all(np.diff(solution.frequencies_hz) > 0.0))

    def test_z_targeted_h_records_unchanged_convergence_gate_result(self) -> None:
        errors = [
            self._targeted_solution(counts)[3]
            for counts in ((16, 48), (32, 96), (64, 192))
        ]
        allowance_16_32 = errors[1][:6] <= 1.01 * errors[0][:6] + 1e-10
        allowance_32_64 = errors[2][:6] <= 1.01 * errors[1][:6] + 1e-10
        np.testing.assert_array_equal(
            np.flatnonzero(~allowance_16_32), np.asarray([0])
        )
        np.testing.assert_array_equal(
            np.flatnonzero(~allowance_32_64), np.asarray([0])
        )
        self.assertTrue(np.all(allowance_16_32[1:]))
        self.assertTrue(np.all(allowance_32_64[1:]))

    def test_za_targeted_i_mode_2_order_is_second_order(self) -> None:
        errors = [
            self._targeted_solution(counts)[3][1]
            for counts in ((16, 48), (32, 96), (64, 192))
        ]
        orders = [math.log(errors[index] / errors[index + 1], 2.0) for index in (0, 1)]
        self.assertTrue(all(1.8 <= order <= 2.2 for order in orders))

    def test_zb_targeted_j_finest_raw_first_three_gate(self) -> None:
        _, _, _, errors = self._targeted_solution((64, 192))
        self.assertLessEqual(np.max(errors[:3]), 1e-5)


if __name__ == "__main__":
    unittest.main()
