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
    joint_basis,
    joint_matrix_book,
    joint_matrix_old,
    joint_virtual_work,
    kinematic_joint_residual,
    moment_joint_residual,
    notation_transform_matrix,
    physical_end_map,
    physical_moment_vector,
    physical_rotation_vector,
    straight_boundary_matrix,
    straight_boundary_matrix_raw,
    straight_right_clamp_matrix,
    two_arm_notation_transform_matrix,
)
from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
    cantilever_clamp_matrix,
    cantilever_geometry,
    find_elastic_roots,
    hms_dx_209_material,
    make_rod_point,
    physical_state_transfer_matrix,
)


class YartsevChapter2CoupledRodsTest(unittest.TestCase):
    @staticmethod
    def _point(theta: float, length: float = 0.2):
        return make_rod_point(
            theta,
            geometry=cantilever_geometry(length),
            material=hms_dx_209_material(),
        )

    def test_notation_transform_matrix(self) -> None:
        expected = np.diag([1.0, -1.0, 1.0, 1.0, -1.0, 1.0])
        np.testing.assert_array_equal(notation_transform_matrix(), expected)
        expected_two = np.zeros((12, 12))
        expected_two[:6, :6] = expected
        expected_two[6:, 6:] = expected
        np.testing.assert_array_equal(two_arm_notation_transform_matrix(), expected_two)

    def test_basis_identities_at_beta_zero_and_ninety(self) -> None:
        zero = joint_basis(0.0)
        np.testing.assert_allclose(zero.t_2, -zero.t_1, atol=1e-15)
        np.testing.assert_allclose(zero.n_2, -zero.n_1, atol=1e-15)
        right = joint_basis(math.pi / 2.0)
        np.testing.assert_allclose(right.t_2, right.n_1, atol=1e-15)
        np.testing.assert_allclose(right.n_2, -right.t_1, atol=1e-15)
        for basis in (zero, right, joint_basis(math.pi / 6.0)):
            np.testing.assert_allclose(np.cross(basis.e_z, basis.t_1), basis.n_1)
            np.testing.assert_allclose(np.cross(basis.t_1, basis.n_1), basis.e_z)
            np.testing.assert_allclose(np.cross(basis.e_z, basis.t_2), basis.n_2)
            np.testing.assert_allclose(np.cross(basis.t_2, basis.n_2), basis.e_z)

    def test_exact_scalar_rows_of_joint_matrix(self) -> None:
        beta = math.radians(30.0)
        c, s = math.cos(beta), math.sin(beta)
        expected = np.array(
            [
                [1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0, 0, -s, c, 0, 0, 0],
                [0, 1, 0, 0, 0, 0, 0, c, s, 0, 0, 0],
                [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
                [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, s, -c],
                [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -c, -s],
            ]
        )
        np.testing.assert_allclose(joint_matrix_book(beta), expected, atol=0.0)

    def test_old_project_sign_gate(self) -> None:
        transform = two_arm_notation_transform_matrix()
        for beta_deg in (0.0, 15.0, 30.0, 90.0, 180.0):
            beta = math.radians(beta_deg)
            with self.subTest(beta_deg=beta_deg):
                np.testing.assert_allclose(
                    joint_matrix_book(beta),
                    joint_matrix_old(beta) @ transform,
                    rtol=0.0,
                    atol=2e-16,
                )

    def test_joint_rank(self) -> None:
        for beta_deg in (0.0, 15.0, 30.0, 90.0, 180.0):
            with self.subTest(beta_deg=beta_deg):
                self.assertEqual(
                    np.linalg.matrix_rank(joint_matrix_book(math.radians(beta_deg))),
                    6,
                )

    def test_random_compatible_rotation_and_moment_residuals(self) -> None:
        rng = np.random.default_rng(20260731)
        for beta_deg in (0.0, 30.0, 90.0):
            basis = joint_basis(math.radians(beta_deg))
            for _ in range(100):
                rotation = rng.normal() * basis.t_1 + rng.normal() * basis.n_1
                moment = rng.normal() * basis.t_1 + rng.normal() * basis.n_1
                self.assertLess(
                    np.linalg.norm(kinematic_joint_residual(math.radians(beta_deg), rotation)),
                    2e-14,
                )
                self.assertLess(
                    np.linalg.norm(moment_joint_residual(math.radians(beta_deg), moment)),
                    2e-14,
                )

    def test_virtual_work_and_book_energy_pairing(self) -> None:
        rng = np.random.default_rng(20260731)
        for beta_deg in (0.0, 30.0, 90.0):
            basis = joint_basis(math.radians(beta_deg))
            for _ in range(100):
                delta_w = float(rng.normal())
                q_1 = float(rng.normal())
                delta_rotation = rng.normal() * basis.t_1 + rng.normal() * basis.n_1
                moment_1 = rng.normal() * basis.t_1 + rng.normal() * basis.n_1
                work = joint_virtual_work(
                    q_1,
                    -q_1,
                    delta_w,
                    delta_w,
                    moment_1,
                    -moment_1,
                    delta_rotation,
                    delta_rotation,
                )
                self.assertLess(abs(work), 2e-14)

                Phi, psi = rng.normal(size=2)
                dPhi, dpsi = rng.normal(size=2)
                M_T, M = rng.normal(size=2)
                vector_pairing = physical_moment_vector(
                    M_T, M, basis.t_1, basis.n_1
                ) @ physical_rotation_vector(dPhi, dpsi, basis.t_1, basis.n_1)
                self.assertAlmostEqual(vector_pairing, M_T * dPhi + M * dpsi)

    def test_book_slope_clamp_and_physical_end_map(self) -> None:
        point = self._point(15.0)
        clamp = cantilever_clamp_matrix(point, "book_slope_clamp")
        shear_stiffness = (
            point.geometry.shear_factor
            * point.properties.Gxz
            * point.geometry.area
        )
        self.assertEqual(clamp.shape, (6, 3))
        np.testing.assert_allclose(clamp[0], 0.0)
        np.testing.assert_allclose(clamp[2], 0.0)
        self.assertAlmostEqual((clamp[1, 0] + clamp[3, 0] / shear_stiffness).real, 0.0)
        omega = 2.0 * math.pi * 500.0
        np.testing.assert_allclose(
            physical_end_map(omega, point),
            physical_state_transfer_matrix(omega, point) @ clamp,
            rtol=2e-14,
            atol=1e-12,
        )

    def test_coupled_and_straight_matrix_shapes(self) -> None:
        point = self._point(15.0)
        omega = 2.0 * math.pi * 500.0
        self.assertEqual(coupled_boundary_matrix(omega, math.pi / 6.0, point, point).shape, (6, 6))
        self.assertEqual(coupled_boundary_matrix_raw(omega, math.pi / 6.0, point, point).shape, (6, 6))
        straight = self._point(15.0, 0.4)
        self.assertEqual(straight_boundary_matrix(omega, straight).shape, (3, 3))
        self.assertEqual(straight_boundary_matrix_raw(omega, straight).shape, (3, 3))
        selector = straight_right_clamp_matrix(straight)
        self.assertEqual(selector.shape, (3, 6))
        shear_stiffness = straight.geometry.shear_factor * straight.properties.Gxz * straight.geometry.area
        np.testing.assert_allclose(selector[1], [0, 1, 0, 1 / shear_stiffness, 0, 0])

    def test_beta_zero_spectrum_equivalence(self) -> None:
        arm = self._point(15.0)
        straight = self._point(15.0, 0.4)

        def coupled_builder(omega, _point, _formulation):
            return coupled_boundary_matrix(omega, 0.0, arm, arm)

        def straight_builder(omega, _point, _formulation):
            return straight_boundary_matrix(omega, straight)

        coupled = find_elastic_roots(
            arm,
            "state_corrected",
            num_roots=7,
            scan_step_hz=10.0,
            initial_max_hz=5000.0,
            boundary_matrix_builder=coupled_builder,
        )
        reference = find_elastic_roots(
            straight,
            "state_corrected",
            num_roots=7,
            scan_step_hz=10.0,
            initial_max_hz=5000.0,
            boundary_matrix_builder=straight_builder,
        )
        np.testing.assert_allclose(
            [item.frequency_hz for item in coupled],
            [item.frequency_hz for item in reference],
            rtol=1e-8,
            atol=2e-7,
        )

    def test_orthotropic_beta_zero_block_separation(self) -> None:
        point = self._point(0.0)
        self.assertLess(abs(point.properties.Sbar16), 1e-24)
        matrix = coupled_boundary_matrix_raw(2.0 * math.pi * 500.0, 0.0, point, point)
        bending_rows = [0, 2, 3, 5]
        torsion_rows = [1, 4]
        bending_columns = [0, 1, 3, 4]
        torsion_columns = [2, 5]
        off_block = np.block(
            [
                [matrix[np.ix_(bending_rows, torsion_columns)].ravel()],
                [matrix[np.ix_(torsion_rows, bending_columns)].ravel()],
            ]
        ).ravel()
        self.assertLessEqual(np.linalg.norm(off_block), 2e-14 * max(np.linalg.norm(matrix), 1.0))
        self.assertEqual(np.linalg.matrix_rank(matrix), 6)

    def test_equal_arm_exchange_preserves_singular_values(self) -> None:
        point = self._point(15.0)
        matrix = coupled_boundary_matrix(2.0 * math.pi * 500.0, math.pi / 6.0, point, point)
        swapped = matrix[:, [3, 4, 5, 0, 1, 2]]
        np.testing.assert_allclose(
            np.linalg.svd(matrix, compute_uv=False),
            np.linalg.svd(swapped, compute_uv=False),
            rtol=2e-14,
            atol=2e-14,
        )


if __name__ == "__main__":
    unittest.main()
