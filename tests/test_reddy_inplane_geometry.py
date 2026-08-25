import ast
import sys
import unittest
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib import in_plane_shape_geometry as legacy_display
from scripts.lib import reddy_inplane_geometry as rlb


ATOL = 5.0e-14


class ReddyInPlaneCoordinateGateTest(unittest.TestCase):
    def assertVectorEqual(self, actual, expected, *, atol=ATOL) -> None:
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=atol)

    def test_global_view_basis_and_k(self) -> None:
        self.assertVectorEqual(np.cross(rlb.E_X, rlb.E_Y), rlb.E_Z)
        self.assertVectorEqual(rlb.k, -rlb.E_Z)
        self.assertVectorEqual(rlb.E_Z, np.array([0.0, 0.0, 1.0]))
        self.assertVectorEqual(rlb.k, np.array([0.0, 0.0, -1.0]))
        self.assertEqual(rlb.STATE_ORDER, ("u", "w", "psi", "N", "Q", "M"))

    def test_explicit_beta_zero_thirty_and_ninety_geometry(self) -> None:
        root_three_over_two = np.sqrt(3.0) / 2.0
        expected = {
            0.0: (
                np.array([1.0, 0.0, 0.0]),
                np.array([-1.0, 0.0, 0.0]),
                np.array([0.0, 1.0, 0.0]),
            ),
            30.0: (
                np.array([root_three_over_two, 0.5, 0.0]),
                np.array([-root_three_over_two, -0.5, 0.0]),
                np.array([-0.5, root_three_over_two, 0.0]),
            ),
            90.0: (
                np.array([0.0, 1.0, 0.0]),
                np.array([0.0, -1.0, 0.0]),
                np.array([-1.0, 0.0, 0.0]),
            ),
        }
        for beta_deg, (g2, t2, n2) in expected.items():
            with self.subTest(beta_deg=beta_deg):
                geometry = rlb.reddy_inplane_geometry(beta_deg)
                self.assertVectorEqual(geometry.g2, g2)
                self.assertVectorEqual(geometry.t2, t2)
                self.assertVectorEqual(geometry.n2, n2)

        beta_zero = rlb.reddy_inplane_geometry(0.0)
        with self.assertRaisesRegex(ValueError, "inconsistent with beta_deg"):
            rlb.ReddyInPlaneGeometry(
                beta_deg=30.0,
                g2=beta_zero.g2,
                arm1=beta_zero.arm1,
                arm2=beta_zero.arm2,
            )

    def test_beta_zero_and_ninety_limits(self) -> None:
        beta_zero = rlb.reddy_inplane_geometry(0.0)
        self.assertVectorEqual(beta_zero.t2, -beta_zero.t1)
        self.assertVectorEqual(beta_zero.n2, -beta_zero.n1)

        beta_ninety = rlb.reddy_inplane_geometry(90.0)
        self.assertVectorEqual(beta_ninety.t2, beta_ninety.n1)
        self.assertVectorEqual(beta_ninety.n2, -beta_ninety.t1)

    def test_arm_bases_and_reddy_triads_are_orthonormal_and_right_handed(self) -> None:
        for beta_deg in (0.0, 30.0, 90.0, -37.0, 143.0):
            geometry = rlb.reddy_inplane_geometry(beta_deg)
            for arm_number, basis in ((1, geometry.arm1), (2, geometry.arm2)):
                with self.subTest(beta_deg=beta_deg, arm=arm_number):
                    self.assertVectorEqual(basis.n, np.cross(rlb.k, basis.t))
                    self.assertAlmostEqual(
                        float(np.dot(basis.t, basis.n)), 0.0, delta=ATOL
                    )
                    self.assertAlmostEqual(
                        float(np.linalg.norm(basis.t)), 1.0, delta=ATOL
                    )
                    self.assertAlmostEqual(
                        float(np.linalg.norm(basis.n)), 1.0, delta=ATOL
                    )
                    self.assertVectorEqual(np.cross(basis.t, basis.n), rlb.k)
                    self.assertVectorEqual(np.cross(basis.e_x, basis.e_y), basis.e_z)
                    self.assertVectorEqual(basis.e_y, -rlb.k)
                    self.assertVectorEqual(
                        basis.reddy_matrix.T @ basis.reddy_matrix,
                        np.eye(3),
                    )
                    self.assertAlmostEqual(
                        float(np.linalg.det(basis.reddy_matrix)), 1.0, delta=ATOL
                    )

        first = rlb.arm1_basis()
        self.assertVectorEqual(first.e_x, rlb.E_X)
        self.assertVectorEqual(first.e_y, rlb.E_Z)
        self.assertVectorEqual(first.e_z, -rlb.E_Y)

    def test_exact_relations_to_legacy_timoshenko_display_geometry(self) -> None:
        new_arm1 = rlb.arm1_basis()
        old_arm1 = legacy_display.rod1_display_basis()
        self.assertVectorEqual(old_arm1.tangent, new_arm1.t[:2])
        self.assertVectorEqual(old_arm1.normal, new_arm1.n[:2])

        for beta_deg in (0.0, 30.0, 90.0, 127.0):
            with self.subTest(beta_deg=beta_deg):
                new_arm2 = rlb.arm2_basis(beta_deg)
                old_arm2 = legacy_display.rod2_display_basis(beta_deg)
                self.assertVectorEqual(old_arm2.tangent, -new_arm2.t[:2])
                self.assertVectorEqual(old_arm2.normal, -new_arm2.n[:2])

    def test_legacy_eb_transverse_field_sign_is_explicitly_different(self) -> None:
        new_arm1 = rlb.arm1_basis()
        old_eb_arm1 = legacy_display.eb_rod1_display_basis()
        self.assertVectorEqual(old_eb_arm1.tangent, new_arm1.t[:2])
        self.assertVectorEqual(old_eb_arm1.normal, -new_arm1.n[:2])

        for beta_deg in (0.0, 30.0, 90.0, 127.0):
            with self.subTest(beta_deg=beta_deg):
                new_arm2 = rlb.arm2_basis(beta_deg)
                old_eb_arm2 = legacy_display.eb_rod2_display_basis(beta_deg)
                self.assertVectorEqual(old_eb_arm2.tangent, -new_arm2.t[:2])
                self.assertVectorEqual(old_eb_arm2.normal, new_arm2.n[:2])

    def test_new_physical_helper_does_not_import_legacy_display_basis(self) -> None:
        helper_path = REPO_ROOT / "scripts" / "lib" / "reddy_inplane_geometry.py"
        source = helper_path.read_text(encoding="utf-8")
        tree = ast.parse(source)
        imported_modules = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported_modules.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module is not None:
                imported_modules.add(node.module)
        self.assertNotIn("scripts.lib.in_plane_shape_geometry", imported_modules)
        self.assertNotIn("in_plane_shape_geometry", imported_modules)
        self.assertNotIn("rod2_display_basis", source)

    def test_seeded_random_vector_reconstruction(self) -> None:
        rng = np.random.default_rng(20260825)
        beta_values = np.concatenate(
            (np.array([0.0, 30.0, 90.0]), rng.uniform(-175.0, 175.0, 32))
        )
        for beta_deg in beta_values:
            geometry = rlb.reddy_inplane_geometry(float(beta_deg))
            for basis in (geometry.arm1, geometry.arm2):
                components = rng.normal(size=3)
                global_vector = basis.reddy_matrix @ components
                recovered = basis.reddy_matrix.T @ global_vector
                self.assertVectorEqual(recovered, components, atol=1.0e-13)

                u, w = rng.normal(size=2)
                displacement = rlb.displacement_vector(basis, float(u), float(w))
                self.assertVectorEqual(
                    basis.inplane_matrix.T @ displacement,
                    np.array([u, w]),
                )
                self.assertAlmostEqual(
                    float(np.dot(displacement, basis.e_y)), 0.0, delta=ATOL
                )

    def test_local_and_global_virtual_work_pairings_are_equal(self) -> None:
        rng = np.random.default_rng(20260825)
        beta_values = np.concatenate(
            (np.array([0.0, 30.0, 90.0]), rng.uniform(-175.0, 175.0, 32))
        )
        for beta_deg in beta_values:
            geometry = rlb.reddy_inplane_geometry(float(beta_deg))
            for basis in (geometry.arm1, geometry.arm2):
                n_force, q_force, moment, delta_u, delta_w, delta_psi = rng.normal(
                    size=6
                )
                force = rlb.force_vector(basis, float(n_force), float(q_force))
                delta_displacement = rlb.displacement_vector(
                    basis, float(delta_u), float(delta_w)
                )
                moment_vector = rlb.moment_vector(basis, float(moment))
                delta_rotation = rlb.rotation_vector(basis, float(delta_psi))

                translational_local = n_force * delta_u + q_force * delta_w
                rotational_local = moment * delta_psi
                self.assertAlmostEqual(
                    float(np.dot(force, delta_displacement)),
                    float(translational_local),
                    delta=1.0e-13,
                )
                self.assertAlmostEqual(
                    float(np.dot(moment_vector, delta_rotation)),
                    float(rotational_local),
                    delta=1.0e-13,
                )
                self.assertAlmostEqual(
                    float(
                        np.dot(force, delta_displacement)
                        + np.dot(moment_vector, delta_rotation)
                    ),
                    float(translational_local + rotational_local),
                    delta=2.0e-13,
                )
                self.assertVectorEqual(delta_rotation, -delta_psi * rlb.k)
                self.assertVectorEqual(moment_vector, -moment * rlb.k)


if __name__ == "__main__":
    unittest.main()
