import math
import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from my_project.analytic.formulas_out_of_plane import (  # noqa: E402
    assemble_out_of_plane_beta0_blocks,
    assemble_out_of_plane_matrix,
    det_out_of_plane,
    fixed_fixed_bending_roots_total_length_two,
    out_of_plane_factors,
    torsion_roots_uniform_eta0_beta0,
)


E_Z = np.array([0.0, 0.0, 1.0], dtype=float)
T1 = np.array([1.0, 0.0, 0.0], dtype=float)
N1 = np.cross(E_Z, T1)


def local_basis(beta: float) -> tuple[np.ndarray, np.ndarray]:
    cb = math.cos(beta)
    sb = math.sin(beta)
    t2 = -T1 * cb + N1 * sb
    n2 = -T1 * sb - N1 * cb
    return t2, n2


def normalized_det(matrix: np.ndarray) -> float:
    matrix = np.asarray(matrix, dtype=float)
    det = float(np.linalg.det(matrix))
    row_scale = float(np.prod(np.maximum(np.linalg.norm(matrix, axis=1), 1.0)))
    if row_scale == 0.0 or not np.isfinite(row_scale):
        return det
    return det / row_scale


class OutOfPlaneFormulaTest(unittest.TestCase):
    def test_local_basis_identities(self) -> None:
        t2, n2 = local_basis(0.0)
        np.testing.assert_allclose(t2, -T1, atol=1e-15)
        np.testing.assert_allclose(n2, -N1, atol=1e-15)

        t2, n2 = local_basis(math.pi / 2.0)
        np.testing.assert_allclose(t2, N1, atol=1e-15)
        np.testing.assert_allclose(n2, -T1, atol=1e-15)

        for beta in (0.0, 0.3, math.pi / 2.0):
            with self.subTest(beta=beta):
                t2, n2 = local_basis(beta)
                self.assertAlmostEqual(float(np.dot(t2, t2)), 1.0, delta=1e-15)
                self.assertAlmostEqual(float(np.dot(n2, n2)), 1.0, delta=1e-15)
                self.assertAlmostEqual(float(np.dot(t2, n2)), 0.0, delta=1e-15)
                np.testing.assert_allclose(np.cross(t2, n2), E_Z, atol=1e-15)

    def test_rotation_compatibility_formulas(self) -> None:
        beta = 0.37
        psi2 = 1.4
        phi2 = -0.6
        t2, n2 = local_basis(beta)

        psi1 = -psi2 * math.cos(beta) - phi2 * math.sin(beta)
        phi1 = psi2 * math.sin(beta) - phi2 * math.cos(beta)
        theta1 = psi1 * T1 + phi1 * N1
        theta2 = psi2 * t2 + phi2 * n2

        np.testing.assert_allclose(theta1, theta2, atol=1e-15)

    def test_moment_compatibility_formulas(self) -> None:
        beta = 0.41
        T2 = 2.3
        M2 = -0.7
        t2, n2 = local_basis(beta)

        T1_component = T2 * math.cos(beta) + M2 * math.sin(beta)
        M1 = M2 * math.cos(beta) - T2 * math.sin(beta)
        m1 = T1_component * T1 + M1 * N1
        m2 = T2 * t2 + M2 * n2

        np.testing.assert_allclose(m1 + m2, np.zeros(3), atol=1e-15)

    def test_beta0_matrix_block_structure(self) -> None:
        matrix = assemble_out_of_plane_matrix(
            Lambda=3.1,
            beta=0.0,
            mu=0.3,
            eta=0.2,
            epsilon=0.0025,
            poisson=0.3,
        )

        np.testing.assert_allclose(matrix[np.ix_([0, 2, 3, 5], [4, 5])], 0.0, atol=1e-12)
        np.testing.assert_allclose(matrix[np.ix_([1, 4], [0, 1, 2, 3])], 0.0, atol=1e-12)

    def test_beta0_determinant_factorization(self) -> None:
        matrix = assemble_out_of_plane_matrix(
            Lambda=3.1,
            beta=0.0,
            mu=0.3,
            eta=0.2,
            epsilon=0.0025,
            poisson=0.3,
        )
        bending, torsion = assemble_out_of_plane_beta0_blocks(
            Lambda=3.1,
            mu=0.3,
            eta=0.2,
            epsilon=0.0025,
            poisson=0.3,
        )

        lhs = abs(float(np.linalg.det(matrix)))
        rhs = abs(float(np.linalg.det(bending) * np.linalg.det(torsion)))
        self.assertAlmostEqual(lhs, rhs, delta=1e-10 * max(1.0, lhs, rhs))

    def test_eta_zero_out_of_plane_factors_are_unity(self) -> None:
        for mu in (0.0, 0.3, -0.4):
            with self.subTest(mu=mu):
                factors = out_of_plane_factors(
                    Lambda=2.7,
                    beta=0.2,
                    mu=mu,
                    epsilon=0.0025,
                    eta=0.0,
                    poisson=0.3,
                )
                for name in (
                    "tau1",
                    "tau2",
                    "a1",
                    "a2",
                    "b1",
                    "b2",
                    "c1",
                    "c2",
                    "e1",
                    "e2",
                ):
                    self.assertAlmostEqual(getattr(factors, name), 1.0, delta=1e-15)

    def test_eta0_beta0_bending_known_roots(self) -> None:
        roots = fixed_fixed_bending_roots_total_length_two()
        for mu in (0.0, 0.3, 0.7):
            for Lambda in roots:
                with self.subTest(mu=mu, Lambda=Lambda):
                    bending, _ = assemble_out_of_plane_beta0_blocks(
                        Lambda=Lambda,
                        mu=mu,
                        epsilon=0.0025,
                        eta=0.0,
                        poisson=0.3,
                    )
                    self.assertLess(abs(normalized_det(bending)), 1e-10)

            bending, _ = assemble_out_of_plane_beta0_blocks(
                Lambda=3.0,
                mu=mu,
                epsilon=0.0025,
                eta=0.0,
                poisson=0.3,
            )
            self.assertGreater(abs(normalized_det(bending)), 1e-6)

    def test_eta0_beta0_torsion_known_roots(self) -> None:
        epsilon = 0.0025
        poisson = 0.3
        for mu in (0.0, 0.3, 0.7):
            for n in (1, 2, 3):
                Lambda = torsion_roots_uniform_eta0_beta0(n, epsilon, poisson)
                with self.subTest(mu=mu, n=n, Lambda=Lambda):
                    _, torsion = assemble_out_of_plane_beta0_blocks(
                        Lambda=Lambda,
                        mu=mu,
                        epsilon=epsilon,
                        eta=0.0,
                        poisson=poisson,
                    )
                    self.assertLess(abs(normalized_det(torsion)), 1e-10)

            _, torsion = assemble_out_of_plane_beta0_blocks(
                Lambda=10.0,
                mu=mu,
                epsilon=epsilon,
                eta=0.0,
                poisson=poisson,
            )
            self.assertGreater(abs(normalized_det(torsion)), 1e-3)

    def test_determinant_not_identically_zero(self) -> None:
        values = []
        for Lambda in (1.3, 2.7, 4.1):
            for beta in (0.0, 0.2):
                values.append(
                    abs(
                        det_out_of_plane(
                            Lambda=Lambda,
                            beta=beta,
                            mu=0.2,
                            eta=0.1,
                            epsilon=0.0025,
                            poisson=0.3,
                        )
                    )
                )

        self.assertTrue(any(value > 1e-12 for value in values))

    def test_invalid_eta_values_are_rejected(self) -> None:
        valid = dict(Lambda=2.0, beta=0.0, mu=0.0, epsilon=0.0025, poisson=0.3)
        for eta in (1.0, -1.0, 1.2):
            params = dict(valid, eta=eta)
            with self.subTest(eta=eta):
                with self.assertRaises(ValueError):
                    assemble_out_of_plane_matrix(**params)

    def test_invalid_poisson_values_are_rejected(self) -> None:
        valid = dict(Lambda=2.0, beta=0.0, mu=0.0, epsilon=0.0025, eta=0.0)
        for poisson in (-1.0, -1.2):
            params = dict(valid, poisson=poisson)
            with self.subTest(poisson=poisson):
                with self.assertRaises(ValueError):
                    assemble_out_of_plane_matrix(**params)

    def test_valid_poisson_ratio_produces_finite_torsion_factors(self) -> None:
        factors = out_of_plane_factors(
            Lambda=2.0,
            beta=0.0,
            mu=0.2,
            epsilon=0.0025,
            eta=0.1,
            poisson=0.3,
        )

        for value in (factors.gamma1, factors.gamma2, factors.chi_T):
            self.assertTrue(math.isfinite(value))
            self.assertGreater(value, 0.0)

    def test_invalid_parameters_are_rejected(self) -> None:
        valid = dict(Lambda=2.0, beta=0.0, mu=0.0, epsilon=0.0025, eta=0.0, poisson=0.3)
        for name, value in (("Lambda", 0.0), ("epsilon", 0.0), ("poisson", -1.0), ("mu", 1.0)):
            params = dict(valid)
            params[name] = value
            with self.subTest(name=name):
                with self.assertRaises(ValueError):
                    assemble_out_of_plane_matrix(**params)

if __name__ == "__main__":
    unittest.main()
