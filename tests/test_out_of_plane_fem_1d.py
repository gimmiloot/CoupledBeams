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
    fixed_fixed_bending_roots_total_length_two,
    torsion_roots_uniform_eta0_beta0,
)
from my_project.analytic.out_of_plane_fem_1d import (  # noqa: E402
    assemble_out_of_plane_fem_1d_matrices,
    bending_element_matrices_phi,
    first_beta0_eta0_torsion_fem_root,
    solve_out_of_plane_fem_1d_frequencies,
    torsion_element_matrices,
)
from my_project.analytic.solvers_out_of_plane import (  # noqa: E402
    find_first_n_roots_out_of_plane,
)


class OutOfPlaneFem1DTest(unittest.TestCase):
    def test_element_matrices_are_symmetric(self) -> None:
        k_b, m_b = bending_element_matrices_phi(
            length=0.25,
            bending_stiffness=0.003,
            mass_per_length=1.2,
        )
        k_t, m_t = torsion_element_matrices(
            length=0.25,
            torsional_stiffness=0.004,
            polar_mass_per_length=0.002,
        )

        np.testing.assert_allclose(k_b, k_b.T, atol=1e-15)
        np.testing.assert_allclose(m_b, m_b.T, atol=1e-15)
        np.testing.assert_allclose(k_t, k_t.T, atol=1e-15)
        np.testing.assert_allclose(m_t, m_t.T, atol=1e-15)

    def test_reduced_global_matrices_are_symmetric_with_positive_mass(self) -> None:
        stiffness, mass = assemble_out_of_plane_fem_1d_matrices(
            beta=0.3,
            mu=0.2,
            epsilon=0.0025,
            eta=0.1,
            poisson=0.3,
            n_elements_per_rod=4,
        )

        np.testing.assert_allclose(stiffness, stiffness.T, atol=1e-12)
        np.testing.assert_allclose(mass, mass.T, atol=1e-12)
        self.assertGreater(float(np.min(np.linalg.eigvalsh(mass))), 0.0)

    def test_beta0_eta0_bending_matches_straight_clamped_clamped_roots(self) -> None:
        epsilon = 0.0025
        fem = solve_out_of_plane_fem_1d_frequencies(
            beta=0.0,
            mu=0.0,
            epsilon=epsilon,
            eta=0.0,
            poisson=0.3,
            n_elements_per_rod=24,
            n_modes=4,
        )
        exact = fixed_fixed_bending_roots_total_length_two()

        for mode_index, Lambda_exact in enumerate(exact):
            with self.subTest(mode_index=mode_index + 1):
                rel_error = abs(float(fem[mode_index]) - float(Lambda_exact)) / float(Lambda_exact)
                self.assertLess(rel_error, 5e-4)

    def test_beta0_eta0_first_torsion_root_visible(self) -> None:
        epsilon = 0.0025
        poisson = 0.3
        fem_lambda = first_beta0_eta0_torsion_fem_root(
            epsilon=epsilon,
            poisson=poisson,
            n_elements_per_rod=32,
        )
        exact_lambda = torsion_roots_uniform_eta0_beta0(1, epsilon, poisson)

        self.assertLess(abs(fem_lambda - exact_lambda) / exact_lambda, 6e-4)

    def test_analytic_vs_fem_small_smoke_comparison(self) -> None:
        epsilon = 0.0025
        beta = math.radians(15.0)
        analytic = np.asarray(
            find_first_n_roots_out_of_plane(
                beta=beta,
                mu=0.3,
                epsilon=epsilon,
                eta=0.1,
                poisson=0.3,
                n_roots=4,
                lambda_max=16.0,
            ),
            dtype=float,
        )
        fem = solve_out_of_plane_fem_1d_frequencies(
            beta=beta,
            mu=0.3,
            epsilon=epsilon,
            eta=0.1,
            poisson=0.3,
            n_elements_per_rod=16,
            n_modes=4,
        )

        self.assertEqual(len(analytic), 4)
        np.testing.assert_allclose(fem[:4], analytic, rtol=5e-4, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
