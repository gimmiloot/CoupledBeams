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

from my_project.analytic.formulas import (  # noqa: E402
    assemble_clamped_coupled_matrix,
    det_clamped_coupled,
)
from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    assemble_clamped_coupled_matrix_eta,
    det_eta,
    local_epsilons,
    thickness_mismatch_factors,
)


class ThicknessMismatchFormulaTest(unittest.TestCase):
    def test_eta_zero_has_unit_thickness_factors(self) -> None:
        for mu in (-0.6, 0.0, 0.3):
            with self.subTest(mu=mu):
                factors = thickness_mismatch_factors(mu, 0.0)
                self.assertAlmostEqual(1.0, factors.tau1, delta=1e-15)
                self.assertAlmostEqual(1.0, factors.tau2, delta=1e-15)
                self.assertAlmostEqual(1.0, factors.mass_factor / 2.0, delta=1e-15)
                self.assertEqual((0.0025, 0.0025), local_epsilons(0.0025, mu, 0.0))
                for value in (
                    factors.tau1 ** (-0.5),
                    factors.tau2 ** (-0.5),
                    factors.tau1**3,
                    factors.tau2**3,
                    factors.tau1 ** 2.5,
                    factors.tau2 ** 2.5,
                    factors.tau1**2,
                    factors.tau2**2,
                ):
                    self.assertAlmostEqual(1.0, value, delta=1e-15)

    def test_eta_zero_matrix_and_determinant_match_baseline(self) -> None:
        epsilon = 0.0025
        for Lambda in (2.5, 4.0, 7.0):
            for beta_deg in (0.0, 15.0):
                for mu in (0.0, 0.3):
                    with self.subTest(Lambda=Lambda, beta_deg=beta_deg, mu=mu):
                        beta = float(np.deg2rad(beta_deg))
                        baseline = assemble_clamped_coupled_matrix(Lambda, beta, mu, epsilon)
                        eta_matrix = assemble_clamped_coupled_matrix_eta(Lambda, beta, mu, epsilon, 0.0)
                        np.testing.assert_allclose(eta_matrix, baseline, rtol=1e-12, atol=1e-12)
                        self.assertAlmostEqual(
                            det_clamped_coupled(Lambda, beta, mu, epsilon),
                            det_eta(Lambda, beta, mu, epsilon, 0.0),
                            delta=1e-8 * max(1.0, abs(det_clamped_coupled(Lambda, beta, mu, epsilon))),
                        )

    def test_rod_exchange_symmetry_of_tau_factors(self) -> None:
        for mu, eta in ((0.3, 0.1), (0.8, 0.5), (-0.2, 0.4)):
            with self.subTest(mu=mu, eta=eta):
                factors = thickness_mismatch_factors(mu, eta)
                swapped = thickness_mismatch_factors(-mu, -eta)
                self.assertAlmostEqual(factors.tau1, swapped.tau2, delta=1e-14)
                self.assertAlmostEqual(factors.tau2, swapped.tau1, delta=1e-14)

    def test_invalid_parameters_are_rejected_and_denominator_is_positive(self) -> None:
        for mu in (-1.0, 1.0, 1.2):
            with self.subTest(mu=mu):
                with self.assertRaises(ValueError):
                    thickness_mismatch_factors(mu, 0.0)
        for eta in (-1.0, 1.0, 1.2):
            with self.subTest(eta=eta):
                with self.assertRaises(ValueError):
                    thickness_mismatch_factors(0.0, eta)

        factors = thickness_mismatch_factors(0.9, -0.9)
        self.assertGreater(factors.denom, 0.0)
        self.assertGreater(factors.tau1, 0.0)
        self.assertGreater(factors.tau2, 0.0)

    def test_single_rod_reference_formula_is_finite_positive(self) -> None:
        alpha = 4.730040744863
        for mu, eta in ((0.0, 0.0), (0.3, 0.5), (0.3, -0.5)):
            with self.subTest(mu=mu, eta=eta):
                factors = thickness_mismatch_factors(mu, eta)
                lambda_ref_1 = alpha * np.sqrt(factors.tau1) / (1.0 - mu)
                lambda_ref_2 = alpha * np.sqrt(factors.tau2) / (1.0 + mu)
                self.assertTrue(np.isfinite(lambda_ref_1))
                self.assertTrue(np.isfinite(lambda_ref_2))
                self.assertGreater(lambda_ref_1, 0.0)
                self.assertGreater(lambda_ref_2, 0.0)


if __name__ == "__main__":
    unittest.main()
