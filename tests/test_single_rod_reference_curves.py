from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.FreqMuNet import roots_clamped_supported
from my_project.analytic.formulas_thickness_mismatch import thickness_mismatch_factors
from my_project.analytic.solvers import fixed_fixed_lambdas


class SingleRodReferenceCurveTests(unittest.TestCase):
    def test_eta_sign_thickness_logic_for_positive_mu(self) -> None:
        mu = 0.3
        positive_eta = thickness_mismatch_factors(mu, 0.1)
        self.assertGreater(positive_eta.tau2, positive_eta.tau1)
        self.assertLess(1.0 - mu, 1.0 + mu)

        negative_eta = thickness_mismatch_factors(mu, -0.1)
        self.assertGreater(negative_eta.tau1, negative_eta.tau2)
        self.assertLess(1.0 - mu, 1.0 + mu)

    def test_reference_formula_uses_sqrt_tau_over_current_length(self) -> None:
        mu = 0.3
        eta = 0.1
        alpha = float(roots_clamped_supported(1)[0])
        factors = thickness_mismatch_factors(mu, eta)

        lambda_rod1 = alpha * np.sqrt(factors.tau1) / (1.0 - mu)
        lambda_rod2 = alpha * np.sqrt(factors.tau2) / (1.0 + mu)

        self.assertAlmostEqual(lambda_rod1, 5.232318691534995, places=12)
        self.assertAlmostEqual(lambda_rod2, 3.114755517530622, places=12)

    def test_boundary_condition_root_equations(self) -> None:
        cp_alpha = float(roots_clamped_supported(1)[0])
        cc_alpha = float(fixed_fixed_lambdas(1)[0])

        self.assertAlmostEqual(float(np.tan(cp_alpha) - np.tanh(cp_alpha)), 0.0, places=12)
        self.assertAlmostEqual(float(np.cosh(cc_alpha) * np.cos(cc_alpha) - 1.0), 0.0, places=12)


if __name__ == "__main__":
    unittest.main()
