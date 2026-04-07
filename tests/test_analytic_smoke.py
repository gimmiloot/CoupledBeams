import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from my_project.analytic.FreqFromAngle import (  # noqa: E402
    BeamParams as AngleBeamParams,
    det_clamped_coupled as angle_det,
    find_first_n_roots as angle_find_first_n_roots,
    fixed_fixed_lambdas,
)
from my_project.analytic.FreqFromMu import (  # noqa: E402
    BeamParams as MuBeamParams,
    det_clamped_coupled as mu_det,
    find_first_n_roots as mu_find_first_n_roots,
)


class AnalyticSmokeTest(unittest.TestCase):
    def test_modules_share_beam_params(self):
        self.assertIs(AngleBeamParams, MuBeamParams)

    def test_det_baseline(self):
        cases = [
            (2.3, 0.4, 0.25, 0.0025),
            (4.1, 0.9, 0.4, 0.0018),
            (3.2, 0.2, 0.0, 0.0025),
        ]
        expected = np.array(
            [
                -2.8568375032250938e01,
                5.5839553362582237e03,
                -4.7214931163951405e01,
            ]
        )
        angle_values = np.array([angle_det(*case) for case in cases])
        mu_values = np.array([mu_det(*case) for case in cases])
        np.testing.assert_allclose(angle_values, expected, rtol=1e-13, atol=1e-13)
        np.testing.assert_allclose(mu_values, expected, rtol=1e-13, atol=1e-13)

    def test_root_baseline(self):
        expected_roots = np.array([3.57209560105, 5.297600200548, 6.335440244269, 8.612310726774])
        angle_roots = angle_find_first_n_roots(beta=0.4, mu=0.2, eps=0.0025, n_roots=4)
        mu_roots = mu_find_first_n_roots(beta=0.4, mu=0.2, eps=0.0025, n_roots=4)
        np.testing.assert_allclose(angle_roots, expected_roots, rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(mu_roots, expected_roots, rtol=1e-12, atol=1e-12)

    def test_fixed_fixed_baseline(self):
        expected = np.array([4.730040744863, 7.853204624096, 10.995607838002, 14.137165491257])
        np.testing.assert_allclose(fixed_fixed_lambdas(4), expected, rtol=1e-12, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
