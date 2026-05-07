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

from scripts.lib.analytic_branch_tracking import (  # noqa: E402
    branch_id_from_base_sorted_index,
    dense_mu_values_for_targets,
    track_mu_sweep,
)


class AnalyticBranchTrackingRegressionTest(unittest.TestCase):
    BETA = 15.0
    EPSILON = 0.0025
    BRANCH_ID = branch_id_from_base_sorted_index(5)
    CHECK_MUS = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 0.9], dtype=float)

    @classmethod
    def setUpClass(cls):
        tracking_mu_values = dense_mu_values_for_targets(cls.CHECK_MUS, mu_steps=91)
        cls.result = track_mu_sweep(
            epsilon=cls.EPSILON,
            beta=cls.BETA,
            mu_values=tracking_mu_values,
            n_track=12,
            n_solve=20,
            shape_metric="full",
            required_branch_ids=[cls.BRANCH_ID],
        )

    def test_desc05_current_sorted_index_is_reproducible(self):
        for mu in self.CHECK_MUS:
            with self.subTest(mu=float(mu)):
                point = self.result.point_at(self.BRANCH_ID, beta=self.BETA, mu=float(mu))
                self.assertEqual(6, int(point.current_sorted_index))

        point_mu08 = self.result.point_at(self.BRANCH_ID, beta=self.BETA, mu=0.8)
        self.assertAlmostEqual(10.90652269, float(point_mu08.Lambda), delta=5e-7)
        self.assertGreater(abs(float(point_mu08.Lambda) - 12.49256427), 0.1)

    def test_branch_lambda_matches_current_sorted_root(self):
        for mu in self.CHECK_MUS:
            with self.subTest(mu=float(mu)):
                point = self.result.point_at(self.BRANCH_ID, beta=self.BETA, mu=float(mu))
                sorted_lambda = float(point.sorted_lambdas[int(point.current_sorted_index) - 1])
                self.assertAlmostEqual(sorted_lambda, float(point.Lambda), delta=1e-9)

    def test_coarse_low_mac_path_is_not_canonical(self):
        coarse_mu_values = dense_mu_values_for_targets([0.8], mu_steps=41)
        with self.assertRaisesRegex(RuntimeError, "Low-MAC analytic branch assignment"):
            track_mu_sweep(
                epsilon=self.EPSILON,
                beta=self.BETA,
                mu_values=coarse_mu_values,
                n_track=12,
            n_solve=20,
            shape_metric="full",
            required_branch_ids=[self.BRANCH_ID],
            max_refinement_depth=0,
        )

        exploratory = track_mu_sweep(
            epsilon=self.EPSILON,
            beta=self.BETA,
            mu_values=coarse_mu_values,
            n_track=12,
            n_solve=20,
            shape_metric="full",
            allow_low_mac=True,
            required_branch_ids=[self.BRANCH_ID],
        )
        point = exploratory.point_at(self.BRANCH_ID, beta=self.BETA, mu=0.8)
        branch_macs = [
            float(branch_point.mac_to_previous)
            for branch_point in exploratory.points_for_branch(self.BRANCH_ID)
            if branch_point.step_type != "base"
        ]
        self.assertLess(min(branch_macs), 0.8)
        self.assertEqual(7, int(point.current_sorted_index))
        self.assertAlmostEqual(12.49256427, float(point.Lambda), delta=5e-7)

    def test_desc05_beta30_continuation_has_no_low_mac_failure(self):
        result = track_mu_sweep(
            epsilon=self.EPSILON,
            beta=30.0,
            mu_values=[0.0],
            n_track=12,
            n_solve=20,
            shape_metric="full",
            required_branch_ids=[self.BRANCH_ID],
        )
        point = result.point_at(self.BRANCH_ID, beta=30.0, mu=0.0)
        self.assertEqual("ok", point.warning_flag)
        branch_macs = [
            float(branch_point.mac_to_previous)
            for branch_point in result.points_for_branch(self.BRANCH_ID)
            if branch_point.step_type != "base"
        ]
        self.assertGreaterEqual(min(branch_macs), 0.8)

    def test_desc05_beta45_adaptive_continuation_does_not_accept_tiny_mac(self):
        result = track_mu_sweep(
            epsilon=self.EPSILON,
            beta=45.0,
            mu_values=[0.0],
            n_track=12,
            n_solve=20,
            shape_metric="full",
            required_branch_ids=[self.BRANCH_ID],
        )
        point = result.point_at(self.BRANCH_ID, beta=45.0, mu=0.0)
        self.assertEqual("ok", point.warning_flag)
        branch_macs = [
            float(branch_point.mac_to_previous)
            for branch_point in result.points_for_branch(self.BRANCH_ID)
            if branch_point.step_type != "base"
        ]
        self.assertGreaterEqual(min(branch_macs), 0.8)


if __name__ == "__main__":
    unittest.main()
