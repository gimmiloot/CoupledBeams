import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from my_project.analytic.formulas_thickness_mismatch import thickness_mismatch_factors  # noqa: E402
from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_eb_validity_vs_timoshenko_stage1 as audit,
)


class EbValidityMetricTest(unittest.TestCase):
    def test_delta_f_uses_lambda_squared(self) -> None:
        metrics = audit.frequency_metrics(2.0, 1.0)
        self.assertAlmostEqual(metrics["delta_Lambda"], 1.0)
        self.assertAlmostEqual(metrics["delta_f"], 3.0)
        self.assertAlmostEqual(metrics["bias_f"], 3.0)

    def test_identical_lambda_has_zero_frequency_divergence(self) -> None:
        metrics = audit.frequency_metrics(4.25, 4.25)
        self.assertAlmostEqual(metrics["delta_Lambda"], 0.0)
        self.assertAlmostEqual(metrics["delta_f"], 0.0)
        self.assertAlmostEqual(metrics["delta_f_symmetric"], 0.0)
        self.assertAlmostEqual(metrics["bias_f"], 0.0)

    def test_consecutive_pass_count_stops_at_first_failure(self) -> None:
        deltas = [0.01, 0.08, 0.12, 0.04]
        self.assertEqual(audit.consecutive_pass_count(deltas, 0.10), 2)
        self.assertEqual(audit.count_passes(deltas, 0.10), 3)
        self.assertEqual(audit.pass_pattern(deltas, 0.10), "Y,Y,N,Y")

    def test_local_epsilon_formulas_include_local_lengths(self) -> None:
        epsilon_0 = 0.02
        mu = 0.3
        eta = 0.2
        factors = thickness_mismatch_factors(mu, eta)
        values = audit.local_thickness_parameters(epsilon_0, mu, eta)
        self.assertAlmostEqual(values["epsilon_1"], epsilon_0 * factors.tau1 / (1.0 - mu))
        self.assertAlmostEqual(values["epsilon_2"], epsilon_0 * factors.tau2 / (1.0 + mu))

    def test_synthetic_monotone_threshold_bracketing_and_refinement(self) -> None:
        xs = [0.0, 1.0, 2.0]
        ys = [0.02, 0.08, 0.12]
        brackets = audit.find_threshold_brackets(xs, ys, 0.10)
        self.assertEqual(len(brackets), 1)
        self.assertEqual(brackets[0][:2], (1.0, 2.0))

        root, value, iterations = audit.bisect_threshold(
            lambda x: 0.03 + 0.05 * x,
            1.0,
            2.0,
            threshold=0.10,
            eps_tol=1.0e-6,
            delta_tol=1.0e-8,
            max_iterations=40,
        )
        self.assertLess(iterations, 40)
        self.assertAlmostEqual(root, 1.4, delta=1.0e-5)
        self.assertAlmostEqual(value, 0.10, delta=1.0e-6)

    def test_synthetic_multiple_crossing_detection(self) -> None:
        xs = [0, 1, 2, 3, 4, 5]
        ys = [0.01, 0.12, 0.05, 0.11, 0.08, 0.13]
        brackets = audit.find_threshold_brackets(xs, ys, 0.10)
        self.assertEqual(len(brackets), 3)
        self.assertTrue(audit.is_nonmonotonic(ys))

    def test_energy_fractions_sum_to_one_for_synthetic_partitions(self) -> None:
        eb = audit.shape_audit.energy_dict(
            model=audit.MODEL_EB,
            U_axial=2.0,
            U_bending=3.0,
            U_shear=0.0,
            disp_u=1.0,
            disp_w=1.0,
        )
        timo = audit.shape_audit.energy_dict(
            model=audit.MODEL_TIMO,
            U_axial=2.0,
            U_bending=3.0,
            U_shear=5.0,
            disp_u=1.0,
            disp_w=1.0,
        )
        self.assertAlmostEqual(eb["axial_energy_fraction"] + eb["bending_energy_fraction"], 1.0)
        self.assertAlmostEqual(
            timo["axial_energy_fraction"]
            + timo["bending_energy_fraction"]
            + timo["shear_energy_fraction"],
            1.0,
        )

    def test_thin_limit_low_frequency_differences_are_small(self) -> None:
        deltas = audit.thin_limit_low_frequency_deltas(n_roots=4)
        self.assertEqual(len(deltas), 4)
        self.assertLess(max(deltas), 0.005)


if __name__ == "__main__":
    unittest.main()
