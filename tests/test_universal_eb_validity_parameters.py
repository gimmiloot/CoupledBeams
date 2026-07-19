import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    analyze_universal_eb_validity_parameters_stage1 as universal,
)


class UniversalEbValidityParameterTest(unittest.TestCase):
    def test_pi_eb_is_invariant_under_eigenvector_normalization(self) -> None:
        row = {
            "epsilon_0": "0.0025",
            "beta_deg": "45.0",
            "eta": "0.0",
            "mu": "0.0",
            "eb_sorted_index": "1",
            "Lambda_EB": "3.9265715971118738",
        }
        result = universal.eb_result_from_mode_row(row, n_shape_points=401)
        base = universal.eb_pi_metrics(result, amplitude_scale=1.0)
        scaled = universal.eb_pi_metrics(result, amplitude_scale=7.5)

        for key in ("Pi_shear", "Pi_rotary", "Pi_EB"):
            self.assertAlmostEqual(base[key], scaled[key], delta=1.0e-12)
            self.assertGreater(base[key], 0.0)

    def test_pi_eb_rod_fractions_sum_to_one(self) -> None:
        row = {
            "epsilon_0": "0.01",
            "beta_deg": "45.0",
            "eta": "0.0",
            "mu": "0.3",
            "eb_sorted_index": "4",
            "Lambda_EB": "6.987337734505156",
        }
        result = universal.eb_result_from_mode_row(row, n_shape_points=401)
        metrics = universal.eb_pi_metrics(result)
        self.assertAlmostEqual(
            metrics["Pi_shear_rod1_fraction"] + metrics["Pi_shear_rod2_fraction"],
            1.0,
            delta=1.0e-12,
        )
        self.assertAlmostEqual(
            metrics["Pi_rotary_rod1_fraction"] + metrics["Pi_rotary_rod2_fraction"],
            1.0,
            delta=1.0e-12,
        )


if __name__ == "__main__":
    unittest.main()
