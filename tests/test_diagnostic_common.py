import csv
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib.diagnostic_common import (  # noqa: E402
    inclusive_grid,
    number_text,
    number_token,
    write_dict_rows_csv,
)


class DiagnosticCommonTest(unittest.TestCase):
    def test_number_token_basic_cases(self) -> None:
        self.assertEqual("0", number_token(0.0))
        self.assertEqual("1p25", number_token(1.25))
        self.assertEqual("m0p5", number_token(-0.5))

    def test_number_text_finite_and_nonfinite_cases(self) -> None:
        self.assertEqual("", number_text(None))
        self.assertEqual("raw", number_text("raw"))
        self.assertEqual("2", number_text(2))
        self.assertEqual("1.25", number_text(1.25))
        self.assertEqual("", number_text(np.nan))
        self.assertEqual("nan", number_text(np.inf, nonfinite_text="nan"))

    def test_inclusive_grid_includes_both_endpoints(self) -> None:
        grid = inclusive_grid(0.0, 1.0, 0.3)
        np.testing.assert_allclose(grid, np.array([0.0, 0.3, 0.6, 0.9, 1.0]))
        self.assertEqual(0.0, float(grid[0]))
        self.assertEqual(1.0, float(grid[-1]))

    def test_inclusive_grid_rejects_nonpositive_step(self) -> None:
        with self.assertRaisesRegex(ValueError, "grid step must be positive"):
            inclusive_grid(0.0, 1.0, 0.0)
        with self.assertRaisesRegex(ValueError, "mu step must be positive"):
            inclusive_grid(0.0, 1.0, -0.1, step_name="mu step")

    def test_inclusive_grid_rejects_stop_before_start(self) -> None:
        with self.assertRaisesRegex(ValueError, "grid step stop must be greater than or equal to start"):
            inclusive_grid(1.0, 0.0, 0.1)

    def test_write_dict_rows_csv_writes_headers_and_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "nested" / "rows.csv"
            returned = write_dict_rows_csv(
                output,
                [{"a": 1, "b": "x", "ignored": "drop"}, {"a": 2, "b": "y"}],
                ["a", "b"],
            )

            self.assertEqual(output, returned)
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(["a", "b"], list(rows[0].keys()))
        self.assertEqual([{"a": "1", "b": "x"}, {"a": "2", "b": "y"}], rows)


if __name__ == "__main__":
    unittest.main()
