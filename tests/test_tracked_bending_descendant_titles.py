import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.lib.tracked_bending_descendant_shapes import infer_bending_descendant_title_label  # noqa: E402


class TrackedBendingDescendantTitleTest(unittest.TestCase):
    def test_bending_descendant_number_label(self):
        cases = {
            "bending_desc_01": "потомок 1-й изгибной ветви",
            "bending_desc_02": "потомок 2-й изгибной ветви",
            "bending_desc_04": "потомок 4-й изгибной ветви",
            "bending_desc_05": "потомок 5-й изгибной ветви",
            "bending_desc_06": "потомок 6-й изгибной ветви",
            "bending_desc_12": "потомок 12-й изгибной ветви",
        }
        for branch_id, expected in cases.items():
            with self.subTest(branch_id=branch_id):
                self.assertEqual(infer_bending_descendant_title_label(branch_id), expected)

    def test_fallback_label(self):
        self.assertEqual(infer_bending_descendant_title_label("foo"), "ветвь foo")


if __name__ == "__main__":
    unittest.main()
