import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts"
    / "analysis"
    / "thickness_mismatch"
    / "audits"
    / "compare_full_spectrum_analytic_vs_3d_fem.py"
)


def _load_module():
    spec = importlib.util.spec_from_file_location("compare_full_spectrum_analytic_vs_3d_fem", SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


MODULE = _load_module()


class FullSpectrumAnalyticVs3DFemHelpersTest(unittest.TestCase):
    def test_analytic_union_preserves_sorted_duplicate_subsystem_roots(self) -> None:
        rows = MODULE.build_analytic_full_spectrum_union(
            in_plane_roots=[1.0, 2.0, 2.0],
            out_of_plane_roots=[1.5, 2.0],
        )

        self.assertEqual([row["analytic_full_index"] for row in rows], [1, 2, 3, 4, 5])
        self.assertEqual([row["Lambda"] for row in rows], [1.0, 1.5, 2.0, 2.0, 2.0])
        self.assertEqual(
            [(row["analytic_subsystem"], row["analytic_subsystem_index"]) for row in rows],
            [
                ("in_plane", 1),
                ("out_of_plane", 1),
                ("in_plane", 2),
                ("in_plane", 3),
                ("out_of_plane", 2),
            ],
        )

    def test_3d_mode_classifier_uses_xy_vs_z_displacement_energy(self) -> None:
        label, in_plane, out_of_plane = MODULE.classify_3d_displacement_mode(
            np.asarray([[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]])
        )
        self.assertEqual(label, "mostly_in_plane")
        self.assertAlmostEqual(in_plane, 1.0)
        self.assertAlmostEqual(out_of_plane, 0.0)

        label, in_plane, out_of_plane = MODULE.classify_3d_displacement_mode(
            np.asarray([[0.0, 0.0, 1.0], [0.0, 0.0, 2.0]])
        )
        self.assertEqual(label, "mostly_out_of_plane")
        self.assertAlmostEqual(in_plane, 0.0)
        self.assertAlmostEqual(out_of_plane, 1.0)

        label, in_plane, out_of_plane = MODULE.classify_3d_displacement_mode(np.asarray([[1.0, 0.0, 1.0]]))
        self.assertEqual(label, "mixed_or_unclassified")
        self.assertAlmostEqual(in_plane, 0.5)
        self.assertAlmostEqual(out_of_plane, 0.5)

    def test_nearest_frequency_matching_keeps_duplicate_analytic_roots(self) -> None:
        analytic_rows = [
            {
                "case_id": "case",
                "analytic_full_index": 1,
                "analytic_subsystem": "in_plane",
                "analytic_subsystem_index": 1,
                "Lambda": 2.0,
            },
            {
                "case_id": "case",
                "analytic_full_index": 2,
                "analytic_subsystem": "out_of_plane",
                "analytic_subsystem_index": 1,
                "Lambda": 2.0,
            },
            {
                "case_id": "case",
                "analytic_full_index": 3,
                "analytic_subsystem": "in_plane",
                "analytic_subsystem_index": 2,
                "Lambda": 3.0,
            },
        ]
        fem_rows = [
            {"case_id": "case", "fem_mode_index": 1, "Lambda_3d_fem": 2.001},
            {"case_id": "case", "fem_mode_index": 2, "Lambda_3d_fem": 2.002},
        ]

        rows = MODULE.nearest_frequency_match_rows(analytic_rows, fem_rows, ambiguity_rel_tol=2e-3)

        self.assertEqual([row["nearest_analytic_full_index"] for row in rows], [1, 2])
        self.assertEqual([row["ambiguity_count_within_tolerance"] for row in rows], [2, 2])
        self.assertTrue(all("ambiguous" in row["warning"] for row in rows))


if __name__ == "__main__":
    unittest.main()
