import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
import sys

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_eb_timo_general_spectrum_completeness as audit,
)
from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    run_eb_epsilon_apriori_pilot as pilot,
)
from scripts.lib import general_spectrum_completeness as complete  # noqa: E402


SYNTHETIC_ROOTS = (1.003, 1.0075, *tuple(float(value) for value in range(2, 13)))


def scalar_root_matrix(roots: tuple[float, ...] = SYNTHETIC_ROOTS):
    def provider(value: float) -> np.ndarray:
        residual = float(np.prod([float(value) - root for root in roots]))
        matrix = np.eye(6, dtype=float)
        matrix[0, 0] = residual
        matrix[0, 1] = 1.0
        return matrix

    return provider


def double_nullity_matrix(roots: tuple[float, ...] = SYNTHETIC_ROOTS):
    def provider(value: float) -> np.ndarray:
        residual = float(np.prod([float(value) - root for root in roots]))
        matrix = np.eye(6, dtype=float)
        matrix[0, :] = 0.0
        matrix[1, :] = 0.0
        matrix[0, 0] = residual
        matrix[0, 2] = 1.0
        matrix[1, 1] = residual
        matrix[1, 3] = 1.0
        return matrix

    return provider


def synthetic_settings(**overrides: object) -> complete.SearchSettings:
    values: dict[str, object] = {
        "requested_roots": 12,
        "candidate_roots": 12,
        "verification_candidate_roots": 13,
        "lambda_min": 0.2,
        "lambda_max": 13.5,
        "scan_step": 0.01,
        "sigma_prefilter": 0.2,
        "adaptive_depth": 1,
    }
    values.update(overrides)
    return complete.SearchSettings(**values)  # type: ignore[arg-type]


class SyntheticCompletenessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = complete.resolve_matrix_spectrum(
            scalar_root_matrix(), settings=synthetic_settings(), model="synthetic"
        )

    def test_two_roots_inside_one_coarse_interval_are_recovered(self) -> None:
        values = self.result.values
        self.assertAlmostEqual(values[0], 1.003, places=7)
        self.assertAlmostEqual(values[1], 1.0075, places=7)

    def test_close_roots_remain_distinct(self) -> None:
        self.assertGreater(self.result.values[1] - self.result.values[0], 0.0015)

    def test_shifted_sign_scan_detects_same_sign_endpoint_pair(self) -> None:
        # Both roots lie in the unshifted [1.00, 1.01] interval, but the
        # half-step grid boundary at 1.005 separates their sign changes.
        first_two = self.result.primary.roots[:2]
        self.assertTrue(all("primary_shifted:sign_change" in root.detection_sources for root in first_two))

    def test_multiple_detection_paths_deduplicate(self) -> None:
        near_two = [root for root in self.result.primary.roots if abs(root.Lambda - 2.0) < 1.0e-5]
        self.assertEqual(len(near_two), 1)
        self.assertGreaterEqual(len(near_two[0].detection_sources), 2)

    def test_shifted_and_half_step_searches_agree(self) -> None:
        self.assertTrue(self.result.independent_agreement)
        self.assertEqual(self.result.spectrum_status, "resolved_complete")

    def test_candidate_12_not_in_k10_target_and_guards_available(self) -> None:
        self.assertEqual(len(self.result.values[:10]), 10)
        self.assertTrue(self.result.root11_available)
        self.assertTrue(self.result.root12_available)
        self.assertFalse(self.result.root12_boundary_warning)

    def test_operation_counts_are_deterministic(self) -> None:
        repeated = complete.resolve_matrix_spectrum(
            scalar_root_matrix(), settings=synthetic_settings(), model="synthetic"
        )
        self.assertEqual(
            self.result.operations.characteristic_matrix_evaluations,
            repeated.operations.characteristic_matrix_evaluations,
        )
        self.assertEqual(
            self.result.operations.full_6x6_SVD_calls,
            repeated.operations.full_6x6_SVD_calls,
        )

    def test_exact_nullity_two_preserves_multiplicity(self) -> None:
        result = complete.resolve_matrix_spectrum(
            double_nullity_matrix((1.0, 1.004, *tuple(float(value) for value in range(2, 13)))),
            settings=synthetic_settings(),
            model="synthetic_double",
        )
        first = [root for root in result.primary.roots if abs(root.Lambda - 1.0) < 1.0e-5]
        self.assertEqual(len(first), 2)
        self.assertTrue(all(root.detected_nullity == 2 for root in first))

    def test_two_continuation_tracks_can_be_marked_coalesced_without_nullity_claim(self) -> None:
        annotated = complete.annotate_coalesced_tracks(
            self.result,
            (1.003, 1.004),
            synthetic_settings(),
        )
        first = annotated.primary.roots[0]
        self.assertEqual(first.multiplicity_status, "coalesced_track_cluster")
        self.assertEqual(first.track_multiplicity, 2)
        self.assertEqual(first.detected_nullity, 1)

    def test_deep_nonzero_sigma_valley_is_rejected(self) -> None:
        roots = tuple(float(value) for value in range(2, 14))

        def provider(value: float) -> np.ndarray:
            residual = ((float(value) - 1.0) ** 2 + 1.0e-3) * float(
                np.prod([float(value) - root for root in roots])
            )
            matrix = np.eye(6)
            matrix[0, 0] = residual
            matrix[0, 1] = 1.0
            return matrix

        result = complete.resolve_matrix_spectrum(provider, settings=synthetic_settings(), model="false_valley")
        self.assertFalse(any(abs(root.Lambda - 1.0) < 0.05 for root in result.primary.roots))

    def test_boundary_minimum_is_not_accepted(self) -> None:
        roots = (0.2, *tuple(float(value) for value in range(1, 13)))
        result = complete.resolve_matrix_spectrum(
            scalar_root_matrix(roots), settings=synthetic_settings(), model="boundary"
        )
        self.assertFalse(any(abs(root.Lambda - 0.2) < 1.0e-8 for root in result.primary.roots))

    def test_candidate_boundary_warning_blocks_resolved_status(self) -> None:
        result = complete.resolve_matrix_spectrum(
            scalar_root_matrix(),
            settings=synthetic_settings(candidate_roots=20, verification_candidate_roots=21),
            model="candidate_boundary",
        )
        self.assertTrue(result.root12_available)
        self.assertTrue(result.root12_boundary_warning)
        self.assertEqual(result.spectrum_status, "unresolved")
        self.assertIn("candidate_boundary_warning", result.exclusion_reason)

    def test_independent_disagreement_blocks_resolution(self) -> None:
        primary = self.result.primary
        shortened = complete.SearchConfigurationResult(
            configuration="verification",
            scan_step=primary.scan_step / 2.0,
            grid_phase=0.5,
            candidate_root_target=13,
            lambda_upper=primary.lambda_upper,
            roots=primary.roots[1:],
            candidates=primary.candidates,
            interval_rows=primary.interval_rows,
            unresolved_intervals=(),
            operations=complete.OperationCounts(),
        )
        _rows, agreement = complete._compare_configurations(primary, shortened, synthetic_settings())
        self.assertFalse(agreement)


class StraightOracleTests(unittest.TestCase):
    def test_r1_r2_r3_oracle_contains_close_pair(self) -> None:
        for name, (epsilon, left, right) in audit.R_CASES.items():
            with self.subTest(case=name):
                geometry = complete.Geometry(epsilon, 0.0, 0.0, 0.0)
                values = complete.straight_oracle_values(complete.MODEL_TIMO, geometry, 12)
                self.assertTrue(any(abs(value - left) < 2.0e-3 for value in values))
                self.assertTrue(any(abs(value - right) < 2.0e-3 for value in values))

    def test_general_r1_matches_factorized_first12(self) -> None:
        epsilon, _left, _right = audit.R_CASES["R1"]
        geometry = complete.Geometry(epsilon, 0.0, 0.0, 0.0)
        settings = synthetic_settings(lambda_max=22.0)
        result = complete.resolve_general_spectrum(complete.MODEL_TIMO, geometry, settings=settings)
        oracle = complete.straight_oracle_values(complete.MODEL_TIMO, geometry, 12)
        np.testing.assert_allclose(result.values, oracle, rtol=0.0, atol=settings.root_match_tol)
        self.assertEqual(result.spectrum_status, "resolved_complete")

    def test_cross_family_multiplicity_is_not_deduplicated_by_oracle(self) -> None:
        epsilon = math.pi / (2.0 * (4.730040744862704 / 2.0) ** 2)
        geometry = complete.Geometry(epsilon, 0.0, 0.0, 0.0)
        values = complete.straight_oracle_values(complete.MODEL_EB, geometry, 12)
        gaps = [abs(right - left) for left, right in zip(values[:-1], values[1:])]
        self.assertLess(min(gaps), 1.0e-7)


class CacheAndIntegrationTests(unittest.TestCase):
    def test_stale_cache_version_is_rejected_and_new_cache_reused(self) -> None:
        geometry = complete.Geometry(0.02, 45.0, 0.5, 0.0)
        settings = synthetic_settings()
        result = complete.resolve_matrix_spectrum(
            scalar_root_matrix(), settings=settings, model=complete.MODEL_EB, geometry=geometry
        )
        with tempfile.TemporaryDirectory() as directory:
            cache = complete.GeneralSpectrumCache(Path(directory))
            path = cache.save(result)
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload["algorithm_version"] = "legacy"
            path.write_text(json.dumps(payload), encoding="utf-8")
            self.assertIsNone(cache.load(complete.MODEL_EB, geometry, settings))
            self.assertEqual(cache.last_load_status, "stale_cache_algorithm_version")
            cache.save(result)
            loaded = cache.load(complete.MODEL_EB, geometry, settings)
            self.assertIsNotNone(loaded)
            self.assertEqual(cache.last_load_status, "hit")

    def test_cached_coalesced_track_annotation_uses_current_seeds(self) -> None:
        geometry = complete.Geometry(0.02, 45.0, 0.5, 0.0)
        settings = synthetic_settings()
        result = complete.resolve_matrix_spectrum(
            scalar_root_matrix(), settings=settings, model=complete.MODEL_EB, geometry=geometry
        )
        annotated = complete.annotate_coalesced_tracks(result, (1.003, 1.004), settings)
        self.assertEqual(annotated.primary.roots[0].multiplicity_status, "coalesced_track_cluster")
        with tempfile.TemporaryDirectory() as directory:
            cache = complete.GeneralSpectrumCache(Path(directory))
            cache.save(annotated)
            reused = cache.resolve(
                complete.MODEL_EB,
                geometry,
                settings,
                continuation_seeds=(1.003,),
            )
        self.assertEqual(reused.cache_status, "hit")
        self.assertEqual(reused.primary.roots[0].multiplicity_status, "simple_root")
        self.assertEqual(reused.primary.roots[0].track_multiplicity, 1)

    def test_pilot_default_spectrum_method_remains_legacy(self) -> None:
        args = pilot.parse_args([])
        self.assertEqual(args.spectrum_method, pilot.SPECTRUM_METHOD_LEGACY)

    def test_auto_mode_requires_first12_and_larger_verification_margin(self) -> None:
        with self.assertRaises(ValueError):
            pilot.parse_args([
                "--spectrum-method", pilot.SPECTRUM_METHOD_AUTO,
                "--n-spectrum-roots", "11",
            ])
        with self.assertRaises(ValueError):
            pilot.parse_args([
                "--spectrum-method", pilot.SPECTRUM_METHOD_AUTO,
                "--verification-candidate-roots", "20",
            ])

    def test_incomplete_corrected_spectrum_is_not_replaced_by_legacy(self) -> None:
        case = pilot.PilotCase("X", 0.02, 45.0, 0.5, 0.0, (), "")
        args = pilot.parse_args([
            "--spectrum-method", pilot.SPECTRUM_METHOD_AUTO,
            "--output-dir", "results/_smoke/test_general_incomplete",
        ])
        entry = pilot._entry_from_values(tuple(float(index) for index in range(1, 21)), source="synthetic", cache_hit=False)
        auto = pilot.AutoSpectrumProducts(
            eb_entry=entry,
            timo_entry=entry,
            eb_status="unresolved",
            timo_status="resolved_complete",
            eb_source=complete.GENERAL_SPECTRUM_ALGORITHM_VERSION,
            timo_source=complete.GENERAL_SPECTRUM_ALGORITHM_VERSION,
            root12_boundary_warning=False,
            operations=complete.OperationCounts(),
        )
        with mock.patch.object(pilot, "auto_complete_spectrum_products", return_value=auto), mock.patch.object(
            pilot.fixed_scan, "point_products"
        ) as point_products:
            point_products.side_effect = RuntimeError("stop after corrected source selection")
            with self.assertRaises(RuntimeError):
                pilot.solve_case(
                    args,
                    case,
                    pilot.fixed_scan.RuntimeTracker(),
                    general_cache=complete.GeneralSpectrumCache(Path("results/_smoke/test_general_incomplete/cache")),
                    continuation_bank={complete.MODEL_EB: [], complete.MODEL_TIMO: []},
                )
            self.assertTrue(point_products.called)

    def test_plot_only_requires_existing_csv_and_performs_no_solver_call(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "audit"
            corrected = Path(directory) / "corrected"
            output.mkdir(parents=True)
            audit.write_csv(
                output / "general_spectrum_completeness_summary.csv",
                [{"case_id": "B01", "model": complete.MODEL_EB, "spectrum_status": "resolved_complete"}],
                audit.SUMMARY_FIELDS,
            )
            args = audit.parse_args([
                "--output-dir", str(output),
                "--corrected-pilot-dir", str(corrected),
                "--plot-only",
                "--smoke",
            ])
            with mock.patch.object(complete, "resolve_general_spectrum") as solve:
                product = audit.plot_only(args)
            solve.assert_not_called()
            self.assertEqual(product["root_calculations"], 0)


if __name__ == "__main__":
    unittest.main()
