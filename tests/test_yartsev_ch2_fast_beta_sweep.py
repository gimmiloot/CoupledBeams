from __future__ import annotations

import math
import tempfile
import unittest
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from scripts.lib.yartsev_ch2_fast_beta_sweep import (
    MANDATORY_GLOBAL_ANCHORS_DEG,
    ExactTransferLRU,
    FastSolverValidationError,
    FastSweepSettings,
    PerformanceCounters,
    SearchInterval,
    connected_clusters,
    load_family_checkpoint,
    predict_sorted_frequencies,
    run_fast_beta_sweep,
    stable_fingerprint,
    write_family_checkpoint,
)


ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class SyntheticRoot:
    frequency_hz: float
    status: str = "accepted"


@dataclass
class SyntheticSpectrum:
    roots: list[SyntheticRoot]


def frequencies(spectrum: SyntheticSpectrum) -> list[float]:
    return [root.frequency_hz for root in spectrum.roots]


def standard_roots(beta: float) -> list[SyntheticRoot]:
    return [SyntheticRoot(100.0 * mode + 2.0 * beta) for mode in range(1, 8)]


def finalize(beta: float, roots: list[SyntheticRoot]) -> SyntheticSpectrum:
    del beta
    if any(root.status.startswith("rejected") for root in roots):
        raise ValueError("rejected local root")
    return SyntheticSpectrum(sorted(roots, key=lambda root: root.frequency_hz))


class FastBetaSweepTest(unittest.TestCase):
    def _run(
        self,
        beta_values: list[float],
        *,
        global_search=None,
        scan_interval=None,
        settings: FastSweepSettings | None = None,
    ):
        global_callback = global_search or (
            lambda beta: SyntheticSpectrum(standard_roots(beta))
        )

        def default_scan(beta: float, interval: SearchInterval):
            actual = standard_roots(beta)
            return [actual[index] for index in interval.indices]

        return run_fast_beta_sweep(
            beta_values,
            global_search=global_callback,
            scan_interval=scan_interval or default_scan,
            finalize_local=finalize,
            spectrum_frequencies=frequencies,
            root_frequency=lambda root: root.frequency_hz,
            settings=settings or FastSweepSettings(),
        )

    def test_first_beta_uses_global_search(self) -> None:
        calls: list[float] = []

        def global_search(beta: float) -> SyntheticSpectrum:
            calls.append(beta)
            return SyntheticSpectrum(standard_roots(beta))

        result = self._run([0.0, 0.5, 1.0], global_search=global_search)
        self.assertEqual(calls, [0.0])
        self.assertEqual(result.data_origins[0.0], "global_anchor")
        self.assertEqual(result.counters.global_anchor_scans, 1)

    def test_predictor_for_constant_and_linear_roots(self) -> None:
        constant = predict_sorted_frequencies(1.0, 0.5, [1.0, 2.0])
        np.testing.assert_array_equal(constant, [1.0, 2.0])
        linear = predict_sorted_frequencies(
            1.0,
            0.5,
            [2.0, 4.0],
            older_beta_deg=0.0,
            older_frequencies_hz=[1.0, 2.0],
        )
        np.testing.assert_array_equal(linear, [3.0, 6.0])

    def test_local_roots_are_sorted_after_cluster_scan(self) -> None:
        base = [100.0, 200.0, 300.0, 300.1, 500.0, 600.0, 700.0]

        def global_search(beta: float) -> SyntheticSpectrum:
            return SyntheticSpectrum([SyntheticRoot(value + beta) for value in base])

        def scan(beta: float, interval: SearchInterval):
            roots = [SyntheticRoot(base[index] + beta) for index in interval.indices]
            return list(reversed(roots))

        result = self._run([0.0, 0.5], global_search=global_search, scan_interval=scan)
        self.assertEqual(
            frequencies(result.spectra[0.5]),
            sorted(frequencies(result.spectra[0.5])),
        )

    def test_connected_cluster_is_merged(self) -> None:
        clusters = connected_clusters(
            [100.0, 200.0, 300.0, 300.1, 300.2, 600.0, 700.0],
            [101.0, 201.0, 301.0, 301.1, 301.2, 601.0, 701.0],
        )
        self.assertEqual(clusters, ((2, 3, 4),))

    def test_cluster_wrong_count_triggers_global_fallback(self) -> None:
        base = [100.0, 200.0, 300.0, 300.1, 500.0, 600.0, 700.0]

        def global_search(beta: float) -> SyntheticSpectrum:
            return SyntheticSpectrum([SyntheticRoot(value + beta) for value in base])

        def scan(beta: float, interval: SearchInterval):
            roots = [SyntheticRoot(base[index] + beta) for index in interval.indices]
            return roots[:-1] if interval.is_cluster else roots

        result = self._run([0.0, 0.5], global_search=global_search, scan_interval=scan)
        self.assertEqual(result.data_origins[0.5], "global_fallback")
        self.assertEqual(result.counters.global_fallback_scans, 1)
        self.assertIn("root count mismatch", result.fallback_reasons[0.5])

    def test_rejected_local_root_triggers_fallback(self) -> None:
        def scan(beta: float, interval: SearchInterval):
            roots = [standard_roots(beta)[index] for index in interval.indices]
            if interval.indices == (2,):
                roots[0] = SyntheticRoot(roots[0].frequency_hz, "rejected_test")
            return roots

        result = self._run([0.0, 0.5], scan_interval=scan)
        self.assertEqual(result.data_origins[0.5], "global_fallback")
        self.assertIn("rejected local root", result.fallback_reasons[0.5])

    def test_mandatory_anchor_set(self) -> None:
        self.assertEqual(
            MANDATORY_GLOBAL_ANCHORS_DEG,
            (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0),
        )

    def test_anchor_mismatch_above_tolerance_fails(self) -> None:
        def global_search(beta: float) -> SyntheticSpectrum:
            del beta
            return SyntheticSpectrum(standard_roots(0.0))

        def scan(beta: float, interval: SearchInterval):
            del beta
            return [
                SyntheticRoot(standard_roots(0.0)[index].frequency_hz + 0.01)
                for index in interval.indices
            ]

        with self.assertRaisesRegex(
            FastSolverValidationError, "FAST_SOLVER_VALIDATION_FAIL"
        ):
            self._run(
                [0.0, 15.0],
                global_search=global_search,
                scan_interval=scan,
            )

    def test_fallback_spectrum_becomes_next_predictor_seed(self) -> None:
        observed_predictions: dict[float, tuple[float, ...]] = {}

        def actual(beta: float) -> list[SyntheticRoot]:
            return [SyntheticRoot(100.0 * mode + 4.0 * beta) for mode in range(1, 8)]

        def global_search(beta: float) -> SyntheticSpectrum:
            return SyntheticSpectrum(actual(beta))

        def scan(beta: float, interval: SearchInterval):
            observed_predictions.setdefault(beta, interval.predicted_hz)
            if beta == 0.5:
                return []
            return [actual(beta)[index] for index in interval.indices]

        result = self._run(
            [0.0, 0.5, 1.0],
            global_search=global_search,
            scan_interval=scan,
            settings=FastSweepSettings(predictor_residual_limit=0.5),
        )
        self.assertEqual(result.data_origins[0.5], "global_fallback")
        self.assertAlmostEqual(observed_predictions[1.0][0], 104.0)

    def test_cache_does_not_round_omega(self) -> None:
        counters = PerformanceCounters()
        cache = ExactTransferLRU(4, counters=counters)
        point = ("point", 1)
        omega_1 = 100.0
        omega_2 = np.nextafter(omega_1, math.inf)
        cache.get("model", omega_1, point, lambda: np.eye(2))
        cache.get("model", omega_2, point, lambda: 2.0 * np.eye(2))
        self.assertEqual(counters.transfer_cache_misses, 2)
        self.assertEqual(counters.transfer_expm_evaluations, 2)

    def test_identical_arms_use_one_transfer(self) -> None:
        counters = PerformanceCounters()
        cache = ExactTransferLRU(4, counters=counters)
        arm_1 = ("same", 0.4, 0.005)
        arm_2 = ("same", 0.4, 0.005)
        cache.get("timo", 10.0, arm_1, lambda: np.eye(2))
        cache.get("timo", 10.0, arm_2, lambda: np.eye(2))
        self.assertEqual(counters.transfer_expm_evaluations, 1)
        self.assertEqual(counters.transfer_cache_hits, 1)

    def test_clamp_variants_reuse_transfer(self) -> None:
        counters = PerformanceCounters()
        cache = ExactTransferLRU(4, counters=counters)
        point = ("same", 0.4, 0.005)
        cache.get("chapter2_timo", 10.0, point, lambda: np.eye(2))
        cache.get("chapter2_timo", 10.0, point, lambda: 2.0 * np.eye(2))
        self.assertEqual(counters.transfer_expm_evaluations, 1)
        self.assertEqual(counters.transfer_cache_hits, 1)

    def test_cache_is_bounded(self) -> None:
        cache = ExactTransferLRU(2)
        for omega in (1.0, 2.0, 3.0):
            cache.get("model", omega, "point", lambda: np.eye(2))
        self.assertEqual(len(cache), 2)

    def test_counters_are_independent(self) -> None:
        counters = PerformanceCounters(
            transfer_expm_evaluations=1,
            boundary_matrix_assemblies=2,
            scaled_quality_evaluations=3,
            raw_quality_evaluations=4,
            local_root_solves=5,
            cluster_local_scans=6,
            global_anchor_scans=7,
            global_fallback_scans=8,
            predictor_failures=9,
        )
        self.assertEqual(
            [
                counters.transfer_expm_evaluations,
                counters.boundary_matrix_assemblies,
                counters.scaled_quality_evaluations,
                counters.raw_quality_evaluations,
                counters.local_root_solves,
                counters.cluster_local_scans,
                counters.global_anchor_scans,
                counters.global_fallback_scans,
                counters.predictor_failures,
            ],
            list(range(1, 10)),
        )

    def test_checkpoint_fingerprint_is_verified(self) -> None:
        rows = [{"beta_deg": 0.0, "mode": 1, "frequency_hz": 100.0}]
        fingerprint = stable_fingerprint({"material": "HMS", "root_count": 7})
        results_root = ROOT / "results"
        results_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=results_root) as temporary:
            directory = Path(temporary)
            write_family_checkpoint(
                directory,
                "family",
                rows,
                fingerprint=fingerprint,
                metadata={"status_detail": "PASS"},
            )
            self.assertIsNotNone(
                load_family_checkpoint(
                    directory,
                    "family",
                    expected_fingerprint=fingerprint,
                    expected_row_count=1,
                )
            )
            self.assertIsNone(
                load_family_checkpoint(
                    directory,
                    "family",
                    expected_fingerprint=stable_fingerprint({"material": "other"}),
                    expected_row_count=1,
                )
            )

    def test_resume_load_avoids_recalculation(self) -> None:
        rows = [{"beta_deg": 0.0, "mode": 1, "frequency_hz": 100.0}]
        fingerprint = stable_fingerprint({"family": "complete"})
        results_root = ROOT / "results"
        results_root.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=results_root) as temporary:
            directory = Path(temporary)
            write_family_checkpoint(
                directory,
                "complete",
                rows,
                fingerprint=fingerprint,
                metadata={},
            )
            calls = 0
            loaded = load_family_checkpoint(
                directory,
                "complete",
                expected_fingerprint=fingerprint,
                expected_row_count=1,
            )
            if loaded is None:
                calls += 1
            self.assertEqual(calls, 0)


if __name__ == "__main__":
    unittest.main()
