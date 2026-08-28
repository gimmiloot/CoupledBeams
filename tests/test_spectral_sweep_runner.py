from __future__ import annotations

import ast
import json
import math
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np

from scripts.analysis.laminated_beams import (
    validate_spectral_sweep_runner as validation_cli,
)
from scripts.lib.spectral_sweep_runner import (
    CheckpointConfig,
    ExactLRU,
    RootRecord,
    SearchInterval,
    SpectrumRecord,
    SweepCallbacks,
    SweepConfigurationError,
    SweepCounters,
    SweepSettings,
    anchor_indices,
    build_search_intervals,
    detect_spikes,
    expanded_interval,
    load_point_transaction,
    predict_sorted_roots,
    run_spectral_sweep,
    stable_fingerprint,
    write_checkpoint_manifest,
    write_point_transaction,
)


ROOT = Path(__file__).resolve().parents[1]


def _values(parameter: float) -> list[float]:
    return [100.0 * position + parameter for position in range(1, 10)]


def _close_values(parameter: float) -> list[float]:
    return [
        100.0 + parameter,
        200.0 + parameter,
        300.0 + parameter,
        400.0 + parameter,
        400.01 + parameter,
        600.0 + parameter,
        700.0 + parameter,
        800.0 + parameter,
        900.0 + parameter,
    ]


def _roots(values: list[float], **changes: object) -> tuple[RootRecord, ...]:
    return tuple(RootRecord(value=value, **changes) for value in values)


def _spectrum(
    parameter: float,
    *,
    values: list[float] | None = None,
    origin: str = "SYNTHETIC_FULL",
) -> SpectrumRecord:
    return SpectrumRecord(
        parameter=parameter,
        roots=_roots(_values(parameter) if values is None else values),
        origin=origin,
    )


def _local_from(
    value_provider,
):
    def local_search(
        parameter: float,
        interval: SearchInterval,
        verification: bool,
    ) -> list[RootRecord]:
        del verification
        actual = _roots(value_provider(parameter))
        return [actual[position - 1] for position in interval.positions]

    return local_search


def _settings(**changes: object) -> SweepSettings:
    return replace(
        SweepSettings(enable_spike_audit=False, anchor_stride=10),
        **changes,
    )


class SpectralSweepRunnerTest(unittest.TestCase):
    def test_closing_validation_contract_is_first8_plus_root9_only(self) -> None:
        self.assertEqual(validation_cli.FAST_FULL_RELATIVE_TOLERANCE, 1.0e-8)
        self.assertEqual(len(validation_cli._fixed_validation_jobs()), 24)
        self.assertEqual(len(validation_cli.RLB2D_PRODUCTION_FIELDS), 59)
        self.assertEqual(
            validation_cli.RLB2D_MISSING_REFERENCE_GROUPS,
            (
                (validation_cli.MODEL_OLD, 0.78),
                (validation_cli.MODEL_OLD, 0.80),
                (validation_cli.MODEL_RLB, 0.80),
            ),
        )

    def test_closing_missing_only_contract_excludes_completed_points(self) -> None:
        self.assertEqual(
            validation_cli.RLB2D_EXPECTED_FAST_MISSING,
            {
                validation_cli.MODEL_OLD: (0.68, 0.70, 0.72, 0.74, 0.76),
                validation_cli.MODEL_RLB: (0.74,),
            },
        )
        validation_path = (
            ROOT
            / "scripts"
            / "analysis"
            / "laminated_beams"
            / "validate_spectral_sweep_runner.py"
        )
        source = validation_path.read_text(encoding="utf-8")
        tree = ast.parse(source)
        main = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "main"
        )
        self.assertIn("run_rlb2d_missing_only", ast.unparse(main))

    def test_required_physical_fast_probe_cannot_be_relabelled_as_anchor(self) -> None:
        self.assertTrue(
            validation_cli._requires_fast_local_evidence(
                "beta", validation_cli.MODEL_RLB, 54.0
            )
        )
        self.assertTrue(
            validation_cli._requires_fast_local_evidence(
                "mu", validation_cli.MODEL_OLD, 0.78
            )
        )
        self.assertFalse(
            validation_cli._requires_fast_local_evidence(
                "beta", validation_cli.MODEL_RLB, 45.0
            )
        )
        with self.assertRaisesRegex(RuntimeError, "did not produce physical"):
            validation_cli._require_fast_local_shard(
                "beta",
                validation_cli.MODEL_RLB,
                54.0,
                {
                    "target_origin": "GLOBAL_ANCHOR",
                    "backend_counters": {
                        "local_primary_scans": 0,
                        "local_verification_scans": 0,
                        "determinant_evaluations": 0,
                    },
                },
            )

    def test_predictor_hold_and_secant_on_nonuniform_grid(self) -> None:
        hold = predict_sorted_roots(0.35, 0.20, [10.0, 20.0])
        np.testing.assert_array_equal(hold, [10.0, 20.0])

        secant = predict_sorted_roots(
            0.50,
            0.20,
            [10.0, 20.0],
            older_parameter=0.05,
            older_roots=[7.0, 14.0],
        )
        np.testing.assert_allclose(secant, [16.0, 32.0], rtol=1.0e-15, atol=0.0)

    def test_intervals_are_nonoverlapping_and_preserve_close_pair(self) -> None:
        previous = _close_values(0.0)
        predicted = _close_values(0.5)
        intervals = build_search_intervals(
            previous,
            predicted,
            older_roots=None,
            settings=_settings(cluster_relative_gap=1.0e-3),
        )

        self.assertEqual(sum(item.expected_count for item in intervals), 9)
        close = [item for item in intervals if item.positions == (4, 5)]
        self.assertEqual(len(close), 1)
        self.assertTrue(close[0].is_cluster)
        for left, right in zip(intervals, intervals[1:], strict=False):
            self.assertLessEqual(left.upper, right.lower)

    def test_inventory_is_sorted_and_root9_is_retained_bit_exactly(self) -> None:
        guard = float(np.nextafter(900.0, math.inf))
        unsorted = [guard, 500.0, 100.0, 800.0, 300.0, 200.0, 700.0, 400.0, 600.0]
        result = run_spectral_sweep(
            [0.0],
            callbacks=SweepCallbacks(
                global_search=lambda parameter: _spectrum(parameter, values=unsorted),
                local_search=_local_from(_values),
            ),
            settings=_settings(),
        )

        values = result.spectra[0.0].values
        np.testing.assert_array_equal(values[:-1], np.arange(100.0, 900.0, 100.0))
        self.assertEqual(float(values[-1]).hex(), guard.hex())
        self.assertEqual(len(values), 9)

    def test_anchor_schedules_include_endpoints_and_periodic_points(self) -> None:
        self.assertEqual(anchor_indices(91, 10), tuple(range(0, 91, 10)))
        self.assertEqual(anchor_indices(41, 5), tuple(range(0, 41, 5)))
        self.assertEqual(anchor_indices(3, 10), (0, 2))

    def test_missing_local_root_triggers_full_fallback_and_next_anchor(self) -> None:
        global_calls: list[float] = []

        def global_search(parameter: float) -> SpectrumRecord:
            global_calls.append(parameter)
            return _spectrum(parameter)

        def local_search(
            parameter: float,
            interval: SearchInterval,
            verification: bool,
        ) -> list[RootRecord]:
            del verification
            if parameter == 1.0 and interval.positions == (3,):
                return []
            actual = _roots(_values(parameter))
            return [actual[position - 1] for position in interval.positions]

        result = run_spectral_sweep(
            [0.0, 1.0, 2.0, 3.0, 4.0],
            callbacks=SweepCallbacks(global_search, local_search),
            settings=_settings(),
        )

        self.assertEqual(result.spectra[1.0].origin, "GLOBAL_FALLBACK")
        self.assertEqual(result.spectra[2.0].origin, "GLOBAL_ANCHOR")
        self.assertIn(1.0, global_calls)
        self.assertIn(2.0, global_calls)
        audit = next(item for item in result.point_audits if item.parameter == 1.0)
        self.assertIn("wrong root count", " ".join(audit.fallback_reasons))

    def test_small_gap_triggers_fallback(self) -> None:
        result = run_spectral_sweep(
            [0.0, 1.0, 2.0],
            callbacks=SweepCallbacks(
                global_search=lambda parameter: _spectrum(
                    parameter, values=_close_values(parameter)
                ),
                local_search=_local_from(_close_values),
            ),
            settings=_settings(
                cluster_relative_gap=1.0e-3,
                fallback_gap_relative=1.0e-4,
            ),
        )

        self.assertEqual(result.spectra[1.0].origin, "GLOBAL_FALLBACK")
        audit = next(item for item in result.point_audits if item.parameter == 1.0)
        self.assertIn("spectral gap", " ".join(audit.fallback_reasons))

    def test_bad_frozen_quality_residual_triggers_fallback(self) -> None:
        def local_search(
            parameter: float,
            interval: SearchInterval,
            verification: bool,
        ) -> list[RootRecord]:
            del verification
            actual = _roots(_values(parameter))
            records = [actual[position - 1] for position in interval.positions]
            if interval.positions == (4,):
                records[0] = replace(records[0], boundary_residual=1.0e-4)
            return records

        quality = lambda root: (
            root.boundary_residual <= 1.0e-9,
            "boundary residual",
        )
        result = run_spectral_sweep(
            [0.0, 1.0, 2.0],
            callbacks=SweepCallbacks(
                lambda parameter: _spectrum(parameter),
                local_search,
                quality_gate=quality,
            ),
            settings=_settings(),
        )

        self.assertEqual(result.spectra[1.0].origin, "GLOBAL_FALLBACK")
        audit = next(item for item in result.point_audits if item.parameter == 1.0)
        self.assertIn("frozen gate", " ".join(audit.fallback_reasons))

    def test_primary_verification_mismatch_triggers_fallback(self) -> None:
        def local_search(
            parameter: float,
            interval: SearchInterval,
            verification: bool,
        ) -> list[RootRecord]:
            actual = _roots(_values(parameter))
            records = [actual[position - 1] for position in interval.positions]
            if verification and interval.positions == (4,):
                records[0] = replace(records[0], value=records[0].value + 1.0e-4)
            return records

        result = run_spectral_sweep(
            [0.0, 1.0, 2.0],
            callbacks=SweepCallbacks(
                lambda parameter: _spectrum(parameter),
                local_search,
            ),
            settings=_settings(),
        )

        self.assertEqual(result.spectra[1.0].origin, "GLOBAL_FALLBACK")
        audit = next(item for item in result.point_audits if item.parameter == 1.0)
        self.assertIn("verification mismatch", " ".join(audit.fallback_reasons))

    def test_detector_disagreement_and_even_root_each_trigger_fallback(self) -> None:
        for condition in ("detector", "even"):
            with self.subTest(condition=condition):
                def local_search(
                    parameter: float,
                    interval: SearchInterval,
                    verification: bool,
                ) -> list[RootRecord]:
                    del verification
                    actual = _roots(_values(parameter))
                    records = [actual[position - 1] for position in interval.positions]
                    if interval.positions == (5,):
                        if condition == "detector":
                            records[0] = replace(records[0], detector_agreement=False)
                        else:
                            records[0] = replace(records[0], possible_even_root=True)
                    return records

                result = run_spectral_sweep(
                    [0.0, 1.0, 2.0],
                    callbacks=SweepCallbacks(
                        lambda parameter: _spectrum(parameter),
                        local_search,
                    ),
                    settings=_settings(),
                )
                self.assertEqual(result.spectra[1.0].origin, "GLOBAL_FALLBACK")

    def test_completeness_guard_failure_triggers_fallback(self) -> None:
        result = run_spectral_sweep(
            [0.0, 1.0, 2.0],
            callbacks=SweepCallbacks(
                lambda parameter: _spectrum(parameter),
                _local_from(_values),
                completeness_guard=lambda parameter, roots: (
                    parameter != 1.0,
                    "synthetic missing guard evidence",
                ),
            ),
            settings=_settings(),
        )
        self.assertEqual(result.spectra[1.0].origin, "GLOBAL_FALLBACK")

    def test_close_pair_is_searched_together_but_never_merged(self) -> None:
        observed: list[tuple[int, ...]] = []

        def local_search(
            parameter: float,
            interval: SearchInterval,
            verification: bool,
        ) -> list[RootRecord]:
            del verification
            observed.append(interval.positions)
            actual = _roots(_close_values(parameter))
            return list(
                reversed([actual[position - 1] for position in interval.positions])
            )

        result = run_spectral_sweep(
            [0.0, 1.0, 2.0],
            callbacks=SweepCallbacks(
                lambda parameter: _spectrum(
                    parameter, values=_close_values(parameter)
                ),
                local_search,
            ),
            settings=_settings(
                cluster_relative_gap=1.0e-3,
                fallback_gap_relative=1.0e-6,
            ),
        )

        roots = result.spectra[1.0].roots
        self.assertEqual(result.spectra[1.0].origin, "FAST_LOCAL")
        self.assertEqual(len(roots), 9)
        self.assertAlmostEqual(roots[3].value, 401.0)
        self.assertAlmostEqual(roots[4].value, 401.01)
        self.assertIn((4, 5), observed)

    def test_stable_fingerprint_preserves_adjacent_float_distinction(self) -> None:
        value = 1.0
        adjacent = float(np.nextafter(value, math.inf))
        self.assertEqual(
            stable_fingerprint({"a": 1, "b": [2, 3]}),
            stable_fingerprint({"b": [2, 3], "a": 1}),
        )
        self.assertNotEqual(stable_fingerprint(value), stable_fingerprint(adjacent))

    def test_point_transaction_is_atomic_and_hash_checked(self) -> None:
        with tempfile.TemporaryDirectory(dir=ROOT) as temporary:
            directory = Path(temporary)
            config = CheckpointConfig(
                directory=directory,
                sweep_id="sweep",
                model_id="model",
                fingerprint=stable_fingerprint({"contract": 1}),
            )
            stale = directory / "points" / "000000.json.tmp"
            stale.parent.mkdir(parents=True)
            stale.write_text("stale", encoding="utf-8")

            target = write_point_transaction(config, 0, _spectrum(0.0))
            self.assertTrue(target.is_file())
            self.assertFalse(stale.exists())
            loaded = load_point_transaction(config, 0, 0.0)
            self.assertIsNotNone(loaded)
            np.testing.assert_array_equal(loaded.values, _values(0.0))

            payload = json.loads(target.read_text(encoding="utf-8"))
            payload["record"]["roots"][0]["value"] = 123.0
            target.write_text(json.dumps(payload), encoding="utf-8")
            self.assertIsNone(load_point_transaction(config, 0, 0.0))

    def test_resume_and_missing_only_never_recalculate_valid_points(self) -> None:
        with tempfile.TemporaryDirectory(dir=ROOT) as temporary:
            config = CheckpointConfig(
                directory=Path(temporary),
                sweep_id="resume",
                model_id="synthetic",
                fingerprint=stable_fingerprint({"grid": [0.0, 1.0, 2.0]}),
            )
            write_point_transaction(config, 0, _spectrum(0.0))
            global_calls: list[float] = []

            def global_search(parameter: float) -> SpectrumRecord:
                global_calls.append(parameter)
                return _spectrum(parameter)

            callbacks = SweepCallbacks(global_search, _local_from(_values))
            first = run_spectral_sweep(
                [0.0, 1.0, 2.0],
                callbacks=callbacks,
                settings=_settings(),
                checkpoint=config,
                resume=True,
            )
            self.assertEqual(first.counters.points_resumed, 1)
            self.assertEqual(first.counters.points_committed, 2)
            self.assertNotIn(0.0, global_calls)

            global_calls.clear()
            second = run_spectral_sweep(
                [0.0, 1.0, 2.0],
                callbacks=callbacks,
                settings=_settings(),
                checkpoint=config,
                missing_only=True,
            )
            self.assertEqual(second.counters.points_resumed, 3)
            self.assertEqual(second.counters.points_committed, 0)
            self.assertEqual(global_calls, [])

    def test_resume_restores_forced_anchor_after_committed_fallback(self) -> None:
        with tempfile.TemporaryDirectory(dir=ROOT) as temporary:
            config = CheckpointConfig(
                directory=Path(temporary),
                sweep_id="forced-anchor-resume",
                model_id="synthetic",
                fingerprint=stable_fingerprint({"grid": [0.0, 1.0, 2.0, 3.0]}),
            )
            phase = "interrupt"
            global_calls: list[tuple[str, float]] = []
            local_calls: list[tuple[str, float]] = []

            def global_search(parameter: float) -> SpectrumRecord:
                global_calls.append((phase, parameter))
                if phase == "interrupt" and parameter == 2.0:
                    raise RuntimeError("synthetic interruption after fallback commit")
                return _spectrum(parameter)

            def local_search(
                parameter: float,
                interval: SearchInterval,
                verification: bool,
            ) -> list[RootRecord]:
                del verification
                local_calls.append((phase, parameter))
                if parameter == 1.0 and interval.positions == (3,):
                    return []
                actual = _roots(_values(parameter))
                return [actual[position - 1] for position in interval.positions]

            callbacks = SweepCallbacks(global_search, local_search)
            with self.assertRaises(RuntimeError):
                run_spectral_sweep(
                    [0.0, 1.0, 2.0, 3.0],
                    callbacks=callbacks,
                    settings=_settings(),
                    checkpoint=config,
                )

            committed_fallback = load_point_transaction(config, 1, 1.0)
            self.assertIsNotNone(committed_fallback)
            self.assertEqual(committed_fallback.origin, "GLOBAL_FALLBACK")

            phase = "resume"
            resumed = run_spectral_sweep(
                [0.0, 1.0, 2.0, 3.0],
                callbacks=callbacks,
                settings=_settings(),
                checkpoint=config,
                resume=True,
            )
            self.assertEqual(resumed.spectra[2.0].origin, "GLOBAL_ANCHOR")
            self.assertIn(("resume", 2.0), global_calls)
            self.assertNotIn(("resume", 2.0), local_calls)

    def test_overlapping_expanded_intervals_trigger_full_fallback(self) -> None:
        settings = _settings(verification_expansion_factor=10.0)
        primary = build_search_intervals(
            _values(0.0),
            _values(1.0),
            older_roots=None,
            settings=settings,
        )
        expanded = [expanded_interval(interval, settings) for interval in primary]
        self.assertTrue(
            any(left.upper > right.lower for left, right in zip(expanded, expanded[1:]))
        )

        result = run_spectral_sweep(
            [0.0, 1.0, 2.0],
            callbacks=SweepCallbacks(
                lambda parameter: _spectrum(parameter),
                _local_from(_values),
            ),
            settings=settings,
        )
        self.assertEqual(result.spectra[1.0].origin, "GLOBAL_FALLBACK")
        audit = next(item for item in result.point_audits if item.parameter == 1.0)
        self.assertIn("overlap", " ".join(audit.fallback_reasons).lower())

    def test_missing_only_recomputes_semantically_incomplete_transaction(self) -> None:
        with tempfile.TemporaryDirectory(dir=ROOT) as temporary:
            config = CheckpointConfig(
                directory=Path(temporary),
                sweep_id="repair",
                model_id="synthetic",
                fingerprint=stable_fingerprint({"grid": [0.0]}),
            )
            incomplete = SpectrumRecord(
                parameter=0.0,
                roots=_roots(_values(0.0)[:8]),
                origin="INCOMPLETE",
            )
            write_point_transaction(config, 0, incomplete)
            calls: list[float] = []

            def global_search(parameter: float) -> SpectrumRecord:
                calls.append(parameter)
                return _spectrum(parameter)

            result = run_spectral_sweep(
                [0.0],
                callbacks=SweepCallbacks(global_search, _local_from(_values)),
                settings=_settings(),
                checkpoint=config,
                missing_only=True,
            )
            self.assertEqual(calls, [0.0])
            self.assertEqual(len(result.spectra[0.0].roots), 9)

    def test_global_search_failure_is_persisted_in_checkpoint_manifest(self) -> None:
        with tempfile.TemporaryDirectory(dir=ROOT) as temporary:
            directory = Path(temporary)
            config = CheckpointConfig(
                directory=directory,
                sweep_id="failed-global",
                model_id="synthetic",
                fingerprint=stable_fingerprint({"grid": [0.0, 1.0]}),
            )

            def global_search(parameter: float) -> SpectrumRecord:
                if parameter == 1.0:
                    raise RuntimeError("synthetic global-search failure")
                return _spectrum(parameter)

            with self.assertRaises(RuntimeError):
                run_spectral_sweep(
                    [0.0, 1.0],
                    callbacks=SweepCallbacks(global_search, _local_from(_values)),
                    settings=_settings(),
                    checkpoint=config,
                )

            manifest = json.loads(
                (directory / "checkpoint_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["completed_points"], [0.0])
            self.assertEqual(manifest["missing_points"], [1.0])
            self.assertEqual(manifest["failed_points"], [1.0])

    def test_resume_rejects_checkpoint_fingerprint_mismatch_before_callbacks(self) -> None:
        with tempfile.TemporaryDirectory(dir=ROOT) as temporary:
            directory = Path(temporary)
            original = CheckpointConfig(
                directory=directory,
                sweep_id="fingerprint",
                model_id="synthetic",
                fingerprint=stable_fingerprint({"contract": "original"}),
            )
            callbacks = SweepCallbacks(
                lambda parameter: _spectrum(parameter),
                _local_from(_values),
            )
            run_spectral_sweep(
                [0.0, 1.0],
                callbacks=callbacks,
                settings=_settings(),
                checkpoint=original,
            )

            changed = replace(
                original,
                fingerprint=stable_fingerprint({"contract": "changed"}),
            )
            callback_calls: list[float] = []

            def forbidden_global(parameter: float) -> SpectrumRecord:
                callback_calls.append(parameter)
                return _spectrum(parameter)

            with self.assertRaises(SweepConfigurationError):
                run_spectral_sweep(
                    [0.0, 1.0],
                    callbacks=SweepCallbacks(forbidden_global, _local_from(_values)),
                    settings=_settings(),
                    checkpoint=changed,
                    resume=True,
                )
            self.assertEqual(callback_calls, [])

    def test_checkpoint_manifest_contains_required_resume_fields(self) -> None:
        with tempfile.TemporaryDirectory(dir=ROOT) as temporary:
            config = CheckpointConfig(
                directory=Path(temporary),
                sweep_id="manifest",
                model_id="synthetic",
                fingerprint=stable_fingerprint("manifest"),
            )
            counters = SweepCounters(
                elapsed_time=4.0,
                global_full_scans=2,
                fallback_count=1,
                spike_audit_count=3,
            )
            target = write_checkpoint_manifest(
                config,
                [0.0, 1.0],
                {0.0: _spectrum(0.0)},
                counters,
                failed_points=[1.0],
                qualified_points=[0.0],
                estimated_remaining_time=2.0,
            )
            payload = json.loads(target.read_text(encoding="utf-8"))
            required = {
                "completed_points",
                "missing_points",
                "failed_points",
                "qualified_points",
                "last_completed_parameter",
                "elapsed_time",
                "estimated_remaining_time",
                "full_scan_count",
                "fallback_count",
                "spike_audit_count",
            }
            self.assertTrue(required.issubset(payload))
            self.assertEqual(payload["completed_points"], [0.0])
            self.assertEqual(payload["missing_points"], [1.0])

    def test_reproduced_physical_kink_is_not_smoothed(self) -> None:
        def kink_values(parameter: float) -> list[float]:
            values = _values(parameter)
            if parameter == 1.0:
                values[0] = 101.9
            return values

        result = run_spectral_sweep(
            [0.0, 1.0, 2.0],
            callbacks=SweepCallbacks(
                lambda parameter: _spectrum(
                    parameter, values=kink_values(parameter)
                ),
                _local_from(kink_values),
            ),
            settings=replace(_settings(), enable_spike_audit=True),
        )

        outcomes = [item for item in result.spike_audits if item.position == 1]
        self.assertEqual([item.outcome for item in outcomes], ["REPRODUCED_BY_FULL_SCAN"])
        self.assertEqual(result.spectra[1.0].roots[0].value, 101.9)
        self.assertNotEqual(result.spectra[1.0].roots[0].value, 101.0)

    def test_false_fast_spike_is_replaced_only_by_full_scan(self) -> None:
        def local_values(parameter: float) -> list[float]:
            values = _values(parameter)
            if parameter == 1.0:
                values[0] = 101.9
            return values

        result = run_spectral_sweep(
            [0.0, 1.0, 2.0],
            callbacks=SweepCallbacks(
                lambda parameter: _spectrum(parameter),
                _local_from(local_values),
            ),
            settings=replace(_settings(), enable_spike_audit=True),
        )

        outcomes = [item for item in result.spike_audits if item.position == 1]
        self.assertEqual([item.outcome for item in outcomes], ["FAST_LOCATOR_CORRECTED"])
        self.assertEqual(outcomes[0].fast_value, 101.9)
        self.assertEqual(outcomes[0].full_value, 101.0)
        self.assertEqual(result.spectra[1.0].roots[0].value, 101.0)

    def test_spike_audit_validates_all_three_full_scans_or_is_unresolved(self) -> None:
        call_counts = {0.0: 0, 1.0: 0, 2.0: 0}

        def kink_values(parameter: float) -> list[float]:
            values = _values(parameter)
            if parameter == 1.0:
                values[0] = 101.9
            return values

        def global_search(parameter: float) -> SpectrumRecord:
            call_counts[parameter] += 1
            record = _spectrum(parameter, values=kink_values(parameter))
            if parameter == 0.0 and call_counts[parameter] == 2:
                damaged = (replace(record.roots[0], accepted=False), *record.roots[1:])
                return replace(record, roots=damaged)
            return record

        result = run_spectral_sweep(
            [0.0, 1.0, 2.0],
            callbacks=SweepCallbacks(global_search, _local_from(kink_values)),
            settings=replace(_settings(), enable_spike_audit=True),
        )

        outcomes = [item for item in result.spike_audits if item.position == 1]
        self.assertEqual([item.outcome for item in outcomes], ["UNRESOLVED"])
        self.assertEqual(call_counts, {0.0: 2, 1.0: 1, 2.0: 2})
        self.assertEqual(result.spectra[1.0].roots[0].value, 101.9)
        self.assertEqual(result.status, "FAIL")

    def test_spike_detector_only_marks_data_and_does_not_mutate_it(self) -> None:
        spectra = {
            0.0: _spectrum(0.0),
            1.0: _spectrum(1.0, values=[101.9, *_values(1.0)[1:]]),
            2.0: _spectrum(2.0),
        }
        before = spectra[1.0].values.copy()
        suspicious = detect_spikes(
            [0.0, 1.0, 2.0],
            spectra,
            replace(_settings(), enable_spike_audit=True),
        )
        self.assertIn(1, [item.position for item in suspicious])
        np.testing.assert_array_equal(spectra[1.0].values, before)

    def test_exact_bounded_cache_matches_uncached_and_evicts(self) -> None:
        cache: ExactLRU[np.ndarray] = ExactLRU(max_entries=2, max_bytes=4096)
        builds = 0

        def build(value: float) -> np.ndarray:
            nonlocal builds
            builds += 1
            return np.asarray([[value, value**2]], dtype=float)

        omega = 10.0
        uncached = build(omega)
        cached_1 = cache.get_or_build(("arm", omega), lambda: build(omega))
        cached_2 = cache.get_or_build(("arm", omega), lambda: build(omega))
        np.testing.assert_allclose(cached_1, uncached, rtol=0.0, atol=0.0)
        np.testing.assert_allclose(cached_2, uncached, rtol=0.0, atol=0.0)
        self.assertIs(cached_1, cached_2)

        adjacent = float(np.nextafter(omega, math.inf))
        cache.get_or_build(("arm", adjacent), lambda: build(adjacent))
        cache.get_or_build(("arm", 12.0), lambda: build(12.0))
        self.assertEqual(len(cache), 2)
        self.assertEqual(cache.evictions, 1)
        self.assertNotIn(("arm", omega), tuple(cache))
        self.assertEqual(cache.diagnostics()["peak_entries"], 2)

    def test_exact_cache_is_independently_bounded_by_payload_bytes(self) -> None:
        cache: ExactLRU[np.ndarray] = ExactLRU(max_entries=10, max_bytes=24)
        first = np.asarray([1.0, 2.0], dtype=np.float64)
        second = np.asarray([3.0, 4.0], dtype=np.float64)
        cache["first"] = first
        cache["second"] = second

        self.assertEqual(tuple(cache), ("second",))
        self.assertLessEqual(cache.current_bytes, cache.max_bytes)
        self.assertEqual(cache.evictions, 1)
        self.assertEqual(cache.diagnostics()["peak_bytes"], second.nbytes)

    def test_runner_rejects_non_strict_or_nonfinite_parameter_grids(self) -> None:
        callbacks = SweepCallbacks(
            lambda parameter: _spectrum(parameter),
            _local_from(_values),
        )
        invalid_grids = (
            [],
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, math.nan],
            [0.0, math.inf],
            [[0.0, 1.0]],
        )
        for grid in invalid_grids:
            with self.subTest(grid=grid):
                with self.assertRaises(SweepConfigurationError):
                    run_spectral_sweep(
                        grid,
                        callbacks=callbacks,
                        settings=_settings(),
                    )

    def test_generic_runner_imports_no_production_physics(self) -> None:
        module_path = ROOT / "scripts" / "lib" / "spectral_sweep_runner.py"
        source = module_path.read_text(encoding="utf-8")
        tree = ast.parse(source)
        imported_roots: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported_roots.update(alias.name.split(".", 1)[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported_roots.add(node.module.split(".", 1)[0])

        allowed = {
            "__future__",
            "collections",
            "dataclasses",
            "hashlib",
            "json",
            "math",
            "numpy",
            "os",
            "pathlib",
            "time",
            "typing",
        }
        self.assertLessEqual(imported_roots, allowed)
        self.assertTrue(
            imported_roots.isdisjoint(
                {"scripts", "src", "scipy", "multiprocessing", "concurrent"}
            )
        )


if __name__ == "__main__":
    unittest.main()
