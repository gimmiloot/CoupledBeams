import csv
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_eb_epsilon_baseline_thresholds as audit,
)


def make_args(output_dir: Path, **overrides: object) -> audit.Args:
    values: dict[str, object] = {
        "epsilon_min": 0.005,
        "epsilon_max": 0.008,
        "primary_epsilon_min": 0.005,
        "primary_epsilon_max": 0.008,
        "coarse_step": 0.001,
        "k_max": 10,
        "n_spectrum_roots": 12,
        "n_candidate_roots": 20,
        "epsilon_tolerance": 1.0e-5,
        "delta_tolerance": 1.0e-6,
        "threshold": 0.10,
        "output_dir": output_dir,
        "cache_dir": output_dir / "cache",
        "reuse_cache": True,
        "force_recompute": False,
        "plot_only": False,
        "smoke": False,
    }
    values.update(overrides)
    return audit.Args(**values)  # type: ignore[arg-type]


def synthetic_evaluation(
    epsilon: float,
    deltas: list[float],
    *,
    mu: float = 0.0,
    eb_signature: tuple[str, ...] | None = None,
    timo_signature: tuple[str, ...] | None = None,
    quality: str = "resolved",
) -> audit.Evaluation:
    eb_roots = tuple(2.0 + index for index in range(12))
    timo_roots = tuple(
        eb_roots[index] / math.sqrt(1.0 + (deltas[index] if index < 10 else 0.02))
        for index in range(12)
    )
    eb_families = eb_signature or tuple("bending_EB" for _ in range(12))
    timo_families = timo_signature or tuple("bending_Timoshenko" for _ in range(12))
    references = tuple(
        audit.ReferenceRoot(value=value, family=eb_families[index], family_index=index + 1)
        for index, value in enumerate(eb_roots)
    )
    gaps_eb = audit.normalized_adjacent_gaps(eb_roots)
    gaps_timo = audit.normalized_adjacent_gaps(timo_roots)
    prefix = audit.running_prefix_maxima(deltas)
    return audit.Evaluation(
        epsilon=epsilon,
        mu=mu,
        eb_roots=eb_roots,
        timo_roots=timo_roots,
        deltas=tuple(deltas),
        prefix_deltas=prefix,
        n_true=audit.true_safe_prefix(deltas),
        eb_families=eb_families,
        timo_families=timo_families,
        family_ambiguity=(False,) * 12,
        eb_references=references,
        axial_references=tuple(audit.analytic_axial_roots(max(epsilon, 1.0e-3), 24)),
        eb_reference_abs_errors=(0.0,) * 12,
        eb_reference_rel_errors=(0.0,) * 12,
        timo_axial_abs_errors=(1.0,) * 12,
        timo_axial_rel_errors=(1.0,) * 12,
        eb_gaps=gaps_eb,
        timo_gaps=gaps_timo,
        eb_singular_values=((0.0, 1.0),) * 12,
        timo_singular_values=((0.0, 1.0),) * 12,
        possible_multiplicity=(False,) * 12,
        quality_status=quality,
        quality_reasons=() if quality == "resolved" else ("synthetic_warning",),
        cache_status="hit",
        eb_found=20,
        timo_found=20,
        eb_warnings=(),
        timo_warnings=(),
        eb_retry_attempted=False,
        timo_retry_attempted=False,
        candidate_boundary_warning=False,
        svd_recovery_attempted=False,
        svd_recovery_changed_roots=False,
        svd_recovery_call_count=0,
        point_types={"synthetic"},
    )


def fake_solve(
    args: audit.Args,
    epsilon: float,
    mu: float,
    **_kwargs: object,
) -> audit.Evaluation:
    crossing = 0.00635
    base = 0.06 + 35.0 * (epsilon - 0.005)
    deltas = [base + 0.002 * index for index in range(10)]
    result = synthetic_evaluation(epsilon, deltas, mu=mu)
    if epsilon > crossing:
        result.deltas = tuple(value + 0.012 for value in result.deltas)
        result.prefix_deltas = audit.running_prefix_maxima(result.deltas)
        result.n_true = audit.true_safe_prefix(result.deltas)
    return result


class EbEpsilonBaselineThresholdsTest(unittest.TestCase):
    def test_01_delta_n_is_running_maximum(self) -> None:
        self.assertEqual(audit.running_prefix_maxima([0.02, 0.04, 0.12, 0.07]), (0.02, 0.04, 0.12, 0.12))

    def test_02_n_true_stops_at_first_failed_mode(self) -> None:
        self.assertEqual(audit.true_safe_prefix([0.02, 0.04, 0.12, 0.07]), 2)

    def test_03_equality_is_safe(self) -> None:
        self.assertEqual(audit.true_safe_prefix([0.10] * 10), 10)

    def test_04_first_loss_uses_first_crossing(self) -> None:
        evaluations = [
            synthetic_evaluation(0.01, [0.08] * 10),
            synthetic_evaluation(0.02, [0.12] * 10),
            synthetic_evaluation(0.03, [0.07] * 10),
            synthetic_evaluation(0.04, [0.13] * 10),
        ]
        status, safe, unsafe, _note = audit.first_loss_bracket(evaluations, 1, 0.10)
        self.assertEqual(status, "resolved")
        self.assertEqual((safe.epsilon, unsafe.epsilon), (0.01, 0.02))  # type: ignore[union-attr]

    def test_05_synthetic_reentry_does_not_change_first_loss(self) -> None:
        evaluations = [
            synthetic_evaluation(0.01, [0.08] * 10),
            synthetic_evaluation(0.02, [0.12] * 10),
            synthetic_evaluation(0.03, [0.07] * 10),
        ]
        self.assertEqual(audit.reentry_count(evaluations, 1, 0.10), 1)
        _status, _safe, unsafe, _note = audit.first_loss_bracket(evaluations, 1, 0.10)
        self.assertEqual(unsafe.epsilon, 0.02)  # type: ignore[union-attr]

    def test_06_refined_bracket_contains_safe_and_unsafe_sides(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            args = make_args(Path(tmp), epsilon_tolerance=1.0e-4)
            manager = audit.EvaluationManager(args)
            with mock.patch.object(audit, "_solve_evaluation", side_effect=fake_solve):
                safe = manager.evaluate(0.005, "coarse")
                unsafe = manager.evaluate(0.006, "coarse")
                safe, unsafe, _rows, _iterations = audit.refine_first_loss(args, manager, 10, safe, unsafe)
            self.assertLessEqual(safe.prefix_deltas[9], args.threshold)
            self.assertGreater(unsafe.prefix_deltas[9], args.threshold)

    def test_07_certified_threshold_is_safe_lower_not_midpoint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            args = make_args(Path(tmp), epsilon_tolerance=1.0e-4)
            manager = audit.EvaluationManager(args)
            with mock.patch.object(audit, "_solve_evaluation", side_effect=fake_solve):
                evaluations = [manager.evaluate(value, "coarse") for value in (0.005, 0.006, 0.007, 0.008)]
                rows, _audit = audit.build_threshold_rows(args, manager, evaluations)
            resolved = next(row for row in rows if row["threshold_status"] == "resolved")
            self.assertEqual(resolved["epsilon_certified_n"], audit.finite_float(resolved, "epsilon_certified_n"))
            self.assertLess(audit.finite_float(resolved, "epsilon_certified_n"), audit.finite_float(resolved, "epsilon_star_estimate"))

    def test_08_epsilon_star_ordering_check(self) -> None:
        rows = [{"prefix_n": index, "epsilon_star_estimate": 0.05 - 0.001 * index} for index in range(1, 11)]
        self.assertTrue(audit.check_threshold_ordering(rows, 1.0e-6))
        self.assertTrue(all(row["monotonicity_status"] == "pass" for row in rows))

    def test_09_simultaneous_threshold_group(self) -> None:
        rows = [
            {"prefix_n": 1, "epsilon_star_estimate": 0.02},
            {"prefix_n": 2, "epsilon_star_estimate": 0.0200004},
            {"prefix_n": 3, "epsilon_star_estimate": 0.03},
        ]
        audit.assign_simultaneous_groups(rows, 1.0e-6)
        self.assertEqual(rows[0]["simultaneous_transition_group"], rows[1]["simultaneous_transition_group"])
        self.assertTrue(rows[0]["simultaneous_transition_group"])

    def test_10_not_reached_status(self) -> None:
        evaluations = [synthetic_evaluation(0.01, [0.02] * 10), synthetic_evaluation(0.02, [0.03] * 10)]
        status, safe, unsafe, _note = audit.first_loss_bracket(evaluations, 10, 0.10)
        self.assertEqual(status, "not_reached_in_scan_range")
        self.assertEqual(safe.epsilon, 0.02)  # type: ignore[union-attr]
        self.assertIsNone(unsafe)

    def test_10b_not_reached_row_certifies_last_safe_scan_point(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            args = make_args(Path(tmp))
            manager = audit.EvaluationManager(args)
            evaluations = [synthetic_evaluation(0.005, [0.02] * 10), synthetic_evaluation(0.008, [0.03] * 10)]
            rows, _audit_rows = audit.build_threshold_rows(args, manager, evaluations)
        row = rows[0]
        self.assertEqual(row["threshold_status"], "not_reached_in_scan_range")
        self.assertEqual(row["epsilon_certified_n"], 0.008)
        self.assertTrue(math.isnan(float(row["epsilon_star_estimate"])))

    def test_11_below_scan_status(self) -> None:
        evaluations = [synthetic_evaluation(0.01, [0.12] * 10), synthetic_evaluation(0.02, [0.13] * 10)]
        status, safe, _unsafe, _note = audit.first_loss_bracket(evaluations, 1, 0.10)
        self.assertEqual(status, "below_scan_range")
        self.assertIsNone(safe)

    def test_12_axial_root_project_normalization(self) -> None:
        epsilon = 0.02
        roots = audit.analytic_axial_roots(epsilon, 3)
        expected = np.sqrt(np.arange(1, 4) * np.pi / (2.0 * epsilon))
        np.testing.assert_allclose(roots, expected, rtol=0.0, atol=1.0e-14)

    def test_13_axial_family_classification_matches_both_numerical_theories(self) -> None:
        epsilon = 0.02
        axial = audit.analytic_axial_roots(epsilon, 8)
        roots = [float(axial[0]), 7.1, float(axial[1])]
        families, ambiguity, errors, _relative = audit.classify_timo_families(
            roots, axial, audit.normalized_adjacent_gaps(roots)
        )
        self.assertEqual(families[0], "axial")
        self.assertEqual(families[2], "axial")
        self.assertFalse(ambiguity[0])
        self.assertAlmostEqual(errors[0], 0.0)

    def test_14_fixed_fixed_bending_equation_and_lambda_normalization(self) -> None:
        roots = audit.analytic_bending_roots(4)
        for root in roots:
            self.assertAlmostEqual(math.cosh(2.0 * root) * math.cos(2.0 * root), 1.0, delta=2.0e-8)
        self.assertAlmostEqual(roots[0], 4.730040744862704 / 2.0, places=12)

    def test_15_analytic_union_contains_axial_and_bending_families(self) -> None:
        union = audit.analytic_eb_union(0.03, 12)
        self.assertEqual(len(union), 12)
        self.assertIn("axial", {item.family for item in union})
        self.assertIn("bending_EB", {item.family for item in union})
        self.assertTrue(all(right.value >= left.value for left, right in zip(union, union[1:])))

    def test_16_mu_invariance_compares_first_twelve_roots(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            args = make_args(Path(tmp))
            manager = audit.EvaluationManager(args)
            thresholds: list[dict[str, object]] = []
            with mock.patch.object(audit, "MU_AUDIT_EPSILONS", (0.006,)), mock.patch.object(audit, "_solve_evaluation", side_effect=fake_solve):
                rows = audit.mu_invariance_rows(args, manager, thresholds)
            self.assertEqual(len(rows), 4 * 12)
            self.assertTrue(all(row["status"] == "pass" for row in rows))

    def test_17_exact_axial_bending_crossing_keeps_two_reference_roots(self) -> None:
        bending = float(audit.analytic_bending_roots(6)[3])
        epsilon = math.pi / (2.0 * bending**2)
        union = audit.analytic_eb_union(epsilon, 12)
        close = [item for item in union if abs(item.value - bending) <= 1.0e-10]
        self.assertEqual(len(close), 2)
        self.assertEqual({item.family for item in close}, {"axial", "bending_EB"})
        self.assertTrue(all(item.ambiguous for item in close))

    def test_17b_near_crossing_general_recovery_keeps_two_roots(self) -> None:
        epsilon = 0.0314375
        case = audit.root_workflow.CaseSpec(mu=0.0, eta=0.0, epsilon=epsilon)
        initial = audit.root_workflow.solve_model(
            case,
            0.0,
            20,
            audit.root_workflow.MODEL_EB,
        )
        entry = audit._entry_from_root_result(initial)
        references = audit.analytic_eb_union(epsilon, 12)
        _absolute, _relative, initial_mismatches = audit._analytic_reference_errors(
            entry.roots, references, 12
        )
        recovered, sign_calls, svd_calls = audit._try_reference_guided_eb_recovery(
            epsilon,
            0.0,
            entry,
            references,
            20,
        )
        _absolute, _relative, recovered_mismatches = audit._analytic_reference_errors(
            recovered.roots, references, 12
        )
        self.assertGreater(initial_mismatches, 0)
        self.assertEqual(recovered_mismatches, 0)
        self.assertGreaterEqual(sign_calls, 1)
        self.assertGreaterEqual(svd_calls, 1)
        self.assertTrue(audit.roots_strictly_increasing(recovered.roots, 12))

    def test_18_ambiguous_family_is_marked(self) -> None:
        roots = [5.0, 5.00001, 8.0]
        families, ambiguity, _absolute, _relative = audit.classify_timo_families(
            roots, [5.0, 9.0], audit.normalized_adjacent_gaps(roots)
        )
        self.assertEqual(families[0], "axial_bending_cluster")
        self.assertTrue(ambiguity[0])

    def test_19_mode_ten_right_gap_uses_root_eleven(self) -> None:
        roots = tuple(float(index) for index in range(1, 13))
        gaps = audit.normalized_adjacent_gaps(roots)
        _left, right, local = audit.sided_gap(gaps, 9)
        expected = 2.0 * (11.0 - 10.0) / (11.0 + 10.0)
        self.assertAlmostEqual(right, expected)
        self.assertLessEqual(local, expected)

    def test_20_root_twelve_is_not_in_k10_target(self) -> None:
        evaluation = synthetic_evaluation(0.02, [0.01] * 10)
        evaluation.timo_roots = (*evaluation.timo_roots[:11], 1.0e9)
        self.assertEqual(len(evaluation.deltas), 10)
        self.assertEqual(evaluation.n_true, 10)

    def test_21_missing_root_eleven_blocks_quality(self) -> None:
        self.assertFalse(audit.roots_strictly_increasing([1.0] * 10 + [float("nan")], 11))

    def test_22_candidate_boundary_blocks_threshold_use(self) -> None:
        evaluation = synthetic_evaluation(0.02, [0.01] * 10, quality="unresolved")
        status, _safe, _unsafe, _note = audit.first_loss_bracket([evaluation], 1, 0.10)
        self.assertEqual(status, "unresolved_root_quality")

    def test_23_event_detector_finds_n_true_transition(self) -> None:
        left = synthetic_evaluation(0.01, [0.08] * 10)
        midpoint = synthetic_evaluation(0.015, [0.11] * 10)
        right = synthetic_evaluation(0.02, [0.12] * 10)
        reasons = audit.interval_event_reasons(left, midpoint, right, threshold=0.10)
        self.assertIn("N_true_transition", reasons)

    def test_24_event_detector_finds_family_reorder(self) -> None:
        left = synthetic_evaluation(0.01, [0.08] * 10)
        changed = ("axial",) + left.eb_families[1:]
        midpoint = synthetic_evaluation(0.015, [0.08] * 10, eb_signature=changed)
        right = synthetic_evaluation(0.02, [0.08] * 10, eb_signature=changed)
        reasons = audit.interval_event_reasons(left, midpoint, right, threshold=0.10)
        self.assertIn("family_reorder", reasons)

    def test_25_event_detector_finds_midpoint_only_peak(self) -> None:
        left = synthetic_evaluation(0.01, [0.08] * 10)
        midpoint = synthetic_evaluation(0.015, [0.105] * 10)
        right = synthetic_evaluation(0.02, [0.08] * 10)
        reasons = audit.interval_event_reasons(left, midpoint, right, threshold=0.10)
        self.assertIn("prefix_1_midpoint_peak", reasons)

    def test_26_refinement_is_deterministic_under_input_order(self) -> None:
        evaluations = [
            synthetic_evaluation(0.01, [0.08] * 10),
            synthetic_evaluation(0.02, [0.12] * 10),
            synthetic_evaluation(0.03, [0.13] * 10),
        ]
        first = audit.first_loss_bracket(evaluations, 1, 0.10)
        second = audit.first_loss_bracket(list(reversed(evaluations)), 1, 0.10)
        self.assertEqual(first[0], second[0])
        self.assertEqual(first[1].epsilon, second[1].epsilon)  # type: ignore[union-attr]
        self.assertEqual(first[2].epsilon, second[2].epsilon)  # type: ignore[union-attr]

    def test_27_force_recompute_verification_is_compared(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            args = make_args(Path(tmp), epsilon_tolerance=1.0e-4)
            threshold_rows = [{
                "prefix_n": 1,
                "threshold_status": "resolved",
                "epsilon_certified_n": 0.005,
                "epsilon_unsafe_upper": 0.006,
                "epsilon_star_estimate": 0.0055,
                "numerical_verification_status": "pending",
                "notes": "",
            }]
            with mock.patch.object(audit, "_solve_evaluation", side_effect=fake_solve):
                manager, _rows = audit.verify_thresholds(args, threshold_rows)
            self.assertFalse(manager.reuse_cache)
            self.assertTrue(manager.force_recompute)
            self.assertIn(threshold_rows[0]["numerical_verification_status"], {"pass", "fail"})

    def test_28_synthetic_cli_creates_all_required_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "baseline"
            with mock.patch.object(audit, "_solve_evaluation", side_effect=fake_solve), mock.patch.object(audit, "MU_AUDIT_EPSILONS", (0.006,)):
                audit.main(
                    [
                        "--epsilon-min", "0.005",
                        "--epsilon-max", "0.008",
                        "--primary-epsilon-min", "0.005",
                        "--primary-epsilon-max", "0.008",
                        "--coarse-step", "0.001",
                        "--epsilon-tolerance", "0.0001",
                        "--output-dir", str(output),
                    ]
                )
            self.assertEqual(set(audit.OUTPUT_NAMES), {path.name for path in output.iterdir() if path.name != "cache"})

    def test_29_plot_only_does_not_solve_roots(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "baseline"
            with mock.patch.object(audit, "_solve_evaluation", side_effect=fake_solve), mock.patch.object(audit, "MU_AUDIT_EPSILONS", (0.006,)):
                audit.main(["--epsilon-min", "0.005", "--epsilon-max", "0.008", "--primary-epsilon-min", "0.005", "--primary-epsilon-max", "0.008", "--coarse-step", "0.001", "--epsilon-tolerance", "0.0001", "--output-dir", str(output)])
            with mock.patch.object(audit, "_solve_evaluation", side_effect=AssertionError("root solve called")):
                result = audit.main(["--epsilon-min", "0.005", "--epsilon-max", "0.008", "--primary-epsilon-min", "0.005", "--primary-epsilon-max", "0.008", "--coarse-step", "0.001", "--epsilon-tolerance", "0.0001", "--output-dir", str(output), "--plot-only"])
            self.assertEqual(result["root_solves"], 0)

    def test_30_cli_validation_requires_root_margin(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(ValueError):
                audit.parse_args(["--n-spectrum-roots", "10", "--output-dir", tmp])

    def test_31_existing_safe_prefix_definition_agrees(self) -> None:
        deltas = [0.02, 0.04, 0.12, 0.07] + [0.01] * 6
        self.assertEqual(audit.true_safe_prefix(deltas), 2)

    def test_32_output_csv_threshold_has_exactly_ten_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "baseline"
            with mock.patch.object(audit, "_solve_evaluation", side_effect=fake_solve), mock.patch.object(audit, "MU_AUDIT_EPSILONS", (0.006,)):
                audit.main(["--epsilon-min", "0.005", "--epsilon-max", "0.008", "--primary-epsilon-min", "0.005", "--primary-epsilon-max", "0.008", "--coarse-step", "0.001", "--epsilon-tolerance", "0.0001", "--output-dir", str(output)])
            with (output / "baseline_critical_prefix_thresholds.csv").open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 10)


if __name__ == "__main__":
    unittest.main()
