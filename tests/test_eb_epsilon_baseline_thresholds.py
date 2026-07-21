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
    factor_roots = tuple(
        audit.factorized.FamilyRoot(
            value=value,
            family="bending_Timoshenko",
            family_index=index + 1,
            detection_source="synthetic",
            block_sigma_min=0.0,
        )
        for index, value in enumerate(timo_roots)
    )
    factor_spectrum = audit.factorized.FactorizedSpectrum(
        epsilon=epsilon,
        mu=mu,
        requested_roots=12,
        roots=factor_roots,
        bending_search=audit.factorized.BendingRootSearch(
            roots=timo_roots,
            sources=("synthetic",) * 12,
            sigma_minima=(0.0,) * 12,
            scan_upper=20.0,
            determinant_evaluations=0,
            svd_refinements=0,
            warnings=(),
        ),
        cache_status="hit",
    )
    return audit.Evaluation(
        epsilon=epsilon,
        mu=mu,
        eb_roots=eb_roots,
        timo_roots=timo_roots,
        factorized_timo=factor_spectrum,
        raw_eb_roots=eb_roots,
        raw_timo_roots=timo_roots,
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
        general_cache_status="hit",
        eb_found=20,
        timo_found=20,
        eb_warnings=(),
        timo_warnings=(),
        raw_eb_warnings=(),
        raw_timo_warnings=(),
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
            self.assertEqual(set(audit.OUTPUT_NAMES), {path.name for path in output.iterdir() if path.is_file()})

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

    def test_33_exact_timoshenko_block_indices_and_zero_off_block(self) -> None:
        matrix, bending, axial, _warnings = audit.factorized.straight_blocks(4.0, 0.3, 0.02)
        self.assertEqual(bending.shape, (4, 4))
        self.assertEqual(axial.shape, (2, 2))
        self.assertEqual(audit.factorized.BENDING_ROWS, (0, 2, 3, 4))
        self.assertEqual(audit.factorized.BENDING_COLUMNS, (0, 1, 3, 4))
        self.assertEqual(audit.factorized.AXIAL_ROWS, (1, 5))
        self.assertEqual(audit.factorized.AXIAL_COLUMNS, (2, 5))
        self.assertLessEqual(audit.factorized.off_block_max_abs(matrix), 1.0e-14)

    def test_34_exact_axial_root_is_singular_in_axial_and_full_blocks(self) -> None:
        epsilon = 0.02
        root = float(audit.factorized.exact_axial_roots(epsilon, 1)[0])
        axial_sigma, _ = audit.factorized.axial_block_singular_values(root, 0.0, epsilon)
        full_sigma, _ = audit.factorized.full_matrix_singular_values(root, 0.0, epsilon)
        self.assertLessEqual(axial_sigma, audit.factorized.AXIAL_BLOCK_SVD_ACCEPT)
        self.assertLessEqual(full_sigma, audit.factorized.FULL_MATRIX_SVD_ACCEPT)

    def assert_regression_pair(self, epsilon: float, bending: float, axial: float) -> None:
        spectrum = audit.factorized.factorized_straight_spectrum(epsilon, 0.0, 12)
        self.assertEqual(len(spectrum.roots), 12)
        bend = min(
            (item for item in spectrum.roots if item.family == "bending_Timoshenko"),
            key=lambda item: abs(item.value - bending),
        )
        axial_item = min(
            (item for item in spectrum.roots if item.family == "axial"),
            key=lambda item: abs(item.value - axial),
        )
        self.assertLess(abs(bend.value - bending), 0.05)
        self.assertLess(abs(axial_item.value - axial), 0.05)
        self.assertNotEqual(spectrum.roots.index(bend), spectrum.roots.index(axial_item))
        self.assertLessEqual(
            audit.factorized.full_matrix_singular_values(bend.value, 0.0, epsilon)[0],
            audit.factorized.FULL_MATRIX_SVD_ACCEPT,
        )
        self.assertLessEqual(
            audit.factorized.full_matrix_singular_values(axial_item.value, 0.0, epsilon)[0],
            audit.factorized.FULL_MATRIX_SVD_ACCEPT,
        )

    def test_35_regression_R1_keeps_bending_and_axial_roots(self) -> None:
        self.assert_regression_pair(0.0228017578125, 8.2932, 8.3000)

    def test_36_regression_R2_keeps_bending_and_axial_roots(self) -> None:
        self.assert_regression_pair(0.0159306640625, 9.9287, 9.9299)

    def test_37_regression_R3_keeps_bending_and_axial_roots(self) -> None:
        self.assert_regression_pair(0.015625, 14.1773, 14.1796)

    def test_38_regression_R1_survives_plus_minus_epsilon_perturbation(self) -> None:
        for offset in (-1.0e-6, 1.0e-6):
            with self.subTest(offset=offset):
                self.assert_regression_pair(0.0228017578125 + offset, 8.2932, 8.3000)

    def test_39_regression_R2_survives_plus_minus_epsilon_perturbation(self) -> None:
        for offset in (-1.0e-6, 1.0e-6):
            with self.subTest(offset=offset):
                self.assert_regression_pair(0.0159306640625 + offset, 9.9287, 9.9299)

    def test_40_regression_R3_survives_plus_minus_epsilon_perturbation(self) -> None:
        for offset in (-1.0e-6, 1.0e-6):
            with self.subTest(offset=offset):
                self.assert_regression_pair(0.015625 + offset, 14.1773, 14.1796)

    def test_41_cross_family_duplicates_are_not_deduplicated(self) -> None:
        epsilon = audit.factorized.brentq(
            lambda value: audit.factorized.bending_block_det(
                float(audit.factorized.exact_axial_roots(value, 1)[0]),
                0.0,
                value,
            ),
            0.01592,
            0.01594,
        )
        bending = float(audit.factorized.exact_axial_roots(epsilon, 1)[0])
        spectrum = audit.factorized.factorized_straight_spectrum(epsilon, 0.0, 12)
        close = [item for item in spectrum.roots if abs(item.value - bending) <= 2.0e-6]
        self.assertEqual(len(close), 2)
        self.assertEqual({item.family for item in close}, {"axial", "bending_Timoshenko"})

    def test_42_general_matching_preserves_multiplicity_and_marks_missing(self) -> None:
        spectrum = audit.factorized.factorized_straight_spectrum(0.0228017578125, 0.0, 12)
        general, _warnings = audit.factorized.TIMO.timo_sorted_roots(
            0.0, 0.0, 0.0228017578125, 12, eta=0.0
        )
        matches = audit.match_general_roots(spectrum.roots, general)
        self.assertGreaterEqual(sum(match[0] is None for match in matches), 2)

    def test_43_cache_key_changes_with_versioned_settings(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache = audit.factorized.FactorizedSpectrumCache(Path(tmp))
            first = cache.path_for(0.02, 0.0, 12)
            second = cache.path_for(0.02, 0.0, 13)
            candidate_changed = cache.path_for(0.02, 0.0, 12, 20)
            self.assertNotEqual(first, second)
            self.assertNotEqual(first, candidate_changed)
            self.assertIn("factorized_straight_spectrum_v2", audit.factorized.ALGORITHM_VERSION)

    def test_44_stale_cache_algorithm_version_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache = audit.factorized.FactorizedSpectrumCache(Path(tmp))
            path = cache.path_for(0.02, 0.0, 12)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text('{"algorithm_version":"factorized_straight_spectrum_v1"}\n', encoding="utf-8")
            result = cache.get(0.02, 0.0, 12)
            self.assertEqual(result.cache_status, "stale_cache_algorithm_version")
            self.assertEqual(result.algorithm_version, audit.factorized.ALGORITHM_VERSION)

    def test_45_factorized_cache_reuse_returns_identical_spectrum(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache = audit.factorized.FactorizedSpectrumCache(Path(tmp))
            first = cache.get(0.02, 0.0, 12)
            second = cache.get(0.02, 0.0, 12)
            self.assertEqual(first.values, second.values)
            self.assertEqual(first.cache_status, "miss")
            self.assertEqual(second.cache_status, "hit")

    def test_46_old_general_cache_file_cannot_contaminate_factorized_cache(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "cache").mkdir()
            (root / "cache" / "legacy.json").write_text("{}\n", encoding="utf-8")
            cache = audit.factorized.FactorizedSpectrumCache(root / "cache_factorized_straight_spectrum_v2")
            result = cache.get(0.02, 0.0, 12)
            self.assertEqual(result.cache_status, "miss")

    def test_47_conservative_floor_is_not_rounding(self) -> None:
        self.assertEqual(audit.conservative_decimal_floor(0.02280999, 5), "0.02280")
        self.assertEqual(audit.conservative_decimal_floor(0.02280999, 4), "0.0228")
        self.assertEqual(audit.presentation_round(0.02280999, 5), "0.02281")

    def test_48_factorized_mu_invariance_for_first_twelve_roots(self) -> None:
        baseline = audit.factorized.factorized_straight_spectrum(0.02, 0.0, 12)
        for mu in (0.3, 0.7, 0.9):
            with self.subTest(mu=mu):
                current = audit.factorized.factorized_straight_spectrum(0.02, mu, 12)
                np.testing.assert_allclose(current.values, baseline.values, rtol=1.0e-7, atol=1.0e-6)
                self.assertEqual(current.families, baseline.families)
                self.assertEqual(
                    tuple(item.multiplicity_group for item in current.roots),
                    tuple(item.multiplicity_group for item in baseline.roots),
                )

    def test_49_off_block_entries_vanish_over_required_parameter_grid(self) -> None:
        maximum = 0.0
        for mu in (0.0, 0.3, 0.7, 0.9):
            for epsilon in (0.01, 0.02, 0.05):
                for Lambda in (2.75, 7.37, 11.11):
                    with self.subTest(mu=mu, epsilon=epsilon, Lambda=Lambda):
                        matrix, bending, axial, _warnings = audit.factorized.straight_blocks(
                            Lambda, mu, epsilon
                        )
                        self.assertEqual(bending.shape, (4, 4))
                        self.assertEqual(axial.shape, (2, 2))
                        maximum = max(maximum, audit.factorized.off_block_max_abs(matrix))
        self.assertLessEqual(maximum, audit.factorized.BLOCK_OFFDIAGONAL_ABS_TOL)

    def test_50_block_determinant_product_matches_full_determinant(self) -> None:
        for Lambda, mu, epsilon in ((4.123, 0.0, 0.02), (7.37, 0.3, 0.01), (11.11, 0.9, 0.05)):
            with self.subTest(Lambda=Lambda, mu=mu, epsilon=epsilon):
                matrix, bending, axial, _warnings = audit.factorized.straight_blocks(
                    Lambda, mu, epsilon
                )
                full = float(np.linalg.det(matrix))
                product = float(np.linalg.det(bending) * np.linalg.det(axial))
                self.assertTrue(np.isclose(full, -product, rtol=2.0e-12, atol=1.0e-14))

    def test_51_exact_axial_roots_agree_for_eb_and_timoshenko(self) -> None:
        for epsilon in (0.01, 0.02, 0.05):
            np.testing.assert_allclose(
                audit.analytic_axial_roots(epsilon, 12),
                audit.factorized.exact_axial_roots(epsilon, 12),
                rtol=0.0,
                atol=1.0e-14,
            )

    def test_52_axial_generation_is_sufficient_for_first_twelve_union(self) -> None:
        spectrum = audit.factorized.factorized_straight_spectrum(0.05, 0.0, 12)
        axial_records = [item for item in spectrum.roots if item.family == "axial"]
        exact = audit.factorized.exact_axial_roots(0.05, 12)
        self.assertTrue(axial_records)
        self.assertTrue(all(item.family_index <= len(exact) for item in axial_records))
        self.assertEqual(len(spectrum.roots), 12)

    def test_53_same_family_candidates_are_deduplicated(self) -> None:
        candidates = [
            (5.0, "sign_change", 1.0e-10),
            (5.0 + 0.5 * audit.factorized.BENDING_DEDUP_TOL, "svd_minimum", 1.0e-12),
            (7.0, "sign_change", 1.0e-11),
        ]
        deduplicated = audit.factorized._deduplicate_candidates(  # noqa: SLF001
            candidates, audit.factorized.BENDING_DEDUP_TOL
        )
        self.assertEqual(len(deduplicated), 2)
        self.assertIn("sign_change", deduplicated[0][1])
        self.assertIn("svd_minimum", deduplicated[0][1])

    def test_54_missing_required_axial_root_is_detected_and_locally_verified(self) -> None:
        epsilon = 0.0228017578125
        spectrum = audit.factorized.factorized_straight_spectrum(epsilon, 0.0, 12)
        general, _warnings = audit.factorized.TIMO.timo_sorted_roots(
            0.0, 0.0, epsilon, 12, eta=0.0
        )
        matches = audit.match_general_roots(spectrum.roots, general)
        missing_axial = [
            index
            for index, (item, match) in enumerate(zip(spectrum.roots, matches))
            if item.family == "axial" and match[0] is None
        ]
        self.assertTrue(missing_axial)
        evaluation = synthetic_evaluation(epsilon, [0.01] * 10)
        evaluation.factorized_timo = spectrum
        evaluation.timo_roots = spectrum.values
        evaluation.raw_timo_roots = tuple(float(value) for value in general)
        _root, sigma, confirmed = audit.local_full_svd_confirmation(
            evaluation, missing_axial[0]
        )
        self.assertTrue(confirmed)
        self.assertLessEqual(sigma, audit.factorized.FULL_MATRIX_SVD_ACCEPT)

    def test_55_factorized_csv_exposes_required_quality_status(self) -> None:
        epsilon = 0.0228017578125
        spectrum = audit.factorized.factorized_straight_spectrum(epsilon, 0.0, 12)
        general, _warnings = audit.factorized.TIMO.timo_sorted_roots(
            0.0, 0.0, epsilon, 12, eta=0.0
        )
        evaluation = synthetic_evaluation(epsilon, [0.01] * 10)
        evaluation.factorized_timo = spectrum
        evaluation.timo_roots = spectrum.values
        evaluation.raw_timo_roots = tuple(float(value) for value in general)
        rows = audit.factorized_spectrum_rows([evaluation])
        flagged = [
            row
            for row in rows
            if row["model"] == "Timoshenko"
            and row["family"] == "axial"
            and row["raw_general_missing"]
        ]
        self.assertTrue(flagged)
        self.assertTrue(
            all(
                row["quality_status"]
                == "missing_required_axial_root_in_general_scan"
                for row in flagged
            )
        )

    def test_56_legacy_comparison_output_is_generated_when_legacy_exists(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "baseline"
            legacy_dir = output / "legacy_pre_factorized_root_fix"
            legacy_dir.mkdir(parents=True)
            legacy_path = legacy_dir / "baseline_critical_prefix_thresholds.csv"
            audit.write_csv(
                legacy_path,
                [{
                    "prefix_n": 1,
                    "threshold_status": "resolved",
                    "epsilon_certified_n": 0.049,
                    "epsilon_star_estimate": 0.05,
                    "triggering_sorted_indices": "1",
                }],
                [
                    "prefix_n",
                    "threshold_status",
                    "epsilon_certified_n",
                    "epsilon_star_estimate",
                    "triggering_sorted_indices",
                ],
            )
            args = make_args(output)
            corrected = [{
                "prefix_n": 1,
                "threshold_status": "resolved",
                "epsilon_certified_n": 0.049,
                "epsilon_star_estimate": 0.05,
                "triggering_sorted_indices": "1",
            }]
            rows = audit.legacy_threshold_comparison_rows(args, corrected, [])
            target = output / "baseline_legacy_threshold_comparison.csv"
            audit.write_csv(target, rows, audit.LEGACY_COMPARISON_FIELDS)
            self.assertTrue(target.exists())
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["scientific_status"], "unchanged_within_tolerance")


if __name__ == "__main__":
    unittest.main()
