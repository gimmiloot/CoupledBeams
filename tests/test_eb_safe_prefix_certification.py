import csv
import math
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_eb_validity_fixed_epsilon_geometry_scan as fixed_scan,
)
from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_eb_validity_vs_timoshenko_stage1 as stage1,
)
from scripts.analysis.thickness_mismatch.postprocess import (  # noqa: E402
    analyze_eb_safe_prefix_certification as certification,
)


def mode_rows(
    pi_values: list[float],
    *,
    shear_values: list[float] | None = None,
    rotary_values: list[float] | None = None,
    classifications: list[str] | None = None,
    pair_gaps: list[float] | None = None,
) -> list[dict[str, object]]:
    count = len(pi_values)
    shear = shear_values or [0.5 * value for value in pi_values]
    rotary = rotary_values or [0.5 * value for value in pi_values]
    characters = classifications or ["bending_dominated"] * count
    gaps = pair_gaps or [0.2] * max(count - 1, 0)
    rows: list[dict[str, object]] = []
    for position in range(count):
        left = gaps[position - 1] if position > 0 else float("nan")
        right = gaps[position] if position < len(gaps) else float("nan")
        rows.append(
            {
                "sorted_index": position + 1,
                "Pi_EB": pi_values[position],
                "Pi_shear_EB": shear[position],
                "Pi_rotary_EB": rotary[position],
                "EB_classification": characters[position],
                "gap_left_EB": left,
                "gap_right_EB": right,
                "gap_EB": min(value for value in (left, right) if math.isfinite(value)) if count > 1 else float("nan"),
            }
        )
    return rows


def record(
    key: tuple[float, float, float, float],
    modes: list[dict[str, object]],
    n_true: int,
    *,
    source: str = certification.SOURCE_STAGE1,
) -> certification.GeometryRecord:
    return certification.GeometryRecord(
        source_dataset=source,
        geometry_id=certification.geometry_identifier(key),
        geometry_key=key,
        epsilon_0=key[0],
        beta_deg=key[1],
        mu=key[2],
        eta=key[3],
        modes=modes,
        source_file=Path("synthetic.csv"),
        complete=True,
        included=True,
        N_true=n_true,
    )


def synthetic_delta_pattern(first_failure: int | None) -> list[float]:
    values = [0.02 + 0.002 * index for index in range(10)]
    if first_failure is not None:
        values[first_failure - 1] = 0.12
        for index in range(first_failure, 10):
            values[index] = 0.03 + 0.001 * index
    return values


def write_source_csv(
    path: Path,
    geometries: list[tuple[tuple[float, float, float, float], list[float]]],
    *,
    fixed_source: bool,
) -> None:
    fields = [
        "comparison_type",
        "epsilon_0",
        "beta_deg",
        "mu",
        "eta",
        "eb_sorted_index",
        "timo_sorted_index",
        "Lambda_EB",
        "Lambda_Timo",
        "delta_f",
        "root_warning",
        "candidate_boundary",
        "tracking_warning",
        "notes",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for key, deltas in geometries:
            epsilon, beta, mu, eta = key
            for index, delta in enumerate(deltas, start=1):
                lambda_eb = 2.0 + 0.6 * index
                lambda_timo = lambda_eb / math.sqrt(1.0 + delta)
                writer.writerow(
                    {
                        "comparison_type": "sorted_index",
                        "epsilon_0": epsilon,
                        "beta_deg": beta,
                        "mu": mu,
                        "eta": eta,
                        "eb_sorted_index": index,
                        "timo_sorted_index": index,
                        "Lambda_EB": f"{lambda_eb:.16e}",
                        "Lambda_Timo": f"{lambda_timo:.16e}",
                        "delta_f": f"{certification.squared_lambda_delta(lambda_eb, lambda_timo):.16e}",
                        "root_warning": "false" if fixed_source else "",
                        "candidate_boundary": "false" if fixed_source else "",
                        "tracking_warning": "",
                        "notes": "synthetic K=10 fixture",
                    }
                )


class EbSafePrefixCertificationTest(unittest.TestCase):
    def test_true_safe_prefix_stops_at_first_failure(self) -> None:
        values = [0.01, 0.08, 0.12, 0.04] + [0.01] * 6
        self.assertEqual(certification.true_safe_prefix(values), 2)

    def test_true_safe_prefix_all_ten_pass(self) -> None:
        self.assertEqual(certification.true_safe_prefix([0.01] * 10), 10)

    def test_true_safe_prefix_first_failure_returns_zero(self) -> None:
        self.assertEqual(certification.true_safe_prefix([0.11] + [0.01] * 9), 0)

    def test_threshold_equality_is_safe(self) -> None:
        self.assertEqual(certification.true_safe_prefix([0.10] * 10), 10)

    def test_frequency_discrepancy_uses_squared_lambdas(self) -> None:
        self.assertAlmostEqual(certification.squared_lambda_delta(2.0, 1.0), 3.0)

    def test_missing_mode_row_is_rejected_in_strict_mode(self) -> None:
        rows = [{"sorted_index": index, "delta_f": 0.01} for index in range(1, 10)]
        with self.assertRaisesRegex(ValueError, "missing sorted mode indices"):
            certification.deltas_from_mode_rows(rows, strict=True)

    def test_duplicate_sorted_index_is_detected(self) -> None:
        rows = [{"sorted_index": index, "delta_f": 0.01} for index in range(1, 11)]
        rows.append({"sorted_index": 4, "delta_f": 0.01})
        with self.assertRaisesRegex(ValueError, "duplicate sorted mode indices"):
            certification.deltas_from_mode_rows(rows, strict=True)

    def test_candidate_boundary_is_an_explicit_source_warning(self) -> None:
        warning, text = certification.source_warning(
            {"root_warning": "false", "candidate_boundary": "true"}
        )
        self.assertTrue(warning)
        self.assertEqual(text, "candidate_boundary")

    def test_rule_a_returns_a_continuous_prefix(self) -> None:
        modes = mode_rows([0.1, 0.3, 0.1])
        result = certification.predict_safe_prefix(modes, "Rule_A", {"T_pi": 0.2})
        self.assertEqual(result["N_hat"], 1)

    def test_rule_b_requires_shear_and_rotary_limits(self) -> None:
        modes = mode_rows(
            [0.1, 0.1],
            shear_values=[0.04, 0.04],
            rotary_values=[0.04, 0.20],
        )
        result = certification.predict_safe_prefix(
            modes,
            "Rule_B",
            {"T_shear": 0.05, "T_rotary": 0.05},
        )
        self.assertEqual(result["N_hat"], 1)
        self.assertEqual(result["trigger_reason"], "Pi_rotary_limit")

    def test_cluster_guard_blocks_the_whole_close_pair(self) -> None:
        modes = mode_rows([0.01, 0.01, 0.01], pair_gaps=[0.2, 0.001])
        result = certification.predict_safe_prefix(
            modes,
            "Rule_C",
            {"T_shear": 1.0, "T_rotary": 1.0, "gap_min": 0.01},
        )
        self.assertEqual(result["N_hat"], 1)
        self.assertEqual(result["trigger_reason"], "close_EB_cluster_guard")

    def test_rule_d_uses_only_eb_classification(self) -> None:
        modes = mode_rows(
            [0.01, 0.01],
            classifications=["bending_dominated", "mixed"],
        )
        modes[0]["Timo_classification"] = "mixed"
        modes[1]["Timo_classification"] = "bending_dominated"
        result = certification.predict_safe_prefix(
            modes,
            "Rule_D",
            {"T_shear": 1.0, "T_rotary": 1.0, "gap_min": 0.0},
        )
        self.assertEqual(result["N_hat"], 1)
        self.assertEqual(result["trigger_reason"], "mixed_EB_mode_guard")

    def test_rule_a_calibration_has_zero_observed_false_safe(self) -> None:
        records = [
            record((0.01, 45.0, 0.0, 0.0), mode_rows([0.05, 0.10, 0.40]), 2),
            record((0.02, 45.0, 0.2, 0.0), mode_rows([0.04, 0.30, 0.50]), 1),
        ]
        result = certification.calibrate_rule(records, "Rule_A", k_max=3, max_candidates_per_axis=4)
        metrics = result["metrics"]
        self.assertEqual(metrics["false_safe_geometry_count"], 0)
        self.assertEqual(result["objective_retained_frequencies"], 3)

    def test_calibrated_threshold_does_not_use_held_out_targets(self) -> None:
        train = [record((0.01, 45.0, 0.0, 0.0), mode_rows([0.05, 0.10, 0.40]), 2)]
        held_out = record((0.02, 45.0, 0.2, 0.0), mode_rows([0.02, 0.03, 0.04]), 3)
        first = certification.calibrate_rule(train, "Rule_A", k_max=3, max_candidates_per_axis=4)
        held_out.N_true = 0
        second = certification.calibrate_rule(train, "Rule_A", k_max=3, max_candidates_per_axis=4)
        self.assertEqual(first["thresholds"], second["thresholds"])

    def test_insufficient_transfer_fold_is_recorded_but_not_evaluated(self) -> None:
        records = [
            record((0.01, 45.0, 0.0, 0.0), mode_rows([0.01]), 1),
            record(
                (0.02, 30.0, 0.2, 0.1),
                mode_rows([0.20]),
                0,
                source=certification.SOURCE_FIXED,
            ),
        ]
        calibration, predictions, _operations = certification.run_calibration_and_validation(
            records,
            k_max=1,
            max_candidates_per_axis=3,
        )
        transfer = next(
            row
            for row in calibration
            if row["split_kind"] == "stage1_to_fixed_epsilon" and row["rule"] == "Rule_A"
        )
        self.assertEqual(transfer["fold_status"], "limited_training_size")
        self.assertEqual(transfer["thresholds_json"], "{}")
        self.assertFalse(
            any(row["split_kind"] == "stage1_to_fixed_epsilon" for row in predictions)
        )

    def test_cross_source_overlap_is_removed_from_validation(self) -> None:
        shared = (0.02, 45.0, 0.2, 0.0)
        train = [record(shared, mode_rows([0.01]), 1)]
        validation = [
            record(shared, mode_rows([0.01]), 1, source=certification.SOURCE_FIXED),
            record((0.02, 0.0, 0.0, 0.0), mode_rows([0.01]), 1, source=certification.SOURCE_FIXED),
        ]
        kept, overlap_count = certification.exclude_transfer_overlaps(train, validation)
        self.assertEqual(overlap_count, 1)
        self.assertEqual(len(kept), 1)

    def test_false_safe_and_conservative_metrics_are_separate(self) -> None:
        predictions = [
            {"N_true": 2, "N_hat": 4},
            {"N_true": 5, "N_hat": 3},
        ]
        metrics = certification.prediction_metrics(predictions, k_max=10)
        self.assertEqual(metrics["false_safe_frequency_count"], 2)
        self.assertEqual(metrics["mean_conservative_loss"], 1.0)
        self.assertEqual(metrics["retained_EB_frequency_count"], 5)

    def test_all_held_out_fold_count_keeps_split_identity(self) -> None:
        predictions = [
            {
                "split_kind": "stage1_to_fixed_epsilon",
                "fold_id": "transfer",
                "fold_status": "ok",
                "rule": "Rule_A",
                "N_true": 1,
                "N_hat": 1,
            },
            {
                "split_kind": "fixed_epsilon_to_stage1",
                "fold_id": "transfer",
                "fold_status": "ok",
                "rule": "Rule_A",
                "N_true": 1,
                "N_hat": 1,
            },
        ]
        rows = certification.validation_summary_rows(predictions, k_max=1)
        aggregate = next(row for row in rows if row["split_kind"] == "all_held_out")
        self.assertEqual(aggregate["fold_count"], 2)

    def test_stage1_reconstruction_uses_lambda_eb(self) -> None:
        key = (0.01, 45.0, 0.0, 0.0)
        lambda_eb = 3.9
        mode = {
            "sorted_index": 1,
            "Lambda_EB": lambda_eb,
            "Lambda_Timo": 99.0,
            "delta_f": certification.squared_lambda_delta(lambda_eb, 99.0),
        }
        item = record(key, [mode], 0)
        operations = certification.OperationCounts()
        certification.reconstruct_predictors(
            [item],
            n_shape_points=51,
            k_max=1,
            threshold=0.10,
            operations=operations,
        )
        expected = fixed_scan.chi_values_eb(lambda_eb, key[0], key[2], key[3])["chi_max_EB"]
        self.assertAlmostEqual(float(item.modes[0]["chi_max_EB"]), expected)

    def test_candidate_grid_is_deterministic(self) -> None:
        values = [0.4, 0.1, 0.3, 0.2, 0.5]
        first = certification.deterministic_candidate_axis(
            values,
            max_candidates=4,
            include_negative_infinity=True,
            priority_values=[0.3],
        )
        second = certification.deterministic_candidate_axis(
            list(reversed(values)),
            max_candidates=4,
            include_negative_infinity=True,
            priority_values=[0.3],
        )
        self.assertEqual(first, second)

    def test_operation_counts_are_reproducible(self) -> None:
        key = (0.01, 45.0, 0.0, 0.0)
        mode = {
            "sorted_index": 1,
            "Lambda_EB": 3.9,
            "Lambda_Timo": 3.8,
            "delta_f": certification.squared_lambda_delta(3.9, 3.8),
        }
        item = record(key, [mode], 0)
        operations = certification.OperationCounts()
        certification.reconstruct_predictors(
            [item],
            n_shape_points=51,
            k_max=1,
            threshold=0.10,
            operations=operations,
        )
        self.assertEqual(operations.EB_mode_reconstructions, 1)
        self.assertEqual(operations.EB_characteristic_matrix_evaluations, 1)
        self.assertEqual(operations.SVD_6x6_calls, 1)
        self.assertEqual(operations.EB_shape_grid_points, 102)
        self.assertEqual(operations.quadrature_calls, 12)
        self.assertEqual(operations.quadrature_point_evaluations, 612)

    def test_legacy_first8_and_firstK_fields_remain_distinct(self) -> None:
        deltas = [0.01] * 8 + [0.20, 0.01]
        stage_summary = stage1.k_aware_sorted_summary_fields(deltas, target_mode_count=10)
        fixed_summary = fixed_scan.k_aware_point_summary_fields(deltas, target_mode_count=10)
        self.assertEqual(stage_summary["n_pass_10_among_first8"], 8)
        self.assertEqual(stage_summary["n_pass_10_among_firstK"], 9)
        self.assertAlmostEqual(stage_summary["max_delta_f_firstK"], 0.20)
        self.assertEqual(fixed_summary["n_pass_10_among_first8"], 8)
        self.assertEqual(fixed_summary["n_pass_10_among_firstK"], 9)
        self.assertAlmostEqual(fixed_summary["max_delta_f_first8"], 0.01)
        self.assertAlmostEqual(fixed_summary["max_delta_f_firstK"], 0.20)

    def test_stage1_report_uses_dynamic_reported_mode_count(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            args = stage1.Args(
                beta_deg=45.0,
                eta=0.0,
                mu_values=(0.0,),
                epsilon_values=(0.0025,),
                n_reported_modes=10,
                n_candidate_roots=16,
                n_shape_points=51,
                output_dir=output,
                cache_dir=output / "cache",
                reuse_cache=True,
                force_recompute=False,
                plot_only=False,
                skip_critical_refinement=True,
                smoke=True,
                max_refinement_iterations=0,
            )
            report = stage1.write_report(
                output_dir=output,
                args=args,
                mode_rows=[],
                summary_rows=[],
                critical_rows=[],
                tracking_rows=[],
                threshold_rows=[],
                chi_fit_rows=[],
                warnings=[],
                plot_paths=[],
            )
            text = report.read_text(encoding="utf-8")
            self.assertIn("base_branch_1..base_branch_10", text)
            self.assertNotIn("base_branch_1..base_branch_8", text)

    def test_stage1_plot_only_recovers_saved_base_pairing_warning(self) -> None:
        rows = [
            {
                "comparison_type": "physical_branch",
                "mu": 0.0,
                "epsilon_0": 0.0025,
                "branch_id": "base_branch_9",
                "eb_sorted_index": 9,
                "timo_sorted_index": 10,
                "tracking_warning": "",
            }
        ]
        warnings = stage1.warnings_from_saved_audits(rows, [])
        self.assertEqual(len(warnings), 1)
        self.assertIn("EB sorted 9 matched Timoshenko sorted 10", warnings[0])

    def test_synthetic_cli_writes_all_required_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stage_dir = root / "stage1"
            fixed_dir = root / "fixed"
            output_dir = root / "output"
            overlap_key = (0.02, 45.0, 0.2, 0.0)
            stage_geometries = [
                ((0.01, 45.0, 0.0, 0.0), synthetic_delta_pattern(None)),
                (overlap_key, synthetic_delta_pattern(7)),
                ((0.03, 45.0, 0.4, 0.0), synthetic_delta_pattern(5)),
            ]
            fixed_geometries = [
                (overlap_key, synthetic_delta_pattern(7)),
                ((0.02, 0.0, 0.0, -0.1), synthetic_delta_pattern(None)),
                ((0.02, 90.0, 0.7, 0.1), synthetic_delta_pattern(4)),
            ]
            write_source_csv(
                stage_dir / "eb_timo_mode_level_metrics.csv",
                stage_geometries,
                fixed_source=False,
            )
            write_source_csv(
                fixed_dir / "fixed_epsilon_mode_metrics.csv",
                fixed_geometries,
                fixed_source=True,
            )
            result = certification.main(
                [
                    "--stage1-dir",
                    str(stage_dir),
                    "--fixed-epsilon-dir",
                    str(fixed_dir),
                    "--output-dir",
                    str(output_dir),
                    "--k-max",
                    "10",
                    "--n-shape-points",
                    "51",
                    "--max-candidates-per-axis",
                    "4",
                ]
            )
            required = {
                "eb_safe_prefix_mode_audit.csv",
                "eb_safe_prefix_geometry_metrics.csv",
                "eb_safe_prefix_rule_calibration.csv",
                "eb_safe_prefix_fold_predictions.csv",
                "eb_safe_prefix_validation_summary.csv",
                "eb_safe_prefix_exclusion_audit.csv",
                "eb_safe_prefix_source_overlap_audit.csv",
                "eb_safe_prefix_predictor_consistency_audit.csv",
                "eb_safe_prefix_operation_counts.csv",
                "eb_safe_prefix_certification_report.md",
            }
            self.assertEqual(required, {path.name for path in output_dir.iterdir()})
            self.assertEqual(len(result["records"]), 6)
            self.assertEqual(len(result["operation_rows"]) > 3, True)
            overlap_rows = certification.read_csv(output_dir / "eb_safe_prefix_source_overlap_audit.csv")
            self.assertEqual(len(overlap_rows), 1)


if __name__ == "__main__":
    unittest.main()
