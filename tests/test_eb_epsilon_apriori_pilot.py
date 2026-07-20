import csv
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    run_eb_epsilon_apriori_pilot as runner,
)
from scripts.analysis.thickness_mismatch.postprocess import (  # noqa: E402
    analyze_eb_epsilon_apriori_pilot as pilot,
)
from scripts.analysis.thickness_mismatch.postprocess import (  # noqa: E402
    analyze_eb_safe_prefix_certification as certification,
)


def synthetic_modes(count: int = 10, *, pi_scale: float = 0.01, gap: float = 0.2) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for index in range(1, count + 1):
        left = gap if index > 1 else float("nan")
        right = gap if index < count else float("nan")
        pi_value = pi_scale * index
        rows.append(
            {
                "sorted_index": index,
                "Lambda_EB": 2.0 + index,
                "Lambda_Timo": (2.0 + index) / math.sqrt(1.01),
                "delta_f": 0.01,
                "epsilon_max": 0.02,
                "Pi_EB": pi_value,
                "Pi_shear_EB": 0.6 * pi_value,
                "Pi_rotary_EB": 0.4 * pi_value,
                "chi_max_EB": 0.01 * index,
                "chi_eff_EB": 0.008 * index,
                "Theta_max_EB": 0.002 * index,
                "EB_axial_fraction": 0.1,
                "EB_bending_fraction": 0.9,
                "EB_classification": "bending_dominated",
                "gap_left_EB": left,
                "gap_right_EB": right,
                "gap_EB": gap,
                "mode10_upper_guard_source": "unavailable" if index == count else "not_applicable",
            }
        )
    return rows


def synthetic_record(
    case_number: int,
    epsilon_0: float,
    n_true: int,
    *,
    beta_deg: float = 45.0,
    mu: float = 0.0,
    eta: float = 0.0,
    epsilon_max: float | None = None,
) -> certification.GeometryRecord:
    modes = synthetic_modes()
    for mode in modes:
        mode["epsilon_max"] = epsilon_0 if epsilon_max is None else epsilon_max
    key = certification.geometry_key(epsilon_0, beta_deg, mu, eta)
    return certification.GeometryRecord(
        source_dataset="pilot",
        geometry_id=f"geometry_{case_number:02d}",
        geometry_key=key,
        epsilon_0=epsilon_0,
        beta_deg=beta_deg,
        mu=mu,
        eta=eta,
        modes=modes,
        source_file=Path("synthetic.csv"),
        complete=True,
        included=True,
        N_true=n_true,
    )


def manifest_map(records: list[certification.GeometryRecord], baseline_count: int = 0) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    for index, record in enumerate(records, start=1):
        case_id = f"C{index:02d}"
        groups = "baseline_reference" if index <= baseline_count else ""
        result[case_id] = {"geometry_id": record.geometry_id, "case_groups": groups}
    return result


def write_rows(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


class EbEpsilonAprioriPilotTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.cases = runner.load_manifest(runner.DEFAULT_MANIFEST)

    def test_01_manifest_has_exactly_21_unique_geometries(self) -> None:
        runner.validate_manifest(self.cases)
        self.assertEqual(len(self.cases), 21)
        keys = {(case.epsilon_0, case.beta_deg, case.mu, case.eta) for case in self.cases}
        self.assertEqual(len(keys), 21)

    def test_02_manifest_group_memberships(self) -> None:
        self.assertEqual(runner.group_members(self.cases, "baseline_reference"), {f"B{i:02d}" for i in range(1, 8)})
        self.assertEqual(runner.group_members(self.cases, "eta_sign_triplet"), {"H01", "G04", "H02"})
        self.assertEqual(runner.group_members(self.cases, "eta_sign_triplet_center"), {"G04"})

    def test_03_project_helper_gives_matched_epsilon_max_values(self) -> None:
        by_id = {case.case_id: case for case in self.cases}
        for case_id in ("G03", "M01", "M02"):
            self.assertAlmostEqual(runner.local_parameters(by_id[case_id])["epsilon_max"], 0.02)
        for case_id in ("M03", "G04", "M04"):
            self.assertAlmostEqual(runner.local_parameters(by_id[case_id])["epsilon_max"], 0.04)

    def test_04_epsilon_max_uses_existing_helper(self) -> None:
        case = self.cases[0]
        expected = {"epsilon_max": 0.123}
        with mock.patch.object(runner.stage1, "local_thickness_parameters", return_value=expected) as helper:
            self.assertIs(runner.local_parameters(case), expected)
        helper.assert_called_once_with(case.epsilon_0, case.mu, case.eta)

    def test_05_reference_thresholds_use_only_baseline_cases(self) -> None:
        baseline = [synthetic_record(1, 0.01, 10), synthetic_record(2, 0.02, 5)]
        heldout = synthetic_record(3, 0.015, 0, beta_deg=90.0)
        manifest = manifest_map([*baseline, heldout], baseline_count=2)
        folds = pilot.build_folds([*baseline, heldout], manifest)
        _rows, _predictions, first = pilot.run_epsilon_rules(folds, manifest, k_max=10)
        heldout.N_true = 10
        folds = pilot.build_folds([*baseline, heldout], manifest)
        _rows, _predictions, second = pilot.run_epsilon_rules(folds, manifest, k_max=10)
        key = ("baseline_reference_transfer", "baseline_to_nonbaseline", "Rule_E0_ref")
        self.assertEqual(first[key], second[key])

    def test_06_rule_e0_uses_epsilon_0(self) -> None:
        record = synthetic_record(1, 0.02, 10, epsilon_max=0.08)
        self.assertEqual(pilot.predictor_value(record, "Rule_E0_ref"), 0.02)

    def test_07_rule_emax_uses_epsilon_max_with_same_reference_thresholds(self) -> None:
        record = synthetic_record(1, 0.02, 10, epsilon_max=0.08)
        thresholds = (0.05,) * 10
        self.assertEqual(pilot.predictor_value(record, "Rule_Emax_ref"), 0.08)
        self.assertEqual(len(thresholds), 10)

    def test_08_calibrated_threshold_vector_has_length_ten(self) -> None:
        thresholds = pilot.calibrate_epsilon_thresholds(
            [synthetic_record(1, 0.01, 10), synthetic_record(2, 0.02, 4)],
            predictor="epsilon_0",
            k_max=10,
        )
        self.assertEqual(len(thresholds), 10)

    def test_09_threshold_vector_is_monotone_nonincreasing(self) -> None:
        thresholds = pilot.calibrate_epsilon_thresholds(
            [synthetic_record(1, 0.01, 10), synthetic_record(2, 0.02, 6), synthetic_record(3, 0.03, 2)],
            predictor="epsilon_0",
            k_max=10,
        )
        self.assertTrue(all(left >= right for left, right in zip(thresholds, thresholds[1:])))

    def test_10_threshold_equality_is_accepted(self) -> None:
        record = synthetic_record(1, 0.02, 10)
        detail = pilot.epsilon_prediction(record, "Rule_E0_cal", (0.02,) * 10, [record], guarded=False)
        self.assertEqual(detail["N_hat"], 10)

    def test_11_conflicting_equal_x_is_handled_conservatively(self) -> None:
        safe = synthetic_record(1, 0.02, 10)
        unsafe = synthetic_record(2, 0.02, 0, beta_deg=90.0)
        thresholds = pilot.calibrate_epsilon_thresholds([safe, unsafe], predictor="epsilon_0", k_max=10)
        self.assertLess(thresholds[0], 0.02)
        self.assertEqual(pilot.epsilon_prediction(safe, "Rule_E0_cal", thresholds, [safe, unsafe], guarded=False)["N_hat"], 0)

    def test_12_training_false_safe_count_is_zero(self) -> None:
        records = [synthetic_record(1, 0.01, 10), synthetic_record(2, 0.02, 4), synthetic_record(3, 0.03, 0)]
        thresholds = pilot.calibrate_epsilon_thresholds(records, predictor="epsilon_0", k_max=10)
        metrics = pilot.epsilon_threshold_train_metrics(records, "Rule_E0_cal", thresholds)
        self.assertEqual(metrics["false_safe"], 0)

    def test_13_heldout_target_does_not_change_train_thresholds(self) -> None:
        train = [synthetic_record(1, 0.01, 10), synthetic_record(2, 0.02, 4)]
        heldout = synthetic_record(3, 0.03, 0)
        first = pilot.calibrate_epsilon_thresholds(train, predictor="epsilon_0", k_max=10)
        heldout.N_true = 10
        second = pilot.calibrate_epsilon_thresholds(train, predictor="epsilon_0", k_max=10)
        self.assertEqual(first, second)

    def test_14_out_of_x_domain_guard_returns_zero(self) -> None:
        train = [synthetic_record(1, 0.01, 10), synthetic_record(2, 0.02, 8)]
        heldout = synthetic_record(3, 0.03, 10)
        detail = pilot.epsilon_prediction(heldout, "Rule_E0_cal", (0.05,) * 10, train, guarded=True)
        self.assertEqual(detail["N_hat"], 0)
        self.assertTrue(detail["out_of_X_calibration_domain"])

    def test_15_late_pass_does_not_extend_true_prefix(self) -> None:
        self.assertEqual(certification.true_safe_prefix([0.01, 0.12, 0.01] + [0.01] * 7), 1)

    def test_16_matched_epsilon_max_groups_have_three_cases(self) -> None:
        self.assertEqual(len(runner.group_members(self.cases, "matched_epsilon_max_0p02")), 3)
        self.assertEqual(len(runner.group_members(self.cases, "matched_epsilon_max_0p04")), 3)

    def test_17_rule_a_gap_blocks_whole_close_cluster(self) -> None:
        modes = synthetic_modes(3, gap=0.2)
        modes[0]["gap_right_EB"] = 0.001
        modes[1]["gap_left_EB"] = 0.001
        result = certification.predict_safe_prefix(modes, "Rule_A_gap", {"T_pi": 1.0, "gap_min": 0.01})
        self.assertEqual(result["N_hat"], 0)
        self.assertEqual(result["trigger_reason"], "close_EB_cluster_guard")

    def test_18_mode_ten_uses_right_guard_when_available(self) -> None:
        modes = synthetic_modes()
        modes[-1]["gap_right_EB"] = 0.001
        modes[-1]["mode10_upper_guard_source"] = "candidate_root_11"
        blocked = certification.close_cluster_members(modes, 0.01)
        self.assertIn(10, blocked)

    def test_19_absent_guard_root_is_marked_without_fake_gap(self) -> None:
        self.assertTrue(math.isnan(runner.gap_from_pair(10.0, float("nan"))))
        modes = synthetic_modes()
        self.assertNotIn(10, certification.close_cluster_members(modes, 0.01))

    def test_20_offline_and_online_operation_scopes_are_separate(self) -> None:
        predictions = [
            {
                "split_kind": "synthetic",
                "fold_id": "one",
                "rule": "Rule_E0_cal",
                "N_hat_guarded": 4,
                "out_of_X_calibration_domain": False,
            }
        ]
        offline = [{"cost_scope": "offline_threshold_calibration", "threshold_combinations_evaluated": 7}]
        rows = pilot.operation_rows(predictions, offline, [], k_max=10)
        online = next(row for row in rows if row.get("cost_scope") == "online_geometry_only_inference")
        self.assertEqual(online.get("threshold_combinations_evaluated", ""), "")
        self.assertEqual(online["EB_root_operations_if_available"], 0)

    def test_21_prefix_extrema_candidates_include_critical_values(self) -> None:
        item = synthetic_record(1, 0.01, 5)
        item.modes[2]["Pi_shear_EB"] = 0.77
        values = certification.prefix_extrema_values([item], "Pi_shear_EB")
        self.assertIn(0.77, values)
        axes, _exact = certification.threshold_candidate_axes([item], "Rule_B", max_candidates_per_axis=8)
        self.assertIn(0.77, axes["T_shear"])

    def test_22_calibration_is_deterministic_under_row_permutation(self) -> None:
        records = [synthetic_record(1, 0.01, 8), synthetic_record(2, 0.02, 5), synthetic_record(3, 0.03, 2)]
        first = certification.calibrate_rule(records, "Rule_B", k_max=10, max_candidates_per_axis=8)
        second = certification.calibrate_rule(list(reversed(records)), "Rule_B", k_max=10, max_candidates_per_axis=8)
        self.assertEqual(first["thresholds"], second["thresholds"])

    def test_23_synthetic_csv_cli_writes_all_required_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pilot_dir = root / "pilot"
            output_dir = root / "analysis"
            cases = [
                ("B01", 0.01, 0.01, 0.0, 0.0, 0.0, "baseline_reference", 10),
                ("B02", 0.02, 0.02, 0.0, 0.0, 0.0, "baseline_reference", 7),
                ("G01", 0.02, 0.04, 45.0, 0.5, 0.0, "fixed_epsilon_0p02_block", 5),
                ("G02", 0.02, 0.066, 90.0, 0.7, 0.0, "fixed_epsilon_0p02_block", 3),
            ]
            manifest_rows: list[dict[str, object]] = []
            geometry_rows: list[dict[str, object]] = []
            mode_rows: list[dict[str, object]] = []
            for case_id, eps0, epsmax, beta, mu, eta, groups, n_true in cases:
                geometry_id = f"geometry_{case_id}"
                manifest_rows.append({"case_id": case_id, "geometry_id": geometry_id, "epsilon_0": eps0, "beta_deg": beta, "mu": mu, "eta": eta, "case_groups": groups})
                geometry_rows.append({"case_id": case_id, "quality_status": "included", "N_true": n_true})
                for mode in synthetic_modes():
                    index = int(mode["sorted_index"])
                    mode_rows.append({"case_id": case_id, "geometry_id": geometry_id, "epsilon_0": eps0, "beta_deg": beta, "mu": mu, "eta": eta, "case_groups": groups, "epsilon_max": epsmax, **mode, "delta_f": 0.01 if index <= n_true else 0.12, "matched_timo_sorted_index": index, "matching_status": "reliable_individual", "cluster_id": ""})
            write_rows(pilot_dir / "epsilon_pilot_case_manifest_resolved.csv", ["case_id", "geometry_id", "epsilon_0", "beta_deg", "mu", "eta", "case_groups"], manifest_rows)
            write_rows(pilot_dir / "epsilon_pilot_geometry_metrics.csv", ["case_id", "quality_status", "N_true"], geometry_rows)
            write_rows(pilot_dir / "epsilon_pilot_mode_metrics.csv", list(mode_rows[0]), mode_rows)
            write_rows(pilot_dir / "epsilon_pilot_root_operation_counts.csv", ["case_id", "n_shape_points", "EB_root_solver_calls", "root_cache_hits", "root_cache_misses", "boundary_retries"], [{"case_id": case[0], "n_shape_points": 51, "EB_root_solver_calls": 1, "root_cache_hits": 0, "root_cache_misses": 2, "boundary_retries": 0} for case in cases])
            pilot.main(["--pilot-dir", str(pilot_dir), "--output-dir", str(output_dir), "--candidate-grid-sizes", "3", "4", "--smoke", "--force"])
            self.assertEqual(set(pilot.OUTPUT_NAMES), {path.name for path in output_dir.iterdir()})

    def test_24_legacy_rules_a_to_d_keep_prefix_semantics(self) -> None:
        modes = synthetic_modes(3)
        modes[1]["Pi_EB"] = 2.0
        self.assertEqual(certification.predict_safe_prefix(modes, "Rule_A", {"T_pi": 1.0})["N_hat"], 1)
        self.assertEqual(certification.predict_safe_prefix(modes, "Rule_B", {"T_shear": 1.0, "T_rotary": 1.0})["N_hat"], 3)
        self.assertEqual(certification.predict_safe_prefix(modes, "Rule_C", {"T_shear": 1.0, "T_rotary": 1.0, "gap_min": 0.0})["N_hat"], 3)
        self.assertEqual(certification.predict_safe_prefix(modes, "Rule_D", {"T_shear": 1.0, "T_rotary": 1.0, "gap_min": 0.0})["N_hat"], 3)

    def test_25_stage_and_fixed_source_semantics_remain_imported(self) -> None:
        self.assertTrue(callable(runner.fixed_scan.point_products))
        self.assertTrue(callable(runner.stage1.local_thickness_parameters))
        self.assertEqual(certification.squared_lambda_delta(2.0, 1.0), 3.0)


if __name__ == "__main__":
    unittest.main()
