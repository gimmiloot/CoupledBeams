from __future__ import annotations

from dataclasses import replace
import hashlib
import inspect
import itertools
from pathlib import Path
import shutil

import pytest

from scripts.analysis.thickness_mismatch.postprocess import analyze_eb_rule_ab_exact_pareto as audit
from scripts.analysis.thickness_mismatch.benchmarks import benchmark_rule_s_cost_break_even as benchmark


def make_record(
    case_id: str,
    *,
    n_true: int,
    shear: list[float],
    rotary: list[float],
    geometry: tuple[float, float, float, float] = (0.01, 0.0, 0.0, 0.0),
    source: str = "synthetic",
    partition: str = audit.PARTITION_DEVELOPMENT,
) -> audit.GeometryData:
    if len(shear) != audit.K_MAX or len(rotary) != audit.K_MAX:
        raise ValueError("synthetic records require ten modes")
    modes = [
        {
            "sorted_index": index,
            "Lambda_EB": float(index + 1),
            "Lambda_Timo": float(index + 1),
            "delta_f": 0.01 if index <= n_true else 0.11,
            "Pi_shear_EB": s_value,
            "Pi_rotary_EB": r_value,
            "Pi_EB": s_value + r_value,
        }
        for index, (s_value, r_value) in enumerate(zip(shear, rotary, strict=True), start=1)
    ]
    return audit.GeometryData(
        source=source,
        case_id=case_id,
        epsilon_0=geometry[0],
        beta_deg=geometry[1],
        mu=geometry[2],
        eta=geometry[3],
        modes=modes,
        N_true=n_true,
        quality_status="quality_approved",
        root_source="synthetic_saved_roots",
        predictor_source="synthetic",
        assigned_partition=partition,
    )


def test_continuous_safe_prefix_stops_at_first_failed_mode() -> None:
    record = make_record(
        "continuous",
        n_true=1,
        shear=[0.05, 0.30, 0.04] + [0.04] * 7,
        rotary=[0.05] * 10,
    )
    n_a, _, _ = audit.predict_rule_a(record, 0.20)
    n_b, _, _ = audit.predict_rule_b(record, 0.20, 0.20)
    n_s, _, _ = audit.predict_rule_s(record, 0.20)
    assert n_a == n_b == n_s == 1


def test_exact_rule_a_search_and_zero_false_safe() -> None:
    first = make_record(
        "A1",
        n_true=2,
        shear=[0.05, 0.10, 0.20] + [0.30] * 7,
        rotary=[0.05, 0.10, 0.10] + [0.10] * 7,
        geometry=(0.01, 0.0, 0.0, 0.0),
    )
    second = make_record(
        "A2",
        n_true=3,
        shear=[0.02, 0.08, 0.15, 0.25] + [0.30] * 6,
        rotary=[0.03, 0.07, 0.10, 0.10] + [0.10] * 6,
        geometry=(0.02, 0.0, 0.0, 0.0),
    )
    result = audit.exact_rule_a_search([first, second])
    assert result.selected_threshold == pytest.approx(0.25)
    assert result.objective == 5
    assert all(
        audit.predict_rule_a(record, result.selected_threshold)[0] <= record.N_true
        for record in (first, second)
    )
    assert result.candidates[0] == float("-inf")
    assert all(bool(row["search_is_exact"]) for row in result.rows)


def rule_b_records() -> list[audit.GeometryData]:
    return [
        make_record(
            "B1",
            n_true=1,
            shear=[0.10, 0.30] + [0.50] * 8,
            rotary=[0.10, 0.10] + [0.50] * 8,
            geometry=(0.01, 10.0, 0.1, 0.0),
        ),
        make_record(
            "B2",
            n_true=2,
            shear=[0.10, 0.20, 0.20] + [0.50] * 7,
            rotary=[0.05, 0.15, 0.30] + [0.50] * 7,
            geometry=(0.02, 20.0, 0.2, 0.0),
        ),
    ]


def rule_s_records() -> list[audit.GeometryData]:
    return [
        make_record(
            "S1",
            n_true=1,
            shear=[0.10, 0.30] + [0.50] * 8,
            rotary=[0.01] * 10,
            geometry=(0.01, 10.0, 0.1, 0.0),
        ),
        make_record(
            "S2",
            n_true=2,
            shear=[0.10, 0.20, 0.30] + [0.50] * 7,
            rotary=[0.01] * 10,
            geometry=(0.02, 20.0, 0.2, 0.0),
        ),
    ]


def independent_rule_b_oracle(records: list[audit.GeometryData]) -> tuple[int, set[tuple[float, float]]]:
    shear, rotary = audit.rule_b_candidate_axes(records)
    best = -1
    pairs: set[tuple[float, float]] = set()
    for pair in itertools.product(shear, rotary):
        predictions = [audit.predict_rule_b(record, *pair)[0] for record in records]
        if any(value > record.N_true for value, record in zip(predictions, records, strict=True)):
            continue
        objective = sum(predictions)
        if objective > best:
            best, pairs = objective, {pair}
        elif objective == best:
            pairs.add(pair)
    return best, pairs


def test_exact_rule_b_cartesian_pareto_tie_break_and_all_optima() -> None:
    records = rule_b_records()
    expected_objective, expected_pairs = independent_rule_b_oracle(records)
    brute = audit.brute_force_rule_b_search(records)
    pareto = audit.exact_rule_b_pareto_search(records)
    assert pareto.objective == brute.objective == expected_objective == 3
    assert set(pareto.all_optimal_pairs) == set(brute.all_optimal_pairs) == expected_pairs
    assert pareto.pareto_optimal_pairs == audit.pareto_nondominated_threshold_pairs(expected_pairs)
    assert pareto.representative == min(pareto.all_optimal_pairs, key=lambda pair: (pair[0], pair[1]))
    assert len(pareto.all_optimal_pairs) >= len(pareto.pareto_optimal_pairs) >= 1
    assert all(
        audit.predict_rule_b(record, *pareto.representative)[0] <= record.N_true
        for record in records
    )
    repeat = audit.exact_rule_b_pareto_search(records)
    assert repeat.all_optimal_pairs == pareto.all_optimal_pairs
    assert repeat.representative == pareto.representative


def test_exact_rule_s_search_is_independent_and_matches_rule_b_shear_threshold() -> None:
    records = rule_s_records()
    rule_s = audit.exact_rule_s_search(records)
    rule_b = audit.exact_rule_b_pareto_search(records)
    assert rule_s.objective == rule_b.objective == 3
    assert rule_s.selected_threshold == rule_b.representative[0]
    assert rule_s.candidates[0] == float("-inf")
    assert all(
        audit.predict_rule_s(record, rule_s.selected_threshold)[0] <= record.N_true
        for record in records
    )


def test_equal_optimum_sensitivity_checks_every_pair_and_hashes_match() -> None:
    development = rule_b_records()
    result = audit.exact_rule_b_pareto_search(development)
    optimal_rows, _pareto_rows = audit.threshold_pair_rows(result)
    common_shear = [0.10, 0.30] + [0.50] * 8
    common_rotary = [0.01] * 10
    validation = [
        make_record("VBASE", n_true=1, shear=common_shear, rotary=common_rotary, geometry=(0.03, 0, 0, 0), partition=audit.PARTITION_BASELINE),
        make_record("VDIR", n_true=1, shear=common_shear, rotary=common_rotary, geometry=(0.04, 45, 0.2, 0), partition=audit.PARTITION_DIRECTED),
        make_record("S3_12", n_true=1, shear=common_shear, rotary=common_rotary, geometry=(0.05, 90, 0.7, 0), partition=audit.PARTITION_HOLDOUT),
        make_record("S3_14", n_true=1, shear=common_shear, rotary=common_rotary, geometry=(0.06, 45, 0.5, -0.1), partition=audit.PARTITION_HOLDOUT),
    ]
    rows, equivalence, invariants, unstable = audit.sensitivity_rows(validation, result, optimal_rows)
    assert len(rows) == len(optimal_rows) * 5
    assert len(equivalence) == len(optimal_rows) * 4
    assert invariants["all_equal_optimum_pairs_checked"]
    assert invariants["all_equal_optimum_prediction_vectors_identical"]
    assert invariants["all_equal_optimum_safety_signatures_identical"]
    assert all(row["prediction_vector_equals_representative"] for row in rows)
    assert all(row["safety_signature_equals_representative"] for row in rows)
    assert not unstable


def test_rule_s_rule_b_equivalence_and_nonbinding_rotary() -> None:
    records = rule_s_records()
    predictions_s = audit.predictions_for_rule_s(records, 0.20)
    predictions_b = audit.predictions_for_rule_b(records, (0.20, 0.15), pair_id="representative")
    rows, summary = audit.rule_s_rule_b_equivalence_rows(
        records,
        predictions_s,
        predictions_b,
        threshold_s=0.20,
        threshold_r=0.15,
    )
    assert len(rows) == 2
    assert summary["Rule_S_equals_Rule_B_on_all_checked_geometries"]
    assert summary["Pi_rotary_nonbinding_on_all_checked_predictions"]
    assert all(not row["rotary_binding"] for row in rows)


def test_dominance_witness_and_crossing_coordinates() -> None:
    rows = [
        {
            "case_id": "safe_dominated",
            "prefix_n": 2,
            "prefix_label": "safe_prefix",
            "S_prefix_max_Pi_shear": 0.20,
            "R_prefix_max_Pi_rotary": 0.20,
        },
        {
            "case_id": "safe_crossing",
            "prefix_n": 1,
            "prefix_label": "safe_prefix",
            "S_prefix_max_Pi_shear": 0.05,
            "R_prefix_max_Pi_rotary": 0.25,
        },
        {
            "case_id": "unsafe",
            "prefix_n": 3,
            "prefix_label": "unsafe_prefix",
            "S_prefix_max_Pi_shear": 0.10,
            "R_prefix_max_Pi_rotary": 0.10,
        },
    ]
    witnesses, summary = audit.dominance_audit(rows)
    assert [(row["safe_case_id"], row["unsafe_case_id"]) for row in witnesses] == [
        ("safe_dominated", "unsafe")
    ]
    assert summary["safe_prefix_point_count"] == 2
    assert summary["unsafe_prefix_point_count"] == 1
    assert summary["unavoidably_rejected_safe_prefix_count"] == 1


def test_rule_s_threshold_margin_calculation() -> None:
    rows = [
        {"assigned_partition": audit.PARTITION_DEVELOPMENT, "case_id": "A", "prefix_n": 1, "prefix_label": "safe_prefix", "S_prefix_max_Pi_shear": 0.18},
        {"assigned_partition": audit.PARTITION_DEVELOPMENT, "case_id": "B", "prefix_n": 2, "prefix_label": "unsafe_prefix", "S_prefix_max_Pi_shear": 0.22},
        {"assigned_partition": audit.PARTITION_DIRECTED, "case_id": "V", "prefix_n": 1, "prefix_label": "safe_prefix", "S_prefix_max_Pi_shear": 0.19},
        {"assigned_partition": audit.PARTITION_DIRECTED, "case_id": "V", "prefix_n": 2, "prefix_label": "unsafe_prefix", "S_prefix_max_Pi_shear": 0.21},
        {"assigned_partition": audit.PARTITION_HOLDOUT, "case_id": "S3_12", "prefix_n": 1, "prefix_label": "safe_prefix", "S_prefix_max_Pi_shear": 0.17},
        {"assigned_partition": audit.PARTITION_HOLDOUT, "case_id": "S3_14", "prefix_n": 2, "prefix_label": "unsafe_prefix", "S_prefix_max_Pi_shear": 0.23},
    ]
    margins = {row["partition"]: row for row in audit.rule_s_threshold_margin_rows(rows, threshold_s=0.20)}
    assert margins[audit.PARTITION_DEVELOPMENT]["accepted_absolute_margin"] == pytest.approx(0.02)
    assert margins[audit.PARTITION_DEVELOPMENT]["unsafe_absolute_margin"] == pytest.approx(0.02)
    assert margins[audit.PARTITION_DIRECTED]["nearest_accepted_case_id"] == "V"
    assert margins["S3_combined"]["nearest_unsafe_above_case_id"] == "S3_14"


def test_rule_s_leave_one_geometry_out_recalibrates_without_changing_frozen_threshold() -> None:
    records = [
        replace(record, case_id=f"L{index}", epsilon_0=0.01 * index)
        for index, record in enumerate([*rule_s_records(), rule_s_records()[1]], start=1)
    ]
    rows, comparisons = audit.rule_s_leave_one_geometry_out_rows(records)
    assert len(rows) == 3
    assert comparisons > 0
    assert all(row["training_geometry_count"] == 2 for row in rows)
    assert all(not row["frozen_full_development_threshold_changed"] for row in rows)


def test_specialized_shear_predictor_matches_existing_pi_implementation() -> None:
    mode = audit.certification.fixed_scan.eb_mode_result(
        epsilon=0.02,
        beta_deg=45.0,
        mu=0.2,
        eta=0.0,
        sorted_index=1,
        Lambda=2.5,
        n_points=101,
        root_warnings=(),
    )
    existing = audit.certification.fixed_scan.eb_predictors(mode)["Pi_shear_EB"]
    specialized = audit.certification.fixed_scan.eb_shear_predictor(mode)["Pi_shear_EB"]
    assert specialized == pytest.approx(existing, rel=1.0e-13, abs=1.0e-15)


def test_retention_conservative_loss_and_worst_overprediction() -> None:
    metrics = audit.metrics_for_predictions(
        [
            {"N_true": 4, "N_hat": 2},
            {"N_true": 3, "N_hat": 5},
            {"N_true": 3, "N_hat": 3},
        ]
    )
    assert metrics["false_safe_geometry_count"] == 1
    assert metrics["false_safe_frequency_count"] == 2
    assert metrics["worst_overprediction"] == 2
    assert metrics["mean_conservative_loss"] == pytest.approx(2 / 3)
    assert metrics["maximum_conservative_loss"] == 2
    assert metrics["usable_frequency_retention"] == pytest.approx(0.8)


def test_partition_lock_deduplication_and_reporting_split() -> None:
    common = [0.01] * 10
    development = [
        make_record("D1", n_true=10, shear=common, rotary=common, geometry=(0.01, 0.0, 0.0, 0.0)),
        make_record("D2", n_true=10, shear=common, rotary=common, geometry=(0.02, 0.0, 0.0, 0.0)),
    ]
    step3a = [
        make_record("SBASE", n_true=10, shear=common, rotary=common, geometry=(0.03, 0.0, 0.0, 0.0), partition=""),
        make_record("SDIR", n_true=10, shear=common, rotary=common, geometry=(0.04, 45.0, 0.2, 0.0), partition=""),
        make_record("DUP", n_true=10, shear=common, rotary=common, geometry=(0.01 + 0.4e-12, 0.0, 0.0, 0.0), partition=""),
        make_record("S3_12", n_true=10, shear=common, rotary=common, geometry=(0.05, 90.0, 0.7, 0.0), partition=""),
        make_record("S3_14", n_true=10, shear=common, rotary=common, geometry=(0.06, 45.0, 0.5, -0.1), partition=""),
    ]
    audit.assign_partitions(development, step3a)
    partitions = {record.case_id: record.assigned_partition for record in [*development, *step3a]}
    assert partitions["SBASE"] == audit.PARTITION_BASELINE
    assert partitions["SDIR"] == audit.PARTITION_DIRECTED
    assert partitions["DUP"] == audit.PARTITION_EXCLUDED
    assert partitions["S3_12"] == partitions["S3_14"] == audit.PARTITION_HOLDOUT
    assert not {"S3_12", "S3_14"}.intersection(record.case_id for record in development)
    assert not {"S3_12", "S3_14"}.intersection(
        record.case_id for record in step3a if record.assigned_partition in {audit.PARTITION_BASELINE, audit.PARTITION_DIRECTED}
    )
    included = [record for record in [*development, *step3a] if record.assigned_partition != audit.PARTITION_EXCLUDED]
    assert not any(
        audit.same_geometry(left, right)
        for index, left in enumerate(included)
        for right in included[index + 1 :]
    )
    predictions_a = audit.predictions_for_rule_a(included, 1.0)
    predictions_b = audit.predictions_for_rule_b(included, (1.0, 1.0), pair_id="synthetic")
    summaries = audit.summary_partition_rows(predictions_a, predictions_b, threshold_a=1.0, threshold_b=(1.0, 1.0))
    counts = {(row["summary_partition"], row["rule"]): row["geometry_count"] for row in summaries}
    assert counts[(audit.PARTITION_BASELINE, "Rule_B")] == 1
    assert counts[(audit.PARTITION_DIRECTED, "Rule_B")] == 1


def test_online_cost_is_sequential_per_geometry() -> None:
    rows = audit.online_operation_rows(
        [
            {"assigned_partition": audit.PARTITION_DEVELOPMENT, "rule": "Rule_B", "case_id": "G1", "N_hat": 4},
            {"assigned_partition": audit.PARTITION_DEVELOPMENT, "rule": "Rule_B", "case_id": "G2", "N_hat": 10},
        ],
        n_shape_points=401,
    )
    geometry = {row["case_id"]: row for row in rows if row["scope"].endswith("_geometry")}
    assert geometry["G1"]["estimated_online_mode_reconstructions"] == 5
    assert geometry["G2"]["estimated_online_mode_reconstructions"] == 10
    assert geometry["G1"]["Timoshenko_suffix_frequencies"] == 6
    assert geometry["G1"]["estimated_online_scalar_comparisons"] == 10


def copy_csv_fixture(source: Path, target: Path, names: tuple[str, ...]) -> None:
    target.mkdir(parents=True)
    for name in names:
        shutil.copy2(source / name, target / name)


def test_csv_only_full_audit_never_calls_roots_and_is_reproducible(tmp_path: Path, monkeypatch) -> None:
    development_dir = tmp_path / "development"
    step3a_dir = tmp_path / "step3a"
    copy_csv_fixture(
        audit.DEFAULT_DEVELOPMENT_DIR,
        development_dir,
        (
            "epsilon_pilot_case_manifest_resolved.csv",
            "epsilon_pilot_geometry_metrics.csv",
            "epsilon_pilot_mode_metrics.csv",
        ),
    )
    copy_csv_fixture(
        audit.DEFAULT_STEP3A_DIR,
        step3a_dir,
        (
            "step3a_manifest_resolved.csv",
            "step3a_case_summary.csv",
            "step3a_case_spectrum_summary.csv",
            "step3a_mode_metrics.csv",
            "step3a_strict_verification_audit.csv",
            "step3a_baseline_control_audit.csv",
        ),
    )

    def forbidden_root_call(*_args, **_kwargs):
        raise AssertionError("root solver was called by CSV-only audit")

    monkeypatch.setattr(audit.certification.fixed_scan, "solve_point_roots", forbidden_root_call)
    monkeypatch.setattr(audit.certification.stage1, "solve_products", forbidden_root_call)
    source = inspect.getsource(audit)
    assert "solve_point_roots(" not in source
    assert "resolve_branch_spectrum(" not in source
    assert "solve_products(" not in source

    original_load_development = audit.load_development

    def prepared_development(args, counts):
        records = original_load_development(args, counts)
        for record in records:
            for mode in record.modes:
                index = int(mode["sorted_index"])
                mode["Pi_shear_EB"] = 0.01 * index
                mode["Pi_rotary_EB"] = 0.005 * index
                mode["Pi_EB"] = 0.015 * index
        return records

    monkeypatch.setattr(audit, "load_development", prepared_development)

    def prepared_predictors(records, *, n_shape_points, counts):
        del n_shape_points, counts
        equivalence = []
        for record in records:
            if record.assigned_partition == audit.PARTITION_EXCLUDED:
                continue
            for mode in record.modes:
                index = int(mode["sorted_index"])
                mode["Pi_shear_EB"] = 0.01 * index
                mode["Pi_rotary_EB"] = 0.005 * index
                mode["Pi_EB"] = 0.015 * index
                equivalence.append(
                    {
                        "case_id": record.case_id,
                        "sorted_index": index,
                        "Pi_shear_existing": mode["Pi_shear_EB"],
                        "Pi_shear_specialized": mode["Pi_shear_EB"],
                        "equivalent": True,
                    }
                )
        return equivalence

    monkeypatch.setattr(audit, "reconstruct_step3a_predictors", prepared_predictors)

    def prepared_plots(output_dir: Path):
        paths = tuple(output_dir / name for name in audit.PLOT_NAMES)
        for path in paths:
            path.write_bytes(b"prepared plot fixture\n")
        return paths

    monkeypatch.setattr(audit, "create_plots", prepared_plots)
    args_one = audit.Args(
        development_dir=development_dir,
        pilot_metadata_dir=tmp_path / "unused_metadata",
        step3a_dir=step3a_dir,
        output_dir=tmp_path / "output_one",
        n_shape_points=401,
        force=False,
        plot_only=False,
    )
    args_two = replace(args_one, output_dir=tmp_path / "output_two")
    first = audit.run_analysis(args_one)
    second = audit.run_analysis(args_two)
    assert first["root_calculations"] == second["root_calculations"] == 0
    for name in ("rule_A_exact_search.csv", "rule_B_exact_search_summary.csv", "rule_B_optimal_threshold_pairs.csv"):
        assert (args_one.output_dir / name).read_bytes() == (args_two.output_dir / name).read_bytes()

    numerical_hashes = {
        name: hashlib.sha256((args_one.output_dir / name).read_bytes()).hexdigest()
        for name in audit.NUMERICAL_OUTPUT_NAMES
    }
    audit.plot_only(replace(args_one, plot_only=True))
    assert numerical_hashes == {
        name: hashlib.sha256((args_one.output_dir / name).read_bytes()).hexdigest()
        for name in audit.NUMERICAL_OUTPUT_NAMES
    }
    with pytest.raises(FileExistsError):
        audit.ensure_output_policy(args_one)
    audit.ensure_output_policy(replace(args_one, force=True))


def test_import_smoke() -> None:
    assert audit.ALGORITHM_VERSION == "eb_rule_ab_exact_pareto_v2"
    assert callable(audit.main)


def test_data_quality_failure_is_explicitly_inconclusive(tmp_path: Path) -> None:
    args = audit.Args(
        development_dir=tmp_path / "missing_development",
        pilot_metadata_dir=tmp_path / "missing_metadata",
        step3a_dir=tmp_path / "missing_step3a",
        output_dir=tmp_path / "output",
        n_shape_points=401,
        force=False,
        plot_only=False,
    )
    report = audit.write_inconclusive_quality_report(args, audit.DataQualityError("missing approved saved roots"))
    text = report.read_text(encoding="utf-8")
    assert "rule_A_inconclusive_data_quality" in text
    assert "rule_B_inconclusive_data_quality" in text
    assert "New root calculations: `0`" in text


def benchmark_online_case(*, clustered: bool = False) -> benchmark.OnlineCase:
    roots = [float(index) for index in range(1, 11)]
    if clustered:
        roots[5] = roots[4] + 0.5 * benchmark.beta_workflow.SVD_RECOVERY_CLUSTER_TOL
        for index in range(6, 10):
            roots[index] = roots[index - 1] + 1.0
    return benchmark.OnlineCase(
        case_id="synthetic",
        epsilon_0=0.02,
        beta_deg=45.0,
        mu=0.5,
        eta=0.0,
        eb_roots=tuple(roots),
    )


def test_cost_proposal_is_exact_and_rule_s_threshold_is_frozen() -> None:
    cases, manifest, audits = benchmark.load_cases(benchmark.DEFAULT_PROPOSAL_CSV)
    assert tuple(case.case_id for case in cases) == benchmark.EXPECTED_CASE_IDS
    assert len({case.frozen_T_s for case in cases}) == 1
    assert cases[0].frozen_T_s == pytest.approx(0.16762413001084248)
    assert [case.expected_N_S for case in cases] == [10, 9, 4, 4, 4]
    assert all(row["online_timo_reference_policy"] == "withheld_until_post_solve_verification" for row in manifest)
    assert all(row["online_bracketing_uses_saved_Timo_references"] is False for row in audits)


def test_online_case_structurally_withholds_saved_timo_references() -> None:
    fields = set(benchmark.OnlineCase.__dataclass_fields__)
    assert "timo_references" not in fields
    assert "N_true" not in fields
    assert "expected_N_S" not in fields


def test_rule_s_selector_stops_at_first_rejected_mode(monkeypatch) -> None:
    values = iter((0.05, 0.10, 0.30, 0.01))
    calls: list[int] = []

    def fake_mode(**kwargs):
        calls.append(int(kwargs["sorted_index"]))
        return object()

    monkeypatch.setattr(benchmark.fixed_scan, "eb_mode_result", fake_mode)
    monkeypatch.setattr(
        benchmark.fixed_scan,
        "eb_shear_predictor",
        lambda _mode: {"Pi_shear_EB": next(values)},
    )
    counts = benchmark.OperationCounts()
    n_safe, predictors = benchmark.run_rule_s_selector(benchmark_online_case(), 0.20, 101, counts)
    assert n_safe == 2
    assert predictors == pytest.approx((0.05, 0.10, 0.30))
    assert calls == [1, 2, 3]
    assert counts.EB_mode_reconstructions == 3
    assert counts.EB_quadrature_calls == 18


def test_local_suffix_roots_remain_sorted_without_fallback(monkeypatch) -> None:
    monkeypatch.setattr(benchmark, "run_rule_s_selector", lambda *_args, **_kwargs: (4, (0.1,) * 5))
    monkeypatch.setattr(
        benchmark,
        "_local_root",
        lambda _case, mode_index, _counts: (float(mode_index) + 0.25, True, f"mode {mode_index}"),
    )

    def forbidden_fallback(*_args, **_kwargs):
        raise AssertionError("full spectrum fallback should not run")

    monkeypatch.setattr(benchmark, "solve_direct_timo", forbidden_fallback)
    result = benchmark.solve_hybrid_timo(benchmark_online_case(), 0.2, 101, benchmark.OperationCounts())
    assert result.N_S == 4
    assert result.suffix_roots == tuple(float(index) + 0.25 for index in range(5, 11))
    assert result.fallback_used is False


def test_cluster_risk_triggers_full_k10_fallback(monkeypatch) -> None:
    monkeypatch.setattr(benchmark, "run_rule_s_selector", lambda *_args, **_kwargs: (4, (0.1,) * 5))
    monkeypatch.setattr(
        benchmark,
        "_local_root",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("local roots must be skipped for cluster")),
    )

    def fake_direct(_case, counts):
        counts.Timo_root_solver_calls += 1
        counts.Timo_roots_requested += 10
        counts.Timo_roots_returned_by_solver_calls += 10
        return tuple(float(index) for index in range(1, 11))

    monkeypatch.setattr(benchmark, "solve_direct_timo", fake_direct)
    counts = benchmark.OperationCounts()
    result = benchmark.solve_hybrid_timo(benchmark_online_case(clustered=True), 0.2, 101, counts)
    assert result.fallback_used is True
    assert result.suffix_roots == tuple(float(index) for index in range(5, 11))
    assert counts.full_spectrum_fallback_count == 1
    assert counts.Timo_root_solver_calls == 1


def test_failed_local_attempt_is_counted_before_fallback(monkeypatch) -> None:
    monkeypatch.setattr(benchmark, "run_rule_s_selector", lambda *_args, **_kwargs: (8, (0.1,) * 9))

    def failed_local(_case, _mode_index, counts):
        counts.Timo_root_solver_calls += 1
        counts.Timo_roots_requested += 1
        counts.Timo_characteristic_determinant_evaluations += 80
        return float("nan"), False, "failed"

    def fake_direct(_case, counts):
        counts.Timo_root_solver_calls += 1
        counts.Timo_roots_requested += 10
        counts.Timo_roots_returned_by_solver_calls += 10
        counts.Timo_characteristic_determinant_evaluations += 100
        return tuple(float(index) for index in range(1, 11))

    monkeypatch.setattr(benchmark, "_local_root", failed_local)
    monkeypatch.setattr(benchmark, "solve_direct_timo", fake_direct)
    counts = benchmark.OperationCounts()
    result = benchmark.solve_hybrid_timo(benchmark_online_case(), 0.2, 101, counts)
    assert result.fallback_used
    assert counts.Timo_root_solver_calls == 2
    assert counts.Timo_roots_requested == 11
    assert counts.Timo_characteristic_determinant_evaluations == 180


def test_postsolve_reference_mismatch_triggers_global_fallback_without_reference_bracketing(monkeypatch) -> None:
    case = benchmark_online_case()
    initial = benchmark.HybridResult(
        N_S=8,
        pi_shear=(0.1,) * 9,
        suffix_roots=(90.0, 100.0),
        suffix_sources=("local_existing_EB_interval", "local_existing_EB_interval"),
        fallback_used=False,
        local_notes=("local solve completed",),
    )
    received_cases: list[benchmark.OnlineCase] = []

    def fake_direct(online, counts):
        received_cases.append(online)
        counts.Timo_root_solver_calls += 1
        counts.Timo_roots_requested += 10
        counts.Timo_roots_returned_by_solver_calls += 10
        return tuple(float(index) for index in range(1, 11))

    monkeypatch.setattr(benchmark, "solve_direct_timo", fake_direct)
    counts = benchmark.OperationCounts()
    result = benchmark.apply_postsolve_verification_fallback(
        case,
        tuple(float(index) for index in range(1, 11)),
        initial,
        counts,
    )
    assert received_cases == [case]
    assert "timo_references" not in received_cases[0].__dataclass_fields__
    assert result.suffix_roots == (9.0, 10.0)
    assert result.fallback_used
    assert counts.full_spectrum_fallback_count == 1


def test_direct_workflow_requests_exactly_k10(monkeypatch) -> None:
    requested: list[tuple[int, bool]] = []

    def fake_solve(_case, _beta, n_roots, *, retry=False, **_kwargs):
        requested.append((n_roots, retry))
        return benchmark.beta_workflow.RootResult(
            roots=tuple(float(index) for index in range(1, 11)),
            warnings=(),
            root_count_found=10,
            lambda_max_used=20.0,
            scan_step_used=0.01,
            retry_attempted=retry,
            retry_changed_value=False,
            notes=(),
        )

    monkeypatch.setattr(benchmark.beta_workflow, "solve_timo", fake_solve)
    counts = benchmark.OperationCounts()
    roots = benchmark.solve_direct_timo(benchmark_online_case(), counts)
    assert requested == [(10, False)]
    assert roots == tuple(float(index) for index in range(1, 11))
    assert counts.Timo_root_solver_calls == 1
    assert counts.Timo_roots_requested == 10


def test_benchmark_operation_aggregation_is_componentwise() -> None:
    rows = [
        benchmark.operation_row("A", "hybrid_Rule_S", benchmark.OperationCounts(
            EB_mode_reconstructions=2,
            Timo_root_solver_calls=1,
            Timo_roots_requested=3,
            Timo_characteristic_determinant_evaluations=11,
            suffix_frequency_count=3,
        )),
        benchmark.operation_row("B", "hybrid_Rule_S", benchmark.OperationCounts(
            EB_mode_reconstructions=4,
            Timo_root_solver_calls=2,
            Timo_roots_requested=12,
            Timo_characteristic_determinant_evaluations=29,
            full_spectrum_fallback_count=1,
            suffix_frequency_count=6,
        )),
    ]
    aggregate = benchmark.aggregate_operation_rows(rows)
    hybrid = next(row for row in aggregate if row["workflow"] == "hybrid_Rule_S")
    assert hybrid["EB_mode_reconstructions"] == 6
    assert hybrid["Timo_root_solver_calls"] == 3
    assert hybrid["Timo_roots_requested"] == 15
    assert hybrid["Timo_characteristic_determinant_evaluations"] == 40
    assert hybrid["full_spectrum_fallback_count"] == 1


def test_benchmark_plot_only_never_calls_roots_and_preserves_csv(tmp_path: Path, monkeypatch) -> None:
    output_dir = tmp_path / "benchmark"
    output_dir.mkdir()
    runtime_rows = []
    operation_rows = []
    for case_id in benchmark.EXPECTED_CASE_IDS:
        for workflow in ("direct_Timoshenko_K10", "hybrid_Rule_S"):
            runtime_rows.append({
                "algorithm_version": benchmark.ALGORITHM_VERSION,
                "case_id": case_id,
                "workflow": workflow,
                "repeat_index": 0,
                "temperature": "cold",
                "elapsed_seconds": 0.1,
                "result_hash": "fixture",
                "N_S": 4,
                "suffix_frequency_count": 6,
                "full_spectrum_fallback_count": 0,
                "Timo_root_solver_calls": 1,
                "Timo_roots_returned_by_solver_calls": 10,
                "Timo_characteristic_determinant_evaluations": 100,
                "EB_mode_reconstructions": 5,
            })
            operation_rows.append(benchmark.operation_row(case_id, workflow, benchmark.OperationCounts(
                Timo_characteristic_determinant_evaluations=100,
            )))
    benchmark.write_csv(
        output_dir / "runtime_measurements.csv",
        runtime_rows,
        benchmark.OUTPUT_FIELDS["runtime_measurements.csv"],
    )
    benchmark.write_csv(
        output_dir / "per_case_operation_counts.csv",
        operation_rows,
        benchmark.OUTPUT_FIELDS["per_case_operation_counts.csv"],
    )
    before = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in output_dir.glob("*.csv")
    }

    def forbidden(*_args, **_kwargs):
        raise AssertionError("root solver called during --plot-only")

    monkeypatch.setattr(benchmark, "solve_direct_timo", forbidden)
    monkeypatch.setattr(benchmark, "solve_hybrid_timo", forbidden)
    assert benchmark.main(["--output-dir", str(output_dir), "--plot-only"]) == 0
    after = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in output_dir.glob("*.csv")
    }
    assert before == after
    assert all((output_dir / name).is_file() for name in benchmark.PLOT_NAMES)


def test_benchmark_force_policy_and_import_smoke(tmp_path: Path) -> None:
    args = benchmark.Args(
        proposal_csv=benchmark.DEFAULT_PROPOSAL_CSV,
        output_dir=tmp_path,
        force=False,
        reuse_cache=False,
        plot_only=False,
        runtime_repeats=2,
        n_shape_points=101,
    )
    (tmp_path / "cost_break_even_summary.csv").write_text("fixture\n", encoding="utf-8")
    with pytest.raises(FileExistsError):
        benchmark.ensure_write_policy(args)
    benchmark.ensure_write_policy(replace(args, force=True))
    assert benchmark.ALGORITHM_VERSION == "rule_s_cost_break_even_v2"
    assert callable(benchmark.main)


def test_benchmark_reuse_cache_is_bound_to_all_numerical_inputs(tmp_path: Path) -> None:
    phase1 = tmp_path / "phase1"
    phase1.mkdir()
    for name in ("rule_B_cost_break_even_proposal.csv", "prefix_predictor_metrics.csv", "rule_S_predictions.csv"):
        shutil.copy2(benchmark.DEFAULT_PHASE1_DIR / name, phase1 / name)
    args = benchmark.Args(
        proposal_csv=phase1 / "rule_B_cost_break_even_proposal.csv",
        output_dir=tmp_path / "output",
        force=False,
        reuse_cache=True,
        plot_only=False,
        runtime_repeats=2,
        n_shape_points=101,
    )
    outputs = {"cost_break_even_summary.csv": [{"cache_reused": False}]}
    cache = tmp_path / "cache.json"
    benchmark.save_cache(cache, outputs, args)
    loaded = benchmark.load_cache(cache, args)
    assert loaded["cost_break_even_summary.csv"][0]["cache_reused"] is True
    prefix = phase1 / "prefix_predictor_metrics.csv"
    prefix.write_text(prefix.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="cache metadata"):
        benchmark.load_cache(cache, args)
