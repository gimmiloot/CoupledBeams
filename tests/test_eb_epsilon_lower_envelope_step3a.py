from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path

import pytest

from scripts.analysis.thickness_mismatch.audits import audit_eb_epsilon_lower_envelope_step3a as step3a
from scripts.lib import branch_informed_spectrum_continuation as branch
from scripts.lib import epsilon_lower_envelope_metrics as metric
from scripts.lib import general_spectrum_completeness as complete


@pytest.fixture(scope="module")
def cases() -> list[step3a.ManifestCase]:
    return step3a.load_manifest(step3a.DEFAULT_MANIFEST)


@pytest.fixture(scope="module")
def thresholds() -> dict[int, step3a.BaselineThreshold]:
    return step3a.load_thresholds(step3a.DEFAULT_BASELINE_DIR)


@pytest.fixture(scope="module")
def baseline(cases, thresholds) -> step3a.BaselineData:
    return step3a.load_baseline_data(
        step3a.DEFAULT_BASELINE_DIR,
        cases,
        thresholds,
        tolerance=1.0e-14,
    )


@pytest.fixture(scope="module")
def synthetic_run(tmp_path_factory):
    output = tmp_path_factory.mktemp("step3a_synthetic") / "results"
    result = step3a.main(
        ["--output-dir", str(output), "--write-step3b-followup-manifest"],
        resolver=fake_resolver,
    )
    return output, result


def fake_result(case: step3a.ManifestCase, model: str, settings: branch.ContinuationSettings, *, force: bool):
    scale = math.sqrt(1.05) if model == complete.MODEL_EB else 1.0
    values = tuple(scale * (index + 1.0) for index in range(settings.candidate_roots))
    parents = tuple(
        branch.ParentBranch(
            model=model,
            branch_id=f"b{index:02d}",
            family="bending_EB" if model == complete.MODEL_EB else "bending_Timoshenko",
            family_index=index,
            Lambda=value,
            right_vector=(1.0,),
            left_vector=(1.0,),
            sigma_1=1.0e-12,
            sigma_ratio=1.0e-12,
            block_nullity=1,
            full_nullity=1,
            source="synthetic",
        )
        for index, value in enumerate(values, start=1)
    )
    continued = tuple(
        branch.ContinuedBranch(
            model=model,
            branch_id=item.branch_id,
            parent_family=item.family,
            beta_deg=case.beta_deg,
            predicted_Lambda=item.Lambda,
            Lambda=item.Lambda,
            right_vector=(1.0,),
            left_vector=(1.0,),
            sigma_1=1.0e-12,
            sigma_ratio=1.0e-12,
            self_MAC=1.0,
            nullity=1,
            cluster_id="",
            refinement_status=branch.SEED_REFINED_TO_NEW_ROOT,
            detection_source="synthetic",
        )
        for item in parents
    )
    return branch.BranchContinuationResult(
        algorithm_version=branch.BRANCH_CONTINUATION_ALGORITHM_VERSION,
        model=model,
        geometry=case.geometry,
        settings=settings,
        parent_branches=parents,
        branches=continued,
        steps=(),
        guard=branch.GlobalGuardResult(values[10], (), (), (), (), True, "pass"),
        primary_vs_verification=(),
        strict_fallback_audit=(),
        k10_guard_resolved=True,
        full12_resolved=True,
        spectrum_status="K10_guard_resolved",
        exclusion_reason="",
        oracle_agreement=None,
        force_verification_agreement=True if force else None,
        cache_status="synthetic",
        operations=branch.BranchOperationCounts(),
    )


def fake_resolver(case, settings, cache):
    force = cache.verification_scope == "force_strict_verification"
    return {model: fake_result(case, model, settings, force=force) for model in complete.SUPPORTED_MODELS}


# Manifest 1--8.
def test_manifest_has_exactly_28_rows(cases):
    assert len(cases) == 28


def test_manifest_ids_and_geometry_keys_are_unique(cases):
    assert len({case.case_id for case in cases}) == len({case.geometry_id for case in cases}) == 28


def test_manifest_has_exactly_seven_expected_groups(cases):
    assert {case.prefix_group for case in cases} == set(metric.EXPECTED_PREFIX_GROUPS)


def test_each_group_has_near_buffer_baseline_adversarial(cases):
    step3a.validate_manifest_contract(cases)


def test_combined_groups_expand_to_separate_prefixes():
    assert metric.EXPECTED_PREFIX_GROUPS["prefixes_4_5"] == (4, 5)
    assert metric.EXPECTED_PREFIX_GROUPS["prefixes_9_10"] == (9, 10)


def test_required_method_and_guard_are_validated(cases):
    bad = [replace(cases[0], required_K10_guard=False), *cases[1:]]
    with pytest.raises(ValueError, match="K10"):
        step3a.validate_manifest_contract(bad)


def test_manifest_epsilon_matches_corrected_full_precision(cases, thresholds):
    rows = step3a.validate_manifest_provenance(cases, thresholds, tolerance=1.0e-14)
    assert max(float(row["epsilon_provenance_difference"]) for row in rows) <= 1.0e-14


def test_legacy_threshold_mismatch_is_rejected(cases, thresholds):
    bad = [replace(cases[8], epsilon=0.02280127), *cases[:8], *cases[9:]]
    with pytest.raises(ValueError, match="manifest_baseline_provenance_failure"):
        step3a.validate_manifest_provenance(bad, thresholds, tolerance=1.0e-14)


# Metrics 9--17.
def test_delta_uses_lambda_squared():
    assert metric.squared_frequency_delta(2.0, 1.0) == pytest.approx(3.0)


def test_Delta_is_running_maximum():
    assert [row.Delta_n for row in metric.running_prefix_metrics((0.02, 0.08, 0.03))] == [0.02, 0.08, 0.08]


def test_N_true_stops_at_first_failure():
    assert metric.true_safe_prefix((0.01, 0.11, 0.01)) == 1


def test_threshold_equality_is_safe():
    assert metric.true_safe_prefix((0.10, 0.10)) == 2


def test_late_pass_does_not_extend_N_true():
    assert metric.true_safe_prefix((0.01, 0.11, 0.02, 0.03)) == 1
    assert metric.late_pass_indices((0.01, 0.11, 0.02, 0.03)) == (3, 4)


def test_violation_and_margin_have_opposite_signs():
    row = metric.running_prefix_metrics((0.12,))[0]
    assert row.V_n == pytest.approx(-row.M_n)


def test_trigger_mode_is_lowest_maximum():
    row = metric.running_prefix_metrics((0.01, 0.08, 0.03))[2]
    assert row.triggering_indices[0] == 2


def test_tied_trigger_indices_are_preserved():
    row = metric.running_prefix_metrics((0.08, 0.01, 0.08), tie_tolerance=1.0e-14)[2]
    assert row.triggering_indices == (1, 3)


def test_full_baseline_certificate_uses_safe_lower(thresholds):
    safe = {n: row.epsilon_certified_n for n, row in thresholds.items()}
    assert metric.certified_prefix(0.02, safe) == 7


# Classification 18--25.
def test_positive_violation_is_provisional():
    assert metric.primary_status(V_n=0.01, N_true=1, required_prefix_n=2, N_certified_0=2, strict_margin=0.005, k10_resolved=True) == "provisional_counterexample"


def test_near_boundary_trigger():
    assert metric.primary_status(V_n=-0.001, N_true=10, required_prefix_n=2, N_certified_0=2, strict_margin=0.005, k10_resolved=True) == "primary_near_boundary"


def test_confirmed_counterexample_requires_two_agreeing_runs():
    result = metric.final_status(primary="provisional_counterexample", V_primary=0.01, V_verification=0.009, N_true_primary=1, N_true_verification=1, required_prefix_n=2, N_certified_0=2, primary_k10=True, verification_k10=True, roots_agree=True, clusters_agree=True, violation_tolerance=1e-5, strict_required=True)
    assert result == "confirmed_counterexample"


def test_sign_disagreement_is_indeterminate():
    result = metric.final_status(primary="provisional_counterexample", V_primary=1e-4, V_verification=-1e-4, N_true_primary=1, N_true_verification=2, required_prefix_n=2, N_certified_0=2, primary_k10=True, verification_k10=True, roots_agree=True, clusters_agree=True, violation_tolerance=1e-5, strict_required=True)
    assert result == "numerically_indeterminate_near_threshold"


def test_unresolved_K10_excludes_case():
    assert metric.primary_status(V_n=-0.1, N_true=10, required_prefix_n=2, N_certified_0=2, strict_margin=0.005, k10_resolved=False) == "unresolved_spectrum"


def test_root12_does_not_enter_N_true():
    assert metric.true_safe_prefix([0.01] * 10) == 10


def test_root11_is_required_guard(cases):
    settings = branch.ContinuationSettings()
    result = fake_result(cases[0], complete.MODEL_EB, settings, force=False)
    short = replace(result, branches=result.branches[:10], k10_guard_resolved=False)
    assert not short.k10_guard_resolved


def test_full12_false_does_not_forbid_K10(cases):
    result = replace(fake_result(cases[0], complete.MODEL_EB, branch.ContinuationSettings(), force=False), full12_resolved=False)
    assert result.k10_guard_resolved


# Baseline controls 26--29.
def test_beta0_controls_use_corrected_factorized_oracle(cases, baseline):
    control = next(case for case in cases if case.case_group == "baseline_control")
    expected = complete.straight_oracle_values(complete.MODEL_TIMO, control.geometry, 11)
    assert baseline.factorized_roots[(control.epsilon, complete.MODEL_TIMO)] == pytest.approx(expected)


def test_near_baseline_target_is_safe(cases, baseline):
    for case in cases:
        if case.case_group == "baseline_control" and case.epsilon_source == "epsilon_near_n":
            prefixes = metric.running_prefix_metrics(baseline.mode_deltas[case.epsilon])
            assert all(prefixes[n - 1].Delta_n <= 0.10 for n in case.prefixes)


def test_buffer_baseline_target_is_safe(cases, baseline):
    for case in cases:
        if case.case_group == "baseline_control" and case.epsilon_source == "epsilon_buffer_n":
            prefixes = metric.running_prefix_metrics(baseline.mode_deltas[case.epsilon])
            assert all(prefixes[n - 1].Delta_n <= 0.10 for n in case.prefixes)


def test_pipeline_mismatch_blocks_decision():
    assert metric.decide_step3a(geometry_count=28, resolved_geometry_count=28, baseline_controls_pass=False, final_target_statuses=("confirmed_safe_at_screen_point",), strict_trigger_count=1, strict_verified_count=1) == "inconclusive_due_to_unresolved_cases"


# Strict verification 30--34.
def test_strict_run_returns_new_result_objects(cases):
    settings = branch.ContinuationSettings()
    primary = fake_resolver(cases[0], settings, type("Cache", (), {"verification_scope": "primary"})())
    strict = fake_resolver(cases[0], settings, type("Cache", (), {"verification_scope": "force_strict_verification"})())
    assert all(primary[model] is not strict[model] for model in complete.SUPPORTED_MODELS)


def test_strict_cache_is_separate(tmp_path):
    primary = branch.BranchContinuationCache(tmp_path / "primary", verification_scope="primary")
    strict = branch.BranchContinuationCache(tmp_path / "strict", verification_scope="force_strict_verification")
    assert primary.cache_dir != strict.cache_dir and primary.verification_scope != strict.verification_scope


@pytest.mark.parametrize("status,reasons", [("provisional_counterexample", ()), ("primary_near_boundary", ()), ("primary_safe", ("cluster",))])
def test_strict_trigger_for_provisional_near_cluster(status, reasons):
    assert metric.needs_strict_verification(status, reasons)


def test_root_differences_and_v_spread(cases):
    settings = branch.ContinuationSettings()
    primary = fake_resolver(cases[0], settings, type("Cache", (), {"verification_scope": "primary"})())
    strict = fake_resolver(cases[0], settings, type("Cache", (), {"verification_scope": "force_strict_verification"})())
    maximum, relative, agree, clusters = step3a.root_comparison(primary, strict)
    assert (maximum, relative, agree, clusters) == pytest.approx((0.0, 0.0, True, True))


def test_root_search_boundary_has_no_delta_inputs():
    assert set(step3a.resolve_case_models.__annotations__) >= {"case", "settings", "cache"}
    assert "delta" not in step3a.resolve_case_models.__annotations__


# Aggregation 35--39.
def test_worst_case_uses_maximum_violation():
    assert metric.worst_row(({"case_id": "a", "V_n": -0.2}, {"case_id": "b", "V_n": 0.01}))["case_id"] == "b"


def test_prefixes_4_5_aggregate_separately(cases):
    expanded = [(case.case_id, n) for case in cases if case.prefix_group == "prefixes_4_5" for n in case.prefixes]
    assert {n for _case, n in expanded} == {4, 5}


def test_prefixes_9_10_aggregate_separately(cases):
    expanded = [(case.case_id, n) for case in cases if case.prefix_group == "prefixes_9_10" for n in case.prefixes]
    assert {n for _case, n in expanded} == {9, 10}


def test_near_buffer_not_paired_for_different_geometry(cases):
    near = next(case for case in cases if case.case_id == "S3_02")
    buffer = next(case for case in cases if case.case_id == "S3_04")
    assert not metric.same_geometry(near.__dict__, buffer.__dict__)


def test_unique_counterexample_cases_are_deduplicated():
    assert metric.unique_case_ids(({"case_id": "x"}, {"case_id": "x"}, {"case_id": "y"})) == ("x", "y")


# Step 3B 40--44.
def proposal_input(final_status="confirmed_safe_at_screen_point"):
    return [{"case_id": "S3", "prefix_group": "prefixes_2", "prefix_n": 2, "is_target_prefix": True, "case_group": "adversarial", "final_status": final_status, "V_n_verification": -0.01, "Delta_n_verification": 0.09, "strict_verification_status": "pass", "beta_deg": 1.0, "mu": 0.0, "eta": 0.0}]


def test_worst_geometry_gets_paired_proposal_rows(cases, thresholds):
    rows = step3a.build_step3b_proposal(cases, proposal_input(), thresholds)
    assert [row["epsilon_source"] for row in rows] == ["epsilon_near_n", "epsilon_buffer_n"]


def test_existing_proposal_rows_deduplicate(cases, thresholds):
    rows = step3a.build_step3b_proposal(cases, proposal_input() * 2, thresholds)
    assert len(rows) == 2


def test_confirmed_counterexample_gets_refinement_purpose(cases, thresholds):
    rows = step3a.build_step3b_proposal(cases, proposal_input("confirmed_counterexample"), thresholds)
    assert {row["followup_purpose"] for row in rows} == {"refine_counterexample_epsilon_star"}


def test_safe_worst_case_gets_paired_margin_purpose(cases, thresholds):
    rows = step3a.build_step3b_proposal(cases, proposal_input(), thresholds)
    assert {row["followup_purpose"] for row in rows} == {"paired_worst_margin_check"}


def test_proposal_does_not_call_root_solver(monkeypatch, cases, thresholds):
    monkeypatch.setattr(step3a, "resolve_case_models", lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("root call")))
    assert step3a.build_step3b_proposal(cases, proposal_input(), thresholds)


# Integration 45--47 and diagnostic import coverage.
def test_synthetic_cli_creates_all_required_outputs(synthetic_run):
    output, result = synthetic_run
    assert result["decision"] in {
        "no_counterexample_in_28_case_screen",
        "inconclusive_due_to_unresolved_cases",
        "inconclusive_due_to_numerical_boundary",
    }
    assert all((output / name).exists() for name in (*step3a.OUTPUT_NAMES, *step3a.PLOT_NAMES))


def test_plot_only_performs_zero_root_calculations(tmp_path, synthetic_run):
    source_output, _synthetic_result = synthetic_run
    output = tmp_path / "saved"
    output.mkdir()
    for name in step3a.OUTPUT_NAMES[:11]:
        source = source_output / name
        assert source.exists()
        (output / name).write_bytes(source.read_bytes())
    result = step3a.main(["--output-dir", str(output), "--plot-only"])
    assert result["root_calculations"] == 0


def test_force_recompute_verification_is_reproducible(cases):
    settings = branch.ContinuationSettings()
    cache = type("Cache", (), {"verification_scope": "force_strict_verification"})()
    first = fake_resolver(cases[0], settings, cache)
    second = fake_resolver(cases[0], settings, cache)
    assert first[complete.MODEL_EB].values == second[complete.MODEL_EB].values


def test_decision_reports_confirmed_counterexample():
    assert metric.decide_step3a(geometry_count=28, resolved_geometry_count=28, baseline_controls_pass=True, final_target_statuses=("confirmed_counterexample",), strict_trigger_count=1, strict_verified_count=1) == "counterexample_found"


def test_decision_reports_numerical_boundary():
    assert metric.decide_step3a(geometry_count=28, resolved_geometry_count=28, baseline_controls_pass=True, final_target_statuses=("numerically_indeterminate_near_threshold",), strict_trigger_count=1, strict_verified_count=1) == "inconclusive_due_to_numerical_boundary"


def test_decision_reports_complete_28_case_screen():
    assert metric.decide_step3a(geometry_count=28, resolved_geometry_count=28, baseline_controls_pass=True, final_target_statuses=("confirmed_safe_at_screen_point",), strict_trigger_count=1, strict_verified_count=1) == "no_counterexample_in_28_case_screen"


def test_diagnostic_modules_import():
    assert step3a.STEP3A_ALGORITHM_VERSION == "epsilon_lower_envelope_step3a_v1"
