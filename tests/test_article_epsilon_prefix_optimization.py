from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import gzip
import json

import pytest

from scripts.analysis.thickness_mismatch.benchmarks import (
    benchmark_article_epsilon_prefix_optimization as benchmark,
)
from scripts.analysis.thickness_mismatch.audits import (
    run_article_epsilon_upper_envelope_grid as runner,
)
from scripts.analysis.thickness_mismatch.postprocess import (
    analyze_article_epsilon_upper_envelope as analyzer,
)
from scripts.lib import article_epsilon_prefix_optimization as prefix
from scripts.lib import article_epsilon_upper_envelope as workflow
from scripts.lib import general_spectrum_completeness as complete


def point() -> workflow.GridPoint:
    return workflow.build_manifest(smoke=True)[0]


def resolved_assessments(_sessions, guard: int):
    assessment = prefix.PrefixGuardAssessment(
        guard, "prefix_guard_resolved", guard - 1, (), False, False
    )
    return {model: assessment for model in complete.SUPPORTED_MODELS}


class FakeSession:
    roots = {
        complete.MODEL_TIMO: tuple(float(index) for index in range(1, 13)),
        complete.MODEL_EB: tuple(float(index) for index in range(1, 13)),
    }
    calls: list[tuple[str, int, bool]] = []
    unresolved_model: str | None = None

    def __init__(self, model, geometry, *, state=None):
        self.model = model
        self.geometry = geometry
        self.highest_resolved_mode = int(state.get("highest_resolved_mode", 0)) if state else 0
        self.highest_guard_mode = int(state.get("highest_guard_mode", 0)) if state else 0
        self.full_reference_completed = bool(state.get("full_reference_completed", False)) if state else False
        self.latest_result = self._result(self.highest_guard_mode) if self.highest_guard_mode else None

    def _result(self, guard):
        values = self.roots[self.model]
        roots = [SimpleNamespace(Lambda=value) for value in values[: max(guard, 1)]]
        return SimpleNamespace(primary=SimpleNamespace(roots=roots))

    def extend_to_guard(
        self,
        guard_index,
        *,
        continuation_seeds=(),
        eb_seed_roots=(),
        full_reference=False,
    ):
        guard = int(guard_index)
        hit = self.highest_guard_mode >= guard and (not full_reference or self.full_reference_completed)
        if not hit:
            self.calls.append((self.model, guard, bool(full_reference)))
            self.highest_guard_mode = guard
            self.highest_resolved_mode = guard - 1
            self.full_reference_completed = self.full_reference_completed or bool(full_reference)
            self.latest_result = self._result(12 if full_reference else guard)
        status = "prefix_guard_unresolved" if self.unresolved_model == self.model else "prefix_guard_resolved"
        assessment = prefix.PrefixGuardAssessment(
            guard,
            status,
            guard - 1 if status.endswith("resolved") and "unresolved" not in status else 0,
            (() if status == "prefix_guard_resolved" else ("synthetic_guard_failure",)),
            False,
            False,
        )
        timing = prefix.StageTiming(self.model, guard, 0.01, 0.02, 0.001, 10, 20, 2, 4)
        return self.latest_result, assessment, timing, hit

    def snapshot(self):
        return {
            "model": self.model,
            "highest_resolved_mode": self.highest_resolved_mode,
            "highest_guard_mode": self.highest_guard_mode,
            "full_reference_completed": self.full_reference_completed,
            "resolved_roots": list(self.roots[self.model][: self.highest_guard_mode]),
            "latest_result": {},
            "primary_evaluator": {},
            "verification_evaluator": {},
            "stage_timings": [],
        }


@pytest.fixture(autouse=True)
def reset_fake() -> None:
    FakeSession.calls = []
    FakeSession.unresolved_model = None
    FakeSession.roots = {
        complete.MODEL_TIMO: tuple(float(index) for index in range(1, 13)),
        complete.MODEL_EB: tuple(float(index) for index in range(1, 13)),
    }


def run_fake(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, **kwargs):
    monkeypatch.setattr(prefix, "IncrementalModelSession", FakeSession)
    monkeypatch.setattr(prefix, "_stage_assessments", resolved_assessments)
    cache = prefix.PartialPointCache(tmp_path / "cache")
    return prefix.run_staged_point(point(), cache=cache, strict_policy="local", **kwargs)


def set_delta_failure(mode: int, delta: float = 0.2) -> None:
    eb = list(FakeSession.roots[complete.MODEL_EB])
    eb[mode - 1] = FakeSession.roots[complete.MODEL_TIMO][mode - 1] * (1.0 + delta) ** 0.5
    FakeSession.roots[complete.MODEL_EB] = tuple(eb)


def test_early_stop_returns_prefix_N_true(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    set_delta_failure(3)
    payload = run_fake(tmp_path, monkeypatch)
    assert payload["execution_status"] == "resolved_prefix_early_stop"
    assert payload["N_true"] == 2
    assert payload["first_failed_mode"] == 3
    assert payload["prefix_guard_resolved_through"] == 3


def test_late_safe_mode_does_not_restore_prefix(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    set_delta_failure(2)
    payload = run_fake(tmp_path, monkeypatch)
    assert payload["N_true"] == 1


def test_N_true_zero(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    set_delta_failure(1)
    payload = run_fake(tmp_path, monkeypatch)
    assert payload["N_true"] == 0


def test_N_true_10_requires_root_11(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    payload = run_fake(tmp_path, monkeypatch)
    assert payload["N_true"] == 10
    assert payload["full_K10_guard_status"] == "full_K10_guard_resolved"
    assert all((model, 11, False) in FakeSession.calls for model in complete.SUPPORTED_MODELS)


def test_failure_at_k_requests_right_guard(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    set_delta_failure(5)
    payload = run_fake(tmp_path, monkeypatch)
    assert payload["N_true"] == 4
    assert payload["prefix_guard_resolved_through"] == 5
    assert all(any(call[0] == model and call[1] >= 6 for call in FakeSession.calls) for model in complete.SUPPORTED_MODELS)


def test_unresolved_guard_has_no_N_true(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(prefix, "IncrementalModelSession", FakeSession)
    FakeSession.unresolved_model = complete.MODEL_TIMO

    def unresolved(sessions, guard):
        result = resolved_assessments(sessions, guard)
        result[complete.MODEL_TIMO] = prefix.PrefixGuardAssessment(
            guard, "prefix_guard_unresolved", 0, ("synthetic_guard_failure",), False, False
        )
        return result

    monkeypatch.setattr(prefix, "_stage_assessments", unresolved)
    payload = prefix.run_staged_point(
        point(), cache=prefix.PartialPointCache(tmp_path / "cache"), strict_policy="local"
    )
    assert payload["execution_status"] == "attempted_unresolved"
    assert payload["N_true"] is None


def test_basis_error_is_not_interpreted_as_physical_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    class BasisErrorSession(FakeSession):
        def extend_to_guard(self, *args, **kwargs):
            if self.model == complete.MODEL_TIMO:
                raise ValueError("unsupported Timoshenko basis regime")
            return super().extend_to_guard(*args, **kwargs)

    monkeypatch.setattr(prefix, "IncrementalModelSession", BasisErrorSession)
    monkeypatch.setattr(prefix, "_stage_assessments", resolved_assessments)
    payload = prefix.run_staged_point(
        point(), cache=prefix.PartialPointCache(tmp_path / "cache"), strict_policy="local"
    )
    assert payload["N_true"] is None
    assert str(payload["execution_status"]).startswith("attempted_")


def test_partial_cache_resumes_from_last_guard(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    set_delta_failure(2)
    cache = prefix.PartialPointCache(tmp_path / "cache")
    monkeypatch.setattr(prefix, "IncrementalModelSession", FakeSession)
    monkeypatch.setattr(prefix, "_stage_assessments", resolved_assessments)
    prefix.run_staged_point(point(), cache=cache, strict_policy="local")
    first_calls = list(FakeSession.calls)
    prefix.run_staged_point(point(), cache=cache, strict_policy="local", full_k10=True)
    later_calls = FakeSession.calls[len(first_calls) :]
    assert later_calls
    assert all(call[1] == 11 and call[2] for call in later_calls)
    assert not any(call[1] == 5 for call in later_calls)


def test_full_k10_after_stop_completes_same_cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    set_delta_failure(2)
    cache = prefix.PartialPointCache(tmp_path / "cache")
    monkeypatch.setattr(prefix, "IncrementalModelSession", FakeSession)
    monkeypatch.setattr(prefix, "_stage_assessments", resolved_assessments)
    prefix.run_staged_point(point(), cache=cache, strict_policy="local")
    full = prefix.run_staged_point(point(), cache=cache, strict_policy="local", full_k10=True)
    assert full["execution_status"] == "resolved_full_K10"
    assert full["N_true"] == 1


@pytest.mark.parametrize("stages", [(4, 7, 10), (2, 6, 10), (1, 5, 9, 10)])
def test_stage_block_size_does_not_change_N_true(
    stages, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    set_delta_failure(6)
    payload = run_fake(tmp_path / str(stages), monkeypatch, stage_safe_ends=stages)
    assert payload["N_true"] == 5


@pytest.mark.parametrize("strategy", prefix.PREFIX_STRATEGIES)
def test_paired_and_full_eb_progressive_timo_agree(
    strategy: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    set_delta_failure(4)
    payload = run_fake(tmp_path / strategy, monkeypatch, strategy=strategy)
    assert payload["N_true"] == 3


def test_local_strict_checks_full_prefix_count() -> None:
    assert prefix.local_strict_modes(5, 6) == (1, 2, 3, 4, 5, 6)


def test_near_threshold_triggers_strict() -> None:
    reasons = prefix.strict_trigger_reasons(
        first_failed_mode=5,
        deltas={4: 0.09, 5: 0.1001},
        assessments={},
    )
    assert "near_delta_threshold" in reasons


def test_far_threshold_clean_case_stays_local_in_auto() -> None:
    reasons = prefix.strict_trigger_reasons(
        first_failed_mode=2,
        deltas={1: 0.01, 2: 0.2},
        assessments={},
    )
    assert prefix.choose_effective_strict_policy("auto", reasons) == "local"


def test_declared_control_case_escalates_even_local_policy() -> None:
    assert prefix.choose_effective_strict_policy("local", ("deterministic_control_case",)) == "full"


def test_control_case_full_strict_failure_cannot_be_accepted(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    set_delta_failure(2)
    monkeypatch.setattr(prefix, "IncrementalModelSession", FakeSession)
    monkeypatch.setattr(prefix, "_stage_assessments", resolved_assessments)
    payload = prefix.run_staged_point(
        point(),
        cache=prefix.PartialPointCache(tmp_path / "cache"),
        strict_policy="local",
        strict_callback=lambda *_args: {"status": "full_strict_failed"},
        force_audit=True,
    )
    assert payload["execution_status"] == "attempted_unresolved"
    assert payload["N_true"] is None


def test_full_strict_failure_above_required_guard_preserves_exact_prefix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    set_delta_failure(2)
    monkeypatch.setattr(prefix, "IncrementalModelSession", FakeSession)
    monkeypatch.setattr(prefix, "_stage_assessments", resolved_assessments)
    payload = prefix.run_staged_point(
        point(),
        cache=prefix.PartialPointCache(tmp_path / "cache"),
        strict_policy="local",
        strict_callback=lambda *_args: {
            "status": "required_prefix_strict_pass_upper_audit_incomplete",
            "required_prefix_strict_status": "pass",
            "optional_upper_spectrum_full_audit_status": "incomplete",
        },
        force_audit=True,
    )
    assert payload["execution_status"] == "resolved_prefix_early_stop"
    assert payload["N_true"] == 1
    assert payload["required_prefix_strict_status"] == "pass"
    assert payload["optional_upper_spectrum_full_audit_status"] == "incomplete"


def test_strict_failure_at_or_below_guard_remains_unresolved(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    set_delta_failure(2)
    monkeypatch.setattr(prefix, "IncrementalModelSession", FakeSession)
    monkeypatch.setattr(prefix, "_stage_assessments", resolved_assessments)
    payload = prefix.run_staged_point(
        point(),
        cache=prefix.PartialPointCache(tmp_path / "cache"),
        strict_policy="local",
        strict_callback=lambda *_args: {
            "status": "full_strict_failed",
            "required_prefix_strict_status": "fail",
            "optional_upper_spectrum_full_audit_status": "incomplete",
        },
        force_audit=True,
    )
    assert payload["execution_status"] == "attempted_unresolved"
    assert payload["N_true"] is None


def test_unresolved_cluster_through_guard_blocks_required_prefix_strict() -> None:
    branches = [
        SimpleNamespace(
            Lambda=float(index),
            sigma_1=0.0,
            sigma_ratio=0.0,
            branch_id=f"b{index}",
        )
        for index in range(1, 5)
    ]
    result = SimpleNamespace(
        branches=branches,
        primary_vs_verification=[{"within_tolerance": True}] * 4,
        settings=SimpleNamespace(
            sigma_accept=1.0e-6,
            sigma_ratio_accept=1.0e-3,
            seed_half_width=1.0e-3,
        ),
        steps=[
            SimpleNamespace(
                accepted=True,
                clusters=[SimpleNamespace(resolved=False, branch_ids=("b2",))],
            )
        ],
        guard=SimpleNamespace(unmatched_candidates=(), unresolved_intervals=()),
        k10_guard_resolved=False,
        exclusion_reason="cluster_unresolved_above_or_within_scan",
    )
    assessment = prefix.assess_branch_strict_through_guard(
        result,
        tuple(float(index) for index in range(1, 5)),
        4,
        1.0e-8,
    )
    assert assessment["agreement"] is False
    assert assessment["cluster_resolved"] is False


def test_four_known_smoke_failures_are_registered_as_full_path_controls() -> None:
    points = workflow.build_manifest(smoke=True)
    controls = [item for item in points if prefix.is_known_full_path_unresolved_smoke(item)]
    assert len(controls) == 4


def test_not_needed_after_first_failure_is_not_unresolved() -> None:
    rows = [{"quality_status": "resolved", "case_id": "a"}]
    assert analyzer.build_unresolved(rows) == []


def test_not_attempted_is_separate_from_unresolved() -> None:
    row = {
        "quality_status": "not_attempted",
        "case_id": "a",
        "case_identity": "i",
        "grid_group": "base",
        "epsilon_0": 0.02,
        "beta_deg": 0.0,
        "mu": 0.0,
        "eta": 0.0,
        "execution_status": "not_attempted",
        "unresolved_reason": "not_attempted",
    }
    assert analyzer.build_unresolved([row]) == []
    assert len(analyzer.build_not_attempted([row])) == 1


def test_envelope_saturation_marks_only_upper_envelope_complete() -> None:
    status = prefix.envelope_saturation_status(10)
    assert status["complete_for_upper_envelope"] is True
    assert status["complete_for_distribution"] is False
    assert status["remaining_status"] == "not_needed_after_envelope_saturation"


def test_full_and_prefix_csv_N_true_semantics_match() -> None:
    p = point()
    manifest = [{key: str(value) for key, value in p.manifest_row().items()}]
    rows = []
    for model in complete.SUPPORTED_MODELS:
        for index in range(1, 7):
            rows.append(
                {
                    "case_id": p.case_id,
                    "model": model,
                    "sorted_index": str(index),
                    "Lambda": str(float(index) * (1.1**0.5 if model == complete.MODEL_EB and index == 5 else 1.0)),
                    "solver_status": "resolved",
                    "root_quality": "pass",
                    "execution_status": "resolved_prefix_early_stop",
                    "N_true_cached": "4",
                    "first_failed_delta_f": "0.1",
                    "prefix_guard_status": "prefix_guard_resolved",
                    "prefix_guard_resolved_through": "5",
                    "full_K10_guard_status": "not_needed_after_first_failure",
                    "strict_verification_status": "local_prefix_count_and_guard_pass",
                }
            )
    summary = analyzer.build_case_summary(manifest, rows)
    assert summary[0]["N_true"] == 4
    assert summary[0]["quality_status"] == "resolved"


def test_incompatible_cache_schema_is_rejected(tmp_path: Path) -> None:
    p = point()
    cache = prefix.PartialPointCache(tmp_path)
    path = cache.path(p, strategy="paired", strict_policy="local")
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        json.dump({"cache_schema_version": "old"}, handle)
    assert cache.load(p, strategy="paired", strict_policy="local") is None
    assert cache.last_load_status == "incompatible_schema"


def test_benchmark_cache_isolated_from_main() -> None:
    assert "optimization_benchmark" in str(benchmark.OUTPUT_DIR)
    assert benchmark.OUTPUT_DIR / "cache" != workflow.MAIN_OUTPUT_DIR / "cache"


def test_validation_manifest_has_24_unique_difficult_cases() -> None:
    cases = benchmark.build_validation_manifest()
    assert len(cases) == 24
    assert len({point.case_id for _category, point in cases}) == 24
    categories = {category for category, _point in cases}
    assert sum(category.startswith("smoke_unresolved") for category in categories) == 4
    assert {"regression_S3_12", "regression_S3_14"}.issubset(categories)


def test_control_sample_is_deterministic_and_at_least_five_percent() -> None:
    rows = [
        {
            "case_id": index,
            "epsilon_0": index % 8,
            "beta_deg": index % 7,
            "mu": index % 5,
            "eta": index % 3,
            "N_true": index % 11,
            "thin_0p1_flag": index % 2,
            "cutoff_regime": index % 4,
        }
        for index in range(101)
    ]
    left = prefix.deterministic_control_sample(rows)
    right = prefix.deterministic_control_sample(rows)
    assert left == right
    assert len(left) >= 6


def test_future_full_k10_policy_selects_argmax_threshold_and_sample() -> None:
    reasons = prefix.full_k10_control_reasons(
        {
            "case_id": "case",
            "N_true": 7,
            "deltas": {"8": 0.1001},
            "execution_status": "resolved_prefix_early_stop",
        },
        n_up_raw=7,
        sampled_case_ids=("case",),
    )
    assert {"defines_N_up_raw", "near_delta_threshold", "stratified_five_percent_control"}.issubset(reasons)


def test_prefix_cache_identity_includes_strategy_policy_and_schema() -> None:
    p = point()
    paired = prefix.prefix_cache_identity(p, strategy="paired", strict_policy="local")
    other = prefix.prefix_cache_identity(
        p, strategy="full-eb-progressive-timo", strict_policy="auto"
    )
    assert paired != other
    assert paired["cache_schema_version"] == prefix.PREFIX_CACHE_SCHEMA_VERSION


def test_full_k10_remains_default_cli_mode() -> None:
    args = runner.parse_args([])
    assert args.prefix_until_failure is False
    assert args.full_k10 is False


def test_prefix_cli_exposes_both_strategies_and_strict_policies() -> None:
    args = runner.parse_args(
        [
            "--prefix-until-failure",
            "--prefix-strategy",
            "full-eb-progressive-timo",
            "--strict-policy",
            "local",
        ]
    )
    assert args.prefix_until_failure is True
    assert args.prefix_strategy == "full-eb-progressive-timo"
    assert args.strict_policy == "local"


def test_envelope_only_requires_prefix_mode() -> None:
    with pytest.raises(SystemExit):
        runner.parse_args(["--envelope-only"])


def status_payload(status: str, reason: str = "", *, n_true: int | None = None, failed: int | None = None):
    return {
        "execution_status": status,
        "unresolved_reason": reason,
        "N_true": n_true,
        "first_failed_mode": failed,
        "models": {},
    }


def summary_case_row(mode: str, *, early: bool, failed: int | None, category: str = "ordinary"):
    return {
        "run_mode": mode,
        "case_id": "case",
        "category": category,
        "execution_status": "resolved_full_K10" if mode == benchmark.MODE_ORDER[0] else "resolved_prefix_early_stop",
        "N_true": 4,
        "first_failed_mode": failed if failed is not None else "",
        "early_stop_used": early,
        "wall_seconds": 1.0,
        "peak_memory_bytes": 0,
        "EB_primary_seconds": 0.1,
        "Timo_primary_seconds": 0.2,
        "EB_strict_seconds": 0.3,
        "Timo_strict_seconds": 0.4,
        "preparation_seconds": 0.0,
        "branch_strict_seconds": 0.0,
        "branch_strict_cache_hits": 0,
        "branch_strict_cache_misses": 0,
        "EB_root_count": 11,
        "Timo_root_count": 11,
        "guard_root_count": 12,
        "cache_hits": 0,
        "warm_cache_seconds": "",
    }


def test_resolved_execution_status_agreement_is_true() -> None:
    payloads = [status_payload("resolved_full_K10"), status_payload("resolved_prefix_early_stop")]
    assert benchmark.execution_status_agreement(payloads) is True


def test_same_unresolved_execution_reason_agrees() -> None:
    payloads = [
        status_payload("attempted_unresolved", "strict_failure_or_disagreement"),
        status_payload("attempted_error", "strict_failure_or_disagreement"),
    ]
    assert benchmark.execution_status_agreement(payloads) is True


def test_resolved_and_unresolved_execution_status_disagree() -> None:
    payloads = [
        status_payload("resolved_full_K10"),
        status_payload("attempted_unresolved", "root_count_failure"),
    ]
    assert benchmark.execution_status_agreement(payloads) is False


def test_full_mode_is_never_counted_as_early_stop() -> None:
    rows = [summary_case_row(benchmark.MODE_ORDER[0], early=True, failed=5)]
    summary = benchmark._run_summary(rows)
    assert summary[0]["early_stopped_points"] == 0


def test_first_failure_is_counted_separately_from_early_stop() -> None:
    rows = [summary_case_row(benchmark.MODE_ORDER[0], early=True, failed=5)]
    summary = benchmark._run_summary(rows)
    assert summary[0]["points_with_first_failure"] == 1


def test_equivalence_and_readiness_gates_are_independent() -> None:
    rows = [
        {"equivalence_case_status": "PASS", "solver_readiness_case_status": "READY"},
        {
            "equivalence_case_status": "NOT_EVALUABLE_REFERENCE_UNRESOLVED",
            "solver_readiness_case_status": "IMPLEMENTATION_LIMIT",
        },
    ]
    assert benchmark.calculate_gates(rows) == ("PASS", "BLOCKED_BY_UNRESOLVED_REFERENCE")


def test_unresolved_reference_cannot_create_false_optimized_pass() -> None:
    category, p = benchmark.build_validation_manifest()[2]
    results = {
        benchmark.MODE_ORDER[0]: {p.case_id: status_payload("attempted_unresolved", "basis regime")},
        benchmark.MODE_ORDER[1]: {p.case_id: status_payload("resolved_prefix_early_stop", n_true=2, failed=3)},
        benchmark.MODE_ORDER[2]: {p.case_id: status_payload("attempted_unresolved", "basis regime")},
    }
    row = benchmark._accuracy_rows([(category, p)], results, [])[0]
    assert row["equivalence_case_status"] == "NOT_EVALUABLE_REFERENCE_UNRESOLVED"
    assert row["local_execution_status"] == "not_accepted_reference_unresolved"
    assert row["local_N_true"] == ""


def test_targeted_audit_implementation_limit_is_preserved_without_basis_word() -> None:
    category, p = benchmark.build_validation_manifest()[4]
    full = status_payload("attempted_unresolved", "unsupported two-trigonometric-pair regime")
    full["implementation_limit"] = True
    results = {
        benchmark.MODE_ORDER[0]: {p.case_id: full},
        benchmark.MODE_ORDER[1]: {},
        benchmark.MODE_ORDER[2]: {},
    }
    row = benchmark._accuracy_rows([(category, p)], results, [])[0]
    assert row["solver_readiness_case_status"] == "IMPLEMENTATION_LIMIT"


def comparison_rows(case_id: str) -> list[dict[str, object]]:
    rows = []
    for mode in benchmark.OPTIMIZED_MODES:
        rows.append({"case_id": case_id, "run_mode": mode, "comparison_status": "PASS"})
        rows.append(
            {
                "case_id": case_id,
                "run_mode": mode,
                "comparison_status": "not_needed_after_first_failure",
            }
        )
    return rows


def test_N_true_mismatch_closes_equivalence_gate() -> None:
    category, p = benchmark.build_validation_manifest()[6]
    results = {
        benchmark.MODE_ORDER[0]: {p.case_id: status_payload("resolved_full_K10", n_true=4, failed=5)},
        benchmark.MODE_ORDER[1]: {p.case_id: status_payload("resolved_prefix_early_stop", n_true=3, failed=5)},
        benchmark.MODE_ORDER[2]: {p.case_id: status_payload("resolved_prefix_early_stop", n_true=4, failed=5)},
    }
    row = benchmark._accuracy_rows([(category, p)], results, comparison_rows(p.case_id))[0]
    assert row["equivalence_case_status"] == "FAIL"
    assert benchmark.calculate_gates([row])[0] == "FAIL_DISAGREEMENT"


def test_not_needed_roots_do_not_count_as_disagreement() -> None:
    category, p = benchmark.build_validation_manifest()[6]
    results = {
        benchmark.MODE_ORDER[0]: {p.case_id: status_payload("resolved_full_K10", n_true=4, failed=5)},
        benchmark.MODE_ORDER[1]: {p.case_id: status_payload("resolved_prefix_early_stop", n_true=4, failed=5)},
        benchmark.MODE_ORDER[2]: {p.case_id: status_payload("resolved_prefix_early_stop", n_true=4, failed=5)},
    }
    row = benchmark._accuracy_rows([(category, p)], results, comparison_rows(p.case_id))[0]
    assert row["root_agreement"] is True
    assert row["equivalence_case_status"] == "PASS"


def test_root_comparison_uses_existing_project_tolerance() -> None:
    assert complete.DEFAULT_ROOT_MATCH_TOL == workflow.primary_settings().root_match_tol


def test_smoke_rows_are_excluded_from_runtime_summary() -> None:
    ordinary = summary_case_row(benchmark.MODE_ORDER[0], early=False, failed=5)
    smoke = summary_case_row(
        benchmark.MODE_ORDER[0], early=False, failed=5, category="smoke_unresolved_control"
    )
    smoke["case_id"] = "smoke"
    summary = benchmark._run_summary([ordinary, smoke])
    assert summary[0]["case_count"] == 1


def test_report_distinguishes_incomplete_and_implementation_limit(tmp_path: Path) -> None:
    rows = [
        {
            "equivalence_case_status": "NOT_ATTEMPTED",
            "solver_readiness_case_status": "IMPLEMENTATION_LIMIT",
        }
    ]
    path = benchmark.write_report(tmp_path, 1, 0, [], rows, [], "paired")
    text = path.read_text(encoding="utf-8")
    assert "INCOMPLETE" in text
    assert "BLOCKED_BY_UNRESOLVED_REFERENCE" in text
    assert "not a contradiction" in text


def test_full_local_auto_cache_directories_do_not_conflict(tmp_path: Path) -> None:
    roots = [tmp_path / mode / "partial" for mode in benchmark.MODE_ORDER]
    assert len(set(roots)) == 3


def test_postprocess_only_is_zero_solve_and_deterministic(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def forbidden(*_args, **_kwargs):
        raise AssertionError("root solver called during postprocess-only")

    monkeypatch.setattr(prefix, "run_staged_point", forbidden)
    args = ["--output-dir", str(tmp_path), "--postprocess-only"]
    first = benchmark.main(args)
    snapshots = {
        name: (tmp_path / name).read_bytes()
        for name in ("benchmark_cases.csv", "benchmark_runs.csv", "root_comparison.csv", "accuracy_gate.csv")
    }
    second = benchmark.main(args)
    assert first["root_calculations"] == second["root_calculations"] == 0
    assert snapshots == {name: (tmp_path / name).read_bytes() for name in snapshots}


def test_case_row_order_is_deterministic() -> None:
    cases = benchmark.build_validation_manifest()[:2]
    rows = [
        {"run_mode": benchmark.MODE_ORDER[2], "case_id": cases[1][1].case_id},
        {"run_mode": benchmark.MODE_ORDER[1], "case_id": cases[0][1].case_id},
        {"run_mode": benchmark.MODE_ORDER[0], "case_id": cases[0][1].case_id},
    ]
    ordered = benchmark._deterministic_case_rows(rows, cases)
    assert [(row["case_id"], row["run_mode"]) for row in ordered] == [
        (cases[0][1].case_id, benchmark.MODE_ORDER[0]),
        (cases[0][1].case_id, benchmark.MODE_ORDER[1]),
        (cases[1][1].case_id, benchmark.MODE_ORDER[2]),
    ]


def test_saved_S3_regression_values_remain_unchanged() -> None:
    path = benchmark.OUTPUT_DIR / "benchmark_cases.csv"
    if not path.exists():
        pytest.skip("saved benchmark artifact is not available")
    rows = workflow.read_csv(path)
    regressions = [row for row in rows if row.get("category") in {"regression_S3_12", "regression_S3_14"}]
    assert regressions
    assert {row["N_true"] for row in regressions} == {"4"}
    assert {row["first_failed_mode"] for row in regressions} == {"5"}
