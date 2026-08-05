from __future__ import annotations

from dataclasses import asdict
import math
from pathlib import Path

import pytest

from scripts.analysis.thickness_mismatch.audits import (
    run_article_epsilon_upper_envelope_grid as runner,
)
from scripts.analysis.thickness_mismatch.benchmarks import (
    benchmark_article_epsilon_prefix_parallelism as parallel_benchmark,
)
from scripts.analysis.thickness_mismatch.postprocess import (
    analyze_article_epsilon_upper_envelope as analyzer,
)
from scripts.lib import article_epsilon_prefix_optimization as prefix
from scripts.lib import article_epsilon_upper_envelope as workflow
from scripts.lib import branch_informed_spectrum_continuation as branch
from scripts.lib import general_spectrum_completeness as complete


@pytest.fixture(scope="module")
def manifest() -> list[workflow.GridPoint]:
    return workflow.build_manifest()


def test_manifest_ids_and_full_precision_identities_are_unique(manifest: list[workflow.GridPoint]) -> None:
    assert len({point.case_id for point in manifest}) == 1554
    assert len({point.case_identity for point in manifest}) == 1554
    assert all(float(point.epsilon_0).hex() in point.case_identity for point in manifest)


def test_expected_grid_sizes(manifest: list[workflow.GridPoint]) -> None:
    assert workflow.group_counts(manifest) == {
        "base": 1400,
        "low_angle": 144,
        "s3_14_sweep": 8,
        "regression": 2,
    }
    assert len(workflow.select_points(manifest, base_only=True)) == 1400
    assert len(workflow.select_points(manifest, low_angle_only=True)) == 144
    assert len(workflow.select_points(manifest, regressions_only=True)) == 2
    assert len(workflow.build_manifest(smoke=True)) == 16


def test_tau_formulas_and_mass_preservation(manifest: list[workflow.GridPoint]) -> None:
    for point in manifest:
        denominator = math.sqrt(1.0 + 2.0 * point.mu * point.eta + point.eta**2)
        assert point.tau1 == pytest.approx((1.0 - point.eta) / denominator, abs=2.0e-15)
        assert point.tau2 == pytest.approx((1.0 + point.eta) / denominator, abs=2.0e-15)
        assert (1.0 - point.mu) * point.tau1**2 + (1.0 + point.mu) * point.tau2**2 == pytest.approx(
            2.0, abs=5.0e-13
        )


@pytest.mark.parametrize(
    "deltas,expected",
    [
        ([0.01] * 10, 10),
        ([0.11] + [0.01] * 9, 0),
        ([0.01, 0.11, 0.01] + [0.01] * 7, 1),
        ([0.10] * 4 + [0.1000001] + [0.01] * 5, 4),
    ],
)
def test_true_safe_prefix_is_contiguous(deltas: list[float], expected: int) -> None:
    assert workflow.true_safe_prefix(deltas) == expected


def test_squared_frequency_delta_uses_timoshenko_squared_denominator() -> None:
    # |2^2 - 1^2| / 1^2 = 3; an EB denominator would give 0.75.
    assert workflow.squared_frequency_delta(2.0, 1.0) == pytest.approx(3.0)


def test_suffix_max_monotone_envelope() -> None:
    assert workflow.suffix_max([4, 3, 5, 2]) == [5, 5, 5, 2]


def test_thin_flag_is_diagnostic_only_and_high_s_points_remain(manifest: list[workflow.GridPoint]) -> None:
    high = [point for point in manifest if point.claim_eligible and not point.thin_0p1_flag]
    assert high
    assert any(point.s_max > 0.1 for point in high)
    assert len([point for point in manifest if point.claim_eligible]) == 1552
    assert len(workflow.select_points(manifest)) == 1554


def test_cache_identity_depends_on_full_precision_geometry_and_solver_configuration(
    manifest: list[workflow.GridPoint],
) -> None:
    point = manifest[0]
    near = workflow._point(
        point.grid_group,
        math.nextafter(point.epsilon_0, math.inf),
        point.beta_deg,
        point.mu,
        point.eta,
    )
    baseline = workflow.case_cache_identity(point)
    changed_geometry = workflow.case_cache_identity(near)
    assert baseline != changed_geometry
    changed_settings = workflow.solver_configuration()
    changed_settings["primary_settings"] = {
        **changed_settings["primary_settings"],  # type: ignore[dict-item]
        "scan_step": 0.005,
    }
    assert workflow.case_cache_identity(point, changed_settings) != baseline


def test_smoke_and_main_outputs_do_not_overlap() -> None:
    assert workflow.SMOKE_OUTPUT_DIR != workflow.MAIN_OUTPUT_DIR
    assert workflow.MAIN_OUTPUT_DIR not in workflow.SMOKE_OUTPUT_DIR.parents
    assert runner.COARSE_GRID_OUTPUT_DIR != workflow.MAIN_OUTPUT_DIR


def test_coarse_runner_accepts_only_bounded_prefix_workers() -> None:
    assert runner.parse_args(["--workers", "1"]).workers == 1
    assert runner.parse_args(["--prefix-until-failure", "--workers", "2"]).workers == 2
    assert runner.parse_args(["--prefix-until-failure", "--workers", "4"]).workers == 4
    with pytest.raises(SystemExit):
        runner.parse_args(["--prefix-until-failure", "--workers", "3"])
    with pytest.raises(SystemExit):
        runner.parse_args(["--full-k10", "--workers", "2"])


def test_main_pass_resume_arguments_are_orchestration_only(tmp_path: Path) -> None:
    deferred = tmp_path / "deferred.csv"
    args = runner.parse_args(
        [
            "--prefix-until-failure",
            "--workers",
            "4",
            "--reuse-cache",
            "--main-pass-only",
            "--skip-existing-unresolved",
            "--skip-interrupted",
            "--defer-case-list",
            str(deferred),
        ]
    )
    assert args.main_pass_only is True
    assert args.skip_existing_unresolved is True
    assert args.skip_interrupted is True
    assert args.defer_case_list == deferred
    assert workflow.solver_configuration()["workflow_version"] == workflow.WORKFLOW_VERSION


def test_main_pass_selects_only_not_attempted_non_deferred(manifest: list[workflow.GridPoint]) -> None:
    points = manifest[:5]
    statuses = {
        points[0].case_id: "resolved_prefix_early_stop",
        points[1].case_id: "resolved_full_K10",
        points[2].case_id: "attempted_unresolved",
        points[3].case_id: "interrupted_incomplete",
        points[4].case_id: "not_attempted",
    }
    assert runner._main_pass_selection(points, statuses, set()) == [points[4]]
    assert runner._main_pass_selection(points, statuses, {points[4].case_id}) == []


@pytest.mark.parametrize(
    "resolved,n_up,upper_complete,distribution_complete,upper_status",
    [
        (194, 7, True, True, "exact"),
        (190, 10, True, False, "exact_by_K_saturation"),
        (190, 8, False, False, "provisional_lower_bound"),
    ],
)
def test_main_pass_completeness_and_k_saturation(
    resolved: int,
    n_up: int,
    upper_complete: bool,
    distribution_complete: bool,
    upper_status: str,
) -> None:
    result = analyzer.main_pass_completeness(
        intended_count=194,
        resolved_count=resolved,
        n_up_observed=n_up,
    )
    assert result["complete_for_upper_envelope"] is upper_complete
    assert result["complete_for_distribution"] is distribution_complete
    assert result["upper_envelope_status"] == upper_status


def _synthetic_spectra_rows(points: list[workflow.GridPoint]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for point in points:
        for model in complete.SUPPORTED_MODELS:
            for index in range(1, 12):
                value = 1.0 + index
                rows.append(
                    {
                        "case_id": point.case_id,
                        "case_identity": point.case_identity,
                        "grid_group": point.grid_group,
                        "regression_label": point.regression_label,
                        "claim_eligible": point.claim_eligible,
                        "epsilon_0": point.epsilon_0,
                        "beta_deg": point.beta_deg,
                        "mu": point.mu,
                        "eta": point.eta,
                        "model": model,
                        "sorted_index": index,
                        "root_role": "reported_K10" if index <= 10 else "K10_right_guard",
                        "Lambda": value,
                        "Lambda_squared": value**2,
                        "solver_status": "resolved",
                        "root_quality": "pass",
                        "guard_status": "root11_resolved",
                        "strict_verification_status": "not_triggered_primary_independent_pass",
                        "strict_trigger_reasons": "",
                        "cluster_id": "",
                        "cluster_member_index": 1,
                        "cluster_size": 1,
                        "detected_nullity": 1,
                        "track_multiplicity": 1,
                        "multiplicity_status": "simple_root",
                        "branch_id": "",
                        "parent_family": "",
                        "branch_reordered": False,
                        "detection_sources": "synthetic",
                        "sigma_1": 0.0,
                        "sigma_ratio": 0.0,
                        "source_path": "synthetic",
                        "cache_source_path": "synthetic",
                        "omega": value,
                        "omega_over_cutoff_1": 0.01,
                        "omega_over_cutoff_2": 0.01,
                        "max_cutoff_ratio": 0.01,
                    }
                )
    return rows


def test_postprocess_only_performs_no_root_solver_call(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    points = workflow.build_manifest(smoke=True)
    workflow.write_csv(tmp_path / "grid_manifest.csv", [point.manifest_row() for point in points], workflow.MANIFEST_FIELDS)
    workflow.write_csv(tmp_path / "spectra_long.csv", _synthetic_spectra_rows(points), runner.SPECTRA_FIELDS)
    monkeypatch.setattr(
        complete,
        "resolve_general_spectrum",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("root solver called")),
    )
    monkeypatch.setattr(
        branch,
        "resolve_branch_spectrum",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("strict root solver called")),
    )
    result = runner.postprocess_only(tmp_path)
    assert result["root_calculations"] == 0
    assert result["resolved_count"] == 16
    assert (tmp_path / "epsilon_summary.csv").exists()


def test_all_group_does_not_filter_high_s_rows() -> None:
    rows = [
        {
            "epsilon_0": 0.02,
            "grid_group": "base",
            "claim_eligible": True,
            "thin_0p1_flag": True,
            "quality_status": "resolved",
            "N_true": 4,
        },
        {
            "epsilon_0": 0.02,
            "grid_group": "base",
            "claim_eligible": True,
            "thin_0p1_flag": False,
            "quality_status": "resolved",
            "N_true": 7,
        },
    ]
    summary = analyzer.build_epsilon_summary(rows)
    all_row = next(row for row in summary if row["grid_group"] == "all")
    assert all_row["intended_case_count"] == 2
    assert all_row["N_up_raw"] == 7


def test_full_controls_gate_is_independent_of_resolved_prefix_rows() -> None:
    rows = [
        {
            "case_id": "a",
            "epsilon_0": 0.02,
            "grid_group": "base",
            "claim_eligible": True,
            "thin_0p1_flag": True,
            "quality_status": "resolved",
            "N_true": 4,
        }
    ]
    blocked = analyzer.build_epsilon_summary(rows, [], controls_required=True)
    assert next(row for row in blocked if row["grid_group"] == "all")["complete_for_claim"] is False
    passed = analyzer.build_epsilon_summary(
        rows,
        [{"case_id": "a", "comparison_status": "PASS"}],
        controls_required=True,
    )
    assert next(row for row in passed if row["grid_group"] == "all")["complete_for_claim"] is True


def test_stratified_control_sample_is_at_least_five_percent_and_covers_levels() -> None:
    rows = [
        {
            "case_id": f"case_{index:04d}",
            "epsilon_0": (0.01, 0.02)[index % 2],
            "beta_deg": (0.0, 90.0)[index % 2],
            "mu": (0.0, 0.9)[(index // 2) % 2],
            "eta": (-0.5, 0.5)[(index // 4) % 2],
            "N_true": index % 11,
            "thin_0p1_flag": bool(index % 2),
            "basis_regimes": ("mixed", "mixed;two_trig")[index % 2],
            "first_failure_stratum": ("early_failure", "middle_prefix", "late_prefix", "N_true_10")[index % 4],
        }
        for index in range(1552)
    ]
    selected = runner._sample_controls(rows, population_size=1552)
    assert len(selected) >= 78
    selected_rows = [row for row in rows if row["case_id"] in selected]
    for field in (
        "epsilon_0",
        "beta_deg",
        "mu",
        "eta",
        "N_true",
        "thin_0p1_flag",
        "basis_regimes",
        "first_failure_stratum",
    ):
        assert {row[field] for row in selected_rows} == {row[field] for row in rows}


def test_prefix_full_comparison_ignores_unneeded_roots_above_guard() -> None:
    def payload(extra: float) -> dict[str, object]:
        models = {}
        for model in complete.SUPPORTED_MODELS:
            roots = [
                {
                    "sorted_index": index,
                    "Lambda": float(index) if index <= 4 else float(index) + extra,
                    "cluster_id": f"root_cluster_{index:03d}",
                    "cluster_member_index": 1,
                    "cluster_size": 1,
                    "detected_nullity": 1,
                    "multiplicity_status": "simple_root",
                }
                for index in range(1, 12)
            ]
            models[model] = {"roots": roots}
        return {
            "case_status": "resolved",
            "execution_status": "resolved_prefix_early_stop" if extra == 0.0 else "resolved_full_K10",
            "N_true": 2,
            "first_failed_mode": 3,
            "prefix_guard_status": "prefix_guard_resolved",
            "full_K10_guard_status": "full_K10_guard_resolved",
            "models": models,
            "accepted_models": models,
            "strict": {"status": "full_strict_pass"},
        }

    comparison = runner._compare_prefix_full(payload(0.0), payload(100.0))
    assert comparison["comparison_status"] == "PASS"
    assert comparison["compared_root_count_EB"] == 4
    assert comparison["compared_root_count_Timo"] == 4


def test_parallel_benchmark_uses_eight_indivisible_geometry_chains() -> None:
    points = parallel_benchmark.benchmark_points()
    chains = parallel_benchmark.group_chains(points)
    assert len(chains) == 8
    assert all(len(chain) == 8 for chain in chains)
    for chain in chains:
        assert len({(point.beta_deg, point.mu, point.eta) for point in chain}) == 1
        assert tuple(point.epsilon_0 for point in chain) == workflow.EPSILON_VALUES


def test_parallel_benchmark_cache_is_isolated_from_main_grid() -> None:
    assert parallel_benchmark.OUTPUT_DIR != parallel_benchmark.MAIN_OUTPUT_DIR
    assert parallel_benchmark.MAIN_OUTPUT_DIR not in parallel_benchmark.OUTPUT_DIR.parents
    for workers in parallel_benchmark.WORKER_COUNTS:
        cache = parallel_benchmark.OUTPUT_DIR / "cache" / f"workers_{workers}"
        assert parallel_benchmark.MAIN_OUTPUT_DIR not in cache.parents


def test_parallel_comparison_requires_identical_worker_results() -> None:
    point = parallel_benchmark.benchmark_points()[0]
    rows = []
    for workers in parallel_benchmark.WORKER_COUNTS:
        rows.append(
            {
                "workers": workers,
                "case_id": point.case_id,
                "epsilon_0": point.epsilon_0,
                "beta_deg": point.beta_deg,
                "mu": point.mu,
                "eta": point.eta,
                "execution_status": "resolved_prefix_early_stop",
                "N_true": 4,
                "first_failed_mode": 5,
                "EB_roots_json": "[1.0,2.0,3.0,4.0,5.0,6.0]",
                "Timoshenko_roots_json": "[1.0,2.0,3.0,4.0,5.0,6.0]",
            }
        )
    comparison = parallel_benchmark._compare_rows(rows)
    assert comparison[0]["comparison_status"] == "PASS"
    rows[-1]["N_true"] = 3
    assert parallel_benchmark._compare_rows(rows)[0]["comparison_status"] == "FAIL"


def test_atomic_prefix_cache_survives_interrupted_temporary_write(tmp_path: Path) -> None:
    point = parallel_benchmark.benchmark_points()[0]
    cache = prefix.PartialPointCache(tmp_path / "isolated", reuse_cache=True)
    payload = prefix.fresh_point_payload(point, strategy="paired", strict_policy="auto")
    payload["execution_status"] = "resolved_prefix_early_stop"
    payload["N_true"] = 4
    path = cache.save(point, payload, strategy="paired", strict_policy="auto")
    path.with_name(path.name + ".tmp").write_bytes(b"interrupted worker payload")
    loaded = cache.load(point, strategy="paired", strict_policy="auto")
    assert loaded is not None
    assert loaded["N_true"] == 4


def test_partial_postprocess_marks_interrupted_cache_without_root_solves(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    point = workflow.build_manifest(smoke=True)[0]
    monkeypatch.setattr(workflow, "build_manifest", lambda *args, **kwargs: [point])
    workflow.write_csv(tmp_path / "runtime_by_case.csv", [], runner.RUNTIME_FIELDS)
    cache = prefix.PartialPointCache(tmp_path / "cache" / "prefix", reuse_cache=True)
    payload = prefix.fresh_point_payload(point, strategy="paired", strict_policy="auto")
    payload["execution_status"] = "attempted_unresolved"
    payload["unresolved_reason"] = "execution_in_progress_or_interrupted"
    cache.save(point, payload, strategy="paired", strict_policy="auto")
    monkeypatch.setattr(
        complete,
        "resolve_general_spectrum",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("root solver called")),
    )
    monkeypatch.setattr(
        branch,
        "resolve_branch_spectrum",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(AssertionError("strict solver called")),
    )
    result = analyzer.analyze_partial_cache(tmp_path)
    assert result["root_calculations"] == 0
    assert result["interrupted_incomplete_count"] == 1
    row = workflow.read_csv(tmp_path / "partial_case_summary.csv")[0]
    assert row["execution_status"] == "interrupted_incomplete"
