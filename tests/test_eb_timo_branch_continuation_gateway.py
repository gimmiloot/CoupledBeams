from __future__ import annotations

from dataclasses import replace
import inspect
from pathlib import Path

import numpy as np
import pytest

from scripts.analysis.thickness_mismatch.audits import (
    audit_eb_timo_branch_continuation_gateway as gateway,
)
from scripts.analysis.thickness_mismatch.audits import run_eb_epsilon_apriori_pilot as pilot
from scripts.lib import branch_informed_spectrum_continuation as branch
from scripts.lib import general_spectrum_completeness as complete


GEOMETRY_0 = complete.Geometry(0.015625, 0.0, 0.0, 0.0)
GEOMETRY_SMALL = complete.Geometry(0.015625, 0.1, 0.0, 0.0)


@pytest.fixture(scope="module")
def beta0_results() -> dict[str, branch.BranchContinuationResult]:
    return {
        model: branch.resolve_branch_spectrum(model, GEOMETRY_0)
        for model in complete.SUPPORTED_MODELS
    }


@pytest.fixture(scope="module")
def small_angle_results() -> dict[str, branch.BranchContinuationResult]:
    return {
        model: branch.resolve_branch_spectrum(model, GEOMETRY_SMALL)
        for model in complete.SUPPORTED_MODELS
    }


def test_eb_block_map_uses_production_ordering() -> None:
    item = branch.BETA0_BLOCK_MAPS[complete.MODEL_EB]
    assert item.bending_rows == (0, 2, 3, 4)
    assert item.bending_columns == (0, 1, 2, 3)
    assert item.axial_rows == (1, 5)
    assert item.axial_columns == (4, 5)


def test_timo_block_map_uses_production_ordering() -> None:
    item = branch.BETA0_BLOCK_MAPS[complete.MODEL_TIMO]
    assert item.bending_rows == (0, 2, 3, 4)
    assert item.bending_columns == (0, 1, 3, 4)
    assert item.axial_rows == (1, 5)
    assert item.axial_columns == (2, 5)


@pytest.mark.parametrize("model", complete.SUPPORTED_MODELS)
@pytest.mark.parametrize("mu,eta", [(0.0, 0.0), (0.4, 0.1), (-0.2, -0.1)])
def test_beta0_off_block_norm_is_small(model: str, mu: float, eta: float) -> None:
    geometry = complete.Geometry(0.02, 0.0, mu, eta)
    matrix = complete.model_matrix_provider(model, geometry)(5.0)
    assert branch._off_block_max_abs(matrix, branch.BETA0_BLOCK_MAPS[model]) <= 1.0e-12


@pytest.mark.parametrize("model", complete.SUPPORTED_MODELS)
def test_parent_families_remain_separate(
    beta0_results: dict[str, branch.BranchContinuationResult], model: str
) -> None:
    families = {item.family for item in beta0_results[model].parent_branches}
    assert families == {"axial", "bending"}


def test_cross_family_multiplicity_is_not_deduplicated() -> None:
    bending = complete.fixed_fixed_eb_roots(1)[0] / 2.0
    epsilon = np.pi / (2.0 * bending**2)
    result = branch.resolve_branch_spectrum(
        complete.MODEL_EB,
        complete.Geometry(float(epsilon), 0.0, 0.0, 0.0),
    )
    close = [item for item in result.parent_branches if abs(item.Lambda - bending) < 2.0e-8]
    assert {item.family for item in close} == {"axial", "bending"}


@pytest.mark.parametrize("model", complete.SUPPORTED_MODELS)
def test_straight_parent_spectrum_matches_factorized_oracle(
    beta0_results: dict[str, branch.BranchContinuationResult], model: str
) -> None:
    result = beta0_results[model]
    oracle = complete.straight_oracle_values(model, GEOMETRY_0, 12)
    assert np.allclose(result.values[:12], oracle, rtol=0.0, atol=2.0e-7)
    assert result.oracle_agreement is True


def test_seed_point_is_not_directly_accepted_in_general_helper() -> None:
    source = inspect.getsource(complete._global_candidates)
    assert "direct_full_matrix_SVD" not in source
    assert complete.GENERAL_SPECTRUM_ALGORITHM_VERSION == "general_complete_svd_v2"


def _synthetic_branch(root: float = 3.0) -> tuple[branch._Evaluator, branch.ContinuedBranch, branch.BranchOperationCounts]:
    def provider(value: float) -> np.ndarray:
        matrix = np.eye(6)
        matrix[:2, :2] = np.array([[value - root + 1.0, 1.0], [1.0, 1.0]])
        return matrix

    operations = branch.BranchOperationCounts()
    evaluator = branch._Evaluator(provider, operations, "local")
    left, right, _singular, sigma, ratio, nullity = evaluator.svd(root - 0.02)
    item = branch.ContinuedBranch(
        model="synthetic",
        branch_id="S001",
        parent_family="synthetic",
        beta_deg=0.0,
        predicted_Lambda=root - 0.02,
        Lambda=root - 0.02,
        right_vector=tuple(right),
        left_vector=tuple(left),
        sigma_1=sigma,
        sigma_ratio=ratio,
        self_MAC=1.0,
        nullity=nullity,
        cluster_id="",
        refinement_status=branch.SEED_REFINED_TO_NEW_ROOT,
        detection_source="synthetic_parent",
    )
    return evaluator, item, operations


def test_seed_is_refined_before_new_root_creation() -> None:
    evaluator, item, operations = _synthetic_branch()
    refined = branch._refine_seed(
        evaluator,
        item,
        2.97,
        branch.ContinuationSettings(lambda_max=10.0),
        operations,
    )
    assert refined.acceptance_status == branch.SEED_REFINED_TO_NEW_ROOT
    assert refined.Lambda != pytest.approx(refined.seed, abs=1.0e-6)
    assert refined.stationary and refined.full_matrix_svd_pass


def test_two_seeds_in_one_basin_create_one_new_root_record() -> None:
    evaluator, item, operations = _synthetic_branch()
    settings = branch.ContinuationSettings(lambda_max=10.0)
    refined = branch._refine_seed(evaluator, item, 2.97, settings, operations)
    duplicate = replace(refined, branch_id="S002", seed=3.03)
    operations.seed_windows_accepted = 2
    merged = branch._mark_duplicate_seed_refinements(
        [refined, duplicate], 2, settings, operations
    )
    assert [item.acceptance_status for item in merged] == [
        branch.SEED_REFINED_TO_NEW_ROOT,
        branch.SEED_REFINED_TO_EXISTING_ROOT,
    ]
    assert operations.seed_windows_accepted == 1
    assert operations.seed_windows_merged == 1


def test_seed_terminal_status_vocabulary_is_exact() -> None:
    assert {
        branch.SEED_REFINED_TO_EXISTING_ROOT,
        branch.SEED_REFINED_TO_NEW_ROOT,
        branch.SEED_REJECTED_NONSTATIONARY,
        branch.SEED_REJECTED_BOUNDARY_MINIMUM,
        branch.SEED_REJECTED_INSUFFICIENT_SINGULARITY,
        branch.SEED_REJECTED_INDEPENDENT_DISAGREEMENT,
    } == {
        "seed_refined_to_existing_root",
        "seed_refined_to_new_root",
        "seed_rejected_nonstationary",
        "seed_rejected_boundary_minimum",
        "seed_rejected_insufficient_singularity",
        "seed_rejected_independent_disagreement",
    }


def test_nonstationary_seed_is_rejected() -> None:
    def provider(value: float) -> np.ndarray:
        return np.diag([1.0 + 0.01 * value, 2.0, 3.0, 4.0, 5.0, 6.0])

    operations = branch.BranchOperationCounts()
    evaluator = branch._Evaluator(provider, operations, "local")
    item = branch.ContinuedBranch(
        "synthetic", "S001", "synthetic", 0.0, 3.0, 3.0,
        (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        1.0, 1.0, 1.0, 1, "", branch.SEED_REFINED_TO_NEW_ROOT, "synthetic",
    )
    refined = branch._refine_seed(
        evaluator, item, 3.0, branch.ContinuationSettings(lambda_max=10.0), operations
    )
    assert refined.acceptance_status != branch.SEED_REFINED_TO_NEW_ROOT


@pytest.mark.parametrize("model", complete.SUPPORTED_MODELS)
def test_isolated_small_angle_continuation_preserves_K10(
    small_angle_results: dict[str, branch.BranchContinuationResult], model: str
) -> None:
    result = small_angle_results[model]
    assert result.k10_guard_resolved
    assert len(result.values) >= 11
    assert result.operations.strict_fallback_runs == 0
    assert all(step.accepted for step in result.steps)


def test_timo_close_pair_uses_dimension_two_cluster(
    small_angle_results: dict[str, branch.BranchContinuationResult],
) -> None:
    clusters = [item for step in small_angle_results[complete.MODEL_TIMO].steps for item in step.clusters]
    assert any(item.expected_dimension == 2 and item.found_dimension == 2 for item in clusters)
    assert all(item.resolved for item in clusters)


def test_subspace_mac_is_permutation_invariant() -> None:
    e1 = (1.0, 0.0, 0.0)
    e2 = (0.0, 1.0, 0.0)
    assert branch._subspace_mac([e1, e2], [e2, e1]) == pytest.approx(1.0)


def test_two_stationary_close_roots_remain_distinct(
    small_angle_results: dict[str, branch.BranchContinuationResult],
) -> None:
    values = small_angle_results[complete.MODEL_TIMO].values
    assert all(right > left for left, right in zip(values, values[1:]))
    assert min(right - left for left, right in zip(values, values[1:])) > 1.0e-6


@pytest.mark.parametrize("model", complete.SUPPORTED_MODELS)
def test_first11_guard_and_full12_are_separate_flags(
    beta0_results: dict[str, branch.BranchContinuationResult], model: str
) -> None:
    result = beta0_results[model]
    assert isinstance(result.k10_guard_resolved, bool)
    assert isinstance(result.full12_resolved, bool)
    assert result.guard.passed


def test_compare_values_ignores_root12_for_K10_agreement() -> None:
    primary = tuple(float(index) for index in range(1, 13))
    verification = (*primary[:11], 99.0)
    rows, agreement = branch._compare_values(primary, verification, 1.0e-8)
    assert agreement
    assert rows[11]["within_tolerance"] is False


@pytest.mark.parametrize(
    "updates",
    [
        {"requested_roots": 11},
        {"candidate_roots": 11},
        {"verification_candidate_roots": 20},
        {"beta_min_step_deg": 1.0},
        {"mac_accept": 1.1},
        {"sigma_accept": -1.0},
    ],
)
def test_continuation_settings_validation_rejects_invalid_values(updates: dict[str, object]) -> None:
    settings = replace(branch.ContinuationSettings(), **updates)
    with pytest.raises(ValueError):
        settings.validate()


def test_cache_identity_includes_scope_and_settings(tmp_path: Path) -> None:
    settings = branch.ContinuationSettings()
    primary = branch.BranchContinuationCache(tmp_path, verification_scope="primary")
    force = branch.BranchContinuationCache(tmp_path, verification_scope="force_strict_verification")
    assert primary.path(complete.MODEL_EB, GEOMETRY_0, settings) != force.path(
        complete.MODEL_EB, GEOMETRY_0, settings
    )
    changed = replace(settings, seed_half_width=0.06)
    assert primary.path(complete.MODEL_EB, GEOMETRY_0, settings) != primary.path(
        complete.MODEL_EB, GEOMETRY_0, changed
    )
    identity = primary.identity(complete.MODEL_EB, GEOMETRY_SMALL, settings, "primary")
    assert identity["beta_path"] == {
        "start_beta_deg": 0.0,
        "target_beta_deg": 0.1,
        "policy": "adaptive_local_continuation",
    }
    assert identity["beta_targets"] == [0.0, 0.1]
    assert identity["strict_fallback_settings"]


def test_cache_identity_normalizes_equivalent_numeric_geometry(tmp_path: Path) -> None:
    settings = branch.ContinuationSettings()
    cache = branch.BranchContinuationCache(tmp_path)
    integer_geometry = complete.Geometry(0.02, 90, 0, 0)
    float_geometry = complete.Geometry(0.02, 90.0, 0.0, 0.0)
    assert cache.path(complete.MODEL_EB, integer_geometry, settings) == cache.path(
        complete.MODEL_EB, float_geometry, settings
    )


def test_cache_round_trip_reuses_result(tmp_path: Path) -> None:
    cache = branch.BranchContinuationCache(tmp_path)
    first = cache.resolve(complete.MODEL_EB, GEOMETRY_0, branch.ContinuationSettings())
    second = cache.resolve(complete.MODEL_EB, GEOMETRY_0, branch.ContinuationSettings())
    assert first.values == second.values
    assert second.cache_status == "hit"
    assert second.operations.cache_hits >= 1


def test_pilot_default_remains_legacy(tmp_path: Path) -> None:
    args = pilot.parse_args(["--output-dir", str(tmp_path / "pilot")])
    assert args.spectrum_method == pilot.SPECTRUM_METHOD_LEGACY


def test_pilot_accepts_branch_method(tmp_path: Path) -> None:
    args = pilot.parse_args(
        [
            "--output-dir", str(tmp_path / "pilot"),
            "--spectrum-method", branch.BRANCH_CONTINUATION_ALGORITHM_VERSION,
        ]
    )
    assert args.spectrum_method == branch.BRANCH_CONTINUATION_ALGORITHM_VERSION


@pytest.mark.parametrize(
    "argv",
    [
        ["--min-beta-step", "0"],
        ["--min-beta-step", "1", "--initial-beta-step", "0.5"],
        ["--n-spectrum-roots", "11"],
        ["--sigma-accept", "-1"],
        ["--mac-accept", "1.1"],
    ],
)
def test_gateway_cli_rejects_invalid_settings(tmp_path: Path, argv: list[str]) -> None:
    common = [
        "--output-dir", str(tmp_path / "gateway"),
        "--corrected-pilot-dir", str(tmp_path / "pilot"),
    ]
    with pytest.raises(ValueError):
        gateway.parse_args([*common, *argv])


def test_gateway_output_contract_lists_all_required_files() -> None:
    assert len(gateway.OUTPUT_NAMES) == 13
    assert "branch_k10_guard_summary.csv" in gateway.OUTPUT_NAMES
    assert "eb_timo_branch_continuation_gateway_report.md" in gateway.OUTPUT_NAMES


def test_gateway_rows_expose_required_audit_columns() -> None:
    result = branch.resolve_branch_spectrum(complete.MODEL_TIMO, GEOMETRY_SMALL)
    tables: dict[str, list[dict[str, object]]] = {
        name: []
        for name in (
            "parent", "steps", "seed", "cluster", "guard", "fallback",
            "summary", "verification", "close", "operations",
        )
    }
    gateway.collect_result("contract", result, tables)
    assert {
        "parent_family", "parent_family_index", "Lambda_beta0", "parent_sorted_index",
        "multiplicity_group", "block_sigma_1", "full_sigma_1", "quality_status",
    } <= set(tables["parent"][0])
    assert {
        "branch_or_cluster_id", "beta_previous", "beta_target", "attempted_step_deg",
        "accepted_step", "predicted_Lambda", "refined_Lambda", "MAC", "MAC_margin",
        "gap", "refinement_status", "step_reduction_reason",
    } <= set(next(row for row in tables["steps"] if row["record_type"] == "branch"))
    assert {
        "seed_source", "seed_Lambda", "window", "stationary_minimum", "refined_Lambda",
        "sigma_values", "merged_new_rejected_status", "related_accepted_root",
    } <= set(tables["seed"][0])
    assert {"found_root_count", "full_roots", "cluster_status"} <= set(tables["cluster"][0])
    assert {
        "roots_1_10_status", "root11_status", "root12_status", "first_disagreement_index",
        "unresolved_interval_below_guard", "seed_only_accepted_count",
        "cluster_ambiguity_count", "strict_fallback_used",
    } <= set(tables["summary"][0])


def test_future_manifest_has_required_schema_and_no_roots(tmp_path: Path) -> None:
    baseline_dir = tmp_path / "baseline"
    threshold_rows = []
    values = {
        2: (0.049091972167968748, 0.048649702148437494),
        3: (0.036973243652343751, 0.036640151367187498),
        4: (0.029675860839843748, 0.029408510742187498),
        5: (0.029675860839843748, 0.029408510742187498),
        6: (0.024798906738281248, 0.0245754931640625),
        7: (0.021305333496093745, 0.021113393554687495),
        8: (0.018677104980468751, 0.018508842773437499),
        9: (0.016627398925781252, 0.0164776025390625),
        10: (0.016627398925781252, 0.0164776025390625),
    }
    for prefix_n, (near, buffer) in values.items():
        threshold_rows.append(
            {
                "prefix_n": prefix_n,
                "threshold_status": "resolved",
                "epsilon_near_n": near,
                "epsilon_buffer_n": buffer,
            }
        )
    gateway.write_csv(baseline_dir / "baseline_critical_prefix_thresholds.csv", threshold_rows)
    path = gateway.write_step3_manifest(tmp_path / "manifest.csv", baseline_dir)
    rows = gateway.read_csv(path)
    assert len(rows) == 28
    required = {
        "case_id", "prefix_group", "epsilon_source", "epsilon", "beta_deg", "mu", "eta",
        "case_group", "adversarial_rationale", "required_spectrum_method",
        "required_K10_guard", "notes",
    }
    assert required == set(rows[0])
    assert all("root" not in key.lower() for key in rows[0])
    assert {row["epsilon_source"] for row in rows} == {"epsilon_near_n", "epsilon_buffer_n"}
    assert {row["prefix_group"] for row in rows} >= {"prefixes_4_5", "prefixes_9_10"}
    assert {float(row["epsilon"]) for row in rows} == {
        value for pair in values.values() for value in pair
    }
    assert any(float(row["beta_deg"]) == 90.0 for row in rows)
    assert any(float(row["mu"]) == 0.7 for row in rows)
    assert {float(row["eta"]) for row in rows} >= {-0.1, 0.1}


def test_operation_counter_categories_are_separate() -> None:
    fields = branch.BranchOperationCounts.__dataclass_fields__
    assert "parent_matrix_evaluations" in fields
    assert "local_matrix_evaluations" in fields
    assert "guard_matrix_evaluations" in fields
    assert "strict_fallback_matrix_evaluations" in fields
    assert "force_verification_matrix_evaluations" in fields


def test_algorithm_version_and_import_contract() -> None:
    assert branch.BRANCH_CONTINUATION_ALGORITHM_VERSION == "branch_informed_continuation_v1"
    assert callable(gateway.main)
    assert callable(pilot.main)
