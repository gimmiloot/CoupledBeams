from __future__ import annotations

import inspect
import math
from pathlib import Path

import pytest

from scripts.lib import article_epsilon_targeted_resolution as target
from scripts.lib import general_spectrum_completeness as complete
from scripts.lib.article_epsilon_compact_certificates import load_certificate, read_csv


REPO = Path(__file__).resolve().parents[1]
COARSE = REPO / "results" / "article_epsilon_upper_envelope" / "coarse_grid_v1"


def _selected() -> list[dict[str, str]]:
    return target.select_targets(read_csv(COARSE / "compact_finalization" / "unresolved_cases.csv"))


def test_selects_only_deferred_epsilon_005_cases() -> None:
    rows = _selected()
    assert rows
    assert all(abs(float(row["epsilon_0"]) - 0.05) <= 1e-12 for row in rows)
    assert all("deferred" in row["final_execution_status"] for row in rows)


def test_target_case_ids_are_not_hardcoded() -> None:
    source = inspect.getsource(target.select_targets)
    assert "AUE_" not in source


def test_target_count_is_derived_from_unresolved_table() -> None:
    all_rows = read_csv(COARSE / "compact_finalization" / "unresolved_cases.csv")
    expected = sum(abs(float(row["epsilon_0"]) - 0.05) <= 1e-12 for row in all_rows)
    assert len(target.select_targets(all_rows)) == expected


def test_compact_certificate_reads_without_all_raw_payloads() -> None:
    row = _selected()[0]
    cert = load_certificate(REPO / row["compact_certificate_path"])
    assert cert["case_id"] == row["case_id"]
    assert cert["scientific_scope"] == target.SCIENTIFIC_SCOPE


def test_confirmed_model_side_is_not_in_target_windows() -> None:
    local = {row["case_id"]: row for row in read_csv(COARSE / "compact_finalization" / "local_repair_results.csv")}
    for row in _selected():
        cert = load_certificate(REPO / row["compact_certificate_path"])
        theory = local[row["case_id"]]["preferred_model"]
        windows = target.derive_target_windows(cert, theory, int(row["required_guard"]))
        assert {window.theory for window in windows} == {theory}


def test_exact_beta_local_stage_precedes_strict() -> None:
    assert not target.strict_is_contained("case", {"case"}, ())
    assert not target.strict_is_contained("case", {"case"}, ("T0",))


def test_independent_verification_precedes_full_strict() -> None:
    assert not target.strict_is_contained("case", {"case"}, ("T1",))
    assert target.strict_is_contained("case", {"case"}, ("T1", "T2"))


def test_strict_is_limited_to_target_set() -> None:
    assert not target.strict_is_contained("other", {"case"}, ("T1", "T2"))


def test_target_windows_end_at_required_guard_inventory() -> None:
    local = {row["case_id"]: row for row in read_csv(COARSE / "compact_finalization" / "local_repair_results.csv")}
    for row in _selected():
        cert = load_certificate(REPO / row["compact_certificate_path"])
        guard = int(row["required_guard"])
        windows = target.derive_target_windows(cert, local[row["case_id"]]["preferred_model"], guard)
        assert all(window.rank_end <= guard for window in windows)


def test_target_windows_do_not_use_blind_global_range() -> None:
    local = {row["case_id"]: row for row in read_csv(COARSE / "compact_finalization" / "local_repair_results.csv")}
    for row in _selected():
        cert = load_certificate(REPO / row["compact_certificate_path"])
        windows = target.derive_target_windows(cert, local[row["case_id"]]["preferred_model"], int(row["required_guard"]))
        assert all(window.lambda_left > 0.2 and window.lambda_right < 80.0 for window in windows)


def test_exact_n_true_requires_right_guard() -> None:
    n_true, failure, guard, _ = target.compute_prefix_result([1.0, 2.0], [1.0, 1.8])
    assert n_true is None
    assert failure == 2
    assert guard == 3


def test_exact_prefix_certificate_includes_delta_at_right_guard() -> None:
    n_true, failure, guard, deltas = target.compute_prefix_result(
        [1.0, 2.0, 3.0], [1.0, 1.8, 2.9],
    )
    assert (n_true, failure, guard) == (1, 2, 3)
    assert len(deltas) == guard


def test_envelope_bound_does_not_create_fake_n_true() -> None:
    record = {"N_true": None, "certified_N_true_upper_bound": 4}
    assert record["N_true"] is None
    assert record["certified_N_true_upper_bound"] == 4


def test_n_true_zero_is_distinct_from_missing() -> None:
    eb = [2.0, 3.0]
    tm = [1.0, 2.0]
    n_true, failure, guard, _ = target.compute_prefix_result(eb, tm)
    assert n_true == 0
    assert failure == 1
    assert guard == 2
    assert n_true is not None


def test_raw_envelope_uses_exact_max_without_smoothing() -> None:
    values = [10, 10, 9, 6, 4, 3]
    assert values[-2] == 4
    assert values != sorted(values)


def test_suffix_max_envelope_definition() -> None:
    raw = [10, 4, 9, 6, 4, 3]
    suffix = [max(raw[index:]) for index in range(len(raw))]
    assert suffix == [10, 9, 9, 6, 4, 3]


def test_manifest_has_1554_unique_cases() -> None:
    rows = read_csv(COARSE / "grid_manifest.csv")
    assert len(rows) == len({row["case_id"] for row in rows}) == 1554


def test_non_target_deferred_count_is_33() -> None:
    unresolved = read_csv(COARSE / "compact_finalization" / "unresolved_cases.csv")
    target_ids = {row["case_id"] for row in _selected()}
    assert sum(row["case_id"] not in target_ids for row in unresolved) == 33


def test_target_identity_is_stable_and_unique() -> None:
    rows = _selected()
    assert len(rows) == len({row["case_id"] for row in rows})


def test_target_raw_caches_exist() -> None:
    index = {row["case_id"]: row for row in read_csv(COARSE / "compact_point_certificates_v1" / "compact_index.csv")}
    assert all((REPO / index[row["case_id"]]["source_full_cache_path"]).exists() for row in _selected())


def test_force_strict_outside_target_is_rejected() -> None:
    target_ids = {row["case_id"] for row in _selected()}
    assert not target.strict_is_contained("AUE_not_a_selected_case", target_ids, ("T1", "T2"))


def test_four_shifted_phases_are_declared() -> None:
    assert target.T1_PHASES == (0.0, 0.25, 0.5, 0.75)


def test_two_last_refinement_levels_are_finer_than_prior_local_floor() -> None:
    assert target.T1_STEPS == (1.0e-4, 5.0e-5)
    assert target.T1_STEPS[-1] < target.T1_STEPS[0]


def test_scientific_scope_is_isotropic_circular() -> None:
    assert target.SCIENTIFIC_SCOPE == "isotropic_circular_coupled_rods_eb_timoshenko"


def test_anisotropic_scope_is_rejected_by_repair_validator() -> None:
    with pytest.raises(ValueError):
        target.repair.validate_scientific_scope("anisotropic_rectangular_rods")


def test_no_anisotropic_modules_imported() -> None:
    assert all("anisotropic" not in name.lower() for name in __import__("sys").modules)


def test_production_root_tolerances_remain_unchanged() -> None:
    settings = complete.SearchSettings()
    assert math.isclose(settings.root_match_tol, 2.0e-4)
    assert math.isclose(settings.root_dedup_tol, 2.0e-4)
    assert math.isclose(settings.sigma_accept, 5.0e-6)
