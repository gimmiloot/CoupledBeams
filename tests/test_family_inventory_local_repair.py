from __future__ import annotations

from collections import defaultdict
import csv
import inspect
import math
from pathlib import Path
import re

import numpy as np

from scripts.analysis.thickness_mismatch.audits import audit_family_inventory_local_repair as audit
from scripts.lib import family_inventory_local_repair as repair
from scripts.lib import general_spectrum_completeness as complete


ROOT = Path(__file__).resolve().parents[1]
SOURCE_CSV = (
    ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "beta_sorted_spectrum_pilot"
    / "beta_sorted_spectrum_pilot.csv"
)
ORACLE_CSV = (
    ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "beta_sorted_spectrum_refined_pilot"
    / "refined_beta_sorted_spectrum.csv"
)
AUDIT_DIR = (
    ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "family_inventory_local_repair_audit"
)


def _thresholds() -> repair.DetectorThresholds:
    return repair.THRESHOLD_PROFILES["nominal"]


def _synthetic_family(*, missing: int = 0, beta: tuple[float, ...] = (0.0, 1.0, 2.0)) -> repair.FamilySpectrum:
    base = np.arange(1.0, 15.0)
    rows = [list((base + 0.02 * index)[:12]) for index in range(3)]
    if missing == 1:
        rows[1] = list((base + 0.02)[[0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12]])
    elif missing == 2:
        rows[1] = list((base + 0.02)[[0, 1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13]])
    return repair.FamilySpectrum(
        family_id="synthetic:Timoshenko",
        case_id="synthetic",
        theory="Timoshenko",
        epsilon_0=0.01,
        mu=0.0,
        eta=0.0,
        beta_values=beta,
        inventories=tuple(tuple(row) for row in rows),
    )


def _window(
    left: float = 0.52,
    right: float = 0.55,
    expected: int = 2,
) -> repair.RepairWindow:
    return repair.RepairWindow(
        event_id="test",
        case_id="test",
        theory="test",
        beta=1.0,
        rank_start=3,
        expected_missing_count=expected,
        lambda_left=left,
        lambda_right=right,
        source="test",
        lower_anchor=left - 0.1,
        upper_anchor=right + 0.1,
        predicted_roots=(),
        margin=0.01,
        beta_probe_required=False,
        status="window_inferred",
    )


def _settings() -> complete.SearchSettings:
    return complete.SearchSettings(
        requested_roots=12,
        candidate_roots=12,
        verification_candidate_roots=13,
        max_upper_growth_tries=1,
    )


def test_one_missing_root_produces_shift_one() -> None:
    events = repair.detect_family_inventory(_synthetic_family(missing=1), _thresholds())
    assert events
    assert {event.best_shift for event in events} == {1}


def test_two_missing_roots_produce_shift_two() -> None:
    events = repair.detect_family_inventory(_synthetic_family(missing=2), _thresholds())
    assert events
    assert {event.best_shift for event in events} == {2}


def test_smooth_crossing_does_not_trigger_missing_root() -> None:
    family = _synthetic_family()
    rows = [list(row) for row in family.inventories]
    rows[1][4] = 5.49
    rows[1][5] = 5.51
    crossing = repair.FamilySpectrum(**{**family.__dict__, "inventories": tuple(tuple(sorted(row)) for row in rows)})
    assert repair.detect_family_inventory(crossing, _thresholds()) == ()


def test_avoided_crossing_does_not_trigger_missing_root() -> None:
    family = _synthetic_family()
    rows = [list(row) for row in family.inventories]
    rows[1][4], rows[1][5] = 5.45, 5.55
    avoided = repair.FamilySpectrum(**{**family.__dict__, "inventories": tuple(tuple(row) for row in rows)})
    assert repair.detect_family_inventory(avoided, _thresholds()) == ()


def test_correct_close_pair_does_not_trigger_repair() -> None:
    family = _synthetic_family()
    rows = [list(row) for row in family.inventories]
    for row in rows:
        row[4], row[5] = 5.0, 5.0005
    close = repair.FamilySpectrum(**{**family.__dict__, "inventories": tuple(tuple(row) for row in rows)})
    assert repair.detect_family_inventory(close, _thresholds()) == ()


def test_solver_unresolved_interval_is_a_direct_trigger() -> None:
    triggered, reasons = repair.solver_diagnostic_trigger(
        unresolved_intervals=("interval_7",), required_guard=8
    )
    assert triggered
    assert reasons == ("unresolved_interval",)


def test_automatic_window_contains_source_inferred_missing_roots() -> None:
    family = _synthetic_family(missing=1)
    event = repair.detect_family_inventory(family, _thresholds())[0]
    window = repair.infer_repair_window(family, event, _thresholds())
    assert window.status == "window_inferred"
    assert window.lambda_left < 5.02 < window.lambda_right


def test_ambiguous_window_is_deferred() -> None:
    family = _synthetic_family(missing=1)
    event = repair.detect_family_inventory(family, _thresholds())[0]
    event = repair.DetectorEvent(
        **{
            **event.__dict__,
            "segment_beta_left": family.beta_values[0],
            "segment_beta_right": family.beta_values[-1],
        }
    )
    window = repair.infer_repair_window(family, event, _thresholds())
    assert window.status == "defer_window_ambiguous"


def test_exact_beta_lambda_stages_precede_any_beta_probe() -> None:
    provider = lambda value: np.diag([value - 0.531, value - 0.537, 1.0])
    result = repair.staged_local_search(provider, _window(), base_settings=_settings())
    assert result.status == "resolved_after_local_repair"
    assert result.stage_roots[0][0] == "L1"
    assert result.beta_probes == ()


def test_sign_change_pair_is_recovered() -> None:
    provider = lambda value: np.diag([value - 0.531, value - 0.537, 1.0])
    result = repair.staged_local_search(provider, _window(), base_settings=_settings())
    assert [round(entry.value, 6) for entry in result.entries] == [0.531, 0.537]


def test_sigma_only_root_is_recovered() -> None:
    center = 0.53743
    provider = lambda value: np.diag([(value - center) ** 2, 1.0, 1.0])
    result = repair.staged_local_search(provider, _window(expected=1), base_settings=_settings())
    assert result.status == "resolved_after_local_repair"
    assert abs(result.entries[0].value - center) < 1.0e-6
    assert any(row["accepted"] and not row["sign_change"] for row in result.candidate_rows)


def test_rejected_sigma_candidate_is_not_accepted() -> None:
    center = 0.53743
    provider = lambda value: np.asarray(
        [[1.0, 1.0, 0.0], [1.0, 1.0 + (value - center) ** 2 + 1.0e-3, 0.0], [0.0, 0.0, 1.0]]
    )
    result = repair.staged_local_search(provider, _window(expected=1), base_settings=_settings())
    assert result.status == "deferred_local_inventory_unstable"
    assert not result.entries
    assert any(not row["accepted"] for row in result.candidate_rows)


def test_confirmed_multiplicity_is_preserved() -> None:
    center = 0.535
    provider = lambda value: np.diag([value - center, value - center, 1.0])
    result = repair.staged_local_search(provider, _window(expected=2), base_settings=_settings())
    assert len(result.entries) == 2
    assert [entry.repeated_root_slot for entry in result.entries] == [1, 2]
    assert all(entry.multiplicity == 2 for entry in result.entries)


def test_beta_zero_distinct_block_pair_keeps_provenance() -> None:
    first, second = 0.531, 0.537
    provider = lambda value: np.diag([value - first, value - second, 1.0])
    blocks = {
        "axial_block": lambda value: np.diag([value - first, 1.0]),
        "bending_block": lambda value: np.diag([value - second, 1.0]),
    }
    result = repair.staged_local_search(
        provider, _window(), base_settings=_settings(), block_providers=blocks
    )
    assert result.block_classification == "resolved_distinct_pair"
    assert {entry.block_family for entry in result.entries} == {"axial_block", "bending_block"}


def test_merge_resorts_inventory_after_lower_recovery() -> None:
    original = tuple(float(index) for index in range(1, 13))
    recovered = (repair.LocalRootEntry(4.5, 1, 1, "full_matrix", 1, "test"),)
    window = _window(4.4, 4.6, expected=1)
    merged = repair.merge_inventory(original, recovered, window, root_dedup_tolerance=2.0e-4)
    assert merged[:7] == (1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 6.0)


def test_repair_scope_is_limited_to_required_guard() -> None:
    assert repair.repair_rank_is_required(6, 6)
    assert not repair.repair_rank_is_required(9, 6)


def test_problem_above_guard_does_not_block_n_true() -> None:
    eb = tuple(float(index) for index in range(1, 12))
    timo = tuple(float(index) for index in range(1, 12))
    n_true, failure, guard, _ = repair.compute_n_true(eb, timo)
    assert (n_true, failure, guard) == (10, None, 11)
    assert not repair.repair_rank_is_required(12, guard)


def _source_families() -> dict[tuple[str, str], repair.FamilySpectrum]:
    _rows, parameters, values, statuses = audit._source_payload()
    return audit._families(parameters, values, statuses)


def test_p2_negative_control_has_no_trigger() -> None:
    families = _source_families()
    assert repair.detect_family_inventory(families[("P2", "Euler-Bernoulli")], _thresholds()) == ()
    assert repair.detect_family_inventory(families[("P2", "Timoshenko")], _thresholds()) == ()


def test_manual_defect_nodes_are_detected_without_window_constants() -> None:
    events = [
        event
        for family in _source_families().values()
        for event in repair.detect_family_inventory(family, _thresholds())
    ]
    assert len(events) == 9
    assert len({(event.case_id, event.theory, event.segment_beta_left, event.segment_beta_right) for event in events}) == 7
    helper_source = inspect.getsource(repair)
    assert re.search(r"\bR[1-7]\b", helper_source) is None
    assert all(token not in helper_source for token in ("P1", "P2", "P3", "P4"))


def test_automated_repair_matches_post_run_oracle() -> None:
    before_after = list(csv.DictReader((AUDIT_DIR / "before_after.csv").open(newline="", encoding="utf-8")))
    assert before_after
    assert max(float(row["oracle_difference"]) for row in before_after) <= complete.DEFAULT_ROOT_MATCH_TOL


def test_completed_audit_has_all_gates_pass_and_distinct_p3_pair() -> None:
    gates = list(csv.DictReader((AUDIT_DIR / "gate_summary.csv").open(newline="", encoding="utf-8")))
    assert len(gates) == 6
    assert all(row["status"] == "PASS" for row in gates)
    report = (AUDIT_DIR / "report.md").read_text(encoding="utf-8")
    assert "P3 beta=0 classification: resolved_distinct_pair" in report
    assert "force/full strict execution: 0" in report


def test_readiness_24_n_true_is_unchanged() -> None:
    rows = list(csv.DictReader((audit.READINESS_DIR / "validation_24_cases.csv").open(newline="", encoding="utf-8")))
    assert len(rows) == 24
    assert all(row["full_N_true"] == row["auto_N_true"] for row in rows)


def test_force_strict_execution_is_structurally_zero() -> None:
    source = inspect.getsource(audit)
    assert "force_strict_verification(" not in source
    assert "resolve_general_spectrum(" not in source
    assert "resolve_matrix_spectrum(" not in source


def test_cache_round_trip_is_deterministic(tmp_path: Path) -> None:
    identity = {"cache": "test"}
    result = repair.LocalSearchResult(
        status="resolved_after_local_repair",
        repair_stage="L2",
        entries=(repair.LocalRootEntry(1.0, 1, 1, "full_matrix", 1, "test"),),
        candidate_rows=(),
        matrix_evaluations=10,
        stage_matrix_evaluations=(("L1", 5), ("L2", 5)),
        stage_roots=(("L1", (1.0,)), ("L2", (1.0,))),
        beta_probes=(),
        block_classification="not_applicable",
    )
    path = tmp_path / "cache" / "point.json"
    repair.save_cache(path, identity, result)
    first = path.read_bytes()
    loaded = repair.load_cache(path, identity)
    repair.save_cache(path, identity, repair.result_without_cache_flag(loaded))
    assert path.read_bytes() == first
    assert repair.load_cache(path, identity).cache_hit


def test_audit_paths_are_isolated_from_source_and_coarse_caches() -> None:
    assert audit.OUTPUT_DIR != audit.SOURCE_DIR
    assert audit.OUTPUT_DIR != audit.ORACLE_DIR
    assert audit.OUTPUT_DIR != audit.COARSE_DIR
    source = inspect.getsource(audit.main)
    assert "_write_csv(SOURCE_DIR" not in source
    assert "_write_csv(ORACLE_DIR" not in source
    assert "_write_csv(COARSE_DIR" not in source


def test_sparse_production_style_grid_is_supported() -> None:
    events = repair.detect_family_inventory(
        _synthetic_family(missing=1, beta=(0.0, 15.0, 30.0)), _thresholds()
    )
    assert events
    assert events[0].beta == 15.0
