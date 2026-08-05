from __future__ import annotations

import csv
import json
from pathlib import Path
from unittest import mock

import pytest

from scripts.analysis.thickness_mismatch.audits import (
    audit_article_epsilon_solver_readiness_v2 as readiness,
)


SEVEN_EXPECTED = {
    "AUE_cc9a93b84d6bd27b0e06": (10, None),
    "AUE_54ec94c09aa621016e48": (2, 3),
    "AUE_dae89daaca36b3b79c16": (1, 2),
    "AUE_a4d07dc121ad182e53e4": (10, None),
    "AUE_e5e288ed487815468c25": (9, 10),
    "AUE_3d054f56e6eaefbd1952": (10, None),
    "AUE_5d57110aceecdc1ef72e": (1, 2),
}


def _payload(n_true: int, first_failure: int | None, *, resolved: bool = True) -> dict[str, object]:
    return {
        "execution_status": "resolved_full_K10" if resolved else "unresolved_interval",
        "N_true": n_true if resolved else None,
        "first_failed_mode": first_failure if resolved else None,
        "prefix_guard_status": "prefix_guard_resolved" if resolved else "prefix_guard_unresolved",
    }


def test_paired_auto_equivalence_contract_covers_all_seven_former_blockers() -> None:
    assert tuple(SEVEN_EXPECTED) == readiness.SEVEN_IDS
    for n_true, first_failure in SEVEN_EXPECTED.values():
        full = _payload(n_true, first_failure)
        auto = _payload(n_true, first_failure)
        checks = readiness.validation_case_checks(full, auto, root_agreement=True)
        assert checks == {
            "execution_status_agreement": True,
            "N_true_agreement": True,
            "first_failed_mode_agreement": True,
            "prefix_guard_agreement": True,
            "case_pass": True,
        }


def test_unresolved_or_root_disagreement_cannot_pass_case_gate() -> None:
    full = _payload(4, 5)
    unresolved = _payload(4, 5, resolved=False)
    assert not readiness.validation_case_checks(full, unresolved, root_agreement=True)["case_pass"]
    assert not readiness.validation_case_checks(full, full, root_agreement=False)["case_pass"]


def test_staged_gates_are_independent_and_full_readiness_requires_all_four() -> None:
    blocked = readiness.calculate_gate_statuses(
        oracle_pass=True,
        cutoff_rank_pass=True,
        old_regime_pass=True,
        seven_references_pass=False,
        seven_optimization_pass=True,
        validation_pass=True,
        full_resolved=False,
    )
    assert blocked["optimization_equivalence_gate"] == "PASS"
    assert blocked["full_grid_solver_readiness_gate"] == "BLOCKED_BY_UNRESOLVED_REFERENCE"

    passed = readiness.calculate_gate_statuses(
        oracle_pass=True,
        cutoff_rank_pass=True,
        old_regime_pass=True,
        seven_references_pass=True,
        seven_optimization_pass=True,
        validation_pass=True,
        full_resolved=True,
    )
    assert passed == {
        "basis_formula_gate": "PASS",
        "old_regime_regression_gate": "PASS",
        "seven_case_reference_gate": "PASS",
        "optimization_recheck_gate": "PASS",
        "optimization_equivalence_gate": "PASS",
        "full_grid_solver_readiness_gate": "PASS",
    }


def test_full_24_point_gate_fails_on_any_single_case_disagreement() -> None:
    statuses = [True] * 24
    statuses[17] = False
    gates = readiness.calculate_gate_statuses(
        oracle_pass=True,
        cutoff_rank_pass=True,
        old_regime_pass=True,
        seven_references_pass=True,
        seven_optimization_pass=True,
        validation_pass=all(statuses),
        full_resolved=True,
    )
    assert gates["optimization_equivalence_gate"] == "FAIL_DISAGREEMENT"
    assert gates["full_grid_solver_readiness_gate"] == "PASS"


def test_postprocess_only_reuses_high_precision_csv_without_root_solve(tmp_path: Path) -> None:
    path = tmp_path / "high_precision_checks.csv"
    fields = (
        "check_id",
        "validation_id",
        "double_precision_root",
        "double_precision_scaled_sigma_min",
        "status",
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerow(
            {
                "check_id": "existing",
                "validation_id": "case",
                "double_precision_root": "3.0",
                "double_precision_scaled_sigma_min": "1e-12",
                "status": "PASS",
            }
        )
    with mock.patch.object(readiness, "high_precision_rows", side_effect=AssertionError("root solve")):
        rows, solves, provenance = readiness.load_or_compute_high_precision_rows(
            postprocess_only=True,
            path=path,
            cases=(),
            results={},
        )
    assert len(rows) == 1
    assert solves == 0
    assert provenance.startswith("reused_existing")


def test_saved_24_point_and_seven_case_scientific_regressions() -> None:
    output = readiness.OUTPUT_DIR
    required = (
        output / "validation_24_cases.csv",
        output / "benchmark_cases.csv",
        output / "run_metadata.json",
        output / "unresolved_cases.csv",
        output / "old_regime_comparison.csv",
        output / "accuracy_gate.csv",
    )
    if not all(path.exists() for path in required):
        pytest.skip("versioned solver-readiness artifacts are not available")

    with required[0].open(encoding="utf-8", newline="") as handle:
        validation = list(csv.DictReader(handle))
    assert len(validation) == 24
    assert all(row["full_execution_status"] == "resolved_full_K10" for row in validation)
    assert all(row["optimization_case_status"] == "PASS" for row in validation)
    assert all(row["root_agreement"] == "True" for row in validation)
    assert all(row["prefix_guard_agreement"] == "True" for row in validation)

    by_id = {row["validation_id"]: row for row in validation}
    for case_id, (n_true, first_failure) in SEVEN_EXPECTED.items():
        row = by_id[case_id]
        assert int(row["full_N_true"]) == int(row["auto_N_true"]) == n_true
        expected_failure = "" if first_failure is None else str(first_failure)
        assert row["first_failed_mode_full"] == row["first_failed_mode_auto"] == expected_failure

    with required[1].open(encoding="utf-8", newline="") as handle:
        cases = list(csv.DictReader(handle))
    seven_full = [
        row
        for row in cases
        if row["run_mode"] == readiness.FULL_MODE and row["case_id"] in SEVEN_EXPECTED
    ]
    assert len(seven_full) == 7
    assert all(row["EB_root_count"] == row["Timo_root_count"] == "11" for row in seven_full)
    assert all(row["strict_status"] == "full_strict_pass" for row in seven_full)

    metadata = json.loads(required[2].read_text(encoding="utf-8"))
    assert metadata["resolved_full_references"] == 24
    assert metadata["optimization_pass_cases"] == 24
    assert metadata["S3_12_delta_f_5"] == pytest.approx(0.11739469908977435, abs=5.0e-10)
    assert metadata["S3_14_delta_f_5"] == pytest.approx(0.10050934855191349, abs=5.0e-10)
    assert metadata["full_grid_solver_readiness_gate"] == "PASS"
    assert metadata["historical_old_regime_comparison_count"] > 0
    assert metadata["historical_old_regime_max_root_difference"] <= 2.0e-4

    with required[3].open(encoding="utf-8", newline="") as handle:
        assert list(csv.DictReader(handle)) == []
    with required[4].open(encoding="utf-8", newline="") as handle:
        assert all(row["status"] == "PASS" for row in csv.DictReader(handle))
    with required[5].open(encoding="utf-8", newline="") as handle:
        gates = list(csv.DictReader(handle))
    assert len(gates) == 24
    assert all(row["optimization_equivalence_gate"] == "PASS" for row in gates)
    assert all(row["full_grid_solver_readiness_gate"] == "PASS" for row in gates)
