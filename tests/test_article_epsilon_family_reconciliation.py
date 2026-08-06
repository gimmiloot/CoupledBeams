from __future__ import annotations

import csv
import inspect
import json
from pathlib import Path
import sys

import pytest

from scripts.analysis.thickness_mismatch.audits import (
    run_article_epsilon_upper_envelope_grid as runner,
)
from scripts.lib import article_epsilon_family_reconciliation as reconciliation


REPO_ROOT = Path(__file__).resolve().parents[1]
COARSE_DIR = (
    REPO_ROOT / "results" / "article_epsilon_upper_envelope" / "coarse_grid_v1"
)
OUTPUT_DIR = COARSE_DIR / reconciliation.OUTPUT_DIRECTORY_NAME


def _rows(name: str) -> list[dict[str, str]]:
    with (OUTPUT_DIR / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _metadata() -> dict[str, object]:
    return json.loads((OUTPUT_DIR / "run_metadata.json").read_text(encoding="utf-8"))


def test_scope_accepts_only_isotropic_circular_rods() -> None:
    assert reconciliation.validate_scope(reconciliation.SCIENTIFIC_SCOPE) == (
        "isotropic_circular_coupled_rods_eb_timoshenko"
    )


@pytest.mark.parametrize(
    "invalid",
    ("anisotropic_rods", "rectangular", "orthotropic", "monoclinic"),
)
def test_anisotropic_or_rectangular_scope_is_rejected(invalid: str) -> None:
    with pytest.raises(ValueError):
        reconciliation.validate_scope(invalid)


@pytest.mark.parametrize(
    "operation",
    ("matrix_evaluator", "point_solver", "local_repair"),
)
def test_reconcile_only_guard_rejects_scientific_operations(operation: str) -> None:
    with pytest.raises(reconciliation.ReconciliationIntegrityError):
        reconciliation.forbid_scientific_operation(operation)


@pytest.mark.parametrize(
    ("gate", "cache", "manifest", "mismatches"),
    (
        (False, "cache", "manifest", ()),
        (True, "wrong", "manifest", ()),
        (True, "cache", "wrong", ()),
        (True, "cache", "manifest", ("artifact",)),
    ),
)
def test_invalid_shadow_contract_blocks_promotion(
    gate: bool, cache: str, manifest: str, mismatches: tuple[str, ...]
) -> None:
    assert not reconciliation.shadow_input_contract_valid(
        required_gates_pass=gate,
        actual_cache_fingerprint=cache,
        expected_cache_fingerprint="cache",
        actual_manifest_sha256=manifest,
        expected_manifest_sha256="manifest",
        artifact_mismatches=mismatches,
    )


def test_zero_is_confirmed_but_nan_is_not() -> None:
    assert reconciliation.confirmed_n_true_value(0)
    assert reconciliation.confirmed_n_true_value("10")
    assert not reconciliation.confirmed_n_true_value("NaN")
    assert not reconciliation.confirmed_n_true_value("")


def test_reconciliation_cli_is_explicit_and_default_runner_is_unchanged() -> None:
    assert not runner.parse_args([]).reconcile_family_local_repair_shadow
    args = runner.parse_args(
        [
            "--reconcile-family-local-repair-shadow", "--reconcile-only",
            "--no-new-point-solves",
        ]
    )
    assert args.reconcile_only and args.no_new_point_solves
    with pytest.raises(SystemExit):
        runner.parse_args(["--reconcile-family-local-repair-shadow"])
    with pytest.raises(SystemExit):
        runner.parse_args(
            [
                "--reconcile-family-local-repair-shadow", "--reconcile-only",
                "--no-new-point-solves", "--prefix-until-failure",
            ]
        )


def test_reconciliation_dispatch_precedes_manifest_or_solver_work(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    called: dict[str, object] = {}

    def fake_reconcile(output_dir: Path, **kwargs: object) -> dict[str, object]:
        called.update(output_dir=output_dir, **kwargs)
        return {"root_calculations": 0}

    monkeypatch.setattr(reconciliation, "reconcile", fake_reconcile)
    monkeypatch.setattr(
        runner.workflow, "build_manifest",
        lambda **_kwargs: pytest.fail("manifest construction reached"),
    )
    result = runner.main(
        [
            "--output-dir", str(tmp_path / "coarse"),
            "--reconcile-family-local-repair-shadow", "--reconcile-only",
            "--no-new-point-solves",
        ]
    )
    assert result["root_calculations"] == 0
    assert called["promotion_policy"] == "verified-only"


def test_pending_family_status_is_not_an_expensive_strict_defer() -> None:
    assert reconciliation.pending_family_status(
        primary_complete=True, family_context_complete=False
    ) == "pending_family_inventory_check"
    assert reconciliation.pending_family_status(
        primary_complete=True, family_context_complete=True,
        poststage_outcome="resolved",
    ) == "resolved_after_local_repair"
    assert reconciliation.pending_family_status(
        primary_complete=True, family_context_complete=True,
        poststage_outcome="deferred",
    ) == "deferred_expensive_strict"


def test_family_poststage_is_independent_of_worker_completion_order() -> None:
    rows = [
        {
            "case_id": "b", "epsilon_0": 0.01, "mu": 0.9, "eta": -0.5,
            "beta": 15.0, "family_context_complete": True,
            "poststage_outcome": "resolved",
        },
        {
            "case_id": "a", "epsilon_0": 0.01, "mu": 0.9, "eta": -0.5,
            "beta": 0.0, "family_context_complete": True,
            "poststage_outcome": "deferred",
        },
    ]
    forward = reconciliation.deterministic_family_poststage(rows)
    reverse = reconciliation.deterministic_family_poststage(list(reversed(rows)))
    assert forward == reverse
    assert [row["case_id"] for row in forward] == ["a", "b"]


def test_artifact_gates_and_operation_counts_pass() -> None:
    gates = {row["gate"]: row["status"] for row in _rows("gate_summary.csv")}
    assert len(gates) == 9
    assert set(gates.values()) == {"PASS"}
    operations = {row["operation"]: int(row["count"]) for row in _rows("operation_counts.csv")}
    for name in (
        "point_solver_calls", "matrix_evaluator_calls", "local_repair_calls",
        "detector_calls", "force_strict_executed",
    ):
        assert operations[name] == 0
    assert operations["processed_not_attempted_points"] == 0


def test_verified_promotions_are_unique_and_provenance_complete() -> None:
    promotions = _rows("promotion_overlay.csv")
    assert len(promotions) == 7
    assert len({row["case_id"] for row in promotions}) == 7
    assert {row["promotion_status"] for row in promotions} == {
        "verified_shadow_promoted"
    }
    assert {row["N_true"] for row in promotions} == {"10"}
    assert all(row["source_cache_hash"] and row["repair_cache_hash"] for row in promotions)
    audit = {row["case_id"]: row for row in _rows("eligibility_audit.csv")}
    assert all(audit[row["case_id"]]["provenance_complete"] == "True" for row in promotions)


def test_deferred_cases_are_not_promoted_and_keep_nan() -> None:
    deferred = _rows("deferred_cases.csv")
    promoted = {row["case_id"] for row in _rows("promotion_overlay.csv")}
    assert len(deferred) == 18
    assert not promoted.intersection(row["case_id"] for row in deferred)
    assert all(row["execution_status"] == "deferred_expensive_strict" for row in deferred)
    assert all(row["N_true"] == "NaN" for row in deferred)
    assert all(row["force_strict_executed"] == "0" for row in deferred)


def test_original_resolved_results_and_reference_sets_are_preserved() -> None:
    comparisons = _rows("original_vs_reconciled_cases.csv")
    unchanged = [row for row in comparisons if row["reconciliation_status"] == "unchanged_original_resolved"]
    assert len(unchanged) == 684
    assert all(row["N_true_changed"] == "False" for row in unchanged)
    assert all(row["first_failed_mode_changed"] == "False" for row in unchanged)
    reference = _metadata()["reference_validation"]
    assert reference["S3_pass"] and reference["readiness_pass"] and reference["former_pass"]
    assert reference["readiness_count"] == 24
    assert reference["former_count"] == 7
    assert reference["S3"] == {
        "S3_12_N_true": 4,
        "S3_12_delta_f_5": 0.11739469908796035,
        "S3_14_N_true": 4,
        "S3_14_delta_f_5": 0.10050934855181458,
    }


def test_resume_plan_has_deterministic_disjoint_statuses() -> None:
    rows = _rows("resume_plan.csv")
    assert len(rows) == 1554
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["resume_status"]] = counts.get(row["resume_status"], 0) + 1
    assert counts == {
        "ready_not_attempted": 845,
        "skipped_deferred": 18,
        "skipped_existing_resolved": 684,
        "skipped_promoted_resolved": 7,
    }
    ready = {row["case_id"] for row in rows if row["resume_status"] == "ready_not_attempted"}
    assert ready == reconciliation.load_ready_resume_ids(COARSE_DIR)
    assert not ready.intersection(row["case_id"] for row in _rows("promotion_overlay.csv"))
    assert not ready.intersection(row["case_id"] for row in _rows("deferred_cases.csv"))


def test_reconciled_table_has_no_final_pending_or_fake_results() -> None:
    rows = _rows("reconciled_case_results.csv")
    assert len(rows) == 1554
    assert len({row["case_id"] for row in rows}) == 1554
    assert all(row["reconciled_execution_status"] != "pending_family_inventory_check" for row in rows)
    absent = [row for row in rows if row["reconciled_N_true"] == "NaN"]
    assert all(
        row["reconciled_execution_status"] in {"not_attempted", "deferred_expensive_strict"}
        for row in absent
    )


def test_source_inputs_are_immutable_and_missing_cache_count_is_explicit() -> None:
    metadata = _metadata()
    assert metadata["source_cache_fingerprint_before"] == metadata["source_cache_fingerprint_after"]
    assert metadata["manifest_sha256_before"] == metadata["manifest_sha256_after"]
    assert metadata["shadow_tree_fingerprint_before"] == metadata["shadow_tree_fingerprint_after"]
    assert metadata["immutable_hashes_before"] == metadata["immutable_hashes_after"]
    assert metadata["counts"]["point_cache_missing"] == 846
    assert metadata["counts"]["ready_not_attempted"] == 845


def test_repeat_hashes_are_idempotent_and_anisotropic_modules_are_absent() -> None:
    metadata = _metadata()
    assert metadata["all_csv_repeat_match"] is True
    assert metadata["gates"]["idempotence_gate"] == "PASS"
    assert not [
        name for name in sys.modules
        if "anisotropic_rods" in name or "yartsev" in name.lower()
    ]
    source = inspect.getsource(reconciliation)
    assert "scripts.analysis.anisotropic_rods" not in source

