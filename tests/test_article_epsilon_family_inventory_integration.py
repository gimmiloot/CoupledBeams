from __future__ import annotations

import inspect
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest

from scripts.analysis.thickness_mismatch.audits import (
    run_article_epsilon_upper_envelope_grid as runner,
)
from scripts.lib import article_epsilon_family_inventory_integration as integration
from scripts.lib import family_inventory_local_repair as repair
from scripts.lib import general_spectrum_completeness as complete


def test_scope_is_explicitly_isotropic_circular_only() -> None:
    assert integration.validate_scope(integration.SCIENTIFIC_SCOPE) == (
        "isotropic_circular_coupled_rods_eb_timoshenko"
    )
    for invalid in (
        "anisotropic", "rectangular", "orthotropic", "monoclinic"
    ):
        with pytest.raises(ValueError):
            integration.validate_scope(invalid)


def test_family_key_is_parameter_tuple_without_beta() -> None:
    left = {"epsilon_0": 0.01, "mu": 0.9, "eta": -0.25, "beta": 0.0}
    right = {**left, "beta": 90.0}
    assert integration.family_key(left, "Timoshenko") == integration.family_key(
        right, "Timoshenko"
    )
    assert integration.family_id(integration.family_key(left, "Timoshenko")) == (
        integration.family_id(integration.family_key(right, "Timoshenko"))
    )


def test_sparse_and_low_angle_beta_nodes_are_supported_deterministically() -> None:
    base = np.arange(1.0, 13.0)
    sparse = repair.FamilySpectrum(
        "sparse", "sparse", "Timoshenko", 0.01, 0.0, 0.0,
        (0.0, 15.0, 90.0),
        (tuple(base), tuple(base + 0.02), tuple(base + 0.04)),
    )
    low_angle = repair.FamilySpectrum(
        "low", "low", "Timoshenko", 0.01, 0.0, 0.0,
        (0.0, 5.0, 10.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0),
        tuple(tuple(base + beta * 1.0e-3) for beta in (0, 5, 10, 15, 30, 45, 60, 75, 90)),
    )
    sparse.validate()
    low_angle.validate()
    assert sparse.beta_values == tuple(sorted(sparse.beta_values))
    assert low_angle.beta_values == tuple(sorted(low_angle.beta_values))


def test_runner_reuses_detector_helper_instead_of_copying_formula() -> None:
    source = inspect.getsource(integration._detect)
    assert "repair.detect_family_inventory" in source
    assert "normalized_mismatch" not in source


def test_runner_default_is_off_and_shadow_is_postprocess_only() -> None:
    assert runner.parse_args([]).family_inventory_policy == "off"
    shadow = runner.parse_args(
        ["--postprocess-only", "--family-inventory-policy", "shadow"]
    )
    assert shadow.family_inventory_policy == "shadow"
    with pytest.raises(SystemExit):
        runner.parse_args(["--family-inventory-policy", "shadow"])


def test_local_repair_policy_requires_main_pass_and_defer() -> None:
    with pytest.raises(SystemExit):
        runner.parse_args(["--family-inventory-policy", "local-repair"])
    args = runner.parse_args(
        [
            "--prefix-until-failure", "--main-pass-only",
            "--family-inventory-policy", "local-repair",
            "--defer-expensive-strict",
        ]
    )
    assert args.main_pass_only and args.defer_expensive_strict


def test_shadow_dispatch_does_not_call_ordinary_postprocess(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    called: list[Path] = []

    def fake_shadow(path: Path) -> dict[str, object]:
        called.append(path)
        return {"root_calculations": 0, "gates": {}}

    monkeypatch.setattr(integration, "run_shadow", fake_shadow)
    monkeypatch.setattr(
        runner, "postprocess_only",
        lambda _path: pytest.fail("ordinary postprocess must not run in shadow mode"),
    )
    result = runner.main(
        [
            "--postprocess-only", "--family-inventory-policy", "shadow",
            "--output-dir", str(tmp_path),
        ]
    )
    assert called == [tmp_path]
    assert result["root_calculations"] == 0


def test_repair_identity_includes_scope_versions_tolerances_and_window() -> None:
    window = repair.RepairWindow(
        "event", "case", "Timoshenko", 15.0, 5, 1, 4.0, 5.0,
        "test", 3.9, 5.1, (4.5,), 0.1, False, "window_inferred",
    )
    identity = repair.cache_identity(
        family_id="family", theory="Timoshenko", beta=15.0,
        source_hash="abc", threshold_profile=repair.THRESHOLD_PROFILES["nominal"],
        window=window,
        base_settings=complete.SearchSettings(requested_roots=12, candidate_roots=12),
    )
    assert identity["scientific_scope"] == integration.SCIENTIFIC_SCOPE
    assert identity["detector_version"] == repair.DETECTOR_VERSION
    assert identity["repair_algorithm_version"] == repair.REPAIR_ALGORITHM_VERSION
    assert identity["tolerance_hash"]


def test_required_guard_limits_repair_and_upper_issue_is_independent() -> None:
    assert repair.repair_rank_is_required(6, 6)
    assert not repair.repair_rank_is_required(7, 6)


def test_no_prohibited_research_modules_are_imported_by_integration() -> None:
    source = inspect.getsource(integration)
    assert "scripts.analysis.anisotropic_rods" not in source
    assert "src.my_project.fem" not in source
    loaded = set(sys.modules)
    assert not any(name.startswith("scripts.analysis.anisotropic_rods") for name in loaded)


def test_standalone_reference_gate_is_available_without_matrix_solve() -> None:
    gate_rows = integration._read_csv(
        integration.standalone.OUTPUT_DIR / "gate_summary.csv"
    )
    assert gate_rows
    assert all(row["status"] == "PASS" for row in gate_rows)


def test_saved_shadow_overlay_preserves_resolved_and_never_assigns_deferred_n() -> None:
    output = (
        integration.standalone.COARSE_DIR / integration.SHADOW_DIRECTORY_NAME
    )
    rows = integration._read_csv(output / "original_vs_shadow_cases.csv")
    assert rows
    for row in rows:
        assert row["scientific_scope"] == integration.SCIENTIFIC_SCOPE
        if row["original_N_true"]:
            assert row["shadow_N_true"] == row["original_N_true"]
        if row["shadow_execution_status"] == "deferred_expensive_strict":
            assert row["shadow_N_true"] == ""


def test_saved_shadow_matches_standalone_and_has_zero_solve_repeat() -> None:
    output = integration.standalone.COARSE_DIR / integration.SHADOW_DIRECTORY_NAME
    rows = integration._read_csv(output / "original_vs_shadow_cases.csv")
    passing, checked, differences = integration._standalone_gate(rows)
    metadata = json.loads((output / "run_metadata.json").read_text(encoding="utf-8"))
    assert passing and checked == 17 and differences == []
    assert metadata["current_invocation"]["root_calculations"] == 0
    assert metadata["operation_counts"]["force_strict_executed"] == 0


def test_saved_shadow_csv_hashes_and_gates_are_deterministic() -> None:
    output = integration.standalone.COARSE_DIR / integration.SHADOW_DIRECTORY_NAME
    metadata = json.loads((output / "run_metadata.json").read_text(encoding="utf-8"))
    for name, expected in metadata["deterministic_csv_hashes"].items():
        actual = hashlib.sha256((output / name).read_bytes()).hexdigest()
        assert actual == expected
    gates = {row["gate"]: row["status"] for row in integration._read_csv(output / "gate_summary.csv")}
    assert gates["production_resume_readiness_gate"] == "PASS"


def test_shadow_repair_cache_paths_are_unique_per_case_and_theory() -> None:
    output = integration.standalone.COARSE_DIR / integration.SHADOW_DIRECTORY_NAME
    paths = sorted((output / "cache" / "coarse_cases").glob("*.json"))
    assert paths
    assert len(paths) == len({path.name for path in paths})
