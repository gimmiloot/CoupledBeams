from __future__ import annotations

import inspect
from pathlib import Path

import pytest
import numpy as np

from scripts.analysis.thickness_mismatch.maps import (
    plot_eb_vs_timoshenko_lambda_beta_cases as beta_map,
)
from scripts.lib import general_spectrum_completeness as complete


def _model_payload(beta: int, *, shift: float = 0.0) -> dict[str, object]:
    values = [float(index + 1) + 0.01 * beta + shift for index in range(12)]
    return {
        "values": values,
        "root_sources": ["sign_change" for _ in values],
        "root_count": 12,
        "candidate_root_count": 12,
        "point_status": "complete_root_inventory",
        "lambda_min": 0.2,
        "lambda_max": 22.0,
        "scan_step": 0.01,
        "matrix_evaluations": 100,
        "unresolved_intervals": [],
        "resolved_cluster_count": 0,
        "independent_verification_runs": 0,
    }


def _synthetic_results() -> dict[str, dict[str, object]]:
    output: dict[str, dict[str, object]] = {}
    for label, epsilon, mu, eta in beta_map.SORTED_PILOT_CASES:
        points = []
        for beta in range(91):
            points.append(
                {
                    "case_id": label,
                    "epsilon_0": epsilon,
                    "mu": mu,
                    "eta": eta,
                    "beta_deg": beta,
                    "models": {
                        beta_map.MODEL_EB: _model_payload(beta),
                        beta_map.MODEL_TIMO: _model_payload(beta, shift=-0.02),
                    },
                }
            )
        output[label] = {"case": (label, epsilon, mu, eta), "points": points, "cache_hits": 0}
    return output


def test_sorted_pilot_preset_is_explicit_and_mutually_exclusive() -> None:
    args = beta_map.parse_args(["--beta-sorted-spectrum-pilot"])
    assert args.beta_sorted_spectrum_pilot is True
    assert args.beta_branch_pilot is False
    with pytest.raises(ValueError, match="mutually exclusive"):
        beta_map.parse_args(["--beta-sorted-spectrum-pilot", "--beta-branch-pilot"])


def test_primary_only_entry_point_passes_no_seeds_or_verification(monkeypatch: pytest.MonkeyPatch) -> None:
    settings = beta_map.sorted_pilot_search_settings()
    geometry = complete.Geometry(epsilon_0=0.01, beta_deg=12.0, mu=0.0, eta=0.0)
    provider = lambda value: value  # noqa: E731
    captured: dict[str, object] = {}

    monkeypatch.setattr(complete, "model_matrix_provider", lambda model, point: provider)

    def fake_run_configuration(matrix_provider, active, **kwargs):
        captured.update(kwargs)
        return "primary-result"

    monkeypatch.setattr(complete, "_run_configuration", fake_run_configuration)
    result = complete.resolve_primary_spectrum(complete.MODEL_EB, geometry, settings=settings)

    assert result == "primary-result"
    assert captured["configuration"] == "primary"
    assert captured["candidate_target"] == 12
    assert captured["phases"] == (0.0, settings.shifted_grid_phase)
    assert captured["seed_roots"] == ()
    assert settings.max_upper_growth_tries == 1


def test_sorted_pilot_rows_and_diagnostics_are_deterministic() -> None:
    results = _synthetic_results()
    rows = beta_map._sorted_pilot_csv_rows(results)
    diagnostics = beta_map._sorted_pilot_step_rows(results)

    assert len(rows) == 4 * 91 * 12
    assert len(diagnostics) == 4 * 2 * 90 * 12
    assert [(rows[index]["case_id"], rows[index]["beta_deg"], rows[index]["sorted_mode_index"]) for index in range(3)] == [
        ("P1", 0, 1),
        ("P1", 0, 2),
        ("P1", 0, 3),
    ]
    assert beta_map._sorted_pilot_suspect_rows(diagnostics) == []


def test_incomplete_inventory_is_nan_and_enters_suspect_list() -> None:
    results = _synthetic_results()
    p2_point = results["P2"]["points"][20]
    for model in beta_map.MODELS:
        payload = p2_point["models"][model]
        payload["values"] = payload["values"][:11]
        payload["root_sources"] = payload["root_sources"][:11]
        payload["root_count"] = 11
        payload["point_status"] = "incomplete_root_inventory"

    rows = beta_map._sorted_pilot_csv_rows(results)
    missing = [
        row
        for row in rows
        if row["case_id"] == "P2" and row["beta_deg"] == 20 and row["sorted_mode_index"] == 12
    ]
    assert len(missing) == 1
    assert str(missing[0]["lambda_eb"]) == "nan"
    assert str(missing[0]["lambda_timo"]) == "nan"

    suspect = beta_map._sorted_pilot_suspect_rows(beta_map._sorted_pilot_step_rows(results))
    assert {(row["case_id"], row["theory"], row["beta_left"], row["beta_right"]) for row in suspect} == {
        ("P2", "EB", 19, 21),
        ("P2", "Timoshenko", 19, 21),
    }
    assert all("incomplete_root_inventory" in str(row["reason"]) for row in suspect)


def test_pointwise_unresolved_interval_enters_suspect_list_without_recalculation() -> None:
    results = _synthetic_results()
    payload = results["P2"]["points"][30]["models"][beta_map.MODEL_TIMO]
    payload["point_status"] = "diagnostic_unresolved_interval"
    payload["unresolved_intervals"] = ["11.8:12.2:primary_shifted:global_sigma_valley"]

    suspect = beta_map._sorted_pilot_suspect_rows(beta_map._sorted_pilot_step_rows(results))
    matches = [row for row in suspect if row["case_id"] == "P2" and row["theory"] == "Timoshenko"]
    assert len(matches) == 1
    assert (matches[0]["beta_left"], matches[0]["beta_right"]) == (29, 31)
    assert "pointwise_unresolved_interval" in str(matches[0]["reason"])
    assert "12" in str(matches[0]["affected_sorted_modes"])


def test_sorted_model_search_has_no_tracking_mac_shape_or_strict_calls() -> None:
    source = inspect.getsource(beta_map._sorted_pilot_model_search)
    for forbidden in (
        "TRACK.",
        "analytic_shape_vector_eta",
        "MAC",
        "force_strict_verification",
        "resolve_general_spectrum",
        "resolve_matrix_spectrum",
    ):
        assert forbidden not in source


def test_output_validation_rejects_pdf_and_accepts_exact_png_set(tmp_path: Path) -> None:
    rows = beta_map._sorted_pilot_csv_rows(_synthetic_results())
    pngs = []
    for label, *_ in beta_map.SORTED_PILOT_CASES:
        path = tmp_path / f"{label}_sorted_lambda_beta_comparison.png"
        path.write_bytes(b"png")
        pngs.append(path)
    validation = beta_map.validate_sorted_pilot_artifacts(rows, pngs, tmp_path)
    assert validation["csv_rows"] == 4368
    assert validation["force_strict_verification_calls"] == 0

    (tmp_path / "forbidden.pdf").write_bytes(b"pdf")
    with pytest.raises(AssertionError, match="forbidden"):
        beta_map.validate_sorted_pilot_artifacts(rows, pngs, tmp_path)


def _local_settings(start: float = 0.4, end: float = 0.7, step: float = 0.001) -> complete.SearchSettings:
    return complete.SearchSettings(
        requested_roots=12,
        candidate_roots=12,
        verification_candidate_roots=13,
        lambda_min=start,
        lambda_max=end,
        scan_step=step,
        max_upper_growth_tries=1,
    )


def test_refined_preset_is_explicit_and_mutually_exclusive() -> None:
    args = beta_map.parse_args(["--beta-sorted-spectrum-refined-pilot"])
    assert args.beta_sorted_spectrum_refined_pilot is True
    assert args.beta_sorted_spectrum_pilot is False
    with pytest.raises(ValueError, match="mutually exclusive"):
        beta_map.parse_args(
            ["--beta-sorted-spectrum-refined-pilot", "--beta-sorted-spectrum-pilot"]
        )


def test_local_search_finds_two_roots_inside_one_old_coarse_interval() -> None:
    provider = lambda value: np.diag([value - 0.531, value - 0.537, 1.0])
    result = beta_map._dense_local_candidate_search(
        provider,
        _local_settings(0.52, 0.55, 0.001),
        candidate_source_prefix="test_local_dense",
        block_family="full_6x6",
    )
    values = sorted({round(float(entry["value"]), 6) for entry in result["entries"]})
    assert values == [0.531, 0.537]


def test_sigma_minimum_without_sign_change_can_be_accepted() -> None:
    center = 0.53743
    provider = lambda value: np.diag([(value - center) ** 2, 1.0, 1.0])
    result = beta_map._dense_local_candidate_search(
        provider,
        _local_settings(0.52, 0.55, 0.001),
        candidate_source_prefix="test_local_dense",
        block_family="full_6x6",
    )
    sigma_only = [
        row
        for row in result["candidate_rows"]
        if row["accepted"] and "sigma" in str(row["candidate_source"]) and not row["sign_change"]
    ]
    assert sigma_only
    assert abs(float(result["entries"][0]["value"]) - center) < 1.0e-6


def test_rejected_sigma_minimum_is_not_added_to_inventory() -> None:
    center = 0.53743
    provider = lambda value: np.asarray(
        [[1.0, 1.0, 0.0], [1.0, 1.0 + (value - center) ** 2 + 1.0e-3, 0.0], [0.0, 0.0, 1.0]]
    )
    result = beta_map._dense_local_candidate_search(
        provider,
        _local_settings(0.52, 0.55, 0.001),
        candidate_source_prefix="test_local_dense",
        block_family="full_6x6",
    )
    assert result["entries"] == []
    assert any(not row["accepted"] for row in result["candidate_rows"])


def test_confirmed_double_root_keeps_two_multiplicity_slots() -> None:
    center = 0.55
    provider = lambda value: np.diag([value - center, value - center, 1.0])
    result = beta_map._dense_local_candidate_search(
        provider,
        _local_settings(0.5, 0.6, 0.01),
        candidate_source_prefix="test_local_dense",
        block_family="full_6x6",
    )
    entries = result["entries"]
    assert len(entries) == 2
    assert [entry["multiplicity"] for entry in entries] == [2, 2]
    assert [entry["repeated_root_slot"] for entry in entries] == [1, 2]


def _entry(value: float, *, multiplicity: int = 1, slot: int = 1, source: str = "primary") -> dict[str, object]:
    return {
        "value": value,
        "source": source,
        "multiplicity": multiplicity,
        "multiplicity_source": "test",
        "block_family": "test",
        "nullity": multiplicity,
        "repeated_root_slot": slot,
    }


def test_merge_preserves_double_slots_and_resorts_inventory() -> None:
    primary = [_entry(float(index)) for index in range(1, 13)]
    local = [
        _entry(2.5, multiplicity=2, slot=1, source="R_local_dense"),
        _entry(2.5, multiplicity=2, slot=2, source="R_local_dense"),
    ]
    merged = beta_map.merge_primary_and_local_inventory(
        primary,
        local,
        lambda_min=2.4,
        lambda_max=2.6,
        settings=_local_settings(2.4, 2.6, 0.001),
    )
    assert [entry["value"] for entry in merged[:5]] == [1.0, 2.0, 2.5, 2.5, 3.0]
    assert [entry["repeated_root_slot"] for entry in merged[2:4]] == [1, 2]


def test_lower_recovered_root_shifts_only_following_sorted_positions() -> None:
    primary = [_entry(float(index)) for index in range(1, 13)]
    local = [_entry(4.5, source="R_local_dense")]
    merged = beta_map.merge_primary_and_local_inventory(
        primary,
        local,
        lambda_min=4.4,
        lambda_max=4.6,
        settings=_local_settings(4.4, 4.6, 0.001),
    )
    assert [entry["value"] for entry in merged[:4]] == [1.0, 2.0, 3.0, 4.0]
    assert [entry["value"] for entry in merged[4:7]] == [4.5, 5.0, 6.0]


def test_unresolved_inventory_uses_nan_without_interpolation() -> None:
    unresolved = beta_map.unresolved_inventory_from_position(
        [_entry(float(index)) for index in range(1, 13)],
        lambda_min=4.5,
    )
    assert [entry["value"] for entry in unresolved[:4]] == [1.0, 2.0, 3.0, 4.0]
    assert all(np.isnan(float(entry["value"])) for entry in unresolved[4:])


def test_p3_beta0_block_provenance_distinguishes_close_roots() -> None:
    full = {"entries": [_entry(8.20612), _entry(8.20700)]}
    bending = {"entries": [_entry(8.20700)]}
    axial = {"entries": [_entry(8.20612)]}
    entries, status = beta_map._annotate_p3_beta0_block_provenance(
        full,
        bending,
        axial,
        beta_map._local_refinement_settings(beta_map.refined_pilot_regions()[2]),
    )
    assert status == "resolved_distinct_pair"
    assert {entry["block_family"] for entry in entries} == {"axial_block", "bending_block"}


def test_refinement_cache_round_trip_is_resumable(tmp_path: Path) -> None:
    identity = {"cache_version": "test", "beta": 1.25}
    point = {"status": "resolved_distinct_pair", "entries": [_entry(1.0)]}
    path = tmp_path / "cache" / "point.json"
    beta_map._save_refined_cache(path, identity, point)
    assert beta_map._load_refined_cache(path, identity) == point
    assert not list(tmp_path.rglob("*.tmp.*"))


def test_refinement_path_contains_no_tracking_mac_shape_or_strict_calls() -> None:
    source = inspect.getsource(beta_map._run_refinement_point)
    for forbidden in (
        "TRACK.",
        "analytic_shape_vector_eta",
        "force_strict_verification",
        "resolve_general_spectrum",
        "resolve_matrix_spectrum",
        "continuation_seeds",
    ):
        assert forbidden not in source
