from __future__ import annotations

import csv
import inspect
import json
from pathlib import Path
import shutil
import sys

import pytest

from scripts.analysis.thickness_mismatch.audits import (
    audit_article_energy_local_minima_grid as grid,
)


def _csv(path: Path, *, delimiter: str = ",") -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


@pytest.fixture(scope="module")
def generated(tmp_path_factory: pytest.TempPathFactory) -> dict[str, object]:
    output_dir = tmp_path_factory.mktemp("article_energy_grid") / "smoke"
    prior_before = grid.prior_fingerprints()
    selection = grid.select_frequency_triplets(
        grid.DEFAULT_SOURCE_DIR,
        grid.DEFAULT_COMPACT_DIR,
        output_dir,
        resume=False,
        smoke=True,
        smoke_cases=4,
    )
    energy = grid.compute_energy_phase(output_dir, resume=True, workers=1)
    summary = grid.summarize_phase(output_dir)
    key_names = (
        "frequency_selected_triplets.csv",
        "unique_selected_modes.csv",
        "mode_energy_801.csv",
        "mode_energy_1601.csv",
        "mode_energy_convergence.csv",
        "energy_triplet_results.csv",
        "summary_overall.csv",
    )
    hashes_before_resume = {
        name: grid.sha256_file(output_dir / name) for name in key_names
    }
    resume_result = grid.compute_energy_phase(output_dir, resume=True, workers=2)
    grid.summarize_phase(output_dir)
    hashes_after_resume = {
        name: grid.sha256_file(output_dir / name) for name in key_names
    }
    return {
        "output_dir": output_dir,
        "selection": selection,
        "energy": energy,
        "summary": summary,
        "resume_result": resume_result,
        "hashes_before_resume": hashes_before_resume,
        "hashes_after_resume": hashes_after_resume,
        "prior_before": prior_before,
        "prior_after": grid.prior_fingerprints(),
    }


def test_scientific_scope_is_isotropic_circular_only() -> None:
    assert grid.SCIENTIFIC_SCOPE == "isotropic_circular_coupled_rods_eb_timoshenko"
    assert grid.K_MAX == 10
    assert grid.LOCAL_MIN_K == tuple(range(2, 10))


def test_anisotropic_modules_are_not_imported() -> None:
    source = Path(grid.__file__).read_text(encoding="utf-8")
    assert "scripts.analysis.anisotropic_rods" not in source
    assert not any("anisotropic_rods" in name for name in sys.modules)


def test_phase_a_selection_source_does_not_read_energy_tables() -> None:
    source = inspect.getsource(grid.select_frequency_triplets)
    discovery_source = inspect.getsource(grid.load_discovery_cases)
    forbidden = (
        "timo_energy_partition",
        "mode_energy_801",
        "mode_energy_1601",
        "energy_triplets.csv",
        "chi_axial_Timo",
        "U_axial_Timo",
    )
    assert all(item not in source for item in forbidden)
    assert all(item not in discovery_source for item in forbidden)


def test_phase_b_verifies_the_fixed_selection_hash(
    generated: dict[str, object], tmp_path: Path
) -> None:
    source = generated["output_dir"]
    shutil.copy(source / "frequency_selected_triplets.csv", tmp_path)
    shutil.copy(source / "frequency_selection_metadata.json", tmp_path)
    selection = tmp_path / "frequency_selected_triplets.csv"
    selection.write_bytes(selection.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="SHA-256"):
        grid.verify_frequency_selection(tmp_path)


def test_selected_k_and_frequency_local_metrics(generated: dict[str, object]) -> None:
    rows = _csv(generated["output_dir"] / "frequency_selected_triplets.csv")
    assert rows
    assert all(2 <= int(row["k"]) <= 9 for row in rows)
    for row in rows:
        left = float(row["delta_f_k_minus_1"])
        center = float(row["delta_f_k"])
        right = float(row["delta_f_k_plus_1"])
        metrics = grid.frequency_local_characteristics(left, center, right)
        assert metrics["strict_local_minimum"] is True
        assert float(row["D_f"]) == pytest.approx(0.5 * (left + right) - center)


def test_energy_local_metrics_are_threshold_free() -> None:
    metrics = grid.energy_local_characteristics(0.10, 0.25, 0.20)
    assert metrics["D_a"] == pytest.approx(0.10)
    assert metrics["axial_above_neighbor_mean"] is True
    assert metrics["axial_strict_local_max"] is True
    not_strict = grid.energy_local_characteristics(0.30, 0.25, 0.10)
    assert not_strict["D_a"] == pytest.approx(0.05)
    assert not_strict["axial_strict_local_max"] is False


def test_overlapping_triplets_deduplicate_modes(generated: dict[str, object]) -> None:
    rows = _csv(generated["output_dir"] / "frequency_selected_triplets.csv")
    by_case: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_case.setdefault(row["case_id"], []).append(row)
    overlapping = next(
        values
        for values in by_case.values()
        if len(values) >= 2 and abs(int(values[0]["k"]) - int(values[1]["k"])) <= 2
    )[:2]
    unique = grid.build_unique_modes(overlapping)
    requested = {
        (row["case_id"], rank)
        for row in overlapping
        for rank in range(int(row["k"]) - 1, int(row["k"]) + 2)
    }
    assert len(unique) == len(requested) == 5


@pytest.mark.parametrize(
    ("overrides", "expected_primary", "expected_named"),
    [
        ({"beta_deg": 0.0}, False, "beta0_control"),
        ({"s_max": 0.1000001}, False, "extended_range"),
        ({"discovery_case": "true"}, False, "discovery_cohort"),
        ({"cluster_flag": "true"}, False, "cluster_affected"),
    ],
)
def test_nonconfirmatory_cases_are_retained_in_named_cohorts(
    overrides: dict[str, object], expected_primary: bool, expected_named: str
) -> None:
    row: dict[str, object] = {
        "beta_deg": 15.0,
        "s_max": 0.09,
        "mu": 0.3,
        "eta": 0.0,
        "cluster_flag": "false",
        "discovery_case": "false",
    }
    row.update(overrides)
    flags = grid.cohort_flags(row, energy_quality_pass=True)
    assert flags["primary_confirmatory"] is expected_primary
    assert flags[expected_named] is True


def test_strict_family_holdout_excludes_complete_mu_eta_family() -> None:
    row = {
        "beta_deg": 15.0,
        "s_max": 0.09,
        "mu": 0.5,
        "eta": 0.5,
        "cluster_flag": "false",
        "discovery_case": "false",
    }
    flags = grid.cohort_flags(row, energy_quality_pass=True)
    assert flags["primary_confirmatory"] is True
    assert flags["strict_family_holdout"] is False


def test_batch_energy_wrapper_matches_scalar_helper(generated: dict[str, object]) -> None:
    mode = grid.build_unique_modes(
        _csv(generated["output_dir"] / "frequency_selected_triplets.csv")[:1]
    )[0]
    coefficients = grid.PILOT.TIMO.timo_mode_coefficients(
        float(mode["Lambda_Timo"]),
        float(mode["beta_deg"]),
        float(mode["mu"]),
        float(mode["epsilon_0"]),
        float(mode["eta"]),
    ).coeff
    batch = grid.evaluate_energy_batch(mode, coefficients, (801,))[801]
    scalar = grid.PILOT.TIMO.timo_energy_partition(
        float(mode["Lambda_Timo"]),
        float(mode["beta_deg"]),
        float(mode["mu"]),
        float(mode["epsilon_0"]),
        float(mode["eta"]),
        coeff=coefficients,
        n_points=801,
    )
    for field in ("U_a_total", "U_b_total", "U_s_total", "U_total"):
        assert batch[field] == pytest.approx(scalar[field], rel=1.0e-14, abs=0.0)


def test_energy_identities_and_grid_checks_are_recorded(generated: dict[str, object]) -> None:
    for name in ("mode_energy_801.csv", "mode_energy_1601.csv"):
        rows = _csv(generated["output_dir"] / name)
        assert rows
        for row in rows:
            chi_a = float(row["chi_a"])
            chi_b = float(row["chi_b"])
            chi_s = float(row["chi_s"])
            assert chi_a + chi_b + chi_s == pytest.approx(1.0, abs=grid.PILOT.ENERGY_SUM_TOL)
            assert float(row["chi_bs"]) == pytest.approx(chi_b + chi_s, abs=1.0e-15)
    convergence = _csv(generated["output_dir"] / "mode_energy_convergence.csv")
    assert convergence
    assert all(row["integration_pass"] == "true" for row in convergence)
    assert all(row["sign_flip_pass"] == "true" for row in convergence)
    assert all(row["svd_pass"] == "true" for row in convergence)


def test_cluster_rows_are_not_silently_dropped(generated: dict[str, object]) -> None:
    all_rows = _csv(generated["output_dir"] / "energy_triplet_results.csv")
    cluster_rows = _csv(generated["output_dir"] / "cluster_affected_results.csv")
    assert [(row["case_id"], row["k"]) for row in cluster_rows] == [
        (row["case_id"], row["k"])
        for row in all_rows
        if row["cluster_affected"] == "true"
    ]


def test_resume_is_deterministic_and_does_not_duplicate_rows(
    generated: dict[str, object],
) -> None:
    assert generated["hashes_before_resume"] == generated["hashes_after_resume"]
    assert generated["resume_result"]["cache_hit_count"] == generated["energy"][
        "unique_mode_count"
    ]
    unique = _csv(generated["output_dir"] / "unique_selected_modes.csv")
    assert len(unique) == len({(row["case_id"], row["sorted_k"]) for row in unique})


def test_prior_pilots_are_immutable_and_operations_are_zero(
    generated: dict[str, object],
) -> None:
    assert generated["prior_before"] == generated["prior_after"]
    metadata = json.loads(
        (generated["output_dir"] / "run_metadata.json").read_text(encoding="utf-8")
    )
    for name in (
        "root_solver_calls",
        "strict_solver_calls",
        "family_detector_calls",
        "local_repair_calls",
        "MAC_calls",
        "shape_plot_calls",
        "FEM_calls",
        "anisotropic_calls",
    ):
        assert metadata["operation_counts"][name] == 0


def test_all_phase_a_inputs_remain_hash_identical(generated: dict[str, object]) -> None:
    metadata = json.loads(
        (generated["output_dir"] / "frequency_selection_metadata.json").read_text(
            encoding="utf-8"
        )
    )
    for raw_path, expected_hash in metadata["input_file_hashes"].items():
        assert grid.sha256_file(Path(raw_path)) == expected_hash


def test_discovery_inventory_and_compact_article_table_are_emitted(
    generated: dict[str, object],
) -> None:
    discovery = _csv(generated["output_dir"] / "discovery_cases_results.csv")
    assert len(discovery) == 5
    assert all(row["first10_coverage_status"] for row in discovery)
    compact = _csv(
        generated["output_dir"] / "article_energy_local_minima_summary_overall.tsv",
        delimiter="\t",
    )
    assert len(compact) == 1
    assert compact[0]["совокупность"] == "Подтверждающая совокупность"
