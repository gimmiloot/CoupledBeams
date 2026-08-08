from __future__ import annotations

import csv
import inspect
import math
from pathlib import Path
import sys

from PIL import Image
import pytest

from scripts.analysis.thickness_mismatch.audits import (
    audit_article_longitudinal_energy_examples as pilot,
)


REPO_ROOT = Path(__file__).resolve().parents[1]


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


@pytest.fixture(scope="module")
def generated(tmp_path_factory: pytest.TempPathFactory) -> dict[str, object]:
    output_dir = tmp_path_factory.mktemp("article_energy_triplet") / "pilot"
    records, selected, input_hashes = pilot.select_candidates(
        pilot.DEFAULT_FINALIZATION_DIR,
        pilot.DEFAULT_COMPACT_DIR,
    )
    result = pilot.run_pilot(
        pilot.DEFAULT_FINALIZATION_DIR,
        pilot.DEFAULT_COMPACT_DIR,
        output_dir,
    )
    after_hashes = {path: pilot._sha256(Path(path)) for path in input_hashes}
    return {
        "output_dir": output_dir,
        "records": records,
        "selected": selected,
        "input_hashes": input_hashes,
        "after_hashes": after_hashes,
        "result": result,
    }


def test_scientific_scope_is_only_isotropic_circular() -> None:
    assert pilot.SCIENTIFIC_SCOPE == "isotropic_circular_coupled_rods_eb_timoshenko"
    assert pilot.PILOT_SPECS == ((0.030, 5), (0.050, 4))


def test_anisotropic_modules_are_not_imported() -> None:
    assert not any("anisotropic_rods" in name for name in sys.modules)
    source = Path(pilot.__file__).read_text(encoding="utf-8")
    assert "scripts.analysis.anisotropic_rods" not in source


def test_candidates_are_resolved_main_grid_cases(generated: dict[str, object]) -> None:
    final_rows = _csv(pilot.DEFAULT_FINALIZATION_DIR / "final_case_certificates.csv")
    compact_rows = _csv(pilot.DEFAULT_COMPACT_DIR / "compact_index.csv")
    final_by_id = {row["case_id"]: row for row in final_rows}
    compact_by_id = {row["case_id"]: row for row in compact_rows}
    records = generated["records"]
    assert records
    for record in records:
        final = final_by_id[record["case_id"]]
        compact = compact_by_id[record["case_id"]]
        assert final["final_execution_status"].startswith("resolved")
        assert final["N_true"]
        assert not final["unresolved_reason"]
        assert compact["n_true_status"] == "exact"
        assert compact["required_guard_confirmed"] == "true"


def test_s3_regression_points_do_not_participate(generated: dict[str, object]) -> None:
    epsilons = {round(float(row["epsilon_0"]), 12) for row in generated["records"]}
    assert epsilons == {0.03, 0.05}
    assert round(0.024798906738281248, 12) not in epsilons
    assert round(0.029408510742187498, 12) not in epsilons


def test_requested_central_indices_are_used(generated: dict[str, object]) -> None:
    selected = {round(float(row["epsilon_0"]), 3): row for row in generated["selected"]}
    assert int(selected[0.030]["central_sorted_index"]) == 5
    assert int(selected[0.050]["central_sorted_index"]) == 4


def test_every_candidate_has_a_strict_local_delta_minimum(generated: dict[str, object]) -> None:
    for row in generated["records"]:
        assert float(row["delta_center"]) < float(row["delta_minus"])
        assert float(row["delta_center"]) < float(row["delta_plus"])


def test_candidate_selection_does_not_use_energy_values() -> None:
    selection_source = inspect.getsource(pilot.select_candidates)
    assert "timo_energy_partition" not in selection_source
    assert "chi_axial" not in selection_source
    assert "U_axial" not in selection_source


def test_selection_is_deterministic_and_has_two_geometries(
    generated: dict[str, object],
) -> None:
    selected_ids = [row["case_id"] for row in generated["selected"]]
    output_ids = [
        row["case_id"]
        for row in _csv(generated["output_dir"] / "selected_cases.csv")
    ]
    assert selected_ids == output_ids
    assert selected_ids == [
        "AUE_1e02916c7f48cdbdd59b",
        "AUE_68c469e17f3058b9bd6c",
    ]


def test_expected_six_energy_rows_and_fraction_identities(generated: dict[str, object]) -> None:
    rows = _csv(generated["output_dir"] / "energy_triplets.csv")
    assert len(rows) == 6
    assert {row["role"] for row in rows} == {
        "left_neighbor",
        "central_dip",
        "right_neighbor",
    }
    for row in rows:
        chi_a = float(row["chi_axial_Timo"])
        chi_b = float(row["chi_bending_Timo"])
        chi_s = float(row["chi_shear_Timo"])
        chi_bs = float(row["chi_bending_shear_Timo"])
        assert chi_a + chi_b + chi_s == pytest.approx(1.0, abs=pilot.ENERGY_SUM_TOL)
        assert chi_bs == pytest.approx(chi_b + chi_s, abs=1.0e-15)


def test_801_1601_integration_convergence_is_checked(generated: dict[str, object]) -> None:
    rows = _csv(generated["output_dir"] / "integration_convergence.csv")
    assert len(rows) == 6
    assert {(int(row["n_points_coarse"]), int(row["n_points_fine"])) for row in rows} == {
        (801, 1601)
    }
    for row in rows:
        assert row["passed"] == "true"
        assert max(
            float(row["delta_chi_a"]),
            float(row["delta_chi_b"]),
            float(row["delta_chi_s"]),
        ) <= pilot.INTEGRATION_TOL


def test_sign_flip_and_joint_residual_audits_pass(generated: dict[str, object]) -> None:
    rows = _csv(generated["output_dir"] / "mode_residual_audit.csv")
    assert len(rows) == 6
    assert sum(row["null_vector_source"].endswith("fallback") for row in rows) == 1
    for row in rows:
        assert row["sign_flip_pass"] == "true"
        assert float(row["sign_flip_max_abs_delta"]) <= pilot.ENERGY_SUM_TOL
        assert row["svd_pass"] == "true"
        assert row["null_residual_pass"] == "true"
        assert row["clamp_pass"] == "true"
        assert row["joint_kinematic_pass"] == "true"
        assert row["joint_force_pass"] == "true"
        assert row["passed"] == "true"


def test_central_energy_contribution_exceeds_both_neighbors(
    generated: dict[str, object],
) -> None:
    rows = _csv(generated["output_dir"] / "energy_triplets.csv")
    for epsilon, central_index in pilot.PILOT_SPECS:
        level = [row for row in rows if math.isclose(float(row["epsilon_0"]), epsilon)]
        by_index = {int(row["sorted_index"]): float(row["chi_axial_Timo"]) for row in level}
        assert by_index[central_index] > by_index[central_index - 1]
        assert by_index[central_index] > by_index[central_index + 1]


def test_forbidden_operation_counters_are_zero(generated: dict[str, object]) -> None:
    rows = _csv(generated["output_dir"] / "operation_counters.csv")
    assert len(rows) == 1
    row = rows[0]
    for field in (
        "root_solver_calls",
        "strict_solver_calls",
        "family_detector_calls",
        "local_repair_calls",
        "shape_plot_calls",
        "MAC_calls",
        "FEM_calls",
        "anisotropic_calls",
    ):
        assert int(row[field]) == 0
    assert int(row["matrix_evaluator_calls"]) == 12
    assert int(row["null_vector_reconstructions"]) == 7
    assert int(row["energy_integral_evaluations"]) == 18


def test_output_contains_one_600_dpi_png_and_no_pdf(generated: dict[str, object]) -> None:
    output_dir = generated["output_dir"]
    pngs = list(output_dir.glob("*.png"))
    assert [path.name for path in pngs] == ["energy_triplet_comparison.png"]
    assert not list(output_dir.rglob("*.pdf"))
    with Image.open(pngs[0]) as image:
        dpi = image.info["dpi"]
    assert dpi[0] == pytest.approx(600, abs=0.1)
    assert dpi[1] == pytest.approx(600, abs=0.1)


def test_all_gates_and_scientific_status_pass(generated: dict[str, object]) -> None:
    result = generated["result"]
    assert result["scientific_status"] == "hypothesis_supported_by_two_examples"
    assert all(result["gates"].values())


def test_compact_inputs_are_unchanged(generated: dict[str, object]) -> None:
    assert generated["input_hashes"] == generated["after_hashes"]
    rows = _csv(generated["output_dir"] / "source_integrity.csv")
    assert rows
    assert all(row["unchanged"] == "true" for row in rows)


@pytest.fixture(scope="module")
def coupled_generated(tmp_path_factory: pytest.TempPathFactory) -> dict[str, object]:
    output_dir = tmp_path_factory.mktemp("article_energy_coupled_triplet") / "pilot"
    beta0_hashes_before = pilot._directory_hashes(pilot.DEFAULT_BETA0_CONTROL_DIR)
    records, selected, input_hashes, context = pilot.select_coupled_candidates(
        pilot.DEFAULT_FINALIZATION_DIR,
        pilot.DEFAULT_COMPACT_DIR,
        pilot.DEFAULT_BETA0_CONTROL_DIR,
    )
    result = pilot.run_coupled_pilot(
        finalization_dir=pilot.DEFAULT_FINALIZATION_DIR,
        compact_dir=pilot.DEFAULT_COMPACT_DIR,
        beta0_control_dir=pilot.DEFAULT_BETA0_CONTROL_DIR,
        output_dir=output_dir,
    )
    return {
        "output_dir": output_dir,
        "records": records,
        "selected": selected,
        "context": context,
        "input_hashes": input_hashes,
        "input_hashes_after": {
            path: pilot._sha256(Path(path)) for path in input_hashes
        },
        "beta0_hashes_before": beta0_hashes_before,
        "beta0_hashes_after": pilot._directory_hashes(pilot.DEFAULT_BETA0_CONTROL_DIR),
        "result": result,
    }


def test_coupled_selection_is_nonzero_angle_frequency_only(
    coupled_generated: dict[str, object],
) -> None:
    records = coupled_generated["records"]
    assert records
    for record in records:
        assert float(record["beta_deg"]) >= pilot.COUPLED_MIN_BETA_DEG
        assert float(record["delta_center"]) < float(record["delta_minus"])
        assert float(record["delta_center"]) < float(record["delta_plus"])
    source = inspect.getsource(pilot.select_coupled_candidates)
    assert "timo_energy_partition" not in source
    assert "chi_axial" not in source
    assert "U_axial" not in source
    assert "AUE_70fc3defd1f14b5706da" not in source
    assert "AUE_a5f9d866334acfbeb71d" not in source
    assert coupled_generated["context"]["selection_uses_energy_values"] is False


def test_coupled_pool_policy_and_regression_ids(
    coupled_generated: dict[str, object],
) -> None:
    selected = {
        round(float(row["epsilon_0"]), 3): row
        for row in coupled_generated["selected"]
    }
    assert len(selected) == 2
    assert selected[0.030]["case_id"] == "AUE_70fc3defd1f14b5706da"
    assert selected[0.050]["case_id"] == "AUE_a5f9d866334acfbeb71d"
    assert selected[0.030]["selection_pool"] == "nonzero_angle_smax_le_0p1"
    assert float(selected[0.030]["s_max"]) <= pilot.THINNESS_LIMIT
    assert selected[0.050]["selection_pool"] == "nonzero_angle_all"
    assert float(selected[0.050]["s_max"]) > pilot.THINNESS_LIMIT
    assert (
        coupled_generated["context"]["main_grid_smax_le_0p1_count_epsilon_0p050"]
        == 0
    )
    assert all(float(row["beta_deg"]) == 15.0 for row in selected.values())


def test_coupled_selection_is_deterministic_and_uses_requested_indices(
    coupled_generated: dict[str, object],
) -> None:
    selected = coupled_generated["selected"]
    output = _csv(coupled_generated["output_dir"] / "selected_coupled_cases.csv")
    assert [row["case_id"] for row in selected] == [row["case_id"] for row in output]
    by_epsilon = {round(float(row["epsilon_0"]), 3): row for row in selected}
    assert int(by_epsilon[0.030]["central_sorted_index"]) == 5
    assert int(by_epsilon[0.050]["central_sorted_index"]) == 4
    top_rows = _csv(coupled_generated["output_dir"] / "coupled_candidate_ranking_top10.csv")
    counts: dict[tuple[str, str], int] = {}
    for row in top_rows:
        key = (row["epsilon_0"], row["selection_pool"])
        counts[key] = counts.get(key, 0) + 1
    assert sorted(counts.values()) == [3, 10, 10]


def test_coupled_energy_has_exactly_six_rows_and_fraction_identities(
    coupled_generated: dict[str, object],
) -> None:
    rows = _csv(coupled_generated["output_dir"] / "coupled_energy_triplets.csv")
    assert len(rows) == 6
    assert {
        (round(float(row["epsilon_0"]), 3), int(row["sorted_index"])) for row in rows
    } == {(0.030, 4), (0.030, 5), (0.030, 6), (0.050, 3), (0.050, 4), (0.050, 5)}
    for row in rows:
        chi_a = float(row["chi_axial_Timo"])
        chi_b = float(row["chi_bending_Timo"])
        chi_s = float(row["chi_shear_Timo"])
        assert chi_a + chi_b + chi_s == pytest.approx(1.0, abs=pilot.ENERGY_SUM_TOL)
        assert float(row["chi_bending_shear_Timo"]) == pytest.approx(
            chi_b + chi_s, abs=1.0e-15
        )
        assert row["article_ready_status"] == "article_ready"


def test_coupled_center_has_higher_axial_contribution_than_neighbors(
    coupled_generated: dict[str, object],
) -> None:
    rows = _csv(coupled_generated["output_dir"] / "coupled_energy_triplets.csv")
    for epsilon, central_index in pilot.PILOT_SPECS:
        level = {
            int(row["sorted_index"]): float(row["chi_axial_Timo"])
            for row in rows
            if math.isclose(float(row["epsilon_0"]), epsilon)
        }
        assert level[central_index] > level[central_index - 1]
        assert level[central_index] > level[central_index + 1]


def test_coupled_convergence_residual_and_sign_checks_pass(
    coupled_generated: dict[str, object],
) -> None:
    convergence = _csv(coupled_generated["output_dir"] / "integration_convergence.csv")
    residuals = _csv(coupled_generated["output_dir"] / "mode_residual_audit.csv")
    assert len(convergence) == len(residuals) == 6
    for row in convergence:
        assert (int(row["n_points_coarse"]), int(row["n_points_fine"])) == (801, 1601)
        assert row["passed"] == "true"
        assert max(float(row[name]) for name in ("delta_chi_a", "delta_chi_b", "delta_chi_s")) <= pilot.INTEGRATION_TOL
    for row in residuals:
        assert all(
            row[name] == "true"
            for name in (
                "svd_pass",
                "null_residual_pass",
                "clamp_pass",
                "joint_kinematic_pass",
                "joint_force_pass",
                "sign_flip_pass",
                "passed",
            )
        )
        assert float(row["sign_flip_max_abs_delta"]) <= pilot.ENERGY_SUM_TOL


def test_beta0_control_is_read_zero_solve_and_preserved(
    coupled_generated: dict[str, object],
) -> None:
    control = _csv(coupled_generated["output_dir"] / "beta0_control_reference.csv")
    assert len(control) == 3
    assert {int(row["sorted_index"]) for row in control} == {4, 5, 6}
    assert all(float(row["beta_deg"]) == 0.0 for row in control)
    assert coupled_generated["beta0_hashes_before"] == coupled_generated["beta0_hashes_after"]


def test_coupled_operation_counts_and_gates(
    coupled_generated: dict[str, object],
) -> None:
    counts = _csv(coupled_generated["output_dir"] / "operation_counts.csv")[0]
    for field in (
        "root_solver_calls",
        "strict_solver_calls",
        "family_detector_calls",
        "local_repair_calls",
        "shape_plot_calls",
        "MAC_calls",
        "FEM_calls",
        "anisotropic_calls",
    ):
        assert int(counts[field]) == 0
    assert int(counts["matrix_evaluator_calls"]) == 12
    assert int(counts["null_vector_reconstructions"]) == 6
    assert int(counts["energy_integral_evaluations"]) == 18
    result = coupled_generated["result"]
    assert result["scientific_status"] == "coupled_hypothesis_supported_by_two_examples"
    assert all(result["gates"].values())


def test_coupled_output_has_one_600_dpi_png_and_no_pdf(
    coupled_generated: dict[str, object],
) -> None:
    output_dir = coupled_generated["output_dir"]
    pngs = list(output_dir.glob("*.png"))
    assert [path.name for path in pngs] == ["coupled_energy_triplet_comparison.png"]
    assert not list(output_dir.rglob("*.pdf"))
    with Image.open(pngs[0]) as image:
        dpi = image.info["dpi"]
    assert dpi[0] == pytest.approx(600, abs=0.1)
    assert dpi[1] == pytest.approx(600, abs=0.1)


def test_coupled_compact_certificates_are_unchanged(
    coupled_generated: dict[str, object],
) -> None:
    assert coupled_generated["input_hashes"] == coupled_generated["input_hashes_after"]


@pytest.fixture(scope="module")
def additional_angle_generated(
    tmp_path_factory: pytest.TempPathFactory,
) -> dict[str, object]:
    base = tmp_path_factory.mktemp("article_energy_additional_angle")
    selection_dir = base / "selection"
    output_dir = base / "pilot"
    beta0_hashes_before = pilot._directory_hashes(pilot.DEFAULT_BETA0_CONTROL_DIR)
    coupled_hashes_before = pilot._directory_hashes(pilot.DEFAULT_COUPLED_CONTROL_DIR)
    candidates, selected, input_hashes, context = pilot.select_additional_angle_candidate(
        pilot.DEFAULT_FINALIZATION_DIR,
        pilot.DEFAULT_COMPACT_DIR,
        pilot.DEFAULT_COUPLED_CONTROL_DIR,
    )
    selection_result = pilot.run_additional_angle_pilot(
        finalization_dir=pilot.DEFAULT_FINALIZATION_DIR,
        compact_dir=pilot.DEFAULT_COMPACT_DIR,
        beta0_control_dir=pilot.DEFAULT_BETA0_CONTROL_DIR,
        coupled_control_dir=pilot.DEFAULT_COUPLED_CONTROL_DIR,
        output_dir=selection_dir,
        selection_only=True,
    )
    result = pilot.run_additional_angle_pilot(
        finalization_dir=pilot.DEFAULT_FINALIZATION_DIR,
        compact_dir=pilot.DEFAULT_COMPACT_DIR,
        beta0_control_dir=pilot.DEFAULT_BETA0_CONTROL_DIR,
        coupled_control_dir=pilot.DEFAULT_COUPLED_CONTROL_DIR,
        output_dir=output_dir,
    )
    return {
        "selection_dir": selection_dir,
        "output_dir": output_dir,
        "candidates": candidates,
        "selected": selected,
        "context": context,
        "selection_result": selection_result,
        "result": result,
        "input_hashes": input_hashes,
        "input_hashes_after": {
            path: pilot._sha256(Path(path)) for path in input_hashes
        },
        "beta0_hashes_before": beta0_hashes_before,
        "beta0_hashes_after": pilot._directory_hashes(pilot.DEFAULT_BETA0_CONTROL_DIR),
        "coupled_hashes_before": coupled_hashes_before,
        "coupled_hashes_after": pilot._directory_hashes(pilot.DEFAULT_COUPLED_CONTROL_DIR),
    }


def test_additional_angle_selection_is_zero_solve_frequency_only(
    additional_angle_generated: dict[str, object],
) -> None:
    selection_result = additional_angle_generated["selection_result"]
    assert selection_result["selection_only"] is True
    assert all(value == 0 for value in selection_result["operation_counts"].values())
    source = inspect.getsource(pilot.select_additional_angle_candidate)
    assert "timo_energy_partition" not in source
    assert "chi_axial" not in source
    assert "U_axial" not in source
    assert "AUE_ac5482dd4f84aa6d7dfc" not in source
    assert additional_angle_generated["context"]["selection_uses_energy_values"] is False


def test_additional_angle_pool_constraints_and_case_id_deduplication(
    additional_angle_generated: dict[str, object],
) -> None:
    candidates = additional_angle_generated["candidates"]
    assert candidates
    assert len({row["case_id"] for row in candidates}) == len(candidates)
    for row in candidates:
        assert math.isclose(float(row["epsilon_0"]), 0.030)
        assert int(row["central_sorted_index"]) == 5
        assert float(row["beta_deg"]) >= 30.0
        assert float(row["s_max"]) <= pilot.THINNESS_LIMIT
        assert float(row["delta_center"]) < float(row["delta_minus"])
        assert float(row["delta_center"]) < float(row["delta_plus"])
    context = additional_angle_generated["context"]
    assert context["ranking_source_row_count"] == 156
    assert context["ranking_unique_case_id_count"] == 153
    assert context["ranking_duplicate_rows_removed"] == 3


def test_additional_angle_regression_selection_and_fixed_geometry(
    additional_angle_generated: dict[str, object],
) -> None:
    selected = additional_angle_generated["selected"]
    assert selected["case_id"] == "AUE_ac5482dd4f84aa6d7dfc"
    assert float(selected["beta_deg"]) == 30.0
    assert float(selected["mu"]) == 0.5
    assert float(selected["eta"]) == 0.5
    combined = _csv(
        additional_angle_generated["output_dir"] / "combined_article_energy_triplets.csv"
    )
    rows15 = [row for row in combined if float(row["epsilon_0"]) == 0.03 and float(row["beta_deg"]) == 15.0]
    rows30 = [row for row in combined if float(row["epsilon_0"]) == 0.03 and float(row["beta_deg"]) == 30.0]
    assert len(rows15) == len(rows30) == 3
    assert {
        (float(row["epsilon_0"]), float(row["mu"]), float(row["eta"])) for row in rows15 + rows30
    } == {(0.03, 0.5, 0.5)}
    assert {float(row["beta_deg"]) for row in rows15 + rows30} == {15.0, 30.0}


def test_additional_angle_has_exactly_three_energy_rows_and_identities(
    additional_angle_generated: dict[str, object],
) -> None:
    rows = _csv(
        additional_angle_generated["output_dir"] / "additional_angle_energy_triplet.csv"
    )
    assert len(rows) == 3
    assert [int(row["sorted_index"]) for row in rows] == [4, 5, 6]
    for row in rows:
        chi_a = float(row["chi_axial_Timo"])
        chi_b = float(row["chi_bending_Timo"])
        chi_s = float(row["chi_shear_Timo"])
        assert chi_a + chi_b + chi_s == pytest.approx(1.0, abs=pilot.ENERGY_SUM_TOL)
        assert float(row["chi_bending_shear_Timo"]) == pytest.approx(
            chi_b + chi_s, abs=1.0e-15
        )
        assert row["article_ready_status"] == "article_ready"


def test_additional_angle_center_exceeds_both_neighbors(
    additional_angle_generated: dict[str, object],
) -> None:
    rows = _csv(
        additional_angle_generated["output_dir"] / "additional_angle_energy_triplet.csv"
    )
    axial = {int(row["sorted_index"]): float(row["chi_axial_Timo"]) for row in rows}
    assert axial[5] > axial[4]
    assert axial[5] > axial[6]


def test_additional_angle_convergence_residual_and_sign_checks_pass(
    additional_angle_generated: dict[str, object],
) -> None:
    output_dir = additional_angle_generated["output_dir"]
    convergence = _csv(output_dir / "integration_convergence.csv")
    residuals = _csv(output_dir / "mode_residual_audit.csv")
    assert len(convergence) == len(residuals) == 3
    for row in convergence:
        assert (int(row["n_points_coarse"]), int(row["n_points_fine"])) == (801, 1601)
        assert row["passed"] == "true"
        assert max(float(row[name]) for name in ("delta_chi_a", "delta_chi_b", "delta_chi_s")) <= pilot.INTEGRATION_TOL
    for row in residuals:
        assert all(
            row[name] == "true"
            for name in (
                "svd_pass",
                "null_residual_pass",
                "clamp_pass",
                "joint_kinematic_pass",
                "joint_force_pass",
                "sign_flip_pass",
                "passed",
            )
        )
        assert float(row["sign_flip_max_abs_delta"]) <= pilot.ENERGY_SUM_TOL


def test_additional_angle_combined_and_figure_source_row_counts(
    additional_angle_generated: dict[str, object],
) -> None:
    output_dir = additional_angle_generated["output_dir"]
    combined = _csv(output_dir / "combined_article_energy_triplets.csv")
    compact_source = _csv(output_dir / "article_energy_angle_comparison_source.csv")
    extended_source = _csv(output_dir / "article_energy_examples_extended_source.csv")
    assert len(combined) == 9
    assert len(compact_source) == 6
    assert len(extended_source) == 9
    assert not any(float(row["beta_deg"]) == 0.0 for row in combined)


def test_additional_angle_output_has_two_600_dpi_pngs_and_no_pdf(
    additional_angle_generated: dict[str, object],
) -> None:
    output_dir = additional_angle_generated["output_dir"]
    pngs = sorted(output_dir.glob("*.png"))
    assert [path.name for path in pngs] == [
        "article_energy_angle_comparison.png",
        "article_energy_examples_extended.png",
    ]
    assert not list(output_dir.rglob("*.pdf"))
    for path in pngs:
        with Image.open(path) as image:
            dpi = image.info["dpi"]
        assert dpi[0] == pytest.approx(600, abs=0.1)
        assert dpi[1] == pytest.approx(600, abs=0.1)


def test_additional_angle_prior_results_and_compact_inputs_are_unchanged(
    additional_angle_generated: dict[str, object],
) -> None:
    assert additional_angle_generated["beta0_hashes_before"] == additional_angle_generated["beta0_hashes_after"]
    assert additional_angle_generated["coupled_hashes_before"] == additional_angle_generated["coupled_hashes_after"]
    assert additional_angle_generated["input_hashes"] == additional_angle_generated["input_hashes_after"]


def test_additional_angle_operation_counts_gates_and_root_immutability(
    additional_angle_generated: dict[str, object],
) -> None:
    output_dir = additional_angle_generated["output_dir"]
    counts = _csv(output_dir / "operation_counts.csv")[0]
    for field in (
        "root_solver_calls",
        "strict_solver_calls",
        "family_detector_calls",
        "local_repair_calls",
        "shape_plot_calls",
        "MAC_calls",
        "FEM_calls",
        "anisotropic_calls",
    ):
        assert int(counts[field]) == 0
    assert int(counts["matrix_evaluator_calls"]) == 6
    assert int(counts["null_vector_reconstructions"]) == 3
    assert int(counts["energy_integral_evaluations"]) == 9
    result = additional_angle_generated["result"]
    assert result["root_immutability_pass"] is True
    assert result["scientific_status"] == "additional_angle_hypothesis_supported"
    assert all(result["gates"].values())
