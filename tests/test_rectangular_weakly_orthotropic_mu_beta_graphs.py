"""Fast contract tests for the finite RLB-2D ``mu``/``beta`` graphs.

The tests deliberately do not execute either full sweep.  They check the
fixed normalization and geometry contracts, the weak four-ply arm factory,
the additive EB adapter at ``tau=0``, and the pure plotting path.  When the
ignored full result directory has been finalized, a final test audits all six
root tables at the requested representative points and verifies both PNGs.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from matplotlib.legend import Legend
import numpy as np
from numpy.testing import assert_allclose
import pytest

from scripts.analysis.laminated_beams import (
    compare_rectangular_weakly_orthotropic_mu_beta_graphs as target,
)


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = (
    ROOT
    / "scripts"
    / "analysis"
    / "laminated_beams"
    / "compare_rectangular_weakly_orthotropic_mu_beta_graphs.py"
)
RESULT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "rectangular_weakly_orthotropic_graphs_mu_beta"
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return [dict(row) for row in csv.DictReader(stream)]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest().upper()


@pytest.fixture(scope="module")
def small_contract() -> dict[str, Any]:
    return target.build_case_contract(
        mu_values=[0.0, 0.4, 0.8],
        beta_values=[0.0, 15.0, 90.0],
    )


def test_old_lambda_normalization_uses_only_the_fixed_reference() -> None:
    expected_scale = target.L_REFERENCE**2 * math.sqrt(
        target.RHO0 * target.A0 / (target.E0 * target.I_Y0)
    )
    assert_allclose(
        expected_scale,
        target.L_REFERENCE**2
        * math.sqrt(12.0 * target.RHO0 / (target.E0 * target.H0**2)),
        rtol=4.0e-16,
        atol=0.0,
    )
    assert_allclose(
        target.reference_omega_scale(), expected_scale, rtol=2.0e-16, atol=0.0
    )

    omega = 0.731
    Omega = target.omega_to_Omega(omega)
    Lambda = target.Omega_to_Lambda(Omega)
    assert_allclose(Omega, omega * expected_scale, rtol=2.0e-16, atol=0.0)
    assert_allclose(Lambda, math.sqrt(Omega), rtol=0.0, atol=0.0)
    assert_allclose(target.old_Lambda(omega), Lambda, rtol=0.0, atol=0.0)
    assert_allclose(Lambda**2, Omega, rtol=2.0e-16, atol=0.0)


def test_exact_default_grids_and_finite_sweep_contract(
    small_contract: dict[str, Any],
) -> None:
    assert_allclose(
        target.mu_grid(), np.arange(41, dtype=float) / 50.0, rtol=0.0, atol=0.0
    )
    assert_allclose(
        target.beta_grid(), np.arange(91, dtype=float), rtol=0.0, atol=0.0
    )
    assert small_contract["sweeps"][target.SWEEP_MU] == {
        "axis": "mu",
        "values": [0.0, 0.4, 0.8],
        "step": 0.4,
        "beta_deg": 15.0,
        "tau": 0.0,
        "requested_step": 0.01,
        "allowed_runtime_fallback_step": 0.02,
        "fallback_used": True,
        "fallback_reason": target.MU_FALLBACK_REASON,
        "fallback_decision_used_spectrum": False,
    }
    assert small_contract["sweeps"][target.SWEEP_BETA] == {
        "axis": "beta_deg",
        "values": [0.0, 15.0, 90.0],
        "step": 15.0,
        "mu": 0.5,
        "tau": 0.2,
    }
    assert small_contract["plotted_sorted_positions"] == 8
    assert small_contract["output_guard_position"] == 9
    assert small_contract["modal_descendant_tracking"] is False
    assert small_contract["inter_model_relative_differences_computed"] is False


def test_reference_and_weak_case_A_material_contract(
    small_contract: dict[str, Any],
) -> None:
    assert small_contract["reference_constants"] == {
        "E0": 1.0,
        "nu0": 0.3,
        "rho0": 1.0,
        "b0": 0.2,
        "h0": 0.05,
        "l": 1.0,
        "L_total": 2.0,
        "K": 5.0 / 6.0,
        "A0": target.A0,
        "I_y0": target.I_Y0,
    }
    assert small_contract["new_RLB_lamina"] == {
        "case_id": "A",
        "delta": 0.1,
        "E1": 1.1,
        "E2": 0.9,
        "nu12": 0.3,
        "G12": 1.0 / 2.6,
        "G13": 1.0 / 2.6,
        "G23": 1.0 / 2.6,
        "rho": 1.0,
    }
    assert small_contract["new_RLB_laminate"] == {
        "number_of_plies_per_arm": 4,
        "stacking_sequence_deg": [0.0, 90.0, 90.0, 0.0],
        "ply_thickness_arm_i": "h_i/4",
        "one_ply_shortcut": False,
        "pipeline": (
            "ply properties->Q->Qbar->A/B/D->shear/mass"
            "->beam reduction->coupled determinant"
        ),
    }
    assert small_contract["old_models"] == {
        "EB": "isotropic rectangular baseline with actual arm sections",
        "old_Timoshenko": (
            "isotropic rectangular baseline with actual arm sections"
        ),
        "equivalent_isotropic_fitting": False,
    }


@pytest.mark.parametrize(
    ("mu", "tau", "expected"),
    [
        (0.0, 0.0, (1.0, 1.0, 0.05, 0.05)),
        (0.8, 0.0, (0.2, 1.8, 0.05, 0.05)),
        (0.5, 0.2, (0.5, 1.5, 0.04, 0.06)),
    ],
)
def test_rectangular_mu_tau_geometry_and_four_equal_plies(
    mu: float,
    tau: float,
    expected: tuple[float, float, float, float],
) -> None:
    geometry = target.geometry_for(mu, tau, beta_deg=15.0)
    assert_allclose(
        [geometry.L1, geometry.L2, geometry.h1, geometry.h2],
        expected,
        rtol=0.0,
        atol=2.0e-16,
    )
    assert geometry.b1 == geometry.b2 == 0.2
    assert_allclose(geometry.L1 + geometry.L2, 2.0, rtol=0.0, atol=2.0e-16)

    objects = target.build_model_objects(geometry)
    checks = target.constitutive_checks(geometry, objects)
    assert checks["status"] == "PASS"
    assert checks["exact_case_A_material"] is True
    for arm, thickness in (
        (objects.arm1, geometry.h1),
        (objects.arm2, geometry.h2),
    ):
        material = arm.weak_laminate.plies[0].material
        assert (
            material.E1,
            material.E2,
            material.nu12,
            material.G12,
            material.G13,
            material.G23,
            material.rho,
        ) == (1.1, 0.9, 0.3, 1.0 / 2.6, 1.0 / 2.6, 1.0 / 2.6, 1.0)
        assert len(arm.weak_laminate.plies) == 4
        assert [ply.angle_deg for ply in arm.weak_laminate.plies] == [
            0.0,
            90.0,
            90.0,
            0.0,
        ]
        assert_allclose(
            [ply.thickness for ply in arm.weak_laminate.plies],
            np.full(4, thickness / 4.0),
            rtol=0.0,
            atol=0.0,
        )


def test_physical_eb_adapter_matches_frozen_tau0_determinant_identity() -> None:
    diagnostic = target.eb_tau0_equivalence_diagnostic()
    assert diagnostic["status"] == "PASS"
    assert diagnostic["maximum_relative_residual"] <= diagnostic["tolerance"]
    assert len(diagnostic["probes"]) >= 3
    assert {float(row["mu"]) for row in diagnostic["probes"]} >= {
        0.0,
        0.4,
        0.8,
    }


def test_analysis_entrypoint_has_no_forbidden_project_imports() -> None:
    tree = ast.parse(SCRIPT_PATH.read_text(encoding="utf-8"))
    project_imports: set[str] = set()
    dynamic_imports: list[ast.Call] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            project_imports.update(
                alias.name
                for alias in node.names
                if alias.name.startswith(("scripts.", "my_project."))
            )
        elif isinstance(node, ast.ImportFrom) and node.module:
            if node.module.startswith(("scripts.", "my_project.")):
                project_imports.add(node.module)
                project_imports.update(
                    f"{node.module}.{alias.name}" for alias in node.names
                )
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name) and node.func.id == "__import__":
                dynamic_imports.append(node)
            if isinstance(node.func, ast.Attribute) and node.func.attr == "import_module":
                dynamic_imports.append(node)

    imported_text = " ".join(sorted(project_imports)).lower()
    forbidden = (
        "ritz",
        "fem",
        "yartsev",
        "circular",
        "shellbuckling",
        "shell_buckling",
        "branch_tracking",
        "layer_order",
    )
    assert not any(fragment in imported_text for fragment in forbidden)
    assert project_imports == {
        "scripts.analysis.laminated_beams",
        "scripts.analysis.laminated_beams.compare_rectangular_weakly_orthotropic_models_vs_beta",
    }
    assert dynamic_imports == []


def test_closing_atomic_csv_and_complete_group_contract(tmp_path: Path) -> None:
    rows = [
        {
            "sweep": target.SWEEP_MU,
            "model": target.MODEL_RLB,
            "mu": 0.74,
            "sorted_position": position,
            "guard_flag": position == 9,
        }
        for position in range(1, 10)
    ]
    assert target._complete_mu_values(rows, target.MODEL_RLB) == [0.74]
    path = tmp_path / "roots.csv"
    target._atomic_write_csv(path, rows)
    assert path.is_file()
    assert len(_read_csv(path)) == 9
    assert not (tmp_path / ".roots.csv.closing.tmp").exists()

    with pytest.raises(RuntimeError, match="Duplicate closing key"):
        target._complete_mu_values([*rows, dict(rows[-1])], target.MODEL_RLB)
    with pytest.raises(RuntimeError, match="Incomplete"):
        target._complete_mu_values(rows[:-1], target.MODEL_RLB)


def test_closing_merge_ignores_only_identical_historical_overlap(
    tmp_path: Path,
) -> None:
    canonical = tmp_path / "canonical.csv"
    shard = tmp_path / "shard.csv"
    rows = [
        {
            "sweep": target.SWEEP_MU,
            "model": target.MODEL_RLB,
            "mu": 0.4,
            "sorted_position": position,
            "guard_flag": position == 9,
        }
        for position in range(1, 10)
    ]
    target._atomic_write_csv(canonical, rows)
    target._atomic_write_csv(shard, rows)

    merged, sources, ignored = target._merge_disjoint_or_identical_mu_rows(
        [canonical, shard], target.MODEL_RLB
    )
    assert len(merged) == 9
    assert len(sources) == 2
    assert ignored == 9

    conflicting = [dict(row) for row in rows]
    conflicting[0]["guard_flag"] = True
    target._atomic_write_csv(shard, conflicting)
    with pytest.raises(RuntimeError, match="Conflicting duplicate"):
        target._merge_disjoint_or_identical_mu_rows(
            [canonical, shard], target.MODEL_RLB
        )


def test_closing_thread_limit_contract_is_explicit() -> None:
    assert target.CLOSING_THREAD_LIMITS == {
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "NUMEXPR_NUM_THREADS": "1",
    }


def test_bounded_root9_policy_preserves_frozen_scan_spacings_and_gates() -> None:
    frozen = target.rlb2c.rlb2b.frozen_root_policy()
    bounded, record = target._bounded_root9_policy(frozen, 230.15893863973855)
    assert bounded.requested_roots == 8
    assert bounded.guard_roots == 1
    assert bounded.required_slots == 9
    assert bounded.Omega_max >= record["requested_Omega_max"]
    assert bounded.Omega_max < record["requested_Omega_max"] + 0.21
    assert_allclose(
        record["primary_scan_spacing"],
        record["frozen_primary_scan_spacing"],
        rtol=32.0 * np.finfo(float).eps,
        atol=0.0,
    )
    assert_allclose(
        record["verification_scan_spacing"],
        record["frozen_verification_scan_spacing"],
        rtol=32.0 * np.finfo(float).eps,
        atol=0.0,
    )
    for field in (
        "sigma_prefilter",
        "root_singular_ratio",
        "nullity_relative_threshold",
        "boundary_null_residual",
        "root_xtol_Omega",
        "root_rtol",
        "dedup_atol_Omega",
        "dedup_rtol",
        "cluster_atol_Omega",
        "cluster_rtol",
        "post_guard_tail_Omega",
        "local_close_pair_guard_subintervals",
        "primary_phases",
        "verification_phases",
    ):
        assert getattr(bounded, field) == getattr(frozen, field)


def _synthetic_attempt_inventory(
    *,
    agreement: bool,
    guard_at_boundary: bool = False,
    slot_count: int = 9,
    unresolved: int = 0,
    finite: bool = True,
) -> tuple[SimpleNamespace, dict[str, Any], Any]:
    policy = target.rlb2c.rlb2b.frozen_root_policy()

    def slot(position: int) -> SimpleNamespace:
        diagnostics = SimpleNamespace(
            finite=finite,
            scaled_sigma_ratio=1.0e-14,
            scaled_null_residual=1.0e-14,
            raw_boundary_null_residual=1.0e-14,
            root_gate_nullity=1,
            detected_nullity=1,
        )
        candidate = SimpleNamespace(diagnostics=diagnostics)
        event = SimpleNamespace(
            candidate=candidate,
            multiplicity=1,
            detected_nullity=1,
            cluster_id="",
            cluster_multiplicity=1,
            cluster_total_nullity=1,
        )
        return SimpleNamespace(sorted_slot=position, event=event)

    slots = [slot(position) for position in range(1, slot_count + 1)]
    inventory = SimpleNamespace(
        primary=SimpleNamespace(slots=slots),
        verification=SimpleNamespace(slots=slots),
        guard_available=slot_count >= 9,
        guard_not_at_scan_boundary=not guard_at_boundary,
    )
    export = {
        "unresolved_candidates_below_export_guard": unresolved,
        "export_primary_verification_agreement": agreement,
    }
    return inventory, export, policy


def test_agreement_only_failure_routes_to_local_adjudication_not_expansion() -> None:
    inventory, export, policy = _synthetic_attempt_inventory(agreement=False)
    assert (
        target._bounded_attempt_disposition(inventory, export, policy)
        == target.ATTEMPT_LOCAL_ADJUDICATION_REQUIRED
    )


@pytest.mark.parametrize(
    ("guard_at_boundary", "slot_count", "unresolved"),
    ((True, 9, 0), (False, 8, 0), (False, 9, 1)),
)
def test_only_completeness_truncation_routes_to_single_range_expansion(
    guard_at_boundary: bool, slot_count: int, unresolved: int
) -> None:
    inventory, export, policy = _synthetic_attempt_inventory(
        agreement=False,
        guard_at_boundary=guard_at_boundary,
        slot_count=slot_count,
        unresolved=unresolved,
    )
    assert (
        target._bounded_attempt_disposition(inventory, export, policy)
        == target.ATTEMPT_RANGE_EXPANSION_REQUIRED
    )


def test_complete_bad_physical_quality_is_hard_fail_not_range_expansion() -> None:
    inventory, export, policy = _synthetic_attempt_inventory(
        agreement=False, finite=False
    )
    assert (
        target._bounded_attempt_disposition(inventory, export, policy)
        == target.ATTEMPT_HARD_FAIL
    )


def test_adjudication_contract_keeps_initial_fail_and_writes_no_failed_rows() -> None:
    failed = {
        "initial_primary_verification_status": "FAIL",
        "local_adjudication_status": "FAIL",
        "final_point_status": "FAIL",
        "refined_root_diagnostics": [],
    }
    assert (
        target._adjudicated_point_rows(
            target.MODEL_RLB,
            target.geometry_for(0.8, target.MU_TAU, target.MU_BETA_DEG),
            [],
            {},
            failed,
        )
        == []
    )
    assert failed["initial_primary_verification_status"] == "FAIL"


def test_saved_adjudication_uses_persisted_and_reconstructed_brackets() -> None:
    source = SCRIPT_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source)
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "_adjudicate_reconstructed_pairs"
    )
    rendered = ast.unparse(function)
    assert "persisted['bracket_left_Omega']" in rendered
    assert "verification_candidate.interval_left_Omega" in rendered
    assert "LOCAL_ADJUDICATION_RELATIVE_TOLERANCE" in rendered
    assert "mean(" not in rendered
    assert "average(" not in rendered


def _synthetic_root9_guard_evidence() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], Any]:
    rows = []
    for position in range(1, 10):
        Omega = 10.0 * position
        rows.append(
            {
                "root_index": position,
                "Omega_primary": Omega,
                "Omega_verification": Omega * (1.0 + 1.0e-10),
                "relative_difference_Omega": 1.0e-10,
                "nearest_adjacent_gap_Omega": 10.0,
                "primary_original_bracket": [Omega - 0.4, Omega + 0.4],
                "verification_original_bracket": [Omega - 0.2, Omega + 0.2],
                "primary_scaled_sigma_ratio": 1.0e-14,
                "verification_scaled_sigma_ratio": 1.0e-14,
                "primary_boundary_residual": 1.0e-14,
                "verification_boundary_residual": 1.0e-14,
            }
        )

    def outcome(position: int, relative: float, passed: bool) -> dict[str, Any]:
        Omega = 10.0 * position
        values = [Omega, Omega + relative * Omega, Omega, Omega]
        path = {
            "local_sign_bracket": [Omega - 0.1, Omega + 0.1],
            "sign_bracket_count": 1,
            "status": "PASS",
            "determinant_converged": True,
            "sigma_converged": True,
            "determinant_quality_pass": True,
            "sigma_quality_pass": True,
            "determinant_scaled_sigma_ratio": 1.0e-14,
            "sigma_scaled_sigma_ratio": 1.0e-14,
            "determinant_boundary_residual": 1.0e-14,
            "sigma_boundary_residual": 1.0e-14,
        }
        return {
            "root_index": position,
            "local_adjudication_status": "PASS" if passed else "FAIL",
            "all_refined_values_Omega": values,
            "final_relative_agreement_Omega": relative,
            "nearest_adjacent_gap_Omega": 10.0,
            "ordering_preserved": True,
            "distinct_physical_candidate_between_estimates": False,
            "canonical_selection": "DETERMINISTIC_PRIMARY_DETERMINANT",
            "primary_path": dict(path),
            "verification_path": dict(path),
        }

    initial = {"root_diagnostics": rows}
    adjudication = {
        "refined_root_diagnostics": [
            outcome(8, 2.0e-9, True),
            outcome(9, 8.0e-8, False),
        ]
    }
    attempt = {
        "root_count": 9,
        "unresolved_candidates_below_root9": 0,
        "effective_Omega_max": 110.0,
        "guard_not_at_scan_boundary": True,
    }
    policy = SimpleNamespace(post_guard_tail_Omega=2.0)
    return initial, adjudication, attempt, policy


def test_isolated_root9_strict_fail_can_receive_guard_pass() -> None:
    evidence = _synthetic_root9_guard_evidence()
    result = target._root9_guard_contract(*evidence)
    assert result["strict_roots_1_to_8_status"] == "PASS"
    assert result["root9_strict_agreement_status"] == "FAIL"
    assert result["root9_guard_status"] == target.ROOT9_GUARD_INTERVAL_PASS
    assert result["point_status"] == target.POINT_PASS_WITH_GUARD_QUALIFICATION
    assert result["root_estimates_averaged"] is False


def test_root9_guard_never_exempts_a_strict_root_in_positions_1_to_8() -> None:
    initial, adjudication, attempt, policy = _synthetic_root9_guard_evidence()
    initial["root_diagnostics"][6]["relative_difference_Omega"] = 2.0e-8
    result = target._root9_guard_contract(
        initial, adjudication, attempt, policy
    )
    assert result["strict_by_position"]["7"] == "FAIL"
    assert result["strict_roots_1_to_8_status"] == "FAIL"
    assert result["root9_guard_status"] == "FAIL"


def test_intersecting_root8_root9_enclosures_fail_guard() -> None:
    initial, adjudication, attempt, policy = _synthetic_root9_guard_evidence()
    root9 = adjudication["refined_root_diagnostics"][1]
    for path in (root9["primary_path"], root9["verification_path"]):
        path["local_sign_bracket"] = [79.95, 80.05]
    root9["all_refined_values_Omega"] = [80.0, 80.01, 80.0, 80.01]
    result = target._root9_guard_contract(
        initial, adjudication, attempt, policy
    )
    assert result["enclosures_intersect"] is True
    assert result["root9_guard_status"] == "FAIL"


def test_ambiguous_root9_ordering_fails_guard() -> None:
    initial, adjudication, attempt, policy = _synthetic_root9_guard_evidence()
    adjudication["refined_root_diagnostics"][1]["ordering_preserved"] = False
    result = target._root9_guard_contract(
        initial, adjudication, attempt, policy
    )
    assert result["ordering_same_in_all_passes"] is False
    assert result["root9_guard_status"] == "FAIL"


def test_root9_initial_fail_is_retained_in_guard_diagnostics() -> None:
    result = target._root9_guard_contract(*_synthetic_root9_guard_evidence())
    assert result["root9_strict_agreement_status"] == "FAIL"
    assert result["root9_guard_status"] == target.ROOT9_GUARD_INTERVAL_PASS
    assert result["root9_plotted"] is False


def test_guard_sentinel_requires_a_strict_first8_certificate() -> None:
    rows = [
        {
            "sweep": target.SWEEP_MU,
            "model": target.MODEL_OLD,
            "mu": 0.8,
            "sorted_position": position,
            "root_status": "PASS",
            "export_range_status": target.POINT_PASS_WITH_GUARD_QUALIFICATION,
            "internal_inventory_status": "NOT_COMPUTED_ABOVE_ROOT9",
            "export_primary_verification_max_relative": 8.0e-8,
        }
        for position in range(1, 10)
    ]
    uncertified = target._point_qualification(rows)
    assert uncertified["PLOTTED_FIRST_8"] != "PASS"
    assert uncertified["ROOT_9_GUARD"] != "PASS"

    certificate = {
        "strict_roots_1_to_8_status": "PASS",
        "root9_strict_agreement_status": "FAIL",
        "root9_guard_status": target.ROOT9_GUARD_INTERVAL_PASS,
        "point_status": target.POINT_PASS_WITH_GUARD_QUALIFICATION,
    }
    certified = target._point_qualification(
        rows, root9_guard_certificate=certificate
    )
    assert certified["PLOTTED_FIRST_8"] == "PASS"
    assert certified["ROOT_9_GUARD"] == "PASS"

    certificate["strict_roots_1_to_8_status"] = "FAIL"
    rejected = target._point_qualification(
        rows, root9_guard_certificate=certificate
    )
    assert rejected["PLOTTED_FIRST_8"] != "PASS"
    assert rejected["ROOT_9_GUARD"] != "PASS"


def test_root9_prediction_uses_two_nearest_complete_points_only() -> None:
    rows = []
    for mu_value, guard in ((0.64, 165.2), (0.66, 165.0), (0.70, 170.0)):
        rows.extend(
            {
                "sweep": target.SWEEP_MU,
                "model": target.MODEL_OLD,
                "mu": mu_value,
                "sorted_position": position,
                "guard_flag": position == 9,
                "Omega": guard if position == 9 else position,
            }
            for position in range(1, 10)
        )
    prediction = target._predict_root9_Omega(rows, target.MODEL_OLD, 0.68)
    assert prediction["source_points"] == [
        {"mu": 0.66, "Omega_9": 165.0},
        {"mu": 0.7, "Omega_9": 170.0},
    ]
    assert_allclose(
        prediction["predicted_Omega_9"], 167.5, rtol=0.0, atol=2.0e-13
    )


def test_atomic_mu_merge_preserves_every_existing_field(tmp_path: Path) -> None:
    path = tmp_path / "roots.csv"
    existing = [
        {
            "sweep": target.SWEEP_MU,
            "model": target.MODEL_OLD,
            "mu": 0.0,
            "sorted_position": position,
            "guard_flag": position == 9,
            "Omega": 10.0 + position,
        }
        for position in range(1, 10)
    ]
    target._atomic_write_csv(path, existing)
    existing_from_csv = _read_csv(path)
    new_rows = [
        {
            **{field: "" for field in existing_from_csv[0]},
            "sweep": target.SWEEP_MU,
            "model": target.MODEL_OLD,
            "mu": "0.8",
            "sorted_position": str(position),
            "guard_flag": str(position == 9),
            "Omega": str(20.0 + position),
        }
        for position in range(1, 10)
    ]
    target._atomic_merge_mu_csv_preserving_existing(
        path, existing_from_csv, [*existing_from_csv, *new_rows]
    )
    merged = _read_csv(path)
    assert merged[:9] == existing_from_csv
    assert len(merged) == 18
    assert not (tmp_path / ".roots.csv.production.tmp").exists()


def _synthetic_plot_rows(
    sweep: str,
) -> dict[str, list[dict[str, float | int]]]:
    axis_name = "mu" if sweep == target.SWEEP_MU else "beta_deg"
    axis_values = (0.0, 0.4, 0.8) if sweep == target.SWEEP_MU else (0.0, 45.0, 90.0)
    return {
        model: [
            {
                axis_name: axis_value,
                "sorted_position": position,
                "Lambda": 3.0 * position + axis_value / 90.0 + 0.01 * model_index,
            }
            for position in range(1, target.OUTPUT_GUARD_POSITION + 1)
            for axis_value in axis_values
        ]
        for model_index, model in enumerate(target.MODELS)
    }


@pytest.mark.parametrize("sweep", target.SWEEPS)
def test_create_plot_has_24_curves_requested_styles_and_one_png(
    sweep: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    original_close = target.plt.close
    monkeypatch.setattr(target.plt, "close", lambda *_args, **_kwargs: None)
    axis_field = "mu" if sweep == target.SWEEP_MU else "beta_deg"
    x_label = "mu" if sweep == target.SWEEP_MU else "beta"
    x_limits = (0.0, 0.8) if sweep == target.SWEEP_MU else (0.0, 90.0)
    output = target.create_plot(
        _synthetic_plot_rows(sweep),
        tmp_path / target.PLOT_FILENAMES[sweep],
        axis_field=axis_field,
        x_label=x_label,
        x_limits=x_limits,
    )
    fig = target.plt.gcf()
    ax = fig.axes[0]
    try:
        assert output.is_file() and output.stat().st_size > 0
        assert len(ax.lines) == 3 * target.PLOTTED_POSITIONS == 24
        for position_zero in range(target.PLOTTED_POSITIONS):
            group = ax.lines[3 * position_zero : 3 * position_zero + 3]
            assert [line.get_linestyle() for line in group] == [
                target.LINE_STYLES[model] for model in target.MODELS
            ]
            assert len({tuple(line.get_color()) for line in group}) == 1
        assert len(
            {tuple(ax.lines[3 * index].get_color()) for index in range(8)}
        ) == 8
        assert len(
            [child for child in ax.get_children() if isinstance(child, Legend)]
        ) == 2
        assert ax.get_ylabel() == "Lambda"
        assert ax.get_xlabel() == x_label
        graphics = [
            path
            for path in tmp_path.iterdir()
            if path.suffix.lower() in {".png", ".pdf", ".svg"}
        ]
        assert graphics == [output]
    finally:
        original_close(fig)


def _assert_representative_inventory(
    rows: list[dict[str, str]],
    *,
    axis_name: str,
    representative_values: tuple[float, ...],
) -> None:
    for axis_value in representative_values:
        selected = [
            row
            for row in rows
            if math.isclose(
                float(row[axis_name]), axis_value, rel_tol=0.0, abs_tol=5.0e-13
            )
        ]
        assert len(selected) == target.OUTPUT_GUARD_POSITION
        selected.sort(key=lambda row: int(row["sorted_position"]))
        assert [int(row["sorted_position"]) for row in selected] == list(range(1, 10))
        omega = np.asarray([float(row["omega"]) for row in selected])
        Omega = np.asarray([float(row["Omega"]) for row in selected])
        Lambda = np.asarray([float(row["Lambda"]) for row in selected])
        assert np.all(omega > 0.0)
        assert np.all(np.diff(omega) >= 0.0)
        assert_allclose(
            Omega,
            omega * target.reference_omega_scale(),
            rtol=target.NORMALIZATION_IDENTITY_TOLERANCE,
            atol=0.0,
        )
        assert_allclose(
            Lambda**2,
            Omega,
            rtol=target.NORMALIZATION_IDENTITY_TOLERANCE,
            atol=0.0,
        )
        assert [row["guard_flag"].lower() == "true" for row in selected] == [
            False
        ] * 8 + [True]


def test_existing_generated_outputs_have_representative_roots_and_two_pngs() -> None:
    run_manifest_path = RESULT_DIR / "run_manifest.json"
    if not run_manifest_path.is_file():
        pytest.skip("ignored full RLB-2D output has not been finalized")

    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    case_contract = json.loads(
        (RESULT_DIR / "case_contract.json").read_text(encoding="utf-8")
    )
    assert case_contract["sweeps"][target.SWEEP_MU]["values"] == [
        float(value) for value in target.mu_grid()
    ]
    assert case_contract["sweeps"][target.SWEEP_BETA]["values"] == [
        float(value) for value in target.beta_grid()
    ]
    assert case_contract["inter_model_relative_differences_computed"] is False
    assert case_contract["sweeps"][target.SWEEP_MU]["beta_deg"] == 15.0
    assert case_contract["sweeps"][target.SWEEP_MU]["tau"] == 0.0
    assert case_contract["sweeps"][target.SWEEP_BETA]["mu"] == 0.5
    assert case_contract["sweeps"][target.SWEEP_BETA]["tau"] == 0.2
    assert case_contract["new_RLB_lamina"]["case_id"] == "A"
    assert case_contract["new_RLB_lamina"]["delta"] == 0.1

    for sweep, expected_points, axis_name, representatives in (
        (target.SWEEP_MU, 41, "mu", (0.0, 0.4, 0.8)),
        (target.SWEEP_BETA, 91, "beta_deg", (0.0, 15.0, 90.0)),
    ):
        for model in target.MODELS:
            rows = _read_csv(RESULT_DIR / target.ROOT_FILENAMES[(sweep, model)])
            assert len(rows) == expected_points * target.OUTPUT_GUARD_POSITION
            assert sum(row["guard_flag"].lower() == "true" for row in rows) == expected_points
            expected_keys = {
                (round(float(value), 12), position)
                for value in (
                    target.mu_grid()
                    if sweep == target.SWEEP_MU
                    else target.beta_grid()
                )
                for position in range(1, 10)
            }
            assert {
                (round(float(row[axis_name]), 12), int(row["sorted_position"]))
                for row in rows
            } == expected_keys
            summary = target.validate_root_rows(
                sweep,
                model,
                rows,
                target.mu_grid()
                if sweep == target.SWEEP_MU
                else target.beta_grid(),
            )
            assert summary["exact_row_structure_passed"] is True
            assert summary["accepted_root_quality_passed"] is True
            assert summary["no_unresolved_below_export_guard_passed"] is True
            _assert_representative_inventory(
                rows,
                axis_name=axis_name,
                representative_values=representatives,
            )

    assert len(_read_csv(RESULT_DIR / "geometry_sanity_checks.csv")) == 3
    assert len(_read_csv(RESULT_DIR / "laminate_properties_summary.csv")) == 6
    geometry_rows = _read_csv(RESULT_DIR / "geometry_sanity_checks.csv")
    beta_geometry = next(
        row for row in geometry_rows if row["case_id"] == "mu_0p5_tau_0p2"
    )
    assert (float(beta_geometry["h1"]), float(beta_geometry["h2"])) == (
        0.04,
        0.06,
    )
    graphics = sorted(
        path.name
        for path in RESULT_DIR.iterdir()
        if path.suffix.lower() in {".png", ".pdf", ".svg"}
    )
    assert graphics == sorted(target.PLOT_FILENAMES.values())
    assert not list(RESULT_DIR.glob("*comparison*.csv"))

    for name, expected_hash in run_manifest["generated_file_hashes"].items():
        assert _sha256(RESULT_DIR / name) == expected_hash
    closing = run_manifest["closing_stage"]
    assert closing["new_points_executed"] == sum(
        len(values) for values in closing["missing_mu_before"].values()
    ) == 9
    assert closing["reused_points"] == 387
    assert closing["ready_points_recalculated"] == 0
    assert closing["parallel_workers_used_in_closing_stage"] == 0
    assert closing["inherited_workers_drained"] == 0
    assert closing["global_restarts"] == 0
    assert closing["thread_limits"] == target.CLOSING_THREAD_LIMITS
    assert closing["missing_mu_after"] == {
        target.MODEL_OLD: [],
        target.MODEL_RLB: [],
    }
    assert closing["preexisting_groups_preserved"] is True
    assert closing["roots_above_9_requested_or_exported"] is False
    assert closing["FAST_LOCAL_point_count"] == 0
    assert closing["beta_plot_reused_without_redraw"] is True
    assert closing["beta_plot_sha256_before"] == closing["beta_plot_sha256_after"]
    assert run_manifest["plotted_sorted_positions"] == 8
    assert run_manifest["output_guard_position"] == 9
    assert run_manifest["root9_plotted"] is False
    assert run_manifest["figures_created_in_closing_stage"] == 1
    assert run_manifest["figures_available"] == 2
    assert isinstance(
        run_manifest["exact_qualifications_affecting_plotted_range"], list
    )
    assert isinstance(
        run_manifest["exact_qualifications_only_above_root9"], list
    )
    assert set(run_manifest["statuses"]) == {
        "RLB-2D-BETA-PLOTTED-FIRST-8",
        "RLB-2D-MU-PLOTTED-FIRST-8",
        "RLB-2D-ROOT9-GUARDS",
        "RLB-2D-INTERNAL-TAIL",
        "RLB-2D-PLOT-GENERATION",
        "SCIENTIFIC_OVERALL",
    }
    assert run_manifest["frozen_models_preserved"] is True
    assert run_manifest["inter_model_relative_differences_computed"] is False
    assert run_manifest["Ritz_run"] is False
    assert run_manifest["FEM_run"] is False
    assert run_manifest["branch_tracking_run"] is False
    assert run_manifest["commit_performed"] is False
    assert run_manifest["push_performed"] is False
    report_text = (RESULT_DIR / "report.md").read_text(encoding="utf-8")
    for name, value in run_manifest["statuses"].items():
        assert f"- {name}: {value}" in report_text
