"""Targeted regression gates for the finite RLB-2A layer-order pilot.

The physical coupled searches are not repeated here.  Their frozen roots are
evaluated directly in the corresponding boundary matrices, while the ignored
generated manifest supplies the local closing audit when it is present.
"""

from __future__ import annotations

import ast
import csv
import json
import math
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose
import pytest

from scripts.analysis.laminated_beams import (
    compare_reddy_cross_ply_layer_order as pilot,
)


ROOT = Path(__file__).resolve().parents[1]
RESULT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_cross_ply_layer_order_pilot"
)

FROZEN_BETA0_ROOTS = {
    "OUTER_0": (
        22.67206118744846,
        50.61286710541817,
        83.29900985509377,
        117.87342669122222,
        153.2219871010395,
        188.6755267183602,
        224.03164905362817,
        226.78523261343506,
        259.2125392941688,
        294.2274722772023,
        329.0885957451078,
        363.8240553255954,
        398.4504559370047,
    ),
    "OUTER_90": (
        12.086284916899912,
        31.004705692928237,
        56.108284534676876,
        85.27410108593436,
        117.12847018649592,
        150.73099154053125,
        185.44728810381324,
        220.84145554832565,
        226.78523261343506,
        256.61618894583944,
        292.5730154522818,
        328.58384626009246,
        364.56922367340457,
    ),
}

FROZEN_BETA30_ROOTS = {
    "OUTER_0": (
        47.63362821858805,
        50.584682507818,
        95.31651389758419,
        117.6972865246834,
        158.44954384296773,
        187.97753248692615,
        224.21358584438644,
        226.17072474234803,
        260.1009930501306,
        289.5226157525381,
        329.33834863440234,
        352.6780780491725,
        398.53975395141117,
    ),
    "OUTER_90": (
        30.998487270561117,
        36.845268377173184,
        75.33518275489152,
        85.22466279953893,
        126.07878889016139,
        150.539821745891,
        188.19811757981418,
        218.28075731217518,
        228.71526158011986,
        254.69955152130456,
        292.90265053105895,
        321.3782615698171,
        364.6961781154533,
    ),
}

FROZEN_BETA30_HASHES = {
    "OUTER_0": "41985174A53C1B120D17C52BEFB691D3D32C54B4CFF2CD68EC20673E7C72F985",
    "OUTER_90": "B69FE2997EC6216F57619867B3CBD825F8F7A1E340169AB80B20C996A68AFFA4",
}


@pytest.fixture(scope="module")
def cases() -> dict[str, pilot.LayerOrderCase]:
    return pilot.build_layer_order_cases()


@pytest.fixture(scope="module")
def direct_families(
    cases: dict[str, pilot.LayerOrderCase],
) -> dict[str, pilot.DirectFamilyInventory]:
    policy = pilot.beta0_pilot.SearchPolicy()
    return {
        stack_id: pilot._direct_family_inventory(case, policy)
        for stack_id, case in cases.items()
    }


def _relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), np.finfo(float).tiny)


def test_analysis_script_imports_only_the_frozen_rlb_workflow() -> None:
    path = (
        ROOT
        / "scripts"
        / "analysis"
        / "laminated_beams"
        / "compare_reddy_cross_ply_layer_order.py"
    )
    tree = ast.parse(path.read_text(encoding="utf-8"))
    project_imports: set[str] = set()
    imported_names: set[str] = set()
    dynamic_imports: list[ast.Call] = []
    locally_defined_functions: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("scripts."):
                    project_imports.add(alias.name)
                    imported_names.add(alias.name)
        elif isinstance(node, ast.ImportFrom) and node.module:
            if node.module.startswith("scripts."):
                project_imports.add(node.module)
                imported_names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.FunctionDef):
            locally_defined_functions.add(node.name)
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name) and node.func.id == "__import__":
                dynamic_imports.append(node)
            if isinstance(node.func, ast.Attribute) and node.func.attr == "import_module":
                dynamic_imports.append(node)

    assert project_imports == {
        "scripts.analysis.laminated_beams",
        "scripts.lib",
    }
    assert imported_names == {
        "pilot_reddy_symmetric_coupled_beams_beta0",
        "reddy_symmetric_coupled_beams",
        "reddy_symmetric_laminated_beam",
    }
    forbidden = (
        "isotropic_rectangular",
        "circular",
        "ritz",
        "yartsev",
        "shellbuckling",
        "euler",
        "fem",
    )
    imported_text = " ".join(sorted(imported_names)).lower()
    assert not any(fragment in imported_text for fragment in forbidden)
    assert dynamic_imports == []
    assert "seed_free_root_inventory" not in locally_defined_functions
    assert "find_bending_roots" not in locally_defined_functions


def test_frozen_material_geometry_K_and_frequency_scale_are_reused(
    cases: dict[str, pilot.LayerOrderCase],
) -> None:
    outer_0 = cases["OUTER_0"]
    material = outer_0.laminate.plies[0].material
    assert material.E1 == 25.0
    assert material.E2 == 1.0
    assert material.G12 == 0.5
    assert material.G13 == 0.5
    assert material.G23 == 0.2
    assert material.nu12 == 0.25
    assert material.rho == 1.0
    assert outer_0.properties.width == 1.0
    assert outer_0.laminate.thickness == 1.0
    assert outer_0.total_length == 20.0
    assert outer_0.length_1 == outer_0.length_2 == 10.0
    assert outer_0.properties.K == pytest.approx(5.0 / 6.0, rel=0.0, abs=0.0)
    assert outer_0.frequency_scale == 400.0
    assert cases["OUTER_90"].frequency_scale == 400.0


def test_exact_four_equal_plies_orientations_and_interfaces(
    cases: dict[str, pilot.LayerOrderCase],
) -> None:
    expected = {
        "OUTER_0": (0.0, 90.0, 90.0, 0.0),
        "OUTER_90": (90.0, 0.0, 0.0, 90.0),
    }
    expected_interfaces = np.array([-0.5, -0.25, 0.0, 0.25, 0.5])
    for stack_id, case in cases.items():
        assert len(case.laminate.plies) == 4
        assert tuple(ply.angle_deg for ply in case.laminate.plies) == expected[stack_id]
        assert_allclose(
            [ply.thickness for ply in case.laminate.plies],
            np.full(4, 0.25),
            rtol=0.0,
            atol=0.0,
        )
        assert_allclose(
            case.laminate.z_interfaces,
            expected_interfaces,
            rtol=0.0,
            atol=0.0,
        )
        assert case.stacking_sequence[0] == case.stacking_sequence[-1]
        assert case.stacking_sequence[1] == case.stacking_sequence[-2]


def test_full_constitutive_analytic_gate_passes(
    cases: dict[str, pilot.LayerOrderCase],
) -> None:
    rows, summary = pilot.constitutive_gate(cases)
    assert len(rows) == 24
    assert summary["pass"] is True
    assert all(row["status"] == "PASS" for row in rows)
    assert summary["maximum_A_formula_or_equality_relative"] <= 1.0e-12
    assert summary["maximum_D_formula_relative"] <= 1.0e-12
    assert summary["maximum_scaled_B_residual"] <= 1.0e-12
    assert summary["Dbeam_ratio_OUTER_0_over_OUTER_90"] == pytest.approx(
        5.5, rel=3.0e-16
    )


def test_A_shear_mass_and_reduced_properties_are_equal_but_D_differs(
    cases: dict[str, pilot.LayerOrderCase],
) -> None:
    outer_0 = cases["OUTER_0"]
    outer_90 = cases["OUTER_90"]
    assert_allclose(outer_0.laminate.A, outer_90.laminate.A, rtol=1.0e-12, atol=0.0)
    assert_allclose(
        outer_0.laminate.shear,
        outer_90.laminate.shear,
        rtol=1.0e-12,
        atol=0.0,
    )
    assert outer_0.laminate.I0 == outer_90.laminate.I0
    assert outer_0.laminate.I1 == outer_90.laminate.I1 == 0.0
    assert outer_0.laminate.I2 == outer_90.laminate.I2
    assert_allclose(
        outer_0.laminate.B,
        np.zeros((3, 3)),
        rtol=0.0,
        atol=1.0e-16,
    )
    assert_allclose(
        outer_90.laminate.B,
        np.zeros((3, 3)),
        rtol=0.0,
        atol=1.0e-16,
    )
    assert not np.allclose(outer_0.laminate.D, outer_90.laminate.D)
    assert outer_0.properties.A == outer_90.properties.A
    assert outer_0.properties.S == outer_90.properties.S
    assert outer_0.properties.m == outer_90.properties.m
    assert outer_0.properties.J == outer_90.properties.J
    assert outer_0.properties.D > outer_90.properties.D


def test_beta0_axial_roots_are_identical(
    cases: dict[str, pilot.LayerOrderCase],
) -> None:
    roots: dict[str, tuple[float, ...]] = {}
    for stack_id, case in cases.items():
        modes = pilot.rlb.exact_axial_modes(
            case.properties, case.total_length, "FF", n_modes=24
        )
        roots[stack_id] = tuple(mode.omega * case.frequency_scale for mode in modes)
    assert_allclose(
        roots["OUTER_0"],
        roots["OUTER_90"],
        rtol=1.0e-10,
        atol=0.0,
    )


def test_beta0_bending_roots_increase_for_outer_zero(
    cases: dict[str, pilot.LayerOrderCase],
    direct_families: dict[str, pilot.DirectFamilyInventory],
) -> None:
    roots_0 = np.array(
        [
            root.omega * cases["OUTER_0"].frequency_scale
            for root in direct_families["OUTER_0"].bending
        ]
    )
    roots_90 = np.array(
        [
            root.omega * cases["OUTER_90"].frequency_scale
            for root in direct_families["OUTER_90"].bending
        ]
    )
    assert roots_0.size == roots_90.size == 24
    allowance = 1.0e-10 * np.maximum.reduce(
        (np.abs(roots_0), np.abs(roots_90), np.ones_like(roots_0))
    )
    assert np.all(roots_0 + allowance > roots_90)
    assert float(np.min((roots_0 - roots_90) / roots_90)) > 0.0


def test_direct_family_through_guard_flags_use_one_based_output_indices(
    cases: dict[str, pilot.LayerOrderCase],
    direct_families: dict[str, pilot.DirectFamilyInventory],
) -> None:
    for stack_id in pilot.STACK_IDS:
        axial_rows = pilot._axial_rows(cases[stack_id], direct_families[stack_id])
        bending_rows = pilot._bending_rows(
            cases[stack_id], direct_families[stack_id]
        )
        represented_axial = [
            int(row["family_index"])
            for row in axial_rows
            if row["represented_through_union_guard"]
        ]
        represented_bending = [
            int(row["family_index"])
            for row in bending_rows
            if row["represented_through_union_guard"]
        ]
        assert represented_axial == [1]
        assert represented_bending == list(range(1, 13))


@pytest.mark.parametrize("stack_id", ["OUTER_0", "OUTER_90"])
def test_beta0_coupled_frozen_roots_match_the_family_union_and_quality_gates(
    stack_id: str,
    cases: dict[str, pilot.LayerOrderCase],
    direct_families: dict[str, pilot.DirectFamilyInventory],
) -> None:
    case = cases[stack_id]
    union = direct_families[stack_id].union_rows
    roots = FROZEN_BETA0_ROOTS[stack_id]
    provider = pilot.beta0_pilot._coupled_provider(pilot.coupled, pilot._case_spec(case))
    assert len(union) == len(roots) == 13
    for row, omega_bar in zip(union, roots, strict=True):
        assert _relative_difference(float(row["omega_bar"]), omega_bar) <= 1.0e-9
        diagnostic = pilot.beta0_pilot.boundary_matrix_diagnostics(
            omega_bar, provider, case.frequency_scale
        )
        assert diagnostic.detected_nullity == 1
        assert diagnostic.scaled_sigma_ratio <= 1.0e-9
        assert diagnostic.raw_boundary_null_residual <= 1.0e-9


def test_beta30_first_twelve_plus_guard_are_ordered_monotonically() -> None:
    roots_0 = np.asarray(FROZEN_BETA30_ROOTS["OUTER_0"])
    roots_90 = np.asarray(FROZEN_BETA30_ROOTS["OUTER_90"])
    assert roots_0.size == roots_90.size == 13
    assert np.all(np.diff(roots_0) > 0.0)
    assert np.all(np.diff(roots_90) > 0.0)
    scale = np.maximum.reduce(
        (np.abs(roots_0), np.abs(roots_90), np.ones_like(roots_0))
    )
    assert np.all(roots_0 + 1.0e-10 * scale >= roots_90)
    assert np.all((roots_0 - roots_90) / roots_90 > 0.0)


@pytest.mark.parametrize("stack_id", ["OUTER_0", "OUTER_90"])
def test_beta30_frozen_roots_pass_physical_boundary_quality(
    stack_id: str,
    cases: dict[str, pilot.LayerOrderCase],
) -> None:
    case = cases[stack_id]
    provider = pilot._beta30_provider(case)
    roots = FROZEN_BETA30_ROOTS[stack_id]
    assert len(roots) == 13
    for omega_bar in roots:
        diagnostic = pilot.beta0_pilot.boundary_matrix_diagnostics(
            omega_bar, provider, case.frequency_scale
        )
        assert diagnostic.detected_nullity == 1
        assert diagnostic.scaled_sigma_ratio <= 1.0e-9
        assert diagnostic.raw_boundary_null_residual <= 1.0e-9


def test_generated_closing_audit_has_no_unresolved_candidates_below_guard() -> None:
    manifest_path = RESULT_DIR / "run_manifest.json"
    if not manifest_path.exists():
        pytest.skip("ignored generated RLB-2A closing output is not present")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["beta30_inventory_hashes_frozen_before_comparison"] == (
        FROZEN_BETA30_HASHES
    )
    assert manifest["one_stack_roots_used_as_other_stack_seeds"] is False
    assert manifest["preservation_pass"] is True
    assert manifest["result_phrase"] == "CONTROLLED_LAYER_ORDER_PILOT_COMPLETED"
    assert set(manifest["statuses"].values()) == {"PASS"}
    quality = manifest["inventory_quality"]
    assert len(quality) == 4
    for item in quality:
        assert item["slot_count"] == 13
        assert item["primary_slot_count"] == 13
        assert item["verification_slot_count"] == 13
        assert item["unresolved_candidates_below_guard"] == 0
        assert item["guard_available"] is True
        assert item["guard_not_at_scan_boundary"] is True
        assert item["maximum_scaled_sigma_ratio"] <= 1.0e-9
        assert item["maximum_boundary_null_residual"] <= 1.0e-9
        assert item["pass"] is True


def test_generated_output_contract_and_row_counts() -> None:
    if not RESULT_DIR.exists():
        pytest.skip("ignored generated RLB-2A output is not present")
    expected_csv_rows = {
        "laminate_properties.csv": 2,
        "analytic_stiffness_checks.csv": 24,
        "beta0_axial_roots.csv": 48,
        "beta0_bending_roots.csv": 48,
        "beta0_coupled_roots.csv": 26,
        "beta0_family_comparison.csv": 74,
        "beta30_roots.csv": 26,
        "beta30_spectrum_shift.csv": 13,
    }
    for filename, expected in expected_csv_rows.items():
        with (RESULT_DIR / filename).open(encoding="utf-8", newline="") as stream:
            assert len(list(csv.DictReader(stream))) == expected
    candidate_case_ids: set[str] = set()
    for filename in ("root_candidates.csv", "rejected_candidates.csv"):
        with (RESULT_DIR / filename).open(encoding="utf-8", newline="") as stream:
            candidate_case_ids.update(row["case_id"] for row in csv.DictReader(stream))
    for filename in ("beta0_coupled_roots.csv", "beta30_roots.csv"):
        with (RESULT_DIR / filename).open(encoding="utf-8", newline="") as stream:
            root_rows = list(csv.DictReader(stream))
        assert all(row["case_id"] in candidate_case_ids for row in root_rows)
    assert not list(RESULT_DIR.glob("*.png"))
    assert not list(RESULT_DIR.glob("*.pdf"))
