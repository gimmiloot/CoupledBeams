"""Fast contract tests for the finite rectangular RLB-2B beta workflow.

The tests below never run the 91-point three-model sweep.  They exercise the
shared dimensional contract, two isolated EB inventories, frozen ISO-01 and
ISO-03 root evidence through direct boundary residuals, and the pure plotting
path.  The ignored generated output is audited only when it already exists.
"""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from matplotlib.legend import Legend
import matplotlib.pyplot as plt
import numpy as np
from numpy.testing import assert_allclose
import pytest

from scripts.analysis.laminated_beams import (
    compare_rectangular_isotropic_models_vs_beta as target,
)


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = (
    ROOT
    / "scripts"
    / "analysis"
    / "laminated_beams"
    / "compare_rectangular_isotropic_models_vs_beta.py"
)
RESULT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "rectangular_isotropic_beta_sweep_comparison"
)

FROZEN_EB_COMMON_LAMBDA = {
    0.0: np.asarray(
        [
            5.593321362015330,
            15.418205716980063,
            30.225847931780976,
            49.964862031800190,
            74.638883824545790,
            104.24769645885272,
            108.82796185405304,
            138.79131189164934,
            178.26972949431268,
        ],
        dtype=float,
    ),
    30.0: np.asarray(
        [
            15.414786047534621,
            19.465013608304375,
            41.244114687393800,
            49.919175647938680,
            77.403659707073830,
            102.91681248389719,
            109.81555610725472,
            136.83812416891540,
            178.38685182943000,
        ],
        dtype=float,
    ),
}

# First eight sorted frequencies plus guard 9 from the committed ISO-01 and
# ISO-03 evidence.  They are used only as bounded residual probes; no scan is
# performed and no root from one model is used to locate a root of the other.
FROZEN_ISO_COMMON_LAMBDA = {
    0.0: {
        target.MODEL_OLD: np.asarray(
            [
                5.569209524503223,
                15.265747505349882,
                29.704528495665777,
                48.647470790117080,
                71.873403276132490,
                99.128226738179170,
                108.82796185405307,
                130.13966282522130,
                164.62650975773600,
            ],
            dtype=float,
        ),
        target.MODEL_RLB: np.asarray(
            [
                5.569209524502683,
                15.265747505349847,
                29.704528495663720,
                48.647470790116960,
                71.873403276051600,
                99.128226738227060,
                108.82796185405304,
                130.13966282486214,
                164.62650975800102,
            ],
            dtype=float,
        ),
    },
    30.0: {
        target.MODEL_OLD: np.asarray(
            [
                15.262445704715402,
                19.224364706350773,
                40.639214325859115,
                48.606172133489030,
                74.863656673453630,
                98.472436174017180,
                109.10012643497338,
                128.79413103922084,
                164.77216353285718,
            ],
            dtype=float,
        ),
        target.MODEL_RLB: np.asarray(
            [
                15.262445704715411,
                19.224364706349820,
                40.639214325867510,
                48.606172133489600,
                74.863656673385580,
                98.472436174004120,
                109.10012643499483,
                128.79413103921087,
                164.77216353301307,
            ],
            dtype=float,
        ),
    },
}


@pytest.fixture(scope="module")
def contract() -> dict[str, Any]:
    return target.build_case_contract([0.0, 30.0])


@pytest.fixture(scope="module")
def model_objects(contract: dict[str, Any]) -> target.ModelObjects:
    return target.build_model_objects(contract)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return [dict(row) for row in csv.DictReader(stream)]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest().upper()


def test_common_lambda_formula_and_historical_eb_square_root_map(
    contract: dict[str, Any],
    model_objects: target.ModelObjects,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    material = contract["material"]
    geometry = contract["geometry"]
    E = float(material["E"])
    rho = float(material["rho"])
    width = float(geometry["width_b"])
    thickness = float(geometry["thickness_h"])
    arm_length = float(geometry["l"])
    area = width * thickness
    inertia = width * thickness**3 / 12.0
    expected_scale = arm_length**2 * math.sqrt(rho * area / (E * inertia))
    rectangular_scale = arm_length**2 * math.sqrt(
        12.0 * rho / (E * thickness**2)
    )
    assert_allclose(
        target.lambda_scale(
            E=E,
            rho=rho,
            width=width,
            thickness=thickness,
            arm_length=arm_length,
        ),
        expected_scale,
        rtol=2.0e-16,
        atol=0.0,
    )
    assert_allclose(expected_scale, rectangular_scale, rtol=4.0e-16, atol=0.0)
    omega = 0.731
    common_value = target.common_lambda(
        omega,
        E=E,
        rho=rho,
        width=width,
        thickness=thickness,
        arm_length=arm_length,
    )
    assert_allclose(common_value, omega * expected_scale, rtol=2.0e-16, atol=0.0)

    captured: dict[str, float] = {}

    def fake_eb_matrix(
        historical_wavenumber: float, beta_rad: float, mu: float, epsilon: float
    ) -> np.ndarray:
        captured.update(
            wavenumber=float(historical_wavenumber),
            beta_rad=float(beta_rad),
            mu=float(mu),
            epsilon=float(epsilon),
        )
        return np.eye(6, dtype=float)

    monkeypatch.setattr(target, "assemble_eb_coupled_matrix", fake_eb_matrix)
    provider, metadata = target.make_matrix_provider(
        target.MODEL_EB, 30.0, contract, model_objects
    )
    assert_allclose(provider(omega), np.eye(6), rtol=0.0, atol=0.0)
    assert_allclose(captured["wavenumber"], math.sqrt(common_value), rtol=0.0, atol=2.0e-15)
    assert_allclose(captured["beta_rad"], math.radians(30.0), rtol=0.0, atol=0.0)
    assert captured["mu"] == 0.0
    assert_allclose(
        captured["epsilon"], thickness / (math.sqrt(12.0) * arm_length)
    )
    assert metadata["matrix_argument"] == "sqrt(Lambda)"


def test_canonical_contract_and_four_equal_isotropic_plies(
    contract: dict[str, Any], model_objects: target.ModelObjects
) -> None:
    assert contract["mu"] == 0.0
    assert contract["tau"] == 0.0
    assert contract["material"] == {
        "E": 1.0,
        "nu": 0.3,
        "rho": 1.0,
        "G": 1.0 / 2.6,
    }
    assert {
        key: contract["geometry"][key]
        for key in (
            "L1",
            "L2",
            "l",
            "L_total",
            "width_b",
            "thickness_h",
        )
    } == {
        "L1": 1.0,
        "L2": 1.0,
        "l": 1.0,
        "L_total": 2.0,
        "width_b": 0.2,
        "thickness_h": 0.05,
    }
    assert_allclose(contract["geometry"]["area_A"], 0.2 * 0.05, rtol=0.0, atol=0.0)
    assert_allclose(
        contract["geometry"]["second_moment_I_y"],
        0.2 * 0.05**3 / 12.0,
        rtol=0.0,
        atol=0.0,
    )
    assert contract["shear_correction"]["K"] == 5.0 / 6.0
    assert contract["new_RLB_laminate"]["number_of_plies"] == 4
    assert contract["new_RLB_laminate"]["one_ply_shortcut"] is False
    assert contract["new_RLB_laminate"]["stacking_sequence_deg"] == [
        0.0,
        90.0,
        90.0,
        0.0,
    ]
    assert len(model_objects.laminate.plies) == 4
    assert [ply.angle_deg for ply in model_objects.laminate.plies] == [
        0.0,
        90.0,
        90.0,
        0.0,
    ]
    assert_allclose(
        [ply.thickness for ply in model_objects.laminate.plies],
        np.full(4, 0.05 / 4.0),
        rtol=0.0,
        atol=0.0,
    )
    first_material = model_objects.laminate.plies[0].material
    assert first_material.E1 == first_material.E2 == 1.0
    assert first_material.nu12 == 0.3
    assert first_material.G12 == first_material.G13 == first_material.G23 == 1.0 / 2.6
    assert first_material.rho == 1.0


def test_old_timoshenko_and_four_ply_rlb_reduce_to_one_section_contract(
    contract: dict[str, Any], model_objects: target.ModelObjects
) -> None:
    checks = target.model_contract_checks(contract, model_objects)
    assert checks["status"] == "PASS"
    assert checks["four_equal_plies"] is True
    assert checks["stacking_sequence_deg"] == [0.0, 90.0, 90.0, 0.0]
    assert checks["maximum_relative_residual"] <= target.CONTRACT_RELATIVE_TOLERANCE
    assert checks["scaled_B_residual"] <= 1.0e-12
    assert checks["scaled_I1_residual"] <= 1.0e-12
    expected = checks["expected"]
    assert_allclose(
        [checks["RLB"][name] for name in ("A", "D", "S", "m", "J")],
        [expected[name] for name in ("A", "D", "S", "m", "J")],
        rtol=target.CONTRACT_RELATIVE_TOLERANCE,
        atol=0.0,
    )
    assert_allclose(
        [checks["old_Timoshenko"][name] for name in ("A", "D", "S", "m", "J")],
        [checks["RLB"][name] for name in ("A", "D", "S", "m", "J")],
        rtol=target.CONTRACT_RELATIVE_TOLERANCE,
        atol=0.0,
    )


def test_default_beta_grid_is_exactly_zero_through_ninety_by_one_degree() -> None:
    values = target.beta_grid()
    assert values.shape == (91,)
    assert_allclose(values, np.arange(91, dtype=float), rtol=0.0, atol=0.0)


def test_analysis_entrypoint_has_only_allowed_model_imports_and_no_root_finder() -> None:
    tree = ast.parse(SCRIPT_PATH.read_text(encoding="utf-8"))
    project_imports: set[str] = set()
    locally_defined_functions: set[str] = set()
    dynamic_imports: list[ast.Call] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            project_imports.update(
                alias.name for alias in node.names if alias.name.startswith(("scripts.", "my_project."))
            )
        elif isinstance(node, ast.ImportFrom) and node.module:
            if node.module.startswith(("scripts.", "my_project.")):
                project_imports.add(node.module)
                project_imports.update(
                    f"{node.module}.{alias.name}" for alias in node.names
                )
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            locally_defined_functions.add(node.name)
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
        "layer_order",
        "branch_tracking",
    )
    assert not any(fragment in imported_text for fragment in forbidden)
    assert project_imports == {
        "my_project.analytic.formulas",
        "my_project.analytic.formulas.assemble_clamped_coupled_matrix",
        "my_project.analytic.solvers",
        "my_project.analytic.solvers.find_first_n_roots",
        "scripts.analysis.laminated_beams",
        "scripts.analysis.laminated_beams.validate_reddy_four_ply_isotropic_limit",
        "scripts.lib",
        "scripts.lib.isotropic_rectangular_timoshenko_coupled_beams",
        "scripts.lib.reddy_symmetric_coupled_beams",
        "scripts.lib.reddy_symmetric_laminated_beam",
    }
    assert dynamic_imports == []
    assert not any(
        "root" in name.lower() and name.lower().startswith(("find", "search", "scan"))
        for name in locally_defined_functions
    )
    assert "seed_free_root_inventory" not in locally_defined_functions


@pytest.mark.parametrize("beta_deg", [0.0, 30.0])
def test_frozen_eb_first_eight_plus_guard_and_quality(
    beta_deg: float, contract: dict[str, Any], model_objects: target.ModelObjects
) -> None:
    rows = target._eb_root_rows(
        beta_deg, contract, model_objects, target.frozen_root_policy()
    )
    assert len(rows) == target.OUTPUT_GUARD_POSITION == 9
    assert [int(row["sorted_position"]) for row in rows] == list(range(1, 10))
    assert_allclose(
        [float(row["Lambda"]) for row in rows],
        FROZEN_EB_COMMON_LAMBDA[beta_deg],
        rtol=2.0e-11,
        atol=2.0e-11,
    )
    assert_allclose(
        [float(row["historical_EB_wavenumber"]) ** 2 for row in rows],
        [float(row["Lambda"]) for row in rows],
        rtol=2.0e-16,
        atol=2.0e-14,
    )
    assert all(row["root_status"] == "PASS" for row in rows)
    assert all(
        row["inventory_status"] == "PASS_SIGN_SCAN_CROSSCHECK" for row in rows
    )
    assert max(
        float(row["primary_verification_max_relative"]) for row in rows
    ) <= target.EB_PRIMARY_VERIFICATION_RELATIVE_TOLERANCE
    assert all(
        row["unresolved_candidates_below_internal_guard"]
        == "NOT_ASSESSED_BY_EB_SIGN_SCAN"
        for row in rows
    )
    assert all(
        row["cluster_semantics"] == "NO_CLUSTER_CLAIM_SIGN_SCAN_ONLY"
        for row in rows
    )
    assert max(float(row["scaled_sigma_ratio"]) for row in rows) <= target.ROOT_SINGULAR_RATIO_TOLERANCE
    assert max(float(row["boundary_null_residual"]) for row in rows) <= target.BOUNDARY_RESIDUAL_TOLERANCE
    assert [bool(row["guard_flag"]) for row in rows] == [False] * 8 + [True]
    assert len({str(row["inventory_sha256"]) for row in rows}) == 1


@pytest.mark.parametrize("beta_deg", [0.0, 30.0])
def test_frozen_old_vs_rlb_roots_match_and_pass_direct_boundary_residuals(
    beta_deg: float, contract: dict[str, Any], model_objects: target.ModelObjects
) -> None:
    old_values = FROZEN_ISO_COMMON_LAMBDA[beta_deg][target.MODEL_OLD]
    rlb_values = FROZEN_ISO_COMMON_LAMBDA[beta_deg][target.MODEL_RLB]
    relative = np.abs(old_values - rlb_values) / np.maximum(
        np.maximum(np.abs(old_values), np.abs(rlb_values)), np.finfo(float).tiny
    )
    assert float(np.max(relative)) <= target.OLD_VS_NEW_RELATIVE_TOLERANCE

    policy = target.frozen_root_policy()
    scale = float(contract["frequency"]["Lambda_per_omega"])
    for model, values in (
        (target.MODEL_OLD, old_values),
        (target.MODEL_RLB, rlb_values),
    ):
        provider, _metadata = target.make_matrix_provider(
            model, beta_deg, contract, model_objects
        )
        diagnostics = [
            target.iso_inventory._boundary_matrix_diagnostics(
                float(value), provider, scale, policy
            )
            for value in values
        ]
        assert max(item.scaled_sigma_ratio for item in diagnostics) <= target.ROOT_SINGULAR_RATIO_TOLERANCE
        assert max(item.raw_boundary_null_residual for item in diagnostics) <= target.BOUNDARY_RESIDUAL_TOLERANCE
        assert all(item.root_gate_nullity >= 1 for item in diagnostics)


def _synthetic_plot_rows() -> dict[str, list[dict[str, float | int]]]:
    return {
        model: [
            {
                "beta_deg": beta_deg,
                "sorted_position": position,
                "Lambda": 10.0 * position + beta_deg / 90.0 + 0.01 * model_index,
            }
            for position in range(1, target.PLOTTED_POSITIONS + 1)
            for beta_deg in (0.0, 30.0, 90.0)
        ]
        for model_index, model in enumerate(target.MODELS)
    }


def test_create_plot_has_24_curves_requested_styles_two_legends_and_one_png(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    original_close = target.plt.close
    monkeypatch.setattr(target.plt, "close", lambda *_args, **_kwargs: None)
    output = target.create_plot(_synthetic_plot_rows(), tmp_path / "plot.png")
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
        assert len([child for child in ax.get_children() if isinstance(child, Legend)]) == 2
        assert ax.get_xlabel() == r"$\beta$, degrees"
        assert ax.get_ylabel() == r"$\Lambda$"
        graphics = [
            path
            for path in tmp_path.iterdir()
            if path.suffix.lower() in {".png", ".pdf", ".svg"}
        ]
        assert graphics == [output]
    finally:
        original_close(fig)


def test_existing_generated_output_has_exact_row_counts_and_manifest_hashes() -> None:
    run_manifest_path = RESULT_DIR / "run_manifest.json"
    if not run_manifest_path.is_file():
        pytest.skip("ignored full RLB-2B output has not been generated")

    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    contract = json.loads((RESULT_DIR / "case_contract.json").read_text(encoding="utf-8"))
    assert contract["beta_grid_deg"] == [float(value) for value in range(91)]
    for model in target.MODELS:
        rows = _read_csv(RESULT_DIR / target.ROOT_FILENAMES[model])
        assert len(rows) == 91 * target.OUTPUT_GUARD_POSITION
        assert sum(row["guard_flag"].lower() == "true" for row in rows) == 91
        assert {(float(row["beta_deg"]), int(row["sorted_position"])) for row in rows} == {
            (float(beta), position)
            for beta in range(91)
            for position in range(1, 10)
        }
    assert len(_read_csv(RESULT_DIR / "old_vs_new_rlb_comparison.csv")) == 91 * 8
    assert len(_read_csv(RESULT_DIR / "spectrum_summary.csv")) == 4
    graphics = [
        path
        for path in RESULT_DIR.iterdir()
        if path.suffix.lower() in {".png", ".pdf", ".svg"}
    ]
    assert [path.name for path in graphics] == ["lambda_vs_beta_plot.png"]
    for name, expected_hash in run_manifest["generated_file_hashes"].items():
        assert _sha256(RESULT_DIR / name) == expected_hash
    assert run_manifest["frozen_models_preserved"] is True
    assert run_manifest["statuses"]["OVERALL"] == "PASS"
