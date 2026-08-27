"""Fast contract tests for the finite RLB-2C beta workflow.

The tests do not run the 91-angle sweep.  They check the historical frequency
normalization, the exact case-A laminate construction, frozen EB/Timoshenko
and weak-RLB root residuals at two declared angles, and the pure plotting path.
The ignored full output is audited only when its final run manifest exists.
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
import numpy as np
from numpy.testing import assert_allclose
import pytest

from scripts.analysis.laminated_beams import (
    compare_rectangular_weakly_orthotropic_models_vs_beta as target,
)


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = (
    ROOT
    / "scripts"
    / "analysis"
    / "laminated_beams"
    / "compare_rectangular_weakly_orthotropic_models_vs_beta.py"
)
RESULT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "rectangular_weakly_orthotropic_beta_sweep_comparison"
)

# These values are the independently frozen RLB-2B isotropic-baseline
# inventories.  RLB-2C intentionally leaves the EB and old-Timoshenko models
# unchanged.  The arrays contain the first eight Omega values and root 9 guard;
# they are residual probes, not seeds for any new search.
FROZEN_BASELINE_OMEGA = {
    0.0: {
        target.MODEL_EB: np.asarray(
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
    },
    30.0: {
        target.MODEL_EB: np.asarray(
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
    },
}

# These RLB-2C case-A values were frozen by the independent full-grid worker
# before adding this regression probe.  As above, they are used only to
# evaluate the weak-RLB boundary matrix at accepted roots; the root search does
# not read them.
FROZEN_WEAK_RLB_OMEGA = {
    0.0: np.asarray(
        [
            5.773555147079765,
            15.820298594035771,
            30.769807621425430,
            50.364903425318470,
            74.364877899141450,
            102.49466408644187,
            108.87120676984792,
            134.46096399155357,
            169.96245721583230,
        ],
        dtype=float,
    ),
    30.0: np.asarray(
        [
            15.816622538984234,
            19.713491482392190,
            41.386749111423940,
            50.318401601844380,
            77.090077330952740,
            101.49339584150292,
            109.50558443453953,
            132.84161514744900,
            170.09579250721820,
        ],
        dtype=float,
    ),
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


def test_old_lambda_formula_and_historical_eb_mapping(
    contract: dict[str, Any],
    model_objects: target.ModelObjects,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    geometry = contract["geometry"]
    width = float(geometry["width_b"])
    thickness = float(geometry["thickness_h"])
    arm_length = float(geometry["l"])
    area = width * thickness
    inertia = width * thickness**3 / 12.0
    expected_scale = arm_length**2 * math.sqrt(
        target.RHO0 * area / (target.E0 * inertia)
    )
    scale = target.omega_scale(
        E_reference=target.E0,
        rho_reference=target.RHO0,
        width=width,
        thickness=thickness,
        arm_length=arm_length,
    )
    assert_allclose(scale, expected_scale, rtol=2.0e-16, atol=0.0)
    assert_allclose(
        scale,
        arm_length**2
        * math.sqrt(12.0 * target.RHO0 / (target.E0 * thickness**2)),
        rtol=4.0e-16,
        atol=0.0,
    )

    omega = 0.731
    Omega = target.omega_to_Omega(omega, scale)
    Lambda = target.Omega_to_Lambda(Omega)
    assert_allclose(Omega, omega * expected_scale, rtol=2.0e-16, atol=0.0)
    assert_allclose(Lambda, math.sqrt(Omega), rtol=0.0, atol=0.0)
    assert_allclose(target.old_Lambda(omega, scale), Lambda, rtol=0.0, atol=0.0)
    assert_allclose(Lambda**2, Omega, rtol=2.0e-16, atol=0.0)

    captured: dict[str, float] = {}

    def fake_eb_matrix(
        historical_wavenumber: float,
        beta_rad: float,
        mu: float,
        epsilon: float,
    ) -> np.ndarray:
        captured.update(
            wavenumber=float(historical_wavenumber),
            beta_rad=float(beta_rad),
            mu=float(mu),
            epsilon=float(epsilon),
        )
        return np.eye(6, dtype=float)

    monkeypatch.setattr(
        target.rlb2b, "assemble_eb_coupled_matrix", fake_eb_matrix
    )
    provider, metadata = target.make_matrix_provider(
        target.MODEL_EB, 30.0, contract, model_objects
    )
    assert_allclose(provider(omega), np.eye(6), rtol=0.0, atol=0.0)
    assert_allclose(captured["wavenumber"], Lambda, rtol=0.0, atol=2.0e-15)
    assert_allclose(captured["beta_rad"], math.radians(30.0), rtol=0.0, atol=0.0)
    assert captured["mu"] == 0.0
    assert_allclose(
        captured["epsilon"], thickness / (math.sqrt(12.0) * arm_length)
    )
    assert metadata["matrix_argument"] == "Lambda=sqrt(Omega)"


def test_root_row_transform_preserves_omega_Omega_Lambda_identity() -> None:
    scale = 69.2820323027551
    source = {
        "model": target.MODEL_EB,
        "beta_deg": 0.0,
        "sorted_position": 1,
        "role": "FIRST_8",
        "omega": 9.0 / scale,
        "Lambda": 9.0,
        "historical_EB_wavenumber": 3.0,
        "bracket_left_Lambda": 4.0,
        "bracket_right_Lambda": 16.0,
        "internal_root13_Lambda": 25.0,
    }
    row = target.transform_root_row(source, target.MODEL_EB, scale)
    assert row["omega"] == source["omega"]
    assert row["Omega"] == 9.0
    assert row["Lambda"] == 3.0
    assert row["bracket_left_Omega"] == 4.0
    assert row["bracket_right_Omega"] == 16.0
    assert row["bracket_left_Lambda"] == 2.0
    assert row["bracket_right_Lambda"] == 4.0
    assert row["internal_root13_Omega"] == 25.0
    assert row["internal_root13_Lambda"] == 5.0
    assert row["normalization_identity_relative_residual"] <= 2.0e-16
    assert row["historical_EB_mapping_relative_residual"] == 0.0


def test_shared_geometry_and_reference_material_contract(
    contract: dict[str, Any], model_objects: target.ModelObjects
) -> None:
    assert contract["mu"] == 0.0
    assert contract["tau"] == 0.0
    assert contract["reference_isotropic_material"] == {
        "E0": 1.0,
        "nu0": 0.3,
        "rho0": 1.0,
        "G0": 1.0 / 2.6,
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
    assert contract["old_models"] == {
        "EB": "isotropic rectangular baseline",
        "old_Timoshenko": "isotropic rectangular baseline",
        "K": 5.0 / 6.0,
    }
    baseline = model_objects.baseline_contract
    assert baseline["material"] == {
        "E": 1.0,
        "nu": 0.3,
        "rho": 1.0,
        "G": 1.0 / 2.6,
    }
    assert model_objects.baseline_objects.legacy_section.K == 5.0 / 6.0


def test_case_A_exact_four_ply_laminate_and_constitutive_sanity(
    contract: dict[str, Any], model_objects: target.ModelObjects
) -> None:
    material_contract = contract["new_RLB_lamina"]
    assert material_contract == {
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
    laminate = model_objects.weak_laminate
    material = laminate.plies[0].material
    assert target.DELTA == 0.1
    assert material.E1 == 1.1
    assert material.E2 == 0.9
    assert material.nu12 == 0.3
    assert material.G12 == material.G13 == material.G23 == 1.0 / 2.6
    assert material.rho == 1.0
    assert len(laminate.plies) == 4
    assert [ply.angle_deg for ply in laminate.plies] == [0.0, 90.0, 90.0, 0.0]
    assert_allclose(
        [ply.thickness for ply in laminate.plies],
        np.full(4, 0.05 / 4.0),
        rtol=0.0,
        atol=0.0,
    )
    assert_allclose(
        laminate.z_interfaces,
        [-0.025, -0.0125, 0.0, 0.0125, 0.025],
        rtol=0.0,
        atol=0.0,
    )

    checks = target.constitutive_checks(contract, model_objects)
    assert checks["status"] == "PASS"
    assert checks["exact_case_A_material"] is True
    assert checks["exact_four_equal_plies"] is True
    assert checks["B_relative"] <= target.SYMMETRY_RELATIVE_TOLERANCE
    assert checks["I1_relative"] <= target.SYMMETRY_RELATIVE_TOLERANCE
    assert checks["mass_max_relative_residual"] <= target.CONTRACT_RELATIVE_TOLERANCE
    assert checks["shear_relative_residual"] <= target.CONTRACT_RELATIVE_TOLERANCE
    assert checks["analytic_A_relative_residual"] <= 1.0e-12
    assert checks["analytic_D_relative_residual"] <= 1.0e-12
    assert checks["stiffnesses_differ_from_isotropic_reference"] is True
    assert_allclose(checks["mass_ratio_to_isotropic"], 1.0, rtol=4.0e-16)
    assert_allclose(
        checks["rotary_inertia_ratio_to_isotropic"], 1.0, rtol=4.0e-16
    )
    assert_allclose(checks["Sbeam_ratio_to_isotropic"], 1.0, rtol=4.0e-16)
    assert_allclose(
        model_objects.weak_properties.A,
        0.01000794896957802,
        rtol=2.0e-15,
        atol=0.0,
    )
    assert_allclose(
        model_objects.weak_properties.D,
        2.240366593286123e-6,
        rtol=2.0e-15,
        atol=0.0,
    )
    assert_allclose(
        model_objects.weak_properties.S,
        0.0032051282051282055,
        rtol=2.0e-15,
        atol=0.0,
    )
    assert_allclose(model_objects.weak_properties.m, 0.01, rtol=4.0e-16)
    assert_allclose(
        model_objects.weak_properties.J,
        2.083333333333334e-6,
        rtol=4.0e-16,
    )


def test_default_beta_grid_is_exactly_zero_through_ninety_by_one_degree() -> None:
    values = target.beta_grid()
    assert values.shape == (91,)
    assert_allclose(values, np.arange(91, dtype=float), rtol=0.0, atol=0.0)


def test_analysis_entrypoint_has_only_allowed_imports_and_no_local_root_finder() -> None:
    tree = ast.parse(SCRIPT_PATH.read_text(encoding="utf-8"))
    project_imports: set[str] = set()
    locally_defined_functions: set[str] = set()
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
        "shellbuckling",
        "shell_buckling",
        "layer_order",
        "branch_tracking",
    )
    assert not any(fragment in imported_text for fragment in forbidden)
    assert project_imports == {
        "scripts.analysis.laminated_beams",
        "scripts.analysis.laminated_beams.compare_rectangular_isotropic_models_vs_beta",
        "scripts.lib",
        "scripts.lib.reddy_symmetric_coupled_beams",
        "scripts.lib.reddy_symmetric_laminated_beam",
    }
    assert dynamic_imports == []
    assert not any(
        "root" in name.lower()
        and name.lower().startswith(("find", "search", "scan"))
        for name in locally_defined_functions
    )
    assert "seed_free_root_inventory" not in locally_defined_functions


@pytest.mark.parametrize("beta_deg", [0.0, 30.0])
@pytest.mark.parametrize("model", target.MODELS)
def test_frozen_first_eight_plus_guard_have_direct_matrix_residuals(
    beta_deg: float,
    model: str,
    contract: dict[str, Any],
    model_objects: target.ModelObjects,
) -> None:
    values = (
        FROZEN_WEAK_RLB_OMEGA[beta_deg]
        if model == target.MODEL_RLB
        else FROZEN_BASELINE_OMEGA[beta_deg][model]
    )
    assert values.shape == (target.OUTPUT_GUARD_POSITION,) == (9,)
    assert np.all(values > 0.0)
    assert np.all(np.diff(values) > 0.0)
    Lambdas = np.sqrt(values)
    assert_allclose(Lambdas**2, values, rtol=2.0e-16, atol=2.0e-14)

    provider, _metadata = target.make_matrix_provider(
        model, beta_deg, contract, model_objects
    )
    policy = target.rlb2b.frozen_root_policy()
    scale = float(contract["frequency"]["Omega_per_omega"])
    diagnostics = [
        target.rlb2b.iso_inventory._boundary_matrix_diagnostics(
            float(Omega), provider, scale, policy
        )
        for Omega in values
    ]
    assert max(item.scaled_sigma_ratio for item in diagnostics) <= target.rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE
    assert max(item.raw_boundary_null_residual for item in diagnostics) <= target.rlb2b.BOUNDARY_RESIDUAL_TOLERANCE
    assert all(item.root_gate_nullity >= 1 for item in diagnostics)


def _synthetic_plot_rows() -> dict[str, list[dict[str, float | int]]]:
    return {
        model: [
            {
                "beta_deg": beta_deg,
                "sorted_position": position,
                "Lambda": 3.0 * position + beta_deg / 90.0 + 0.01 * model_index,
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
        assert len(
            [child for child in ax.get_children() if isinstance(child, Legend)]
        ) == 2
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


def test_existing_generated_output_has_exact_rows_hashes_and_statuses() -> None:
    run_manifest_path = RESULT_DIR / "run_manifest.json"
    if not run_manifest_path.is_file():
        pytest.skip("ignored full RLB-2C output has not been finalized")

    run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
    contract = json.loads(
        (RESULT_DIR / "case_contract.json").read_text(encoding="utf-8")
    )
    assert contract["beta_grid_deg"] == [float(value) for value in range(91)]
    assert contract["frequency"]["Lambda_definition"] == "sqrt(Omega)"
    assert contract["new_RLB_lamina"]["case_id"] == "A"
    assert contract["new_RLB_lamina"]["delta"] == 0.1
    for model in target.MODELS:
        rows = _read_csv(RESULT_DIR / target.ROOT_FILENAMES[model])
        assert len(rows) == 91 * target.OUTPUT_GUARD_POSITION
        assert sum(row["guard_flag"].lower() == "true" for row in rows) == 91
        assert {
            (float(row["beta_deg"]), int(row["sorted_position"]))
            for row in rows
        } == {
            (float(beta), position)
            for beta in range(91)
            for position in range(1, 10)
        }
        assert max(
            float(row["normalization_identity_relative_residual"])
            for row in rows
        ) <= target.NORMALIZATION_IDENTITY_TOLERANCE
        assert_allclose(
            [float(row["Lambda"]) ** 2 for row in rows],
            [float(row["Omega"]) for row in rows],
            rtol=target.NORMALIZATION_IDENTITY_TOLERANCE,
            atol=0.0,
        )
    assert len(_read_csv(RESULT_DIR / "laminate_properties.csv")) == 2
    assert len(_read_csv(RESULT_DIR / "old_vs_new_rlb_comparison.csv")) == 91 * 8
    assert len(_read_csv(RESULT_DIR / "eb_vs_new_rlb_comparison.csv")) == 91 * 8
    assert len(_read_csv(RESULT_DIR / "spectrum_summary.csv")) == 5
    graphics = [
        path
        for path in RESULT_DIR.iterdir()
        if path.suffix.lower() in {".png", ".pdf", ".svg"}
    ]
    assert [path.name for path in graphics] == ["lambda_vs_beta_plot.png"]
    for name, expected_hash in run_manifest["generated_file_hashes"].items():
        assert _sha256(RESULT_DIR / name) == expected_hash
    assert run_manifest["frozen_models_preserved"] is True
    assert run_manifest["delta_values_run"] == [0.1]
    assert run_manifest["case_ids_run"] == ["A"]
    assert run_manifest["mu"] == 0.0
    assert run_manifest["tau"] == 0.0
    assert run_manifest["Ritz_run"] is False
    assert run_manifest["FEM_run"] is False
    assert run_manifest["branch_tracking_run"] is False
    statuses = run_manifest["statuses"]
    assert set(statuses) == {
        "OVERALL",
        "RLB-2C-LAMBDA-DEFINITION",
        "RLB-2C-PLOT-GENERATION",
        "RLB-2C-ROOT-DATA",
        "RLB-2C-WEAK-ORTHOTROPIC-RLB-RUN",
    }
    assert statuses["RLB-2C-LAMBDA-DEFINITION"] == "PASS"
    assert statuses["RLB-2C-PLOT-GENERATION"] == "PASS"
    assert statuses["RLB-2C-ROOT-DATA"] in {"PASS", "PARTIAL_PASS"}
    assert statuses["RLB-2C-WEAK-ORTHOTROPIC-RLB-RUN"] in {
        "PASS",
        "PARTIAL_PASS",
    }
    expected_overall = (
        "PASS"
        if all(
            value == "PASS"
            for name, value in statuses.items()
            if name != "OVERALL"
        )
        else "FAIL"
    )
    assert statuses["OVERALL"] == expected_overall
    assert run_manifest["summary"]["exact_full_grid_passed"] is True
    assert run_manifest["summary"]["comparison_counts_passed"] is True
    assert run_manifest["frequency_coordinates"]["root_search_by_model"] == {
        target.MODEL_EB: "historical Lambda=sqrt(Omega)",
        target.MODEL_OLD: "Omega",
        target.MODEL_RLB: "Omega",
    }
