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
            for position in range(1, target.PLOTTED_POSITIONS + 1)
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

    for sweep, expected_points, axis_name, representatives in (
        (target.SWEEP_MU, 41, "mu", (0.0, 0.4, 0.8)),
        (target.SWEEP_BETA, 91, "beta_deg", (0.0, 15.0, 90.0)),
    ):
        for model in target.MODELS:
            rows = _read_csv(RESULT_DIR / target.ROOT_FILENAMES[(sweep, model)])
            assert len(rows) == expected_points * target.OUTPUT_GUARD_POSITION
            assert sum(row["guard_flag"].lower() == "true" for row in rows) == expected_points
            _assert_representative_inventory(
                rows,
                axis_name=axis_name,
                representative_values=representatives,
            )

    assert len(_read_csv(RESULT_DIR / "geometry_sanity_checks.csv")) == 3
    assert len(_read_csv(RESULT_DIR / "laminate_properties_summary.csv")) == 6
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
    assert closing["new_points_executed"] == 11
    assert closing["reused_points"] == 385
    assert closing["ready_points_recalculated"] == 0
    assert closing["parallel_workers_used_in_closing_stage"] == 0
    assert closing["inherited_workers_drained"] == 2
    assert closing["global_restarts"] == 0
    assert closing["thread_limits"] == target.CLOSING_THREAD_LIMITS
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
