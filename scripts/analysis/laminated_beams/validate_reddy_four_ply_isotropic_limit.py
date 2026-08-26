"""RLB-1C-ISO four-ply isotropic-limit validation driver.

This script deliberately separates the two determinant inventories into fresh
subprocesses.  The RLB worker imports only the frozen RLB model; the legacy
worker imports only the independent closed-form rectangular comparator.  Both
workers read the same frozen JSON case contract, perform seed-free primary and
verification scans, and freeze their inventory digest before any future
cross-model comparison.

The current CLI implements the manifest, preliminary algebraic/local gates,
and independent root inventories.  Spectrum/mode comparison and report prose
are intentionally left for the next implementation increment.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
from datetime import datetime, timezone
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq, minimize_scalar


FloatArray = NDArray[np.float64]
MatrixProvider = Callable[[float], FloatArray]

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))
DEFAULT_CONTRACT_PATH = (
    REPOSITORY_ROOT / "tests" / "data" / "reddy_four_ply_isotropic_limit_cases.json"
)
DEFAULT_OUTPUT_DIRECTORY = (
    REPOSITORY_ROOT
    / "results"
    / "laminated_beams"
    / "reddy_four_ply_isotropic_limit_validation"
)
STAGE_ID = "RLB-1C-ISO"
IMPLEMENTATION_SCOPE = "manifest_preliminary_gates_and_independent_inventories_v1"

FROZEN_REFERENCE_PATHS = (
    "scripts/lib/reddy_symmetric_laminated_beam.py",
    "scripts/lib/reddy_inplane_geometry.py",
    "scripts/lib/reddy_symmetric_coupled_beams.py",
    "scripts/lib/reddy_symmetric_coupled_beams_ritz.py",
    "scripts/analysis/laminated_beams/pilot_reddy_symmetric_coupled_beams_beta0.py",
    "scripts/analysis/laminated_beams/validate_reddy_symmetric_coupled_beams_nonzero_beta.py",
    "scripts/lib/variable_length_timoshenko.py",
    "src/my_project/analytic/formulas.py",
    "src/my_project/analytic/solvers.py",
    "scripts/lib/analytic_coupled_rods_shapes.py",
    "scripts/lib/in_plane_shape_geometry.py",
    "tests/test_reddy_symmetric_laminated_beam.py",
    "tests/test_reddy_inplane_geometry.py",
    "tests/test_reddy_symmetric_coupled_beams_beta0.py",
    "tests/test_reddy_symmetric_coupled_beams_ritz.py",
    "tests/test_reddy_table_4_3_3_discrepancy_audit.py",
    "tests/test_timoshenko_basis_regime_extension.py",
    "tests/test_common_thickness_scaling.py",
    "docs/laminated_beams/reddy_ch4_source_contract.md",
    "docs/laminated_beams/reddy_symmetric_single_beam_validation.md",
    "docs/laminated_beams/reddy_table_4_3_3_discrepancy_audit.md",
    "docs/laminated_beams/reddy_inplane_coordinate_contract.md",
    "docs/laminated_beams/reddy_symmetric_rigid_joint.md",
    "docs/laminated_beams/reddy_symmetric_coupled_beta0_validation.md",
    "docs/laminated_beams/reddy_symmetric_coupled_nonzero_beta_validation.md",
    "results/laminated_beams/reddy_symmetric_coupled_beta0_pilot/model_manifest.json",
    "results/laminated_beams/reddy_symmetric_coupled_beta0_pilot/run_manifest.json",
    "results/laminated_beams/reddy_symmetric_coupled_nonzero_beta_validation/model_manifest.json",
    "results/laminated_beams/reddy_symmetric_coupled_nonzero_beta_validation/run_manifest.json",
)

PRELIMINARY_MODEL_PATHS = (
    "scripts/lib/reddy_symmetric_laminated_beam.py",
    "scripts/lib/reddy_inplane_geometry.py",
    "scripts/lib/reddy_symmetric_coupled_beams.py",
    "scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py",
    "scripts/lib/variable_length_timoshenko.py",
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest().upper()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def _json_value(value: Any) -> Any:
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float):
        if math.isfinite(value):
            return value
        if math.isnan(value):
            return "NaN"
        return "Infinity" if value > 0.0 else "-Infinity"
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_json_value(payload), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _csv_value(value: Any) -> Any:
    if isinstance(value, (tuple, list, dict, np.ndarray)):
        return json.dumps(_json_value(value), ensure_ascii=False, separators=(",", ":"))
    if isinstance(value, np.generic):
        return value.item()
    return value


def _write_csv(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    data = [dict(row) for row in rows]
    fields: list[str] = []
    for row in data:
        for key in row:
            if key not in fields:
                fields.append(key)
    if not fields:
        fields = ["status"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in data:
            writer.writerow({key: _csv_value(row.get(key, "")) for key in fields})


def _load_contract(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    expected: dict[str, Any] = {
        "contract_version": "rlb_1c_iso_v1",
        "material": {"E": 1.0, "rho": 1.0, "nu": 0.3, "K": 5.0 / 6.0},
        "lengths": {"L_total": 2.0, "L_ref": 1.0},
        "geometries": {
            "G20": {
                "width": 0.2,
                "thickness": 0.05,
                "width_to_thickness": 4.0,
                "L_ref_to_thickness": 20.0,
            },
            "G10": {
                "width": 0.4,
                "thickness": 0.1,
                "width_to_thickness": 4.0,
                "L_ref_to_thickness": 10.0,
            },
        },
        "spectral_stack_deg": [0.0, 90.0, 90.0, 0.0],
        "algebraic_control_stacks_deg": [
            [0.0, 0.0, 0.0, 0.0],
            [17.0, -38.0, -38.0, 17.0],
        ],
        "cases": [
            {"case_id": "ISO-01", "geometry": "G20", "L1": 1.0, "L2": 1.0, "beta_deg": 0.0},
            {"case_id": "ISO-02", "geometry": "G20", "L1": 0.7, "L2": 1.3, "beta_deg": 0.0},
            {"case_id": "ISO-03", "geometry": "G20", "L1": 1.0, "L2": 1.0, "beta_deg": 30.0},
            {"case_id": "ISO-04", "geometry": "G20", "L1": 1.0, "L2": 1.0, "beta_deg": 90.0},
            {"case_id": "ISO-05", "geometry": "G20", "L1": 1.0, "L2": 1.0, "beta_deg": -30.0},
            {"case_id": "ISO-06", "geometry": "G20", "L1": 0.7, "L2": 1.3, "beta_deg": 30.0},
            {"case_id": "ISO-07", "geometry": "G20", "L1": 1.3, "L2": 0.7, "beta_deg": 30.0},
            {"case_id": "ISO-08", "geometry": "G10", "L1": 1.0, "L2": 1.0, "beta_deg": 30.0},
        ],
        "mode_positions": {"ISO-03": [1, 2, 3, 6], "ISO-04": [1, 2, 3]},
        "local_solution_space_policy": {
            "L_local": 0.25,
            "cutoff_factors": [0.23, 0.5, 0.73, 0.99, 1.0, 1.01, 1.5, 2.2],
            "geometries": ["G20", "G10"],
            "rationale": "conditioning-controlled local gate",
        },
        "thresholds": {
            "qbar_matrix_analytic": 1.0e-12,
            "section_reduction": 1.0e-11,
            "legacy_circular_regular": 1.0e-11,
            "legacy_circular_cutoff": 1.0e-8,
            "local_solution_space_regular": 1.0e-10,
            "local_solution_space_cutoff": 1.0e-8,
            "isolated_frequency_comparison": 1.0e-8,
            "cluster_center": 1.0e-8,
            "symmetry_spectrum": 1.0e-9,
            "root_singular_ratio": 1.0e-9,
            "boundary_null_residual": 1.0e-9,
            "outer_clamp_residual": 1.0e-10,
            "joint_compatibility": 1.0e-10,
            "joint_equilibrium": 1.0e-9,
            "energy_identity": 1.0e-8,
            "isolated_one_minus_MAC": 1.0e-6,
            "grid_convergence": 1.0e-8,
        },
        "root_policy": {
            "requested_roots": 12,
            "guard_roots": 1,
            "Omega_min": 1.0e-8,
            "Omega_max": 1600.0,
            "primary_scan_points": 8001,
            "verification_scan_points": 16001,
            "primary_phase": 0.0,
            "verification_phase": 0.5,
            "sigma_prefilter": 1.0e-5,
            "nullity_relative_threshold": 1.0e-12,
            "root_xtol_Omega": 1.0e-11,
            "dedup_atol_Omega": 5.0e-10,
            "dedup_rtol": 5.0e-12,
            "cluster_atol_Omega": 1.0e-10,
            "cluster_rtol": 1.0e-10,
            "post_guard_tail_Omega": 2.0,
            "local_close_pair_guard_subintervals": 32,
        },
    }
    expected_keys = set(expected) | {"explicit_exclusions"}
    if set(payload) != expected_keys:
        raise ValueError(
            "The frozen contract top-level schema changed: "
            f"expected {sorted(expected_keys)}, got {sorted(payload)}."
        )
    for field, expected_value in expected.items():
        if payload.get(field) != expected_value:
            raise ValueError(
                f"Frozen RLB-1C-ISO contract field {field!r} changed: "
                f"expected {expected_value!r}, got {payload.get(field)!r}."
            )
    expected_exclusions = {
        "one_ply_numerical_case",
        "Rayleigh_Ritz",
        "Euler_Bernoulli_comparison",
        "FEM",
        "torsion",
        "damping",
        "complex_frequencies",
        "nonsymmetric_laminate",
        "B_nonzero",
        "joint_mass_or_compliance",
        "parameter_sweep",
        "branch_tracking",
        "article_work",
        "commit",
        "push",
    }
    if set(payload["explicit_exclusions"]) != expected_exclusions:
        raise ValueError("The frozen explicit-exclusions set changed.")
    return payload


def _git_output(*arguments: str) -> str:
    result = subprocess.run(
        ["git", *arguments],
        cwd=REPOSITORY_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def _git_snapshot() -> dict[str, Any]:
    return {
        "repository_root": _git_output("rev-parse", "--show-toplevel"),
        "branch": _git_output("branch", "--show-current"),
        "head": _git_output("rev-parse", "HEAD"),
        "latest_commit": _git_output("log", "-1", "--oneline"),
        "status_short": _git_output("status", "--short").splitlines(),
    }


def _frozen_reference_hashes() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for relative in FROZEN_REFERENCE_PATHS:
        path = REPOSITORY_ROOT / relative
        rows.append(
            {
                "path": relative,
                "exists": path.is_file(),
                "sha256": _sha256_file(path) if path.is_file() else None,
            }
        )
    return rows


def _preliminary_model_hashes() -> dict[str, str]:
    return {
        relative: _sha256_file(REPOSITORY_ROOT / relative)
        for relative in PRELIMINARY_MODEL_PATHS
    }


def _base_model_manifest(contract_path: Path, contract: Mapping[str, Any]) -> dict[str, Any]:
    comparator_path = REPOSITORY_ROOT / "scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py"
    return {
        "stage": STAGE_ID,
        "implementation_scope": IMPLEMENTATION_SCOPE,
        "generated_utc": _utc_now(),
        "git_at_generation": _git_snapshot(),
        "contract_source": str(contract_path.relative_to(REPOSITORY_ROOT)).replace("\\", "/"),
        "contract_sha256": _sha256_file(contract_path),
        "thresholds_frozen_before_calculation": contract["thresholds"],
        "root_policy_frozen_before_calculation": contract["root_policy"],
        "material_contract": contract["material"],
        "geometry_contract": contract["geometries"],
        "section_orientation": {
            "width_b": "local Reddy e_y=-k, perpendicular to construction plane",
            "thickness_h": "local Reddy e_z=n_i, in construction plane",
            "bending_axis": "e_y",
            "area": "b*h",
            "I_y": "b*h^3/12",
        },
        "four_ply_numerical_contract": {
            "number_of_plies": 4,
            "equal_ply_thickness": "h/4",
            "spectral_stack_deg": contract["spectral_stack_deg"],
            "algebraic_control_stacks_deg": contract["algebraic_control_stacks_deg"],
            "one_ply_numerical_case": False,
            "isotropic_shortcut_in_RLB_pipeline": False,
        },
        "workflow_separation": {
            "inventory_processes": ["rlb", "legacy_rectangular"],
            "fresh_subprocess_per_workflow": True,
            "cross_seeding": False,
            "comparison_implemented_in_this_increment": False,
            "rlb_worker_imports_legacy_comparator": False,
            "legacy_worker_imports_RLB": False,
        },
        "independent_comparator": {
            "path": "scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py",
            "sha256": _sha256_file(comparator_path) if comparator_path.is_file() else None,
            "uses_matrix_exponential": False,
            "uses_RLB_properties_or_matrices": False,
        },
        "frozen_reference_hashes": _frozen_reference_hashes(),
        "historical_RLB_1C_Ritz_status_preserved": "FAIL_CONVERGENCE_AT_N_LE_16",
        "explicit_exclusions": contract["explicit_exclusions"],
        "preliminary_statuses": {
            "RLB-ISO-4PLY-CONSTITUTIVE": "NOT_RUN",
            "RLB-ISO-SECTION-REDUCTION": "NOT_RUN",
            "RLB-ISO-LEGACY-RECTANGULAR-ADAPTER": "NOT_RUN",
            "RLB-ISO-LOCAL-ARM-EQUIVALENCE": "NOT_RUN",
        },
    }


def write_manifest_only(contract_path: Path, output_directory: Path) -> dict[str, Any]:
    contract = _load_contract(contract_path)
    output_directory.mkdir(parents=True, exist_ok=True)
    manifest = _base_model_manifest(contract_path, contract)
    _write_json(output_directory / "model_manifest.json", manifest)
    _write_json(output_directory / "case_contract.json", contract)
    return manifest


def _relative_residual(actual: Any, expected: Any) -> float:
    left = np.asarray(actual, dtype=float)
    right = np.asarray(expected, dtype=float)
    numerator = float(np.linalg.norm(left - right))
    denominator = max(
        float(np.linalg.norm(left)),
        float(np.linalg.norm(right)),
        np.finfo(float).tiny,
    )
    return numerator / denominator


def _scaled_zero_residual(actual: Any, scale: float) -> float:
    return float(np.linalg.norm(np.asarray(actual, dtype=float))) / max(
        abs(float(scale)), np.finfo(float).tiny
    )


def _stack_label(angles: Sequence[float]) -> str:
    return "[" + "/".join(f"{float(angle):g}" for angle in angles) + "]"


def _isotropic_material(single_beam: Any, material: Mapping[str, Any]) -> Any:
    E = float(material["E"])
    nu = float(material["nu"])
    shear = E / (2.0 * (1.0 + nu))
    return single_beam.LaminaMaterial(
        E1=E,
        E2=E,
        nu12=nu,
        G12=shear,
        G13=shear,
        G23=shear,
        rho=float(material["rho"]),
        name="RLB-1C-ISO isotropic lamina",
    )


def _four_ply_properties(
    single_beam: Any,
    contract: Mapping[str, Any],
    geometry_id: str,
    angles: Sequence[float],
) -> tuple[Any, Any, Any]:
    geometry = contract["geometries"][geometry_id]
    h = float(geometry["thickness"])
    material = _isotropic_material(single_beam, contract["material"])
    plies = tuple(
        single_beam.Ply(material, float(angle), h / 4.0, label=f"ply-{index}")
        for index, angle in enumerate(angles, start=1)
    )
    laminate = single_beam.integrate_laminate(plies)
    properties = single_beam.reduce_to_beam_properties(
        laminate,
        width=float(geometry["width"]),
        K=float(contract["material"]["K"]),
        symmetry_tolerance=1.0e-12,
        reduction_tolerance=1.0e-12,
    )
    return material, laminate, properties


def _orthonormal_projector(matrix: FloatArray) -> tuple[FloatArray, int]:
    u, singular, _vh = np.linalg.svd(np.asarray(matrix, dtype=float), full_matrices=False)
    threshold = max(matrix.shape) * np.finfo(float).eps * float(singular[0])
    rank = int(np.count_nonzero(singular > threshold))
    basis = u[:, :rank]
    return basis @ basis.T, rank


def _legacy_first_order_state_matrix(omega: float, section: Any) -> FloatArray:
    w2 = float(omega) ** 2
    return np.array(
        [
            [0.0, 0.0, 0.0, 1.0 / section.EA, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 1.0 / section.KGA, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / section.EI],
            [-section.rhoA * w2, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, -section.rhoA * w2, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, -section.rhoI * w2, 0.0, -1.0, 0.0],
        ],
        dtype=float,
    )


def _row_normalized_determinant(matrix: FloatArray) -> float:
    array = np.asarray(matrix, dtype=float)
    norms = np.linalg.norm(array, axis=1)
    scaled = np.array(array, copy=True)
    nonzero = norms > 0.0
    scaled[nonzero] /= norms[nonzero, None]
    return float(np.linalg.det(scaled))


def _legacy_circular_first_six_root_rows(
    comparator: Any,
    frozen_circular: Any,
    generic_circle: Any,
    *,
    epsilon: float,
    beta_deg: float,
    tolerance: float,
) -> list[dict[str, Any]]:
    """Independently bracket the generic adapter's first six circular roots."""

    scan_start = 0.2
    scan_stop = 45.0
    scan_step = 0.005
    required = 6

    def generic_determinant(lambda_value: float) -> float:
        omega = frozen_circular.project_omega(float(lambda_value), epsilon)
        return comparator.legacy_coupled_characteristic_determinant(
            omega,
            generic_circle,
            frozen_circular.L_SEGMENT,
            generic_circle,
            frozen_circular.L_SEGMENT,
            beta_deg=beta_deg,
            scaled=True,
        )

    grid = np.arange(scan_start, scan_stop + 0.5 * scan_step, scan_step)
    determinant_values = np.asarray(
        [generic_determinant(float(value)) for value in grid], dtype=float
    )
    generic_roots: list[float] = []
    generic_brackets: list[tuple[float, float]] = []
    for left, right, value_left, value_right in zip(
        grid[:-1], grid[1:], determinant_values[:-1], determinant_values[1:]
    ):
        if not (math.isfinite(float(value_left)) and math.isfinite(float(value_right))):
            continue
        candidate: float | None = None
        if value_left == 0.0:
            candidate = float(left)
        elif np.signbit(value_left) != np.signbit(value_right):
            candidate = float(
                brentq(
                    generic_determinant,
                    float(left),
                    float(right),
                    xtol=1.0e-13,
                    rtol=1.0e-13,
                    maxiter=160,
                )
            )
        if candidate is None:
            continue
        if generic_roots and abs(candidate - generic_roots[-1]) <= 1.0e-7:
            continue
        generic_roots.append(candidate)
        generic_brackets.append((float(left), float(right)))
        if len(generic_roots) == required:
            break
    if len(generic_roots) != required:
        raise RuntimeError(
            f"Generic circular adapter found {len(generic_roots)} roots; expected {required}."
        )

    frozen_roots, warnings = frozen_circular.timo_sorted_roots(
        beta_deg,
        0.0,
        epsilon,
        required,
        eta=0.0,
    )
    frozen_roots = np.asarray(frozen_roots, dtype=float)
    if frozen_roots.shape != (required,) or not np.all(np.isfinite(frozen_roots)):
        raise RuntimeError("Frozen circular reference did not return six finite roots.")

    rows: list[dict[str, Any]] = []
    for root_index, (old_root, new_root, bracket) in enumerate(
        zip(frozen_roots, generic_roots, generic_brackets, strict=True), start=1
    ):
        relative = _relative_difference(float(old_root), float(new_root))
        rows.append(
            {
                "check_kind": "first_six_root_backcompat",
                "epsilon": epsilon,
                "beta_deg": beta_deg,
                "root_index": root_index,
                "old_Lambda": float(old_root),
                "new_Lambda": float(new_root),
                "old_omega": frozen_circular.project_omega(float(old_root), epsilon),
                "new_omega": frozen_circular.project_omega(float(new_root), epsilon),
                "generic_scan_start_Lambda": scan_start,
                "generic_scan_stop_Lambda": scan_stop,
                "generic_scan_step_Lambda": scan_step,
                "generic_bracket_left_Lambda": bracket[0],
                "generic_bracket_right_Lambda": bracket[1],
                "cross_seeded": False,
                "frozen_reference_warnings": warnings,
                "root_relative_residual": relative,
                "threshold": tolerance,
                "maximum_residual": relative,
                "status": "PASS" if relative <= tolerance else "FAIL",
            }
        )
    return rows


def run_preflight_snapshot_worker(
    worker: str, contract_path: Path, output_directory: Path
) -> dict[str, Any]:
    """Freeze one model's preliminary values in a fresh, isolated process."""

    if worker not in {"rlb", "legacy"}:
        raise ValueError("preflight worker must be 'rlb' or 'legacy'.")
    _assert_fresh_worker_import_boundary(worker, after_import=False)
    contract = _load_contract(contract_path)
    material_data = contract["material"]
    local_policy = contract["local_solution_space_policy"]
    local_length = float(local_policy["L_local"])
    factors = tuple(float(value) for value in local_policy["cutoff_factors"])
    snapshot_rows: list[dict[str, Any]] = []
    if worker == "rlb":
        single_beam = importlib.import_module(
            "scripts.lib.reddy_symmetric_laminated_beam"
        )
        coupled = importlib.import_module("scripts.lib.reddy_symmetric_coupled_beams")
        material = _isotropic_material(single_beam, material_data)
        for angle in (0.0, 17.0, 45.0, 90.0, -38.0):
            snapshot_rows.append(
                {
                    "kind": "Qbar",
                    "angle_deg": angle,
                    "Qbar": single_beam.transformed_reduced_stiffness(material, angle),
                    "shear_Qbar": single_beam.transformed_transverse_shear_stiffness(
                        material, angle
                    ),
                }
            )
        stacks = [
            tuple(float(value) for value in contract["spectral_stack_deg"]),
            *(
                tuple(float(value) for value in stack)
                for stack in contract["algebraic_control_stacks_deg"]
            ),
        ]
        for geometry_id in contract["geometries"]:
            for stack in stacks:
                _material, laminate, properties = _four_ply_properties(
                    single_beam, contract, geometry_id, stack
                )
                snapshot_rows.append(
                    {
                        "kind": "four_ply_section",
                        "geometry": geometry_id,
                        "stack_deg": stack,
                        "z_interfaces": laminate.z_interfaces,
                        "A": laminate.A,
                        "B": laminate.B,
                        "D": laminate.D,
                        "shear": laminate.shear,
                        "I0": laminate.I0,
                        "I1": laminate.I1,
                        "I2": laminate.I2,
                        "beam_properties": {
                            "A": properties.A,
                            "D": properties.D,
                            "S": properties.S,
                            "m": properties.m,
                            "J": properties.J,
                        },
                    }
                )
            _material, _laminate, properties = _four_ply_properties(
                single_beam,
                contract,
                geometry_id,
                tuple(float(value) for value in contract["spectral_stack_deg"]),
            )
            cutoff = math.sqrt(properties.S / properties.J)
            for factor in factors:
                omega = factor * cutoff
                snapshot_rows.append(
                    {
                        "kind": "local_arm",
                        "geometry": geometry_id,
                        "cutoff_factor": factor,
                        "omega": omega,
                        "state_matrix": single_beam.combined_state_matrix(
                            omega, properties
                        ),
                        "clamp_to_endpoint": coupled.arm_clamp_map(
                            omega, local_length, properties
                        ),
                    }
                )
    else:
        comparator = importlib.import_module(
            "scripts.lib.isotropic_rectangular_timoshenko_coupled_beams"
        )
        frozen_circular = importlib.import_module(
            "scripts.lib.variable_length_timoshenko"
        )
        epsilon = 0.05
        old_section = frozen_circular.section_from_epsilon(epsilon)
        circle = comparator.circular_section(
            E=frozen_circular.E,
            nu=frozen_circular.NU,
            rho=frozen_circular.RHO,
            radius=old_section.radius,
            K=old_section.kappa,
        )
        cutoff_lambda = frozen_circular.lambda_cutoff(epsilon, old_section)
        for factor in (0.8, 1.0, 1.2):
            Lambda = factor * cutoff_lambda
            omega = frozen_circular.project_omega(Lambda, epsilon)
            basis = comparator.timoshenko_spatial_basis(omega, circle)
            for beta_deg in (0.0, 30.0):
                snapshot_rows.append(
                    {
                        "kind": "circular_backcompat",
                        "cutoff_factor": factor,
                        "Lambda": Lambda,
                        "omega": omega,
                        "beta_deg": beta_deg,
                        "basis": basis.__dict__,
                        "boundary_matrix": comparator.legacy_coupled_boundary_matrix_raw(
                            omega,
                            circle,
                            1.0,
                            circle,
                            1.0,
                            beta_deg=beta_deg,
                        ),
                    }
                )
        for geometry_id in local_policy["geometries"]:
            geometry = contract["geometries"][geometry_id]
            section = comparator.rectangular_section(
                E=float(material_data["E"]),
                nu=float(material_data["nu"]),
                rho=float(material_data["rho"]),
                width=float(geometry["width"]),
                thickness=float(geometry["thickness"]),
                K=float(material_data["K"]),
            )
            cutoff = comparator.omega_cutoff(section)
            for factor in factors:
                omega = factor * cutoff
                basis = comparator.timoshenko_spatial_basis(omega, section)
                for coordinate_sign in (1.0, -1.0):
                    endpoint = comparator.clamped_endpoint_columns(
                        coordinate_sign * local_length,
                        omega,
                        section,
                        basis=basis,
                    )
                    snapshot_rows.append(
                        {
                            "kind": "local_arm",
                            "geometry": geometry_id,
                            "cutoff_factor": factor,
                            "omega": omega,
                            "coordinate_sign": coordinate_sign,
                            "state_matrix": _legacy_first_order_state_matrix(
                                omega, section
                            ),
                            "endpoint_columns": {
                                field: endpoint[field]
                                for field in ("u", "w", "psi", "N", "Q", "M")
                            },
                        }
                    )
    _assert_fresh_worker_import_boundary(worker, after_import=True)
    semantic = {
        "stage": STAGE_ID,
        "worker": worker,
        "contract_sha256": _sha256_file(contract_path),
        "rows": snapshot_rows,
    }
    semantic_sha = _sha256_bytes(
        json.dumps(
            _json_value(semantic), sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
    )
    payload = {
        **semantic,
        "generated_utc": _utc_now(),
        "fresh_process_import_boundary_passed": True,
        "semantic_snapshot_sha256": semantic_sha,
    }
    path = output_directory / f"{worker}_preliminary_snapshot.json"
    _write_json(path, payload)
    return {
        "worker": worker,
        "path": path.name,
        "contract_sha256": payload["contract_sha256"],
        "semantic_snapshot_sha256": semantic_sha,
        "file_sha256": _sha256_file(path),
        "row_count": len(snapshot_rows),
    }


def _freeze_isolated_preliminary_snapshots(
    contract_path: Path, output_directory: Path
) -> dict[str, Any]:
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    summaries: dict[str, Any] = {}
    for worker in ("rlb", "legacy"):
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--preflight-worker",
            worker,
            "--contract",
            str(contract_path),
            "--output",
            str(output_directory),
        ]
        completed = subprocess.run(
            command,
            cwd=REPOSITORY_ROOT,
            env=environment,
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"Isolated {worker} preliminary snapshot failed: {completed.stderr}"
            )
        snapshot_path = output_directory / f"{worker}_preliminary_snapshot.json"
        snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))
        if snapshot.get("contract_sha256") != _sha256_file(contract_path):
            raise RuntimeError(f"{worker} preliminary snapshot contract mismatch.")
        summaries[worker] = {
            "path": snapshot_path.name,
            "semantic_snapshot_sha256": snapshot["semantic_snapshot_sha256"],
            "file_sha256": _sha256_file(snapshot_path),
            "fresh_subprocess": True,
        }
    return summaries


def run_preliminary_gates(
    contract_path: Path, output_directory: Path
) -> dict[str, Any]:
    started = time.perf_counter()
    contract = _load_contract(contract_path)
    manifest = write_manifest_only(contract_path, output_directory)
    thresholds = contract["thresholds"]

    # Freeze each model's preliminary values in a fresh one-model subprocess
    # before this mixed audit layer is permitted to import both implementations.
    preflight_snapshots = _freeze_isolated_preliminary_snapshots(
        contract_path, output_directory
    )

    # Imports are intentionally local.  Inventory workers start in fresh
    # subprocesses and therefore never inherit this mixed audit process.
    single_beam = importlib.import_module("scripts.lib.reddy_symmetric_laminated_beam")
    coupled = importlib.import_module("scripts.lib.reddy_symmetric_coupled_beams")
    comparator = importlib.import_module(
        "scripts.lib.isotropic_rectangular_timoshenko_coupled_beams"
    )
    frozen_circular = importlib.import_module("scripts.lib.variable_length_timoshenko")

    material_contract = contract["material"]
    E = float(material_contract["E"])
    rho = float(material_contract["rho"])
    nu = float(material_contract["nu"])
    K = float(material_contract["K"])
    G = E / (2.0 * (1.0 + nu))
    material = _isotropic_material(single_beam, material_contract)
    analytic_Q = np.array(
        [
            [E / (1.0 - nu**2), nu * E / (1.0 - nu**2), 0.0],
            [nu * E / (1.0 - nu**2), E / (1.0 - nu**2), 0.0],
            [0.0, 0.0, G],
        ],
        dtype=float,
    )
    analytic_shear = np.diag([G, G])
    qbar_rows: list[dict[str, Any]] = []
    for angle in (0.0, 17.0, 45.0, 90.0, -38.0):
        qbar = single_beam.transformed_reduced_stiffness(material, angle)
        shear_qbar = single_beam.transformed_transverse_shear_stiffness(material, angle)
        q_residual = _relative_residual(qbar, analytic_Q)
        shear_residual = _relative_residual(shear_qbar, analytic_shear)
        row_pass = max(q_residual, shear_residual) <= float(
            thresholds["qbar_matrix_analytic"]
        )
        qbar_rows.append(
            {
                "angle_deg": angle,
                "nu21": material.nu21,
                "nu21_minus_nu": material.nu21 - nu,
                "Qbar_relative_residual": q_residual,
                "shear_Qbar_relative_residual": shear_residual,
                "Qbar_symmetric": bool(np.array_equal(qbar, qbar.T)),
                "Qbar_min_eigenvalue": float(np.linalg.eigvalsh(qbar)[0]),
                "status": "PASS" if row_pass else "FAIL",
            }
        )

    stack_contracts = [
        tuple(float(value) for value in contract["spectral_stack_deg"]),
        *(
            tuple(float(value) for value in stack)
            for stack in contract["algebraic_control_stacks_deg"]
        ),
    ]
    laminate_rows: list[dict[str, Any]] = []
    section_rows: list[dict[str, Any]] = []
    primary_by_geometry: dict[str, tuple[Any, Any, Any]] = {}
    for geometry_id, geometry in contract["geometries"].items():
        h = float(geometry["thickness"])
        b = float(geometry["width"])
        analytic_A = h * analytic_Q
        analytic_D = h**3 * analytic_Q / 12.0
        analytic_integrated_shear = h * analytic_shear
        analytic_I0 = rho * h
        analytic_I2 = rho * h**3 / 12.0
        analytic_properties = {
            "A": E * b * h,
            "D": E * b * h**3 / 12.0,
            "S": K * G * b * h,
            "m": rho * b * h,
            "J": rho * b * h**3 / 12.0,
        }
        for stack_index, stack in enumerate(stack_contracts):
            stack_label = _stack_label(stack)
            stack_material, laminate, properties = _four_ply_properties(
                single_beam, contract, geometry_id, stack
            )
            if stack_index == 0:
                primary_by_geometry[geometry_id] = (stack_material, laminate, properties)
            A_residual = _relative_residual(laminate.A, analytic_A)
            D_residual = _relative_residual(laminate.D, analytic_D)
            shear_residual = _relative_residual(laminate.shear, analytic_integrated_shear)
            B_residual = _scaled_zero_residual(laminate.B, np.linalg.norm(analytic_A) * h)
            I0_residual = _relative_residual(laminate.I0, analytic_I0)
            I2_residual = _relative_residual(laminate.I2, analytic_I2)
            I1_residual = abs(float(laminate.I1)) / max(abs(analytic_I0) * h, np.finfo(float).tiny)
            interface_expected = np.linspace(-0.5 * h, 0.5 * h, 5)
            interface_residual = _relative_residual(laminate.z_interfaces, interface_expected)
            matrix_max = max(
                A_residual,
                D_residual,
                shear_residual,
                B_residual,
                I0_residual,
                I2_residual,
                I1_residual,
                interface_residual,
            )
            matrix_pass = matrix_max <= float(thresholds["qbar_matrix_analytic"])
            laminate_rows.append(
                {
                    "geometry": geometry_id,
                    "stack_deg": stack_label,
                    "number_of_plies": len(laminate.plies),
                    "ply_thicknesses": [ply.thickness for ply in laminate.plies],
                    "z_interfaces": laminate.z_interfaces,
                    "A_relative_residual": A_residual,
                    "scaled_B_residual": B_residual,
                    "D_relative_residual": D_residual,
                    "shear_relative_residual": shear_residual,
                    "I0_relative_residual": I0_residual,
                    "scaled_I1_residual": I1_residual,
                    "I2_relative_residual": I2_residual,
                    "interface_relative_residual": interface_residual,
                    "maximum_residual": matrix_max,
                    "status": "PASS" if matrix_pass else "FAIL",
                }
            )
            property_residuals = {
                name: _relative_residual(getattr(properties, name), expected)
                for name, expected in analytic_properties.items()
            }
            reduction_route_max = max(
                properties.axial_reduction.relative_difference,
                properties.bending_reduction.relative_difference,
                properties.shear_reduction_before_K.relative_difference,
            )
            section_max = max(*property_residuals.values(), reduction_route_max)
            section_pass = section_max <= float(thresholds["section_reduction"])
            section_rows.append(
                {
                    "geometry": geometry_id,
                    "stack_deg": stack_label,
                    "width": b,
                    "thickness": h,
                    "area": b * h,
                    "I_y": b * h**3 / 12.0,
                    "EA_relative_residual": property_residuals["A"],
                    "EI_relative_residual": property_residuals["D"],
                    "KGA_relative_residual": property_residuals["S"],
                    "rhoA_relative_residual": property_residuals["m"],
                    "rhoI_relative_residual": property_residuals["J"],
                    "compliance_Schur_max_relative": reduction_route_max,
                    "maximum_residual": section_max,
                    "status": "PASS" if section_pass else "FAIL",
                }
            )

    circular_rows: list[dict[str, Any]] = []
    epsilon = 0.05
    old_section = frozen_circular.section_from_epsilon(epsilon)
    generic_circle = comparator.circular_section(
        E=frozen_circular.E,
        nu=frozen_circular.NU,
        rho=frozen_circular.RHO,
        radius=old_section.radius,
        K=old_section.kappa,
    )
    lambda_cutoff = frozen_circular.lambda_cutoff(epsilon, old_section)
    for cutoff_factor in (0.8, 1.0, 1.2):
        Lambda = cutoff_factor * lambda_cutoff
        omega = frozen_circular.project_omega(Lambda, epsilon)
        old_basis = frozen_circular.timo_basis(Lambda, epsilon, old_section)
        new_basis = comparator.timoshenko_spatial_basis(omega, generic_circle)
        basis_residual = max(
            _relative_residual(getattr(new_basis, name), getattr(old_basis, name))
            for name in ("a", "b", "h", "q", "z_a", "z_b", "alpha")
        )
        tolerance = float(
            thresholds[
                "legacy_circular_cutoff"
                if cutoff_factor == 1.0
                else "legacy_circular_regular"
            ]
        )
        for beta_deg in (0.0, 30.0):
            old_matrix, _warnings = frozen_circular.timo_coupling_matrix(
                Lambda, beta_deg, 0.0, epsilon, 0.0
            )
            new_matrix = comparator.legacy_coupled_boundary_matrix_raw(
                omega,
                generic_circle,
                1.0,
                generic_circle,
                1.0,
                beta_deg=beta_deg,
            )
            matrix_residual = _relative_residual(new_matrix, old_matrix)
            determinant_residual = _relative_residual(
                _row_normalized_determinant(new_matrix),
                frozen_circular.normalized_det(old_matrix),
            )
            endpoint_residual = 0.0
            for coordinate in (1.0, -1.0):
                old_endpoint = frozen_circular.endpoint_columns(
                    coordinate, omega, old_basis, old_section
                )
                new_endpoint = comparator.clamped_endpoint_columns(
                    coordinate, omega, generic_circle, basis=new_basis
                )
                endpoint_residual = max(
                    endpoint_residual,
                    *(
                        _relative_residual(new_endpoint[field], old_endpoint[field])
                        for field in ("u", "w", "psi", "N", "Q", "M")
                    ),
                )
            maximum = max(basis_residual, endpoint_residual, matrix_residual, determinant_residual)
            circular_rows.append(
                {
                    "check_kind": "basis_endpoint_matrix_backcompat",
                    "epsilon": epsilon,
                    "cutoff_factor": cutoff_factor,
                    "Lambda": Lambda,
                    "omega": omega,
                    "beta_deg": beta_deg,
                    "old_regime": old_basis.regime,
                    "new_regime": new_basis.regime,
                    "basis_relative_residual": basis_residual,
                    "endpoint_relative_residual": endpoint_residual,
                    "coupled_matrix_relative_residual": matrix_residual,
                    "normalized_determinant_relative_residual": determinant_residual,
                    "threshold": tolerance,
                    "first_six_legacy_roots": "RECORDED_IN_SEPARATE_ROOT_ROWS",
                    "maximum_residual": maximum,
                    "status": (
                        "PASS"
                        if old_basis.regime == new_basis.regime and maximum <= tolerance
                        else "FAIL"
                    ),
                }
            )
    circular_rows.extend(
        _legacy_circular_first_six_root_rows(
            comparator,
            frozen_circular,
            generic_circle,
            epsilon=epsilon,
            beta_deg=30.0,
            tolerance=float(thresholds["legacy_circular_regular"]),
        )
    )

    local_rows: list[dict[str, Any]] = []
    maps = (
        ("arm1", 1.0, np.diag([1.0, 1.0, -1.0, 1.0, 1.0, -1.0])),
        ("arm2", -1.0, np.diag([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0])),
    )
    local_policy = contract["local_solution_space_policy"]
    local_length = float(local_policy["L_local"])
    local_cutoff_factors = tuple(
        float(value) for value in local_policy["cutoff_factors"]
    )
    for geometry_id in local_policy["geometries"]:
        geometry = contract["geometries"][geometry_id]
        legacy_rectangle = comparator.rectangular_section(
            E=E,
            nu=nu,
            rho=rho,
            width=float(geometry["width"]),
            thickness=float(geometry["thickness"]),
            K=K,
        )
        _mat, _laminate, rlb_properties = primary_by_geometry[geometry_id]
        omega_cutoff = comparator.omega_cutoff(legacy_rectangle)
        state_scale = single_beam.combined_state_scale(rlb_properties, local_length)
        state_scale_inverse = np.diag(1.0 / np.diag(state_scale))
        for cutoff_factor in local_cutoff_factors:
            omega = cutoff_factor * omega_cutoff
            old_state = _legacy_first_order_state_matrix(omega, legacy_rectangle)
            rlb_state = single_beam.combined_state_matrix(omega, rlb_properties)
            basis = comparator.timoshenko_spatial_basis(omega, legacy_rectangle)
            cutoff_neighbourhood = abs(cutoff_factor - 1.0) <= 0.0100000000001
            tolerance = float(
                thresholds[
                    "local_solution_space_cutoff"
                    if cutoff_neighbourhood
                    else "local_solution_space_regular"
                ]
            )
            for arm_id, coordinate_sign, state_map in maps:
                similarity = coordinate_sign * state_map @ old_state @ state_map
                state_residual = _relative_residual(similarity, rlb_state)
                old_endpoint = comparator.clamped_endpoint_columns(
                    coordinate_sign * local_length,
                    omega,
                    legacy_rectangle,
                    basis=basis,
                )
                old_columns = np.vstack(
                    [
                        old_endpoint[field]
                        for field in ("u", "w", "psi", "N", "Q", "M")
                    ]
                )
                mapped_columns = state_map @ old_columns
                rlb_columns = coupled.arm_clamp_map(
                    omega, local_length, rlb_properties
                )
                scaled_legacy = state_scale_inverse @ mapped_columns
                scaled_rlb = state_scale_inverse @ rlb_columns
                legacy_projector, legacy_rank = _orthonormal_projector(scaled_legacy)
                rlb_projector, rlb_rank = _orthonormal_projector(scaled_rlb)
                subspace_residual = _relative_residual(
                    legacy_projector, rlb_projector
                )
                best_column_map, _residuals, column_rank, _singular = np.linalg.lstsq(
                    scaled_legacy, scaled_rlb, rcond=None
                )
                column_map_residual = _relative_residual(
                    scaled_legacy @ best_column_map, scaled_rlb
                )
                column_map_condition = float(np.linalg.cond(best_column_map))
                kinematic_map = state_map[:3, :3]
                resultant_map = state_map[3:, 3:]
                work_orientation = coordinate_sign
                work_residual = _relative_residual(
                    kinematic_map.T @ resultant_map,
                    work_orientation * np.eye(3),
                )
                maximum = max(
                    state_residual,
                    subspace_residual,
                    column_map_residual,
                    work_residual,
                )
                local_rows.append(
                    {
                        "geometry": geometry_id,
                        "arm": arm_id,
                        "local_length": local_length,
                        "cutoff_factor": cutoff_factor,
                        "cutoff_neighbourhood": cutoff_neighbourhood,
                        "omega": omega,
                        "regime": basis.regime,
                        "x_old_over_x_RLB": coordinate_sign,
                        "old_to_RLB_state_diagonal": np.diag(state_map),
                        "boundary_work_orientation_factor": work_orientation,
                        "state_matrix_similarity_relative_residual": state_residual,
                        "endpoint_solution_space_relative_residual": subspace_residual,
                        "best_column_map_relative_residual": column_map_residual,
                        "best_column_map_rank": int(column_rank),
                        "best_column_map_condition_number": column_map_condition,
                        "virtual_work_map_relative_residual": work_residual,
                        "legacy_endpoint_rank": legacy_rank,
                        "RLB_endpoint_rank": rlb_rank,
                        "threshold": tolerance,
                        "maximum_residual": maximum,
                        "status": (
                            "PASS"
                            if legacy_rank == 3
                            and rlb_rank == 3
                            and int(column_rank) == 3
                            and math.isfinite(column_map_condition)
                            and maximum <= tolerance
                            else "FAIL"
                        ),
                    }
                )

    _write_csv(output_directory / "isotropic_qbar_checks.csv", qbar_rows)
    _write_csv(output_directory / "four_ply_laminate_checks.csv", laminate_rows)
    _write_csv(output_directory / "section_reduction_checks.csv", section_rows)
    _write_csv(output_directory / "legacy_circular_backcompat.csv", circular_rows)
    _write_csv(output_directory / "local_solution_space_checks.csv", local_rows)

    statuses = {
        "RLB-ISO-4PLY-CONSTITUTIVE": (
            "PASS"
            if all(row["status"] == "PASS" for row in (*qbar_rows, *laminate_rows))
            else "FAIL"
        ),
        "RLB-ISO-SECTION-REDUCTION": (
            "PASS" if all(row["status"] == "PASS" for row in section_rows) else "FAIL"
        ),
        "RLB-ISO-LEGACY-RECTANGULAR-ADAPTER": (
            "PASS" if all(row["status"] == "PASS" for row in circular_rows) else "FAIL"
        ),
        "RLB-ISO-LOCAL-ARM-EQUIVALENCE": (
            "PASS" if all(row["status"] == "PASS" for row in local_rows) else "FAIL"
        ),
    }
    preliminary_status = "PASS" if all(value == "PASS" for value in statuses.values()) else "FAIL"
    gate_payload = {
        "stage": STAGE_ID,
        "generated_utc": _utc_now(),
        "contract_sha256": _sha256_file(contract_path),
        "gated_model_sha256": _preliminary_model_hashes(),
        "preliminary_status": preliminary_status,
        "statuses": statuses,
        "isolated_preliminary_snapshots": preflight_snapshots,
        "root_inventories_launched": False,
        "elapsed_seconds": time.perf_counter() - started,
    }
    _write_json(output_directory / "preliminary_gate.json", gate_payload)
    manifest["preliminary_statuses"] = statuses
    manifest["preliminary_gate"] = preliminary_status
    manifest["isolated_preliminary_snapshots"] = preflight_snapshots
    manifest["generated_utc"] = _utc_now()
    _write_json(output_directory / "model_manifest.json", manifest)
    _write_json(
        output_directory / "run_manifest.json",
        {
            **gate_payload,
            "workflow_state": (
                "PRELIMINARY_PASS_INVENTORIES_NOT_STARTED"
                if preliminary_status == "PASS"
                else "PRELIMINARY_FAIL_STOP_BEFORE_INVENTORIES"
            ),
            "comparison_implemented": False,
            "mode_reconstruction_implemented": False,
        },
    )
    return gate_payload


@dataclass(frozen=True)
class SearchPolicy:
    requested_roots: int
    guard_roots: int
    Omega_min: float
    Omega_max: float
    primary_scan_points: int
    verification_scan_points: int
    primary_phases: tuple[float, ...]
    verification_phases: tuple[float, ...]
    sigma_prefilter: float
    root_singular_ratio: float
    nullity_relative_threshold: float
    boundary_null_residual: float
    root_xtol_Omega: float
    root_rtol: float
    dedup_atol_Omega: float
    dedup_rtol: float
    cluster_atol_Omega: float
    cluster_rtol: float
    post_guard_tail_Omega: float
    local_close_pair_guard_subintervals: int = 0

    @classmethod
    def from_contract(cls, contract: Mapping[str, Any]) -> "SearchPolicy":
        source = contract["root_policy"]
        thresholds = contract["thresholds"]
        policy = cls(
            requested_roots=int(source["requested_roots"]),
            guard_roots=int(source["guard_roots"]),
            Omega_min=float(source["Omega_min"]),
            Omega_max=float(source["Omega_max"]),
            primary_scan_points=int(source["primary_scan_points"]),
            verification_scan_points=int(source["verification_scan_points"]),
            primary_phases=(float(source["primary_phase"]),),
            verification_phases=(float(source["verification_phase"]),),
            sigma_prefilter=float(source["sigma_prefilter"]),
            root_singular_ratio=float(thresholds["root_singular_ratio"]),
            nullity_relative_threshold=float(source["nullity_relative_threshold"]),
            boundary_null_residual=float(thresholds["boundary_null_residual"]),
            root_xtol_Omega=float(source["root_xtol_Omega"]),
            root_rtol=8.0 * np.finfo(float).eps,
            dedup_atol_Omega=float(source["dedup_atol_Omega"]),
            dedup_rtol=float(source["dedup_rtol"]),
            cluster_atol_Omega=float(source["cluster_atol_Omega"]),
            cluster_rtol=float(source["cluster_rtol"]),
            post_guard_tail_Omega=float(source["post_guard_tail_Omega"]),
            local_close_pair_guard_subintervals=int(
                source["local_close_pair_guard_subintervals"]
            ),
        )
        if policy.requested_roots != 12 or policy.guard_roots != 1:
            raise ValueError("The frozen inventory contract is 12 roots plus one guard.")
        if not 0.0 < policy.Omega_min < policy.Omega_max:
            raise ValueError("Invalid positive Omega scan interval.")
        if policy.verification_scan_points <= policy.primary_scan_points:
            raise ValueError("Verification scan must be finer than primary scan.")
        if policy.local_close_pair_guard_subintervals < 2:
            raise ValueError("The close-pair guard requires at least two subintervals.")
        return policy

    @property
    def required_slots(self) -> int:
        return self.requested_roots + self.guard_roots


@dataclass(frozen=True)
class PositiveScaling:
    scaled_matrix: FloatArray
    row_factors: FloatArray
    column_factors: FloatArray


@dataclass(frozen=True)
class MatrixDiagnostics:
    Omega: float
    omega: float
    raw_matrix: FloatArray
    scaled_matrix: FloatArray
    row_factors: FloatArray
    column_factors: FloatArray
    raw_determinant: float
    raw_det_sign: float
    raw_logabsdet: float
    scaled_determinant: float
    scaled_det_sign: float
    scaled_logabsdet: float
    raw_sigma_min: float
    raw_sigma_max: float
    raw_sigma_ratio: float
    scaled_sigma_min: float
    scaled_sigma_max: float
    scaled_sigma_ratio: float
    raw_condition_number: float
    scaled_condition_number: float
    detected_nullity: int
    root_gate_nullity: int
    scaled_null_residual: float
    raw_boundary_null_residual: float
    raw_right_null_vector: FloatArray
    finite: bool


@dataclass(frozen=True)
class RootCandidate:
    case_id: str
    builder_id: str
    scan_id: str
    Omega: float
    detection_sources: tuple[str, ...]
    interval_left_Omega: float
    interval_right_Omega: float
    interior_minimum: bool
    diagnostics: MatrixDiagnostics
    accepted: bool
    rejection_reason: str
    canonical: bool = True
    merge_group_id: str = ""


@dataclass(frozen=True)
class RootEvent:
    event_id: str
    Omega: float
    omega: float
    multiplicity: int
    detected_nullity: int
    candidate: RootCandidate
    cluster_id: str = ""
    cluster_semantics: str = "ISOLATED"
    cluster_multiplicity: int = 1
    cluster_total_nullity: int = 1
    cluster_center_Omega: float = math.nan


@dataclass(frozen=True)
class SpectrumSlot:
    sorted_slot: int
    role: str
    repeated_root_slot: int
    event: RootEvent


@dataclass(frozen=True)
class ScanResult:
    scan_id: str
    candidates: tuple[RootCandidate, ...]
    rejected_candidates: tuple[RootCandidate, ...]
    events: tuple[RootEvent, ...]
    slots: tuple[SpectrumSlot, ...]


@dataclass(frozen=True)
class RootInventory:
    case_id: str
    builder_id: str
    frequency_scale: float
    primary: ScanResult
    verification: ScanResult
    slots: tuple[SpectrumSlot, ...]
    independent_agreement: bool
    maximum_primary_verification_relative: float
    unresolved_low_sigma_count: int
    guard_available: bool
    guard_not_at_scan_boundary: bool
    status: str
    inventory_sha256: str


def _positive_equilibrate_model_neutral(matrix: FloatArray) -> PositiveScaling:
    """Positive finite row/column scaling shared only as numerical policy."""

    raw = np.asarray(matrix, dtype=float)
    if raw.ndim != 2 or raw.shape[0] != raw.shape[1]:
        raise ValueError("Boundary matrix must be square.")
    if not np.all(np.isfinite(raw)):
        raise ArithmeticError("Boundary matrix contains non-finite values.")
    row_norms = np.linalg.norm(raw, axis=1)
    row_factors = np.ones(raw.shape[0], dtype=float)
    row_reference = float(np.max(row_norms)) if row_norms.size else 0.0
    if row_reference > 0.0:
        relative_floor = float(np.finfo(float).eps ** 0.25)
        row_factors = 1.0 / np.maximum(row_norms, relative_floor * row_reference)
    row_scaled = row_factors[:, None] * raw
    column_norms = np.linalg.norm(row_scaled, axis=0)
    column_factors = 1.0 / np.maximum(column_norms, 1.0)
    scaled = row_scaled * column_factors[None, :]
    if not (
        np.all(np.isfinite(row_factors))
        and np.all(row_factors > 0.0)
        and np.all(np.isfinite(column_factors))
        and np.all(column_factors > 0.0)
        and np.all(np.isfinite(scaled))
    ):
        raise ArithmeticError("Positive equilibration failed.")
    return PositiveScaling(
        scaled_matrix=np.asarray(scaled, dtype=float),
        row_factors=np.asarray(row_factors, dtype=float),
        column_factors=np.asarray(column_factors, dtype=float),
    )


def _boundary_matrix_diagnostics(
    Omega: float,
    matrix_provider: MatrixProvider,
    frequency_scale: float,
    policy: SearchPolicy,
) -> MatrixDiagnostics:
    Omega_value = float(Omega)
    omega = Omega_value / float(frequency_scale)
    raw = np.asarray(matrix_provider(omega), dtype=float)
    scaling = _positive_equilibrate_model_neutral(raw)
    scaled = scaling.scaled_matrix
    raw_sign, raw_log = np.linalg.slogdet(raw)
    scaled_sign, scaled_log = np.linalg.slogdet(scaled)
    raw_singular = np.linalg.svd(raw, compute_uv=False)
    _left, scaled_singular, scaled_vh = np.linalg.svd(scaled, full_matrices=False)
    scaled_null = np.asarray(scaled_vh[-1], dtype=float)
    raw_null = scaling.column_factors * scaled_null
    raw_null_norm = float(np.linalg.norm(raw_null))
    scaled_null_norm = float(np.linalg.norm(scaled_null))
    if raw_null_norm == 0.0 or scaled_null_norm == 0.0:
        raise ArithmeticError("Recovered null candidate is zero.")
    raw_null = raw_null / raw_null_norm
    scaled_null = scaled_null / scaled_null_norm
    pivot = int(np.argmax(np.abs(raw_null)))
    if raw_null[pivot] < 0.0:
        raw_null = -raw_null
        scaled_null = -scaled_null
    raw_max = float(raw_singular[0])
    raw_min = float(raw_singular[-1])
    scaled_max = float(scaled_singular[0])
    scaled_min = float(scaled_singular[-1])
    scaled_residual = float(
        np.linalg.norm(scaled @ scaled_null)
        / max(scaled_max * np.linalg.norm(scaled_null), np.finfo(float).tiny)
    )
    raw_residual = float(
        np.linalg.norm(raw @ raw_null)
        / max(raw_max * np.linalg.norm(raw_null), np.finfo(float).tiny)
    )
    null_threshold = policy.nullity_relative_threshold * scaled_max
    detected_nullity = int(np.count_nonzero(scaled_singular <= null_threshold))
    root_gate_nullity = int(
        np.count_nonzero(scaled_singular <= policy.root_singular_ratio * scaled_max)
    )
    return MatrixDiagnostics(
        Omega=Omega_value,
        omega=omega,
        raw_matrix=raw,
        scaled_matrix=scaled,
        row_factors=scaling.row_factors,
        column_factors=scaling.column_factors,
        raw_determinant=float(np.linalg.det(raw)),
        raw_det_sign=float(raw_sign),
        raw_logabsdet=float(raw_log),
        scaled_determinant=float(np.linalg.det(scaled)),
        scaled_det_sign=float(scaled_sign),
        scaled_logabsdet=float(scaled_log),
        raw_sigma_min=raw_min,
        raw_sigma_max=raw_max,
        raw_sigma_ratio=raw_min / raw_max if raw_max > 0.0 else 0.0,
        scaled_sigma_min=scaled_min,
        scaled_sigma_max=scaled_max,
        scaled_sigma_ratio=scaled_min / scaled_max if scaled_max > 0.0 else 0.0,
        raw_condition_number=raw_max / raw_min if raw_min > 0.0 else math.inf,
        scaled_condition_number=scaled_max / scaled_min if scaled_min > 0.0 else math.inf,
        detected_nullity=detected_nullity,
        root_gate_nullity=root_gate_nullity,
        scaled_null_residual=scaled_residual,
        raw_boundary_null_residual=raw_residual,
        raw_right_null_vector=raw_null,
        finite=bool(
            np.all(np.isfinite(raw))
            and np.all(np.isfinite(scaled))
            and np.all(np.isfinite(scaled_singular))
        ),
    )


class _DiagnosticEvaluator:
    def __init__(
        self, provider: MatrixProvider, frequency_scale: float, policy: SearchPolicy
    ) -> None:
        self.provider = provider
        self.frequency_scale = float(frequency_scale)
        self.policy = policy
        self.cache: dict[float, MatrixDiagnostics] = {}

    def diagnostics(self, Omega: float) -> MatrixDiagnostics:
        key = float(Omega)
        if key not in self.cache:
            self.cache[key] = _boundary_matrix_diagnostics(
                key, self.provider, self.frequency_scale, self.policy
            )
        return self.cache[key]

    def determinant_objective(self, Omega: float) -> float:
        diagnostic = self.diagnostics(Omega)
        if diagnostic.scaled_det_sign == 0.0:
            return 0.0
        if not math.isfinite(diagnostic.scaled_logabsdet):
            return math.nan
        dimension = diagnostic.scaled_matrix.shape[0]
        return float(
            diagnostic.scaled_det_sign
            * math.exp(diagnostic.scaled_logabsdet / dimension)
        )

    def sigma_ratio(self, Omega: float) -> float:
        return self.diagnostics(Omega).scaled_sigma_ratio


def _candidate_quality(
    diagnostic: MatrixDiagnostics, policy: SearchPolicy
) -> tuple[bool, str]:
    if not diagnostic.finite:
        return False, "NONFINITE_MATRIX"
    if diagnostic.scaled_sigma_ratio > policy.root_singular_ratio:
        return False, "ROOT_SINGULAR_RATIO_FAIL"
    if diagnostic.scaled_null_residual > policy.boundary_null_residual:
        return False, "SCALED_NULL_RESIDUAL_FAIL"
    if diagnostic.raw_boundary_null_residual > policy.boundary_null_residual:
        return False, "RAW_BOUNDARY_NULL_RESIDUAL_FAIL"
    if diagnostic.root_gate_nullity < 1:
        return False, "ROOT_GATE_NULLITY_UNRESOLVED"
    if diagnostic.detected_nullity < 1:
        return False, "STRICT_NULLITY_UNRESOLVED"
    return True, ""


def _scan_grid(policy: SearchPolicy, points: int, phase: float) -> FloatArray:
    step = (policy.Omega_max - policy.Omega_min) / (points - 1)
    start = policy.Omega_min + float(phase) * step
    count = int(math.floor((policy.Omega_max - start) / step)) + 1
    return start + step * np.arange(count, dtype=float)


def _root_candidate(
    *,
    evaluator: _DiagnosticEvaluator,
    policy: SearchPolicy,
    case_id: str,
    builder_id: str,
    scan_id: str,
    source: str,
    left: float,
    right: float,
    Omega: float,
    interior: bool,
) -> RootCandidate:
    polished = float(Omega)
    local_half_width = max(
        1.0e-4 * (right - left),
        64.0 * policy.root_xtol_Omega,
        256.0 * abs(float(np.spacing(polished))),
    )
    delta_left = max(left - polished, -local_half_width)
    delta_right = min(right - polished, local_half_width)
    if delta_left < delta_right:
        try:
            fit = minimize_scalar(
                lambda delta: evaluator.sigma_ratio(polished + float(delta)) ** 2,
                bounds=(delta_left, delta_right),
                method="bounded",
                options={
                    "xatol": max(np.finfo(float).tiny, policy.root_xtol_Omega * 1.0e-4),
                    "maxiter": 240,
                },
            )
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            fit = None
        if fit is not None and fit.success and math.isfinite(float(fit.x)):
            fitted = polished + float(fit.x)
            if evaluator.sigma_ratio(fitted) <= evaluator.sigma_ratio(polished):
                polished = fitted
    for _iteration in range(8):
        neighbours = (
            polished,
            float(np.nextafter(polished, -math.inf)),
            float(np.nextafter(polished, math.inf)),
        )
        bounded = [value for value in neighbours if left <= value <= right]
        best = min(bounded, key=evaluator.sigma_ratio)
        if best == polished:
            break
        polished = best
    diagnostic = evaluator.diagnostics(polished)
    accepted, reason = _candidate_quality(diagnostic, policy)
    if not interior:
        accepted, reason = False, "BOUNDARY_MINIMUM"
    return RootCandidate(
        case_id=case_id,
        builder_id=builder_id,
        scan_id=scan_id,
        Omega=polished,
        detection_sources=(source,),
        interval_left_Omega=float(left),
        interval_right_Omega=float(right),
        interior_minimum=bool(interior),
        diagnostics=diagnostic,
        accepted=accepted,
        rejection_reason=reason,
    )


def _failure_candidate(
    *,
    case_id: str,
    builder_id: str,
    scan_id: str,
    source: str,
    left: float,
    right: float,
    diagnostic: MatrixDiagnostics,
    reason: str,
) -> RootCandidate:
    return RootCandidate(
        case_id=case_id,
        builder_id=builder_id,
        scan_id=scan_id,
        Omega=diagnostic.Omega,
        detection_sources=(source,),
        interval_left_Omega=float(left),
        interval_right_Omega=float(right),
        interior_minimum=True,
        diagnostics=diagnostic,
        accepted=False,
        rejection_reason=reason,
    )


def _scan_candidates(
    matrix_provider: MatrixProvider,
    frequency_scale: float,
    policy: SearchPolicy,
    *,
    case_id: str,
    builder_id: str,
    scan_id: str,
    points: int,
    phases: Sequence[float],
) -> list[RootCandidate]:
    evaluator = _DiagnosticEvaluator(matrix_provider, frequency_scale, policy)
    candidates: list[RootCandidate] = []
    for phase in phases:
        grid = _scan_grid(policy, points, float(phase))
        diagnostics = [evaluator.diagnostics(float(value)) for value in grid]
        determinants = np.asarray(
            [evaluator.determinant_objective(float(value)) for value in grid], dtype=float
        )
        sigma = np.asarray([item.scaled_sigma_ratio for item in diagnostics], dtype=float)
        source_prefix = f"{scan_id}_phase_{phase:g}"
        for index in range(grid.size - 1):
            left, right = float(grid[index]), float(grid[index + 1])
            f_left, f_right = float(determinants[index]), float(determinants[index + 1])
            if not (math.isfinite(f_left) and math.isfinite(f_right)):
                candidates.append(
                    _failure_candidate(
                        case_id=case_id,
                        builder_id=builder_id,
                        scan_id=scan_id,
                        source=source_prefix + ":determinant_interval",
                        left=left,
                        right=right,
                        diagnostic=diagnostics[index],
                        reason="NONFINITE_DETERMINANT_INTERVAL",
                    )
                )
                continue
            if f_left == 0.0:
                candidates.append(
                    _root_candidate(
                        evaluator=evaluator,
                        policy=policy,
                        case_id=case_id,
                        builder_id=builder_id,
                        scan_id=scan_id,
                        source=source_prefix + ":determinant_grid_zero",
                        left=left,
                        right=right,
                        Omega=left,
                        interior=left > policy.Omega_min + policy.root_xtol_Omega,
                    )
                )
            elif f_left * f_right < 0.0:
                try:
                    root = float(
                        brentq(
                            evaluator.determinant_objective,
                            left,
                            right,
                            xtol=policy.root_xtol_Omega,
                            rtol=policy.root_rtol,
                            maxiter=180,
                        )
                    )
                except (RuntimeError, ValueError, FloatingPointError):
                    candidates.append(
                        _failure_candidate(
                            case_id=case_id,
                            builder_id=builder_id,
                            scan_id=scan_id,
                            source=source_prefix + ":determinant_bracket",
                            left=left,
                            right=right,
                            diagnostic=diagnostics[index],
                            reason="BRENT_FAILURE",
                        )
                    )
                else:
                    candidates.append(
                        _root_candidate(
                            evaluator=evaluator,
                            policy=policy,
                            case_id=case_id,
                            builder_id=builder_id,
                            scan_id=scan_id,
                            source=source_prefix + ":determinant_bracket",
                            left=left,
                            right=right,
                            Omega=root,
                            interior=True,
                        )
                    )
        for index in range(1, grid.size - 1):
            current = float(sigma[index])
            if not math.isfinite(current):
                continue
            if current > float(sigma[index - 1]) or current > float(sigma[index + 1]):
                continue
            left, right = float(grid[index - 1]), float(grid[index + 1])
            if policy.local_close_pair_guard_subintervals > 0:
                local_grid = np.linspace(
                    left,
                    right,
                    policy.local_close_pair_guard_subintervals + 1,
                    dtype=float,
                )
                local_determinants = [
                    evaluator.determinant_objective(float(value))
                    for value in local_grid
                ]
                for local_index in range(local_grid.size - 1):
                    local_left = float(local_grid[local_index])
                    local_right = float(local_grid[local_index + 1])
                    local_f_left = float(local_determinants[local_index])
                    local_f_right = float(local_determinants[local_index + 1])
                    local_source = source_prefix + ":determinant_close_pair_guard"
                    if not (
                        math.isfinite(local_f_left)
                        and math.isfinite(local_f_right)
                    ):
                        candidates.append(
                            _failure_candidate(
                                case_id=case_id,
                                builder_id=builder_id,
                                scan_id=scan_id,
                                source=local_source,
                                left=local_left,
                                right=local_right,
                                diagnostic=evaluator.diagnostics(local_left),
                                reason="NONFINITE_DETERMINANT_INTERVAL",
                            )
                        )
                        continue
                    if local_f_left == 0.0:
                        candidates.append(
                            _root_candidate(
                                evaluator=evaluator,
                                policy=policy,
                                case_id=case_id,
                                builder_id=builder_id,
                                scan_id=scan_id,
                                source=local_source,
                                left=local_left,
                                right=local_right,
                                Omega=local_left,
                                interior=(
                                    local_left
                                    > policy.Omega_min + policy.root_xtol_Omega
                                ),
                            )
                        )
                    elif local_f_left * local_f_right < 0.0:
                        try:
                            local_root = float(
                                brentq(
                                    evaluator.determinant_objective,
                                    local_left,
                                    local_right,
                                    xtol=policy.root_xtol_Omega,
                                    rtol=policy.root_rtol,
                                    maxiter=180,
                                )
                            )
                        except (RuntimeError, ValueError, FloatingPointError):
                            candidates.append(
                                _failure_candidate(
                                    case_id=case_id,
                                    builder_id=builder_id,
                                    scan_id=scan_id,
                                    source=local_source,
                                    left=local_left,
                                    right=local_right,
                                    diagnostic=evaluator.diagnostics(local_left),
                                    reason="BRENT_FAILURE",
                                )
                            )
                        else:
                            candidates.append(
                                _root_candidate(
                                    evaluator=evaluator,
                                    policy=policy,
                                    case_id=case_id,
                                    builder_id=builder_id,
                                    scan_id=scan_id,
                                    source=local_source,
                                    left=local_left,
                                    right=local_right,
                                    Omega=local_root,
                                    interior=True,
                                )
                            )
            try:
                first = minimize_scalar(
                    lambda value: evaluator.sigma_ratio(float(value)) ** 2,
                    bounds=(left, right),
                    method="bounded",
                    options={"xatol": policy.root_xtol_Omega, "maxiter": 180},
                )
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                candidates.append(
                    _failure_candidate(
                        case_id=case_id,
                        builder_id=builder_id,
                        scan_id=scan_id,
                        source=source_prefix + ":sigma_ratio_minimum",
                        left=left,
                        right=right,
                        diagnostic=diagnostics[index],
                        reason="MINIMIZER_EXCEPTION",
                    )
                )
                continue
            if not first.success or not math.isfinite(float(first.x)):
                candidates.append(
                    _failure_candidate(
                        case_id=case_id,
                        builder_id=builder_id,
                        scan_id=scan_id,
                        source=source_prefix + ":sigma_ratio_minimum",
                        left=left,
                        right=right,
                        diagnostic=diagnostics[index],
                        reason="MINIMIZER_FAILURE",
                    )
                )
                continue
            root = float(first.x)
            half = max((right - left) / 8.0, 8.0 * policy.root_xtol_Omega)
            narrow_left, narrow_right = max(left, root - half), min(right, root + half)
            try:
                second = minimize_scalar(
                    lambda delta: evaluator.sigma_ratio(root + float(delta)) ** 2,
                    bounds=(narrow_left - root, narrow_right - root),
                    method="bounded",
                    options={"xatol": policy.root_xtol_Omega, "maxiter": 180},
                )
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                candidates.append(
                    _failure_candidate(
                        case_id=case_id,
                        builder_id=builder_id,
                        scan_id=scan_id,
                        source=source_prefix + ":sigma_ratio_minimum",
                        left=left,
                        right=right,
                        diagnostic=evaluator.diagnostics(root),
                        reason="NESTED_MINIMIZER_EXCEPTION",
                    )
                )
                continue
            if second.success and math.isfinite(float(second.fun)) and float(second.fun) <= float(first.fun):
                root += float(second.x)
            edge = max(4.0 * policy.root_xtol_Omega, 1.0e-5 * (right - left))
            candidate = _root_candidate(
                evaluator=evaluator,
                policy=policy,
                case_id=case_id,
                builder_id=builder_id,
                scan_id=scan_id,
                source=source_prefix + ":sigma_ratio_minimum",
                left=left,
                right=right,
                Omega=root,
                interior=left + edge < root < right - edge,
            )
            if not candidate.accepted and candidate.rejection_reason == "ROOT_SINGULAR_RATIO_FAIL":
                candidate = replace(candidate, rejection_reason="FALSE_SIGMA_VALLEY")
            candidates.append(candidate)
    return candidates


def _candidate_close(
    left: RootCandidate, right: RootCandidate, policy: SearchPolicy
) -> bool:
    tolerance = policy.dedup_atol_Omega + policy.dedup_rtol * max(
        abs(left.Omega), abs(right.Omega)
    )
    return abs(left.Omega - right.Omega) <= tolerance


def _same_root_detection(
    left: RootCandidate, right: RootCandidate, policy: SearchPolicy
) -> bool:
    """Recognize duplicate detectors without broadening the frozen root tolerance.

    A determinant bracket and a sigma-minimum search can converge to slightly
    different floating-point representatives when the boundary matrix is very
    ill-conditioned, even though both recover the same one-dimensional null
    vector.  The ordinary frozen Omega tolerance remains the primary rule.  A
    second rule is allowed only for candidates whose search intervals overlap,
    contain both estimates, and whose normalized raw null vectors are parallel
    to roundoff.  Distinct close roots with independent null directions
    therefore remain distinct events, including roots recovered by the local
    close-pair guard.
    """

    if _candidate_close(left, right, policy):
        return True
    overlap_left = max(left.interval_left_Omega, right.interval_left_Omega)
    overlap_right = min(left.interval_right_Omega, right.interval_right_Omega)
    if overlap_left > overlap_right:
        return False
    if not (
        overlap_left <= left.Omega <= overlap_right
        and overlap_left <= right.Omega <= overlap_right
    ):
        return False
    null_overlap = abs(
        float(
            np.dot(
                left.diagnostics.raw_right_null_vector,
                right.diagnostics.raw_right_null_vector,
            )
        )
    )
    return 1.0 - min(null_overlap, 1.0) <= 256.0 * np.finfo(float).eps


def _merge_candidates(
    candidates: Sequence[RootCandidate], policy: SearchPolicy
) -> tuple[list[RootCandidate], list[RootCandidate]]:
    accepted = sorted(
        (candidate for candidate in candidates if candidate.accepted),
        key=lambda item: (item.Omega, item.diagnostics.scaled_sigma_ratio),
    )
    canonical: list[RootCandidate] = []
    duplicates: list[RootCandidate] = []
    for candidate in accepted:
        match_index = next(
            (
                index
                for index, previous in enumerate(canonical)
                if _same_root_detection(previous, candidate, policy)
            ),
            None,
        )
        if match_index is None:
            canonical.append(candidate)
            continue
        previous = canonical[match_index]
        best, duplicate = (
            (candidate, previous)
            if candidate.diagnostics.scaled_sigma_ratio
            < previous.diagnostics.scaled_sigma_ratio
            else (previous, candidate)
        )
        merge_id = f"merge_{match_index + 1:04d}"
        canonical[match_index] = replace(
            best,
            detection_sources=tuple(
                sorted(set(previous.detection_sources + candidate.detection_sources))
            ),
            merge_group_id=merge_id,
            canonical=True,
        )
        duplicates.append(
            replace(
                duplicate,
                accepted=False,
                rejection_reason="DUPLICATE_DETECTION_SAME_ROOT",
                canonical=False,
                merge_group_id=merge_id,
            )
        )
    rejected = [candidate for candidate in candidates if not candidate.accepted]
    rejected.extend(duplicates)
    canonical.sort(key=lambda item: item.Omega)
    return canonical, rejected


def _events_and_slots(
    canonical: Sequence[RootCandidate], policy: SearchPolicy
) -> tuple[tuple[RootEvent, ...], tuple[SpectrumSlot, ...]]:
    raw_events = [
        RootEvent(
            event_id=f"event_{index:04d}",
            Omega=candidate.Omega,
            omega=candidate.diagnostics.omega,
            multiplicity=candidate.diagnostics.detected_nullity,
            detected_nullity=candidate.diagnostics.detected_nullity,
            candidate=candidate,
            cluster_center_Omega=candidate.Omega,
        )
        for index, candidate in enumerate(canonical, start=1)
    ]
    groups: list[list[int]] = []
    for index, event in enumerate(raw_events):
        if not groups:
            groups.append([index])
            continue
        previous = raw_events[groups[-1][-1]]
        tolerance = policy.cluster_atol_Omega + policy.cluster_rtol * max(
            abs(previous.Omega), abs(event.Omega)
        )
        if event.Omega - previous.Omega <= tolerance:
            groups[-1].append(index)
        else:
            groups.append([index])
    events = list(raw_events)
    cluster_counter = 0
    for indices in groups:
        multiplicity = sum(events[index].multiplicity for index in indices)
        if multiplicity <= 1:
            continue
        cluster_counter += 1
        cluster_id = f"cluster_{cluster_counter:04d}"
        centre = math.fsum(
            events[index].Omega * events[index].multiplicity for index in indices
        ) / multiplicity
        total_nullity = sum(events[index].detected_nullity for index in indices)
        exact = len(indices) == 1 and events[indices[0]].multiplicity > 1
        semantics = "EXACT_DEGENERATE_SUBSPACE" if exact else "NEAR_DEGENERATE_CLUSTER"
        for index in indices:
            events[index] = replace(
                events[index],
                cluster_id=cluster_id,
                cluster_semantics=semantics,
                cluster_multiplicity=multiplicity,
                cluster_total_nullity=total_nullity,
                cluster_center_Omega=centre,
            )
    slots: list[SpectrumSlot] = []
    for event in events:
        for repeated in range(1, event.multiplicity + 1):
            slot_number = len(slots) + 1
            role = "FIRST_12" if slot_number <= policy.requested_roots else "ROOT_13_GUARD"
            slots.append(
                SpectrumSlot(
                    sorted_slot=slot_number,
                    role=role,
                    repeated_root_slot=repeated,
                    event=event,
                )
            )
    retained = slots[: policy.required_slots]
    if retained and retained[-1].event.cluster_id:
        final_cluster = retained[-1].event.cluster_id
        retained = [
            replace(
                slot,
                role=(
                    "GUARD_CLUSTER_COMPLETION"
                    if slot.sorted_slot > policy.required_slots
                    else slot.role
                ),
            )
            for slot in slots
            if slot.sorted_slot <= policy.required_slots or slot.event.cluster_id == final_cluster
        ]
    return tuple(events), tuple(retained)


def _run_scan(
    matrix_provider: MatrixProvider,
    frequency_scale: float,
    policy: SearchPolicy,
    *,
    case_id: str,
    builder_id: str,
    scan_id: str,
    points: int,
    phases: Sequence[float],
) -> ScanResult:
    raw = _scan_candidates(
        matrix_provider,
        frequency_scale,
        policy,
        case_id=case_id,
        builder_id=builder_id,
        scan_id=scan_id,
        points=points,
        phases=phases,
    )
    canonical, rejected = _merge_candidates(raw, policy)
    events, slots = _events_and_slots(canonical, policy)
    return ScanResult(
        scan_id=scan_id,
        candidates=tuple(canonical),
        rejected_candidates=tuple(rejected),
        events=events,
        slots=slots,
    )


def _relative_difference(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / max(
        abs(float(left)), abs(float(right)), np.finfo(float).tiny
    )


def _inventory_digest_payload(slots: Sequence[SpectrumSlot]) -> list[dict[str, Any]]:
    return [
        {
            "sorted_slot": slot.sorted_slot,
            "Omega": format(slot.event.Omega, ".17g"),
            "omega": format(slot.event.omega, ".17g"),
            "multiplicity": slot.event.multiplicity,
            "nullity": slot.event.detected_nullity,
            "cluster_id": slot.event.cluster_id,
        }
        for slot in slots
    ]


def _digest_float(value: float) -> str:
    return format(float(value), ".17g")


def _candidate_digest_payload(candidate: RootCandidate) -> dict[str, Any]:
    diagnostic = candidate.diagnostics
    return {
        "Omega": _digest_float(candidate.Omega),
        "scan_id": candidate.scan_id,
        "detection_sources": list(candidate.detection_sources),
        "interval_left_Omega": _digest_float(candidate.interval_left_Omega),
        "interval_right_Omega": _digest_float(candidate.interval_right_Omega),
        "interior_minimum": candidate.interior_minimum,
        "accepted": candidate.accepted,
        "rejection_reason": candidate.rejection_reason,
        "canonical": candidate.canonical,
        "merge_group_id": candidate.merge_group_id,
        "raw_det_sign": _digest_float(diagnostic.raw_det_sign),
        "raw_logabsdet": _digest_float(diagnostic.raw_logabsdet),
        "scaled_det_sign": _digest_float(diagnostic.scaled_det_sign),
        "scaled_logabsdet": _digest_float(diagnostic.scaled_logabsdet),
        "raw_sigma_ratio": _digest_float(diagnostic.raw_sigma_ratio),
        "scaled_sigma_ratio": _digest_float(diagnostic.scaled_sigma_ratio),
        "detected_nullity": diagnostic.detected_nullity,
        "root_gate_nullity": diagnostic.root_gate_nullity,
        "scaled_null_residual": _digest_float(diagnostic.scaled_null_residual),
        "raw_boundary_null_residual": _digest_float(
            diagnostic.raw_boundary_null_residual
        ),
    }


def _scan_digest_payload(scan: ScanResult) -> dict[str, Any]:
    return {
        "scan_id": scan.scan_id,
        "candidates": [_candidate_digest_payload(item) for item in scan.candidates],
        "rejected_candidates": [
            _candidate_digest_payload(item) for item in scan.rejected_candidates
        ],
        "events": [
            {
                "event_id": item.event_id,
                "Omega": _digest_float(item.Omega),
                "omega": _digest_float(item.omega),
                "multiplicity": item.multiplicity,
                "detected_nullity": item.detected_nullity,
                "cluster_id": item.cluster_id,
                "cluster_semantics": item.cluster_semantics,
                "cluster_multiplicity": item.cluster_multiplicity,
                "cluster_total_nullity": item.cluster_total_nullity,
                "cluster_center_Omega": _digest_float(item.cluster_center_Omega),
            }
            for item in scan.events
        ],
        "slots": _inventory_digest_payload(scan.slots),
    }


def _root_inventory_digest_payload(
    *,
    case_id: str,
    builder_id: str,
    contract_sha256: str,
    frequency_scale: float,
    policy: SearchPolicy,
    primary: ScanResult,
    verification: ScanResult,
    independent_agreement: bool,
    maximum_primary_verification_relative: float,
    unresolved_low_sigma_count: int,
    guard_available: bool,
    guard_not_at_scan_boundary: bool,
    status: str,
) -> dict[str, Any]:
    return {
        "case_id": case_id,
        "builder_id": builder_id,
        "contract_sha256": contract_sha256,
        "frequency_scale": _digest_float(frequency_scale),
        "policy": _json_value(policy.__dict__),
        "primary": _scan_digest_payload(primary),
        "verification": _scan_digest_payload(verification),
        "certification": {
            "independent_agreement": independent_agreement,
            "maximum_primary_verification_relative": _digest_float(
                maximum_primary_verification_relative
            ),
            "unresolved_low_sigma_count": unresolved_low_sigma_count,
            "guard_available": guard_available,
            "guard_not_at_scan_boundary": guard_not_at_scan_boundary,
            "status": status,
        },
    }


def seed_free_root_inventory(
    matrix_provider: MatrixProvider,
    frequency_scale: float,
    policy: SearchPolicy,
    *,
    case_id: str,
    builder_id: str,
    contract_sha256: str = "",
) -> RootInventory:
    primary = _run_scan(
        matrix_provider,
        frequency_scale,
        policy,
        case_id=case_id,
        builder_id=builder_id,
        scan_id="primary",
        points=policy.primary_scan_points,
        phases=policy.primary_phases,
    )
    verification = _run_scan(
        matrix_provider,
        frequency_scale,
        policy,
        case_id=case_id,
        builder_id=builder_id,
        scan_id="verification",
        points=policy.verification_scan_points,
        phases=policy.verification_phases,
    )
    required = policy.required_slots
    primary_slots = primary.slots[:required]
    verification_slots = verification.slots[:required]
    agreement = len(primary_slots) >= required and len(verification_slots) >= required
    maximum_relative = 0.0
    if agreement:
        for left, right in zip(primary_slots, verification_slots, strict=True):
            left_clustered = bool(left.event.cluster_id)
            right_clustered = bool(right.event.cluster_id)
            comparison_left = (
                left.event.cluster_center_Omega if left_clustered else left.event.Omega
            )
            comparison_right = (
                right.event.cluster_center_Omega if right_clustered else right.event.Omega
            )
            relative = _relative_difference(comparison_left, comparison_right)
            maximum_relative = max(maximum_relative, relative)
            if (
                relative > 1.0e-8
                or left.event.multiplicity != right.event.multiplicity
                or left.event.detected_nullity != right.event.detected_nullity
                or left_clustered != right_clustered
                or left.event.cluster_multiplicity != right.event.cluster_multiplicity
                or left.event.cluster_total_nullity != right.event.cluster_total_nullity
            ):
                agreement = False
    guard_available = len(primary_slots) >= required
    guard_not_boundary = bool(
        guard_available
        and primary_slots[required - 1].event.Omega
        <= policy.Omega_max - policy.post_guard_tail_Omega
    )
    guard_frequency = (
        primary_slots[required - 1].event.Omega if guard_available else policy.Omega_max
    )
    accepted_events = (*primary.events, *verification.events)

    def genuinely_unresolved(candidate: RootCandidate) -> bool:
        if candidate.interval_left_Omega > guard_frequency:
            return False
        if candidate.rejection_reason == "DUPLICATE_DETECTION_SAME_ROOT":
            return False
        if any(
            abs(candidate.Omega - event.Omega)
            <= policy.dedup_atol_Omega
            + policy.dedup_rtol * max(abs(candidate.Omega), abs(event.Omega))
            for event in accepted_events
        ):
            return False
        if candidate.rejection_reason in {
            "NONFINITE_MATRIX",
            "NONFINITE_DETERMINANT_INTERVAL",
            "BRENT_FAILURE",
            "MINIMIZER_EXCEPTION",
            "MINIMIZER_FAILURE",
            "NESTED_MINIMIZER_EXCEPTION",
        }:
            return True
        return candidate.diagnostics.scaled_sigma_ratio <= policy.sigma_prefilter

    unresolved = sum(
        genuinely_unresolved(candidate)
        for candidate in (*primary.rejected_candidates, *verification.rejected_candidates)
    )
    if not guard_available or len(verification_slots) < required or not agreement:
        status = "FAIL"
    elif guard_not_boundary and unresolved == 0:
        status = "PASS"
    else:
        status = "PARTIAL_PASS"
    semantic_payload = _root_inventory_digest_payload(
        case_id=case_id,
        builder_id=builder_id,
        contract_sha256=contract_sha256,
        frequency_scale=frequency_scale,
        policy=policy,
        primary=primary,
        verification=verification,
        independent_agreement=agreement,
        maximum_primary_verification_relative=maximum_relative,
        unresolved_low_sigma_count=int(unresolved),
        guard_available=guard_available,
        guard_not_at_scan_boundary=guard_not_boundary,
        status=status,
    )
    digest = _sha256_bytes(
        json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode(
            "utf-8"
        )
    )
    return RootInventory(
        case_id=case_id,
        builder_id=builder_id,
        frequency_scale=float(frequency_scale),
        primary=primary,
        verification=verification,
        slots=primary.slots,
        independent_agreement=agreement,
        maximum_primary_verification_relative=maximum_relative,
        unresolved_low_sigma_count=int(unresolved),
        guard_available=guard_available,
        guard_not_at_scan_boundary=guard_not_boundary,
        status=status,
        inventory_sha256=digest,
    )


def _candidate_row(candidate: RootCandidate) -> dict[str, Any]:
    diagnostic = candidate.diagnostics
    return {
        "case_id": candidate.case_id,
        "builder_id": candidate.builder_id,
        "scan_id": candidate.scan_id,
        "Omega": candidate.Omega,
        "omega": diagnostic.omega,
        "detection_sources": candidate.detection_sources,
        "interval_left_Omega": candidate.interval_left_Omega,
        "interval_right_Omega": candidate.interval_right_Omega,
        "interior_minimum": candidate.interior_minimum,
        "accepted": candidate.accepted,
        "rejection_reason": candidate.rejection_reason,
        "canonical": candidate.canonical,
        "merge_group_id": candidate.merge_group_id,
        "raw_determinant": diagnostic.raw_determinant,
        "raw_det_sign": diagnostic.raw_det_sign,
        "raw_logabsdet": diagnostic.raw_logabsdet,
        "scaled_determinant": diagnostic.scaled_determinant,
        "scaled_det_sign": diagnostic.scaled_det_sign,
        "scaled_logabsdet": diagnostic.scaled_logabsdet,
        "raw_sigma_min": diagnostic.raw_sigma_min,
        "raw_sigma_max": diagnostic.raw_sigma_max,
        "raw_sigma_ratio": diagnostic.raw_sigma_ratio,
        "scaled_sigma_min": diagnostic.scaled_sigma_min,
        "scaled_sigma_max": diagnostic.scaled_sigma_max,
        "scaled_sigma_ratio": diagnostic.scaled_sigma_ratio,
        "raw_condition_number": diagnostic.raw_condition_number,
        "scaled_condition_number": diagnostic.scaled_condition_number,
        "detected_nullity": diagnostic.detected_nullity,
        "root_gate_nullity": diagnostic.root_gate_nullity,
        "scaled_null_residual": diagnostic.scaled_null_residual,
        "raw_boundary_null_residual": diagnostic.raw_boundary_null_residual,
        "row_factors": diagnostic.row_factors,
        "column_factors": diagnostic.column_factors,
    }


def _root_row(inventory: RootInventory, slot: SpectrumSlot) -> dict[str, Any]:
    event = slot.event
    diagnostic = event.candidate.diagnostics
    return {
        "case_id": inventory.case_id,
        "builder_id": inventory.builder_id,
        "sorted_slot": slot.sorted_slot,
        "role": slot.role,
        "repeated_root_slot": slot.repeated_root_slot,
        "event_id": event.event_id,
        "omega": event.omega,
        "Omega": event.Omega,
        "multiplicity": event.multiplicity,
        "detected_nullity": event.detected_nullity,
        "root_gate_nullity": diagnostic.root_gate_nullity,
        "cluster_id": event.cluster_id,
        "cluster_semantics": event.cluster_semantics,
        "cluster_multiplicity": event.cluster_multiplicity,
        "cluster_total_nullity": event.cluster_total_nullity,
        "cluster_center_Omega": event.cluster_center_Omega,
        "raw_determinant": diagnostic.raw_determinant,
        "raw_det_sign": diagnostic.raw_det_sign,
        "raw_logabsdet": diagnostic.raw_logabsdet,
        "scaled_determinant": diagnostic.scaled_determinant,
        "scaled_det_sign": diagnostic.scaled_det_sign,
        "scaled_logabsdet": diagnostic.scaled_logabsdet,
        "raw_sigma_min": diagnostic.raw_sigma_min,
        "raw_sigma_max": diagnostic.raw_sigma_max,
        "raw_sigma_ratio": diagnostic.raw_sigma_ratio,
        "scaled_sigma_min": diagnostic.scaled_sigma_min,
        "scaled_sigma_max": diagnostic.scaled_sigma_max,
        "scaled_sigma_ratio": diagnostic.scaled_sigma_ratio,
        "raw_condition_number": diagnostic.raw_condition_number,
        "scaled_condition_number": diagnostic.scaled_condition_number,
        "scaled_null_residual": diagnostic.scaled_null_residual,
        "raw_boundary_null_residual": diagnostic.raw_boundary_null_residual,
        "row_factors": diagnostic.row_factors,
        "column_factors": diagnostic.column_factors,
        "detection_sources": event.candidate.detection_sources,
        "bracket_left_Omega": event.candidate.interval_left_Omega,
        "bracket_right_Omega": event.candidate.interval_right_Omega,
        "inventory_status": inventory.status,
        "case_inventory_sha256": inventory.inventory_sha256,
    }


def _frequency_scale(contract: Mapping[str, Any], geometry_id: str) -> float:
    material = contract["material"]
    geometry = contract["geometries"][geometry_id]
    E = float(material["E"])
    rho = float(material["rho"])
    b = float(geometry["width"])
    h = float(geometry["thickness"])
    area = b * h
    inertia = b * h**3 / 12.0
    L_ref = float(contract["lengths"]["L_ref"])
    return L_ref**2 * math.sqrt(rho * area / (E * inertia))


def _scan_manifest_summary(scan: ScanResult) -> dict[str, Any]:
    return {
        "scan_id": scan.scan_id,
        "candidate_count": len(scan.candidates),
        "rejected_candidate_count": len(scan.rejected_candidates),
        "event_count": len(scan.events),
        "slot_count_including_guard_cluster_completion": len(scan.slots),
        "strict_nullities": [event.detected_nullity for event in scan.events],
        "roots": [
            {
                "event_id": event.event_id,
                "Omega": event.Omega,
                "omega": event.omega,
                "multiplicity": event.multiplicity,
                "detected_nullity": event.detected_nullity,
                "cluster_id": event.cluster_id,
                "cluster_semantics": event.cluster_semantics,
                "cluster_multiplicity": event.cluster_multiplicity,
                "cluster_total_nullity": event.cluster_total_nullity,
                "scaled_sigma_ratio": event.candidate.diagnostics.scaled_sigma_ratio,
                "raw_boundary_null_residual": (
                    event.candidate.diagnostics.raw_boundary_null_residual
                ),
            }
            for event in scan.events
        ],
    }


def _assert_fresh_worker_import_boundary(worker: str, *, after_import: bool) -> None:
    if worker == "rlb":
        forbidden = (
            "scripts.lib.isotropic_rectangular_timoshenko_coupled_beams",
            "scripts.lib.variable_length_timoshenko",
        )
    else:
        forbidden = (
            "scripts.lib.reddy_symmetric_laminated_beam",
            "scripts.lib.reddy_inplane_geometry",
            "scripts.lib.reddy_symmetric_coupled_beams",
            "scripts.lib.reddy_symmetric_coupled_beams_ritz",
        )
    loaded = [name for name in forbidden if name in sys.modules]
    if loaded:
        stage = "after model import" if after_import else "before model import"
        raise RuntimeError(f"{worker} worker isolation failed {stage}: {loaded}")


def _worker_case_provider(
    worker: str,
    contract: Mapping[str, Any],
    case: Mapping[str, Any],
    *,
    rlb_arm_cache: dict[tuple[str, float, float], FloatArray] | None = None,
) -> tuple[MatrixProvider, dict[str, Any]]:
    geometry_id = str(case["geometry"])
    geometry = contract["geometries"][geometry_id]
    material = contract["material"]
    length_1 = float(case["L1"])
    length_2 = float(case["L2"])
    beta_deg = float(case["beta_deg"])
    if worker == "rlb":
        single_beam = importlib.import_module(
            "scripts.lib.reddy_symmetric_laminated_beam"
        )
        coupled = importlib.import_module("scripts.lib.reddy_symmetric_coupled_beams")
        _material, laminate, properties = _four_ply_properties(
            single_beam,
            contract,
            geometry_id,
            tuple(float(value) for value in contract["spectral_stack_deg"]),
        )
        if len(laminate.plies) != 4:
            raise RuntimeError("RLB worker must use exactly four plies.")

        cache = rlb_arm_cache if rlb_arm_cache is not None else {}
        beta_rad = math.radians(beta_deg)
        frozen_joint_matrix = np.asarray(coupled.joint_matrix(beta_rad), dtype=float)

        def cached_clamp_map(omega: float, length: float) -> FloatArray:
            key = (geometry_id, float(length), float(omega))
            if key not in cache:
                cache[key] = np.asarray(
                    coupled.arm_clamp_map(float(omega), float(length), properties),
                    dtype=float,
                )
            return cache[key]

        def provider(omega: float) -> FloatArray:
            combined = np.zeros((12, 6), dtype=float)
            combined[:6, :3] = cached_clamp_map(float(omega), length_1)
            combined[6:, 3:] = cached_clamp_map(float(omega), length_2)
            return np.asarray(frozen_joint_matrix @ combined, dtype=float)

        cache_equivalence = 0.0
        for omega_check in (0.0, 0.731, 3.217):
            direct = np.asarray(
                coupled.coupled_boundary_matrix(
                    omega_check,
                    beta_rad,
                    length_1,
                    properties,
                    length_2,
                    properties,
                ),
                dtype=float,
            )
            cache_equivalence = max(
                cache_equivalence,
                float(np.max(np.abs(provider(omega_check) - direct))),
            )
        if cache_equivalence > 16.0 * np.finfo(float).eps:
            raise RuntimeError("Cached RLB arm assembly differs from the frozen builder.")

        metadata = {
            "model_module": "scripts.lib.reddy_symmetric_coupled_beams",
            "arm_module": "scripts.lib.reddy_symmetric_laminated_beam",
            "number_of_plies": len(laminate.plies),
            "stack_deg": [ply.angle_deg for ply in laminate.plies],
            "cached_arm_assembly": True,
            "cached_vs_frozen_builder_max_abs": cache_equivalence,
            "beam_properties": {
                "A": properties.A,
                "D": properties.D,
                "S": properties.S,
                "m": properties.m,
                "J": properties.J,
                "K": properties.K,
            },
        }
        return provider, metadata

    comparator = importlib.import_module(
        "scripts.lib.isotropic_rectangular_timoshenko_coupled_beams"
    )
    section = comparator.rectangular_section(
        E=float(material["E"]),
        nu=float(material["nu"]),
        rho=float(material["rho"]),
        width=float(geometry["width"]),
        thickness=float(geometry["thickness"]),
        K=float(material["K"]),
    )

    def provider(omega: float) -> FloatArray:
        return np.asarray(
            comparator.legacy_coupled_boundary_matrix_raw(
                float(omega),
                section,
                length_1,
                section,
                length_2,
                beta_deg=beta_deg,
            ),
            dtype=float,
        )

    metadata = {
        "model_module": "scripts.lib.isotropic_rectangular_timoshenko_coupled_beams",
        "uses_matrix_exponential": False,
        "section_properties": {
            "E": section.E,
            "nu": section.nu,
            "rho": section.rho,
            "width": section.width,
            "thickness": section.thickness,
            "area": section.area,
            "I_y": section.I_y,
            "EA": section.EA,
            "EI": section.EI,
            "KGA": section.KGA,
            "rhoA": section.rhoA,
            "rhoI": section.rhoI,
            "K": section.K,
        },
    }
    return provider, metadata


def run_inventory_worker(
    worker: str, contract_path: Path, output_directory: Path
) -> dict[str, Any]:
    if worker not in {"rlb", "legacy"}:
        raise ValueError("worker must be 'rlb' or 'legacy'.")
    _assert_fresh_worker_import_boundary(worker, after_import=False)
    contract = _load_contract(contract_path)
    gate_path = output_directory / "preliminary_gate.json"
    if not gate_path.is_file():
        raise RuntimeError("Preliminary gate is absent; refusing spectral search.")
    gate = json.loads(gate_path.read_text(encoding="utf-8"))
    if gate.get("preliminary_status") != "PASS":
        raise RuntimeError("Preliminary gate is not PASS; refusing spectral search.")
    if gate.get("contract_sha256") != _sha256_file(contract_path):
        raise RuntimeError("Preliminary gate and worker contract hashes differ.")
    if gate.get("gated_model_sha256") != _preliminary_model_hashes():
        raise RuntimeError(
            "Scientific model files changed after the preliminary PASS gate."
        )

    started = time.perf_counter()
    policy = SearchPolicy.from_contract(contract)
    builder_id = "RLB_FOUR_PLY_TRANSFER" if worker == "rlb" else "LEGACY_RECT_CLOSED_FORM"
    prefix = "rlb" if worker == "rlb" else "legacy_rectangular"
    inventories: list[RootInventory] = []
    case_metadata: list[dict[str, Any]] = []
    rlb_arm_cache: dict[tuple[str, float, float], FloatArray] = {}
    contract_sha256 = _sha256_file(contract_path)
    for case in contract["cases"]:
        case_started = time.perf_counter()
        provider, metadata = _worker_case_provider(
            worker,
            contract,
            case,
            rlb_arm_cache=rlb_arm_cache if worker == "rlb" else None,
        )
        _assert_fresh_worker_import_boundary(worker, after_import=True)
        scale = _frequency_scale(contract, str(case["geometry"]))
        print(
            f"[{worker}] {case['case_id']}: independent primary+verification scan started",
            flush=True,
        )
        inventory = seed_free_root_inventory(
            provider,
            scale,
            policy,
            case_id=str(case["case_id"]),
            builder_id=builder_id,
            contract_sha256=contract_sha256,
        )
        inventories.append(inventory)
        case_metadata.append(
            {
                "case_id": case["case_id"],
                "geometry": case["geometry"],
                "L1": case["L1"],
                "L2": case["L2"],
                "beta_deg": case["beta_deg"],
                "frequency_scale_Omega_per_omega": scale,
                "inventory_status": inventory.status,
                "inventory_sha256": inventory.inventory_sha256,
                "primary": _scan_manifest_summary(inventory.primary),
                "verification": _scan_manifest_summary(inventory.verification),
                "primary_verification_independent_agreement": (
                    inventory.independent_agreement
                ),
                "maximum_primary_verification_relative": (
                    inventory.maximum_primary_verification_relative
                ),
                "unresolved_low_sigma_count": inventory.unresolved_low_sigma_count,
                "guard_available": inventory.guard_available,
                "guard_not_at_scan_boundary": inventory.guard_not_at_scan_boundary,
                "elapsed_seconds": time.perf_counter() - case_started,
                "model_metadata": metadata,
            }
        )
        print(
            f"[{worker}] {case['case_id']}: {inventory.status}, "
            f"slots={len(inventory.slots)}, sha256={inventory.inventory_sha256}",
            flush=True,
        )

    candidate_rows = [
        _candidate_row(candidate)
        for inventory in inventories
        for scan in (inventory.primary, inventory.verification)
        for candidate in scan.candidates
    ]
    rejected_rows = [
        _candidate_row(candidate)
        for inventory in inventories
        for scan in (inventory.primary, inventory.verification)
        for candidate in scan.rejected_candidates
    ]
    root_rows = [
        _root_row(inventory, slot)
        for inventory in inventories
        for slot in inventory.slots
    ]
    candidate_path = output_directory / f"{prefix}_root_candidates.csv"
    rejected_path = output_directory / f"{prefix}_rejected_candidates.csv"
    roots_path = output_directory / f"{prefix}_roots.csv"
    _write_csv(candidate_path, candidate_rows)
    _write_csv(rejected_path, rejected_rows)
    _write_csv(roots_path, root_rows)

    semantic_payload = {
        "contract_sha256": contract_sha256,
        "builder_id": builder_id,
        "root_policy": _json_value(policy.__dict__),
        "cases": [
            {
                "case_id": inventory.case_id,
                "case_inventory_sha256": inventory.inventory_sha256,
                "status": inventory.status,
                "primary_slot_count": len(inventory.primary.slots),
                "verification_slot_count": len(inventory.verification.slots),
            }
            for inventory in inventories
        ],
    }
    overall_inventory_sha = _sha256_bytes(
        json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode(
            "utf-8"
        )
    )
    manifest_path = output_directory / f"{prefix}_inventory_manifest.json"
    worker_manifest = {
        "stage": STAGE_ID,
        "worker": worker,
        "builder_id": builder_id,
        "generated_utc": _utc_now(),
        "process_id": os.getpid(),
        "contract_sha256": contract_sha256,
        "root_policy": _json_value(policy.__dict__),
        "seed_free": True,
        "cross_model_roots_read": False,
        "fresh_process_import_boundary_passed": True,
        "rlb_arm_cache_entries": len(rlb_arm_cache) if worker == "rlb" else 0,
        "cases": case_metadata,
        "all_case_statuses_pass": all(item.status == "PASS" for item in inventories),
        "semantic_inventory_sha256_before_comparison": overall_inventory_sha,
        "output_sha256": {
            candidate_path.name: _sha256_file(candidate_path),
            rejected_path.name: _sha256_file(rejected_path),
            roots_path.name: _sha256_file(roots_path),
        },
        "elapsed_seconds": time.perf_counter() - started,
        "comparison_performed": False,
    }
    _write_json(manifest_path, worker_manifest)
    return worker_manifest


def _read_preliminary_gate(output_directory: Path) -> dict[str, Any]:
    path = output_directory / "preliminary_gate.json"
    if not path.is_file():
        raise RuntimeError("Preliminary gate file is missing.")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("preliminary_status") != "PASS":
        raise RuntimeError("Preliminary gate is not PASS.")
    return payload


def _verify_worker_output_hashes(
    output_directory: Path, manifest: Mapping[str, Any]
) -> dict[str, str]:
    recorded = {
        str(name): str(digest) for name, digest in manifest["output_sha256"].items()
    }
    actual: dict[str, str] = {}
    for name, expected in recorded.items():
        path = output_directory / name
        if not path.is_file():
            raise RuntimeError(f"Frozen inventory output is missing: {path}")
        actual[name] = _sha256_file(path)
        if actual[name] != expected:
            raise RuntimeError(
                f"Frozen inventory output hash mismatch for {name}: "
                f"expected {expected}, got {actual[name]}."
            )
    return actual


def run_isolated_inventory_subprocesses(
    contract_path: Path, output_directory: Path
) -> dict[str, Any]:
    """Run RLB then legacy workers without passing roots between processes."""

    started = time.perf_counter()
    # The gate is cheap relative to the inventories.  Re-run it unconditionally
    # so an old PASS file can never authorize searches after a model edit.
    run_preliminary_gates(contract_path, output_directory)
    preliminary = _read_preliminary_gate(output_directory)
    if preliminary.get("contract_sha256") != _sha256_file(contract_path):
        raise RuntimeError("Preliminary and requested contract hashes differ.")
    commands: list[dict[str, Any]] = []
    environment = dict(os.environ)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    for worker in ("rlb", "legacy"):
        command = [
            sys.executable,
            str(Path(__file__).resolve()),
            "--worker",
            worker,
            "--contract",
            str(contract_path),
            "--output",
            str(output_directory),
        ]
        worker_started = time.perf_counter()
        completed = subprocess.run(
            command,
            cwd=REPOSITORY_ROOT,
            env=environment,
            check=False,
            text=True,
        )
        commands.append(
            {
                "worker": worker,
                "command": command,
                "exit_code": completed.returncode,
                "elapsed_seconds": time.perf_counter() - worker_started,
            }
        )
        if completed.returncode != 0:
            raise RuntimeError(f"Independent {worker} worker failed with exit code {completed.returncode}.")

    rlb_manifest = json.loads(
        (output_directory / "rlb_inventory_manifest.json").read_text(encoding="utf-8")
    )
    legacy_manifest = json.loads(
        (output_directory / "legacy_rectangular_inventory_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    if rlb_manifest.get("contract_sha256") != _sha256_file(contract_path):
        raise RuntimeError("RLB inventory manifest contract hash mismatch.")
    if legacy_manifest.get("contract_sha256") != _sha256_file(contract_path):
        raise RuntimeError("Legacy inventory manifest contract hash mismatch.")
    rlb_output_hashes = _verify_worker_output_hashes(output_directory, rlb_manifest)
    legacy_output_hashes = _verify_worker_output_hashes(
        output_directory, legacy_manifest
    )
    run_manifest = {
        "stage": STAGE_ID,
        "generated_utc": _utc_now(),
        "workflow_state": "BOTH_INVENTORIES_FROZEN_COMPARISON_NOT_IMPLEMENTED",
        "preliminary_gate": preliminary,
        "commands": commands,
        "rlb_inventory_sha256_before_comparison": rlb_manifest[
            "semantic_inventory_sha256_before_comparison"
        ],
        "legacy_inventory_sha256_before_comparison": legacy_manifest[
            "semantic_inventory_sha256_before_comparison"
        ],
        "rlb_output_sha256_before_comparison": rlb_output_hashes,
        "legacy_output_sha256_before_comparison": legacy_output_hashes,
        "all_worker_case_statuses_pass": bool(
            rlb_manifest.get("all_case_statuses_pass")
            and legacy_manifest.get("all_case_statuses_pass")
        ),
        "cross_seeding": False,
        "comparison_performed": False,
        "mode_reconstruction_performed": False,
        "elapsed_seconds": time.perf_counter() - started,
    }
    _write_json(output_directory / "run_manifest.json", run_manifest)
    return run_manifest


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--manifest-only",
        action="store_true",
        help="write only model_manifest.json and case_contract.json",
    )
    mode.add_argument(
        "--preliminary-only",
        action="store_true",
        help="run algebraic, adapter, and local-arm gates without spectra",
    )
    mode.add_argument(
        "--run-inventories",
        action="store_true",
        help="after PASS preliminary gate, launch both isolated root workers",
    )
    parser.add_argument(
        "--worker",
        choices=("rlb", "legacy"),
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--preflight-worker",
        choices=("rlb", "legacy"),
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT_PATH)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_DIRECTORY)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = _build_parser().parse_args(argv)
    contract_path = arguments.contract.resolve()
    output_directory = arguments.output.resolve()
    if arguments.worker and arguments.preflight_worker:
        raise SystemExit("--worker and --preflight-worker are mutually exclusive.")
    internal_worker = arguments.worker or arguments.preflight_worker
    if internal_worker and (
        arguments.manifest_only
        or arguments.preliminary_only
        or arguments.run_inventories
    ):
        raise SystemExit("An internal worker cannot be combined with a parent workflow mode.")
    if arguments.preflight_worker:
        payload = run_preflight_snapshot_worker(
            arguments.preflight_worker, contract_path, output_directory
        )
    elif arguments.worker:
        payload = run_inventory_worker(
            arguments.worker, contract_path, output_directory
        )
    elif arguments.manifest_only:
        payload = write_manifest_only(contract_path, output_directory)
    elif arguments.run_inventories:
        payload = run_isolated_inventory_subprocesses(
            contract_path, output_directory
        )
    else:
        # Safe default: stop after preliminary gates.  A full inventory requires
        # the explicit --run-inventories opt-in.
        payload = run_preliminary_gates(contract_path, output_directory)
    printable = payload
    if arguments.worker:
        printable = {
            "stage": payload["stage"],
            "worker": payload["worker"],
            "all_case_statuses_pass": payload["all_case_statuses_pass"],
            "semantic_inventory_sha256_before_comparison": payload[
                "semantic_inventory_sha256_before_comparison"
            ],
            "output_sha256": payload["output_sha256"],
            "elapsed_seconds": payload["elapsed_seconds"],
        }
    print(json.dumps(_json_value(printable), ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
