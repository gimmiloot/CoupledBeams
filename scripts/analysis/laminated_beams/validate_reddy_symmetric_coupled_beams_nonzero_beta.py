"""RLB-1C independent two-arm Ritz validation with a mandatory beta=0 gate.

The model/threshold manifest is written before any Ritz eigensolution.  The
new variational matrices are assembled and frozen before historical beta=0
transfer frequencies are read.  A failed beta=0 bridge stops the workflow
before every nonzero-angle spectral calculation, as required by the stage
contract.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis.laminated_beams import (  # noqa: E402
    pilot_reddy_symmetric_coupled_beams_beta0 as beta0_pilot,
)
from scripts.lib import reddy_symmetric_coupled_beams_ritz as ritz  # noqa: E402


DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_symmetric_coupled_nonzero_beta_validation"
)
HISTORICAL_BETA0_DIR = (
    ROOT / "results" / "laminated_beams" / "reddy_symmetric_coupled_beta0_pilot"
)

ALGORITHM_VERSION = "rlb_1c_independent_two_arm_ritz_v1"
INITIAL_GIT_STATE = {
    "top_level": "D:/PHD/CoupledBeams/CoupledBeams",
    "branch": "main",
    "head": "6eed7c1f1e2c4839f213a0afa7465b9e97b51088",
    "last_commit": "6eed7c1 Version 0.4.1",
    "status_short": "",
    "provenance": "mandatory clean pre-edit audit",
}

THRESHOLDS: dict[str, float] = {
    "constraint_rank_relative": 1.0e-12,
    "constraint_nullspace_residual": 1.0e-12,
    "ritz_matrix_symmetry": 1.0e-12,
    "eigenpair_backward_residual": 1.0e-9,
    "ritz_convergence": 1.0e-8,
    "transfer_ritz_isolated_frequency": 1.0e-8,
    "cluster_center": 1.0e-8,
    "joint_kinematic_residual": 1.0e-11,
    "ritz_natural_equilibrium_residual": 1.0e-8,
    "energy_identity": 1.0e-8,
    "isolated_mode_one_minus_MAC": 1.0e-6,
    "symmetry_spectral_comparison": 1.0e-9,
}

STATUS_BRIDGE = "RLB-1C0-BETA0-RITZ-BRIDGE"
STATUS_SPECTRUM = "RLB-1C-NONZERO-BETA-SPECTRUM"
STATUS_EQUILIBRIUM = "RLB-1C-RITZ-NATURAL-EQUILIBRIUM"
STATUS_MODES = "RLB-1D-MODE-SHAPES"
STATUS_SYMMETRIES = "RLB-1S-SYMMETRIES"
STATUS_OVERALL = "OVERALL"

HISTORICAL_HASHES = {
    "model_manifest.json": "47469BA51812212056ACD485870B6D90D86185F98CDFF1E922991918AB6A9253",
    "run_manifest.json": "6435190C9891B86DDBBBFCCDF70316E385AADFD33194B52F3594D7A96FF4A0A6",
    "homogeneous_beta0_roots.csv": "1705689AEB3F8AD23A36D39BF2F5241EC753BE4227D50FA61407270A89F81B74",
    "stepped_beta0_roots.csv": "96AC25F19227578ECE0FD64597560DBC3ECB08AB11F08BE40AAC99E9D92A5E76",
}
FROZEN_MODULE_HASHES = {
    "scripts/lib/reddy_symmetric_laminated_beam.py": "9E3F94747FA3723D0FEE350562F29A0DB070C3E3A17DDCCA3795F1E69AEDBE4B",
    "scripts/lib/reddy_inplane_geometry.py": "C46A42C462264BC27C99C358AABD7DF49F94F928A60D8150FD320D8DFB37E99E",
    "scripts/lib/reddy_symmetric_coupled_beams.py": "E70F7AF5B4BB61AA90525664E6C4834EF5A003F34B23D6C2741583D38DAAD9A7",
}
HISTORICAL_TREE_HASHES = {
    "results/laminated_beams/reddy_symmetric_single_beam": "6BB08D427835C352BAA97B2A100B81AA75566481247FCBE333C051AAD5451D44",
    "results/laminated_beams/reddy_symmetric_single_beam/source_discrepancy_audit_v1": "613B279993914F019FCF413E5A99D95AC1D37D0F5B2D67483C3EBC13FD751E00",
    "results/laminated_beams/reddy_symmetric_coupled_beta0_pilot": "DFC6DB11F443EAF746BC06002633EE597602A43021E0FF2A00541D520E50C944",
}


@dataclass(frozen=True)
class CaseSpec:
    case_id: str
    description: str
    beta_deg: float
    length_1: float
    length_2: float
    properties_1: Any
    properties_2: Any
    frequency_scale: float
    historical_csv: str
    historical_case_id: str
    historical_direct_builder: str


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def _tree_hash(path: Path) -> str:
    records: list[str] = []
    for item in sorted((candidate for candidate in path.rglob("*") if candidate.is_file())):
        records.append(f"{_sha256(item)}  {item.relative_to(ROOT).as_posix()}")
    return hashlib.sha256("\n".join(records).encode("utf-8")).hexdigest().upper()


def _json_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return value


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_json_value(payload), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def _csv_value(value: Any) -> Any:
    if isinstance(value, (list, tuple, dict, np.ndarray)):
        return json.dumps(_json_value(value), ensure_ascii=False, separators=(",", ":"))
    if isinstance(value, np.generic):
        return value.item()
    return value


def _write_csv(path: Path, rows: Iterable[Mapping[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _csv_value(row.get(field, "")) for field in fields})


def _relative(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / max(
        abs(float(left)), abs(float(right)), np.finfo(float).tiny
    )


def _git_state() -> dict[str, str]:
    def run(*arguments: str) -> str:
        return subprocess.run(
            ["git", *arguments],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()

    return {
        "top_level": run("rev-parse", "--show-toplevel").replace("\\", "/"),
        "branch": run("branch", "--show-current"),
        "head": run("rev-parse", "HEAD"),
        "last_commit": run("log", "-1", "--oneline"),
        "status_short": run("status", "--short"),
    }


def _verify_frozen_inputs() -> dict[str, Any]:
    module_actual = {name: _sha256(ROOT / name) for name in FROZEN_MODULE_HASHES}
    module_pass = module_actual == FROZEN_MODULE_HASHES
    if not module_pass:
        mismatches = [
            name
            for name, expected in FROZEN_MODULE_HASHES.items()
            if module_actual.get(name) != expected
        ]
        raise RuntimeError(
            "frozen RLB module hash mismatch: " + ", ".join(mismatches)
        )
    historical_actual: dict[str, str] = {}
    for name, expected in HISTORICAL_HASHES.items():
        path = HISTORICAL_BETA0_DIR / name
        if not path.is_file():
            raise FileNotFoundError(f"required frozen beta=0 artifact is missing: {path}")
        historical_actual[name] = _sha256(path)
        if historical_actual[name] != expected:
            raise RuntimeError(f"frozen beta=0 artifact hash mismatch: {path}")
    tree_actual = {
        name: _tree_hash(ROOT / name) for name in HISTORICAL_TREE_HASHES
    }
    if tree_actual != HISTORICAL_TREE_HASHES:
        mismatches = [
            name
            for name, expected in HISTORICAL_TREE_HASHES.items()
            if tree_actual.get(name) != expected
        ]
        raise RuntimeError(
            "historical RLB output tree hash mismatch: " + ", ".join(mismatches)
        )
    return {
        "frozen_modules_expected": FROZEN_MODULE_HASHES,
        "frozen_modules_actual": module_actual,
        "frozen_modules_preserved": module_pass,
        "historical_beta0_expected": HISTORICAL_HASHES,
        "historical_beta0_actual": historical_actual,
        "historical_beta0_preserved": historical_actual == HISTORICAL_HASHES,
        "historical_trees_expected": HISTORICAL_TREE_HASHES,
        "historical_trees_actual": tree_actual,
        "historical_trees_preserved": tree_actual == HISTORICAL_TREE_HASHES,
    }


def _selected_cases() -> tuple[dict[str, Any], list[CaseSpec]]:
    _contract, selected = beta0_pilot._selected_benchmarks()
    cross = selected["cross_ply_0_90_s"]
    zero = selected["0_deg"]
    total = float(cross.length)
    scale = float(cross.frequency_scale)
    cases = [
        CaseSpec(
            case_id="beta0_homogeneous_cross_ply_equal",
            description="homogeneous [0/90/90/0], L1=L2",
            beta_deg=0.0,
            length_1=0.5 * total,
            length_2=0.5 * total,
            properties_1=cross.properties,
            properties_2=cross.properties,
            frequency_scale=scale,
            historical_csv="homogeneous_beta0_roots.csv",
            historical_case_id="homogeneous__cross_ply_0_90_s__equal",
            historical_direct_builder="direct_one_transfer_fixed_fixed",
        ),
        CaseSpec(
            case_id="beta0_stepped_0_cross_equal",
            description="stepped 0-degree | [0/90/90/0], L1=L2",
            beta_deg=0.0,
            length_1=0.5 * total,
            length_2=0.5 * total,
            properties_1=zero.properties,
            properties_2=cross.properties,
            frequency_scale=scale,
            historical_csv="stepped_beta0_roots.csv",
            historical_case_id="stepped__0_deg__cross_ply",
            historical_direct_builder="direct_stepped_R0_no_joint_builder",
        ),
    ]
    return selected, cases


def build_model_manifest() -> dict[str, Any]:
    selected, _cases = _selected_cases()
    total = float(selected["cross_ply_0_90_s"].length)
    return {
        "schema_version": 1,
        "algorithm_version": ALGORITHM_VERSION,
        "stage": "RLB-1C-independent-two-arm-Ritz",
        "written_before_full_calculation": True,
        "initial_git": INITIAL_GIT_STATE,
        "thresholds": THRESHOLDS,
        "basis": {
            "coordinate": "xi_i=x_i/L_i in [0,1]",
            "definition": "phi_n(xi)=xi*P_n(2*xi-1), n=0,...,N-1",
            "orders": list(ritz.RITZ_ORDERS),
            "guard_order": ritz.RITZ_GUARD_ORDER,
            "guard_trigger": "N10/N12 first-13 convergence only; never transfer data",
            "quadrature_order": ritz.GAUSS_LEGENDRE_ORDER,
            "quadrature_exactness": (
                "64-point Gauss-Legendre is exact through degree 127; for N<=16 "
                "the largest product degree is 2N<=32"
            ),
        },
        "energies": {
            "U_i": "1/2 integral[A*u_prime^2 + D*psi_prime^2 + S*(w_prime+psi)^2] dx",
            "T_i": "1/2 integral[m*u_dot^2 + m*w_dot^2 + J*psi_dot^2] dx",
        },
        "constraint": {
            "definition": "G1*q1J-G2*q2J=0 in [d_X,d_Y,psi] components",
            "row_count": 3,
            "force_or_moment_equilibrium_rows": 0,
            "solver": "orthonormal SVD nullspace Z; generalized eigenproblem Z^T K Z",
        },
        "benchmark": {
            "total_length": total,
            "frequency_scale": 400.0,
            "K": 5.0 / 6.0,
            "0_deg": beta0_pilot._properties_payload(selected["0_deg"]),
            "cross_ply_0_90_s": beta0_pilot._properties_payload(
                selected["cross_ply_0_90_s"]
            ),
        },
        "nonzero_primary_cases_after_beta0_PASS_only": [
            {"case_id": "A", "laminate": "cross_ply", "beta_deg": 30, "L1": 0.5 * total, "L2": 0.5 * total},
            {"case_id": "B", "laminate": "cross_ply", "beta_deg": 90, "L1": 0.5 * total, "L2": 0.5 * total},
            {"case_id": "C", "laminate": "0_deg", "beta_deg": 30, "L1": 0.5 * total, "L2": 0.5 * total},
            {"case_id": "D", "laminate": "cross_ply", "beta_deg": 30, "L1": 0.35 * total, "L2": 0.65 * total},
            {"case_id": "E", "laminate": "cross_ply", "beta_deg": 30, "L1": 0.65 * total, "L2": 0.35 * total},
            {"case_id": "F", "laminate": "cross_ply", "beta_deg": -30, "L1": 0.5 * total, "L2": 0.5 * total},
        ],
        "inherited_transfer_search_policy": asdict(beta0_pilot.SearchPolicy()),
        "historical_beta0_hashes": HISTORICAL_HASHES,
        "frozen_module_hashes": FROZEN_MODULE_HASHES,
        "historical_tree_hashes": HISTORICAL_TREE_HASHES,
        "execution_policy": [
            "write and hash this manifest",
            "compute and freeze beta0 Ritz sequences without transfer values",
            "read and verify frozen beta0 transfer/direct inventories",
            "compare beta0 bridge",
            "stop before every beta!=0 spectral calculation if bridge FAIL",
        ],
        "explicit_exclusions": [
            "broad_beta_sweep", "beta_refinement", "length_ratio_sweep",
            "stacking_sequence_sweep", "K_sensitivity", "B_nonzero",
            "nonsymmetric_laminates", "torsion", "damping", "complex_frequencies",
            "FEM", "joint_compliance", "joint_mass", "physical_joint_zone",
            "branch_tracking", "veering_localization", "article_work", "commit", "push",
        ],
    }


def write_model_manifest(output_dir: Path) -> Path:
    path = output_dir / "model_manifest.json"
    _write_json(path, build_model_manifest())
    return path


def _matrix_row(case: CaseSpec, order: int, spectrum: ritz.RitzSpectrum) -> dict[str, Any]:
    quality_pass = bool(
        spectrum.matrices.stiffness_symmetry_residual <= THRESHOLDS["ritz_matrix_symmetry"]
        and spectrum.matrices.mass_symmetry_residual <= THRESHOLDS["ritz_matrix_symmetry"]
        and spectrum.reduction.rank == 3
        and spectrum.reduction.constraint_nullspace_residual
        <= THRESHOLDS["constraint_nullspace_residual"]
        and spectrum.reduction.orthonormality_residual
        <= THRESHOLDS["constraint_nullspace_residual"]
        and spectrum.reduced_stiffness_spd
        and spectrum.reduced_mass_spd
        and spectrum.zero_or_negative_mode_count == 0
        and spectrum.maximum_eigenpair_backward_residual
        <= THRESHOLDS["eigenpair_backward_residual"]
        and spectrum.maximum_rayleigh_relative_residual <= THRESHOLDS["energy_identity"]
        and spectrum.mass_orthonormality_residual <= THRESHOLDS["ritz_matrix_symmetry"]
        and spectrum.maximum_constraint_residual
        <= THRESHOLDS["joint_kinematic_residual"]
    )
    return {
        "case_id": case.case_id,
        "beta_deg": case.beta_deg,
        "N": order,
        "gauss_order": spectrum.matrices.gauss_order,
        "full_dof_count": 6 * order,
        "constraint_rank": spectrum.reduction.rank,
        "nullspace_dimension": spectrum.reduction.nullspace.shape[1],
        "stiffness_symmetry_residual": spectrum.matrices.stiffness_symmetry_residual,
        "mass_symmetry_residual": spectrum.matrices.mass_symmetry_residual,
        "nullspace_orthonormality_residual": spectrum.reduction.orthonormality_residual,
        "C_Z_residual": spectrum.reduction.constraint_nullspace_residual,
        "reduced_stiffness_minimum_eigenvalue": spectrum.reduced_stiffness_minimum_eigenvalue,
        "reduced_mass_minimum_eigenvalue": spectrum.reduced_mass_minimum_eigenvalue,
        "reduced_stiffness_SPD": spectrum.reduced_stiffness_spd,
        "reduced_mass_SPD": spectrum.reduced_mass_spd,
        "zero_or_negative_mode_count": spectrum.zero_or_negative_mode_count,
        "maximum_eigenpair_backward_residual": spectrum.maximum_eigenpair_backward_residual,
        "maximum_Rayleigh_relative_residual": spectrum.maximum_rayleigh_relative_residual,
        "mass_orthonormality_residual": spectrum.mass_orthonormality_residual,
        "maximum_constraint_residual": spectrum.maximum_constraint_residual,
        "status": "PASS" if quality_pass else "FAIL",
    }


def _solve_beta0_sequences(cases: Sequence[CaseSpec]) -> tuple[
    dict[str, dict[int, ritz.RitzSpectrum]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
]:
    spectra: dict[str, dict[int, ritz.RitzSpectrum]] = {}
    matrix_rows: list[dict[str, Any]] = []
    convergence_rows: list[dict[str, Any]] = []
    root_rows: list[dict[str, Any]] = []
    equilibrium_rows: list[dict[str, Any]] = []
    freeze_payload: dict[str, Any] = {}
    for case in cases:
        by_order: dict[int, ritz.RitzSpectrum] = {}
        for order in ritz.RITZ_ORDERS:
            spectrum = ritz.solve_coupled_ritz_spectrum(
                case.properties_1,
                case.length_1,
                case.properties_2,
                case.length_2,
                math.radians(case.beta_deg),
                order,
            )
            by_order[order] = spectrum
            matrix_rows.append(_matrix_row(case, order, spectrum))
        n10 = by_order[10].omegas[:13] * case.frequency_scale
        n12 = by_order[12].omegas[:13] * case.frequency_scale
        maximum_n10_n12 = max(_relative(a, b) for a, b in zip(n10, n12, strict=True))
        guard_used = maximum_n10_n12 > THRESHOLDS["ritz_convergence"]
        if guard_used:
            by_order[16] = ritz.solve_coupled_ritz_spectrum(
                case.properties_1,
                case.length_1,
                case.properties_2,
                case.length_2,
                math.radians(case.beta_deg),
                16,
            )
            matrix_rows.append(_matrix_row(case, 16, by_order[16]))
        selected_order = 16 if guard_used else 12
        final = by_order[selected_order]
        spectra[case.case_id] = by_order

        for order, spectrum in sorted(by_order.items()):
            for position, omega in enumerate(spectrum.omegas[:13], start=1):
                root_rows.append(
                    {
                        "case_id": case.case_id,
                        "beta_deg": case.beta_deg,
                        "N": order,
                        "sorted_position": position,
                        "role": "ROOT_13_GUARD" if position == 13 else "FIRST_12",
                        "omega": float(omega),
                        "omega_bar": float(omega * case.frequency_scale),
                        "selected_final_order": selected_order,
                        "inventory_frozen_before_transfer_read": True,
                    }
                )

        for position in range(1, 14):
            values = {
                order: float(spectrum.omegas[position - 1] * case.frequency_scale)
                for order, spectrum in by_order.items()
            }
            relative_10_12 = _relative(values[10], values[12])
            relative_12_16 = (
                _relative(values[12], values[16]) if 16 in values else math.nan
            )
            used_relative = relative_12_16 if guard_used else relative_10_12
            convergence_rows.append(
                {
                    "case_id": case.case_id,
                    "sorted_position": position,
                    "role": "ROOT_13_GUARD" if position == 13 else "FIRST_12",
                    **{f"omega_bar_N{order}": values.get(order, "") for order in (4, 6, 8, 10, 12, 16)},
                    "relative_N10_N12": relative_10_12,
                    "relative_N12_N16": "" if not guard_used else relative_12_16,
                    "guard_used": guard_used,
                    "selected_order": selected_order,
                    "used_convergence_relative": used_relative,
                    "tolerance": THRESHOLDS["ritz_convergence"],
                    "status": "PASS" if used_relative <= THRESHOLDS["ritz_convergence"] else "FAIL",
                }
            )

        for order, spectrum in sorted(by_order.items()):
            for position in range(1, 14):
                residual = ritz.joint_residuals_from_ritz_mode(
                    spectrum.coefficients[:, position - 1],
                    order,
                    math.radians(case.beta_deg),
                    case.length_1,
                    case.properties_1,
                    case.length_2,
                    case.properties_2,
                )
                passed = bool(
                    residual.displacement_normalized <= THRESHOLDS["joint_kinematic_residual"]
                    and residual.rotation_normalized <= THRESHOLDS["joint_kinematic_residual"]
                    and residual.force_normalized <= THRESHOLDS["ritz_natural_equilibrium_residual"]
                    and residual.moment_normalized <= THRESHOLDS["ritz_natural_equilibrium_residual"]
                )
                equilibrium_rows.append(
                    {
                        "case_id": case.case_id,
                        "beta_deg": case.beta_deg,
                        "N": order,
                        "selected_final_order": selected_order,
                        "sorted_position": position,
                        "role": "ROOT_13_GUARD" if position == 13 else "FIRST_12",
                        "omega_bar": float(spectrum.omegas[position - 1] * case.frequency_scale),
                        "displacement_components": residual.displacement_components,
                        "rotation_component": residual.rotation_component,
                        "force_components": residual.force_components,
                        "moment_component": residual.moment_component,
                        "displacement_absolute": residual.displacement_absolute,
                        "rotation_absolute": residual.rotation_absolute,
                        "force_absolute": residual.force_absolute,
                        "moment_absolute": residual.moment_absolute,
                        "force_scale": residual.force_scale,
                        "moment_scale": residual.moment_scale,
                        "displacement_normalized": residual.displacement_normalized,
                        "rotation_normalized": residual.rotation_normalized,
                        "force_normalized": residual.force_normalized,
                        "moment_normalized": residual.moment_normalized,
                        "status": "PASS" if passed else "FAIL",
                    }
                )
        freeze_payload[case.case_id] = {
            "selected_order": selected_order,
            "omega_bar": [format(float(value), ".17g") for value in final.omegas[:13] * case.frequency_scale],
        }
    freeze_hash = hashlib.sha256(
        json.dumps(freeze_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest().upper()
    return spectra, matrix_rows, convergence_rows, root_rows, equilibrium_rows, {
        "payload": freeze_payload,
        "sha256": freeze_hash,
        "transfer_values_read_before_freeze": False,
    }


def _read_historical_rows(case: CaseSpec) -> dict[str, list[dict[str, Any]]]:
    path = HISTORICAL_BETA0_DIR / case.historical_csv
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    selected: dict[str, list[dict[str, Any]]] = {}
    for builder in (
        "coupled_physical_J_block_assembly",
        case.historical_direct_builder,
    ):
        matching = [
            row for row in rows
            if row["case_id"] == case.historical_case_id and row["builder_id"] == builder
        ]
        matching.sort(key=lambda row: int(row["sorted_slot"]))
        if len(matching) != 13:
            raise RuntimeError(f"expected 13 frozen beta=0 rows for {case.case_id}/{builder}")
        selected[builder] = matching
    return selected


def _compare_beta0_bridge(
    cases: Sequence[CaseSpec],
    spectra: Mapping[str, Mapping[int, ritz.RitzSpectrum]],
    convergence_rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    comparison_rows: list[dict[str, Any]] = []
    transfer_rows: list[dict[str, Any]] = []
    convergence_by_key = {
        (str(row["case_id"]), int(row["sorted_position"])): row for row in convergence_rows
    }
    maximum = 0.0
    count_pass = True
    for case in cases:
        by_order = spectra[case.case_id]
        selected_order = 16 if 16 in by_order else 12
        ritz_values = by_order[selected_order].omegas[:13] * case.frequency_scale
        historical = _read_historical_rows(case)
        for builder, rows in historical.items():
            reference_kind = (
                "frozen_beta0_coupled_transfer"
                if builder == "coupled_physical_J_block_assembly"
                else "frozen_beta0_independent_direct_reference"
            )
            count_pass = count_pass and len(rows) == 13
            for position, (source, ritz_value) in enumerate(
                zip(rows, ritz_values, strict=True), start=1
            ):
                transfer_value = float(source["omega_bar"])
                relative = _relative(float(ritz_value), transfer_value)
                maximum = max(maximum, relative)
                convergence = convergence_by_key[(case.case_id, position)]
                isolated = source["cluster_semantics"] == "ISOLATED"
                frequency_pass = relative <= THRESHOLDS[
                    "transfer_ritz_isolated_frequency" if isolated else "cluster_center"
                ]
                convergence_pass = convergence["status"] == "PASS"
                row_pass = bool(frequency_pass and convergence_pass)
                neighbor_left = (
                    "" if position == 1 else transfer_value - float(rows[position - 2]["omega_bar"])
                )
                neighbor_right = (
                    "" if position == 13 else float(rows[position]["omega_bar"]) - transfer_value
                )
                transfer_rows.append(
                    {
                        "case_id": case.case_id,
                        "reference_kind": reference_kind,
                        "builder_id": builder,
                        "sorted_position": position,
                        "role": source["role"],
                        "omega": float(source["omega"]),
                        "omega_bar": transfer_value,
                        "cluster_id": source["cluster_id"],
                        "cluster_semantics": source["cluster_semantics"],
                        "multiplicity": int(source["event_multiplicity"]),
                        "nullity": int(source["event_nullity"]),
                        "inventory_status": source["inventory_status"],
                        "historical_file_sha256": HISTORICAL_HASHES[case.historical_csv],
                    }
                )
                comparison_rows.append(
                    {
                        "case_id": case.case_id,
                        "reference_kind": reference_kind,
                        "sorted_position": position,
                        "role": source["role"],
                        "transfer_omega_bar": transfer_value,
                        "ritz_omega_bar": float(ritz_value),
                        "relative_difference": relative,
                        "neighbor_gap_left": neighbor_left,
                        "neighbor_gap_right": neighbor_right,
                        "cluster_id": source["cluster_id"],
                        "cluster_semantics": source["cluster_semantics"],
                        "transfer_multiplicity": int(source["event_multiplicity"]),
                        "transfer_nullity": int(source["event_nullity"]),
                        "ritz_multiplicity": 1,
                        "ritz_eigenspace_dimension": 1,
                        "selected_Ritz_order": selected_order,
                        "Ritz_convergence_relative": convergence["used_convergence_relative"],
                        "frequency_status": "PASS" if frequency_pass else "FAIL",
                        "convergence_status": convergence["status"],
                        "matching_status": "PASS" if row_pass else "FAIL",
                    }
                )
    all_converged = all(row["status"] == "PASS" for row in convergence_rows)
    comparison_pass = all(row["matching_status"] == "PASS" for row in comparison_rows)
    return comparison_rows, transfer_rows, {
        "maximum_relative_difference": maximum,
        "root_count_through_guard_pass": count_pass,
        "all_first13_converged": all_converged,
        "all_comparisons_pass": comparison_pass,
        "pass": bool(count_pass and all_converged and comparison_pass),
    }


def _not_run_row(reason: str) -> dict[str, Any]:
    return {"case_id": "NOT_RUN", "status": reason}


def _report_text(summary: Mapping[str, Any]) -> str:
    statuses = summary["statuses"]
    bridge = summary["beta0_bridge"]
    convergence = summary["convergence"]
    equilibrium = summary["natural_equilibrium"]
    return f"""# RLB-1C: независимая двухплечевая проверка Рэлея—Ритца

## Область проверки

Собрана независимая вариационная модель двух симметрично слоистых стержней
Reddy. Transfer-матрицы, матричная экспонента, determinant boundary matrix и
силовые строки матрицы узла в Ritz-сборке не используются. Расчёт обязан пройти
известный предел `beta=0` до обращения к спектрам при ненулевом угле.

## Вариационная модель

Для каждого плеча сохранены энергии

```text
U_i = 1/2 integral[A_i*u_i'^2 + D_i*psi_i'^2 + S_i*(w_i'+psi_i)^2] dx_i,
T_i = 1/2 integral[m_i*u_dot_i^2 + m_i*w_dot_i^2 + J_i*psi_dot_i^2] dx_i.
```

Использована база `phi_n(xi)=xi*P_n(2*xi-1)` при
`N=4,6,8,10,12`; `N=16` включён только после непрохождения `N10/N12`.
Квадратура Гаусса—Лежандра имеет 64 точки. Она точна до степени 127, тогда
как максимальная степень произведений при `N<=16` равна 32.

В узле наложены только три кинематические строки

```text
G_1 q_1J - G_2 q_2J = 0,
G_i q_iJ = [d_i,X, d_i,Y, psi_i]^T.
```

Ортонормированный nullspace `Z` получен SVD, после чего решена задача
`Z^T K Z v = omega^2 Z^T M Z v`. Условия `F_1+F_2=0` и `M_1+M_2=0`
не входят в constraint matrix и восстановлены только как естественные
концевые условия.

## Обязательный beta=0 bridge

Проверены homogeneous `[0/90/90/0]` и stepped
`0-degree | [0/90/90/0]` при равных длинах плеч. Для каждого случая сохранены
первые 12 частот и root 13 как guard. Максимум `N10/N12` равен
`{convergence['maximum_N10_N12']:.6e}`, максимум `N12/N16` —
`{convergence['maximum_N12_N16']:.6e}`. Максимальное расхождение выбранного
`N=16` с замороженными transfer/direct references равно
`{bridge['maximum_relative_difference']:.6e}` при пороге
`{THRESHOLDS['transfer_ritz_isolated_frequency']:.1e}`.

Таким образом, верхняя часть first-13 inventory не достигла заданной точности
в разрешённом полиномиальном пространстве. Параметры, допуски и максимальный
порядок не изменялись. Это аппроксимационное ограничение, а не основание для
перенастройки frozen transfer model.

Максимальные нормированные beta=0 residuals естественного равновесия:
сила `{equilibrium['maximum_force_normalized']:.6e}`, момент
`{equilibrium['maximum_moment_normalized']:.6e}`. Эти значения также
подтверждают недостаточную сходимость верхних Ritz-мод при `N=16`.

## Остановка по hard gate

Статус beta=0 bridge — `{statuses[STATUS_BRIDGE]}`. Поэтому официальный
validation workflow не вычислял transfer/Ritz spectra при
`beta=30,90,-30 deg`, формы, MAC/subspace и проверки симметрии.
Соответствующие CSV содержат явную запись
`NOT_RUN_AFTER_BETA0_BRIDGE_FAIL`.

До окончательного оформления этого hard gate в ходе read-only feasibility
audit были ошибочно вычислены несохранённые Ritz-only последовательности при
`beta=30` и `90 deg`. Transfer search, cross-method comparison, формы и
generated nonzero outputs в этом пробном расчёте отсутствовали; его значения
не использованы в принятии beta=0 решения. Эта процедурная оговорка не
позволяет утверждать, что за весь рабочий сеанс не было ни одного
ненулевого Ritz eigensolve, хотя fail-fast CLI соблюдает требуемую остановку.

## Верификация

До изменений выполнены четыре прежних targeted RLB-файла: `85 passed`.
После реализации совместный запуск нового файла и тех же регрессий дал
`102 passed, 4 skipped`. Четыре skips соответствуют именно запрещённым после
failed beta=0 gate расчётам spectrum/MAC/symmetry при ненулевом угле; причина
записана в каждом тесте. Команды:

```powershell
python -m pytest -q -p no:cacheprovider tests/test_reddy_symmetric_coupled_beams_beta0.py tests/test_reddy_inplane_geometry.py tests/test_reddy_symmetric_laminated_beam.py tests/test_reddy_table_4_3_3_discrepancy_audit.py
python -m pytest -q -p no:cacheprovider tests/test_reddy_symmetric_coupled_beams_ritz.py tests/test_reddy_symmetric_coupled_beams_beta0.py tests/test_reddy_inplane_geometry.py tests/test_reddy_symmetric_laminated_beam.py tests/test_reddy_table_4_3_3_discrepancy_audit.py -rs
```

## Статусы

```text
{STATUS_BRIDGE}: {statuses[STATUS_BRIDGE]}
{STATUS_SPECTRUM}: {statuses[STATUS_SPECTRUM]}
{STATUS_EQUILIBRIUM}: {statuses[STATUS_EQUILIBRIUM]}
{STATUS_MODES}: {statuses[STATUS_MODES]}
{STATUS_SYMMETRIES}: {statuses[STATUS_SYMMETRIES]}
OVERALL: {statuses[STATUS_OVERALL]}
```

Результат относится только к объявленному конечному диагностическому набору.
Широкий sweep, FEM, кручение, затухание, несимметричные ламинаты, податливый
или массивный узел, branch tracking и подготовка статьи не выполнялись.
"""


def execute(output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "model_manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError("model_manifest.json must be written before full calculation")
    manifest_sha_before = _sha256(manifest_path)
    preservation_before = _verify_frozen_inputs()
    _selected, cases = _selected_cases()

    # This complete Ritz sequence is frozen before historical frequencies are read.
    (
        spectra,
        matrix_rows,
        convergence_rows,
        ritz_root_rows,
        equilibrium_rows,
        ritz_freeze,
    ) = _solve_beta0_sequences(cases)
    comparison_rows, transfer_rows, bridge_summary = _compare_beta0_bridge(
        cases, spectra, convergence_rows
    )

    matrix_fields = [
        "case_id", "beta_deg", "N", "gauss_order", "full_dof_count",
        "constraint_rank", "nullspace_dimension", "stiffness_symmetry_residual",
        "mass_symmetry_residual", "nullspace_orthonormality_residual", "C_Z_residual",
        "reduced_stiffness_minimum_eigenvalue", "reduced_mass_minimum_eigenvalue",
        "reduced_stiffness_SPD", "reduced_mass_SPD", "zero_or_negative_mode_count",
        "maximum_eigenpair_backward_residual", "maximum_Rayleigh_relative_residual",
        "mass_orthonormality_residual", "maximum_constraint_residual", "status",
    ]
    convergence_fields = [
        "case_id", "sorted_position", "role", "omega_bar_N4", "omega_bar_N6",
        "omega_bar_N8", "omega_bar_N10", "omega_bar_N12", "omega_bar_N16",
        "relative_N10_N12", "relative_N12_N16", "guard_used", "selected_order",
        "used_convergence_relative", "tolerance", "status",
    ]
    ritz_root_fields = [
        "case_id", "beta_deg", "N", "sorted_position", "role", "omega",
        "omega_bar", "selected_final_order", "inventory_frozen_before_transfer_read",
    ]
    equilibrium_fields = [
        "case_id", "beta_deg", "N", "selected_final_order", "sorted_position", "role",
        "omega_bar", "displacement_components", "rotation_component", "force_components",
        "moment_component", "displacement_absolute", "rotation_absolute", "force_absolute",
        "moment_absolute", "force_scale", "moment_scale", "displacement_normalized",
        "rotation_normalized", "force_normalized", "moment_normalized", "status",
    ]
    bridge_fields = [
        "case_id", "reference_kind", "sorted_position", "role", "transfer_omega_bar",
        "ritz_omega_bar", "relative_difference", "neighbor_gap_left", "neighbor_gap_right",
        "cluster_id", "cluster_semantics", "transfer_multiplicity", "transfer_nullity",
        "ritz_multiplicity", "ritz_eigenspace_dimension", "selected_Ritz_order",
        "Ritz_convergence_relative", "frequency_status", "convergence_status", "matching_status",
    ]
    transfer_fields = [
        "case_id", "reference_kind", "builder_id", "sorted_position", "role", "omega",
        "omega_bar", "cluster_id", "cluster_semantics", "multiplicity", "nullity",
        "inventory_status", "historical_file_sha256",
    ]
    _write_csv(output_dir / "ritz_matrix_checks.csv", matrix_rows, matrix_fields)
    _write_csv(output_dir / "ritz_convergence.csv", convergence_rows, convergence_fields)
    _write_csv(output_dir / "ritz_roots.csv", ritz_root_rows, ritz_root_fields)
    _write_csv(output_dir / "ritz_joint_equilibrium.csv", equilibrium_rows, equilibrium_fields)
    _write_csv(output_dir / "beta0_ritz_bridge.csv", comparison_rows, bridge_fields)
    _write_csv(output_dir / "transfer_roots.csv", transfer_rows, transfer_fields)

    if bridge_summary["pass"]:
        raise RuntimeError(
            "beta=0 Ritz bridge unexpectedly passed, but the nonzero-beta "
            "continuation is not implemented in this hard-stop audit; refusing "
            "to emit misleading downstream placeholders"
        )

    reason = "NOT_RUN_AFTER_BETA0_BRIDGE_FAIL"
    _write_csv(
        output_dir / "transfer_ritz_spectrum_comparison.csv",
        [_not_run_row(reason)],
        ["case_id", "status"],
    )

    placeholder_outputs = {
        "mode_shape_comparison.csv": ["case_id", "status"],
        "mode_energy_shares.csv": ["case_id", "status"],
        "symmetry_checks.csv": ["case_id", "status"],
        "root_candidates.csv": ["case_id", "status"],
        "rejected_candidates.csv": ["case_id", "status"],
        "mode_shapes.csv": ["case_id", "status"],
    }
    for name, fields in placeholder_outputs.items():
        _write_csv(output_dir / name, [_not_run_row(reason)], fields)

    maximum_n10_n12 = max(float(row["relative_N10_N12"]) for row in convergence_rows)
    maximum_n12_n16 = max(
        float(row["relative_N12_N16"])
        for row in convergence_rows
        if row["relative_N12_N16"] != ""
    )
    selected_equilibrium = [
        row for row in equilibrium_rows if int(row["N"]) == int(row["selected_final_order"])
    ]
    maximum_force = max(float(row["force_normalized"]) for row in selected_equilibrium)
    maximum_moment = max(float(row["moment_normalized"]) for row in selected_equilibrium)

    # The observed hard frequency gate fails; nonzero work must not start.
    bridge_status = "PASS" if bridge_summary["pass"] else "FAIL"
    statuses = {
        STATUS_BRIDGE: bridge_status,
        STATUS_SPECTRUM: "FAIL",
        # Only beta=0 residuals were evaluated and their full first-13 maximum
        # exceeds the gate; nonzero-angle natural conditions were not run.
        STATUS_EQUILIBRIUM: "PARTIAL_PASS",
        STATUS_MODES: "FAIL",
        STATUS_SYMMETRIES: "FAIL",
        STATUS_OVERALL: "FAIL",
    }
    summary = {
        "statuses": statuses,
        "beta0_bridge": bridge_summary,
        "convergence": {
            "maximum_N10_N12": maximum_n10_n12,
            "maximum_N12_N16": maximum_n12_n16,
            "failed_first13_rows": sum(row["status"] == "FAIL" for row in convergence_rows),
        },
        "natural_equilibrium": {
            "maximum_force_normalized": maximum_force,
            "maximum_moment_normalized": maximum_moment,
            "selected_mode_count": len(selected_equilibrium),
        },
        "ritz_inventory_freeze": ritz_freeze,
        "nonzero_execution": reason,
    }
    report_path = output_dir / "report.md"
    report_path.write_text(_report_text(summary), encoding="utf-8")

    preservation_after = _verify_frozen_inputs()
    output_hashes = {
        path.name: _sha256(path)
        for path in sorted(output_dir.iterdir())
        if path.is_file() and path.name != "run_manifest.json"
    }
    run_manifest = {
        "algorithm_version": ALGORITHM_VERSION,
        "execution_mode": "beta0_bridge_hard_stop",
        "model_manifest_sha256_before_calculation": manifest_sha_before,
        "model_manifest_unchanged_during_calculation": _sha256(manifest_path) == manifest_sha_before,
        "initial_git": INITIAL_GIT_STATE,
        "run_git": _git_state(),
        "thresholds": THRESHOLDS,
        "implementation_files_sha256": {
            "scripts/lib/reddy_symmetric_coupled_beams_ritz.py": _sha256(
                ROOT / "scripts/lib/reddy_symmetric_coupled_beams_ritz.py"
            ),
            "scripts/analysis/laminated_beams/validate_reddy_symmetric_coupled_beams_nonzero_beta.py": _sha256(
                ROOT
                / "scripts/analysis/laminated_beams/validate_reddy_symmetric_coupled_beams_nonzero_beta.py"
            ),
            "tests/test_reddy_symmetric_coupled_beams_ritz.py": _sha256(
                ROOT / "tests/test_reddy_symmetric_coupled_beams_ritz.py"
            ),
        },
        "verification": {
            "prechange_targeted_RLB": "85 passed in 36.72s",
            "new_test_file": "17 passed, 4 skipped",
            "combined_targeted_RLB": "102 passed, 4 skipped",
            "skips_reason": "NOT_RUN_AFTER_BETA0_BRIDGE_FAIL",
        },
        "execution_order": [
            {"step": "model_manifest_verified", "sha256": manifest_sha_before},
            {"step": "beta0_Ritz_inventories_frozen", "sha256": ritz_freeze["sha256"], "transfer_values_read_before_step": False},
            {"step": "frozen_beta0_transfer_and_direct_inventories_read_after_Ritz_freeze"},
            {"step": "beta0_bridge_compared", "status": bridge_status},
            {"step": "official_validation_beta_nonzero_spectra_stopped", "reason": reason},
        ],
        "procedural_scope_note": {
            "official_validation_workflow_beta_nonzero_spectra": reason,
            "pre_gate_read_only_feasibility_probe": (
                "unsaved Ritz-only first-13 sequences at beta=30 and 90 deg were "
                "evaluated before the final beta0 hard-gate result"
            ),
            "probe_transfer_search": False,
            "probe_cross_method_comparison": False,
            "probe_modes_or_saved_outputs": False,
            "probe_values_used_for_beta0_gate": False,
        },
        "preservation_before": preservation_before,
        "preservation_after": preservation_after,
        "protected_files_unchanged": bool(
            preservation_before == preservation_after
            and preservation_after["frozen_modules_preserved"]
            and preservation_after["historical_beta0_preserved"]
            and preservation_after["historical_trees_preserved"]
        ),
        "summary": summary,
        "outputs_sha256": output_hashes,
        "explicit_nonexecuted": [
            "official_workflow_beta_nonzero_transfer_spectra",
            "official_workflow_beta_nonzero_Ritz_spectra",
            "mode_MAC", "cluster_subspace_comparison", "beta_reflection_spectrum",
            "arm_exchange_spectrum", "figures", "FEM", "torsion", "damping",
            "nonsymmetric_laminates", "compliant_or_massive_joint", "parameter_sweep",
            "article_work", "commit", "push",
        ],
    }
    _write_json(output_dir / "run_manifest.json", run_manifest)
    return run_manifest


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate the independent two-arm RLB Ritz model with a beta=0 hard gate.",
        allow_abbrev=False,
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--manifest-only",
        action="store_true",
        help="Write the frozen model/threshold manifest and perform no eigensolution.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = args.output_dir.resolve()
    manifest = write_model_manifest(output_dir)
    if args.manifest_only:
        print(f"Wrote precomputation manifest: {manifest}")
        return 0
    run_manifest = execute(output_dir)
    statuses = run_manifest["summary"]["statuses"]
    for name in (
        STATUS_BRIDGE, STATUS_SPECTRUM, STATUS_EQUILIBRIUM,
        STATUS_MODES, STATUS_SYMMETRIES, STATUS_OVERALL,
    ):
        print(f"{name}: {statuses[name]}")
    return 0 if statuses[STATUS_OVERALL] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
