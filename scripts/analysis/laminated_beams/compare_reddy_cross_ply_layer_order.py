"""RLB-2A controlled pilot for four-ply cross-ply layer order.

This diagnostic workflow is not a new solver and is not a parameter sweep.  It
reuses the frozen RLB-0 laminate reduction, the RLB-1 coupled boundary matrix,
and the RLB-1 seed-free root inventory.  Its distinct input/output contract is
the full-block analytic constitutive check for exactly ``[0/90/90/0]`` and
``[90/0/0/90]``, followed by the declared ``beta=0`` and ``beta=30 deg``
finite spectral comparisons.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis.laminated_beams import (  # noqa: E402
    pilot_reddy_symmetric_coupled_beams_beta0 as beta0_pilot,
)
from scripts.lib import reddy_symmetric_coupled_beams as coupled  # noqa: E402
from scripts.lib import reddy_symmetric_laminated_beam as rlb  # noqa: E402


FloatArray = NDArray[np.float64]

DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_cross_ply_layer_order_pilot"
)

ALGORITHM_VERSION = "rlb_2a_layer_order_pilot_v1"
TASK_INITIAL_GIT_STATE = {
    "top_level": "D:/PHD/CoupledBeams/CoupledBeams",
    "branch": "main",
    "head": "5d63a4c519e1384230722143324045de3d3fe44c",
    "last_commit": "5d63a4c Version 0.4.3",
    "status_short": "",
    "provenance": "mandatory clean pre-edit audit",
}

STACK_OUTER_0 = (0.0, 90.0, 90.0, 0.0)
STACK_OUTER_90 = (90.0, 0.0, 0.0, 90.0)
STACK_IDS = ("OUTER_0", "OUTER_90")

STATUS_CONSTITUTIVE = "RLB-2A-CONSTITUTIVE-LAYER-ORDER"
STATUS_BETA0 = "RLB-2A-BETA0-FAMILY-CHECK"
STATUS_BETA30 = "RLB-2A-BETA30-COUPLED-PILOT"
STATUS_INVENTORY = "RLB-2A-ROOT-INVENTORY"
STATUS_OVERALL = "OVERALL"

THRESHOLDS: dict[str, float] = {
    "A_equality_relative": 1.0e-12,
    "B_scaled_residual": 1.0e-12,
    "shear_equality_relative": 1.0e-12,
    "mass_equality_relative": 1.0e-12,
    "analytic_D_formula_relative": 1.0e-12,
    "reduced_property_equality_relative": 1.0e-11,
    "beta0_axial_relative": 1.0e-10,
    "beta0_bending_order_allowance_relative": 1.0e-10,
    "beta0_coupled_union_relative": 1.0e-9,
    "beta30_ordered_frequency_allowance_relative": 1.0e-10,
    "root_singular_ratio": 1.0e-9,
    "boundary_residual": 1.0e-9,
}

FROZEN_PATHS = (
    ROOT / "scripts" / "lib" / "reddy_symmetric_laminated_beam.py",
    ROOT / "scripts" / "lib" / "reddy_inplane_geometry.py",
    ROOT / "scripts" / "lib" / "reddy_symmetric_coupled_beams.py",
    ROOT
    / "scripts"
    / "analysis"
    / "laminated_beams"
    / "pilot_reddy_symmetric_coupled_beams_beta0.py",
)

HISTORICAL_RESULT_DIRS = (
    ROOT / "results" / "laminated_beams" / "reddy_symmetric_single_beam",
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_symmetric_coupled_beta0_pilot",
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_symmetric_coupled_nonzero_beta_validation",
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_four_ply_isotropic_limit_validation",
)


@dataclass(frozen=True)
class LayerOrderCase:
    """One declared stack with the frozen RLB baseline realization."""

    stack_id: str
    stacking_sequence: tuple[float, ...]
    laminate: rlb.LaminateSection
    properties: rlb.BeamProperties
    total_length: float
    frequency_scale: float
    K_provenance: str

    @property
    def length_1(self) -> float:
        return 0.5 * self.total_length

    @property
    def length_2(self) -> float:
        return 0.5 * self.total_length


@dataclass(frozen=True)
class DirectFamilyInventory:
    """Existing RLB-0 axial/bending roots plus their RLB-1 union."""

    axial: tuple[Any, ...]
    bending: tuple[Any, ...]
    union_rows: tuple[dict[str, Any], ...]
    reconciliation: Mapping[str, int]


def _json_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value).replace("\\", "/")
    if isinstance(value, tuple):
        return [_json_value(item) for item in value]
    if isinstance(value, list):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, float) and not math.isfinite(value):
        return "inf" if value > 0.0 else ("-inf" if value < 0.0 else "nan")
    return value


def _csv_value(value: Any) -> Any:
    if isinstance(value, (tuple, list, dict, np.ndarray)):
        return json.dumps(
            _json_value(value), ensure_ascii=False, separators=(",", ":")
        )
    if isinstance(value, np.generic):
        return value.item()
    return value


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_json_value(dict(payload)), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _write_csv(
    path: Path,
    rows: Iterable[Mapping[str, Any]],
    fields: Sequence[str] | None = None,
) -> None:
    data = [dict(row) for row in rows]
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        ordered: list[str] = []
        for row in data:
            for key in row:
                if key not in ordered:
                    ordered.append(key)
        fields = ordered or ("status",)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(fields),
            extrasaction="raise",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in data:
            writer.writerow({key: _csv_value(row.get(key, "")) for key in fields})


def _read_csv(path: Path) -> tuple[list[dict[str, str]], tuple[str, ...]]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        rows = [dict(row) for row in reader]
        fields = tuple(reader.fieldnames or ())
    return rows, fields


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _tree_hash(path: Path) -> str:
    """Return a deterministic path-and-content digest without writing there."""

    digest = hashlib.sha256()
    if not path.exists():
        digest.update(b"MISSING")
        return digest.hexdigest().upper()
    for item in sorted(candidate for candidate in path.rglob("*") if candidate.is_file()):
        relative = item.relative_to(path).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(bytes.fromhex(_sha256(item)))
    return digest.hexdigest().upper()


def _preservation_hashes() -> dict[str, str]:
    payload = {
        path.relative_to(ROOT).as_posix(): _sha256(path) for path in FROZEN_PATHS
    }
    payload.update(
        {
            f"{path.relative_to(ROOT).as_posix()}/**": _tree_hash(path)
            for path in HISTORICAL_RESULT_DIRS
        }
    )
    return payload


def _git_state() -> dict[str, str]:
    def command(*arguments: str) -> str:
        completed = subprocess.run(
            ["git", *arguments],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
        )
        return completed.stdout.strip().replace("\\", "/")

    return {
        "top_level": command("rev-parse", "--show-toplevel"),
        "branch": command("branch", "--show-current"),
        "head": command("rev-parse", "HEAD"),
        "last_commit": command("log", "-1", "--oneline"),
        "status_short": command("status", "--short", "--untracked-files=all"),
    }


def _relative_difference(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / max(
        abs(float(left)), abs(float(right)), np.finfo(float).tiny
    )


def _relative_matrix_residual(actual: FloatArray, expected: FloatArray) -> float:
    numerator = float(np.linalg.norm(np.asarray(actual) - np.asarray(expected), ord="fro"))
    denominator = max(
        float(np.linalg.norm(np.asarray(expected), ord="fro")), np.finfo(float).tiny
    )
    return numerator / denominator


def _scaled_B_residual(laminate: rlb.LaminateSection) -> float:
    scale = max(
        float(np.linalg.norm(laminate.A, ord="fro")) * laminate.thickness,
        np.finfo(float).tiny,
    )
    return float(np.linalg.norm(laminate.B, ord="fro")) / scale


def _baseline_case() -> Any:
    _contract, selected = beta0_pilot._selected_benchmarks()
    baseline = selected["cross_ply_0_90_s"]
    actual_stack = tuple(float(ply.angle_deg) for ply in baseline.laminate.plies)
    if actual_stack != STACK_OUTER_0:
        raise RuntimeError(
            f"Frozen cross-ply factory changed: expected {STACK_OUTER_0}, got {actual_stack}"
        )
    return baseline


def build_layer_order_cases(baseline: Any | None = None) -> dict[str, LayerOrderCase]:
    """Build the two declared stacks from one frozen baseline material."""

    source = _baseline_case() if baseline is None else baseline
    material = source.laminate.plies[0].material
    thickness = float(source.laminate.thickness)
    ply_thickness = thickness / 4.0
    outer_90_laminate = rlb.integrate_laminate(
        [rlb.Ply(material, angle, ply_thickness) for angle in STACK_OUTER_90]
    )
    outer_90_properties = rlb.reduce_to_beam_properties(
        outer_90_laminate,
        width=float(source.properties.width),
        K=float(source.properties.K),
    )
    outer_90_scale = float(
        source.length**2
        * math.sqrt(
            outer_90_laminate.I0
            / outer_90_laminate.thickness**3
        )
    )
    cases = {
        "OUTER_0": LayerOrderCase(
            stack_id="OUTER_0",
            stacking_sequence=STACK_OUTER_0,
            laminate=source.laminate,
            properties=source.properties,
            total_length=float(source.length),
            frequency_scale=float(source.frequency_scale),
            K_provenance=str(source.K_provenance),
        ),
        "OUTER_90": LayerOrderCase(
            stack_id="OUTER_90",
            stacking_sequence=STACK_OUTER_90,
            laminate=outer_90_laminate,
            properties=outer_90_properties,
            total_length=float(source.length),
            frequency_scale=outer_90_scale,
            K_provenance=str(source.K_provenance),
        ),
    }
    return cases


def build_model_manifest(
    baseline: Any | None = None,
    policy: beta0_pilot.SearchPolicy | None = None,
) -> dict[str, Any]:
    """Freeze the declared finite contract before the new-stack gates."""

    source = _baseline_case() if baseline is None else baseline
    active = beta0_pilot.SearchPolicy() if policy is None else policy
    material = source.laminate.plies[0].material
    return {
        "schema_version": 1,
        "algorithm_version": ALGORITHM_VERSION,
        "stage": "RLB-2A",
        "scientific_scope": "controlled four-ply cross-ply layer-order pilot",
        "result_phrase_if_all_gates_pass": "CONTROLLED_LAYER_ORDER_PILOT_COMPLETED",
        "task_initial_git_state": TASK_INITIAL_GIT_STATE,
        "thresholds_frozen_before_full_calculation": True,
        "thresholds": THRESHOLDS,
        "root_search_policy_reused_without_change": True,
        "root_search_policy": asdict(active),
        "root_search_contract_sha256": beta0_pilot._active_contract_hash(active),
        "material": {
            "name": material.name,
            "E1": material.E1,
            "E2": material.E2,
            "E1_over_E2": material.E1 / material.E2,
            "G12": material.G12,
            "G13": material.G13,
            "G23": material.G23,
            "nu12": material.nu12,
            "rho": material.rho,
            "source": "frozen RLB-1 beta=0 cross-ply factory",
        },
        "geometry": {
            "width": source.properties.width,
            "total_thickness": source.laminate.thickness,
            "total_length": source.length,
            "arm_lengths": [0.5 * source.length, 0.5 * source.length],
            "a_over_h": source.a_over_h,
            "frequency_scale_omega_bar_per_omega": source.frequency_scale,
            "outer_clamps": "u_i(0)=w_i(0)=psi_i(0)=0",
        },
        "shear_correction": {
            "K": source.properties.K,
            "provenance": source.K_provenance,
        },
        "stacks": {
            "OUTER_0": {
                "stacking_sequence_bottom_to_top_deg": STACK_OUTER_0,
                "outer_ply_angle_deg": 0.0,
                "inner_ply_angle_deg": 90.0,
            },
            "OUTER_90": {
                "stacking_sequence_bottom_to_top_deg": STACK_OUTER_90,
                "outer_ply_angle_deg": 90.0,
                "inner_ply_angle_deg": 0.0,
            },
            "equal_ply_thickness": source.laminate.thickness / 4.0,
            "z_interfaces": [
                -0.5 * source.laminate.thickness,
                -0.25 * source.laminate.thickness,
                0.0,
                0.25 * source.laminate.thickness,
                0.5 * source.laminate.thickness,
            ],
        },
        "spectral_cases": [
            {"beta_deg": 0.0, "stack_id": stack_id} for stack_id in STACK_IDS
        ]
        + [{"beta_deg": 30.0, "stack_id": stack_id} for stack_id in STACK_IDS],
        "method_independence": {
            "one_stack_roots_used_as_other_stack_seeds": False,
            "comparison_starts_after_both_beta30_inventories_are_frozen": True,
            "primary_and_verification_scans_are_seed_free": True,
        },
        "frozen_input_hashes": _preservation_hashes(),
        "new_solver_created": False,
        "script_role": "diagnostic-only stable entry point",
        "script_proliferation_control": {
            "why_not_existing_parameter": (
                "new full-block constitutive layer-order gate and freeze-before-compare "
                "output contract"
            ),
            "reused_helpers": [
                "reddy_symmetric_laminated_beam",
                "reddy_symmetric_coupled_beams",
                "pilot_reddy_symmetric_coupled_beams_beta0.SearchPolicy",
                "pilot_reddy_symmetric_coupled_beams_beta0.seed_free_root_inventory",
            ],
        },
        "explicit_exclusions": [
            "beta_90",
            "beta_minus_30",
            "unequal_arm_lengths",
            "other_stacking_sequences",
            "material_or_K_sensitivity",
            "Rayleigh_Ritz",
            "Euler_Bernoulli",
            "FEM",
            "torsion",
            "damping",
            "mode_shapes_or_MAC",
            "branch_tracking",
            "figures",
            "article_work",
            "commit",
            "push",
        ],
    }


def laminate_property_rows(cases: Mapping[str, LayerOrderCase]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for stack_id in STACK_IDS:
        case = cases[stack_id]
        laminate = case.laminate
        properties = case.properties
        rows.append(
            {
                "stack_id": stack_id,
                "stacking_sequence": case.stacking_sequence,
                "stack_reading_direction": "bottom_to_top_z_project",
                "ply_count": len(laminate.plies),
                "ply_thicknesses": [ply.thickness for ply in laminate.plies],
                "z_interfaces": laminate.z_interfaces,
                "outer_ply_angle_deg": case.stacking_sequence[0],
                "inner_ply_angle_deg": case.stacking_sequence[1],
                "A_matrix": laminate.A,
                "B_matrix": laminate.B,
                "D_matrix": laminate.D,
                "transverse_shear_matrix_yz_xz": laminate.shear,
                "I0": laminate.I0,
                "I1": laminate.I1,
                "I2": laminate.I2,
                "Abeam": properties.A,
                "Dbeam": properties.D,
                "Sbeam": properties.S,
                "mass_per_length": properties.m,
                "rotary_inertia_per_length": properties.J,
                "width": properties.width,
                "K": properties.K,
                "K_provenance": case.K_provenance,
                "total_length": case.total_length,
                "L1": case.length_1,
                "L2": case.length_2,
                "frequency_scale": case.frequency_scale,
            }
        )
    return rows


def constitutive_gate(
    cases: Mapping[str, LayerOrderCase],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Check the complete analytic A/B/D, shear, mass, and reduction contract."""

    outer_0 = cases["OUTER_0"]
    outer_90 = cases["OUTER_90"]
    material = outer_0.laminate.plies[0].material
    h = outer_0.laminate.thickness
    Q0 = rlb.transformed_reduced_stiffness(material, 0.0)
    Q90 = rlb.transformed_reduced_stiffness(material, 90.0)
    expected_A = 0.5 * h * (Q0 + Q90)
    expected_D_outer_0 = h**3 / 96.0 * (7.0 * Q0 + Q90)
    expected_D_outer_90 = h**3 / 96.0 * (Q0 + 7.0 * Q90)
    expected_interfaces = np.array([-0.5, -0.25, 0.0, 0.25, 0.5]) * h

    rows: list[dict[str, Any]] = []

    def record(
        check_id: str,
        quantity: str,
        residual: float,
        tolerance: float,
        *,
        relation: str = "residual<=tolerance",
        passed: bool | None = None,
        details: Any = "",
    ) -> None:
        row_pass = residual <= tolerance if passed is None else bool(passed)
        rows.append(
            {
                "check_id": check_id,
                "quantity": quantity,
                "relation": relation,
                "residual_or_metric": residual,
                "tolerance": tolerance,
                "details": details,
                "status": "PASS" if row_pass else "FAIL",
            }
        )

    for stack_id in STACK_IDS:
        case = cases[stack_id]
        expected_stack = STACK_OUTER_0 if stack_id == "OUTER_0" else STACK_OUTER_90
        actual_stack = tuple(float(ply.angle_deg) for ply in case.laminate.plies)
        thickness_residual = max(
            abs(ply.thickness - h / 4.0) for ply in case.laminate.plies
        )
        record(
            f"{stack_id}_FOUR_EQUAL_PLIES",
            "ply_count_and_thickness",
            thickness_residual,
            64.0 * np.finfo(float).eps * h,
            passed=len(case.laminate.plies) == 4
            and thickness_residual <= 64.0 * np.finfo(float).eps * h,
            details={"actual_stack": actual_stack},
        )
        record(
            f"{stack_id}_ORIENTATION",
            "outer_and_inner_orientations",
            0.0 if actual_stack == expected_stack else math.inf,
            0.0,
            passed=actual_stack == expected_stack,
            details={"expected_stack": expected_stack},
        )
        interface_residual = float(
            np.max(np.abs(case.laminate.z_interfaces - expected_interfaces))
        )
        record(
            f"{stack_id}_Z_INTERFACES",
            "z_interfaces",
            interface_residual,
            64.0 * np.finfo(float).eps * h,
            details={"expected": expected_interfaces},
        )
        A_formula_residual = _relative_matrix_residual(
            case.laminate.A, expected_A
        )
        record(
            f"{stack_id}_A_ANALYTIC",
            "full_A_matrix",
            A_formula_residual,
            THRESHOLDS["A_equality_relative"],
        )
        B_residual = _scaled_B_residual(case.laminate)
        record(
            f"{stack_id}_B_ZERO",
            "scaled_B_matrix",
            B_residual,
            THRESHOLDS["B_scaled_residual"],
        )
        D_expected = expected_D_outer_0 if stack_id == "OUTER_0" else expected_D_outer_90
        D_formula_residual = _relative_matrix_residual(case.laminate.D, D_expected)
        record(
            f"{stack_id}_D_ANALYTIC",
            "full_D_matrix",
            D_formula_residual,
            THRESHOLDS["analytic_D_formula_relative"],
        )
        I1_scale = max(abs(case.laminate.I0) * h, np.finfo(float).tiny)
        record(
            f"{stack_id}_I1_ZERO",
            "I1_scaled",
            abs(case.laminate.I1) / I1_scale,
            THRESHOLDS["mass_equality_relative"],
        )

    record(
        "A_STACK_EQUALITY",
        "full_A_matrix",
        _relative_matrix_residual(outer_0.laminate.A, outer_90.laminate.A),
        THRESHOLDS["A_equality_relative"],
    )
    record(
        "SHEAR_STACK_EQUALITY",
        "full_transverse_shear_matrix",
        _relative_matrix_residual(
            outer_0.laminate.shear, outer_90.laminate.shear
        ),
        THRESHOLDS["shear_equality_relative"],
    )
    for name in ("I0", "I2"):
        record(
            f"{name}_STACK_EQUALITY",
            name,
            _relative_difference(
                getattr(outer_0.laminate, name), getattr(outer_90.laminate, name)
            ),
            THRESHOLDS["mass_equality_relative"],
        )
    for name, attribute in (
        ("ABEAM", "A"),
        ("SBEAM", "S"),
        ("MASS_PER_LENGTH", "m"),
        ("ROTARY_INERTIA", "J"),
    ):
        record(
            f"{name}_STACK_EQUALITY",
            attribute,
            _relative_difference(
                getattr(outer_0.properties, attribute),
                getattr(outer_90.properties, attribute),
            ),
            THRESHOLDS["reduced_property_equality_relative"],
        )

    D_difference = float(
        np.linalg.norm(outer_0.laminate.D - outer_90.laminate.D, ord="fro")
    )
    record(
        "D_MATRICES_DIFFER",
        "full_D_matrix_difference_norm",
        D_difference,
        0.0,
        relation="metric>0",
        passed=D_difference > 0.0,
    )
    Dbeam_difference = outer_0.properties.D - outer_90.properties.D
    record(
        "DBEAM_OUTER_0_GREATER",
        "Dbeam_OUTER_0-Dbeam_OUTER_90",
        Dbeam_difference,
        0.0,
        relation="metric>0",
        passed=Dbeam_difference > 0.0,
    )

    passed = all(row["status"] == "PASS" for row in rows)
    return rows, {
        "pass": passed,
        "Q0": Q0,
        "Q90": Q90,
        "expected_A": expected_A,
        "expected_D_OUTER_0": expected_D_outer_0,
        "expected_D_OUTER_90": expected_D_outer_90,
        "maximum_A_formula_or_equality_relative": max(
            float(row["residual_or_metric"])
            for row in rows
            if row["check_id"].endswith("A_ANALYTIC")
            or row["check_id"] == "A_STACK_EQUALITY"
        ),
        "maximum_D_formula_relative": max(
            float(row["residual_or_metric"])
            for row in rows
            if row["check_id"].endswith("D_ANALYTIC")
        ),
        "maximum_scaled_B_residual": max(
            float(row["residual_or_metric"])
            for row in rows
            if row["check_id"].endswith("B_ZERO")
        ),
        "Dbeam_OUTER_0": outer_0.properties.D,
        "Dbeam_OUTER_90": outer_90.properties.D,
        "Dbeam_ratio_OUTER_0_over_OUTER_90": (
            outer_0.properties.D / outer_90.properties.D
        ),
    }


def _direct_family_inventory(
    case: LayerOrderCase,
    policy: beta0_pilot.SearchPolicy,
) -> DirectFamilyInventory:
    """Call the frozen RLB-0 family solvers and RLB-1 cluster helpers once."""

    scale = case.frequency_scale
    axial = tuple(
        rlb.exact_axial_modes(
            case.properties, case.total_length, "FF", n_modes=24
        )
    )
    bending_detections = rlb.find_bending_roots(
        case.properties,
        case.total_length,
        "CC",
        omega_max=policy.omega_bar_max / scale,
        n_roots=None,
        omega_min=0.0,
        scan_points=policy.verification_scan_points,
        sigma_ratio_tolerance=policy.sigma_ratio_tolerance,
        root_xtol=max(policy.root_xtol_bar / scale, 1.0e-15),
        dedup_rtol=1.0e-12,
    )
    bending_groups, reconciliation = beta0_pilot._reconcile_reference_detector_detections(
        bending_detections,
        frequency_scale=scale,
        policy=policy,
    )
    bending = tuple(
        min(
            group,
            key=lambda item: (
                0 if item.detection.startswith("determinant") else 1,
                item.boundary_residual,
                item.sigma_ratio,
            ),
        )
        for group in bending_groups
    )[:24]
    clusters = rlb.union_subsystem_spectra(
        [item.omega for item in axial],
        [item.omega for item in bending],
        atol=policy.cluster_atol_bar / scale,
        rtol=policy.cluster_rtol,
    )
    union_rows: list[dict[str, Any]] = []
    sorted_slot = 0
    for cluster_index, cluster in enumerate(clusters, start=1):
        if sorted_slot >= policy.required_slots and cluster.multiplicity == 1:
            break
        semantics = (
            "ISOLATED"
            if cluster.multiplicity == 1
            else (
                "EXACT_DEGENERATE_SUBSPACE"
                if cluster.exact_degeneracy
                else "NEAR_DEGENERATE_CLUSTER"
            )
        )
        family_content = "+".join(
            sorted(member.subsystem for member in cluster.members)
        )
        for cluster_slot, member in enumerate(cluster.members, start=1):
            sorted_slot += 1
            union_rows.append(
                {
                    "sorted_slot": sorted_slot,
                    "role": (
                        "FIRST_12"
                        if sorted_slot <= 12
                        else (
                            "ROOT_13_GUARD"
                            if sorted_slot == 13
                            else "GUARD_CLUSTER_COMPLETION"
                        )
                    ),
                    "cluster_index": cluster_index,
                    "cluster_semantics": semantics,
                    "cluster_multiplicity": cluster.multiplicity,
                    "cluster_total_nullity": cluster.multiplicity,
                    "cluster_center_omega_bar": cluster.representative_omega * scale,
                    "cluster_member_slot": cluster_slot,
                    "omega": member.omega,
                    "omega_bar": member.omega * scale,
                    "family": member.subsystem,
                    "family_index": member.subsystem_index,
                    "family_content": family_content,
                }
            )
        if sorted_slot >= policy.required_slots:
            break
    return DirectFamilyInventory(
        axial=axial,
        bending=bending,
        union_rows=tuple(union_rows),
        reconciliation=reconciliation,
    )


def _axial_rows(
    case: LayerOrderCase, family: DirectFamilyInventory
) -> list[dict[str, Any]]:
    through_guard = {
        int(row["family_index"]) + 1
        for row in family.union_rows
        if row["family"] == "axial"
    }
    return [
        {
            "stack_id": case.stack_id,
            "family": "axial",
            "family_index": mode.n,
            "omega": mode.omega,
            "omega_bar": mode.omega * case.frequency_scale,
            "wavenumber": mode.wavenumber,
            "boundary_condition": mode.boundary_condition,
            "represented_through_union_guard": mode.n in through_guard,
            "solver": "rlb.exact_axial_modes",
        }
        for mode in family.axial
    ]


def _bending_rows(
    case: LayerOrderCase, family: DirectFamilyInventory
) -> list[dict[str, Any]]:
    through_guard = {
        int(row["family_index"]) + 1
        for row in family.union_rows
        if row["family"] == "bending"
    }
    rows: list[dict[str, Any]] = []
    for index, root in enumerate(family.bending, start=1):
        rows.append(
            {
                "stack_id": case.stack_id,
                "family": "bending",
                "family_index": index,
                "omega": root.omega,
                "omega_bar": root.omega * case.frequency_scale,
                "determinant": root.determinant,
                "sigma_min": root.sigma_min,
                "sigma_max": root.sigma_max,
                "sigma_min_over_sigma_max": root.sigma_ratio,
                "boundary_null_residual": root.boundary_residual,
                "detector_type": root.detection,
                "root_status": "PASS" if root.accepted else "FAIL",
                "represented_through_union_guard": index in through_guard,
                "solver": "rlb.find_bending_roots",
            }
        )
    return rows


def _case_spec(case: LayerOrderCase) -> beta0_pilot.CoupledCaseSpec:
    return beta0_pilot.CoupledCaseSpec(
        case_id=f"{case.stack_id.lower()}__equal",
        category="homogeneous_layer_order",
        split_id="equal",
        order_id=f"{case.stack_id}|{case.stack_id}",
        length_1=case.length_1,
        properties_1=case.properties,
        length_2=case.length_2,
        properties_2=case.properties,
        total_length=case.total_length,
        frequency_scale=case.frequency_scale,
    )


def _beta30_provider(case: LayerOrderCase) -> Any:
    beta_rad = math.radians(30.0)

    def provider(omega: float) -> FloatArray:
        return np.asarray(
            coupled.coupled_boundary_matrix(
                omega,
                beta_rad,
                case.length_1,
                case.properties,
                case.length_2,
                case.properties,
            ),
            dtype=float,
        )

    return provider


def _root_row(
    case: LayerOrderCase,
    beta_deg: float,
    inventory: beta0_pilot.RootInventory,
    slot: beta0_pilot.SpectrumSlot,
) -> dict[str, Any]:
    spec = _case_spec(case)
    row = beta0_pilot._slot_row(spec, inventory, slot)
    candidate = slot.event.candidate
    diagnostic = candidate.diagnostics
    row.update(
        {
            "case_id": inventory.case_id,
            "stack_id": case.stack_id,
            "stacking_sequence": case.stacking_sequence,
            "beta_deg": beta_deg,
            "project_nondimensional_frequency": slot.event.omega_bar,
            "bracket_left_omega_bar": candidate.interval_left_bar,
            "bracket_right_omega_bar": candidate.interval_right_bar,
            "detector_type": candidate.detection_sources,
            "root_status": (
                "PASS"
                if diagnostic.scaled_sigma_ratio
                <= THRESHOLDS["root_singular_ratio"]
                and diagnostic.raw_boundary_null_residual
                <= THRESHOLDS["boundary_residual"]
                else "FAIL"
            ),
            "guard_flag": slot.role != "FIRST_12",
            "inventory_sha256": inventory.inventory_sha256,
            "primary_verification_max_relative": (
                inventory.maximum_primary_verification_relative
            ),
            "unresolved_candidates_below_guard": (
                inventory.unresolved_low_sigma_count
            ),
        }
    )
    return row


def _inventory_quality(
    inventory: beta0_pilot.RootInventory,
    policy: beta0_pilot.SearchPolicy,
) -> dict[str, Any]:
    slots = inventory.slots
    maximum_sigma = max(
        (slot.event.candidate.diagnostics.scaled_sigma_ratio for slot in slots),
        default=math.inf,
    )
    maximum_boundary = max(
        (
            slot.event.candidate.diagnostics.raw_boundary_null_residual
            for slot in slots
        ),
        default=math.inf,
    )
    primary_count = len(inventory.primary.slots)
    verification_count = len(inventory.verification.slots)
    passed = bool(
        inventory.status == "PASS"
        and len(slots) == policy.required_slots
        and primary_count == policy.required_slots
        and verification_count == policy.required_slots
        and inventory.guard_available
        and inventory.guard_not_at_scan_boundary
        and inventory.unresolved_low_sigma_count == 0
        and maximum_sigma <= THRESHOLDS["root_singular_ratio"]
        and maximum_boundary <= THRESHOLDS["boundary_residual"]
    )
    return {
        "case_id": inventory.case_id,
        "builder_id": inventory.builder_id,
        "inventory_sha256": inventory.inventory_sha256,
        "inventory_status": inventory.status,
        "slot_count": len(slots),
        "primary_slot_count": primary_count,
        "verification_slot_count": verification_count,
        "guard_available": inventory.guard_available,
        "guard_not_at_scan_boundary": inventory.guard_not_at_scan_boundary,
        "unresolved_candidates_below_guard": inventory.unresolved_low_sigma_count,
        "maximum_primary_verification_relative": (
            inventory.maximum_primary_verification_relative
        ),
        "maximum_scaled_sigma_ratio": maximum_sigma,
        "maximum_boundary_null_residual": maximum_boundary,
        "pass": passed,
    }


def _family_comparison_rows(
    cases: Mapping[str, LayerOrderCase],
    families: Mapping[str, DirectFamilyInventory],
    beta0_inventories: Mapping[str, beta0_pilot.RootInventory],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    outer_0 = cases["OUTER_0"]
    outer_90 = cases["OUTER_90"]
    family_0 = families["OUTER_0"]
    family_90 = families["OUTER_90"]

    axial_pass = len(family_0.axial) == len(family_90.axial) >= 13
    maximum_axial_relative = 0.0
    for left, right in zip(family_0.axial, family_90.axial, strict=True):
        omega_0 = left.omega * outer_0.frequency_scale
        omega_90 = right.omega * outer_90.frequency_scale
        relative = _relative_difference(omega_0, omega_90)
        maximum_axial_relative = max(maximum_axial_relative, relative)
        passed = relative <= THRESHOLDS["beta0_axial_relative"]
        axial_pass = axial_pass and passed
        rows.append(
            {
                "comparison_kind": "OUTER_0_vs_OUTER_90_physical_family",
                "family": "axial",
                "family_index": left.n,
                "omega_bar_OUTER_0": omega_0,
                "omega_bar_OUTER_90": omega_90,
                "relative_difference": relative,
                "delta_OUTER_0_over_OUTER_90": (omega_0 - omega_90) / omega_90,
                "relation": "equal",
                "tolerance": THRESHOLDS["beta0_axial_relative"],
                "status": "PASS" if passed else "FAIL",
            }
        )

    bending_pass = len(family_0.bending) == len(family_90.bending) >= 13
    minimum_bending_shift = math.inf
    maximum_bending_shift = -math.inf
    for index, (left, right) in enumerate(
        zip(family_0.bending, family_90.bending, strict=True), start=1
    ):
        omega_0 = left.omega * outer_0.frequency_scale
        omega_90 = right.omega * outer_90.frequency_scale
        scale = max(abs(omega_0), abs(omega_90), 1.0)
        allowance = THRESHOLDS["beta0_bending_order_allowance_relative"] * scale
        delta = (omega_0 - omega_90) / omega_90
        passed = omega_0 + allowance > omega_90
        bending_pass = bending_pass and passed
        minimum_bending_shift = min(minimum_bending_shift, delta)
        maximum_bending_shift = max(maximum_bending_shift, delta)
        rows.append(
            {
                "comparison_kind": "OUTER_0_vs_OUTER_90_physical_family",
                "family": "bending",
                "family_index": index,
                "omega_bar_OUTER_0": omega_0,
                "omega_bar_OUTER_90": omega_90,
                "relative_difference": _relative_difference(omega_0, omega_90),
                "delta_OUTER_0_over_OUTER_90": delta,
                "relation": "OUTER_0>OUTER_90_with_numerical_allowance",
                "tolerance": THRESHOLDS[
                    "beta0_bending_order_allowance_relative"
                ],
                "status": "PASS" if passed else "FAIL",
            }
        )

    union_pass = True
    union_maximum = 0.0
    union_summaries: dict[str, Any] = {}
    for stack_id in STACK_IDS:
        comparison, summary = beta0_pilot._compare_inventory_to_union(
            beta0_inventories[stack_id],
            families[stack_id].union_rows,
            case_id=f"{stack_id.lower()}__equal",
        )
        union_summaries[stack_id] = summary
        union_pass = union_pass and bool(summary["pass"])
        union_maximum = max(
            union_maximum, float(summary["maximum_relative_difference"])
        )
        for row in comparison:
            converted = dict(row)
            converted["stack_id"] = stack_id
            converted["comparison_kind"] = "coupled_beta0_vs_axial_bending_union"
            converted["tolerance"] = THRESHOLDS["beta0_coupled_union_relative"]
            rows.append(converted)

    return rows, {
        "pass": axial_pass and bending_pass and union_pass,
        "axial_pass": axial_pass,
        "bending_pass": bending_pass,
        "coupled_union_pass": union_pass,
        "maximum_axial_relative_difference": maximum_axial_relative,
        "minimum_bending_relative_shift": minimum_bending_shift,
        "maximum_bending_relative_shift": maximum_bending_shift,
        "maximum_coupled_union_relative_difference": union_maximum,
        "coupled_union": union_summaries,
    }


def _neighbor_gap(slots: Sequence[beta0_pilot.SpectrumSlot], index: int) -> tuple[Any, Any]:
    left = (
        ""
        if index == 0
        else slots[index].event.omega_bar - slots[index - 1].event.omega_bar
    )
    right = (
        ""
        if index + 1 >= len(slots)
        else slots[index + 1].event.omega_bar - slots[index].event.omega_bar
    )
    return left, right


def _beta30_shift_rows(
    inventories: Mapping[str, beta0_pilot.RootInventory],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    outer_0 = inventories["OUTER_0"].slots
    outer_90 = inventories["OUTER_90"].slots
    rows: list[dict[str, Any]] = []
    passed = len(outer_0) == len(outer_90) == 13
    minimum_shift = math.inf
    maximum_shift = -math.inf
    maximum_quality_ratio = 0.0
    for index in range(max(len(outer_0), len(outer_90))):
        if index >= len(outer_0) or index >= len(outer_90):
            rows.append(
                {
                    "sorted_position": index + 1,
                    "status": "FAIL",
                    "reason": "MISSING_SORTED_SLOT",
                }
            )
            passed = False
            continue
        slot_0 = outer_0[index]
        slot_90 = outer_90[index]
        omega_0 = slot_0.event.omega_bar
        omega_90 = slot_90.event.omega_bar
        absolute = omega_0 - omega_90
        relative = _relative_difference(omega_0, omega_90)
        delta = absolute / omega_90
        scale = max(abs(omega_0), abs(omega_90), 1.0)
        allowance = THRESHOLDS["beta30_ordered_frequency_allowance_relative"] * scale
        monotonic = omega_0 + allowance >= omega_90
        quality_0 = slot_0.event.candidate.diagnostics
        quality_90 = slot_90.event.candidate.diagnostics
        quality_ratio = max(quality_0.scaled_sigma_ratio, quality_90.scaled_sigma_ratio)
        boundary = max(
            quality_0.raw_boundary_null_residual,
            quality_90.raw_boundary_null_residual,
        )
        row_pass = bool(
            monotonic
            and quality_ratio <= THRESHOLDS["root_singular_ratio"]
            and boundary <= THRESHOLDS["boundary_residual"]
        )
        passed = passed and row_pass
        minimum_shift = min(minimum_shift, delta)
        maximum_shift = max(maximum_shift, delta)
        maximum_quality_ratio = max(maximum_quality_ratio, quality_ratio)
        left_0, right_0 = _neighbor_gap(outer_0, index)
        left_90, right_90 = _neighbor_gap(outer_90, index)
        rows.append(
            {
                "sorted_position": index + 1,
                "role": slot_0.role,
                "omega_bar_OUTER_0": omega_0,
                "omega_bar_OUTER_90": omega_90,
                "absolute_difference_OUTER_0_minus_OUTER_90": absolute,
                "relative_difference": relative,
                "delta_k": delta,
                "left_neighbor_gap_OUTER_0": left_0,
                "right_neighbor_gap_OUTER_0": right_0,
                "left_neighbor_gap_OUTER_90": left_90,
                "right_neighbor_gap_OUTER_90": right_90,
                "root_quality_scaled_sigma_ratio_OUTER_0": (
                    quality_0.scaled_sigma_ratio
                ),
                "root_quality_scaled_sigma_ratio_OUTER_90": (
                    quality_90.scaled_sigma_ratio
                ),
                "boundary_null_residual_OUTER_0": (
                    quality_0.raw_boundary_null_residual
                ),
                "boundary_null_residual_OUTER_90": (
                    quality_90.raw_boundary_null_residual
                ),
                "cluster_id_OUTER_0": slot_0.event.cluster_id,
                "cluster_id_OUTER_90": slot_90.event.cluster_id,
                "guard_flag": index == 12,
                "comparison_semantics": (
                    "ordered_eigenvalue_only_not_modal_descendant"
                ),
                "monotonicity_allowance_relative": THRESHOLDS[
                    "beta30_ordered_frequency_allowance_relative"
                ],
                "status": "PASS" if row_pass else "FAIL",
            }
        )
    return rows, {
        "pass": passed,
        "minimum_relative_shift": minimum_shift,
        "maximum_relative_shift": maximum_shift,
        "maximum_root_quality_scaled_sigma_ratio": maximum_quality_ratio,
        "slot_count_OUTER_0": len(outer_0),
        "slot_count_OUTER_90": len(outer_90),
    }


def _status(
    scientific_pass: bool,
    inventory_pass: bool,
    *,
    allow_partial: bool = True,
) -> str:
    if scientific_pass and inventory_pass:
        return "PASS"
    if scientific_pass and allow_partial:
        return "PARTIAL_PASS"
    return "FAIL"


def _report_text(summary: Mapping[str, Any]) -> str:
    constitutive = summary["constitutive"]
    beta0 = summary["beta0"]
    beta30 = summary["beta30"]
    quality = summary["root_inventory"]
    statuses = summary["statuses"]
    return f"""# RLB-2A: порядок слоёв четырёхслойного cross-ply-ламината

## Объём проверки

Сопоставлены ровно две симметричные укладки из одного и того же однонаправленного материала: `OUTER_0 = [0/90/90/0]` и `OUTER_90 = [90/0/0/90]`. Все четыре слоя имеют толщину `h/4`; использованы `E2=rho=b=h=1`, `E1/E2=25`, `G12=G13=0.5 E2`, `G23=0.2 E2`, `nu12=0.25`, `L_total=20`, `L1=L2=10` и прежний `K=5/6`. Рассмотрены только `beta=0` и `beta=30 deg`.

## Constitutive gate

Границы слоёв равны `[-h/2,-h/4,0,h/4,h/2]`. Для полного in-plane блока численно проверены формулы

`A = h/2 (Q0+Q90)`,

`D_OUTER_0 = h^3/96 (7 Q0+Q90)`,

`D_OUTER_90 = h^3/96 (Q0+7 Q90)`.

Максимальная относительная невязка аналитической формулы `D` составила `{constitutive['maximum_D_formula_relative']:.6e}`, максимальная масштабированная невязка `B=0` — `{constitutive['maximum_scaled_B_residual']:.6e}`. Матрицы `A`, transverse-shear stiffness и массовые моменты совпадают. После free-edge reduction совпадают `Abeam`, `Sbeam`, `m` и `J`; различается только `Dbeam`. Получено

`Dbeam_OUTER_0 = {constitutive['Dbeam_OUTER_0']:.15g}`,

`Dbeam_OUTER_90 = {constitutive['Dbeam_OUTER_90']:.15g}`,

`Dbeam_OUTER_0 / Dbeam_OUTER_90 = {constitutive['Dbeam_ratio_OUTER_0_over_OUTER_90']:.15g}`.

## Прямой предел beta=0

Для каждой укладки объединены точное продольное fixed--fixed семейство и проверенное изгибное clamped--clamped семейство RLB-0. Это объединение сопоставлено с coupled determinant inventory RLB-1. Максимальное расхождение coupled/union равно `{beta0['maximum_coupled_union_relative_difference']:.6e}`. Максимальное расхождение продольных частот двух укладок равно `{beta0['maximum_axial_relative_difference']:.6e}`. Для проверенных изгибных членов относительный прирост при `OUTER_0` лежит в диапазоне от `{beta0['minimum_bending_relative_shift']:.6e}` до `{beta0['maximum_bending_relative_shift']:.6e}`.

## Угловой pilot beta=30 deg

Оба seed-free inventories вычислены независимо и заморожены до сопоставления; частоты одной укладки не использовались как seeds другой. Для каждой укладки сохранены первые 12 положительных roots и root 13 guard. Минимальный и максимальный ordered shift `(omega_OUTER_0-omega_OUTER_90)/omega_OUTER_90` равны `{beta30['minimum_relative_shift']:.6e}` и `{beta30['maximum_relative_shift']:.6e}`. Все 13 упорядоченных частот `OUTER_0` не ниже соответствующих частот `OUTER_90` с заранее заданным относительным allowance `1e-10`. Совпадающий sorted index трактуется только как позиция упорядоченного собственного значения, а не как modal descendant.

## Полнота и качество roots

Проверены четыре coupled inventories: два при `beta=0` и два при `beta=30 deg`. Во всех случаях требуются ровно 13 slots, согласованные primary/verification inventories, отсутствие неразрешённых кандидатов ниже guard, `sigma_min/sigma_max <= 1e-9` и boundary null residual `<=1e-9`. Максимальные значения этих диагностик составили соответственно `{quality['maximum_scaled_sigma_ratio']:.6e}` и `{quality['maximum_boundary_null_residual']:.6e}`.

## Статусы

- `{STATUS_CONSTITUTIVE}: {statuses[STATUS_CONSTITUTIVE]}`
- `{STATUS_BETA0}: {statuses[STATUS_BETA0]}`
- `{STATUS_BETA30}: {statuses[STATUS_BETA30]}`
- `{STATUS_INVENTORY}: {statuses[STATUS_INVENTORY]}`
- `{STATUS_OVERALL}: {statuses[STATUS_OVERALL]}`

Итоговая формулировка при полном прохождении gates: `{summary['result_phrase']}`.

## Ограничения

Результат относится только к объявленному конечному pilot для двух четырёхслойных симметричных cross-ply-укладок и углов `0` и `30 deg`. Он не является parameter study и не переносится автоматически на произвольный `beta`, другие материалы, толщины, значения `K` или несимметричные укладки. Формы, MAC, branch tracking, veering, localization, Ritz, FEM, Euler--Bernoulli, torsion, damping и figures не вычислялись.
"""


def repair_generated_metadata(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    """Repair two CSV provenance fields without any physical calculation.

    The repair maps zero-based subsystem indices to the one-based public family
    indices and makes root-table case IDs identical to the candidate-table
    inventory IDs.  A digest excluding only those two metadata fields proves
    that frequencies and all root diagnostics are unchanged.
    """

    target = Path(output_dir)
    comparison_rows, _comparison_fields = _read_csv(
        target / "beta0_family_comparison.csv"
    )
    union_frequencies: dict[tuple[str, str], list[float]] = {}
    for row in comparison_rows:
        if row.get("comparison_kind") != "coupled_beta0_vs_axial_bending_union":
            continue
        key = (str(row["stack_id"]), str(row["family_content"]))
        union_frequencies.setdefault(key, []).append(
            float(row["union_omega_bar_or_center"])
        )

    affected = (
        "beta0_axial_roots.csv",
        "beta0_bending_roots.csv",
        "beta0_coupled_roots.csv",
        "beta30_roots.csv",
    )

    def scientific_digest() -> str:
        payload: list[dict[str, Any]] = []
        for filename in affected:
            rows, _fields = _read_csv(target / filename)
            for row in rows:
                payload.append(
                    {
                        "file": filename,
                        **{
                            key: value
                            for key, value in row.items()
                            if key not in {"case_id", "represented_through_union_guard"}
                        },
                    }
                )
        return hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode(
                "utf-8"
            )
        ).hexdigest().upper()

    before_digest = scientific_digest()
    repair_counts = {"family_guard_flags": 0, "root_case_ids": 0}
    for filename, family in (
        ("beta0_axial_roots.csv", "axial"),
        ("beta0_bending_roots.csv", "bending"),
    ):
        rows, fields = _read_csv(target / filename)
        for row in rows:
            omega_bar = float(row["omega_bar"])
            references = union_frequencies[(str(row["stack_id"]), family)]
            represented = any(
                _relative_difference(omega_bar, reference)
                <= THRESHOLDS["beta0_coupled_union_relative"]
                for reference in references
            )
            expected = "True" if represented else "False"
            if row["represented_through_union_guard"] != expected:
                repair_counts["family_guard_flags"] += 1
            row["represented_through_union_guard"] = expected
        _write_csv(target / filename, rows, fields)

    for filename in ("beta0_coupled_roots.csv", "beta30_roots.csv"):
        rows, fields = _read_csv(target / filename)
        for row in rows:
            beta_label = int(round(float(row["beta_deg"])))
            expected = f"beta{beta_label}__{str(row['stack_id']).lower()}"
            if row["case_id"] != expected:
                repair_counts["root_case_ids"] += 1
            row["case_id"] = expected
        _write_csv(target / filename, rows, fields)

    after_digest = scientific_digest()
    if before_digest != after_digest:
        raise RuntimeError("metadata repair changed scientific CSV values")

    manifest_path = target / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["closing_metadata_repair"] = {
        "status": "PASS",
        "reason": [
            "ZERO_BASED_SUBSYSTEM_INDEX_TO_ONE_BASED_OUTPUT_FLAG",
            "ROOT_TABLE_CASE_ID_TO_INVENTORY_CASE_ID",
        ],
        "repair_counts": repair_counts,
        "scientific_payload_sha256_before": before_digest,
        "scientific_payload_sha256_after": after_digest,
        "scientific_values_changed": False,
        "global_root_searches_run": 0,
        "direct_family_searches_run": 0,
        "physical_matrices_evaluated": 0,
    }
    manifest["final_output_hashes_after_metadata_repair"] = {
        filename: _sha256(target / filename) for filename in affected
    }
    _write_json(manifest_path, manifest)
    return manifest["closing_metadata_repair"]


def run_pilot(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    policy: beta0_pilot.SearchPolicy | None = None,
) -> dict[str, Any]:
    active = beta0_pilot.SearchPolicy() if policy is None else policy
    target = Path(output_dir)
    baseline = _baseline_case()
    model_manifest = build_model_manifest(baseline, active)
    _write_json(target / "model_manifest.json", model_manifest)

    run_start_git = _git_state()
    preservation_before = _preservation_hashes()
    cases = build_layer_order_cases(baseline)
    property_rows = laminate_property_rows(cases)
    constitutive_rows, constitutive_summary = constitutive_gate(cases)
    _write_csv(target / "laminate_properties.csv", property_rows)
    _write_csv(target / "analytic_stiffness_checks.csv", constitutive_rows)

    if not constitutive_summary["pass"]:
        raise RuntimeError(
            "RLB-2A constitutive gate failed; spectral calculations were not started"
        )

    families = {
        stack_id: _direct_family_inventory(cases[stack_id], active)
        for stack_id in STACK_IDS
    }
    _write_csv(
        target / "beta0_axial_roots.csv",
        [
            row
            for stack_id in STACK_IDS
            for row in _axial_rows(cases[stack_id], families[stack_id])
        ],
    )
    _write_csv(
        target / "beta0_bending_roots.csv",
        [
            row
            for stack_id in STACK_IDS
            for row in _bending_rows(cases[stack_id], families[stack_id])
        ],
    )

    beta0_inventories: dict[str, beta0_pilot.RootInventory] = {}
    for stack_id in STACK_IDS:
        case = cases[stack_id]
        spec = _case_spec(case)
        beta0_inventories[stack_id] = beta0_pilot.seed_free_root_inventory(
            beta0_pilot._coupled_provider(coupled, spec),
            case.frequency_scale,
            active,
            case_id=f"beta0__{stack_id.lower()}",
            builder_id="frozen_RLB1_coupled_beta0",
        )
    beta0_root_rows = [
        _root_row(case=cases[stack_id], beta_deg=0.0, inventory=inventory, slot=slot)
        for stack_id, inventory in beta0_inventories.items()
        for slot in inventory.slots
    ]
    _write_csv(target / "beta0_coupled_roots.csv", beta0_root_rows)
    beta0_comparison_rows, beta0_summary = _family_comparison_rows(
        cases, families, beta0_inventories
    )
    _write_csv(target / "beta0_family_comparison.csv", beta0_comparison_rows)

    beta30_inventories: dict[str, beta0_pilot.RootInventory] = {}
    for stack_id in STACK_IDS:
        case = cases[stack_id]
        beta30_inventories[stack_id] = beta0_pilot.seed_free_root_inventory(
            _beta30_provider(case),
            case.frequency_scale,
            active,
            case_id=f"beta30__{stack_id.lower()}",
            builder_id="frozen_RLB1_coupled_beta30",
        )
    beta30_inventory_hashes = {
        stack_id: inventory.inventory_sha256
        for stack_id, inventory in beta30_inventories.items()
    }
    beta30_root_rows = [
        _root_row(case=cases[stack_id], beta_deg=30.0, inventory=inventory, slot=slot)
        for stack_id, inventory in beta30_inventories.items()
        for slot in inventory.slots
    ]
    _write_csv(target / "beta30_roots.csv", beta30_root_rows)
    beta30_roots_file_sha256_before_comparison = _sha256(target / "beta30_roots.csv")

    beta30_shift_rows, beta30_summary = _beta30_shift_rows(beta30_inventories)
    _write_csv(target / "beta30_spectrum_shift.csv", beta30_shift_rows)

    all_inventories = [
        *[beta0_inventories[stack_id] for stack_id in STACK_IDS],
        *[beta30_inventories[stack_id] for stack_id in STACK_IDS],
    ]
    candidate_rows, rejected_rows = beta0_pilot._all_inventory_candidates(
        all_inventories
    )
    _write_csv(target / "root_candidates.csv", candidate_rows)
    _write_csv(target / "rejected_candidates.csv", rejected_rows)

    quality_rows = [_inventory_quality(inventory, active) for inventory in all_inventories]
    inventory_pass = all(bool(row["pass"]) for row in quality_rows)
    maximum_sigma = max(float(row["maximum_scaled_sigma_ratio"]) for row in quality_rows)
    maximum_boundary = max(
        float(row["maximum_boundary_null_residual"]) for row in quality_rows
    )
    beta0_inventory_pass = all(
        bool(_inventory_quality(inventory, active)["pass"])
        for inventory in beta0_inventories.values()
    )
    beta30_inventory_pass = all(
        bool(_inventory_quality(inventory, active)["pass"])
        for inventory in beta30_inventories.values()
    )
    statuses = {
        STATUS_CONSTITUTIVE: "PASS" if constitutive_summary["pass"] else "FAIL",
        STATUS_BETA0: _status(
            bool(beta0_summary["pass"]), beta0_inventory_pass
        ),
        STATUS_BETA30: _status(
            bool(beta30_summary["pass"]), beta30_inventory_pass
        ),
        STATUS_INVENTORY: "PASS" if inventory_pass else "FAIL",
    }
    statuses[STATUS_OVERALL] = (
        "PASS"
        if all(statuses[name] == "PASS" for name in (
            STATUS_CONSTITUTIVE,
            STATUS_BETA0,
            STATUS_BETA30,
            STATUS_INVENTORY,
        ))
        else (
            "PARTIAL_PASS"
            if not any(statuses[name] == "FAIL" for name in (
                STATUS_CONSTITUTIVE,
                STATUS_BETA0,
                STATUS_BETA30,
                STATUS_INVENTORY,
            ))
            else "FAIL"
        )
    )
    result_phrase = (
        "CONTROLLED_LAYER_ORDER_PILOT_COMPLETED"
        if statuses[STATUS_OVERALL] == "PASS"
        else "CONTROLLED_LAYER_ORDER_PILOT_NOT_COMPLETED"
    )

    preservation_after = _preservation_hashes()
    preservation_pass = preservation_before == preservation_after
    summary = {
        "constitutive": constitutive_summary,
        "beta0": beta0_summary,
        "beta30": beta30_summary,
        "root_inventory": {
            "pass": inventory_pass,
            "case_count": len(quality_rows),
            "cases": quality_rows,
            "maximum_scaled_sigma_ratio": maximum_sigma,
            "maximum_boundary_null_residual": maximum_boundary,
        },
        "statuses": statuses,
        "result_phrase": result_phrase,
    }
    (target / "report.md").write_text(_report_text(summary), encoding="utf-8")

    run_manifest = {
        "schema_version": 1,
        "algorithm_version": ALGORITHM_VERSION,
        "stage": "RLB-2A",
        "task_initial_git_state": TASK_INITIAL_GIT_STATE,
        "execution_start_git_state": run_start_git,
        "execution_end_git_state": _git_state(),
        "model_manifest_sha256": _sha256(target / "model_manifest.json"),
        "thresholds": THRESHOLDS,
        "root_search_policy": asdict(active),
        "root_search_policy_changed": False,
        "beta0_inventory_hashes": {
            stack_id: inventory.inventory_sha256
            for stack_id, inventory in beta0_inventories.items()
        },
        "beta30_inventory_hashes_frozen_before_comparison": beta30_inventory_hashes,
        "beta30_roots_file_sha256_frozen_before_comparison": (
            beta30_roots_file_sha256_before_comparison
        ),
        "one_stack_roots_used_as_other_stack_seeds": False,
        "inventory_quality": quality_rows,
        "direct_family_reconciliation": {
            stack_id: dict(families[stack_id].reconciliation)
            for stack_id in STACK_IDS
        },
        "summary": summary,
        "statuses": statuses,
        "result_phrase": result_phrase,
        "preservation_before": preservation_before,
        "preservation_after": preservation_after,
        "preservation_pass": preservation_pass,
        "generated_files": [
            "model_manifest.json",
            "laminate_properties.csv",
            "analytic_stiffness_checks.csv",
            "beta0_axial_roots.csv",
            "beta0_bending_roots.csv",
            "beta0_coupled_roots.csv",
            "beta0_family_comparison.csv",
            "beta30_roots.csv",
            "beta30_spectrum_shift.csv",
            "root_candidates.csv",
            "rejected_candidates.csv",
            "report.md",
            "run_manifest.json",
        ],
        "figures_created": 0,
        "new_solver_modules_created": 0,
        "restricted_operations": {
            "parameter_sweeps": 0,
            "Rayleigh_Ritz_solves": 0,
            "FEM_solves": 0,
            "Euler_Bernoulli_solves": 0,
            "torsion_or_damping_solves": 0,
            "mode_shape_or_MAC_calculations": 0,
            "article_work": 0,
            "commits": 0,
            "pushes": 0,
        },
    }
    _write_json(target / "run_manifest.json", run_manifest)
    return run_manifest


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--repair-metadata-only",
        action="store_true",
        help="repair existing CSV provenance fields without running any roots",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_args(argv)
    if arguments.repair_metadata_only:
        repair = repair_generated_metadata(arguments.output_dir)
        print(json.dumps(repair, ensure_ascii=False, indent=2))
        return 0
    manifest = run_pilot(arguments.output_dir)
    print(json.dumps(manifest["statuses"], ensure_ascii=False, indent=2))
    return 0 if manifest["statuses"][STATUS_OVERALL] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
