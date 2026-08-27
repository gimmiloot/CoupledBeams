"""RLB-2C finite beta-grid comparison in the historical frequency notation.

The script reuses the frozen RLB-2B EB and old-Timoshenko workflows. Only the
new RLB material is changed: four equal weakly orthotropic plies form the
symmetric stack [0/90/90/0]. The old-Timoshenko and RLB determinant/SVD
searches remain in the frozen Omega coordinate,

    Omega = omega*l**2*sqrt(rho0*A/(E0*I_y)),

while the frozen EB sign scan uses its historical Lambda=sqrt(Omega)
wavenumber directly. Tables and the plot use Lambda. Curves join pointwise
sorted spectral positions; no modal continuation is performed.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from numpy.typing import NDArray


ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = ROOT / "src"
for import_root in (ROOT, SRC_ROOT):
    if str(import_root) not in sys.path:
        sys.path.insert(0, str(import_root))

from scripts.analysis.laminated_beams import (  # noqa: E402
    compare_rectangular_isotropic_models_vs_beta as rlb2b,
)
from scripts.lib import reddy_symmetric_coupled_beams as rlb_coupled  # noqa: E402
from scripts.lib import reddy_symmetric_laminated_beam as rlb_beam  # noqa: E402


FloatArray = NDArray[np.float64]

STAGE_ID = "RLB-2C"
ALGORITHM_VERSION = "rectangular_weak_orthotropic_old_lambda_beta_grid_v1"
MODEL_EB = rlb2b.MODEL_EB
MODEL_OLD = rlb2b.MODEL_OLD
MODEL_RLB = "NEW_WEAKLY_ORTHOTROPIC_RLB"
MODELS = (MODEL_EB, MODEL_OLD, MODEL_RLB)
ROOT_FILENAMES = {
    MODEL_EB: "eb_roots.csv",
    MODEL_OLD: "old_timoshenko_roots.csv",
    MODEL_RLB: "new_rlb_roots.csv",
}
MODEL_LABELS = {
    MODEL_EB: "Euler--Bernoulli",
    MODEL_OLD: "old rectangular Timoshenko",
    MODEL_RLB: "new weakly orthotropic RLB",
}
LINE_STYLES = {MODEL_EB: "-", MODEL_OLD: "--", MODEL_RLB: "-."}

DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "rectangular_weakly_orthotropic_beta_sweep_comparison"
)

DELTA = 0.1
E0 = 1.0
NU0 = 0.3
RHO0 = 1.0
SHEAR0 = E0 / (2.0 * (1.0 + NU0))
STACK_DEG = (0.0, 90.0, 90.0, 0.0)
PLOTTED_POSITIONS = rlb2b.PLOTTED_POSITIONS
OUTPUT_GUARD_POSITION = rlb2b.OUTPUT_GUARD_POSITION
CONTRACT_RELATIVE_TOLERANCE = 1.0e-11
SYMMETRY_RELATIVE_TOLERANCE = 1.0e-12
NORMALIZATION_IDENTITY_TOLERANCE = 1.0e-12
INVENTORY_VERIFICATION_TOLERANCE = (
    rlb2b.TIMO_PRIMARY_VERIFICATION_RELATIVE_TOLERANCE
)

TASK_INITIAL_GIT_STATE = {
    "top_level": "D:/PHD/CoupledBeams/CoupledBeams",
    "branch": "main",
    "head": "c4719d8853f6ccc4f71e67b35abfe794212603b5",
    "last_commit": "c4719d8 Version 0.4.3.2",
    "status_short": "",
}

FROZEN_MODEL_PATHS = tuple(
    dict.fromkeys(
        (
            *rlb2b.FROZEN_MODEL_PATHS,
            "scripts/analysis/laminated_beams/"
            "compare_rectangular_isotropic_models_vs_beta.py",
            "tests/test_rectangular_isotropic_models_vs_beta.py",
            "docs/laminated_beams/rectangular_isotropic_models_vs_beta_note.md",
        )
    )
)


@dataclass(frozen=True)
class ModelObjects:
    """Baseline EB/Timoshenko objects and the separate weak RLB laminate."""

    baseline_contract: Mapping[str, Any]
    baseline_objects: Any
    weak_laminate: Any
    weak_properties: Any


def _relative(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / max(
        abs(float(left)), abs(float(right)), np.finfo(float).tiny
    )


def _matrix_relative(left: Any, right: Any) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    denominator = max(
        np.linalg.norm(a, ord="fro"),
        np.linalg.norm(b, ord="fro"),
        np.finfo(float).tiny,
    )
    return float(np.linalg.norm(a - b, ord="fro") / denominator)


def _frozen_hashes() -> dict[str, str]:
    return {
        relative: rlb2b._sha256(ROOT / relative)
        for relative in FROZEN_MODEL_PATHS
    }


def beta_grid(
    beta_min_deg: float = rlb2b.DEFAULT_BETA_MIN_DEG,
    beta_max_deg: float = rlb2b.DEFAULT_BETA_MAX_DEG,
    beta_step_deg: float = rlb2b.DEFAULT_BETA_STEP_DEG,
) -> np.ndarray:
    return rlb2b.beta_grid(beta_min_deg, beta_max_deg, beta_step_deg)


def omega_scale(
    *,
    E_reference: float,
    rho_reference: float,
    width: float,
    thickness: float,
    arm_length: float,
) -> float:
    """Return Omega/omega for the fixed isotropic reference normalization."""

    return rlb2b.lambda_scale(
        E=E_reference,
        rho=rho_reference,
        width=width,
        thickness=thickness,
        arm_length=arm_length,
    )


def omega_to_Omega(omega: float, scale: float) -> float:
    return float(omega) * float(scale)


def Omega_to_Lambda(Omega: float) -> float:
    value = float(Omega)
    if value < 0.0 or not math.isfinite(value):
        raise ValueError("Omega must be finite and nonnegative.")
    return math.sqrt(value)


def old_Lambda(omega: float, scale: float) -> float:
    return Omega_to_Lambda(omega_to_Omega(omega, scale))


def build_case_contract(
    betas_deg: Sequence[float] | None = None,
) -> dict[str, Any]:
    active_betas = (
        beta_grid() if betas_deg is None else np.asarray(betas_deg, dtype=float)
    )
    baseline = rlb2b.build_case_contract(active_betas)
    geometry = baseline["geometry"]
    scale = omega_scale(
        E_reference=E0,
        rho_reference=RHO0,
        width=float(geometry["width_b"]),
        thickness=float(geometry["thickness_h"]),
        arm_length=float(geometry["l"]),
    )
    return {
        "schema_version": 1,
        "stage": STAGE_ID,
        "source_contract": baseline["source_contract"],
        "source_contract_sha256": baseline["source_contract_sha256"],
        "mu": 0.0,
        "tau": 0.0,
        "geometry": dict(geometry),
        "reference_isotropic_material": {
            "E0": E0,
            "nu0": NU0,
            "rho0": RHO0,
            "G0": SHEAR0,
        },
        "old_models": {
            "EB": "isotropic rectangular baseline",
            "old_Timoshenko": "isotropic rectangular baseline",
            "K": float(baseline["shear_correction"]["K"]),
        },
        "new_RLB_lamina": {
            "case_id": "A",
            "delta": DELTA,
            "E1": E0 * (1.0 + DELTA),
            "E2": E0 * (1.0 - DELTA),
            "nu12": NU0,
            "G12": SHEAR0,
            "G13": SHEAR0,
            "G23": SHEAR0,
            "rho": RHO0,
        },
        "new_RLB_laminate": {
            "number_of_plies": 4,
            "stacking_sequence_deg": list(STACK_DEG),
            "equal_ply_thickness": float(geometry["thickness_h"]) / 4.0,
            "one_ply_shortcut": False,
            "pipeline": (
                "ply properties->Q->Qbar->A/B/D->shear/mass"
                "->beam reduction->coupled determinant"
            ),
        },
        "frequency": {
            "Omega_definition": "omega*l^2*sqrt(rho0*A/(E0*I_y))",
            "Lambda_definition": "sqrt(Omega)",
            "Lambda_squared_definition": (
                "omega*l^2*sqrt(rho0*A/(E0*I_y))"
            ),
            "Omega_per_omega": scale,
            "historical_EB_matrix_variable": "Lambda",
            "root_search_coordinate": "Omega",
            "plot_coordinate": "Lambda",
        },
        "beta_grid_deg": [float(value) for value in active_betas],
        "beta_step_deg": (
            float(active_betas[1] - active_betas[0])
            if len(active_betas) > 1
            else None
        ),
        "plotted_sorted_positions": PLOTTED_POSITIONS,
        "output_guard_position": OUTPUT_GUARD_POSITION,
        "modal_descendant_tracking": False,
        "curves_meaning": "pointwise independently sorted spectral positions",
    }


def _baseline_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    baseline = rlb2b.build_case_contract(contract["beta_grid_deg"])
    material = baseline["material"]
    geometry = baseline["geometry"]
    if (
        float(material["E"]) != E0
        or float(material["nu"]) != NU0
        or float(material["rho"]) != RHO0
        or float(geometry["L1"]) != 1.0
        or float(geometry["L2"]) != 1.0
        or float(geometry["L_total"]) != 2.0
        or float(geometry["width_b"]) != 0.20
        or float(geometry["thickness_h"]) != 0.05
        or float(baseline["shear_correction"]["K"]) != 5.0 / 6.0
        or any(
            float(geometry[name]) != float(contract["geometry"][name])
            for name in (
                "L1",
                "L2",
                "l",
                "L_total",
                "width_b",
                "thickness_h",
            )
        )
    ):
        raise RuntimeError("The RLB-2B baseline differs from the RLB-2C contract.")
    return baseline


def build_model_objects(contract: Mapping[str, Any]) -> ModelObjects:
    baseline_contract = _baseline_contract(contract)
    baseline_objects = rlb2b.build_model_objects(baseline_contract)
    material_data = contract["new_RLB_lamina"]
    geometry = contract["geometry"]
    material = rlb_beam.LaminaMaterial(
        E1=float(material_data["E1"]),
        E2=float(material_data["E2"]),
        nu12=float(material_data["nu12"]),
        G12=float(material_data["G12"]),
        G13=float(material_data["G13"]),
        G23=float(material_data["G23"]),
        rho=float(material_data["rho"]),
        name="RLB-2C weakly orthotropic case A",
    )
    ply_thickness = float(
        contract["new_RLB_laminate"]["equal_ply_thickness"]
    )
    laminate = rlb_beam.integrate_laminate(
        [rlb_beam.Ply(material, angle, ply_thickness) for angle in STACK_DEG]
    )
    properties = rlb_beam.reduce_to_beam_properties(
        laminate,
        width=float(geometry["width_b"]),
        K=float(contract["old_models"]["K"]),
    )
    return ModelObjects(
        baseline_contract=baseline_contract,
        baseline_objects=baseline_objects,
        weak_laminate=laminate,
        weak_properties=properties,
    )


def constitutive_checks(
    contract: Mapping[str, Any], objects: ModelObjects
) -> dict[str, Any]:
    baseline_check = rlb2b.model_contract_checks(
        objects.baseline_contract, objects.baseline_objects
    )
    material = objects.weak_laminate.plies[0].material
    thickness = float(contract["geometry"]["thickness_h"])
    expected_interfaces = np.asarray(
        [
            -0.5 * thickness,
            -0.25 * thickness,
            0.0,
            0.25 * thickness,
            0.5 * thickness,
        ],
        dtype=float,
    )
    symmetry = rlb_beam.check_laminate_symmetry(
        objects.weak_laminate, tolerance=SYMMETRY_RELATIVE_TOLERANCE
    )
    q0 = rlb_beam.transformed_reduced_stiffness(material, 0.0)
    q90 = rlb_beam.transformed_reduced_stiffness(material, 90.0)
    analytic_A = 0.5 * thickness * (q0 + q90)
    analytic_D = thickness**3 / 96.0 * (7.0 * q0 + q90)
    A_formula_residual = _matrix_relative(objects.weak_laminate.A, analytic_A)
    D_formula_residual = _matrix_relative(objects.weak_laminate.D, analytic_D)
    reference = objects.baseline_objects.rlb_properties
    reference_laminate = objects.baseline_objects.laminate
    weak = objects.weak_properties
    mass_residual = max(
        _relative(weak.m, reference.m), _relative(weak.J, reference.J)
    )
    shear_residual = _relative(weak.S, reference.S)
    shear_matrix_residual = _matrix_relative(
        objects.weak_laminate.shear, reference_laminate.shear
    )
    mass_moment_residual = max(
        _relative(objects.weak_laminate.I0, reference_laminate.I0),
        _relative(objects.weak_laminate.I2, reference_laminate.I2),
    )
    reduction_route_residual = max(
        weak.axial_reduction.relative_difference,
        weak.bending_reduction.relative_difference,
        weak.shear_reduction_before_K.relative_difference,
    )
    exact_material = bool(
        material.E1 == 1.1
        and material.E2 == 0.9
        and material.nu12 == 0.3
        and material.G12 == SHEAR0
        and material.G13 == SHEAR0
        and material.G23 == SHEAR0
        and material.rho == 1.0
    )
    exact_stack = bool(
        len(objects.weak_laminate.plies) == 4
        and tuple(ply.angle_deg for ply in objects.weak_laminate.plies)
        == STACK_DEG
        and all(
            ply.thickness == thickness / 4.0
            for ply in objects.weak_laminate.plies
        )
        and np.array_equal(
            objects.weak_laminate.z_interfaces, expected_interfaces
        )
    )
    stiffness_changed = bool(
        _relative(weak.A, reference.A) > 1.0e-6
        and _relative(weak.D, reference.D) > 1.0e-6
    )
    passed = bool(
        baseline_check["status"] == "PASS"
        and exact_material
        and exact_stack
        and symmetry.is_symmetric
        and symmetry.B_relative <= SYMMETRY_RELATIVE_TOLERANCE
        and symmetry.I1_relative <= SYMMETRY_RELATIVE_TOLERANCE
        and mass_residual <= CONTRACT_RELATIVE_TOLERANCE
        and shear_residual <= CONTRACT_RELATIVE_TOLERANCE
        and shear_matrix_residual <= CONTRACT_RELATIVE_TOLERANCE
        and mass_moment_residual <= CONTRACT_RELATIVE_TOLERANCE
        and A_formula_residual <= 1.0e-12
        and D_formula_residual <= 1.0e-12
        and reduction_route_residual <= CONTRACT_RELATIVE_TOLERANCE
        and stiffness_changed
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "baseline_contract_status": baseline_check["status"],
        "exact_case_A_material": exact_material,
        "exact_four_equal_plies": exact_stack,
        "stacking_sequence_deg": [
            ply.angle_deg for ply in objects.weak_laminate.plies
        ],
        "z_interfaces": objects.weak_laminate.z_interfaces,
        "B_relative": symmetry.B_relative,
        "I1_relative": symmetry.I1_relative,
        "analytic_A_relative_residual": A_formula_residual,
        "analytic_D_relative_residual": D_formula_residual,
        "mass_max_relative_residual": mass_residual,
        "shear_relative_residual": shear_residual,
        "shear_matrix_relative_residual": shear_matrix_residual,
        "mass_moment_max_relative_residual": mass_moment_residual,
        "reduction_route_max_relative": reduction_route_residual,
        "stiffnesses_differ_from_isotropic_reference": stiffness_changed,
        "Abeam_ratio_to_isotropic": weak.A / reference.A,
        "Dbeam_ratio_to_isotropic": weak.D / reference.D,
        "Sbeam_ratio_to_isotropic": weak.S / reference.S,
        "mass_ratio_to_isotropic": weak.m / reference.m,
        "rotary_inertia_ratio_to_isotropic": weak.J / reference.J,
    }


def laminate_property_rows(
    contract: Mapping[str, Any], objects: ModelObjects
) -> list[dict[str, Any]]:
    weak_material = objects.weak_laminate.plies[0].material
    reference_laminate = objects.baseline_objects.laminate
    reference_material = reference_laminate.plies[0].material
    reference_props = objects.baseline_objects.rlb_properties

    def row(
        model: str,
        laminate: Any,
        properties: Any,
        material: Any,
    ) -> dict[str, Any]:
        return {
            "model": model,
            "case_id": (
                "A_DELTA_0.1"
                if model == MODEL_RLB
                else "ISOTROPIC_REFERENCE"
            ),
            "delta": DELTA if model == MODEL_RLB else 0.0,
            "stacking_sequence_deg": [
                ply.angle_deg for ply in laminate.plies
            ],
            "stack_reading_direction": "bottom_to_top_in_z_project",
            "one_ply_shortcut": False,
            "ply_count": len(laminate.plies),
            "equal_ply_thickness": laminate.plies[0].thickness,
            "z_interfaces": laminate.z_interfaces,
            "E1": material.E1,
            "E2": material.E2,
            "nu12": material.nu12,
            "nu21": material.nu21,
            "G12": material.G12,
            "G13": material.G13,
            "G23": material.G23,
            "rho": material.rho,
            "A_matrix": laminate.A,
            "B_matrix": laminate.B,
            "D_matrix": laminate.D,
            "shear_matrix": laminate.shear,
            "I0": laminate.I0,
            "I1": laminate.I1,
            "I2": laminate.I2,
            "A_beam": properties.A,
            "D_beam": properties.D,
            "S_beam": properties.S,
            "mass_per_length": properties.m,
            "rotary_inertia_per_length": properties.J,
            "K": properties.K,
            "B_relative": rlb_beam.check_laminate_symmetry(
                laminate, tolerance=SYMMETRY_RELATIVE_TOLERANCE
            ).B_relative,
            "I1_relative": rlb_beam.check_laminate_symmetry(
                laminate, tolerance=SYMMETRY_RELATIVE_TOLERANCE
            ).I1_relative,
            "A_beam_ratio_to_isotropic": properties.A / reference_props.A,
            "D_beam_ratio_to_isotropic": properties.D / reference_props.D,
            "S_beam_ratio_to_isotropic": properties.S / reference_props.S,
            "mass_ratio_to_isotropic": properties.m / reference_props.m,
            "rotary_inertia_ratio_to_isotropic": (
                properties.J / reference_props.J
            ),
        }

    return [
        row(
            MODEL_RLB,
            objects.weak_laminate,
            objects.weak_properties,
            weak_material,
        ),
        row(
            "ISOTROPIC_REFERENCE",
            reference_laminate,
            reference_props,
            reference_material,
        ),
    ]


def validate_laminate_property_file(
    path: Path,
    contract: Mapping[str, Any],
    objects: ModelObjects,
) -> bool:
    if not path.is_file():
        return False
    actual = rlb2b._read_csv(path)
    expected = laminate_property_rows(contract, objects)
    if len(actual) != len(expected):
        return False
    for actual_row, expected_row in zip(actual, expected, strict=True):
        if set(actual_row) != set(expected_row):
            return False
        for key, value in expected_row.items():
            serialized = str(rlb2b._csv_value(value))
            if actual_row[key] != serialized:
                return False
    return True


def make_matrix_provider(
    model: str,
    beta_deg: float,
    contract: Mapping[str, Any],
    objects: ModelObjects,
    *,
    rlb_arm_cache: dict[float, FloatArray] | None = None,
) -> tuple[Any, dict[str, Any]]:
    if model in (MODEL_EB, MODEL_OLD):
        provider, metadata = rlb2b.make_matrix_provider(
            model,
            beta_deg,
            objects.baseline_contract,
            objects.baseline_objects,
        )
        metadata = dict(metadata)
        if model == MODEL_EB:
            metadata["matrix_argument"] = "Lambda=sqrt(Omega)"
        return provider, metadata
    if model != MODEL_RLB:
        raise ValueError(f"Unknown model: {model!r}.")
    cache = {} if rlb_arm_cache is None else rlb_arm_cache
    beta_rad = math.radians(float(beta_deg))
    length = float(contract["geometry"]["l"])
    properties = objects.weak_properties
    joint = np.asarray(rlb_coupled.joint_matrix(beta_rad), dtype=float)

    def arm_map(omega: float) -> FloatArray:
        key = float(omega)
        if key not in cache:
            cache[key] = np.asarray(
                rlb_coupled.arm_clamp_map(key, length, properties), dtype=float
            )
        return cache[key]

    def provider(omega: float) -> FloatArray:
        combined = np.zeros((12, 6), dtype=float)
        endpoint = arm_map(float(omega))
        combined[:6, :3] = endpoint
        combined[6:, 3:] = endpoint
        return np.asarray(joint @ combined, dtype=float)

    direct_residual = 0.0
    for check_omega in (0.731, 3.217):
        direct = np.asarray(
            rlb_coupled.coupled_boundary_matrix(
                check_omega,
                beta_rad,
                length,
                properties,
                length,
                properties,
            ),
            dtype=float,
        )
        direct_residual = max(
            direct_residual,
            float(np.max(np.abs(provider(check_omega) - direct))),
        )
    if direct_residual > 16.0 * np.finfo(float).eps:
        raise RuntimeError(
            "Cached weak-RLB arm assembly differs from the frozen builder."
        )
    return provider, {
        "matrix_source": "scripts/lib/reddy_symmetric_coupled_beams.py",
        "arm_source": "scripts/lib/reddy_symmetric_laminated_beam.py",
        "case_id": "A",
        "delta": DELTA,
        "number_of_plies": len(objects.weak_laminate.plies),
        "cached_arm_map": True,
        "cached_vs_public_builder_max_abs": direct_residual,
    }


def _optional_sqrt(value: Any) -> Any:
    if value in ("", None):
        return ""
    return Omega_to_Lambda(float(value))


def transform_root_row(
    source: Mapping[str, Any],
    model: str,
    scale: float,
) -> dict[str, Any]:
    raw = dict(source)
    Omega = float(raw.pop("Lambda"))
    omega = float(raw.pop("omega"))
    Lambda = Omega_to_Lambda(Omega)
    bracket_left_Omega = raw.pop("bracket_left_Lambda")
    bracket_right_Omega = raw.pop("bracket_right_Lambda")
    internal_root13_Omega = raw.pop("internal_root13_Lambda")
    identity_residual = max(
        _relative(Omega, omega_to_Omega(omega, scale)),
        _relative(Omega, Lambda**2),
    )
    historical = raw.get("historical_EB_wavenumber", "")
    eb_mapping_residual: Any = ""
    if model == MODEL_EB:
        eb_mapping_residual = _relative(Lambda, float(historical))
    return {
        "model": model,
        "beta_deg": float(raw.pop("beta_deg")),
        "sorted_position": int(raw.pop("sorted_position")),
        "role": raw.pop("role"),
        "omega": omega,
        "Omega": Omega,
        "Lambda": Lambda,
        "normalization_identity_relative_residual": identity_residual,
        "historical_EB_mapping_relative_residual": eb_mapping_residual,
        **raw,
        "bracket_left_Omega": bracket_left_Omega,
        "bracket_right_Omega": bracket_right_Omega,
        "bracket_left_Lambda": _optional_sqrt(bracket_left_Omega),
        "bracket_right_Lambda": _optional_sqrt(bracket_right_Omega),
        "internal_root13_Omega": internal_root13_Omega,
        "internal_root13_Lambda": _optional_sqrt(internal_root13_Omega),
    }


def solve_model_grid(
    model: str,
    betas_deg: Sequence[float],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, Any]]:
    if model not in MODELS:
        raise ValueError(f"Unknown model: {model!r}.")
    contract = build_case_contract(betas_deg)
    objects = build_model_objects(contract)
    checks = constitutive_checks(contract, objects)
    if checks["status"] != "PASS":
        raise RuntimeError(
            "Material/geometry/normalization contract failed before roots."
        )
    policy = rlb2b.frozen_root_policy()
    scale = float(contract["frequency"]["Omega_per_omega"])
    prepared_contract_path = Path(output_dir) / "case_contract.json"
    if not prepared_contract_path.is_file():
        raise RuntimeError("Run --mode prepare before a model worker.")
    prepared_contract = json.loads(
        prepared_contract_path.read_text(encoding="utf-8")
    )
    if prepared_contract != rlb2b._json_value(contract):
        raise RuntimeError("Worker contract differs from the prepared contract.")
    prepared_contract_sha256 = rlb2b._sha256(prepared_contract_path)
    rows: list[dict[str, Any]] = []
    cache: dict[float, FloatArray] = {}
    target = Path(output_dir) / ROOT_FILENAMES[model]
    for beta_index, beta_deg in enumerate(betas_deg, start=1):
        if model == MODEL_EB:
            raw_rows = rlb2b._eb_root_rows(
                float(beta_deg),
                objects.baseline_contract,
                objects.baseline_objects,
                policy,
            )
        else:
            provider, _metadata = make_matrix_provider(
                model,
                float(beta_deg),
                contract,
                objects,
                rlb_arm_cache=cache,
            )
            inventory = rlb2b.iso_inventory.seed_free_root_inventory(
                provider,
                scale,
                policy,
                case_id=f"{model.lower()}__beta_{float(beta_deg):06.2f}",
                builder_id=f"RLB2C_{model}",
                contract_sha256=prepared_contract_sha256,
            )
            if len(inventory.slots) < policy.required_slots:
                raise RuntimeError(
                    f"{model}, beta={float(beta_deg):g}: only "
                    f"{len(inventory.slots)} internal slots were accepted."
                )
            raw_rows = [
                rlb2b._root_row(
                    model, float(beta_deg), inventory, slot
                )
                for slot in inventory.slots[:OUTPUT_GUARD_POSITION]
            ]
        beta_rows = [
            transform_root_row(row, model, scale) for row in raw_rows
        ]
        rows.extend(beta_rows)
        rlb2b._write_csv(target, rows)
        print(
            f"[{model}] beta={float(beta_deg):g} deg "
            f"({beta_index}/{len(betas_deg)}): "
            f"{beta_rows[0]['inventory_status']}, "
            f"Lambda_9={float(beta_rows[-1]['Lambda']):.12g}",
            flush=True,
        )
    return rows


def prepare_output(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    betas_deg: Sequence[float] | None = None,
) -> dict[str, Any]:
    active_betas = (
        beta_grid() if betas_deg is None else np.asarray(betas_deg, dtype=float)
    )
    target = Path(output_dir)
    contract = build_case_contract(active_betas)
    objects = build_model_objects(contract)
    checks = constitutive_checks(contract, objects)
    contract_path = target / "case_contract.json"
    rlb2b._write_json(contract_path, contract)
    rlb2b._write_csv(
        target / "laminate_properties.csv",
        laminate_property_rows(contract, objects),
    )
    manifest = {
        "schema_version": 1,
        "algorithm_version": ALGORITHM_VERSION,
        "stage": STAGE_ID,
        "scientific_role": "finite controlled applicability comparison",
        "task_initial_git_state": TASK_INITIAL_GIT_STATE,
        "git_at_preparation": rlb2b._git_state(),
        "case_contract_sha256": rlb2b._sha256(contract_path),
        "constitutive_checks": checks,
        "common_export_coordinate": "Omega",
        "root_search_coordinates_by_model": {
            MODEL_EB: "historical Lambda=sqrt(Omega)",
            MODEL_OLD: "Omega",
            MODEL_RLB: "Omega",
        },
        "plotted_frequency_coordinate": "Lambda=sqrt(Omega)",
        "models": {
            MODEL_EB: {
                "physics": "frozen isotropic rectangular EB baseline",
                "matrix": "src/my_project/analytic/formulas.py",
            },
            MODEL_OLD: {
                "physics": "frozen isotropic rectangular Timoshenko baseline",
                "matrix": (
                    "scripts/lib/isotropic_rectangular_"
                    "timoshenko_coupled_beams.py"
                ),
            },
            MODEL_RLB: {
                "physics": "frozen laminated RLB with case-A ply inputs",
                "matrix": "scripts/lib/reddy_symmetric_coupled_beams.py",
                "delta": DELTA,
                "plies": 4,
                "one_ply_shortcut": False,
            },
        },
        "thresholds": {
            "shared_contract_relative": CONTRACT_RELATIVE_TOLERANCE,
            "symmetry_relative": SYMMETRY_RELATIVE_TOLERANCE,
            "normalization_identity_relative": (
                NORMALIZATION_IDENTITY_TOLERANCE
            ),
            "root_singular_ratio": rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE,
            "boundary_null_residual": (
                rlb2b.BOUNDARY_RESIDUAL_TOLERANCE
            ),
            "EB_primary_verification_relative": (
                rlb2b.EB_PRIMARY_VERIFICATION_RELATIVE_TOLERANCE
            ),
            "Timoshenko_RLB_primary_verification_relative": (
                INVENTORY_VERIFICATION_TOLERANCE
            ),
        },
        "frozen_root_policy": asdict(rlb2b.frozen_root_policy()),
        "frozen_model_hashes_before": _frozen_hashes(),
        "root_frequencies_cross_seeded_between_models": False,
        "comparison_differences_are_diagnostic_not_gates": True,
        "new_production_solver": False,
        "explicit_exclusions": [
            "delta_values_other_than_0.1",
            "equivalent_isotropic_fitting",
            "new_stacking_sequences",
            "nonzero_mu_or_tau",
            "unequal_lengths",
            "Rayleigh_Ritz",
            "FEM",
            "MAC",
            "modal_branch_tracking",
            "veering",
            "torsion",
            "damping",
            "article_preparation",
            "commit",
            "push",
        ],
    }
    rlb2b._write_json(target / "model_manifest.json", manifest)
    return manifest


def _float(row: Mapping[str, str], key: str) -> float:
    return float(row[key])


def validate_root_rows(
    model: str,
    rows: Sequence[Mapping[str, str]],
    betas_deg: Sequence[float],
) -> dict[str, Any]:
    expected_rows = len(betas_deg) * OUTPUT_GUARD_POSITION
    expected_keys = {
        (round(float(beta), 12), position)
        for beta in betas_deg
        for position in range(1, OUTPUT_GUARD_POSITION + 1)
    }
    actual_keys = {
        (round(_float(row, "beta_deg"), 12), int(row["sorted_position"]))
        for row in rows
    }
    finite = all(
        math.isfinite(_float(row, field))
        for row in rows
        for field in (
            "omega",
            "Omega",
            "Lambda",
            "normalization_identity_relative_residual",
            "scaled_sigma_ratio",
            "boundary_null_residual",
            "primary_verification_max_relative",
        )
    )
    positive = all(
        _float(row, field) > 0.0
        for row in rows
        for field in ("omega", "Omega", "Lambda")
    )
    sorted_positive_positions = True
    for beta in betas_deg:
        beta_rows = sorted(
            (
                row
                for row in rows
                if round(_float(row, "beta_deg"), 12)
                == round(float(beta), 12)
            ),
            key=lambda row: int(row["sorted_position"]),
        )
        sorted_positive_positions = sorted_positive_positions and all(
            _float(right, "Omega") >= _float(left, "Omega")
            for left, right in zip(beta_rows, beta_rows[1:])
        )
    verification_tolerance = (
        rlb2b.EB_PRIMARY_VERIFICATION_RELATIVE_TOLERANCE
        if model == MODEL_EB
        else INVENTORY_VERIFICATION_TOLERANCE
    )
    accepted_root_quality = all(
        row["root_status"] == "PASS"
        and _float(row, "scaled_sigma_ratio")
        <= rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE
        and _float(row, "boundary_null_residual")
        <= rlb2b.BOUNDARY_RESIDUAL_TOLERANCE
        and _float(row, "normalization_identity_relative_residual")
        <= NORMALIZATION_IDENTITY_TOLERANCE
        for row in rows
    )
    if model == MODEL_EB:
        accepted_root_quality = accepted_root_quality and all(
            _float(row, "historical_EB_mapping_relative_residual")
            <= NORMALIZATION_IDENTITY_TOLERANCE
            for row in rows
        )
    accepted_inventory = (
        {"PASS_SIGN_SCAN_CROSSCHECK"} if model == MODEL_EB else {"PASS"}
    )
    inventory_status = all(
        row["inventory_status"] in accepted_inventory for row in rows
    )
    primary_verification = all(
        _float(row, "primary_verification_max_relative")
        <= verification_tolerance
        for row in rows
    )
    no_unresolved = (
        True
        if model == MODEL_EB
        else all(
            int(row["unresolved_candidates_below_internal_guard"]) == 0
            for row in rows
        )
    )
    guard_structure = all(
        (row["guard_flag"].lower() == "true")
        == (int(row["sorted_position"]) == OUTPUT_GUARD_POSITION)
        and row["role"]
        == (
            "ROOT_9_GUARD"
            if int(row["sorted_position"]) == OUTPUT_GUARD_POSITION
            else "FIRST_8"
        )
        for row in rows
    )
    hashes_by_beta = {
        round(float(beta), 12): {
            row["inventory_sha256"]
            for row in rows
            if round(_float(row, "beta_deg"), 12) == round(float(beta), 12)
        }
        for beta in betas_deg
    }
    inventory_hash_structure = bool(
        all(len(values) == 1 for values in hashes_by_beta.values())
        and len({next(iter(values)) for values in hashes_by_beta.values()})
        == len(betas_deg)
    )
    exact_structure = bool(
        len(rows) == expected_rows
        and actual_keys == expected_keys
        and finite
        and positive
        and sorted_positive_positions
        and guard_structure
        and inventory_hash_structure
    )
    verification_inventory = bool(
        inventory_status and primary_verification and no_unresolved
    )
    verification_failure_betas = sorted(
        {
            _float(row, "beta_deg")
            for row in rows
            if row["inventory_status"] not in accepted_inventory
            or _float(row, "primary_verification_max_relative")
            > verification_tolerance
            or (
                model != MODEL_EB
                and int(row["unresolved_candidates_below_internal_guard"])
                != 0
            )
        }
    )
    maximum_verification_row = max(
        rows,
        key=lambda row: _float(row, "primary_verification_max_relative"),
        default=None,
    )
    passed = bool(
        exact_structure
        and accepted_root_quality
        and verification_inventory
    )
    partial = bool(
        not passed and exact_structure and accepted_root_quality
    )
    return {
        "summary_kind": "MODEL_ROOT_INVENTORY",
        "model": model,
        "beta_count": len(betas_deg),
        "expected_rows": expected_rows,
        "actual_rows": len(rows),
        "first_frequency_count": len(betas_deg) * PLOTTED_POSITIONS,
        "guard_count": sum(
            row["guard_flag"].lower() == "true" for row in rows
        ),
        "maximum_scaled_sigma_ratio": max(
            (_float(row, "scaled_sigma_ratio") for row in rows),
            default=math.inf,
        ),
        "maximum_boundary_null_residual": max(
            (_float(row, "boundary_null_residual") for row in rows),
            default=math.inf,
        ),
        "maximum_primary_verification_relative": max(
            (
                _float(row, "primary_verification_max_relative")
                for row in rows
            ),
            default=math.inf,
        ),
        "beta_at_maximum_primary_verification_relative": (
            _float(maximum_verification_row, "beta_deg")
            if maximum_verification_row is not None
            else math.nan
        ),
        "verification_failure_betas_deg": verification_failure_betas,
        "maximum_normalization_identity_relative_residual": max(
            (
                _float(row, "normalization_identity_relative_residual")
                for row in rows
            ),
            default=math.inf,
        ),
        "inventory_hash_count": len(
            {row["inventory_sha256"] for row in rows}
        ),
        "unresolved_candidate_sum": (
            "NOT_ASSESSED_BY_EB_SIGN_SCAN"
            if model == MODEL_EB
            else sum(
                int(row["unresolved_candidates_below_internal_guard"])
                for row in rows
            )
        ),
        "guard_structure_passed": guard_structure,
        "positive_frequencies_passed": positive,
        "sorted_positions_nondecreasing_passed": sorted_positive_positions,
        "inventory_hash_structure_passed": inventory_hash_structure,
        "exact_row_structure_passed": exact_structure,
        "accepted_root_quality_passed": accepted_root_quality,
        "primary_verification_passed": primary_verification,
        "inventory_status_passed": inventory_status,
        "verification_inventory_passed": verification_inventory,
        "status": "PASS" if passed else "PARTIAL_PASS" if partial else "FAIL",
    }


def comparison_rows(
    reference_rows: Sequence[Mapping[str, str]],
    rlb_rows: Sequence[Mapping[str, str]],
    *,
    reference_tag: str,
) -> list[dict[str, Any]]:
    reference = {
        (round(_float(row, "beta_deg"), 12), int(row["sorted_position"])): row
        for row in reference_rows
        if int(row["sorted_position"]) <= PLOTTED_POSITIONS
    }
    new = {
        (round(_float(row, "beta_deg"), 12), int(row["sorted_position"])): row
        for row in rlb_rows
        if int(row["sorted_position"]) <= PLOTTED_POSITIONS
    }
    if set(reference) != set(new):
        raise RuntimeError(f"{reference_tag} and weak RLB position keys differ.")
    rows: list[dict[str, Any]] = []
    for beta_deg, position in sorted(reference):
        left = reference[(beta_deg, position)]
        right = new[(beta_deg, position)]
        Lambda_left = _float(left, "Lambda")
        Lambda_right = _float(right, "Lambda")
        Omega_left = _float(left, "Omega")
        Omega_right = _float(right, "Omega")
        omega_left = _float(left, "omega")
        omega_right = _float(right, "omega")

        def differences(a: float, b: float) -> tuple[float, float]:
            absolute = abs(a - b)
            relative = absolute / max(
                abs(a), abs(b), np.finfo(float).tiny
            )
            return absolute, relative

        Lambda_abs, Lambda_rel = differences(Lambda_left, Lambda_right)
        Omega_abs, Omega_rel = differences(Omega_left, Omega_right)
        omega_abs, omega_rel = differences(omega_left, omega_right)
        rows.append(
            {
                "comparison": f"{reference_tag}_vs_NEW_WEAK_RLB",
                "beta_deg": beta_deg,
                "sorted_position": position,
                f"Lambda_{reference_tag}": Lambda_left,
                "Lambda_new_rlb": Lambda_right,
                "abs_difference_Lambda": Lambda_abs,
                "rel_difference_Lambda": Lambda_rel,
                f"Omega_{reference_tag}": Omega_left,
                "Omega_new_rlb": Omega_right,
                "abs_difference_Omega": Omega_abs,
                "rel_difference_Omega": Omega_rel,
                f"omega_{reference_tag}": omega_left,
                "omega_new_rlb": omega_right,
                "abs_difference_omega": omega_abs,
                "rel_difference_omega": omega_rel,
                "comparison_semantics": "INDEPENDENTLY_SORTED_POSITIONS",
                "scientific_role": (
                    "DIAGNOSTIC_DIFFERENCE_NOT_AGREEMENT_GATE"
                ),
            }
        )
    return rows


def comparison_summary(
    rows: Sequence[Mapping[str, Any]], model: str
) -> dict[str, Any]:
    max_Lambda = max(
        rows, key=lambda row: float(row["rel_difference_Lambda"])
    )
    max_Omega = max(
        rows, key=lambda row: float(row["rel_difference_Omega"])
    )
    max_omega = max(
        rows, key=lambda row: float(row["rel_difference_omega"])
    )
    return {
        "summary_kind": "SORTED_POSITION_COMPARISON",
        "model": model,
        "comparison_count": len(rows),
        "maximum_relative_difference_Lambda": float(
            max_Lambda["rel_difference_Lambda"]
        ),
        "beta_at_max_relative_Lambda": float(max_Lambda["beta_deg"]),
        "position_at_max_relative_Lambda": int(
            max_Lambda["sorted_position"]
        ),
        "maximum_relative_difference_Omega": float(
            max_Omega["rel_difference_Omega"]
        ),
        "beta_at_max_relative_Omega": float(max_Omega["beta_deg"]),
        "position_at_max_relative_Omega": int(
            max_Omega["sorted_position"]
        ),
        "maximum_relative_difference_omega": float(
            max_omega["rel_difference_omega"]
        ),
        "beta_at_max_relative_omega": float(max_omega["beta_deg"]),
        "position_at_max_relative_omega": int(
            max_omega["sorted_position"]
        ),
        "status": "DIAGNOSTIC_ONLY",
    }


def create_plot(
    rows_by_model: Mapping[str, Sequence[Mapping[str, Any]]],
    output_path: Path,
) -> Path:
    fig, ax = plt.subplots(figsize=(11.0, 6.6), constrained_layout=True)
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, PLOTTED_POSITIONS))
    for position in range(1, PLOTTED_POSITIONS + 1):
        for model in MODELS:
            selected = sorted(
                (
                    row
                    for row in rows_by_model[model]
                    if int(row["sorted_position"]) == position
                ),
                key=lambda row: float(row["beta_deg"]),
            )
            ax.plot(
                [float(row["beta_deg"]) for row in selected],
                [float(row["Lambda"]) for row in selected],
                color=colors[position - 1],
                linestyle=LINE_STYLES[model],
                linewidth=1.55 if model != MODEL_RLB else 1.35,
                alpha=0.96,
            )
    ax.set_xlabel(r"$\beta$, degrees")
    ax.set_ylabel(r"$\Lambda$")
    ax.set_xlim(0.0, 90.0)
    ax.grid(True, alpha=0.24)
    model_legend = ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color="black",
                linestyle=LINE_STYLES[model],
                linewidth=1.7,
                label=MODEL_LABELS[model].replace("--", "–"),
            )
            for model in MODELS
        ],
        loc="upper left",
        frameon=False,
        fontsize=9,
    )
    ax.add_artist(model_legend)
    ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color=colors[position - 1],
                linestyle="-",
                linewidth=1.8,
                label=f"sorted position {position}",
            )
            for position in range(1, PLOTTED_POSITIONS + 1)
        ],
        loc="upper right",
        ncols=2,
        frameon=False,
        fontsize=8,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=190)
    plt.close(fig)
    return output_path


def _report_text(
    contract: Mapping[str, Any],
    checks: Mapping[str, Any],
    root_summaries: Sequence[Mapping[str, Any]],
    old_summary: Mapping[str, Any],
    eb_summary: Mapping[str, Any],
    statuses: Mapping[str, str],
) -> str:
    material = contract["new_RLB_lamina"]
    geometry = contract["geometry"]
    lines = [
        "# Слабо ортотропный RLB и старые прямоугольные модели",
        "",
        "## Область расчёта",
        "",
        (
            "Рассматривается конечная сетка $\\beta=0,1,\\ldots,90^\\circ$ "
            "при $\\mu=\\tau=0$. Модели Эйлера--Бернулли и Тимошенко "
            "сохраняют изотропный baseline. Новая RLB использует weakly "
            "orthotropic case A."
        ),
        "",
        (
            "Линии соединяют независимо отсортированные спектральные позиции. "
            "Они не задают продолжение модальных форм."
        ),
        "",
        "## Нормировка",
        "",
        "$$",
        "\\Omega=\\omega l^2\\sqrt{\\frac{\\rho_0 A}{E_0 I_y}},",
        "\\qquad \\Lambda=\\sqrt{\\Omega},",
        "\\qquad A=bh,\\quad I_y=\\frac{bh^3}{12}.",
        "$$",
        "",
        (
            f"Использованы $E_0=1$, $\\nu_0=0.3$, $\\rho_0=1$, "
            f"$b={float(geometry['width_b']):g}$, "
            f"$h={float(geometry['thickness_h']):g}$ и "
            f"$l={float(geometry['l']):g}$. Для old Timoshenko и RLB "
            "determinant/SVD search выполнялся в $\\Omega$. Замороженный "
            "EB sign scan выполнялся непосредственно в историческом "
            "$\\Lambda=\\sqrt{\\Omega}$."
        ),
        "",
        "## Слабо ортотропный ламинат",
        "",
        (
            f"Стопка $[0/90/90/0]$ содержит четыре слоя толщиной $h/4$. "
            f"Приняты $\\delta={float(material['delta']):g}$, "
            f"$E_1={float(material['E1']):g}$, "
            f"$E_2={float(material['E2']):g}$, "
            f"$\\nu_{{12}}={float(material['nu12']):g}$, "
            "$G_{12}=G_{13}=G_{23}=1/2.6$, $\\rho=1$ и $K=5/6$. "
            "One-ply shortcut не использовался: четыре слоя прошли цепочку "
            "$Q\\rightarrow\\overline Q\\rightarrow A/B/D$, shear/mass "
            "integration и beam reduction."
        ),
        "",
        (
            f"Scaled $B$-residual равен "
            f"{float(checks['B_relative']):.6e}. "
            "Отношения редуцированных жёсткостей к изотропному baseline "
            f"равны {float(checks['Abeam_ratio_to_isotropic']):.9g} "
            "для растяжения и "
            f"{float(checks['Dbeam_ratio_to_isotropic']):.9g} "
            "для изгиба. Сдвиговая жёсткость и массовые характеристики "
            "не изменились."
        ),
        "",
        "## Root inventory",
        "",
        (
            "| Модель | строк | guards | max singular ratio | "
            "max boundary residual | status |"
        ),
        "|---|---:|---:|---:|---:|---|",
    ]
    for summary in root_summaries:
        lines.append(
            f"| {summary['model']} | {int(summary['actual_rows'])} | "
            f"{int(summary['guard_count'])} | "
            f"{float(summary['maximum_scaled_sigma_ratio']):.6e} | "
            f"{float(summary['maximum_boundary_null_residual']):.6e} | "
            f"{summary['status']} |"
        )
    for summary in root_summaries:
        if summary["status"] == "PARTIAL_PASS":
            betas = ", ".join(
                f"{float(value):g}"
                for value in summary["verification_failure_betas_deg"]
            )
            lines.extend(
                [
                    "",
                    (
                        f"Для {summary['model']} экспортированные корни "
                        "прошли индивидуальные singular/boundary gates, но "
                        "primary/verification inventory не совпали в "
                        f"пределах frozen threshold при $\\beta={betas}^\\circ$. "
                        "Поэтому root status ограничен PARTIAL_PASS; "
                        "расхождение не скрывается и policy не меняется."
                    ),
                ]
            )
    lines.extend(
        [
            "",
            "## Диагностическое сопоставление",
            "",
            (
                "Для old Timoshenko и weak RLB максимальная относительная "
                f"разность $\\Lambda$ равна "
                f"{float(old_summary['maximum_relative_difference_Lambda']):.6e} "
                f"при $\\beta={float(old_summary['beta_at_max_relative_Lambda']):g}^\\circ$, "
                f"позиция {int(old_summary['position_at_max_relative_Lambda'])}. "
                "Максимальная относительная разность $\\Omega$ равна "
                f"{float(old_summary['maximum_relative_difference_Omega']):.6e} "
                f"при $\\beta={float(old_summary['beta_at_max_relative_Omega']):g}^\\circ$, "
                f"позиция {int(old_summary['position_at_max_relative_Omega'])}."
            ),
            "",
            (
                "Для EB и weak RLB максимальная относительная разность "
                f"$\\Lambda$ равна "
                f"{float(eb_summary['maximum_relative_difference_Lambda']):.6e} "
                f"при $\\beta={float(eb_summary['beta_at_max_relative_Lambda']):g}^\\circ$, "
                f"позиция {int(eb_summary['position_at_max_relative_Lambda'])}. "
                "Максимальная относительная разность $\\Omega$ равна "
                f"{float(eb_summary['maximum_relative_difference_Omega']):.6e} "
                f"при $\\beta={float(eb_summary['beta_at_max_relative_Omega']):g}^\\circ$, "
                f"позиция {int(eb_summary['position_at_max_relative_Omega'])}."
            ),
            "",
            (
                "Относительная разность определена как "
                "$|a-b|/\\max(|a|,|b|)$. Эти величины описывают различие "
                "независимо отсортированных позиций. Они не используются "
                "как agreement gate и не устанавливают идентичность "
                "модальных потомков."
            ),
            "",
            "## График",
            "",
            (
                "Файл [lambda_vs_beta_plot.png](lambda_vs_beta_plot.png) "
                "содержит 24 кривые. Цвет задаёт отсортированную позицию "
                "$k=1,\\ldots,8$; EB показана сплошной линией, old "
                "Timoshenko -- пунктирной, weak RLB -- штрихпунктирной."
            ),
            "",
            "## Статусы",
            "",
        ]
    )
    lines.extend(f"- {name}: {value}" for name, value in statuses.items())
    lines.extend(
        [
            "",
            "## Ограничения",
            "",
            (
                "Результат относится только к case A при $\\delta=0.1$, "
                "одной геометрии и принятой угловой сетке. Не выполнялись "
                "equivalent isotropic fitting, MAC, branch tracking, Ritz "
                "и FEM. Не изучались veering, кручение и демпфирование; "
                "статья не подготавливалась."
            ),
        ]
    )
    return "\n".join(lines) + "\n"


def finalize_output(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    betas_deg: Sequence[float] | None = None,
) -> dict[str, Any]:
    target = Path(output_dir)
    active_betas = (
        beta_grid() if betas_deg is None else np.asarray(betas_deg, dtype=float)
    )
    contract_path = target / "case_contract.json"
    manifest_path = target / "model_manifest.json"
    if not contract_path.is_file() or not manifest_path.is_file():
        raise RuntimeError("Run --mode prepare before finalization.")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    model_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    model_manifest.pop("root_search_coordinate", None)
    model_manifest["common_export_coordinate"] = "Omega"
    model_manifest["root_search_coordinates_by_model"] = {
        MODEL_EB: "historical Lambda=sqrt(Omega)",
        MODEL_OLD: "Omega",
        MODEL_RLB: "Omega",
    }
    model_manifest["case_contract_coordinate_note"] = (
        "frequency.root_search_coordinate=Omega names the common exported "
        "inventory coordinate; the frozen EB sign scan itself uses the "
        "historical Lambda=sqrt(Omega) wavenumber"
    )
    rlb2b._write_json(manifest_path, model_manifest)
    if not np.array_equal(
        active_betas, np.asarray(contract["beta_grid_deg"], dtype=float)
    ):
        raise RuntimeError(
            "Finalization grid differs from the prepared contract."
        )
    contract_hash = rlb2b._sha256(contract_path)
    if model_manifest["case_contract_sha256"] != contract_hash:
        raise RuntimeError("Prepared case-contract hash no longer matches.")
    preservation_after = _frozen_hashes()
    if model_manifest["frozen_model_hashes_before"] != preservation_after:
        raise RuntimeError("A frozen RLB-2B or production file changed.")

    rows_by_model = {
        model: rlb2b._read_csv(target / ROOT_FILENAMES[model])
        for model in MODELS
    }
    root_summaries = [
        validate_root_rows(model, rows_by_model[model], active_betas)
        for model in MODELS
    ]
    old_comparison = comparison_rows(
        rows_by_model[MODEL_OLD],
        rows_by_model[MODEL_RLB],
        reference_tag="old_timo",
    )
    eb_comparison = comparison_rows(
        rows_by_model[MODEL_EB],
        rows_by_model[MODEL_RLB],
        reference_tag="EB",
    )
    rlb2b._write_csv(
        target / "old_vs_new_rlb_comparison.csv", old_comparison
    )
    rlb2b._write_csv(
        target / "eb_vs_new_rlb_comparison.csv", eb_comparison
    )
    old_summary = comparison_summary(
        old_comparison, "OLD_TIMOSHENKO_vs_NEW_WEAK_RLB"
    )
    eb_summary = comparison_summary(
        eb_comparison, "EB_vs_NEW_WEAK_RLB"
    )
    rlb2b._write_csv(
        target / "spectrum_summary.csv",
        [*root_summaries, old_summary, eb_summary],
    )
    required_betas = beta_grid()
    exact_full_grid_pass = bool(
        np.array_equal(active_betas, required_betas)
        and np.array_equal(
            np.asarray(contract["beta_grid_deg"], dtype=float),
            required_betas,
        )
    )
    comparison_counts_pass = bool(
        len(old_comparison) == len(required_betas) * PLOTTED_POSITIONS
        and len(eb_comparison) == len(required_betas) * PLOTTED_POSITIONS
    )
    plot_path = create_plot(
        rows_by_model, target / "lambda_vs_beta_plot.png"
    )

    objects = build_model_objects(contract)
    checks = constitutive_checks(contract, objects)
    laminate_properties_pass = validate_laminate_property_file(
        target / "laminate_properties.csv", contract, objects
    )
    scale = float(contract["frequency"]["Omega_per_omega"])
    expected_scale = omega_scale(
        E_reference=E0,
        rho_reference=RHO0,
        width=float(contract["geometry"]["width_b"]),
        thickness=float(contract["geometry"]["thickness_h"]),
        arm_length=float(contract["geometry"]["l"]),
    )
    normalization_pass = bool(
        math.isfinite(scale)
        and scale > 0.0
        and _relative(scale, expected_scale)
        <= NORMALIZATION_IDENTITY_TOLERANCE
        and _relative(
            old_Lambda(0.731, scale) ** 2,
            omega_to_Omega(0.731, scale),
        )
        <= NORMALIZATION_IDENTITY_TOLERANCE
    )
    root_summary_statuses = [summary["status"] for summary in root_summaries]
    if (
        exact_full_grid_pass
        and comparison_counts_pass
        and all(value == "PASS" for value in root_summary_statuses)
    ):
        root_status = "PASS"
    elif (
        exact_full_grid_pass
        and comparison_counts_pass
        and all(
            value in {"PASS", "PARTIAL_PASS"}
            for value in root_summary_statuses
        )
    ):
        root_status = "PARTIAL_PASS"
    else:
        root_status = "FAIL"
    weak_summary = next(
        summary
        for summary in root_summaries
        if summary["model"] == MODEL_RLB
    )
    if (
        checks["status"] == "PASS"
        and laminate_properties_pass
        and weak_summary["status"] == "PASS"
    ):
        weak_status = "PASS"
    elif (
        checks["status"] == "PASS"
        and laminate_properties_pass
        and weak_summary["status"] == "PARTIAL_PASS"
    ):
        weak_status = "PARTIAL_PASS"
    else:
        weak_status = "FAIL"
    plot_pass = plot_path.is_file() and plot_path.stat().st_size > 0
    statuses = {
        "RLB-2C-LAMBDA-DEFINITION": (
            "PASS" if normalization_pass else "FAIL"
        ),
        "RLB-2C-ROOT-DATA": root_status,
        "RLB-2C-WEAK-ORTHOTROPIC-RLB-RUN": weak_status,
        "RLB-2C-PLOT-GENERATION": "PASS" if plot_pass else "FAIL",
    }
    statuses["OVERALL"] = (
        "PASS"
        if all(value == "PASS" for value in statuses.values())
        else "FAIL"
    )
    (target / "report.md").write_text(
        _report_text(
            contract,
            checks,
            root_summaries,
            old_summary,
            eb_summary,
            statuses,
        ),
        encoding="utf-8",
    )
    generated_names = [
        "model_manifest.json",
        "case_contract.json",
        "laminate_properties.csv",
        "eb_roots.csv",
        "old_timoshenko_roots.csv",
        "new_rlb_roots.csv",
        "old_vs_new_rlb_comparison.csv",
        "eb_vs_new_rlb_comparison.csv",
        "spectrum_summary.csv",
        "lambda_vs_beta_plot.png",
        "report.md",
    ]
    run_manifest = {
        "schema_version": 1,
        "algorithm_version": ALGORITHM_VERSION,
        "stage": STAGE_ID,
        "task_initial_git_state": TASK_INITIAL_GIT_STATE,
        "finalization_git_state": rlb2b._git_state(),
        "case_contract_sha256": contract_hash,
        "model_manifest_sha256": rlb2b._sha256(manifest_path),
        "frequency_coordinates": {
            "common_export": "Omega",
            "root_search_by_model": {
                MODEL_EB: "historical Lambda=sqrt(Omega)",
                MODEL_OLD: "Omega",
                MODEL_RLB: "Omega",
            },
            "plot_and_primary_CSV": "Lambda=sqrt(Omega)",
            "Omega_per_omega": scale,
            "case_contract_coordinate_note": (
                "frequency.root_search_coordinate=Omega names the common "
                "exported inventory coordinate, not the raw EB scan variable"
            ),
        },
        "beta_grid_deg": [float(value) for value in active_betas],
        "models": list(MODELS),
        "plotted_sorted_positions": PLOTTED_POSITIONS,
        "output_guard_position": OUTPUT_GUARD_POSITION,
        "root_frequencies_cross_seeded_between_models": False,
        "comparison_differences_used_as_agreement_gate": False,
        "summary": {
            "constitutive": checks,
            "laminate_properties_semantic_check_passed": (
                laminate_properties_pass
            ),
            "root_data": root_summaries,
            "exact_full_grid_passed": exact_full_grid_pass,
            "comparison_counts_passed": comparison_counts_pass,
            "old_vs_new": old_summary,
            "eb_vs_new": eb_summary,
        },
        "statuses": statuses,
        "frozen_model_hashes_before": (
            model_manifest["frozen_model_hashes_before"]
        ),
        "frozen_model_hashes_after": preservation_after,
        "frozen_models_preserved": (
            model_manifest["frozen_model_hashes_before"]
            == preservation_after
        ),
        "generated_file_hashes": {
            name: rlb2b._sha256(target / name)
            for name in generated_names
        },
        "figures_created": 1,
        "new_production_solver_created": False,
        "delta_values_run": [DELTA],
        "case_ids_run": ["A"],
        "mu": 0.0,
        "tau": 0.0,
        "Ritz_run": False,
        "FEM_run": False,
        "branch_tracking_run": False,
        "commit_performed": False,
        "push_performed": False,
    }
    rlb2b._write_json(target / "run_manifest.json", run_manifest)
    return run_manifest


def run_all(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    betas_deg: Sequence[float] | None = None,
) -> dict[str, Any]:
    active_betas = (
        beta_grid() if betas_deg is None else np.asarray(betas_deg, dtype=float)
    )
    prepare_output(output_dir, active_betas)
    for model in MODELS:
        solve_model_grid(model, active_betas, output_dir)
    return finalize_output(output_dir, active_betas)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=(
            "prepare",
            "eb",
            "old_timoshenko",
            "new_rlb",
            "finalize",
            "all",
        ),
        default="all",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR
    )
    parser.add_argument("--beta-min-deg", type=float, default=0.0)
    parser.add_argument("--beta-max-deg", type=float, default=90.0)
    parser.add_argument("--beta-step-deg", type=float, default=1.0)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    betas = beta_grid(
        args.beta_min_deg, args.beta_max_deg, args.beta_step_deg
    )
    output_dir = Path(args.output_dir)
    if args.mode == "prepare":
        prepare_output(output_dir, betas)
    elif args.mode == "eb":
        solve_model_grid(MODEL_EB, betas, output_dir)
    elif args.mode == "old_timoshenko":
        solve_model_grid(MODEL_OLD, betas, output_dir)
    elif args.mode == "new_rlb":
        solve_model_grid(MODEL_RLB, betas, output_dir)
    elif args.mode == "finalize":
        finalize_output(output_dir, betas)
    else:
        run_all(output_dir, betas)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
