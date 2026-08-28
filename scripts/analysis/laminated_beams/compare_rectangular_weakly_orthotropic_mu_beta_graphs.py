"""RLB-2D plots of pointwise sorted frequencies over ``mu`` and ``beta``.

The stage preserves the RLB-2C material and frequency contracts.  Euler--
Bernoulli (EB) and the old rectangular Timoshenko model remain isotropic,
whereas the new RLB arm is the four-ply weakly orthotropic case A.  Every
accepted dimensional angular frequency ``omega`` is exported together with

    Omega = omega*l**2*sqrt(rho0*A0/(E0*I_y0)),
    Lambda = sqrt(Omega).

Only independently sorted positions 1--8 are plotted.  Position 9 is retained
as an output guard.  No modal continuation or inter-model difference is
computed.

The frozen EB sign scan accepts unequal lengths but only one section
slenderness.  It is therefore used unchanged for the ``tau=0`` mu sweep.  The
``tau=0.2`` beta sweep uses an additive physical-state EB matrix adapter for
the two actual rectangular sections.  Its determinant is passed to the frozen
EB sign-scan algorithm; only the determinant provider changes.  At ``tau=0``
an exact determinant identity checks this adapter against the frozen
closed-form EB matrix.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Mapping, MutableMapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from numpy.typing import NDArray
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = ROOT / "src"
for import_root in (ROOT, SRC_ROOT):
    if str(import_root) not in sys.path:
        sys.path.insert(0, str(import_root))

from scripts.analysis.laminated_beams import (  # noqa: E402
    compare_rectangular_weakly_orthotropic_models_vs_beta as rlb2c,
)


FloatArray = NDArray[np.float64]

STAGE_ID = "RLB-2D"
ALGORITHM_VERSION = "rectangular_weak_orthotropic_mu_beta_graphs_v1"

MODEL_EB = rlb2c.MODEL_EB
MODEL_OLD = rlb2c.MODEL_OLD
MODEL_RLB = rlb2c.MODEL_RLB
MODELS = (MODEL_EB, MODEL_OLD, MODEL_RLB)
MODEL_LABELS = dict(rlb2c.MODEL_LABELS)
LINE_STYLES = dict(rlb2c.LINE_STYLES)

SWEEP_MU = "mu"
SWEEP_BETA = "beta"
SWEEPS = (SWEEP_MU, SWEEP_BETA)

ROOT_FILENAMES = {
    (SWEEP_MU, MODEL_EB): "mu_sweep_eb_roots.csv",
    (SWEEP_MU, MODEL_OLD): "mu_sweep_old_timoshenko_roots.csv",
    (SWEEP_MU, MODEL_RLB): "mu_sweep_new_rlb_roots.csv",
    (SWEEP_BETA, MODEL_EB): "beta_sweep_eb_roots.csv",
    (SWEEP_BETA, MODEL_OLD): "beta_sweep_old_timoshenko_roots.csv",
    (SWEEP_BETA, MODEL_RLB): "beta_sweep_new_rlb_roots.csv",
}
PLOT_FILENAMES = {
    SWEEP_MU: "lambda_vs_mu_plot.png",
    SWEEP_BETA: "lambda_vs_beta_plot.png",
}
CLOSING_CHECKPOINT_FILENAME = "closing_checkpoint.json"
CLOSING_THREAD_LIMITS = {
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}

DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "rectangular_weakly_orthotropic_graphs_mu_beta"
)

E0 = rlb2c.E0
NU0 = rlb2c.NU0
RHO0 = rlb2c.RHO0
G0 = rlb2c.SHEAR0
DELTA = rlb2c.DELTA
STACK_DEG = rlb2c.STACK_DEG
B0 = 0.20
H0 = 0.05
L_REFERENCE = 1.0
L_TOTAL = 2.0
K = 5.0 / 6.0
A0 = B0 * H0
I_Y0 = B0 * H0**3 / 12.0

MU_BETA_DEG = 15.0
MU_TAU = 0.0
DEFAULT_MU_MIN = 0.0
DEFAULT_MU_MAX = 0.8
REQUESTED_MU_STEP = 0.01
ALLOWED_MU_FALLBACK_STEP = 0.02
DEFAULT_MU_STEP = ALLOWED_MU_FALLBACK_STEP
MU_FALLBACK_REASON = "EXTREME_MU_FULL_FROZEN_INVENTORY_OVER_5_MINUTES"

BETA_MU = 0.5
BETA_TAU = 0.2
DEFAULT_BETA_MIN_DEG = 0.0
DEFAULT_BETA_MAX_DEG = 90.0
DEFAULT_BETA_STEP_DEG = 1.0

PLOTTED_POSITIONS = 8
OUTPUT_GUARD_POSITION = 9
NORMALIZATION_IDENTITY_TOLERANCE = 1.0e-12
CONTRACT_RELATIVE_TOLERANCE = 1.0e-11
EB_DETERMINANT_IDENTITY_TOLERANCE = 1.0e-10
INVENTORY_VERIFICATION_TOLERANCE = (
    rlb2c.INVENTORY_VERIFICATION_TOLERANCE
)

RLB2C_SCRIPT = (
    "scripts/analysis/laminated_beams/"
    "compare_rectangular_weakly_orthotropic_models_vs_beta.py"
)
RLB2C_TEST = "tests/test_rectangular_weakly_orthotropic_models_vs_beta.py"
RLB2C_NOTE = (
    "docs/laminated_beams/"
    "rectangular_weakly_orthotropic_models_vs_beta_note.md"
)
FROZEN_MODEL_PATHS = tuple(
    dict.fromkeys(
        (*rlb2c.FROZEN_MODEL_PATHS, RLB2C_SCRIPT, RLB2C_TEST, RLB2C_NOTE)
    )
)


@dataclass(frozen=True)
class GeometryPoint:
    """Actual arm geometry at one independently solved sweep point."""

    mu: float
    tau: float
    beta_deg: float
    L1: float
    L2: float
    h1: float
    h2: float
    b1: float
    b2: float


@dataclass(frozen=True)
class ArmObjects:
    """Section and four-ply objects for one arm."""

    length: float
    width: float
    thickness: float
    legacy_section: Any
    weak_laminate: Any
    weak_properties: Any


@dataclass(frozen=True)
class ModelObjects:
    """Arm-specific objects shared by the three matrix providers."""

    arm1: ArmObjects
    arm2: ArmObjects


@dataclass
class ArmPairCache:
    """Frequency caches whose identity is restricted to one arm pair."""

    arm1: MutableMapping[float, FloatArray]
    arm2: MutableMapping[float, FloatArray]

    @classmethod
    def empty(cls) -> "ArmPairCache":
        return cls(arm1={}, arm2={})


def _relative(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / max(
        abs(float(left)), abs(float(right)), np.finfo(float).tiny
    )


def _frozen_hashes() -> dict[str, str]:
    return {
        relative: rlb2c.rlb2b._sha256(ROOT / relative)
        for relative in FROZEN_MODEL_PATHS
    }


def _inclusive_grid(lower: float, upper: float, step: float) -> np.ndarray:
    lower_value = float(lower)
    upper_value = float(upper)
    step_value = float(step)
    if not all(
        math.isfinite(value)
        for value in (lower_value, upper_value, step_value)
    ):
        raise ValueError("Grid inputs must be finite.")
    if step_value <= 0.0 or upper_value < lower_value:
        raise ValueError("Require step>0 and upper>=lower.")
    count = int(
        math.floor((upper_value - lower_value) / step_value + 1.0e-12)
    ) + 1
    values = lower_value + step_value * np.arange(count, dtype=float)
    if values[-1] < upper_value - 1.0e-12:
        values = np.append(values, upper_value)
    return np.round(values, 12)


def mu_grid(
    mu_min: float = DEFAULT_MU_MIN,
    mu_max: float = DEFAULT_MU_MAX,
    mu_step: float = DEFAULT_MU_STEP,
) -> np.ndarray:
    """Return the endpoint-inclusive RLB-2D length-asymmetry grid."""

    return _inclusive_grid(mu_min, mu_max, mu_step)


def beta_grid(
    beta_min_deg: float = DEFAULT_BETA_MIN_DEG,
    beta_max_deg: float = DEFAULT_BETA_MAX_DEG,
    beta_step_deg: float = DEFAULT_BETA_STEP_DEG,
) -> np.ndarray:
    """Return the endpoint-inclusive RLB-2D angle grid in degrees."""

    return _inclusive_grid(beta_min_deg, beta_max_deg, beta_step_deg)


def geometry_for(
    mu: float, tau: float, beta_deg: float = 0.0
) -> GeometryPoint:
    """Apply the rectangular ``mu``/``tau`` parameterization exactly."""

    mu_value = float(mu)
    tau_value = float(tau)
    beta_value = float(beta_deg)
    if not all(
        math.isfinite(value) for value in (mu_value, tau_value, beta_value)
    ):
        raise ValueError("mu, tau, and beta must be finite.")
    if not (0.0 <= mu_value < 1.0):
        raise ValueError("mu must satisfy 0<=mu<1 for positive arm lengths.")
    if not (-1.0 < tau_value < 1.0):
        raise ValueError("tau must satisfy -1<tau<1 for positive thicknesses.")
    # Decimal sweep inputs are part of the exported contract; remove only
    # binary representation noise from their elementary affine formulas.
    length_1 = round(L_REFERENCE * (1.0 - mu_value), 15)
    length_2 = round(L_REFERENCE * (1.0 + mu_value), 15)
    thickness_1 = round(H0 * (1.0 - tau_value), 15)
    thickness_2 = round(H0 * (1.0 + tau_value), 15)
    return GeometryPoint(
        mu=mu_value,
        tau=tau_value,
        beta_deg=beta_value,
        L1=length_1,
        L2=length_2,
        h1=thickness_1,
        h2=thickness_2,
        b1=B0,
        b2=B0,
    )


def reference_omega_scale() -> float:
    """Return the fixed reference value ``Omega/omega``."""

    return rlb2c.omega_scale(
        E_reference=E0,
        rho_reference=RHO0,
        width=B0,
        thickness=H0,
        arm_length=L_REFERENCE,
    )


def omega_to_Omega(omega: float, scale: float | None = None) -> float:
    return rlb2c.omega_to_Omega(
        float(omega), reference_omega_scale() if scale is None else float(scale)
    )


def Omega_to_Lambda(Omega: float) -> float:
    return rlb2c.Omega_to_Lambda(float(Omega))


def old_Lambda(omega: float, scale: float | None = None) -> float:
    return rlb2c.old_Lambda(
        float(omega), reference_omega_scale() if scale is None else float(scale)
    )


def build_case_contract(
    mu_values: Sequence[float] | None = None,
    beta_values: Sequence[float] | None = None,
) -> dict[str, Any]:
    """Build one immutable contract for both required finite sweeps."""

    active_mu = (
        mu_grid() if mu_values is None else np.asarray(mu_values, dtype=float)
    )
    active_beta = (
        beta_grid()
        if beta_values is None
        else np.asarray(beta_values, dtype=float)
    )
    if active_mu.ndim != 1 or active_mu.size == 0:
        raise ValueError("mu_values must be a nonempty one-dimensional grid.")
    if active_beta.ndim != 1 or active_beta.size == 0:
        raise ValueError("beta_values must be a nonempty one-dimensional grid.")
    for value in active_mu:
        geometry_for(float(value), MU_TAU, MU_BETA_DEG)
    for value in active_beta:
        geometry_for(BETA_MU, BETA_TAU, float(value))

    source = rlb2c.build_case_contract([0.0])
    scale = reference_omega_scale()
    mu_step = (
        float(active_mu[1] - active_mu[0]) if active_mu.size > 1 else None
    )
    beta_step = (
        float(active_beta[1] - active_beta[0])
        if active_beta.size > 1
        else None
    )
    fallback_used = bool(
        mu_step is not None
        and abs(mu_step - REQUESTED_MU_STEP) > 1.0e-14
    )
    return {
        "schema_version": 1,
        "stage": STAGE_ID,
        "source_stage": rlb2c.STAGE_ID,
        "source_contract": source["source_contract"],
        "source_contract_sha256": source["source_contract_sha256"],
        "reference_constants": {
            "E0": E0,
            "nu0": NU0,
            "rho0": RHO0,
            "b0": B0,
            "h0": H0,
            "l": L_REFERENCE,
            "L_total": L_TOTAL,
            "K": K,
            "A0": A0,
            "I_y0": I_Y0,
        },
        "geometry_parameterization": {
            "L1": "l*(1-mu)",
            "L2": "l*(1+mu)",
            "h1": "h0*(1-tau)",
            "h2": "h0*(1+tau)",
            "b1": "b0",
            "b2": "b0",
        },
        "old_models": {
            "EB": "isotropic rectangular baseline with actual arm sections",
            "old_Timoshenko": (
                "isotropic rectangular baseline with actual arm sections"
            ),
            "equivalent_isotropic_fitting": False,
        },
        "new_RLB_lamina": {
            "case_id": "A",
            "delta": DELTA,
            "E1": E0 * (1.0 + DELTA),
            "E2": E0 * (1.0 - DELTA),
            "nu12": NU0,
            "G12": G0,
            "G13": G0,
            "G23": G0,
            "rho": RHO0,
        },
        "new_RLB_laminate": {
            "number_of_plies_per_arm": 4,
            "stacking_sequence_deg": list(STACK_DEG),
            "ply_thickness_arm_i": "h_i/4",
            "one_ply_shortcut": False,
            "pipeline": (
                "ply properties->Q->Qbar->A/B/D->shear/mass"
                "->beam reduction->coupled determinant"
            ),
        },
        "normalization": {
            "Omega_definition": (
                "omega*l^2*sqrt(rho0*A0/(E0*I_y0))"
            ),
            "Lambda_definition": "sqrt(Omega)",
            "Lambda_squared_definition": (
                "omega*l^2*sqrt(rho0*A0/(E0*I_y0))"
            ),
            "Omega_per_omega": scale,
            "reference_only": True,
        },
        "sweeps": {
            SWEEP_MU: {
                "axis": "mu",
                "values": [float(value) for value in active_mu],
                "step": mu_step,
                "beta_deg": MU_BETA_DEG,
                "tau": MU_TAU,
                "requested_step": REQUESTED_MU_STEP,
                "allowed_runtime_fallback_step": ALLOWED_MU_FALLBACK_STEP,
                "fallback_used": fallback_used,
                "fallback_reason": (
                    MU_FALLBACK_REASON if fallback_used else None
                ),
                "fallback_decision_used_spectrum": False,
            },
            SWEEP_BETA: {
                "axis": "beta_deg",
                "values": [float(value) for value in active_beta],
                "step": beta_step,
                "mu": BETA_MU,
                "tau": BETA_TAU,
            },
        },
        "plotted_sorted_positions": PLOTTED_POSITIONS,
        "output_guard_position": OUTPUT_GUARD_POSITION,
        "modal_descendant_tracking": False,
        "curves_meaning": (
            "pointwise independently sorted spectral positions"
        ),
        "inter_model_relative_differences_computed": False,
    }


def _weak_material() -> Any:
    data = {
        "E1": E0 * (1.0 + DELTA),
        "E2": E0 * (1.0 - DELTA),
        "nu12": NU0,
        "G12": G0,
        "G13": G0,
        "G23": G0,
        "rho": RHO0,
    }
    return rlb2c.rlb_beam.LaminaMaterial(
        **data, name="RLB-2D weakly orthotropic case A"
    )


def _build_arm(length: float, width: float, thickness: float) -> ArmObjects:
    legacy_section = rlb2c.rlb2b.legacy_timoshenko.rectangular_section(
        E=E0,
        nu=NU0,
        rho=RHO0,
        width=float(width),
        thickness=float(thickness),
        K=K,
    )
    material = _weak_material()
    weak_laminate = rlb2c.rlb_beam.integrate_laminate(
        [
            rlb2c.rlb_beam.Ply(
                material, float(angle), float(thickness) / 4.0
            )
            for angle in STACK_DEG
        ]
    )
    weak_properties = rlb2c.rlb_beam.reduce_to_beam_properties(
        weak_laminate, width=float(width), K=K
    )
    return ArmObjects(
        length=float(length),
        width=float(width),
        thickness=float(thickness),
        legacy_section=legacy_section,
        weak_laminate=weak_laminate,
        weak_properties=weak_properties,
    )


def build_model_objects(geometry: GeometryPoint) -> ModelObjects:
    """Build both actual arm sections without an equivalent-fit shortcut."""

    return ModelObjects(
        arm1=_build_arm(geometry.L1, geometry.b1, geometry.h1),
        arm2=_build_arm(geometry.L2, geometry.b2, geometry.h2),
    )


def constitutive_checks(
    geometry: GeometryPoint, objects: ModelObjects
) -> dict[str, Any]:
    """Check arm geometry, four-ply construction, and reduced properties."""

    checks: list[dict[str, Any]] = []
    exact_material = True
    for arm_index, (arm, expected_length, expected_width, expected_h) in enumerate(
        (
            (objects.arm1, geometry.L1, geometry.b1, geometry.h1),
            (objects.arm2, geometry.L2, geometry.b2, geometry.h2),
        ),
        start=1,
    ):
        material = arm.weak_laminate.plies[0].material
        symmetry = rlb2c.rlb_beam.check_laminate_symmetry(
            arm.weak_laminate,
            tolerance=rlb2c.SYMMETRY_RELATIVE_TOLERANCE,
        )
        exact_material = exact_material and bool(
            material.E1 == 1.1
            and material.E2 == 0.9
            and material.nu12 == 0.3
            and material.G12 == G0
            and material.G13 == G0
            and material.G23 == G0
            and material.rho == 1.0
        )
        exact_stack = bool(
            len(arm.weak_laminate.plies) == 4
            and tuple(
                float(ply.angle_deg) for ply in arm.weak_laminate.plies
            )
            == STACK_DEG
            and all(
                ply.thickness == float(expected_h) / 4.0
                for ply in arm.weak_laminate.plies
            )
        )
        legacy = arm.legacy_section
        expected_area = float(expected_width) * float(expected_h)
        expected_inertia = (
            float(expected_width) * float(expected_h) ** 3 / 12.0
        )
        residual = max(
            _relative(arm.length, expected_length),
            _relative(arm.width, expected_width),
            _relative(arm.thickness, expected_h),
            _relative(legacy.area, expected_area),
            _relative(legacy.second_moment, expected_inertia),
            _relative(legacy.EA, E0 * expected_area),
            _relative(legacy.EI, E0 * expected_inertia),
            _relative(legacy.rhoA, RHO0 * expected_area),
        )
        arm_pass = bool(
            exact_stack
            and symmetry.is_symmetric
            and residual <= CONTRACT_RELATIVE_TOLERANCE
        )
        checks.append(
            {
                "arm": arm_index,
                "status": "PASS" if arm_pass else "FAIL",
                "exact_four_equal_plies": exact_stack,
                "B_relative": symmetry.B_relative,
                "I1_relative": symmetry.I1_relative,
                "geometry_and_isotropic_section_max_relative": residual,
            }
        )
    passed = exact_material and all(item["status"] == "PASS" for item in checks)
    return {
        "status": "PASS" if passed else "FAIL",
        "exact_case_A_material": exact_material,
        "arms": checks,
    }


def geometry_sanity_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for case_id, mu, tau in (
        ("mu_0_tau_0", 0.0, 0.0),
        ("mu_0p8_tau_0", 0.8, 0.0),
        ("mu_0p5_tau_0p2", 0.5, 0.2),
    ):
        geometry = geometry_for(mu, tau)
        rows.append(
            {
                "case_id": case_id,
                "mu": geometry.mu,
                "tau": geometry.tau,
                "L1": geometry.L1,
                "L2": geometry.L2,
                "h1": geometry.h1,
                "h2": geometry.h2,
                "b1": geometry.b1,
                "b2": geometry.b2,
            }
        )
    return rows


def laminate_property_rows() -> list[dict[str, Any]]:
    """Return weak-RLB arm summaries for the three sanity geometries."""

    rows: list[dict[str, Any]] = []
    for geometry_row in geometry_sanity_rows():
        geometry = geometry_for(
            float(geometry_row["mu"]), float(geometry_row["tau"])
        )
        objects = build_model_objects(geometry)
        for arm_index, arm in enumerate((objects.arm1, objects.arm2), start=1):
            props = arm.weak_properties
            rows.append(
                {
                    "case_id": geometry_row["case_id"],
                    "arm": arm_index,
                    "mu": geometry.mu,
                    "tau": geometry.tau,
                    "length": arm.length,
                    "width": arm.width,
                    "thickness": arm.thickness,
                    "ply_count": len(arm.weak_laminate.plies),
                    "stacking_sequence_deg": [
                        ply.angle_deg for ply in arm.weak_laminate.plies
                    ],
                    "ply_thicknesses": [
                        ply.thickness for ply in arm.weak_laminate.plies
                    ],
                    "A_beam": props.A,
                    "D_beam": props.D,
                    "S_beam": props.S,
                    "m": props.m,
                    "J": props.J,
                    "B_relative": rlb2c.rlb_beam.check_laminate_symmetry(
                        arm.weak_laminate
                    ).B_relative,
                    "one_ply_shortcut": False,
                }
            )
    return rows


def _eb_state_matrix(omega: float, section: Any) -> FloatArray:
    """Return physical EB state equations in canonical RLB joint ordering."""

    frequency = float(omega)
    if not math.isfinite(frequency) or frequency < 0.0:
        raise ValueError("omega must be finite and nonnegative.")
    omega2 = frequency * frequency
    return np.asarray(
        [
            [0.0, 0.0, 0.0, 1.0 / section.EA, 0.0, 0.0],
            [0.0, 0.0, -1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / section.EI],
            [-section.rhoA * omega2, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, -section.rhoA * omega2, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        ],
        dtype=float,
    )


def _cached_eb_arm_map(
    omega: float,
    arm: ArmObjects,
    cache: MutableMapping[float, FloatArray],
) -> FloatArray:
    key = float(omega)
    if key not in cache:
        cache[key] = np.asarray(
            expm(_eb_state_matrix(key, arm.legacy_section) * arm.length)
            @ rlb2c.rlb_coupled.CLAMP_BASIS,
            dtype=float,
        )
    return cache[key]


def _physical_eb_provider(
    geometry: GeometryPoint,
    objects: ModelObjects,
    cache: ArmPairCache,
) -> Any:
    joint = np.asarray(
        rlb2c.rlb_coupled.joint_matrix(math.radians(geometry.beta_deg)),
        dtype=float,
    )

    def provider(omega: float) -> FloatArray:
        combined = np.zeros((12, 6), dtype=float)
        combined[:6, :3] = _cached_eb_arm_map(
            float(omega), objects.arm1, cache.arm1
        )
        combined[6:, 3:] = _cached_eb_arm_map(
            float(omega), objects.arm2, cache.arm2
        )
        return np.asarray(joint @ combined, dtype=float)

    return provider


def make_matrix_provider(
    model: str,
    geometry: GeometryPoint,
    objects: ModelObjects,
    *,
    cache: ArmPairCache | None = None,
    force_physical_eb: bool = False,
) -> tuple[Any, dict[str, Any]]:
    """Build an arm-specific matrix provider without changing production APIs."""

    if model not in MODELS:
        raise ValueError(f"Unknown model: {model!r}.")
    active_cache = ArmPairCache.empty() if cache is None else cache
    beta_rad = math.radians(geometry.beta_deg)
    if model == MODEL_EB:
        if geometry.tau == 0.0 and not force_physical_eb:
            epsilon = H0 / (math.sqrt(12.0) * L_REFERENCE)
            scale = reference_omega_scale()

            def provider(omega: float) -> FloatArray:
                Omega = omega_to_Omega(float(omega), scale)
                historical_wavenumber = Omega_to_Lambda(Omega)
                return np.asarray(
                    rlb2c.rlb2b.assemble_eb_coupled_matrix(
                        historical_wavenumber,
                        beta_rad,
                        geometry.mu,
                        epsilon,
                    ),
                    dtype=float,
                )

            return provider, {
                "matrix_source": "src/my_project/analytic/formulas.py",
                "matrix_argument": "historical Lambda=sqrt(Omega)",
                "root_workflow": "frozen EB sign scan",
                "mu": geometry.mu,
                "tau": geometry.tau,
                "epsilon": epsilon,
            }
        return _physical_eb_provider(geometry, objects, active_cache), {
            "matrix_source": "RLB-2D additive physical EB state adapter",
            "state_order": ["u", "w", "psi", "N", "Q", "M"],
            "root_workflow": "frozen EB sign scan with physical determinant",
            "reason": "frozen EB API has no unequal-thickness argument",
            "actual_arm_thicknesses": [geometry.h1, geometry.h2],
        }

    if model == MODEL_OLD:

        def provider(omega: float) -> FloatArray:
            return np.asarray(
                rlb2c.rlb2b.legacy_timoshenko.legacy_coupled_boundary_matrix_raw(
                    float(omega),
                    objects.arm1.legacy_section,
                    geometry.L1,
                    objects.arm2.legacy_section,
                    geometry.L2,
                    beta_deg=geometry.beta_deg,
                ),
                dtype=float,
            )

        return provider, {
            "matrix_source": (
                "scripts/lib/isotropic_rectangular_"
                "timoshenko_coupled_beams.py"
            ),
            "actual_arm_sections": True,
        }

    joint = np.asarray(rlb2c.rlb_coupled.joint_matrix(beta_rad), dtype=float)

    def rlb_arm_map(
        omega: float,
        arm: ArmObjects,
        arm_cache: MutableMapping[float, FloatArray],
    ) -> FloatArray:
        key = float(omega)
        if key not in arm_cache:
            arm_cache[key] = np.asarray(
                rlb2c.rlb_coupled.arm_clamp_map(
                    key, arm.length, arm.weak_properties
                ),
                dtype=float,
            )
        return arm_cache[key]

    def provider(omega: float) -> FloatArray:
        combined = np.zeros((12, 6), dtype=float)
        combined[:6, :3] = rlb_arm_map(
            float(omega), objects.arm1, active_cache.arm1
        )
        combined[6:, 3:] = rlb_arm_map(
            float(omega), objects.arm2, active_cache.arm2
        )
        return np.asarray(joint @ combined, dtype=float)

    direct_residual = 0.0
    for check_omega in (0.731, 3.217):
        direct = np.asarray(
            rlb2c.rlb_coupled.coupled_boundary_matrix(
                check_omega,
                beta_rad,
                geometry.L1,
                objects.arm1.weak_properties,
                geometry.L2,
                objects.arm2.weak_properties,
            ),
            dtype=float,
        )
        direct_residual = max(
            direct_residual,
            float(np.max(np.abs(provider(check_omega) - direct))),
        )
    if direct_residual > 16.0 * np.finfo(float).eps:
        raise RuntimeError("Cached arm-pair RLB assembly differs from public builder.")
    return provider, {
        "matrix_source": "scripts/lib/reddy_symmetric_coupled_beams.py",
        "case_id": "A",
        "delta": DELTA,
        "number_of_plies_per_arm": [
            len(objects.arm1.weak_laminate.plies),
            len(objects.arm2.weak_laminate.plies),
        ],
        "two_arm_caches": True,
        "cached_vs_public_builder_max_abs": direct_residual,
    }


def eb_tau0_equivalence_diagnostic() -> dict[str, Any]:
    """Check the exact determinant identity of the additive EB adapter.

    For equal rectangular arm sections the physical-state and frozen EB
    matrices satisfy

    ``det(D_phys) = -l^7 det(D_frozen)/(16 (E I_y)^3 Lambda^7)``.

    The check uses ordinary non-root probes and never supplies search seeds.
    """

    epsilon = H0 / (math.sqrt(12.0) * L_REFERENCE)
    scale = reference_omega_scale()
    rows: list[dict[str, Any]] = []
    maximum = 0.0
    for mu, beta_deg, historical_Lambda in (
        (0.0, 15.0, 2.3),
        (0.4, 37.0, 4.7),
        # Keep the long-arm probe below the raw hyperbolic determinant's
        # cancellation range; accepted-root SVD gates cover larger values.
        (0.8, 15.0, 3.1),
    ):
        geometry = geometry_for(mu, 0.0, beta_deg)
        objects = build_model_objects(geometry)
        physical, _ = make_matrix_provider(
            MODEL_EB,
            geometry,
            objects,
            force_physical_eb=True,
        )
        Omega = historical_Lambda**2
        omega = Omega / scale
        frozen = np.asarray(
            rlb2c.rlb2b.assemble_eb_coupled_matrix(
                historical_Lambda,
                math.radians(beta_deg),
                mu,
                epsilon,
            ),
            dtype=float,
        )
        actual = float(np.linalg.det(physical(omega)))
        expected = float(
            -L_REFERENCE**7
            * np.linalg.det(frozen)
            / (
                16.0
                * (E0 * I_Y0) ** 3
                * historical_Lambda**7
            )
        )
        residual = _relative(actual, expected)
        maximum = max(maximum, residual)
        rows.append(
            {
                "mu": mu,
                "beta_deg": beta_deg,
                "historical_Lambda": historical_Lambda,
                "omega": omega,
                "actual_physical_determinant": actual,
                "identity_determinant": expected,
                "relative_residual": residual,
            }
        )
    return {
        "status": (
            "PASS"
            if maximum <= EB_DETERMINANT_IDENTITY_TOLERANCE
            else "FAIL"
        ),
        "identity": (
            "det(D_phys)=-l^7*det(D_frozen)/"
            "(16*(E0*I_y0)^3*Lambda^7)"
        ),
        "tolerance": EB_DETERMINANT_IDENTITY_TOLERANCE,
        "maximum_relative_residual": maximum,
        "probes": rows,
    }


def _inventory_export_diagnostics(
    inventory: Any, policy: Any
) -> dict[str, Any]:
    """Qualify primary/verification agreement only through root 9."""

    primary = inventory.primary.slots[:OUTPUT_GUARD_POSITION]
    verification = inventory.verification.slots[:OUTPUT_GUARD_POSITION]
    agreement = bool(
        len(primary) == OUTPUT_GUARD_POSITION
        and len(verification) == OUTPUT_GUARD_POSITION
    )
    maximum_relative = 0.0
    relative_by_position: list[float] = []
    agreement_by_position: list[bool] = []
    if agreement:
        for left, right in zip(primary, verification, strict=True):
            left_clustered = bool(left.event.cluster_id)
            right_clustered = bool(right.event.cluster_id)
            left_value = (
                left.event.cluster_center_Omega
                if left_clustered
                else left.event.Omega
            )
            right_value = (
                right.event.cluster_center_Omega
                if right_clustered
                else right.event.Omega
            )
            relative = _relative(left_value, right_value)
            position_agreement = bool(
                relative <= INVENTORY_VERIFICATION_TOLERANCE
                and left.event.multiplicity == right.event.multiplicity
                and left.event.detected_nullity == right.event.detected_nullity
                and left_clustered == right_clustered
                and left.event.cluster_multiplicity
                == right.event.cluster_multiplicity
                and left.event.cluster_total_nullity
                == right.event.cluster_total_nullity
            )
            relative_by_position.append(relative)
            agreement_by_position.append(position_agreement)
            maximum_relative = max(maximum_relative, relative)
            if not position_agreement:
                agreement = False

    export_guard = (
        float(primary[-1].event.Omega) if len(primary) == 9 else -math.inf
    )
    accepted_events = (*inventory.primary.events, *inventory.verification.events)

    def genuinely_unresolved(candidate: Any) -> bool:
        if candidate.interval_left_Omega > export_guard:
            return False
        if candidate.rejection_reason == "DUPLICATE_DETECTION_SAME_ROOT":
            return False
        if any(
            abs(candidate.Omega - event.Omega)
            <= policy.dedup_atol_Omega
            + policy.dedup_rtol
            * max(abs(candidate.Omega), abs(event.Omega))
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

    unresolved_below_export = sum(
        genuinely_unresolved(candidate)
        for candidate in (
            *inventory.primary.rejected_candidates,
            *inventory.verification.rejected_candidates,
        )
    )
    passed = bool(agreement and unresolved_below_export == 0)
    plotted_agreement = bool(
        len(agreement_by_position) == OUTPUT_GUARD_POSITION
        and all(agreement_by_position[:PLOTTED_POSITIONS])
        and unresolved_below_export == 0
    )
    guard_agreement = bool(
        len(agreement_by_position) == OUTPUT_GUARD_POSITION
        and agreement_by_position[OUTPUT_GUARD_POSITION - 1]
        and unresolved_below_export == 0
    )
    return {
        "export_guard_available": len(primary) == OUTPUT_GUARD_POSITION,
        "export_primary_verification_agreement": agreement,
        "export_primary_verification_max_relative": maximum_relative,
        "unresolved_candidates_below_export_guard": unresolved_below_export,
        "export_range_status": "PASS" if passed else "FAIL",
        "plotted_first8_primary_verification_agreement": plotted_agreement,
        "plotted_first8_primary_verification_max_relative": (
            max(relative_by_position[:PLOTTED_POSITIONS], default=math.inf)
        ),
        "root9_primary_verification_agreement": guard_agreement,
        "root9_primary_verification_relative": (
            relative_by_position[OUTPUT_GUARD_POSITION - 1]
            if len(relative_by_position) == OUTPUT_GUARD_POSITION
            else math.inf
        ),
        "primary_verification_relative_by_position": relative_by_position,
        "internal_inventory_status": inventory.status,
        "internal_primary_verification_max_relative": (
            inventory.maximum_primary_verification_relative
        ),
        "internal_unresolved_candidate_count": (
            inventory.unresolved_low_sigma_count
        ),
    }


def transform_root_row(
    source: Mapping[str, Any],
    model: str,
    geometry: GeometryPoint,
    sweep: str,
    export_diagnostics: Mapping[str, Any],
) -> dict[str, Any]:
    """Apply the RLB-2C coordinate transform and add actual geometry."""

    transformed = rlb2c.transform_root_row(
        source, model, reference_omega_scale()
    )
    return {
        "sweep": sweep,
        "model": transformed.pop("model"),
        "mu": geometry.mu,
        "tau": geometry.tau,
        "beta_deg": transformed.pop("beta_deg"),
        "L1": geometry.L1,
        "L2": geometry.L2,
        "h1": geometry.h1,
        "h2": geometry.h2,
        "b1": geometry.b1,
        "b2": geometry.b2,
        **transformed,
        **dict(export_diagnostics),
    }


def _eb_sign_scan_rows(
    geometry: GeometryPoint,
    objects: ModelObjects,
    policy: Any,
    *,
    cache: ArmPairCache | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run the frozen EB sign-scan algorithm for the actual arm sections.

    At ``tau=0`` the frozen closed-form determinant is used directly.  For
    unequal thicknesses the same solver function is temporarily supplied with
    the positively scaled physical-state determinant.  Positive row/column
    factors preserve its zeros and sign; no root-search code is copied or
    modified.
    """

    epsilon = H0 / (math.sqrt(12.0) * L_REFERENCE)
    physical_determinant = geometry.tau != 0.0
    provider, _ = make_matrix_provider(
        MODEL_EB,
        geometry,
        objects,
        cache=cache,
        force_physical_eb=physical_determinant,
    )
    scale = reference_omega_scale()

    def determinant_adapter(
        historical_wavenumber: float,
        _beta_rad: float,
        _mu: float,
        _epsilon: float,
    ) -> float:
        Omega = float(historical_wavenumber) ** 2
        matrix = provider(Omega / scale)
        scaled = rlb2c.rlb_coupled.positively_equilibrate_matrix(
            matrix
        ).scaled_matrix
        return float(np.linalg.det(scaled))

    def scan(step: float) -> np.ndarray:
        return np.asarray(
            rlb2c.rlb2b.find_eb_roots(
                beta=math.radians(geometry.beta_deg),
                mu=geometry.mu,
                eps=epsilon,
                n_roots=OUTPUT_GUARD_POSITION,
                Lmin=0.2,
                Lmax0=35.0,
                scan_step=float(step),
                grow_factor=1.35,
                max_tries=8,
            ),
            dtype=float,
        )

    solver_globals = rlb2c.rlb2b.find_eb_roots.__globals__
    original_determinant = solver_globals["det_clamped_coupled"]
    if physical_determinant:
        solver_globals["det_clamped_coupled"] = determinant_adapter
    try:
        roots = scan(rlb2c.rlb2b.EB_PRIMARY_SCAN_STEP)
        verification = scan(rlb2c.rlb2b.EB_VERIFICATION_SCAN_STEP)
    finally:
        if physical_determinant:
            solver_globals["det_clamped_coupled"] = original_determinant
    for label, values in (("primary", roots), ("verification", verification)):
        if values.shape != (OUTPUT_GUARD_POSITION,) or not np.all(
            np.isfinite(values)
        ):
            raise RuntimeError(
                f"EB {geometry}: the {label} sign scan did not return "
                f"{OUTPUT_GUARD_POSITION} finite roots."
            )
    primary_Omega = roots**2
    verification_Omega = verification**2
    verification_relative = float(
        np.max(
            np.abs(primary_Omega - verification_Omega)
            / np.maximum.reduce(
                (
                    np.abs(primary_Omega),
                    np.abs(verification_Omega),
                    np.full(primary_Omega.shape, np.finfo(float).tiny),
                )
            )
        )
    )
    raw_rows: list[dict[str, Any]] = []
    digest_rows: list[dict[str, Any]] = []
    for position, historical_wavenumber in enumerate(roots, start=1):
        Omega = float(historical_wavenumber**2)
        diagnostic = rlb2c.rlb2b.iso_inventory._boundary_matrix_diagnostics(
            Omega, provider, scale, policy
        )
        quality_pass = bool(
            diagnostic.scaled_sigma_ratio
            <= rlb2c.rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE
            and diagnostic.raw_boundary_null_residual
            <= rlb2c.rlb2b.BOUNDARY_RESIDUAL_TOLERANCE
        )
        row = {
            "model": MODEL_EB,
            "beta_deg": geometry.beta_deg,
            "sorted_position": position,
            "role": (
                "FIRST_8" if position <= PLOTTED_POSITIONS else "ROOT_9_GUARD"
            ),
            "omega": Omega / scale,
            "Lambda": Omega,
            "historical_EB_wavenumber": float(historical_wavenumber),
            "multiplicity": "NOT_ASSESSED_BY_SIGN_SCAN",
            "detected_nullity": diagnostic.detected_nullity,
            "cluster_id": "",
            "cluster_semantics": "NO_CLUSTER_CLAIM_SIGN_SCAN_ONLY",
            "cluster_multiplicity": "NOT_ASSESSED_BY_SIGN_SCAN",
            "cluster_total_nullity": "NOT_ASSESSED_BY_SIGN_SCAN",
            "raw_determinant": diagnostic.raw_determinant,
            "scaled_determinant": diagnostic.scaled_determinant,
            "raw_sigma_min": diagnostic.raw_sigma_min,
            "raw_sigma_max": diagnostic.raw_sigma_max,
            "raw_sigma_ratio": diagnostic.raw_sigma_ratio,
            "scaled_sigma_min": diagnostic.scaled_sigma_min,
            "scaled_sigma_max": diagnostic.scaled_sigma_max,
            "scaled_sigma_ratio": diagnostic.scaled_sigma_ratio,
            "boundary_null_residual": diagnostic.raw_boundary_null_residual,
            "bracket_left_Lambda": "",
            "bracket_right_Lambda": "",
            "detector_type": (
                (
                    "frozen_EB_sign_scan_on_physical_determinant_primary"
                    if physical_determinant
                    else "frozen_EB_sign_scan_bisection_primary"
                ),
                (
                    "finer_EB_sign_scan_on_physical_determinant_verification"
                    if physical_determinant
                    else "finer_sign_scan_bisection_verification"
                ),
            ),
            "root_status": "PASS" if quality_pass else "FAIL",
            "inventory_status": "PENDING",
            "inventory_sha256": "PENDING",
            "primary_slot_count_internal": OUTPUT_GUARD_POSITION,
            "verification_slot_count_internal": OUTPUT_GUARD_POSITION,
            "internal_requested_roots": PLOTTED_POSITIONS,
            "internal_guard_position": OUTPUT_GUARD_POSITION,
            "internal_root13_Lambda": "",
            "primary_verification_max_relative": verification_relative,
            "unresolved_candidates_below_internal_guard": (
                "NOT_ASSESSED_BY_EB_SIGN_SCAN"
            ),
            "guard_flag": position == OUTPUT_GUARD_POSITION,
        }
        raw_rows.append(row)
        digest_rows.append(
            {
                "position": position,
                "Omega": format(Omega, ".17g"),
                "verification_Omega": format(
                    float(verification_Omega[position - 1]), ".17g"
                ),
                "scaled_sigma_ratio": format(
                    diagnostic.scaled_sigma_ratio, ".17g"
                ),
                "boundary_residual": format(
                    diagnostic.raw_boundary_null_residual, ".17g"
                ),
            }
        )
    inventory_status = (
        "PASS_SIGN_SCAN_CROSSCHECK"
        if all(row["root_status"] == "PASS" for row in raw_rows)
        and verification_relative
        <= rlb2c.rlb2b.EB_PRIMARY_VERIFICATION_RELATIVE_TOLERANCE
        else "FAIL"
    )
    digest = hashlib.sha256(
        json.dumps(
            {
                "stage": STAGE_ID,
                "model": MODEL_EB,
                "mu": geometry.mu,
                "tau": geometry.tau,
                "beta_deg": geometry.beta_deg,
                "primary_scan_step": rlb2c.rlb2b.EB_PRIMARY_SCAN_STEP,
                "verification_scan_step": (
                    rlb2c.rlb2b.EB_VERIFICATION_SCAN_STEP
                ),
                "roots": digest_rows,
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest().upper()
    for row in raw_rows:
        row["inventory_status"] = inventory_status
        row["inventory_sha256"] = digest
    export = {
        "export_guard_available": True,
        "export_primary_verification_agreement": (
            verification_relative
            <= rlb2c.rlb2b.EB_PRIMARY_VERIFICATION_RELATIVE_TOLERANCE
        ),
        "export_primary_verification_max_relative": verification_relative,
        "unresolved_candidates_below_export_guard": (
            "NOT_ASSESSED_BY_EB_SIGN_SCAN"
        ),
        "export_range_status": (
            "PASS" if inventory_status == "PASS_SIGN_SCAN_CROSSCHECK" else "FAIL"
        ),
        "internal_inventory_status": inventory_status,
        "internal_primary_verification_max_relative": verification_relative,
        "internal_unresolved_candidate_count": (
            "NOT_ASSESSED_BY_EB_SIGN_SCAN"
        ),
    }
    return raw_rows, export


def _inventory_rows(
    model: str,
    geometry: GeometryPoint,
    objects: ModelObjects,
    policy: Any,
    contract_sha256: str,
    cache: ArmPairCache,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    provider, _ = make_matrix_provider(
        model, geometry, objects, cache=cache
    )
    inventory = rlb2c.rlb2b.iso_inventory.seed_free_root_inventory(
        provider,
        reference_omega_scale(),
        policy,
        case_id=(
            f"rlb2d__{model.lower()}__mu_{geometry.mu:.4f}__"
            f"tau_{geometry.tau:.4f}__beta_{geometry.beta_deg:.4f}"
        ),
        builder_id=f"RLB2D_{model}",
        contract_sha256=contract_sha256,
    )
    if len(inventory.slots) < policy.required_slots:
        raise RuntimeError(
            f"{model}, geometry={geometry}: only {len(inventory.slots)} "
            "internal slots were accepted."
        )
    rows = [
        rlb2c.rlb2b._root_row(
            model, geometry.beta_deg, inventory, slot
        )
        for slot in inventory.slots[:OUTPUT_GUARD_POSITION]
    ]
    return rows, _inventory_export_diagnostics(inventory, policy)


def _worker_contract(
    output_dir: Path,
    mu_values: Sequence[float],
    beta_values: Sequence[float],
) -> tuple[dict[str, Any], str]:
    path = Path(output_dir) / "case_contract.json"
    if not path.is_file():
        raise RuntimeError("Run --mode prepare before a model worker.")
    expected = build_case_contract(mu_values, beta_values)
    prepared = json.loads(path.read_text(encoding="utf-8"))
    if prepared != rlb2c.rlb2b._json_value(expected):
        raise RuntimeError("Worker grids or contract differ from preparation.")
    return expected, rlb2c.rlb2b._sha256(path)


def solve_sweep_model(
    sweep: str,
    model: str,
    *,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    mu_values: Sequence[float] | None = None,
    beta_values: Sequence[float] | None = None,
    point_start: int = 0,
    point_stop: int | None = None,
) -> list[dict[str, Any]]:
    """Run one of the six independent root workers.

    ``point_start``/``point_stop`` provide a sequential, resumable half-open
    slice of the already prepared full grid.  They do not change the contract
    or final CSV semantics.  Concurrent slices must not write the same CSV.
    """

    if sweep not in SWEEPS:
        raise ValueError(f"Unknown sweep: {sweep!r}.")
    if model not in MODELS:
        raise ValueError(f"Unknown model: {model!r}.")
    active_mu = (
        mu_grid() if mu_values is None else np.asarray(mu_values, dtype=float)
    )
    active_beta = (
        beta_grid()
        if beta_values is None
        else np.asarray(beta_values, dtype=float)
    )
    _contract, contract_sha256 = _worker_contract(
        Path(output_dir), active_mu, active_beta
    )
    policy = rlb2c.rlb2b.frozen_root_policy()
    target = Path(output_dir) / ROOT_FILENAMES[(sweep, model)]

    if sweep == SWEEP_MU:
        points = [
            geometry_for(float(value), MU_TAU, MU_BETA_DEG)
            for value in active_mu
        ]
        shared_objects: ModelObjects | None = None
        shared_cache: ArmPairCache | None = None
    else:
        points = [
            geometry_for(BETA_MU, BETA_TAU, float(value))
            for value in active_beta
        ]
        shared_objects = build_model_objects(points[0])
        if constitutive_checks(points[0], shared_objects)["status"] != "PASS":
            raise RuntimeError("Arm contract failed before the beta worker.")
        shared_cache = ArmPairCache.empty()

    start = int(point_start)
    stop = len(points) if point_stop is None else int(point_stop)
    if start < 0 or stop < start or stop > len(points):
        raise ValueError(
            f"Require 0<=point_start<=point_stop<={len(points)}."
        )
    full_worker = start == 0 and stop == len(points)
    rows = (
        []
        if start == 0 or full_worker or not target.is_file()
        else [dict(row) for row in rlb2c.rlb2b._read_csv(target)]
    )
    selected_points = points[start:stop]

    for local_index, geometry in enumerate(selected_points, start=1):
        point_index = start + local_index
        if shared_objects is None:
            objects = build_model_objects(geometry)
            if constitutive_checks(geometry, objects)["status"] != "PASS":
                raise RuntimeError(
                    f"Arm contract failed before roots at {geometry}."
                )
            point_cache = ArmPairCache.empty()
        else:
            objects = shared_objects
            point_cache = shared_cache
            if point_cache is None:  # pragma: no cover - defensive invariant
                raise RuntimeError("Missing shared beta-sweep arm cache.")

        if model == MODEL_EB:
            raw_rows, export = _eb_sign_scan_rows(
                geometry, objects, policy, cache=point_cache
            )
        else:
            raw_rows, export = _inventory_rows(
                model,
                geometry,
                objects,
                policy,
                contract_sha256,
                point_cache,
            )
        point_rows = [
            transform_root_row(row, model, geometry, sweep, export)
            for row in raw_rows
        ]
        axis_field = "mu" if sweep == SWEEP_MU else "beta_deg"
        axis_value = geometry.mu if sweep == SWEEP_MU else geometry.beta_deg
        rows = [
            row
            for row in rows
            if round(float(row[axis_field]), 12)
            != round(float(axis_value), 12)
        ]
        rows.extend(point_rows)
        rows.sort(
            key=lambda row: (
                float(row[axis_field]), int(row["sorted_position"])
            )
        )
        rlb2c.rlb2b._write_csv(target, rows)
        print(
            f"[{sweep}/{model}] {sweep}={axis_value:g} "
            f"({point_index}/{len(points)}; slice {start}:{stop}): "
            f"export={export['export_range_status']}, "
            f"internal={export['internal_inventory_status']}, "
            f"Lambda_9={float(point_rows[-1]['Lambda']):.12g}",
            flush=True,
        )
    return rows


def _atomic_write_csv(
    path: Path, rows: Sequence[Mapping[str, Any]]
) -> None:
    """Write a complete CSV beside its target, then replace atomically."""

    target = Path(path)
    temporary = target.with_name(f".{target.name}.closing.tmp")
    try:
        rlb2c.rlb2b._write_csv(temporary, rows)
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    """Write JSON beside its target, then replace atomically."""

    target = Path(path)
    temporary = target.with_name(f".{target.name}.closing.tmp")
    try:
        rlb2c.rlb2b._write_json(temporary, payload)
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def _complete_mu_values(
    rows: Sequence[Mapping[str, Any]], model: str
) -> list[float]:
    """Validate unique 1--8 plus root-9 groups and return their mu values."""

    groups: dict[float, list[Mapping[str, Any]]] = {}
    keys: set[tuple[str, str, float, int]] = set()
    for row in rows:
        mu_value = round(float(row["mu"]), 12)
        position = int(row["sorted_position"])
        key = (str(row["sweep"]), str(row["model"]), mu_value, position)
        if key in keys:
            raise RuntimeError(f"Duplicate closing key: {key}.")
        keys.add(key)
        if row["sweep"] != SWEEP_MU or row["model"] != model:
            raise RuntimeError(f"Unexpected row contract for key {key}.")
        groups.setdefault(mu_value, []).append(row)
    for mu_value, point_rows in groups.items():
        positions = sorted(int(row["sorted_position"]) for row in point_rows)
        if positions != list(range(1, OUTPUT_GUARD_POSITION + 1)):
            raise RuntimeError(
                f"Incomplete mu={mu_value:g}, model={model}: {positions}."
            )
        guard_rows = [row for row in point_rows if str(row["guard_flag"]).lower() == "true"]
        if len(guard_rows) != 1 or int(guard_rows[0]["sorted_position"]) != 9:
            raise RuntimeError(
                f"Invalid root-9 guard at mu={mu_value:g}, model={model}."
            )
    return sorted(groups)


def _canonical_mu_rows(output_dir: Path, model: str) -> list[dict[str, Any]]:
    path = Path(output_dir) / ROOT_FILENAMES[(SWEEP_MU, model)]
    return (
        [dict(row) for row in rlb2c.rlb2b._read_csv(path)]
        if path.is_file()
        else []
    )


def _merge_disjoint_or_identical_mu_rows(
    paths: Sequence[Path],
    model: str,
) -> tuple[list[dict[str, Any]], dict[str, int], int]:
    """Merge point files while ignoring byte-equivalent historical overlaps.

    Closing shards may remain after their rows have already been consolidated
    into the canonical CSV.  Such overlaps are provenance duplicates, not new
    spectral points.  A conflicting duplicate remains a hard error.
    """

    merged: dict[tuple[str, str, float, int], dict[str, Any]] = {}
    sources: dict[str, int] = {}
    identical_overlap_count = 0
    common_parent = Path(paths[0]).parent if paths else Path(".")
    for path in paths:
        if not path.is_file():
            continue
        data = [dict(row) for row in rlb2c.rlb2b._read_csv(path)]
        try:
            label = str(path.relative_to(common_parent))
        except ValueError:
            label = str(path)
        sources[label] = len(data)
        for row in data:
            key = (
                str(row["sweep"]),
                str(row["model"]),
                round(float(row["mu"]), 12),
                int(row["sorted_position"]),
            )
            if key in merged:
                if merged[key] != row:
                    raise RuntimeError(
                        f"Conflicting duplicate closing key: {key}."
                    )
                identical_overlap_count += 1
                continue
            if row["sweep"] != SWEEP_MU or row["model"] != model:
                raise RuntimeError(f"Unexpected row contract for key {key}.")
            merged[key] = row
    rows = sorted(
        merged.values(),
        key=lambda row: (float(row["mu"]), int(row["sorted_position"])),
    )
    return rows, sources, identical_overlap_count


def closing_mu_audit(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    *,
    include_rlb_shards: bool = False,
) -> dict[str, Any]:
    """Audit existing mu rows without running a determinant calculation."""

    target = Path(output_dir)
    grid = [round(float(value), 12) for value in mu_grid()]
    result: dict[str, Any] = {}
    for model in (MODEL_OLD, MODEL_RLB):
        paths = [target / ROOT_FILENAMES[(SWEEP_MU, model)]]
        if model == MODEL_RLB and include_rlb_shards:
            paths.extend(
                sorted(
                    (target / "_shards").glob(
                        "mu_rlb_*/mu_sweep_new_rlb_roots.csv"
                    )
                )
            )
        rows, sources, identical_overlaps = _merge_disjoint_or_identical_mu_rows(
            paths, model
        )
        complete = _complete_mu_values(rows, model)
        result[model] = {
            "source_rows": sources,
            "row_count": len(rows),
            "complete_point_count": len(complete),
            "complete_mu": complete,
            "missing_mu": [value for value in grid if value not in complete],
            "duplicate_key_count": 0,
            "ignored_identical_overlap_count": identical_overlaps,
            "incomplete_group_count": 0,
        }
    return result


def consolidate_existing_rlb_mu_rows(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> dict[str, Any]:
    """Atomically merge disjoint completed RLB shards into the canonical CSV."""

    target = Path(output_dir)
    canonical = target / ROOT_FILENAMES[(SWEEP_MU, MODEL_RLB)]
    paths = [canonical]
    paths.extend(
        sorted(
            (target / "_shards").glob(
                "mu_rlb_*/mu_sweep_new_rlb_roots.csv"
            )
        )
    )
    rows, source_rows, identical_overlaps = _merge_disjoint_or_identical_mu_rows(
        paths, MODEL_RLB
    )
    complete = _complete_mu_values(rows, MODEL_RLB)
    _atomic_write_csv(canonical, rows)
    return {
        "source_rows": source_rows,
        "ignored_identical_overlap_count": identical_overlaps,
        "merged_row_count": len(rows),
        "merged_point_count": len(complete),
        "merged_mu": complete,
        "canonical_sha256": rlb2c.rlb2b._sha256(canonical),
    }


def _point_qualification(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    first = rows[0]
    export_pass = str(first["export_range_status"]) == "PASS"
    first8_quality = all(
        str(row["root_status"]) == "PASS"
        for row in rows
        if int(row["sorted_position"]) <= PLOTTED_POSITIONS
    )
    guard_quality = all(
        str(row["root_status"]) == "PASS"
        for row in rows
        if int(row["sorted_position"]) == OUTPUT_GUARD_POSITION
    )
    detailed = str(
        first.get("plotted_first8_primary_verification_agreement", "")
    ) not in {"", "None"}
    if detailed:
        first8 = (
            "PASS"
            if first8_quality
            and str(first["plotted_first8_primary_verification_agreement"]).lower()
            == "true"
            else "FAIL"
        )
        guard = (
            "PASS"
            if guard_quality
            and str(first["root9_primary_verification_agreement"]).lower()
            == "true"
            else "FAIL"
        )
    elif export_pass and first8_quality and guard_quality:
        first8 = "PASS"
        guard = "PASS"
    elif not first8_quality or not guard_quality:
        first8 = "FAIL" if not first8_quality else "UNRESOLVED_AGGREGATE_POSITIONS_1_TO_9"
        guard = "FAIL" if not guard_quality else "UNRESOLVED_AGGREGATE_POSITIONS_1_TO_9"
    else:
        first8 = "UNRESOLVED_AGGREGATE_POSITIONS_1_TO_9"
        guard = "UNRESOLVED_AGGREGATE_POSITIONS_1_TO_9"
    internal = str(first["internal_inventory_status"])
    model = str(first["model"])
    sweep = str(first["sweep"])
    parameter_name = "mu" if sweep == SWEEP_MU else "beta_deg"
    return {
        "sweep": sweep,
        "parameter_name": parameter_name,
        "parameter_value": float(first[parameter_name]),
        "mu": float(first["mu"]),
        "model": model,
        "PLOTTED_FIRST_8": first8,
        "ROOT_9_GUARD": guard,
        "INTERNAL_TAIL": (
            "NOT_APPLICABLE_SIGN_SCAN"
            if model == MODEL_EB
            else "PASS"
            if internal == "PASS"
            else "QUALIFIED"
        ),
        "aggregate_export_range_status": str(first["export_range_status"]),
        "aggregate_primary_verification_max_relative": float(
            first["export_primary_verification_max_relative"]
        ),
        "internal_inventory_status": internal,
    }


def update_closing_checkpoint(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    *,
    inherited_workers_drained: int | None = None,
) -> dict[str, Any]:
    """Atomically refresh the closing checkpoint from canonical mu CSVs."""

    target = Path(output_dir)
    path = target / CLOSING_CHECKPOINT_FILENAME
    checkpoint = (
        json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {}
    )
    audit = closing_mu_audit(target, include_rlb_shards=False)
    initial_missing = checkpoint.get("initial_missing_mu", {})
    completed: dict[str, list[float]] = {}
    qualifications: list[dict[str, Any]] = []
    for model in (MODEL_OLD, MODEL_RLB):
        initial_values = initial_missing.get(model)
        if initial_values is None and model == MODEL_RLB:
            # The first closing checkpoint used the short provenance label
            # ``NEW_RLB``.  Accept it once and normalize all later writes to
            # the canonical model label used by the CSV inventory.
            initial_values = initial_missing.get("NEW_RLB", [])
        initial = [round(float(value), 12) for value in (initial_values or [])]
        missing = set(audit[model]["missing_mu"])
        completed[model] = [value for value in initial if value not in missing]
        rows = _canonical_mu_rows(target, model)
        groups: dict[float, list[Mapping[str, Any]]] = {}
        for row in rows:
            groups.setdefault(round(float(row["mu"]), 12), []).append(row)
        for value in completed[model]:
            qualifications.append(_point_qualification(groups[value]))
    checkpoint.update(
        {
            "initial_missing_mu": {
                MODEL_OLD: [
                    round(float(value), 12)
                    for value in initial_missing.get(MODEL_OLD, [])
                ],
                MODEL_RLB: [
                    round(float(value), 12)
                    for value in (
                        initial_missing.get(MODEL_RLB)
                        or initial_missing.get("NEW_RLB", [])
                    )
                ],
            },
            "completed_new_mu_points": completed,
            "current_missing_mu": {
                model: audit[model]["missing_mu"]
                for model in (MODEL_OLD, MODEL_RLB)
            },
            "new_points_executed": sum(len(values) for values in completed.values()),
            "ready_points_recalculated": 0,
            "parallel_workers_used_in_closing_stage": 0,
            "global_restarts": 0,
            "latest_canonical_audit": audit,
            "new_point_qualifications": qualifications,
        }
    )
    if inherited_workers_drained is not None:
        checkpoint["inherited_workers_drained"] = int(inherited_workers_drained)
    _atomic_write_json(path, checkpoint)
    return checkpoint


def solve_one_missing_mu_point(
    model: str,
    point_index: int,
    *,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> dict[str, Any]:
    """Solve exactly one absent mu point and atomically merge its nine rows."""

    if model not in {MODEL_OLD, MODEL_RLB}:
        raise ValueError("Closing point model must be OLD_TIMOSHENKO or NEW_RLB.")
    bad_limits = {
        name: os.environ.get(name)
        for name, expected in CLOSING_THREAD_LIMITS.items()
        if os.environ.get(name) != expected
    }
    if bad_limits:
        raise RuntimeError(f"Closing thread limits are not fixed at one: {bad_limits}.")
    active_mu = mu_grid()
    if point_index < 0 or point_index >= len(active_mu):
        raise ValueError(f"Invalid mu point index: {point_index}.")
    mu_value = round(float(active_mu[point_index]), 12)
    target_dir = Path(output_dir)
    existing = _canonical_mu_rows(target_dir, model)
    complete = _complete_mu_values(existing, model)
    if mu_value in complete:
        raise RuntimeError(
            f"Refusing to recalculate ready point model={model}, mu={mu_value:g}."
        )
    _contract, contract_sha256 = _worker_contract(
        target_dir, active_mu, beta_grid()
    )
    geometry = geometry_for(mu_value, MU_TAU, MU_BETA_DEG)
    objects = build_model_objects(geometry)
    if constitutive_checks(geometry, objects)["status"] != "PASS":
        raise RuntimeError(f"Arm contract failed before roots at {geometry}.")
    policy = rlb2c.rlb2b.frozen_root_policy()
    raw_rows, export = _inventory_rows(
        model,
        geometry,
        objects,
        policy,
        contract_sha256,
        ArmPairCache.empty(),
    )
    point_rows = [
        transform_root_row(row, model, geometry, SWEEP_MU, export)
        for row in raw_rows
    ]
    _complete_mu_values(point_rows, model)
    merged = [*existing, *point_rows]
    _complete_mu_values(merged, model)
    merged.sort(key=lambda row: (float(row["mu"]), int(row["sorted_position"])))
    target = target_dir / ROOT_FILENAMES[(SWEEP_MU, model)]
    _atomic_write_csv(target, merged)
    checkpoint = update_closing_checkpoint(target_dir)
    qualification = _point_qualification(point_rows)
    print(
        f"[closing/{model}] mu={mu_value:g}: "
        f"export={qualification['aggregate_export_range_status']}, "
        f"internal={qualification['INTERNAL_TAIL']}; checkpoint updated",
        flush=True,
    )
    return {
        "model": model,
        "mu": mu_value,
        "row_count": len(point_rows),
        "qualification": qualification,
        "checkpoint_new_points": checkpoint["new_points_executed"],
    }


def prepare_output(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    mu_values: Sequence[float] | None = None,
    beta_values: Sequence[float] | None = None,
) -> dict[str, Any]:
    """Write the shared contract, sanity tables, and frozen-input manifest."""

    active_mu = (
        mu_grid() if mu_values is None else np.asarray(mu_values, dtype=float)
    )
    active_beta = (
        beta_grid()
        if beta_values is None
        else np.asarray(beta_values, dtype=float)
    )
    target = Path(output_dir)
    contract = build_case_contract(active_mu, active_beta)
    reference_geometry = geometry_for(0.0, 0.0, MU_BETA_DEG)
    reference_objects = build_model_objects(reference_geometry)
    checks = constitutive_checks(reference_geometry, reference_objects)
    equivalence = eb_tau0_equivalence_diagnostic()
    if checks["status"] != "PASS" or equivalence["status"] != "PASS":
        raise RuntimeError(
            "Material/geometry or additive EB equivalence failed before roots."
        )
    contract_path = target / "case_contract.json"
    rlb2c.rlb2b._write_json(contract_path, contract)
    rlb2c.rlb2b._write_csv(
        target / "geometry_sanity_checks.csv", geometry_sanity_rows()
    )
    rlb2c.rlb2b._write_csv(
        target / "laminate_properties_summary.csv", laminate_property_rows()
    )
    manifest = {
        "schema_version": 1,
        "algorithm_version": ALGORITHM_VERSION,
        "stage": STAGE_ID,
        "scientific_role": (
            "two finite plots of independently sorted spectral positions"
        ),
        "git_at_preparation": rlb2c.rlb2b._git_state(),
        "case_contract_sha256": rlb2c.rlb2b._sha256(contract_path),
        "constitutive_checks": checks,
        "eb_tau0_equivalence": equivalence,
        "models": {
            MODEL_EB: {
                "physics": "isotropic rectangular Euler--Bernoulli",
                "mu_sweep_root_workflow": "frozen EB sign scan",
                "beta_sweep_root_workflow": (
                    "physical unequal-section determinant plus frozen EB sign scan"
                ),
            },
            MODEL_OLD: {
                "physics": "frozen isotropic rectangular Timoshenko",
                "matrix": (
                    "scripts/lib/isotropic_rectangular_"
                    "timoshenko_coupled_beams.py"
                ),
            },
            MODEL_RLB: {
                "physics": "frozen weakly orthotropic four-ply RLB case A",
                "matrix": "scripts/lib/reddy_symmetric_coupled_beams.py",
                "delta": DELTA,
                "one_ply_shortcut": False,
            },
        },
        "root_search_coordinates": {
            "mu_EB": "historical Lambda=sqrt(Omega)",
            "beta_EB": "Omega",
            "old_Timoshenko": "Omega",
            "weak_RLB": "Omega",
        },
        "frozen_root_policy": asdict(rlb2c.rlb2b.frozen_root_policy()),
        "thresholds": {
            "normalization_identity_relative": (
                NORMALIZATION_IDENTITY_TOLERANCE
            ),
            "root_singular_ratio": (
                rlb2c.rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE
            ),
            "boundary_null_residual": (
                rlb2c.rlb2b.BOUNDARY_RESIDUAL_TOLERANCE
            ),
            "inventory_primary_verification_relative": (
                INVENTORY_VERIFICATION_TOLERANCE
            ),
            "EB_determinant_identity_relative": (
                EB_DETERMINANT_IDENTITY_TOLERANCE
            ),
        },
        "frozen_model_hashes_before": _frozen_hashes(),
        "root_frequencies_cross_seeded_between_models": False,
        "inter_model_relative_differences_computed": False,
        "new_root_solver": False,
        "new_production_solver": False,
        "explicit_exclusions": [
            "relative_model_differences",
            "equivalent_isotropic_fitting",
            "delta_values_other_than_0.1",
            "other_stacking_sequences",
            "branch_tracking",
            "MAC",
            "veering",
            "Ritz",
            "FEM",
            "torsion",
            "damping",
            "article_preparation",
            "commit",
            "push",
        ],
    }
    rlb2c.rlb2b._write_json(target / "model_manifest.json", manifest)
    return manifest


def _float(row: Mapping[str, Any], key: str) -> float:
    return float(row[key])


def _bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text == "true":
        return True
    if text == "false":
        return False
    raise ValueError(f"Expected serialized boolean; got {value!r}.")


def validate_root_rows(
    sweep: str,
    model: str,
    rows: Sequence[Mapping[str, Any]],
    axis_values: Sequence[float],
) -> dict[str, Any]:
    """Validate the exported first eight roots and root-9 guard."""

    if sweep not in SWEEPS or model not in MODELS:
        raise ValueError("Unknown sweep/model pair.")
    axis_field = "mu" if sweep == SWEEP_MU else "beta_deg"
    expected_axis = [round(float(value), 12) for value in axis_values]
    expected_keys = {
        (axis, position)
        for axis in expected_axis
        for position in range(1, OUTPUT_GUARD_POSITION + 1)
    }
    actual_keys = {
        (
            round(_float(row, axis_field), 12),
            int(row["sorted_position"]),
        )
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
            "export_primary_verification_max_relative",
        )
    )
    positive = all(
        _float(row, field) > 0.0
        for row in rows
        for field in ("omega", "Omega", "Lambda")
    )
    sorted_positions = True
    for axis in expected_axis:
        point_rows = sorted(
            (
                row
                for row in rows
                if round(_float(row, axis_field), 12) == axis
            ),
            key=lambda row: int(row["sorted_position"]),
        )
        sorted_positions = sorted_positions and len(point_rows) == 9 and all(
            _float(right, "Omega") >= _float(left, "Omega")
            for left, right in zip(point_rows, point_rows[1:])
        )
    guard_structure = all(
        (_bool(row["guard_flag"]))
        == (int(row["sorted_position"]) == OUTPUT_GUARD_POSITION)
        and row["role"]
        == (
            "ROOT_9_GUARD"
            if int(row["sorted_position"]) == OUTPUT_GUARD_POSITION
            else "FIRST_8"
        )
        for row in rows
    )
    accepted_quality = all(
        row["root_status"] == "PASS"
        and _float(row, "scaled_sigma_ratio")
        <= rlb2c.rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE
        and _float(row, "boundary_null_residual")
        <= rlb2c.rlb2b.BOUNDARY_RESIDUAL_TOLERANCE
        and _float(row, "normalization_identity_relative_residual")
        <= NORMALIZATION_IDENTITY_TOLERANCE
        for row in rows
    )
    if model == MODEL_EB:
        accepted_quality = accepted_quality and all(
            _float(row, "historical_EB_mapping_relative_residual")
            <= NORMALIZATION_IDENTITY_TOLERANCE
            for row in rows
        )
    export_agreement = all(
        _bool(row["export_primary_verification_agreement"])
        and row["export_range_status"] == "PASS"
        for row in rows
    )
    export_issue_axes = [
        axis
        for axis in expected_axis
        if not all(
            _bool(row["export_primary_verification_agreement"])
            and row["export_range_status"] == "PASS"
            for row in rows
            if round(_float(row, axis_field), 12) == axis
        )
    ]
    no_unresolved_below_export = all(
        str(row["unresolved_candidates_below_export_guard"])
        == "NOT_ASSESSED_BY_EB_SIGN_SCAN"
        or int(row["unresolved_candidates_below_export_guard"]) == 0
        for row in rows
    )
    hashes_by_axis = {
        axis: {
            str(row["inventory_sha256"])
            for row in rows
            if round(_float(row, axis_field), 12) == axis
        }
        for axis in expected_axis
    }
    inventory_hash_structure = bool(
        all(len(values) == 1 for values in hashes_by_axis.values())
        and len(
            {
                next(iter(values))
                for values in hashes_by_axis.values()
                if values
            }
        )
        == len(expected_axis)
    )
    expected_rows = len(expected_axis) * OUTPUT_GUARD_POSITION
    exact_structure = bool(
        len(rows) == expected_rows
        and actual_keys == expected_keys
        and finite
        and positive
        and sorted_positions
        and guard_structure
        and inventory_hash_structure
    )
    export_pass = bool(
        exact_structure
        and accepted_quality
        and export_agreement
        and no_unresolved_below_export
    )
    point_internal_statuses = {
        axis: next(
            (
                str(row["internal_inventory_status"])
                for row in rows
                if round(_float(row, axis_field), 12) == axis
            ),
            "MISSING",
        )
        for axis in expected_axis
    }
    accepted_internal = {"PASS", "PASS_SIGN_SCAN_CROSSCHECK"}
    tail_issue_axes = [
        axis
        for axis, status in point_internal_statuses.items()
        if status not in accepted_internal
    ]
    usable_prefix = bool(
        exact_structure and accepted_quality and no_unresolved_below_export
    )
    status = (
        "PASS"
        if export_pass and not tail_issue_axes
        else "PARTIAL_PASS"
        if usable_prefix
        else "FAIL"
    )
    return {
        "summary_kind": "SWEEP_MODEL_ROOT_INVENTORY",
        "sweep": sweep,
        "model": model,
        "axis_count": len(expected_axis),
        "expected_rows": expected_rows,
        "actual_rows": len(rows),
        "first_frequency_count": len(expected_axis) * PLOTTED_POSITIONS,
        "guard_count": sum(_bool(row["guard_flag"]) for row in rows),
        "maximum_scaled_sigma_ratio": max(
            (_float(row, "scaled_sigma_ratio") for row in rows),
            default=math.inf,
        ),
        "maximum_boundary_null_residual": max(
            (_float(row, "boundary_null_residual") for row in rows),
            default=math.inf,
        ),
        "maximum_normalization_identity_relative_residual": max(
            (
                _float(row, "normalization_identity_relative_residual")
                for row in rows
            ),
            default=math.inf,
        ),
        "maximum_export_primary_verification_relative": max(
            (
                _float(row, "export_primary_verification_max_relative")
                for row in rows
            ),
            default=math.inf,
        ),
        "export_issue_axis_values": export_issue_axes,
        "export_issue_count": len(export_issue_axes),
        "tail_issue_axis_values": tail_issue_axes,
        "tail_issue_count": len(tail_issue_axes),
        "exact_row_structure_passed": exact_structure,
        "accepted_root_quality_passed": accepted_quality,
        "export_primary_verification_passed": export_agreement,
        "no_unresolved_below_export_guard_passed": (
            no_unresolved_below_export
        ),
        "export_range_status": "PASS" if export_pass else "FAIL",
        "status": status,
    }


def create_plot(
    rows_by_model: Mapping[str, Sequence[Mapping[str, Any]]],
    output_path: Path,
    *,
    axis_field: str,
    x_label: str,
    x_limits: tuple[float, float],
) -> Path:
    """Create one 24-curve sorted-position plot."""

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
                key=lambda row: float(row[axis_field]),
            )
            ax.plot(
                [float(row[axis_field]) for row in selected],
                [float(row["Lambda"]) for row in selected],
                color=colors[position - 1],
                linestyle=LINE_STYLES[model],
                linewidth=1.55 if model != MODEL_RLB else 1.35,
                alpha=0.96,
            )
    ax.set_xlabel(x_label)
    ax.set_ylabel("Lambda")
    ax.set_xlim(*x_limits)
    ax.grid(True, alpha=0.24)
    model_legend = ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                color="black",
                linestyle=LINE_STYLES[model],
                linewidth=1.7,
                label=MODEL_LABELS[model],
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


def collect_point_qualifications(
    rows_by_pair: Mapping[
        tuple[str, str], Sequence[Mapping[str, Any]]
    ],
) -> list[dict[str, Any]]:
    """Classify the plotted positions, root 9, and internal tail separately."""

    qualifications: list[dict[str, Any]] = []
    for sweep in SWEEPS:
        axis_field = "mu" if sweep == SWEEP_MU else "beta_deg"
        for model in MODELS:
            groups: dict[float, list[Mapping[str, Any]]] = {}
            for row in rows_by_pair[(sweep, model)]:
                groups.setdefault(
                    round(float(row[axis_field]), 12), []
                ).append(row)
            for axis_value in sorted(groups):
                point_rows = sorted(
                    groups[axis_value],
                    key=lambda row: int(row["sorted_position"]),
                )
                if [int(row["sorted_position"]) for row in point_rows] != list(
                    range(1, OUTPUT_GUARD_POSITION + 1)
                ):
                    raise RuntimeError(
                        f"Incomplete qualification group {sweep}/{model}/"
                        f"{axis_value:g}."
                    )
                qualifications.append(_point_qualification(point_rows))
    return qualifications


def _closing_statuses(
    qualifications: Sequence[Mapping[str, Any]], plot_pass: bool
) -> dict[str, str]:
    beta = [item for item in qualifications if item["sweep"] == SWEEP_BETA]
    mu = [item for item in qualifications if item["sweep"] == SWEEP_MU]
    first8_beta_pass = bool(
        beta and all(item["PLOTTED_FIRST_8"] == "PASS" for item in beta)
    )
    first8_mu_pass = bool(
        mu and all(item["PLOTTED_FIRST_8"] == "PASS" for item in mu)
    )
    root9_pass = bool(
        qualifications
        and all(item["ROOT_9_GUARD"] == "PASS" for item in qualifications)
    )
    tail_qualified = any(
        item["INTERNAL_TAIL"] == "QUALIFIED"
        and item["aggregate_export_range_status"] == "PASS"
        for item in qualifications
    )
    statuses = {
        "RLB-2D-BETA-PLOTTED-FIRST-8": (
            "PASS" if first8_beta_pass else "FAIL"
        ),
        "RLB-2D-MU-PLOTTED-FIRST-8": (
            "PASS" if first8_mu_pass else "PARTIAL_PASS"
        ),
        "RLB-2D-ROOT9-GUARDS": (
            "PASS" if root9_pass else "PARTIAL_PASS"
        ),
        "RLB-2D-INTERNAL-TAIL": "QUALIFIED" if tail_qualified else "PASS",
        "RLB-2D-PLOT-GENERATION": "PASS" if plot_pass else "FAIL",
    }
    statuses["SCIENTIFIC_OVERALL"] = (
        "PASS_WITH_INTERNAL_TAIL_QUALIFICATIONS"
        if first8_beta_pass
        and first8_mu_pass
        and root9_pass
        and plot_pass
        else "PARTIAL_PASS_PLOTTED_RANGE_QUALIFICATIONS"
    )
    return statuses


def _sweep_status(
    summaries: Sequence[Mapping[str, Any]], sweep: str
) -> str:
    selected = [summary for summary in summaries if summary["sweep"] == sweep]
    if len(selected) != len(MODELS) or any(
        summary["status"] == "FAIL" for summary in selected
    ):
        return "FAIL"
    return (
        "PASS"
        if all(summary["status"] == "PASS" for summary in selected)
        else "PARTIAL_PASS"
    )


def _report_text(
    contract: Mapping[str, Any],
    summaries: Sequence[Mapping[str, Any]],
    statuses: Mapping[str, str],
    qualifications: Sequence[Mapping[str, Any]],
    closing_checkpoint: Mapping[str, Any],
) -> str:
    reference = contract["reference_constants"]
    mu_spec = contract["sweeps"][SWEEP_MU]
    beta_spec = contract["sweeps"][SWEEP_BETA]
    lines = [
        "# Сравнение спектральных позиций трёх моделей",
        "",
        "## Расчётные случаи",
        "",
        (
            "Рассматриваются две конечные параметрические сетки. На первой "
            f"сетке $\\beta={MU_BETA_DEG:g}^\\circ$, $\\tau=0$, а $\\mu$ "
            f"изменяется от {float(mu_spec['values'][0]):g} до "
            f"{float(mu_spec['values'][-1]):g} с шагом "
            f"{float(mu_spec['step']):g}. На второй сетке $\\mu=0.5$, "
            "$\\tau=0.2$, а $\\beta$ изменяется от $0$ до $90^\\circ$ "
            f"с шагом {float(beta_spec['step']):g}$^\\circ$."
        ),
        "",
    ]
    if mu_spec["fallback_used"]:
        lines.extend(
            [
                (
                    "Для сетки по $\\mu$ использован разрешённый шаг 0.02 "
                    "вместо исходного шага 0.01. Причина зафиксирована как "
                    f"`{mu_spec['fallback_reason']}`. Пороговые значения и "
                    "root policy не изменялись; выбор шага сделан до анализа "
                    "спектра."
                ),
                "",
            ]
        )
    lines.extend(
        [
            (
                "Линии соединяют независимо отсортированные спектральные "
                "позиции. Они не задают продолжение форм колебаний и не "
                "используются как идентификаторы модальных ветвей."
            ),
            "",
            "## Нормировка",
            "",
            "$$",
            "\\Omega=\\omega l^2\\sqrt{\\frac{\\rho_0 A_0}{E_0 I_{y0}}},",
            "\\qquad \\Lambda=\\sqrt{\\Omega},",
            "\\qquad A_0=b_0h_0,\\quad I_{y0}=\\frac{b_0h_0^3}{12}.",
            "$$",
            "",
            (
                f"Использованы $E_0={float(reference['E0']):g}$, "
                f"$\\nu_0={float(reference['nu0']):g}$, "
                f"$\\rho_0={float(reference['rho0']):g}$, "
                f"$b_0={float(reference['b0']):g}$, "
                f"$h_0={float(reference['h0']):g}$ и "
                f"$l={float(reference['l']):g}$. Одна reference normalization "
                "применена ко всем моделям и всем точкам."
            ),
            "",
            "## Геометрия и weakly orthotropic case A",
            "",
            (
                "Фактическая параметризация имеет вид "
                "$L_1=l(1-\\mu)$, $L_2=l(1+\\mu)$, "
                "$h_1=h_0(1-\\tau)$, $h_2=h_0(1+\\tau)$ и "
                "$b_1=b_2=b_0$."
            ),
            "",
            (
                "Каждое плечо RLB содержит четыре равных слоя "
                "$[0/90/90/0]$ толщиной $h_i/4$. Приняты "
                "$\\delta=0.1$, $E_1=1.1$, $E_2=0.9$, "
                "$\\nu_{12}=0.3$, "
                "$G_{12}=G_{13}=G_{23}=1/2.6$ и $\\rho=1$. "
                "One-ply shortcut и equivalent-isotropic fitting не "
                "использовались."
            ),
            "",
            "## Root data",
            "",
            (
                "| sweep | модель | точек | строк | first 8 | guards | "
                "tail issues | status |"
            ),
            "|---|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for summary in summaries:
        lines.append(
            f"| {summary['sweep']} | {summary['model']} | "
            f"{int(summary['axis_count'])} | {int(summary['actual_rows'])} | "
            f"{int(summary['first_frequency_count'])} | "
            f"{int(summary['guard_count'])} | "
            f"{int(summary['tail_issue_count'])} | {summary['status']} |"
        )
    tail_summaries = [
        summary for summary in summaries if summary["tail_issue_count"]
    ]
    export_issue_summaries = [
        summary
        for summary in summaries
        if summary["export_range_status"] != "PASS"
    ]
    if export_issue_summaries:
        lines.extend(
            [
                "",
                (
                    "Следующие inventories имеют отдельное замечание к "
                    "qualification экспортируемого диапазона 1--9:"
                ),
            ]
        )
        for summary in export_issue_summaries:
            lines.append(
                f"- `{summary['sweep']}/{summary['model']}`: "
                "primary/verification agreement = "
                f"{summary['export_primary_verification_passed']}; "
                "maximum relative difference = "
                f"{float(summary['maximum_export_primary_verification_relative']):.6e}; "
                "axis values = "
                f"{summary['export_issue_axis_values']}."
            )
    if tail_summaries:
        lines.extend(
            [
                "",
                (
                    "Отмеченные internal tail issues относятся к полной "
                    "qualification до внутреннего root 13. Статус "
                    "экспортируемого диапазона устанавливается отдельно по "
                    "первым восьми корням и root 9 guard."
                ),
            ]
        )
    range_qualifications = [
        item
        for item in qualifications
        if item["PLOTTED_FIRST_8"] != "PASS"
        or item["ROOT_9_GUARD"] != "PASS"
    ]
    tail_only_qualifications = [
        item
        for item in qualifications
        if item["aggregate_export_range_status"] == "PASS"
        and item["INTERNAL_TAIL"] == "QUALIFIED"
    ]
    lines.extend(
        [
            "",
            "## Раздельная qualification",
            "",
        ]
    )
    if range_qualifications:
        lines.append(
            "Следующие сохранённые inventories имеют qualification внутри "
            "агрегированного диапазона positions 1--9:"
        )
        lines.append("")
        for item in range_qualifications:
            lines.append(
                f"- `{item['sweep']}/{item['model']}`, "
                f"{item['parameter_name']}={float(item['parameter_value']):g}: "
                f"first 8 = `{item['PLOTTED_FIRST_8']}`, root 9 = "
                f"`{item['ROOT_9_GUARD']}`, maximum relative = "
                f"{float(item['aggregate_primary_verification_max_relative']):.6e}."
            )
        lines.extend(
            [
                "",
                (
                    "Для строк, рассчитанных до closing stage, исходный CSV "
                    "сохраняет только aggregate maximum по positions 1--9. "
                    "Поэтому FAIL этого aggregate нельзя без повторного "
                    "расчёта отнести только к first 8 или только к root 9. "
                    "Такие точки консервативно имеют статус "
                    "`UNRESOLVED_AGGREGATE_POSITIONS_1_TO_9`."
                ),
            ]
        )
    else:
        lines.append("Positions 1--8 и root 9 прошли во всех точках.")
    lines.extend(
        [
            "",
            (
                "Qualifications строго выше root 9 сохранены отдельно: "
                f"{len(tail_only_qualifications)} point inventories. Они не "
                "изменяют статусы plotted range."
            ),
            "",
            "## Closing stage",
            "",
            (
                f"Повторно использовано {int(closing_checkpoint.get('initial_reused_point_count', 0))} "
                "готовых point inventories. В closing stage рассчитано "
                f"{int(closing_checkpoint.get('new_points_executed', 0))} "
                "отсутствовавших points; готовых points повторно не "
                "рассчитывали. После завершения двух унаследованных workers "
                "каждый следующий point вычислялся отдельным последовательным "
                "Python-процессом с BLAS/OpenMP thread limits, равными единице."
            ),
        ]
    )
    lines.extend(
        [
            "",
            "## Графики",
            "",
            (
                "Файлы [lambda_vs_mu_plot.png](lambda_vs_mu_plot.png) и "
                "[lambda_vs_beta_plot.png](lambda_vs_beta_plot.png) содержат "
                "по 24 кривые. Цвет задаёт sorted position от 1 до 8. "
                "Сплошные линии соответствуют EB, пунктирные линии — old "
                "Timoshenko, штрихпунктирные линии — weak RLB."
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
                "Результат относится только к двум указанным конечным "
                "сеткам, case A и $\\delta=0.1$. Относительные расхождения "
                "между моделями не вычислялись. Не выполнялись branch "
                "tracking, MAC, veering analysis, Ritz и FEM."
            ),
            "",
            (
                "Для beta sweep с $\\tau=0.2$ frozen EB API не содержит "
                "двух arm-specific thickness arguments. Поэтому применён "
                "additive physical-state EB adapter и frozen EB sign scan. "
                "При $\\tau=0$ adapter проверен по точному "
                "determinant identity с frozen EB matrix. Production physics "
                "modules и root thresholds не изменялись."
            ),
            "",
            (
                "До closing stage часть mu inventory weak RLB вычислялась "
                "неперекрывающимися slices. Closing audit подтвердил "
                "отсутствие дубликатов и неполных групп. Незавершённые "
                "resource-control attempts не входят в итоговые inventories."
            ),
        ]
    )
    return "\n".join(lines) + "\n"


def finalize_output(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    mu_values: Sequence[float] | None = None,
    beta_values: Sequence[float] | None = None,
) -> dict[str, Any]:
    """Validate six root files and create exactly the two required PNGs."""

    active_mu = (
        mu_grid() if mu_values is None else np.asarray(mu_values, dtype=float)
    )
    active_beta = (
        beta_grid()
        if beta_values is None
        else np.asarray(beta_values, dtype=float)
    )
    target = Path(output_dir)
    contract_path = target / "case_contract.json"
    manifest_path = target / "model_manifest.json"
    if not contract_path.is_file() or not manifest_path.is_file():
        raise RuntimeError("Run --mode prepare before finalization.")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    model_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected_contract = build_case_contract(active_mu, active_beta)
    if contract != rlb2c.rlb2b._json_value(expected_contract):
        raise RuntimeError("Finalization grids differ from preparation.")
    contract_hash = rlb2c.rlb2b._sha256(contract_path)
    if model_manifest["case_contract_sha256"] != contract_hash:
        raise RuntimeError("Prepared case-contract hash no longer matches.")
    frozen_after = _frozen_hashes()
    if model_manifest["frozen_model_hashes_before"] != frozen_after:
        raise RuntimeError("A frozen production or predecessor file changed.")

    rows_by_pair = {
        (sweep, model): rlb2c.rlb2b._read_csv(
            target / ROOT_FILENAMES[(sweep, model)]
        )
        for sweep in SWEEPS
        for model in MODELS
    }
    summaries = [
        validate_root_rows(
            sweep,
            model,
            rows_by_pair[(sweep, model)],
            active_mu if sweep == SWEEP_MU else active_beta,
        )
        for sweep in SWEEPS
        for model in MODELS
    ]
    scale = float(contract["normalization"]["Omega_per_omega"])
    normalization_pass = bool(
        math.isfinite(scale)
        and scale > 0.0
        and _relative(scale, reference_omega_scale())
        <= NORMALIZATION_IDENTITY_TOLERANCE
        and _relative(
            old_Lambda(0.731, scale) ** 2,
            omega_to_Omega(0.731, scale),
        )
        <= NORMALIZATION_IDENTITY_TOLERANCE
        and model_manifest["eb_tau0_equivalence"]["status"] == "PASS"
    )
    if not normalization_pass:
        raise RuntimeError("The frozen old-Lambda normalization gate failed.")
    checkpoint = update_closing_checkpoint(
        target, inherited_workers_drained=2
    )
    if any(checkpoint["current_missing_mu"].values()):
        raise RuntimeError(
            f"Cannot finalize with missing mu points: "
            f"{checkpoint['current_missing_mu']}."
        )
    mu_plot = create_plot(
        {model: rows_by_pair[(SWEEP_MU, model)] for model in MODELS},
        target / PLOT_FILENAMES[SWEEP_MU],
        axis_field="mu",
        x_label="mu",
        x_limits=(float(active_mu[0]), float(active_mu[-1])),
    )
    beta_plot = create_plot(
        {model: rows_by_pair[(SWEEP_BETA, model)] for model in MODELS},
        target / PLOT_FILENAMES[SWEEP_BETA],
        axis_field="beta_deg",
        x_label="beta",
        x_limits=(float(active_beta[0]), float(active_beta[-1])),
    )
    plot_pass = all(
        path.is_file() and path.stat().st_size > 0
        for path in (mu_plot, beta_plot)
    )
    all_export_pass = all(
        summary["export_range_status"] == "PASS" for summary in summaries
    )
    exact_mu_grid = any(
        np.array_equal(active_mu, mu_grid(mu_step=step))
        for step in (REQUESTED_MU_STEP, ALLOWED_MU_FALLBACK_STEP)
    )
    exact_beta_grid = np.array_equal(active_beta, beta_grid())
    if not exact_mu_grid or not exact_beta_grid:
        raise RuntimeError("Finalization grid differs from the fixed contract.")
    qualifications = collect_point_qualifications(rows_by_pair)
    statuses = _closing_statuses(qualifications, plot_pass)
    report_path = target / "report.md"
    report_path.write_text(
        _report_text(
            contract,
            summaries,
            statuses,
            qualifications,
            checkpoint,
        ),
        encoding="utf-8",
    )
    generated_names = [
        "model_manifest.json",
        "case_contract.json",
        "geometry_sanity_checks.csv",
        "laminate_properties_summary.csv",
        *[ROOT_FILENAMES[(sweep, model)] for sweep in SWEEPS for model in MODELS],
        PLOT_FILENAMES[SWEEP_MU],
        PLOT_FILENAMES[SWEEP_BETA],
        CLOSING_CHECKPOINT_FILENAME,
        "report.md",
    ]
    plotted_range_qualifications = [
        item
        for item in qualifications
        if item["PLOTTED_FIRST_8"] != "PASS"
        or item["ROOT_9_GUARD"] != "PASS"
    ]
    tail_only_qualifications = [
        item
        for item in qualifications
        if item["aggregate_export_range_status"] == "PASS"
        and item["INTERNAL_TAIL"] == "QUALIFIED"
    ]
    run_manifest = {
        "schema_version": 1,
        "algorithm_version": ALGORITHM_VERSION,
        "stage": STAGE_ID,
        "finalization_git_state": rlb2c.rlb2b._git_state(),
        "case_contract_sha256": contract_hash,
        "model_manifest_sha256": rlb2c.rlb2b._sha256(manifest_path),
        "reference_normalization": contract["normalization"],
        "mu_grid": contract["sweeps"][SWEEP_MU],
        "beta_grid": contract["sweeps"][SWEEP_BETA],
        "models": list(MODELS),
        "plotted_sorted_positions": PLOTTED_POSITIONS,
        "output_guard_position": OUTPUT_GUARD_POSITION,
        "summary": {
            "root_data": summaries,
            "all_export_ranges_passed": all_export_pass,
            "exact_mu_grid_passed": exact_mu_grid,
            "exact_beta_grid_passed": exact_beta_grid,
        },
        "closing_stage": {
            "new_points_executed": checkpoint["new_points_executed"],
            "reused_points": checkpoint["initial_reused_point_count"],
            "ready_points_recalculated": 0,
            "parallel_workers_used_in_closing_stage": 0,
            "inherited_workers_drained": checkpoint[
                "inherited_workers_drained"
            ],
            "inherited_points_allowed_to_finish": checkpoint[
                "inherited_points_allowed_to_finish"
            ],
            "global_restarts": 0,
            "thread_limits": CLOSING_THREAD_LIMITS,
        },
        "exact_qualifications_affecting_plotted_range": (
            plotted_range_qualifications
        ),
        "exact_qualifications_only_above_root9": tail_only_qualifications,
        "statuses": statuses,
        "frozen_model_hashes_before": (
            model_manifest["frozen_model_hashes_before"]
        ),
        "frozen_model_hashes_after": frozen_after,
        "frozen_models_preserved": (
            model_manifest["frozen_model_hashes_before"] == frozen_after
        ),
        "generated_file_hashes": {
            name: rlb2c.rlb2b._sha256(target / name)
            for name in generated_names
        },
        "figures_created": 2,
        "inter_model_relative_differences_computed": False,
        "delta_values_run": [DELTA],
        "case_ids_run": ["A"],
        "root_frequencies_cross_seeded_between_models": False,
        "pre_closing_execution_history": {
            "mu_sweep_new_rlb": {
                "used": True,
                "worker_slices_half_open": [
                    [0, 20],
                    [20, 27],
                    [27, 34],
                    [34, 37],
                    [37, 38],
                    [38, 41],
                ],
                "overlap": False,
                "physics_or_root_policy_changed": False,
                "resource_control": {
                    "unwritten_parallel_attempts_stopped": 3,
                    "minimum_observed_available_memory_MB": 107,
                    "partial_rows_written_by_stopped_attempts": 0,
                    "reason": "PREVENT_OUT_OF_MEMORY",
                },
            }
        },
        "production_physics_changed": False,
        "new_root_solver_created": False,
        "Ritz_run": False,
        "FEM_run": False,
        "branch_tracking_run": False,
        "MAC_run": False,
        "veering_analysis_run": False,
        "commit_performed": False,
        "push_performed": False,
    }
    _atomic_write_json(target / "run_manifest.json", run_manifest)
    return run_manifest


def run_all(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    mu_values: Sequence[float] | None = None,
    beta_values: Sequence[float] | None = None,
) -> dict[str, Any]:
    active_mu = (
        mu_grid() if mu_values is None else np.asarray(mu_values, dtype=float)
    )
    active_beta = (
        beta_grid()
        if beta_values is None
        else np.asarray(beta_values, dtype=float)
    )
    prepare_output(output_dir, active_mu, active_beta)
    for sweep in SWEEPS:
        for model in MODELS:
            solve_sweep_model(
                sweep,
                model,
                output_dir=output_dir,
                mu_values=active_mu,
                beta_values=active_beta,
            )
    return finalize_output(output_dir, active_mu, active_beta)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=(
            "prepare",
            "mu_eb",
            "mu_old_timoshenko",
            "mu_new_rlb",
            "beta_eb",
            "beta_old_timoshenko",
            "beta_new_rlb",
            "closing_audit",
            "closing_consolidate_rlb",
            "closing_mu_old_timoshenko",
            "closing_mu_new_rlb",
            "plot_beta",
            "plot_mu",
            "finalize",
            "all",
        ),
        default="all",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR
    )
    parser.add_argument("--mu-min", type=float, default=DEFAULT_MU_MIN)
    parser.add_argument("--mu-max", type=float, default=DEFAULT_MU_MAX)
    parser.add_argument("--mu-step", type=float, default=DEFAULT_MU_STEP)
    parser.add_argument(
        "--beta-min-deg", type=float, default=DEFAULT_BETA_MIN_DEG
    )
    parser.add_argument(
        "--beta-max-deg", type=float, default=DEFAULT_BETA_MAX_DEG
    )
    parser.add_argument(
        "--beta-step-deg", type=float, default=DEFAULT_BETA_STEP_DEG
    )
    parser.add_argument(
        "--point-start",
        type=int,
        default=0,
        help=(
            "Worker-only zero-based start index on the prepared full grid; "
            "use sequentially for resumable chunks."
        ),
    )
    parser.add_argument(
        "--point-stop",
        type=int,
        default=None,
        help="Worker-only exclusive stop index on the prepared full grid.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    active_mu = mu_grid(args.mu_min, args.mu_max, args.mu_step)
    active_beta = beta_grid(
        args.beta_min_deg, args.beta_max_deg, args.beta_step_deg
    )
    output_dir = Path(args.output_dir)
    worker_modes = {
        "mu_eb": (SWEEP_MU, MODEL_EB),
        "mu_old_timoshenko": (SWEEP_MU, MODEL_OLD),
        "mu_new_rlb": (SWEEP_MU, MODEL_RLB),
        "beta_eb": (SWEEP_BETA, MODEL_EB),
        "beta_old_timoshenko": (SWEEP_BETA, MODEL_OLD),
        "beta_new_rlb": (SWEEP_BETA, MODEL_RLB),
    }
    if args.mode == "prepare":
        prepare_output(output_dir, active_mu, active_beta)
    elif args.mode in worker_modes:
        sweep, model = worker_modes[args.mode]
        solve_sweep_model(
            sweep,
            model,
            output_dir=output_dir,
            mu_values=active_mu,
            beta_values=active_beta,
            point_start=args.point_start,
            point_stop=args.point_stop,
        )
    elif args.mode == "closing_audit":
        print(
            json.dumps(
                closing_mu_audit(output_dir, include_rlb_shards=True),
                indent=2,
                sort_keys=True,
            )
        )
    elif args.mode == "closing_consolidate_rlb":
        print(
            json.dumps(
                consolidate_existing_rlb_mu_rows(output_dir),
                indent=2,
                sort_keys=True,
            )
        )
    elif args.mode == "closing_mu_old_timoshenko":
        solve_one_missing_mu_point(
            MODEL_OLD, args.point_start, output_dir=output_dir
        )
    elif args.mode == "closing_mu_new_rlb":
        solve_one_missing_mu_point(
            MODEL_RLB, args.point_start, output_dir=output_dir
        )
    elif args.mode in {"plot_beta", "plot_mu"}:
        sweep = SWEEP_BETA if args.mode == "plot_beta" else SWEEP_MU
        rows = {
            model: rlb2c.rlb2b._read_csv(
                output_dir / ROOT_FILENAMES[(sweep, model)]
            )
            for model in MODELS
        }
        expected = active_beta if sweep == SWEEP_BETA else active_mu
        for model in MODELS:
            summary = validate_root_rows(sweep, model, rows[model], expected)
            if not summary["exact_row_structure_passed"]:
                raise RuntimeError(
                    f"Cannot plot incomplete {sweep}/{model}: {summary}."
                )
        create_plot(
            rows,
            output_dir / PLOT_FILENAMES[sweep],
            axis_field="beta_deg" if sweep == SWEEP_BETA else "mu",
            x_label="beta" if sweep == SWEEP_BETA else "mu",
            x_limits=(float(expected[0]), float(expected[-1])),
        )
    elif args.mode == "finalize":
        finalize_output(output_dir, active_mu, active_beta)
    else:
        run_all(output_dir, active_mu, active_beta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
