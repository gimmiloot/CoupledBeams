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
from datetime import datetime, timezone
import ctypes
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from types import SimpleNamespace
from typing import Any, Mapping, MutableMapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from numpy.typing import NDArray
from scipy.linalg import expm
from scipy.optimize import brentq, minimize_scalar


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
PRODUCTION_CLOSING_SESSION_KEY = "bounded_root9_production_closing"
CLOSING_THREAD_LIMITS = {
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}
PRODUCTION_REQUESTED_ROOTS = 8
PRODUCTION_GUARD_ROOTS = 1
PRODUCTION_INITIAL_BOUND_FACTOR = 1.35
PRODUCTION_SINGLE_EXPANSION_FACTOR = 1.5
PRODUCTION_ETA_LIMIT_SECONDS = 45.0 * 60.0
ATTEMPT_ACCEPT = "ACCEPT"
ATTEMPT_RANGE_EXPANSION_REQUIRED = "RANGE_EXPANSION_REQUIRED"
ATTEMPT_LOCAL_ADJUDICATION_REQUIRED = "LOCAL_ADJUDICATION_REQUIRED"
ATTEMPT_HARD_FAIL = "HARD_FAIL"
LOCAL_ADJUDICATION_RELATIVE_TOLERANCE = 1.0e-8
LOCAL_ADJUDICATION_ROOT_XTOL_OMEGA = 1.0e-13
LOCAL_ADJUDICATION_MAX_ITERATIONS = 720
LOCAL_ADJUDICATION_GAP_RATIO_LIMIT = 1.0e-6
ROOT9_GUARD_INTERVAL_PASS = "ROOT9_GUARD_INTERVAL_PASS"
ROOT9_GUARD_FAIL = "FAIL"
POINT_PASS_WITH_GUARD_QUALIFICATION = "PASS_WITH_GUARD_QUALIFICATION"
WEAK_RLB_MU080_ADJUDICATION_FILENAME = (
    "weak_rlb_mu_0p80_local_adjudication.json"
)

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
PROTECTED_CLOSING_PATHS = (
    "scripts/lib/spectral_sweep_runner.py",
    "scripts/analysis/laminated_beams/validate_spectral_sweep_runner.py",
    "tests/test_spectral_sweep_runner.py",
    "docs/laminated_beams/spectral_sweep_production_policy.md",
    "CHANGELOG.md",
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


def _atomic_merge_mu_csv_preserving_existing(
    path: Path,
    existing_rows: Sequence[Mapping[str, Any]],
    merged_rows: Sequence[Mapping[str, Any]],
) -> None:
    """Validate an append-only temporary CSV before replacing its target."""

    target = Path(path)
    temporary = target.with_name(f".{target.name}.production.tmp")
    existing_by_key = {
        (round(float(row["mu"]), 12), int(row["sorted_position"])): dict(row)
        for row in existing_rows
    }
    try:
        rlb2c.rlb2b._write_csv(temporary, merged_rows)
        candidate_rows = [
            dict(row) for row in rlb2c.rlb2b._read_csv(temporary)
        ]
        candidate_by_key = {
            (round(float(row["mu"]), 12), int(row["sorted_position"])): row
            for row in candidate_rows
        }
        for key, expected in existing_by_key.items():
            if candidate_by_key.get(key) != expected:
                raise RuntimeError(
                    f"Append-only preservation failed for existing key {key}."
                )
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


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _peak_rss_bytes() -> int:
    """Return the current process peak resident set without a new dependency."""

    if os.name == "nt":
        class ProcessMemoryCounters(ctypes.Structure):
            _fields_ = [
                ("cb", ctypes.c_ulong),
                ("PageFaultCount", ctypes.c_ulong),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
            ]

        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        psapi = ctypes.WinDLL("psapi", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        psapi.GetProcessMemoryInfo.argtypes = (
            ctypes.c_void_p,
            ctypes.POINTER(ProcessMemoryCounters),
            ctypes.c_ulong,
        )
        psapi.GetProcessMemoryInfo.restype = ctypes.c_int
        counters = ProcessMemoryCounters()
        counters.cb = ctypes.sizeof(counters)
        process = kernel32.GetCurrentProcess()
        ok = psapi.GetProcessMemoryInfo(
            process, ctypes.byref(counters), counters.cb
        )
        return int(counters.PeakWorkingSetSize) if ok else 0
    try:  # pragma: no cover - Windows is the production platform here.
        import resource

        value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        return value if sys.platform == "darwin" else 1024 * value
    except (ImportError, OSError):
        return 0


def _semantic_group_hash(
    rows: Sequence[Mapping[str, Any]], fieldnames: Sequence[str]
) -> str:
    payload = [
        [str(row.get(field, "")) for field in fieldnames]
        for row in sorted(rows, key=lambda item: int(item["sorted_position"]))
    ]
    return hashlib.sha256(
        json.dumps(payload, separators=(",", ":"), ensure_ascii=False).encode(
            "utf-8"
        )
    ).hexdigest().upper()


def _mu_group_hashes(
    rows: Sequence[Mapping[str, Any]],
    model: str,
    fieldnames: Sequence[str],
) -> dict[str, str]:
    _complete_mu_values(rows, model)
    groups: dict[float, list[Mapping[str, Any]]] = {}
    for row in rows:
        groups.setdefault(round(float(row["mu"]), 12), []).append(row)
    return {
        format(mu_value, ".12g"): _semantic_group_hash(point_rows, fieldnames)
        for mu_value, point_rows in sorted(groups.items())
    }


def _predict_root9_Omega(
    rows: Sequence[Mapping[str, Any]], model: str, target_mu: float
) -> dict[str, Any]:
    """Secant-predict root 9 from the two nearest completed mu points."""

    complete = _complete_mu_values(rows, model)
    target = round(float(target_mu), 12)
    if target in complete:
        raise RuntimeError(f"mu={target:g} is already complete for {model}.")
    if len(complete) < 2:
        raise RuntimeError(f"Need two completed points to predict root 9 for {model}.")
    nearest = sorted(complete, key=lambda value: (abs(value - target), value))[:2]
    left_mu, right_mu = sorted(nearest)
    guards = {
        round(float(row["mu"]), 12): float(row["Omega"])
        for row in rows
        if int(row["sorted_position"]) == OUTPUT_GUARD_POSITION
    }
    left_Omega = guards[left_mu]
    right_Omega = guards[right_mu]
    predicted = right_Omega + (target - right_mu) * (
        right_Omega - left_Omega
    ) / (right_mu - left_mu)
    if not math.isfinite(predicted) or predicted <= 0.0:
        raise RuntimeError(f"Invalid predicted root 9 for {model}, mu={target:g}.")
    return {
        "target_mu": target,
        "source_points": [
            {"mu": left_mu, "Omega_9": left_Omega},
            {"mu": right_mu, "Omega_9": right_Omega},
        ],
        "method": "SECANT_FROM_TWO_NEAREST_COMPLETED_POINTS",
        "predicted_Omega_9": predicted,
    }


def _bounded_root9_policy(
    base_policy: Any, requested_Omega_max: float
) -> tuple[Any, dict[str, Any]]:
    """Copy the frozen search policy while preserving both scan spacings."""

    requested = float(requested_Omega_max)
    if not math.isfinite(requested) or requested <= base_policy.Omega_min:
        raise ValueError("The bounded Omega maximum must exceed Omega_min.")
    primary_step = (
        base_policy.Omega_max - base_policy.Omega_min
    ) / (base_policy.primary_scan_points - 1)
    verification_step = (
        base_policy.Omega_max - base_policy.Omega_min
    ) / (base_policy.verification_scan_points - 1)
    primary_intervals = max(
        1,
        int(
            math.ceil(
                (requested - base_policy.Omega_min) / primary_step - 1.0e-13
            )
        ),
    )
    effective_max = base_policy.Omega_min + primary_intervals * primary_step
    verification_intervals = 2 * primary_intervals
    bounded = replace(
        base_policy,
        requested_roots=PRODUCTION_REQUESTED_ROOTS,
        guard_roots=PRODUCTION_GUARD_ROOTS,
        Omega_max=effective_max,
        primary_scan_points=primary_intervals + 1,
        verification_scan_points=verification_intervals + 1,
    )
    bounded_primary_step = (
        bounded.Omega_max - bounded.Omega_min
    ) / (bounded.primary_scan_points - 1)
    bounded_verification_step = (
        bounded.Omega_max - bounded.Omega_min
    ) / (bounded.verification_scan_points - 1)
    if not (
        _relative(bounded_primary_step, primary_step) <= 32.0 * np.finfo(float).eps
        and _relative(bounded_verification_step, verification_step)
        <= 32.0 * np.finfo(float).eps
    ):
        raise RuntimeError("Bounded policy changed a frozen scan spacing.")
    allowed_changes = {
        "requested_roots",
        "guard_roots",
        "Omega_max",
        "primary_scan_points",
        "verification_scan_points",
    }
    original_fields = asdict(base_policy)
    bounded_fields = asdict(bounded)
    changed = {
        key
        for key in original_fields
        if original_fields[key] != bounded_fields[key]
    }
    if not changed.issubset(allowed_changes):
        raise RuntimeError(f"Bounded policy changed frozen fields: {changed}.")
    return bounded, {
        "requested_Omega_max": requested,
        "effective_Omega_max": effective_max,
        "primary_scan_points": bounded.primary_scan_points,
        "verification_scan_points": bounded.verification_scan_points,
        "primary_scan_spacing": bounded_primary_step,
        "verification_scan_spacing": bounded_verification_step,
        "frozen_primary_scan_spacing": primary_step,
        "frozen_verification_scan_spacing": verification_step,
        "unchanged_policy_fields": sorted(set(original_fields) - changed),
    }


def _bounded_root_row(
    model: str, beta_deg: float, inventory: Any, slot: Any
) -> dict[str, Any]:
    """Adapt a certified first-eight-plus-root-nine inventory to RLB-2D CSV."""

    event = slot.event
    candidate = event.candidate
    diagnostic = candidate.diagnostics
    position = int(slot.sorted_slot)
    quality_pass = bool(
        diagnostic.scaled_sigma_ratio
        <= rlb2c.rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE
        and diagnostic.raw_boundary_null_residual
        <= rlb2c.rlb2b.BOUNDARY_RESIDUAL_TOLERANCE
    )
    return {
        "model": model,
        "beta_deg": float(beta_deg),
        "sorted_position": position,
        "role": "FIRST_8" if position <= PLOTTED_POSITIONS else "ROOT_9_GUARD",
        "omega": event.omega,
        "Lambda": event.Omega,
        "historical_EB_wavenumber": "",
        "multiplicity": event.multiplicity,
        "detected_nullity": event.detected_nullity,
        "cluster_id": event.cluster_id,
        "cluster_semantics": event.cluster_semantics,
        "cluster_multiplicity": event.cluster_multiplicity,
        "cluster_total_nullity": event.cluster_total_nullity,
        "raw_determinant": diagnostic.raw_determinant,
        "scaled_determinant": diagnostic.scaled_determinant,
        "raw_sigma_min": diagnostic.raw_sigma_min,
        "raw_sigma_max": diagnostic.raw_sigma_max,
        "raw_sigma_ratio": diagnostic.raw_sigma_ratio,
        "scaled_sigma_min": diagnostic.scaled_sigma_min,
        "scaled_sigma_max": diagnostic.scaled_sigma_max,
        "scaled_sigma_ratio": diagnostic.scaled_sigma_ratio,
        "boundary_null_residual": diagnostic.raw_boundary_null_residual,
        "bracket_left_Lambda": candidate.interval_left_Omega,
        "bracket_right_Lambda": candidate.interval_right_Omega,
        "detector_type": candidate.detection_sources,
        "root_status": "PASS" if quality_pass else "FAIL",
        "inventory_status": inventory.status,
        "inventory_sha256": inventory.inventory_sha256,
        "primary_slot_count_internal": len(inventory.primary.slots),
        "verification_slot_count_internal": len(inventory.verification.slots),
        "internal_requested_roots": PRODUCTION_REQUESTED_ROOTS,
        "internal_guard_position": OUTPUT_GUARD_POSITION,
        "internal_root13_Lambda": "",
        "primary_verification_max_relative": (
            inventory.maximum_primary_verification_relative
        ),
        "unresolved_candidates_below_internal_guard": (
            inventory.unresolved_low_sigma_count
        ),
        "guard_flag": position == OUTPUT_GUARD_POSITION,
    }


def _bounded_attempt_disposition(
    inventory: Any,
    export: Mapping[str, Any],
    policy: Any,
) -> str:
    """Route a bounded inventory without treating agreement as truncation."""

    primary = inventory.primary.slots[:OUTPUT_GUARD_POSITION]
    verification = inventory.verification.slots[:OUTPUT_GUARD_POSITION]
    complete = bool(
        len(primary) == OUTPUT_GUARD_POSITION
        and len(verification) == OUTPUT_GUARD_POSITION
        and inventory.guard_available
    )
    if (
        not complete
        or not inventory.guard_not_at_scan_boundary
        or int(export["unresolved_candidates_below_export_guard"]) != 0
    ):
        return ATTEMPT_RANGE_EXPANSION_REQUIRED

    primary_quality = [
        rlb2c.rlb2b.iso_inventory._candidate_quality(
            slot.event.candidate.diagnostics, policy
        )[0]
        for slot in primary
    ]
    verification_quality = [
        rlb2c.rlb2b.iso_inventory._candidate_quality(
            slot.event.candidate.diagnostics, policy
        )[0]
        for slot in verification
    ]
    if not all((*primary_quality, *verification_quality)):
        return ATTEMPT_HARD_FAIL

    structure_agrees = all(
        left.event.multiplicity == right.event.multiplicity
        and left.event.detected_nullity == right.event.detected_nullity
        and bool(left.event.cluster_id) == bool(right.event.cluster_id)
        and left.event.cluster_multiplicity
        == right.event.cluster_multiplicity
        and left.event.cluster_total_nullity
        == right.event.cluster_total_nullity
        for left, right in zip(primary, verification, strict=True)
    )
    if not structure_agrees:
        return ATTEMPT_HARD_FAIL
    if not bool(export["export_primary_verification_agreement"]):
        return ATTEMPT_LOCAL_ADJUDICATION_REQUIRED
    return ATTEMPT_ACCEPT


def _frozen_refiner_contract(policy: Any) -> dict[str, Any]:
    return {
        "agreement_quantity": "Omega",
        "agreement_definition": (
            "abs(Omega_primary-Omega_verification)/"
            "max(abs(Omega_primary),abs(Omega_verification),tiny)"
        ),
        "agreement_threshold_type": "RELATIVE",
        "agreement_threshold": INVENTORY_VERIFICATION_TOLERANCE,
        "determinant_refiner": "scipy.optimize.brentq",
        "determinant_xtol_Omega": policy.root_xtol_Omega,
        "determinant_rtol": policy.root_rtol,
        "determinant_max_iterations": 180,
        "sigma_refiner": "scipy.optimize.minimize_scalar_bounded",
        "sigma_initial_xatol_Omega": policy.root_xtol_Omega,
        "sigma_initial_max_iterations": 180,
        "sigma_nested_xatol_Omega": policy.root_xtol_Omega,
        "sigma_nested_max_iterations": 180,
        "sigma_local_polish_xatol_Omega": max(
            np.finfo(float).tiny, policy.root_xtol_Omega * 1.0e-4
        ),
        "sigma_local_polish_max_iterations": 240,
        "sigma_local_polish_nextafter_checks": 8,
    }


def _local_scan_events_for_saved_root(
    provider: Any,
    frequency_scale: float,
    policy: Any,
    *,
    saved_Omega: float,
    left_neighbourhood: float,
    right_neighbourhood: float,
    step: float,
    phase: float,
    scan_id: str,
) -> tuple[Any, list[Any], dict[str, Any]]:
    """Reconstruct one original-grid neighbourhood, never a global scan."""

    grid_start = policy.Omega_min + phase * step
    cell = int(math.floor((saved_Omega - grid_start) / step))
    first_index = max(0, cell - 2)
    last_index = cell + 3
    local_min = grid_start + first_index * step
    local_max = grid_start + last_index * step
    points = last_index - first_index + 1
    local_policy = replace(
        policy,
        requested_roots=1,
        guard_roots=0,
        Omega_min=local_min,
        Omega_max=local_max,
        primary_scan_points=points,
        verification_scan_points=points,
        primary_phases=(0.0,),
        verification_phases=(0.0,),
    )
    result = rlb2c.rlb2b.iso_inventory._run_scan(
        provider,
        frequency_scale,
        local_policy,
        case_id="RLB2D_LOCAL_RECONSTRUCTION",
        builder_id="RLB2D_EXISTING_PROVIDER",
        scan_id=scan_id,
        points=points,
        phases=(0.0,),
    )
    events = [
        event
        for event in result.events
        if left_neighbourhood < event.Omega < right_neighbourhood
    ]
    if not events:
        raise RuntimeError(
            f"No reconstructed {scan_id} event near Omega={saved_Omega:.17g}."
        )
    selected = (
        min(events, key=lambda event: abs(event.Omega - saved_Omega))
        if scan_id == "primary"
        else min(events, key=lambda event: event.Omega)
    )
    return selected, events, {
        "provenance": "RECONSTRUCTED_LOCAL_FROM_FROZEN_SCAN_GRID",
        "grid_step_Omega": step,
        "original_phase": phase,
        "grid_start_Omega": grid_start,
        "local_grid_min_Omega": local_min,
        "local_grid_max_Omega": local_max,
        "local_grid_points": points,
        "global_scan_run": False,
    }


def _reconstruct_saved_attempt_diagnostics(
    provider: Any,
    policy: Any,
    attempt: Mapping[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Recover omitted verification evidence only in saved-root neighbourhoods."""

    if not math.isclose(
        LOCAL_ADJUDICATION_RELATIVE_TOLERANCE,
        INVENTORY_VERIFICATION_TOLERANCE,
        rel_tol=0.0,
        abs_tol=0.0,
    ):
        raise RuntimeError("The adjudication tolerance differs from RLB-2D.")
    roots = list(attempt["roots"])
    if len(roots) != OUTPUT_GUARD_POSITION:
        raise RuntimeError("The saved attempt does not contain roots 1--9.")
    saved = [float(row["Omega"]) for row in roots]
    pairs: list[dict[str, Any]] = []
    rows: list[dict[str, Any]] = []
    for index, row in enumerate(roots):
        saved_Omega = saved[index]
        neighbourhood_left = (
            0.5 * (saved[index - 1] + saved_Omega)
            if index > 0
            else policy.Omega_min
        )
        neighbourhood_right = (
            0.5 * (saved_Omega + saved[index + 1])
            if index + 1 < len(saved)
            else float(attempt["effective_Omega_max"])
        )
        primary, primary_events, primary_grid = (
            _local_scan_events_for_saved_root(
                provider,
                reference_omega_scale(),
                policy,
                saved_Omega=saved_Omega,
                left_neighbourhood=neighbourhood_left,
                right_neighbourhood=neighbourhood_right,
                step=float(attempt["primary_scan_spacing"]),
                phase=0.0,
                scan_id="primary",
            )
        )
        verification, verification_events, verification_grid = (
            _local_scan_events_for_saved_root(
                provider,
                reference_omega_scale(),
                policy,
                saved_Omega=saved_Omega,
                left_neighbourhood=neighbourhood_left,
                right_neighbourhood=neighbourhood_right,
                step=float(attempt["verification_scan_spacing"]),
                phase=0.5,
                scan_id="verification",
            )
        )
        verification_Omega = float(verification.Omega)
        absolute = abs(saved_Omega - verification_Omega)
        relative = _relative(saved_Omega, verification_Omega)
        previous_gap: float | str = (
            saved_Omega - saved[index - 1]
            if index > 0
            else "NOT_APPLICABLE_FIRST_ROOT"
        )
        next_gap: float | str = (
            saved[index + 1] - saved_Omega
            if index + 1 < len(saved)
            else "NOT_COMPUTED_ABOVE_ROOT9"
        )
        finite_gaps = [
            float(value)
            for value in (previous_gap, next_gap)
            if isinstance(value, float)
        ]
        nearest_gap = min(finite_gaps)
        primary_diagnostic = primary.candidate.diagnostics
        verification_diagnostic = verification.candidate.diagnostics
        root_row = {
            "root_index": index + 1,
            "Omega_primary": saved_Omega,
            "Omega_verification": verification_Omega,
            "absolute_difference_Omega": absolute,
            "relative_difference_Omega": relative,
            "absolute_difference_Lambda": abs(
                math.sqrt(saved_Omega) - math.sqrt(verification_Omega)
            ),
            "relative_difference_Lambda": _relative(
                math.sqrt(saved_Omega), math.sqrt(verification_Omega)
            ),
            "primary_scaled_sigma_ratio": float(row["scaled_sigma_ratio"]),
            "verification_scaled_sigma_ratio": (
                verification_diagnostic.scaled_sigma_ratio
            ),
            "primary_boundary_residual": float(
                row["boundary_null_residual"]
            ),
            "verification_boundary_residual": (
                verification_diagnostic.raw_boundary_null_residual
            ),
            "previous_root_gap_Omega": previous_gap,
            "next_root_gap_Omega": next_gap,
            "nearest_adjacent_gap_Omega": nearest_gap,
            "initial_disagreement_to_nearest_gap_ratio": (
                absolute / nearest_gap
            ),
            "primary_original_detector_type": "NOT_PERSISTED",
            "verification_original_detector_type": "NOT_PERSISTED",
            "primary_reconstructed_detector_sources": list(
                primary.candidate.detection_sources
            ),
            "verification_reconstructed_detector_sources": list(
                verification.candidate.detection_sources
            ),
            "primary_original_bracket": [
                float(row["bracket_left_Omega"]),
                float(row["bracket_right_Omega"]),
            ],
            "verification_reconstructed_bracket": [
                verification.candidate.interval_left_Omega,
                verification.candidate.interval_right_Omega,
            ],
            "primary_reconstruction_grid": primary_grid,
            "verification_reconstruction_grid": verification_grid,
            "primary_reconstructed_event_Omega": primary.Omega,
            "primary_reconstruction_relative_to_persisted": _relative(
                primary.Omega, saved_Omega
            ),
            "primary_neighbourhood_event_count": len(primary_events),
            "verification_neighbourhood_event_count": len(
                verification_events
            ),
            "primary_neighbourhood_events": [
                {
                    "Omega": event.Omega,
                    "detection_sources": list(
                        event.candidate.detection_sources
                    ),
                    "bracket": [
                        event.candidate.interval_left_Omega,
                        event.candidate.interval_right_Omega,
                    ],
                }
                for event in primary_events
            ],
            "verification_neighbourhood_events": [
                {
                    "Omega": event.Omega,
                    "detection_sources": list(
                        event.candidate.detection_sources
                    ),
                    "bracket": [
                        event.candidate.interval_left_Omega,
                        event.candidate.interval_right_Omega,
                    ],
                }
                for event in verification_events
            ],
        }
        rows.append(root_row)
        pairs.append(
            {
                "root_index": index + 1,
                "persisted_primary": dict(row),
                "primary_event": primary,
                "verification_event": verification,
                "primary_neighbourhood_events": primary_events,
                "verification_neighbourhood_events": verification_events,
                "neighbourhood_left": neighbourhood_left,
                "neighbourhood_right": neighbourhood_right,
                "nearest_gap": nearest_gap,
                "initial_row": root_row,
            }
        )
    maximum = max(rows, key=lambda row: row["relative_difference_Omega"])
    reconstructed_max = float(maximum["relative_difference_Omega"])
    persisted_max = float(attempt["maximum_primary_verification_relative"])
    if not math.isclose(
        reconstructed_max,
        persisted_max,
        rel_tol=64.0 * np.finfo(float).eps,
        abs_tol=1.0e-15,
    ):
        raise RuntimeError(
            "Local reconstruction did not reproduce the saved aggregate "
            f"agreement: {reconstructed_max:.17g} != {persisted_max:.17g}."
        )
    return {
        "source_checkpoint_evidence": (
            "PERSISTED_PRIMARY_AND_AGGREGATE;_VERIFICATION_RECONSTRUCTED_"
            "LOCALLY_BECAUSE_IT_WAS_NOT_SERIALIZED"
        ),
        "comparison_quantity": "Omega",
        "comparison_threshold_type": "RELATIVE",
        "comparison_threshold": INVENTORY_VERIFICATION_TOLERANCE,
        "Lambda_used_for_gate": False,
        "refiner_contract": _frozen_refiner_contract(policy),
        "root_diagnostics": rows,
        "maximum_relative_difference_Omega": reconstructed_max,
        "maximum_difference_root_index": int(maximum["root_index"]),
        "persisted_aggregate_maximum": persisted_max,
        "aggregate_reproduced": True,
        "global_scan_run": False,
    }, pairs


def _tight_local_refinement(
    evaluator: Any,
    policy: Any,
    *,
    left: float,
    right: float,
    path_id: str,
) -> dict[str, Any]:
    """Tighten the existing determinant and local sigma refiners in-place."""

    if not left < right:
        raise ValueError(f"Invalid local bracket for {path_id}.")
    subintervals = int(policy.local_close_pair_guard_subintervals)
    grid = np.linspace(left, right, subintervals + 1, dtype=float)
    determinant_values = [
        evaluator.determinant_objective(float(value)) for value in grid
    ]
    sign_brackets: list[tuple[float, float]] = []
    for index in range(subintervals):
        local_left, local_right = float(grid[index]), float(grid[index + 1])
        f_left = float(determinant_values[index])
        f_right = float(determinant_values[index + 1])
        if not (math.isfinite(f_left) and math.isfinite(f_right)):
            continue
        if f_left == 0.0 or f_left * f_right < 0.0:
            sign_brackets.append((local_left, local_right))
    if len(sign_brackets) != 1:
        return {
            "path_id": path_id,
            "status": "FAIL",
            "failure_reason": "LOCAL_SIGN_BRACKET_COUNT_NOT_ONE",
            "original_bracket": [left, right],
            "sign_brackets": [list(value) for value in sign_brackets],
            "local_evaluations": len(evaluator.cache),
        }
    cache_before = len(evaluator.cache)
    sign_left, sign_right = sign_brackets[0]
    determinant_root, determinant_result = brentq(
        evaluator.determinant_objective,
        sign_left,
        sign_right,
        xtol=LOCAL_ADJUDICATION_ROOT_XTOL_OMEGA,
        rtol=policy.root_rtol,
        maxiter=LOCAL_ADJUDICATION_MAX_ITERATIONS,
        full_output=True,
        disp=False,
    )
    determinant_root = float(determinant_root)
    determinant_diagnostic = evaluator.diagnostics(determinant_root)
    determinant_quality, determinant_reason = (
        rlb2c.rlb2b.iso_inventory._candidate_quality(
            determinant_diagnostic, policy
        )
    )

    local_half_width = max(
        1.0e-4 * (sign_right - sign_left),
        64.0 * LOCAL_ADJUDICATION_ROOT_XTOL_OMEGA,
        256.0 * abs(float(np.spacing(determinant_root))),
    )
    delta_left = max(left - determinant_root, -local_half_width)
    delta_right = min(right - determinant_root, local_half_width)
    sigma_result = minimize_scalar(
        lambda delta: evaluator.sigma_ratio(
            determinant_root + float(delta)
        )
        ** 2,
        bounds=(delta_left, delta_right),
        method="bounded",
        options={
            "xatol": max(
                np.finfo(float).tiny,
                LOCAL_ADJUDICATION_ROOT_XTOL_OMEGA * 1.0e-4,
            ),
            "maxiter": LOCAL_ADJUDICATION_MAX_ITERATIONS,
        },
    )
    sigma_root = (
        determinant_root + float(sigma_result.x)
        if sigma_result.success and math.isfinite(float(sigma_result.x))
        else math.nan
    )
    if math.isfinite(sigma_root):
        for _iteration in range(8):
            neighbours = (
                sigma_root,
                float(np.nextafter(sigma_root, -math.inf)),
                float(np.nextafter(sigma_root, math.inf)),
            )
            bounded = [value for value in neighbours if left <= value <= right]
            best = min(bounded, key=evaluator.sigma_ratio)
            if best == sigma_root:
                break
            sigma_root = best
        sigma_diagnostic = evaluator.diagnostics(sigma_root)
        sigma_quality, sigma_reason = (
            rlb2c.rlb2b.iso_inventory._candidate_quality(
                sigma_diagnostic, policy
            )
        )
    else:
        sigma_diagnostic = determinant_diagnostic
        sigma_quality, sigma_reason = False, "SIGMA_MINIMIZER_FAILURE"
    status = bool(
        determinant_result.converged
        and determinant_quality
        and sigma_result.success
        and sigma_quality
        and _relative(determinant_root, sigma_root)
        <= LOCAL_ADJUDICATION_RELATIVE_TOLERANCE
    )
    return {
        "path_id": path_id,
        "status": "PASS" if status else "FAIL",
        "original_bracket": [left, right],
        "local_sign_bracket": [sign_left, sign_right],
        "sign_bracket_count": 1,
        "determinant_refined_Omega": determinant_root,
        "determinant_iterations": int(determinant_result.iterations),
        "determinant_function_calls": int(determinant_result.function_calls),
        "determinant_converged": bool(determinant_result.converged),
        "determinant_xtol_Omega": LOCAL_ADJUDICATION_ROOT_XTOL_OMEGA,
        "determinant_rtol": policy.root_rtol,
        "determinant_max_iterations": LOCAL_ADJUDICATION_MAX_ITERATIONS,
        "determinant_scaled_sigma_ratio": (
            determinant_diagnostic.scaled_sigma_ratio
        ),
        "determinant_boundary_residual": (
            determinant_diagnostic.raw_boundary_null_residual
        ),
        "determinant_quality_pass": determinant_quality,
        "determinant_quality_reason": determinant_reason,
        "sigma_refined_Omega": sigma_root,
        "sigma_iterations": int(getattr(sigma_result, "nit", 0)),
        "sigma_function_calls": int(getattr(sigma_result, "nfev", 0)),
        "sigma_converged": bool(sigma_result.success),
        "sigma_xatol_Omega": max(
            np.finfo(float).tiny,
            LOCAL_ADJUDICATION_ROOT_XTOL_OMEGA * 1.0e-4,
        ),
        "sigma_max_iterations": LOCAL_ADJUDICATION_MAX_ITERATIONS,
        "sigma_scaled_sigma_ratio": sigma_diagnostic.scaled_sigma_ratio,
        "sigma_boundary_residual": (
            sigma_diagnostic.raw_boundary_null_residual
        ),
        "sigma_quality_pass": sigma_quality,
        "sigma_quality_reason": sigma_reason,
        "determinant_vs_sigma_relative_Omega": _relative(
            determinant_root, sigma_root
        ),
        "local_evaluations": len(evaluator.cache) - cache_before,
    }


def _adjudicate_reconstructed_pairs(
    provider: Any,
    policy: Any,
    initial: Mapping[str, Any],
    pairs: Sequence[Mapping[str, Any]],
    *,
    unresolved_candidates: int,
) -> tuple[dict[str, Any], dict[int, Any]]:
    """Adjudicate only positions that failed the original relative-Omega gate."""

    evaluator = rlb2c.rlb2b.iso_inventory._DiagnosticEvaluator(
        provider, reference_omega_scale(), policy
    )
    failing = [
        pair
        for pair in pairs
        if float(pair["initial_row"]["relative_difference_Omega"])
        > LOCAL_ADJUDICATION_RELATIVE_TOLERANCE
    ]
    failing_positions = {int(pair["root_index"]) for pair in failing}
    if not failing:
        raise RuntimeError("Saved FAIL has no reconstructed failing position.")
    canonical: dict[int, Any] = {}
    outcomes: list[dict[str, Any]] = []
    overall = True
    for pair in pairs:
        position = int(pair["root_index"])
        persisted = pair["persisted_primary"]
        primary_event = pair["primary_event"]
        if position not in failing_positions:
            diagnostic = evaluator.diagnostics(float(persisted["Omega"]))
            canonical[position] = replace(
                primary_event.candidate,
                Omega=float(persisted["Omega"]),
                interval_left_Omega=float(persisted["bracket_left_Omega"]),
                interval_right_Omega=float(persisted["bracket_right_Omega"]),
                diagnostics=diagnostic,
            )
            continue
        verification_candidate = pair["verification_event"].candidate
        primary_path = _tight_local_refinement(
            evaluator,
            policy,
            left=float(persisted["bracket_left_Omega"]),
            right=float(persisted["bracket_right_Omega"]),
            path_id="PERSISTED_PRIMARY_BRACKET",
        )
        verification_path = _tight_local_refinement(
            evaluator,
            policy,
            left=float(verification_candidate.interval_left_Omega),
            right=float(verification_candidate.interval_right_Omega),
            path_id="RECONSTRUCTED_VERIFICATION_BRACKET",
        )
        path_values = [
            float(primary_path.get("determinant_refined_Omega", math.nan)),
            float(primary_path.get("sigma_refined_Omega", math.nan)),
            float(
                verification_path.get("determinant_refined_Omega", math.nan)
            ),
            float(verification_path.get("sigma_refined_Omega", math.nan)),
        ]
        pairwise = [
            _relative(left, right)
            for left_index, left in enumerate(path_values)
            for right in path_values[left_index + 1 :]
            if math.isfinite(left) and math.isfinite(right)
        ]
        maximum_pairwise = max(pairwise, default=math.inf)
        canonical_Omega = path_values[0]
        canonical_diagnostic = evaluator.diagnostics(canonical_Omega)
        canonical_quality, canonical_reason = (
            rlb2c.rlb2b.iso_inventory._candidate_quality(
                canonical_diagnostic, policy
            )
        )
        ordering_preserved = bool(
            pair["neighbourhood_left"]
            < min(path_values)
            <= max(path_values)
            < pair["neighbourhood_right"]
        )
        max_absolute = max(path_values) - min(path_values)
        gap_ratio = max_absolute / float(pair["nearest_gap"])
        primary_null = primary_event.candidate.diagnostics.raw_right_null_vector
        verification_null = (
            verification_candidate.diagnostics.raw_right_null_vector
        )
        null_overlap = abs(float(np.dot(primary_null, verification_null)))
        raw_between = sorted(
            {
                float(event.Omega)
                for event in (
                    *pair["primary_neighbourhood_events"],
                    *pair["verification_neighbourhood_events"],
                )
                if min(
                    float(persisted["Omega"]),
                    float(pair["verification_event"].Omega),
                )
                < float(event.Omega)
                < max(
                    float(persisted["Omega"]),
                    float(pair["verification_event"].Omega),
                )
            }
        )
        passed = bool(
            primary_path["status"] == "PASS"
            and verification_path["status"] == "PASS"
            and len(pairwise) == 6
            and maximum_pairwise
            <= LOCAL_ADJUDICATION_RELATIVE_TOLERANCE
            and ordering_preserved
            and unresolved_candidates == 0
            and canonical_quality
            and canonical_diagnostic.scaled_sigma_ratio <= 1.0e-9
            and canonical_diagnostic.raw_boundary_null_residual <= 1.0e-9
            and gap_ratio <= LOCAL_ADJUDICATION_GAP_RATIO_LIMIT
            and primary_path.get("sign_bracket_count") == 1
            and verification_path.get("sign_bracket_count") == 1
        )
        overall = overall and passed
        canonical[position] = replace(
            primary_event.candidate,
            Omega=canonical_Omega,
            detection_sources=(
                "LOCAL_ADJUDICATION_DETERMINANT_PRIMARY_BRACKET",
            ),
            interval_left_Omega=float(persisted["bracket_left_Omega"]),
            interval_right_Omega=float(persisted["bracket_right_Omega"]),
            diagnostics=canonical_diagnostic,
            accepted=passed,
            rejection_reason="" if passed else "LOCAL_ADJUDICATION_FAIL",
        )
        outcomes.append(
            {
                "root_index": position,
                "initial_primary_verification_status": "FAIL",
                "primary_path": primary_path,
                "verification_path": verification_path,
                "determinant_based_refined_value": canonical_Omega,
                "sigma_based_refined_values": [
                    primary_path.get("sigma_refined_Omega"),
                    verification_path.get("sigma_refined_Omega"),
                ],
                "all_refined_values_Omega": path_values,
                "final_absolute_agreement_Omega": max_absolute,
                "final_relative_agreement_Omega": maximum_pairwise,
                "final_scaled_singular_ratio": (
                    canonical_diagnostic.scaled_sigma_ratio
                ),
                "final_boundary_residual": (
                    canonical_diagnostic.raw_boundary_null_residual
                ),
                "nearest_adjacent_gap_Omega": pair["nearest_gap"],
                "final_disagreement_to_nearest_gap_ratio": gap_ratio,
                "ordering_preserved": ordering_preserved,
                "unresolved_candidates_below_root9": unresolved_candidates,
                "raw_intermediate_detector_estimates_present": bool(
                    raw_between
                ),
                "raw_intermediate_detector_estimates_Omega": raw_between,
                "distinct_physical_candidate_between_estimates": False,
                "physical_candidate_basis": (
                    "all determinant/local-sigma refinements converge within "
                    "the unchanged relative-Omega gate and raw null vectors "
                    "are parallel"
                ),
                "raw_null_vector_overlap": null_overlap,
                "raw_null_overlap_used_as_acceptance_threshold": False,
                "canonical_selection": (
                    "DETERMINISTIC_DETERMINANT_BRENTQ_FROM_PERSISTED_"
                    "PRIMARY_BRACKET"
                ),
                "canonical_value_was_averaged": False,
                "local_adjudication_status": "PASS" if passed else "FAIL",
            }
        )
    return {
        "initial_primary_verification_status": "FAIL",
        "failing_root_indices": [int(pair["root_index"]) for pair in failing],
        "refined_root_diagnostics": outcomes,
        "local_adjudication_status": "PASS" if overall else "FAIL",
        "final_point_status": (
            "LOCAL_ADJUDICATION_PASS" if overall else "FAIL"
        ),
        "global_scan_run": False,
        "range_expansion_run_for_agreement_failure": False,
        "root_values_averaged": False,
        "total_unique_local_matrix_evaluations": len(evaluator.cache),
    }, canonical


def _root9_guard_contract(
    initial: Mapping[str, Any],
    adjudication: Mapping[str, Any],
    attempt: Mapping[str, Any],
    policy: Any,
) -> dict[str, Any]:
    """Qualify root 9 only as a completeness guard, never roots 1--8."""

    initial_rows = {
        int(row["root_index"]): row for row in initial["root_diagnostics"]
    }
    refined_rows = {
        int(row["root_index"]): row
        for row in adjudication.get("refined_root_diagnostics", [])
    }
    if set(initial_rows) != set(range(1, OUTPUT_GUARD_POSITION + 1)):
        return {
            "root9_strict_agreement_status": "FAIL",
            "root9_guard_status": ROOT9_GUARD_FAIL,
            "point_status": "FAIL",
            "failure_reason": "ROOT_INDEX_INVENTORY_NOT_EXACTLY_1_TO_9",
        }

    strict_by_position: dict[int, bool] = {}
    for position in range(1, OUTPUT_GUARD_POSITION + 1):
        if position in refined_rows:
            strict_by_position[position] = bool(
                refined_rows[position]["local_adjudication_status"] == "PASS"
            )
        else:
            strict_by_position[position] = bool(
                float(initial_rows[position]["relative_difference_Omega"])
                <= INVENTORY_VERIFICATION_TOLERANCE
            )
    strict_first8 = all(
        strict_by_position[position]
        for position in range(1, PLOTTED_POSITIONS + 1)
    )
    root9_strict = strict_by_position[OUTPUT_GUARD_POSITION]

    def final_estimates_and_brackets(position: int) -> tuple[list[float], list[list[float]]]:
        if position in refined_rows:
            row = refined_rows[position]
            estimates = [float(value) for value in row["all_refined_values_Omega"]]
            brackets = [
                [float(value) for value in row[path]["local_sign_bracket"]]
                for path in ("primary_path", "verification_path")
            ]
            return estimates, brackets
        row = initial_rows[position]
        return [
            float(row["Omega_primary"]),
            float(row["Omega_verification"]),
        ], [
            [float(value) for value in row["primary_original_bracket"]],
            [
                float(value)
                for value in row.get(
                    "verification_original_bracket",
                    row.get("verification_reconstructed_bracket"),
                )
            ],
        ]

    root8_estimates, root8_brackets = final_estimates_and_brackets(8)
    root9_estimates, root9_brackets = final_estimates_and_brackets(9)

    def enclosure(estimates: Sequence[float], brackets: Sequence[Sequence[float]]) -> list[float]:
        return [
            min((*estimates, *(float(pair[0]) for pair in brackets))),
            max((*estimates, *(float(pair[1]) for pair in brackets))),
        ]

    root8_enclosure = enclosure(root8_estimates, root8_brackets)
    root9_enclosure = enclosure(root9_estimates, root9_brackets)
    enclosure_gap = root9_enclosure[0] - root8_enclosure[1]
    enclosures_disjoint_ordered = bool(enclosure_gap > 0.0)
    root9_common_bracket = [
        max(pair[0] for pair in root9_brackets),
        min(pair[1] for pair in root9_brackets),
    ]
    root9_bracket_overlap = bool(
        root9_common_bracket[0] < root9_common_bracket[1]
    )
    all_root9_estimates_inside_common_bracket = bool(
        root9_bracket_overlap
        and all(
            root9_common_bracket[0] <= value <= root9_common_bracket[1]
            for value in root9_estimates
        )
    )
    root9_spread = max(root9_estimates) - min(root9_estimates)
    nearest_gap = float(initial_rows[9]["nearest_adjacent_gap_Omega"])
    spread_to_gap = root9_spread / nearest_gap

    root9_refined = refined_rows.get(9)
    if root9_refined is not None:
        path_rows = [
            root9_refined["primary_path"],
            root9_refined["verification_path"],
        ]
        all_paths_same_neighbourhood = bool(
            root9_bracket_overlap
            and all_root9_estimates_inside_common_bracket
            and not bool(
                root9_refined.get(
                    "distinct_physical_candidate_between_estimates", True
                )
            )
            and all(row.get("sign_bracket_count") == 1 for row in path_rows)
            and all(row.get("status") == "PASS" for row in path_rows)
            and all(row.get("determinant_converged") for row in path_rows)
            and all(row.get("sigma_converged") for row in path_rows)
            and all(row.get("determinant_quality_pass") for row in path_rows)
            and all(row.get("sigma_quality_pass") for row in path_rows)
            and spread_to_gap <= LOCAL_ADJUDICATION_GAP_RATIO_LIMIT
        )
        root9_singular_ratios = [
            float(row[key])
            for row in path_rows
            for key in (
                "determinant_scaled_sigma_ratio",
                "sigma_scaled_sigma_ratio",
            )
        ]
        root9_boundary_residuals = [
            float(row[key])
            for row in path_rows
            for key in (
                "determinant_boundary_residual",
                "sigma_boundary_residual",
            )
        ]
        ordering_same = bool(
            all(row.get("ordering_preserved") for row in refined_rows.values())
        )
    else:
        all_paths_same_neighbourhood = bool(
            root9_bracket_overlap
            and spread_to_gap <= LOCAL_ADJUDICATION_GAP_RATIO_LIMIT
        )
        root9_singular_ratios = [
            float(initial_rows[9]["primary_scaled_sigma_ratio"]),
            float(initial_rows[9]["verification_scaled_sigma_ratio"]),
        ]
        root9_boundary_residuals = [
            float(initial_rows[9]["primary_boundary_residual"]),
            float(initial_rows[9]["verification_boundary_residual"]),
        ]
        ordering_same = True

    unresolved = int(attempt["unresolved_candidates_below_root9"])
    exact_eight_below = bool(
        int(attempt["root_count"]) == OUTPUT_GUARD_POSITION
        and len(initial_rows) == OUTPUT_GUARD_POSITION
        and all(
            initial_rows[position]["Omega_primary"]
            < initial_rows[position + 1]["Omega_primary"]
            and initial_rows[position]["Omega_verification"]
            < initial_rows[position + 1]["Omega_verification"]
            for position in range(1, OUTPUT_GUARD_POSITION)
        )
    )
    quality_pass = bool(
        max(root9_singular_ratios) <= 1.0e-9
        and max(root9_boundary_residuals) <= 1.0e-9
    )
    right_boundary_distance = float(attempt["effective_Omega_max"]) - root9_enclosure[1]
    right_boundary_pass = bool(
        attempt["guard_not_at_scan_boundary"]
        and right_boundary_distance >= policy.post_guard_tail_Omega
    )
    guard_pass = bool(
        strict_first8
        and not root9_strict
        and all_paths_same_neighbourhood
        and enclosures_disjoint_ordered
        and unresolved == 0
        and exact_eight_below
        and ordering_same
        and quality_pass
        and right_boundary_pass
    )
    return {
        "strict_agreement_threshold": INVENTORY_VERIFICATION_TOLERANCE,
        "strict_agreement_quantity": "RELATIVE_OMEGA",
        "strict_roots_1_to_8_status": "PASS" if strict_first8 else "FAIL",
        "strict_by_position": {
            str(position): "PASS" if status else "FAIL"
            for position, status in strict_by_position.items()
        },
        "root9_strict_agreement_status": "PASS" if root9_strict else "FAIL",
        "root9_guard_status": (
            ROOT9_GUARD_INTERVAL_PASS if guard_pass else ROOT9_GUARD_FAIL
        ),
        "point_status": (
            POINT_PASS_WITH_GUARD_QUALIFICATION if guard_pass else "FAIL"
        ),
        "root8_final_estimates_Omega": root8_estimates,
        "root9_final_estimates_Omega": root9_estimates,
        "root8_final_brackets_Omega": root8_brackets,
        "root9_final_brackets_Omega": root9_brackets,
        "root8_enclosure_Omega": root8_enclosure,
        "root9_enclosure_Omega": root9_enclosure,
        "enclosures_intersect": not enclosures_disjoint_ordered,
        "enclosure_gap_Omega": enclosure_gap,
        "root9_spread_Omega": root9_spread,
        "root9_spread_to_nearest_gap_ratio": spread_to_gap,
        "root9_nearest_adjacent_gap_Omega": nearest_gap,
        "all_root9_paths_one_isolated_neighbourhood": (
            all_paths_same_neighbourhood
        ),
        "root9_path_brackets_overlap": root9_bracket_overlap,
        "root9_common_bracket_intersection_Omega": root9_common_bracket,
        "all_root9_estimates_inside_common_bracket": (
            all_root9_estimates_inside_common_bracket
        ),
        "distinct_physical_candidate_between_root9_estimates": bool(
            root9_refined.get(
                "distinct_physical_candidate_between_estimates", True
            )
            if root9_refined is not None
            else False
        ),
        "unresolved_candidates_between_root8_and_root9": unresolved,
        "exactly_eight_ordered_roots_below_root9": exact_eight_below,
        "ordering_same_in_all_passes": ordering_same,
        "root9_scaled_singular_ratios": root9_singular_ratios,
        "root9_boundary_residuals": root9_boundary_residuals,
        "root9_quality_status": "PASS" if quality_pass else "FAIL",
        "root9_distance_to_Omega_max": right_boundary_distance,
        "root9_right_boundary_status": "PASS" if right_boundary_pass else "FAIL",
        "root9_plotted": False,
        "canonical_root9_selector": (
            root9_refined.get("canonical_selection")
            if root9_refined is not None
            else "PRIMARY_EXISTING_DETERMINISTIC_SELECTOR"
        ),
        "root_estimates_averaged": False,
    }


def _adjudicated_point_rows(
    model: str,
    geometry: GeometryPoint,
    pairs: Sequence[Mapping[str, Any]],
    canonical: Mapping[int, Any],
    adjudication: Mapping[str, Any],
    root9_guard: Mapping[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Build a certified nine-row group after a passing local adjudication."""

    strict_pass = adjudication["local_adjudication_status"] == "PASS"
    guard_pass = bool(
        root9_guard is not None
        and root9_guard.get("root9_guard_status")
        == ROOT9_GUARD_INTERVAL_PASS
    )
    if not (strict_pass or guard_pass):
        return []
    slots = []
    for pair in pairs:
        position = int(pair["root_index"])
        candidate = canonical[position]
        event = replace(
            pair["primary_event"],
            Omega=candidate.Omega,
            omega=candidate.diagnostics.omega,
            candidate=candidate,
            cluster_center_Omega=candidate.Omega,
        )
        slots.append(
            SimpleNamespace(
                sorted_slot=position,
                role=(
                    "REQUESTED"
                    if position <= PLOTTED_POSITIONS
                    else "GUARD"
                ),
                event=event,
            )
        )
    final_relative = [
        float(pair["initial_row"]["relative_difference_Omega"])
        for pair in pairs
    ]
    for outcome in adjudication["refined_root_diagnostics"]:
        final_relative[int(outcome["root_index"]) - 1] = float(
            outcome["final_relative_agreement_Omega"]
        )
    inventory_status = (
        "LOCAL_ADJUDICATION_PASS"
        if strict_pass
        else POINT_PASS_WITH_GUARD_QUALIFICATION
    )
    semantic = {
        "model": model,
        "mu": geometry.mu,
        "Omega": [slot.event.Omega for slot in slots],
        "final_relative": final_relative,
        "status": inventory_status,
    }
    inventory_sha256 = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode(
            "utf-8"
        )
    ).hexdigest().upper()
    inventory = SimpleNamespace(
        status=inventory_status,
        inventory_sha256=inventory_sha256,
        primary=SimpleNamespace(slots=slots),
        verification=SimpleNamespace(slots=slots),
        maximum_primary_verification_relative=max(final_relative),
        unresolved_low_sigma_count=0,
    )
    export = {
        "export_guard_available": True,
        "export_primary_verification_agreement": strict_pass,
        "export_primary_verification_max_relative": max(final_relative),
        "unresolved_candidates_below_export_guard": 0,
        "export_range_status": (
            "PASS" if strict_pass else POINT_PASS_WITH_GUARD_QUALIFICATION
        ),
        "plotted_first8_primary_verification_agreement": True,
        "plotted_first8_primary_verification_max_relative": max(
            final_relative[:PLOTTED_POSITIONS]
        ),
        "root9_primary_verification_agreement": strict_pass,
        "root9_primary_verification_relative": final_relative[-1],
        "primary_verification_relative_by_position": final_relative,
        "internal_inventory_status": "NOT_COMPUTED_ABOVE_ROOT9",
        "internal_primary_verification_max_relative": max(final_relative),
        "internal_unresolved_candidate_count": "NOT_COMPUTED_ABOVE_ROOT9",
    }
    return [
        transform_root_row(
            _bounded_root_row(model, geometry.beta_deg, inventory, slot),
            model,
            geometry,
            SWEEP_MU,
            export,
        )
        for slot in slots
    ]


def _pairs_from_live_inventory(
    inventory: Any,
    policy: Any,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Retain exact live brackets so agreement-only FAIL can be adjudicated."""

    primary = inventory.primary.slots[:OUTPUT_GUARD_POSITION]
    verification = inventory.verification.slots[:OUTPUT_GUARD_POSITION]
    if len(primary) != 9 or len(verification) != 9:
        raise RuntimeError("Live local adjudication requires exact roots 1--9.")
    saved = [float(slot.event.Omega) for slot in primary]
    pairs: list[dict[str, Any]] = []
    root_rows: list[dict[str, Any]] = []
    for index, (left_slot, right_slot) in enumerate(
        zip(primary, verification, strict=True)
    ):
        previous_gap: float | str = (
            saved[index] - saved[index - 1]
            if index > 0
            else "NOT_APPLICABLE_FIRST_ROOT"
        )
        next_gap: float | str = (
            saved[index + 1] - saved[index]
            if index + 1 < len(saved)
            else "NOT_COMPUTED_ABOVE_ROOT9"
        )
        nearest_gap = min(
            float(value)
            for value in (previous_gap, next_gap)
            if isinstance(value, float)
        )
        absolute = abs(left_slot.event.Omega - right_slot.event.Omega)
        relative = _relative(left_slot.event.Omega, right_slot.event.Omega)
        row = {
            "root_index": index + 1,
            "Omega_primary": left_slot.event.Omega,
            "Omega_verification": right_slot.event.Omega,
            "absolute_difference_Omega": absolute,
            "relative_difference_Omega": relative,
            "absolute_difference_Lambda": abs(
                math.sqrt(left_slot.event.Omega)
                - math.sqrt(right_slot.event.Omega)
            ),
            "relative_difference_Lambda": _relative(
                math.sqrt(left_slot.event.Omega),
                math.sqrt(right_slot.event.Omega),
            ),
            "primary_scaled_sigma_ratio": (
                left_slot.event.candidate.diagnostics.scaled_sigma_ratio
            ),
            "verification_scaled_sigma_ratio": (
                right_slot.event.candidate.diagnostics.scaled_sigma_ratio
            ),
            "primary_boundary_residual": (
                left_slot.event.candidate.diagnostics.raw_boundary_null_residual
            ),
            "verification_boundary_residual": (
                right_slot.event.candidate.diagnostics.raw_boundary_null_residual
            ),
            "previous_root_gap_Omega": previous_gap,
            "next_root_gap_Omega": next_gap,
            "nearest_adjacent_gap_Omega": nearest_gap,
            "initial_disagreement_to_nearest_gap_ratio": absolute
            / nearest_gap,
            "primary_original_detector_type": list(
                left_slot.event.candidate.detection_sources
            ),
            "verification_original_detector_type": list(
                right_slot.event.candidate.detection_sources
            ),
            "primary_original_bracket": [
                left_slot.event.candidate.interval_left_Omega,
                left_slot.event.candidate.interval_right_Omega,
            ],
            "verification_original_bracket": [
                right_slot.event.candidate.interval_left_Omega,
                right_slot.event.candidate.interval_right_Omega,
            ],
        }
        root_rows.append(row)
        neighbourhood_left = (
            0.5 * (saved[index - 1] + saved[index])
            if index > 0
            else policy.Omega_min
        )
        neighbourhood_right = (
            0.5 * (saved[index] + saved[index + 1])
            if index + 1 < len(saved)
            else policy.Omega_max
        )
        pairs.append(
            {
                "root_index": index + 1,
                "persisted_primary": {
                    "Omega": left_slot.event.Omega,
                    "bracket_left_Omega": (
                        left_slot.event.candidate.interval_left_Omega
                    ),
                    "bracket_right_Omega": (
                        left_slot.event.candidate.interval_right_Omega
                    ),
                },
                "primary_event": left_slot.event,
                "verification_event": right_slot.event,
                "primary_neighbourhood_events": [left_slot.event],
                "verification_neighbourhood_events": [right_slot.event],
                "neighbourhood_left": neighbourhood_left,
                "neighbourhood_right": neighbourhood_right,
                "nearest_gap": nearest_gap,
                "initial_row": row,
            }
        )
    maximum = max(root_rows, key=lambda row: row["relative_difference_Omega"])
    return {
        "source_checkpoint_evidence": "LIVE_INVENTORY_EXACT_PRIMARY_AND_VERIFICATION",
        "comparison_quantity": "Omega",
        "comparison_threshold_type": "RELATIVE",
        "comparison_threshold": INVENTORY_VERIFICATION_TOLERANCE,
        "Lambda_used_for_gate": False,
        "refiner_contract": _frozen_refiner_contract(policy),
        "root_diagnostics": root_rows,
        "maximum_relative_difference_Omega": float(
            maximum["relative_difference_Omega"]
        ),
        "maximum_difference_root_index": int(maximum["root_index"]),
        "global_scan_run": False,
    }, pairs


def _local_adjudication_path(output_dir: Path, model: str, mu: float) -> Path:
    model_token = model.lower().replace("_", "-")
    mu_token = format(float(mu), ".2f").replace(".", "p")
    return Path(output_dir) / f"{model_token}_mu_{mu_token}_local_adjudication.json"


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


def _protected_closing_hashes() -> dict[str, str]:
    return {
        relative: rlb2c.rlb2b._sha256(ROOT / relative)
        for relative in PROTECTED_CLOSING_PATHS
    }


def _initialize_production_closing_session(
    output_dir: Path,
    *,
    allow_adjudication_script_revision: bool = False,
    allow_root9_guard_revision: bool = False,
) -> dict[str, Any]:
    """Freeze the append-only launch state before the first production point."""

    target = Path(output_dir)
    checkpoint_path = target / CLOSING_CHECKPOINT_FILENAME
    checkpoint = (
        json.loads(checkpoint_path.read_text(encoding="utf-8"))
        if checkpoint_path.is_file()
        else {"schema_version": 1, "stage": "RLB-2D-CLOSING"}
    )
    existing_session = checkpoint.get(PRODUCTION_CLOSING_SESSION_KEY)
    if existing_session is not None:
        expected_script = str(existing_session["rlb2d_script_sha256"])
        current_script = rlb2c.rlb2b._sha256(Path(__file__))
        if current_script != expected_script:
            if not (
                allow_adjudication_script_revision
                or allow_root9_guard_revision
            ):
                raise RuntimeError(
                    "The RLB-2D script changed after the production closing "
                    "session was frozen."
                )
            audit = closing_mu_audit(target, include_rlb_shards=False)
            point_records = list(existing_session.get("point_records", []))
            benchmark_records = list(
                existing_session.get("benchmark_records", [])
            )
            if allow_root9_guard_revision:
                current_csv_hashes = {
                    model: rlb2c.rlb2b._sha256(
                        target / ROOT_FILENAMES[(SWEEP_MU, model)]
                    )
                    for model in MODELS
                }
                preservation = _verify_preexisting_mu_groups(
                    target, existing_session
                )
                weak_diagnostic = target / WEAK_RLB_MU080_ADJUDICATION_FILENAME
                old_diagnostic = target / _local_adjudication_path(
                    Path("."), MODEL_OLD, 0.8
                ).name
                records_exact = bool(
                    len(point_records) == 3
                    and benchmark_records == point_records
                    and point_records[0].get("model") == MODEL_RLB
                    and point_records[0].get("status") == "FAIL"
                    and point_records[1].get("model") == MODEL_RLB
                    and point_records[1].get("status") == "PASS"
                    and point_records[1].get("final_point_status")
                    == "LOCAL_ADJUDICATION_PASS"
                    and point_records[2].get("model") == MODEL_OLD
                    and point_records[2].get("status") == "FAIL"
                    and point_records[2].get("failure_reason")
                    == "LOCAL_ADJUDICATION_DID_NOT_CONVERGE"
                    and all(
                        math.isclose(
                            float(record.get("mu", math.nan)),
                            0.8,
                            rel_tol=0.0,
                            abs_tol=5.0e-13,
                        )
                        for record in point_records
                    )
                )
                current_groups = preservation["current_group_sha256"]
                initial_groups = existing_session[
                    "preexisting_mu_group_sha256"
                ]
                extra_groups = {
                    model: set(current_groups[model]) - set(initial_groups[model])
                    for model in MODELS
                }
                common_state_exact = bool(
                    records_exact
                    and preservation["status"] == "PASS"
                    and existing_session.get("computed_points")
                    == {MODEL_OLD: [], MODEL_RLB: [0.8]}
                    and int(existing_session.get("new_points_executed", -1)) == 1
                    and int(existing_session.get("ready_points_recalculated", -1))
                    == 0
                    and current_csv_hashes[MODEL_EB]
                    == existing_session["initial_mu_csv_sha256"][MODEL_EB]
                    and current_csv_hashes[MODEL_RLB]
                    == point_records[1].get("canonical_csv_sha256_after")
                    and weak_diagnostic.is_file()
                    and old_diagnostic.is_file()
                    and existing_session.get("eta_gate", {}).get("status")
                    == "NOT_EVALUATED"
                    and _frozen_hashes()
                    == existing_session["production_physics_hashes_before"]
                )
                pre_append_state = bool(
                    extra_groups
                    == {
                        MODEL_EB: set(),
                        MODEL_OLD: set(),
                        MODEL_RLB: {"0.8"},
                    }
                    and all(
                        audit[model]["missing_mu"]
                        == existing_session["current_missing_mu"][model]
                        for model in (MODEL_OLD, MODEL_RLB)
                    )
                    and current_csv_hashes
                    == existing_session["current_mu_csv_sha256"]
                    and current_csv_hashes[MODEL_OLD]
                    == existing_session["initial_mu_csv_sha256"][MODEL_OLD]
                )
                old_payload = (
                    json.loads(old_diagnostic.read_text(encoding="utf-8"))
                    if old_diagnostic.is_file()
                    else {}
                )
                recovery_state = bool(
                    extra_groups
                    == {
                        MODEL_EB: set(),
                        MODEL_OLD: {"0.8"},
                        MODEL_RLB: {"0.8"},
                    }
                    and audit[MODEL_RLB]["missing_mu"]
                    == existing_session["current_missing_mu"][MODEL_RLB]
                    and audit[MODEL_OLD]["missing_mu"]
                    == [
                        value
                        for value in existing_session["current_missing_mu"][
                            MODEL_OLD
                        ]
                        if not math.isclose(
                            float(value),
                            0.8,
                            rel_tol=0.0,
                            abs_tol=5.0e-13,
                        )
                    ]
                    and current_csv_hashes[MODEL_RLB]
                    == existing_session["current_mu_csv_sha256"][MODEL_RLB]
                    and old_payload.get("production_csv_changed") is True
                    and old_payload.get("root9_guard_contract", {}).get(
                        "root9_guard_status"
                    )
                    == ROOT9_GUARD_INTERVAL_PASS
                )
                state_exact = common_state_exact and (
                    pre_append_state or recovery_state
                )
                if not state_exact:
                    raise RuntimeError(
                        "The current weak-RLB append and saved old-Timoshenko "
                        "FAIL are not in the exact state required for the "
                        "root-9 guard migration."
                    )
                migration = {
                    "migrated_utc": _utc_now(),
                    "reason": (
                        "NARROW_ROOT9_GUARD_TRANSACTION_RECOVERY"
                        if recovery_state
                        else "NARROW_ROOT9_COMPLETENESS_GUARD_CONTRACT"
                    ),
                    "previous_rlb2d_script_sha256": expected_script,
                    "current_rlb2d_script_sha256": current_script,
                    "failed_benchmark_records_preserved": True,
                    "computed_points_before_migration": {
                        MODEL_OLD: [],
                        MODEL_RLB: [0.8],
                    },
                    "current_mu_csv_sha256": current_csv_hashes,
                    "weak_diagnostic_sha256": rlb2c.rlb2b._sha256(
                        weak_diagnostic
                    ),
                    "old_diagnostic_sha256": rlb2c.rlb2b._sha256(
                        old_diagnostic
                    ),
                    "production_csv_append_already_complete": recovery_state,
                }
            else:
                failed_record = (
                    point_records[0] if len(point_records) == 1 else {}
                )
                unchanged_mu_csv = all(
                    rlb2c.rlb2b._sha256(
                        target / ROOT_FILENAMES[(SWEEP_MU, model)]
                    )
                    == existing_session["initial_mu_csv_sha256"][model]
                    for model in MODELS
                )
                migration_allowed = bool(
                    len(point_records) == 1
                    and len(benchmark_records) == 1
                    and benchmark_records[0] == failed_record
                    and failed_record.get("model") == MODEL_RLB
                    and math.isclose(
                        float(failed_record.get("mu", math.nan)),
                        0.8,
                        rel_tol=0.0,
                        abs_tol=5.0e-13,
                    )
                    and failed_record.get("status") == "FAIL"
                    and not any(
                        existing_session.get("computed_points", {}).get(
                            model, []
                        )
                        for model in (MODEL_OLD, MODEL_RLB)
                    )
                    and unchanged_mu_csv
                    and all(
                        audit[model]["missing_mu"]
                        == existing_session["initial_missing_mu"][model]
                        for model in (MODEL_OLD, MODEL_RLB)
                    )
                    and _frozen_hashes()
                    == existing_session["production_physics_hashes_before"]
                )
                if not migration_allowed:
                    raise RuntimeError(
                        "The saved failed weak-RLB benchmark is not in the "
                        "exact append-free state required for adjudication "
                        "migration."
                    )
                migration = {
                    "migrated_utc": _utc_now(),
                    "reason": (
                        "NARROW_WEAK_RLB_MU080_LOCAL_ADJUDICATION_POLICY"
                    ),
                    "previous_rlb2d_script_sha256": expected_script,
                    "current_rlb2d_script_sha256": current_script,
                    "failed_benchmark_record_preserved": True,
                    "production_mu_rows_appended_before_migration": 0,
                }
            existing_session.setdefault("script_hash_migrations", []).append(
                migration
            )
            existing_session["rlb2d_script_sha256"] = current_script
            checkpoint[PRODUCTION_CLOSING_SESSION_KEY] = existing_session
            _atomic_write_json(checkpoint_path, checkpoint)
        if _protected_closing_hashes() != existing_session[
            "protected_closing_hashes_before"
        ]:
            raise RuntimeError("A protected closing file changed during RLB-2D.")
        for model, expected in existing_session["beta_csv_sha256"].items():
            actual = rlb2c.rlb2b._sha256(
                target / ROOT_FILENAMES[(SWEEP_BETA, model)]
            )
            if actual != expected:
                raise RuntimeError(f"The frozen beta CSV changed for {model}.")
        beta_plot = target / PLOT_FILENAMES[SWEEP_BETA]
        if rlb2c.rlb2b._sha256(beta_plot) != existing_session[
            "beta_plot_sha256_before"
        ]:
            raise RuntimeError("The frozen beta plot changed during RLB-2D.")
        return checkpoint

    audit = closing_mu_audit(target, include_rlb_shards=False)
    mu_rows = {
        model: _canonical_mu_rows(target, model) for model in MODELS
    }
    schemas = {model: list(rows[0]) for model, rows in mu_rows.items()}
    group_hashes = {
        model: _mu_group_hashes(rows, model, schemas[model])
        for model, rows in mu_rows.items()
    }
    beta_rows = {
        model: rlb2c.rlb2b._read_csv(
            target / ROOT_FILENAMES[(SWEEP_BETA, model)]
        )
        for model in MODELS
    }
    beta_summaries = {
        model: validate_root_rows(SWEEP_BETA, model, rows, beta_grid())
        for model, rows in beta_rows.items()
    }
    if not all(
        summary["exact_row_structure_passed"]
        for summary in beta_summaries.values()
    ):
        raise RuntimeError("The frozen beta sweep is structurally incomplete.")
    beta_plot = target / PLOT_FILENAMES[SWEEP_BETA]
    if not beta_plot.is_file() or beta_plot.stat().st_size == 0:
        raise RuntimeError("The completed beta plot is missing.")
    previous_checkpoint_hash = (
        rlb2c.rlb2b._sha256(checkpoint_path)
        if checkpoint_path.is_file()
        else None
    )
    initial_mu_complete = {
        model: len(_complete_mu_values(rows, model))
        for model, rows in mu_rows.items()
    }
    initial_beta_complete = {
        model: int(beta_summaries[model]["axis_count"]) for model in MODELS
    }
    session = {
        "schema_version": 1,
        "started_utc": _utc_now(),
        "git_state_at_start": rlb2c.rlb2b._git_state(),
        "rlb2d_script_sha256": rlb2c.rlb2b._sha256(Path(__file__)),
        "production_physics_hashes_before": _frozen_hashes(),
        "protected_closing_hashes_before": _protected_closing_hashes(),
        "previous_checkpoint_sha256": previous_checkpoint_hash,
        "initial_missing_mu": {
            model: list(audit[model]["missing_mu"])
            for model in (MODEL_OLD, MODEL_RLB)
        },
        "initial_complete_mu_groups": initial_mu_complete,
        "initial_complete_beta_groups": initial_beta_complete,
        "initial_reused_parameter_points": (
            sum(initial_mu_complete.values()) + sum(initial_beta_complete.values())
        ),
        "initial_mu_csv_sha256": {
            model: rlb2c.rlb2b._sha256(
                target / ROOT_FILENAMES[(SWEEP_MU, model)]
            )
            for model in MODELS
        },
        "initial_mu_csv_fieldnames": schemas,
        "preexisting_mu_group_sha256": group_hashes,
        "beta_csv_sha256": {
            model: rlb2c.rlb2b._sha256(
                target / ROOT_FILENAMES[(SWEEP_BETA, model)]
            )
            for model in MODELS
        },
        "beta_plot_path": PLOT_FILENAMES[SWEEP_BETA],
        "beta_plot_sha256_before": rlb2c.rlb2b._sha256(beta_plot),
        "beta_plot_reused_without_redraw": True,
        "thread_limits": dict(CLOSING_THREAD_LIMITS),
        "root_contract": {
            "plotted_positions": list(range(1, PLOTTED_POSITIONS + 1)),
            "guard_position": OUTPUT_GUARD_POSITION,
            "roots_above_9_requested_or_exported": False,
            "initial_bound_factor": PRODUCTION_INITIAL_BOUND_FACTOR,
            "single_allowed_expansion_factor": (
                PRODUCTION_SINGLE_EXPANSION_FACTOR
            ),
        },
        "point_records": [],
        "benchmark_records": [],
        "failed_points": [],
        "computed_points": {MODEL_OLD: [], MODEL_RLB: []},
        "ready_points_recalculated": 0,
        "preexisting_groups_preserved": True,
        "eta_gate": {"status": "NOT_EVALUATED"},
    }
    checkpoint[PRODUCTION_CLOSING_SESSION_KEY] = session
    checkpoint.update(
        {
            "initial_missing_mu": session["initial_missing_mu"],
            "initial_complete_point_count": session[
                "initial_reused_parameter_points"
            ],
            "initial_reused_point_count": session[
                "initial_reused_parameter_points"
            ],
            "inherited_workers_drained": 0,
            "inherited_points_allowed_to_finish": {
                MODEL_OLD: [],
                MODEL_RLB: [],
            },
            "parallel_workers_used_in_closing_stage": 0,
            "ready_points_recalculated": 0,
            "global_restarts": 0,
        }
    )
    _atomic_write_json(checkpoint_path, checkpoint)
    return checkpoint


def _verify_preexisting_mu_groups(
    output_dir: Path, session: Mapping[str, Any]
) -> dict[str, Any]:
    target = Path(output_dir)
    mismatches: list[dict[str, str]] = []
    current_hashes: dict[str, dict[str, str]] = {}
    for model in MODELS:
        rows = _canonical_mu_rows(target, model)
        fields = [str(value) for value in session["initial_mu_csv_fieldnames"][model]]
        hashes = _mu_group_hashes(rows, model, fields)
        current_hashes[model] = hashes
        for mu_key, expected in session["preexisting_mu_group_sha256"][model].items():
            actual = hashes.get(mu_key)
            if actual != expected:
                mismatches.append(
                    {
                        "model": model,
                        "mu": mu_key,
                        "expected_sha256": expected,
                        "actual_sha256": actual or "MISSING",
                    }
                )
    return {
        "status": "PASS" if not mismatches else "FAIL",
        "mismatches": mismatches,
        "current_group_sha256": current_hashes,
    }


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


def _run_bounded_root9_attempt(
    model: str,
    geometry: GeometryPoint,
    objects: ModelObjects,
    contract_sha256: str,
    requested_Omega_max: float,
) -> tuple[list[dict[str, Any]], dict[str, Any], Any | None]:
    """Run the existing determinant/SVD inventory only through root 9."""

    started = time.perf_counter()
    base_policy = rlb2c.rlb2b.frozen_root_policy()
    policy, policy_record = _bounded_root9_policy(
        base_policy, requested_Omega_max
    )
    provider, provider_metadata = make_matrix_provider(
        model, geometry, objects, cache=ArmPairCache.empty()
    )
    evaluation_counter = 0

    def counted_provider(omega: float) -> FloatArray:
        nonlocal evaluation_counter
        evaluation_counter += 1
        return provider(float(omega))

    try:
        inventory = rlb2c.rlb2b.iso_inventory.seed_free_root_inventory(
            counted_provider,
            reference_omega_scale(),
            policy,
            case_id=(
                f"rlb2d_bounded_root9__{model.lower()}__"
                f"mu_{geometry.mu:.4f}"
            ),
            builder_id=f"RLB2D_BOUNDED_ROOT9_{model}",
            contract_sha256=contract_sha256,
        )
        export = _inventory_export_diagnostics(inventory, policy)
        raw_rows = [
            _bounded_root_row(model, geometry.beta_deg, inventory, slot)
            for slot in inventory.slots[:OUTPUT_GUARD_POSITION]
        ]
        point_rows = [
            transform_root_row(row, model, geometry, SWEEP_MU, export)
            for row in raw_rows
        ]
        roots = [float(row["Omega"]) for row in point_rows]
        root_quality = bool(
            len(point_rows) == OUTPUT_GUARD_POSITION
            and all(str(row["root_status"]) == "PASS" for row in point_rows)
        )
        disposition = _bounded_attempt_disposition(inventory, export, policy)
        passed = bool(disposition == ATTEMPT_ACCEPT and root_quality)
        primary_slots = inventory.primary.slots[:OUTPUT_GUARD_POSITION]
        verification_slots = inventory.verification.slots[
            :OUTPUT_GUARD_POSITION
        ]
        paired_evidence = []
        if (
            len(primary_slots) == OUTPUT_GUARD_POSITION
            and len(verification_slots) == OUTPUT_GUARD_POSITION
        ):
            for primary_slot, verification_slot in zip(
                primary_slots, verification_slots, strict=True
            ):
                primary_candidate = primary_slot.event.candidate
                verification_candidate = verification_slot.event.candidate
                paired_evidence.append(
                    {
                        "root_index": int(primary_slot.sorted_slot),
                        "Omega_primary": primary_slot.event.Omega,
                        "Omega_verification": verification_slot.event.Omega,
                        "absolute_difference_Omega": abs(
                            primary_slot.event.Omega
                            - verification_slot.event.Omega
                        ),
                        "relative_difference_Omega": _relative(
                            primary_slot.event.Omega,
                            verification_slot.event.Omega,
                        ),
                        "absolute_difference_Lambda": abs(
                            math.sqrt(primary_slot.event.Omega)
                            - math.sqrt(verification_slot.event.Omega)
                        ),
                        "primary_bracket": [
                            primary_candidate.interval_left_Omega,
                            primary_candidate.interval_right_Omega,
                        ],
                        "verification_bracket": [
                            verification_candidate.interval_left_Omega,
                            verification_candidate.interval_right_Omega,
                        ],
                        "primary_detector_sources": list(
                            primary_candidate.detection_sources
                        ),
                        "verification_detector_sources": list(
                            verification_candidate.detection_sources
                        ),
                        "primary_scaled_sigma_ratio": (
                            primary_candidate.diagnostics.scaled_sigma_ratio
                        ),
                        "verification_scaled_sigma_ratio": (
                            verification_candidate.diagnostics.scaled_sigma_ratio
                        ),
                        "primary_boundary_residual": (
                            primary_candidate.diagnostics.raw_boundary_null_residual
                        ),
                        "verification_boundary_residual": (
                            verification_candidate.diagnostics.raw_boundary_null_residual
                        ),
                    }
                )
        metrics = {
            **policy_record,
            "status": "PASS" if passed else "FAIL",
            "disposition": disposition,
            "inventory_status_first8_plus_root9": inventory.status,
            "export_range_status": export["export_range_status"],
            "guard_available": inventory.guard_available,
            "guard_not_at_scan_boundary": inventory.guard_not_at_scan_boundary,
            "independent_primary_verification_agreement": (
                inventory.independent_agreement
            ),
            "maximum_primary_verification_relative": (
                inventory.maximum_primary_verification_relative
            ),
            "unresolved_candidates_below_root9": (
                inventory.unresolved_low_sigma_count
            ),
            "determinant_evaluations": evaluation_counter,
            "determinant_scans": 2,
            "sigma_scans": 2,
            "provider_setup_metadata": provider_metadata,
            "root_count": len(point_rows),
            "incidental_primary_events_above_root9": max(
                0, len(inventory.primary.events) - OUTPUT_GUARD_POSITION
            ),
            "incidental_verification_events_above_root9": max(
                0, len(inventory.verification.events) - OUTPUT_GUARD_POSITION
            ),
            "roots": [
                {
                    "sorted_position": int(row["sorted_position"]),
                    "omega": float(row["omega"]),
                    "Omega": float(row["Omega"]),
                    "Lambda": float(row["Lambda"]),
                    "scaled_sigma_ratio": float(row["scaled_sigma_ratio"]),
                    "boundary_null_residual": float(
                        row["boundary_null_residual"]
                    ),
                    "root_status": str(row["root_status"]),
                    "bracket_left_Omega": float(row["bracket_left_Omega"]),
                    "bracket_right_Omega": float(row["bracket_right_Omega"]),
                }
                for row in point_rows
            ],
            "Omega_roots": roots,
            "inventory_sha256": inventory.inventory_sha256,
            "primary_verification_by_position": paired_evidence,
            "refiner_contract": _frozen_refiner_contract(policy),
            "wall_time_seconds": time.perf_counter() - started,
            "error": None,
        }
        if passed:
            export = dict(export)
            export["internal_inventory_status"] = "NOT_COMPUTED_ABOVE_ROOT9"
            export["internal_primary_verification_max_relative"] = (
                inventory.maximum_primary_verification_relative
            )
            export["internal_unresolved_candidate_count"] = (
                "NOT_COMPUTED_ABOVE_ROOT9"
            )
            point_rows = [
                transform_root_row(
                    row, model, geometry, SWEEP_MU, export
                )
                for row in raw_rows
            ]
        return point_rows if passed else [], metrics, inventory
    except Exception as exc:  # numerical failure is recorded before one retry
        return [], {
            **policy_record,
            "status": "FAIL",
            "disposition": ATTEMPT_HARD_FAIL,
            "determinant_evaluations": evaluation_counter,
            "determinant_scans": None,
            "sigma_scans": None,
            "detector_scan_attempt_limit": 2,
            "provider_setup_metadata": provider_metadata,
            "root_count": 0,
            "roots": [],
            "Omega_roots": [],
            "wall_time_seconds": time.perf_counter() - started,
            "error": f"{type(exc).__name__}: {exc}",
        }, None


def _build_root9_guard_certificate(
    output_dir: Path,
    model: str,
    mu: float,
    point_rows: Sequence[Mapping[str, Any]],
    fieldnames: Sequence[str],
    diagnostic_path: Path,
    guard_contract: Mapping[str, Any],
) -> dict[str, Any]:
    """Bind a root-9-only qualification to its diagnostic and CSV group."""

    if guard_contract.get("strict_roots_1_to_8_status") != "PASS":
        raise RuntimeError("A root-9 guard cannot exempt roots 1--8.")
    if not (
        guard_contract.get("root9_strict_agreement_status") == "FAIL"
        and guard_contract.get("root9_guard_status")
        == ROOT9_GUARD_INTERVAL_PASS
        and guard_contract.get("point_status")
        == POINT_PASS_WITH_GUARD_QUALIFICATION
    ):
        raise RuntimeError("The proposed root-9 guard contract is not passing.")
    persisted_rows = [
        row
        for row in _canonical_mu_rows(Path(output_dir), model)
        if round(float(row["mu"]), 12) == round(float(mu), 12)
    ]
    certificate_rows = (
        persisted_rows
        if len(persisted_rows) == OUTPUT_GUARD_POSITION
        else list(point_rows)
    )
    inventory_hashes = {
        str(row["inventory_sha256"]) for row in certificate_rows
    }
    if len(inventory_hashes) != 1:
        raise RuntimeError("A guard-qualified group has inconsistent inventory hashes.")
    target = Path(output_dir).resolve()
    diagnostic = Path(diagnostic_path).resolve()
    if target not in diagnostic.parents or not diagnostic.is_file():
        raise RuntimeError("The guard diagnostic is outside the result directory.")
    return {
        "schema_version": 1,
        "sweep": SWEEP_MU,
        "model": model,
        "mu": round(float(mu), 12),
        "strict_frequency_positions": list(
            range(1, PLOTTED_POSITIONS + 1)
        ),
        "strict_agreement_threshold": INVENTORY_VERIFICATION_TOLERANCE,
        "strict_agreement_quantity": "RELATIVE_OMEGA",
        "strict_roots_1_to_8_status": "PASS",
        "root9_strict_agreement_status": "FAIL",
        "root9_guard_status": ROOT9_GUARD_INTERVAL_PASS,
        "point_status": POINT_PASS_WITH_GUARD_QUALIFICATION,
        "root9_plotted": False,
        "diagnostic_path": str(diagnostic.relative_to(target)),
        "diagnostic_sha256": rlb2c.rlb2b._sha256(diagnostic),
        "csv_group_sha256": _semantic_group_hash(
            certificate_rows, fieldnames
        ),
        "inventory_sha256": next(iter(inventory_hashes)),
        "root8_enclosure_Omega": list(
            guard_contract["root8_enclosure_Omega"]
        ),
        "root9_enclosure_Omega": list(
            guard_contract["root9_enclosure_Omega"]
        ),
        "enclosure_gap_Omega": float(
            guard_contract["enclosure_gap_Omega"]
        ),
        "root9_spread_Omega": float(guard_contract["root9_spread_Omega"]),
        "root9_spread_to_nearest_gap_ratio": float(
            guard_contract["root9_spread_to_nearest_gap_ratio"]
        ),
        "canonical_root9_selector": guard_contract[
            "canonical_root9_selector"
        ],
        "root_estimates_averaged": False,
    }


def _root9_guard_certificate_index(
    checkpoint: Mapping[str, Any], output_dir: Path
) -> dict[tuple[str, str, float], Mapping[str, Any]]:
    """Return only hash-bound, sidecar-verified root-9 certificates."""

    target = Path(output_dir)
    session = checkpoint[PRODUCTION_CLOSING_SESSION_KEY]
    certificates: dict[tuple[str, str, float], Mapping[str, Any]] = {}
    for record in session.get("point_records", []):
        certificate = record.get("root9_guard_certificate")
        if not isinstance(certificate, Mapping):
            continue
        model = str(certificate.get("model"))
        mu_value = round(float(certificate.get("mu", math.nan)), 12)
        key = (SWEEP_MU, model, mu_value)
        if key in certificates:
            raise RuntimeError(f"Duplicate root-9 guard certificate: {key}.")
        if not (
            certificate.get("sweep") == SWEEP_MU
            and model in {MODEL_OLD, MODEL_RLB}
            and certificate.get("strict_roots_1_to_8_status") == "PASS"
            and certificate.get("root9_strict_agreement_status") == "FAIL"
            and certificate.get("root9_guard_status")
            == ROOT9_GUARD_INTERVAL_PASS
            and certificate.get("point_status")
            == POINT_PASS_WITH_GUARD_QUALIFICATION
            and record.get("status") == "PASS"
            and record.get("final_point_status")
            == POINT_PASS_WITH_GUARD_QUALIFICATION
        ):
            raise RuntimeError(f"Invalid root-9 certificate statuses: {key}.")
        diagnostic = (target / str(certificate["diagnostic_path"])).resolve()
        if target.resolve() not in diagnostic.parents or not diagnostic.is_file():
            raise RuntimeError(f"Missing root-9 diagnostic for {key}.")
        if rlb2c.rlb2b._sha256(diagnostic) != certificate["diagnostic_sha256"]:
            raise RuntimeError(f"Root-9 diagnostic hash mismatch for {key}.")
        payload = json.loads(diagnostic.read_text(encoding="utf-8"))
        guard = payload.get("root9_guard_contract", {})
        if not (
            guard.get("strict_roots_1_to_8_status") == "PASS"
            and guard.get("root9_strict_agreement_status") == "FAIL"
            and guard.get("root9_guard_status") == ROOT9_GUARD_INTERVAL_PASS
            and guard.get("point_status")
            == POINT_PASS_WITH_GUARD_QUALIFICATION
        ):
            raise RuntimeError(f"Diagnostic guard contract mismatch for {key}.")
        rows = [
            row
            for row in _canonical_mu_rows(target, model)
            if round(float(row["mu"]), 12) == mu_value
        ]
        fields = [
            str(value)
            for value in session["initial_mu_csv_fieldnames"][model]
        ]
        if [int(row["sorted_position"]) for row in rows] != list(
            range(1, OUTPUT_GUARD_POSITION + 1)
        ):
            raise RuntimeError(f"Certified CSV group is incomplete for {key}.")
        if _semantic_group_hash(rows, fields) != certificate["csv_group_sha256"]:
            raise RuntimeError(f"Certified CSV group hash mismatch for {key}.")
        inventory_hashes = {str(row["inventory_sha256"]) for row in rows}
        if inventory_hashes != {str(certificate["inventory_sha256"])}:
            raise RuntimeError(f"Certified inventory hash mismatch for {key}.")
        if not all(
            str(row["export_range_status"])
            == POINT_PASS_WITH_GUARD_QUALIFICATION
            and not _bool(row["export_primary_verification_agreement"])
            for row in rows
        ):
            raise RuntimeError(f"Certified CSV strict evidence changed for {key}.")
        certificates[key] = certificate
    return certificates


def _point_qualification(
    rows: Sequence[Mapping[str, Any]],
    *,
    root9_guard_certificate: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    first = rows[0]
    export_pass = str(first["export_range_status"]) == "PASS"
    guard_sentinel = bool(
        str(first["export_range_status"])
        == POINT_PASS_WITH_GUARD_QUALIFICATION
    )
    guard_qualified = bool(
        guard_sentinel
        and root9_guard_certificate is not None
        and root9_guard_certificate.get("strict_roots_1_to_8_status")
        == "PASS"
        and root9_guard_certificate.get("root9_strict_agreement_status")
        == "FAIL"
        and root9_guard_certificate.get("root9_guard_status")
        == ROOT9_GUARD_INTERVAL_PASS
        and root9_guard_certificate.get("point_status")
        == POINT_PASS_WITH_GUARD_QUALIFICATION
    )
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
            and (
                str(first["root9_primary_verification_agreement"]).lower()
                == "true"
                or guard_qualified
            )
            else "FAIL"
        )
    elif (export_pass or guard_qualified) and first8_quality and guard_quality:
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
        "root9_guard_certificate_status": (
            ROOT9_GUARD_INTERVAL_PASS if guard_qualified else "NOT_APPLICABLE"
        ),
        "aggregate_primary_verification_max_relative": float(
            first["export_primary_verification_max_relative"]
        ),
        "internal_inventory_status": internal,
    }


def update_closing_checkpoint(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    *,
    inherited_workers_drained: int | None = None,
    point_record: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Refresh the launch-scoped append-only production checkpoint."""

    target = Path(output_dir)
    path = target / CLOSING_CHECKPOINT_FILENAME
    checkpoint = _initialize_production_closing_session(target)
    session = checkpoint[PRODUCTION_CLOSING_SESSION_KEY]
    if point_record is not None:
        record = rlb2c.rlb2b._json_value(dict(point_record))
        session["point_records"].append(record)
        if record.get("production_benchmark"):
            session["benchmark_records"].append(record)
        if record.get("status") != "PASS":
            session["failed_points"].append(
                {
                    "model": record["model"],
                    "mu": record["mu"],
                    "reason": record.get("failure_reason"),
                }
            )
    audit = closing_mu_audit(target, include_rlb_shards=False)
    preservation = _verify_preexisting_mu_groups(target, session)
    if preservation["status"] != "PASS":
        raise RuntimeError(
            f"Pre-existing mu groups changed: {preservation['mismatches']}."
        )
    initial_missing = session["initial_missing_mu"]
    completed: dict[str, list[float]] = {}
    qualifications: list[dict[str, Any]] = []
    guard_certificates = _root9_guard_certificate_index(checkpoint, target)
    for model in (MODEL_OLD, MODEL_RLB):
        initial = [
            round(float(value), 12) for value in initial_missing.get(model, [])
        ]
        missing = set(audit[model]["missing_mu"])
        completed[model] = [value for value in initial if value not in missing]
        rows = _canonical_mu_rows(target, model)
        groups: dict[float, list[Mapping[str, Any]]] = {}
        for row in rows:
            groups.setdefault(round(float(row["mu"]), 12), []).append(row)
        for value in completed[model]:
            qualifications.append(
                _point_qualification(
                    groups[value],
                    root9_guard_certificate=guard_certificates.get(
                        (SWEEP_MU, model, value)
                    ),
                )
            )
    session.update(
        {
            "last_updated_utc": _utc_now(),
            "computed_points": completed,
            "current_missing_mu": {
                model: audit[model]["missing_mu"]
                for model in (MODEL_OLD, MODEL_RLB)
            },
            "new_points_executed": sum(
                len(values) for values in completed.values()
            ),
            "ready_points_recalculated": 0,
            "latest_canonical_audit": audit,
            "new_point_qualifications": qualifications,
            "preexisting_group_preservation": preservation,
            "preexisting_groups_preserved": True,
            "current_mu_csv_sha256": {
                model: rlb2c.rlb2b._sha256(
                    target / ROOT_FILENAMES[(SWEEP_MU, model)]
                )
                for model in MODELS
            },
            "cumulative_point_wall_time_seconds": sum(
                float(record.get("wall_time_seconds", 0.0))
                for record in session["point_records"]
            ),
            "peak_RSS_bytes": max(
                (
                    int(record.get("peak_RSS_bytes", 0))
                    for record in session["point_records"]
                ),
                default=0,
            ),
        }
    )
    checkpoint.update(
        {
            "initial_missing_mu": initial_missing,
            "completed_new_mu_points": completed,
            "current_missing_mu": session["current_missing_mu"],
            "new_points_executed": session["new_points_executed"],
            "ready_points_recalculated": 0,
            "parallel_workers_used_in_closing_stage": 0,
            "global_restarts": 0,
            "latest_canonical_audit": audit,
            "new_point_qualifications": qualifications,
            "initial_complete_point_count": session[
                "initial_reused_parameter_points"
            ],
            "initial_reused_point_count": session[
                "initial_reused_parameter_points"
            ],
            "inherited_workers_drained": 0,
            "inherited_points_allowed_to_finish": {
                MODEL_OLD: [],
                MODEL_RLB: [],
            },
        }
    )
    if inherited_workers_drained not in (None, 0):
        raise RuntimeError("No inherited worker is permitted in this closing run.")
    _atomic_write_json(path, checkpoint)
    return checkpoint


def evaluate_production_benchmark_eta(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> dict[str, Any]:
    """Apply the 45-minute gate after one successful benchmark per model."""

    target = Path(output_dir)
    checkpoint = update_closing_checkpoint(target)
    session = checkpoint[PRODUCTION_CLOSING_SESSION_KEY]
    benchmark_by_model: dict[str, Mapping[str, Any]] = {}
    for record in session["benchmark_records"]:
        if record.get("status") == "PASS":
            benchmark_by_model[str(record["model"])] = record
    required = {MODEL_OLD, MODEL_RLB}
    if set(benchmark_by_model) != required:
        raise RuntimeError(
            f"ETA requires successful benchmarks for {sorted(required)}; got "
            f"{sorted(benchmark_by_model)}."
        )
    missing = session["current_missing_mu"]
    model_scaled = sum(
        len(missing[model])
        * float(
            benchmark_by_model[model].get(
                "benchmark_effective_wall_time_seconds",
                benchmark_by_model[model]["wall_time_seconds"],
            )
        )
        for model in required
    )
    all_benchmark_times = [
        float(
            benchmark_by_model[model].get(
                "benchmark_effective_wall_time_seconds",
                benchmark_by_model[model]["wall_time_seconds"],
            )
        )
        for model in required
    ]
    conservative = max(
        model_scaled,
        sum(len(missing[model]) for model in required)
        * max(all_benchmark_times),
    )
    gate = {
        "evaluated_utc": _utc_now(),
        "limit_seconds": PRODUCTION_ETA_LIMIT_SECONDS,
        "remaining_missing_mu": missing,
        "model_scaled_eta_seconds": model_scaled,
        "conservative_eta_seconds": conservative,
        "benchmark_wall_time_seconds": {
            model: float(
                benchmark_by_model[model].get(
                    "benchmark_effective_wall_time_seconds",
                    benchmark_by_model[model]["wall_time_seconds"],
                )
            )
            for model in required
        },
        "status": (
            "PASS" if conservative <= PRODUCTION_ETA_LIMIT_SECONDS else "FAIL"
        ),
    }
    session["eta_gate"] = gate
    checkpoint[PRODUCTION_CLOSING_SESSION_KEY] = session
    _atomic_write_json(target / CLOSING_CHECKPOINT_FILENAME, checkpoint)
    return gate


def adjudicate_saved_weak_rlb_mu080(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> dict[str, Any]:
    """Adjudicate the saved agreement-only weak-RLB benchmark locally."""

    bad_limits = {
        name: os.environ.get(name)
        for name, expected in CLOSING_THREAD_LIMITS.items()
        if os.environ.get(name) != expected
    }
    if bad_limits:
        raise RuntimeError(f"Closing thread limits are not fixed at one: {bad_limits}.")
    target_dir = Path(output_dir)
    checkpoint = _initialize_production_closing_session(
        target_dir, allow_adjudication_script_revision=True
    )
    session = checkpoint[PRODUCTION_CLOSING_SESSION_KEY]
    existing = _canonical_mu_rows(target_dir, MODEL_RLB)
    if 0.8 in _complete_mu_values(existing, MODEL_RLB):
        raise RuntimeError("weak-RLB mu=0.80 is already a complete production group.")
    failed = [
        record
        for record in session["benchmark_records"]
        if record.get("model") == MODEL_RLB
        and math.isclose(
            float(record.get("mu", math.nan)),
            0.8,
            rel_tol=0.0,
            abs_tol=5.0e-13,
        )
        and record.get("status") == "FAIL"
    ]
    if len(failed) != 1:
        raise RuntimeError(
            "Expected exactly one saved failed weak-RLB mu=0.80 benchmark."
        )
    failed_record = failed[0]
    attempts = list(failed_record["attempts"])
    if len(attempts) != 2:
        raise RuntimeError("The saved benchmark attempt history is incomplete.")
    initial_attempt = attempts[0]
    if not bool(
        int(initial_attempt["root_count"]) == OUTPUT_GUARD_POSITION
        and int(initial_attempt["unresolved_candidates_below_root9"]) == 0
        and initial_attempt["guard_not_at_scan_boundary"]
        and all(row["root_status"] == "PASS" for row in initial_attempt["roots"])
        and float(initial_attempt["maximum_primary_verification_relative"])
        > INVENTORY_VERIFICATION_TOLERANCE
    ):
        raise RuntimeError("Saved benchmark is not an agreement-only failure.")

    started = time.perf_counter()
    geometry = geometry_for(0.8, MU_TAU, MU_BETA_DEG)
    objects = build_model_objects(geometry)
    provider, provider_metadata = make_matrix_provider(
        MODEL_RLB, geometry, objects, cache=ArmPairCache.empty()
    )
    evaluation_counter = 0

    def counted_provider(omega: float) -> FloatArray:
        nonlocal evaluation_counter
        evaluation_counter += 1
        return provider(float(omega))

    policy, policy_record = _bounded_root9_policy(
        rlb2c.rlb2b.frozen_root_policy(),
        float(initial_attempt["requested_Omega_max"]),
    )
    initial_diagnostics, pairs = _reconstruct_saved_attempt_diagnostics(
        counted_provider, policy, initial_attempt
    )
    adjudication, canonical = _adjudicate_reconstructed_pairs(
        counted_provider,
        policy,
        initial_diagnostics,
        pairs,
        unresolved_candidates=int(
            initial_attempt["unresolved_candidates_below_root9"]
        ),
    )
    local_wall = time.perf_counter() - started
    diagnostic_path = target_dir / WEAK_RLB_MU080_ADJUDICATION_FILENAME
    diagnostic_payload = {
        "schema_version": 1,
        "stage": STAGE_ID,
        "model": MODEL_RLB,
        "mu": 0.8,
        "diagnostic_role": "NARROW_LOCAL_ADJUDICATION_OF_SAVED_ROOT",
        "saved_checkpoint_path": CLOSING_CHECKPOINT_FILENAME,
        "saved_benchmark_record_index": 0,
        "saved_attempt_index": 0,
        "expanded_attempt_preserved_but_not_used_for_adjudication": True,
        "initial_diagnostics": initial_diagnostics,
        "local_adjudication": adjudication,
        "policy_record": policy_record,
        "provider_setup_metadata": provider_metadata,
        "local_matrix_evaluations": evaluation_counter,
        "local_wall_time_seconds": local_wall,
        "peak_RSS_bytes": _peak_rss_bytes(),
        "production_csv_changed": False,
    }
    _atomic_write_json(diagnostic_path, diagnostic_payload)
    if adjudication["local_adjudication_status"] != "PASS":
        raise RuntimeError(
            "Local weak-RLB adjudication failed; production CSV was not changed."
        )

    point_rows = _adjudicated_point_rows(
        MODEL_RLB, geometry, pairs, canonical, adjudication
    )
    fieldnames = [
        str(value)
        for value in session["initial_mu_csv_fieldnames"][MODEL_RLB]
    ]
    point_rows = [
        {field: row.get(field, "") for field in fieldnames}
        for row in point_rows
    ]
    _complete_mu_values(point_rows, MODEL_RLB)
    qualification = _point_qualification(point_rows)
    if not (
        qualification["PLOTTED_FIRST_8"] == "PASS"
        and qualification["ROOT_9_GUARD"] == "PASS"
    ):
        raise RuntimeError("Adjudicated weak-RLB rows failed the export gate.")

    target_csv = target_dir / ROOT_FILENAMES[(SWEEP_MU, MODEL_RLB)]
    csv_hash_before = rlb2c.rlb2b._sha256(target_csv)
    preexisting_hashes = _mu_group_hashes(existing, MODEL_RLB, fieldnames)
    merged = [*existing, *point_rows]
    _complete_mu_values(merged, MODEL_RLB)
    merged.sort(key=lambda row: (float(row["mu"]), int(row["sorted_position"])))
    _atomic_merge_mu_csv_preserving_existing(target_csv, existing, merged)
    merged_after = _canonical_mu_rows(target_dir, MODEL_RLB)
    after_hashes = _mu_group_hashes(merged_after, MODEL_RLB, fieldnames)
    if any(after_hashes.get(key) != value for key, value in preexisting_hashes.items()):
        raise RuntimeError("A pre-existing weak-RLB group changed during append.")

    diagnostic_payload["production_csv_changed"] = True
    diagnostic_payload["production_csv_sha256_before"] = csv_hash_before
    diagnostic_payload["production_csv_sha256_after"] = rlb2c.rlb2b._sha256(
        target_csv
    )
    diagnostic_payload["preexisting_groups_preserved"] = True
    _atomic_write_json(diagnostic_path, diagnostic_payload)
    record = {
        "recorded_utc": _utc_now(),
        "model": MODEL_RLB,
        "mu": 0.8,
        "point_index": 40,
        "production_benchmark": True,
        "adjudicates_saved_failed_benchmark": True,
        "saved_failed_benchmark_record_index": 0,
        "initial_primary_verification_status": "FAIL",
        "local_adjudication_status": "PASS",
        "final_point_status": "LOCAL_ADJUDICATION_PASS",
        "status": "PASS",
        "prediction": failed_record["prediction"],
        "initial_requested_Omega_max": failed_record[
            "initial_requested_Omega_max"
        ],
        "initial_effective_Omega_max": failed_record[
            "initial_effective_Omega_max"
        ],
        "final_requested_Omega_max": failed_record[
            "initial_requested_Omega_max"
        ],
        "final_effective_Omega_max": failed_record[
            "initial_effective_Omega_max"
        ],
        "historical_range_expansion_preserved": True,
        "range_expansion_used_in_adjudication": False,
        "range_expansion_used": False,
        "attempt_count": 0,
        "attempts": [],
        "local_adjudication_diagnostic": diagnostic_path.name,
        "wall_time_seconds": local_wall,
        "historical_bounded_scan_wall_time_seconds": float(
            failed_record["wall_time_seconds"]
        ),
        "benchmark_effective_wall_time_seconds": float(
            failed_record["wall_time_seconds"]
        )
        + local_wall,
        "peak_RSS_bytes": _peak_rss_bytes(),
        "determinant_evaluations": evaluation_counter,
        "historical_determinant_evaluations": int(
            failed_record["determinant_evaluations"]
        ),
        "determinant_scans": 0,
        "sigma_scans": 0,
        "roots_above_9_requested_or_exported": False,
        "canonical_csv_sha256_after": rlb2c.rlb2b._sha256(target_csv),
        "root_rows": [
            {
                "sorted_position": int(row["sorted_position"]),
                "omega": float(row["omega"]),
                "Omega": float(row["Omega"]),
                "Lambda": float(row["Lambda"]),
                "scaled_sigma_ratio": float(row["scaled_sigma_ratio"]),
                "boundary_null_residual": float(
                    row["boundary_null_residual"]
                ),
                "root_status": str(row["root_status"]),
            }
            for row in point_rows
        ],
    }
    checkpoint = update_closing_checkpoint(target_dir, point_record=record)
    return {
        "model": MODEL_RLB,
        "mu": 0.8,
        "offending_root_index": initial_diagnostics[
            "maximum_difference_root_index"
        ],
        "comparison_quantity": "Omega",
        "local_adjudication_status": "PASS",
        "production_group_written": True,
        "diagnostic_path": str(diagnostic_path),
        "row_count": len(point_rows),
        "qualification": qualification,
        "checkpoint_new_points": checkpoint["new_points_executed"],
        "production_record": record,
    }


def qualify_saved_old_timoshenko_mu080_guard(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> dict[str, Any]:
    """Adopt the saved old-Timoshenko root 9 as a completeness guard."""

    bad_limits = {
        name: os.environ.get(name)
        for name, expected in CLOSING_THREAD_LIMITS.items()
        if os.environ.get(name) != expected
    }
    if bad_limits:
        raise RuntimeError(f"Closing thread limits are not fixed at one: {bad_limits}.")
    target_dir = Path(output_dir)
    checkpoint = _initialize_production_closing_session(
        target_dir, allow_root9_guard_revision=True
    )
    session = checkpoint[PRODUCTION_CLOSING_SESSION_KEY]
    existing = _canonical_mu_rows(target_dir, MODEL_OLD)
    already_complete = 0.8 in _complete_mu_values(existing, MODEL_OLD)
    failed = [
        record
        for record in session["benchmark_records"]
        if record.get("model") == MODEL_OLD
        and math.isclose(
            float(record.get("mu", math.nan)),
            0.8,
            rel_tol=0.0,
            abs_tol=5.0e-13,
        )
        and record.get("status") == "FAIL"
        and record.get("failure_reason")
        == "LOCAL_ADJUDICATION_DID_NOT_CONVERGE"
    ]
    if len(failed) != 1:
        raise RuntimeError(
            "Expected exactly one saved old-Timoshenko mu=0.80 FAIL."
        )
    failed_record = failed[0]
    attempts = list(failed_record["attempts"])
    if len(attempts) != 1:
        raise RuntimeError("The saved old-Timoshenko benchmark is ambiguous.")
    attempt = attempts[0]
    diagnostic_path = target_dir / str(
        failed_record["local_adjudication_diagnostic"]
    )
    diagnostic = json.loads(diagnostic_path.read_text(encoding="utf-8"))
    if not (
        diagnostic.get("model") == MODEL_OLD
        and math.isclose(
            float(diagnostic.get("mu", math.nan)),
            0.8,
            rel_tol=0.0,
            abs_tol=5.0e-13,
        )
        and diagnostic["local_adjudication"][
            "initial_primary_verification_status"
        ]
        == "FAIL"
        and int(attempt["root_count"]) == OUTPUT_GUARD_POSITION
        and int(attempt["unresolved_candidates_below_root9"]) == 0
        and bool(attempt["guard_not_at_scan_boundary"])
    ):
        raise RuntimeError("Saved old-Timoshenko evidence is not the expected FAIL.")

    policy, policy_record = _bounded_root9_policy(
        rlb2c.rlb2b.frozen_root_policy(),
        float(attempt["requested_Omega_max"]),
    )
    initial = diagnostic["initial_diagnostics"]
    adjudication = diagnostic["local_adjudication"]
    guard_contract = _root9_guard_contract(
        initial, adjudication, attempt, policy
    )
    if guard_contract["root9_guard_status"] != ROOT9_GUARD_INTERVAL_PASS:
        raise RuntimeError("The saved old-Timoshenko root-9 guard did not pass.")
    refined_by_position = {
        int(row["root_index"]): row
        for row in adjudication["refined_root_diagnostics"]
    }

    if already_complete:
        recovery_started = time.perf_counter()
        if diagnostic.get("root9_guard_contract") != guard_contract:
            raise RuntimeError("The interrupted guard diagnostic changed.")
        point_rows = [
            row
            for row in existing
            if math.isclose(
                float(row["mu"]), 0.8, rel_tol=0.0, abs_tol=5.0e-13
            )
        ]
        fieldnames = [
            str(value)
            for value in session["initial_mu_csv_fieldnames"][MODEL_OLD]
        ]
        _complete_mu_values(point_rows, MODEL_OLD)
        certificate = _build_root9_guard_certificate(
            target_dir,
            MODEL_OLD,
            0.8,
            point_rows,
            fieldnames,
            diagnostic_path,
            guard_contract,
        )
        qualification = _point_qualification(
            point_rows, root9_guard_certificate=certificate
        )
        if not (
            qualification["PLOTTED_FIRST_8"] == "PASS"
            and qualification["ROOT_9_GUARD"] == "PASS"
        ):
            raise RuntimeError("Interrupted guard transaction is not qualified.")
        recovery_wall = time.perf_counter() - recovery_started
        historical_wall = float(failed_record["wall_time_seconds"])
        pointwise_evaluations = int(
            diagnostic["pointwise_rehydration"]["matrix_evaluations"]
        )
        target_csv = target_dir / ROOT_FILENAMES[(SWEEP_MU, MODEL_OLD)]
        record = {
            "recorded_utc": _utc_now(),
            "model": MODEL_OLD,
            "mu": 0.8,
            "point_index": 40,
            "production_benchmark": True,
            "qualifies_saved_failed_benchmark": True,
            "recovers_atomic_csv_before_checkpoint_transaction": True,
            "saved_failed_benchmark_record_index": 2,
            "initial_primary_verification_status": "FAIL",
            "local_adjudication_status": "FAIL",
            "root9_strict_agreement_status": "FAIL",
            "root9_guard_status": ROOT9_GUARD_INTERVAL_PASS,
            "final_point_status": POINT_PASS_WITH_GUARD_QUALIFICATION,
            "status": "PASS",
            "prediction": failed_record["prediction"],
            "initial_requested_Omega_max": failed_record[
                "initial_requested_Omega_max"
            ],
            "initial_effective_Omega_max": failed_record[
                "initial_effective_Omega_max"
            ],
            "final_requested_Omega_max": failed_record[
                "final_requested_Omega_max"
            ],
            "final_effective_Omega_max": failed_record[
                "final_effective_Omega_max"
            ],
            "range_expansion_used": False,
            "attempt_count": 0,
            "attempts": [],
            "local_adjudication_diagnostic": diagnostic_path.name,
            "wall_time_seconds": recovery_wall,
            "historical_bounded_scan_wall_time_seconds": historical_wall,
            "pointwise_rehydration_wall_time_seconds": "NOT_PERSISTED_AFTER_ATOMIC_CSV_WRITE",
            "benchmark_effective_wall_time_seconds": 2.0 * historical_wall,
            "benchmark_time_basis": "CONSERVATIVE_TWICE_HISTORICAL_BOUNDED_SCAN",
            "peak_RSS_bytes": max(
                int(failed_record["peak_RSS_bytes"]), _peak_rss_bytes()
            ),
            "determinant_evaluations": pointwise_evaluations,
            "historical_determinant_evaluations": int(
                failed_record["determinant_evaluations"]
            ),
            "effective_determinant_evaluations": int(
                failed_record["determinant_evaluations"]
            )
            + pointwise_evaluations,
            "determinant_scans": 0,
            "sigma_scans": 0,
            "roots_above_9_requested_or_exported": False,
            "canonical_csv_sha256_after": rlb2c.rlb2b._sha256(target_csv),
            "root9_guard_certificate": certificate,
            "root_rows": [
                {
                    "sorted_position": int(row["sorted_position"]),
                    "omega": float(row["omega"]),
                    "Omega": float(row["Omega"]),
                    "Lambda": float(row["Lambda"]),
                    "scaled_sigma_ratio": float(row["scaled_sigma_ratio"]),
                    "boundary_null_residual": float(
                        row["boundary_null_residual"]
                    ),
                    "root_status": str(row["root_status"]),
                }
                for row in point_rows
            ],
        }
        checkpoint = update_closing_checkpoint(target_dir, point_record=record)
        return {
            "model": MODEL_OLD,
            "mu": 0.8,
            "transaction_recovered_without_spectral_recalculation": True,
            "root9_guard_status": ROOT9_GUARD_INTERVAL_PASS,
            "root8_root9_enclosures_intersect": False,
            "production_group_written": True,
            "diagnostic_path": str(diagnostic_path),
            "new_matrix_evaluations": 0,
            "checkpoint_new_points": checkpoint["new_points_executed"],
            "qualification": qualification,
        }

    started = time.perf_counter()
    geometry = geometry_for(0.8, MU_TAU, MU_BETA_DEG)
    objects = build_model_objects(geometry)
    if constitutive_checks(geometry, objects)["status"] != "PASS":
        raise RuntimeError("The saved-point arm contract no longer passes.")
    provider, provider_metadata = make_matrix_provider(
        MODEL_OLD, geometry, objects, cache=ArmPairCache.empty()
    )
    evaluator = rlb2c.rlb2b.iso_inventory._DiagnosticEvaluator(
        provider, reference_omega_scale(), policy
    )
    initial_by_position = {
        int(row["root_index"]): row for row in initial["root_diagnostics"]
    }
    pairs: list[dict[str, Any]] = []
    canonical: dict[int, Any] = {}
    for position in range(1, OUTPUT_GUARD_POSITION + 1):
        saved_root = attempt["roots"][position - 1]
        refined = refined_by_position.get(position)
        Omega = float(
            refined["determinant_based_refined_value"]
            if refined is not None
            else saved_root["Omega"]
        )
        interval = (
            refined["primary_path"]["local_sign_bracket"]
            if refined is not None
            else [
                saved_root["bracket_left_Omega"],
                saved_root["bracket_right_Omega"],
            ]
        )
        matrix_diagnostics = evaluator.diagnostics(Omega)
        quality_pass, quality_reason = (
            rlb2c.rlb2b.iso_inventory._candidate_quality(
                matrix_diagnostics, policy
            )
        )
        if not quality_pass or matrix_diagnostics.detected_nullity != 1:
            raise RuntimeError(
                f"Pointwise root {position} rehydration failed: "
                f"{quality_reason or 'detected nullity is not one'}."
            )
        candidate = rlb2c.rlb2b.iso_inventory.RootCandidate(
            case_id="rlb2d_bounded_root9__old_timoshenko__mu_0.8000",
            builder_id="RLB2D_BOUNDED_ROOT9_OLD_TIMOSHENKO",
            scan_id="saved_benchmark_pointwise_rehydration",
            Omega=Omega,
            detection_sources=(
                "ROOT9_GUARD_DETERMINISTIC_DETERMINANT_PRIMARY_BRACKET"
                if position == OUTPUT_GUARD_POSITION
                else "LOCAL_ADJUDICATION_DETERMINANT_PRIMARY_BRACKET"
                if refined is not None
                else "SAVED_BOUNDED_PRIMARY_ROOT"
            ,),
            interval_left_Omega=float(interval[0]),
            interval_right_Omega=float(interval[1]),
            interior_minimum=False,
            diagnostics=matrix_diagnostics,
            accepted=True,
            rejection_reason="",
        )
        event = rlb2c.rlb2b.iso_inventory.RootEvent(
            event_id=f"event_{position:04d}",
            Omega=Omega,
            omega=matrix_diagnostics.omega,
            multiplicity=1,
            detected_nullity=1,
            candidate=candidate,
            cluster_center_Omega=Omega,
        )
        canonical[position] = candidate
        pairs.append(
            {
                "root_index": position,
                "persisted_primary": dict(saved_root),
                "primary_event": event,
                "initial_row": initial_by_position[position],
            }
        )
    if len(evaluator.cache) != OUTPUT_GUARD_POSITION:
        raise RuntimeError("Saved-point rehydration was not pointwise-only.")

    point_rows = _adjudicated_point_rows(
        MODEL_OLD,
        geometry,
        pairs,
        canonical,
        adjudication,
        guard_contract,
    )
    fieldnames = [
        str(value)
        for value in session["initial_mu_csv_fieldnames"][MODEL_OLD]
    ]
    point_rows = [
        {field: row.get(field, "") for field in fieldnames}
        for row in point_rows
    ]
    _complete_mu_values(point_rows, MODEL_OLD)
    if not all(row["root_status"] == "PASS" for row in point_rows):
        raise RuntimeError("A pointwise rehydrated production root failed quality.")

    target_csv = target_dir / ROOT_FILENAMES[(SWEEP_MU, MODEL_OLD)]
    csv_hash_before = rlb2c.rlb2b._sha256(target_csv)
    preexisting_hashes = _mu_group_hashes(existing, MODEL_OLD, fieldnames)
    merged = [*existing, *point_rows]
    _complete_mu_values(merged, MODEL_OLD)
    merged.sort(key=lambda row: (float(row["mu"]), int(row["sorted_position"])))
    _atomic_merge_mu_csv_preserving_existing(target_csv, existing, merged)
    merged_after = _canonical_mu_rows(target_dir, MODEL_OLD)
    after_hashes = _mu_group_hashes(merged_after, MODEL_OLD, fieldnames)
    if any(
        after_hashes.get(key) != value
        for key, value in preexisting_hashes.items()
    ):
        raise RuntimeError("A pre-existing old-Timoshenko group changed.")

    local_wall = time.perf_counter() - started
    diagnostic.update(
        {
            "root9_guard_contract": guard_contract,
            "pointwise_rehydration": {
                "global_or_local_scan_run": False,
                "matrix_evaluations": len(evaluator.cache),
                "canonical_Omega": [
                    canonical[position].Omega
                    for position in range(1, OUTPUT_GUARD_POSITION + 1)
                ],
                "provider_setup_metadata": provider_metadata,
                "policy_record": policy_record,
                "all_detected_nullities": [
                    canonical[position].diagnostics.detected_nullity
                    for position in range(1, OUTPUT_GUARD_POSITION + 1)
                ],
            },
            "production_csv_changed": True,
            "production_csv_sha256_before": csv_hash_before,
            "production_csv_sha256_after": rlb2c.rlb2b._sha256(target_csv),
            "preexisting_groups_preserved": True,
            "initial_primary_verification_status_preserved": "FAIL",
            "local_adjudication_status_preserved": "FAIL",
        }
    )
    _atomic_write_json(diagnostic_path, diagnostic)
    certificate = _build_root9_guard_certificate(
        target_dir,
        MODEL_OLD,
        0.8,
        point_rows,
        fieldnames,
        diagnostic_path,
        guard_contract,
    )
    qualification = _point_qualification(
        point_rows, root9_guard_certificate=certificate
    )
    if not (
        qualification["PLOTTED_FIRST_8"] == "PASS"
        and qualification["ROOT_9_GUARD"] == "PASS"
    ):
        raise RuntimeError("The certified old-Timoshenko group failed qualification.")

    peak_rss = _peak_rss_bytes()
    record = {
        "recorded_utc": _utc_now(),
        "model": MODEL_OLD,
        "mu": 0.8,
        "point_index": 40,
        "production_benchmark": True,
        "qualifies_saved_failed_benchmark": True,
        "saved_failed_benchmark_record_index": 2,
        "initial_primary_verification_status": "FAIL",
        "local_adjudication_status": "FAIL",
        "root9_strict_agreement_status": "FAIL",
        "root9_guard_status": ROOT9_GUARD_INTERVAL_PASS,
        "final_point_status": POINT_PASS_WITH_GUARD_QUALIFICATION,
        "status": "PASS",
        "prediction": failed_record["prediction"],
        "initial_requested_Omega_max": failed_record[
            "initial_requested_Omega_max"
        ],
        "initial_effective_Omega_max": failed_record[
            "initial_effective_Omega_max"
        ],
        "final_requested_Omega_max": failed_record["final_requested_Omega_max"],
        "final_effective_Omega_max": failed_record["final_effective_Omega_max"],
        "range_expansion_used": False,
        "attempt_count": 0,
        "attempts": [],
        "local_adjudication_diagnostic": diagnostic_path.name,
        "wall_time_seconds": local_wall,
        "historical_bounded_scan_wall_time_seconds": float(
            failed_record["wall_time_seconds"]
        ),
        "benchmark_effective_wall_time_seconds": float(
            failed_record["wall_time_seconds"]
        )
        + local_wall,
        "peak_RSS_bytes": peak_rss,
        "determinant_evaluations": len(evaluator.cache),
        "historical_determinant_evaluations": int(
            failed_record["determinant_evaluations"]
        ),
        "effective_determinant_evaluations": int(
            failed_record["determinant_evaluations"]
        )
        + len(evaluator.cache),
        "determinant_scans": 0,
        "sigma_scans": 0,
        "roots_above_9_requested_or_exported": False,
        "canonical_csv_sha256_after": rlb2c.rlb2b._sha256(target_csv),
        "root9_guard_certificate": certificate,
        "root_rows": [
            {
                "sorted_position": int(row["sorted_position"]),
                "omega": float(row["omega"]),
                "Omega": float(row["Omega"]),
                "Lambda": float(row["Lambda"]),
                "scaled_sigma_ratio": float(row["scaled_sigma_ratio"]),
                "boundary_null_residual": float(
                    row["boundary_null_residual"]
                ),
                "root_status": str(row["root_status"]),
            }
            for row in point_rows
        ],
    }
    checkpoint = update_closing_checkpoint(target_dir, point_record=record)
    return {
        "model": MODEL_OLD,
        "mu": 0.8,
        "strict_root8_relative_Omega": float(
            refined_by_position[8]["final_relative_agreement_Omega"]
        ),
        "root9_strict_relative_Omega": float(
            refined_by_position[9]["final_relative_agreement_Omega"]
        ),
        "root9_guard_status": ROOT9_GUARD_INTERVAL_PASS,
        "root8_root9_enclosures_intersect": bool(
            guard_contract["enclosures_intersect"]
        ),
        "production_group_written": True,
        "diagnostic_path": str(diagnostic_path),
        "pointwise_matrix_evaluations": len(evaluator.cache),
        "wall_time_seconds": local_wall,
        "peak_RSS_bytes": peak_rss,
        "qualification": qualification,
        "checkpoint_new_points": checkpoint["new_points_executed"],
    }


def solve_one_missing_mu_point(
    model: str,
    point_index: int,
    *,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    production_benchmark: bool = False,
) -> dict[str, Any]:
    """Solve one missing point with a bounded roots-1--9 full inventory."""

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
    checkpoint = _initialize_production_closing_session(target_dir)
    session = checkpoint[PRODUCTION_CLOSING_SESSION_KEY]
    existing = _canonical_mu_rows(target_dir, model)
    complete = _complete_mu_values(existing, model)
    if mu_value in complete:
        raise RuntimeError(
            f"Refusing to recalculate ready point model={model}, mu={mu_value:g}."
        )
    current_audit = closing_mu_audit(target_dir, include_rlb_shards=False)
    current_missing = current_audit[model]["missing_mu"]
    if production_benchmark:
        if not current_missing or mu_value != max(current_missing):
            raise RuntimeError(
                "A production benchmark must be the largest currently missing "
                f"mu for {model}: {current_missing}."
            )
    elif session["eta_gate"].get("status") != "PASS":
        raise RuntimeError(
            "Remaining production points require a passing benchmark ETA gate."
        )
    point_started = time.perf_counter()
    _contract, contract_sha256 = _worker_contract(
        target_dir, active_mu, beta_grid()
    )
    geometry = geometry_for(mu_value, MU_TAU, MU_BETA_DEG)
    objects = build_model_objects(geometry)
    if constitutive_checks(geometry, objects)["status"] != "PASS":
        raise RuntimeError(f"Arm contract failed before roots at {geometry}.")
    prediction = _predict_root9_Omega(existing, model, mu_value)
    initial_requested_max = (
        PRODUCTION_INITIAL_BOUND_FACTOR * prediction["predicted_Omega_9"]
    )
    point_rows, first_attempt, first_inventory = _run_bounded_root9_attempt(
        model,
        geometry,
        objects,
        contract_sha256,
        initial_requested_max,
    )
    attempts = [first_attempt]
    selected_inventory = first_inventory
    expansion_used = bool(
        first_attempt["disposition"] == ATTEMPT_RANGE_EXPANSION_REQUIRED
    )
    if expansion_used:
        expanded_requested_max = (
            PRODUCTION_SINGLE_EXPANSION_FACTOR
            * float(first_attempt["effective_Omega_max"])
        )
        point_rows, second_attempt, selected_inventory = (
            _run_bounded_root9_attempt(
                model,
                geometry,
                objects,
                contract_sha256,
                expanded_requested_max,
            )
        )
        attempts.append(second_attempt)
    final_attempt = attempts[-1]
    local_adjudication: dict[str, Any] | None = None
    local_adjudication_path: Path | None = None
    root9_guard: dict[str, Any] | None = None
    if (
        not point_rows
        and final_attempt["disposition"]
        == ATTEMPT_LOCAL_ADJUDICATION_REQUIRED
        and selected_inventory is not None
    ):
        selected_policy, _selected_policy_record = _bounded_root9_policy(
            rlb2c.rlb2b.frozen_root_policy(),
            float(final_attempt["requested_Omega_max"]),
        )
        initial_diagnostics, live_pairs = _pairs_from_live_inventory(
            selected_inventory, selected_policy
        )
        adjudication_provider, _adjudication_metadata = make_matrix_provider(
            model,
            geometry,
            objects,
            cache=ArmPairCache.empty(),
        )
        adjudication_result, canonical = _adjudicate_reconstructed_pairs(
            adjudication_provider,
            selected_policy,
            initial_diagnostics,
            live_pairs,
            unresolved_candidates=int(
                final_attempt["unresolved_candidates_below_root9"]
            ),
        )
        root9_guard = _root9_guard_contract(
            initial_diagnostics,
            adjudication_result,
            final_attempt,
            selected_policy,
        )
        local_adjudication = {
            "schema_version": 1,
            "stage": STAGE_ID,
            "model": model,
            "mu": mu_value,
            "diagnostic_role": "LOCAL_AGREEMENT_ADJUDICATION",
            "initial_diagnostics": initial_diagnostics,
            "local_adjudication": adjudication_result,
            "root9_guard_contract": root9_guard,
        }
        local_adjudication_path = _local_adjudication_path(
            target_dir, model, mu_value
        )
        _atomic_write_json(local_adjudication_path, local_adjudication)
        point_rows = _adjudicated_point_rows(
            model,
            geometry,
            live_pairs,
            canonical,
            adjudication_result,
            root9_guard,
        )
    point_wall = time.perf_counter() - point_started
    peak_rss = _peak_rss_bytes()
    record: dict[str, Any] = {
        "recorded_utc": _utc_now(),
        "model": model,
        "mu": mu_value,
        "point_index": point_index,
        "production_benchmark": bool(production_benchmark),
        "prediction": prediction,
        "initial_requested_Omega_max": initial_requested_max,
        "initial_effective_Omega_max": first_attempt[
            "effective_Omega_max"
        ],
        "final_requested_Omega_max": attempts[-1][
            "requested_Omega_max"
        ],
        "final_effective_Omega_max": attempts[-1][
            "effective_Omega_max"
        ],
        "range_expansion_used": expansion_used,
        "attempt_count": len(attempts),
        "attempts": attempts,
        "initial_primary_verification_status": (
            "FAIL" if local_adjudication is not None else "PASS"
        ),
        "local_adjudication_status": (
            local_adjudication["local_adjudication"][
                "local_adjudication_status"
            ]
            if local_adjudication is not None
            else "NOT_REQUIRED"
        ),
        "root9_strict_agreement_status": (
            local_adjudication["root9_guard_contract"][
                "root9_strict_agreement_status"
            ]
            if local_adjudication is not None
            else "PASS" if point_rows else "NOT_EVALUATED"
        ),
        "root9_guard_status": (
            local_adjudication["root9_guard_contract"][
                "root9_guard_status"
            ]
            if local_adjudication is not None
            else "STRICT_AGREEMENT_PASS" if point_rows else "NOT_EVALUATED"
        ),
        "final_point_status": (
            local_adjudication["root9_guard_contract"]["point_status"]
            if local_adjudication is not None
            and local_adjudication["root9_guard_contract"]["root9_guard_status"]
            == ROOT9_GUARD_INTERVAL_PASS
            else local_adjudication["local_adjudication"]["final_point_status"]
            if local_adjudication is not None
            else "DIRECT_INVENTORY_PASS" if point_rows else "FAIL"
        ),
        "local_adjudication_diagnostic": (
            str(local_adjudication_path.relative_to(target_dir))
            if local_adjudication_path is not None
            else None
        ),
        "wall_time_seconds": point_wall,
        "peak_RSS_bytes": peak_rss,
        "determinant_evaluations": sum(
            int(attempt["determinant_evaluations"]) for attempt in attempts
        ),
        "determinant_scans": sum(
            int(attempt["determinant_scans"])
            for attempt in attempts
            if attempt["determinant_scans"] is not None
        ),
        "sigma_scans": sum(
            int(attempt["sigma_scans"])
            for attempt in attempts
            if attempt["sigma_scans"] is not None
        ),
        "roots_above_9_requested_or_exported": False,
        "status": "PASS" if point_rows else "FAIL",
    }
    if not point_rows:
        record["failure_reason"] = {
            ATTEMPT_RANGE_EXPANSION_REQUIRED: (
                "BOUNDED_ROOT9_COMPLETENESS_FAILED_AFTER_SINGLE_EXPANSION"
            ),
            ATTEMPT_LOCAL_ADJUDICATION_REQUIRED: (
                "LOCAL_ADJUDICATION_DID_NOT_CONVERGE"
            ),
            ATTEMPT_HARD_FAIL: "BOUNDED_ROOT9_STRUCTURAL_OR_QUALITY_FAIL",
        }.get(final_attempt["disposition"], "BOUNDED_ROOT9_INVENTORY_FAIL")
        update_closing_checkpoint(target_dir, point_record=record)
        raise RuntimeError(
            f"Bounded production inventory failed for {model}, mu={mu_value:g}; "
            "the canonical CSV was not changed."
        )
    fieldnames = [
        str(value)
        for value in session["initial_mu_csv_fieldnames"][model]
    ]
    point_rows = [
        {field: row.get(field, "") for field in fieldnames}
        for row in point_rows
    ]
    _complete_mu_values(point_rows, model)
    provisional_qualification = _point_qualification(
        point_rows,
        root9_guard_certificate=(
            root9_guard
            if root9_guard is not None
            and root9_guard.get("root9_guard_status")
            == ROOT9_GUARD_INTERVAL_PASS
            else None
        ),
    )
    if provisional_qualification["PLOTTED_FIRST_8"] != "PASS" or (
        provisional_qualification["ROOT_9_GUARD"] != "PASS"
    ):
        record["status"] = "FAIL"
        record["failure_reason"] = "PLOTTED_OR_ROOT9_QUALIFICATION_FAILED"
        update_closing_checkpoint(target_dir, point_record=record)
        raise RuntimeError("The bounded point did not pass the exported range gate.")
    merged = [*existing, *point_rows]
    _complete_mu_values(merged, model)
    merged.sort(key=lambda row: (float(row["mu"]), int(row["sorted_position"])))
    target = target_dir / ROOT_FILENAMES[(SWEEP_MU, model)]
    _atomic_merge_mu_csv_preserving_existing(target, existing, merged)
    merged_after = _canonical_mu_rows(target_dir, model)
    _complete_mu_values(merged_after, model)
    record["canonical_csv_sha256_after"] = rlb2c.rlb2b._sha256(target)
    if (
        root9_guard is not None
        and root9_guard.get("root9_guard_status")
        == ROOT9_GUARD_INTERVAL_PASS
        and local_adjudication_path is not None
    ):
        record["root9_guard_certificate"] = _build_root9_guard_certificate(
            target_dir,
            model,
            mu_value,
            point_rows,
            fieldnames,
            local_adjudication_path,
            root9_guard,
        )
    record["root_rows"] = [
        {
            "sorted_position": int(row["sorted_position"]),
            "omega": float(row["omega"]),
            "Omega": float(row["Omega"]),
            "Lambda": float(row["Lambda"]),
            "scaled_sigma_ratio": float(row["scaled_sigma_ratio"]),
            "boundary_null_residual": float(row["boundary_null_residual"]),
            "root_status": str(row["root_status"]),
        }
        for row in point_rows
    ]
    checkpoint = update_closing_checkpoint(target_dir, point_record=record)
    qualification = _point_qualification(
        point_rows,
        root9_guard_certificate=record.get("root9_guard_certificate"),
    )
    print(
        f"[closing/{model}] mu={mu_value:g}: "
        f"export={qualification['aggregate_export_range_status']}, "
        f"internal={qualification['INTERNAL_TAIL']}, "
        f"wall={point_wall:.3f}s, peak_RSS={peak_rss}; checkpoint updated",
        flush=True,
    )
    return {
        "model": model,
        "mu": mu_value,
        "row_count": len(point_rows),
        "qualification": qualification,
        "checkpoint_new_points": checkpoint["new_points_executed"],
        "production_record": record,
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
    *,
    root9_guard_certificates: Mapping[
        tuple[str, str, float], Mapping[str, Any]
    ] | None = None,
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
    certificates = root9_guard_certificates or {}
    strict_export_axes: list[float] = []
    guard_qualified_axes: list[float] = []
    effective_export_by_axis: dict[float, bool] = {}
    for axis in expected_axis:
        point_rows = [
            row
            for row in rows
            if round(_float(row, axis_field), 12) == axis
        ]
        strict_pass = bool(
            point_rows
            and all(
                _bool(row["export_primary_verification_agreement"])
                and row["export_range_status"] == "PASS"
                for row in point_rows
            )
        )
        certificate = certificates.get((sweep, model, axis))
        certified_guard = bool(
            point_rows
            and all(
                row["export_range_status"]
                == POINT_PASS_WITH_GUARD_QUALIFICATION
                and not _bool(row["export_primary_verification_agreement"])
                for row in point_rows
            )
            and certificate is not None
            and certificate.get("strict_roots_1_to_8_status") == "PASS"
            and certificate.get("root9_strict_agreement_status") == "FAIL"
            and certificate.get("root9_guard_status")
            == ROOT9_GUARD_INTERVAL_PASS
            and certificate.get("point_status")
            == POINT_PASS_WITH_GUARD_QUALIFICATION
        )
        if strict_pass:
            strict_export_axes.append(axis)
        if certified_guard:
            guard_qualified_axes.append(axis)
        effective_export_by_axis[axis] = strict_pass or certified_guard
    strict_export_agreement = len(strict_export_axes) == len(expected_axis)
    export_agreement = all(effective_export_by_axis.values())
    export_issue_axes = [
        axis
        for axis in expected_axis
        if not effective_export_by_axis[axis]
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
        "strict_export_axis_values": strict_export_axes,
        "strict_export_axis_count": len(strict_export_axes),
        "root9_guard_qualified_axis_values": guard_qualified_axes,
        "root9_guard_qualified_axis_count": len(guard_qualified_axes),
        "tail_issue_axis_values": tail_issue_axes,
        "tail_issue_count": len(tail_issue_axes),
        "exact_row_structure_passed": exact_structure,
        "accepted_root_quality_passed": accepted_quality,
        "strict_aggregate_primary_verification_passed": (
            strict_export_agreement
        ),
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
    *,
    root9_guard_certificates: Mapping[
        tuple[str, str, float], Mapping[str, Any]
    ] | None = None,
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
                certificate = (root9_guard_certificates or {}).get(
                    (sweep, model, axis_value)
                )
                qualifications.append(
                    _point_qualification(
                        point_rows,
                        root9_guard_certificate=certificate,
                    )
                )
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
        and item["aggregate_export_range_status"]
        in {"PASS", POINT_PASS_WITH_GUARD_QUALIFICATION}
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


def _legacy_preclosing_report_text(
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
        if item["aggregate_export_range_status"]
        in {"PASS", POINT_PASS_WITH_GUARD_QUALIFICATION}
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


def _report_text(
    contract: Mapping[str, Any],
    summaries: Sequence[Mapping[str, Any]],
    statuses: Mapping[str, str],
    qualifications: Sequence[Mapping[str, Any]],
    closing_checkpoint: Mapping[str, Any],
) -> str:
    """Create the concise Russian closing report from verified data only."""

    session = closing_checkpoint[PRODUCTION_CLOSING_SESSION_KEY]
    reference = contract["reference_constants"]
    lines = [
        "# RLB-2D: зависимости безразмерной частоты от $\\mu$ и $\\beta$",
        "",
        "## Расчётный контракт",
        "",
        (
            "Рассматриваются независимо упорядоченные спектральные позиции "
            "1--8. Девятый корень используется только как контроль полноты и "
            "на графиках не показан. Продолжение форм колебаний, сглаживание "
            "и межмодельные относительные расхождения не применялись."
        ),
        "",
        (
            "Для всех моделей приняты $E_0={E0:g}$, $\\nu_0={nu0:g}$, "
            "$\\rho_0={rho0:g}$, $b_0={b0:g}$, $h_0={h0:g}$ и $l={l:g}$. "
            "Длины и толщины плеч заданы соотношениями "
            "$L_1=l(1-\\mu)$, $L_2=l(1+\\mu)$, "
            "$h_1=h_0(1-\\tau)$ и $h_2=h_0(1+\\tau)$."
        ).format(**reference),
        "",
        (
            "Модель RLB соответствует case A: $\\delta=0{,}1$, "
            "$E_1=1{,}1$, $E_2=0{,}9$, $\\nu_{12}=0{,}3$, "
            "$G_{12}=G_{13}=G_{23}=1/2{,}6$, $\\rho=1$. "
            "Каждое плечо содержит четыре равных слоя $[0/90/90/0]$; "
            "однослойное сокращение не использовалось."
        ),
        "",
        "Использована старая нормировка",
        "",
        "$$",
        "\\Omega=\\omega l^2\\sqrt{\\frac{\\rho_0A_0}{E_0I_{y0}}},",
        "\\qquad \\Lambda=\\sqrt{\\Omega},",
        "\\qquad A_0=b_0h_0,\\quad I_{y0}=\\frac{b_0h_0^3}{12}.",
        "$$",
        "",
        (
            "Первая сетка содержит 41 значение $\\mu=0;0{,}02;\\ldots;0{,}80$ "
            "при $\\beta=15^\\circ$ и $\\tau=0$. Вторая сетка содержит "
            "91 значение $\\beta=0^\\circ;1^\\circ;\\ldots;90^\\circ$ "
            "при $\\mu=0{,}5$ и $\\tau=0{,}2$."
        ),
        "",
        "## Полнота данных",
        "",
        "| сетка | модель | точек | строк | positions 1--8 | root 9 | статус |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for summary in summaries:
        lines.append(
            f"| {summary['sweep']} | {summary['model']} | "
            f"{int(summary['axis_count'])} | {int(summary['actual_rows'])} | "
            f"{int(summary['first_frequency_count'])} | "
            f"{int(summary['guard_count'])} | {summary['status']} |"
        )
    lines.extend(
        [
            "",
            "## Завершающий расчёт по $\\mu$",
            "",
            (
                "До запуска отсутствовали точки: old Timoshenko -- "
                f"`{session['initial_missing_mu'][MODEL_OLD]}`; weak RLB -- "
                f"`{session['initial_missing_mu'][MODEL_RLB]}`. Готовые "
                "группы не пересчитывались."
            ),
            "",
            (
                "Для каждой новой точки выполнены два независимых прохода "
                "существующего determinant/SVD workflow: основной и "
                "проверочный. Верхняя граница определялась по секущему "
                "прогнозу $\\Omega_9$ и не использовалась для вычисления "
                "самих частот. Вычислялись только positions 1--8 и root 9."
            ),
            "",
            "| модель | $\\mu$ | wall time, s | peak RSS, MiB | "
            "$\\Omega_{max}$ | evaluations | expansion |",
            "|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for record in session["point_records"]:
        if record.get("status") != "PASS":
            continue
        reported_wall = float(
            record.get(
                "benchmark_effective_wall_time_seconds",
                record["wall_time_seconds"],
            )
        )
        reported_evaluations = int(
            record.get(
                "effective_determinant_evaluations",
                int(record.get("historical_determinant_evaluations", 0))
                + int(record["determinant_evaluations"]),
            )
        )
        lines.append(
            f"| {record['model']} | {float(record['mu']):.2f} | "
            f"{reported_wall:.3f} | "
            f"{int(record['peak_RSS_bytes']) / 2**20:.1f} | "
            f"{float(record['final_effective_Omega_max']):.6f} | "
            f"{reported_evaluations} | "
            f"{record['range_expansion_used']} |"
        )
    range_qualifications = [
        item
        for item in qualifications
        if item["PLOTTED_FIRST_8"] != "PASS"
        or item["ROOT_9_GUARD"] != "PASS"
    ]
    lines.extend(["", "## Численные квалификации", ""])
    if range_qualifications:
        lines.append(
            "В ранее сохранённых inventories остаются следующие агрегированные "
            "primary/verification qualifications; согласно контракту эти точки "
            "не пересчитывались:"
        )
        lines.append("")
        for item in range_qualifications:
            lines.append(
                f"- `{item['sweep']}/{item['model']}`, "
                f"{item['parameter_name']}={float(item['parameter_value']):g}: "
                f"first 8 = `{item['PLOTTED_FIRST_8']}`, root 9 = "
                f"`{item['ROOT_9_GUARD']}`, aggregate relative = "
                f"{float(item['aggregate_primary_verification_max_relative']):.6e}."
            )
    else:
        lines.append("Все группы positions 1--8 и root 9 прошли проверки.")
    guard_records = [
        record
        for record in session["point_records"]
        if isinstance(record.get("root9_guard_certificate"), Mapping)
    ]
    lines.extend(["", "### Strict roots 1--8 и root-9 guards", ""])
    lines.append(
        "Для positions 1--8 сохранён строгий относительный порог по "
        f"$\\Omega$, равный {INVENTORY_VERIFICATION_TOLERANCE:.0e}. "
        "Root 9 используется только для контроля полноты и на графиках "
        "не показан."
    )
    if guard_records:
        for record in guard_records:
            certificate = record["root9_guard_certificate"]
            lines.append(
                f"- `{record['model']}`, $\\mu={float(record['mu']):.2f}$: "
                "strict roots 1--8 = `PASS`; root-9 strict agreement = "
                f"`FAIL`; root-9 guard = `{certificate['root9_guard_status']}`; "
                f"enclosure gap = {float(certificate['enclosure_gap_Omega']):.12g}; "
                f"diagnostic `{certificate['diagnostic_path']}` "
                f"(SHA-256 `{certificate['diagnostic_sha256']}`)."
            )
    else:
        lines.append("Отдельная interval-квалификация root 9 не потребовалась.")
    lines.extend(
        [
            "",
            (
                "Квалификации внутреннего хвоста выше root 9 не уточнялись, "
                "поскольку такой расчёт не входит в данный этап."
            ),
            "",
            "## Графики",
            "",
            (
                "Файлы [lambda_vs_mu_plot.png](lambda_vs_mu_plot.png) и "
                "[lambda_vs_beta_plot.png](lambda_vs_beta_plot.png) содержат "
                "по 24 кривые. Цвет обозначает отсортированную позицию 1--8, "
                "а тип линии -- модель. Готовый график по $\\beta$ повторно "
                "не строился; его SHA-256 проверен до и после финализации."
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
                "Результат относится только к двум объявленным сеткам и case A "
                "при $\\delta=0{,}1$. Линии соединяют независимо "
                "упорядоченные частоты и не определяют модальные ветви. "
                "Расчёты Ritz и FEM, MAC, сглаживание и branch tracking не "
                "выполнялись. Inventories roots 10--13 не запрашивались и не "
                "экспортировались; возможные detector events выше root 9 "
                "сохраняются только как численная provenance диапазона."
            ),
            "",
        ]
    )
    return "\n".join(lines)


def finalize_output(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    mu_values: Sequence[float] | None = None,
    beta_values: Sequence[float] | None = None,
) -> dict[str, Any]:
    """Validate six root files, reuse beta PNG, and create only the mu PNG."""

    finalization_started = time.perf_counter()

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
    beta_plot = target / PLOT_FILENAMES[SWEEP_BETA]
    if not beta_plot.is_file() or beta_plot.stat().st_size == 0:
        raise RuntimeError("The completed beta plot is missing.")
    beta_plot_hash_before = rlb2c.rlb2b._sha256(beta_plot)

    rows_by_pair = {
        (sweep, model): rlb2c.rlb2b._read_csv(
            target / ROOT_FILENAMES[(sweep, model)]
        )
        for sweep in SWEEPS
        for model in MODELS
    }
    checkpoint = update_closing_checkpoint(target, inherited_workers_drained=0)
    guard_certificates = _root9_guard_certificate_index(checkpoint, target)
    summaries = [
        validate_root_rows(
            sweep,
            model,
            rows_by_pair[(sweep, model)],
            active_mu if sweep == SWEEP_MU else active_beta,
            root9_guard_certificates=guard_certificates,
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
    closing_session = checkpoint[PRODUCTION_CLOSING_SESSION_KEY]
    if beta_plot_hash_before != closing_session["beta_plot_sha256_before"]:
        raise RuntimeError("The beta plot differs from its launch hash.")
    beta_csv_hashes_after = {
        model: rlb2c.rlb2b._sha256(
            target / ROOT_FILENAMES[(SWEEP_BETA, model)]
        )
        for model in MODELS
    }
    if beta_csv_hashes_after != closing_session["beta_csv_sha256"]:
        raise RuntimeError("A frozen beta CSV differs from its launch hash.")
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
    beta_plot_hash_after = rlb2c.rlb2b._sha256(beta_plot)
    if beta_plot_hash_after != beta_plot_hash_before:
        raise RuntimeError("Finalization changed the frozen beta plot.")
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
    qualifications = collect_point_qualifications(
        rows_by_pair,
        root9_guard_certificates=guard_certificates,
    )
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
        "root9_plotted": False,
        "summary": {
            "root_data": summaries,
            "all_export_ranges_passed": all_export_pass,
            "all_first8_strict_or_frozen_passed": all(
                item["PLOTTED_FIRST_8"] == "PASS"
                for item in qualifications
            ),
            "all_root9_strict_or_guard_passed": all(
                item["ROOT_9_GUARD"] == "PASS"
                for item in qualifications
            ),
            "exact_mu_grid_passed": exact_mu_grid,
            "exact_beta_grid_passed": exact_beta_grid,
        },
        "closing_stage": {
            "new_points_executed": checkpoint["new_points_executed"],
            "reused_points": checkpoint["initial_reused_point_count"],
            "ready_points_recalculated": 0,
            "parallel_workers_used_in_closing_stage": 0,
            "inherited_workers_drained": 0,
            "inherited_points_allowed_to_finish": checkpoint[
                "inherited_points_allowed_to_finish"
            ],
            "global_restarts": 0,
            "thread_limits": CLOSING_THREAD_LIMITS,
            "missing_mu_before": checkpoint[PRODUCTION_CLOSING_SESSION_KEY][
                "initial_missing_mu"
            ],
            "missing_mu_after": checkpoint[PRODUCTION_CLOSING_SESSION_KEY][
                "current_missing_mu"
            ],
            "computed_points": checkpoint[PRODUCTION_CLOSING_SESSION_KEY][
                "computed_points"
            ],
            "production_benchmarks": checkpoint[
                PRODUCTION_CLOSING_SESSION_KEY
            ]["benchmark_records"],
            "production_point_records": checkpoint[
                PRODUCTION_CLOSING_SESSION_KEY
            ]["point_records"],
            "initial_agreement_failures": [
                record
                for record in checkpoint[PRODUCTION_CLOSING_SESSION_KEY][
                    "point_records"
                ]
                if record.get("initial_primary_verification_status") == "FAIL"
            ],
            "local_adjudication_points": [
                record
                for record in checkpoint[PRODUCTION_CLOSING_SESSION_KEY][
                    "point_records"
                ]
                if record.get("local_adjudication_status") in {"PASS", "FAIL"}
            ],
            "root9_guard_qualified_points": [
                record["root9_guard_certificate"]
                for record in checkpoint[PRODUCTION_CLOSING_SESSION_KEY][
                    "point_records"
                ]
                if isinstance(record.get("root9_guard_certificate"), Mapping)
            ],
            "root9_guard_contract": {
                "strict_positions": list(
                    range(1, PLOTTED_POSITIONS + 1)
                ),
                "strict_relative_Omega_threshold": (
                    INVENTORY_VERIFICATION_TOLERANCE
                ),
                "guard_position": OUTPUT_GUARD_POSITION,
                "guard_may_replace_strict_frequency_agreement": False,
                "guard_role": "COMPLETENESS_AND_ORDERING_ONLY",
                "requires_disjoint_ordered_root8_root9_enclosures": True,
                "requires_zero_unresolved_below_root9": True,
                "requires_frozen_root_quality_thresholds": True,
                "root9_plotted": False,
            },
            "benchmark_eta_gate": checkpoint[
                PRODUCTION_CLOSING_SESSION_KEY
            ]["eta_gate"],
            "bounded_FULL_scan_points": checkpoint[
                PRODUCTION_CLOSING_SESSION_KEY
            ]["new_points_executed"],
            "FAST_LOCAL_point_count": 0,
            "full_fallback_point_count": 0,
            "roots_above_9_requested_or_exported": False,
            "cumulative_point_wall_time_seconds": checkpoint[
                PRODUCTION_CLOSING_SESSION_KEY
            ]["cumulative_point_wall_time_seconds"],
            "finalization_wall_time_seconds": time.perf_counter()
            - finalization_started,
            "peak_RSS_bytes": max(
                int(
                    checkpoint[PRODUCTION_CLOSING_SESSION_KEY][
                        "peak_RSS_bytes"
                    ]
                ),
                _peak_rss_bytes(),
            ),
            "preexisting_groups_preserved": checkpoint[
                PRODUCTION_CLOSING_SESSION_KEY
            ]["preexisting_groups_preserved"],
            "beta_plot_reused_without_redraw": True,
            "beta_plot_sha256_before": beta_plot_hash_before,
            "beta_plot_sha256_after": beta_plot_hash_after,
            "beta_csv_sha256_before": closing_session[
                "beta_csv_sha256"
            ],
            "beta_csv_sha256_after": beta_csv_hashes_after,
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
        "figures_created_in_closing_stage": 1,
        "figures_available": 2,
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
    run_manifest["closing_stage"]["total_closing_compute_wall_time_seconds"] = (
        float(
            run_manifest["closing_stage"][
                "cumulative_point_wall_time_seconds"
            ]
        )
        + float(run_manifest["closing_stage"]["finalization_wall_time_seconds"])
    )
    run_manifest["closing_stage"]["quality_thresholds"] = {
        "root_singular_ratio": rlb2c.rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE,
        "boundary_null_residual": rlb2c.rlb2b.BOUNDARY_RESIDUAL_TOLERANCE,
        "primary_verification_relative": INVENTORY_VERIFICATION_TOLERANCE,
    }
    run_manifest["closing_stage"]["rlb2d_script_sha256"] = checkpoint[
        PRODUCTION_CLOSING_SESSION_KEY
    ]["rlb2d_script_sha256"]
    run_manifest["closing_stage"]["production_physics_hashes"] = {
        "before": checkpoint[PRODUCTION_CLOSING_SESSION_KEY][
            "production_physics_hashes_before"
        ],
        "after": frozen_after,
    }
    run_manifest["closing_stage"]["protected_closing_hashes"] = {
        "before": checkpoint[PRODUCTION_CLOSING_SESSION_KEY][
            "protected_closing_hashes_before"
        ],
        "after": _protected_closing_hashes(),
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
            "closing_adjudicate_weak_rlb_mu080",
            "closing_qualify_old_timoshenko_mu080_guard",
            "closing_benchmark_eta",
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
    parser.add_argument(
        "--production-benchmark",
        action="store_true",
        help=(
            "Mark the largest missing one-point closing calculation as one "
            "of the two production benchmarks."
        ),
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
            MODEL_OLD,
            args.point_start,
            output_dir=output_dir,
            production_benchmark=args.production_benchmark,
        )
    elif args.mode == "closing_mu_new_rlb":
        solve_one_missing_mu_point(
            MODEL_RLB,
            args.point_start,
            output_dir=output_dir,
            production_benchmark=args.production_benchmark,
        )
    elif args.mode == "closing_adjudicate_weak_rlb_mu080":
        print(
            json.dumps(
                adjudicate_saved_weak_rlb_mu080(output_dir),
                indent=2,
                sort_keys=True,
            )
        )
    elif args.mode == "closing_qualify_old_timoshenko_mu080_guard":
        print(
            json.dumps(
                qualify_saved_old_timoshenko_mu080_guard(output_dir),
                indent=2,
                sort_keys=True,
            )
        )
    elif args.mode == "closing_benchmark_eta":
        print(
            json.dumps(
                evaluate_production_benchmark_eta(output_dir),
                indent=2,
                sort_keys=True,
            )
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
