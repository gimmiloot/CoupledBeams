"""RLB-2F signed one-arm layered-contrast frequency map.

The analysis reuses the frozen RLB-2E determinant/SVD search machinery but
defines a separate, narrow physical-case and I/O contract.  Arm 1 contains
four equal zero-degree plies with signed outer/inner stiffness scaling; arm 2
is the homogeneous four-ply M0 reference.  Independently sorted positions
1--8 are plotted and position 9 is retained only as a completeness guard.

``--plot-only`` imports no production-physics module and performs no root
calculation.  Predictors are numerical locators only; every exported root is
refined from the frozen characteristic matrix.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Any, Mapping, Sequence


for _thread_variable in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import numpy as np
from numpy.typing import NDArray


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import sweep_grid_policy  # noqa: E402
from scripts.analysis.laminated_beams import (  # noqa: E402
    sweep_reddy_stiffness_layout_contrast as rlb2e,
)


FloatArray = NDArray[np.float64]

STAGE_ID = "RLB-2F"
ALGORITHM_VERSION = "one_arm_layered_signed_contrast_fast_plot_v1"
POLICY_ID = "frequency-map-v1"
CONFIGURATION_ID = "ONE_ARM_LAYERED"

OUTER_PLY_MATERIAL = "OUTER_PLY_MATERIAL"
INNER_PLY_MATERIAL = "INNER_PLY_MATERIAL"
HOMOGENEOUS_M0 = "HOMOGENEOUS_M0"
LAYERED_LAYOUT = (
    OUTER_PLY_MATERIAL,
    INNER_PLY_MATERIAL,
    INNER_PLY_MATERIAL,
    OUTER_PLY_MATERIAL,
)
HOMOGENEOUS_LAYOUT = (HOMOGENEOUS_M0,) * 4
PLY_REGIONS = ("OUTER", "INNER", "INNER", "OUTER")

E1_0 = rlb2e.E1_0
E2_0 = rlb2e.E2_0
NU12_0 = rlb2e.NU12_0
G0 = rlb2e.G0
RHO_0 = rlb2e.RHO_0
DELTA = rlb2e.DELTA

MU = rlb2e.MU
TAU = rlb2e.TAU
BETA_DEG = rlb2e.BETA_DEG
L_REFERENCE = rlb2e.L_REFERENCE
L1 = rlb2e.L1
L2 = rlb2e.L2
L_TOTAL = rlb2e.L_TOTAL
WIDTH = rlb2e.WIDTH
THICKNESS = rlb2e.THICKNESS
PLY_THICKNESS = rlb2e.PLY_THICKNESS
K = rlb2e.K

XI_MIN = -0.8
XI_MAX = 0.8
XI_STEP = 0.02
LOCAL_XI_STEP = rlb2e.LOCAL_CHI_STEP
K_PLOT = rlb2e.K_PLOT
K_GUARD = rlb2e.K_GUARD

OMEGA_TO_OMEGA_SCALE = rlb2e.OMEGA_TO_OMEGA_SCALE
MATRIX_RELATIVE_TOLERANCE = rlb2e.MATRIX_RELATIVE_TOLERANCE
SYMMETRY_RELATIVE_TOLERANCE = rlb2e.SYMMETRY_RELATIVE_TOLERANCE
REDUCED_PROPERTY_TOLERANCE = rlb2e.REDUCED_PROPERTY_TOLERANCE
REDUCTION_ROUTE_TOLERANCE = rlb2e.REDUCTION_ROUTE_TOLERANCE
ROOT_SINGULAR_RATIO_TOLERANCE = rlb2e.ROOT_SINGULAR_RATIO_TOLERANCE
BOUNDARY_RESIDUAL_TOLERANCE = rlb2e.BOUNDARY_RESIDUAL_TOLERANCE
ARM_SWAP_RELATIVE_TOLERANCE = rlb2e.ARM_SWAP_RELATIVE_TOLERANCE
ROOT9_RIGHT_TAIL_OMEGA = rlb2e.ROOT9_RIGHT_TAIL_OMEGA
ETA_LIMIT_SECONDS = rlb2e.ETA_LIMIT_SECONDS
NEIGHBOUR_MAD_MULTIPLIER = rlb2e.NEIGHBOUR_MAD_MULTIPLIER
NEIGHBOUR_ABSOLUTE_TRIGGER = rlb2e.NEIGHBOUR_ABSOLUTE_TRIGGER

DEFAULT_OUTPUT_DIR = (
    ROOT / "results" / "laminated_beams" / "reddy_one_arm_layered_contrast_sweep"
)
SPECTRUM_FILENAME = "spectrum_roots.csv"
SECTION_FILENAME = "section_properties.csv"
AUDIT_FILENAME = "neighbour_audit.csv"
BENCHMARK_FILENAME = "benchmark.json"
CHECKPOINT_FILENAME = "checkpoint.json"
ARM_SWAP_FILENAME = "arm_swap_diagnostics.json"
PLOT_FILENAME = "lambda_vs_xi_one_arm_layered.png"
REPORT_FILENAME = "report.md"
MANIFEST_FILENAME = "run_manifest.json"
MANDATORY_COMPLETED_OUTPUTS = frozenset(
    {
        SPECTRUM_FILENAME,
        SECTION_FILENAME,
        AUDIT_FILENAME,
        BENCHMARK_FILENAME,
        CHECKPOINT_FILENAME,
        ARM_SWAP_FILENAME,
        PLOT_FILENAME,
        REPORT_FILENAME,
    }
)

PRODUCTION_PHYSICS_PATHS = rlb2e.PRODUCTION_PHYSICS_PATHS
RLB2E_SCRIPT_PATH = (
    "scripts/analysis/laminated_beams/sweep_reddy_stiffness_layout_contrast.py"
)
RLB2E_RESULT_DIR = (
    ROOT / "results" / "laminated_beams" / "reddy_stiffness_layout_contrast_sweep"
)

FREQUENCY_MAP_POLICY = {
    "frequency_map_policy": POLICY_ID,
    "calculation_mode": "fast_plot",
    "spectrum_semantics": "sorted_positions",
    "sweep_parameter": "xi",
    "parameter_grid": "-0.80:0.02:0.80",
    "continuation_anchor": "0.00",
    "continuation_paths": ["0.00:-0.02:-0.80", "0.00:0.02:0.80"],
    "K_plot": K_PLOT,
    "K_guard": K_GUARD,
    "guard_root_role": "completeness_only",
    "neighbour_audit": "enabled",
    "local_repair_policy": "triggered_only",
    "strict_audit_default": False,
    "branch_tracking": False,
    "mac": False,
    "mode_shapes": False,
    "energy_analysis": False,
}

SPECTRUM_FIELDS = (
    "row_id",
    "configuration_id",
    "xi",
    "grid_kind",
    "continuation_leg",
    "solve_id",
    "transaction_id",
    "sorted_position",
    "root_role",
    "guard_flag",
    "omega",
    "Omega",
    "Lambda",
    "predictor_Omega",
    "predictor_used_as_final",
    "locator_interval_left_Omega",
    "locator_interval_right_Omega",
    "root_interval_left_Omega",
    "root_interval_right_Omega",
    "detector_refiner_provenance",
    "raw_determinant",
    "scaled_determinant",
    "raw_sigma_ratio",
    "scaled_sigma_ratio",
    "boundary_null_residual",
    "detected_nullity",
    "unresolved_candidates_below_root9",
    "search_right_Omega",
    "root9_right_margin_Omega",
    "solve_mode",
    "fallback_used",
    "quality_status",
    "point_status",
    "is_canonical_plot_source",
    "supersedes_row_id",
    "repair_id",
    "roots_above_9_computed",
)


@dataclass(frozen=True)
class PointSolution:
    xi: float
    rows: tuple[dict[str, Any], ...]
    wall_time_seconds: float
    peak_rss_bytes: int
    determinant_evaluations: int
    sigma_evaluations: int
    search_left_Omega: float
    search_right_Omega: float
    local_refinements: int
    solve_mode: str
    fallback_used: bool
    unresolved_candidates_below_root9: int
    continuation_leg: str
    swapped_arms: bool = False
    candidate_count_total: int = K_GUARD
    accepted_candidates_above_root9: int = 0
    retained_slots_above_root9: int = 0
    roots_above_9_computed: bool = False


ROOT_CALCULATION_COUNT = 0


def xi_grid() -> FloatArray:
    """Return the exact 81-point signed grid."""

    return np.asarray(
        sweep_grid_policy.parameter_grid(XI_MIN, XI_MAX, XI_STEP), dtype=float
    )


def continuation_paths() -> tuple[FloatArray, FloatArray]:
    """Return the declared zero-to-negative and zero-to-positive legs."""

    negative = np.asarray(
        sweep_grid_policy.parameter_grid(0.0, abs(XI_MIN), XI_STEP), dtype=float
    )
    negative = -negative
    negative[0] = 0.0
    positive = np.asarray(
        sweep_grid_policy.parameter_grid(0.0, XI_MAX, XI_STEP), dtype=float
    )
    return negative, positive


def omega_to_Omega(omega: float) -> float:
    return rlb2e.omega_to_Omega(omega)


def Omega_to_Lambda(Omega: float) -> float:
    return rlb2e.Omega_to_Lambda(Omega)


def base_material_contract() -> dict[str, float]:
    return rlb2e.base_material_contract()


def signed_materials(xi: float) -> dict[str, Any]:
    """Return signed outer/inner materials and the unchanged M0 reference."""

    value = float(xi)
    if not math.isfinite(value) or not XI_MIN <= value <= XI_MAX:
        raise ValueError("xi must lie in [-0.8, 0.8].")
    contrast = rlb2e.contrast_materials(abs(value))
    outer_key, inner_key = (("H", "L") if value >= 0.0 else ("L", "H"))
    reference = rlb2e.contrast_materials(0.0)["H"]
    return {
        OUTER_PLY_MATERIAL: contrast[outer_key],
        INNER_PLY_MATERIAL: contrast[inner_key],
        HOMOGENEOUS_M0: reference,
    }


def material_scales(xi: float) -> dict[str, float]:
    value = float(xi)
    if not XI_MIN <= value <= XI_MAX:
        raise ValueError("xi must lie in [-0.8, 0.8].")
    return {
        OUTER_PLY_MATERIAL: 1.0 + value,
        INNER_PLY_MATERIAL: 1.0 - value,
        HOMOGENEOUS_M0: 1.0,
    }


def build_layered_section(xi: float) -> Any:
    beam, _coupled = rlb2e._physics_modules()
    materials = signed_materials(xi)
    laminate = beam.integrate_laminate(
        [
            beam.Ply(
                materials[label],
                angle_deg=0.0,
                thickness=PLY_THICKNESS,
                label=label,
            )
            for label in LAYERED_LAYOUT
        ]
    )
    properties = beam.reduce_to_beam_properties(
        laminate,
        width=WIDTH,
        K=K,
        symmetry_tolerance=SYMMETRY_RELATIVE_TOLERANCE,
        reduction_tolerance=REDUCTION_ROUTE_TOLERANCE,
    )
    return rlb2e.SectionObjects(LAYERED_LAYOUT, laminate, properties)


def build_homogeneous_section() -> Any:
    baseline = rlb2e._baseline_section()
    return rlb2e.SectionObjects(HOMOGENEOUS_LAYOUT, baseline.laminate, baseline.properties)


def build_arm_sections(xi: float, *, swapped: bool = False) -> tuple[Any, Any]:
    layered = build_layered_section(xi)
    homogeneous = build_homogeneous_section()
    return (homogeneous, layered) if swapped else (layered, homogeneous)


def constitutive_gate(
    values: Sequence[float] = (-0.8, -0.4, 0.0, 0.4, 0.8),
) -> dict[str, Any]:
    """Verify the signed analytic section contract before any root search."""

    baseline = build_homogeneous_section()
    checks: list[dict[str, Any]] = []
    maxima = {
        "D_matrix_formula_relative": 0.0,
        "Dbeam_formula_relative": 0.0,
        "A_matrix_relative": 0.0,
        "shear_matrix_relative": 0.0,
        "reduced_invariant_relative": 0.0,
        "symmetry_relative": 0.0,
        "reduction_route_relative": 0.0,
    }
    full_grid_positive = True
    for xi in xi_grid():
        materials = signed_materials(float(xi))
        full_grid_positive = full_grid_positive and all(
            getattr(material, name) > 0.0
            for material in materials.values()
            for name in ("E1", "E2", "G12", "G13", "G23")
        )
    passed = full_grid_positive
    for raw_xi in values:
        xi = float(raw_xi)
        layered = build_layered_section(xi)
        expected = 1.0 + 0.75 * xi
        matrix_residual = rlb2e._matrix_relative(
            layered.laminate.D, expected * baseline.laminate.D
        )
        beam_residual = rlb2e._relative(
            layered.properties.D / baseline.properties.D, expected
        )
        A_residual = rlb2e._matrix_relative(
            layered.laminate.A, baseline.laminate.A
        )
        shear_residual = rlb2e._matrix_relative(
            layered.laminate.shear, baseline.laminate.shear
        )
        mass_residual = max(
            rlb2e._relative(layered.laminate.I0, baseline.laminate.I0),
            rlb2e._relative(layered.laminate.I2, baseline.laminate.I2),
        )
        reduced_residual = max(
            rlb2e._relative(
                getattr(layered.properties, name),
                getattr(baseline.properties, name),
            )
            for name in ("A", "S", "m", "J")
        )
        symmetry_residual = max(
            rlb2e._scaled_B(layered.laminate),
            rlb2e._scaled_B(baseline.laminate),
            rlb2e._scaled_I1(layered.laminate),
            rlb2e._scaled_I1(baseline.laminate),
        )
        route_residual = max(
            rlb2e._reduction_max_residual(layered.properties),
            rlb2e._reduction_max_residual(baseline.properties),
        )
        point_pass = bool(
            matrix_residual <= MATRIX_RELATIVE_TOLERANCE
            and beam_residual <= REDUCED_PROPERTY_TOLERANCE
            and A_residual <= MATRIX_RELATIVE_TOLERANCE
            and shear_residual <= MATRIX_RELATIVE_TOLERANCE
            and mass_residual <= MATRIX_RELATIVE_TOLERANCE
            and reduced_residual <= REDUCED_PROPERTY_TOLERANCE
            and symmetry_residual <= SYMMETRY_RELATIVE_TOLERANCE
            and route_residual <= REDUCTION_ROUTE_TOLERANCE
        )
        passed = passed and point_pass
        maxima["D_matrix_formula_relative"] = max(
            maxima["D_matrix_formula_relative"], matrix_residual
        )
        maxima["Dbeam_formula_relative"] = max(
            maxima["Dbeam_formula_relative"], beam_residual
        )
        maxima["A_matrix_relative"] = max(maxima["A_matrix_relative"], A_residual)
        maxima["shear_matrix_relative"] = max(
            maxima["shear_matrix_relative"], shear_residual
        )
        maxima["reduced_invariant_relative"] = max(
            maxima["reduced_invariant_relative"], mass_residual, reduced_residual
        )
        maxima["symmetry_relative"] = max(
            maxima["symmetry_relative"], symmetry_residual
        )
        maxima["reduction_route_relative"] = max(
            maxima["reduction_route_relative"], route_residual
        )
        checks.append(
            {
                "xi": xi,
                "status": "PASS" if point_pass else "FAIL",
                "D_layered_over_D0": layered.properties.D / baseline.properties.D,
                "expected_D_layered_over_D0": expected,
                "D_total_over_D0": 1.0 + layered.properties.D / baseline.properties.D,
                "expected_D_total_over_D0": 2.0 + 0.75 * xi,
                "D_matrix_formula_relative": matrix_residual,
                "Dbeam_formula_relative": beam_residual,
                "A_matrix_relative": A_residual,
                "shear_matrix_relative": shear_residual,
                "mass_relative": mass_residual,
                "reduced_invariant_relative": reduced_residual,
                "symmetry_relative": symmetry_residual,
                "reduction_route_relative": route_residual,
            }
        )
    return {
        "status": "PASS" if passed else "FAIL",
        "passed": bool(passed),
        "full_grid_moduli_positive": bool(full_grid_positive),
        "checks": checks,
        "maximum_residuals": maxima,
        "D0_matrix": baseline.laminate.D,
        "Dbeam0": baseline.properties.D,
        "Abeam0": baseline.properties.A,
        "Sbeam0": baseline.properties.S,
        "m0": baseline.properties.m,
        "J0": baseline.properties.J,
        "tolerances": {
            "matrix_relative": MATRIX_RELATIVE_TOLERANCE,
            "symmetry_relative": SYMMETRY_RELATIVE_TOLERANCE,
            "reduced_property_relative": REDUCED_PROPERTY_TOLERANCE,
            "reduction_route_relative": REDUCTION_ROUTE_TOLERANCE,
        },
    }


def section_property_rows() -> list[dict[str, Any]]:
    """Return 81 signed values times two arm rows."""

    baseline = build_homogeneous_section()
    rows: list[dict[str, Any]] = []
    for xi_value in xi_grid():
        xi = float(xi_value)
        scales = material_scales(xi)
        for arm_id, arm_role, section in (
            (1, "LAYERED", build_layered_section(xi)),
            (2, "HOMOGENEOUS_REFERENCE", baseline),
        ):
            layered = arm_role == "LAYERED"
            expected = 1.0 + 0.75 * xi if layered else 1.0
            ply_scales = (
                [
                    scales[OUTER_PLY_MATERIAL],
                    scales[INNER_PLY_MATERIAL],
                    scales[INNER_PLY_MATERIAL],
                    scales[OUTER_PLY_MATERIAL],
                ]
                if layered
                else [1.0] * 4
            )
            props = section.properties
            laminate = section.laminate
            actual = props.D / baseline.properties.D
            D_matrix_residual = rlb2e._matrix_relative(
                laminate.D, expected * baseline.laminate.D
            )
            A_matrix_residual = rlb2e._matrix_relative(
                laminate.A, baseline.laminate.A
            )
            shear_matrix_residual = rlb2e._matrix_relative(
                laminate.shear, baseline.laminate.shear
            )
            B_relative = rlb2e._scaled_B(laminate)
            I1_relative = rlb2e._scaled_I1(laminate)
            D_formula_residual = rlb2e._relative(actual, expected)
            invariant_residual = max(
                rlb2e._relative(props.A, baseline.properties.A),
                rlb2e._relative(props.S, baseline.properties.S),
                rlb2e._relative(props.m, baseline.properties.m),
                rlb2e._relative(props.J, baseline.properties.J),
                rlb2e._relative(laminate.I0, baseline.laminate.I0),
                rlb2e._relative(laminate.I2, baseline.laminate.I2),
            )
            reduction_residual = rlb2e._reduction_max_residual(props)
            constitutive_pass = bool(
                D_matrix_residual <= MATRIX_RELATIVE_TOLERANCE
                and A_matrix_residual <= MATRIX_RELATIVE_TOLERANCE
                and shear_matrix_residual <= MATRIX_RELATIVE_TOLERANCE
                and B_relative <= SYMMETRY_RELATIVE_TOLERANCE
                and I1_relative <= SYMMETRY_RELATIVE_TOLERANCE
                and D_formula_residual <= REDUCED_PROPERTY_TOLERANCE
                and invariant_residual <= REDUCED_PROPERTY_TOLERANCE
                and reduction_residual <= REDUCTION_ROUTE_TOLERANCE
            )
            rows.append(
                {
                    "xi": xi,
                    "arm_id": arm_id,
                    "arm_role": arm_role,
                    "stack_bottom_to_top": list(section.layout),
                    "ply_index": [1, 2, 3, 4],
                    "ply_region": list(PLY_REGIONS),
                    "material_scale": ply_scales,
                    "angle_deg": [0.0] * 4,
                    "ply_thickness": [PLY_THICKNESS] * 4,
                    "z_interfaces": laminate.z_interfaces,
                    "A_matrix": laminate.A,
                    "B_matrix": laminate.B,
                    "D_matrix": laminate.D,
                    "shear_matrix_yz_xz": laminate.shear,
                    "I0": laminate.I0,
                    "I1": laminate.I1,
                    "I2": laminate.I2,
                    "A_beam": props.A,
                    "D_beam": props.D,
                    "S_beam": props.S,
                    "m": props.m,
                    "J": props.J,
                    "B_relative": B_relative,
                    "I1_relative": I1_relative,
                    "D_beam_over_D0": actual,
                    "expected_D_beam_over_D0": expected,
                    "D_matrix_formula_residual": D_matrix_residual,
                    "A_matrix_invariance_residual": A_matrix_residual,
                    "shear_matrix_invariance_residual": shear_matrix_residual,
                    "D_formula_residual": D_formula_residual,
                    "beam_mass_invariance_residual": invariant_residual,
                    "outer_material_scale": (
                        scales[OUTER_PLY_MATERIAL] if layered else 1.0
                    ),
                    "inner_material_scale": (
                        scales[INNER_PLY_MATERIAL] if layered else 1.0
                    ),
                    "reduction_route_max_relative": reduction_residual,
                    "constitutive_status": "PASS" if constitutive_pass else "FAIL",
                }
            )
    return rows


def make_matrix_provider(xi: float, *, swapped: bool = False) -> tuple[Any, dict[str, Any]]:
    """Build the frozen two-arm RLB boundary matrix for one signed point."""

    _beam, coupled = rlb2e._physics_modules()
    arm1, arm2 = build_arm_sections(xi, swapped=swapped)
    joint = np.asarray(coupled.joint_matrix(math.radians(BETA_DEG)), dtype=float)
    identical = all(
        rlb2e._relative(
            getattr(arm1.properties, name), getattr(arm2.properties, name)
        )
        <= 8.0 * np.finfo(float).eps
        for name in ("A", "D", "S", "m", "J")
    )

    def provider(omega: float) -> FloatArray:
        map1 = np.asarray(
            coupled.arm_clamp_map(float(omega), L1, arm1.properties), dtype=float
        )
        map2 = map1
        if not identical:
            map2 = np.asarray(
                coupled.arm_clamp_map(float(omega), L2, arm2.properties), dtype=float
            )
        combined = np.zeros((12, 6), dtype=float)
        combined[:6, :3] = map1
        combined[6:, 3:] = map2
        return np.asarray(joint @ combined, dtype=float)

    direct_residual = 0.0
    for probe in (0.731, 3.217):
        direct = coupled.coupled_boundary_matrix(
            probe,
            math.radians(BETA_DEG),
            L1,
            arm1.properties,
            L2,
            arm2.properties,
        )
        direct_residual = max(
            direct_residual,
            float(np.max(np.abs(provider(probe) - np.asarray(direct)))),
        )
    if direct_residual > 16.0 * np.finfo(float).eps:
        raise RuntimeError("RLB-2F provider differs from the frozen public builder.")
    return provider, {
        "xi": float(xi),
        "swapped_arms": bool(swapped),
        "beta_deg": BETA_DEG,
        "cached_vs_public_builder_max_abs": direct_residual,
        "production_modules_only": True,
    }


def _root_rows(
    xi: float,
    slots: Sequence[Any],
    *,
    solve_id: str,
    transaction_id: str,
    solve_mode: str,
    fallback_used: bool,
    predicted: Sequence[float] | None,
    search_right: float,
    unresolved: int,
    continuation_leg: str,
    guard_locator_right_width: float = 1.2,
    grid_kind: str = "BASE",
    repair_id: str = "",
) -> tuple[dict[str, Any], ...]:
    windows = (
        None
        if predicted is None
        else rlb2e.local_search_windows(
            predicted, guard_right_width=guard_locator_right_width
        )
    )
    guard = float(slots[K_GUARD - 1].event.omega_bar)
    rows: list[dict[str, Any]] = []
    for position, slot in enumerate(slots[:K_GUARD], start=1):
        event = slot.event
        candidate = event.candidate
        diagnostic = candidate.diagnostics
        Omega = float(event.omega_bar)
        locator = (
            (float(candidate.interval_left_bar), float(candidate.interval_right_bar))
            if windows is None
            else windows[position - 1]
        )
        row_id = (
            f"{CONFIGURATION_ID}__{float(xi):+.6f}__{grid_kind}__"
            f"p{position:02d}__{repair_id or 'base'}"
        )
        rows.append(
            {
                "row_id": row_id,
                "configuration_id": CONFIGURATION_ID,
                "xi": float(xi),
                "grid_kind": grid_kind,
                "continuation_leg": continuation_leg,
                "solve_id": solve_id,
                "transaction_id": transaction_id,
                "sorted_position": position,
                "root_role": "PLOTTED" if position <= K_PLOT else "ROOT_9_GUARD",
                "guard_flag": position == K_GUARD,
                "omega": float(event.omega),
                "Omega": Omega,
                "Lambda": Omega_to_Lambda(Omega),
                "predictor_Omega": (
                    "" if predicted is None else float(predicted[position - 1])
                ),
                "predictor_used_as_final": False,
                "locator_interval_left_Omega": locator[0],
                "locator_interval_right_Omega": locator[1],
                "root_interval_left_Omega": candidate.interval_left_bar,
                "root_interval_right_Omega": candidate.interval_right_bar,
                "detector_refiner_provenance": candidate.detection_sources,
                "raw_determinant": diagnostic.raw_determinant,
                "scaled_determinant": diagnostic.scaled_determinant,
                "raw_sigma_ratio": diagnostic.raw_sigma_ratio,
                "scaled_sigma_ratio": diagnostic.scaled_sigma_ratio,
                "boundary_null_residual": diagnostic.raw_boundary_null_residual,
                "detected_nullity": diagnostic.detected_nullity,
                "unresolved_candidates_below_root9": unresolved,
                "search_right_Omega": search_right,
                "root9_right_margin_Omega": search_right - guard,
                "solve_mode": solve_mode,
                "fallback_used": fallback_used,
                "quality_status": "PASS",
                "point_status": "PASS",
                "is_canonical_plot_source": True,
                "supersedes_row_id": "",
                "repair_id": repair_id,
                "roots_above_9_computed": False,
            }
        )
    return tuple(rows)


def _guard_detector_duplicate(slots: Sequence[Any], policy: Any) -> bool:
    """Identify two overlapping estimates of one root-9 neighbourhood.

    This is not a frequency matching tolerance and does not merge a close
    cluster.  It only triggers a fresh dense local SVD/determinant refinement
    when the inherited primary/verification detector reconciliation tolerance,
    overlapping enclosures, and unit nullities all identify one neighbourhood.
    """

    if len(slots) != K_GUARD + 1:
        return False
    left_slot, right_slot = slots[K_GUARD - 1 : K_GUARD + 1]
    left = left_slot.event.candidate
    right = right_slot.event.candidate
    left_value = float(left_slot.event.omega_bar)
    right_value = float(right_slot.event.omega_bar)
    tolerance = float(policy.reference_detector_reconciliation_atol_bar) + float(
        policy.reference_detector_reconciliation_rtol
    ) * max(abs(left_value), abs(right_value), 1.0)
    enclosures_overlap = bool(
        max(float(left.interval_left_bar), float(right.interval_left_bar))
        <= min(float(left.interval_right_bar), float(right.interval_right_bar))
    )
    unit_nullities = bool(
        int(left.diagnostics.detected_nullity) == 1
        and int(right.diagnostics.detected_nullity) == 1
    )
    return bool(
        0.0 < right_value - left_value <= tolerance
        and enclosures_overlap
        and unit_nullities
    )


def _dense_guard_reconciliation(
    counted: Any,
    policy: Any,
    slots: Sequence[Any],
    *,
    solve_id: str,
) -> tuple[list[Any], list[Any], list[Any], float, int, FloatArray]:
    """Repeat one identified root-9 neighbourhood with the frozen local tools."""

    predicted = np.asarray(
        [float(slot.event.omega_bar) for slot in slots[:K_GUARD]], dtype=float
    )
    canonical, rejected, refined_slots, search_right, refinements = (
        rlb2e._local_candidates(
            counted,
            policy,
            predicted,
            solve_id=solve_id + "_guard_reconciliation",
            dense=True,
            dense_positions=(K_GUARD,),
            guard_right_width=0.2,
        )
    )
    return (
        canonical,
        rejected,
        refined_slots,
        search_right,
        refinements,
        predicted,
    )


def solve_point(
    xi: float,
    *,
    previous: tuple[float, Sequence[float]] | None = None,
    second_previous: tuple[float, Sequence[float]] | None = None,
    force_anchor: bool = False,
    dense_local: bool = False,
    dense_positions: Sequence[int] | None = None,
    grid_kind: str = "BASE",
    repair_id: str = "",
    continuation_leg: str = "ANCHOR",
    swapped: bool = False,
    guard_locator_right_width: float = 1.2,
) -> PointSolution:
    """Solve one signed point from the frozen characteristic matrix."""

    global ROOT_CALCULATION_COUNT
    ROOT_CALCULATION_COUNT += 1
    started = time.perf_counter()
    solve_id = (
        f"{CONFIGURATION_ID}__xi_{float(xi):+.6f}"
        + ("__ARM_SWAP" if swapped else "")
    )
    transaction_id = hashlib.sha256(
        f"{STAGE_ID}|{solve_id}|{grid_kind}|{repair_id}".encode("utf-8")
    ).hexdigest()[:20].upper()
    provider, _metadata = make_matrix_provider(xi, swapped=swapped)
    counted = rlb2e.CountedProvider(provider)
    policy = rlb2e._rlb2e_search_policy()
    predicted: FloatArray | None = None
    fallback_used = False
    used_guard_locator_right_width = guard_locator_right_width
    solve_mode = ""
    if not force_anchor and previous is not None:
        predicted = rlb2e.hold_secant_predictor(
            xi,
            previous[0],
            previous[1],
            None if second_previous is None else second_previous[0],
            None if second_previous is None else second_previous[1],
        )
        predicted = np.sort(predicted)
        try:
            canonical, rejected, slots, search_right, refinements = (
                rlb2e._local_candidates(
                    counted,
                    policy,
                    predicted,
                    solve_id=solve_id,
                    dense=dense_local,
                    dense_positions=dense_positions,
                    guard_right_width=guard_locator_right_width,
                )
            )
            accepted, quality = rlb2e._point_is_acceptable(
                canonical, rejected, slots, search_right, policy
            )
        except (ValueError, RuntimeError, ArithmeticError, np.linalg.LinAlgError):
            accepted = False
            quality = {}
            canonical, rejected, slots = [], [], []
            search_right = math.nan
            refinements = 0
        if not accepted and _guard_detector_duplicate(slots, policy):
            (
                canonical,
                rejected,
                slots,
                search_right,
                extra_refinements,
                predicted,
            ) = _dense_guard_reconciliation(
                counted, policy, slots, solve_id=solve_id
            )
            used_guard_locator_right_width = 0.2
            refinements += extra_refinements
            accepted, quality = rlb2e._point_is_acceptable(
                canonical, rejected, slots, search_right, policy
            )
            if accepted:
                solve_mode = "FAST_LOCAL_DENSE_GUARD_RECONCILIATION"
        if not accepted:
            fallback_used = True
            canonical, rejected, slots, search_right, refinements = (
                rlb2e._progressive_anchor_candidates(
                    counted, policy, solve_id=solve_id + "_fallback"
                )
            )
            accepted, quality = rlb2e._point_is_acceptable(
                canonical, rejected, slots, search_right, policy
            )
            solve_mode = "BOUNDED_GLOBAL_RECOVERY"
            if not accepted and _guard_detector_duplicate(slots, policy):
                (
                    canonical,
                    rejected,
                    slots,
                    search_right,
                    extra_refinements,
                    predicted,
                ) = _dense_guard_reconciliation(
                    counted, policy, slots, solve_id=solve_id + "_fallback"
                )
                used_guard_locator_right_width = 0.2
                refinements += extra_refinements
                accepted, quality = rlb2e._point_is_acceptable(
                    canonical, rejected, slots, search_right, policy
                )
                if accepted:
                    solve_mode = "GLOBAL_RECOVERY_DENSE_GUARD_RECONCILIATION"
        elif not solve_mode:
            solve_mode = "FAST_LOCAL"
    else:
        canonical, rejected, slots, search_right, refinements = (
            rlb2e._progressive_anchor_candidates(counted, policy, solve_id=solve_id)
        )
        accepted, quality = rlb2e._point_is_acceptable(
            canonical, rejected, slots, search_right, policy
        )
        solve_mode = "PROGRESSIVE_ANCHOR"
        if not accepted and _guard_detector_duplicate(slots, policy):
            # The two inherited detector phases have produced overlapping
            # estimates of the same unit-nullity guard neighbourhood.  Repeat
            # roots 1...9 with the existing local detector and a dense root-9
            # SVD window.  No tenth-root window is scanned or retained.
            (
                canonical,
                rejected,
                slots,
                search_right,
                local_refinements,
                predicted,
            ) = _dense_guard_reconciliation(
                counted, policy, slots, solve_id=solve_id
            )
            used_guard_locator_right_width = 0.2
            refinements += local_refinements
            accepted, quality = rlb2e._point_is_acceptable(
                canonical, rejected, slots, search_right, policy
            )
            if accepted:
                solve_mode = "PROGRESSIVE_ANCHOR_DENSE_GUARD_RECONCILIATION"
    if not accepted:
        raise RuntimeError(f"{solve_id}: first-eight-plus-root9 failed: {quality}")
    exported_prediction = (
        predicted
        if solve_mode
        in {
            "FAST_LOCAL",
            "FAST_LOCAL_DENSE_GUARD_RECONCILIATION",
            "GLOBAL_RECOVERY_DENSE_GUARD_RECONCILIATION",
            "PROGRESSIVE_ANCHOR_DENSE_GUARD_RECONCILIATION",
        }
        else None
    )
    rows = _root_rows(
        xi,
        slots,
        solve_id=solve_id,
        transaction_id=transaction_id,
        solve_mode=solve_mode,
        fallback_used=fallback_used,
        predicted=exported_prediction,
        search_right=search_right,
        unresolved=int(quality["unresolved_candidates_below_root9"]),
        continuation_leg=continuation_leg,
        guard_locator_right_width=used_guard_locator_right_width,
        grid_kind=grid_kind,
        repair_id=repair_id,
    )
    return PointSolution(
        xi=float(xi),
        rows=rows,
        wall_time_seconds=time.perf_counter() - started,
        peak_rss_bytes=rlb2e._peak_rss_bytes(),
        determinant_evaluations=counted.evaluations,
        sigma_evaluations=counted.evaluations,
        search_left_Omega=(
            1.0e-8
            if exported_prediction is None
            else rlb2e.local_search_windows(exported_prediction)[0][0]
        ),
        search_right_Omega=search_right,
        local_refinements=refinements,
        solve_mode=solve_mode,
        fallback_used=fallback_used,
        unresolved_candidates_below_root9=int(
            quality["unresolved_candidates_below_root9"]
        ),
        continuation_leg=continuation_leg,
        swapped_arms=swapped,
        candidate_count_total=int(quality["candidate_count_total"]),
        accepted_candidates_above_root9=int(
            quality["accepted_candidates_above_root9"]
        ),
        retained_slots_above_root9=int(quality["retained_slots_above_root9"]),
        roots_above_9_computed=bool(quality["roots_above_9_computed"]),
    )


def _as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _base_group_index(
    rows: Sequence[Mapping[str, Any]],
) -> dict[float, list[Mapping[str, Any]]]:
    groups: dict[float, list[Mapping[str, Any]]] = {}
    for row in rows:
        if str(row.get("grid_kind")) != "BASE":
            continue
        groups.setdefault(round(float(row["xi"]), 10), []).append(row)
    return groups


def _base_group_has_exact_positions(group: Sequence[Mapping[str, Any]]) -> bool:
    try:
        positions = [int(row["sorted_position"]) for row in group]
    except (KeyError, TypeError, ValueError):
        return False
    return len(positions) == K_GUARD and sorted(positions) == list(
        range(1, K_GUARD + 1)
    )


def _base_group_is_acceptable(group: Sequence[Mapping[str, Any]]) -> bool:
    if not _base_group_has_exact_positions(group):
        return False
    try:
        ordered = sorted(group, key=lambda row: int(row["sorted_position"]))
        Omegas = np.asarray([float(row["Omega"]) for row in ordered], dtype=float)
        omegas = np.asarray([float(row["omega"]) for row in ordered], dtype=float)
        Lambdas = np.asarray([float(row["Lambda"]) for row in ordered], dtype=float)
        roles_ok = all(
            str(row["root_role"])
            == ("PLOTTED" if position <= K_PLOT else "ROOT_9_GUARD")
            and _as_bool(row["guard_flag"]) == (position == K_GUARD)
            for position, row in enumerate(ordered, start=1)
        )
        numerical_ok = bool(
            np.all(np.isfinite(Omegas))
            and np.all(np.isfinite(omegas))
            and np.all(np.isfinite(Lambdas))
            and np.all(Omegas > 0.0)
            and np.all(omegas > 0.0)
            and np.all(Lambdas > 0.0)
            and np.all(np.diff(Omegas) > 0.0)
            and np.allclose(
                Omegas,
                omegas * OMEGA_TO_OMEGA_SCALE,
                rtol=2.0e-12,
                atol=2.0e-12,
            )
            and np.allclose(
                Lambdas * Lambdas,
                Omegas,
                rtol=2.0e-12,
                atol=2.0e-12,
            )
        )
        quality_ok = all(
            str(row["quality_status"]) == "PASS"
            and int(row["unresolved_candidates_below_root9"]) == 0
            and float(row["scaled_sigma_ratio"])
            <= ROOT_SINGULAR_RATIO_TOLERANCE
            and float(row["boundary_null_residual"])
            <= BOUNDARY_RESIDUAL_TOLERANCE
            and not _as_bool(row["predictor_used_as_final"])
            and not _as_bool(row["roots_above_9_computed"])
            for row in ordered
        )
        guard_ok = (
            float(ordered[-1]["root9_right_margin_Omega"])
            >= ROOT9_RIGHT_TAIL_OMEGA - 1.0e-10
        )
    except (KeyError, TypeError, ValueError, OverflowError):
        return False
    return bool(roles_ok and numerical_ok and quality_ok and guard_ok)


def _complete_base_group_index(
    rows: Sequence[Mapping[str, Any]],
) -> dict[float, list[Mapping[str, Any]]]:
    return {
        key: group
        for key, group in _base_group_index(rows).items()
        if _base_group_is_acceptable(group)
    }


def _canonical_group(
    rows: Sequence[Mapping[str, Any]], xi: float
) -> list[Mapping[str, Any]]:
    selected = [
        row
        for row in rows
        if round(float(row["xi"]), 10) == round(float(xi), 10)
        and _as_bool(row.get("is_canonical_plot_source", True))
    ]
    selected.sort(key=lambda row: int(row["sorted_position"]))
    if [int(row["sorted_position"]) for row in selected] != list(
        range(1, K_GUARD + 1)
    ):
        raise RuntimeError(f"Incomplete canonical group at xi={xi:+.2f}.")
    return selected


def _rows_for_roots(rows: Sequence[Mapping[str, Any]], xi: float) -> FloatArray:
    return np.asarray(
        [float(row["Omega"]) for row in _canonical_group(rows, xi)], dtype=float
    )


def audit_spectrum_rows(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    groups = _base_group_index(rows)
    complete = _complete_base_group_index(rows)
    expected = {round(float(value), 10) for value in xi_grid()}
    duplicates: list[str] = []
    incomplete: list[str] = []
    quality_failures: list[str] = []
    for xi, group in groups.items():
        positions = [int(row["sorted_position"]) for row in group]
        if len(positions) != len(set(positions)):
            duplicates.append(f"{xi:+.2f}")
        if not _base_group_has_exact_positions(group):
            incomplete.append(f"{xi:+.2f}")
        elif not _base_group_is_acceptable(group):
            quality_failures.append(f"{xi:+.2f}")
    base = [row for row in rows if str(row.get("grid_kind")) == "BASE"]
    above = [row for row in rows if int(row["sorted_position"]) > K_GUARD]
    duplicate_row_ids = len(rows) - len({str(row["row_id"]) for row in rows})
    canonical_counts: dict[tuple[float, int], int] = {}
    for row in rows:
        if _as_bool(row.get("is_canonical_plot_source", False)):
            key = (round(float(row["xi"]), 10), int(row["sorted_position"]))
            canonical_counts[key] = canonical_counts.get(key, 0) + 1
    expected_canonical = {
        (round(float(xi), 10), position)
        for xi in xi_grid()
        for position in range(1, K_GUARD + 1)
    }
    canonical_failures = [
        f"{key[0]:+.2f}:p{key[1]}"
        for key in sorted(expected_canonical | set(canonical_counts))
        if canonical_counts.get(key, 0) != 1
    ]
    missing = sorted(expected - set(groups))
    extra = sorted(set(groups) - expected)
    passed = bool(
        not duplicates
        and not incomplete
        and not quality_failures
        and not missing
        and not extra
        and not above
        and duplicate_row_ids == 0
        and not canonical_failures
        and len(base) == len(xi_grid()) * K_GUARD
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "base_group_count": len(complete),
        "base_row_count": len(base),
        "missing_groups": [f"{value:+.2f}" for value in missing],
        "extra_groups": [f"{value:+.2f}" for value in extra],
        "duplicate_groups": duplicates,
        "incomplete_groups": incomplete,
        "physical_quality_failures": quality_failures,
        "duplicate_row_id_count": duplicate_row_ids,
        "canonical_source_failures": canonical_failures,
        "roots_above_guard_count": len(above),
    }


def canonical_plot_rows(
    rows: Sequence[Mapping[str, Any]],
) -> list[Mapping[str, Any]]:
    return [
        row
        for row in rows
        if int(row["sorted_position"]) <= K_PLOT
        and _as_bool(row.get("is_canonical_plot_source", True))
    ]


def audit_plot_rows(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    selected = canonical_plot_rows(rows)
    expected = {
        (round(float(xi), 10), position)
        for xi in xi_grid()
        for position in range(1, K_PLOT + 1)
    }
    counts: dict[tuple[float, int], int] = {}
    invalid: list[str] = []
    for row in selected:
        key = (round(float(row["xi"]), 10), int(row["sorted_position"]))
        counts[key] = counts.get(key, 0) + 1
        value = float(row["Lambda"])
        gap = math.isnan(value) and str(row.get("point_status", "")) == (
            "UNRESOLVED_AFTER_LOCAL_REPAIR"
        )
        if not ((math.isfinite(value) and value > 0.0) or gap):
            invalid.append(f"{key[0]:+.2f}:p{key[1]}")
    failures = [
        f"{key[0]:+.2f}:p{key[1]}"
        for key in sorted(expected | set(counts))
        if counts.get(key, 0) != 1
    ]
    return {
        "status": "PASS" if not failures and not invalid else "FAIL",
        "row_count": len(selected),
        "key_failures": failures,
        "invalid_values": invalid,
    }


def _sort_rows(rows: list[dict[str, Any]]) -> None:
    rows.sort(
        key=lambda row: (
            float(row["xi"]),
            0 if str(row["grid_kind"]) == "BASE" else 1,
            int(row["sorted_position"]),
            str(row.get("repair_id", "")),
        )
    )


def _write_point_transaction(
    output_dir: Path,
    existing_rows: Sequence[Mapping[str, Any]],
    solution: PointSolution,
) -> list[dict[str, Any]]:
    rows = [dict(row) for row in existing_rows]
    group = [
        row
        for row in rows
        if str(row["grid_kind"]) == "BASE"
        and round(float(row["xi"]), 10) == round(solution.xi, 10)
    ]
    if group:
        if _base_group_is_acceptable(group):
            return rows
        rows = [
            row
            for row in rows
            if not (
                round(float(row["xi"]), 10) == round(solution.xi, 10)
                and str(row["grid_kind"]) in {"BASE", "LOCAL_REFINEMENT"}
            )
        ]
    rows.extend(dict(row) for row in solution.rows)
    _sort_rows(rows)
    rlb2e._atomic_write_csv(
        Path(output_dir) / SPECTRUM_FILENAME, rows, SPECTRUM_FIELDS
    )
    return rows


def _solution_record(solution: PointSolution, *, benchmark: bool) -> dict[str, Any]:
    return {
        "xi": solution.xi,
        "benchmark": benchmark,
        "continuation_leg": solution.continuation_leg,
        "swapped_arms": solution.swapped_arms,
        "wall_time_seconds": solution.wall_time_seconds,
        "peak_rss_bytes": solution.peak_rss_bytes,
        "determinant_evaluations": solution.determinant_evaluations,
        "sigma_evaluations": solution.sigma_evaluations,
        "initial_Omega_max": solution.search_right_Omega,
        "actual_Omega_max": solution.search_right_Omega,
        "local_refinements": solution.local_refinements,
        "solve_mode": solution.solve_mode,
        "fallback_used": solution.fallback_used,
        "unresolved_candidates_below_root9": (
            solution.unresolved_candidates_below_root9
        ),
        "candidate_count_total": solution.candidate_count_total,
        "accepted_candidates_above_root9": (
            solution.accepted_candidates_above_root9
        ),
        "retained_slots_above_root9": solution.retained_slots_above_root9,
        "roots_above_9_computed": solution.roots_above_9_computed,
        "roots": [
            {
                "sorted_position": int(row["sorted_position"]),
                "Omega": float(row["Omega"]),
                "Lambda": float(row["Lambda"]),
                "singular_ratio": float(row["scaled_sigma_ratio"]),
                "boundary_residual": float(row["boundary_null_residual"]),
                "quality_status": str(row["quality_status"]),
            }
            for row in solution.rows
        ],
    }


def contract_payload() -> dict[str, Any]:
    return {
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "frequency_map_policy": FREQUENCY_MAP_POLICY,
        "geometry": {
            "mu": MU,
            "tau": TAU,
            "beta_deg": BETA_DEG,
            "l": L_REFERENCE,
            "L1": L1,
            "L2": L2,
            "L_total": L_TOTAL,
            "b1": WIDTH,
            "b2": WIDTH,
            "h1": THICKNESS,
            "h2": THICKNESS,
            "ply_thickness": PLY_THICKNESS,
            "K": K,
            "outer_clamps": True,
            "joint": "frozen ideal rigid RLB joint",
        },
        "material_M0": base_material_contract(),
        "signed_xi_definition": {
            "range": [XI_MIN, XI_MAX],
            "outer_stiffness_scale": "1+xi",
            "inner_stiffness_scale": "1-xi",
            "nu12_fixed": NU12_0,
            "rho_fixed": RHO_0,
            "xi_negative_meaning": "stiffer inner plies",
            "xi_positive_meaning": "stiffer outer plies",
            "all_ply_angles_deg": 0.0,
        },
        "stacks_bottom_to_top": {
            "arm1_LAYERED": list(LAYERED_LAYOUT),
            "arm2_HOMOGENEOUS_REFERENCE": list(HOMOGENEOUS_LAYOUT),
        },
        "xi_grid": [float(value) for value in xi_grid()],
        "continuation": {
            "anchor": 0.0,
            "negative_leg": [float(value) for value in continuation_paths()[0]],
            "positive_leg": [float(value) for value in continuation_paths()[1]],
            "xi0_calculated_once": True,
        },
        "normalization": {
            "Omega": "omega*l^2*sqrt(rho0*A0/(E0*Iy0))",
            "Lambda": "sqrt(Omega)",
            "E0": 1.0,
            "rho0": RHO_0,
            "b0": WIDTH,
            "h0": THICKNESS,
            "l": L_REFERENCE,
            "Omega_per_omega": OMEGA_TO_OMEGA_SCALE,
        },
        "thresholds": {
            "matrix_relative": MATRIX_RELATIVE_TOLERANCE,
            "symmetry_relative": SYMMETRY_RELATIVE_TOLERANCE,
            "reduced_property_relative": REDUCED_PROPERTY_TOLERANCE,
            "reduction_route_relative": REDUCTION_ROUTE_TOLERANCE,
            "root_singular_ratio": ROOT_SINGULAR_RATIO_TOLERANCE,
            "boundary_residual": BOUNDARY_RESIDUAL_TOLERANCE,
            "arm_swap_relative": ARM_SWAP_RELATIVE_TOLERANCE,
            "neighbour_MAD_multiplier": NEIGHBOUR_MAD_MULTIPLIER,
            "neighbour_absolute_trigger": NEIGHBOUR_ABSOLUTE_TRIGGER,
        },
        "root_search_reuse": {
            "source": RLB2E_SCRIPT_PATH,
            "requested_roots": K_PLOT,
            "guard_roots": 1,
            "required_slots": K_GUARD,
            "physics_formulas_copied": False,
        },
        "explicit_exclusions": [
            "roots_10_and_above",
            "spectral_sweep_runner",
            "certified_audit",
            "branch_tracking",
            "MAC",
            "mode_shapes",
            "energy_analysis",
            "Ritz",
            "FEM",
            "smoothing",
            "interpolation_based_frequencies",
            "commit",
            "push",
        ],
    }


def contract_hash() -> str:
    payload = json.dumps(
        rlb2e._json_value(contract_payload()),
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest().upper()


def _checkpoint_payload(
    rows: Sequence[Mapping[str, Any]],
    point_records: Sequence[Mapping[str, Any]],
    *,
    constitutive: Mapping[str, Any],
    started_at: str,
    benchmark_status: str,
) -> dict[str, Any]:
    groups = _complete_base_group_index(rows)
    missing = [
        float(value)
        for value in xi_grid()
        if round(float(value), 10) not in groups
    ]
    return {
        "schema_version": 1,
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "contract_sha256": contract_hash(),
        "started_at_utc": started_at,
        "updated_at_utc": rlb2e._utc_now(),
        "constitutive_gate": constitutive["status"],
        "benchmark_status": benchmark_status,
        "completed_points": len(groups),
        "completed_base_rows": len(groups) * K_GUARD,
        "missing_points": missing,
        "failed_points": [],
        "qualified_points": [],
        "last_completed_parameter": (
            point_records[-1].get("xi") if point_records else None
        ),
        "point_records": list(point_records),
        "root_calculation_count": ROOT_CALCULATION_COUNT,
        "parallel_workers_used": 0,
        "thread_limits": {
            name: os.environ.get(name, "")
            for name in (
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            )
        },
    }


def _write_failure_checkpoint(
    output_dir: Path,
    rows: Sequence[Mapping[str, Any]],
    point_records: Sequence[Mapping[str, Any]],
    *,
    constitutive: Mapping[str, Any],
    started_at: str,
    benchmark_status: str,
    xi: float,
    error: BaseException,
) -> None:
    payload = _checkpoint_payload(
        rows,
        point_records,
        constitutive=constitutive,
        started_at=started_at,
        benchmark_status=benchmark_status,
    )
    payload["failed_points"] = [
        {
            "xi": float(xi),
            "error": f"{type(error).__name__}: {error}",
            "recorded_at_utc": rlb2e._utc_now(),
        }
    ]
    rlb2e._atomic_write_json(Path(output_dir) / CHECKPOINT_FILENAME, payload)


def _benchmark_payload(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    endpoint_times = [
        float(item["wall_time_seconds"])
        for item in records
        if abs(float(item.get("xi", 0.0))) == 0.8
    ]
    remaining = len(xi_grid()) - 3
    maximum = max(endpoint_times, default=0.0)
    eta = 1.5 * maximum * remaining
    return {
        "schema_version": 1,
        "stage_id": STAGE_ID,
        "anchors": [dict(item) for item in records],
        "total_unique_base_points": len(xi_grid()),
        "declared_anchor_points": 3,
        "remaining_base_points": remaining,
        "endpoint_max_wall_time_seconds": maximum,
        "eta_formula": "1.5*max(time_xi_minus_0.8,time_xi_plus_0.8)*78",
        "conservative_eta_seconds": eta,
        "eta_limit_seconds": ETA_LIMIT_SECONDS,
        "production_run_permitted": eta <= ETA_LIMIT_SECONDS,
        "xi0_calculated_once": True,
    }


def _existing_rows(output_dir: Path) -> list[dict[str, Any]]:
    path = Path(output_dir) / SPECTRUM_FILENAME
    return [] if not path.is_file() else rlb2e._read_csv(path)


def _existing_point_records(output_dir: Path) -> list[dict[str, Any]]:
    path = Path(output_dir) / CHECKPOINT_FILENAME
    if not path.is_file():
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("contract_sha256") != contract_hash():
        raise RuntimeError("Existing RLB-2F checkpoint has an incompatible contract.")
    return [dict(item) for item in payload.get("point_records", [])]


def _recover_interrupted_finalization(
    output_dir: Path, rows: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    """Discard only non-durable repair rows after an interrupted finalization.

    BASE point transactions are already durable and are never recalculated.
    LOCAL_REFINEMENT provenance becomes durable only in the final manifest, so
    an incomplete finalization must repeat those triggered checks rather than
    retain two competing canonical sources.
    """

    if not any(str(row.get("grid_kind")) == "LOCAL_REFINEMENT" for row in rows):
        return [dict(row) for row in rows]
    recovered = [
        dict(row) for row in rows if str(row.get("grid_kind")) == "BASE"
    ]
    for row in recovered:
        row["is_canonical_plot_source"] = True
        if str(row.get("point_status")) != "PASS":
            row["point_status"] = "PASS"
    _sort_rows(recovered)
    rlb2e._atomic_write_csv(
        Path(output_dir) / SPECTRUM_FILENAME, recovered, SPECTRUM_FIELDS
    )
    return recovered


def _reconcile_nonbenchmark_point_records(
    rows: Sequence[Mapping[str, Any]], point_records: list[dict[str, Any]]
) -> None:
    """Record durable non-anchor CSV groups whose metrics sidecar was interrupted."""

    recorded = {
        round(float(item.get("xi", math.nan)), 10) for item in point_records
    }
    anchors = {round(value, 10) for value in (0.0, -0.8, 0.8)}
    for xi in sorted(_complete_base_group_index(rows)):
        if xi in recorded or xi in anchors:
            continue
        point_records.append(
            {
                "xi": xi,
                "benchmark": False,
                "reused_existing_without_metrics": True,
                "metrics_available": False,
                "wall_time_seconds": 0.0,
                "peak_rss_bytes": 0,
                "determinant_evaluations": 0,
                "sigma_evaluations": 0,
                "local_refinements": 0,
                "solve_mode": "RECOVERED_DURABLE_BASE_TRANSACTION",
                "fallback_used": False,
                "swapped_arms": False,
            }
        )


def run_benchmarks(
    output_dir: Path,
    rows: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    started_at: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Calculate or reuse xi=0,-0.8,+0.8 and enforce the measured ETA."""

    target = Path(output_dir)
    existing_index = {
        round(float(item.get("xi", math.nan)), 10): item
        for item in point_records
        if bool(item.get("benchmark", False))
    }
    benchmark_records: list[dict[str, Any]] = []
    for xi, leg in ((0.0, "ANCHOR"), (-0.8, "NEGATIVE"), (0.8, "POSITIVE")):
        groups = _complete_base_group_index(rows)
        key = round(xi, 10)
        if key in groups:
            if key not in existing_index:
                raise RuntimeError(
                    f"Benchmark xi={xi:+.2f} exists without measured metrics."
                )
            benchmark_records.append(dict(existing_index[key]))
            continue
        try:
            solution = solve_point(
                xi, force_anchor=True, continuation_leg=leg
            )
        except Exception as exc:
            _write_failure_checkpoint(
                target,
                rows,
                point_records,
                constitutive=constitutive,
                started_at=started_at,
                benchmark_status="RUNNING",
                xi=xi,
                error=exc,
            )
            raise
        rows = _write_point_transaction(target, rows, solution)
        record = _solution_record(solution, benchmark=True)
        point_records.append(record)
        benchmark_records.append(record)
        rlb2e._atomic_write_json(
            target / BENCHMARK_FILENAME, _benchmark_payload(benchmark_records)
        )
        rlb2e._atomic_write_json(
            target / CHECKPOINT_FILENAME,
            _checkpoint_payload(
                rows,
                point_records,
                constitutive=constitutive,
                started_at=started_at,
                benchmark_status="RUNNING",
            ),
        )
    benchmark = _benchmark_payload(benchmark_records)
    status = "PASS" if benchmark["production_run_permitted"] else "STOPPED_BY_ETA_GATE"
    rlb2e._atomic_write_json(target / BENCHMARK_FILENAME, benchmark)
    rlb2e._atomic_write_json(
        target / CHECKPOINT_FILENAME,
        _checkpoint_payload(
            rows,
            point_records,
            constitutive=constitutive,
            started_at=started_at,
            benchmark_status=status,
        ),
    )
    return rows, benchmark


def complete_missing_points(
    output_dir: Path,
    rows: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    started_at: str,
) -> list[dict[str, Any]]:
    """Traverse both declared paths from xi=0 and preserve completed groups."""

    target = Path(output_dir)
    for path, leg in zip(continuation_paths(), ("NEGATIVE", "POSITIVE"), strict=True):
        history: list[float] = []
        for xi_value in path:
            xi = round(float(xi_value), 10)
            groups = _complete_base_group_index(rows)
            if xi in groups:
                history.append(xi)
                continue
            previous = None
            second = None
            if history:
                previous = (history[-1], _rows_for_roots(rows, history[-1]))
            if len(history) >= 2:
                second = (history[-2], _rows_for_roots(rows, history[-2]))
            try:
                solution = solve_point(
                    xi,
                    previous=previous,
                    second_previous=second,
                    force_anchor=previous is None,
                    continuation_leg=leg,
                )
            except Exception as exc:
                _write_failure_checkpoint(
                    target,
                    rows,
                    point_records,
                    constitutive=constitutive,
                    started_at=started_at,
                    benchmark_status="PASS",
                    xi=xi,
                    error=exc,
                )
                raise
            rows = _write_point_transaction(target, rows, solution)
            point_records.append(_solution_record(solution, benchmark=False))
            history.append(xi)
            rlb2e._atomic_write_json(
                target / CHECKPOINT_FILENAME,
                _checkpoint_payload(
                    rows,
                    point_records,
                    constitutive=constitutive,
                    started_at=started_at,
                    benchmark_status="PASS",
                ),
            )
    return rows


def neighbour_audit_rows(
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Apply the unchanged RLB-2E median+8*MAD neighbour trigger."""

    grid = [round(float(value), 10) for value in xi_grid()]
    spectra = {
        xi: np.sqrt(_rows_for_roots(rows, xi)[:K_PLOT]) for xi in grid
    }
    gap_flags: set[tuple[float, int]] = set()
    for lower_position in range(1, K_PLOT):
        gaps = np.asarray(
            [
                spectra[xi][lower_position] - spectra[xi][lower_position - 1]
                for xi in grid
            ],
            dtype=float,
        )
        residuals = np.asarray(
            [
                rlb2e.centred_secant_residual(
                    gaps[index - 1], gaps[index], gaps[index + 1]
                )
                for index in range(1, len(grid) - 1)
            ],
            dtype=float,
        )
        median = float(np.median(residuals))
        mad = float(np.median(np.abs(residuals - median)))
        threshold = median + NEIGHBOUR_MAD_MULTIPLIER * mad
        for offset, index in enumerate(range(1, len(grid) - 1)):
            if (
                float(residuals[offset]) > threshold
                and float(residuals[offset]) > NEIGHBOUR_ABSOLUTE_TRIGGER
            ):
                gap_flags.add((grid[index], lower_position))
                gap_flags.add((grid[index], lower_position + 1))

    result: list[dict[str, Any]] = []
    for position in range(1, K_PLOT + 1):
        values = {xi: float(spectra[xi][position - 1]) for xi in grid}
        residuals = [
            rlb2e.centred_secant_residual(
                values[grid[index - 1]],
                values[grid[index]],
                values[grid[index + 1]],
            )
            for index in range(1, len(grid) - 1)
        ]
        median = float(np.median(residuals))
        mad = float(np.median(np.abs(np.asarray(residuals) - median)))
        robust_threshold = median + NEIGHBOUR_MAD_MULTIPLIER * mad
        for offset, index in enumerate(range(1, len(grid) - 1)):
            xi = grid[index]
            group = _base_group_index(rows)[xi]
            ordered = sorted(group, key=lambda row: int(row["sorted_position"]))
            root_row = ordered[position - 1]
            residual = float(residuals[offset])
            statistical_flag = bool(
                residual > robust_threshold
                and residual > NEIGHBOUR_ABSOLUTE_TRIGGER
            )
            root_count_warning = len(ordered) != K_GUARD
            ordering_warning = any(
                float(left["Omega"]) >= float(right["Omega"])
                for left, right in zip(ordered[:-1], ordered[1:], strict=True)
            )
            unresolved_warning = int(
                root_row["unresolved_candidates_below_root9"]
            ) != 0
            bad_residual_warning = bool(
                float(root_row["scaled_sigma_ratio"])
                > ROOT_SINGULAR_RATIO_TOLERANCE
                or float(root_row["boundary_null_residual"])
                > BOUNDARY_RESIDUAL_TOLERANCE
            )
            gap_warning = (xi, position) in gap_flags
            flagged = bool(
                statistical_flag
                or root_count_warning
                or ordering_warning
                or unresolved_warning
                or bad_residual_warning
                or gap_warning
            )
            result.append(
                {
                    "configuration_id": CONFIGURATION_ID,
                    "sorted_position": position,
                    "xi_left": grid[index - 1],
                    "xi": xi,
                    "xi_right": grid[index + 1],
                    "Lambda_left": values[grid[index - 1]],
                    "Lambda_center": values[xi],
                    "Lambda_right": values[grid[index + 1]],
                    "centered_predictor_Lambda": 0.5
                    * (values[grid[index - 1]] + values[grid[index + 1]]),
                    "centered_secant_residual": residual,
                    "median_residual_for_position": median,
                    "MAD_residual_for_position": mad,
                    "robust_threshold": robust_threshold,
                    "absolute_trigger": NEIGHBOUR_ABSOLUTE_TRIGGER,
                    "statistical_flag": statistical_flag,
                    "root_count_warning": root_count_warning,
                    "ordering_warning": ordering_warning,
                    "unresolved_candidate_warning": unresolved_warning,
                    "bad_residual_warning": bad_residual_warning,
                    "gap_jump_warning": gap_warning,
                    "flagged": flagged,
                    "repair_id": "",
                    "repair_status": "PENDING" if flagged else "NOT_REQUIRED",
                    "local_xi_values": [],
                    "smoothing_applied": False,
                }
            )
    return result


def flagged_repair_points(
    audit_rows: Sequence[Mapping[str, Any]],
) -> list[float]:
    return sorted(
        {
            round(float(row["xi"]), 10)
            for row in audit_rows
            if _as_bool(row["flagged"])
        }
    )


def apply_local_repairs(
    rows: list[dict[str, Any]], audit_rows: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Re-refine only flagged root neighbourhoods without smoothing."""

    records: list[dict[str, Any]] = []
    for repair_index, xi in enumerate(flagged_repair_points(audit_rows), start=1):
        repair_id = f"repair_{repair_index:04d}"
        positions = sorted(
            {
                int(row["sorted_position"])
                for row in audit_rows
                if round(float(row["xi"]), 10) == xi and _as_bool(row["flagged"])
            }
        )
        existing = [
            row
            for row in rows
            if str(row["grid_kind"]) == "LOCAL_REFINEMENT"
            and round(float(row["xi"]), 10) == xi
            and str(row.get("repair_id", "")) == repair_id
        ]
        if sorted(int(row["sorted_position"]) for row in existing) == list(
            range(1, K_GUARD + 1)
        ) and all(
            str(row["point_status"])
            in {
                "REPRODUCED_AFTER_LOCAL_REPAIR",
                "LOCATOR_CORRECTED_AFTER_LOCAL_REPAIR",
            }
            for row in existing
        ):
            original = np.asarray(
                [
                    float(row["Omega"])
                    for row in sorted(
                        _base_group_index(rows)[xi],
                        key=lambda item: int(item["sorted_position"]),
                    )
                ]
            )
            refined = np.asarray(
                [
                    float(row["Omega"])
                    for row in sorted(existing, key=lambda item: int(item["sorted_position"]))
                ]
            )
            relative = float(
                np.max(
                    np.abs(original - refined)
                    / np.maximum.reduce(
                        (
                            np.abs(original),
                            np.abs(refined),
                            np.full(K_GUARD, np.finfo(float).tiny),
                        )
                    )
                )
            )
            status = str(existing[0]["point_status"])
            for audit in audit_rows:
                if round(float(audit["xi"]), 10) == xi and _as_bool(audit["flagged"]):
                    audit["repair_id"] = repair_id
                    audit["repair_status"] = status
            records.append(
                {
                    "repair_id": repair_id,
                    "xi": xi,
                    "status": status,
                    "affected_positions": positions,
                    "maximum_relative_Omega_change": relative,
                    "wall_time_seconds": 0.0,
                    "determinant_evaluations": 0,
                    "sigma_evaluations": 0,
                    "reused_existing_repair": True,
                    "smoothing_applied": False,
                    "predictor_used_as_final": False,
                }
            )
            continue
        original = _rows_for_roots(rows, xi)
        try:
            solution = solve_point(
                xi,
                previous=(xi, original),
                dense_local=True,
                dense_positions=positions,
                grid_kind="LOCAL_REFINEMENT",
                repair_id=repair_id,
                continuation_leg="LOCAL_REPAIR",
            )
        except (RuntimeError, ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            gaps: list[dict[str, Any]] = []
            for row in rows:
                if (
                    str(row["grid_kind"]) == "BASE"
                    and round(float(row["xi"]), 10) == xi
                    and int(row["sorted_position"]) in positions
                ):
                    row["is_canonical_plot_source"] = False
                    row["point_status"] = "UNRESOLVED_AFTER_LOCAL_REPAIR"
                    gap = dict(row)
                    gap["row_id"] = (
                        f"{CONFIGURATION_ID}__{xi:+.6f}__LOCAL_REFINEMENT__"
                        f"p{int(row['sorted_position']):02d}__{repair_id}_gap"
                    )
                    gap["grid_kind"] = "LOCAL_REFINEMENT"
                    gap["Lambda"] = math.nan
                    gap["quality_status"] = "UNRESOLVED_AFTER_LOCAL_REPAIR"
                    gap["is_canonical_plot_source"] = True
                    gap["supersedes_row_id"] = row["row_id"]
                    gap["repair_id"] = repair_id
                    gaps.append(gap)
            rows.extend(gaps)
            for audit in audit_rows:
                if round(float(audit["xi"]), 10) == xi and _as_bool(audit["flagged"]):
                    audit["repair_id"] = repair_id
                    audit["repair_status"] = "UNRESOLVED_AFTER_LOCAL_REPAIR"
            records.append(
                {
                    "repair_id": repair_id,
                    "xi": xi,
                    "status": "UNRESOLVED",
                    "affected_positions": positions,
                    "error": f"{type(exc).__name__}: {exc}",
                    "wall_time_seconds": 0.0,
                    "determinant_evaluations": 0,
                    "sigma_evaluations": 0,
                    "smoothing_applied": False,
                    "predictor_used_as_final": False,
                }
            )
            continue
        refined = np.asarray([float(row["Omega"]) for row in solution.rows])
        relative = float(
            np.max(
                np.abs(original - refined)
                / np.maximum.reduce(
                    (
                        np.abs(original),
                        np.abs(refined),
                        np.full(K_GUARD, np.finfo(float).tiny),
                    )
                )
            )
        )
        status = (
            "REPRODUCED_AFTER_LOCAL_REPAIR"
            if relative <= 1.0e-8
            else "LOCATOR_CORRECTED_AFTER_LOCAL_REPAIR"
        )
        for row in rows:
            if str(row["grid_kind"]) == "BASE" and round(float(row["xi"]), 10) == xi:
                row["is_canonical_plot_source"] = False
                row["point_status"] = status
        for row in solution.rows:
            row["grid_kind"] = "LOCAL_REFINEMENT"
            row["repair_id"] = repair_id
            row["is_canonical_plot_source"] = True
            row["supersedes_row_id"] = next(
                str(item["row_id"])
                for item in rows
                if str(item["grid_kind"]) == "BASE"
                and round(float(item["xi"]), 10) == xi
                and int(item["sorted_position"]) == int(row["sorted_position"])
            )
            row["point_status"] = status
        rows.extend(dict(row) for row in solution.rows)
        for audit in audit_rows:
            if round(float(audit["xi"]), 10) == xi and _as_bool(audit["flagged"]):
                audit["repair_id"] = repair_id
                audit["repair_status"] = status
        records.append(
            {
                "repair_id": repair_id,
                "xi": xi,
                "status": status,
                "affected_positions": positions,
                "maximum_relative_Omega_change": relative,
                "wall_time_seconds": solution.wall_time_seconds,
                "peak_rss_bytes": solution.peak_rss_bytes,
                "determinant_evaluations": solution.determinant_evaluations,
                "sigma_evaluations": solution.sigma_evaluations,
                "smoothing_applied": False,
                "predictor_used_as_final": False,
            }
        )
    _sort_rows(rows)
    return rows, audit_rows, records


def create_plot_from_csv(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    """Render one panel from completed CSV without importing physics."""

    started = time.perf_counter()
    rows = rlb2e._read_csv(Path(output_dir) / SPECTRUM_FILENAME)
    spectrum_audit = audit_spectrum_rows(rows)
    if spectrum_audit["status"] != "PASS":
        raise RuntimeError(f"plot_only rejected root inventory: {spectrum_audit}")
    plot_audit = audit_plot_rows(rows)
    if plot_audit["status"] != "PASS":
        raise RuntimeError(f"plot_only rejected plot data: {plot_audit}")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(8.6, 5.2))
    colors = plt.get_cmap("tab10").colors[:K_PLOT]
    selected_rows = canonical_plot_rows(rows)
    for position in range(1, K_PLOT + 1):
        selected = [
            row
            for row in selected_rows
            if int(row["sorted_position"]) == position
        ]
        selected.sort(key=lambda row: float(row["xi"]))
        axis.plot(
            [float(row["xi"]) for row in selected],
            [float(row["Lambda"]) for row in selected],
            color=colors[position - 1],
            linestyle="-",
            linewidth=1.3,
            label=f"k={position}",
        )
    axis.axvline(0.0, color="0.45", linewidth=0.8, alpha=0.65)
    axis.set_xlabel(r"$\xi$")
    axis.set_ylabel(r"$\Lambda$")
    axis.grid(True, alpha=0.22, linewidth=0.5)
    axis.legend(ncol=4, frameon=False)
    figure.tight_layout()
    target = Path(output_dir) / PLOT_FILENAME
    temporary = target.with_name(target.stem + ".tmp" + target.suffix)
    figure.savefig(temporary, dpi=180, bbox_inches="tight")
    plt.close(figure)
    os.replace(temporary, target)
    return {
        "path": target.as_posix(),
        "wall_time_seconds": time.perf_counter() - started,
        "panel_count": 1,
        "spectral_line_count": K_PLOT,
        "reference_line_count": 1,
        "plotted_positions": list(range(1, K_PLOT + 1)),
        "root9_plotted": False,
        "root_calculation_count": 0,
        "spectrum_data_audit": spectrum_audit,
        "plot_data_audit": plot_audit,
    }


def arm_swap_checks(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    for xi in (-0.8, -0.4, 0.4, 0.8):
        reference = _rows_for_roots(rows, xi)
        solution = solve_point(
            xi,
            previous=(xi, reference),
            continuation_leg="ARM_SWAP_DIAGNOSTIC",
            swapped=True,
            guard_locator_right_width=0.2,
            grid_kind="ARM_SWAP_DIAGNOSTIC",
        )
        swapped = np.asarray([float(row["Omega"]) for row in solution.rows])
        relative = np.abs(reference - swapped) / np.maximum.reduce(
            (np.abs(reference), np.abs(swapped), np.full(K_GUARD, np.finfo(float).tiny))
        )
        maximum = float(np.max(relative))
        checks.append(
            {
                "xi": xi,
                "reference_arm_roles": ["LAYERED", "HOMOGENEOUS_REFERENCE"],
                "swapped_arm_roles": ["HOMOGENEOUS_REFERENCE", "LAYERED"],
                "reference_Omega": reference,
                "swapped_Omega": swapped,
                "relative_differences": relative,
                "maximum_relative_Omega": maximum,
                "tolerance": ARM_SWAP_RELATIVE_TOLERANCE,
                "status": "PASS" if maximum <= ARM_SWAP_RELATIVE_TOLERANCE else "FAIL",
                "root_count": len(swapped),
                "roots_above_9_computed": False,
                "wall_time_seconds": solution.wall_time_seconds,
                "peak_rss_bytes": solution.peak_rss_bytes,
                "determinant_evaluations": solution.determinant_evaluations,
                "sigma_evaluations": solution.sigma_evaluations,
                "root_diagnostics": [dict(row) for row in solution.rows],
            }
        )
    return {
        "stage_id": STAGE_ID,
        "status": "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL",
        "tolerance": ARM_SWAP_RELATIVE_TOLERANCE,
        "checks": checks,
        "plotted_configuration_count": 1,
        "diagnostic_configuration_exported_to_spectrum_csv": False,
    }


def baseline_consistency(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Check xi=0 matrix identity and optional frozen RLB-2E spectrum evidence."""

    provider, _metadata = make_matrix_provider(0.0)
    _beam, coupled = rlb2e._physics_modules()
    homogeneous = build_homogeneous_section()

    def homogeneous_provider(omega: float) -> FloatArray:
        return np.asarray(
            coupled.coupled_boundary_matrix(
                float(omega),
                math.radians(BETA_DEG),
                L1,
                homogeneous.properties,
                L2,
                homogeneous.properties,
            ),
            dtype=float,
        )

    matrix_residual = max(
        float(np.max(np.abs(provider(probe) - homogeneous_provider(probe))))
        for probe in (0.731, 3.217, 7.553)
    )
    external: dict[str, Any] = {
        "status": "NOT_AVAILABLE",
        "used_for_reproducibility": False,
    }
    manifest_path = RLB2E_RESULT_DIR / rlb2e.MANIFEST_FILENAME
    spectrum_path = RLB2E_RESULT_DIR / rlb2e.SPECTRUM_FILENAME
    if manifest_path.is_file() and spectrum_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        hashes = manifest.get("output_hashes", {})
        hash_failures = [
            name
            for name, expected in hashes.items()
            if not (RLB2E_RESULT_DIR / name).is_file()
            or rlb2e._sha256(RLB2E_RESULT_DIR / name) != expected
        ]
        if not hash_failures:
            predecessor_rows = rlb2e._read_csv(spectrum_path)
            predecessor = rlb2e._rows_for_roots(
                predecessor_rows, rlb2e.CONFIG_BOTH_OUTER, 0.0
            )
            current = _rows_for_roots(rows, 0.0)
            relative = np.abs(current - predecessor) / np.maximum.reduce(
                (
                    np.abs(current),
                    np.abs(predecessor),
                    np.full(K_GUARD, np.finfo(float).tiny),
                )
            )
            maximum = float(np.max(relative))
            external = {
                "status": "PASS" if maximum <= 1.0e-8 else "FAIL",
                "used_for_reproducibility": False,
                "evidence_role": "optional predecessor regression",
                "result_directory": RLB2E_RESULT_DIR.as_posix(),
                "manifest_sha256": rlb2e._sha256(manifest_path),
                "spectrum_sha256": rlb2e._sha256(spectrum_path),
                "output_hash_failures": [],
                "maximum_relative_Omega": maximum,
                "tolerance": 1.0e-8,
            }
        else:
            external = {
                "status": "UNTRUSTED_OPTIONAL_EVIDENCE",
                "used_for_reproducibility": False,
                "evidence_role": "optional predecessor regression",
                "output_hash_failures": hash_failures,
            }
    internal_pass = matrix_residual <= 16.0 * np.finfo(float).eps
    external_verified_failure = external.get("status") == "FAIL"
    overall = internal_pass and not external_verified_failure
    return {
        "status": "PASS" if overall else "FAIL",
        "xi0_layered_vs_homogeneous_matrix_max_abs": matrix_residual,
        "matrix_tolerance": 16.0 * np.finfo(float).eps,
        "current_spectrum_reproducible_without_predecessor_tree": True,
        "predecessor_evidence": external,
    }


def _minimum_adjacent_gaps(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for xi in xi_grid():
        roots = np.sqrt(_rows_for_roots(rows, float(xi))[:K_PLOT])
        gaps = np.diff(roots)
        index = int(np.argmin(gaps))
        records.append(
            {
                "xi": float(xi),
                "minimum_adjacent_Lambda_gap": float(gaps[index]),
                "between_sorted_positions": [index + 1, index + 2],
            }
        )
    records.sort(key=lambda item: float(item["minimum_adjacent_Lambda_gap"]))
    for item in records:
        item["interpretation"] = (
            "candidate interval only; no branch, crossing, veering, or shape claim"
        )
    return records


def _endpoint_changes(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    left = np.sqrt(_rows_for_roots(rows, XI_MIN)[:K_PLOT])
    right = np.sqrt(_rows_for_roots(rows, XI_MAX)[:K_PLOT])
    return [
        {
            "sorted_position": position,
            "Lambda_xi_minus_0.8": float(left[position - 1]),
            "Lambda_xi_plus_0.8": float(right[position - 1]),
            "relative_endpoint_change": float(
                (right[position - 1] - left[position - 1])
                / max(abs(left[position - 1]), np.finfo(float).tiny)
            ),
        }
        for position in range(1, K_PLOT + 1)
    ]


def _monotonicity(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    grid = [float(value) for value in xi_grid()]
    for position in range(1, K_PLOT + 1):
        values = np.asarray(
            [
                float(_canonical_group(rows, xi)[position - 1]["Lambda"])
                for xi in grid
            ],
            dtype=float,
        )
        differences = np.diff(values)
        tolerance = 1.0e-12 * max(float(np.max(np.abs(values))), 1.0)
        if np.all(differences >= -tolerance):
            label = "NONDECREASING"
        elif np.all(differences <= tolerance):
            label = "NONINCREASING"
        else:
            label = "NON_MONOTONE"
        result.append(
            {
                "sorted_position": position,
                "classification": label,
                "classification_tolerance_Lambda": tolerance,
                "minimum_step_change_Lambda": float(np.min(differences)),
                "maximum_step_change_Lambda": float(np.max(differences)),
                "branch_identity_claimed": False,
            }
        )
    return result


def _root_quality_summary(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    base = [row for row in rows if str(row["grid_kind"]) == "BASE"]
    guards = [row for row in base if int(row["sorted_position"]) == K_GUARD]
    return {
        "base_quality_failures": sum(
            str(row["quality_status"]) != "PASS" for row in base
        ),
        "maximum_base_scaled_sigma_ratio": max(
            float(row["scaled_sigma_ratio"]) for row in base
        ),
        "maximum_base_boundary_null_residual": max(
            float(row["boundary_null_residual"]) for row in base
        ),
        "maximum_root9_scaled_sigma_ratio": max(
            float(row["scaled_sigma_ratio"]) for row in guards
        ),
        "maximum_root9_boundary_null_residual": max(
            float(row["boundary_null_residual"]) for row in guards
        ),
        "minimum_root9_right_margin_Omega": min(
            float(row["root9_right_margin_Omega"]) for row in guards
        ),
        "maximum_unresolved_candidates_below_root9": max(
            int(row["unresolved_candidates_below_root9"]) for row in base
        ),
    }


def _runtime_summary(
    point_records: Sequence[Mapping[str, Any]],
    repair_records: Sequence[Mapping[str, Any]],
    arm_swap: Mapping[str, Any],
    plot: Mapping[str, Any],
) -> dict[str, Any]:
    swap_records = arm_swap["checks"]
    base_seconds = sum(float(item.get("wall_time_seconds", 0.0)) for item in point_records)
    repair_seconds = sum(float(item.get("wall_time_seconds", 0.0)) for item in repair_records)
    swap_seconds = sum(float(item.get("wall_time_seconds", 0.0)) for item in swap_records)
    plot_seconds = float(plot.get("wall_time_seconds", 0.0))
    determinant = (
        sum(int(item.get("determinant_evaluations", 0)) for item in point_records)
        + sum(int(item.get("determinant_evaluations", 0)) for item in repair_records)
        + sum(int(item.get("determinant_evaluations", 0)) for item in swap_records)
    )
    sigma = (
        sum(int(item.get("sigma_evaluations", 0)) for item in point_records)
        + sum(int(item.get("sigma_evaluations", 0)) for item in repair_records)
        + sum(int(item.get("sigma_evaluations", 0)) for item in swap_records)
    )
    peak = max(
        [rlb2e._peak_rss_bytes()]
        + [int(item.get("peak_rss_bytes", 0)) for item in point_records]
        + [int(item.get("peak_rss_bytes", 0)) for item in repair_records]
        + [int(item.get("peak_rss_bytes", 0)) for item in swap_records]
    )
    spectrum_seconds = base_seconds + repair_seconds + swap_seconds
    return {
        "base_point_wall_time_sum_seconds": base_seconds,
        "local_repair_wall_time_sum_seconds": repair_seconds,
        "arm_swap_wall_time_sum_seconds": swap_seconds,
        "spectrum_runtime_seconds": spectrum_seconds,
        "plot_only_seconds": plot_seconds,
        "total_measured_workflow_seconds": spectrum_seconds + plot_seconds,
        "peak_rss_bytes": peak,
        "determinant_evaluations": determinant,
        "sigma_evaluations": sigma,
        "base_root_solve_count": sum(
            int(item.get("determinant_evaluations", 0)) > 0
            for item in point_records
        ),
        "local_repair_solve_count": sum(
            int(item.get("determinant_evaluations", 0)) > 0
            for item in repair_records
        ),
        "arm_swap_solve_count": len(swap_records),
        "parallel_spectral_workers": 0,
    }


def _output_hashes(output_dir: Path) -> dict[str, str]:
    return {
        path.name: rlb2e._sha256(path)
        for path in sorted(Path(output_dir).iterdir())
        if path.is_file() and path.name != MANIFEST_FILENAME
    }


def _report_text(manifest: Mapping[str, Any]) -> str:
    gate = manifest["constitutive_gate"]
    counts = manifest["counts"]
    audit = manifest["neighbour_audit"]
    runtime = manifest["runtime"]
    quality = manifest["root_quality_summary"]
    endpoints = "\n".join(
        "- position {sorted_position}: $\\Delta_{{end}}={relative_endpoint_change:.6g}$.".format(
            **item
        )
        for item in manifest["endpoint_changes"]
    )
    monotonic = "\n".join(
        f"- position {item['sorted_position']}: `{item['classification']}`."
        for item in manifest["sorted_position_monotonicity"]
    )
    gap = manifest["minimum_adjacent_sorted_gaps"][0]
    return rf"""# RLB-2F: знаковый контраст слоёв в одном плече

## Цель и область применимости

Рассматривается изменение первых восьми независимо упорядоченных частот при
перераспределении жёсткости между наружной и внутренней парами слоёв первого
плеча. Второе плечо сохраняет однородный материал $M_0$. Идентичность
модальных ветвей, формы, MAC, локализация и распределение энергии не
определяются.

## Фиксированный контракт

Использованы $\mu=\tau=0$, $\beta=30^\circ$, $L_1=L_2=l=1$,
$b=0.20$, $h=0.05$ и $K=5/6$. Каждое плечо содержит четыре равных
слоя толщины 0.0125. Ориентация материала во всех слоях равна $0^\circ$.

Базовый материал имеет $E_1=1.1$, $E_2=0.9$, $\nu_{{12}}=0.3$,
$G_{{12}}=G_{{13}}=G_{{23}}=1/2.6$ и $\rho=1$. В слоистом плече
наружные модули умножены на $1+\xi$, внутренние -- на $1-\xi$.
Плотность и коэффициент Пуассона не менялись. При $\xi<0$ внутренние слои
жёстче наружных; при $\xi>0$ жёстче наружные слои.

## Constitutive gate

Статус: **{gate['status']}**. Подтверждены соотношения

$$\frac{{D_{{layered}}}}{{D_0}}=1+\frac{{3\xi}}{{4}},\qquad
\frac{{D_{{total}}}}{{D_0}}=2+\frac{{3\xi}}{{4}}.$$

При $\xi=-0.8$, 0 и 0.8 первое отношение равно соответственно 0.4, 1 и
1.6. Значения $A_{{beam}}$, $S_{{beam}}$, $m$ и $J$ совпадают с
однородным плечом. Максимальный residual полной формулы для $D$ равен
{gate['maximum_residuals']['D_matrix_formula_relative']:.3e}; $B$ и $I_1$
остаются нулевыми в пределах объявленных допусков.

## Frequency-map policy и численная сборка

Локальный instance использует `frequency-map-v1`, режим `fast_plot` и
семантику `sorted_positions`. Сетка содержит 81 значение
$\xi=-0.80,-0.78,\ldots,0.80$. Расчёт начат при $\xi=0$ и продолжен
отдельно к отрицательному и положительному концам.

Позиции 1--8 показаны на рисунке, root 9 служит только проверкой полноты.
Корни с позициями 10 и выше не вычислялись. Predictor использовался только
для локализации; каждое сохранённое значение получено из замороженной
characteristic matrix с существующими determinant/SVD detector и refiner.

Использована нормировка

$$\Omega=\omega l^2\sqrt{{\rho_0A_0/(E_0I_{{y0}})}},\qquad
\Lambda=\sqrt{{\Omega}}.$$

## Результаты и диагностика

Получено {counts['base_points']}/81 base-точек и
{counts['base_rows']}/729 base-строк. Полных root-9 guards:
{counts['root9_guards']}/81. Сохранено {counts['local_refinement_rows']}
строк `LOCAL_REFINEMENT`.

Production benchmark дал консервативную ETA
{manifest['benchmark']['conservative_eta_seconds']:.1f} s для 78 точек при
лимите 2700 s. Полное измеренное время составило
{runtime['total_measured_workflow_seconds']:.1f} s, peak RSS --
{runtime['peak_rss_bytes'] / 2**20:.1f} MiB. Выполнено
{runtime['determinant_evaluations']} determinant и
{runtime['sigma_evaluations']} SVD/sigma evaluations.

Neighbour audit отметил {audit['flagged_point_count']} точек. Выполнено
{audit['repair_count']} локальных проверок; unresolved points:
{audit['unresolved_point_count']}. Сглаживание не применялось. Максимальные
base residuals равны {quality['maximum_base_scaled_sigma_ratio']:.3e} для
$\sigma_{{min}}/\sigma_{{max}}$ и
{quality['maximum_base_boundary_null_residual']:.3e} для boundary residual.

Проверка перестановки плеч имеет статус **{manifest['arm_swap']['status']}**.
Сопоставление при $\xi=0$ с однородной конструкцией имеет статус
**{manifest['baseline_consistency']['status']}**.

Относительные изменения между концами сетки:

{endpoints}

Классификация упорядоченных позиций на принятой сетке:

{monotonic}

Минимальный соседний интервал равен
$\Delta\Lambda={float(gap['minimum_adjacent_Lambda_gap']):.6g}$ при
$\xi={float(gap['xi']):.2f}$ между позициями
{int(gap['between_sorted_positions'][0])} и
{int(gap['between_sorted_positions'][1])}. Это только кандидат для
последующего анализа форм. По частотной карте нельзя установить crossing,
veering, обмен модальным характером или локализацию.

## Статус и ограничения

**RLB-2F: {manifest['scientific_status']}.** Результат относится только к
фиксированной геометрии, одной конфигурации и конечной сетке $\xi$. Все слои
имеют ориентацию $0^\circ$, а плотность постоянна. Branch tracking, MAC,
формы, energy analysis, Ritz, FEM, damping, complex roots и certified audit
не выполнялись. Production physics не изменялась.
"""


def finalize_outputs(
    output_dir: Path,
    rows: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    benchmark: Mapping[str, Any],
    point_records: Sequence[Mapping[str, Any]],
    started_at: str,
) -> dict[str, Any]:
    target = Path(output_dir)
    audit_rows = neighbour_audit_rows(rows)
    rows, audit_rows, repair_records = apply_local_repairs(rows, audit_rows)
    rlb2e._atomic_write_csv(target / SPECTRUM_FILENAME, rows, SPECTRUM_FIELDS)
    rlb2e._atomic_write_csv(target / AUDIT_FILENAME, audit_rows)
    spectrum_audit = audit_spectrum_rows(rows)
    if spectrum_audit["status"] != "PASS":
        raise RuntimeError(f"Final spectrum audit failed: {spectrum_audit}")
    arm_swap = arm_swap_checks(rows)
    rlb2e._atomic_write_json(target / ARM_SWAP_FILENAME, arm_swap)
    if arm_swap["status"] != "PASS":
        raise RuntimeError(f"Arm-swap diagnostic failed: {arm_swap}")
    baseline = baseline_consistency(rows)
    if baseline["status"] != "PASS":
        raise RuntimeError(f"Homogeneous baseline consistency failed: {baseline}")
    plot = create_plot_from_csv(target)
    unresolved = sum(record["status"] == "UNRESOLVED" for record in repair_records)
    scientific_status = "PASS" if unresolved == 0 else "PARTIAL_PASS"
    counts = {
        "base_points": spectrum_audit["base_group_count"],
        "base_rows": spectrum_audit["base_row_count"],
        "local_refinement_rows": sum(
            str(row["grid_kind"]) == "LOCAL_REFINEMENT" for row in rows
        ),
        "root9_guards": sum(
            str(row["grid_kind"]) == "BASE"
            and int(row["sorted_position"]) == K_GUARD
            for row in rows
        ),
        "computed_base_points": sum(
            int(record.get("determinant_evaluations", 0)) > 0
            and not bool(record.get("swapped_arms", False))
            for record in point_records
        ),
        "reused_base_points": len(xi_grid())
        - sum(
            int(record.get("determinant_evaluations", 0)) > 0
            and not bool(record.get("swapped_arms", False))
            for record in point_records
        ),
        "locally_repaired_points": len(repair_records),
        "unresolved_points": unresolved,
    }
    runtime = _runtime_summary(point_records, repair_records, arm_swap, plot)
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "scientific_status": scientific_status,
        "completed_at_utc": rlb2e._utc_now(),
        "git": rlb2e._git_state(),
        "contract": contract_payload(),
        "contract_sha256": contract_hash(),
        "production_physics_hashes": {
            path: rlb2e._sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS
        },
        "analysis_script_sha256": rlb2e._sha256(Path(__file__)),
        "reused_RLB2E_script_sha256": rlb2e._sha256(ROOT / RLB2E_SCRIPT_PATH),
        "constitutive_gate": dict(constitutive),
        "counts": counts,
        "spectrum_audit": spectrum_audit,
        "root_quality_summary": _root_quality_summary(rows),
        "benchmark": dict(benchmark),
        "point_records": list(point_records),
        "neighbour_audit": {
            "criterion": (
                "centered secant residual > median+8*MAD and >1e-3; "
                "plus root-count, ordering, unresolved-candidate, residual, "
                "and adjacent-gap triggers"
            ),
            "row_count": len(audit_rows),
            "flagged_point_count": len(flagged_repair_points(audit_rows)),
            "repair_count": len(repair_records),
            "repair_records": repair_records,
            "unresolved_point_count": unresolved,
            "smoothing_applied": False,
        },
        "arm_swap": arm_swap,
        "baseline_consistency": baseline,
        "endpoint_changes": _endpoint_changes(rows),
        "sorted_position_monotonicity": _monotonicity(rows),
        "minimum_adjacent_sorted_gaps": _minimum_adjacent_gaps(rows),
        "plot": plot,
        "runtime": runtime,
        "root_contract": {
            "plotted_positions": list(range(1, K_PLOT + 1)),
            "guard_position": K_GUARD,
            "root9_role": "completeness_only",
            "search_policy_requested_roots": K_PLOT,
            "search_policy_guard_roots": 1,
            "root9_plotted": False,
            "roots_above_9_computed": False,
            "accepted_candidates_above_root9": 0,
            "branch_tracking": False,
        },
        "exclusions_confirmed": {
            "spectral_sweep_runner_used": False,
            "certified_audit_run": False,
            "full_inventory_run": False,
            "parallel_spectral_workers": 0,
            "interpolation_based_frequencies": False,
            "smoothing": False,
            "branch_tracking": False,
            "MAC": False,
            "mode_shapes": False,
            "energy_analysis": False,
            "Ritz": False,
            "FEM": False,
            "commit": False,
            "push": False,
        },
    }
    checkpoint = _checkpoint_payload(
        rows,
        point_records,
        constitutive=constitutive,
        started_at=started_at,
        benchmark_status="PASS",
    )
    checkpoint["scientific_status"] = scientific_status
    checkpoint["completed_at_utc"] = rlb2e._utc_now()
    checkpoint["local_repair_count"] = len(repair_records)
    checkpoint["arm_swap_check_count"] = len(arm_swap["checks"])
    checkpoint["terminal_unresolved_points"] = [
        {"xi": record["xi"], "repair_id": record["repair_id"]}
        for record in repair_records
        if record["status"] == "UNRESOLVED"
    ]
    rlb2e._atomic_write_json(target / CHECKPOINT_FILENAME, checkpoint)
    rlb2e._atomic_write_text(target / REPORT_FILENAME, _report_text(manifest))
    manifest["output_hashes"] = _output_hashes(target)
    rlb2e._atomic_write_json(target / MANIFEST_FILENAME, manifest)
    return manifest


def _completed_manifest(output_dir: Path) -> dict[str, Any] | None:
    target = Path(output_dir)
    manifest_path = target / MANIFEST_FILENAME
    spectrum_path = target / SPECTRUM_FILENAME
    if not manifest_path.is_file() or not spectrum_path.is_file():
        return None
    rows = rlb2e._read_csv(spectrum_path)
    if audit_spectrum_rows(rows)["status"] != "PASS":
        return None
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    recorded_hashes = manifest.get("output_hashes", {})
    expected_physics = {
        path: rlb2e._sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS
    }
    compatible = bool(
        manifest.get("stage_id") == STAGE_ID
        and manifest.get("algorithm_version") == ALGORITHM_VERSION
        and manifest.get("contract_sha256") == contract_hash()
        and manifest.get("analysis_script_sha256") == rlb2e._sha256(Path(__file__))
        and manifest.get("production_physics_hashes") == expected_physics
        and manifest.get("reused_RLB2E_script_sha256")
        == rlb2e._sha256(ROOT / RLB2E_SCRIPT_PATH)
        and isinstance(recorded_hashes, dict)
        and MANDATORY_COMPLETED_OUTPUTS.issubset(recorded_hashes)
        and all(
            (target / name).is_file()
            and rlb2e._sha256(target / name) == expected
            for name, expected in recorded_hashes.items()
        )
    )
    return manifest if compatible else None


def run_workflow(
    output_dir: Path = DEFAULT_OUTPUT_DIR, *, missing_only: bool = False
) -> dict[str, Any]:
    """Run the finite RLB-2F map or resume only incomplete transactions."""

    global ROOT_CALCULATION_COUNT
    ROOT_CALCULATION_COUNT = 0
    target = Path(output_dir)
    completed = _completed_manifest(target)
    if completed is not None:
        result = dict(completed)
        result["invocation"] = {
            "missing_only": bool(missing_only),
            "root_calculation_count": 0,
            "outputs_modified": False,
        }
        return result

    target.mkdir(parents=True, exist_ok=True)
    started_at = rlb2e._utc_now()
    constitutive = constitutive_gate()
    section_rows = section_property_rows()
    rlb2e._atomic_write_csv(target / SECTION_FILENAME, section_rows)
    if constitutive["status"] != "PASS":
        raise RuntimeError(f"RLB-2F constitutive gate failed: {constitutive}")
    failed_sections = [
        f"xi={float(row['xi']):+.2f},arm={int(row['arm_id'])}"
        for row in section_rows
        if str(row["constitutive_status"]) != "PASS"
    ]
    if failed_sections:
        raise RuntimeError(
            "RLB-2F section-row constitutive gate failed: "
            + ", ".join(failed_sections)
        )
    rows = _existing_rows(target)
    rows = _recover_interrupted_finalization(target, rows)
    point_records = _existing_point_records(target)
    _reconcile_nonbenchmark_point_records(rows, point_records)
    rows, benchmark = run_benchmarks(
        target, rows, point_records, constitutive, started_at
    )
    if not benchmark["production_run_permitted"]:
        raise RuntimeError(
            "RLB-2F ETA exceeds 45 minutes: "
            f"{benchmark['conservative_eta_seconds']:.1f} s."
        )
    rows = complete_missing_points(
        target, rows, point_records, constitutive, started_at
    )
    return finalize_outputs(
        target, rows, constitutive, benchmark, point_records, started_at
    )


def manifest_only(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    """Validate and return the contract without calculating roots or writing."""

    return {
        "stage_id": STAGE_ID,
        "contract": contract_payload(),
        "contract_sha256": contract_hash(),
        "root_calculation_count": 0,
        "outputs_modified": False,
        "output_dir": Path(output_dir).as_posix(),
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--missing-only", action="store_true")
    mode.add_argument("--plot-only", action="store_true")
    mode.add_argument("--manifest-only", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.manifest_only:
        result = manifest_only(args.output_dir)
    elif args.plot_only:
        result = create_plot_from_csv(args.output_dir)
    else:
        result = run_workflow(args.output_dir, missing_only=args.missing_only)
    print(json.dumps(rlb2e._json_value(result), indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
