"""RLB-2G density-layout frequency maps for two coupled Reddy beams.

The workflow is the mass-layout analogue of RLB-2E and RLB-2F.  Every ply
retains the same elastic material; only its density changes.  The palindromic
four-ply layouts preserve the static stiffnesses and translational mass while
changing only the rotary inertia ``J`` of the reduced beam.

Two local ``frequency-map-v1`` instances are produced.  Experiment A compares
three two-arm density layouts on ``eta=0...0.8``.  Experiment B changes the
density layout of arm 1 on the signed ``xi_rho=-0.8...0.8`` grid while arm 2
remains homogeneous.  Independently sorted positions 1--8 are plotted and
position 9 is retained only as a completeness guard.

The characteristic matrix and root-search machinery are imported from the
verified RLB-2E/RLB-2F analysis stack.  Predictors are numerical locators only;
every exported frequency is refined from the frozen characteristic matrix.
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
    sweep_reddy_one_arm_layered_contrast as rlb2f,
)
from scripts.analysis.laminated_beams import (  # noqa: E402
    sweep_reddy_stiffness_layout_contrast as rlb2e,
)


FloatArray = NDArray[np.float64]

STAGE_ID = "RLB-2G"
ALGORITHM_VERSION = "mass_layout_duality_fast_plot_v1"
POLICY_ID = "frequency-map-v1"

EXPERIMENT_A = "A_THREE_CONFIGURATION_DENSITY_LAYOUT"
EXPERIMENT_B = "B_ONE_ARM_SIGNED_DENSITY_LAYOUT"

CONFIG_BOTH_OUTER = "BOTH_OUTER_HEAVY"
CONFIG_BOTH_INNER = "BOTH_INNER_HEAVY"
CONFIG_ANTI_PHASE = "ANTI_PHASE"
CONFIG_ONE_ARM = "ONE_ARM_DENSITY_LAYOUT"
CONFIGURATIONS_A = (
    CONFIG_BOTH_OUTER,
    CONFIG_BOTH_INNER,
    CONFIG_ANTI_PHASE,
)

HEAVY = "H"
LIGHT = "L"
HOMOGENEOUS = "M0"
OUTER_HEAVY_LAYOUT = (HEAVY, LIGHT, LIGHT, HEAVY)
INNER_HEAVY_LAYOUT = (LIGHT, HEAVY, HEAVY, LIGHT)
HOMOGENEOUS_LAYOUT = (HOMOGENEOUS,) * 4
CONFIGURATION_LAYOUTS_A = {
    CONFIG_BOTH_OUTER: (OUTER_HEAVY_LAYOUT, OUTER_HEAVY_LAYOUT),
    CONFIG_BOTH_INNER: (INNER_HEAVY_LAYOUT, INNER_HEAVY_LAYOUT),
    CONFIG_ANTI_PHASE: (INNER_HEAVY_LAYOUT, OUTER_HEAVY_LAYOUT),
}

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

ETA_MIN = 0.0
ETA_MAX = 0.8
ETA_STEP = 0.02
XI_RHO_MIN = -0.8
XI_RHO_MAX = 0.8
XI_RHO_STEP = 0.02
K_PLOT = rlb2e.K_PLOT
K_GUARD = rlb2e.K_GUARD

A_REFERENCE = WIDTH * THICKNESS
IY_REFERENCE = WIDTH * THICKNESS**3 / 12.0
OMEGA_TO_OMEGA_SCALE = rlb2e.OMEGA_TO_OMEGA_SCALE
M0 = RHO_0 * WIDTH * THICKNESS
J0 = RHO_0 * WIDTH * THICKNESS**3 / 12.0

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
    ROOT / "results" / "laminated_beams" / "reddy_mass_layout_duality"
)
SPECTRUM_A_FILENAME = "three_configuration_spectrum_roots.csv"
SPECTRUM_B_FILENAME = "one_arm_spectrum_roots.csv"
SECTION_FILENAME = "section_properties.csv"
AUDIT_FILENAME = "neighbour_audit.csv"
ARM_SWAP_FILENAME = "arm_swap_diagnostics.json"
BENCHMARK_FILENAME = "benchmark.json"
CHECKPOINT_FILENAME = "checkpoint.json"
PLOT_A_FILENAME = "lambda_vs_eta_three_density_configurations.png"
PLOT_B_FILENAME = "lambda_vs_xi_rho_one_arm_density.png"
PLOT_J_FILENAME = "J_over_J0_vs_contrast.png"
REPORT_FILENAME = "report.md"
MANIFEST_FILENAME = "run_manifest.json"

PRODUCTION_PHYSICS_PATHS = rlb2e.PRODUCTION_PHYSICS_PATHS
RLB2E_SCRIPT_PATH = (
    "scripts/analysis/laminated_beams/sweep_reddy_stiffness_layout_contrast.py"
)
RLB2F_SCRIPT_PATH = (
    "scripts/analysis/laminated_beams/sweep_reddy_one_arm_layered_contrast.py"
)
RLB2E_RESULT_DIR = (
    ROOT / "results" / "laminated_beams" / "reddy_stiffness_layout_contrast_sweep"
)
RLB2F_RESULT_DIR = (
    ROOT / "results" / "laminated_beams" / "reddy_one_arm_layered_contrast_sweep"
)

POLICY_EXPERIMENT_A = {
    "frequency_map_policy": POLICY_ID,
    "calculation_mode": "fast_plot",
    "spectrum_semantics": "sorted_positions",
    "sweep_parameter": "eta",
    "parameter_grid": "0.00:0.02:0.80",
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

POLICY_EXPERIMENT_B = {
    "frequency_map_policy": POLICY_ID,
    "calculation_mode": "fast_plot",
    "spectrum_semantics": "sorted_positions",
    "sweep_parameter": "xi_rho",
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
    "experiment_id",
    "configuration_id",
    "parameter_name",
    "parameter_value",
    "eta",
    "xi_rho",
    "grid_kind",
    "continuation_leg",
    "solve_id",
    "physical_solve_id",
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
    "shared_zero_contrast_anchor_reused",
    "shared_zero_contrast_source",
    "is_canonical_plot_source",
    "supersedes_row_id",
    "repair_id",
    "roots_above_9_computed",
)

SECTION_FIELDS = (
    "experiment_id",
    "parameter_name",
    "parameter_value",
    "configuration_id",
    "arm_id",
    "arm_role",
    "stack_bottom_to_top",
    "ply_densities",
    "angle_deg",
    "ply_thickness",
    "z_interfaces",
    "A_matrix",
    "B_matrix",
    "D_matrix",
    "shear_matrix_yz_xz",
    "I0",
    "I1",
    "I2",
    "A_beam",
    "D_beam",
    "S_beam",
    "m",
    "J",
    "J_over_J0",
    "expected_J_over_J0",
    "J_formula_residual",
    "B_relative",
    "I1_relative",
    "A_matrix_invariance_residual",
    "D_matrix_invariance_residual",
    "shear_matrix_invariance_residual",
    "A_beam_invariance_residual",
    "D_beam_invariance_residual",
    "S_beam_invariance_residual",
    "m_invariance_residual",
    "reduction_route_max_relative",
    "elastic_properties_identical_across_plies",
    "A_invariant",
    "D_invariant",
    "S_invariant",
    "m_invariant",
    "only_J_changes",
    "constitutive_status",
)

MANDATORY_COMPLETED_OUTPUTS = frozenset(
    {
        SPECTRUM_A_FILENAME,
        SPECTRUM_B_FILENAME,
        SECTION_FILENAME,
        AUDIT_FILENAME,
        ARM_SWAP_FILENAME,
        BENCHMARK_FILENAME,
        CHECKPOINT_FILENAME,
        PLOT_A_FILENAME,
        PLOT_B_FILENAME,
        PLOT_J_FILENAME,
        REPORT_FILENAME,
    }
)


@dataclass(frozen=True)
class SectionObjects:
    layout: tuple[str, str, str, str]
    ply_densities: tuple[float, float, float, float]
    laminate: Any
    properties: Any


@dataclass(frozen=True)
class PointSolution:
    experiment_id: str
    configuration_id: str
    parameter_name: str
    parameter_value: float
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


def eta_grid() -> FloatArray:
    return np.asarray(
        sweep_grid_policy.parameter_grid(ETA_MIN, ETA_MAX, ETA_STEP), dtype=float
    )


def xi_rho_grid() -> FloatArray:
    return np.asarray(
        sweep_grid_policy.parameter_grid(XI_RHO_MIN, XI_RHO_MAX, XI_RHO_STEP),
        dtype=float,
    )


def xi_rho_continuation_paths() -> tuple[FloatArray, FloatArray]:
    negative = -np.asarray(
        sweep_grid_policy.parameter_grid(0.0, abs(XI_RHO_MIN), XI_RHO_STEP),
        dtype=float,
    )
    negative[0] = 0.0
    positive = np.asarray(
        sweep_grid_policy.parameter_grid(0.0, XI_RHO_MAX, XI_RHO_STEP),
        dtype=float,
    )
    return negative, positive


def omega_to_Omega(omega: float) -> float:
    return rlb2e.omega_to_Omega(omega)


def Omega_to_Lambda(Omega: float) -> float:
    return rlb2e.Omega_to_Lambda(Omega)


def base_material_contract() -> dict[str, float]:
    return {
        "delta": DELTA,
        "E1": E1_0,
        "E2": E2_0,
        "nu12": NU12_0,
        "G12": G0,
        "G13": G0,
        "G23": G0,
        "rho0": RHO_0,
    }


def _density_material(rho: float, *, label: str) -> Any:
    beam, _coupled = rlb2e._physics_modules()
    return beam.OrthotropicLamina(
        E1=E1_0,
        E2=E2_0,
        nu12=NU12_0,
        G12=G0,
        G13=G0,
        G23=G0,
        rho=float(rho),
        name=f"RLB-2G {label}",
    )


def density_materials(eta: float) -> dict[str, Any]:
    value = float(eta)
    if not math.isfinite(value) or not ETA_MIN <= value <= ETA_MAX:
        raise ValueError("eta must lie in [0, 0.8].")
    return {
        HEAVY: _density_material((1.0 + value) * RHO_0, label=f"H eta={value:.2f}"),
        LIGHT: _density_material((1.0 - value) * RHO_0, label=f"L eta={value:.2f}"),
        HOMOGENEOUS: _density_material(RHO_0, label="M0"),
    }


def signed_density_values(xi_rho: float) -> tuple[float, float]:
    value = float(xi_rho)
    if not math.isfinite(value) or not XI_RHO_MIN <= value <= XI_RHO_MAX:
        raise ValueError("xi_rho must lie in [-0.8, 0.8].")
    return (1.0 + value) * RHO_0, (1.0 - value) * RHO_0


def _build_section(
    layout: Sequence[str], densities: Sequence[float]
) -> SectionObjects:
    labels = tuple(str(item) for item in layout)
    rho_values = tuple(float(value) for value in densities)
    if len(labels) != 4 or len(rho_values) != 4:
        raise ValueError("RLB-2G sections require exactly four plies.")
    if any(value <= 0.0 or not math.isfinite(value) for value in rho_values):
        raise ValueError("Every RLB-2G ply density must be finite and positive.")
    beam, _coupled = rlb2e._physics_modules()
    plies = [
        beam.Ply(
            _density_material(rho, label=label),
            angle_deg=0.0,
            thickness=PLY_THICKNESS,
            label=label,
        )
        for label, rho in zip(labels, rho_values, strict=True)
    ]
    laminate = beam.integrate_laminate(plies)
    properties = beam.reduce_to_beam_properties(
        laminate,
        width=WIDTH,
        K=K,
        symmetry_tolerance=SYMMETRY_RELATIVE_TOLERANCE,
        reduction_tolerance=REDUCTION_ROUTE_TOLERANCE,
    )
    return SectionObjects(labels, rho_values, laminate, properties)


def build_homogeneous_section() -> SectionObjects:
    return _build_section(HOMOGENEOUS_LAYOUT, (RHO_0,) * 4)


def build_eta_layout_section(layout: Sequence[str], eta: float) -> SectionObjects:
    labels = tuple(str(item) for item in layout)
    if labels not in (OUTER_HEAVY_LAYOUT, INNER_HEAVY_LAYOUT):
        raise ValueError(f"Unknown density layout: {labels!r}.")
    materials = density_materials(eta)
    return _build_section(labels, tuple(materials[label].rho for label in labels))


def build_experiment_A_sections(
    configuration_id: str, eta: float, *, swapped: bool = False
) -> tuple[SectionObjects, SectionObjects]:
    try:
        layouts = CONFIGURATION_LAYOUTS_A[str(configuration_id)]
    except KeyError as exc:
        raise ValueError(f"Unknown Experiment A configuration: {configuration_id!r}.") from exc
    sections = tuple(build_eta_layout_section(layout, eta) for layout in layouts)
    return (sections[1], sections[0]) if swapped else sections  # type: ignore[return-value]


def build_one_arm_layered_section(xi_rho: float) -> SectionObjects:
    rho_outer, rho_inner = signed_density_values(xi_rho)
    return _build_section(
        ("OUTER", "INNER", "INNER", "OUTER"),
        (rho_outer, rho_inner, rho_inner, rho_outer),
    )


def build_experiment_B_sections(
    xi_rho: float, *, swapped: bool = False
) -> tuple[SectionObjects, SectionObjects]:
    layered = build_one_arm_layered_section(xi_rho)
    homogeneous = build_homogeneous_section()
    return (homogeneous, layered) if swapped else (layered, homogeneous)


def _elastic_tuple(material: Any) -> tuple[float, ...]:
    return tuple(
        float(getattr(material, field))
        for field in ("E1", "E2", "nu12", "G12", "G13", "G23")
    )


def elastic_properties_identical(section: SectionObjects) -> bool:
    values = [_elastic_tuple(ply.material) for ply in section.laminate.plies]
    return all(value == values[0] for value in values[1:])


def _properties_only_J_can_change(left: Any, right: Any) -> bool:
    scalar_names = ("A", "D", "S", "m", "K", "width")
    scalars_equal = all(
        rlb2e._relative(getattr(left, name), getattr(right, name))
        <= REDUCED_PROPERTY_TOLERANCE
        for name in scalar_names
    )
    reductions_equal = all(
        rlb2e._relative(
            getattr(left, reduction_name).value,
            getattr(right, reduction_name).value,
        )
        <= REDUCTION_ROUTE_TOLERANCE
        for reduction_name in (
            "axial_reduction",
            "bending_reduction",
            "shear_reduction_before_K",
        )
    )
    return bool(scalars_equal and reductions_equal)


def _section_check(
    section: SectionObjects,
    baseline: SectionObjects,
    expected_J_ratio: float,
) -> dict[str, Any]:
    laminate = section.laminate
    properties = section.properties
    reference_laminate = baseline.laminate
    reference = baseline.properties
    A_matrix_residual = rlb2e._matrix_relative(laminate.A, reference_laminate.A)
    D_matrix_residual = rlb2e._matrix_relative(laminate.D, reference_laminate.D)
    shear_matrix_residual = rlb2e._matrix_relative(
        laminate.shear, reference_laminate.shear
    )
    A_residual = rlb2e._relative(properties.A, reference.A)
    D_residual = rlb2e._relative(properties.D, reference.D)
    S_residual = rlb2e._relative(properties.S, reference.S)
    m_residual = rlb2e._relative(properties.m, reference.m)
    J_ratio = float(properties.J / reference.J)
    J_formula_residual = rlb2e._relative(J_ratio, expected_J_ratio)
    B_relative = rlb2e._scaled_B(laminate)
    I1_relative = rlb2e._scaled_I1(laminate)
    I0_residual = rlb2e._relative(laminate.I0, reference_laminate.I0)
    route_residual = rlb2e._reduction_max_residual(properties)
    elastic_identical = elastic_properties_identical(section)
    only_J_changes = _properties_only_J_can_change(properties, reference)
    passed = bool(
        elastic_identical
        and A_matrix_residual <= MATRIX_RELATIVE_TOLERANCE
        and D_matrix_residual <= MATRIX_RELATIVE_TOLERANCE
        and shear_matrix_residual <= MATRIX_RELATIVE_TOLERANCE
        and max(A_residual, D_residual, S_residual, m_residual)
        <= REDUCED_PROPERTY_TOLERANCE
        and I0_residual <= MATRIX_RELATIVE_TOLERANCE
        and B_relative <= SYMMETRY_RELATIVE_TOLERANCE
        and I1_relative <= SYMMETRY_RELATIVE_TOLERANCE
        and J_formula_residual <= REDUCED_PROPERTY_TOLERANCE
        and route_residual <= REDUCTION_ROUTE_TOLERANCE
        and only_J_changes
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "passed": passed,
        "elastic_properties_identical_across_plies": elastic_identical,
        "A_matrix_invariance_residual": A_matrix_residual,
        "D_matrix_invariance_residual": D_matrix_residual,
        "shear_matrix_invariance_residual": shear_matrix_residual,
        "A_beam_invariance_residual": A_residual,
        "D_beam_invariance_residual": D_residual,
        "S_beam_invariance_residual": S_residual,
        "m_invariance_residual": m_residual,
        "I0_invariance_residual": I0_residual,
        "B_relative": B_relative,
        "I1_relative": I1_relative,
        "J_over_J0": J_ratio,
        "expected_J_over_J0": expected_J_ratio,
        "J_formula_residual": J_formula_residual,
        "reduction_route_max_relative": route_residual,
        "only_J_changes": only_J_changes,
    }


def constitutive_gate() -> dict[str, Any]:
    """Verify the density-only mass-moment contract before root searches."""

    baseline = build_homogeneous_section()
    checks: list[dict[str, Any]] = []
    full_grid_positive = all(
        min((1.0 + float(value)) * RHO_0, (1.0 - float(value)) * RHO_0) > 0.0
        for value in np.concatenate((eta_grid(), xi_rho_grid()))
    )
    for eta in (0.0, 0.4, 0.8):
        for layout_id, layout, expected in (
            ("OUTER_HEAVY", OUTER_HEAVY_LAYOUT, 1.0 + 0.75 * eta),
            ("INNER_HEAVY", INNER_HEAVY_LAYOUT, 1.0 - 0.75 * eta),
        ):
            check = _section_check(
                build_eta_layout_section(layout, eta), baseline, expected
            )
            checks.append(
                {
                    "experiment_id": EXPERIMENT_A,
                    "configuration_id": layout_id,
                    "parameter_name": "eta",
                    "parameter_value": eta,
                    **check,
                }
            )
        inner = build_eta_layout_section(INNER_HEAVY_LAYOUT, eta)
        outer = build_eta_layout_section(OUTER_HEAVY_LAYOUT, eta)
        sum_residual = rlb2e._relative(
            inner.properties.J + outer.properties.J, 2.0 * baseline.properties.J
        )
        checks.append(
            {
                "experiment_id": EXPERIMENT_A,
                "configuration_id": CONFIG_ANTI_PHASE,
                "parameter_name": "eta",
                "parameter_value": eta,
                "status": (
                    "PASS" if sum_residual <= REDUCED_PROPERTY_TOLERANCE else "FAIL"
                ),
                "passed": sum_residual <= REDUCED_PROPERTY_TOLERANCE,
                "anti_phase_J_sum_relative": sum_residual,
            }
        )
    for xi_rho in (-0.8, -0.4, 0.0, 0.4, 0.8):
        check = _section_check(
            build_one_arm_layered_section(xi_rho),
            baseline,
            1.0 + 0.75 * xi_rho,
        )
        checks.append(
            {
                "experiment_id": EXPERIMENT_B,
                "configuration_id": CONFIG_ONE_ARM,
                "parameter_name": "xi_rho",
                "parameter_value": xi_rho,
                **check,
            }
        )
    maxima: dict[str, float] = {}
    for name in (
        "A_matrix_invariance_residual",
        "D_matrix_invariance_residual",
        "shear_matrix_invariance_residual",
        "A_beam_invariance_residual",
        "D_beam_invariance_residual",
        "S_beam_invariance_residual",
        "m_invariance_residual",
        "I0_invariance_residual",
        "B_relative",
        "I1_relative",
        "J_formula_residual",
        "reduction_route_max_relative",
        "anti_phase_J_sum_relative",
    ):
        maxima[name] = max(float(item.get(name, 0.0)) for item in checks)
    passed = bool(full_grid_positive and all(bool(item["passed"]) for item in checks))
    return {
        "status": "PASS" if passed else "FAIL",
        "passed": passed,
        "full_grid_density_positive": full_grid_positive,
        "checks": checks,
        "maximum_residuals": maxima,
        "Abeam0": baseline.properties.A,
        "Dbeam0": baseline.properties.D,
        "Sbeam0": baseline.properties.S,
        "m0": baseline.properties.m,
        "J0": baseline.properties.J,
        "I0": baseline.laminate.I0,
        "I2": baseline.laminate.I2,
        "tolerances": {
            "matrix_relative": MATRIX_RELATIVE_TOLERANCE,
            "symmetry_relative": SYMMETRY_RELATIVE_TOLERANCE,
            "reduced_property_relative": REDUCED_PROPERTY_TOLERANCE,
            "reduction_route_relative": REDUCTION_ROUTE_TOLERANCE,
        },
    }


def _section_row(
    *,
    experiment_id: str,
    parameter_name: str,
    parameter_value: float,
    configuration_id: str,
    arm_id: int,
    arm_role: str,
    section: SectionObjects,
    baseline: SectionObjects,
    expected_J_ratio: float,
) -> dict[str, Any]:
    check = _section_check(section, baseline, expected_J_ratio)
    laminate = section.laminate
    properties = section.properties
    return {
        "experiment_id": experiment_id,
        "parameter_name": parameter_name,
        "parameter_value": float(parameter_value),
        "configuration_id": configuration_id,
        "arm_id": int(arm_id),
        "arm_role": arm_role,
        "stack_bottom_to_top": list(section.layout),
        "ply_densities": list(section.ply_densities),
        "angle_deg": [float(ply.angle_deg) for ply in laminate.plies],
        "ply_thickness": [float(ply.thickness) for ply in laminate.plies],
        "z_interfaces": laminate.z_interfaces,
        "A_matrix": laminate.A,
        "B_matrix": laminate.B,
        "D_matrix": laminate.D,
        "shear_matrix_yz_xz": laminate.shear,
        "I0": laminate.I0,
        "I1": laminate.I1,
        "I2": laminate.I2,
        "A_beam": properties.A,
        "D_beam": properties.D,
        "S_beam": properties.S,
        "m": properties.m,
        "J": properties.J,
        "J_over_J0": check["J_over_J0"],
        "expected_J_over_J0": expected_J_ratio,
        "J_formula_residual": check["J_formula_residual"],
        "B_relative": check["B_relative"],
        "I1_relative": check["I1_relative"],
        "A_matrix_invariance_residual": check["A_matrix_invariance_residual"],
        "D_matrix_invariance_residual": check["D_matrix_invariance_residual"],
        "shear_matrix_invariance_residual": check[
            "shear_matrix_invariance_residual"
        ],
        "A_beam_invariance_residual": check["A_beam_invariance_residual"],
        "D_beam_invariance_residual": check["D_beam_invariance_residual"],
        "S_beam_invariance_residual": check["S_beam_invariance_residual"],
        "m_invariance_residual": check["m_invariance_residual"],
        "reduction_route_max_relative": check["reduction_route_max_relative"],
        "elastic_properties_identical_across_plies": check[
            "elastic_properties_identical_across_plies"
        ],
        "A_invariant": check["A_beam_invariance_residual"]
        <= REDUCED_PROPERTY_TOLERANCE,
        "D_invariant": check["D_beam_invariance_residual"]
        <= REDUCED_PROPERTY_TOLERANCE,
        "S_invariant": check["S_beam_invariance_residual"]
        <= REDUCED_PROPERTY_TOLERANCE,
        "m_invariant": check["m_invariance_residual"]
        <= REDUCED_PROPERTY_TOLERANCE,
        "only_J_changes": check["only_J_changes"],
        "constitutive_status": check["status"],
    }


def section_property_rows() -> list[dict[str, Any]]:
    """Return 408 logical arm rows for both experiments."""

    baseline = build_homogeneous_section()
    rows: list[dict[str, Any]] = []
    for configuration_id in CONFIGURATIONS_A:
        for eta_value in eta_grid():
            eta = float(eta_value)
            sections = build_experiment_A_sections(configuration_id, eta)
            layouts = CONFIGURATION_LAYOUTS_A[configuration_id]
            for arm_id, (layout, section) in enumerate(
                zip(layouts, sections, strict=True), start=1
            ):
                outer = layout == OUTER_HEAVY_LAYOUT
                expected = 1.0 + (0.75 * eta if outer else -0.75 * eta)
                rows.append(
                    _section_row(
                        experiment_id=EXPERIMENT_A,
                        parameter_name="eta",
                        parameter_value=eta,
                        configuration_id=configuration_id,
                        arm_id=arm_id,
                        arm_role="OUTER_HEAVY" if outer else "INNER_HEAVY",
                        section=section,
                        baseline=baseline,
                        expected_J_ratio=expected,
                    )
                )
    for xi_value in xi_rho_grid():
        xi_rho = float(xi_value)
        layered, homogeneous = build_experiment_B_sections(xi_rho)
        rows.append(
            _section_row(
                experiment_id=EXPERIMENT_B,
                parameter_name="xi_rho",
                parameter_value=xi_rho,
                configuration_id=CONFIG_ONE_ARM,
                arm_id=1,
                arm_role="LAYERED_DENSITY",
                section=layered,
                baseline=baseline,
                expected_J_ratio=1.0 + 0.75 * xi_rho,
            )
        )
        rows.append(
            _section_row(
                experiment_id=EXPERIMENT_B,
                parameter_name="xi_rho",
                parameter_value=xi_rho,
                configuration_id=CONFIG_ONE_ARM,
                arm_id=2,
                arm_role="HOMOGENEOUS_REFERENCE",
                section=homogeneous,
                baseline=baseline,
                expected_J_ratio=1.0,
            )
        )
    return rows


def audit_section_property_rows(
    rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Validate the exact logical-arm inventory used by the J diagnostic."""

    expected = {
        (
            EXPERIMENT_A,
            configuration_id,
            round(float(value), 10),
            arm_id,
        )
        for configuration_id in CONFIGURATIONS_A
        for value in eta_grid()
        for arm_id in (1, 2)
    }
    expected.update(
        {
            (EXPERIMENT_B, CONFIG_ONE_ARM, round(float(value), 10), arm_id)
            for value in xi_rho_grid()
            for arm_id in (1, 2)
        }
    )
    keys: list[tuple[str, str, float, int]] = []
    invalid: list[str] = []
    for row in rows:
        try:
            key = (
                str(row["experiment_id"]),
                str(row["configuration_id"]),
                round(float(row["parameter_value"]), 10),
                int(row["arm_id"]),
            )
            keys.append(key)
            if str(row["constitutive_status"]) != "PASS":
                invalid.append(f"{key}:constitutive_status")
            for field in ("J", "J_over_J0", "expected_J_over_J0"):
                if not math.isfinite(float(row[field])):
                    invalid.append(f"{key}:{field}")
        except (KeyError, TypeError, ValueError) as exc:
            invalid.append(f"malformed-row:{type(exc).__name__}")
    duplicate_keys = sorted({key for key in keys if keys.count(key) > 1})
    actual = set(keys)
    missing = sorted(expected - actual)
    unexpected = sorted(actual - expected)
    passed = not duplicate_keys and not missing and not unexpected and not invalid
    return {
        "status": "PASS" if passed else "FAIL",
        "row_count": len(rows),
        "expected_row_count": len(expected),
        "duplicate_keys": duplicate_keys,
        "missing_keys": missing,
        "unexpected_keys": unexpected,
        "invalid_rows": invalid,
    }


if not (rlb2f.K_PLOT == rlb2e.K_PLOT == K_PLOT == 8):
    raise RuntimeError("RLB-2G requires the frozen plotted-root count K=8.")
if not (rlb2f.K_GUARD == rlb2e.K_GUARD == K_GUARD == 9):
    raise RuntimeError("RLB-2G requires the frozen root-9 guard contract.")


def contract_payload() -> dict[str, Any]:
    return {
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "frequency_map_policy_instances": {
            "experiment_A": POLICY_EXPERIMENT_A,
            "experiment_B": POLICY_EXPERIMENT_B,
        },
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
            "plies_per_arm": 4,
            "all_ply_angles_deg": 0.0,
            "K": K,
            "outer_clamps": True,
            "joint": "frozen ideal rigid RLB joint",
        },
        "material_M0": base_material_contract(),
        "density_contrast": {
            "rho_H": "(1+eta)*rho0",
            "rho_L": "(1-eta)*rho0",
            "rho_outer_signed": "(1+xi_rho)*rho0",
            "rho_inner_signed": "(1-xi_rho)*rho0",
            "elastic_properties_scaled": False,
            "minimum_grid_density": 0.2 * RHO_0,
        },
        "mass_moment_contract": {
            "m0": M0,
            "J0": J0,
            "J_outer_over_J0": "1+3*eta/4",
            "J_inner_over_J0": "1-3*eta/4",
            "J_layered_over_J0": "1+3*xi_rho/4",
            "anti_phase_J1_plus_J2": "2*J0",
        },
        "experiment_A": {
            "parameter": "eta",
            "grid": [float(value) for value in eta_grid()],
            "configurations": {
                key: [list(layout) for layout in value]
                for key, value in CONFIGURATION_LAYOUTS_A.items()
            },
        },
        "experiment_B": {
            "parameter": "xi_rho",
            "grid": [float(value) for value in xi_rho_grid()],
            "continuation_anchor": 0.0,
            "continuation_paths": [
                [float(value) for value in path]
                for path in xi_rho_continuation_paths()
            ],
            "arm1": "layered density [outer,inner,inner,outer]",
            "arm2": "homogeneous rho0 reference",
        },
        "shared_zero_contrast_anchor": {
            "logical_group_count": 4,
            "physical_solve_count": 1,
            "clone_count": 3,
        },
        "normalization": {
            "Omega": "omega*l^2*sqrt(rho0*A0/(E0*Iy0))",
            "Lambda": "sqrt(Omega)",
            "E0": 1.0,
            "rho0": RHO_0,
            "b0": WIDTH,
            "h0": THICKNESS,
            "l": L_REFERENCE,
            "A0": A_REFERENCE,
            "Iy0": IY_REFERENCE,
            "Omega_per_omega": OMEGA_TO_OMEGA_SCALE,
            "density_or_J_dependent": False,
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
            "RLB2E_source": RLB2E_SCRIPT_PATH,
            "RLB2F_guard_source": RLB2F_SCRIPT_PATH,
            "requested_roots": K_PLOT,
            "guard_roots": 1,
            "required_slots": K_GUARD,
            "physics_formulas_copied": False,
        },
        "explicit_exclusions": [
            "roots_10_and_above",
            "universal_spectral_runner",
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


def _sections_for_point(
    experiment_id: str,
    configuration_id: str,
    parameter_value: float,
    *,
    swapped: bool = False,
) -> tuple[SectionObjects, SectionObjects]:
    if experiment_id == EXPERIMENT_A:
        return build_experiment_A_sections(
            configuration_id, parameter_value, swapped=swapped
        )
    if experiment_id == EXPERIMENT_B and configuration_id == CONFIG_ONE_ARM:
        return build_experiment_B_sections(parameter_value, swapped=swapped)
    raise ValueError(
        f"Unknown RLB-2G point: {experiment_id!r}, {configuration_id!r}."
    )


def make_matrix_provider(
    experiment_id: str,
    configuration_id: str,
    parameter_value: float,
    *,
    swapped: bool = False,
) -> tuple[Any, dict[str, Any]]:
    """Build the frozen coupled boundary matrix from density-only sections."""

    _beam, coupled = rlb2e._physics_modules()
    arm1, arm2 = _sections_for_point(
        experiment_id, configuration_id, parameter_value, swapped=swapped
    )
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
        raise RuntimeError("RLB-2G provider differs from the frozen public builder.")
    return provider, {
        "experiment_id": experiment_id,
        "configuration_id": configuration_id,
        "parameter_value": float(parameter_value),
        "swapped_arms": bool(swapped),
        "identical_arm_map_reused": identical,
        "cached_vs_public_builder_max_abs": direct_residual,
        "production_modules_only": True,
    }


def _parameter_name(experiment_id: str) -> str:
    if experiment_id == EXPERIMENT_A:
        return "eta"
    if experiment_id == EXPERIMENT_B:
        return "xi_rho"
    raise ValueError(f"Unknown experiment: {experiment_id!r}.")


def _root_rows(
    experiment_id: str,
    configuration_id: str,
    parameter_value: float,
    slots: Sequence[Any],
    *,
    solve_id: str,
    physical_solve_id: str,
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
    shared_zero: bool = False,
) -> tuple[dict[str, Any], ...]:
    windows = (
        None
        if predicted is None
        else rlb2e.local_search_windows(
            predicted, guard_right_width=guard_locator_right_width
        )
    )
    parameter_name = _parameter_name(experiment_id)
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
            f"{experiment_id}__{configuration_id}__{parameter_name}_"
            f"{float(parameter_value):+.6f}__{grid_kind}__p{position:02d}__"
            f"{repair_id or 'base'}"
        )
        rows.append(
            {
                "row_id": row_id,
                "experiment_id": experiment_id,
                "configuration_id": configuration_id,
                "parameter_name": parameter_name,
                "parameter_value": float(parameter_value),
                "eta": float(parameter_value) if experiment_id == EXPERIMENT_A else "",
                "xi_rho": (
                    float(parameter_value) if experiment_id == EXPERIMENT_B else ""
                ),
                "grid_kind": grid_kind,
                "continuation_leg": continuation_leg,
                "solve_id": solve_id,
                "physical_solve_id": physical_solve_id,
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
                "shared_zero_contrast_anchor_reused": shared_zero,
                "shared_zero_contrast_source": (
                    "RLB-2G_SHARED_ZERO" if shared_zero else ""
                ),
                "is_canonical_plot_source": True,
                "supersedes_row_id": "",
                "repair_id": repair_id,
                "roots_above_9_computed": False,
            }
        )
    return tuple(rows)


def solve_point(
    experiment_id: str,
    configuration_id: str,
    parameter_value: float,
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
    """Solve one density-layout point with the frozen RLB detector/refiner."""

    global ROOT_CALCULATION_COUNT
    ROOT_CALCULATION_COUNT += 1
    started = time.perf_counter()
    parameter_name = _parameter_name(experiment_id)
    solve_id = (
        f"{experiment_id}__{configuration_id}__{parameter_name}_"
        f"{float(parameter_value):+.6f}"
        + ("__ARM_SWAP" if swapped else "")
    )
    transaction_id = hashlib.sha256(
        f"{STAGE_ID}|{solve_id}|{grid_kind}|{repair_id}".encode("utf-8")
    ).hexdigest()[:20].upper()
    provider, _metadata = make_matrix_provider(
        experiment_id, configuration_id, parameter_value, swapped=swapped
    )
    counted = rlb2e.CountedProvider(provider)
    policy = rlb2e._rlb2e_search_policy()
    predicted: FloatArray | None = None
    fallback_used = False
    used_guard_locator_right_width = guard_locator_right_width
    solve_mode = ""
    if not force_anchor and previous is not None:
        predicted = rlb2e.hold_secant_predictor(
            parameter_value,
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
        if not accepted and rlb2f._guard_detector_duplicate(slots, policy):
            (
                canonical,
                rejected,
                slots,
                search_right,
                extra_refinements,
                predicted,
            ) = rlb2f._dense_guard_reconciliation(
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
            if not accepted and rlb2f._guard_detector_duplicate(slots, policy):
                (
                    canonical,
                    rejected,
                    slots,
                    search_right,
                    extra_refinements,
                    predicted,
                ) = rlb2f._dense_guard_reconciliation(
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
            rlb2e._progressive_anchor_candidates(
                counted, policy, solve_id=solve_id
            )
        )
        accepted, quality = rlb2e._point_is_acceptable(
            canonical, rejected, slots, search_right, policy
        )
        solve_mode = "PROGRESSIVE_ANCHOR"
        if not accepted and rlb2f._guard_detector_duplicate(slots, policy):
            (
                canonical,
                rejected,
                slots,
                search_right,
                extra_refinements,
                predicted,
            ) = rlb2f._dense_guard_reconciliation(
                counted, policy, slots, solve_id=solve_id
            )
            used_guard_locator_right_width = 0.2
            refinements += extra_refinements
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
        experiment_id,
        configuration_id,
        parameter_value,
        slots,
        solve_id=solve_id,
        physical_solve_id=solve_id,
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
        experiment_id=experiment_id,
        configuration_id=configuration_id,
        parameter_name=parameter_name,
        parameter_value=float(parameter_value),
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
) -> dict[tuple[str, float], list[Mapping[str, Any]]]:
    groups: dict[tuple[str, float], list[Mapping[str, Any]]] = {}
    for row in rows:
        if str(row.get("grid_kind")) != "BASE":
            continue
        key = (
            str(row["configuration_id"]),
            round(float(row["parameter_value"]), 10),
        )
        groups.setdefault(key, []).append(row)
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
) -> dict[tuple[str, float], list[Mapping[str, Any]]]:
    return {
        key: group
        for key, group in _base_group_index(rows).items()
        if _base_group_is_acceptable(group)
    }


def _expected_keys(experiment_id: str) -> set[tuple[str, float]]:
    if experiment_id == EXPERIMENT_A:
        return {
            (configuration_id, round(float(eta), 10))
            for configuration_id in CONFIGURATIONS_A
            for eta in eta_grid()
        }
    if experiment_id == EXPERIMENT_B:
        return {
            (CONFIG_ONE_ARM, round(float(xi_rho), 10))
            for xi_rho in xi_rho_grid()
        }
    raise ValueError(f"Unknown experiment: {experiment_id!r}.")


def audit_spectrum_rows(
    rows: Sequence[Mapping[str, Any]], experiment_id: str
) -> dict[str, Any]:
    groups = _base_group_index(rows)
    complete = _complete_base_group_index(rows)
    expected = _expected_keys(experiment_id)
    duplicates: list[str] = []
    incomplete: list[str] = []
    quality_failures: list[str] = []
    for key, group in groups.items():
        positions = [int(row["sorted_position"]) for row in group]
        label = f"{key[0]}:{key[1]:+.2f}"
        if len(positions) != len(set(positions)):
            duplicates.append(label)
        if not _base_group_has_exact_positions(group):
            incomplete.append(label)
        elif not _base_group_is_acceptable(group):
            quality_failures.append(label)
    base = [row for row in rows if str(row.get("grid_kind")) == "BASE"]
    above = [row for row in rows if int(row["sorted_position"]) > K_GUARD]
    duplicate_row_ids = len(rows) - len({str(row["row_id"]) for row in rows})
    canonical_counts: dict[tuple[str, float, int], int] = {}
    for row in rows:
        if _as_bool(row.get("is_canonical_plot_source", False)):
            key = (
                str(row["configuration_id"]),
                round(float(row["parameter_value"]), 10),
                int(row["sorted_position"]),
            )
            canonical_counts[key] = canonical_counts.get(key, 0) + 1
    expected_canonical = {
        (configuration_id, value, position)
        for configuration_id, value in expected
        for position in range(1, K_GUARD + 1)
    }
    canonical_failures = [
        f"{key[0]}:{key[1]:+.2f}:p{key[2]}"
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
        and len(base) == len(expected) * K_GUARD
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "base_group_count": len(complete),
        "base_row_count": len(base),
        "missing_groups": [f"{key[0]}:{key[1]:+.2f}" for key in missing],
        "extra_groups": [f"{key[0]}:{key[1]:+.2f}" for key in extra],
        "duplicate_groups": duplicates,
        "incomplete_groups": incomplete,
        "physical_quality_failures": quality_failures,
        "duplicate_row_id_count": duplicate_row_ids,
        "canonical_source_failures": canonical_failures,
        "roots_above_guard_count": len(above),
    }


def _canonical_group(
    rows: Sequence[Mapping[str, Any]], configuration_id: str, parameter_value: float
) -> list[Mapping[str, Any]]:
    selected = [
        row
        for row in rows
        if str(row["configuration_id"]) == configuration_id
        and round(float(row["parameter_value"]), 10)
        == round(float(parameter_value), 10)
        and _as_bool(row.get("is_canonical_plot_source", True))
    ]
    selected.sort(key=lambda row: int(row["sorted_position"]))
    if [int(row["sorted_position"]) for row in selected] != list(
        range(1, K_GUARD + 1)
    ):
        raise RuntimeError(
            f"Incomplete canonical group at {configuration_id}, "
            f"parameter={parameter_value:+.2f}."
        )
    return selected


def _roots_for(
    rows: Sequence[Mapping[str, Any]], configuration_id: str, parameter_value: float
) -> FloatArray:
    return np.asarray(
        [
            float(row["Omega"])
            for row in _canonical_group(rows, configuration_id, parameter_value)
        ],
        dtype=float,
    )


def _sort_rows(rows: list[dict[str, Any]]) -> None:
    rows.sort(
        key=lambda row: (
            str(row["configuration_id"]),
            float(row["parameter_value"]),
            0 if str(row["grid_kind"]) == "BASE" else 1,
            int(row["sorted_position"]),
            str(row.get("repair_id", "")),
        )
    )


def _spectrum_path(output_dir: Path, experiment_id: str) -> Path:
    if experiment_id == EXPERIMENT_A:
        return Path(output_dir) / SPECTRUM_A_FILENAME
    if experiment_id == EXPERIMENT_B:
        return Path(output_dir) / SPECTRUM_B_FILENAME
    raise ValueError(f"Unknown experiment: {experiment_id!r}.")


def _read_spectrum(output_dir: Path, experiment_id: str) -> list[dict[str, Any]]:
    path = _spectrum_path(output_dir, experiment_id)
    return [] if not path.is_file() else rlb2e._read_csv(path)


def _write_spectrum(
    output_dir: Path, experiment_id: str, rows: Sequence[Mapping[str, Any]]
) -> None:
    rlb2e._atomic_write_csv(
        _spectrum_path(output_dir, experiment_id), rows, SPECTRUM_FIELDS
    )


def _write_point_transaction(
    output_dir: Path,
    existing_rows: Sequence[Mapping[str, Any]],
    solution: PointSolution,
) -> list[dict[str, Any]]:
    rows = [dict(row) for row in existing_rows]
    key = (solution.configuration_id, round(solution.parameter_value, 10))
    group = _base_group_index(rows).get(key, [])
    if group:
        if _base_group_is_acceptable(group):
            return rows
        rows = [
            row
            for row in rows
            if not (
                str(row["configuration_id"]) == solution.configuration_id
                and round(float(row["parameter_value"]), 10) == key[1]
                and str(row["grid_kind"]) in {"BASE", "LOCAL_REFINEMENT"}
            )
        ]
    rows.extend(dict(row) for row in solution.rows)
    _sort_rows(rows)
    _write_spectrum(output_dir, solution.experiment_id, rows)
    return rows


def _clone_shared_zero_rows(
    source_rows: Sequence[Mapping[str, Any]],
    experiment_id: str,
    configuration_id: str,
) -> tuple[dict[str, Any], ...]:
    parameter_name = _parameter_name(experiment_id)
    transaction_id = hashlib.sha256(
        f"{STAGE_ID}|SHARED_ZERO|{experiment_id}|{configuration_id}".encode("utf-8")
    ).hexdigest()[:20].upper()
    result: list[dict[str, Any]] = []
    for source in sorted(source_rows, key=lambda row: int(row["sorted_position"])):
        row = dict(source)
        position = int(row["sorted_position"])
        row.update(
            {
                "row_id": (
                    f"{experiment_id}__{configuration_id}__{parameter_name}_"
                    f"+0.000000__BASE__p{position:02d}__base"
                ),
                "experiment_id": experiment_id,
                "configuration_id": configuration_id,
                "parameter_name": parameter_name,
                "parameter_value": 0.0,
                "eta": 0.0 if experiment_id == EXPERIMENT_A else "",
                "xi_rho": 0.0 if experiment_id == EXPERIMENT_B else "",
                "grid_kind": "BASE",
                "continuation_leg": "SHARED_ZERO_ANCHOR",
                "solve_id": f"{experiment_id}__{configuration_id}__SHARED_ZERO_CLONE",
                "physical_solve_id": "RLB-2G_SHARED_ZERO_PHYSICAL_SOLVE",
                "transaction_id": transaction_id,
                "shared_zero_contrast_anchor_reused": True,
                "shared_zero_contrast_source": "RLB-2G_SHARED_ZERO",
                "is_canonical_plot_source": True,
                "supersedes_row_id": "",
                "repair_id": "",
            }
        )
        result.append(row)
    return tuple(result)


def _solution_record(solution: PointSolution, *, benchmark: bool) -> dict[str, Any]:
    return {
        "experiment_id": solution.experiment_id,
        "configuration_id": solution.configuration_id,
        "parameter_name": solution.parameter_name,
        "parameter_value": solution.parameter_value,
        "benchmark": benchmark,
        "continuation_leg": solution.continuation_leg,
        "swapped_arms": solution.swapped_arms,
        "wall_time_seconds": solution.wall_time_seconds,
        "peak_rss_bytes": solution.peak_rss_bytes,
        "determinant_evaluations": solution.determinant_evaluations,
        "sigma_evaluations": solution.sigma_evaluations,
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


def _point_record_key(record: Mapping[str, Any]) -> tuple[str, str, float, bool]:
    return (
        str(record.get("experiment_id", "")),
        str(record.get("configuration_id", "")),
        round(float(record.get("parameter_value", math.nan)), 10),
        bool(record.get("swapped_arms", False)),
    )


def _checkpoint_payload(
    rows_A: Sequence[Mapping[str, Any]],
    rows_B: Sequence[Mapping[str, Any]],
    point_records: Sequence[Mapping[str, Any]],
    *,
    constitutive: Mapping[str, Any],
    started_at: str,
    benchmark_status: str,
) -> dict[str, Any]:
    complete_A = _complete_base_group_index(rows_A)
    complete_B = _complete_base_group_index(rows_B)
    missing_A = [
        {"configuration_id": configuration_id, "eta": float(eta)}
        for configuration_id, eta in sorted(_expected_keys(EXPERIMENT_A))
        if (configuration_id, eta) not in complete_A
    ]
    missing_B = [
        float(value)
        for configuration_id, value in sorted(_expected_keys(EXPERIMENT_B))
        if (configuration_id, value) not in complete_B
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
        "experiment_A": {
            "completed_groups": len(complete_A),
            "missing_groups": missing_A,
        },
        "experiment_B": {
            "completed_groups": len(complete_B),
            "missing_points": missing_B,
            "continuation_paths": [
                [float(value) for value in path]
                for path in xi_rho_continuation_paths()
            ],
        },
        "logical_completed_groups": len(complete_A) + len(complete_B),
        "logical_completed_base_rows": (
            len(complete_A) + len(complete_B)
        )
        * K_GUARD,
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


def _write_checkpoint(
    output_dir: Path,
    rows_A: Sequence[Mapping[str, Any]],
    rows_B: Sequence[Mapping[str, Any]],
    point_records: Sequence[Mapping[str, Any]],
    *,
    constitutive: Mapping[str, Any],
    started_at: str,
    benchmark_status: str,
) -> None:
    rlb2e._atomic_write_json(
        Path(output_dir) / CHECKPOINT_FILENAME,
        _checkpoint_payload(
            rows_A,
            rows_B,
            point_records,
            constitutive=constitutive,
            started_at=started_at,
            benchmark_status=benchmark_status,
        ),
    )


def _existing_point_records(output_dir: Path) -> list[dict[str, Any]]:
    path = Path(output_dir) / CHECKPOINT_FILENAME
    if not path.is_file():
        return []
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("contract_sha256") != contract_hash():
        raise RuntimeError("Existing RLB-2G checkpoint has an incompatible contract.")
    return [dict(item) for item in payload.get("point_records", [])]


def _append_record_once(
    point_records: list[dict[str, Any]], record: Mapping[str, Any]
) -> None:
    key = _point_record_key(record)
    if key not in {_point_record_key(item) for item in point_records}:
        point_records.append(dict(record))


def ensure_shared_zero_anchor(
    output_dir: Path,
    rows_A: list[dict[str, Any]],
    rows_B: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Solve the homogeneous spectrum once and bind four logical groups."""

    existing_groups: list[list[Mapping[str, Any]]] = []
    for configuration_id in CONFIGURATIONS_A:
        group = _base_group_index(rows_A).get((configuration_id, 0.0), [])
        if _base_group_is_acceptable(group):
            existing_groups.append(group)
    group_B = _base_group_index(rows_B).get((CONFIG_ONE_ARM, 0.0), [])
    if _base_group_is_acceptable(group_B):
        existing_groups.append(group_B)
    solution: PointSolution | None = None
    if existing_groups:
        source = sorted(existing_groups[0], key=lambda row: int(row["sorted_position"]))
        reference = np.asarray([float(row["Omega"]) for row in source])
        for group in existing_groups[1:]:
            candidate = np.asarray(
                [
                    float(row["Omega"])
                    for row in sorted(group, key=lambda row: int(row["sorted_position"]))
                ]
            )
            if not np.array_equal(reference, candidate):
                raise RuntimeError("Existing RLB-2G zero-anchor groups disagree.")
        source_rows = source
        reused_existing = True
        wall_time = 0.0
        evaluations = 0
        peak = rlb2e._peak_rss_bytes()
    else:
        solution = solve_point(
            EXPERIMENT_A,
            CONFIG_BOTH_OUTER,
            0.0,
            force_anchor=True,
            continuation_leg="SHARED_ZERO_ANCHOR",
        )
        source_rows = list(solution.rows)
        reused_existing = False
        wall_time = solution.wall_time_seconds
        evaluations = solution.determinant_evaluations
        peak = solution.peak_rss_bytes
        record = _solution_record(solution, benchmark=True)
        record["shared_zero_physical_solve"] = True
        _append_record_once(point_records, record)
    for configuration_id in CONFIGURATIONS_A:
        key = (configuration_id, 0.0)
        if not _base_group_is_acceptable(_base_group_index(rows_A).get(key, [])):
            rows_A = [
                row
                for row in rows_A
                if not (
                    str(row.get("configuration_id")) == configuration_id
                    and round(float(row.get("parameter_value", math.nan)), 10) == 0.0
                )
            ]
            rows_A.extend(
                _clone_shared_zero_rows(source_rows, EXPERIMENT_A, configuration_id)
            )
    if not _base_group_is_acceptable(
        _base_group_index(rows_B).get((CONFIG_ONE_ARM, 0.0), [])
    ):
        rows_B = [
            row
            for row in rows_B
            if round(float(row.get("parameter_value", math.nan)), 10) != 0.0
        ]
        rows_B.extend(
            _clone_shared_zero_rows(source_rows, EXPERIMENT_B, CONFIG_ONE_ARM)
        )
    _sort_rows(rows_A)
    _sort_rows(rows_B)
    _write_spectrum(output_dir, EXPERIMENT_A, rows_A)
    _write_spectrum(output_dir, EXPERIMENT_B, rows_B)
    return rows_A, rows_B, {
        "experiment_id": "SHARED_ZERO_ANCHOR",
        "configuration_id": "HOMOGENEOUS",
        "parameter_name": "zero_contrast",
        "parameter_value": 0.0,
        "benchmark": True,
        "shared_zero_physical_solve": True,
        "logical_group_count": 4,
        "physical_solve_count": 0 if reused_existing else 1,
        "wall_time_seconds": wall_time,
        "peak_rss_bytes": peak,
        "determinant_evaluations": evaluations,
        "sigma_evaluations": evaluations,
        "actual_Omega_max": float(source_rows[-1]["search_right_Omega"]),
        "quality_status": "PASS",
        "reused_existing": reused_existing,
        "roots": [
            {
                "sorted_position": int(row["sorted_position"]),
                "Omega": float(row["Omega"]),
                "Lambda": float(row["Lambda"]),
            }
            for row in source_rows
        ],
    }


def _ensure_benchmark_point(
    output_dir: Path,
    rows: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    *,
    experiment_id: str,
    configuration_id: str,
    parameter_value: float,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    key = (configuration_id, round(float(parameter_value), 10))
    group = _complete_base_group_index(rows).get(key)
    if group is not None:
        existing = next(
            (
                record
                for record in point_records
                if _point_record_key(record)
                == (experiment_id, configuration_id, key[1], False)
            ),
            None,
        )
        if (
            existing is not None
            and math.isfinite(float(existing.get("wall_time_seconds", math.nan)))
            and float(existing.get("wall_time_seconds", 0.0)) > 0.0
            and int(existing.get("determinant_evaluations", 0)) > 0
        ):
            return rows, dict(existing)
        raise RuntimeError(
            "A durable nonzero benchmark group exists without measured timing "
            "and determinant-evaluation evidence; conservative ETA cannot be "
            f"reconstructed for {experiment_id}:{configuration_id}:"
            f"{float(parameter_value):+.6f}. No output was rewritten."
        )
    solution = solve_point(
        experiment_id,
        configuration_id,
        parameter_value,
        force_anchor=True,
        continuation_leg="BENCHMARK_ANCHOR",
    )
    rows = _write_point_transaction(output_dir, rows, solution)
    record = _solution_record(solution, benchmark=True)
    _append_record_once(point_records, record)
    return rows, record


def run_benchmarks(
    output_dir: Path,
    rows_A: list[dict[str, Any]],
    rows_B: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    started_at: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    rows_A, rows_B, zero_record = ensure_shared_zero_anchor(
        output_dir, rows_A, rows_B, point_records
    )
    records: list[dict[str, Any]] = [zero_record]
    for experiment_id, configuration_id, parameter_value in (
        (EXPERIMENT_A, CONFIG_BOTH_OUTER, 0.8),
        (EXPERIMENT_A, CONFIG_ANTI_PHASE, 0.8),
        (EXPERIMENT_B, CONFIG_ONE_ARM, -0.8),
        (EXPERIMENT_B, CONFIG_ONE_ARM, 0.8),
    ):
        if experiment_id == EXPERIMENT_A:
            rows_A, record = _ensure_benchmark_point(
                output_dir,
                rows_A,
                point_records,
                experiment_id=experiment_id,
                configuration_id=configuration_id,
                parameter_value=parameter_value,
            )
        else:
            rows_B, record = _ensure_benchmark_point(
                output_dir,
                rows_B,
                point_records,
                experiment_id=experiment_id,
                configuration_id=configuration_id,
                parameter_value=parameter_value,
            )
        records.append(record)
        _write_checkpoint(
            output_dir,
            rows_A,
            rows_B,
            point_records,
            constitutive=constitutive,
            started_at=started_at,
            benchmark_status="RUNNING",
        )
    nonzero_times = [
        float(record.get("wall_time_seconds", 0.0)) for record in records[1:]
    ]
    maximum = max(nonzero_times, default=0.0)
    remaining = 201 - 5
    eta_seconds = 1.5 * maximum * remaining
    payload = {
        "schema_version": 1,
        "stage_id": STAGE_ID,
        "physical_anchors": records,
        "logical_anchor_bindings": {
            "shared_zero": 4,
            "nonzero_benchmarks": 4,
        },
        "total_logical_groups": 204,
        "total_unique_physical_base_points": 201,
        "physical_anchor_count": 5,
        "logical_anchor_group_count": 8,
        "remaining_unique_physical_points": remaining,
        "maximum_nonzero_anchor_wall_time_seconds": maximum,
        "eta_formula": "1.5*max(nonzero_anchor_wall_time)*196",
        "conservative_eta_seconds": eta_seconds,
        "eta_limit_seconds": ETA_LIMIT_SECONDS,
        "production_run_permitted": eta_seconds <= ETA_LIMIT_SECONDS,
        "shared_zero_calculated_once": True,
    }
    rlb2e._atomic_write_json(Path(output_dir) / BENCHMARK_FILENAME, payload)
    _write_checkpoint(
        output_dir,
        rows_A,
        rows_B,
        point_records,
        constitutive=constitutive,
        started_at=started_at,
        benchmark_status="PASS" if payload["production_run_permitted"] else "FAIL",
    )
    return rows_A, rows_B, payload


def _record_recovered_groups(
    rows: Sequence[Mapping[str, Any]],
    experiment_id: str,
    point_records: list[dict[str, Any]],
) -> None:
    recorded = {_point_record_key(record) for record in point_records}
    for (configuration_id, parameter_value), group in sorted(
        _complete_base_group_index(rows).items()
    ):
        key = (experiment_id, configuration_id, parameter_value, False)
        if key in recorded or parameter_value == 0.0:
            continue
        ordered = sorted(group, key=lambda row: int(row["sorted_position"]))
        record = {
            "experiment_id": experiment_id,
            "configuration_id": configuration_id,
            "parameter_name": _parameter_name(experiment_id),
            "parameter_value": parameter_value,
            "benchmark": False,
            "continuation_leg": "RECOVERED_DURABLE_BASE_TRANSACTION",
            "swapped_arms": False,
            "wall_time_seconds": 0.0,
            "peak_rss_bytes": rlb2e._peak_rss_bytes(),
            "determinant_evaluations": 0,
            "sigma_evaluations": 0,
            "actual_Omega_max": float(ordered[-1]["search_right_Omega"]),
            "solve_mode": "RECOVERED_DURABLE_BASE_TRANSACTION",
            "fallback_used": False,
            "unresolved_candidates_below_root9": 0,
            "roots_above_9_computed": False,
        }
        _append_record_once(point_records, record)
        recorded.add(key)


def _complete_path(
    output_dir: Path,
    rows: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    *,
    experiment_id: str,
    configuration_id: str,
    path: Sequence[float],
    continuation_leg: str,
    constitutive: Mapping[str, Any],
    started_at: str,
    rows_A: list[dict[str, Any]],
    rows_B: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    history: list[tuple[float, FloatArray]] = []
    for raw_value in path:
        value = float(raw_value)
        key = (configuration_id, round(value, 10))
        complete = _complete_base_group_index(rows)
        if key in complete:
            roots = _roots_for(rows, configuration_id, value)
            history.append((value, roots))
            history = history[-2:]
            continue
        previous = history[-1] if history else None
        second_previous = history[-2] if len(history) > 1 else None
        solution = solve_point(
            experiment_id,
            configuration_id,
            value,
            previous=previous,
            second_previous=second_previous,
            continuation_leg=continuation_leg,
        )
        rows = _write_point_transaction(output_dir, rows, solution)
        _append_record_once(point_records, _solution_record(solution, benchmark=False))
        roots = np.asarray([float(row["Omega"]) for row in solution.rows])
        history.append((value, roots))
        history = history[-2:]
        if experiment_id == EXPERIMENT_A:
            rows_A = rows
        else:
            rows_B = rows
        _write_checkpoint(
            output_dir,
            rows_A,
            rows_B,
            point_records,
            constitutive=constitutive,
            started_at=started_at,
            benchmark_status="PASS",
        )
    return rows, rows_A, rows_B


def complete_missing_points(
    output_dir: Path,
    rows_A: list[dict[str, Any]],
    rows_B: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    started_at: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    for configuration_id in CONFIGURATIONS_A:
        rows_A, rows_A, rows_B = _complete_path(
            output_dir,
            rows_A,
            point_records,
            experiment_id=EXPERIMENT_A,
            configuration_id=configuration_id,
            path=[float(value) for value in eta_grid()],
            continuation_leg=f"ETA_POSITIVE__{configuration_id}",
            constitutive=constitutive,
            started_at=started_at,
            rows_A=rows_A,
            rows_B=rows_B,
        )
    negative, positive = xi_rho_continuation_paths()
    for leg_name, leg in (("XI_RHO_NEGATIVE", negative), ("XI_RHO_POSITIVE", positive)):
        rows_B, rows_A, rows_B = _complete_path(
            output_dir,
            rows_B,
            point_records,
            experiment_id=EXPERIMENT_B,
            configuration_id=CONFIG_ONE_ARM,
            path=[float(value) for value in leg],
            continuation_leg=leg_name,
            constitutive=constitutive,
            started_at=started_at,
            rows_A=rows_A,
            rows_B=rows_B,
        )
    return rows_A, rows_B


def _recover_interrupted_finalization(
    output_dir: Path,
    experiment_id: str,
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    if not any(str(row.get("grid_kind")) == "LOCAL_REFINEMENT" for row in rows):
        return [dict(row) for row in rows]
    recovered = [dict(row) for row in rows if str(row.get("grid_kind")) == "BASE"]
    for row in recovered:
        row["is_canonical_plot_source"] = True
        row["point_status"] = "PASS"
    _sort_rows(recovered)
    _write_spectrum(output_dir, experiment_id, recovered)
    return recovered


def _curve_specs(
    rows_A: Sequence[Mapping[str, Any]], rows_B: Sequence[Mapping[str, Any]]
) -> list[tuple[str, str, str, list[float], Sequence[Mapping[str, Any]]]]:
    result = [
        (
            EXPERIMENT_A,
            configuration_id,
            "eta",
            [float(value) for value in eta_grid()],
            rows_A,
        )
        for configuration_id in CONFIGURATIONS_A
    ]
    result.append(
        (
            EXPERIMENT_B,
            CONFIG_ONE_ARM,
            "xi_rho",
            [float(value) for value in xi_rho_grid()],
            rows_B,
        )
    )
    return result


def neighbour_audit_rows(
    rows_A: Sequence[Mapping[str, Any]], rows_B: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    """Apply the unchanged RLB-2E median-plus-eight-MAD trigger."""

    result: list[dict[str, Any]] = []
    for experiment_id, configuration_id, parameter_name, raw_grid, rows in _curve_specs(
        rows_A, rows_B
    ):
        grid = [round(float(value), 10) for value in raw_grid]
        spectra = {
            value: np.sqrt(_roots_for(rows, configuration_id, value)[:K_PLOT])
            for value in grid
        }
        gap_flags: set[tuple[float, int]] = set()
        for lower_position in range(1, K_PLOT):
            gaps = np.asarray(
                [
                    spectra[value][lower_position]
                    - spectra[value][lower_position - 1]
                    for value in grid
                ],
                dtype=float,
            )
            gap_residuals = np.asarray(
                [
                    rlb2e.centred_secant_residual(
                        gaps[index - 1], gaps[index], gaps[index + 1]
                    )
                    for index in range(1, len(grid) - 1)
                ],
                dtype=float,
            )
            median = float(np.median(gap_residuals))
            mad = float(np.median(np.abs(gap_residuals - median)))
            threshold = median + NEIGHBOUR_MAD_MULTIPLIER * mad
            for offset, index in enumerate(range(1, len(grid) - 1)):
                if (
                    float(gap_residuals[offset]) > threshold
                    and float(gap_residuals[offset]) > NEIGHBOUR_ABSOLUTE_TRIGGER
                ):
                    gap_flags.add((grid[index], lower_position))
                    gap_flags.add((grid[index], lower_position + 1))
        groups = _base_group_index(rows)
        for position in range(1, K_PLOT + 1):
            values = {
                value: float(spectra[value][position - 1]) for value in grid
            }
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
                parameter_value = grid[index]
                group = groups[(configuration_id, parameter_value)]
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
                unresolved_warning = (
                    int(root_row["unresolved_candidates_below_root9"]) != 0
                )
                bad_residual_warning = bool(
                    float(root_row["scaled_sigma_ratio"])
                    > ROOT_SINGULAR_RATIO_TOLERANCE
                    or float(root_row["boundary_null_residual"])
                    > BOUNDARY_RESIDUAL_TOLERANCE
                )
                gap_warning = (parameter_value, position) in gap_flags
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
                        "experiment_id": experiment_id,
                        "configuration_id": configuration_id,
                        "parameter_name": parameter_name,
                        "parameter_left": grid[index - 1],
                        "parameter_value": parameter_value,
                        "parameter_right": grid[index + 1],
                        "sorted_position": position,
                        "Lambda_left": values[grid[index - 1]],
                        "Lambda_center": values[parameter_value],
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
                        "smoothing_applied": False,
                    }
                )
    return result


def flagged_repair_points(
    audit_rows: Sequence[Mapping[str, Any]],
) -> list[tuple[str, str, float]]:
    return sorted(
        {
            (
                str(row["experiment_id"]),
                str(row["configuration_id"]),
                round(float(row["parameter_value"]), 10),
            )
            for row in audit_rows
            if _as_bool(row["flagged"])
        }
    )


def apply_local_repairs(
    rows_A: list[dict[str, Any]],
    rows_B: list[dict[str, Any]],
    audit_rows: list[dict[str, Any]],
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    records: list[dict[str, Any]] = []
    for repair_index, (experiment_id, configuration_id, parameter_value) in enumerate(
        flagged_repair_points(audit_rows), start=1
    ):
        repair_id = f"repair_{repair_index:04d}"
        target_rows = rows_A if experiment_id == EXPERIMENT_A else rows_B
        positions = sorted(
            {
                int(row["sorted_position"])
                for row in audit_rows
                if str(row["experiment_id"]) == experiment_id
                and str(row["configuration_id"]) == configuration_id
                and round(float(row["parameter_value"]), 10) == parameter_value
                and _as_bool(row["flagged"])
            }
        )
        original = _roots_for(target_rows, configuration_id, parameter_value)
        try:
            solution = solve_point(
                experiment_id,
                configuration_id,
                parameter_value,
                previous=(parameter_value, original),
                dense_local=True,
                dense_positions=positions,
                grid_kind="LOCAL_REFINEMENT",
                repair_id=repair_id,
                continuation_leg="LOCAL_REPAIR",
            )
        except (RuntimeError, ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            gaps: list[dict[str, Any]] = []
            for row in target_rows:
                if (
                    str(row["grid_kind"]) == "BASE"
                    and str(row["configuration_id"]) == configuration_id
                    and round(float(row["parameter_value"]), 10) == parameter_value
                    and int(row["sorted_position"]) in positions
                ):
                    row["is_canonical_plot_source"] = False
                    row["point_status"] = "UNRESOLVED_AFTER_LOCAL_REPAIR"
                    gap = dict(row)
                    gap["row_id"] = (
                        f"{experiment_id}__{configuration_id}__{parameter_value:+.6f}__"
                        f"LOCAL_REFINEMENT__p{int(row['sorted_position']):02d}__"
                        f"{repair_id}_gap"
                    )
                    gap["grid_kind"] = "LOCAL_REFINEMENT"
                    gap["Lambda"] = math.nan
                    gap["quality_status"] = "UNRESOLVED_AFTER_LOCAL_REPAIR"
                    gap["is_canonical_plot_source"] = True
                    gap["supersedes_row_id"] = row["row_id"]
                    gap["repair_id"] = repair_id
                    gaps.append(gap)
            target_rows.extend(gaps)
            status = "UNRESOLVED"
            records.append(
                {
                    "repair_id": repair_id,
                    "experiment_id": experiment_id,
                    "configuration_id": configuration_id,
                    "parameter_value": parameter_value,
                    "status": status,
                    "affected_positions": positions,
                    "error": f"{type(exc).__name__}: {exc}",
                    "wall_time_seconds": 0.0,
                    "determinant_evaluations": 0,
                    "sigma_evaluations": 0,
                    "smoothing_applied": False,
                    "predictor_used_as_final": False,
                }
            )
        else:
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
            original_by_position = {
                int(row["sorted_position"]): row
                for row in target_rows
                if str(row["grid_kind"]) == "BASE"
                and str(row["configuration_id"]) == configuration_id
                and round(float(row["parameter_value"]), 10) == parameter_value
            }
            for row in original_by_position.values():
                row["is_canonical_plot_source"] = False
                row["point_status"] = status
            for row in solution.rows:
                row["grid_kind"] = "LOCAL_REFINEMENT"
                row["repair_id"] = repair_id
                row["is_canonical_plot_source"] = True
                row["supersedes_row_id"] = original_by_position[
                    int(row["sorted_position"])
                ]["row_id"]
                row["point_status"] = status
            target_rows.extend(dict(row) for row in solution.rows)
            records.append(
                {
                    "repair_id": repair_id,
                    "experiment_id": experiment_id,
                    "configuration_id": configuration_id,
                    "parameter_value": parameter_value,
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
        for audit in audit_rows:
            if (
                str(audit["experiment_id"]) == experiment_id
                and str(audit["configuration_id"]) == configuration_id
                and round(float(audit["parameter_value"]), 10) == parameter_value
                and _as_bool(audit["flagged"])
            ):
                audit["repair_id"] = repair_id
                audit["repair_status"] = status
        if experiment_id == EXPERIMENT_A:
            rows_A = target_rows
        else:
            rows_B = target_rows
    _sort_rows(rows_A)
    _sort_rows(rows_B)
    return rows_A, rows_B, audit_rows, records


def canonical_plot_rows(
    rows: Sequence[Mapping[str, Any]],
) -> list[Mapping[str, Any]]:
    return [
        row
        for row in rows
        if int(row["sorted_position"]) <= K_PLOT
        and _as_bool(row.get("is_canonical_plot_source", True))
    ]


def audit_plot_rows(
    rows: Sequence[Mapping[str, Any]], experiment_id: str
) -> dict[str, Any]:
    selected = canonical_plot_rows(rows)
    expected = {
        (configuration_id, value, position)
        for configuration_id, value in _expected_keys(experiment_id)
        for position in range(1, K_PLOT + 1)
    }
    counts: dict[tuple[str, float, int], int] = {}
    invalid: list[str] = []
    for row in selected:
        key = (
            str(row["configuration_id"]),
            round(float(row["parameter_value"]), 10),
            int(row["sorted_position"]),
        )
        counts[key] = counts.get(key, 0) + 1
        value = float(row["Lambda"])
        gap = math.isnan(value) and str(row.get("point_status", "")) == (
            "UNRESOLVED_AFTER_LOCAL_REPAIR"
        )
        if not ((math.isfinite(value) and value > 0.0) or gap):
            invalid.append(f"{key[0]}:{key[1]:+.2f}:p{key[2]}")
    failures = [
        f"{key[0]}:{key[1]:+.2f}:p{key[2]}"
        for key in sorted(expected | set(counts))
        if counts.get(key, 0) != 1
    ]
    return {
        "status": "PASS" if not failures and not invalid else "FAIL",
        "row_count": len(selected),
        "key_failures": failures,
        "invalid_values": invalid,
    }


def _atomic_save_figure(figure: Any, path: Path, *, dpi: int = 190) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.stem + ".tmp" + target.suffix)
    try:
        figure.savefig(temporary, dpi=dpi, bbox_inches="tight")
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def create_plots_from_csv(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    """Render all three figures without physics, determinant, SVD or roots."""

    global ROOT_CALCULATION_COUNT
    before = ROOT_CALCULATION_COUNT
    target = Path(output_dir)
    rows_A = rlb2e._read_csv(target / SPECTRUM_A_FILENAME)
    rows_B = rlb2e._read_csv(target / SPECTRUM_B_FILENAME)
    section_rows = rlb2e._read_csv(target / SECTION_FILENAME)
    spectrum_A = audit_spectrum_rows(rows_A, EXPERIMENT_A)
    spectrum_B = audit_spectrum_rows(rows_B, EXPERIMENT_B)
    section_audit = audit_section_property_rows(section_rows)
    if (
        spectrum_A["status"] != "PASS"
        or spectrum_B["status"] != "PASS"
        or section_audit["status"] != "PASS"
    ):
        raise RuntimeError(
            "Plot-only requires complete first-eight-plus-root9 and section "
            "inventories: "
            f"A={spectrum_A}, B={spectrum_B}, sections={section_audit}"
        )
    audit_A = audit_plot_rows(rows_A, EXPERIMENT_A)
    audit_B = audit_plot_rows(rows_B, EXPERIMENT_B)
    if audit_A["status"] != "PASS" or audit_B["status"] != "PASS":
        raise RuntimeError(f"Plot-data audit failed: A={audit_A}, B={audit_B}")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    colours = [plt.get_cmap("tab10")(index) for index in range(K_PLOT)]
    figure_A, axes = plt.subplots(1, 3, figsize=(15.0, 4.8), sharey=True)
    base_A = canonical_plot_rows(rows_A)
    for axis, configuration_id, panel in zip(
        axes, CONFIGURATIONS_A, ("(a)", "(b)", "(c)"), strict=True
    ):
        for position in range(1, K_PLOT + 1):
            selected = sorted(
                (
                    row
                    for row in base_A
                    if str(row["configuration_id"]) == configuration_id
                    and int(row["sorted_position"]) == position
                ),
                key=lambda row: float(row["parameter_value"]),
            )
            axis.plot(
                [float(row["parameter_value"]) for row in selected],
                [float(row["Lambda"]) for row in selected],
                color=colours[position - 1],
                linewidth=1.35,
                linestyle="-",
            )
        axis.set_xlabel(r"$\eta$")
        axis.set_title(f"{panel} {configuration_id}", fontsize=10)
        axis.grid(alpha=0.22)
    axes[0].set_ylabel(r"$\Lambda$")
    handles = [
        Line2D([0], [0], color=colours[index], lw=1.6, label=f"position {index+1}")
        for index in range(K_PLOT)
    ]
    figure_A.legend(handles=handles, loc="upper center", ncol=8, frameon=False)
    figure_A.subplots_adjust(top=0.82, bottom=0.20)
    figure_A.text(
        0.5,
        0.02,
        "H: higher density; L: lower density; elastic properties are identical.",
        ha="center",
        fontsize=9,
    )
    _atomic_save_figure(figure_A, target / PLOT_A_FILENAME)
    plt.close(figure_A)

    figure_B, axis_B = plt.subplots(figsize=(8.8, 5.4))
    base_B = canonical_plot_rows(rows_B)
    for position in range(1, K_PLOT + 1):
        selected = sorted(
            (row for row in base_B if int(row["sorted_position"]) == position),
            key=lambda row: float(row["parameter_value"]),
        )
        axis_B.plot(
            [float(row["parameter_value"]) for row in selected],
            [float(row["Lambda"]) for row in selected],
            color=colours[position - 1],
            linewidth=1.35,
            linestyle="-",
            label=f"position {position}",
        )
    axis_B.axvline(0.0, color="0.35", linewidth=0.9, linestyle=":")
    axis_B.set_xlabel(r"$\xi_\rho$")
    axis_B.set_ylabel(r"$\Lambda$")
    axis_B.grid(alpha=0.22)
    axis_B.legend(ncol=2, frameon=False)
    figure_B.subplots_adjust(bottom=0.20)
    figure_B.text(
        0.5,
        0.02,
        r"$\xi_\rho<0$: heavier inner plies; $\xi_\rho>0$: heavier outer plies.",
        ha="center",
        fontsize=9,
    )
    _atomic_save_figure(figure_B, target / PLOT_B_FILENAME)
    plt.close(figure_B)

    def property_series(
        experiment_id: str, configuration_id: str, arm_id: int
    ) -> tuple[list[float], list[float]]:
        selected = sorted(
            (
                row
                for row in section_rows
                if str(row["experiment_id"]) == experiment_id
                and str(row["configuration_id"]) == configuration_id
                and int(row["arm_id"]) == arm_id
            ),
            key=lambda row: float(row["parameter_value"]),
        )
        return (
            [float(row["parameter_value"]) for row in selected],
            [float(row["J_over_J0"]) for row in selected],
        )

    eta_outer, J_outer = property_series(EXPERIMENT_A, CONFIG_BOTH_OUTER, 1)
    eta_inner, J_inner = property_series(EXPERIMENT_A, CONFIG_BOTH_INNER, 1)
    xi_values, J_layered = property_series(EXPERIMENT_B, CONFIG_ONE_ARM, 1)
    figure_J, (axis_eta, axis_xi) = plt.subplots(1, 2, figsize=(10.5, 4.4))
    axis_eta.plot(eta_outer, J_outer, label=r"$J_{outer}/J_0$", linewidth=1.6)
    axis_eta.plot(eta_inner, J_inner, label=r"$J_{inner}/J_0$", linewidth=1.6)
    axis_eta.set_xlabel(r"$\eta$")
    axis_eta.set_ylabel(r"$J/J_0$")
    axis_eta.grid(alpha=0.22)
    axis_eta.legend(frameon=False)
    axis_xi.plot(xi_values, J_layered, color="C2", label=r"$J_{layered}/J_0$")
    axis_xi.axvline(0.0, color="0.35", linewidth=0.9, linestyle=":")
    axis_xi.set_xlabel(r"$\xi_\rho$")
    axis_xi.set_ylabel(r"$J/J_0$")
    axis_xi.grid(alpha=0.22)
    axis_xi.legend(frameon=False)
    _atomic_save_figure(figure_J, target / PLOT_J_FILENAME)
    plt.close(figure_J)
    if ROOT_CALCULATION_COUNT != before:
        raise RuntimeError("plot-only unexpectedly calculated spectral roots.")
    return {
        "status": "PASS",
        "plot_A": (target / PLOT_A_FILENAME).as_posix(),
        "plot_B": (target / PLOT_B_FILENAME).as_posix(),
        "plot_J": (target / PLOT_J_FILENAME).as_posix(),
        "experiment_A_panels": 3,
        "experiment_B_panels": 1,
        "plotted_positions": list(range(1, K_PLOT + 1)),
        "root9_plotted": False,
        "frequency_line_count_A": 3 * K_PLOT,
        "frequency_line_count_B": K_PLOT,
        "property_line_count": 3,
        "root_calculation_count": ROOT_CALCULATION_COUNT - before,
        "matrix_assembly_count": 0,
        "determinant_evaluations": 0,
        "SVD_evaluations": 0,
        "smoothing_applied": False,
    }


def arm_swap_checks(
    rows_A: Sequence[Mapping[str, Any]], rows_B: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []
    declared = (
        (EXPERIMENT_A, CONFIG_ANTI_PHASE, 0.4),
        (EXPERIMENT_A, CONFIG_ANTI_PHASE, 0.8),
        (EXPERIMENT_B, CONFIG_ONE_ARM, -0.8),
        (EXPERIMENT_B, CONFIG_ONE_ARM, -0.4),
        (EXPERIMENT_B, CONFIG_ONE_ARM, 0.4),
        (EXPERIMENT_B, CONFIG_ONE_ARM, 0.8),
    )
    for experiment_id, configuration_id, parameter_value in declared:
        source = rows_A if experiment_id == EXPERIMENT_A else rows_B
        reference = _roots_for(source, configuration_id, parameter_value)
        solution = solve_point(
            experiment_id,
            configuration_id,
            parameter_value,
            previous=(parameter_value, reference),
            continuation_leg="ARM_SWAP_DIAGNOSTIC",
            swapped=True,
            guard_locator_right_width=0.2,
            grid_kind="ARM_SWAP_DIAGNOSTIC",
        )
        swapped = np.asarray([float(row["Omega"]) for row in solution.rows])
        relative = np.abs(reference - swapped) / np.maximum.reduce(
            (
                np.abs(reference),
                np.abs(swapped),
                np.full(K_GUARD, np.finfo(float).tiny),
            )
        )
        maximum = float(np.max(relative))
        checks.append(
            {
                "experiment_id": experiment_id,
                "configuration_id": configuration_id,
                "parameter_value": parameter_value,
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
        "diagnostic_configuration_exported_to_spectrum_csv": False,
    }


def _validated_predecessor_zero(
    result_dir: Path, spectrum_filename: str, selector: Any
) -> dict[str, Any]:
    manifest_path = Path(result_dir) / "run_manifest.json"
    spectrum_path = Path(result_dir) / spectrum_filename
    if not manifest_path.is_file() or not spectrum_path.is_file():
        return {"status": "NOT_AVAILABLE", "reason": "missing files"}
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        recorded = manifest.get("output_hashes", {}).get(spectrum_filename)
        actual = rlb2e._sha256(spectrum_path)
        if recorded != actual:
            return {
                "status": "NOT_AVAILABLE_HASH_MISMATCH",
                "path": spectrum_path.as_posix(),
                "recorded_sha256": recorded,
                "actual_sha256": actual,
            }
        rows = rlb2e._read_csv(spectrum_path)
        selected = selector(rows)
        selected.sort(key=lambda row: int(row["sorted_position"]))
        if [int(row["sorted_position"]) for row in selected] != list(
            range(1, K_GUARD + 1)
        ):
            return {"status": "NOT_AVAILABLE", "reason": "incomplete zero group"}
        return {
            "status": "AVAILABLE",
            "path": spectrum_path.as_posix(),
            "sha256": actual,
            "Omega": [float(row["Omega"]) for row in selected],
        }
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        return {"status": "NOT_AVAILABLE", "reason": f"{type(exc).__name__}: {exc}"}


def baseline_consistency(
    rows_A: Sequence[Mapping[str, Any]], rows_B: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    groups = [
        _roots_for(rows_A, configuration_id, 0.0)
        for configuration_id in CONFIGURATIONS_A
    ]
    groups.append(_roots_for(rows_B, CONFIG_ONE_ARM, 0.0))
    reference = groups[0]
    logical_max = max(
        float(
            np.max(
                np.abs(reference - group)
                / np.maximum.reduce(
                    (
                        np.abs(reference),
                        np.abs(group),
                        np.full(K_GUARD, np.finfo(float).tiny),
                    )
                )
            )
        )
        for group in groups[1:]
    )
    predecessor_A = _validated_predecessor_zero(
        RLB2E_RESULT_DIR,
        "spectrum_roots.csv",
        lambda rows: [
            row
            for row in rows
            if str(row.get("configuration_id")) == rlb2e.CONFIG_BOTH_OUTER
            and round(float(row.get("chi", math.nan)), 10) == 0.0
            and str(row.get("grid_kind")) == "BASE"
        ],
    )
    predecessor_B = _validated_predecessor_zero(
        RLB2F_RESULT_DIR,
        "spectrum_roots.csv",
        lambda rows: [
            row
            for row in rows
            if round(float(row.get("xi", math.nan)), 10) == 0.0
            and str(row.get("grid_kind")) == "BASE"
        ],
    )
    comparisons: list[dict[str, Any]] = []
    for name, evidence in (("RLB-2E", predecessor_A), ("RLB-2F", predecessor_B)):
        if evidence.get("status") != "AVAILABLE":
            comparisons.append({"source": name, **evidence})
            continue
        values = np.asarray(evidence["Omega"], dtype=float)
        relative = float(
            np.max(
                np.abs(reference - values)
                / np.maximum.reduce(
                    (
                        np.abs(reference),
                        np.abs(values),
                        np.full(K_GUARD, np.finfo(float).tiny),
                    )
                )
            )
        )
        comparisons.append(
            {
                "source": name,
                **evidence,
                "maximum_relative_Omega": relative,
                "comparison_status": (
                    "PASS" if relative <= ARM_SWAP_RELATIVE_TOLERANCE else "FAIL"
                ),
            }
        )
    predecessor_fail = any(
        item.get("status") == "AVAILABLE" and item.get("comparison_status") != "PASS"
        for item in comparisons
    )
    passed = logical_max <= ARM_SWAP_RELATIVE_TOLERANCE and not predecessor_fail
    return {
        "status": "PASS" if passed else "FAIL",
        "shared_logical_group_count": 4,
        "shared_physical_solve_count": 1,
        "maximum_logical_relative_Omega": logical_max,
        "tolerance": ARM_SWAP_RELATIVE_TOLERANCE,
        "predecessor_evidence_optional": True,
        "predecessor_comparisons": comparisons,
    }


def _endpoint_changes(
    rows_A: Sequence[Mapping[str, Any]], rows_B: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for experiment_id, configuration_id, left, right, rows in (
        *[
            (EXPERIMENT_A, config, 0.0, 0.8, rows_A)
            for config in CONFIGURATIONS_A
        ],
        (EXPERIMENT_B, CONFIG_ONE_ARM, -0.8, 0.8, rows_B),
    ):
        left_values = np.sqrt(_roots_for(rows, configuration_id, left)[:K_PLOT])
        right_values = np.sqrt(_roots_for(rows, configuration_id, right)[:K_PLOT])
        for position, (left_value, right_value) in enumerate(
            zip(left_values, right_values, strict=True), start=1
        ):
            result.append(
                {
                    "experiment_id": experiment_id,
                    "configuration_id": configuration_id,
                    "sorted_position": position,
                    "parameter_left": left,
                    "parameter_right": right,
                    "Lambda_left": left_value,
                    "Lambda_right": right_value,
                    "relative_endpoint_change": (
                        float(right_value - left_value)
                        / max(abs(float(left_value)), np.finfo(float).tiny)
                    ),
                }
            )
    return result


def _monotonicity(
    rows_A: Sequence[Mapping[str, Any]], rows_B: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for experiment_id, configuration_id, _name, grid, rows in _curve_specs(
        rows_A, rows_B
    ):
        for position in range(1, K_PLOT + 1):
            values = np.asarray(
                [
                    math.sqrt(
                        _roots_for(rows, configuration_id, value)[position - 1]
                    )
                    for value in grid
                ],
                dtype=float,
            )
            differences = np.diff(values)
            scale = max(float(np.max(np.abs(values))), 1.0)
            tolerance = 1.0e-10 * scale
            if np.all(differences >= -tolerance):
                classification = "NONDECREASING"
            elif np.all(differences <= tolerance):
                classification = "NONINCREASING"
            else:
                classification = "NONMONOTONE_ON_FINITE_GRID"
            result.append(
                {
                    "experiment_id": experiment_id,
                    "configuration_id": configuration_id,
                    "sorted_position": position,
                    "classification": classification,
                    "minimum_step_change": float(np.min(differences)),
                    "maximum_step_change": float(np.max(differences)),
                }
            )
    return result


def _minimum_adjacent_gaps(
    rows_A: Sequence[Mapping[str, Any]], rows_B: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for experiment_id, configuration_id, parameter_name, grid, rows in _curve_specs(
        rows_A, rows_B
    ):
        best: dict[str, Any] | None = None
        for value in grid:
            lambdas = np.sqrt(_roots_for(rows, configuration_id, value)[:K_PLOT])
            gaps = np.diff(lambdas)
            index = int(np.argmin(gaps))
            candidate = {
                "experiment_id": experiment_id,
                "configuration_id": configuration_id,
                "parameter_name": parameter_name,
                "parameter_value": value,
                "between_sorted_positions": [index + 1, index + 2],
                "minimum_adjacent_Lambda_gap": float(gaps[index]),
                "interpretation": "candidate_interval_only",
            }
            if best is None or candidate["minimum_adjacent_Lambda_gap"] < best[
                "minimum_adjacent_Lambda_gap"
            ]:
                best = candidate
        if best is not None:
            result.append(best)
    return result


def _root_quality_summary(
    rows_A: Sequence[Mapping[str, Any]], rows_B: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    base = [
        row
        for row in [*rows_A, *rows_B]
        if str(row.get("grid_kind")) == "BASE"
    ]
    return {
        "maximum_base_scaled_sigma_ratio": max(
            float(row["scaled_sigma_ratio"]) for row in base
        ),
        "maximum_base_boundary_null_residual": max(
            float(row["boundary_null_residual"]) for row in base
        ),
        "maximum_unresolved_candidates_below_root9": max(
            int(row["unresolved_candidates_below_root9"]) for row in base
        ),
        "minimum_root9_right_margin_Omega": min(
            float(row["root9_right_margin_Omega"])
            for row in base
            if int(row["sorted_position"]) == K_GUARD
        ),
        "roots_above_9_computed_count": sum(
            _as_bool(row["roots_above_9_computed"]) for row in base
        ),
    }


def _runtime_summary(
    point_records: Sequence[Mapping[str, Any]],
    repair_records: Sequence[Mapping[str, Any]],
    arm_swap: Mapping[str, Any],
    plot: Mapping[str, Any],
) -> dict[str, Any]:
    swap_records = arm_swap.get("checks", [])
    base_seconds = sum(float(item.get("wall_time_seconds", 0.0)) for item in point_records)
    repair_seconds = sum(
        float(item.get("wall_time_seconds", 0.0)) for item in repair_records
    )
    swap_seconds = sum(float(item.get("wall_time_seconds", 0.0)) for item in swap_records)
    determinant = sum(
        int(item.get("determinant_evaluations", 0))
        for item in [*point_records, *repair_records, *swap_records]
    )
    sigma = sum(
        int(item.get("sigma_evaluations", 0))
        for item in [*point_records, *repair_records, *swap_records]
    )
    peak = max(
        [rlb2e._peak_rss_bytes()]
        + [int(item.get("peak_rss_bytes", 0)) for item in point_records]
        + [int(item.get("peak_rss_bytes", 0)) for item in repair_records]
        + [int(item.get("peak_rss_bytes", 0)) for item in swap_records]
    )
    plot_seconds = float(plot.get("wall_time_seconds", 0.0))
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
        "unique_base_root_solve_count": sum(
            int(item.get("determinant_evaluations", 0)) > 0
            and not bool(item.get("swapped_arms", False))
            for item in point_records
        ),
        "local_repair_solve_count": len(repair_records),
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
    runtime = manifest["runtime"]
    audit = manifest["neighbour_audit"]
    endpoint_lines = "\n".join(
        "- {experiment_id}, {configuration_id}, position {sorted_position}: "
        "$\\Delta_{{end}}={relative_endpoint_change:.6g}$.".format(**item)
        for item in manifest["endpoint_changes"]
    )
    monotonic_lines = "\n".join(
        "- {experiment_id}, {configuration_id}, position {sorted_position}: "
        "`{classification}`.".format(**item)
        for item in manifest["sorted_position_monotonicity"]
    )
    gap_lines = "\n".join(
        "- {experiment_id}, {configuration_id}: $\\Delta\\Lambda_{{min}}="
        "{minimum_adjacent_Lambda_gap:.6g}$ at {parameter_name}="
        "{parameter_value:.2f} between positions {between_sorted_positions[0]} "
        "and {between_sorted_positions[1]}.".format(**item)
        for item in manifest["minimum_adjacent_sorted_gaps"]
    )
    endpoints_by_configuration: dict[str, dict[int, float]] = {}
    for item in manifest["endpoint_changes"]:
        endpoints_by_configuration.setdefault(str(item["configuration_id"]), {})[
            int(item["sorted_position"])
        ] = float(item["relative_endpoint_change"])

    def absolute_extrema(configuration_id: str) -> tuple[tuple[int, float], tuple[int, float]]:
        values = endpoints_by_configuration[configuration_id]
        ordered = list(values.items())
        return (
            min(ordered, key=lambda pair: abs(pair[1])),
            max(ordered, key=lambda pair: abs(pair[1])),
        )

    outer_min, outer_max = absolute_extrema(CONFIG_BOTH_OUTER)
    inner_min, inner_max = absolute_extrema(CONFIG_BOTH_INNER)
    anti_min, anti_max = absolute_extrema(CONFIG_ANTI_PHASE)
    one_arm_min, one_arm_max = absolute_extrema(CONFIG_ONE_ARM)
    symmetric_maximum = max(
        abs(value)
        for configuration_id in (CONFIG_BOTH_OUTER, CONFIG_BOTH_INNER)
        for value in endpoints_by_configuration[configuration_id].values()
    )
    anti_phase_reduction = symmetric_maximum / abs(anti_max[1])
    one_arm_outer_difference = max(
        abs(
            endpoints_by_configuration[CONFIG_ONE_ARM][position]
            - endpoints_by_configuration[CONFIG_BOTH_OUTER][position]
        )
        for position in range(1, K_PLOT + 1)
    )

    def latex_scientific(value: float) -> str:
        exponent = int(math.floor(math.log10(abs(value))))
        mantissa = value / (10.0**exponent)
        return f"{mantissa:.6g}\\cdot10^{{{exponent}}}"

    sensitivity_summary = (
        "Для `BOTH_OUTER_HEAVY` и `BOTH_INNER_HEAVY` наименьший по модулю endpoint\n"
        f"shift получен в position {outer_min[0]} "
        f"($|\\Delta_{{end}}|={abs(outer_min[1]):.6g}$ и "
        f"${abs(inner_min[1]):.6g}$),\n"
        f"а наибольший — в position {outer_max[0]} "
        f"($|\\Delta_{{end}}|={abs(outer_max[1]):.6g}$ и "
        f"${abs(inner_max[1]):.6g}$).\n"
        "Position 7 остаётся близкой по чувствительности к position 2, поэтому\n"
        "чувствительность не растёт монотонно с номером независимо упорядоченной\n"
        "позиции.\n\n"
        "Для `ANTI_PHASE` минимальный и максимальный по модулю сдвиги относятся к\n"
        f"positions {anti_min[0]} и {anti_max[0]}; максимальный модуль равен "
        f"${latex_scientific(abs(anti_max[1]))}$, что в\n"
        f"{anti_phase_reduction:.1f} раза меньше максимального симметричного "
        "сдвига. Постоянство\n"
        "$J_1+J_2$ существенно уменьшает, но не устраняет изменения упорядоченных\n"
        "частот.\n\n"
        f"В one-arm experiment наименее и наиболее чувствительны positions "
        f"{one_arm_min[0]} и {one_arm_max[0]}.\n"
        "Максимальная разность его endpoint shifts и соответствующих shifts\n"
        "`BOTH_OUTER_HEAVY` равна "
        f"${latex_scientific(one_arm_outer_difference)}$. Эта близость относится\n"
        "только к данному конечному сопоставлению и не устанавливает эквивалентность\n"
        "конфигураций или идентичность модальных ветвей."
    )
    return rf"""# RLB-2G: двойственность компоновок по массе

## Цель и область расчёта

Рассматриваются четырёхслойные сопряжённые стержни с разным распределением
массы по толщине. Сопоставляются первые восемь независимо упорядоченных частот. Все
слои имеют одинаковые упругие свойства. Изменяется только плотность слоёв при
неизменной погонной массе каждого плеча. Термин mass-layout duality обозначает
два density-аналога экспериментов RLB-2E и RLB-2F, а не общую математическую
двойственность жёсткости и массы.

## Фиксированный контракт

Использованы $\mu=\tau=0$, $\beta=30^\circ$, $L_1=L_2=l=1$,
$b=0.20$, $h=0.05$ и $K=5/6$. Каждое плечо содержит четыре равных слоя
толщины 0.0125. Ориентация материала во всех слоях равна $0^\circ$.

Базовый материал имеет $E_1=1.1$, $E_2=0.9$, $\nu_{{12}}=0.3$,
$G_{{12}}=G_{{13}}=G_{{23}}=1/2.6$ и $\rho_0=1$. Для положительного
контраста $\rho_H=(1+\eta)\rho_0$ и $\rho_L=(1-\eta)\rho_0$.
В знаковом эксперименте плотности наружных и внутренних слоёв первого плеча
равны соответственно $(1+\xi_\rho)\rho_0$ и
$(1-\xi_\rho)\rho_0$. Второе плечо остаётся однородным.

## Constitutive gate

Статус: **{gate['status']}**. Получены

$$A_{{beam}}={gate['Abeam0']:.16g},\quad
D_{{beam}}={gate['Dbeam0']:.16g},\quad
S_{{beam}}={gate['Sbeam0']:.16g},$$

$$m={gate['m0']:.16g},\qquad J_0={gate['J0']:.16g}.$$

Для всех объявленных компоновок $A$, $D$, $S$ и $m$ сохраняются. Изменяется
только $J$. Максимальная невязка аналитической формулы для $J/J_0$ равна
{gate['maximum_residuals']['J_formula_residual']:.3e}. Проверены соотношения

$$\frac{{J_{{outer}}}}{{J_0}}=1+\frac{{3\eta}}{{4}},\qquad
\frac{{J_{{inner}}}}{{J_0}}=1-\frac{{3\eta}}{{4}},\qquad
\frac{{J_{{layered}}}}{{J_0}}=1+\frac{{3\xi_\rho}}{{4}}.$$

Поэтому внутри текущей одномерной модели RLB спектральные различия
рассмотренных компоновок связаны только с изменением вращательной инерции.

## Frequency-map policy и численная сборка

Оба локальных instance используют `frequency-map-v1`, режим `fast_plot` и
семантику `sorted_positions`. Эксперимент A содержит 41 значение
$\eta=0,0.02,\ldots,0.8$ для трёх конфигураций. Эксперимент B содержит 81
значение $\xi_\rho=-0.8,-0.78,\ldots,0.8$ и две continuation legs от нуля.

Позиции 1--8 показаны на рисунках. Root 9 служит только проверкой полноты.
Корни с позициями 10 и выше не вычислялись. Predictor использовался только
для локализации, а каждое сохранённое значение получено из замороженной
characteristic matrix с существующими determinant/SVD detector и refiner.
Нормировка не зависит от $J$:

$$\Omega=\omega l^2\sqrt{{\rho_0A_0/(E_0I_{{y0}})}},\qquad
\Lambda=\sqrt{{\Omega}}.$$

## Результаты и диагностика

Эксперимент A содержит {counts['experiment_A_groups']}/123 групп и
{counts['experiment_A_base_rows']}/1107 строк `BASE`. Эксперимент B содержит
{counts['experiment_B_groups']}/81 групп и
{counts['experiment_B_base_rows']}/729 строк `BASE`. Общий однородный спектр
рассчитан один раз и использован для четырёх логических групп.

Production benchmark дал консервативную ETA
{manifest['benchmark']['conservative_eta_seconds']:.1f} s при лимите 2700 s.
Сумма измеренных времён workflow равна
{runtime['total_measured_workflow_seconds']:.1f} s, peak RSS --
{runtime['peak_rss_bytes'] / 2**20:.1f} MiB. Выполнено
{runtime['determinant_evaluations']} determinant и
{runtime['sigma_evaluations']} SVD/sigma evaluations.

Neighbour audit отметил {audit['flagged_point_count']} точек. Выполнено
{audit['repair_count']} локальных проверок; unresolved points:
{audit['unresolved_point_count']}. Сглаживание не применялось. Проверка
перестановки плеч имеет статус **{manifest['arm_swap']['status']}**.

Относительные изменения между концами объявленных сеток:

{endpoint_lines}

{sensitivity_summary}

Классификация упорядоченных позиций на конечных сетках:

{monotonic_lines}

Минимальные соседние интервалы:

{gap_lines}

Эти интервалы являются только кандидатами для отдельного анализа форм.
По частотным картам нельзя установить crossing, veering, обмен модальным
характером или локализацию.

## Статус и ограничения

**RLB-2G: {manifest['scientific_status']}.** Результат относится только к
фиксированной геометрии, материалу $M_0$, двум объявленным схемам изменения
плотности и конечным сеткам параметров. Branch tracking, MAC, формы, анализ
энергии, Ritz, FEM, damping, complex roots и certified audit не выполнялись.
Production physics не изменялась.
"""


def finalize_outputs(
    output_dir: Path,
    rows_A: list[dict[str, Any]],
    rows_B: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    benchmark: Mapping[str, Any],
    point_records: Sequence[Mapping[str, Any]],
    started_at: str,
) -> dict[str, Any]:
    target = Path(output_dir)
    audit_rows = neighbour_audit_rows(rows_A, rows_B)
    rows_A, rows_B, audit_rows, repair_records = apply_local_repairs(
        rows_A, rows_B, audit_rows
    )
    _write_spectrum(target, EXPERIMENT_A, rows_A)
    _write_spectrum(target, EXPERIMENT_B, rows_B)
    rlb2e._atomic_write_csv(target / AUDIT_FILENAME, audit_rows)
    audit_A = audit_spectrum_rows(rows_A, EXPERIMENT_A)
    audit_B = audit_spectrum_rows(rows_B, EXPERIMENT_B)
    if audit_A["status"] != "PASS" or audit_B["status"] != "PASS":
        raise RuntimeError(f"Final spectrum audit failed: A={audit_A}, B={audit_B}")
    arm_swap = arm_swap_checks(rows_A, rows_B)
    rlb2e._atomic_write_json(target / ARM_SWAP_FILENAME, arm_swap)
    if arm_swap["status"] != "PASS":
        raise RuntimeError(f"Arm-swap diagnostic failed: {arm_swap}")
    baseline = baseline_consistency(rows_A, rows_B)
    if baseline["status"] != "PASS":
        raise RuntimeError(f"Shared-zero baseline consistency failed: {baseline}")
    plot_started = time.perf_counter()
    plot = create_plots_from_csv(target)
    plot["wall_time_seconds"] = time.perf_counter() - plot_started
    unresolved = sum(record["status"] == "UNRESOLVED" for record in repair_records)
    scientific_status = "PASS" if unresolved == 0 else "PARTIAL_PASS"
    counts = {
        "experiment_A_groups": audit_A["base_group_count"],
        "experiment_A_base_rows": audit_A["base_row_count"],
        "experiment_A_root9_guards": sum(
            str(row["grid_kind"]) == "BASE"
            and int(row["sorted_position"]) == K_GUARD
            for row in rows_A
        ),
        "experiment_B_groups": audit_B["base_group_count"],
        "experiment_B_base_rows": audit_B["base_row_count"],
        "experiment_B_root9_guards": sum(
            str(row["grid_kind"]) == "BASE"
            and int(row["sorted_position"]) == K_GUARD
            for row in rows_B
        ),
        "logical_group_count": audit_A["base_group_count"]
        + audit_B["base_group_count"],
        "logical_base_row_count": audit_A["base_row_count"]
        + audit_B["base_row_count"],
        "unique_physical_base_points": 201,
        "shared_zero_logical_groups": 4,
        "shared_zero_physical_solves": 1,
        "shared_zero_clones": 3,
        "section_property_rows": 408,
        "local_refinement_rows": sum(
            str(row["grid_kind"]) == "LOCAL_REFINEMENT"
            for row in [*rows_A, *rows_B]
        ),
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
        "reused_RLB2F_script_sha256": rlb2e._sha256(ROOT / RLB2F_SCRIPT_PATH),
        "constitutive_gate": dict(constitutive),
        "counts": counts,
        "spectrum_audits": {"experiment_A": audit_A, "experiment_B": audit_B},
        "root_quality_summary": _root_quality_summary(rows_A, rows_B),
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
        "endpoint_changes": _endpoint_changes(rows_A, rows_B),
        "sorted_position_monotonicity": _monotonicity(rows_A, rows_B),
        "minimum_adjacent_sorted_gaps": _minimum_adjacent_gaps(rows_A, rows_B),
        "plots": plot,
        "runtime": runtime,
        "root_contract": {
            "plotted_positions": list(range(1, K_PLOT + 1)),
            "guard_position": K_GUARD,
            "root9_role": "completeness_only",
            "root9_plotted": False,
            "roots_above_9_computed": False,
            "accepted_candidates_above_root9": 0,
            "branch_tracking": False,
        },
        "exclusions_confirmed": {
            "universal_spectral_runner_used": False,
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
        rows_A,
        rows_B,
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
        {
            "experiment_id": record["experiment_id"],
            "configuration_id": record["configuration_id"],
            "parameter_value": record["parameter_value"],
        }
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
    if not manifest_path.is_file():
        return None
    if not all((target / name).is_file() for name in MANDATORY_COMPLETED_OUTPUTS):
        return None
    rows_A = rlb2e._read_csv(target / SPECTRUM_A_FILENAME)
    rows_B = rlb2e._read_csv(target / SPECTRUM_B_FILENAME)
    audit_A = audit_spectrum_rows(rows_A, EXPERIMENT_A)
    audit_B = audit_spectrum_rows(rows_B, EXPERIMENT_B)
    if audit_A["status"] != "PASS" or audit_B["status"] != "PASS":
        raise RuntimeError(
            "A materially complete RLB-2G result tree has an invalid spectrum "
            f"inventory and is not resumable: A={audit_A}, B={audit_B}"
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    recorded_hashes = manifest.get("output_hashes", {})
    expected_physics = {
        path: rlb2e._sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS
    }
    mismatches: list[str] = []
    expected_scalars = {
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "contract_sha256": contract_hash(),
        "analysis_script_sha256": rlb2e._sha256(Path(__file__)),
        "reused_RLB2E_script_sha256": rlb2e._sha256(ROOT / RLB2E_SCRIPT_PATH),
        "reused_RLB2F_script_sha256": rlb2e._sha256(ROOT / RLB2F_SCRIPT_PATH),
    }
    for field, expected in expected_scalars.items():
        if manifest.get(field) != expected:
            mismatches.append(field)
    if manifest.get("production_physics_hashes") != expected_physics:
        mismatches.append("production_physics_hashes")
    if not isinstance(recorded_hashes, dict):
        mismatches.append("output_hashes:not-a-mapping")
    else:
        missing_hashes = sorted(MANDATORY_COMPLETED_OUTPUTS - recorded_hashes.keys())
        if missing_hashes:
            mismatches.append("output_hashes:missing=" + ",".join(missing_hashes))
        for name, expected in recorded_hashes.items():
            path = target / name
            if not path.is_file() or rlb2e._sha256(path) != expected:
                mismatches.append(f"output_hash:{name}")
    if mismatches:
        raise RuntimeError(
            "A materially complete RLB-2G result tree is hash-incompatible and "
            "must not be mutated or recomputed: " + "; ".join(mismatches)
        )
    return manifest


def run_workflow(
    output_dir: Path = DEFAULT_OUTPUT_DIR, *, missing_only: bool = False
) -> dict[str, Any]:
    """Run both finite RLB-2G maps or resume incomplete transactions."""

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
    rlb2e._atomic_write_csv(target / SECTION_FILENAME, section_rows, SECTION_FIELDS)
    if constitutive["status"] != "PASS":
        raise RuntimeError(f"RLB-2G constitutive gate failed: {constitutive}")
    failed_sections = [
        (
            f"{row['experiment_id']}:{row['configuration_id']}:"
            f"{float(row['parameter_value']):+.2f}:arm{int(row['arm_id'])}"
        )
        for row in section_rows
        if str(row["constitutive_status"]) != "PASS"
    ]
    if failed_sections:
        raise RuntimeError(
            "RLB-2G section-row constitutive gate failed: "
            + ", ".join(failed_sections)
        )
    rows_A = _recover_interrupted_finalization(
        target, EXPERIMENT_A, _read_spectrum(target, EXPERIMENT_A)
    )
    rows_B = _recover_interrupted_finalization(
        target, EXPERIMENT_B, _read_spectrum(target, EXPERIMENT_B)
    )
    point_records = _existing_point_records(target)
    _record_recovered_groups(rows_A, EXPERIMENT_A, point_records)
    _record_recovered_groups(rows_B, EXPERIMENT_B, point_records)
    rows_A, rows_B, benchmark = run_benchmarks(
        target, rows_A, rows_B, point_records, constitutive, started_at
    )
    if not benchmark["production_run_permitted"]:
        raise RuntimeError(
            "RLB-2G ETA exceeds 45 minutes: "
            f"{benchmark['conservative_eta_seconds']:.1f} s."
        )
    rows_A, rows_B = complete_missing_points(
        target, rows_A, rows_B, point_records, constitutive, started_at
    )
    return finalize_outputs(
        target,
        rows_A,
        rows_B,
        constitutive,
        benchmark,
        point_records,
        started_at,
    )


def manifest_only(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
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
        result = create_plots_from_csv(args.output_dir)
    else:
        result = run_workflow(args.output_dir, missing_only=args.missing_only)
    print(json.dumps(rlb2e._json_value(result), indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
