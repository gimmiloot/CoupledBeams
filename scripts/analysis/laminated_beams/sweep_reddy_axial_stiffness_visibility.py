"""RLB-2H axial-stiffness visibility under geometric coupling.

The workflow reuses the frozen RLB-2E/RLB-2F determinant/SVD search
machinery, but changes only the physical case: all four plies stay at
``angle_deg = 0`` and only the in-plane material stiffness distribution
between the outer and inner plies varies.  The reduction is designed so that
only the reduced axial stiffness ``A`` changes, while ``D``, ``S``, ``m`` and
``J`` remain fixed.

Two finite sorted-position maps are produced at ``beta = 0 deg`` and
``beta = 30 deg``.  Positions 1--8 are plotted and position 9 is retained
only as a completeness guard.  Predictors are numerical locators only; every
exported root is refined from the frozen characteristic matrix.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time
from types import SimpleNamespace
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
    pilot_reddy_symmetric_coupled_beams_beta0 as beta0_tools,
)
from scripts.analysis.laminated_beams import (  # noqa: E402
    sweep_reddy_stiffness_layout_contrast as rlb2e,
)


FloatArray = NDArray[np.float64]

STAGE_ID = "RLB-2H"
ALGORITHM_VERSION = "axial_stiffness_visibility_fast_plot_v1"
POLICY_ID = "frequency-map-v1"
REQUIRED_STATUS_GATE_NAMES = {
    "RLB-2H-CONSTITUTIVE-A-ONLY",
    "RLB-2H-BETA0-DIRECT-SUBSYSTEM-REFERENCE",
    "RLB-2H-FREQUENCY-MAP",
    "RLB-2H-ROOT-INVENTORY",
    "OVERALL",
}

BETA_VALUES_DEG = (0.0, 30.0)
BETA_REFERENCE = 30.0
MU = 0.0
TAU = 0.0
L_REFERENCE = 1.0
L1 = 1.0
L2 = 1.0
L_TOTAL = 2.0
WIDTH = 0.20
THICKNESS = 0.05
PLY_THICKNESS = THICKNESS / 4.0
K = 5.0 / 6.0

E1_0 = rlb2e.E1_0
E2_0 = rlb2e.E2_0
NU12_0 = rlb2e.NU12_0
G12_0 = rlb2e.G0
G13_0 = rlb2e.G0
G23_0 = rlb2e.G0
RHO_0 = rlb2e.RHO_0

ALPHA_A_MIN = 0.70
ALPHA_A_MAX = 1.30
ALPHA_A_STEP = 0.02
ALPHA_A_ANCHOR = 1.00

K_PLOT = rlb2e.K_PLOT
K_GUARD = rlb2e.K_GUARD
LOCAL_ALPHA_STEP = rlb2e.LOCAL_CHI_STEP

MATRIX_RELATIVE_TOLERANCE = rlb2e.MATRIX_RELATIVE_TOLERANCE
SYMMETRY_RELATIVE_TOLERANCE = rlb2e.SYMMETRY_RELATIVE_TOLERANCE
REDUCED_PROPERTY_TOLERANCE = rlb2e.REDUCED_PROPERTY_TOLERANCE
REDUCTION_ROUTE_TOLERANCE = rlb2e.REDUCTION_ROUTE_TOLERANCE
ROOT_SINGULAR_RATIO_TOLERANCE = rlb2e.ROOT_SINGULAR_RATIO_TOLERANCE
BOUNDARY_RESIDUAL_TOLERANCE = rlb2e.BOUNDARY_RESIDUAL_TOLERANCE
ROOT9_RIGHT_TAIL_OMEGA = rlb2e.ROOT9_RIGHT_TAIL_OMEGA
ETA_LIMIT_SECONDS = 45.0 * 60.0
NEIGHBOUR_MAD_MULTIPLIER = rlb2e.NEIGHBOUR_MAD_MULTIPLIER
NEIGHBOUR_ABSOLUTE_TRIGGER = rlb2e.NEIGHBOUR_ABSOLUTE_TRIGGER
DIRECT_REFERENCE_RELATIVE_TOLERANCE = (
    beta0_tools.THRESHOLDS["beta0_isolated_spectral_relative"]
)
CLUSTER_REFERENCE_RELATIVE_TOLERANCE = (
    beta0_tools.THRESHOLDS["cluster_center_relative"]
)

REFERENCE_AREA = WIDTH * THICKNESS
IY_REFERENCE = WIDTH * THICKNESS**3 / 12.0
OMEGA_TO_OMEGA_SCALE = (
    L_REFERENCE**2 * math.sqrt(RHO_0 * REFERENCE_AREA / (1.0 * IY_REFERENCE))
)

OUTER = "OUTER"
INNER = "INNER"
LAYOUT = (OUTER, INNER, INNER, OUTER)

DEFAULT_OUTPUT_DIR = (
    ROOT / "results" / "laminated_beams" / "reddy_axial_stiffness_visibility"
)
SPECTRUM_FILENAME = "spectrum_roots.csv"
SECTION_FILENAME = "section_properties.csv"
AUDIT_FILENAME = "neighbour_audit.csv"
REFERENCE_FILENAME = "beta0_subsystem_reference.csv"
BENCHMARK_FILENAME = "benchmark.json"
CHECKPOINT_FILENAME = "checkpoint.json"
PLOT_FILENAME = "lambda_vs_A_ratio_beta0_beta30.png"
BETA0_PLOT_FILENAME = "beta0_decoupled_reference.png"
REPORT_FILENAME = "report.md"
MANIFEST_FILENAME = "run_manifest.json"

PRODUCTION_PHYSICS_PATHS = rlb2e.PRODUCTION_PHYSICS_PATHS
EXPECTED_PRODUCTION_PHYSICS_HASHES = {
    "scripts/lib/reddy_symmetric_laminated_beam.py": (
        "9E3F94747FA3723D0FEE350562F29A0DB070C3E3A17DDCCA3795F1E69AEDBE4B"
    ),
    "scripts/lib/reddy_symmetric_coupled_beams.py": (
        "E70F7AF5B4BB61AA90525664E6C4834EF5A003F34B23D6C2741583D38DAAD9A7"
    ),
    "scripts/lib/reddy_inplane_geometry.py": (
        "C46A42C462264BC27C99C358AABD7DF49F94F928A60D8150FD320D8DFB37E99E"
    ),
}
PREDECESSOR_RESULT_DIRS = {
    "RLB-2A": ROOT / "results" / "laminated_beams" / "reddy_cross_ply_layer_order_pilot",
    "RLB-2D": ROOT
    / "results"
    / "laminated_beams"
    / "rectangular_weakly_orthotropic_graphs_mu_beta",
    "RLB-2E": ROOT / "results" / "laminated_beams" / "reddy_stiffness_layout_contrast_sweep",
    "RLB-2F": ROOT / "results" / "laminated_beams" / "reddy_one_arm_layered_contrast_sweep",
    "RLB-2G": ROOT / "results" / "laminated_beams" / "reddy_mass_layout_duality",
}
INITIAL_PREDECESSOR_TREE_HASHES = {
    "RLB-2A": "07D0B115FE8B42AC4EF11A32B3B37A075A5477203869BDB4A343F46C79533106",
    "RLB-2D": "86D34750EB13CE6039D8FFA18D9FE15A4CC518FCD7921A5646F3ADB0129F0250",
    "RLB-2E": "57E9FFCFD518FADF198F30C84F04EE181F1C645814A7BCFA834FCC920426B008",
    "RLB-2F": "10C80E8136AF917BCC5EFCB351FAD2FBE6665856A91C1F71241BC650372046C5",
    "RLB-2G": "4DA662EB77240C59B78017CCCD38136522561F2F3BE48BAD2FA50AACBB059CC1",
}

FREQUENCY_MAP_POLICY_BASE = {
    "frequency_map_policy": POLICY_ID,
    "calculation_mode": "fast_plot",
    "spectrum_semantics": "sorted_positions",
    "sweep_parameter": "alpha_A",
    "parameter_grid": "0.70:0.02:1.30",
    "continuation_anchor": "1.00",
    "continuation_paths": ["1.00:0.98:0.70", "1.00:1.02:1.30"],
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
    "beta_deg",
    "alpha_A",
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
    "cluster_id",
    "cluster_semantics",
    "cluster_multiplicity",
    "cluster_total_nullity",
    "cluster_center_Omega",
    "cluster_metadata_source",
    "unresolved_candidates_below_root9",
    "search_right_Omega",
    "root9_right_margin_Omega",
    "guard_cluster_multiplicity",
    "guard_cluster_extends_beyond_export",
    "solve_mode",
    "fallback_used",
    "quality_status",
    "point_status",
    "is_canonical_plot_source",
    "supersedes_row_id",
    "repair_id",
    "decoupled_reference_classification",
    "roots_above_9_computed",
)

SECTION_FIELDS = (
    "alpha_A",
    "s_outer",
    "s_inner",
    "stack_bottom_to_top",
    "properties_identical_between_arms",
    "A_beam",
    "A_beam_over_A0",
    "expected_A_beam_over_A0",
    "D_beam",
    "D_beam_over_D0",
    "S_beam",
    "S_beam_over_S0",
    "m",
    "m_over_m0",
    "J",
    "J_over_J0",
    "B_relative",
    "I1_relative",
    "A_matrix_scale_residual",
    "D_matrix_invariance_residual",
    "A_formula_residual",
    "D_formula_residual",
    "only_A_field_max_invariance_residual",
    "reduction_route_max_relative",
    "z_interfaces",
    "A_matrix",
    "B_matrix",
    "D_matrix",
    "shear_matrix_yz_xz",
    "constitutive_status",
)


@dataclass(frozen=True)
class SectionObjects:
    layout: tuple[str, str, str, str]
    laminate: Any
    properties: Any


@dataclass(frozen=True)
class PointSolution:
    beta_deg: float
    alpha_A: float
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
    candidate_count_total: int = K_GUARD
    accepted_candidates_above_root9: int = 0
    retained_slots_above_root9: int = 0
    roots_above_9_computed: bool = False


ROOT_CALCULATION_COUNT = 0


def alpha_grid() -> FloatArray:
    point_count = int(round((ALPHA_A_MAX - ALPHA_A_MIN) / ALPHA_A_STEP)) + 1
    return np.linspace(ALPHA_A_MIN, ALPHA_A_MAX, point_count, dtype=float)


def continuation_paths() -> tuple[FloatArray, FloatArray]:
    lower = np.asarray(
        sweep_grid_policy.parameter_grid(0.0, ALPHA_A_ANCHOR - ALPHA_A_MIN, ALPHA_A_STEP),
        dtype=float,
    )
    lower = ALPHA_A_ANCHOR - lower
    lower[0] = ALPHA_A_ANCHOR
    upper = np.asarray(
        sweep_grid_policy.parameter_grid(ALPHA_A_ANCHOR, ALPHA_A_MAX, ALPHA_A_STEP),
        dtype=float,
    )
    return lower, upper


def omega_to_Omega(omega: float) -> float:
    return float(omega) * OMEGA_TO_OMEGA_SCALE


def Omega_to_Lambda(Omega: float) -> float:
    value = float(Omega)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError("Omega must be finite and nonnegative.")
    return math.sqrt(value)


def base_material_contract() -> dict[str, float]:
    return {
        "E1_0": E1_0,
        "E2_0": E2_0,
        "nu12_0": NU12_0,
        "G12_0": G12_0,
        "G13_0": G13_0,
        "G23_0": G23_0,
        "rho_0": RHO_0,
    }


def _physics_modules() -> tuple[Any, Any]:
    return rlb2e._physics_modules()


def alpha_scales(alpha_A: float) -> tuple[float, float]:
    value = float(alpha_A)
    if not math.isfinite(value) or not ALPHA_A_MIN <= value <= ALPHA_A_MAX:
        raise ValueError("alpha_A must lie in [0.70, 1.30].")
    s_outer = (4.0 - value) / 3.0
    s_inner = (7.0 * value - 4.0) / 3.0
    if s_outer <= 0.0 or s_inner <= 0.0:
        raise ValueError("Outer and inner multipliers must remain positive.")
    return s_outer, s_inner


def build_scaled_material(scale: float, label: str) -> Any:
    beam, _coupled = _physics_modules()
    factor = float(scale)
    return beam.OrthotropicLamina(
        E1=factor * E1_0,
        E2=factor * E2_0,
        nu12=NU12_0,
        G12=factor * G12_0,
        G13=G13_0,
        G23=G23_0,
        rho=RHO_0,
        name=f"RLB-2H {label}(scale={factor:.6f})",
    )


def _baseline_section() -> SectionObjects:
    beam, _coupled = _physics_modules()
    material = build_scaled_material(1.0, "M0")
    laminate = beam.integrate_laminate(
        [beam.Ply(material, 0.0, PLY_THICKNESS, label="M0") for _ in range(4)]
    )
    properties = beam.reduce_to_beam_properties(
        laminate,
        width=WIDTH,
        K=K,
        symmetry_tolerance=SYMMETRY_RELATIVE_TOLERANCE,
        reduction_tolerance=REDUCTION_ROUTE_TOLERANCE,
    )
    return SectionObjects(("M0", "M0", "M0", "M0"), laminate, properties)


def build_layered_section(alpha_A: float) -> SectionObjects:
    beam, _coupled = _physics_modules()
    s_outer, s_inner = alpha_scales(alpha_A)
    materials = {
        OUTER: build_scaled_material(s_outer, OUTER),
        INNER: build_scaled_material(s_inner, INNER),
    }
    laminate = beam.integrate_laminate(
        [
            beam.Ply(
                materials[label],
                angle_deg=0.0,
                thickness=PLY_THICKNESS,
                label=label,
            )
            for label in LAYOUT
        ]
    )
    properties = beam.reduce_to_beam_properties(
        laminate,
        width=WIDTH,
        K=K,
        symmetry_tolerance=SYMMETRY_RELATIVE_TOLERANCE,
        reduction_tolerance=REDUCTION_ROUTE_TOLERANCE,
    )
    return SectionObjects(LAYOUT, laminate, properties)


def _relative(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / max(
        abs(float(left)), abs(float(right)), np.finfo(float).tiny
    )


def _matrix_relative(left: Any, right: Any) -> float:
    return rlb2e._matrix_relative(left, right)


def _scaled_B(laminate: Any) -> float:
    return rlb2e._scaled_B(laminate)


def _scaled_I1(laminate: Any) -> float:
    return rlb2e._scaled_I1(laminate)


def _reduction_max_residual(properties: Any) -> float:
    return rlb2e._reduction_max_residual(properties)


def constitutive_gate() -> dict[str, Any]:
    baseline = _baseline_section()
    beam, _coupled = _physics_modules()
    Q0 = beam.lamina_reduced_stiffness(build_scaled_material(1.0, "Q0"))
    shear0 = beam.lamina_transverse_shear_stiffness(build_scaled_material(1.0, "Q0"))
    checks: list[dict[str, Any]] = []
    maxima = {
        "A_matrix_scale_residual": 0.0,
        "D_matrix_invariance_residual": 0.0,
        "A_formula_residual": 0.0,
        "D_formula_residual": 0.0,
        "S_formula_residual": 0.0,
        "m_formula_residual": 0.0,
        "J_formula_residual": 0.0,
        "B_relative": 0.0,
        "I1_relative": 0.0,
        "only_A_field_max_invariance_residual": 0.0,
        "reduction_route_max_relative": 0.0,
        "Q_scaling_residual": 0.0,
        "shear_Q_invariance_residual": 0.0,
    }
    passed = True
    for alpha_A_value in alpha_grid():
        alpha_A = float(alpha_A_value)
        section = build_layered_section(alpha_A)
        s_outer, s_inner = alpha_scales(alpha_A)
        outer_material = build_scaled_material(s_outer, "OUTER")
        inner_material = build_scaled_material(s_inner, "INNER")
        Q_outer = beam.lamina_reduced_stiffness(outer_material)
        Q_inner = beam.lamina_reduced_stiffness(inner_material)
        shear_outer = beam.lamina_transverse_shear_stiffness(outer_material)
        shear_inner = beam.lamina_transverse_shear_stiffness(inner_material)
        A_matrix_scale_residual = _matrix_relative(
            section.laminate.A, alpha_A * baseline.laminate.A
        )
        D_matrix_invariance_residual = _matrix_relative(
            section.laminate.D, baseline.laminate.D
        )
        A_formula_residual = _relative(
            section.properties.A / baseline.properties.A, alpha_A
        )
        D_formula_residual = _relative(
            section.properties.D / baseline.properties.D, 1.0
        )
        S_formula_residual = _relative(
            section.properties.S / baseline.properties.S, 1.0
        )
        m_formula_residual = _relative(section.properties.m / baseline.properties.m, 1.0)
        J_formula_residual = _relative(section.properties.J / baseline.properties.J, 1.0)
        B_relative = _scaled_B(section.laminate)
        I1_relative = _scaled_I1(section.laminate)
        only_A_field_max_invariance_residual = max(
            _relative(section.properties.D, baseline.properties.D),
            _relative(section.properties.S, baseline.properties.S),
            _relative(section.properties.m, baseline.properties.m),
            _relative(section.properties.J, baseline.properties.J),
            _relative(section.properties.K, baseline.properties.K),
            _relative(section.properties.width, baseline.properties.width),
        )
        reduction_route_max_relative = _reduction_max_residual(section.properties)
        Q_scaling_residual = max(
            _matrix_relative(Q_outer, s_outer * Q0),
            _matrix_relative(Q_inner, s_inner * Q0),
        )
        shear_Q_invariance_residual = max(
            _matrix_relative(shear_outer, shear0),
            _matrix_relative(shear_inner, shear0),
        )
        point_pass = bool(
            A_matrix_scale_residual <= MATRIX_RELATIVE_TOLERANCE
            and D_matrix_invariance_residual <= MATRIX_RELATIVE_TOLERANCE
            and A_formula_residual <= REDUCED_PROPERTY_TOLERANCE
            and D_formula_residual <= REDUCED_PROPERTY_TOLERANCE
            and S_formula_residual <= REDUCED_PROPERTY_TOLERANCE
            and m_formula_residual <= REDUCED_PROPERTY_TOLERANCE
            and J_formula_residual <= REDUCED_PROPERTY_TOLERANCE
            and B_relative <= SYMMETRY_RELATIVE_TOLERANCE
            and I1_relative <= SYMMETRY_RELATIVE_TOLERANCE
            and only_A_field_max_invariance_residual <= REDUCED_PROPERTY_TOLERANCE
            and reduction_route_max_relative <= REDUCTION_ROUTE_TOLERANCE
            and Q_scaling_residual <= MATRIX_RELATIVE_TOLERANCE
            and shear_Q_invariance_residual <= MATRIX_RELATIVE_TOLERANCE
        )
        passed = passed and point_pass
        for key, value in (
            ("A_matrix_scale_residual", A_matrix_scale_residual),
            ("D_matrix_invariance_residual", D_matrix_invariance_residual),
            ("A_formula_residual", A_formula_residual),
            ("D_formula_residual", D_formula_residual),
            ("S_formula_residual", S_formula_residual),
            ("m_formula_residual", m_formula_residual),
            ("J_formula_residual", J_formula_residual),
            ("B_relative", B_relative),
            ("I1_relative", I1_relative),
            ("only_A_field_max_invariance_residual", only_A_field_max_invariance_residual),
            ("reduction_route_max_relative", reduction_route_max_relative),
            ("Q_scaling_residual", Q_scaling_residual),
            ("shear_Q_invariance_residual", shear_Q_invariance_residual),
        ):
            maxima[key] = max(maxima[key], float(value))
        checks.append(
            {
                "alpha_A": alpha_A,
                "s_outer": s_outer,
                "s_inner": s_inner,
                "status": "PASS" if point_pass else "FAIL",
                "A_beam_over_A0": section.properties.A / baseline.properties.A,
                "D_beam_over_D0": section.properties.D / baseline.properties.D,
                "S_beam_over_S0": section.properties.S / baseline.properties.S,
                "m_over_m0": section.properties.m / baseline.properties.m,
                "J_over_J0": section.properties.J / baseline.properties.J,
                "A_matrix_scale_residual": A_matrix_scale_residual,
                "D_matrix_invariance_residual": D_matrix_invariance_residual,
                "A_formula_residual": A_formula_residual,
                "D_formula_residual": D_formula_residual,
                "Q_scaling_residual": Q_scaling_residual,
                "shear_Q_invariance_residual": shear_Q_invariance_residual,
                "B_relative": B_relative,
                "I1_relative": I1_relative,
                "only_A_field_max_invariance_residual": (
                    only_A_field_max_invariance_residual
                ),
                "reduction_route_max_relative": reduction_route_max_relative,
            }
        )
    return {
        "status": "PASS" if passed else "FAIL",
        "passed": bool(passed),
        "checks": checks,
        "maximum_residuals": maxima,
        "A0": baseline.properties.A,
        "D0": baseline.properties.D,
        "S0": baseline.properties.S,
        "m0": baseline.properties.m,
        "J0": baseline.properties.J,
        "reference_section": {
            "A_matrix": baseline.laminate.A,
            "D_matrix": baseline.laminate.D,
            "shear_matrix": baseline.laminate.shear,
        },
        "tolerances": {
            "matrix_relative": MATRIX_RELATIVE_TOLERANCE,
            "symmetry_relative": SYMMETRY_RELATIVE_TOLERANCE,
            "reduced_property_relative": REDUCED_PROPERTY_TOLERANCE,
            "reduction_route_relative": REDUCTION_ROUTE_TOLERANCE,
        },
    }


def section_property_rows() -> list[dict[str, Any]]:
    baseline = _baseline_section()
    rows: list[dict[str, Any]] = []
    for alpha_A_value in alpha_grid():
        alpha_A = float(alpha_A_value)
        section = build_layered_section(alpha_A)
        s_outer, s_inner = alpha_scales(alpha_A)
        B_relative = _scaled_B(section.laminate)
        I1_relative = _scaled_I1(section.laminate)
        row = {
            "alpha_A": alpha_A,
            "s_outer": s_outer,
            "s_inner": s_inner,
            "stack_bottom_to_top": list(section.layout),
            "properties_identical_between_arms": True,
            "A_beam": section.properties.A,
            "A_beam_over_A0": section.properties.A / baseline.properties.A,
            "expected_A_beam_over_A0": alpha_A,
            "D_beam": section.properties.D,
            "D_beam_over_D0": section.properties.D / baseline.properties.D,
            "S_beam": section.properties.S,
            "S_beam_over_S0": section.properties.S / baseline.properties.S,
            "m": section.properties.m,
            "m_over_m0": section.properties.m / baseline.properties.m,
            "J": section.properties.J,
            "J_over_J0": section.properties.J / baseline.properties.J,
            "B_relative": B_relative,
            "I1_relative": I1_relative,
            "A_matrix_scale_residual": _matrix_relative(
                section.laminate.A, alpha_A * baseline.laminate.A
            ),
            "D_matrix_invariance_residual": _matrix_relative(
                section.laminate.D, baseline.laminate.D
            ),
            "A_formula_residual": _relative(
                section.properties.A / baseline.properties.A, alpha_A
            ),
            "D_formula_residual": _relative(
                section.properties.D / baseline.properties.D, 1.0
            ),
            "only_A_field_max_invariance_residual": max(
                _relative(section.properties.D, baseline.properties.D),
                _relative(section.properties.S, baseline.properties.S),
                _relative(section.properties.m, baseline.properties.m),
                _relative(section.properties.J, baseline.properties.J),
                _relative(section.properties.K, baseline.properties.K),
                _relative(section.properties.width, baseline.properties.width),
            ),
            "reduction_route_max_relative": _reduction_max_residual(
                section.properties
            ),
            "z_interfaces": section.laminate.z_interfaces,
            "A_matrix": section.laminate.A,
            "B_matrix": section.laminate.B,
            "D_matrix": section.laminate.D,
            "shear_matrix_yz_xz": section.laminate.shear,
            "constitutive_status": "PASS",
        }
        rows.append(row)
    return rows


def _policy_root9() -> Any:
    return rlb2e._rlb2e_search_policy()


def _truncate_inventory_to_root9(
    canonical: Sequence[Any], slots: Sequence[Any], policy: Any
) -> tuple[list[Any], list[Any]]:
    """Retain nine export slots while preserving full guard-cluster metadata."""
    if len(slots) <= K_GUARD:
        return list(canonical), list(slots)
    guard_event = slots[K_GUARD - 1].event
    guard_cluster = guard_event.cluster_id or guard_event.event_id
    distinct_extra_slots = [
        slot
        for slot in slots[K_GUARD:]
        if (slot.event.cluster_id or slot.event.event_id) != guard_cluster
    ]
    if distinct_extra_slots:
        raise RuntimeError(
            "The bounded first-eight-plus-root9 search localized a distinct "
            "root above the guard; silently truncating it is forbidden."
        )
    guard = float(slots[K_GUARD - 1].event.omega_bar)
    tolerance = policy.dedup_atol_bar + policy.dedup_rtol * abs(guard)
    trimmed_canonical = [
        candidate
        for candidate in canonical
        if float(candidate.omega_bar) <= guard + tolerance
    ]
    return trimmed_canonical, list(slots[:K_GUARD])


def _contract_instance(beta_deg: float) -> dict[str, Any]:
    payload = dict(FREQUENCY_MAP_POLICY_BASE)
    payload["beta_deg"] = float(beta_deg)
    return payload


def make_matrix_provider(beta_deg: float, alpha_A: float) -> tuple[Any, dict[str, Any]]:
    _beam, coupled = _physics_modules()
    section = build_layered_section(alpha_A)
    joint = np.asarray(coupled.joint_matrix(math.radians(beta_deg)), dtype=float)

    def provider(omega: float) -> FloatArray:
        arm_map = np.asarray(
            coupled.arm_clamp_map(float(omega), L1, section.properties), dtype=float
        )
        combined = np.zeros((12, 6), dtype=float)
        combined[:6, :3] = arm_map
        combined[6:, 3:] = arm_map
        return np.asarray(joint @ combined, dtype=float)

    direct_residual = 0.0
    for probe in (0.731, 3.217):
        direct = coupled.coupled_boundary_matrix(
            probe,
            math.radians(beta_deg),
            L1,
            section.properties,
            L2,
            section.properties,
        )
        direct_residual = max(
            direct_residual,
            float(np.max(np.abs(provider(probe) - np.asarray(direct)))),
        )
    if direct_residual > 16.0 * np.finfo(float).eps:
        raise RuntimeError("RLB-2H provider differs from public coupled builder.")
    return provider, {
        "beta_deg": float(beta_deg),
        "alpha_A": float(alpha_A),
        "cached_vs_public_builder_max_abs": direct_residual,
        "production_modules_only": True,
        "identical_arms": True,
    }


def _root_rows(
    beta_deg: float,
    alpha_A: float,
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
    grid_kind: str = "BASE",
    repair_id: str = "",
) -> tuple[dict[str, Any], ...]:
    windows = (
        None if predicted is None else rlb2e.local_search_windows(predicted)
    )
    guard_event = slots[K_GUARD - 1].event
    guard = float(guard_event.omega_bar)
    guard_cluster_id = guard_event.cluster_id or guard_event.event_id
    represented_guard_slots = sum(
        (slot.event.cluster_id or slot.event.event_id) == guard_cluster_id
        for slot in slots[:K_GUARD]
    )
    guard_cluster_extends = bool(
        int(guard_event.cluster_multiplicity) > represented_guard_slots
    )
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
        rows.append(
            {
                "row_id": (
                    f"beta{int(beta_deg):02d}__alpha_{alpha_A:.6f}__{grid_kind}"
                    f"__p{position:02d}__{repair_id or 'base'}"
                ),
                "beta_deg": float(beta_deg),
                "alpha_A": float(alpha_A),
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
                "cluster_id": event.cluster_id or event.event_id,
                "cluster_semantics": event.cluster_semantics,
                "cluster_multiplicity": event.cluster_multiplicity,
                "cluster_total_nullity": event.cluster_total_nullity,
                "cluster_center_Omega": event.cluster_center_omega_bar,
                "cluster_metadata_source": "DETERMINANT_SVD_EVENT",
                "unresolved_candidates_below_root9": unresolved,
                "search_right_Omega": search_right,
                "root9_right_margin_Omega": search_right - guard,
                "guard_cluster_multiplicity": (
                    guard_event.cluster_multiplicity if position == K_GUARD else ""
                ),
                "guard_cluster_extends_beyond_export": (
                    guard_cluster_extends if position == K_GUARD else False
                ),
                "solve_mode": solve_mode,
                "fallback_used": fallback_used,
                "quality_status": "PASS",
                "point_status": "PASS",
                "is_canonical_plot_source": True,
                "supersedes_row_id": "",
                "repair_id": repair_id,
                "decoupled_reference_classification": "",
                "roots_above_9_computed": False,
            }
        )
    return tuple(rows)


def _point_is_acceptable_with_multiplicity(
    canonical: Sequence[Any],
    rejected: Sequence[Any],
    slots: Sequence[Any],
    search_right: float,
    policy: Any,
) -> tuple[bool, dict[str, Any]]:
    """Apply frozen quality gates while retaining a proven repeated slot."""

    predecessor_pass, predecessor_quality = rlb2e._point_is_acceptable(
        canonical, rejected, slots, search_right, policy
    )
    if predecessor_pass:
        quality = dict(predecessor_quality)
        quality["multiplicity_aware_order"] = True
        guard_event = slots[K_GUARD - 1].event
        guard_id = guard_event.cluster_id or guard_event.event_id
        represented = sum(
            (slot.event.cluster_id or slot.event.event_id) == guard_id
            for slot in slots[:K_GUARD]
        )
        quality["guard_cluster_multiplicity"] = int(
            guard_event.cluster_multiplicity
        )
        quality["guard_cluster_represented_slots"] = represented
        quality["guard_cluster_extends_beyond_export"] = bool(
            int(guard_event.cluster_multiplicity) > represented
        )
        return True, quality

    tools = rlb2e._root_tools()
    first = list(slots[:K_GUARD])
    values = [float(slot.event.omega_bar) for slot in first]
    order_tolerance = policy.dedup_atol_bar + policy.dedup_rtol * max(
        [1.0, *[abs(value) for value in values]]
    )
    nondecreasing = bool(
        len(first) == K_GUARD
        and all(left <= right + order_tolerance for left, right in zip(values[:-1], values[1:]))
    )
    repeated_pairs_are_clusters = True
    for left, right in zip(first[:-1], first[1:]):
        if float(left.event.omega_bar) < float(right.event.omega_bar):
            continue
        left_id = left.event.cluster_id or left.event.event_id
        right_id = right.event.cluster_id or right.event.event_id
        repeated_pairs_are_clusters = repeated_pairs_are_clusters and bool(
            left_id == right_id
            and int(left.event.cluster_multiplicity) >= 2
            and int(right.event.cluster_multiplicity) >= 2
            and int(left.event.cluster_total_nullity) >= 2
            and int(right.event.cluster_total_nullity) >= 2
        )
    qualities = [
        bool(tools._candidate_quality(slot.event.candidate.diagnostics, policy)[0])
        for slot in first
    ]
    guard = values[-1] if len(values) == K_GUARD else math.inf
    tolerance = policy.dedup_atol_bar + policy.dedup_rtol * abs(guard)
    accepted_above_guard = [
        candidate
        for candidate in canonical
        if float(candidate.omega_bar) > guard + tolerance
    ]
    unresolved = (
        rlb2e._unresolved_candidates_below(rejected, canonical, guard, policy)
        if len(first) == K_GUARD
        else -1
    )
    margin = float(search_right) - guard
    guard_cluster_id = (
        first[-1].event.cluster_id or first[-1].event.event_id
        if first
        else ""
    )
    represented_guard_slots = sum(
        (slot.event.cluster_id or slot.event.event_id) == guard_cluster_id
        for slot in first
    )
    guard_cluster_multiplicity = (
        int(first[-1].event.cluster_multiplicity) if first else 0
    )
    passed = bool(
        len(first) == K_GUARD
        and len(slots) == K_GUARD
        and not accepted_above_guard
        and nondecreasing
        and repeated_pairs_are_clusters
        and all(qualities)
        and unresolved == 0
        and margin >= ROOT9_RIGHT_TAIL_OMEGA - 1.0e-10
    )
    quality = {
        "root_count": len(first),
        "strict_order": all(left < right for left, right in zip(values[:-1], values[1:])),
        "multiplicity_aware_order": nondecreasing and repeated_pairs_are_clusters,
        "all_quality_pass": bool(qualities and all(qualities)),
        "unresolved_candidates_below_root9": unresolved,
        "root9_right_margin_Omega": margin,
        "candidate_count_total": len(canonical),
        "accepted_candidates_above_root9": len(accepted_above_guard),
        "retained_slots_above_root9": max(0, len(slots) - K_GUARD),
        "roots_above_9_computed": bool(accepted_above_guard or len(slots) > K_GUARD),
        "guard_cluster_multiplicity": guard_cluster_multiplicity,
        "guard_cluster_represented_slots": represented_guard_slots,
        "guard_cluster_extends_beyond_export": bool(
            guard_cluster_multiplicity > represented_guard_slots
        ),
    }
    return passed, quality


def solve_point(
    beta_deg: float,
    alpha_A: float,
    *,
    previous: tuple[float, Sequence[float]] | None = None,
    second_previous: tuple[float, Sequence[float]] | None = None,
    force_anchor: bool = False,
    dense_local: bool = False,
    dense_positions: Sequence[int] | None = None,
    grid_kind: str = "BASE",
    repair_id: str = "",
    continuation_leg: str = "ANCHOR",
) -> PointSolution:
    global ROOT_CALCULATION_COUNT
    ROOT_CALCULATION_COUNT += 1
    started = time.perf_counter()
    solve_id = f"beta{int(beta_deg):02d}__alpha_{alpha_A:.6f}"
    transaction_id = hashlib.sha256(
        f"{STAGE_ID}|{solve_id}|{grid_kind}|{repair_id}".encode("utf-8")
    ).hexdigest()[:20].upper()
    provider, _metadata = make_matrix_provider(beta_deg, alpha_A)
    counted = rlb2e.CountedProvider(provider)
    policy = _policy_root9()
    predicted: FloatArray | None = None
    fallback_used = False
    if not force_anchor and previous is not None:
        predicted = rlb2e.hold_secant_predictor(
            alpha_A,
            previous[0],
            previous[1],
            None if second_previous is None else second_previous[0],
            None if second_previous is None else second_previous[1],
        )
        predicted = np.sort(predicted)
        try:
            canonical, rejected, slots, search_right, refinements = rlb2e._local_candidates(
                counted,
                policy,
                predicted,
                solve_id=solve_id,
                dense=dense_local,
                dense_positions=dense_positions,
            )
            canonical, slots = _truncate_inventory_to_root9(canonical, slots, policy)
            accepted, quality = _point_is_acceptable_with_multiplicity(
                canonical, rejected, slots, search_right, policy
            )
        except (ValueError, RuntimeError, ArithmeticError, np.linalg.LinAlgError):
            accepted = False
            quality = {}
        if not accepted:
            fallback_used = True
            canonical, rejected, slots, search_right, refinements = (
                rlb2e._progressive_anchor_candidates(
                    counted, policy, solve_id=solve_id + "_fallback"
                )
            )
            canonical, slots = _truncate_inventory_to_root9(canonical, slots, policy)
            accepted, quality = _point_is_acceptable_with_multiplicity(
                canonical, rejected, slots, search_right, policy
            )
            solve_mode = "BOUNDED_GLOBAL_RECOVERY"
        else:
            solve_mode = "FAST_LOCAL"
    else:
        canonical, rejected, slots, search_right, refinements = (
            rlb2e._progressive_anchor_candidates(counted, policy, solve_id=solve_id)
        )
        canonical, slots = _truncate_inventory_to_root9(canonical, slots, policy)
        accepted, quality = _point_is_acceptable_with_multiplicity(
            canonical, rejected, slots, search_right, policy
        )
        solve_mode = "PROGRESSIVE_ANCHOR"
    if not accepted:
        raise RuntimeError(f"{solve_id}: first-eight-plus-root9 quality failed: {quality}")
    rows = _root_rows(
        beta_deg,
        alpha_A,
        slots,
        solve_id=solve_id,
        transaction_id=transaction_id,
        solve_mode=solve_mode,
        fallback_used=fallback_used,
        predicted=None if fallback_used else predicted,
        search_right=search_right,
        unresolved=int(quality["unresolved_candidates_below_root9"]),
        continuation_leg=continuation_leg,
        grid_kind=grid_kind,
        repair_id=repair_id,
    )
    return PointSolution(
        beta_deg=float(beta_deg),
        alpha_A=float(alpha_A),
        rows=rows,
        wall_time_seconds=time.perf_counter() - started,
        peak_rss_bytes=rlb2e._peak_rss_bytes(),
        determinant_evaluations=counted.evaluations,
        sigma_evaluations=counted.evaluations,
        search_left_Omega=(
            1.0e-8
            if predicted is None or fallback_used
            else rlb2e.local_search_windows(predicted)[0][0]
        ),
        search_right_Omega=search_right,
        local_refinements=refinements,
        solve_mode=solve_mode,
        fallback_used=fallback_used,
        unresolved_candidates_below_root9=int(
            quality["unresolved_candidates_below_root9"]
        ),
        continuation_leg=continuation_leg,
        candidate_count_total=int(quality["candidate_count_total"]),
        accepted_candidates_above_root9=int(
            quality["accepted_candidates_above_root9"]
        ),
        retained_slots_above_root9=int(quality["retained_slots_above_root9"]),
        roots_above_9_computed=bool(quality["roots_above_9_computed"]),
    )


def _as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _group_key(beta_deg: float, alpha_A: float) -> tuple[float, float]:
    return (round(float(beta_deg), 10), round(float(alpha_A), 10))


def _base_group_index(
    rows: Sequence[Mapping[str, Any]],
) -> dict[tuple[float, float], list[Mapping[str, Any]]]:
    groups: dict[tuple[float, float], list[Mapping[str, Any]]] = {}
    for row in rows:
        if str(row.get("grid_kind")) != "BASE":
            continue
        key = _group_key(float(row["beta_deg"]), float(row["alpha_A"]))
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
            and np.all(Omegas > 0.0)
            and np.all(np.isfinite(omegas))
            and np.all(omegas > 0.0)
            and np.all(np.isfinite(Lambdas))
            and np.all(Lambdas > 0.0)
            and np.all(np.diff(Omegas) >= 0.0)
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
        multiplicity_ok = all(
            float(left["Omega"]) < float(right["Omega"])
            or (
                str(left.get("cluster_id", ""))
                == str(right.get("cluster_id", ""))
                and bool(str(left.get("cluster_id", "")))
                and int(left.get("cluster_multiplicity", 1)) >= 2
                and int(right.get("cluster_multiplicity", 1)) >= 2
                and int(left.get("cluster_total_nullity", 1)) >= 2
                and int(right.get("cluster_total_nullity", 1)) >= 2
            )
            for left, right in zip(ordered[:-1], ordered[1:])
        )
        quality_ok = all(
            str(row["quality_status"]) == "PASS"
            and int(row["unresolved_candidates_below_root9"]) == 0
            and float(row["scaled_sigma_ratio"]) <= ROOT_SINGULAR_RATIO_TOLERANCE
            and float(row["boundary_null_residual"]) <= BOUNDARY_RESIDUAL_TOLERANCE
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
    return bool(roles_ok and numerical_ok and multiplicity_ok and quality_ok and guard_ok)


def _complete_base_group_index(
    rows: Sequence[Mapping[str, Any]],
) -> dict[tuple[float, float], list[Mapping[str, Any]]]:
    return {
        key: group
        for key, group in _base_group_index(rows).items()
        if _base_group_is_acceptable(group)
    }


def _canonical_group(
    rows: Sequence[Mapping[str, Any]], beta_deg: float, alpha_A: float
) -> list[Mapping[str, Any]]:
    selected = [
        row
        for row in rows
        if _group_key(float(row["beta_deg"]), float(row["alpha_A"]))
        == _group_key(beta_deg, alpha_A)
        and _as_bool(row.get("is_canonical_plot_source", True))
    ]
    selected.sort(key=lambda row: int(row["sorted_position"]))
    if [int(row["sorted_position"]) for row in selected] != list(
        range(1, K_GUARD + 1)
    ):
        raise RuntimeError(
            f"Incomplete canonical group at beta={beta_deg:.0f}, alpha_A={alpha_A:.2f}."
        )
    return selected


def _rows_for_roots(
    rows: Sequence[Mapping[str, Any]], beta_deg: float, alpha_A: float
) -> FloatArray:
    return np.asarray(
        [
            float(row["Omega"])
            for row in _canonical_group(rows, beta_deg, alpha_A)
        ],
        dtype=float,
    )


def _sort_rows(rows: list[dict[str, Any]]) -> None:
    rows.sort(
        key=lambda row: (
            float(row["beta_deg"]),
            float(row["alpha_A"]),
            0 if str(row["grid_kind"]) == "BASE" else 1,
            int(row["sorted_position"]),
            str(row.get("repair_id", "")),
        )
    )


def _spectrum_row_domain_failures(
    rows: Sequence[Mapping[str, Any]],
) -> list[str]:
    expected_points = {
        _group_key(beta_deg, alpha_A)
        for beta_deg in BETA_VALUES_DEG
        for alpha_A in alpha_grid()
    }
    failures: list[str] = []
    for row in rows:
        row_id = str(row.get("row_id", "<missing-row-id>"))
        try:
            key = _group_key(float(row["beta_deg"]), float(row["alpha_A"]))
            position = int(row["sorted_position"])
            grid_kind = str(row["grid_kind"])
        except (KeyError, TypeError, ValueError):
            failures.append(f"{row_id}:MALFORMED_ROW_DOMAIN")
            continue
        if key not in expected_points:
            failures.append(f"{row_id}:UNEXPECTED_PARAMETER_POINT")
        if position not in range(1, K_GUARD + 1):
            failures.append(f"{row_id}:INVALID_SORTED_POSITION")
        if grid_kind not in {"BASE", "LOCAL_REFINEMENT"}:
            failures.append(f"{row_id}:INVALID_GRID_KIND")
        if grid_kind == "LOCAL_REFINEMENT" and not str(
            row.get("repair_id", "")
        ).strip():
            failures.append(f"{row_id}:MISSING_REPAIR_ID")
    return failures


def audit_spectrum_rows(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    groups = _base_group_index(rows)
    complete = _complete_base_group_index(rows)
    expected = {
        _group_key(beta_deg, alpha)
        for beta_deg in BETA_VALUES_DEG
        for alpha in alpha_grid()
    }
    duplicates: list[str] = []
    incomplete: list[str] = []
    quality_failures: list[str] = []
    for key, group in groups.items():
        positions = [int(row["sorted_position"]) for row in group]
        label = f"beta={key[0]:.0f},alpha={key[1]:.2f}"
        if len(positions) != len(set(positions)):
            duplicates.append(label)
        if not _base_group_has_exact_positions(group):
            incomplete.append(label)
        elif not _base_group_is_acceptable(group):
            quality_failures.append(label)
    base = [row for row in rows if str(row.get("grid_kind")) == "BASE"]
    above = [row for row in rows if int(row["sorted_position"]) > K_GUARD]
    duplicate_row_ids = len(rows) - len({str(row["row_id"]) for row in rows})
    canonical_counts: dict[tuple[float, float, int], int] = {}
    for row in rows:
        if _as_bool(row.get("is_canonical_plot_source", False)):
            key = (
                round(float(row["beta_deg"]), 10),
                round(float(row["alpha_A"]), 10),
                int(row["sorted_position"]),
            )
            canonical_counts[key] = canonical_counts.get(key, 0) + 1
    expected_canonical = {
        (round(float(beta), 10), round(float(alpha), 10), position)
        for beta in BETA_VALUES_DEG
        for alpha in alpha_grid()
        for position in range(1, K_GUARD + 1)
    }
    canonical_failures = [
        f"beta={key[0]:.0f},alpha={key[1]:.2f}:p{key[2]}"
        for key in sorted(expected_canonical | set(canonical_counts))
        if key not in expected_canonical or canonical_counts.get(key, 0) != 1
    ]
    missing = sorted(expected - set(groups))
    extra = sorted(set(groups) - expected)
    row_domain_failures = _spectrum_row_domain_failures(rows)
    passed = bool(
        not duplicates
        and not incomplete
        and not quality_failures
        and not missing
        and not extra
        and not above
        and duplicate_row_ids == 0
        and not canonical_failures
        and not row_domain_failures
        and len(base) == len(BETA_VALUES_DEG) * len(alpha_grid()) * K_GUARD
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "base_group_count": len(complete),
        "base_row_count": len(base),
        "missing_groups": [f"beta={key[0]:.0f},alpha={key[1]:.2f}" for key in missing],
        "extra_groups": [f"beta={key[0]:.0f},alpha={key[1]:.2f}" for key in extra],
        "duplicate_groups": duplicates,
        "incomplete_groups": incomplete,
        "physical_quality_failures": quality_failures,
        "duplicate_row_id_count": duplicate_row_ids,
        "canonical_source_failures": canonical_failures,
        "row_domain_failures": row_domain_failures,
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
        (round(float(beta), 10), round(float(alpha), 10), position)
        for beta in BETA_VALUES_DEG
        for alpha in alpha_grid()
        for position in range(1, K_PLOT + 1)
    }
    counts: dict[tuple[float, float, int], int] = {}
    invalid: list[str] = []
    for row in selected:
        key = (
            round(float(row["beta_deg"]), 10),
            round(float(row["alpha_A"]), 10),
            int(row["sorted_position"]),
        )
        counts[key] = counts.get(key, 0) + 1
        value = float(row["Lambda"])
        gap = math.isnan(value) and str(row.get("point_status", "")) == (
            "UNRESOLVED_AFTER_LOCAL_REPAIR"
        )
        if not ((math.isfinite(value) and value > 0.0) or gap):
            invalid.append(f"beta={key[0]:.0f},alpha={key[1]:.2f}:p{key[2]}")
    failures = [
        f"beta={key[0]:.0f},alpha={key[1]:.2f}:p{key[2]}"
        for key in sorted(expected | set(counts))
        if key not in expected or counts.get(key, 0) != 1
    ]
    row_domain_failures = _spectrum_row_domain_failures(rows)
    return {
        "status": (
            "PASS"
            if not failures and not invalid and not row_domain_failures
            else "FAIL"
        ),
        "row_count": len(selected),
        "key_failures": failures,
        "invalid_values": invalid,
        "row_domain_failures": row_domain_failures,
    }


def audit_final_spectrum_schema(
    rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Check additive provenance fields required by a completed result tree."""

    missing_fields: list[str] = []
    blank_cluster_sources: list[str] = []
    incomplete_guards: list[str] = []
    invalid_beta0_classification: list[str] = []
    roots_above_guard_claims: list[str] = []
    required = set(SPECTRUM_FIELDS)
    allowed_beta0 = {"AXIAL", "BENDING", "DEGENERATE"}
    for row in rows:
        row_id = str(row.get("row_id", "<missing-row-id>"))
        absent = sorted(required - set(row))
        if absent:
            missing_fields.append(f"{row_id}:{','.join(absent)}")
            continue
        if not str(row.get("cluster_metadata_source", "")).strip():
            blank_cluster_sources.append(row_id)
        if int(row["sorted_position"]) == K_GUARD and (
            not str(row.get("guard_cluster_multiplicity", "")).strip()
            or not str(row.get("guard_cluster_extends_beyond_export", "")).strip()
        ):
            incomplete_guards.append(row_id)
        if (
            abs(float(row["beta_deg"])) <= 1.0e-12
            and _as_bool(row.get("is_canonical_plot_source", True))
            and str(row.get("decoupled_reference_classification", ""))
            not in allowed_beta0
        ):
            invalid_beta0_classification.append(row_id)
        if _as_bool(row.get("roots_above_9_computed", False)):
            roots_above_guard_claims.append(row_id)
    passed = not (
        missing_fields
        or blank_cluster_sources
        or incomplete_guards
        or invalid_beta0_classification
        or roots_above_guard_claims
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "missing_field_rows": missing_fields,
        "blank_cluster_metadata_source_rows": blank_cluster_sources,
        "incomplete_guard_metadata_rows": incomplete_guards,
        "invalid_beta0_classification_rows": invalid_beta0_classification,
        "roots_above_guard_claim_rows": roots_above_guard_claims,
    }


def _existing_rows(output_dir: Path) -> list[dict[str, Any]]:
    path = Path(output_dir) / SPECTRUM_FILENAME
    return [] if not path.is_file() else [dict(row) for row in rlb2e._read_csv(path)]


def _write_point_transaction(
    output_dir: Path,
    existing_rows: Sequence[Mapping[str, Any]],
    solution: PointSolution,
) -> list[dict[str, Any]]:
    rows = [dict(row) for row in existing_rows]
    key = _group_key(solution.beta_deg, solution.alpha_A)
    group = [
        row
        for row in rows
        if str(row["grid_kind"]) == "BASE"
        and _group_key(float(row["beta_deg"]), float(row["alpha_A"])) == key
    ]
    if group:
        if _base_group_is_acceptable(group):
            return rows
        rows = [
            row
            for row in rows
            if not (
                _group_key(float(row["beta_deg"]), float(row["alpha_A"])) == key
                and str(row["grid_kind"]) in {"BASE", "LOCAL_REFINEMENT"}
            )
        ]
    rows.extend(dict(row) for row in solution.rows)
    _sort_rows(rows)
    rlb2e._atomic_write_csv(Path(output_dir) / SPECTRUM_FILENAME, rows, SPECTRUM_FIELDS)
    return rows


def _solution_record(solution: PointSolution, *, benchmark: bool) -> dict[str, Any]:
    return {
        "beta_deg": solution.beta_deg,
        "alpha_A": solution.alpha_A,
        "benchmark": benchmark,
        "continuation_leg": solution.continuation_leg,
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
        "accepted_candidates_above_root9": solution.accepted_candidates_above_root9,
        "retained_slots_above_root9": solution.retained_slots_above_root9,
        "roots_above_9_computed": solution.roots_above_9_computed,
        "guard_cluster_multiplicity": int(
            solution.rows[K_GUARD - 1]["guard_cluster_multiplicity"]
        ),
        "guard_cluster_extends_beyond_export": _as_bool(
            solution.rows[K_GUARD - 1]["guard_cluster_extends_beyond_export"]
        ),
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


def _benchmark_payload(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    measured = [
        float(item["wall_time_seconds"])
        for item in records
        if int(item.get("determinant_evaluations", 0)) > 0
    ]
    remaining = max(0, len(BETA_VALUES_DEG) * len(alpha_grid()) - len(records))
    max_time = max(measured, default=0.0)
    eta = remaining * max_time
    return {
        "anchors": [dict(item) for item in records],
        "anchor_count": len(records),
        "remaining_unique_root_points": remaining,
        "conservative_eta_seconds": eta,
        "eta_limit_seconds": ETA_LIMIT_SECONDS,
        "production_run_permitted": eta <= ETA_LIMIT_SECONDS,
        "peak_rss_bytes": max((int(item.get("peak_rss_bytes", 0)) for item in records), default=0),
    }


def _checkpoint_payload(
    rows: Sequence[Mapping[str, Any]],
    point_records: Sequence[Mapping[str, Any]],
    *,
    constitutive: Mapping[str, Any],
    started_at: str,
    benchmark_status: str,
) -> dict[str, Any]:
    groups = _complete_base_group_index(rows)
    expected = [
        {"beta_deg": beta_deg, "alpha_A": float(alpha)}
        for beta_deg in BETA_VALUES_DEG
        for alpha in alpha_grid()
    ]
    completed = [
        item for item in expected if _group_key(item["beta_deg"], item["alpha_A"]) in groups
    ]
    missing = [
        item for item in expected if _group_key(item["beta_deg"], item["alpha_A"]) not in groups
    ]
    return {
        "schema_version": 1,
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "started_at_utc": started_at,
        "updated_at_utc": rlb2e._utc_now(),
        "contract_sha256": contract_hash(),
        "completed_points": completed,
        "missing_points": missing,
        "point_records": list(point_records),
        "benchmark_status": benchmark_status,
        "constitutive_gate_status": constitutive["status"],
        "root_calculation_count": ROOT_CALCULATION_COUNT,
        "parallel_workers_used": 0,
        "thread_limits": {
            name: os.environ[name]
            for name in (
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
            )
        },
    }


def _validated_checkpoint_point_records(
    checkpoint: Mapping[str, Any],
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    if (
        checkpoint.get("schema_version") != 1
        or checkpoint.get("stage_id") != STAGE_ID
        or checkpoint.get("algorithm_version") != ALGORITHM_VERSION
        or checkpoint.get("contract_sha256") != contract_hash()
    ):
        raise RuntimeError("Checkpoint contract differs from RLB-2H.")
    raw_records = checkpoint.get("point_records")
    if not isinstance(raw_records, list):
        raise RuntimeError("Checkpoint point_records must be a list.")
    expected = {
        _group_key(beta_deg, alpha_A)
        for beta_deg in BETA_VALUES_DEG
        for alpha_A in alpha_grid()
    }
    complete = _complete_base_group_index(rows)
    base_groups = _base_group_index(rows)
    records: list[dict[str, Any]] = []
    seen: set[tuple[float, float]] = set()
    for item in raw_records:
        if not isinstance(item, Mapping):
            raise RuntimeError("Checkpoint contains a malformed point record.")
        try:
            key = _group_key(float(item["beta_deg"]), float(item["alpha_A"]))
        except (KeyError, TypeError, ValueError) as exc:
            raise RuntimeError("Checkpoint contains a malformed point key.") from exc
        if key not in expected or key not in complete or key in seen:
            raise RuntimeError(
                "Checkpoint point_records disagree with completed spectrum groups."
            )
        required_record_fields = {
            "benchmark",
            "continuation_leg",
            "wall_time_seconds",
            "peak_rss_bytes",
            "determinant_evaluations",
            "sigma_evaluations",
            "solve_mode",
            "fallback_used",
            "unresolved_candidates_below_root9",
            "roots",
        }
        if not required_record_fields <= set(item):
            raise RuntimeError("Checkpoint point record is missing required fields.")
        try:
            numeric_metrics = (
                float(item["wall_time_seconds"]),
                int(item["peak_rss_bytes"]),
                int(item["determinant_evaluations"]),
                int(item["sigma_evaluations"]),
                int(item["unresolved_candidates_below_root9"]),
            )
        except (TypeError, ValueError) as exc:
            raise RuntimeError("Checkpoint point record has invalid metrics.") from exc
        if any(value < 0 for value in numeric_metrics):
            raise RuntimeError("Checkpoint point record has negative metrics.")
        roots = item.get("roots")
        if (
            not isinstance(roots, list)
            or len(roots) != K_GUARD
            or not all(isinstance(root, Mapping) for root in roots)
            or [int(root.get("sorted_position", -1)) for root in roots]
            != list(range(1, K_GUARD + 1))
        ):
            raise RuntimeError("Checkpoint point record has an invalid root inventory.")
        base_rows = sorted(
            base_groups[key], key=lambda row: int(row["sorted_position"])
        )
        for stored, spectrum in zip(roots, base_rows, strict=True):
            try:
                consistent = bool(
                    _relative(float(stored["Omega"]), float(spectrum["Omega"]))
                    <= 1.0e-14
                    and _relative(
                        float(stored["Lambda"]), float(spectrum["Lambda"])
                    )
                    <= 1.0e-14
                    and _relative(
                        float(stored["singular_ratio"]),
                        float(spectrum["scaled_sigma_ratio"]),
                    )
                    <= 1.0e-14
                    and _relative(
                        float(stored["boundary_residual"]),
                        float(spectrum["boundary_null_residual"]),
                    )
                    <= 1.0e-14
                    and str(stored["quality_status"])
                    == str(spectrum["quality_status"])
                )
            except (KeyError, TypeError, ValueError) as exc:
                raise RuntimeError(
                    "Checkpoint point root has malformed diagnostics."
                ) from exc
            if not consistent:
                raise RuntimeError(
                    "Checkpoint point root disagrees with the stored BASE spectrum."
                )
        seen.add(key)
        records.append(dict(item))
    return records


def prepare_constitutive_outputs(output_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    gate = constitutive_gate()
    rows = section_property_rows()
    rlb2e._atomic_write_csv(target / SECTION_FILENAME, rows, SECTION_FIELDS)
    if gate["status"] != "PASS":
        raise RuntimeError(f"Constitutive gate failed before roots: {gate}")
    return gate, rows


def _complete_point_sequence(
    target: Path,
    rows: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    *,
    beta_deg: float,
    values: Sequence[float],
    leg_name: str,
    constitutive: Mapping[str, Any],
    started_at: str,
) -> list[dict[str, Any]]:
    previous: tuple[float, Sequence[float]] | None = None
    second_previous: tuple[float, Sequence[float]] | None = None
    groups = _complete_base_group_index(rows)
    for raw_alpha in values:
        alpha_A = float(raw_alpha)
        key = _group_key(beta_deg, alpha_A)
        if key in groups:
            roots = _rows_for_roots(rows, beta_deg, alpha_A)
            second_previous = previous
            previous = (alpha_A, roots)
            continue
        solution = solve_point(
            beta_deg,
            alpha_A,
            previous=previous,
            second_previous=second_previous,
            force_anchor=(previous is None),
            continuation_leg=leg_name,
        )
        rows = _write_point_transaction(target, rows, solution)
        point_records.append(_solution_record(solution, benchmark=False))
        groups = _complete_base_group_index(rows)
        second_previous = previous
        previous = (alpha_A, _rows_for_roots(rows, beta_deg, alpha_A))
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


def run_benchmarks(
    target: Path,
    rows: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    started_at: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    required = (
        (0.0, 1.0, "ANCHOR"),
        (0.0, 0.70, "LOWER"),
        (0.0, 1.30, "UPPER"),
        (30.0, 1.0, "ANCHOR"),
        (30.0, 0.70, "LOWER"),
        (30.0, 1.30, "UPPER"),
    )
    groups = _complete_base_group_index(rows)
    benchmark_records: list[dict[str, Any]] = []
    existing_index = {
        (round(float(item["beta_deg"]), 10), round(float(item["alpha_A"]), 10)): item
        for item in point_records
        if bool(item.get("benchmark", False))
    }
    for beta_deg, alpha_A, leg in required:
        key = _group_key(beta_deg, alpha_A)
        if key in existing_index:
            benchmark_records.append(dict(existing_index[key]))
            continue
        if key in groups:
            record = {
                "beta_deg": beta_deg,
                "alpha_A": alpha_A,
                "benchmark": True,
                "continuation_leg": leg,
                "wall_time_seconds": 0.0,
                "peak_rss_bytes": 0,
                "determinant_evaluations": 0,
                "sigma_evaluations": 0,
                "initial_Omega_max": "",
                "actual_Omega_max": "",
                "local_refinements": 0,
                "solve_mode": "REUSED_EXISTING_COMPLETE_GROUP",
                "fallback_used": False,
                "unresolved_candidates_below_root9": 0,
                "candidate_count_total": K_GUARD,
                "accepted_candidates_above_root9": 0,
                "retained_slots_above_root9": 0,
                "roots_above_9_computed": False,
                "roots": [
                    {
                        "sorted_position": int(row["sorted_position"]),
                        "Omega": float(row["Omega"]),
                        "Lambda": float(row["Lambda"]),
                        "singular_ratio": float(row["scaled_sigma_ratio"]),
                        "boundary_residual": float(row["boundary_null_residual"]),
                        "quality_status": str(row["quality_status"]),
                    }
                    for row in _canonical_group(rows, beta_deg, alpha_A)
                ],
            }
            benchmark_records.append(record)
            point_records.append(record)
            continue
        solution = solve_point(
            beta_deg,
            alpha_A,
            force_anchor=True,
            continuation_leg=leg,
        )
        rows = _write_point_transaction(target, rows, solution)
        record = _solution_record(solution, benchmark=True)
        benchmark_records.append(record)
        point_records.append(record)
        rlb2e._atomic_write_json(target / BENCHMARK_FILENAME, _benchmark_payload(benchmark_records))
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
    payload = _benchmark_payload(benchmark_records)
    rlb2e._atomic_write_json(target / BENCHMARK_FILENAME, payload)
    rlb2e._atomic_write_json(
        target / CHECKPOINT_FILENAME,
        _checkpoint_payload(
            rows,
            point_records,
            constitutive=constitutive,
            started_at=started_at,
            benchmark_status=(
                "PASS" if payload["production_run_permitted"] else "STOPPED_BY_ETA_GATE"
            ),
        ),
    )
    return rows, payload


def complete_missing_points(
    target: Path,
    rows: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    started_at: str,
) -> list[dict[str, Any]]:
    lower, upper = continuation_paths()
    for beta_deg in BETA_VALUES_DEG:
        rows = _complete_point_sequence(
            target,
            rows,
            point_records,
            beta_deg=beta_deg,
            values=lower,
            leg_name="LOWER",
            constitutive=constitutive,
            started_at=started_at,
        )
        rows = _complete_point_sequence(
            target,
            rows,
            point_records,
            beta_deg=beta_deg,
            values=upper[1:],
            leg_name="UPPER",
            constitutive=constitutive,
            started_at=started_at,
        )
    return rows


def neighbour_audit_rows(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    grid = [round(float(value), 10) for value in alpha_grid()]
    groups = _base_group_index(rows)
    for beta_deg in BETA_VALUES_DEG:
        spectra = {
            alpha: np.sqrt(_rows_for_roots(rows, beta_deg, alpha)[:K_PLOT]) for alpha in grid
        }
        gap_flags: set[tuple[float, int]] = set()
        for lower_position in range(1, K_PLOT):
            gaps = np.asarray(
                [
                    spectra[alpha][lower_position] - spectra[alpha][lower_position - 1]
                    for alpha in grid
                ],
                dtype=float,
            )
            residuals = np.asarray(
                [
                    rlb2e.centred_secant_residual(gaps[index - 1], gaps[index], gaps[index + 1])
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
        for position in range(1, K_PLOT + 1):
            values = {alpha: float(spectra[alpha][position - 1]) for alpha in grid}
            residuals = [
                rlb2e.centred_secant_residual(
                    values[grid[index - 1]], values[grid[index]], values[grid[index + 1]]
                )
                for index in range(1, len(grid) - 1)
            ]
            median = float(np.median(residuals))
            mad = float(np.median(np.abs(np.asarray(residuals) - median)))
            robust_threshold = median + NEIGHBOUR_MAD_MULTIPLIER * mad
            for offset, index in enumerate(range(1, len(grid) - 1)):
                alpha = grid[index]
                group = groups[_group_key(beta_deg, alpha)]
                ordered = sorted(group, key=lambda row: int(row["sorted_position"]))
                root_row = ordered[position - 1]
                residual = float(residuals[offset])
                statistical_flag = bool(
                    residual > robust_threshold
                    and residual > NEIGHBOUR_ABSOLUTE_TRIGGER
                )
                root_count_warning = len(ordered) != K_GUARD
                ordering_warning = any(
                    float(left["Omega"]) > float(right["Omega"])
                    for left, right in zip(ordered[:-1], ordered[1:], strict=True)
                )
                unresolved_warning = int(root_row["unresolved_candidates_below_root9"]) != 0
                bad_residual_warning = bool(
                    float(root_row["scaled_sigma_ratio"]) > ROOT_SINGULAR_RATIO_TOLERANCE
                    or float(root_row["boundary_null_residual"]) > BOUNDARY_RESIDUAL_TOLERANCE
                )
                gap_warning = (alpha, position) in gap_flags
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
                        "beta_deg": beta_deg,
                        "sorted_position": position,
                        "alpha_left": grid[index - 1],
                        "alpha_A": alpha,
                        "alpha_right": grid[index + 1],
                        "Lambda_left": values[grid[index - 1]],
                        "Lambda_center": values[alpha],
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
                        "local_alpha_values": [],
                        "smoothing_applied": False,
                    }
                )
    return result


def flagged_repair_points(
    audit_rows: Sequence[Mapping[str, Any]],
) -> list[tuple[float, float]]:
    return sorted(
        {
            _group_key(float(row["beta_deg"]), float(row["alpha_A"]))
            for row in audit_rows
            if _as_bool(row["flagged"])
        }
    )


def apply_local_repairs(
    rows: list[dict[str, Any]], audit_rows: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    records: list[dict[str, Any]] = []
    for repair_index, (beta_deg, alpha_A) in enumerate(flagged_repair_points(audit_rows), start=1):
        repair_id = f"repair_{repair_index:04d}"
        positions = sorted(
            {
                int(row["sorted_position"])
                for row in audit_rows
                if _group_key(float(row["beta_deg"]), float(row["alpha_A"]))
                == _group_key(beta_deg, alpha_A)
                and _as_bool(row["flagged"])
            }
        )
        original = _rows_for_roots(rows, beta_deg, alpha_A)
        try:
            solution = solve_point(
                beta_deg,
                alpha_A,
                previous=(alpha_A, original),
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
                    and _group_key(float(row["beta_deg"]), float(row["alpha_A"]))
                    == _group_key(beta_deg, alpha_A)
                    and int(row["sorted_position"]) in positions
                ):
                    row["is_canonical_plot_source"] = False
                    row["point_status"] = "UNRESOLVED_AFTER_LOCAL_REPAIR"
                    gap = dict(row)
                    gap["row_id"] = (
                        f"beta{int(beta_deg):02d}__alpha_{alpha_A:.6f}__LOCAL_REFINEMENT"
                        f"__p{int(row['sorted_position']):02d}__{repair_id}_gap"
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
                if _group_key(float(audit["beta_deg"]), float(audit["alpha_A"])) == _group_key(beta_deg, alpha_A) and _as_bool(audit["flagged"]):
                    audit["repair_id"] = repair_id
                    audit["repair_status"] = "UNRESOLVED_AFTER_LOCAL_REPAIR"
            records.append(
                {
                    "repair_id": repair_id,
                    "beta_deg": beta_deg,
                    "alpha_A": alpha_A,
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
                    (np.abs(original), np.abs(refined), np.full(K_GUARD, np.finfo(float).tiny))
                )
            )
        )
        status = (
            "REPRODUCED_AFTER_LOCAL_REPAIR"
            if relative <= 1.0e-8
            else "LOCATOR_CORRECTED_AFTER_LOCAL_REPAIR"
        )
        for row in rows:
            if (
                str(row["grid_kind"]) == "BASE"
                and _group_key(float(row["beta_deg"]), float(row["alpha_A"]))
                == _group_key(beta_deg, alpha_A)
            ):
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
                and _group_key(float(item["beta_deg"]), float(item["alpha_A"]))
                == _group_key(beta_deg, alpha_A)
                and int(item["sorted_position"]) == int(row["sorted_position"])
            )
            row["point_status"] = status
        rows.extend(dict(row) for row in solution.rows)
        for audit in audit_rows:
            if _group_key(float(audit["beta_deg"]), float(audit["alpha_A"])) == _group_key(beta_deg, alpha_A) and _as_bool(audit["flagged"]):
                audit["repair_id"] = repair_id
                audit["repair_status"] = status
        records.append(
            {
                "repair_id": repair_id,
                "beta_deg": beta_deg,
                "alpha_A": alpha_A,
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


def _beta0_reference_case(alpha_A: float) -> Any:
    section = build_layered_section(alpha_A)
    return SimpleNamespace(
        properties=section.properties,
        length=L_TOTAL,
        frequency_scale=OMEGA_TO_OMEGA_SCALE,
    )


def _beta0_family_reference_data(
    spectrum_rows: Sequence[Mapping[str, Any]],
    reference_Omega_max: float,
) -> dict[str, Any]:
    """Compute one bending reference and exact axial families for beta=0."""

    started = time.perf_counter()
    beam, _coupled = _physics_modules()
    baseline = _beta0_reference_case(1.0)
    axial_base = beam.exact_axial_modes(
        baseline.properties, L_TOTAL, "FF", n_modes=K_GUARD
    )
    alpha1_rows = _canonical_group(spectrum_rows, 0.0, 1.0)
    axial_Omega = [
        float(mode.omega) * OMEGA_TO_OMEGA_SCALE for mode in axial_base
    ]
    bending: list[Any] = []
    axial_matches: list[dict[str, Any]] = []
    for row in alpha1_rows:
        Omega = float(row["Omega"])
        nearest_index = int(np.argmin(np.abs(np.asarray(axial_Omega) - Omega)))
        relative = _relative(Omega, axial_Omega[nearest_index])
        if relative <= DIRECT_REFERENCE_RELATIVE_TOLERANCE:
            axial_matches.append(
                {
                    "coupled_sorted_position": int(row["sorted_position"]),
                    "axial_family_index": nearest_index + 1,
                    "relative_difference": relative,
                }
            )
            continue
        bending.append(
            SimpleNamespace(
                omega=float(row["omega"]),
                detection=str(row["detector_refiner_provenance"]),
                boundary_residual=float(row["boundary_null_residual"]),
                sigma_ratio=float(row["scaled_sigma_ratio"]),
                omega_bar=Omega,
            )
        )
    reconciliation = {
        "source": "COMPLETED_ALPHA1_COUPLED_CHARACTERISTIC_MATRIX_GROUP",
        "source_root_count": len(alpha1_rows),
        "exact_axial_matches": axial_matches,
        "bending_root_count": len(bending),
        "bounded_reference_max_Omega": float(reference_Omega_max),
        "root_values_selected_by": "existing_completed_alpha1_determinant_SVD_refiner",
        "predictor_or_interpolation_used_as_frequency": False,
        "completed_sweep_group_recomputed": False,
        "direct_fixed_fixed_matrix_used_for_independent_singularity_check": True,
        "two_by_two_bending_minimum_used_for_root_localization": False,
    }
    if len(bending) < K_PLOT:
        raise RuntimeError("The single beta=0 bending reference did not retain eight roots.")

    bending_matrix_invariance = 0.0
    axial_scaling_error = 0.0
    for alpha_A in (0.70, 1.00, 1.30):
        section = _beta0_reference_case(alpha_A)
        axial = beam.exact_axial_modes(
            section.properties, L_TOTAL, "FF", n_modes=K_GUARD
        )
        for reference_mode, mode in zip(axial_base, axial, strict=True):
            expected = reference_mode.omega * math.sqrt(alpha_A)
            axial_scaling_error = max(
                axial_scaling_error, _relative(mode.omega, expected)
            )
        for detection in bending:
            baseline_matrix = beam.bending_boundary_matrix(
                detection.omega, L_TOTAL, baseline.properties, "CC"
            )
            current_matrix = beam.bending_boundary_matrix(
                detection.omega, L_TOTAL, section.properties, "CC"
            )
            bending_matrix_invariance = max(
                bending_matrix_invariance,
                _matrix_relative(current_matrix, baseline_matrix),
            )
    return {
        "axial_base": axial_base,
        "bending": tuple(bending),
        "reference_detector_reconciliation": reconciliation,
        "reference_Omega_max": float(reference_Omega_max),
        "bending_matrix_invariance_error": bending_matrix_invariance,
        "axial_scaling_error": axial_scaling_error,
        "wall_time_seconds": time.perf_counter() - started,
    }


def _union_reference_rows_for_alpha(
    alpha_A: float, family_data: Mapping[str, Any]
) -> list[dict[str, Any]]:
    beam, _coupled = _physics_modules()
    policy = _policy_root9()
    axial = beam.exact_axial_modes(
        build_layered_section(alpha_A).properties,
        L_TOTAL,
        "FF",
        n_modes=K_GUARD,
    )
    bending = tuple(family_data["bending"])
    clusters = beam.union_subsystem_spectra(
        [mode.omega for mode in axial],
        [item.omega for item in bending],
        atol=policy.cluster_atol_bar / OMEGA_TO_OMEGA_SCALE,
        rtol=policy.cluster_rtol,
    )
    rows: list[dict[str, Any]] = []
    slot = 0
    for cluster_index, cluster in enumerate(clusters, start=1):
        cluster_id = f"alpha_{alpha_A:.6f}__cluster_{cluster_index:03d}"
        semantics = (
            "ISOLATED"
            if cluster.multiplicity == 1
            else (
                "EXACT_DEGENERATE_SUBSPACE"
                if cluster.exact_degeneracy
                else "NEAR_DEGENERATE_CLUSTER"
            )
        )
        for member_index, member in enumerate(cluster.members, start=1):
            if slot >= K_GUARD:
                break
            slot += 1
            Omega = float(member.omega) * OMEGA_TO_OMEGA_SCALE
            family_index = int(member.subsystem_index) + 1
            rows.append(
                {
                    "alpha_A": float(alpha_A),
                    "sorted_position": slot,
                    "root_role": "PLOTTED" if slot <= K_PLOT else "ROOT_9_GUARD",
                    "family": str(member.subsystem),
                    "family_index": family_index,
                    "classification": (
                        "DEGENERATE"
                        if cluster.multiplicity > 1
                        else str(member.subsystem).upper()
                    ),
                    "omega": float(member.omega),
                    "Omega": Omega,
                    "Lambda": Omega_to_Lambda(Omega),
                    "cluster_id": cluster_id,
                    "cluster_semantics": semantics,
                    "cluster_multiplicity": cluster.multiplicity,
                    "cluster_total_nullity": cluster.multiplicity,
                    "cluster_center_Omega": (
                        float(cluster.representative_omega) * OMEGA_TO_OMEGA_SCALE
                    ),
                    "cluster_member_slot": member_index,
                }
            )
        if slot >= K_GUARD:
            break
    if len(rows) != K_GUARD:
        raise RuntimeError(
            f"beta=0 subsystem union at alpha_A={alpha_A:.2f} has {len(rows)} slots."
        )
    return rows


def beta0_reference_rows(
    spectrum_rows: Sequence[Mapping[str, Any]],
    *,
    include_internal: bool = False,
) -> Any:
    root9_values = [
        float(row["Omega"])
        for row in spectrum_rows
        if str(row.get("grid_kind")) == "BASE"
        and abs(float(row["beta_deg"])) <= 1.0e-12
        and int(row["sorted_position"]) == K_GUARD
    ]
    if len(root9_values) != len(alpha_grid()):
        raise RuntimeError("beta=0 reference requires all 31 completed root-9 guards.")
    reference_Omega_max = max(root9_values) + ROOT9_RIGHT_TAIL_OMEGA
    family_data = _beta0_family_reference_data(spectrum_rows, reference_Omega_max)
    rows = [
        row
        for alpha_A in alpha_grid()
        for row in _union_reference_rows_for_alpha(float(alpha_A), family_data)
    ]
    axial_scaling_error = float(family_data["axial_scaling_error"])
    bending_invariance_error = float(
        family_data["bending_matrix_invariance_error"]
    )
    axial_base = {
        int(row["family_index"]): float(row["Lambda"])
        for row in rows
        if abs(float(row["alpha_A"]) - 1.0) < 1.0e-14
        and str(row["family"]) == "axial"
    }
    bending_base = {
        int(row["family_index"]): float(row["Lambda"])
        for row in rows
        if abs(float(row["alpha_A"]) - 1.0) < 1.0e-14
        and str(row["family"]) == "bending"
    }
    for row in rows:
        family_index = int(row["family_index"])
        alpha_A = float(row["alpha_A"])
        if str(row["family"]) == "axial":
            expected = axial_base.get(family_index)
            if expected is not None:
                expected *= alpha_A**0.25
                axial_scaling_error = max(
                    axial_scaling_error,
                    _relative(float(row["Lambda"]), expected),
                )
        else:
            expected = bending_base.get(family_index)
            if expected is not None:
                bending_invariance_error = max(
                    bending_invariance_error,
                    _relative(float(row["Lambda"]), expected),
                )
        row["expected_Lambda_from_alpha1"] = "" if expected is None else expected
        row["relative_error_vs_reference"] = (
            "" if expected is None else _relative(float(row["Lambda"]), expected)
        )
    passed = bool(
        bending_invariance_error <= DIRECT_REFERENCE_RELATIVE_TOLERANCE
        and axial_scaling_error <= DIRECT_REFERENCE_RELATIVE_TOLERANCE
    )
    summary = {
        "status": "PASS" if passed else "FAIL",
        "bending_reference_solve_count": 0,
        "bending_reference_source_group_count": 1,
        "bending_reference_construction_count": 1,
        "bending_reference_root_count": len(family_data["bending"]),
        "reference_row_count": len(rows),
        "bending_invariance_error": bending_invariance_error,
        "bending_matrix_invariance_error": float(
            family_data["bending_matrix_invariance_error"]
        ),
        "axial_scaling_error": axial_scaling_error,
        "reference_detector_reconciliation": family_data[
            "reference_detector_reconciliation"
        ],
        "reference_Omega_max": family_data["reference_Omega_max"],
        "wall_time_seconds": family_data["wall_time_seconds"],
        "root_solver_invocations": 0,
        "determinant_evaluations": 0,
        "sigma_evaluations": 0,
        "completed_alpha1_root_group_reused": True,
        "roots_above_global_guard_exported": False,
    }
    if include_internal:
        return rows, summary, {}
    return rows, summary


def _cluster_groups(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    ordered = sorted(rows, key=lambda row: int(row["sorted_position"]))
    groups: list[dict[str, Any]] = []
    for row in ordered:
        cluster_id = str(row.get("cluster_id", "")) or (
            f"slot_{int(row['sorted_position']):02d}"
        )
        if groups and groups[-1]["cluster_id"] == cluster_id:
            groups[-1]["last_slot"] = int(row["sorted_position"])
            groups[-1]["represented_values"].append(float(row["Omega"]))
            stored_center = row.get("cluster_center_Omega", "")
            if str(stored_center).strip():
                groups[-1]["stored_centers"].append(float(stored_center))
            continue
        stored_center = row.get("cluster_center_Omega", "")
        groups.append(
            {
                "cluster_id": cluster_id,
                "first_slot": int(row["sorted_position"]),
                "last_slot": int(row["sorted_position"]),
                "semantics": str(row.get("cluster_semantics") or "ISOLATED"),
                "multiplicity": int(row.get("cluster_multiplicity") or 1),
                "total_nullity": int(row.get("cluster_total_nullity") or 1),
                "represented_values": [float(row["Omega"])],
                "stored_centers": (
                    [] if not str(stored_center).strip() else [float(stored_center)]
                ),
            }
        )
    for group in groups:
        stored_centers = group.pop("stored_centers")
        represented = group.pop("represented_values")
        if group["multiplicity"] > 1 and stored_centers:
            spread = max(stored_centers) - min(stored_centers)
            scale = max(1.0, *(abs(value) for value in stored_centers))
            if spread > 1.0e-12 * scale:
                raise RuntimeError("Inconsistent stored cluster-center metadata.")
            group["center_Omega"] = float(np.mean(stored_centers))
            group["cluster_center_source"] = "STORED_FULL_CLUSTER_CENTER"
        else:
            group["center_Omega"] = float(np.mean(represented))
            group["cluster_center_source"] = "REPRESENTED_SLOT_VALUE"
        group["represented_slot_count"] = len(represented)
    return groups


def _compare_cluster_groups(
    left_rows: Sequence[Mapping[str, Any]],
    right_rows: Sequence[Mapping[str, Any]],
    *,
    alpha_A: float,
    comparison_kind: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    left_groups = _cluster_groups(left_rows)
    right_groups = _cluster_groups(right_rows)
    records: list[dict[str, Any]] = []
    maximum = 0.0
    passed = len(left_groups) == len(right_groups)
    for index in range(max(len(left_groups), len(right_groups))):
        left = left_groups[index] if index < len(left_groups) else None
        right = right_groups[index] if index < len(right_groups) else None
        if left is None or right is None:
            records.append(
                {
                    "alpha_A": alpha_A,
                    "comparison_kind": comparison_kind,
                    "comparison_index": index + 1,
                    "status": "FAIL",
                    "reason": "MISSING_CLUSTER",
                }
            )
            passed = False
            continue
        clustered = left["multiplicity"] > 1 or right["multiplicity"] > 1
        tolerance = (
            CLUSTER_REFERENCE_RELATIVE_TOLERANCE
            if clustered
            else DIRECT_REFERENCE_RELATIVE_TOLERANCE
        )
        relative = _relative(left["center_Omega"], right["center_Omega"])
        maximum = max(maximum, relative)
        row_pass = bool(
            relative <= tolerance
            and left["multiplicity"] == right["multiplicity"]
            and left["total_nullity"] == right["total_nullity"]
            and left["first_slot"] == right["first_slot"]
            and left["last_slot"] == right["last_slot"]
        )
        passed = passed and row_pass
        records.append(
            {
                "alpha_A": alpha_A,
                "comparison_kind": comparison_kind,
                "comparison_index": index + 1,
                "comparison_unit": "CLUSTER" if clustered else "ISOLATED",
                "left_first_slot": left["first_slot"],
                "left_last_slot": left["last_slot"],
                "right_first_slot": right["first_slot"],
                "right_last_slot": right["last_slot"],
                "left_center_Omega": left["center_Omega"],
                "right_center_Omega": right["center_Omega"],
                "relative_difference": relative,
                "tolerance": tolerance,
                "left_multiplicity": left["multiplicity"],
                "right_multiplicity": right["multiplicity"],
                "left_total_nullity": left["total_nullity"],
                "right_total_nullity": right["total_nullity"],
                "status": "PASS" if row_pass else "FAIL",
            }
        )
    return records, {
        "status": "PASS" if passed else "FAIL",
        "maximum_relative_difference": maximum,
        "left_group_count": len(left_groups),
        "right_group_count": len(right_groups),
    }


def _direct_reference_rows(
    alpha_A: float,
    union_rows: Sequence[Mapping[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    started = time.perf_counter()
    _beam, coupled = _physics_modules()
    section = build_layered_section(alpha_A)

    def provider(omega: float) -> FloatArray:
        return np.asarray(
            coupled.direct_fixed_fixed_boundary_matrix(
                float(omega), L_TOTAL, section.properties
            ),
            dtype=float,
        )

    policy = _policy_root9()
    evaluator = beta0_tools._DiagnosticEvaluator(
        provider,
        OMEGA_TO_OMEGA_SCALE,
        policy,
    )
    rows: list[dict[str, Any]] = []
    quality_failures = 0
    for row in sorted(union_rows, key=lambda item: int(item["sorted_position"])):
        diagnostic = evaluator.diagnostics(float(row["Omega"]))
        quality_ok, reason = beta0_tools._candidate_quality(diagnostic, policy)
        if not quality_ok:
            quality_failures += 1
        rows.append(
            {
                "sorted_position": int(row["sorted_position"]),
                "Omega": float(row["Omega"]),
                "cluster_id": str(row["cluster_id"]),
                "cluster_semantics": str(row["cluster_semantics"]),
                "cluster_multiplicity": int(row["cluster_multiplicity"]),
                "cluster_total_nullity": int(diagnostic.detected_nullity),
                "cluster_center_Omega": float(row["cluster_center_Omega"]),
                "cluster_metadata_source": (
                    "UNION_SLOT_GROUPING_WITH_DIRECT_FIXED_FIXED_DETECTED_NULLITY"
                ),
                "direct_detected_nullity": int(diagnostic.detected_nullity),
                "scaled_sigma_ratio": float(diagnostic.scaled_sigma_ratio),
                "boundary_null_residual": float(
                    diagnostic.raw_boundary_null_residual
                ),
                "quality_status": "PASS" if quality_ok else reason,
            }
        )
    return rows, {
        "alpha_A": alpha_A,
        "wall_time_seconds": time.perf_counter() - started,
        "peak_rss_bytes": rlb2e._peak_rss_bytes(),
        "determinant_evaluations": len(rows),
        "sigma_evaluations": len(rows),
        "maximum_diagnostic_Omega": max(float(row["Omega"]) for row in union_rows),
        "local_refinements": 0,
        "diagnostic_frequency_count": len(rows),
        "root_count_localized_by_direct_check": 0,
        "root9_right_margin_Omega": "NOT_APPLICABLE_NO_DIRECT_SCAN",
        "maximum_scaled_sigma_ratio": max(
            float(row["scaled_sigma_ratio"]) for row in rows
        ),
        "maximum_boundary_residual": max(
            float(row["boundary_null_residual"]) for row in rows
        ),
        "unresolved_candidates_below_root9": "NOT_APPLICABLE_NO_DIRECT_SCAN",
        "roots_above_9_computed": "NOT_APPLICABLE_NO_DIRECT_SCAN",
        "frequencies_supplied_by": "EXACT_AXIAL_PLUS_SINGLE_BENDING_REFERENCE_UNION",
        "used_for_root_localization": False,
        "scientific_role": "DIRECT_FIXED_FIXED_BOUNDARY_SINGULARITY_CHECK",
        "quality_failure_count": quality_failures,
        "detected_nullities": [
            int(row["direct_detected_nullity"]) for row in rows
        ],
        "status": "PASS" if quality_failures == 0 else "FAIL",
    }


def coupled_beta0_direct_checks(
    rows: Sequence[Mapping[str, Any]],
    reference_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    comparisons: list[dict[str, Any]] = []
    direct_records: list[dict[str, Any]] = []
    maxima = {
        "coupled_vs_direct": 0.0,
        "coupled_vs_subsystem_union": 0.0,
        "direct_vs_subsystem_union": 0.0,
    }
    passed = True
    for alpha_A in (0.70, 1.00, 1.30):
        coupled_rows = _canonical_group(rows, 0.0, alpha_A)
        union_rows = [
            row
            for row in reference_rows
            if abs(float(row["alpha_A"]) - alpha_A) < 1.0e-12
        ]
        direct_rows, direct_record = _direct_reference_rows(alpha_A, union_rows)
        direct_records.append(direct_record)
        for left, right, kind, key in (
            (
                coupled_rows,
                direct_rows,
                "coupled_vs_direct_fixed_fixed",
                "coupled_vs_direct",
            ),
            (
                coupled_rows,
                union_rows,
                "coupled_vs_axial_bending_union",
                "coupled_vs_subsystem_union",
            ),
            (
                direct_rows,
                union_rows,
                "direct_fixed_fixed_vs_axial_bending_union",
                "direct_vs_subsystem_union",
            ),
        ):
            comparison_rows, summary = _compare_cluster_groups(
                left,
                right,
                alpha_A=alpha_A,
                comparison_kind=kind,
            )
            comparisons.extend(comparison_rows)
            maxima[key] = max(
                maxima[key], float(summary["maximum_relative_difference"])
            )
            passed = passed and summary["status"] == "PASS"
        passed = passed and direct_record["status"] == "PASS"
    return {
        "status": "PASS" if passed else "FAIL",
        "maximum_relative_difference": maxima["coupled_vs_direct"],
        "maximum_coupled_vs_subsystem_union_relative_difference": maxima[
            "coupled_vs_subsystem_union"
        ],
        "maximum_direct_vs_subsystem_union_relative_difference": maxima[
            "direct_vs_subsystem_union"
        ],
        "comparisons": comparisons,
        "direct_reference_records": direct_records,
        "coupled_roots_recomputed": False,
        "direct_reference_check_count": len(direct_records),
        "direct_reference_solve_count": 0,
        "direct_reference_frequencies_localized_independently": False,
        "direct_reference_frequencies_supplied_by": (
            "EXACT_AXIAL_PLUS_SINGLE_ALPHA1_BENDING_REFERENCE_UNION"
        ),
    }


def beta0_all_grid_subsystem_check(
    rows: Sequence[Mapping[str, Any]],
    reference_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Compare every beta=0 grid point with the decoupled family union."""

    records: list[dict[str, Any]] = []
    maximum = 0.0
    passed = True
    for alpha_A in alpha_grid():
        value = float(alpha_A)
        coupled_rows = _canonical_group(rows, 0.0, value)
        union_rows = [
            row
            for row in reference_rows
            if abs(float(row["alpha_A"]) - value) < 1.0e-12
        ]
        _comparisons, summary = _compare_cluster_groups(
            coupled_rows,
            union_rows,
            alpha_A=value,
            comparison_kind="coupled_vs_axial_bending_union_all_grid",
        )
        maximum = max(maximum, float(summary["maximum_relative_difference"]))
        passed = passed and summary["status"] == "PASS"
        records.append(
            {
                "alpha_A": value,
                "status": summary["status"],
                "maximum_relative_difference": summary[
                    "maximum_relative_difference"
                ],
                "coupled_group_count": summary["left_group_count"],
                "reference_group_count": summary["right_group_count"],
            }
        )
    return {
        "status": "PASS" if passed else "FAIL",
        "point_count": len(records),
        "maximum_relative_difference": maximum,
        "point_records": records,
    }


def apply_beta0_reference_classification(
    rows: list[dict[str, Any]], reference_rows: Sequence[Mapping[str, Any]]
) -> None:
    reference = {
        (
            round(float(row["alpha_A"]), 10),
            int(row["sorted_position"]),
        ): row
        for row in reference_rows
    }
    for row in rows:
        if abs(float(row["beta_deg"])) > 1.0e-12:
            continue
        item = reference.get(
            (
                round(float(row["alpha_A"]), 10),
                int(row["sorted_position"]),
            )
        )
        if item is None:
            row["decoupled_reference_classification"] = "UNCLASSIFIED"
            continue
        relative = _relative(float(row["Omega"]), float(item["Omega"]))
        row["decoupled_reference_classification"] = (
            str(item["classification"])
            if relative <= DIRECT_REFERENCE_RELATIVE_TOLERANCE
            else "UNCLASSIFIED"
        )
        if relative <= DIRECT_REFERENCE_RELATIVE_TOLERANCE:
            previous_cluster = {
                field: row.get(field, "")
                for field in (
                    "cluster_id",
                    "cluster_semantics",
                    "cluster_multiplicity",
                    "cluster_total_nullity",
                    "cluster_center_Omega",
                )
            }
            for field in (
                "cluster_id",
                "cluster_semantics",
                "cluster_multiplicity",
                "cluster_total_nullity",
                "cluster_center_Omega",
            ):
                row[field] = item[field]
            changed = any(
                str(previous_cluster[field]) != str(item[field])
                for field in previous_cluster
            )
            row["cluster_metadata_source"] = (
                "BETA0_EXACT_SUBSYSTEM_REFERENCE_OVERRIDE"
                if changed
                else "BETA0_EXACT_SUBSYSTEM_REFERENCE_CONFIRMED"
            )


def backfill_resumed_spectrum_metadata(
    rows: list[dict[str, Any]],
) -> dict[str, int]:
    """Fill additive resume-schema fields without changing any root value."""

    source_backfills = 0
    guard_backfills = 0
    for row in rows:
        if not str(row.get("cluster_metadata_source", "")).strip():
            row["cluster_metadata_source"] = "DETERMINANT_SVD_EVENT_RESUMED"
            source_backfills += 1
    groups: dict[tuple[float, float, str, str], list[dict[str, Any]]] = {}
    for row in rows:
        key = (
            round(float(row["beta_deg"]), 10),
            round(float(row["alpha_A"]), 10),
            str(row["grid_kind"]),
            str(row.get("repair_id", "")),
        )
        groups.setdefault(key, []).append(row)
    for group in groups.values():
        guards = [row for row in group if int(row["sorted_position"]) == K_GUARD]
        if len(guards) != 1:
            continue
        guard = guards[0]
        cluster_id = str(guard.get("cluster_id", ""))
        multiplicity = int(guard.get("cluster_multiplicity") or 1)
        represented = sum(
            str(row.get("cluster_id", "")) == cluster_id for row in group
        )
        if not str(guard.get("guard_cluster_multiplicity", "")).strip():
            guard["guard_cluster_multiplicity"] = multiplicity
            guard_backfills += 1
        guard["guard_cluster_extends_beyond_export"] = multiplicity > represented
    return {
        "cluster_metadata_source_backfilled_rows": source_backfills,
        "guard_metadata_backfilled_groups": guard_backfills,
        "root_values_modified": 0,
    }


def backfill_resumed_point_record_metadata(
    point_records: Sequence[dict[str, Any]],
    rows: Sequence[Mapping[str, Any]],
) -> int:
    """Synchronize additive root-9 provenance in durable point records."""

    groups = _base_group_index(rows)
    updated = 0
    for record in point_records:
        key = _group_key(float(record["beta_deg"]), float(record["alpha_A"]))
        guard = next(
            (
                row
                for row in groups.get(key, [])
                if int(row["sorted_position"]) == K_GUARD
            ),
            None,
        )
        if guard is None:
            continue
        before = (
            record.get("guard_cluster_multiplicity"),
            record.get("guard_cluster_extends_beyond_export"),
        )
        record["guard_cluster_multiplicity"] = int(
            guard["guard_cluster_multiplicity"]
        )
        record["guard_cluster_extends_beyond_export"] = _as_bool(
            guard["guard_cluster_extends_beyond_export"]
        )
        after = (
            record["guard_cluster_multiplicity"],
            record["guard_cluster_extends_beyond_export"],
        )
        updated += before != after
    return updated


def create_plot_from_csv(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    started = time.perf_counter()
    rows = rlb2e._read_csv(Path(output_dir) / SPECTRUM_FILENAME)
    spectrum_audit = audit_spectrum_rows(rows)
    if spectrum_audit["status"] != "PASS":
        raise RuntimeError(f"plot_only rejected root inventory: {spectrum_audit}")
    plot_audit = audit_plot_rows(rows)
    if plot_audit["status"] != "PASS":
        raise RuntimeError(f"plot_only rejected plot data: {plot_audit}")
    selected_rows = canonical_plot_rows(rows)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = plt.get_cmap("tab10").colors[:K_PLOT]
    figure, axes = plt.subplots(1, 2, figsize=(12.8, 5.0), sharex=True, sharey=True)
    titles = ("(a) beta = 0 deg", "(b) beta = 30 deg")
    for axis, beta_deg, title in zip(axes, BETA_VALUES_DEG, titles, strict=True):
        for position in range(1, K_PLOT + 1):
            panel_rows = [
                row
                for row in selected_rows
                if round(float(row["beta_deg"]), 10) == round(beta_deg, 10)
                and int(row["sorted_position"]) == position
            ]
            panel_rows.sort(key=lambda row: float(row["alpha_A"]))
            axis.plot(
                [float(row["alpha_A"]) for row in panel_rows],
                [float(row["Lambda"]) for row in panel_rows],
                color=colors[position - 1],
                linestyle="-",
                linewidth=1.25,
                label=f"k={position}",
            )
        axis.axvline(1.0, color="0.55", linewidth=0.8, linestyle=":")
        axis.set_title(title)
        axis.set_xlabel(r"$A/A_0$")
        axis.grid(True, alpha=0.22, linewidth=0.5)
    axes[0].set_ylabel(r"$\Lambda$")
    figure.legend(
        handles=axes[0].lines[:K_PLOT],
        labels=[f"k={index}" for index in range(1, K_PLOT + 1)],
        loc="upper center",
        ncol=8,
        frameon=False,
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.91))
    target = Path(output_dir) / PLOT_FILENAME
    temporary = target.with_name(target.stem + ".tmp" + target.suffix)
    figure.savefig(temporary, dpi=180, bbox_inches="tight")
    plt.close(figure)
    os.replace(temporary, target)
    return {
        "path": target.as_posix(),
        "wall_time_seconds": time.perf_counter() - started,
        "panel_count": 2,
        "lines_per_panel": K_PLOT,
        "plotted_positions": list(range(1, K_PLOT + 1)),
        "root9_plotted": False,
        "root_calculation_count": 0,
        "spectrum_data_audit": spectrum_audit,
        "plot_data_audit": plot_audit,
    }


def create_beta0_reference_plot(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    started = time.perf_counter()
    reference_path = Path(output_dir) / REFERENCE_FILENAME
    rows = rlb2e._read_csv(reference_path)
    if not rows:
        raise RuntimeError("beta0_subsystem_reference.csv is missing or empty.")
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(7.0, 4.8))
    alpha_values = sorted({float(row["alpha_A"]) for row in rows})
    axial_groups: dict[int, list[tuple[float, float]]] = {}
    bending_groups: dict[int, list[tuple[float, float]]] = {}
    for row in rows:
        family_index = int(row["family_index"])
        pair = (float(row["alpha_A"]), float(row["Lambda"]))
        if str(row["family"]) == "axial":
            axial_groups.setdefault(family_index, []).append(pair)
        elif str(row["family"]) == "bending":
            bending_groups.setdefault(family_index, []).append(pair)
    for family_index, data in sorted(bending_groups.items()):
        data.sort()
        axis.plot(
            [item[0] for item in data],
            [item[1] for item in data],
            color="0.35",
            linestyle="--",
            linewidth=1.0,
            label="bending reference" if family_index == 1 else None,
        )
    for family_index, data in sorted(axial_groups.items()):
        data.sort()
        axis.plot(
            [item[0] for item in data],
            [item[1] for item in data],
            color="tab:blue",
            linestyle="-.",
            linewidth=1.0,
            label="axial exact reference" if family_index == 1 else None,
        )
    axis.set_xlabel(r"$A/A_0$")
    axis.set_ylabel(r"$\Lambda$")
    axis.axvline(1.0, color="0.55", linewidth=0.8, linestyle=":")
    axis.grid(True, alpha=0.22, linewidth=0.5)
    axis.legend(frameon=False)
    target = Path(output_dir) / BETA0_PLOT_FILENAME
    temporary = target.with_name(target.stem + ".tmp" + target.suffix)
    figure.tight_layout()
    figure.savefig(temporary, dpi=180, bbox_inches="tight")
    plt.close(figure)
    os.replace(temporary, target)
    return {
        "path": target.as_posix(),
        "wall_time_seconds": time.perf_counter() - started,
        "root_calculation_count": 0,
    }


def _minimum_adjacent_gaps(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for beta_deg in BETA_VALUES_DEG:
        best = (math.inf, math.nan, -1)
        for alpha_A in alpha_grid():
            roots = _rows_for_roots(rows, beta_deg, float(alpha_A))
            gaps = np.diff(np.sqrt(roots[:K_PLOT]))
            index = int(np.argmin(gaps))
            if float(gaps[index]) < best[0]:
                best = (float(gaps[index]), float(alpha_A), index + 1)
        records.append(
            {
                "beta_deg": beta_deg,
                "minimum_adjacent_Lambda_gap": best[0],
                "alpha_A": best[1],
                "between_sorted_positions": [best[2], best[2] + 1],
                "interpretation": "candidate interval only; no branch or shape claim",
            }
        )
    return records


def _beta30_shift_summary(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    summary: list[dict[str, Any]] = []
    for position in range(1, K_PLOT + 1):
        values = []
        minimum_gap = math.inf
        minimum_gap_alpha = math.nan
        minimum_gap_neighbour = -1
        for alpha_A in alpha_grid():
            ordered = _canonical_group(rows, 30.0, float(alpha_A))
            values.append(float(ordered[position - 1]["Lambda"]))
            neighbours = []
            if position > 1:
                neighbours.append(position - 1)
            if position < K_GUARD:
                neighbours.append(position + 1)
            for neighbour in neighbours:
                gap = abs(
                    float(ordered[position - 1]["Lambda"])
                    - float(ordered[neighbour - 1]["Lambda"])
                )
                if gap < minimum_gap:
                    minimum_gap = gap
                    minimum_gap_alpha = float(alpha_A)
                    minimum_gap_neighbour = neighbour
        low, mid, high = values[0], values[15], values[-1]
        changes = np.diff(values)
        if np.all(changes >= -1.0e-12):
            monotonicity = "NONDECREASING"
        elif np.all(changes <= 1.0e-12):
            monotonicity = "NONINCREASING"
        else:
            monotonicity = "NONMONOTONE"
        summary.append(
            {
                "sorted_position": position,
                "Lambda_alpha_0p70": low,
                "Lambda_alpha_1p00": mid,
                "Lambda_alpha_1p30": high,
                "relative_change_0p70_to_1p30": (high - low)
                / max(abs(low), np.finfo(float).tiny),
                "maximum_relative_deviation_from_alpha_1": max(
                    abs(value - mid) / max(abs(mid), np.finfo(float).tiny) for value in values
                ),
                "monotonicity_class": monotonicity,
                "minimum_adjacent_Lambda_gap": minimum_gap,
                "minimum_gap_alpha_A": minimum_gap_alpha,
                "minimum_gap_neighbour_position": minimum_gap_neighbour,
                "minimum_gap_neighbour_role": (
                    "ROOT_9_GUARD"
                    if minimum_gap_neighbour == K_GUARD
                    else "PLOTTED"
                ),
            }
        )
    return summary


def _root_quality_summary(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    base = [row for row in rows if str(row["grid_kind"]) == "BASE"]
    guards = [row for row in base if int(row["sorted_position"]) == K_GUARD]
    canonical = [
        row
        for row in rows
        if _as_bool(row.get("is_canonical_plot_source", False))
        and int(row["sorted_position"]) <= K_GUARD
    ]
    unresolved_canonical = [
        row
        for row in canonical
        if str(row.get("point_status", "")) == "UNRESOLVED_AFTER_LOCAL_REPAIR"
    ]
    resolved_canonical = [row for row in canonical if row not in unresolved_canonical]
    canonical_quality_failures = [
        str(row.get("row_id", ""))
        for row in resolved_canonical
        if not (
            str(row.get("quality_status", "")) == "PASS"
            and math.isfinite(float(row["Omega"]))
            and float(row["Omega"]) > 0.0
            and float(row["scaled_sigma_ratio"]) <= ROOT_SINGULAR_RATIO_TOLERANCE
            and float(row["boundary_null_residual"]) <= BOUNDARY_RESIDUAL_TOLERANCE
            and int(row["unresolved_candidates_below_root9"]) == 0
            and not _as_bool(row["predictor_used_as_final"])
        )
    ]
    return {
        "status": (
            "FAIL"
            if canonical_quality_failures
            else ("PARTIAL_PASS" if unresolved_canonical else "PASS")
        ),
        "base_quality_failures": sum(str(row["quality_status"]) != "PASS" for row in base),
        "maximum_base_scaled_sigma_ratio": max(float(row["scaled_sigma_ratio"]) for row in base),
        "maximum_base_boundary_null_residual": max(
            float(row["boundary_null_residual"]) for row in base
        ),
        "maximum_root9_scaled_sigma_ratio": max(float(row["scaled_sigma_ratio"]) for row in guards),
        "maximum_root9_boundary_null_residual": max(
            float(row["boundary_null_residual"]) for row in guards
        ),
        "minimum_root9_right_margin_Omega": min(
            float(row["root9_right_margin_Omega"]) for row in guards
        ),
        "maximum_unresolved_candidates_below_root9": max(
            int(row["unresolved_candidates_below_root9"]) for row in base
        ),
        "canonical_row_count": len(canonical),
        "canonical_resolved_row_count": len(resolved_canonical),
        "canonical_unresolved_row_count": len(unresolved_canonical),
        "canonical_quality_failure_row_ids": canonical_quality_failures,
        "maximum_canonical_scaled_sigma_ratio": max(
            (float(row["scaled_sigma_ratio"]) for row in resolved_canonical),
            default=math.nan,
        ),
        "maximum_canonical_boundary_null_residual": max(
            (float(row["boundary_null_residual"]) for row in resolved_canonical),
            default=math.nan,
        ),
    }


def _section_property_range_summary(output_dir: Path) -> dict[str, Any]:
    rows = rlb2e._read_csv(Path(output_dir) / SECTION_FILENAME)
    if len(rows) != len(alpha_grid()):
        raise RuntimeError("section_properties.csv must contain exactly 31 rows.")
    fields = (
        "A_beam",
        "A_beam_over_A0",
        "D_beam",
        "D_beam_over_D0",
        "S_beam",
        "S_beam_over_S0",
        "m",
        "m_over_m0",
        "J",
        "J_over_J0",
        "s_outer",
        "s_inner",
    )
    summary: dict[str, Any] = {"row_count": len(rows)}
    for field in fields:
        values = [float(row[field]) for row in rows]
        summary[field] = {"minimum": min(values), "maximum": max(values)}
    return summary


def _recover_existing_repair_records(
    rows: Sequence[Mapping[str, Any]],
    audit_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    groups: dict[tuple[float, float, str], list[Mapping[str, Any]]] = {}
    for row in rows:
        if str(row.get("grid_kind")) != "LOCAL_REFINEMENT":
            continue
        repair_id = str(row.get("repair_id", ""))
        key = (float(row["beta_deg"]), float(row["alpha_A"]), repair_id)
        groups.setdefault(key, []).append(row)
    records: list[dict[str, Any]] = []
    for (beta_deg, alpha_A, repair_id), group in sorted(groups.items()):
        positions = sorted(int(row["sorted_position"]) for row in group)
        if not repair_id or len(positions) != len(set(positions)):
            raise RuntimeError("Stored LOCAL_REFINEMENT transaction is malformed.")
        statuses = {str(row.get("point_status", "")) for row in group}
        if len(statuses) != 1:
            raise RuntimeError("Stored LOCAL_REFINEMENT transaction has mixed statuses.")
        stored_status = statuses.pop()
        unresolved = stored_status == "UNRESOLVED_AFTER_LOCAL_REPAIR"
        if not unresolved and positions != list(range(1, K_GUARD + 1)):
            raise RuntimeError("Stored resolved LOCAL_REFINEMENT transaction is incomplete.")
        affected = sorted(
            {
                int(row["sorted_position"])
                for row in audit_rows
                if str(row.get("repair_id", "")) == repair_id
                and _group_key(float(row["beta_deg"]), float(row["alpha_A"]))
                == _group_key(beta_deg, alpha_A)
            }
        )
        if not affected:
            raise RuntimeError("Stored local repair has no matching neighbour-audit evidence.")
        if unresolved and positions != affected:
            raise RuntimeError(
                "Stored unresolved LOCAL_REFINEMENT gaps do not match the "
                "flagged neighbour-audit positions."
            )
        records.append(
            {
                "repair_id": repair_id,
                "beta_deg": beta_deg,
                "alpha_A": alpha_A,
                "status": "UNRESOLVED" if unresolved else stored_status,
                "affected_positions": affected,
                "wall_time_seconds": None,
                "peak_rss_bytes": None,
                "determinant_evaluations": None,
                "sigma_evaluations": None,
                "metrics_provenance": (
                    "RECOVERED_AFTER_INTERRUPTED_PREFINALIZATION; numerical repair "
                    "metrics were not persisted"
                ),
                "smoothing_applied": False,
                "predictor_used_as_final": False,
            }
        )
    return records


def _runtime_summary(
    point_records: Sequence[Mapping[str, Any]],
    repair_records: Sequence[Mapping[str, Any]],
    direct_records: Sequence[Mapping[str, Any]],
    beta0_reference: Mapping[str, Any],
    plot: Mapping[str, Any],
    beta0_plot: Mapping[str, Any],
    *,
    finalization_process_seconds: float,
) -> dict[str, Any]:
    base_seconds = sum(float(item.get("wall_time_seconds", 0.0)) for item in point_records)
    repair_seconds = sum(
        float(item.get("wall_time_seconds") or 0.0) for item in repair_records
    )
    direct_seconds = sum(
        float(item.get("wall_time_seconds") or 0.0) for item in direct_records
    )
    reference_seconds = float(beta0_reference.get("wall_time_seconds") or 0.0)
    plot_seconds = float(plot.get("wall_time_seconds", 0.0)) + float(
        beta0_plot.get("wall_time_seconds", 0.0)
    )
    missing_repair_metrics = sum(
        item.get("determinant_evaluations") is None for item in repair_records
    )
    determinant_evaluations = sum(
        int(item.get("determinant_evaluations", 0)) for item in point_records
    ) + sum(int(item.get("determinant_evaluations") or 0) for item in repair_records)
    determinant_evaluations += sum(
        int(item.get("determinant_evaluations") or 0) for item in direct_records
    )
    sigma_evaluations = sum(
        int(item.get("sigma_evaluations", 0)) for item in point_records
    ) + sum(int(item.get("sigma_evaluations") or 0) for item in repair_records)
    sigma_evaluations += sum(
        int(item.get("sigma_evaluations") or 0) for item in direct_records
    )
    return {
        "base_point_wall_time_sum_seconds": base_seconds,
        "local_repair_wall_time_sum_seconds": repair_seconds,
        "direct_reference_wall_time_sum_seconds": direct_seconds,
        "beta0_bending_reference_wall_time_seconds": reference_seconds,
        "figure_rendering_seconds": plot_seconds,
        "total_measured_workflow_seconds_lower_bound": (
            base_seconds
            + repair_seconds
            + direct_seconds
            + reference_seconds
            + plot_seconds
        ),
        "total_measured_workflow_seconds": (
            base_seconds
            + repair_seconds
            + direct_seconds
            + reference_seconds
            + plot_seconds
        ),
        "timing_is_complete": missing_repair_metrics == 0,
        "timing_qualification": (
            "COMPLETE"
            if missing_repair_metrics == 0
            else "LOWER_BOUND; LOCAL_REPAIR_TIMES_NOT_PERSISTED_BEFORE_INTERRUPTED_FINALIZATION"
        ),
        "finalization_process_seconds": float(finalization_process_seconds),
        "peak_rss_bytes": max(
            rlb2e._peak_rss_bytes(),
            max(
                [int(item.get("peak_rss_bytes") or 0) for item in point_records]
                + [int(item.get("peak_rss_bytes") or 0) for item in direct_records]
                + [int(item.get("peak_rss_bytes") or 0) for item in repair_records]
                + [0]
            ),
        ),
        "determinant_evaluations": determinant_evaluations,
        "sigma_evaluations": sigma_evaluations,
        "base_root_solve_count": sum(
            int(item.get("determinant_evaluations", 0)) > 0 for item in point_records
        ),
        "local_repair_solve_count": sum(
            int(item.get("determinant_evaluations") or 0) > 0 for item in repair_records
        ),
        "direct_reference_solve_count": 0,
        "direct_boundary_diagnostic_point_count": len(direct_records),
        "beta0_bending_reference_root_solve_count": int(
            beta0_reference.get("root_solver_invocations", 0)
        ),
        "total_root_solve_count": len(point_records)
        + len(repair_records)
        + int(beta0_reference.get("root_solver_invocations", 0)),
        "root_calculations_in_current_closing_invocation": ROOT_CALCULATION_COUNT,
        "determinant_evaluation_count_is_lower_bound": missing_repair_metrics > 0,
        "sigma_evaluation_count_is_lower_bound": missing_repair_metrics > 0,
        "local_repair_metrics_missing_after_interrupted_prefinalization": missing_repair_metrics,
        "parallel_spectral_workers": 0,
    }


def _sha_tree(path: Path) -> str | None:
    if not path.exists():
        return None
    digest = hashlib.sha256()
    for item in sorted(candidate for candidate in path.rglob("*") if candidate.is_file()):
        relative = item.relative_to(path).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(bytes.fromhex(rlb2e._sha256(item)))
    return digest.hexdigest().upper()


def contract_payload() -> dict[str, Any]:
    lower, upper = continuation_paths()
    return {
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "frequency_map_policy_instances": {
            f"beta_{int(beta_deg)}": _contract_instance(beta_deg) for beta_deg in BETA_VALUES_DEG
        },
        "geometry": {
            "mu": MU,
            "tau": TAU,
            "beta_values_deg": list(BETA_VALUES_DEG),
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
            "all_ply_angles_deg": 0.0,
        },
        "material_M0": base_material_contract(),
        "alpha_A_definition": {
            "range": [ALPHA_A_MIN, ALPHA_A_MAX],
            "step": ALPHA_A_STEP,
            "anchor": ALPHA_A_ANCHOR,
            "s_outer": "(4-alpha_A)/3",
            "s_inner": "(7*alpha_A-4)/3",
            "in_plane_scaled": ["E1", "E2", "G12"],
            "in_plane_fixed": ["nu12"],
            "transverse_shear_fixed": ["G13", "G23"],
            "density_fixed": True,
            "expected_only_changed_reduced_field": "A",
        },
        "stack_bottom_to_top": list(LAYOUT),
        "alpha_A_grid": [float(value) for value in alpha_grid()],
        "continuation": {
            "anchor": ALPHA_A_ANCHOR,
            "lower_leg": [float(value) for value in lower],
            "upper_leg": [float(value) for value in upper],
        },
        "normalization": {
            "reference_area": REFERENCE_AREA,
            "I_y0": IY_REFERENCE,
            "Omega": "omega*l^2*sqrt(rho0*reference_area/(E0*I_y0))",
            "Lambda": "sqrt(Omega)",
            "E0": 1.0,
            "rho0": RHO_0,
            "b0": WIDTH,
            "h0": THICKNESS,
            "l": L_REFERENCE,
            "Omega_per_omega": OMEGA_TO_OMEGA_SCALE,
        },
        "root_contract": {
            "plotted_positions": list(range(1, K_PLOT + 1)),
            "guard_position": K_GUARD,
            "root9_role": "completeness_only",
            "search_policy_requested_roots": K_PLOT,
            "search_policy_guard_roots": 1,
            "root9_plotted": False,
            "roots_above_9_computed": False,
            "branch_tracking": False,
        },
        "thresholds": {
            "matrix_relative": MATRIX_RELATIVE_TOLERANCE,
            "symmetry_relative": SYMMETRY_RELATIVE_TOLERANCE,
            "reduced_property_relative": REDUCED_PROPERTY_TOLERANCE,
            "reduction_route_relative": REDUCTION_ROUTE_TOLERANCE,
            "root_singular_ratio": ROOT_SINGULAR_RATIO_TOLERANCE,
            "boundary_residual": BOUNDARY_RESIDUAL_TOLERANCE,
            "beta0_direct_reference_relative": DIRECT_REFERENCE_RELATIVE_TOLERANCE,
            "cluster_reference_relative": CLUSTER_REFERENCE_RELATIVE_TOLERANCE,
        },
        "exclusions": {
            "new_physical_model": False,
            "governing_equation_changes": False,
            "joint_matrix_changes": False,
            "beta_values_beyond_declared": False,
            "mu_tau_sweep": False,
            "roots_10_plus": False,
            "branch_tracking": False,
            "MAC": False,
            "mode_shapes": False,
            "energy_analysis": False,
            "Ritz": False,
            "FEM": False,
            "certified_audit": False,
            "universal_runner": False,
        },
    }


def contract_hash() -> str:
    return hashlib.sha256(
        json.dumps(
            rlb2e._json_value(contract_payload()),
            sort_keys=True,
            ensure_ascii=False,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest().upper()


def _output_hashes(output_dir: Path) -> dict[str, str]:
    return {
        path.name: rlb2e._sha256(path)
        for path in sorted(Path(output_dir).iterdir())
        if path.is_file() and path.name != MANIFEST_FILENAME
    }


def _report_text(manifest: Mapping[str, Any]) -> str:
    counts = manifest["counts"]
    constitutive = manifest["constitutive_gate"]
    reference = manifest["beta0_reference"]
    benchmark = manifest["benchmark"]
    audit = manifest["neighbour_audit"]
    quality = manifest["root_quality_summary"]
    runtime = manifest["runtime"]
    beta30 = manifest["beta30_shift_summary"]
    gap_lines = "\n".join(
        "- `beta={beta:.0f}`: $\\Delta\\Lambda_{{min}}={gap:.6g}$ при "
        "$A/A_0={alpha:.2f}$ между позициями {left} и {right}.".format(
            beta=float(item["beta_deg"]),
            gap=float(item["minimum_adjacent_Lambda_gap"]),
            alpha=float(item["alpha_A"]),
            left=int(item["between_sorted_positions"][0]),
            right=int(item["between_sorted_positions"][1]),
        )
        for item in manifest["minimum_adjacent_sorted_gaps"]
    )
    beta30_lines = "\n".join(
        "- `k={k}`: $\\Lambda(0.70)={lo:.6f}$, $\\Lambda(1.00)={mid:.6f}$, "
        "$\\Lambda(1.30)={hi:.6f}$, относительное изменение {chg:+.6e}, "
        "класс `{cls}`.".format(
            k=int(item["sorted_position"]),
            lo=float(item["Lambda_alpha_0p70"]),
            mid=float(item["Lambda_alpha_1p00"]),
            hi=float(item["Lambda_alpha_1p30"]),
            chg=float(item["relative_change_0p70_to_1p30"]),
            cls=str(item["monotonicity_class"]),
        )
        for item in beta30
    )
    return f"""# RLB-2H: axial-stiffness visibility under geometric coupling

## Scope

Этап строит ordinary frequency map первых восьми independently sorted
частот при изменении только reduced axial stiffness $A_{{beam}}$. Два
конечных случая: $\\beta=0^\\circ$ и $\\beta=30^\\circ$.

## Frozen contract

Использованы $\\mu=\\tau=0$, $L_1=L_2=1$, $L_{{total}}=2$,
$b=0.20$, $h=0.05$, $K=5/6$ и четыре равных 0°-слоя толщины 0.0125.

Базовый материал $M_0$:
$E_1=1.1$, $E_2=0.9$, $\\nu_{{12}}=0.3$,
$G_{{12}}=G_{{13}}=G_{{23}}=1/2.6$, $\\rho=1$.

Sweep parameter:
$\\alpha_A=A/A_0=0.70,0.72,\\ldots,1.30$.
Множители слоёв:
$s_{{outer}}=(4-\\alpha_A)/3$,
$s_{{inner}}=(7\\alpha_A-4)/3$.

## Constitutive gate

Статус: **{constitutive['status']}**.
Подтверждено, что меняется только $A_{{beam}}$.
Максимальные residuals:

- $A$-matrix scale: {constitutive['maximum_residuals']['A_matrix_scale_residual']:.3e}
- $D$-matrix invariance: {constitutive['maximum_residuals']['D_matrix_invariance_residual']:.3e}
- $A/A_0$ formula: {constitutive['maximum_residuals']['A_formula_residual']:.3e}
- $D$, $S$, $m$, $J$ invariants: {constitutive['maximum_residuals']['only_A_field_max_invariance_residual']:.3e}
- $B$: {constitutive['maximum_residuals']['B_relative']:.3e}
- $I_1$: {constitutive['maximum_residuals']['I1_relative']:.3e}

## beta = 0 reference

Direct fixed-fixed and subsystem reference status:
**{reference['status']}**.

- Bending invariance error: {reference['bending_invariance_error']:.3e}
- Axial $\\Lambda\\propto\\alpha_A^{{1/4}}$ error:
  {reference['axial_scaling_error']:.3e}
- Coupled/direct fixed-fixed maximum relative difference:
  {manifest['beta0_direct_coupled_check']['maximum_relative_difference']:.3e}

## Frequency map

Получено {counts['beta0_points_complete']}/31 point(s) при $\\beta=0^\\circ$
и {counts['beta30_points_complete']}/31 point(s) при $\\beta=30^\\circ$.
BASE rows: {counts['base_rows']}/558. Root-9 guards:
{counts['root9_guards']}/62.

Neighbour flags: {audit['flagged_point_count']}; local repairs:
{audit['repair_count']}; unresolved: {audit['unresolved_point_count']}.

Minimum adjacent sorted gaps:

{gap_lines}

### beta = 30 endpoint response

{beta30_lines}

## Runtime

Benchmark ETA:
{benchmark['conservative_eta_seconds']:.1f} s for
{benchmark['remaining_unique_root_points']} remaining point(s) after six
production anchors.

Measured workflow time:
{runtime['total_measured_workflow_seconds']:.1f} s.
Peak RSS:
{runtime['peak_rss_bytes'] / 2**20:.1f} MiB.
Determinant evaluations:
{runtime['determinant_evaluations']}.
Sigma/SVD evaluations:
{runtime['sigma_evaluations']}.

Maximum base $\\sigma_{{min}}/\\sigma_{{max}}$:
{quality['maximum_base_scaled_sigma_ratio']:.3e}.
Maximum base boundary residual:
{quality['maximum_base_boundary_null_residual']:.3e}.
Minimum root-9 right margin:
{quality['minimum_root9_right_margin_Omega']:.3g}.

## Outputs and status

Main plot:
`{manifest['plot']['path']}`.

Scientific status:
**{manifest['scientific_status']}**.

This result concerns only the declared finite grid, independently sorted
positions 1--8, and root 9 as completeness guard. No branch identity,
crossing, veering, localization, MAC, mode-shape, Ritz, FEM, smoothing, or
certified-audit claim is made.
"""


def _report_text_v2(manifest: Mapping[str, Any]) -> str:
    """Render the closing scientific report without modal-identity claims."""

    counts = manifest["counts"]
    constitutive = manifest["constitutive_gate"]
    residuals = constitutive["maximum_residuals"]
    reference = manifest["beta0_reference"]
    direct = manifest["beta0_direct_coupled_check"]
    beta0_all_grid = manifest["beta0_all_grid_coupled_subsystem_check"]
    benchmark = manifest["benchmark"]
    audit = manifest["neighbour_audit"]
    quality = manifest["root_quality_summary"]
    runtime = manifest["runtime"]
    ranges = manifest["section_property_ranges"]
    status_lines = "\n".join(
        f"- `{name}`: **{status}**"
        for name, status in manifest["status_gates"].items()
    )
    gap_lines = "\n".join(
        "- $\\beta={beta:.0f}^\\circ$: $\\Delta\\Lambda_{{min}}={gap:.6g}$ "
        "при $A/A_0={alpha:.2f}$ между позициями {left} и {right}.".format(
            beta=float(item["beta_deg"]),
            gap=float(item["minimum_adjacent_Lambda_gap"]),
            alpha=float(item["alpha_A"]),
            left=int(item["between_sorted_positions"][0]),
            right=int(item["between_sorted_positions"][1]),
        )
        for item in manifest["minimum_adjacent_sorted_gaps"]
    )
    beta30_lines = "\n".join(
        "- `k={k}`: $\\Lambda_{{0.70}}={lo:.6f}$, "
        "$\\Lambda_{{1.00}}={mid:.6f}$, $\\Lambda_{{1.30}}={hi:.6f}$; "
        "$\\Delta_{{0.70\\to1.30}}={chg:+.6e}$; "
        "$\\Delta_{{30,k}}={dev:.6e}$; `{cls}`; минимальный соседний "
        "gap {gap:.6g} при $A/A_0={gap_alpha:.2f}$ с позицией {neighbour} "
        "(`{role}`).".format(
            k=int(item["sorted_position"]),
            lo=float(item["Lambda_alpha_0p70"]),
            mid=float(item["Lambda_alpha_1p00"]),
            hi=float(item["Lambda_alpha_1p30"]),
            chg=float(item["relative_change_0p70_to_1p30"]),
            dev=float(item["maximum_relative_deviation_from_alpha_1"]),
            cls=str(item["monotonicity_class"]),
            gap=float(item["minimum_adjacent_Lambda_gap"]),
            gap_alpha=float(item["minimum_gap_alpha_A"]),
            neighbour=int(item["minimum_gap_neighbour_position"]),
            role=str(item["minimum_gap_neighbour_role"]),
        )
        for item in manifest["beta30_shift_summary"]
    )
    timing_note = (
        "полная"
        if runtime["timing_is_complete"]
        else (
            "является нижней границей: времена 14 уже завершённых local-repair "
            "транзакций не были сохранены до прерванной первой финализации"
        )
    )
    return f"""# RLB-2H — axial-stiffness visibility under geometric coupling

## Scope and frozen contract

Проверен конечный ordinary frequency map двух идентичных четырёхслойных
Reddy-стержней. Во всех слоях `angle_deg=0`; геометрия фиксирована:
$\\mu=\\tau=0$, $L_1=L_2=1$, $b=0.20$, $h=0.05$, $K=5/6$.
Рассмотрены только $\\beta=0^\\circ$ и $30^\\circ$ и сетка
$\\alpha_A=A/A_0=0.70:0.02:1.30$.

Материал $M_0$: $E_1=1.1$, $E_2=0.9$, $\\nu_{{12}}=0.3$,
$G_{{12}}=G_{{13}}=G_{{23}}=1/2.6$, $\\rho=1$.
Для наружных и внутренних слоёв использованы
$s_{{outer}}=(4-\\alpha_A)/3$ и
$s_{{inner}}=(7\\alpha_A-4)/3$. Диапазоны множителей:
`s_outer=[{ranges['s_outer']['minimum']:.6g}, {ranges['s_outer']['maximum']:.6g}]`,
`s_inner=[{ranges['s_inner']['minimum']:.6g}, {ranges['s_inner']['maximum']:.6g}]`.

## Constitutive A-only gate

Статус: **{constitutive['status']}**. Получены диапазоны:

- $A/A_0=[{ranges['A_beam_over_A0']['minimum']:.6g},
  {ranges['A_beam_over_A0']['maximum']:.6g}]$;
- $D/D_0=[{ranges['D_beam_over_D0']['minimum']:.16g},
  {ranges['D_beam_over_D0']['maximum']:.16g}]$;
- $S/S_0=[{ranges['S_beam_over_S0']['minimum']:.16g},
  {ranges['S_beam_over_S0']['maximum']:.16g}]$;
- $m/m_0=[{ranges['m_over_m0']['minimum']:.16g},
  {ranges['m_over_m0']['maximum']:.16g}]$;
- $J/J_0=[{ranges['J_over_J0']['minimum']:.16g},
  {ranges['J_over_J0']['maximum']:.16g}]$.

Максимальные residuals: $A$-matrix scaling
{residuals['A_matrix_scale_residual']:.3e}, $D$-matrix invariance
{residuals['D_matrix_invariance_residual']:.3e}, reduced $A/A_0$ formula
{residuals['A_formula_residual']:.3e}, совокупная инвариантность
$D,S,m,J,K,b$ {residuals['only_A_field_max_invariance_residual']:.3e},
$B$ {residuals['B_relative']:.3e}, $I_1$ {residuals['I1_relative']:.3e}.

## Beta=0 reference mechanism

Статус subsystem reference: **{reference['status']}**; direct/coupled check:
**{direct['status']}**. Единственная завершённая coupled-группа при
$\\alpha_A=1$ повторно использована как источник восьми bending roots;
отдельный bending root solve при финализации не выполнялся. Ошибка
инвариантности bending matrix равна
{reference['bending_matrix_invariance_error']:.3e}; ошибка закона
$\\Lambda_{{axial}}\\propto\\alpha_A^{{1/4}}$ —
{reference['axial_scaling_error']:.3e}. Максимальная cluster-aware разность
coupled spectrum и exact axial+bending union в трёх direct-control points:
{direct['maximum_coupled_vs_subsystem_union_relative_difference']:.3e}.
На всей 31-точечной сетке maximum coupled/subsystem difference равна
{beta0_all_grid['maximum_relative_difference']:.3e}; unclassified canonical
slots: {beta0_all_grid['unclassified_canonical_slot_count']}.

В трёх точках $\\alpha_A=0.70,1.00,1.30$ дополнительно проверена
сингулярность direct fixed-fixed boundary matrix общей длины 2 на частотах
exact axial plus single-$\\alpha_A=1$ bending union. Direct builder независим
от joint matrix, однако частоты для этой проверки отдельно не локализовались:
они поданы из subsystem union и использованы только для boundary-singularity
diagnostics. Они не использовались для локализации coupled roots. Максимальная
разность coupled/direct fixed-fixed diagnostic-frequency равна
{direct['maximum_relative_difference']:.3e}.
В coupled/subsystem comparison repeated/near-degenerate roots сравнивались по
multiplicity, total nullity и cluster centre. Direct diagnostic независимо
проверял boundary singularity и detected nullity на поданных частотах, но не
локализовал roots и не устанавливал cluster grouping. Sorted position не
трактуется как modal branch.

## Frequency map and neighbour audit

Полнота: $\\beta=0^\\circ$ — {counts['beta0_points_complete']}/31,
$\\beta=30^\\circ$ — {counts['beta30_points_complete']}/31; BASE rows
{counts['base_rows']}/558; root-9 guards {counts['root9_guards']}/62.
Сохранено {counts['local_refinement_rows']} LOCAL_REFINEMENT rows для
{counts['locally_repaired_points']} точек. Audit отметил
{audit['flagged_point_count']} точек; все они пересчитаны локально;
unresolved points: {audit['unresolved_point_count']}. Сглаживание и подмена
частот predictor/interpolation не применялись.

Минимальные gaps во всём independently sorted spectrum:

{gap_lines}

### Beta=30 finite-grid response

$\\Delta_{{30,k}}=\\max_{{\\alpha_A}}|\\Lambda_k(\\alpha_A)-
\\Lambda_k(1)|/\\Lambda_k(1)$.

{beta30_lines}

Эти значения описывают independently sorted positions на объявленной сетке,
а не прослеженные модальные ветви.

## Numerical quality and performance

Шесть production anchors дали conservative ETA
{benchmark['conservative_eta_seconds']:.1f} s для
{benchmark['remaining_unique_root_points']} оставшихся base points при лимите
{benchmark['eta_limit_seconds']:.0f} s; peak RSS benchmark —
{benchmark['peak_rss_bytes']/2**20:.1f} MiB.

Сумма сохранённых measured times равна
{runtime['total_measured_workflow_seconds_lower_bound']:.1f} s и {timing_note}.
Peak RSS всего workflow: {runtime['peak_rss_bytes']/2**20:.1f} MiB.
Известная нижняя граница determinant/SVD evaluations:
{runtime['determinant_evaluations']}/{runtime['sigma_evaluations']}.
Отсутствуют per-repair evaluation counters для
{runtime['local_repair_metrics_missing_after_interrupted_prefinalization']}
восстановленных repair-транзакций; они не пересчитывались ради метрик.

Максимальный BASE $\\sigma_{{min}}/\\sigma_{{max}}$:
{quality['maximum_base_scaled_sigma_ratio']:.3e}; максимальный boundary
residual: {quality['maximum_base_boundary_null_residual']:.3e}; минимальный
right margin root 9: {quality['minimum_root9_right_margin_Omega']:.6g}.

## Statuses and limitations

{status_lines}

Основной рисунок: `{manifest['plot']['path']}`.
Диаграмма beta=0 mechanism: `{manifest['beta0_diagnostic_plot']['path']}`.

Вывод относится только к четырём равным 0°-слоям, двум идентичным плечам,
двум объявленным углам и конечной сетке $\\alpha_A$. Позиции 1--8 показаны
на рисунке; root 9 служит только completeness guard. Отдельный intentional
поиск или audit spectral tail выше root 9 не выполнялся, и rows выше guard не
экспортировались. Текущая версия скрипта отклоняет distinct pre-trim slot выше
guard; исторические завершённые группы после добавления этого gate не
пересчитывались. Не выполнялись branch tracking, MAC, mode shapes, energy
partition, Ritz, FEM, smoothing, certified audit, изменение production physics,
commit или push.
"""


def finalize_outputs(
    output_dir: Path,
    rows: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    benchmark: Mapping[str, Any],
    point_records: Sequence[Mapping[str, Any]],
    started_at: str,
    total_started: float,
) -> dict[str, Any]:
    target = Path(output_dir)
    spectrum_path = target / SPECTRUM_FILENAME
    resumed_pre_finalization_hash = (
        rlb2e._sha256(spectrum_path) if spectrum_path.is_file() else None
    )
    existing_local_rows = [
        row for row in rows if str(row.get("grid_kind")) == "LOCAL_REFINEMENT"
    ]
    if existing_local_rows:
        audit_rows = rlb2e._read_csv(target / AUDIT_FILENAME)
        repair_records = _recover_existing_repair_records(rows, audit_rows)
        repair_resume_mode = "REUSED_COMPLETED_LOCAL_REPAIR_TRANSACTIONS"
    else:
        audit_rows = neighbour_audit_rows(rows)
        rows, audit_rows, repair_records = apply_local_repairs(rows, audit_rows)
        repair_resume_mode = "COMPUTED_IN_CURRENT_INVOCATION"
    reference_rows, beta0_reference = beta0_reference_rows(rows)
    if beta0_reference["status"] != "PASS":
        raise RuntimeError(f"beta=0 subsystem reference failed: {beta0_reference}")
    # This gate must see determinant/SVD cluster metadata before the optional
    # subsystem classification annotates beta=0 rows.
    beta0_all_grid = beta0_all_grid_subsystem_check(rows, reference_rows)
    if beta0_all_grid["status"] != "PASS":
        raise RuntimeError(
            f"beta=0 all-grid subsystem comparison failed: {beta0_all_grid}"
        )
    apply_beta0_reference_classification(rows, reference_rows)
    unclassified_beta0_slots = sum(
        1
        for row in rows
        if abs(float(row["beta_deg"])) <= 1.0e-12
        and _as_bool(row.get("is_canonical_plot_source", True))
        and str(row.get("decoupled_reference_classification", ""))
        == "UNCLASSIFIED"
    )
    beta0_all_grid["unclassified_canonical_slot_count"] = (
        unclassified_beta0_slots
    )
    if unclassified_beta0_slots:
        beta0_all_grid["status"] = "FAIL"
        raise RuntimeError(
            "beta=0 subsystem classification left canonical slots unclassified."
        )
    spectrum_metadata_backfill = backfill_resumed_spectrum_metadata(rows)
    point_record_backfills = backfill_resumed_point_record_metadata(
        point_records, rows
    )
    final_schema_audit = audit_final_spectrum_schema(rows)
    if final_schema_audit["status"] != "PASS":
        raise RuntimeError(
            f"Final spectrum provenance/schema audit failed: {final_schema_audit}"
        )
    rlb2e._atomic_write_csv(target / SPECTRUM_FILENAME, rows, SPECTRUM_FIELDS)
    rlb2e._atomic_write_csv(target / AUDIT_FILENAME, audit_rows)
    rlb2e._atomic_write_csv(target / REFERENCE_FILENAME, reference_rows)
    spectrum_audit = audit_spectrum_rows(rows)
    if spectrum_audit["status"] != "PASS":
        raise RuntimeError(f"Final spectrum audit failed: {spectrum_audit}")
    root_quality = _root_quality_summary(rows)
    if root_quality["status"] == "FAIL":
        raise RuntimeError(f"Canonical root quality failed: {root_quality}")
    beta0_direct = coupled_beta0_direct_checks(rows, reference_rows)
    if beta0_direct["status"] != "PASS":
        raise RuntimeError(f"beta=0 coupled/direct reference failed: {beta0_direct}")
    plot = create_plot_from_csv(target)
    beta0_plot = create_beta0_reference_plot(target)
    unresolved = sum(1 for record in repair_records if record["status"] == "UNRESOLVED")
    physics_hashes = {
        path: rlb2e._sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS
    }
    physics_preserved = physics_hashes == EXPECTED_PRODUCTION_PHYSICS_HASHES
    predecessor_hashes = {
        name: _sha_tree(path) for name, path in PREDECESSOR_RESULT_DIRS.items()
    }
    predecessor_preserved = predecessor_hashes == INITIAL_PREDECESSOR_TREE_HASHES
    if not physics_preserved:
        raise RuntimeError("Frozen production-physics hashes changed during RLB-2H.")
    if not predecessor_preserved:
        raise RuntimeError("A frozen predecessor result tree changed during RLB-2H.")
    scientific_status = (
        "PASS"
        if unresolved == 0 and root_quality["status"] == "PASS"
        else "PARTIAL_PASS"
    )
    frequency_map_status = (
        "PASS"
        if unresolved == 0 and spectrum_audit["status"] == "PASS"
        else "PARTIAL_PASS"
    )
    status_gates = {
        "RLB-2H-CONSTITUTIVE-A-ONLY": str(constitutive["status"]),
        "RLB-2H-BETA0-DIRECT-SUBSYSTEM-REFERENCE": (
            "PASS"
            if beta0_direct["status"] == "PASS"
            and beta0_all_grid["status"] == "PASS"
            else "FAIL"
        ),
        "RLB-2H-FREQUENCY-MAP": frequency_map_status,
        "RLB-2H-ROOT-INVENTORY": str(root_quality["status"]),
        "OVERALL": scientific_status,
    }
    counts = {
        "base_rows": spectrum_audit["base_row_count"],
        "base_groups": spectrum_audit["base_group_count"],
        "beta0_points_complete": sum(
            _group_key(0.0, float(alpha)) in _complete_base_group_index(rows)
            for alpha in alpha_grid()
        ),
        "beta30_points_complete": sum(
            _group_key(30.0, float(alpha)) in _complete_base_group_index(rows)
            for alpha in alpha_grid()
        ),
        "local_refinement_rows": sum(str(row["grid_kind"]) == "LOCAL_REFINEMENT" for row in rows),
        "root9_guards": sum(
            str(row["grid_kind"]) == "BASE" and int(row["sorted_position"]) == K_GUARD
            for row in rows
        ),
        "locally_repaired_points": len(repair_records),
        "unresolved_points": unresolved,
    }
    runtime = _runtime_summary(
        point_records,
        repair_records,
        beta0_direct["direct_reference_records"],
        beta0_reference,
        plot,
        beta0_plot,
        finalization_process_seconds=time.perf_counter() - total_started,
    )
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "scientific_status": scientific_status,
        "completed_at_utc": rlb2e._utc_now(),
        "git": rlb2e._git_state(),
        "contract": contract_payload(),
        "contract_sha256": contract_hash(),
        "production_physics_hashes": physics_hashes,
        "expected_production_physics_hashes": EXPECTED_PRODUCTION_PHYSICS_HASHES,
        "production_physics_preserved": physics_preserved,
        "analysis_script_sha256": rlb2e._sha256(Path(__file__)),
        "constitutive_gate": dict(constitutive),
        "section_property_ranges": _section_property_range_summary(target),
        "counts": counts,
        "spectrum_audit": spectrum_audit,
        "final_spectrum_schema_audit": final_schema_audit,
        "root_quality_summary": root_quality,
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
        "beta0_reference": beta0_reference,
        "beta0_direct_coupled_check": beta0_direct,
        "beta0_all_grid_coupled_subsystem_check": beta0_all_grid,
        "minimum_adjacent_sorted_gaps": _minimum_adjacent_gaps(rows),
        "beta30_shift_summary": _beta30_shift_summary(rows),
        "plot": plot,
        "beta0_diagnostic_plot": beta0_plot,
        "status_gates": status_gates,
        "runtime": runtime,
        "resume_and_repair_provenance": {
            "mode": repair_resume_mode,
            "resumed_pre_finalization_spectrum_sha256": resumed_pre_finalization_hash,
            "existing_local_refinement_rows_reused": len(existing_local_rows),
            "existing_local_repair_points_reused": len(repair_records),
            "completed_base_points_recomputed": 0,
            "spectrum_metadata_backfill": spectrum_metadata_backfill,
            "point_records_with_metadata_backfill": point_record_backfills,
            "pretrim_distinct_root_above_guard_gate_in_current_script": True,
            "historical_base_groups_recomputed_after_gate_addition": 0,
        },
        "root_above_guard_provenance": {
            "exported_rows_above_guard": 0,
            "intentional_tail_search_or_audit_run": False,
            "current_script_rejects_distinct_pretrim_slot_above_guard": True,
            "same_guard_cluster_continuation_is_retained_as_multiplicity_metadata": True,
            "historical_pretrim_distinct_candidate_count_persisted": False,
            "historical_completed_groups_recomputed_for_this_gate": 0,
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
        "initial_predecessor_result_tree_hashes": INITIAL_PREDECESSOR_TREE_HASHES,
        "predecessor_result_tree_hashes": predecessor_hashes,
        "predecessor_result_trees_preserved": predecessor_preserved,
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
    checkpoint["root_calculation_count"] = runtime["total_root_solve_count"]
    checkpoint["local_repair_count"] = len(repair_records)
    checkpoint["terminal_unresolved_points"] = [
        {
            "beta_deg": record["beta_deg"],
            "alpha_A": record["alpha_A"],
            "repair_id": record["repair_id"],
        }
        for record in repair_records
        if record["status"] == "UNRESOLVED"
    ]
    rlb2e._atomic_write_json(target / CHECKPOINT_FILENAME, checkpoint)
    report = _report_text_v2(manifest)
    rlb2e._atomic_write_text(target / REPORT_FILENAME, report)
    manifest["output_hashes"] = _output_hashes(target)
    rlb2e._atomic_write_json(target / MANIFEST_FILENAME, manifest)
    return manifest


def _results_are_complete(rows: Sequence[Mapping[str, Any]], target: Path) -> bool:
    required = {
        SECTION_FILENAME,
        SPECTRUM_FILENAME,
        AUDIT_FILENAME,
        REFERENCE_FILENAME,
        BENCHMARK_FILENAME,
        CHECKPOINT_FILENAME,
        PLOT_FILENAME,
        BETA0_PLOT_FILENAME,
        REPORT_FILENAME,
        MANIFEST_FILENAME,
    }
    if (
        audit_spectrum_rows(rows)["status"] != "PASS"
        or audit_final_spectrum_schema(rows)["status"] != "PASS"
        or not all((target / name).is_file() for name in required)
    ):
        return False
    try:
        manifest = json.loads(
            (target / MANIFEST_FILENAME).read_text(encoding="utf-8")
        )
    except (OSError, UnicodeError, json.JSONDecodeError):
        return False
    if not (
        manifest.get("stage_id") == STAGE_ID
        and manifest.get("contract_sha256") == contract_hash()
        and manifest.get("analysis_script_sha256") == rlb2e._sha256(Path(__file__))
        and manifest.get("scientific_status") == "PASS"
        and set(manifest.get("status_gates", {})) == REQUIRED_STATUS_GATE_NAMES
        and set(manifest.get("status_gates", {}).values()) == {"PASS"}
        and manifest.get("production_physics_preserved") is True
        and manifest.get("predecessor_result_trees_preserved") is True
    ):
        return False
    current_physics_hashes = {
        path: rlb2e._sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS
    }
    current_predecessor_hashes = {
        name: _sha_tree(path) for name, path in PREDECESSOR_RESULT_DIRS.items()
    }
    if (
        current_physics_hashes != EXPECTED_PRODUCTION_PHYSICS_HASHES
        or current_predecessor_hashes != INITIAL_PREDECESSOR_TREE_HASHES
        or manifest.get("production_physics_hashes") != current_physics_hashes
        or manifest.get("predecessor_result_tree_hashes")
        != current_predecessor_hashes
    ):
        return False
    hashes = manifest.get("output_hashes")
    required_hashed = required - {MANIFEST_FILENAME}
    if not isinstance(hashes, dict) or not required_hashed <= set(hashes):
        return False
    return all(
        (target / name).is_file()
        and rlb2e._sha256(target / name) == str(expected)
        for name, expected in hashes.items()
    )


def run_workflow(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    *,
    missing_only: bool = False,
) -> dict[str, Any]:
    total_started = time.perf_counter()
    started_at = rlb2e._utc_now()
    target = Path(output_dir)
    rows = _existing_rows(target)
    if rows and _results_are_complete(rows, target):
        manifest = json.loads((target / MANIFEST_FILENAME).read_text(encoding="utf-8"))
        manifest["resume_root_calculation_count"] = 0
        manifest["resume_outputs_modified"] = False
        return manifest
    checkpoint_path = target / CHECKPOINT_FILENAME
    if checkpoint_path.is_file():
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        point_records = _validated_checkpoint_point_records(checkpoint, rows)
    else:
        point_records = []
    constitutive, _section_rows = prepare_constitutive_outputs(target)
    rows, benchmark = run_benchmarks(target, rows, point_records, constitutive, started_at)
    if not benchmark["production_run_permitted"]:
        return {
            "stage_id": STAGE_ID,
            "scientific_status": "NOT_EVALUATED",
            "workflow_status": "STOPPED_BY_ETA_GATE",
            "benchmark": benchmark,
        }
    rows = complete_missing_points(target, rows, point_records, constitutive, started_at)
    if missing_only and audit_spectrum_rows(rows)["status"] == "PASS":
        return finalize_outputs(
            target,
            rows,
            constitutive,
            benchmark,
            point_records,
            started_at,
            total_started,
        )
    return finalize_outputs(
        target,
        rows,
        constitutive,
        benchmark,
        point_records,
        started_at,
        total_started,
    )


def manifest_only(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    gate = constitutive_gate()
    rows = section_property_rows()
    return {
        "stage_id": STAGE_ID,
        "mode": "manifest_only",
        "scientific_status": "NOT_EVALUATED",
        "root_calculation_count": 0,
        "contract": contract_payload(),
        "constitutive_gate": gate,
        "section_property_row_count": len(rows),
        "output_directory": Path(output_dir).as_posix(),
        "resume_artifacts_modified": False,
        "created_at_utc": rlb2e._utc_now(),
    }


def plot_only_workflow(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    """Redraw both figures from immutable, manifest-verified CSV inputs."""

    target = Path(output_dir)
    manifest_path = target / MANIFEST_FILENAME
    if not manifest_path.is_file():
        raise RuntimeError("plot_only requires an existing completed manifest.")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        manifest.get("stage_id") != STAGE_ID
        or manifest.get("contract_sha256") != contract_hash()
        or manifest.get("analysis_script_sha256")
        != rlb2e._sha256(Path(__file__))
        or manifest.get("scientific_status") != "PASS"
        or set(manifest.get("status_gates", {})) != REQUIRED_STATUS_GATE_NAMES
        or set(manifest.get("status_gates", {}).values()) != {"PASS"}
    ):
        raise RuntimeError("plot_only rejected a foreign or stale RLB-2H contract.")
    current_physics_hashes = {
        path: rlb2e._sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS
    }
    current_predecessor_hashes = {
        name: _sha_tree(path) for name, path in PREDECESSOR_RESULT_DIRS.items()
    }
    if (
        current_physics_hashes != EXPECTED_PRODUCTION_PHYSICS_HASHES
        or current_predecessor_hashes != INITIAL_PREDECESSOR_TREE_HASHES
        or manifest.get("production_physics_hashes") != current_physics_hashes
        or manifest.get("predecessor_result_tree_hashes")
        != current_predecessor_hashes
    ):
        raise RuntimeError("plot_only rejected changed frozen inputs.")
    hashes = manifest.get("output_hashes", {})
    if not isinstance(hashes, dict):
        raise RuntimeError("plot_only rejected malformed output hashes.")
    actual_files = {path.name for path in target.iterdir() if path.is_file()}
    expected_files = set(hashes) | {MANIFEST_FILENAME}
    if actual_files != expected_files:
        raise RuntimeError("plot_only rejected an unexpected result-tree file set.")
    mutable_figures = {PLOT_FILENAME, BETA0_PLOT_FILENAME}
    for name, expected in hashes.items():
        if name in mutable_figures:
            continue
        path = target / name
        if not path.is_file() or rlb2e._sha256(path) != str(expected):
            raise RuntimeError(f"plot_only rejected modified evidence: {name}")
    # Render both figures completely before replacing either production PNG.
    with tempfile.TemporaryDirectory(prefix=".rlb2h_plot_only_", dir=target) as raw:
        temporary_dir = Path(raw)
        shutil.copyfile(
            target / SPECTRUM_FILENAME,
            temporary_dir / SPECTRUM_FILENAME,
        )
        shutil.copyfile(
            target / REFERENCE_FILENAME,
            temporary_dir / REFERENCE_FILENAME,
        )
        plot = create_plot_from_csv(temporary_dir)
        beta0_plot = create_beta0_reference_plot(temporary_dir)
        os.replace(temporary_dir / PLOT_FILENAME, target / PLOT_FILENAME)
        os.replace(
            temporary_dir / BETA0_PLOT_FILENAME,
            target / BETA0_PLOT_FILENAME,
        )
    plot["path"] = (target / PLOT_FILENAME).as_posix()
    beta0_plot["path"] = (target / BETA0_PLOT_FILENAME).as_posix()
    manifest["plot"] = plot
    manifest["beta0_diagnostic_plot"] = beta0_plot
    manifest["plot_only_last_rendered_at_utc"] = rlb2e._utc_now()
    manifest["plot_only_root_calculation_count"] = 0
    manifest["output_hashes"] = _output_hashes(target)
    rlb2e._atomic_write_json(manifest_path, manifest)
    return {
        "mode": "plot_only",
        "plot": plot,
        "beta0_diagnostic_plot": beta0_plot,
        "root_calculation_count": 0,
        "matrix_assembly_count": 0,
        "determinant_evaluation_count": 0,
        "SVD_evaluation_count": 0,
        "manifest_updated": True,
    }


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--missing-only", action="store_true")
    mode.add_argument("--plot-only", action="store_true")
    mode.add_argument("--manifest-only", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.manifest_only:
        result = manifest_only(args.output_dir)
    elif args.plot_only:
        result = plot_only_workflow(args.output_dir)
    else:
        result = run_workflow(args.output_dir, missing_only=args.missing_only)
    print(json.dumps(rlb2e._json_value(result), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
