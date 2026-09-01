"""RLB-2J six-ply pairwise stiffness-transfer frequency maps.

This analysis-level entry point implements a local ``frequency-map-v1``
instance for three symmetric six-ply, all-zero-degree laminates.  The three
families preserve membrane, transverse-shear and mass properties while their
bending stiffness changes with a signed transfer parameter ``xi``.  Spectral
positions are independently sorted at every point: positions 1--8 are plotted
and root 9 is retained only as a completeness guard.

The characteristic matrix and every accepted root come from the frozen RLB
production modules and the established RLB-2E determinant/SVD machinery.
Continuation values are numerical locators only.  ``--plot-only`` imports no
production physics and performs no laminate, matrix or root calculation.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Iterable, Mapping, Sequence

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


FloatArray = NDArray[np.float64]

STAGE_ID = "RLB-2J"
ALGORITHM_VERSION = "six_ply_pairwise_stiffness_transfer_fast_plot_v1"
POLICY_ID = "frequency-map-v1"

CENTER_MIDDLE_TRANSFER = "CENTER_MIDDLE_TRANSFER"
MIDDLE_OUTER_TRANSFER = "MIDDLE_OUTER_TRANSFER"
CENTER_OUTER_TRANSFER = "CENTER_OUTER_TRANSFER"
CONFIGURATIONS = (
    CENTER_MIDDLE_TRANSFER,
    MIDDLE_OUTER_TRANSFER,
    CENTER_OUTER_TRANSFER,
)

PAIR_CENTER = "CENTER"
PAIR_MIDDLE = "MIDDLE"
PAIR_OUTER = "OUTER"
PAIR_ORDER = (PAIR_CENTER, PAIR_MIDDLE, PAIR_OUTER)
STACK_BOTTOM_TO_TOP = (
    PAIR_OUTER,
    PAIR_MIDDLE,
    PAIR_CENTER,
    PAIR_CENTER,
    PAIR_MIDDLE,
    PAIR_OUTER,
)
PAIR_A_WEIGHTS = {
    PAIR_CENTER: Fraction(1, 3),
    PAIR_MIDDLE: Fraction(1, 3),
    PAIR_OUTER: Fraction(1, 3),
}
PAIR_D_WEIGHTS = {
    PAIR_CENTER: Fraction(1, 27),
    PAIR_MIDDLE: Fraction(7, 27),
    PAIR_OUTER: Fraction(19, 27),
}
TRANSFER_VECTORS = {
    CENTER_MIDDLE_TRANSFER: (-1, 1, 0),
    MIDDLE_OUTER_TRANSFER: (0, -1, 1),
    CENTER_OUTER_TRANSFER: (-1, 0, 1),
}
TRANSFER_LEVERS = {
    CENTER_MIDDLE_TRANSFER: 1,
    MIDDLE_OUTER_TRANSFER: 2,
    CENTER_OUTER_TRANSFER: 3,
}
TRANSFER_METADATA = {
    CENTER_MIDDLE_TRANSFER: {
        "fixed_pair": PAIR_OUTER,
        "donor_pair_for_positive_xi": PAIR_CENTER,
        "receiver_pair_for_positive_xi": PAIR_MIDDLE,
    },
    MIDDLE_OUTER_TRANSFER: {
        "fixed_pair": PAIR_CENTER,
        "donor_pair_for_positive_xi": PAIR_MIDDLE,
        "receiver_pair_for_positive_xi": PAIR_OUTER,
    },
    CENTER_OUTER_TRANSFER: {
        "fixed_pair": PAIR_MIDDLE,
        "donor_pair_for_positive_xi": PAIR_CENTER,
        "receiver_pair_for_positive_xi": PAIR_OUTER,
    },
}

MU = 0.0
TAU = 0.0
BETA_DEG = 30.0
L_REFERENCE = 1.0
L1 = 1.0
L2 = 1.0
L_TOTAL = 2.0
WIDTH = 0.20
THICKNESS = 0.05
PLY_THICKNESS = THICKNESS / 6.0
K = 5.0 / 6.0

E1_0 = 1.1
E2_0 = 0.9
NU12_0 = 0.3
G12_0 = 1.0 / 2.6
G13_0 = 1.0 / 2.6
G23_0 = 1.0 / 2.6
RHO_0 = 1.0

XI_MIN = -0.80
XI_MAX = 0.80
XI_STEP = 0.02
K_PLOT = 8
K_GUARD = 9

REFERENCE_AREA = WIDTH * THICKNESS
IY_REFERENCE = WIDTH * THICKNESS**3 / 12.0
OMEGA_TO_OMEGA_SCALE = L_REFERENCE**2 * math.sqrt(
    RHO_0 * REFERENCE_AREA / IY_REFERENCE
)

EXACT_IDENTITY_TOLERANCE = 1.0e-14
MATRIX_RELATIVE_TOLERANCE = 1.0e-12
SYMMETRY_RELATIVE_TOLERANCE = 1.0e-12
REDUCED_PROPERTY_TOLERANCE = 1.0e-11
REDUCTION_ROUTE_TOLERANCE = 1.0e-11
SPECTRAL_RELATIVE_TOLERANCE = 1.0e-8
ROOT_SINGULAR_RATIO_TOLERANCE = 1.0e-9
BOUNDARY_RESIDUAL_TOLERANCE = 1.0e-9
ROOT9_RIGHT_TAIL_OMEGA = 2.0
ETA_LIMIT_SECONDS = 45.0 * 60.0
NEIGHBOUR_MAD_MULTIPLIER = 8.0
NEIGHBOUR_ABSOLUTE_TRIGGER = 1.0e-3
SLOPE_INSENSITIVE_ABSOLUTE = 1.0e-10
SLOPE_DIAGNOSTIC_RELATIVE = 5.0e-4

DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_six_ply_pairwise_stiffness_transfer"
)
SPECTRUM_FILENAME = "spectrum_roots.csv"
SECTION_FILENAME = "section_properties.csv"
AUDIT_FILENAME = "neighbour_audit.csv"
MATCHED_D_FILENAME = "matched_D_collapse.csv"
SLOPE_FILENAME = "initial_slope_check.csv"
BENCHMARK_FILENAME = "benchmark.json"
CHECKPOINT_FILENAME = "checkpoint.json"
XI_PLOT_FILENAME = "lambda_vs_xi_three_pair_transfers.png"
MASTER_PLOT_FILENAME = "lambda_vs_D_ratio_master_collapse.png"
REPORT_FILENAME = "report.md"
MANIFEST_FILENAME = "run_manifest.json"
FOUR_VS_SIX_FILENAME = "four_vs_six_ply_equivalence.csv"

PRODUCTION_PHYSICS_PATHS = (
    "scripts/lib/reddy_symmetric_laminated_beam.py",
    "scripts/lib/reddy_symmetric_coupled_beams.py",
    "scripts/lib/reddy_inplane_geometry.py",
)
EXPECTED_PRODUCTION_PHYSICS_HASHES = {
    "scripts/lib/reddy_symmetric_laminated_beam.py": "9E3F94747FA3723D0FEE350562F29A0DB070C3E3A17DDCCA3795F1E69AEDBE4B",
    "scripts/lib/reddy_symmetric_coupled_beams.py": "E70F7AF5B4BB61AA90525664E6C4834EF5A003F34B23D6C2741583D38DAAD9A7",
    "scripts/lib/reddy_inplane_geometry.py": "C46A42C462264BC27C99C358AABD7DF49F94F928A60D8150FD320D8DFB37E99E",
}
PREDECESSOR_RESULT_DIRS = {
    "RLB-2A": ROOT / "results/laminated_beams/reddy_cross_ply_layer_order_pilot",
    "RLB-2D": ROOT / "results/laminated_beams/rectangular_weakly_orthotropic_graphs_mu_beta",
    "RLB-2E": ROOT / "results/laminated_beams/reddy_stiffness_layout_contrast_sweep",
    "RLB-2F": ROOT / "results/laminated_beams/reddy_one_arm_layered_contrast_sweep",
    "RLB-2G": ROOT / "results/laminated_beams/reddy_mass_layout_duality",
    "RLB-2H": ROOT / "results/laminated_beams/reddy_axial_stiffness_visibility",
    "RLB-2I": ROOT / "results/laminated_beams/reddy_six_ply_equivalent_laminates",
}
EXPECTED_PREDECESSOR_TREE_HASHES = {
    "RLB-2A": "07D0B115FE8B42AC4EF11A32B3B37A075A5477203869BDB4A343F46C79533106",
    "RLB-2D": "86D34750EB13CE6039D8FFA18D9FE15A4CC518FCD7921A5646F3ADB0129F0250",
    "RLB-2E": "57E9FFCFD518FADF198F30C84F04EE181F1C645814A7BCFA834FCC920426B008",
    "RLB-2F": "10C80E8136AF917BCC5EFCB351FAD2FBE6665856A91C1F71241BC650372046C5",
    "RLB-2G": "4DA662EB77240C59B78017CCCD38136522561F2F3BE48BAD2FA50AACBB059CC1",
    "RLB-2H": "8669AFC4D8BF4831D1DBC295C4A6EABA4FA16A864B77810DAC5F652030ECA874",
    "RLB-2I": "A14CAD4FA1A548BEAA2761EB037888E05581DF12D23B5A1785214B223D3A5085",
}

FREQUENCY_MAP_POLICY = {
    "frequency_map_policy": POLICY_ID,
    "calculation_mode": "fast_plot",
    "spectrum_semantics": "sorted_positions",
    "sweep_parameter": "xi",
    "parameter_grid": "-0.80:0.02:0.80",
    "continuation_anchor": 0.0,
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
    "row_id", "configuration_id", "xi", "xi_index", "s_C", "s_M", "s_O",
    "D_over_D0", "D_key_q", "grid_kind", "physical_solve_id", "transaction_id",
    "continuation_leg", "sorted_position", "root_role", "guard_flag", "omega",
    "Omega", "Lambda", "predictor_Omega", "predictor_used_as_final",
    "locator_interval_left_Omega", "locator_interval_right_Omega",
    "root_interval_left_Omega", "root_interval_right_Omega",
    "detector_refiner_provenance", "raw_determinant", "scaled_determinant",
    "raw_sigma_ratio", "scaled_sigma_ratio", "boundary_null_residual",
    "detected_nullity", "cluster_id", "cluster_multiplicity", "cluster_total_nullity",
    "unresolved_candidates_below_root9", "search_right_Omega",
    "root9_right_margin_Omega", "solve_mode", "fallback_used", "quality_status",
    "point_status", "shared_xi0_anchor_reused", "shared_xi0_source_configuration",
    "is_canonical_plot_source", "supersedes_row_id", "repair_id",
    "roots_above_9_computed",
)

SECTION_FIELDS = (
    "configuration_id", "xi", "xi_index", "transfer_vector_C_M_O", "fixed_pair",
    "donor_pair_for_positive_xi", "receiver_pair_for_positive_xi", "s_C", "s_M",
    "s_O", "sum_multipliers", "D_key_q", "A_beam", "A_over_A0", "D_beam",
    "D_over_D0", "expected_D_over_D0", "S_beam", "S_over_S0", "m", "m_over_m0",
    "J", "J_over_J0", "K", "width", "stack_bottom_to_top", "ply_angles_deg",
    "ply_thicknesses", "z_interfaces", "A_matrix", "B_matrix", "D_matrix",
    "shear_matrix_yz_xz", "I0", "I1", "I2", "B_relative", "I1_relative",
    "A_invariance_relative", "A_matrix_invariance_relative", "D_formula_relative",
    "D_matrix_formula_relative", "shear_matrix_invariance_relative",
    "S_invariance_relative", "m_invariance_relative", "J_invariance_relative",
    "reduction_route_max_relative", "constitutive_status",
)


@dataclass(frozen=True)
class SectionObjects:
    laminate: Any
    properties: Any
    multipliers: Mapping[str, float]


@dataclass(frozen=True)
class PointSolution:
    configuration_id: str
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
    candidate_count_total: int = K_GUARD
    accepted_candidates_above_root9: int = 0
    retained_slots_above_root9: int = 0
    roots_above_9_computed: bool = False


ROOT_CALCULATION_COUNT = 0
LAMINATE_CALCULATION_COUNT = 0
MATRIX_ASSEMBLY_COUNT = 0


def _base() -> Any:
    from scripts.analysis.laminated_beams import (
        sweep_reddy_stiffness_layout_contrast as base,
    )
    return base


def _physics_modules() -> tuple[Any, Any]:
    from scripts.lib import reddy_symmetric_coupled_beams as coupled
    from scripts.lib import reddy_symmetric_laminated_beam as beam
    return beam, coupled


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _sha_tree(path: Path) -> str | None:
    root = Path(path)
    if not root.is_dir():
        return None
    digest = hashlib.sha256()
    for item in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
        relative = item.relative_to(root).as_posix().encode("utf-8")
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(bytes.fromhex(_sha256(item)))
    return digest.hexdigest().upper()


def _json_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, Fraction):
        return f"{value.numerator}/{value.denominator}"
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _csv_value(value: Any) -> Any:
    converted = _json_value(value)
    if isinstance(converted, (dict, list)):
        return json.dumps(converted, ensure_ascii=False, separators=(",", ":"))
    if isinstance(converted, bool):
        return "true" if converted else "false"
    return converted


def _atomic_write_text(path: Path, text: str) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".tmp")
    try:
        temporary.write_text(text, encoding="utf-8")
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    _atomic_write_text(
        path,
        json.dumps(_json_value(payload), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )


def _atomic_write_csv(
    path: Path, rows: Iterable[Mapping[str, Any]], fields: Sequence[str] | None = None
) -> None:
    records = [dict(row) for row in rows]
    fieldnames = list(fields or (records[0].keys() if records else ()))
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="raise")
            writer.writeheader()
            for row in records:
                writer.writerow({name: _csv_value(row.get(name, "")) for name in fieldnames})
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _git_state() -> dict[str, str]:
    def run(*args: str) -> str:
        return subprocess.run(
            ["git", *args], cwd=ROOT, check=True, capture_output=True, text=True
        ).stdout.strip()
    return {
        "top_level": run("rev-parse", "--show-toplevel").replace("\\", "/"),
        "branch": run("branch", "--show-current"),
        "head": run("rev-parse", "HEAD"),
        "last_commit": run("log", "-1", "--oneline"),
        "status_short": run("status", "--short"),
    }


def _peak_rss_bytes() -> int:
    if os.name == "nt":
        class Counters(ctypes.Structure):
            _fields_ = [
                ("cb", ctypes.c_ulong), ("PageFaultCount", ctypes.c_ulong),
                ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t), ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t),
            ]
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        psapi = ctypes.WinDLL("psapi", use_last_error=True)
        kernel32.GetCurrentProcess.restype = ctypes.c_void_p
        psapi.GetProcessMemoryInfo.argtypes = (
            ctypes.c_void_p,
            ctypes.POINTER(Counters),
            ctypes.c_ulong,
        )
        psapi.GetProcessMemoryInfo.restype = ctypes.c_int
        counters = Counters()
        counters.cb = ctypes.sizeof(counters)
        ok = psapi.GetProcessMemoryInfo(
            kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb
        )
        return int(counters.PeakWorkingSetSize) if ok else 0
    try:
        import resource
        value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        return value if sys.platform == "darwin" else value * 1024
    except (ImportError, OSError):
        return 0


def _relative(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / max(
        abs(float(left)), abs(float(right)), np.finfo(float).tiny
    )


def _matrix_relative(left: Any, right: Any) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    scale = max(float(np.linalg.norm(a, ord="fro")), float(np.linalg.norm(b, ord="fro")), np.finfo(float).tiny)
    return float(np.linalg.norm(a - b, ord="fro") / scale)


def _scaled_B(laminate: Any) -> float:
    scale = max(float(np.linalg.norm(laminate.A, ord="fro")) * laminate.thickness, np.finfo(float).tiny)
    return float(np.linalg.norm(laminate.B, ord="fro") / scale)


def _scaled_I1(laminate: Any) -> float:
    scale = max(abs(float(laminate.I0)) * laminate.thickness, np.finfo(float).tiny)
    return abs(float(laminate.I1)) / scale


def _reduction_route_max(properties: Any) -> float:
    return max(
        float(item.relative_difference)
        for item in (
            properties.axial_reduction,
            properties.bending_reduction,
            properties.shear_reduction_before_K,
        )
        if item is not None
    )


def xi_grid() -> FloatArray:
    return np.asarray([index / 50.0 for index in range(-40, 41)], dtype=float)


def xi_index(xi: float) -> int:
    index = int(round(float(xi) * 50.0))
    if index < -40 or index > 40 or abs(float(xi) - index / 50.0) > 1.0e-12:
        raise ValueError("xi must lie on the exact -0.80:0.02:0.80 grid.")
    return index


def stiffness_multipliers(configuration_id: str, xi: float) -> dict[str, float]:
    value = xi_index(xi) / 50.0
    if configuration_id == CENTER_MIDDLE_TRANSFER:
        result = {PAIR_CENTER: 1.0 - value, PAIR_MIDDLE: 1.0 + value, PAIR_OUTER: 1.0}
    elif configuration_id == MIDDLE_OUTER_TRANSFER:
        result = {PAIR_CENTER: 1.0, PAIR_MIDDLE: 1.0 - value, PAIR_OUTER: 1.0 + value}
    elif configuration_id == CENTER_OUTER_TRANSFER:
        result = {PAIR_CENTER: 1.0 - value, PAIR_MIDDLE: 1.0, PAIR_OUTER: 1.0 + value}
    else:
        raise ValueError(f"Unknown transfer configuration: {configuration_id!r}.")
    if min(result.values()) <= 0.0:
        raise ValueError("All pair stiffness multipliers must remain positive.")
    return result


def D_key_q(configuration_id: str, xi: float) -> int:
    return TRANSFER_LEVERS[configuration_id] * xi_index(xi)


def exact_D_ratio(configuration_id: str, xi: float) -> Fraction:
    return Fraction(225 + D_key_q(configuration_id, xi), 225)


def expected_D_ratio(configuration_id: str, xi: float) -> float:
    return float(exact_D_ratio(configuration_id, xi))


def base_material_contract() -> dict[str, float]:
    return {
        "E1_0": E1_0, "E2_0": E2_0, "nu12_0": NU12_0,
        "G12_0": G12_0, "G13_0": G13_0, "G23_0": G23_0, "rho_0": RHO_0,
    }


def build_scaled_material(scale: float, pair_id: str) -> Any:
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
        name=f"RLB-2J {pair_id}(scale={factor:.6f})",
    )


def build_six_ply_section(configuration_id: str, xi: float) -> SectionObjects:
    global LAMINATE_CALCULATION_COUNT
    LAMINATE_CALCULATION_COUNT += 1
    beam, _coupled = _physics_modules()
    multipliers = stiffness_multipliers(configuration_id, xi)
    materials = {
        pair_id: build_scaled_material(multiplier, pair_id)
        for pair_id, multiplier in multipliers.items()
    }
    laminate = beam.integrate_laminate(
        [
            beam.Ply(materials[pair_id], 0.0, PLY_THICKNESS, label=pair_id)
            for pair_id in STACK_BOTTOM_TO_TOP
        ]
    )
    properties = beam.reduce_to_beam_properties(
        laminate,
        width=WIDTH,
        K=K,
        symmetry_tolerance=SYMMETRY_RELATIVE_TOLERANCE,
        reduction_tolerance=REDUCTION_ROUTE_TOLERANCE,
    )
    return SectionObjects(laminate, properties, multipliers)


def _baseline_section() -> SectionObjects:
    return build_six_ply_section(CENTER_OUTER_TRANSFER, 0.0)


def compute_constitutive_data() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    baseline = _baseline_section()
    rows: list[dict[str, Any]] = []
    maxima = {
        "sum_multiplier_absolute": 0.0,
        "A_matrix_invariance_relative": 0.0,
        "A_beam_invariance_relative": 0.0,
        "D_matrix_formula_relative": 0.0,
        "D_beam_formula_relative": 0.0,
        "shear_matrix_invariance_relative": 0.0,
        "S_beam_invariance_relative": 0.0,
        "mass_invariance_relative": 0.0,
        "rotary_inertia_invariance_relative": 0.0,
        "B_scaled_residual": 0.0,
        "I1_scaled_residual": 0.0,
        "reduction_route_relative": 0.0,
    }
    passed = True
    for configuration_id in CONFIGURATIONS:
        for xi_value in xi_grid():
            xi = float(xi_value)
            section = build_six_ply_section(configuration_id, xi)
            p = section.properties
            lam = section.laminate
            expected = expected_D_ratio(configuration_id, xi)
            values = {
                "sum_multiplier_absolute": abs(sum(section.multipliers.values()) - 3.0),
                "A_matrix_invariance_relative": _matrix_relative(lam.A, baseline.laminate.A),
                "A_beam_invariance_relative": _relative(p.A, baseline.properties.A),
                "D_matrix_formula_relative": _matrix_relative(lam.D, expected * baseline.laminate.D),
                "D_beam_formula_relative": _relative(p.D / baseline.properties.D, expected),
                "shear_matrix_invariance_relative": _matrix_relative(lam.shear, baseline.laminate.shear),
                "S_beam_invariance_relative": _relative(p.S, baseline.properties.S),
                "mass_invariance_relative": _relative(p.m, baseline.properties.m),
                "rotary_inertia_invariance_relative": _relative(p.J, baseline.properties.J),
                "B_scaled_residual": _scaled_B(lam),
                "I1_scaled_residual": _scaled_I1(lam),
                "reduction_route_relative": _reduction_route_max(p),
            }
            for name, value in values.items():
                maxima[name] = max(maxima[name], float(value))
            point_pass = bool(
                values["sum_multiplier_absolute"] <= EXACT_IDENTITY_TOLERANCE
                and values["A_matrix_invariance_relative"] <= MATRIX_RELATIVE_TOLERANCE
                and values["A_beam_invariance_relative"] <= REDUCED_PROPERTY_TOLERANCE
                and values["D_matrix_formula_relative"] <= MATRIX_RELATIVE_TOLERANCE
                and values["D_beam_formula_relative"] <= REDUCED_PROPERTY_TOLERANCE
                and values["shear_matrix_invariance_relative"] <= MATRIX_RELATIVE_TOLERANCE
                and max(values["S_beam_invariance_relative"], values["mass_invariance_relative"], values["rotary_inertia_invariance_relative"]) <= REDUCED_PROPERTY_TOLERANCE
                and max(values["B_scaled_residual"], values["I1_scaled_residual"]) <= SYMMETRY_RELATIVE_TOLERANCE
                and values["reduction_route_relative"] <= REDUCTION_ROUTE_TOLERANCE
            )
            passed = passed and point_pass
            meta = TRANSFER_METADATA[configuration_id]
            rows.append({
                "configuration_id": configuration_id,
                "xi": xi,
                "xi_index": xi_index(xi),
                "transfer_vector_C_M_O": TRANSFER_VECTORS[configuration_id],
                "fixed_pair": meta["fixed_pair"],
                "donor_pair_for_positive_xi": meta["donor_pair_for_positive_xi"],
                "receiver_pair_for_positive_xi": meta["receiver_pair_for_positive_xi"],
                "s_C": section.multipliers[PAIR_CENTER],
                "s_M": section.multipliers[PAIR_MIDDLE],
                "s_O": section.multipliers[PAIR_OUTER],
                "sum_multipliers": sum(section.multipliers.values()),
                "D_key_q": D_key_q(configuration_id, xi),
                "A_beam": p.A, "A_over_A0": p.A / baseline.properties.A,
                "D_beam": p.D, "D_over_D0": p.D / baseline.properties.D,
                "expected_D_over_D0": expected,
                "S_beam": p.S, "S_over_S0": p.S / baseline.properties.S,
                "m": p.m, "m_over_m0": p.m / baseline.properties.m,
                "J": p.J, "J_over_J0": p.J / baseline.properties.J,
                "K": p.K, "width": p.width,
                "stack_bottom_to_top": STACK_BOTTOM_TO_TOP,
                "ply_angles_deg": [0.0] * 6,
                "ply_thicknesses": [PLY_THICKNESS] * 6,
                "z_interfaces": lam.z_interfaces,
                "A_matrix": lam.A, "B_matrix": lam.B, "D_matrix": lam.D,
                "shear_matrix_yz_xz": lam.shear,
                "I0": lam.I0, "I1": lam.I1, "I2": lam.I2,
                "B_relative": values["B_scaled_residual"],
                "I1_relative": values["I1_scaled_residual"],
                "A_invariance_relative": values["A_beam_invariance_relative"],
                "A_matrix_invariance_relative": values["A_matrix_invariance_relative"],
                "D_formula_relative": values["D_beam_formula_relative"],
                "D_matrix_formula_relative": values["D_matrix_formula_relative"],
                "shear_matrix_invariance_relative": values["shear_matrix_invariance_relative"],
                "S_invariance_relative": values["S_beam_invariance_relative"],
                "m_invariance_relative": values["mass_invariance_relative"],
                "J_invariance_relative": values["rotary_inertia_invariance_relative"],
                "reduction_route_max_relative": values["reduction_route_relative"],
                "constitutive_status": "PASS" if point_pass else "FAIL",
            })
    vectors_ok = (
        np.array_equal(
            np.asarray(TRANSFER_VECTORS[CENTER_OUTER_TRANSFER]),
            np.asarray(TRANSFER_VECTORS[CENTER_MIDDLE_TRANSFER]) + np.asarray(TRANSFER_VECTORS[MIDDLE_OUTER_TRANSFER]),
        )
        and [sum(w * v for w, v in zip((1, 7, 19), TRANSFER_VECTORS[c], strict=True)) for c in CONFIGURATIONS] == [6, 12, 18]
    )
    passed = bool(passed and vectors_ok)
    return {
        "status": "PASS" if passed else "FAIL",
        "passed": passed,
        "row_count": len(rows),
        "maximum_residuals": maxima,
        "baseline_beam_properties": {field: getattr(baseline.properties, field) for field in ("A", "D", "S", "m", "J", "K", "width")},
        "pair_A_weights": {key: float(value) for key, value in PAIR_A_WEIGHTS.items()},
        "pair_D_weights": {key: float(value) for key, value in PAIR_D_WEIGHTS.items()},
        "D_projections": {key: sum(w * v for w, v in zip((1, 7, 19), TRANSFER_VECTORS[key], strict=True)) for key in CONFIGURATIONS},
        "neutral_RLB2I_direction": [2, -3, 1],
        "vector_sum_identity": vectors_ok,
        "thresholds": contract_payload()["thresholds"],
    }, rows


def make_matrix_provider(configuration_id: str, xi: float) -> tuple[Any, dict[str, Any]]:
    """Build the frozen coupled matrix with one cached identical-arm map."""

    global MATRIX_ASSEMBLY_COUNT
    _beam, coupled = _physics_modules()
    section = build_six_ply_section(configuration_id, xi)
    joint = np.asarray(coupled.joint_matrix(math.radians(BETA_DEG)), dtype=float)

    def provider(omega: float) -> FloatArray:
        global MATRIX_ASSEMBLY_COUNT
        MATRIX_ASSEMBLY_COUNT += 1
        arm_map = np.asarray(
            coupled.arm_clamp_map(float(omega), L1, section.properties), dtype=float
        )
        combined = np.zeros((12, 6), dtype=float)
        combined[:6, :3] = arm_map
        combined[6:, 3:] = arm_map
        return np.asarray(joint @ combined, dtype=float)

    public_residual = 0.0
    for probe in (0.731, 3.217):
        direct = np.asarray(
            coupled.coupled_boundary_matrix(
                probe,
                math.radians(BETA_DEG),
                L1,
                section.properties,
                L2,
                section.properties,
            ),
            dtype=float,
        )
        public_residual = max(public_residual, float(np.max(np.abs(provider(probe) - direct))))
    if public_residual > 16.0 * np.finfo(float).eps:
        raise RuntimeError("RLB-2J cached provider differs from the public coupled builder.")
    return provider, {
        "configuration_id": configuration_id,
        "xi": float(xi),
        "beta_deg": BETA_DEG,
        "identical_arms": True,
        "identical_arm_map_reused": True,
        "cached_vs_public_builder_max_abs": public_residual,
        "production_modules_only": True,
    }


def Omega_to_Lambda(Omega: float) -> float:
    value = float(Omega)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError("Omega must be finite and nonnegative.")
    return math.sqrt(value)


def continuation_leg(xi: float) -> str:
    if xi < 0.0:
        return "NEGATIVE"
    if xi > 0.0:
        return "POSITIVE"
    return "ZERO_ANCHOR"


def _root_rows(
    configuration_id: str,
    xi: float,
    slots: Sequence[Any],
    *,
    physical_solve_id: str,
    transaction_id: str,
    solve_mode: str,
    fallback_used: bool,
    predicted: Sequence[float] | None,
    search_right: float,
    unresolved: int,
    grid_kind: str = "BASE",
    repair_id: str = "",
    shared_anchor_reused: bool = False,
    shared_anchor_source: str = "",
) -> tuple[dict[str, Any], ...]:
    base = _base()
    windows = None if predicted is None else base.local_search_windows(predicted)
    multipliers = stiffness_multipliers(configuration_id, xi)
    guard = float(slots[K_GUARD - 1].event.omega_bar)
    rows: list[dict[str, Any]] = []
    for position, slot in enumerate(slots[:K_GUARD], start=1):
        event = slot.event
        candidate = event.candidate
        diagnostic = candidate.diagnostics
        locator = (
            (float(candidate.interval_left_bar), float(candidate.interval_right_bar))
            if windows is None
            else windows[position - 1]
        )
        rows.append({
            "row_id": f"{configuration_id}__{xi:+.6f}__{grid_kind}__p{position:02d}__{repair_id or 'base'}",
            "configuration_id": configuration_id,
            "xi": float(xi),
            "xi_index": xi_index(xi),
            "s_C": multipliers[PAIR_CENTER],
            "s_M": multipliers[PAIR_MIDDLE],
            "s_O": multipliers[PAIR_OUTER],
            "D_over_D0": expected_D_ratio(configuration_id, xi),
            "D_key_q": D_key_q(configuration_id, xi),
            "grid_kind": grid_kind,
            "physical_solve_id": physical_solve_id,
            "transaction_id": transaction_id,
            "continuation_leg": continuation_leg(xi),
            "sorted_position": position,
            "root_role": "PLOTTED" if position <= K_PLOT else "ROOT_9_GUARD",
            "guard_flag": position == K_GUARD,
            "omega": float(event.omega),
            "Omega": float(event.omega_bar),
            "Lambda": Omega_to_Lambda(float(event.omega_bar)),
            "predictor_Omega": "" if predicted is None else float(predicted[position - 1]),
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
            "cluster_multiplicity": event.cluster_multiplicity,
            "cluster_total_nullity": event.cluster_total_nullity,
            "unresolved_candidates_below_root9": unresolved,
            "search_right_Omega": search_right,
            "root9_right_margin_Omega": search_right - guard,
            "solve_mode": solve_mode,
            "fallback_used": fallback_used,
            "quality_status": "PASS",
            "point_status": "PASS",
            "shared_xi0_anchor_reused": shared_anchor_reused,
            "shared_xi0_source_configuration": shared_anchor_source,
            "is_canonical_plot_source": True,
            "supersedes_row_id": "",
            "repair_id": repair_id,
            "roots_above_9_computed": False,
        })
    return tuple(rows)


def solve_point(
    configuration_id: str,
    xi: float,
    *,
    previous: tuple[float, Sequence[float]] | None = None,
    second_previous: tuple[float, Sequence[float]] | None = None,
    force_anchor: bool = False,
    dense_local: bool = False,
    dense_positions: Sequence[int] | None = None,
    grid_kind: str = "BASE",
    repair_id: str = "",
) -> PointSolution:
    """Compute one nine-root transaction; predictors are locators only."""

    global ROOT_CALCULATION_COUNT
    ROOT_CALCULATION_COUNT += 1
    started = time.perf_counter()
    base = _base()
    from scripts.analysis.laminated_beams import (
        sweep_reddy_axial_stiffness_visibility as rlb2h,
    )
    physical_solve_id = f"{configuration_id}__xi_{float(xi):+.6f}"
    transaction_id = hashlib.sha256(
        f"{STAGE_ID}|{physical_solve_id}|{grid_kind}|{repair_id}".encode("utf-8")
    ).hexdigest()[:20].upper()
    provider, _metadata = make_matrix_provider(configuration_id, xi)
    counted = base.CountedProvider(provider)
    policy = base._rlb2e_search_policy()
    predicted: FloatArray | None = None
    fallback_used = False
    if not force_anchor and previous is not None:
        predicted = base.hold_secant_predictor(
            xi,
            previous[0],
            previous[1],
            None if second_previous is None else second_previous[0],
            None if second_previous is None else second_previous[1],
        )
        predicted = np.sort(predicted)
        try:
            canonical, rejected, slots, search_right, refinements = base._local_candidates(
                counted,
                policy,
                predicted,
                solve_id=physical_solve_id,
                dense=dense_local,
                dense_positions=dense_positions,
            )
            canonical, slots = rlb2h._truncate_inventory_to_root9(canonical, slots, policy)
            accepted, quality = rlb2h._point_is_acceptable_with_multiplicity(
                canonical, rejected, slots, search_right, policy
            )
        except (ValueError, RuntimeError, ArithmeticError, np.linalg.LinAlgError):
            accepted = False
            quality = {}
        if not accepted:
            fallback_used = True
            canonical, rejected, slots, search_right, refinements = base._progressive_anchor_candidates(
                counted, policy, solve_id=physical_solve_id + "_fallback"
            )
            canonical, slots = rlb2h._truncate_inventory_to_root9(canonical, slots, policy)
            accepted, quality = rlb2h._point_is_acceptable_with_multiplicity(
                canonical, rejected, slots, search_right, policy
            )
            solve_mode = "BOUNDED_GLOBAL_RECOVERY"
        else:
            solve_mode = "FAST_LOCAL"
    else:
        canonical, rejected, slots, search_right, refinements = base._progressive_anchor_candidates(
            counted, policy, solve_id=physical_solve_id
        )
        canonical, slots = rlb2h._truncate_inventory_to_root9(canonical, slots, policy)
        accepted, quality = rlb2h._point_is_acceptable_with_multiplicity(
            canonical, rejected, slots, search_right, policy
        )
        solve_mode = "PROGRESSIVE_ANCHOR"
    if not accepted:
        raise RuntimeError(
            f"{physical_solve_id}: first-eight-plus-root9 quality failed: {quality}"
        )
    rows = _root_rows(
        configuration_id,
        xi,
        slots,
        physical_solve_id=physical_solve_id,
        transaction_id=transaction_id,
        solve_mode=solve_mode,
        fallback_used=fallback_used,
        predicted=None if fallback_used else predicted,
        search_right=search_right,
        unresolved=int(quality["unresolved_candidates_below_root9"]),
        grid_kind=grid_kind,
        repair_id=repair_id,
    )
    return PointSolution(
        configuration_id=configuration_id,
        xi=float(xi),
        rows=rows,
        wall_time_seconds=time.perf_counter() - started,
        peak_rss_bytes=_peak_rss_bytes(),
        determinant_evaluations=counted.evaluations,
        sigma_evaluations=counted.evaluations,
        search_left_Omega=(
            1.0e-8
            if predicted is None or fallback_used
            else base.local_search_windows(predicted)[0][0]
        ),
        search_right_Omega=search_right,
        local_refinements=refinements,
        solve_mode=solve_mode,
        fallback_used=fallback_used,
        unresolved_candidates_below_root9=int(quality["unresolved_candidates_below_root9"]),
        candidate_count_total=int(quality["candidate_count_total"]),
        accepted_candidates_above_root9=int(quality["accepted_candidates_above_root9"]),
        retained_slots_above_root9=int(quality["retained_slots_above_root9"]),
        roots_above_9_computed=bool(quality["roots_above_9_computed"]),
    )


def _as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _base_group_index(rows: Sequence[Mapping[str, Any]]) -> dict[tuple[str, int], list[Mapping[str, Any]]]:
    groups: dict[tuple[str, int], list[Mapping[str, Any]]] = {}
    for row in rows:
        if str(row.get("grid_kind")) != "BASE":
            continue
        key = (str(row["configuration_id"]), int(row["xi_index"]))
        groups.setdefault(key, []).append(row)
    return groups


def _base_group_is_acceptable(group: Sequence[Mapping[str, Any]]) -> bool:
    try:
        ordered = sorted(group, key=lambda row: int(row["sorted_position"]))
        positions = [int(row["sorted_position"]) for row in ordered]
        Omegas = np.asarray([float(row["Omega"]) for row in ordered], dtype=float)
        multiplicities = [int(row.get("cluster_multiplicity", 1)) for row in ordered]
        cluster_ids = [str(row.get("cluster_id", "")) for row in ordered]
        order_ok = all(
            right > left
            or (
                abs(right - left) <= 5.0e-10 + 5.0e-12 * max(abs(left), abs(right))
                and cluster_ids[index] == cluster_ids[index + 1]
                and multiplicities[index] >= 2
                and multiplicities[index + 1] >= 2
            )
            for index, (left, right) in enumerate(zip(Omegas[:-1], Omegas[1:], strict=True))
        )
        return bool(
            len(group) == K_GUARD
            and positions == list(range(1, K_GUARD + 1))
            and np.all(np.isfinite(Omegas))
            and np.all(Omegas > 0.0)
            and order_ok
            and all(str(row["quality_status"]) == "PASS" for row in ordered)
            and all(int(row["unresolved_candidates_below_root9"]) == 0 for row in ordered)
            and all(float(row["scaled_sigma_ratio"]) <= ROOT_SINGULAR_RATIO_TOLERANCE for row in ordered)
            and all(float(row["boundary_null_residual"]) <= BOUNDARY_RESIDUAL_TOLERANCE for row in ordered)
            and all(not _as_bool(row["predictor_used_as_final"]) for row in ordered)
            and all(not _as_bool(row["roots_above_9_computed"]) for row in ordered)
            and float(ordered[-1]["root9_right_margin_Omega"]) >= ROOT9_RIGHT_TAIL_OMEGA - 1.0e-10
        )
    except (KeyError, TypeError, ValueError, OverflowError):
        return False


def _complete_base_groups(rows: Sequence[Mapping[str, Any]]) -> dict[tuple[str, int], list[Mapping[str, Any]]]:
    return {key: group for key, group in _base_group_index(rows).items() if _base_group_is_acceptable(group)}


def audit_spectrum_rows(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    groups = _base_group_index(rows)
    complete = _complete_base_groups(rows)
    expected = {(configuration, index) for configuration in CONFIGURATIONS for index in range(-40, 41)}
    duplicates = []
    invalid = []
    for key, group in groups.items():
        positions = [int(row["sorted_position"]) for row in group]
        if len(positions) != len(set(positions)):
            duplicates.append(f"{key[0]}:{key[1]}")
        if not _base_group_is_acceptable(group):
            invalid.append(f"{key[0]}:{key[1]}")
    base_rows = [row for row in rows if str(row.get("grid_kind")) == "BASE"]
    row_ids = [str(row["row_id"]) for row in rows]
    roots_above = [row for row in rows if int(row["sorted_position"]) > K_GUARD]
    missing = sorted(expected - set(groups), key=lambda item: (CONFIGURATIONS.index(item[0]), item[1]))
    extra = sorted(set(groups) - expected)
    status = "PASS" if not duplicates and not invalid and not missing and not extra and len(row_ids) == len(set(row_ids)) and not roots_above and len(base_rows) == 2187 else "FAIL"
    return {
        "status": status,
        "base_group_count": len(complete),
        "base_row_count": len(base_rows),
        "missing_groups": [f"{c}:{i}" for c, i in missing],
        "extra_groups": [f"{c}:{i}" for c, i in extra],
        "duplicate_groups": duplicates,
        "invalid_groups": invalid,
        "duplicate_row_id_count": len(row_ids) - len(set(row_ids)),
        "roots_above_guard_count": len(roots_above),
    }


def _rows_for_roots(rows: Sequence[Mapping[str, Any]], configuration_id: str, xi: float) -> FloatArray:
    group = _complete_base_groups(rows).get((configuration_id, xi_index(xi)), [])
    ordered = sorted(group, key=lambda row: int(row["sorted_position"]))
    if len(ordered) != K_GUARD:
        raise RuntimeError(f"Incomplete root group {configuration_id}, xi={xi}.")
    return np.asarray([float(row["Omega"]) for row in ordered], dtype=float)


def _write_point_transaction(
    output_dir: Path, existing_rows: Sequence[Mapping[str, Any]], solution: PointSolution
) -> list[dict[str, Any]]:
    rows = [dict(row) for row in existing_rows]
    key = (solution.configuration_id, xi_index(solution.xi))
    current = _base_group_index(rows).get(key, [])
    if current and _base_group_is_acceptable(current):
        return rows
    rows = [
        row for row in rows
        if not (str(row["configuration_id"]) == key[0] and int(row["xi_index"]) == key[1] and str(row["grid_kind"]) == "BASE")
    ]
    rows.extend(dict(row) for row in solution.rows)
    rows.sort(key=lambda row: (CONFIGURATIONS.index(str(row["configuration_id"])), int(row["xi_index"]), 0 if str(row["grid_kind"]) == "BASE" else 1, int(row["sorted_position"])))
    _atomic_write_csv(Path(output_dir) / SPECTRUM_FILENAME, rows, SPECTRUM_FIELDS)
    return rows


def _solution_record(solution: PointSolution, *, benchmark: bool) -> dict[str, Any]:
    return {
        "configuration_id": solution.configuration_id,
        "xi": solution.xi,
        "benchmark": benchmark,
        "wall_time_seconds": solution.wall_time_seconds,
        "peak_rss_bytes": solution.peak_rss_bytes,
        "determinant_evaluations": solution.determinant_evaluations,
        "sigma_evaluations": solution.sigma_evaluations,
        "search_left_Omega": solution.search_left_Omega,
        "search_right_Omega": solution.search_right_Omega,
        "local_refinements": solution.local_refinements,
        "solve_mode": solution.solve_mode,
        "fallback_used": solution.fallback_used,
        "roots_above_9_computed": solution.roots_above_9_computed,
        "roots": [{
            "sorted_position": int(row["sorted_position"]),
            "Omega": float(row["Omega"]),
            "Lambda": float(row["Lambda"]),
            "scaled_sigma_ratio": float(row["scaled_sigma_ratio"]),
            "boundary_null_residual": float(row["boundary_null_residual"]),
        } for row in solution.rows],
    }


def contract_payload() -> dict[str, Any]:
    return {
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "frequency_map_policy": FREQUENCY_MAP_POLICY,
        "geometry": {
            "mu": MU, "tau": TAU, "beta_deg": BETA_DEG,
            "l": L_REFERENCE, "L1": L1, "L2": L2, "L_total": L_TOTAL,
            "width": WIDTH, "thickness": THICKNESS,
            "ply_thickness": PLY_THICKNESS, "K": K,
            "outer_clamps": True, "joint": "frozen ideal rigid RLB joint",
            "identical_arms": True,
        },
        "material_M0": base_material_contract(),
        "six_ply_stack": {
            "stack_bottom_to_top": STACK_BOTTOM_TO_TOP,
            "pair_order_midplane_to_surface": PAIR_ORDER,
            "ply_angles_deg": [0.0] * 6,
            "ply_thicknesses": [PLY_THICKNESS] * 6,
            "scaled_fields": ["E1", "E2", "G12"],
            "fixed_fields": ["nu12", "G13", "G23", "rho"],
        },
        "transfer_definitions": {
            key: {
                "vector_C_M_O": TRANSFER_VECTORS[key],
                "lever_index": TRANSFER_LEVERS[key],
                "D_over_D0": f"1+{TRANSFER_LEVERS[key]}*xi/4.5",
                **TRANSFER_METADATA[key],
            }
            for key in CONFIGURATIONS
        },
        "pair_weights": {
            "A_C_M_O": [1.0 / 3.0] * 3,
            "D_C_M_O": [1.0 / 27.0, 7.0 / 27.0, 19.0 / 27.0],
            "D_integer_weights_C_M_O": [1, 7, 19],
            "pairwise_D_projections": [6, 12, 18],
            "bending_lever_ratio": [1, 2, 3],
            "RLB2I_neutral_direction_C_M_O": [2, -3, 1],
        },
        "normalization": {
            "Omega": "omega*l^2*sqrt(rho0*A0/(E0*Iy0))",
            "Lambda": "sqrt(Omega)",
            "E0": 1.0, "rho0": RHO_0, "b0": WIDTH,
            "h0": THICKNESS, "l": L_REFERENCE,
            "Omega_per_omega": OMEGA_TO_OMEGA_SCALE,
        },
        "xi_grid": [float(value) for value in xi_grid()],
        "exact_D_grid_key": "xi=n/50; q=lever_index*n; D/D0=1+q/225",
        "thresholds": {
            "exact_identity_absolute": EXACT_IDENTITY_TOLERANCE,
            "matrix_relative": MATRIX_RELATIVE_TOLERANCE,
            "symmetry_relative": SYMMETRY_RELATIVE_TOLERANCE,
            "reduced_property_relative": REDUCED_PROPERTY_TOLERANCE,
            "reduction_route_relative": REDUCTION_ROUTE_TOLERANCE,
            "matched_D_spectral_relative": SPECTRAL_RELATIVE_TOLERANCE,
            "root_singular_ratio": ROOT_SINGULAR_RATIO_TOLERANCE,
            "boundary_residual": BOUNDARY_RESIDUAL_TOLERANCE,
            "slope_diagnostic_relative": SLOPE_DIAGNOSTIC_RELATIVE,
            "slope_insensitive_absolute": SLOPE_INSENSITIVE_ABSOLUTE,
            "neighbour_MAD_multiplier": NEIGHBOUR_MAD_MULTIPLIER,
            "neighbour_absolute_trigger": NEIGHBOUR_ABSOLUTE_TRIGGER,
        },
        "root_contract": {
            "plotted_positions": list(range(1, K_PLOT + 1)),
            "guard_position": K_GUARD,
            "guard_role": "completeness_only",
            "roots_above_9": "not_computed",
        },
        "explicit_exclusions": [
            "roots_10_and_above", "branch_tracking", "MAC", "mode_shapes",
            "energy_or_stress_reconstruction", "Ritz", "FEM", "smoothing",
            "spectral_sweep_runner", "certified_audit", "arm_asymmetry",
            "beta_mu_tau_sweeps", "commit", "push",
        ],
    }


def contract_hash() -> str:
    encoded = json.dumps(
        _json_value(contract_payload()), sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest().upper()


def _checkpoint_payload(
    rows: Sequence[Mapping[str, Any]],
    point_records: Sequence[Mapping[str, Any]],
    *,
    constitutive: Mapping[str, Any],
    started_at: str,
    benchmark_status: str,
) -> dict[str, Any]:
    groups = _complete_base_groups(rows)
    expected = [(configuration, index) for configuration in CONFIGURATIONS for index in range(-40, 41)]
    missing = [
        {"configuration_id": configuration, "xi": index / 50.0}
        for configuration, index in expected
        if (configuration, index) not in groups
    ]
    return {
        "schema_version": 1,
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "contract_sha256": contract_hash(),
        "started_at_utc": started_at,
        "updated_at_utc": _utc_now(),
        "constitutive_gate": constitutive["status"],
        "benchmark_status": benchmark_status,
        "completed_base_groups": len(groups),
        "completed_base_rows": sum(len(group) for group in groups.values()),
        "missing_points": missing,
        "failed_points": [],
        "last_completed_parameter": point_records[-1].get("xi") if point_records else None,
        "point_records": list(point_records),
        "invocation_root_calculation_count": ROOT_CALCULATION_COUNT,
        "parallel_workers_used": 0,
        "thread_limits": {name: os.environ.get(name, "") for name in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
    }


def _write_checkpoint(
    target: Path,
    rows: Sequence[Mapping[str, Any]],
    point_records: Sequence[Mapping[str, Any]],
    constitutive: Mapping[str, Any],
    started_at: str,
    benchmark_status: str,
) -> None:
    _atomic_write_json(
        target / CHECKPOINT_FILENAME,
        _checkpoint_payload(
            rows,
            point_records,
            constitutive=constitutive,
            started_at=started_at,
            benchmark_status=benchmark_status,
        ),
    )


def _clone_shared_xi0(
    source: PointSolution, configuration_id: str
) -> PointSolution:
    cloned_rows: list[dict[str, Any]] = []
    multipliers = stiffness_multipliers(configuration_id, 0.0)
    for row in source.rows:
        clone = dict(row)
        clone.update({
            "row_id": str(clone["row_id"]).replace(source.configuration_id, configuration_id, 1),
            "configuration_id": configuration_id,
            "s_C": multipliers[PAIR_CENTER], "s_M": multipliers[PAIR_MIDDLE], "s_O": multipliers[PAIR_OUTER],
            "D_over_D0": 1.0, "D_key_q": 0,
            "shared_xi0_anchor_reused": True,
            "shared_xi0_source_configuration": source.configuration_id,
            "solve_mode": "SHARED_XI0_REUSE",
        })
        cloned_rows.append(clone)
    return replace(
        source,
        configuration_id=configuration_id,
        rows=tuple(cloned_rows),
        wall_time_seconds=0.0,
        determinant_evaluations=0,
        sigma_evaluations=0,
        solve_mode="SHARED_XI0_REUSE",
    )


def _benchmark_payload(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    endpoint_times = [
        float(item["wall_time_seconds"])
        for item in records
        if str(item.get("configuration_id")) == CENTER_OUTER_TRANSFER
        and abs(float(item.get("xi", 0.0))) == 0.8
        and int(item.get("determinant_evaluations", 0)) > 0
    ]
    remaining = 241 - 3
    maximum = max(endpoint_times, default=math.inf)
    eta = 1.5 * maximum * remaining
    return {
        "schema_version": 1,
        "stage_id": STAGE_ID,
        "anchors": [dict(item) for item in records],
        "anchor_count": 3,
        "shared_xi0_anchor_calculated_once": True,
        "total_unique_physical_base_solves": 241,
        "remaining_unique_root_points": remaining,
        "measured_endpoint_max_seconds": maximum,
        "eta_formula": "1.5*max(CENTER_OUTER endpoint anchor times)*238",
        "conservative_eta_seconds": eta,
        "eta_limit_seconds": ETA_LIMIT_SECONDS,
        "production_run_permitted": bool(eta <= ETA_LIMIT_SECONDS),
    }


def run_benchmarks(
    target: Path,
    rows: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    started_at: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    groups = _complete_base_groups(rows)
    existing_benchmark_path = target / BENCHMARK_FILENAME
    preserved = (
        json.loads(existing_benchmark_path.read_text(encoding="utf-8"))
        if existing_benchmark_path.is_file()
        else {}
    )
    index = {
        (str(item.get("configuration_id")), xi_index(float(item.get("xi", 0.0)))): dict(item)
        for item in [*point_records, *preserved.get("anchors", [])]
        if bool(item.get("benchmark", False))
    }
    benchmark_records: list[dict[str, Any]] = []
    source_key = (CENTER_OUTER_TRANSFER, 0)
    if source_key not in groups:
        source = solve_point(CENTER_OUTER_TRANSFER, 0.0, force_anchor=True)
        rows = _write_point_transaction(target, rows, source)
        record = _solution_record(source, benchmark=True)
        point_records.append(record)
        benchmark_records.append(record)
    else:
        if source_key not in index:
            raise RuntimeError("Shared xi=0 spectrum exists without benchmark provenance.")
        benchmark_records.append(index[source_key])
        source_rows = sorted(groups[source_key], key=lambda row: int(row["sorted_position"]))
        source = PointSolution(
            CENTER_OUTER_TRANSFER, 0.0, tuple(dict(row) for row in source_rows),
            0.0, 0, 0, 0, 1.0e-8, float(source_rows[-1]["search_right_Omega"]),
            0, "SHARED_XI0_SOURCE_REUSE", False, 0,
        )
    groups = _complete_base_groups(rows)
    for configuration_id in (CENTER_MIDDLE_TRANSFER, MIDDLE_OUTER_TRANSFER):
        key = (configuration_id, 0)
        if key in groups:
            continue
        clone = _clone_shared_xi0(source, configuration_id)
        rows = _write_point_transaction(target, rows, clone)
        point_records.append(_solution_record(clone, benchmark=True))
        groups = _complete_base_groups(rows)
    _write_checkpoint(target, rows, point_records, constitutive, started_at, "RUNNING")

    for endpoint in (-0.8, 0.8):
        key = (CENTER_OUTER_TRANSFER, xi_index(endpoint))
        groups = _complete_base_groups(rows)
        if key in groups:
            if key not in index:
                raise RuntimeError(f"Benchmark endpoint xi={endpoint} lacks durable metrics.")
            benchmark_records.append(index[key])
            continue
        solution = solve_point(CENTER_OUTER_TRANSFER, endpoint, force_anchor=True)
        rows = _write_point_transaction(target, rows, solution)
        record = _solution_record(solution, benchmark=True)
        point_records.append(record)
        benchmark_records.append(record)
        _atomic_write_json(target / BENCHMARK_FILENAME, _benchmark_payload(benchmark_records))
        _write_checkpoint(target, rows, point_records, constitutive, started_at, "RUNNING")
    benchmark = _benchmark_payload(benchmark_records)
    _atomic_write_json(target / BENCHMARK_FILENAME, benchmark)
    _write_checkpoint(
        target,
        rows,
        point_records,
        constitutive,
        started_at,
        "PASS" if benchmark["production_run_permitted"] else "STOPPED_BY_ETA_GATE",
    )
    return rows, benchmark


def _complete_leg(
    target: Path,
    rows: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    *,
    configuration_id: str,
    values: Sequence[float],
    constitutive: Mapping[str, Any],
    started_at: str,
) -> list[dict[str, Any]]:
    previous: tuple[float, Sequence[float]] | None = None
    second: tuple[float, Sequence[float]] | None = None
    groups = _complete_base_groups(rows)
    for raw_value in values:
        xi = float(raw_value)
        key = (configuration_id, xi_index(xi))
        if key in groups:
            roots = _rows_for_roots(rows, configuration_id, xi)
            second = previous
            previous = (xi, roots)
            continue
        solution = solve_point(
            configuration_id,
            xi,
            previous=previous,
            second_previous=second,
            force_anchor=previous is None,
        )
        rows = _write_point_transaction(target, rows, solution)
        point_records.append(_solution_record(solution, benchmark=False))
        groups = _complete_base_groups(rows)
        second = previous
        previous = (xi, _rows_for_roots(rows, configuration_id, xi))
        _write_checkpoint(target, rows, point_records, constitutive, started_at, "PASS")
    return rows


def complete_missing_points(
    target: Path,
    rows: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    started_at: str,
) -> list[dict[str, Any]]:
    negative = [0.0, *[index / 50.0 for index in range(-1, -41, -1)]]
    positive = [0.0, *[index / 50.0 for index in range(1, 41)]]
    for configuration_id in CONFIGURATIONS:
        rows = _complete_leg(
            target, rows, point_records,
            configuration_id=configuration_id, values=negative,
            constitutive=constitutive, started_at=started_at,
        )
        rows = _complete_leg(
            target, rows, point_records,
            configuration_id=configuration_id, values=positive,
            constitutive=constitutive, started_at=started_at,
        )
    return rows


def centred_secant_residual(left: float, center: float, right: float) -> float:
    predictor = 0.5 * (float(left) + float(right))
    return abs(float(center) - predictor) / max(
        abs(float(center)), abs(predictor), np.finfo(float).tiny
    )


def _canonical_root_group(
    rows: Sequence[Mapping[str, Any]], configuration_id: str, xi: float
) -> list[Mapping[str, Any]]:
    selected = [
        row for row in rows
        if str(row["configuration_id"]) == configuration_id
        and int(row["xi_index"]) == xi_index(xi)
        and _as_bool(row.get("is_canonical_plot_source", True))
    ]
    selected.sort(key=lambda row: int(row["sorted_position"]))
    if [int(row["sorted_position"]) for row in selected] != list(range(1, K_GUARD + 1)):
        raise RuntimeError(f"Canonical root group incomplete: {configuration_id}, xi={xi}.")
    return selected


def neighbour_audit_rows(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Apply the unchanged RLB-2E median+8*MAD local anomaly trigger."""

    grid = [index / 50.0 for index in range(-40, 41)]
    audit: list[dict[str, Any]] = []
    gap_flags: set[tuple[str, int, int]] = set()
    for configuration_id in CONFIGURATIONS:
        spectra = {
            xi_index(xi): np.asarray(
                [float(row["Lambda"]) for row in _canonical_root_group(rows, configuration_id, xi)[:K_PLOT]],
                dtype=float,
            )
            for xi in grid
        }
        for lower in range(K_PLOT - 1):
            gaps = np.asarray([spectra[index][lower + 1] - spectra[index][lower] for index in range(-40, 41)])
            residuals = np.asarray([
                centred_secant_residual(gaps[offset - 1], gaps[offset], gaps[offset + 1])
                for offset in range(1, 80)
            ])
            median = float(np.median(residuals))
            mad = float(np.median(np.abs(residuals - median)))
            threshold = median + NEIGHBOUR_MAD_MULTIPLIER * mad
            for offset, index in enumerate(range(-39, 40)):
                if float(residuals[offset]) > threshold and float(residuals[offset]) > NEIGHBOUR_ABSOLUTE_TRIGGER:
                    gap_flags.add((configuration_id, index, lower + 1))
                    gap_flags.add((configuration_id, index, lower + 2))
        for position in range(1, K_PLOT + 1):
            values = {index: spectra[index][position - 1] for index in range(-40, 41)}
            residuals = [
                centred_secant_residual(values[index - 1], values[index], values[index + 1])
                for index in range(-39, 40)
            ]
            median = float(np.median(residuals))
            mad = float(np.median(np.abs(np.asarray(residuals) - median)))
            threshold = median + NEIGHBOUR_MAD_MULTIPLIER * mad
            for offset, index in enumerate(range(-39, 40)):
                group = _canonical_root_group(rows, configuration_id, index / 50.0)
                row = group[position - 1]
                residual = float(residuals[offset])
                statistical = residual > threshold and residual > NEIGHBOUR_ABSOLUTE_TRIGGER
                bad_quality = bool(
                    float(row["scaled_sigma_ratio"]) > ROOT_SINGULAR_RATIO_TOLERANCE
                    or float(row["boundary_null_residual"]) > BOUNDARY_RESIDUAL_TOLERANCE
                    or int(row["unresolved_candidates_below_root9"]) != 0
                )
                gap = (configuration_id, index, position) in gap_flags
                flagged = bool(statistical or bad_quality or gap)
                audit.append({
                    "configuration_id": configuration_id,
                    "sorted_position": position,
                    "xi_left": (index - 1) / 50.0,
                    "xi": index / 50.0,
                    "xi_right": (index + 1) / 50.0,
                    "Lambda_left": values[index - 1],
                    "Lambda_center": values[index],
                    "Lambda_right": values[index + 1],
                    "centered_predictor_Lambda": 0.5 * (values[index - 1] + values[index + 1]),
                    "centered_secant_residual": residual,
                    "median_residual_for_position": median,
                    "MAD_residual_for_position": mad,
                    "robust_threshold": threshold,
                    "absolute_trigger": NEIGHBOUR_ABSOLUTE_TRIGGER,
                    "statistical_flag": statistical,
                    "bad_residual_or_unresolved_warning": bad_quality,
                    "gap_jump_warning": gap,
                    "flagged": flagged,
                    "repair_id": "",
                    "repair_status": "PENDING" if flagged else "NOT_REQUIRED",
                    "smoothing_applied": False,
                })
    return audit


def flagged_repair_points(audit_rows: Sequence[Mapping[str, Any]]) -> list[tuple[str, float]]:
    return sorted(
        {
            (str(row["configuration_id"]), float(row["xi"]))
            for row in audit_rows if bool(row["flagged"])
        },
        key=lambda item: (CONFIGURATIONS.index(item[0]), item[1]),
    )


def apply_local_repairs(
    rows: list[dict[str, Any]], audit_rows: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Recompute only flagged points; never smooth or interpolate."""

    repair_records: list[dict[str, Any]] = []
    for repair_number, (configuration_id, xi) in enumerate(flagged_repair_points(audit_rows), start=1):
        repair_id = f"repair_{repair_number:04d}"
        positions = sorted({
            int(row["sorted_position"]) for row in audit_rows
            if str(row["configuration_id"]) == configuration_id
            and xi_index(float(row["xi"])) == xi_index(xi)
            and bool(row["flagged"])
        })
        original_group = _canonical_root_group(rows, configuration_id, xi)
        original = np.asarray([float(row["Omega"]) for row in original_group])
        existing_local = sorted(
            (
                row for row in rows
                if str(row["configuration_id"]) == configuration_id
                and int(row["xi_index"]) == xi_index(xi)
                and str(row["grid_kind"]) == "LOCAL_REFINEMENT"
                and str(row.get("repair_id", "")) == repair_id
            ),
            key=lambda row: int(row["sorted_position"]),
        )
        if len(existing_local) == K_GUARD and [int(row["sorted_position"]) for row in existing_local] == list(range(1, K_GUARD + 1)):
            refined = np.asarray([float(row["Omega"]) for row in existing_local])
            relative = float(np.max(np.abs(original - refined) / np.maximum.reduce((np.abs(original), np.abs(refined), np.full(K_GUARD, np.finfo(float).tiny)))))
            status = "LOCATOR_CORRECTED" if any(str(row["point_status"]) == "LOCATOR_CORRECTED_AFTER_LOCAL_REPAIR" for row in existing_local) else "REPRODUCED"
            for audit_row in audit_rows:
                if str(audit_row["configuration_id"]) == configuration_id and xi_index(float(audit_row["xi"])) == xi_index(xi) and bool(audit_row["flagged"]):
                    audit_row["repair_id"] = repair_id
                    audit_row["repair_status"] = "LOCATOR_CORRECTED_AFTER_LOCAL_REPAIR" if status == "LOCATOR_CORRECTED" else "REPRODUCED_AFTER_LOCAL_REPAIR"
            repair_records.append({
                "repair_id": repair_id,
                "configuration_id": configuration_id,
                "xi": xi,
                "affected_positions": positions,
                "status": status,
                "maximum_relative_Omega_difference": relative,
                "tolerance": SPECTRAL_RELATIVE_TOLERANCE,
                "wall_time_seconds": 0.0,
                "peak_rss_bytes": 0,
                "determinant_evaluations": 0,
                "sigma_evaluations": 0,
                "error": "",
                "reused_existing_local_repair": True,
                "smoothing_applied": False,
            })
            continue
        try:
            solution = solve_point(
                configuration_id,
                xi,
                previous=(xi, original),
                force_anchor=False,
                dense_local=True,
                dense_positions=positions,
                grid_kind="LOCAL_REFINEMENT",
                repair_id=repair_id,
            )
        except Exception as exc:
            status = "UNRESOLVED"
            relative = math.inf
            solution = None
            error = f"{type(exc).__name__}: {exc}"
        else:
            refined = np.asarray([float(row["Omega"]) for row in solution.rows])
            relative = float(np.max(np.abs(original - refined) / np.maximum.reduce((np.abs(original), np.abs(refined), np.full(K_GUARD, np.finfo(float).tiny)))))
            status = "REPRODUCED" if relative <= SPECTRAL_RELATIVE_TOLERANCE else "LOCATOR_CORRECTED"
            error = ""
            local_rows = [dict(row) for row in solution.rows]
            for row in local_rows:
                row["repair_id"] = repair_id
                row["grid_kind"] = "LOCAL_REFINEMENT"
                row["point_status"] = (
                    "REPRODUCED_AFTER_LOCAL_REPAIR"
                    if status == "REPRODUCED"
                    else "LOCATOR_CORRECTED_AFTER_LOCAL_REPAIR"
                )
                row["is_canonical_plot_source"] = status == "LOCATOR_CORRECTED"
                if status == "LOCATOR_CORRECTED":
                    matching = next(item for item in original_group if int(item["sorted_position"]) == int(row["sorted_position"]))
                    row["supersedes_row_id"] = matching["row_id"]
            if status == "LOCATOR_CORRECTED":
                for row in rows:
                    if str(row["configuration_id"]) == configuration_id and int(row["xi_index"]) == xi_index(xi) and _as_bool(row.get("is_canonical_plot_source", True)):
                        row["is_canonical_plot_source"] = False
                        row["point_status"] = "SUPERSEDED_BY_LOCAL_REPAIR"
            rows.extend(local_rows)
        for audit_row in audit_rows:
            if str(audit_row["configuration_id"]) == configuration_id and xi_index(float(audit_row["xi"])) == xi_index(xi) and bool(audit_row["flagged"]):
                audit_row["repair_id"] = repair_id
                audit_row["repair_status"] = (
                    "REPRODUCED_AFTER_LOCAL_REPAIR" if status == "REPRODUCED"
                    else "LOCATOR_CORRECTED_AFTER_LOCAL_REPAIR" if status == "LOCATOR_CORRECTED"
                    else "UNRESOLVED_AFTER_LOCAL_REPAIR"
                )
        repair_records.append({
            "repair_id": repair_id,
            "configuration_id": configuration_id,
            "xi": xi,
            "affected_positions": positions,
            "status": status,
            "maximum_relative_Omega_difference": relative,
            "tolerance": SPECTRAL_RELATIVE_TOLERANCE,
            "wall_time_seconds": 0.0 if solution is None else solution.wall_time_seconds,
            "peak_rss_bytes": 0 if solution is None else solution.peak_rss_bytes,
            "determinant_evaluations": 0 if solution is None else solution.determinant_evaluations,
            "sigma_evaluations": 0 if solution is None else solution.sigma_evaluations,
            "error": error,
            "smoothing_applied": False,
        })
    rows.sort(key=lambda row: (CONFIGURATIONS.index(str(row["configuration_id"])), int(row["xi_index"]), 0 if str(row["grid_kind"]) == "BASE" else 1, int(row["sorted_position"])))
    return rows, audit_rows, repair_records


def matched_D_rows(
    spectrum_rows: Sequence[Mapping[str, Any]],
    section_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Compare exact rational matched-D points without interpolation."""

    groups: dict[int, dict[str, int]] = {}
    for configuration_id in CONFIGURATIONS:
        for index in range(-40, 41):
            groups.setdefault(TRANSFER_LEVERS[configuration_id] * index, {})[configuration_id] = index
    section_index = {
        (str(row["configuration_id"]), int(row["xi_index"])): row
        for row in section_rows
    }
    result: list[dict[str, Any]] = []
    for q in sorted(key for key, values in groups.items() if len(values) >= 2):
        available = groups[q]
        configurations = [item for item in CONFIGURATIONS if item in available]
        for left_index in range(len(configurations) - 1):
            for right_index in range(left_index + 1, len(configurations)):
                left_configuration = configurations[left_index]
                right_configuration = configurations[right_index]
                left_xi_index = available[left_configuration]
                right_xi_index = available[right_configuration]
                left_roots = _canonical_root_group(spectrum_rows, left_configuration, left_xi_index / 50.0)
                right_roots = _canonical_root_group(spectrum_rows, right_configuration, right_xi_index / 50.0)
                left_section = section_index[(left_configuration, left_xi_index)]
                right_section = section_index[(right_configuration, right_xi_index)]
                beam_relative = max(
                    _relative(float(left_section[field]), float(right_section[field]))
                    for field in ("A_beam", "D_beam", "S_beam", "m", "J", "K", "width")
                )
                for position, (left, right) in enumerate(zip(left_roots, right_roots, strict=True), start=1):
                    left_Omega = float(left["Omega"])
                    right_Omega = float(right["Omega"])
                    relative = _relative(left_Omega, right_Omega)
                    result.append({
                        "D_key_q": q,
                        "D_over_D0_exact": str(Fraction(225 + q, 225)),
                        "D_over_D0": float(Fraction(225 + q, 225)),
                        "common_to_all_three": len(available) == 3,
                        "left_configuration": left_configuration,
                        "right_configuration": right_configuration,
                        "left_xi": left_xi_index / 50.0,
                        "right_xi": right_xi_index / 50.0,
                        "sorted_position": position,
                        "root_role": "PLOTTED" if position <= K_PLOT else "ROOT_9_GUARD",
                        "left_Omega": left_Omega,
                        "right_Omega": right_Omega,
                        "absolute_Omega_difference": abs(left_Omega - right_Omega),
                        "relative_Omega_difference": relative,
                        "left_Lambda": float(left["Lambda"]),
                        "right_Lambda": float(right["Lambda"]),
                        "relative_Lambda_difference": _relative(float(left["Lambda"]), float(right["Lambda"])),
                        "maximum_scaled_sigma_ratio": max(float(left["scaled_sigma_ratio"]), float(right["scaled_sigma_ratio"])),
                        "maximum_boundary_residual": max(float(left["boundary_null_residual"]), float(right["boundary_null_residual"])),
                        "beam_properties_max_relative": beam_relative,
                        "interpolation_used": False,
                        "shared_xi0_provenance": q == 0,
                        "status": "PASS" if beam_relative <= REDUCED_PROPERTY_TOLERANCE and (position == K_GUARD or relative <= SPECTRAL_RELATIVE_TOLERANCE) else "FAIL",
                    })
    return result


def initial_slope_rows(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    root_groups = {
        (configuration, index): _canonical_root_group(rows, configuration, index / 50.0)
        for configuration in CONFIGURATIONS for index in (-1, 0, 1)
    }
    baseline = root_groups[(CENTER_OUTER_TRANSFER, 0)]
    for position in range(1, K_PLOT + 1):
        adjacent_gaps = []
        if position > 1:
            adjacent_gaps.append(float(baseline[position - 1]["Omega"]) - float(baseline[position - 2]["Omega"]))
        if position < K_GUARD:
            adjacent_gaps.append(float(baseline[position]["Omega"]) - float(baseline[position - 1]["Omega"]))
        close = min(adjacent_gaps, default=math.inf) <= 1.0e-8 * max(float(baseline[position - 1]["Omega"]), 1.0)
        slopes: dict[str, float] = {}
        scaled: dict[str, float] = {}
        for configuration in CONFIGURATIONS:
            minus = float(root_groups[(configuration, -1)][position - 1]["Lambda"])
            plus = float(root_groups[(configuration, 1)][position - 1]["Lambda"])
            slopes[configuration] = (plus - minus) / (2.0 * XI_STEP)
            scaled[configuration] = slopes[configuration] / (2.0 * TRANSFER_LEVERS[configuration] / 9.0)
        scale = max(abs(value) for value in scaled.values())
        spread = max(scaled.values()) - min(scaled.values())
        relative_spread = abs(spread) / max(scale, np.finfo(float).tiny)
        if close:
            status = "NOT_APPLICABLE_CLOSE_OR_MULTIPLE"
        elif max(abs(value) for value in slopes.values()) <= SLOPE_INSENSITIVE_ABSOLUTE:
            status = "INSENSITIVE_AT_BASELINE"
        else:
            status = "PASS" if relative_spread <= SLOPE_DIAGNOSTIC_RELATIVE else "DIAGNOSTIC_DEVIATION"
        common_D_secants: dict[str, float] = {}
        for configuration, offset in (
            (CENTER_MIDDLE_TRANSFER, 6),
            (MIDDLE_OUTER_TRANSFER, 3),
            (CENTER_OUTER_TRANSFER, 2),
        ):
            minus_row = _canonical_root_group(rows, configuration, -offset / 50.0)[position - 1]
            plus_row = _canonical_root_group(rows, configuration, offset / 50.0)[position - 1]
            common_D_secants[configuration] = (
                float(plus_row["Lambda"]) - float(minus_row["Lambda"])
            ) / (12.0 / 225.0)
        exact_secant_spread = max(common_D_secants.values()) - min(common_D_secants.values())
        exact_secant_relative = abs(exact_secant_spread) / max(
            max(abs(value) for value in common_D_secants.values()), np.finfo(float).tiny
        )
        result.append({
            "sorted_position": position,
            "baseline_Omega": float(baseline[position - 1]["Omega"]),
            "nearest_adjacent_Omega_gap": min(adjacent_gaps, default=math.inf),
            "dLambda_CM_dxi": slopes[CENTER_MIDDLE_TRANSFER],
            "dLambda_MO_dxi": slopes[MIDDLE_OUTER_TRANSFER],
            "dLambda_CO_dxi": slopes[CENTER_OUTER_TRANSFER],
            "scaled_q_CM": scaled[CENTER_MIDDLE_TRANSFER],
            "scaled_q_MO": scaled[MIDDLE_OUTER_TRANSFER],
            "scaled_q_CO": scaled[CENTER_OUTER_TRANSFER],
            "scaled_q_relative_spread": relative_spread,
            "common_D_secant_CM": common_D_secants[CENTER_MIDDLE_TRANSFER],
            "common_D_secant_MO": common_D_secants[MIDDLE_OUTER_TRANSFER],
            "common_D_secant_CO": common_D_secants[CENTER_OUTER_TRANSFER],
            "common_D_secant_relative_spread": exact_secant_relative,
            "common_D_key_excursion": 6,
            "expected_lever_ratio": "1:2:3",
            "finite_difference_step_xi": XI_STEP,
            "diagnostic_tolerance": SLOPE_DIAGNOSTIC_RELATIVE,
            "hard_spectral_gate": False,
            "status": status,
        })
    return result


def audit_plot_rows(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    canonical = [
        row for row in rows
        if _as_bool(row.get("is_canonical_plot_source", True))
        and int(row["sorted_position"]) <= K_PLOT
    ]
    keys = [
        (str(row["configuration_id"]), int(row["xi_index"]), int(row["sorted_position"]))
        for row in canonical
    ]
    expected = {
        (configuration, index, position)
        for configuration in CONFIGURATIONS
        for index in range(-40, 41)
        for position in range(1, K_PLOT + 1)
    }
    return {
        "status": "PASS" if len(keys) == len(set(keys)) == len(expected) and set(keys) == expected else "FAIL",
        "row_count": len(canonical),
        "expected_row_count": 1944,
        "duplicate_key_count": len(keys) - len(set(keys)),
        "missing_key_count": len(expected - set(keys)),
        "root9_plotted": any(int(row["sorted_position"]) == K_GUARD for row in canonical),
    }


def _save_figure_atomic(figure: Any, path: Path) -> None:
    target = Path(path)
    temporary = target.with_name(target.stem + ".tmp" + target.suffix)
    figure.savefig(temporary, dpi=180, bbox_inches="tight")
    os.replace(temporary, target)


def create_plots_from_csv(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    """Render both figures from CSV only, without importing root physics."""

    started = time.perf_counter()
    target = Path(output_dir)
    rows = _read_csv(target / SPECTRUM_FILENAME)
    spectrum_audit = audit_spectrum_rows(rows)
    plot_audit = audit_plot_rows(rows)
    if spectrum_audit["status"] != "PASS" or plot_audit["status"] != "PASS":
        raise RuntimeError(f"plot-only rejected incomplete data: {spectrum_audit}; {plot_audit}")
    plot_rows = [
        row for row in rows
        if _as_bool(row.get("is_canonical_plot_source", True))
        and int(row["sorted_position"]) <= K_PLOT
    ]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    colors = plt.get_cmap("tab10").colors[:K_PLOT]
    titles = (
        r"(a) CENTER $\leftrightarrow$ MIDDLE",
        r"(b) MIDDLE $\leftrightarrow$ OUTER",
        r"(c) CENTER $\leftrightarrow$ OUTER",
    )
    figure, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), sharex=True, sharey=True)
    for axis, configuration_id, title in zip(axes, CONFIGURATIONS, titles, strict=True):
        for position in range(1, K_PLOT + 1):
            selected = sorted(
                (row for row in plot_rows if str(row["configuration_id"]) == configuration_id and int(row["sorted_position"]) == position),
                key=lambda row: int(row["xi_index"]),
            )
            axis.plot(
                [float(row["xi"]) for row in selected],
                [float(row["Lambda"]) for row in selected],
                color=colors[position - 1], linewidth=1.25, linestyle="-",
            )
        axis.axvline(0.0, color="0.45", linewidth=0.7, alpha=0.7)
        axis.set_title(title)
        axis.set_xlabel(r"$\xi$")
        axis.grid(True, alpha=0.22, linewidth=0.5)
    axes[0].set_ylabel(r"$\Lambda$")
    position_handles = [Line2D([0], [0], color=colors[index - 1], linewidth=1.5, label=f"k={index}") for index in range(1, K_PLOT + 1)]
    figure.legend(handles=position_handles, loc="upper center", ncol=8, frameon=False)
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.91))
    _save_figure_atomic(figure, target / XI_PLOT_FILENAME)
    plt.close(figure)

    by_q: dict[int, dict[str, list[Mapping[str, Any]]]] = {}
    for row in plot_rows:
        by_q.setdefault(int(row["D_key_q"]), {}).setdefault(str(row["configuration_id"]), []).append(row)
    figure, axis = plt.subplots(figsize=(8.8, 5.8))
    priority = (CENTER_OUTER_TRANSFER, MIDDLE_OUTER_TRANSFER, CENTER_MIDDLE_TRANSFER)
    for position in range(1, K_PLOT + 1):
        master_x: list[float] = []
        master_y: list[float] = []
        for q in sorted(by_q):
            candidates = by_q[q]
            chosen = next(configuration for configuration in priority if configuration in candidates)
            row = next(item for item in candidates[chosen] if int(item["sorted_position"]) == position)
            master_x.append(float(row["D_over_D0"]))
            master_y.append(float(row["Lambda"]))
        axis.plot(master_x, master_y, color=colors[position - 1], linewidth=1.25)
    marker_map = {
        CENTER_MIDDLE_TRANSFER: "o",
        MIDDLE_OUTER_TRANSFER: "s",
        CENTER_OUTER_TRANSFER: "^",
    }
    for configuration_id in CONFIGURATIONS:
        for position in range(1, K_PLOT + 1):
            selected = sorted(
                (row for row in plot_rows if str(row["configuration_id"]) == configuration_id and int(row["sorted_position"]) == position),
                key=lambda row: int(row["D_key_q"]),
            )[::4]
            axis.plot(
                [float(row["D_over_D0"]) for row in selected],
                [float(row["Lambda"]) for row in selected],
                linestyle="none", marker=marker_map[configuration_id], markersize=2.8,
                markerfacecolor="none", markeredgewidth=0.65, color=colors[position - 1],
            )
    axis.set_xlabel(r"$D/D_0$")
    axis.set_ylabel(r"$\Lambda$")
    axis.grid(True, alpha=0.22, linewidth=0.5)
    position_legend = axis.legend(handles=position_handles, loc="upper left", ncol=2, frameon=False, title="sorted position")
    axis.add_artist(position_legend)
    marker_handles = [
        Line2D([0], [0], marker=marker_map[configuration], linestyle="none", markerfacecolor="none", color="0.2", label=configuration.replace("_TRANSFER", ""))
        for configuration in CONFIGURATIONS
    ]
    axis.legend(handles=marker_handles, loc="lower right", frameon=False, title="transfer")
    figure.tight_layout()
    _save_figure_atomic(figure, target / MASTER_PLOT_FILENAME)
    plt.close(figure)
    return {
        "wall_time_seconds": time.perf_counter() - started,
        "xi_plot": (target / XI_PLOT_FILENAME).as_posix(),
        "master_plot": (target / MASTER_PLOT_FILENAME).as_posix(),
        "xi_panel_count": 3,
        "spectral_lines_per_xi_panel": K_PLOT,
        "master_solid_lines": K_PLOT,
        "root9_plotted": False,
        "interpolation_used": False,
        "smoothing_applied": False,
        "root_calculation_count": 0,
        "spectrum_audit": spectrum_audit,
        "plot_audit": plot_audit,
    }


def four_vs_six_ply_rows(
    rows: Sequence[Mapping[str, Any]], section_rows: Sequence[Mapping[str, Any]]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Use frozen RLB-2E CSV only at exact reduced-property matches."""

    directory = PREDECESSOR_RESULT_DIRS["RLB-2E"]
    manifest_path = directory / "run_manifest.json"
    spectrum_path = directory / "spectrum_roots.csv"
    section_path = directory / "section_properties.csv"
    if not (manifest_path.is_file() and spectrum_path.is_file() and section_path.is_file()):
        return [], {"status": "NOT_AVAILABLE", "blocking": False}
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    recorded = manifest.get("output_hashes", {})
    hashes_ok = bool(
        manifest.get("scientific_status") == "PASS"
        and recorded.get(spectrum_path.name) == _sha256(spectrum_path)
        and recorded.get(section_path.name) == _sha256(section_path)
    )
    if not hashes_ok:
        return [], {"status": "FROZEN_HASH_MISMATCH", "blocking": False}
    old_spectrum = _read_csv(spectrum_path)
    old_sections = _read_csv(section_path)
    old_root_index: dict[tuple[str, int], list[Mapping[str, Any]]] = {}
    for configuration in ("BOTH_OUTER_STIFF", "BOTH_INNER_STIFF"):
        for m in range(0, 41):
            group = [
                row for row in old_spectrum
                if str(row["configuration_id"]) == configuration
                and round(float(row["chi"]), 10) == round(m / 50.0, 10)
                and str(row["grid_kind"]) == "BASE"
                and _as_bool(row.get("is_canonical_plot_source", True))
            ]
            group.sort(key=lambda row: int(row["sorted_position"]))
            if len(group) == K_GUARD:
                old_root_index[(configuration, m)] = group
    old_section_index: dict[tuple[str, int], Mapping[str, Any]] = {}
    for row in old_sections:
        if int(row["arm_id"]) != 1:
            continue
        old_section_index[(str(row["configuration_id"]), int(round(float(row["chi"]) * 50.0)))] = row
    six_section_index = {(str(row["configuration_id"]), int(row["xi_index"])): row for row in section_rows}
    result: list[dict[str, Any]] = []
    for configuration in CONFIGURATIONS:
        for index in range(-40, 41):
            q = TRANSFER_LEVERS[configuration] * index
            if q == 0:
                old_configuration, m = "BOTH_OUTER_STIFF", 0
            else:
                numerator = 8 * abs(q)
                if numerator % 27:
                    continue
                m = numerator // 27
                if m > 40:
                    continue
                old_configuration = "BOTH_OUTER_STIFF" if q > 0 else "BOTH_INNER_STIFF"
            old_roots = old_root_index.get((old_configuration, m))
            old_section = old_section_index.get((old_configuration, m))
            if old_roots is None or old_section is None:
                continue
            six_roots = _canonical_root_group(rows, configuration, index / 50.0)
            six_section = six_section_index[(configuration, index)]
            property_relative = max(
                _relative(float(six_section[six_name]), float(old_section[old_name]))
                for six_name, old_name in (
                    ("A_beam", "A_beam"), ("D_beam", "D_beam"),
                    ("S_beam", "S_beam"), ("m", "m"), ("J", "J"),
                )
            )
            for position, (six, old) in enumerate(zip(six_roots, old_roots, strict=True), start=1):
                relative = _relative(float(six["Omega"]), float(old["Omega"]))
                result.append({
                    "six_configuration": configuration,
                    "six_xi": index / 50.0,
                    "RLB2E_configuration": old_configuration,
                    "RLB2E_chi": m / 50.0,
                    "D_key_q": q,
                    "D_over_D0": float(Fraction(225 + q, 225)),
                    "sorted_position": position,
                    "root_role": "PLOTTED" if position <= K_PLOT else "ROOT_9_GUARD",
                    "six_Omega": float(six["Omega"]),
                    "four_Omega": float(old["Omega"]),
                    "relative_Omega_difference": relative,
                    "beam_properties_max_relative": property_relative,
                    "interpolation_used": False,
                    "status": "PASS" if property_relative <= REDUCED_PROPERTY_TOLERANCE and (position == K_GUARD or relative <= SPECTRAL_RELATIVE_TOLERANCE) else "FAIL",
                })
    return result, {
        "status": "PASS" if result and all(row["status"] == "PASS" for row in result) else "NOT_AVAILABLE" if not result else "FAIL",
        "blocking": False,
        "row_count": len(result),
        "maximum_plotted_relative_Omega_difference": max((float(row["relative_Omega_difference"]) for row in result if int(row["sorted_position"]) <= K_PLOT), default=None),
        "RLB2E_spectrum_sha256": _sha256(spectrum_path),
        "RLB2E_section_sha256": _sha256(section_path),
        "RLB2E_result_tree_sha256": _sha_tree(directory),
    }


def _minimum_adjacent_gaps(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for configuration in CONFIGURATIONS:
        best = (math.inf, math.nan, -1)
        for xi in xi_grid():
            group = _canonical_root_group(rows, configuration, float(xi))
            values = np.asarray([float(row["Lambda"]) for row in group[:K_PLOT]])
            gaps = np.diff(values)
            index = int(np.argmin(gaps))
            if float(gaps[index]) < best[0]:
                best = (float(gaps[index]), float(xi), index + 1)
        result.append({
            "configuration_id": configuration,
            "minimum_adjacent_Lambda_gap": best[0],
            "xi": best[1],
            "between_sorted_positions": [best[2], best[2] + 1],
            "interpretation": "candidate interval only; no modal-identity claim",
        })
    return result


def _endpoint_and_monotonicity(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for configuration in CONFIGURATIONS:
        for position in range(1, K_PLOT + 1):
            values = [
                float(_canonical_root_group(rows, configuration, index / 50.0)[position - 1]["Lambda"])
                for index in range(-40, 41)
            ]
            differences = np.diff(np.asarray(values))
            allowance = 1.0e-10 * max(max(abs(value) for value in values), 1.0)
            classification = "NONDECREASING" if float(np.min(differences)) >= -allowance else "NOT_MONOTONE_ON_SORTED_GRID"
            result.append({
                "configuration_id": configuration,
                "sorted_position": position,
                "Lambda_xi_minus_0p8": values[0],
                "Lambda_xi_zero": values[40],
                "Lambda_xi_plus_0p8": values[-1],
                "relative_endpoint_change": (values[-1] - values[0]) / max(abs(values[0]), np.finfo(float).tiny),
                "minimum_forward_Lambda_increment": float(np.min(differences)),
                "numerical_allowance": allowance,
                "classification": classification,
                "spectrum_semantics": "sorted_positions",
            })
    return result


def _root_quality_summary(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    base_rows = [row for row in rows if str(row["grid_kind"]) == "BASE"]
    guards = [row for row in base_rows if int(row["sorted_position"]) == K_GUARD]
    return {
        "maximum_base_scaled_sigma_ratio": max(float(row["scaled_sigma_ratio"]) for row in base_rows),
        "maximum_base_boundary_null_residual": max(float(row["boundary_null_residual"]) for row in base_rows),
        "maximum_root9_scaled_sigma_ratio": max(float(row["scaled_sigma_ratio"]) for row in guards),
        "maximum_root9_boundary_null_residual": max(float(row["boundary_null_residual"]) for row in guards),
        "minimum_root9_right_margin_Omega": min(float(row["root9_right_margin_Omega"]) for row in guards),
        "maximum_unresolved_candidates_below_root9": max(int(row["unresolved_candidates_below_root9"]) for row in base_rows),
        "root9_guard_count": len(guards),
        "roots_above_9_computed": any(_as_bool(row["roots_above_9_computed"]) for row in rows),
    }


def _output_hashes(output_dir: Path) -> dict[str, str]:
    return {
        path.name: _sha256(path)
        for path in sorted(Path(output_dir).iterdir())
        if path.is_file() and path.name != MANIFEST_FILENAME
    }


def _runtime_metadata(
    point_records: Sequence[Mapping[str, Any]],
    repair_records: Sequence[Mapping[str, Any]],
    *,
    plot_only_seconds: float,
    finalization_invocation_seconds: float,
    peak_rss_bytes: int,
) -> dict[str, Any]:
    """Build honest timing/evaluation provenance from retained transactions.

    A completed resume can reuse local-repair rows produced by an earlier
    invocation.  If that invocation ended before its manifest was written,
    the scientific repair result is retained but its wall time and evaluation
    counters are not.  Those missing counters must be reported as unavailable,
    not silently replaced by zero or by the duration of the final resume.
    """

    base_seconds = sum(float(record.get("wall_time_seconds", 0.0)) for record in point_records)
    base_evaluations = sum(int(record.get("determinant_evaluations", 0)) for record in point_records)
    base_solve_count = sum(int(record.get("determinant_evaluations", 0)) > 0 for record in point_records)
    repair_metrics_complete = all(
        not _as_bool(record.get("reused_existing_local_repair", False))
        for record in repair_records
    )
    retained_repair_seconds = sum(float(record.get("wall_time_seconds", 0.0)) for record in repair_records)
    retained_repair_evaluations = sum(
        int(record.get("determinant_evaluations", 0)) for record in repair_records
    )
    recorded_lower_bound = base_seconds + retained_repair_seconds + float(plot_only_seconds)
    provenance = (
        "COMPLETE_RETAINED_TRANSACTION_METRICS"
        if repair_metrics_complete
        else "QUALIFIED_LOCAL_REPAIR_METRICS_NOT_RETAINED_AFTER_INTERRUPTED_PRE_MANIFEST_FINALIZATION"
    )
    return {
        "measurement_scope_version": 2,
        "recorded_base_point_wall_time_sum_seconds": base_seconds,
        "base_point_wall_time_sum_seconds": base_seconds,
        "local_repair_wall_time_sum_seconds": retained_repair_seconds if repair_metrics_complete else None,
        "local_repair_wall_time_provenance": provenance,
        "plot_generation_seconds": float(plot_only_seconds),
        "plot_only_seconds": float(plot_only_seconds),
        "latest_finalization_invocation_wall_time_seconds": float(finalization_invocation_seconds),
        "finalization_invocation_seconds": float(finalization_invocation_seconds),
        "recorded_workflow_wall_time_lower_bound_seconds": recorded_lower_bound,
        "campaign_total_wall_time_seconds": recorded_lower_bound if repair_metrics_complete else None,
        "campaign_total_wall_time_available": repair_metrics_complete,
        "total_measured_workflow_seconds": recorded_lower_bound if repair_metrics_complete else None,
        "total_measured_workflow_status": provenance,
        "recorded_peak_rss_bytes_lower_bound": int(peak_rss_bytes),
        "campaign_peak_rss_complete": repair_metrics_complete,
        "peak_rss_bytes": int(peak_rss_bytes),
        "recorded_base_determinant_evaluations": base_evaluations,
        "base_determinant_evaluations": base_evaluations,
        "recorded_local_repair_determinant_evaluations": (
            retained_repair_evaluations if repair_metrics_complete else None
        ),
        "local_repair_determinant_evaluations": (
            retained_repair_evaluations if repair_metrics_complete else None
        ),
        "determinant_evaluations_recorded_lower_bound": base_evaluations + retained_repair_evaluations,
        "determinant_evaluations": (
            base_evaluations + retained_repair_evaluations if repair_metrics_complete else None
        ),
        "sigma_evaluations_recorded_lower_bound": base_evaluations + retained_repair_evaluations,
        "sigma_evaluations": (
            base_evaluations + retained_repair_evaluations if repair_metrics_complete else None
        ),
        "base_root_solve_count": base_solve_count,
        "accepted_local_refinement_group_count": len(repair_records),
        "local_repair_solve_count": len(repair_records),
        "local_refinement_groups_executed_in_latest_finalization": sum(
            int(record.get("determinant_evaluations", 0)) > 0 for record in repair_records
        ),
        "local_repair_solve_count_with_retained_metrics": sum(
            int(record.get("determinant_evaluations", 0)) > 0 for record in repair_records
        ),
        "total_root_solve_count": base_solve_count + len(repair_records),
        "parallel_spectral_workers": 0,
    }


def _report_text(manifest: Mapping[str, Any]) -> str:
    status = manifest["scientific_status"]
    gate = manifest["constitutive_gate"]
    counts = manifest["counts"]
    quality = manifest["root_quality_summary"]
    matched = manifest["matched_D_summary"]
    neighbour = manifest["neighbour_audit"]
    runtime = manifest["runtime"]
    benchmark = manifest["benchmark"]
    slope = manifest["initial_slope_summary"]
    endpoints = manifest["endpoint_and_monotonicity"]
    optional = manifest["four_vs_six_ply_equivalence"]
    minimum_gaps = manifest["minimum_adjacent_sorted_gaps"]
    endpoint_lines = "\n".join(
        f"- `{item['configuration_id']}`, position {item['sorted_position']}: "
        f"$\\Delta\\Lambda/\\Lambda_-={float(item['relative_endpoint_change']):.6g}$, "
        f"`{item['classification']}`."
        for item in endpoints
    )
    gap_lines = "\n".join(
        f"- `{item['configuration_id']}`: {float(item['minimum_adjacent_Lambda_gap']):.9g} "
        f"при $\\xi={float(item['xi']):.2f}$ между positions "
        f"{item['between_sorted_positions'][0]}–{item['between_sorted_positions'][1]}."
        for item in minimum_gaps
    )
    if runtime["total_measured_workflow_seconds"] is None:
        runtime_lines = (
            f"Сумма сохранённых wall times BASE root transactions равна "
            f"{runtime['base_point_wall_time_sum_seconds']:.1f} s; построение рисунков — "
            f"{runtime['plot_only_seconds']:.1f} s. Нижняя граница измеренного времени "
            f"составляет {runtime['recorded_workflow_wall_time_lower_bound_seconds']:.1f} s. "
            "Времена и evaluation counters 25 ранее выполненных local repairs не были "
            "сохранены после прерванной pre-manifest finalization; repairs повторно не "
            "рассчитывались и отсутствующие метрики не заменялись нулями."
        )
        evaluation_lines = (
            f"Для BASE points сохранено {runtime['base_determinant_evaluations']} "
            "determinant/SVD evaluations; полное число с учётом local repairs "
            "квалифицировано как недоступное."
        )
    else:
        runtime_lines = (
            "Сумма сохранённых wall times root transactions и plot rendering равна "
            f"{runtime['total_measured_workflow_seconds']:.1f} s."
        )
        evaluation_lines = (
            f"Выполнено {runtime['determinant_evaluations']} determinant/SVD evaluations."
        )
    return rf"""# RLB-2J: попарный перенос жёсткости в шестислойном ламинате

## Цель и расчётный контракт

Рассмотрены три способа попарного перераспределения внутриплоскостной
жёсткости между зеркальными парами слоёв симметричного ламината
O/M/C/C/M/O. Все шесть слоёв имеют одинаковую толщину $h/6$ и ориентацию
$0^\circ$. Геометрия конструкции фиксирована: $\mu=\tau=0$,
$\beta=30^\circ$, $L_1=L_2=1$, $b=0.20$, $h=0.05$, $K=5/6$.

Базовый материал задан величинами $E_1=1.1$, $E_2=0.9$,
$\nu_{{12}}=0.3$, $G_{{12}}=G_{{13}}=G_{{23}}=1/2.6$ и $\rho=1$.
Множитель пары изменяет только $E_1$, $E_2$ и $G_{{12}}$; $G_{{13}}$,
$G_{{23}}$, плотность и коэффициент Пуассона остаются постоянными.

## Аналитические соотношения

В порядке C, M, O веса мембранной и изгибной жёсткости равны

$$w_A=(1,1,1)/3,\qquad w_D=(1,7,19)/27.$$

Для направлений $v_{{CM}}=(-1,1,0)$, $v_{{MO}}=(0,-1,1)$ и
$v_{{CO}}=(-1,0,1)$ получены

$$\frac{{D_{{CM}}}}{{D_0}}=1+\frac{{2\xi}}9,\qquad
  \frac{{D_{{MO}}}}{{D_0}}=1+\frac{{4\xi}}9,\qquad
  \frac{{D_{{CO}}}}{{D_0}}=1+\frac{{2\xi}}3.$$

Проекции на целочисленные веса $(1,7,19)$ равны 6, 12 и 18, поэтому
изгибные рычаги относятся как 1:2:3. Направление RLB-2I
$v_N=(2,-3,1)$, напротив, ортогонально как мембранным, так и изгибным
весам и сохраняет одновременно $A$ и $D$.

Constitutive gate имеет статус **{gate['status']}**. Максимальные residuals:
$A$ — {gate['maximum_residuals']['A_matrix_invariance_relative']:.3e},
$D$ — {gate['maximum_residuals']['D_matrix_formula_relative']:.3e},
поперечный сдвиг — {gate['maximum_residuals']['shear_matrix_invariance_relative']:.3e};
максимум для reduced $A,S,m,J$ равен
{max(gate['maximum_residuals'][name] for name in ('A_beam_invariance_relative','S_beam_invariance_relative','mass_invariance_relative','rotary_inertia_invariance_relative')):.3e}.

## Frequency-map workflow

Использован локальный instance `frequency-map-v1`: `fast_plot`,
`sorted_positions`, сетка $\xi=-0.80,-0.78,\ldots,0.80$. Продолжение
выполнялось отдельно в отрицательную и положительную стороны от общего
$\xi=0$. Прогнозы задавали только интервалы локализации; все частоты
получены из characteristic matrix. Позиции 1–8 показаны на рисунках, root 9
использован только как completeness guard. Roots 10 и выше не вычислялись.

Получено {counts['base_groups']}/243 logical groups и
{counts['base_rows']}/2187 BASE rows. Общий $\xi=0$ рассчитан один раз;
число уникальных физических BASE solves равно {counts['unique_physical_base_solves']}.
Neighbour audit отметил {neighbour['flagged_point_count']} point(s), локально
проверено {neighbour['repair_count']}, unresolved point(s):
{neighbour['unresolved_point_count']}.

Benchmark дал ETA {benchmark['conservative_eta_seconds']:.1f} s при лимите
2700 s. {runtime_lines} Доступная нижняя оценка максимума RSS —
{runtime['peak_rss_bytes']/2**20:.1f} MiB. {evaluation_lines}

Максимальные root-quality diagnostics составляют
{quality['maximum_base_scaled_sigma_ratio']:.3e} для
$\sigma_{{min}}/\sigma_{{max}}$ и
{quality['maximum_base_boundary_null_residual']:.3e} для boundary residual.
Все {quality['root9_guard_count']} guards имеют минимальный правый запас
{quality['minimum_root9_right_margin_Omega']:.3g} по $\Omega$.

## Matched-D и начальные наклоны

Точное объединение выполнено по целому ключу $q$, где
$D/D_0=1+q/225$; интерполяция частот не применялась. Найдено
{matched['unique_matched_D_values']} значений $D/D_0$, общих хотя бы двум
семействам, в том числе {matched['common_all_three_D_values']} значений,
общих всем трём. Максимальная относительная разность positions 1–8 равна
{matched['maximum_plotted_relative_Omega_difference']:.3e}; статус collapse:
**{matched['status']}**.

Диагностика начальных наклонов использует центральные разности на основной
сетке. Масштабированные чувствительности имеют максимум относительного
разброса {slope['maximum_scaled_q_relative_spread']:.3e}. Поскольку три
разности соответствуют разным шагам по $D/D_0$, это не жёсткий spectral
gate. Дополнительные секущие на одном и том же точном интервале по $D/D_0$
имеют максимальный разброс {slope['maximum_common_D_secant_relative_spread']:.3e}
и непосредственно подтверждают отношение рычагов 1:2:3 для применимых
простых sorted positions.

Дополнительное точное сравнение с frozen RLB-2E содержит
{optional['row_count']} rows и имеет статус **{optional['status']}**;
максимальная относительная разность plotted $\Omega$ равна
{optional['maximum_plotted_relative_Omega_difference']:.3e}. RLB-2E при
этом не пересчитывался.

Минимальные соседние gaps в $\Lambda$:

{gap_lines}

## Изменения частот на концах сетки

{endpoint_lines}

Эти классификации относятся только к независимо отсортированным значениям и
не устанавливают modal identity, crossing, veering или localization.

## Статус и ограничения

**RLB-2J: {status}.** Для фиксированной конструкции и объявленной конечной
сетки подтверждено: при постоянных $A,S,m,J$ частотный эффект попарного
переноса определяется изменением $D$ с рычагами 1:2:3. После перехода к
фактическому $D/D_0$ три семейства совпадают в пределах заявленного допуска,
что проверяет структуру текущей одномерной редуцированной модели.

Branch tracking, MAC, формы, анализ энергии или напряжений, Ritz, FEM,
сглаживание и certified audit не выполнялись. Производственные physics-модули
и predecessor result trees не изменялись.
"""


def _refresh_completed_metadata(
    target: Path, manifest: Mapping[str, Any]
) -> dict[str, Any]:
    """Refresh completed-run provenance without touching numerical outputs."""

    refreshed = json.loads(json.dumps(_json_value(manifest)))
    repair_records = [
        dict(record)
        for record in refreshed.get("neighbour_audit", {}).get("repair_records", [])
    ]
    previous_runtime = dict(refreshed.get("runtime", {}))
    finalization_seconds = float(
        previous_runtime.get(
            "finalization_invocation_seconds",
            previous_runtime.get("total_measured_workflow_seconds") or 0.0,
        )
    )
    refreshed["runtime"] = _runtime_metadata(
        [dict(record) for record in refreshed.get("point_records", [])],
        repair_records,
        plot_only_seconds=float(previous_runtime.get("plot_only_seconds", 0.0)),
        finalization_invocation_seconds=finalization_seconds,
        peak_rss_bytes=int(previous_runtime.get("peak_rss_bytes", 0)),
    )
    refreshed["analysis_script_sha256"] = _sha256(Path(__file__))
    refreshed["git"] = _git_state()
    refreshed["contract"] = _json_value(contract_payload())
    refreshed["contract_sha256"] = contract_hash()
    refreshed["production_physics_hashes"] = {
        path: _sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS
    }
    refreshed["production_physics_hashes_match_frozen"] = (
        refreshed["production_physics_hashes"] == EXPECTED_PRODUCTION_PHYSICS_HASHES
    )
    refreshed["predecessor_result_tree_hashes"] = {
        stage: _sha_tree(path) for stage, path in PREDECESSOR_RESULT_DIRS.items()
    }
    refreshed["predecessor_result_trees_preserved"] = all(
        refreshed["predecessor_result_tree_hashes"].get(stage) == expected
        for stage, expected in EXPECTED_PREDECESSOR_TREE_HASHES.items()
    )
    report_text = _report_text(refreshed)
    current_report = (target / REPORT_FILENAME).read_text(encoding="utf-8")
    comparison = dict(refreshed)
    comparison["output_hashes"] = manifest.get("output_hashes", {})
    if comparison == dict(manifest) and report_text == current_report:
        return dict(manifest)
    refreshed["metadata_refreshed_at_utc"] = _utc_now()
    refreshed["metadata_refresh_root_calculation_count"] = 0
    refreshed["metadata_refresh_numerical_outputs_changed"] = False
    _atomic_write_text(target / REPORT_FILENAME, report_text)
    refreshed["output_hashes"] = _output_hashes(target)
    _atomic_write_json(target / MANIFEST_FILENAME, refreshed)
    return refreshed


def finalize_outputs(
    target: Path,
    rows: list[dict[str, Any]],
    section_rows: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    benchmark: Mapping[str, Any],
    point_records: Sequence[Mapping[str, Any]],
    started_at: str,
    workflow_started: float,
) -> dict[str, Any]:
    audit_rows = neighbour_audit_rows(rows)
    rows, audit_rows, repair_records = apply_local_repairs(rows, audit_rows)
    _atomic_write_csv(target / SPECTRUM_FILENAME, rows, SPECTRUM_FIELDS)
    _atomic_write_csv(target / AUDIT_FILENAME, audit_rows)
    spectrum_audit = audit_spectrum_rows(rows)
    if spectrum_audit["status"] != "PASS":
        raise RuntimeError(f"Final BASE spectrum audit failed: {spectrum_audit}")
    matched = matched_D_rows(rows, section_rows)
    _atomic_write_csv(target / MATCHED_D_FILENAME, matched)
    slope = initial_slope_rows(rows)
    _atomic_write_csv(target / SLOPE_FILENAME, slope)
    optional_rows, optional_summary = four_vs_six_ply_rows(rows, section_rows)
    if optional_rows:
        _atomic_write_csv(target / FOUR_VS_SIX_FILENAME, optional_rows)
    plots = create_plots_from_csv(target)
    unresolved = sum(record["status"] == "UNRESOLVED" for record in repair_records)
    matched_plotted = [row for row in matched if int(row["sorted_position"]) <= K_PLOT]
    matched_summary = {
        "status": "PASS" if matched and all(row["status"] == "PASS" for row in matched) else "FAIL",
        "row_count": len(matched),
        "expected_pairwise_row_count": 855,
        "unique_matched_D_values": len({int(row["D_key_q"]) for row in matched}),
        "common_all_three_D_values": len({int(row["D_key_q"]) for row in matched if _as_bool(row["common_to_all_three"])}),
        "maximum_plotted_relative_Omega_difference": max(float(row["relative_Omega_difference"]) for row in matched_plotted),
        "maximum_guard_relative_Omega_difference": max(float(row["relative_Omega_difference"]) for row in matched if int(row["sorted_position"]) == K_GUARD),
        "maximum_beam_properties_relative": max(float(row["beam_properties_max_relative"]) for row in matched),
        "interpolation_used": False,
    }
    slope_summary = {
        "row_count": len(slope),
        "PASS_count": sum(row["status"] == "PASS" for row in slope),
        "insensitive_count": sum(row["status"] == "INSENSITIVE_AT_BASELINE" for row in slope),
        "not_applicable_close_or_multiple_count": sum(row["status"] == "NOT_APPLICABLE_CLOSE_OR_MULTIPLE" for row in slope),
        "diagnostic_deviation_count": sum(row["status"] == "DIAGNOSTIC_DEVIATION" for row in slope),
        "maximum_scaled_q_relative_spread": max(float(row["scaled_q_relative_spread"]) for row in slope),
        "maximum_common_D_secant_relative_spread": max(float(row["common_D_secant_relative_spread"]) for row in slope),
        "hard_gate": False,
    }
    endpoint = _endpoint_and_monotonicity(rows)
    base_solve_count = sum(int(record.get("determinant_evaluations", 0)) > 0 for record in point_records)
    peak = max(
        [_peak_rss_bytes(), *[int(record.get("peak_rss_bytes", 0)) for record in point_records], *[int(record.get("peak_rss_bytes", 0)) for record in repair_records]]
    )
    status = (
        "PASS"
        if constitutive["status"] == "PASS"
        and spectrum_audit["status"] == "PASS"
        and matched_summary["status"] == "PASS"
        and unresolved == 0
        and plots["plot_audit"]["status"] == "PASS"
        else "PARTIAL_PASS" if constitutive["status"] == "PASS" and spectrum_audit["base_group_count"] > 0 else "FAIL"
    )
    predecessor_hashes = {stage: _sha_tree(path) for stage, path in PREDECESSOR_RESULT_DIRS.items()}
    production_hashes = {path: _sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS}
    runtime = _runtime_metadata(
        point_records,
        repair_records,
        plot_only_seconds=float(plots["wall_time_seconds"]),
        finalization_invocation_seconds=time.perf_counter() - workflow_started,
        peak_rss_bytes=peak,
    )
    counts = {
        "base_groups": spectrum_audit["base_group_count"],
        "base_rows": spectrum_audit["base_row_count"],
        "plotted_base_rows": 243 * K_PLOT,
        "root9_guards": 243,
        "local_refinement_rows": sum(str(row["grid_kind"]) == "LOCAL_REFINEMENT" for row in rows),
        "unique_physical_base_solves": base_solve_count,
        "shared_xi0_reused_logical_groups": 2,
        "logical_base_groups": 243,
        "unique_physical_base_solve_contract": 241,
    }
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "scientific_status": status,
        "completed_at_utc": _utc_now(),
        "git": _git_state(),
        "contract": contract_payload(),
        "contract_sha256": contract_hash(),
        "analysis_script_sha256": _sha256(Path(__file__)),
        "production_physics_hashes": production_hashes,
        "production_physics_hashes_match_frozen": production_hashes == EXPECTED_PRODUCTION_PHYSICS_HASHES,
        "predecessor_result_tree_hashes": predecessor_hashes,
        "predecessor_result_trees_preserved": all(predecessor_hashes.get(key) == value for key, value in EXPECTED_PREDECESSOR_TREE_HASHES.items()),
        "constitutive_gate": dict(constitutive),
        "counts": counts,
        "spectrum_audit": spectrum_audit,
        "root_quality_summary": _root_quality_summary(rows),
        "benchmark": dict(benchmark),
        "point_records": list(point_records),
        "neighbour_audit": {
            "criterion": "centered secant residual > median+8*MAD and >1e-3, plus gap/root-quality triggers",
            "row_count": len(audit_rows),
            "flagged_point_count": len(flagged_repair_points(audit_rows)),
            "repair_count": len(repair_records),
            "repair_records": repair_records,
            "unresolved_point_count": unresolved,
            "smoothing_applied": False,
        },
        "matched_D_summary": matched_summary,
        "initial_slope_summary": slope_summary,
        "endpoint_and_monotonicity": endpoint,
        "four_vs_six_ply_equivalence": optional_summary,
        "minimum_adjacent_sorted_gaps": _minimum_adjacent_gaps(rows),
        "plots": plots,
        "runtime": runtime,
        "exclusions_confirmed": {
            "roots_above_9_computed": False,
            "spectral_sweep_runner_used": False,
            "certified_audit_run": False,
            "full_inventory_run": False,
            "parallel_spectral_workers": 0,
            "interpolation_based_frequencies": False,
            "smoothing": False,
            "branch_tracking": False,
            "MAC": False,
            "mode_shapes": False,
            "energy_or_stress_reconstruction": False,
            "Ritz": False,
            "FEM": False,
            "commit": False,
            "push": False,
        },
    }
    checkpoint = _checkpoint_payload(rows, point_records, constitutive=constitutive, started_at=started_at, benchmark_status="PASS")
    checkpoint.update({
        "scientific_status": status,
        "completed_at_utc": _utc_now(),
        "local_repair_count": len(repair_records),
        "terminal_unresolved_points": [
            {"configuration_id": record["configuration_id"], "xi": record["xi"]}
            for record in repair_records if record["status"] == "UNRESOLVED"
        ],
    })
    _atomic_write_json(target / CHECKPOINT_FILENAME, checkpoint)
    _atomic_write_text(target / REPORT_FILENAME, _report_text(manifest))
    manifest["output_hashes"] = _output_hashes(target)
    _atomic_write_json(target / MANIFEST_FILENAME, manifest)
    return manifest


def prepare_constitutive_outputs(output_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    gate, rows = compute_constitutive_data()
    _atomic_write_csv(target / SECTION_FILENAME, rows, SECTION_FIELDS)
    if gate["status"] != "PASS":
        raise RuntimeError(f"Constitutive gate failed before spectra: {gate}")
    return gate, rows


def _existing_rows(output_dir: Path) -> list[dict[str, Any]]:
    path = Path(output_dir) / SPECTRUM_FILENAME
    return [] if not path.is_file() else [dict(row) for row in _read_csv(path)]


def run_workflow(
    output_dir: Path = DEFAULT_OUTPUT_DIR, *, missing_only: bool = False
) -> dict[str, Any]:
    del missing_only
    global ROOT_CALCULATION_COUNT
    invocation_start_count = ROOT_CALCULATION_COUNT
    workflow_started = time.perf_counter()
    started_at = _utc_now()
    target = Path(output_dir)
    rows = _existing_rows(target)
    required = (
        SPECTRUM_FILENAME, SECTION_FILENAME, AUDIT_FILENAME, MATCHED_D_FILENAME,
        SLOPE_FILENAME, BENCHMARK_FILENAME, CHECKPOINT_FILENAME,
        XI_PLOT_FILENAME, MASTER_PLOT_FILENAME, REPORT_FILENAME, MANIFEST_FILENAME,
    )
    if rows and audit_spectrum_rows(rows)["status"] == "PASS" and all((target / name).is_file() for name in required):
        manifest = json.loads((target / MANIFEST_FILENAME).read_text(encoding="utf-8"))
        manifest = _refresh_completed_metadata(target, manifest)
        manifest["invocation_root_calculation_count"] = ROOT_CALCULATION_COUNT - invocation_start_count
        manifest["completed_missing_only_numerical_outputs_unchanged"] = True
        return manifest
    constitutive, section_rows = prepare_constitutive_outputs(target)
    checkpoint_path = target / CHECKPOINT_FILENAME
    if checkpoint_path.is_file():
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        if checkpoint.get("contract_sha256") != contract_hash():
            raise RuntimeError("Checkpoint contract differs from RLB-2J.")
        point_records = [dict(item) for item in checkpoint.get("point_records", [])]
    else:
        point_records = []
    rows, benchmark = run_benchmarks(target, rows, point_records, constitutive, started_at)
    if not benchmark["production_run_permitted"]:
        return {
            "stage_id": STAGE_ID,
            "scientific_status": "NOT_EVALUATED",
            "workflow_status": "STOPPED_BY_ETA_GATE",
            "benchmark": benchmark,
        }
    rows = complete_missing_points(target, rows, point_records, constitutive, started_at)
    return finalize_outputs(
        target, rows, section_rows, constitutive, benchmark,
        point_records, started_at, workflow_started,
    )


def manifest_only(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    gate, rows = compute_constitutive_data()
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
    if args.plot_only:
        result = create_plots_from_csv(args.output_dir)
    elif args.manifest_only:
        result = manifest_only(args.output_dir)
    else:
        result = run_workflow(args.output_dir, missing_only=args.missing_only)
    print(json.dumps(_json_value(result), ensure_ascii=False, indent=2))
    return 1 if result.get("scientific_status") == "FAIL" else 0


if __name__ == "__main__":
    raise SystemExit(main())
