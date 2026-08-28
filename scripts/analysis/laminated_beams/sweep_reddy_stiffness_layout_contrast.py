"""RLB-2E stiffness-layout contrast map for two coupled Reddy beams.

The entry point implements one local ``frequency-map-v1`` instance.  It
varies only the scalar material contrast ``chi`` in four equal, zero-degree
plies.  Curves are independently sorted spectral positions, not tracked modal
branches.  Positions 1--8 are plotted; position 9 is a completeness guard.

The production physics is imported lazily.  Consequently ``--plot-only``
reads the completed CSV and renders the figure without importing laminate,
coupled-beam, determinant, SVD, or root-refinement code.
"""

from __future__ import annotations

import argparse
import copy
import csv
import ctypes
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
from typing import Any, Iterable, Mapping, Sequence

# The production CLI is intentionally single-threaded.  These variables must
# be fixed before NumPy (and therefore any BLAS runtime) is imported.
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


FloatArray = NDArray[np.float64]

STAGE_ID = "RLB-2E"
ALGORITHM_VERSION = "stiffness_layout_contrast_fast_plot_v2"
POLICY_ID = "frequency-map-v1"

CONFIG_BOTH_OUTER = "BOTH_OUTER_STIFF"
CONFIG_BOTH_INNER = "BOTH_INNER_STIFF"
CONFIG_ANTI_PHASE = "ANTI_PHASE"
CONFIGURATIONS = (
    CONFIG_BOTH_OUTER,
    CONFIG_BOTH_INNER,
    CONFIG_ANTI_PHASE,
)
OUTER_LAYOUT = ("H", "L", "L", "H")
INNER_LAYOUT = ("L", "H", "H", "L")
CONFIGURATION_LAYOUTS = {
    CONFIG_BOTH_OUTER: (OUTER_LAYOUT, OUTER_LAYOUT),
    CONFIG_BOTH_INNER: (INNER_LAYOUT, INNER_LAYOUT),
    CONFIG_ANTI_PHASE: (INNER_LAYOUT, OUTER_LAYOUT),
}

E1_0 = 1.1
E2_0 = 0.9
NU12_0 = 0.3
G0 = 1.0 / 2.6
RHO_0 = 1.0
DELTA = 0.1

MU = 0.0
TAU = 0.0
BETA_DEG = 30.0
L_REFERENCE = 1.0
L1 = 1.0
L2 = 1.0
L_TOTAL = 2.0
WIDTH = 0.20
THICKNESS = 0.05
PLY_THICKNESS = THICKNESS / 4.0
K = 5.0 / 6.0

CHI_MIN = 0.0
CHI_MAX = 0.8
CHI_STEP = 0.02
LOCAL_CHI_STEP = 0.005
K_PLOT = 8
K_GUARD = 9

A_REFERENCE = WIDTH * THICKNESS
IY_REFERENCE = WIDTH * THICKNESS**3 / 12.0
OMEGA_TO_OMEGA_SCALE = L_REFERENCE**2 * math.sqrt(
    RHO_0 * A_REFERENCE / IY_REFERENCE
)

MATRIX_RELATIVE_TOLERANCE = 1.0e-12
SYMMETRY_RELATIVE_TOLERANCE = 1.0e-12
REDUCED_PROPERTY_TOLERANCE = 1.0e-11
REDUCTION_ROUTE_TOLERANCE = 1.0e-11
ROOT_SINGULAR_RATIO_TOLERANCE = 1.0e-9
BOUNDARY_RESIDUAL_TOLERANCE = 1.0e-9
ARM_SWAP_RELATIVE_TOLERANCE = 1.0e-8
ROOT9_RIGHT_TAIL_OMEGA = 2.0
ANCHOR_BLOCK_WIDTH_OMEGA = 1.0
ANCHOR_SCAN_SPACING_OMEGA = 0.1
LOCAL_SCAN_MAX_SPACING_OMEGA = 0.1
MAX_ANCHOR_OMEGA = 500.0
ETA_LIMIT_SECONDS = 45.0 * 60.0
NEIGHBOUR_MAD_MULTIPLIER = 8.0
NEIGHBOUR_ABSOLUTE_TRIGGER = 1.0e-3

DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_stiffness_layout_contrast_sweep"
)
SPECTRUM_FILENAME = "spectrum_roots.csv"
SECTION_FILENAME = "section_properties.csv"
AUDIT_FILENAME = "neighbour_audit.csv"
BENCHMARK_FILENAME = "benchmark.json"
CHECKPOINT_FILENAME = "checkpoint.json"
PLOT_FILENAME = "lambda_vs_chi_three_configurations.png"
REPORT_FILENAME = "report.md"
MANIFEST_FILENAME = "run_manifest.json"

PRODUCTION_PHYSICS_PATHS = (
    "scripts/lib/reddy_symmetric_laminated_beam.py",
    "scripts/lib/reddy_symmetric_coupled_beams.py",
    "scripts/lib/reddy_inplane_geometry.py",
)

FREQUENCY_MAP_POLICY = {
    "frequency_map_policy": POLICY_ID,
    "calculation_mode": "fast_plot",
    "spectrum_semantics": "sorted_positions",
    "sweep_parameter": "chi",
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

SPECTRUM_FIELDS = (
    "row_id",
    "configuration_id",
    "chi",
    "grid_kind",
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
    "shared_chi0_anchor_reused",
    "shared_chi0_source_configuration",
    "is_canonical_plot_source",
    "supersedes_row_id",
    "repair_id",
    "roots_above_9_computed",
)


@dataclass(frozen=True)
class SectionObjects:
    """One four-ply laminate and its frozen one-dimensional reduction."""

    layout: tuple[str, str, str, str]
    laminate: Any
    properties: Any


@dataclass(frozen=True)
class PointSolution:
    """One accepted first-eight-plus-guard transaction."""

    configuration_id: str
    chi: float
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


class CountedProvider:
    """Count full matrix evaluations; each computes determinant and SVD."""

    def __init__(self, provider: Any) -> None:
        self.provider = provider
        self.evaluations = 0

    def __call__(self, omega: float) -> FloatArray:
        self.evaluations += 1
        return np.asarray(self.provider(float(omega)), dtype=float)


ROOT_CALCULATION_COUNT = 0


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _json_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _csv_value(value: Any) -> Any:
    if isinstance(value, (tuple, list, dict, np.ndarray)):
        return json.dumps(
            _json_value(value), ensure_ascii=False, separators=(",", ":")
        )
    if isinstance(value, np.generic):
        return value.item()
    return value


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".tmp")
    try:
        temporary.write_text(
            json.dumps(_json_value(payload), indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


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


def _atomic_write_csv(
    path: Path,
    rows: Iterable[Mapping[str, Any]],
    fieldnames: Sequence[str] | None = None,
) -> None:
    records = [dict(row) for row in rows]
    fields = list(fieldnames or (records[0].keys() if records else ()))
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + ".tmp")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="raise")
            writer.writeheader()
            for row in records:
                writer.writerow({key: _csv_value(row.get(key, "")) for key in fields})
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8", newline="") as stream:
        return [dict(row) for row in csv.DictReader(stream)]


def _git_state() -> dict[str, str]:
    def run(*arguments: str) -> str:
        result = subprocess.run(
            ["git", *arguments],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()

    return {
        "top_level": run("rev-parse", "--show-toplevel").replace("\\", "/"),
        "branch": run("branch", "--show-current"),
        "head": run("rev-parse", "HEAD"),
        "last_commit": run("log", "-1", "--oneline"),
        "status_short": run("status", "--short"),
    }


def _peak_rss_bytes() -> int:
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
        ok = psapi.GetProcessMemoryInfo(
            kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb
        )
        return int(counters.PeakWorkingSetSize) if ok else 0
    try:  # pragma: no cover - Windows is the production platform.
        import resource

        value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        return value if sys.platform == "darwin" else 1024 * value
    except (ImportError, OSError):
        return 0


def chi_grid() -> FloatArray:
    """Return the exact 41-point project grid."""

    return np.asarray(
        sweep_grid_policy.parameter_grid(CHI_MIN, CHI_MAX, CHI_STEP),
        dtype=float,
    )


def omega_to_Omega(omega: float) -> float:
    return float(omega) * OMEGA_TO_OMEGA_SCALE


def Omega_to_Lambda(Omega: float) -> float:
    value = float(Omega)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError("Omega must be finite and nonnegative.")
    return math.sqrt(value)


def base_material_contract() -> dict[str, float]:
    return {
        "delta": DELTA,
        "E1": E1_0,
        "E2": E2_0,
        "nu12": NU12_0,
        "G12": G0,
        "G13": G0,
        "G23": G0,
        "rho": RHO_0,
    }


def _physics_modules() -> tuple[Any, Any]:
    from scripts.lib import reddy_symmetric_coupled_beams as coupled
    from scripts.lib import reddy_symmetric_laminated_beam as beam

    return beam, coupled


def _root_tools() -> Any:
    from scripts.analysis.laminated_beams import (
        pilot_reddy_symmetric_coupled_beams_beta0 as tools,
    )

    return tools


def contrast_materials(chi: float) -> dict[str, Any]:
    """Build H/L from M0 by one scalar stiffness multiplier."""

    value = float(chi)
    if not math.isfinite(value) or not CHI_MIN <= value <= CHI_MAX:
        raise ValueError("chi must lie in [0, 0.8].")
    beam, _coupled = _physics_modules()
    materials: dict[str, Any] = {}
    for label, factor in (("H", 1.0 + value), ("L", 1.0 - value)):
        materials[label] = beam.OrthotropicLamina(
            E1=factor * E1_0,
            E2=factor * E2_0,
            nu12=NU12_0,
            G12=factor * G0,
            G13=factor * G0,
            G23=factor * G0,
            rho=RHO_0,
            name=f"RLB-2E {label}(chi={value:.2f})",
        )
    return materials


def build_layout_section(layout: Sequence[str], chi: float) -> SectionObjects:
    """Build one palindromic four-ply, all-zero-degree section."""

    labels = tuple(str(item) for item in layout)
    if labels not in (OUTER_LAYOUT, INNER_LAYOUT):
        raise ValueError(f"Unknown RLB-2E layout: {labels!r}.")
    beam, _coupled = _physics_modules()
    materials = contrast_materials(chi)
    laminate = beam.integrate_laminate(
        [
            beam.Ply(
                materials[label],
                angle_deg=0.0,
                thickness=PLY_THICKNESS,
                label=label,
            )
            for label in labels
        ]
    )
    properties = beam.reduce_to_beam_properties(
        laminate,
        width=WIDTH,
        K=K,
        symmetry_tolerance=SYMMETRY_RELATIVE_TOLERANCE,
        reduction_tolerance=REDUCTION_ROUTE_TOLERANCE,
    )
    return SectionObjects(labels, laminate, properties)


def build_configuration_sections(
    configuration_id: str, chi: float
) -> tuple[SectionObjects, SectionObjects]:
    try:
        layouts = CONFIGURATION_LAYOUTS[str(configuration_id)]
    except KeyError as exc:
        raise ValueError(f"Unknown configuration: {configuration_id!r}.") from exc
    return tuple(build_layout_section(layout, chi) for layout in layouts)  # type: ignore[return-value]


def _relative(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / max(
        abs(float(left)), abs(float(right)), np.finfo(float).tiny
    )


def _matrix_relative(left: Any, right: Any) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    scale = max(
        float(np.linalg.norm(a, ord="fro")),
        float(np.linalg.norm(b, ord="fro")),
        np.finfo(float).tiny,
    )
    return float(np.linalg.norm(a - b, ord="fro") / scale)


def _scaled_B(laminate: Any) -> float:
    scale = max(
        float(np.linalg.norm(laminate.A, ord="fro")) * laminate.thickness,
        np.finfo(float).tiny,
    )
    return float(np.linalg.norm(laminate.B, ord="fro") / scale)


def _scaled_I1(laminate: Any) -> float:
    scale = max(abs(float(laminate.I0)) * laminate.thickness, np.finfo(float).tiny)
    return abs(float(laminate.I1)) / scale


def _reduction_max_residual(properties: Any) -> float:
    routes = (
        properties.axial_reduction,
        properties.bending_reduction,
        properties.shear_reduction_before_K,
    )
    return max(float(route.relative_difference) for route in routes)


def _baseline_section() -> SectionObjects:
    beam, _coupled = _physics_modules()
    material = beam.OrthotropicLamina(
        E1=E1_0,
        E2=E2_0,
        nu12=NU12_0,
        G12=G0,
        G13=G0,
        G23=G0,
        rho=RHO_0,
        name="RLB-2E M0",
    )
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


def constitutive_gate(
    values: Sequence[float] = (0.0, 0.4, 0.8),
) -> dict[str, Any]:
    """Check the complete analytic layout contrast before any root search."""

    baseline = _baseline_section()
    records: list[dict[str, Any]] = []
    maximums = {
        "D_matrix_formula_relative": 0.0,
        "Dbeam_formula_relative": 0.0,
        "A_matrix_layout_relative": 0.0,
        "shear_matrix_layout_relative": 0.0,
        "reduced_invariant_relative": 0.0,
        "symmetry_relative": 0.0,
        "reduction_route_relative": 0.0,
    }
    passed = True
    for raw_chi in values:
        chi = float(raw_chi)
        outer = build_layout_section(OUTER_LAYOUT, chi)
        inner = build_layout_section(INNER_LAYOUT, chi)
        expected_outer = 1.0 + 0.75 * chi
        expected_inner = 1.0 - 0.75 * chi
        D_outer_matrix_residual = _matrix_relative(
            outer.laminate.D, expected_outer * baseline.laminate.D
        )
        D_inner_matrix_residual = _matrix_relative(
            inner.laminate.D, expected_inner * baseline.laminate.D
        )
        D_outer_beam_residual = _relative(
            outer.properties.D / baseline.properties.D, expected_outer
        )
        D_inner_beam_residual = _relative(
            inner.properties.D / baseline.properties.D, expected_inner
        )
        sum_matrix_residual = _matrix_relative(
            outer.laminate.D + inner.laminate.D,
            2.0 * baseline.laminate.D,
        )
        A_residual = _matrix_relative(outer.laminate.A, inner.laminate.A)
        shear_residual = _matrix_relative(
            outer.laminate.shear, inner.laminate.shear
        )
        reduced_residual = max(
            _relative(
                getattr(outer.properties, name), getattr(inner.properties, name)
            )
            for name in ("A", "S", "m", "J")
        )
        mass_residual = max(
            _relative(outer.laminate.I0, inner.laminate.I0),
            _relative(outer.laminate.I2, inner.laminate.I2),
        )
        symmetry_residual = max(
            _scaled_B(outer.laminate),
            _scaled_B(inner.laminate),
            _scaled_I1(outer.laminate),
            _scaled_I1(inner.laminate),
        )
        route_residual = max(
            _reduction_max_residual(outer.properties),
            _reduction_max_residual(inner.properties),
        )
        moduli_positive = all(
            getattr(material, field) > 0.0
            for material in contrast_materials(chi).values()
            for field in ("E1", "E2", "G12", "G13", "G23")
        )
        point_pass = bool(
            moduli_positive
            and max(
                D_outer_matrix_residual,
                D_inner_matrix_residual,
                sum_matrix_residual,
            )
            <= MATRIX_RELATIVE_TOLERANCE
            and max(D_outer_beam_residual, D_inner_beam_residual)
            <= REDUCED_PROPERTY_TOLERANCE
            and A_residual <= MATRIX_RELATIVE_TOLERANCE
            and shear_residual <= MATRIX_RELATIVE_TOLERANCE
            and mass_residual <= MATRIX_RELATIVE_TOLERANCE
            and reduced_residual <= REDUCED_PROPERTY_TOLERANCE
            and symmetry_residual <= SYMMETRY_RELATIVE_TOLERANCE
            and route_residual <= REDUCTION_ROUTE_TOLERANCE
        )
        passed = passed and point_pass
        maximums["D_matrix_formula_relative"] = max(
            maximums["D_matrix_formula_relative"],
            D_outer_matrix_residual,
            D_inner_matrix_residual,
            sum_matrix_residual,
        )
        maximums["Dbeam_formula_relative"] = max(
            maximums["Dbeam_formula_relative"],
            D_outer_beam_residual,
            D_inner_beam_residual,
        )
        maximums["A_matrix_layout_relative"] = max(
            maximums["A_matrix_layout_relative"], A_residual
        )
        maximums["shear_matrix_layout_relative"] = max(
            maximums["shear_matrix_layout_relative"], shear_residual
        )
        maximums["reduced_invariant_relative"] = max(
            maximums["reduced_invariant_relative"], reduced_residual, mass_residual
        )
        maximums["symmetry_relative"] = max(
            maximums["symmetry_relative"], symmetry_residual
        )
        maximums["reduction_route_relative"] = max(
            maximums["reduction_route_relative"], route_residual
        )
        records.append(
            {
                "chi": chi,
                "status": "PASS" if point_pass else "FAIL",
                "D_outer_matrix_over_D0_expected": expected_outer,
                "D_inner_matrix_over_D0_expected": expected_inner,
                "D_outer_beam_over_Dbeam0": (
                    outer.properties.D / baseline.properties.D
                ),
                "D_inner_beam_over_Dbeam0": (
                    inner.properties.D / baseline.properties.D
                ),
                "D_outer_matrix_formula_relative": D_outer_matrix_residual,
                "D_inner_matrix_formula_relative": D_inner_matrix_residual,
                "D_sum_formula_relative": sum_matrix_residual,
                "D_outer_beam_formula_relative": D_outer_beam_residual,
                "D_inner_beam_formula_relative": D_inner_beam_residual,
                "A_layout_relative": A_residual,
                "shear_layout_relative": shear_residual,
                "mass_layout_relative": mass_residual,
                "reduced_invariant_relative": reduced_residual,
                "symmetry_relative": symmetry_residual,
                "reduction_route_relative": route_residual,
                "all_moduli_positive": moduli_positive,
            }
        )
    return {
        "status": "PASS" if passed else "FAIL",
        "passed": passed,
        "checks": records,
        "maximum_residuals": maximums,
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
    """Return all 246 configuration/arm/base-grid section rows."""

    baseline = _baseline_section()
    rows: list[dict[str, Any]] = []
    for configuration_id in CONFIGURATIONS:
        for chi in chi_grid():
            sections = build_configuration_sections(configuration_id, float(chi))
            for arm_id, section in enumerate(sections, start=1):
                layout_kind = "OUTER_STIFF" if section.layout == OUTER_LAYOUT else "INNER_STIFF"
                expected_ratio = 1.0 + 0.75 * float(chi)
                if layout_kind == "INNER_STIFF":
                    expected_ratio = 1.0 - 0.75 * float(chi)
                props = section.properties
                laminate = section.laminate
                actual_ratio = props.D / baseline.properties.D
                rows.append(
                    {
                        "configuration_id": configuration_id,
                        "arm_id": arm_id,
                        "chi": float(chi),
                        "grid_kind": "BASE",
                        "layout_kind": layout_kind,
                        "stack_bottom_to_top": section.layout,
                        "material_labels": [ply.label for ply in laminate.plies],
                        "ply_angles_deg": [ply.angle_deg for ply in laminate.plies],
                        "ply_thicknesses": [ply.thickness for ply in laminate.plies],
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
                        "B_relative": _scaled_B(laminate),
                        "I1_relative": _scaled_I1(laminate),
                        "Dbeam0": baseline.properties.D,
                        "D_beam_over_Dbeam0": actual_ratio,
                        "analytic_D_beam_over_Dbeam0": expected_ratio,
                        "analytic_D_beam_ratio_residual": _relative(
                            actual_ratio, expected_ratio
                        ),
                        "reduction_route_max_relative": _reduction_max_residual(
                            props
                        ),
                    }
                )
    return rows


def make_matrix_provider(
    configuration_id: str, chi: float
) -> tuple[Any, dict[str, Any]]:
    """Assemble only the frozen RLB arm maps and frozen rigid joint."""

    _beam, coupled = _physics_modules()
    arm1, arm2 = build_configuration_sections(configuration_id, chi)
    joint = np.asarray(coupled.joint_matrix(math.radians(BETA_DEG)), dtype=float)
    identical = bool(
        arm1.layout == arm2.layout
        and all(
            _relative(getattr(arm1.properties, name), getattr(arm2.properties, name))
            <= 8.0 * np.finfo(float).eps
            for name in ("A", "D", "S", "m", "J")
        )
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
        raise RuntimeError("RLB-2E provider differs from the public coupled builder.")
    return provider, {
        "configuration_id": configuration_id,
        "chi": float(chi),
        "beta_deg": BETA_DEG,
        "identical_arm_map_reused": identical,
        "cached_vs_public_builder_max_abs": direct_residual,
        "production_modules_only": True,
    }


def hold_secant_predictor(
    parameter: float,
    previous_parameter: float,
    previous_roots: Sequence[float],
    second_previous_parameter: float | None = None,
    second_previous_roots: Sequence[float] | None = None,
) -> FloatArray:
    """Return hold/secant locators; never an accepted root value."""

    last = np.asarray(previous_roots, dtype=float)
    if last.shape != (K_GUARD,):
        raise ValueError("The predictor requires exactly roots 1...9.")
    if second_previous_parameter is None or second_previous_roots is None:
        return last.copy()
    older = np.asarray(second_previous_roots, dtype=float)
    if older.shape != last.shape:
        raise ValueError("Previous root arrays must have identical shape.")
    denominator = float(previous_parameter) - float(second_previous_parameter)
    if denominator == 0.0:
        raise ValueError("Predictor parameter points must be distinct.")
    ratio = (float(parameter) - float(previous_parameter)) / denominator
    return last + ratio * (last - older)


def local_search_windows(
    predicted_roots: Sequence[float], *, guard_right_width: float = 1.2
) -> tuple[tuple[float, float], ...]:
    """Partition the expected first-nine range by predicted midpoints."""

    roots = np.asarray(predicted_roots, dtype=float)
    if roots.shape != (K_GUARD,) or not np.all(np.isfinite(roots)):
        raise ValueError("Nine finite predicted roots are required.")
    if np.any(np.diff(roots) <= 0.0):
        raise ValueError("Predicted roots must be strictly increasing.")
    midpoints = 0.5 * (roots[:-1] + roots[1:])
    first_width = max(2.0, 0.5 * (roots[1] - roots[0]))
    # Root 9 is a completeness guard.  Do not widen its locator by half of a
    # potentially large root-8/root-9 gap: that can unnecessarily enter the
    # root-10 neighbourhood.  Continuation error is handled by the existing
    # bounded fallback, while the guard keeps its declared right margin.
    last_width = float(guard_right_width)
    if not math.isfinite(last_width) or last_width <= 0.0:
        raise ValueError("guard_right_width must be finite and positive.")
    edges = np.concatenate(
        ([max(1.0e-8, roots[0] - first_width)], midpoints, [roots[-1] + last_width])
    )
    return tuple((float(edges[index]), float(edges[index + 1])) for index in range(K_GUARD))


def _rlb2e_search_policy() -> Any:
    """Adapt frozen RLB thresholds to the local first-eight-plus-guard count."""

    tools = _root_tools()
    # The frozen predecessor policy intentionally validates a 12+1 inventory
    # in its constructor.  This analysis-local shallow copy changes only the
    # inventory length; all physics-independent quality thresholds stay byte-
    # for-byte equal to the predecessor and are asserted in targeted tests.
    policy = copy.copy(tools.SearchPolicy())
    object.__setattr__(policy, "requested_roots", K_PLOT)
    object.__setattr__(policy, "guard_roots", 1)
    if policy.required_slots != K_GUARD:
        raise RuntimeError("RLB-2E search policy must retain exactly roots 1...9.")
    return policy


def _scan_policy(base: Any, left: float, right: float, spacing: float) -> Any:
    width = float(right) - float(left)
    points = max(33, int(math.ceil(width / float(spacing))) + 1)
    policy = copy.copy(base)
    object.__setattr__(policy, "omega_bar_min", float(left))
    object.__setattr__(policy, "omega_bar_max", float(right))
    object.__setattr__(policy, "primary_scan_points", points)
    object.__setattr__(
        policy,
        "verification_scan_points",
        max(points + 1, 2 * points - 1),
    )
    return policy


def _unresolved_candidates_below(
    rejected: Sequence[Any], canonical: Sequence[Any], guard_Omega: float, policy: Any
) -> int:
    failure_reasons = {
        "NONFINITE_MATRIX",
        "NONFINITE_DETERMINANT_INTERVAL",
        "BRENT_FAILURE",
        "MINIMIZER_EXCEPTION",
        "MINIMIZER_FAILURE",
        "NESTED_MINIMIZER_EXCEPTION",
    }
    count = 0
    for candidate in rejected:
        if float(candidate.interval_left_bar) > float(guard_Omega):
            continue
        if candidate.rejection_reason == "DUPLICATE_DETECTION_SAME_ROOT":
            continue
        if any(
            abs(float(candidate.omega_bar) - float(root.omega_bar))
            <= policy.dedup_atol_bar
            + policy.dedup_rtol
            * max(abs(float(candidate.omega_bar)), abs(float(root.omega_bar)))
            for root in canonical
        ):
            continue
        if candidate.rejection_reason in failure_reasons:
            count += 1
            continue
        if candidate.diagnostics.scaled_sigma_ratio <= policy.sigma_prefilter:
            count += 1
    return count


def _canonical_slots(raw_candidates: Sequence[Any], policy: Any) -> tuple[list[Any], list[Any], list[Any]]:
    tools = _root_tools()
    canonical, rejected = tools._merge_candidates(raw_candidates, policy)
    _events, slots = tools._events_and_slots(canonical, policy)
    return list(canonical), list(rejected), list(slots)


def _scan_candidate_window(
    counted: CountedProvider,
    base_policy: Any,
    left: float,
    right: float,
    *,
    scan_id: str,
    phases: Sequence[float] = (0.0,),
    spacing: float = LOCAL_SCAN_MAX_SPACING_OMEGA,
) -> list[Any]:
    tools = _root_tools()
    policy = _scan_policy(base_policy, left, right, spacing)
    candidates, _evaluator = tools._scan_candidates(
        counted,
        OMEGA_TO_OMEGA_SCALE,
        policy,
        case_id=scan_id,
        builder_id="RLB2E_FROZEN_RLB",
        scan_id=scan_id,
        points=policy.primary_scan_points,
        phases=tuple(float(value) for value in phases),
    )
    return candidates


def _progressive_anchor_candidates(
    counted: CountedProvider,
    base_policy: Any,
    *,
    solve_id: str,
) -> tuple[list[Any], list[Any], list[Any], float, int]:
    raw: list[Any] = []
    lower = max(1.0e-8, float(base_policy.omega_bar_min))
    upper = lower
    segment_index = 0
    local_refinements = 0
    while upper < MAX_ANCHOR_OMEGA:
        segment_index += 1
        segment_left = lower if segment_index == 1 else max(lower, upper - 0.25)
        segment_right = min(MAX_ANCHOR_OMEGA, upper + ANCHOR_BLOCK_WIDTH_OMEGA)
        raw.extend(
            _scan_candidate_window(
                counted,
                base_policy,
                segment_left,
                segment_right,
                scan_id=f"{solve_id}_anchor_{segment_index:03d}",
                phases=(0.0, 0.5),
                spacing=ANCHOR_SCAN_SPACING_OMEGA,
            )
        )
        local_refinements += 1
        canonical, rejected, slots = _canonical_slots(raw, base_policy)
        upper = segment_right
        if len(slots) >= K_GUARD:
            guard = float(slots[K_GUARD - 1].event.omega_bar)
            search_right = max(upper, guard + ROOT9_RIGHT_TAIL_OMEGA)
            # One matrix diagnostic at the declared right boundary establishes
            # the guard margin without locating/refining root 10.
            counted(search_right / OMEGA_TO_OMEGA_SCALE)
            return canonical, rejected, slots, search_right, local_refinements
    raise RuntimeError(f"{solve_id}: root 9 was not found below Omega=500.")


def _local_candidates(
    counted: CountedProvider,
    base_policy: Any,
    predicted: Sequence[float],
    *,
    solve_id: str,
    dense: bool = False,
    dense_positions: Sequence[int] | None = None,
    guard_right_width: float = 1.2,
) -> tuple[list[Any], list[Any], list[Any], float, int]:
    raw: list[Any] = []
    windows = local_search_windows(predicted, guard_right_width=guard_right_width)
    dense_set = (
        set(range(1, K_GUARD + 1))
        if dense and dense_positions is None
        else {int(value) for value in (dense_positions or ())}
    )
    if any(position < 1 or position > K_GUARD for position in dense_set):
        raise ValueError("Dense root positions must lie in 1...9.")
    for index, (left, right) in enumerate(windows, start=1):
        raw.extend(
            _scan_candidate_window(
                counted,
                base_policy,
                left,
                right,
                scan_id=f"{solve_id}_local_{index:02d}",
                phases=(0.0,),
                spacing=(
                    0.05 if index in dense_set else LOCAL_SCAN_MAX_SPACING_OMEGA
                ),
            )
        )
    canonical, rejected, slots = _canonical_slots(raw, base_policy)
    search_right = windows[-1][1]
    if len(slots) >= K_GUARD:
        guard = float(slots[K_GUARD - 1].event.omega_bar)
        search_right = guard + ROOT9_RIGHT_TAIL_OMEGA
        # Establish the declared guard-to-boundary margin without scanning or
        # refining a tenth root neighbourhood.
        counted(search_right / OMEGA_TO_OMEGA_SCALE)
    return canonical, rejected, slots, search_right, len(windows)


def _point_is_acceptable(
    canonical: Sequence[Any],
    rejected: Sequence[Any],
    slots: Sequence[Any],
    search_right: float,
    policy: Any,
) -> tuple[bool, dict[str, Any]]:
    tools = _root_tools()
    first = list(slots[:K_GUARD])
    strict_order = bool(
        len(first) == K_GUARD
        and all(
            float(first[index].event.omega_bar)
            < float(first[index + 1].event.omega_bar)
            for index in range(K_GUARD - 1)
        )
    )
    qualities: list[bool] = []
    if len(first) == K_GUARD:
        qualities = [
            bool(tools._candidate_quality(slot.event.candidate.diagnostics, policy)[0])
            for slot in first
        ]
    guard = float(first[-1].event.omega_bar) if first else math.inf
    accepted_above_guard = []
    if len(first) == K_GUARD:
        tolerance = policy.dedup_atol_bar + policy.dedup_rtol * abs(guard)
        accepted_above_guard = [
            candidate
            for candidate in canonical
            if float(candidate.omega_bar) > guard + tolerance
        ]
    unresolved = (
        _unresolved_candidates_below(rejected, canonical, guard, policy)
        if len(first) == K_GUARD
        else -1
    )
    margin = float(search_right) - guard
    passed = bool(
        len(first) == K_GUARD
        and len(slots) == K_GUARD
        and not accepted_above_guard
        and strict_order
        and all(qualities)
        and unresolved == 0
        and margin >= ROOT9_RIGHT_TAIL_OMEGA - 1.0e-10
    )
    return passed, {
        "root_count": len(first),
        "strict_order": strict_order,
        "all_quality_pass": bool(qualities and all(qualities)),
        "unresolved_candidates_below_root9": unresolved,
        "root9_right_margin_Omega": margin,
        "candidate_count_total": len(canonical),
        "accepted_candidates_above_root9": len(accepted_above_guard),
        "retained_slots_above_root9": max(0, len(slots) - K_GUARD),
        "roots_above_9_computed": bool(
            accepted_above_guard or len(slots) > K_GUARD
        ),
    }


def _root_rows(
    configuration_id: str,
    chi: float,
    slots: Sequence[Any],
    *,
    solve_id: str,
    transaction_id: str,
    solve_mode: str,
    fallback_used: bool,
    predicted: Sequence[float] | None,
    search_right: float,
    unresolved: int,
    grid_kind: str = "BASE",
    repair_id: str = "",
    shared_chi0_anchor_reused: bool = False,
) -> tuple[dict[str, Any], ...]:
    windows = None if predicted is None else local_search_windows(predicted)
    result: list[dict[str, Any]] = []
    guard_Omega = float(slots[K_GUARD - 1].event.omega_bar)
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
            f"{configuration_id}__{float(chi):.6f}__{grid_kind}__"
            f"p{position:02d}__{repair_id or 'base'}"
        )
        result.append(
            {
                "row_id": row_id,
                "configuration_id": configuration_id,
                "chi": float(chi),
                "grid_kind": grid_kind,
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
                "root9_right_margin_Omega": search_right - guard_Omega,
                "solve_mode": solve_mode,
                "fallback_used": fallback_used,
                "quality_status": "PASS",
                "point_status": "PASS",
                "shared_chi0_anchor_reused": shared_chi0_anchor_reused,
                "shared_chi0_source_configuration": (
                    CONFIG_BOTH_OUTER if shared_chi0_anchor_reused else ""
                ),
                "is_canonical_plot_source": True,
                "supersedes_row_id": "",
                "repair_id": repair_id,
                "roots_above_9_computed": False,
            }
        )
    return tuple(result)


def solve_point(
    configuration_id: str,
    chi: float,
    *,
    previous: tuple[float, Sequence[float]] | None = None,
    second_previous: tuple[float, Sequence[float]] | None = None,
    force_anchor: bool = False,
    dense_local: bool = False,
    dense_positions: Sequence[int] | None = None,
    grid_kind: str = "BASE",
    repair_id: str = "",
    guard_locator_right_width: float = 1.2,
) -> PointSolution:
    """Solve one point from the characteristic matrix, never interpolation."""

    global ROOT_CALCULATION_COUNT
    ROOT_CALCULATION_COUNT += 1
    started = time.perf_counter()
    solve_id = f"{configuration_id}__chi_{float(chi):.6f}"
    transaction_id = hashlib.sha256(
        f"{STAGE_ID}|{solve_id}|{grid_kind}|{repair_id}".encode("utf-8")
    ).hexdigest()[:20].upper()
    provider, _metadata = make_matrix_provider(configuration_id, chi)
    counted = CountedProvider(provider)
    base_policy = _rlb2e_search_policy()
    predicted: FloatArray | None = None
    fallback_used = False
    if not force_anchor and previous is not None:
        predicted = hold_secant_predictor(
            chi,
            previous[0],
            previous[1],
            None if second_previous is None else second_previous[0],
            None if second_previous is None else second_previous[1],
        )
        # The map semantics are independently sorted positions.  A secant
        # locator may reverse two neighbouring predictions near an approach;
        # sorting the locators restores disjoint numerical windows without
        # claiming modal identity or changing any final root.
        predicted = np.sort(predicted)
        try:
            canonical, rejected, slots, search_right, refinements = _local_candidates(
                counted,
                base_policy,
                predicted,
                solve_id=solve_id,
                dense=dense_local,
                dense_positions=dense_positions,
                guard_right_width=guard_locator_right_width,
            )
            accepted, quality = _point_is_acceptable(
                canonical, rejected, slots, search_right, base_policy
            )
        except (ValueError, RuntimeError, ArithmeticError, np.linalg.LinAlgError):
            accepted = False
            quality = {}
        if not accepted:
            fallback_used = True
            canonical, rejected, slots, search_right, refinements = (
                _progressive_anchor_candidates(
                    counted, base_policy, solve_id=solve_id + "_fallback"
                )
            )
            accepted, quality = _point_is_acceptable(
                canonical, rejected, slots, search_right, base_policy
            )
            solve_mode = "BOUNDED_GLOBAL_RECOVERY"
        else:
            solve_mode = "FAST_LOCAL"
    else:
        canonical, rejected, slots, search_right, refinements = (
            _progressive_anchor_candidates(counted, base_policy, solve_id=solve_id)
        )
        accepted, quality = _point_is_acceptable(
            canonical, rejected, slots, search_right, base_policy
        )
        solve_mode = "PROGRESSIVE_ANCHOR"
    if not accepted:
        raise RuntimeError(f"{solve_id}: first-eight-plus-root9 quality failed: {quality}")
    rows = _root_rows(
        configuration_id,
        chi,
        slots,
        solve_id=solve_id,
        transaction_id=transaction_id,
        solve_mode=solve_mode,
        fallback_used=fallback_used,
        # A rejected/non-monotone predictor remains only a fallback trigger.
        # The accepted global-recovery roots must export their actual scan
        # intervals and must never rebuild locator windows from that predictor.
        predicted=None if fallback_used else predicted,
        search_right=search_right,
        unresolved=int(quality["unresolved_candidates_below_root9"]),
        grid_kind=grid_kind,
        repair_id=repair_id,
    )
    return PointSolution(
        configuration_id=configuration_id,
        chi=float(chi),
        rows=rows,
        wall_time_seconds=time.perf_counter() - started,
        peak_rss_bytes=_peak_rss_bytes(),
        determinant_evaluations=counted.evaluations,
        sigma_evaluations=counted.evaluations,
        search_left_Omega=(
            1.0e-8
            if predicted is None or fallback_used
            else local_search_windows(predicted)[0][0]
        ),
        search_right_Omega=search_right,
        local_refinements=refinements,
        solve_mode=solve_mode,
        fallback_used=fallback_used,
        unresolved_candidates_below_root9=int(
            quality["unresolved_candidates_below_root9"]
        ),
        candidate_count_total=int(quality["candidate_count_total"]),
        accepted_candidates_above_root9=int(
            quality["accepted_candidates_above_root9"]
        ),
        retained_slots_above_root9=int(quality["retained_slots_above_root9"]),
        roots_above_9_computed=bool(quality["roots_above_9_computed"]),
    )


def _base_group_index(rows: Sequence[Mapping[str, Any]]) -> dict[tuple[str, float], list[Mapping[str, Any]]]:
    groups: dict[tuple[str, float], list[Mapping[str, Any]]] = {}
    for row in rows:
        if str(row.get("grid_kind")) != "BASE":
            continue
        key = (str(row["configuration_id"]), round(float(row["chi"]), 10))
        groups.setdefault(key, []).append(row)
    return groups


def _as_bool(value: Any) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def _base_group_has_exact_positions(group: Sequence[Mapping[str, Any]]) -> bool:
    try:
        positions = [int(row["sorted_position"]) for row in group]
    except (KeyError, TypeError, ValueError):
        return False
    return len(positions) == K_GUARD and sorted(positions) == list(
        range(1, K_GUARD + 1)
    )


def _base_group_is_acceptable(group: Sequence[Mapping[str, Any]]) -> bool:
    """Return whether a saved BASE transaction may be reused without solving."""

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
) -> dict[tuple[str, float], list[Mapping[str, Any]]]:
    return {
        key: group
        for key, group in _base_group_index(rows).items()
        if _base_group_is_acceptable(group)
    }


def audit_spectrum_rows(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Validate the immutable 1107-row BASE contract and repair additions."""

    groups = _base_group_index(rows)
    complete_groups = _complete_base_group_index(rows)
    duplicates: list[str] = []
    incomplete: list[str] = []
    physical_failures: list[str] = []
    for key, group in groups.items():
        positions = [int(row["sorted_position"]) for row in group]
        if len(positions) != len(set(positions)):
            duplicates.append(f"{key[0]}:{key[1]:.2f}")
        if not _base_group_has_exact_positions(group):
            incomplete.append(f"{key[0]}:{key[1]:.2f}")
            continue
        if not _base_group_is_acceptable(group):
            physical_failures.append(f"{key[0]}:{key[1]:.2f}")
    expected = {
        (configuration, round(float(chi), 10))
        for configuration in CONFIGURATIONS
        for chi in chi_grid()
    }
    missing = sorted(expected - set(groups))
    extra = sorted(set(groups) - expected)
    base_rows = [row for row in rows if str(row.get("grid_kind")) == "BASE"]
    roots_above_guard = [
        row for row in rows if int(row["sorted_position"]) > K_GUARD
    ]
    duplicate_row_ids = len(rows) - len({str(row["row_id"]) for row in rows})
    canonical_counts: dict[tuple[str, float, int], int] = {}
    for row in rows:
        if _as_bool(row.get("is_canonical_plot_source", False)):
            key = (
                str(row["configuration_id"]),
                round(float(row["chi"]), 10),
                int(row["sorted_position"]),
            )
            canonical_counts[key] = canonical_counts.get(key, 0) + 1
    canonical_failures = sum(value != 1 for value in canonical_counts.values())
    expected_canonical = len(CONFIGURATIONS) * len(chi_grid()) * K_GUARD
    canonical_failures += max(0, expected_canonical - len(canonical_counts))
    return {
        "status": "PASS"
        if not duplicates
        and not incomplete
        and not missing
        and not extra
        and not roots_above_guard
        and not physical_failures
        and duplicate_row_ids == 0
        and canonical_failures == 0
        and len(base_rows) == len(CONFIGURATIONS) * len(chi_grid()) * K_GUARD
        else "FAIL",
        "base_group_count": len(complete_groups),
        "base_row_count": len(base_rows),
        "missing_groups": [f"{item[0]}:{item[1]:.2f}" for item in missing],
        "extra_groups": [f"{item[0]}:{item[1]:.2f}" for item in extra],
        "duplicate_groups": duplicates,
        "incomplete_groups": incomplete,
        "physical_quality_failures": physical_failures,
        "duplicate_row_id_count": duplicate_row_ids,
        "canonical_source_failure_count": canonical_failures,
        "roots_above_guard_count": len(roots_above_guard),
    }


def _write_point_transaction(
    output_dir: Path,
    existing_rows: Sequence[Mapping[str, Any]],
    solution: PointSolution,
) -> list[dict[str, Any]]:
    rows = [dict(row) for row in existing_rows]
    key = (solution.configuration_id, round(solution.chi, 10), "BASE")
    complete = [
        row
        for row in rows
        if (
            str(row["configuration_id"]),
            round(float(row["chi"]), 10),
            str(row["grid_kind"]),
        )
        == key
    ]
    if complete:
        if _base_group_is_acceptable(complete):
            return rows
        # A structurally or physically damaged transaction is not complete.
        # Replace only that group after a fresh nine-root solution has passed.
        rows = [
            row
            for row in rows
            if not (
                str(row["configuration_id"]) == solution.configuration_id
                and round(float(row["chi"]), 10) == round(solution.chi, 10)
                and str(row["grid_kind"]) == "BASE"
            )
        ]
    rows.extend(dict(row) for row in solution.rows)
    rows.sort(
        key=lambda row: (
            CONFIGURATIONS.index(str(row["configuration_id"])),
            float(row["chi"]),
            0 if str(row["grid_kind"]) == "BASE" else 1,
            int(row["sorted_position"]),
        )
    )
    _atomic_write_csv(output_dir / SPECTRUM_FILENAME, rows, SPECTRUM_FIELDS)
    return rows


def _rows_for_roots(
    rows: Sequence[Mapping[str, Any]], configuration_id: str, chi: float
) -> FloatArray:
    group = [
        row
        for row in rows
        if str(row["configuration_id"]) == configuration_id
        and round(float(row["chi"]), 10) == round(float(chi), 10)
        and _as_bool(row.get("is_canonical_plot_source", True))
    ]
    group.sort(key=lambda row: int(row["sorted_position"]))
    if [int(row["sorted_position"]) for row in group] != list(
        range(1, K_GUARD + 1)
    ):
        raise RuntimeError(f"Incomplete root group {configuration_id}, chi={chi}.")
    return np.asarray([float(row["Omega"]) for row in group], dtype=float)


def _checkpoint_payload(
    rows: Sequence[Mapping[str, Any]],
    point_records: Sequence[Mapping[str, Any]],
    *,
    constitutive: Mapping[str, Any],
    started_at: str,
    benchmark_status: str,
) -> dict[str, Any]:
    audit = audit_spectrum_rows(rows)
    groups = _complete_base_group_index(rows)
    expected = [
        (configuration, round(float(chi), 10))
        for configuration in CONFIGURATIONS
        for chi in chi_grid()
    ]
    missing = [
        {"configuration_id": configuration, "chi": chi}
        for configuration, chi in expected
        if (configuration, chi) not in groups
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
        "completed_base_rows": audit["base_row_count"],
        "missing_points": missing,
        "failed_points": [],
        "terminal_unresolved_points": [],
        "last_completed_parameter": (
            point_records[-1].get("chi") if point_records else None
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
    configuration_id: str,
    chi: float,
    error: BaseException,
) -> None:
    """Persist a failed transaction without publishing an incomplete group."""

    payload = _checkpoint_payload(
        rows,
        point_records,
        constitutive=constitutive,
        started_at=started_at,
        benchmark_status=benchmark_status,
    )
    payload["failed_points"] = [
        {
            "configuration_id": configuration_id,
            "chi": float(chi),
            "error": f"{type(error).__name__}: {error}",
            "recorded_at_utc": _utc_now(),
        }
    ]
    _atomic_write_json(Path(output_dir) / CHECKPOINT_FILENAME, payload)


def contract_payload() -> dict[str, Any]:
    return {
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "frequency_map_policy": FREQUENCY_MAP_POLICY,
        "material_M0": base_material_contract(),
        "contrast": {
            "H": "(1+chi) times every E and G of M0",
            "L": "(1-chi) times every E and G of M0",
            "nu12_fixed": NU12_0,
            "rho_fixed": RHO_0,
            "all_ply_angles_deg": 0.0,
        },
        "geometry": {
            "mu": MU,
            "tau": TAU,
            "beta_deg": BETA_DEG,
            "l": L_REFERENCE,
            "L1": L1,
            "L2": L2,
            "L_total": L_TOTAL,
            "width": WIDTH,
            "thickness": THICKNESS,
            "ply_thickness": PLY_THICKNESS,
            "K": K,
            "outer_clamps": True,
            "joint": "frozen ideal rigid RLB joint",
        },
        "configurations": {
            key: {
                "arm1": {
                    "arm_id": 1,
                    "stack_bottom_to_top": value[0],
                    "material_labels": value[0],
                    "ply_angles_deg": [0.0] * 4,
                    "ply_thicknesses": [PLY_THICKNESS] * 4,
                },
                "arm2": {
                    "arm_id": 2,
                    "stack_bottom_to_top": value[1],
                    "material_labels": value[1],
                    "ply_angles_deg": [0.0] * 4,
                    "ply_thicknesses": [PLY_THICKNESS] * 4,
                },
            }
            for key, value in CONFIGURATION_LAYOUTS.items()
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
        "chi_grid": [float(value) for value in chi_grid()],
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
        "root_search_adapter": {
            "source": "RLB-1 frozen determinant/SVD thresholds and refiners",
            "requested_roots": K_PLOT,
            "guard_roots": 1,
            "required_slots": K_GUARD,
            "source_inventory_length_not_inherited": True,
        },
        "explicit_exclusions": [
            "roots_10_and_above",
            "branch_tracking",
            "MAC",
            "mode_shapes",
            "energy_analysis",
            "Ritz",
            "FEM",
            "smoothing",
            "spectral_sweep_runner",
            "certified_audit",
            "commit",
            "push",
        ],
    }


def contract_hash() -> str:
    payload = json.dumps(
        _json_value(contract_payload()), sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest().upper()


def _solution_record(solution: PointSolution, *, benchmark: bool) -> dict[str, Any]:
    return {
        "configuration_id": solution.configuration_id,
        "chi": solution.chi,
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
                "scaled_sigma_ratio": float(row["scaled_sigma_ratio"]),
                "boundary_null_residual": float(row["boundary_null_residual"]),
                "quality_status": row["quality_status"],
            }
            for row in solution.rows
        ],
    }


def _benchmark_payload(
    benchmark_records: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Build the pre-production ETA from the three declared anchor solves.

    The estimate describes the original production decision, so resume must
    not silently replace its denominator by the smaller number of points that
    happen to remain at a later invocation.
    """

    measured_nonzero = [
        float(item["wall_time_seconds"])
        for item in benchmark_records
        if float(item.get("chi", 0.0)) > 0.0
        and not bool(item.get("reused_existing_without_benchmark_record", False))
    ]
    total_unique_base_solves = 1 + len(CONFIGURATIONS) * (len(chi_grid()) - 1)
    declared_anchor_solves = 3  # shared chi=0 plus two nonzero endpoint anchors
    remaining = total_unique_base_solves - declared_anchor_solves
    max_time = max(measured_nonzero, default=0.0)
    eta = 1.5 * max_time * remaining
    return {
        "schema_version": 1,
        "stage_id": STAGE_ID,
        "anchors": [dict(item) for item in benchmark_records],
        "measured_nonzero_max_seconds": max_time,
        "total_unique_base_solves": total_unique_base_solves,
        "declared_anchor_solves": declared_anchor_solves,
        "remaining_unique_root_points": remaining,
        "eta_formula": (
            "1.5*max(measured_nonzero_point_time)*remaining_root_points"
        ),
        "conservative_eta_seconds": eta,
        "eta_limit_seconds": ETA_LIMIT_SECONDS,
        "production_run_permitted": eta <= ETA_LIMIT_SECONDS,
        "shared_chi0_calculated_once": True,
    }


def run_benchmarks(
    output_dir: Path,
    rows: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    started_at: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run/reuse the three declared production anchors and enforce ETA."""

    target = Path(output_dir)
    benchmark_records: list[dict[str, Any]] = []
    groups = _complete_base_group_index(rows)
    existing_benchmark_path = target / BENCHMARK_FILENAME
    existing_benchmark = (
        json.loads(existing_benchmark_path.read_text(encoding="utf-8"))
        if existing_benchmark_path.is_file()
        else {}
    )
    existing_anchor_index: dict[tuple[str, float], Mapping[str, Any]] = {}
    for item in point_records:
        if bool(item.get("benchmark", False)):
            existing_anchor_index[
                (
                    str(item.get("configuration_id")),
                    round(float(item.get("chi", -1.0)), 10),
                )
            ] = item
    for item in existing_benchmark.get("anchors", []):
        existing_anchor_index[
            (
                str(item.get("configuration_id")),
                round(float(item.get("chi", -1.0)), 10),
            )
        ] = item

    shared_solution: PointSolution | None = None
    if (CONFIG_BOTH_OUTER, 0.0) not in groups:
        try:
            shared = solve_point(CONFIG_BOTH_OUTER, 0.0, force_anchor=True)
        except Exception as exc:
            _write_failure_checkpoint(
                target,
                rows,
                point_records,
                constitutive=constitutive,
                started_at=started_at,
                benchmark_status="RUNNING",
                configuration_id=CONFIG_BOTH_OUTER,
                chi=0.0,
                error=exc,
            )
            raise
        rows = _write_point_transaction(target, rows, shared)
        point_records.append(_solution_record(shared, benchmark=True))
        benchmark_records.append(_solution_record(shared, benchmark=True))
        shared_solution = shared
        _atomic_write_json(
            target / BENCHMARK_FILENAME,
            _benchmark_payload(benchmark_records),
        )
        _atomic_write_json(
            target / CHECKPOINT_FILENAME,
            _checkpoint_payload(
                rows,
                point_records,
                constitutive=constitutive,
                started_at=started_at,
                benchmark_status="RUNNING",
            ),
        )
    else:
        preserved = existing_anchor_index.get((CONFIG_BOTH_OUTER, 0.0))
        if preserved is None:
            raise RuntimeError(
                "The shared chi=0 anchor exists without its benchmark metrics; "
                "ETA cannot be reconstructed safely."
            )
        benchmark_records.append(dict(preserved))

    source_rows = [
        dict(row)
        for row in _complete_base_group_index(rows)[(CONFIG_BOTH_OUTER, 0.0)]
    ]
    source_rows.sort(key=lambda row: int(row["sorted_position"]))
    if shared_solution is None:
        shared_solution = PointSolution(
            configuration_id=CONFIG_BOTH_OUTER,
            chi=0.0,
            rows=tuple(source_rows),
            wall_time_seconds=0.0,
            peak_rss_bytes=0,
            determinant_evaluations=0,
            sigma_evaluations=0,
            search_left_Omega=float(source_rows[0]["locator_interval_left_Omega"]),
            search_right_Omega=float(source_rows[-1]["search_right_Omega"]),
            local_refinements=0,
            solve_mode="SHARED_CHI0_SOURCE_REUSE",
            fallback_used=False,
            unresolved_candidates_below_root9=0,
        )
    for configuration_id in (CONFIG_BOTH_INNER, CONFIG_ANTI_PHASE):
        groups = _complete_base_group_index(rows)
        if (configuration_id, 0.0) in groups:
            continue
        reused_rows: list[dict[str, Any]] = []
        for row in source_rows:
            clone = dict(row)
            clone["configuration_id"] = configuration_id
            clone["row_id"] = str(clone["row_id"]).replace(
                CONFIG_BOTH_OUTER, configuration_id, 1
            )
            clone["shared_chi0_anchor_reused"] = True
            clone["shared_chi0_source_configuration"] = CONFIG_BOTH_OUTER
            reused_rows.append(clone)
        reused = replace(
            shared_solution,
            configuration_id=configuration_id,
            rows=tuple(reused_rows),
            wall_time_seconds=0.0,
            determinant_evaluations=0,
            sigma_evaluations=0,
            solve_mode="SHARED_CHI0_REUSE",
        )
        rows = _write_point_transaction(target, rows, reused)
        point_records.append(_solution_record(reused, benchmark=True))

    # Make the shared-anchor reuse transactions durable before either
    # nonzero benchmark starts.  A restart can recover the measured anchor
    # record from checkpoint point_records even if final benchmark assembly
    # did not run.
    _atomic_write_json(
        target / CHECKPOINT_FILENAME,
        _checkpoint_payload(
            rows,
            point_records,
            constitutive=constitutive,
            started_at=started_at,
            benchmark_status="RUNNING",
        ),
    )

    for configuration_id in (CONFIG_BOTH_OUTER, CONFIG_ANTI_PHASE):
        groups = _complete_base_group_index(rows)
        if (configuration_id, 0.8) in groups:
            preserved = existing_anchor_index.get((configuration_id, 0.8))
            if preserved is None:
                raise RuntimeError(
                    f"The {configuration_id} chi=0.8 anchor exists without "
                    "benchmark metrics; ETA cannot be reconstructed safely."
                )
            benchmark_records.append(dict(preserved))
            continue
        try:
            solution = solve_point(configuration_id, 0.8, force_anchor=True)
        except Exception as exc:
            _write_failure_checkpoint(
                target,
                rows,
                point_records,
                constitutive=constitutive,
                started_at=started_at,
                benchmark_status="RUNNING",
                configuration_id=configuration_id,
                chi=0.8,
                error=exc,
            )
            raise
        rows = _write_point_transaction(target, rows, solution)
        record = _solution_record(solution, benchmark=True)
        point_records.append(record)
        benchmark_records.append(record)
        _atomic_write_json(
            target / BENCHMARK_FILENAME,
            _benchmark_payload(benchmark_records),
        )
        _atomic_write_json(
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
    eta = float(benchmark["conservative_eta_seconds"])
    _atomic_write_json(target / BENCHMARK_FILENAME, benchmark)
    _atomic_write_json(
        target / CHECKPOINT_FILENAME,
        _checkpoint_payload(
            rows,
            point_records,
            constitutive=constitutive,
            started_at=started_at,
            benchmark_status=("PASS" if eta <= ETA_LIMIT_SECONDS else "STOPPED_BY_ETA_GATE"),
        ),
    )
    return rows, benchmark


def _existing_rows(output_dir: Path) -> list[dict[str, Any]]:
    path = Path(output_dir) / SPECTRUM_FILENAME
    return [] if not path.is_file() else [dict(row) for row in _read_csv(path)]


def complete_missing_points(
    output_dir: Path,
    rows: list[dict[str, Any]],
    point_records: list[dict[str, Any]],
    constitutive: Mapping[str, Any],
    started_at: str,
) -> list[dict[str, Any]]:
    """Walk each configuration sequentially and never recompute a full group."""

    target = Path(output_dir)
    for configuration_id in CONFIGURATIONS:
        for chi_value in chi_grid():
            chi = float(chi_value)
            groups = _complete_base_group_index(rows)
            if (configuration_id, round(chi, 10)) in groups:
                continue
            earlier = sorted(
                value
                for (config, value) in groups
                if config == configuration_id and value < round(chi, 10)
            )
            previous = None
            second = None
            if earlier:
                last = earlier[-1]
                previous = (last, _rows_for_roots(rows, configuration_id, last))
            if len(earlier) >= 2:
                older = earlier[-2]
                second = (older, _rows_for_roots(rows, configuration_id, older))
            try:
                solution = solve_point(
                    configuration_id,
                    chi,
                    previous=previous,
                    second_previous=second,
                    force_anchor=previous is None,
                )
            except Exception as exc:
                _write_failure_checkpoint(
                    target,
                    rows,
                    point_records,
                    constitutive=constitutive,
                    started_at=started_at,
                    benchmark_status="PASS",
                    configuration_id=configuration_id,
                    chi=chi,
                    error=exc,
                )
                raise
            rows = _write_point_transaction(target, rows, solution)
            point_records.append(_solution_record(solution, benchmark=False))
            _atomic_write_json(
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


def centred_secant_residual(left: float, center: float, right: float) -> float:
    predictor = 0.5 * (float(left) + float(right))
    return abs(float(center) - predictor) / max(
        abs(float(center)), abs(predictor), np.finfo(float).tiny
    )


def neighbour_audit_rows(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Apply the declared median+8*MAD anomaly trigger to plotted positions."""

    audit_rows: list[dict[str, Any]] = []
    gap_flags: set[tuple[str, float, int]] = set()
    base_grid = [round(float(value), 10) for value in chi_grid()]
    for configuration_id in CONFIGURATIONS:
        spectra = {
            chi: np.sqrt(_rows_for_roots(rows, configuration_id, chi)[:K_PLOT])
            for chi in base_grid
        }
        for lower_position in range(1, K_PLOT):
            gaps = np.asarray(
                [
                    spectra[chi][lower_position]
                    - spectra[chi][lower_position - 1]
                    for chi in base_grid
                ],
                dtype=float,
            )
            residuals = np.asarray(
                [
                    centred_secant_residual(
                        gaps[index - 1], gaps[index], gaps[index + 1]
                    )
                    for index in range(1, len(base_grid) - 1)
                ],
                dtype=float,
            )
            median = float(np.median(residuals))
            mad = float(np.median(np.abs(residuals - median)))
            threshold = median + NEIGHBOUR_MAD_MULTIPLIER * mad
            for offset, index in enumerate(range(1, len(base_grid) - 1)):
                if (
                    float(residuals[offset]) > threshold
                    and float(residuals[offset]) > NEIGHBOUR_ABSOLUTE_TRIGGER
                ):
                    for position in (lower_position, lower_position + 1):
                        gap_flags.add(
                            (configuration_id, base_grid[index], position)
                        )
    for configuration_id in CONFIGURATIONS:
        for position in range(1, K_PLOT + 1):
            values: dict[float, float] = {}
            for row in rows:
                if (
                    str(row["configuration_id"]) == configuration_id
                    and int(row["sorted_position"]) == position
                    and _as_bool(row.get("is_canonical_plot_source", True))
                ):
                    values[round(float(row["chi"]), 10)] = float(row["Lambda"])
            grid = base_grid
            residuals = [
                centred_secant_residual(
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
                chi = grid[index]
                residual = residuals[offset]
                statistical_flag = bool(
                    residual > robust_threshold
                    and residual > NEIGHBOUR_ABSOLUTE_TRIGGER
                )
                center_group = [
                    row
                    for row in rows
                    if str(row["configuration_id"]) == configuration_id
                    and str(row["grid_kind"]) == "BASE"
                    and round(float(row["chi"]), 10) == chi
                ]
                root_row = next(
                    row
                    for row in center_group
                    if int(row["sorted_position"]) == position
                )
                root_count_warning = len(center_group) != K_GUARD
                ordering_warning = any(
                    float(left["Omega"]) >= float(right["Omega"])
                    for left, right in zip(
                        sorted(center_group, key=lambda row: int(row["sorted_position"]))[:-1],
                        sorted(center_group, key=lambda row: int(row["sorted_position"]))[1:],
                        strict=True,
                    )
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
                gap_warning = (
                    configuration_id, chi, position
                ) in gap_flags
                flagged = bool(
                    statistical_flag
                    or root_count_warning
                    or ordering_warning
                    or unresolved_warning
                    or bad_residual_warning
                    or gap_warning
                )
                audit_rows.append(
                    {
                        "configuration_id": configuration_id,
                        "sorted_position": position,
                        "chi_left": grid[index - 1],
                        "chi": chi,
                        "chi_right": grid[index + 1],
                        "Lambda_left": values[grid[index - 1]],
                        "Lambda_center": values[chi],
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
                        "local_chi_values": [],
                        "smoothing_applied": False,
                    }
                )
    return audit_rows


def flagged_repair_points(audit_rows: Sequence[Mapping[str, Any]]) -> list[tuple[str, float]]:
    return sorted(
        {
            (str(row["configuration_id"]), round(float(row["chi"]), 10))
            for row in audit_rows
            if bool(row["flagged"])
        },
        key=lambda item: (CONFIGURATIONS.index(item[0]), item[1]),
    )


def apply_local_repairs(
    rows: list[dict[str, Any]], audit_rows: list[dict[str, Any]]
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    """Re-solve only flagged base points; never smooth or interpolate."""

    repair_records: list[dict[str, Any]] = []
    for repair_index, (configuration_id, chi) in enumerate(
        flagged_repair_points(audit_rows), start=1
    ):
        repair_id = f"repair_{repair_index:04d}"
        affected_positions = sorted(
            {
                int(item["sorted_position"])
                for item in audit_rows
                if str(item["configuration_id"]) == configuration_id
                and round(float(item["chi"]), 10) == round(chi, 10)
                and bool(item["flagged"])
            }
        )
        existing_repair = [
            row
            for row in rows
            if str(row["configuration_id"]) == configuration_id
            and round(float(row["chi"]), 10) == round(chi, 10)
            and str(row["grid_kind"]) == "LOCAL_REFINEMENT"
            and str(row.get("repair_id", "")) == repair_id
        ]
        if sorted(int(row["sorted_position"]) for row in existing_repair) == list(
            range(1, K_GUARD + 1)
        ) and all(
            str(row["point_status"])
            in {
                "REPRODUCED_AFTER_LOCAL_REPAIR",
                "LOCATOR_CORRECTED_AFTER_LOCAL_REPAIR",
            }
            for row in existing_repair
        ):
            original = np.asarray(
                [
                    float(row["Omega"])
                    for row in sorted(
                        (
                            item
                            for item in rows
                            if str(item["configuration_id"]) == configuration_id
                            and round(float(item["chi"]), 10) == round(chi, 10)
                            and str(item["grid_kind"]) == "BASE"
                        ),
                        key=lambda item: int(item["sorted_position"]),
                    )
                ],
                dtype=float,
            )
            refined = np.asarray(
                [
                    float(row["Omega"])
                    for row in sorted(
                        existing_repair,
                        key=lambda item: int(item["sorted_position"]),
                    )
                ],
                dtype=float,
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
            status = str(existing_repair[0]["point_status"])
            for audit in audit_rows:
                if (
                    str(audit["configuration_id"]) == configuration_id
                    and round(float(audit["chi"]), 10) == round(chi, 10)
                    and bool(audit["flagged"])
                ):
                    audit["repair_id"] = repair_id
                    audit["repair_status"] = status
            repair_records.append(
                {
                    "repair_id": repair_id,
                    "configuration_id": configuration_id,
                    "chi": chi,
                    "status": status,
                    "affected_positions": affected_positions,
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
        previous = (
            chi - CHI_STEP,
            _rows_for_roots(rows, configuration_id, chi - CHI_STEP),
        )
        second = (
            chi - 2.0 * CHI_STEP,
            _rows_for_roots(rows, configuration_id, chi - 2.0 * CHI_STEP),
        ) if chi >= 2.0 * CHI_STEP else None
        try:
            repaired = solve_point(
                configuration_id,
                chi,
                previous=previous,
                second_previous=second,
                dense_local=True,
                dense_positions=affected_positions,
                grid_kind="LOCAL_REFINEMENT",
                repair_id=repair_id,
            )
        except (RuntimeError, ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            gap_rows: list[dict[str, Any]] = []
            for row in rows:
                if (
                    str(row["configuration_id"]) == configuration_id
                    and round(float(row["chi"]), 10) == round(chi, 10)
                    and str(row["grid_kind"]) == "BASE"
                    and int(row["sorted_position"]) in affected_positions
                ):
                    row["is_canonical_plot_source"] = False
                    row["point_status"] = "UNRESOLVED_AFTER_LOCAL_REPAIR"
                    gap = dict(row)
                    gap["row_id"] = (
                        f"{configuration_id}__{float(chi):.6f}__"
                        f"LOCAL_REFINEMENT__p{int(row['sorted_position']):02d}__"
                        f"{repair_id}_gap"
                    )
                    gap["grid_kind"] = "LOCAL_REFINEMENT"
                    gap["Lambda"] = math.nan
                    gap["quality_status"] = "UNRESOLVED_AFTER_LOCAL_REPAIR"
                    gap["is_canonical_plot_source"] = True
                    gap["supersedes_row_id"] = row["row_id"]
                    gap["repair_id"] = repair_id
                    gap_rows.append(gap)
            rows.extend(gap_rows)
            for audit in audit_rows:
                if (
                    str(audit["configuration_id"]) == configuration_id
                    and round(float(audit["chi"]), 10) == round(chi, 10)
                    and bool(audit["flagged"])
                ):
                    audit["repair_id"] = repair_id
                    audit["repair_status"] = "UNRESOLVED_AFTER_LOCAL_REPAIR"
            repair_records.append(
                {
                    "repair_id": repair_id,
                    "configuration_id": configuration_id,
                    "chi": chi,
                    "status": "UNRESOLVED",
                    "affected_positions": affected_positions,
                    "error": f"{type(exc).__name__}: {exc}",
                    "wall_time_seconds": 0.0,
                    "determinant_evaluations": 0,
                    "sigma_evaluations": 0,
                    "smoothing_applied": False,
                    "predictor_used_as_final": False,
                }
            )
            continue
        original = _rows_for_roots(rows, configuration_id, chi)
        refined = np.asarray([float(row["Omega"]) for row in repaired.rows])
        relative = float(
            np.max(
                np.abs(original - refined)
                / np.maximum.reduce(
                    (np.abs(original), np.abs(refined), np.full(9, np.finfo(float).tiny))
                )
            )
        )
        status = "REPRODUCED_AFTER_LOCAL_REPAIR"
        if relative > 1.0e-8:
            status = "LOCATOR_CORRECTED_AFTER_LOCAL_REPAIR"
        for row in rows:
            if (
                str(row["configuration_id"]) == configuration_id
                and round(float(row["chi"]), 10) == round(chi, 10)
                and str(row["grid_kind"]) == "BASE"
            ):
                row["is_canonical_plot_source"] = False
                row["point_status"] = status
        repair_rows: list[dict[str, Any]] = []
        for row in repaired.rows:
            clone = dict(row)
            clone["grid_kind"] = "LOCAL_REFINEMENT"
            clone["repair_id"] = repair_id
            clone["is_canonical_plot_source"] = True
            clone["supersedes_row_id"] = next(
                str(item["row_id"])
                for item in rows
                if str(item["configuration_id"]) == configuration_id
                and round(float(item["chi"]), 10) == round(chi, 10)
                and int(item["sorted_position"]) == int(clone["sorted_position"])
                and str(item["grid_kind"]) == "BASE"
            )
            clone["point_status"] = status
            repair_rows.append(clone)
        rows.extend(repair_rows)
        for audit in audit_rows:
            if (
                str(audit["configuration_id"]) == configuration_id
                and round(float(audit["chi"]), 10) == round(chi, 10)
                and bool(audit["flagged"])
            ):
                audit["repair_id"] = repair_id
                audit["repair_status"] = status
        repair_records.append(
            {
                "repair_id": repair_id,
                "configuration_id": configuration_id,
                "chi": chi,
                "status": status,
                "affected_positions": affected_positions,
                "maximum_relative_Omega_change": relative,
                "wall_time_seconds": repaired.wall_time_seconds,
                "determinant_evaluations": repaired.determinant_evaluations,
                "sigma_evaluations": repaired.sigma_evaluations,
                "smoothing_applied": False,
                "predictor_used_as_final": False,
            }
        )
    return rows, audit_rows, repair_records


def canonical_plot_rows(rows: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    return [
        row
        for row in rows
        if int(row["sorted_position"]) <= K_PLOT
        and str(row.get("is_canonical_plot_source", "True")).lower()
        in {"true", "1", "yes"}
    ]


def audit_plot_rows(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Require one canonical plotted value for every declared map key."""

    plot_rows = canonical_plot_rows(rows)
    expected = {
        (configuration, round(float(chi), 10), position)
        for configuration in CONFIGURATIONS
        for chi in chi_grid()
        for position in range(1, K_PLOT + 1)
    }
    counts: dict[tuple[str, float, int], int] = {}
    invalid_values: list[str] = []
    for row in plot_rows:
        key = (
            str(row["configuration_id"]),
            round(float(row["chi"]), 10),
            int(row["sorted_position"]),
        )
        counts[key] = counts.get(key, 0) + 1
        value = float(row["Lambda"])
        unresolved_gap = (
            math.isnan(value)
            and str(row.get("point_status", ""))
            == "UNRESOLVED_AFTER_LOCAL_REPAIR"
        )
        if not ((math.isfinite(value) and value > 0.0) or unresolved_gap):
            invalid_values.append(f"{key[0]}:{key[1]:.2f}:p{key[2]}")
    duplicate = sorted(key for key, count in counts.items() if count != 1)
    missing = sorted(expected - set(counts))
    extra = sorted(set(counts) - expected)
    return {
        "status": (
            "PASS"
            if not duplicate and not missing and not extra and not invalid_values
            else "FAIL"
        ),
        "row_count": len(plot_rows),
        "missing_keys": [f"{key[0]}:{key[1]:.2f}:p{key[2]}" for key in missing],
        "duplicate_keys": [
            f"{key[0]}:{key[1]:.2f}:p{key[2]}" for key in duplicate
        ],
        "extra_keys": [f"{key[0]}:{key[1]:.2f}:p{key[2]}" for key in extra],
        "invalid_values": invalid_values,
    }


def create_plot_from_csv(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    """Render only from CSV; this function imports no root-calculation path."""

    started = time.perf_counter()
    rows = _read_csv(Path(output_dir) / SPECTRUM_FILENAME)
    spectrum_audit = audit_spectrum_rows(rows)
    if spectrum_audit["status"] != "PASS":
        raise RuntimeError(
            f"plot_only rejected incomplete root-guard data: {spectrum_audit}"
        )
    plot_audit = audit_plot_rows(rows)
    if plot_audit["status"] != "PASS":
        raise RuntimeError(f"plot_only rejected incomplete CSV data: {plot_audit}")
    plot_rows = canonical_plot_rows(rows)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = plt.get_cmap("tab10").colors[:K_PLOT]
    figure, axes = plt.subplots(1, 3, figsize=(14.5, 4.8), sharex=True, sharey=True)
    titles = (
        "(a) BOTH_OUTER_STIFF",
        "(b) BOTH_INNER_STIFF",
        "(c) ANTI_PHASE",
    )
    for axis, configuration_id, title in zip(
        axes, CONFIGURATIONS, titles, strict=True
    ):
        for position in range(1, K_PLOT + 1):
            selected = [
                row
                for row in plot_rows
                if str(row["configuration_id"]) == configuration_id
                and int(row["sorted_position"]) == position
            ]
            selected.sort(key=lambda row: float(row["chi"]))
            axis.plot(
                [float(row["chi"]) for row in selected],
                [float(row["Lambda"]) for row in selected],
                color=colors[position - 1],
                linestyle="-",
                linewidth=1.25,
                label=f"k={position}",
            )
        axis.set_title(title)
        axis.set_xlabel(r"$\chi$")
        axis.grid(True, alpha=0.22, linewidth=0.5)
    axes[0].set_ylabel(r"$\Lambda$")
    figure.legend(
        handles=axes[0].lines,
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
        "panel_count": 3,
        "lines_per_panel": K_PLOT,
        "plotted_positions": list(range(1, K_PLOT + 1)),
        "root9_plotted": False,
        "root_calculation_count": 0,
        "spectrum_data_audit": spectrum_audit,
        "plot_data_audit": plot_audit,
    }


def arm_swap_checks(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Check the unplotted ANTI_PHASE arm permutation at chi=.4 and .8."""

    checks: list[dict[str, Any]] = []
    for chi in (0.4, 0.8):
        reference = _rows_for_roots(rows, CONFIG_ANTI_PHASE, chi)
        # Swapping ANTI_PHASE arms is represented by the same frozen provider
        # with the explicit reversed layout order, not a fourth configuration.
        original = CONFIGURATION_LAYOUTS[CONFIG_ANTI_PHASE]
        CONFIGURATION_LAYOUTS[CONFIG_ANTI_PHASE] = (original[1], original[0])
        try:
            solution = solve_point(
                CONFIG_ANTI_PHASE,
                chi,
                previous=(chi, reference),
                force_anchor=False,
                guard_locator_right_width=0.2,
            )
        finally:
            CONFIGURATION_LAYOUTS[CONFIG_ANTI_PHASE] = original
        swapped = np.asarray([float(row["Omega"]) for row in solution.rows])
        relatives = np.abs(reference - swapped) / np.maximum.reduce(
            (np.abs(reference), np.abs(swapped), np.full(9, np.finfo(float).tiny))
        )
        checks.append(
            {
                "chi": chi,
                "maximum_relative_Omega": float(np.max(relatives)),
                "tolerance": ARM_SWAP_RELATIVE_TOLERANCE,
                "status": (
                    "PASS"
                    if float(np.max(relatives)) <= ARM_SWAP_RELATIVE_TOLERANCE
                    else "FAIL"
                ),
                "root_count": len(swapped),
                "roots_above_9_computed": False,
                "wall_time_seconds": solution.wall_time_seconds,
                "determinant_evaluations": solution.determinant_evaluations,
                "sigma_evaluations": solution.sigma_evaluations,
            }
        )
    return checks


def _output_hashes(output_dir: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for path in sorted(Path(output_dir).iterdir()):
        if path.is_file() and path.name != MANIFEST_FILENAME:
            result[path.name] = _sha256(path)
    return result


def _minimum_adjacent_gaps(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for configuration_id in CONFIGURATIONS:
        best = (math.inf, math.nan, -1)
        for chi in chi_grid():
            roots = _rows_for_roots(rows, configuration_id, float(chi))
            gaps = np.diff(np.sqrt(roots[:K_PLOT]))
            index = int(np.argmin(gaps))
            if float(gaps[index]) < best[0]:
                best = (float(gaps[index]), float(chi), index + 1)
        records.append(
            {
                "configuration_id": configuration_id,
                "minimum_adjacent_Lambda_gap": best[0],
                "chi": best[1],
                "between_sorted_positions": [best[2], best[2] + 1],
                "interpretation": "candidate interval only; no branch or shape claim",
            }
        )
    return records


def _spectrum_quality_summary(
    rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
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
    arm_swap: Sequence[Mapping[str, Any]],
    plot: Mapping[str, Any],
    *,
    finalization_process_seconds: float,
) -> dict[str, Any]:
    base_seconds = sum(float(item.get("wall_time_seconds", 0.0)) for item in point_records)
    repair_seconds = sum(
        float(item.get("wall_time_seconds", 0.0)) for item in repair_records
    )
    arm_swap_seconds = sum(
        float(item.get("wall_time_seconds", 0.0)) for item in arm_swap
    )
    plot_seconds = float(plot.get("wall_time_seconds", 0.0))
    determinant_evaluations = (
        sum(int(item.get("determinant_evaluations", 0)) for item in point_records)
        + sum(int(item.get("determinant_evaluations", 0)) for item in repair_records)
        + sum(int(item.get("determinant_evaluations", 0)) for item in arm_swap)
    )
    sigma_evaluations = (
        sum(int(item.get("sigma_evaluations", 0)) for item in point_records)
        + sum(int(item.get("sigma_evaluations", 0)) for item in repair_records)
        + sum(int(item.get("sigma_evaluations", 0)) for item in arm_swap)
    )
    recorded_peak = max(
        (int(item.get("peak_rss_bytes", 0)) for item in point_records),
        default=0,
    )
    base_root_solve_count = sum(
        int(item.get("determinant_evaluations", 0)) > 0 for item in point_records
    )
    local_repair_solve_count = sum(
        int(item.get("determinant_evaluations", 0)) > 0
        for item in repair_records
    )
    reused_local_repair_count = len(repair_records) - local_repair_solve_count
    return {
        "base_point_wall_time_sum_seconds": base_seconds,
        "local_repair_wall_time_sum_seconds": repair_seconds,
        "arm_swap_wall_time_sum_seconds": arm_swap_seconds,
        "plot_only_seconds": plot_seconds,
        "total_measured_workflow_seconds": (
            base_seconds + repair_seconds + arm_swap_seconds + plot_seconds
        ),
        "finalization_process_seconds": float(finalization_process_seconds),
        "peak_rss_bytes": max(_peak_rss_bytes(), recorded_peak),
        "determinant_evaluations": determinant_evaluations,
        "sigma_evaluations": sigma_evaluations,
        "base_root_solve_count": base_root_solve_count,
        "local_repair_solve_count": local_repair_solve_count,
        "reused_local_repair_count": reused_local_repair_count,
        "arm_swap_solve_count": len(arm_swap),
        "total_root_solve_count": (
            base_root_solve_count + local_repair_solve_count + len(arm_swap)
        ),
        "parallel_spectral_workers": 0,
    }


def _report_text(manifest: Mapping[str, Any]) -> str:
    counts = manifest["counts"]
    gate = manifest["constitutive_gate"]
    benchmark = manifest["benchmark"]
    repairs = manifest["neighbour_audit"]
    runtime = manifest["runtime"]
    quality = manifest["root_quality_summary"]
    status = manifest["scientific_status"]
    gap_lines = "\n".join(
        "- `{configuration_id}`: $\\Delta\\Lambda_{{min}}={gap:.6g}$ при "
        "$\\chi={chi:.2f}$ между позициями {left} и {right}.".format(
            configuration_id=item["configuration_id"],
            gap=float(item["minimum_adjacent_Lambda_gap"]),
            chi=float(item["chi"]),
            left=int(item["between_sorted_positions"][0]),
            right=int(item["between_sorted_positions"][1]),
        )
        for item in manifest["minimum_adjacent_sorted_gaps"]
    )
    return f"""# RLB-2E: контраст расположения жёстких и мягких слоёв

## Цель и область применимости

Этап строит конечную карту первых восьми независимо отсортированных частот
двух жёстко сопряжённых симметрично слоистых стержней при изменении
контраста жёсткости $\\chi$. Идентичность модальных ветвей, формы, MAC и
распределения энергии не определяются.

## Фиксированный контракт

Использованы $\\mu=\\tau=0$, $\\beta=30^\\circ$, $L_1=L_2=1$,
$b=0.20$, $h=0.05$ и $K=5/6$. Каждый стержень содержит четыре равных
слоя толщины 0.0125; угол материала во всех слоях равен $0^\\circ$.

Базовый материал имеет $E_1=1.1$, $E_2=0.9$, $\\nu_{{12}}=0.3$,
$G_{{12}}=G_{{13}}=G_{{23}}=1/2.6$ и $\\rho=1$. Для материала H все
модули упругости умножены на $1+\\chi$, для материала L — на
$1-\\chi$; коэффициент Пуассона и плотность не менялись.

Рассмотрены конфигурации `BOTH_OUTER_STIFF`, `BOTH_INNER_STIFF` и
`ANTI_PHASE`. Стопки задавались снизу вверх как H/L/L/H или L/H/H/L.

## Constitutive gate

Статус: **{gate['status']}**. Для четырёх равных слоёв подтверждены

$$D_{{outer}}/D_0=1+3\\chi/4,\\qquad
D_{{inner}}/D_0=1-3\\chi/4,$$

а также $D_{{outer}}+D_{{inner}}=2D_0$. Перестановка материалов не
изменила $A$, поперечную сдвиговую жёсткость, $I_0$, $I_2$, $A_{{beam}}$,
$S_{{beam}}$, $m$ и $J$; $B$ и $I_1$ остались нулевыми в пределах
объявленных допусков. Максимальный residual полной D-формулы равен
{gate['maximum_residuals']['D_matrix_formula_relative']:.3e}. При
$\\chi=0.8$ получены $D_{{outer}}/D_0=1.6$ и
$D_{{inner}}/D_0=0.4$.

## Frequency-map policy и численная сборка

Локальный instance — `frequency-map-v1`, режим `fast_plot`, семантика
`sorted_positions`, сетка $\\chi=0.00,0.02,\\ldots,0.80$.

Позиции 1–8 показаны на рисунке, root 9 используется только как проверка
полноты.

Прогноз соседними точками определял лишь локальные интервалы. Каждая итоговая
частота вычислена из замороженной coupled RLB matrix с существующими
determinant/SVD detector и refiner.

Локальный adapter сохранял ровно восемь plotted roots и root 9. Принятые
roots 10 и выше отсутствуют.

Нормировка сохранена от RLB-2D:

$$\\Omega=\\omega l^2\\sqrt{{\\rho_0A_0/(E_0I_{{y0}})}},\\qquad
\\Lambda=\\sqrt{{\\Omega}}.$$

## Результаты и диагностика

Получено {counts['base_configuration_points']}/123 base-точек и
{counts['base_rows']}/1107 base-строк. Полных root-9 guards:
{counts['root9_guards']}/123. При $\\chi=0$ один расчёт повторно использован
для трёх конфигураций. После локального контроля сохранено
{counts['local_refinement_rows']} строк `LOCAL_REFINEMENT`.

Production benchmark дал консервативную ETA
{benchmark['conservative_eta_seconds']:.1f} s для исходных
{benchmark['remaining_unique_root_points']} расчётов после трёх anchors при
лимите 2700 s.

Сумма измеренных времён расчёта, локального контроля, проверки перестановки
плеч и построения рисунка равна
{runtime['total_measured_workflow_seconds']:.1f} s; peak RSS составил
{runtime['peak_rss_bytes'] / 2**20:.1f} MiB.

Выполнено
{runtime['determinant_evaluations']} determinant и
{runtime['sigma_evaluations']} SVD/sigma evaluations.

Neighbour audit отметил {repairs['flagged_point_count']} point(s); выполнено
{repairs['repair_count']} локальных repair(s), unresolved point(s):
{repairs['unresolved_point_count']}. Максимальные base residuals равны
{quality['maximum_base_scaled_sigma_ratio']:.3e} для
$\\sigma_{{min}}/\\sigma_{{max}}$ и
{quality['maximum_base_boundary_null_residual']:.3e} для boundary null
residual; минимальный запас root 9 до правой границы равен
{quality['minimum_root9_right_margin_Omega']:.3g} по $\\Omega$.

Проверка перестановки плеч `ANTI_PHASE` при $\\chi=0.4$ и 0.8 имеет
статус **{manifest['arm_swap_status']}**. Минимальные соседние интервалы:

{gap_lines}

Эти интервалы являются только кандидатами для возможного будущего анализа
форм и не доказывают пересечение, veering или сохранение идентичности моды.

## Статус и ограничения

**RLB-2E: {status}.** Вывод относится только к указанной конечной сетке,
фиксированной геометрии, четырём 0°-слоям и трём объявленным раскладкам.
Branch tracking, MAC, формы, energy analysis, Ritz, FEM, smoothing и
certified audit не выполнялись. Производственные физические модули не
изменялись.
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
    audit_rows = neighbour_audit_rows(rows)
    rows, audit_rows, repair_records = apply_local_repairs(rows, audit_rows)
    _atomic_write_csv(target / SPECTRUM_FILENAME, rows, SPECTRUM_FIELDS)
    _atomic_write_csv(target / AUDIT_FILENAME, audit_rows)
    spectrum_audit = audit_spectrum_rows(rows)
    if spectrum_audit["status"] != "PASS":
        raise RuntimeError(f"Final spectrum audit failed: {spectrum_audit}")
    arm_swap = arm_swap_checks(rows)
    if any(item["status"] != "PASS" for item in arm_swap):
        raise RuntimeError(f"ANTI_PHASE arm-swap check failed: {arm_swap}")
    plot = create_plot_from_csv(target)
    unresolved = sum(
        1 for record in repair_records if record["status"] == "UNRESOLVED"
    )
    scientific_status = "PASS" if unresolved == 0 else "PARTIAL_PASS"
    unique_base_solves = sum(
        int(record.get("determinant_evaluations", 0)) > 0
        for record in point_records
    )
    counts = {
        "base_configuration_points": spectrum_audit["base_group_count"],
        "base_rows": spectrum_audit["base_row_count"],
        "local_refinement_rows": sum(
            str(row["grid_kind"]) == "LOCAL_REFINEMENT" for row in rows
        ),
        "root9_guards": sum(
            str(row["grid_kind"]) == "BASE"
            and int(row["sorted_position"]) == K_GUARD
            for row in rows
        ),
        "unique_physics_base_solves": unique_base_solves,
        "computed_base_configuration_points": unique_base_solves,
        "reused_base_configuration_points": 2,
        "shared_chi0_reused_configuration_groups": 2,
        "locally_repaired_points": len(repair_records),
        "unresolved_points": unresolved,
    }
    runtime = _runtime_summary(
        point_records,
        repair_records,
        arm_swap,
        plot,
        finalization_process_seconds=time.perf_counter() - total_started,
    )
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "scientific_status": scientific_status,
        "completed_at_utc": _utc_now(),
        "git": _git_state(),
        "contract": contract_payload(),
        "contract_sha256": contract_hash(),
        "production_physics_hashes": {
            path: _sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS
        },
        "analysis_script_sha256": _sha256(Path(__file__)),
        "constitutive_gate": dict(constitutive),
        "counts": counts,
        "spectrum_audit": spectrum_audit,
        "root_quality_summary": _spectrum_quality_summary(rows),
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
        "arm_swap_checks": arm_swap,
        "arm_swap_status": "PASS",
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
    checkpoint["completed_at_utc"] = _utc_now()
    checkpoint["root_calculation_count"] = runtime["total_root_solve_count"]
    checkpoint["local_repair_count"] = len(repair_records)
    checkpoint["arm_swap_check_count"] = len(arm_swap)
    checkpoint["terminal_unresolved_points"] = [
        {
            "configuration_id": record["configuration_id"],
            "chi": record["chi"],
            "repair_id": record["repair_id"],
        }
        for record in repair_records
        if record["status"] == "UNRESOLVED"
    ]
    _atomic_write_json(target / CHECKPOINT_FILENAME, checkpoint)
    report = _report_text(manifest)
    _atomic_write_text(target / REPORT_FILENAME, report)
    manifest["output_hashes"] = _output_hashes(target)
    _atomic_write_json(target / MANIFEST_FILENAME, manifest)
    return manifest


def refresh_completed_outputs(
    output_dir: Path,
    rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Refresh final metadata without determinant, SVD, repair, or plotting work."""

    target = Path(output_dir)
    manifest_path = target / MANIFEST_FILENAME
    if not manifest_path.is_file():
        raise RuntimeError("A completed spectrum requires an existing run manifest.")
    spectrum_audit = audit_spectrum_rows(rows)
    if spectrum_audit["status"] != "PASS":
        raise RuntimeError(f"Completed-output refresh rejected: {spectrum_audit}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    root_contract = manifest.get("root_contract", {})
    compatible = bool(
        manifest.get("stage_id") == STAGE_ID
        and manifest.get("algorithm_version") == ALGORITHM_VERSION
        and manifest.get("contract_sha256") == contract_hash()
        and root_contract.get("search_policy_requested_roots") == K_PLOT
        and root_contract.get("search_policy_guard_roots") == 1
        and root_contract.get("guard_position") == K_GUARD
        and root_contract.get("root9_plotted") is False
    )
    if not compatible:
        raise RuntimeError(
            "Completed-output refresh refused incompatible legacy provenance."
        )
    expected_physics = {
        path: _sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS
    }
    if manifest.get("production_physics_hashes") != expected_physics:
        raise RuntimeError(
            "Completed-output refresh refused changed production physics."
        )
    recorded_hashes = manifest.get("output_hashes", {})
    if not recorded_hashes or any(
        not (target / name).is_file()
        or _sha256(target / name) != expected
        for name, expected in recorded_hashes.items()
    ):
        raise RuntimeError(
            "Completed-output refresh refused modified or incomplete evidence."
        )
    current_git = _git_state()
    current_script_hash = _sha256(Path(__file__))
    if (
        manifest.get("analysis_script_sha256") == current_script_hash
        and manifest.get("git") == current_git
    ):
        # Fully current completed outputs are byte-stable under resume.
        return manifest
    benchmark_path = target / BENCHMARK_FILENAME
    benchmark_existing = json.loads(benchmark_path.read_text(encoding="utf-8"))
    benchmark = _benchmark_payload(benchmark_existing["anchors"])
    _atomic_write_json(benchmark_path, benchmark)

    point_records = [dict(item) for item in manifest.get("point_records", [])]
    repairs = [
        dict(item)
        for item in manifest.get("neighbour_audit", {}).get("repair_records", [])
    ]
    arm_swap = [dict(item) for item in manifest.get("arm_swap_checks", [])]
    plot = dict(manifest["plot"])
    prior_runtime = dict(manifest.get("runtime", {}))
    runtime = _runtime_summary(
        point_records,
        repairs,
        arm_swap,
        plot,
        finalization_process_seconds=float(
            prior_runtime.get("finalization_process_seconds", 0.0)
        ),
    )
    runtime["peak_rss_bytes"] = max(
        int(runtime["peak_rss_bytes"]),
        int(prior_runtime.get("peak_rss_bytes", 0)),
    )
    manifest.update(
        {
            "algorithm_version": ALGORITHM_VERSION,
            "metadata_refreshed_at_utc": _utc_now(),
            "git": current_git,
            "contract": contract_payload(),
            "contract_sha256": contract_hash(),
            "production_physics_hashes": expected_physics,
            "analysis_script_sha256": current_script_hash,
            "benchmark": benchmark,
            "spectrum_audit": spectrum_audit,
            "root_quality_summary": _spectrum_quality_summary(rows),
            "runtime": runtime,
        }
    )
    counts = dict(manifest["counts"])
    base_solves = int(runtime["base_root_solve_count"])
    counts.update(
        {
            "unique_physics_base_solves": base_solves,
            "computed_base_configuration_points": base_solves,
            "reused_base_configuration_points": 2,
            "locally_repaired_points": len(repairs),
            "unresolved_points": int(
                manifest["neighbour_audit"]["unresolved_point_count"]
            ),
        }
    )
    manifest["counts"] = counts

    checkpoint_path = target / CHECKPOINT_FILENAME
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    checkpoint.update(
        {
            "algorithm_version": ALGORITHM_VERSION,
            "contract_sha256": contract_hash(),
            "updated_at_utc": _utc_now(),
            "benchmark_status": "PASS",
            "scientific_status": manifest["scientific_status"],
            "root_calculation_count": runtime["total_root_solve_count"],
            "local_repair_count": len(repairs),
            "arm_swap_check_count": len(arm_swap),
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
    )
    _atomic_write_json(checkpoint_path, checkpoint)
    _atomic_write_text(target / REPORT_FILENAME, _report_text(manifest))
    manifest["output_hashes"] = _output_hashes(target)
    _atomic_write_json(manifest_path, manifest)
    return manifest


def prepare_constitutive_outputs(output_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)
    gate = constitutive_gate()
    rows = section_property_rows()
    _atomic_write_csv(target / SECTION_FILENAME, rows)
    if gate["status"] != "PASS":
        raise RuntimeError(f"Constitutive gate failed before roots: {gate}")
    return gate, rows


def run_workflow(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    *,
    missing_only: bool = False,
) -> dict[str, Any]:
    """Calculate missing transactions, then audit and render the map."""

    del missing_only  # both modes preserve and calculate only incomplete groups
    total_started = time.perf_counter()
    started_at = _utc_now()
    target = Path(output_dir)
    rows = _existing_rows(target)
    if (
        rows
        and not (target / MANIFEST_FILENAME).is_file()
        and any(str(row.get("grid_kind")) == "LOCAL_REFINEMENT" for row in rows)
    ):
        # A prior finalization stopped after writing repair rows but before the
        # manifest made their metrics durable.  Preserve every BASE
        # transaction and repeat only those triggered local repairs.
        rows = [row for row in rows if str(row.get("grid_kind")) == "BASE"]
        for row in rows:
            row["is_canonical_plot_source"] = True
            if str(row.get("point_status")) in {
                "REPRODUCED_AFTER_LOCAL_REPAIR",
                "LOCATOR_CORRECTED_AFTER_LOCAL_REPAIR",
            }:
                row["point_status"] = "PASS"
        _atomic_write_csv(target / SPECTRUM_FILENAME, rows, SPECTRUM_FIELDS)
    required_final = (
        BENCHMARK_FILENAME,
        CHECKPOINT_FILENAME,
        SECTION_FILENAME,
        AUDIT_FILENAME,
        PLOT_FILENAME,
        REPORT_FILENAME,
        MANIFEST_FILENAME,
    )
    if (
        rows
        and audit_spectrum_rows(rows)["status"] == "PASS"
        and all((target / name).is_file() for name in required_final)
    ):
        return refresh_completed_outputs(target, rows)
    constitutive, _section_rows = prepare_constitutive_outputs(target)
    checkpoint_path = target / CHECKPOINT_FILENAME
    if checkpoint_path.is_file():
        checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
        if checkpoint.get("contract_sha256") != contract_hash():
            raise RuntimeError("Checkpoint contract differs from RLB-2E.")
        point_records = [dict(item) for item in checkpoint.get("point_records", [])]
    else:
        point_records = []
    rows, benchmark = run_benchmarks(
        target, rows, point_records, constitutive, started_at
    )
    if not benchmark["production_run_permitted"]:
        return {
            "stage_id": STAGE_ID,
            "scientific_status": "NOT_EVALUATED",
            "workflow_status": "STOPPED_BY_ETA_GATE",
            "benchmark": benchmark,
        }
    rows = complete_missing_points(
        target, rows, point_records, constitutive, started_at
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
    """Return the pre-root contract without mutating resume artifacts."""

    started = _utc_now()
    gate = constitutive_gate()
    rows = section_property_rows()
    payload = {
        "stage_id": STAGE_ID,
        "mode": "manifest_only",
        "scientific_status": "NOT_EVALUATED",
        "root_calculation_count": 0,
        "contract": contract_payload(),
        "constitutive_gate": gate,
        "section_property_row_count": len(rows),
        "output_directory": Path(output_dir).as_posix(),
        "resume_artifacts_modified": False,
        "created_at_utc": started,
    }
    return payload


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
        result = create_plot_from_csv(args.output_dir)
    elif args.manifest_only:
        result = manifest_only(args.output_dir)
    else:
        result = run_workflow(args.output_dir, missing_only=args.missing_only)
    print(json.dumps(_json_value(result), indent=2, ensure_ascii=False))
    return 0 if result.get("scientific_status", "PASS") != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
