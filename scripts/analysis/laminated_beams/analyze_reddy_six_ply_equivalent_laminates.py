"""RLB-2I six-ply globally equivalent laminate diagnostic.

The analysis constructs a one-parameter family of symmetric six-ply,
zero-degree laminates.  Their full membrane, bending, transverse-shear and
mass integrals are invariant, while the layerwise stress and strain-energy
distributions vary.  Only six coupled spectral spot checks are performed:
three values of ``zeta`` at ``beta = 0`` and ``30 deg``, roots 1--8 plus the
root-9 completeness guard.

The implementation is analysis-level orchestration.  Production laminate,
beam, geometry and coupled-joint modules are imported lazily and are never
modified.  ``--plot-only`` reads completed CSV files and performs no laminate
integration, matrix assembly, determinant/SVD evaluation or root search.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
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

STAGE_ID = "RLB-2I"
ALGORITHM_VERSION = "six_ply_equivalent_laminates_v1"

MU = 0.0
TAU = 0.0
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

PAIR_CENTER = "CENTER"
PAIR_MIDDLE = "MIDDLE"
PAIR_OUTER = "OUTER"
STACK_BOTTOM_TO_TOP = (
    PAIR_OUTER,
    PAIR_MIDDLE,
    PAIR_CENTER,
    PAIR_CENTER,
    PAIR_MIDDLE,
    PAIR_OUTER,
)
PAIR_ORDER = (PAIR_CENTER, PAIR_MIDDLE, PAIR_OUTER)
PAIR_A_WEIGHTS = {PAIR_CENTER: 1.0 / 3.0, PAIR_MIDDLE: 1.0 / 3.0, PAIR_OUTER: 1.0 / 3.0}
PAIR_D_WEIGHTS = {PAIR_CENTER: 1.0 / 27.0, PAIR_MIDDLE: 7.0 / 27.0, PAIR_OUTER: 19.0 / 27.0}

ZETA_MIN = -0.25
ZETA_MAX = 0.25
ZETA_STEP = 0.01
SPOT_ZETA_VALUES = (-0.25, 0.0, 0.25)
BETA_VALUES_DEG = (0.0, 30.0)
BOUNDARY_PROBE_OMEGA = (0.0, 0.731, 3.217)
K_PLOT = 8
K_GUARD = 9

MATRIX_RELATIVE_TOLERANCE = 1.0e-12
SYMMETRY_RELATIVE_TOLERANCE = 1.0e-12
REDUCED_PROPERTY_TOLERANCE = 1.0e-11
REDUCTION_ROUTE_TOLERANCE = 1.0e-11
RESULTANT_RECONSTRUCTION_TOLERANCE = 1.0e-12
ENERGY_IDENTITY_TOLERANCE = 1.0e-12
ENERGY_FRACTION_TOLERANCE = 1.0e-12
STRESS_PROPORTIONALITY_TOLERANCE = 1.0e-12
SPECTRAL_RELATIVE_TOLERANCE = 1.0e-9
ROOT_SINGULAR_RATIO_TOLERANCE = 1.0e-9
BOUNDARY_RESIDUAL_TOLERANCE = 1.0e-9

REFERENCE_AREA = WIDTH * THICKNESS
IY_REFERENCE = WIDTH * THICKNESS**3 / 12.0
OMEGA_TO_OMEGA_SCALE = L_REFERENCE**2 * math.sqrt(
    RHO_0 * REFERENCE_AREA / IY_REFERENCE
)

DEFAULT_OUTPUT_DIR = (
    ROOT / "results" / "laminated_beams" / "reddy_six_ply_equivalent_laminates"
)
SECTION_FILENAME = "section_equivalence.csv"
PLY_FILENAME = "ply_response.csv"
SPECTRAL_FILENAME = "spectral_spot_check.csv"
BOUNDARY_FILENAME = "boundary_matrix_equivalence.csv"
ENERGY_PLOT_FILENAME = "pair_energy_fractions_vs_zeta.png"
STRESS_PLOT_FILENAME = "pair_peak_sigma_vs_zeta.png"
PROFILE_PLOT_FILENAME = "bending_stress_profiles_selected_zeta.png"
REPORT_FILENAME = "report.md"
MANIFEST_FILENAME = "run_manifest.json"

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
}
EXPECTED_PREDECESSOR_TREE_HASHES = {
    "RLB-2A": "07D0B115FE8B42AC4EF11A32B3B37A075A5477203869BDB4A343F46C79533106",
    "RLB-2D": "86D34750EB13CE6039D8FFA18D9FE15A4CC518FCD7921A5646F3ADB0129F0250",
    "RLB-2E": "57E9FFCFD518FADF198F30C84F04EE181F1C645814A7BCFA834FCC920426B008",
    "RLB-2F": "10C80E8136AF917BCC5EFCB351FAD2FBE6665856A91C1F71241BC650372046C5",
    "RLB-2G": "4DA662EB77240C59B78017CCCD38136522561F2F3BE48BAD2FA50AACBB059CC1",
    "RLB-2H": "8669AFC4D8BF4831D1DBC295C4A6EABA4FA16A864B77810DAC5F652030ECA874",
}

SECTION_FIELDS = (
    "zeta",
    "s_center",
    "s_middle",
    "s_outer",
    "sum_A_identity",
    "sum_D_identity",
    "A_beam",
    "D_beam",
    "S_beam",
    "m",
    "J",
    "K",
    "width",
    "A_normalized",
    "D_normalized",
    "S_normalized",
    "m_normalized",
    "J_normalized",
    "A_matrix_relative",
    "D_matrix_relative",
    "shear_matrix_relative",
    "B_scaled_residual",
    "I1_scaled_residual",
    "beam_properties_max_relative",
    "reduction_route_max_relative",
    "z_interfaces",
    "A_matrix_per_unit_width",
    "B_matrix_per_unit_width",
    "D_matrix_per_unit_width",
    "shear_matrix_per_unit_width",
    "constitutive_status",
)

PLY_FIELDS = (
    "zeta",
    "load_case",
    "ply_index",
    "pair_id",
    "z_lower",
    "z_mid",
    "z_upper",
    "material_multiplier",
    "epsilon_lower",
    "epsilon_mid",
    "epsilon_upper",
    "sigma_lower",
    "sigma_mid",
    "sigma_upper",
    "ply_energy",
    "ply_energy_fraction",
    "pair_energy",
    "pair_energy_fraction",
    "analytical_pair_energy_fraction",
    "pair_fraction_residual",
    "maximum_abs_sigma_x",
    "maximum_sigma_x_coordinate",
    "maximum_sigma_x_coordinate_semantics",
    "pair_peak_abs_sigma_x",
    "pair_peak_normalized",
    "hotspot_flag",
    "canonical_hotspot_flag",
    "resultant_N_reconstructed",
    "resultant_M_reconstructed",
    "resultant_absolute_residual",
    "resultant_relative_residual",
    "energy_identity_rhs",
    "energy_identity_absolute_residual",
    "energy_identity_relative_residual",
)

SPECTRAL_FIELDS = (
    "beta_deg",
    "zeta",
    "sorted_position",
    "root_role",
    "guard_flag",
    "omega",
    "Omega",
    "Lambda",
    "reference_omega",
    "relative_frequency_difference",
    "raw_determinant",
    "scaled_determinant",
    "raw_sigma_ratio",
    "scaled_sigma_ratio",
    "boundary_null_residual",
    "detected_nullity",
    "cluster_id",
    "cluster_multiplicity",
    "cluster_total_nullity",
    "root_interval_left_Omega",
    "root_interval_right_Omega",
    "detector_refiner_provenance",
    "unresolved_candidates_below_root9",
    "root9_right_margin_Omega",
    "roots_above_9_computed",
    "quality_status",
)

BOUNDARY_FIELDS = (
    "beta_deg",
    "zeta",
    "omega",
    "state_matrix_max_abs_difference",
    "state_matrix_relative_difference",
    "boundary_matrix_max_abs_difference",
    "boundary_matrix_relative_difference",
    "reference_zeta",
    "status",
)


@dataclass(frozen=True)
class SectionObjects:
    laminate: Any
    properties: Any
    multipliers: Mapping[str, float]


ROOT_CALCULATION_COUNT = 0
LAMINATE_CALCULATION_COUNT = 0
MATRIX_ASSEMBLY_COUNT = 0


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
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
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
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", delete=False, dir=target.parent,
        prefix=f".{target.name}.", suffix=".tmp"
    ) as stream:
        stream.write(text)
        temporary = Path(stream.name)
    os.replace(temporary, target)


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> None:
    _atomic_write_text(
        path,
        json.dumps(_json_value(payload), ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )


def _atomic_write_csv(
    path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]
) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", newline="", delete=False, dir=target.parent,
        prefix=f".{target.name}.", suffix=".tmp"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fields), extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _csv_value(row.get(field, "")) for field in fields})
        temporary = Path(stream.name)
    os.replace(temporary, target)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _git_state() -> dict[str, str]:
    def run(*args: str) -> str:
        return subprocess.run(
            ["git", *args], cwd=ROOT, check=True, capture_output=True, text=True
        ).stdout.strip()

    return {
        "toplevel": run("rev-parse", "--show-toplevel"),
        "branch": run("branch", "--show-current"),
        "head": run("rev-parse", "HEAD"),
        "last_commit": run("log", "-1", "--oneline"),
        "status_short": run("status", "--short"),
    }


def _peak_rss_bytes() -> int:
    try:
        import psutil

        return int(psutil.Process().memory_info().rss)
    except (ImportError, OSError):
        return 0


def _physics_modules() -> tuple[Any, Any]:
    from scripts.lib import reddy_symmetric_coupled_beams as coupled
    from scripts.lib import reddy_symmetric_laminated_beam as beam

    return beam, coupled


def _root_helpers() -> tuple[Any, Any]:
    from scripts.analysis.laminated_beams import (
        sweep_reddy_axial_stiffness_visibility as rlb2h,
    )
    from scripts.analysis.laminated_beams import (
        sweep_reddy_stiffness_layout_contrast as rlb2e,
    )

    return rlb2e, rlb2h


def zeta_grid() -> FloatArray:
    return np.round(np.linspace(ZETA_MIN, ZETA_MAX, 51, dtype=float), 2)


def stiffness_multipliers(zeta: float) -> dict[str, float]:
    value = float(zeta)
    if not math.isfinite(value) or not ZETA_MIN <= value <= ZETA_MAX:
        raise ValueError("zeta must lie in [-0.25, 0.25].")
    values = {
        PAIR_CENTER: 1.0 + 2.0 * value,
        PAIR_MIDDLE: 1.0 - 3.0 * value,
        PAIR_OUTER: 1.0 + value,
    }
    if min(values.values()) <= 0.0:
        raise ValueError("All six-ply material multipliers must be positive.")
    return values


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
        name=f"RLB-2I {pair_id}(scale={factor:.6f})",
    )


def build_six_ply_section(zeta: float) -> SectionObjects:
    global LAMINATE_CALCULATION_COUNT
    LAMINATE_CALCULATION_COUNT += 1
    beam, _coupled = _physics_modules()
    multipliers = stiffness_multipliers(zeta)
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


def _relative(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / max(
        abs(float(left)), abs(float(right)), np.finfo(float).tiny
    )


def _matrix_relative(left: Any, right: Any) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    return float(np.linalg.norm(a - b, ord="fro") / max(
        np.linalg.norm(a, ord="fro"), np.linalg.norm(b, ord="fro"), np.finfo(float).tiny
    ))


def _scaled_B(laminate: Any) -> float:
    scale = max(
        float(np.linalg.norm(laminate.A, ord="fro")) * laminate.thickness,
        np.finfo(float).tiny,
    )
    return float(np.linalg.norm(laminate.B, ord="fro") / scale)


def _scaled_I1(laminate: Any) -> float:
    scale = max(abs(float(laminate.I0)) * laminate.thickness, np.finfo(float).tiny)
    return abs(float(laminate.I1)) / scale


def _reduction_route_max(properties: Any) -> float:
    reductions = (
        properties.axial_reduction,
        properties.bending_reduction,
        properties.shear_reduction_before_K,
    )
    return max(float(item.relative_difference) for item in reductions if item is not None)


def _beam_field_max_relative(left: Any, right: Any) -> float:
    return max(
        _relative(getattr(left, field), getattr(right, field))
        for field in ("A", "D", "S", "m", "J", "K", "width")
    )


def compute_constitutive_data() -> tuple[dict[str, Any], list[dict[str, Any]], dict[float, SectionObjects]]:
    sections = {float(value): build_six_ply_section(float(value)) for value in zeta_grid()}
    baseline = sections[0.0]
    beam, _coupled = _physics_modules()
    Q0 = np.asarray(beam.transformed_reduced_stiffness(build_scaled_material(1.0, "M0"), 0.0))
    shear0 = np.asarray(beam.transformed_transverse_shear_stiffness(build_scaled_material(1.0, "M0"), 0.0))
    expected_A0 = Q0 * THICKNESS
    expected_D0 = Q0 * THICKNESS**3 / 12.0
    expected_shear0 = shear0 * THICKNESS
    rows: list[dict[str, Any]] = []
    maxima = {
        "A_matrix_relative": 0.0,
        "D_matrix_relative": 0.0,
        "shear_matrix_relative": 0.0,
        "A_analytic_relative": 0.0,
        "D_analytic_relative": 0.0,
        "shear_analytic_relative": 0.0,
        "B_scaled_residual": 0.0,
        "I1_scaled_residual": 0.0,
        "beam_properties_max_relative": 0.0,
        "reduction_route_max_relative": 0.0,
        "exact_A_identity_absolute": 0.0,
        "exact_D_identity_absolute": 0.0,
    }
    for zeta, section in sections.items():
        s = section.multipliers
        identity_A = s[PAIR_CENTER] + s[PAIR_MIDDLE] + s[PAIR_OUTER]
        identity_D = s[PAIR_CENTER] + 7.0 * s[PAIR_MIDDLE] + 19.0 * s[PAIR_OUTER]
        values = {
            "A_matrix_relative": _matrix_relative(section.laminate.A, baseline.laminate.A),
            "D_matrix_relative": _matrix_relative(section.laminate.D, baseline.laminate.D),
            "shear_matrix_relative": _matrix_relative(section.laminate.shear, baseline.laminate.shear),
            "A_analytic_relative": _matrix_relative(section.laminate.A, expected_A0),
            "D_analytic_relative": _matrix_relative(section.laminate.D, expected_D0),
            "shear_analytic_relative": _matrix_relative(section.laminate.shear, expected_shear0),
            "B_scaled_residual": _scaled_B(section.laminate),
            "I1_scaled_residual": _scaled_I1(section.laminate),
            "beam_properties_max_relative": _beam_field_max_relative(section.properties, baseline.properties),
            "reduction_route_max_relative": _reduction_route_max(section.properties),
            "exact_A_identity_absolute": abs(identity_A - 3.0),
            "exact_D_identity_absolute": abs(identity_D - 27.0),
        }
        for name, value in values.items():
            maxima[name] = max(maxima[name], float(value))
        p = section.properties
        rows.append(
            {
                "zeta": zeta,
                "s_center": s[PAIR_CENTER],
                "s_middle": s[PAIR_MIDDLE],
                "s_outer": s[PAIR_OUTER],
                "sum_A_identity": identity_A,
                "sum_D_identity": identity_D,
                "A_beam": p.A,
                "D_beam": p.D,
                "S_beam": p.S,
                "m": p.m,
                "J": p.J,
                "K": p.K,
                "width": p.width,
                "A_normalized": p.A / baseline.properties.A,
                "D_normalized": p.D / baseline.properties.D,
                "S_normalized": p.S / baseline.properties.S,
                "m_normalized": p.m / baseline.properties.m,
                "J_normalized": p.J / baseline.properties.J,
                "A_matrix_relative": values["A_matrix_relative"],
                "D_matrix_relative": values["D_matrix_relative"],
                "shear_matrix_relative": values["shear_matrix_relative"],
                "B_scaled_residual": values["B_scaled_residual"],
                "I1_scaled_residual": values["I1_scaled_residual"],
                "beam_properties_max_relative": values["beam_properties_max_relative"],
                "reduction_route_max_relative": values["reduction_route_max_relative"],
                "z_interfaces": section.laminate.z_interfaces,
                "A_matrix_per_unit_width": section.laminate.A,
                "B_matrix_per_unit_width": section.laminate.B,
                "D_matrix_per_unit_width": section.laminate.D,
                "shear_matrix_per_unit_width": section.laminate.shear,
                "constitutive_status": "PASS",
            }
        )
    passed = bool(
        maxima["exact_A_identity_absolute"] <= 1.0e-14
        and maxima["exact_D_identity_absolute"] <= 1.0e-14
        and maxima["A_matrix_relative"] <= MATRIX_RELATIVE_TOLERANCE
        and maxima["D_matrix_relative"] <= MATRIX_RELATIVE_TOLERANCE
        and maxima["shear_matrix_relative"] <= MATRIX_RELATIVE_TOLERANCE
        and maxima["A_analytic_relative"] <= MATRIX_RELATIVE_TOLERANCE
        and maxima["D_analytic_relative"] <= MATRIX_RELATIVE_TOLERANCE
        and maxima["shear_analytic_relative"] <= MATRIX_RELATIVE_TOLERANCE
        and maxima["B_scaled_residual"] <= SYMMETRY_RELATIVE_TOLERANCE
        and maxima["I1_scaled_residual"] <= SYMMETRY_RELATIVE_TOLERANCE
        and maxima["beam_properties_max_relative"] <= REDUCED_PROPERTY_TOLERANCE
        and maxima["reduction_route_max_relative"] <= REDUCTION_ROUTE_TOLERANCE
    )
    for row in rows:
        row["constitutive_status"] = "PASS" if passed else "FAIL"
    gate = {
        "status": "PASS" if passed else "FAIL",
        "zeta_point_count": len(sections),
        "maxima": maxima,
        "baseline_beam_properties": {
            field: getattr(baseline.properties, field)
            for field in ("A", "D", "S", "m", "J", "K", "width")
        },
        "A_pair_weights": PAIR_A_WEIGHTS,
        "D_pair_weights": PAIR_D_WEIGHTS,
        "neutral_direction": [2, -3, 1],
    }
    return gate, rows, sections


def _vector_relative_residual(actual: FloatArray, expected: FloatArray) -> tuple[float, float]:
    absolute = float(np.linalg.norm(np.asarray(actual) - np.asarray(expected)))
    scale = max(float(np.linalg.norm(expected)), np.finfo(float).tiny)
    return absolute, absolute / scale


def _analytical_pair_fraction(
    load_case: str, pair_id: str, multipliers: Mapping[str, float]
) -> float:
    if load_case == "UNIT_AXIAL_RESULTANT":
        return float(multipliers[pair_id] * PAIR_A_WEIGHTS[pair_id])
    if load_case == "UNIT_BENDING_MOMENT":
        return float(multipliers[pair_id] * PAIR_D_WEIGHTS[pair_id])
    raise ValueError(f"Unsupported load case: {load_case}")


def _load_generalized_fields(
    section: SectionObjects, load_case: str
) -> tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """Return full-width resultants and generalized strains.

    Production ``LaminateSection.A`` and ``D`` are per unit width.  The unit
    vectors below are total beam resultants, so the unreduced full-width
    matrices ``b*A`` and ``b*D`` are used.  This makes the requested
    ``b/2`` ply-energy integral consistent with ``0.5*N.T@epsilon0`` and
    ``0.5*M.T@kappa``.
    """

    N = np.zeros(3, dtype=float)
    M = np.zeros(3, dtype=float)
    epsilon0 = np.zeros(3, dtype=float)
    kappa = np.zeros(3, dtype=float)
    if load_case == "UNIT_AXIAL_RESULTANT":
        N[0] = 1.0
        epsilon0 = np.linalg.solve(WIDTH * np.asarray(section.laminate.A), N)
    elif load_case == "UNIT_BENDING_MOMENT":
        M[0] = 1.0
        kappa = np.linalg.solve(WIDTH * np.asarray(section.laminate.D), M)
    else:
        raise ValueError(f"Unsupported load case: {load_case}")
    return N, M, epsilon0, kappa


def _ply_energy_exact(
    Qbar: FloatArray,
    epsilon0: FloatArray,
    kappa: FloatArray,
    z_lower: float,
    z_upper: float,
) -> float:
    dz1 = z_upper - z_lower
    dz2 = z_upper**2 - z_lower**2
    dz3 = z_upper**3 - z_lower**3
    eQe = float(epsilon0 @ Qbar @ epsilon0)
    eQk = float(epsilon0 @ Qbar @ kappa)
    kQk = float(kappa @ Qbar @ kappa)
    integral = eQe * dz1 + eQk * dz2 + (kQk / 3.0) * dz3
    return 0.5 * WIDTH * integral


def _ply_resultant_integrals(
    Qbar: FloatArray,
    epsilon0: FloatArray,
    kappa: FloatArray,
    z_lower: float,
    z_upper: float,
) -> tuple[FloatArray, FloatArray]:
    dz1 = z_upper - z_lower
    dz2 = z_upper**2 - z_lower**2
    dz3 = z_upper**3 - z_lower**3
    N = WIDTH * (Qbar @ epsilon0 * dz1 + 0.5 * Qbar @ kappa * dz2)
    M = WIDTH * (0.5 * Qbar @ epsilon0 * dz2 + (1.0 / 3.0) * Qbar @ kappa * dz3)
    return np.asarray(N, dtype=float), np.asarray(M, dtype=float)


def _one_load_response(
    zeta: float, section: SectionObjects, load_case: str
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    beam, _coupled = _physics_modules()
    N_target, M_target, epsilon0, kappa = _load_generalized_fields(section, load_case)
    provisional: list[dict[str, Any]] = []
    N_reconstructed = np.zeros(3, dtype=float)
    M_reconstructed = np.zeros(3, dtype=float)
    for index, ply in enumerate(section.laminate.plies):
        z_lower = float(section.laminate.z_interfaces[index])
        z_upper = float(section.laminate.z_interfaces[index + 1])
        z_mid = 0.5 * (z_lower + z_upper)
        Qbar = np.asarray(
            beam.transformed_reduced_stiffness(ply.material, ply.angle_deg), dtype=float
        )
        points = (z_lower, z_mid, z_upper)
        strains = [np.asarray(epsilon0 + z * kappa, dtype=float) for z in points]
        stresses = [np.asarray(Qbar @ strain, dtype=float) for strain in strains]
        sigma_abs = [abs(float(stress[0])) for stress in stresses]
        maximum = max(sigma_abs)
        if max(sigma_abs) - min(sigma_abs) <= 32.0 * np.finfo(float).eps * max(1.0, maximum):
            z_at_max = z_mid
            coordinate_semantics = "PLY_WIDE_CONSTANT"
        else:
            local_index = int(np.argmax(sigma_abs))
            z_at_max = points[local_index]
            coordinate_semantics = "ENDPOINT_MAXIMUM"
        energy = _ply_energy_exact(Qbar, epsilon0, kappa, z_lower, z_upper)
        N_part, M_part = _ply_resultant_integrals(
            Qbar, epsilon0, kappa, z_lower, z_upper
        )
        N_reconstructed += N_part
        M_reconstructed += M_part
        provisional.append(
            {
                "zeta": float(zeta),
                "load_case": load_case,
                "ply_index": index + 1,
                "pair_id": str(ply.label),
                "z_lower": z_lower,
                "z_mid": z_mid,
                "z_upper": z_upper,
                "material_multiplier": section.multipliers[str(ply.label)],
                "epsilon_lower": strains[0],
                "epsilon_mid": strains[1],
                "epsilon_upper": strains[2],
                "sigma_lower": stresses[0],
                "sigma_mid": stresses[1],
                "sigma_upper": stresses[2],
                "ply_energy": energy,
                "maximum_abs_sigma_x": maximum,
                "maximum_sigma_x_coordinate": z_at_max,
                "maximum_sigma_x_coordinate_semantics": coordinate_semantics,
            }
        )

    total_energy = math.fsum(float(row["ply_energy"]) for row in provisional)
    pair_energies = {
        pair_id: math.fsum(
            float(row["ply_energy"]) for row in provisional if row["pair_id"] == pair_id
        )
        for pair_id in PAIR_ORDER
    }
    pair_peaks = {
        pair_id: max(
            float(row["maximum_abs_sigma_x"])
            for row in provisional
            if row["pair_id"] == pair_id
        )
        for pair_id in PAIR_ORDER
    }
    target_energy = 0.5 * float(N_target @ epsilon0 + M_target @ kappa)
    energy_absolute = abs(total_energy - target_energy)
    energy_relative = energy_absolute / max(abs(target_energy), np.finfo(float).tiny)
    N_abs, _N_rel = _vector_relative_residual(N_reconstructed, N_target)
    M_abs, _M_rel = _vector_relative_residual(M_reconstructed, M_target)
    resultant_absolute = float(
        np.linalg.norm(
            np.concatenate((N_reconstructed - N_target, M_reconstructed - M_target))
        )
    )
    target_scale = max(
        float(np.linalg.norm(np.concatenate((N_target, M_target)))),
        np.finfo(float).tiny,
    )
    resultant_relative = resultant_absolute / target_scale
    global_peak = max(pair_peaks.values())
    tie_tolerance = 64.0 * np.finfo(float).eps * max(1.0, global_peak)
    hotspot_rows = [
        row
        for row in provisional
        if abs(float(row["maximum_abs_sigma_x"]) - global_peak) <= tie_tolerance
    ]
    canonical_index = min(int(row["ply_index"]) for row in hotspot_rows)

    rows: list[dict[str, Any]] = []
    for row in provisional:
        pair_id = str(row["pair_id"])
        pair_fraction = pair_energies[pair_id] / total_energy
        analytical_fraction = _analytical_pair_fraction(
            load_case, pair_id, section.multipliers
        )
        is_hotspot = abs(float(row["maximum_abs_sigma_x"]) - global_peak) <= tie_tolerance
        rows.append(
            {
                **row,
                "ply_energy_fraction": float(row["ply_energy"]) / total_energy,
                "pair_energy": pair_energies[pair_id],
                "pair_energy_fraction": pair_fraction,
                "analytical_pair_energy_fraction": analytical_fraction,
                "pair_fraction_residual": abs(pair_fraction - analytical_fraction),
                "pair_peak_abs_sigma_x": pair_peaks[pair_id],
                "pair_peak_normalized": 0.0,
                "hotspot_flag": is_hotspot,
                "canonical_hotspot_flag": is_hotspot and int(row["ply_index"]) == canonical_index,
                "resultant_N_reconstructed": N_reconstructed,
                "resultant_M_reconstructed": M_reconstructed,
                "resultant_absolute_residual": resultant_absolute,
                "resultant_relative_residual": resultant_relative,
                "energy_identity_rhs": target_energy,
                "energy_identity_absolute_residual": energy_absolute,
                "energy_identity_relative_residual": energy_relative,
            }
        )
    winning_pairs = [
        pair_id for pair_id in PAIR_ORDER if abs(pair_peaks[pair_id] - global_peak) <= tie_tolerance
    ]
    summary = {
        "zeta": float(zeta),
        "load_case": load_case,
        "epsilon0": epsilon0,
        "kappa": kappa,
        "N_reconstructed": N_reconstructed,
        "M_reconstructed": M_reconstructed,
        "N_absolute_residual": N_abs,
        "M_absolute_residual": M_abs,
        "resultant_absolute_residual": resultant_absolute,
        "resultant_relative_residual": resultant_relative,
        "total_energy": total_energy,
        "energy_identity_rhs": target_energy,
        "energy_identity_absolute_residual": energy_absolute,
        "energy_identity_relative_residual": energy_relative,
        "pair_energy_fractions": {
            pair_id: pair_energies[pair_id] / total_energy for pair_id in PAIR_ORDER
        },
        "pair_peaks": pair_peaks,
        "hotspot_pairs": winning_pairs,
        "hotspot_ply_indices": [int(row["ply_index"]) for row in hotspot_rows],
        "hotspot_z_coordinates": [float(row["maximum_sigma_x_coordinate"]) for row in hotspot_rows],
        "canonical_hotspot_ply_index": canonical_index,
        "hotspot_tie_count": len(hotspot_rows),
    }
    return rows, summary


def compute_ply_response_data(
    sections: Mapping[float, SectionObjects]
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for zeta in (float(value) for value in zeta_grid()):
        for load_case in ("UNIT_AXIAL_RESULTANT", "UNIT_BENDING_MOMENT"):
            load_rows, summary = _one_load_response(zeta, sections[zeta], load_case)
            rows.extend(load_rows)
            summaries.append(summary)

    reference_pair_peaks: dict[str, float] = {}
    for load_case in ("UNIT_AXIAL_RESULTANT", "UNIT_BENDING_MOMENT"):
        zero_rows = [
            row for row in rows if float(row["zeta"]) == 0.0 and row["load_case"] == load_case
        ]
        reference_pair_peaks[load_case] = max(
            float(row["maximum_abs_sigma_x"]) for row in zero_rows
        )
    for row in rows:
        row["pair_peak_normalized"] = (
            float(row["pair_peak_abs_sigma_x"])
            / reference_pair_peaks[str(row["load_case"])]
        )

    maximums = {
        "resultant_absolute_residual": max(
            float(item["resultant_absolute_residual"]) for item in summaries
        ),
        "resultant_relative_residual": max(
            float(item["resultant_relative_residual"]) for item in summaries
        ),
        "energy_identity_absolute_residual": max(
            float(item["energy_identity_absolute_residual"]) for item in summaries
        ),
        "energy_identity_relative_residual": max(
            float(item["energy_identity_relative_residual"]) for item in summaries
        ),
        "pair_fraction_residual": max(float(row["pair_fraction_residual"]) for row in rows),
        "pair_fraction_sum_residual": max(
            abs(math.fsum(item["pair_energy_fractions"].values()) - 1.0)
            for item in summaries
        ),
    }
    stress_proportionality_residual = 0.0
    for summary in summaries:
        zeta = float(summary["zeta"])
        multipliers = stiffness_multipliers(zeta)
        load_case = str(summary["load_case"])
        for pair_id in PAIR_ORDER:
            actual = float(summary["pair_peaks"][pair_id]) / reference_pair_peaks[load_case]
            expected = multipliers[pair_id]
            if load_case == "UNIT_BENDING_MOMENT":
                expected *= {PAIR_CENTER: 1.0 / 3.0, PAIR_MIDDLE: 2.0 / 3.0, PAIR_OUTER: 1.0}[pair_id]
            stress_proportionality_residual = max(
                stress_proportionality_residual, abs(actual - expected)
            )
    maximums["stress_proportionality_residual"] = stress_proportionality_residual

    axial_winners = {
        float(item["zeta"]): list(item["hotspot_pairs"])
        for item in summaries
        if item["load_case"] == "UNIT_AXIAL_RESULTANT"
    }
    bending_winners = {
        float(item["zeta"]): list(item["hotspot_pairs"])
        for item in summaries
        if item["load_case"] == "UNIT_BENDING_MOMENT"
    }
    axial_transition_pass = all(
        winners == ([PAIR_MIDDLE] if zeta < 0.0 else [PAIR_CENTER] if zeta > 0.0 else list(PAIR_ORDER))
        for zeta, winners in axial_winners.items()
    )
    bending_transition = -1.0 / 9.0
    bending_transition_pass = bool(
        bending_winners[-0.12] == [PAIR_MIDDLE]
        and bending_winners[-0.11] == [PAIR_OUTER]
        and all(PAIR_CENTER not in winners for winners in bending_winners.values())
    )
    passed = bool(
        maximums["resultant_relative_residual"] <= RESULTANT_RECONSTRUCTION_TOLERANCE
        and maximums["energy_identity_relative_residual"] <= ENERGY_IDENTITY_TOLERANCE
        and maximums["pair_fraction_residual"] <= ENERGY_FRACTION_TOLERANCE
        and maximums["pair_fraction_sum_residual"] <= ENERGY_FRACTION_TOLERANCE
        and maximums["stress_proportionality_residual"] <= STRESS_PROPORTIONALITY_TOLERANCE
        and axial_transition_pass
        and bending_transition_pass
    )
    gate = {
        "status": "PASS" if passed else "FAIL",
        "response_row_count": len(rows),
        "load_state_count": len(summaries),
        "full_width_resultant_convention": {
            "A_section": "width * LaminateSection.A",
            "D_section": "width * LaminateSection.D",
            "N_and_M": "total beam resultants",
            "ply_energy": "width/2 * integral(epsilon.T @ Qbar @ epsilon, z)",
        },
        "maximum_residuals": maximums,
        "reference_pair_peaks": reference_pair_peaks,
        "axial_hotspot_transition": {
            "parameter": 0.0,
            "status": "PASS" if axial_transition_pass else "FAIL",
            "negative_side_pair": PAIR_MIDDLE,
            "positive_side_pair": PAIR_CENTER,
            "zero_pair_tie": list(PAIR_ORDER),
        },
        "bending_hotspot_transition": {
            "analytical_parameter": bending_transition,
            "grid_bracket": [-0.12, -0.11],
            "lower_side_pair": PAIR_MIDDLE,
            "upper_side_pair": PAIR_OUTER,
            "status": "PASS" if bending_transition_pass else "FAIL",
        },
        "summaries": summaries,
    }
    return gate, rows, summaries


def make_matrix_provider(beta_deg: float, section: SectionObjects) -> Any:
    _beam, coupled = _physics_modules()
    beta_rad = math.radians(float(beta_deg))

    def provider(omega: float) -> FloatArray:
        global MATRIX_ASSEMBLY_COUNT
        MATRIX_ASSEMBLY_COUNT += 1
        return np.asarray(
            coupled.coupled_boundary_matrix(
                float(omega),
                beta_rad,
                L1,
                section.properties,
                L2,
                section.properties,
            ),
            dtype=float,
        )

    return provider


def boundary_matrix_equivalence_rows(
    sections: Mapping[float, SectionObjects]
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    beam, coupled = _physics_modules()
    baseline = sections[0.0]
    rows: list[dict[str, Any]] = []
    for beta_deg in BETA_VALUES_DEG:
        beta_rad = math.radians(beta_deg)
        for zeta in SPOT_ZETA_VALUES:
            section = sections[float(zeta)]
            for omega in BOUNDARY_PROBE_OMEGA:
                state = np.asarray(beam.combined_state_matrix(omega, section.properties))
                state0 = np.asarray(beam.combined_state_matrix(omega, baseline.properties))
                boundary = np.asarray(
                    coupled.coupled_boundary_matrix(
                        omega,
                        beta_rad,
                        L1,
                        section.properties,
                        L2,
                        section.properties,
                    )
                )
                boundary0 = np.asarray(
                    coupled.coupled_boundary_matrix(
                        omega,
                        beta_rad,
                        L1,
                        baseline.properties,
                        L2,
                        baseline.properties,
                    )
                )
                state_relative = _matrix_relative(state, state0)
                boundary_relative = _matrix_relative(boundary, boundary0)
                rows.append(
                    {
                        "beta_deg": beta_deg,
                        "zeta": zeta,
                        "omega": omega,
                        "state_matrix_max_abs_difference": float(np.max(np.abs(state - state0))),
                        "state_matrix_relative_difference": state_relative,
                        "boundary_matrix_max_abs_difference": float(np.max(np.abs(boundary - boundary0))),
                        "boundary_matrix_relative_difference": boundary_relative,
                        "reference_zeta": 0.0,
                        "status": (
                            "PASS"
                            if max(state_relative, boundary_relative)
                            <= MATRIX_RELATIVE_TOLERANCE
                            else "FAIL"
                        ),
                    }
                )
    maximum_state = max(float(row["state_matrix_relative_difference"]) for row in rows)
    maximum_boundary = max(float(row["boundary_matrix_relative_difference"]) for row in rows)
    passed = bool(
        maximum_state <= MATRIX_RELATIVE_TOLERANCE
        and maximum_boundary <= MATRIX_RELATIVE_TOLERANCE
        and all(row["status"] == "PASS" for row in rows)
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "row_count": len(rows),
        "probe_omega_values": list(BOUNDARY_PROBE_OMEGA),
        "maximum_state_matrix_relative_difference": maximum_state,
        "maximum_boundary_matrix_relative_difference": maximum_boundary,
    }, rows


def _solve_spectral_case(
    beta_deg: float, zeta: float, section: SectionObjects
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    global ROOT_CALCULATION_COUNT
    ROOT_CALCULATION_COUNT += 1
    rlb2e, rlb2h = _root_helpers()
    provider = make_matrix_provider(beta_deg, section)
    counted = rlb2e.CountedProvider(provider)
    policy = rlb2e._rlb2e_search_policy()
    solve_id = f"RLB2I_beta{int(beta_deg):02d}_zeta{float(zeta):+.2f}"
    started = time.perf_counter()
    canonical, rejected, slots, search_right, refinements = (
        rlb2e._progressive_anchor_candidates(counted, policy, solve_id=solve_id)
    )
    canonical, slots = rlb2h._truncate_inventory_to_root9(canonical, slots, policy)
    accepted, quality = rlb2h._point_is_acceptable_with_multiplicity(
        canonical, rejected, slots, search_right, policy
    )
    if not accepted:
        raise RuntimeError(f"{solve_id}: first-eight-plus-root9 gate failed: {quality}")
    rows: list[dict[str, Any]] = []
    for position, slot in enumerate(slots, start=1):
        event = slot.event
        candidate = event.candidate
        diagnostic = candidate.diagnostics
        rows.append(
            {
                "beta_deg": float(beta_deg),
                "zeta": float(zeta),
                "sorted_position": position,
                "root_role": "PLOTTED" if position <= K_PLOT else "ROOT_9_GUARD",
                "guard_flag": position == K_GUARD,
                "omega": float(event.omega),
                "Omega": float(event.omega_bar),
                "Lambda": math.sqrt(float(event.omega_bar)),
                "reference_omega": 0.0,
                "relative_frequency_difference": 0.0,
                "raw_determinant": diagnostic.raw_determinant,
                "scaled_determinant": diagnostic.scaled_determinant,
                "raw_sigma_ratio": diagnostic.raw_sigma_ratio,
                "scaled_sigma_ratio": diagnostic.scaled_sigma_ratio,
                "boundary_null_residual": diagnostic.raw_boundary_null_residual,
                "detected_nullity": diagnostic.detected_nullity,
                "cluster_id": event.cluster_id or event.event_id,
                "cluster_multiplicity": event.cluster_multiplicity,
                "cluster_total_nullity": event.cluster_total_nullity,
                "root_interval_left_Omega": candidate.interval_left_bar,
                "root_interval_right_Omega": candidate.interval_right_bar,
                "detector_refiner_provenance": candidate.detection_sources,
                "unresolved_candidates_below_root9": quality[
                    "unresolved_candidates_below_root9"
                ],
                "root9_right_margin_Omega": quality["root9_right_margin_Omega"],
                "roots_above_9_computed": quality["roots_above_9_computed"],
                "quality_status": "PASS",
            }
        )
    return rows, {
        "beta_deg": float(beta_deg),
        "zeta": float(zeta),
        "root_count": len(rows),
        "wall_time_seconds": time.perf_counter() - started,
        "peak_rss_bytes": _peak_rss_bytes(),
        "determinant_evaluations": counted.evaluations,
        "SVD_evaluations": counted.evaluations,
        "local_refinement_blocks": refinements,
        "search_right_Omega": search_right,
        "quality": quality,
    }


def spectral_spot_check(
    sections: Mapping[float, SectionObjects]
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    records: list[dict[str, Any]] = []
    for beta_deg in BETA_VALUES_DEG:
        for zeta in SPOT_ZETA_VALUES:
            case_rows, record = _solve_spectral_case(
                beta_deg, zeta, sections[float(zeta)]
            )
            rows.extend(case_rows)
            records.append(record)
    references = {
        (float(row["beta_deg"]), int(row["sorted_position"])): row
        for row in rows
        if float(row["zeta"]) == 0.0
    }
    maximum_relative = 0.0
    cluster_metadata_match = True
    for row in rows:
        reference = references[(float(row["beta_deg"]), int(row["sorted_position"]))]
        row["reference_omega"] = reference["omega"]
        difference = _relative(float(row["omega"]), float(reference["omega"]))
        row["relative_frequency_difference"] = difference
        maximum_relative = max(maximum_relative, difference)
        cluster_metadata_match = cluster_metadata_match and bool(
            int(row["cluster_multiplicity"]) == int(reference["cluster_multiplicity"])
            and int(row["cluster_total_nullity"])
            == int(reference["cluster_total_nullity"])
        )
    maximum_sigma = max(float(row["scaled_sigma_ratio"]) for row in rows)
    maximum_boundary = max(float(row["boundary_null_residual"]) for row in rows)
    group_keys = {(float(row["beta_deg"]), float(row["zeta"])) for row in rows}
    exact_groups = all(
        sorted(
            int(row["sorted_position"])
            for row in rows
            if (float(row["beta_deg"]), float(row["zeta"])) == key
        )
        == list(range(1, K_GUARD + 1))
        for key in group_keys
    )
    passed = bool(
        len(rows) == len(BETA_VALUES_DEG) * len(SPOT_ZETA_VALUES) * K_GUARD
        and len(group_keys) == 6
        and exact_groups
        and maximum_relative <= SPECTRAL_RELATIVE_TOLERANCE
        and maximum_sigma <= ROOT_SINGULAR_RATIO_TOLERANCE
        and maximum_boundary <= BOUNDARY_RESIDUAL_TOLERANCE
        and cluster_metadata_match
        and not any(bool(row["roots_above_9_computed"]) for row in rows)
        and all(int(row["unresolved_candidates_below_root9"]) == 0 for row in rows)
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "case_count": len(group_keys),
        "row_count": len(rows),
        "root_calculation_count": ROOT_CALCULATION_COUNT,
        "maximum_relative_frequency_difference": maximum_relative,
        "maximum_scaled_sigma_ratio": maximum_sigma,
        "maximum_boundary_null_residual": maximum_boundary,
        "cluster_metadata_match": cluster_metadata_match,
        "roots_above_9_computed": False,
        "case_records": records,
    }, rows


def _parse_json_vector(value: str) -> FloatArray:
    return np.asarray(json.loads(value), dtype=float)


def _audit_ply_csv_rows(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    keys = [
        (round(float(row["zeta"]), 2), str(row["load_case"]), int(row["ply_index"]))
        for row in rows
    ]
    expected = {
        (round(float(zeta), 2), load_case, ply_index)
        for zeta in zeta_grid()
        for load_case in ("UNIT_AXIAL_RESULTANT", "UNIT_BENDING_MOMENT")
        for ply_index in range(1, 7)
    }
    return {
        "row_count": len(rows),
        "duplicate_key_count": len(keys) - len(set(keys)),
        "missing_key_count": len(expected - set(keys)),
        "unexpected_key_count": len(set(keys) - expected),
        "status": (
            "PASS"
            if len(rows) == 612 and len(keys) == len(set(keys)) and set(keys) == expected
            else "FAIL"
        ),
    }


def _save_figure_atomic(figure: Any, path: Path) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.parent / f".{target.name}.{os.getpid()}.tmp.png"
    figure.savefig(temporary, dpi=180, bbox_inches="tight")
    os.replace(temporary, target)


def create_figures_from_csv(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> dict[str, Any]:
    target = Path(output_dir)
    rows = _read_csv(target / PLY_FILENAME)
    audit = _audit_ply_csv_rows(rows)
    if audit["status"] != "PASS":
        raise RuntimeError(f"plot-only rejected incomplete ply response CSV: {audit}")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    colors = {PAIR_CENTER: "#0072B2", PAIR_MIDDLE: "#D55E00", PAIR_OUTER: "#009E73"}
    labels = {PAIR_CENTER: "CENTER", PAIR_MIDDLE: "MIDDLE", PAIR_OUTER: "OUTER"}
    unique_pair_rows: dict[tuple[float, str, str], Mapping[str, Any]] = {}
    for row in rows:
        key = (round(float(row["zeta"]), 2), str(row["load_case"]), str(row["pair_id"]))
        unique_pair_rows.setdefault(key, row)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), sharex=True, sharey=True)
    for axis, load_case, title in zip(
        axes,
        ("UNIT_AXIAL_RESULTANT", "UNIT_BENDING_MOMENT"),
        ("(a) Unit axial resultant", "(b) Unit bending moment"),
    ):
        for pair_id in PAIR_ORDER:
            selected = sorted(
                (
                    (key[0], float(row["pair_energy_fraction"]))
                    for key, row in unique_pair_rows.items()
                    if key[1] == load_case and key[2] == pair_id
                ),
                key=lambda item: item[0],
            )
            axis.plot(
                [item[0] for item in selected],
                [item[1] for item in selected],
                color=colors[pair_id],
                label=labels[pair_id],
                linewidth=1.8,
            )
        axis.set_title(title)
        axis.set_xlabel(r"$\zeta$")
        axis.grid(True, alpha=0.25)
    axes[0].set_ylabel("Energy fraction")
    axes[1].legend(frameon=False, loc="best")
    fig.tight_layout()
    _save_figure_atomic(fig, target / ENERGY_PLOT_FILENAME)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), sharex=True)
    for axis, load_case, title in zip(
        axes,
        ("UNIT_AXIAL_RESULTANT", "UNIT_BENDING_MOMENT"),
        ("(a) Unit axial resultant", "(b) Unit bending moment"),
    ):
        for pair_id in PAIR_ORDER:
            selected = sorted(
                (
                    (key[0], float(row["pair_peak_normalized"]))
                    for key, row in unique_pair_rows.items()
                    if key[1] == load_case and key[2] == pair_id
                ),
                key=lambda item: item[0],
            )
            axis.plot(
                [item[0] for item in selected],
                [item[1] for item in selected],
                color=colors[pair_id],
                label=labels[pair_id],
                linewidth=1.8,
            )
        axis.axvline(0.0, color="0.45", linewidth=0.8, linestyle=":")
        if load_case == "UNIT_BENDING_MOMENT":
            axis.axvline(-1.0 / 9.0, color="0.2", linewidth=0.8, linestyle="--")
        axis.set_title(title)
        axis.set_xlabel(r"$\zeta$")
        axis.grid(True, alpha=0.25)
    axes[0].set_ylabel(r"Normalized pair peak $|\sigma_x|$")
    axes[1].legend(frameon=False, loc="best")
    fig.tight_layout()
    _save_figure_atomic(fig, target / STRESS_PLOT_FILENAME)
    plt.close(fig)

    selected_profiles: dict[float, list[Mapping[str, Any]]] = {}
    all_sigma: list[float] = []
    for zeta in SPOT_ZETA_VALUES:
        profile = sorted(
            (
                row
                for row in rows
                if round(float(row["zeta"]), 2) == round(zeta, 2)
                and row["load_case"] == "UNIT_BENDING_MOMENT"
            ),
            key=lambda row: int(row["ply_index"]),
        )
        selected_profiles[zeta] = profile
        for row in profile:
            all_sigma.extend(
                [
                    float(_parse_json_vector(row["sigma_lower"])[0]),
                    float(_parse_json_vector(row["sigma_upper"])[0]),
                ]
            )
    sigma_limit = 1.05 * max(abs(value) for value in all_sigma)
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.5), sharex=True, sharey=True)
    for axis, zeta, panel in zip(axes, SPOT_ZETA_VALUES, ("(a)", "(b)", "(c)")):
        for row in selected_profiles[zeta]:
            pair_id = str(row["pair_id"])
            sigma_lower = float(_parse_json_vector(row["sigma_lower"])[0])
            sigma_upper = float(_parse_json_vector(row["sigma_upper"])[0])
            axis.plot(
                [sigma_lower, sigma_upper],
                [float(row["z_lower"]), float(row["z_upper"])],
                color=colors[pair_id],
                linewidth=2.0,
            )
        for interface in np.linspace(-THICKNESS / 2.0, THICKNESS / 2.0, 7):
            axis.axhline(interface, color="0.75", linewidth=0.55, linestyle=":")
        axis.axvline(0.0, color="0.45", linewidth=0.7)
        axis.set_xlim(-sigma_limit, sigma_limit)
        axis.set_ylim(-THICKNESS / 2.0, THICKNESS / 2.0)
        axis.set_title(fr"{panel} $\zeta={zeta:+.2f}$")
        axis.set_xlabel(r"$\sigma_x$")
        axis.grid(True, alpha=0.18)
    axes[0].set_ylabel(r"$z_{\mathrm{project}}$")
    handles = [Line2D([0], [0], color=colors[pair], lw=2, label=labels[pair]) for pair in PAIR_ORDER]
    axes[-1].legend(handles=handles, frameon=False, loc="best")
    fig.tight_layout()
    _save_figure_atomic(fig, target / PROFILE_PLOT_FILENAME)
    plt.close(fig)

    return {
        "status": "PASS",
        "input_audit": audit,
        "figures": {
            ENERGY_PLOT_FILENAME: {"panels": 2},
            STRESS_PLOT_FILENAME: {"panels": 2},
            PROFILE_PLOT_FILENAME: {"panels": 3, "piecewise_segments": 18},
        },
        "root_calculation_count": 0,
        "laminate_recalculation_count": 0,
        "matrix_assembly_count": 0,
        "determinant_evaluation_count": 0,
        "SVD_evaluation_count": 0,
        "smoothing_applied": False,
    }


def contract_payload() -> dict[str, Any]:
    return {
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "geometry": {
            "mu": MU,
            "tau": TAU,
            "l": L_REFERENCE,
            "L1": L1,
            "L2": L2,
            "L_total": L_TOTAL,
            "width": WIDTH,
            "thickness": THICKNESS,
            "ply_thickness": PLY_THICKNESS,
            "K": K,
            "spectral_beta_values_deg": list(BETA_VALUES_DEG),
            "coordinate": "z_project=-z_Reddy",
        },
        "material_M0": base_material_contract(),
        "six_ply_stack": {
            "bottom_to_top": list(STACK_BOTTOM_TO_TOP),
            "pair_order_from_midplane": list(PAIR_ORDER),
            "ply_angles_deg": [0.0] * 6,
            "equal_ply_thicknesses": [PLY_THICKNESS] * 6,
        },
        "zeta": {
            "grid": [float(value) for value in zeta_grid()],
            "range": [ZETA_MIN, ZETA_MAX],
            "step": ZETA_STEP,
            "s_center": "1+2*zeta",
            "s_middle": "1-3*zeta",
            "s_outer": "1+zeta",
            "neutral_direction_center_middle_outer": [2, -3, 1],
        },
        "material_scaling": {
            "scaled": ["E1", "E2", "G12"],
            "fixed": ["nu12", "G13", "G23", "rho"],
            "all_Qbar_angles_deg": 0.0,
        },
        "exact_equivalence": {
            "A_pair_weights": PAIR_A_WEIGHTS,
            "D_pair_weights": PAIR_D_WEIGHTS,
            "A_identity": "s_center+s_middle+s_outer=3",
            "D_identity": "s_center+7*s_middle+19*s_outer=27",
        },
        "load_cases": {
            "UNIT_AXIAL_RESULTANT": {"N": [1.0, 0.0, 0.0], "M": [0.0, 0.0, 0.0]},
            "UNIT_BENDING_MOMENT": {"N": [0.0, 0.0, 0.0], "M": [1.0, 0.0, 0.0]},
            "resultant_semantics": "total beam resultants",
            "full_width_matrices": {
                "A_section": "width*LaminateSection.A",
                "D_section": "width*LaminateSection.D",
            },
            "reason": (
                "LaminateSection matrices are per unit width; multiplying by width "
                "keeps the requested width-weighted energy identity dimensionally consistent."
            ),
        },
        "spectral_spot_check": {
            "zeta_values": list(SPOT_ZETA_VALUES),
            "beta_values_deg": list(BETA_VALUES_DEG),
            "positions": list(range(1, K_PLOT + 1)),
            "guard_position": K_GUARD,
            "guard_role": "completeness_only",
            "root_calculation_case_count": 6,
            "wide_zeta_frequency_sweep": False,
            "roots_10_plus": False,
        },
        "boundary_matrix_probes": {
            "physical_omega_values": list(BOUNDARY_PROBE_OMEGA),
            "beta_values_deg": list(BETA_VALUES_DEG),
            "zeta_values": list(SPOT_ZETA_VALUES),
        },
        "thresholds": {
            "matrix_relative": MATRIX_RELATIVE_TOLERANCE,
            "symmetry_relative": SYMMETRY_RELATIVE_TOLERANCE,
            "reduced_property_relative": REDUCED_PROPERTY_TOLERANCE,
            "reduction_route_relative": REDUCTION_ROUTE_TOLERANCE,
            "resultant_reconstruction_relative": RESULTANT_RECONSTRUCTION_TOLERANCE,
            "energy_identity_relative": ENERGY_IDENTITY_TOLERANCE,
            "energy_fraction_absolute": ENERGY_FRACTION_TOLERANCE,
            "stress_proportionality_absolute": STRESS_PROPORTIONALITY_TOLERANCE,
            "spectral_relative": SPECTRAL_RELATIVE_TOLERANCE,
            "root_singular_ratio": ROOT_SINGULAR_RATIO_TOLERANCE,
            "boundary_residual": BOUNDARY_RESIDUAL_TOLERANCE,
        },
        "exclusions": {
            "wide_frequency_sweep": False,
            "roots_10_plus": False,
            "new_physical_solver": False,
            "production_physics_changes": False,
            "branch_tracking": False,
            "MAC": False,
            "coupled_mode_shapes": False,
            "Ritz": False,
            "FEM": False,
            "failure_criteria": False,
            "damage_modelling": False,
            "delamination": False,
            "stress_smoothing": False,
            "interpolation_based_roots": False,
            "commit": False,
            "push": False,
        },
    }


def contract_hash() -> str:
    return hashlib.sha256(
        json.dumps(
            _json_value(contract_payload()),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest().upper()


def _endpoint_fraction_table(
    summaries: Sequence[Mapping[str, Any]]
) -> list[dict[str, Any]]:
    table: list[dict[str, Any]] = []
    for zeta in SPOT_ZETA_VALUES:
        record: dict[str, Any] = {"zeta": zeta}
        for load_case, prefix in (
            ("UNIT_AXIAL_RESULTANT", "axial"),
            ("UNIT_BENDING_MOMENT", "bending"),
        ):
            summary = next(
                item
                for item in summaries
                if float(item["zeta"]) == zeta and item["load_case"] == load_case
            )
            record[prefix] = dict(summary["pair_energy_fractions"])
            record[f"{prefix}_hotspot_pairs"] = list(summary["hotspot_pairs"])
        table.append(record)
    return table


def _stress_jump_summary(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    for zeta in SPOT_ZETA_VALUES:
        profile = sorted(
            (
                row
                for row in rows
                if float(row["zeta"]) == zeta
                and row["load_case"] == "UNIT_BENDING_MOMENT"
            ),
            key=lambda row: int(row["ply_index"]),
        )
        for left, right in zip(profile[:-1], profile[1:]):
            left_sigma = float(np.asarray(left["sigma_upper"])[0])
            right_sigma = float(np.asarray(right["sigma_lower"])[0])
            records.append(
                {
                    "zeta": zeta,
                    "z_interface": float(left["z_upper"]),
                    "left_ply": int(left["ply_index"]),
                    "right_ply": int(right["ply_index"]),
                    "left_sigma_x": left_sigma,
                    "right_sigma_x": right_sigma,
                    "jump": right_sigma - left_sigma,
                }
            )
    nonzero = [record for record in records if abs(float(record["jump"])) > 1.0e-12]
    return {
        "selected_zeta_values": list(SPOT_ZETA_VALUES),
        "interface_record_count": len(records),
        "nonzero_jump_count": len(nonzero),
        "maximum_absolute_jump": max(abs(float(record["jump"])) for record in records),
        "piecewise_segments_retained": True,
        "smoothing_applied": False,
        "records": records,
    }


def _output_hashes(output_dir: Path) -> dict[str, str]:
    return {
        path.name: _sha256(path)
        for path in sorted(Path(output_dir).iterdir())
        if path.is_file() and path.name != MANIFEST_FILENAME
    }


def _report_text(manifest: Mapping[str, Any]) -> str:
    constitutive = manifest["constitutive_gate"]
    response = manifest["layerwise_response_gate"]
    spectral = manifest["spectral_spot_check"]
    boundary = manifest["boundary_matrix_equivalence"]
    base = constitutive["baseline_beam_properties"]
    endpoint = manifest["endpoint_energy_and_hotspot_summary"]

    def fractions(record: Mapping[str, Any], key: str) -> str:
        values = record[key]
        return ", ".join(
            f"{pair}={float(values[pair]):.8f}" for pair in PAIR_ORDER
        )

    endpoint_lines = []
    for record in endpoint:
        endpoint_lines.append(
            f"| {float(record['zeta']):+.2f} | {fractions(record, 'axial')} | "
            f"{fractions(record, 'bending')} |"
        )
    endpoint_table = "\n".join(endpoint_lines)
    hotspot = response["bending_hotspot_transition"]
    maxima = response["maximum_residuals"]
    statuses = manifest["status_gates"]
    status_lines = "\n".join(f"- `{name}`: **{value}**" for name, value in statuses.items())
    return rf"""# RLB-2I: globally equivalent six-ply laminates

## Scope

Рассматривается однопараметрическое семейство симметричных шестислойных
ламинатов с одинаковыми редуцированными свойствами. Три зеркальные пары дают
три независимых множителя жёсткости, а условия сохранения матриц $A$ и $D$
задают две независимые связи. Поэтому шесть слоёв являются минимальным
нетривиальным симметричным случаем с одним свободным направлением.

Все слои ориентированы под углом $0^\circ$ и имеют толщину $h/6$. Порядок
снизу вверх равен `OUTER/MIDDLE/CENTER/CENTER/MIDDLE/OUTER`. Использована
координата $z_{{\mathrm{{project}}}}=-z_{{\mathrm{{Reddy}}}}$.

## Exact equivalence

Множители имеют вид

\[
s_C=1+2\zeta,\qquad s_M=1-3\zeta,\qquad s_O=1+\zeta,
\]

где $\zeta=-0.25:0.01:0.25$. Веса пар в матрице $A$ равны
$(1,1,1)/3$, а в матрице $D$ -- $(1,7,19)/27$. Направление
$(2,-3,1)$ сохраняет обе матрицы, поскольку

\[
s_C+s_M+s_O=3,\qquad s_C+7s_M+19s_O=27.
\]

Constitutive gate получил статус `{constitutive['status']}`. Максимальные
relative residuals матриц $A$, $D$ и transverse shear равны соответственно
{constitutive['maxima']['A_matrix_relative']:.3e},
{constitutive['maxima']['D_matrix_relative']:.3e} и
{constitutive['maxima']['shear_matrix_relative']:.3e}. Максимальные scaled
residuals $B$ и $I_1$ равны
{constitutive['maxima']['B_scaled_residual']:.3e} и
{constitutive['maxima']['I1_scaled_residual']:.3e}.

Редуцированные свойства во всём семействе совпадают:

- $A_{{\rm beam}}={base['A']:.16g}$;
- $D_{{\rm beam}}={base['D']:.16g}$;
- $S_{{\rm beam}}={base['S']:.16g}$;
- $m={base['m']:.16g}$;
- $J={base['J']:.16g}$.

## Layerwise reconstruction

Матрицы `LaminateSection.A/D` заданы на единицу ширины. Единичные $N_x$ и
$M_x$ в этом расчёте трактуются как полные результирующие стержня. Поэтому
использованы полные матрицы $A_{{\rm sec}}=bA$ и $D_{{\rm sec}}=bD$. Такая
конвенция согласует послойную энергию
$b\int\varepsilon^T\bar Q\varepsilon\,dz/2$ с тождеством
$U=(N^T\varepsilon_0+M^T\kappa)/2$.

Максимальный normalized residual восстановления результирующих равен
{maxima['resultant_relative_residual']:.3e}. Для энергетического тождества
получено {maxima['energy_identity_relative_residual']:.3e}. Максимальная
разность вычисленной и аналитической доли энергии пары равна
{maxima['pair_fraction_residual']:.3e}.

| $\zeta$ | Axial pair fractions | Bending pair fractions |
| ---: | --- | --- |
{endpoint_table}

При осевом нагружении максимальное $|\sigma_x|$ находится в паре `MIDDLE`
при $\zeta<0$ и в паре `CENTER` при $\zeta>0$. При $\zeta=0$ все пары имеют
одинаковый уровень. Для изгиба переход между `MIDDLE` и `OUTER` происходит
при $\zeta=-1/9={float(hotspot['analytical_parameter']):.12g}$; на принятой
сетке он заключён между -0.12 и -0.11. Следовательно, при неизменной
глобальной изгибной жёсткости место максимального изгибного напряжения может
перемещаться вследствие перераспределения жёсткости между слоями.

Из-за зеркальной симметрии максимум $|\sigma_x|$ не является уникальным по
слою. В CSV отмечены все связанные слои, а один deterministic representative
хранится только как служебный selector. Межслойные скачки напряжения сохранены
как односторонние значения; сглаживание не применялось.

## Spectral regression

Широкая частотная карта по $\zeta$ не строилась. Выполнены только шесть
spot-check cases: $\zeta\in\{{-0.25,0,+0.25\}}$ при
$\beta\in\{{0^\circ,30^\circ\}}$. В каждом случае сохранены positions 1--8
и root 9 как completeness guard. Корни 10 и выше не вычислялись.

Максимальная relative difference относительно $\zeta=0$ равна
{spectral['maximum_relative_frequency_difference']:.3e}. Максимальные
$\sigma_{{\min}}/\sigma_{{\max}}$ и boundary residual равны
{spectral['maximum_scaled_sigma_ratio']:.3e} и
{spectral['maximum_boundary_null_residual']:.3e}. Boundary matrices на
фиксированных частотах совпали с максимальным relative residual
{boundary['maximum_boundary_matrix_relative_difference']:.3e}.

## Figures

- `{ENERGY_PLOT_FILENAME}` -- доли энергии трёх зеркальных пар;
- `{STRESS_PLOT_FILENAME}` -- нормированные максимумы $|\sigma_x|$;
- `{PROFILE_PLOT_FILENAME}` -- кусочно-линейные изгибные профили для трёх
  выбранных значений $\zeta$.

## Status and limitations

{status_lines}

Главный вывод ограничен объявленным шестислойным семейством, принятой
одномерной RLB-моделью, двумя load-controlled состояниями и шестью
спектральными spot checks. Полученные данные не задают критерий прочности и не
позволяют делать выводы о начале разрушения, усталостной долговечности,
расслоении или повреждаемости. Branch tracking, MAC, coupled mode shapes,
Ritz и FEM не выполнялись.
"""


def manifest_only(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    return {
        "stage_id": STAGE_ID,
        "mode": "manifest_only",
        "scientific_status": "NOT_EVALUATED",
        "contract": contract_payload(),
        "contract_sha256": contract_hash(),
        "output_directory": Path(output_dir).as_posix(),
        "root_calculation_count": 0,
        "laminate_calculation_count": 0,
        "matrix_assembly_count": 0,
        "created_at_utc": _utc_now(),
    }


def plot_only_workflow(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    before_roots = ROOT_CALCULATION_COUNT
    before_laminates = LAMINATE_CALCULATION_COUNT
    before_matrices = MATRIX_ASSEMBLY_COUNT
    result = create_figures_from_csv(output_dir)
    if (
        ROOT_CALCULATION_COUNT != before_roots
        or LAMINATE_CALCULATION_COUNT != before_laminates
        or MATRIX_ASSEMBLY_COUNT != before_matrices
    ):
        raise RuntimeError("plot-only unexpectedly entered a calculation path.")
    return {"mode": "plot_only", **result}


def run_workflow(output_dir: Path = DEFAULT_OUTPUT_DIR) -> dict[str, Any]:
    total_started = time.perf_counter()
    started_at = _utc_now()
    target = Path(output_dir)
    target.mkdir(parents=True, exist_ok=True)

    physics_hashes = {path: _sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS}
    if physics_hashes != EXPECTED_PRODUCTION_PHYSICS_HASHES:
        raise RuntimeError("Frozen production-physics hash mismatch before RLB-2I.")
    predecessor_hashes = {
        name: _sha_tree(path) for name, path in PREDECESSOR_RESULT_DIRS.items()
    }
    if predecessor_hashes != EXPECTED_PREDECESSOR_TREE_HASHES:
        raise RuntimeError("Predecessor result-tree hash mismatch before RLB-2I.")

    constitutive_gate, section_rows, sections = compute_constitutive_data()
    _atomic_write_csv(target / SECTION_FILENAME, section_rows, SECTION_FIELDS)
    if constitutive_gate["status"] != "PASS":
        raise RuntimeError("RLB-2I constitutive equivalence gate failed before roots.")

    response_gate, ply_rows, response_summaries = compute_ply_response_data(sections)
    _atomic_write_csv(target / PLY_FILENAME, ply_rows, PLY_FIELDS)
    if response_gate["status"] != "PASS":
        raise RuntimeError("RLB-2I layerwise resultant/energy gate failed.")

    boundary_gate, boundary_rows = boundary_matrix_equivalence_rows(sections)
    _atomic_write_csv(target / BOUNDARY_FILENAME, boundary_rows, BOUNDARY_FIELDS)
    if boundary_gate["status"] != "PASS":
        raise RuntimeError("RLB-2I boundary-matrix equivalence gate failed before roots.")

    spectral_gate, spectral_rows = spectral_spot_check(sections)
    _atomic_write_csv(target / SPECTRAL_FILENAME, spectral_rows, SPECTRAL_FIELDS)

    figures = create_figures_from_csv(target)
    all_primary_gates_pass = bool(
        constitutive_gate["status"] == "PASS"
        and response_gate["status"] == "PASS"
        and boundary_gate["status"] == "PASS"
        and spectral_gate["status"] == "PASS"
        and figures["status"] == "PASS"
    )
    status_gates = {
        "RLB-2I-EXACT-EQUIVALENCE": constitutive_gate["status"],
        "RLB-2I-LAYERWISE-RECONSTRUCTION": response_gate["status"],
        "RLB-2I-BOUNDARY-MATRIX-EQUIVALENCE": boundary_gate["status"],
        "RLB-2I-SPECTRAL-SPOT-CHECK": spectral_gate["status"],
        "RLB-2I-FIGURES": figures["status"],
        "OVERALL": "PASS" if all_primary_gates_pass else "FAIL",
    }
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "stage_id": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "scientific_status": "PASS" if all_primary_gates_pass else "FAIL",
        "permitted_conclusion": (
            "Distinct six-ply laminates can have identical reduced properties and "
            "the same spectrum in the current RLB model while retaining different "
            "layerwise stresses and energy distributions."
        ),
        "started_at_utc": started_at,
        "completed_at_utc": _utc_now(),
        "git": _git_state(),
        "contract": contract_payload(),
        "contract_sha256": contract_hash(),
        "analysis_script_sha256": _sha256(Path(__file__)),
        "production_physics_hashes": physics_hashes,
        "production_physics_preserved": physics_hashes == EXPECTED_PRODUCTION_PHYSICS_HASHES,
        "predecessor_result_tree_hashes": predecessor_hashes,
        "predecessor_result_trees_preserved": predecessor_hashes == EXPECTED_PREDECESSOR_TREE_HASHES,
        "constitutive_gate": constitutive_gate,
        "layerwise_response_gate": response_gate,
        "boundary_matrix_equivalence": boundary_gate,
        "spectral_spot_check": spectral_gate,
        "stress_jump_diagnostic": _stress_jump_summary(ply_rows),
        "endpoint_energy_and_hotspot_summary": _endpoint_fraction_table(response_summaries),
        "figures": figures,
        "status_gates": status_gates,
        "counts": {
            "section_equivalence_rows": len(section_rows),
            "ply_response_rows": len(ply_rows),
            "spectral_spot_check_rows": len(spectral_rows),
            "boundary_matrix_equivalence_rows": len(boundary_rows),
            "spectral_cases": 6,
            "figures": 3,
        },
        "runtime": {
            "wall_time_seconds": time.perf_counter() - total_started,
            "peak_rss_bytes": _peak_rss_bytes(),
            "root_calculation_count": ROOT_CALCULATION_COUNT,
            "laminate_calculation_count": LAMINATE_CALCULATION_COUNT,
            "boundary_matrix_assembly_count": MATRIX_ASSEMBLY_COUNT,
            "determinant_evaluations": sum(
                int(record["determinant_evaluations"])
                for record in spectral_gate["case_records"]
            ),
            "SVD_evaluations": sum(
                int(record["SVD_evaluations"])
                for record in spectral_gate["case_records"]
            ),
            "parallel_workers": 0,
        },
        "explicit_confirmations": {
            "wide_51_point_frequency_sweep_run": False,
            "spectral_spot_check_case_count": 6,
            "roots_above_9_computed": False,
            "new_physical_solver_created": False,
            "production_physics_changed": False,
            "branch_tracking_run": False,
            "MAC_run": False,
            "coupled_mode_shapes_reconstructed": False,
            "Ritz_run": False,
            "FEM_run": False,
            "failure_or_damage_model_run": False,
            "stress_smoothing_applied": False,
            "interpolation_based_root_values": False,
            "commit_performed": False,
            "push_performed": False,
        },
    }
    _atomic_write_text(target / REPORT_FILENAME, _report_text(manifest))
    manifest["output_hashes"] = _output_hashes(target)
    _atomic_write_json(target / MANIFEST_FILENAME, manifest)
    return manifest


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--plot-only", action="store_true")
    mode.add_argument("--manifest-only", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    if args.plot_only:
        result = plot_only_workflow(args.output_dir)
    elif args.manifest_only:
        result = manifest_only(args.output_dir)
    else:
        result = run_workflow(args.output_dir)
    print(json.dumps(_json_value(result), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
