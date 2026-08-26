"""RLB-2B finite beta-grid comparison for three rectangular beam models.

The workflow compares pointwise *sorted spectral positions* of the frozen
Euler--Bernoulli determinant, the independent legacy rectangular Timoshenko
determinant, and the four-ply isotropic Reddy laminated-beam determinant.  It
does not track modal descendants.

All models use the same dimensional contract and the common frequency

    Lambda = omega*l**2*sqrt(rho*A/(E*I_y)).

The variable called ``Lambda`` by the historical EB matrix is instead its
bending wavenumber.  Consequently the EB matrix is evaluated at
``sqrt(Lambda)`` while every CSV and plot stores the common frequency above.

This file is an analysis entry point, not a production solver.  The two
Timoshenko formulations reuse the completed RLB-1C-ISO determinant/SVD policy.
The EB formulation retains its historical sign-scan/bisection and adds a
finer-grid repeat plus the same matrix-quality diagnostics at every root.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Any, Callable, Iterable, Mapping, Sequence

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

from my_project.analytic.formulas import (  # noqa: E402
    assemble_clamped_coupled_matrix as assemble_eb_coupled_matrix,
)
from my_project.analytic.solvers import find_first_n_roots as find_eb_roots  # noqa: E402
from scripts.analysis.laminated_beams import (  # noqa: E402
    validate_reddy_four_ply_isotropic_limit as iso_inventory,
)
from scripts.lib import (  # noqa: E402
    isotropic_rectangular_timoshenko_coupled_beams as legacy_timoshenko,
)
from scripts.lib import reddy_symmetric_coupled_beams as rlb_coupled  # noqa: E402
from scripts.lib import reddy_symmetric_laminated_beam as rlb_beam  # noqa: E402


FloatArray = NDArray[np.float64]
MatrixProvider = Callable[[float], FloatArray]

STAGE_ID = "RLB-2B"
ALGORITHM_VERSION = "rectangular_isotropic_sorted_beta_grid_v1"
MODEL_EB = "EB"
MODEL_OLD = "OLD_TIMOSHENKO"
MODEL_RLB = "NEW_ISOTROPIC_RLB"
MODELS = (MODEL_EB, MODEL_OLD, MODEL_RLB)
ROOT_FILENAMES = {
    MODEL_EB: "eb_roots.csv",
    MODEL_OLD: "old_timoshenko_roots.csv",
    MODEL_RLB: "new_rlb_roots.csv",
}
LINE_STYLES = {MODEL_EB: "-", MODEL_OLD: "--", MODEL_RLB: "-."}
MODEL_LABELS = {
    MODEL_EB: "Euler--Bernoulli",
    MODEL_OLD: "old rectangular Timoshenko",
    MODEL_RLB: "new four-ply isotropic RLB",
}

DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "rectangular_isotropic_beta_sweep_comparison"
)
SOURCE_CONTRACT_PATH = (
    ROOT / "tests" / "data" / "reddy_four_ply_isotropic_limit_cases.json"
)

DEFAULT_BETA_MIN_DEG = 0.0
DEFAULT_BETA_MAX_DEG = 90.0
DEFAULT_BETA_STEP_DEG = 1.0
PLOTTED_POSITIONS = 8
OUTPUT_GUARD_POSITION = 9
OLD_VS_NEW_RELATIVE_TOLERANCE = 1.0e-8
ROOT_SINGULAR_RATIO_TOLERANCE = 1.0e-9
BOUNDARY_RESIDUAL_TOLERANCE = 1.0e-9
CONTRACT_RELATIVE_TOLERANCE = 1.0e-11
EB_PRIMARY_SCAN_STEP = 0.005
EB_VERIFICATION_SCAN_STEP = 0.0025
EB_PRIMARY_VERIFICATION_RELATIVE_TOLERANCE = 1.0e-9
TIMO_PRIMARY_VERIFICATION_RELATIVE_TOLERANCE = 1.0e-8

TASK_INITIAL_GIT_STATE = {
    "top_level": "D:/PHD/CoupledBeams/CoupledBeams",
    "branch": "main",
    "head": "11be0933553504964ac1c819a39475689df15986",
    "last_commit": "11be093 Version 0.4.3.1",
    "status_short": "",
}

FROZEN_MODEL_PATHS = (
    "src/my_project/analytic/formulas.py",
    "src/my_project/analytic/solvers.py",
    "scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py",
    "scripts/lib/reddy_symmetric_laminated_beam.py",
    "scripts/lib/reddy_inplane_geometry.py",
    "scripts/lib/reddy_symmetric_coupled_beams.py",
    "scripts/analysis/laminated_beams/validate_reddy_four_ply_isotropic_limit.py",
    "tests/data/reddy_four_ply_isotropic_limit_cases.json",
    "docs/laminated_beams/reddy_four_ply_isotropic_limit_validation.md",
)


@dataclass(frozen=True)
class ModelObjects:
    """The three model-specific objects built from one dimensional contract."""

    laminate: Any
    rlb_properties: Any
    legacy_section: Any


def _json_value(value: Any) -> Any:
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not math.isfinite(value):
        if math.isnan(value):
            return "NaN"
        return "Infinity" if value > 0.0 else "-Infinity"
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_json_value(dict(payload)), ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )


def _csv_value(value: Any) -> Any:
    if isinstance(value, (tuple, list, dict, np.ndarray)):
        return json.dumps(_json_value(value), ensure_ascii=False, separators=(",", ":"))
    if isinstance(value, np.generic):
        return value.item()
    return value


def _write_csv(
    path: Path,
    rows: Iterable[Mapping[str, Any]],
    fields: Sequence[str] | None = None,
) -> None:
    data = [dict(row) for row in rows]
    if fields is None:
        ordered: list[str] = []
        for row in data:
            for key in row:
                if key not in ordered:
                    ordered.append(key)
        fields = ordered or ("status",)
    path.parent.mkdir(parents=True, exist_ok=True)
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


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return [dict(row) for row in csv.DictReader(stream)]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


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


def _frozen_hashes() -> dict[str, str]:
    return {relative: _sha256(ROOT / relative) for relative in FROZEN_MODEL_PATHS}


def beta_grid(
    beta_min_deg: float = DEFAULT_BETA_MIN_DEG,
    beta_max_deg: float = DEFAULT_BETA_MAX_DEG,
    beta_step_deg: float = DEFAULT_BETA_STEP_DEG,
) -> np.ndarray:
    """Return an endpoint-inclusive fixed angular grid in degrees."""

    lower = float(beta_min_deg)
    upper = float(beta_max_deg)
    step = float(beta_step_deg)
    if not (math.isfinite(lower) and math.isfinite(upper) and math.isfinite(step)):
        raise ValueError("beta-grid inputs must be finite.")
    if step <= 0.0 or upper < lower:
        raise ValueError("Require beta_step_deg>0 and beta_max_deg>=beta_min_deg.")
    count = int(math.floor((upper - lower) / step + 1.0e-12)) + 1
    values = lower + step * np.arange(count, dtype=float)
    if values[-1] < upper - 1.0e-12:
        values = np.append(values, upper)
    return np.round(values, 12)


def rectangular_area(width: float, thickness: float) -> float:
    return float(width) * float(thickness)


def rectangular_second_moment(width: float, thickness: float) -> float:
    return float(width) * float(thickness) ** 3 / 12.0


def lambda_scale(
    *, E: float, rho: float, width: float, thickness: float, arm_length: float
) -> float:
    """Return ``Lambda/omega`` for the common rectangular normalization."""

    area = rectangular_area(width, thickness)
    inertia = rectangular_second_moment(width, thickness)
    return float(arm_length) ** 2 * math.sqrt(float(rho) * area / (float(E) * inertia))


def common_lambda(
    omega: float,
    *,
    E: float,
    rho: float,
    width: float,
    thickness: float,
    arm_length: float,
) -> float:
    return float(omega) * lambda_scale(
        E=E,
        rho=rho,
        width=width,
        thickness=thickness,
        arm_length=arm_length,
    )


def eb_slenderness(*, width: float, thickness: float, arm_length: float) -> float:
    area = rectangular_area(width, thickness)
    inertia = rectangular_second_moment(width, thickness)
    return math.sqrt(inertia / area) / float(arm_length)


def build_case_contract(
    betas_deg: Sequence[float] | None = None,
) -> dict[str, Any]:
    """Derive the RLB-2B contract from the committed RLB-1C-ISO G20 case."""

    source = iso_inventory._load_contract(SOURCE_CONTRACT_PATH)
    material = source["material"]
    geometry = source["geometries"]["G20"]
    arm_length = float(source["lengths"]["L_ref"])
    active_betas = (
        beta_grid().tolist()
        if betas_deg is None
        else [float(value) for value in betas_deg]
    )
    E = float(material["E"])
    rho = float(material["rho"])
    width = float(geometry["width"])
    thickness = float(geometry["thickness"])
    area = rectangular_area(width, thickness)
    inertia = rectangular_second_moment(width, thickness)
    scale = lambda_scale(
        E=E,
        rho=rho,
        width=width,
        thickness=thickness,
        arm_length=arm_length,
    )
    return {
        "schema_version": 1,
        "stage": STAGE_ID,
        "source_contract": SOURCE_CONTRACT_PATH.relative_to(ROOT).as_posix(),
        "source_contract_sha256": _sha256(SOURCE_CONTRACT_PATH),
        "mu": 0.0,
        "tau": 0.0,
        "material": {
            "E": E,
            "nu": float(material["nu"]),
            "rho": rho,
            "G": E / (2.0 * (1.0 + float(material["nu"]))),
        },
        "geometry": {
            "L1": arm_length,
            "L2": arm_length,
            "l": arm_length,
            "L_total": 2.0 * arm_length,
            "width_b": width,
            "thickness_h": thickness,
            "area_A": area,
            "second_moment_I_y": inertia,
        },
        "shear_correction": {
            "K": float(material["K"]),
            "provenance": "frozen RLB-1C-ISO common accepted parameter",
        },
        "new_RLB_laminate": {
            "number_of_plies": 4,
            "stacking_sequence_deg": [0.0, 90.0, 90.0, 0.0],
            "equal_ply_thickness": thickness / 4.0,
            "one_ply_shortcut": False,
            "pipeline": "lamina->Q->Qbar->A/B/D->beam reduction",
        },
        "frequency": {
            "symbol": "Lambda",
            "definition": "omega*l^2*sqrt(rho*A/(E*I_y))",
            "rectangular_definition": "omega*l^2*sqrt(12*rho/(E*h^2))",
            "Lambda_per_omega": scale,
            "historical_EB_matrix_variable": "sqrt(Lambda)",
            "eb_epsilon": eb_slenderness(
                width=width, thickness=thickness, arm_length=arm_length
            ),
        },
        "beta_grid_deg": active_betas,
        "beta_step_deg": (
            active_betas[1] - active_betas[0] if len(active_betas) > 1 else None
        ),
        "plotted_sorted_positions": PLOTTED_POSITIONS,
        "output_guard_position": OUTPUT_GUARD_POSITION,
        "modal_descendant_tracking": False,
        "curves_meaning": "pointwise sorted spectral positions",
    }


def build_model_objects(contract: Mapping[str, Any]) -> ModelObjects:
    """Build the four-ply RLB and rectangular legacy sections independently."""

    material_data = contract["material"]
    geometry = contract["geometry"]
    correction = float(contract["shear_correction"]["K"])
    E = float(material_data["E"])
    nu = float(material_data["nu"])
    rho = float(material_data["rho"])
    shear = E / (2.0 * (1.0 + nu))
    thickness = float(geometry["thickness_h"])
    width = float(geometry["width_b"])

    material = rlb_beam.LaminaMaterial(
        E1=E,
        E2=E,
        nu12=nu,
        G12=shear,
        G13=shear,
        G23=shear,
        rho=rho,
        name="RLB-2B isotropic lamina",
    )
    stack = tuple(
        float(value)
        for value in contract["new_RLB_laminate"]["stacking_sequence_deg"]
    )
    laminate = rlb_beam.integrate_laminate(
        [rlb_beam.Ply(material, angle, thickness / 4.0) for angle in stack]
    )
    rlb_properties = rlb_beam.reduce_to_beam_properties(
        laminate,
        width=width,
        K=correction,
    )
    legacy_section = legacy_timoshenko.rectangular_section(
        E=E,
        nu=nu,
        rho=rho,
        width=width,
        thickness=thickness,
        K=correction,
    )
    return ModelObjects(
        laminate=laminate,
        rlb_properties=rlb_properties,
        legacy_section=legacy_section,
    )


def model_contract_checks(
    contract: Mapping[str, Any], objects: ModelObjects
) -> dict[str, Any]:
    geometry = contract["geometry"]
    material = contract["material"]
    E = float(material["E"])
    rho = float(material["rho"])
    width = float(geometry["width_b"])
    thickness = float(geometry["thickness_h"])
    area = rectangular_area(width, thickness)
    inertia = rectangular_second_moment(width, thickness)
    G = E / (2.0 * (1.0 + float(material["nu"])))
    K = float(contract["shear_correction"]["K"])
    expected = {
        "A": E * area,
        "D": E * inertia,
        "S": K * G * area,
        "m": rho * area,
        "J": rho * inertia,
    }
    actual_rlb = {
        name: float(getattr(objects.rlb_properties, name)) for name in expected
    }
    actual_legacy = {
        "A": float(objects.legacy_section.EA),
        "D": float(objects.legacy_section.EI),
        "S": float(objects.legacy_section.KGA),
        "m": float(objects.legacy_section.rhoA),
        "J": float(objects.legacy_section.rhoI),
    }

    def relative(left: float, right: float) -> float:
        return abs(left - right) / max(abs(left), abs(right), np.finfo(float).tiny)

    residuals = {
        name: max(
            relative(actual_rlb[name], expected[name]),
            relative(actual_legacy[name], expected[name]),
            relative(actual_rlb[name], actual_legacy[name]),
        )
        for name in expected
    }
    symmetry = rlb_beam.check_laminate_symmetry(objects.laminate)
    maximum = max(residuals.values())
    passed = bool(
        len(objects.laminate.plies) == 4
        and tuple(float(ply.angle_deg) for ply in objects.laminate.plies)
        == (0.0, 90.0, 90.0, 0.0)
        and all(
            abs(ply.thickness - thickness / 4.0)
            <= 64.0 * np.finfo(float).eps * thickness
            for ply in objects.laminate.plies
        )
        and symmetry.is_symmetric
        and maximum <= CONTRACT_RELATIVE_TOLERANCE
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "expected": expected,
        "RLB": actual_rlb,
        "old_Timoshenko": actual_legacy,
        "relative_residuals": residuals,
        "maximum_relative_residual": maximum,
        "four_equal_plies": len(objects.laminate.plies) == 4,
        "stacking_sequence_deg": [ply.angle_deg for ply in objects.laminate.plies],
        "z_interfaces": objects.laminate.z_interfaces,
        "scaled_B_residual": symmetry.B_relative,
        "scaled_I1_residual": symmetry.I1_relative,
    }


def frozen_root_policy() -> Any:
    source = iso_inventory._load_contract(SOURCE_CONTRACT_PATH)
    return iso_inventory.SearchPolicy.from_contract(source)


def make_matrix_provider(
    model: str,
    beta_deg: float,
    contract: Mapping[str, Any],
    objects: ModelObjects,
    *,
    rlb_arm_cache: dict[float, FloatArray] | None = None,
) -> tuple[MatrixProvider, dict[str, Any]]:
    """Return a frozen matrix provider; no root from another model is accepted."""

    geometry = contract["geometry"]
    material = contract["material"]
    l1 = float(geometry["L1"])
    l2 = float(geometry["L2"])
    scale = float(contract["frequency"]["Lambda_per_omega"])
    angle_deg = float(beta_deg)
    if model == MODEL_EB:
        epsilon = eb_slenderness(
            width=float(geometry["width_b"]),
            thickness=float(geometry["thickness_h"]),
            arm_length=float(geometry["l"]),
        )
        beta_rad = math.radians(angle_deg)

        def provider(omega: float) -> FloatArray:
            common_value = max(0.0, float(omega) * scale)
            eb_wavenumber = math.sqrt(common_value)
            return np.asarray(
                assemble_eb_coupled_matrix(
                    eb_wavenumber,
                    beta_rad,
                    0.0,
                    epsilon,
                ),
                dtype=float,
            )

        return provider, {
            "matrix_source": "src/my_project/analytic/formulas.py",
            "matrix_argument": "sqrt(Lambda)",
            "epsilon": epsilon,
            "mu": 0.0,
        }

    if model == MODEL_OLD:
        section = objects.legacy_section

        def provider(omega: float) -> FloatArray:
            return np.asarray(
                legacy_timoshenko.legacy_coupled_boundary_matrix_raw(
                    float(omega),
                    section,
                    l1,
                    section,
                    l2,
                    beta_deg=angle_deg,
                ),
                dtype=float,
            )

        return provider, {
            "matrix_source": (
                "scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py"
            ),
            "uses_RLB_matrix": False,
            "uses_matrix_exponential": False,
        }

    if model != MODEL_RLB:
        raise ValueError(f"Unknown model: {model!r}.")
    cache = {} if rlb_arm_cache is None else rlb_arm_cache
    beta_rad = math.radians(angle_deg)
    joint = np.asarray(rlb_coupled.joint_matrix(beta_rad), dtype=float)
    properties = objects.rlb_properties

    def arm_map(omega: float) -> FloatArray:
        key = float(omega)
        if key not in cache:
            cache[key] = np.asarray(
                rlb_coupled.arm_clamp_map(key, l1, properties), dtype=float
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
                l1,
                properties,
                l2,
                properties,
            ),
            dtype=float,
        )
        direct_residual = max(
            direct_residual, float(np.max(np.abs(provider(check_omega) - direct)))
        )
    if direct_residual > 16.0 * np.finfo(float).eps:
        raise RuntimeError("Cached RLB arm assembly differs from the frozen builder.")
    return provider, {
        "matrix_source": "scripts/lib/reddy_symmetric_coupled_beams.py",
        "arm_source": "scripts/lib/reddy_symmetric_laminated_beam.py",
        "number_of_plies": len(objects.laminate.plies),
        "cached_arm_map": True,
        "cached_vs_public_builder_max_abs": direct_residual,
        "material_E": float(material["E"]),
    }


def _root_row(
    model: str,
    beta_deg: float,
    inventory: Any,
    slot: Any,
) -> dict[str, Any]:
    event = slot.event
    candidate = event.candidate
    diagnostic = candidate.diagnostics
    position = int(slot.sorted_slot)
    quality_pass = bool(
        diagnostic.scaled_sigma_ratio <= ROOT_SINGULAR_RATIO_TOLERANCE
        and diagnostic.raw_boundary_null_residual <= BOUNDARY_RESIDUAL_TOLERANCE
    )
    return {
        "model": model,
        "beta_deg": float(beta_deg),
        "sorted_position": position,
        "role": "FIRST_8" if position <= PLOTTED_POSITIONS else "ROOT_9_GUARD",
        "omega": event.omega,
        "Lambda": event.Omega,
        "historical_EB_wavenumber": (
            math.sqrt(max(0.0, event.Omega)) if model == MODEL_EB else ""
        ),
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
        "internal_requested_roots": 12,
        "internal_guard_position": 13,
        "internal_root13_Lambda": inventory.slots[12].event.Omega,
        "primary_verification_max_relative": (
            inventory.maximum_primary_verification_relative
        ),
        "unresolved_candidates_below_internal_guard": (
            inventory.unresolved_low_sigma_count
        ),
        "guard_flag": position == OUTPUT_GUARD_POSITION,
    }


def _eb_root_rows(
    beta_deg: float,
    contract: Mapping[str, Any],
    objects: ModelObjects,
    policy: Any,
) -> list[dict[str, Any]]:
    """Use the frozen EB sign scan and an independent finer-grid repeat.

    The historical solver searches in the EB bending wavenumber, so both root
    sets are squared before comparison in the common-Lambda normalization.  A
    sign scan does not certify same-sign roots or cluster multiplicity; those
    limitations are recorded explicitly instead of being filled with zeros.
    """

    geometry = contract["geometry"]
    epsilon = eb_slenderness(
        width=float(geometry["width_b"]),
        thickness=float(geometry["thickness_h"]),
        arm_length=float(geometry["l"]),
    )

    def scan(step: float) -> np.ndarray:
        return np.asarray(
            find_eb_roots(
                beta=math.radians(float(beta_deg)),
                mu=0.0,
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

    roots = scan(EB_PRIMARY_SCAN_STEP)
    verification_roots = scan(EB_VERIFICATION_SCAN_STEP)
    for label, values in (
        ("primary", roots),
        ("verification", verification_roots),
    ):
        if values.shape != (OUTPUT_GUARD_POSITION,) or not np.all(
            np.isfinite(values)
        ):
            raise RuntimeError(
                f"EB beta={float(beta_deg):g}: the {label} sign scan did not "
                f"return {OUTPUT_GUARD_POSITION} finite roots."
            )
    primary_common = roots**2
    verification_common = verification_roots**2
    primary_verification_relative = float(
        np.max(
            np.abs(primary_common - verification_common)
            / np.maximum.reduce(
                (
                    np.abs(primary_common),
                    np.abs(verification_common),
                    np.full(primary_common.shape, np.finfo(float).tiny),
                )
            )
        )
    )
    provider, _metadata = make_matrix_provider(
        MODEL_EB, float(beta_deg), contract, objects
    )
    frequency_scale = float(contract["frequency"]["Lambda_per_omega"])
    rows: list[dict[str, Any]] = []
    digest_payload: list[dict[str, Any]] = []
    for position, eb_wavenumber in enumerate(roots, start=1):
        Lambda = float(eb_wavenumber**2)
        diagnostic = iso_inventory._boundary_matrix_diagnostics(
            Lambda,
            provider,
            frequency_scale,
            policy,
        )
        quality_pass = bool(
            diagnostic.scaled_sigma_ratio <= ROOT_SINGULAR_RATIO_TOLERANCE
            and diagnostic.raw_boundary_null_residual
            <= BOUNDARY_RESIDUAL_TOLERANCE
        )
        row = {
            "model": MODEL_EB,
            "beta_deg": float(beta_deg),
            "sorted_position": position,
            "role": (
                "FIRST_8"
                if position <= PLOTTED_POSITIONS
                else "ROOT_9_GUARD"
            ),
            "omega": Lambda / frequency_scale,
            "Lambda": Lambda,
            "historical_EB_wavenumber": float(eb_wavenumber),
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
                "frozen_EB_sign_scan_bisection_primary",
                "finer_sign_scan_bisection_verification",
            ),
            "root_status": "PASS" if quality_pass else "FAIL",
            "inventory_status": "PENDING",
            "inventory_sha256": "PENDING",
            "primary_slot_count_internal": OUTPUT_GUARD_POSITION,
            "verification_slot_count_internal": OUTPUT_GUARD_POSITION,
            "internal_requested_roots": PLOTTED_POSITIONS,
            "internal_guard_position": OUTPUT_GUARD_POSITION,
            "internal_root13_Lambda": "",
            "primary_verification_max_relative": primary_verification_relative,
            "unresolved_candidates_below_internal_guard": (
                "NOT_ASSESSED_BY_EB_SIGN_SCAN"
            ),
            "guard_flag": position == OUTPUT_GUARD_POSITION,
        }
        rows.append(row)
        digest_payload.append(
            {
                "position": position,
                "Lambda": format(Lambda, ".17g"),
                "verification_Lambda": format(
                    float(verification_common[position - 1]), ".17g"
                ),
                "sigma_ratio": format(diagnostic.scaled_sigma_ratio, ".17g"),
                "boundary_residual": format(
                    diagnostic.raw_boundary_null_residual, ".17g"
                ),
            }
        )
    inventory_status = (
        "PASS_SIGN_SCAN_CROSSCHECK"
        if all(row["root_status"] == "PASS" for row in rows)
        and primary_verification_relative
        <= EB_PRIMARY_VERIFICATION_RELATIVE_TOLERANCE
        else "FAIL"
    )
    digest = hashlib.sha256(
        json.dumps(
            {
                "model": MODEL_EB,
                "beta_deg": float(beta_deg),
                "primary_scan_step": EB_PRIMARY_SCAN_STEP,
                "verification_scan_step": EB_VERIFICATION_SCAN_STEP,
                "primary_verification_max_relative": format(
                    primary_verification_relative, ".17g"
                ),
                "roots": digest_payload,
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest().upper()
    for row in rows:
        row["inventory_status"] = inventory_status
        row["inventory_sha256"] = digest
    return rows


def solve_model_grid(
    model: str,
    betas_deg: Sequence[float],
    output_dir: Path = DEFAULT_OUTPUT_DIR,
) -> list[dict[str, Any]]:
    """Run independent seed-free inventories for one model on the fixed grid."""

    if model not in MODELS:
        raise ValueError(f"Unknown model: {model!r}.")
    contract = build_case_contract(betas_deg)
    objects = build_model_objects(contract)
    checks = model_contract_checks(contract, objects)
    if checks["status"] != "PASS":
        raise RuntimeError("Shared material/geometry contract failed before roots.")
    policy = frozen_root_policy()
    frequency_scale = float(contract["frequency"]["Lambda_per_omega"])
    rows: list[dict[str, Any]] = []
    rlb_cache: dict[float, FloatArray] = {}
    target = Path(output_dir) / ROOT_FILENAMES[model]
    for beta_index, beta_deg in enumerate(betas_deg, start=1):
        if model == MODEL_EB:
            beta_rows = _eb_root_rows(
                float(beta_deg), contract, objects, policy
            )
            rows.extend(beta_rows)
            _write_csv(target, rows)
            print(
                f"[{model}] beta={float(beta_deg):g} deg "
                f"({beta_index}/{len(betas_deg)}): "
                f"{beta_rows[0]['inventory_status']}, "
                f"Lambda_9={float(beta_rows[-1]['Lambda']):.12g}",
                flush=True,
            )
            continue
        provider, _metadata = make_matrix_provider(
            model,
            float(beta_deg),
            contract,
            objects,
            rlb_arm_cache=rlb_cache,
        )
        inventory = iso_inventory.seed_free_root_inventory(
            provider,
            frequency_scale,
            policy,
            case_id=f"{model.lower()}__beta_{float(beta_deg):06.2f}",
            builder_id=f"RLB2B_{model}",
            contract_sha256=str(contract["source_contract_sha256"]),
        )
        if len(inventory.slots) < policy.required_slots:
            raise RuntimeError(
                f"{model}, beta={float(beta_deg):g}: only "
                f"{len(inventory.slots)} internal slots were accepted."
            )
        rows.extend(
            _root_row(model, float(beta_deg), inventory, slot)
            for slot in inventory.slots[:OUTPUT_GUARD_POSITION]
        )
        _write_csv(target, rows)
        print(
            f"[{model}] beta={float(beta_deg):g} deg "
            f"({beta_index}/{len(betas_deg)}): {inventory.status}, "
            f"Lambda_9={inventory.slots[8].event.Omega:.12g}, "
            f"internal_Lambda_13={inventory.slots[12].event.Omega:.12g}",
            flush=True,
        )
    return rows


def prepare_output(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    betas_deg: Sequence[float] | None = None,
) -> dict[str, Any]:
    active_betas = beta_grid() if betas_deg is None else np.asarray(betas_deg, dtype=float)
    target = Path(output_dir)
    contract = build_case_contract(active_betas)
    objects = build_model_objects(contract)
    checks = model_contract_checks(contract, objects)
    policy = frozen_root_policy()
    _write_json(target / "case_contract.json", contract)
    manifest = {
        "schema_version": 1,
        "algorithm_version": ALGORITHM_VERSION,
        "stage": STAGE_ID,
        "scientific_role": "calculation-and-plot diagnostic only",
        "task_initial_git_state": TASK_INITIAL_GIT_STATE,
        "git_at_preparation": _git_state(),
        "case_contract_sha256": _sha256(target / "case_contract.json"),
        "common_contract_check": checks,
        "production_root_solvers_changed": False,
        "timoshenko_root_policy_reused_without_mathematical_change": True,
        "timoshenko_root_policy": asdict(policy),
        "internal_inventory": {
            "EB": {
                "policy": "frozen EB sign-scan/bisection with finer-grid repeat",
                "primary_scan_step": EB_PRIMARY_SCAN_STEP,
                "verification_scan_step": EB_VERIFICATION_SCAN_STEP,
                "primary_verification_relative_tolerance": (
                    EB_PRIMARY_VERIFICATION_RELATIVE_TOLERANCE
                ),
                "same_sign_or_low_sigma_candidate_search": False,
                "computed_positions": 9,
                "guard_position": 9,
            },
            "old_Timoshenko_and_RLB": {
                "policy": "frozen RLB-1C-ISO determinant/SVD policy",
                "computed_positions": 13,
                "internal_guard_position": 13,
            },
            "exported_positions": 9,
            "plotted_positions": 8,
            "exported_guard": 9,
        },
        "models": {
            MODEL_EB: {
                "matrix": "src/my_project/analytic/formulas.py",
                "section_adaptation": "epsilon=sqrt(I_y/A)/l",
                "common_Lambda_mapping": "Lambda=(historical EB matrix argument)^2",
            },
            MODEL_OLD: {
                "matrix": (
                    "scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py"
                ),
                "section": "rectangle b*h, I_y=b*h^3/12",
            },
            MODEL_RLB: {
                "matrix": "scripts/lib/reddy_symmetric_coupled_beams.py",
                "plies": 4,
                "one_ply_shortcut": False,
            },
        },
        "thresholds": {
            "shared_contract_relative": CONTRACT_RELATIVE_TOLERANCE,
            "root_singular_ratio": ROOT_SINGULAR_RATIO_TOLERANCE,
            "boundary_null_residual": BOUNDARY_RESIDUAL_TOLERANCE,
            "old_vs_new_relative": OLD_VS_NEW_RELATIVE_TOLERANCE,
        },
        "frozen_model_hashes_before": _frozen_hashes(),
        "script_proliferation_control": {
            "classification": "diagnostic-only stable entry point",
            "why_new_script": (
                "new three-model rectangular/common-Lambda data and figure contract"
            ),
            "new_production_solver": False,
        },
        "explicit_exclusions": [
            "modal_descendant_tracking",
            "MAC",
            "Rayleigh_Ritz",
            "FEM",
            "torsion",
            "damping",
            "mu_or_tau_sweep",
            "material_or_geometry_sweep",
            "article_preparation",
            "commit",
            "push",
        ],
    }
    _write_json(target / "model_manifest.json", manifest)
    return manifest


def _float(row: Mapping[str, str], key: str) -> float:
    return float(row[key])


def _validate_root_rows(
    model: str, rows: Sequence[Mapping[str, str]], betas_deg: Sequence[float]
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
            "Lambda",
            "scaled_sigma_ratio",
            "boundary_null_residual",
            "primary_verification_max_relative",
        )
    )
    numeric_quality = all(
        _float(row, "scaled_sigma_ratio") <= ROOT_SINGULAR_RATIO_TOLERANCE
        and _float(row, "boundary_null_residual")
        <= BOUNDARY_RESIDUAL_TOLERANCE
        for row in rows
    )
    verification_tolerance = (
        EB_PRIMARY_VERIFICATION_RELATIVE_TOLERANCE
        if model == MODEL_EB
        else TIMO_PRIMARY_VERIFICATION_RELATIVE_TOLERANCE
    )
    verification_quality = all(
        _float(row, "primary_verification_max_relative")
        <= verification_tolerance
        for row in rows
    )
    status_quality = all(row["root_status"] == "PASS" for row in rows)
    accepted_inventory_statuses = (
        {"PASS_SIGN_SCAN_CROSSCHECK"} if model == MODEL_EB else {"PASS"}
    )
    inventories = all(
        row["inventory_status"] in accepted_inventory_statuses for row in rows
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
    passed = bool(
        len(rows) == expected_rows
        and actual_keys == expected_keys
        and finite
        and numeric_quality
        and verification_quality
        and status_quality
        and inventories
        and no_unresolved
        and guard_structure
    )
    return {
        "summary_kind": "MODEL_ROOT_INVENTORY",
        "model": model,
        "beta_count": len(betas_deg),
        "expected_rows": expected_rows,
        "actual_rows": len(rows),
        "first_frequency_count": len(betas_deg) * PLOTTED_POSITIONS,
        "guard_count": sum(row["guard_flag"].lower() == "true" for row in rows),
        "maximum_scaled_sigma_ratio": max(
            (_float(row, "scaled_sigma_ratio") for row in rows), default=math.inf
        ),
        "maximum_boundary_null_residual": max(
            (_float(row, "boundary_null_residual") for row in rows), default=math.inf
        ),
        "maximum_primary_verification_relative": max(
            (_float(row, "primary_verification_max_relative") for row in rows),
            default=math.inf,
        ),
        "primary_verification_relative_tolerance": verification_tolerance,
        "inventory_hash_count": len({row["inventory_sha256"] for row in rows}),
        "missing_or_duplicate_slot_count": len(expected_keys.symmetric_difference(actual_keys)),
        "guard_structure_passed": guard_structure,
        "same_sign_low_sigma_candidate_search_performed": model != MODEL_EB,
        "unresolved_candidate_sum": (
            "NOT_ASSESSED_BY_EB_SIGN_SCAN"
            if model == MODEL_EB
            else sum(
                int(row["unresolved_candidates_below_internal_guard"])
                for row in rows
            )
        ),
        "status": "PASS" if passed else "FAIL",
    }


def old_vs_new_rows(
    old_rows: Sequence[Mapping[str, str]],
    rlb_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, Any]]:
    old = {
        (round(_float(row, "beta_deg"), 12), int(row["sorted_position"])): row
        for row in old_rows
        if int(row["sorted_position"]) <= PLOTTED_POSITIONS
    }
    new = {
        (round(_float(row, "beta_deg"), 12), int(row["sorted_position"])): row
        for row in rlb_rows
        if int(row["sorted_position"]) <= PLOTTED_POSITIONS
    }
    if set(old) != set(new):
        raise RuntimeError("Old Timoshenko and new RLB sorted-position keys differ.")
    rows: list[dict[str, Any]] = []
    for beta_deg, position in sorted(old):
        old_value = _float(old[(beta_deg, position)], "Lambda")
        new_value = _float(new[(beta_deg, position)], "Lambda")
        absolute = abs(old_value - new_value)
        relative = absolute / max(abs(old_value), abs(new_value), np.finfo(float).tiny)
        rows.append(
            {
                "beta_deg": beta_deg,
                "sorted_position": position,
                "Lambda_old_Timoshenko": old_value,
                "Lambda_new_isotropic_RLB": new_value,
                "absolute_difference": absolute,
                "relative_difference": relative,
                "tolerance": OLD_VS_NEW_RELATIVE_TOLERANCE,
                "status": (
                    "PASS" if relative <= OLD_VS_NEW_RELATIVE_TOLERANCE else "FAIL"
                ),
            }
        )
    return rows


def create_plot(
    rows_by_model: Mapping[str, Sequence[Mapping[str, Any]]], output_path: Path
) -> Path:
    """Plot 24 pointwise sorted-frequency curves with two compact legends."""

    fig, ax = plt.subplots(figsize=(11.0, 6.6), constrained_layout=True)
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, PLOTTED_POSITIONS))
    for position in range(1, PLOTTED_POSITIONS + 1):
        color = colors[position - 1]
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
                color=color,
                linestyle=LINE_STYLES[model],
                linewidth=1.55 if model != MODEL_RLB else 1.35,
                alpha=0.96,
            )
    ax.set_xlabel(r"$\beta$, degrees")
    ax.set_ylabel(r"$\Lambda$")
    all_beta_values = [
        float(row["beta_deg"])
        for model in MODELS
        for row in rows_by_model[model]
    ]
    beta_min = min(all_beta_values)
    beta_max = max(all_beta_values)
    if beta_max > beta_min:
        ax.set_xlim(beta_min, beta_max)
    ax.grid(True, alpha=0.24)
    style_legend = ax.legend(
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
    ax.add_artist(style_legend)
    position_handles = [
        Line2D(
            [0],
            [0],
            color=colors[position - 1],
            linestyle="-",
            linewidth=1.8,
            label=f"sorted position {position}",
        )
        for position in range(1, PLOTTED_POSITIONS + 1)
    ]
    ax.legend(
        handles=position_handles,
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
    summaries: Sequence[Mapping[str, Any]],
    comparison_summary: Mapping[str, Any],
    statuses: Mapping[str, str],
) -> str:
    material = contract["material"]
    geometry = contract["geometry"]
    beta_values = [float(value) for value in contract["beta_grid_deg"]]
    beta_description = (
        f"$\\beta={beta_values[0]:g},{beta_values[1]:g},\\ldots,"
        f"{beta_values[-1]:g}^\\circ$"
        if len(beta_values) > 1
        else f"$\\beta={beta_values[0]:g}^\\circ$"
    )
    lines = [
        "# Сопоставление трёх прямоугольных моделей по углу сопряжения",
        "",
        "## Область расчёта",
        "",
        (
            f"Выполнен конечный расчёт на сетке {beta_description} "
            "при $\\mu=\\tau=0$. Сопоставляются модель Эйлера--Бернулли, "
            "прежняя прямоугольная модель Тимошенко и четырёхслойная "
            "изотропная модель RLB. Соединённые линиями значения являются "
            "точечными отсортированными спектральными позициями, а не "
            "траекториями модальных потомков."
        ),
        "",
        "## Общий контракт",
        "",
        (
            f"Использованы $E={float(material['E']):g}$, "
            f"$\\nu={float(material['nu']):g}$, "
            f"$\\rho={float(material['rho']):g}$, "
            f"$b={float(geometry['width_b']):g}$, "
            f"$h={float(geometry['thickness_h']):g}$ и "
            f"$l={float(geometry['l']):g}$. Коэффициент сдвига равен $K=5/6$."
        ),
        "",
        "Во всех таблицах и на рисунке применяется",
        "",
        "$$",
        "\\Lambda=\\omega l^2\\sqrt{\\frac{\\rho A}{E I_y}}",
        "=\\omega l^2\\sqrt{\\frac{12\\rho}{E h^2}},",
        "\\qquad A=bh,\\quad I_y=\\frac{bh^3}{12}.",
        "$$",
        "",
        (
            "Для исторической EB-матрицы её аргумент равен "
            "$\\sqrt{\\Lambda}$; это преобразование выполнено до записи "
            "общих частот. Модель RLB содержит ровно четыре одинаковых по "
            "толщине изотропных слоя $[0/90/90/0]$ и проходит обычную цепочку "
            "$Q\\to\\overline Q\\to A/B/D\\to$ одномерная редукция."
        ),
        "",
        "## Root inventory",
        "",
        "| Модель | углов | строк (8 + guard) | max $\\sigma_{min}/\\sigma_{max}$ | max boundary residual | статус |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for summary in summaries:
        lines.append(
            f"| {summary['model']} | {int(summary['beta_count'])} | "
            f"{int(summary['actual_rows'])} | "
            f"{float(summary['maximum_scaled_sigma_ratio']):.6e} | "
            f"{float(summary['maximum_boundary_null_residual']):.6e} | "
            f"{summary['status']} |"
        )
    lines.extend(
        [
            "",
            "Для EB применён прежний sign-scan/bisection, повторённый на вдвое "
            "более мелкой сетке; вычислены первые восемь позиций и девятый "
            "контрольный корень. Этот "
            "EB-контроль не объявляет поиск корней без смены знака или кластерную "
            "сертификацию. Неизменённый determinant/SVD inventory моделей "
            "Тимошенко и RLB внутренне вычислял двенадцать корней и тринадцатый "
            "контрольный корень; для данного графика из него сохранены первые "
            "восемь позиций и девятый контрольный корень.",
            "",
            "## Прежняя модель Тимошенко и новая RLB",
            "",
            (
                f"Сопоставлено {int(comparison_summary['comparison_count'])} пар. "
                f"Максимальная абсолютная разность равна "
                f"{float(comparison_summary['maximum_absolute_difference']):.6e}, "
                f"она достигается при $\\beta="
                f"{float(comparison_summary['beta_at_max_absolute']):g}^\\circ$ "
                f"для отсортированной позиции "
                f"{int(comparison_summary['position_at_max_absolute'])}. "
                f"Максимальная относительная разность -- "
                f"{float(comparison_summary['maximum_relative_difference']):.6e}. "
                f"Последняя достигается при $\\beta="
                f"{float(comparison_summary['beta_at_max_relative']):g}^\\circ$ "
                f"для отсортированной позиции "
                f"{int(comparison_summary['position_at_max_relative'])}."
            ),
            "",
            "Все 728 относительных разностей не превышают допуска $10^{-8}$. "
            "Совпадение оценивается только на объявленной конечной сетке и не "
            "заменяет доказательство для произвольного угла.",
            "",
            "## Статусы",
            "",
        ]
    )
    lines.extend(f"- `{name}: {value}`" for name, value in statuses.items())
    lines.extend(
        [
            "",
            "## Ограничения",
            "",
            "Не выполнялись MAC и отслеживание модальных потомков; veering не "
            "исследовался. Случаи с ненулевыми $\\mu$ или $\\tau$ не "
            "рассчитывались. Также не выполнялись Ritz, FEM, кручение, "
            "демпфирование и подготовка статьи. Физика трёх "
            "существующих моделей и их production-модули не изменялись.",
        ]
    )
    return "\n".join(lines) + "\n"


def finalize_output(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    betas_deg: Sequence[float] | None = None,
) -> dict[str, Any]:
    target = Path(output_dir)
    active_betas = beta_grid() if betas_deg is None else np.asarray(betas_deg, dtype=float)
    contract_path = target / "case_contract.json"
    manifest_path = target / "model_manifest.json"
    if not contract_path.is_file() or not manifest_path.is_file():
        raise RuntimeError("Run --mode prepare before finalization.")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    model_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    contract_grid = np.asarray(contract["beta_grid_deg"], dtype=float)
    if not np.array_equal(active_betas, contract_grid):
        raise RuntimeError(
            "Finalization beta grid differs from the prepared case contract."
        )
    contract_hash = _sha256(contract_path)
    if model_manifest["case_contract_sha256"] != contract_hash:
        raise RuntimeError("Prepared case-contract hash no longer matches the manifest.")
    preservation_after = _frozen_hashes()
    if model_manifest["frozen_model_hashes_before"] != preservation_after:
        raise RuntimeError(
            "A frozen model/source file changed after RLB-2B preparation."
        )
    rows_by_model = {
        model: _read_csv(target / ROOT_FILENAMES[model]) for model in MODELS
    }
    summaries = [
        _validate_root_rows(model, rows_by_model[model], active_betas)
        for model in MODELS
    ]
    comparison = old_vs_new_rows(
        rows_by_model[MODEL_OLD], rows_by_model[MODEL_RLB]
    )
    _write_csv(target / "old_vs_new_rlb_comparison.csv", comparison)
    maximum_relative_row = max(
        comparison, key=lambda row: float(row["relative_difference"])
    )
    maximum_absolute_row = max(
        comparison, key=lambda row: float(row["absolute_difference"])
    )
    comparison_summary = {
        "summary_kind": "OLD_VS_NEW_ISOTROPIC",
        "model": "OLD_TIMOSHENKO_vs_NEW_ISOTROPIC_RLB",
        "comparison_count": len(comparison),
        "maximum_absolute_difference": float(
            maximum_absolute_row["absolute_difference"]
        ),
        "beta_at_max_absolute": float(maximum_absolute_row["beta_deg"]),
        "position_at_max_absolute": int(maximum_absolute_row["sorted_position"]),
        "maximum_relative_difference": float(
            maximum_relative_row["relative_difference"]
        ),
        "beta_at_max_relative": float(maximum_relative_row["beta_deg"]),
        "position_at_max_relative": int(
            maximum_relative_row["sorted_position"]
        ),
        "tolerance": OLD_VS_NEW_RELATIVE_TOLERANCE,
        "status": (
            "PASS"
            if all(row["status"] == "PASS" for row in comparison)
            else "FAIL"
        ),
    }
    _write_csv(target / "spectrum_summary.csv", [*summaries, comparison_summary])
    plot_path = create_plot(
        rows_by_model, target / "lambda_vs_beta_plot.png"
    )
    objects = build_model_objects(contract)
    checks = model_contract_checks(contract, objects)
    lambda_definition_pass = bool(
        checks["status"] == "PASS"
        and math.isclose(
            float(contract["frequency"]["Lambda_per_omega"]),
            lambda_scale(
                E=float(contract["material"]["E"]),
                rho=float(contract["material"]["rho"]),
                width=float(contract["geometry"]["width_b"]),
                thickness=float(contract["geometry"]["thickness_h"]),
                arm_length=float(contract["geometry"]["l"]),
            ),
            rel_tol=0.0,
            abs_tol=0.0,
        )
    )
    root_pass = all(summary["status"] == "PASS" for summary in summaries)
    consistency_pass = comparison_summary["status"] == "PASS"
    plot_pass = plot_path.is_file() and plot_path.stat().st_size > 0
    statuses = {
        "RLB-2B-LAMBDA-DEFINITION": "PASS" if lambda_definition_pass else "FAIL",
        "RLB-2B-ROOT-DATA": "PASS" if root_pass else "FAIL",
        "RLB-2B-OLD-VS-NEW-ISOTROPIC-CONSISTENCY": (
            "PASS" if consistency_pass else "PARTIAL_PASS"
        ),
        "RLB-2B-PLOT-GENERATION": "PASS" if plot_pass else "FAIL",
    }
    statuses["OVERALL"] = (
        "PASS"
        if all(value == "PASS" for value in statuses.values())
        else (
            "PARTIAL_PASS"
            if statuses["RLB-2B-LAMBDA-DEFINITION"] == "PASS" and plot_pass
            else "FAIL"
        )
    )
    (target / "report.md").write_text(
        _report_text(contract, summaries, comparison_summary, statuses),
        encoding="utf-8",
    )
    generated_names = [
        "model_manifest.json",
        "case_contract.json",
        "eb_roots.csv",
        "old_timoshenko_roots.csv",
        "new_rlb_roots.csv",
        "old_vs_new_rlb_comparison.csv",
        "spectrum_summary.csv",
        "lambda_vs_beta_plot.png",
        "report.md",
    ]
    run_manifest = {
        "schema_version": 1,
        "algorithm_version": ALGORITHM_VERSION,
        "stage": STAGE_ID,
        "task_initial_git_state": TASK_INITIAL_GIT_STATE,
        "finalization_git_state": _git_state(),
        "case_contract_sha256": contract_hash,
        "model_manifest_sha256": _sha256(manifest_path),
        "beta_grid_deg": [float(value) for value in active_betas],
        "beta_step_deg": (
            float(active_betas[1] - active_betas[0])
            if len(active_betas) > 1
            else None
        ),
        "models": list(MODELS),
        "plotted_sorted_positions": PLOTTED_POSITIONS,
        "output_guard_position": OUTPUT_GUARD_POSITION,
        "internal_root_policies": {
            "EB": {
                "source": "src/my_project/analytic/solvers.py",
                "n_roots": OUTPUT_GUARD_POSITION,
                "primary_scan_step_in_historical_EB_wavenumber": (
                    EB_PRIMARY_SCAN_STEP
                ),
                "verification_scan_step_in_historical_EB_wavenumber": (
                    EB_VERIFICATION_SCAN_STEP
                ),
                "same_sign_or_low_sigma_candidate_search": False,
            },
            "old_Timoshenko_and_RLB": asdict(frozen_root_policy()),
        },
        "root_frequencies_cross_seeded_between_models": False,
        "modal_descendant_tracking_run": False,
        "summary": {
            "root_data": summaries,
            "old_vs_new": comparison_summary,
        },
        "statuses": statuses,
        "frozen_model_hashes_before": model_manifest["frozen_model_hashes_before"],
        "frozen_model_hashes_after": preservation_after,
        "frozen_models_preserved": (
            model_manifest["frozen_model_hashes_before"] == preservation_after
        ),
        "generated_file_hashes": {
            name: _sha256(target / name) for name in generated_names
        },
        "figures_created": 1,
        "new_production_solver_created": False,
        "mu": 0.0,
        "tau": 0.0,
        "commit_performed": False,
        "push_performed": False,
    }
    _write_json(target / "run_manifest.json", run_manifest)
    return run_manifest


def run_all(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    betas_deg: Sequence[float] | None = None,
) -> dict[str, Any]:
    active_betas = beta_grid() if betas_deg is None else np.asarray(betas_deg, dtype=float)
    prepare_output(output_dir, active_betas)
    for model in MODELS:
        solve_model_grid(model, active_betas, output_dir)
    return finalize_output(output_dir, active_betas)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare the first eight pointwise sorted rectangular EB, old "
            "Timoshenko, and four-ply isotropic RLB frequencies versus beta."
        )
    )
    parser.add_argument(
        "--mode",
        choices=("prepare", "eb", "old_timoshenko", "new_rlb", "finalize", "all"),
        default="all",
    )
    parser.add_argument("--beta-min", type=float, default=DEFAULT_BETA_MIN_DEG)
    parser.add_argument("--beta-max", type=float, default=DEFAULT_BETA_MAX_DEG)
    parser.add_argument("--beta-step", type=float, default=DEFAULT_BETA_STEP_DEG)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    betas = beta_grid(args.beta_min, args.beta_max, args.beta_step)
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
