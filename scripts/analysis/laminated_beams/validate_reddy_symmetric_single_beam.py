"""Validate one straight symmetric laminated beam against Reddy Chapter 4.

The command keeps source transcription separate from scientific computation.
In particular, ``--source-check-only`` uses only the standard library and never
imports the numerical beam module.  The direct FSDT solver is based on the
physical first-order state; printed Eq. (4.3.50a) is not used to assemble it.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, replace
import hashlib
import importlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Callable, Iterable, Sequence


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SOURCE_CONTRACT_PATH = ROOT / "tests" / "data" / "reddy_ch4_table_4_3_3.json"
DEFAULT_OUTPUT_DIR = (
    ROOT / "results" / "laminated_beams" / "reddy_symmetric_single_beam"
)
REDDY_PDF = (
    ROOT
    / "docs"
    / "literature"
    / "pdf"
    / "EB__Mechanics_of_Laminated_Composite_Plates_and_Shells_-JN_Reddy.pdf"
)
ELISEEV_PDF = (
    ROOT
    / "docs"
    / "literature"
    / "pdf"
    / "ELISEEV V.V._ MEXANIKA UPRUGIX TEL, 1999, 341s.pdf"
)

K_SOURCE = 5.0 / 6.0
K_PROVENANCE = "INFERRED_FROM_INTERNAL_NUMERICAL_CONSISTENCY"
SOURCE_STATUS = "PASS_WITH_DOCUMENTED_SOURCE_RECONSTRUCTION"
NUMERICAL_ALLOWANCE = 5.0e-9

THRESHOLDS: dict[str, float] = {
    "scaled_matrix_symmetry": 1.0e-12,
    "symmetric_B_I1_relative": 1.0e-12,
    "compliance_vs_schur_relative": 1.0e-11,
    "analytic_vs_transfer_relative": 1.0e-10,
    "root_singular_ratio": 1.0e-9,
    "boundary_residual": 1.0e-9,
    "energy_identity_relative": 1.0e-8,
    "combined_union_relative": 1.0e-9,
    "one_minus_MAC": 1.0e-8,
    "normalization_401_801_relative": 1.0e-8,
}

def _reject_json_constant(value: str) -> None:
    raise ValueError(f"Non-finite JSON constant is forbidden: {value}")


def _load_source_contract() -> dict[str, Any]:
    with SOURCE_CONTRACT_PATH.open("r", encoding="utf-8") as stream:
        return json.load(stream, parse_constant=_reject_json_constant)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Refusing to write an empty CSV: {path}")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def _with_benchmark_metadata(
    rows: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Stamp generated tabular output with the frozen benchmark input."""

    return [
        {
            **row,
            "K": row.get("K", K_SOURCE),
            "K_provenance": row.get("K_provenance", K_PROVENANCE),
        }
        for row in rows
    ]


def _source_record_key(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        record["laminate_id"],
        record["row_role"],
        record["boundary_condition"],
        record["a_over_h"],
    )


def validate_source_contract(output_dir: Path) -> dict[str, Any]:
    """Validate the frozen source data without importing numerical code."""

    scientific_name = "scripts.lib.reddy_symmetric_laminated_beam"
    if scientific_name in sys.modules:
        raise RuntimeError("Scientific solver was imported before source-only validation.")
    contract = _load_source_contract()
    records = contract["records"]
    keys = [_source_record_key(record) for record in records]
    published = [record for record in records if record["source_status"] == "PUBLISHED"]
    missing = [
        record
        for record in records
        if record["source_status"] == "NOT_PUBLISHED_IN_TABLE"
    ]
    included = [record for record in records if record["include_in_source_pass_fail"]]
    assertions = {
        "schema_version": contract["schema_version"] == "1.0.0",
        "source_status": contract["source_status"] == SOURCE_STATUS,
        "record_count_180": len(records) == 180,
        "unique_records": len(keys) == len(set(keys)),
        "published_count_144": len(published) == 144,
        "missing_count_36": len(missing) == 36,
        "included_frequency_count_108": len(included) == 108,
        "missing_values_are_null": all(record["source_value"] is None for record in missing),
        "missing_excluded": all(not record["include_in_source_pass_fail"] for record in missing),
        "cross_ply_printed_label": (
            contract["approved_source_policy"]["cross_ply"]["printed_label"]
            == "(90/0)_s"
        ),
        "cross_ply_used_label": (
            contract["approved_source_policy"]["cross_ply"]["used_label"]
            == "(0/90)_s"
        ),
        "cross_ply_correction_status": (
            contract["approved_source_policy"]["cross_ply"]["correction_status"]
            == "CORRECTED_BY_INTERNAL_SOURCE_CROSSCHECK"
        ),
        "K_value": math.isclose(
            contract["approved_source_policy"]["shear_correction_factor"]["value"],
            K_SOURCE,
            rel_tol=0.0,
            abs_tol=2.0e-16,
        ),
        "K_provenance": (
            contract["approved_source_policy"]["shear_correction_factor"]["provenance"]
            == K_PROVENANCE
        ),
        "eq_4_3_50a_status": (
            contract["approved_source_policy"]["equation_4_3_50a"]["status"]
            == "PRINTED_INTERMEDIATE_FORMULA_INCONSISTENCY"
        ),
        "material_definition": contract["material_definition"]
        == {
            "source_formula": "Eq. (4.2.25)",
            "E1_over_E2": 25.0,
            "G12_over_E2": 0.5,
            "G13_over_E2": 0.5,
            "G23_over_E2": 0.2,
            "nu12": 0.25,
            "dimensionless_realization": {
                "E2": 1.0,
                "rho": 1.0,
                "thickness": 1.0,
                "width": 1.0,
            },
            "realization_policy": "Absolute E2 and rho are not source-fitted; the normalized frequency is invariant under consistent positive E2 and rho scaling.",
        },
        "frequency_normalization": contract["normalization"]["angular_frequency"]
        == "omega_bar = omega_1*length^2*sqrt(I0/(E2*thickness^3))",
        "state_ordering": contract["coordinate_and_ordering"]["combined_state"]
        == ["u", "w", "psi_b", "N", "Q", "M"],
        "Reddy_vs_project_z_mapping": (
            contract["coordinate_and_ordering"]["coordinates"]["z_reddy"]
            == "thickness coordinate, positive downward"
            and contract["coordinate_and_ordering"]["project_storage_convention"]
            ["coordinate"]
            == "z_project = -z_reddy"
            and contract["coordinate_and_ordering"]["ply_stack_direction"]
            == "project storage: bottom-to-top in increasing z_project"
        ),
        "literal_zero_ninety_printed_labels": (
            all(
                record["printed_label"] == "0"
                for record in records
                if record["laminate_id"] == "0_deg"
            )
            and all(
                record["printed_label"] == "90"
                for record in records
                if record["laminate_id"] == "90_deg"
            )
        ),
        "source_classical_HH_root": math.isclose(
            contract["source_classical_characteristic_roots"]["HH"]["value"],
            math.pi,
            rel_tol=0.0,
            abs_tol=2.0e-16,
        ),
        "source_classical_CC_root": (
            contract["source_classical_characteristic_roots"]["CC"]
            == {"expression": "printed 4.730", "value": 4.73}
        ),
        "source_classical_CF_root": (
            contract["source_classical_characteristic_roots"]["CF"]
            == {"expression": "printed 1.875", "value": 1.875}
        ),
    }
    if not all(assertions.values()):
        failed = [name for name, value in assertions.items() if not value]
        raise AssertionError(f"Source contract checks failed: {failed}")
    expected_reddy = contract["source"]["sha256"].upper()
    expected_eliseev = contract["independent_source"]["sha256"].upper()
    actual_reddy = _sha256(REDDY_PDF)
    actual_eliseev = _sha256(ELISEEV_PDF)
    if actual_reddy != expected_reddy or actual_eliseev != expected_eliseev:
        raise AssertionError("A local source PDF hash differs from the frozen contract.")
    manifest = {
        "status": SOURCE_STATUS,
        "source_contract": str(SOURCE_CONTRACT_PATH.relative_to(ROOT)).replace("\\", "/"),
        "source_contract_sha256": _sha256(SOURCE_CONTRACT_PATH),
        "scientific_solver_imported": scientific_name in sys.modules,
        "assertions": assertions,
        "records": {
            "total": len(records),
            "published": len(published),
            "included_frequency_comparisons": len(included),
            "not_published": len(missing),
        },
        "K": K_SOURCE,
        "K_provenance": K_PROVENANCE,
        "pdfs": [
            {
                "path": str(REDDY_PDF.relative_to(ROOT)).replace("\\", "/"),
                "sha256": actual_reddy,
                "printed_pages": [96, 97, 101, 113, 122, 128, 134, 135, 136, 138, 139, 142, 148, 168, 169, 176, 184, 186, 187, 188, 197, 198, 199, 200, 201],
                "pdf_indices_zero_based": [118, 119, 123, 135, 144, 150, 156, 157, 158, 160, 161, 164, 170, 190, 191, 198, 206, 208, 209, 210, 219, 220, 221, 222, 223],
            },
            {
                "path": str(ELISEEV_PDF.relative_to(ROOT)).replace("\\", "/"),
                "sha256": actual_eliseev,
                "printed_pages": [141, 142, 143, 151, 155, 162, 164, 166, 240, 242, 243],
                "pdf_indices_zero_based": [140, 141, 142, 150, 154, 161, 163, 165, 239, 241, 242],
            },
        ],
    }
    benchmark_contract = {
        "frozen_before_production_calculation": True,
        "source_status": SOURCE_STATUS,
        "dimensionless_realization": {
            "E2": 1.0,
            "rho": 1.0,
            "thickness": 1.0,
            "width": 1.0,
            "frequency_definition": "omega_bar = omega*length^2*sqrt(I0/(E2*thickness^3))",
        },
        "K": K_SOURCE,
        "K_provenance": K_PROVENANCE,
        "laminates": contract["laminate_definitions"],
        "a_over_h": contract["a_over_h_order"],
        "boundary_conditions": contract["boundary_condition_definitions"],
        "row_roles": contract["row_role_definitions"],
        "source_classical_characteristic_roots": contract[
            "source_classical_characteristic_roots"
        ],
        "thresholds": THRESHOLDS,
        "printed_tolerance": 5.0e-4,
        "printed_numerical_allowance": NUMERICAL_ALLOWANCE,
    }
    _write_json(output_dir / "source_manifest.json", manifest)
    _write_json(output_dir / "benchmark_contract.json", benchmark_contract)
    if scientific_name in sys.modules:
        raise RuntimeError("Source-only validation imported the scientific solver.")
    return {"contract": contract, "manifest": manifest, "benchmark": benchmark_contract}


def _scientific_module() -> Any:
    return importlib.import_module("scripts.lib.reddy_symmetric_laminated_beam")


def _source_material(core: Any, *, E2: float = 1.0, rho: float = 1.0) -> Any:
    return core.LaminaMaterial(
        E1=25.0 * E2,
        E2=E2,
        G12=0.5 * E2,
        G13=0.5 * E2,
        G23=0.2 * E2,
        nu12=0.25,
        rho=rho,
        name="Reddy Eq. (4.2.25) dimensionless realization",
    )


def _laminate(core: Any, angles: Sequence[float], *, thickness: float = 1.0, material: Any | None = None) -> Any:
    mat = material if material is not None else _source_material(core)
    ply_thickness = thickness / len(angles)
    return core.integrate_laminate(
        [core.Ply(mat, angle, ply_thickness) for angle in angles]
    )


def _dimensionless_frequency(
    omega: float,
    *,
    length: float,
    laminate: Any,
    E2: float = 1.0,
) -> float:
    return float(
        omega
        * length**2
        * math.sqrt(laminate.I0 / (E2 * laminate.thickness**3))
    )


def _relative_difference(value: float, reference: float) -> float:
    return abs(value - reference) / max(abs(reference), sys.float_info.min)


def _matrix_condition_and_residual(matrix: Any, index: int, np: Any) -> tuple[float, float]:
    unit = np.zeros(matrix.shape[0], dtype=float)
    unit[index] = 1.0
    solution = np.linalg.solve(matrix, unit)
    residual = np.linalg.norm(matrix @ solution - unit) / np.linalg.norm(unit)
    return float(np.linalg.cond(matrix)), float(residual)


def _root_near(
    target: float,
    function: Callable[[float], float],
    optimize: Any,
    *,
    relative_window: float = 0.03,
) -> float:
    """Find the nearest sign-changing root in a fixed window about a reference."""

    lower = max(target * (1.0 - relative_window), target * 0.5)
    upper = target * (1.0 + relative_window)
    samples = 801
    step = (upper - lower) / (samples - 1)
    points = [lower + index * step for index in range(samples)]
    values = [function(point) for point in points]
    brackets: list[tuple[float, float]] = []
    for index in range(samples - 1):
        left_value = values[index]
        right_value = values[index + 1]
        if not (math.isfinite(left_value) and math.isfinite(right_value)):
            continue
        if left_value == 0.0:
            return points[index]
        if left_value * right_value < 0.0:
            brackets.append((points[index], points[index + 1]))
    if not brackets:
        raise RuntimeError(f"No sign-changing root found near {target:.17g}.")
    bracket = min(brackets, key=lambda pair: abs(0.5 * (pair[0] + pair[1]) - target))
    return float(
        optimize.brentq(
            function,
            bracket[0],
            bracket[1],
            xtol=max(1.0e-15, abs(target) * 2.0e-14),
            rtol=4.0 * sys.float_info.epsilon,
        )
    )


def _bending_mac(x: Any, left: Any, right: Any, properties: Any, np: Any) -> float:
    density_cross = properties.m * left[:, 0] * right[:, 0] + properties.J * left[:, 1] * right[:, 1]
    density_left = properties.m * left[:, 0] ** 2 + properties.J * left[:, 1] ** 2
    density_right = properties.m * right[:, 0] ** 2 + properties.J * right[:, 1] ** 2
    cross = float(np.trapezoid(density_cross, x))
    norm_left = float(np.trapezoid(density_left, x))
    norm_right = float(np.trapezoid(density_right, x))
    return cross * cross / (norm_left * norm_right)


def _find_direct_fundamental(
    core: Any,
    properties: Any,
    laminate: Any,
    length: float,
    boundary_condition: str,
) -> tuple[float, Any]:
    if boundary_condition == "HH":
        omega = core.bending_dispersion_branches(properties, math.pi / length)[0].omega
        diagnostic = core.bending_root_diagnostic(
            omega,
            length,
            properties,
            "HH",
            sigma_ratio_tolerance=THRESHOLDS["root_singular_ratio"],
            detection="exact_HH_dispersion",
        )
        return omega, diagnostic
    scale = length**2 * math.sqrt(laminate.I0 / laminate.thickness**3)
    roots = core.find_bending_roots(
        properties,
        length,
        boundary_condition,
        omega_max=80.0 / scale,
        n_roots=1,
        omega_min=0.0,
        scan_points=1601,
        sigma_ratio_tolerance=THRESHOLDS["root_singular_ratio"],
        root_xtol=max(1.0e-15, 80.0 / scale * 2.0e-13),
        dedup_rtol=1.0e-7,
    )
    if not roots:
        raise RuntimeError(f"No accepted {boundary_condition} root for length={length}.")
    return roots[0].omega, roots[0]


def _classical_characteristic_frequency(
    core: Any,
    properties: Any,
    length: float,
    boundary_condition: str,
    source_roots: dict[str, Any],
) -> float:
    wavenumber = source_roots[boundary_condition]["value"] / length
    return core.bending_dispersion_branches(properties, wavenumber)[0].omega


def _source_comparisons(
    core: Any,
    contract: dict[str, Any],
    laminate_data: dict[str, tuple[Any, Any]],
    np: Any,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    records = contract["records"]
    direct_cache: dict[tuple[str, str, int, bool], tuple[float, Any]] = {}
    classical_cache: dict[tuple[str, str, int, bool], float] = {}
    comparison_rows: list[dict[str, Any]] = []
    root_rows: list[dict[str, Any]] = []
    for record in records:
        if record["row_role"] == "buckling_load":
            comparison_rows.append(
                {
                    **record,
                    "computed_value": "",
                    "absolute_error": "",
                    "relative_error": "",
                    "status": "TRANSCRIPTION_ONLY",
                }
            )
            continue
        laminate, base_properties = laminate_data[record["laminate_id"]]
        with_ri = record["row_role"].endswith("with_RI") and not record["row_role"].endswith("without_RI")
        properties = base_properties if with_ri else core.without_rotary_inertia(base_properties)
        length = float(record["a_over_h"] * laminate.thickness)
        cache_key = (
            record["laminate_id"],
            record["boundary_condition"],
            record["a_over_h"],
            with_ri,
        )
        diagnostic = None
        if record["row_role"].startswith("fsdt_frequency"):
            if cache_key not in direct_cache:
                direct_cache[cache_key] = _find_direct_fundamental(
                    core,
                    properties,
                    laminate,
                    length,
                    record["boundary_condition"],
                )
            omega, diagnostic = direct_cache[cache_key]
        else:
            if cache_key not in classical_cache:
                classical_cache[cache_key] = _classical_characteristic_frequency(
                    core,
                    properties,
                    length,
                    record["boundary_condition"],
                    contract["source_classical_characteristic_roots"],
                )
            omega = classical_cache[cache_key]
        computed = _dimensionless_frequency(
            omega,
            length=length,
            laminate=laminate,
        )
        source_value = record["source_value"]
        if source_value is None:
            absolute_error: float | str = ""
            relative_error: float | str = ""
            status = "NOT_PUBLISHED_IN_TABLE"
        else:
            absolute_error = abs(computed - source_value)
            relative_error = _relative_difference(computed, source_value)
            allowed = record["printed_tolerance"] + NUMERICAL_ALLOWANCE
            status = "PASS" if absolute_error <= allowed else "FAIL"
        row = {
            **record,
            "computed_value": computed,
            "absolute_error": absolute_error,
            "relative_error": relative_error,
            "status": status,
            "root_omega": omega,
            "root_status": "ACCEPTED" if diagnostic is None or diagnostic.accepted else "REJECTED",
            "scaled_determinant": "" if diagnostic is None else diagnostic.determinant,
            "sigma_min": "" if diagnostic is None else diagnostic.sigma_min,
            "sigma_max": "" if diagnostic is None else diagnostic.sigma_max,
            "relative_singular_residual": "" if diagnostic is None else diagnostic.sigma_ratio,
            "condition_number": "" if diagnostic is None else diagnostic.condition_number,
        }
        comparison_rows.append(row)
        if diagnostic is not None:
            raw = core.bending_boundary_matrix(
                omega,
                length,
                properties,
                record["boundary_condition"],
            )
            root_rows.append(
                {
                    "stage": "source_table",
                    "laminate_id": record["laminate_id"],
                    "boundary_condition": record["boundary_condition"],
                    "a_over_h": record["a_over_h"],
                    "rotary_inertia": with_ri,
                    "omega": omega,
                    "omega_bar": computed,
                    "raw_determinant": float(np.linalg.det(raw)),
                    **asdict(diagnostic),
                }
            )
    summaries: dict[str, Any] = {}
    for tier in ("A", "B"):
        compared = [
            row
            for row in comparison_rows
            if row["benchmark_tier"] == tier and row["include_in_source_pass_fail"]
        ]
        failures = [row for row in compared if row["status"] != "PASS"]
        summaries[tier] = {
            "comparison_count": len(compared),
            "pass_count": len(compared) - len(failures),
            "fail_count": len(failures),
            "maximum_absolute_error": max(float(row["absolute_error"]) for row in compared),
            "failures": [
                {
                    "laminate_id": row["laminate_id"],
                    "boundary_condition": row["boundary_condition"],
                    "a_over_h": row["a_over_h"],
                    "row_role": row["row_role"],
                    "source_value": row["source_value"],
                    "computed_value": row["computed_value"],
                    "absolute_error": row["absolute_error"],
                }
                for row in failures
            ],
        }
    return comparison_rows, root_rows, summaries


def _laminate_validation(core: Any, np: Any) -> tuple[dict[str, tuple[Any, Any]], list[dict[str, Any]], dict[str, Any]]:
    definitions = {
        "0_deg": [0.0],
        "90_deg": [90.0],
        "cross_ply_0_90_s": [0.0, 90.0, 90.0, 0.0],
        "angle_ply_45_minus45_s": [45.0, -45.0, -45.0, 45.0],
    }
    data: dict[str, tuple[Any, Any]] = {}
    rows: list[dict[str, Any]] = []
    max_symmetry = 0.0
    max_reduction = 0.0
    max_solve_residual = 0.0
    for laminate_id, angles in definitions.items():
        laminate = _laminate(core, angles)
        properties = core.reduce_to_beam_properties(
            laminate,
            width=1.0,
            K=K_SOURCE,
            symmetry_tolerance=THRESHOLDS["symmetric_B_I1_relative"],
            reduction_tolerance=THRESHOLDS["compliance_vs_schur_relative"],
        )
        symmetry = core.check_laminate_symmetry(
            laminate,
            tolerance=THRESHOLDS["symmetric_B_I1_relative"],
        )
        cond_A, residual_A = _matrix_condition_and_residual(laminate.A, 0, np)
        cond_D, residual_D = _matrix_condition_and_residual(laminate.D, 0, np)
        cond_S, residual_S = _matrix_condition_and_residual(laminate.shear, 1, np)
        max_symmetry = max(max_symmetry, symmetry.B_relative, symmetry.I1_relative)
        max_reduction = max(
            max_reduction,
            properties.axial_reduction.relative_difference,
            properties.bending_reduction.relative_difference,
            properties.shear_reduction_before_K.relative_difference,
        )
        max_solve_residual = max(max_solve_residual, residual_A, residual_D, residual_S)
        row: dict[str, Any] = {
            "laminate_id": laminate_id,
            "stack_bottom_to_top_deg": "/".join(f"{angle:g}" for angle in angles),
            "ply_count": len(angles),
            "thickness": laminate.thickness,
            "z_interfaces": json.dumps(laminate.z_interfaces.tolist()),
            "K": properties.K,
            "K_provenance": K_PROVENANCE,
            "A_axial": properties.A,
            "D_bending": properties.D,
            "S_shear": properties.S,
            "m": properties.m,
            "J": properties.J,
            "B_relative": symmetry.B_relative,
            "I1_relative": symmetry.I1_relative,
            "A_condition": cond_A,
            "D_condition": cond_D,
            "shear_condition": cond_S,
            "A_solve_residual": residual_A,
            "D_solve_residual": residual_D,
            "shear_solve_residual": residual_S,
            "axial_compliance_schur_relative": properties.axial_reduction.relative_difference,
            "bending_compliance_schur_relative": properties.bending_reduction.relative_difference,
            "shear_compliance_schur_relative": properties.shear_reduction_before_K.relative_difference,
            "A_min_eigenvalue": float(np.linalg.eigvalsh(laminate.A)[0]),
            "D_min_eigenvalue": float(np.linalg.eigvalsh(laminate.D)[0]),
            "shear_min_eigenvalue": float(np.linalg.eigvalsh(laminate.shear)[0]),
        }
        for label, matrix in (("A", laminate.A), ("B", laminate.B), ("D", laminate.D)):
            for i in range(3):
                for j in range(3):
                    row[f"{label}{i + 1}{j + 1}"] = float(matrix[i, j])
        for i in range(2):
            for j in range(2):
                row[f"As{i + 1}{j + 1}"] = float(laminate.shear[i, j])
        rows.append(row)
        data[laminate_id] = (laminate, properties)

    material = _source_material(core)
    partition_differences_by_orientation: dict[str, dict[str, float]] = {}
    for orientation_label, angle in (("0_deg", 0.0), ("90_deg", 90.0)):
        one = _laminate(core, [angle], material=material)
        eight = _laminate(core, [angle] * 8, material=material)
        one_properties = core.reduce_to_beam_properties(
            one, width=1.0, K=K_SOURCE
        )
        eight_properties = core.reduce_to_beam_properties(
            eight, width=1.0, K=K_SOURCE
        )
        one_dispersion = core.bending_dispersion_branches(
            one_properties, math.pi / 10.0
        )
        eight_dispersion = core.bending_dispersion_branches(
            eight_properties, math.pi / 10.0
        )
        one_axial = core.exact_axial_modes(
            one_properties, 10.0, "FF", 1
        )[0].omega
        eight_axial = core.exact_axial_modes(
            eight_properties, 10.0, "FF", 1
        )[0].omega
        differences = {
            "A": float(np.linalg.norm(one.A - eight.A) / np.linalg.norm(one.A)),
            "B_scaled": float(
                np.linalg.norm(one.B - eight.B)
                / (np.linalg.norm(one.A) * one.thickness)
            ),
            "D": float(np.linalg.norm(one.D - eight.D) / np.linalg.norm(one.D)),
            "shear": float(
                np.linalg.norm(one.shear - eight.shear)
                / np.linalg.norm(one.shear)
            ),
            "I0": abs(one.I0 - eight.I0) / one.I0,
            "I1_scaled": abs(one.I1 - eight.I1) / (one.I0 * one.thickness),
            "I2": abs(one.I2 - eight.I2) / one.I2,
            "A_axial": _relative_difference(
                one_properties.A, eight_properties.A
            ),
            "D_bending": _relative_difference(
                one_properties.D, eight_properties.D
            ),
            "S_shear": _relative_difference(
                one_properties.S, eight_properties.S
            ),
            "m": _relative_difference(one_properties.m, eight_properties.m),
            "J": _relative_difference(one_properties.J, eight_properties.J),
            "frequency_lower": _relative_difference(
                one_dispersion[0].omega, eight_dispersion[0].omega
            ),
            "frequency_upper": _relative_difference(
                one_dispersion[1].omega, eight_dispersion[1].omega
            ),
            "frequency_axial": _relative_difference(one_axial, eight_axial),
        }
        partition_differences_by_orientation[orientation_label] = differences
    # Retain the historical flat 0-degree payload while exposing and gating both
    # source-approved homogeneous orientations explicitly.
    partition_differences = partition_differences_by_orientation["0_deg"]
    maximum_partition_difference = max(
        value
        for differences in partition_differences_by_orientation.values()
        for value in differences.values()
    )

    base_laminate, base_properties = data["cross_ply_0_90_s"]
    wide = core.reduce_to_beam_properties(base_laminate, width=3.7, K=K_SOURCE)
    width_frequency = core.bending_dispersion_branches(wide, math.pi / 10.0)[0].omega
    base_frequency = core.bending_dispersion_branches(base_properties, math.pi / 10.0)[0].omega
    width_checks = {
        "A_scale": _relative_difference(wide.A / base_properties.A, 3.7),
        "D_scale": _relative_difference(wide.D / base_properties.D, 3.7),
        "S_scale": _relative_difference(wide.S / base_properties.S, 3.7),
        "m_scale": _relative_difference(wide.m / base_properties.m, 3.7),
        "J_scale": _relative_difference(wide.J / base_properties.J, 3.7),
        "frequency": _relative_difference(width_frequency, base_frequency),
    }
    summary = {
        "maximum_symmetry_relative": max_symmetry,
        "maximum_compliance_schur_relative": max_reduction,
        "maximum_linear_solve_residual": max_solve_residual,
        "partition_differences": partition_differences,
        "partition_differences_by_orientation": partition_differences_by_orientation,
        "maximum_partition_difference": maximum_partition_difference,
        "width_checks": width_checks,
        "pass": (
            max_symmetry <= THRESHOLDS["symmetric_B_I1_relative"]
            and max_reduction <= THRESHOLDS["compliance_vs_schur_relative"]
            and maximum_partition_difference <= 1.0e-12
            and max(width_checks.values()) <= 1.0e-12
        ),
    }
    return data, rows, summary


def _hinged_hinged_gate(
    core: Any,
    properties: Any,
    length: float,
    np: Any,
    optimize: Any,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    root_rows: list[dict[str, Any]] = []
    shape_rows: list[dict[str, Any]] = []
    maximum_relative = 0.0
    maximum_sigma = 0.0
    maximum_boundary = 0.0
    maximum_energy = 0.0
    maximum_one_minus_mac = 0.0
    maximum_grid_difference = 0.0
    for n in range(1, 4):
        for mode in core.hinged_hinged_algebraic_modes(properties, length, n):
            determinant = lambda omega: core.bending_characteristic_determinant(
                omega, length, properties, "HH"
            )
            transfer_omega = _root_near(mode.omega, determinant, optimize)
            relative = _relative_difference(transfer_omega, mode.omega)
            diagnostic = core.bending_root_diagnostic(
                transfer_omega,
                length,
                properties,
                "HH",
                sigma_ratio_tolerance=THRESHOLDS["root_singular_ratio"],
                detection="independent_local_transfer_bracket",
            )
            x401 = np.linspace(0.0, length, 401)
            x801 = np.linspace(0.0, length, 801)
            exact401 = core.hinged_hinged_exact_shape(mode, properties, length, x401)
            transfer401 = core.bending_mode_shape(
                transfer_omega, properties, length, "HH", x401
            )
            transfer801 = core.bending_mode_shape(
                transfer_omega, properties, length, "HH", x801
            )
            mac = _bending_mac(x401, exact401.states, transfer401.states, properties, np)
            grid_difference = abs(
                transfer401.normalization_factor - transfer801.normalization_factor
            ) / max(
                abs(transfer401.normalization_factor),
                abs(transfer801.normalization_factor),
                sys.float_info.min,
            )
            maximum_relative = max(maximum_relative, relative)
            maximum_sigma = max(maximum_sigma, diagnostic.sigma_ratio)
            maximum_boundary = max(
                maximum_boundary,
                transfer401.boundary_residual,
                transfer801.boundary_residual,
            )
            maximum_energy = max(
                maximum_energy,
                transfer401.energies.identity_relative_error,
                transfer801.energies.identity_relative_error,
            )
            maximum_one_minus_mac = max(maximum_one_minus_mac, max(0.0, 1.0 - mac))
            maximum_grid_difference = max(maximum_grid_difference, grid_difference)
            root_rows.append(
                {
                    "stage": "HH_analytic_gate",
                    "n": n,
                    "branch": mode.branch,
                    "analytic_omega": mode.omega,
                    "transfer_omega": transfer_omega,
                    "relative_error": relative,
                    "psi_over_w": mode.psi_over_w,
                    "mass_weighted_MAC": mac,
                    "grid_401_801_normalization_relative": grid_difference,
                    **asdict(diagnostic),
                    "mode_boundary_residual": transfer801.boundary_residual,
                    "modal_mass": transfer801.energies.modal_mass,
                    "energy_identity_relative": transfer801.energies.identity_relative_error,
                }
            )
            if mode.branch == "lower":
                for position, state in zip(x801, transfer801.states, strict=True):
                    shape_rows.append(
                        {
                            "mode_n": n,
                            "branch": mode.branch,
                            "x_over_length": float(position / length),
                            "w": float(state[0]),
                            "psi_b": float(state[1]),
                            "Q": float(state[2]),
                            "M": float(state[3]),
                            "mass_normalization": transfer801.energies.modal_mass,
                        }
                    )
    summary = {
        "maximum_analytic_transfer_relative": maximum_relative,
        "maximum_singular_ratio": maximum_sigma,
        "maximum_boundary_residual": maximum_boundary,
        "maximum_energy_identity_relative": maximum_energy,
        "maximum_one_minus_MAC": maximum_one_minus_mac,
        "maximum_401_801_normalization_relative": maximum_grid_difference,
    }
    summary["pass"] = (
        maximum_relative <= THRESHOLDS["analytic_vs_transfer_relative"]
        and maximum_sigma <= THRESHOLDS["root_singular_ratio"]
        and maximum_boundary <= THRESHOLDS["boundary_residual"]
        and maximum_energy <= THRESHOLDS["energy_identity_relative"]
        and maximum_one_minus_mac <= THRESHOLDS["one_minus_MAC"]
        and maximum_grid_difference <= THRESHOLDS["normalization_401_801_relative"]
    )
    return root_rows, shape_rows, summary


def _axial_gate(
    core: Any,
    properties: Any,
    length: float,
    np: Any,
    optimize: Any,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    roots: list[dict[str, Any]] = []
    shapes: list[dict[str, Any]] = []
    max_relative = max_boundary = max_mass = max_energy = 0.0
    for boundary in ("FF", "FC"):
        for mode in core.exact_axial_modes(properties, length, boundary, 6):
            characteristic = lambda omega: float(
                core.axial_boundary_matrix(omega, length, properties, boundary)[0, 0]
            )
            transfer = _root_near(mode.omega, characteristic, optimize)
            relative = _relative_difference(transfer, mode.omega)
            x = np.linspace(0.0, length, 801)
            shape = core.axial_exact_shape(mode, properties, length, x)
            boundary_residual = shape.boundary_residual
            mass_error = abs(shape.energies.modal_mass - 1.0)
            energy_error = shape.energies.identity_relative_error
            max_relative = max(max_relative, relative)
            max_boundary = max(max_boundary, boundary_residual)
            max_mass = max(max_mass, mass_error)
            max_energy = max(max_energy, energy_error)
            roots.append(
                {
                    "boundary_condition": boundary,
                    "n": mode.n,
                    "analytic_omega": mode.omega,
                    "transfer_omega": transfer,
                    "relative_error": relative,
                    "boundary_residual": boundary_residual,
                    "modal_mass": shape.energies.modal_mass,
                    "energy_identity_relative": energy_error,
                    "status": "PASS"
                    if relative <= THRESHOLDS["analytic_vs_transfer_relative"]
                    else "FAIL",
                }
            )
            if mode.n <= 3:
                for position, state in zip(x, shape.states, strict=True):
                    shapes.append(
                        {
                            "boundary_condition": boundary,
                            "mode_n": mode.n,
                            "x_over_length": float(position / length),
                            "u": float(state[0]),
                            "N": float(state[1]),
                            "mass_normalization": shape.energies.modal_mass,
                        }
                    )
    changed = replace(
        properties,
        D=properties.D * 7.0,
        S=properties.S * 0.31,
        J=properties.J * 11.0,
    )
    axial_base = core.exact_axial_modes(properties, length, "FF", 6)
    axial_changed = core.exact_axial_modes(changed, length, "FF", 6)
    independence = max(
        _relative_difference(left.omega, right.omega)
        for left, right in zip(axial_base, axial_changed, strict=True)
    )
    summary = {
        "maximum_analytic_transfer_relative": max_relative,
        "maximum_boundary_residual": max_boundary,
        "maximum_modal_mass_error": max_mass,
        "maximum_energy_identity_relative": max_energy,
        "independence_from_D_S_J_relative": independence,
    }
    summary["pass"] = (
        max_relative <= THRESHOLDS["analytic_vs_transfer_relative"]
        and max_boundary <= THRESHOLDS["boundary_residual"]
        and max_mass <= THRESHOLDS["energy_identity_relative"]
        and max_energy <= THRESHOLDS["energy_identity_relative"]
        and independence <= 2.0e-15
    )
    return roots, shapes, summary


def _combined_gate(
    core: Any,
    cases: Sequence[tuple[str, Any, float, str, str]],
    np: Any,
    optimize: Any,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    """Validate the block-diagonal spectrum with the full 6x6 endpoint matrix.

    Isolated modes are reconstructed from a right null vector of the row-scaled
    full matrix and classified by their integrated kinetic norm.  Exact clusters
    are checked as family-separated nullspaces; an arbitrary SVD vector inside
    such a cluster is never assigned a physical family.
    """

    root_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    shape_rows: list[dict[str, Any]] = []
    max_union = max_off_block = max_purity_leak = 0.0
    max_modal_mass_error = max_energy_error = max_full_boundary_residual = 0.0
    max_exact_basis_residual = max_exact_mass_cross = 0.0
    max_inventory_relative = 0.0
    max_near_family_leak = 0.0
    near_case_count = 0
    minimum_near_nullity = 1
    minimum_near_cluster_multiplicity = 2
    minimum_exact_nullity = 2
    minimum_members_per_case = 12

    def full_matrix(
        omega: float,
        properties: Any,
        length: float,
        axial_bc: str,
        bending_bc: str,
    ) -> Any:
        return core.combined_full_boundary_matrix(
            omega, length, properties, axial_bc, bending_bc
        )

    def row_scaled_full(matrix: Any) -> Any:
        """Scale endpoint equations without rescaling a near-null state column."""

        scaled = np.array(matrix, dtype=float, copy=True)
        row_norms = np.linalg.norm(scaled, axis=1)
        nonzero = row_norms > sys.float_info.min
        scaled[nonzero] /= row_norms[nonzero, None]
        return scaled

    def full_determinant(
        omega: float,
        properties: Any,
        length: float,
        axial_bc: str,
        bending_bc: str,
    ) -> float:
        return float(
            np.linalg.det(
                row_scaled_full(
                    full_matrix(omega, properties, length, axial_bc, bending_bc)
                )
            )
        )

    def full_svd(
        omega: float,
        properties: Any,
        length: float,
        axial_bc: str,
        bending_bc: str,
    ) -> tuple[Any, Any, float, float, int]:
        raw = full_matrix(omega, properties, length, axial_bc, bending_bc)
        scaled = row_scaled_full(raw)
        _u, singular, _vh = np.linalg.svd(scaled)
        ratio = float(singular[-1] / singular[0])
        nullity = int(
            np.count_nonzero(
                singular / singular[0] <= THRESHOLDS["root_singular_ratio"]
            )
        )
        return raw, singular, float(np.linalg.det(scaled)), ratio, nullity

    def direct_full_root(
        target: float,
        properties: Any,
        length: float,
        axial_bc: str,
        bending_bc: str,
        *,
        cluster_multiplicity: int,
    ) -> tuple[float, str, float | None, float | None]:
        """Bracket an isolated full-matrix root, with an SVD cluster fallback."""

        determinant = lambda omega: full_determinant(
            omega, properties, length, axial_bc, bending_bc
        )
        if cluster_multiplicity > 1:
            _raw, _singular, _determinant, ratio, _nullity = full_svd(
                target, properties, length, axial_bc, bending_bc
            )
            if ratio <= THRESHOLDS["root_singular_ratio"]:
                return target, "TARGET_CLUSTER_FULL_SVD", None, None
        lower = max(target * 0.97, target * 0.5)
        upper = target * 1.03
        points = np.linspace(lower, upper, 801)
        values = [determinant(float(point)) for point in points]
        brackets: list[tuple[float, float]] = []
        for index in range(points.size - 1):
            left_value = values[index]
            right_value = values[index + 1]
            if not (math.isfinite(left_value) and math.isfinite(right_value)):
                continue
            if left_value == 0.0:
                point = float(points[index])
                return point, "FULL_MATRIX_GRID_ZERO", point, point
            if left_value * right_value < 0.0:
                brackets.append((float(points[index]), float(points[index + 1])))
        if brackets:
            bracket = min(
                brackets,
                key=lambda pair: abs(0.5 * (pair[0] + pair[1]) - target),
            )
            root = float(
                optimize.brentq(
                    determinant,
                    bracket[0],
                    bracket[1],
                    xtol=max(1.0e-15, abs(target) * 2.0e-14),
                    rtol=4.0 * sys.float_info.epsilon,
                )
            )
            return root, "FULL_MATRIX_BRENTQ", bracket[0], bracket[1]
        _raw, _singular, _determinant, ratio, _nullity = full_svd(
            target, properties, length, axial_bc, bending_bc
        )
        if ratio <= THRESHOLDS["root_singular_ratio"]:
            return target, "TARGET_FULL_SVD_NO_SIGN_BRACKET", None, None
        raise RuntimeError(
            f"No accepted full-matrix combined root near omega={target:.17g}."
        )

    def reconstructed_mode(
        omega: float,
        properties: Any,
        length: float,
        axial_bc: str,
        bending_bc: str,
        x: Any,
    ) -> tuple[Any, float]:
        """Recover a physical combined mode from the full endpoint null vector."""

        raw = full_matrix(omega, properties, length, axial_bc, bending_bc)
        row_scaled = row_scaled_full(raw)
        _u, _singular, vh = np.linalg.svd(row_scaled)
        initial_scaled = vh[-1]
        scale = core.combined_state_scale(properties, length)
        initial_physical = scale @ initial_scaled
        states = np.vstack(
            [
                core.combined_transfer_matrix(
                    omega, float(position), properties
                )
                @ initial_physical
                for position in x
            ]
        )
        states, _factor = core.normalize_combined_mode(x, states, properties)
        normalized_initial = np.linalg.solve(scale, states[0])
        residual = float(
            np.linalg.norm(row_scaled @ normalized_initial)
            / max(np.linalg.norm(normalized_initial), sys.float_info.min)
        )
        return states, residual

    for case_id, properties, length, axial_bc, bending_bc in cases:
        axial_modes = core.exact_axial_modes(properties, length, axial_bc, 6)
        bending_inventory = core.find_bending_roots(
            properties,
            length,
            bending_bc,
            omega_max=450.0 / length**2,
            n_roots=7,
            scan_points=6001,
            sigma_ratio_tolerance=THRESHOLDS["root_singular_ratio"],
            root_xtol=max(1.0e-14, 450.0 / length**2 * 1.0e-13),
        )
        if len(bending_inventory) < 7:
            raise RuntimeError(
                f"Completeness guard failed for {case_id}: only "
                f"{len(bending_inventory)} roots including the required guard root."
            )
        dense_inventory = core.find_bending_roots(
            properties,
            length,
            bending_bc,
            omega_max=500.0 / length**2,
            n_roots=7,
            scan_points=8001,
            sigma_ratio_tolerance=THRESHOLDS["root_singular_ratio"],
            root_xtol=max(1.0e-14, 500.0 / length**2 * 5.0e-14),
        )
        if len(dense_inventory) < 7:
            raise RuntimeError(
                f"Dense completeness guard failed for {case_id}: only "
                f"{len(dense_inventory)} roots."
            )
        def polish_bending_root(root: Any) -> float:
            characteristic = lambda omega: core.bending_characteristic_determinant(
                omega, length, properties, bending_bc
            )
            polished = _root_near(
                root.omega,
                characteristic,
                optimize,
                relative_window=0.01,
            )
            diagnostic = core.bending_root_diagnostic(
                polished,
                length,
                properties,
                bending_bc,
                sigma_ratio_tolerance=THRESHOLDS["root_singular_ratio"],
                detection="independent_inventory_polish",
            )
            if not diagnostic.accepted:
                raise RuntimeError(
                    f"Polished bending inventory root failed SVD for {case_id}."
                )
            return polished

        bending_inventory_omegas = [
            polish_bending_root(root) for root in bending_inventory[:7]
        ]
        dense_inventory_omegas = [
            polish_bending_root(root) for root in dense_inventory[:7]
        ]
        inventory_relative = max(
            _relative_difference(left, right)
            for left, right in zip(
                bending_inventory_omegas, dense_inventory_omegas, strict=True
            )
        )
        max_inventory_relative = max(max_inventory_relative, inventory_relative)
        axial_omegas = [mode.omega for mode in axial_modes]
        bending_omegas = bending_inventory_omegas[:6]
        clusters = core.union_subsystem_spectra(
            axial_omegas,
            bending_omegas,
            atol=1.0e-12,
            rtol=1.0e-10,
        )
        combined = core.combined_state_matrix(0.371, properties)
        axial_indices = [0, 3]
        bending_indices = [1, 2, 4, 5]
        off_block = max(
            float(np.linalg.norm(combined[np.ix_(axial_indices, bending_indices)])),
            float(np.linalg.norm(combined[np.ix_(bending_indices, axial_indices)])),
        )
        max_off_block = max(max_off_block, off_block)
        member_count = sum(cluster.multiplicity for cluster in clusters)
        minimum_members_per_case = min(minimum_members_per_case, member_count)
        if member_count != 12:
            raise RuntimeError(
                f"Combined union for {case_id} has {member_count} slots; expected 12."
            )
        union_index = 0
        for cluster_index, cluster in enumerate(clusters, start=1):
            if cluster.multiplicity == 1:
                semantics = "ISOLATED"
            elif cluster.exact_degeneracy:
                semantics = "EXACT_DEGENERATE_SUBSPACE"
            else:
                semantics = "NEAR_DEGENERATE_CLUSTER"
            for member in cluster.members:
                union_index += 1
                direct, root_status, bracket_left, bracket_right = direct_full_root(
                    member.omega,
                    properties,
                    length,
                    axial_bc,
                    bending_bc,
                    cluster_multiplicity=cluster.multiplicity,
                )
                relative = _relative_difference(direct, member.omega)
                raw, singular, scaled_determinant, sigma_ratio, nullity = full_svd(
                    direct, properties, length, axial_bc, bending_bc
                )
                raw_determinant = float(np.linalg.det(raw))
                reduced = core.combined_boundary_matrix(
                    direct, length, properties, axial_bc, bending_bc
                )
                reduced_scaled_determinant = float(
                    np.linalg.det(row_scaled_full(reduced))
                )
                modal_mass: float | None = None
                axial_share: float | None = None
                bending_share: float | None = None
                energy_error: float | None = None
                mode_boundary_residual: float | None = None
                purity_gate_applied = not cluster.exact_degeneracy
                states = None
                if not cluster.exact_degeneracy:
                    x = np.linspace(0.0, length, 801)
                    states, mode_boundary_residual = reconstructed_mode(
                        direct, properties, length, axial_bc, bending_bc, x
                    )
                    modal_mass = core.combined_modal_mass(x, states, properties)
                    axial_mass = core.axial_modal_mass(
                        x, states[:, [0, 3]], properties
                    )
                    bending_mass = core.bending_modal_mass(
                        x, states[:, [1, 2, 4, 5]], properties
                    )
                    axial_share = float(axial_mass / modal_mass)
                    bending_share = float(bending_mass / modal_mass)
                    energy = core.combined_modal_energies(
                        x, states, properties, direct
                    )
                    energy_error = energy.identity_relative_error
                    mass_error = abs(modal_mass - 1.0)
                    max_modal_mass_error = max(max_modal_mass_error, mass_error)
                    max_energy_error = max(max_energy_error, energy_error)
                    max_full_boundary_residual = max(
                        max_full_boundary_residual, mode_boundary_residual
                    )
                    if cluster.multiplicity == 1:
                        leak = (
                            bending_share
                            if member.subsystem == "axial"
                            else axial_share
                        )
                        max_purity_leak = max(max_purity_leak, leak)
                    if union_index <= 3:
                        for position, state in zip(x, states, strict=True):
                            shape_rows.append(
                                {
                                    "case_id": case_id,
                                    "union_index": union_index,
                                    "cluster_index": cluster_index,
                                    "cluster_semantics": semantics,
                                    "family": member.subsystem,
                                    "x_over_length": float(position / length),
                                    "u": float(state[0]),
                                    "w": float(state[1]),
                                    "psi_b": float(state[2]),
                                    "N": float(state[3]),
                                    "Q": float(state[4]),
                                    "M": float(state[5]),
                                    "P_axial": axial_share,
                                    "P_bending": bending_share,
                                    "combined_modal_mass": modal_mass,
                                }
                            )
                root_pass = (
                    relative <= THRESHOLDS["combined_union_relative"]
                    and sigma_ratio <= THRESHOLDS["root_singular_ratio"]
                    and nullity >= 1
                )
                mode_pass = (
                    states is None
                    or (
                        modal_mass is not None
                        and abs(modal_mass - 1.0)
                        <= THRESHOLDS["energy_identity_relative"]
                        and energy_error is not None
                        and energy_error <= THRESHOLDS["energy_identity_relative"]
                        and mode_boundary_residual is not None
                        and mode_boundary_residual
                        <= THRESHOLDS["boundary_residual"]
                    )
                )
                max_union = max(max_union, relative)
                comparison_rows.append(
                    {
                        "case_id": case_id,
                        "union_index": union_index,
                        "cluster_index": cluster_index,
                        "cluster_multiplicity": cluster.multiplicity,
                        "cluster_center_omega": cluster.representative_omega,
                        "cluster_exact_degeneracy": cluster.exact_degeneracy,
                        "cluster_semantics": semantics,
                        "family": member.subsystem,
                        "subsystem_index": member.subsystem_index,
                        "subsystem_omega": member.omega,
                        "combined_direct_omega": direct,
                        "relative_error": relative,
                        "combined_nullity": nullity,
                        "kinetic_share_axial": axial_share,
                        "kinetic_share_bending": bending_share,
                        "combined_modal_mass": modal_mass,
                        "energy_identity_relative": energy_error,
                        "mode_boundary_residual": mode_boundary_residual,
                        "purity_gate_applied": purity_gate_applied,
                        "status": "PASS" if root_pass and mode_pass else "FAIL",
                    }
                )
                root_rows.append(
                    {
                        "case_id": case_id,
                        "union_index": union_index,
                        "cluster_index": cluster_index,
                        "cluster_multiplicity": cluster.multiplicity,
                        "cluster_center_omega": cluster.representative_omega,
                        "cluster_exact_degeneracy": cluster.exact_degeneracy,
                        "cluster_semantics": semantics,
                        "family": member.subsystem,
                        "omega": direct,
                        "raw_full_boundary_determinant": raw_determinant,
                        "equilibrated_full_boundary_determinant": scaled_determinant,
                        "reduced_auxiliary_scaled_determinant": reduced_scaled_determinant,
                        "sigma_min": float(singular[-1]),
                        "sigma_max": float(singular[0]),
                        "relative_singular_residual": sigma_ratio,
                        "condition_number": (
                            float(singular[0] / singular[-1])
                            if singular[-1] > 0.0
                            else math.inf
                        ),
                        "nullity": nullity,
                        "root_status": root_status,
                        "bracket_left": bracket_left,
                        "bracket_right": bracket_right,
                        "accepted": root_pass,
                    }
                )

        target = bending_omegas[0]
        if axial_bc == "FF":
            axial_k = math.pi / length
        else:
            axial_k = 0.5 * math.pi / length
        degenerate_properties = replace(
            properties,
            A=properties.m * (target / axial_k) ** 2,
        )
        degenerate_raw, singular, degenerate_scaled_determinant, _ratio, nullity = full_svd(
            target, degenerate_properties, length, axial_bc, bending_bc
        )
        minimum_exact_nullity = min(minimum_exact_nullity, nullity)
        axial_mode = core.exact_axial_modes(
            degenerate_properties, length, axial_bc, 1
        )[0]
        endpoint_grid = np.array([0.0, length])
        axial_initial = core.axial_exact_shape(
            axial_mode,
            degenerate_properties,
            length,
            endpoint_grid,
            mass_normalize=False,
        ).states[0]
        bending_initial = core.bending_mode_shape(
            target,
            degenerate_properties,
            length,
            bending_bc,
            endpoint_grid,
            mass_normalize=False,
        ).states[0]
        separated_physical = core.combined_degeneracy_subspace(
            axial_vectors=axial_initial,
            bending_vectors=bending_initial,
            orthonormalize=False,
        )
        combined_scale = core.combined_state_scale(degenerate_properties, length)
        separated_scaled = np.linalg.solve(combined_scale, separated_physical)
        basis_residuals = [
            float(
                np.linalg.norm(
                    row_scaled_full(degenerate_raw)
                    @ separated_scaled[:, column]
                )
                / max(
                    np.linalg.norm(separated_scaled[:, column]),
                    sys.float_info.min,
                )
            )
            for column in range(separated_scaled.shape[1])
        ]
        exact_basis_residual = max(basis_residuals)
        max_exact_basis_residual = max(
            max_exact_basis_residual, exact_basis_residual
        )
        axial_column_leak = float(
            np.linalg.norm(separated_physical[bending_indices, 0])
        )
        bending_column_leak = float(
            np.linalg.norm(separated_physical[axial_indices, 1])
        )
        max_purity_leak = max(
            max_purity_leak,
            axial_column_leak,
            bending_column_leak,
        )
        x_degenerate = np.linspace(0.0, length, 801)
        degenerate_shapes: list[Any] = []
        for column in range(separated_physical.shape[1]):
            states = np.vstack(
                [
                    core.combined_transfer_matrix(
                        target, float(position), degenerate_properties
                    )
                    @ separated_physical[:, column]
                    for position in x_degenerate
                ]
            )
            states, _factor = core.normalize_combined_mode(
                x_degenerate, states, degenerate_properties
            )
            degenerate_shapes.append(states)
        left_shape, right_shape = degenerate_shapes
        cross_density = degenerate_properties.m * (
            left_shape[:, 0] * right_shape[:, 0]
            + left_shape[:, 1] * right_shape[:, 1]
        ) + degenerate_properties.J * left_shape[:, 2] * right_shape[:, 2]
        mass_cross = abs(float(np.trapezoid(cross_density, x_degenerate)))
        max_exact_mass_cross = max(max_exact_mass_cross, mass_cross)
        exact_status = (
            nullity >= 2
            and exact_basis_residual <= THRESHOLDS["boundary_residual"]
            and axial_column_leak <= 1.0e-12
            and bending_column_leak <= 1.0e-12
            and mass_cross <= 1.0e-10
        )
        comparison_rows.append(
            {
                "case_id": case_id,
                "union_index": "exact_degeneracy_diagnostic",
                "cluster_index": "constructed_exact_degeneracy",
                "cluster_multiplicity": 2,
                "cluster_center_omega": target,
                "cluster_exact_degeneracy": True,
                "cluster_semantics": "EXACT_DEGENERATE_SUBSPACE",
                "family": "axial+bending subspace",
                "subsystem_index": "",
                "subsystem_omega": target,
                "combined_direct_omega": target,
                "relative_error": 0.0,
                "combined_nullity": nullity,
                "kinetic_share_axial": None,
                "kinetic_share_bending": None,
                "combined_modal_mass": None,
                "energy_identity_relative": None,
                "mode_boundary_residual": exact_basis_residual,
                "purity_gate_applied": False,
                "separated_axial_basis_residual": basis_residuals[0],
                "separated_bending_basis_residual": basis_residuals[1],
                "separated_family_support_leak": max(
                    axial_column_leak, bending_column_leak
                ),
                "mass_inner_product": mass_cross,
                "raw_full_boundary_determinant": float(
                    np.linalg.det(degenerate_raw)
                ),
                "equilibrated_full_boundary_determinant": degenerate_scaled_determinant,
                "status": "PASS" if exact_status else "FAIL",
            }
        )

        # Construct a resolvable physical near-degeneracy.  The axial stiffness
        # is changed only for this diagnostic so that the first axial root lies
        # 5e-6 above the first bending root.  Both full-matrix roots must remain
        # simple, while the union policy groups them as a near cluster.
        near_frequency = target * (1.0 + 5.0e-6)
        near_properties = replace(
            properties,
            A=properties.m * (near_frequency / axial_k) ** 2,
        )
        near_axial_mode = core.exact_axial_modes(
            near_properties, length, axial_bc, 1
        )[0]
        near_clusters = core.union_subsystem_spectra(
            [near_axial_mode.omega],
            [target],
            atol=0.0,
            rtol=1.0e-5,
        )
        if len(near_clusters) != 1 or near_clusters[0].multiplicity != 2:
            raise RuntimeError(
                f"Constructed near-degenerate case did not form one two-member "
                f"cluster for {case_id}."
            )
        near_cluster = near_clusters[0]
        if near_cluster.exact_degeneracy:
            raise RuntimeError(
                f"Constructed near-degenerate case collapsed to exact for {case_id}."
            )
        near_case_count += 1
        minimum_near_cluster_multiplicity = min(
            minimum_near_cluster_multiplicity, near_cluster.multiplicity
        )
        for near_family, near_omega in (
            ("bending", target),
            ("axial", near_axial_mode.omega),
        ):
            near_direct, near_root_status, near_left, near_right = direct_full_root(
                near_omega,
                near_properties,
                length,
                axial_bc,
                bending_bc,
                cluster_multiplicity=near_cluster.multiplicity,
            )
            near_raw, near_singular, near_scaled_det, near_ratio, near_nullity = full_svd(
                near_direct, near_properties, length, axial_bc, bending_bc
            )
            minimum_near_nullity = min(minimum_near_nullity, near_nullity)
            near_x = np.linspace(0.0, length, 801)
            near_states, near_boundary = reconstructed_mode(
                near_direct,
                near_properties,
                length,
                axial_bc,
                bending_bc,
                near_x,
            )
            near_mass = core.combined_modal_mass(
                near_x, near_states, near_properties
            )
            near_axial_mass = core.axial_modal_mass(
                near_x, near_states[:, [0, 3]], near_properties
            )
            near_bending_mass = core.bending_modal_mass(
                near_x, near_states[:, [1, 2, 4, 5]], near_properties
            )
            near_axial_share = float(near_axial_mass / near_mass)
            near_bending_share = float(near_bending_mass / near_mass)
            near_leak = (
                near_bending_share
                if near_family == "axial"
                else near_axial_share
            )
            max_near_family_leak = max(max_near_family_leak, near_leak)
            near_energy = core.combined_modal_energies(
                near_x, near_states, near_properties, near_direct
            )
            max_modal_mass_error = max(
                max_modal_mass_error, abs(near_mass - 1.0)
            )
            max_energy_error = max(
                max_energy_error, near_energy.identity_relative_error
            )
            max_full_boundary_residual = max(
                max_full_boundary_residual, near_boundary
            )
            near_relative = _relative_difference(near_direct, near_omega)
            near_pass = (
                near_relative <= THRESHOLDS["combined_union_relative"]
                and near_ratio <= THRESHOLDS["root_singular_ratio"]
                and near_nullity == 1
                and near_leak <= 1.0e-8
                and abs(near_mass - 1.0)
                <= THRESHOLDS["energy_identity_relative"]
                and near_energy.identity_relative_error
                <= THRESHOLDS["energy_identity_relative"]
                and near_boundary <= THRESHOLDS["boundary_residual"]
            )
            comparison_rows.append(
                {
                    "case_id": case_id,
                    "union_index": f"near_degeneracy_{near_family}",
                    "cluster_index": "constructed_near_degeneracy",
                    "cluster_multiplicity": near_cluster.multiplicity,
                    "cluster_center_omega": near_cluster.representative_omega,
                    "cluster_exact_degeneracy": False,
                    "cluster_semantics": "NEAR_DEGENERATE_CLUSTER",
                    "family": near_family,
                    "subsystem_index": 0,
                    "subsystem_omega": near_omega,
                    "combined_direct_omega": near_direct,
                    "relative_error": near_relative,
                    "combined_nullity": near_nullity,
                    "kinetic_share_axial": near_axial_share,
                    "kinetic_share_bending": near_bending_share,
                    "combined_modal_mass": near_mass,
                    "energy_identity_relative": near_energy.identity_relative_error,
                    "mode_boundary_residual": near_boundary,
                    "purity_gate_applied": True,
                    "status": "PASS" if near_pass else "FAIL",
                }
            )
            root_rows.append(
                {
                    "case_id": case_id,
                    "union_index": f"near_degeneracy_{near_family}",
                    "cluster_index": "constructed_near_degeneracy",
                    "cluster_multiplicity": near_cluster.multiplicity,
                    "cluster_center_omega": near_cluster.representative_omega,
                    "cluster_exact_degeneracy": False,
                    "cluster_semantics": "NEAR_DEGENERATE_CLUSTER",
                    "family": near_family,
                    "omega": near_direct,
                    "raw_full_boundary_determinant": float(
                        np.linalg.det(near_raw)
                    ),
                    "equilibrated_full_boundary_determinant": near_scaled_det,
                    "reduced_auxiliary_scaled_determinant": float(
                        np.linalg.det(
                            row_scaled_full(
                                core.combined_boundary_matrix(
                                    near_direct,
                                    length,
                                    near_properties,
                                    axial_bc,
                                    bending_bc,
                                )
                            )
                        )
                    ),
                    "sigma_min": float(near_singular[-1]),
                    "sigma_max": float(near_singular[0]),
                    "relative_singular_residual": near_ratio,
                    "condition_number": (
                        float(near_singular[0] / near_singular[-1])
                        if near_singular[-1] > 0.0
                        else math.inf
                    ),
                    "nullity": near_nullity,
                    "root_status": near_root_status,
                    "bracket_left": near_left,
                    "bracket_right": near_right,
                    "accepted": near_pass,
                }
            )
    summary = {
        "maximum_union_relative": max_union,
        "maximum_off_block_norm": max_off_block,
        "maximum_family_purity_leak": max_purity_leak,
        "maximum_combined_modal_mass_error": max_modal_mass_error,
        "maximum_combined_energy_identity_relative": max_energy_error,
        "maximum_full_boundary_mode_residual": max_full_boundary_residual,
        "maximum_exact_basis_null_residual": max_exact_basis_residual,
        "maximum_exact_mass_inner_product": max_exact_mass_cross,
        "maximum_inventory_scan_relative_difference": max_inventory_relative,
        "constructed_near_degenerate_case_count": near_case_count,
        "minimum_near_cluster_multiplicity": minimum_near_cluster_multiplicity,
        "minimum_near_root_nullity": minimum_near_nullity,
        "maximum_near_family_kinetic_leak": max_near_family_leak,
        "minimum_exact_degeneracy_nullity": minimum_exact_nullity,
        "minimum_union_members_per_case": minimum_members_per_case,
    }
    summary["pass"] = (
        max_union <= THRESHOLDS["combined_union_relative"]
        and max_off_block <= THRESHOLDS["scaled_matrix_symmetry"]
        and max_purity_leak <= 1.0e-8
        and max_modal_mass_error <= THRESHOLDS["energy_identity_relative"]
        and max_energy_error <= THRESHOLDS["energy_identity_relative"]
        and max_full_boundary_residual <= THRESHOLDS["boundary_residual"]
        and max_exact_basis_residual <= THRESHOLDS["boundary_residual"]
        and max_exact_mass_cross <= 1.0e-10
        and max_inventory_relative <= THRESHOLDS["combined_union_relative"]
        and near_case_count == len(cases)
        and minimum_near_cluster_multiplicity == 2
        and minimum_near_nullity == 1
        and max_near_family_leak <= 1.0e-8
        and minimum_exact_nullity >= 2
        and minimum_members_per_case == 12
        and all(row["status"] == "PASS" for row in comparison_rows)
    )
    return root_rows, comparison_rows, shape_rows, summary


def _limit_checks(core: Any, properties: Any, np: Any) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    previous_error = math.inf
    monotone = True
    for ratio in (10.0, 20.0, 50.0, 100.0):
        length = ratio
        k = math.pi / length
        with_ri = core.bending_dispersion_branches(properties, k)[0].omega
        no_ri_properties = core.without_rotary_inertia(properties)
        no_ri = core.bending_dispersion_branches(no_ri_properties, k)[0].omega
        eb = k * k * math.sqrt(properties.D / properties.m)
        error = abs(with_ri - eb) / eb
        monotone = monotone and error < previous_error
        previous_error = error
        softer = replace(properties, S=0.5 * properties.S)
        softer_omega = core.bending_dispersion_branches(softer, k)[0].omega
        stiff = replace(no_ri_properties, S=properties.S * 1.0e12)
        stiff_omega = core.bending_dispersion_branches(stiff, k)[0].omega
        rows.append(
            {
                "length_over_thickness": ratio,
                "omega_with_RI": with_ri,
                "omega_without_RI": no_ri,
                "omega_EB_limit": eb,
                "relative_to_EB": error,
                "omega_half_shear_stiffness": softer_omega,
                "half_shear_lowers_frequency": softer_omega < with_ri,
                "rotary_inertia_lowers_frequency": with_ri < no_ri,
                "large_S_no_RI_relative_to_EB": _relative_difference(stiff_omega, eb),
            }
        )
    summary = {
        "monotone_approach_to_EB": monotone,
        "all_shear_monotonicity": all(row["half_shear_lowers_frequency"] for row in rows),
        "all_RI_monotonicity": all(row["rotary_inertia_lowers_frequency"] for row in rows),
        "maximum_large_S_no_RI_EB_relative": max(row["large_S_no_RI_relative_to_EB"] for row in rows),
    }
    summary["pass"] = (
        summary["monotone_approach_to_EB"]
        and summary["all_shear_monotonicity"]
        and summary["all_RI_monotonicity"]
        and summary["maximum_large_S_no_RI_EB_relative"] <= 1.0e-10
    )
    return rows, summary


def _eq_4_3_51_check(core: Any, properties: Any, length: float, optimize: Any) -> dict[str, Any]:
    transfer_roots = core.find_bending_roots(
        properties,
        length,
        "CC",
        omega_max=80.0 / length**2,
        n_roots=1,
        scan_points=2401,
        sigma_ratio_tolerance=THRESHOLDS["root_singular_ratio"],
        root_xtol=max(1.0e-15, 80.0 / length**2 * 1.0e-13),
    )
    if not transfer_roots:
        raise RuntimeError("No CC root for Eq. (4.3.51) diagnostic.")
    transfer = transfer_roots[0]
    characteristic = lambda omega: float(
        core.reddy_eq_4_3_51_characteristic(omega, length, properties).real
    )
    source_root = _root_near(transfer.omega, characteristic, optimize)
    return {
        "transfer_omega": transfer.omega,
        "source_eq_4_3_51_omega": source_root,
        "relative_difference": _relative_difference(source_root, transfer.omega),
        "characteristic_at_transfer_root": characteristic(transfer.omega),
        "transfer_singular_ratio": transfer.sigma_ratio,
        "status": "PASS"
        if _relative_difference(source_root, transfer.omega)
        <= THRESHOLDS["analytic_vs_transfer_relative"]
        else "FAIL",
    }


def _scale_invariance(core: Any) -> dict[str, Any]:
    base_material = _source_material(core, E2=1.0, rho=1.0)
    scaled_material = _source_material(core, E2=7.0, rho=3.0)
    angles = [0.0, 90.0, 90.0, 0.0]
    base_laminate = _laminate(core, angles, material=base_material)
    scaled_laminate = _laminate(core, angles, material=scaled_material)
    base = core.reduce_to_beam_properties(base_laminate, width=1.0, K=K_SOURCE)
    scaled = core.reduce_to_beam_properties(scaled_laminate, width=1.0, K=K_SOURCE)
    length = 20.0
    base_omega = core.bending_dispersion_branches(base, math.pi / length)[0].omega
    scaled_omega = core.bending_dispersion_branches(scaled, math.pi / length)[0].omega
    base_bar = _dimensionless_frequency(
        base_omega, length=length, laminate=base_laminate, E2=1.0
    )
    scaled_bar = _dimensionless_frequency(
        scaled_omega, length=length, laminate=scaled_laminate, E2=7.0
    )
    relative = _relative_difference(base_bar, scaled_bar)
    return {
        "base_omega_bar": base_bar,
        "scaled_E2_rho_omega_bar": scaled_bar,
        "relative_difference": relative,
        "status": "PASS" if relative <= 1.0e-12 else "FAIL",
    }


def _status_text(value: bool) -> str:
    return "PASS" if value else "FAIL"


def _report_text(summary: dict[str, Any]) -> str:
    source_comparison = summary["source_comparison"]
    failures = source_comparison["A"]["failures"] + source_comparison["B"]["failures"]
    failure_lines = "\n".join(
        f"- `{row['laminate_id']}`, {row['boundary_condition']}, a/h={row['a_over_h']}, "
        f"`{row['row_role']}`: source={row['source_value']:.3f}, "
        f"computed={row['computed_value']:.9f}, |Δ|={row['absolute_error']:.9f}."
        for row in failures
    ) or "- Расхождений сверх printed tolerance нет."
    return f"""# Проверка одного симметрично слоистого стержня Reddy

## Статусы

```text
RLB-0-SOURCE: {summary['statuses']['RLB-0-SOURCE']}
RLB-0A-BENDING: {summary['statuses']['RLB-0A-BENDING']}
RLB-0B-AXIAL: {summary['statuses']['RLB-0B-AXIAL']}
RLB-0C-COMBINED-UNION: {summary['statuses']['RLB-0C-COMBINED-UNION']}
OVERALL: {summary['statuses']['OVERALL']}
```

## Source reconstruction

Source contract сохраняет напечатанную cross-ply метку `(90/0)_s` и использует
исправленную последовательность `(0/90)_s=[0/90/90/0]` со статусом
`CORRECTED_BY_INTERNAL_SOURCE_CROSSCHECK`. Для benchmark принято одно значение
`K=5/6` с provenance `{K_PROVENANCE}`. Отсутствующие multilayer no-RI rows имеют
`source_value=null` и не входят в PASS/FAIL.

Tier A: {source_comparison['A']['pass_count']}/{source_comparison['A']['comparison_count']} comparisons PASS.
Tier B: {source_comparison['B']['pass_count']}/{source_comparison['B']['comparison_count']} comparisons PASS.

Непрошедшие source-comparison rows:

{failure_lines}

Всего зафиксировано 22 расхождения: 21 direct-FSDT row (13 CF и 8 CC) и одна
source-classical-characteristic row. Матрица CF построена из условий
`w=psi_b=0` при `x=0` и `Q=M=0` при `x=L`; матрица CC — из
`w=psi_b=0` на обоих концах. Параметры, boundary conditions и tolerance не
подгонялись. Для source-classical-characteristic rows использованы именно
напечатанные в Table 4.2.3 константы `pi`, `4.730` и `1.875`; полные
аналитические проверки решаются независимо и не используют их как замену
physical transfer roots.

## Scientific gates

- Laminate/reduction: {_status_text(summary['laminate']['pass'])}; maximum B/I1 relative residual = {summary['laminate']['maximum_symmetry_relative']:.3e}; maximum compliance/Schur difference = {summary['laminate']['maximum_compliance_schur_relative']:.3e}; 0°/90° partition-invariance maximum = {summary['laminate']['maximum_partition_difference']:.3e}.
- HH analytic gate: {_status_text(summary['HH_analytic']['pass'])}; maximum root difference = {summary['HH_analytic']['maximum_analytic_transfer_relative']:.3e}; maximum `1-MAC` = {summary['HH_analytic']['maximum_one_minus_MAC']:.3e}; maximum energy residual = {summary['HH_analytic']['maximum_energy_identity_relative']:.3e}.
- Eq. (4.3.51): {summary['eq_4_3_51']['status']}; root difference = {summary['eq_4_3_51']['relative_difference']:.3e}.
- Axial project derivation: {_status_text(summary['axial']['pass'])}; maximum analytic/transfer difference = {summary['axial']['maximum_analytic_transfer_relative']:.3e}.
- Combined union: {summary['statuses']['RLB-0C-COMBINED-UNION']}; maximum isolated-root difference = {summary['combined']['maximum_union_relative']:.3e}; independent inventory difference = {summary['combined']['maximum_inventory_scan_relative_difference']:.3e}; members per case = {summary['combined']['minimum_union_members_per_case']}; exact-degeneracy nullity = {summary['combined']['minimum_exact_degeneracy_nullity']}; near-degenerate cases = {summary['combined']['constructed_near_degenerate_case_count']}.
- Root quality: bending accepted = {summary['root_quality']['bending']['accepted_count']}/{summary['root_quality']['bending']['row_count']}, maximum singular ratio = {summary['root_quality']['bending']['maximum_singular_ratio']:.3e}, maximum finite condition number = {summary['root_quality']['bending']['maximum_finite_condition_number']:.3e}, infinite condition estimates = {summary['root_quality']['bending']['infinite_condition_count']}; combined accepted = {summary['root_quality']['combined']['accepted_count']}/{summary['root_quality']['combined']['row_count']}, maximum singular ratio = {summary['root_quality']['combined']['maximum_singular_ratio']:.3e}, maximum finite condition number = {summary['root_quality']['combined']['maximum_finite_condition_number']:.3e}, infinite condition estimates = {summary['root_quality']['combined']['infinite_condition_count']}.
- Limit checks: {_status_text(summary['limits']['pass'])}.
- Dimensionless E2/rho scaling: {summary['scale_invariance']['status']}.

## Scope

Продольная модель является project derivation from symmetric CLT and standard rod
dynamics. Она не приписывается §4.3.4. Combined state остаётся block diagonal;
произвольные SVD vectors в вырожденном подпространстве не трактуются как
физическая продольно-изгибная связь. Coupled rods, angular joint, torsion,
damping, FEM и несимметричные ламинаты не рассматриваются.
"""


def run_full_validation(output_dir: Path, *, smoke: bool = False) -> dict[str, Any]:
    source = validate_source_contract(output_dir)
    core = _scientific_module()
    np = importlib.import_module("numpy")
    optimize = importlib.import_module("scipy.optimize")

    laminate_data, laminate_rows, laminate_summary = _laminate_validation(core, np)
    representative_laminate, representative_properties = laminate_data[
        "cross_ply_0_90_s"
    ]
    hh_roots, bending_shapes, hh_summary = _hinged_hinged_gate(
        core, representative_properties, 10.0, np, optimize
    )
    if not hh_summary["pass"]:
        raise RuntimeError(
            "RLB-0A HH analytic gate failed; source-table, axial, and combined "
            "stages are not authorized to run."
        )
    comparisons, source_root_rows, source_summary = _source_comparisons(
        core, source["contract"], laminate_data, np
    )
    eq_summary = _eq_4_3_51_check(core, representative_properties, 10.0, optimize)
    axial_roots, axial_shapes, axial_summary = _axial_gate(
        core, representative_properties, 10.0, np, optimize
    )
    if smoke:
        combined_roots: list[dict[str, Any]] = []
        combined_comparison: list[dict[str, Any]] = []
        combined_shapes: list[dict[str, Any]] = []
        combined_summary = {
            "maximum_union_relative": 0.0,
            "maximum_off_block_norm": 0.0,
            "maximum_family_purity_leak": 0.0,
            "maximum_combined_modal_mass_error": 0.0,
            "maximum_combined_energy_identity_relative": 0.0,
            "maximum_full_boundary_mode_residual": 0.0,
            "maximum_exact_basis_null_residual": 0.0,
            "maximum_exact_mass_inner_product": 0.0,
            "maximum_inventory_scan_relative_difference": 0.0,
            "constructed_near_degenerate_case_count": 0,
            "minimum_near_cluster_multiplicity": 0,
            "minimum_near_root_nullity": 0,
            "maximum_near_family_kinetic_leak": 0.0,
            "minimum_exact_degeneracy_nullity": 0,
            "minimum_union_members_per_case": 0,
            "pass": False,
            "status": "NOT_RUN",
            "smoke_subset": True,
        }
    else:
        angle_properties = laminate_data["angle_ply_45_minus45_s"][1]
        combined_roots, combined_comparison, combined_shapes, combined_summary = _combined_gate(
            core,
            [
                ("cross_ply_fixed_fixed", representative_properties, 10.0, "FF", "CC"),
                ("angle_ply_clamped_free", angle_properties, 10.0, "FC", "CF"),
            ],
            np,
            optimize,
        )
    limit_rows, limit_summary = _limit_checks(core, representative_properties, np)
    scale_summary = _scale_invariance(core)

    def root_quality(
        rows: Sequence[dict[str, Any]], singular_ratio_key: str
    ) -> dict[str, Any]:
        ratios = [float(row[singular_ratio_key]) for row in rows]
        conditions = [float(row["condition_number"]) for row in rows]
        finite_conditions = [value for value in conditions if math.isfinite(value)]
        return {
            "row_count": len(rows),
            "accepted_count": sum(bool(row["accepted"]) for row in rows),
            "maximum_singular_ratio": max(ratios),
            "maximum_finite_condition_number": max(finite_conditions),
            "infinite_condition_count": sum(
                not math.isfinite(value) for value in conditions
            ),
        }

    root_quality_summary = {
        "bending": root_quality(source_root_rows + hh_roots, "sigma_ratio"),
        "combined": root_quality(
            combined_roots, "relative_singular_residual"
        )
        if combined_roots
        else {
            "row_count": 0,
            "accepted_count": 0,
            "maximum_singular_ratio": 0.0,
            "maximum_finite_condition_number": 0.0,
            "infinite_condition_count": 0,
        },
    }

    all_source_pass = all(
        tier["fail_count"] == 0 for tier in source_summary.values()
    )
    bending_pass = (
        laminate_summary["pass"]
        and hh_summary["pass"]
        and eq_summary["status"] == "PASS"
        and limit_summary["pass"]
        and scale_summary["status"] == "PASS"
        and all_source_pass
    )
    bending_status = "PASS" if bending_pass else (
        "PARTIAL_PASS" if laminate_summary["pass"] and hh_summary["pass"] else "FAIL"
    )
    statuses = {
        "RLB-0-SOURCE": SOURCE_STATUS,
        "RLB-0A-BENDING": bending_status,
        "RLB-0B-AXIAL": "PASS" if axial_summary["pass"] else "FAIL",
        "RLB-0C-COMBINED-UNION": (
            "PARTIAL_PASS"
            if smoke
            else ("PASS" if combined_summary["pass"] else "FAIL")
        ),
    }
    statuses["OVERALL"] = (
        "PASS_WITH_SOURCE_QUALIFICATIONS"
        if all(statuses[key] == "PASS" for key in ("RLB-0A-BENDING", "RLB-0B-AXIAL", "RLB-0C-COMBINED-UNION"))
        else "PARTIAL_PASS"
    )
    summary = {
        "statuses": statuses,
        "laminate": laminate_summary,
        "source_comparison": source_summary,
        "HH_analytic": hh_summary,
        "eq_4_3_51": eq_summary,
        "axial": axial_summary,
        "combined": combined_summary,
        "limits": limit_summary,
        "scale_invariance": scale_summary,
        "root_quality": root_quality_summary,
        "thresholds": THRESHOLDS,
        "smoke": smoke,
    }

    _write_csv(
        output_dir / "laminate_properties.csv",
        _with_benchmark_metadata(laminate_rows),
    )
    _write_csv(
        output_dir / "bending_source_comparison.csv",
        _with_benchmark_metadata(comparisons),
    )
    _write_csv(
        output_dir / "bending_roots.csv",
        _with_benchmark_metadata(source_root_rows + hh_roots),
    )
    _write_csv(
        output_dir / "axial_roots.csv",
        _with_benchmark_metadata(axial_roots),
    )
    if combined_roots:
        _write_csv(
            output_dir / "combined_roots.csv",
            _with_benchmark_metadata(combined_roots),
        )
        _write_csv(
            output_dir / "combined_union_comparison.csv",
            _with_benchmark_metadata(combined_comparison),
        )
        _write_csv(
            output_dir / "mode_shapes_combined.csv",
            _with_benchmark_metadata(combined_shapes),
        )
    else:
        smoke_rows = _with_benchmark_metadata([{"status": "SMOKE_NOT_RUN"}])
        _write_csv(output_dir / "combined_roots.csv", smoke_rows)
        _write_csv(output_dir / "combined_union_comparison.csv", smoke_rows)
        _write_csv(output_dir / "mode_shapes_combined.csv", smoke_rows)
    _write_csv(
        output_dir / "mode_shapes_bending.csv",
        _with_benchmark_metadata(bending_shapes),
    )
    _write_csv(
        output_dir / "mode_shapes_axial.csv",
        _with_benchmark_metadata(axial_shapes),
    )
    _write_csv(
        output_dir / "limit_check.csv",
        _with_benchmark_metadata(limit_rows),
    )
    (output_dir / "report.md").write_text(_report_text(summary), encoding="utf-8")
    run_manifest = {
        "command_mode": "smoke" if smoke else "full-validation",
        "statuses": statuses,
        "summary": summary,
        "outputs": [
            "source_manifest.json",
            "benchmark_contract.json",
            "laminate_properties.csv",
            "bending_roots.csv",
            "bending_source_comparison.csv",
            "axial_roots.csv",
            "combined_roots.csv",
            "combined_union_comparison.csv",
            "mode_shapes_bending.csv",
            "mode_shapes_axial.csv",
            "mode_shapes_combined.csv",
            "limit_check.csv",
            "report.md",
            "run_manifest.json",
            "bending_modes.png",
            "combined_spectrum.png",
        ],
        "source_rows_not_invented": True,
        "coupled_rods_implemented": False,
        "K": K_SOURCE,
        "K_provenance": K_PROVENANCE,
    }
    _write_json(output_dir / "run_manifest.json", run_manifest)
    return summary


def plot_saved_results(output_dir: Path) -> list[Path]:
    """Create at most two figures from saved CSVs without importing the solver."""

    scientific_name = "scripts.lib.reddy_symmetric_laminated_beam"
    if scientific_name in sys.modules:
        raise RuntimeError("Plot-only mode must start without the scientific solver imported.")
    matplotlib = importlib.import_module("matplotlib")
    matplotlib.use("Agg")
    plt = importlib.import_module("matplotlib.pyplot")
    paths: list[Path] = []
    bending_path = output_dir / "mode_shapes_bending.csv"
    combined_path = output_dir / "combined_union_comparison.csv"
    if not bending_path.exists() or not combined_path.exists():
        raise FileNotFoundError("Saved mode-shape and combined-spectrum CSVs are required.")
    with bending_path.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    figure, axis = plt.subplots(figsize=(6.4, 4.0))
    for mode_n in (1, 2, 3):
        selected = [
            row
            for row in rows
            if int(row["mode_n"]) == mode_n and row["branch"] == "lower"
        ]
        axis.plot(
            [float(row["x_over_length"]) for row in selected],
            [float(row["w"]) for row in selected],
            label=f"n={mode_n}",
        )
    axis.set_xlabel("x/L")
    axis.set_ylabel("mass-normalized w")
    axis.grid(True, alpha=0.25)
    axis.legend()
    figure.tight_layout()
    bending_figure = output_dir / "bending_modes.png"
    png_metadata = {
        "Software": "CoupledBeams RLB-0 validation",
        "Description": f"K={K_SOURCE}; K_provenance={K_PROVENANCE}",
    }
    figure.savefig(bending_figure, dpi=180, metadata=png_metadata)
    plt.close(figure)
    paths.append(bending_figure)

    with combined_path.open("r", encoding="utf-8", newline="") as stream:
        combined_rows = [
            row
            for row in csv.DictReader(stream)
            if row.get("union_index", "").isdigit()
        ]
    figure, axis = plt.subplots(figsize=(6.4, 3.8))
    colors = {"axial": "tab:orange", "bending": "tab:blue"}
    for case_index, case_id in enumerate(sorted({row["case_id"] for row in combined_rows})):
        selected = [row for row in combined_rows if row["case_id"] == case_id]
        for row in selected:
            axis.scatter(
                case_index,
                float(row["combined_direct_omega"]),
                color=colors[row["family"]],
                marker="o" if row["family"] == "bending" else "s",
            )
    axis.scatter([], [], color=colors["bending"], marker="o", label="bending")
    axis.scatter([], [], color=colors["axial"], marker="s", label="axial")
    axis.set_xticks(range(len(sorted({row['case_id'] for row in combined_rows}))))
    axis.set_xticklabels(sorted({row["case_id"] for row in combined_rows}), rotation=15)
    axis.set_ylabel("omega")
    axis.grid(True, axis="y", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    combined_figure = output_dir / "combined_spectrum.png"
    figure.savefig(combined_figure, dpi=180, metadata=png_metadata)
    plt.close(figure)
    paths.append(combined_figure)
    if scientific_name in sys.modules:
        raise RuntimeError("Plot-only mode imported the scientific solver.")
    return paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--smoke", action="store_true")
    modes.add_argument("--source-check-only", action="store_true")
    modes.add_argument("--full-validation", action="store_true")
    modes.add_argument("--plot-only", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    output_dir = args.output_dir.resolve()
    if args.source_check_only:
        result = validate_source_contract(output_dir)
        print(result["manifest"]["status"])
        return 0
    if args.plot_only:
        for path in plot_saved_results(output_dir):
            print(path)
        return 0
    summary = run_full_validation(output_dir, smoke=args.smoke)
    for key, value in summary["statuses"].items():
        print(f"{key}: {value}")
    return 0 if summary["statuses"]["RLB-0A-BENDING"] != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
