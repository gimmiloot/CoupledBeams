from __future__ import annotations

"""Article-facing energy pilot for selected isotropic circular sorted modes.

The selection phase reads only finalized compact certificates.  It fixes one
geometry for each requested epsilon level from the local frequency-discrepancy
criterion before any mode reconstruction is performed.  The reconstruction
phase evaluates only the six selected Timoshenko modes and never searches for
roots, tracks branches, evaluates MAC, calls FEM, or creates shape figures.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.analysis.thickness_mismatch.audits.audit_longitudinal_suspect_modes_eb_timo import (  # noqa: E402
    ENERGY_SUM_TOL,
    JOINT_FORCE_TOL,
    JOINT_KINEMATIC_TOL,
    timo_joint_continuity_row,
)
from scripts.lib import general_spectrum_completeness as COMPLETE  # noqa: E402
from scripts.lib import variable_length_timoshenko as TIMO  # noqa: E402


SCIENTIFIC_SCOPE = "isotropic_circular_coupled_rods_eb_timoshenko"
CERTIFICATE_SCHEMA = "article_epsilon_compact_certificate_v1"
DEFAULT_FINALIZATION_DIR = (
    REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "coarse_grid_v1"
    / "compact_finalization_epsilon_005_resolved"
)
DEFAULT_COMPACT_DIR = (
    REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "coarse_grid_v1"
    / "compact_point_certificates_v1"
)
DEFAULT_OUTPUT_DIR = (
    REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "energy_triplet_pilot_v1"
)
DEFAULT_COUPLED_OUTPUT_DIR = (
    REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "energy_triplet_coupled_pilot_v2"
)
DEFAULT_BETA0_CONTROL_DIR = DEFAULT_OUTPUT_DIR
DEFAULT_ADDITIONAL_ANGLE_OUTPUT_DIR = (
    REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "energy_triplet_angle_extension_v3"
)
DEFAULT_COUPLED_CONTROL_DIR = DEFAULT_COUPLED_OUTPUT_DIR

PILOT_SPECS = ((0.030, 5), (0.050, 4))
INTEGRATION_GRIDS = (801, 1601)
MODE_RESIDUAL_TOL = min(JOINT_KINEMATIC_TOL, JOINT_FORCE_TOL)
INTEGRATION_TOL = ENERGY_SUM_TOL
COUPLED_MIN_BETA_DEG = 15.0
COUPLED_SELECTED_BETA_DEG = 15.0
THINNESS_LIMIT = 0.1
CONTRAST_NUMERICAL_FLOOR = 1.0e-30
ADDITIONAL_ANGLE_EPSILON = 0.030
ADDITIONAL_ANGLE_CENTRAL_INDEX = 5
ADDITIONAL_ANGLE_MIN_BETA_DEG = 30.0

CANDIDATE_FIELDS = (
    "epsilon_0",
    "case_id",
    "beta_deg",
    "mu",
    "eta",
    "central_sorted_index",
    "delta_minus",
    "delta_center",
    "delta_plus",
    "R_k",
    "D_k",
    "selected",
    "s_max",
    "thinness_flag",
)

ENERGY_FIELDS = (
    "epsilon_0",
    "case_id",
    "beta_deg",
    "mu",
    "eta",
    "sorted_index",
    "role",
    "Lambda_EB",
    "Lambda_Timo",
    "delta_f",
    "U_axial_Timo",
    "U_bending_Timo",
    "U_shear_Timo",
    "U_total_Timo",
    "chi_axial_Timo",
    "chi_bending_Timo",
    "chi_shear_Timo",
    "chi_bending_shear_Timo",
    "chi_axial_EB",
    "s_max",
    "thinness_flag",
    "article_ready_status",
)

CONVERGENCE_FIELDS = (
    "epsilon_0",
    "case_id",
    "sorted_index",
    "role",
    "n_points_coarse",
    "n_points_fine",
    "chi_axial_801",
    "chi_bending_801",
    "chi_shear_801",
    "chi_axial_1601",
    "chi_bending_1601",
    "chi_shear_1601",
    "delta_chi_a",
    "delta_chi_b",
    "delta_chi_s",
    "integration_tolerance",
    "passed",
)

RESIDUAL_FIELDS = (
    "epsilon_0",
    "case_id",
    "sorted_index",
    "role",
    "Lambda_Timo",
    "null_vector_source",
    "timo_mode_coefficients_smallest_singular_value",
    "timo_mode_coefficients_relative_residual_on_certificate_scaling",
    "smallest_singular_value",
    "second_smallest_singular_value",
    "singular_value_ratio",
    "relative_null_residual",
    "rod1_clamp_u",
    "rod1_clamp_w",
    "rod1_clamp_psi",
    "rod2_clamp_u",
    "rod2_clamp_w",
    "rod2_clamp_psi",
    "max_abs_clamp_gap",
    "gap_w",
    "gap_u",
    "gap_psi",
    "gap_M",
    "gap_Q",
    "gap_N",
    "max_abs_kinematic_gap",
    "max_abs_force_gap",
    "sign_flip_delta_chi_a",
    "sign_flip_delta_chi_b",
    "sign_flip_delta_chi_s",
    "sign_flip_max_abs_delta",
    "mode_residual_tolerance",
    "energy_tolerance",
    "svd_pass",
    "null_residual_pass",
    "clamp_pass",
    "joint_kinematic_pass",
    "joint_force_pass",
    "sign_flip_pass",
    "passed",
    "warnings",
)

OPERATION_COUNTER_FIELDS = (
    "root_solver_calls",
    "strict_solver_calls",
    "family_detector_calls",
    "local_repair_calls",
    "matrix_evaluator_calls",
    "null_vector_reconstructions",
    "energy_integral_evaluations",
    "shape_plot_calls",
    "MAC_calls",
    "FEM_calls",
    "anisotropic_calls",
)

COUPLED_SELECTION_FIELDS = (
    "epsilon_0",
    "case_id",
    "beta_deg",
    "mu",
    "eta",
    "central_sorted_index",
    "delta_minus",
    "delta_center",
    "delta_plus",
    "R_k",
    "D_k",
    "s_max",
    "thinness_flag",
    "s_max_le_0p1",
    "selection_pool",
    "selected_for_energy",
    "selection_reason",
)

COUPLED_ENERGY_FIELDS = tuple(field for field in ENERGY_FIELDS if field != "chi_axial_EB")

COUPLED_SUMMARY_FIELDS = (
    "example_role",
    "epsilon_0",
    "case_id",
    "beta_deg",
    "mu",
    "eta",
    "central_sorted_index",
    "delta_minus",
    "delta_center",
    "delta_plus",
    "R_k",
    "D_k",
    "chi_axial_left",
    "chi_axial_center",
    "chi_axial_right",
    "chi_bending_shear_left",
    "chi_bending_shear_center",
    "chi_bending_shear_right",
    "Q_left",
    "Q_right",
    "contrast_numerical_floor",
    "s_max",
    "thinness_flag",
    "center_exceeds_both_neighbors",
    "notes",
)

ADDITIONAL_SELECTION_FIELDS = (
    "epsilon_0",
    "case_id",
    "beta_deg",
    "mu",
    "eta",
    "central_sorted_index",
    "delta_minus",
    "delta_center",
    "delta_plus",
    "R_k",
    "D_k",
    "s_max",
    "thinness_flag",
    "selection_pool",
    "selected_for_energy",
    "selection_reason",
)

ANGLE_COMPARISON_FIELDS = (
    "epsilon_0",
    "mu",
    "eta",
    "case_id_beta15",
    "case_id_beta30",
    "beta15_deg",
    "beta30_deg",
    "sorted_index",
    "delta_f_beta15",
    "delta_f_beta30",
    "delta_f_change_beta30_minus_beta15",
    "chi_axial_beta15",
    "chi_axial_beta30",
    "chi_axial_change_beta30_minus_beta15",
    "chi_bending_shear_beta15",
    "chi_bending_shear_beta30",
    "chi_bending_shear_change_beta30_minus_beta15",
    "Q_left_beta15",
    "Q_right_beta15",
    "Q_left_beta30",
    "Q_right_beta30",
    "Q_left_change_beta30_minus_beta15",
    "Q_right_change_beta30_minus_beta15",
)


def operation_counters() -> dict[str, int]:
    return {field: 0 for field in OPERATION_COUNTER_FIELDS}


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _atomic_csv(
    path: Path,
    fields: Sequence[str],
    rows: Iterable[Mapping[str, object]],
    *,
    delimiter: str = ",",
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter=delimiter,
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)
    return path


def _atomic_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    os.replace(temporary, path)
    return path


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> Path:
    return _atomic_text(
        path,
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True, default=str) + "\n",
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _directory_hashes(directory: Path) -> dict[str, str]:
    directory = directory.resolve()
    if not directory.is_dir():
        raise FileNotFoundError(f"immutable control directory is missing: {directory}")
    return {
        path.relative_to(directory).as_posix(): _sha256(path)
        for path in sorted(directory.rglob("*"))
        if path.is_file()
    }


def _resolve_repo_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path


def _same_float(left: str | float, right: float, tolerance: float = 1.0e-12) -> bool:
    return abs(float(left) - float(right)) <= tolerance


def _bool_text(value: object) -> str:
    return "true" if bool(value) else "false"


def _load_certificate(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("scientific_scope") != SCIENTIFIC_SCOPE:
        raise ValueError(f"mixed scientific scope in compact certificate: {path}")
    if payload.get("schema_version") != CERTIFICATE_SCHEMA:
        raise ValueError(f"unexpected compact certificate schema: {path}")
    return payload


def _target_level_row(rows: Sequence[dict[str, str]], epsilon: float) -> dict[str, str]:
    matches = [row for row in rows if _same_float(row["epsilon_0"], epsilon)]
    if len(matches) != 1:
        raise ValueError(f"epsilon_0={epsilon:.3f}: expected one level summary, found {len(matches)}")
    return matches[0]


def _quality_is_resolved(model_payload: Mapping[str, Any], ranks: set[int]) -> bool:
    quality_by_rank = {
        int(row["rank"]): row for row in model_payload.get("quality_by_rank", [])
    }
    for rank in ranks:
        row = quality_by_rank.get(rank)
        if row is None or row.get("status") != "matrix_confirmed_resolved":
            return False
        if str(row.get("cluster_id", "")).strip():
            return False
    clustered_ranks = {
        int(row["sorted_index"])
        for row in model_payload.get("clusters", [])
        if str(row.get("sorted_index", "")).strip()
    }
    return not bool(clustered_ranks & ranks)


def _candidate_output_row(record: Mapping[str, Any], *, selected: bool) -> dict[str, object]:
    return {
        "epsilon_0": float(record["epsilon_0"]),
        "case_id": record["case_id"],
        "beta_deg": float(record["beta_deg"]),
        "mu": float(record["mu"]),
        "eta": float(record["eta"]),
        "central_sorted_index": int(record["central_sorted_index"]),
        "delta_minus": float(record["delta_minus"]),
        "delta_center": float(record["delta_center"]),
        "delta_plus": float(record["delta_plus"]),
        "R_k": float(record["R_k"]),
        "D_k": float(record["D_k"]),
        "selected": _bool_text(selected),
        "s_max": float(record["s_max"]),
        "thinness_flag": _bool_text(record["thinness_flag"]),
    }


def select_candidates(
    finalization_dir: Path,
    compact_dir: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, str]]:
    """Select the two geometries using compact frequency data only."""

    finalization_dir = finalization_dir.resolve()
    compact_dir = compact_dir.resolve()
    source_paths = (
        finalization_dir / "final_case_certificates.csv",
        finalization_dir / "epsilon_level_summary.csv",
        finalization_dir / "unresolved_cases.csv",
        compact_dir / "compact_index.csv",
    )
    missing = [str(path) for path in source_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing compact pilot inputs: " + ", ".join(missing))

    final_rows = _read_csv(source_paths[0])
    level_rows = _read_csv(source_paths[1])
    unresolved_rows = _read_csv(source_paths[2])
    compact_rows = _read_csv(source_paths[3])
    if {row["scientific_scope"] for row in final_rows} != {SCIENTIFIC_SCOPE}:
        raise ValueError("final case table contains a scientific scope outside the pilot")

    unresolved_ids = {row["case_id"] for row in unresolved_rows}
    final_by_id = {row["case_id"]: row for row in final_rows}
    compact_by_id = {row["case_id"]: row for row in compact_rows}
    if len(final_by_id) != len(final_rows) or len(compact_by_id) != len(compact_rows):
        raise ValueError("duplicate case ID in compact inputs")

    for epsilon, _central_index in PILOT_SPECS:
        level = _target_level_row(level_rows, epsilon)
        if str(level.get("envelope_finalizable", "")).lower() != "true":
            raise ValueError(f"epsilon_0={epsilon:.3f}: finalized level is not finalizable")
        main_cases = [row for row in final_rows if _same_float(row["epsilon_0"], epsilon)]
        if len(main_cases) != 194:
            raise ValueError(
                f"epsilon_0={epsilon:.3f}: expected 194 main-grid cases, found {len(main_cases)}"
            )

    input_hashes = {path.resolve().as_posix(): _sha256(path) for path in source_paths}
    records: list[dict[str, Any]] = []
    certificates_read: list[Path] = []

    for epsilon, central_index in PILOT_SPECS:
        ranks = {central_index - 1, central_index, central_index + 1}
        level_records: list[dict[str, Any]] = []
        level_final_rows = [row for row in final_rows if _same_float(row["epsilon_0"], epsilon)]
        for final_row in level_final_rows:
            case_id = final_row["case_id"]
            if case_id in unresolved_ids:
                continue
            if not final_row.get("N_true", "").strip():
                continue
            if not final_row.get("final_execution_status", "").startswith("resolved"):
                continue
            if final_row.get("unresolved_reason", "").strip():
                continue

            compact_row = compact_by_id.get(case_id)
            if compact_row is None:
                continue
            if compact_row.get("scientific_scope") != SCIENTIFIC_SCOPE:
                continue
            if compact_row.get("n_true_status") != "exact":
                continue
            if compact_row.get("required_guard_confirmed", "").lower() != "true":
                continue
            certificate_path = _resolve_repo_path(compact_row["certificate_path"])
            certificate = _load_certificate(certificate_path)
            certificates_read.append(certificate_path)
            result = certificate.get("result", {})
            if result.get("n_true_status") != "exact":
                continue
            if not bool(result.get("required_guard_confirmed")):
                continue
            if str(result.get("unresolved_reason", "")).strip():
                continue
            if not bool(certificate.get("diagnostics", {}).get("unresolved_cluster_below_guard")):
                continue

            spectra = certificate.get("spectra", {})
            eb = spectra.get("Euler-Bernoulli", {})
            timo = spectra.get("Timoshenko", {})
            if not _quality_is_resolved(eb, ranks) or not _quality_is_resolved(timo, ranks):
                continue
            eb_roots = [float(value) for value in eb.get("roots", [])]
            timo_roots = [float(value) for value in timo.get("roots", [])]
            if len(eb_roots) < max(ranks) or len(timo_roots) < max(ranks):
                continue
            roots_needed = [eb_roots[rank - 1] for rank in sorted(ranks)] + [
                timo_roots[rank - 1] for rank in sorted(ranks)
            ]
            if not all(math.isfinite(value) and value > 0.0 for value in roots_needed):
                continue
            delta_payload = spectra.get("delta_f", {})
            try:
                delta_minus = float(delta_payload[str(central_index - 1)])
                delta_center = float(delta_payload[str(central_index)])
                delta_plus = float(delta_payload[str(central_index + 1)])
            except (KeyError, TypeError, ValueError):
                continue
            if not all(
                math.isfinite(value) and value >= 0.0
                for value in (delta_minus, delta_center, delta_plus)
            ):
                continue
            if not (delta_center < delta_minus and delta_center < delta_plus):
                continue

            neighbor_mean = 0.5 * (delta_minus + delta_plus)
            if neighbor_mean <= 0.0:
                continue
            geometry = certificate.get("geometry", {})
            if any(
                not _same_float(geometry[name], final_row[source_name])
                for name, source_name in (
                    ("epsilon_0", "epsilon_0"),
                    ("beta_deg", "beta"),
                    ("mu", "mu"),
                    ("eta", "eta"),
                )
            ):
                raise ValueError(f"geometry mismatch for {case_id}")
            s_max = float(final_row["s_max"])
            level_records.append(
                {
                    "epsilon_0": epsilon,
                    "case_id": case_id,
                    "beta_deg": float(final_row["beta"]),
                    "mu": float(final_row["mu"]),
                    "eta": float(final_row["eta"]),
                    "central_sorted_index": central_index,
                    "delta_minus": delta_minus,
                    "delta_center": delta_center,
                    "delta_plus": delta_plus,
                    "R_k": delta_center / neighbor_mean,
                    "D_k": neighbor_mean - delta_center,
                    "s_max": s_max,
                    "thinness_flag": s_max > 0.1,
                    "certificate_path": certificate_path,
                    "certificate": certificate,
                }
            )

        level_records.sort(
            key=lambda row: (
                float(row["R_k"]),
                -float(row["D_k"]),
                float(row["beta_deg"]),
                float(row["mu"]),
                float(row["eta"]),
                str(row["case_id"]),
            )
        )
        records.extend(level_records)

    selected: list[dict[str, Any]] = []
    for epsilon, _central_index in PILOT_SPECS:
        matches = [row for row in records if _same_float(row["epsilon_0"], epsilon)]
        if matches:
            selected.append(matches[0])

    for path in sorted(set(certificates_read)):
        input_hashes[path.resolve().as_posix()] = _sha256(path)
    return records, selected, input_hashes


def write_selection_outputs(
    output_dir: Path,
    records: Sequence[Mapping[str, Any]],
    selected: Sequence[Mapping[str, Any]],
) -> None:
    selected_ids = {str(row["case_id"]) for row in selected}
    ranking_rows = [
        _candidate_output_row(row, selected=str(row["case_id"]) in selected_ids)
        for row in records
    ]
    top_ten_rows: list[dict[str, object]] = []
    for epsilon, _central_index in PILOT_SPECS:
        level = [row for row in ranking_rows if _same_float(row["epsilon_0"], epsilon)]
        top_ten_rows.extend(level[:10])
    selected_rows = [
        _candidate_output_row(row, selected=True)
        for row in selected
    ]
    _atomic_csv(output_dir / "candidate_ranking.csv", CANDIDATE_FIELDS, ranking_rows)
    _atomic_csv(output_dir / "candidate_ranking_top10.csv", CANDIDATE_FIELDS, top_ten_rows)
    _atomic_csv(output_dir / "selected_cases.csv", CANDIDATE_FIELDS, selected_rows)


def _candidate_sort_key(row: Mapping[str, Any]) -> tuple[float, float, float, float, float, str]:
    return (
        float(row["R_k"]),
        -float(row["D_k"]),
        float(row["beta_deg"]),
        float(row["mu"]),
        float(row["eta"]),
        str(row["case_id"]),
    )


def _coupled_selection_output_row(record: Mapping[str, Any]) -> dict[str, object]:
    return {
        "epsilon_0": float(record["epsilon_0"]),
        "case_id": record["case_id"],
        "beta_deg": float(record["beta_deg"]),
        "mu": float(record["mu"]),
        "eta": float(record["eta"]),
        "central_sorted_index": int(record["central_sorted_index"]),
        "delta_minus": float(record["delta_minus"]),
        "delta_center": float(record["delta_center"]),
        "delta_plus": float(record["delta_plus"]),
        "R_k": float(record["R_k"]),
        "D_k": float(record["D_k"]),
        "s_max": float(record["s_max"]),
        "thinness_flag": _bool_text(record["thinness_flag"]),
        "s_max_le_0p1": _bool_text(float(record["s_max"]) <= THINNESS_LIMIT),
        "selection_pool": record["selection_pool"],
        "selected_for_energy": _bool_text(record["selected_for_energy"]),
        "selection_reason": record["selection_reason"],
    }


def select_coupled_candidates(
    finalization_dir: Path,
    compact_dir: Path,
    beta0_control_dir: Path,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, str],
    dict[str, object],
]:
    """Build nonzero-angle pools from frequency-only finalized evidence."""

    beta0_control_dir = beta0_control_dir.resolve()
    ranking_path = beta0_control_dir / "candidate_ranking.csv"
    if not ranking_path.is_file():
        raise FileNotFoundError(f"beta=0 pilot candidate ranking is missing: {ranking_path}")

    regenerated, _old_selected, input_hashes = select_candidates(
        finalization_dir,
        compact_dir,
    )
    regenerated_by_id = {str(row["case_id"]): row for row in regenerated}
    ranking_rows = _read_csv(ranking_path)
    if len(ranking_rows) != len(regenerated_by_id):
        raise ValueError("immutable candidate ranking does not match compact frequency candidates")
    metric_fields = (
        "epsilon_0",
        "beta_deg",
        "mu",
        "eta",
        "central_sorted_index",
        "delta_minus",
        "delta_center",
        "delta_plus",
        "R_k",
        "D_k",
        "s_max",
    )
    ranked_records: list[dict[str, Any]] = []
    for ranking_row in ranking_rows:
        case_id = ranking_row["case_id"]
        regenerated_row = regenerated_by_id.get(case_id)
        if regenerated_row is None:
            raise ValueError(f"candidate ranking contains a non-finalized case: {case_id}")
        for field in metric_fields:
            if not math.isclose(
                float(ranking_row[field]),
                float(regenerated_row[field]),
                rel_tol=1.0e-12,
                abs_tol=1.0e-15,
            ):
                raise ValueError(f"candidate ranking mismatch for {case_id}: {field}")
        ranked_records.append(regenerated_row)

    input_hashes[ranking_path.resolve().as_posix()] = _sha256(ranking_path)
    nonzero_records = [
        row for row in ranked_records if float(row["beta_deg"]) >= COUPLED_MIN_BETA_DEG
    ]
    all_by_epsilon = {
        epsilon: sorted(
            (row for row in nonzero_records if _same_float(row["epsilon_0"], epsilon)),
            key=_candidate_sort_key,
        )
        for epsilon, _central_index in PILOT_SPECS
    }
    thin_003 = [
        row for row in all_by_epsilon[0.030] if float(row["s_max"]) <= THINNESS_LIMIT
    ]
    if not thin_003:
        raise ValueError("epsilon_0=0.030 has no beta>=15 candidate with s_max<=0.1")
    if not all_by_epsilon[0.050]:
        raise ValueError("epsilon_0=0.050 has no beta>=15 strict-local-minimum candidate")

    final_rows = _read_csv(finalization_dir.resolve() / "final_case_certificates.csv")
    main_005_rows = [row for row in final_rows if _same_float(row["epsilon_0"], 0.050)]
    if len(main_005_rows) != 194:
        raise ValueError("epsilon_0=0.050 main-grid inventory must contain 194 cases")
    thin_main_005 = [row for row in main_005_rows if float(row["s_max"]) <= THINNESS_LIMIT]
    if thin_main_005:
        raise ValueError(
            "epsilon_0=0.050 unexpectedly contains main-grid geometries with s_max<=0.1"
        )

    selected_003 = dict(thin_003[0])
    selected_003.update(
        {
            "selection_pool": "nonzero_angle_smax_le_0p1",
            "selected_for_energy": True,
            "selection_reason": "minimum_R_k_then_maximum_D_k_in_nonzero_angle_smax_le_0p1",
        }
    )
    selected_005 = dict(all_by_epsilon[0.050][0])
    selected_005.update(
        {
            "selection_pool": "nonzero_angle_all",
            "selected_for_energy": True,
            "selection_reason": (
                "no_epsilon_0_0p050_main_grid_case_with_s_max_le_0p1;"
                "minimum_R_k_then_maximum_D_k_in_nonzero_angle_all"
            ),
        }
    )
    selected = [selected_003, selected_005]
    selected_keys = {
        (str(row["case_id"]), str(row["selection_pool"])) for row in selected
    }

    pool_records: list[dict[str, Any]] = []
    pool_specs = (
        (0.030, "nonzero_angle_all", all_by_epsilon[0.030]),
        (0.030, "nonzero_angle_smax_le_0p1", thin_003),
        (0.050, "nonzero_angle_all", all_by_epsilon[0.050]),
    )
    for epsilon, pool_name, pool in pool_specs:
        for record in pool:
            item = dict(record)
            key = (str(item["case_id"]), pool_name)
            item["selection_pool"] = pool_name
            item["selected_for_energy"] = key in selected_keys
            if key in selected_keys:
                selected_item = next(row for row in selected if str(row["case_id"]) == key[0])
                item["selection_reason"] = selected_item["selection_reason"]
            elif (
                _same_float(epsilon, 0.030)
                and pool_name == "nonzero_angle_all"
                and str(item["case_id"]) == str(selected_003["case_id"])
            ):
                item["selection_reason"] = (
                    "eligible_in_nonzero_angle_all;energy_selection_uses_smax_le_0p1_pool"
                )
            else:
                item["selection_reason"] = "eligible_not_selected"
            pool_records.append(item)

    context: dict[str, object] = {
        "nonzero_angle_all_candidate_count_epsilon_0p030": len(all_by_epsilon[0.030]),
        "nonzero_angle_smax_le_0p1_candidate_count_epsilon_0p030": len(thin_003),
        "nonzero_angle_all_candidate_count_epsilon_0p050": len(all_by_epsilon[0.050]),
        "main_grid_smax_le_0p1_count_epsilon_0p050": len(thin_main_005),
        "selection_uses_energy_values": False,
    }
    return pool_records, selected, input_hashes, context


def write_coupled_selection_outputs(
    output_dir: Path,
    pool_records: Sequence[Mapping[str, Any]],
    selected: Sequence[Mapping[str, Any]],
) -> None:
    ranking_rows = [_coupled_selection_output_row(row) for row in pool_records]
    top_ten_rows: list[dict[str, object]] = []
    for epsilon, pool_name in (
        (0.030, "nonzero_angle_all"),
        (0.030, "nonzero_angle_smax_le_0p1"),
        (0.050, "nonzero_angle_all"),
    ):
        rows = [
            row
            for row in ranking_rows
            if _same_float(row["epsilon_0"], epsilon) and row["selection_pool"] == pool_name
        ]
        top_ten_rows.extend(rows[:10])
    _atomic_csv(
        output_dir / "coupled_candidate_ranking.csv",
        COUPLED_SELECTION_FIELDS,
        ranking_rows,
    )
    _atomic_csv(
        output_dir / "coupled_candidate_ranking_top10.csv",
        COUPLED_SELECTION_FIELDS,
        top_ten_rows,
    )
    _atomic_csv(
        output_dir / "selected_coupled_cases.csv",
        COUPLED_SELECTION_FIELDS,
        (_coupled_selection_output_row(row) for row in selected),
    )


def _additional_selection_output_row(record: Mapping[str, Any]) -> dict[str, object]:
    return {
        "epsilon_0": float(record["epsilon_0"]),
        "case_id": record["case_id"],
        "beta_deg": float(record["beta_deg"]),
        "mu": float(record["mu"]),
        "eta": float(record["eta"]),
        "central_sorted_index": int(record["central_sorted_index"]),
        "delta_minus": float(record["delta_minus"]),
        "delta_center": float(record["delta_center"]),
        "delta_plus": float(record["delta_plus"]),
        "R_k": float(record["R_k"]),
        "D_k": float(record["D_k"]),
        "s_max": float(record["s_max"]),
        "thinness_flag": _bool_text(record["thinness_flag"]),
        "selection_pool": record["selection_pool"],
        "selected_for_energy": _bool_text(record["selected_for_energy"]),
        "selection_reason": record["selection_reason"],
    }


def select_additional_angle_candidate(
    finalization_dir: Path,
    compact_dir: Path,
    coupled_control_dir: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any], dict[str, str], dict[str, object]]:
    """Select the beta>=30 extension after case-ID deduplication, without energy data."""

    coupled_control_dir = coupled_control_dir.resolve()
    ranking_path = coupled_control_dir / "coupled_candidate_ranking.csv"
    if not ranking_path.is_file():
        raise FileNotFoundError(f"coupled candidate ranking is missing: {ranking_path}")

    regenerated, _selected, input_hashes = select_candidates(
        finalization_dir,
        compact_dir,
    )
    regenerated_by_id = {str(row["case_id"]): row for row in regenerated}
    ranking_rows = _read_csv(ranking_path)
    metric_fields = (
        "epsilon_0",
        "beta_deg",
        "mu",
        "eta",
        "central_sorted_index",
        "delta_minus",
        "delta_center",
        "delta_plus",
        "R_k",
        "D_k",
        "s_max",
    )
    deduplicated: dict[str, dict[str, str]] = {}
    for ranking_row in ranking_rows:
        case_id = str(ranking_row["case_id"])
        previous = deduplicated.get(case_id)
        if previous is not None:
            if any(
                not math.isclose(
                    float(previous[field]),
                    float(ranking_row[field]),
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-15,
                )
                for field in metric_fields
            ):
                raise ValueError(f"conflicting duplicate ranking rows for {case_id}")
            continue
        deduplicated[case_id] = ranking_row

    validated: list[dict[str, Any]] = []
    for case_id, ranking_row in deduplicated.items():
        regenerated_row = regenerated_by_id.get(case_id)
        if regenerated_row is None:
            raise ValueError(f"additional-angle ranking case is not finalized: {case_id}")
        for field in metric_fields:
            if not math.isclose(
                float(ranking_row[field]),
                float(regenerated_row[field]),
                rel_tol=1.0e-12,
                abs_tol=1.0e-15,
            ):
                raise ValueError(f"additional-angle ranking mismatch for {case_id}: {field}")
        validated.append(regenerated_row)

    pool = [
        row
        for row in validated
        if _same_float(row["epsilon_0"], ADDITIONAL_ANGLE_EPSILON)
        and int(row["central_sorted_index"]) == ADDITIONAL_ANGLE_CENTRAL_INDEX
        and float(row["beta_deg"]) >= ADDITIONAL_ANGLE_MIN_BETA_DEG
        and float(row["s_max"]) <= THINNESS_LIMIT
        and float(row["delta_center"]) < float(row["delta_minus"])
        and float(row["delta_center"]) < float(row["delta_plus"])
    ]
    pool.sort(key=_candidate_sort_key)
    if not pool:
        raise ValueError("no epsilon_0=0.030, beta>=30, s_max<=0.1 candidate exists")

    pool_name = "epsilon_0p030_beta_ge_30_smax_le_0p1"
    selection_reason = "minimum_R_5_then_maximum_D_5_then_beta_mu_eta_case_id"
    selected = dict(pool[0])
    selected.update(
        {
            "selection_pool": pool_name,
            "selected_for_energy": True,
            "selection_reason": selection_reason,
        }
    )
    selected_id = str(selected["case_id"])
    candidates: list[dict[str, Any]] = []
    for row in pool:
        item = dict(row)
        is_selected = str(item["case_id"]) == selected_id
        item.update(
            {
                "selection_pool": pool_name,
                "selected_for_energy": is_selected,
                "selection_reason": selection_reason if is_selected else "eligible_not_selected",
            }
        )
        candidates.append(item)

    input_hashes[ranking_path.resolve().as_posix()] = _sha256(ranking_path)
    context: dict[str, object] = {
        "ranking_source_row_count": len(ranking_rows),
        "ranking_unique_case_id_count": len(deduplicated),
        "ranking_duplicate_rows_removed": len(ranking_rows) - len(deduplicated),
        "additional_angle_candidate_count": len(candidates),
        "selection_pool": pool_name,
        "selection_uses_energy_values": False,
    }
    return candidates, selected, input_hashes, context


def write_additional_angle_selection_outputs(
    output_dir: Path,
    candidates: Sequence[Mapping[str, Any]],
    selected: Mapping[str, Any],
) -> None:
    _atomic_csv(
        output_dir / "additional_angle_candidate_ranking.csv",
        ADDITIONAL_SELECTION_FIELDS,
        (_additional_selection_output_row(row) for row in candidates),
    )
    _atomic_csv(
        output_dir / "selected_additional_angle_case.csv",
        ADDITIONAL_SELECTION_FIELDS,
        [_additional_selection_output_row(selected)],
    )


def _energy_fractions(energy: Mapping[str, float]) -> tuple[float, float, float]:
    return (
        float(energy["axial_fraction"]),
        float(energy["bending_fraction"]),
        float(energy["shear_fraction"]),
    )


def _role(index: int, central_index: int) -> str:
    return {
        central_index - 1: "left_neighbor",
        central_index: "central_dip",
        central_index + 1: "right_neighbor",
    }[index]


def reconstruct_selected_modes(
    selected: Sequence[Mapping[str, Any]],
    counters: dict[str, int],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    energy_rows: list[dict[str, object]] = []
    convergence_rows: list[dict[str, object]] = []
    residual_rows: list[dict[str, object]] = []

    for selected_case in selected:
        epsilon = float(selected_case["epsilon_0"])
        beta_deg = float(selected_case["beta_deg"])
        mu = float(selected_case["mu"])
        eta = float(selected_case["eta"])
        central_index = int(selected_case["central_sorted_index"])
        certificate = selected_case["certificate"]
        spectra = certificate["spectra"]
        eb_roots = [float(value) for value in spectra["Euler-Bernoulli"]["roots"]]
        timo_roots = [float(value) for value in spectra["Timoshenko"]["roots"]]

        for sorted_index in range(central_index - 1, central_index + 2):
            role = _role(sorted_index, central_index)
            lambda_eb = eb_roots[sorted_index - 1]
            lambda_timo = timo_roots[sorted_index - 1]
            delta_f = float(spectra["delta_f"][str(sorted_index)])

            mode = TIMO.timo_mode_coefficients(
                lambda_timo,
                beta_deg,
                mu,
                epsilon,
                eta,
            )
            counters["matrix_evaluator_calls"] += 1
            counters["null_vector_reconstructions"] += 1

            # Compact roots were certified with the existing floored diagnostic
            # row scaling.  At an exact straight-geometry axial root a
            # structurally zero row may contain only round-off.  The older
            # timo_mode_coefficients scaling divides that row by its own tiny
            # norm, so retain its vector when it is valid on the certificate
            # scaling and use the certificate-consistent SVD only when needed.
            matrix, matrix_warnings = TIMO.timo_coupling_matrix(
                lambda_timo,
                beta_deg,
                mu,
                epsilon,
                eta,
            )
            counters["matrix_evaluator_calls"] += 1
            certificate_scaled_matrix = COMPLETE.diagnostic_scaled_matrix(matrix)
            singular_values = np.linalg.svd(certificate_scaled_matrix, compute_uv=False)
            mode_coeff = np.asarray(mode.coeff, dtype=float)
            matrix_norm = float(np.linalg.norm(certificate_scaled_matrix, ord=2))
            mode_denominator = float(matrix_norm * np.linalg.norm(mode_coeff))
            mode_relative_residual = (
                float(np.linalg.norm(certificate_scaled_matrix @ mode_coeff)) / mode_denominator
                if mode_denominator > 1.0e-30
                else float("nan")
            )
            if math.isfinite(mode_relative_residual) and mode_relative_residual <= MODE_RESIDUAL_TOL:
                coeff = mode_coeff
                null_vector_source = "timo_mode_coefficients"
            else:
                _u, _singular_values, vh = np.linalg.svd(
                    certificate_scaled_matrix,
                    full_matrices=True,
                )
                coeff = np.asarray(vh[-1, :], dtype=float)
                coeff_norm = float(np.linalg.norm(coeff))
                if coeff_norm <= 0.0 or not math.isfinite(coeff_norm):
                    raise ValueError("certificate-scaled Timoshenko null vector is not finite")
                coeff = coeff / coeff_norm
                pivot = int(np.argmax(np.abs(coeff)))
                if coeff[pivot] < 0.0:
                    coeff = -coeff
                counters["null_vector_reconstructions"] += 1
                null_vector_source = "certificate_floored_row_scaling_fallback"

            energies: dict[int, dict[str, float]] = {}
            for n_points in INTEGRATION_GRIDS:
                energies[n_points] = TIMO.timo_energy_partition(
                    lambda_timo,
                    beta_deg,
                    mu,
                    epsilon,
                    eta,
                    coeff=coeff,
                    n_points=n_points,
                )
                counters["energy_integral_evaluations"] += 1

            sign_energy = TIMO.timo_energy_partition(
                lambda_timo,
                beta_deg,
                mu,
                epsilon,
                eta,
                coeff=-coeff,
                n_points=INTEGRATION_GRIDS[0],
            )
            counters["energy_integral_evaluations"] += 1

            coarse = energies[INTEGRATION_GRIDS[0]]
            fine = energies[INTEGRATION_GRIDS[1]]
            coarse_fractions = _energy_fractions(coarse)
            fine_fractions = _energy_fractions(fine)
            sign_fractions = _energy_fractions(sign_energy)
            convergence_deltas = tuple(
                abs(fine_value - coarse_value)
                for fine_value, coarse_value in zip(fine_fractions, coarse_fractions)
            )
            sign_deltas = tuple(
                abs(sign_value - coarse_value)
                for sign_value, coarse_value in zip(sign_fractions, coarse_fractions)
            )
            convergence_pass = all(value <= INTEGRATION_TOL for value in convergence_deltas)
            sign_flip_pass = all(value <= ENERGY_SUM_TOL for value in sign_deltas)

            fields = TIMO.timo_mode_fields(
                lambda_timo,
                beta_deg,
                mu,
                epsilon,
                eta,
                coeff=coeff,
                n_points=INTEGRATION_GRIDS[1],
            )
            denominator = float(matrix_norm * np.linalg.norm(coeff))
            relative_null_residual = (
                float(np.linalg.norm(certificate_scaled_matrix @ coeff)) / denominator
                if denominator > 1.0e-30
                else float("nan")
            )
            smallest = float(singular_values[-1])
            second_smallest = float(singular_values[-2])
            singular_ratio = (
                smallest / second_smallest if second_smallest > 1.0e-30 else float("nan")
            )

            rod1 = fields["rod1"]
            rod2 = fields["rod2"]
            clamp_values = (
                float(np.asarray(rod1["u"], dtype=float)[0]),
                float(np.asarray(rod1["w"], dtype=float)[0]),
                float(np.asarray(rod1["psi"], dtype=float)[0]),
                float(np.asarray(rod2["u"], dtype=float)[0]),
                float(np.asarray(rod2["w"], dtype=float)[0]),
                float(np.asarray(rod2["psi"], dtype=float)[0]),
            )
            max_clamp = max(abs(value) for value in clamp_values)
            joint = timo_joint_continuity_row(
                epsilon=epsilon,
                beta_deg=beta_deg,
                eta=eta,
                mu=mu,
                sorted_index=sorted_index,
                Lambda=lambda_timo,
                fields=fields,
            )

            fine_sum_error = abs(sum(fine_fractions) - 1.0)
            coarse_sum_error = abs(sum(coarse_fractions) - 1.0)
            finite_energy = all(
                math.isfinite(float(value))
                for value in (
                    fine["U_a_total"],
                    fine["U_b_total"],
                    fine["U_s_total"],
                    fine["U_total"],
                    *fine_fractions,
                )
            ) and float(fine["U_total"]) > 0.0
            energy_sum_pass = (
                fine_sum_error <= ENERGY_SUM_TOL and coarse_sum_error <= ENERGY_SUM_TOL
            )
            svd_pass = math.isfinite(smallest) and smallest <= MODE_RESIDUAL_TOL
            null_pass = (
                math.isfinite(relative_null_residual)
                and relative_null_residual <= MODE_RESIDUAL_TOL
            )
            clamp_pass = max_clamp <= MODE_RESIDUAL_TOL
            joint_kinematic_pass = bool(joint["pass_kinematic"])
            joint_force_pass = bool(joint["pass_force"])
            residual_pass = all(
                (
                    svd_pass,
                    null_pass,
                    clamp_pass,
                    joint_kinematic_pass,
                    joint_force_pass,
                    sign_flip_pass,
                )
            )
            ready = all(
                (
                    finite_energy,
                    energy_sum_pass,
                    convergence_pass,
                    residual_pass,
                )
            )
            readiness_reasons: list[str] = []
            if not finite_energy:
                readiness_reasons.append("nonfinite_or_nonpositive_energy")
            if not energy_sum_pass:
                readiness_reasons.append("energy_fraction_sum_failed")
            if not convergence_pass:
                readiness_reasons.append("integration_convergence_failed")
            if not residual_pass:
                readiness_reasons.append("mode_residual_or_sign_check_failed")
            status = "article_ready" if ready else "not_article_ready:" + ";".join(readiness_reasons)

            chi_a, chi_b, chi_s = fine_fractions
            energy_rows.append(
                {
                    "epsilon_0": epsilon,
                    "case_id": selected_case["case_id"],
                    "beta_deg": beta_deg,
                    "mu": mu,
                    "eta": eta,
                    "sorted_index": sorted_index,
                    "role": role,
                    "Lambda_EB": lambda_eb,
                    "Lambda_Timo": lambda_timo,
                    "delta_f": delta_f,
                    "U_axial_Timo": float(fine["U_a_total"]),
                    "U_bending_Timo": float(fine["U_b_total"]),
                    "U_shear_Timo": float(fine["U_s_total"]),
                    "U_total_Timo": float(fine["U_total"]),
                    "chi_axial_Timo": chi_a,
                    "chi_bending_Timo": chi_b,
                    "chi_shear_Timo": chi_s,
                    "chi_bending_shear_Timo": chi_b + chi_s,
                    "chi_axial_EB": "",
                    "s_max": float(selected_case["s_max"]),
                    "thinness_flag": _bool_text(selected_case["thinness_flag"]),
                    "article_ready_status": status,
                }
            )
            convergence_rows.append(
                {
                    "epsilon_0": epsilon,
                    "case_id": selected_case["case_id"],
                    "sorted_index": sorted_index,
                    "role": role,
                    "n_points_coarse": INTEGRATION_GRIDS[0],
                    "n_points_fine": INTEGRATION_GRIDS[1],
                    "chi_axial_801": coarse_fractions[0],
                    "chi_bending_801": coarse_fractions[1],
                    "chi_shear_801": coarse_fractions[2],
                    "chi_axial_1601": fine_fractions[0],
                    "chi_bending_1601": fine_fractions[1],
                    "chi_shear_1601": fine_fractions[2],
                    "delta_chi_a": convergence_deltas[0],
                    "delta_chi_b": convergence_deltas[1],
                    "delta_chi_s": convergence_deltas[2],
                    "integration_tolerance": INTEGRATION_TOL,
                    "passed": _bool_text(convergence_pass),
                }
            )
            warning_items = [*mode.warnings, *matrix_warnings, *fields["warnings"]]
            if null_vector_source.endswith("fallback"):
                warning_items.append(
                    "timo_mode_coefficients_roundoff_row_scaling_replaced_by_certificate_scaling"
                )
            warnings = tuple(dict.fromkeys(warning_items))
            residual_rows.append(
                {
                    "epsilon_0": epsilon,
                    "case_id": selected_case["case_id"],
                    "sorted_index": sorted_index,
                    "role": role,
                    "Lambda_Timo": lambda_timo,
                    "null_vector_source": null_vector_source,
                    "timo_mode_coefficients_smallest_singular_value": mode.smallest_singular_value,
                    "timo_mode_coefficients_relative_residual_on_certificate_scaling": mode_relative_residual,
                    "smallest_singular_value": smallest,
                    "second_smallest_singular_value": second_smallest,
                    "singular_value_ratio": singular_ratio,
                    "relative_null_residual": relative_null_residual,
                    "rod1_clamp_u": clamp_values[0],
                    "rod1_clamp_w": clamp_values[1],
                    "rod1_clamp_psi": clamp_values[2],
                    "rod2_clamp_u": clamp_values[3],
                    "rod2_clamp_w": clamp_values[4],
                    "rod2_clamp_psi": clamp_values[5],
                    "max_abs_clamp_gap": max_clamp,
                    "gap_w": joint["gap_w"],
                    "gap_u": joint["gap_u"],
                    "gap_psi": joint["gap_psi"],
                    "gap_M": joint["gap_M"],
                    "gap_Q": joint["gap_Q"],
                    "gap_N": joint["gap_N"],
                    "max_abs_kinematic_gap": joint["max_abs_kinematic_gap"],
                    "max_abs_force_gap": joint["max_abs_force_gap"],
                    "sign_flip_delta_chi_a": sign_deltas[0],
                    "sign_flip_delta_chi_b": sign_deltas[1],
                    "sign_flip_delta_chi_s": sign_deltas[2],
                    "sign_flip_max_abs_delta": max(sign_deltas),
                    "mode_residual_tolerance": MODE_RESIDUAL_TOL,
                    "energy_tolerance": ENERGY_SUM_TOL,
                    "svd_pass": _bool_text(svd_pass),
                    "null_residual_pass": _bool_text(null_pass),
                    "clamp_pass": _bool_text(clamp_pass),
                    "joint_kinematic_pass": _bool_text(joint_kinematic_pass),
                    "joint_force_pass": _bool_text(joint_force_pass),
                    "sign_flip_pass": _bool_text(sign_flip_pass),
                    "passed": _bool_text(residual_pass),
                    "warnings": "; ".join(warnings),
                }
            )

    return energy_rows, convergence_rows, residual_rows


def _percent_text(delta_f: float) -> str:
    percent = 100.0 * float(delta_f)
    if abs(percent) < 1.0e-3:
        return f"{percent:.2e}"
    return f"{percent:.3f}"


def write_figure(output_dir: Path, energy_rows: Sequence[Mapping[str, object]]) -> Path:
    source_fields = (
        "epsilon_0",
        "case_id",
        "beta_deg",
        "mu",
        "eta",
        "sorted_index",
        "role",
        "delta_f",
        "delta_f_percent",
        "chi_axial_Timo",
        "chi_bending_shear_Timo",
    )
    source_rows = [
        {
            field: (100.0 * float(row["delta_f"]) if field == "delta_f_percent" else row[field])
            for field in source_fields
        }
        for row in energy_rows
    ]
    _atomic_csv(
        output_dir / "energy_triplet_comparison_source.csv",
        source_fields,
        source_rows,
    )

    x_values = np.array([0.0, 1.0, 2.0, 4.0, 5.0, 6.0], dtype=float)
    chi_a = np.array([float(row["chi_axial_Timo"]) for row in energy_rows], dtype=float)
    chi_bs = np.array(
        [float(row["chi_bending_shear_Timo"]) for row in energy_rows], dtype=float
    )
    central = [str(row["role"]) == "central_dip" for row in energy_rows]
    line_widths = [1.9 if is_central else 0.7 for is_central in central]
    edge_colors = ["#222222" if is_central else "#5f5f5f" for is_central in central]

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.bar(
        x_values,
        chi_a,
        width=0.72,
        color="#4C78A8",
        edgecolor=edge_colors,
        linewidth=line_widths,
        label="продольная деформация",
        zorder=3,
    )
    ax.bar(
        x_values,
        chi_bs,
        width=0.72,
        bottom=chi_a,
        color="#F2B66D",
        edgecolor=edge_colors,
        linewidth=line_widths,
        label="изгиб и поперечный сдвиг",
        zorder=3,
    )
    labels = [
        f"$\\varepsilon_0=0{{,}}030$\n$k={int(row['sorted_index'])}$"
        if index < 3
        else f"$\\varepsilon_0=0{{,}}050$\n$k={int(row['sorted_index'])}$"
        for index, row in enumerate(energy_rows)
    ]
    ax.set_xticks(x_values, labels)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Относительный вклад в потенциальную энергию")
    ax.grid(axis="y", color="0.88", linewidth=0.6, zorder=0)
    ax.grid(axis="x", visible=False)
    for x_value, row in zip(x_values, energy_rows):
        ax.text(
            x_value,
            0.985,
            f"$\\delta_f={_percent_text(float(row['delta_f']))}\\%$",
            ha="center",
            va="top",
            fontsize=8.2,
        )
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 0.015), frameon=False, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    output = output_dir / "energy_triplet_comparison.png"
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return output


def _compact_table_rows(energy_rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    return [
        {
            "epsilon_0": f"{float(row['epsilon_0']):.3f}",
            "beta": f"{float(row['beta_deg']):g}",
            "mu": f"{float(row['mu']):g}",
            "eta": f"{float(row['eta']):g}",
            "k": int(row["sorted_index"]),
            "delta_f, %": _percent_text(float(row["delta_f"])),
            "chi_a": f"{float(row['chi_axial_Timo']):.6f}",
            "chi_b+s": f"{float(row['chi_bending_shear_Timo']):.6f}",
        }
        for row in energy_rows
    ]


def write_tables(output_dir: Path, energy_rows: Sequence[Mapping[str, object]]) -> None:
    rows = _compact_table_rows(energy_rows)
    fields = tuple(rows[0])
    _atomic_csv(output_dir / "energy_triplets_table.tsv", fields, rows, delimiter="\t")

    markdown = [
        "| " + " | ".join(fields) + " |",
        "| " + " | ".join("---" for _ in fields) + " |",
    ]
    markdown.extend("| " + " | ".join(str(row[field]) for field in fields) + " |" for row in rows)
    _atomic_text(output_dir / "energy_triplets_table.md", "\n".join(markdown) + "\n")

    tex_headers = (
        "$\\epsilon_0$",
        "$\\beta$, град",
        "$\\mu$",
        "$\\eta$",
        "$k$",
        "$\\delta_f$, \\%",
        "$\\chi_a$",
        "$\\chi_{b+s}$",
    )
    tex = [
        "\\begin{tabular}{rrrrrrrr}",
        "\\hline",
        " & ".join(tex_headers) + " \\\\",
        "\\hline",
    ]
    tex.extend(" & ".join(str(row[field]) for field in fields) + " \\\\" for row in rows)
    tex.extend(("\\hline", "\\end{tabular}"))
    _atomic_text(output_dir / "energy_triplets_table.tex", "\n".join(tex) + "\n")


def evaluate_gates(
    selected: Sequence[Mapping[str, Any]],
    energy_rows: Sequence[Mapping[str, object]],
    convergence_rows: Sequence[Mapping[str, object]],
    residual_rows: Sequence[Mapping[str, object]],
    counters: Mapping[str, int],
) -> tuple[dict[str, bool], str, dict[float, bool]]:
    scope_ok = (
        len(selected) == 2
        and all(row["certificate"]["scientific_scope"] == SCIENTIFIC_SCOPE for row in selected)
        and not any("anisotropic_rods" in name for name in sys.modules)
    )
    selection_ok = len(selected) == 2 and all(
        float(row["delta_center"]) < float(row["delta_minus"])
        and float(row["delta_center"]) < float(row["delta_plus"])
        for row in selected
    )
    energy_ok = len(energy_rows) == 6 and all(
        str(row["article_ready_status"]).startswith("article_ready")
        and all(
            math.isfinite(float(row[field]))
            for field in (
                "U_axial_Timo",
                "U_bending_Timo",
                "U_shear_Timo",
                "U_total_Timo",
                "chi_axial_Timo",
                "chi_bending_Timo",
                "chi_shear_Timo",
                "chi_bending_shear_Timo",
            )
        )
        for row in energy_rows
    )
    convergence_ok = len(convergence_rows) == 6 and all(
        row["passed"] == "true" for row in convergence_rows
    )
    residual_ok = len(residual_rows) == 6 and all(
        row["passed"] == "true" for row in residual_rows
    )
    zero_solve_ok = all(
        counters[name] == 0
        for name in (
            "root_solver_calls",
            "strict_solver_calls",
            "family_detector_calls",
            "local_repair_calls",
        )
    )

    support_by_epsilon: dict[float, bool] = {}
    for epsilon, central_index in PILOT_SPECS:
        level = [row for row in energy_rows if _same_float(row["epsilon_0"], epsilon)]
        if len(level) != 3:
            support_by_epsilon[epsilon] = False
            continue
        by_index = {int(row["sorted_index"]): row for row in level}
        support_by_epsilon[epsilon] = (
            float(by_index[central_index]["chi_axial_Timo"])
            > float(by_index[central_index - 1]["chi_axial_Timo"])
            and float(by_index[central_index]["chi_axial_Timo"])
            > float(by_index[central_index + 1]["chi_axial_Timo"])
        )
    interpretation_ok = len(support_by_epsilon) == 2 and all(support_by_epsilon.values())
    gates = {
        "Gate A: scope_isolation_gate": scope_ok,
        "Gate B: candidate_selection_gate": selection_ok,
        "Gate C: energy_reconstruction_gate": energy_ok,
        "Gate D: integration_convergence_gate": convergence_ok,
        "Gate E: mode_residual_gate": residual_ok,
        "Gate F: zero_root_solve_gate": zero_solve_ok,
        "Gate G: interpretation_gate": interpretation_ok,
    }
    gates["Gate H: pilot_readiness_gate"] = all(gates.values())
    supported_count = sum(support_by_epsilon.values())
    scientific_status = (
        "hypothesis_supported_by_two_examples"
        if supported_count == 2
        else "partially_supported"
        if supported_count == 1
        else "hypothesis_not_supported_by_selected_examples"
    )
    return gates, scientific_status, support_by_epsilon


def _markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> list[str]:
    return [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
        *("| " + " | ".join(row) + " |" for row in rows),
    ]


def write_report(
    output_dir: Path,
    records: Sequence[Mapping[str, Any]],
    selected: Sequence[Mapping[str, Any]],
    energy_rows: Sequence[Mapping[str, object]],
    convergence_rows: Sequence[Mapping[str, object]],
    residual_rows: Sequence[Mapping[str, object]],
    counters: Mapping[str, int],
    gates: Mapping[str, bool],
    scientific_status: str,
    support_by_epsilon: Mapping[float, bool],
) -> Path:
    selection_rows = [
        (
            f"{float(row['epsilon_0']):.3f}",
            str(row["case_id"]),
            f"{float(row['beta_deg']):g}",
            f"{float(row['mu']):g}",
            f"{float(row['eta']):g}",
            f"{float(row['delta_minus']):.9g}",
            f"{float(row['delta_center']):.9g}",
            f"{float(row['delta_plus']):.9g}",
            f"{float(row['R_k']):.9g}",
            f"{float(row['D_k']):.9g}",
        )
        for row in selected
    ]
    energy_table_rows = [
        (
            f"{float(row['epsilon_0']):.3f}",
            str(row["sorted_index"]),
            str(row["role"]),
            f"{float(row['delta_f']):.9g}",
            f"{float(row['U_axial_Timo']):.9g}",
            f"{float(row['U_bending_Timo']):.9g}",
            f"{float(row['U_shear_Timo']):.9g}",
            f"{float(row['chi_axial_Timo']):.9g}",
            f"{float(row['chi_bending_shear_Timo']):.9g}",
            str(row["article_ready_status"]),
        )
        for row in energy_rows
    ]
    lines = [
        "# Энергетический пилот для локальных минимумов расхождения частот",
        "",
        f"Scientific scope: `{SCIENTIFIC_SCOPE}`.",
        "",
        "Рассматриваются только одинаковые номера в упорядоченных спектрах двух сопряжённых изотропных стержней круглого сечения. Физические ветви не отслеживаются; MAC, формы колебаний, FEM и поиск новых корней не используются.",
        "",
        "## Выбор геометрий",
        "",
        "Геометрии выбраны до восстановления форм только по минимальному `R_k`; при равенстве использовались максимальный `D_k`, затем `beta`, `mu`, `eta` и `case_id`. Энергетические величины в выборе не участвуют. Полный ранжированный список и первые десять кандидатов каждого уровня сохранены отдельно.",
        "",
        *_markdown_table(
            ("epsilon_0", "case_id", "beta", "mu", "eta", "delta_-", "delta_k", "delta_+", "R_k", "D_k"),
            selection_rows,
        ),
        "",
        f"Допустимых кандидатов: epsilon_0=0.030 — {sum(_same_float(row['epsilon_0'], 0.030) for row in records)}; epsilon_0=0.050 — {sum(_same_float(row['epsilon_0'], 0.050) for row in records)}.",
        f"Обе выбранные геометрии имеют `thinness_flag=true` (`s_max`={float(selected[0]['s_max']):.6g} и {float(selected[1]['s_max']):.6g}); поэтому вывод ограничен заданным аналитическим пилотом и не заменяет отдельную пространственную валидацию этих геометрий.",
        "",
        "## Потенциальная энергия форм Тимошенко",
        "",
        "Основные величины — относительный вклад продольной деформации в потенциальную энергию `chi_a` и совокупный относительный вклад изгиба и поперечного сдвига `chi_b+s`. Классификация по пороговому значению не применяется. Вспомогательная величина `chi_axial_EB` не вычислялась.",
        "",
        *_markdown_table(
            ("epsilon_0", "k", "роль", "delta_f", "U_a", "U_b", "U_s", "chi_a", "chi_b+s", "статус"),
            energy_table_rows,
        ),
        "",
        "## Численные проверки",
        "",
        f"Допуск суммы энергетических вкладов и сходимости интегрирования: `{ENERGY_SUM_TOL:.1e}`. Допуск невязок формы: `{MODE_RESIDUAL_TOL:.1e}`. Для каждой формы сопоставлены сетки 801 и 1601 точек и отдельно проверена инвариантность при смене знака null vector.",
        "",
        f"Максимальные |delta chi| между сетками: chi_a={max(float(row['delta_chi_a']) for row in convergence_rows):.6e}, chi_b={max(float(row['delta_chi_b']) for row in convergence_rows):.6e}, chi_s={max(float(row['delta_chi_s']) for row in convergence_rows):.6e}.",
        f"Максимумы для шести форм: sigma_min={max(float(row['smallest_singular_value']) for row in residual_rows):.6e}, relative null residual={max(float(row['relative_null_residual']) for row in residual_rows):.6e}, clamp={max(float(row['max_abs_clamp_gap']) for row in residual_rows):.6e}, joint kinematic={max(float(row['max_abs_kinematic_gap']) for row in residual_rows):.6e}, joint force={max(float(row['max_abs_force_gap']) for row in residual_rows):.6e}, sign-flip delta={max(float(row['sign_flip_max_abs_delta']) for row in residual_rows):.6e}.",
        f"Для {sum(str(row['null_vector_source']).endswith('fallback') for row in residual_rows)} из шести форм потребовалась SVD-нормировка с floored row scaling, совпадающая с provenance compact solver; сохранённая частота при этом не изменялась.",
        "",
        "## Gates",
        "",
        *[f"- {name}: {'PASS' if passed else 'FAIL'}" for name, passed in gates.items()],
        "",
        f"Научный статус: `{scientific_status}`.",
        "",
        "## Интерпретация",
        "",
    ]
    if all(support_by_epsilon.values()) and len(support_by_epsilon) == 2:
        lines.append(
            "В выбранных примерах частота, для которой расхождение моделей Эйлера—Бернулли и Тимошенко существенно меньше, чем для соседних частот, соответствует форме с более высоким относительным вкладом продольной деформации в потенциальную энергию."
        )
    else:
        failed = ", ".join(
            f"epsilon_0={epsilon:.3f}"
            for epsilon, supported in support_by_epsilon.items()
            if not supported
        )
        lines.append(
            f"Для {failed} центральная форма не имеет более высокого относительного вклада продольной деформации, чем обе соседние формы; в этом примере гипотеза не подтверждена."
        )
    lines.extend(
        (
            "",
            "Приведённые случаи служат иллюстрацией механизма и не устанавливают универсальной зависимости для всех геометрий и всех номеров спектра.",
            "",
            "## Operation counters",
            "",
            *[f"- `{name}`: {value}" for name, value in counters.items()],
            "",
            "Исходные compact certificates использованы только для чтения. Scientific formulas, сохранённые корни и частоты не изменялись.",
        )
    )
    return _atomic_text(output_dir / "report.md", "\n".join(lines) + "\n")


def _write_source_integrity(
    output_dir: Path,
    input_hashes: Mapping[str, str],
) -> bool:
    rows = []
    passed = True
    for path_text, before_hash in sorted(input_hashes.items()):
        path = Path(path_text)
        after_hash = _sha256(path)
        unchanged = before_hash == after_hash
        passed = passed and unchanged
        rows.append(
            {
                "path": path.resolve().as_posix(),
                "sha256_before": before_hash,
                "sha256_after": after_hash,
                "unchanged": _bool_text(unchanged),
            }
        )
    _atomic_csv(
        output_dir / "source_integrity.csv",
        ("path", "sha256_before", "sha256_after", "unchanged"),
        rows,
    )
    return passed


def run_pilot(
    finalization_dir: Path = DEFAULT_FINALIZATION_DIR,
    compact_dir: Path = DEFAULT_COMPACT_DIR,
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    *,
    selection_only: bool = False,
) -> dict[str, object]:
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    records, selected, input_hashes = select_candidates(finalization_dir, compact_dir)

    # This write is deliberately completed before reconstruction starts.
    write_selection_outputs(output_dir, records, selected)
    counters = operation_counters()
    if selection_only:
        _atomic_csv(output_dir / "operation_counters.csv", OPERATION_COUNTER_FIELDS, [counters])
        source_integrity = _write_source_integrity(output_dir, input_hashes)
        return {
            "selection_only": True,
            "selected_case_ids": [row["case_id"] for row in selected],
            "candidate_counts": {
                f"{epsilon:.3f}": sum(_same_float(row["epsilon_0"], epsilon) for row in records)
                for epsilon, _central_index in PILOT_SPECS
            },
            "source_integrity": source_integrity,
            "operation_counters": counters,
        }

    energy_rows, convergence_rows, residual_rows = reconstruct_selected_modes(selected, counters)
    _atomic_csv(output_dir / "energy_triplets.csv", ENERGY_FIELDS, energy_rows)
    _atomic_csv(
        output_dir / "integration_convergence.csv",
        CONVERGENCE_FIELDS,
        convergence_rows,
    )
    _atomic_csv(output_dir / "mode_residual_audit.csv", RESIDUAL_FIELDS, residual_rows)
    _atomic_csv(output_dir / "operation_counters.csv", OPERATION_COUNTER_FIELDS, [counters])
    write_figure(output_dir, energy_rows)
    write_tables(output_dir, energy_rows)

    gates, scientific_status, support_by_epsilon = evaluate_gates(
        selected,
        energy_rows,
        convergence_rows,
        residual_rows,
        counters,
    )
    _atomic_csv(
        output_dir / "gate_summary.csv",
        ("gate", "status"),
        (
            {"gate": name, "status": "PASS" if passed else "FAIL"}
            for name, passed in gates.items()
        ),
    )
    write_report(
        output_dir,
        records,
        selected,
        energy_rows,
        convergence_rows,
        residual_rows,
        counters,
        gates,
        scientific_status,
        support_by_epsilon,
    )
    source_integrity = _write_source_integrity(output_dir, input_hashes)
    if not source_integrity:
        raise RuntimeError("compact input integrity check failed")
    return {
        "selection_only": False,
        "selected_case_ids": [row["case_id"] for row in selected],
        "energy_row_count": len(energy_rows),
        "scientific_status": scientific_status,
        "gates": gates,
        "source_integrity": source_integrity,
        "operation_counters": counters,
        "output_dir": output_dir.as_posix(),
    }


def read_beta0_control(beta0_control_dir: Path) -> list[dict[str, object]]:
    """Read the immutable epsilon=0.030 beta=0 control without reconstruction."""

    path = beta0_control_dir.resolve() / "energy_triplets.csv"
    if not path.is_file():
        raise FileNotFoundError(f"beta=0 control energy table is missing: {path}")
    rows = [
        row
        for row in _read_csv(path)
        if _same_float(row["epsilon_0"], 0.030)
        and _same_float(row["beta_deg"], 0.0)
        and int(row["sorted_index"]) in (4, 5, 6)
    ]
    rows.sort(key=lambda row: int(row["sorted_index"]))
    if [int(row["sorted_index"]) for row in rows] != [4, 5, 6]:
        raise ValueError("immutable beta=0 control must contain epsilon_0=0.030, k=4,5,6")
    return [{field: row.get(field, "") for field in ENERGY_FIELDS} for row in rows]


def _summary_row(
    *,
    example_role: str,
    rows: Sequence[Mapping[str, object]],
    selected_case: Mapping[str, Any] | None,
    notes: str,
) -> dict[str, object]:
    ordered = sorted(rows, key=lambda row: int(row["sorted_index"]))
    if len(ordered) != 3:
        raise ValueError(f"{example_role}: expected three consecutive energy rows")
    left, center, right = ordered
    delta_minus = float(left["delta_f"])
    delta_center = float(center["delta_f"])
    delta_plus = float(right["delta_f"])
    neighbor_mean = 0.5 * (delta_minus + delta_plus)
    chi_left = float(left["chi_axial_Timo"])
    chi_center = float(center["chi_axial_Timo"])
    chi_right = float(right["chi_axial_Timo"])
    if selected_case is None:
        mu = float(center["mu"])
        eta = float(center["eta"])
        s_max = float(center["s_max"])
        thinness_flag = str(center["thinness_flag"])
    else:
        mu = float(selected_case["mu"])
        eta = float(selected_case["eta"])
        s_max = float(selected_case["s_max"])
        thinness_flag = _bool_text(selected_case["thinness_flag"])
    return {
        "example_role": example_role,
        "epsilon_0": float(center["epsilon_0"]),
        "case_id": center["case_id"],
        "beta_deg": float(center["beta_deg"]),
        "mu": mu,
        "eta": eta,
        "central_sorted_index": int(center["sorted_index"]),
        "delta_minus": delta_minus,
        "delta_center": delta_center,
        "delta_plus": delta_plus,
        "R_k": delta_center / neighbor_mean,
        "D_k": neighbor_mean - delta_center,
        "chi_axial_left": chi_left,
        "chi_axial_center": chi_center,
        "chi_axial_right": chi_right,
        "chi_bending_shear_left": float(left["chi_bending_shear_Timo"]),
        "chi_bending_shear_center": float(center["chi_bending_shear_Timo"]),
        "chi_bending_shear_right": float(right["chi_bending_shear_Timo"]),
        "Q_left": chi_center / max(chi_left, CONTRAST_NUMERICAL_FLOOR),
        "Q_right": chi_center / max(chi_right, CONTRAST_NUMERICAL_FLOOR),
        "contrast_numerical_floor": CONTRAST_NUMERICAL_FLOOR,
        "s_max": s_max,
        "thinness_flag": thinness_flag,
        "center_exceeds_both_neighbors": _bool_text(
            chi_center > chi_left and chi_center > chi_right
        ),
        "notes": notes,
    }


def build_coupled_vs_beta0_summary(
    beta0_rows: Sequence[Mapping[str, object]],
    coupled_rows: Sequence[Mapping[str, object]],
    selected: Sequence[Mapping[str, Any]],
) -> list[dict[str, object]]:
    summary = [
        _summary_row(
            example_role="beta0_control_reference",
            rows=beta0_rows,
            selected_case=None,
            notes=(
                "immutable zero-solve control; exact separation of axial and bending motions "
                "in the straight beta=0 limit"
            ),
        )
    ]
    selected_by_epsilon = {float(row["epsilon_0"]): row for row in selected}
    for epsilon, _central_index in PILOT_SPECS:
        rows = [row for row in coupled_rows if _same_float(row["epsilon_0"], epsilon)]
        selected_case = selected_by_epsilon[epsilon]
        notes = (
            "beta=15 coupled example within s_max<=0.1"
            if epsilon == 0.030
            else "beta=15 coupled example in the extended one-dimensional-model range"
        )
        summary.append(
            _summary_row(
                example_role=(
                    "beta15_coupled_smax_le_0p1"
                    if epsilon == 0.030
                    else "beta15_coupled_extended_range"
                ),
                rows=rows,
                selected_case=selected_case,
                notes=notes,
            )
        )
    return summary


def write_coupled_figure(
    output_dir: Path,
    energy_rows: Sequence[Mapping[str, object]],
) -> Path:
    source_fields = (
        "epsilon_0",
        "case_id",
        "beta_deg",
        "mu",
        "eta",
        "sorted_index",
        "role",
        "delta_f",
        "delta_f_percent",
        "chi_axial_Timo",
        "chi_bending_shear_Timo",
        "s_max",
    )
    source_rows = [
        {
            field: (100.0 * float(row["delta_f"]) if field == "delta_f_percent" else row[field])
            for field in source_fields
        }
        for row in energy_rows
    ]
    _atomic_csv(
        output_dir / "coupled_energy_triplet_comparison_source.csv",
        source_fields,
        source_rows,
    )

    x_values = np.array([0.0, 1.0, 2.0, 4.0, 5.0, 6.0], dtype=float)
    chi_a = np.array([float(row["chi_axial_Timo"]) for row in energy_rows], dtype=float)
    chi_bs = np.array(
        [float(row["chi_bending_shear_Timo"]) for row in energy_rows], dtype=float
    )
    central = [str(row["role"]) == "central_dip" for row in energy_rows]
    line_widths = [1.9 if is_central else 0.7 for is_central in central]
    edge_colors = ["#222222" if is_central else "#5f5f5f" for is_central in central]

    fig, ax = plt.subplots(figsize=(8.6, 5.1))
    ax.bar(
        x_values,
        chi_a,
        width=0.72,
        color="#4C78A8",
        edgecolor=edge_colors,
        linewidth=line_widths,
        label="Продольная деформация",
        zorder=3,
    )
    ax.bar(
        x_values,
        chi_bs,
        width=0.72,
        bottom=chi_a,
        color="#F2B66D",
        edgecolor=edge_colors,
        linewidth=line_widths,
        label="Изгиб и поперечный сдвиг",
        zorder=3,
    )
    labels = [f"$k={int(row['sorted_index'])}$" for row in energy_rows]
    ax.set_xticks(x_values, labels)
    ax.text(
        1.0,
        -0.13,
        "$\\varepsilon_0=0.030,\\ \\beta=15^\\circ$",
        ha="center",
        va="top",
        transform=ax.get_xaxis_transform(),
    )
    ax.text(
        5.0,
        -0.13,
        "$\\varepsilon_0=0.050,\\ \\beta=15^\\circ$",
        ha="center",
        va="top",
        transform=ax.get_xaxis_transform(),
    )
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Относительный вклад в потенциальную энергию")
    ax.grid(axis="y", color="0.88", linewidth=0.6, zorder=0)
    ax.grid(axis="x", visible=False)
    for x_value, row in zip(x_values, energy_rows):
        ax.text(
            x_value,
            0.985,
            f"$\\delta_f={_percent_text(float(row['delta_f']))}\\%$",
            ha="center",
            va="top",
            fontsize=8.0,
        )
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.005), frameon=False, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.subplots_adjust(left=0.14, right=0.99, bottom=0.19, top=0.88)
    output = output_dir / "coupled_energy_triplet_comparison.png"
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return output


def _coupled_table_rows(
    energy_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    return [
        {
            "epsilon_0": f"{float(row['epsilon_0']):.3f}",
            "beta": f"{float(row['beta_deg']):g}",
            "mu": f"{float(row['mu']):g}",
            "eta": f"{float(row['eta']):g}",
            "k": int(row["sorted_index"]),
            "delta_f, %": _percent_text(float(row["delta_f"])),
            "chi_a": f"{float(row['chi_axial_Timo']):.6f}",
            "chi_b+s": f"{float(row['chi_bending_shear_Timo']):.6f}",
            "s_max": f"{float(row['s_max']):.6f}",
        }
        for row in energy_rows
    ]


def write_coupled_tables(
    output_dir: Path,
    energy_rows: Sequence[Mapping[str, object]],
) -> None:
    rows = _coupled_table_rows(energy_rows)
    fields = tuple(rows[0])
    _atomic_csv(
        output_dir / "coupled_energy_triplets_table.tsv",
        fields,
        rows,
        delimiter="\t",
    )
    note = (
        "Примечание: пример epsilon_0=0.030 удовлетворяет s_max<=0.1; "
        "пример epsilon_0=0.050 относится к расширенному диапазону одномерной модели."
    )
    markdown = [
        "| " + " | ".join(fields) + " |",
        "| " + " | ".join("---" for _ in fields) + " |",
        *("| " + " | ".join(str(row[field]) for field in fields) + " |" for row in rows),
        "",
        note,
    ]
    _atomic_text(
        output_dir / "coupled_energy_triplets_table.md",
        "\n".join(markdown) + "\n",
    )
    tex_headers = (
        "$\\epsilon_0$",
        "$\\beta$, град",
        "$\\mu$",
        "$\\eta$",
        "$k$",
        "$\\delta_f$, \\%",
        "$\\chi_a$",
        "$\\chi_{b+s}$",
        "$s_{\\max}$",
    )
    tex = [
        "\\begin{tabular}{rrrrrrrrr}",
        "\\hline",
        " & ".join(tex_headers) + " \\\\",
        "\\hline",
        *(" & ".join(str(row[field]) for field in fields) + " \\\\" for row in rows),
        "\\hline",
        "\\end{tabular}",
        "",
        "\\par\\smallskip\\noindent\\footnotesize "
        "Примечание: пример $\\epsilon_0=0.030$ удовлетворяет $s_{\\max}\\leq0.1$; "
        "пример $\\epsilon_0=0.050$ относится к расширенному диапазону одномерной модели.",
    ]
    _atomic_text(
        output_dir / "coupled_energy_triplets_table.tex",
        "\n".join(tex) + "\n",
    )


def evaluate_coupled_gates(
    selected: Sequence[Mapping[str, Any]],
    energy_rows: Sequence[Mapping[str, object]],
    convergence_rows: Sequence[Mapping[str, object]],
    residual_rows: Sequence[Mapping[str, object]],
    summary_rows: Sequence[Mapping[str, object]],
    beta0_rows: Sequence[Mapping[str, object]],
    counters: Mapping[str, int],
    selection_context: Mapping[str, object],
    *,
    compact_inputs_unchanged: bool,
    beta0_control_unchanged: bool,
) -> tuple[dict[str, bool], str]:
    scope_ok = (
        len(selected) == 2
        and all(row["certificate"]["scientific_scope"] == SCIENTIFIC_SCOPE for row in selected)
        and not any("anisotropic_rods" in name for name in sys.modules)
    )
    selection_ok = (
        len(selected) == 2
        and compact_inputs_unchanged
        and not bool(selection_context["selection_uses_energy_values"])
        and int(selection_context["main_grid_smax_le_0p1_count_epsilon_0p050"]) == 0
        and all(_same_float(row["beta_deg"], COUPLED_SELECTED_BETA_DEG) for row in selected)
        and all(
            float(row["delta_center"]) < float(row["delta_minus"])
            and float(row["delta_center"]) < float(row["delta_plus"])
            for row in selected
        )
        and selected[0]["selection_pool"] == "nonzero_angle_smax_le_0p1"
        and float(selected[0]["s_max"]) <= THINNESS_LIMIT
        and selected[1]["selection_pool"] == "nonzero_angle_all"
    )
    zero_solve_ok = all(
        counters[name] == 0
        for name in (
            "root_solver_calls",
            "strict_solver_calls",
            "family_detector_calls",
            "local_repair_calls",
            "shape_plot_calls",
            "MAC_calls",
            "FEM_calls",
            "anisotropic_calls",
        )
    )
    energy_ok = len(energy_rows) == 6 and all(
        str(row["article_ready_status"]) == "article_ready"
        and all(
            math.isfinite(float(row[field]))
            for field in (
                "U_axial_Timo",
                "U_bending_Timo",
                "U_shear_Timo",
                "U_total_Timo",
                "chi_axial_Timo",
                "chi_bending_Timo",
                "chi_shear_Timo",
                "chi_bending_shear_Timo",
            )
        )
        for row in energy_rows
    )
    convergence_ok = len(convergence_rows) == 6 and all(
        row["passed"] == "true" for row in convergence_rows
    )
    residual_ok = len(residual_rows) == 6 and all(
        row["passed"] == "true" for row in residual_rows
    )
    coupled_summary = [row for row in summary_rows if str(row["example_role"]).startswith("beta15")]
    interpretation_ok = len(coupled_summary) == 2 and all(
        row["center_exceeds_both_neighbors"] == "true" for row in coupled_summary
    )
    control_ok = (
        beta0_control_unchanged
        and len(beta0_rows) == 3
        and [int(row["sorted_index"]) for row in beta0_rows] == [4, 5, 6]
        and all(_same_float(row["beta_deg"], 0.0) for row in beta0_rows)
    )
    gates = {
        "Gate A: scope_isolation_gate": scope_ok,
        "Gate B: coupled_candidate_selection_gate": selection_ok,
        "Gate C: zero_root_solve_gate": zero_solve_ok,
        "Gate D: energy_reconstruction_gate": energy_ok,
        "Gate E: integration_convergence_gate": convergence_ok,
        "Gate F: mode_residual_gate": residual_ok,
        "Gate G: coupled_interpretation_gate": interpretation_ok,
        "Gate H: beta0_control_preservation_gate": control_ok,
    }
    gates["Gate I: coupled_energy_pilot_readiness_gate"] = all(gates.values())
    supported_count = sum(row["center_exceeds_both_neighbors"] == "true" for row in coupled_summary)
    scientific_status = (
        "coupled_hypothesis_supported_by_two_examples"
        if supported_count == 2
        else "coupled_hypothesis_partially_supported"
        if supported_count == 1
        else "coupled_hypothesis_not_supported"
    )
    return gates, scientific_status


def write_coupled_report(
    output_dir: Path,
    selected: Sequence[Mapping[str, Any]],
    energy_rows: Sequence[Mapping[str, object]],
    convergence_rows: Sequence[Mapping[str, object]],
    residual_rows: Sequence[Mapping[str, object]],
    summary_rows: Sequence[Mapping[str, object]],
    counters: Mapping[str, int],
    gates: Mapping[str, bool],
    scientific_status: str,
    selection_context: Mapping[str, object],
) -> Path:
    selection_table = [
        (
            f"{float(row['epsilon_0']):.3f}",
            str(row["case_id"]),
            f"{float(row['beta_deg']):g}",
            f"{float(row['mu']):g}",
            f"{float(row['eta']):g}",
            str(row["selection_pool"]),
            f"{float(row['delta_minus']):.9g}",
            f"{float(row['delta_center']):.9g}",
            f"{float(row['delta_plus']):.9g}",
            f"{float(row['R_k']):.9g}",
            f"{float(row['D_k']):.9g}",
            f"{float(row['s_max']):.9g}",
        )
        for row in selected
    ]
    energy_table = [
        (
            f"{float(row['epsilon_0']):.3f}",
            str(row["sorted_index"]),
            f"{100.0 * float(row['delta_f']):.6g}",
            f"{float(row['U_axial_Timo']):.8g}",
            f"{float(row['U_bending_Timo']):.8g}",
            f"{float(row['U_shear_Timo']):.8g}",
            f"{float(row['chi_axial_Timo']):.8g}",
            f"{float(row['chi_bending_shear_Timo']):.8g}",
        )
        for row in energy_rows
    ]
    comparison_table = [
        (
            str(row["example_role"]),
            f"{float(row['epsilon_0']):.3f}",
            f"{float(row['beta_deg']):g}",
            f"{float(row['chi_axial_left']):.7g}",
            f"{float(row['chi_axial_center']):.7g}",
            f"{float(row['chi_axial_right']):.7g}",
            f"{float(row['Q_left']):.7g}",
            f"{float(row['Q_right']):.7g}",
        )
        for row in summary_rows
    ]
    lines = [
        "# Энергетический пилот для углового сопряжения beta=15 deg",
        "",
        f"Scientific scope: `{SCIENTIFIC_SCOPE}`.",
        "",
        "Рассматриваются одинаковые номера в упорядоченных спектрах двух сопряжённых изотропных стержней круглого сечения. Новые собственные значения не вычислялись, а физические ветви не отслеживались.",
        "",
        "## Выбор геометрий",
        "",
        "Геометрии зафиксированы до восстановления форм. Выбор использует только сохранённые `delta_f`, показатели `R_k`, `D_k` и геометрические ограничения; энергетические величины в ранжировании не участвуют.",
        "",
        *_markdown_table(
            ("epsilon_0", "case_id", "beta", "mu", "eta", "pool", "delta_-", "delta_k", "delta_+", "R_k", "D_k", "s_max"),
            selection_table,
        ),
        "",
        f"Для epsilon_0=0.030 получено {selection_context['nonzero_angle_all_candidate_count_epsilon_0p030']} допустимых beta>=15 кандидатов, из них {selection_context['nonzero_angle_smax_le_0p1_candidate_count_epsilon_0p030']} удовлетворяют s_max<=0.1. Первый пример выбран из этого ограниченного pool.",
        f"Для epsilon_0=0.050 получено {selection_context['nonzero_angle_all_candidate_count_epsilon_0p050']} допустимых beta>=15 кандидатов. Среди всех 194 геометрий основной сетки нет ни одной с s_max<=0.1, поэтому второй пример выбран из полного nonzero-angle pool и относится к расширенному исследованию одномерных моделей.",
        "",
        "## Потенциальная энергия форм Тимошенко",
        "",
        "Основные показатели — относительный вклад продольной деформации в потенциальную энергию `chi_a` и совокупный относительный вклад изгиба и поперечного сдвига `chi_b+s`. Пороговая классификация форм не используется.",
        "",
        *_markdown_table(
            ("epsilon_0", "k", "delta_f, %", "U_a", "U_b", "U_s", "chi_a", "chi_b+s"),
            energy_table,
        ),
        "",
        "## Сопоставление с контрольным случаем beta=0",
        "",
        *_markdown_table(
            ("роль", "epsilon_0", "beta", "chi_a(left)", "chi_a(center)", "chi_a(right)", "Q_left", "Q_right"),
            comparison_table,
        ),
        "",
        "Контрольная тройка beta=0 прочитана из неизменённого предыдущего пилота без повторного восстановления форм. В этом предельном случае продольные и изгибные движения разделяются практически точно. В отличие от контрольного случая beta=0, при beta=15 deg продольное и поперечное движения связаны условиями сопряжения.",
        "",
        "Показатели `Q_left` и `Q_right` характеризуют только численный контраст и не вводятся как новый физический критерий.",
        "",
        "## Численные проверки",
        "",
        f"Максимальные различия между сетками 801 и 1601 точек: delta_chi_a={max(float(row['delta_chi_a']) for row in convergence_rows):.6e}, delta_chi_b={max(float(row['delta_chi_b']) for row in convergence_rows):.6e}, delta_chi_s={max(float(row['delta_chi_s']) for row in convergence_rows):.6e}; допуск равен {INTEGRATION_TOL:.1e}.",
        f"Максимумы для шести форм: sigma_min={max(float(row['smallest_singular_value']) for row in residual_rows):.6e}, relative null residual={max(float(row['relative_null_residual']) for row in residual_rows):.6e}, clamp={max(float(row['max_abs_clamp_gap']) for row in residual_rows):.6e}, joint kinematic={max(float(row['max_abs_kinematic_gap']) for row in residual_rows):.6e}, joint force={max(float(row['max_abs_force_gap']) for row in residual_rows):.6e}, sign-flip delta={max(float(row['sign_flip_max_abs_delta']) for row in residual_rows):.6e}.",
        "",
        "## Gates",
        "",
        *[f"- {name}: {'PASS' if passed else 'FAIL'}" for name, passed in gates.items()],
        "",
        f"Научный статус: `{scientific_status}`.",
        "",
        "## Интерпретация",
        "",
    ]
    if scientific_status == "coupled_hypothesis_supported_by_two_examples":
        lines.append(
            "В рассмотренных угловых системах частота, для которой расхождение моделей Эйлера—Бернулли и Тимошенко существенно меньше, чем для соседних частот, соответствует форме с повышенным относительным вкладом продольной деформации в потенциальную энергию."
        )
    elif scientific_status == "coupled_hypothesis_partially_supported":
        lines.append(
            "Повышенный относительный вклад продольной деформации центральной формы наблюдается только в одном из двух выбранных угловых примеров."
        )
    else:
        lines.append(
            "В выбранных угловых примерах центральные формы не имеют повышенного относительного вклада продольной деформации по сравнению с обеими соседними формами."
        )
    lines.extend(
        (
            "",
            "Приведённые примеры иллюстрируют возможный механизм локального уменьшения расхождения частот и не устанавливают универсальной зависимости для всех геометрий и номеров спектра.",
            "",
            "## Operation counters",
            "",
            *[f"- `{name}`: {value}" for name, value in counters.items()],
            "",
            "Compact certificates и beta=0 control использованы только для чтения. Scientific formulas, basis, tolerances, сохранённые корни и частоты не изменялись.",
        )
    )
    return _atomic_text(output_dir / "report.md", "\n".join(lines) + "\n")


def run_coupled_pilot(
    finalization_dir: Path = DEFAULT_FINALIZATION_DIR,
    compact_dir: Path = DEFAULT_COMPACT_DIR,
    beta0_control_dir: Path = DEFAULT_BETA0_CONTROL_DIR,
    output_dir: Path = DEFAULT_COUPLED_OUTPUT_DIR,
    *,
    selection_only: bool = False,
) -> dict[str, object]:
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    pool_records, selected, input_hashes, selection_context = select_coupled_candidates(
        finalization_dir,
        compact_dir,
        beta0_control_dir,
    )

    # Frequency-only selection is durably written before any energy table is opened.
    write_coupled_selection_outputs(output_dir, pool_records, selected)
    counters = operation_counters()
    if selection_only:
        _atomic_csv(output_dir / "operation_counts.csv", OPERATION_COUNTER_FIELDS, [counters])
        compact_inputs_unchanged = all(
            _sha256(Path(path_text)) == before_hash
            for path_text, before_hash in input_hashes.items()
        )
        _atomic_json(
            output_dir / "run_metadata.json",
            {
                "scientific_scope": SCIENTIFIC_SCOPE,
                "run_mode": "coupled_angle_selection_only",
                "selected_case_ids": [row["case_id"] for row in selected],
                "selection_context": selection_context,
                "compact_inputs_unchanged": compact_inputs_unchanged,
                "beta0_energy_data_opened": False,
                "operation_counts": counters,
            },
        )
        return {
            "selection_only": True,
            "selected_case_ids": [row["case_id"] for row in selected],
            "selection_context": selection_context,
            "compact_inputs_unchanged": compact_inputs_unchanged,
            "operation_counts": counters,
            "output_dir": output_dir.as_posix(),
        }

    beta0_control_dir = beta0_control_dir.resolve()
    beta0_hashes_before = _directory_hashes(beta0_control_dir)
    energy_rows, convergence_rows, residual_rows = reconstruct_selected_modes(selected, counters)
    coupled_energy_rows = [
        {field: row[field] for field in COUPLED_ENERGY_FIELDS}
        for row in energy_rows
    ]
    _atomic_csv(
        output_dir / "coupled_energy_triplets.csv",
        COUPLED_ENERGY_FIELDS,
        coupled_energy_rows,
    )
    _atomic_csv(
        output_dir / "integration_convergence.csv",
        CONVERGENCE_FIELDS,
        convergence_rows,
    )
    _atomic_csv(
        output_dir / "mode_residual_audit.csv",
        RESIDUAL_FIELDS,
        residual_rows,
    )
    _atomic_csv(output_dir / "operation_counts.csv", OPERATION_COUNTER_FIELDS, [counters])

    beta0_rows = read_beta0_control(beta0_control_dir)
    _atomic_csv(
        output_dir / "beta0_control_reference.csv",
        ENERGY_FIELDS,
        beta0_rows,
    )
    summary_rows = build_coupled_vs_beta0_summary(beta0_rows, coupled_energy_rows, selected)
    _atomic_csv(
        output_dir / "coupled_vs_beta0_summary.csv",
        COUPLED_SUMMARY_FIELDS,
        summary_rows,
    )
    write_coupled_figure(output_dir, coupled_energy_rows)
    write_coupled_tables(output_dir, coupled_energy_rows)

    compact_inputs_unchanged = all(
        _sha256(Path(path_text)) == before_hash
        for path_text, before_hash in input_hashes.items()
    )
    beta0_hashes_after = _directory_hashes(beta0_control_dir)
    beta0_control_unchanged = beta0_hashes_before == beta0_hashes_after
    gates, scientific_status = evaluate_coupled_gates(
        selected,
        coupled_energy_rows,
        convergence_rows,
        residual_rows,
        summary_rows,
        beta0_rows,
        counters,
        selection_context,
        compact_inputs_unchanged=compact_inputs_unchanged,
        beta0_control_unchanged=beta0_control_unchanged,
    )
    _atomic_csv(
        output_dir / "gate_summary.csv",
        ("gate", "status"),
        (
            {"gate": name, "status": "PASS" if passed else "FAIL"}
            for name, passed in gates.items()
        ),
    )
    write_coupled_report(
        output_dir,
        selected,
        coupled_energy_rows,
        convergence_rows,
        residual_rows,
        summary_rows,
        counters,
        gates,
        scientific_status,
        selection_context,
    )
    generated_files = sorted(
        path.name
        for path in output_dir.iterdir()
        if path.is_file() and path.name != "run_metadata.json"
    )
    _atomic_json(
        output_dir / "run_metadata.json",
        {
            "scientific_scope": SCIENTIFIC_SCOPE,
            "run_mode": "coupled_angle_pilot",
            "selected_case_ids": [row["case_id"] for row in selected],
            "selection_context": selection_context,
            "energy_row_count": len(coupled_energy_rows),
            "integration_grids": list(INTEGRATION_GRIDS),
            "energy_sum_tolerance": ENERGY_SUM_TOL,
            "integration_tolerance": INTEGRATION_TOL,
            "mode_residual_tolerance": MODE_RESIDUAL_TOL,
            "contrast_numerical_floor": CONTRAST_NUMERICAL_FLOOR,
            "compact_inputs_unchanged": compact_inputs_unchanged,
            "beta0_control_unchanged": beta0_control_unchanged,
            "beta0_control_file_hashes_before": beta0_hashes_before,
            "beta0_control_file_hashes_after": beta0_hashes_after,
            "operation_counts": counters,
            "gates": gates,
            "scientific_status": scientific_status,
            "generated_files": generated_files,
        },
    )
    if not compact_inputs_unchanged or not beta0_control_unchanged:
        raise RuntimeError("immutable compact/control input integrity check failed")
    return {
        "selection_only": False,
        "selected_case_ids": [row["case_id"] for row in selected],
        "energy_row_count": len(coupled_energy_rows),
        "scientific_status": scientific_status,
        "gates": gates,
        "compact_inputs_unchanged": compact_inputs_unchanged,
        "beta0_control_unchanged": beta0_control_unchanged,
        "operation_counts": counters,
        "output_dir": output_dir.as_posix(),
    }


def read_coupled_control_energy(coupled_control_dir: Path) -> list[dict[str, object]]:
    """Read the six immutable beta=15 rows only after additional-case selection."""

    path = coupled_control_dir.resolve() / "coupled_energy_triplets.csv"
    if not path.is_file():
        raise FileNotFoundError(f"beta=15 coupled energy table is missing: {path}")
    rows = _read_csv(path)
    if len(rows) != 6 or not all(_same_float(row["beta_deg"], 15.0) for row in rows):
        raise ValueError("immutable beta=15 control must contain exactly six beta=15 rows")
    expected = {
        (0.030, 4),
        (0.030, 5),
        (0.030, 6),
        (0.050, 3),
        (0.050, 4),
        (0.050, 5),
    }
    actual = {(round(float(row["epsilon_0"]), 3), int(row["sorted_index"])) for row in rows}
    if actual != expected:
        raise ValueError("immutable beta=15 control has an unexpected sorted-index inventory")
    return [{field: row[field] for field in COUPLED_ENERGY_FIELDS} for row in rows]


def build_combined_article_energy_rows(
    beta15_rows: Sequence[Mapping[str, object]],
    beta30_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    beta15_003 = sorted(
        (row for row in beta15_rows if _same_float(row["epsilon_0"], 0.030)),
        key=lambda row: int(row["sorted_index"]),
    )
    beta15_005 = sorted(
        (row for row in beta15_rows if _same_float(row["epsilon_0"], 0.050)),
        key=lambda row: int(row["sorted_index"]),
    )
    beta30_003 = sorted(beta30_rows, key=lambda row: int(row["sorted_index"]))
    if [int(row["sorted_index"]) for row in beta15_003] != [4, 5, 6]:
        raise ValueError("beta=15 epsilon_0=0.030 control must contain k=4,5,6")
    if [int(row["sorted_index"]) for row in beta15_005] != [3, 4, 5]:
        raise ValueError("beta=15 epsilon_0=0.050 control must contain k=3,4,5")
    if [int(row["sorted_index"]) for row in beta30_003] != [4, 5, 6]:
        raise ValueError("additional beta=30 example must contain k=4,5,6")
    reference = beta15_003[0]
    extension = beta30_003[0]
    for field in ("epsilon_0", "mu", "eta"):
        if not _same_float(reference[field], float(extension[field])):
            raise ValueError(f"beta=15/beta=30 geometry mismatch in {field}")
    if not _same_float(reference["beta_deg"], 15.0) or not _same_float(
        extension["beta_deg"], 30.0
    ):
        raise ValueError("angle comparison must change beta from 15 to 30 degrees only")
    ordered = [*beta15_003, *beta30_003, *beta15_005]
    return [{field: row[field] for field in COUPLED_ENERGY_FIELDS} for row in ordered]


def build_beta15_vs_beta30_comparison(
    combined_rows: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], dict[str, object], dict[str, object]]:
    rows15 = [
        row
        for row in combined_rows
        if _same_float(row["epsilon_0"], 0.030) and _same_float(row["beta_deg"], 15.0)
    ]
    rows30 = [
        row
        for row in combined_rows
        if _same_float(row["epsilon_0"], 0.030) and _same_float(row["beta_deg"], 30.0)
    ]
    summary15 = _summary_row(
        example_role="beta15_fixed_geometry_reference",
        rows=rows15,
        selected_case=None,
        notes="immutable beta=15 coupled reference",
    )
    summary30 = _summary_row(
        example_role="beta30_additional_angle_example",
        rows=rows30,
        selected_case=None,
        notes="new beta=30 additional-angle example",
    )
    by15 = {int(row["sorted_index"]): row for row in rows15}
    by30 = {int(row["sorted_index"]): row for row in rows30}
    comparison: list[dict[str, object]] = []
    for sorted_index in (4, 5, 6):
        left = by15[sorted_index]
        right = by30[sorted_index]
        comparison.append(
            {
                "epsilon_0": float(right["epsilon_0"]),
                "mu": float(right["mu"]),
                "eta": float(right["eta"]),
                "case_id_beta15": left["case_id"],
                "case_id_beta30": right["case_id"],
                "beta15_deg": float(left["beta_deg"]),
                "beta30_deg": float(right["beta_deg"]),
                "sorted_index": sorted_index,
                "delta_f_beta15": float(left["delta_f"]),
                "delta_f_beta30": float(right["delta_f"]),
                "delta_f_change_beta30_minus_beta15": float(right["delta_f"])
                - float(left["delta_f"]),
                "chi_axial_beta15": float(left["chi_axial_Timo"]),
                "chi_axial_beta30": float(right["chi_axial_Timo"]),
                "chi_axial_change_beta30_minus_beta15": float(right["chi_axial_Timo"])
                - float(left["chi_axial_Timo"]),
                "chi_bending_shear_beta15": float(left["chi_bending_shear_Timo"]),
                "chi_bending_shear_beta30": float(right["chi_bending_shear_Timo"]),
                "chi_bending_shear_change_beta30_minus_beta15": float(
                    right["chi_bending_shear_Timo"]
                )
                - float(left["chi_bending_shear_Timo"]),
                "Q_left_beta15": float(summary15["Q_left"]),
                "Q_right_beta15": float(summary15["Q_right"]),
                "Q_left_beta30": float(summary30["Q_left"]),
                "Q_right_beta30": float(summary30["Q_right"]),
                "Q_left_change_beta30_minus_beta15": float(summary30["Q_left"])
                - float(summary15["Q_left"]),
                "Q_right_change_beta30_minus_beta15": float(summary30["Q_right"])
                - float(summary15["Q_right"]),
            }
        )
    return comparison, summary15, summary30


ARTICLE_FIGURE_SOURCE_FIELDS = (
    "epsilon_0",
    "case_id",
    "beta_deg",
    "mu",
    "eta",
    "sorted_index",
    "role",
    "delta_f",
    "delta_f_percent",
    "chi_axial_Timo",
    "chi_bending_shear_Timo",
    "s_max",
)


def _write_article_energy_group_figure(
    output_dir: Path,
    energy_rows: Sequence[Mapping[str, object]],
    *,
    figure_name: str,
    source_name: str,
    group_labels: Sequence[str],
    figsize: tuple[float, float],
) -> Path:
    if len(energy_rows) != 3 * len(group_labels):
        raise ValueError("article energy figure requires three bars per group")
    source_rows = [
        {
            field: (
                100.0 * float(row["delta_f"])
                if field == "delta_f_percent"
                else row[field]
            )
            for field in ARTICLE_FIGURE_SOURCE_FIELDS
        }
        for row in energy_rows
    ]
    _atomic_csv(output_dir / source_name, ARTICLE_FIGURE_SOURCE_FIELDS, source_rows)

    x_values = np.array(
        [4.0 * group_index + within for group_index in range(len(group_labels)) for within in range(3)],
        dtype=float,
    )
    chi_a = np.array([float(row["chi_axial_Timo"]) for row in energy_rows], dtype=float)
    chi_bs = np.array(
        [float(row["chi_bending_shear_Timo"]) for row in energy_rows], dtype=float
    )
    central = [str(row["role"]) == "central_dip" for row in energy_rows]
    line_widths = [1.9 if is_central else 0.7 for is_central in central]
    edge_colors = ["#222222" if is_central else "#5f5f5f" for is_central in central]

    fig, ax = plt.subplots(figsize=figsize, facecolor="white")
    ax.set_facecolor("white")
    ax.bar(
        x_values,
        chi_a,
        width=0.72,
        color="#4C78A8",
        edgecolor=edge_colors,
        linewidth=line_widths,
        label="Продольная деформация",
        zorder=3,
    )
    ax.bar(
        x_values,
        chi_bs,
        width=0.72,
        bottom=chi_a,
        color="#F2B66D",
        edgecolor=edge_colors,
        linewidth=line_widths,
        label="Изгиб и поперечный сдвиг",
        zorder=3,
    )
    ax.set_xticks(
        x_values,
        [f"$k={int(row['sorted_index'])}$" for row in energy_rows],
    )
    for group_index, label in enumerate(group_labels):
        ax.text(
            4.0 * group_index + 1.0,
            -0.13,
            label,
            ha="center",
            va="top",
            transform=ax.get_xaxis_transform(),
        )
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Относительный вклад в потенциальную энергию")
    ax.grid(axis="y", color="0.88", linewidth=0.6, zorder=0)
    ax.grid(axis="x", visible=False)
    for x_value, row in zip(x_values, energy_rows):
        ax.text(
            x_value,
            0.985,
            f"$\\delta_f={_percent_text(float(row['delta_f']))}\\%$",
            ha="center",
            va="top",
            fontsize=7.8,
        )
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.005), frameon=False, ncol=2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.19, top=0.88)
    output = output_dir / figure_name
    fig.savefig(output, dpi=600, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return output


def write_additional_angle_figures(
    output_dir: Path,
    combined_rows: Sequence[Mapping[str, object]],
) -> tuple[Path, Path]:
    compact_rows = [
        row
        for row in combined_rows
        if _same_float(row["epsilon_0"], 0.030)
        and (_same_float(row["beta_deg"], 15.0) or _same_float(row["beta_deg"], 30.0))
    ]
    compact = _write_article_energy_group_figure(
        output_dir,
        compact_rows,
        figure_name="article_energy_angle_comparison.png",
        source_name="article_energy_angle_comparison_source.csv",
        group_labels=("$\\beta=15^\\circ$", "$\\beta=30^\\circ$"),
        figsize=(8.6, 5.1),
    )
    extended = _write_article_energy_group_figure(
        output_dir,
        combined_rows,
        figure_name="article_energy_examples_extended.png",
        source_name="article_energy_examples_extended_source.csv",
        group_labels=(
            "$\\varepsilon_0=0.030,\\ \\beta=15^\\circ$",
            "$\\varepsilon_0=0.030,\\ \\beta=30^\\circ$",
            "$\\varepsilon_0=0.050,\\ \\beta=15^\\circ$",
        ),
        figsize=(11.0, 5.1),
    )
    return compact, extended


def write_additional_angle_tables(
    output_dir: Path,
    combined_rows: Sequence[Mapping[str, object]],
) -> None:
    compact_energy = [
        row for row in combined_rows if _same_float(row["epsilon_0"], 0.030)
    ]
    compact_rows = [
        {
            "beta": f"{float(row['beta_deg']):g}",
            "k": int(row["sorted_index"]),
            "delta_f, %": _percent_text(float(row["delta_f"])),
            "chi_a": f"{float(row['chi_axial_Timo']):.6f}",
            "chi_bs": f"{float(row['chi_bending_shear_Timo']):.6f}",
        }
        for row in compact_energy
    ]
    compact_fields = tuple(compact_rows[0])
    _atomic_csv(
        output_dir / "article_energy_angle_comparison_table.tsv",
        compact_fields,
        compact_rows,
        delimiter="\t",
    )
    markdown = [
        "| " + " | ".join(compact_fields) + " |",
        "| " + " | ".join("---" for _ in compact_fields) + " |",
        *("| " + " | ".join(str(row[field]) for field in compact_fields) + " |" for row in compact_rows),
        "",
        "Параметры обеих групп: epsilon_0=0.030, mu=0.5, eta=0.5.",
    ]
    _atomic_text(
        output_dir / "article_energy_angle_comparison_table.md",
        "\n".join(markdown) + "\n",
    )
    tex_headers = ("$\\beta$, град", "$k$", "$\\delta_f$, \\%", "$\\chi_a$", "$\\chi_{b+s}$")
    tex = [
        "\\begin{tabular}{rrrrr}",
        "\\hline",
        " & ".join(tex_headers) + " \\\\",
        "\\hline",
        *(" & ".join(str(row[field]) for field in compact_fields) + " \\\\" for row in compact_rows),
        "\\hline",
        "\\end{tabular}",
        "",
        "\\par\\smallskip\\noindent\\footnotesize "
        "Параметры обеих групп: $\\epsilon_0=0.030$, $\\mu=0.5$, $\\eta=0.5$.",
    ]
    _atomic_text(
        output_dir / "article_energy_angle_comparison_table.tex",
        "\n".join(tex) + "\n",
    )

    extended_rows = _coupled_table_rows(combined_rows)
    extended_fields = tuple(extended_rows[0])
    _atomic_csv(
        output_dir / "article_energy_examples_extended_table.csv",
        extended_fields,
        extended_rows,
    )
    _atomic_csv(
        output_dir / "article_energy_examples_extended_table.tsv",
        extended_fields,
        extended_rows,
        delimiter="\t",
    )


def _article_figure_source_consistent(
    path: Path,
    expected_rows: Sequence[Mapping[str, object]],
) -> bool:
    actual = _read_csv(path)
    if len(actual) != len(expected_rows):
        return False
    for source, expected in zip(actual, expected_rows):
        if str(source["case_id"]) != str(expected["case_id"]):
            return False
        if int(source["sorted_index"]) != int(expected["sorted_index"]):
            return False
        for field in ("epsilon_0", "beta_deg", "mu", "eta", "delta_f", "chi_axial_Timo", "chi_bending_shear_Timo"):
            if not math.isclose(
                float(source[field]),
                float(expected[field]),
                rel_tol=1.0e-12,
                abs_tol=1.0e-15,
            ):
                return False
    return True


def evaluate_additional_angle_gates(
    selected: Mapping[str, Any],
    energy_rows: Sequence[Mapping[str, object]],
    convergence_rows: Sequence[Mapping[str, object]],
    residual_rows: Sequence[Mapping[str, object]],
    summary30: Mapping[str, object],
    counters: Mapping[str, int],
    selection_context: Mapping[str, object],
    *,
    compact_inputs_unchanged: bool,
    prior_results_unchanged: bool,
    root_immutability_pass: bool,
    figure_consistency_pass: bool,
) -> tuple[dict[str, bool], str]:
    scope_ok = (
        selected["certificate"]["scientific_scope"] == SCIENTIFIC_SCOPE
        and not any("anisotropic_rods" in name for name in sys.modules)
    )
    selection_ok = (
        compact_inputs_unchanged
        and not bool(selection_context["selection_uses_energy_values"])
        and int(selection_context["ranking_duplicate_rows_removed"]) >= 1
        and _same_float(selected["epsilon_0"], ADDITIONAL_ANGLE_EPSILON)
        and int(selected["central_sorted_index"]) == ADDITIONAL_ANGLE_CENTRAL_INDEX
        and float(selected["beta_deg"]) >= ADDITIONAL_ANGLE_MIN_BETA_DEG
        and float(selected["s_max"]) <= THINNESS_LIMIT
        and float(selected["delta_center"]) < float(selected["delta_minus"])
        and float(selected["delta_center"]) < float(selected["delta_plus"])
    )
    zero_solve_ok = all(
        counters[name] == 0
        for name in (
            "root_solver_calls",
            "strict_solver_calls",
            "family_detector_calls",
            "local_repair_calls",
            "shape_plot_calls",
            "MAC_calls",
            "FEM_calls",
            "anisotropic_calls",
        )
    )
    energy_ok = root_immutability_pass and len(energy_rows) == 3 and all(
        str(row["article_ready_status"]) == "article_ready"
        and all(
            math.isfinite(float(row[field]))
            for field in (
                "U_axial_Timo",
                "U_bending_Timo",
                "U_shear_Timo",
                "U_total_Timo",
                "chi_axial_Timo",
                "chi_bending_Timo",
                "chi_shear_Timo",
                "chi_bending_shear_Timo",
            )
        )
        for row in energy_rows
    )
    convergence_ok = len(convergence_rows) == 3 and all(
        row["passed"] == "true" for row in convergence_rows
    )
    residual_ok = len(residual_rows) == 3 and all(
        row["passed"] == "true" for row in residual_rows
    )
    interpretation_ok = summary30["center_exceeds_both_neighbors"] == "true"
    gates = {
        "Gate A: scope_isolation_gate": scope_ok,
        "Gate B: additional_angle_selection_gate": selection_ok,
        "Gate C: zero_root_solve_gate": zero_solve_ok,
        "Gate D: energy_reconstruction_gate": energy_ok,
        "Gate E: integration_convergence_gate": convergence_ok,
        "Gate F: mode_residual_gate": residual_ok,
        "Gate G: additional_angle_interpretation_gate": interpretation_ok,
        "Gate H: prior_results_preservation_gate": prior_results_unchanged,
        "Gate I: article_figure_consistency_gate": figure_consistency_pass,
    }
    gates["Gate J: angle_extension_readiness_gate"] = all(gates.values())
    status = (
        "additional_angle_hypothesis_supported"
        if interpretation_ok
        else "additional_angle_hypothesis_not_supported"
    )
    return gates, status


def write_figure_selection_notes(
    output_dir: Path,
    *,
    scientific_status: str,
) -> Path:
    diagnostic_note = (
        ""
        if scientific_status == "additional_angle_hypothesis_supported"
        else " Поскольку дополнительный пример не прошёл интерпретационный gate, расширенный рисунок следует считать диагностическим."
    )
    lines = [
        "# Внутренние заметки по выбору рисунка",
        "",
        "`article_energy_angle_comparison.png` изолирует влияние угла: epsilon_0=0.030, mu=0.5 и eta=0.5 одинаковы, а beta изменяется от 15 до 30 градусов.",
        "",
        "`article_energy_examples_extended.png` показывает оба угловых примера вместе с ранее рассчитанным epsilon_0=0.050, beta=15 примером и поэтому даёт более широкий контекст, но смешивает изменение угла с изменением epsilon_0 и геометрии.",
        "",
        "Случай beta=0 сохраняется отдельно как контрольный предел разделения движений и не включён ни в один основной рисунок.",
        "",
        "Пример epsilon_0=0.050 имеет s_max>0.1 и относится к расширенному диапазону одномерной модели." + diagnostic_note,
        "",
        "После визуальной проверки для статьи рекомендуется компактный `article_energy_angle_comparison.png`, поскольку он представляет контролируемое сравнение при изменении только beta. Расширенный вариант целесообразно сохранить как альтернативу для более широкого контекста.",
    ]
    return _atomic_text(output_dir / "figure_selection_notes.md", "\n".join(lines) + "\n")


def write_additional_angle_report(
    output_dir: Path,
    selected: Mapping[str, Any],
    energy_rows: Sequence[Mapping[str, object]],
    comparison_rows: Sequence[Mapping[str, object]],
    convergence_rows: Sequence[Mapping[str, object]],
    residual_rows: Sequence[Mapping[str, object]],
    summary30: Mapping[str, object],
    counters: Mapping[str, int],
    gates: Mapping[str, bool],
    scientific_status: str,
    selection_context: Mapping[str, object],
) -> Path:
    lines = [
        "# Внутренний отчёт: дополнительный угловой энергетический пример",
        "",
        "Этот файл фиксирует вычислительный пилот и не является текстом раздела статьи.",
        "",
        f"Scientific scope: `{SCIENTIFIC_SCOPE}`.",
        "",
        "## Frequency-only selection",
        "",
        f"Выбран `{selected['case_id']}`: epsilon_0={float(selected['epsilon_0']):.3f}, beta={float(selected['beta_deg']):g} deg, mu={float(selected['mu']):g}, eta={float(selected['eta']):g}, k={int(selected['central_sorted_index'])}.",
        f"delta_4={float(selected['delta_minus']):.16g}, delta_5={float(selected['delta_center']):.16g}, delta_6={float(selected['delta_plus']):.16g}, R_5={float(selected['R_k']):.16g}, D_5={float(selected['D_k']):.16g}, s_max={float(selected['s_max']):.16g}.",
        f"Исходный ranking содержал {selection_context['ranking_source_row_count']} строк; после дедупликации по case_id осталось {selection_context['ranking_unique_case_id_count']}, удалено повторов: {selection_context['ranking_duplicate_rows_removed']}. Допустимых beta>=30, s_max<=0.1 кандидатов: {selection_context['additional_angle_candidate_count']}.",
        "",
        "## Три сохранённые формы",
        "",
        *_markdown_table(
            ("k", "delta_f, %", "U_a", "U_b", "U_s", "U_total", "chi_a", "chi_b+s"),
            [
                (
                    str(row["sorted_index"]),
                    f"{100.0 * float(row['delta_f']):.6g}",
                    f"{float(row['U_axial_Timo']):.8g}",
                    f"{float(row['U_bending_Timo']):.8g}",
                    f"{float(row['U_shear_Timo']):.8g}",
                    f"{float(row['U_total_Timo']):.8g}",
                    f"{float(row['chi_axial_Timo']):.8g}",
                    f"{float(row['chi_bending_shear_Timo']):.8g}",
                )
                for row in energy_rows
            ],
        ),
        "",
        f"Контраст beta=30: Q_left={float(summary30['Q_left']):.8g}, Q_right={float(summary30['Q_right']):.8g}.",
        "",
        "## Изменение beta=15 -> 30",
        "",
        *_markdown_table(
            ("k", "Delta delta_f", "Delta chi_a", "Delta chi_b+s"),
            [
                (
                    str(row["sorted_index"]),
                    f"{float(row['delta_f_change_beta30_minus_beta15']):.8g}",
                    f"{float(row['chi_axial_change_beta30_minus_beta15']):.8g}",
                    f"{float(row['chi_bending_shear_change_beta30_minus_beta15']):.8g}",
                )
                for row in comparison_rows
            ],
        ),
        "",
        "## Numerical checks",
        "",
        f"Максимумы 801/1601: delta_chi_a={max(float(row['delta_chi_a']) for row in convergence_rows):.6e}, delta_chi_b={max(float(row['delta_chi_b']) for row in convergence_rows):.6e}, delta_chi_s={max(float(row['delta_chi_s']) for row in convergence_rows):.6e}; tolerance={INTEGRATION_TOL:.1e}.",
        f"Максимумы residual: sigma_min={max(float(row['smallest_singular_value']) for row in residual_rows):.6e}, relative null={max(float(row['relative_null_residual']) for row in residual_rows):.6e}, clamp={max(float(row['max_abs_clamp_gap']) for row in residual_rows):.6e}, joint kinematic={max(float(row['max_abs_kinematic_gap']) for row in residual_rows):.6e}, joint force={max(float(row['max_abs_force_gap']) for row in residual_rows):.6e}, sign flip={max(float(row['sign_flip_max_abs_delta']) for row in residual_rows):.6e}.",
        "",
        "## Gates",
        "",
        *[f"- {name}: {'PASS' if passed else 'FAIL'}" for name, passed in gates.items()],
        "",
        f"Научный статус: `{scientific_status}`.",
        "",
        "## Operation counters",
        "",
        *[f"- `{name}`: {value}" for name, value in counters.items()],
        "",
        "Новые корни не искались и не уточнялись. Scientific formulas, basis и tolerances не изменялись.",
    ]
    return _atomic_text(output_dir / "report.md", "\n".join(lines) + "\n")


def run_additional_angle_pilot(
    finalization_dir: Path = DEFAULT_FINALIZATION_DIR,
    compact_dir: Path = DEFAULT_COMPACT_DIR,
    beta0_control_dir: Path = DEFAULT_BETA0_CONTROL_DIR,
    coupled_control_dir: Path = DEFAULT_COUPLED_CONTROL_DIR,
    output_dir: Path = DEFAULT_ADDITIONAL_ANGLE_OUTPUT_DIR,
    *,
    selection_only: bool = False,
) -> dict[str, object]:
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    beta0_control_dir = beta0_control_dir.resolve()
    coupled_control_dir = coupled_control_dir.resolve()
    beta0_hashes_before = _directory_hashes(beta0_control_dir)
    coupled_hashes_before = _directory_hashes(coupled_control_dir)

    candidates, selected, input_hashes, selection_context = select_additional_angle_candidate(
        finalization_dir,
        compact_dir,
        coupled_control_dir,
    )
    # The selected case is durably fixed before any prior energy table is parsed.
    write_additional_angle_selection_outputs(output_dir, candidates, selected)
    counters = operation_counters()
    if selection_only:
        compact_inputs_unchanged = all(
            _sha256(Path(path_text)) == before_hash
            for path_text, before_hash in input_hashes.items()
        )
        beta0_unchanged = beta0_hashes_before == _directory_hashes(beta0_control_dir)
        coupled_unchanged = coupled_hashes_before == _directory_hashes(coupled_control_dir)
        _atomic_csv(output_dir / "operation_counts.csv", OPERATION_COUNTER_FIELDS, [counters])
        _atomic_json(
            output_dir / "run_metadata.json",
            {
                "scientific_scope": SCIENTIFIC_SCOPE,
                "run_mode": "additional_angle_selection_only",
                "selected_case_id": selected["case_id"],
                "selection_context": selection_context,
                "energy_data_opened": False,
                "compact_inputs_unchanged": compact_inputs_unchanged,
                "beta0_control_unchanged": beta0_unchanged,
                "coupled_control_unchanged": coupled_unchanged,
                "operation_counts": counters,
            },
        )
        return {
            "selection_only": True,
            "selected_case_id": selected["case_id"],
            "selection_context": selection_context,
            "compact_inputs_unchanged": compact_inputs_unchanged,
            "prior_results_unchanged": beta0_unchanged and coupled_unchanged,
            "operation_counts": counters,
            "output_dir": output_dir.as_posix(),
        }

    energy_rows_raw, convergence_rows, residual_rows = reconstruct_selected_modes(
        [selected], counters
    )
    energy_rows = [
        {field: row[field] for field in COUPLED_ENERGY_FIELDS}
        for row in energy_rows_raw
    ]
    _atomic_csv(
        output_dir / "additional_angle_energy_triplet.csv",
        COUPLED_ENERGY_FIELDS,
        energy_rows,
    )
    _atomic_csv(
        output_dir / "integration_convergence.csv",
        CONVERGENCE_FIELDS,
        convergence_rows,
    )
    _atomic_csv(
        output_dir / "mode_residual_audit.csv",
        RESIDUAL_FIELDS,
        residual_rows,
    )
    certificate = selected["certificate"]
    eb_roots = [float(value) for value in certificate["spectra"]["Euler-Bernoulli"]["roots"]]
    timo_roots = [float(value) for value in certificate["spectra"]["Timoshenko"]["roots"]]
    root_immutability_pass = all(
        float(row["Lambda_EB"]) == eb_roots[int(row["sorted_index"]) - 1]
        and float(row["Lambda_Timo"]) == timo_roots[int(row["sorted_index"]) - 1]
        for row in energy_rows
    )

    # Prior energy values are opened only after selected_additional_angle_case.csv exists.
    beta15_rows = read_coupled_control_energy(coupled_control_dir)
    combined_rows = build_combined_article_energy_rows(beta15_rows, energy_rows)
    _atomic_csv(
        output_dir / "combined_article_energy_triplets.csv",
        COUPLED_ENERGY_FIELDS,
        combined_rows,
    )
    comparison_rows, summary15, summary30 = build_beta15_vs_beta30_comparison(combined_rows)
    _atomic_csv(
        output_dir / "beta15_vs_beta30_comparison.csv",
        ANGLE_COMPARISON_FIELDS,
        comparison_rows,
    )
    compact_figure, extended_figure = write_additional_angle_figures(
        output_dir, combined_rows
    )
    write_additional_angle_tables(output_dir, combined_rows)

    compact_rows = [row for row in combined_rows if _same_float(row["epsilon_0"], 0.030)]
    compact_source = output_dir / "article_energy_angle_comparison_source.csv"
    extended_source = output_dir / "article_energy_examples_extended_source.csv"
    figure_consistency_pass = (
        compact_figure.is_file()
        and extended_figure.is_file()
        and _article_figure_source_consistent(compact_source, compact_rows)
        and _article_figure_source_consistent(extended_source, combined_rows)
        and not any(output_dir.rglob("*.pdf"))
    )

    compact_inputs_unchanged = all(
        _sha256(Path(path_text)) == before_hash
        for path_text, before_hash in input_hashes.items()
    )
    beta0_hashes_after = _directory_hashes(beta0_control_dir)
    coupled_hashes_after = _directory_hashes(coupled_control_dir)
    beta0_unchanged = beta0_hashes_before == beta0_hashes_after
    coupled_unchanged = coupled_hashes_before == coupled_hashes_after
    prior_results_unchanged = beta0_unchanged and coupled_unchanged
    gates, scientific_status = evaluate_additional_angle_gates(
        selected,
        energy_rows,
        convergence_rows,
        residual_rows,
        summary30,
        counters,
        selection_context,
        compact_inputs_unchanged=compact_inputs_unchanged,
        prior_results_unchanged=prior_results_unchanged,
        root_immutability_pass=root_immutability_pass,
        figure_consistency_pass=figure_consistency_pass,
    )
    _atomic_csv(output_dir / "operation_counts.csv", OPERATION_COUNTER_FIELDS, [counters])
    _atomic_csv(
        output_dir / "gate_summary.csv",
        ("gate", "status"),
        ({"gate": name, "status": "PASS" if passed else "FAIL"} for name, passed in gates.items()),
    )
    write_figure_selection_notes(output_dir, scientific_status=scientific_status)
    write_additional_angle_report(
        output_dir,
        selected,
        energy_rows,
        comparison_rows,
        convergence_rows,
        residual_rows,
        summary30,
        counters,
        gates,
        scientific_status,
        selection_context,
    )
    generated_files = sorted(
        path.name for path in output_dir.iterdir() if path.is_file() and path.name != "run_metadata.json"
    )
    _atomic_json(
        output_dir / "run_metadata.json",
        {
            "scientific_scope": SCIENTIFIC_SCOPE,
            "run_mode": "additional_angle_pilot",
            "selected_case_id": selected["case_id"],
            "selection_context": selection_context,
            "energy_row_count": len(energy_rows),
            "combined_energy_row_count": len(combined_rows),
            "integration_grids": list(INTEGRATION_GRIDS),
            "energy_sum_tolerance": ENERGY_SUM_TOL,
            "integration_tolerance": INTEGRATION_TOL,
            "mode_residual_tolerance": MODE_RESIDUAL_TOL,
            "contrast_numerical_floor": CONTRAST_NUMERICAL_FLOOR,
            "root_immutability_pass": root_immutability_pass,
            "compact_inputs_unchanged": compact_inputs_unchanged,
            "beta0_control_unchanged": beta0_unchanged,
            "coupled_control_unchanged": coupled_unchanged,
            "beta0_control_file_hashes_before": beta0_hashes_before,
            "beta0_control_file_hashes_after": beta0_hashes_after,
            "coupled_control_file_hashes_before": coupled_hashes_before,
            "coupled_control_file_hashes_after": coupled_hashes_after,
            "beta15_summary": summary15,
            "beta30_summary": summary30,
            "operation_counts": counters,
            "gates": gates,
            "scientific_status": scientific_status,
            "generated_files": generated_files,
        },
    )
    if not compact_inputs_unchanged or not prior_results_unchanged:
        raise RuntimeError("immutable compact/prior-result integrity check failed")
    return {
        "selection_only": False,
        "selected_case_id": selected["case_id"],
        "energy_row_count": len(energy_rows),
        "combined_energy_row_count": len(combined_rows),
        "scientific_status": scientific_status,
        "gates": gates,
        "compact_inputs_unchanged": compact_inputs_unchanged,
        "prior_results_unchanged": prior_results_unchanged,
        "root_immutability_pass": root_immutability_pass,
        "operation_counts": counters,
        "output_dir": output_dir.as_posix(),
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Select two finalized isotropic-circular EB/Timoshenko frequency triplets "
            "and reconstruct only their six Timoshenko energy partitions."
        )
    )
    parser.add_argument("--finalization-dir", type=Path, default=DEFAULT_FINALIZATION_DIR)
    parser.add_argument("--compact-dir", type=Path, default=DEFAULT_COMPACT_DIR)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--coupled-angle-pilot",
        action="store_true",
        help=(
            "Run the beta>=15 degree frequency-only selection and reconstruct the "
            "two fixed coupled-angle triplets; the default beta=0 pilot is unchanged."
        ),
    )
    parser.add_argument(
        "--additional-angle-pilot",
        action="store_true",
        help=(
            "Select the epsilon_0=0.030, beta>=30, s_max<=0.1 extension and "
            "reconstruct only its saved k=4,5,6 Timoshenko modes."
        ),
    )
    parser.add_argument(
        "--beta0-control-dir",
        type=Path,
        default=DEFAULT_BETA0_CONTROL_DIR,
        help="Immutable beta=0 pilot directory used only as a zero-solve control reference.",
    )
    parser.add_argument(
        "--coupled-control-dir",
        type=Path,
        default=DEFAULT_COUPLED_CONTROL_DIR,
        help="Immutable beta=15 coupled pilot used after additional-angle selection.",
    )
    parser.add_argument("--selection-only", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.additional_angle_pilot and args.coupled_angle_pilot:
        raise ValueError("choose only one of --coupled-angle-pilot and --additional-angle-pilot")
    if args.additional_angle_pilot:
        output_dir = args.output_dir or DEFAULT_ADDITIONAL_ANGLE_OUTPUT_DIR
        result = run_additional_angle_pilot(
            finalization_dir=args.finalization_dir,
            compact_dir=args.compact_dir,
            beta0_control_dir=args.beta0_control_dir,
            coupled_control_dir=args.coupled_control_dir,
            output_dir=output_dir,
            selection_only=args.selection_only,
        )
    elif args.coupled_angle_pilot:
        output_dir = args.output_dir or DEFAULT_COUPLED_OUTPUT_DIR
        result = run_coupled_pilot(
            finalization_dir=args.finalization_dir,
            compact_dir=args.compact_dir,
            beta0_control_dir=args.beta0_control_dir,
            output_dir=output_dir,
            selection_only=args.selection_only,
        )
    else:
        output_dir = args.output_dir or DEFAULT_OUTPUT_DIR
        result = run_pilot(
            args.finalization_dir,
            args.compact_dir,
            output_dir,
            selection_only=args.selection_only,
        )
    print(json.dumps(result, ensure_ascii=False, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
