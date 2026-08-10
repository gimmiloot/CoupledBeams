from __future__ import annotations

"""Finite-grid energy check for local minima of same-index EB/Timo mismatch.

Phase A reads only finalized compact frequency certificates and durably fixes
the selected triplets. Phase B verifies that selection-file hash before it
reconstructs the deduplicated saved-root Timoshenko modes. No root search,
branch tracking, MAC, FEM, or anisotropic implementation is used.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import statistics
import sys
import time
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_article_longitudinal_energy_examples as PILOT,
)


SCIENTIFIC_SCOPE = "isotropic_circular_coupled_rods_eb_timoshenko"
SCHEMA_VERSION = "article_energy_local_minima_grid_v1"
BASE_EPSILON_LEVELS = (0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.060)
K_MAX = 10
LOCAL_MIN_K = tuple(range(2, 10))
DELTA_0 = 0.10
THINNESS_LIMIT = 0.10
INTEGRATION_GRIDS = (801, 1601)
NUMERICAL_FLOOR = 1.0e-30
DEFAULT_SOURCE_DIR = (
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
    / "energy_local_minima_grid_v1"
)
PRIOR_RESULT_DIRS = {
    "energy_triplet_pilot_v1": REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "energy_triplet_pilot_v1",
    "energy_triplet_coupled_pilot_v2": REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "energy_triplet_coupled_pilot_v2",
    "energy_triplet_angle_extension_v3": REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "energy_triplet_angle_extension_v3",
}
DISCOVERY_SELECTION_FILES = (
    PRIOR_RESULT_DIRS["energy_triplet_pilot_v1"] / "selected_cases.csv",
    PRIOR_RESULT_DIRS["energy_triplet_coupled_pilot_v2"] / "selected_coupled_cases.csv",
    PRIOR_RESULT_DIRS["energy_triplet_angle_extension_v3"]
    / "selected_additional_angle_case.csv",
)

FREQUENCY_FIELDS = (
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "s_max",
    "k",
    "delta_f_k_minus_1",
    "delta_f_k",
    "delta_f_k_plus_1",
    "D_f",
    "threshold_any",
    "threshold_both",
    "eb_first10_quality_flag",
    "timo_first10_quality_flag",
    "eb_triplet_quality_flags",
    "timo_triplet_quality_flags",
    "cluster_flag",
    "eb_cluster_flag",
    "timo_cluster_flag",
    "discovery_case",
    "strict_family_excluded",
    "source_frequency_file",
    "source_frequency_sha256",
)

COVERAGE_FIELDS = (
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "s_max",
    "coverage_status",
    "coverage_reason",
    "eb_root_count",
    "timo_root_count",
    "eb_confirmed_first10",
    "timo_confirmed_first10",
    "delta_count_first10",
    "source_frequency_file",
    "source_frequency_sha256",
)

DISCOVERY_FIELDS = (
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "central_sorted_index",
    "source_dataset",
    "source_selection_file",
    "source_selection_sha256",
)

DISCOVERY_RESULT_FIELDS = (
    *DISCOVERY_FIELDS,
    "first10_coverage_status",
    "first10_coverage_reason",
    "frequency_selected_triplet_count",
    "energy_quality_pass_triplet_count",
)

UNIQUE_MODE_FIELDS = (
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "s_max",
    "sorted_k",
    "Lambda_EB",
    "Lambda_Timo",
    "delta_f",
    "source_frequency_file",
    "source_frequency_sha256",
    "mode_cache_key",
)

MODE_ENERGY_FIELDS = (
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "s_max",
    "sorted_k",
    "Lambda_EB",
    "Lambda_Timo",
    "delta_f",
    "n_points",
    "U_a",
    "U_b",
    "U_s",
    "U_total",
    "chi_a",
    "chi_b",
    "chi_s",
    "chi_bs",
    "U_a_rod1",
    "U_b_rod1",
    "U_s_rod1",
    "U_a_rod2",
    "U_b_rod2",
    "U_s_rod2",
    "rod1_energy_fraction",
    "rod2_energy_fraction",
    "energy_sum_residual",
    "energy_sum_pass",
)

MODE_CONVERGENCE_FIELDS = (
    "case_id",
    "sorted_k",
    "chi_a_801",
    "chi_a_1601",
    "delta_chi_a_abs",
    "delta_chi_a_rel",
    "delta_chi_b_abs",
    "delta_chi_s_abs",
    "integration_tolerance",
    "integration_pass",
    "smallest_singular_value",
    "second_smallest_singular_value",
    "singular_value_ratio",
    "relative_null_residual",
    "max_abs_clamp_gap",
    "max_abs_kinematic_gap",
    "max_abs_force_gap",
    "sign_flip_max_abs_delta",
    "svd_pass",
    "null_residual_pass",
    "clamp_pass",
    "joint_kinematic_pass",
    "joint_force_pass",
    "sign_flip_pass",
    "quality_pass",
    "null_vector_source",
    "warnings",
)

TRIPLET_FIELDS = (
    *FREQUENCY_FIELDS,
    "chi_a_k_minus_1_801",
    "chi_a_k_801",
    "chi_a_k_plus_1_801",
    "chi_a_k_minus_1",
    "chi_a_k",
    "chi_a_k_plus_1",
    "chi_bs_k_minus_1",
    "chi_bs_k",
    "chi_bs_k_plus_1",
    "D_a_801",
    "D_a",
    "axial_above_neighbor_mean",
    "axial_strict_local_max_801",
    "axial_strict_local_max",
    "grid_classification_changed",
    "energy_quality_pass",
    "energy_failure_reason",
    "primary_confirmatory",
    "strict_family_holdout",
    "beta0_control",
    "extended_range",
    "cluster_affected",
    "discovery_cohort",
)

CONTRADICTING_FIELDS = (*TRIPLET_FIELDS, "contradiction_reasons")

SUMMARY_FIELDS = (
    "cohort",
    "group_field",
    "group_value",
    "n_loc",
    "n_geometries",
    "n_mu_eta",
    "n_D_a_positive",
    "n_axial_strict_local_max",
    "ratio_D_a_positive_percent",
    "ratio_axial_strict_local_max_percent",
    "n_threshold_any",
    "n_threshold_any_axial_strict_local_max",
    "n_threshold_both",
    "n_threshold_both_axial_strict_local_max",
    "D_f_min",
    "D_f_q1",
    "D_f_median",
    "D_f_q3",
    "D_f_max",
    "D_a_min",
    "D_a_q1",
    "D_a_median",
    "D_a_q3",
    "D_a_max",
    "strongest_support_case_id",
    "strongest_support_k",
    "strongest_support_D_a",
    "strongest_contradiction_case_id",
    "strongest_contradiction_k",
    "strongest_contradiction_D_a",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def directory_fingerprint(path: Path) -> str:
    root = path.resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"immutable prior-result directory is missing: {root}")
    digest = hashlib.sha256()
    for item in sorted(candidate for candidate in root.rglob("*") if candidate.is_file()):
        digest.update(item.relative_to(root).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256_file(item).encode("ascii"))
        digest.update(b"\0")
        digest.update(str(item.stat().st_size).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def prior_fingerprints() -> dict[str, str]:
    return {name: directory_fingerprint(path) for name, path in PRIOR_RESULT_DIRS.items()}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def csv_bytes(fields: Sequence[str], rows: Iterable[Mapping[str, object]], *, delimiter: str = ",") -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream,
        fieldnames=fields,
        delimiter=delimiter,
        lineterminator="\n",
        extrasaction="raise",
    )
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def atomic_bytes(path: Path, payload: bytes) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(payload)
    os.replace(temporary, path)
    return path


def atomic_csv(
    path: Path,
    fields: Sequence[str],
    rows: Iterable[Mapping[str, object]],
    *,
    delimiter: str = ",",
) -> Path:
    return atomic_bytes(path, csv_bytes(fields, rows, delimiter=delimiter))


def atomic_text(path: Path, text: str) -> Path:
    return atomic_bytes(path, text.encode("utf-8"))


def atomic_json(path: Path, payload: Mapping[str, object]) -> Path:
    return atomic_text(
        path,
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True, default=str) + "\n",
    )


def bool_text(value: object) -> str:
    return "true" if bool(value) else "false"


def true_text(value: object) -> bool:
    return str(value).strip().lower() == "true"


def same_float(left: object, right: float, tolerance: float = 1.0e-12) -> bool:
    return abs(float(left) - float(right)) <= tolerance


def optional_float_text(value: object, format_spec: str) -> str:
    if value is None or str(value) == "":
        return "n/a"
    return format(float(value), format_spec)


def frequency_local_characteristics(
    left: float,
    center: float,
    right: float,
) -> dict[str, object]:
    """Return the frequency-only triplet characteristics used by Phase A."""

    return {
        "strict_local_minimum": center < left and center < right,
        "D_f": 0.5 * (left + right) - center,
        "threshold_any": (
            center <= DELTA_0 and (left > DELTA_0 or right > DELTA_0)
        ),
        "threshold_both": center <= DELTA_0 and left > DELTA_0 and right > DELTA_0,
    }


def energy_local_characteristics(
    left: float,
    center: float,
    right: float,
) -> dict[str, object]:
    """Return threshold-free local characteristics of the axial contribution."""

    d_a = center - 0.5 * (left + right)
    return {
        "D_a": d_a,
        "axial_above_neighbor_mean": d_a > 0.0,
        "axial_strict_local_max": center > max(left, right),
    }


def cohort_flags(
    selected: Mapping[str, object],
    *,
    energy_quality_pass: bool,
) -> dict[str, bool]:
    """Apply the fixed confirmatory/control cohort policy to one triplet."""

    beta0 = same_float(selected["beta_deg"], 0.0)
    extended = float(selected["s_max"]) > THINNESS_LIMIT
    cluster = true_text(selected["cluster_flag"])
    discovery = true_text(selected["discovery_case"])
    primary = (
        energy_quality_pass
        and not beta0
        and not extended
        and not cluster
        and not discovery
    )
    strict = primary and not (
        same_float(selected["mu"], 0.5) and same_float(selected["eta"], 0.5)
    )
    return {
        "primary_confirmatory": primary,
        "strict_family_holdout": strict,
        "beta0_control": beta0 and energy_quality_pass,
        "extended_range": extended and energy_quality_pass,
        "cluster_affected": cluster,
        "discovery_cohort": discovery,
    }


def resolve_repo_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path


def load_certificate(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("scientific_scope") != SCIENTIFIC_SCOPE:
        raise ValueError(f"non-isotropic scientific scope in {path}")
    if payload.get("schema_version") != PILOT.CERTIFICATE_SCHEMA:
        raise ValueError(f"unexpected compact certificate schema in {path}")
    return payload


def load_discovery_cases() -> tuple[list[dict[str, object]], dict[str, str]]:
    rows: list[dict[str, object]] = []
    hashes: dict[str, str] = {}
    seen: set[tuple[str, str]] = set()
    for path in DISCOVERY_SELECTION_FILES:
        if not path.is_file():
            continue
        source_hash = sha256_file(path)
        hashes[path.resolve().as_posix()] = source_hash
        for source in read_csv(path):
            case_id = str(source["case_id"])
            source_dataset = path.parent.name
            key = (case_id, source_dataset)
            if key in seen:
                continue
            seen.add(key)
            rows.append(
                {
                    "case_id": case_id,
                    "epsilon_0": float(source["epsilon_0"]),
                    "beta_deg": float(source["beta_deg"]),
                    "mu": float(source["mu"]),
                    "eta": float(source["eta"]),
                    "central_sorted_index": int(float(source["central_sorted_index"])),
                    "source_dataset": source_dataset,
                    "source_selection_file": path.resolve().as_posix(),
                    "source_selection_sha256": source_hash,
                }
            )
    if not rows:
        rows = [
            {
                "case_id": "",
                "epsilon_0": 0.030,
                "beta_deg": beta,
                "mu": 0.5,
                "eta": 0.5,
                "central_sorted_index": 5,
                "source_dataset": "fallback_geometry",
                "source_selection_file": "",
                "source_selection_sha256": "",
            }
            for beta in (15.0, 30.0)
        ]
    rows.sort(
        key=lambda row: (
            float(row["epsilon_0"]),
            float(row["beta_deg"]),
            float(row["mu"]),
            float(row["eta"]),
            str(row["case_id"]),
            str(row["source_dataset"]),
        )
    )
    return rows, hashes


def discovery_case_ids(rows: Sequence[Mapping[str, object]]) -> set[str]:
    return {str(row["case_id"]) for row in rows if str(row["case_id"]).strip()}


def quality_by_rank(model: Mapping[str, Any]) -> dict[int, Mapping[str, Any]]:
    return {int(row["rank"]): row for row in model.get("quality_by_rank", [])}


def confirmed_first_ten(model: Mapping[str, Any]) -> tuple[bool, int, str]:
    by_rank = quality_by_rank(model)
    confirmed = sum(
        by_rank.get(rank, {}).get("status") == "matrix_confirmed_resolved"
        for rank in range(1, K_MAX + 1)
    )
    flag = ";".join(
        f"{rank}:{by_rank.get(rank, {}).get('status', 'missing')}"
        for rank in range(1, K_MAX + 1)
    )
    return confirmed == K_MAX, confirmed, flag


def cluster_ranks(model: Mapping[str, Any]) -> set[int]:
    ranks = {
        int(row["sorted_index"])
        for row in model.get("clusters", [])
        if str(row.get("sorted_index", "")).strip()
    }
    ranks.update(
        int(rank)
        for rank, row in quality_by_rank(model).items()
        if str(row.get("cluster_id", "")).strip()
    )
    return ranks


def source_paths(source_dir: Path, compact_dir: Path) -> tuple[Path, Path, Path]:
    final_path = source_dir.resolve() / "final_case_certificates.csv"
    final_metadata = source_dir.resolve() / "run_metadata.json"
    compact_index = compact_dir.resolve() / "compact_index.csv"
    for path in (final_path, final_metadata, compact_index):
        if not path.is_file():
            raise FileNotFoundError(f"required compact source is missing: {path}")
    metadata = json.loads(final_metadata.read_text(encoding="utf-8"))
    if metadata.get("scientific_scope") != SCIENTIFIC_SCOPE:
        raise ValueError("final compact source has the wrong scientific scope")
    if sha256_file(compact_index) != metadata.get("compact_index_sha256"):
        raise ValueError("compact index SHA does not match finalization metadata")
    return final_path, final_metadata, compact_index


def _triplet_quality_text(model: Mapping[str, Any], ranks: Sequence[int]) -> str:
    by_rank = quality_by_rank(model)
    return ";".join(
        f"{rank}:{by_rank.get(rank, {}).get('status', 'missing')}"
        for rank in ranks
    )


def select_frequency_triplets(
    source_dir: Path,
    compact_dir: Path,
    output_dir: Path,
    *,
    resume: bool,
    smoke: bool,
    smoke_cases: int,
) -> dict[str, object]:
    """Phase A: select strict local minima without reading any energy table."""

    final_path, final_metadata_path, compact_index_path = source_paths(source_dir, compact_dir)
    final_rows = read_csv(final_path)
    compact_rows = read_csv(compact_index_path)
    if {row["scientific_scope"] for row in final_rows} != {SCIENTIFIC_SCOPE}:
        raise ValueError("mixed scientific scopes in final case table")
    compact_by_id = {row["case_id"]: row for row in compact_rows}
    if len(compact_by_id) != len(compact_rows):
        raise ValueError("duplicate case IDs in compact index")
    base_rows = [
        row
        for row in final_rows
        if any(same_float(row["epsilon_0"], epsilon) for epsilon in BASE_EPSILON_LEVELS)
    ]
    if len(base_rows) != 1552 or len({row["case_id"] for row in base_rows}) != 1552:
        raise ValueError("main article grid must contain 1552 unique cases")

    discovery_rows, discovery_hashes = load_discovery_cases()
    discovery_ids = discovery_case_ids(discovery_rows)
    selected: list[dict[str, object]] = []
    coverage: list[dict[str, object]] = []
    certificate_hashes: dict[str, str] = {}
    eligible_case_ids: list[str] = []

    base_rows.sort(
        key=lambda row: (
            float(row["epsilon_0"]),
            float(row["beta"]),
            float(row["mu"]),
            float(row["eta"]),
            str(row["case_id"]),
        )
    )
    for final_row in base_rows:
        case_id = str(final_row["case_id"])
        compact_row = compact_by_id.get(case_id)
        reason = ""
        certificate_path: Path | None = None
        certificate_sha = ""
        eb_count = 0
        timo_count = 0
        eb_confirmed = 0
        timo_confirmed = 0
        delta_count = 0
        certificate: dict[str, Any] | None = None
        if compact_row is None:
            reason = "missing_compact_index_row"
        elif compact_row.get("scientific_scope") != SCIENTIFIC_SCOPE:
            reason = "wrong_scientific_scope"
        else:
            certificate_path = resolve_repo_path(compact_row["certificate_path"]).resolve()
            if not certificate_path.is_file():
                reason = "missing_compact_certificate"
            else:
                certificate_sha = sha256_file(certificate_path)
                certificate_hashes[certificate_path.as_posix()] = certificate_sha
                certificate = load_certificate(certificate_path)
                spectra = certificate.get("spectra", {})
                eb = spectra.get("Euler-Bernoulli", {})
                timo = spectra.get("Timoshenko", {})
                eb_roots = [float(value) for value in eb.get("roots", [])]
                timo_roots = [float(value) for value in timo.get("roots", [])]
                eb_count = len(eb_roots)
                timo_count = len(timo_roots)
                eb_ok, eb_confirmed, eb_flag = confirmed_first_ten(eb)
                timo_ok, timo_confirmed, timo_flag = confirmed_first_ten(timo)
                delta_payload = spectra.get("delta_f", {})
                delta_count = sum(str(rank) in delta_payload for rank in range(1, K_MAX + 1))
                geometry = certificate.get("geometry", {})
                if any(
                    not same_float(geometry[name], final_row[source_name])
                    for name, source_name in (
                        ("epsilon_0", "epsilon_0"),
                        ("beta_deg", "beta"),
                        ("mu", "mu"),
                        ("eta", "eta"),
                    )
                ):
                    reason = "geometry_mismatch"
                elif eb_count < K_MAX:
                    reason = "fewer_than_10_EB_roots"
                elif timo_count < K_MAX:
                    reason = "fewer_than_10_Timoshenko_roots"
                elif not eb_ok:
                    reason = "EB_first10_quality_not_confirmed"
                elif not timo_ok:
                    reason = "Timoshenko_first10_quality_not_confirmed"
                elif delta_count < K_MAX:
                    reason = "fewer_than_10_saved_delta_f_values"
                elif not bool(certificate.get("diagnostics", {}).get("unresolved_cluster_below_guard")):
                    reason = "unresolved_inventory_below_guard"
                else:
                    for rank in range(1, K_MAX + 1):
                        stored = float(delta_payload[str(rank)])
                        recalculated = abs(eb_roots[rank - 1] ** 2 - timo_roots[rank - 1] ** 2) / (
                            timo_roots[rank - 1] ** 2
                        )
                        if not math.isclose(stored, recalculated, rel_tol=1.0e-12, abs_tol=1.0e-14):
                            raise ValueError(
                                f"saved delta_f formula mismatch for {case_id}, k={rank}"
                            )

        eligible = reason == ""
        coverage.append(
            {
                "case_id": case_id,
                "epsilon_0": float(final_row["epsilon_0"]),
                "beta_deg": float(final_row["beta"]),
                "mu": float(final_row["mu"]),
                "eta": float(final_row["eta"]),
                "s_max": float(final_row["s_max"]),
                "coverage_status": "eligible_first10" if eligible else "quality_or_incomplete",
                "coverage_reason": reason,
                "eb_root_count": eb_count,
                "timo_root_count": timo_count,
                "eb_confirmed_first10": eb_confirmed,
                "timo_confirmed_first10": timo_confirmed,
                "delta_count_first10": delta_count,
                "source_frequency_file": certificate_path.as_posix() if certificate_path else "",
                "source_frequency_sha256": certificate_sha,
            }
        )
        if not eligible or certificate is None or certificate_path is None:
            continue
        eligible_case_ids.append(case_id)
        spectra = certificate["spectra"]
        eb = spectra["Euler-Bernoulli"]
        timo = spectra["Timoshenko"]
        delta_payload = spectra["delta_f"]
        eb_clusters = cluster_ranks(eb)
        timo_clusters = cluster_ranks(timo)
        eb_first10_flag = confirmed_first_ten(eb)[2]
        timo_first10_flag = confirmed_first_ten(timo)[2]
        deltas = {rank: float(delta_payload[str(rank)]) for rank in range(1, K_MAX + 1)}
        for k in LOCAL_MIN_K:
            frequency_metrics = frequency_local_characteristics(
                deltas[k - 1], deltas[k], deltas[k + 1]
            )
            if not bool(frequency_metrics["strict_local_minimum"]):
                continue
            ranks = (k - 1, k, k + 1)
            eb_cluster = bool(set(ranks) & eb_clusters)
            timo_cluster = bool(set(ranks) & timo_clusters)
            selected.append(
                {
                    "case_id": case_id,
                    "epsilon_0": float(final_row["epsilon_0"]),
                    "beta_deg": float(final_row["beta"]),
                    "mu": float(final_row["mu"]),
                    "eta": float(final_row["eta"]),
                    "s_max": float(final_row["s_max"]),
                    "k": k,
                    "delta_f_k_minus_1": deltas[k - 1],
                    "delta_f_k": deltas[k],
                    "delta_f_k_plus_1": deltas[k + 1],
                    "D_f": frequency_metrics["D_f"],
                    "threshold_any": bool_text(frequency_metrics["threshold_any"]),
                    "threshold_both": bool_text(frequency_metrics["threshold_both"]),
                    "eb_first10_quality_flag": eb_first10_flag,
                    "timo_first10_quality_flag": timo_first10_flag,
                    "eb_triplet_quality_flags": _triplet_quality_text(eb, ranks),
                    "timo_triplet_quality_flags": _triplet_quality_text(timo, ranks),
                    "cluster_flag": bool_text(eb_cluster or timo_cluster),
                    "eb_cluster_flag": bool_text(eb_cluster),
                    "timo_cluster_flag": bool_text(timo_cluster),
                    "discovery_case": bool_text(case_id in discovery_ids),
                    "strict_family_excluded": bool_text(
                        same_float(final_row["mu"], 0.5) and same_float(final_row["eta"], 0.5)
                    ),
                    "source_frequency_file": certificate_path.as_posix(),
                    "source_frequency_sha256": certificate_sha,
                }
            )

    selected.sort(
        key=lambda row: (
            float(row["epsilon_0"]),
            float(row["beta_deg"]),
            float(row["mu"]),
            float(row["eta"]),
            str(row["case_id"]),
            int(row["k"]),
        )
    )
    if smoke:
        smoke_ids: list[str] = []
        for row in selected:
            case_id = str(row["case_id"])
            if case_id not in smoke_ids:
                smoke_ids.append(case_id)
            if len(smoke_ids) >= smoke_cases:
                break
        selected = [row for row in selected if str(row["case_id"]) in set(smoke_ids)]

    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    selection_path = output_dir / "frequency_selected_triplets.csv"
    selection_payload = csv_bytes(FREQUENCY_FIELDS, selected)
    if selection_path.exists():
        if not resume:
            raise FileExistsError(
                f"frequency selection already exists; use --resume to verify it: {selection_path}"
            )
        if selection_path.read_bytes() != selection_payload:
            raise ValueError("existing frequency selection differs from deterministic Phase-A result")
    else:
        atomic_bytes(selection_path, selection_payload)
    selection_sha = sha256_file(selection_path)
    atomic_csv(output_dir / "discovery_case_ids.csv", DISCOVERY_FIELDS, discovery_rows)
    incomplete = [row for row in coverage if row["coverage_status"] != "eligible_first10"]
    atomic_csv(output_dir / "quality_and_incomplete_cases.csv", COVERAGE_FIELDS, incomplete)

    input_hashes = {
        final_path.resolve().as_posix(): sha256_file(final_path),
        final_metadata_path.resolve().as_posix(): sha256_file(final_metadata_path),
        compact_index_path.resolve().as_posix(): sha256_file(compact_index_path),
        **discovery_hashes,
        **certificate_hashes,
    }
    coverage_by_epsilon = {
        f"{epsilon:.3f}": {
            "total": sum(same_float(row["epsilon_0"], epsilon) for row in coverage),
            "eligible_first10": sum(
                same_float(row["epsilon_0"], epsilon)
                and row["coverage_status"] == "eligible_first10"
                for row in coverage
            ),
        }
        for epsilon in BASE_EPSILON_LEVELS
    }
    metadata = {
        "schema_version": SCHEMA_VERSION,
        "scientific_scope": SCIENTIFIC_SCOPE,
        "phase": "frequency_selection_fixed",
        "energy_fields_read": False,
        "source_dir": source_dir.resolve().as_posix(),
        "compact_dir": compact_dir.resolve().as_posix(),
        "base_epsilon_levels": list(BASE_EPSILON_LEVELS),
        "k_max": K_MAX,
        "local_minimum_k": list(LOCAL_MIN_K),
        "delta_0": DELTA_0,
        "selection_rule": "strict delta_f[k] < delta_f[k-1] and delta_f[k] < delta_f[k+1]",
        "D_f_rule": "0.5*(delta_f[k-1]+delta_f[k+1])-delta_f[k]",
        "main_grid_case_count": len(base_rows),
        "eligible_first10_case_count": len(eligible_case_ids),
        "quality_or_incomplete_case_count": len(incomplete),
        "frequency_selected_triplet_count": len(selected),
        "frequency_selected_geometry_count": len({row["case_id"] for row in selected}),
        "coverage_by_epsilon": coverage_by_epsilon,
        "smoke": smoke,
        "smoke_case_limit": smoke_cases if smoke else None,
        "frequency_selected_triplets_sha256": selection_sha,
        "input_file_hashes": input_hashes,
        "prior_result_fingerprints": prior_fingerprints(),
    }
    atomic_json(output_dir / "frequency_selection_metadata.json", metadata)
    return metadata


def verify_frequency_selection(output_dir: Path) -> tuple[list[dict[str, str]], dict[str, Any]]:
    selection_path = output_dir.resolve() / "frequency_selected_triplets.csv"
    metadata_path = output_dir.resolve() / "frequency_selection_metadata.json"
    if not selection_path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError("Phase B requires the fixed Phase-A selection and metadata")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("scientific_scope") != SCIENTIFIC_SCOPE:
        raise ValueError("frequency-selection metadata has the wrong scientific scope")
    actual_hash = sha256_file(selection_path)
    expected_hash = str(metadata.get("frequency_selected_triplets_sha256", ""))
    if actual_hash != expected_hash:
        raise ValueError("frequency-selected triplet SHA-256 verification failed")
    rows = read_csv(selection_path)
    if any(int(row["k"]) not in LOCAL_MIN_K for row in rows):
        raise ValueError("frequency selection contains k outside 2,...,9")
    return rows, metadata


def build_unique_modes(
    selected_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, object]]:
    certificates: dict[str, dict[str, Any]] = {}
    unique: dict[tuple[str, int], dict[str, object]] = {}
    for triplet in selected_rows:
        case_id = str(triplet["case_id"])
        certificate_path = Path(str(triplet["source_frequency_file"]))
        if not certificate_path.is_file():
            raise FileNotFoundError(f"selected certificate is missing: {certificate_path}")
        if sha256_file(certificate_path) != str(triplet["source_frequency_sha256"]):
            raise ValueError(f"selected certificate hash changed: {certificate_path}")
        certificate = certificates.setdefault(case_id, load_certificate(certificate_path))
        spectra = certificate["spectra"]
        eb_roots = [float(value) for value in spectra["Euler-Bernoulli"]["roots"]]
        timo_roots = [float(value) for value in spectra["Timoshenko"]["roots"]]
        k = int(triplet["k"])
        for sorted_k in (k - 1, k, k + 1):
            key = (case_id, sorted_k)
            if key in unique:
                continue
            payload = {
                "case_id": case_id,
                "epsilon_0": float(triplet["epsilon_0"]),
                "beta_deg": float(triplet["beta_deg"]),
                "mu": float(triplet["mu"]),
                "eta": float(triplet["eta"]),
                "s_max": float(triplet["s_max"]),
                "sorted_k": sorted_k,
                "Lambda_EB": eb_roots[sorted_k - 1],
                "Lambda_Timo": timo_roots[sorted_k - 1],
                "delta_f": float(spectra["delta_f"][str(sorted_k)]),
                "source_frequency_file": certificate_path.resolve().as_posix(),
                "source_frequency_sha256": str(triplet["source_frequency_sha256"]),
            }
            cache_basis = json.dumps(payload, sort_keys=True, separators=(",", ":"))
            payload["mode_cache_key"] = hashlib.sha256(cache_basis.encode("utf-8")).hexdigest()
            unique[key] = payload
    return [unique[key] for key in sorted(unique)]


def evaluate_energy_batch(
    mode: Mapping[str, object],
    coeff: np.ndarray,
    grids: Sequence[int],
) -> dict[int, dict[str, float]]:
    """Batch wrapper around the unchanged, array-vectorized project helper."""

    return {
        int(n_points): PILOT.TIMO.timo_energy_partition(
            float(mode["Lambda_Timo"]),
            float(mode["beta_deg"]),
            float(mode["mu"]),
            float(mode["epsilon_0"]),
            float(mode["eta"]),
            coeff=np.asarray(coeff, dtype=float),
            n_points=int(n_points),
        )
        for n_points in grids
    }


def _energy_row(
    mode: Mapping[str, object],
    energy: Mapping[str, float],
    n_points: int,
) -> dict[str, object]:
    chi_a = float(energy["axial_fraction"])
    chi_b = float(energy["bending_fraction"])
    chi_s = float(energy["shear_fraction"])
    sum_residual = abs(chi_a + chi_b + chi_s - 1.0)
    return {
        "case_id": mode["case_id"],
        "epsilon_0": float(mode["epsilon_0"]),
        "beta_deg": float(mode["beta_deg"]),
        "mu": float(mode["mu"]),
        "eta": float(mode["eta"]),
        "s_max": float(mode["s_max"]),
        "sorted_k": int(mode["sorted_k"]),
        "Lambda_EB": float(mode["Lambda_EB"]),
        "Lambda_Timo": float(mode["Lambda_Timo"]),
        "delta_f": float(mode["delta_f"]),
        "n_points": int(n_points),
        "U_a": float(energy["U_a_total"]),
        "U_b": float(energy["U_b_total"]),
        "U_s": float(energy["U_s_total"]),
        "U_total": float(energy["U_total"]),
        "chi_a": chi_a,
        "chi_b": chi_b,
        "chi_s": chi_s,
        "chi_bs": chi_b + chi_s,
        "U_a_rod1": float(energy["U_a1"]),
        "U_b_rod1": float(energy["U_b1"]),
        "U_s_rod1": float(energy["U_s1"]),
        "U_a_rod2": float(energy["U_a2"]),
        "U_b_rod2": float(energy["U_b2"]),
        "U_s_rod2": float(energy["U_s2"]),
        "rod1_energy_fraction": float(energy["rod1_energy_fraction"]),
        "rod2_energy_fraction": float(energy["rod2_energy_fraction"]),
        "energy_sum_residual": sum_residual,
        "energy_sum_pass": bool_text(sum_residual <= PILOT.ENERGY_SUM_TOL),
    }


def reconstruct_unique_mode(mode: Mapping[str, object]) -> dict[str, object]:
    lambda_timo = float(mode["Lambda_Timo"])
    beta_deg = float(mode["beta_deg"])
    mu = float(mode["mu"])
    epsilon = float(mode["epsilon_0"])
    eta = float(mode["eta"])
    sorted_k = int(mode["sorted_k"])
    try:
        mode_coefficients = PILOT.TIMO.timo_mode_coefficients(
            lambda_timo, beta_deg, mu, epsilon, eta
        )
        matrix, matrix_warnings = PILOT.TIMO.timo_coupling_matrix(
            lambda_timo, beta_deg, mu, epsilon, eta
        )
        scaled_matrix = PILOT.COMPLETE.diagnostic_scaled_matrix(matrix)
        singular_values = np.linalg.svd(scaled_matrix, compute_uv=False)
        matrix_norm = float(np.linalg.norm(scaled_matrix, ord=2))
        mode_coeff = np.asarray(mode_coefficients.coeff, dtype=float)
        denominator = float(matrix_norm * np.linalg.norm(mode_coeff))
        coefficient_residual = (
            float(np.linalg.norm(scaled_matrix @ mode_coeff)) / denominator
            if denominator > NUMERICAL_FLOOR
            else float("nan")
        )
        if math.isfinite(coefficient_residual) and coefficient_residual <= PILOT.MODE_RESIDUAL_TOL:
            coeff = mode_coeff
            null_source = "timo_mode_coefficients"
        else:
            _u, _s, vh = np.linalg.svd(scaled_matrix, full_matrices=True)
            coeff = np.asarray(vh[-1, :], dtype=float)
            coeff_norm = float(np.linalg.norm(coeff))
            if coeff_norm <= 0.0 or not math.isfinite(coeff_norm):
                raise ValueError("certificate-scaled null vector is not finite")
            coeff /= coeff_norm
            pivot = int(np.argmax(np.abs(coeff)))
            if coeff[pivot] < 0.0:
                coeff = -coeff
            null_source = "certificate_floored_row_scaling_fallback"

        energies = evaluate_energy_batch(mode, coeff, INTEGRATION_GRIDS)
        sign_energy = evaluate_energy_batch(mode, -coeff, (INTEGRATION_GRIDS[0],))[
            INTEGRATION_GRIDS[0]
        ]
        coarse = energies[INTEGRATION_GRIDS[0]]
        fine = energies[INTEGRATION_GRIDS[1]]
        coarse_fractions = (
            float(coarse["axial_fraction"]),
            float(coarse["bending_fraction"]),
            float(coarse["shear_fraction"]),
        )
        fine_fractions = (
            float(fine["axial_fraction"]),
            float(fine["bending_fraction"]),
            float(fine["shear_fraction"]),
        )
        sign_fractions = (
            float(sign_energy["axial_fraction"]),
            float(sign_energy["bending_fraction"]),
            float(sign_energy["shear_fraction"]),
        )
        deltas = tuple(abs(a - b) for a, b in zip(fine_fractions, coarse_fractions))
        sign_deltas = tuple(abs(a - b) for a, b in zip(sign_fractions, coarse_fractions))
        integration_pass = all(value <= PILOT.INTEGRATION_TOL for value in deltas)
        sign_pass = all(value <= PILOT.ENERGY_SUM_TOL for value in sign_deltas)

        fields = PILOT.TIMO.timo_mode_fields(
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
            float(np.linalg.norm(scaled_matrix @ coeff)) / denominator
            if denominator > NUMERICAL_FLOOR
            else float("nan")
        )
        smallest = float(singular_values[-1])
        second_smallest = float(singular_values[-2])
        singular_ratio = (
            smallest / second_smallest if second_smallest > NUMERICAL_FLOOR else float("nan")
        )
        rod1 = fields["rod1"]
        rod2 = fields["rod2"]
        clamp_values = (
            float(np.asarray(rod1["u"])[0]),
            float(np.asarray(rod1["w"])[0]),
            float(np.asarray(rod1["psi"])[0]),
            float(np.asarray(rod2["u"])[0]),
            float(np.asarray(rod2["w"])[0]),
            float(np.asarray(rod2["psi"])[0]),
        )
        max_clamp = max(abs(value) for value in clamp_values)
        joint = PILOT.timo_joint_continuity_row(
            epsilon=epsilon,
            beta_deg=beta_deg,
            eta=eta,
            mu=mu,
            sorted_index=sorted_k,
            Lambda=lambda_timo,
            fields=fields,
        )
        row801 = _energy_row(mode, coarse, INTEGRATION_GRIDS[0])
        row1601 = _energy_row(mode, fine, INTEGRATION_GRIDS[1])
        finite_energy = all(
            math.isfinite(float(row1601[field]))
            for field in ("U_a", "U_b", "U_s", "U_total", "chi_a", "chi_b", "chi_s")
        ) and float(row1601["U_total"]) > 0.0
        energy_sum_pass = true_text(row801["energy_sum_pass"]) and true_text(
            row1601["energy_sum_pass"]
        )
        svd_pass = math.isfinite(smallest) and smallest <= PILOT.MODE_RESIDUAL_TOL
        null_pass = (
            math.isfinite(relative_null_residual)
            and relative_null_residual <= PILOT.MODE_RESIDUAL_TOL
        )
        clamp_pass = max_clamp <= PILOT.MODE_RESIDUAL_TOL
        joint_kinematic_pass = bool(joint["pass_kinematic"])
        joint_force_pass = bool(joint["pass_force"])
        quality_pass = all(
            (
                finite_energy,
                energy_sum_pass,
                integration_pass,
                svd_pass,
                null_pass,
                clamp_pass,
                joint_kinematic_pass,
                joint_force_pass,
                sign_pass,
            )
        )
        warnings = tuple(
            dict.fromkeys(
                [
                    *mode_coefficients.warnings,
                    *matrix_warnings,
                    *fields["warnings"],
                ]
            )
        )
        convergence = {
            "case_id": mode["case_id"],
            "sorted_k": sorted_k,
            "chi_a_801": coarse_fractions[0],
            "chi_a_1601": fine_fractions[0],
            "delta_chi_a_abs": deltas[0],
            "delta_chi_a_rel": deltas[0] / max(abs(fine_fractions[0]), NUMERICAL_FLOOR),
            "delta_chi_b_abs": deltas[1],
            "delta_chi_s_abs": deltas[2],
            "integration_tolerance": PILOT.INTEGRATION_TOL,
            "integration_pass": bool_text(integration_pass),
            "smallest_singular_value": smallest,
            "second_smallest_singular_value": second_smallest,
            "singular_value_ratio": singular_ratio,
            "relative_null_residual": relative_null_residual,
            "max_abs_clamp_gap": max_clamp,
            "max_abs_kinematic_gap": float(joint["max_abs_kinematic_gap"]),
            "max_abs_force_gap": float(joint["max_abs_force_gap"]),
            "sign_flip_max_abs_delta": max(sign_deltas),
            "svd_pass": bool_text(svd_pass),
            "null_residual_pass": bool_text(null_pass),
            "clamp_pass": bool_text(clamp_pass),
            "joint_kinematic_pass": bool_text(joint_kinematic_pass),
            "joint_force_pass": bool_text(joint_force_pass),
            "sign_flip_pass": bool_text(sign_pass),
            "quality_pass": bool_text(quality_pass),
            "null_vector_source": null_source,
            "warnings": "; ".join(warnings),
        }
        return {
            "schema_version": SCHEMA_VERSION,
            "mode_cache_key": mode["mode_cache_key"],
            "mode": dict(mode),
            "status": "complete" if quality_pass else "quality_failed",
            "error": "" if quality_pass else "one_or_more_energy_or_residual_checks_failed",
            "energy_801": row801,
            "energy_1601": row1601,
            "convergence": convergence,
        }
    except Exception as exc:
        return {
            "schema_version": SCHEMA_VERSION,
            "mode_cache_key": mode["mode_cache_key"],
            "mode": dict(mode),
            "status": "reconstruction_failed",
            "error": f"{type(exc).__name__}: {exc}",
            "energy_801": None,
            "energy_1601": None,
            "convergence": None,
        }


def _mode_cache_path(output_dir: Path, mode: Mapping[str, object]) -> Path:
    return (
        output_dir.resolve()
        / "cache"
        / "modes"
        / f"{mode['case_id']}_k{int(mode['sorted_k']):02d}.json"
    )


def _load_mode_cache(path: Path, mode: Mapping[str, object]) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if (
        payload.get("schema_version") != SCHEMA_VERSION
        or payload.get("mode_cache_key") != mode["mode_cache_key"]
    ):
        return None
    return payload


def _write_mode_cache(path: Path, payload: Mapping[str, object]) -> None:
    atomic_json(path, payload)


def compute_energy_phase(
    output_dir: Path,
    *,
    resume: bool,
    workers: int,
) -> dict[str, object]:
    """Phase B: verify fixed selection SHA and reconstruct unique saved-root modes."""

    output_dir = output_dir.resolve()
    selected_rows, selection_metadata = verify_frequency_selection(output_dir)
    unique_modes = build_unique_modes(selected_rows)
    unique_path = output_dir / "unique_selected_modes.csv"
    unique_payload = csv_bytes(UNIQUE_MODE_FIELDS, unique_modes)
    if unique_path.exists() and unique_path.read_bytes() != unique_payload:
        raise ValueError("existing unique mode inventory differs from deterministic inventory")
    atomic_bytes(unique_path, unique_payload)

    results_by_key: dict[tuple[str, int], dict[str, Any]] = {}
    missing: list[dict[str, object]] = []
    for mode in unique_modes:
        key = (str(mode["case_id"]), int(mode["sorted_k"]))
        cached = _load_mode_cache(_mode_cache_path(output_dir, mode), mode) if resume else None
        if cached is None:
            missing.append(mode)
        else:
            results_by_key[key] = cached
    cached_before_compute = len(results_by_key)
    missing_before_compute = len(missing)

    benchmark_path = output_dir / "energy_benchmark.json"
    benchmark_count = min(3, len(missing))
    benchmark_computed = 0
    if benchmark_count and not benchmark_path.exists():
        start = time.perf_counter()
        for mode in missing[:benchmark_count]:
            payload = reconstruct_unique_mode(mode)
            _write_mode_cache(_mode_cache_path(output_dir, mode), payload)
            results_by_key[(str(mode["case_id"]), int(mode["sorted_k"]))] = payload
        seconds = time.perf_counter() - start
        atomic_json(
            benchmark_path,
            {
                "schema_version": SCHEMA_VERSION,
                "scientific_scope": SCIENTIFIC_SCOPE,
                "unique_mode_count": len(unique_modes),
                "cached_before_benchmark": len(results_by_key) - benchmark_count,
                "benchmark_mode_count": benchmark_count,
                "benchmark_wall_seconds": seconds,
                "seconds_per_mode": seconds / benchmark_count,
                "estimated_full_seconds_at_scalar_rate": seconds / benchmark_count * len(unique_modes),
                "workers_requested_for_remaining": workers,
            },
        )
        missing = missing[benchmark_count:]
        benchmark_computed = benchmark_count

    if missing:
        if workers <= 1:
            computed = map(reconstruct_unique_mode, missing)
        else:
            executor = ThreadPoolExecutor(max_workers=workers)
            computed = executor.map(reconstruct_unique_mode, missing)
        try:
            for mode, payload in zip(missing, computed):
                _write_mode_cache(_mode_cache_path(output_dir, mode), payload)
                results_by_key[(str(mode["case_id"]), int(mode["sorted_k"]))] = payload
        finally:
            if workers > 1:
                executor.shutdown(wait=True)

    rows801: list[dict[str, object]] = []
    rows1601: list[dict[str, object]] = []
    convergence: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    for mode in unique_modes:
        key = (str(mode["case_id"]), int(mode["sorted_k"]))
        payload = results_by_key.get(key)
        if payload is None:
            payload = _load_mode_cache(_mode_cache_path(output_dir, mode), mode)
        if payload is None:
            raise RuntimeError(f"missing mode result after Phase B: {key}")
        if payload.get("energy_801") is not None:
            rows801.append(dict(payload["energy_801"]))
            rows1601.append(dict(payload["energy_1601"]))
            convergence.append(dict(payload["convergence"]))
        if payload.get("status") != "complete":
            failures.append(
                {
                    "case_id": key[0],
                    "epsilon_0": mode["epsilon_0"],
                    "beta_deg": mode["beta_deg"],
                    "mu": mode["mu"],
                    "eta": mode["eta"],
                    "s_max": mode["s_max"],
                    "coverage_status": "energy_quality_failed",
                    "coverage_reason": payload.get("error", "unknown_energy_failure"),
                    "eb_root_count": "",
                    "timo_root_count": "",
                    "eb_confirmed_first10": "",
                    "timo_confirmed_first10": "",
                    "delta_count_first10": "",
                    "source_frequency_file": mode["source_frequency_file"],
                    "source_frequency_sha256": mode["source_frequency_sha256"],
                }
            )

    atomic_csv(output_dir / "mode_energy_801.csv", MODE_ENERGY_FIELDS, rows801)
    atomic_csv(output_dir / "mode_energy_1601.csv", MODE_ENERGY_FIELDS, rows1601)
    atomic_csv(output_dir / "mode_energy_convergence.csv", MODE_CONVERGENCE_FIELDS, convergence)
    existing = [
        row
        for row in read_csv(output_dir / "quality_and_incomplete_cases.csv")
        if row.get("coverage_status") != "energy_quality_failed"
    ]
    atomic_csv(
        output_dir / "quality_and_incomplete_cases.csv",
        COVERAGE_FIELDS,
        [*existing, *failures],
    )
    return {
        "scientific_scope": SCIENTIFIC_SCOPE,
        "selection_sha256_verified": selection_metadata[
            "frequency_selected_triplets_sha256"
        ],
        "selected_triplet_count": len(selected_rows),
        "unique_mode_count": len(unique_modes),
        "cache_hit_count": cached_before_compute,
        "computed_mode_count": missing_before_compute,
        "benchmark_computed_mode_count": benchmark_computed,
        "quality_failed_mode_count": len(failures),
        "workers": workers,
    }


def assemble_triplet_results(output_dir: Path) -> list[dict[str, object]]:
    selected_rows, _metadata = verify_frequency_selection(output_dir)
    rows801 = read_csv(output_dir / "mode_energy_801.csv")
    rows1601 = read_csv(output_dir / "mode_energy_1601.csv")
    convergence_rows = read_csv(output_dir / "mode_energy_convergence.csv")
    energy801 = {(row["case_id"], int(row["sorted_k"])): row for row in rows801}
    energy1601 = {(row["case_id"], int(row["sorted_k"])): row for row in rows1601}
    convergence = {
        (row["case_id"], int(row["sorted_k"])): row for row in convergence_rows
    }
    results: list[dict[str, object]] = []
    for selected in selected_rows:
        case_id = selected["case_id"]
        k = int(selected["k"])
        keys = [(case_id, rank) for rank in (k - 1, k, k + 1)]
        missing = [key for key in keys if key not in energy801 or key not in energy1601]
        quality_failures = [
            key
            for key in keys
            if key in convergence and not true_text(convergence[key]["quality_pass"])
        ]
        quality_pass = not missing and not quality_failures
        failure_reasons: list[str] = []
        if missing:
            failure_reasons.append(
                "missing_energy_modes:" + ";".join(f"{case}:k{rank}" for case, rank in missing)
            )
        if quality_failures:
            failure_reasons.append(
                "quality_failed_modes:"
                + ";".join(f"{case}:k{rank}" for case, rank in quality_failures)
            )
        row: dict[str, object] = {field: selected[field] for field in FREQUENCY_FIELDS}
        if quality_pass:
            coarse = [float(energy801[key]["chi_a"]) for key in keys]
            fine = [float(energy1601[key]["chi_a"]) for key in keys]
            fine_bs = [float(energy1601[key]["chi_bs"]) for key in keys]
            coarse_metrics = energy_local_characteristics(*coarse)
            fine_metrics = energy_local_characteristics(*fine)
            row.update(
                {
                    "chi_a_k_minus_1_801": coarse[0],
                    "chi_a_k_801": coarse[1],
                    "chi_a_k_plus_1_801": coarse[2],
                    "chi_a_k_minus_1": fine[0],
                    "chi_a_k": fine[1],
                    "chi_a_k_plus_1": fine[2],
                    "chi_bs_k_minus_1": fine_bs[0],
                    "chi_bs_k": fine_bs[1],
                    "chi_bs_k_plus_1": fine_bs[2],
                    "D_a_801": coarse_metrics["D_a"],
                    "D_a": fine_metrics["D_a"],
                    "axial_above_neighbor_mean": bool_text(
                        fine_metrics["axial_above_neighbor_mean"]
                    ),
                    "axial_strict_local_max_801": bool_text(
                        coarse_metrics["axial_strict_local_max"]
                    ),
                    "axial_strict_local_max": bool_text(
                        fine_metrics["axial_strict_local_max"]
                    ),
                    "grid_classification_changed": bool_text(
                        coarse_metrics["axial_strict_local_max"]
                        != fine_metrics["axial_strict_local_max"]
                    ),
                }
            )
        else:
            row.update(
                {
                    "chi_a_k_minus_1_801": "",
                    "chi_a_k_801": "",
                    "chi_a_k_plus_1_801": "",
                    "chi_a_k_minus_1": "",
                    "chi_a_k": "",
                    "chi_a_k_plus_1": "",
                    "chi_bs_k_minus_1": "",
                    "chi_bs_k": "",
                    "chi_bs_k_plus_1": "",
                    "D_a_801": "",
                    "D_a": "",
                    "axial_above_neighbor_mean": "false",
                    "axial_strict_local_max_801": "false",
                    "axial_strict_local_max": "false",
                    "grid_classification_changed": "false",
                }
            )
        cohorts = cohort_flags(selected, energy_quality_pass=quality_pass)
        row.update(
            {
                "energy_quality_pass": bool_text(quality_pass),
                "energy_failure_reason": ";".join(failure_reasons),
                **{name: bool_text(value) for name, value in cohorts.items()},
            }
        )
        results.append(row)
    return results


def _quantiles(values: Sequence[float]) -> tuple[object, object, object, object, object]:
    if not values:
        return "", "", "", "", ""
    array = np.asarray(values, dtype=float)
    q1, median, q3 = np.quantile(array, [0.25, 0.5, 0.75], method="linear")
    return float(np.min(array)), float(q1), float(median), float(q3), float(np.max(array))


def summarize_subset(
    rows: Sequence[Mapping[str, object]],
    *,
    cohort: str,
    group_field: str,
    group_value: str,
) -> dict[str, object]:
    valid = [row for row in rows if true_text(row["energy_quality_pass"])]
    d_f = [float(row["D_f"]) for row in valid]
    d_a = [float(row["D_a"]) for row in valid]
    d_f_stats = _quantiles(d_f)
    d_a_stats = _quantiles(d_a)
    support = max(valid, key=lambda row: float(row["D_a"])) if valid else None
    contradict = min(valid, key=lambda row: float(row["D_a"])) if valid else None
    n_loc = len(valid)
    n_d_a = sum(float(row["D_a"]) > 0.0 for row in valid)
    n_ax = sum(true_text(row["axial_strict_local_max"]) for row in valid)
    n_threshold_any = sum(true_text(row["threshold_any"]) for row in valid)
    n_threshold_any_ax = sum(
        true_text(row["threshold_any"]) and true_text(row["axial_strict_local_max"])
        for row in valid
    )
    n_threshold_both = sum(true_text(row["threshold_both"]) for row in valid)
    n_threshold_both_ax = sum(
        true_text(row["threshold_both"]) and true_text(row["axial_strict_local_max"])
        for row in valid
    )
    return {
        "cohort": cohort,
        "group_field": group_field,
        "group_value": group_value,
        "n_loc": n_loc,
        "n_geometries": len({row["case_id"] for row in valid}),
        "n_mu_eta": len({(float(row["mu"]), float(row["eta"])) for row in valid}),
        "n_D_a_positive": n_d_a,
        "n_axial_strict_local_max": n_ax,
        "ratio_D_a_positive_percent": 100.0 * n_d_a / n_loc if n_loc else "",
        "ratio_axial_strict_local_max_percent": 100.0 * n_ax / n_loc if n_loc else "",
        "n_threshold_any": n_threshold_any,
        "n_threshold_any_axial_strict_local_max": n_threshold_any_ax,
        "n_threshold_both": n_threshold_both,
        "n_threshold_both_axial_strict_local_max": n_threshold_both_ax,
        "D_f_min": d_f_stats[0],
        "D_f_q1": d_f_stats[1],
        "D_f_median": d_f_stats[2],
        "D_f_q3": d_f_stats[3],
        "D_f_max": d_f_stats[4],
        "D_a_min": d_a_stats[0],
        "D_a_q1": d_a_stats[1],
        "D_a_median": d_a_stats[2],
        "D_a_q3": d_a_stats[3],
        "D_a_max": d_a_stats[4],
        "strongest_support_case_id": support["case_id"] if support else "",
        "strongest_support_k": support["k"] if support else "",
        "strongest_support_D_a": support["D_a"] if support else "",
        "strongest_contradiction_case_id": contradict["case_id"] if contradict else "",
        "strongest_contradiction_k": contradict["k"] if contradict else "",
        "strongest_contradiction_D_a": contradict["D_a"] if contradict else "",
    }


def _group_rows(
    rows: Sequence[Mapping[str, object]],
    field: str,
) -> list[tuple[str, list[Mapping[str, object]]]]:
    groups: dict[str, list[Mapping[str, object]]] = {}
    for row in rows:
        if field == "mu_eta":
            value = f"mu={float(row['mu']):g},eta={float(row['eta']):g}"
        elif field == "k":
            value = str(int(row["k"]))
        else:
            value = f"{float(row[field]):g}"
        groups.setdefault(value, []).append(row)
    return sorted(
        groups.items(),
        key=lambda item: tuple(
            float(part.split("=")[-1])
            for part in item[0].split(",")
        ),
    )


def write_article_summary_table(
    output_dir: Path,
    primary: Sequence[Mapping[str, object]],
) -> None:
    rows: list[dict[str, object]] = []
    for epsilon in BASE_EPSILON_LEVELS:
        subset = [row for row in primary if same_float(row["epsilon_0"], epsilon)]
        n_loc = len(subset)
        n_ax = sum(true_text(row["axial_strict_local_max"]) for row in subset)
        threshold_any = sum(true_text(row["threshold_any"]) for row in subset)
        threshold_any_ax = sum(
            true_text(row["threshold_any"]) and true_text(row["axial_strict_local_max"])
            for row in subset
        )
        rows.append(
            {
                "epsilon_0": f"{epsilon:.3f}",
                "n_loc": n_loc,
                "n_ax": n_ax,
                "n_ax/n_loc, %": f"{100.0*n_ax/n_loc:.1f}" if n_loc else "—",
                "threshold_any": threshold_any,
                "threshold_any и локальный максимум chi_a": threshold_any_ax,
            }
        )
    n_loc = len(primary)
    n_ax = sum(true_text(row["axial_strict_local_max"]) for row in primary)
    threshold_any = sum(true_text(row["threshold_any"]) for row in primary)
    threshold_any_ax = sum(
        true_text(row["threshold_any"]) and true_text(row["axial_strict_local_max"])
        for row in primary
    )
    rows.append(
        {
            "epsilon_0": "Всего",
            "n_loc": n_loc,
            "n_ax": n_ax,
            "n_ax/n_loc, %": f"{100.0*n_ax/n_loc:.1f}" if n_loc else "—",
            "threshold_any": threshold_any,
            "threshold_any и локальный максимум chi_a": threshold_any_ax,
        }
    )
    fields = tuple(rows[0])
    atomic_csv(
        output_dir / "article_energy_local_minima_summary_table.tsv",
        fields,
        rows,
        delimiter="\t",
    )
    tex = [
        "\\begin{tabular}{lrrrrr}",
        "\\hline",
        "$\\epsilon_0$ & $n_{loc}$ & $n_{ax}$ & $n_{ax}/n_{loc}$, \\% & "
        "$n_{thr}$ & $n_{thr,ax}$ \\\\",
        "\\hline",
        *(
            f"{row['epsilon_0']} & {row['n_loc']} & {row['n_ax']} & "
            f"{row['n_ax/n_loc, %']} & {row['threshold_any']} & "
            f"{row['threshold_any и локальный максимум chi_a']} \\\\"
            for row in rows
        ),
        "\\hline",
        "\\end{tabular}",
    ]
    atomic_text(
        output_dir / "article_energy_local_minima_summary_table.tex",
        "\n".join(tex) + "\n",
    )

    overall = rows[-1]
    compact = [
        {
            "совокупность": "Подтверждающая совокупность",
            "n_loc": overall["n_loc"],
            "n_ax": overall["n_ax"],
            "n_ax/n_loc, %": overall["n_ax/n_loc, %"],
            "threshold_any": overall["threshold_any"],
            "threshold_any и локальный максимум chi_a": overall[
                "threshold_any и локальный максимум chi_a"
            ],
        }
    ]
    compact_fields = tuple(compact[0])
    atomic_csv(
        output_dir / "article_energy_local_minima_summary_overall.tsv",
        compact_fields,
        compact,
        delimiter="\t",
    )
    compact_tex = [
        "\\begin{tabular}{lrrrrr}",
        "\\hline",
        "Совокупность & $n_{loc}$ & $n_{ax}$ & $n_{ax}/n_{loc}$, \\% & "
        "$n_{thr}$ & $n_{thr,ax}$ \\\\",
        "\\hline",
        f"Подтверждающая & {overall['n_loc']} & {overall['n_ax']} & "
        f"{overall['n_ax/n_loc, %']} & {overall['threshold_any']} & "
        f"{overall['threshold_any и локальный максимум chi_a']} \\\\",
        "\\hline",
        "\\end{tabular}",
    ]
    atomic_text(
        output_dir / "article_energy_local_minima_summary_overall.tex",
        "\n".join(compact_tex) + "\n",
    )


def write_article_text(
    output_dir: Path,
    primary_summary: Mapping[str, object],
    holdout_summary: Mapping[str, object],
) -> None:
    n_loc = int(primary_summary["n_loc"])
    n_ax = int(primary_summary["n_axial_strict_local_max"])
    h_loc = int(holdout_summary["n_loc"])
    h_ax = int(holdout_summary["n_axial_strict_local_max"])
    if n_loc and 2 * n_ax > n_loc and h_loc and 2 * h_ax > h_loc:
        core = (
            "На принятой конечной сетке большинство строгих локальных минимумов "
            "расхождения частот сопровождалось локальным повышением относительного "
            "вклада продольной деформации по сравнению с соседними формами."
        )
    elif n_ax > 0:
        core = (
            "Повышенный относительный вклад продольной деформации сопровождает часть "
            "наблюдаемых локальных минимумов расхождения, но не является их "
            "универсальным признаком."
        )
    else:
        core = (
            "На принятой конечной сетке устойчивая связь между строгими локальными "
            "минимумами расхождения частот и локальным повышением относительного "
            "вклада продольной деформации не обнаружена; два ранее рассмотренных "
            "случая следует сохранять только как иллюстративные примеры."
        )
    text = (
        core
        + "\n\n"
        + f"В подтверждающей совокупности локальный максимум наблюдался в {n_ax} из "
        f"{n_loc} случаев; после исключения всей семьи $\\mu=0.5$, $\\eta=0.5$ — "
        f"в {h_ax} из {h_loc} случаев.\n\n"
        + "Полученный результат относится к принятой конечной сетке геометрических "
        "параметров и не устанавливает универсальной зависимости для всех геометрий "
        "и номеров спектра.\n"
    )
    atomic_text(output_dir / "article_energy_local_minima_text_ru.md", text)


def write_hash_inventory(output_dir: Path, input_hashes: Mapping[str, str]) -> Path:
    inventory_path = output_dir / "input_output_sha256.csv"
    rows = [
        {"kind": "input", "path": path, "sha256": value}
        for path, value in sorted(input_hashes.items())
    ]
    rows.extend(
        {
            "kind": "output",
            "path": path.resolve().as_posix(),
            "sha256": sha256_file(path),
        }
        for path in sorted(output_dir.iterdir())
        if path.is_file() and path != inventory_path
    )
    return atomic_csv(inventory_path, ("kind", "path", "sha256"), rows)


def summarize_phase(output_dir: Path) -> dict[str, object]:
    output_dir = output_dir.resolve()
    selected_rows, selection_metadata = verify_frequency_selection(output_dir)
    for name in ("mode_energy_801.csv", "mode_energy_1601.csv", "mode_energy_convergence.csv"):
        if not (output_dir / name).is_file():
            raise FileNotFoundError(f"summarize-only requires {name}")
    triplets = assemble_triplet_results(output_dir)
    atomic_csv(output_dir / "energy_triplet_results.csv", TRIPLET_FIELDS, triplets)
    primary = [row for row in triplets if true_text(row["primary_confirmatory"])]
    holdout = [row for row in triplets if true_text(row["strict_family_holdout"])]
    beta0 = [row for row in triplets if true_text(row["beta0_control"])]
    extended = [row for row in triplets if true_text(row["extended_range"])]
    clusters = [row for row in triplets if true_text(row["cluster_affected"])]
    discovery = [row for row in triplets if true_text(row["discovery_cohort"])]
    atomic_csv(output_dir / "primary_confirmatory_results.csv", TRIPLET_FIELDS, primary)
    atomic_csv(output_dir / "strict_family_holdout_results.csv", TRIPLET_FIELDS, holdout)
    atomic_csv(output_dir / "beta0_control_results.csv", TRIPLET_FIELDS, beta0)
    atomic_csv(output_dir / "extended_range_results.csv", TRIPLET_FIELDS, extended)
    atomic_csv(output_dir / "cluster_affected_results.csv", TRIPLET_FIELDS, clusters)

    contradicting: list[dict[str, object]] = []
    for row in triplets:
        reasons: list[str] = []
        if not true_text(row["energy_quality_pass"]):
            reasons.append("energy_quality_failed")
        elif float(row["D_a"]) <= 0.0:
            reasons.append("D_a_not_positive")
        if not true_text(row["axial_strict_local_max"]):
            reasons.append("axial_strict_local_max_false")
        if true_text(row["grid_classification_changed"]):
            reasons.append("801_1601_classification_changed")
        if true_text(row["cluster_affected"]):
            reasons.append("cluster_affected")
        if reasons:
            item = dict(row)
            item["contradiction_reasons"] = ";".join(reasons)
            contradicting.append(item)
    atomic_csv(output_dir / "contradicting_cases.csv", CONTRADICTING_FIELDS, contradicting)

    cohorts = {
        "primary_confirmatory": primary,
        "strict_family_holdout": holdout,
        "beta0_control": beta0,
        "extended_range": extended,
        "cluster_affected": clusters,
        "discovery_cases": discovery,
    }
    overall = [
        summarize_subset(rows, cohort=name, group_field="all", group_value="all")
        for name, rows in cohorts.items()
    ]
    atomic_csv(output_dir / "summary_overall.csv", SUMMARY_FIELDS, overall)
    breakdown_specs = (
        ("epsilon_0", "summary_by_epsilon.csv"),
        ("k", "summary_by_k.csv"),
        ("beta_deg", "summary_by_beta.csv"),
        ("mu", "summary_by_mu.csv"),
        ("eta", "summary_by_eta.csv"),
        ("mu_eta", "summary_by_mu_eta.csv"),
    )
    for field, filename in breakdown_specs:
        summary_rows: list[dict[str, object]] = []
        for cohort in ("primary_confirmatory", "strict_family_holdout"):
            for value, rows in _group_rows(cohorts[cohort], field):
                summary_rows.append(
                    summarize_subset(
                        rows,
                        cohort=cohort,
                        group_field=field,
                        group_value=value,
                    )
                )
        atomic_csv(output_dir / filename, SUMMARY_FIELDS, summary_rows)

    primary_summary = next(row for row in overall if row["cohort"] == "primary_confirmatory")
    holdout_summary = next(row for row in overall if row["cohort"] == "strict_family_holdout")
    write_article_summary_table(output_dir, primary)
    write_article_text(output_dir, primary_summary, holdout_summary)

    convergence = read_csv(output_dir / "mode_energy_convergence.csv")
    quality_incomplete = read_csv(output_dir / "quality_and_incomplete_cases.csv")
    discovery_inventory = read_csv(output_dir / "discovery_case_ids.csv")
    coverage_by_id: dict[str, dict[str, str]] = {}
    for coverage_row in quality_incomplete:
        if coverage_row["coverage_status"] != "energy_quality_failed":
            coverage_by_id.setdefault(coverage_row["case_id"], coverage_row)
    discovery_results: list[dict[str, object]] = []
    for item in discovery_inventory:
        case_id = item["case_id"]
        coverage_row = coverage_by_id.get(case_id, {})
        case_triplets = [row for row in triplets if row["case_id"] == case_id]
        discovery_results.append(
            {
                **item,
                "first10_coverage_status": coverage_row.get(
                    "coverage_status", "missing_main_grid_coverage_row"
                ),
                "first10_coverage_reason": coverage_row.get("coverage_reason", ""),
                "frequency_selected_triplet_count": len(case_triplets),
                "energy_quality_pass_triplet_count": sum(
                    true_text(row["energy_quality_pass"]) for row in case_triplets
                ),
            }
        )
    atomic_csv(
        output_dir / "discovery_cases_results.csv",
        DISCOVERY_RESULT_FIELDS,
        discovery_results,
    )

    convergence_abs_stats = _quantiles(
        [float(row["delta_chi_a_abs"]) for row in convergence]
    )
    convergence_rel_stats = _quantiles(
        [float(row["delta_chi_a_rel"]) for row in convergence]
    )
    quality_check_failure_counts = {
        field: sum(not true_text(row[field]) for row in convergence)
        for field in (
            "integration_pass",
            "svd_pass",
            "null_residual_pass",
            "clamp_pass",
            "joint_kinematic_pass",
            "joint_force_pass",
            "sign_flip_pass",
            "quality_pass",
        )
    }
    prior_after = prior_fingerprints()
    prior_before = selection_metadata["prior_result_fingerprints"]
    prior_unchanged = prior_after == prior_before
    if not prior_unchanged:
        raise RuntimeError("immutable v1/v2/v3 energy result fingerprint changed")
    unique_mode_count = len(read_csv(output_dir / "unique_selected_modes.csv"))
    fallback_count = sum(
        str(row["null_vector_source"]).endswith("fallback") for row in convergence
    )
    operation_counts = {
        "root_solver_calls": 0,
        "strict_solver_calls": 0,
        "family_detector_calls": 0,
        "local_repair_calls": 0,
        "MAC_calls": 0,
        "shape_plot_calls": 0,
        "FEM_calls": 0,
        "anisotropic_calls": 0,
        "matrix_evaluator_calls": 2 * unique_mode_count,
        "null_vector_reconstructions": unique_mode_count + fallback_count,
        "energy_integral_evaluations": 3 * unique_mode_count,
    }
    metadata = {
        "schema_version": SCHEMA_VERSION,
        "scientific_scope": SCIENTIFIC_SCOPE,
        "selection_sha256_verified": selection_metadata[
            "frequency_selected_triplets_sha256"
        ],
        "main_grid_case_count": selection_metadata["main_grid_case_count"],
        "eligible_first10_case_count": selection_metadata["eligible_first10_case_count"],
        "quality_or_incomplete_case_count": selection_metadata[
            "quality_or_incomplete_case_count"
        ],
        "frequency_selected_triplet_count": len(selected_rows),
        "unique_mode_count": unique_mode_count,
        "energy_complete_mode_count": len(read_csv(output_dir / "mode_energy_1601.csv")),
        "energy_quality_failed_mode_count": sum(
            not true_text(row["quality_pass"]) for row in convergence
        ),
        "primary_confirmatory_triplet_count": len(primary),
        "strict_family_holdout_triplet_count": len(holdout),
        "beta0_control_triplet_count": len(beta0),
        "extended_range_triplet_count": len(extended),
        "cluster_affected_triplet_count": len(clusters),
        "discovery_triplet_count": len(discovery),
        "discovery_case_inventory_count": len(discovery_inventory),
        "contradicting_case_count": len(contradicting),
        "energy_quality_failed_triplet_count": sum(
            not true_text(row["energy_quality_pass"]) for row in triplets
        ),
        "grid_classification_changed_triplet_count": sum(
            true_text(row["grid_classification_changed"]) for row in triplets
        ),
        "quality_and_incomplete_row_count": len(quality_incomplete),
        "integration_grids": list(INTEGRATION_GRIDS),
        "integration_tolerance": PILOT.INTEGRATION_TOL,
        "energy_sum_tolerance": PILOT.ENERGY_SUM_TOL,
        "mode_residual_tolerance": PILOT.MODE_RESIDUAL_TOL,
        "chi_a_801_1601_absolute_change": {
            "min": convergence_abs_stats[0],
            "q1": convergence_abs_stats[1],
            "median": convergence_abs_stats[2],
            "q3": convergence_abs_stats[3],
            "max": convergence_abs_stats[4],
        },
        "chi_a_801_1601_relative_change": {
            "min": convergence_rel_stats[0],
            "q1": convergence_rel_stats[1],
            "median": convergence_rel_stats[2],
            "q3": convergence_rel_stats[3],
            "max": convergence_rel_stats[4],
        },
        "quality_check_failure_counts": quality_check_failure_counts,
        "prior_result_fingerprints_before": prior_before,
        "prior_result_fingerprints_after": prior_after,
        "prior_results_unchanged": prior_unchanged,
        "operation_counts": operation_counts,
        "no_new_roots": True,
        "same_index_sorted_spectra": True,
    }
    atomic_json(output_dir / "run_metadata.json", metadata)

    primary_ratio = (
        100.0
        * int(primary_summary["n_axial_strict_local_max"])
        / int(primary_summary["n_loc"])
        if int(primary_summary["n_loc"])
        else float("nan")
    )
    holdout_ratio = (
        100.0
        * int(holdout_summary["n_axial_strict_local_max"])
        / int(holdout_summary["n_loc"])
        if int(holdout_summary["n_loc"])
        else float("nan")
    )
    report = [
        "# Finite-grid local-minimum energy check",
        "",
        f"Scientific scope: `{SCIENTIFIC_SCOPE}`.",
        "",
        "Triplets were fixed from saved frequency data before energy reconstruction. "
        "The comparison uses identical positions in the ordered EB and Timoshenko spectra; "
        "physical descendants are not tracked.",
        "",
        "## Coverage",
        "",
        f"- Main-grid geometries: {selection_metadata['main_grid_case_count']}",
        f"- Complete, quality-confirmed first-ten spectra: {selection_metadata['eligible_first10_case_count']}",
        f"- Quality/incomplete geometries: {selection_metadata['quality_or_incomplete_case_count']}",
        f"- Frequency-selected triplets: {len(selected_rows)}",
        f"- Unique reconstructed forms: {metadata['unique_mode_count']}",
        "- Eligible first-ten spectra by epsilon_0: "
        + ", ".join(
            f"{epsilon}={counts['eligible_first10']}/{counts['total']}"
            for epsilon, counts in selection_metadata["coverage_by_epsilon"].items()
        ),
        "",
        "## Confirmatory cohorts",
        "",
        f"- PRIMARY CONFIRMATORY: {primary_summary['n_axial_strict_local_max']}/{primary_summary['n_loc']} ({primary_ratio:.3f}%) strict local maxima of chi_a.",
        f"- STRICT FAMILY HOLDOUT: {holdout_summary['n_axial_strict_local_max']}/{holdout_summary['n_loc']} ({holdout_ratio:.3f}%) strict local maxima of chi_a.",
        f"- beta=0 control triplets: {len(beta0)}",
        f"- extended s_max>0.1 triplets: {len(extended)}",
        f"- cluster-affected triplets: {len(clusters)}",
        f"- prior discovery cases in inventory: {len(discovery_inventory)}; "
        f"frequency-selected full-first-ten discovery triplets: {len(discovery)}",
        f"- rows requiring contradictory/caution reporting: {len(contradicting)}",
        f"- PRIMARY D_a>0: {primary_summary['n_D_a_positive']}/{primary_summary['n_loc']} "
        f"({optional_float_text(primary_summary['ratio_D_a_positive_percent'], '.3f')}%).",
        f"- strongest PRIMARY support: {primary_summary['strongest_support_case_id']} "
        f"k={primary_summary['strongest_support_k'] or 'n/a'}, "
        f"D_a={optional_float_text(primary_summary['strongest_support_D_a'], '.6e')}.",
        f"- strongest PRIMARY contradiction: {primary_summary['strongest_contradiction_case_id']} "
        f"k={primary_summary['strongest_contradiction_k'] or 'n/a'}, "
        f"D_a={optional_float_text(primary_summary['strongest_contradiction_D_a'], '.6e')}.",
        "",
        "## Numerical checks",
        "",
        f"Absolute 801/1601 chi_a change (q1/median/q3/max): "
        f"{float(convergence_abs_stats[1]):.6e} / {float(convergence_abs_stats[2]):.6e} / "
        f"{float(convergence_abs_stats[3]):.6e} / {float(convergence_abs_stats[4]):.6e}.",
        f"Relative 801/1601 chi_a change (q1/median/q3/max): "
        f"{float(convergence_rel_stats[1]):.6e} / {float(convergence_rel_stats[2]):.6e} / "
        f"{float(convergence_rel_stats[3]):.6e} / {float(convergence_rel_stats[4]):.6e}.",
        f"Integration failures: {quality_check_failure_counts['integration_pass']}; "
        f"grid-classification changes: {metadata['grid_classification_changed_triplet_count']}.",
        f"Energy/residual quality failures among reconstructed modes: {metadata['energy_quality_failed_mode_count']}.",
        f"Residual failure counts (kinematic/force): "
        f"{quality_check_failure_counts['joint_kinematic_pass']}/"
        f"{quality_check_failure_counts['joint_force_pass']}.",
        "",
        "No new roots, strict solve, branch tracking, MAC, FEM, or anisotropic calculation was performed.",
    ]
    atomic_text(output_dir / "report.md", "\n".join(report) + "\n")
    write_hash_inventory(output_dir, selection_metadata["input_file_hashes"])
    return metadata


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Two-phase finite-grid isotropic-circular energy check for strict local "
            "minima of same-index EB/Timoshenko frequency mismatch."
        ),
    )
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--select-only", action="store_true")
    modes.add_argument("--compute-energy", action="store_true")
    modes.add_argument("--summarize-only", action="store_true")
    modes.add_argument("--run-all", action="store_true")
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--compact-dir", type=Path, default=DEFAULT_COMPACT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--smoke-cases", type=int, default=4)
    parser.add_argument("--workers", type=int, default=1)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.workers < 1:
        raise ValueError("--workers must be positive")
    if args.smoke_cases < 1:
        raise ValueError("--smoke-cases must be positive")
    output_dir = args.output_dir
    if args.smoke and output_dir.resolve() == DEFAULT_OUTPUT_DIR.resolve():
        output_dir = REPO_ROOT / "results" / "_smoke" / "energy_local_minima_grid_v1"
    result: dict[str, object] = {
        "scientific_scope": SCIENTIFIC_SCOPE,
        "output_dir": output_dir.resolve().as_posix(),
    }
    if args.select_only or args.run_all:
        result["selection"] = select_frequency_triplets(
            args.source_dir,
            args.compact_dir,
            output_dir,
            resume=args.resume,
            smoke=args.smoke,
            smoke_cases=args.smoke_cases,
        )
    if args.compute_energy or args.run_all:
        result["energy"] = compute_energy_phase(
            output_dir,
            resume=args.resume or args.run_all,
            workers=args.workers,
        )
    if args.summarize_only or args.run_all:
        result["summary"] = summarize_phase(output_dir)
    console_result: dict[str, object] = {
        "scientific_scope": result["scientific_scope"],
        "output_dir": result["output_dir"],
    }
    if "selection" in result:
        selection = result["selection"]
        assert isinstance(selection, Mapping)
        console_result["selection"] = {
            key: selection[key]
            for key in (
                "main_grid_case_count",
                "eligible_first10_case_count",
                "quality_or_incomplete_case_count",
                "frequency_selected_triplet_count",
                "frequency_selected_geometry_count",
                "frequency_selected_triplets_sha256",
            )
        }
    if "energy" in result:
        console_result["energy"] = result["energy"]
    if "summary" in result:
        summary = result["summary"]
        assert isinstance(summary, Mapping)
        console_result["summary"] = {
            key: summary[key]
            for key in (
                "frequency_selected_triplet_count",
                "unique_mode_count",
                "energy_quality_failed_mode_count",
                "primary_confirmatory_triplet_count",
                "strict_family_holdout_triplet_count",
                "beta0_control_triplet_count",
                "extended_range_triplet_count",
                "cluster_affected_triplet_count",
                "contradicting_case_count",
                "prior_results_unchanged",
            )
        }
    print(json.dumps(console_result, ensure_ascii=False, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
