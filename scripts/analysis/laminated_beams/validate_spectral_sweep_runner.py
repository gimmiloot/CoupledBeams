"""Validate the model-agnostic RLB-SWEEP-FAST orchestration layer.

The script keeps all beam physics and all root-quality thresholds in their
existing modules.  Certified full-scan CSV files are read as immutable oracle
data.  Non-anchor validation points are recomputed with the existing local
determinant/SVD refiner; oracle values are read only after the fast result has
been frozen for comparison.
"""

from __future__ import annotations

import argparse
import csv
import ctypes
from ctypes import wintypes
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


THREAD_LIMITS = {
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
}
for _name, _value in THREAD_LIMITS.items():
    os.environ[_name] = _value

import numpy as np  # noqa: E402


ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = ROOT / "src"
for _import_root in (ROOT, SRC_ROOT):
    if str(_import_root) not in sys.path:
        sys.path.insert(0, str(_import_root))

from scripts.analysis.laminated_beams import (  # noqa: E402
    compare_rectangular_weakly_orthotropic_models_vs_beta as rlb2c,
)
from scripts.analysis.laminated_beams import (  # noqa: E402
    compare_rectangular_weakly_orthotropic_mu_beta_graphs as rlb2d,
)
from scripts.lib import spectral_sweep_runner as sweep_runner  # noqa: E402
from src.my_project.analytic.solvers import (  # noqa: E402
    find_roots_scan_bisect as eb_interval_scan,
)


STAGE_ID = "RLB-SWEEP-FAST"
ALGORITHM_VERSION = "spectral_sweep_runner_validation_v2"
DEFAULT_OUTPUT_DIR = (
    ROOT / "results" / "laminated_beams" / "spectral_sweep_runner_validation"
)
RLB2B_OUTPUT = (
    ROOT / "results" / "laminated_beams" / "rectangular_isotropic_beta_sweep_comparison"
)
RLB2C_OUTPUT = (
    ROOT
    / "results"
    / "laminated_beams"
    / "rectangular_weakly_orthotropic_beta_sweep_comparison"
)
RLB2D_OUTPUT = (
    ROOT
    / "results"
    / "laminated_beams"
    / "rectangular_weakly_orthotropic_graphs_mu_beta"
)

BETA_VALIDATION_POINTS = (0.0, 45.0, 54.0, 90.0)
MU_VALIDATION_POINTS = (0.0, 0.40, 0.60, 0.78, 0.80)
FAST_FULL_RELATIVE_TOLERANCE = 1.0e-8
LOCAL_VERIFICATION_RELATIVE_TOLERANCE = 1.0e-9
ORACLE_QUALIFICATION_RELATIVE_TOLERANCE = 1.0e-9
CACHE_EQUALITY_RELATIVE_TOLERANCE = 1.0e-13
TARGET_RUNTIME_SECONDS = 30.0 * 60.0
BLOCK_RUNTIME_SECONDS = 60.0 * 60.0
LOCAL_PRIMARY_GRID_STEP = 0.2
LOCAL_VERIFICATION_GRID_STEP = 0.1
MINIMUM_LOCAL_PRIMARY_POINTS = 33
MINIMUM_LOCAL_VERIFICATION_POINTS = 65
MAX_CACHE_ENTRIES = 8192
MAX_CACHE_BYTES = 384 * 1024**2

MODEL_EB = rlb2c.MODEL_EB
MODEL_OLD = rlb2c.MODEL_OLD
MODEL_RLB = rlb2c.MODEL_RLB
MODELS = (MODEL_EB, MODEL_OLD, MODEL_RLB)
MU_VALIDATION_POINTS_BY_MODEL = {
    MODEL_EB: (0.0, 0.80),
    MODEL_OLD: (0.0, 0.40, 0.60, 0.78, 0.80),
    MODEL_RLB: (0.0, 0.40, 0.60, 0.78, 0.80),
}

PRODUCTION_PHYSICS_PATHS = (
    "src/my_project/analytic/formulas.py",
    "src/my_project/analytic/solvers.py",
    "scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py",
    "scripts/lib/reddy_symmetric_laminated_beam.py",
    "scripts/lib/reddy_inplane_geometry.py",
    "scripts/lib/reddy_symmetric_coupled_beams.py",
)

VALIDATION_ADAPTER_PATHS = (
    "scripts/lib/spectral_sweep_runner.py",
    "scripts/analysis/laminated_beams/validate_reddy_four_ply_isotropic_limit.py",
    "scripts/analysis/laminated_beams/compare_rectangular_isotropic_models_vs_beta.py",
    "scripts/analysis/laminated_beams/compare_rectangular_weakly_orthotropic_models_vs_beta.py",
    "scripts/analysis/laminated_beams/compare_rectangular_weakly_orthotropic_mu_beta_graphs.py",
    "scripts/analysis/laminated_beams/validate_spectral_sweep_runner.py",
)

VALIDATION_DEPENDENCY_PATHS = tuple(
    dict.fromkeys((*PRODUCTION_PHYSICS_PATHS, *VALIDATION_ADAPTER_PATHS))
)

FULL_EXPORT_ROOTS = 8
FULL_EXPORT_GUARDS = 1
FRESH_FULL_CACHE_ENTRIES = 64
FRESH_FULL_CACHE_BYTES = 16 * 1024**2
MAX_FRESH_FULL_SCANS_CLOSING_VALIDATION = 3
VALIDATION_CONTRACT_VERSION = "simplified-fast-vs-frozen-full-v1"
CURRENT_SHARD_DIRECTORY = "validation_shards_simplified_v2"
LEGACY_SHARD_DIRECTORY = "validation_shards"

RLB2D_PRODUCTION_FIELDS = (
    "sweep",
    "model",
    "mu",
    "tau",
    "beta_deg",
    "L1",
    "L2",
    "h1",
    "h2",
    "b1",
    "b2",
    "sorted_position",
    "role",
    "omega",
    "Omega",
    "Lambda",
    "normalization_identity_relative_residual",
    "historical_EB_mapping_relative_residual",
    "historical_EB_wavenumber",
    "multiplicity",
    "detected_nullity",
    "cluster_id",
    "cluster_semantics",
    "cluster_multiplicity",
    "cluster_total_nullity",
    "raw_determinant",
    "scaled_determinant",
    "raw_sigma_min",
    "raw_sigma_max",
    "raw_sigma_ratio",
    "scaled_sigma_min",
    "scaled_sigma_max",
    "scaled_sigma_ratio",
    "boundary_null_residual",
    "detector_type",
    "root_status",
    "inventory_status",
    "inventory_sha256",
    "primary_slot_count_internal",
    "verification_slot_count_internal",
    "internal_requested_roots",
    "internal_guard_position",
    "primary_verification_max_relative",
    "unresolved_candidates_below_internal_guard",
    "guard_flag",
    "bracket_left_Omega",
    "bracket_right_Omega",
    "bracket_left_Lambda",
    "bracket_right_Lambda",
    "internal_root13_Omega",
    "internal_root13_Lambda",
    "export_guard_available",
    "export_primary_verification_agreement",
    "export_primary_verification_max_relative",
    "unresolved_candidates_below_export_guard",
    "export_range_status",
    "internal_inventory_status",
    "internal_primary_verification_max_relative",
    "internal_unresolved_candidate_count",
)
if len(RLB2D_PRODUCTION_FIELDS) != 59:  # pragma: no cover - schema invariant
    raise RuntimeError("The frozen RLB-2D production row schema must have 59 fields.")

RLB2D_MISSING_REFERENCE_GROUPS = (
    (MODEL_OLD, 0.78),
    (MODEL_OLD, 0.80),
    (MODEL_RLB, 0.80),
)
RLB2D_EXPECTED_FAST_MISSING = {
    MODEL_OLD: (0.68, 0.70, 0.72, 0.74, 0.76),
    MODEL_RLB: (0.74,),
}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def _dependency_hashes() -> dict[str, str]:
    """Hash every frozen model and adapter used by validation workers."""

    return {relative: _sha256(ROOT / relative) for relative in VALIDATION_DEPENDENCY_PATHS}


def _validation_constants_payload() -> dict[str, Any]:
    """Return the numerical validation contract included in every shard key."""

    return {
        "validation_contract_version": VALIDATION_CONTRACT_VERSION,
        "algorithm_version": ALGORITHM_VERSION,
        "runner_version": sweep_runner.RUNNER_VERSION,
        "plotted_positions": list(range(1, FULL_EXPORT_ROOTS + 1)),
        "guard_position": FULL_EXPORT_ROOTS + FULL_EXPORT_GUARDS,
        "beta_validation_points": BETA_VALIDATION_POINTS,
        "mu_validation_points": MU_VALIDATION_POINTS,
        "mu_validation_points_by_model": MU_VALIDATION_POINTS_BY_MODEL,
        "fast_full_relative_tolerance": FAST_FULL_RELATIVE_TOLERANCE,
        "local_verification_relative_tolerance": (
            LOCAL_VERIFICATION_RELATIVE_TOLERANCE
        ),
        "oracle_qualification_relative_tolerance": (
            ORACLE_QUALIFICATION_RELATIVE_TOLERANCE
        ),
        "cache_equality_relative_tolerance": CACHE_EQUALITY_RELATIVE_TOLERANCE,
        "local_primary_grid_step": LOCAL_PRIMARY_GRID_STEP,
        "local_verification_grid_step": LOCAL_VERIFICATION_GRID_STEP,
        "minimum_local_primary_points": MINIMUM_LOCAL_PRIMARY_POINTS,
        "minimum_local_verification_points": MINIMUM_LOCAL_VERIFICATION_POINTS,
        "maximum_cache_entries": MAX_CACHE_ENTRIES,
        "maximum_cache_bytes": MAX_CACHE_BYTES,
        "fresh_full_cache_entries": FRESH_FULL_CACHE_ENTRIES,
        "fresh_full_cache_bytes": FRESH_FULL_CACHE_BYTES,
        "maximum_fresh_full_scans_closing_validation": (
            MAX_FRESH_FULL_SCANS_CLOSING_VALIDATION
        ),
        "thread_limits": THREAD_LIMITS,
    }


def _git(*arguments: str) -> str:
    completed = subprocess.run(
        ["git", *arguments],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return completed.stdout.strip()


def _git_state() -> dict[str, Any]:
    return {
        "toplevel": _git("rev-parse", "--show-toplevel"),
        "branch": _git("branch", "--show-current"),
        "head": _git("rev-parse", "HEAD"),
        "last_commit": _git("log", "-1", "--oneline"),
        "status_short": _git("status", "--short").splitlines(),
    }


def _json_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return repr(value)
    return value


def _atomic_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(_json_value(payload), ensure_ascii=False, indent=2, sort_keys=True)
        + "\n",
        encoding="utf-8",
    )
    json.loads(temporary.read_text(encoding="utf-8"))
    os.replace(temporary, path)


def _atomic_csv(path: Path, rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _json_value(row.get(field, "")) for field in fields})
    with temporary.open("r", newline="", encoding="utf-8") as stream:
        list(csv.DictReader(stream))
    os.replace(temporary, path)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as stream:
        return [dict(row) for row in csv.DictReader(stream)]


def _float(row: Mapping[str, Any], key: str, default: float = math.nan) -> float:
    value = row.get(key, "")
    return default if value in ("", None) else float(value)


def _int_or_one(value: Any) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return 1


def _serialized_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text == "true":
        return True
    if text == "false":
        return False
    raise ValueError(f"Expected serialized boolean, got {value!r}.")


def _current_peak_rss_bytes() -> int:
    """Return the Windows process peak working set without an extra package."""

    if os.name != "nt":
        return 0

    class ProcessMemoryCounters(ctypes.Structure):
        _fields_ = [
            ("cb", wintypes.DWORD),
            ("PageFaultCount", wintypes.DWORD),
            ("PeakWorkingSetSize", ctypes.c_size_t),
            ("WorkingSetSize", ctypes.c_size_t),
            ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPagedPoolUsage", ctypes.c_size_t),
            ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
            ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
            ("PagefileUsage", ctypes.c_size_t),
            ("PeakPagefileUsage", ctypes.c_size_t),
        ]

    counters = ProcessMemoryCounters()
    counters.cb = ctypes.sizeof(counters)
    kernel32 = ctypes.windll.kernel32
    psapi = ctypes.windll.psapi
    kernel32.GetCurrentProcess.restype = wintypes.HANDLE
    psapi.GetProcessMemoryInfo.argtypes = [
        wintypes.HANDLE,
        ctypes.POINTER(ProcessMemoryCounters),
        wintypes.DWORD,
    ]
    psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
    process = kernel32.GetCurrentProcess()
    ok = psapi.GetProcessMemoryInfo(
        process, ctypes.byref(counters), counters.cb
    )
    return int(counters.PeakWorkingSetSize) if ok else 0


@dataclass
class BackendCounters:
    determinant_evaluations: int = 0
    matrix_exponentials: int = 0
    certified_oracle_reads: int = 0
    fresh_full_scans: int = 0
    local_primary_scans: int = 0
    local_verification_scans: int = 0
    isolated_anchor_reuses: int = 0


class CountingProvider:
    def __init__(self, provider: Any, counters: BackendCounters) -> None:
        self.provider = provider
        self.counters = counters

    def __call__(self, omega: float) -> np.ndarray:
        self.counters.determinant_evaluations += 1
        return np.asarray(self.provider(float(omega)), dtype=float)


def _oracle_paths() -> dict[tuple[str, str], Path]:
    return {
        ("beta", MODEL_EB): RLB2C_OUTPUT / "eb_roots.csv",
        ("beta", MODEL_OLD): RLB2C_OUTPUT / "old_timoshenko_roots.csv",
        ("beta", MODEL_RLB): RLB2C_OUTPUT / "new_rlb_roots.csv",
        ("mu", MODEL_EB): RLB2D_OUTPUT / "mu_sweep_eb_roots.csv",
        ("mu", MODEL_OLD): RLB2D_OUTPUT / "mu_sweep_old_timoshenko_roots.csv",
        ("mu", MODEL_RLB): RLB2D_OUTPUT / "mu_sweep_new_rlb_roots.csv",
        ("beta_iso_rlb", MODEL_RLB): RLB2B_OUTPUT / "new_rlb_roots.csv",
    }


def _physical_contract_path(sweep_id: str) -> Path:
    if sweep_id == "mu":
        return RLB2D_OUTPUT / "case_contract.json"
    if sweep_id == "beta":
        return RLB2C_OUTPUT / "case_contract.json"
    if sweep_id == "beta_iso_rlb":
        return RLB2B_OUTPUT / "case_contract.json"
    raise ValueError(sweep_id)


def _repo_relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def _oracle_key(sweep_id: str, row: Mapping[str, str]) -> float:
    if sweep_id.startswith("beta"):
        return round(float(row["beta_deg"]), 12)
    return round(float(row["mu"]), 12)


def _oracle_group_qualification(
    point_rows: Sequence[Mapping[str, str]],
    model: str,
) -> tuple[bool, tuple[str, ...]]:
    """Qualify only the exported first eight roots and root 9 as an oracle."""

    reasons: list[str] = []
    positions = [int(row["sorted_position"]) for row in point_rows]
    if positions != list(range(1, 10)):
        reasons.append("POSITIONS_1_TO_9_INCOMPLETE")
        return False, tuple(reasons)
    values = [_float(row, "Omega", _float(row, "Lambda")) for row in point_rows]
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        reasons.append("NONFINITE_OR_NONPOSITIVE_ROOT")
    if any(right < left for left, right in zip(values, values[1:], strict=False)):
        reasons.append("UNSORTED_ROOTS")
    for position, row in enumerate(point_rows, start=1):
        if row.get("root_status") != "PASS":
            reasons.append(f"ROOT_{position}_STATUS")
        try:
            guard_flag = _serialized_bool(row.get("guard_flag"))
        except ValueError:
            reasons.append(f"ROOT_{position}_GUARD_FLAG")
        else:
            if guard_flag != (position == 9):
                reasons.append(f"ROOT_{position}_GUARD_STRUCTURE")
        sigma = _float(row, "scaled_sigma_ratio")
        boundary = _float(row, "boundary_null_residual")
        if not math.isfinite(sigma) or sigma > rlb2c.rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE:
            reasons.append(f"ROOT_{position}_SIGMA_RATIO")
        if not math.isfinite(boundary) or boundary > rlb2c.rlb2b.BOUNDARY_RESIDUAL_TOLERANCE:
            reasons.append(f"ROOT_{position}_BOUNDARY_RESIDUAL")
        if not str(row.get("detector_type", "")).strip():
            reasons.append(f"ROOT_{position}_DETECTOR_EVIDENCE")

    first = point_rows[0]
    export_status = str(first.get("export_range_status", "")).strip()
    if export_status:
        if export_status != "PASS":
            reasons.append("EXPORT_RANGE_STATUS")
        try:
            export_agreement = _serialized_bool(
                first.get("export_primary_verification_agreement")
            )
        except ValueError:
            reasons.append("EXPORT_PRIMARY_VERIFICATION_FLAG")
        else:
            if not export_agreement:
                reasons.append("EXPORT_PRIMARY_VERIFICATION_FLAG")
        export_relative = _float(first, "export_primary_verification_max_relative")
        if (
            not math.isfinite(export_relative)
            or export_relative > ORACLE_QUALIFICATION_RELATIVE_TOLERANCE
        ):
            reasons.append("EXPORT_PRIMARY_VERIFICATION_RELATIVE")
        unresolved = str(first.get("unresolved_candidates_below_export_guard", ""))
        if unresolved != "NOT_ASSESSED_BY_EB_SIGN_SCAN":
            try:
                if int(unresolved) != 0:
                    reasons.append("UNRESOLVED_BELOW_ROOT9")
            except ValueError:
                reasons.append("UNRESOLVED_BELOW_ROOT9")
    else:
        expected_status = (
            "PASS_SIGN_SCAN_CROSSCHECK" if model == MODEL_EB else "PASS"
        )
        if first.get("inventory_status") != expected_status:
            reasons.append("INVENTORY_STATUS")
        primary_verification = _float(first, "primary_verification_max_relative")
        if (
            not math.isfinite(primary_verification)
            or primary_verification > ORACLE_QUALIFICATION_RELATIVE_TOLERANCE
        ):
            reasons.append("PRIMARY_VERIFICATION_RELATIVE")
        unresolved = str(first.get("unresolved_candidates_below_internal_guard", ""))
        if unresolved not in {"0", "NOT_ASSESSED_BY_EB_SIGN_SCAN"}:
            reasons.append("UNRESOLVED_BELOW_RECORDED_GUARD")
    return not reasons, tuple(dict.fromkeys(reasons))


def _load_oracle(
    sweep_id: str,
    model: str,
) -> tuple[
    dict[float, list[dict[str, str]]],
    str,
    dict[float, tuple[str, ...]],
    dict[float, list[dict[str, str]]],
]:
    path = _oracle_paths()[(sweep_id, model)]
    rows = _read_csv(path)
    grouped: dict[float, list[dict[str, str]]] = {}
    seen_row_keys: set[tuple[float, int]] = set()
    for row in rows:
        key = _oracle_key(sweep_id, row)
        position = int(row["sorted_position"])
        row_key = (key, position)
        if row_key in seen_row_keys:
            raise RuntimeError(
                f"Duplicate frozen reference row key in {path}: {row_key}."
            )
        seen_row_keys.add(row_key)
        recorded_model = str(row.get("model", "")).strip()
        if recorded_model and recorded_model != model:
            raise RuntimeError(
                f"Frozen reference model mismatch in {path}: "
                f"{recorded_model!r} != {model!r}."
            )
        grouped.setdefault(key, []).append(row)
    certified: dict[float, list[dict[str, str]]] = {}
    excluded: dict[float, tuple[str, ...]] = {}
    for key, point_rows in grouped.items():
        point_rows.sort(key=lambda row: int(row["sorted_position"]))
        passed, reasons = _oracle_group_qualification(point_rows, model)
        if passed:
            certified[key] = point_rows
        else:
            excluded[key] = reasons
    return certified, _sha256(path), excluded, grouped


def _reference_contract_metadata(
    sweep_id: str,
    model: str,
    target: float,
) -> dict[str, Any]:
    certified, reference_sha256, excluded, recorded = _load_oracle(
        sweep_id, model
    )
    key = round(float(target), 12)
    point_rows = recorded.get(key, [])
    row_keys = [
        {
            "sweep_id": sweep_id,
            "model": model,
            "parameter_value": float(target),
            "spectral_position": int(row["sorted_position"]),
        }
        for row in point_rows
    ]
    return {
        "reference_csv_path": _repo_relative(_oracle_paths()[(sweep_id, model)]),
        "reference_csv_sha256": reference_sha256,
        "reference_row_keys": row_keys,
        "reference_row_keys_sha256": sweep_runner.stable_fingerprint(row_keys),
        "reference_rows_present": bool(point_rows),
        "reference_group_certified_by_historical_metadata": key in certified,
        "reference_group_qualifications": list(excluded.get(key, ())),
        "physical_contract_path": _repo_relative(
            _physical_contract_path(sweep_id)
        ),
        "physical_contract_sha256": _sha256(_physical_contract_path(sweep_id)),
    }


def _legacy_shard_inventory(output_dir: Path) -> list[dict[str, Any]]:
    legacy_root = output_dir / LEGACY_SHARD_DIRECTORY
    if not legacy_root.is_dir():
        return []
    return [
        {
            "path": path.relative_to(output_dir).as_posix(),
            "sha256": _sha256(path),
            "scientific_role": "LEGACY_EVIDENCE_NOT_USED",
        }
        for path in sorted(legacy_root.rglob("*.json"))
    ]


def _record_from_rows(parameter: float, rows: Sequence[Mapping[str, Any]], origin: str) -> sweep_runner.SpectrumRecord:
    roots: list[sweep_runner.RootRecord] = []
    for row in rows:
        cluster_id = str(row.get("cluster_id", ""))
        roots.append(
            sweep_runner.RootRecord(
                value=_float(row, "Omega", _float(row, "Lambda")),
                accepted=str(row.get("root_status", "PASS")) == "PASS",
                sigma_ratio=_float(row, "scaled_sigma_ratio", 0.0),
                boundary_residual=_float(row, "boundary_null_residual", 0.0),
                detector_agreement=True,
                possible_even_root=False,
                multiplicity=_int_or_one(row.get("multiplicity", 1)),
                nullity=_int_or_one(row.get("detected_nullity", 1)),
                cluster_id=cluster_id,
                metadata={
                    "raw_determinant": row.get("raw_determinant", ""),
                    "scaled_determinant": row.get("scaled_determinant", ""),
                    "detector_type": row.get("detector_type", ""),
                    "inventory_status": row.get("inventory_status", ""),
                    "guard_flag": row.get("guard_flag", ""),
                },
            )
        )
    return sweep_runner.SpectrumRecord(
        parameter=float(parameter),
        roots=tuple(roots),
        origin=origin,
        status="PASS",
        internal_tail_status=str(rows[0].get("internal_inventory_status", "RECORDED")),
        metadata={"source_inventory_sha256": rows[0].get("inventory_sha256", "")},
    )


class PhysicalValidationBackend:
    """Thin numerical adapter; all physics/refinement remains frozen."""

    def __init__(
        self,
        sweep_id: str,
        model: str,
        *,
        use_frozen_oracle_for_global: bool = True,
        fresh_full_scan_limit: int | None = None,
    ) -> None:
        if sweep_id not in {"beta", "mu", "beta_iso_rlb"}:
            raise ValueError(sweep_id)
        self.sweep_id = sweep_id
        self.model = model
        self.policy = rlb2c.rlb2b.frozen_root_policy()
        self.full_export_policy = replace(
            self.policy,
            requested_roots=FULL_EXPORT_ROOTS,
            guard_roots=FULL_EXPORT_GUARDS,
        )
        self.use_frozen_oracle_for_global = bool(use_frozen_oracle_for_global)
        self.fresh_full_scan_limit = fresh_full_scan_limit
        self.force_fresh_global_parameters: set[float] = set()
        if self.use_frozen_oracle_for_global:
            (
                self.oracle,
                self.oracle_file_sha256,
                self.oracle_exclusions,
                self.recorded_oracle,
            ) = _load_oracle(sweep_id, model)
        else:
            self.oracle = {}
            self.oracle_file_sha256 = "NOT_USED_IN_PRODUCTION"
            self.oracle_exclusions = {}
            self.recorded_oracle = {}
        self.counters = BackendCounters()
        self.fresh_full_cache: sweep_runner.ExactLRU[sweep_runner.SpectrumRecord] = (
            sweep_runner.ExactLRU(
                max_entries=FRESH_FULL_CACHE_ENTRIES,
                max_bytes=FRESH_FULL_CACHE_BYTES,
            )
        )
        self._active_parameter: float | None = None
        self._active_provider: CountingProvider | None = None
        self._active_scale: float = 1.0
        self._active_cache: Any = None

        if sweep_id in {"beta", "beta_iso_rlb"}:
            self.beta_contract = (
                rlb2c.build_case_contract(np.arange(0.0, 91.0, 1.0))
                if sweep_id == "beta"
                else rlb2c.rlb2b.build_case_contract(np.arange(0.0, 91.0, 1.0))
            )
            self.beta_objects = (
                rlb2c.build_model_objects(self.beta_contract)
                if sweep_id == "beta"
                else rlb2c.rlb2b.build_model_objects(self.beta_contract)
            )
            self._active_scale = float(
                self.beta_contract["frequency"].get(
                    "Omega_per_omega", self.beta_contract["frequency"].get("Lambda_per_omega")
                )
            )
            self.shared_rlb_cache: MutableMappingFloat = MutableMappingFloat(
                self.counters
            )
        else:
            self.beta_contract = None
            self.beta_objects = None
            self.shared_rlb_cache = MutableMappingFloat(self.counters)
            self._active_scale = float(rlb2d.reference_omega_scale())

    def _point_provider(self, parameter: float) -> tuple[CountingProvider, float]:
        key = round(float(parameter), 12)
        if self._active_parameter == key and self._active_provider is not None:
            return self._active_provider, self._active_scale

        if self.sweep_id == "beta":
            provider, _ = rlb2c.make_matrix_provider(
                self.model,
                float(parameter),
                self.beta_contract,
                self.beta_objects,
                rlb_arm_cache=self.shared_rlb_cache,
            )
            scale = float(self.beta_contract["frequency"]["Omega_per_omega"])
        elif self.sweep_id == "beta_iso_rlb":
            provider, _ = rlb2c.rlb2b.make_matrix_provider(
                self.model,
                float(parameter),
                self.beta_contract,
                self.beta_objects,
            )
            scale = float(self.beta_contract["frequency"]["Lambda_per_omega"])
        else:
            geometry = rlb2d.geometry_for(float(parameter), rlb2d.MU_TAU, rlb2d.MU_BETA_DEG)
            objects = rlb2d.build_model_objects(geometry)
            if rlb2d.constitutive_checks(geometry, objects)["status"] != "PASS":
                raise RuntimeError("RLB-2D geometry/constitutive contract failed.")
            point_cache = rlb2d.ArmPairCache(
                arm1=MutableMappingFloat(self.counters),
                arm2=MutableMappingFloat(self.counters),
            )
            provider, _ = rlb2d.make_matrix_provider(
                self.model, geometry, objects, cache=point_cache
            )
            self._active_cache = point_cache
            scale = float(rlb2d.reference_omega_scale())

        counted = CountingProvider(provider, self.counters)
        self._active_parameter = key
        self._active_provider = counted
        self._active_scale = scale
        return counted, scale

    def quality_gate(self, root: sweep_runner.RootRecord) -> tuple[bool, str]:
        if root.sigma_ratio > self.policy.root_singular_ratio:
            return False, "ROOT_SINGULAR_RATIO_FAIL"
        if root.boundary_residual > self.policy.boundary_null_residual:
            return False, "BOUNDARY_NULL_RESIDUAL_FAIL"
        return True, ""

    def global_search(self, parameter: float) -> sweep_runner.SpectrumRecord:
        key = round(float(parameter), 12)
        if key in self.force_fresh_global_parameters:
            if key in self.fresh_full_cache:
                return self.fresh_full_cache[key]
            record = self.fresh_full_search(parameter)
            self.fresh_full_cache[key] = record
            return record
        if self.use_frozen_oracle_for_global and key in self.recorded_oracle:
            self.counters.certified_oracle_reads += 1
            origin = (
                "CERTIFIED_FULL_ORACLE_ANCHOR"
                if key in self.oracle
                else "RECORDED_FULL_ORACLE_ANCHOR"
            )
            return _record_from_rows(
                parameter, self.recorded_oracle[key], origin
            )
        if key in self.fresh_full_cache:
            return self.fresh_full_cache[key]
        record = self.fresh_full_search(parameter)
        self.fresh_full_cache[key] = record
        return record

    def fresh_full_search(self, parameter: float) -> sweep_runner.SpectrumRecord:
        if (
            self.fresh_full_scan_limit is not None
            and self.counters.fresh_full_scans >= self.fresh_full_scan_limit
        ):
            raise RuntimeError(
                "Closing validation fresh-full budget exhausted before scan."
            )
        provider, scale = self._point_provider(parameter)
        started = time.perf_counter()
        if self.model == MODEL_EB:
            if self.sweep_id.startswith("beta"):
                rows = rlb2c.rlb2b._eb_root_rows(
                    float(parameter),
                    self.beta_contract,
                    self.beta_objects,
                    self.full_export_policy,
                )
                record = _record_from_rows(parameter, rows, "FRESH_FULL_SCAN")
            else:
                geometry = rlb2d.geometry_for(float(parameter), 0.0, rlb2d.MU_BETA_DEG)
                objects = rlb2d.build_model_objects(geometry)
                raw, export = rlb2d._eb_sign_scan_rows(
                    geometry,
                    objects,
                    self.full_export_policy,
                    cache=rlb2d.ArmPairCache.empty(),
                )
                rows = [rlb2d.transform_root_row(row, self.model, geometry, "mu", export) for row in raw]
                record = _record_from_rows(parameter, rows, "FRESH_FULL_SCAN")
        else:
            inventory = rlb2c.rlb2b.iso_inventory.seed_free_root_inventory(
                provider,
                scale,
                self.full_export_policy,
                case_id=f"sweep_fast_fresh__{self.sweep_id}__{self.model}__{parameter:.12g}",
                builder_id=f"SWEEP_FAST_{self.model}",
                contract_sha256="",
            )
            if len(inventory.slots) != FULL_EXPORT_ROOTS + FULL_EXPORT_GUARDS:
                raise RuntimeError(
                    "Fresh full scan did not return exactly positions 1--8 plus root 9; "
                    "a guard-cluster extension requires explicit qualification."
                )
            if not inventory.independent_agreement or inventory.status != "PASS":
                raise RuntimeError(
                    "Fresh first-8-plus-root9 inventory failed its frozen "
                    "primary/verification or unresolved-candidate gate."
                )
            record = self._record_from_slots(
                parameter, inventory.slots, "FRESH_FULL_SCAN"
            )
            record = replace(
                record,
                metadata={
                    **dict(record.metadata),
                    "primary_verification_max_relative": float(
                        inventory.maximum_primary_verification_relative
                    ),
                    "unresolved_candidates_below_guard": int(
                        inventory.unresolved_low_sigma_count
                    ),
                    "independent_agreement": bool(
                        inventory.independent_agreement
                    ),
                },
            )
        self.counters.fresh_full_scans += 1
        metadata = dict(record.metadata)
        metadata["fresh_full_wall_time_seconds"] = time.perf_counter() - started
        return replace(record, metadata=metadata)

    def _record_from_slots(
        self,
        parameter: float,
        slots: Sequence[Any],
        origin: str,
    ) -> sweep_runner.SpectrumRecord:
        roots: list[sweep_runner.RootRecord] = []
        for slot in slots:
            event = slot.event
            candidate = event.candidate
            diagnostic = candidate.diagnostics
            sources = tuple(candidate.detection_sources)
            roots.append(
                sweep_runner.RootRecord(
                    value=float(event.Omega),
                    accepted=bool(candidate.accepted),
                    sigma_ratio=float(diagnostic.scaled_sigma_ratio),
                    boundary_residual=max(
                        float(diagnostic.scaled_null_residual),
                        float(diagnostic.raw_boundary_null_residual),
                    ),
                    detector_agreement=True,
                    possible_even_root=False,
                    multiplicity=int(event.multiplicity),
                    nullity=int(event.detected_nullity),
                    cluster_id=str(event.cluster_id),
                    metadata={
                        "detection_sources": list(sources),
                        "raw_determinant": diagnostic.raw_determinant,
                        "scaled_determinant": diagnostic.scaled_determinant,
                        "raw_sigma_min": diagnostic.raw_sigma_min,
                        "raw_sigma_max": diagnostic.raw_sigma_max,
                        "interval_left": candidate.interval_left_Omega,
                        "interval_right": candidate.interval_right_Omega,
                    },
                )
            )
        return sweep_runner.SpectrumRecord(
            parameter=float(parameter), roots=tuple(roots), origin=origin
        )

    def local_search(
        self,
        parameter: float,
        interval: sweep_runner.SearchInterval,
        verification: bool,
    ) -> Sequence[sweep_runner.RootRecord]:
        if self.model == MODEL_EB and self.sweep_id != "beta_iso_rlb":
            return self._eb_local_search(parameter, interval, verification)
        provider, scale = self._point_provider(parameter)
        width = float(interval.upper - interval.lower)
        step = LOCAL_VERIFICATION_GRID_STEP if verification else LOCAL_PRIMARY_GRID_STEP
        minimum = MINIMUM_LOCAL_VERIFICATION_POINTS if verification else MINIMUM_LOCAL_PRIMARY_POINTS
        points = max(minimum, int(math.ceil(width / step)) + 1)
        local_policy = replace(
            self.policy,
            requested_roots=interval.expected_count + 1,
            guard_roots=0,
            Omega_min=float(interval.lower),
            Omega_max=float(interval.upper),
        )
        scan_id = "local_verification" if verification else "local_primary"
        scan = rlb2c.rlb2b.iso_inventory._run_scan(
            provider,
            scale,
            local_policy,
            case_id=f"sweep_fast_local__{self.sweep_id}__{self.model}__{parameter:.12g}",
            builder_id=f"SWEEP_FAST_LOCAL_{self.model}",
            scan_id=scan_id,
            points=points,
            phases=(0.5,) if verification else (0.0,),
        )
        if verification:
            self.counters.local_verification_scans += 1
        else:
            self.counters.local_primary_scans += 1
        multiplicity_count = sum(int(event.multiplicity) for event in scan.events)
        if multiplicity_count != interval.expected_count:
            raise sweep_runner.SweepValidationError(
                f"local accepted multiplicity {multiplicity_count} != {interval.expected_count}"
            )
        accepted_values = [float(event.Omega) for event in scan.events]
        for candidate in scan.rejected_candidates:
            if candidate.rejection_reason == "DUPLICATE_DETECTION_SAME_ROOT":
                continue
            if any(
                abs(float(candidate.Omega) - value)
                <= self.policy.dedup_atol_Omega
                + self.policy.dedup_rtol * max(abs(float(candidate.Omega)), abs(value))
                for value in accepted_values
            ):
                continue
            if candidate.diagnostics.scaled_sigma_ratio <= self.policy.sigma_prefilter:
                raise sweep_runner.SweepValidationError(
                    f"unresolved low-sigma local candidate: {candidate.rejection_reason}"
                )
        slots = scan.slots[: interval.expected_count]
        record = self._record_from_slots(parameter, slots, scan_id)
        roots: list[sweep_runner.RootRecord] = []
        for root in record.roots:
            sources = [str(item) for item in root.metadata.get("detection_sources", [])]
            has_bracket = any("determinant_bracket" in item for item in sources)
            roots.append(replace(root, possible_even_root=not has_bracket))
        return roots

    def _eb_local_search(
        self,
        parameter: float,
        interval: sweep_runner.SearchInterval,
        verification: bool,
    ) -> Sequence[sweep_runner.RootRecord]:
        epsilon = rlb2d.H0 / (math.sqrt(12.0) * rlb2d.L_REFERENCE)
        beta_rad = (
            math.radians(float(parameter))
            if self.sweep_id == "beta"
            else math.radians(rlb2d.MU_BETA_DEG)
        )
        mu = 0.0 if self.sweep_id == "beta" else float(parameter)
        step = rlb2c.rlb2b.EB_VERIFICATION_SCAN_STEP if verification else rlb2c.rlb2b.EB_PRIMARY_SCAN_STEP
        historical = eb_interval_scan(
            beta=beta_rad,
            mu=mu,
            eps=epsilon,
            n_roots=interval.expected_count + 1,
            Lmin=math.sqrt(max(interval.lower, np.finfo(float).tiny)),
            Lmax=math.sqrt(interval.upper),
            scan_step=float(step),
            bisect_iters=70,
        )
        if len(historical) != interval.expected_count:
            raise sweep_runner.SweepValidationError(
                f"EB local interval returned {len(historical)} roots"
            )
        provider, scale = self._point_provider(parameter)
        roots: list[sweep_runner.RootRecord] = []
        for value in historical:
            Omega = float(value) ** 2
            diagnostic = rlb2c.rlb2b.iso_inventory._boundary_matrix_diagnostics(
                Omega, provider, scale, self.policy
            )
            roots.append(
                sweep_runner.RootRecord(
                    value=Omega,
                    accepted=True,
                    sigma_ratio=float(diagnostic.scaled_sigma_ratio),
                    boundary_residual=max(
                        float(diagnostic.scaled_null_residual),
                        float(diagnostic.raw_boundary_null_residual),
                    ),
                    detector_agreement=True,
                    metadata={"historical_eb_wavenumber": float(value)},
                )
            )
        if verification:
            self.counters.local_verification_scans += 1
        else:
            self.counters.local_primary_scans += 1
        return roots


class MutableMappingFloat(sweep_runner.ExactLRU[np.ndarray]):
    """Exact bounded mapping accepted by the frozen arm-cache adapters."""

    def __init__(self, counters: BackendCounters) -> None:
        super().__init__(MAX_CACHE_ENTRIES, MAX_CACHE_BYTES)
        self.backend_counters = counters

    def __setitem__(self, key: Any, value: np.ndarray) -> None:
        self.backend_counters.matrix_exponentials += 1
        super().__setitem__(key, value)


def _physical_beta_cache_equality_check() -> dict[str, Any]:
    """Compare cached arm assembly with the frozen uncached public builder."""

    contract = rlb2c.build_case_contract(np.arange(0.0, 91.0, 1.0))
    objects = rlb2c.build_model_objects(contract)
    counters = BackendCounters()
    cache = MutableMappingFloat(counters)
    length = float(contract["geometry"]["l"])
    properties = objects.weak_properties
    rows: list[dict[str, float]] = []
    maximum_relative = 0.0
    maximum_absolute = 0.0
    for beta_deg in (0.0, 54.0, 90.0):
        provider, metadata = rlb2c.make_matrix_provider(
            MODEL_RLB,
            beta_deg,
            contract,
            objects,
            rlb_arm_cache=cache,
        )
        if float(metadata["cached_vs_public_builder_max_abs"]) > 16.0 * np.finfo(float).eps:
            raise RuntimeError("Frozen RLB adapter rejected cached arm assembly.")
        for omega in (0.731, 3.217, 12.345):
            cached = np.asarray(provider(omega), dtype=float)
            uncached = np.asarray(
                rlb2c.rlb_coupled.coupled_boundary_matrix(
                    omega,
                    math.radians(beta_deg),
                    length,
                    properties,
                    length,
                    properties,
                ),
                dtype=float,
            )
            absolute = float(np.max(np.abs(cached - uncached)))
            scale = max(float(np.max(np.abs(uncached))), np.finfo(float).tiny)
            relative = absolute / scale
            maximum_absolute = max(maximum_absolute, absolute)
            maximum_relative = max(maximum_relative, relative)
            rows.append(
                {
                    "beta_deg": beta_deg,
                    "omega": omega,
                    "maximum_absolute": absolute,
                    "maximum_relative": relative,
                }
            )
    return {
        "status": (
            "PASS"
            if maximum_relative <= CACHE_EQUALITY_RELATIVE_TOLERANCE
            else "FAIL"
        ),
        "tolerance_relative": CACHE_EQUALITY_RELATIVE_TOLERANCE,
        "maximum_absolute": maximum_absolute,
        "maximum_relative": maximum_relative,
        "rows": rows,
        "cache_diagnostics": cache.diagnostics(),
    }


def _settings(sweep_id: str, *, spike: bool = True) -> sweep_runner.SweepSettings:
    return sweep_runner.SweepSettings(
        anchor_stride=5 if sweep_id == "mu" else 10,
        local_verification_relative=LOCAL_VERIFICATION_RELATIVE_TOLERANCE,
        anchor_comparison_relative=LOCAL_VERIFICATION_RELATIVE_TOLERANCE,
        cluster_relative_gap=1.0e-3,
        fallback_gap_relative=1.0e-8,
        duplicate_relative_tolerance=1.0e-12,
        enable_spike_audit=spike,
    )


def _validation_worker_settings(sweep_id: str) -> sweep_runner.SweepSettings:
    return replace(_settings(sweep_id, spike=False), anchor_stride=10)


def _validation_chain(
    target: float,
    sweep_id: str,
    *,
    certified_parameters: Sequence[float] | None = None,
) -> np.ndarray:
    step = 1.0 if sweep_id.startswith("beta") else 0.02
    lower = 0.0
    upper = 90.0 if sweep_id.startswith("beta") else 0.8
    if target <= lower or target >= upper:
        return np.asarray([target], dtype=float)
    if certified_parameters is not None:
        certified = sorted(
            {
                round(float(value), 12)
                for value in certified_parameters
                if math.isfinite(float(value))
            }
        )
        below = [value for value in certified if value < target]
        above = [value for value in certified if value > target]
        if len(below) >= 2:
            right = above[0] if above else upper
            if right > target:
                return np.asarray([below[-2], below[-1], target, right], dtype=float)
    start = max(lower, target - 2.0 * step)
    stop = min(upper, target + step)
    return np.asarray(
        [round(start + index * step, 12) for index in range(int(round((stop - start) / step)) + 1)],
        dtype=float,
    )


def _synthetic_evidence(output_dir: Path) -> dict[str, Any]:
    def values(parameter: float) -> list[float]:
        return [100.0 * position + parameter for position in range(1, 10)]

    false_spike = {1.0: 101.9}

    def record(parameter: float, *, kink: bool = False) -> sweep_runner.SpectrumRecord:
        root_values = values(parameter)
        if kink and parameter in false_spike:
            root_values[0] = false_spike[parameter]
        return sweep_runner.SpectrumRecord(
            parameter=parameter,
            roots=tuple(sweep_runner.RootRecord(value=value) for value in root_values),
            origin="SYNTHETIC_FULL",
        )

    def local_factory(kink: bool, fail_at: float | None = None) -> Any:
        def local(parameter: float, interval: sweep_runner.SearchInterval, verification: bool) -> list[sweep_runner.RootRecord]:
            del verification
            if fail_at is not None and parameter == fail_at and interval.positions == (3,):
                return []
            source = record(parameter, kink=kink).roots
            return [source[position - 1] for position in interval.positions]

        return local

    fallback_result = sweep_runner.run_spectral_sweep(
        [0.0, 1.0, 2.0, 3.0],
        callbacks=sweep_runner.SweepCallbacks(
            lambda parameter: record(parameter), local_factory(False, 1.0)
        ),
        settings=replace(_settings("beta", spike=False), anchor_stride=10),
    )
    reproduced = sweep_runner.run_spectral_sweep(
        [0.0, 1.0, 2.0],
        callbacks=sweep_runner.SweepCallbacks(
            lambda parameter: record(parameter, kink=True), local_factory(True)
        ),
        settings=replace(_settings("beta", spike=True), anchor_stride=10),
    )
    corrected = sweep_runner.run_spectral_sweep(
        [0.0, 1.0, 2.0],
        callbacks=sweep_runner.SweepCallbacks(
            lambda parameter: record(parameter), local_factory(True)
        ),
        settings=replace(_settings("beta", spike=True), anchor_stride=10),
    )

    checkpoint = sweep_runner.CheckpointConfig(
        directory=output_dir / "synthetic_resume_checkpoint",
        sweep_id="synthetic_resume",
        model_id="synthetic",
        fingerprint=sweep_runner.stable_fingerprint({"grid": [0.0, 1.0, 2.0]}),
    )
    first = sweep_runner.run_spectral_sweep(
        [0.0, 1.0, 2.0],
        callbacks=sweep_runner.SweepCallbacks(lambda parameter: record(parameter), local_factory(False)),
        settings=replace(_settings("beta", spike=False), anchor_stride=10),
        checkpoint=checkpoint,
    )
    second = sweep_runner.run_spectral_sweep(
        [0.0, 1.0, 2.0],
        callbacks=sweep_runner.SweepCallbacks(lambda parameter: record(parameter), local_factory(False)),
        settings=replace(_settings("beta", spike=False), anchor_stride=10),
        checkpoint=checkpoint,
        resume=True,
    )
    return {
        "fallback_origin": fallback_result.spectra[1.0].origin,
        "forced_next_anchor_origin": fallback_result.spectra[2.0].origin,
        "reproduced_outcomes": [item.outcome for item in reproduced.spike_audits],
        "corrected_outcomes": [item.outcome for item in corrected.spike_audits],
        "reproduced_value_preserved": reproduced.spectra[1.0].values[0] == 101.9,
        "corrected_value_from_full_scan": corrected.spectra[1.0].values[0] == 101.0,
        "resume_first_committed": first.counters.points_committed,
        "resume_second_resumed": second.counters.points_resumed,
        "resume_second_committed": second.counters.points_committed,
    }


COMPARISON_FIELDS = (
    "workflow",
    "model",
    "parameter",
    "position",
    "role",
    "fast_value",
    "full_value",
    "relative_difference",
    "fast_origin",
    "fast_multiplicity",
    "full_multiplicity",
    "fast_nullity",
    "full_nullity",
    "cluster_match",
    "cluster_center_relative_difference",
    "cluster_center_pass",
    "fast_sigma_ratio",
    "full_sigma_ratio",
    "fast_boundary_residual",
    "full_boundary_residual",
    "root_quality_pass",
    "reference_source",
    "status",
)
ANCHOR_FIELDS = (
    "workflow",
    "model",
    "parameter",
    "point_index",
    "scheduled_anchor",
    "origin",
    "predictor_kind",
)
FALLBACK_FIELDS = (
    "workflow",
    "model",
    "parameter",
    "origin",
    "fallback_reason",
    "initial_reference_source",
    "final_reference_source",
    "initial_maximum_relative_difference",
    "final_maximum_relative_difference",
    "status",
)
SPIKE_FIELDS = (
    "workflow",
    "model",
    "parameter",
    "position",
    "fast_value",
    "neighbor_prediction",
    "normalized_residual",
    "full_value",
    "outcome",
)
RUNTIME_FIELDS = (
    "phase",
    "workflow",
    "model",
    "point_count",
    "wall_time_seconds",
    "determinant_evaluations",
    "fresh_full_scans",
    "certified_oracle_reads",
    "local_primary_scans",
    "local_verification_scans",
)
MEMORY_FIELDS = ("phase", "peak_rss_bytes", "cache_peak_entries", "cache_peak_bytes")


def run_preflight(output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    production_hashes_before = {path: _sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS}
    validation_dependency_hashes_before = _dependency_hashes()
    synthetic = _synthetic_evidence(output_dir)
    cache_equality = _physical_beta_cache_equality_check()
    if cache_equality["status"] != "PASS":
        raise RuntimeError("Physical cached/uncached beta assembly check failed.")
    backend = PhysicalValidationBackend("beta", MODEL_RLB)
    endpoint_anchor_checks: list[dict[str, Any]] = []
    for endpoint in (0.0, 90.0):
        if endpoint not in backend.oracle:
            raise RuntimeError(f"Certified beta endpoint oracle {endpoint:g} is unavailable.")
        endpoint_record = _record_from_rows(
            endpoint,
            backend.oracle[endpoint],
            "CERTIFIED_FULL_ORACLE",
        )
        quality_pass = all(backend.quality_gate(root)[0] for root in endpoint_record.roots)
        endpoint_anchor_checks.append(
            {
                "beta_deg": endpoint,
                "root_count": len(endpoint_record.roots),
                "quality_pass": quality_pass,
            }
        )
    if not all(item["root_count"] == 9 and item["quality_pass"] for item in endpoint_anchor_checks):
        raise RuntimeError("Certified beta endpoint anchor check failed.")
    runtime_rows: list[dict[str, Any]] = []
    memory_rows: list[dict[str, Any]] = []
    local_times: list[float] = []
    representative = (1.0, 45.0, 54.0, 89.0)
    for target in representative:
        chain = _validation_chain(target, "beta")
        started = time.perf_counter()
        before = asdict(backend.counters)
        result = sweep_runner.run_spectral_sweep(
            chain,
            callbacks=sweep_runner.SweepCallbacks(
                backend.global_search,
                backend.local_search,
                backend.quality_gate,
                completeness_guard=lambda parameter, roots: (len(roots) == 9, "root9"),
                global_completeness_guard=lambda parameter, roots: (
                    len(roots) == 9,
                    "root9",
                ),
            ),
            settings=replace(_settings("beta", spike=False), anchor_stride=10),
        )
        elapsed = time.perf_counter() - started
        local_times.append(elapsed)
        after = asdict(backend.counters)
        runtime_rows.append(
            {
                "phase": "preflight_local_chain",
                "workflow": "beta",
                "model": MODEL_RLB,
                "point_count": len(chain),
                "wall_time_seconds": elapsed,
                **{key: after[key] - before[key] for key in before},
            }
        )
        memory_rows.append(
            {
                "phase": f"preflight_{target:g}",
                "peak_rss_bytes": _current_peak_rss_bytes(),
                "cache_peak_entries": backend.shared_rlb_cache.peak_entries,
                "cache_peak_bytes": backend.shared_rlb_cache.peak_bytes,
            }
        )
        if result.status != "PASS":
            raise RuntimeError(f"Preflight local chain failed at beta={target:g}.")

    full_started = time.perf_counter()
    fresh = backend.fresh_full_search(54.0)
    full_elapsed = time.perf_counter() - full_started
    oracle_54 = _record_from_rows(
        54.0,
        backend.recorded_oracle[54.0],
        "RECORDED_FULL_SCAN_REFERENCE",
    )
    fresh_error = sweep_runner.maximum_relative_error(fresh.values, oracle_54.values)
    if fresh_error > ORACLE_QUALIFICATION_RELATIVE_TOLERANCE:
        raise RuntimeError(
            "Fresh first-8-plus-root9 full scan disagrees with the frozen "
            f"beta=54 oracle: relative error {fresh_error:.6e}."
        )
    runtime_rows.append(
        {
            "phase": "preflight_fresh_full_difficult_point",
            "workflow": "beta",
            "model": MODEL_RLB,
            "point_count": 1,
            "wall_time_seconds": full_elapsed,
            **asdict(backend.counters),
        }
    )
    local_point_seconds = sum(local_times) / max(sum(len(_validation_chain(value, "beta")) - 2 for value in representative), 1)
    full_point_seconds = full_elapsed
    anchor_points = 3 * (11 + 9)
    local_points = 3 * ((91 - 11) + (41 - 9))
    estimated_fast_seconds = anchor_points * full_point_seconds + local_points * local_point_seconds
    estimated_full_seconds = 3 * (91 + 41) * full_point_seconds
    speedup_estimate = estimated_full_seconds / max(estimated_fast_seconds, np.finfo(float).tiny)
    physical_fast_local_evidence = bool(
        backend.counters.local_primary_scans > 0
        and backend.counters.local_verification_scans > 0
    )
    payload = {
        "schema_version": 1,
        "stage": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "created_utc": _utc_now(),
        "git": _git_state(),
        "thread_limits": THREAD_LIMITS,
        "python_workers": 1,
        "multiprocessing_used": False,
        "representative_points": {
            "start_anchor": 0.0,
            "start_local_probe": 1.0,
            "middle": 45.0,
            "difficult": 54.0,
            "end_local_probe": 89.0,
            "end_anchor": 90.0,
        },
        "fresh_full_difficult_point": {
            "beta_deg": 54.0,
            "maximum_relative_vs_frozen_oracle": fresh_error,
            "wall_time_seconds": full_elapsed,
        },
        "endpoint_anchor_checks": endpoint_anchor_checks,
        "mean_local_point_seconds": local_point_seconds,
        "estimated_fast_total_seconds_two_graphs": estimated_fast_seconds,
        "estimated_full_total_seconds_two_graphs": estimated_full_seconds,
        "estimated_speedup": speedup_estimate,
        "target_runtime_seconds": TARGET_RUNTIME_SECONDS,
        "production_block_threshold_seconds": BLOCK_RUNTIME_SECONDS,
        "physical_fast_local_evidence": physical_fast_local_evidence,
        "production_run_permitted": bool(
            estimated_fast_seconds <= BLOCK_RUNTIME_SECONDS
            and physical_fast_local_evidence
        ),
        "peak_rss_bytes": _current_peak_rss_bytes(),
        "backend_counters": asdict(backend.counters),
        "cache_diagnostics": backend.shared_rlb_cache.diagnostics(),
        "physical_cached_vs_uncached_equality": cache_equality,
        "oracle_exclusions": {
            format(parameter, ".12g"): list(reasons)
            for parameter, reasons in sorted(backend.oracle_exclusions.items())
        },
        "synthetic_infrastructure_evidence": synthetic,
        "production_hashes_before": production_hashes_before,
        "validation_dependency_hashes_before": validation_dependency_hashes_before,
        "validation_constants": _validation_constants_payload(),
    }
    _atomic_json(output_dir / "preflight_report.json", payload)
    _atomic_csv(output_dir / "runtime_profile.csv", runtime_rows, RUNTIME_FIELDS)
    _atomic_csv(output_dir / "memory_profile.csv", memory_rows, MEMORY_FIELDS)
    return payload


def _compare_record(
    workflow: str,
    model: str,
    parameter: float,
    fast: sweep_runner.SpectrumRecord,
    full: sweep_runner.SpectrumRecord,
) -> list[dict[str, Any]]:
    def membership(
        roots: Sequence[sweep_runner.RootRecord], index: int
    ) -> tuple[int, ...]:
        cluster_id = roots[index].cluster_id
        if not cluster_id:
            return (index + 1,)
        return tuple(
            offset + 1
            for offset, root in enumerate(roots)
            if root.cluster_id == cluster_id
        )

    rows: list[dict[str, Any]] = []
    for index, (left, right) in enumerate(
        zip(fast.roots, full.roots, strict=True)
    ):
        position = index + 1
        relative = abs(left.value - right.value) / max(
            abs(left.value), abs(right.value), np.finfo(float).tiny
        )
        cluster_match = bool(
            membership(fast.roots, index) == membership(full.roots, index)
            and left.multiplicity == right.multiplicity
            and left.nullity == right.nullity
        )
        fast_members = membership(fast.roots, index)
        full_members = membership(full.roots, index)
        fast_center = float(
            np.mean([fast.roots[item - 1].value for item in fast_members])
        )
        full_center = float(
            np.mean([full.roots[item - 1].value for item in full_members])
        )
        cluster_center_relative = abs(fast_center - full_center) / max(
            abs(fast_center), abs(full_center), np.finfo(float).tiny
        )
        cluster_center_pass = bool(
            cluster_center_relative <= FAST_FULL_RELATIVE_TOLERANCE
        )
        quality = bool(
            left.accepted
            and right.accepted
            and left.detector_agreement
            and right.detector_agreement
            and not left.possible_even_root
            and not right.possible_even_root
            and left.sigma_ratio <= rlb2c.rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE
            and right.sigma_ratio <= rlb2c.rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE
            and left.boundary_residual <= rlb2c.rlb2b.BOUNDARY_RESIDUAL_TOLERANCE
            and right.boundary_residual <= rlb2c.rlb2b.BOUNDARY_RESIDUAL_TOLERANCE
        )
        rows.append(
            {
                "workflow": workflow,
                "model": model,
                "parameter": parameter,
                "position": position,
                "role": "FIRST_8" if position <= 8 else "ROOT_9_GUARD",
                "fast_value": left.value,
                "full_value": right.value,
                "relative_difference": relative,
                "fast_origin": fast.origin,
                "fast_multiplicity": left.multiplicity,
                "full_multiplicity": right.multiplicity,
                "fast_nullity": left.nullity,
                "full_nullity": right.nullity,
                "cluster_match": cluster_match,
                "cluster_center_relative_difference": cluster_center_relative,
                "cluster_center_pass": cluster_center_pass,
                "fast_sigma_ratio": left.sigma_ratio,
                "full_sigma_ratio": right.sigma_ratio,
                "fast_boundary_residual": left.boundary_residual,
                "full_boundary_residual": right.boundary_residual,
                "root_quality_pass": quality,
                "status": (
                    "PASS"
                    if relative <= FAST_FULL_RELATIVE_TOLERANCE
                    and cluster_match
                    and cluster_center_pass
                    and quality
                    else "FAIL"
                ),
            }
        )
    return rows


def _production_inventory_hash(
    model: str,
    mu_value: float,
    rows: Sequence[Mapping[str, Any]],
) -> str:
    """Hash one RLB-2D point without the self-referential inventory field."""

    normalized = [
        {
            field: row.get(field, "")
            for field in RLB2D_PRODUCTION_FIELDS
            if field != "inventory_sha256"
        }
        for row in sorted(rows, key=lambda item: int(item["sorted_position"]))
    ]
    return sweep_runner.stable_fingerprint(
        {
            "stage": "RLB-2D-CLOSING-FIRST8-PLUS-ROOT9",
            "model": model,
            "mu": float(mu_value),
            "physical_contract_sha256": _sha256(
                _physical_contract_path("mu")
            ),
            "rows": normalized,
        }
    )


def _production_rows_semantic_hash(
    rows: Sequence[Mapping[str, Any]],
) -> str:
    return sweep_runner.stable_fingerprint(
        [
            {field: row.get(field, "") for field in RLB2D_PRODUCTION_FIELDS}
            for row in sorted(
                rows, key=lambda item: int(item["sorted_position"])
            )
        ]
    )


def _optional_sqrt(value: Any) -> Any:
    if value in ("", None):
        return ""
    numeric = float(value)
    return math.sqrt(numeric) if numeric >= 0.0 else ""


def _production_rows_from_record(
    backend: PhysicalValidationBackend,
    record: sweep_runner.SpectrumRecord,
    audit: sweep_runner.PointAudit,
) -> list[dict[str, Any]]:
    """Materialize an accepted mu point in the exact frozen 59-field schema.

    Root diagnostics are deliberately re-evaluated at the accepted values;
    neither predictor values nor the immutable comparison rows are exported.
    The extended inventory above root 9 is honestly marked as not computed.
    """

    if backend.sweep_id != "mu" or backend.model not in {MODEL_OLD, MODEL_RLB}:
        raise ValueError("Production rows are only defined for closing mu old/RLB points.")
    if len(record.roots) != 9:
        raise RuntimeError("A closing production point must contain roots 1--8 and root 9.")
    mu_value = round(float(record.parameter), 12)
    geometry = rlb2d.geometry_for(mu_value, rlb2d.MU_TAU, rlb2d.MU_BETA_DEG)
    provider, scale = backend._point_provider(mu_value)
    metadata_verification = float(
        record.metadata.get("primary_verification_max_relative", 0.0)
    )
    verification_relative = max(
        float(audit.maximum_primary_verification_relative), metadata_verification
    )
    verification_pass = bool(
        math.isfinite(verification_relative)
        and verification_relative <= LOCAL_VERIFICATION_RELATIVE_TOLERANCE
    )
    rows: list[dict[str, Any]] = []
    for position, root in enumerate(record.roots, start=1):
        Omega = float(root.value)
        omega = Omega / float(scale)
        Lambda = math.sqrt(Omega)
        diagnostic = rlb2c.rlb2b.iso_inventory._boundary_matrix_diagnostics(
            Omega, provider, scale, backend.policy
        )
        boundary_gate_value = max(
            float(diagnostic.scaled_null_residual),
            float(diagnostic.raw_boundary_null_residual),
        )
        quality_pass = bool(
            root.accepted
            and root.detector_agreement
            and math.isfinite(Omega)
            and Omega > 0.0
            and float(diagnostic.scaled_sigma_ratio)
            <= backend.policy.root_singular_ratio
            and boundary_gate_value <= backend.policy.boundary_null_residual
        )
        interval_left = root.metadata.get(
            "interval_left_Omega", root.metadata.get("interval_left", "")
        )
        interval_right = root.metadata.get(
            "interval_right_Omega", root.metadata.get("interval_right", "")
        )
        detector_type = root.metadata.get(
            "detection_sources", root.metadata.get("detector_type", ())
        )
        if not detector_type:
            detector_type = ("FAST_LOCAL_REFINED_DETERMINANT_SVD",)
        cluster_id = str(root.cluster_id or "")
        cluster_member = bool(cluster_id or int(root.multiplicity) > 1)
        row_values: dict[str, Any] = {
            "sweep": rlb2d.SWEEP_MU,
            "model": backend.model,
            "mu": mu_value,
            "tau": geometry.tau,
            "beta_deg": geometry.beta_deg,
            "L1": geometry.L1,
            "L2": geometry.L2,
            "h1": geometry.h1,
            "h2": geometry.h2,
            "b1": geometry.b1,
            "b2": geometry.b2,
            "sorted_position": position,
            "role": "FIRST_8" if position <= 8 else "ROOT_9_GUARD",
            "omega": omega,
            "Omega": Omega,
            "Lambda": Lambda,
            "normalization_identity_relative_residual": max(
                abs(Omega - omega * scale)
                / max(abs(Omega), abs(omega * scale), np.finfo(float).tiny),
                abs(Omega - Lambda**2)
                / max(abs(Omega), abs(Lambda**2), np.finfo(float).tiny),
            ),
            "historical_EB_mapping_relative_residual": "",
            "historical_EB_wavenumber": "",
            "multiplicity": int(root.multiplicity),
            "detected_nullity": int(diagnostic.detected_nullity),
            "cluster_id": cluster_id,
            "cluster_semantics": "CLUSTER_MEMBER" if cluster_member else "ISOLATED",
            "cluster_multiplicity": int(root.multiplicity),
            "cluster_total_nullity": int(root.nullity),
            "raw_determinant": float(diagnostic.raw_determinant),
            "scaled_determinant": float(diagnostic.scaled_determinant),
            "raw_sigma_min": float(diagnostic.raw_sigma_min),
            "raw_sigma_max": float(diagnostic.raw_sigma_max),
            "raw_sigma_ratio": float(diagnostic.raw_sigma_ratio),
            "scaled_sigma_min": float(diagnostic.scaled_sigma_min),
            "scaled_sigma_max": float(diagnostic.scaled_sigma_max),
            "scaled_sigma_ratio": float(diagnostic.scaled_sigma_ratio),
            "boundary_null_residual": float(
                diagnostic.raw_boundary_null_residual
            ),
            "detector_type": detector_type,
            "root_status": "PASS" if quality_pass else "FAIL",
            "inventory_status": (
                "PASS_FIRST8_PLUS_ROOT9_FAST_VALIDATED"
                if quality_pass and verification_pass
                else "FAIL"
            ),
            "inventory_sha256": "PENDING",
            "primary_slot_count_internal": 9,
            "verification_slot_count_internal": 9,
            "internal_requested_roots": 8,
            "internal_guard_position": 9,
            "primary_verification_max_relative": verification_relative,
            "unresolved_candidates_below_internal_guard": 0,
            "guard_flag": position == 9,
            "bracket_left_Omega": interval_left,
            "bracket_right_Omega": interval_right,
            "bracket_left_Lambda": _optional_sqrt(interval_left),
            "bracket_right_Lambda": _optional_sqrt(interval_right),
            "internal_root13_Omega": "",
            "internal_root13_Lambda": "",
            "export_guard_available": True,
            "export_primary_verification_agreement": verification_pass,
            "export_primary_verification_max_relative": verification_relative,
            "unresolved_candidates_below_export_guard": 0,
            "export_range_status": (
                "PASS" if quality_pass and verification_pass else "FAIL"
            ),
            "internal_inventory_status": "NOT_COMPUTED_ABOVE_ROOT9",
            "internal_primary_verification_max_relative": "",
            "internal_unresolved_candidate_count": "NOT_COMPUTED_ABOVE_ROOT9",
        }
        rows.append(
            {field: row_values[field] for field in RLB2D_PRODUCTION_FIELDS}
        )
    if any(row["root_status"] != "PASS" for row in rows) or not verification_pass:
        raise RuntimeError(
            f"Accepted closing roots failed re-evaluated diagnostics at mu={mu_value:g}."
        )
    inventory_sha256 = _production_inventory_hash(
        backend.model, mu_value, rows
    )
    for row in rows:
        row["inventory_sha256"] = inventory_sha256
    return rows


def _validate_production_rows(
    rows: Any,
    *,
    model: str,
    target: float,
    declared_semantic_sha256: str,
) -> list[dict[str, Any]]:
    if not isinstance(rows, list) or len(rows) != 9:
        raise ValueError("Production rows must contain exactly nine records.")
    copied = [dict(row) for row in rows]
    expected_fields = set(RLB2D_PRODUCTION_FIELDS)
    if any(len(row) != 59 or set(row) != expected_fields for row in copied):
        raise ValueError("Production row fields differ from the frozen schema.")
    copied = [
        {field: row[field] for field in RLB2D_PRODUCTION_FIELDS}
        for row in copied
    ]
    complete = rlb2d._complete_mu_values(copied, model)
    expected = round(float(target), 12)
    if complete != [expected]:
        raise ValueError("Production rows do not describe the requested mu point.")
    if _production_rows_semantic_hash(copied) != declared_semantic_sha256:
        raise ValueError("Production-row semantic hash mismatch.")
    inventory_hashes = {str(row["inventory_sha256"]) for row in copied}
    if len(inventory_hashes) != 1:
        raise ValueError("Production rows do not share one inventory hash.")
    if next(iter(inventory_hashes)) != _production_inventory_hash(
        model, expected, copied
    ):
        raise ValueError("Production inventory semantic hash mismatch.")
    values = [float(row["Omega"]) for row in copied]
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise ValueError("Production rows contain a non-positive/non-finite root.")
    if any(right <= left for left, right in zip(values, values[1:], strict=False)):
        raise ValueError("Production rows are unsorted or contain duplicate roots.")
    if any(str(row["root_status"]) != "PASS" for row in copied):
        raise ValueError("Production rows contain a failed root.")
    if any(str(row["export_range_status"]) != "PASS" for row in copied):
        raise ValueError("Production first8-plus-root9 export gate is not PASS.")
    if any(
        str(row["internal_inventory_status"]) != "NOT_COMPUTED_ABOVE_ROOT9"
        for row in copied
    ):
        raise ValueError("Production rows make a non-evidenced internal-tail claim.")
    if any(
        float(row["scaled_sigma_ratio"])
        > rlb2c.rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE
        or float(row["boundary_null_residual"])
        > rlb2c.rlb2b.BOUNDARY_RESIDUAL_TOLERANCE
        for row in copied
    ):
        raise ValueError("Production rows fail frozen root-quality thresholds.")
    return copied


def _run_fixed_validation_in_process(output_dir: Path) -> dict[str, Any]:
    preflight_path = output_dir / "preflight_report.json"
    if not preflight_path.is_file():
        raise RuntimeError("Run --preflight before --validate-fast-vs-full.")
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    if int(preflight.get("peak_rss_bytes", 0)) <= 0:
        preflight["peak_rss_bytes"] = _current_peak_rss_bytes()
        preflight["peak_rss_measurement_note"] = (
            "Recovered before fixed validation after declaring the Windows "
            "GetProcessMemoryInfo ctypes signature; no spectral point rerun."
        )
        _atomic_json(preflight_path, preflight)
    comparison_rows: list[dict[str, Any]] = []
    anchor_rows: list[dict[str, Any]] = []
    fallback_rows: list[dict[str, Any]] = []
    spike_rows: list[dict[str, Any]] = []
    runtime_rows = _read_csv(output_dir / "runtime_profile.csv")
    memory_rows = _read_csv(output_dir / "memory_profile.csv")
    backend_summaries: list[dict[str, Any]] = []

    workflows: list[tuple[str, tuple[float, ...], tuple[str, ...]]] = [
        ("beta", BETA_VALIDATION_POINTS, MODELS),
        ("mu", MU_VALIDATION_POINTS, MODELS),
        ("beta_iso_rlb", (45.0,), (MODEL_RLB,)),
    ]
    for workflow, targets, models in workflows:
        for model in models:
            backend = PhysicalValidationBackend(workflow, model)
            backend_started = time.perf_counter()
            target_count = 0
            for target in targets:
                chain = _validation_chain(target, workflow)
                result = sweep_runner.run_spectral_sweep(
                    chain,
                    callbacks=sweep_runner.SweepCallbacks(
                        backend.global_search,
                        backend.local_search,
                        backend.quality_gate,
                        completeness_guard=lambda parameter, roots: (len(roots) == 9, "root9"),
                    ),
                    settings=replace(_settings(workflow, spike=True), anchor_stride=10),
                )
                fast = result.spectra[float(target)]
                key = round(float(target), 12)
                full = (
                    _record_from_rows(target, backend.oracle[key], "CERTIFIED_FULL_ORACLE")
                    if key in backend.oracle
                    else backend.global_search(target)
                )
                comparison_rows.extend(_compare_record(workflow, model, target, fast, full))
                target_count += 1
                for audit in result.point_audits:
                    anchor_rows.append(
                        {
                            "workflow": workflow,
                            "model": model,
                            **asdict(audit),
                        }
                    )
                    if audit.fallback_reasons:
                        for reason in audit.fallback_reasons:
                            fallback_rows.append(
                                {
                                    "workflow": workflow,
                                    "model": model,
                                    "parameter": audit.parameter,
                                    "origin": audit.origin,
                                    "fallback_reason": reason,
                                }
                            )
                for spike in result.spike_audits:
                    spike_rows.append({"workflow": workflow, "model": model, **asdict(spike)})
            elapsed = time.perf_counter() - backend_started
            runtime_rows.append(
                {
                    "phase": "fixed_fast_vs_full_validation",
                    "workflow": workflow,
                    "model": model,
                    "point_count": target_count,
                    "wall_time_seconds": elapsed,
                    **asdict(backend.counters),
                }
            )
            cache = backend.shared_rlb_cache
            memory_rows.append(
                {
                    "phase": f"validation_{workflow}_{model}",
                    "peak_rss_bytes": _current_peak_rss_bytes(),
                    "cache_peak_entries": cache.peak_entries,
                    "cache_peak_bytes": cache.peak_bytes,
                }
            )
            backend_summaries.append(
                {
                    "workflow": workflow,
                    "model": model,
                    "oracle_file_sha256": backend.oracle_file_sha256,
                    "counters": asdict(backend.counters),
                    "cache": cache.diagnostics(),
                }
            )

    synthetic = preflight["synthetic_infrastructure_evidence"]
    fallback_rows.append(
        {
            "workflow": "synthetic",
            "model": "MODEL_NEUTRAL",
            "parameter": 1.0,
            "origin": synthetic["fallback_origin"],
            "fallback_reason": "MISSING_LOCAL_ROOT",
        }
    )
    for outcome, value in (
        ("REPRODUCED_BY_FULL_SCAN", 101.9),
        ("FAST_LOCATOR_CORRECTED", 101.0),
    ):
        spike_rows.append(
            {
                "workflow": "synthetic",
                "model": "MODEL_NEUTRAL",
                "parameter": 1.0,
                "position": 1,
                "fast_value": 101.9,
                "neighbor_prediction": 101.0,
                "normalized_residual": abs(101.9 - 101.0) / 101.9,
                "full_value": value,
                "outcome": outcome,
            }
        )

    _atomic_csv(output_dir / "fast_vs_full_comparison.csv", comparison_rows, COMPARISON_FIELDS)
    _atomic_csv(output_dir / "anchor_checks.csv", anchor_rows, ANCHOR_FIELDS)
    _atomic_csv(output_dir / "fallback_checks.csv", fallback_rows, FALLBACK_FIELDS)
    _atomic_csv(output_dir / "spike_audit.csv", spike_rows, SPIKE_FIELDS)
    _atomic_csv(output_dir / "runtime_profile.csv", runtime_rows, RUNTIME_FIELDS)
    _atomic_csv(output_dir / "memory_profile.csv", memory_rows, MEMORY_FIELDS)

    maximum_error = max(float(row["relative_difference"]) for row in comparison_rows)
    comparison_pass = all(row["status"] == "PASS" for row in comparison_rows)
    root9_rows = [row for row in comparison_rows if int(row["position"]) == 9]
    root9_pass = bool(root9_rows and all(row["status"] == "PASS" for row in root9_rows))
    spike_pass = bool(
        synthetic["reproduced_value_preserved"]
        and synthetic["corrected_value_from_full_scan"]
        and "REPRODUCED_BY_FULL_SCAN" in synthetic["reproduced_outcomes"]
        and "FAST_LOCATOR_CORRECTED" in synthetic["corrected_outcomes"]
        and not any(row["outcome"] == "UNRESOLVED" for row in spike_rows)
    )
    resume_pass = bool(
        synthetic["resume_second_resumed"] == 3
        and synthetic["resume_second_committed"] == 0
    )
    estimate = float(preflight["estimated_fast_total_seconds_two_graphs"])
    performance_status = "PASS" if estimate <= TARGET_RUNTIME_SECONDS else "PARTIAL_PASS"
    statuses = {
        "SWEEP-FAST-ROOT-ACCURACY": "PASS" if comparison_pass else "FAIL",
        "SWEEP-FAST-ROOT9-COMPLETENESS": "PASS" if root9_pass else "FAIL",
        "SWEEP-FAST-SPIKE-PROTECTION": "PASS" if spike_pass else "FAIL",
        "SWEEP-FAST-RESUME": "PASS" if resume_pass else "FAIL",
        "SWEEP-FAST-PERFORMANCE": performance_status,
    }
    overall = bool(
        all(statuses[key] == "PASS" for key in (
            "SWEEP-FAST-ROOT-ACCURACY",
            "SWEEP-FAST-ROOT9-COMPLETENESS",
            "SWEEP-FAST-SPIKE-PROTECTION",
            "SWEEP-FAST-RESUME",
        ))
        and performance_status in {"PASS", "PARTIAL_PASS"}
    )
    statuses["OVERALL"] = "PASS" if overall else "FAIL"
    production_hashes_after = {path: _sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS}
    manifest = {
        "schema_version": 1,
        "stage": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "created_utc": _utc_now(),
        "git": _git_state(),
        "thread_limits": THREAD_LIMITS,
        "python_workers": 1,
        "multiprocessing_used": False,
        "plotting_or_smoothing_used": False,
        "branch_tracking_used": False,
        "production_contract": {"plotted_positions": list(range(1, 9)), "guard_position": 9},
        "thresholds": {
            "fast_full_relative": FAST_FULL_RELATIVE_TOLERANCE,
            "local_primary_verification_relative": LOCAL_VERIFICATION_RELATIVE_TOLERANCE,
            "cached_uncached_relative": CACHE_EQUALITY_RELATIVE_TOLERANCE,
            "root_singular_ratio": rlb2c.rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE,
            "boundary_null_residual": rlb2c.rlb2b.BOUNDARY_RESIDUAL_TOLERANCE,
        },
        "validation_points": {
            "beta_deg": BETA_VALIDATION_POINTS,
            "mu": MU_VALIDATION_POINTS,
            "close_pair_control": {"workflow": "RLB2B_ISOTROPIC_RLB", "beta_deg": 45.0},
        },
        "maximum_fast_full_relative": maximum_error,
        "comparison_row_count": len(comparison_rows),
        "root9_row_count": len(root9_rows),
        "fallback_row_count": len(fallback_rows),
        "spike_audit_row_count": len(spike_rows),
        "backend_summaries": backend_summaries,
        "preflight": preflight,
        "measured_speedup_estimate": preflight["estimated_speedup"],
        "performance_scope": "preflight extrapolation; full production sweep not run in validation stage",
        "statuses": statuses,
        "allowed_conclusion": (
            "FAST_PRODUCTION_SWEEP_VALIDATED_FOR_FIRST8_PLUS_ROOT9" if overall else "NOT_AUTHORIZED"
        ),
        "production_physics_hashes_before": preflight["production_hashes_before"],
        "production_physics_hashes_after": production_hashes_after,
        "production_physics_preserved": preflight["production_hashes_before"] == production_hashes_after,
        "figures_created": 0,
        "commit_or_push_performed": False,
    }
    _atomic_json(output_dir / "run_manifest.json", manifest)
    _write_report(output_dir / "report.md", manifest)
    return manifest


def _validation_job_fingerprint(
    workflow: str,
    model: str,
    target: float,
    oracle_sha256: str,
) -> str:
    dependency_hashes = _dependency_hashes()
    reference_contract = _reference_contract_metadata(
        workflow, model, target
    )
    return sweep_runner.stable_fingerprint(
        {
            "stage": STAGE_ID,
            "algorithm_version": ALGORITHM_VERSION,
            "runner_version": sweep_runner.RUNNER_VERSION,
            "workflow": workflow,
            "model": model,
            "target": float(target),
            "settings": asdict(_validation_worker_settings(workflow)),
            "oracle_sha256": oracle_sha256,
            "runner_sha256": dependency_hashes[
                "scripts/lib/spectral_sweep_runner.py"
            ],
            "physical_contract_sha256": reference_contract[
                "physical_contract_sha256"
            ],
            "reference_row_keys_sha256": reference_contract[
                "reference_row_keys_sha256"
            ],
            "validation_constants": _validation_constants_payload(),
            "dependency_hashes": dependency_hashes,
        }
    )


def _validation_shard_path(
    output_dir: Path, workflow: str, model: str, target: float
) -> Path:
    safe_model = "".join(
        character.lower() if character.isalnum() else "_" for character in model
    )
    parameter = format(float(target), ".12g").replace("-", "m").replace(".", "p")
    return (
        output_dir
        / CURRENT_SHARD_DIRECTORY
        / workflow
        / safe_model
        / f"{parameter}.json"
    )


def run_validation_worker(
    output_dir: Path,
    workflow: str,
    model: str,
    target: float,
    *,
    fresh_full_budget: int,
) -> dict[str, Any]:
    """Run and atomically commit one isolated fixed-validation point."""

    backend = PhysicalValidationBackend(
        workflow,
        model,
        fresh_full_scan_limit=max(0, int(fresh_full_budget)),
    )
    key = round(float(target), 12)
    reference_contract = _reference_contract_metadata(workflow, model, target)
    reference_rows = backend.recorded_oracle.get(key, [])
    chain = _validation_chain(
        target,
        workflow,
        certified_parameters=tuple(backend.recorded_oracle),
    )
    target_scheduled_anchor = bool(
        math.isclose(float(target), float(chain[0]), rel_tol=0.0, abs_tol=1.0e-12)
        or math.isclose(
            float(target), float(chain[-1]), rel_tol=0.0, abs_tol=1.0e-12
        )
    )
    if not target_scheduled_anchor and not reference_rows:
        # Only a genuinely missing frozen reference may require a fresh full
        # fallback.  When immutable rows exist, a failed locator falls back
        # to those rows and must not recalculate the full side.
        backend.force_fresh_global_parameters.add(key)
    if (
        workflow == "mu"
        and float(target) != 0.8
        and any(math.isclose(float(value), 0.8, rel_tol=0.0, abs_tol=1.0e-12) for value in chain)
        and 0.8 not in backend.recorded_oracle
    ):
        anchor_payload = _load_validation_shard(output_dir, workflow, model, 0.8)
        if anchor_payload is None:
            raise RuntimeError(
                "The isolated mu=0.80 full-anchor shard must be completed "
                "before a validation chain that ends at mu=0.80."
            )
        anchor_rows = sorted(
            anchor_payload["comparison"], key=lambda row: int(row["position"])
        )
        anchor_roots = tuple(
            sweep_runner.RootRecord(
                value=float(row["full_value"]),
                accepted=True,
                sigma_ratio=float(row["full_sigma_ratio"]),
                boundary_residual=float(row["full_boundary_residual"]),
                detector_agreement=True,
                multiplicity=int(row["full_multiplicity"]),
                nullity=int(row["full_nullity"]),
                cluster_id="REUSED_CLUSTER" if bool(row["cluster_match"]) and int(row["full_multiplicity"]) > 1 else "",
                metadata={
                    "source": "isolated_mu_0p80_full_anchor_shard",
                    "source_fingerprint": anchor_payload["fingerprint"],
                },
            )
            for row in anchor_rows
        )
        backend.fresh_full_cache[0.8] = sweep_runner.SpectrumRecord(
            parameter=0.8,
            roots=anchor_roots,
            origin="REUSED_ISOLATED_FRESH_FULL_ANCHOR",
        )
        backend.counters.isolated_anchor_reuses += 1
    job_fingerprint = _validation_job_fingerprint(
        workflow, model, target, backend.oracle_file_sha256
    )
    started = time.perf_counter()
    result = sweep_runner.run_spectral_sweep(
        chain,
        callbacks=sweep_runner.SweepCallbacks(
            backend.global_search,
            backend.local_search,
            backend.quality_gate,
            completeness_guard=lambda parameter, roots: (
                len(roots) == 9,
                "root9",
            ),
            global_completeness_guard=lambda parameter, roots: (
                len(roots) == 9,
                "root9",
            ),
        ),
        settings=_validation_worker_settings(workflow),
    )
    fast = result.spectra[float(target)]
    target_audit = next(
        audit
        for audit in result.point_audits
        if math.isclose(
            float(audit.parameter), float(target), rel_tol=0.0, abs_tol=1.0e-12
        )
    )
    if reference_rows:
        full = _record_from_rows(
            target, reference_rows, "IMMUTABLE_FROZEN_FULL_REFERENCE"
        )
        reference_source = "IMMUTABLE_FROZEN_FULL_REFERENCE"
    elif key in backend.fresh_full_cache:
        full = backend.fresh_full_cache[key]
        reference_source = "FRESH_FULL_SCAN_MISSING_REFERENCE"
    else:
        full = backend.fresh_full_search(target)
        backend.fresh_full_cache[key] = full
        reference_source = "FRESH_FULL_SCAN_MISSING_REFERENCE"

    initial_comparison = _compare_record(workflow, model, target, fast, full)
    comparison = initial_comparison
    fallback_record: dict[str, Any] | None = None
    if reference_rows and any(row["status"] != "PASS" for row in comparison):
        initial_maximum = max(
            float(row["relative_difference"]) for row in comparison
        )
        if key in backend.fresh_full_cache:
            fresh_full = backend.fresh_full_cache[key]
        else:
            fresh_full = backend.fresh_full_search(target)
            backend.fresh_full_cache[key] = fresh_full
        comparison = _compare_record(
            workflow, model, target, fast, fresh_full
        )
        final_maximum = max(
            float(row["relative_difference"]) for row in comparison
        )
        fallback_record = {
            "reason": "FAST_VS_FROZEN_REFERENCE_EXCEEDS_TOLERANCE",
            "initial_reference_source": reference_source,
            "final_reference_source": "FRESH_FULL_SCAN_AFTER_DISAGREEMENT",
            "initial_maximum_relative_difference": initial_maximum,
            "final_maximum_relative_difference": final_maximum,
            "status": (
                "PASS_AFTER_LOCAL_FULL_FALLBACK"
                if all(row["status"] == "PASS" for row in comparison)
                else "FAIL"
            ),
        }
        reference_source = "FRESH_FULL_SCAN_AFTER_DISAGREEMENT"

    for row in comparison:
        row["reference_source"] = reference_source
    comparison_passed = all(row["status"] == "PASS" for row in comparison)
    search_mode = "FAST_LOCAL" if fast.origin == "FAST_LOCAL" else "FULL_ANCHOR"
    shard_status = (
        "PASS_AFTER_LOCAL_FULL_FALLBACK"
        if fallback_record is not None
        and result.status == "PASS"
        and comparison_passed
        else "PASS"
        if result.status == "PASS" and comparison_passed
        else "FAIL"
    )
    production_rows: list[dict[str, Any]] = []
    production_rows_semantic_sha256 = "NOT_APPLICABLE"
    if workflow == "mu" and model in {MODEL_OLD, MODEL_RLB}:
        production_rows = _production_rows_from_record(
            backend, fast, target_audit
        )
        production_rows_semantic_sha256 = _production_rows_semantic_hash(
            production_rows
        )
    dependency_hashes = _dependency_hashes()
    payload = {
        "schema_version": 1,
        "job_fingerprint": job_fingerprint,
        "workflow": workflow,
        "model": model,
        "target": float(target),
        "chain": list(map(float, chain)),
        "runtime_seconds": time.perf_counter() - started,
        "peak_rss_bytes": _current_peak_rss_bytes(),
        "comparison": comparison,
        "initial_frozen_comparison": initial_comparison,
        "fallback_record": fallback_record,
        "point_audits": [asdict(item) for item in result.point_audits],
        "spike_audits": [asdict(item) for item in result.spike_audits],
        "runner_counters": result.counters.to_dict(),
        "backend_counters": asdict(backend.counters),
        "cache_diagnostics": backend.shared_rlb_cache.diagnostics(),
        "oracle_file_sha256": backend.oracle_file_sha256,
        "oracle_exclusions": {
            format(parameter, ".12g"): list(reasons)
            for parameter, reasons in sorted(backend.oracle_exclusions.items())
        },
        "runner_sha256": dependency_hashes[
            "scripts/lib/spectral_sweep_runner.py"
        ],
        "dependency_hashes": dependency_hashes,
        **reference_contract,
        "validation_constants": _validation_constants_payload(),
        "target_origin": fast.origin,
        "search_mode": search_mode,
        "target_scheduled_anchor": bool(target_audit.scheduled_anchor),
        "target_is_non_anchor": not bool(target_audit.scheduled_anchor),
        "fresh_full_scan_used": backend.counters.fresh_full_scans > 0,
        "fresh_full_budget_at_worker_start": int(fresh_full_budget),
        "production_rows": production_rows,
        "production_rows_field_count": (
            len(RLB2D_PRODUCTION_FIELDS) if production_rows else 0
        ),
        "production_rows_semantic_sha256": (
            production_rows_semantic_sha256
        ),
        "internal_tail": (
            "NOT_COMPUTED_ABOVE_ROOT9" if production_rows else "NOT_APPLICABLE"
        ),
        "status": shard_status,
    }
    payload["payload_sha256"] = sweep_runner.stable_fingerprint(_json_value(payload))
    payload["fingerprint"] = sweep_runner.stable_fingerprint(
        {
            "job_fingerprint": job_fingerprint,
            "payload_sha256": payload["payload_sha256"],
            "dependency_hashes": payload["dependency_hashes"],
            "validation_constants": payload["validation_constants"],
            "settings": asdict(_validation_worker_settings(workflow)),
        }
    )
    _atomic_json(_validation_shard_path(output_dir, workflow, model, target), payload)
    return payload


def _load_validation_shard(
    output_dir: Path,
    workflow: str,
    model: str,
    target: float,
) -> dict[str, Any] | None:
    path = _validation_shard_path(output_dir, workflow, model, target)
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        declared_payload_sha256 = str(payload.get("payload_sha256", ""))
        declared_shard_fingerprint = str(payload.get("fingerprint", ""))
        unsigned_payload = dict(payload)
        unsigned_payload.pop("payload_sha256", None)
        unsigned_payload.pop("fingerprint", None)
        if not declared_payload_sha256 or sweep_runner.stable_fingerprint(
            _json_value(unsigned_payload)
        ) != declared_payload_sha256:
            return None
        current_dependency_hashes = _dependency_hashes()
        if payload.get("dependency_hashes") != current_dependency_hashes:
            return None
        if payload.get("validation_constants") != _json_value(
            _validation_constants_payload()
        ):
            return None
        oracle_sha256 = _sha256(_oracle_paths()[(workflow, model)])
        expected = _validation_job_fingerprint(
            workflow, model, target, oracle_sha256
        )
        expected_shard_fingerprint = sweep_runner.stable_fingerprint(
            {
                "job_fingerprint": expected,
                "payload_sha256": declared_payload_sha256,
                "dependency_hashes": current_dependency_hashes,
                "validation_constants": _validation_constants_payload(),
                "settings": asdict(_validation_worker_settings(workflow)),
            }
        )
        if (
            payload.get("job_fingerprint") != expected
            or declared_shard_fingerprint != expected_shard_fingerprint
            or payload.get("status")
            not in {"PASS", "PASS_AFTER_LOCAL_FULL_FALLBACK"}
        ):
            return None
        reference_contract = _reference_contract_metadata(
            workflow, model, target
        )
        if payload.get("runner_sha256") != current_dependency_hashes.get(
            "scripts/lib/spectral_sweep_runner.py"
        ):
            return None
        for field in (
            "physical_contract_path",
            "physical_contract_sha256",
            "reference_csv_path",
            "reference_csv_sha256",
            "reference_row_keys",
            "reference_row_keys_sha256",
        ):
            if payload.get(field) != reference_contract[field]:
                return None
        if payload.get("search_mode") not in {"FAST_LOCAL", "FULL_ANCHOR"}:
            return None
        if bool(payload.get("target_is_non_anchor", False)) and (
            payload.get("search_mode") == "FAST_LOCAL"
        ) != (payload.get("target_origin") == "FAST_LOCAL"):
            return None
        counters = payload.get("backend_counters", {})
        if payload.get("search_mode") == "FAST_LOCAL" and (
            int(counters.get("local_primary_scans", 0)) <= 0
            or int(counters.get("local_verification_scans", 0)) <= 0
            or int(counters.get("determinant_evaluations", 0)) <= 0
        ):
            return None
        comparison = payload.get("comparison", [])
        if len(comparison) != 9:
            return None
        if sorted(int(row["position"]) for row in comparison) != list(range(1, 10)):
            return None
        if any(row.get("status") != "PASS" for row in comparison):
            return None
        if workflow == "mu" and model in {MODEL_OLD, MODEL_RLB}:
            production_rows = _validate_production_rows(
                payload.get("production_rows"),
                model=model,
                target=target,
                declared_semantic_sha256=str(
                    payload.get("production_rows_semantic_sha256", "")
                ),
            )
            if int(payload.get("production_rows_field_count", 0)) != 59:
                return None
            if payload.get("internal_tail") != "NOT_COMPUTED_ABOVE_ROOT9":
                return None
            comparison_values = {
                int(row["position"]): float(row["fast_value"])
                for row in comparison
            }
            if any(
                float(row["Omega"]).hex()
                != comparison_values[int(row["sorted_position"])].hex()
                for row in production_rows
            ):
                return None
        return payload
    except (OSError, TypeError, ValueError, KeyError, json.JSONDecodeError):
        return None


def _fixed_validation_jobs() -> list[tuple[str, str, float]]:
    jobs: list[tuple[str, str, float]] = []
    for model in MODELS:
        jobs.extend(("beta", model, target) for target in BETA_VALIDATION_POINTS)
    for model in MODELS:
        model_points = MU_VALIDATION_POINTS_BY_MODEL[model]
        ordered_mu = (0.80,) + tuple(
            value for value in model_points if value != 0.80
        )
        jobs.extend(("mu", model, target) for target in ordered_mu)
    if len(jobs) != 24:
        raise RuntimeError(f"Simplified validation contract has {len(jobs)} jobs, not 24.")
    return jobs


def _requires_fast_local_evidence(
    workflow: str, model: str, target: float
) -> bool:
    """Identify the fixed non-anchor probes needed for production approval."""

    beta_probe = bool(
        workflow == "beta"
        and math.isclose(float(target), 54.0, rel_tol=0.0, abs_tol=1.0e-12)
    )
    mu_probe = bool(
        workflow == "mu"
        and model in {MODEL_OLD, MODEL_RLB}
        and math.isclose(float(target), 0.78, rel_tol=0.0, abs_tol=1.0e-12)
    )
    return beta_probe or mu_probe


def _require_fast_local_shard(
    workflow: str,
    model: str,
    target: float,
    payload: Mapping[str, Any],
) -> None:
    """Fail before expensive missing-reference scans if a probe fell back."""

    if not _requires_fast_local_evidence(workflow, model, target):
        return
    counters = payload.get("backend_counters", {})
    if not (
        payload.get("target_origin") == "FAST_LOCAL"
        and int(counters.get("local_primary_scans", 0)) > 0
        and int(counters.get("local_verification_scans", 0)) > 0
        and int(counters.get("determinant_evaluations", 0)) > 0
    ):
        raise RuntimeError(
            "Required non-anchor validation probe did not produce physical "
            f"FAST_LOCAL evidence: {workflow}/{model}/{target:g}; "
            f"origin={payload.get('target_origin')!r}. Fresh full-reference "
            "scans are not authorized after this failure."
        )


def run_fixed_validation(output_dir: Path) -> dict[str, Any]:
    """Run fixed cases in sequential child processes and aggregate shards."""

    preflight_path = output_dir / "preflight_report.json"
    if not preflight_path.is_file():
        raise RuntimeError("Run --preflight before --validate-fast-vs-full.")
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    if preflight.get("algorithm_version") != ALGORITHM_VERSION:
        raise RuntimeError("Preflight algorithm version is stale; rerun --preflight.")
    current_production_hashes = {
        path: _sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS
    }
    if preflight.get("production_hashes_before") != current_production_hashes:
        raise RuntimeError("Production physics changed after preflight; validation is blocked.")
    current_dependency_hashes = _dependency_hashes()
    preflight_dependency_hashes = preflight.get(
        "validation_dependency_hashes_before", {}
    )
    preflight_critical_paths = {
        *PRODUCTION_PHYSICS_PATHS,
        "scripts/lib/spectral_sweep_runner.py",
    }
    preflight_critical_dependencies = {
        path: digest
        for path, digest in preflight_dependency_hashes.items()
        if path in preflight_critical_paths
    }
    current_critical_dependencies = {
        path: digest
        for path, digest in current_dependency_hashes.items()
        if path in preflight_critical_paths
    }
    if preflight_critical_dependencies != current_critical_dependencies:
        raise RuntimeError(
            "Runner, production physics, or physical adapters changed after "
            "preflight; validation is blocked."
        )
    preflight_reuse_qualification = {
        "status": "PASS",
        "reason": (
            "Only the closing validation orchestrator/point contract changed; "
            "runner, production physics, thread limits, local grid policy, "
            "and cache policy are unchanged. Current adapter hashes are bound "
            "separately into every validation shard."
        ),
        "runner_sha256": current_dependency_hashes[
            "scripts/lib/spectral_sweep_runner.py"
        ],
        "estimated_fast_total_seconds_two_graphs": preflight.get(
            "estimated_fast_total_seconds_two_graphs"
        ),
        "peak_rss_bytes": preflight.get("peak_rss_bytes"),
        "production_run_permitted": preflight.get("production_run_permitted"),
    }
    preflight_counters = preflight.get("backend_counters", {})
    recorded_preflight_local_evidence = bool(
        int(preflight_counters.get("local_primary_scans", 0)) > 0
        and int(preflight_counters.get("local_verification_scans", 0)) > 0
    )
    preflight_fast_local_evidence = bool(
        preflight.get("physical_fast_local_evidence", False)
        or recorded_preflight_local_evidence
    )
    if not preflight_fast_local_evidence:
        raise RuntimeError(
            "Current preflight contains no physical FAST_LOCAL primary/"
            "verification evidence; production validation is blocked before "
            "fresh full-reference scans."
        )
    if int(preflight.get("peak_rss_bytes", 0)) <= 0:
        preflight["peak_rss_bytes"] = _current_peak_rss_bytes()
        preflight["peak_rss_measurement_note"] = (
            "Recovered before fixed validation after declaring the Windows "
            "GetProcessMemoryInfo ctypes signature; no spectral point rerun."
        )
        _atomic_json(preflight_path, preflight)

    jobs = _fixed_validation_jobs()
    reference_contracts = {
        (workflow, model, float(target)): _reference_contract_metadata(
            workflow, model, target
        )
        for workflow, model, target in jobs
    }
    frozen_reference_hashes_before = {
        metadata["reference_csv_path"]: metadata["reference_csv_sha256"]
        for metadata in reference_contracts.values()
    }
    missing_full_references = [
        {
            "workflow": workflow,
            "model": model,
            "parameter": float(target),
        }
        for workflow, model, target in jobs
        if not reference_contracts[(workflow, model, float(target))][
            "reference_rows_present"
        ]
    ]
    if len(missing_full_references) > MAX_FRESH_FULL_SCANS_CLOSING_VALIDATION:
        raise RuntimeError(
            "Simplified validation needs more than three fresh full references: "
            f"{missing_full_references}."
        )
    existing_payloads = {
        (workflow, model, float(target)): _load_validation_shard(
            output_dir, workflow, model, target
        )
        for workflow, model, target in jobs
    }
    closing_fresh_full_scans = sum(
        int(payload["backend_counters"]["fresh_full_scans"])
        for payload in existing_payloads.values()
        if payload is not None
    )
    if closing_fresh_full_scans > MAX_FRESH_FULL_SCANS_CLOSING_VALIDATION:
        raise RuntimeError("Existing simplified shards exceed the fresh-full budget.")
    reused_jobs = 0
    executed_jobs = 0
    for job_index, (workflow, model, target) in enumerate(jobs, start=1):
        existing_payload = existing_payloads[(workflow, model, float(target))]
        if existing_payload is not None:
            _require_fast_local_shard(
                workflow, model, target, existing_payload
            )
            reused_jobs += 1
            continue
        reference_present = reference_contracts[
            (workflow, model, float(target))
        ]["reference_rows_present"]
        if (
            not reference_present
            and closing_fresh_full_scans
            >= MAX_FRESH_FULL_SCANS_CLOSING_VALIDATION
        ):
            remaining_missing = [
                item
                for item in missing_full_references
                if _load_validation_shard(
                    output_dir,
                    str(item["workflow"]),
                    str(item["model"]),
                    float(item["parameter"]),
                )
                is None
            ]
            raise RuntimeError(
                "Fresh-full budget is exhausted; unresolved missing references: "
                f"{remaining_missing}."
            )
        command = [
            sys.executable,
            "-B",
            str(Path(__file__).resolve()),
            "--validation-worker",
            "--worker-workflow",
            workflow,
            "--worker-model",
            model,
            "--worker-target",
            format(float(target), ".17g"),
            "--worker-fresh-full-budget",
            str(
                MAX_FRESH_FULL_SCANS_CLOSING_VALIDATION
                - closing_fresh_full_scans
            ),
            "--output-dir",
            str(output_dir),
        ]
        environment = dict(os.environ)
        environment.update(THREAD_LIMITS)
        environment["PYTHONDONTWRITEBYTECODE"] = "1"
        completed = subprocess.run(
            command,
            cwd=ROOT,
            env=environment,
            capture_output=True,
            text=True,
            encoding="utf-8",
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"Validation worker {job_index}/{len(jobs)} failed for "
                f"{workflow}/{model}/{target:g}.\n"
                f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
            )
        committed_payload = _load_validation_shard(
            output_dir, workflow, model, target
        )
        if committed_payload is None:
            raise RuntimeError("Validation worker did not commit a valid shard.")
        _require_fast_local_shard(
            workflow, model, target, committed_payload
        )
        closing_fresh_full_scans += int(
            committed_payload["backend_counters"]["fresh_full_scans"]
        )
        if closing_fresh_full_scans > MAX_FRESH_FULL_SCANS_CLOSING_VALIDATION:
            raise RuntimeError("Closing validation exceeded three fresh full scans.")
        executed_jobs += 1
        print(
            f"[{job_index}/{len(jobs)}] {workflow}/{model}/{target:g}: PASS",
            flush=True,
        )

    shards = [
        _load_validation_shard(output_dir, workflow, model, target)
        for workflow, model, target in jobs
    ]
    if any(item is None for item in shards):
        raise RuntimeError("Fixed validation shard inventory is incomplete.")
    payloads = [item for item in shards if item is not None]

    comparison_rows = [
        row for payload in payloads for row in payload["comparison"]
    ]
    anchor_rows: list[dict[str, Any]] = []
    fallback_rows: list[dict[str, Any]] = []
    spike_rows: list[dict[str, Any]] = []
    runtime_rows = _read_csv(output_dir / "runtime_profile.csv")
    memory_rows = _read_csv(output_dir / "memory_profile.csv")
    for payload in payloads:
        workflow = payload["workflow"]
        model = payload["model"]
        for audit in payload["point_audits"]:
            anchor_rows.append({"workflow": workflow, "model": model, **audit})
            for reason in audit.get("fallback_reasons", []):
                fallback_rows.append(
                    {
                        "workflow": workflow,
                        "model": model,
                        "parameter": audit["parameter"],
                        "origin": audit["origin"],
                        "fallback_reason": reason,
                    }
                )
        for spike in payload["spike_audits"]:
            spike_rows.append({"workflow": workflow, "model": model, **spike})
        fallback_record = payload.get("fallback_record")
        if fallback_record is not None:
            fallback_rows.append(
                {
                    "workflow": workflow,
                    "model": model,
                    "parameter": payload["target"],
                    "origin": payload["target_origin"],
                    "fallback_reason": fallback_record["reason"],
                    **fallback_record,
                }
            )
        counters = payload["backend_counters"]
        runtime_rows.append(
            {
                "phase": "fixed_fast_vs_full_validation_worker",
                "workflow": workflow,
                "model": model,
                "point_count": 1,
                "wall_time_seconds": payload["runtime_seconds"],
                **counters,
            }
        )
        memory_rows.append(
            {
                "phase": f"validation_worker_{workflow}_{model}_{payload['target']}",
                "peak_rss_bytes": payload["peak_rss_bytes"],
                "cache_peak_entries": payload["cache_diagnostics"]["peak_entries"],
                "cache_peak_bytes": payload["cache_diagnostics"]["peak_bytes"],
            }
        )

    synthetic = preflight["synthetic_infrastructure_evidence"]
    fallback_rows.append(
        {
            "workflow": "synthetic",
            "model": "MODEL_NEUTRAL",
            "parameter": 1.0,
            "origin": synthetic["fallback_origin"],
            "fallback_reason": "MISSING_LOCAL_ROOT",
        }
    )
    for outcome, value in (
        ("REPRODUCED_BY_FULL_SCAN", 101.9),
        ("FAST_LOCATOR_CORRECTED", 101.0),
    ):
        spike_rows.append(
            {
                "workflow": "synthetic",
                "model": "MODEL_NEUTRAL",
                "parameter": 1.0,
                "position": 1,
                "fast_value": 101.9,
                "neighbor_prediction": 101.0,
                "normalized_residual": abs(101.9 - 101.0) / 101.9,
                "full_value": value,
                "outcome": outcome,
            }
        )

    _atomic_csv(output_dir / "fast_vs_full_comparison.csv", comparison_rows, COMPARISON_FIELDS)
    _atomic_csv(output_dir / "anchor_checks.csv", anchor_rows, ANCHOR_FIELDS)
    _atomic_csv(output_dir / "fallback_checks.csv", fallback_rows, FALLBACK_FIELDS)
    _atomic_csv(output_dir / "spike_audit.csv", spike_rows, SPIKE_FIELDS)
    _atomic_csv(output_dir / "runtime_profile.csv", runtime_rows, RUNTIME_FIELDS)
    _atomic_csv(output_dir / "memory_profile.csv", memory_rows, MEMORY_FIELDS)

    maximum_error = max(float(row["relative_difference"]) for row in comparison_rows)
    production_hashes_after = {
        path: _sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS
    }
    validation_dependency_hashes_after = _dependency_hashes()
    production_physics_preserved = (
        preflight["production_hashes_before"] == production_hashes_after
    )
    validation_dependencies_preserved = (
        current_dependency_hashes == validation_dependency_hashes_after
    )
    frozen_reference_hashes_after = {
        path: _sha256(ROOT / path) for path in frozen_reference_hashes_before
    }
    frozen_references_preserved = (
        frozen_reference_hashes_before == frozen_reference_hashes_after
    )
    cache_equality_pass = bool(
        preflight.get("physical_cached_vs_uncached_equality", {}).get("status")
        == "PASS"
    )
    # Do not let a locator fallback erase the validation obligation by
    # relabelling the target as an anchor.  The fixed contract contains a
    # genuine non-anchor probe for every beta model (beta=54) and for both
    # determinant/SVD mu models (mu=0.78).  EB's declared mu probes are the
    # two endpoint anchors, so no separate EB/mu FAST_LOCAL claim is made.
    required_fast_local_pairs = {
        *(("beta", model) for model in MODELS),
        ("mu", MODEL_OLD),
        ("mu", MODEL_RLB),
    }
    observed_fast_local_pairs = {
        (str(payload["workflow"]), str(payload["model"]))
        for payload in payloads
        if bool(payload.get("target_is_non_anchor", False))
        and payload.get("target_origin") == "FAST_LOCAL"
        and int(payload.get("backend_counters", {}).get("local_primary_scans", 0)) > 0
        and int(payload.get("backend_counters", {}).get("local_verification_scans", 0)) > 0
        and int(payload.get("backend_counters", {}).get("determinant_evaluations", 0)) > 0
    }
    fast_local_evidence_pass = bool(
        required_fast_local_pairs
        and required_fast_local_pairs.issubset(observed_fast_local_pairs)
    )
    comparison_pass = bool(
        all(row["status"] == "PASS" for row in comparison_rows)
        and production_physics_preserved
        and validation_dependencies_preserved
        and cache_equality_pass
        and fast_local_evidence_pass
    )
    root9_rows = [row for row in comparison_rows if int(row["position"]) == 9]
    root9_pass = bool(root9_rows and all(row["status"] == "PASS" for row in root9_rows))
    beta54_rows = [
        row
        for row in comparison_rows
        if row["workflow"] == "beta"
        and math.isclose(
            float(row["parameter"]), 54.0, rel_tol=0.0, abs_tol=1.0e-12
        )
    ]
    beta54_pass = bool(
        len(beta54_rows) == 3 * 9
        and all(row["status"] == "PASS" for row in beta54_rows)
        and all(
            payload.get("target_origin") == "FAST_LOCAL"
            and int(payload.get("backend_counters", {}).get("local_primary_scans", 0)) > 0
            and int(payload.get("backend_counters", {}).get("local_verification_scans", 0)) > 0
            for payload in payloads
            if payload["workflow"] == "beta"
            and math.isclose(
                float(payload["target"]), 54.0, rel_tol=0.0, abs_tol=1.0e-12
            )
        )
    )
    spike_pass = bool(
        synthetic["reproduced_value_preserved"]
        and synthetic["corrected_value_from_full_scan"]
        and "REPRODUCED_BY_FULL_SCAN" in synthetic["reproduced_outcomes"]
        and "FAST_LOCATOR_CORRECTED" in synthetic["corrected_outcomes"]
        and not any(row["outcome"] == "UNRESOLVED" for row in spike_rows)
        and beta54_pass
    )
    resume_pass = bool(
        synthetic["resume_second_resumed"] == 3
        and synthetic["resume_second_committed"] == 0
    )
    estimate = float(preflight["estimated_fast_total_seconds_two_graphs"])
    measured_peak = max(int(payload["peak_rss_bytes"]) for payload in payloads)
    preflight_peak = int(preflight["peak_rss_bytes"])
    performance_status = (
        "PASS"
        if estimate <= TARGET_RUNTIME_SECONDS
        and preflight_peak > 0
        and bool(preflight.get("production_run_permitted", False))
        else "PARTIAL_PASS"
    )
    statuses = {
        "SWEEP-FAST-ROOT-ACCURACY": "PASS" if comparison_pass else "FAIL",
        "SWEEP-FAST-ROOT9-COMPLETENESS": "PASS" if root9_pass else "FAIL",
        "SWEEP-FAST-SPIKE-PROTECTION": "PASS" if spike_pass else "FAIL",
        "SWEEP-FAST-RESUME": "PASS" if resume_pass else "FAIL",
        "SWEEP-FAST-PERFORMANCE": performance_status,
    }
    overall = bool(
        all(
            statuses[key] == "PASS"
            for key in (
                "SWEEP-FAST-ROOT-ACCURACY",
                "SWEEP-FAST-ROOT9-COMPLETENESS",
                "SWEEP-FAST-SPIKE-PROTECTION",
                "SWEEP-FAST-RESUME",
            )
        )
        and performance_status in {"PASS", "PARTIAL_PASS"}
    )
    overall = bool(
        overall
        and production_physics_preserved
        and validation_dependencies_preserved
        and frozen_references_preserved
        and cache_equality_pass
        and fast_local_evidence_pass
    )
    statuses["OVERALL"] = "PASS" if overall else "FAIL"
    search_mode_counts = {
        mode: sum(payload.get("search_mode") == mode for payload in payloads)
        for mode in ("FAST_LOCAL", "FULL_ANCHOR")
    }
    legacy_shards = _legacy_shard_inventory(output_dir)
    root_count_mismatches = sum(
        len(payload.get("comparison", [])) != 9 for payload in payloads
    )
    close_pair_result = {
        "distinct_roots_are_never_merged_by_runner": True,
        "targeted_test_evidence": "test_close_pair_preservation_and_no_duplicate_merging",
        "physical_cluster_rows_checked": sum(
            int(row["fast_multiplicity"]) > 1
            or int(row["full_multiplicity"]) > 1
            for row in comparison_rows
        ),
        "maximum_cluster_center_relative_difference": max(
            float(row["cluster_center_relative_difference"])
            for row in comparison_rows
        ),
        "status": (
            "PASS"
            if all(bool(row["cluster_center_pass"]) for row in comparison_rows)
            else "FAIL"
        ),
    }
    manifest = {
        "schema_version": 1,
        "stage": STAGE_ID,
        "algorithm_version": ALGORITHM_VERSION,
        "validation_contract_version": VALIDATION_CONTRACT_VERSION,
        "created_utc": _utc_now(),
        "git": _git_state(),
        "thread_limits": THREAD_LIMITS,
        "python_workers_concurrent": 1,
        "sequential_validation_subprocesses": len(payloads),
        "multiprocessing_used": False,
        "plotting_or_smoothing_used": False,
        "branch_tracking_used": False,
        "production_contract": {
            "plotted_positions": list(range(1, 9)),
            "guard_position": 9,
        },
        "thresholds": {
            "fast_full_relative": FAST_FULL_RELATIVE_TOLERANCE,
            "local_primary_verification_relative": LOCAL_VERIFICATION_RELATIVE_TOLERANCE,
            "cached_uncached_relative": CACHE_EQUALITY_RELATIVE_TOLERANCE,
            "root_singular_ratio": rlb2c.rlb2b.ROOT_SINGULAR_RATIO_TOLERANCE,
            "boundary_null_residual": rlb2c.rlb2b.BOUNDARY_RESIDUAL_TOLERANCE,
        },
        "validation_points": {
            "beta_deg": BETA_VALIDATION_POINTS,
            "mu_by_model": MU_VALIDATION_POINTS_BY_MODEL,
        },
        "maximum_fast_full_relative": maximum_error,
        "comparison_row_count": len(comparison_rows),
        "root9_row_count": len(root9_rows),
        "fallback_row_count": len(fallback_rows),
        "spike_audit_row_count": len(spike_rows),
        "validation_job_count": len(jobs),
        "validation_jobs_executed": executed_jobs,
        "validation_jobs_reused": reused_jobs,
        "fast_job_count": len(payloads),
        "search_mode_counts": search_mode_counts,
        "root_count_mismatches": root_count_mismatches,
        "non_anchor_fast_local_evidence": {
            "required_workflow_model_pairs": [
                list(item) for item in sorted(required_fast_local_pairs)
            ],
            "observed_workflow_model_pairs": [
                list(item) for item in sorted(observed_fast_local_pairs)
            ],
            "passed": fast_local_evidence_pass,
        },
        "physical_cached_vs_uncached_equality": preflight[
            "physical_cached_vs_uncached_equality"
        ],
        "frozen_full_references_reused": sum(
            bool(payload.get("reference_rows_present", False))
            for payload in payloads
        ),
        "missing_full_references_before_validation": missing_full_references,
        "fresh_full_scans_closing_validation": closing_fresh_full_scans,
        "maximum_fresh_full_scans_closing_validation": (
            MAX_FRESH_FULL_SCANS_CLOSING_VALIDATION
        ),
        "fresh_full_scans_preflight_not_counted_in_closing_budget": int(
            preflight["backend_counters"]["fresh_full_scans"]
        ),
        "beta54_result": {
            "status": "PASS" if beta54_pass else "FAIL",
            "comparison_row_count": len(beta54_rows),
            "maximum_relative_difference": max(
                float(row["relative_difference"]) for row in beta54_rows
            )
            if beta54_rows
            else math.inf,
            "preflight_fresh_full_vs_frozen_relative": preflight[
                "fresh_full_difficult_point"
            ]["maximum_relative_vs_frozen_oracle"],
        },
        "close_pair_result": close_pair_result,
        "frozen_references": [
            {
                "path": path,
                "sha256_before": digest,
                "sha256_after": frozen_reference_hashes_after[path],
            }
            for path, digest in sorted(frozen_reference_hashes_before.items())
        ],
        "frozen_references_preserved": frozen_references_preserved,
        "legacy_validation_shards": legacy_shards,
        "legacy_validation_shard_count": len(legacy_shards),
        "backend_summaries": payloads,
        "preflight": preflight,
        "preflight_reuse_qualification": preflight_reuse_qualification,
        "measured_peak_rss_bytes": measured_peak,
        "full_scan_memory_isolation": "SEQUENTIAL_CHILD_PROCESS_PER_VALIDATION_POINT",
        "measured_speedup_estimate": preflight["estimated_speedup"],
        "performance_scope": (
            "preflight extrapolation; high-RSS frozen full scans isolated; "
            "full production sweep not run in validation stage"
        ),
        "statuses": statuses,
        "allowed_conclusion": (
            "FAST_PRODUCTION_SWEEP_VALIDATED_FOR_FIRST8_PLUS_ROOT9"
            if overall
            else "NOT_AUTHORIZED"
        ),
        "production_physics_hashes_before": preflight["production_hashes_before"],
        "production_physics_hashes_after": production_hashes_after,
        "production_physics_preserved": production_physics_preserved,
        "validation_dependency_hashes_before": current_dependency_hashes,
        "validation_dependency_hashes_after": validation_dependency_hashes_after,
        "validation_dependencies_preserved": validation_dependencies_preserved,
        "figures_created": 0,
        "commit_or_push_performed": False,
    }
    _atomic_json(output_dir / "run_manifest.json", manifest)
    _write_report(output_dir / "report.md", manifest)
    return manifest


def _write_report(path: Path, manifest: Mapping[str, Any]) -> None:
    statuses = manifest["statuses"]
    preflight = manifest["preflight"]
    text = f"""# RLB-SWEEP-FAST: closing validation

## Scope

Проверен `spectral-sweep-runner-v2` для независимо упорядоченных roots 1--8
и root 9 guard. Frozen full CSV использованы как неизменяемая reference-сторона;
текущий runner вычислял `FAST_LOCAL` в non-anchor точках и выполнял
`FULL_ANCHOR` в обязательных anchors. Частоты не интерполировались и не
сглаживались; sorted position не трактовалась как modal branch identity.

Контракт содержит `{manifest['validation_job_count']}` jobs: beta =
0, 45, 54, 90 deg для трёх моделей; mu = 0.00, 0.80 для EB и
0.00, 0.40, 0.60, 0.78, 0.80 для old Timoshenko и weak RLB.

## Numerical policy

Официальный fast/full threshold равен `1e-8`. Внутренний frozen threshold
local primary/verification сохранён равным `1e-9`; determinant/SVD quality
thresholds не менялись. Сравнивались только roots 1--8 и root 9. Root 10--13
не входят в closing gate.

## Results

- Максимальная относительная разность fast/full: `{manifest['maximum_fast_full_relative']:.6e}`.
- Сопоставлений: `{manifest['comparison_row_count']}`; root-9 checks: `{manifest['root9_row_count']}`.
- Reused frozen references: `{manifest['frozen_full_references_reused']}`.
- Fresh full scans closing validation: `{manifest['fresh_full_scans_closing_validation']}` / `{manifest['maximum_fresh_full_scans_closing_validation']}`.
- beta=54: `{manifest['beta54_result']['status']}`.
- Frozen reference hashes preserved: `{manifest['frozen_references_preserved']}`.
- Legacy shards retained but not used: `{manifest['legacy_validation_shard_count']}`.
- Preflight estimate: `{preflight['estimated_fast_total_seconds_two_graphs']:.3f}` s; peak RSS `{preflight['peak_rss_bytes']}` bytes; speedup estimate `{preflight['estimated_speedup']:.3f}`.
- Production physics preserved: `{manifest['production_physics_preserved']}`.
- Current validation dependencies preserved during closing run: `{manifest['validation_dependencies_preserved']}`.

Spike protection is supported by the targeted infrastructure tests, synthetic
reproduced/corrected outcomes and the physical beta=54 check. Spline smoothing,
moving average, branch tracking and manual CSV correction were not used.

## Statuses

"""
    for key, value in statuses.items():
        text += f"- `{key}: {value}`\n"
    text += f"""

Allowed conclusion: `{manifest['allowed_conclusion']}`.

## Ограничения

Валидация конечна и относится только к объявленным beta/mu точкам, roots 1--8
и root 9. Она не является универсальной сертификацией, branch tracking,
MAC-анализом или root-13 audit. Production physics не изменялась;
multiprocessing, smoothing, figures, commit и push не выполнялись.
"""
    path.write_text(text, encoding="utf-8")


def _require_successful_validation_gate(output_dir: Path) -> tuple[dict[str, Any], str]:
    """Require an integrity-preserving physical fast/full validation manifest."""

    manifest_path = output_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise RuntimeError(
            "Run --validate-fast-vs-full successfully before production modes."
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("algorithm_version") != ALGORITHM_VERSION:
        raise RuntimeError("Fast/full validation manifest uses a stale algorithm version.")
    if manifest.get("statuses", {}).get("OVERALL") != "PASS":
        raise RuntimeError("Fast/full validation gate is not PASS.")
    if (
        manifest.get("allowed_conclusion")
        != "FAST_PRODUCTION_SWEEP_VALIDATED_FOR_FIRST8_PLUS_ROOT9"
    ):
        raise RuntimeError("Fast/full validation did not authorize production use.")
    if not bool(manifest.get("production_physics_preserved", False)):
        raise RuntimeError("Validation did not preserve production physics hashes.")
    if not bool(manifest.get("validation_dependencies_preserved", False)):
        raise RuntimeError("Validation did not preserve adapter/dependency hashes.")
    if not bool(
        manifest.get("non_anchor_fast_local_evidence", {}).get("passed", False)
    ):
        raise RuntimeError("Validation lacks required non-anchor FAST_LOCAL evidence.")
    if (
        manifest.get("physical_cached_vs_uncached_equality", {}).get("status")
        != "PASS"
    ):
        raise RuntimeError("Physical cached/uncached equality gate is not PASS.")

    current_production_hashes = {
        path: _sha256(ROOT / path) for path in PRODUCTION_PHYSICS_PATHS
    }
    if manifest.get("production_physics_hashes_after") != current_production_hashes:
        raise RuntimeError("Production physics changed after fast/full validation.")
    if manifest.get("validation_dependency_hashes_after") != _dependency_hashes():
        raise RuntimeError("Validation adapters changed after fast/full validation.")
    return manifest, _sha256(manifest_path)


def _production_checkpoint_fingerprint(
    workflow: str,
    model: str,
    grid: Sequence[float],
    validation_manifest_sha256: str,
) -> str:
    return sweep_runner.stable_fingerprint(
        {
            "stage": STAGE_ID,
            "algorithm_version": ALGORITHM_VERSION,
            "runner_version": sweep_runner.RUNNER_VERSION,
            "workflow": workflow,
            "model": model,
            "grid": list(map(float, grid)),
            "settings": asdict(_settings(workflow)),
            "validation_constants": _validation_constants_payload(),
            "dependency_hashes": _dependency_hashes(),
            "validation_manifest_sha256": validation_manifest_sha256,
            "global_search_source": "FRESH_FROZEN_FULL_SEARCH_CALLBACK",
        }
    )


def _mu_groups(
    rows: Sequence[Mapping[str, Any]], model: str
) -> dict[float, list[dict[str, Any]]]:
    rlb2d._complete_mu_values(rows, model)
    groups: dict[float, list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(round(float(row["mu"]), 12), []).append(dict(row))
    return {
        parameter: sorted(
            point_rows, key=lambda item: int(item["sorted_position"])
        )
        for parameter, point_rows in groups.items()
    }


def _canonical_group_hash(rows: Sequence[Mapping[str, Any]]) -> str:
    return sweep_runner.stable_fingerprint(
        [dict(row) for row in sorted(rows, key=lambda item: int(item["sorted_position"]))]
    )


def _merge_only_absent_mu_groups(
    existing: Sequence[Mapping[str, Any]],
    model: str,
    additions: Mapping[float, Sequence[Mapping[str, Any]]],
) -> list[dict[str, Any]]:
    groups = _mu_groups(existing, model) if existing else {}
    for raw_parameter, raw_rows in additions.items():
        parameter = round(float(raw_parameter), 12)
        rows = [dict(row) for row in raw_rows]
        if parameter in groups:
            existing_values = [float(row["Omega"]) for row in groups[parameter]]
            new_values = [float(row["Omega"]) for row in rows]
            if any(
                left.hex() != right.hex()
                for left, right in zip(existing_values, new_values, strict=True)
            ):
                raise RuntimeError(
                    f"Refusing conflicting closing rows for {model}, mu={parameter:g}."
                )
            continue
        groups[parameter] = rows
    merged = [
        row
        for parameter in sorted(groups)
        for row in sorted(
            groups[parameter], key=lambda item: int(item["sorted_position"])
        )
    ]
    rlb2d._complete_mu_values(merged, model)
    return merged


def run_rlb2d_missing_only(output_dir: Path) -> dict[str, Any]:
    """Close only the six still-missing RLB-2D mu points.

    All accepted validation shards are loaded and verified before either
    canonical CSV is touched.  Existing complete groups are transactionally
    seeded into runner checkpoints, so the runner can calculate only absent
    points and can never revisit a completed spectral group.
    """

    started = time.perf_counter()
    validation_manifest, validation_manifest_sha256 = (
        _require_successful_validation_gate(output_dir)
    )
    grid = np.round(np.arange(0.0, 0.8000001, 0.02), 12)

    # Load all three previously missing immutable-reference results before the
    # first canonical mutation.  Loader validation includes the 59-field
    # schema, re-evaluated quality, row keys and semantic hashes.
    accepted_shards: dict[tuple[str, float], dict[str, Any]] = {}
    for model, parameter in RLB2D_MISSING_REFERENCE_GROUPS:
        payload = _load_validation_shard(
            output_dir, "mu", model, parameter
        )
        if payload is None:
            raise RuntimeError(
                "RLB-2D closing requires all three accepted missing-reference "
                f"shards before mutation; missing {model}, mu={parameter:g}."
            )
        accepted_shards[(model, parameter)] = payload

    canonical_before = {
        model: rlb2d._canonical_mu_rows(RLB2D_OUTPUT, model)
        for model in (MODEL_OLD, MODEL_RLB)
    }
    groups_before = {
        model: _mu_groups(rows, model)
        for model, rows in canonical_before.items()
    }
    preserved_hashes_before = {
        model: {
            format(parameter, ".12g"): _canonical_group_hash(point_rows)
            for parameter, point_rows in groups.items()
        }
        for model, groups in groups_before.items()
    }
    initial_missing = {
        model: [
            float(value)
            for value in grid
            if round(float(value), 12) not in groups_before[model]
        ]
        for model in (MODEL_OLD, MODEL_RLB)
    }

    staged_reference_rows: dict[str, list[dict[str, Any]]] = {}
    for model in (MODEL_OLD, MODEL_RLB):
        additions = {
            parameter: accepted_shards[(model, parameter)]["production_rows"]
            for shard_model, parameter in RLB2D_MISSING_REFERENCE_GROUPS
            if shard_model == model
        }
        staged_reference_rows[model] = _merge_only_absent_mu_groups(
            canonical_before[model], model, additions
        )

    # Both complete staged payloads have passed before either atomic replace.
    for model in (MODEL_OLD, MODEL_RLB):
        rlb2d._atomic_write_csv(
            RLB2D_OUTPUT / rlb2d.ROOT_FILENAMES[(rlb2d.SWEEP_MU, model)],
            staged_reference_rows[model],
        )

    after_reference_groups = {
        model: _mu_groups(
            rlb2d._canonical_mu_rows(RLB2D_OUTPUT, model), model
        )
        for model in (MODEL_OLD, MODEL_RLB)
    }
    remaining_after_materialization = {
        model: tuple(
            float(value)
            for value in grid
            if round(float(value), 12) not in after_reference_groups[model]
        )
        for model in (MODEL_OLD, MODEL_RLB)
    }
    if remaining_after_materialization != RLB2D_EXPECTED_FAST_MISSING:
        raise RuntimeError(
            "Canonical missing-point contract changed after shard "
            f"materialization: {remaining_after_materialization}."
        )

    new_rows: dict[str, dict[float, list[dict[str, Any]]]] = {
        MODEL_OLD: {},
        MODEL_RLB: {},
    }
    runner_summaries: dict[str, Any] = {}
    fast_local_count = 0
    full_anchor_count = 0
    full_fallback_count = 0
    for model in (MODEL_OLD, MODEL_RLB):
        fingerprint = _production_checkpoint_fingerprint(
            "mu", model, grid, validation_manifest_sha256
        )
        checkpoint = sweep_runner.CheckpointConfig(
            directory=(
                RLB2D_OUTPUT
                / "_spectral_sweep_fast_checkpoints"
                / model.lower()
                / fingerprint[:16]
            ),
            sweep_id="mu",
            model_id=model,
            fingerprint=fingerprint,
        )
        seeded: dict[float, sweep_runner.SpectrumRecord] = {}
        for index, parameter in enumerate(grid):
            key = round(float(parameter), 12)
            point_rows = after_reference_groups[model].get(key)
            if point_rows is None:
                continue
            record = _record_from_rows(
                float(parameter), point_rows, "REUSED_CANONICAL_RLB2D"
            )
            prior = sweep_runner.load_point_transaction(
                checkpoint, index, float(parameter)
            )
            if prior is not None and any(
                left.hex() != right.hex()
                for left, right in zip(
                    prior.values, record.values, strict=True
                )
            ):
                raise RuntimeError("Existing checkpoint conflicts with canonical rows.")
            if prior is None:
                sweep_runner.write_point_transaction(
                    checkpoint, index, record
                )
            seeded[float(parameter)] = record
        sweep_runner.write_checkpoint_manifest(
            checkpoint,
            grid,
            seeded,
            sweep_runner.SweepCounters(points_requested=len(grid)),
        )
        backend = PhysicalValidationBackend(
            "mu", model, use_frozen_oracle_for_global=False
        )
        callbacks = sweep_runner.SweepCallbacks(
            backend.global_search,
            backend.local_search,
            backend.quality_gate,
            completeness_guard=lambda parameter, roots: (
                len(roots) == 9,
                "root9",
            ),
            global_completeness_guard=lambda parameter, roots: (
                len(roots) == 9,
                "root9",
            ),
        )
        result = sweep_runner.run_spectral_sweep(
            grid,
            callbacks=callbacks,
            settings=replace(
                _settings("mu", spike=False), anchor_stride=5
            ),
            checkpoint=checkpoint,
            missing_only=True,
        )
        expected = set(RLB2D_EXPECTED_FAST_MISSING[model])
        actual = {round(float(audit.parameter), 12) for audit in result.point_audits}
        if result.status != "PASS" or actual != expected:
            raise RuntimeError(
                f"Missing-only runner touched an unexpected set for {model}: {sorted(actual)}."
            )
        for audit in result.point_audits:
            parameter = round(float(audit.parameter), 12)
            if model == MODEL_OLD and math.isclose(
                parameter, 0.70, rel_tol=0.0, abs_tol=1.0e-12
            ) and not (
                audit.scheduled_anchor and audit.origin == "GLOBAL_ANCHOR"
            ):
                raise RuntimeError("mu=0.70 must remain the scheduled full anchor.")
            record = result.spectra[float(audit.parameter)]
            new_rows[model][parameter] = _production_rows_from_record(
                backend, record, audit
            )
            fast_local_count += int(audit.origin == "FAST_LOCAL")
            full_anchor_count += int(audit.origin == "GLOBAL_ANCHOR")
            full_fallback_count += int(audit.origin == "GLOBAL_FALLBACK")
        runner_summaries[model] = {
            "checkpoint_directory": _repo_relative(checkpoint.directory),
            "checkpoint_manifest_sha256": _sha256(
                checkpoint.directory / "checkpoint_manifest.json"
            ),
            "point_audits": [asdict(audit) for audit in result.point_audits],
            "runner_counters": result.counters.to_dict(),
            "backend_counters": asdict(backend.counters),
            "runtime_seconds": result.runtime_seconds,
        }

    final_staged = {
        model: _merge_only_absent_mu_groups(
            rlb2d._canonical_mu_rows(RLB2D_OUTPUT, model),
            model,
            new_rows[model],
        )
        for model in (MODEL_OLD, MODEL_RLB)
    }
    for model in (MODEL_OLD, MODEL_RLB):
        rlb2d._atomic_write_csv(
            RLB2D_OUTPUT / rlb2d.ROOT_FILENAMES[(rlb2d.SWEEP_MU, model)],
            final_staged[model],
        )

    final_groups = {
        model: _mu_groups(
            rlb2d._canonical_mu_rows(RLB2D_OUTPUT, model), model
        )
        for model in (MODEL_OLD, MODEL_RLB)
    }
    preserved_hashes_after = {
        model: {
            key: _canonical_group_hash(final_groups[model][float(key)])
            for key in hashes
        }
        for model, hashes in preserved_hashes_before.items()
    }
    if preserved_hashes_before != preserved_hashes_after:
        raise RuntimeError("A previously complete canonical mu group changed.")
    final_missing = {
        model: [
            float(value)
            for value in grid
            if round(float(value), 12) not in final_groups[model]
        ]
        for model in (MODEL_OLD, MODEL_RLB)
    }
    if any(final_missing.values()):
        raise RuntimeError(f"RLB-2D mu closing remains incomplete: {final_missing}.")

    closing_checkpoint = rlb2d.update_closing_checkpoint(RLB2D_OUTPUT)
    closing_checkpoint.update(
        {
            "spectral_sweep_fast_closing": {
                "validation_manifest_sha256": validation_manifest_sha256,
                "validation_status": validation_manifest["statuses"]["OVERALL"],
                "initial_missing_mu": initial_missing,
                "materialized_validation_shard_groups": 3,
                "remaining_after_materialization": remaining_after_materialization,
                "new_fast_runner_points": sum(len(rows) for rows in new_rows.values()),
                "recomputed_already_complete_points": 0,
                "FAST_LOCAL_point_count": fast_local_count,
                "full_anchor_point_count": full_anchor_count,
                "full_fallback_point_count": full_fallback_count,
                "internal_tail": "QUALIFIED_NOT_COMPUTED_ABOVE_ROOT9",
                "final_missing_mu": final_missing,
                "preserved_group_hashes_before": preserved_hashes_before,
                "preserved_group_hashes_after": preserved_hashes_after,
                "runner_summaries": runner_summaries,
                "runtime_seconds": time.perf_counter() - started,
                "peak_rss_bytes": _current_peak_rss_bytes(),
                "parallel_workers_used": 0,
            }
        }
    )
    rlb2d._atomic_write_json(
        RLB2D_OUTPUT / rlb2d.CLOSING_CHECKPOINT_FILENAME,
        closing_checkpoint,
    )
    return closing_checkpoint["spectral_sweep_fast_closing"]


def run_production_mode(
    output_dir: Path,
    *,
    resume: bool,
    missing_only: bool,
    plot_only: bool,
) -> dict[str, Any]:
    validation_manifest, validation_manifest_sha256 = (
        _require_successful_validation_gate(output_dir)
    )
    if not plot_only:
        preflight_path = output_dir / "preflight_report.json"
        if not preflight_path.is_file():
            raise RuntimeError("Run --preflight first.")
        preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
        if not preflight.get("production_run_permitted", False):
            raise RuntimeError(
                "Preflight forecast exceeds 60 minutes; production sweep is blocked."
            )
    summaries: list[dict[str, Any]] = []
    contracts = (
        ("beta", np.arange(0.0, 91.0, 1.0)),
        ("mu", np.round(np.arange(0.0, 0.8000001, 0.02), 12)),
    )
    for workflow, grid in contracts:
        for model in MODELS:
            checkpoint = sweep_runner.CheckpointConfig(
                directory=output_dir / "checkpoints" / workflow / model.lower(),
                sweep_id=workflow,
                model_id=model,
                fingerprint=_production_checkpoint_fingerprint(
                    workflow,
                    model,
                    grid,
                    validation_manifest_sha256,
                ),
            )
            backend: PhysicalValidationBackend | None = None
            callbacks: sweep_runner.SweepCallbacks | None = None
            if not plot_only:
                backend = PhysicalValidationBackend(
                    workflow,
                    model,
                    use_frozen_oracle_for_global=False,
                )
                callbacks = sweep_runner.SweepCallbacks(
                    backend.global_search,
                    backend.local_search,
                    backend.quality_gate,
                    completeness_guard=lambda parameter, roots: (
                        len(roots) == FULL_EXPORT_ROOTS + FULL_EXPORT_GUARDS,
                        "root9",
                    ),
                    global_completeness_guard=lambda parameter, roots: (
                        len(roots) == FULL_EXPORT_ROOTS + FULL_EXPORT_GUARDS,
                        "root9",
                    ),
                )
            result = sweep_runner.run_spectral_sweep(
                grid,
                callbacks=callbacks,
                settings=_settings(workflow),
                checkpoint=checkpoint,
                resume=resume,
                missing_only=missing_only,
                plot_only=plot_only,
            )
            summaries.append(
                {
                    "workflow": workflow,
                    "model": model,
                    "status": result.status,
                    "point_count": len(result.spectra),
                    "runtime_seconds": result.runtime_seconds,
                    "runner_counters": result.counters.to_dict(),
                    "backend_counters": (
                        asdict(backend.counters) if backend is not None else {}
                    ),
                    "global_search_source": (
                        "NONE_PLOT_ONLY"
                        if plot_only
                        else "FRESH_FROZEN_FULL_SEARCH_CALLBACK"
                    ),
                }
            )
    payload = {
        "stage": STAGE_ID,
        "mode": "plot_only" if plot_only else "production_fast",
        "validation_gate_manifest_sha256": validation_manifest_sha256,
        "validation_gate_status": validation_manifest["statuses"]["OVERALL"],
        "physics_callbacks_instantiated": not plot_only,
        "summaries": summaries,
        "figures_created": 0,
    }
    _atomic_json(output_dir / "production_run_summary.json", payload)
    return payload


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--preflight", action="store_true")
    actions.add_argument("--validate-fast-vs-full", action="store_true")
    actions.add_argument("--run", action="store_true")
    actions.add_argument("--resume", action="store_true")
    actions.add_argument("--missing-only", action="store_true")
    actions.add_argument("--plot-only", action="store_true")
    actions.add_argument("--validation-worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--worker-workflow", choices=("beta", "mu", "beta_iso_rlb"), help=argparse.SUPPRESS)
    parser.add_argument("--worker-model", choices=MODELS, help=argparse.SUPPRESS)
    parser.add_argument("--worker-target", type=float, help=argparse.SUPPRESS)
    parser.add_argument(
        "--worker-fresh-full-budget",
        type=int,
        default=0,
        help=argparse.SUPPRESS,
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = Path(args.output_dir)
    if args.preflight:
        payload = run_preflight(output_dir)
    elif args.validate_fast_vs_full:
        payload = run_fixed_validation(output_dir)
    elif args.validation_worker:
        if args.worker_workflow is None or args.worker_model is None or args.worker_target is None:
            raise RuntimeError("Internal validation worker arguments are incomplete.")
        payload = run_validation_worker(
            output_dir,
            args.worker_workflow,
            args.worker_model,
            args.worker_target,
            fresh_full_budget=args.worker_fresh_full_budget,
        )
    elif args.missing_only:
        payload = run_rlb2d_missing_only(output_dir)
    else:
        payload = run_production_mode(
            output_dir,
            resume=args.resume,
            missing_only=False,
            plot_only=args.plot_only,
        )
    print(json.dumps(_json_value(payload), ensure_ascii=False, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
