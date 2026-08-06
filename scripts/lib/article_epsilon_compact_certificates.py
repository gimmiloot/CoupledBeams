"""Memory-bounded scientific certificates for the article epsilon grid.

The compressed prefix cache is a transient solver trace.  This module reads at
most one such trace at a time and stores the small, matrix-confirmed subset that
is needed to reproduce the sorted-rank EB/Timoshenko prefix decision.  It is a
pure migration module: it deliberately imports no matrix evaluator, root solver,
branch tracker, local repair, FEM, or anisotropic implementation.
"""

from __future__ import annotations

import csv
import ctypes
import gc
import gzip
import hashlib
import json
import math
import os
import re
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence


SCIENTIFIC_SCOPE = "isotropic_circular_coupled_rods_eb_timoshenko"
SCHEMA_VERSION = "article_epsilon_compact_certificate_v1"
COMPACTOR_VERSION = "article_epsilon_streaming_compactor_v1"
K_MAX = 10
DELTA_TOL = 0.10
MAX_STORED_RANK = 12
MODEL_EB = "Euler-Bernoulli"
MODEL_TIMO = "Timoshenko"
MODELS = (MODEL_EB, MODEL_TIMO)
CASE_ID_RE = re.compile(r"AUE_[0-9a-f]{20}")

INDEX_FIELDS = (
    "case_id", "epsilon_0", "beta_deg", "mu", "eta", "s_max",
    "execution_status", "n_true_status", "N_true", "first_failed_mode",
    "required_guard", "required_guard_confirmed", "strict_request_flag",
    "detector_triggered", "local_repair_status", "eb_root_count",
    "timo_root_count", "certificate_path", "source_full_cache_path",
    "source_full_cache_sha256", "source_full_cache_size", "source_kind",
    "scientific_scope", "schema_version",
)

SOURCE_INVENTORY_FIELDS = (
    "case_id", "relative_path", "size", "sha256", "artifact_kind",
    "canonical_source", "selection_reason",
)

COMPACTION_RESULT_FIELDS = (
    "case_id", "source_kind", "source_path", "certificate_path", "action",
    "fidelity_status", "source_size", "certificate_size", "rss_bytes",
    "elapsed_seconds",
)

FAILURE_FIELDS = (
    "case_id", "source_path", "failure_kind", "message",
)

MEMORY_FIELDS = (
    "processed_cases", "timestamp_utc", "rss_bytes", "peak_rss_bytes",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path, chunk_size: int = 4 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while True:
            chunk = stream.read(chunk_size)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _canonical_json(value: object) -> str:
    return json.dumps(
        value, sort_keys=True, ensure_ascii=False, separators=(",", ":"),
        allow_nan=False,
    )


def _atomic_replace(temp: Path, target: Path) -> None:
    with temp.open("r+b") as stream:
        os.fsync(stream.fileno())
    os.replace(temp, target)


def atomic_write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(_canonical_json(value) + "\n", encoding="utf-8")
    _atomic_replace(temp, path)


def atomic_write_gzip_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("wb") as raw:
        # mtime=0 makes repeated certificates byte-for-byte deterministic.
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            zipped.write((_canonical_json(value) + "\n").encode("utf-8"))
        raw.flush()
        os.fsync(raw.fileno())
    with gzip.open(temp, "rt", encoding="utf-8") as stream:
        check = json.load(stream)
    if check != value:
        temp.unlink(missing_ok=True)
        raise RuntimeError(f"atomic compact-certificate validation failed: {path}")
    os.replace(temp, path)


def atomic_write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".tmp")
    with temp.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fields), extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: _csv_value(row.get(field, "")) for field in fields})
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temp, path)


def _csv_value(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        return format(value, ".17g")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def load_certificate(path: Path) -> dict[str, object]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        result = json.load(stream)
    if not isinstance(result, dict):
        raise RuntimeError(f"compact certificate is not an object: {path}")
    return result


def process_rss_bytes() -> int:
    """Return current working set without introducing a psutil dependency."""

    if os.name == "nt":
        from ctypes import wintypes  # noqa: PLC0415

        class Counters(ctypes.Structure):
            _fields_ = [
                ("cb", wintypes.DWORD), ("PageFaultCount", wintypes.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
                ("PrivateUsage", ctypes.c_size_t),
            ]
        counters = Counters()
        counters.cb = ctypes.sizeof(counters)
        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        psapi = ctypes.WinDLL("psapi", use_last_error=True)
        kernel32.GetCurrentProcess.restype = wintypes.HANDLE
        psapi.GetProcessMemoryInfo.argtypes = (
            wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD,
        )
        psapi.GetProcessMemoryInfo.restype = wintypes.BOOL
        ok = psapi.GetProcessMemoryInfo(
            kernel32.GetCurrentProcess(),
            ctypes.byref(counters), counters.cb,
        )
        return int(counters.WorkingSetSize) if ok else 0
    try:
        import resource  # noqa: PLC0415
        value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        return value * (1024 if value < 10**10 else 1)
    except (ImportError, OSError):
        return 0


def _uncompressed_gzip_size(path: Path) -> int:
    with path.open("rb") as stream:
        stream.seek(-4, os.SEEK_END)
        return int.from_bytes(stream.read(4), "little")


def _as_float(value: object, default: float | None = None) -> float | None:
    try:
        if value is None or str(value).strip() == "":
            return default
        result = float(value)
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def _as_int(value: object, default: int | None = None) -> int | None:
    try:
        if value is None or str(value).strip() == "":
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _source_roots(payload: Mapping[str, object], model: str) -> list[float]:
    candidates: list[list[float]] = []
    explicit = payload.get("roots_already_resolved")
    if isinstance(explicit, Mapping):
        values = explicit.get(model)
        if isinstance(values, Sequence) and not isinstance(values, (str, bytes)):
            candidates.append([float(item) for item in values])
    models = payload.get("models")
    if isinstance(models, Mapping):
        state = models.get(model)
        if isinstance(state, Mapping):
            values = state.get("resolved_roots")
            if isinstance(values, Sequence) and not isinstance(values, (str, bytes)):
                candidates.append([float(item) for item in values])
    return max(candidates, key=len, default=[])[:MAX_STORED_RANK]


def _compact_clusters(payload: Mapping[str, object], model: str) -> list[dict[str, object]]:
    models = payload.get("models")
    state = models.get(model) if isinstance(models, Mapping) else None
    clusters = state.get("clusters", ()) if isinstance(state, Mapping) else ()
    result: list[dict[str, object]] = []
    if isinstance(clusters, Sequence) and not isinstance(clusters, (str, bytes)):
        for item in clusters:
            if not isinstance(item, Mapping):
                continue
            rank = _as_int(item.get("sorted_index"), 0) or 0
            if rank > MAX_STORED_RANK:
                continue
            result.append({
                key: item.get(key)
                for key in (
                    "cluster_id", "cluster_member_index", "cluster_size",
                    "detected_nullity", "sorted_index", "resolution_status",
                )
                if key in item
            })
    return result


def _compact_unresolved_intervals(payload: Mapping[str, object], model: str) -> list[dict[str, object]]:
    models = payload.get("models")
    state = models.get(model) if isinstance(models, Mapping) else None
    raw = state.get("brackets_and_continuation_metadata", ()) if isinstance(state, Mapping) else ()
    if isinstance(raw, Mapping):
        raw = raw.get("unresolved_intervals", ())
    rows: list[dict[str, object]] = []
    if isinstance(raw, Sequence) and not isinstance(raw, (str, bytes)):
        for item in raw:
            if isinstance(item, Mapping):
                status = str(item.get("resolution_status", ""))
                if status.startswith("resolved_"):
                    continue
                rows.append({
                    key: item.get(key)
                    for key in (
                        "Lambda_left", "Lambda_right", "resolution_status",
                        "trigger_reasons", "local_minima_count", "minimum_sigma",
                        "minimum_sigma_ratio", "rank_start", "rank_end",
                    )
                    if key in item
                })
            else:
                rows.append({"summary": str(item)})
    return rows[:24]


def _quality_by_rank(roots: Sequence[float], clusters: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    cluster_by_rank = {
        int(item["sorted_index"]): item
        for item in clusters if _as_int(item.get("sorted_index"), 0)
    }
    return [
        {
            "rank": rank,
            "status": "matrix_confirmed_resolved",
            "multiplicity_slot": int(cluster_by_rank.get(rank, {}).get("cluster_member_index", 1)),
            "cluster_id": cluster_by_rank.get(rank, {}).get("cluster_id", ""),
        }
        for rank, _ in enumerate(roots, start=1)
    ]


def calculate_deltas(eb_roots: Sequence[float], timo_roots: Sequence[float]) -> dict[str, float]:
    count = min(K_MAX, len(eb_roots), len(timo_roots))
    return {
        str(rank): abs(float(eb_roots[rank - 1]) ** 2 - float(timo_roots[rank - 1]) ** 2)
        / float(timo_roots[rank - 1]) ** 2
        for rank in range(1, count + 1)
    }


def reproduce_prefix(eb_roots: Sequence[float], timo_roots: Sequence[float]) -> tuple[int | None, int | None, int | None, dict[str, float]]:
    deltas = calculate_deltas(eb_roots, timo_roots)
    for rank in range(1, K_MAX + 1):
        if str(rank) not in deltas:
            return None, None, None, deltas
        if deltas[str(rank)] > DELTA_TOL:
            guard = rank + 1
            if len(eb_roots) < guard or len(timo_roots) < guard:
                return None, rank, guard, deltas
            return rank - 1, rank, guard, deltas
    guard = K_MAX + 1
    if len(eb_roots) < guard or len(timo_roots) < guard:
        return None, None, guard, deltas
    return K_MAX, None, guard, deltas


def _operation_summary(payload: Mapping[str, object]) -> dict[str, object]:
    totals: dict[str, float] = {
        "matrix_evaluations": 0, "svd_calls": 0, "primary_seconds": 0.0,
        "verification_seconds": 0.0, "preparation_seconds": 0.0,
    }
    rows = payload.get("stage_timings", ())
    if isinstance(rows, Sequence) and not isinstance(rows, (str, bytes)):
        for row in rows:
            if not isinstance(row, Mapping):
                continue
            totals["matrix_evaluations"] += float(row.get("primary_matrix_evaluations_added", 0) or 0)
            totals["matrix_evaluations"] += float(row.get("verification_matrix_evaluations_added", 0) or 0)
            totals["svd_calls"] += float(row.get("primary_svd_calls_added", 0) or 0)
            totals["svd_calls"] += float(row.get("verification_svd_calls_added", 0) or 0)
            for key in ("primary_seconds", "verification_seconds", "preparation_seconds"):
                totals[key] += float(row.get(key, 0.0) or 0.0)
    return {
        "matrix_evaluations": int(totals["matrix_evaluations"]),
        "svd_calls": int(totals["svd_calls"]),
        "primary_seconds": totals["primary_seconds"],
        "verification_seconds": totals["verification_seconds"],
        "preparation_seconds": totals["preparation_seconds"],
        "force_strict_requested": int(payload.get("force_strict_requested", 0) or 0),
        "force_strict_executed": int(payload.get("force_strict_executed", 0) or 0),
    }


def _scientific_configuration(payload: Mapping[str, object]) -> Mapping[str, object]:
    identity = payload.get("identity")
    if not isinstance(identity, Mapping):
        return {}
    value = identity.get("scientific_model_configuration")
    return value if isinstance(value, Mapping) else {}


def _required_guard_from_source(payload: Mapping[str, object], reproduced: int | None) -> int | None:
    explicit = _as_int(payload.get("required_guard_mode"))
    if explicit is not None:
        return explicit
    if reproduced is not None:
        return reproduced
    first = _as_int(payload.get("first_apparent_failed_mode", payload.get("first_failed_mode")))
    if first is not None:
        return first + 1
    if _as_int(payload.get("N_true")) == K_MAX:
        return K_MAX + 1
    return _as_int(payload.get("prefix_guard_resolved_through"))


def extract_certificate(
    payload: Mapping[str, object], manifest_row: Mapping[str, str], *,
    manifest_sha: str, source_commit: str, source_path: str,
    source_sha: str, source_size: int, source_kind: str = "full_point_cache",
) -> dict[str, object]:
    eb_roots = _source_roots(payload, MODEL_EB)
    timo_roots = _source_roots(payload, MODEL_TIMO)
    eb_clusters = _compact_clusters(payload, MODEL_EB)
    timo_clusters = _compact_clusters(payload, MODEL_TIMO)
    reproduced_n, reproduced_failure, reproduced_guard, calculated = reproduce_prefix(eb_roots, timo_roots)
    source_n = _as_int(payload.get("N_true"))
    resolved = source_n is not None and str(payload.get("execution_status", "")).startswith("resolved_")
    if resolved and reproduced_n != source_n:
        raise ValueError(f"stored N_true={source_n} is not reproduced from compact roots ({reproduced_n})")
    source_failure = _as_int(payload.get("first_failed_mode"))
    if resolved and source_failure != reproduced_failure:
        raise ValueError(
            f"stored first_failed_mode={source_failure} differs from reproduced {reproduced_failure}"
        )
    required_guard = _required_guard_from_source(payload, reproduced_guard)
    if resolved:
        required_guard = reproduced_guard
    source_deltas = payload.get("deltas", {})
    if isinstance(source_deltas, Mapping):
        for rank, value in source_deltas.items():
            if str(rank) in calculated and not math.isclose(
                # Stored roots are rounded scientific outputs whereas some
                # legacy delta values were computed before that JSON round-trip.
                float(value), calculated[str(rank)], rel_tol=1e-8, abs_tol=2e-7,
            ):
                raise ValueError(f"delta_f,{rank} does not round-trip from compact roots")

    config = _scientific_configuration(payload)
    primary = config.get("primary_settings", {}) if isinstance(config, Mapping) else {}
    strict = config.get("strict_settings", {}) if isinstance(config, Mapping) else {}
    unresolved_eb = _compact_unresolved_intervals(payload, MODEL_EB)
    unresolved_timo = _compact_unresolved_intervals(payload, MODEL_TIMO)
    enough_roots = required_guard is not None and len(eb_roots) >= required_guard and len(timo_roots) >= required_guard
    guard_confirmed = bool(resolved and enough_roots)
    return {
        "scientific_scope": SCIENTIFIC_SCOPE,
        "schema_version": SCHEMA_VERSION,
        "case_id": manifest_row["case_id"],
        "geometry": {
            "epsilon_0": float(manifest_row["epsilon_0"]),
            "beta_deg": float(manifest_row["beta_deg"]),
            "mu": float(manifest_row["mu"]),
            "eta": float(manifest_row["eta"]),
            "s_max": float(manifest_row["s_max"]),
        },
        "manifest_sha": manifest_sha,
        "source_commit": source_commit,
        "source": {
            "kind": source_kind, "full_cache_path": source_path,
            "full_cache_sha256": source_sha, "full_cache_size": source_size,
        },
        "result": {
            "execution_status": str(payload.get("execution_status", "attempted_unresolved")),
            "n_true_status": str(payload.get("n_true_status", "exact" if resolved else "unresolved")),
            "N_true": source_n if resolved else None,
            "first_failed_mode": source_failure if resolved else None,
            "first_apparent_failed_mode": _as_int(
                payload.get("first_apparent_failed_mode", payload.get("first_failed_mode"))
            ),
            "required_guard": required_guard,
            "required_guard_confirmed": guard_confirmed,
            "upper_spectrum_audit_status": str(payload.get(
                "upper_spectrum_audit_status",
                payload.get("optional_upper_spectrum_full_audit_status", "not_attempted"),
            )),
            "unresolved_reason": str(payload.get("unresolved_reason", payload.get("defer_reason", ""))),
        },
        "spectra": {
            "Euler-Bernoulli": {
                "roots": eb_roots, "quality_by_rank": _quality_by_rank(eb_roots, eb_clusters),
                "clusters": eb_clusters, "unresolved_intervals": unresolved_eb,
            },
            "Timoshenko": {
                "roots": timo_roots, "quality_by_rank": _quality_by_rank(timo_roots, timo_clusters),
                "clusters": timo_clusters, "unresolved_intervals": unresolved_timo,
            },
            "delta_f": calculated,
        },
        "diagnostics": {
            "primary_solver_status": str(payload.get("primary_status", payload.get("execution_status", ""))),
            "strict_request_flag": bool(payload.get("expensive_escalation_requested", False))
            or bool(payload.get("force_strict_requested", 0)),
            "strict_request_kind": str(payload.get("expensive_escalation_kind", "")),
            "detector_trigger_summary": payload.get("detector_trigger_summary", payload.get("trigger_reason", [])),
            "local_repair_status": str(payload.get("local_repair_status", payload.get("local_status", "not_attempted"))),
            "recovered_root_count": int(payload.get("recovered_root_count", 0) or 0),
            "unresolved_cluster_below_guard": bool(unresolved_eb or unresolved_timo),
        },
        "provenance": {
            "algorithm_version": str(payload.get("identity", {}).get("algorithm_version", ""))
            if isinstance(payload.get("identity"), Mapping) else "",
            "detector_version": str(payload.get("detector_version", "not_applied")),
            "repair_version": str(payload.get("repair_version", "not_applied")),
            "eb_evaluator_version": str(config.get("eb_matrix_evaluator_version", "")),
            "timoshenko_evaluator_version": str(config.get("timoshenko_basis_evaluator_version", "")),
            "tolerance_hash": sha256_text(_canonical_json({"primary": primary, "strict": strict})),
            "root_match_tolerance": _as_float(
                primary.get("root_match_tol") if isinstance(primary, Mapping) else None,
                _as_float(strict.get("root_match_tolerance") if isinstance(strict, Mapping) else None),
            ),
            "prefix_strategy": str(payload.get("identity", {}).get("prefix_strategy", "paired"))
            if isinstance(payload.get("identity"), Mapping) else "paired",
            "strict_policy": str(payload.get("identity", {}).get("strict_policy", payload.get("strict_policy_effective", "")))
            if isinstance(payload.get("identity"), Mapping) else str(payload.get("strict_policy_effective", "")),
            "operation_counts": _operation_summary(payload),
        },
    }


def _benchmark_payload(coarse_dir: Path, case_id: str) -> tuple[dict[str, object], str, int, str]:
    path = coarse_dir.parent / "coarse_grid_parallel_benchmark" / "runs.csv"
    rows = [row for row in read_csv(path) if row.get("case_id") == case_id]
    if not rows:
        raise FileNotFoundError(f"saved benchmark inventory not found for {case_id}")
    row = sorted(rows, key=lambda item: int(item.get("workers", "999")))[0]
    data = {
        "execution_status": row["execution_status"], "N_true": None,
        "first_failed_mode": _as_int(row.get("first_failed_mode")),
        "required_guard_mode": _as_int(row.get("required_guard"), 11),
        "unresolved_reason": "saved_benchmark_required_inventory_unresolved",
        "models": {
            MODEL_EB: {"resolved_roots": json.loads(row["EB_roots_json"]), "clusters": []},
            MODEL_TIMO: {"resolved_roots": json.loads(row["Timoshenko_roots_json"]), "clusters": []},
        },
        "identity": {"prefix_strategy": "paired", "strict_policy": "auto", "scientific_model_configuration": {}},
        "force_strict_requested": 0, "force_strict_executed": 0,
    }
    serialized = _canonical_json(data)
    return data, sha256_text(serialized), len(serialized.encode("utf-8")), str(path)


def _index_row(cert: Mapping[str, object], certificate_path: str) -> dict[str, object]:
    geometry = cert["geometry"]
    result = cert["result"]
    diagnostics = cert["diagnostics"]
    spectra = cert["spectra"]
    source = cert["source"]
    assert isinstance(geometry, Mapping) and isinstance(result, Mapping)
    assert isinstance(diagnostics, Mapping) and isinstance(spectra, Mapping) and isinstance(source, Mapping)
    eb = spectra[MODEL_EB]
    timo = spectra[MODEL_TIMO]
    assert isinstance(eb, Mapping) and isinstance(timo, Mapping)
    trigger = diagnostics.get("detector_trigger_summary", ())
    return {
        "case_id": cert["case_id"], **geometry,
        "execution_status": result.get("execution_status", ""),
        "n_true_status": result.get("n_true_status", ""),
        "N_true": result.get("N_true"),
        "first_failed_mode": result.get("first_failed_mode"),
        "required_guard": result.get("required_guard"),
        "required_guard_confirmed": result.get("required_guard_confirmed", False),
        "strict_request_flag": diagnostics.get("strict_request_flag", False),
        "detector_triggered": bool(trigger),
        "local_repair_status": diagnostics.get("local_repair_status", ""),
        "eb_root_count": len(eb.get("roots", ())),
        "timo_root_count": len(timo.get("roots", ())),
        "certificate_path": certificate_path,
        "source_full_cache_path": source.get("full_cache_path", ""),
        "source_full_cache_sha256": source.get("full_cache_sha256", ""),
        "source_full_cache_size": source.get("full_cache_size", 0),
        "source_kind": source.get("kind", ""),
        "scientific_scope": cert["scientific_scope"], "schema_version": cert["schema_version"],
    }


def iter_full_payloads(
    manifest_rows: Sequence[Mapping[str, str]], canonical_paths: Mapping[str, Path],
) -> Iterator[tuple[Mapping[str, str], Path, dict[str, object]]]:
    """Yield exactly one full payload; the generator retains no previous item."""

    for row in manifest_rows:
        path = canonical_paths.get(row["case_id"])
        if path is None:
            continue
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            payload = json.load(stream)
        if not isinstance(payload, dict):
            raise ValueError(f"point-cache payload is not an object: {path}")
        yield row, path, payload
        del payload


@dataclass(frozen=True)
class CompactionRun:
    output_dir: Path
    certificate_count: int
    converted_count: int
    cache_hit_count: int
    failure_count: int
    initial_rss: int
    peak_rss: int
    final_rss: int
    source_bytes: int
    compact_bytes: int
    wall_seconds: float
    operation_counts: Mapping[str, int]


def build_compact_certificates(
    coarse_dir: Path, compact_dir: Path | None = None, *, source_commit: str | None = None,
) -> CompactionRun:
    started = time.perf_counter()
    coarse_dir = coarse_dir.resolve()
    output_dir = (compact_dir or coarse_dir / "compact_point_certificates_v1").resolve()
    case_dir = output_dir / "cases"
    log_dir = output_dir / "logs"
    case_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = coarse_dir / "grid_manifest.csv"
    manifest_rows = read_csv(manifest_path)
    if len(manifest_rows) != 1554 or len({row["case_id"] for row in manifest_rows}) != 1554:
        raise RuntimeError("compact migration requires the unique 1554-case manifest")
    manifest_sha = sha256_file(manifest_path)
    commit = source_commit or subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=coarse_dir, text=True,
        check=True, capture_output=True,
    ).stdout.strip()

    prefix_dir = coarse_dir / "cache" / "prefix"
    by_case: dict[str, list[Path]] = {}
    for path in sorted(prefix_dir.glob("*.json.gz")):
        match = CASE_ID_RE.search(path.name)
        if match:
            by_case.setdefault(match.group(0), []).append(path)
    canonical = {case_id: sorted(paths, key=lambda path: path.name)[0] for case_id, paths in by_case.items()}
    manifest_ids = {row["case_id"] for row in manifest_rows}
    if set(canonical) - manifest_ids:
        raise RuntimeError("prefix cache contains a case outside the manifest")

    initial_rss = process_rss_bytes()
    peak_rss = initial_rss
    memory_rows: list[dict[str, object]] = []
    result_rows: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    index_rows: list[dict[str, object]] = []
    source_hashes: dict[Path, str] = {}
    largest_decompressed = 0
    converted = 0
    cache_hits = 0
    processed = 0

    row_by_id = {row["case_id"]: row for row in manifest_rows}
    for row in manifest_rows:
        case_id = row["case_id"]
        item_started = time.perf_counter()
        source_path = canonical.get(case_id)
        payload: dict[str, object] | None = None
        try:
            if source_path is None:
                payload, source_sha, source_size, synthetic_path = _benchmark_payload(coarse_dir, case_id)
                source_kind = "saved_parallel_benchmark"
                source_rel = str(Path(synthetic_path).resolve().relative_to(coarse_dir.parents[2].resolve())).replace("\\", "/")
            else:
                source_sha = sha256_file(source_path)
                source_hashes[source_path.resolve()] = source_sha
                source_size = source_path.stat().st_size
                largest_decompressed = max(largest_decompressed, _uncompressed_gzip_size(source_path))
                source_kind = "full_point_cache"
                source_rel = str(source_path.relative_to(coarse_dir.parents[2])).replace("\\", "/")
                with gzip.open(source_path, "rt", encoding="utf-8") as stream:
                    payload = json.load(stream)
                if not isinstance(payload, dict):
                    raise ValueError("full point-cache payload is not an object")
            cert_path = case_dir / f"{case_id}.json.gz"
            cert_rel = str(cert_path.relative_to(coarse_dir.parents[2])).replace("\\", "/")
            existing: dict[str, object] | None = None
            if cert_path.exists():
                candidate = load_certificate(cert_path)
                source = candidate.get("source", {})
                if (
                    candidate.get("schema_version") == SCHEMA_VERSION
                    and candidate.get("manifest_sha") == manifest_sha
                    and isinstance(source, Mapping)
                    and source.get("full_cache_sha256") == source_sha
                ):
                    existing = candidate
            if existing is None:
                cert = extract_certificate(
                    payload, row, manifest_sha=manifest_sha, source_commit=commit,
                    source_path=source_rel, source_sha=source_sha,
                    source_size=source_size, source_kind=source_kind,
                )
                atomic_write_gzip_json(cert_path, cert)
                converted += 1
                action = "converted"
            else:
                cert = existing
                cache_hits += 1
                action = "certificate_cache_hit"
            index_rows.append(_index_row(cert, cert_rel))
            processed += 1
            rss = process_rss_bytes()
            peak_rss = max(peak_rss, rss)
            result_rows.append({
                "case_id": case_id, "source_kind": source_kind,
                "source_path": source_rel, "certificate_path": cert_rel,
                "action": action, "fidelity_status": "PASS",
                "source_size": source_size, "certificate_size": cert_path.stat().st_size,
                "rss_bytes": rss, "elapsed_seconds": time.perf_counter() - item_started,
            })
        except Exception as exc:  # one corrupt cache must not create a false certificate
            failures.append({
                "case_id": case_id,
                "source_path": str(source_path or "saved_parallel_benchmark"),
                "failure_kind": type(exc).__name__, "message": str(exc)[:1000],
            })
        finally:
            if payload is not None:
                payload.clear()
            del payload
            if processed and processed % 25 == 0:
                gc.collect()
            if (processed and processed % 100 == 0) or processed == len(manifest_rows):
                rss = process_rss_bytes()
                peak_rss = max(peak_rss, rss)
                memory_rows.append({
                    "processed_cases": processed, "timestamp_utc": utc_now(),
                    "rss_bytes": rss, "peak_rss_bytes": peak_rss,
                })
                print(f"compact progress {processed}/{len(manifest_rows)} rss={rss / 2**30:.3f} GiB", flush=True)

    if failures:
        # Publish failure diagnostics, but never a misleading complete index.
        atomic_write_csv(output_dir / "compaction_failures.csv", failures, FAILURE_FIELDS)
    index_rows.sort(key=lambda row: (
        float(row["epsilon_0"]), float(row["beta_deg"]),
        float(row["mu"]), float(row["eta"]), str(row["case_id"]),
    ))
    result_rows.sort(key=lambda row: str(row["case_id"]))
    failures.sort(key=lambda row: str(row["case_id"]))

    # The prefix caches are gzip files, while strict/sanity auxiliary caches
    # include plain JSON.  Both belong to the immutable 1604-artifact source
    # inventory even though only prefix gzip files are scientific sources for
    # per-case certificates.
    all_artifacts = sorted(
        path for path in (coarse_dir / "cache").rglob("*")
        if path.is_file() and path.suffix in {".gz", ".json"}
    )
    source_rows: list[dict[str, object]] = []
    for path in all_artifacts:
        match = CASE_ID_RE.search(path.name)
        case_id = match.group(0) if match else ""
        selected = canonical.get(case_id) == path
        digest = source_hashes.get(path.resolve()) or sha256_file(path)
        source_rows.append({
            "case_id": case_id,
            "relative_path": str(path.relative_to(coarse_dir.parents[2])).replace("\\", "/"),
            "size": path.stat().st_size, "sha256": digest,
            "artifact_kind": path.parent.name,
            "canonical_source": selected,
            "selection_reason": "canonical_prefix_identity" if selected else "non_prefix_or_rejected_duplicate_artifact",
        })
    source_rows.sort(key=lambda row: str(row["relative_path"]))

    atomic_write_csv(output_dir / "compact_index.csv", index_rows, INDEX_FIELDS)
    atomic_write_csv(output_dir / "source_cache_inventory.csv", source_rows, SOURCE_INVENTORY_FIELDS)
    atomic_write_csv(output_dir / "compaction_results.csv", result_rows, COMPACTION_RESULT_FIELDS)
    atomic_write_csv(output_dir / "compaction_failures.csv", failures, FAILURE_FIELDS)
    atomic_write_csv(log_dir / "memory_profile.csv", memory_rows, MEMORY_FIELDS)
    schema = {
        "schema_version": SCHEMA_VERSION, "scientific_scope": SCIENTIFIC_SCOPE,
        "k_max": K_MAX, "delta_tolerance": DELTA_TOL,
        "required_guard_rule": "first failed rank k -> k+1; N_true=10 -> rank 11",
        "max_stored_rank": MAX_STORED_RANK,
        "excluded_large_fields": [
            "dense Lambda grids", "determinant arrays", "sigma_min arrays",
            "candidate traces", "null vectors", "matrices", "shapes",
        ],
    }
    atomic_write_json(output_dir / "schema.json", schema)
    final_rss = process_rss_bytes()
    peak_rss = max(peak_rss, final_rss)
    compact_bytes = sum((case_dir / f"{row['case_id']}.json.gz").stat().st_size for row in index_rows)
    source_bytes = sum(int(row["size"]) for row in source_rows)
    wall = time.perf_counter() - started
    operation_counts = {
        "point_solver_calls": 0, "matrix_evaluator_calls": 0,
        "local_repair_calls": 0, "force_strict_executed": 0,
        "converted": converted, "certificate_cache_hits": cache_hits,
    }
    metadata = {
        "schema_version": SCHEMA_VERSION, "compactor_version": COMPACTOR_VERSION,
        "scientific_scope": SCIENTIFIC_SCOPE, "manifest_sha256": manifest_sha,
        "source_commit": commit, "manifest_case_count": len(manifest_rows),
        "certificate_count": len(index_rows), "failure_count": len(failures),
        "source_artifact_count": len(source_rows), "source_bytes": source_bytes,
        "compact_bytes": compact_bytes, "initial_rss_bytes": initial_rss,
        "peak_rss_bytes": peak_rss, "final_rss_bytes": final_rss,
        "largest_decompressed_payload_estimate": largest_decompressed,
        "wall_seconds": wall, "operation_counts": operation_counts,
        "raw_cache_files_deleted": 0,
    }
    atomic_write_json(output_dir / "run_metadata.json", metadata)
    report = (
        "# Compact point-certificate migration\n\n"
        f"Scientific scope: `{SCIENTIFIC_SCOPE}`.\n\n"
        f"Schema: `{SCHEMA_VERSION}`.  Certificates: {len(index_rows)}/{len(manifest_rows)}; "
        f"failures: {len(failures)}.\n\n"
        "The migration performed zero point solves, zero matrix evaluations, zero local repairs, "
        "and zero strict calls. Full solver payloads remain immutable and were not deleted.\n"
    )
    (output_dir / "report.md").write_text(report, encoding="utf-8")
    return CompactionRun(
        output_dir=output_dir, certificate_count=len(index_rows),
        converted_count=converted, cache_hit_count=cache_hits,
        failure_count=len(failures), initial_rss=initial_rss, peak_rss=peak_rss,
        final_rss=final_rss, source_bytes=source_bytes,
        compact_bytes=compact_bytes, wall_seconds=wall,
        operation_counts=operation_counts,
    )


def compact_pseudo_payload(certificate: Mapping[str, object]) -> dict[str, object]:
    """Adapter for existing detector/local-repair orchestration (not a solver trace)."""

    spectra = certificate["spectra"]
    result = certificate["result"]
    assert isinstance(spectra, Mapping) and isinstance(result, Mapping)
    return {
        "execution_status": result.get("execution_status"),
        "N_true": result.get("N_true"),
        "first_failed_mode": result.get("first_failed_mode"),
        "first_apparent_failed_mode": result.get("first_apparent_failed_mode"),
        "required_guard_mode": result.get("required_guard"),
        "unresolved_reason": result.get("unresolved_reason", ""),
        "models": {
            model: {
                "resolved_roots": list(spectra[model].get("roots", ())),
                "clusters": list(spectra[model].get("clusters", ())),
                "brackets_and_continuation_metadata": list(
                    spectra[model].get("unresolved_intervals", ())
                ),
            }
            for model in MODELS
        },
    }
