"""Diagnostic sequential beta continuation with mandatory global inventories.

This module contains no rod equations, determinant, material model, or boundary
condition.  Callers provide the legacy global search, interval scan, local
spectrum finalizer, and frequency accessors.  The helper only coordinates
sorted-spectrum predictors, close-root clusters, global anchors, fallback,
bounded exact transfer caching, counters, and atomic family checkpoints.

The result at every beta is an independently sorted spectrum.  Predictor
indices are numerical window hints and are not modal-descendant identities.
This is an internal diagnostic helper, not a production API.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import time
from collections import OrderedDict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Callable, Generic, Mapping, Sequence, TypeVar

import numpy as np
from numpy.typing import NDArray


FAST_SOLVER_VERSION = "yartsev-ch2-fast-beta-v1"
MANDATORY_GLOBAL_ANCHORS_DEG = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)
ANCHOR_RELATIVE_TOLERANCE = 1.0e-8
CLUSTER_RELATIVE_GAP = 1.0e-3
DUPLICATE_RELATIVE_TOLERANCE = 2.0e-8

RootT = TypeVar("RootT")
SpectrumT = TypeVar("SpectrumT")


class FastSolverValidationError(RuntimeError):
    """Raised when a mandatory fast/global validation contract fails."""


@dataclass(frozen=True)
class FastSweepSettings:
    root_count: int = 7
    mandatory_anchors_deg: tuple[float, ...] = MANDATORY_GLOBAL_ANCHORS_DEG
    anchor_relative_tolerance: float = ANCHOR_RELATIVE_TOLERANCE
    cluster_relative_gap: float = CLUSTER_RELATIVE_GAP
    duplicate_relative_tolerance: float = DUPLICATE_RELATIVE_TOLERANCE
    minimum_window_hz: float = 4.0
    relative_window_half_width: float = 3.0e-3
    motion_window_factor: float = 1.25
    predictor_residual_limit: float = 2.0e-2
    cache_maxsize: int = 65536
    solver_version: str = FAST_SOLVER_VERSION


@dataclass
class PerformanceCounters:
    transfer_expm_evaluations: int = 0
    transfer_cache_hits: int = 0
    transfer_cache_misses: int = 0
    boundary_matrix_assemblies: int = 0
    scaled_quality_evaluations: int = 0
    raw_quality_evaluations: int = 0
    local_root_solves: int = 0
    cluster_local_scans: int = 0
    global_inventory_checks: int = 0
    global_anchor_scans: int = 0
    global_fallback_scans: int = 0
    predictor_failures: int = 0
    family_runtime_s: float = 0.0
    total_scientific_runtime_s: float = 0.0

    def copy(self) -> "PerformanceCounters":
        return PerformanceCounters(**asdict(self))

    def difference(self, earlier: "PerformanceCounters") -> "PerformanceCounters":
        result: dict[str, int | float] = {}
        for name, value in asdict(self).items():
            result[name] = value - getattr(earlier, name)
        return PerformanceCounters(**result)

    @property
    def transfer_cache_requests(self) -> int:
        return self.transfer_cache_hits + self.transfer_cache_misses

    @property
    def transfer_cache_hit_rate(self) -> float:
        requests = self.transfer_cache_requests
        return self.transfer_cache_hits / requests if requests else 0.0

    def to_dict(self) -> dict[str, int | float]:
        result = asdict(self)
        result["transfer_cache_requests"] = self.transfer_cache_requests
        result["transfer_cache_hit_rate"] = self.transfer_cache_hit_rate
        return result


def _fingerprint_value(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            return repr(value)
        return {"float_hex": value.hex()}
    if isinstance(value, complex):
        return {
            "complex_real_hex": float(value.real).hex(),
            "complex_imag_hex": float(value.imag).hex(),
        }
    if isinstance(value, Mapping):
        return {
            str(key): _fingerprint_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_fingerprint_value(item) for item in value]
    if hasattr(value, "__dataclass_fields__"):
        return {
            "dataclass": f"{type(value).__module__}.{type(value).__qualname__}",
            "fields": _fingerprint_value(asdict(value)),
        }
    return {
        "type": f"{type(value).__module__}.{type(value).__qualname__}",
        "repr": repr(value),
    }


def stable_fingerprint(payload: Any) -> str:
    """Return a deterministic SHA-256 fingerprint without decimal rounding."""

    encoded = json.dumps(
        _fingerprint_value(payload),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


class ExactTransferLRU:
    """Bounded LRU for exact transfer-matrix calls.

    Keys contain the model type, the exact real/imaginary IEEE-754 values of
    omega, and a complete immutable point fingerprint.  Omega is never rounded.
    """

    def __init__(
        self,
        maxsize: int = 65536,
        *,
        counters: PerformanceCounters | None = None,
    ) -> None:
        if maxsize < 1:
            raise ValueError("cache maxsize must be positive")
        self.maxsize = int(maxsize)
        self.counters = counters or PerformanceCounters()
        self._values: OrderedDict[
            tuple[str, str, str, str], NDArray[np.complex128]
        ] = OrderedDict()
        self._point_fingerprints: dict[int, str] = {}

    def point_fingerprint(self, point: Any) -> str:
        identity = id(point)
        if identity not in self._point_fingerprints:
            self._point_fingerprints[identity] = stable_fingerprint(point)
        return self._point_fingerprints[identity]

    @staticmethod
    def _omega_key(omega: complex) -> tuple[str, str]:
        exact = complex(omega)
        return float(exact.real).hex(), float(exact.imag).hex()

    def get(
        self,
        model_type: str,
        omega: complex,
        point: Any,
        builder: Callable[[], NDArray[np.complex128]],
    ) -> NDArray[np.complex128]:
        real_hex, imag_hex = self._omega_key(omega)
        key = (str(model_type), real_hex, imag_hex, self.point_fingerprint(point))
        if key in self._values:
            self.counters.transfer_cache_hits += 1
            value = self._values.pop(key)
            self._values[key] = value
            return value
        self.counters.transfer_cache_misses += 1
        self.counters.transfer_expm_evaluations += 1
        value = np.asarray(builder(), dtype=np.complex128)
        self._values[key] = value
        if len(self._values) > self.maxsize:
            self._values.popitem(last=False)
        return value

    def __len__(self) -> int:
        return len(self._values)

    def clear(self) -> None:
        self._values.clear()


@dataclass(frozen=True)
class SearchInterval:
    lower_hz: float
    upper_hz: float
    indices: tuple[int, ...]
    predicted_hz: tuple[float, ...]
    is_cluster: bool

    @property
    def expected_count(self) -> int:
        return len(self.indices)


@dataclass
class FastSweepResult(Generic[SpectrumT]):
    spectra: dict[float, SpectrumT]
    data_origins: dict[float, str]
    fallback_reasons: dict[float, str]
    anchor_relative_errors: dict[float, float]
    counters: PerformanceCounters
    runtime_seconds: float
    status: str = "PASS"


def predict_sorted_frequencies(
    beta_deg: float,
    previous_beta_deg: float,
    previous_frequencies_hz: Sequence[float],
    *,
    older_beta_deg: float | None = None,
    older_frequencies_hz: Sequence[float] | None = None,
) -> NDArray[np.float64]:
    """Return the one-point or linear two-point numerical predictor."""

    previous = np.asarray(previous_frequencies_hz, dtype=float)
    if older_beta_deg is None or older_frequencies_hz is None:
        return previous.copy()
    older = np.asarray(older_frequencies_hz, dtype=float)
    denominator = previous_beta_deg - older_beta_deg
    if denominator <= 0.0:
        raise ValueError("predictor beta values must be strictly increasing")
    factor = (float(beta_deg) - previous_beta_deg) / denominator
    return previous + factor * (previous - older)


def connected_clusters(
    previous_frequencies_hz: Sequence[float],
    predicted_frequencies_hz: Sequence[float],
    *,
    threshold: float = CLUSTER_RELATIVE_GAP,
) -> tuple[tuple[int, ...], ...]:
    """Return connected close-root index sets from previous or predicted gaps."""

    previous = np.asarray(previous_frequencies_hz, dtype=float)
    predicted = np.asarray(predicted_frequencies_hz, dtype=float)
    if previous.shape != predicted.shape:
        raise ValueError("previous and predicted spectra must have equal shape")
    if previous.size < 2:
        return ()
    previous_gap = np.diff(previous) / np.maximum(previous[:-1], previous[1:])
    predicted_gap = np.diff(predicted) / np.maximum(predicted[:-1], predicted[1:])
    connected_edges = np.logical_or(previous_gap < threshold, predicted_gap < threshold)
    clusters: list[tuple[int, ...]] = []
    start: int | None = None
    for edge, is_connected in enumerate(connected_edges):
        if is_connected and start is None:
            start = edge
        if start is not None and (not is_connected or edge == len(connected_edges) - 1):
            stop = edge + 1 if is_connected and edge == len(connected_edges) - 1 else edge
            clusters.append(tuple(range(start, stop + 1)))
            start = None
    return tuple(clusters)


def _cluster_edges(
    frequencies_hz: Sequence[float], threshold: float
) -> tuple[bool, ...]:
    frequencies = np.asarray(frequencies_hz, dtype=float)
    if frequencies.size < 2:
        return ()
    gaps = np.diff(frequencies) / np.maximum(frequencies[:-1], frequencies[1:])
    return tuple(bool(value < threshold) for value in gaps)


def build_search_intervals(
    previous_frequencies_hz: Sequence[float],
    predicted_frequencies_hz: Sequence[float],
    *,
    older_frequencies_hz: Sequence[float] | None,
    settings: FastSweepSettings,
) -> list[SearchInterval]:
    """Build non-overlapping singleton/cluster search windows."""

    previous = np.asarray(previous_frequencies_hz, dtype=float)
    predicted = np.asarray(predicted_frequencies_hz, dtype=float)
    older = previous if older_frequencies_hz is None else np.asarray(
        older_frequencies_hz, dtype=float
    )
    if previous.shape != predicted.shape or previous.shape != older.shape:
        raise ValueError("window spectra must have equal shape")
    if np.any(~np.isfinite(predicted)) or np.any(predicted <= 0.0):
        raise ValueError("predictor produced non-positive or non-finite roots")
    if np.any(np.diff(predicted) <= 0.0):
        raise ValueError("predictor ordering crossed or became duplicate")

    clusters = connected_clusters(
        previous, predicted, threshold=settings.cluster_relative_gap
    )
    cluster_by_index = {index: cluster for cluster in clusters for index in cluster}
    groups: list[tuple[int, ...]] = []
    index = 0
    while index < len(predicted):
        group = cluster_by_index.get(index, (index,))
        groups.append(group)
        index = group[-1] + 1

    raw_windows: list[SearchInterval] = []
    for group in groups:
        motion = max(
            max(abs(predicted[item] - previous[item]) for item in group),
            max(abs(previous[item] - older[item]) for item in group),
        )
        half_width = max(
            settings.minimum_window_hz,
            settings.relative_window_half_width
            * max(abs(predicted[item]) for item in group),
            settings.motion_window_factor * motion,
        )
        lower = max(1.0e-6, float(min(predicted[item] for item in group) - half_width))
        upper = float(max(predicted[item] for item in group) + half_width)
        raw_windows.append(
            SearchInterval(
                lower_hz=lower,
                upper_hz=upper,
                indices=group,
                predicted_hz=tuple(float(predicted[item]) for item in group),
                is_cluster=len(group) > 1,
            )
        )

    # Partition any generous raw overlap at predictor midpoints.  These are
    # still local intervals, but overlap is impossible by construction.  A
    # predictor crossing or cluster-topology change is rejected before this
    # partition and therefore still forces a global inventory.
    partitioned: list[SearchInterval] = []
    for index, interval in enumerate(raw_windows):
        lower = interval.lower_hz
        upper = interval.upper_hz
        if index:
            separator = 0.5 * (
                raw_windows[index - 1].predicted_hz[-1]
                + interval.predicted_hz[0]
            )
            lower = max(lower, float(separator))
        if index + 1 < len(raw_windows):
            separator = 0.5 * (
                interval.predicted_hz[-1]
                + raw_windows[index + 1].predicted_hz[0]
            )
            upper = min(upper, float(separator))
        if lower >= upper:
            raise ValueError("partitioned local search window is empty")
        partitioned.append(
            SearchInterval(
                lower_hz=lower,
                upper_hz=upper,
                indices=interval.indices,
                predicted_hz=interval.predicted_hz,
                is_cluster=interval.is_cluster,
            )
        )
    return partitioned


def maximum_relative_frequency_error(
    left_hz: Sequence[float], right_hz: Sequence[float]
) -> float:
    left = np.asarray(left_hz, dtype=float)
    right = np.asarray(right_hz, dtype=float)
    if left.shape != right.shape:
        return math.inf
    return float(np.max(np.abs(left - right) / np.maximum(np.abs(right), 1.0)))


def _validate_numeric_inventory(
    frequencies_hz: Sequence[float], settings: FastSweepSettings
) -> None:
    frequencies = np.asarray(frequencies_hz, dtype=float)
    if len(frequencies) != settings.root_count:
        raise ValueError(
            f"local root count {len(frequencies)} != {settings.root_count}"
        )
    if np.any(~np.isfinite(frequencies)) or np.any(frequencies <= 0.0):
        raise ValueError("local spectrum contains invalid positive roots")
    if np.any(np.diff(frequencies) <= 0.0):
        raise ValueError("local spectrum is unsorted or contains duplicates")
    gaps = np.diff(frequencies) / np.maximum(frequencies[:-1], frequencies[1:])
    if np.any(gaps <= settings.duplicate_relative_tolerance):
        raise ValueError("local spectrum contains duplicate/indistinguishable roots")


def run_fast_beta_sweep(
    beta_values_deg: Sequence[float],
    *,
    global_search: Callable[[float], SpectrumT],
    scan_interval: Callable[[float, SearchInterval], Sequence[RootT]],
    finalize_local: Callable[[float, Sequence[RootT]], SpectrumT],
    spectrum_frequencies: Callable[[SpectrumT], Sequence[float]],
    root_frequency: Callable[[RootT], float],
    fallback_search: Callable[[float, SpectrumT, SpectrumT | None], SpectrumT]
    | None = None,
    settings: FastSweepSettings = FastSweepSettings(),
    counters: PerformanceCounters | None = None,
) -> FastSweepResult[SpectrumT]:
    """Run sequential sorted-spectrum continuation with mandatory fallback."""

    beta_values = np.asarray(beta_values_deg, dtype=float)
    if beta_values.ndim != 1 or beta_values.size < 1:
        raise ValueError("beta grid must be a non-empty one-dimensional sequence")
    if np.any(np.diff(beta_values) <= 0.0):
        raise ValueError("beta grid must be strictly increasing")
    counts = counters or PerformanceCounters()
    started = time.perf_counter()
    spectra: dict[float, SpectrumT] = {}
    origins: dict[float, str] = {}
    fallbacks: dict[float, str] = {}
    anchor_errors: dict[float, float] = {}
    anchor_set = {round(float(value), 12) for value in settings.mandatory_anchors_deg}

    for index, beta_value in enumerate(beta_values):
        beta = float(beta_value)
        is_anchor = round(beta, 12) in anchor_set
        if index == 0:
            counts.global_anchor_scans += 1
            spectrum = global_search(beta)
            frequencies = np.asarray(spectrum_frequencies(spectrum), dtype=float)
            _validate_numeric_inventory(frequencies, settings)
            spectra[beta] = spectrum
            origins[beta] = "global_anchor"
            anchor_errors[beta] = 0.0
            continue

        previous_beta = float(beta_values[index - 1])
        previous_spectrum = spectra[previous_beta]
        previous = np.asarray(spectrum_frequencies(previous_spectrum), dtype=float)
        older_beta = float(beta_values[index - 2]) if index >= 2 else None
        older = (
            np.asarray(spectrum_frequencies(spectra[older_beta]), dtype=float)
            if older_beta is not None
            else None
        )
        local_spectrum: SpectrumT | None = None
        local_failure: str | None = None
        try:
            predicted = predict_sorted_frequencies(
                beta,
                previous_beta,
                previous,
                older_beta_deg=older_beta,
                older_frequencies_hz=older,
            )
            if np.any(np.diff(predicted) <= 0.0):
                raise ValueError("predicted sorted ordering crossed")
            if older is not None and _cluster_edges(
                previous, settings.cluster_relative_gap
            ) != _cluster_edges(predicted, settings.cluster_relative_gap):
                raise ValueError("cluster topology changed")
            intervals = build_search_intervals(
                previous,
                predicted,
                older_frequencies_hz=older,
                settings=settings,
            )
            candidates: list[RootT] = []
            for interval in intervals:
                counts.local_root_solves += 1
                if interval.is_cluster:
                    counts.cluster_local_scans += 1
                interval_roots = list(scan_interval(beta, interval))
                if len(interval_roots) != interval.expected_count:
                    raise ValueError(
                        "cluster/local interval root count mismatch: "
                        f"expected {interval.expected_count}, got {len(interval_roots)}"
                    )
                for root in interval_roots:
                    frequency = float(root_frequency(root))
                    if not (
                        interval.lower_hz * (1.0 - 1.0e-12)
                        <= frequency
                        <= interval.upper_hz * (1.0 + 1.0e-12)
                    ):
                        raise ValueError("local root escaped its search interval")
                candidates.extend(interval_roots)
            candidates.sort(key=root_frequency)
            frequencies = np.asarray(
                [float(root_frequency(root)) for root in candidates], dtype=float
            )
            _validate_numeric_inventory(frequencies, settings)
            predictor_residual = maximum_relative_frequency_error(
                frequencies, predicted
            )
            if predictor_residual > settings.predictor_residual_limit:
                raise ValueError(
                    "sharp predictor residual "
                    f"{predictor_residual:.6e} > {settings.predictor_residual_limit:.6e}"
                )
            local_spectrum = finalize_local(beta, candidates)
            _validate_numeric_inventory(
                spectrum_frequencies(local_spectrum), settings
            )
        except Exception as error:  # fallback is an explicit part of the algorithm
            counts.predictor_failures += 1
            local_failure = f"{type(error).__name__}: {error}"

        if is_anchor:
            counts.global_anchor_scans += 1
            global_spectrum = global_search(beta)
            global_frequencies = np.asarray(
                spectrum_frequencies(global_spectrum), dtype=float
            )
            _validate_numeric_inventory(global_frequencies, settings)
            if local_spectrum is None:
                counts.global_fallback_scans += 1
                fallbacks[beta] = local_failure or "local path unavailable at anchor"
                spectrum = (
                    fallback_search(
                        beta,
                        previous_spectrum,
                        spectra[older_beta] if older_beta is not None else None,
                    )
                    if fallback_search is not None
                    else global_spectrum
                )
                _validate_numeric_inventory(
                    spectrum_frequencies(spectrum), settings
                )
                origins[beta] = "global_fallback"
                anchor_errors[beta] = 0.0
            else:
                error = maximum_relative_frequency_error(
                    spectrum_frequencies(local_spectrum), global_frequencies
                )
                anchor_errors[beta] = error
                if error > settings.anchor_relative_tolerance:
                    raise FastSolverValidationError(
                        "FAST_SOLVER_VALIDATION_FAIL: anchor mismatch at "
                        f"beta={beta:g} deg: {error:.6e} > "
                        f"{settings.anchor_relative_tolerance:.6e}"
                    )
                spectrum = global_spectrum
                origins[beta] = "global_anchor"
        elif local_spectrum is None:
            counts.global_fallback_scans += 1
            fallbacks[beta] = local_failure or "local path unavailable"
            spectrum = (
                fallback_search(
                    beta,
                    previous_spectrum,
                    spectra[older_beta] if older_beta is not None else None,
                )
                if fallback_search is not None
                else global_search(beta)
            )
            _validate_numeric_inventory(spectrum_frequencies(spectrum), settings)
            origins[beta] = "global_fallback"
        else:
            spectrum = local_spectrum
            origins[beta] = "new_fast_solve"

        spectra[beta] = spectrum

    elapsed = time.perf_counter() - started
    counts.family_runtime_s += elapsed
    counts.total_scientific_runtime_s += elapsed
    return FastSweepResult(
        spectra=spectra,
        data_origins=origins,
        fallback_reasons=fallbacks,
        anchor_relative_errors=anchor_errors,
        counters=counts,
        runtime_seconds=elapsed,
    )


def checkpoint_paths(directory: Path, family_id: str) -> tuple[Path, Path]:
    safe = "".join(character if character.isalnum() or character in "-_" else "_" for character in family_id)
    return directory / f"{safe}.csv", directory / f"{safe}.manifest.json"


def write_family_checkpoint(
    directory: Path,
    family_id: str,
    rows: Sequence[Mapping[str, Any]],
    *,
    fingerprint: str,
    metadata: Mapping[str, Any],
) -> tuple[Path, Path]:
    """Atomically write and validate one completed spectral family."""

    if not rows:
        raise ValueError("cannot checkpoint an empty family")
    directory.mkdir(parents=True, exist_ok=True)
    csv_path, manifest_path = checkpoint_paths(directory, family_id)
    csv_tmp = csv_path.with_suffix(csv_path.suffix + ".tmp")
    manifest_tmp = manifest_path.with_suffix(manifest_path.suffix + ".tmp")
    fieldnames = list(rows[0].keys())
    with csv_tmp.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            if list(row.keys()) != fieldnames:
                raise ValueError("checkpoint rows do not share one schema")
            writer.writerow(row)
    with csv_tmp.open("r", newline="", encoding="utf-8") as stream:
        validated_rows = list(csv.DictReader(stream))
    if len(validated_rows) != len(rows):
        raise RuntimeError("temporary family checkpoint row count mismatch")
    csv_sha256 = hashlib.sha256(csv_tmp.read_bytes()).hexdigest()
    payload = {
        "family_id": family_id,
        "fingerprint": fingerprint,
        "row_count": len(rows),
        "csv_sha256": csv_sha256,
        "status": "COMPLETE",
        **dict(metadata),
    }
    manifest_tmp.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(csv_tmp, csv_path)
    os.replace(manifest_tmp, manifest_path)
    return csv_path, manifest_path


def load_family_checkpoint(
    directory: Path,
    family_id: str,
    *,
    expected_fingerprint: str,
    expected_row_count: int,
) -> tuple[list[dict[str, str]], dict[str, Any]] | None:
    """Load one complete checkpoint only when its full fingerprint matches."""

    csv_path, manifest_path = checkpoint_paths(directory, family_id)
    if not csv_path.is_file() or not manifest_path.is_file():
        return None
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (
        manifest.get("status") != "COMPLETE"
        or manifest.get("fingerprint") != expected_fingerprint
        or int(manifest.get("row_count", -1)) != expected_row_count
    ):
        return None
    if hashlib.sha256(csv_path.read_bytes()).hexdigest() != manifest.get("csv_sha256"):
        return None
    with csv_path.open("r", newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != expected_row_count:
        return None
    return rows, manifest


__all__ = [
    "ANCHOR_RELATIVE_TOLERANCE",
    "CLUSTER_RELATIVE_GAP",
    "ExactTransferLRU",
    "FAST_SOLVER_VERSION",
    "FastSolverValidationError",
    "FastSweepResult",
    "FastSweepSettings",
    "MANDATORY_GLOBAL_ANCHORS_DEG",
    "PerformanceCounters",
    "SearchInterval",
    "build_search_intervals",
    "checkpoint_paths",
    "connected_clusters",
    "load_family_checkpoint",
    "maximum_relative_frequency_error",
    "predict_sorted_frequencies",
    "run_fast_beta_sweep",
    "stable_fingerprint",
    "write_family_checkpoint",
]
