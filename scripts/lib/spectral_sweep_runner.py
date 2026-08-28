"""Model-agnostic continuation runner for sorted spectral sweeps.

The runner contains no beam equations, determinant, boundary condition, or
root refiner.  Model adapters provide global and local searches together with
the frozen acceptance evidence.  Continuation is used only to locate local
search intervals: every exported value must be returned by a supplied root
search callback.

The integer position of a root is an independently sorted numerical locator.
It is not a modal ``branch_id`` and it does not establish modal ancestry.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import time
from collections import OrderedDict
from collections.abc import Callable, Iterator, Mapping, MutableMapping, Sequence
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Any, Generic, TypeVar

import numpy as np
from numpy.typing import NDArray


RUNNER_VERSION = "spectral-sweep-runner-v2"
PLOT_COUNT = 8
GUARD_POSITION = 9

ValueT = TypeVar("ValueT")


class SweepConfigurationError(ValueError):
    """Raised when a sweep/checkpoint contract is internally inconsistent."""


class SweepValidationError(RuntimeError):
    """Raised when neither a local solve nor the mandatory fallback is valid."""


@dataclass(frozen=True)
class RootRecord:
    """One accepted root and model-owned quality evidence.

    ``value`` is the canonical positive spectral coordinate selected by the
    adapter.  The generic runner never transforms or interpolates it.
    """

    value: float
    accepted: bool = True
    sigma_ratio: float = 0.0
    boundary_residual: float = 0.0
    detector_agreement: bool = True
    possible_even_root: bool = False
    multiplicity: int = 1
    nullity: int = 1
    cluster_id: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class SpectrumRecord:
    """Exactly positions 1--8 and root 9 for one parameter point."""

    parameter: float
    roots: tuple[RootRecord, ...]
    origin: str
    status: str = "PASS"
    internal_tail_status: str = "NOT_COMPUTED"
    metadata: Mapping[str, Any] = field(default_factory=dict)

    @property
    def values(self) -> NDArray[np.float64]:
        return np.asarray([root.value for root in self.roots], dtype=float)


@dataclass(frozen=True)
class SearchInterval:
    """Non-overlapping interval for one root or one multi-root locator group.

    A multi-root locator group is numerical search ownership only.  It does
    not assert physical multiplicity, a modal branch, or a merged root.
    """

    lower: float
    upper: float
    positions: tuple[int, ...]
    predicted: tuple[float, ...]
    is_cluster: bool = False
    verification: bool = False
    partition_lower: float | None = None
    partition_upper: float | None = None

    @property
    def expected_count(self) -> int:
        return len(self.positions)


@dataclass(frozen=True)
class SweepSettings:
    """Frozen numerical orchestration settings, independent of physics."""

    plot_count: int = PLOT_COUNT
    guard_position: int = GUARD_POSITION
    anchor_stride: int = 10
    local_verification_relative: float = 1.0e-9
    anchor_comparison_relative: float = 1.0e-9
    locator_relative_limit: float = 5.0e-2
    minimum_window: float = 1.0e-8
    spectral_minimum: float = 1.0e-12
    relative_window_half_width: float = 2.5e-2
    motion_window_factor: float = 2.0
    verification_expansion_factor: float = 1.25
    cluster_relative_gap: float = 1.0e-6
    fallback_gap_relative: float = 1.0e-7
    duplicate_relative_tolerance: float = 1.0e-12
    spike_relative_threshold: float = 5.0e-3
    spike_neighbor_factor: float = 5.0
    spike_absolute_floor: float = 1.0e-12
    enable_spike_audit: bool = True
    force_anchor_after_fallback: bool = True
    runner_version: str = RUNNER_VERSION

    def __post_init__(self) -> None:
        if self.plot_count != 8 or self.guard_position != 9:
            raise SweepConfigurationError(
                "The production contract is positions 1--8 plus root 9."
            )
        if self.anchor_stride < 1:
            raise SweepConfigurationError("anchor_stride must be positive")
        positive = (
            self.local_verification_relative,
            self.anchor_comparison_relative,
            self.locator_relative_limit,
            self.minimum_window,
            self.spectral_minimum,
            self.relative_window_half_width,
            self.motion_window_factor,
            self.verification_expansion_factor,
            self.cluster_relative_gap,
            self.fallback_gap_relative,
            self.duplicate_relative_tolerance,
            self.spike_relative_threshold,
            self.spike_neighbor_factor,
            self.spike_absolute_floor,
        )
        if any(not math.isfinite(value) or value <= 0.0 for value in positive):
            raise SweepConfigurationError("runner thresholds must be positive")
        if self.verification_expansion_factor <= 1.0:
            raise SweepConfigurationError(
                "verification_expansion_factor must exceed one"
            )


@dataclass
class SweepCounters:
    points_requested: int = 0
    points_committed: int = 0
    points_resumed: int = 0
    local_intervals: int = 0
    cluster_intervals: int = 0
    local_primary_searches: int = 0
    local_verification_searches: int = 0
    global_full_scans: int = 0
    global_anchor_scans: int = 0
    global_fallback_scans: int = 0
    fallback_count: int = 0
    predictor_failures: int = 0
    spike_audit_count: int = 0
    spike_full_scans: int = 0
    checkpoint_transactions: int = 0
    checkpoint_rejects: int = 0
    elapsed_time: float = 0.0

    def to_dict(self) -> dict[str, int | float]:
        return asdict(self)


@dataclass(frozen=True)
class PointAudit:
    parameter: float
    point_index: int
    scheduled_anchor: bool
    origin: str
    predictor_kind: str
    fallback_reasons: tuple[str, ...] = ()
    maximum_primary_verification_relative: float = 0.0
    maximum_predictor_relative: float = 0.0
    minimum_relative_gap: float = math.inf


@dataclass(frozen=True)
class SpikeAuditRecord:
    point_index: int
    parameter: float
    position: int
    fast_value: float
    neighbor_prediction: float
    normalized_residual: float
    outcome: str
    full_value: float | None = None
    triplet_parameters: tuple[float, float, float] | None = None
    corrected_parameters: tuple[float, ...] = ()
    full_triplet_max_relative: float | None = None


@dataclass(frozen=True)
class SweepCallbacks:
    """Model-owned root operations used by :func:`run_spectral_sweep`."""

    global_search: Callable[[float], SpectrumRecord]
    local_search: Callable[[float, SearchInterval, bool], Sequence[RootRecord]]
    quality_gate: Callable[[RootRecord], tuple[bool, str]] | None = None
    completeness_guard: Callable[
        [float, Sequence[RootRecord]], tuple[bool, str]
    ] | None = None
    global_completeness_guard: Callable[
        [float, Sequence[RootRecord]], tuple[bool, str]
    ] | None = None


@dataclass(frozen=True)
class CheckpointConfig:
    directory: Path
    sweep_id: str
    model_id: str
    fingerprint: str


@dataclass
class SpectralSweepResult:
    spectra: dict[float, SpectrumRecord]
    point_audits: list[PointAudit]
    spike_audits: list[SpikeAuditRecord]
    counters: SweepCounters
    status: str
    runtime_seconds: float
    failed_points: tuple[float, ...] = ()
    qualified_points: tuple[float, ...] = ()


def _fingerprint_value(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return {"float_hex": value.hex()} if math.isfinite(value) else repr(value)
    if isinstance(value, Mapping):
        return {
            str(key): _fingerprint_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_fingerprint_value(item) for item in value]
    if hasattr(value, "__dataclass_fields__"):
        return _fingerprint_value(asdict(value))
    return {"type": f"{type(value).__module__}.{type(value).__qualname__}", "repr": repr(value)}


def stable_fingerprint(payload: Any) -> str:
    """Return a deterministic SHA-256 without decimal rounding."""

    encoded = json.dumps(
        _fingerprint_value(payload),
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest().upper()


def _estimated_bytes(value: Any) -> int:
    if isinstance(value, np.ndarray):
        return int(value.nbytes)
    if isinstance(value, (bytes, bytearray, str)):
        return len(value)
    if isinstance(value, Mapping):
        return sum(_estimated_bytes(key) + _estimated_bytes(item) for key, item in value.items())
    if isinstance(value, (tuple, list)):
        return sum(_estimated_bytes(item) for item in value)
    return 64


class ExactLRU(MutableMapping[Any, ValueT], Generic[ValueT]):
    """Exact-key LRU bounded by both entry count and estimated payload bytes."""

    def __init__(self, max_entries: int = 4096, max_bytes: int = 256 * 1024**2) -> None:
        if max_entries < 1 or max_bytes < 1:
            raise ValueError("cache bounds must be positive")
        self.max_entries = int(max_entries)
        self.max_bytes = int(max_bytes)
        self._values: OrderedDict[Any, tuple[ValueT, int]] = OrderedDict()
        self.current_bytes = 0
        self.peak_entries = 0
        self.peak_bytes = 0
        self.hits = 0
        self.misses = 0
        self.evictions = 0

    def __getitem__(self, key: Any) -> ValueT:
        if key not in self._values:
            self.misses += 1
            raise KeyError(key)
        self.hits += 1
        value, size = self._values.pop(key)
        self._values[key] = (value, size)
        return value

    def __setitem__(self, key: Any, value: ValueT) -> None:
        size = _estimated_bytes(value)
        if key in self._values:
            _, old_size = self._values.pop(key)
            self.current_bytes -= old_size
        self._values[key] = (value, size)
        self.current_bytes += size
        while len(self._values) > self.max_entries or self.current_bytes > self.max_bytes:
            _, (_, removed_size) = self._values.popitem(last=False)
            self.current_bytes -= removed_size
            self.evictions += 1
        self.peak_entries = max(self.peak_entries, len(self._values))
        self.peak_bytes = max(self.peak_bytes, self.current_bytes)

    def __delitem__(self, key: Any) -> None:
        _, size = self._values.pop(key)
        self.current_bytes -= size

    def __iter__(self) -> Iterator[Any]:
        return iter(self._values)

    def __len__(self) -> int:
        return len(self._values)

    def clear(self) -> None:
        self._values.clear()
        self.current_bytes = 0

    def get_or_build(self, key: Any, builder: Callable[[], ValueT]) -> ValueT:
        try:
            return self[key]
        except KeyError:
            value = builder()
            self[key] = value
            return value

    def diagnostics(self) -> dict[str, int | float]:
        requests = self.hits + self.misses
        return {
            "entries": len(self),
            "current_bytes": self.current_bytes,
            "peak_entries": self.peak_entries,
            "peak_bytes": self.peak_bytes,
            "hits": self.hits,
            "misses": self.misses,
            "evictions": self.evictions,
            "hit_rate": self.hits / requests if requests else 0.0,
        }


def predict_sorted_roots(
    parameter: float,
    previous_parameter: float,
    previous_roots: Sequence[float],
    *,
    older_parameter: float | None = None,
    older_roots: Sequence[float] | None = None,
) -> NDArray[np.float64]:
    """Return hold or nonuniform-grid secant predictions as locators only."""

    previous = np.asarray(previous_roots, dtype=float)
    if older_parameter is None or older_roots is None:
        return previous.copy()
    older = np.asarray(older_roots, dtype=float)
    if previous.shape != older.shape:
        raise SweepConfigurationError("predictor spectra have different shapes")
    denominator = float(previous_parameter) - float(older_parameter)
    if not math.isfinite(denominator) or denominator <= 0.0:
        raise SweepConfigurationError("predictor parameters must increase")
    factor = (float(parameter) - float(previous_parameter)) / denominator
    return previous + factor * (previous - older)


def anchor_indices(point_count: int, stride: int) -> tuple[int, ...]:
    """Return endpoint-inclusive periodic anchor indices."""

    if point_count < 1 or stride < 1:
        raise SweepConfigurationError("invalid point_count or anchor stride")
    return tuple(sorted({0, point_count - 1, *range(0, point_count, stride)}))


def _relative_gaps(values: Sequence[float]) -> NDArray[np.float64]:
    array = np.asarray(values, dtype=float)
    if len(array) < 2:
        return np.asarray([], dtype=float)
    return np.diff(array) / np.maximum.reduce(
        (np.abs(array[:-1]), np.abs(array[1:]), np.full(len(array) - 1, np.finfo(float).tiny))
    )


def _cluster_groups(
    previous: NDArray[np.float64],
    predicted: NDArray[np.float64],
    threshold: float,
) -> list[tuple[int, ...]]:
    connected = np.logical_or(
        _relative_gaps(previous) < threshold,
        _relative_gaps(predicted) < threshold,
    )
    groups: list[tuple[int, ...]] = []
    start = 0
    for edge, joined in enumerate(connected):
        if not joined:
            groups.append(tuple(range(start, edge + 1)))
            start = edge + 1
    groups.append(tuple(range(start, len(predicted))))
    return groups


def build_search_intervals(
    previous_roots: Sequence[float],
    predicted_roots: Sequence[float],
    *,
    older_roots: Sequence[float] | None,
    settings: SweepSettings,
) -> list[SearchInterval]:
    """Build midpoint-partitioned singleton or close-pair local windows."""

    previous = np.asarray(previous_roots, dtype=float)
    predicted = np.asarray(predicted_roots, dtype=float)
    older = previous if older_roots is None else np.asarray(older_roots, dtype=float)
    expected_shape = (settings.guard_position,)
    if previous.shape != expected_shape or predicted.shape != expected_shape or older.shape != expected_shape:
        raise SweepConfigurationError("interval spectra must contain exactly root 1..9")
    if np.any(~np.isfinite(predicted)) or np.any(predicted <= 0.0):
        raise SweepConfigurationError("prediction is nonfinite or nonpositive")
    groups = _cluster_groups(previous, predicted, settings.cluster_relative_gap)

    def raw_interval(group: tuple[int, ...]) -> tuple[float, float, tuple[int, ...]]:
        motion = max(
            max(abs(predicted[i] - previous[i]) for i in group),
            max(abs(previous[i] - older[i]) for i in group),
        )
        scale = max(abs(predicted[i]) for i in group)
        half = max(
            settings.minimum_window,
            settings.relative_window_half_width * scale,
            settings.motion_window_factor * motion,
        )
        return (
            max(np.finfo(float).tiny, min(predicted[i] for i in group) - half),
            max(predicted[i] for i in group) + half,
            group,
        )

    # A predictor safety margin may overlap its neighbour even when the two
    # roots remain distinct and well resolved.  Search every connected set of
    # overlapping margins as one multi-root locator interval.  This changes
    # only the numerical search window: the callback must still return the
    # exact expected number of separate RootRecord objects, and all topology,
    # duplicate, quality, and gap gates remain active below.
    merged_raw: list[tuple[float, float, tuple[int, ...]]] = []
    for initial_group in groups:
        candidate = raw_interval(initial_group)
        while merged_raw and candidate[0] <= merged_raw[-1][1]:
            previous_group = merged_raw.pop()[2]
            candidate = raw_interval(previous_group + candidate[2])
        merged_raw.append(candidate)
    raw = merged_raw

    intervals: list[SearchInterval] = []
    for index, (raw_lower, raw_upper, group) in enumerate(raw):
        partition_lower = None if index == 0 else 0.5 * (
            max(predicted[i] for i in raw[index - 1][2])
            + min(predicted[i] for i in group)
        )
        partition_upper = None if index + 1 == len(raw) else 0.5 * (
            max(predicted[i] for i in group)
            + min(predicted[i] for i in raw[index + 1][2])
        )
        # Primary windows retain the predictor safety margin.  The midpoint
        # partition is stored separately as the hard non-overlap envelope for
        # the independently expanded verification search.
        lower = max(
            settings.spectral_minimum,
            raw_lower,
            partition_lower if partition_lower is not None else settings.spectral_minimum,
        )
        upper = min(
            raw_upper,
            partition_upper if partition_upper is not None else raw_upper,
        )
        if not lower < upper:
            raise SweepConfigurationError("local intervals overlap or are empty")
        intervals.append(SearchInterval(
            lower=float(lower),
            upper=float(upper),
            positions=tuple(i + 1 for i in group),
            predicted=tuple(float(predicted[i]) for i in group),
            is_cluster=len(group) > 1,
            partition_lower=(
                None if partition_lower is None else float(partition_lower)
            ),
            partition_upper=(
                None if partition_upper is None else float(partition_upper)
            ),
        ))
    for left, right in zip(intervals, intervals[1:]):
        if left.upper > right.lower:
            raise SweepConfigurationError("partitioned local intervals overlap")
    return intervals


def expanded_interval(interval: SearchInterval, settings: SweepSettings) -> SearchInterval:
    minimum_prediction = min(interval.predicted)
    maximum_prediction = max(interval.predicted)
    clearance = max(
        settings.minimum_window,
        min(minimum_prediction - interval.lower, interval.upper - maximum_prediction),
    )
    padding = (settings.verification_expansion_factor - 1.0) * clearance
    lower = max(settings.spectral_minimum, interval.lower - padding)
    upper = interval.upper + padding
    if interval.partition_lower is not None:
        lower = max(lower, interval.partition_lower)
    if interval.partition_upper is not None:
        upper = min(upper, interval.partition_upper)
    if not lower < upper:
        raise SweepConfigurationError("expanded local interval is empty")
    return replace(
        interval,
        lower=lower,
        upper=upper,
        verification=True,
    )


def maximum_relative_error(left: Sequence[float], right: Sequence[float]) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    if a.shape != b.shape or not a.size:
        return math.inf
    return float(np.max(np.abs(a - b) / np.maximum.reduce((np.abs(a), np.abs(b), np.full(a.shape, np.finfo(float).tiny)))))


def _root_topology_matches(left: RootRecord, right: RootRecord) -> bool:
    return bool(
        left.multiplicity == right.multiplicity
        and left.nullity == right.nullity
        and bool(left.cluster_id) == bool(right.cluster_id)
    )


def _topology_signature(roots: Sequence[RootRecord]) -> tuple[tuple[int, ...], ...]:
    """Describe the contiguous cluster/nullity partition without label names."""

    signature: list[tuple[int, ...]] = []
    index = 0
    while index < len(roots):
        root = roots[index]
        if not root.cluster_id:
            signature.append((1, root.multiplicity, root.nullity, 0))
            index += 1
            continue
        stop = index + 1
        while stop < len(roots) and roots[stop].cluster_id == root.cluster_id:
            stop += 1
        members = roots[index:stop]
        signature.append((
            len(members),
            *tuple(item.multiplicity for item in members),
            *tuple(item.nullity for item in members),
            1,
        ))
        index = stop
    return tuple(signature)


def _validate_roots(
    parameter: float,
    roots: Sequence[RootRecord],
    settings: SweepSettings,
    quality_gate: Callable[[RootRecord], tuple[bool, str]] | None,
) -> tuple[RootRecord, ...]:
    ordered = tuple(sorted(roots, key=lambda root: root.value))
    if len(ordered) != settings.guard_position:
        raise SweepValidationError(f"expected exactly root 1..9, got {len(ordered)}")
    values = np.asarray([root.value for root in ordered], dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise SweepValidationError("root inventory is nonfinite or nonpositive")
    if np.any(np.diff(values) < 0.0):
        raise SweepValidationError("root inventory is not sorted")
    for index, gap in enumerate(_relative_gaps(values)):
        if gap > settings.duplicate_relative_tolerance:
            continue
        left = ordered[index]
        right = ordered[index + 1]
        proven_multiplicity = bool(
            left.cluster_id
            and left.cluster_id == right.cluster_id
            and left.multiplicity > 1
            and right.multiplicity == left.multiplicity
        )
        if not proven_multiplicity:
            raise SweepValidationError("duplicate representation of one root")
    for position, root in enumerate(ordered, start=1):
        if not root.accepted:
            raise SweepValidationError(f"root {position} has rejected status")
        if not root.detector_agreement:
            raise SweepValidationError(f"root {position} detector disagreement")
        if root.possible_even_root:
            raise SweepValidationError(
                f"root {position} has unresolved even/sigma-only evidence"
            )
        if quality_gate is not None:
            passed, reason = quality_gate(root)
            if not passed:
                raise SweepValidationError(
                    f"root {position} failed frozen quality gate: {reason}"
                )
    return ordered


def _validate_spectrum(
    parameter: float,
    record: SpectrumRecord,
    callbacks: SweepCallbacks,
    settings: SweepSettings,
    origin: str,
) -> SpectrumRecord:
    if record.status not in {"PASS", "QUALIFIED", "PARTIAL_PASS_INTERNAL_TAIL"}:
        raise SweepValidationError(
            f"global inventory has rejected status {record.status!r}"
        )
    roots = _validate_roots(
        parameter,
        record.roots,
        settings,
        callbacks.quality_gate,
    )
    if callbacks.global_completeness_guard is not None:
        passed, reason = callbacks.global_completeness_guard(parameter, roots)
        if not passed:
            raise SweepValidationError(f"global completeness guard failed: {reason}")
    return _spectrum_with_origin(replace(record, roots=roots), parameter, origin)


def _record_payload(record: SpectrumRecord) -> dict[str, Any]:
    return {
        "parameter": record.parameter,
        "origin": record.origin,
        "status": record.status,
        "internal_tail_status": record.internal_tail_status,
        "metadata": dict(record.metadata),
        "roots": [
            {
                **{key: value for key, value in asdict(root).items() if key != "metadata"},
                "metadata": dict(root.metadata),
            }
            for root in record.roots
        ],
    }


def _record_from_payload(payload: Mapping[str, Any]) -> SpectrumRecord:
    return SpectrumRecord(
        parameter=float(payload["parameter"]),
        roots=tuple(RootRecord(**root) for root in payload["roots"]),
        origin=str(payload["origin"]),
        status=str(payload.get("status", "PASS")),
        internal_tail_status=str(payload.get("internal_tail_status", "NOT_COMPUTED")),
        metadata=dict(payload.get("metadata", {})),
    )


def _point_path(config: CheckpointConfig, point_index: int) -> Path:
    return Path(config.directory) / "points" / f"{point_index:06d}.json"


def write_point_transaction(
    config: CheckpointConfig,
    point_index: int,
    record: SpectrumRecord,
) -> Path:
    """Write one verified point through temporary-file validation and rename."""

    target = _point_path(config, point_index)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "runner_version": RUNNER_VERSION,
        "sweep_id": config.sweep_id,
        "model_id": config.model_id,
        "fingerprint": config.fingerprint,
        "point_index": int(point_index),
        "record": _record_payload(record),
    }
    payload["record_sha256"] = stable_fingerprint(payload["record"])
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    checked = json.loads(temporary.read_text(encoding="utf-8"))
    if checked != payload or stable_fingerprint(checked["record"]) != checked["record_sha256"]:
        raise SweepValidationError("temporary point transaction failed validation")
    os.replace(temporary, target)
    return target


def load_point_transaction(
    config: CheckpointConfig,
    point_index: int,
    expected_parameter: float,
) -> SpectrumRecord | None:
    target = _point_path(config, point_index)
    if not target.is_file():
        return None
    try:
        payload = json.loads(target.read_text(encoding="utf-8"))
        if (
            payload.get("fingerprint") != config.fingerprint
            or payload.get("sweep_id") != config.sweep_id
            or payload.get("model_id") != config.model_id
            or int(payload.get("point_index", -1)) != point_index
            or stable_fingerprint(payload["record"]) != payload.get("record_sha256")
        ):
            return None
        record = _record_from_payload(payload["record"])
        if float(record.parameter).hex() != float(expected_parameter).hex():
            return None
        return record
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None


def write_checkpoint_manifest(
    config: CheckpointConfig,
    parameters: Sequence[float],
    records: Mapping[float, SpectrumRecord],
    counters: SweepCounters,
    *,
    failed_points: Sequence[float] = (),
    qualified_points: Sequence[float] = (),
    estimated_remaining_time: float = 0.0,
) -> Path:
    completed = [float(value) for value in parameters if float(value) in records]
    missing = [float(value) for value in parameters if float(value) not in records]
    payload = {
        "schema_version": 1,
        "runner_version": RUNNER_VERSION,
        "sweep_id": config.sweep_id,
        "model_id": config.model_id,
        "fingerprint": config.fingerprint,
        "completed_points": completed,
        "missing_points": missing,
        "failed_points": [float(value) for value in failed_points],
        "qualified_points": [float(value) for value in qualified_points],
        "last_completed_parameter": completed[-1] if completed else None,
        "elapsed_time": counters.elapsed_time,
        "estimated_remaining_time": float(estimated_remaining_time),
        "full_scan_count": counters.global_full_scans,
        "fallback_count": counters.fallback_count,
        "spike_audit_count": counters.spike_audit_count,
        "counters": counters.to_dict(),
    }
    target = Path(config.directory) / "checkpoint_manifest.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    json.loads(temporary.read_text(encoding="utf-8"))
    os.replace(temporary, target)
    return target


def _require_checkpoint_identity(config: CheckpointConfig) -> None:
    """Reject reuse of a checkpoint produced under a different contract."""

    target = Path(config.directory) / "checkpoint_manifest.json"
    if not target.is_file():
        return
    try:
        payload = json.loads(target.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise SweepConfigurationError(
            "checkpoint manifest is unreadable; use an empty directory for a new contract"
        ) from error
    identity = (
        payload.get("sweep_id"),
        payload.get("model_id"),
        payload.get("fingerprint"),
    )
    expected = (config.sweep_id, config.model_id, config.fingerprint)
    if identity != expected:
        raise SweepConfigurationError(
            "checkpoint identity/fingerprint does not match the requested sweep"
        )


def _spectrum_with_origin(record: SpectrumRecord, parameter: float, origin: str) -> SpectrumRecord:
    return replace(record, parameter=float(parameter), origin=origin)


def detect_spikes(
    parameters: Sequence[float],
    spectra: Mapping[float, SpectrumRecord],
    settings: SweepSettings,
) -> list[SpikeAuditRecord]:
    """Flag isolated points; this function never changes spectral data."""

    values = np.asarray(parameters, dtype=float)
    records: list[SpikeAuditRecord] = []
    for index in range(1, len(values) - 1):
        left_p, point_p, right_p = map(float, values[index - 1 : index + 2])
        if any(parameter not in spectra for parameter in (left_p, point_p, right_p)):
            continue
        left = spectra[left_p].values
        point = spectra[point_p].values
        right = spectra[right_p].values
        factor = (point_p - left_p) / (right_p - left_p)
        prediction = left + factor * (right - left)
        residual = np.abs(point - prediction)
        scale = np.maximum.reduce((np.abs(point), np.abs(prediction), np.full(point.shape, settings.spike_absolute_floor)))
        normalized = residual / scale
        isolated = residual > np.maximum(
            np.abs(right - left) / settings.spike_neighbor_factor,
            settings.spike_absolute_floor,
        )
        for position in np.flatnonzero(
            np.logical_and(normalized > settings.spike_relative_threshold, isolated)
        ):
            records.append(SpikeAuditRecord(
                point_index=index,
                parameter=point_p,
                position=int(position + 1),
                fast_value=float(point[position]),
                neighbor_prediction=float(prediction[position]),
                normalized_residual=float(normalized[position]),
                outcome="SUSPICIOUS",
            ))
    return records


def _audit_spikes_with_full_scans(
    parameters: Sequence[float],
    spectra: dict[float, SpectrumRecord],
    suspicious: Sequence[SpikeAuditRecord],
    callbacks: SweepCallbacks,
    settings: SweepSettings,
    counters: SweepCounters,
) -> list[SpikeAuditRecord]:
    outcomes: list[SpikeAuditRecord] = []
    full_cache: dict[float, SpectrumRecord] = {}
    raw_full_cache: dict[float, SpectrumRecord] = {}
    original_spectra = dict(spectra)
    for item in suspicious:
        index = item.point_index
        triplet = tuple(
            float(parameters[offset]) for offset in (index - 1, index, index + 1)
        )
        try:
            for parameter in triplet:
                if parameter not in raw_full_cache:
                    raw_full_cache[parameter] = callbacks.global_search(parameter)
                    counters.global_full_scans += 1
                    counters.spike_full_scans += 1
            # Run all three scans before accepting any member.  This makes the
            # audit evidence complete even when one full record is defective.
            for parameter in triplet:
                if parameter not in full_cache:
                    full_cache[parameter] = _validate_spectrum(
                        parameter,
                        raw_full_cache[parameter],
                        callbacks,
                        settings,
                        "GLOBAL_SPIKE_AUDIT",
                    )
            full_current = full_cache[item.parameter]
            triplet_errors: dict[float, float] = {}
            for parameter in triplet:
                fast = original_spectra[parameter]
                full = full_cache[parameter]
                topology_matches = _topology_signature(fast.roots) == _topology_signature(full.roots)
                triplet_errors[parameter] = (
                    maximum_relative_error(fast.values, full.values)
                    if topology_matches
                    else math.inf
                )
            corrected = tuple(
                parameter
                for parameter in triplet
                if triplet_errors[parameter] > settings.anchor_comparison_relative
            )
            if corrected:
                for parameter in corrected:
                    spectra[parameter] = full_cache[parameter]
                outcome = "FAST_LOCATOR_CORRECTED"
            else:
                full_suspicious = detect_spikes(
                    triplet,
                    {parameter: full_cache[parameter] for parameter in triplet},
                    settings,
                )
                reproduced = any(
                    candidate.point_index == 1
                    and candidate.position == item.position
                    for candidate in full_suspicious
                )
                outcome = "REPRODUCED_BY_FULL_SCAN" if reproduced else "UNRESOLVED"
            outcomes.append(replace(
                item,
                outcome=outcome,
                full_value=float(full_current.values[item.position - 1]),
                triplet_parameters=triplet,
                corrected_parameters=corrected,
                full_triplet_max_relative=max(triplet_errors.values()),
            ))
        except Exception:
            outcomes.append(replace(
                item,
                outcome="UNRESOLVED",
                triplet_parameters=triplet,
            ))
    return outcomes


def run_spectral_sweep(
    parameters: Sequence[float],
    *,
    callbacks: SweepCallbacks | None = None,
    settings: SweepSettings = SweepSettings(),
    checkpoint: CheckpointConfig | None = None,
    resume: bool = False,
    missing_only: bool = False,
    plot_only: bool = False,
) -> SpectralSweepResult:
    """Run one sequential fast sweep with mandatory global fallback.

    ``resume`` and ``missing_only`` never recalculate a valid committed point.
    ``plot_only`` loads transactions and performs no numerical callback.
    """

    grid = np.asarray(parameters, dtype=float)
    if grid.ndim != 1 or not grid.size or np.any(~np.isfinite(grid)) or np.any(np.diff(grid) <= 0.0):
        raise SweepConfigurationError("parameter grid must be finite and increasing")
    if plot_only and checkpoint is None:
        raise SweepConfigurationError("plot_only requires checkpoint data")
    if callbacks is None and not plot_only:
        raise SweepConfigurationError("numerical sweep modes require callbacks")
    if checkpoint is not None and (resume or missing_only or plot_only):
        _require_checkpoint_identity(checkpoint)
    started = time.perf_counter()
    counters = SweepCounters(points_requested=len(grid))
    spectra: dict[float, SpectrumRecord] = {}
    audits: list[PointAudit] = []
    scheduled = set(anchor_indices(len(grid), settings.anchor_stride))
    forced_anchors: set[int] = set()
    failed_points: set[float] = set()
    qualified_points: set[float] = set()

    if checkpoint is not None and (resume or missing_only or plot_only):
        for index, raw_parameter in enumerate(grid):
            parameter = float(raw_parameter)
            record = load_point_transaction(checkpoint, index, parameter)
            if record is not None:
                try:
                    validated = _validate_roots(
                        parameter,
                        record.roots,
                        settings,
                        None if plot_only else callbacks.quality_gate,
                    )
                except SweepValidationError:
                    counters.checkpoint_rejects += 1
                else:
                    spectra[parameter] = replace(record, roots=validated)
                    counters.points_resumed += 1
                    if (
                        settings.force_anchor_after_fallback
                        and record.origin == "GLOBAL_FALLBACK"
                        and index + 1 < len(grid)
                    ):
                        forced_anchors.add(index + 1)
            elif plot_only:
                counters.checkpoint_rejects += 1
        if plot_only and len(spectra) != len(grid):
            raise SweepValidationError("plot_only checkpoint is incomplete or corrupt")
    if plot_only:
        elapsed = time.perf_counter() - started
        counters.elapsed_time = elapsed
        return SpectralSweepResult(
            spectra,
            audits,
            [],
            counters,
            "PASS",
            elapsed,
            (),
            (),
        )

    assert callbacks is not None

    def persist_failure(parameter: float) -> None:
        failed_points.add(parameter)
        if checkpoint is None:
            return
        counters.elapsed_time = time.perf_counter() - started
        remaining = len(grid) - len(spectra)
        estimate = counters.elapsed_time / max(len(spectra), 1) * remaining
        write_checkpoint_manifest(
            checkpoint,
            grid,
            spectra,
            counters,
            failed_points=sorted(failed_points),
            qualified_points=sorted(qualified_points),
            estimated_remaining_time=estimate,
        )

    for index, raw_parameter in enumerate(grid):
        parameter = float(raw_parameter)
        if parameter in spectra:
            continue
        is_anchor = index in scheduled or index in forced_anchors
        fallback_reasons: list[str] = []
        predictor_kind = "NONE"
        maximum_verification = 0.0
        maximum_predictor = 0.0
        minimum_gap = math.inf
        origin = ""
        record: SpectrumRecord | None = None

        if is_anchor or index == 0:
            counters.global_full_scans += 1
            counters.global_anchor_scans += 1
            try:
                global_record = callbacks.global_search(parameter)
                record = _validate_spectrum(
                    parameter,
                    global_record,
                    callbacks,
                    settings,
                    "GLOBAL_ANCHOR",
                )
            except Exception:
                persist_failure(parameter)
                raise
            origin = record.origin
        else:
            try:
                previous_parameter = float(grid[index - 1])
                if previous_parameter not in spectra:
                    raise SweepValidationError("previous point unavailable for continuation")
                previous = spectra[previous_parameter].values
                older_parameter = float(grid[index - 2]) if index >= 2 else None
                older = spectra[older_parameter].values if older_parameter is not None and older_parameter in spectra else None
                predictor_kind = "SECANT" if older is not None else "HOLD"
                predicted = predict_sorted_roots(
                    parameter,
                    previous_parameter,
                    previous,
                    older_parameter=older_parameter,
                    older_roots=older,
                )
                intervals = build_search_intervals(
                    previous,
                    predicted,
                    older_roots=older,
                    settings=settings,
                )
                expanded_intervals = [
                    expanded_interval(interval, settings) for interval in intervals
                ]
                if any(
                    left.upper > right.lower
                    for left, right in zip(
                        expanded_intervals,
                        expanded_intervals[1:],
                        strict=False,
                    )
                ):
                    raise SweepValidationError(
                        "expanded local verification intervals overlap"
                    )
                primary: list[RootRecord] = []
                verification: list[RootRecord] = []
                for interval, expanded in zip(
                    intervals, expanded_intervals, strict=True
                ):
                    counters.local_intervals += 1
                    counters.local_primary_searches += 1
                    if interval.is_cluster:
                        counters.cluster_intervals += 1
                    left = list(callbacks.local_search(parameter, interval, False))
                    if len(left) != interval.expected_count:
                        raise SweepValidationError("local interval returned wrong root count")
                    for root in left:
                        if not interval.lower <= root.value <= interval.upper:
                            raise SweepValidationError("local root escaped its safety interval")
                    counters.local_verification_searches += 1
                    right = list(callbacks.local_search(parameter, expanded, True))
                    if len(right) != interval.expected_count:
                        raise SweepValidationError("local verification returned wrong root count")
                    for root in right:
                        if not expanded.lower <= root.value <= expanded.upper:
                            raise SweepValidationError(
                                "local verification root escaped its safety interval"
                            )
                    for label, roots_here in (("primary", left), ("verification", right)):
                        for local_position, root in enumerate(roots_here, start=1):
                            if not root.accepted or not root.detector_agreement:
                                raise SweepValidationError(
                                    f"{label} local root {local_position} is rejected"
                                )
                            if root.possible_even_root:
                                raise SweepValidationError(
                                    f"{label} local root has unresolved even/sigma-only evidence"
                                )
                            if callbacks.quality_gate is not None:
                                passed, reason = callbacks.quality_gate(root)
                                if not passed:
                                    raise SweepValidationError(
                                        f"{label} local root failed frozen gate: {reason}"
                                    )
                    left.sort(key=lambda root: root.value)
                    right.sort(key=lambda root: root.value)
                    error = maximum_relative_error(
                        [root.value for root in left],
                        [root.value for root in right],
                    )
                    maximum_verification = max(maximum_verification, error)
                    if error > settings.local_verification_relative:
                        raise SweepValidationError("primary/local verification mismatch")
                    if any(not _root_topology_matches(a, b) for a, b in zip(left, right, strict=True)):
                        raise SweepValidationError("local cluster/nullity mismatch")
                    if _topology_signature(left) != _topology_signature(right):
                        raise SweepValidationError("local cluster partition mismatch")
                    if any(root.possible_even_root for root in (*left, *right)):
                        raise SweepValidationError("possible even-multiplicity or sigma-only root")
                    primary.extend(left)
                    verification.extend(right)
                roots = _validate_roots(parameter, primary, settings, callbacks.quality_gate)
                maximum_predictor = maximum_relative_error(
                    [root.value for root in roots], np.sort(predicted)
                )
                if maximum_predictor > settings.locator_relative_limit:
                    raise SweepValidationError("root moved beyond frozen locator threshold")
                gaps = _relative_gaps([root.value for root in roots])
                minimum_gap = float(np.min(gaps)) if gaps.size else math.inf
                if minimum_gap < settings.fallback_gap_relative:
                    raise SweepValidationError("neighboring spectral gap triggered full fallback")
                if callbacks.completeness_guard is not None:
                    passed, reason = callbacks.completeness_guard(parameter, roots)
                    if not passed:
                        raise SweepValidationError(f"completeness guard failed: {reason}")
                record = SpectrumRecord(
                    parameter=parameter,
                    roots=roots,
                    origin="FAST_LOCAL",
                    status="PASS",
                    metadata={
                        "predictor_kind": predictor_kind,
                        "predicted_values": [float(value) for value in predicted],
                        "maximum_primary_verification_relative": maximum_verification,
                    },
                )
                origin = record.origin
            except Exception as error:
                counters.predictor_failures += 1
                fallback_reasons.append(f"{type(error).__name__}: {error}")

        if record is None:
            counters.global_full_scans += 1
            counters.global_fallback_scans += 1
            counters.fallback_count += 1
            try:
                global_record = callbacks.global_search(parameter)
                record = _validate_spectrum(
                    parameter,
                    global_record,
                    callbacks,
                    settings,
                    "GLOBAL_FALLBACK",
                )
            except Exception:
                persist_failure(parameter)
                raise
            origin = record.origin
            if settings.force_anchor_after_fallback and index + 1 < len(grid):
                forced_anchors.add(index + 1)

        spectra[parameter] = record
        counters.points_committed += 1
        if checkpoint is not None:
            write_point_transaction(checkpoint, index, record)
            counters.checkpoint_transactions += 1
            counters.elapsed_time = time.perf_counter() - started
            remaining = len(grid) - len(spectra)
            estimate = counters.elapsed_time / max(len(spectra), 1) * remaining
            write_checkpoint_manifest(
                checkpoint,
                grid,
                spectra,
                counters,
                failed_points=sorted(failed_points),
                qualified_points=sorted(qualified_points),
                estimated_remaining_time=estimate,
            )
        audits.append(PointAudit(
            parameter=parameter,
            point_index=index,
            scheduled_anchor=is_anchor,
            origin=origin,
            predictor_kind=predictor_kind,
            fallback_reasons=tuple(fallback_reasons),
            maximum_primary_verification_relative=maximum_verification,
            maximum_predictor_relative=maximum_predictor,
            minimum_relative_gap=minimum_gap,
        ))

    suspicious = detect_spikes(grid, spectra, settings) if settings.enable_spike_audit else []
    counters.spike_audit_count = len(suspicious)
    spike_outcomes = _audit_spikes_with_full_scans(
        grid, spectra, suspicious, callbacks, settings, counters
    ) if suspicious else []
    qualified_points.update(
        item.parameter for item in spike_outcomes if item.outcome == "UNRESOLVED"
    )
    if checkpoint is not None:
        corrected_parameters = {
            parameter
            for outcome in spike_outcomes
            for parameter in outcome.corrected_parameters
        }
        index_by_parameter = {
            float(parameter): index for index, parameter in enumerate(grid)
        }
        for parameter in sorted(corrected_parameters):
            write_point_transaction(
                checkpoint,
                index_by_parameter[parameter],
                spectra[parameter],
            )
            counters.checkpoint_transactions += 1
        counters.elapsed_time = time.perf_counter() - started
        write_checkpoint_manifest(
            checkpoint,
            grid,
            spectra,
            counters,
            failed_points=sorted(failed_points),
            qualified_points=sorted(qualified_points),
        )
    elapsed = time.perf_counter() - started
    counters.elapsed_time = elapsed
    status = "FAIL" if any(item.outcome == "UNRESOLVED" for item in spike_outcomes) else "PASS"
    return SpectralSweepResult(
        spectra=spectra,
        point_audits=audits,
        spike_audits=spike_outcomes,
        counters=counters,
        status=status,
        runtime_seconds=elapsed,
        failed_points=tuple(sorted(failed_points)),
        qualified_points=tuple(sorted(qualified_points)),
    )


__all__ = [
    "CheckpointConfig",
    "ExactLRU",
    "GUARD_POSITION",
    "PLOT_COUNT",
    "PointAudit",
    "RUNNER_VERSION",
    "RootRecord",
    "SearchInterval",
    "SpectralSweepResult",
    "SpectrumRecord",
    "SpikeAuditRecord",
    "SweepCallbacks",
    "SweepConfigurationError",
    "SweepCounters",
    "SweepSettings",
    "SweepValidationError",
    "anchor_indices",
    "build_search_intervals",
    "detect_spikes",
    "expanded_interval",
    "load_point_transaction",
    "maximum_relative_error",
    "predict_sorted_roots",
    "run_spectral_sweep",
    "stable_fingerprint",
    "write_checkpoint_manifest",
    "write_point_transaction",
]
