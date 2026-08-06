from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np

from scripts.lib import general_spectrum_completeness as COMPLETE


DETECTOR_VERSION = "sorted_family_inventory_detector_v2"
REPAIR_ALGORITHM_VERSION = "staged_local_matrix_repair_v3_timed_cache"
CACHE_SCHEMA_VERSION = "family_inventory_local_repair_cache_v1"
SCIENTIFIC_SCOPE = "isotropic_circular_coupled_rods_eb_timoshenko"
SIGNED_SHIFT_CANDIDATES = (-2, -1, 1, 2)
MISSING_COUNTS = tuple(sorted({abs(value) for value in SIGNED_SHIFT_CANDIDATES}))


@dataclass(frozen=True)
class DetectorThresholds:
    name: str
    minimum_tail_length: int = 3
    noise_mad_multiplier: float = 3.0
    normalized_noise_floor: float = 1.0e-4
    same_rank_noise_multiplier: float = 5.0
    shifted_noise_multiplier: float = 3.0
    minimum_improvement_ratio: float = 2.5
    minimum_tail_consistency: float = 1.0
    short_buffer_minimum_jump_count: int = 4
    close_gap_fraction: float = 0.08
    window_margin_fraction: float = 0.35
    maximum_window_scale: float = 3.0

    def validate(self) -> None:
        if self.minimum_tail_length < 3:
            raise ValueError("minimum_tail_length must be at least three")
        positive = (
            self.noise_mad_multiplier,
            self.normalized_noise_floor,
            self.same_rank_noise_multiplier,
            self.shifted_noise_multiplier,
            self.minimum_improvement_ratio,
            self.short_buffer_minimum_jump_count,
            self.close_gap_fraction,
            self.window_margin_fraction,
            self.maximum_window_scale,
        )
        if any(not math.isfinite(float(item)) or float(item) <= 0.0 for item in positive):
            raise ValueError("detector thresholds must be finite and positive")
        if not 0.5 <= self.minimum_tail_consistency <= 1.0:
            raise ValueError("minimum_tail_consistency must lie in [0.5, 1]")


THRESHOLD_PROFILES: dict[str, DetectorThresholds] = {
    "permissive": DetectorThresholds(
        name="permissive",
        noise_mad_multiplier=3.0,
        same_rank_noise_multiplier=4.5,
        shifted_noise_multiplier=3.25,
        minimum_improvement_ratio=2.25,
        minimum_tail_consistency=1.0,
        short_buffer_minimum_jump_count=4,
        close_gap_fraction=0.10,
        window_margin_fraction=0.40,
    ),
    "nominal": DetectorThresholds(name="nominal"),
    "conservative": DetectorThresholds(
        name="conservative",
        noise_mad_multiplier=4.0,
        same_rank_noise_multiplier=8.0,
        shifted_noise_multiplier=2.5,
        minimum_improvement_ratio=4.0,
        minimum_tail_consistency=1.0,
        short_buffer_minimum_jump_count=4,
        close_gap_fraction=0.06,
        window_margin_fraction=0.30,
    ),
}


@dataclass(frozen=True)
class FamilySpectrum:
    family_id: str
    case_id: str
    theory: str
    epsilon_0: float
    mu: float
    eta: float
    beta_values: tuple[float, ...]
    inventories: tuple[tuple[float, ...], ...]
    point_statuses: tuple[str, ...] = ()
    required_guards: tuple[int, ...] = ()

    def validate(self) -> None:
        if not self.beta_values or len(self.beta_values) != len(self.inventories):
            raise ValueError("family beta/inventory lengths are inconsistent")
        if tuple(sorted(self.beta_values)) != self.beta_values:
            raise ValueError("family beta values must be strictly ordered")
        if len(set(self.beta_values)) != len(self.beta_values):
            raise ValueError("family beta values must be unique")
        width = len(self.inventories[0])
        if width < 3 or any(len(row) != width for row in self.inventories):
            raise ValueError("family inventories must have a common width >= 3")
        if self.point_statuses and len(self.point_statuses) != len(self.beta_values):
            raise ValueError("point status length does not match beta grid")
        if self.required_guards and len(self.required_guards) != len(self.beta_values):
            raise ValueError("required guard length does not match beta grid")

    @property
    def width(self) -> int:
        return len(self.inventories[0])


@dataclass(frozen=True)
class BoundaryEvent:
    beta_left: float
    beta_right: float
    missing_side: str
    trigger_type: str
    tail_start_rank: int
    missing_count: int
    affected_rank_count: int
    same_rank_score: float
    shifted_score: float
    improvement_ratio: float
    consistency_fraction: float
    robust_noise_scale: float


@dataclass(frozen=True)
class DetectorEvent:
    event_id: str
    family_id: str
    case_id: str
    theory: str
    beta: float
    segment_beta_left: float
    segment_beta_right: float
    trigger_types: tuple[str, ...]
    tail_start_rank: int
    best_shift: int
    affected_rank_count: int
    same_rank_score: float
    shifted_score: float
    improvement_ratio: float
    robust_noise_scale: float
    threshold_profile: str
    required_guard: int
    detector_status: str


@dataclass(frozen=True)
class RepairWindow:
    event_id: str
    case_id: str
    theory: str
    beta: float
    rank_start: int
    expected_missing_count: int
    lambda_left: float
    lambda_right: float
    source: str
    lower_anchor: float
    upper_anchor: float
    predicted_roots: tuple[float, ...]
    margin: float
    beta_probe_required: bool
    status: str


@dataclass(frozen=True)
class LocalRootEntry:
    value: float
    multiplicity: int
    repeated_root_slot: int
    block_family: str
    nullity: int
    source: str


@dataclass(frozen=True)
class LocalSearchResult:
    status: str
    repair_stage: str
    entries: tuple[LocalRootEntry, ...]
    candidate_rows: tuple[dict[str, object], ...]
    matrix_evaluations: int
    stage_matrix_evaluations: tuple[tuple[str, int], ...]
    stage_roots: tuple[tuple[str, tuple[float, ...]], ...]
    beta_probes: tuple[float, ...]
    block_classification: str
    cache_hit: bool = False
    wall_seconds: float = 0.0


def threshold_metadata() -> dict[str, dict[str, object]]:
    return {name: asdict(profile) for name, profile in THRESHOLD_PROFILES.items()}


def solver_diagnostic_trigger(
    *,
    point_status: str = "",
    unresolved_intervals: Sequence[str] = (),
    unresolved_clusters: Sequence[str] = (),
    primary_agreement: bool = True,
    root_count: int | None = None,
    required_guard: int = 11,
    guard_warning: bool = False,
    rejected_candidate_ranks: Sequence[int] = (),
) -> tuple[bool, tuple[str, ...]]:
    reasons: list[str] = []
    if unresolved_intervals:
        reasons.append("unresolved_interval")
    if unresolved_clusters:
        reasons.append("unresolved_cluster")
    if not primary_agreement:
        reasons.append("primary_disagreement")
    if root_count is not None and int(root_count) < int(required_guard):
        reasons.append("insufficient_root_count")
    if guard_warning:
        reasons.append("guard_warning")
    if any(int(rank) <= int(required_guard) for rank in rejected_candidate_ranks):
        reasons.append("rejected_candidate_below_guard")
    if "unresolved" in str(point_status).lower() and not reasons:
        reasons.append("solver_unresolved_status")
    return bool(reasons), tuple(sorted(set(reasons)))


def repair_rank_is_required(rank_start: int, required_guard: int) -> bool:
    return int(rank_start) <= int(required_guard)


def _finite_array(values: Sequence[float]) -> np.ndarray:
    return np.asarray([float(item) for item in values], dtype=float)


def normalized_mismatch(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    denominator = np.maximum(0.5 * (np.abs(left) + np.abs(right)), 1.0e-12)
    return np.abs(left - right) / denominator


def robust_noise_scale(family: FamilySpectrum, thresholds: DetectorThresholds) -> float:
    values = np.asarray(family.inventories, dtype=float)
    adjacent = normalized_mismatch(values[:-1, :], values[1:, :]).ravel()
    finite = adjacent[np.isfinite(adjacent)]
    if finite.size == 0:
        return float(thresholds.normalized_noise_floor)
    # Use the lower, internally consistent part of the rankwise residuals.
    # This stays well-defined on a three-node sparse family even when a missing
    # root shifts more than half of the stored upper ranks at both boundaries.
    cutoff = float(np.quantile(finite, 0.25))
    core = finite[finite <= cutoff]
    if core.size < 3:
        core = finite
    median = float(np.median(core))
    mad = 1.4826 * float(np.median(np.abs(core - median)))
    return max(float(thresholds.normalized_noise_floor), median + thresholds.noise_mad_multiplier * mad)


def _tail_candidate(
    observed: np.ndarray,
    reference: np.ndarray,
    *,
    missing_count: int,
    noise: float,
    thresholds: DetectorThresholds,
) -> tuple[int, int, float, float, float, float] | None:
    width = len(observed)
    accepted: list[tuple[int, int, float, float, float, float]] = []
    for start in range(width):
        indices = np.arange(start, width - missing_count, dtype=int)
        if len(indices) < thresholds.minimum_tail_length:
            continue
        same = normalized_mismatch(observed[indices], reference[indices])
        shifted = normalized_mismatch(observed[indices], reference[indices + missing_count])
        finite = np.isfinite(same) & np.isfinite(shifted)
        if int(np.count_nonzero(finite)) < thresholds.minimum_tail_length:
            continue
        same = same[finite]
        shifted = shifted[finite]
        same_score = float(np.median(same))
        shifted_score = float(np.median(shifted))
        consistency = float(np.mean(shifted < same))
        improvement = same_score / max(shifted_score, 0.05 * noise)
        if (
            same_score >= thresholds.same_rank_noise_multiplier * noise
            and shifted_score <= thresholds.shifted_noise_multiplier * noise
            and improvement >= thresholds.minimum_improvement_ratio
            and consistency + 1.0e-12 >= thresholds.minimum_tail_consistency
        ):
            accepted.append((start + 1, len(indices), same_score, shifted_score, improvement, consistency))
    return min(accepted, key=lambda item: (item[0], -item[1])) if accepted else None


def _large_jump_tail_start(
    observed: np.ndarray,
    reference: np.ndarray,
    *,
    noise: float,
    thresholds: DetectorThresholds,
) -> int | None:
    residuals = normalized_mismatch(observed, reference)
    large = residuals >= thresholds.same_rank_noise_multiplier * noise
    width = len(residuals)
    for start in range(width):
        count = width - start
        if count < thresholds.short_buffer_minimum_jump_count:
            continue
        if bool(np.all(large[start:])):
            return start
    return None


def _short_buffer_shift(
    observed: np.ndarray,
    reference: np.ndarray,
    *,
    start: int,
    noise: float,
    thresholds: DetectorThresholds,
) -> tuple[int, int, float, float, float, float] | None:
    options: list[tuple[int, int, float, float, float, float]] = []
    for missing_count in MISSING_COUNTS:
        indices = np.arange(start, len(observed) - missing_count, dtype=int)
        if len(indices) < 2:
            continue
        same = normalized_mismatch(observed[indices], reference[indices])
        shifted = normalized_mismatch(observed[indices], reference[indices + missing_count])
        same_score = float(np.median(same))
        shifted_score = float(np.median(shifted))
        improvement = same_score / max(shifted_score, 0.05 * noise)
        consistency = float(np.mean(shifted < same))
        if (
            shifted_score <= thresholds.shifted_noise_multiplier * noise
            and improvement >= thresholds.minimum_improvement_ratio
            and consistency >= thresholds.minimum_tail_consistency
        ):
            options.append((missing_count, len(indices), same_score, shifted_score, improvement, consistency))
    return min(options, key=lambda item: (item[3], -item[0])) if options else None


def boundary_events(
    family: FamilySpectrum,
    thresholds: DetectorThresholds,
) -> tuple[BoundaryEvent, ...]:
    family.validate()
    thresholds.validate()
    noise = robust_noise_scale(family, thresholds)
    rows = [np.asarray(item, dtype=float) for item in family.inventories]
    output: list[BoundaryEvent] = []
    for index in range(len(rows) - 1):
        left = rows[index]
        right = rows[index + 1]
        beta_left = float(family.beta_values[index])
        beta_right = float(family.beta_values[index + 1])
        for missing_side, observed, reference in (
            ("right", right, left),
            ("left", left, right),
        ):
            candidates: list[BoundaryEvent] = []
            signed_candidates = (
                (shift for shift in SIGNED_SHIFT_CANDIDATES if shift > 0)
                if missing_side == "right"
                else (shift for shift in SIGNED_SHIFT_CANDIDATES if shift < 0)
            )
            for signed_shift in signed_candidates:
                missing_count = abs(signed_shift)
                item = _tail_candidate(
                    observed,
                    reference,
                    missing_count=missing_count,
                    noise=noise,
                    thresholds=thresholds,
                )
                if item is None:
                    continue
                rank, affected, same_score, shifted_score, improvement, consistency = item
                candidates.append(
                    BoundaryEvent(
                        beta_left,
                        beta_right,
                        missing_side,
                        "tail_shift",
                        rank,
                        missing_count,
                        affected,
                        same_score,
                        shifted_score,
                        improvement,
                        consistency,
                        noise,
                    )
                )
            if candidates:
                # The earliest rank whose complete upper tail is explained by
                # a +1/+2 displacement is the first affected sorted slot.
                # Choosing only the smallest residual can incorrectly move
                # the boundary one slot upward when both descriptions fit.
                output.append(min(candidates, key=lambda event: (event.tail_start_rank, event.shifted_score)))
                continue
            start = _large_jump_tail_start(observed, reference, noise=noise, thresholds=thresholds)
            if start is None:
                continue
            short = _short_buffer_shift(
                observed,
                reference,
                start=start,
                noise=noise,
                thresholds=thresholds,
            )
            if short is None:
                continue
            missing_count, affected, same_score, shifted_score, improvement, consistency = short
            output.append(
                BoundaryEvent(
                    beta_left,
                    beta_right,
                    missing_side,
                    "gap_anomaly_short_buffer",
                    start + 1,
                    missing_count,
                    max(affected, len(observed) - start),
                    same_score,
                    shifted_score,
                    improvement,
                    consistency,
                    noise,
                )
            )
    output.sort(
        key=lambda event: (
            event.beta_left,
            event.beta_right,
            event.missing_side,
            event.tail_start_rank,
            event.missing_count,
        )
    )
    return tuple(output)


def _two_sided_support(family: FamilySpectrum, beta: float, rank: int, missing_count: int) -> bool:
    beta_to_index = {float(value): index for index, value in enumerate(family.beta_values)}
    index = beta_to_index[float(beta)]
    if index <= 0 or index >= len(family.beta_values) - 1:
        return False
    current = np.asarray(family.inventories[index], dtype=float)
    predicted = 0.5 * (
        np.asarray(family.inventories[index - 1], dtype=float)
        + np.asarray(family.inventories[index + 1], dtype=float)
    )
    start = max(0, int(rank) - 1)
    indices = np.arange(start, len(current) - int(missing_count), dtype=int)
    if len(indices) < 2:
        return False
    same = normalized_mismatch(current[indices], predicted[indices])
    shifted = normalized_mismatch(current[indices], predicted[indices + int(missing_count)])
    return float(np.mean(shifted < same)) >= 0.75


def _cluster_boundary_events(
    family: FamilySpectrum,
    boundaries: Sequence[BoundaryEvent],
    thresholds: DetectorThresholds,
) -> tuple[DetectorEvent, ...]:
    if not boundaries:
        return ()
    beta_min = float(family.beta_values[0])
    beta_max = float(family.beta_values[-1])
    grouped: dict[tuple[int, int], list[BoundaryEvent]] = {}
    for boundary in boundaries:
        grouped.setdefault((boundary.tail_start_rank, boundary.missing_count), []).append(boundary)
    segments: list[tuple[float, float, int, int, tuple[BoundaryEvent, ...]]] = []
    for (rank, missing_count), items in sorted(grouped.items()):
        items = sorted(items, key=lambda item: (item.beta_left, item.beta_right, item.missing_side))
        open_start: float | None = None
        support: list[BoundaryEvent] = []
        for item in items:
            if item.missing_side == "right":
                if open_start is not None:
                    segments.append((open_start, float(item.beta_left), rank, missing_count, tuple(support)))
                open_start = float(item.beta_right)
                support = [item]
            else:
                if open_start is None:
                    open_start = beta_min
                    support = []
                support.append(item)
                segments.append((open_start, float(item.beta_left), rank, missing_count, tuple(support)))
                open_start = None
                support = []
        if open_start is not None:
            segments.append((open_start, beta_max, rank, missing_count, tuple(support)))
    beta_array = np.asarray(family.beta_values, dtype=float)
    events: list[DetectorEvent] = []
    counter = 0
    for beta_left, beta_right, rank, missing_count, support in sorted(segments):
        selected = beta_array[(beta_array >= beta_left - 1.0e-12) & (beta_array <= beta_right + 1.0e-12)]
        if selected.size == 0:
            continue
        same_score = max((item.same_rank_score for item in support), default=float("nan"))
        shifted_score = min((item.shifted_score for item in support), default=float("nan"))
        improvement = max((item.improvement_ratio for item in support), default=float("nan"))
        affected = max((item.affected_rank_count for item in support), default=missing_count)
        noise = max((item.robust_noise_scale for item in support), default=robust_noise_scale(family, thresholds))
        trigger_types = sorted({item.trigger_type for item in support})
        for beta in selected:
            counter += 1
            point_index = int(np.where(beta_array == beta)[0][0])
            guard = int(family.required_guards[point_index]) if family.required_guards else family.width
            point_triggers = list(trigger_types)
            if _two_sided_support(family, float(beta), rank, missing_count):
                point_triggers.append("two_sided_prediction")
            local_values = np.asarray(family.inventories[point_index], dtype=float)
            gaps = np.diff(local_values)
            finite_gaps = gaps[np.isfinite(gaps) & (gaps > 0.0)]
            if finite_gaps.size:
                median_gap = float(np.median(finite_gaps))
                start = max(0, rank - 2)
                stop = min(len(gaps), rank + missing_count)
                if any(float(gap) <= thresholds.close_gap_fraction * median_gap for gap in gaps[start:stop]):
                    point_triggers.append("close_pair_neighborhood")
            status = "repair_trigger" if rank <= guard else "upper_spectrum_audit_trigger"
            event_id = f"{family.case_id}_{family.theory}_{float(beta):09.4f}_{counter:03d}"
            events.append(
                DetectorEvent(
                    event_id=event_id,
                    family_id=family.family_id,
                    case_id=family.case_id,
                    theory=family.theory,
                    beta=float(beta),
                    segment_beta_left=float(beta_left),
                    segment_beta_right=float(beta_right),
                    trigger_types=tuple(sorted(set(point_triggers))),
                    tail_start_rank=int(rank),
                    best_shift=int(missing_count),
                    affected_rank_count=int(affected),
                    same_rank_score=float(same_score),
                    shifted_score=float(shifted_score),
                    improvement_ratio=float(improvement),
                    robust_noise_scale=float(noise),
                    threshold_profile=thresholds.name,
                    required_guard=guard,
                    detector_status=status,
                )
            )
    events.sort(key=lambda event: (event.case_id, event.theory, event.beta, event.tail_start_rank))
    return tuple(events)


def detect_family_inventory(
    family: FamilySpectrum,
    thresholds: DetectorThresholds,
) -> tuple[DetectorEvent, ...]:
    return _cluster_boundary_events(family, boundary_events(family, thresholds), thresholds)


def infer_repair_window(
    family: FamilySpectrum,
    event: DetectorEvent,
    thresholds: DetectorThresholds,
) -> RepairWindow:
    beta_array = np.asarray(family.beta_values, dtype=float)
    matches = np.where(np.isclose(beta_array, event.beta, atol=1.0e-12, rtol=0.0))[0]
    if matches.size != 1:
        raise ValueError("event beta is not a unique family node")
    point_index = int(matches[0])
    current = np.asarray(family.inventories[point_index], dtype=float)
    rank_index = int(event.tail_start_rank) - 1
    missing_count = int(event.best_shift)
    clean_indices = [
        index
        for index, beta in enumerate(beta_array)
        if beta < event.segment_beta_left - 1.0e-12 or beta > event.segment_beta_right + 1.0e-12
    ]
    left_candidates = [index for index in clean_indices if index < point_index]
    right_candidates = [index for index in clean_indices if index > point_index]
    neighbor_indices: list[int] = []
    if left_candidates:
        neighbor_indices.append(max(left_candidates))
    if right_candidates:
        neighbor_indices.append(min(right_candidates))
    if not neighbor_indices:
        return RepairWindow(
            event.event_id,
            event.case_id,
            event.theory,
            event.beta,
            event.tail_start_rank,
            missing_count,
            float("nan"),
            float("nan"),
            "no_clean_neighbor",
            float("nan"),
            float("nan"),
            (),
            float("nan"),
            False,
            "defer_window_ambiguous",
        )
    neighbor_rows = np.asarray([family.inventories[index] for index in neighbor_indices], dtype=float)
    # A displaced tail can occasionally be recognized one slot above the
    # actual inventory hole: the first shifted rank is then itself a member
    # of a very close pair.  Resolve that purely from sorted values by finding
    # the nearest lower insertion slot whose neighboring prediction is fully
    # bracketed by consecutive roots of the defective inventory.
    for candidate_rank_index in range(rank_index, max(-1, rank_index - missing_count), -1):
        candidate_prediction = np.median(
            neighbor_rows[:, candidate_rank_index : candidate_rank_index + missing_count], axis=0
        )
        if len(candidate_prediction) != missing_count or not np.all(np.isfinite(candidate_prediction)):
            continue
        lower = current[candidate_rank_index - 1] if candidate_rank_index > 0 else -np.inf
        upper = current[candidate_rank_index] if candidate_rank_index < len(current) else np.inf
        tolerance = COMPLETE.DEFAULT_ROOT_MATCH_TOL
        if float(np.min(candidate_prediction)) >= float(lower) - tolerance and float(
            np.max(candidate_prediction)
        ) <= float(upper) + tolerance:
            rank_index = int(candidate_rank_index)
            break
    predicted = np.median(neighbor_rows[:, rank_index : rank_index + missing_count], axis=0)
    finite_predicted = predicted[np.isfinite(predicted)]
    if len(finite_predicted) != missing_count:
        status = "defer_window_ambiguous"
        predicted_values: tuple[float, ...] = tuple(float(item) for item in finite_predicted)
        return RepairWindow(
            event.event_id,
            event.case_id,
            event.theory,
            event.beta,
            event.tail_start_rank,
            missing_count,
            float("nan"),
            float("nan"),
            "incomplete_neighbor_prediction",
            float("nan"),
            float("nan"),
            predicted_values,
            float("nan"),
            False,
            status,
        )
    local_gaps: list[float] = []
    for row in neighbor_rows:
        gaps = np.diff(row)
        start = max(0, rank_index - 2)
        stop = min(len(gaps), rank_index + missing_count + 1)
        local_gaps.extend(float(item) for item in gaps[start:stop] if math.isfinite(float(item)) and float(item) > 0.0)
    local_scale = float(np.median(local_gaps)) if local_gaps else max(1.0e-3, abs(float(np.median(finite_predicted))) * 0.05)
    margin = thresholds.window_margin_fraction * local_scale
    if rank_index > 0 and math.isfinite(float(current[rank_index - 1])):
        lower_anchor = float(current[rank_index - 1])
    else:
        lower_anchor = float(np.min(finite_predicted) - local_scale)
    if rank_index < len(current) and math.isfinite(float(current[rank_index])):
        upper_anchor = float(current[rank_index])
    else:
        upper_anchor = float(np.max(finite_predicted) + local_scale)
    left = max(lower_anchor + 1.0e-8, float(np.min(finite_predicted) - margin))
    right = min(upper_anchor - 1.0e-8, float(np.max(finite_predicted) + margin))
    if right <= left:
        left = float(np.min(finite_predicted) - margin)
        right = float(np.max(finite_predicted) + margin)
    width = right - left
    ambiguous = (
        not math.isfinite(width)
        or width <= 20.0 * COMPLETE.DEFAULT_BRENT_XTOL
        or width > thresholds.maximum_window_scale * local_scale
        or np.min(finite_predicted) < left - COMPLETE.DEFAULT_ROOT_MATCH_TOL
        or np.max(finite_predicted) > right + COMPLETE.DEFAULT_ROOT_MATCH_TOL
    )
    return RepairWindow(
        event.event_id,
        event.case_id,
        event.theory,
        event.beta,
        rank_index + 1,
        missing_count,
        float(left),
        float(right),
        "neighbor_sorted_anchors",
        lower_anchor,
        upper_anchor,
        tuple(float(item) for item in finite_predicted),
        float(margin),
        False,
        "defer_window_ambiguous" if ambiguous else "window_inferred",
    )


def compute_n_true(
    eb_roots: Sequence[float],
    timo_roots: Sequence[float],
    *,
    k_max: int = 10,
    delta_tol: float = 0.10,
) -> tuple[int | None, int | None, int, tuple[float, ...]]:
    eb = [float(item) for item in eb_roots]
    timo = [float(item) for item in timo_roots]
    if min(len(eb), len(timo)) < k_max:
        return None, None, k_max + 1, ()
    deltas: list[float] = []
    first_failed: int | None = None
    for index in range(k_max):
        if not (math.isfinite(eb[index]) and math.isfinite(timo[index]) and timo[index] != 0.0):
            return None, None, k_max + 1, tuple(deltas)
        delta = abs(eb[index] ** 2 - timo[index] ** 2) / (timo[index] ** 2)
        deltas.append(float(delta))
        if first_failed is None and delta > delta_tol:
            first_failed = index + 1
    n_true = k_max if first_failed is None else first_failed - 1
    guard = 11 if n_true == k_max else first_failed + 1
    if min(len(eb), len(timo)) < guard:
        return None, first_failed, guard, tuple(deltas)
    return int(n_true), first_failed, int(guard), tuple(deltas)


def _candidate_nullity(candidate: COMPLETE.RootCandidate, settings: COMPLETE.SearchSettings) -> int:
    diag = candidate.diagnostics
    return 2 if (
        diag.sigma_2 <= settings.nullity_sigma
        and math.isfinite(diag.sigma_3)
        and diag.sigma_3 > 0.0
        and diag.sigma_2 / diag.sigma_3 <= settings.sigma_ratio_accept
    ) else 1


def _phase_support(candidate: COMPLETE.RootCandidate) -> tuple[bool, bool]:
    sources = "+".join(candidate.detection_sources)
    return "unshifted" in sources, "shifted" in sources


def _dense_search(
    provider: Callable[[float], np.ndarray],
    *,
    lambda_left: float,
    lambda_right: float,
    scan_step: float,
    stage: str,
    block_family: str,
    base_settings: COMPLETE.SearchSettings,
    phases: Sequence[float] = (0.0, 0.5),
) -> tuple[tuple[LocalRootEntry, ...], tuple[dict[str, object], ...], int]:
    settings = replace(
        base_settings,
        lambda_min=float(lambda_left),
        lambda_max=float(lambda_right),
        scan_step=float(scan_step),
        max_upper_growth_tries=1,
    )
    operations = COMPLETE.OperationCounts()
    evaluator = COMPLETE._MatrixEvaluator(provider, operations)
    candidates, _intervals, _unresolved = COMPLETE._global_candidates(
        evaluator,
        settings,
        configuration="primary",
        scan_step=float(scan_step),
        phases=tuple(float(item) for item in phases),
        upper=float(lambda_right),
        seed_roots=(),
        seed_source="automatic_local_repair_no_seed",
    )
    merged = COMPLETE._merge_candidates(candidates, settings)
    records = COMPLETE._root_records(merged, settings)
    entries: list[LocalRootEntry] = []
    slots: dict[tuple[float, int], int] = {}
    for record in records:
        multiplicity = max(1, int(record.detected_nullity))
        key = (round(float(record.Lambda), 10), multiplicity)
        slots[key] = slots.get(key, 0) + 1
        entries.append(
            LocalRootEntry(
                value=float(record.Lambda),
                multiplicity=multiplicity,
                repeated_root_slot=slots[key] if multiplicity > 1 else 1,
                block_family=block_family,
                nullity=multiplicity,
                source=stage + ":" + "+".join(record.detection_sources),
            )
        )
    rows: list[dict[str, object]] = []
    for candidate in candidates:
        unshifted, shifted = _phase_support(candidate)
        rows.append(
            {
                "stage": stage,
                "lambda_candidate": float(candidate.Lambda),
                "source": "+".join(candidate.detection_sources),
                "sign_change": any("sign_change" in item or "grid_zero" in item for item in candidate.detection_sources),
                "sigma_min": float(candidate.diagnostics.sigma_1),
                "sigma_ratio": float(candidate.diagnostics.sigma_ratio),
                "residual": float(candidate.diagnostics.sigma_1),
                "accepted": candidate.acceptance_status == "accepted_full_matrix_svd",
                "rejection_reason": "" if candidate.acceptance_status == "accepted_full_matrix_svd" else candidate.acceptance_status,
                "multiplicity": _candidate_nullity(candidate, settings),
                "block_family": block_family,
                "unshifted_phase_support": unshifted,
                "shifted_phase_support": shifted,
                "bracket_left": float(candidate.interval_left),
                "bracket_right": float(candidate.interval_right),
            }
        )
    return tuple(entries), tuple(rows), int(operations.characteristic_matrix_evaluations)


def root_sets_stable(
    left: Sequence[LocalRootEntry],
    right: Sequence[LocalRootEntry],
    *,
    tolerance: float,
) -> bool:
    if len(left) != len(right):
        return False
    return all(
        abs(float(a.value) - float(b.value)) <= tolerance
        and int(a.multiplicity) == int(b.multiplicity)
        for a, b in zip(left, right)
    )


def _annotate_block_provenance(
    full_entries: Sequence[LocalRootEntry],
    block_entries: Mapping[str, Sequence[LocalRootEntry]],
    *,
    root_match_tolerance: float,
) -> tuple[tuple[LocalRootEntry, ...], str]:
    annotated: list[LocalRootEntry] = []
    for entry in full_entries:
        families = [
            name
            for name, roots in block_entries.items()
            if any(abs(float(entry.value) - float(root.value)) <= root_match_tolerance for root in roots)
        ]
        annotated.append(replace(entry, block_family="+".join(sorted(families)) if families else "full_matrix_only"))
    family_entries = [entry for entry in annotated if entry.block_family != "full_matrix_only"]
    double = any(entry.multiplicity > 1 and "+" in entry.block_family for entry in family_entries)
    distinct = any(
        right.value - left.value > root_match_tolerance
        for left, right in zip(family_entries, family_entries[1:])
    )
    classification = "resolved_double_root" if double else ("resolved_distinct_pair" if distinct else "block_pair_not_present")
    return tuple(annotated), classification


def staged_local_search(
    provider: Callable[[float], np.ndarray],
    window: RepairWindow,
    *,
    base_settings: COMPLETE.SearchSettings,
    block_providers: Mapping[str, Callable[[float], np.ndarray]] | None = None,
) -> LocalSearchResult:
    if window.status != "window_inferred":
        return LocalSearchResult(
            status="deferred_window_ambiguous",
            repair_stage="",
            entries=(),
            candidate_rows=(),
            matrix_evaluations=0,
            stage_matrix_evaluations=(),
            stage_roots=(),
            beta_probes=(),
            block_classification="not_attempted",
        )
    width = float(window.lambda_right - window.lambda_left)
    stage_steps = (
        ("L1", min(1.0e-3, max(1.0e-4, width / 400.0))),
        ("L2", min(2.5e-4, max(1.0e-4, width / 1200.0))),
        ("L3", 1.0e-4),
    )
    all_rows: list[dict[str, object]] = []
    stage_entries: list[tuple[str, tuple[LocalRootEntry, ...]]] = []
    stage_counts: list[tuple[str, int]] = []
    total_evaluations = 0
    final_entries: tuple[LocalRootEntry, ...] = ()
    final_stage = ""
    for stage, step in stage_steps:
        entries, rows, evaluations = _dense_search(
            provider,
            lambda_left=window.lambda_left,
            lambda_right=window.lambda_right,
            scan_step=step,
            stage=stage,
            block_family="full_matrix",
            base_settings=base_settings,
        )
        all_rows.extend(rows)
        stage_entries.append((stage, entries))
        stage_counts.append((stage, evaluations))
        total_evaluations += evaluations
        if len(stage_entries) >= 2:
            previous = stage_entries[-2][1]
            stable = root_sets_stable(previous, entries, tolerance=base_settings.root_match_tol)
            enough = len(entries) >= int(window.expected_missing_count)
            if stable and enough:
                final_entries = entries
                final_stage = stage
                break
    if not final_entries:
        return LocalSearchResult(
            status="deferred_local_inventory_unstable",
            repair_stage=stage_entries[-1][0] if stage_entries else "",
            entries=stage_entries[-1][1] if stage_entries else (),
            candidate_rows=tuple(all_rows),
            matrix_evaluations=total_evaluations,
            stage_matrix_evaluations=tuple(stage_counts),
            stage_roots=tuple((stage, tuple(entry.value for entry in entries)) for stage, entries in stage_entries),
            beta_probes=(),
            block_classification="not_resolved",
        )
    block_classification = "not_applicable"
    if block_providers:
        block_results: dict[str, tuple[LocalRootEntry, ...]] = {}
        final_step = next(step for stage, step in stage_steps if stage == final_stage)
        for block_family, block_provider in sorted(block_providers.items()):
            block_entries, block_rows, evaluations = _dense_search(
                block_provider,
                lambda_left=window.lambda_left,
                lambda_right=window.lambda_right,
                scan_step=final_step,
                stage=final_stage + "_block",
                block_family=block_family,
                base_settings=base_settings,
            )
            block_results[block_family] = block_entries
            all_rows.extend(block_rows)
            total_evaluations += evaluations
            stage_counts.append((final_stage + "_" + block_family, evaluations))
        final_entries, block_classification = _annotate_block_provenance(
            final_entries,
            block_results,
            root_match_tolerance=base_settings.root_match_tol,
        )
    return LocalSearchResult(
        status="resolved_after_local_repair",
        repair_stage=final_stage,
        entries=tuple(final_entries),
        candidate_rows=tuple(all_rows),
        matrix_evaluations=total_evaluations,
        stage_matrix_evaluations=tuple(stage_counts),
        stage_roots=tuple((stage, tuple(entry.value for entry in entries)) for stage, entries in stage_entries),
        beta_probes=(),
        block_classification=block_classification,
    )


def merge_inventory(
    original: Sequence[float],
    recovered: Sequence[LocalRootEntry],
    window: RepairWindow,
    *,
    root_dedup_tolerance: float,
    limit: int = 12,
) -> tuple[float, ...]:
    kept = [
        float(value)
        for value in original
        if math.isfinite(float(value))
        and not (
            window.lambda_left - root_dedup_tolerance
            <= float(value)
            <= window.lambda_right + root_dedup_tolerance
        )
    ]
    combined: list[tuple[float, int, int]] = [(value, 1, 1) for value in kept]
    combined.extend((entry.value, entry.multiplicity, entry.repeated_root_slot) for entry in recovered)
    combined.sort(key=lambda item: (item[0], item[2]))
    output: list[tuple[float, int, int]] = []
    for item in combined:
        if not output:
            output.append(item)
            continue
        previous = output[-1]
        preserve_slots = item[1] > 1 and previous[1] > 1 and item[2] != previous[2]
        if abs(item[0] - previous[0]) <= root_dedup_tolerance and not preserve_slots:
            continue
        output.append(item)
    return tuple(float(item[0]) for item in output[:limit])


def cache_identity(
    *,
    family_id: str,
    theory: str,
    beta: float,
    source_hash: str,
    threshold_profile: DetectorThresholds,
    window: RepairWindow,
    base_settings: COMPLETE.SearchSettings,
    beta_probes: Sequence[float] = (),
    scientific_scope: str = SCIENTIFIC_SCOPE,
) -> dict[str, object]:
    validate_scientific_scope(scientific_scope)
    tolerance_payload = {
        key: value
        for key, value in asdict(base_settings).items()
        if "tol" in key or "sigma" in key or "gap" in key
    }
    tolerance_hash = hashlib.sha256(
        json.dumps(tolerance_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    identity = {
        "cache_schema_version": CACHE_SCHEMA_VERSION,
        "scientific_scope": scientific_scope,
        "detector_version": DETECTOR_VERSION,
        "repair_algorithm_version": REPAIR_ALGORITHM_VERSION,
        "family_id": family_id,
        "theory": theory,
        "beta": float(beta),
        "source_hash": source_hash,
        "threshold_profile": asdict(threshold_profile),
        "inferred_window": asdict(window),
        "lambda_steps": [1.0e-3, 2.5e-4, 1.0e-4],
        "beta_probes": [float(item) for item in beta_probes],
        "general_spectrum_evaluator_version": COMPLETE.GENERAL_SPECTRUM_ALGORITHM_VERSION,
        "eb_evaluator_version": COMPLETE.EB_MATRIX_EVALUATOR_VERSION,
        "timo_evaluator_version": COMPLETE.TIMO.TIMOSHENKO_BASIS_EVALUATOR_VERSION,
        "tolerance_hash": tolerance_hash,
        "force_strict_enabled": False,
    }
    # Cache files pass through JSON, where tuples become lists.  Returning the
    # same canonical JSON-native structure prevents a false identity miss on
    # every resume.
    return json.loads(json.dumps(identity, sort_keys=True, separators=(",", ":")))


def validate_scientific_scope(scientific_scope: str) -> str:
    """Accept only the isotropic circular EB/Timoshenko research scope.

    The family repair helper is intentionally not a generic constitutive-model
    dispatcher.  Keeping this boundary explicit prevents a diagnostic cache or
    orchestration preset from being reused for a different section/material
    model merely because the numerical arrays happen to have a similar shape.
    """

    if scientific_scope != SCIENTIFIC_SCOPE:
        raise ValueError(
            "family inventory repair supports only "
            f"{SCIENTIFIC_SCOPE!r}; received {scientific_scope!r}"
        )
    return scientific_scope


def load_cache(path: Path, identity: Mapping[str, object]) -> LocalSearchResult | None:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return None
    if payload.get("identity") != dict(identity):
        return None
    raw = payload.get("result")
    if not isinstance(raw, dict):
        return None
    return LocalSearchResult(
        status=str(raw["status"]),
        repair_stage=str(raw["repair_stage"]),
        entries=tuple(LocalRootEntry(**item) for item in raw["entries"]),
        candidate_rows=tuple(dict(item) for item in raw["candidate_rows"]),
        matrix_evaluations=int(raw["matrix_evaluations"]),
        stage_matrix_evaluations=tuple((str(item[0]), int(item[1])) for item in raw["stage_matrix_evaluations"]),
        stage_roots=tuple((str(item[0]), tuple(float(value) for value in item[1])) for item in raw["stage_roots"]),
        beta_probes=tuple(float(item) for item in raw["beta_probes"]),
        block_classification=str(raw["block_classification"]),
        cache_hit=True,
        wall_seconds=float(raw.get("wall_seconds", 0.0)),
    )


def save_cache(path: Path, identity: Mapping[str, object], result: LocalSearchResult) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    payload = {
        "identity": dict(identity),
        "result": {
            **asdict(result),
            "cache_hit": False,
        },
    }
    temporary.write_text(json.dumps(payload, sort_keys=True, indent=2), encoding="utf-8")
    os.replace(temporary, path)


def result_without_cache_flag(result: LocalSearchResult) -> LocalSearchResult:
    return replace(result, cache_hit=False)


def direct_window_from_inventory(
    *,
    event_id: str,
    case_id: str,
    theory: str,
    beta: float,
    roots: Sequence[float],
    rank_start: int,
    expected_root_count: int,
    thresholds: DetectorThresholds,
) -> RepairWindow:
    values = [float(item) for item in roots if math.isfinite(float(item))]
    index = max(0, min(len(values) - 1, int(rank_start) - 1))
    count = max(1, min(2, int(expected_root_count)))
    selected = values[index : min(len(values), index + count)]
    if not selected:
        return RepairWindow(
            event_id, case_id, theory, beta, rank_start, count,
            float("nan"), float("nan"), "direct_diagnostic_no_anchor",
            float("nan"), float("nan"), (), float("nan"), False, "defer_window_ambiguous",
        )
    gaps = [values[i + 1] - values[i] for i in range(max(0, index - 2), min(len(values) - 1, index + count + 1))]
    finite_gaps = [gap for gap in gaps if math.isfinite(gap) and gap > 0.0]
    local_scale = float(np.median(finite_gaps)) if finite_gaps else max(1.0e-3, abs(selected[0]) * 0.05)
    margin = thresholds.window_margin_fraction * local_scale
    lower_anchor = values[index - 1] if index > 0 else selected[0] - local_scale
    upper_index = min(len(values) - 1, index + count)
    upper_anchor = values[upper_index] if upper_index > index + count - 1 else selected[-1] + local_scale
    left = max(lower_anchor + 1.0e-8, min(selected) - margin)
    right = min(upper_anchor - 1.0e-8, max(selected) + margin)
    if right <= left:
        left = min(selected) - margin
        right = max(selected) + margin
    status = "window_inferred" if right > left else "defer_window_ambiguous"
    return RepairWindow(
        event_id=event_id,
        case_id=case_id,
        theory=theory,
        beta=float(beta),
        rank_start=int(rank_start),
        expected_missing_count=count,
        lambda_left=float(left),
        lambda_right=float(right),
        source="solver_diagnostic_local_anchors",
        lower_anchor=float(lower_anchor),
        upper_anchor=float(upper_anchor),
        predicted_roots=tuple(selected),
        margin=float(margin),
        beta_probe_required=False,
        status=status,
    )


__all__ = [
    "CACHE_SCHEMA_VERSION",
    "DETECTOR_VERSION",
    "REPAIR_ALGORITHM_VERSION",
    "DetectorEvent",
    "DetectorThresholds",
    "FamilySpectrum",
    "LocalRootEntry",
    "LocalSearchResult",
    "RepairWindow",
    "THRESHOLD_PROFILES",
    "boundary_events",
    "cache_identity",
    "SCIENTIFIC_SCOPE",
    "validate_scientific_scope",
    "compute_n_true",
    "detect_family_inventory",
    "direct_window_from_inventory",
    "infer_repair_window",
    "load_cache",
    "merge_inventory",
    "normalized_mismatch",
    "repair_rank_is_required",
    "result_without_cache_flag",
    "root_sets_stable",
    "save_cache",
    "solver_diagnostic_trigger",
    "staged_local_search",
    "threshold_metadata",
]
