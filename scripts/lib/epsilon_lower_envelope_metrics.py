from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Mapping, Sequence


EXPECTED_PREFIX_GROUPS: dict[str, tuple[int, ...]] = {
    "prefixes_2": (2,),
    "prefixes_3": (3,),
    "prefixes_4_5": (4, 5),
    "prefixes_6": (6,),
    "prefixes_7": (7,),
    "prefixes_8": (8,),
    "prefixes_9_10": (9, 10),
}


@dataclass(frozen=True)
class PrefixMetric:
    prefix_n: int
    Delta_n: float
    V_n: float
    M_n: float
    triggering_indices: tuple[int, ...]


def squared_frequency_delta(lambda_eb: float, lambda_timo: float) -> float:
    eb = float(lambda_eb)
    timo = float(lambda_timo)
    if not math.isfinite(eb) or not math.isfinite(timo) or abs(timo) <= 1.0e-30:
        return float("nan")
    return abs(eb * eb - timo * timo) / (timo * timo)


def running_prefix_metrics(
    deltas: Sequence[float],
    *,
    threshold: float = 0.10,
    tie_tolerance: float = 1.0e-12,
) -> tuple[PrefixMetric, ...]:
    if not 0.0 < float(threshold) < 1.0:
        raise ValueError("threshold must lie in (0, 1)")
    if tie_tolerance < 0.0 or not math.isfinite(float(tie_tolerance)):
        raise ValueError("tie_tolerance must be finite and nonnegative")
    metrics: list[PrefixMetric] = []
    prefix: list[float] = []
    for index, raw in enumerate(deltas, start=1):
        value = float(raw)
        prefix.append(value)
        if not all(math.isfinite(item) for item in prefix):
            maximum = float("nan")
            triggers: tuple[int, ...] = ()
        else:
            maximum = max(prefix)
            triggers = tuple(
                position
                for position, item in enumerate(prefix, start=1)
                if abs(item - maximum) <= tie_tolerance
            )
        violation = maximum - threshold if math.isfinite(maximum) else float("nan")
        metrics.append(
            PrefixMetric(
                prefix_n=index,
                Delta_n=maximum,
                V_n=violation,
                M_n=-violation,
                triggering_indices=triggers,
            )
        )
    return tuple(metrics)


def true_safe_prefix(deltas: Sequence[float], *, threshold: float = 0.10) -> int:
    """Return the continuous safe prefix; a late pass cannot repair a failure."""

    count = 0
    for raw in deltas:
        value = float(raw)
        if not math.isfinite(value) or value > threshold:
            break
        count += 1
    return count


def late_pass_indices(deltas: Sequence[float], *, threshold: float = 0.10) -> tuple[int, ...]:
    first_failure_seen = False
    late: list[int] = []
    for index, raw in enumerate(deltas, start=1):
        value = float(raw)
        if not math.isfinite(value) or value > threshold:
            first_failure_seen = True
        elif first_failure_seen:
            late.append(index)
    return tuple(late)


def certified_prefix(
    epsilon: float,
    safe_lower_by_prefix: Mapping[int, float],
    *,
    k_max: int = 10,
) -> int:
    epsilon_f = float(epsilon)
    certified = 0
    for prefix_n in range(1, int(k_max) + 1):
        safe = float(safe_lower_by_prefix[prefix_n])
        if math.isfinite(safe) and epsilon_f <= safe:
            certified = prefix_n
    return certified


def primary_status(
    *,
    V_n: float,
    N_true: int,
    required_prefix_n: int,
    N_certified_0: int,
    strict_margin: float,
    k10_resolved: bool,
    baseline_control_failure: bool = False,
) -> str:
    if baseline_control_failure:
        return "baseline_control_failure"
    if not k10_resolved or not math.isfinite(float(V_n)):
        return "unresolved_spectrum"
    certificate_probe = int(required_prefix_n) == int(N_certified_0)
    if (
        float(V_n) > 0.0
        or int(N_true) < int(required_prefix_n)
        or (certificate_probe and int(N_true) < int(N_certified_0))
    ):
        return "provisional_counterexample"
    if -float(strict_margin) <= float(V_n) <= 0.0:
        return "primary_near_boundary"
    return "primary_safe"


def needs_strict_verification(primary: str, quality_triggers: Sequence[str] = ()) -> bool:
    return primary in {"provisional_counterexample", "primary_near_boundary"} or bool(quality_triggers)


def final_status(
    *,
    primary: str,
    V_primary: float,
    V_verification: float | None,
    N_true_primary: int,
    N_true_verification: int | None,
    required_prefix_n: int,
    N_certified_0: int,
    primary_k10: bool,
    verification_k10: bool | None,
    roots_agree: bool | None,
    clusters_agree: bool | None,
    violation_tolerance: float,
    strict_required: bool,
) -> str:
    if primary == "baseline_control_failure":
        return "unresolved_spectrum"
    if not primary_k10:
        return "unresolved_spectrum"
    if not strict_required:
        return "screen_safe_not_strictly_recomputed"
    if (
        V_verification is None
        or N_true_verification is None
        or verification_k10 is not True
        or roots_agree is not True
        or clusters_agree is not True
    ):
        return "unresolved_spectrum"
    v_primary = float(V_primary)
    v_verification = float(V_verification)
    if not math.isfinite(v_primary) or not math.isfinite(v_verification):
        return "unresolved_spectrum"
    if (v_primary < 0.0 < v_verification) or (v_verification < 0.0 < v_primary):
        return "numerically_indeterminate_near_threshold"
    failure_primary = (
        int(N_true_primary) < int(required_prefix_n)
        or int(N_true_primary) < int(N_certified_0)
    )
    failure_verification = (
        int(N_true_verification) < int(required_prefix_n)
        or int(N_true_verification) < int(N_certified_0)
    )
    if min(v_primary, v_verification) > float(violation_tolerance) and failure_primary and failure_verification:
        return "confirmed_counterexample"
    spread = abs(v_primary - v_verification)
    if max(abs(v_primary), abs(v_verification)) <= max(float(violation_tolerance), spread):
        return "numerically_indeterminate_near_threshold"
    if max(v_primary, v_verification) <= 0.0:
        return "confirmed_safe_at_screen_point"
    return "numerically_indeterminate_near_threshold"


def same_geometry(left: Mapping[str, object], right: Mapping[str, object], *, tolerance: float = 1.0e-14) -> bool:
    return all(
        math.isclose(float(left[key]), float(right[key]), rel_tol=0.0, abs_tol=tolerance)
        for key in ("beta_deg", "mu", "eta")
    )


def worst_row(rows: Sequence[Mapping[str, object]]) -> Mapping[str, object] | None:
    finite = [row for row in rows if math.isfinite(float(row.get("V_n", float("nan"))))]
    if not finite:
        return None
    return max(finite, key=lambda row: (float(row["V_n"]), str(row.get("case_id", ""))))


def unique_case_ids(rows: Sequence[Mapping[str, object]]) -> tuple[str, ...]:
    return tuple(sorted({str(row["case_id"]) for row in rows if row.get("case_id")}))


def decide_step3a(
    *,
    geometry_count: int,
    resolved_geometry_count: int,
    baseline_controls_pass: bool,
    final_target_statuses: Sequence[str],
    strict_trigger_count: int,
    strict_verified_count: int,
) -> str:
    statuses = tuple(str(value) for value in final_target_statuses)
    if "confirmed_counterexample" in statuses:
        return "counterexample_found"
    if resolved_geometry_count < geometry_count or "unresolved_spectrum" in statuses or not baseline_controls_pass:
        return "inconclusive_due_to_unresolved_cases"
    if "numerically_indeterminate_near_threshold" in statuses:
        return "inconclusive_due_to_numerical_boundary"
    if (
        geometry_count == 28
        and resolved_geometry_count == 28
        and baseline_controls_pass
        and strict_verified_count == strict_trigger_count
    ):
        return "no_counterexample_in_28_case_screen"
    return "inconclusive_due_to_unresolved_cases"


__all__ = [
    "EXPECTED_PREFIX_GROUPS",
    "PrefixMetric",
    "certified_prefix",
    "decide_step3a",
    "final_status",
    "late_pass_indices",
    "needs_strict_verification",
    "primary_status",
    "running_prefix_metrics",
    "same_geometry",
    "squared_frequency_delta",
    "true_safe_prefix",
    "unique_case_ids",
    "worst_row",
]
