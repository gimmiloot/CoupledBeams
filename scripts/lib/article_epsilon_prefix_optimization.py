"""Staged, value-only spectrum orchestration for the article epsilon study.

The module deliberately reuses the unchanged general completeness search and
branch-informed strict gateway.  It contains no characteristic matrix, model
equation, mode-shape, energy, or predictor implementation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import time
from typing import Callable, Mapping, Sequence

from scripts.lib import article_epsilon_upper_envelope as workflow
from scripts.lib import branch_informed_spectrum_continuation as branch
from scripts.lib import general_spectrum_completeness as complete


PREFIX_ALGORITHM_VERSION = "article_prefix_staged_v5_floored_row_scaling"
PREFIX_CACHE_SCHEMA_VERSION = "article_prefix_cache_v6_gzip"
DEFAULT_STAGE_SAFE_ENDS = (4, 7, 10)
PREFIX_STRATEGIES = ("paired", "full-eb-progressive-timo")
STRICT_POLICIES = ("full", "local", "auto", "main-pass")
KNOWN_FULL_PATH_UNRESOLVED_SMOKE = (
    (0.020, 0.0, 0.7, -0.5),
    (0.050, 0.0, 0.7, -0.5),
    (0.050, 90.0, 0.0, 0.0),
    (0.020, 90.0, 0.7, -0.5),
)
EXECUTION_STATUSES = (
    "not_attempted",
    "resolved_prefix_early_stop",
    "resolved_full_K10",
    "attempted_unresolved",
    "attempted_error",
    "deferred_expensive_strict",
)


@dataclass(frozen=True)
class PrefixGuardAssessment:
    guard_index: int
    status: str
    resolved_through: int
    reasons: tuple[str, ...]
    close_cluster: bool
    root_ordering_warning: bool


@dataclass(frozen=True)
class StageTiming:
    model: str
    guard_index: int
    primary_seconds: float
    verification_seconds: float
    preparation_seconds: float
    primary_matrix_evaluations_added: int
    verification_matrix_evaluations_added: int
    primary_svd_calls_added: int
    verification_svd_calls_added: int


def validate_stage_safe_ends(values: Sequence[int]) -> tuple[int, ...]:
    stages = tuple(int(value) for value in values)
    if not stages or stages[-1] != workflow.K_MAX:
        raise ValueError("the final staged block must end at K=10")
    if any(value < 1 or value > workflow.K_MAX for value in stages):
        raise ValueError("staged safe-mode endpoints must lie in [1, 10]")
    if tuple(sorted(set(stages))) != stages:
        raise ValueError("staged safe-mode endpoints must be strictly increasing")
    return stages


def is_known_full_path_unresolved_smoke(point: workflow.GridPoint) -> bool:
    key = (point.epsilon_0, point.beta_deg, point.mu, point.eta)
    return key in KNOWN_FULL_PATH_UNRESOLVED_SMOKE


def stage_guard_indices(stage_safe_ends: Sequence[int] = DEFAULT_STAGE_SAFE_ENDS) -> tuple[int, ...]:
    return tuple(value + 1 for value in validate_stage_safe_ends(stage_safe_ends))


def staged_settings(guard_index: int) -> complete.SearchSettings:
    """Return unchanged tolerances with a smaller, adaptively grown range."""

    guard = int(guard_index)
    if not 2 <= guard <= workflow.ROOT_GUARD_INDEX:
        raise ValueError("guard_index must lie in [2, 11]")
    base = workflow.primary_settings()
    candidate_target = guard + 2
    return replace(
        base,
        requested_roots=guard,
        candidate_roots=candidate_target,
        verification_candidate_roots=candidate_target + 2,
        # This is only the first adaptive upper bound.  _run_configuration
        # grows it with the unchanged factor until its candidate target is met.
        lambda_max=max(6.0, 1.05 * float(candidate_target + 2)),
        max_upper_growth_tries=max(8, base.max_upper_growth_tries),
    )


def prefix_cache_identity(
    point: workflow.GridPoint,
    *,
    strategy: str,
    strict_policy: str,
    stage_safe_ends: Sequence[int] = DEFAULT_STAGE_SAFE_ENDS,
) -> dict[str, object]:
    if strategy not in PREFIX_STRATEGIES:
        raise ValueError(f"unsupported prefix strategy: {strategy}")
    if strict_policy not in STRICT_POLICIES:
        raise ValueError(f"unsupported strict policy: {strict_policy}")
    stages = validate_stage_safe_ends(stage_safe_ends)
    return {
        "cache_schema_version": PREFIX_CACHE_SCHEMA_VERSION,
        "algorithm_version": PREFIX_ALGORITHM_VERSION,
        "case_identity": point.case_identity,
        "full_precision_geometry": {
            "epsilon_0": float(point.epsilon_0).hex(),
            "beta_deg": float(point.beta_deg).hex(),
            "mu": float(point.mu).hex(),
            "eta": float(point.eta).hex(),
        },
        "scientific_model_configuration": workflow.solver_configuration(),
        "prefix_strategy": strategy,
        "strict_policy": strict_policy,
        "stage_safe_ends": list(stages),
    }


def _operations_copy(operations: complete.OperationCounts) -> complete.OperationCounts:
    return complete.OperationCounts(**asdict(operations))


def _operations_delta(
    after: complete.OperationCounts,
    before: complete.OperationCounts,
    name: str,
) -> int:
    return int(getattr(after, name)) - int(getattr(before, name))


def _evaluator_snapshot(evaluator: object) -> dict[str, object]:
    diagnostics = getattr(evaluator, "_diagnostic_cache")
    determinants = getattr(evaluator, "_determinant_cache")
    return {
        "operations": asdict(getattr(evaluator, "operations")),
        "diagnostics": [
            {"cache_key": float(key).hex(), "value": asdict(value)}
            for key, value in sorted(diagnostics.items())
        ],
        "determinants": [
            {"cache_key": float(key).hex(), "value": float(value)}
            for key, value in sorted(determinants.items())
        ],
    }


def _restore_evaluator(provider: complete.MatrixProvider, payload: Mapping[str, object]) -> object:
    operations = complete._operations_from_dict(payload.get("operations", {}))  # type: ignore[attr-defined,arg-type]
    evaluator = complete._MatrixEvaluator(provider, operations)  # type: ignore[attr-defined]
    diagnostic_cache = {
        float.fromhex(str(item["cache_key"])): complete._matrix_diagnostics_from_dict(item["value"])  # type: ignore[attr-defined,index,arg-type]
        for item in payload.get("diagnostics", ())  # type: ignore[union-attr]
    }
    determinant_cache = {
        float.fromhex(str(item["cache_key"])): float(item["value"])
        for item in payload.get("determinants", ())  # type: ignore[union-attr,index]
    }
    setattr(evaluator, "_diagnostic_cache", diagnostic_cache)
    setattr(evaluator, "_determinant_cache", determinant_cache)
    return evaluator


def _result_payload(result: complete.CompleteSpectrumResult) -> dict[str, object]:
    return {"result": asdict(result)}


def _result_from_payload(payload: Mapping[str, object]) -> complete.CompleteSpectrumResult:
    return complete._result_from_payload(payload, "partial_cache")  # type: ignore[attr-defined]


def _unresolved_below_guard(result: complete.CompleteSpectrumResult, guard_index: int) -> bool:
    roots = result.primary.roots
    if len(roots) < int(guard_index):
        return True
    upper = roots[int(guard_index) - 1].Lambda + max(
        result.settings.seed_half_width,
        2.0 * result.settings.scan_step,
    )
    for entry in (*result.primary.unresolved_intervals, *result.verification.unresolved_intervals):
        try:
            left = float(str(entry).split(":", 1)[0])
        except (TypeError, ValueError):
            return True
        if left <= upper:
            return True
    return False


def _cluster_contract(roots: Sequence[complete.RootRecord], guard_index: int) -> tuple[bool, bool]:
    relevant_ids = {root.root_cluster_id for root in roots[:guard_index] if root.root_cluster_id}
    close_cluster = bool(relevant_ids)
    for cluster_id in relevant_ids:
        members = [root for root in roots if root.root_cluster_id == cluster_id]
        if not members:
            return False, close_cluster
        expected_size = members[0].cluster_size
        if len(members) < expected_size:
            return False, close_cluster
        if sorted(root.cluster_member_index for root in members[:expected_size]) != list(range(1, expected_size + 1)):
            return False, close_cluster
        if any(root.cluster_size != expected_size for root in members[:expected_size]):
            return False, close_cluster
    return True, close_cluster


def assess_prefix_guard(
    result: complete.CompleteSpectrumResult,
    guard_index: int,
) -> PrefixGuardAssessment:
    guard = int(guard_index)
    reasons: list[str] = []
    primary = result.primary.roots
    verification = result.verification.roots
    if len(primary) < guard:
        reasons.append(f"primary_found_only_{len(primary)}_of_{guard}")
    if len(verification) < guard:
        reasons.append(f"verification_found_only_{len(verification)}_of_{guard}")
    comparisons = result.primary_vs_verification[:guard]
    if len(comparisons) < guard or any(row.get("status") != "pass" for row in comparisons):
        reasons.append("primary_strict_prefix_disagreement")
    relevant = primary[:guard]
    if any(
        root.acceptance_status != "accepted_full_matrix_svd"
        or not math.isfinite(root.sigma_1)
        or root.sigma_1 > result.settings.sigma_accept
        for root in relevant
    ):
        reasons.append("root_quality_failure")
    cluster_ok, close_cluster = _cluster_contract(primary, guard)
    verification_cluster_ok, verification_close_cluster = _cluster_contract(verification, guard)
    if not cluster_ok or not verification_cluster_ok:
        reasons.append("unresolved_cluster_multiplicity")
    if _unresolved_below_guard(result, guard):
        reasons.append("unresolved_interval_below_prefix_guard")
    ordering_warning = any(row.get("status") == "disagreement" for row in comparisons)
    resolved = not reasons
    return PrefixGuardAssessment(
        guard_index=guard,
        status="prefix_guard_resolved" if resolved else "prefix_guard_unresolved",
        resolved_through=guard - 1 if resolved else 0,
        reasons=tuple(dict.fromkeys(reasons)),
        close_cluster=bool(close_cluster or verification_close_cluster),
        root_ordering_warning=ordering_warning,
    )


class IncrementalModelSession:
    """Persistent in-memory SVD/determinant cache for one model and point."""

    def __init__(
        self,
        model: str,
        geometry: complete.Geometry,
        *,
        state: Mapping[str, object] | None = None,
    ) -> None:
        if model not in complete.SUPPORTED_MODELS:
            raise ValueError(f"unsupported model: {model}")
        self.model = model
        self.geometry = geometry
        self.provider = complete.model_matrix_provider(model, geometry)
        if state and isinstance(state.get("primary_evaluator"), Mapping):
            self.primary_evaluator = _restore_evaluator(self.provider, state["primary_evaluator"])  # type: ignore[arg-type]
        else:
            self.primary_evaluator = complete._MatrixEvaluator(  # type: ignore[attr-defined]
                self.provider, complete.OperationCounts()
            )
        if state and isinstance(state.get("verification_evaluator"), Mapping):
            self.verification_evaluator = _restore_evaluator(self.provider, state["verification_evaluator"])  # type: ignore[arg-type]
        else:
            self.verification_evaluator = complete._MatrixEvaluator(  # type: ignore[attr-defined]
                self.provider, complete.OperationCounts()
            )
        self.latest_result: complete.CompleteSpectrumResult | None = None
        if state and isinstance(state.get("latest_result"), Mapping):
            self.latest_result = _result_from_payload(state["latest_result"])  # type: ignore[arg-type]
        self.highest_resolved_mode = int(state.get("highest_resolved_mode", 0)) if state else 0
        self.highest_guard_mode = int(state.get("highest_guard_mode", 0)) if state else 0
        self.full_reference_completed = bool(state.get("full_reference_completed", False)) if state else False
        self.stage_timings: list[dict[str, object]] = list(state.get("stage_timings", ())) if state else []  # type: ignore[arg-type]

    @property
    def values(self) -> tuple[float, ...]:
        if self.latest_result is None:
            return ()
        return tuple(root.Lambda for root in self.latest_result.primary.roots)

    def snapshot(self) -> dict[str, object]:
        return {
            "model": self.model,
            "highest_resolved_mode": self.highest_resolved_mode,
            "highest_guard_mode": self.highest_guard_mode,
            "full_reference_completed": self.full_reference_completed,
            "resolved_roots": [root.Lambda for root in (self.latest_result.roots if self.latest_result else ())],
            "clusters": [
                {
                    "sorted_index": root.sorted_index,
                    "cluster_id": root.root_cluster_id,
                    "cluster_member_index": root.cluster_member_index,
                    "cluster_size": root.cluster_size,
                    "detected_nullity": root.detected_nullity,
                }
                for root in (self.latest_result.primary.roots if self.latest_result else ())
                if root.root_cluster_id or root.detected_nullity > 1
            ],
            "brackets_and_continuation_metadata": [
                dict(row)
                for row in (self.latest_result.primary.interval_rows if self.latest_result else ())
            ],
            "latest_result": _result_payload(self.latest_result) if self.latest_result else {},
            "primary_evaluator": _evaluator_snapshot(self.primary_evaluator),
            "verification_evaluator": _evaluator_snapshot(self.verification_evaluator),
            "stage_timings": list(self.stage_timings),
        }

    def extend_to_guard(
        self,
        guard_index: int,
        *,
        continuation_seeds: Sequence[float] = (),
        eb_seed_roots: Sequence[float] = (),
        full_reference: bool = False,
    ) -> tuple[complete.CompleteSpectrumResult, PrefixGuardAssessment, StageTiming, bool]:
        guard = int(guard_index)
        if (
            self.latest_result is not None
            and self.highest_guard_mode >= guard
            and (not full_reference or self.full_reference_completed)
        ):
            assessment = assess_prefix_guard(self.latest_result, guard)
            return (
                self.latest_result,
                assessment,
                StageTiming(self.model, guard, 0.0, 0.0, 0.0, 0, 0, 0, 0),
                True,
            )
        prepared_at = time.perf_counter()
        settings = workflow.primary_settings() if full_reference else staged_settings(guard)
        seeds = [float(value) for value in continuation_seeds if math.isfinite(float(value))]
        seed_source = "continuation_seed_window"
        if abs(self.geometry.beta_deg) <= 1.0e-14 and abs(self.geometry.eta) <= 1.0e-14:
            seeds.extend(complete.straight_oracle_values(self.model, self.geometry, settings.candidate_roots))
            seed_source = "straight_oracle_seed"
        if self.model == complete.MODEL_TIMO and eb_seed_roots:
            seeds.extend(float(value) for value in eb_seed_roots if math.isfinite(float(value)))
            seed_source = "EB_seed_window_for_Timo"
        if self.latest_result is not None:
            seeds.extend(root.Lambda for root in self.latest_result.primary.roots)
        preparation_seconds = time.perf_counter() - prepared_at

        primary_before = _operations_copy(self.primary_evaluator.operations)
        primary_at = time.perf_counter()
        primary = complete._run_configuration(  # type: ignore[attr-defined]
            self.provider,
            settings,
            configuration="primary",
            candidate_target=settings.candidate_roots,
            scan_step=settings.scan_step,
            phases=(0.0, settings.shifted_grid_phase),
            seed_roots=tuple(seeds),
            seed_source=seed_source,
            evaluator=self.primary_evaluator,
        )
        primary_seconds = time.perf_counter() - primary_at

        verification_before = _operations_copy(self.verification_evaluator.operations)
        verification_at = time.perf_counter()
        verification = complete._run_configuration(  # type: ignore[attr-defined]
            self.provider,
            settings,
            configuration="verification",
            candidate_target=settings.verification_candidate_roots,
            scan_step=0.5 * settings.scan_step,
            phases=(settings.shifted_grid_phase,),
            seed_roots=(),
            seed_source="independent_no_primary_objects",
            evaluator=self.verification_evaluator,
        )
        primary, verification = complete.reconcile_independent_configurations(
            primary,
            verification,
            self.primary_evaluator,
            self.verification_evaluator,
            settings,
        )
        verification_seconds = time.perf_counter() - verification_at
        comparison, agreement = complete._compare_configurations(primary, verification, settings)  # type: ignore[attr-defined]
        roots = tuple(primary.roots[:guard])
        operations = complete.OperationCounts()
        operations.add(primary.operations)
        operations.add(verification.operations)
        requested_count = settings.requested_roots
        roots = tuple(primary.roots[:requested_count])
        result = complete.CompleteSpectrumResult(
            algorithm_version=f"{complete.GENERAL_SPECTRUM_ALGORITHM_VERSION}+{PREFIX_ALGORITHM_VERSION}",
            model=self.model,
            geometry=self.geometry,
            settings=settings,
            primary=primary,
            verification=verification,
            roots=roots,
            primary_vs_verification=comparison,
            independent_agreement=agreement,
            root11_available=len(roots) >= 11,
            root12_available=len(roots) >= 12,
            root12_boundary_warning=(
                len(primary.roots) < settings.candidate_roots
                or len(verification.roots) < settings.verification_candidate_roots
            ),
            spectrum_status="resolved_complete" if agreement and len(roots) >= requested_count else "unresolved",
            exclusion_reason="" if agreement and len(roots) >= requested_count else "staged_prefix_incomplete",
            cache_status="miss",
            operations=operations,
        )
        result = complete.annotate_coalesced_tracks(result, continuation_seeds, settings)
        assessment = assess_prefix_guard(result, guard)
        if assessment.status == "prefix_guard_resolved":
            self.highest_resolved_mode = max(self.highest_resolved_mode, guard - 1)
            self.highest_guard_mode = max(self.highest_guard_mode, guard)
        if full_reference:
            self.full_reference_completed = True
        self.latest_result = result
        timing = StageTiming(
            model=self.model,
            guard_index=guard,
            primary_seconds=primary_seconds,
            verification_seconds=verification_seconds,
            preparation_seconds=preparation_seconds,
            primary_matrix_evaluations_added=_operations_delta(
                self.primary_evaluator.operations, primary_before, "characteristic_matrix_evaluations"
            ),
            verification_matrix_evaluations_added=_operations_delta(
                self.verification_evaluator.operations, verification_before, "characteristic_matrix_evaluations"
            ),
            primary_svd_calls_added=_operations_delta(
                self.primary_evaluator.operations, primary_before, "full_6x6_SVD_calls"
            ),
            verification_svd_calls_added=_operations_delta(
                self.verification_evaluator.operations, verification_before, "full_6x6_SVD_calls"
            ),
        )
        self.stage_timings.append(asdict(timing))
        return result, assessment, timing, False


class PartialPointCache:
    def __init__(self, cache_dir: Path, *, reuse_cache: bool = True, force: bool = False) -> None:
        self.cache_dir = Path(cache_dir)
        self.reuse_cache = bool(reuse_cache)
        self.force = bool(force)
        self.last_load_status = "not_checked"

    def path(
        self,
        point: workflow.GridPoint,
        *,
        strategy: str,
        strict_policy: str,
        stage_safe_ends: Sequence[int] = DEFAULT_STAGE_SAFE_ENDS,
    ) -> Path:
        identity = prefix_cache_identity(
            point,
            strategy=strategy,
            strict_policy=strict_policy,
            stage_safe_ends=stage_safe_ends,
        )
        digest = hashlib.sha256(
            json.dumps(identity, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()[:24]
        return self.cache_dir / f"prefix_case_{point.case_id}_{digest}.json.gz"

    def load(
        self,
        point: workflow.GridPoint,
        *,
        strategy: str,
        strict_policy: str,
        stage_safe_ends: Sequence[int] = DEFAULT_STAGE_SAFE_ENDS,
    ) -> dict[str, object] | None:
        if self.force or not self.reuse_cache:
            self.last_load_status = "force" if self.force else "disabled"
            return None
        path = self.path(
            point,
            strategy=strategy,
            strict_policy=strict_policy,
            stage_safe_ends=stage_safe_ends,
        )
        if not path.exists():
            self.last_load_status = "miss"
            return None
        try:
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                payload = json.load(handle)
        except (OSError, EOFError, json.JSONDecodeError):
            self.last_load_status = "invalid_payload"
            return None
        identity = prefix_cache_identity(
            point,
            strategy=strategy,
            strict_policy=strict_policy,
            stage_safe_ends=stage_safe_ends,
        )
        if payload.get("cache_schema_version") != PREFIX_CACHE_SCHEMA_VERSION:
            self.last_load_status = "incompatible_schema"
            return None
        if payload.get("identity") != identity:
            self.last_load_status = "incompatible_identity"
            return None
        self.last_load_status = "hit"
        return payload

    def save(
        self,
        point: workflow.GridPoint,
        payload: Mapping[str, object],
        *,
        strategy: str,
        strict_policy: str,
        stage_safe_ends: Sequence[int] = DEFAULT_STAGE_SAFE_ENDS,
    ) -> Path:
        path = self.path(
            point,
            strategy=strategy,
            strict_policy=strict_policy,
            stage_safe_ends=stage_safe_ends,
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_name(path.name + ".tmp")
        encoded = json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=True,
        ).encode("utf-8")
        temporary.write_bytes(gzip.compress(encoded, compresslevel=3))
        os.replace(temporary, path)
        return path


StrictCallback = Callable[
    [workflow.GridPoint, Mapping[str, IncrementalModelSession], int],
    Mapping[str, object],
]


def _persist_sessions(
    payload: dict[str, object],
    sessions: Mapping[str, IncrementalModelSession],
) -> None:
    payload["models"] = {model: session.snapshot() for model, session in sessions.items()}
    payload["highest_resolved_mode_eb"] = sessions[complete.MODEL_EB].highest_resolved_mode
    payload["highest_guard_mode_eb"] = sessions[complete.MODEL_EB].highest_guard_mode
    payload["highest_resolved_mode_timo"] = sessions[complete.MODEL_TIMO].highest_resolved_mode
    payload["highest_guard_mode_timo"] = sessions[complete.MODEL_TIMO].highest_guard_mode


def _save_partial(
    cache: PartialPointCache,
    point: workflow.GridPoint,
    payload: dict[str, object],
    sessions: Mapping[str, IncrementalModelSession],
    *,
    strategy: str,
    strict_policy: str,
    stage_safe_ends: Sequence[int],
) -> None:
    _persist_sessions(payload, sessions)
    cache.save(
        point,
        payload,
        strategy=strategy,
        strict_policy=strict_policy,
        stage_safe_ends=stage_safe_ends,
    )


def _stage_assessments(
    sessions: Mapping[str, IncrementalModelSession],
    guard_index: int,
) -> dict[str, PrefixGuardAssessment]:
    result: dict[str, PrefixGuardAssessment] = {}
    for model, session in sessions.items():
        if session.latest_result is None:
            result[model] = PrefixGuardAssessment(
                int(guard_index),
                "prefix_guard_unresolved",
                0,
                ("model_not_attempted",),
                False,
                False,
            )
        else:
            result[model] = assess_prefix_guard(session.latest_result, int(guard_index))
    return result


def _append_timing(payload: dict[str, object], timing: StageTiming, cache_hit: bool) -> None:
    rows = payload.setdefault("stage_timings", [])
    if isinstance(rows, list):
        rows.append({**asdict(timing), "cache_hit": bool(cache_hit)})


def run_staged_point(
    point: workflow.GridPoint,
    *,
    cache: PartialPointCache,
    strategy: str = "paired",
    strict_policy: str = "local",
    full_k10: bool = False,
    stage_safe_ends: Sequence[int] = DEFAULT_STAGE_SAFE_ENDS,
    continuation_roots: Mapping[str, Sequence[float]] | None = None,
    strict_callback: StrictCallback | None = None,
    force_audit: bool = False,
) -> dict[str, object]:
    """Resolve one point, persisting after every successful model stage.

    ``full_k10=True`` deliberately ignores an earlier cached stop and extends
    the same evaluator caches through root 11.  It therefore supplies the
    requested energy/article follow-up continuation contract without silently
    recomputing the already cached prefix.
    """

    point_started = time.perf_counter()
    if strategy not in PREFIX_STRATEGIES:
        raise ValueError(f"unsupported prefix strategy: {strategy}")
    if strict_policy not in STRICT_POLICIES:
        raise ValueError(f"unsupported strict policy: {strict_policy}")
    stages = validate_stage_safe_ends(stage_safe_ends)
    loaded = cache.load(
        point,
        strategy=strategy,
        strict_policy=strict_policy,
        stage_safe_ends=stages,
    )
    if loaded is None:
        payload = fresh_point_payload(
            point,
            strategy=strategy,
            strict_policy=strict_policy,
            stage_safe_ends=stages,
        )
    else:
        payload = dict(loaded)
        payload["cache_hits"] = int(payload.get("cache_hits", 0)) + 1
        if payload.get("execution_status") == "deferred_expensive_strict":
            return payload
        if (
            not full_k10
            and payload.get("execution_status") == "resolved_prefix_early_stop"
            and payload.get("N_true") is not None
        ):
            return payload
        if full_k10 and payload.get("execution_status") == "resolved_full_K10":
            return payload

    model_states = payload.get("models", {})
    sessions = {
        model: IncrementalModelSession(
            model,
            point.geometry,
            state=(model_states.get(model) if isinstance(model_states, Mapping) else None),  # type: ignore[arg-type,union-attr]
        )
        for model in complete.SUPPORTED_MODELS
    }
    continuation = continuation_roots or {}
    # Persisted stage snapshots written before the final decision must never
    # look like a point that was not sent to a solver.  If the process is
    # interrupted, postprocessing will classify this as an attempted case.
    payload["execution_status"] = "attempted_unresolved"
    payload["unresolved_reason"] = "execution_in_progress_or_interrupted"
    payload["run_mode"] = "full-k10" if full_k10 else "prefix-until-failure"
    deltas = {int(key): float(value) for key, value in dict(payload.get("deltas", {})).items()}
    first_failed_mode: int | None = None
    first_failed_delta: float | None = None
    last_assessments: dict[str, PrefixGuardAssessment] = {}
    last_guard = 0

    try:
        if full_k10:
            eb_result, _eb_assessment, eb_timing, eb_hit = sessions[complete.MODEL_EB].extend_to_guard(
                workflow.ROOT_GUARD_INDEX,
                continuation_seeds=continuation.get(complete.MODEL_EB, ()),
                full_reference=True,
            )
            _append_timing(payload, eb_timing, eb_hit)
            _save_partial(
                cache,
                point,
                payload,
                sessions,
                strategy=strategy,
                strict_policy=strict_policy,
                stage_safe_ends=stages,
            )
            tm_result, _tm_assessment, tm_timing, tm_hit = sessions[complete.MODEL_TIMO].extend_to_guard(
                workflow.ROOT_GUARD_INDEX,
                continuation_seeds=continuation.get(complete.MODEL_TIMO, ()),
                eb_seed_roots=tuple(root.Lambda for root in eb_result.primary.roots),
                full_reference=True,
            )
            _append_timing(payload, tm_timing, tm_hit)
            last_guard = workflow.ROOT_GUARD_INDEX
            last_assessments = _stage_assessments(sessions, last_guard)
            _save_partial(
                cache,
                point,
                payload,
                sessions,
                strategy=strategy,
                strict_policy=strict_policy,
                stage_safe_ends=stages,
            )
            if any(item.status != "prefix_guard_resolved" for item in last_assessments.values()):
                payload.update(
                    {
                        "execution_status": "attempted_unresolved",
                        "N_true": None,
                        "prefix_guard_status": "prefix_guard_unresolved",
                        "full_K10_guard_status": "full_K10_guard_unresolved",
                        "unresolved_reason": ";".join(
                            f"{model}:{reason}"
                            for model, item in last_assessments.items()
                            for reason in item.reasons
                        ),
                    }
                )
                _save_partial(
                    cache,
                    point,
                    payload,
                    sessions,
                    strategy=strategy,
                    strict_policy=strict_policy,
                    stage_safe_ends=stages,
                )
                return payload
            eb_values = tuple(root.Lambda for root in eb_result.primary.roots)
            tm_values = tuple(root.Lambda for root in tm_result.primary.roots)
            deltas = {
                mode: workflow.squared_frequency_delta(eb_values[mode - 1], tm_values[mode - 1])
                for mode in range(1, workflow.K_MAX + 1)
            }
            payload["deltas"] = {str(key): value for key, value in sorted(deltas.items())}

        if strategy == "full-eb-progressive-timo" and not full_k10:
            eb_result, eb_assessment, timing, hit = sessions[complete.MODEL_EB].extend_to_guard(
                workflow.ROOT_GUARD_INDEX,
                continuation_seeds=continuation.get(complete.MODEL_EB, ()),
            )
            _append_timing(payload, timing, hit)
            _save_partial(
                cache,
                point,
                payload,
                sessions,
                strategy=strategy,
                strict_policy=strict_policy,
                stage_safe_ends=stages,
            )
            if eb_assessment.status != "prefix_guard_resolved":
                payload.update(
                    {
                        "execution_status": "attempted_unresolved",
                        "prefix_guard_status": "prefix_guard_unresolved",
                        "unresolved_reason": ";".join(eb_assessment.reasons),
                    }
                )
                _save_partial(
                    cache,
                    point,
                    payload,
                    sessions,
                    strategy=strategy,
                    strict_policy=strict_policy,
                    stage_safe_ends=stages,
                )
                return payload

        previous_safe_end = 0
        for safe_end in (() if full_k10 else stages):
            guard = safe_end + 1
            last_guard = guard
            if strategy == "paired":
                eb_result, _eb_assessment, eb_timing, eb_hit = sessions[complete.MODEL_EB].extend_to_guard(
                    guard,
                    continuation_seeds=continuation.get(complete.MODEL_EB, ()),
                )
                _append_timing(payload, eb_timing, eb_hit)
                _save_partial(
                    cache,
                    point,
                    payload,
                    sessions,
                    strategy=strategy,
                    strict_policy=strict_policy,
                    stage_safe_ends=stages,
                )
            else:
                eb_result = sessions[complete.MODEL_EB].latest_result
                if eb_result is None:  # pragma: no cover - defensive contract
                    raise RuntimeError("full EB stage was not prepared")
            tm_result, _tm_assessment, tm_timing, tm_hit = sessions[complete.MODEL_TIMO].extend_to_guard(
                guard,
                continuation_seeds=continuation.get(complete.MODEL_TIMO, ()),
                eb_seed_roots=tuple(root.Lambda for root in eb_result.primary.roots),
            )
            _append_timing(payload, tm_timing, tm_hit)
            _save_partial(
                cache,
                point,
                payload,
                sessions,
                strategy=strategy,
                strict_policy=strict_policy,
                stage_safe_ends=stages,
            )
            last_assessments = _stage_assessments(sessions, guard)
            if any(item.status != "prefix_guard_resolved" for item in last_assessments.values()):
                reasons = [
                    f"{model}:{reason}"
                    for model, item in last_assessments.items()
                    for reason in item.reasons
                ]
                payload.update(
                    {
                        "execution_status": "attempted_unresolved",
                        "N_true": None,
                        "prefix_guard_status": "prefix_guard_unresolved",
                        "prefix_guard_resolved_through": 0,
                        "full_K10_guard_status": (
                            "full_K10_guard_unresolved" if guard == workflow.ROOT_GUARD_INDEX else "not_attempted"
                        ),
                        "unresolved_reason": ";".join(reasons),
                    }
                )
                _save_partial(
                    cache,
                    point,
                    payload,
                    sessions,
                    strategy=strategy,
                    strict_policy=strict_policy,
                    stage_safe_ends=stages,
                )
                return payload

            eb_values = tuple(root.Lambda for root in eb_result.primary.roots)
            tm_values = tuple(root.Lambda for root in tm_result.primary.roots)
            for mode in range(previous_safe_end + 1, safe_end + 1):
                deltas[mode] = workflow.squared_frequency_delta(eb_values[mode - 1], tm_values[mode - 1])
                if first_failed_mode is None and deltas[mode] > workflow.DELTA_TOL:
                    first_failed_mode = mode
                    first_failed_delta = deltas[mode]
                    break
            payload["deltas"] = {str(key): value for key, value in sorted(deltas.items())}
            if first_failed_mode is not None and not full_k10:
                failure_guard = first_failed_mode + 1
                failure_assessments = _stage_assessments(sessions, failure_guard)
                if any(item.status != "prefix_guard_resolved" for item in failure_assessments.values()):
                    payload.update(
                        {
                            "execution_status": "attempted_unresolved",
                            "N_true": None,
                            "prefix_guard_status": "prefix_guard_unresolved",
                            "unresolved_reason": "failure_right_guard_unresolved",
                        }
                    )
                    _save_partial(
                        cache,
                        point,
                        payload,
                        sessions,
                        strategy=strategy,
                        strict_policy=strict_policy,
                        stage_safe_ends=stages,
                    )
                    return payload
                last_assessments = failure_assessments
                last_guard = failure_guard
                break
            previous_safe_end = safe_end

        if full_k10:
            if len(deltas) < workflow.K_MAX:
                payload.update(
                    {
                        "execution_status": "attempted_unresolved",
                        "N_true": None,
                        "full_K10_guard_status": "full_K10_guard_unresolved",
                        "unresolved_reason": "not_all_K10_deltas_resolved",
                    }
                )
                _save_partial(
                    cache,
                    point,
                    payload,
                    sessions,
                    strategy=strategy,
                    strict_policy=strict_policy,
                    stage_safe_ends=stages,
                )
                return payload
            first_failed_mode = next(
                (mode for mode in range(1, workflow.K_MAX + 1) if deltas[mode] > workflow.DELTA_TOL),
                None,
            )
            first_failed_delta = deltas.get(first_failed_mode) if first_failed_mode else None

        n_true = workflow.K_MAX if first_failed_mode is None else first_failed_mode - 1
        trigger_reasons = strict_trigger_reasons(
            first_failed_mode=first_failed_mode,
            deltas=deltas,
            assessments=last_assessments,
            regression_label=point.regression_label,
            force_audit=force_audit,
        )
        effective_policy = choose_effective_strict_policy(strict_policy, trigger_reasons)
        verified_modes = local_strict_modes(first_failed_mode, last_guard)
        strict_status = "local_prefix_count_and_guard_pass"
        required_prefix_strict_status = "pass"
        optional_upper_spectrum_full_audit_status = "not_requested"
        strict_payload: Mapping[str, object] = {}
        if effective_policy == "full":
            payload["force_strict_requested"] = int(payload.get("force_strict_requested", 0)) + 1
            if strict_policy == "main-pass":
                roots_resolved: dict[str, list[float]] = {}
                cluster_status: dict[str, object] = {}
                primary_status: dict[str, object] = {}
                for model, session in sessions.items():
                    state = session.latest_result
                    roots_resolved[model] = (
                        [float(root.Lambda) for root in state.primary.roots]
                        if state is not None
                        else []
                    )
                    assessment = last_assessments.get(model)
                    cluster_status[model] = {
                        "close_cluster": bool(assessment.close_cluster) if assessment else False,
                        "root_ordering_warning": (
                            bool(assessment.root_ordering_warning) if assessment else False
                        ),
                        "guard_status": assessment.status if assessment else "not_available",
                    }
                    primary_status[model] = (
                        getattr(state, "spectrum_status", "available")
                        if state is not None
                        else "not_available"
                    )
                payload.update(
                    {
                        "execution_status": "deferred_expensive_strict",
                        "n_true_status": "unresolved_pending_complex_pass",
                        "N_true": None,
                        "expensive_escalation_requested": True,
                        "expensive_escalation_kind": "force_strict_verification",
                        "trigger_reason": list(trigger_reasons),
                        "first_apparent_failed_mode": first_failed_mode,
                        "first_failed_mode": first_failed_mode,
                        "first_failed_delta_f": first_failed_delta,
                        "required_guard_mode": int(last_guard),
                        "roots_already_resolved": roots_resolved,
                        "clusters_already_resolved": cluster_status,
                        "primary_status": primary_status,
                        "local_status": "local_prefix_count_and_guard_pass",
                        "elapsed_time_before_escalation": time.perf_counter() - point_started,
                        "strict_policy_effective": "deferred_before_full",
                        "strict_trigger_reasons": list(trigger_reasons),
                        "strict_verification_status": "not_executed_deferred",
                        "required_prefix_strict_status": "unresolved_without_expensive_strict",
                        "prefix_guard_status": "unresolved_without_expensive_strict",
                        "required_prefix_guard_status": "unresolved_without_expensive_strict",
                        "upper_spectrum_audit_status": "not_attempted",
                        "optional_upper_spectrum_full_audit_status": "not_attempted",
                        "full_K10_control_status": "not_attempted",
                        "full_K10_guard_status": "not_attempted",
                        "defer_reason": "force_strict_verification_required",
                        "deferred_before_force_strict": int(
                            payload.get("deferred_before_force_strict", 0)
                        )
                        + 1,
                        "force_strict_executed": int(payload.get("force_strict_executed", 0)),
                        "unresolved_reason": "force_strict_verification_required",
                    }
                )
                _save_partial(
                    cache,
                    point,
                    payload,
                    sessions,
                    strategy=strategy,
                    strict_policy=strict_policy,
                    stage_safe_ends=stages,
                )
                return payload
            if strict_callback is None:
                strict_status = "full_strict_not_available"
                required_prefix_strict_status = "fail"
                optional_upper_spectrum_full_audit_status = "not_available"
            else:
                payload["force_strict_executed"] = int(
                    payload.get("force_strict_executed", 0)
                ) + 1
                strict_payload = strict_callback(point, sessions, last_guard)
                strict_status = str(strict_payload.get("status", "full_strict_failed"))
                required_prefix_strict_status = str(
                    strict_payload.get(
                        "required_prefix_strict_status",
                        "pass" if strict_status in {"pass", "full_strict_pass"} else "fail",
                    )
                )
                optional_upper_spectrum_full_audit_status = str(
                    strict_payload.get(
                        "optional_upper_spectrum_full_audit_status",
                        "pass" if strict_status in {"pass", "full_strict_pass"} else "fail",
                    )
                )
            if required_prefix_strict_status != "pass":
                payload.update(
                    {
                        "execution_status": "attempted_unresolved",
                        "N_true": None,
                        "prefix_guard_status": "prefix_guard_unresolved",
                        "strict_verification_status": strict_status,
                        "required_prefix_strict_status": required_prefix_strict_status,
                        "optional_upper_spectrum_full_audit_status": (
                            optional_upper_spectrum_full_audit_status
                        ),
                        "strict_details": dict(strict_payload),
                        "unresolved_reason": "strict_failure_or_disagreement",
                    }
                )
                _save_partial(
                    cache,
                    point,
                    payload,
                    sessions,
                    strategy=strategy,
                    strict_policy=strict_policy,
                    stage_safe_ends=stages,
                )
                return payload
            verified_modes = tuple(range(1, last_guard + 1))

        payload.update(
            {
                "execution_status": (
                    "resolved_full_K10" if full_k10 or n_true == workflow.K_MAX else "resolved_prefix_early_stop"
                ),
                "N_true": int(n_true),
                "first_failed_mode": first_failed_mode,
                "first_failed_delta_f": first_failed_delta,
                "prefix_guard_status": "prefix_guard_resolved",
                "prefix_guard_resolved_through": (
                    workflow.K_MAX if first_failed_mode is None else int(first_failed_mode)
                ),
                "full_K10_guard_status": (
                    "full_K10_guard_resolved"
                    if full_k10 or n_true == workflow.K_MAX
                    else "not_needed_after_first_failure"
                ),
                "early_stop_used": bool(not full_k10 and n_true < workflow.K_MAX),
                "early_stop_reason": (
                    "confirmed_first_delta_failure" if not full_k10 and n_true < workflow.K_MAX else ""
                ),
                "strict_policy_effective": effective_policy,
                "strict_trigger_reasons": list(trigger_reasons),
                "strict_verification_status": strict_status,
                "required_prefix_strict_status": required_prefix_strict_status,
                "optional_upper_spectrum_full_audit_status": (
                    optional_upper_spectrum_full_audit_status
                ),
                "strict_verified_modes": list(verified_modes),
                "strict_details": dict(strict_payload),
                "unresolved_reason": "",
            }
        )
        _save_partial(
            cache,
            point,
            payload,
            sessions,
            strategy=strategy,
            strict_policy=strict_policy,
            stage_safe_ends=stages,
        )
        return payload
    except Exception as exc:
        payload.update(
            {
                "execution_status": "attempted_error",
                "N_true": None,
                "prefix_guard_status": "prefix_guard_unresolved",
                "unresolved_reason": f"{type(exc).__name__}: {exc}",
            }
        )
        _save_partial(
            cache,
            point,
            payload,
            sessions,
            strategy=strategy,
            strict_policy=strict_policy,
            stage_safe_ends=stages,
        )
        return payload


def fresh_point_payload(
    point: workflow.GridPoint,
    *,
    strategy: str,
    strict_policy: str,
    stage_safe_ends: Sequence[int] = DEFAULT_STAGE_SAFE_ENDS,
) -> dict[str, object]:
    identity = prefix_cache_identity(
        point,
        strategy=strategy,
        strict_policy=strict_policy,
        stage_safe_ends=stage_safe_ends,
    )
    return {
        "cache_schema_version": PREFIX_CACHE_SCHEMA_VERSION,
        "identity": identity,
        "solver_configuration_hash": hashlib.sha256(
            json.dumps(identity, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "execution_status": "not_attempted",
        "models": {},
        "first_failed_mode": None,
        "first_failed_delta_f": None,
        "N_true": None,
        "deltas": {},
        "prefix_guard_status": "not_attempted",
        "prefix_guard_resolved_through": 0,
        "full_K10_guard_status": "not_attempted",
        "early_stop_used": False,
        "early_stop_reason": "",
        "strict_verification_status": "not_attempted",
        "required_prefix_strict_status": "not_attempted",
        "optional_upper_spectrum_full_audit_status": "not_attempted",
        "strict_verified_modes": [],
        "stage_timings": [],
        "cache_hits": 0,
        "force_strict_requested": 0,
        "force_strict_executed": 0,
        "deferred_before_force_strict": 0,
    }


def local_strict_modes(first_failed_mode: int | None, guard_index: int) -> tuple[int, ...]:
    """Local strict covers the complete root count plus k-1/k/k+1 context."""

    if first_failed_mode is None:
        return tuple(range(1, int(guard_index) + 1))
    k = int(first_failed_mode)
    local = set(range(1, k + 1))
    local.update(range(max(1, k - 1), min(int(guard_index), k + 1) + 1))
    return tuple(sorted(local))


def assess_branch_strict_through_guard(
    result: branch.BranchContinuationResult | None,
    primary_values: Sequence[float],
    guard_index: int,
    tolerance: float,
) -> dict[str, object]:
    """Assess strict evidence only through the scientifically required guard.

    The branch solver may continue beyond ``guard_index`` for its own full-K10
    audit.  A localized failure to the right of that guard remains visible in
    the full-audit status, but it must not invalidate an independently complete
    prefix.
    """

    guard = int(guard_index)
    if result is None:
        return {
            "cache_status": "miss",
            "agreement": False,
            "first_disagreement": None,
            "cluster_resolved": False,
            "quality_resolved": False,
            "boundary_resolved": False,
            "reasons": ("strict_result_missing",),
        }
    ordered = tuple(sorted(result.branches, key=lambda item: item.Lambda))
    strict_values = tuple(item.Lambda for item in ordered)
    disagreements: list[int] = []
    reasons: list[str] = []
    if len(primary_values) < guard or len(strict_values) < guard:
        disagreements.append(min(len(primary_values), len(strict_values)) + 1)
        reasons.append("insufficient_root_count_through_guard")
    for index in range(min(guard, len(primary_values), len(strict_values))):
        if abs(float(primary_values[index]) - strict_values[index]) > float(tolerance):
            disagreements.append(index + 1)
            reasons.append("primary_strict_root_disagreement")
    comparisons = tuple(result.primary_vs_verification[:guard])
    for index, comparison in enumerate(comparisons, start=1):
        if not bool(comparison.get("within_tolerance")):
            disagreements.append(index)
            reasons.append("strict_independent_root_disagreement")
    if len(comparisons) < guard:
        disagreements.append(len(comparisons) + 1)
        reasons.append("strict_independent_comparison_incomplete")
    quality_failure = next(
        (
            index
            for index, item in enumerate(ordered[:guard], start=1)
            if item.sigma_1 > result.settings.sigma_accept
            or item.sigma_ratio > result.settings.sigma_ratio_accept
        ),
        None,
    )
    if quality_failure is not None:
        disagreements.append(quality_failure)
        reasons.append("strict_root_quality_failure")
    relevant_ids = {item.branch_id for item in ordered[:guard]}
    cluster_resolved = all(
        cluster.resolved
        for step in result.steps
        if step.accepted
        for cluster in step.clusters
        if relevant_ids.intersection(cluster.branch_ids)
    )
    if not cluster_resolved:
        disagreements.append(guard)
        reasons.append("strict_cluster_multiplicity_unresolved")
    guard_lambda = strict_values[guard - 1] if len(strict_values) >= guard else float("nan")
    boundary_resolved = math.isfinite(guard_lambda)
    if boundary_resolved:
        margin = max(float(tolerance), result.settings.seed_half_width)
        if any(float(value) <= guard_lambda + margin for value in result.guard.unmatched_candidates):
            boundary_resolved = False
        for entry in result.guard.unresolved_intervals:
            try:
                left = float(str(entry).split(":", 1)[0])
            except (TypeError, ValueError):
                boundary_resolved = False
                break
            if left <= guard_lambda + margin:
                boundary_resolved = False
                break
    if not boundary_resolved:
        disagreements.append(guard)
        reasons.append("strict_boundary_unresolved_through_guard")
    first = min(disagreements) if disagreements else None
    return {
        "cache_status": "hit",
        "agreement": first is None,
        "first_disagreement": first,
        "cluster_resolved": cluster_resolved,
        "quality_resolved": quality_failure is None,
        "boundary_resolved": boundary_resolved,
        "reasons": tuple(dict.fromkeys(reasons)),
        "k10_guard_resolved": result.k10_guard_resolved,
        "exclusion_reason": result.exclusion_reason,
    }


def strict_trigger_reasons(
    *,
    first_failed_mode: int | None,
    deltas: Mapping[int, float],
    assessments: Mapping[str, PrefixGuardAssessment],
    regression_label: str = "",
    force_audit: bool = False,
) -> tuple[str, ...]:
    reasons: list[str] = []
    if first_failed_mode is not None:
        reasons.append("first_failed_mode")
        value = deltas.get(int(first_failed_mode), float("nan"))
        if math.isfinite(value) and 0.095 <= value <= 0.105:
            reasons.append("near_delta_threshold")
        if first_failed_mode > 1:
            previous = deltas.get(first_failed_mode - 1, float("nan"))
            if math.isfinite(previous) and 0.095 <= previous <= 0.105:
                reasons.append("last_safe_mode_near_threshold")
    if any(item.close_cluster for item in assessments.values()):
        reasons.append("close_root_cluster")
    if any(item.root_ordering_warning for item in assessments.values()):
        reasons.append("root_ordering_warning")
    if regression_label:
        reasons.append("exact_regression_point")
    if force_audit:
        reasons.append("deterministic_control_case")
    return tuple(dict.fromkeys(reasons))


def choose_effective_strict_policy(
    requested: str,
    trigger_reasons: Sequence[str],
) -> str:
    if requested not in STRICT_POLICIES:
        raise ValueError(f"unsupported strict policy: {requested}")
    # A declared validation/control case remains a full audit even in a local
    # benchmark arm.  Otherwise a known full-spectrum failure to the right of
    # an early prefix guard could be hidden by a scientifically valid local
    # stop, contradicting the validation contract.
    if "deterministic_control_case" in trigger_reasons:
        return "full"
    if requested not in {"auto", "main-pass"}:
        return requested
    full_reasons = {
        "exact_regression_point",
        "deterministic_control_case",
        "close_root_cluster",
        "root_ordering_warning",
    }
    return "full" if full_reasons.intersection(trigger_reasons) else "local"


def envelope_saturation_status(n_up_raw: int | None) -> dict[str, object]:
    saturated = n_up_raw == workflow.K_MAX
    return {
        "complete_for_upper_envelope": bool(saturated),
        "complete_for_distribution": False,
        "remaining_status": "not_needed_after_envelope_saturation" if saturated else "not_attempted",
    }


def deterministic_control_sample(
    rows: Sequence[Mapping[str, object]],
    *,
    fraction: float = 0.05,
    seed: int = 20260802,
) -> list[Mapping[str, object]]:
    """Deterministic >=5% audit sample, stratified by declared diagnostics."""

    if not 0.0 < float(fraction) <= 1.0:
        raise ValueError("fraction must lie in (0, 1]")
    if not rows:
        return []
    strata_fields = (
        "epsilon_0",
        "beta_deg",
        "mu",
        "eta",
        "N_true",
        "thin_0p1_flag",
        "cutoff_regime",
    )
    ranked = sorted(
        rows,
        key=lambda row: hashlib.sha256(
            (str(seed) + "|" + "|".join(str(row.get(field, "")) for field in strata_fields)).encode("utf-8")
        ).hexdigest(),
    )
    required = max(1, int(math.ceil(float(fraction) * len(rows))))
    chosen: list[Mapping[str, object]] = []
    seen: set[tuple[object, ...]] = set()
    for row in ranked:
        stratum = tuple(row.get(field) for field in strata_fields)
        if stratum not in seen:
            chosen.append(row)
            seen.add(stratum)
        if len(chosen) >= required:
            break
    if len(chosen) < required:
        chosen_ids = {id(row) for row in chosen}
        chosen.extend(row for row in ranked if id(row) not in chosen_ids)  # pragma: no branch
    return chosen[:required]


def full_k10_control_reasons(
    row: Mapping[str, object],
    *,
    n_up_raw: int | None,
    sampled_case_ids: Sequence[str] = (),
) -> tuple[str, ...]:
    """Return the mandatory future full-K10 audit provenance for one row."""

    reasons: list[str] = []
    case_id = str(row.get("case_id", ""))
    if row.get("regression_label") in {"S3_12", "S3_14"}:
        reasons.append("exact_regression_point")
    n_true = row.get("N_true")
    try:
        n_value = int(n_true)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        n_value = None
    if n_up_raw is not None and n_value == int(n_up_raw):
        reasons.append("defines_N_up_raw")
    deltas = row.get("deltas", {})
    def near_threshold(value: object) -> bool:
        try:
            number = float(value)
        except (TypeError, ValueError):
            return False
        return math.isfinite(number) and 0.095 <= number <= 0.105

    if isinstance(deltas, Mapping) and any(near_threshold(value) for value in deltas.values()):
        reasons.append("near_delta_threshold")
    if row.get("execution_status") in {"attempted_unresolved", "attempted_error"}:
        reasons.append("unresolved_or_error_case")
    if bool(row.get("close_cluster")) or bool(row.get("root_ordering_warning")):
        reasons.append("cluster_or_ordering_case")
    if bool(row.get("selected_for_article")):
        reasons.append("selected_article_case")
    if case_id and case_id in set(str(value) for value in sampled_case_ids):
        reasons.append("stratified_five_percent_control")
    return tuple(dict.fromkeys(reasons))


__all__ = [
    "DEFAULT_STAGE_SAFE_ENDS",
    "EXECUTION_STATUSES",
    "IncrementalModelSession",
    "KNOWN_FULL_PATH_UNRESOLVED_SMOKE",
    "PREFIX_ALGORITHM_VERSION",
    "PREFIX_CACHE_SCHEMA_VERSION",
    "PREFIX_STRATEGIES",
    "PartialPointCache",
    "PrefixGuardAssessment",
    "STRICT_POLICIES",
    "StageTiming",
    "assess_prefix_guard",
    "choose_effective_strict_policy",
    "deterministic_control_sample",
    "envelope_saturation_status",
    "fresh_point_payload",
    "full_k10_control_reasons",
    "local_strict_modes",
    "is_known_full_path_unresolved_smoke",
    "prefix_cache_identity",
    "run_staged_point",
    "stage_guard_indices",
    "staged_settings",
    "strict_trigger_reasons",
    "validate_stage_safe_ends",
]
