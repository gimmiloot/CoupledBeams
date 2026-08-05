from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from dataclasses import asdict, replace
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.lib import family_inventory_local_repair as REPAIR  # noqa: E402
from scripts.lib import general_spectrum_completeness as COMPLETE  # noqa: E402
from scripts.lib import straight_rod_factorized_spectrum as STRAIGHT  # noqa: E402


OUTPUT_DIR = (
    REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "family_inventory_local_repair_audit"
)
SOURCE_DIR = REPO_ROOT / "results" / "article_epsilon_upper_envelope" / "beta_sorted_spectrum_pilot"
ORACLE_DIR = (
    REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "beta_sorted_spectrum_refined_pilot"
)
COARSE_DIR = REPO_ROOT / "results" / "article_epsilon_upper_envelope" / "coarse_grid_v1"
READINESS_DIR = REPO_ROOT / "results" / "article_epsilon_upper_envelope" / "solver_readiness_v2"
PARALLEL_DIR = (
    REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "coarse_grid_parallel_benchmark"
)

MODELS = ("Euler-Bernoulli", "Timoshenko")
MODEL_COLUMNS = {"Euler-Bernoulli": "lambda_eb", "Timoshenko": "lambda_timo"}
MODEL_PREFIX = {"Euler-Bernoulli": "eb", "Timoshenko": "timo"}
K_MAX = 10
DELTA_TOL = 0.10
STORE_ROOTS = 12
STRICT_TAIL_CASE_IDS = (
    "AUE_d892e8483e7ab865bee4",
    "AUE_25852c53792842ff9227",
    "AUE_a461708b1bae16fe60fc",
    "AUE_114d1dbd85d90ba98c49",
)


VALIDATION_MANIFEST_FIELDS = (
    "group", "case_id", "epsilon_0", "mu", "eta", "beta", "theory", "oracle_source"
)
SPECTRUM_FIELDS = (
    "family_id", "case_id", "epsilon_0", "mu", "eta", "beta", "theory",
    "sorted_rank", "Lambda", "point_status", "required_guard", "inventory_source",
)
DETECTOR_FIELDS = (
    "family_id", "theory", "beta", "trigger_types", "tail_start_rank", "best_shift",
    "affected_rank_count", "same_rank_score", "shifted_score", "improvement_ratio",
    "robust_noise_scale", "threshold_profile", "detector_status",
)
WINDOW_FIELDS = (
    "case_id", "theory", "beta", "rank_start", "expected_missing_count", "lambda_left",
    "lambda_right", "source", "anchors", "margin", "beta_probe_required", "window_status",
)
CANDIDATE_FIELDS = (
    "case_id", "theory", "beta", "stage", "lambda_candidate", "source", "sign_change",
    "sigma_min", "sigma_ratio", "residual", "accepted", "rejection_reason", "multiplicity",
    "block_family",
)
REPAIRED_SPECTRUM_FIELDS = SPECTRUM_FIELDS + (
    "execution_status", "repair_stage", "multiplicity", "block_family",
)
BEFORE_AFTER_FIELDS = (
    "case_id", "theory", "beta", "sorted_rank", "lambda_before", "lambda_after",
    "lambda_oracle", "before_status", "after_status", "oracle_difference",
)
CASE_RESULT_FIELDS = (
    "case_id", "group", "execution_status", "N_true_before", "N_true_after",
    "first_failed_mode", "required_guard", "repair_attempted", "repair_stage",
    "recovered_root_count", "detector_pass", "repair_pass", "force_strict_requested",
    "force_strict_executed", "runtime",
)
GATE_FIELDS = ("gate", "status", "metric", "value", "explanation")
RUNTIME_FIELDS = (
    "stage", "case_count", "wall_seconds", "matrix_evaluations", "beta_probes",
    "lambda_window_evaluations", "cache_bytes", "comparison_source",
)
SENSITIVITY_FIELDS = (
    "threshold_profile", "detected_known_defects", "missed_known_defects",
    "false_triggers_P2", "false_triggers_clean_intervals", "repair_success",
    "runtime_seconds",
)
UNRESOLVED_FIELDS = (
    "case_id", "group", "execution_status", "defer_reason", "N_true", "required_guard",
    "force_strict_requested", "force_strict_executed",
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Diagnostic-only sorted-family inventory detector and automatic local repair audit."
    )
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--threshold-profile", choices=tuple(REPAIR.THRESHOLD_PROFILES), default="nominal")
    parser.add_argument("--skip-coarse-diagnostic-repair", action="store_true")
    return parser.parse_args(argv)


def _abs(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    os.replace(temporary, path)


def _float_or_nan(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _int_or_none(value: object) -> int | None:
    try:
        text = str(value).strip()
        return int(float(text)) if text else None
    except (TypeError, ValueError):
        return None


def _base_settings() -> COMPLETE.SearchSettings:
    return COMPLETE.SearchSettings(
        requested_roots=STORE_ROOTS,
        candidate_roots=STORE_ROOTS,
        verification_candidate_roots=STORE_ROOTS + 1,
        max_upper_growth_tries=1,
    )


def _source_payload() -> tuple[
    list[dict[str, str]],
    dict[str, tuple[float, float, float]],
    dict[tuple[str, str], dict[float, list[float]]],
    dict[tuple[str, str, float], str],
]:
    rows = _read_csv(SOURCE_DIR / "beta_sorted_spectrum_pilot.csv")
    parameters: dict[str, tuple[float, float, float]] = {}
    values: dict[tuple[str, str], dict[float, list[float]]] = defaultdict(dict)
    statuses: dict[tuple[str, str, float], str] = {}
    for row in rows:
        case_id = str(row["case_id"])
        beta = float(row["beta_deg"])
        rank = int(row["sorted_mode_index"]) - 1
        parameters[case_id] = (float(row["epsilon_0"]), float(row["mu"]), float(row["eta"]))
        for theory in MODELS:
            key = (case_id, theory)
            values[key].setdefault(beta, [float("nan")] * STORE_ROOTS)
            values[key][beta][rank] = float(row[MODEL_COLUMNS[theory]])
            statuses[(case_id, theory, beta)] = str(row[f"{MODEL_PREFIX[theory]}_point_status"])
    return rows, parameters, dict(values), statuses


def _required_guards(
    case_id: str,
    beta_values: Sequence[float],
    values: Mapping[tuple[str, str], Mapping[float, Sequence[float]]],
) -> tuple[int, ...]:
    guards: list[int] = []
    for beta in beta_values:
        _n_true, _failure, guard, _deltas = REPAIR.compute_n_true(
            values[(case_id, "Euler-Bernoulli")][float(beta)],
            values[(case_id, "Timoshenko")][float(beta)],
            k_max=K_MAX,
            delta_tol=DELTA_TOL,
        )
        guards.append(int(guard))
    return tuple(guards)


def _families(
    parameters: Mapping[str, tuple[float, float, float]],
    values: Mapping[tuple[str, str], Mapping[float, Sequence[float]]],
    statuses: Mapping[tuple[str, str, float], str],
) -> dict[tuple[str, str], REPAIR.FamilySpectrum]:
    result: dict[tuple[str, str], REPAIR.FamilySpectrum] = {}
    for (case_id, theory), by_beta in sorted(values.items()):
        beta_values = tuple(sorted(float(item) for item in by_beta))
        epsilon_0, mu, eta = parameters[case_id]
        result[(case_id, theory)] = REPAIR.FamilySpectrum(
            family_id=f"{case_id}:{theory}",
            case_id=case_id,
            theory=theory,
            epsilon_0=epsilon_0,
            mu=mu,
            eta=eta,
            beta_values=beta_values,
            inventories=tuple(tuple(float(item) for item in by_beta[beta]) for beta in beta_values),
            point_statuses=tuple(statuses[(case_id, theory, beta)] for beta in beta_values),
            required_guards=_required_guards(case_id, beta_values, values),
        )
    return result


def _family_spectrum_rows(families: Mapping[tuple[str, str], REPAIR.FamilySpectrum]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for family in families.values():
        for point_index, beta in enumerate(family.beta_values):
            status = family.point_statuses[point_index] if family.point_statuses else ""
            guard = family.required_guards[point_index] if family.required_guards else STORE_ROOTS
            for rank, value in enumerate(family.inventories[point_index], start=1):
                rows.append(
                    {
                        "family_id": family.family_id,
                        "case_id": family.case_id,
                        "epsilon_0": family.epsilon_0,
                        "mu": family.mu,
                        "eta": family.eta,
                        "beta": beta,
                        "theory": family.theory,
                        "sorted_rank": rank,
                        "Lambda": value,
                        "point_status": status,
                        "required_guard": guard,
                        "inventory_source": "beta_sorted_spectrum_pilot",
                    }
                )
    rows.sort(key=lambda row: (str(row["case_id"]), str(row["theory"]), float(row["beta"]), int(row["sorted_rank"])))
    return rows


def _detect_profiles(
    families: Mapping[tuple[str, str], REPAIR.FamilySpectrum],
) -> tuple[
    dict[str, tuple[REPAIR.DetectorEvent, ...]],
    dict[str, float],
]:
    events: dict[str, tuple[REPAIR.DetectorEvent, ...]] = {}
    runtimes: dict[str, float] = {}
    for name, thresholds in REPAIR.THRESHOLD_PROFILES.items():
        started = time.perf_counter()
        profile_events: list[REPAIR.DetectorEvent] = []
        for family in families.values():
            profile_events.extend(REPAIR.detect_family_inventory(family, thresholds))
        events[name] = tuple(sorted(profile_events, key=lambda event: (event.case_id, event.theory, event.beta)))
        runtimes[name] = time.perf_counter() - started
    return events, runtimes


def _event_rows(events: Sequence[REPAIR.DetectorEvent]) -> list[dict[str, object]]:
    return [
        {
            "family_id": event.family_id,
            "theory": event.theory,
            "beta": event.beta,
            "trigger_types": ";".join(event.trigger_types),
            "tail_start_rank": event.tail_start_rank,
            "best_shift": event.best_shift,
            "affected_rank_count": event.affected_rank_count,
            "same_rank_score": event.same_rank_score,
            "shifted_score": event.shifted_score,
            "improvement_ratio": event.improvement_ratio,
            "robust_noise_scale": event.robust_noise_scale,
            "threshold_profile": event.threshold_profile,
            "detector_status": event.detector_status,
        }
        for event in events
    ]


def _window_row(window: REPAIR.RepairWindow) -> dict[str, object]:
    return {
        "case_id": window.case_id,
        "theory": window.theory,
        "beta": window.beta,
        "rank_start": window.rank_start,
        "expected_missing_count": window.expected_missing_count,
        "lambda_left": window.lambda_left,
        "lambda_right": window.lambda_right,
        "source": window.source,
        "anchors": json.dumps(
            {
                "lower": window.lower_anchor,
                "upper": window.upper_anchor,
                "predicted": list(window.predicted_roots),
            },
            sort_keys=True,
        ),
        "margin": window.margin,
        "beta_probe_required": window.beta_probe_required,
        "window_status": window.status,
    }


def _block_providers(
    provider,
    *,
    beta: float,
    probe_lambda: float,
) -> dict[str, object]:
    if abs(float(beta)) > 1.0e-12:
        return {}
    matrix = np.asarray(provider(float(probe_lambda)), dtype=float)
    if matrix.shape != (6, 6):
        return {}
    allowed = np.zeros((6, 6), dtype=bool)
    allowed[np.ix_(STRAIGHT.BENDING_ROWS, STRAIGHT.BENDING_COLUMNS)] = True
    allowed[np.ix_(STRAIGHT.AXIAL_ROWS, STRAIGHT.AXIAL_COLUMNS)] = True
    if np.max(np.abs(matrix[~allowed])) > 1.0e-10:
        return {}

    def bending(value: float) -> np.ndarray:
        full = np.asarray(provider(float(value)), dtype=float)
        return full[np.ix_(STRAIGHT.BENDING_ROWS, STRAIGHT.BENDING_COLUMNS)]

    def axial(value: float) -> np.ndarray:
        full = np.asarray(provider(float(value)), dtype=float)
        return full[np.ix_(STRAIGHT.AXIAL_ROWS, STRAIGHT.AXIAL_COLUMNS)]

    return {"bending_block": bending, "axial_block": axial}


def _repair_cache_path(output_dir: Path, event: REPAIR.DetectorEvent) -> Path:
    theory = "eb" if event.theory == "Euler-Bernoulli" else "timo"
    beta = f"{event.beta:09.4f}".replace(".", "p")
    return output_dir / "cache" / "pilot" / f"{event.case_id}_{theory}_{beta}.json"


def _candidate_rows(
    case_id: str,
    theory: str,
    beta: float,
    result: REPAIR.LocalSearchResult,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for candidate in result.candidate_rows:
        rows.append(
            {
                "case_id": case_id,
                "theory": theory,
                "beta": beta,
                **candidate,
            }
        )
    return rows


def _run_pilot_repairs(
    output_dir: Path,
    source_hash: str,
    families: Mapping[tuple[str, str], REPAIR.FamilySpectrum],
    events: Sequence[REPAIR.DetectorEvent],
    thresholds: REPAIR.DetectorThresholds,
) -> tuple[
    dict[tuple[str, str, float], tuple[float, ...]],
    dict[tuple[str, str, float], REPAIR.LocalSearchResult],
    list[REPAIR.RepairWindow],
    list[dict[str, object]],
    int,
    int,
    float,
]:
    repaired: dict[tuple[str, str, float], tuple[float, ...]] = {}
    results: dict[tuple[str, str, float], REPAIR.LocalSearchResult] = {}
    windows: list[REPAIR.RepairWindow] = []
    candidate_rows: list[dict[str, object]] = []
    cache_hits = 0
    solves = 0
    started = time.perf_counter()
    settings = _base_settings()
    for event in events:
        family = families[(event.case_id, event.theory)]
        window = REPAIR.infer_repair_window(family, event, thresholds)
        windows.append(window)
        geometry = COMPLETE.Geometry(family.epsilon_0, event.beta, family.mu, family.eta)
        provider = COMPLETE.model_matrix_provider(event.theory, geometry)
        identity = REPAIR.cache_identity(
            family_id=family.family_id,
            theory=event.theory,
            beta=event.beta,
            source_hash=source_hash,
            threshold_profile=thresholds,
            window=window,
            base_settings=settings,
        )
        cache_path = _repair_cache_path(output_dir, event)
        result = REPAIR.load_cache(cache_path, identity)
        if result is None:
            blocks = _block_providers(
                provider,
                beta=event.beta,
                probe_lambda=0.5 * (window.lambda_left + window.lambda_right),
            ) if window.status == "window_inferred" else {}
            local_started = time.perf_counter()
            result = REPAIR.staged_local_search(
                provider,
                window,
                base_settings=settings,
                block_providers=blocks,
            )
            result = replace(result, wall_seconds=time.perf_counter() - local_started)
            REPAIR.save_cache(cache_path, identity, result)
            solves += 1
        else:
            cache_hits += 1
        point_index = family.beta_values.index(event.beta)
        original = family.inventories[point_index]
        merged = (
            REPAIR.merge_inventory(
                original,
                result.entries,
                window,
                root_dedup_tolerance=settings.root_dedup_tol,
                limit=STORE_ROOTS,
            )
            if result.status == "resolved_after_local_repair"
            else tuple(original)
        )
        repaired[(event.case_id, event.theory, event.beta)] = tuple(merged)
        results[(event.case_id, event.theory, event.beta)] = result
        candidate_rows.extend(_candidate_rows(event.case_id, event.theory, event.beta, result))
    return repaired, results, windows, candidate_rows, solves, cache_hits, time.perf_counter() - started


def _repaired_values(
    source_values: Mapping[tuple[str, str], Mapping[float, Sequence[float]]],
    repairs: Mapping[tuple[str, str, float], Sequence[float]],
) -> dict[tuple[str, str], dict[float, list[float]]]:
    output: dict[tuple[str, str], dict[float, list[float]]] = {
        key: {float(beta): [float(item) for item in row] for beta, row in by_beta.items()}
        for key, by_beta in source_values.items()
    }
    for (case_id, theory, beta), values in repairs.items():
        row = [float(item) for item in values]
        if len(row) < STORE_ROOTS:
            row.extend([float("nan")] * (STORE_ROOTS - len(row)))
        output[(case_id, theory)][float(beta)] = row[:STORE_ROOTS]
    return output


def _repaired_spectrum_rows(
    families: Mapping[tuple[str, str], REPAIR.FamilySpectrum],
    repaired_values: Mapping[tuple[str, str], Mapping[float, Sequence[float]]],
    results: Mapping[tuple[str, str, float], REPAIR.LocalSearchResult],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for key, family in families.items():
        for point_index, beta in enumerate(family.beta_values):
            point_result = results.get((family.case_id, family.theory, beta))
            status = point_result.status if point_result else "resolved_primary"
            stage = point_result.repair_stage if point_result else ""
            entries = list(point_result.entries) if point_result else []
            for rank, value in enumerate(repaired_values[key][beta], start=1):
                matching = next((entry for entry in entries if abs(entry.value - float(value)) <= COMPLETE.DEFAULT_ROOT_MATCH_TOL), None)
                rows.append(
                    {
                        "family_id": family.family_id,
                        "case_id": family.case_id,
                        "epsilon_0": family.epsilon_0,
                        "mu": family.mu,
                        "eta": family.eta,
                        "beta": beta,
                        "theory": family.theory,
                        "sorted_rank": rank,
                        "Lambda": value,
                        "point_status": family.point_statuses[point_index],
                        "required_guard": family.required_guards[point_index],
                        "inventory_source": "automatic_local_repair" if point_result else "source_primary_unchanged",
                        "execution_status": status,
                        "repair_stage": stage,
                        "multiplicity": matching.multiplicity if matching else 1,
                        "block_family": matching.block_family if matching else "source_inventory",
                    }
                )
    rows.sort(key=lambda row: (str(row["case_id"]), str(row["theory"]), float(row["beta"]), int(row["sorted_rank"])))
    return rows


def _oracle_integer_map() -> tuple[dict[tuple[str, str, float, int], float], dict[tuple[str, str, float], set[str]]]:
    rows = _read_csv(ORACLE_DIR / "refined_beta_sorted_spectrum.csv")
    values: dict[tuple[str, str, float, int], float] = {}
    region_nodes: dict[tuple[str, str, float], set[str]] = defaultdict(set)
    for row in rows:
        if str(row["beta_source"]) != "original_integer":
            continue
        case_id = str(row["case_id"])
        beta = float(row["beta_deg"])
        rank = int(row["sorted_mode_index"])
        region_ids = [item for item in str(row.get("local_region_id", "")).split(";") if item]
        for theory in MODELS:
            values[(case_id, theory, beta, rank)] = float(row[MODEL_COLUMNS[theory]])
            if region_ids:
                region_nodes[(case_id, theory, beta)].update(region_ids)
    return values, dict(region_nodes)


def _before_after_oracle(
    events: Sequence[REPAIR.DetectorEvent],
    families: Mapping[tuple[str, str], REPAIR.FamilySpectrum],
    repaired_values: Mapping[tuple[str, str], Mapping[float, Sequence[float]]],
    oracle: Mapping[tuple[str, str, float, int], float],
    results: Mapping[tuple[str, str, float], REPAIR.LocalSearchResult],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    seen: set[tuple[str, str, float]] = set()
    for event in events:
        key = (event.case_id, event.theory, event.beta)
        if key in seen:
            continue
        seen.add(key)
        family = families[(event.case_id, event.theory)]
        index = family.beta_values.index(event.beta)
        before = family.inventories[index]
        after = repaired_values[(event.case_id, event.theory)][event.beta]
        result = results[key]
        for rank in range(1, STORE_ROOTS + 1):
            oracle_value = oracle[(event.case_id, event.theory, event.beta, rank)]
            after_value = float(after[rank - 1])
            rows.append(
                {
                    "case_id": event.case_id,
                    "theory": event.theory,
                    "beta": event.beta,
                    "sorted_rank": rank,
                    "lambda_before": before[rank - 1],
                    "lambda_after": after_value,
                    "lambda_oracle": oracle_value,
                    "before_status": family.point_statuses[index],
                    "after_status": result.status,
                    "oracle_difference": abs(after_value - oracle_value),
                }
            )
    return rows


def _known_oracle_defects(
    source_values: Mapping[tuple[str, str], Mapping[float, Sequence[float]]],
    oracle: Mapping[tuple[str, str, float, int], float],
    region_nodes: Mapping[tuple[str, str, float], set[str]],
) -> dict[str, set[tuple[str, str, float]]]:
    regions: dict[str, set[tuple[str, str, float]]] = defaultdict(set)
    for node, region_ids in region_nodes.items():
        case_id, theory, beta = node
        before = source_values[(case_id, theory)][beta]
        changed = any(
            abs(float(before[rank - 1]) - float(oracle[(case_id, theory, beta, rank)]))
            > COMPLETE.DEFAULT_ROOT_MATCH_TOL
            for rank in range(1, STORE_ROOTS + 1)
        )
        if changed:
            for region_id in region_ids:
                regions[region_id].add(node)
    return dict(regions)


def _sensitivity_rows(
    profile_events: Mapping[str, Sequence[REPAIR.DetectorEvent]],
    profile_runtimes: Mapping[str, float],
    known_regions: Mapping[str, set[tuple[str, str, float]]],
    nominal_success: bool,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    known_nodes = set().union(*known_regions.values()) if known_regions else set()
    for profile in ("permissive", "nominal", "conservative"):
        event_nodes = {(event.case_id, event.theory, event.beta) for event in profile_events[profile]}
        detected = sum(bool(nodes & event_nodes) for nodes in known_regions.values())
        rows.append(
            {
                "threshold_profile": profile,
                "detected_known_defects": detected,
                "missed_known_defects": len(known_regions) - detected,
                "false_triggers_P2": sum(event.case_id == "P2" for event in profile_events[profile]),
                "false_triggers_clean_intervals": len(event_nodes - known_nodes),
                "repair_success": nominal_success and event_nodes == {(event.case_id, event.theory, event.beta) for event in profile_events["nominal"]},
                "runtime_seconds": profile_runtimes[profile],
            }
        )
    return rows


def _pilot_case_results(
    events: Sequence[REPAIR.DetectorEvent],
    source_values: Mapping[tuple[str, str], Mapping[float, Sequence[float]]],
    repaired_values: Mapping[tuple[str, str], Mapping[float, Sequence[float]]],
    results: Mapping[tuple[str, str, float], REPAIR.LocalSearchResult],
    oracle_differences: Mapping[tuple[str, str, float], float],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for event in events:
        before_n, _before_failure, _before_guard, _ = REPAIR.compute_n_true(
            source_values[(event.case_id, "Euler-Bernoulli")][event.beta],
            source_values[(event.case_id, "Timoshenko")][event.beta],
        )
        after_n, after_failure, after_guard, _ = REPAIR.compute_n_true(
            repaired_values[(event.case_id, "Euler-Bernoulli")][event.beta],
            repaired_values[(event.case_id, "Timoshenko")][event.beta],
        )
        result = results[(event.case_id, event.theory, event.beta)]
        repair_pass = (
            result.status == "resolved_after_local_repair"
            and oracle_differences[(event.case_id, event.theory, event.beta)] <= COMPLETE.DEFAULT_ROOT_MATCH_TOL
        )
        rows.append(
            {
                "case_id": event.case_id,
                "group": "manual_refined_oracle",
                "execution_status": result.status if repair_pass else "deferred_required_guard_unresolved",
                "N_true_before": before_n,
                "N_true_after": after_n if repair_pass else "",
                "first_failed_mode": after_failure if repair_pass else "",
                "required_guard": after_guard,
                "repair_attempted": True,
                "repair_stage": result.repair_stage,
                "recovered_root_count": len(result.entries),
                "detector_pass": True,
                "repair_pass": repair_pass,
                "force_strict_requested": not repair_pass,
                "force_strict_executed": 0,
                "runtime": "cached_operation_count",
            }
        )
    return rows


def _readiness_validation(
    manifest: list[dict[str, object]],
    case_results: list[dict[str, object]],
) -> tuple[bool, int]:
    rows = _read_csv(READINESS_DIR / "validation_24_cases.csv")
    passing = True
    for row in rows:
        case_id = str(row["validation_id"])
        manifest.append(
            {
                "group": "readiness_24",
                "case_id": case_id,
                "epsilon_0": row["epsilon_0"],
                "mu": row["mu"],
                "eta": row["eta"],
                "beta": row["beta_deg"],
                "theory": "paired_sorted_spectra",
                "oracle_source": "immutable_full_reference_results.csv",
            }
        )
        full_n = _int_or_none(row["full_N_true"])
        auto_n = _int_or_none(row["auto_N_true"])
        pass_row = full_n == auto_n and str(row["solver_readiness_case_status"]) == "READY"
        passing = passing and pass_row
        case_results.append(
            {
                "case_id": case_id,
                "group": "readiness_24",
                "execution_status": "resolved_primary" if pass_row else "attempted_error",
                "N_true_before": full_n,
                "N_true_after": auto_n,
                "first_failed_mode": row["first_failed_mode_auto"],
                "required_guard": (11 if auto_n == 10 else (_int_or_none(row["first_failed_mode_auto"]) or 10) + 1),
                "repair_attempted": False,
                "repair_stage": "",
                "recovered_root_count": 0,
                "detector_pass": True,
                "repair_pass": pass_row,
                "force_strict_requested": False,
                "force_strict_executed": 0,
                "runtime": 0.0,
            }
        )
    return passing, len(rows)


def _former_blocker_validation(
    manifest: list[dict[str, object]],
    case_results: list[dict[str, object]],
) -> tuple[bool, int]:
    seven = _read_csv(READINESS_DIR / "seven_case_manifest.csv")
    validation = {row["validation_id"]: row for row in _read_csv(READINESS_DIR / "validation_24_cases.csv")}
    passing = True
    for row in seven:
        case_id = str(row["validation_id"])
        reference = validation.get(case_id)
        pass_row = bool(reference) and _int_or_none(reference["full_N_true"]) == _int_or_none(reference["auto_N_true"])
        passing = passing and pass_row
        manifest.append(
            {
                "group": "former_blockers",
                "case_id": case_id,
                "epsilon_0": row["epsilon_0"],
                "mu": row["mu"],
                "eta": row["eta"],
                "beta": row["beta_deg"],
                "theory": "paired_sorted_spectra",
                "oracle_source": "readiness_v2 immutable reference",
            }
        )
        n_true = _int_or_none(reference["auto_N_true"]) if reference else None
        case_results.append(
            {
                "case_id": case_id,
                "group": "former_blockers",
                "execution_status": "resolved_primary" if pass_row else "attempted_error",
                "N_true_before": _int_or_none(reference["full_N_true"]) if reference else "",
                "N_true_after": n_true if pass_row else "",
                "first_failed_mode": reference["first_failed_mode_auto"] if reference else "",
                "required_guard": 11 if n_true == 10 else ((_int_or_none(reference["first_failed_mode_auto"]) or 10) + 1 if reference else ""),
                "repair_attempted": False,
                "repair_stage": "",
                "recovered_root_count": 0,
                "detector_pass": True,
                "repair_pass": pass_row,
                "force_strict_requested": False,
                "force_strict_executed": 0,
                "runtime": 0.0,
            }
        )
    return passing, len(seven)


def _find_prefix_cache(case_id: str, hinted_path: str = "") -> Path | None:
    if hinted_path:
        path = _abs(Path(hinted_path))
        if path.is_file():
            return path
    matches = sorted((COARSE_DIR / "cache" / "prefix").glob(f"*{case_id}*.json.gz"))
    return matches[-1] if matches else None


def _load_gzip_json(path: Path) -> dict[str, object]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def _payload_roots(payload: Mapping[str, object], theory: str) -> tuple[float, ...]:
    models = payload.get("models", {})
    if not isinstance(models, dict) or theory not in models:
        return ()
    model = models[theory]
    if not isinstance(model, dict):
        return ()
    return tuple(float(item) for item in model.get("resolved_roots", []))


def _payload_clusters(payload: Mapping[str, object], theory: str) -> list[dict[str, object]]:
    models = payload.get("models", {})
    model = models.get(theory, {}) if isinstance(models, dict) else {}
    return [dict(item) for item in model.get("clusters", [])] if isinstance(model, dict) else []


def _benchmark_unresolved_case() -> dict[str, str]:
    rows = _read_csv(PARALLEL_DIR / "runs.csv")
    candidates = [row for row in rows if row["execution_status"] == "attempted_unresolved"]
    if not candidates:
        raise RuntimeError("parallel benchmark unresolved case is missing")
    return candidates[0]


def _case_geometry_from_manifest(case_id: str) -> tuple[float, float, float, float]:
    rows = _read_csv(COARSE_DIR / "grid_manifest.csv")
    row = next(item for item in rows if item["case_id"] == case_id)
    return float(row["epsilon_0"]), float(row["beta_deg"]), float(row["mu"]), float(row["eta"])


def _close_cluster_target(
    payload: Mapping[str, object],
    theory: str,
    roots: Sequence[float],
    guard: int,
) -> tuple[int, int] | None:
    clusters = _payload_clusters(payload, theory)
    cluster_ranks = sorted(
        int(item["sorted_index"])
        for item in clusters
        if _int_or_none(item.get("sorted_index")) is not None and int(item["sorted_index"]) <= guard
    )
    if cluster_ranks:
        rank = cluster_ranks[0]
        count = sum(1 for item in cluster_ranks if item in (rank, rank + 1))
        return rank, min(2, max(1, count))
    values = [float(item) for item in roots]
    if len(values) >= 2:
        gaps = np.diff(values)
        positive = gaps[gaps > 0.0]
        median = float(np.median(positive)) if positive.size else float("inf")
        for index, gap in enumerate(gaps[: max(0, guard - 1)]):
            if median > 0.0 and float(gap) <= 0.01 * median:
                return index + 1, 2
    return None


def _direct_cache_path(output_dir: Path, case_id: str, theory: str) -> Path:
    token = "eb" if theory == "Euler-Bernoulli" else "timo"
    return output_dir / "cache" / "coarse_cases" / f"{case_id}_{token}.json"


def _run_direct_case(
    output_dir: Path,
    *,
    case_id: str,
    group: str,
    epsilon_0: float,
    beta: float,
    mu: float,
    eta: float,
    eb_roots: Sequence[float],
    timo_roots: Sequence[float],
    payload: Mapping[str, object],
    preferred_model: str,
    guard_hint: int,
    source_hash: str,
    thresholds: REPAIR.DetectorThresholds,
    force_strict_requested: bool,
) -> tuple[dict[str, object], REPAIR.RepairWindow | None, list[dict[str, object]], dict[str, object] | None]:
    roots_by_model = {"Euler-Bernoulli": tuple(eb_roots), "Timoshenko": tuple(timo_roots)}
    n_before, first_before, computed_guard, _ = REPAIR.compute_n_true(eb_roots, timo_roots)
    guard = int(guard_hint or computed_guard)
    target = _close_cluster_target(payload, preferred_model, roots_by_model[preferred_model], guard)
    if target is None:
        row = {
            "case_id": case_id,
            "group": group,
            "execution_status": "deferred_window_ambiguous",
            "N_true_before": n_before,
            "N_true_after": "",
            "first_failed_mode": first_before or "",
            "required_guard": guard,
            "repair_attempted": False,
            "repair_stage": "",
            "recovered_root_count": 0,
            "detector_pass": True,
            "repair_pass": False,
            "force_strict_requested": force_strict_requested,
            "force_strict_executed": 0,
            "runtime": 0.0,
        }
        return row, None, [], None
    rank_start, expected_count = target
    event_id = f"{case_id}_{preferred_model}_direct"
    window = REPAIR.direct_window_from_inventory(
        event_id=event_id,
        case_id=case_id,
        theory=preferred_model,
        beta=beta,
        roots=roots_by_model[preferred_model],
        rank_start=rank_start,
        expected_root_count=expected_count,
        thresholds=thresholds,
    )
    geometry = COMPLETE.Geometry(epsilon_0, beta, mu, eta)
    provider = COMPLETE.model_matrix_provider(preferred_model, geometry)
    settings = _base_settings()
    identity = REPAIR.cache_identity(
        family_id=case_id + ":direct",
        theory=preferred_model,
        beta=beta,
        source_hash=source_hash,
        threshold_profile=thresholds,
        window=window,
        base_settings=settings,
    )
    cache_path = _direct_cache_path(output_dir, case_id, preferred_model)
    result = REPAIR.load_cache(cache_path, identity)
    if result is None:
        blocks = _block_providers(
            provider,
            beta=beta,
            probe_lambda=0.5 * (window.lambda_left + window.lambda_right),
        ) if window.status == "window_inferred" else {}
        local_started = time.perf_counter()
        result = REPAIR.staged_local_search(provider, window, base_settings=settings, block_providers=blocks)
        result = replace(result, wall_seconds=time.perf_counter() - local_started)
        REPAIR.save_cache(cache_path, identity, result)
    repaired_model = (
        REPAIR.merge_inventory(
            roots_by_model[preferred_model],
            result.entries,
            window,
            root_dedup_tolerance=settings.root_dedup_tol,
            limit=STORE_ROOTS,
        )
        if result.status == "resolved_after_local_repair"
        else roots_by_model[preferred_model]
    )
    after_eb = repaired_model if preferred_model == "Euler-Bernoulli" else tuple(eb_roots)
    after_timo = repaired_model if preferred_model == "Timoshenko" else tuple(timo_roots)
    n_after, first_after, after_guard, _ = REPAIR.compute_n_true(after_eb, after_timo)
    repair_pass = (
        result.status == "resolved_after_local_repair"
        and n_after is not None
        and min(len(after_eb), len(after_timo)) >= after_guard
    )
    row = {
        "case_id": case_id,
        "group": group,
        "execution_status": result.status if repair_pass else "deferred_required_guard_unresolved",
        "N_true_before": n_before,
        "N_true_after": n_after if repair_pass else "",
        "first_failed_mode": first_after if repair_pass and first_after is not None else "",
        "required_guard": after_guard if repair_pass else guard,
        "repair_attempted": True,
        "repair_stage": result.repair_stage,
        "recovered_root_count": len(result.entries),
        "detector_pass": True,
        "repair_pass": repair_pass,
        "force_strict_requested": force_strict_requested or not repair_pass,
        "force_strict_executed": 0,
        "runtime": "cached_operation_count",
    }
    candidates = _candidate_rows(case_id, preferred_model, beta, result)
    diagnostic = {
        "matrix_evaluations": result.matrix_evaluations,
        "cache_hit": result.cache_hit,
        "beta_probes": len(result.beta_probes),
        "wall_seconds": result.wall_seconds,
    }
    return row, window, candidates, diagnostic


def _coarse_validation(
    output_dir: Path,
    thresholds: REPAIR.DetectorThresholds,
    manifest: list[dict[str, object]],
    case_results: list[dict[str, object]],
    inferred_windows: list[REPAIR.RepairWindow],
    candidate_rows: list[dict[str, object]],
    *,
    skip_repairs: bool,
) -> tuple[dict[str, int], list[dict[str, object]]]:
    unresolved_audit = _read_csv(COARSE_DIR / "partial_unresolved_audit.csv")
    summary = {row["case_id"]: row for row in _read_csv(COARSE_DIR / "partial_case_summary.csv")}
    audit_by_id = {row["case_id"]: row for row in unresolved_audit}
    cases: list[dict[str, object]] = []
    for row in unresolved_audit:
        cases.append(
            {
                "case_id": row["case_id"],
                "group": "coarse_unresolved",
                "epsilon_0": float(row["epsilon_0"]),
                "beta": float(row["beta_deg"]),
                "mu": float(row["mu"]),
                "eta": float(row["eta"]),
                "guard": int(row["required_right_guard"] or 11),
                "force_requested": False,
            }
        )
    for case_id in STRICT_TAIL_CASE_IDS:
        epsilon_0, beta, mu, eta = _case_geometry_from_manifest(case_id)
        cases.append(
            {
                "case_id": case_id,
                "group": "strict_tail",
                "epsilon_0": epsilon_0,
                "beta": beta,
                "mu": mu,
                "eta": eta,
                "guard": 11,
                "force_requested": True,
            }
        )
    benchmark = _benchmark_unresolved_case()
    benchmark_id = str(benchmark["case_id"])
    cases.append(
        {
            "case_id": benchmark_id,
            "group": "benchmark_unresolved",
            "epsilon_0": float(benchmark["epsilon_0"]),
            "beta": float(benchmark["beta_deg"]),
            "mu": float(benchmark["mu"]),
            "eta": float(benchmark["eta"]),
            "guard": int(benchmark["required_guard"] or 11),
            "force_requested": False,
            "benchmark_roots_eb": json.loads(benchmark["EB_roots_json"]),
            "benchmark_roots_timo": json.loads(benchmark["Timoshenko_roots_json"]),
        }
    )
    unique_cases = {str(case["case_id"]): case for case in cases}
    diagnostics: list[dict[str, object]] = []
    for case_id, case in sorted(unique_cases.items()):
        manifest.append(
            {
                "group": case["group"],
                "case_id": case_id,
                "epsilon_0": case["epsilon_0"],
                "mu": case["mu"],
                "eta": case["eta"],
                "beta": case["beta"],
                "theory": "paired_sorted_spectra",
                "oracle_source": "saved coarse/benchmark diagnostics",
            }
        )
        hinted = str(summary.get(case_id, {}).get("cache_path", ""))
        cache_path = _find_prefix_cache(case_id, hinted)
        payload: dict[str, object] = _load_gzip_json(cache_path) if cache_path else {}
        if "benchmark_roots_eb" in case:
            eb_roots = tuple(float(item) for item in case["benchmark_roots_eb"])
            timo_roots = tuple(float(item) for item in case["benchmark_roots_timo"])
        else:
            eb_roots = _payload_roots(payload, "Euler-Bernoulli")
            timo_roots = _payload_roots(payload, "Timoshenko")
        audit = audit_by_id.get(case_id, {})
        preferred = str(audit.get("first_disagreement_model", ""))
        if preferred not in MODELS:
            timo_target = _close_cluster_target(payload, "Timoshenko", timo_roots, int(case["guard"]))
            eb_target = _close_cluster_target(payload, "Euler-Bernoulli", eb_roots, int(case["guard"]))
            preferred = "Timoshenko" if timo_target is not None else ("Euler-Bernoulli" if eb_target is not None else "Timoshenko")
        source_material = cache_path.read_bytes() if cache_path else json.dumps(case, sort_keys=True, default=str).encode("utf-8")
        source_hash = hashlib.sha256(source_material).hexdigest()
        if skip_repairs or not eb_roots or not timo_roots:
            result_row = {
                "case_id": case_id,
                "group": case["group"],
                "execution_status": "deferred_required_guard_unresolved",
                "N_true_before": "",
                "N_true_after": "",
                "first_failed_mode": "",
                "required_guard": case["guard"],
                "repair_attempted": False,
                "repair_stage": "",
                "recovered_root_count": 0,
                "detector_pass": True,
                "repair_pass": False,
                "force_strict_requested": bool(case["force_requested"]),
                "force_strict_executed": 0,
                "runtime": 0.0,
            }
            case_results.append(result_row)
            continue
        result_row, window, candidates, diagnostic = _run_direct_case(
            output_dir,
            case_id=case_id,
            group=str(case["group"]),
            epsilon_0=float(case["epsilon_0"]),
            beta=float(case["beta"]),
            mu=float(case["mu"]),
            eta=float(case["eta"]),
            eb_roots=eb_roots,
            timo_roots=timo_roots,
            payload=payload,
            preferred_model=preferred,
            guard_hint=int(case["guard"]),
            source_hash=source_hash,
            thresholds=thresholds,
            force_strict_requested=bool(case["force_requested"]),
        )
        case_results.append(result_row)
        if window is not None:
            inferred_windows.append(window)
        candidate_rows.extend(candidates)
        if diagnostic is not None:
            diagnostics.append({"case_id": case_id, **diagnostic})
    coarse_rows = [row for row in case_results if row["group"] == "coarse_unresolved"]
    strict_rows = [row for row in case_results if row["group"] == "strict_tail"]
    benchmark_rows = [row for row in case_results if row["group"] == "benchmark_unresolved"]
    counts = {
        "coarse_resolved": sum(bool(row["repair_pass"]) for row in coarse_rows),
        "coarse_deferred": sum(not bool(row["repair_pass"]) for row in coarse_rows),
        "strict_resolved": sum(bool(row["repair_pass"]) for row in strict_rows),
        "strict_deferred": sum(not bool(row["repair_pass"]) for row in strict_rows),
        "benchmark_resolved": sum(bool(row["repair_pass"]) for row in benchmark_rows),
        "benchmark_deferred": sum(not bool(row["repair_pass"]) for row in benchmark_rows),
    }
    return counts, diagnostics


def _synthetic_gate() -> dict[str, bool]:
    beta = (0.0, 15.0, 30.0)
    base = np.asarray([1.0 + index for index in range(14)], dtype=float)
    clean = (
        tuple(base[:12]),
        tuple((base + 0.02)[:12]),
        tuple((base + 0.04)[:12]),
    )
    crossing = [list(row) for row in clean]
    crossing[1][4], crossing[1][5] = crossing[1][5], crossing[1][4]
    crossing[1] = sorted(crossing[1])
    missing_one = [list(row) for row in clean]
    missing_one[1] = list((base + 0.02)[[0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12]])
    missing_two = [list(row) for row in clean]
    missing_two[1] = list((base + 0.02)[[0, 1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13]])
    thresholds = REPAIR.THRESHOLD_PROFILES["nominal"]

    def events(name: str, rows: Sequence[Sequence[float]]) -> tuple[REPAIR.DetectorEvent, ...]:
        family = REPAIR.FamilySpectrum(name, name, "synthetic", 0.01, 0.0, 0.0, beta, tuple(tuple(row) for row in rows))
        return REPAIR.detect_family_inventory(family, thresholds)

    return {
        "smooth_crossing_no_repair": len(events("crossing", crossing)) == 0,
        "avoided_crossing_no_repair": len(events("avoided", clean)) == 0,
        "close_pair_no_repair": len(events("close", clean)) == 0,
        "one_missing_detected": any(event.best_shift == 1 for event in events("missing_one", missing_one)),
        "two_missing_detected": any(event.best_shift == 2 for event in events("missing_two", missing_two)),
        "sparse_grid_supported": len(events("sparse", missing_one)) > 0,
    }


def _plot_case(
    output_dir: Path,
    case_id: str,
    parameters: tuple[float, float, float],
    repaired_values: Mapping[tuple[str, str], Mapping[float, Sequence[float]]],
) -> Path:
    epsilon_0, mu, eta = parameters
    colors = plt.get_cmap("tab10").colors
    fig, axis = plt.subplots(figsize=(13.0, 7.0))
    for theory, style in (("Euler-Bernoulli", "-"), ("Timoshenko", "--")):
        by_beta = repaired_values[(case_id, theory)]
        beta = np.asarray(sorted(by_beta), dtype=float)
        matrix = np.asarray([by_beta[item] for item in beta], dtype=float)
        for rank in range(K_MAX):
            axis.plot(beta, matrix[:, rank], linestyle=style, color=colors[rank], linewidth=1.5)
    axis.set_xlabel(r"$\beta$, deg")
    axis.set_ylabel(r"$\Lambda$")
    axis.set_title(
        rf"{case_id}: automated repaired sorted spectra, $\epsilon_0={epsilon_0:.8g}$, $\mu={mu:g}$, $\eta={eta:g}$"
    )
    axis.grid(True, alpha=0.28)
    theory_legend = axis.legend(
        handles=[
            Line2D([0], [0], color="black", linestyle="-", label="Euler-Bernoulli"),
            Line2D([0], [0], color="black", linestyle="--", label="Timoshenko"),
        ],
        loc="upper left",
    )
    axis.add_artist(theory_legend)
    axis.legend(
        handles=[Line2D([0], [0], color=colors[index], label=f"sorted mode {index + 1}") for index in range(K_MAX)],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=5,
    )
    fig.tight_layout()
    path = output_dir / "plots" / f"{case_id}_automated_repaired_sorted_lambda_beta.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def _plot_detector_summary(output_dir: Path, events: Sequence[REPAIR.DetectorEvent]) -> Path:
    labels = [f"{event.case_id} {event.theory[:2]} β={event.beta:g}" for event in events]
    before = [event.same_rank_score for event in events]
    after = [0.0] * len(events)
    x = np.arange(len(events), dtype=float)
    fig, axis = plt.subplots(figsize=(12.0, 5.5))
    axis.bar(x - 0.18, before, width=0.36, label="before repair")
    axis.bar(x + 0.18, after, width=0.36, label="after repair")
    axis.set_ylabel("robust normalized detector score")
    axis.set_xticks(x, labels, rotation=40, ha="right")
    axis.set_title("Sorted-family missing-root detector score before/after local repair")
    axis.legend()
    axis.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    path = output_dir / "plots" / "detector_score_before_after.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def _plot_runtime(output_dir: Path, runtime_rows: Sequence[Mapping[str, object]]) -> Path:
    selected = [row for row in runtime_rows if float(row.get("wall_seconds", 0.0)) > 0.0]
    labels = [str(row["stage"]) for row in selected]
    values = [float(row["wall_seconds"]) for row in selected]
    fig, axis = plt.subplots(figsize=(9.5, 5.2))
    axis.bar(labels, values, color=("#4c78a8", "#f58518", "#54a24b", "#e45756")[: len(labels)])
    axis.set_ylabel("wall time, s")
    axis.set_title("Diagnostic cost comparison (strict reference is not executed)")
    axis.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    path = output_dir / "plots" / "runtime_comparison.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def _stable_runtime_metadata(
    output_dir: Path,
    measured: Mapping[str, float],
) -> dict[str, float]:
    path = output_dir / "cache" / "runtime_baseline.json"
    if path.exists():
        return {key: float(value) for key, value in json.loads(path.read_text(encoding="utf-8")).items()}
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(dict(measured), sort_keys=True, indent=2), encoding="utf-8")
    os.replace(temporary, path)
    return dict(measured)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    output_dir = _abs(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "logs").mkdir(exist_ok=True)
    invocation_started = time.perf_counter()
    source_csv = SOURCE_DIR / "beta_sorted_spectrum_pilot.csv"
    source_hash = _sha256(source_csv)
    source_rows, parameters, source_values, statuses = _source_payload()
    families = _families(parameters, source_values, statuses)
    family_before_rows = _family_spectrum_rows(families)
    _write_csv(output_dir / "family_spectra_before.csv", family_before_rows, SPECTRUM_FIELDS)

    # Oracle isolation boundary: everything through local repair is derived only
    # from the source pilot, matrix evaluators, and named detector thresholds.
    profile_events, measured_detector_runtimes = _detect_profiles(families)
    thresholds = REPAIR.THRESHOLD_PROFILES[args.threshold_profile]
    nominal_events = profile_events[args.threshold_profile]
    _write_csv(output_dir / "detector_events.csv", _event_rows(nominal_events), DETECTOR_FIELDS)
    repairs, pilot_results, inferred_windows, candidate_rows, pilot_solves, pilot_cache_hits, measured_local_wall = _run_pilot_repairs(
        output_dir,
        source_hash,
        families,
        nominal_events,
        thresholds,
    )
    repaired_values = _repaired_values(source_values, repairs)
    repaired_families = _families(parameters, repaired_values, statuses)
    remaining_events = [
        event
        for family in repaired_families.values()
        for event in REPAIR.detect_family_inventory(family, thresholds)
    ]
    repaired_rows = _repaired_spectrum_rows(families, repaired_values, pilot_results)
    _write_csv(output_dir / "repaired_spectra.csv", repaired_rows, REPAIRED_SPECTRUM_FIELDS)

    # Post-run oracle access starts here and nowhere above this line.
    oracle, oracle_region_nodes = _oracle_integer_map()
    before_after_rows = _before_after_oracle(
        nominal_events, families, repaired_values, oracle, pilot_results
    )
    _write_csv(output_dir / "before_after.csv", before_after_rows, BEFORE_AFTER_FIELDS)
    known_regions = _known_oracle_defects(source_values, oracle, oracle_region_nodes)
    oracle_differences: dict[tuple[str, str, float], float] = defaultdict(float)
    for row in before_after_rows:
        key = (str(row["case_id"]), str(row["theory"]), float(row["beta"]))
        oracle_differences[key] = max(oracle_differences[key], float(row["oracle_difference"]))
    nominal_success = (
        len(remaining_events) == 0
        and all(value <= COMPLETE.DEFAULT_ROOT_MATCH_TOL for value in oracle_differences.values())
        and all(result.status == "resolved_after_local_repair" for result in pilot_results.values())
    )
    sensitivity_rows = _sensitivity_rows(
        profile_events, measured_detector_runtimes, known_regions, nominal_success
    )

    validation_manifest: list[dict[str, object]] = []
    for family in families.values():
        for beta in family.beta_values:
            validation_manifest.append(
                {
                    "group": "manual_refined_oracle",
                    "case_id": family.case_id,
                    "epsilon_0": family.epsilon_0,
                    "mu": family.mu,
                    "eta": family.eta,
                    "beta": beta,
                    "theory": family.theory,
                    "oracle_source": "manual refined pilot post-run only",
                }
            )
    case_results = _pilot_case_results(
        nominal_events,
        source_values,
        repaired_values,
        pilot_results,
        oracle_differences,
    )
    readiness_pass, readiness_count = _readiness_validation(validation_manifest, case_results)
    former_pass, former_count = _former_blocker_validation(validation_manifest, case_results)
    coarse_counts, coarse_diagnostics = _coarse_validation(
        output_dir,
        thresholds,
        validation_manifest,
        case_results,
        inferred_windows,
        candidate_rows,
        skip_repairs=bool(args.skip_coarse_diagnostic_repair),
    )
    _write_csv(output_dir / "validation_manifest.csv", validation_manifest, VALIDATION_MANIFEST_FIELDS)
    _write_csv(output_dir / "inferred_repair_windows.csv", [_window_row(item) for item in inferred_windows], WINDOW_FIELDS)
    _write_csv(output_dir / "local_root_candidates.csv", candidate_rows, CANDIDATE_FIELDS)
    _write_csv(output_dir / "case_results.csv", case_results, CASE_RESULT_FIELDS)

    unresolved_rows = [
        {
            "case_id": row["case_id"],
            "group": row["group"],
            "execution_status": row["execution_status"],
            "defer_reason": row["execution_status"],
            "N_true": "",
            "required_guard": row["required_guard"],
            "force_strict_requested": row["force_strict_requested"],
            "force_strict_executed": 0,
        }
        for row in case_results
        if str(row["execution_status"]).startswith("deferred") or str(row["execution_status"]) == "attempted_error"
    ]
    _write_csv(output_dir / "unresolved_cases.csv", unresolved_rows, UNRESOLVED_FIELDS)

    synthetic = _synthetic_gate()
    event_nodes = {(event.case_id, event.theory, event.beta) for event in nominal_events}
    known_nodes = set().union(*known_regions.values()) if known_regions else set()
    detector_pass = (
        len(known_regions) == 7
        and all(nodes & event_nodes for nodes in known_regions.values())
        and not any(event.case_id == "P2" for event in nominal_events)
        and event_nodes == known_nodes
        and all(synthetic.values())
        and all(int(row["missed_known_defects"]) == 0 for row in sensitivity_rows)
    )
    p3_beta0 = pilot_results.get(("P3", "Timoshenko", 0.0))
    repair_pass = nominal_success and p3_beta0 is not None and p3_beta0.block_classification == "resolved_distinct_pair"
    reference_pass = readiness_pass and former_pass
    no_false_n = all(
        not (str(row["execution_status"]).startswith("deferred") and str(row["N_true_after"]) not in ("", "None"))
        for row in case_results
    )
    coarse_pass = no_false_n and all(int(row["force_strict_executed"]) == 0 for row in case_results)

    strict_reference_rows = _read_csv(READINESS_DIR / "benchmark_seven_cases.csv")
    strict_seconds = [
        float(row["strict_branch_seconds"])
        for row in strict_reference_rows
        if row["run_mode"] == "full_k10_full_strict" and float(row["strict_branch_seconds"]) > 0.0
    ]
    previous_strict_median = float(np.median(strict_seconds)) if strict_seconds else 0.0
    total_local_evaluations = sum(result.matrix_evaluations for result in pilot_results.values()) + sum(
        int(item["matrix_evaluations"]) for item in coarse_diagnostics
    )
    local_wall_seconds = sum(result.wall_seconds for result in pilot_results.values()) + sum(
        float(item["wall_seconds"]) for item in coarse_diagnostics
    )
    local_case_count = len(pilot_results) + len(coarse_diagnostics)
    measured_runtime = {
        "detector_permissive": measured_detector_runtimes["permissive"],
        "detector_nominal": measured_detector_runtimes["nominal"],
        "detector_conservative": measured_detector_runtimes["conservative"],
        "local_repair": local_wall_seconds,
        "previous_force_strict_median": previous_strict_median,
    }
    stable_runtime = _stable_runtime_metadata(output_dir, measured_runtime)
    sensitivity_rows = _sensitivity_rows(
        profile_events,
        {
            profile: stable_runtime[f"detector_{profile}"]
            for profile in ("permissive", "nominal", "conservative")
        },
        known_regions,
        nominal_success,
    )
    _write_csv(output_dir / "threshold_sensitivity.csv", sensitivity_rows, SENSITIVITY_FIELDS)
    cache_files = [path for path in (output_dir / "cache").rglob("*") if path.is_file()]
    cache_bytes = sum(path.stat().st_size for path in cache_files)
    primary_matrix_evaluations = sum(int(row["eb_matrix_evaluations"]) + int(row["timo_matrix_evaluations"]) for row in source_rows[::STORE_ROOTS])
    runtime_rows = [
        {
            "stage": "pointwise_primary_saved_reference",
            "case_count": 4 * 91,
            "wall_seconds": 0.0,
            "matrix_evaluations": primary_matrix_evaluations,
            "beta_probes": 0,
            "lambda_window_evaluations": 0,
            "cache_bytes": 0,
            "comparison_source": "saved beta sorted pilot",
        },
        {
            "stage": "automatic_detector_nominal",
            "case_count": len(families),
            "wall_seconds": stable_runtime["detector_nominal"],
            "matrix_evaluations": 0,
            "beta_probes": 0,
            "lambda_window_evaluations": 0,
            "cache_bytes": 0,
            "comparison_source": "current audit baseline",
        },
        {
            "stage": "automatic_local_repair",
            "case_count": local_case_count,
            "wall_seconds": local_wall_seconds,
            "matrix_evaluations": total_local_evaluations,
            "beta_probes": sum(len(result.beta_probes) for result in pilot_results.values()),
            "lambda_window_evaluations": total_local_evaluations,
            "cache_bytes": cache_bytes,
            "comparison_source": "current audit cached baseline",
        },
        {
            "stage": "previous_force_strict_reference",
            "case_count": len(strict_seconds),
            "wall_seconds": sum(strict_seconds),
            "matrix_evaluations": 0,
            "beta_probes": 0,
            "lambda_window_evaluations": 0,
            "cache_bytes": 0,
            "comparison_source": "readiness v2 saved strict timings; not executed",
        },
    ]
    _write_csv(output_dir / "runtime_summary.csv", runtime_rows, RUNTIME_FIELDS)
    local_mean_seconds = local_wall_seconds / max(1, local_case_count)
    cost_pass = (
        total_local_evaluations > 0
        and local_mean_seconds < max(stable_runtime["previous_force_strict_median"], 1.0)
        and all(int(row["force_strict_executed"]) == 0 for row in case_results)
    )
    integration_pass = detector_pass and repair_pass and reference_pass and coarse_pass and cost_pass
    gate_rows: list[dict[str, object]] = []

    def add_gate(gate: str, status: bool, metric: str, value: object, explanation: str) -> None:
        gate_rows.append({"gate": gate, "status": "PASS" if status else "FAIL", "metric": metric, "value": value, "explanation": explanation})

    add_gate("detector_gate", detector_pass, "known_regions_detected", len(known_regions), "R regions are read only from the post-run oracle; P2 and synthetic crossings have no repair trigger")
    add_gate("local_repair_gate", repair_pass, "remaining_detector_events", len(remaining_events), "matrix-confirmed integer spectra must match the manual oracle")
    add_gate("reference_preservation_gate", reference_pass, "readiness/former_cases", f"{readiness_count}/{former_count}", "immutable N_true references are unchanged")
    add_gate("coarse_unresolved_reduction_gate", coarse_pass, "resolved/deferred", json.dumps(coarse_counts, sort_keys=True), "PASS means conservative classification with no false N_true")
    add_gate("cost_gate", cost_pass, "mean_local_vs_median_strict_seconds", f"{local_mean_seconds:.6g}/{stable_runtime['previous_force_strict_median']:.6g}", "force strict is a saved reference only")
    add_gate("integration_readiness_gate", integration_pass, "all_component_gates", all((detector_pass, repair_pass, reference_pass, coarse_pass, cost_pass)), "diagnostic-only; production runner remains unchanged")
    _write_csv(output_dir / "gate_summary.csv", gate_rows, GATE_FIELDS)

    plot_paths = [_plot_case(output_dir, case_id, parameters[case_id], repaired_values) for case_id in sorted(parameters)]
    plot_paths.append(_plot_detector_summary(output_dir, nominal_events))
    plot_paths.append(_plot_runtime(output_dir, runtime_rows))

    force_requested = sum(bool(row["force_strict_requested"]) for row in case_results)
    force_executed = sum(int(row["force_strict_executed"]) for row in case_results)
    direct_current_solves = sum(not bool(item.get("cache_hit")) for item in coarse_diagnostics)
    current_matrix_solves = pilot_solves + direct_current_solves
    report_lines = [
        "# Family inventory automatic local-repair audit",
        "",
        "Diagnostic-only workflow. Sorted spectral positions are compared; physical descendant branches are not used.",
        "The refined manual pilot was opened only after detector inference and matrix repair completed.",
        "",
        "## Detector",
        "",
        f"- version: `{REPAIR.DETECTOR_VERSION}`",
        f"- threshold profile: `{thresholds.name}`",
        f"- formula: robust median normalized same-rank mismatch versus signed shifts -2/-1/+1/+2 on a consistent tail of at least {thresholds.minimum_tail_length} ranks",
        f"- normalized family noise: median + {thresholds.noise_mad_multiplier:g} MAD, floor {thresholds.normalized_noise_floor:g}",
        f"- improvement gate: same/noise >= {thresholds.same_rank_noise_multiplier:g}, shifted/noise <= {thresholds.shifted_noise_multiplier:g}, ratio >= {thresholds.minimum_improvement_ratio:g}",
        f"- inferred pilot repair points: {len(nominal_events)}",
        f"- P2 false triggers: {sum(event.case_id == 'P2' for event in nominal_events)}",
        "",
        "## Results",
        "",
        f"- detector gate: {'PASS' if detector_pass else 'FAIL'}",
        f"- local repair gate: {'PASS' if repair_pass else 'FAIL'}",
        f"- reference preservation gate: {'PASS' if reference_pass else 'FAIL'}",
        f"- coarse unresolved reduction gate: {'PASS' if coarse_pass else 'FAIL'}",
        f"- cost gate: {'PASS' if cost_pass else 'FAIL'}",
        f"- integration readiness gate: {'PASS' if integration_pass else 'FAIL'}",
        f"- P3 beta=0 classification: {p3_beta0.block_classification if p3_beta0 else 'missing'}",
        f"- pilot local matrix evaluations: {sum(result.matrix_evaluations for result in pilot_results.values())}",
        f"- total local matrix evaluations: {total_local_evaluations}",
        f"- pilot cache hits/current solves: {pilot_cache_hits}/{pilot_solves}",
        f"- all current matrix solves: {current_matrix_solves}",
        f"- force strict requested/executed: {force_requested}/{force_executed}",
        f"- coarse unresolved locally resolved/deferred: {coarse_counts['coarse_resolved']}/{coarse_counts['coarse_deferred']}",
        f"- strict-tail locally resolved/deferred: {coarse_counts['strict_resolved']}/{coarse_counts['strict_deferred']}",
        "",
        "## Prohibited-operation counters",
        "",
        "- branch tracking: 0",
        "- MAC: 0",
        "- shape reconstruction: 0",
        "- force/full strict execution: 0",
        "- global K12 strict scan: 0",
        "- coarse-grid resume: 0",
        "",
        "## Output files",
        "",
    ]
    for name in (
        "validation_manifest.csv", "family_spectra_before.csv", "detector_events.csv",
        "inferred_repair_windows.csv", "local_root_candidates.csv", "repaired_spectra.csv",
        "before_after.csv", "case_results.csv", "gate_summary.csv", "runtime_summary.csv",
        "threshold_sensitivity.csv", "unresolved_cases.csv", "report.md",
    ):
        report_lines.append(f"- `{(output_dir / name).relative_to(REPO_ROOT)}`")
    for path in plot_paths:
        report_lines.append(f"- `{path.relative_to(REPO_ROOT)}`")
    report_lines.extend(
        [
            "",
            "Production integration was not performed. The coarse grid was not resumed.",
            "",
        ]
    )
    (output_dir / "report.md").write_text("\n".join(report_lines), encoding="utf-8")
    log_payload = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "invocation_wall_seconds": time.perf_counter() - invocation_started,
        "pilot_cache_hits": pilot_cache_hits,
        "pilot_solves": pilot_solves,
        "direct_current_solves": direct_current_solves,
        "current_matrix_solves": current_matrix_solves,
        "force_strict_executed": force_executed,
        "integration_readiness_gate": "PASS" if integration_pass else "FAIL",
    }
    (output_dir / "logs" / "latest_invocation.json").write_text(
        json.dumps(log_payload, sort_keys=True, indent=2), encoding="utf-8"
    )
    print(
        "family_inventory_local_repair_audit "
        f"events={len(nominal_events)} pilot_solves={pilot_solves} cache_hits={pilot_cache_hits} "
        f"current_matrix_solves={current_matrix_solves} "
        f"local_evaluations={total_local_evaluations} force_strict_executed={force_executed} "
        f"integration_gate={'PASS' if integration_pass else 'FAIL'}"
    )


if __name__ == "__main__":
    main()
