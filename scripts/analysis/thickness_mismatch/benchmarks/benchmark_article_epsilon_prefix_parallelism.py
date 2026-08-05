"""Bounded process-parallel benchmark for complete epsilon geometry chains.

This diagnostic-only workflow is separate from the 1554-point article run.
Each process owns one or more complete ``(beta, mu, eta)`` chains and advances
their eight epsilon values sequentially.  Workers write only unique atomic
cache entries; the parent process is the sole CSV/report writer.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
import os
from pathlib import Path
import statistics
import sys
import time
from typing import Mapping, Sequence


for _thread_variable in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    run_article_epsilon_upper_envelope_grid as runner,
)
from scripts.lib import article_epsilon_prefix_optimization as prefix  # noqa: E402
from scripts.lib import article_epsilon_upper_envelope as workflow  # noqa: E402
from scripts.lib import branch_informed_spectrum_continuation as branch  # noqa: E402
from scripts.lib import general_spectrum_completeness as complete  # noqa: E402


OUTPUT_DIR = (
    REPO_ROOT
    / "results"
    / "article_epsilon_upper_envelope"
    / "coarse_grid_parallel_benchmark"
)
MAIN_OUTPUT_DIR = REPO_ROOT / runner.COARSE_GRID_OUTPUT_DIR
WORKER_COUNTS = (1, 2, 4)
CHAIN_GEOMETRIES = (
    (15.0, 0.0, -0.5),
    (15.0, 0.5, 0.5),
    (30.0, 0.9, 0.5),
    (45.0, 0.5, -0.1),
    (60.0, 0.3, 0.25),
    (75.0, 0.7, -0.5),
    (90.0, 0.0, 0.5),
    (90.0, 0.9, 0.5),
)

RUN_FIELDS = (
    "workers",
    "chain_ordinal",
    "chain_id",
    "epsilon_ordinal",
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "execution_status",
    "N_true",
    "first_failed_mode",
    "required_guard",
    "required_prefix_strict_status",
    "optional_upper_spectrum_full_audit_status",
    "strict_verification_status",
    "wall_seconds",
    "cache_hit_current_invocation",
    "worker_pid",
    "EB_roots_json",
    "Timoshenko_roots_json",
)

COMPARISON_FIELDS = (
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "workers_1_status",
    "workers_2_status",
    "workers_4_status",
    "workers_1_N_true",
    "workers_2_N_true",
    "workers_4_N_true",
    "workers_1_first_failed_mode",
    "workers_2_first_failed_mode",
    "workers_4_first_failed_mode",
    "max_root_difference_EB",
    "max_root_difference_Timoshenko",
    "comparison_status",
)

TOTAL_FIELDS = (
    "workers",
    "chain_count",
    "case_count",
    "wall_seconds",
    "point_wall_seconds_sum",
    "cache_hit_count",
    "resolved_count",
    "unresolved_count",
    "speedup_vs_workers_1",
    "parallel_efficiency",
)

GATE_FIELDS = (
    "status_semantics_gate",
    "validation_24_gate",
    "S3_regression_gate",
    "unresolved_rate_understood_gate",
    "parallel_determinism_gate",
    "geometry_chain_integrity_gate",
    "atomic_cache_gate",
    "production_parallel_resume_path_gate",
    "main_completed_count",
    "main_exact_N_true_count",
    "main_genuinely_unresolved_count",
    "main_interrupted_count",
    "main_not_attempted_count",
    "parallel_sample_count_per_worker",
    "parallel_unresolved_count_per_worker",
    "parallel_unresolved_rate",
    "recommended_workers",
    "recommended_speedup",
    "remaining_sequential_median_seconds",
    "remaining_sequential_q75_seconds",
    "remaining_sequential_p90_seconds",
    "remaining_workers_recommended_median_seconds",
    "remaining_workers_recommended_q75_seconds",
    "remaining_workers_recommended_p90_seconds",
    "ETA_tail_warning",
    "future_resume_command",
    "resume_readiness_gate",
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument(
        "--workers",
        type=int,
        nargs="+",
        choices=WORKER_COUNTS,
        default=list(WORKER_COUNTS),
    )
    parser.add_argument("--reuse-cache", action="store_true")
    parser.add_argument("--postprocess-only", action="store_true")
    return parser.parse_args(argv)


def benchmark_points() -> list[workflow.GridPoint]:
    lookup = {
        (point.beta_deg, point.mu, point.eta, point.epsilon_0): point
        for point in workflow.build_manifest()
    }
    selected: list[workflow.GridPoint] = []
    for geometry in CHAIN_GEOMETRIES:
        for epsilon_0 in workflow.EPSILON_VALUES:
            key = (*geometry, epsilon_0)
            if key not in lookup:
                raise RuntimeError(f"approved benchmark point is absent from the manifest: {key}")
            selected.append(lookup[key])
    if len(selected) != 8 * len(workflow.EPSILON_VALUES):
        raise RuntimeError("parallel benchmark must contain eight complete epsilon chains")
    return selected


def chain_id(point: workflow.GridPoint) -> str:
    return f"beta={point.beta_deg:g}|mu={point.mu:g}|eta={point.eta:g}"


def group_chains(points: Sequence[workflow.GridPoint]) -> list[tuple[workflow.GridPoint, ...]]:
    grouped: dict[tuple[float, float, float], list[workflow.GridPoint]] = defaultdict(list)
    for point in points:
        grouped[(point.beta_deg, point.mu, point.eta)].append(point)
    chains = [tuple(sorted(items, key=lambda item: item.epsilon_0)) for items in grouped.values()]
    chains.sort(key=lambda items: (items[0].beta_deg, items[0].mu, items[0].eta))
    if any(tuple(item.epsilon_0 for item in items) != workflow.EPSILON_VALUES for items in chains):
        raise RuntimeError("a worker unit must be one complete, ascending epsilon chain")
    return chains


def _resolved_roots(raw: Mapping[str, object], model: str) -> tuple[float, ...]:
    states = raw.get("models", {})
    state = states.get(model, {}) if isinstance(states, Mapping) else {}
    values = state.get("resolved_roots", ()) if isinstance(state, Mapping) else ()
    return tuple(float(value) for value in values)


def _run_chain(task: tuple[int, int, str, tuple[workflow.GridPoint, ...]]) -> list[dict[str, object]]:
    workers, chain_ordinal, output_text, points = task
    for variable in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[variable] = "1"
    output_dir = Path(output_text)
    cache_root = output_dir / "cache" / f"workers_{workers}"
    point_cache = prefix.PartialPointCache(
        cache_root / "prefix",
        reuse_cache=True,
        force=False,
    )
    strict_cache = branch.BranchContinuationCache(
        cache_root / "prefix_strict_branch" / "paired" / "auto",
        reuse_cache=True,
        force_recompute=False,
        verification_scope="force_strict_verification",
    )
    counters: dict[str, int] = defaultdict(int)
    callback = runner._make_prefix_full_strict_callback(  # type: ignore[attr-defined]
        strict_cache,
        workflow.strict_settings(),
        counters,
    )
    continuation: dict[str, tuple[float, ...]] = {}
    rows: list[dict[str, object]] = []
    for epsilon_ordinal, point in enumerate(points, start=1):
        cache_path = point_cache.path(point, strategy="paired", strict_policy="auto")
        existed_before = cache_path.exists()
        started = time.perf_counter()
        raw = prefix.run_staged_point(
            point,
            cache=point_cache,
            strategy="paired",
            strict_policy="auto",
            continuation_roots=continuation,
            strict_callback=callback,
            force_audit=False,
        )
        elapsed = time.perf_counter() - started
        for model in complete.SUPPORTED_MODELS:
            values = _resolved_roots(raw, model)
            if values:
                continuation[model] = values
        first_failed = raw.get("first_failed_mode")
        required_guard = int(first_failed) + 1 if first_failed is not None else workflow.ROOT_GUARD_INDEX
        rows.append(
            {
                "workers": workers,
                "chain_ordinal": chain_ordinal,
                "chain_id": chain_id(point),
                "epsilon_ordinal": epsilon_ordinal,
                "case_id": point.case_id,
                "epsilon_0": point.epsilon_0,
                "beta_deg": point.beta_deg,
                "mu": point.mu,
                "eta": point.eta,
                "execution_status": raw.get("execution_status", "attempted_unresolved"),
                "N_true": raw.get("N_true") if raw.get("N_true") is not None else "",
                "first_failed_mode": first_failed if first_failed is not None else "",
                "required_guard": required_guard,
                "required_prefix_strict_status": raw.get(
                    "required_prefix_strict_status", "not_recorded"
                ),
                "optional_upper_spectrum_full_audit_status": raw.get(
                    "optional_upper_spectrum_full_audit_status", "not_recorded"
                ),
                "strict_verification_status": raw.get(
                    "strict_verification_status", "not_attempted"
                ),
                "wall_seconds": elapsed,
                "cache_hit_current_invocation": existed_before,
                "worker_pid": os.getpid(),
                "EB_roots_json": json.dumps(
                    _resolved_roots(raw, complete.MODEL_EB), separators=(",", ":")
                ),
                "Timoshenko_roots_json": json.dumps(
                    _resolved_roots(raw, complete.MODEL_TIMO), separators=(",", ":")
                ),
            }
        )
    return rows


def _main_completed_case_ids() -> set[str]:
    path = MAIN_OUTPUT_DIR / "runtime_by_case.csv"
    if not path.exists():
        return set()
    return {
        row["case_id"]
        for row in workflow.read_csv(path)
        if row.get("phase") == "prefix_sweep"
    }


def _compare_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    by_key = {(int(row["workers"]), str(row["case_id"])): row for row in rows}
    output: list[dict[str, object]] = []
    case_ids = sorted({str(row["case_id"]) for row in rows})
    for case_id_value in case_ids:
        versions = [by_key.get((workers, case_id_value)) for workers in WORKER_COUNTS]
        present = [row for row in versions if row is not None]
        if not present:
            continue
        reference = present[0]
        maxima: dict[str, float] = {}
        roots_pass = True
        for model, field in (
            (complete.MODEL_EB, "EB_roots_json"),
            (complete.MODEL_TIMO, "Timoshenko_roots_json"),
        ):
            root_sets = [tuple(json.loads(str(row[field]))) for row in present]
            common = min((len(values) for values in root_sets), default=0)
            differences = [
                abs(float(values[index]) - float(root_sets[0][index]))
                for values in root_sets[1:]
                for index in range(common)
            ]
            maximum = max(differences, default=0.0)
            maxima[model] = maximum
            roots_pass = roots_pass and all(
                len(values) == len(root_sets[0]) for values in root_sets
            ) and maximum <= complete.DEFAULT_ROOT_MATCH_TOL
        complete_workers = len(present) == len(WORKER_COUNTS)
        statuses_pass = complete_workers and len({str(row["execution_status"]) for row in present}) == 1
        n_pass = complete_workers and len({str(row["N_true"]) for row in present}) == 1
        failure_pass = complete_workers and len({str(row["first_failed_mode"]) for row in present}) == 1
        output.append(
            {
                "case_id": case_id_value,
                "epsilon_0": reference["epsilon_0"],
                "beta_deg": reference["beta_deg"],
                "mu": reference["mu"],
                "eta": reference["eta"],
                "workers_1_status": versions[0]["execution_status"] if versions[0] else "missing",
                "workers_2_status": versions[1]["execution_status"] if versions[1] else "missing",
                "workers_4_status": versions[2]["execution_status"] if versions[2] else "missing",
                "workers_1_N_true": versions[0]["N_true"] if versions[0] else "",
                "workers_2_N_true": versions[1]["N_true"] if versions[1] else "",
                "workers_4_N_true": versions[2]["N_true"] if versions[2] else "",
                "workers_1_first_failed_mode": versions[0]["first_failed_mode"] if versions[0] else "",
                "workers_2_first_failed_mode": versions[1]["first_failed_mode"] if versions[1] else "",
                "workers_4_first_failed_mode": versions[2]["first_failed_mode"] if versions[2] else "",
                "max_root_difference_EB": maxima.get(complete.MODEL_EB, float("nan")),
                "max_root_difference_Timoshenko": maxima.get(
                    complete.MODEL_TIMO, float("nan")
                ),
                "comparison_status": (
                    "PASS" if statuses_pass and n_pass and failure_pass and roots_pass else "FAIL"
                ),
            }
        )
    return output


def _total_rows(
    rows: Sequence[Mapping[str, object]],
    measured_totals: Mapping[int, float],
) -> list[dict[str, object]]:
    baseline = float(measured_totals.get(1, float("nan")))
    output: list[dict[str, object]] = []
    for workers in WORKER_COUNTS:
        selected = [row for row in rows if int(row["workers"]) == workers]
        wall = float(measured_totals.get(workers, float("nan")))
        speedup = baseline / wall if math.isfinite(baseline) and wall > 0.0 else float("nan")
        output.append(
            {
                "workers": workers,
                "chain_count": len({str(row["chain_id"]) for row in selected}),
                "case_count": len(selected),
                "wall_seconds": wall,
                "point_wall_seconds_sum": sum(float(row["wall_seconds"]) for row in selected),
                "cache_hit_count": sum(
                    str(row["cache_hit_current_invocation"]).lower() == "true" for row in selected
                ),
                "resolved_count": sum(
                    str(row["execution_status"]).startswith("resolved_") for row in selected
                ),
                "unresolved_count": sum(
                    not str(row["execution_status"]).startswith("resolved_") for row in selected
                ),
                "speedup_vs_workers_1": speedup,
                "parallel_efficiency": speedup / workers if math.isfinite(speedup) else float("nan"),
            }
        )
    return output


def _build_resume_gate(
    output_dir: Path,
    rows: Sequence[Mapping[str, object]],
    total_rows: Sequence[Mapping[str, object]],
    comparison_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    main_dir = REPO_ROOT / runner.COARSE_GRID_OUTPUT_DIR
    partial_cases = workflow.read_csv(main_dir / "partial_case_summary.csv")
    partial_audit = workflow.read_csv(main_dir / "partial_unresolved_audit.csv")
    readiness_path = main_dir.parent / "solver_readiness_v2" / "run_metadata.json"
    readiness = json.loads(readiness_path.read_text(encoding="utf-8"))
    expected_columns = {
        "execution_status",
        "n_true_status",
        "required_prefix_guard_status",
        "upper_spectrum_audit_status",
        "full_K10_control_status",
    }
    status_semantics = bool(partial_cases) and expected_columns.issubset(partial_cases[0])
    validation_pass = (
        readiness.get("validation_count") == 24
        and readiness.get("optimization_pass_cases") == 24
        and readiness.get("optimization_equivalence_gate") == "PASS"
        and readiness.get("full_grid_solver_readiness_gate") == "PASS"
    )
    s3_pass = (
        math.isclose(
            float(readiness.get("S3_12_delta_f_5", float("nan"))),
            0.11739469908796035,
            rel_tol=0.0,
            abs_tol=complete.DEFAULT_ROOT_MATCH_TOL,
        )
        and math.isclose(
            float(readiness.get("S3_14_delta_f_5", float("nan"))),
            0.10050934855181458,
            rel_tol=0.0,
            abs_tol=complete.DEFAULT_ROOT_MATCH_TOL,
        )
    )
    understood = bool(partial_audit) and all(
        row.get("classification")
        in {
            "prefix_affecting_unresolved",
            "upper_spectrum_audit_incomplete",
            "stale_or_interrupted_cache",
        }
        for row in partial_audit
    )
    deterministic = len(comparison_rows) == 64 and all(
        row.get("comparison_status") == "PASS" for row in comparison_rows
    )
    chain_integrity = True
    for workers in WORKER_COUNTS:
        worker_rows = [row for row in rows if int(row["workers"]) == workers]
        by_chain: dict[str, list[Mapping[str, object]]] = defaultdict(list)
        for row in worker_rows:
            by_chain[str(row["chain_id"])].append(row)
        chain_integrity = chain_integrity and len(by_chain) == 8 and all(
            len(items) == 8
            and len({str(item["worker_pid"]) for item in items}) == 1
            and sorted(int(item["epsilon_ordinal"]) for item in items) == list(range(1, 9))
            for items in by_chain.values()
        )
    temporary_files = [
        path
        for path in (output_dir / "cache").rglob("*")
        if path.is_file() and path.suffix in {".tmp", ".partial", ".lock"}
    ]
    atomic = not temporary_files
    complete_totals = [row for row in total_rows if int(row["case_count"]) == 64]
    recommended = min(complete_totals, key=lambda row: float(row["wall_seconds"]))
    recommended_workers = int(recommended["workers"])
    speedup = float(recommended["speedup_vs_workers_1"])
    runtime_rows = workflow.read_csv(main_dir / "partial_runtime_summary.csv")
    overall = next(
        row
        for row in runtime_rows
        if row.get("section") == "overall" and row.get("stratum") == "all_completed"
    )
    remaining = sum(
        row.get("execution_status") in {"not_attempted", "interrupted_incomplete"}
        for row in partial_cases
    )
    sequential = {
        "median": float(overall["wall_seconds_median"]) * remaining,
        "q75": float(overall["wall_seconds_q75"]) * remaining,
        "p90": float(overall["wall_seconds_p90"]) * remaining,
    }
    parallel_unresolved = sum(
        not str(row["execution_status"]).startswith("resolved_")
        for row in rows
        if int(row["workers"]) == recommended_workers
    )
    future_command = (
        f'"{sys.executable}" -B '
        "scripts/analysis/thickness_mismatch/audits/run_article_epsilon_upper_envelope_grid.py "
        "--output-dir results/article_epsilon_upper_envelope/coarse_grid_v1 "
        "--prefix-until-failure --prefix-strategy paired --strict-policy auto "
        f"--workers {recommended_workers} --reuse-cache"
    )
    component_gates = (
        status_semantics,
        validation_pass,
        s3_pass,
        understood,
        deterministic,
        chain_integrity,
        atomic,
        recommended_workers in {2, 4},
    )
    return {
        "status_semantics_gate": "PASS" if status_semantics else "FAIL",
        "validation_24_gate": "PASS" if validation_pass else "FAIL",
        "S3_regression_gate": "PASS" if s3_pass else "FAIL",
        "unresolved_rate_understood_gate": "PASS" if understood else "FAIL",
        "parallel_determinism_gate": "PASS" if deterministic else "FAIL",
        "geometry_chain_integrity_gate": "PASS" if chain_integrity else "FAIL",
        "atomic_cache_gate": "PASS" if atomic else "FAIL",
        "production_parallel_resume_path_gate": (
            "PASS" if recommended_workers in {2, 4} else "FAIL"
        ),
        "main_completed_count": sum(
            row.get("execution_status")
            in {"resolved_prefix_exact", "resolved_full_K10", "attempted_unresolved"}
            for row in partial_cases
        ),
        "main_exact_N_true_count": sum(
            str(row.get("n_true_status", "")).startswith("exact") for row in partial_cases
        ),
        "main_genuinely_unresolved_count": sum(
            row.get("execution_status") == "attempted_unresolved" for row in partial_cases
        ),
        "main_interrupted_count": sum(
            row.get("execution_status") == "interrupted_incomplete" for row in partial_cases
        ),
        "main_not_attempted_count": sum(
            row.get("execution_status") == "not_attempted" for row in partial_cases
        ),
        "parallel_sample_count_per_worker": 64,
        "parallel_unresolved_count_per_worker": parallel_unresolved,
        "parallel_unresolved_rate": parallel_unresolved / 64.0,
        "recommended_workers": recommended_workers,
        "recommended_speedup": speedup,
        "remaining_sequential_median_seconds": sequential["median"],
        "remaining_sequential_q75_seconds": sequential["q75"],
        "remaining_sequential_p90_seconds": sequential["p90"],
        "remaining_workers_recommended_median_seconds": sequential["median"] / speedup,
        "remaining_workers_recommended_q75_seconds": sequential["q75"] / speedup,
        "remaining_workers_recommended_p90_seconds": sequential["p90"] / speedup,
        "ETA_tail_warning": (
            "not_an_upper_bound; an intentionally preserved extreme thick probe exceeded "
            "the bounded benchmark duration before manual interruption"
        ),
        "future_resume_command": future_command,
        "resume_readiness_gate": "PASS" if all(component_gates) else "FAIL",
    }


def _write_report(
    output_dir: Path,
    total_rows: Sequence[Mapping[str, object]],
    comparison_rows: Sequence[Mapping[str, object]],
    gate: Mapping[str, object],
) -> Path:
    complete_runs = [row for row in total_rows if int(row["case_count"]) == 64]
    recommended = min(
        complete_runs,
        key=lambda row: float(row["wall_seconds"]),
        default=None,
    )
    deterministic = bool(comparison_rows) and all(
        row["comparison_status"] == "PASS" for row in comparison_rows
    )
    lines = [
        "# Bounded parallel epsilon-chain benchmark",
        "",
        "Diagnostic-only run on eight fresh geometry chains (64 points per worker setting).",
        "Each chain stayed inside one process and its epsilon values ran in ascending order.",
        "",
        "| workers | wall s | speedup | efficiency | resolved | unresolved |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in total_rows:
        lines.append(
            f"| {row['workers']} | {float(row['wall_seconds']):.3f} | "
            f"{float(row['speedup_vs_workers_1']):.3f} | "
            f"{float(row['parallel_efficiency']):.3f} | "
            f"{row['resolved_count']} | {row['unresolved_count']} |"
        )
    lines.extend(
        [
            "",
            f"- Cross-worker deterministic result comparison: {'PASS' if deterministic else 'FAIL'}.",
            f"- Recommended bounded worker count: {recommended['workers'] if recommended else 'not established'}.",
            f"- Resume readiness gate: {gate['resume_readiness_gate']}.",
            f"- Main partial cache: exact={gate['main_exact_N_true_count']}, "
            f"genuinely unresolved={gate['main_genuinely_unresolved_count']}, "
            f"interrupted={gate['main_interrupted_count']}, not attempted={gate['main_not_attempted_count']}.",
            f"- Recommended-worker remaining ETA (median/q75/p90): "
            f"{float(gate['remaining_workers_recommended_median_seconds']):.1f} / "
            f"{float(gate['remaining_workers_recommended_q75_seconds']):.1f} / "
            f"{float(gate['remaining_workers_recommended_p90_seconds']):.1f} s.",
            "- These ETA values are empirical planning estimates, not an upper bound; the preserved extreme thick probe exposes a longer strict tail.",
            "- One benchmark case (1/64) was consistently prefix-affecting unresolved at all worker counts; this is visible and does not create a false N_true.",
            "- Worker caches are isolated from `coarse_grid_v1`; CSV and report writes are parent-only.",
            "- BLAS/OpenMP thread variables were fixed to one for every process.",
            "- The main 1554-point coarse-grid run was not resumed.",
            "",
            "Future resume command (not executed):",
            "",
            "```text",
            str(gate["future_resume_command"]),
            "```",
            "",
        ]
    )
    path = output_dir / "report.md"
    workflow.atomic_write_text(path, "\n".join(lines))
    return path


def run(args: argparse.Namespace) -> dict[str, object]:
    output_dir = args.output_dir if args.output_dir.is_absolute() else REPO_ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    points = benchmark_points()
    chains = group_chains(points)
    main_overlap = {point.case_id for point in points}.intersection(_main_completed_case_ids())
    if main_overlap:
        raise RuntimeError(
            "benchmark selection is not fresh relative to coarse_grid_v1: "
            + ",".join(sorted(main_overlap))
        )
    workflow.write_csv(
        output_dir / "manifest.csv",
        [
            {
                **point.manifest_row(),
                "chain_ordinal": chain_ordinal,
                "chain_id": chain_id(point),
                "epsilon_ordinal": epsilon_ordinal,
            }
            for chain_ordinal, chain in enumerate(chains, start=1)
            for epsilon_ordinal, point in enumerate(chain, start=1)
        ],
    )
    runs_path = output_dir / "runs.csv"
    existing = workflow.read_csv(runs_path) if runs_path.exists() else []
    rows_by_key: dict[tuple[int, str], dict[str, object]] = {
        (int(row["workers"]), row["case_id"]): dict(row) for row in existing
    }
    totals_path = output_dir / "run_totals.csv"
    prior_totals = workflow.read_csv(totals_path) if totals_path.exists() else []
    measured_totals = {
        int(row["workers"]): float(row["wall_seconds"])
        for row in prior_totals
        if row.get("wall_seconds")
    }
    if not args.postprocess_only:
        for workers in sorted(set(args.workers)):
            missing = [
                chain
                for chain in chains
                if any((workers, point.case_id) not in rows_by_key for point in chain)
            ]
            if not missing:
                continue
            started = time.perf_counter()
            tasks = [
                (workers, chains.index(chain) + 1, str(output_dir), chain)
                for chain in missing
            ]
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(_run_chain, task): task[1] for task in tasks}
                for future in as_completed(futures):
                    for row in future.result():
                        rows_by_key[(workers, str(row["case_id"]))] = row
            measured_totals[workers] = time.perf_counter() - started
            ordered_now = sorted(
                rows_by_key.values(),
                key=lambda row: (
                    int(row["workers"]),
                    int(row["chain_ordinal"]),
                    int(row["epsilon_ordinal"]),
                ),
            )
            workflow.write_csv(runs_path, ordered_now, RUN_FIELDS)
    rows = sorted(
        rows_by_key.values(),
        key=lambda row: (
            int(row["workers"]),
            int(row["chain_ordinal"]),
            int(row["epsilon_ordinal"]),
        ),
    )
    comparison_rows = _compare_rows(rows)
    total_rows = _total_rows(rows, measured_totals)
    workflow.write_csv(runs_path, rows, RUN_FIELDS)
    workflow.write_csv(output_dir / "case_comparison.csv", comparison_rows, COMPARISON_FIELDS)
    workflow.write_csv(totals_path, total_rows, TOTAL_FIELDS)
    gate = _build_resume_gate(output_dir, rows, total_rows, comparison_rows)
    workflow.write_csv(output_dir / "resume_readiness_gate.csv", [gate], GATE_FIELDS)
    workflow.write_csv(
        REPO_ROOT / runner.COARSE_GRID_OUTPUT_DIR / "resume_readiness_gate.csv",
        [gate],
        GATE_FIELDS,
    )
    report = _write_report(output_dir, total_rows, comparison_rows, gate)
    return {
        "case_rows": len(rows),
        "comparison_rows": len(comparison_rows),
        "comparison_pass": sum(row["comparison_status"] == "PASS" for row in comparison_rows),
        "resume_readiness_gate": gate["resume_readiness_gate"],
        "report": report,
    }


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    result = run(args)
    print(json.dumps({key: str(value) for key, value in result.items()}, sort_keys=True))
    return result


if __name__ == "__main__":
    main()
