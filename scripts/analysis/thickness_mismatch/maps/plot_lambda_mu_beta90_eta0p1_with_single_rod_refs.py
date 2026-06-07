from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import sys
from typing import Sequence

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

from my_project.analytic.FreqMuNet import roots_clamped_supported  # noqa: E402
from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    find_first_n_roots_eta,
    thickness_mismatch_factors,
)
from my_project.analytic.solvers import fixed_fixed_lambdas  # noqa: E402


DEFAULT_BETA_DEG = 90.0
DEFAULT_EPSILON = 0.0025
DEFAULT_ETA = 0.1
DEFAULT_MU_MIN = 0.0
DEFAULT_MU_MAX = 0.8
DEFAULT_MU_STEP = 0.01
DEFAULT_N_SYSTEM_ROOTS = 6
DEFAULT_N_REFERENCE_ROOTS = 6
DEFAULT_N_REFERENCE_ROOTS_PLOTTED = 4
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "lambda_mu_beta90_eta0p1_single_rod_refs"
SMOKE_OUTPUT_DIR = REPO_ROOT / "results" / "_smoke" / "lambda_mu_beta90_eta0p1_single_rod_refs"

DEFAULT_ROOT_SCAN_STEP = 0.02
RETRY_ROOT_SCAN_STEP = 0.01
DEFAULT_ROOT_LMAX0 = 35.0
RETRY_ROOT_LMAX0 = 45.0

AMBIGUOUS_ABS_TOL = 1.0e-3
AMBIGUOUS_REL_TOL = 1.0e-3

SYSTEM_CSV_NAME = "lambda_mu_beta90_eta0p1_system_in_plane_sorted.csv"
REFERENCE_CSV_NAME = "lambda_mu_beta90_eta0p1_single_rod_references.csv"
MATCHES_CSV_NAME = "lambda_mu_beta90_eta0p1_nearest_reference_matches.csv"
ALL_CANDIDATES_CSV_NAME = "lambda_mu_beta90_eta0p1_all_candidate_matches.csv"
BRANCH_SUMMARY_CSV_NAME = "lambda_mu_beta90_eta0p1_branch_family_summary.csv"
SPIKE_AUDIT_CSV_NAME = "lambda_mu_beta90_eta0p1_spike_audit.csv"
MAIN_PNG_NAME = "lambda_mu_beta90_eta0p1_in_plane_vs_single_rod_refs.png"
COMPACT_PNG_NAME = "lambda_mu_beta90_eta0p1_in_plane_vs_single_rod_refs_compact.png"
REPORT_NAME = "lambda_mu_beta90_eta0p1_report.md"

SYSTEM_FIELDS = [
    "mu",
    "beta_deg",
    "epsilon",
    "eta",
    "sorted_index",
    "Lambda_system",
    "root_solver_warning",
    "notes",
]

REFERENCE_FIELDS = [
    "mu",
    "beta_deg",
    "epsilon",
    "eta",
    "rod_id",
    "rod_label",
    "length_factor",
    "tau",
    "boundary_condition",
    "reference_mode_index",
    "alpha_root",
    "Lambda_reference",
]

MATCH_FIELDS = [
    "mu",
    "system_sorted_index",
    "Lambda_system",
    "best_rod_id",
    "best_rod_label",
    "best_boundary_condition",
    "best_reference_mode_index",
    "Lambda_reference",
    "abs_diff",
    "rel_diff_system",
    "rel_diff_reference",
    "ambiguous",
    "notes",
]

ALL_CANDIDATE_FIELDS = [
    "mu",
    "system_sorted_index",
    "Lambda_system",
    "candidate_rank_by_abs_diff",
    "rod_id",
    "rod_label",
    "boundary_condition",
    "reference_mode_index",
    "Lambda_reference",
    "abs_diff",
    "rel_diff_system",
    "rel_diff_reference",
    "notes",
]

BRANCH_SUMMARY_FIELDS = [
    "system_sorted_index",
    "dominant_reference_family",
    "fraction_of_mu_points",
    "number_of_family_switches",
    "switch_mu_locations",
    "notes",
]

SPIKE_AUDIT_FIELDS = [
    "mu",
    "sorted_index",
    "Lambda_raw",
    "Lambda_retry",
    "action",
    "notes",
]


@dataclass(frozen=True)
class Args:
    beta_deg: float
    epsilon: float
    eta: float
    mu_min: float
    mu_max: float
    mu_step: float
    n_system_roots: int
    n_reference_roots: int
    n_reference_roots_plotted: int
    output_dir: Path
    smoke: bool


@dataclass(frozen=True)
class ReferenceValue:
    mu: float
    rod_id: int
    rod_label: str
    length_factor: float
    tau: float
    boundary_condition: str
    reference_mode_index: int
    alpha_root: float
    Lambda_reference: float

    @property
    def family(self) -> str:
        return f"{self.rod_label}_{self.boundary_condition}"


def _fmt(value: object) -> str:
    try:
        value_f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not isfinite(value_f):
        return "nan"
    return f"{value_f:.16g}"


def _output_dir(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def _write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def mu_grid(args: Args) -> np.ndarray:
    if args.smoke:
        return np.array([0.0, 0.2, 0.4, 0.6, 0.8], dtype=float)
    if args.mu_step <= 0.0:
        raise ValueError("mu-step must be positive.")
    values = np.arange(args.mu_min, args.mu_max + 0.5 * args.mu_step, args.mu_step, dtype=float)
    if values.size == 0 or abs(float(values[-1]) - args.mu_max) > 1.0e-10:
        values = np.append(values, args.mu_max)
    values[0] = args.mu_min
    values[-1] = args.mu_max
    return np.unique(np.round(values, 12))


def validate_args(args: Args, mu_values: np.ndarray) -> None:
    if not isfinite(args.beta_deg):
        raise ValueError("beta-deg must be finite.")
    if not isfinite(args.epsilon) or args.epsilon <= 0.0:
        raise ValueError("epsilon must be positive and finite.")
    if not (-1.0 < args.eta < 1.0):
        raise ValueError("eta must lie inside (-1, 1).")
    if args.n_system_roots <= 0:
        raise ValueError("n-system-roots must be positive.")
    if args.n_reference_roots <= 0:
        raise ValueError("n-reference-roots must be positive.")
    if args.n_reference_roots_plotted <= 0:
        raise ValueError("n-reference-roots-plotted must be positive.")
    if args.n_reference_roots_plotted > args.n_reference_roots:
        raise ValueError("n-reference-roots-plotted must not exceed n-reference-roots.")
    if mu_values.size == 0:
        raise ValueError("mu grid is empty.")
    if np.any(mu_values <= -1.0) or np.any(mu_values >= 1.0):
        raise ValueError("mu values must lie inside (-1, 1).")
    if np.any(np.diff(mu_values) <= 0.0):
        raise ValueError("mu grid must be strictly increasing.")
    for mu in mu_values:
        thickness_mismatch_factors(float(mu), args.eta)


def rod_state(mu: float, eta: float, rod_id: int) -> tuple[str, float, float]:
    factors = thickness_mismatch_factors(float(mu), float(eta))
    length1 = 1.0 - float(mu)
    length2 = 1.0 + float(mu)
    if rod_id == 1:
        if length1 <= length2:
            return "short", length1, factors.tau1
        return "long", length1, factors.tau1
    if rod_id == 2:
        if length2 >= length1:
            return "long", length2, factors.tau2
        return "short", length2, factors.tau2
    raise ValueError(f"unsupported rod id {rod_id!r}")


def alpha_roots(n_roots: int) -> dict[str, np.ndarray]:
    roots = {
        "clamped_pinned": np.asarray(roots_clamped_supported(int(n_roots)), dtype=float),
        "clamped_clamped": np.asarray(fixed_fixed_lambdas(int(n_roots)), dtype=float),
    }
    for name, values in roots.items():
        if values.size < int(n_roots) or np.any(~np.isfinite(values)) or np.any(np.diff(values) <= 0.0):
            raise RuntimeError(f"{name} alpha roots must be finite and increasing.")
    return roots


def build_reference_values(args: Args, mu_values: np.ndarray) -> tuple[list[ReferenceValue], dict[int, list[ReferenceValue]]]:
    roots = alpha_roots(args.n_reference_roots)
    values: list[ReferenceValue] = []
    by_mu_index: dict[int, list[ReferenceValue]] = {}
    for mu_index, mu in enumerate(mu_values):
        local: list[ReferenceValue] = []
        for rod_id in (1, 2):
            rod_label, length_factor, tau = rod_state(float(mu), args.eta, rod_id)
            if not (isfinite(length_factor) and length_factor > 0.0 and isfinite(tau) and tau > 0.0):
                raise ValueError(f"invalid reference length/tau at mu={float(mu):g}, rod_id={rod_id}.")
            for boundary_condition, alphas in roots.items():
                lambdas = alphas * np.sqrt(tau) / length_factor
                if np.any(~np.isfinite(lambdas)) or np.any(lambdas <= 0.0) or np.any(np.diff(lambdas) <= 0.0):
                    raise RuntimeError(
                        "reference frequencies must be finite, positive, and increasing "
                        f"at mu={float(mu):g}, rod_id={rod_id}, {boundary_condition}."
                    )
                for mode_index, (alpha, value) in enumerate(zip(alphas, lambdas), start=1):
                    item = ReferenceValue(
                        mu=float(mu),
                        rod_id=int(rod_id),
                        rod_label=rod_label,
                        length_factor=float(length_factor),
                        tau=float(tau),
                        boundary_condition=boundary_condition,
                        reference_mode_index=int(mode_index),
                        alpha_root=float(alpha),
                        Lambda_reference=float(value),
                    )
                    values.append(item)
                    local.append(item)
        by_mu_index[int(mu_index)] = local
    return values, by_mu_index


def solve_roots_once(args: Args, mu: float, *, scan_step: float, lmax0: float) -> np.ndarray:
    return np.asarray(
        find_first_n_roots_eta(
            float(np.deg2rad(args.beta_deg)),
            float(mu),
            args.epsilon,
            args.eta,
            args.n_system_roots,
            scan_step=float(scan_step),
            Lmax0=float(lmax0),
        ),
        dtype=float,
    )


def solve_system_roots(args: Args, mu_values: np.ndarray) -> tuple[np.ndarray, list[list[str]], list[list[str]]]:
    roots_grid = np.full((len(mu_values), args.n_system_roots), np.nan, dtype=float)
    warnings = [["none" for _ in range(args.n_system_roots)] for _ in range(len(mu_values))]
    notes = [["sorted in-plane Euler-Bernoulli root" for _ in range(args.n_system_roots)] for _ in range(len(mu_values))]
    for mu_index, mu in enumerate(mu_values):
        try:
            roots = solve_roots_once(args, float(mu), scan_step=DEFAULT_ROOT_SCAN_STEP, lmax0=DEFAULT_ROOT_LMAX0)
        except Exception as exc:  # pragma: no cover - diagnostic CSV path.
            roots = np.full(args.n_system_roots, np.nan, dtype=float)
            for mode_index in range(args.n_system_roots):
                warnings[mu_index][mode_index] = "root_solver_exception"
                notes[mu_index][mode_index] = f"default root solve failed: {exc}"
        if roots.size < args.n_system_roots:
            padded = np.full(args.n_system_roots, np.nan, dtype=float)
            padded[: roots.size] = roots
            roots = padded
        roots = np.sort(roots[: args.n_system_roots])
        if np.count_nonzero(np.isfinite(roots)) < args.n_system_roots:
            try:
                retry = solve_roots_once(args, float(mu), scan_step=RETRY_ROOT_SCAN_STEP, lmax0=RETRY_ROOT_LMAX0)
                if retry.size < args.n_system_roots:
                    padded = np.full(args.n_system_roots, np.nan, dtype=float)
                    padded[: retry.size] = retry
                    retry = padded
                retry = np.sort(retry[: args.n_system_roots])
                if np.count_nonzero(np.isfinite(retry)) > np.count_nonzero(np.isfinite(roots)):
                    roots = retry
                    for mode_index in range(args.n_system_roots):
                        warnings[mu_index][mode_index] = "default_missing_retry_success"
                        notes[mu_index][mode_index] = "stricter root search used after missing default roots"
            except Exception as exc:  # pragma: no cover - diagnostic CSV path.
                for mode_index in range(args.n_system_roots):
                    warnings[mu_index][mode_index] = "missing_root_retry_exception"
                    notes[mu_index][mode_index] = f"retry root solve failed: {exc}"
        for mode_index, value in enumerate(roots[: args.n_system_roots]):
            if not (isfinite(float(value)) and float(value) > 0.0):
                warnings[mu_index][mode_index] = (
                    "missing_root_after_retry" if warnings[mu_index][mode_index] == "none" else warnings[mu_index][mode_index]
                )
                notes[mu_index][mode_index] = "missing or nonpositive root; plot uses NaN"
        finite = roots[np.isfinite(roots)]
        if finite.size > 1 and np.any(np.diff(finite) <= 0.0):
            for mode_index in range(args.n_system_roots):
                warnings[mu_index][mode_index] = "nonincreasing_sorted_roots"
                notes[mu_index][mode_index] = "root list is not strictly increasing where finite"
        roots_grid[mu_index] = roots[: args.n_system_roots]
    return roots_grid, warnings, notes


def suspicious_isolated_value(prev_value: float, value: float, next_value: float) -> bool:
    if not all(isfinite(x) for x in (prev_value, value, next_value)):
        return False
    scale = max(1.0, abs(prev_value), abs(value), abs(next_value))
    neighbors_close = abs(next_value - prev_value) <= max(0.75, 0.10 * scale)
    current_far = min(abs(value - prev_value), abs(value - next_value)) >= max(1.25, 0.22 * scale)
    return bool(neighbors_close and current_far)


def audit_and_repair_spikes(
    args: Args,
    mu_values: np.ndarray,
    roots_grid: np.ndarray,
    warnings: list[list[str]],
    notes: list[list[str]],
) -> tuple[np.ndarray, list[dict[str, object]]]:
    plot_grid = np.array(roots_grid, dtype=float)
    audit_rows: list[dict[str, object]] = []
    for mode_index in range(args.n_system_roots):
        for mu_index in range(1, len(mu_values) - 1):
            prev_value = float(roots_grid[mu_index - 1, mode_index])
            value = float(roots_grid[mu_index, mode_index])
            next_value = float(roots_grid[mu_index + 1, mode_index])
            if not suspicious_isolated_value(prev_value, value, next_value):
                continue
            retry_roots = solve_roots_once(
                args,
                float(mu_values[mu_index]),
                scan_step=RETRY_ROOT_SCAN_STEP,
                lmax0=RETRY_ROOT_LMAX0,
            )
            if retry_roots.size < args.n_system_roots:
                padded = np.full(args.n_system_roots, np.nan, dtype=float)
                padded[: retry_roots.size] = retry_roots
                retry_roots = padded
            retry_value = float(np.sort(retry_roots[: args.n_system_roots])[mode_index])
            if not suspicious_isolated_value(prev_value, retry_value, next_value):
                roots_grid[mu_index, mode_index] = retry_value
                plot_grid[mu_index, mode_index] = retry_value
                warnings[mu_index][mode_index] = "suspicious_spike_retry_success"
                notes[mu_index][mode_index] = "raw value replaced after stricter root-search retry"
                action = "replaced_with_retry_value"
            else:
                plot_grid[mu_index, mode_index] = np.nan
                warnings[mu_index][mode_index] = "unresolved_suspicious_spike"
                notes[mu_index][mode_index] = "raw value retained in CSV but omitted from plot"
                action = "plot_nan_keep_raw_csv"
            audit_rows.append(
                {
                    "mu": _fmt(mu_values[mu_index]),
                    "sorted_index": mode_index + 1,
                    "Lambda_raw": _fmt(value),
                    "Lambda_retry": _fmt(retry_value),
                    "action": action,
                    "notes": notes[mu_index][mode_index],
                }
            )
    return plot_grid, audit_rows


def system_rows(
    args: Args,
    mu_values: np.ndarray,
    roots_grid: np.ndarray,
    warnings: Sequence[Sequence[str]],
    notes: Sequence[Sequence[str]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for mu_index, mu in enumerate(mu_values):
        for mode_index in range(args.n_system_roots):
            rows.append(
                {
                    "mu": _fmt(mu),
                    "beta_deg": _fmt(args.beta_deg),
                    "epsilon": _fmt(args.epsilon),
                    "eta": _fmt(args.eta),
                    "sorted_index": mode_index + 1,
                    "Lambda_system": _fmt(roots_grid[mu_index, mode_index]),
                    "root_solver_warning": warnings[mu_index][mode_index],
                    "notes": notes[mu_index][mode_index],
                }
            )
    return rows


def reference_rows(args: Args, values: Sequence[ReferenceValue]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in values:
        rows.append(
            {
                "mu": _fmt(item.mu),
                "beta_deg": _fmt(args.beta_deg),
                "epsilon": _fmt(args.epsilon),
                "eta": _fmt(args.eta),
                "rod_id": item.rod_id,
                "rod_label": item.rod_label,
                "length_factor": _fmt(item.length_factor),
                "tau": _fmt(item.tau),
                "boundary_condition": item.boundary_condition,
                "reference_mode_index": item.reference_mode_index,
                "alpha_root": _fmt(item.alpha_root),
                "Lambda_reference": _fmt(item.Lambda_reference),
            }
        )
    return rows


def compare_candidates(
    args: Args,
    mu_values: np.ndarray,
    roots_grid: np.ndarray,
    references_by_mu: dict[int, list[ReferenceValue]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    match_rows: list[dict[str, object]] = []
    all_rows: list[dict[str, object]] = []
    for mu_index, mu in enumerate(mu_values):
        references = references_by_mu[int(mu_index)]
        for mode_index in range(args.n_system_roots):
            system_value = float(roots_grid[mu_index, mode_index])
            candidate_rows: list[dict[str, object]] = []
            for reference in references:
                if isfinite(system_value) and system_value > 0.0:
                    abs_diff = abs(system_value - reference.Lambda_reference)
                    rel_diff_system = abs_diff / system_value
                    rel_diff_reference = abs_diff / reference.Lambda_reference
                    notes = "frequency-proximity diagnostic only"
                else:
                    abs_diff = float("nan")
                    rel_diff_system = float("nan")
                    rel_diff_reference = float("nan")
                    notes = "comparison unavailable because system root is nonfinite"
                candidate_rows.append(
                    {
                        "mu": _fmt(mu),
                        "system_sorted_index": mode_index + 1,
                        "Lambda_system": _fmt(system_value),
                        "rod_id": reference.rod_id,
                        "rod_label": reference.rod_label,
                        "boundary_condition": reference.boundary_condition,
                        "reference_mode_index": reference.reference_mode_index,
                        "Lambda_reference": _fmt(reference.Lambda_reference),
                        "abs_diff": _fmt(abs_diff),
                        "rel_diff_system": _fmt(rel_diff_system),
                        "rel_diff_reference": _fmt(rel_diff_reference),
                        "notes": notes,
                        "_abs_diff": abs_diff,
                        "_rel_diff_system": rel_diff_system,
                    }
                )
            candidate_rows.sort(
                key=lambda row: (
                    float(row["_abs_diff"]) if isfinite(float(row["_abs_diff"])) else float("inf"),
                    str(row["rod_label"]),
                    str(row["boundary_condition"]),
                    int(row["reference_mode_index"]),
                )
            )
            for rank, row in enumerate(candidate_rows, start=1):
                clean = dict(row)
                clean["candidate_rank_by_abs_diff"] = rank
                clean.pop("_abs_diff", None)
                clean.pop("_rel_diff_system", None)
                all_rows.append(clean)
            finite_rows = [row for row in candidate_rows if isfinite(float(row["_abs_diff"]))]
            if not finite_rows:
                match_rows.append(
                    {
                        "mu": _fmt(mu),
                        "system_sorted_index": mode_index + 1,
                        "Lambda_system": _fmt(system_value),
                        "best_rod_id": "",
                        "best_rod_label": "",
                        "best_boundary_condition": "",
                        "best_reference_mode_index": "",
                        "Lambda_reference": "nan",
                        "abs_diff": "nan",
                        "rel_diff_system": "nan",
                        "rel_diff_reference": "nan",
                        "ambiguous": "yes",
                        "notes": "no finite system/reference comparison",
                    }
                )
                continue
            best = finite_rows[0]
            best_abs = float(best["_abs_diff"])
            best_rel = float(best["_rel_diff_system"])
            ambiguous = [best]
            for row in finite_rows[1:]:
                if abs(float(row["_abs_diff"]) - best_abs) < AMBIGUOUS_ABS_TOL or abs(
                    float(row["_rel_diff_system"]) - best_rel
                ) < AMBIGUOUS_REL_TOL:
                    ambiguous.append(row)
                else:
                    break
            if len(ambiguous) > 1:
                ambiguous_flag = "yes"
                labels = ", ".join(
                    f"{row['rod_label']}_{row['boundary_condition']}_m{row['reference_mode_index']}"
                    for row in ambiguous[:4]
                )
                notes = f"ambiguous nearest candidates: {labels}"
            else:
                ambiguous_flag = "no"
                notes = "unique closest candidate by absolute difference"
            match_rows.append(
                {
                    "mu": _fmt(mu),
                    "system_sorted_index": mode_index + 1,
                    "Lambda_system": best["Lambda_system"],
                    "best_rod_id": best["rod_id"],
                    "best_rod_label": best["rod_label"],
                    "best_boundary_condition": best["boundary_condition"],
                    "best_reference_mode_index": best["reference_mode_index"],
                    "Lambda_reference": best["Lambda_reference"],
                    "abs_diff": best["abs_diff"],
                    "rel_diff_system": best["rel_diff_system"],
                    "rel_diff_reference": best["rel_diff_reference"],
                    "ambiguous": ambiguous_flag,
                    "notes": notes,
                }
            )
    return match_rows, all_rows


def family_from_match(row: dict[str, object]) -> str:
    rod_label = str(row.get("best_rod_label", ""))
    bc = str(row.get("best_boundary_condition", ""))
    if not rod_label or not bc:
        return ""
    return f"{rod_label}_{bc}"


def branch_family_summary(args: Args, mu_values: np.ndarray, match_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    by_mode: dict[int, list[dict[str, object]]] = {mode: [] for mode in range(1, args.n_system_roots + 1)}
    for row in match_rows:
        by_mode[int(row["system_sorted_index"])].append(row)
    summary: list[dict[str, object]] = []
    for mode in range(1, args.n_system_roots + 1):
        rows = sorted(by_mode[mode], key=lambda row: float(row["mu"]))
        families = [family_from_match(row) for row in rows if family_from_match(row)]
        if not families:
            summary.append(
                {
                    "system_sorted_index": mode,
                    "dominant_reference_family": "",
                    "fraction_of_mu_points": "nan",
                    "number_of_family_switches": 0,
                    "switch_mu_locations": "",
                    "notes": "no finite nearest-family assignments",
                }
            )
            continue
        counts = Counter(families)
        dominant, count = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
        switches: list[str] = []
        previous_family = ""
        previous_mu = float("nan")
        for row in rows:
            family = family_from_match(row)
            if not family:
                continue
            mu = float(row["mu"])
            if previous_family and family != previous_family:
                switches.append(f"between {_fmt(previous_mu)} and {_fmt(mu)}")
            previous_family = family
            previous_mu = mu
        summary.append(
            {
                "system_sorted_index": mode,
                "dominant_reference_family": dominant,
                "fraction_of_mu_points": _fmt(count / len(families)),
                "number_of_family_switches": len(switches),
                "switch_mu_locations": "; ".join(switches),
                "notes": "sorted-index curve summary; not descendant branch tracking",
            }
        )
    return summary


def reference_curve_arrays(
    args: Args,
    mu_values: np.ndarray,
    references_by_mu: dict[int, list[ReferenceValue]],
) -> dict[tuple[str, str], np.ndarray]:
    arrays: dict[tuple[str, str], np.ndarray] = {}
    for rod_label in ("long", "short"):
        for boundary_condition in ("clamped_pinned", "clamped_clamped"):
            arrays[(rod_label, boundary_condition)] = np.full(
                (args.n_reference_roots, len(mu_values)),
                np.nan,
                dtype=float,
            )
    for mu_index in range(len(mu_values)):
        for item in references_by_mu[int(mu_index)]:
            arrays[(item.rod_label, item.boundary_condition)][item.reference_mode_index - 1, mu_index] = (
                item.Lambda_reference
            )
    return arrays


def family_style(rod_label: str, boundary_condition: str) -> dict[str, object]:
    styles = {
        ("long", "clamped_pinned"): {"color": "#1f77b4", "linestyle": (0, (6.0, 2.2)), "alpha": 0.58},
        ("long", "clamped_clamped"): {"color": "#2ca02c", "linestyle": (0, (7.0, 2.0, 1.4, 2.0)), "alpha": 0.60},
        ("short", "clamped_pinned"): {"color": "#ff7f0e", "linestyle": (0, (3.0, 2.0)), "alpha": 0.52},
        ("short", "clamped_clamped"): {"color": "#d62728", "linestyle": ":", "alpha": 0.60},
    }
    return styles[(rod_label, boundary_condition)]


def plot_curves(
    path: Path,
    args: Args,
    mu_values: np.ndarray,
    system_plot_grid: np.ndarray,
    reference_arrays: dict[tuple[str, str], np.ndarray],
    *,
    n_reference_roots_plotted: int,
    compact: bool,
) -> None:
    fig, ax = plt.subplots(figsize=(10.8, 6.3))
    system_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for mode_index in range(args.n_system_roots):
        ax.plot(
            mu_values,
            system_plot_grid[:, mode_index],
            color=system_colors[mode_index % len(system_colors)],
            linewidth=2.1,
            linestyle="-",
            label=f"system sorted {mode_index + 1}" if mode_index == 0 else None,
            zorder=5,
        )
    for (rod_label, boundary_condition), values in reference_arrays.items():
        style = family_style(rod_label, boundary_condition)
        for ref_index in range(min(int(n_reference_roots_plotted), values.shape[0])):
            plot_values = np.array(values[ref_index], dtype=float)
            if compact:
                finite_system = system_plot_grid[np.isfinite(system_plot_grid)]
                if finite_system.size:
                    y_limit = 1.18 * float(np.max(finite_system))
                    plot_values = np.where(plot_values <= y_limit, plot_values, np.nan)
            ax.plot(
                mu_values,
                plot_values,
                color=style["color"],
                linewidth=1.25,
                linestyle=style["linestyle"],
                alpha=float(style["alpha"]),
                zorder=2,
            )
    finite_system = system_plot_grid[np.isfinite(system_plot_grid)]
    finite_refs = np.concatenate([values[:n_reference_roots_plotted].ravel() for values in reference_arrays.values()])
    finite_refs = finite_refs[np.isfinite(finite_refs)]
    if compact and finite_system.size:
        ax.set_ylim(0.0, 1.18 * float(np.max(finite_system)))
    elif finite_system.size or finite_refs.size:
        max_value = max(
            float(np.max(finite_system)) if finite_system.size else 0.0,
            float(np.max(finite_refs)) if finite_refs.size else 0.0,
        )
        ax.set_ylim(0.0, 1.08 * max_value)
    ax.set_xlabel("mu")
    ax.set_ylabel("Lambda")
    ref_count = int(n_reference_roots_plotted)
    title_suffix = "compact zoom" if compact else "full reference scale"
    ax.set_title(
        "Lambda(mu), beta = 90 deg, epsilon = 0.0025, eta = 0.1\n"
        f"solid: system in-plane sorted frequencies; dashed: single-rod references ({ref_count} modes plotted, {title_suffix})"
    )
    ax.grid(True, color="0.88", linewidth=0.65)
    handles = [
        Line2D([0], [0], color="black", lw=2.1, ls="-", label=f"system sorted roots, first {args.n_system_roots}"),
    ]
    for rod_label, boundary_condition in (
        ("long", "clamped_pinned"),
        ("long", "clamped_clamped"),
        ("short", "clamped_pinned"),
        ("short", "clamped_clamped"),
    ):
        style = family_style(rod_label, boundary_condition)
        label = f"{rod_label} {'CP' if boundary_condition == 'clamped_pinned' else 'CC'} refs"
        handles.append(
            Line2D(
                [0],
                [0],
                color=style["color"],
                lw=1.4,
                ls=style["linestyle"],
                label=label,
                alpha=style["alpha"],
            )
        )
    ax.legend(handles=handles, loc="upper left", fontsize=8.5, frameon=False, handlelength=3.6)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> list[str]:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return lines


def write_report(
    path: Path,
    args: Args,
    mu_values: np.ndarray,
    system_grid: np.ndarray,
    system_warnings: Sequence[Sequence[str]],
    spike_rows: Sequence[dict[str, object]],
    branch_summary_rows: Sequence[dict[str, object]],
    output_paths: Sequence[Path],
) -> None:
    finite_system = system_grid[np.isfinite(system_grid)]
    missing_count = int(np.count_nonzero(~np.isfinite(system_grid)))
    warning_count = sum(1 for row in system_warnings for item in row if item != "none")
    low_modes = [row for row in branch_summary_rows if int(row["system_sorted_index"]) <= min(3, args.n_system_roots)]
    high_modes = [row for row in branch_summary_rows if int(row["system_sorted_index"]) > min(3, args.n_system_roots)]
    long_low = sum(1 for row in low_modes if str(row["dominant_reference_family"]).startswith("long_"))
    short_high = sum(1 for row in high_modes if str(row["dominant_reference_family"]).startswith("short_"))
    lines: list[str] = [
        "# Lambda(mu) Beta 90 Eta 0.1 Single-Rod Reference Diagnostic",
        "",
        "This is a diagnostic-only analytic plot and proximity comparison. It is not an article figure.",
        "",
        "## Parameters",
        "",
        f"- beta_deg: {_fmt(args.beta_deg)}",
        f"- epsilon: {_fmt(args.epsilon)}",
        f"- eta: {_fmt(args.eta)}",
        f"- mu_min: {_fmt(float(mu_values[0]))}",
        f"- mu_max: {_fmt(float(mu_values[-1]))}",
        f"- mu grid points: {len(mu_values)}",
        f"- nominal mu step: {_fmt(args.mu_step if not args.smoke else 0.2)}",
        f"- system sorted roots per mu: {args.n_system_roots}",
        f"- reference roots saved per family: {args.n_reference_roots}",
        f"- reference roots plotted in main figure per family: {args.n_reference_roots_plotted}",
        "",
        "## Formulas",
        "",
        "- system roots: existing analytic in-plane Euler-Bernoulli thickness-mismatch root finder",
        "- reference CP/FP roots: tan(alpha) = tanh(alpha)",
        "- reference CC/FF roots: cosh(alpha) cos(alpha) = 1",
        "- reference scaling: Lambda_ref_i,n(mu) = alpha_n * sqrt(tau_i(mu, eta)) / L_i(mu)",
        "- axial frequencies are not included",
        "",
        "## Root And Spike Checks",
        "",
        f"- finite system roots: {int(finite_system.size)} / {system_grid.size}",
        f"- missing system roots: {missing_count}",
        f"- system root warning entries: {warning_count}",
        f"- suspicious spike audit rows: {len(spike_rows)}",
        "",
        "## Dominant Nearest Reference Families",
        "",
    ]
    lines.extend(
        markdown_table(
            [
                "sorted index",
                "dominant family",
                "fraction",
                "switches",
                "switch locations",
            ],
            [
                [
                    row["system_sorted_index"],
                    row["dominant_reference_family"],
                    row["fraction_of_mu_points"],
                    row["number_of_family_switches"],
                    row["switch_mu_locations"] or "none",
                ]
                for row in branch_summary_rows
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Qualitative Summary",
            "",
            f"- lower sorted-index curves with long-rod dominant references: {long_low} / {len(low_modes)}",
            f"- higher sorted-index curves with short-rod dominant references: {short_high} / {len(high_modes)}",
            "- Family switches are nearest-frequency assignment changes only; they are not descendant branch changes.",
            "- This is a frequency-proximity diagnostic only, not a proof that a system mode is physically identical to an isolated-rod mode.",
            "- No crossing or no-crossing claim is made.",
            "",
            "## Output Files",
            "",
        ]
    )
    for output_path in output_paths:
        lines.append(f"- {output_path.name}")
    lines.extend(
        [
            "",
            "## Protected Scope",
            "",
            "The script uses sorted in-plane analytic frequencies only. It does not use descendant tracking, FEM, 3D FEM, Gmsh, CalculiX, article workspaces, old determinant edits, old solver edits, baseline-result edits, or analytic formula changes.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        description="Plot diagnostic Lambda(mu) in-plane sorted roots at beta=90 with single-rod references."
    )
    parser.add_argument("--beta-deg", type=float, default=DEFAULT_BETA_DEG)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--mu-min", type=float, default=DEFAULT_MU_MIN)
    parser.add_argument("--mu-max", type=float, default=DEFAULT_MU_MAX)
    parser.add_argument("--mu-step", type=float, default=DEFAULT_MU_STEP)
    parser.add_argument("--n-system-roots", type=int, default=DEFAULT_N_SYSTEM_ROOTS)
    parser.add_argument("--n-reference-roots", type=int, default=DEFAULT_N_REFERENCE_ROOTS)
    parser.add_argument("--n-reference-roots-plotted", type=int, default=DEFAULT_N_REFERENCE_ROOTS_PLOTTED)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(argv)
    n_system_roots = int(ns.n_system_roots)
    n_reference_roots = int(ns.n_reference_roots)
    n_reference_roots_plotted = int(ns.n_reference_roots_plotted)
    output_dir = Path(ns.output_dir)
    if bool(ns.smoke):
        n_system_roots = 4
        n_reference_roots = 4
        n_reference_roots_plotted = 3
        if output_dir == DEFAULT_OUTPUT_DIR:
            output_dir = SMOKE_OUTPUT_DIR
    return Args(
        beta_deg=float(ns.beta_deg),
        epsilon=float(ns.epsilon),
        eta=float(ns.eta),
        mu_min=float(ns.mu_min),
        mu_max=float(ns.mu_max),
        mu_step=float(ns.mu_step),
        n_system_roots=n_system_roots,
        n_reference_roots=n_reference_roots,
        n_reference_roots_plotted=n_reference_roots_plotted,
        output_dir=_output_dir(output_dir),
        smoke=bool(ns.smoke),
    )


def run(args: Args) -> dict[str, object]:
    mu_values = mu_grid(args)
    validate_args(args, mu_values)
    references, references_by_mu = build_reference_values(args, mu_values)
    system_grid, warnings, notes = solve_system_roots(args, mu_values)
    system_plot_grid, spike_rows = audit_and_repair_spikes(args, mu_values, system_grid, warnings, notes)
    matches, all_candidates = compare_candidates(args, mu_values, system_grid, references_by_mu)
    branch_summary = branch_family_summary(args, mu_values, matches)
    reference_arrays = reference_curve_arrays(args, mu_values, references_by_mu)

    output_dir = args.output_dir
    system_csv = output_dir / SYSTEM_CSV_NAME
    reference_csv = output_dir / REFERENCE_CSV_NAME
    matches_csv = output_dir / MATCHES_CSV_NAME
    all_candidates_csv = output_dir / ALL_CANDIDATES_CSV_NAME
    branch_summary_csv = output_dir / BRANCH_SUMMARY_CSV_NAME
    spike_audit_csv = output_dir / SPIKE_AUDIT_CSV_NAME
    main_png = output_dir / MAIN_PNG_NAME
    compact_png = output_dir / COMPACT_PNG_NAME
    report_md = output_dir / REPORT_NAME

    _write_csv(system_csv, system_rows(args, mu_values, system_grid, warnings, notes), SYSTEM_FIELDS)
    _write_csv(reference_csv, reference_rows(args, references), REFERENCE_FIELDS)
    _write_csv(matches_csv, matches, MATCH_FIELDS)
    _write_csv(all_candidates_csv, all_candidates, ALL_CANDIDATE_FIELDS)
    _write_csv(branch_summary_csv, branch_summary, BRANCH_SUMMARY_FIELDS)
    if spike_rows:
        _write_csv(spike_audit_csv, spike_rows, SPIKE_AUDIT_FIELDS)
    plot_curves(
        main_png,
        args,
        mu_values,
        system_plot_grid,
        reference_arrays,
        n_reference_roots_plotted=args.n_reference_roots_plotted,
        compact=False,
    )
    compact_plotted = max(1, min(args.n_reference_roots_plotted, 3))
    plot_curves(
        compact_png,
        args,
        mu_values,
        system_plot_grid,
        reference_arrays,
        n_reference_roots_plotted=compact_plotted,
        compact=True,
    )
    output_paths = [system_csv, reference_csv, matches_csv, all_candidates_csv, branch_summary_csv, main_png, compact_png]
    if spike_rows:
        output_paths.append(spike_audit_csv)
    output_paths.append(report_md)
    write_report(report_md, args, mu_values, system_grid, warnings, spike_rows, branch_summary, output_paths)
    return {
        "mu_values": mu_values,
        "system_grid": system_grid,
        "system_plot_grid": system_plot_grid,
        "warnings": warnings,
        "spike_rows": spike_rows,
        "branch_summary": branch_summary,
        "outputs": output_paths,
        "system_csv": system_csv,
        "reference_csv": reference_csv,
        "matches_csv": matches_csv,
        "all_candidates_csv": all_candidates_csv,
        "branch_summary_csv": branch_summary_csv,
        "main_png": main_png,
        "compact_png": compact_png,
        "report_md": report_md,
        "spike_audit_csv": spike_audit_csv if spike_rows else None,
    }


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    result = run(args)
    mu_values = np.asarray(result["mu_values"], dtype=float)
    system_grid = np.asarray(result["system_grid"], dtype=float)
    warning_count = sum(1 for row in result["warnings"] for item in row if item != "none")
    missing_count = int(np.count_nonzero(~np.isfinite(system_grid)))
    print("diagnostic-only Lambda(mu) map with single-rod references")
    print(f"beta_deg={args.beta_deg:g}, epsilon={args.epsilon:g}, eta={args.eta:g}")
    print(
        f"mu grid: {len(mu_values)} points from {float(mu_values[0]):g} to {float(mu_values[-1]):g}; "
        f"step={'0.2 smoke grid' if args.smoke else f'{args.mu_step:g}'}"
    )
    print(f"system roots per mu={args.n_system_roots}, reference roots saved={args.n_reference_roots}")
    print(f"missing system roots={missing_count}, root warning entries={warning_count}")
    print(f"suspicious spike audit rows={len(result['spike_rows'])}")
    print("dominant nearest reference families:")
    for row in result["branch_summary"]:
        print(
            f"  sorted {row['system_sorted_index']}: {row['dominant_reference_family']} "
            f"(fraction={row['fraction_of_mu_points']}, switches={row['number_of_family_switches']}, "
            f"locations={row['switch_mu_locations'] or 'none'})"
        )
    for path in result["outputs"]:
        print(f"saved: {path}")
    print("solid curves: system in-plane sorted frequencies; dashed curves: single-rod references")
    print("no descendant tracking; no FEM/Gmsh/CalculiX; no article/formula/baseline changes")
    return result


if __name__ == "__main__":
    main()
