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
if str(SCRIPT_PATH.parent) not in sys.path:
    sys.path.insert(0, str(SCRIPT_PATH.parent))

import plot_lambda_mu_beta90_eta0p1_with_single_rod_refs as single_mu  # noqa: E402
from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    find_first_n_roots_eta,
    thickness_mismatch_factors,
)


DEFAULT_BETA_DEG = 90.0
DEFAULT_EPSILON = 0.0025
DEFAULT_MU_VALUES = (0.0, 0.3, 0.6)
DEFAULT_ETA_MIN = -0.5
DEFAULT_ETA_MAX = 0.5
DEFAULT_ETA_STEP = 0.01
DEFAULT_N_SYSTEM_ROOTS = 6
DEFAULT_N_REFERENCE_ROOTS = 6
DEFAULT_READABLE_YMIN = 0.0
DEFAULT_READABLE_YMAX = 13.0
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "lambda_eta_beta90_mu_scan_single_rod_refs"
SMOKE_OUTPUT_DIR = REPO_ROOT / "results" / "_smoke" / "lambda_eta_beta90_mu_scan_single_rod_refs"

DEFAULT_ROOT_SCAN_STEP = 0.02
RETRY_ROOT_SCAN_STEP = 0.01
DEFAULT_ROOT_LMAX0 = 35.0
RETRY_ROOT_LMAX0 = 45.0

AMBIGUOUS_ABS_TOL = 1.0e-3
AMBIGUOUS_REL_TOL = 1.0e-3

COMBINED_SUMMARY_NAME = "mu_scan_summary.csv"
COMBINED_REPORT_NAME = "mu_scan_report.md"

SYSTEM_FIELDS = [
    "eta",
    "mu",
    "beta_deg",
    "epsilon",
    "sorted_index",
    "Lambda_system",
    "root_solver_warning",
    "notes",
]

REFERENCE_FIELDS = [
    "eta",
    "mu",
    "beta_deg",
    "epsilon",
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
    "eta",
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
    "eta",
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
    "mu",
    "system_sorted_index",
    "dominant_reference_family",
    "fraction_of_eta_points",
    "number_of_family_switches",
    "switch_eta_locations",
    "notes",
]

COMBINED_SUMMARY_FIELDS = [
    "mu",
    "system_sorted_index",
    "dominant_reference_family",
    "fraction_of_eta_points",
    "number_of_family_switches",
    "switch_eta_locations",
    "notes",
]


@dataclass(frozen=True)
class Args:
    beta_deg: float
    epsilon: float
    mu_values: tuple[float, ...]
    eta_min: float
    eta_max: float
    eta_step: float
    n_system_roots: int
    n_reference_roots: int
    readable_ymin: float
    readable_ymax: float
    output_dir: Path
    smoke: bool


@dataclass(frozen=True)
class ReferenceValue:
    eta_index: int
    eta: float
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


@dataclass(frozen=True)
class MuCaseResult:
    mu: float
    mu_token: str
    output_dir: Path
    eta_values: np.ndarray
    system_grid: np.ndarray
    system_plot_grid: np.ndarray
    warnings: Sequence[Sequence[str]]
    spike_rows: Sequence[dict[str, object]]
    branch_summary: Sequence[dict[str, object]]
    main_plot_stats: dict[str, object]
    full_plot_stats: dict[str, object]
    output_paths: Sequence[Path]


def _fmt(value: object) -> str:
    try:
        value_f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not isfinite(value_f):
        return "nan"
    return f"{value_f:.16g}"


def decimal_token(value: float) -> str:
    value_f = float(value)
    prefix = "m" if value_f < 0.0 else ""
    text = f"{abs(value_f):.10g}"
    if "e" in text or "E" in text:
        text = f"{abs(value_f):.10f}".rstrip("0").rstrip(".")
    if "." not in text:
        text = f"{text}.0"
    return prefix + text.replace(".", "p")


def angle_token(value: float) -> str:
    value_f = float(value)
    prefix = "m" if value_f < 0.0 else ""
    text = f"{abs(value_f):.10g}"
    if text.endswith(".0"):
        text = text[:-2]
    return prefix + text.replace(".", "p")


def mu_token(value: float) -> str:
    return f"mu_{decimal_token(float(value))}"


def file_prefix(args: Args, mu: float) -> str:
    return f"lambda_eta_beta{angle_token(args.beta_deg)}_mu_{decimal_token(float(mu))}"


def output_dir(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def eta_grid(args: Args) -> np.ndarray:
    if args.smoke:
        return np.array([args.eta_min, -0.25, 0.0, 0.25, args.eta_max], dtype=float)
    if args.eta_step <= 0.0:
        raise ValueError("eta-step must be positive.")
    values = np.arange(args.eta_min, args.eta_max + 0.5 * args.eta_step, args.eta_step, dtype=float)
    values = values[values <= args.eta_max + 1.0e-12]
    if values.size == 0 or not np.isclose(values[0], args.eta_min, rtol=0.0, atol=1.0e-12):
        values = np.insert(values, 0, args.eta_min)
    if not np.isclose(values[-1], args.eta_max, rtol=0.0, atol=1.0e-12):
        values = np.append(values, args.eta_max)
    values[0] = args.eta_min
    values[-1] = args.eta_max
    return np.unique(np.round(values, 12))


def validate_args(args: Args, eta_values: np.ndarray) -> None:
    if not isfinite(args.beta_deg):
        raise ValueError("beta-deg must be finite.")
    if not isfinite(args.epsilon) or args.epsilon <= 0.0:
        raise ValueError("epsilon must be positive and finite.")
    if not args.mu_values:
        raise ValueError("mu-values must not be empty.")
    if any(not (-1.0 < mu < 1.0) for mu in args.mu_values):
        raise ValueError("all mu values must lie inside (-1, 1).")
    if eta_values.size == 0:
        raise ValueError("eta grid is empty.")
    if np.any(eta_values <= -1.0) or np.any(eta_values >= 1.0):
        raise ValueError("eta values must lie inside (-1, 1).")
    if args.eta_max < args.eta_min:
        raise ValueError("eta-max must be greater than or equal to eta-min.")
    if args.n_system_roots <= 0:
        raise ValueError("n-system-roots must be positive.")
    if args.n_reference_roots <= 0:
        raise ValueError("n-reference-roots must be positive.")
    if not isfinite(args.readable_ymin) or not isfinite(args.readable_ymax):
        raise ValueError("readable y-limits must be finite.")
    if args.readable_ymax <= args.readable_ymin:
        raise ValueError("readable-ymax must be greater than readable-ymin.")
    for mu in args.mu_values:
        for eta in eta_values:
            thickness_mismatch_factors(float(mu), float(eta))


def solve_roots_once(args: Args, mu: float, eta: float, *, scan_step: float, lmax0: float) -> np.ndarray:
    return np.asarray(
        find_first_n_roots_eta(
            float(np.deg2rad(args.beta_deg)),
            float(mu),
            float(args.epsilon),
            float(eta),
            int(args.n_system_roots),
            scan_step=float(scan_step),
            Lmax0=float(lmax0),
        ),
        dtype=float,
    )


def solve_system_roots(
    args: Args,
    mu: float,
    eta_values: np.ndarray,
) -> tuple[np.ndarray, list[list[str]], list[list[str]]]:
    roots_grid = np.full((len(eta_values), args.n_system_roots), np.nan, dtype=float)
    warnings = [["none" for _ in range(args.n_system_roots)] for _ in range(len(eta_values))]
    notes = [["sorted in-plane Euler-Bernoulli root" for _ in range(args.n_system_roots)] for _ in range(len(eta_values))]
    for eta_index, eta in enumerate(eta_values):
        try:
            roots = solve_roots_once(args, mu, float(eta), scan_step=DEFAULT_ROOT_SCAN_STEP, lmax0=DEFAULT_ROOT_LMAX0)
        except Exception as exc:  # pragma: no cover - diagnostic CSV path.
            roots = np.full(args.n_system_roots, np.nan, dtype=float)
            for mode_index in range(args.n_system_roots):
                warnings[eta_index][mode_index] = "root_solver_exception"
                notes[eta_index][mode_index] = f"default root solve failed: {exc}"
        if roots.size < args.n_system_roots:
            padded = np.full(args.n_system_roots, np.nan, dtype=float)
            padded[: roots.size] = roots
            roots = padded
        roots = np.sort(roots[: args.n_system_roots])
        if np.count_nonzero(np.isfinite(roots)) < args.n_system_roots:
            try:
                retry = solve_roots_once(args, mu, float(eta), scan_step=RETRY_ROOT_SCAN_STEP, lmax0=RETRY_ROOT_LMAX0)
                if retry.size < args.n_system_roots:
                    padded = np.full(args.n_system_roots, np.nan, dtype=float)
                    padded[: retry.size] = retry
                    retry = padded
                retry = np.sort(retry[: args.n_system_roots])
                if np.count_nonzero(np.isfinite(retry)) > np.count_nonzero(np.isfinite(roots)):
                    roots = retry
                    for mode_index in range(args.n_system_roots):
                        warnings[eta_index][mode_index] = "default_missing_retry_success"
                        notes[eta_index][mode_index] = "stricter root search used after missing default roots"
            except Exception as exc:  # pragma: no cover - diagnostic CSV path.
                for mode_index in range(args.n_system_roots):
                    warnings[eta_index][mode_index] = "missing_root_retry_exception"
                    notes[eta_index][mode_index] = f"retry root solve failed: {exc}"
        for mode_index, value in enumerate(roots[: args.n_system_roots]):
            if not (isfinite(float(value)) and float(value) > 0.0):
                warnings[eta_index][mode_index] = (
                    "missing_root_after_retry"
                    if warnings[eta_index][mode_index] == "none"
                    else warnings[eta_index][mode_index]
                )
                notes[eta_index][mode_index] = "missing or nonpositive root; plot uses NaN"
        finite = roots[np.isfinite(roots)]
        if finite.size > 1 and np.any(np.diff(finite) <= 0.0):
            for mode_index in range(args.n_system_roots):
                warnings[eta_index][mode_index] = "nonincreasing_sorted_roots"
                notes[eta_index][mode_index] = "root list is not strictly increasing where finite"
        roots_grid[eta_index] = roots[: args.n_system_roots]
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
    mu: float,
    eta_values: np.ndarray,
    roots_grid: np.ndarray,
    warnings: list[list[str]],
    notes: list[list[str]],
) -> tuple[np.ndarray, list[dict[str, object]]]:
    plot_grid = np.array(roots_grid, dtype=float)
    spike_rows: list[dict[str, object]] = []
    for mode_index in range(args.n_system_roots):
        for eta_index in range(1, len(eta_values) - 1):
            prev_value = float(roots_grid[eta_index - 1, mode_index])
            value = float(roots_grid[eta_index, mode_index])
            next_value = float(roots_grid[eta_index + 1, mode_index])
            if not suspicious_isolated_value(prev_value, value, next_value):
                continue
            retry_roots = solve_roots_once(
                args,
                float(mu),
                float(eta_values[eta_index]),
                scan_step=RETRY_ROOT_SCAN_STEP,
                lmax0=RETRY_ROOT_LMAX0,
            )
            if retry_roots.size < args.n_system_roots:
                padded = np.full(args.n_system_roots, np.nan, dtype=float)
                padded[: retry_roots.size] = retry_roots
                retry_roots = padded
            retry_value = float(np.sort(retry_roots[: args.n_system_roots])[mode_index])
            if not suspicious_isolated_value(prev_value, retry_value, next_value):
                roots_grid[eta_index, mode_index] = retry_value
                plot_grid[eta_index, mode_index] = retry_value
                warnings[eta_index][mode_index] = "suspicious_spike_retry_success"
                notes[eta_index][mode_index] = "raw value replaced after stricter root-search retry"
                action = "replaced_with_retry_value"
            else:
                plot_grid[eta_index, mode_index] = np.nan
                warnings[eta_index][mode_index] = "unresolved_suspicious_spike"
                notes[eta_index][mode_index] = "raw value retained in CSV but omitted from plot"
                action = "plot_nan_keep_raw_csv"
            spike_rows.append(
                {
                    "eta": _fmt(eta_values[eta_index]),
                    "sorted_index": mode_index + 1,
                    "Lambda_raw": _fmt(value),
                    "Lambda_retry": _fmt(retry_value),
                    "action": action,
                    "notes": notes[eta_index][mode_index],
                }
            )
    return plot_grid, spike_rows


def build_reference_values(
    args: Args,
    mu: float,
    eta_values: np.ndarray,
) -> tuple[list[ReferenceValue], dict[int, list[ReferenceValue]]]:
    alpha = single_mu.alpha_roots(args.n_reference_roots)
    values: list[ReferenceValue] = []
    by_eta_index: dict[int, list[ReferenceValue]] = {}
    for eta_index, eta in enumerate(eta_values):
        local: list[ReferenceValue] = []
        for rod_id in (1, 2):
            rod_label, length_factor, tau = single_mu.rod_state(float(mu), float(eta), rod_id)
            if not (isfinite(length_factor) and length_factor > 0.0 and isfinite(tau) and tau > 0.0):
                raise ValueError(f"invalid reference length/tau at mu={mu:g}, eta={float(eta):g}, rod_id={rod_id}.")
            for boundary_condition, alphas in alpha.items():
                lambdas = alphas * np.sqrt(tau) / length_factor
                if np.any(~np.isfinite(lambdas)) or np.any(lambdas <= 0.0) or np.any(np.diff(lambdas) <= 0.0):
                    raise RuntimeError(
                        "reference frequencies must be finite, positive, and increasing "
                        f"at mu={mu:g}, eta={float(eta):g}, rod_id={rod_id}, {boundary_condition}."
                    )
                for mode_index, (alpha_root, value) in enumerate(zip(alphas, lambdas), start=1):
                    item = ReferenceValue(
                        eta_index=int(eta_index),
                        eta=float(eta),
                        mu=float(mu),
                        rod_id=int(rod_id),
                        rod_label=rod_label,
                        length_factor=float(length_factor),
                        tau=float(tau),
                        boundary_condition=boundary_condition,
                        reference_mode_index=int(mode_index),
                        alpha_root=float(alpha_root),
                        Lambda_reference=float(value),
                    )
                    values.append(item)
                    local.append(item)
        by_eta_index[int(eta_index)] = local
    return values, by_eta_index


def system_rows(
    args: Args,
    mu: float,
    eta_values: np.ndarray,
    roots_grid: np.ndarray,
    warnings: Sequence[Sequence[str]],
    notes: Sequence[Sequence[str]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for eta_index, eta in enumerate(eta_values):
        for mode_index in range(args.n_system_roots):
            rows.append(
                {
                    "eta": _fmt(eta),
                    "mu": _fmt(mu),
                    "beta_deg": _fmt(args.beta_deg),
                    "epsilon": _fmt(args.epsilon),
                    "sorted_index": mode_index + 1,
                    "Lambda_system": _fmt(roots_grid[eta_index, mode_index]),
                    "root_solver_warning": warnings[eta_index][mode_index],
                    "notes": notes[eta_index][mode_index],
                }
            )
    return rows


def reference_rows(args: Args, values: Sequence[ReferenceValue]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for item in values:
        rows.append(
            {
                "eta": _fmt(item.eta),
                "mu": _fmt(item.mu),
                "beta_deg": _fmt(args.beta_deg),
                "epsilon": _fmt(args.epsilon),
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
    mu: float,
    eta_values: np.ndarray,
    roots_grid: np.ndarray,
    references_by_eta: dict[int, list[ReferenceValue]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    match_rows: list[dict[str, object]] = []
    all_rows: list[dict[str, object]] = []
    for eta_index, eta in enumerate(eta_values):
        references = references_by_eta[int(eta_index)]
        for mode_index in range(args.n_system_roots):
            system_value = float(roots_grid[eta_index, mode_index])
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
                        "eta": _fmt(eta),
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
                        "_abs_diff_value": float(abs_diff),
                    }
                )
            ranked = sorted(
                candidate_rows,
                key=lambda row: (not isfinite(float(row["_abs_diff_value"])), float(row["_abs_diff_value"])),
            )
            for rank, row in enumerate(ranked, start=1):
                row = dict(row)
                row["candidate_rank_by_abs_diff"] = rank
                all_rows.append(row)
            best = ranked[0]
            if len(ranked) > 1 and isfinite(float(best["_abs_diff_value"])):
                second = ranked[1]
                tie_abs = abs(float(second["_abs_diff_value"]) - float(best["_abs_diff_value"]))
                scale = max(1.0, abs(float(system_value)))
                ambiguous = tie_abs <= max(AMBIGUOUS_ABS_TOL, AMBIGUOUS_REL_TOL * scale)
            else:
                ambiguous = False
            note_parts = [str(best["notes"])]
            if ambiguous:
                note_parts.append("nearest reference is ambiguous within tolerance")
            match_rows.append(
                {
                    "eta": _fmt(eta),
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
                    "ambiguous": "yes" if ambiguous else "no",
                    "notes": "; ".join(note_parts),
                }
            )
    return match_rows, all_rows


def family_from_match(row: dict[str, object]) -> str:
    rod_label = str(row.get("best_rod_label", ""))
    bc = str(row.get("best_boundary_condition", ""))
    if not rod_label or not bc:
        return ""
    return f"{rod_label}_{bc}"


def branch_family_summary(args: Args, mu: float, match_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    by_mode: dict[int, list[dict[str, object]]] = {mode: [] for mode in range(1, args.n_system_roots + 1)}
    for row in match_rows:
        by_mode[int(row["system_sorted_index"])].append(row)
    rows_out: list[dict[str, object]] = []
    for mode in range(1, args.n_system_roots + 1):
        rows = sorted(by_mode[mode], key=lambda row: float(row["eta"]))
        families = [family_from_match(row) for row in rows if family_from_match(row)]
        if not families:
            rows_out.append(
                {
                    "mu": _fmt(mu),
                    "system_sorted_index": mode,
                    "dominant_reference_family": "",
                    "fraction_of_eta_points": "nan",
                    "number_of_family_switches": 0,
                    "switch_eta_locations": "",
                    "notes": "no finite nearest-family assignments",
                }
            )
            continue
        counts = Counter(families)
        dominant, count = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0]
        switches: list[str] = []
        previous_family = ""
        previous_eta = float("nan")
        for row in rows:
            family = family_from_match(row)
            if not family:
                continue
            eta = float(row["eta"])
            if previous_family and family != previous_family:
                switches.append(f"between {_fmt(previous_eta)} and {_fmt(eta)}")
            previous_family = family
            previous_eta = eta
        rows_out.append(
            {
                "mu": _fmt(mu),
                "system_sorted_index": mode,
                "dominant_reference_family": dominant,
                "fraction_of_eta_points": _fmt(count / len(families)),
                "number_of_family_switches": len(switches),
                "switch_eta_locations": "; ".join(switches),
                "notes": "sorted-index curve summary; not descendant branch tracking",
            }
        )
    return rows_out


def reference_curve_arrays(
    args: Args,
    eta_values: np.ndarray,
    references_by_eta: dict[int, list[ReferenceValue]],
) -> dict[tuple[str, str], np.ndarray]:
    arrays: dict[tuple[str, str], np.ndarray] = {}
    for rod_label in ("long", "short"):
        for boundary_condition in ("clamped_pinned", "clamped_clamped"):
            arrays[(rod_label, boundary_condition)] = np.full(
                (args.n_reference_roots, len(eta_values)),
                np.nan,
                dtype=float,
            )
    for eta_index in range(len(eta_values)):
        for item in references_by_eta[int(eta_index)]:
            arrays[(item.rod_label, item.boundary_condition)][item.reference_mode_index - 1, eta_index] = (
                item.Lambda_reference
            )
    return arrays


def family_style(rod_label: str, boundary_condition: str) -> dict[str, object]:
    return single_mu.family_style(rod_label, boundary_condition)


def line_label(rod_label: str, boundary_condition: str) -> str:
    bc = "CP" if boundary_condition == "clamped_pinned" else "CC"
    return f"{rod_label} {bc} refs"


def plot_mu_case(
    path: Path,
    *,
    args: Args,
    mu: float,
    eta_values: np.ndarray,
    system_plot_grid: np.ndarray,
    reference_arrays: dict[tuple[str, str], np.ndarray],
    full_range: bool,
) -> dict[str, object]:
    fig, ax = plt.subplots(figsize=(10.8, 6.35))
    system_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for mode_index in range(args.n_system_roots):
        ax.plot(
            eta_values,
            system_plot_grid[:, mode_index],
            color=system_colors[mode_index % len(system_colors)],
            linewidth=2.15,
            linestyle="-",
            zorder=5,
        )

    short_curve_plotted = 0
    short_curve_omitted = 0
    short_curve_clipped = 0
    long_curve_plotted = 0
    y_min = 0.0 if full_range else float(args.readable_ymin)
    y_max = float(args.readable_ymax)

    for (rod_label, boundary_condition), values in reference_arrays.items():
        style = family_style(rod_label, boundary_condition)
        for ref_index in range(min(args.n_reference_roots, values.shape[0])):
            raw = np.asarray(values[ref_index], dtype=float)
            finite = raw[np.isfinite(raw)]
            if not finite.size:
                continue
            if full_range:
                should_plot = True
                plot_values = raw
            elif rod_label == "short":
                should_plot = bool(np.nanmin(finite) <= y_max and np.nanmax(finite) >= y_min)
                plot_values = np.where((raw >= y_min) & (raw <= y_max), raw, np.nan)
            else:
                should_plot = True
                plot_values = raw
            if not should_plot:
                short_curve_omitted += 1
                continue
            if rod_label == "short":
                short_curve_plotted += 1
                if np.any(np.isfinite(raw) & ((raw > y_max) | (raw < y_min))):
                    short_curve_clipped += 1
            else:
                long_curve_plotted += 1
            ax.plot(
                eta_values,
                plot_values,
                color=style["color"],
                linewidth=1.22,
                linestyle=style["linestyle"],
                alpha=float(style["alpha"]),
                zorder=2,
            )

    if full_range:
        finite_arrays = [system_plot_grid[np.isfinite(system_plot_grid)]]
        finite_arrays.extend(arr[np.isfinite(arr)] for arr in reference_arrays.values())
        finite_arrays = [arr for arr in finite_arrays if arr.size]
        y_upper = 1.06 * float(np.max(np.concatenate(finite_arrays))) if finite_arrays else 1.0
        ax.set_ylim(0.0, y_upper)
        scope = "full reference range"
    else:
        y_upper = y_max
        ax.set_ylim(y_min, y_upper)
        scope = f"readable plot clipped to {_fmt(y_min)} <= Lambda <= {_fmt(y_upper)}"

    ax.set_xlabel("eta")
    ax.set_ylabel("Lambda")
    ax.set_title(
        f"Lambda(eta), beta = {args.beta_deg:g} deg, epsilon = {args.epsilon:g}, mu = {mu:g}\n"
        f"solid: system sorted in-plane frequencies; dashed: single-rod references ({scope})"
    )
    ax.grid(True, color="0.88", linewidth=0.65)
    handles = [
        Line2D([0], [0], color="black", lw=2.15, ls="-", label=f"system sorted roots, first {args.n_system_roots}"),
    ]
    for rod_label, boundary_condition in (
        ("long", "clamped_pinned"),
        ("long", "clamped_clamped"),
        ("short", "clamped_pinned"),
        ("short", "clamped_clamped"),
    ):
        style = family_style(rod_label, boundary_condition)
        handles.append(
            Line2D(
                [0],
                [0],
                color=style["color"],
                lw=1.4,
                ls=style["linestyle"],
                alpha=style["alpha"],
                label=line_label(rod_label, boundary_condition),
            )
        )
    ax.legend(handles=handles, loc="upper left", fontsize=8.5, frameon=False, handlelength=3.6)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return {
        "y_visible_min": y_min if not full_range else 0.0,
        "y_visible_max": y_upper,
        "long_reference_curves_plotted": long_curve_plotted,
        "short_reference_curves_plotted": short_curve_plotted,
        "short_reference_curves_omitted": short_curve_omitted,
        "short_reference_curves_clipped": short_curve_clipped,
        "full_range": "yes" if full_range else "no",
    }


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> list[str]:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return lines


def write_mu_report(
    path: Path,
    *,
    args: Args,
    result: MuCaseResult,
) -> None:
    finite_system = result.system_grid[np.isfinite(result.system_grid)]
    missing = int(np.count_nonzero(~np.isfinite(result.system_grid)))
    warning_count = sum(1 for row in result.warnings for item in row if item != "none")
    stats = result.main_plot_stats
    step_text = "smoke grid" if args.smoke else _fmt(args.eta_step)
    lines = [
        f"# Lambda(eta) Single-Rod Reference Diagnostic, mu = {_fmt(result.mu)}",
        "",
        "This is a diagnostic-only analytic plot and nearest-frequency comparison.",
        "",
        "## Parameters",
        "",
        f"- beta_deg: {_fmt(args.beta_deg)}",
        f"- epsilon: {_fmt(args.epsilon)}",
        f"- mu: {_fmt(result.mu)}",
        f"- eta grid: {_fmt(float(result.eta_values[0]))} to {_fmt(float(result.eta_values[-1]))}, {len(result.eta_values)} points",
        f"- nominal eta step: {step_text}",
        f"- system sorted roots per eta: {args.n_system_roots}",
        f"- reference roots saved for all four families: {args.n_reference_roots}",
        f"- readable-ymin: {_fmt(args.readable_ymin)}",
        f"- readable-ymax: {_fmt(args.readable_ymax)}",
        "",
        "## Readable Plot Window",
        "",
        "- The readable PNG uses the fixed configured y-range; by default this is 0 <= Lambda <= 13.",
        "- Full-range PNGs keep the complete plotted reference range.",
        "- CSV data, nearest matches, and candidate tables are not clipped.",
        "- Long-rod CP and CC reference modes are plotted for every saved reference mode.",
        "- High short-rod reference curves are clipped or omitted from the readable view only.",
        f"- readable y_visible_min: {_fmt(stats['y_visible_min'])}",
        f"- readable y_visible_max: {_fmt(stats['y_visible_max'])}",
        f"- long-reference curves plotted in readable PNG: {stats['long_reference_curves_plotted']}",
        f"- short-reference curves plotted in readable PNG: {stats['short_reference_curves_plotted']}",
        f"- short-reference curves omitted from readable PNG: {stats['short_reference_curves_omitted']}",
        f"- short-reference curves clipped by readable y-window: {stats['short_reference_curves_clipped']}",
        "",
        "## Root And Spike Checks",
        "",
        f"- finite system roots: {int(finite_system.size)} / {result.system_grid.size}",
        f"- missing system roots: {missing}",
        f"- system root warning entries: {warning_count}",
        f"- suspicious spike audit rows: {len(result.spike_rows)}",
        "",
        "## Dominant Nearest Reference Families",
        "",
    ]
    lines.extend(
        markdown_table(
            ["sorted index", "dominant family", "fraction", "switches", "switch locations"],
            [
                [
                    row["system_sorted_index"],
                    row["dominant_reference_family"],
                    row["fraction_of_eta_points"],
                    row["number_of_family_switches"],
                    row["switch_eta_locations"] or "none",
                ]
                for row in result.branch_summary
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Interpretation Scope",
            "",
            "- This is frequency-proximity diagnostic bookkeeping only.",
            "- Family switches are nearest-reference changes for sorted frequencies, not descendant branch switches.",
            "- No crossing/no-crossing or strict positive-gap claim is made.",
            "",
            "## Output Files",
            "",
        ]
    )
    for output_path in result.output_paths:
        lines.append(f"- {output_path.name}")
    lines.extend(
        [
            "",
            "## Protected Scope",
            "",
            "No descendant tracking, FEM, 3D FEM, Gmsh, CalculiX, article workspace edits, old determinant edits, old solver edits, baseline-result edits, or analytic formula changes are part of this workflow.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_combined_report(path: Path, args: Args, results: Sequence[MuCaseResult]) -> None:
    lines = [
        "# Beta 90 Mu Scan Lambda(eta) Single-Rod Reference Diagnostic",
        "",
        "This report summarizes diagnostic-only sorted in-plane Lambda(eta) plots against isolated-rod references.",
        "",
        "## Parameters",
        "",
        f"- beta_deg: {_fmt(args.beta_deg)}",
        f"- epsilon: {_fmt(args.epsilon)}",
        f"- mu values: {', '.join(_fmt(result.mu) for result in results)}",
        f"- eta grid: {_fmt(float(results[0].eta_values[0]))} to {_fmt(float(results[0].eta_values[-1]))}, {len(results[0].eta_values)} points",
        f"- nominal eta step: {'smoke grid' if args.smoke else _fmt(args.eta_step)}",
        f"- system roots per eta: {args.n_system_roots}",
        f"- reference roots saved per family: {args.n_reference_roots}",
        f"- readable y-range: {_fmt(args.readable_ymin)} <= Lambda <= {_fmt(args.readable_ymax)}",
        "- full-range plots keep the complete plotted reference range",
        "- CSV data and matching tables are not clipped by the readable y-range",
        "",
        "## Combined Dominant Families",
        "",
    ]
    table_rows: list[list[object]] = []
    for result in results:
        for row in result.branch_summary:
            table_rows.append(
                [
                    _fmt(result.mu),
                    row["system_sorted_index"],
                    row["dominant_reference_family"],
                    row["fraction_of_eta_points"],
                    row["number_of_family_switches"],
                    row["switch_eta_locations"] or "none",
                ]
            )
    lines.extend(
        markdown_table(
            ["mu", "sorted index", "dominant family", "fraction", "switches", "switch locations"],
            table_rows,
        )
    )
    lines.extend(
        [
            "",
            "## Scope",
            "",
            "System curves are solid sorted in-plane Euler-Bernoulli frequencies. Reference curves are dashed isolated-rod CP/FP and CC/FF bending frequencies. The readable PNGs use the configured fixed y-range only; high short-rod reference curves may be clipped or omitted from that view, while CSV data and full-range PNGs remain complete. This mu scan uses no descendant tracking, FEM, 3D FEM, Gmsh, CalculiX, article artifacts, determinant edits, old-solver edits, baseline-result edits, or strict gap/crossing verification.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def output_paths_for(case_dir: Path, args: Args, mu: float) -> dict[str, Path]:
    prefix = file_prefix(args, mu)
    return {
        "system": case_dir / f"{prefix}_system_in_plane_sorted.csv",
        "references": case_dir / f"{prefix}_single_rod_references.csv",
        "matches": case_dir / f"{prefix}_nearest_reference_matches.csv",
        "all_candidates": case_dir / f"{prefix}_all_candidate_matches.csv",
        "branch_summary": case_dir / f"{prefix}_branch_family_summary.csv",
        "main_png": case_dir / f"{prefix}_in_plane_vs_single_rod_refs.png",
        "full_png": case_dir / f"{prefix}_in_plane_vs_single_rod_refs_full.png",
        "report": case_dir / f"{prefix}_report.md",
    }


def run_mu_case(args: Args, mu: float, eta_values: np.ndarray) -> MuCaseResult:
    case_dir = args.output_dir / mu_token(mu)
    references, references_by_eta = build_reference_values(args, mu, eta_values)
    system_grid, warnings, notes = solve_system_roots(args, mu, eta_values)
    system_plot_grid, spike_rows = audit_and_repair_spikes(args, mu, eta_values, system_grid, warnings, notes)
    matches, all_candidates = compare_candidates(args, mu, eta_values, system_grid, references_by_eta)
    branch_summary = branch_family_summary(args, mu, matches)
    reference_arrays = reference_curve_arrays(args, eta_values, references_by_eta)
    paths = output_paths_for(case_dir, args, mu)

    write_csv(paths["system"], system_rows(args, mu, eta_values, system_grid, warnings, notes), SYSTEM_FIELDS)
    write_csv(paths["references"], reference_rows(args, references), REFERENCE_FIELDS)
    write_csv(paths["matches"], matches, MATCH_FIELDS)
    write_csv(paths["all_candidates"], all_candidates, ALL_CANDIDATE_FIELDS)
    write_csv(paths["branch_summary"], branch_summary, BRANCH_SUMMARY_FIELDS)
    main_stats = plot_mu_case(
        paths["main_png"],
        args=args,
        mu=mu,
        eta_values=eta_values,
        system_plot_grid=system_plot_grid,
        reference_arrays=reference_arrays,
        full_range=False,
    )
    full_stats = plot_mu_case(
        paths["full_png"],
        args=args,
        mu=mu,
        eta_values=eta_values,
        system_plot_grid=system_plot_grid,
        reference_arrays=reference_arrays,
        full_range=True,
    )
    output_paths = [
        paths["system"],
        paths["references"],
        paths["matches"],
        paths["all_candidates"],
        paths["branch_summary"],
        paths["main_png"],
        paths["full_png"],
        paths["report"],
    ]
    pending = MuCaseResult(
        mu=float(mu),
        mu_token=mu_token(mu),
        output_dir=case_dir,
        eta_values=eta_values,
        system_grid=system_grid,
        system_plot_grid=system_plot_grid,
        warnings=warnings,
        spike_rows=spike_rows,
        branch_summary=branch_summary,
        main_plot_stats=main_stats,
        full_plot_stats=full_stats,
        output_paths=output_paths,
    )
    write_mu_report(paths["report"], args=args, result=pending)
    return pending


def combined_summary_rows(results: Sequence[MuCaseResult]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        for row in result.branch_summary:
            rows.append(
                {
                    "mu": _fmt(result.mu),
                    "system_sorted_index": row["system_sorted_index"],
                    "dominant_reference_family": row["dominant_reference_family"],
                    "fraction_of_eta_points": row["fraction_of_eta_points"],
                    "number_of_family_switches": row["number_of_family_switches"],
                    "switch_eta_locations": row["switch_eta_locations"],
                    "notes": row["notes"],
                }
            )
    return rows


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        description="Plot diagnostic Lambda(eta) sorted in-plane roots at beta=90 with single-rod references."
    )
    parser.add_argument("--beta-deg", type=float, default=DEFAULT_BETA_DEG)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--mu-values", type=float, nargs="+", default=list(DEFAULT_MU_VALUES))
    parser.add_argument("--eta-min", type=float, default=DEFAULT_ETA_MIN)
    parser.add_argument("--eta-max", type=float, default=DEFAULT_ETA_MAX)
    parser.add_argument("--eta-step", type=float, default=DEFAULT_ETA_STEP)
    parser.add_argument("--n-system-roots", type=int, default=DEFAULT_N_SYSTEM_ROOTS)
    parser.add_argument("--n-reference-roots", type=int, default=DEFAULT_N_REFERENCE_ROOTS)
    parser.add_argument("--readable-ymin", type=float, default=DEFAULT_READABLE_YMIN)
    parser.add_argument("--readable-ymax", type=float, default=DEFAULT_READABLE_YMAX)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(argv)
    mu_values = tuple(float(value) for value in ns.mu_values)
    n_system_roots = int(ns.n_system_roots)
    n_reference_roots = int(ns.n_reference_roots)
    out_dir = output_dir(Path(ns.output_dir))
    if bool(ns.smoke):
        n_system_roots = min(n_system_roots, 4)
        n_reference_roots = min(n_reference_roots, 4)
        if Path(ns.output_dir) == DEFAULT_OUTPUT_DIR:
            out_dir = SMOKE_OUTPUT_DIR
    args = Args(
        beta_deg=float(ns.beta_deg),
        epsilon=float(ns.epsilon),
        mu_values=mu_values,
        eta_min=float(ns.eta_min),
        eta_max=float(ns.eta_max),
        eta_step=float(ns.eta_step),
        n_system_roots=n_system_roots,
        n_reference_roots=n_reference_roots,
        readable_ymin=float(ns.readable_ymin),
        readable_ymax=float(ns.readable_ymax),
        output_dir=out_dir,
        smoke=bool(ns.smoke),
    )
    validate_args(args, eta_grid(args))
    return args


def run(args: Args) -> dict[str, object]:
    values = eta_grid(args)
    validate_args(args, values)
    results = [run_mu_case(args, mu, values) for mu in args.mu_values]
    combined_csv = args.output_dir / COMBINED_SUMMARY_NAME
    combined_report = args.output_dir / COMBINED_REPORT_NAME
    write_csv(combined_csv, combined_summary_rows(results), COMBINED_SUMMARY_FIELDS)
    write_combined_report(combined_report, args, results)
    return {
        "eta_values": values,
        "results": results,
        "combined_csv": combined_csv,
        "combined_report": combined_report,
    }


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    result = run(args)
    eta_values = np.asarray(result["eta_values"], dtype=float)
    results: Sequence[MuCaseResult] = result["results"]
    print("diagnostic-only Lambda(eta) mu scan with single-rod references")
    print(f"beta_deg={args.beta_deg:g}, epsilon={args.epsilon:g}")
    print(
        f"eta grid: {len(eta_values)} points from {float(eta_values[0]):g} to {float(eta_values[-1]):g}; "
        f"step={'smoke grid' if args.smoke else f'{args.eta_step:g}'}"
    )
    print(f"system roots per eta={args.n_system_roots}, reference roots saved={args.n_reference_roots}")
    print(f"readable y-range={args.readable_ymin:g}..{args.readable_ymax:g}")
    for case in results:
        missing = int(np.count_nonzero(~np.isfinite(case.system_grid)))
        warning_count = sum(1 for row in case.warnings for item in row if item != "none")
        print(
            f"mu={case.mu:g}: missing roots={missing}; root warnings={warning_count}; "
            f"spike rows={len(case.spike_rows)}"
        )
        print(
            "  readable reference handling: "
            f"long curves plotted={case.main_plot_stats['long_reference_curves_plotted']}; "
            f"short curves plotted={case.main_plot_stats['short_reference_curves_plotted']}; "
            f"short curves omitted={case.main_plot_stats['short_reference_curves_omitted']}; "
            f"short curves clipped={case.main_plot_stats['short_reference_curves_clipped']}"
        )
        for row in case.branch_summary:
            print(
                f"  sorted {row['system_sorted_index']}: {row['dominant_reference_family']} "
                f"(fraction={row['fraction_of_eta_points']}, switches={row['number_of_family_switches']})"
            )
        for output_path in case.output_paths:
            print(f"  saved: {output_path}")
    print(f"combined summary: {result['combined_csv']}")
    print(f"combined report: {result['combined_report']}")
    print("solid curves: sorted in-plane system roots; dashed curves: single-rod references")
    print("no descendant tracking; no FEM/Gmsh/CalculiX; no article/formula/baseline changes")
    return result


if __name__ == "__main__":
    main()
