from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
DEFAULT_MU = 0.3
DEFAULT_ETA = 0.1
DEFAULT_N_SYSTEM_ROOTS = 6
DEFAULT_N_REFERENCE_ROOTS = 10
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "beta90_in_plane_single_rod_reference_comparison"
SMOKE_OUTPUT_DIR = REPO_ROOT / "results" / "_smoke" / "beta90_in_plane_single_rod_reference_comparison"

SUMMARY_CSV_NAME = "beta90_in_plane_first6_reference_matches.csv"
ALL_CANDIDATES_CSV_NAME = "beta90_in_plane_first6_reference_all_candidates.csv"
REFERENCE_CSV_NAME = "beta90_single_rod_reference_frequencies.csv"
REPORT_NAME = "beta90_in_plane_first6_reference_matches_report.md"
PLOT_NAME = "beta90_in_plane_first6_vs_single_rod_references.png"

AMBIGUOUS_ABS_TOL = 1.0e-3
AMBIGUOUS_REL_TOL = 1.0e-3

SUMMARY_FIELDS = [
    "system_sorted_index",
    "Lambda_system",
    "best_rod_id",
    "best_rod_label",
    "best_length_factor",
    "best_tau",
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
    "system_sorted_index",
    "Lambda_system",
    "candidate_rank_by_abs_diff",
    "rod_id",
    "rod_label",
    "length_factor",
    "tau",
    "boundary_condition",
    "reference_mode_index",
    "alpha_root",
    "Lambda_reference",
    "abs_diff",
    "rel_diff_system",
    "rel_diff_reference",
    "notes",
]

REFERENCE_FIELDS = [
    "rod_id",
    "rod_label",
    "length_factor",
    "tau",
    "boundary_condition",
    "reference_mode_index",
    "alpha_root",
    "Lambda_reference",
]


@dataclass(frozen=True)
class AuditParameters:
    beta_deg: float
    beta_rad: float
    epsilon: float
    mu: float
    eta: float
    n_system_roots: int
    n_reference_roots: int
    output_dir: Path


@dataclass(frozen=True)
class ReferenceCandidate:
    rod_id: int
    rod_label: str
    length_factor: float
    tau: float
    boundary_condition: str
    reference_mode_index: int
    alpha_root: float
    Lambda_reference: float


@dataclass(frozen=True)
class SystemRoot:
    sorted_index: int
    Lambda_system: float
    root_solver_warning: str
    notes: str


def _fmt(value: object) -> str:
    try:
        value_f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not isfinite(value_f):
        return "nan"
    return f"{value_f:.16g}"


def _resolve_output_dir(path: Path) -> Path:
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def _write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _validate_args(params: AuditParameters) -> None:
    if not isfinite(params.beta_deg) or not isfinite(params.beta_rad):
        raise ValueError("beta must be finite.")
    if not (-1.0 < params.mu < 1.0):
        raise ValueError("mu must lie inside (-1, 1).")
    if params.epsilon <= 0.0 or not isfinite(params.epsilon):
        raise ValueError("epsilon must be positive and finite.")
    if params.n_system_roots <= 0:
        raise ValueError("n-system-roots must be positive.")
    if params.n_reference_roots <= 0:
        raise ValueError("n-reference-roots must be positive.")
    thickness_mismatch_factors(params.mu, params.eta)


def _rod_label(length_factor: float, other_length_factor: float) -> str:
    return "short" if float(length_factor) <= float(other_length_factor) else "long"


def _rod_states(mu: float, eta: float) -> dict[int, tuple[str, float, float]]:
    factors = thickness_mismatch_factors(float(mu), float(eta))
    l1 = 1.0 - float(mu)
    l2 = 1.0 + float(mu)
    return {
        1: (_rod_label(l1, l2), l1, factors.tau1),
        2: (_rod_label(l2, l1), l2, factors.tau2),
    }


def _validate_lengths(states: dict[int, tuple[str, float, float]]) -> None:
    lengths = [state[1] for state in states.values()]
    if not all(isfinite(length) and length > 0.0 for length in lengths):
        raise ValueError("single-rod length factors must be finite and positive.")
    if not (min(lengths) < max(lengths)):
        raise ValueError("this diagnostic expects one short and one long rod; choose nonzero mu.")
    short = [rod_id for rod_id, (label, _, _) in states.items() if label == "short"]
    long = [rod_id for rod_id, (label, _, _) in states.items() if label == "long"]
    if len(short) != 1 or len(long) != 1:
        raise ValueError("rod labels must identify exactly one short and one long rod.")


def _increasing(values: Sequence[float]) -> bool:
    arr = np.asarray(values, dtype=float)
    return bool(np.all(np.isfinite(arr)) and np.all(np.diff(arr) > 0.0))


def _alpha_roots(n_reference_roots: int) -> dict[str, np.ndarray]:
    return {
        "clamped_pinned": np.asarray(roots_clamped_supported(int(n_reference_roots)), dtype=float),
        "clamped_clamped": np.asarray(fixed_fixed_lambdas(int(n_reference_roots)), dtype=float),
    }


def build_reference_candidates(params: AuditParameters) -> list[ReferenceCandidate]:
    states = _rod_states(params.mu, params.eta)
    _validate_lengths(states)
    alpha_roots = _alpha_roots(params.n_reference_roots)
    candidates: list[ReferenceCandidate] = []
    for boundary_condition, alphas in alpha_roots.items():
        if len(alphas) < params.n_reference_roots or not _increasing(alphas):
            raise RuntimeError(f"{boundary_condition} alpha roots are not a finite increasing sequence.")
        for rod_id in (1, 2):
            rod_label, length_factor, tau = states[rod_id]
            if not (isfinite(tau) and tau > 0.0):
                raise ValueError(f"tau{rod_id} must be finite and positive.")
            lambdas = np.asarray(alphas, dtype=float) * np.sqrt(tau) / length_factor
            if not _increasing(lambdas):
                raise RuntimeError(
                    f"reference frequencies are not increasing for rod {rod_id}, {boundary_condition}."
                )
            for index, (alpha, value) in enumerate(zip(alphas, lambdas), start=1):
                candidates.append(
                    ReferenceCandidate(
                        rod_id=int(rod_id),
                        rod_label=rod_label,
                        length_factor=float(length_factor),
                        tau=float(tau),
                        boundary_condition=boundary_condition,
                        reference_mode_index=int(index),
                        alpha_root=float(alpha),
                        Lambda_reference=float(value),
                    )
                )
    return candidates


def solve_system_roots(params: AuditParameters) -> list[SystemRoot]:
    roots = np.asarray(
        find_first_n_roots_eta(
            params.beta_rad,
            params.mu,
            params.epsilon,
            params.eta,
            params.n_system_roots,
        ),
        dtype=float,
    )
    rows: list[SystemRoot] = []
    for index in range(1, params.n_system_roots + 1):
        value = float(roots[index - 1]) if index - 1 < roots.size else float("nan")
        if isfinite(value) and value > 0.0:
            warning = "none"
            notes = "sorted in-plane Euler-Bernoulli root"
        else:
            warning = "missing_or_nonfinite_root"
            notes = "missing root from find_first_n_roots_eta"
        rows.append(
            SystemRoot(
                sorted_index=index,
                Lambda_system=value,
                root_solver_warning=warning,
                notes=notes,
            )
        )
    finite_values = [root.Lambda_system for root in rows if isfinite(root.Lambda_system)]
    if len(finite_values) >= 2 and any(b <= a for a, b in zip(finite_values, finite_values[1:])):
        raise RuntimeError("system roots are not strictly increasing where finite.")
    return rows


def _candidate_comparison(system: SystemRoot, candidate: ReferenceCandidate) -> dict[str, object]:
    system_value = system.Lambda_system
    reference_value = candidate.Lambda_reference
    if isfinite(system_value) and isfinite(reference_value) and system_value > 0.0 and reference_value > 0.0:
        abs_diff = abs(system_value - reference_value)
        rel_diff_system = abs_diff / system_value
        rel_diff_reference = abs_diff / reference_value
        notes = "frequency-proximity diagnostic only"
    else:
        abs_diff = float("nan")
        rel_diff_system = float("nan")
        rel_diff_reference = float("nan")
        notes = "comparison unavailable because a frequency is nonfinite"
    return {
        "system_sorted_index": system.sorted_index,
        "Lambda_system": _fmt(system_value),
        "rod_id": candidate.rod_id,
        "rod_label": candidate.rod_label,
        "length_factor": _fmt(candidate.length_factor),
        "tau": _fmt(candidate.tau),
        "boundary_condition": candidate.boundary_condition,
        "reference_mode_index": candidate.reference_mode_index,
        "alpha_root": _fmt(candidate.alpha_root),
        "Lambda_reference": _fmt(reference_value),
        "abs_diff": _fmt(abs_diff),
        "rel_diff_system": _fmt(rel_diff_system),
        "rel_diff_reference": _fmt(rel_diff_reference),
        "notes": notes,
        "_abs_diff_value": abs_diff,
        "_rel_diff_system_value": rel_diff_system,
    }


def build_candidate_rows(
    system_roots: Sequence[SystemRoot],
    candidates: Sequence[ReferenceCandidate],
) -> dict[int, list[dict[str, object]]]:
    rows_by_system: dict[int, list[dict[str, object]]] = {}
    for system in system_roots:
        rows = [_candidate_comparison(system, candidate) for candidate in candidates]
        rows.sort(
            key=lambda row: (
                float(row["_abs_diff_value"]) if isfinite(float(row["_abs_diff_value"])) else float("inf"),
                str(row["rod_label"]),
                str(row["boundary_condition"]),
                int(row["reference_mode_index"]),
            )
        )
        for rank, row in enumerate(rows, start=1):
            row["candidate_rank_by_abs_diff"] = rank
        rows_by_system[int(system.sorted_index)] = rows
    return rows_by_system


def _ambiguous_candidates(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    finite_rows = [row for row in rows if isfinite(float(row["_abs_diff_value"]))]
    if len(finite_rows) < 2:
        return []
    best = finite_rows[0]
    best_abs = float(best["_abs_diff_value"])
    best_rel = float(best["_rel_diff_system_value"])
    ambiguous = [best]
    for row in finite_rows[1:]:
        abs_gap = abs(float(row["_abs_diff_value"]) - best_abs)
        rel_gap = abs(float(row["_rel_diff_system_value"]) - best_rel)
        if abs_gap < AMBIGUOUS_ABS_TOL or rel_gap < AMBIGUOUS_REL_TOL:
            ambiguous.append(row)
        else:
            break
    return ambiguous if len(ambiguous) > 1 else []


def _candidate_label(row: dict[str, object]) -> str:
    return (
        f"{row['rod_label']} rod, {row['boundary_condition']}, "
        f"mode {row['reference_mode_index']}"
    )


def build_summary_rows(
    system_roots: Sequence[SystemRoot],
    rows_by_system: dict[int, list[dict[str, object]]],
) -> list[dict[str, object]]:
    summary: list[dict[str, object]] = []
    for system in system_roots:
        rows = rows_by_system[int(system.sorted_index)]
        finite_rows = [row for row in rows if isfinite(float(row["_abs_diff_value"]))]
        if not finite_rows:
            summary.append(
                {
                    "system_sorted_index": system.sorted_index,
                    "Lambda_system": _fmt(system.Lambda_system),
                    "best_rod_id": "",
                    "best_rod_label": "",
                    "best_length_factor": "",
                    "best_tau": "",
                    "best_boundary_condition": "",
                    "best_reference_mode_index": "",
                    "Lambda_reference": "nan",
                    "abs_diff": "nan",
                    "rel_diff_system": "nan",
                    "rel_diff_reference": "nan",
                    "ambiguous": "yes",
                    "notes": f"{system.root_solver_warning}; no finite comparison candidates",
                }
            )
            continue
        best = finite_rows[0]
        ambiguous = _ambiguous_candidates(finite_rows)
        if ambiguous:
            ambiguous_labels = ", ".join(_candidate_label(row) for row in ambiguous[:4])
            ambiguous_text = f"ambiguous with {ambiguous_labels}"
            ambiguous_flag = "yes"
        else:
            ambiguous_text = "unique closest candidate by absolute difference"
            ambiguous_flag = "no"
        notes = f"{system.notes}; root_solver_warning={system.root_solver_warning}; {ambiguous_text}"
        summary.append(
            {
                "system_sorted_index": system.sorted_index,
                "Lambda_system": best["Lambda_system"],
                "best_rod_id": best["rod_id"],
                "best_rod_label": best["rod_label"],
                "best_length_factor": best["length_factor"],
                "best_tau": best["tau"],
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
    return summary


def clean_candidate_rows(rows_by_system: dict[int, list[dict[str, object]]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for system_index in sorted(rows_by_system):
        for row in rows_by_system[system_index]:
            clean = dict(row)
            clean.pop("_abs_diff_value", None)
            clean.pop("_rel_diff_system_value", None)
            out.append(clean)
    return out


def reference_rows(candidates: Sequence[ReferenceCandidate]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    ordered = sorted(
        candidates,
        key=lambda item: (
            0 if item.rod_label == "long" else 1,
            item.boundary_condition,
            item.reference_mode_index,
        ),
    )
    for candidate in ordered:
        rows.append(
            {
                "rod_id": candidate.rod_id,
                "rod_label": candidate.rod_label,
                "length_factor": _fmt(candidate.length_factor),
                "tau": _fmt(candidate.tau),
                "boundary_condition": candidate.boundary_condition,
                "reference_mode_index": candidate.reference_mode_index,
                "alpha_root": _fmt(candidate.alpha_root),
                "Lambda_reference": _fmt(candidate.Lambda_reference),
            }
        )
    return rows


def _markdown_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> list[str]:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return lines


def write_report(
    path: Path,
    params: AuditParameters,
    system_roots: Sequence[SystemRoot],
    candidates: Sequence[ReferenceCandidate],
    rows_by_system: dict[int, list[dict[str, object]]],
    summary_rows: Sequence[dict[str, object]],
    plot_path: Path,
) -> None:
    states = _rod_states(params.mu, params.eta)
    factors = thickness_mismatch_factors(params.mu, params.eta)
    root_warnings = [root for root in system_roots if root.root_solver_warning != "none"]
    ambiguous_rows = [row for row in summary_rows if str(row["ambiguous"]) == "yes"]
    lines: list[str] = [
        "# Beta 90 In-Plane Single-Rod Reference Comparison",
        "",
        "This is a diagnostic-only sorted-frequency comparison. It is not an article figure, not a FEM run, and not descendant branch tracking.",
        "",
        "## Parameters",
        "",
        f"- beta_deg: {_fmt(params.beta_deg)}",
        f"- beta_rad: {_fmt(params.beta_rad)}",
        f"- epsilon: {_fmt(params.epsilon)}",
        f"- mu: {_fmt(params.mu)}",
        f"- eta: {_fmt(params.eta)}",
        f"- n_system_roots: {params.n_system_roots}",
        f"- n_reference_roots: {params.n_reference_roots}",
        "",
        "## Length And Thickness Factors",
        "",
        f"- L1 = 1 - mu = {_fmt(states[1][1])}; rod 1 label: {states[1][0]}",
        f"- L2 = 1 + mu = {_fmt(states[2][1])}; rod 2 label: {states[2][0]}",
        f"- tau1 = {_fmt(factors.tau1)}",
        f"- tau2 = {_fmt(factors.tau2)}",
        f"- mass factor = {_fmt(factors.mass_factor)}",
        "",
        "## System Frequencies",
        "",
    ]
    lines.extend(
        _markdown_table(
            ["sorted index", "Lambda_system", "root solver warning", "notes"],
            [
                [
                    root.sorted_index,
                    _fmt(root.Lambda_system),
                    root.root_solver_warning,
                    root.notes,
                ]
                for root in system_roots
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Reference Formula",
            "",
            "Reference alpha roots are isolated Euler-Bernoulli bending roots.",
            "",
            "- clamped_pinned: tan(alpha) = tanh(alpha)",
            "- clamped_clamped: cosh(alpha) cos(alpha) = 1",
            "- Lambda_ref_i,n = alpha_n * sqrt(tau_i) / L_i",
            "",
            "Axial frequencies are not included.",
            "",
            "## Best Matches",
            "",
        ]
    )
    lines.extend(
        _markdown_table(
            [
                "system mode",
                "Lambda_system",
                "closest family",
                "ref mode",
                "Lambda_ref",
                "abs diff",
                "rel diff system",
                "ambiguous",
            ],
            [
                [
                    row["system_sorted_index"],
                    row["Lambda_system"],
                    f"{row['best_rod_label']} {row['best_boundary_condition']}",
                    row["best_reference_mode_index"],
                    row["Lambda_reference"],
                    row["abs_diff"],
                    row["rel_diff_system"],
                    row["ambiguous"],
                ]
                for row in summary_rows
            ],
        )
    )
    lines.extend(["", "## Top 3 Candidates Per System Mode", ""])
    for system in system_roots:
        lines.append(f"### System Mode {system.sorted_index}")
        top_rows = rows_by_system[int(system.sorted_index)][:3]
        lines.extend(
            _markdown_table(
                ["rank", "family", "ref mode", "Lambda_ref", "abs diff", "rel diff system"],
                [
                    [
                        row["candidate_rank_by_abs_diff"],
                        _candidate_label(row),
                        row["reference_mode_index"],
                        row["Lambda_reference"],
                        row["abs_diff"],
                        row["rel_diff_system"],
                    ]
                    for row in top_rows
                ],
            )
        )
        lines.append("")
    lines.extend(
        [
            "## Interpretation",
            "",
            "Closest references are selected by minimum absolute difference in Lambda.",
            "This is a frequency-proximity diagnostic only, not a proof that the system mode is physically identical to an isolated-rod mode.",
            "The comparison uses analytic in-plane sorted frequencies only, not descendant branches.",
            "No crossing or no-crossing claim is made here.",
            "",
            "## Warnings",
            "",
        ]
    )
    if root_warnings:
        for root in root_warnings:
            lines.append(f"- system mode {root.sorted_index}: {root.root_solver_warning}")
    else:
        lines.append("- no system root warnings")
    if ambiguous_rows:
        for row in ambiguous_rows:
            lines.append(f"- system mode {row['system_sorted_index']}: {row['notes']}")
    else:
        lines.append("- no ambiguity warnings by the configured proximity thresholds")
    lines.extend(
        [
            "",
            "## Output Files",
            "",
            f"- {SUMMARY_CSV_NAME}",
            f"- {ALL_CANDIDATES_CSV_NAME}",
            f"- {REFERENCE_CSV_NAME}",
            f"- {REPORT_NAME}",
            f"- {plot_path.name}",
            "",
            "## Protected Scope",
            "",
            "No FEM, Gmsh, CalculiX, article workspace, old determinant, old solvers, baseline results, or analytic formulas are modified by this audit.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_plot(
    path: Path,
    params: AuditParameters,
    system_roots: Sequence[SystemRoot],
    candidates: Sequence[ReferenceCandidate],
    rows_by_system: dict[int, list[dict[str, object]]],
) -> None:
    visible_reference_keys: set[tuple[int, str, int]] = set()
    for rows in rows_by_system.values():
        for row in rows[:3]:
            visible_reference_keys.add(
                (
                    int(row["rod_id"]),
                    str(row["boundary_condition"]),
                    int(row["reference_mode_index"]),
                )
            )
    visible_candidates = [
        candidate
        for candidate in candidates
        if (candidate.rod_id, candidate.boundary_condition, candidate.reference_mode_index) in visible_reference_keys
    ]
    if not visible_candidates:
        visible_candidates = list(candidates)

    fig, ax = plt.subplots(figsize=(10.0, 5.8))
    family_order = [
        ("system", "system"),
        ("long clamped_pinned", "long\nCP"),
        ("long clamped_clamped", "long\nCC"),
        ("short clamped_pinned", "short\nCP"),
        ("short clamped_clamped", "short\nCC"),
    ]
    x_positions = {key: idx for idx, (key, _) in enumerate(family_order)}
    colors = {
        "system": "#111111",
        "clamped_pinned": "#1f77b4",
        "clamped_clamped": "#d62728",
    }

    system_values = [root.Lambda_system for root in system_roots if isfinite(root.Lambda_system)]
    ax.scatter(
        [x_positions["system"]] * len(system_values),
        system_values,
        s=54,
        color=colors["system"],
        zorder=5,
        label="system sorted roots",
    )
    for root in system_roots:
        if isfinite(root.Lambda_system):
            ax.annotate(
                str(root.sorted_index),
                (x_positions["system"], root.Lambda_system),
                xytext=(6, 2),
                textcoords="offset points",
                fontsize=8.5,
            )

    for candidate in visible_candidates:
        key = f"{candidate.rod_label} {candidate.boundary_condition}"
        x = x_positions[key]
        ax.scatter(
            [x],
            [candidate.Lambda_reference],
            s=38,
            marker="_",
            linewidths=2.0,
            color=colors[candidate.boundary_condition],
            zorder=3,
        )
        ax.annotate(
            str(candidate.reference_mode_index),
            (x, candidate.Lambda_reference),
            xytext=(5, -3),
            textcoords="offset points",
            fontsize=7.5,
            color="0.25",
        )

    y_values = system_values + [candidate.Lambda_reference for candidate in visible_candidates]
    y_max = max(y_values) * 1.12 if y_values else 1.0
    ax.set_xlim(-0.55, len(family_order) - 0.45)
    ax.set_ylim(0.0, y_max)
    ax.set_xticks([idx for idx, _ in enumerate(family_order)])
    ax.set_xticklabels([label for _, label in family_order])
    ax.set_ylabel("Lambda")
    ax.set_title(
        "Diagnostic beta=90 sorted in-plane roots vs nearest isolated-rod references\n"
        f"epsilon={params.epsilon:g}, mu={params.mu:g}, eta={params.eta:g}; reference markers show top-3 candidates per system root"
    )
    ax.grid(True, axis="y", color="0.88", linewidth=0.7)
    ax.grid(False, axis="x")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnostic-only sorted in-plane beta=90 frequency comparison against "
            "isolated single-rod bending references."
        )
    )
    parser.add_argument("--beta-deg", type=float, default=DEFAULT_BETA_DEG)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--mu", type=float, default=DEFAULT_MU)
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--n-system-roots", type=int, default=DEFAULT_N_SYSTEM_ROOTS)
    parser.add_argument("--n-reference-roots", type=int, default=DEFAULT_N_REFERENCE_ROOTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Run a tiny wiring check and write outputs under results/_smoke/ by default.",
    )
    return parser.parse_args(argv)


def params_from_args(args: argparse.Namespace) -> AuditParameters:
    n_system_roots = int(args.n_system_roots)
    n_reference_roots = int(args.n_reference_roots)
    output_dir = Path(args.output_dir)
    if bool(args.smoke):
        n_system_roots = 3
        n_reference_roots = 4
        if output_dir == DEFAULT_OUTPUT_DIR:
            output_dir = SMOKE_OUTPUT_DIR
    beta_deg = float(args.beta_deg)
    params = AuditParameters(
        beta_deg=beta_deg,
        beta_rad=float(np.deg2rad(beta_deg)),
        epsilon=float(args.epsilon),
        mu=float(args.mu),
        eta=float(args.eta),
        n_system_roots=n_system_roots,
        n_reference_roots=n_reference_roots,
        output_dir=_resolve_output_dir(output_dir),
    )
    _validate_args(params)
    return params


def run(params: AuditParameters) -> dict[str, object]:
    candidates = build_reference_candidates(params)
    system_roots = solve_system_roots(params)
    rows_by_system = build_candidate_rows(system_roots, candidates)
    summary = build_summary_rows(system_roots, rows_by_system)

    output_dir = params.output_dir
    summary_csv = output_dir / SUMMARY_CSV_NAME
    all_candidates_csv = output_dir / ALL_CANDIDATES_CSV_NAME
    reference_csv = output_dir / REFERENCE_CSV_NAME
    report_md = output_dir / REPORT_NAME
    plot_png = output_dir / PLOT_NAME

    _write_csv(summary_csv, summary, SUMMARY_FIELDS)
    _write_csv(all_candidates_csv, clean_candidate_rows(rows_by_system), ALL_CANDIDATE_FIELDS)
    _write_csv(reference_csv, reference_rows(candidates), REFERENCE_FIELDS)
    write_plot(plot_png, params, system_roots, candidates, rows_by_system)
    write_report(report_md, params, system_roots, candidates, rows_by_system, summary, plot_png)

    return {
        "system_roots": system_roots,
        "summary_rows": summary,
        "candidates": candidates,
        "summary_csv": summary_csv,
        "all_candidates_csv": all_candidates_csv,
        "reference_csv": reference_csv,
        "report_md": report_md,
        "plot_png": plot_png,
    }


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    params = params_from_args(args)
    result = run(params)
    factors = thickness_mismatch_factors(params.mu, params.eta)
    l1 = 1.0 - params.mu
    l2 = 1.0 + params.mu
    print("diagnostic-only analytic in-plane sorted-frequency comparison")
    print(f"beta_deg={params.beta_deg:g}, epsilon={params.epsilon:g}, mu={params.mu:g}, eta={params.eta:g}")
    print(f"L1={l1:.16g}, L2={l2:.16g}, tau1={factors.tau1:.16g}, tau2={factors.tau2:.16g}")
    print("system roots:")
    for root in result["system_roots"]:
        print(
            f"  sorted {root.sorted_index}: Lambda={root.Lambda_system:.16g}, "
            f"root_solver_warning={root.root_solver_warning}"
        )
    print("best matches:")
    for row in result["summary_rows"]:
        print(
            "  system mode "
            f"{row['system_sorted_index']}: closest to {row['best_rod_label']} rod, "
            f"{row['best_boundary_condition']}, mode {row['best_reference_mode_index']} "
            f"(Lambda_ref={row['Lambda_reference']}, abs_diff={row['abs_diff']}, "
            f"ambiguous={row['ambiguous']})"
        )
    print(f"saved summary CSV: {result['summary_csv']}")
    print(f"saved all-candidates CSV: {result['all_candidates_csv']}")
    print(f"saved reference CSV: {result['reference_csv']}")
    print(f"saved report: {result['report_md']}")
    print(f"saved plot: {result['plot_png']}")
    print("analytic in-plane sorted roots only; no FEM/Gmsh/CalculiX; no article or formula changes")
    return result


if __name__ == "__main__":
    main()
