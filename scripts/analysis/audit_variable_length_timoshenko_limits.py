from __future__ import annotations

import csv
import math
from pathlib import Path
import sys
from typing import Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.solvers import find_first_n_roots  # noqa: E402
from scripts.analysis import compare_coupled_equal_rods_eb_timoshenko as equal_timo  # noqa: E402
from scripts.lib.variable_length_timoshenko import timo_sorted_roots  # noqa: E402


N_MODES = 6
TEST_A_EPSILONS = [0.0025, 0.01, 0.025, 0.05]
TEST_B_BETAS = [0.0, 15.0, 45.0]
TEST_B_MUS = [0.0, 0.3, 0.6, 0.75]
TEST_B_EPSILONS = [0.02, 0.01, 0.005, 0.0025]
TEST_C_EPSILONS = [0.0025, 0.01, 0.025]
TEST_C_MUS = [0.0, 0.3, 0.6, 0.75, 0.9]

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_CSV = OUTPUT_DIR / "variable_length_timoshenko_limits_audit.csv"
OUTPUT_REPORT = OUTPUT_DIR / "variable_length_timoshenko_limits_audit.md"

A_ABS_TOL = 1.0e-7
A_REL_TOL = 1.0e-8
C_EB_ABS_TOL = 1.0e-8
C_EB_REL_TOL = 1.0e-8
C_TIMO_ABS_TOL = 1.0e-6
C_TIMO_REL_TOL = 1.0e-7

CSV_FIELDNAMES = [
    "test_name",
    "beta",
    "mu",
    "epsilon",
    "mode",
    "Lambda_EB",
    "Lambda_Timo_variable",
    "Lambda_Timo_reference",
    "Lambda_reference",
    "abs_diff",
    "rel_diff",
    "status",
    "notes",
]

FORMULA_AUDIT_LINES = [
    "Scope: eta=0 only; tau scaling for eta!=0 is not verified by this audit.",
    "Unknown ordering in the variable-length Timoshenko matrix: (A1, B1, P1, A2, B2, P2).",
    "Rod lengths: l1 = 1 - mu, l2 = 1 + mu.",
    "Endpoint arguments: rod 1 uses x = +l1; rod 2 uses x = -l2.",
    "theta = epsilon*Lambda^2.",
    "w_i(x) = A_i*(cosh(a*x) - cos(b*x)) + B_i*(sinh(a*x) - (h/q)*sin(b*x)).",
    "u_i(x) = P_i*sin(theta*x).",
    "psi_i(x) = A_i*(h*sinh(a*x) + q*sin(b*x)) + B_i*h*(cosh(a*x) - cos(b*x)).",
    "Displacement rows: w1 = w2*cos(beta) - u2*sin(beta).",
    "Displacement rows: u1 = w2*sin(beta) + u2*cos(beta).",
    "Rotation row: psi1 - psi2 = 0.",
    "Moment row: M1 - M2 = 0.",
    "Transverse-force row: -Q1 + Q2*cos(beta) - N2*sin(beta) = 0.",
    "Axial-force row: N1 - Q2*sin(beta) - N2*cos(beta) = 0.",
    "Constitutive definitions: M = E*I*psi'.",
    "Constitutive definitions: Q = kappa*G*A*(w' - psi).",
    "Constitutive definitions: N = E*A*u'.",
]


def finite_abs_rel(value: float, reference: float) -> tuple[float, float]:
    if not (math.isfinite(float(value)) and math.isfinite(float(reference))):
        return float("nan"), float("nan")
    abs_diff = abs(float(value) - float(reference))
    rel_diff = abs_diff / abs(float(reference)) if abs(float(reference)) > 0.0 else float("nan")
    return abs_diff, rel_diff


def status_from_tolerance(abs_diff: float, rel_diff: float, abs_tol: float, rel_tol: float) -> str:
    if not (math.isfinite(float(abs_diff)) and math.isfinite(float(rel_diff))):
        return "fail"
    return "pass" if abs_diff <= abs_tol and rel_diff <= rel_tol else "fail"


def eb_sorted_roots(beta_deg: float, mu: float, epsilon: float) -> np.ndarray:
    return find_first_n_roots(
        beta=float(np.deg2rad(float(beta_deg))),
        mu=float(mu),
        eps=float(epsilon),
        n_roots=N_MODES,
        Lmin=0.2,
        Lmax0=35.0,
        scan_step=0.005,
        grow_factor=1.35,
        max_tries=8,
    )


def max_finite(values: Sequence[float]) -> float:
    finite = [abs(float(value)) for value in values if math.isfinite(float(value))]
    return max(finite) if finite else float("nan")


def classify_decrease(values: Sequence[float]) -> str:
    finite_values = [float(value) for value in values if math.isfinite(float(value))]
    if len(finite_values) != len(values):
        return "fail_missing_value"
    monotone = all(
        finite_values[index + 1] <= finite_values[index] * (1.0 + 1.0e-9) + 1.0e-12
        for index in range(len(finite_values) - 1)
    )
    if monotone:
        return "pass_monotone_decrease"
    final_is_smaller = finite_values[-1] <= 0.1 * finite_values[0]
    final_is_lowest = finite_values[-1] <= min(finite_values[:-1]) * (1.0 + 0.05)
    if final_is_smaller and final_is_lowest:
        return "pass_general_decrease"
    return "fail_no_clear_decrease"


def compute_test_a(rows: list[dict[str, object]], warnings: list[str]) -> dict[str, object]:
    max_abs = 0.0
    max_rel = 0.0
    fail_count = 0
    worst: dict[str, object] | None = None
    for epsilon in TEST_A_EPSILONS:
        eb_reference = equal_timo.eb_roots(float(epsilon))
        timo_reference, reference_warnings = equal_timo.timo_roots(float(epsilon), eb_reference)
        timo_variable, variable_warnings = timo_sorted_roots(15.0, 0.0, float(epsilon), N_MODES)
        warnings.extend(f"test A epsilon={epsilon:g}: reference: {item}" for item in reference_warnings)
        warnings.extend(f"test A epsilon={epsilon:g}: variable: {item}" for item in variable_warnings)
        for mode_index in range(N_MODES):
            lambda_ref = float(timo_reference[mode_index])
            lambda_var = float(timo_variable[mode_index])
            abs_diff, rel_diff = finite_abs_rel(lambda_var, lambda_ref)
            status = status_from_tolerance(abs_diff, rel_diff, A_ABS_TOL, A_REL_TOL)
            if status != "pass":
                fail_count += 1
            if math.isfinite(abs_diff) and abs_diff >= max_abs:
                max_abs = abs_diff
                worst = {
                    "epsilon": float(epsilon),
                    "mode": mode_index + 1,
                    "abs_diff": abs_diff,
                    "rel_diff": rel_diff,
                }
            if math.isfinite(rel_diff):
                max_rel = max(max_rel, rel_diff)
            rows.append(
                {
                    "test_name": "A_mu0_equal_rods_timoshenko_consistency",
                    "beta": 15.0,
                    "mu": 0.0,
                    "epsilon": float(epsilon),
                    "mode": mode_index + 1,
                    "Lambda_EB": float(eb_reference[mode_index]),
                    "Lambda_Timo_variable": lambda_var,
                    "Lambda_Timo_reference": lambda_ref,
                    "Lambda_reference": lambda_ref,
                    "abs_diff": abs_diff,
                    "rel_diff": rel_diff,
                    "status": status,
                    "notes": f"abs_tol={A_ABS_TOL:g}; rel_tol={A_REL_TOL:g}",
                }
            )
    return {
        "status": "pass" if fail_count == 0 else "fail",
        "fail_count": fail_count,
        "max_abs": max_abs,
        "max_rel": max_rel,
        "worst": worst,
    }


def compute_test_b(rows: list[dict[str, object]], warnings: list[str]) -> dict[str, object]:
    summaries: list[dict[str, object]] = []
    fail_count = 0
    for beta in TEST_B_BETAS:
        for mu in TEST_B_MUS:
            case_rows: list[dict[str, object]] = []
            max_rel_by_epsilon: list[float] = []
            max_abs_by_epsilon: list[float] = []
            for epsilon in TEST_B_EPSILONS:
                eb_values = eb_sorted_roots(float(beta), float(mu), float(epsilon))
                timo_values, root_warnings = timo_sorted_roots(float(beta), float(mu), float(epsilon), N_MODES)
                warnings.extend(
                    f"test B beta={beta:g}, mu={mu:g}, epsilon={epsilon:g}: {item}" for item in root_warnings
                )
                rel_values: list[float] = []
                abs_values: list[float] = []
                for mode_index in range(N_MODES):
                    lambda_eb = float(eb_values[mode_index])
                    lambda_timo = float(timo_values[mode_index])
                    abs_diff, rel_diff = finite_abs_rel(lambda_timo, lambda_eb)
                    abs_values.append(abs_diff)
                    rel_values.append(rel_diff)
                    case_rows.append(
                        {
                            "test_name": "B_epsilon_to_zero_EB_limit",
                            "beta": float(beta),
                            "mu": float(mu),
                            "epsilon": float(epsilon),
                            "mode": mode_index + 1,
                            "Lambda_EB": lambda_eb,
                            "Lambda_Timo_variable": lambda_timo,
                            "Lambda_Timo_reference": "",
                            "Lambda_reference": lambda_eb,
                            "abs_diff": abs_diff,
                            "rel_diff": rel_diff,
                            "status": "pending",
                            "notes": "sorted roots; not descendant branches",
                        }
                    )
                max_rel_by_epsilon.append(max_finite(rel_values))
                max_abs_by_epsilon.append(max_finite(abs_values))
            trend = classify_decrease(max_rel_by_epsilon)
            case_status = "pass" if trend.startswith("pass_") else "fail"
            if case_status != "pass":
                fail_count += 1
            for row in case_rows:
                row["status"] = case_status
                row["notes"] = f"{row['notes']}; trend={trend}"
            rows.extend(case_rows)
            summaries.append(
                {
                    "beta": float(beta),
                    "mu": float(mu),
                    "max_rel_by_epsilon": tuple(max_rel_by_epsilon),
                    "max_abs_by_epsilon": tuple(max_abs_by_epsilon),
                    "trend": trend,
                    "status": case_status,
                }
            )
    return {
        "status": "pass" if fail_count == 0 else "fail",
        "fail_count": fail_count,
        "summaries": summaries,
    }


def compute_test_c(rows: list[dict[str, object]], warnings: list[str]) -> dict[str, object]:
    summaries: list[dict[str, object]] = []
    fail_count = 0
    for epsilon in TEST_C_EPSILONS:
        eb_baseline = eb_sorted_roots(0.0, 0.0, float(epsilon))
        timo_baseline, baseline_warnings = timo_sorted_roots(0.0, 0.0, float(epsilon), N_MODES)
        warnings.extend(f"test C epsilon={epsilon:g}: baseline: {item}" for item in baseline_warnings)
        eb_abs_diffs: list[float] = []
        eb_rel_diffs: list[float] = []
        timo_abs_diffs: list[float] = []
        timo_rel_diffs: list[float] = []
        for mu in TEST_C_MUS:
            eb_values = eb_sorted_roots(0.0, float(mu), float(epsilon))
            timo_values, root_warnings = timo_sorted_roots(0.0, float(mu), float(epsilon), N_MODES)
            warnings.extend(f"test C epsilon={epsilon:g}, mu={mu:g}: {item}" for item in root_warnings)
            for mode_index in range(N_MODES):
                lambda_eb = float(eb_values[mode_index])
                lambda_eb_ref = float(eb_baseline[mode_index])
                eb_abs, eb_rel = finite_abs_rel(lambda_eb, lambda_eb_ref)
                eb_status = status_from_tolerance(eb_abs, eb_rel, C_EB_ABS_TOL, C_EB_REL_TOL)
                if eb_status != "pass":
                    fail_count += 1
                eb_abs_diffs.append(eb_abs)
                eb_rel_diffs.append(eb_rel)
                rows.append(
                    {
                        "test_name": "C_beta0_straight_rod_EB_mu_invariance",
                        "beta": 0.0,
                        "mu": float(mu),
                        "epsilon": float(epsilon),
                        "mode": mode_index + 1,
                        "Lambda_EB": lambda_eb,
                        "Lambda_Timo_variable": "",
                        "Lambda_Timo_reference": "",
                        "Lambda_reference": lambda_eb_ref,
                        "abs_diff": eb_abs,
                        "rel_diff": eb_rel,
                        "status": eb_status,
                        "notes": f"reference=beta0_mu0_EB; abs_tol={C_EB_ABS_TOL:g}; rel_tol={C_EB_REL_TOL:g}",
                    }
                )

                lambda_timo = float(timo_values[mode_index])
                lambda_timo_ref = float(timo_baseline[mode_index])
                timo_abs, timo_rel = finite_abs_rel(lambda_timo, lambda_timo_ref)
                timo_status = status_from_tolerance(timo_abs, timo_rel, C_TIMO_ABS_TOL, C_TIMO_REL_TOL)
                if timo_status != "pass":
                    fail_count += 1
                timo_abs_diffs.append(timo_abs)
                timo_rel_diffs.append(timo_rel)
                rows.append(
                    {
                        "test_name": "C_beta0_straight_rod_Timoshenko_mu_invariance",
                        "beta": 0.0,
                        "mu": float(mu),
                        "epsilon": float(epsilon),
                        "mode": mode_index + 1,
                        "Lambda_EB": "",
                        "Lambda_Timo_variable": lambda_timo,
                        "Lambda_Timo_reference": lambda_timo_ref,
                        "Lambda_reference": lambda_timo_ref,
                        "abs_diff": timo_abs,
                        "rel_diff": timo_rel,
                        "status": timo_status,
                        "notes": f"reference=beta0_mu0_Timoshenko; abs_tol={C_TIMO_ABS_TOL:g}; rel_tol={C_TIMO_REL_TOL:g}",
                    }
                )
        eb_max_abs = max_finite(eb_abs_diffs)
        eb_max_rel = max_finite(eb_rel_diffs)
        timo_max_abs = max_finite(timo_abs_diffs)
        timo_max_rel = max_finite(timo_rel_diffs)
        summaries.append(
            {
                "epsilon": float(epsilon),
                "model": "EB",
                "max_abs": eb_max_abs,
                "max_rel": eb_max_rel,
                "status": "pass" if eb_max_abs <= C_EB_ABS_TOL and eb_max_rel <= C_EB_REL_TOL else "fail",
            }
        )
        summaries.append(
            {
                "epsilon": float(epsilon),
                "model": "Timoshenko",
                "max_abs": timo_max_abs,
                "max_rel": timo_max_rel,
                "status": "pass" if timo_max_abs <= C_TIMO_ABS_TOL and timo_max_rel <= C_TIMO_REL_TOL else "fail",
            }
        )
    return {
        "status": "pass" if fail_count == 0 else "fail",
        "fail_count": fail_count,
        "summaries": summaries,
    }


def write_csv(rows: Sequence[dict[str, object]]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def format_float(value: object, digits: int = 6) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(number):
        return "nan"
    return f"{number:.{digits}e}"


def test_b_table(summaries: Sequence[dict[str, object]]) -> list[str]:
    lines = [
        "| beta deg | mu | max rel @ eps=0.02 | max rel @ eps=0.01 | max rel @ eps=0.005 | max rel @ eps=0.0025 | trend | status |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for summary in summaries:
        rels = list(summary["max_rel_by_epsilon"])
        lines.append(
            "| "
            f"{float(summary['beta']):g} | "
            f"{float(summary['mu']):g} | "
            f"{format_float(rels[0])} | "
            f"{format_float(rels[1])} | "
            f"{format_float(rels[2])} | "
            f"{format_float(rels[3])} | "
            f"{summary['trend']} | "
            f"{summary['status']} |"
        )
    return lines


def test_c_table(summaries: Sequence[dict[str, object]]) -> list[str]:
    lines = [
        "| epsilon | model | max abs drift across mu | max relative drift across mu | status |",
        "| ---: | --- | ---: | ---: | --- |",
    ]
    for summary in summaries:
        lines.append(
            "| "
            f"{float(summary['epsilon']):g} | "
            f"{summary['model']} | "
            f"{format_float(summary['max_abs'])} | "
            f"{format_float(summary['max_rel'])} | "
            f"{summary['status']} |"
        )
    return lines


def warning_lines(warnings: Sequence[str]) -> list[str]:
    unique = sorted(set(warnings))
    if not unique:
        return ["No root-search or basis warnings were recorded."]
    lines = ["Warnings:"]
    for warning in unique[:40]:
        lines.append(f"- {warning}")
    if len(unique) > 40:
        lines.append(f"- ... {len(unique) - 40} more warnings")
    return lines


def write_report(
    test_a: dict[str, object],
    test_b: dict[str, object],
    test_c: dict[str, object],
    warnings: Sequence[str],
) -> None:
    overall_pass = all(result["status"] == "pass" for result in (test_a, test_b, test_c))
    validation_line = (
        "The variable-length Timoshenko implementation is diagnostic-verified for eta=0 sorted-root limit checks."
        if overall_pass
        else "The variable-length Timoshenko implementation remains not validated by this audit."
    )
    worst_a = test_a.get("worst") or {}
    lines: list[str] = [
        "# Variable-Length Timoshenko Limits Audit",
        "",
        "## Purpose",
        "",
        "Diagnostic-only verification of the variable-length Timoshenko determinant used by the corrected "
        "`Lambda(mu)` EB/Timoshenko/FEM diagnostic. This audit checks limiting cases only and does not "
        "change article files, article figures, the old determinant, `formulas.py`, old solvers, baseline "
        "results, or the existing FEM physical model.",
        "",
        "The checks use sorted roots. They are determinant-limit sanity checks, not descendant branch tracking.",
        "",
        "## Implemented Formula Audit",
        "",
        "```text",
        *FORMULA_AUDIT_LINES,
        "```",
        "",
        "## Test A: mu=0 Consistency",
        "",
        "Comparison: old equal-rods Timoshenko diagnostic at `beta=15 deg`, `mu=0`, `eta=0` versus the "
        "variable-length Timoshenko determinant evaluated at `mu=0`.",
        "",
        f"- status: {test_a['status']}",
        f"- max abs diff: {format_float(test_a['max_abs'])}",
        f"- max relative diff: {format_float(test_a['max_rel'])}",
        f"- failures: {test_a['fail_count']}",
    ]
    if worst_a:
        lines.append(
            f"- worst row: epsilon={float(worst_a['epsilon']):g}, mode={int(worst_a['mode'])}, "
            f"abs={format_float(worst_a['abs_diff'])}, rel={format_float(worst_a['rel_diff'])}"
        )
    lines.extend(
        [
            "",
            "## Test B: epsilon -> 0 EB Limit",
            "",
            "Comparison: existing EB determinant roots from `find_first_n_roots` versus variable-length "
            "Timoshenko sorted roots for the first six roots.",
            "",
            *test_b_table(test_b["summaries"]),
            "",
            f"- status: {test_b['status']}",
            f"- failed beta/mu trend cases: {test_b['fail_count']}",
            "",
            "Root-tracking caveat: these rows compare sorted roots at each parameter point. Near crossings or "
            "order exchanges would require shape-based descendant tracking before any branch-identity claim.",
            "",
            "## Test C: beta=0 Straight-Rod mu-Invariance",
            "",
            "For `beta=0`, `eta=0`, and equal material/thickness, the two collinear rods should be equivalent "
            "to one straight rod of total length 2, so sorted spectra should be independent of `mu`.",
            "",
            *test_c_table(test_c["summaries"]),
            "",
            f"- status: {test_c['status']}",
            f"- failed drift rows: {test_c['fail_count']}",
            "",
            "## Overall Status",
            "",
            validation_line,
            "",
            "The diagnostic-verified label is conditional on all three checks passing:",
            "",
            "- Test A `mu=0` consistency passes.",
            "- Test B `epsilon -> 0` EB limit passes.",
            "- Test C `beta=0` straight-rod `mu`-invariance passes.",
            "",
            "## What Remains Unverified",
            "",
            "- `eta != 0` tau scaling is not verified here.",
            "- Descendant branch identity is not verified here; this audit deliberately uses sorted roots.",
            "- FEM validation and point-joint solid constraints are not exercised here.",
            "- High-frequency behavior near or above the Timoshenko cut-off is not promoted by this audit.",
            "- This is numerical diagnostic verification of limiting cases, not a symbolic proof.",
            "",
            "## Outputs",
            "",
            f"- CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
            f"- report: `{OUTPUT_REPORT.relative_to(REPO_ROOT)}`",
            "",
            "## Warnings",
            "",
            *warning_lines(warnings),
            "",
        ]
    )
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    rows: list[dict[str, object]] = []
    warnings: list[str] = []
    test_a = compute_test_a(rows, warnings)
    test_b = compute_test_b(rows, warnings)
    test_c = compute_test_c(rows, warnings)
    write_csv(rows)
    write_report(test_a, test_b, test_c, warnings)

    overall_pass = all(result["status"] == "pass" for result in (test_a, test_b, test_c))
    print("Variable-length Timoshenko limits audit")
    print("Formula audit:")
    for line in FORMULA_AUDIT_LINES:
        print(f"  {line}")
    print(f"Test A status: {test_a['status']}; max abs={test_a['max_abs']:.6e}; max rel={test_a['max_rel']:.6e}")
    print(f"Test B status: {test_b['status']}; failed trend cases={test_b['fail_count']}")
    print(f"Test C status: {test_c['status']}; failed drift rows={test_c['fail_count']}")
    print(f"Overall status: {'pass' if overall_pass else 'fail'}")
    print(f"Wrote {OUTPUT_CSV.relative_to(REPO_ROOT)}")
    print(f"Wrote {OUTPUT_REPORT.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
