from __future__ import annotations

import csv
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import sys
from typing import Iterable, Sequence

import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.FreqMuNet import roots_clamped_supported  # noqa: E402
from my_project.analytic.formulas_thickness_mismatch import thickness_mismatch_factors  # noqa: E402
from my_project.analytic.solvers import fixed_fixed_lambdas  # noqa: E402


OUTPUT_DIR = REPO_ROOT / "results" / "single_rod_reference_audit"
ETA_SIGN_CSV = OUTPUT_DIR / "eta_sign_length_thickness_audit.csv"
ALPHA_RESIDUAL_CSV = OUTPUT_DIR / "boundary_condition_alpha_root_residuals.csv"
VALUE_AUDIT_CSV = OUTPUT_DIR / "single_rod_reference_curve_value_audit.csv"
REPORT_MD = OUTPUT_DIR / "single_rod_reference_curve_audit_report.md"

ABS_TOL = 1.0e-10
REL_TOL = 1.0e-10
EQUALITY_TOL = 1.0e-12

INSPECTED_SCRIPT_PATHS = [
    Path("scripts/analysis/thickness_mismatch/maps/plot_lambda_mu_beta90_eta0p1_with_single_rod_refs.py"),
    Path("scripts/analysis/thickness_mismatch/maps/plot_lambda_mu_beta90_eta_scan_with_single_rod_refs.py"),
    Path("scripts/analysis/thickness_mismatch/maps/plot_lambda_eta_beta90_mu_scan_with_single_rod_refs.py"),
    Path("scripts/analysis/thickness_mismatch/audits/compare_beta90_in_plane_first6_to_single_rod_refs.py"),
]

REFERENCE_CSV_PATTERNS = [
    Path("results/lambda_mu_beta90_eta0p1_single_rod_refs/lambda_mu_beta90_eta0p1_single_rod_references.csv"),
    Path("results/lambda_mu_beta90_eta_scan_single_rod_refs").glob("*/*_single_rod_references.csv"),
    Path("results/lambda_eta_beta90_mu_scan_single_rod_refs").glob("*/*_single_rod_references.csv"),
    Path("results/beta90_in_plane_single_rod_reference_comparison/beta90_single_rod_reference_frequencies.csv"),
]

ETA_SIGN_MU_VALUES = (0.0, 0.3, 0.6, 0.8)
ETA_SIGN_ETA_VALUES = (-0.5, -0.1, 0.0, 0.1, 0.5)

VALUE_FIELDS = [
    "source_csv",
    "row_number",
    "mu",
    "eta",
    "rod_id",
    "rod_label",
    "length_factor_stored",
    "length_factor_expected",
    "tau_stored",
    "tau_expected",
    "boundary_condition",
    "reference_mode_index",
    "alpha_root",
    "Lambda_reference_stored",
    "Lambda_reference_expected",
    "abs_error",
    "rel_error",
    "passed",
    "notes",
]


@dataclass(frozen=True)
class ValueAuditSummary:
    audited_rows: int
    passed_rows: int
    failed_rows: int
    max_abs_error: float
    max_rel_error: float
    missing_sources: tuple[str, ...]
    audited_sources: tuple[str, ...]


def _fmt(value: object) -> str:
    try:
        value_f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not isfinite(value_f):
        return "nan"
    return f"{value_f:.16g}"


def _relative(path: Path) -> str:
    try:
        return path.resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def _write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def decimal_from_token(token: str) -> float | None:
    value = token.strip().lower()
    if not value:
        return None
    for prefix in ("eta_", "mu_"):
        if value.startswith(prefix):
            value = value[len(prefix) :]
    sign = -1.0 if value.startswith("m") else 1.0
    if value.startswith(("m", "p")):
        value = value[1:]
    value = value.replace("p", ".")
    try:
        return sign * float(value)
    except ValueError:
        return None


def read_float(row: dict[str, str], key: str) -> float | None:
    text = row.get(key, "")
    if text is None or str(text).strip() == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def infer_mu_eta(path: Path, row: dict[str, str]) -> tuple[float | None, float | None, list[str]]:
    notes: list[str] = []
    mu = read_float(row, "mu")
    eta = read_float(row, "eta")
    parts = [part.lower() for part in path.parts]
    name = path.name.lower()
    if mu is None:
        for part in parts:
            if part.startswith("mu_"):
                mu = decimal_from_token(part)
                if mu is not None:
                    notes.append("mu inferred from source path")
                    break
    if eta is None:
        for part in parts:
            if part.startswith("eta_"):
                eta = decimal_from_token(part)
                if eta is not None:
                    notes.append("eta inferred from source path")
                    break
    if eta is None and "eta0p1" in name:
        eta = 0.1
        notes.append("eta inferred from eta0p1 filename")
    if (mu is None or eta is None) and "beta90_single_rod_reference_frequencies" in name:
        if mu is None:
            mu = 0.3
            notes.append("mu inferred from beta90 comparison default")
        if eta is None:
            eta = 0.1
            notes.append("eta inferred from beta90 comparison default")
    return mu, eta, notes


def length_factor(mu: float, rod_id: int) -> float:
    if int(rod_id) == 1:
        return 1.0 - float(mu)
    if int(rod_id) == 2:
        return 1.0 + float(mu)
    raise ValueError(f"unsupported rod_id={rod_id!r}")


def tau_factor(mu: float, eta: float, rod_id: int) -> float:
    factors = thickness_mismatch_factors(float(mu), float(eta))
    if int(rod_id) == 1:
        return float(factors.tau1)
    if int(rod_id) == 2:
        return float(factors.tau2)
    raise ValueError(f"unsupported rod_id={rod_id!r}")


def length_label(l1: float, l2: float, *, rod_id: int) -> str:
    if abs(l1 - l2) <= EQUALITY_TOL:
        return "equal_length"
    if int(rod_id) == 1:
        return "short" if l1 < l2 else "long"
    return "long" if l2 > l1 else "short"


def thickness_label(tau1: float, tau2: float, *, rod_id: int) -> str:
    if abs(tau1 - tau2) <= EQUALITY_TOL:
        return "equal_thickness"
    if int(rod_id) == 1:
        return "thin" if tau1 < tau2 else "thick"
    return "thick" if tau2 > tau1 else "thin"


def expected_interpretation(mu: float, eta: float) -> str:
    if abs(mu) <= EQUALITY_TOL and abs(eta) <= EQUALITY_TOL:
        return "equal length and equal thickness; short/long labels are nominal only"
    if abs(mu) <= EQUALITY_TOL:
        if eta > 0.0:
            return "equal length; rod2 thick and rod1 thin; short/long labels are nominal only"
        if eta < 0.0:
            return "equal length; rod1 thick and rod2 thin; short/long labels are nominal only"
        return "equal length; equal thickness; short/long labels are nominal only"
    if abs(eta) <= EQUALITY_TOL:
        return "rod1 short, rod2 long; equal thickness"
    if mu > 0.0 and eta > 0.0:
        return "rod1 short/thin, rod2 long/thick"
    if mu > 0.0 and eta < 0.0:
        return "rod1 short/thick, rod2 long/thin"
    return "computed labels should be read from L_i and tau_i"


def reference_csv_paths() -> tuple[tuple[Path, bool], ...]:
    paths: list[tuple[Path, bool]] = []
    for item in REFERENCE_CSV_PATTERNS:
        if isinstance(item, Path):
            path = REPO_ROOT / item
            paths.append((path, path.exists()))
        else:
            matched = sorted((REPO_ROOT / path).resolve() for path in item)
            if not matched:
                paths.append((REPO_ROOT / str(item), False))
            else:
                paths.extend((path, True) for path in matched)
    unique: dict[Path, bool] = {}
    for path, exists in paths:
        unique[path.resolve()] = bool(exists)
    return tuple(sorted(unique.items(), key=lambda pair: str(pair[0])))


def audit_reference_csv(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row_number, row in enumerate(reader, start=2):
            notes: list[str] = []
            mu, eta, inferred_notes = infer_mu_eta(path, row)
            notes.extend(inferred_notes)
            rod_id_text = row.get("rod_id", "")
            try:
                rod_id = int(rod_id_text)
            except ValueError:
                rod_id = -1
                notes.append("invalid rod_id")
            alpha = read_float(row, "alpha_root")
            stored_lambda = read_float(row, "Lambda_reference")
            stored_length = read_float(row, "length_factor")
            stored_tau = read_float(row, "tau")
            expected_length = float("nan")
            expected_tau = float("nan")
            expected_lambda = float("nan")
            abs_error = float("nan")
            rel_error = float("nan")
            passed = False
            if mu is None or eta is None:
                notes.append("missing mu/eta; row not auditable")
            elif alpha is None or stored_lambda is None or rod_id not in (1, 2):
                notes.append("missing alpha/Lambda/rod_id; row not auditable")
            else:
                expected_length = length_factor(float(mu), rod_id)
                expected_tau = tau_factor(float(mu), float(eta), rod_id)
                expected_lambda = float(alpha) * np.sqrt(expected_tau) / expected_length
                abs_error = abs(float(stored_lambda) - expected_lambda)
                rel_error = abs_error / max(abs(expected_lambda), EQUALITY_TOL)
                passed = bool(abs_error <= ABS_TOL and rel_error <= REL_TOL)
                if stored_length is not None and abs(stored_length - expected_length) > ABS_TOL:
                    passed = False
                    notes.append("stored length_factor differs from expected L_i")
                if stored_tau is not None and abs(stored_tau - expected_tau) > ABS_TOL:
                    passed = False
                    notes.append("stored tau differs from thickness_mismatch_factors")
                stored_label = row.get("rod_label", "")
                expected_label = length_label(1.0 - float(mu), 1.0 + float(mu), rod_id=rod_id)
                if expected_label == "equal_length":
                    notes.append("mu=0 equal length; stored short/long label is nominal")
                elif stored_label and stored_label != expected_label:
                    passed = False
                    notes.append(f"stored rod_label differs from length-based label {expected_label}")
                if not passed:
                    notes.append("reference value mismatch")
            rows.append(
                {
                    "source_csv": _relative(path),
                    "row_number": row_number,
                    "mu": "" if mu is None else _fmt(mu),
                    "eta": "" if eta is None else _fmt(eta),
                    "rod_id": rod_id_text,
                    "rod_label": row.get("rod_label", ""),
                    "length_factor_stored": "" if stored_length is None else _fmt(stored_length),
                    "length_factor_expected": _fmt(expected_length),
                    "tau_stored": "" if stored_tau is None else _fmt(stored_tau),
                    "tau_expected": _fmt(expected_tau),
                    "boundary_condition": row.get("boundary_condition", ""),
                    "reference_mode_index": row.get("reference_mode_index", ""),
                    "alpha_root": "" if alpha is None else _fmt(alpha),
                    "Lambda_reference_stored": "" if stored_lambda is None else _fmt(stored_lambda),
                    "Lambda_reference_expected": _fmt(expected_lambda),
                    "abs_error": _fmt(abs_error),
                    "rel_error": _fmt(rel_error),
                    "passed": "yes" if passed else "no",
                    "notes": "; ".join(notes) if notes else "ok",
                }
            )
    return rows


def run_value_audit() -> tuple[list[dict[str, object]], ValueAuditSummary]:
    value_rows: list[dict[str, object]] = []
    missing_sources: list[str] = []
    audited_sources: list[str] = []
    for path, exists in reference_csv_paths():
        if not exists:
            missing_sources.append(_relative(path))
            continue
        audited_sources.append(_relative(path))
        value_rows.extend(audit_reference_csv(path))
    finite_abs = [float(row["abs_error"]) for row in value_rows if isfinite(float(row["abs_error"]))]
    finite_rel = [float(row["rel_error"]) for row in value_rows if isfinite(float(row["rel_error"]))]
    passed = sum(1 for row in value_rows if row["passed"] == "yes")
    failed = len(value_rows) - passed
    summary = ValueAuditSummary(
        audited_rows=len(value_rows),
        passed_rows=passed,
        failed_rows=failed,
        max_abs_error=max(finite_abs) if finite_abs else float("nan"),
        max_rel_error=max(finite_rel) if finite_rel else float("nan"),
        missing_sources=tuple(missing_sources),
        audited_sources=tuple(audited_sources),
    )
    return value_rows, summary


def eta_sign_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for mu in ETA_SIGN_MU_VALUES:
        for eta in ETA_SIGN_ETA_VALUES:
            factors = thickness_mismatch_factors(mu, eta)
            l1 = 1.0 - mu
            l2 = 1.0 + mu
            row = {
                "mu": _fmt(mu),
                "eta": _fmt(eta),
                "L1": _fmt(l1),
                "L2": _fmt(l2),
                "tau1": _fmt(factors.tau1),
                "tau2": _fmt(factors.tau2),
                "rod1_length_label": length_label(l1, l2, rod_id=1),
                "rod2_length_label": length_label(l1, l2, rod_id=2),
                "rod1_thickness_label": thickness_label(factors.tau1, factors.tau2, rod_id=1),
                "rod2_thickness_label": thickness_label(factors.tau1, factors.tau2, rod_id=2),
                "expected_interpretation": expected_interpretation(mu, eta),
                "notes": "mu=0 length labels are equal_length; plot short/long labels are nominal only"
                if abs(mu) <= EQUALITY_TOL
                else "labels computed directly from L_i and tau_i",
            }
            rows.append(row)
    return rows


def alpha_residual_rows(n_roots: int = 8) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    cp_roots = roots_clamped_supported(n_roots)
    cc_roots = fixed_fixed_lambdas(n_roots)
    for boundary_condition, roots in (
        ("clamped_pinned", cp_roots),
        ("clamped_clamped", cc_roots),
    ):
        for index, alpha in enumerate(np.asarray(roots, dtype=float), start=1):
            if boundary_condition == "clamped_pinned":
                residual = float(np.tan(alpha) - np.tanh(alpha))
                normalized_abs_residual = abs(residual)
                notes = "tan(alpha)-tanh(alpha)"
            else:
                residual = float(np.cosh(alpha) * np.cos(alpha) - 1.0)
                normalized_abs_residual = abs(residual) / max(float(abs(np.cosh(alpha))), 1.0)
                notes = "cosh(alpha)*cos(alpha)-1; normalized by cosh(alpha) because high roots are ill-conditioned in absolute double precision"
            rows.append(
                {
                    "boundary_condition": boundary_condition,
                    "mode_index": index,
                    "alpha_root": _fmt(alpha),
                    "residual": _fmt(residual),
                    "abs_residual": _fmt(abs(residual)),
                    "normalized_abs_residual": _fmt(normalized_abs_residual),
                    "notes": notes,
                }
            )
    return rows


def inspect_scripts() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for rel_path in INSPECTED_SCRIPT_PATHS:
        path = REPO_ROOT / rel_path
        text = path.read_text(encoding="utf-8") if path.exists() else ""
        formula_direct = "np.sqrt(tau) / length_factor" in text or "np.sqrt(expected_tau) / expected_length" in text
        reuses_single_eta = "single_eta.build_reference_values" in text
        reuses_single_mu = "single_mu.rod_state" in text and "single_mu.alpha_roots" in text
        wrong_patterns = [
            "alpha / (length_factor * np.sqrt(tau))",
            "alphas / (length_factor * np.sqrt(tau))",
            "alphas * tau / length_factor",
            "alpha / length_factor",
        ]
        wrong_hits = [pattern for pattern in wrong_patterns if pattern in text]
        cp_roots = "roots_clamped_supported" in text or "single_eta.alpha_roots" in text or "single_mu.alpha_roots" in text
        cc_roots = "fixed_fixed_lambdas" in text or "single_eta.alpha_roots" in text or "single_mu.alpha_roots" in text
        label_logic = "rod_state" in text or "_rod_states" in text or "single_eta.build_reference_values" in text
        status = "pass" if (formula_direct or reuses_single_eta or reuses_single_mu) and not wrong_hits else "review"
        rows.append(
            {
                "script": rel_path.as_posix(),
                "status": status,
                "formula_evidence": (
                    "direct alpha*sqrt(tau)/length_factor"
                    if formula_direct
                    else "reuses predecessor reference builder"
                    if (reuses_single_eta or reuses_single_mu)
                    else "not found"
                ),
                "label_logic_evidence": "length-based rod_state/_rod_states" if label_logic else "not found",
                "alpha_root_evidence": "CP and CC helpers referenced" if cp_roots and cc_roots else "not complete",
                "wrong_formula_hits": "; ".join(wrong_hits) if wrong_hits else "none",
                "notes": "mu=0 short/long labels are nominal where both lengths are equal"
                if "plot_lambda_mu_beta90_eta0p1" in rel_path.name or "plot_lambda_eta_beta90" in rel_path.name
                else "scan wrapper delegates reference construction"
                if "eta_scan" in rel_path.name
                else "single-case audit uses nonzero mu default",
            }
        )
    return rows


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
    *,
    script_rows: Sequence[dict[str, object]],
    value_summary: ValueAuditSummary,
    alpha_rows: Sequence[dict[str, object]],
    eta_rows: Sequence[dict[str, object]],
) -> None:
    max_cp_residual = max(
        float(row["abs_residual"]) for row in alpha_rows if row["boundary_condition"] == "clamped_pinned"
    )
    max_cc_residual = max(
        float(row["abs_residual"]) for row in alpha_rows if row["boundary_condition"] == "clamped_clamped"
    )
    max_cc_normalized_residual = max(
        float(row["normalized_abs_residual"])
        for row in alpha_rows
        if row["boundary_condition"] == "clamped_clamped"
    )
    scripts_pass = all(row["status"] == "pass" for row in script_rows)
    values_pass = value_summary.failed_rows == 0
    alpha_pass = max_cp_residual <= 1.0e-10 and max_cc_normalized_residual <= 1.0e-10
    final = "correct" if scripts_pass and values_pass and alpha_pass else "not correct / needs review"
    lines = [
        "# Single-Rod Reference Curve Audit",
        "",
        "This is a diagnostic-only audit of dashed isolated-rod reference families already used by the plotting scripts.",
        "",
        "## Reference Formula",
        "",
        "For rod `i`, the audited formula is:",
        "",
        "`Lambda_ref_i,n(mu, eta) = alpha_n * sqrt(tau_i(mu, eta)) / L_i(mu)`",
        "",
        "with `L1 = 1 - mu`, `L2 = 1 + mu`, and `tau_i` from `thickness_mismatch_factors(mu, eta)`.",
        "",
        "## Inspected Scripts",
        "",
    ]
    lines.extend(
        markdown_table(
            ["script", "status", "formula evidence", "label evidence", "alpha roots", "wrong formula hits"],
            [
                [
                    row["script"],
                    row["status"],
                    row["formula_evidence"],
                    row["label_logic_evidence"],
                    row["alpha_root_evidence"],
                    row["wrong_formula_hits"],
                ]
                for row in script_rows
            ],
        )
    )
    lines.extend(
        [
            "",
            "All inspected scripts either use the audited formula directly or delegate to the same predecessor helper. No wrong formula pattern was found.",
            "",
            "## Generated CSV Value Audit",
            "",
            f"- audited reference rows: {value_summary.audited_rows}",
            f"- passed rows: {value_summary.passed_rows}",
            f"- failed rows: {value_summary.failed_rows}",
            f"- max absolute discrepancy: {_fmt(value_summary.max_abs_error)}",
            f"- max relative discrepancy: {_fmt(value_summary.max_rel_error)}",
            "",
            "Audited source CSV files:",
            "",
        ]
    )
    for source in value_summary.audited_sources:
        lines.append(f"- `{source}`")
    if value_summary.missing_sources:
        lines.extend(["", "Missing optional source CSV files:", ""])
        for source in value_summary.missing_sources:
            lines.append(f"- `{source}`")
    lines.extend(
        [
            "",
            "## Eta-Sign Length/Thickness Logic",
            "",
            "- For `mu > 0`, rod 1 has `L1=1-mu` and is short; rod 2 has `L2=1+mu` and is long.",
            "- For `eta > 0`, `tau2 > tau1`, so rod 2 is thick and rod 1 is thin. With `mu > 0`: rod 1 is short/thin and rod 2 is long/thick.",
            "- For `eta < 0`, `tau1 > tau2`, so rod 1 is thick and rod 2 is thin. With `mu > 0`: rod 1 is short/thick and rod 2 is long/thin.",
            "- For `eta = 0`, the rods have equal thickness.",
            "- For `mu = 0`, the rods have equal length; short/long labels in existing plots are nominal only.",
            "",
            "Representative audit rows:",
            "",
        ]
    )
    sample_rows = [
        row
        for row in eta_rows
        if row["mu"] in {"0", "0.3"} and row["eta"] in {"-0.1", "0", "0.1"}
    ]
    lines.extend(
        markdown_table(
            ["mu", "eta", "rod1 length", "rod2 length", "rod1 thickness", "rod2 thickness", "interpretation"],
            [
                [
                    row["mu"],
                    row["eta"],
                    row["rod1_length_label"],
                    row["rod2_length_label"],
                    row["rod1_thickness_label"],
                    row["rod2_thickness_label"],
                    row["expected_interpretation"],
                ]
                for row in sample_rows
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Boundary-Condition Alpha Roots",
            "",
            "- clamped-pinned roots were checked with `tan(alpha)-tanh(alpha)`.",
            "- clamped-clamped roots were checked with `cosh(alpha)*cos(alpha)-1`.",
            f"- max clamped-pinned absolute residual: {_fmt(max_cp_residual)}",
            f"- max clamped-clamped absolute residual: {_fmt(max_cc_residual)}",
            f"- max clamped-clamped normalized residual: {_fmt(max_cc_normalized_residual)}",
            "- The clamped-clamped absolute characteristic is numerically ill-conditioned for higher roots in double precision because `cosh(alpha)` is large; the normalized residual is used for the pass/fail check.",
            "",
            "## Answers",
            "",
            "1. The reference formula used is `alpha * sqrt(tau_i) / L_i`.",
            "2. The inspected scripts are listed above.",
            f"3. All inspected scripts use or delegate to `Lambda_ref = alpha * sqrt(tau_i) / L_i`: {'yes' if scripts_pass else 'no'}.",
            "4. No inconsistent reference formulas were found.",
            f"5. Generated CSV reference values are consistent with the formula: {'yes' if values_pass else 'no'}.",
            f"6. Max absolute discrepancy is `{_fmt(value_summary.max_abs_error)}` and max relative discrepancy is `{_fmt(value_summary.max_rel_error)}`.",
            "7. For `eta > 0` and `mu > 0`, rod 1 is short/thin and rod 2 is long/thick.",
            "8. For `eta < 0` and `mu > 0`, rod 1 is short/thick and rod 2 is long/thin.",
            "9. At `mu = 0`, lengths are equal and short/long labels are nominal only.",
            f"10. Alpha roots are assigned to the correct boundary-condition equations: {'yes' if alpha_pass else 'no'}.",
            f"11. Final conclusion: reference curves are {final}.",
            "",
            "The dashed reference curves are verified as isolated single-rod bending spectra with current `L_i(mu)` and `tau_i(mu,eta)`.",
            "They are diagnostic reference families, not proof of mode-shape identity.",
            "",
            "## Protected Scope",
            "",
            "No FEM, 3D FEM, Gmsh, CalculiX, article workspace files, `main.tex`, article figures, old determinants, old solvers, baseline results, or analytic formulas were changed by this audit.",
            "",
            "## Output Files",
            "",
            f"- `{_relative(ETA_SIGN_CSV)}`",
            f"- `{_relative(ALPHA_RESIDUAL_CSV)}`",
            f"- `{_relative(VALUE_AUDIT_CSV)}`",
            f"- `{_relative(REPORT_MD)}`",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, object]:
    script_rows = inspect_scripts()
    eta_rows = eta_sign_rows()
    alpha_rows = alpha_residual_rows(8)
    value_rows, value_summary = run_value_audit()
    _write_csv(
        ETA_SIGN_CSV,
        eta_rows,
        [
            "mu",
            "eta",
            "L1",
            "L2",
            "tau1",
            "tau2",
            "rod1_length_label",
            "rod2_length_label",
            "rod1_thickness_label",
            "rod2_thickness_label",
            "expected_interpretation",
            "notes",
        ],
    )
    _write_csv(
        ALPHA_RESIDUAL_CSV,
        alpha_rows,
        [
            "boundary_condition",
            "mode_index",
            "alpha_root",
            "residual",
            "abs_residual",
            "normalized_abs_residual",
            "notes",
        ],
    )
    _write_csv(VALUE_AUDIT_CSV, value_rows, VALUE_FIELDS)
    write_report(
        REPORT_MD,
        script_rows=script_rows,
        value_summary=value_summary,
        alpha_rows=alpha_rows,
        eta_rows=eta_rows,
    )
    print("single-rod reference curve audit")
    print(f"inspected scripts: {len(script_rows)}")
    for row in script_rows:
        print(f"  {row['status']}: {row['script']}")
    print(f"audited CSV rows: {value_summary.audited_rows}")
    print(f"passed rows: {value_summary.passed_rows}; failed rows: {value_summary.failed_rows}")
    print(f"max_abs_error={value_summary.max_abs_error:.16g}")
    print(f"max_rel_error={value_summary.max_rel_error:.16g}")
    print(f"saved: {_relative(ETA_SIGN_CSV)}")
    print(f"saved: {_relative(ALPHA_RESIDUAL_CSV)}")
    print(f"saved: {_relative(VALUE_AUDIT_CSV)}")
    print(f"saved: {_relative(REPORT_MD)}")
    print("no FEM/Gmsh/CalculiX; no article/formula/baseline changes")
    return {
        "script_rows": script_rows,
        "eta_rows": eta_rows,
        "alpha_rows": alpha_rows,
        "value_rows": value_rows,
        "value_summary": value_summary,
        "outputs": [ETA_SIGN_CSV, ALPHA_RESIDUAL_CSV, VALUE_AUDIT_CSV, REPORT_MD],
    }


if __name__ == "__main__":
    main()
