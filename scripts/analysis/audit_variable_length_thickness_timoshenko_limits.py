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

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    find_first_n_roots_eta,
    thickness_mismatch_factors,
)
from scripts.lib.variable_length_timoshenko import (  # noqa: E402
    endpoint_columns,
    find_roots_by_sign_scan,
    lambda_cutoff,
    normalized_det,
    omega_cutoff,
    project_omega,
    section_from_epsilon,
    section_from_epsilon_tau,
    segment_lengths,
    tau_factors,
    timo_basis,
    timo_coupling_matrix,
    timo_sorted_roots,
)


N_MODES = 6
OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_CSV = OUTPUT_DIR / "variable_length_thickness_timoshenko_limits_audit.csv"
OUTPUT_REPORT = OUTPUT_DIR / "variable_length_thickness_timoshenko_limits_audit.md"

E_BETAS = [0.0, 15.0, 45.0]
E_MUS = [0.0, 0.3, 0.6]
E_EPSILONS = [0.0025, 0.01, 0.025]

F_ETAS = [-1.0e-3, -1.0e-4, 1.0e-4, 1.0e-3]
F_BETAS = [15.0]
F_MUS = [0.0, 0.3, 0.6]
F_EPSILONS = [0.01, 0.025]

G_BETAS = [15.0, 45.0]
G_EPSILONS = [0.01, 0.025]
G_PAIRS = [(0.3, 0.5), (0.6, 0.5), (0.3, -0.5)]

H_BETAS = [15.0, 45.0]
H_MUS = [0.0, 0.3, 0.6]
H_ETAS = [-0.5, 0.5]
H_EPSILONS = [0.02, 0.01, 0.005, 0.0025]

I_EPSILONS = [0.01, 0.025]
I_ETAS = [-0.5, 0.5]
I_MUS = [-0.6, -0.3, 0.0, 0.3, 0.6]

E_ABS_TOL = 1.0e-8
E_REL_TOL = 1.0e-9
G_ABS_TOL = 2.0e-6
G_REL_TOL = 2.0e-7
I_SWAP_ABS_TOL = 2.0e-6
I_SWAP_REL_TOL = 2.0e-7
I_ETA0_ABS_TOL = 1.0e-6
I_ETA0_REL_TOL = 1.0e-7

CSV_FIELDNAMES = [
    "test_name",
    "beta",
    "mu",
    "eta",
    "mu_reference",
    "eta_reference",
    "epsilon",
    "mode",
    "Lambda_Timo",
    "Lambda_reference",
    "Lambda_EB",
    "abs_diff",
    "rel_diff",
    "status",
    "notes",
]

FORMULA_AUDIT_LINES = [
    "eta = (r2 - r1)/(r1 + r2).",
    "denom = sqrt(1 + 2*mu*eta + eta^2).",
    "tau1 = (1 - eta)/denom, tau2 = (1 + eta)/denom.",
    "Mass conservation: (1 - mu)*tau1^2 + (1 + mu)*tau2^2 = 2.",
    "section_from_epsilon_tau(epsilon, tau) uses r_i = 2*epsilon*tau_i with L_segment = 1.",
    "Circular A_i = A0*tau_i^2 and I_i = I0*tau_i^4.",
    "K_i = kappa*G*A_i, B_i = E*I_i, m_i = rho*A_i, j_i = rho*I_i.",
    "Each rod has its own Timoshenko characteristic equation using B_i, K_i, m_i, and j_i.",
    "The frequency normalization stays Omega = epsilon*Lambda^2, based on the reference radius r0.",
    "In the EB limit, the bending wave number tends to Lambda/sqrt(tau_i).",
    "The axial basis uses theta = epsilon*Lambda^2; tau enters axial rows through N_i = E*A_i*u_i'.",
    "Moment rows use M_i = E*I_i*psi_i'.",
    "Shear rows use Q_i = kappa*G*A_i*(w_i' - psi_i).",
    "Cut-off per arm is Omega_c,i = sqrt(K_i/j_i), Lambda_c,i = sqrt(Omega_c,i/epsilon).",
    "The current audit uses sorted roots only and makes no descendant branch claim.",
]

ENERGY_PREP_LINES = [
    "For a future Timoshenko mode-shape energy diagnostic, meaningful eta-aware quantities are:",
    "U_b,i = 1/2 int E*I_i*(psi_i')^2 dx.",
    "U_s,i = 1/2 int kappa*G*A_i*(w_i' - psi_i)^2 dx.",
    "U_a,i = 1/2 int E*A_i*(u_i')^2 dx.",
    "For eta != 0, these must use the same tau-scaled A_i and I_i as the determinant.",
]


def finite_abs_rel(value: float, reference: float) -> tuple[float, float]:
    if not (math.isfinite(float(value)) and math.isfinite(float(reference))):
        return float("nan"), float("nan")
    abs_diff = abs(float(value) - float(reference))
    denom = max(abs(float(reference)), 1.0e-15)
    return abs_diff, abs_diff / denom


def status_from_tolerance(abs_diff: float, rel_diff: float, abs_tol: float, rel_tol: float) -> str:
    if not (math.isfinite(float(abs_diff)) and math.isfinite(float(rel_diff))):
        return "fail"
    return "pass" if abs_diff <= abs_tol or rel_diff <= rel_tol else "fail"


def max_finite(values: Sequence[float]) -> float:
    finite = [abs(float(value)) for value in values if math.isfinite(float(value))]
    return max(finite) if finite else float("nan")


def format_float(value: object, digits: int = 6) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(number):
        return "nan"
    return f"{number:.{digits}e}"


def format_roots(values: Sequence[float]) -> str:
    return ", ".join("nan" if not math.isfinite(float(value)) else f"{float(value):.8g}" for value in values)


def classify_decrease(values: Sequence[float]) -> str:
    finite_values = [float(value) for value in values if math.isfinite(float(value))]
    if len(finite_values) != len(values):
        return "fail_missing_value"
    monotone = all(
        finite_values[index + 1] <= finite_values[index] * (1.0 + 1.0e-8) + 1.0e-12
        for index in range(len(finite_values) - 1)
    )
    if monotone:
        return "pass_monotone_decrease"
    final_is_smaller = finite_values[-1] <= 0.25 * finite_values[0]
    final_is_lowest = finite_values[-1] <= min(finite_values[:-1]) * (1.0 + 0.05)
    if final_is_smaller and final_is_lowest:
        return "pass_general_decrease"
    return "fail_no_clear_decrease"


def legacy_timo_coupling_matrix(
    Lambda: float,
    beta_deg: float,
    mu: float,
    epsilon: float,
) -> tuple[np.ndarray, tuple[str, ...]]:
    """Eta-zero reference path copied from the original variable-length helper."""
    section = section_from_epsilon(float(epsilon))
    basis = timo_basis(float(Lambda), float(epsilon), section)
    theta = project_omega(float(Lambda), float(epsilon))
    l1, l2 = segment_lengths(float(mu))
    rod1 = endpoint_columns(l1, theta, basis, section)
    rod2 = endpoint_columns(-l2, theta, basis, section)
    cb = math.cos(math.radians(float(beta_deg)))
    sb = math.sin(math.radians(float(beta_deg)))
    matrix = np.zeros((6, 6), dtype=float)
    matrix[0, 0:3] = rod1["w"]
    matrix[0, 3:6] = -cb * rod2["w"] + sb * rod2["u"]
    matrix[1, 0:3] = rod1["u"]
    matrix[1, 3:6] = -sb * rod2["w"] - cb * rod2["u"]
    matrix[2, 0:3] = rod1["psi"]
    matrix[2, 3:6] = -rod2["psi"]
    matrix[3, 0:3] = rod1["M"]
    matrix[3, 3:6] = -rod2["M"]
    matrix[4, 0:3] = -rod1["Q"]
    matrix[4, 3:6] = cb * rod2["Q"] - sb * rod2["N"]
    matrix[5, 0:3] = rod1["N"]
    matrix[5, 3:6] = -sb * rod2["Q"] - cb * rod2["N"]
    return matrix, basis.warnings


def legacy_timo_det_for_scan(Lambda: float, beta_deg: float, mu: float, epsilon: float) -> float:
    try:
        matrix, _warnings = legacy_timo_coupling_matrix(
            float(Lambda),
            float(beta_deg),
            float(mu),
            float(epsilon),
        )
    except (FloatingPointError, ValueError, OverflowError):
        return float("nan")
    if not np.all(np.isfinite(matrix)):
        return float("nan")
    return normalized_det(matrix)


def legacy_timo_sorted_roots(
    beta_deg: float,
    mu: float,
    epsilon: float,
    n_roots: int,
) -> tuple[np.ndarray, list[str]]:
    l1, _l2 = segment_lengths(float(mu))
    upper = max(18.0, 1.35 * math.pi * (float(n_roots) + 4) / max(l1, 0.08))
    return find_roots_by_sign_scan(
        lambda value: legacy_timo_det_for_scan(value, float(beta_deg), float(mu), float(epsilon)),
        n_roots,
        start=0.2,
        upper=upper,
        scan_step=0.01,
    )


def eb_eta_roots(beta_deg: float, mu: float, epsilon: float, eta: float) -> np.ndarray:
    return find_first_n_roots_eta(
        beta=float(np.deg2rad(float(beta_deg))),
        mu=float(mu),
        epsilon=float(epsilon),
        eta=float(eta),
        n_roots=N_MODES,
        Lmin=0.2,
        Lmax0=45.0,
        scan_step=0.005,
        grow_factor=1.35,
        max_tries=8,
    )


def append_compare_row(
    rows: list[dict[str, object]],
    *,
    test_name: str,
    beta: float,
    mu: float,
    eta: float,
    mu_reference: float | str,
    eta_reference: float | str,
    epsilon: float,
    mode: int,
    value: float,
    reference: float,
    eb_value: float | str = "",
    status: str,
    notes: str,
) -> None:
    abs_diff, rel_diff = finite_abs_rel(value, reference)
    rows.append(
        {
            "test_name": test_name,
            "beta": float(beta),
            "mu": float(mu),
            "eta": float(eta),
            "mu_reference": mu_reference,
            "eta_reference": eta_reference,
            "epsilon": float(epsilon),
            "mode": int(mode),
            "Lambda_Timo": float(value),
            "Lambda_reference": float(reference),
            "Lambda_EB": eb_value,
            "abs_diff": abs_diff,
            "rel_diff": rel_diff,
            "status": status,
            "notes": notes,
        }
    )


def compute_formula_audit() -> dict[str, float]:
    mass_errors: list[float] = []
    tau_errors: list[float] = []
    cutoff_errors: list[float] = []
    for mu in [-0.6, -0.3, 0.0, 0.3, 0.6]:
        for eta in [-0.5, -1.0e-3, 0.0, 1.0e-3, 0.5]:
            local = tau_factors(mu, eta)
            reference = thickness_mismatch_factors(mu, eta)
            mass_errors.append(abs(local.mass_factor - 2.0))
            tau_errors.append(abs(local.tau1 - reference.tau1))
            tau_errors.append(abs(local.tau2 - reference.tau2))
            for epsilon in [0.01, 0.025]:
                section1 = section_from_epsilon_tau(epsilon, local.tau1)
                section2 = section_from_epsilon_tau(epsilon, local.tau2)
                base = section_from_epsilon(epsilon)
                expected1 = lambda_cutoff(epsilon, base) / math.sqrt(local.tau1)
                expected2 = lambda_cutoff(epsilon, base) / math.sqrt(local.tau2)
                cutoff_errors.append(abs(lambda_cutoff(epsilon, section1) - expected1))
                cutoff_errors.append(abs(lambda_cutoff(epsilon, section2) - expected2))
    return {
        "max_mass_error": max_finite(mass_errors),
        "max_tau_helper_abs_diff": max_finite(tau_errors),
        "max_cutoff_scaling_abs_diff": max_finite(cutoff_errors),
    }


def compute_test_e(rows: list[dict[str, object]], warnings: list[str]) -> dict[str, object]:
    fail_count = 0
    abs_values: list[float] = []
    rel_values: list[float] = []
    worst: dict[str, object] | None = None
    for beta in E_BETAS:
        for mu in E_MUS:
            for epsilon in E_EPSILONS:
                legacy, legacy_warnings = legacy_timo_sorted_roots(beta, mu, epsilon, N_MODES)
                current, current_warnings = timo_sorted_roots(beta, mu, epsilon, N_MODES, eta=0.0)
                warnings.extend(
                    f"test E beta={beta:g}, mu={mu:g}, epsilon={epsilon:g}: legacy: {item}"
                    for item in legacy_warnings
                )
                warnings.extend(
                    f"test E beta={beta:g}, mu={mu:g}, epsilon={epsilon:g}: current: {item}"
                    for item in current_warnings
                )
                for mode_index, (left, right) in enumerate(zip(current, legacy, strict=True), start=1):
                    abs_diff, rel_diff = finite_abs_rel(float(left), float(right))
                    status = status_from_tolerance(abs_diff, rel_diff, E_ABS_TOL, E_REL_TOL)
                    if status != "pass":
                        fail_count += 1
                    abs_values.append(abs_diff)
                    rel_values.append(rel_diff)
                    if worst is None or (math.isfinite(abs_diff) and abs_diff > float(worst["abs_diff"])):
                        worst = {
                            "beta": beta,
                            "mu": mu,
                            "epsilon": epsilon,
                            "mode": mode_index,
                            "abs_diff": abs_diff,
                            "rel_diff": rel_diff,
                        }
                    append_compare_row(
                        rows,
                        test_name="E_eta0_regression_against_legacy_variable_length_timo",
                        beta=beta,
                        mu=mu,
                        eta=0.0,
                        mu_reference=mu,
                        eta_reference=0.0,
                        epsilon=epsilon,
                        mode=mode_index,
                        value=float(left),
                        reference=float(right),
                        status=status,
                        notes=f"abs_tol={E_ABS_TOL:g}; rel_tol={E_REL_TOL:g}",
                    )
    return {
        "status": "pass" if fail_count == 0 else "fail",
        "fail_count": fail_count,
        "max_abs": max_finite(abs_values),
        "max_rel": max_finite(rel_values),
        "worst": worst or {},
    }


def compute_test_f(rows: list[dict[str, object]], warnings: list[str]) -> dict[str, object]:
    fail_count = 0
    continuity_by_abs_eta: dict[float, dict[str, float]] = {}
    diff_by_key: dict[tuple[float, float, float, float, int], float] = {}
    mu0_plus_minus_abs: list[float] = []
    all_abs: list[float] = []
    all_rel: list[float] = []

    for beta in F_BETAS:
        for mu in F_MUS:
            for epsilon in F_EPSILONS:
                baseline, baseline_warnings = timo_sorted_roots(beta, mu, epsilon, N_MODES, eta=0.0)
                warnings.extend(
                    f"test F beta={beta:g}, mu={mu:g}, epsilon={epsilon:g}: eta=0: {item}"
                    for item in baseline_warnings
                )
                roots_by_eta: dict[float, np.ndarray] = {}
                for eta in F_ETAS:
                    roots, root_warnings = timo_sorted_roots(beta, mu, epsilon, N_MODES, eta=eta)
                    roots_by_eta[eta] = roots
                    warnings.extend(
                        f"test F beta={beta:g}, mu={mu:g}, epsilon={epsilon:g}, eta={eta:g}: {item}"
                        for item in root_warnings
                    )
                    abs_eta = abs(float(eta))
                    continuity_by_abs_eta.setdefault(abs_eta, {"max_abs": 0.0, "max_rel": 0.0})
                    for mode_index, (value, reference) in enumerate(zip(roots, baseline, strict=True), start=1):
                        abs_diff, rel_diff = finite_abs_rel(float(value), float(reference))
                        all_abs.append(abs_diff)
                        all_rel.append(rel_diff)
                        continuity_by_abs_eta[abs_eta]["max_abs"] = max(
                            continuity_by_abs_eta[abs_eta]["max_abs"],
                            abs_diff if math.isfinite(abs_diff) else float("inf"),
                        )
                        continuity_by_abs_eta[abs_eta]["max_rel"] = max(
                            continuity_by_abs_eta[abs_eta]["max_rel"],
                            rel_diff if math.isfinite(rel_diff) else float("inf"),
                        )
                        diff_by_key[(beta, mu, epsilon, eta, mode_index)] = abs_diff
                        status = "pass" if math.isfinite(abs_diff) and math.isfinite(rel_diff) else "fail"
                        if status != "pass":
                            fail_count += 1
                        append_compare_row(
                            rows,
                            test_name="F_eta_to_zero_continuity",
                            beta=beta,
                            mu=mu,
                            eta=eta,
                            mu_reference=mu,
                            eta_reference=0.0,
                            epsilon=epsilon,
                            mode=mode_index,
                            value=float(value),
                            reference=float(reference),
                            status=status,
                            notes="continuity trend checked by |eta| groups",
                        )
                if mu == 0.0:
                    for abs_eta in [1.0e-4, 1.0e-3]:
                        plus = roots_by_eta[abs_eta]
                        minus = roots_by_eta[-abs_eta]
                        mu0_plus_minus_abs.extend(abs(float(a) - float(b)) for a, b in zip(plus, minus, strict=True))

    small = continuity_by_abs_eta.get(1.0e-4, {"max_abs": float("nan"), "max_rel": float("nan")})
    large = continuity_by_abs_eta.get(1.0e-3, {"max_abs": float("nan"), "max_rel": float("nan")})
    abs_trend_ok = (
        math.isfinite(small["max_abs"])
        and math.isfinite(large["max_abs"])
        and small["max_abs"] <= 0.25 * large["max_abs"] + 1.0e-9
    )
    rel_trend_ok = (
        math.isfinite(small["max_rel"])
        and math.isfinite(large["max_rel"])
        and small["max_rel"] <= 0.25 * large["max_rel"] + 1.0e-9
    )
    mu0_symmetry = max_finite(mu0_plus_minus_abs)
    mu0_symmetry_ok = math.isfinite(mu0_symmetry) and mu0_symmetry <= 2.0e-6
    status = "pass" if fail_count == 0 and abs_trend_ok and rel_trend_ok and mu0_symmetry_ok else "fail"
    return {
        "status": status,
        "fail_count": fail_count + (0 if abs_trend_ok and rel_trend_ok and mu0_symmetry_ok else 1),
        "by_abs_eta": continuity_by_abs_eta,
        "max_abs": max_finite(all_abs),
        "max_rel": max_finite(all_rel),
        "abs_trend_ok": abs_trend_ok,
        "rel_trend_ok": rel_trend_ok,
        "mu0_plus_minus_max_abs": mu0_symmetry,
        "mu0_plus_minus_status": "pass" if mu0_symmetry_ok else "fail",
    }


def compute_test_g(rows: list[dict[str, object]], warnings: list[str]) -> dict[str, object]:
    fail_count = 0
    abs_values: list[float] = []
    rel_values: list[float] = []
    summaries: list[dict[str, object]] = []
    for beta in G_BETAS:
        for epsilon in G_EPSILONS:
            for mu, eta in G_PAIRS:
                left, left_warnings = timo_sorted_roots(beta, mu, epsilon, N_MODES, eta=eta)
                right, right_warnings = timo_sorted_roots(beta, -mu, epsilon, N_MODES, eta=-eta)
                warnings.extend(
                    f"test G beta={beta:g}, epsilon={epsilon:g}, mu={mu:g}, eta={eta:g}: left: {item}"
                    for item in left_warnings
                )
                warnings.extend(
                    f"test G beta={beta:g}, epsilon={epsilon:g}, mu={mu:g}, eta={eta:g}: right: {item}"
                    for item in right_warnings
                )
                case_abs: list[float] = []
                case_rel: list[float] = []
                for mode_index, (value, reference) in enumerate(zip(left, right, strict=True), start=1):
                    abs_diff, rel_diff = finite_abs_rel(float(value), float(reference))
                    status = status_from_tolerance(abs_diff, rel_diff, G_ABS_TOL, G_REL_TOL)
                    if status != "pass":
                        fail_count += 1
                    abs_values.append(abs_diff)
                    rel_values.append(rel_diff)
                    case_abs.append(abs_diff)
                    case_rel.append(rel_diff)
                    append_compare_row(
                        rows,
                        test_name="G_swap_symmetry",
                        beta=beta,
                        mu=mu,
                        eta=eta,
                        mu_reference=-mu,
                        eta_reference=-eta,
                        epsilon=epsilon,
                        mode=mode_index,
                        value=float(value),
                        reference=float(reference),
                        status=status,
                        notes=f"compare (mu,eta) to (-mu,-eta); abs_tol={G_ABS_TOL:g}; rel_tol={G_REL_TOL:g}",
                    )
                summaries.append(
                    {
                        "beta": beta,
                        "epsilon": epsilon,
                        "mu": mu,
                        "eta": eta,
                        "max_abs": max_finite(case_abs),
                        "max_rel": max_finite(case_rel),
                    }
                )
    return {
        "status": "pass" if fail_count == 0 else "fail",
        "fail_count": fail_count,
        "max_abs": max_finite(abs_values),
        "max_rel": max_finite(rel_values),
        "summaries": summaries,
    }


def compute_test_h(rows: list[dict[str, object]], warnings: list[str]) -> dict[str, object]:
    fail_count = 0
    summaries: list[dict[str, object]] = []
    for beta in H_BETAS:
        for mu in H_MUS:
            for eta in H_ETAS:
                case_rows: list[dict[str, object]] = []
                max_rel_by_epsilon: list[float] = []
                max_abs_by_epsilon: list[float] = []
                for epsilon in H_EPSILONS:
                    eb_values = eb_eta_roots(beta, mu, epsilon, eta)
                    timo_values, root_warnings = timo_sorted_roots(beta, mu, epsilon, N_MODES, eta=eta)
                    warnings.extend(
                        f"test H beta={beta:g}, mu={mu:g}, eta={eta:g}, epsilon={epsilon:g}: {item}"
                        for item in root_warnings
                    )
                    rel_values: list[float] = []
                    abs_values: list[float] = []
                    for mode_index, (timo_value, eb_value) in enumerate(
                        zip(timo_values, eb_values, strict=True),
                        start=1,
                    ):
                        abs_diff, rel_diff = finite_abs_rel(float(timo_value), float(eb_value))
                        abs_values.append(abs_diff)
                        rel_values.append(rel_diff)
                        case_rows.append(
                            {
                                "test_name": "H_epsilon_to_zero_EB_thickness_mismatch_limit",
                                "beta": float(beta),
                                "mu": float(mu),
                                "eta": float(eta),
                                "mu_reference": float(mu),
                                "eta_reference": float(eta),
                                "epsilon": float(epsilon),
                                "mode": mode_index,
                                "Lambda_Timo": float(timo_value),
                                "Lambda_reference": float(eb_value),
                                "Lambda_EB": float(eb_value),
                                "abs_diff": abs_diff,
                                "rel_diff": rel_diff,
                                "status": "pending",
                                "notes": "sorted roots only; trend checked across epsilon",
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
                        "beta": beta,
                        "mu": mu,
                        "eta": eta,
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


def compute_test_i(rows: list[dict[str, object]], warnings: list[str]) -> dict[str, object]:
    root_cache: dict[tuple[float, float, float], np.ndarray] = {}
    root_table: list[dict[str, object]] = []
    for epsilon in I_EPSILONS:
        for eta in I_ETAS:
            for mu in I_MUS:
                roots, root_warnings = timo_sorted_roots(0.0, mu, epsilon, N_MODES, eta=eta)
                root_cache[(epsilon, mu, eta)] = roots
                warnings.extend(
                    f"test I beta=0, epsilon={epsilon:g}, mu={mu:g}, eta={eta:g}: {item}"
                    for item in root_warnings
                )
                root_table.append({"epsilon": epsilon, "mu": mu, "eta": eta, "roots": tuple(float(x) for x in roots)})
                for mode_index, value in enumerate(roots, start=1):
                    rows.append(
                        {
                            "test_name": "I_beta0_composite_roots",
                            "beta": 0.0,
                            "mu": float(mu),
                            "eta": float(eta),
                            "mu_reference": "",
                            "eta_reference": "",
                            "epsilon": float(epsilon),
                            "mode": mode_index,
                            "Lambda_Timo": float(value),
                            "Lambda_reference": "",
                            "Lambda_EB": "",
                            "abs_diff": "",
                            "rel_diff": "",
                            "status": "info",
                            "notes": "eta!=0 straight composite rod; mu-invariance is not expected",
                        }
                    )

    swap_abs: list[float] = []
    swap_rel: list[float] = []
    swap_fail_count = 0
    for epsilon in I_EPSILONS:
        for eta in I_ETAS:
            for mu in [0.0, 0.3, 0.6]:
                left = root_cache[(epsilon, mu, eta)]
                right = root_cache[(epsilon, -mu, -eta)]
                for mode_index, (value, reference) in enumerate(zip(left, right, strict=True), start=1):
                    abs_diff, rel_diff = finite_abs_rel(float(value), float(reference))
                    status = status_from_tolerance(abs_diff, rel_diff, I_SWAP_ABS_TOL, I_SWAP_REL_TOL)
                    if status != "pass":
                        swap_fail_count += 1
                    swap_abs.append(abs_diff)
                    swap_rel.append(rel_diff)
                    append_compare_row(
                        rows,
                        test_name="I_beta0_swap_symmetry",
                        beta=0.0,
                        mu=mu,
                        eta=eta,
                        mu_reference=-mu,
                        eta_reference=-eta,
                        epsilon=epsilon,
                        mode=mode_index,
                        value=float(value),
                        reference=float(reference),
                        status=status,
                        notes=(
                            "beta=0 straight composite swap check; "
                            f"abs_tol={I_SWAP_ABS_TOL:g}; rel_tol={I_SWAP_REL_TOL:g}"
                        ),
                    )

    eta0_abs: list[float] = []
    eta0_rel: list[float] = []
    eta0_fail_count = 0
    for epsilon in I_EPSILONS:
        baseline, baseline_warnings = timo_sorted_roots(0.0, 0.0, epsilon, N_MODES, eta=0.0)
        warnings.extend(f"test I beta=0 eta=0, epsilon={epsilon:g}: baseline: {item}" for item in baseline_warnings)
        for mu in I_MUS:
            roots, root_warnings = timo_sorted_roots(0.0, mu, epsilon, N_MODES, eta=0.0)
            warnings.extend(f"test I beta=0 eta=0, epsilon={epsilon:g}, mu={mu:g}: {item}" for item in root_warnings)
            for mode_index, (value, reference) in enumerate(zip(roots, baseline, strict=True), start=1):
                abs_diff, rel_diff = finite_abs_rel(float(value), float(reference))
                status = status_from_tolerance(abs_diff, rel_diff, I_ETA0_ABS_TOL, I_ETA0_REL_TOL)
                if status != "pass":
                    eta0_fail_count += 1
                eta0_abs.append(abs_diff)
                eta0_rel.append(rel_diff)
                append_compare_row(
                    rows,
                    test_name="I_beta0_eta0_mu_invariance_regression",
                    beta=0.0,
                    mu=mu,
                    eta=0.0,
                    mu_reference=0.0,
                    eta_reference=0.0,
                    epsilon=epsilon,
                    mode=mode_index,
                    value=float(value),
                    reference=float(reference),
                    status=status,
                    notes=(
                        "eta=0 straight rod mu-invariance regression; "
                        f"abs_tol={I_ETA0_ABS_TOL:g}; rel_tol={I_ETA0_REL_TOL:g}"
                    ),
                )

    status = "pass" if swap_fail_count == 0 and eta0_fail_count == 0 else "fail"
    return {
        "status": status,
        "fail_count": swap_fail_count + eta0_fail_count,
        "swap_fail_count": swap_fail_count,
        "eta0_fail_count": eta0_fail_count,
        "swap_max_abs": max_finite(swap_abs),
        "swap_max_rel": max_finite(swap_rel),
        "eta0_max_abs": max_finite(eta0_abs),
        "eta0_max_rel": max_finite(eta0_rel),
        "root_table": root_table,
    }


def write_csv(rows: Sequence[dict[str, object]]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def test_f_table(summary: dict[str, object]) -> list[str]:
    by_abs_eta = summary["by_abs_eta"]
    lines = [
        "| |eta| | max abs difference from eta=0 | max relative difference from eta=0 |",
        "| ---: | ---: | ---: |",
    ]
    for abs_eta in sorted(by_abs_eta):
        item = by_abs_eta[abs_eta]
        lines.append(
            f"| {abs_eta:.0e} | {format_float(item['max_abs'])} | {format_float(item['max_rel'])} |"
        )
    return lines


def test_g_table(summaries: Sequence[dict[str, object]]) -> list[str]:
    lines = [
        "| beta deg | epsilon | mu | eta | max abs | max rel |",
        "| ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for item in summaries:
        lines.append(
            "| "
            f"{float(item['beta']):g} | "
            f"{float(item['epsilon']):g} | "
            f"{float(item['mu']):g} | "
            f"{float(item['eta']):g} | "
            f"{format_float(item['max_abs'])} | "
            f"{format_float(item['max_rel'])} |"
        )
    return lines


def test_h_table(summaries: Sequence[dict[str, object]]) -> list[str]:
    lines = [
        "| beta deg | mu | eta | max rel @ eps=0.02 | max rel @ eps=0.01 | max rel @ eps=0.005 | max rel @ eps=0.0025 | trend | status |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    for item in summaries:
        rels = list(item["max_rel_by_epsilon"])
        lines.append(
            "| "
            f"{float(item['beta']):g} | "
            f"{float(item['mu']):g} | "
            f"{float(item['eta']):g} | "
            f"{format_float(rels[0])} | "
            f"{format_float(rels[1])} | "
            f"{format_float(rels[2])} | "
            f"{format_float(rels[3])} | "
            f"{item['trend']} | "
            f"{item['status']} |"
        )
    return lines


def test_i_roots_table(root_table: Sequence[dict[str, object]]) -> list[str]:
    lines = [
        "| epsilon | eta | mu | first 6 sorted Timoshenko roots |",
        "| ---: | ---: | ---: | --- |",
    ]
    for item in root_table:
        lines.append(
            "| "
            f"{float(item['epsilon']):g} | "
            f"{float(item['eta']):g} | "
            f"{float(item['mu']):g} | "
            f"{format_roots(item['roots'])} |"
        )
    return lines


def warning_lines(warnings: Sequence[str]) -> list[str]:
    unique = sorted(set(warnings))
    if not unique:
        return ["No root-search or basis warnings were recorded."]
    lines = ["Warnings:"]
    for warning in unique[:60]:
        lines.append(f"- {warning}")
    if len(unique) > 60:
        lines.append(f"- ... {len(unique) - 60} more warnings")
    return lines


def write_report(
    formula_audit: dict[str, float],
    test_e: dict[str, object],
    test_f: dict[str, object],
    test_g: dict[str, object],
    test_h: dict[str, object],
    test_i: dict[str, object],
    warnings: Sequence[str],
) -> None:
    overall_pass = all(
        item["status"] == "pass"
        for item in [test_e, test_f, test_g, test_h, test_i]
    )
    verification_line = (
        "The tau-aware eta!=0 variable-length Timoshenko implementation is diagnostic-verified for these sorted-root checks."
        if overall_pass
        else "The tau-aware eta!=0 variable-length Timoshenko implementation is not diagnostic-verified by this run."
    )
    worst_e = test_e.get("worst") or {}
    lines: list[str] = [
        "# Variable-Length Thickness-Mismatch Timoshenko Limits Audit",
        "",
        "## Purpose",
        "",
        "Diagnostic-only sorted-root verification of the tau-aware variable-length Timoshenko implementation. "
        "This audit does not change article files, article figures, the old determinant, `formulas.py`, old "
        "solvers, baseline results, the existing FEM physical model, or any Gmsh/CalculiX workflow.",
        "",
        "The checks below are determinant-limit and symmetry checks only. They are not descendant branch tracking.",
        "",
        "## Implemented Formula Audit",
        "",
        "```text",
        *FORMULA_AUDIT_LINES,
        "```",
        "",
        "- max mass-conservation error: " + format_float(formula_audit["max_mass_error"]),
        "- max tau-helper difference versus `formulas_thickness_mismatch.py`: "
        + format_float(formula_audit["max_tau_helper_abs_diff"]),
        "- max cut-off scaling residual: " + format_float(formula_audit["max_cutoff_scaling_abs_diff"]),
        "",
        "## Test E: eta=0 Regression",
        "",
        "Comparison: the new tau-aware matrix at `eta=0` versus the copied legacy eta=0 variable-length "
        "Timoshenko matrix for beta values `0, 15, 45`, mu values `0, 0.3, 0.6`, epsilon values "
        "`0.0025, 0.01, 0.025`, and the first six sorted roots.",
        "",
        f"- status: {test_e['status']}",
        f"- failures: {test_e['fail_count']}",
        f"- max abs diff: {format_float(test_e['max_abs'])}",
        f"- max relative diff: {format_float(test_e['max_rel'])}",
    ]
    if worst_e:
        lines.append(
            "- worst row: "
            f"beta={float(worst_e['beta']):g}, mu={float(worst_e['mu']):g}, "
            f"epsilon={float(worst_e['epsilon']):g}, mode={int(worst_e['mode'])}, "
            f"abs={format_float(worst_e['abs_diff'])}, rel={format_float(worst_e['rel_diff'])}"
        )
    lines.extend(
        [
            "",
            "## Test F: eta -> 0 Continuity",
            "",
            "Comparison: small nonzero eta values against eta=0 at beta `15 deg`, mu values `0, 0.3, 0.6`, "
            "epsilon values `0.01, 0.025`, and the first six sorted roots.",
            "",
            *test_f_table(test_f),
            "",
            f"- status: {test_f['status']}",
            f"- finite-row failures: {test_f['fail_count']}",
            f"- |eta| trend abs check: {test_f['abs_trend_ok']}",
            f"- |eta| trend relative check: {test_f['rel_trend_ok']}",
            f"- mu=0 plus/minus eta max abs difference: {format_float(test_f['mu0_plus_minus_max_abs'])}",
            f"- mu=0 plus/minus eta symmetry status: {test_f['mu0_plus_minus_status']}",
            "",
            "For fixed positive mu, `+eta` and `-eta` are different physical thickness assignments, so only "
            "the `mu=0` plus/minus comparison is treated as a symmetry check here.",
            "",
            "## Test G: Swap Symmetry",
            "",
            "Comparison: `(mu, eta)` versus `(-mu, -eta)` at beta values `15, 45 deg`, epsilon values "
            "`0.01, 0.025`, and the first six sorted roots.",
            "",
            *test_g_table(test_g["summaries"]),
            "",
            f"- status: {test_g['status']}",
            f"- failures: {test_g['fail_count']}",
            f"- max abs diff: {format_float(test_g['max_abs'])}",
            f"- max relative diff: {format_float(test_g['max_rel'])}",
            "",
            "## Test H: epsilon -> 0 EB Thickness-Mismatch Limit",
            "",
            "Comparison: tau-aware Timoshenko sorted roots versus the existing EB thickness-mismatch determinant "
            "roots from `src/my_project/analytic/formulas_thickness_mismatch.py`.",
            "",
            *test_h_table(test_h["summaries"]),
            "",
            f"- status: {test_h['status']}",
            f"- failed trend cases: {test_h['fail_count']}",
            "",
            "This test compares sorted roots only. It does not make descendant branch claims.",
            "",
            "## Test I: beta=0 Straight Composite Rod",
            "",
            "For `beta=0` and `eta!=0`, the model is a straight two-segment rod with different radii and lengths. "
            "No mu-invariance is expected in general. The audit checks swap symmetry and also confirms that "
            "`eta=0` remains mu-invariant.",
            "",
            *test_i_roots_table(test_i["root_table"]),
            "",
            f"- status: {test_i['status']}",
            f"- swap failures: {test_i['swap_fail_count']}",
            f"- eta=0 mu-invariance failures: {test_i['eta0_fail_count']}",
            f"- beta=0 swap max abs diff: {format_float(test_i['swap_max_abs'])}",
            f"- beta=0 swap max relative diff: {format_float(test_i['swap_max_rel'])}",
            f"- eta=0 mu-invariance max abs diff: {format_float(test_i['eta0_max_abs'])}",
            f"- eta=0 mu-invariance max relative diff: {format_float(test_i['eta0_max_rel'])}",
            "",
            "## Energy Formula Preparation",
            "",
            "```text",
            *ENERGY_PREP_LINES,
            "```",
            "",
            "Energy plots were not run in this audit.",
            "",
            "## Overall Status",
            "",
            verification_line,
            "",
            "The eta=0 variable-length Timoshenko implementation remains diagnostic-verified by the earlier "
            "`results/variable_length_timoshenko_limits_audit.md` sorted-root limit checks. This audit extends "
            "the diagnostic verification to tau-aware eta!=0 sorted-root checks only when all tests E-I pass.",
            "",
            "## What Remains Unverified",
            "",
            "- Descendant branch identity and branch tracking for eta!=0.",
            "- Energy partition plots and energy-based localization claims.",
            "- FEM validation for eta!=0.",
            "- Symbolic proof of the tau-aware Timoshenko determinant.",
            "- High-frequency behavior near or above per-arm Timoshenko cut-off.",
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
    formula_audit = compute_formula_audit()
    test_e = compute_test_e(rows, warnings)
    test_f = compute_test_f(rows, warnings)
    test_g = compute_test_g(rows, warnings)
    test_h = compute_test_h(rows, warnings)
    test_i = compute_test_i(rows, warnings)
    write_csv(rows)
    write_report(formula_audit, test_e, test_f, test_g, test_h, test_i, warnings)

    overall_pass = all(item["status"] == "pass" for item in [test_e, test_f, test_g, test_h, test_i])
    print("Variable-length thickness-mismatch Timoshenko limits audit")
    print(f"Formula audit max mass error: {formula_audit['max_mass_error']:.6e}")
    print(f"Test E status: {test_e['status']}; max abs={test_e['max_abs']:.6e}; max rel={test_e['max_rel']:.6e}")
    print(
        "Test F status: "
        f"{test_f['status']}; max abs={test_f['max_abs']:.6e}; max rel={test_f['max_rel']:.6e}; "
        f"mu0 +/- max abs={test_f['mu0_plus_minus_max_abs']:.6e}"
    )
    print(f"Test G status: {test_g['status']}; max abs={test_g['max_abs']:.6e}; max rel={test_g['max_rel']:.6e}")
    print(f"Test H status: {test_h['status']}; failed trend cases={test_h['fail_count']}")
    print(
        "Test I status: "
        f"{test_i['status']}; swap max abs={test_i['swap_max_abs']:.6e}; "
        f"eta0 mu-invariance max abs={test_i['eta0_max_abs']:.6e}"
    )
    print(f"Overall status: {'pass' if overall_pass else 'fail'}")
    print(f"Wrote {OUTPUT_CSV.relative_to(REPO_ROOT)}")
    print(f"Wrote {OUTPUT_REPORT.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
