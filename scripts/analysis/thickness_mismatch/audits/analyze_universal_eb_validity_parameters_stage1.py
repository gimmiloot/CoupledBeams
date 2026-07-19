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
from scipy.stats import pearsonr, spearmanr


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_eb_validity_vs_timoshenko_stage1 as stage1,
)
from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_timoshenko_shape_construction as shape_audit,
)
from scripts.lib import variable_length_timoshenko as TIMO  # noqa: E402


DEFAULT_OUTPUT_DIR = Path("results") / "eb_validity_vs_timoshenko_stage1"
DEFAULT_N_SHAPE_POINTS = stage1.DEFAULT_SHAPE_SAMPLE_COUNT
PREDICTORS = ("epsilon_max", "chi_max", "chi_eff", "Pi_EB")
SUBGROUPS = ("all_physical", "bending_dominated", "longitudinal_mixed")
MAIN_THRESHOLD = 0.10
NUMERICAL_FLOOR = 1.0e-12

PI_FIELDS = [
    "Pi_shear",
    "Pi_rotary",
    "Pi_EB",
    "U_shear_star",
    "U_EB",
    "T_rot_star",
    "T_trans_star",
    "U_shear_star_rod1",
    "U_shear_star_rod2",
    "T_rot_star_rod1",
    "T_rot_star_rod2",
    "Pi_shear_rod1_fraction",
    "Pi_shear_rod2_fraction",
    "Pi_rotary_rod1_fraction",
    "Pi_rotary_rod2_fraction",
]

FIT_FIELDS = [
    "subgroup",
    "predictor",
    "row_count",
    "pearson_r",
    "pearson_p",
    "spearman_r",
    "spearman_p",
    "log_fit_count",
    "C",
    "p",
    "R2",
    "RMSE",
    "MAE",
    "threshold_10",
    "critical_count",
    "critical_mean",
    "critical_std",
    "critical_cv",
    "critical_cv_mu_means",
    "critical_cv_branch_means",
    "critical_cv_character_means",
    "notes",
]

CV_FIELDS = [
    "subgroup",
    "predictor",
    "held_out_mu",
    "train_count",
    "test_count",
    "C",
    "p",
    "threshold_10",
    "RMSE",
    "MAE",
    "R2",
    "false_safe_count",
    "false_unsafe_count",
    "false_safe_rate",
    "false_unsafe_rate",
    "balanced_accuracy",
    "notes",
]

CLASSIFICATION_FIELDS = [
    "subgroup",
    "predictor",
    "threshold_10",
    "row_count",
    "actual_safe_count",
    "actual_unsafe_count",
    "predicted_safe_count",
    "predicted_unsafe_count",
    "true_safe_count",
    "true_unsafe_count",
    "false_safe_count",
    "false_unsafe_count",
    "false_safe_rate",
    "false_unsafe_rate",
    "balanced_accuracy",
    "notes",
]


@dataclass(frozen=True)
class Args:
    output_dir: Path
    n_shape_points: int


def repo_path(path: Path) -> Path:
    path_obj = Path(path)
    return path_obj if path_obj.is_absolute() else REPO_ROOT / path_obj


def rel(path: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Post-process Stage-1 EB/Timoshenko data to compare universal applicability predictors.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-shape-points", type=int, default=DEFAULT_N_SHAPE_POINTS)
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    if int(ns.n_shape_points) < 51:
        raise ValueError("--n-shape-points must be at least 51")
    return Args(output_dir=repo_path(Path(ns.output_dir)), n_shape_points=int(ns.n_shape_points))


def fmt(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        value_f = float(value)
        if not isfinite(value_f):
            return "nan"
        return f"{value_f:.16e}"
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})
    return path


def finite_float(row: dict[str, object] | dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, float("nan")))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return float("nan")


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if abs(float(denominator)) > 1.0e-30 else float("nan")


def trapz(values: np.ndarray, coordinates: np.ndarray) -> float:
    return shape_audit.trapz(np.asarray(values, dtype=float), np.asarray(coordinates, dtype=float))


def eb_rod_arrays(
    result: shape_audit.ModeResult,
    *,
    rod_index: int,
    amplitude_scale: float,
) -> dict[str, np.ndarray]:
    coeff = np.asarray(result.coeff, dtype=float) * float(amplitude_scale)
    if rod_index == 1:
        A, B, _a2, _b2, P, _p2 = [float(value) for value in coeff]
        x = np.asarray(result.rod1.x_local, dtype=float)
    elif rod_index == 2:
        _a1, _b1, A, B, _p1, P = [float(value) for value in coeff]
        x = np.asarray(result.rod2.x_local, dtype=float)
    else:
        raise ValueError("rod_index must be 1 or 2")
    lam = float(result.Lambda)
    eps = float(result.epsilon)
    z = lam * x
    theta = eps * lam**2 * x
    u = P * np.sin(theta)
    w = A * (np.cos(z) - np.cosh(z)) + B * (np.sin(z) - np.sinh(z))
    u_prime = P * eps * lam**2 * np.cos(theta)
    w_prime = lam * (A * (-np.sin(z) - np.sinh(z)) + B * (np.cos(z) - np.cosh(z)))
    w_second = lam**2 * (A * (-np.cos(z) - np.cosh(z)) + B * (-np.sin(z) - np.sinh(z)))
    w_third = lam**3 * (A * (np.sin(z) - np.sinh(z)) + B * (-np.cos(z) - np.cosh(z)))
    return {
        "x": x,
        "u": u,
        "w": w,
        "u_prime": u_prime,
        "w_prime": w_prime,
        "w_second": w_second,
        "w_third": w_third,
    }


def eb_pi_metrics(result: shape_audit.ModeResult, *, amplitude_scale: float = 1.0) -> dict[str, float]:
    if result.model != stage1.MODEL_EB:
        raise ValueError("Pi_EB must be computed from an Euler-Bernoulli mode result")
    factors = TIMO.tau_factors(float(result.mu), float(result.eta))
    sections = {
        1: TIMO.section_from_epsilon_tau(float(result.epsilon), factors.tau1),
        2: TIMO.section_from_epsilon_tau(float(result.epsilon), factors.tau2),
    }
    omega = TIMO.project_omega(float(result.Lambda), float(result.epsilon))

    U_axial = 0.0
    U_bending = 0.0
    U_shear_star: dict[int, float] = {}
    T_trans: dict[int, float] = {}
    T_rot: dict[int, float] = {}

    for rod_index in (1, 2):
        arrays = eb_rod_arrays(result, rod_index=rod_index, amplitude_scale=float(amplitude_scale))
        section = sections[rod_index]
        coordinate = np.abs(np.asarray(arrays["x"], dtype=float))
        u = np.asarray(arrays["u"], dtype=float)
        w = np.asarray(arrays["w"], dtype=float)
        u_prime = np.asarray(arrays["u_prime"], dtype=float)
        w_prime = np.asarray(arrays["w_prime"], dtype=float)
        w_second = np.asarray(arrays["w_second"], dtype=float)
        w_third = np.asarray(arrays["w_third"], dtype=float)

        U_axial += 0.5 * TIMO.E * section.area * trapz(u_prime**2, coordinate)
        U_bending += 0.5 * section.bending_stiffness * trapz(w_second**2, coordinate)
        q_eb = -section.bending_stiffness * w_third
        U_shear_star[rod_index] = 0.5 * trapz(q_eb**2 / section.shear_stiffness, coordinate)
        T_trans[rod_index] = (
            0.5
            * omega**2
            * section.mass_per_length
            * trapz(u**2 + w**2, coordinate)
        )
        T_rot[rod_index] = (
            0.5
            * omega**2
            * section.rotary_inertia_per_length
            * trapz(w_prime**2, coordinate)
        )

    U_total = max(float(U_axial + U_bending), 0.0)
    U_shear_total = max(float(U_shear_star[1] + U_shear_star[2]), 0.0)
    T_trans_total = max(float(T_trans[1] + T_trans[2]), 0.0)
    T_rot_total = max(float(T_rot[1] + T_rot[2]), 0.0)
    pi_shear = safe_ratio(U_shear_total, U_total)
    pi_rotary = safe_ratio(T_rot_total, T_trans_total)
    return {
        "Pi_shear": pi_shear,
        "Pi_rotary": pi_rotary,
        "Pi_EB": pi_shear + pi_rotary if isfinite(pi_shear) and isfinite(pi_rotary) else float("nan"),
        "U_shear_star": U_shear_total,
        "U_EB": U_total,
        "T_rot_star": T_rot_total,
        "T_trans_star": T_trans_total,
        "U_shear_star_rod1": float(U_shear_star[1]),
        "U_shear_star_rod2": float(U_shear_star[2]),
        "T_rot_star_rod1": float(T_rot[1]),
        "T_rot_star_rod2": float(T_rot[2]),
        "Pi_shear_rod1_fraction": safe_ratio(float(U_shear_star[1]), U_shear_total),
        "Pi_shear_rod2_fraction": safe_ratio(float(U_shear_star[2]), U_shear_total),
        "Pi_rotary_rod1_fraction": safe_ratio(float(T_rot[1]), T_rot_total),
        "Pi_rotary_rod2_fraction": safe_ratio(float(T_rot[2]), T_rot_total),
    }


def eb_result_from_mode_row(row: dict[str, str], *, n_shape_points: int) -> shape_audit.ModeResult:
    return shape_audit.eb_mode_result(
        epsilon=finite_float(row, "epsilon_0"),
        beta_deg=finite_float(row, "beta_deg"),
        eta=finite_float(row, "eta"),
        mu=finite_float(row, "mu"),
        sorted_index=int(finite_float(row, "eb_sorted_index")),
        Lambda=finite_float(row, "Lambda_EB"),
        n_points=int(n_shape_points),
        case_role="universal_parameter_postprocess",
        root_source="stage1_mode_level_csv",
        root_warnings=(),
    )


def build_mode_metrics(output_dir: Path, *, n_shape_points: int) -> tuple[list[dict[str, object]], list[str]]:
    source_path = output_dir / "eb_timo_mode_level_metrics.csv"
    source_rows = read_csv(source_path)
    rows: list[dict[str, object]] = []
    warnings: list[str] = []
    for row in source_rows:
        if row.get("comparison_type") != "physical_branch":
            continue
        try:
            eb_result = eb_result_from_mode_row(row, n_shape_points=int(n_shape_points))
            pi = eb_pi_metrics(eb_result)
        except (FloatingPointError, ValueError, OverflowError) as exc:
            pi = {field: float("nan") for field in PI_FIELDS}
            warnings.append(
                "Pi_EB failed for "
                f"mu={row.get('mu')}, epsilon_0={row.get('epsilon_0')}, "
                f"branch={row.get('branch_id')}: {exc}"
            )
        out: dict[str, object] = dict(row)
        out.update(pi)
        rows.append(out)
    return rows, warnings


def subgroup_rows(rows: Sequence[dict[str, object]], subgroup: str) -> list[dict[str, object]]:
    if subgroup == "all_physical":
        return list(rows)
    if subgroup == "bending_dominated":
        return [row for row in rows if str(row.get("Timo_classification", "")) == "bending_dominated"]
    if subgroup == "longitudinal_mixed":
        return [row for row in rows if str(row.get("Timo_classification", "")) != "bending_dominated"]
    raise ValueError(f"unknown subgroup: {subgroup}")


def finite_xy(rows: Sequence[dict[str, object]], predictor: str) -> tuple[np.ndarray, np.ndarray]:
    x_values: list[float] = []
    y_values: list[float] = []
    for row in rows:
        x = finite_float(row, predictor)
        y = finite_float(row, "delta_f")
        if isfinite(x) and isfinite(y):
            x_values.append(float(x))
            y_values.append(float(y))
    return np.asarray(x_values, dtype=float), np.asarray(y_values, dtype=float)


def log_fit_arrays(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > NUMERICAL_FLOOR)
    return np.asarray(x[mask], dtype=float), np.asarray(y[mask], dtype=float)


def fit_power_law(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    xf, yf = log_fit_arrays(x, y)
    if len(xf) < 3:
        return {
            "log_fit_count": len(xf),
            "C": float("nan"),
            "p": float("nan"),
            "R2": float("nan"),
            "RMSE": float("nan"),
            "MAE": float("nan"),
            "threshold_10": float("nan"),
        }
    logx = np.log(xf)
    logy = np.log(yf)
    p, logc = np.polyfit(logx, logy, 1)
    pred_log = p * logx + logc
    c = float(np.exp(logc))
    pred = c * xf**p
    ss_res = float(np.sum((logy - pred_log) ** 2))
    ss_tot = float(np.sum((logy - float(np.mean(logy))) ** 2))
    threshold = float((MAIN_THRESHOLD / c) ** (1.0 / p)) if c > 0.0 and p > 0.0 else float("nan")
    return {
        "log_fit_count": len(xf),
        "C": c,
        "p": float(p),
        "R2": 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan"),
        "RMSE": float(np.sqrt(np.mean((pred - yf) ** 2))),
        "MAE": float(np.mean(np.abs(pred - yf))),
        "threshold_10": threshold,
    }


def correlation_values(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    if len(x) < 3 or float(np.std(x)) <= 0.0 or float(np.std(y)) <= 0.0:
        return {
            "pearson_r": float("nan"),
            "pearson_p": float("nan"),
            "spearman_r": float("nan"),
            "spearman_p": float("nan"),
        }
    pearson = pearsonr(x, y)
    spearman = spearmanr(x, y)
    return {
        "pearson_r": float(pearson.statistic),
        "pearson_p": float(pearson.pvalue),
        "spearman_r": float(spearman.statistic),
        "spearman_p": float(spearman.pvalue),
    }


def threshold_classification(
    rows: Sequence[dict[str, object]],
    predictor: str,
    threshold: float,
) -> dict[str, float | int | str]:
    if not isfinite(float(threshold)):
        return {
            "threshold_10": float("nan"),
            "row_count": 0,
            "actual_safe_count": 0,
            "actual_unsafe_count": 0,
            "predicted_safe_count": 0,
            "predicted_unsafe_count": 0,
            "true_safe_count": 0,
            "true_unsafe_count": 0,
            "false_safe_count": 0,
            "false_unsafe_count": 0,
            "false_safe_rate": float("nan"),
            "false_unsafe_rate": float("nan"),
            "balanced_accuracy": float("nan"),
            "notes": "non-finite threshold",
        }
    finite_rows = [
        row
        for row in rows
        if isfinite(finite_float(row, predictor)) and isfinite(finite_float(row, "delta_f"))
    ]
    actual_safe = [row for row in finite_rows if finite_float(row, "delta_f") <= MAIN_THRESHOLD]
    actual_unsafe = [row for row in finite_rows if finite_float(row, "delta_f") > MAIN_THRESHOLD]
    predicted_safe = [row for row in finite_rows if finite_float(row, predictor) <= float(threshold)]
    predicted_unsafe = [row for row in finite_rows if finite_float(row, predictor) > float(threshold)]
    true_safe = [
        row
        for row in finite_rows
        if finite_float(row, "delta_f") <= MAIN_THRESHOLD and finite_float(row, predictor) <= float(threshold)
    ]
    true_unsafe = [
        row
        for row in finite_rows
        if finite_float(row, "delta_f") > MAIN_THRESHOLD and finite_float(row, predictor) > float(threshold)
    ]
    false_safe = [
        row
        for row in finite_rows
        if finite_float(row, "delta_f") > MAIN_THRESHOLD and finite_float(row, predictor) <= float(threshold)
    ]
    false_unsafe = [
        row
        for row in finite_rows
        if finite_float(row, "delta_f") <= MAIN_THRESHOLD and finite_float(row, predictor) > float(threshold)
    ]
    safe_rate = safe_ratio(len(true_safe), len(actual_safe))
    unsafe_rate = safe_ratio(len(true_unsafe), len(actual_unsafe))
    return {
        "threshold_10": float(threshold),
        "row_count": len(finite_rows),
        "actual_safe_count": len(actual_safe),
        "actual_unsafe_count": len(actual_unsafe),
        "predicted_safe_count": len(predicted_safe),
        "predicted_unsafe_count": len(predicted_unsafe),
        "true_safe_count": len(true_safe),
        "true_unsafe_count": len(true_unsafe),
        "false_safe_count": len(false_safe),
        "false_unsafe_count": len(false_unsafe),
        "false_safe_rate": safe_ratio(len(false_safe), len(actual_unsafe)),
        "false_unsafe_rate": safe_ratio(len(false_unsafe), len(actual_safe)),
        "balanced_accuracy": 0.5 * (safe_rate + unsafe_rate) if isfinite(safe_rate) and isfinite(unsafe_rate) else float("nan"),
        "notes": "safe means predictor <= inferred 10% threshold",
    }


def coefficient_of_variation(values: Sequence[float]) -> float:
    finite = np.asarray([float(value) for value in values if isfinite(float(value))], dtype=float)
    if finite.size < 2:
        return float("nan")
    mean = float(np.mean(finite))
    if abs(mean) <= 1.0e-30:
        return float("nan")
    return float(np.std(finite, ddof=1) / abs(mean))


def grouped_cv(rows: Sequence[dict[str, object]], predictor: str, group_key: str) -> float:
    grouped: dict[str, list[float]] = {}
    for row in rows:
        value = finite_float(row, predictor)
        if isfinite(value):
            grouped.setdefault(str(row.get(group_key, "")), []).append(value)
    means = [float(np.mean(values)) for values in grouped.values() if values]
    return coefficient_of_variation(means)


def stage1_args_for_cache(output_dir: Path, *, n_shape_points: int) -> stage1.Args:
    return stage1.Args(
        beta_deg=stage1.DEFAULT_BETA_DEG,
        eta=stage1.DEFAULT_ETA,
        mu_values=tuple(stage1.DEFAULT_MU_VALUES),
        epsilon_values=tuple(stage1.DEFAULT_EPSILON_VALUES),
        n_reported_modes=stage1.DEFAULT_N_REPORTED_MODES,
        n_candidate_roots=stage1.DEFAULT_N_CANDIDATE_ROOTS,
        n_shape_points=int(n_shape_points),
        output_dir=Path(output_dir),
        cache_dir=Path(output_dir) / stage1.DEFAULT_CACHE_SUBDIR,
        reuse_cache=True,
        force_recompute=False,
        plot_only=False,
        skip_critical_refinement=False,
        smoke=False,
        max_refinement_iterations=stage1.DEFAULT_MAX_REFINEMENT_ITERATIONS,
    )


def critical_rows_with_pi(
    output_dir: Path,
    *,
    n_shape_points: int,
) -> tuple[list[dict[str, object]], list[str]]:
    critical_path = output_dir / "eb_timo_critical_thickness_by_branch.csv"
    critical_rows = read_csv(critical_path)
    cache = stage1.RootCache(stage1_args_for_cache(output_dir, n_shape_points=int(n_shape_points)))
    rows: list[dict[str, object]] = []
    warnings: list[str] = []
    for row in critical_rows:
        if row.get("critical_status") not in {"bracketed", "below_scan_min"}:
            continue
        epsilon = finite_float(row, "epsilon_0_crit_10")
        mu = finite_float(row, "mu")
        sorted_index = finite_float(row, "EB_sorted_index_at_crit")
        if not (isfinite(epsilon) and isfinite(mu) and isfinite(sorted_index)):
            continue
        root_entry = cache.load(model=stage1.MODEL_EB, epsilon=epsilon, mu=mu)
        if root_entry is None:
            warnings.append(
                "critical Pi_EB skipped because Stage-1 root cache entry is absent: "
                f"mu={mu:g}, epsilon_0={epsilon:g}, branch={row.get('branch_id')}"
            )
            pi_value = float("nan")
        else:
            roots = root_entry.roots
            root_index = int(sorted_index) - 1
            if root_index < 0 or root_index >= len(roots) or not isfinite(float(roots[root_index])):
                warnings.append(
                    "critical Pi_EB skipped because cached root index is unavailable: "
                    f"mu={mu:g}, epsilon_0={epsilon:g}, branch={row.get('branch_id')}"
                )
                pi_value = float("nan")
            else:
                result = shape_audit.eb_mode_result(
                    epsilon=epsilon,
                    beta_deg=finite_float(row, "beta_deg"),
                    eta=finite_float(row, "eta"),
                    mu=mu,
                    sorted_index=int(sorted_index),
                    Lambda=float(roots[root_index]),
                    n_points=int(n_shape_points),
                    case_role="universal_parameter_critical_postprocess",
                    root_source="stage1_refinement_root_cache",
                    root_warnings=root_entry.warnings,
                )
                pi_value = eb_pi_metrics(result)["Pi_EB"]
        rows.append(
            {
                "mu": mu,
                "branch_id": row.get("branch_id", ""),
                "mode_character": row.get("mode_character_at_crit", ""),
                "epsilon_max": finite_float(row, "epsilon_max_at_crit"),
                "chi_max": finite_float(row, "chi_max_at_crit"),
                "chi_eff": finite_float(row, "chi_eff_at_crit"),
                "Pi_EB": pi_value,
            }
        )
    return rows, warnings


def fit_summary_rows(
    mode_rows: Sequence[dict[str, object]],
    critical_rows: Sequence[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    fit_rows: list[dict[str, object]] = []
    class_rows: list[dict[str, object]] = []
    for subgroup in SUBGROUPS:
        rows = subgroup_rows(mode_rows, subgroup)
        critical_subset = (
            critical_rows
            if subgroup == "all_physical"
            else [
                row
                for row in critical_rows
                if (
                    str(row.get("mode_character", "")) == "bending_dominated"
                    if subgroup == "bending_dominated"
                    else str(row.get("mode_character", "")) != "bending_dominated"
                )
            ]
        )
        for predictor in PREDICTORS:
            x, y = finite_xy(rows, predictor)
            corr = correlation_values(x, y)
            fit = fit_power_law(x, y)
            threshold_info = threshold_classification(rows, predictor, fit["threshold_10"])
            critical_values = [finite_float(row, predictor) for row in critical_subset]
            fit_rows.append(
                {
                    "subgroup": subgroup,
                    "predictor": predictor,
                    "row_count": len(x),
                    **corr,
                    **fit,
                    "critical_count": sum(1 for value in critical_values if isfinite(value)),
                    "critical_mean": float(np.mean([value for value in critical_values if isfinite(value)]))
                    if any(isfinite(value) for value in critical_values)
                    else float("nan"),
                    "critical_std": float(np.std([value for value in critical_values if isfinite(value)], ddof=1))
                    if sum(1 for value in critical_values if isfinite(value)) >= 2
                    else float("nan"),
                    "critical_cv": coefficient_of_variation(critical_values),
                    "critical_cv_mu_means": grouped_cv(critical_subset, predictor, "mu"),
                    "critical_cv_branch_means": grouped_cv(critical_subset, predictor, "branch_id"),
                    "critical_cv_character_means": grouped_cv(critical_subset, predictor, "mode_character"),
                    "notes": "R2 is log-space; RMSE/MAE are in delta_f space",
                }
            )
            class_rows.append(
                {
                    "subgroup": subgroup,
                    "predictor": predictor,
                    **threshold_info,
                }
            )
    return fit_rows, class_rows


def cross_validation_rows(mode_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    rows_out: list[dict[str, object]] = []
    for subgroup in SUBGROUPS:
        rows = subgroup_rows(mode_rows, subgroup)
        mus = sorted({finite_float(row, "mu") for row in rows if isfinite(finite_float(row, "mu"))})
        for predictor in PREDICTORS:
            for held_mu in mus:
                train = [row for row in rows if abs(finite_float(row, "mu") - held_mu) > 1.0e-12]
                test = [row for row in rows if abs(finite_float(row, "mu") - held_mu) <= 1.0e-12]
                x_train, y_train = finite_xy(train, predictor)
                fit = fit_power_law(x_train, y_train)
                x_test, y_test = finite_xy(test, predictor)
                if (
                    len(x_test) > 0
                    and isfinite(fit["C"])
                    and isfinite(fit["p"])
                    and fit["C"] > 0.0
                ):
                    pred = fit["C"] * np.asarray(x_test, dtype=float) ** fit["p"]
                    rmse = float(np.sqrt(np.mean((pred - y_test) ** 2)))
                    mae = float(np.mean(np.abs(pred - y_test)))
                    ss_res = float(np.sum((pred - y_test) ** 2))
                    ss_tot = float(np.sum((y_test - float(np.mean(y_test))) ** 2))
                    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")
                else:
                    rmse = mae = r2 = float("nan")
                cls = threshold_classification(test, predictor, fit["threshold_10"])
                rows_out.append(
                    {
                        "subgroup": subgroup,
                        "predictor": predictor,
                        "held_out_mu": held_mu,
                        "train_count": len(x_train),
                        "test_count": len(x_test),
                        "C": fit["C"],
                        "p": fit["p"],
                        "threshold_10": fit["threshold_10"],
                        "RMSE": rmse,
                        "MAE": mae,
                        "R2": r2,
                        "false_safe_count": cls["false_safe_count"],
                        "false_unsafe_count": cls["false_unsafe_count"],
                        "false_safe_rate": cls["false_safe_rate"],
                        "false_unsafe_rate": cls["false_unsafe_rate"],
                        "balanced_accuracy": cls["balanced_accuracy"],
                        "notes": "leave-one-mu-out log-log fit; test errors in delta_f space",
                    }
                )
    return rows_out


def y_limits_for_plots(rows: Sequence[dict[str, object]]) -> tuple[float, float]:
    values = [finite_float(row, "delta_f") for row in rows if finite_float(row, "delta_f") > 0.0]
    if not values:
        return 1.0e-4, 1.0
    ymin = max(min(values) * 0.7, 1.0e-5)
    ymax = max(max(values) * 1.25, MAIN_THRESHOLD * 1.5)
    return ymin, ymax


def predictor_label(name: str) -> str:
    return {
        "epsilon_max": "epsilon_max",
        "chi_max": "chi_max",
        "chi_eff": "chi_eff",
        "Pi_EB": "Pi_EB",
    }[name]


def plot_predictor_scatter(
    rows: Sequence[dict[str, object]],
    fit_rows: Sequence[dict[str, object]],
    *,
    predictor: str,
    path: Path,
    y_limits: tuple[float, float],
    ax: plt.Axes | None = None,
) -> plt.Axes:
    own_fig = ax is None
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 5.0), constrained_layout=True)
    colors = {
        "bending_dominated": "#1f77b4",
        "mixed": "#d55e00",
        "longitudinal_dominated": "#009e73",
    }
    for character, color in colors.items():
        subset = [row for row in rows if str(row.get("Timo_classification", "")) == character]
        x = [finite_float(row, predictor) for row in subset if finite_float(row, predictor) > 0.0 and finite_float(row, "delta_f") > 0.0]
        y = [finite_float(row, "delta_f") for row in subset if finite_float(row, predictor) > 0.0 and finite_float(row, "delta_f") > 0.0]
        if x:
            ax.scatter(x, y, s=16, alpha=0.55, label=character.replace("_", " "), color=color)
    threshold_row = next(
        (
            row
            for row in fit_rows
            if row.get("subgroup") == "all_physical" and row.get("predictor") == predictor
        ),
        None,
    )
    if threshold_row is not None and isfinite(finite_float(threshold_row, "threshold_10")):
        ax.axvline(finite_float(threshold_row, "threshold_10"), color="black", lw=1.0, ls=":")
    ax.axhline(MAIN_THRESHOLD, color="black", lw=1.0, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylim(*y_limits)
    ax.set_xlabel(predictor_label(predictor))
    ax.set_ylabel("delta_f")
    ax.set_title(f"delta_f vs {predictor_label(predictor)}")
    ax.grid(True, which="both", alpha=0.25)
    if own_fig:
        ax.legend(fontsize=8)
        path.parent.mkdir(parents=True, exist_ok=True)
        ax.figure.savefig(path, dpi=180)
        plt.close(ax.figure)
    return ax


def plot_outputs(
    output_dir: Path,
    mode_rows: Sequence[dict[str, object]],
    fit_rows: Sequence[dict[str, object]],
    classification_rows: Sequence[dict[str, object]],
) -> list[Path]:
    y_limits = y_limits_for_plots(mode_rows)
    paths: list[Path] = []
    for predictor in PREDICTORS:
        path = output_dir / f"delta_f_vs_{predictor}.png"
        plot_predictor_scatter(mode_rows, fit_rows, predictor=predictor, path=path, y_limits=y_limits)
        paths.append(path)

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.0), constrained_layout=True)
    for ax, predictor in zip(axes.ravel(), PREDICTORS, strict=True):
        plot_predictor_scatter(
            mode_rows,
            fit_rows,
            predictor=predictor,
            path=output_dir / "unused.png",
            y_limits=y_limits,
            ax=ax,
        )
    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, fontsize=8)
    comparison = output_dir / "universal_parameter_collapse_comparison.png"
    fig.savefig(comparison, dpi=180)
    plt.close(fig)
    paths.append(comparison)

    class_all = [row for row in classification_rows if row.get("subgroup") == "all_physical"]
    class_bending = [row for row in classification_rows if row.get("subgroup") == "bending_dominated"]
    x = np.arange(len(PREDICTORS), dtype=float)
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.2, 5.0), constrained_layout=True)
    all_rates = [
        finite_float(next(row for row in class_all if row["predictor"] == predictor), "false_safe_rate")
        for predictor in PREDICTORS
    ]
    bending_rates = [
        finite_float(next(row for row in class_bending if row["predictor"] == predictor), "false_safe_rate")
        for predictor in PREDICTORS
    ]
    ax.bar(x - width / 2.0, all_rates, width=width, label="all physical")
    ax.bar(x + width / 2.0, bending_rates, width=width, label="bending")
    ax.set_xticks(x)
    ax.set_xticklabels([predictor_label(value) for value in PREDICTORS], rotation=20, ha="right")
    ax.set_ylabel("false-safe rate")
    ax.set_title("10% threshold false-safe comparison")
    ax.set_ylim(0.0, max([*all_rates, *bending_rates, 0.05]) * 1.25)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    false_safe = output_dir / "universal_parameter_false_safe_comparison.png"
    fig.savefig(false_safe, dpi=180)
    plt.close(fig)
    paths.append(false_safe)
    return paths


def best_row(
    rows: Sequence[dict[str, object]],
    *,
    subgroup: str,
    metric: str,
    minimize: bool,
) -> dict[str, object] | None:
    subset = [
        row
        for row in rows
        if row.get("subgroup") == subgroup and isfinite(finite_float(row, metric))
    ]
    if not subset:
        return None
    return sorted(subset, key=lambda row: finite_float(row, metric), reverse=not minimize)[0]


def row_text(row: dict[str, object] | None, metric: str) -> str:
    if row is None:
        return "not available"
    return f"{row.get('predictor')} ({metric}={finite_float(row, metric):.4g})"


def write_report(
    output_dir: Path,
    *,
    mode_rows: Sequence[dict[str, object]],
    fit_rows: Sequence[dict[str, object]],
    classification_rows: Sequence[dict[str, object]],
    cv_rows: Sequence[dict[str, object]],
    warnings: Sequence[str],
    plot_paths: Sequence[Path],
) -> Path:
    report = output_dir / "universal_parameter_comparison_report.md"
    bending_best = best_row(fit_rows, subgroup="bending_dominated", metric="R2", minimize=False)
    all_best = best_row(fit_rows, subgroup="all_physical", metric="R2", minimize=False)
    fs_best_all = best_row(classification_rows, subgroup="all_physical", metric="false_safe_rate", minimize=True)
    fs_best_bending = best_row(classification_rows, subgroup="bending_dominated", metric="false_safe_rate", minimize=True)
    epsilon_all = next(
        row
        for row in classification_rows
        if row["subgroup"] == "all_physical" and row["predictor"] == "epsilon_max"
    )
    pi_all = next(
        row
        for row in classification_rows
        if row["subgroup"] == "all_physical" and row["predictor"] == "Pi_EB"
    )
    chi_eff_bending = next(
        row
        for row in fit_rows
        if row["subgroup"] == "bending_dominated" and row["predictor"] == "chi_eff"
    )
    pi_bending = next(
        row
        for row in fit_rows
        if row["subgroup"] == "bending_dominated" and row["predictor"] == "Pi_EB"
    )
    lines = [
        "# Universal Parameter Comparison For Stage-1 EB/Timoshenko Validity",
        "",
        "## Scope",
        "",
        "This is post-processing of the completed Stage-1 EB-vs-Timoshenko study. It reuses the saved physical-branch rows and reconstructs EB mode fields from saved `Lambda_EB` values for energy integration. It does not recompute roots, change formulas, change determinants, change root solvers, change branch identities, or change the Timoshenko shear coefficient `k'`.",
        "",
        "The comparison remains a divergence between Euler--Bernoulli and Timoshenko beam theories, not an error against exact 3D elasticity.",
        "",
        "## Pi_EB",
        "",
        "`Pi_EB = Pi_shear + Pi_rotary`, where `Pi_shear = U_shear_star / U_EB` and `Pi_rotary = T_rot_star / T_trans_star`. The EB shear estimate uses `Q_EB=-EJ w'''` and the current repository `k'`. The ratios are invariant under eigenvector scaling; this is covered by a unit test.",
        "",
        "## Data",
        "",
        f"- Physical branch rows: `{len(mode_rows)}`",
        f"- Bending-dominated rows: `{len(subgroup_rows(mode_rows, 'bending_dominated'))}`",
        f"- Longitudinal/mixed rows: `{len(subgroup_rows(mode_rows, 'longitudinal_mixed'))}`",
        f"- Leave-one-mu-out rows: `{len(cv_rows)}`",
        "",
        "## Answers",
        "",
        f"1. Best collapse for bending-dominated modes by log-fit R2: `{row_text(bending_best, 'R2')}`.",
        f"2. Best collapse for the complete physical spectrum by log-fit R2: `{row_text(all_best, 'R2')}`.",
        f"3. Smallest all-spectrum false-safe rate at the inferred 10% threshold: `{row_text(fs_best_all, 'false_safe_rate')}`. For bending rows only: `{row_text(fs_best_bending, 'false_safe_rate')}`.",
        f"4. `epsilon_max` as a geometry-only rule has all-spectrum false-safe rate `{finite_float(epsilon_all, 'false_safe_rate'):.4g}` and false-unsafe rate `{finite_float(epsilon_all, 'false_unsafe_rate'):.4g}` at its fitted threshold. It is useful as a pre-solution screen, but it is not automatically the sharpest single modal criterion.",
        f"5. `Pi_EB` has all-spectrum false-safe rate `{finite_float(pi_all, 'false_safe_rate'):.4g}`. Its value is that it estimates the missing shear-flexibility and rotary-inertia effects from the EB mode itself.",
        f"6. A two-level recommendation is better than a single universal scalar at this stage: use `epsilon_max` for geometry-only screening, then use the better mode-level parameter (`Pi_EB` or `chi_eff`) for branches near the boundary. For bending rows, `chi_eff` has R2 `{finite_float(chi_eff_bending, 'R2'):.4g}` and `Pi_EB` has R2 `{finite_float(pi_bending, 'R2'):.4g}`.",
        "",
        "## Fit And Threshold Summary",
        "",
    ]
    for subgroup in SUBGROUPS:
        class_map = {
            row["predictor"]: row
            for row in classification_rows
            if row["subgroup"] == subgroup
        }
        lines.extend(
            [
                f"### {subgroup}",
                "",
                "| predictor | Pearson r | Spearman r | R2 | RMSE | MAE | threshold 10 | false-safe rate | false-unsafe rate | balanced accuracy | critical CV |",
                "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
            ]
        )
        for predictor in PREDICTORS:
            fit = next(row for row in fit_rows if row["subgroup"] == subgroup and row["predictor"] == predictor)
            cls = class_map[predictor]
            lines.append(
                "| "
                f"{predictor} | "
                f"{finite_float(fit, 'pearson_r'):.4g} | "
                f"{finite_float(fit, 'spearman_r'):.4g} | "
                f"{finite_float(fit, 'R2'):.4g} | "
                f"{finite_float(fit, 'RMSE'):.4g} | "
                f"{finite_float(fit, 'MAE'):.4g} | "
                f"{finite_float(fit, 'threshold_10'):.4g} | "
                f"{finite_float(cls, 'false_safe_rate'):.4g} | "
                f"{finite_float(cls, 'false_unsafe_rate'):.4g} | "
                f"{finite_float(cls, 'balanced_accuracy'):.4g} | "
                f"{finite_float(fit, 'critical_cv'):.4g} |"
            )
        lines.append("")
    lines.extend(
        [
            "## Cross-Validation",
            "",
            "Leave-one-mu-out results are in `universal_parameter_cross_validation.csv`. The reported RMSE/MAE are in `delta_f` space; threshold classification uses the training-set inferred 10% predictor threshold.",
            "",
            "## Plots",
            "",
        ]
    )
    for path in plot_paths:
        lines.append(f"- `{rel(path)}`")
    if warnings:
        lines.extend(["", "## Warnings", ""])
        for warning in warnings[:80]:
            lines.append(f"- {warning}")
        if len(warnings) > 80:
            lines.append(f"- ... {len(warnings) - 80} additional warnings omitted.")
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    mode_rows, warnings = build_mode_metrics(args.output_dir, n_shape_points=args.n_shape_points)
    critical_rows, critical_warnings = critical_rows_with_pi(args.output_dir, n_shape_points=args.n_shape_points)
    warnings.extend(critical_warnings)
    fit_rows, classification_rows = fit_summary_rows(mode_rows, critical_rows)
    cv_rows = cross_validation_rows(mode_rows)
    mode_fields = list(read_csv(args.output_dir / "eb_timo_mode_level_metrics.csv")[0].keys())
    mode_output_fields = mode_fields + [field for field in PI_FIELDS if field not in mode_fields]
    paths = [
        write_csv(args.output_dir / "universal_parameter_mode_metrics.csv", mode_rows, mode_output_fields),
        write_csv(args.output_dir / "universal_parameter_fit_summary.csv", fit_rows, FIT_FIELDS),
        write_csv(args.output_dir / "universal_parameter_cross_validation.csv", cv_rows, CV_FIELDS),
        write_csv(
            args.output_dir / "universal_parameter_threshold_classification.csv",
            classification_rows,
            CLASSIFICATION_FIELDS,
        ),
    ]
    plot_paths = plot_outputs(args.output_dir, mode_rows, fit_rows, classification_rows)
    report = write_report(
        args.output_dir,
        mode_rows=mode_rows,
        fit_rows=fit_rows,
        classification_rows=classification_rows,
        cv_rows=cv_rows,
        warnings=warnings,
        plot_paths=plot_paths,
    )
    print("generated universal EB applicability parameter outputs:")
    for path in [report, *paths, *plot_paths]:
        print(f"  {rel(path)}")
    print(f"physical branch rows: {len(mode_rows)}")
    print(f"warnings: {len(warnings)}")
    return {
        "mode_rows": mode_rows,
        "fit_rows": fit_rows,
        "classification_rows": classification_rows,
        "cv_rows": cv_rows,
        "paths": paths,
        "plot_paths": plot_paths,
        "report": report,
        "warnings": warnings,
    }


if __name__ == "__main__":
    main()
