from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
import sys
from typing import Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]

RESULTS_DIR = REPO_ROOT / "results"
EPSILON_TAG = "eps0p0025"

INPUT_HEATMAP = RESULTS_DIR / "diagnostic_eta_mu_beta_heatmap_metrics_eps0p0025.csv"
INPUT_LAMBDA_ETA = RESULTS_DIR / "diagnostic_lambda_eta_beta0_mu0_eps0p0025.csv"
INPUT_LAMBDA_BETA = RESULTS_DIR / "diagnostic_lambda_beta_eps0p0025_mu_eta_slices.csv"

OUTPUT_BY_MU = RESULTS_DIR / "diagnostic_eta_mu_beta_global_trends_by_mu_eps0p0025.csv"
OUTPUT_BY_ETA = RESULTS_DIR / "diagnostic_eta_mu_beta_global_trends_by_eta_eps0p0025.csv"
OUTPUT_ETA_SIGN_ABS = RESULTS_DIR / "diagnostic_eta_mu_beta_global_trends_eta_sign_abs_eps0p0025.csv"
OUTPUT_PAIR_TRANSITION = RESULTS_DIR / "diagnostic_eta_mu_beta_pair_transition_by_mu_eps0p0025.csv"
OUTPUT_BETA_HISTOGRAM = RESULTS_DIR / "diagnostic_eta_mu_beta_beta_at_gmin_histogram_eps0p0025.csv"
OUTPUT_BRANCH_DOMINANCE = RESULTS_DIR / "diagnostic_eta_mu_beta_sensitivity_branch_dominance_eps0p0025.csv"
OUTPUT_SURROGATE_FITS = RESULTS_DIR / "diagnostic_eta_mu_beta_global_surrogate_fits_eps0p0025.csv"
OUTPUT_REPORT = RESULTS_DIR / "diagnostic_eta_mu_beta_global_trends_eps0p0025_report.md"

FIG_GMIN_MU = RESULTS_DIR / "diagnostic_global_trend_gmin_vs_mu_quantiles_eps0p0025.png"
FIG_SMEAN_MU = RESULTS_DIR / "diagnostic_global_trend_Smean_vs_mu_quantiles_eps0p0025.png"
FIG_SMEANREL_MU = RESULTS_DIR / "diagnostic_global_trend_Smeanrel_vs_mu_quantiles_eps0p0025.png"
FIG_GMIN_ETA = RESULTS_DIR / "diagnostic_global_trend_gmin_vs_eta_quantiles_eps0p0025.png"
FIG_SMEAN_ETA = RESULTS_DIR / "diagnostic_global_trend_Smean_vs_eta_quantiles_eps0p0025.png"
FIG_PAIR_TRANSITION = RESULTS_DIR / "diagnostic_global_pair_transition_by_mu_eps0p0025.png"
FIG_BETA_HISTOGRAM = RESULTS_DIR / "diagnostic_global_beta_at_gmin_histogram_eps0p0025.png"
FIG_ETA_SIGN_ASYMMETRY = RESULTS_DIR / "diagnostic_global_eta_sign_asymmetry_eps0p0025.png"

PAIR_LABELS = ("1-2", "2-3", "3-4", "4-5", "5-6", "6-7")
METRIC_COLUMNS = ("g_min", "S_mean", "S_max", "S_mean_rel")
QUANTILE_METRICS = ("g_min", "S_mean", "S_max", "S_mean_rel")
BETA_BINS = (
    (0.0, 5.0),
    (5.0, 10.0),
    (10.0, 15.0),
    (15.0, 20.0),
    (20.0, 30.0),
    (30.0, 45.0),
    (45.0, 60.0),
    (60.0, 90.0),
)


def parse_float(value: str | float | int | None) -> float:
    if value is None:
        return float("nan")
    if isinstance(value, (float, int)):
        return float(value)
    text = str(value).strip()
    if not text:
        return float("nan")
    try:
        return float(text)
    except ValueError:
        return float("nan")


def number_text(value: float | int | str | None) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    value_f = float(value)
    if not np.isfinite(value_f):
        return ""
    return f"{value_f:.12g}"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_cell(value) for key, value in row.items()})


def format_cell(value: object) -> object:
    if isinstance(value, float):
        return number_text(value)
    return value


def numeric_values(rows: Sequence[Mapping[str, object]], key: str) -> np.ndarray:
    values = np.array([parse_float(row.get(key)) for row in rows], dtype=float)
    return values[np.isfinite(values)]


def grouped(rows: Sequence[dict[str, object]], key: str) -> dict[float, list[dict[str, object]]]:
    groups: dict[float, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        value = round(float(row[key]), 12)
        groups[value].append(row)
    return dict(sorted(groups.items(), key=lambda item: item[0]))


def quantile_stats(values: np.ndarray, *, include_min_max: bool) -> dict[str, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        stats = {
            "mean": float("nan"),
            "median": float("nan"),
            "q10": float("nan"),
            "q90": float("nan"),
        }
        if include_min_max:
            stats.update({"min": float("nan"), "max": float("nan")})
        return stats
    stats = {
        "mean": float(np.mean(finite)),
        "median": float(np.median(finite)),
        "q10": float(np.quantile(finite, 0.10)),
        "q90": float(np.quantile(finite, 0.90)),
    }
    if include_min_max:
        stats.update({"min": float(np.min(finite)), "max": float(np.max(finite))})
    return stats


def pair_distribution(rows: Sequence[Mapping[str, object]]) -> tuple[str, str, str]:
    counter = Counter(str(row.get("pair_at_g_min", "")).strip() for row in rows)
    total = sum(counter.get(pair, 0) for pair in PAIR_LABELS)
    if total == 0:
        return "", "", ""
    dominant = max(PAIR_LABELS, key=lambda pair: counter.get(pair, 0))
    counts = ";".join(f"{pair}:{counter.get(pair, 0)}" for pair in PAIR_LABELS)
    percents = ";".join(f"{pair}:{100.0 * counter.get(pair, 0) / total:.6g}" for pair in PAIR_LABELS)
    return dominant, counts, percents


def grouped_summary(rows: Sequence[dict[str, object]], group_key: str) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for value, local_rows in grouped(rows, group_key).items():
        row: dict[str, object] = {group_key: value, "count": len(local_rows)}
        for metric in QUANTILE_METRICS:
            stats = quantile_stats(
                numeric_values(local_rows, metric),
                include_min_max=(metric == "g_min"),
            )
            for stat_name, stat_value in stats.items():
                row[f"{metric}_{stat_name}"] = stat_value
        dominant, counts, percents = pair_distribution(local_rows)
        row["dominant_pair_at_g_min"] = dominant
        row["pair_at_g_min_distribution_counts"] = counts
        row["pair_at_g_min_distribution_percent"] = percents
        beta_values = numeric_values(local_rows, "beta_at_g_min")
        row["beta_at_g_min_mean"] = float(np.mean(beta_values)) if beta_values.size else float("nan")
        row["beta_at_g_min_median"] = float(np.median(beta_values)) if beta_values.size else float("nan")
        output.append(row)
    return output


def summary_fieldnames(group_key: str) -> list[str]:
    fields = [group_key, "count"]
    for metric in QUANTILE_METRICS:
        stats = ["mean", "median", "q10", "q90"]
        if metric == "g_min":
            stats = ["mean", "median", "min", "max", "q10", "q90"]
        fields.extend(f"{metric}_{stat}" for stat in stats)
    fields.extend(
        [
            "dominant_pair_at_g_min",
            "pair_at_g_min_distribution_counts",
            "pair_at_g_min_distribution_percent",
            "beta_at_g_min_mean",
            "beta_at_g_min_median",
        ]
    )
    return fields


def eta_sign(value: float) -> str:
    if value > 1.0e-12:
        return "positive"
    if value < -1.0e-12:
        return "negative"
    return "zero"


def eta_sign_abs_rows(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    sign_groups: dict[str, list[dict[str, object]]] = {"positive": [], "negative": [], "zero": []}
    sign_abs_groups: dict[tuple[str, float], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        sign = eta_sign(float(row["eta"]))
        abs_eta = round(abs(float(row["eta"])), 12)
        sign_groups[sign].append(row)
        sign_abs_groups[(sign, abs_eta)].append(row)

    for sign, local_rows in sign_groups.items():
        output.append(metric_summary_row("eta_sign_overall", sign, float("nan"), local_rows))
    for (sign, abs_eta), local_rows in sorted(sign_abs_groups.items(), key=lambda item: (item[0][1], item[0][0])):
        output.append(metric_summary_row("eta_sign_abs", sign, abs_eta, local_rows))

    paired_rows = paired_eta_difference_rows(rows)
    output.extend(paired_rows)
    return output


def metric_summary_row(row_type: str, sign: str, abs_eta: float, rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    row: dict[str, object] = {
        "row_type": row_type,
        "eta_sign": sign,
        "abs_eta": abs_eta,
        "count": len(rows),
    }
    for metric in METRIC_COLUMNS:
        values = numeric_values(rows, metric)
        row[f"{metric}_mean"] = float(np.mean(values)) if values.size else float("nan")
        row[f"{metric}_median"] = float(np.median(values)) if values.size else float("nan")
    return row


def paired_eta_difference_rows(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    lookup = {
        (round(float(row["mu"]), 12), round(float(row["eta"]), 12)): row
        for row in rows
    }
    abs_values = sorted({round(abs(float(row["eta"])), 12) for row in rows if abs(float(row["eta"])) > 1.0e-12})
    mu_values = sorted({round(float(row["mu"]), 12) for row in rows})
    output: list[dict[str, object]] = []
    global_diffs: dict[str, list[float]] = {
        "S_mean": [],
        "S_mean_rel": [],
        "g_min": [],
    }
    for abs_eta in abs_values:
        diffs: dict[str, list[float]] = {"S_mean": [], "S_mean_rel": [], "g_min": []}
        for mu in mu_values:
            plus = lookup.get((mu, abs_eta))
            minus = lookup.get((mu, -abs_eta))
            if plus is None or minus is None:
                continue
            for metric in diffs:
                value = parse_float(plus.get(metric)) - parse_float(minus.get(metric))
                if np.isfinite(value):
                    diffs[metric].append(float(value))
                    global_diffs[metric].append(float(value))
        output.append(paired_difference_summary("paired_abs_eta", abs_eta, diffs))
    output.append(paired_difference_summary("paired_global", float("nan"), global_diffs))
    return output


def paired_difference_summary(row_type: str, abs_eta: float, diffs: Mapping[str, Sequence[float]]) -> dict[str, object]:
    row: dict[str, object] = {
        "row_type": row_type,
        "eta_sign": "positive_minus_negative",
        "abs_eta": abs_eta,
        "count": len(next(iter(diffs.values()))) if diffs else 0,
    }
    for metric, values in diffs.items():
        array = np.asarray(values, dtype=float)
        finite = array[np.isfinite(array)]
        row[f"{metric}_paired_mean_diff"] = float(np.mean(finite)) if finite.size else float("nan")
        row[f"{metric}_paired_median_diff"] = float(np.median(finite)) if finite.size else float("nan")
        row[f"{metric}_paired_mean_abs_diff"] = float(np.mean(np.abs(finite))) if finite.size else float("nan")
        row[f"{metric}_paired_max_abs_diff"] = float(np.max(np.abs(finite))) if finite.size else float("nan")
    return row


def eta_sign_abs_fieldnames() -> list[str]:
    fields = ["row_type", "eta_sign", "abs_eta", "count"]
    for metric in METRIC_COLUMNS:
        fields.extend([f"{metric}_mean", f"{metric}_median"])
    for metric in ("S_mean", "S_mean_rel", "g_min"):
        fields.extend(
            [
                f"{metric}_paired_mean_diff",
                f"{metric}_paired_median_diff",
                f"{metric}_paired_mean_abs_diff",
                f"{metric}_paired_max_abs_diff",
            ]
        )
    return fields


def pair_transition_rows(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for mu, local_rows in grouped(rows, "mu").items():
        counter = Counter(str(row.get("pair_at_g_min", "")).strip() for row in local_rows)
        total = sum(counter.get(pair, 0) for pair in PAIR_LABELS)
        percentages = {
            pair: (100.0 * counter.get(pair, 0) / total if total else 0.0)
            for pair in PAIR_LABELS
        }
        probabilities = np.array([percentages[pair] / 100.0 for pair in PAIR_LABELS], dtype=float)
        nonzero = probabilities[probabilities > 0.0]
        entropy = float(-np.sum(nonzero * np.log(nonzero))) if nonzero.size else 0.0
        normalized_entropy = float(entropy / np.log(len(PAIR_LABELS))) if len(PAIR_LABELS) > 1 else 0.0
        dominant = max(PAIR_LABELS, key=lambda pair: percentages[pair]) if total else ""
        concentration = float(max(percentages.values())) if percentages else 0.0
        row: dict[str, object] = {
            "mu": mu,
            "count": total,
            "dominant_pair": dominant,
            "entropy": entropy,
            "normalized_entropy": normalized_entropy,
            "concentration_percent": concentration,
        }
        for pair in PAIR_LABELS:
            row[f"pair_{pair.replace('-', '_')}_percent"] = percentages[pair]
            row[f"pair_{pair.replace('-', '_')}_count"] = counter.get(pair, 0)
        output.append(row)
    return output


def pair_transition_fieldnames() -> list[str]:
    fields = ["mu", "count", "dominant_pair", "entropy", "normalized_entropy", "concentration_percent"]
    for pair in PAIR_LABELS:
        token = pair.replace("-", "_")
        fields.extend([f"pair_{token}_percent", f"pair_{token}_count"])
    return fields


def beta_histogram_rows(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    values = numeric_values(rows, "beta_at_g_min")
    total = len(values)
    output: list[dict[str, object]] = []
    for idx, (left, right) in enumerate(BETA_BINS):
        if idx == len(BETA_BINS) - 1:
            mask = (values >= left) & (values <= right)
        else:
            mask = (values >= left) & (values < right)
        count = int(np.count_nonzero(mask))
        output.append(
            {
                "bin_left_deg": left,
                "bin_right_deg": right,
                "bin_label": f"{left:g}-{right:g}",
                "count": count,
                "percent": 100.0 * count / total if total else float("nan"),
            }
        )
    return output


def branch_dominance_rows(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    absolute_counter: Counter[int] = Counter()
    relative_counter: Counter[int] = Counter()
    for row in rows:
        s_values = np.array([parse_float(row.get(f"S_branch{idx}")) for idx in range(1, 7)], dtype=float)
        rel_values = np.array([parse_float(row.get(f"S_rel_branch{idx}")) for idx in range(1, 7)], dtype=float)
        abs_branch = int(np.nanargmax(s_values)) + 1 if np.any(np.isfinite(s_values)) else -1
        rel_branch = int(np.nanargmax(rel_values)) + 1 if np.any(np.isfinite(rel_values)) else -1
        if abs_branch > 0:
            absolute_counter[abs_branch] += 1
        if rel_branch > 0:
            relative_counter[rel_branch] += 1
        output.append(
            {
                "row_type": "point",
                "eta": row["eta"],
                "mu": row["mu"],
                "absolute_sensitivity_dominant_branch": abs_branch if abs_branch > 0 else "",
                "relative_sensitivity_dominant_branch": rel_branch if rel_branch > 0 else "",
                "count": "",
                "percent": "",
            }
        )

    total = len(rows)
    for branch in range(1, 7):
        output.append(
            {
                "row_type": "absolute_count",
                "branch": branch,
                "count": absolute_counter.get(branch, 0),
                "percent": 100.0 * absolute_counter.get(branch, 0) / total if total else float("nan"),
            }
        )
    for branch in range(1, 7):
        output.append(
            {
                "row_type": "relative_count",
                "branch": branch,
                "count": relative_counter.get(branch, 0),
                "percent": 100.0 * relative_counter.get(branch, 0) / total if total else float("nan"),
            }
        )
    return output


def branch_dominance_fieldnames() -> list[str]:
    return [
        "row_type",
        "eta",
        "mu",
        "absolute_sensitivity_dominant_branch",
        "relative_sensitivity_dominant_branch",
        "branch",
        "count",
        "percent",
    ]


def design_matrix(mu: np.ndarray, eta: np.ndarray) -> np.ndarray:
    return np.column_stack(
        [
            np.ones_like(mu),
            mu,
            mu**2,
            eta,
            eta**2,
            mu * eta,
        ]
    )


def fit_surrogate(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    mu = numeric_values(rows, "mu")
    eta = numeric_values(rows, "eta")
    if len(mu) != len(rows) or len(eta) != len(rows):
        raise RuntimeError("mu/eta columns must be finite for surrogate fits")
    output: list[dict[str, object]] = []
    for metric in ("g_min", "S_mean", "S_max", "S_mean_rel", "beta_at_g_min"):
        y = np.array([parse_float(row.get(metric)) for row in rows], dtype=float)
        output.append(fit_metric(metric, "identity", y, mu, eta))
        if metric == "g_min" and np.all(np.isfinite(y)) and np.all(y > 0.0):
            output.append(fit_metric(metric, "log", np.log(y), mu, eta))
    return output


def fit_metric(metric: str, transform: str, y: np.ndarray, mu: np.ndarray, eta: np.ndarray) -> dict[str, object]:
    mask = np.isfinite(y) & np.isfinite(mu) & np.isfinite(eta)
    y_fit = y[mask]
    x_fit = design_matrix(mu[mask], eta[mask])
    coeffs, *_rest = np.linalg.lstsq(x_fit, y_fit, rcond=None)
    prediction = x_fit @ coeffs
    residual = y_fit - prediction
    ss_res = float(np.sum(residual**2))
    ss_tot = float(np.sum((y_fit - np.mean(y_fit)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")
    rmse = float(np.sqrt(np.mean(residual**2))) if y_fit.size else float("nan")
    return {
        "metric": metric,
        "transform": transform,
        "n": int(y_fit.size),
        "intercept": float(coeffs[0]),
        "mu": float(coeffs[1]),
        "mu2": float(coeffs[2]),
        "eta": float(coeffs[3]),
        "eta2": float(coeffs[4]),
        "mu_eta": float(coeffs[5]),
        "r2": r2,
        "rmse": rmse,
    }


def one_variable_r2(rows: Sequence[dict[str, object]], metric: str, variable: str) -> float:
    x = np.array([parse_float(row.get(variable)) for row in rows], dtype=float)
    y = np.array([parse_float(row.get(metric)) for row in rows], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if x.size < 3:
        return float("nan")
    X = np.column_stack([np.ones_like(x), x, x**2])
    coeffs, *_rest = np.linalg.lstsq(X, y, rcond=None)
    prediction = X @ coeffs
    ss_res = float(np.sum((y - prediction) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")


def median_slope(summary_rows: Sequence[Mapping[str, object]], x_key: str, metric_key: str) -> float:
    x = np.array([parse_float(row.get(x_key)) for row in summary_rows], dtype=float)
    y = np.array([parse_float(row.get(metric_key)) for row in summary_rows], dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(mask) < 2:
        return float("nan")
    slope, _intercept = np.polyfit(x[mask], y[mask], 1)
    return float(slope)


def trend_word(slope: float, *, positive_word: str = "increases", negative_word: str = "decreases") -> str:
    if not np.isfinite(slope):
        return "is unclear for"
    if abs(slope) < 1.0e-12:
        return "is nearly flat for"
    return positive_word if slope > 0.0 else negative_word


def plot_quantile_trend(
    rows: Sequence[Mapping[str, object]],
    x_key: str,
    metric: str,
    output: Path,
    *,
    xlabel: str,
    ylabel: str,
    title: str,
) -> None:
    sorted_rows = sorted(rows, key=lambda row: parse_float(row.get(x_key)))
    x = np.array([parse_float(row.get(x_key)) for row in sorted_rows], dtype=float)
    mean = np.array([parse_float(row.get(f"{metric}_mean")) for row in sorted_rows], dtype=float)
    median = np.array([parse_float(row.get(f"{metric}_median")) for row in sorted_rows], dtype=float)
    q10 = np.array([parse_float(row.get(f"{metric}_q10")) for row in sorted_rows], dtype=float)
    q90 = np.array([parse_float(row.get(f"{metric}_q90")) for row in sorted_rows], dtype=float)

    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    ax.fill_between(x, q10, q90, color="#9ecae1", alpha=0.45, label="q10..q90")
    ax.plot(x, median, color="#08519c", lw=1.8, label="median")
    ax.plot(x, mean, color="#cb181d", lw=1.3, ls="--", label="mean")
    finite = median[np.isfinite(median) & (median > 0.0)]
    if metric == "g_min" and finite.size and float(np.max(finite) / np.min(finite)) > 20.0:
        ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, color="0.88", linewidth=0.6, which="both")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=230, bbox_inches="tight")
    plt.close(fig)


def plot_pair_transition(rows: Sequence[Mapping[str, object]], output: Path) -> None:
    sorted_rows = sorted(rows, key=lambda row: parse_float(row.get("mu")))
    mu = np.array([parse_float(row.get("mu")) for row in sorted_rows], dtype=float)
    values = np.vstack(
        [
            np.array([parse_float(row.get(f"pair_{pair.replace('-', '_')}_percent")) for row in sorted_rows], dtype=float)
            for pair in PAIR_LABELS
        ]
    )
    fig, ax = plt.subplots(figsize=(9.2, 5.3))
    ax.stackplot(mu, values, labels=PAIR_LABELS, alpha=0.86)
    ax.set_xlabel("mu")
    ax.set_ylabel("percent")
    ax.set_ylim(0.0, 100.0)
    ax.set_title("Dominant adjacent-gap pair distribution by mu")
    ax.grid(True, color="0.9", linewidth=0.5)
    ax.legend(title="pair", loc="upper left", ncol=3, fontsize=8, frameon=False)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=230, bbox_inches="tight")
    plt.close(fig)


def plot_beta_histogram(rows: Sequence[Mapping[str, object]], output: Path) -> None:
    labels = [str(row["bin_label"]) for row in rows]
    percents = np.array([parse_float(row.get("percent")) for row in rows], dtype=float)
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    ax.bar(labels, percents, color="#31a354", edgecolor="white")
    ax.set_xlabel("beta_at_g_min bin (deg)")
    ax.set_ylabel("percent")
    ax.set_title("Distribution of beta at minimum adjacent gap")
    ax.grid(True, axis="y", color="0.88", linewidth=0.6)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=230, bbox_inches="tight")
    plt.close(fig)


def plot_eta_sign_asymmetry(rows: Sequence[Mapping[str, object]], output: Path) -> None:
    paired = [
        row for row in rows
        if row.get("row_type") == "paired_abs_eta" and parse_float(row.get("abs_eta")) > 0.0
    ]
    paired = sorted(paired, key=lambda row: parse_float(row.get("abs_eta")))
    abs_eta = np.array([parse_float(row.get("abs_eta")) for row in paired], dtype=float)
    metrics = ("g_min", "S_mean", "S_mean_rel")
    fig, axes = plt.subplots(3, 1, figsize=(8.8, 8.0), sharex=True, constrained_layout=True)
    for ax, metric in zip(axes, metrics, strict=True):
        mean_diff = np.array([parse_float(row.get(f"{metric}_paired_mean_diff")) for row in paired], dtype=float)
        mean_abs = np.array([parse_float(row.get(f"{metric}_paired_mean_abs_diff")) for row in paired], dtype=float)
        ax.axhline(0.0, color="0.45", lw=0.8)
        ax.plot(abs_eta, mean_diff, marker="o", lw=1.4, color="#08519c", label="+eta - -eta")
        ax.plot(abs_eta, mean_abs, marker="s", lw=1.2, color="#cb181d", label="mean abs diff")
        ax.set_ylabel(metric)
        ax.grid(True, color="0.88", linewidth=0.6)
    axes[-1].set_xlabel("|eta|")
    axes[0].set_title("Fixed-mu eta sign asymmetry")
    axes[0].legend(frameon=False, fontsize=8)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=230, bbox_inches="tight")
    plt.close(fig)


def dominant_summary(rows: Sequence[Mapping[str, object]], row_type: str) -> tuple[int, float]:
    candidates = [row for row in rows if row.get("row_type") == row_type]
    if not candidates:
        return -1, float("nan")
    best = max(candidates, key=lambda row: parse_float(row.get("count")))
    return int(parse_float(best.get("branch"))), parse_float(best.get("percent"))


def best_surrogate_r2(rows: Sequence[Mapping[str, object]], metric: str, transform: str = "identity") -> float:
    for row in rows:
        if row.get("metric") == metric and row.get("transform") == transform:
            return parse_float(row.get("r2"))
    return float("nan")


def report_lines(
    *,
    heatmap_rows: Sequence[dict[str, object]],
    lambda_eta_rows: Sequence[dict[str, str]],
    lambda_beta_rows: Sequence[dict[str, str]],
    by_mu_rows: Sequence[dict[str, object]],
    by_eta_rows: Sequence[dict[str, object]],
    eta_sign_rows: Sequence[dict[str, object]],
    pair_rows: Sequence[dict[str, object]],
    beta_hist_rows: Sequence[dict[str, object]],
    branch_rows: Sequence[dict[str, object]],
    surrogate_rows: Sequence[dict[str, object]],
) -> list[str]:
    driver = {
        metric: {
            "mu_r2": one_variable_r2(heatmap_rows, metric, "mu"),
            "eta_r2": one_variable_r2(heatmap_rows, metric, "eta"),
        }
        for metric in ("g_min", "S_mean", "S_mean_rel")
    }
    gmin_mu_slope = median_slope(by_mu_rows, "mu", "g_min_median")
    smean_mu_slope = median_slope(by_mu_rows, "mu", "S_mean_median")
    smeanrel_mu_slope = median_slope(by_mu_rows, "mu", "S_mean_rel_median")
    dominant_pair_sequence = [(row["mu"], row["dominant_pair"]) for row in pair_rows]
    pair_changes = sum(
        1
        for left, right in zip(dominant_pair_sequence, dominant_pair_sequence[1:])
        if left[1] != right[1]
    )
    beta_best_bin = max(beta_hist_rows, key=lambda row: parse_float(row.get("count")))
    abs_branch, abs_percent = dominant_summary(branch_rows, "absolute_count")
    rel_branch, rel_percent = dominant_summary(branch_rows, "relative_count")
    paired_global = next(row for row in eta_sign_rows if row.get("row_type") == "paired_global")

    smooth_notes = []
    for metric in ("g_min", "S_mean", "S_max", "S_mean_rel", "beta_at_g_min"):
        r2 = best_surrogate_r2(surrogate_rows, metric)
        smooth_notes.append(f"{metric}: R2={r2:.3f}" if np.isfinite(r2) else f"{metric}: R2 unavailable")
    log_r2 = best_surrogate_r2(surrogate_rows, "g_min", transform="log")

    stronger_g = "mu" if driver["g_min"]["mu_r2"] > driver["g_min"]["eta_r2"] else "eta"
    stronger_s = "mu" if driver["S_mean"]["mu_r2"] > driver["S_mean"]["eta_r2"] else "eta"
    stronger_sr = "mu" if driver["S_mean_rel"]["mu_r2"] > driver["S_mean_rel"]["eta_r2"] else "eta"

    return [
        "# Diagnostic Eta-Mu-Beta Global Trends",
        "",
        "## Scope",
        "",
        "Global diagnostic post-processing of already generated sorted-frequency",
        "CSV files. This script does not recompute roots, does not do local",
        "candidate refinement, does not run strict positive-gap verification, and",
        "does not make article figures.",
        "",
        "## Input Files",
        "",
        f"- heatmap metrics rows: `{len(heatmap_rows)}` from `{INPUT_HEATMAP.relative_to(REPO_ROOT)}`",
        f"- Lambda(eta) rows: `{len(lambda_eta_rows)}` from `{INPUT_LAMBDA_ETA.relative_to(REPO_ROOT)}`",
        f"- Lambda(beta) slice rows: `{len(lambda_beta_rows)}` from `{INPUT_LAMBDA_BETA.relative_to(REPO_ROOT)}`",
        "",
        "## Main Answers",
        "",
        f"- Stronger global driver of `g_min`: `{stronger_g}` by one-variable quadratic R2 "
        f"(mu `{driver['g_min']['mu_r2']:.3f}`, eta `{driver['g_min']['eta_r2']:.3f}`).",
        f"- Stronger global driver of `S_mean`: `{stronger_s}` by one-variable quadratic R2 "
        f"(mu `{driver['S_mean']['mu_r2']:.3f}`, eta `{driver['S_mean']['eta_r2']:.3f}`).",
        f"- Stronger global driver of `S_mean_rel`: `{stronger_sr}` by one-variable quadratic R2 "
        f"(mu `{driver['S_mean_rel']['mu_r2']:.3f}`, eta `{driver['S_mean_rel']['eta_r2']:.3f}`).",
        f"- Increasing mu generally {trend_word(gmin_mu_slope)} median `g_min` over the global eta grid "
        f"(linear slope of mu-group medians `{gmin_mu_slope:.6g}`).",
        f"- Increasing mu generally {trend_word(smean_mu_slope)} median `S_mean` "
        f"(slope `{smean_mu_slope:.6g}`) and {trend_word(smeanrel_mu_slope)} median `S_mean_rel` "
        f"(slope `{smeanrel_mu_slope:.6g}`).",
        f"- Positive and negative eta are not symmetric at fixed positive mu on this diagnostic map: global mean absolute paired differences are "
        f"`g_min={parse_float(paired_global.get('g_min_paired_mean_abs_diff')):.6g}`, "
        f"`S_mean={parse_float(paired_global.get('S_mean_paired_mean_abs_diff')):.6g}`, and "
        f"`S_mean_rel={parse_float(paired_global.get('S_mean_rel_paired_mean_abs_diff')):.6g}`.",
        f"- The pair responsible for `g_min` migrates with mu: the dominant pair changes `{pair_changes}` times across the mu grid.",
        f"- `beta_at_g_min` values concentrate most in `{beta_best_bin['bin_label']}` deg "
        f"({parse_float(beta_best_bin.get('percent')):.3g}% of eta-mu points).",
        f"- Absolute sensitivity is most often dominated by branch `{abs_branch}` "
        f"({abs_percent:.3g}% of eta-mu points).",
        f"- Relative sensitivity is most often dominated by branch `{rel_branch}` "
        f"({rel_percent:.3g}% of eta-mu points).",
        "- Low-order surrogate smoothness check: " + "; ".join(smooth_notes)
        + (f"; log(g_min): R2={log_r2:.3f}" if np.isfinite(log_r2) else ""),
        "",
        "These trends are descriptive global summaries only. Small-gap regions",
        "remain candidate close approaches and require strict local verification",
        "before any crossing/no-crossing or article claim.",
        "",
        "## Output Files",
        "",
        f"- `{OUTPUT_BY_MU.relative_to(REPO_ROOT)}`",
        f"- `{OUTPUT_BY_ETA.relative_to(REPO_ROOT)}`",
        f"- `{OUTPUT_ETA_SIGN_ABS.relative_to(REPO_ROOT)}`",
        f"- `{OUTPUT_PAIR_TRANSITION.relative_to(REPO_ROOT)}`",
        f"- `{OUTPUT_BETA_HISTOGRAM.relative_to(REPO_ROOT)}`",
        f"- `{OUTPUT_BRANCH_DOMINANCE.relative_to(REPO_ROOT)}`",
        f"- `{OUTPUT_SURROGATE_FITS.relative_to(REPO_ROOT)}`",
        f"- `{FIG_GMIN_MU.relative_to(REPO_ROOT)}`",
        f"- `{FIG_SMEAN_MU.relative_to(REPO_ROOT)}`",
        f"- `{FIG_SMEANREL_MU.relative_to(REPO_ROOT)}`",
        f"- `{FIG_GMIN_ETA.relative_to(REPO_ROOT)}`",
        f"- `{FIG_SMEAN_ETA.relative_to(REPO_ROOT)}`",
        f"- `{FIG_PAIR_TRANSITION.relative_to(REPO_ROOT)}`",
        f"- `{FIG_BETA_HISTOGRAM.relative_to(REPO_ROOT)}`",
        f"- `{FIG_ETA_SIGN_ASYMMETRY.relative_to(REPO_ROOT)}`",
        "",
    ]


def normalize_heatmap_rows(raw_rows: Sequence[dict[str, str]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for raw in raw_rows:
        row: dict[str, object] = dict(raw)
        for key in (
            "eta",
            "mu",
            "epsilon",
            "g_min",
            "beta_at_g_min",
            "S_mean",
            "S_max",
            "S_mean_rel",
            *[f"S_branch{idx}" for idx in range(1, 7)],
            *[f"S_rel_branch{idx}" for idx in range(1, 7)],
        ):
            row[key] = parse_float(raw.get(key))
        rows.append(row)
    return rows


def main() -> dict[str, object]:
    heatmap_raw = read_csv(INPUT_HEATMAP)
    lambda_eta_rows = read_csv(INPUT_LAMBDA_ETA)
    lambda_beta_rows = read_csv(INPUT_LAMBDA_BETA)
    heatmap_rows = normalize_heatmap_rows(heatmap_raw)

    by_mu_rows = grouped_summary(heatmap_rows, "mu")
    by_eta_rows = grouped_summary(heatmap_rows, "eta")
    eta_sign_rows = eta_sign_abs_rows(heatmap_rows)
    pair_rows = pair_transition_rows(heatmap_rows)
    beta_hist_rows = beta_histogram_rows(heatmap_rows)
    branch_rows = branch_dominance_rows(heatmap_rows)
    surrogate_rows = fit_surrogate(heatmap_rows)

    write_csv(OUTPUT_BY_MU, summary_fieldnames("mu"), by_mu_rows)
    write_csv(OUTPUT_BY_ETA, summary_fieldnames("eta"), by_eta_rows)
    write_csv(OUTPUT_ETA_SIGN_ABS, eta_sign_abs_fieldnames(), eta_sign_rows)
    write_csv(OUTPUT_PAIR_TRANSITION, pair_transition_fieldnames(), pair_rows)
    write_csv(OUTPUT_BETA_HISTOGRAM, ["bin_left_deg", "bin_right_deg", "bin_label", "count", "percent"], beta_hist_rows)
    write_csv(OUTPUT_BRANCH_DOMINANCE, branch_dominance_fieldnames(), branch_rows)
    write_csv(
        OUTPUT_SURROGATE_FITS,
        ["metric", "transform", "n", "intercept", "mu", "mu2", "eta", "eta2", "mu_eta", "r2", "rmse"],
        surrogate_rows,
    )

    plot_quantile_trend(by_mu_rows, "mu", "g_min", FIG_GMIN_MU, xlabel="mu", ylabel="g_min", title="g_min quantiles by mu")
    plot_quantile_trend(by_mu_rows, "mu", "S_mean", FIG_SMEAN_MU, xlabel="mu", ylabel="S_mean", title="S_mean quantiles by mu")
    plot_quantile_trend(by_mu_rows, "mu", "S_mean_rel", FIG_SMEANREL_MU, xlabel="mu", ylabel="S_mean_rel", title="S_mean_rel quantiles by mu")
    plot_quantile_trend(by_eta_rows, "eta", "g_min", FIG_GMIN_ETA, xlabel="eta", ylabel="g_min", title="g_min quantiles by eta")
    plot_quantile_trend(by_eta_rows, "eta", "S_mean", FIG_SMEAN_ETA, xlabel="eta", ylabel="S_mean", title="S_mean quantiles by eta")
    plot_pair_transition(pair_rows, FIG_PAIR_TRANSITION)
    plot_beta_histogram(beta_hist_rows, FIG_BETA_HISTOGRAM)
    plot_eta_sign_asymmetry(eta_sign_rows, FIG_ETA_SIGN_ASYMMETRY)

    OUTPUT_REPORT.write_text(
        "\n".join(
            report_lines(
                heatmap_rows=heatmap_rows,
                lambda_eta_rows=lambda_eta_rows,
                lambda_beta_rows=lambda_beta_rows,
                by_mu_rows=by_mu_rows,
                by_eta_rows=by_eta_rows,
                eta_sign_rows=eta_sign_rows,
                pair_rows=pair_rows,
                beta_hist_rows=beta_hist_rows,
                branch_rows=branch_rows,
                surrogate_rows=surrogate_rows,
            )
        ),
        encoding="utf-8",
    )

    print("diagnostic eta-mu-beta global trend post-processing")
    print(f"read heatmap rows: {len(heatmap_rows)}")
    print(f"read Lambda(eta) rows: {len(lambda_eta_rows)}")
    print(f"read Lambda(beta) rows: {len(lambda_beta_rows)}")
    print(f"saved by-mu summary: {OUTPUT_BY_MU}")
    print(f"saved by-eta summary: {OUTPUT_BY_ETA}")
    print(f"saved report: {OUTPUT_REPORT}")
    return {
        "heatmap_rows": len(heatmap_rows),
        "by_mu": OUTPUT_BY_MU,
        "by_eta": OUTPUT_BY_ETA,
        "report": OUTPUT_REPORT,
    }


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
