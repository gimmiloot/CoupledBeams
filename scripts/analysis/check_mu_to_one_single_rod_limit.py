from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
from pathlib import Path
import sys
from typing import Sequence

import matplotlib
import mpmath as mp

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import BeamParams  # noqa: E402
from my_project.analytic.solvers import find_first_n_roots  # noqa: E402
from my_project.fem import python_fem as fem  # noqa: E402
from scripts.lib.analytic_branch_tracking import track_mu_sweep  # noqa: E402


DEFAULT_EPSILON = 0.0025
DEFAULT_BETAS = [5.0, 10.0, 15.0, 30.0]
DEFAULT_MUS = [0.8, 0.9, 0.95, 0.98, 0.99]
DEFAULT_NUM_ROOTS = 12
DEFAULT_L_BASE = 1.0

REFERENCE_ROOT_COUNT = 12
KNOWN_FP_PREFIX = np.array([3.9266, 7.0686, 10.2102], dtype=float)
KNOWN_FF_PREFIX = np.array([4.730040744862704, 7.853204624095838, 10.99560783800167], dtype=float)
BRANCH_FOCUS_BETA = 15.0
BRANCH_FOCUS_MUS = [0.9, 0.95]
BRANCH_FOCUS_IDS = ["bending_desc_01", "bending_desc_02", "bending_desc_04"]

RESULTS_DIR = REPO_ROOT / "results"
ROOTS_CSV = RESULTS_DIR / "mu_to_one_single_rod_limit_roots.csv"
FEM_ROOTS_CSV = RESULTS_DIR / "mu_to_one_single_rod_limit_fem_roots.csv"
BRANCHES_CSV = RESULTS_DIR / "mu_to_one_single_rod_limit_branches.csv"
REPORT_MD = RESULTS_DIR / "mu_to_one_single_rod_limit_report.md"

ROOT_FIELDNAMES = [
    "beta",
    "mu",
    "epsilon",
    "root_index",
    "Lambda_system",
    "nearest_FP_finite_mode",
    "nearest_FP_finite_Lambda",
    "abs_diff_FP_finite",
    "rel_diff_FP_finite",
    "nearest_FF_finite_mode",
    "nearest_FF_finite_Lambda",
    "abs_diff_FF_finite",
    "rel_diff_FF_finite",
    "nearest_FP_limit_mode",
    "nearest_FP_limit_Lambda",
    "abs_diff_FP_limit",
    "rel_diff_FP_limit",
    "nearest_FF_limit_mode",
    "nearest_FF_limit_Lambda",
    "abs_diff_FF_limit",
    "rel_diff_FF_limit",
    "best_match_type",
    "best_match_mode",
    "best_match_abs_diff",
    "best_match_rel_diff",
]

FEM_FIELDNAMES = [
    "beta",
    "mu",
    "epsilon",
    "root_index",
    "Lambda_fem",
    "Lambda_analytic_same_index",
    "rel_diff_to_analytic_same_index",
    "nearest_FP_finite_mode",
    "nearest_FP_finite_Lambda",
    "abs_diff_FP_finite",
    "rel_diff_FP_finite",
    "nearest_FF_finite_mode",
    "nearest_FF_finite_Lambda",
    "abs_diff_FF_finite",
    "rel_diff_FF_finite",
    "nearest_FP_limit_mode",
    "nearest_FP_limit_Lambda",
    "abs_diff_FP_limit",
    "rel_diff_FP_limit",
    "nearest_FF_limit_mode",
    "nearest_FF_limit_Lambda",
    "abs_diff_FF_limit",
    "rel_diff_FF_limit",
    "best_match_type",
    "best_match_mode",
    "best_match_abs_diff",
    "best_match_rel_diff",
]

BRANCH_FIELDNAMES = [
    "beta",
    "mu",
    "branch_id",
    "current_sorted_index",
    "Lambda_branch",
    "best_match_type",
    "best_match_mode",
    "best_match_rel_diff",
    "rel_diff_FP_finite",
    "rel_diff_FF_finite",
    "rel_diff_FP_limit",
    "rel_diff_FF_limit",
]

FEM_STATE_KEYS = (
    "r",
    "E",
    "rho",
    "L_tot",
    "ell",
    "A",
    "I",
    "EI",
    "EA",
    "rhoA",
    "eps",
    "scale",
    "EI_nd",
    "rhoA_nd",
    "EA_nd",
)


def parse_float_list(values: Sequence[float | str]) -> list[float]:
    return [float(value) for value in values]


def finite_rel_diff(value: float, reference: float) -> float:
    if not np.isfinite(value) or not np.isfinite(reference) or reference == 0.0:
        return float("nan")
    return abs(value - reference) / abs(reference)


def bisect_root(func, left: float, right: float, *, iterations: int = 80) -> float:
    f_left = float(func(left))
    f_right = float(func(right))
    if not (np.isfinite(f_left) and np.isfinite(f_right)):
        raise RuntimeError(f"Non-finite bracket values for [{left:g}, {right:g}].")
    if f_left == 0.0:
        return float(left)
    if f_right == 0.0:
        return float(right)
    if f_left * f_right > 0.0:
        raise RuntimeError(
            f"Root bracket does not change sign: [{left:g}, {right:g}], "
            f"f(left)={f_left:g}, f(right)={f_right:g}."
        )
    a = float(left)
    b = float(right)
    fa = f_left
    for _ in range(iterations):
        mid = 0.5 * (a + b)
        fm = float(func(mid))
        if fm == 0.0:
            return mid
        if fa * fm <= 0.0:
            b = mid
        else:
            a = mid
            fa = fm
    return 0.5 * (a + b)


def clamped_pinned_characteristic(alpha: float) -> float:
    # Clamped-pinned Euler-Bernoulli beam: tan(alpha) = tanh(alpha),
    # equivalently sin(alpha) cosh(alpha) - cos(alpha) sinh(alpha) = 0.
    return float(np.sin(alpha) * np.cosh(alpha) - np.cos(alpha) * np.sinh(alpha))


def clamped_clamped_characteristic(alpha: float) -> float:
    # Clamped-clamped Euler-Bernoulli beam: cosh(alpha) cos(alpha) = 1.
    return float(np.cosh(alpha) * np.cos(alpha) - 1.0)


def clamped_pinned_roots(n_roots: int) -> np.ndarray:
    pad = 1e-9
    roots = []
    for mode in range(1, int(n_roots) + 1):
        left = mode * np.pi + pad
        right = (mode + 0.5) * np.pi - pad
        roots.append(bisect_root(clamped_pinned_characteristic, left, right))
    return np.asarray(roots, dtype=float)


def clamped_clamped_roots(n_roots: int) -> np.ndarray:
    mp.mp.dps = 80

    def characteristic(alpha):
        return mp.cosh(alpha) * mp.cos(alpha) - 1

    roots = []
    for mode in range(1, int(n_roots) + 1):
        roots.append(float(mp.findroot(characteristic, (mode + 0.5) * mp.pi)))
    return np.asarray(roots, dtype=float)


def validate_reference_roots(fp_alphas: np.ndarray, ff_alphas: np.ndarray) -> None:
    if not np.allclose(fp_alphas[: len(KNOWN_FP_PREFIX)], KNOWN_FP_PREFIX, rtol=0.0, atol=5e-4):
        raise RuntimeError(
            "Computed FP roots failed validation against known prefix: "
            f"computed={fp_alphas[:len(KNOWN_FP_PREFIX)]}, expected={KNOWN_FP_PREFIX}"
        )
    if not np.allclose(ff_alphas[: len(KNOWN_FF_PREFIX)], KNOWN_FF_PREFIX, rtol=0.0, atol=5e-11):
        raise RuntimeError(
            "Computed FF roots failed validation against known prefix: "
            f"computed={ff_alphas[:len(KNOWN_FF_PREFIX)]}, expected={KNOWN_FF_PREFIX}"
        )


FP_ALPHAS = clamped_pinned_roots(REFERENCE_ROOT_COUNT)
FF_ALPHAS = clamped_clamped_roots(REFERENCE_ROOT_COUNT)
validate_reference_roots(FP_ALPHAS, FF_ALPHAS)


def beam_params_from_epsilon(epsilon: float, l_base: float) -> BeamParams:
    return BeamParams(E=2.1e11, rho=7800.0, r=2.0 * float(epsilon) * float(l_base), L_total=2.0 * float(l_base))


def reference_lambdas(
    mu: float,
    *,
    fp_alphas: np.ndarray = FP_ALPHAS,
    ff_alphas: np.ndarray = FF_ALPHAS,
) -> dict[str, np.ndarray]:
    l2 = 1.0 + float(mu)
    return {
        "FP_finite": np.asarray(fp_alphas, dtype=float) / l2,
        "FF_finite": np.asarray(ff_alphas, dtype=float) / l2,
        "FP_limit": np.asarray(fp_alphas, dtype=float) / 2.0,
        "FF_limit": np.asarray(ff_alphas, dtype=float) / 2.0,
    }


def nearest_reference(Lambda: float, refs: np.ndarray) -> dict[str, float | int]:
    if not np.isfinite(Lambda):
        return {"mode": -1, "Lambda": float("nan"), "abs_diff": float("nan"), "rel_diff": float("nan")}
    diffs = np.abs(refs - float(Lambda))
    idx = int(np.nanargmin(diffs))
    ref = float(refs[idx])
    abs_diff = float(diffs[idx])
    return {
        "mode": idx + 1,
        "Lambda": ref,
        "abs_diff": abs_diff,
        "rel_diff": finite_rel_diff(float(Lambda), ref),
    }


def match_references(
    Lambda: float,
    mu: float,
    *,
    fp_alphas: np.ndarray = FP_ALPHAS,
    ff_alphas: np.ndarray = FF_ALPHAS,
) -> dict[str, dict[str, float | int] | str | float | int]:
    refs = reference_lambdas(mu, fp_alphas=fp_alphas, ff_alphas=ff_alphas)
    matches = {name: nearest_reference(float(Lambda), values) for name, values in refs.items()}
    best_type = min(matches, key=lambda name: float(matches[name]["rel_diff"]))
    best = matches[best_type]
    out: dict[str, dict[str, float | int] | str | float | int] = dict(matches)
    out["best_match_type"] = str(best_type)
    out["best_match_mode"] = int(best["mode"])
    out["best_match_abs_diff"] = float(best["abs_diff"])
    out["best_match_rel_diff"] = float(best["rel_diff"])
    return out


def root_row(*, beta: float, mu: float, epsilon: float, root_index: int, Lambda: float) -> dict[str, float | int | str]:
    match = match_references(float(Lambda), float(mu))
    row: dict[str, float | int | str] = {
        "beta": float(beta),
        "mu": float(mu),
        "epsilon": float(epsilon),
        "root_index": int(root_index),
        "Lambda_system": float(Lambda),
    }
    for key in ("FP_finite", "FF_finite", "FP_limit", "FF_limit"):
        info = match[key]
        assert isinstance(info, dict)
        row[f"nearest_{key}_mode"] = int(info["mode"])
        row[f"nearest_{key}_Lambda"] = float(info["Lambda"])
        row[f"abs_diff_{key}"] = float(info["abs_diff"])
        row[f"rel_diff_{key}"] = float(info["rel_diff"])
    row["best_match_type"] = str(match["best_match_type"])
    row["best_match_mode"] = int(match["best_match_mode"])
    row["best_match_abs_diff"] = float(match["best_match_abs_diff"])
    row["best_match_rel_diff"] = float(match["best_match_rel_diff"])
    return row


def solve_analytic_roots(*, beta: float, mu: float, epsilon: float, num_roots: int) -> np.ndarray:
    roots = find_first_n_roots(
        beta=float(np.deg2rad(beta)),
        mu=float(mu),
        eps=float(epsilon),
        n_roots=int(num_roots),
        Lmin=0.2,
        Lmax0=55.0,
        scan_step=0.01,
        grow_factor=1.35,
        max_tries=8,
    )
    return np.asarray(roots, dtype=float)


def _apply_fem_params(params: BeamParams) -> None:
    fem.r = params.r
    fem.E = params.E
    fem.rho = params.rho
    fem.L_tot = params.L_total
    fem.ell = params.L_total / 2.0
    fem.A = np.pi * params.r**2
    fem.I = np.pi * params.r**4 / 4.0
    fem.EI = fem.E * fem.I
    fem.EA = fem.E * fem.A
    fem.rhoA = fem.rho * fem.A
    fem.eps = np.sqrt(fem.I / fem.A) / fem.ell
    fem.scale = np.sqrt(fem.EI / fem.rhoA) / (2.0 * np.pi * fem.ell**2)
    fem.EI_nd = 1.0
    fem.rhoA_nd = 1.0
    fem.EA_nd = 1.0 / fem.eps**2


@contextmanager
def fem_parameter_override(params: BeamParams):
    saved = {name: getattr(fem, name) for name in FEM_STATE_KEYS}
    try:
        _apply_fem_params(params)
        yield
    finally:
        for name, value in saved.items():
            setattr(fem, name, value)


def solve_fem_lambdas(*, params: BeamParams, beta: float, mu: float, num_roots: int) -> np.ndarray:
    with fem_parameter_override(params):
        omega, _ = fem.fem_solve(mu=float(mu), beta_deg=float(beta), n_modes=int(num_roots))
    return np.sqrt(np.asarray(omega, dtype=float))


def collect_root_rows(*, betas: Sequence[float], mus: Sequence[float], epsilon: float, num_roots: int) -> tuple[list[dict[str, float | int | str]], dict[tuple[float, float], np.ndarray]]:
    rows: list[dict[str, float | int | str]] = []
    roots_by_case: dict[tuple[float, float], np.ndarray] = {}
    for beta in betas:
        for mu in mus:
            roots = solve_analytic_roots(beta=float(beta), mu=float(mu), epsilon=float(epsilon), num_roots=int(num_roots))
            roots_by_case[(float(beta), float(mu))] = roots
            for idx, Lambda in enumerate(roots, start=1):
                rows.append(root_row(beta=float(beta), mu=float(mu), epsilon=float(epsilon), root_index=idx, Lambda=float(Lambda)))
    return rows, roots_by_case


def collect_fem_rows(
    *,
    betas: Sequence[float],
    mus: Sequence[float],
    epsilon: float,
    l_base: float,
    num_roots: int,
    analytic_roots: dict[tuple[float, float], np.ndarray],
) -> list[dict[str, float | int | str]]:
    params = beam_params_from_epsilon(epsilon, l_base)
    rows: list[dict[str, float | int | str]] = []
    for beta in betas:
        for mu in mus:
            fem_lambdas = solve_fem_lambdas(params=params, beta=float(beta), mu=float(mu), num_roots=int(num_roots))
            analytic = analytic_roots.get((float(beta), float(mu)), np.full(int(num_roots), np.nan))
            for idx, Lambda in enumerate(fem_lambdas, start=1):
                match = match_references(float(Lambda), float(mu))
                same_index_analytic = float(analytic[idx - 1]) if idx - 1 < len(analytic) else float("nan")
                row: dict[str, float | int | str] = {
                    "beta": float(beta),
                    "mu": float(mu),
                    "epsilon": float(epsilon),
                    "root_index": int(idx),
                    "Lambda_fem": float(Lambda),
                    "Lambda_analytic_same_index": same_index_analytic,
                    "rel_diff_to_analytic_same_index": finite_rel_diff(float(Lambda), same_index_analytic),
                }
                for key in ("FP_finite", "FF_finite", "FP_limit", "FF_limit"):
                    info = match[key]
                    assert isinstance(info, dict)
                    row[f"nearest_{key}_mode"] = int(info["mode"])
                    row[f"nearest_{key}_Lambda"] = float(info["Lambda"])
                    row[f"abs_diff_{key}"] = float(info["abs_diff"])
                    row[f"rel_diff_{key}"] = float(info["rel_diff"])
                row["best_match_type"] = str(match["best_match_type"])
                row["best_match_mode"] = int(match["best_match_mode"])
                row["best_match_abs_diff"] = float(match["best_match_abs_diff"])
                row["best_match_rel_diff"] = float(match["best_match_rel_diff"])
                rows.append(row)
    return rows


def collect_branch_rows(*, epsilon: float, mus: Sequence[float], num_roots: int) -> tuple[list[dict[str, float | int | str]], list[str], str]:
    target_mus = [mu for mu in BRANCH_FOCUS_MUS if any(abs(float(mu) - float(item)) <= 1e-12 for item in mus)]
    if not target_mus:
        return [], ["Branch-focused mu values were not included in this run."], "skipped"

    n_track = max(4, max(int(branch.rsplit("_", 1)[1]) for branch in BRANCH_FOCUS_IDS))
    n_solve = max(int(num_roots), 12, n_track)
    common_kwargs = {
        "epsilon": float(epsilon),
        "beta": BRANCH_FOCUS_BETA,
        "mu_values": target_mus,
        "n_track": n_track,
        "n_solve": n_solve,
        "beta_steps": 31,
        "required_branch_ids": BRANCH_FOCUS_IDS,
    }

    warnings_out: list[str] = []
    status = "strict"
    try:
        tracking = track_mu_sweep(**common_kwargs, allow_low_mac=False)
    except Exception as exc:
        warnings_out.append(f"Strict branch tracking failed: {exc}")
        status = "allow_low_mac"
        try:
            tracking = track_mu_sweep(**common_kwargs, allow_low_mac=True)
        except Exception as fallback_exc:
            warnings_out.append(f"Exploratory branch tracking also failed: {fallback_exc}")
            return [], warnings_out, "failed"

    if tracking.warnings:
        warnings_out.extend(str(item) for item in tracking.warnings)

    rows: list[dict[str, float | int | str]] = []
    for mu in target_mus:
        for branch_id in BRANCH_FOCUS_IDS:
            point = tracking.point_at(branch_id, beta=BRANCH_FOCUS_BETA, mu=float(mu))
            match = match_references(float(point.Lambda), float(mu))
            row = {
                "beta": float(BRANCH_FOCUS_BETA),
                "mu": float(mu),
                "branch_id": str(branch_id),
                "current_sorted_index": int(point.current_sorted_index),
                "Lambda_branch": float(point.Lambda),
                "best_match_type": str(match["best_match_type"]),
                "best_match_mode": int(match["best_match_mode"]),
                "best_match_rel_diff": float(match["best_match_rel_diff"]),
                "rel_diff_FP_finite": float(match["FP_finite"]["rel_diff"]),  # type: ignore[index]
                "rel_diff_FF_finite": float(match["FF_finite"]["rel_diff"]),  # type: ignore[index]
                "rel_diff_FP_limit": float(match["FP_limit"]["rel_diff"]),  # type: ignore[index]
                "rel_diff_FF_limit": float(match["FF_limit"]["rel_diff"]),  # type: ignore[index]
            }
            rows.append(row)
    return rows, warnings_out, status


def write_rows(path: Path, rows: list[dict[str, float | int | str]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def beta_token(beta: float) -> str:
    if abs(float(beta) - round(float(beta))) <= 1e-12:
        integer = int(round(float(beta)))
        if 0 <= integer < 100:
            return f"{integer:02d}"
        return str(integer).replace("-", "m")
    return f"{float(beta):g}".replace("-", "m").replace(".", "p")


def plot_beta_roots(*, beta: float, mus: Sequence[float], roots_by_case: dict[tuple[float, float], np.ndarray], num_roots: int) -> Path:
    output = RESULTS_DIR / f"mu_to_one_single_rod_limit_beta{beta_token(beta)}.png"
    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    mu_array = np.asarray(mus, dtype=float)

    for idx in range(int(num_roots)):
        values = [float(roots_by_case[(float(beta), float(mu))][idx]) for mu in mu_array]
        ax.plot(mu_array, values, marker="o", linewidth=1.5, markersize=3.5, label=f"root {idx + 1}")

    for mode_idx, alpha in enumerate(FP_ALPHAS, start=1):
        ax.plot(mu_array, alpha / (1.0 + mu_array), color="black", linestyle="--", linewidth=0.9, alpha=0.45)
        if mode_idx == 1:
            ax.plot([], [], color="black", linestyle="--", linewidth=0.9, alpha=0.45, label="FP finite L2")
    for mode_idx, alpha in enumerate(FF_ALPHAS, start=1):
        ax.plot(mu_array, alpha / (1.0 + mu_array), color="tab:red", linestyle=":", linewidth=1.0, alpha=0.55)
        if mode_idx == 1:
            ax.plot([], [], color="tab:red", linestyle=":", linewidth=1.0, alpha=0.55, label="FF finite L2")

    ax.set_xlabel("mu")
    ax.set_ylabel("Lambda")
    ax.set_title(f"mu -> 1 single-rod limit diagnostic, beta={float(beta):g} deg")
    ax.grid(True, alpha=0.3)
    ax.legend(ncols=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output


def mean_of(rows: Sequence[dict[str, float | int | str]], key: str) -> float:
    values = [float(row[key]) for row in rows if np.isfinite(float(row[key]))]
    return float(np.mean(values)) if values else float("nan")


def count_by(rows: Sequence[dict[str, float | int | str]], key: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        label = str(row[key])
        counts[label] = counts.get(label, 0) + 1
    return counts


def count_finite_better(rows: Sequence[dict[str, float | int | str]]) -> tuple[int, int]:
    finite = 0
    limit = 0
    for row in rows:
        label = str(row["best_match_type"])
        if label.endswith("_finite"):
            finite += 1
        elif label.endswith("_limit"):
            limit += 1
    return finite, limit


def family_means(rows: Sequence[dict[str, float | int | str]]) -> dict[str, float]:
    return {
        "FP_finite": mean_of(rows, "rel_diff_FP_finite"),
        "FF_finite": mean_of(rows, "rel_diff_FF_finite"),
        "FP_limit": mean_of(rows, "rel_diff_FP_limit"),
        "FF_limit": mean_of(rows, "rel_diff_FF_limit"),
    }


def conclusion_text(rows: Sequence[dict[str, float | int | str]]) -> str:
    high_mu_rows = [row for row in rows if float(row["mu"]) >= 0.95]
    focus = high_mu_rows if high_mu_rows else list(rows)
    means = family_means(focus)
    finite_fp = means["FP_finite"]
    finite_ff = means["FF_finite"]
    finite_count, limit_count = count_finite_better(focus)
    if np.isfinite(finite_fp) and np.isfinite(finite_ff):
        if finite_fp < finite_ff:
            family = "FP"
        elif finite_ff < finite_fp:
            family = "FF"
        else:
            family = "unresolved"
    else:
        family = "unresolved"

    length_part = "finite-L2 references fit better than L=2 limit references"
    if limit_count > finite_count:
        length_part = "L=2 limit references fit better than finite-L2 references"
    elif limit_count == finite_count:
        length_part = "finite-L2 and L=2 limit references are tied by best-match counts"

    if family == "FP":
        return (
            "The sorted-root evidence is closer to the clamped-pinned (FP) single-rod family "
            f"over the high-mu subset, and {length_part}."
        )
    if family == "FF":
        return (
            "The sorted-root evidence is closer to the clamped-clamped (FF) single-rod family "
            f"over the high-mu subset, and {length_part}."
        )
    return f"The FP-vs-FF limit is unresolved by this run; {length_part}."


def markdown_counts(counts: dict[str, int]) -> str:
    if not counts:
        return "none"
    return ", ".join(f"{key}: {value}" for key, value in sorted(counts.items()))


def match_label(match: dict[str, dict[str, float | int] | str | float | int]) -> str:
    return f"{match['best_match_type']} {int(match['best_match_mode'])}"


def row_match_label(row: dict[str, float | int | str]) -> str:
    return f"{row['best_match_type']} {int(row['best_match_mode'])}"


def count_fp_ff(rows: Sequence[dict[str, float | int | str]]) -> tuple[int, int]:
    fp = sum(1 for row in rows if str(row["best_match_type"]).startswith("FP"))
    ff = sum(1 for row in rows if str(row["best_match_type"]).startswith("FF"))
    return fp, ff


def old_five_root_match(row: dict[str, float | int | str]) -> dict[str, dict[str, float | int] | str | float | int]:
    return match_references(
        float(row["Lambda_system"]),
        float(row["mu"]),
        fp_alphas=FP_ALPHAS[:5],
        ff_alphas=FF_ALPHAS[:5],
    )


def count_old_five_matches(rows: Sequence[dict[str, float | int | str]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        label = str(old_five_root_match(row)["best_match_type"])
        counts[label] = counts.get(label, 0) + 1
    return counts


def write_report(
    *,
    rows: list[dict[str, float | int | str]],
    fem_rows: list[dict[str, float | int | str]],
    branch_rows: list[dict[str, float | int | str]],
    branch_warnings: list[str],
    branch_status: str,
    plot_paths: Sequence[Path],
    args: argparse.Namespace,
) -> None:
    lines: list[str] = []
    lines.append("# mu -> 1 single-rod limit diagnostic")
    lines.append("")
    lines.append("Diagnostic only. No determinant, formulas.py, solvers.py, python_fem.py, or FEM physical model changes.")
    lines.append("")
    lines.append("## Configuration")
    lines.append("")
    lines.append(f"- epsilon: {float(args.epsilon):g}")
    lines.append(f"- betas: {', '.join(f'{float(beta):g}' for beta in args.betas)} deg")
    lines.append(f"- mus: {', '.join(f'{float(mu):g}' for mu in args.mus)}")
    lines.append(f"- num roots: {int(args.num_roots)}")
    lines.append(f"- include FEM: {bool(args.include_fem)}")
    lines.append(
        f"- FP roots: first {len(FP_ALPHAS)} roots computed from tan(alpha) = tanh(alpha); "
        f"prefix validation max abs error {float(np.max(np.abs(FP_ALPHAS[:len(KNOWN_FP_PREFIX)] - KNOWN_FP_PREFIX))):.6g}"
    )
    lines.append(
        f"- FF roots: first {len(FF_ALPHAS)} roots computed from cosh(alpha) cos(alpha) = 1; "
        f"prefix validation max abs error {float(np.max(np.abs(FF_ALPHAS[:len(KNOWN_FF_PREFIX)] - KNOWN_FF_PREFIX))):.6g}"
    )
    lines.append("- finite reference: Lambda = alpha / (1 + mu)")
    lines.append("- limiting reference: Lambda = alpha / 2")
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    lines.append(f"- analytic sorted-root CSV: `{ROOTS_CSV.relative_to(REPO_ROOT)}`")
    if fem_rows:
        lines.append(f"- FEM sorted-root CSV: `{FEM_ROOTS_CSV.relative_to(REPO_ROOT)}`")
    if branch_rows:
        lines.append(f"- branch-focused CSV: `{BRANCHES_CSV.relative_to(REPO_ROOT)}`")
    lines.extend(f"- plot: `{path.relative_to(REPO_ROOT)}`" for path in plot_paths)
    lines.append("")
    lines.append("## Sorted-root summary")
    lines.append("")
    lines.append(f"- best-match counts: {markdown_counts(count_by(rows, 'best_match_type'))}")
    finite_count, limit_count = count_finite_better(rows)
    lines.append(f"- finite-L2 best matches: {finite_count}")
    lines.append(f"- limiting-L=2 best matches: {limit_count}")
    means = family_means(rows)
    lines.append(
        "- mean relative differences over all sorted roots: "
        f"FP finite {means['FP_finite']:.6g}, FF finite {means['FF_finite']:.6g}, "
        f"FP limit {means['FP_limit']:.6g}, FF limit {means['FF_limit']:.6g}"
    )
    lines.append("")
    lines.append("| beta | best-match counts | finite vs limit | mean rel FP finite | mean rel FF finite |")
    lines.append("| --- | --- | --- | ---: | ---: |")
    for beta in args.betas:
        beta_rows = [row for row in rows if abs(float(row["beta"]) - float(beta)) <= 1e-12]
        beta_finite, beta_limit = count_finite_better(beta_rows)
        beta_means = family_means(beta_rows)
        lines.append(
            f"| {float(beta):g} | {markdown_counts(count_by(beta_rows, 'best_match_type'))} | "
            f"{beta_finite} / {beta_limit} | {beta_means['FP_finite']:.6g} | {beta_means['FF_finite']:.6g} |"
        )
    lines.append("")
    lines.append("## Effect of extending reference roots")
    lines.append("")
    high_rows = [row for row in rows if int(row["root_index"]) >= 6]
    high_label = f"roots 6-{int(args.num_roots)}" if int(args.num_roots) > 6 else "root 6"
    changed_high_rows = []
    for row in high_rows:
        old_match = old_five_root_match(row)
        if str(old_match["best_match_type"]) != str(row["best_match_type"]) or int(old_match["best_match_mode"]) != int(row["best_match_mode"]):
            changed_high_rows.append((row, old_match))
    low_mid_rows = [row for row in rows if int(row["root_index"]) <= 5]
    low_mid_fp, low_mid_ff = count_fp_ff(low_mid_rows)
    high_fp, high_ff = count_fp_ff(high_rows)
    lines.append(f"- {high_label} changed best reference in {len(changed_high_rows)} of {len(high_rows)} rows relative to the old first-5-only comparison.")
    lines.append(f"- old first-5-only best-match counts for {high_label}: {markdown_counts(count_old_five_matches(high_rows))}")
    lines.append(f"- extended-reference best-match counts for {high_label}: {markdown_counts(count_by(high_rows, 'best_match_type'))}")
    lines.append(f"- roots 1-5 remain predominantly FF by best-match family: FP {low_mid_fp}, FF {low_mid_ff}.")
    lines.append(f"- {high_label} after extension: FP {high_fp}, FF {high_ff}.")
    if changed_high_rows:
        lines.append("")
        lines.append("| beta | mu | root | Lambda | old first-5 best | extended best |")
        lines.append("| --- | --- | ---: | ---: | --- | --- |")
        for row, old_match in changed_high_rows[:16]:
            lines.append(
                f"| {float(row['beta']):g} | {float(row['mu']):g} | {int(row['root_index'])} | "
                f"{float(row['Lambda_system']):.8g} | {match_label(old_match)} | {row_match_label(row)} |"
            )
        if len(changed_high_rows) > 16:
            lines.append(f"| ... | ... | ... | ... | ... | {len(changed_high_rows) - 16} more changed rows |")
    fp_exception_rows = [row for row in rows if str(row["best_match_type"]).startswith("FP")]
    lines.append("")
    if fp_exception_rows:
        lines.append(f"Possible FP exceptions by best match: {len(fp_exception_rows)} rows.")
        lines.append("| beta | mu | root | Lambda | best match | rel diff |")
        lines.append("| --- | --- | ---: | ---: | --- | ---: |")
        for row in fp_exception_rows[:16]:
            lines.append(
                f"| {float(row['beta']):g} | {float(row['mu']):g} | {int(row['root_index'])} | "
                f"{float(row['Lambda_system']):.8g} | {row_match_label(row)} | {float(row['best_match_rel_diff']):.6g} |"
            )
        if len(fp_exception_rows) > 16:
            lines.append(f"| ... | ... | ... | ... | ... | {len(fp_exception_rows) - 16} more FP-best rows |")
    else:
        lines.append("No FP-best exceptions were found in the extended-reference run.")
    lines.append("")
    lines.append("| mu | finite-L2 best | limit-L=2 best | best-match counts |")
    lines.append("| --- | ---: | ---: | --- |")
    for mu in (0.9, 0.95, 0.99):
        mu_rows = [row for row in rows if abs(float(row["mu"]) - float(mu)) <= 1e-12]
        mu_finite, mu_limit = count_finite_better(mu_rows)
        lines.append(f"| {mu:g} | {mu_finite} | {mu_limit} | {markdown_counts(count_by(mu_rows, 'best_match_type'))} |")
    lines.append("")
    focus_rows = [
        row
        for row in rows
        if abs(float(row["beta"]) - BRANCH_FOCUS_BETA) <= 1e-12
        and any(abs(float(row["mu"]) - float(mu)) <= 1e-12 for mu in BRANCH_FOCUS_MUS)
    ]
    if focus_rows:
        focus_finite, focus_limit = count_finite_better(focus_rows)
        lines.append("## beta=15, mu=0.9/0.95 sorted-root focus")
        lines.append("")
        lines.append(f"- finite-L2 best matches: {focus_finite}")
        lines.append(f"- limiting-L=2 best matches: {focus_limit}")
        lines.append("| mu | root | Lambda | best match | rel FP finite | rel FF finite | rel FP limit | rel FF limit |")
        lines.append("| --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |")
        for row in focus_rows:
            lines.append(
                f"| {float(row['mu']):g} | {int(row['root_index'])} | {float(row['Lambda_system']):.8g} | "
                f"{row['best_match_type']} {int(row['best_match_mode'])} | "
                f"{float(row['rel_diff_FP_finite']):.6g} | {float(row['rel_diff_FF_finite']):.6g} | "
                f"{float(row['rel_diff_FP_limit']):.6g} | {float(row['rel_diff_FF_limit']):.6g} |"
            )
        lines.append("")
    lines.append("## Branch-focused check")
    lines.append("")
    if branch_rows:
        lines.append(f"- status: {branch_status}")
        lines.append("| beta | mu | branch_id | sorted index | Lambda | best match | rel diff |")
        lines.append("| --- | --- | --- | ---: | ---: | --- | ---: |")
        for row in branch_rows:
            lines.append(
                f"| {float(row['beta']):g} | {float(row['mu']):g} | {row['branch_id']} | "
                f"{int(row['current_sorted_index'])} | {float(row['Lambda_branch']):.8g} | "
                f"{row['best_match_type']} {int(row['best_match_mode'])} | "
                f"{float(row['best_match_rel_diff']):.6g} |"
            )
    else:
        lines.append(f"- branch-focused CSV was not generated; status: {branch_status}")
    if branch_warnings:
        lines.append("")
        lines.append("Branch-tracking warnings:")
        for warning in branch_warnings[:8]:
            lines.append(f"- {warning}")
        if len(branch_warnings) > 8:
            lines.append(f"- ... {len(branch_warnings) - 8} more warnings")
    lines.append("")
    lines.append("## FEM cross-check")
    lines.append("")
    if fem_rows:
        fem_rel = mean_of(fem_rows, "rel_diff_to_analytic_same_index")
        fem_max = max(
            [float(row["rel_diff_to_analytic_same_index"]) for row in fem_rows if np.isfinite(float(row["rel_diff_to_analytic_same_index"]))],
            default=float("nan"),
        )
        near_zero = sum(1 for row in fem_rows if abs(float(row["Lambda_fem"])) <= 1e-10)
        lines.append(
            f"FEM roots were computed with the current production FEM and compared by sorted index to analytic roots. "
            f"Mean relative analytic/FEM Lambda difference is {fem_rel:.6g}; max is {fem_max:.6g}."
        )
        if near_zero:
            lines.append(
                f"The FEM CSV contains {near_zero} near-zero Lambda entries in the highly asymmetric cases; "
                "treat the FEM sorted-index comparison as an auxiliary numerical check near mu=1."
            )
    else:
        lines.append("FEM roots were not computed in this run.")
    lines.append("")
    lines.append("## Conclusion")
    lines.append("")
    lines.append(conclusion_text(rows))
    if focus_rows:
        focus_finite, focus_limit = count_finite_better(focus_rows)
        if focus_limit > focus_finite:
            lines.append(
                "For beta=15 at mu=0.9/0.95, the limiting L=2 references are closer by best-match count than "
                "the finite-L2 references, despite the actual finite long-arm length being L2=1+mu."
            )
        elif focus_finite > focus_limit:
            lines.append(
                "For beta=15 at mu=0.9/0.95, the finite-L2 references are closer by best-match count than "
                "the limiting L=2 references."
            )
        else:
            lines.append(
                "For beta=15 at mu=0.9/0.95, finite-L2 and limiting-L=2 references are tied by best-match count."
            )
    lines.append(
        f"The reference set now covers {len(FP_ALPHAS)} FP and {len(FF_ALPHAS)} FF modes, so the requested "
        f"system roots 1-{int(args.num_roots)} can be compared without the old first-five truncation."
    )
    lines.append("")
    REPORT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Check the mu -> 1 coupled-rods frequency limit against FP and FF single-rod references."
    )
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--betas", nargs="+", type=float, default=DEFAULT_BETAS)
    parser.add_argument("--mus", nargs="+", type=float, default=DEFAULT_MUS)
    parser.add_argument("--num-roots", type=int, default=DEFAULT_NUM_ROOTS)
    parser.add_argument("--l-base", type=float, default=DEFAULT_L_BASE)
    fem_group = parser.add_mutually_exclusive_group()
    fem_group.add_argument("--include-fem", dest="include_fem", action="store_true")
    fem_group.add_argument("--no-include-fem", dest="include_fem", action="store_false")
    parser.set_defaults(include_fem=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    args.betas = parse_float_list(args.betas)
    args.mus = parse_float_list(args.mus)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    root_rows, roots_by_case = collect_root_rows(
        betas=args.betas,
        mus=args.mus,
        epsilon=float(args.epsilon),
        num_roots=int(args.num_roots),
    )
    write_rows(ROOTS_CSV, root_rows, ROOT_FIELDNAMES)

    fem_rows: list[dict[str, float | int | str]] = []
    if bool(args.include_fem):
        fem_rows = collect_fem_rows(
            betas=args.betas,
            mus=args.mus,
            epsilon=float(args.epsilon),
            l_base=float(args.l_base),
            num_roots=int(args.num_roots),
            analytic_roots=roots_by_case,
        )
        write_rows(FEM_ROOTS_CSV, fem_rows, FEM_FIELDNAMES)

    branch_rows, branch_warnings, branch_status = collect_branch_rows(
        epsilon=float(args.epsilon),
        mus=args.mus,
        num_roots=int(args.num_roots),
    )
    if branch_rows:
        write_rows(BRANCHES_CSV, branch_rows, BRANCH_FIELDNAMES)

    plot_paths = [
        plot_beta_roots(beta=float(beta), mus=args.mus, roots_by_case=roots_by_case, num_roots=int(args.num_roots))
        for beta in args.betas
    ]

    write_report(
        rows=root_rows,
        fem_rows=fem_rows,
        branch_rows=branch_rows,
        branch_warnings=branch_warnings,
        branch_status=branch_status,
        plot_paths=plot_paths,
        args=args,
    )

    print(f"saved analytic roots: {ROOTS_CSV}")
    if fem_rows:
        print(f"saved FEM roots: {FEM_ROOTS_CSV}")
    if branch_rows:
        print(f"saved branch roots: {BRANCHES_CSV} ({branch_status})")
    else:
        print(f"branch-focused roots not saved ({branch_status})")
    for path in plot_paths:
        print(f"saved plot: {path}")
    print(f"saved report: {REPORT_MD}")
    print(conclusion_text(root_rows))


if __name__ == "__main__":
    main()
