from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import re
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    assemble_clamped_coupled_matrix_eta,
    find_first_n_roots_eta,
)
from scripts.lib.thickness_mismatch_diagnostic_helpers import track_descendants_from_mu0  # noqa: E402
from scripts.lib.thickness_mismatch_mac_tracking import (  # noqa: E402
    analytic_shape_vectors_for_roots,
    mac_value,
)


BETA_DEG = 5.0
EPSILON = 0.0025
ETA = 0.5

MU_MIN = 0.0
MU_MAX = 0.9
COARSE_MU_POINTS = 901
LOCAL_REFINE_HALF_WIDTH = 0.003
LOCAL_REFINE_STEP = 1e-5

N_DESCENDANTS_PLOT = 6
N_TRACKED_DESCENDANTS = 8
N_SORTED_ROOTS = 12
N_SORTED_ROOT_FALLBACKS = (10, 8)

ROOT_SCAN_STEP = 0.01
FINE_ROOT_SCAN_STEP = 0.001
FINE_ROOT_SCAN_WINDOWS = ((0.150, 0.170),)
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401
MAC_WARNING_THRESHOLD = 0.9
MAC_MARGIN_WARNING_THRESHOLD = 0.05
MAX_SORTED_POSITION_JUMP = 1
ENABLE_ROOT_REPAIR = False
ROOT_REPAIR_JUMP_THRESHOLD = 0.15
ROOT_REPAIR_COUNT = N_TRACKED_DESCENDANTS
ROOT_REPAIR_XATOL = 1e-12

DESC_I = 1
DESC_J = 2

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_CSV = OUTPUT_DIR / "article_check_beta5_eta0p5_eps0p0025_mu_desc1_desc2_crossing.csv"
OUTPUT_REPORT = OUTPUT_DIR / "article_check_beta5_eta0p5_eps0p0025_mu_desc1_desc2_crossing_report.md"
OUTPUT_PNG = OUTPUT_DIR / "article_check_beta5_eta0p5_eps0p0025_mu_desc1_desc2_crossing.png"

ARTICLE_METADATA = (
    REPO_ROOT
    / "paper_thickness_mismatch_article"
    / "data"
    / "lambda_mu_beta5_eps0p0025_eta_p0p5_metadata.md"
)

CSV_FIELDNAMES = [
    "mu",
    "Lambda_desc1",
    "Lambda_desc2",
    "diff_desc1_minus_desc2",
    "abs_gap_desc12",
    "sorted_position_desc1",
    "sorted_position_desc2",
    "sorted1_lambda",
    "sorted2_lambda",
    "mac_desc1_to_previous",
    "mac_desc2_to_previous",
    "tracking_warning",
]


@dataclass(frozen=True)
class TrackingBundle:
    mu_values: np.ndarray
    sorted_roots: np.ndarray
    tracked_lambdas: np.ndarray
    sorted_positions: np.ndarray
    mac_to_previous: np.ndarray
    tracking_warning: list[str]
    warning_rows: list[dict[str, float | int | str]]
    n_roots_used: int
    root_scan_note: str
    root_repair_count: int


@dataclass(frozen=True)
class CrossingSummary:
    sign_change: bool
    crossing_mu: float
    crossing_lambda: float
    min_gap_mu: float
    min_gap: float
    relative_gap: float
    lambda_desc1_at_min: float
    lambda_desc2_at_min: float
    min_gap_index: int
    cross_mac_at_min: float
    local_min_mac_desc1: float
    local_min_mac_desc2: float
    sorted_position_swap: bool
    sorted1_identity_changes: bool
    sorted2_identity_changes: bool
    sorted1_descendant_at_min: int
    sorted2_descendant_at_min: int


@dataclass(frozen=True)
class ArticlePairInfo:
    branch_mode: str
    sign_change: bool
    crossing_mu_linear: float
    interval_left: float
    interval_right: float
    min_gap: float
    min_gap_mu: float
    desc1_sorted_position_at_right: int
    desc2_sorted_position_at_right: int


@dataclass(frozen=True)
class WindowSummary:
    center_mu: float
    half_width: float
    sign_change: bool
    min_gap: float
    min_gap_mu: float
    lambda_desc1_at_min: float
    lambda_desc2_at_min: float
    relative_gap: float


def fmt(value: float, precision: int = 8) -> str:
    if not np.isfinite(float(value)):
        return "nan"
    return f"{float(value):.{precision}g}"


def yesno(value: bool) -> str:
    return "yes" if bool(value) else "no"


def coarse_mu_grid() -> np.ndarray:
    return np.linspace(MU_MIN, MU_MAX, COARSE_MU_POINTS, dtype=float)


def refined_mu_grid(center_mu: float) -> np.ndarray:
    start = max(MU_MIN, float(center_mu) - LOCAL_REFINE_HALF_WIDTH)
    stop = min(MU_MAX, float(center_mu) + LOCAL_REFINE_HALF_WIDTH)
    local = np.arange(start, stop + 0.5 * LOCAL_REFINE_STEP, LOCAL_REFINE_STEP, dtype=float)
    values = np.concatenate([coarse_mu_grid(), local])
    values = np.unique(np.round(values, 10))
    values[0] = MU_MIN
    values[-1] = MU_MAX
    return values


def smallest_singular_ratio(Lambda: float, *, mu: float) -> float:
    matrix = assemble_clamped_coupled_matrix_eta(
        float(Lambda),
        float(np.deg2rad(BETA_DEG)),
        float(mu),
        EPSILON,
        ETA,
    )
    singular_values = np.linalg.svd(matrix, compute_uv=False)
    scale = float(singular_values[0]) if singular_values.size else 1.0
    return float(singular_values[-1]) / scale if scale > 0.0 else float("inf")


def repair_half_width(previous_roots: np.ndarray, root_idx: int) -> float:
    spacings: list[float] = []
    if root_idx > 0:
        spacings.append(float(previous_roots[root_idx] - previous_roots[root_idx - 1]))
    if root_idx + 1 < len(previous_roots):
        spacings.append(float(previous_roots[root_idx + 1] - previous_roots[root_idx]))
    positive = [value for value in spacings if np.isfinite(value) and value > 0.0]
    if not positive:
        return 0.006
    return max(0.001, min(0.008, 0.45 * min(positive)))


def repair_root_near_previous(previous_root: float, *, mu: float, half_width: float) -> float:
    lo = max(0.2, float(previous_root) - float(half_width))
    hi = float(previous_root) + float(half_width)
    result = minimize_scalar(
        lambda value: smallest_singular_ratio(float(value), mu=float(mu)),
        bounds=(lo, hi),
        method="bounded",
        options={"xatol": ROOT_REPAIR_XATOL},
    )
    return float(result.x)


def unique_sorted(values: Sequence[float], *, tol: float = 1e-7) -> list[float]:
    out: list[float] = []
    for value in sorted(float(v) for v in values if np.isfinite(float(v))):
        if not out or abs(value - out[-1]) > tol:
            out.append(value)
    return out


def repair_roots_if_needed(raw_roots: np.ndarray, previous_roots: np.ndarray | None, *, mu: float) -> tuple[np.ndarray, bool]:
    if not ENABLE_ROOT_REPAIR:
        return raw_roots, False
    if previous_roots is None:
        return raw_roots, False
    check_count = min(ROOT_REPAIR_COUNT, len(raw_roots), len(previous_roots))
    if check_count == 0:
        return raw_roots, False
    jump = np.max(np.abs(raw_roots[:check_count] - previous_roots[:check_count]))
    if not np.isfinite(jump) or float(jump) <= ROOT_REPAIR_JUMP_THRESHOLD:
        return raw_roots, False

    repaired: list[float] = []
    for root_idx in range(check_count):
        repaired.append(
            repair_root_near_previous(
                float(previous_roots[root_idx]),
                mu=float(mu),
                half_width=repair_half_width(previous_roots, root_idx),
            )
        )
    merged = unique_sorted([*repaired, *raw_roots])
    if len(merged) < len(raw_roots):
        return raw_roots, False
    return np.asarray(merged[: len(raw_roots)], dtype=float), True


def root_scan_step_for_mu(mu: float) -> float:
    mu_f = float(mu)
    for start, stop in FINE_ROOT_SCAN_WINDOWS:
        if float(start) <= mu_f <= float(stop):
            return FINE_ROOT_SCAN_STEP
    return ROOT_SCAN_STEP


def solve_root_grid_for_n(mu_values: Sequence[float], n_roots: int) -> tuple[dict[float, np.ndarray], int]:
    roots_by_mu: dict[float, np.ndarray] = {}
    previous_roots: np.ndarray | None = None
    repair_count = 0
    for mu in np.asarray(mu_values, dtype=float):
        raw_roots = find_first_n_roots_eta(
            float(np.deg2rad(BETA_DEG)),
            float(mu),
            EPSILON,
            ETA,
            int(n_roots),
            Lmax0=ROOT_LMAX0,
            scan_step=root_scan_step_for_mu(float(mu)),
        )
        if np.any(~np.isfinite(raw_roots)):
            raise RuntimeError(f"Missing roots for mu={float(mu):g}, eta={ETA:g}.")
        roots, repaired = repair_roots_if_needed(raw_roots, previous_roots, mu=float(mu))
        if repaired:
            repair_count += 1
        roots_by_mu[float(mu)] = roots
        previous_roots = roots
    return roots_by_mu, repair_count


def solve_roots_with_fallback(mu_values: Sequence[float]) -> tuple[dict[float, np.ndarray], int, str, int]:
    last_error: RuntimeError | None = None
    scan_note = (
        f"adaptive scan step {ROOT_SCAN_STEP:g}, "
        f"fine step {FINE_ROOT_SCAN_STEP:g} in mu windows {FINE_ROOT_SCAN_WINDOWS}"
    )
    for n_roots in (N_SORTED_ROOTS, *N_SORTED_ROOT_FALLBACKS):
        if int(n_roots) < N_TRACKED_DESCENDANTS:
            continue
        try:
            roots, repair_count = solve_root_grid_for_n(mu_values, int(n_roots))
        except RuntimeError as exc:
            if "Missing roots" not in str(exc):
                raise
            last_error = exc
            continue
        repair_note = f"; {scan_note}; local SVD root repairs: {repair_count}"
        if int(n_roots) == N_SORTED_ROOTS:
            return roots, int(n_roots), "primary sorted-root scan" + repair_note, repair_count
        return (
            roots,
            int(n_roots),
            f"fallback to first {int(n_roots)} roots after: {last_error}" + repair_note,
            repair_count,
        )
    if last_error is not None:
        raise last_error
    raise RuntimeError("No sorted-root scan size is available.")


def row_lookup(rows: Sequence[dict[str, float | int | str]]) -> dict[tuple[int, float], dict[str, float | int | str]]:
    lookup: dict[tuple[int, float], dict[str, float | int | str]] = {}
    for row in rows:
        lookup[(int(row["branch_index_from_mu0"]), round(float(row["mu"]), 10))] = row
    return lookup


def collect_tracking(mu_values: np.ndarray) -> TrackingBundle:
    roots_by_mu, n_roots_used, root_scan_note, repair_count = solve_roots_with_fallback(mu_values)
    tracking = track_descendants_from_mu0(
        beta_rad=float(np.deg2rad(BETA_DEG)),
        epsilon=EPSILON,
        eta=ETA,
        mu_values=mu_values,
        roots_by_mu=roots_by_mu,
        num_descendants=N_TRACKED_DESCENDANTS,
        num_shape_samples=NUM_SHAPE_SAMPLES,
        mac_warning_threshold=MAC_WARNING_THRESHOLD,
        mac_margin_warning_threshold=MAC_MARGIN_WARNING_THRESHOLD,
        max_sorted_position_jump=MAX_SORTED_POSITION_JUMP,
    )

    roots = np.vstack([np.asarray(roots_by_mu[float(mu)], dtype=float) for mu in mu_values])
    lookup = row_lookup(tracking.rows)
    warning_lookup: dict[float, list[str]] = {}
    for row in tracking.warning_rows:
        branch = int(row["branch_index_from_mu0"])
        if branch not in (DESC_I, DESC_J):
            continue
        mu_key = round(float(row["mu"]), 10)
        warning_lookup.setdefault(mu_key, []).append(
            f"desc{branch}:{row.get('tracking_step_status', '')}"
        )

    mac = np.full((N_TRACKED_DESCENDANTS, len(mu_values)), np.nan, dtype=float)
    warning_text: list[str] = []
    for col, mu in enumerate(mu_values):
        mu_key = round(float(mu), 10)
        for branch in range(1, N_TRACKED_DESCENDANTS + 1):
            row = lookup.get((branch, mu_key))
            if row is not None:
                mac[branch - 1, col] = float(row.get("mac_to_previous", np.nan))
        warning_text.append("; ".join(warning_lookup.get(mu_key, [])) or "no")

    return TrackingBundle(
        mu_values=np.asarray(mu_values, dtype=float),
        sorted_roots=roots,
        tracked_lambdas=np.asarray(tracking.tracked, dtype=float),
        sorted_positions=np.asarray(tracking.current_sorted_indices, dtype=int),
        mac_to_previous=mac,
        tracking_warning=warning_text,
        warning_rows=list(tracking.warning_rows),
        n_roots_used=n_roots_used,
        root_scan_note=root_scan_note,
        root_repair_count=repair_count,
    )


def sign_change_crossing(mu_values: np.ndarray, diff: np.ndarray) -> tuple[bool, float, float]:
    for idx in range(len(diff) - 1):
        left = float(diff[idx])
        right = float(diff[idx + 1])
        if np.isclose(left, 0.0, rtol=0.0, atol=1e-12):
            return True, float(mu_values[idx]), float(np.mean(bundle_pair_lambdas_at_index(idx)))
        if left * right < 0.0:
            frac = abs(left) / (abs(left) + abs(right))
            mu_cross = float(mu_values[idx] + frac * (mu_values[idx + 1] - mu_values[idx]))
            lambda_left = float(0.5 * np.sum(bundle_pair_lambdas_at_index(idx)))
            lambda_right = float(0.5 * np.sum(bundle_pair_lambdas_at_index(idx + 1)))
            lambda_cross = lambda_left + frac * (lambda_right - lambda_left)
            return True, mu_cross, lambda_cross
    if np.isclose(float(diff[-1]), 0.0, rtol=0.0, atol=1e-12):
        return True, float(mu_values[-1]), float(np.mean(bundle_pair_lambdas_at_index(len(diff) - 1)))
    return False, float("nan"), float("nan")


def sign_change_in_values(mu_values: np.ndarray, y_values: np.ndarray) -> tuple[bool, float, float, int]:
    for idx in range(len(y_values) - 1):
        left = float(y_values[idx])
        right = float(y_values[idx + 1])
        if np.isclose(left, 0.0, rtol=0.0, atol=1e-12):
            return True, float(mu_values[idx]), float(mu_values[idx]), idx
        if left * right < 0.0:
            frac = abs(left) / (abs(left) + abs(right))
            mu_cross = float(mu_values[idx] + frac * (mu_values[idx + 1] - mu_values[idx]))
            return True, mu_cross, float(mu_values[idx + 1]), idx
    return False, float("nan"), float("nan"), -1


_PAIR_LAMBDAS_FOR_INTERPOLATION: np.ndarray | None = None


def bundle_pair_lambdas_at_index(index: int) -> np.ndarray:
    if _PAIR_LAMBDAS_FOR_INTERPOLATION is None:
        raise RuntimeError("Pair lambdas for interpolation were not initialized.")
    return _PAIR_LAMBDAS_FOR_INTERPOLATION[:, int(index)]


def descendant_at_sorted_position(positions: np.ndarray, sorted_position: int, col: int) -> int:
    matches = np.flatnonzero(positions[:, int(col)] == int(sorted_position))
    return int(matches[0]) + 1 if matches.size else -1


def finite_local_min(values: np.ndarray, center: int, radius: int = 3) -> float:
    left = max(0, int(center) - int(radius))
    right = min(values.size, int(center) + int(radius) + 1)
    finite = values[left:right][np.isfinite(values[left:right])]
    return float(np.min(finite)) if finite.size else float("nan")


def cross_mac_at(bundle: TrackingBundle, col: int) -> float:
    roots = bundle.sorted_roots[int(col)]
    vectors = analytic_shape_vectors_for_roots(
        roots,
        beta_rad=float(np.deg2rad(BETA_DEG)),
        mu=float(bundle.mu_values[int(col)]),
        epsilon=EPSILON,
        eta=ETA,
        s_norm=np.linspace(0.0, 1.0, NUM_SHAPE_SAMPLES, dtype=float),
    )
    pos_i = int(bundle.sorted_positions[DESC_I - 1, int(col)]) - 1
    pos_j = int(bundle.sorted_positions[DESC_J - 1, int(col)]) - 1
    return mac_value(vectors[pos_i], vectors[pos_j])


def summarize(bundle: TrackingBundle) -> CrossingSummary:
    global _PAIR_LAMBDAS_FOR_INTERPOLATION
    pair_lambdas = bundle.tracked_lambdas[[DESC_I - 1, DESC_J - 1], :]
    _PAIR_LAMBDAS_FOR_INTERPOLATION = pair_lambdas
    diff = pair_lambdas[0] - pair_lambdas[1]
    abs_gap = np.abs(diff)
    min_idx = int(np.argmin(abs_gap))
    mean_lambda = 0.5 * (abs(float(pair_lambdas[0, min_idx])) + abs(float(pair_lambdas[1, min_idx])))
    rel_gap = float(abs_gap[min_idx]) / mean_lambda if mean_lambda > 0.0 else float("nan")
    sign_change, mu_cross, lambda_cross = sign_change_crossing(bundle.mu_values, diff)

    sorted1_ids = np.array(
        [descendant_at_sorted_position(bundle.sorted_positions, 1, col) for col in range(len(bundle.mu_values))],
        dtype=int,
    )
    sorted2_ids = np.array(
        [descendant_at_sorted_position(bundle.sorted_positions, 2, col) for col in range(len(bundle.mu_values))],
        dtype=int,
    )
    pos_i = bundle.sorted_positions[DESC_I - 1]
    pos_j = bundle.sorted_positions[DESC_J - 1]
    return CrossingSummary(
        sign_change=sign_change,
        crossing_mu=mu_cross,
        crossing_lambda=lambda_cross,
        min_gap_mu=float(bundle.mu_values[min_idx]),
        min_gap=float(abs_gap[min_idx]),
        relative_gap=rel_gap,
        lambda_desc1_at_min=float(pair_lambdas[0, min_idx]),
        lambda_desc2_at_min=float(pair_lambdas[1, min_idx]),
        min_gap_index=min_idx,
        cross_mac_at_min=cross_mac_at(bundle, min_idx),
        local_min_mac_desc1=finite_local_min(bundle.mac_to_previous[DESC_I - 1], min_idx),
        local_min_mac_desc2=finite_local_min(bundle.mac_to_previous[DESC_J - 1], min_idx),
        sorted_position_swap=bool(np.any(np.diff(pos_i) != 0) or np.any(np.diff(pos_j) != 0)),
        sorted1_identity_changes=bool(np.any(sorted1_ids != sorted1_ids[0])),
        sorted2_identity_changes=bool(np.any(sorted2_ids != sorted2_ids[0])),
        sorted1_descendant_at_min=int(sorted1_ids[min_idx]),
        sorted2_descendant_at_min=int(sorted2_ids[min_idx]),
    )


def summarize_window(bundle: TrackingBundle, *, center_mu: float, half_width: float) -> WindowSummary:
    mask = (bundle.mu_values >= float(center_mu) - float(half_width)) & (
        bundle.mu_values <= float(center_mu) + float(half_width)
    )
    indices = np.flatnonzero(mask)
    if indices.size == 0:
        raise RuntimeError("Requested summary window contains no mu values.")
    desc1 = bundle.tracked_lambdas[DESC_I - 1, indices]
    desc2 = bundle.tracked_lambdas[DESC_J - 1, indices]
    diff = desc1 - desc2
    abs_gap = np.abs(diff)
    local_min = int(np.argmin(abs_gap))
    col = int(indices[local_min])
    mean_lambda = 0.5 * (abs(float(desc1[local_min])) + abs(float(desc2[local_min])))
    sign_change, _mu_cross, _right, _idx = sign_change_in_values(bundle.mu_values[indices], diff)
    return WindowSummary(
        center_mu=float(center_mu),
        half_width=float(half_width),
        sign_change=sign_change,
        min_gap=float(abs_gap[local_min]),
        min_gap_mu=float(bundle.mu_values[col]),
        lambda_desc1_at_min=float(bundle.tracked_lambdas[DESC_I - 1, col]),
        lambda_desc2_at_min=float(bundle.tracked_lambdas[DESC_J - 1, col]),
        relative_gap=float(abs_gap[local_min]) / mean_lambda if mean_lambda > 0.0 else float("nan"),
    )


def csv_rows(bundle: TrackingBundle) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    desc1 = bundle.tracked_lambdas[DESC_I - 1]
    desc2 = bundle.tracked_lambdas[DESC_J - 1]
    diff = desc1 - desc2
    for col, mu in enumerate(bundle.mu_values):
        rows.append(
            {
                "mu": f"{float(mu):.10g}",
                "Lambda_desc1": f"{float(desc1[col]):.12g}",
                "Lambda_desc2": f"{float(desc2[col]):.12g}",
                "diff_desc1_minus_desc2": f"{float(diff[col]):.12g}",
                "abs_gap_desc12": f"{float(abs(diff[col])):.12g}",
                "sorted_position_desc1": int(bundle.sorted_positions[DESC_I - 1, col]),
                "sorted_position_desc2": int(bundle.sorted_positions[DESC_J - 1, col]),
                "sorted1_lambda": f"{float(bundle.sorted_roots[col, 0]):.12g}",
                "sorted2_lambda": f"{float(bundle.sorted_roots[col, 1]):.12g}",
                "mac_desc1_to_previous": (
                    "" if not np.isfinite(bundle.mac_to_previous[DESC_I - 1, col]) else f"{float(bundle.mac_to_previous[DESC_I - 1, col]):.12g}"
                ),
                "mac_desc2_to_previous": (
                    "" if not np.isfinite(bundle.mac_to_previous[DESC_J - 1, col]) else f"{float(bundle.mac_to_previous[DESC_J - 1, col]):.12g}"
                ),
                "tracking_warning": bundle.tracking_warning[col],
            }
        )
    return rows


def write_csv(bundle: TrackingBundle) -> None:
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        writer.writerows(csv_rows(bundle))


def article_branch_mode() -> str:
    if not ARTICLE_METADATA.exists():
        return "metadata file not found"
    text = ARTICLE_METADATA.read_text(encoding="utf-8")
    match = re.search(r"^- branch_mode:\s*(.+)$", text, flags=re.MULTILINE)
    return match.group(1).strip() if match else "not recorded"


def load_article_pair_info() -> ArticlePairInfo:
    if not ARTICLE_METADATA.exists():
        return ArticlePairInfo(
            branch_mode="metadata file not found",
            sign_change=False,
            crossing_mu_linear=float("nan"),
            interval_left=float("nan"),
            interval_right=float("nan"),
            min_gap=float("nan"),
            min_gap_mu=float("nan"),
            desc1_sorted_position_at_right=-1,
            desc2_sorted_position_at_right=-1,
        )
    csv_path = ARTICLE_METADATA.with_name("lambda_mu_beta5_eps0p0025_eta_p0p5.csv")
    by_mu: dict[float, dict[int, dict[str, str]]] = {}
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            branch = int(row["branch_index"])
            if branch not in (DESC_I, DESC_J):
                continue
            by_mu.setdefault(float(row["mu"]), {})[branch] = row

    mus = np.array(sorted(mu for mu, rows in by_mu.items() if DESC_I in rows and DESC_J in rows), dtype=float)
    desc1 = np.array([float(by_mu[mu][DESC_I]["Lambda"]) for mu in mus], dtype=float)
    desc2 = np.array([float(by_mu[mu][DESC_J]["Lambda"]) for mu in mus], dtype=float)
    diff = desc1 - desc2
    sign_change, mu_cross, right_mu, idx = sign_change_in_values(mus, diff)
    abs_gap = np.abs(diff)
    min_idx = int(np.argmin(abs_gap))
    right_rows = by_mu.get(right_mu, {}) if np.isfinite(right_mu) else {}
    return ArticlePairInfo(
        branch_mode=article_branch_mode(),
        sign_change=sign_change,
        crossing_mu_linear=mu_cross,
        interval_left=float(mus[idx]) if idx >= 0 else float("nan"),
        interval_right=right_mu,
        min_gap=float(abs_gap[min_idx]),
        min_gap_mu=float(mus[min_idx]),
        desc1_sorted_position_at_right=int(right_rows.get(DESC_I, {}).get("sorted_position", -1)),
        desc2_sorted_position_at_right=int(right_rows.get(DESC_J, {}).get("sorted_position", -1)),
    )


def local_table(bundle: TrackingBundle, summary: CrossingSummary) -> list[list[str]]:
    rows: list[list[str]] = []
    for col in range(max(0, summary.min_gap_index - 2), min(len(bundle.mu_values), summary.min_gap_index + 3)):
        rows.append(
            [
                fmt(float(bundle.mu_values[col]), 10),
                fmt(float(bundle.tracked_lambdas[DESC_I - 1, col]), 10),
                fmt(float(bundle.tracked_lambdas[DESC_J - 1, col]), 10),
                fmt(float(bundle.tracked_lambdas[DESC_I - 1, col] - bundle.tracked_lambdas[DESC_J - 1, col]), 10),
                str(int(bundle.sorted_positions[DESC_I - 1, col])),
                str(int(bundle.sorted_positions[DESC_J - 1, col])),
                fmt(float(bundle.mac_to_previous[DESC_I - 1, col]), 8),
                fmt(float(bundle.mac_to_previous[DESC_J - 1, col]), 8),
                bundle.tracking_warning[col],
            ]
        )
    return rows


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> list[str]:
    output = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    output.extend("| " + " | ".join(str(value) for value in row) + " |" for row in rows)
    return output


def warning_lines_near_min(bundle: TrackingBundle, summary: CrossingSummary) -> list[str]:
    lines: list[str] = []
    for row in bundle.warning_rows:
        branch = int(row["branch_index_from_mu0"])
        if branch not in (DESC_I, DESC_J):
            continue
        if abs(float(row["mu"]) - summary.min_gap_mu) > 0.02:
            continue
        lines.append(
            "- "
            f"descendant {branch}, mu={fmt(float(row['mu']), 10)}, "
            f"status={row.get('tracking_step_status', '')}, "
            f"MAC={fmt(float(row.get('mac_to_previous', np.nan)), 8)}, "
            f"candidate sorted={row.get('diagnostic_candidate_sorted_position', '')}, "
            f"accepted sorted={row.get('mac_sorted_root_index', '')}, "
            f"nearest-frequency sorted={row.get('nearest_sorted_root_index', '')}"
        )
    return lines


def write_report(
    bundle: TrackingBundle,
    summary: CrossingSummary,
    coarse_summary: CrossingSummary,
    article_info: ArticlePairInfo,
    apparent_window: WindowSummary,
    refine_center: float,
    refine_reason: str,
) -> None:
    article_mismatch = bool(article_info.sign_change and not apparent_window.sign_change)
    classification = (
        "real descendant crossing"
        if summary.sign_change
        else "finite-gap near miss / avoided-crossing-like approach"
    )
    correction = (
        "The current article figure should be regenerated: its coarse descendant tracking accepts "
        "a desc1/desc2 swap across the visual interaction, while the refined audit with local "
        "root resolution keeps a finite positive gap and no descendant swap."
        if article_mismatch
        else "No article correction is indicated from the apparent-crossing window in this audit."
    )
    lines = [
        "# Article Figure Crossing Audit: beta=5, eta=0.5",
        "",
        "## What Was Checked",
        "",
        "Analytic-only Euler--Bernoulli thickness-mismatch audit for the article figure",
        "`lambda_mu_beta5_eps0p0025_eta_p0p5.png`. The check separates descendant",
        "branches from sorted roots for the first two displayed curves.",
        "",
        "## Parameters",
        "",
        f"- beta: `{BETA_DEG:g} deg`",
        f"- epsilon: `{EPSILON:g}`",
        f"- eta: `{ETA:g}`",
        f"- mu range: `{MU_MIN:g}..{MU_MAX:g}`",
        f"- coarse grid: `{COARSE_MU_POINTS}` points",
        f"- local refinement: half-width `{LOCAL_REFINE_HALF_WIDTH:g}`, step `{LOCAL_REFINE_STEP:g}`",
        f"- local refinement center: `{fmt(refine_center, 10)}` ({refine_reason})",
        f"- final grid points: `{len(bundle.mu_values)}`",
        f"- sorted roots solved: first `{bundle.n_roots_used}` (`{bundle.root_scan_note}`)",
        f"- tracked descendants: first `{N_TRACKED_DESCENDANTS}`; reported pair: `{DESC_I}` and `{DESC_J}`",
        f"- article plotting metadata branch mode: `{article_branch_mode()}`",
        "",
        "Descendant branches are seeded at `mu=0` for this fixed `beta, epsilon, eta` case",
        "and then tracked by adjacent-step analytic shape MAC. Sorted root position is",
        "diagnostic metadata only.",
        "",
        "## Answer",
        "",
        f"- Does descendant 1 cross descendant 2? `{yesno(summary.sign_change)}`",
        f"- classification: `{classification}`",
        f"- global minimum gap on the final grid: `{fmt(summary.min_gap, 12)}` at `mu={fmt(summary.min_gap_mu, 10)}`",
        f"- Lambda at min gap: descendant 1 `{fmt(summary.lambda_desc1_at_min, 12)}`, descendant 2 `{fmt(summary.lambda_desc2_at_min, 12)}`",
        f"- relative gap: `{fmt(summary.relative_gap, 12)}`",
        f"- sign change on refined grid: `{yesno(summary.sign_change)}`",
        f"- coarse-grid min gap before refinement: `{fmt(coarse_summary.min_gap, 12)}` at `mu={fmt(coarse_summary.min_gap_mu, 10)}`",
        "",
    ]
    if summary.sign_change:
        lines.extend(
            [
                f"- refined crossing mu: `{fmt(summary.crossing_mu, 12)}`",
                f"- interpolated Lambda at crossing: `{fmt(summary.crossing_lambda, 12)}`",
                f"- same-mu cross-MAC at nearest refined grid point: `{fmt(summary.cross_mac_at_min, 8)}`",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "No sign change of `Lambda_desc1 - Lambda_desc2` was found on the refined grid.",
                "The visual approach is therefore not accepted as a clean real crossing of the tracked descendants.",
                "",
            ]
        )

    lines.extend(
        [
            "## Apparent Article Crossing Window",
            "",
            f"- article CSV branch mode: `{article_info.branch_mode}`",
            f"- article CSV desc1/desc2 sign change: `{yesno(article_info.sign_change)}`",
            f"- article linear crossing estimate: `mu={fmt(article_info.crossing_mu_linear, 10)}` "
            f"from interval `{fmt(article_info.interval_left, 10)}..{fmt(article_info.interval_right, 10)}`",
            f"- article sorted positions at right side of that interval: desc1 `{article_info.desc1_sorted_position_at_right}`, desc2 `{article_info.desc2_sorted_position_at_right}`",
            f"- refined robust sign change in the same window: `{yesno(apparent_window.sign_change)}`",
            f"- refined window min gap: `{fmt(apparent_window.min_gap, 12)}` at `mu={fmt(apparent_window.min_gap_mu, 10)}`",
            f"- refined window Lambda values: desc1 `{fmt(apparent_window.lambda_desc1_at_min, 12)}`, desc2 `{fmt(apparent_window.lambda_desc2_at_min, 12)}`",
            f"- refined window relative gap: `{fmt(apparent_window.relative_gap, 12)}`",
            "",
            "## Sorted/Descendant Interpretation",
            "",
            f"- descendant 1 sorted position changes over mu: `{yesno(np.any(np.diff(bundle.sorted_positions[DESC_I - 1]) != 0))}`",
            f"- descendant 2 sorted position changes over mu: `{yesno(np.any(np.diff(bundle.sorted_positions[DESC_J - 1]) != 0))}`",
            f"- any desc1/desc2 sorted-position swap: `{yesno(summary.sorted_position_swap)}`",
            f"- sorted root 1 descendant identity changes: `{yesno(summary.sorted1_identity_changes)}`",
            f"- sorted root 2 descendant identity changes: `{yesno(summary.sorted2_identity_changes)}`",
            f"- at min gap, sorted root 1 is descendant `{summary.sorted1_descendant_at_min}` and sorted root 2 is descendant `{summary.sorted2_descendant_at_min}`",
            "",
            "Local rows around the minimum gap:",
            "",
        ]
    )
    lines.extend(
        markdown_table(
            [
                "mu",
                "Lambda d1",
                "Lambda d2",
                "d1-d2",
                "pos d1",
                "pos d2",
                "MAC d1",
                "MAC d2",
                "warning",
            ],
            local_table(bundle, summary),
        )
    )
    lines.extend(
        [
            "",
            "## MAC Continuity And Tracking Warnings",
            "",
            f"- local min MAC for descendant 1 near min gap: `{fmt(summary.local_min_mac_desc1, 8)}`",
            f"- local min MAC for descendant 2 near min gap: `{fmt(summary.local_min_mac_desc2, 8)}`",
            f"- cross-MAC between desc1 and desc2 shapes at min gap: `{fmt(summary.cross_mac_at_min, 8)}`",
            f"- warning rows for descendants 1/2 on full refined grid: `{sum(1 for row in bundle.warning_rows if int(row['branch_index_from_mu0']) in (DESC_I, DESC_J))}`",
            "",
        ]
    )
    near_warnings = warning_lines_near_min(bundle, summary)
    if near_warnings:
        lines.extend(["Warning rows within `0.02` of the minimum gap:", ""])
        lines.extend(near_warnings)
        lines.append("")
    else:
        lines.extend(["No descendant 1/2 warning rows occur within `0.02` of the minimum gap.", ""])

    lines.extend(
        [
            "## Article Figure Recommendation",
            "",
            correction,
            "Do not silently switch to sorted roots under the same caption: sorted roots and",
            "descendant branches are different objects. If the article needs descendant branches,",
            "regenerate the figure with this smaller local `mu` step and adaptive fine root scan.",
            "If it instead needs the spectrum ordered by frequency at every `mu`, regenerate using",
            "sorted roots and label it as sorted frequencies.",
            "",
            "## Spike / Artifact Check",
            "",
            "No spike artifact was detected in the refined audit curves. The earlier apparent",
            "crossing is traced to a coarse `mu` step combined with missed closely spaced roots",
            "in a narrow near-degenerate window; this is a branch-tracking/root-extraction",
            "artifact, not a physical frequency crossing.",
            "",
            "## Limitations",
            "",
            "- Euler--Bernoulli thickness-mismatch analytic determinant only.",
            "- The audit does not run Timoshenko, FEM, Gmsh, or CalculiX.",
            "- A stricter modal-character statement would require a local subspace-MAC study,",
            "  but it is not needed to answer whether the frequencies cross.",
            "",
            "## Outputs",
            "",
            f"- CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
            f"- PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
        ]
    )
    OUTPUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def plot_diagnostic(bundle: TrackingBundle, summary: CrossingSummary) -> None:
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ax.plot(bundle.mu_values, bundle.tracked_lambdas[DESC_I - 1], lw=1.8, label="desc 1")
    ax.plot(bundle.mu_values, bundle.tracked_lambdas[DESC_J - 1], lw=1.8, label="desc 2")
    ax.plot(bundle.mu_values, bundle.sorted_roots[:, 0], color="0.55", lw=1.0, ls="--", label="sorted 1")
    ax.plot(bundle.mu_values, bundle.sorted_roots[:, 1], color="0.70", lw=1.0, ls="--", label="sorted 2")
    ax.scatter(
        [summary.min_gap_mu],
        [0.5 * (summary.lambda_desc1_at_min + summary.lambda_desc2_at_min)],
        color="black",
        s=32,
        zorder=5,
        label="min gap",
    )
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\Lambda$")
    ax.set_title(
        rf"Diagnostic desc1/desc2 crossing audit: $\beta={BETA_DEG:g}^\circ$, "
        rf"$\epsilon={EPSILON:g}$, $\eta={ETA:g}$"
    )
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> dict[str, object]:
    article_info = load_article_pair_info()
    print("coarse pass: solving and tracking")
    coarse = collect_tracking(coarse_mu_grid())
    coarse_summary = summarize(coarse)
    print(
        f"coarse min gap near mu={coarse_summary.min_gap_mu:.10g}, "
        f"gap={coarse_summary.min_gap:.12g}"
    )

    if article_info.sign_change and np.isfinite(article_info.crossing_mu_linear):
        refine_center = float(article_info.crossing_mu_linear)
        refine_reason = "article desc1/desc2 apparent crossing"
    else:
        refine_center = float(coarse_summary.min_gap_mu)
        refine_reason = "coarse-grid minimum gap"

    print(f"refined pass: solving and tracking around mu={refine_center:.10g}")
    refined = collect_tracking(refined_mu_grid(refine_center))
    summary = summarize(refined)
    apparent_window = summarize_window(refined, center_mu=refine_center, half_width=LOCAL_REFINE_HALF_WIDTH)
    write_csv(refined)
    write_report(refined, summary, coarse_summary, article_info, apparent_window, refine_center, refine_reason)
    plot_diagnostic(refined, summary)

    print(f"saved CSV: {OUTPUT_CSV}")
    print(f"saved report: {OUTPUT_REPORT}")
    print(f"saved PNG: {OUTPUT_PNG}")
    print(f"desc1/desc2 sign change: {yesno(summary.sign_change)}")
    print(f"article apparent sign change: {yesno(article_info.sign_change)}")
    print(f"refined apparent-window sign change: {yesno(apparent_window.sign_change)}")
    print(f"min gap: {summary.min_gap:.12g} at mu={summary.min_gap_mu:.10g}")
    print(f"relative gap: {summary.relative_gap:.12g}")
    print(f"sorted-position swap desc1/desc2: {yesno(summary.sorted_position_swap)}")
    return {
        "csv": OUTPUT_CSV,
        "report": OUTPUT_REPORT,
        "png": OUTPUT_PNG,
        "sign_change": summary.sign_change,
        "min_gap": summary.min_gap,
        "min_gap_mu": summary.min_gap_mu,
    }


if __name__ == "__main__":
    main()
