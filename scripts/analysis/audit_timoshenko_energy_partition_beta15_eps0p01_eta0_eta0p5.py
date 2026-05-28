from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path
import sys
from typing import Callable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment


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
from scripts.lib.variable_length_timoshenko import (  # noqa: E402
    E,
    omega_cutoff,
    project_omega,
    section_from_epsilon_tau,
    segment_lengths,
    tau_factors,
    timo_energy_partition,
    timo_mode_coefficients,
    timo_mode_fields,
    timo_sorted_roots,
)


BETA_DEG = 15.0
EPSILON = 0.01
ETA_VALUES = [0.0, 0.5]
MU_VALUES = np.round(np.linspace(0.0, 0.9, 91), 10)
N_BRANCHES = 6
N_SOLVE = 18
BETA_STEPS = 16
NUM_SHAPE_SAMPLES = 401
ENERGY_QUAD_POINTS = 801
MAC_WARNING_THRESHOLD = 0.8
FREQ_WEIGHT = 0.03
THIN_ROD_LIMIT = 0.1

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_STEM = "timoshenko_energy_partition_beta15_eps0p01_eta0_eta0p5"
OUTPUT_CSV = OUTPUT_DIR / f"{OUTPUT_STEM}.csv"
OUTPUT_REPORT = OUTPUT_DIR / f"{OUTPUT_STEM}_report.md"
OUTPUT_SHEAR_PNG = OUTPUT_DIR / f"{OUTPUT_STEM}_shear_fraction.png"
OUTPUT_ROD_PNG = OUTPUT_DIR / f"{OUTPUT_STEM}_rod_energy_fraction.png"


@dataclass(frozen=True)
class ModeState:
    sorted_index: int
    Lambda: float
    vector: np.ndarray
    coeff: np.ndarray
    smallest_singular_value: float
    singular_value_ratio: float
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class TrackedPoint:
    model: str
    eta: float
    branch_index: int
    branch_id: str
    base_sorted_index: int
    beta_deg: float
    mu: float
    current_sorted_index: int
    Lambda: float
    mac_to_previous: float
    coeff: np.ndarray
    vector: np.ndarray
    smallest_singular_value: float
    singular_value_ratio: float


@dataclass(frozen=True)
class TrackingResult:
    model: str
    eta: float
    points: dict[str, dict[float, TrackedPoint]]
    warnings: tuple[str, ...]


def branch_id_from_index(index: int) -> str:
    return f"bending_desc_{int(index):02d}"


def mu_key(value: float) -> float:
    return round(float(value), 10)


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if abs(float(denominator)) > 1.0e-30 else float("nan")


def format_float(value: object, digits: int = 6) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(number):
        return "nan"
    return f"{number:.{digits}g}"


def null_vector(matrix: np.ndarray) -> tuple[np.ndarray, float, float]:
    _u, singular_values, vh = np.linalg.svd(np.asarray(matrix, dtype=float))
    coeff = vh[-1, :].astype(float)
    if coeff.size:
        pivot = int(np.argmax(np.abs(coeff)))
        if coeff[pivot] < 0.0:
            coeff = -coeff
    smallest = float(singular_values[-1])
    ratio = (
        float(singular_values[-1] / singular_values[-2])
        if len(singular_values) >= 2 and abs(float(singular_values[-2])) > 1.0e-14
        else float("nan")
    )
    return coeff, smallest, ratio


def normalize_vector(vector: np.ndarray) -> np.ndarray:
    out = np.asarray(vector, dtype=float)
    scale = float(np.linalg.norm(out))
    if scale <= 1.0e-28 or not math.isfinite(scale):
        return np.full_like(out, np.nan, dtype=float)
    return out / scale


def global_centerline_vector(
    beta_deg: float,
    u1: np.ndarray,
    w1: np.ndarray,
    u2: np.ndarray,
    w2: np.ndarray,
) -> np.ndarray:
    beta = math.radians(float(beta_deg))
    cb = math.cos(beta)
    sb = math.sin(beta)
    parts: list[np.ndarray] = []
    ux1 = np.asarray(u1, dtype=float)
    uy1 = np.asarray(w1, dtype=float)
    rod1 = np.empty(2 * len(ux1), dtype=float)
    rod1[0::2] = ux1
    rod1[1::2] = uy1
    parts.append(rod1)

    u2_arr = np.asarray(u2, dtype=float)
    w2_arr = np.asarray(w2, dtype=float)
    ux2 = u2_arr * cb - w2_arr * sb
    uy2 = u2_arr * sb + w2_arr * cb
    rod2 = np.empty(2 * len(ux2), dtype=float)
    rod2[0::2] = ux2
    rod2[1::2] = uy2
    parts.append(rod2)
    return normalize_vector(np.concatenate(parts))


def centerline_mac(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    denominator = float(np.dot(a, a) * np.dot(b, b))
    if denominator <= 1.0e-28 or not math.isfinite(denominator):
        return float("nan")
    dot_value = float(np.dot(a, b))
    return float((dot_value * dot_value) / denominator)


def eb_centerline_vector(
    Lambda: float,
    beta_deg: float,
    mu: float,
    epsilon: float,
    eta: float,
    coeff: np.ndarray,
) -> np.ndarray:
    factors = tau_factors(float(mu), float(eta))
    A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]
    xi = np.linspace(0.0, 1.0, NUM_SHAPE_SAMPLES, dtype=float)
    l1, l2 = segment_lengths(factors.mu)
    z1 = float(Lambda) * l1 * xi / math.sqrt(factors.tau1)
    theta1 = float(epsilon) * float(Lambda) ** 2 * l1 * xi
    z2 = -float(Lambda) * l2 * xi / math.sqrt(factors.tau2)
    theta2 = -float(epsilon) * float(Lambda) ** 2 * l2 * xi
    w1 = A1 * (np.cos(z1) - np.cosh(z1)) + B1 * (np.sin(z1) - np.sinh(z1))
    u1 = P1 * np.sin(theta1)
    w2 = A2 * (np.cos(z2) - np.cosh(z2)) + B2 * (np.sin(z2) - np.sinh(z2))
    u2 = P2 * np.sin(theta2)
    return global_centerline_vector(float(beta_deg), u1, w1, u2, w2)


def timo_centerline_vector(
    Lambda: float,
    beta_deg: float,
    mu: float,
    epsilon: float,
    eta: float,
    coeff: np.ndarray,
) -> tuple[np.ndarray, tuple[str, ...]]:
    fields = timo_mode_fields(
        float(Lambda),
        float(beta_deg),
        float(mu),
        float(epsilon),
        float(eta),
        coeff=np.asarray(coeff, dtype=float),
        n_points=NUM_SHAPE_SAMPLES,
    )
    rod1 = fields["rod1"]
    rod2 = fields["rod2"]
    vector = global_centerline_vector(
        float(beta_deg),
        np.asarray(rod1["u"], dtype=float),  # type: ignore[index]
        np.asarray(rod1["w"], dtype=float),  # type: ignore[index]
        np.asarray(rod2["u"], dtype=float),  # type: ignore[index]
        np.asarray(rod2["w"], dtype=float),  # type: ignore[index]
    )
    return vector, tuple(fields["warnings"])  # type: ignore[arg-type]


def eb_sorted_roots(beta_deg: float, mu: float, epsilon: float, eta: float, n_roots: int) -> np.ndarray:
    return find_first_n_roots_eta(
        beta=math.radians(float(beta_deg)),
        mu=float(mu),
        epsilon=float(epsilon),
        eta=float(eta),
        n_roots=int(n_roots),
        Lmin=0.2,
        Lmax0=45.0,
        scan_step=0.01,
        grow_factor=1.35,
        max_tries=8,
    )


def solve_eb_states(beta_deg: float, mu: float, eta: float, n_roots: int) -> tuple[list[ModeState], list[str]]:
    roots = eb_sorted_roots(float(beta_deg), float(mu), EPSILON, float(eta), int(n_roots))
    states: list[ModeState] = []
    warnings: list[str] = []
    for sorted_index, root in enumerate(roots, start=1):
        if not math.isfinite(float(root)):
            warnings.append(
                f"EB eta={eta:g}, beta={beta_deg:g}, mu={mu:g}, sorted={sorted_index}: missing root"
            )
            continue
        matrix = assemble_clamped_coupled_matrix_eta(
            float(root),
            math.radians(float(beta_deg)),
            float(mu),
            EPSILON,
            float(eta),
        )
        coeff, smallest, ratio = null_vector(matrix)
        vector = eb_centerline_vector(float(root), float(beta_deg), float(mu), EPSILON, float(eta), coeff)
        states.append(
            ModeState(
                sorted_index=int(sorted_index),
                Lambda=float(root),
                vector=vector,
                coeff=coeff,
                smallest_singular_value=smallest,
                singular_value_ratio=ratio,
                warnings=(),
            )
        )
    return states, warnings


def solve_timo_states(beta_deg: float, mu: float, eta: float, n_roots: int) -> tuple[list[ModeState], list[str]]:
    roots, root_warnings = timo_sorted_roots(float(beta_deg), float(mu), EPSILON, int(n_roots), eta=float(eta))
    states: list[ModeState] = []
    warnings = [f"Timoshenko eta={eta:g}, beta={beta_deg:g}, mu={mu:g}: {item}" for item in root_warnings]
    for sorted_index, root in enumerate(roots, start=1):
        if not math.isfinite(float(root)):
            warnings.append(
                f"Timoshenko eta={eta:g}, beta={beta_deg:g}, mu={mu:g}, sorted={sorted_index}: missing root"
            )
            continue
        try:
            mode = timo_mode_coefficients(float(root), float(beta_deg), float(mu), EPSILON, float(eta))
            vector, shape_warnings = timo_centerline_vector(
                float(root),
                float(beta_deg),
                float(mu),
                EPSILON,
                float(eta),
                mode.coeff,
            )
        except (FloatingPointError, ValueError, OverflowError, np.linalg.LinAlgError) as exc:
            warnings.append(
                f"Timoshenko eta={eta:g}, beta={beta_deg:g}, mu={mu:g}, sorted={sorted_index}: shape failed: {exc}"
            )
            continue
        warnings.extend(
            f"Timoshenko eta={eta:g}, beta={beta_deg:g}, mu={mu:g}, sorted={sorted_index}: {item}"
            for item in (*mode.warnings, *shape_warnings)
        )
        states.append(
            ModeState(
                sorted_index=int(sorted_index),
                Lambda=float(root),
                vector=vector,
                coeff=mode.coeff,
                smallest_singular_value=mode.smallest_singular_value,
                singular_value_ratio=mode.singular_value_ratio,
                warnings=tuple(mode.warnings),
            )
        )
    return states, warnings


def assign_previous_to_current(
    branch_ids: Sequence[str],
    previous_points: dict[str, TrackedPoint],
    current_states: Sequence[ModeState],
) -> dict[str, tuple[ModeState, float]]:
    if not branch_ids or not current_states:
        return {}
    cost = np.ones((len(branch_ids), len(current_states)), dtype=float)
    mac = np.zeros_like(cost)
    for row, branch_id in enumerate(branch_ids):
        previous = previous_points[branch_id]
        previous_lambda = max(abs(float(previous.Lambda)), 1.0e-14)
        for col, state in enumerate(current_states):
            mac_value = centerline_mac(previous.vector, state.vector)
            mac[row, col] = mac_value if math.isfinite(mac_value) else 0.0
            relative_frequency_step = abs(float(state.Lambda) - float(previous.Lambda)) / previous_lambda
            cost[row, col] = 1.0 - mac[row, col] + FREQ_WEIGHT * relative_frequency_step
    rows, cols = linear_sum_assignment(cost)
    assigned: dict[str, tuple[ModeState, float]] = {}
    for row, col in zip(rows, cols, strict=True):
        assigned[branch_ids[int(row)]] = (current_states[int(col)], float(mac[int(row), int(col)]))
    return assigned


def tracking_path(mu_values: Sequence[float]) -> list[tuple[float, float]]:
    path: list[tuple[float, float]] = [(float(beta), 0.0) for beta in np.linspace(0.0, BETA_DEG, BETA_STEPS)]
    for mu in mu_values:
        if abs(float(mu)) <= 1.0e-12:
            continue
        path.append((BETA_DEG, float(mu)))
    return path


def track_descendants(
    *,
    model: str,
    eta: float,
    solver: Callable[[float, float, float, int], tuple[list[ModeState], list[str]]],
) -> TrackingResult:
    warnings: list[str] = []
    branch_ids = [branch_id_from_index(index) for index in range(1, N_BRANCHES + 1)]
    first_states, first_warnings = solver(0.0, 0.0, float(eta), N_SOLVE)
    warnings.extend(first_warnings)
    if len(first_states) < N_BRANCHES:
        raise RuntimeError(f"{model} eta={eta:g}: only {len(first_states)} roots found at beta=0, mu=0.")

    previous: dict[str, TrackedPoint] = {}
    tracked: dict[str, dict[float, TrackedPoint]] = {branch_id: {} for branch_id in branch_ids}
    for branch_index, (branch_id, state) in enumerate(zip(branch_ids, first_states[:N_BRANCHES], strict=True), start=1):
        previous[branch_id] = TrackedPoint(
            model=model,
            eta=float(eta),
            branch_index=int(branch_index),
            branch_id=branch_id,
            base_sorted_index=int(branch_index),
            beta_deg=0.0,
            mu=0.0,
            current_sorted_index=int(state.sorted_index),
            Lambda=float(state.Lambda),
            mac_to_previous=1.0,
            coeff=state.coeff,
            vector=state.vector,
            smallest_singular_value=state.smallest_singular_value,
            singular_value_ratio=state.singular_value_ratio,
        )

    target_mu_keys = {mu_key(value) for value in MU_VALUES}
    for beta_deg, mu in tracking_path(MU_VALUES)[1:]:
        states, state_warnings = solver(float(beta_deg), float(mu), float(eta), N_SOLVE)
        warnings.extend(state_warnings)
        assignment = assign_previous_to_current(branch_ids, previous, states)
        if len(assignment) < len(branch_ids):
            warnings.append(f"{model} eta={eta:g}: incomplete assignment at beta={beta_deg:g}, mu={mu:g}.")
        for branch_index, branch_id in enumerate(branch_ids, start=1):
            if branch_id not in assignment:
                continue
            state, mac_value = assignment[branch_id]
            point = TrackedPoint(
                model=model,
                eta=float(eta),
                branch_index=int(branch_index),
                branch_id=branch_id,
                base_sorted_index=int(branch_index),
                beta_deg=float(beta_deg),
                mu=float(mu),
                current_sorted_index=int(state.sorted_index),
                Lambda=float(state.Lambda),
                mac_to_previous=float(mac_value),
                coeff=state.coeff,
                vector=state.vector,
                smallest_singular_value=state.smallest_singular_value,
                singular_value_ratio=state.singular_value_ratio,
            )
            previous[branch_id] = point
            if abs(float(beta_deg) - BETA_DEG) <= 1.0e-10 and mu_key(mu) in target_mu_keys:
                tracked[branch_id][mu_key(mu)] = point
                if math.isfinite(mac_value) and mac_value < MAC_WARNING_THRESHOLD:
                    warnings.append(
                        f"{model} eta={eta:g}, branch={branch_id}, mu={mu:g}: low MAC {mac_value:.6g}"
                    )
    return TrackingResult(model=model, eta=float(eta), points=tracked, warnings=tuple(sorted(set(warnings))))


def applicability_terms(mu: float, eta: float, Lambda: float) -> dict[str, object]:
    factors = tau_factors(float(mu), float(eta))
    l1, l2 = segment_lengths(factors.mu)
    ratio1 = 4.0 * EPSILON * factors.tau1 / l1
    ratio2 = 4.0 * EPSILON * factors.tau2 / l2
    section1 = section_from_epsilon_tau(EPSILON, factors.tau1)
    section2 = section_from_epsilon_tau(EPSILON, factors.tau2)
    omega = project_omega(float(Lambda), EPSILON)
    cutoff_ratio1 = omega / omega_cutoff(section1)
    cutoff_ratio2 = omega / omega_cutoff(section2)
    max_cutoff = max(cutoff_ratio1, cutoff_ratio2)
    return {
        "tau1": factors.tau1,
        "tau2": factors.tau2,
        "l1": l1,
        "l2": l2,
        "diameter_to_length_rod1": ratio1,
        "diameter_to_length_rod2": ratio2,
        "rod1_EB_thin_valid": bool(ratio1 <= THIN_ROD_LIMIT + 1.0e-15),
        "rod2_EB_thin_valid": bool(ratio2 <= THIN_ROD_LIMIT + 1.0e-15),
        "EB_valid_all_rods": bool(max(ratio1, ratio2) <= THIN_ROD_LIMIT + 1.0e-15),
        "max_diameter_to_length": max(ratio1, ratio2),
        "Omega": omega,
        "Omega_c_rod1": omega_cutoff(section1),
        "Omega_c_rod2": omega_cutoff(section2),
        "Omega_over_cutoff_rod1": cutoff_ratio1,
        "Omega_over_cutoff_rod2": cutoff_ratio2,
        "max_Omega_over_cutoff": max_cutoff,
        "cutoff_warning": bool(max_cutoff >= 0.8),
        "cutoff_violation": bool(max_cutoff >= 1.0),
    }


def build_rows(
    timo_tracking: dict[float, TrackingResult],
    eb_tracking: dict[float, TrackingResult],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for eta in ETA_VALUES:
        for branch_index in range(1, N_BRANCHES + 1):
            branch_id = branch_id_from_index(branch_index)
            for mu in MU_VALUES:
                key = mu_key(float(mu))
                timo_point = timo_tracking[eta].points.get(branch_id, {}).get(key)
                eb_point = eb_tracking[eta].points.get(branch_id, {}).get(key)
                if timo_point is None:
                    continue
                energy = timo_energy_partition(
                    timo_point.Lambda,
                    BETA_DEG,
                    float(mu),
                    EPSILON,
                    float(eta),
                    coeff=timo_point.coeff,
                    n_points=ENERGY_QUAD_POINTS,
                )
                applicability = applicability_terms(float(mu), float(eta), timo_point.Lambda)
                lambda_eb = eb_point.Lambda if eb_point is not None else float("nan")
                signed_diff = timo_point.Lambda - lambda_eb if math.isfinite(float(lambda_eb)) else float("nan")
                abs_diff = abs(signed_diff) if math.isfinite(float(signed_diff)) else float("nan")
                rel_diff = signed_diff / lambda_eb if math.isfinite(float(lambda_eb)) and abs(lambda_eb) > 1.0e-14 else float("nan")
                abs_rel_diff = abs_diff / abs(lambda_eb) if math.isfinite(float(lambda_eb)) and abs(lambda_eb) > 1.0e-14 else float("nan")
                thick_energy = energy["rod2_energy_fraction"] if eta > 0.0 else float("nan")
                thick_shear = energy["rod2_shear_fraction"] if eta > 0.0 else float("nan")
                thin_energy = energy["rod1_energy_fraction"] if eta > 0.0 else float("nan")
                rows.append(
                    {
                        "eta": float(eta),
                        "mu": float(mu),
                        "branch_index": int(branch_index),
                        "branch_id": branch_id,
                        "Lambda_Timoshenko": float(timo_point.Lambda),
                        "Lambda_EB": lambda_eb,
                        "signed_diff": signed_diff,
                        "abs_diff": abs_diff,
                        "rel_diff": rel_diff,
                        "abs_rel_diff": abs_rel_diff,
                        "timo_current_sorted_index": int(timo_point.current_sorted_index),
                        "eb_current_sorted_index": eb_point.current_sorted_index if eb_point is not None else "",
                        "timo_mac_to_previous": float(timo_point.mac_to_previous),
                        "eb_mac_to_previous": eb_point.mac_to_previous if eb_point is not None else "",
                        "timo_svd_smallest": float(timo_point.smallest_singular_value),
                        "timo_svd_ratio": float(timo_point.singular_value_ratio),
                        "eb_svd_smallest": eb_point.smallest_singular_value if eb_point is not None else "",
                        "eb_svd_ratio": eb_point.singular_value_ratio if eb_point is not None else "",
                        "thick_rod": "rod2" if eta > 0.0 else "",
                        "thick_rod_energy_fraction": thick_energy,
                        "thick_rod_shear_fraction": thick_shear,
                        "thin_rod_energy_fraction": thin_energy,
                        **energy,
                        **applicability,
                    }
                )
    return rows


def write_csv(rows: Sequence[dict[str, object]]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def finite_row_value(row: dict[str, object], key: str) -> float:
    try:
        value = float(row.get(key, float("nan")))
    except (TypeError, ValueError):
        return float("nan")
    return value if math.isfinite(value) else float("nan")


def rows_for(rows: Sequence[dict[str, object]], *, eta: float | None = None, branch_index: int | None = None) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in rows:
        if eta is not None and abs(finite_row_value(row, "eta") - float(eta)) > 1.0e-12:
            continue
        if branch_index is not None and int(row["branch_index"]) != int(branch_index):
            continue
        out.append(row)
    return out


def mean_finite(values: Sequence[float]) -> float:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    return float(np.mean(finite)) if finite else float("nan")


def max_finite(values: Sequence[float]) -> float:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    return max(finite) if finite else float("nan")


def pearson(x_values: Sequence[float], y_values: Sequence[float]) -> float:
    pairs = [
        (float(x), float(y))
        for x, y in zip(x_values, y_values, strict=True)
        if math.isfinite(float(x)) and math.isfinite(float(y))
    ]
    if len(pairs) < 3:
        return float("nan")
    x = np.array([item[0] for item in pairs], dtype=float)
    y = np.array([item[1] for item in pairs], dtype=float)
    if float(np.std(x)) <= 1.0e-30 or float(np.std(y)) <= 1.0e-30:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def summarize_h1(rows: Sequence[dict[str, object]]) -> tuple[list[dict[str, object]], str, list[int]]:
    summary: list[dict[str, object]] = []
    failing: list[int] = []
    for branch_index in range(1, N_BRANCHES + 1):
        subset = [
            row
            for row in rows_for(rows, eta=0.0, branch_index=branch_index)
            if finite_row_value(row, "mu") > 0.6 + 1.0e-12
        ]
        rod1_values = [finite_row_value(row, "rod1_energy_fraction") for row in subset]
        shear_values = [finite_row_value(row, "shear_fraction") for row in subset]
        mean_rod1 = mean_finite(rod1_values)
        max_rod1 = max_finite(rod1_values)
        mean_shear = mean_finite(shear_values)
        if math.isfinite(max_rod1) and max_rod1 >= 0.5:
            failing.append(branch_index)
        summary.append(
            {
                "branch_index": branch_index,
                "mean_short_rod1_energy_fraction_mu_gt_0p6": mean_rod1,
                "max_short_rod1_energy_fraction_mu_gt_0p6": max_rod1,
                "mean_shear_fraction_mu_gt_0p6": mean_shear,
            }
        )
    if failing:
        status = "does not support"
    elif any(math.isfinite(float(item["max_short_rod1_energy_fraction_mu_gt_0p6"])) for item in summary):
        status = "supports"
    else:
        status = "inconclusive"
    return summary, status, failing


def summarize_h2(rows: Sequence[dict[str, object]]) -> tuple[list[dict[str, object]], str]:
    summary: list[dict[str, object]] = []
    positive_energy_deltas = 0
    finite_deltas = 0
    for branch_index in range(1, N_BRANCHES + 1):
        eta05 = rows_for(rows, eta=0.5, branch_index=branch_index)
        eta0 = rows_for(rows, eta=0.0, branch_index=branch_index)
        by_mu_eta0 = {mu_key(finite_row_value(row, "mu")): row for row in eta0}
        energy_deltas: list[float] = []
        shear_deltas: list[float] = []
        for row in eta05:
            mu = mu_key(finite_row_value(row, "mu"))
            base = by_mu_eta0.get(mu)
            if base is None:
                continue
            delta_energy = finite_row_value(row, "thick_rod_energy_fraction") - finite_row_value(
                base,
                "rod2_energy_fraction",
            )
            delta_shear = finite_row_value(row, "thick_rod_shear_fraction") - finite_row_value(
                base,
                "rod2_shear_fraction",
            )
            if math.isfinite(delta_energy):
                finite_deltas += 1
                if delta_energy > 0.0:
                    positive_energy_deltas += 1
            energy_deltas.append(delta_energy)
            shear_deltas.append(delta_shear)
        summary.append(
            {
                "branch_index": branch_index,
                "mean_thick_rod_energy_fraction_eta0p5": mean_finite(
                    [finite_row_value(row, "thick_rod_energy_fraction") for row in eta05]
                ),
                "max_thick_rod_energy_fraction_eta0p5": max_finite(
                    [finite_row_value(row, "thick_rod_energy_fraction") for row in eta05]
                ),
                "mean_thick_rod_shear_fraction_eta0p5": mean_finite(
                    [finite_row_value(row, "thick_rod_shear_fraction") for row in eta05]
                ),
                "max_thick_rod_shear_fraction_eta0p5": max_finite(
                    [finite_row_value(row, "thick_rod_shear_fraction") for row in eta05]
                ),
                "mean_delta_vs_eta0_rod2_energy_fraction": mean_finite(energy_deltas),
                "mean_delta_vs_eta0_rod2_shear_fraction": mean_finite(shear_deltas),
            }
        )
    if finite_deltas == 0:
        status = "inconclusive"
    elif positive_energy_deltas > 0.6 * finite_deltas:
        status = "supports"
    elif positive_energy_deltas < 0.4 * finite_deltas:
        status = "does not support"
    else:
        status = "inconclusive"
    return summary, status


def top_rows(rows: Sequence[dict[str, object]], key: str, count: int = 10) -> list[dict[str, object]]:
    return sorted(
        [row for row in rows if math.isfinite(finite_row_value(row, key))],
        key=lambda row: finite_row_value(row, key),
        reverse=True,
    )[: int(count)]


def row_identity(row: dict[str, object]) -> tuple[float, int, float]:
    return (finite_row_value(row, "eta"), int(row["branch_index"]), mu_key(finite_row_value(row, "mu")))


def table_from_rows(rows: Sequence[dict[str, object]], columns: Sequence[tuple[str, str]]) -> list[str]:
    header = "| " + " | ".join(label for label, _key in columns) + " |"
    divider = "| " + " | ".join("---" for _ in columns) + " |"
    lines = [header, divider]
    for row in rows:
        cells: list[str] = []
        for _label, key in columns:
            value = row.get(key, "")
            if isinstance(value, float):
                cells.append(format_float(value, 5))
            else:
                try:
                    cells.append(format_float(value, 5))
                except Exception:
                    cells.append(str(value))
        lines.append("| " + " | ".join(cells) + " |")
    return lines


def h_summary_table(summary: Sequence[dict[str, object]], columns: Sequence[tuple[str, str]]) -> list[str]:
    return table_from_rows(summary, columns)


def applicability_summary(rows: Sequence[dict[str, object]]) -> dict[str, object]:
    eb_invalid = [row for row in rows if not bool(row.get("EB_valid_all_rods"))]
    cutoff_warnings = [row for row in rows if bool(row.get("cutoff_warning"))]
    cutoff_violations = [row for row in rows if bool(row.get("cutoff_violation"))]
    max_diameter_row = max(rows, key=lambda row: finite_row_value(row, "max_diameter_to_length"))
    max_cutoff_row = max(rows, key=lambda row: finite_row_value(row, "max_Omega_over_cutoff"))
    first_invalid_by_eta: dict[float, float] = {}
    for eta in ETA_VALUES:
        invalid_mus = [
            finite_row_value(row, "mu")
            for row in rows_for(rows, eta=eta)
            if not bool(row.get("EB_valid_all_rods"))
        ]
        first_invalid_by_eta[eta] = min(invalid_mus) if invalid_mus else float("nan")
    return {
        "eb_invalid_count": len(eb_invalid),
        "cutoff_warning_count": len(cutoff_warnings),
        "cutoff_violation_count": len(cutoff_violations),
        "max_diameter_row": max_diameter_row,
        "max_cutoff_row": max_cutoff_row,
        "first_invalid_by_eta": first_invalid_by_eta,
    }


def plot_shear_fraction(rows: Sequence[dict[str, object]]) -> None:
    colors = plt.cm.viridis(np.linspace(0.05, 0.9, N_BRANCHES))
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), sharey=True, constrained_layout=True)
    for ax, eta in zip(axes, ETA_VALUES, strict=True):
        for branch_index in range(1, N_BRANCHES + 1):
            subset = rows_for(rows, eta=eta, branch_index=branch_index)
            ax.plot(
                [finite_row_value(row, "mu") for row in subset],
                [finite_row_value(row, "shear_fraction") for row in subset],
                color=colors[branch_index - 1],
                linewidth=1.5,
                label=f"branch {branch_index}",
            )
        ax.set_title(f"eta={eta:g}")
        ax.set_xlabel("mu")
        ax.grid(True, color="0.88", linewidth=0.8)
    axes[0].set_ylabel("shear energy fraction")
    axes[1].legend(loc="upper left", fontsize=8)
    fig.savefig(OUTPUT_SHEAR_PNG, dpi=220)
    plt.close(fig)


def plot_rod_energy_fraction(rows: Sequence[dict[str, object]]) -> None:
    colors = plt.cm.plasma(np.linspace(0.05, 0.9, N_BRANCHES))
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), sharey=True, constrained_layout=True)
    panels = [
        (0.0, "rod1_energy_fraction", "eta=0: short rod1 fraction"),
        (0.5, "thick_rod_energy_fraction", "eta=0.5: thick rod2 fraction"),
    ]
    for ax, (eta, key, title) in zip(axes, panels, strict=True):
        for branch_index in range(1, N_BRANCHES + 1):
            subset = rows_for(rows, eta=eta, branch_index=branch_index)
            ax.plot(
                [finite_row_value(row, "mu") for row in subset],
                [finite_row_value(row, key) for row in subset],
                color=colors[branch_index - 1],
                linewidth=1.5,
                label=f"branch {branch_index}",
            )
        ax.set_title(title)
        ax.set_xlabel("mu")
        ax.grid(True, color="0.88", linewidth=0.8)
    axes[0].set_ylabel("potential energy fraction")
    axes[1].legend(loc="upper left", fontsize=8)
    fig.savefig(OUTPUT_ROD_PNG, dpi=220)
    plt.close(fig)


def write_report(rows: Sequence[dict[str, object]], tracking_warnings: Sequence[str]) -> dict[str, object]:
    h1_summary, h1_status, h1_failures = summarize_h1(rows)
    h2_summary, h2_status = summarize_h2(rows)
    app = applicability_summary(rows)

    finite_diff_rows = [row for row in rows if math.isfinite(finite_row_value(row, "abs_diff"))]
    shear_corr_abs = pearson(
        [finite_row_value(row, "shear_fraction") for row in finite_diff_rows],
        [finite_row_value(row, "abs_diff") for row in finite_diff_rows],
    )
    shear_corr_rel = pearson(
        [finite_row_value(row, "shear_fraction") for row in finite_diff_rows],
        [finite_row_value(row, "abs_rel_diff") for row in finite_diff_rows],
    )
    eta05_rows = [row for row in finite_diff_rows if abs(finite_row_value(row, "eta") - 0.5) <= 1.0e-12]
    thick_shear_corr_abs = pearson(
        [finite_row_value(row, "thick_rod_shear_fraction") for row in eta05_rows],
        [finite_row_value(row, "abs_diff") for row in eta05_rows],
    )
    thick_shear_corr_rel = pearson(
        [finite_row_value(row, "thick_rod_shear_fraction") for row in eta05_rows],
        [finite_row_value(row, "abs_rel_diff") for row in eta05_rows],
    )

    top_shear = top_rows(rows, "shear_fraction", 10)
    top_abs = top_rows(rows, "abs_diff", 10)
    top_rel = top_rows(rows, "abs_rel_diff", 10)
    shear_ids = {row_identity(row) for row in top_shear}
    overlap_abs = len(shear_ids & {row_identity(row) for row in top_abs})
    overlap_rel = len(shear_ids & {row_identity(row) for row in top_rel})
    if (math.isfinite(shear_corr_abs) and shear_corr_abs > 0.5) or overlap_abs >= 3 or overlap_rel >= 3:
        h3_status = "supports"
    elif math.isfinite(shear_corr_abs) and shear_corr_abs < 0.1 and overlap_abs == 0 and overlap_rel == 0:
        h3_status = "does not support"
    else:
        h3_status = "inconclusive"

    max_diam = app["max_diameter_row"]
    max_cutoff = app["max_cutoff_row"]
    warnings_unique = sorted(set(tracking_warnings))
    lines: list[str] = [
        "# Timoshenko Energy Partition Audit",
        "",
        "## Purpose",
        "",
        "Diagnostic-only energy partition audit for tau-aware variable-length Timoshenko coupled rods. "
        "The goal is to check whether low descendant branches localize mostly on the long/thin arm and "
        "therefore explain why Euler--Bernoulli and Timoshenko frequencies can remain close even when "
        "the Euler--Bernoulli diameter criterion fails locally.",
        "",
        "## Parameters",
        "",
        "| parameter | value |",
        "| --- | ---: |",
        f"| beta | {BETA_DEG:g} deg |",
        f"| epsilon | {EPSILON:g} |",
        f"| eta values | {', '.join(format_float(value, 3) for value in ETA_VALUES)} |",
        f"| mu grid | {len(MU_VALUES)} points from {MU_VALUES[0]:g} to {MU_VALUES[-1]:g} |",
        f"| descendant branches | {N_BRANCHES} |",
        f"| sorted roots solved per step | {N_SOLVE} |",
        f"| energy quadrature points per rod | {ENERGY_QUAD_POINTS} |",
        "",
        "## Verification Reminder",
        "",
        "The variable-length/tau-aware Timoshenko helper is diagnostic-verified for sorted-root limit checks "
        "in the existing eta=0 and eta!=0 audits. This energy audit adds descendant energy diagnostics and "
        "frequency comparisons; it is not a symbolic proof and does not use FEM.",
        "",
        "## Energy Formulas",
        "",
        "```text",
        "U_b,i = 1/2 int E*I_i*(psi_i')^2 dx",
        "U_s,i = 1/2 int kappa*G*A_i*(w_i' - psi_i)^2 dx",
        "U_a,i = 1/2 int E*A_i*(u_i')^2 dx",
        "T_trans,i proportional to 1/2 int rho*A_i*w_i^2 dx",
        "T_axial,i proportional to 1/2 int rho*A_i*u_i^2 dx",
        "T_rot,i proportional to 1/2 int rho*I_i*psi_i^2 dx",
        "```",
        "",
        "The quadrature uses the same signed local arm coordinates as the Timoshenko determinant "
        "(`x in [0,l1]` for rod 1 and `x in [0,-l2]` for rod 2) and integrates over `|dx|`.",
        "",
        "## Branch Tracking Convention",
        "",
        "For each eta separately, both Timoshenko and EB descendants are seeded from sorted roots at "
        "`beta=0, mu=0` for that eta. They are then continued by centerline shape-MAC assignment from "
        "`beta=0` to `beta=15 deg` at `mu=0`, followed by the `mu=0..0.9` sweep. "
        "The branch id is the base sorted index; current sorted index is reported only as metadata.",
        "",
        "EB frequencies use the existing diagnostic thickness-mismatch Euler--Bernoulli determinant with "
        "the same beta-then-mu tracking convention. The CSV reports both signed and absolute "
        "Timoshenko-minus-EB differences.",
        "",
        "## Applicability Summary",
        "",
        f"- EB diameter criterion: `2*r_i/l_i <= {THIN_ROD_LIMIT:g}`.",
        f"- rows with at least one EB diameter violation: {app['eb_invalid_count']} / {len(rows)}.",
        "- first violating mu by eta: "
        + ", ".join(
            f"eta={eta:g}: {format_float(app['first_invalid_by_eta'][eta], 4)}"
            for eta in ETA_VALUES
        ),
        "- max diameter-to-length ratio: "
        f"{format_float(finite_row_value(max_diam, 'max_diameter_to_length'), 6)} "
        f"(eta={format_float(max_diam['eta'], 3)}, branch={max_diam['branch_index']}, "
        f"mu={format_float(max_diam['mu'], 3)})",
        f"- Timoshenko cutoff warnings (Omega/Omega_c >= 0.8): {app['cutoff_warning_count']}.",
        f"- Timoshenko cutoff violations (Omega/Omega_c >= 1.0): {app['cutoff_violation_count']}.",
        "- max Omega/Omega_c: "
        f"{format_float(finite_row_value(max_cutoff, 'max_Omega_over_cutoff'), 6)} "
        f"(eta={format_float(max_cutoff['eta'], 3)}, branch={max_cutoff['branch_index']}, "
        f"mu={format_float(max_cutoff['mu'], 3)})",
        "",
        "## H1: eta=0 Short-Rod Energy For mu>0.6",
        "",
        f"Status: {h1_status}. For eta=0 and mu>0, rod 1 is the short rod and rod 2 is the long rod.",
        "",
        *h_summary_table(
            h1_summary,
            [
                ("branch", "branch_index"),
                ("mean rod1 energy", "mean_short_rod1_energy_fraction_mu_gt_0p6"),
                ("max rod1 energy", "max_short_rod1_energy_fraction_mu_gt_0p6"),
                ("mean shear", "mean_shear_fraction_mu_gt_0p6"),
            ],
        ),
        "",
        "A branch is counted as failing the 'most energy in the long rod' form of H1 when "
        "`max rod1 energy fraction >= 0.5` on `mu>0.6`.",
        f"Failing branches: {', '.join(str(item) for item in h1_failures) if h1_failures else 'none'}.",
        "",
        "## H2: eta=0.5 Thick-Rod Involvement",
        "",
        f"Status: {h2_status}. For eta=0.5, rod 2 is the thick rod. The delta columns compare eta=0.5 "
        "rod-2 fractions against eta=0 rod-2 fractions at matching branch and mu.",
        "",
        *h_summary_table(
            h2_summary,
            [
                ("branch", "branch_index"),
                ("mean thick energy", "mean_thick_rod_energy_fraction_eta0p5"),
                ("max thick energy", "max_thick_rod_energy_fraction_eta0p5"),
                ("mean thick shear", "mean_thick_rod_shear_fraction_eta0p5"),
                ("max thick shear", "max_thick_rod_shear_fraction_eta0p5"),
                ("mean energy delta", "mean_delta_vs_eta0_rod2_energy_fraction"),
                ("mean shear delta", "mean_delta_vs_eta0_rod2_shear_fraction"),
            ],
        ),
        "",
        "## H3: Frequency Differences Versus Shear Energy",
        "",
        f"Status: {h3_status}.",
        f"- Pearson corr(shear_fraction, abs_diff): {format_float(shear_corr_abs, 5)}",
        f"- Pearson corr(shear_fraction, abs_rel_diff): {format_float(shear_corr_rel, 5)}",
        f"- eta=0.5 Pearson corr(thick_rod_shear_fraction, abs_diff): {format_float(thick_shear_corr_abs, 5)}",
        f"- eta=0.5 Pearson corr(thick_rod_shear_fraction, abs_rel_diff): {format_float(thick_shear_corr_rel, 5)}",
        f"- overlap(top 10 shear, top 10 abs diff): {overlap_abs}",
        f"- overlap(top 10 shear, top 10 relative diff): {overlap_rel}",
        "",
        "Top 10 rows by shear fraction:",
        "",
        *table_from_rows(
            top_shear,
            [
                ("eta", "eta"),
                ("branch", "branch_index"),
                ("mu", "mu"),
                ("shear", "shear_fraction"),
                ("rod1 energy", "rod1_energy_fraction"),
                ("rod2 energy", "rod2_energy_fraction"),
                ("abs diff", "abs_diff"),
                ("abs rel diff", "abs_rel_diff"),
            ],
        ),
        "",
        "Top 10 rows by absolute EB/Timoshenko Lambda difference:",
        "",
        *table_from_rows(
            top_abs,
            [
                ("eta", "eta"),
                ("branch", "branch_index"),
                ("mu", "mu"),
                ("abs diff", "abs_diff"),
                ("abs rel diff", "abs_rel_diff"),
                ("shear", "shear_fraction"),
                ("thick shear", "thick_rod_shear_fraction"),
            ],
        ),
        "",
        "Top 10 rows by relative EB/Timoshenko Lambda difference:",
        "",
        *table_from_rows(
            top_rel,
            [
                ("eta", "eta"),
                ("branch", "branch_index"),
                ("mu", "mu"),
                ("abs rel diff", "abs_rel_diff"),
                ("abs diff", "abs_diff"),
                ("shear", "shear_fraction"),
                ("thick shear", "thick_rod_shear_fraction"),
            ],
        ),
        "",
        "## Limitations",
        "",
        "- Descendant tracking is still diagnostic and MAC-based.",
        "- No FEM was used.",
        "- There is no eta!=0 FEM validation in this audit.",
        "- High-frequency or near-cutoff rows are not promoted; cutoff margins are only diagnostics.",
        "- The energy partition uses numerically reconstructed Timoshenko modes, not a symbolic proof.",
        "",
        "## Next Recommendation",
        "",
        "Use this CSV/report as an explanatory diagnostic for EB/Timoshenko closeness. If article-level "
        "claims are needed later, follow with targeted branch-identity review and eta!=0 FEM validation "
        "rather than promoting this audit directly.",
        "",
        "## Outputs",
        "",
        f"- CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
        f"- report: `{OUTPUT_REPORT.relative_to(REPO_ROOT)}`",
        f"- shear diagnostic PNG: `{OUTPUT_SHEAR_PNG.relative_to(REPO_ROOT)}`",
        f"- rod-energy diagnostic PNG: `{OUTPUT_ROD_PNG.relative_to(REPO_ROOT)}`",
        "",
        "## Tracking Warnings",
        "",
    ]
    if warnings_unique:
        lines.extend(f"- {item}" for item in warnings_unique[:80])
        if len(warnings_unique) > 80:
            lines.append(f"- ... {len(warnings_unique) - 80} more warnings")
    else:
        lines.append("- none")
    OUTPUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "h1_status": h1_status,
        "h1_failures": tuple(h1_failures),
        "h2_status": h2_status,
        "h3_status": h3_status,
        "max_cutoff": finite_row_value(max_cutoff, "max_Omega_over_cutoff"),
        "cutoff_warnings": app["cutoff_warning_count"],
        "cutoff_violations": app["cutoff_violation_count"],
        "shear_corr_abs": shear_corr_abs,
        "thick_shear_corr_abs": thick_shear_corr_abs,
    }


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    tracking_warnings: list[str] = []
    timo_tracking: dict[float, TrackingResult] = {}
    eb_tracking: dict[float, TrackingResult] = {}
    for eta in ETA_VALUES:
        timo = track_descendants(model="Timoshenko", eta=float(eta), solver=solve_timo_states)
        eb = track_descendants(model="Euler-Bernoulli", eta=float(eta), solver=solve_eb_states)
        timo_tracking[float(eta)] = timo
        eb_tracking[float(eta)] = eb
        tracking_warnings.extend(timo.warnings)
        tracking_warnings.extend(eb.warnings)

    rows = build_rows(timo_tracking, eb_tracking)
    write_csv(rows)
    plot_shear_fraction(rows)
    plot_rod_energy_fraction(rows)
    summary = write_report(rows, tracking_warnings)

    print("Timoshenko energy partition audit")
    print(f"Rows: {len(rows)}")
    print(f"H1 status: {summary['h1_status']}; failing branches: {summary['h1_failures'] or 'none'}")
    print(f"H2 status: {summary['h2_status']}")
    print(
        "H3 status: "
        f"{summary['h3_status']}; corr(shear, abs diff)={summary['shear_corr_abs']:.6g}; "
        f"corr(thick shear, abs diff)={summary['thick_shear_corr_abs']:.6g}"
    )
    print(
        "Cutoff: "
        f"max Omega/Omega_c={summary['max_cutoff']:.6g}; "
        f"warnings={summary['cutoff_warnings']}; violations={summary['cutoff_violations']}"
    )
    print(f"Wrote {OUTPUT_CSV.relative_to(REPO_ROOT)}")
    print(f"Wrote {OUTPUT_REPORT.relative_to(REPO_ROOT)}")
    print(f"Wrote {OUTPUT_SHEAR_PNG.relative_to(REPO_ROOT)}")
    print(f"Wrote {OUTPUT_ROD_PNG.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
