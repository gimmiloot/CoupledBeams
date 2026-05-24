from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.linalg import eigh


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    assemble_clamped_coupled_matrix_eta,
    find_first_n_roots_eta,
    thickness_mismatch_factors,
    thickness_to_length_ratios,
)
from my_project.fem import python_fem as baseline_fem  # noqa: E402
from scripts.lib.thickness_mismatch_mac_tracking import track_mu_branches_shape_mac  # noqa: E402


# =========================
# User-editable parameters
# =========================
BETA_DEG = 15.0
EPSILON = 0.0025
ETA = 0.5
MU_VALUES = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90)

NUM_MODES = 10
NUM_TRACKED_BRANCHES = 10
NUM_SORTED_ROOTS = 14
TRACK_MU_STEP = 0.01
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0

FEM_ELEMENTS_PER_ROD = 160
FEM_CONVERGENCE_ELEMENTS = (80, 160)
THICKNESS_RATIO_LIMIT = 0.1

L_TOTAL = 2.0
L_BASE = L_TOTAL / 2.0
E = 2.1e11
RHO = 7800.0

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_CSV = OUTPUT_DIR / "thickness_mismatch_fem_check_beta15_eps0p0025_eta_p0p5.csv"
OUTPUT_PNG = OUTPUT_DIR / "thickness_mismatch_fem_check_beta15_eps0p0025_eta_p0p5.png"
OUTPUT_REPORT = OUTPUT_DIR / "thickness_mismatch_fem_check_beta15_eps0p0025_eta_p0p5_report.md"

NUM_ANALYTIC_SHAPE_SAMPLES = 401
MAC_WARNING_THRESHOLD = 0.9
VALID_TOL = 1e-12
NEAR_ZERO_NORM = 1e-14
MAC_BRANCHES = NUM_TRACKED_BRANCHES


@dataclass(frozen=True)
class Geometry:
    mu: float
    eta: float
    tau1: float
    tau2: float
    length1_nd: float
    length2_nd: float
    r0: float
    r1: float
    r2: float
    area1_factor: float
    area2_factor: float
    inertia1_factor: float
    inertia2_factor: float
    mass_factor: float
    thickness_ratio_1: float
    thickness_ratio_2: float

    @property
    def validity_status(self) -> str:
        if (
            self.thickness_ratio_1 <= THICKNESS_RATIO_LIMIT + VALID_TOL
            and self.thickness_ratio_2 <= THICKNESS_RATIO_LIMIT + VALID_TOL
        ):
            return "valid"
        return "criterion_violated"


@dataclass(frozen=True)
class FemSolveResult:
    elements_per_rod: int
    omega_nd: np.ndarray
    lambda_values: np.ndarray
    vecs_free: np.ndarray
    vecs_full: np.ndarray
    free_dofs: np.ndarray


def number_tag(value: float) -> str:
    return f"{float(value):.10g}".replace("-", "m").replace(".", "p")


def geometry_for(mu: float) -> Geometry:
    factors = thickness_mismatch_factors(float(mu), ETA)
    r0 = 2.0 * EPSILON * L_BASE
    thickness_ratio_1, thickness_ratio_2 = thickness_to_length_ratios(
        EPSILON,
        float(mu),
        ETA,
    )
    return Geometry(
        mu=float(mu),
        eta=float(ETA),
        tau1=float(factors.tau1),
        tau2=float(factors.tau2),
        length1_nd=1.0 - float(mu),
        length2_nd=1.0 + float(mu),
        r0=float(r0),
        r1=float(r0 * factors.tau1),
        r2=float(r0 * factors.tau2),
        area1_factor=float(factors.tau1**2),
        area2_factor=float(factors.tau2**2),
        inertia1_factor=float(factors.tau1**4),
        inertia2_factor=float(factors.tau2**4),
        mass_factor=float(factors.mass_factor),
        thickness_ratio_1=float(thickness_ratio_1),
        thickness_ratio_2=float(thickness_ratio_2),
    )


def elem_k_param(Le: float, *, ea_nd: float, ei_nd: float) -> np.ndarray:
    """Same Euler-Bernoulli frame element form as src/my_project/fem/python_fem.py."""
    K = np.zeros((6, 6), dtype=float)

    ea = float(ea_nd) / float(Le)
    K[0, 0] += ea
    K[0, 3] -= ea
    K[3, 0] -= ea
    K[3, 3] += ea

    c = float(ei_nd) / float(Le) ** 3
    Kb = c * np.array(
        [
            [12.0, 6.0 * Le, -12.0, 6.0 * Le],
            [6.0 * Le, 4.0 * Le**2, -6.0 * Le, 2.0 * Le**2],
            [-12.0, -6.0 * Le, 12.0, -6.0 * Le],
            [6.0 * Le, 2.0 * Le**2, -6.0 * Le, 4.0 * Le**2],
        ],
        dtype=float,
    )
    for i, di in enumerate([1, 2, 4, 5]):
        for j, dj in enumerate([1, 2, 4, 5]):
            K[di, dj] += Kb[i, j]
    return K


def elem_m_param(Le: float, *, rhoa_nd: float) -> np.ndarray:
    """Same consistent axial + bending mass form as src/my_project/fem/python_fem.py."""
    M = np.zeros((6, 6), dtype=float)

    ma = float(rhoa_nd) * float(Le) / 6.0
    M[0, 0] += 2.0 * ma
    M[0, 3] += ma
    M[3, 0] += ma
    M[3, 3] += 2.0 * ma

    c = float(rhoa_nd) * float(Le) / 420.0
    Mb = c * np.array(
        [
            [156.0, 22.0 * Le, 54.0, -13.0 * Le],
            [22.0 * Le, 4.0 * Le**2, 13.0 * Le, -3.0 * Le**2],
            [54.0, 13.0 * Le, 156.0, -22.0 * Le],
            [-13.0 * Le, -3.0 * Le**2, -22.0 * Le, 4.0 * Le**2],
        ],
        dtype=float,
    )
    for i, di in enumerate([1, 2, 4, 5]):
        for j, dj in enumerate([1, 2, 4, 5]):
            M[di, dj] += Mb[i, j]
    return M


def assemble_mismatch_fem(mu: float, elements_per_rod: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    geom = geometry_for(float(mu))
    n1 = n2 = int(elements_per_rod)
    le1 = geom.length1_nd / n1
    le2 = geom.length2_nd / n2
    ndof = 3 * (n1 + n2 + 1)
    K = np.zeros((ndof, ndof), dtype=float)
    M = np.zeros((ndof, ndof), dtype=float)

    def assemble(dofs: Sequence[int], Ke: np.ndarray, Me: np.ndarray) -> None:
        for r_idx, dr in enumerate(dofs):
            for c_idx, dc in enumerate(dofs):
                K[dr, dc] += Ke[r_idx, c_idx]
                M[dr, dc] += Me[r_idx, c_idx]

    ea1 = geom.area1_factor / EPSILON**2
    ea2 = geom.area2_factor / EPSILON**2
    ke1 = elem_k_param(le1, ea_nd=ea1, ei_nd=geom.inertia1_factor)
    me1 = elem_m_param(le1, rhoa_nd=geom.area1_factor)
    ke2_local = elem_k_param(le2, ea_nd=ea2, ei_nd=geom.inertia2_factor)
    me2_local = elem_m_param(le2, rhoa_nd=geom.area2_factor)
    transform = baseline_fem.rotation_matrix_6x6(np.deg2rad(BETA_DEG))
    ke2 = transform @ ke2_local @ transform.T
    me2 = transform @ me2_local @ transform.T

    for i in range(n1):
        dofs = [3 * i, 3 * i + 1, 3 * i + 2, 3 * (i + 1), 3 * (i + 1) + 1, 3 * (i + 1) + 2]
        assemble(dofs, ke1, me1)

    for j in range(n2):
        dofs = [
            3 * (n1 + j),
            3 * (n1 + j) + 1,
            3 * (n1 + j) + 2,
            3 * (n1 + j + 1),
            3 * (n1 + j + 1) + 1,
            3 * (n1 + j + 1) + 2,
        ]
        assemble(dofs, ke2, me2)

    fixed = {0, 1, 2, 3 * (n1 + n2), 3 * (n1 + n2) + 1, 3 * (n1 + n2) + 2}
    free = np.array(sorted(set(range(ndof)) - fixed), dtype=int)
    return K[np.ix_(free, free)], M[np.ix_(free, free)], free


def solve_mismatch_fem(mu: float, elements_per_rod: int, n_modes: int) -> FemSolveResult:
    Kf, Mf, free = assemble_mismatch_fem(float(mu), int(elements_per_rod))
    eigenvalues, eigenvectors = eigh(Kf, Mf, subset_by_index=[0, int(n_modes) - 1])
    omega_nd = np.sqrt(np.maximum(eigenvalues, 0.0))
    lambda_values = np.sqrt(np.maximum(omega_nd, 0.0))

    ndof = 3 * (2 * int(elements_per_rod) + 1)
    full = np.zeros((ndof, int(n_modes)), dtype=float)
    full[free, :] = eigenvectors
    return FemSolveResult(
        elements_per_rod=int(elements_per_rod),
        omega_nd=omega_nd,
        lambda_values=lambda_values,
        vecs_free=eigenvectors,
        vecs_full=full,
        free_dofs=free,
    )


def sorted_roots_for(mu: float, n_roots: int = NUM_SORTED_ROOTS) -> np.ndarray:
    roots = find_first_n_roots_eta(
        float(np.deg2rad(BETA_DEG)),
        float(mu),
        EPSILON,
        ETA,
        int(n_roots),
        Lmax0=ROOT_LMAX0,
        scan_step=ROOT_SCAN_STEP,
    )
    if np.any(~np.isfinite(roots[:NUM_MODES])):
        raise RuntimeError(f"Missing analytic roots for mu={mu:g}, eta={ETA:g}.")
    return roots


def tracking_mu_grid() -> np.ndarray:
    target_max = max(float(mu) for mu in MU_VALUES)
    grid = np.arange(0.0, target_max + 0.5 * TRACK_MU_STEP, TRACK_MU_STEP, dtype=float)
    values = sorted({round(float(value), 10) for value in grid} | {round(float(mu), 10) for mu in MU_VALUES} | {0.0})
    return np.asarray(values, dtype=float)


def track_analytic_branches() -> tuple[dict[float, np.ndarray], dict[float, np.ndarray], list[dict[str, float | int | str]]]:
    mu_grid = tracking_mu_grid()
    roots_by_mu = {float(mu): sorted_roots_for(float(mu)) for mu in mu_grid}
    result = track_mu_branches_shape_mac(
        beta_rad=float(np.deg2rad(BETA_DEG)),
        epsilon=EPSILON,
        eta=ETA,
        mu_values=mu_grid,
        roots_by_mu=roots_by_mu,
        num_tracked_branches=NUM_TRACKED_BRANCHES,
        num_shape_samples=NUM_ANALYTIC_SHAPE_SAMPLES,
        mac_warning_threshold=MAC_WARNING_THRESHOLD,
    )
    tracked = result.tracked
    sorted_indices = result.current_sorted_indices
    rows: list[dict[str, float | int | str]] = []
    for row in result.rows:
        rows.append(
            {
                "mu": float(row["mu"]),
                "branch_index": int(row["branch_index_from_mu0"]),
                "Lambda_tracked": float(row["Lambda_tracked"]),
                "current_sorted_root_index": int(row["mac_sorted_root_index"]),
                "raw_mac_sorted_root_index": int(row["raw_mac_sorted_root_index"]),
                "frequency_nearest_sorted_root_index": int(row["nearest_sorted_root_index"]),
                # Backward-compatible alias for report code that expects a sorted-index field.
                "nearest_sorted_root_index": int(row["mac_sorted_root_index"]),
                "tracking_step_status": str(row["tracking_step_status"]),
                "mac_to_previous": float(row["mac_to_previous"]),
                "raw_mac_to_previous": float(row["raw_mac_to_previous"]),
                "second_best_mac": float(row["second_best_mac"]),
                "frequency_mac_disagreement": str(row["frequency_mac_disagreement"]),
                "used_frequency_fallback": str(row["used_frequency_fallback"]),
                "assignment_margin": float(row["frequency_assignment_margin"]),
            }
        )

    tracked_targets: dict[float, np.ndarray] = {}
    sorted_index_targets: dict[float, np.ndarray] = {}
    for target_mu in MU_VALUES:
        matches = np.where(np.isclose(mu_grid, float(target_mu), rtol=0.0, atol=1e-12))[0]
        if len(matches) != 1:
            raise RuntimeError(f"Tracking grid does not contain target mu={float(target_mu):g}.")
        col = int(matches[0])
        tracked_targets[float(target_mu)] = tracked[:, col].copy()
        sorted_index_targets[float(target_mu)] = sorted_indices[:, col].copy()
    return tracked_targets, sorted_index_targets, rows


def analytic_null_vector_eta(Lambda: float, mu: float) -> np.ndarray:
    matrix = assemble_clamped_coupled_matrix_eta(float(Lambda), np.deg2rad(BETA_DEG), float(mu), EPSILON, ETA)
    _, _singular_values, vh = np.linalg.svd(matrix)
    return vh[-1, :].astype(float)


def reconstruct_components_eta(
    Lambda: float,
    *,
    mu: float,
    coeff: np.ndarray,
    s_norm: np.ndarray,
) -> dict[str, np.ndarray]:
    factors = thickness_mismatch_factors(float(mu), ETA)
    A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]
    xi = np.asarray(s_norm, dtype=float)

    z1 = float(Lambda) * (1.0 - float(mu)) * xi / np.sqrt(factors.tau1)
    theta1 = EPSILON * float(Lambda) ** 2 * (1.0 - float(mu)) * xi
    w_left = A1 * (np.cos(z1) - np.cosh(z1)) + B1 * (np.sin(z1) - np.sinh(z1))
    u_left = P1 * np.sin(theta1)

    z2_external_to_joint = -float(Lambda) * (1.0 + float(mu)) * xi / np.sqrt(factors.tau2)
    theta2_external_to_joint = -EPSILON * float(Lambda) ** 2 * (1.0 + float(mu)) * xi
    w_right = A2 * (np.cos(z2_external_to_joint) - np.cosh(z2_external_to_joint)) + B2 * (
        np.sin(z2_external_to_joint) - np.sinh(z2_external_to_joint)
    )
    u_right = P2 * np.sin(theta2_external_to_joint)

    return {
        "u_left": u_left,
        "w_left": w_left,
        "u_right": u_right[::-1],
        "w_right": w_right[::-1],
    }


def vector_from_components(components: dict[str, np.ndarray]) -> np.ndarray:
    return np.concatenate(
        [
            np.asarray(components["u_left"], dtype=float),
            np.asarray(components["w_left"], dtype=float),
            np.asarray(components["u_right"], dtype=float),
            np.asarray(components["w_right"], dtype=float),
        ]
    )


def analytic_shape_vector(Lambda: float, mu: float, elements_per_rod: int) -> np.ndarray:
    coeff = analytic_null_vector_eta(float(Lambda), float(mu))
    s_norm = np.linspace(0.0, 1.0, int(elements_per_rod) + 1)
    return vector_from_components(reconstruct_components_eta(float(Lambda), mu=float(mu), coeff=coeff, s_norm=s_norm))


def fem_shape_vector(result: FemSolveResult, mode_index: int) -> np.ndarray:
    n = int(result.elements_per_rod)
    full = result.vecs_full[:, int(mode_index) - 1]
    left = full[: 3 * (n + 1)].reshape((n + 1, 3))
    right_global = full[3 * n : 3 * (2 * n + 1)].reshape((n + 1, 3))

    c = np.cos(np.deg2rad(BETA_DEG))
    s = np.sin(np.deg2rad(BETA_DEG))
    rotation = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)
    right_local = right_global @ rotation

    return np.concatenate([left[:, 0], left[:, 1], right_local[:, 0], right_local[:, 1]])


def mac_value(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    numerator = abs(float(np.dot(a, b))) ** 2
    denominator = float(np.dot(a, a)) * float(np.dot(b, b))
    return numerator / denominator if denominator > NEAR_ZERO_NORM else 0.0


def mac_matrix_for_mu(
    *,
    mu: float,
    tracked_lambdas: np.ndarray,
    fem_result: FemSolveResult,
) -> np.ndarray:
    analytic_vectors = [
        analytic_shape_vector(float(tracked_lambdas[idx]), float(mu), fem_result.elements_per_rod)
        for idx in range(MAC_BRANCHES)
    ]
    fem_vectors = [fem_shape_vector(fem_result, mode_index=idx + 1) for idx in range(NUM_MODES)]
    mac = np.zeros((MAC_BRANCHES, NUM_MODES), dtype=float)
    for row, analytic_vector in enumerate(analytic_vectors):
        for col, fem_vector in enumerate(fem_vectors):
            mac[row, col] = mac_value(analytic_vector, fem_vector)
    return mac


def relative_difference(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / abs(float(right)) if abs(float(right)) > NEAR_ZERO_NORM else np.nan


def write_csv_rows(path: Path, rows: Iterable[dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_rows(
    *,
    sorted_roots: dict[float, np.ndarray],
    tracked_roots: dict[float, np.ndarray],
    tracked_sorted_indices: dict[float, np.ndarray],
    fem_by_mesh: dict[int, dict[float, FemSolveResult]],
    mac_by_mu: dict[float, np.ndarray],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    fine_mesh = int(FEM_ELEMENTS_PER_ROD)
    coarse_mesh = int(FEM_CONVERGENCE_ELEMENTS[0])
    for mu in MU_VALUES:
        geom = geometry_for(float(mu))
        sorted_values = sorted_roots[float(mu)]
        tracked_values = tracked_roots[float(mu)]
        tracked_indices = tracked_sorted_indices[float(mu)]
        fem_fine = fem_by_mesh[fine_mesh][float(mu)].lambda_values
        fem_coarse = fem_by_mesh[coarse_mesh][float(mu)].lambda_values
        mac = mac_by_mu[float(mu)]
        best_branch_for_fem = np.argmax(mac, axis=0) + 1
        best_branch_mac = np.max(mac, axis=0)
        best_fem_for_branch = np.argmax(mac, axis=1) + 1
        best_fem_mac = np.max(mac, axis=1)

        for mode_idx in range(1, NUM_MODES + 1):
            sorted_lambda = float(sorted_values[mode_idx - 1])
            tracked_lambda = float(tracked_values[mode_idx - 1]) if mode_idx <= NUM_TRACKED_BRANCHES else np.nan
            fem_lambda = float(fem_fine[mode_idx - 1])
            coarse_lambda = float(fem_coarse[mode_idx - 1])
            same_index_mac = float(mac[mode_idx - 1, mode_idx - 1]) if mode_idx <= MAC_BRANCHES else np.nan
            best_fem_mode = int(best_fem_for_branch[mode_idx - 1]) if mode_idx <= MAC_BRANCHES else ""
            best_fem_mode_mac = float(best_fem_mac[mode_idx - 1]) if mode_idx <= MAC_BRANCHES else np.nan
            rows.append(
                {
                    "beta_deg": BETA_DEG,
                    "epsilon": EPSILON,
                    "eta": ETA,
                    "mu": float(mu),
                    "mode_index": int(mode_idx),
                    "Lambda_analytic_sorted": sorted_lambda,
                    "Lambda_analytic_tracked_if_available": tracked_lambda,
                    "Lambda_fem": fem_lambda,
                    "abs_diff_sorted": abs(sorted_lambda - fem_lambda),
                    "rel_diff_sorted": relative_difference(sorted_lambda, fem_lambda),
                    "abs_diff_tracked_if_available": abs(tracked_lambda - fem_lambda)
                    if np.isfinite(tracked_lambda)
                    else np.nan,
                    "rel_diff_tracked_if_available": relative_difference(tracked_lambda, fem_lambda)
                    if np.isfinite(tracked_lambda)
                    else np.nan,
                    "thickness_ratio_1": geom.thickness_ratio_1,
                    "thickness_ratio_2": geom.thickness_ratio_2,
                    "validity_status": geom.validity_status,
                    "fem_mode_index": int(mode_idx),
                    "analytic_branch_index": int(mode_idx) if mode_idx <= MAC_BRANCHES else "",
                    "MAC": same_index_mac,
                    "best_match_branch": int(best_branch_for_fem[mode_idx - 1]),
                    "best_match_branch_MAC": float(best_branch_mac[mode_idx - 1]),
                    "best_match_fem_mode_for_tracked": best_fem_mode,
                    "best_match_fem_mode_MAC": best_fem_mode_mac,
                    "Lambda_fem_coarse": coarse_lambda,
                    "Lambda_fem_refined": fem_lambda,
                    "fem_mesh_abs_change": abs(fem_lambda - coarse_lambda),
                    "fem_mesh_rel_change": relative_difference(fem_lambda, coarse_lambda),
                    "tracked_mac_sorted_root_index": int(tracked_indices[mode_idx - 1])
                    if mode_idx <= NUM_TRACKED_BRANCHES
                    else "",
                    "tracked_nearest_sorted_root_index": int(tracked_indices[mode_idx - 1])
                    if mode_idx <= NUM_TRACKED_BRANCHES
                    else "",
                }
            )
    return rows


CSV_FIELDNAMES = [
    "beta_deg",
    "epsilon",
    "eta",
    "mu",
    "mode_index",
    "Lambda_analytic_sorted",
    "Lambda_analytic_tracked_if_available",
    "Lambda_fem",
    "abs_diff_sorted",
    "rel_diff_sorted",
    "abs_diff_tracked_if_available",
    "rel_diff_tracked_if_available",
    "thickness_ratio_1",
    "thickness_ratio_2",
    "validity_status",
    "fem_mode_index",
    "analytic_branch_index",
    "MAC",
    "best_match_branch",
    "best_match_branch_MAC",
    "best_match_fem_mode_for_tracked",
    "best_match_fem_mode_MAC",
    "Lambda_fem_coarse",
    "Lambda_fem_refined",
    "fem_mesh_abs_change",
    "fem_mesh_rel_change",
    "tracked_mac_sorted_root_index",
    "tracked_nearest_sorted_root_index",
]


def plot_comparison(
    *,
    tracked_roots: dict[float, np.ndarray],
    fem_by_mesh: dict[int, dict[float, FemSolveResult]],
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    mu_values = np.asarray(MU_VALUES, dtype=float)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, ax = plt.subplots(figsize=(9.8, 5.5))

    for branch_idx in range(6):
        linewidth = 2.8 if branch_idx in (4, 5) else 1.8
        linestyle = "-" if branch_idx not in (4, 5) else ("-" if branch_idx == 4 else "--")
        label = f"tracked branch {branch_idx + 1}"
        ax.plot(
            mu_values,
            [tracked_roots[float(mu)][branch_idx] for mu in MU_VALUES],
            color=colors[branch_idx % len(colors)],
            lw=linewidth,
            ls=linestyle,
            label=label,
        )

    fine = fem_by_mesh[int(FEM_ELEMENTS_PER_ROD)]
    for mode_idx in range(1, NUM_MODES + 1):
        if mode_idx <= 6:
            color = colors[(mode_idx - 1) % len(colors)]
            marker = "o" if mode_idx not in (5, 6) else ("s" if mode_idx == 5 else "D")
            size = 34 if mode_idx not in (5, 6) else 52
            zorder = 3 if mode_idx not in (5, 6) else 5
        else:
            color = "0.55"
            marker = "x"
            size = 30
            zorder = 2
        ax.scatter(
            mu_values,
            [fine[float(mu)].lambda_values[mode_idx - 1] for mu in MU_VALUES],
            color=color,
            marker=marker,
            s=size,
            zorder=zorder,
        )

    style_handles = [
        Line2D([0], [0], color="black", lw=2.0, label="analytic MAC-tracked"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=5.5, label="FEM sorted modes"),
        Line2D([0], [0], color=colors[4], lw=2.8, label="branch/mode 5 highlighted"),
        Line2D([0], [0], color=colors[5], lw=2.8, ls="--", label="branch/mode 6 highlighted"),
    ]

    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\Lambda$")
    ax.set_title(rf"Thickness-mismatch FEM check, $\eta={ETA:g}$, $\beta={BETA_DEG:g}^\circ$, $\epsilon={EPSILON:g}$")
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(handles=style_handles, loc="best", fontsize=9, frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def violation_summary() -> tuple[list[tuple[float, tuple[int, ...]]], float, int, float]:
    violations: list[tuple[float, tuple[int, ...]]] = []
    max_ratio = -np.inf
    max_rod = 1
    max_mu = float(MU_VALUES[0])
    for mu in MU_VALUES:
        geom = geometry_for(float(mu))
        ratios = (geom.thickness_ratio_1, geom.thickness_ratio_2)
        rods = tuple(idx + 1 for idx, ratio in enumerate(ratios) if ratio > THICKNESS_RATIO_LIMIT + VALID_TOL)
        if rods:
            violations.append((float(mu), rods))
        local_rod = int(np.argmax(ratios)) + 1
        local_max = float(max(ratios))
        if local_max > max_ratio:
            max_ratio = local_max
            max_rod = local_rod
            max_mu = float(mu)
    return violations, float(max_ratio), int(max_rod), float(max_mu)


def max_error(rows: Sequence[dict[str, object]], key: str) -> tuple[float, dict[str, object]]:
    finite_rows = [row for row in rows if np.isfinite(float(row[key]))]
    row = max(finite_rows, key=lambda item: float(item[key]))
    return float(row[key]), row


def branch_pair_mac_summary(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    summary = []
    for mu in MU_VALUES:
        local = [
            row
            for row in rows
            if abs(float(row["mu"]) - float(mu)) < 1e-12 and int(row["mode_index"]) in (5, 6, 7)
        ]
        for row in local:
            summary.append(
                {
                    "mu": float(mu),
                    "branch": int(row["mode_index"]),
                    "tracked_lambda": float(row["Lambda_analytic_tracked_if_available"]),
                    "same_index_fem_lambda": float(row["Lambda_fem"]),
                    "same_index_mac": float(row["MAC"]),
                    "best_fem_mode": row["best_match_fem_mode_for_tracked"],
                    "best_fem_mac": float(row["best_match_fem_mode_MAC"]),
                    "best_branch_for_same_fem": int(row["best_match_branch"]),
                    "best_branch_mac_for_same_fem": float(row["best_match_branch_MAC"]),
                }
            )
    return summary


def branch56_swap_conclusion(pair_summary: Sequence[dict[str, object]]) -> str:
    high_mu = [row for row in pair_summary if float(row["mu"]) >= 0.70 and int(row["branch"]) in (5, 6)]
    if not high_mu:
        return "No high-mu MAC rows were available."
    by_mu = {}
    for row in pair_summary:
        by_mu.setdefault(float(row["mu"]), {})[int(row["branch"])] = int(row["best_fem_mode"])

    pair_swap_mus = [
        mu
        for mu, mapping in by_mu.items()
        if mu >= 0.70 and mapping.get(5) == 6 and mapping.get(6) == 5
    ]
    cyclic_mus = [
        mu
        for mu, mapping in by_mu.items()
        if mu >= 0.70 and mapping.get(5) == 7 and mapping.get(6) == 5 and mapping.get(7) == 6
    ]
    if pair_swap_mus and cyclic_mus:
        return (
            "The shape MAC diagnostic suggests possible assignment ambiguity among "
            f"branches 5--7 near mu={', '.join(f'{mu:g}' for mu in sorted(set(pair_swap_mus + cyclic_mus)))}. "
            "This must be checked on a refined local mu grid before being interpreted "
            "as a real change of sorted position."
        )
    if pair_swap_mus:
        return (
            "The shape MAC diagnostic suggests a possible local branch-5/branch-6 "
            f"assignment ambiguity at mu={', '.join(f'{mu:g}' for mu in pair_swap_mus)}. "
            "Interpret it only after a refined local mu-grid check."
        )
    if cyclic_mus:
        return (
            "The shape MAC diagnostic suggests possible ambiguity involving branches "
            f"5--7 at mu={', '.join(f'{mu:g}' for mu in cyclic_mus)}. "
            "Do not interpret this as automatic branch renumbering without a refined "
            "local mu-grid audit."
        )
    branch5 = [row for row in high_mu if int(row["branch"]) == 5]
    branch6 = [row for row in high_mu if int(row["branch"]) == 6]
    branch5_to_5 = sum(1 for row in branch5 if int(row["best_fem_mode"]) == 5)
    branch6_to_6 = sum(1 for row in branch6 if int(row["best_fem_mode"]) == 6)
    if branch5_to_5 >= len(branch5) / 2 and branch6_to_6 >= len(branch6) / 2:
        return (
            "The shape MAC diagnostic does not support a persistent branch-5/branch-6 swap "
            "on this grid; same-number FEM modes are usually the best matches."
        )
    return (
        "The shape MAC diagnostic is mixed for branches 5 and 6; branch identity should be "
        "reviewed with a denser shape-tracking run before drawing a conclusion."
    )


def format_table(rows: Sequence[Sequence[str]]) -> list[str]:
    if not rows:
        return []
    widths = [max(len(row[col]) for row in rows) for col in range(len(rows[0]))]
    lines = []
    for idx, row in enumerate(rows):
        lines.append("| " + " | ".join(value.ljust(widths[col]) for col, value in enumerate(row)) + " |")
        if idx == 0:
            lines.append("| " + " | ".join("-" * widths[col] for col in range(len(row))) + " |")
    return lines


def build_report(
    *,
    rows: Sequence[dict[str, object]],
    tracking_rows: Sequence[dict[str, float | int | str]],
    pair_summary: Sequence[dict[str, object]],
) -> None:
    OUTPUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    max_abs_sorted, max_abs_sorted_row = max_error(rows, "abs_diff_sorted")
    max_rel_sorted, max_rel_sorted_row = max_error(rows, "rel_diff_sorted")
    max_abs_tracked, max_abs_tracked_row = max_error(rows, "abs_diff_tracked_if_available")
    max_rel_tracked, max_rel_tracked_row = max_error(rows, "rel_diff_tracked_if_available")
    max_mesh_change, max_mesh_row = max_error(rows, "fem_mesh_abs_change")
    max_mesh_rel, max_mesh_rel_row = max_error(rows, "fem_mesh_rel_change")
    violations, max_ratio, max_ratio_rod, max_ratio_mu = violation_summary()
    low_mac = [row for row in tracking_rows if "low_mac" in str(row["tracking_step_status"])]
    frequency_warnings = [row for row in tracking_rows if str(row["frequency_mac_disagreement"]) == "yes"]
    switches = [
        row
        for row in tracking_rows
        if int(row["current_sorted_root_index"]) != int(row["branch_index"])
    ]

    lines = [
        "# Thickness-Mismatch FEM Check: eta=0.5, beta=15",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g}",
        f"- mu values: {', '.join(f'{float(mu):g}' for mu in MU_VALUES)}",
        f"- modes compared: first {NUM_MODES}",
        f"- tracked analytic branches: first {NUM_TRACKED_BRANCHES} from mu=0 by analytic shape MAC with low-MAC warnings",
        "",
        "## FEM Formulation",
        "",
        "This diagnostic reuses the existing planar Euler-Bernoulli frame-element",
        "form from `src/my_project/fem/python_fem.py`: axial bar stiffness,",
        "Euler-Bernoulli bending stiffness, consistent axial and bending mass,",
        "shared joint DOFs, and clamped external ends. The right arm uses the",
        "same local-to-global convention `q_global = T @ q_local` and assembly",
        "`K_global = T @ K_local @ T.T`, `M_global = T @ M_local @ T.T`.",
        "",
        "The baseline FEM module is not modified. This script only parameterizes",
        "the same element matrices by arm-specific circular-section factors:",
        "`EA_i ~ tau_i^2 / epsilon^2`, `EI_i ~ tau_i^4`, `rhoA_i ~ tau_i^2`.",
        "",
        "## Mesh",
        "",
        f"- primary mesh: {FEM_ELEMENTS_PER_ROD} elements per rod",
        f"- convergence meshes: {', '.join(str(value) for value in FEM_CONVERGENCE_ELEMENTS)} elements per rod",
        f"- max coarse-to-refined abs Lambda change: {max_mesh_change:.6e} "
        f"at mu={float(max_mesh_row['mu']):g}, mode={int(max_mesh_row['mode_index'])}",
        f"- max coarse-to-refined rel Lambda change: {max_mesh_rel:.6e} "
        f"at mu={float(max_mesh_rel_row['mu']):g}, mode={int(max_mesh_rel_row['mode_index'])}",
        "",
        "## Lambda Normalization",
        "",
        "The dimensionless FEM eigenproblem is scaled with the equal-radius",
        "reference section `S0=pi*r0^2`, `J0=pi*r0^4/4`, and `l=L_total/2`.",
        "Its returned angular-frequency parameter is",
        "`omega_nd = Omega*l^2*sqrt(rho*S0/(E*J0)) = Lambda^2`, so the plotted",
        "FEM value is `Lambda_fem = sqrt(omega_nd)`. Equivalently, if the",
        "physical angular frequency is `Omega`, then",
        "`Lambda_fem = sqrt(Omega)*l*(rho*S0/(E*J0))^(1/4)`.",
        "",
        "## Frequency Comparison",
        "",
        f"- max sorted abs difference: {max_abs_sorted:.6e} "
        f"at mu={float(max_abs_sorted_row['mu']):g}, mode={int(max_abs_sorted_row['mode_index'])}",
        f"- max sorted rel difference: {max_rel_sorted:.6e} "
        f"at mu={float(max_rel_sorted_row['mu']):g}, mode={int(max_rel_sorted_row['mode_index'])}",
        f"- max tracked abs difference: {max_abs_tracked:.6e} "
        f"at mu={float(max_abs_tracked_row['mu']):g}, branch={int(max_abs_tracked_row['mode_index'])}",
        f"- max tracked rel difference: {max_rel_tracked:.6e} "
        f"at mu={float(max_rel_tracked_row['mu']):g}, branch={int(max_rel_tracked_row['mode_index'])}",
        f"- analytic tracking low-MAC rows on full 0..0.9 grid: {len(low_mac)}",
        f"- analytic MAC current-sorted-index switch rows on full 0..0.9 grid: {len(switches)}",
        f"- nearest-frequency vs MAC disagreement rows on full 0..0.9 grid: {len(frequency_warnings)}",
        "",
        "## Branches 5, 6, and Nearby Mode 7",
        "",
        branch56_swap_conclusion(pair_summary),
        "",
    ]

    pair_table = [["mu", "branch", "tracked Lambda", "best FEM mode", "best MAC", "same-index MAC"]]
    for row in pair_summary:
        pair_table.append(
            [
                f"{float(row['mu']):g}",
                str(int(row["branch"])),
                f"{float(row['tracked_lambda']):.6g}",
                str(row["best_fem_mode"]),
                f"{float(row['best_fem_mac']):.6f}",
                f"{float(row['same_index_mac']):.6f}",
            ]
        )
    lines.extend(format_table(pair_table))

    lines.extend(
        [
            "",
            "## Thickness Criterion",
            "",
            f"- criterion: `2r_i/l_i <= {THICKNESS_RATIO_LIMIT:g}` for both rods",
            f"- max ratio on checked grid: {max_ratio:.6g} "
            f"(rod {max_ratio_rod}, mu={max_ratio_mu:g})",
        ]
    )
    if violations:
        lines.append("- warning: criterion violations were found:")
        for mu, rods in violations:
            lines.append(f"  - mu={mu:g}: rod(s) {', '.join(str(rod) for rod in rods)}")
    else:
        lines.append("- no violations of `2r_i/l_i <= 0.1` were found on this FEM-check grid.")

    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
            f"- PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
            "",
            "This is a diagnostic-only check. It does not change the article files,",
            "article figures, baseline equal-radius determinant, old solvers, or",
            "the existing FEM physical model.",
            "",
        ]
    )
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, object]:
    sorted_roots = {float(mu): sorted_roots_for(float(mu)) for mu in MU_VALUES}
    tracked_roots, tracked_sorted_indices, tracking_rows = track_analytic_branches()

    fem_by_mesh: dict[int, dict[float, FemSolveResult]] = {}
    for mesh in FEM_CONVERGENCE_ELEMENTS:
        fem_by_mesh[int(mesh)] = {}
        for mu in MU_VALUES:
            fem_by_mesh[int(mesh)][float(mu)] = solve_mismatch_fem(float(mu), int(mesh), NUM_MODES)

    fine_mesh = int(FEM_ELEMENTS_PER_ROD)
    mac_by_mu = {
        float(mu): mac_matrix_for_mu(
            mu=float(mu),
            tracked_lambdas=tracked_roots[float(mu)],
            fem_result=fem_by_mesh[fine_mesh][float(mu)],
        )
        for mu in MU_VALUES
    }

    rows = build_rows(
        sorted_roots=sorted_roots,
        tracked_roots=tracked_roots,
        tracked_sorted_indices=tracked_sorted_indices,
        fem_by_mesh=fem_by_mesh,
        mac_by_mu=mac_by_mu,
    )
    write_csv_rows(OUTPUT_CSV, rows, CSV_FIELDNAMES)
    plot_comparison(tracked_roots=tracked_roots, fem_by_mesh=fem_by_mesh, output=OUTPUT_PNG)
    pair_summary = branch_pair_mac_summary(rows)
    build_report(rows=rows, tracking_rows=tracking_rows, pair_summary=pair_summary)

    max_abs_sorted, max_abs_sorted_row = max_error(rows, "abs_diff_sorted")
    max_rel_sorted, max_rel_sorted_row = max_error(rows, "rel_diff_sorted")
    max_mesh_change, max_mesh_row = max_error(rows, "fem_mesh_abs_change")
    violations, max_ratio, max_ratio_rod, max_ratio_mu = violation_summary()

    print(f"saved CSV: {OUTPUT_CSV}")
    print(f"saved PNG: {OUTPUT_PNG}")
    print(f"saved report: {OUTPUT_REPORT}")
    print(
        f"max sorted abs diff: {max_abs_sorted:.6e} "
        f"(mu={float(max_abs_sorted_row['mu']):g}, mode={int(max_abs_sorted_row['mode_index'])})"
    )
    print(
        f"max sorted rel diff: {max_rel_sorted:.6e} "
        f"(mu={float(max_rel_sorted_row['mu']):g}, mode={int(max_rel_sorted_row['mode_index'])})"
    )
    print(
        f"max FEM mesh abs change: {max_mesh_change:.6e} "
        f"(mu={float(max_mesh_row['mu']):g}, mode={int(max_mesh_row['mode_index'])})"
    )
    print(branch56_swap_conclusion(pair_summary))
    if violations:
        print("warning: thickness-ratio criterion violations detected")
        for mu, rods in violations:
            print(f"  mu={mu:g}: rod(s) {', '.join(str(rod) for rod in rods)}")
    else:
        print(
            "no thickness-ratio violations detected; "
            f"max 2*r_i/l_i={max_ratio:.6g} at mu={max_ratio_mu:g}, rod={max_ratio_rod}"
        )

    return {
        "csv": OUTPUT_CSV,
        "png": OUTPUT_PNG,
        "report": OUTPUT_REPORT,
        "max_abs_sorted": max_abs_sorted,
        "max_rel_sorted": max_rel_sorted,
        "max_mesh_change": max_mesh_change,
        "max_thickness_ratio": max_ratio,
    }


if __name__ == "__main__":
    main()
