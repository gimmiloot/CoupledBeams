from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


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
)
from scripts.analysis.track_lambda_eta_thickness_mismatch import unique_nearest_assignment  # noqa: E402
from scripts.lib.analytic_coupled_rods_shapes import (  # noqa: E402
    analytic_null_vector,
    normalize_components,
)


# =========================
# User-editable parameters
# =========================
BETA_DEG = 15.0
EPSILON = 0.0025
ETA_VALUES = (-0.1, 0.0, 0.1)
MU_VALUES = (0.0, 0.1)
BRANCH_INDEX = 5

TRACK_MU_STEP = 0.01
NUM_TRACKED_BRANCHES = 6
NUM_SORTED_ROOTS = 12
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0

NUM_SAMPLES = 401
L_TOTAL = 2.0
MODE_SCALE = 0.16
NORMALIZE = "max-full"
DPI = 240
FIGSIZE = (10.8, 3.8)
OUTPUT_DIR = REPO_ROOT / "results"


RIGHT_COORDINATE_FOR_PLOTTING = "joint-to-external"


@dataclass(frozen=True)
class TrackedRoot:
    eta: float
    mu: float
    branch_index_from_mu0: int
    Lambda: float
    current_sorted_index: int
    tracking_step_status: str


@dataclass(frozen=True)
class ShapeCase:
    eta: float
    mu: float
    Lambda: float
    current_sorted_index: int
    tracking_step_status: str
    components: dict[str, np.ndarray]
    matrix_residual_max: float
    smallest_singular_value: float
    singular_value_ratio: float


def number_tag(value: float, *, force_one_decimal: bool = False) -> str:
    if force_one_decimal:
        text = f"{float(value):.1f}"
    else:
        text = f"{float(value):.10g}"
    return text.replace("-", "m").replace(".", "p")


def output_png_path(mu: float) -> Path:
    return OUTPUT_DIR / (
        f"thickness_mismatch_branch{int(BRANCH_INDEX)}_shapes_"
        f"beta{number_tag(BETA_DEG)}_eps{number_tag(EPSILON)}_"
        f"mu{number_tag(mu, force_one_decimal=True)}.png"
    )


def roots_for(mu: float, eta: float) -> np.ndarray:
    roots = find_first_n_roots_eta(
        float(np.deg2rad(BETA_DEG)),
        float(mu),
        EPSILON,
        float(eta),
        NUM_SORTED_ROOTS,
        Lmax0=ROOT_LMAX0,
        scan_step=ROOT_SCAN_STEP,
    )
    if np.any(~np.isfinite(roots)):
        raise RuntimeError(f"Missing roots for eta={eta:g}, mu={mu:g}.")
    return roots


def tracking_mu_grid(target_mu: float) -> np.ndarray:
    if float(target_mu) < -1e-12:
        raise ValueError("This diagnostic tracks from mu=0 to non-negative target mu values.")
    grid = np.arange(0.0, float(target_mu) + 0.5 * TRACK_MU_STEP, TRACK_MU_STEP, dtype=float)
    values = sorted({round(float(value), 10) for value in grid} | {0.0, round(float(target_mu), 10)})
    return np.asarray(values, dtype=float)


def track_branch_for_eta(eta: float, target_mus: Sequence[float]) -> dict[float, TrackedRoot]:
    max_mu = max(float(value) for value in target_mus)
    mu_grid = tracking_mu_grid(max_mu)
    roots_by_mu = {float(mu): roots_for(float(mu), float(eta)) for mu in mu_grid}

    tracked = np.full((NUM_TRACKED_BRANCHES, len(mu_grid)), np.nan, dtype=float)
    tracked[:, 0] = roots_by_mu[0.0][:NUM_TRACKED_BRANCHES]
    current_sorted_indices = np.arange(1, NUM_TRACKED_BRANCHES + 1, dtype=int)
    status_by_col = ["seed_mu0"]

    previous = tracked[:, 0].copy()
    sorted_indices_by_col: list[np.ndarray] = [current_sorted_indices.copy()]
    for col in range(1, len(mu_grid)):
        mu = float(mu_grid[col])
        roots = roots_by_mu[mu]
        assignments = unique_nearest_assignment(previous, roots)
        current = np.full(NUM_TRACKED_BRANCHES, np.nan, dtype=float)
        current_sorted_indices = np.full(NUM_TRACKED_BRANCHES, -1, dtype=int)
        status = "ok"
        for assignment in assignments:
            branch_row = assignment.branch_index - 1
            root_index = assignment.root_index
            current[branch_row] = float(roots[root_index - 1])
            current_sorted_indices[branch_row] = int(root_index)
            if assignment.assignment_margin <= 1e-5:
                status = "ambiguous"
        tracked[:, col] = current
        sorted_indices_by_col.append(current_sorted_indices.copy())
        status_by_col.append(status)
        previous = current

    result: dict[float, TrackedRoot] = {}
    for target_mu in target_mus:
        matches = np.where(np.isclose(mu_grid, float(target_mu), rtol=0.0, atol=1e-12))[0]
        if len(matches) != 1:
            raise RuntimeError(f"Tracking grid does not contain target mu={float(target_mu):g}.")
        col = int(matches[0])
        branch_row = int(BRANCH_INDEX) - 1
        result[float(target_mu)] = TrackedRoot(
            eta=float(eta),
            mu=float(target_mu),
            branch_index_from_mu0=int(BRANCH_INDEX),
            Lambda=float(tracked[branch_row, col]),
            current_sorted_index=int(sorted_indices_by_col[col][branch_row]),
            tracking_step_status=str(status_by_col[col]),
        )
    return result


def reconstruct_components_eta(
    Lambda: float,
    *,
    mu: float,
    eta: float,
    epsilon: float,
    coeff: np.ndarray,
    s_norm: np.ndarray,
) -> dict[str, np.ndarray]:
    factors = thickness_mismatch_factors(float(mu), float(eta))
    A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]
    xi = np.asarray(s_norm, dtype=float)

    z1 = float(Lambda) * (1.0 - float(mu)) * xi / np.sqrt(factors.tau1)
    theta1 = float(epsilon) * float(Lambda) ** 2 * (1.0 - float(mu)) * xi
    w_left = A1 * (np.cos(z1) - np.cosh(z1)) + B1 * (np.sin(z1) - np.sinh(z1))
    u_left = P1 * np.sin(theta1)

    z2 = -float(Lambda) * (1.0 + float(mu)) * xi / np.sqrt(factors.tau2)
    theta2 = -float(epsilon) * float(Lambda) ** 2 * (1.0 + float(mu)) * xi
    w_right = A2 * (np.cos(z2) - np.cosh(z2)) + B2 * (np.sin(z2) - np.sinh(z2))
    u_right = P2 * np.sin(theta2)

    if RIGHT_COORDINATE_FOR_PLOTTING == "joint-to-external":
        u_right = u_right[::-1]
        w_right = w_right[::-1]

    return {
        "u_left": u_left,
        "w_left": w_left,
        "u_right": u_right,
        "w_right": w_right,
    }


def shape_case(root: TrackedRoot, s_norm: np.ndarray) -> ShapeCase:
    matrix = assemble_clamped_coupled_matrix_eta(root.Lambda, np.deg2rad(BETA_DEG), root.mu, EPSILON, root.eta)
    coeff, smallest, ratio = analytic_null_vector(matrix)
    components_raw = reconstruct_components_eta(
        root.Lambda,
        mu=root.mu,
        eta=root.eta,
        epsilon=EPSILON,
        coeff=coeff,
        s_norm=s_norm,
    )
    components, _scale = normalize_components(components_raw, plot_kind="full", normalize=NORMALIZE)
    matrix_residual = matrix @ coeff
    return ShapeCase(
        eta=root.eta,
        mu=root.mu,
        Lambda=root.Lambda,
        current_sorted_index=root.current_sorted_index,
        tracking_step_status=root.tracking_step_status,
        components=components,
        matrix_residual_max=float(np.max(np.abs(matrix_residual))),
        smallest_singular_value=float(smallest),
        singular_value_ratio=float(ratio),
    )


def align_components_to_reference(
    components: dict[str, np.ndarray],
    reference: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    keys = ("u_left", "w_left", "u_right", "w_right")
    vector = np.concatenate([np.asarray(components[key], dtype=float) for key in keys])
    reference_vector = np.concatenate([np.asarray(reference[key], dtype=float) for key in keys])
    if float(np.dot(vector, reference_vector)) >= 0.0:
        return components
    return {key: -np.asarray(value, dtype=float) for key, value in components.items()}


def arm_lengths(mu: float) -> tuple[float, float]:
    ell = float(L_TOTAL) / 2.0
    return ell * (1.0 - float(mu)), ell * (1.0 + float(mu))


def base_coordinates(s_norm: np.ndarray, *, beta_rad: float, mu: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    left_length, right_length = arm_lengths(float(mu))
    xi = np.asarray(s_norm, dtype=float)
    x_left = left_length * xi
    y_left = np.zeros_like(x_left)
    x_right = left_length + right_length * xi * np.cos(beta_rad)
    y_right = right_length * xi * np.sin(beta_rad)
    return x_left, y_left, x_right, y_right


def deformed_coordinates(
    components: dict[str, np.ndarray],
    s_norm: np.ndarray,
    *,
    beta_rad: float,
    mu: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x_left, y_left, x_right, y_right = base_coordinates(s_norm, beta_rad=beta_rad, mu=mu)
    tangent_right = np.array([np.cos(beta_rad), np.sin(beta_rad)], dtype=float)
    normal_right = np.array([-np.sin(beta_rad), np.cos(beta_rad)], dtype=float)
    u_left = np.asarray(components["u_left"], dtype=float)
    w_left = np.asarray(components["w_left"], dtype=float)
    u_right = np.asarray(components["u_right"], dtype=float)
    w_right = np.asarray(components["w_right"], dtype=float)
    x_left_def = x_left + MODE_SCALE * u_left
    y_left_def = y_left + MODE_SCALE * w_left
    x_right_def = x_right + MODE_SCALE * (u_right * tangent_right[0] + w_right * normal_right[0])
    y_right_def = y_right + MODE_SCALE * (u_right * tangent_right[1] + w_right * normal_right[1])
    return x_left_def, y_left_def, x_right_def, y_right_def


def shared_axis_limits(cases: Sequence[ShapeCase], s_norm: np.ndarray, *, beta_rad: float, mu: float) -> tuple[tuple[float, float], tuple[float, float]]:
    x_left_base, y_left_base, x_right_base, y_right_base = base_coordinates(s_norm, beta_rad=beta_rad, mu=mu)
    x_values = [x_left_base, x_right_base]
    y_values = [y_left_base, y_right_base]
    for case in cases:
        x_left, y_left, x_right, y_right = deformed_coordinates(case.components, s_norm, beta_rad=beta_rad, mu=mu)
        x_values.extend([x_left, x_right])
        y_values.extend([y_left, y_right])
    x_all = np.concatenate([np.asarray(value, dtype=float) for value in x_values])
    y_all = np.concatenate([np.asarray(value, dtype=float) for value in y_values])
    x_span = max(float(np.max(x_all) - np.min(x_all)), 1.0)
    y_span = max(float(np.max(y_all) - np.min(y_all)), 0.5)
    return (
        (float(np.min(x_all) - 0.06 * x_span), float(np.max(x_all) + 0.06 * x_span)),
        (float(np.min(y_all) - 0.12 * y_span), float(np.max(y_all) + 0.12 * y_span)),
    )


def plot_cases_for_mu(mu: float, cases: Sequence[ShapeCase], s_norm: np.ndarray, output: Path) -> None:
    beta_rad = float(np.deg2rad(BETA_DEG))
    x_left_base, y_left_base, x_right_base, y_right_base = base_coordinates(s_norm, beta_rad=beta_rad, mu=mu)
    x_limits, y_limits = shared_axis_limits(cases, s_norm, beta_rad=beta_rad, mu=mu)
    colors = {-0.1: "#1f77b4", 0.0: "#ff7f0e", 0.1: "#2ca02c"}
    fallback_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig, axes = plt.subplots(1, len(cases), figsize=FIGSIZE, sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    for idx, (ax, case) in enumerate(zip(axes, cases, strict=True)):
        color = colors.get(round(float(case.eta), 10), fallback_colors[idx % len(fallback_colors)])
        ax.plot(x_left_base, y_left_base, color="0.72", linestyle="--", linewidth=1.0)
        ax.plot(x_right_base, y_right_base, color="0.72", linestyle="--", linewidth=1.0)
        x_left, y_left, x_right, y_right = deformed_coordinates(case.components, s_norm, beta_rad=beta_rad, mu=mu)
        ax.plot(x_left, y_left, color=color, linewidth=2.0)
        ax.plot(x_right, y_right, color=color, linewidth=2.0)
        ax.scatter([x_left_base[-1]], [y_left_base[-1]], color="black", s=12, zorder=5)
        ax.set_title(
            rf"$\eta={case.eta:g}$" + "\n"
            rf"$\Lambda={case.Lambda:.6g}$, sorted {case.current_sorted_index}",
            fontsize=10,
        )
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(*x_limits)
        ax.set_ylim(*y_limits)
        ax.grid(True, color="0.88", linewidth=0.6)
        ax.set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.suptitle(
        rf"Thickness mismatch branch {BRANCH_INDEX}: "
        rf"$\beta={BETA_DEG:g}^\circ$, $\epsilon={EPSILON:g}$, $\mu={mu:g}$",
        fontsize=12,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def build_shape_cases() -> dict[float, list[ShapeCase]]:
    s_norm = np.linspace(0.0, 1.0, NUM_SAMPLES)
    tracked_by_eta: dict[float, dict[float, TrackedRoot]] = {
        float(eta): track_branch_for_eta(float(eta), MU_VALUES) for eta in ETA_VALUES
    }
    cases_by_mu: dict[float, list[ShapeCase]] = {}
    for mu in MU_VALUES:
        cases: list[ShapeCase] = []
        reference_components: dict[str, np.ndarray] | None = None
        for eta in ETA_VALUES:
            case = shape_case(tracked_by_eta[float(eta)][float(mu)], s_norm)
            components = case.components
            if reference_components is None:
                reference_components = components
            else:
                components = align_components_to_reference(components, reference_components)
                case = ShapeCase(
                    eta=case.eta,
                    mu=case.mu,
                    Lambda=case.Lambda,
                    current_sorted_index=case.current_sorted_index,
                    tracking_step_status=case.tracking_step_status,
                    components=components,
                    matrix_residual_max=case.matrix_residual_max,
                    smallest_singular_value=case.smallest_singular_value,
                    singular_value_ratio=case.singular_value_ratio,
                )
            cases.append(case)
        cases_by_mu[float(mu)] = cases
    return cases_by_mu


def main() -> dict[float, list[ShapeCase]]:
    if BRANCH_INDEX <= 0:
        raise ValueError("BRANCH_INDEX must be positive.")
    if NUM_TRACKED_BRANCHES < BRANCH_INDEX:
        raise ValueError("NUM_TRACKED_BRANCHES must be at least BRANCH_INDEX.")
    if NUM_SORTED_ROOTS < NUM_TRACKED_BRANCHES:
        raise ValueError("NUM_SORTED_ROOTS must be at least NUM_TRACKED_BRANCHES.")

    s_norm = np.linspace(0.0, 1.0, NUM_SAMPLES)
    cases_by_mu = build_shape_cases()
    for mu in MU_VALUES:
        output = output_png_path(float(mu))
        plot_cases_for_mu(float(mu), cases_by_mu[float(mu)], s_norm, output)
        print(f"saved shape plot for mu={float(mu):g}: {output}")

    print("branch tracking and shape reconstruction summary:")
    print("mu, eta, branch_from_mu0, current_sorted_index, Lambda, tracking_status, max_matrix_residual")
    for mu in MU_VALUES:
        for case in cases_by_mu[float(mu)]:
            print(
                f"{float(mu):g}, {case.eta:g}, {BRANCH_INDEX}, "
                f"{case.current_sorted_index}, {case.Lambda:.10g}, "
                f"{case.tracking_step_status}, {case.matrix_residual_max:.6e}"
            )
    return cases_by_mu


if __name__ == "__main__":
    main()
