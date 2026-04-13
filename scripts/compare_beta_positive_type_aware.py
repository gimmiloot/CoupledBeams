from __future__ import annotations

import csv
from math import ceil
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.optimize import linear_sum_assignment


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.compare_beta0_analytic_vs_fem import (  # noqa: E402
    BASE_PARAMS,
    N_BENDING,
    N_FEM_TYPED,
    fem_axial_fractions,
    fem_parameter_override,
    first_axial_seed_from_mu0,
    relative_error,
)
from my_project.analytic.formulas import lambdas_to_frequencies  # noqa: E402
from my_project.analytic.solvers import find_first_n_roots, find_roots_scan_bisect  # noqa: E402
from my_project.fem import python_fem as fem  # noqa: E402


BETA_VALUES = (2.0, 5.0, 10.0, 15.0)
MU_VALUES = np.linspace(0.0, 0.9, 101)
BETA_CONTINUATION_STEP = 0.5
N_ANALYTIC_TRACK = 24
ANALYTIC_BENDING_CANDIDATES = 12
N_TRACKED_BRANCHES = N_BENDING + 1
BRANCH_ROWS = tuple(
    [{"analytic_branch_id": f"bending_{idx + 1}", "mode_type": "bending"} for idx in range(N_BENDING)]
    + [{"analytic_branch_id": "axial_1", "mode_type": "axial"}]
)
SEED_MODE_IDS_BETA0 = tuple(range(1, N_BENDING + 1)) + (16,)
FREQ_WEIGHT = 0.03

RESULTS_DIR = REPO_ROOT / "results"
MATCH_TABLE_PATH = RESULTS_DIR / "beta_positive_type_aware_matches.csv"
SUMMARY_TABLE_PATH = RESULTS_DIR / "beta_positive_type_aware_summary.csv"
BENDING_PLOT_PATH = RESULTS_DIR / "beta_positive_bending_vs_mu.png"
AXIAL_PLOT_PATH = RESULTS_DIR / "beta_positive_axial_vs_mu.png"
AXIAL_FRACTION_PLOT_PATH = RESULTS_DIR / "beta_positive_axial_fraction_vs_mu.png"


def solve_fem_modes_beta(
    mu: float,
    beta_deg: float,
    n_modes: int,
) -> tuple[np.ndarray, np.ndarray]:
    with fem_parameter_override(BASE_PARAMS):
        omega, vecs = fem.fem_solve(mu=mu, beta_deg=beta_deg, n_modes=n_modes)
        return omega * fem.scale, vecs


def beta_continuation_grid(beta_values: tuple[float, ...]) -> np.ndarray:
    beta_max = max(beta_values)
    grid = np.arange(0.0, beta_max + 0.5 * BETA_CONTINUATION_STEP, BETA_CONTINUATION_STEP)
    return np.asarray(grid, dtype=float)


def beta_key(beta_deg: float) -> float:
    return round(float(beta_deg), 6)


def write_csv_rows(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def mac_matrix(prev_vecs: np.ndarray, cur_vecs: np.ndarray) -> np.ndarray:
    out = np.zeros((prev_vecs.shape[1], cur_vecs.shape[1]), dtype=float)
    for row in range(prev_vecs.shape[1]):
        for col in range(cur_vecs.shape[1]):
            out[row, col] = fem.mac_value(prev_vecs[:, row], cur_vecs[:, col])
    return out


def assign_by_mac_and_frequency(
    prev_vecs: np.ndarray,
    prev_freqs: np.ndarray,
    cur_vecs: np.ndarray,
    cur_freqs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    mac = mac_matrix(prev_vecs, cur_vecs)
    freq_penalty = np.abs(prev_freqs[:, None] - cur_freqs[None, :]) / np.maximum(prev_freqs[:, None], 1.0)
    cost = (1.0 - mac) + FREQ_WEIGHT * freq_penalty
    rows, cols = linear_sum_assignment(cost)
    assignment = np.full(prev_vecs.shape[1], -1, dtype=int)
    assigned_mac = np.full(prev_vecs.shape[1], np.nan, dtype=float)
    for row, col in zip(rows, cols):
        assignment[row] = col
        assigned_mac[row] = mac[row, col]
    if np.any(assignment < 0):
        raise RuntimeError("Failed to assign all tracked branches.")
    return assignment, assigned_mac


def analytic_sorted_frequencies(beta_deg: float, mu: float, n_modes: int) -> np.ndarray:
    roots = find_first_n_roots(
        beta=np.deg2rad(beta_deg),
        mu=float(mu),
        eps=BASE_PARAMS.eps,
        n_roots=n_modes,
        Lmin=0.2,
        Lmax0=70.0,
        scan_step=0.01,
        grow_factor=1.35,
        max_tries=8,
    )
    return lambdas_to_frequencies(roots, BASE_PARAMS)


def analytic_tracked_frequencies(beta_deg: float, mu_values: np.ndarray, n_modes: int) -> np.ndarray:
    raise NotImplementedError


def assign_analytic_branches(
    prev_bending_mu0: np.ndarray,
    prev_axial_mu0: float,
    current_mu0: np.ndarray,
) -> tuple[np.ndarray, int]:
    candidate_count = min(ANALYTIC_BENDING_CANDIDATES, len(current_mu0))
    candidate_idx = np.arange(candidate_count, dtype=int)
    cost = np.abs(prev_bending_mu0[:, None] - current_mu0[candidate_idx][None, :])
    rows, cols = linear_sum_assignment(cost)
    bending_idx = np.full(len(prev_bending_mu0), -1, dtype=int)
    for row, col in zip(rows, cols):
        bending_idx[row] = candidate_idx[col]
    if np.any(bending_idx < 0):
        raise RuntimeError("Failed to assign all analytic bending branches.")

    axial_cost = np.abs(current_mu0 - prev_axial_mu0)
    axial_idx = int(np.argmin(axial_cost))
    return bending_idx, axial_idx


def track_fem_family(
    beta_values: tuple[float, ...],
    mu_values: np.ndarray,
) -> dict[float, dict[str, np.ndarray]]:
    seed_freqs_beta0, seed_vecs_beta0 = solve_fem_modes_beta(mu=0.0, beta_deg=0.0, n_modes=N_FEM_TYPED)
    axial_seed_index_beta0 = SEED_MODE_IDS_BETA0[-1] - 1

    prev_beta_vecs = seed_vecs_beta0[:, [axial_seed_index_beta0]]
    prev_beta_freqs = seed_freqs_beta0[[axial_seed_index_beta0]]
    target_keys = {beta_key(beta): beta for beta in beta_values}
    mu0_seeds: dict[float, dict[str, np.ndarray]] = {}

    for beta_deg in beta_continuation_grid(beta_values)[1:]:
        freqs_mu0, vecs_mu0 = solve_fem_modes_beta(mu=0.0, beta_deg=float(beta_deg), n_modes=N_FEM_TYPED)
        fracs_mu0 = fem_axial_fractions(vecs_mu0)
        assign_mu0, mac_mu0 = assign_by_mac_and_frequency(
            prev_vecs=prev_beta_vecs,
            prev_freqs=prev_beta_freqs,
            cur_vecs=vecs_mu0,
            cur_freqs=freqs_mu0,
        )
        current_vecs = vecs_mu0[:, assign_mu0]
        current_freqs = freqs_mu0[assign_mu0]
        current_ids = assign_mu0 + 1
        current_fracs = fracs_mu0[assign_mu0]

        key = beta_key(beta_deg)
        if key in target_keys:
            mu0_seeds[target_keys[key]] = {
                "vec": current_vecs[:, 0].copy(),
                "freq": float(current_freqs[0]),
                "mode_id": int(current_ids[0]),
                "axial_fraction": float(current_fracs[0]),
                "assignment_mac": float(mac_mu0[0]),
            }

        prev_beta_vecs = current_vecs
        prev_beta_freqs = current_freqs

    data: dict[float, dict[str, np.ndarray]] = {}

    for beta_deg in beta_values:
        branch_freqs = np.full((N_TRACKED_BRANCHES, len(mu_values)), np.nan, dtype=float)
        branch_ids = np.full((N_TRACKED_BRANCHES, len(mu_values)), -1, dtype=int)
        branch_fracs = np.full((N_TRACKED_BRANCHES, len(mu_values)), np.nan, dtype=float)
        branch_macs = np.full((N_TRACKED_BRANCHES, len(mu_values)), np.nan, dtype=float)

        seed = mu0_seeds[beta_deg]
        freqs_mu0, vecs_mu0 = solve_fem_modes_beta(mu=0.0, beta_deg=beta_deg, n_modes=N_FEM_TYPED)
        fracs_mu0 = fem_axial_fractions(vecs_mu0)
        branch_freqs[:N_BENDING, 0] = freqs_mu0[:N_BENDING]
        branch_ids[:N_BENDING, 0] = np.arange(1, N_BENDING + 1, dtype=int)
        branch_fracs[:N_BENDING, 0] = fracs_mu0[:N_BENDING]
        branch_macs[:N_BENDING, 0] = np.nan
        branch_freqs[-1, 0] = seed["freq"]
        branch_ids[-1, 0] = seed["mode_id"]
        branch_fracs[-1, 0] = seed["axial_fraction"]
        branch_macs[-1, 0] = seed["assignment_mac"]
        prev_axial_vec = np.asarray(seed["vec"], dtype=float)
        prev_axial_freq = float(seed["freq"])

        for col, mu in enumerate(mu_values[1:], start=1):
            cur_freqs, cur_vecs = solve_fem_modes_beta(mu=float(mu), beta_deg=beta_deg, n_modes=N_FEM_TYPED)
            cur_fracs = fem_axial_fractions(cur_vecs)
            branch_freqs[:N_BENDING, col] = cur_freqs[:N_BENDING]
            branch_ids[:N_BENDING, col] = np.arange(1, N_BENDING + 1, dtype=int)
            branch_fracs[:N_BENDING, col] = cur_fracs[:N_BENDING]
            branch_macs[:N_BENDING, col] = np.nan

            assign, mac_values = assign_by_mac_and_frequency(
                prev_vecs=prev_axial_vec[:, None],
                prev_freqs=np.asarray([prev_axial_freq], dtype=float),
                cur_vecs=cur_vecs,
                cur_freqs=cur_freqs,
            )
            axial_idx = int(assign[0])
            branch_freqs[-1, col] = cur_freqs[axial_idx]
            branch_ids[-1, col] = axial_idx + 1
            branch_fracs[-1, col] = cur_fracs[axial_idx]
            branch_macs[-1, col] = mac_values[0]
            prev_axial_vec = cur_vecs[:, axial_idx]
            prev_axial_freq = cur_freqs[axial_idx]

        data[beta_deg] = {
            "freqs": branch_freqs,
            "mode_ids": branch_ids,
            "axial_fractions": branch_fracs,
            "assignment_macs": branch_macs,
        }

    return data


def track_analytic_family(
    beta_values: tuple[float, ...],
    mu_values: np.ndarray,
) -> dict[float, dict[str, np.ndarray]]:
    beta0_axial_roots = find_first_n_roots(
        beta=0.0,
        mu=0.0,
        eps=BASE_PARAMS.eps,
        n_roots=SEED_MODE_IDS_BETA0[-1],
        Lmin=0.2,
        Lmax0=70.0,
        scan_step=0.01,
        grow_factor=1.35,
        max_tries=8,
    )
    prev_axial_lambda = float(beta0_axial_roots[SEED_MODE_IDS_BETA0[-1] - 1])
    target_keys = {beta_key(beta): beta for beta in beta_values}
    axial_seed_lambda: dict[float, float] = {}

    for beta_deg in beta_continuation_grid(beta_values)[1:]:
        roots = find_first_n_roots(
            beta=np.deg2rad(float(beta_deg)),
            mu=0.0,
            eps=BASE_PARAMS.eps,
            n_roots=N_ANALYTIC_TRACK,
            Lmin=0.2,
            Lmax0=70.0,
            scan_step=0.01,
            grow_factor=1.35,
            max_tries=8,
        )
        axial_idx = int(np.nanargmin(np.abs(roots - prev_axial_lambda)))
        prev_axial_lambda = float(roots[axial_idx])

        key = beta_key(beta_deg)
        if key in target_keys:
            axial_seed_lambda[target_keys[key]] = prev_axial_lambda

    data: dict[float, dict[str, np.ndarray]] = {}
    for beta_deg in beta_values:
        branch_freqs = np.full((N_TRACKED_BRANCHES, len(mu_values)), np.nan, dtype=float)
        branch_indices = np.full(N_TRACKED_BRANCHES, -1, dtype=int)
        lambda_prev = axial_seed_lambda[beta_deg]

        for col, mu in enumerate(mu_values):
            bending_freqs = analytic_sorted_frequencies(beta_deg=beta_deg, mu=float(mu), n_modes=N_BENDING)
            branch_freqs[:N_BENDING, col] = bending_freqs
            branch_indices[:N_BENDING] = np.arange(N_BENDING, dtype=int)

            roots_local: list[float] = []
            for window, step in ((1.5, 0.003), (2.5, 0.005), (4.0, 0.01)):
                roots_local = find_roots_scan_bisect(
                    beta=np.deg2rad(beta_deg),
                    mu=float(mu),
                    eps=BASE_PARAMS.eps,
                    n_roots=8,
                    Lmin=max(0.2, lambda_prev - window),
                    Lmax=lambda_prev + window,
                    scan_step=step,
                    bisect_iters=80,
                )
                if roots_local:
                    break
            if not roots_local:
                raise RuntimeError(f"Failed to continue analytic axial branch for beta={beta_deg}, mu={mu}")

            lambda_prev = float(min(roots_local, key=lambda x: abs(x - lambda_prev)))
            branch_freqs[-1, col] = lambdas_to_frequencies(np.asarray([lambda_prev], dtype=float), BASE_PARAMS)[0]
            branch_indices[-1] = 15

        data[beta_deg] = {"freqs": branch_freqs, "branch_indices": branch_indices + 1}

    return data


def build_match_rows(
    beta_values: tuple[float, ...],
    mu_values: np.ndarray,
    analytic_data: dict[float, dict[str, np.ndarray]],
    fem_data: dict[float, dict[str, np.ndarray]],
) -> tuple[list[dict[str, float | int | str]], list[dict[str, float | int | str]]]:
    match_rows: list[dict[str, float | int | str]] = []
    summary_rows: list[dict[str, float | int | str]] = []

    for beta_deg in beta_values:
        analytic_freqs = analytic_data[beta_deg]["freqs"]
        fem_freqs = fem_data[beta_deg]["freqs"]
        fem_ids = fem_data[beta_deg]["mode_ids"]
        fem_fracs = fem_data[beta_deg]["axial_fractions"]
        fem_macs = fem_data[beta_deg]["assignment_macs"]

        for branch_pos, branch_info in enumerate(BRANCH_ROWS):
            rel = np.divide(
                np.abs(analytic_freqs[branch_pos] - fem_freqs[branch_pos]),
                np.abs(fem_freqs[branch_pos]),
                out=np.full_like(fem_freqs[branch_pos], np.nan),
                where=np.abs(fem_freqs[branch_pos]) > 0,
            )
            for col, mu in enumerate(mu_values):
                match_rows.append(
                    {
                        "beta": beta_deg,
                        "mu": float(mu),
                        "analytic_branch_id": branch_info["analytic_branch_id"],
                        "fem_mode_id": int(fem_ids[branch_pos, col]),
                        "mode_type": branch_info["mode_type"],
                        "analytic_hz": float(analytic_freqs[branch_pos, col]),
                        "fem_hz": float(fem_freqs[branch_pos, col]),
                        "relative_error": float(rel[col]),
                        "fem_axial_fraction": float(fem_fracs[branch_pos, col]),
                    }
                )

            summary_rows.append(
                {
                    "beta": beta_deg,
                    "analytic_branch_id": branch_info["analytic_branch_id"],
                    "mode_type": branch_info["mode_type"],
                    "fem_mode_id_min": int(np.min(fem_ids[branch_pos])),
                    "fem_mode_id_max": int(np.max(fem_ids[branch_pos])),
                    "relative_error_max": float(np.nanmax(rel)),
                    "relative_error_mean": float(np.nanmean(rel)),
                    "fem_axial_fraction_mu0": float(fem_fracs[branch_pos, 0]),
                    "fem_axial_fraction_min": float(np.nanmin(fem_fracs[branch_pos])),
                    "fem_axial_fraction_max": float(np.nanmax(fem_fracs[branch_pos])),
                    "assignment_mac_min": (
                        float(np.nanmin(fem_macs[branch_pos]))
                        if np.any(np.isfinite(fem_macs[branch_pos]))
                        else np.nan
                    ),
                    "assignment_mac_mean": (
                        float(np.nanmean(fem_macs[branch_pos]))
                        if np.any(np.isfinite(fem_macs[branch_pos]))
                        else np.nan
                    ),
                }
            )

    return match_rows, summary_rows


def plot_bending_grid(
    beta_values: tuple[float, ...],
    mu_values: np.ndarray,
    analytic_data: dict[float, dict[str, np.ndarray]],
    fem_data: dict[float, dict[str, np.ndarray]],
) -> None:
    ncols = 2
    nrows = ceil(len(beta_values) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(13.5, 4.8 * nrows), squeeze=False)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for idx, beta_deg in enumerate(beta_values):
        ax = axes[idx // ncols][idx % ncols]
        analytic_freqs = analytic_data[beta_deg]["freqs"][:N_BENDING]
        fem_freqs = fem_data[beta_deg]["freqs"][:N_BENDING]

        for branch_idx in range(N_BENDING):
            color = colors[branch_idx % len(colors)]
            ax.plot(mu_values, analytic_freqs[branch_idx], color=color, linewidth=2.0)
            ax.plot(
                mu_values,
                fem_freqs[branch_idx],
                linestyle="None",
                marker="o",
                markersize=3.2,
                color=color,
                markerfacecolor=color,
                markeredgecolor=color,
            )

        ax.set_title(f"beta = {beta_deg:.1f} deg")
        ax.set_xlabel("mu")
        ax.set_ylabel("Frequency (Hz)")
        ax.grid(True, alpha=0.3)

    style_handles = [
        Line2D([0], [0], color="black", lw=2.0, label="analytic"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=5.0, label="FEM"),
    ]
    mode_handles = [
        Line2D([0], [0], color=colors[i % len(colors)], lw=2.0, label=f"bending_{i + 1}")
        for i in range(N_BENDING)
    ]

    for idx in range(len(beta_values), nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle("Type-aware comparison for beta > 0: bending descendants vs mu", y=0.995)
    fig.legend(handles=style_handles + mode_handles, ncols=4, loc="upper center", bbox_to_anchor=(0.5, 0.98))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(BENDING_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_axial_grid(
    beta_values: tuple[float, ...],
    mu_values: np.ndarray,
    analytic_data: dict[float, dict[str, np.ndarray]],
    fem_data: dict[float, dict[str, np.ndarray]],
) -> None:
    ncols = 2
    nrows = ceil(len(beta_values) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(13.0, 4.6 * nrows), squeeze=False)

    for idx, beta_deg in enumerate(beta_values):
        ax = axes[idx // ncols][idx % ncols]
        analytic_axial = analytic_data[beta_deg]["freqs"][-1]
        fem_axial = fem_data[beta_deg]["freqs"][-1]
        ax.plot(mu_values, analytic_axial, linewidth=2.0, label="analytic axial_1")
        ax.plot(mu_values, fem_axial, linestyle="None", marker="o", markersize=3.5, label="FEM axial_1")
        ax.set_title(f"beta = {beta_deg:.1f} deg")
        ax.set_xlabel("mu")
        ax.set_ylabel("Frequency (Hz)")
        ax.grid(True, alpha=0.3)
        ax.legend()

    for idx in range(len(beta_values), nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle("Type-aware comparison for beta > 0: first axial descendant vs mu", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(AXIAL_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_axial_fraction_grid(
    beta_values: tuple[float, ...],
    mu_values: np.ndarray,
    fem_data: dict[float, dict[str, np.ndarray]],
) -> None:
    ncols = 2
    nrows = ceil(len(beta_values) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(13.0, 4.6 * nrows), squeeze=False)

    for idx, beta_deg in enumerate(beta_values):
        ax = axes[idx // ncols][idx % ncols]
        fem_axial_fraction = fem_data[beta_deg]["axial_fractions"][-1]
        ax.plot(mu_values, fem_axial_fraction, linewidth=2.0)
        ax.axhline(0.9, linestyle="--", linewidth=1.0, alpha=0.5)
        ax.axhline(0.5, linestyle=":", linewidth=1.0, alpha=0.5)
        ax.set_title(f"beta = {beta_deg:.1f} deg")
        ax.set_xlabel("mu")
        ax.set_ylabel("FEM axial fraction")
        ax.set_ylim(0.0, 1.0)
        ax.grid(True, alpha=0.3)

    for idx in range(len(beta_values), nrows * ncols):
        axes[idx // ncols][idx % ncols].axis("off")

    fig.suptitle("Mixing diagnostic for the tracked axial descendant", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(AXIAL_FRACTION_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)

    analytic_data = track_analytic_family(BETA_VALUES, MU_VALUES)
    fem_data = track_fem_family(BETA_VALUES, MU_VALUES)
    match_rows, summary_rows = build_match_rows(
        beta_values=BETA_VALUES,
        mu_values=MU_VALUES,
        analytic_data=analytic_data,
        fem_data=fem_data,
    )

    write_csv_rows(MATCH_TABLE_PATH, match_rows)
    write_csv_rows(SUMMARY_TABLE_PATH, summary_rows)
    plot_bending_grid(BETA_VALUES, MU_VALUES, analytic_data, fem_data)
    plot_axial_grid(BETA_VALUES, MU_VALUES, analytic_data, fem_data)
    plot_axial_fraction_grid(BETA_VALUES, MU_VALUES, fem_data)

    print(f"betas analysed: {BETA_VALUES}")
    print(f"mu grid: {MU_VALUES[0]:.6f} .. {MU_VALUES[-1]:.6f} ({len(MU_VALUES)} points)")
    print(f"saved table: {MATCH_TABLE_PATH}")
    print(f"saved summary: {SUMMARY_TABLE_PATH}")
    print(f"saved plot: {BENDING_PLOT_PATH}")
    print(f"saved plot: {AXIAL_PLOT_PATH}")
    print(f"saved plot: {AXIAL_FRACTION_PLOT_PATH}")

    for beta_deg in BETA_VALUES:
        branch_rows = [row for row in summary_rows if float(row["beta"]) == beta_deg]
        axial_row = next(row for row in branch_rows if row["analytic_branch_id"] == "axial_1")
        bending_rows = [row for row in branch_rows if str(row["mode_type"]) == "bending"]
        print(
            f"beta={beta_deg:.1f} deg: "
            f"axial mode id range {int(axial_row['fem_mode_id_min'])}-{int(axial_row['fem_mode_id_max'])}, "
            f"axial fraction mu0={float(axial_row['fem_axial_fraction_mu0']):.3f}, "
            f"axial fraction min={float(axial_row['fem_axial_fraction_min']):.3f}, "
            f"axial rel err max={float(axial_row['relative_error_max']):.3e}, "
            f"axial min MAC={float(axial_row['assignment_mac_min']):.3f}, "
            f"bending rel err max={max(float(row['relative_error_max']) for row in bending_rows):.3e}"
        )


if __name__ == "__main__":
    main()
