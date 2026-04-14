from __future__ import annotations

from contextlib import contextmanager
import csv
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.FreqFromMu import BeamParams  # noqa: E402
from my_project.analytic.formulas import lambdas_to_frequencies  # noqa: E402
from my_project.analytic.solvers import find_first_n_roots, tracked_lambdas_vs_mu  # noqa: E402
from my_project.fem import python_fem as fem  # noqa: E402
from scripts.sweep_grid_policy import ANALYSIS_MU_STEP, analysis_mu_grid  # noqa: E402


BASE_PARAMS = BeamParams(E=2.1e11, rho=7800.0, r=0.005, L_total=2.0)
BETA_DEG = 0.0
N_BENDING = 6
N_CLASSIFY = 20
N_FEM_TYPED = 30
RADIUS_LOW_VALUES = (0.0045, 0.0050, 0.0055)
RADIUS_AXIAL_SWEEP = np.linspace(0.0040, 0.0062, 23)

RESULTS_DIR = REPO_ROOT / "results"
FEM_CSV_PATH = RESULTS_DIR / "fem_spectrum.csv"
LEGACY_BETA0_PLOT_PATH = RESULTS_DIR / "beta0_analytic_vs_fem.png"
BENDING_PLOT_PATH = RESULTS_DIR / "beta0_bending_modes_vs_mu.png"
AXIAL_MU_PLOT_PATH = RESULTS_DIR / "beta0_axial_mode_vs_mu.png"
AXIAL_RADIUS_PLOT_PATH = RESULTS_DIR / "beta0_axial_mode_vs_radius.png"
MU0_TABLE_PATH = RESULTS_DIR / "beta0_mu0_comparison.csv"
MODE_TYPES_PATH = RESULTS_DIR / "beta0_mu0_mode_types.csv"
RADIUS_LOW_PATH = RESULTS_DIR / "beta0_radius_low_modes.csv"
RADIUS_AXIAL_PATH = RESULTS_DIR / "beta0_radius_axial_summary.csv"
TYPE_MATCH_MU_PATH = RESULTS_DIR / "beta0_type_aware_matches_vs_mu.csv"
TYPE_MATCH_RADIUS_PATH = RESULTS_DIR / "beta0_type_aware_matches_vs_radius.csv"
AXIAL_RADIUS_SWEEP_PATH = RESULTS_DIR / "beta0_axial_mode_vs_radius.csv"

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


def build_params(radius: float) -> BeamParams:
    return BeamParams(E=BASE_PARAMS.E, rho=BASE_PARAMS.rho, r=radius, L_total=BASE_PARAMS.L_total)


def relative_error(analytic_hz: float, fem_hz: float) -> float:
    return abs(analytic_hz - fem_hz) / abs(fem_hz) if fem_hz else np.nan


def analytic_sorted_frequencies(
    params: BeamParams,
    mu: float,
    n_modes: int,
    Lmax0: float,
    scan_step: float,
) -> np.ndarray:
    roots = find_first_n_roots(
        beta=np.deg2rad(BETA_DEG),
        mu=float(mu),
        eps=params.eps,
        n_roots=n_modes,
        Lmin=0.2,
        Lmax0=Lmax0,
        scan_step=scan_step,
        grow_factor=1.35,
        max_tries=8,
    )
    return lambdas_to_frequencies(roots, params)


def analytic_mu0_frequencies(params: BeamParams, n_modes: int) -> np.ndarray:
    return analytic_sorted_frequencies(params=params, mu=0.0, n_modes=n_modes, Lmax0=60.0, scan_step=0.01)


def analytic_bending_mu_sweep(params: BeamParams, mu_values: np.ndarray, n_modes: int) -> np.ndarray:
    lambdas_tr = tracked_lambdas_vs_mu(
        params=params,
        beta_deg=BETA_DEG,
        mu_values=np.asarray(mu_values, dtype=float),
        n_modes=n_modes,
        Lmin=0.2,
        Lmax0=55.0,
        scan_step=0.02,
        grow_factor=1.35,
        max_tries=8,
        tracking_method="auto",
    )
    return lambdas_to_frequencies(lambdas_tr, params)


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


def solve_fem_modes(params: BeamParams, mu: float, n_modes: int) -> tuple[np.ndarray, np.ndarray]:
    with fem_parameter_override(params):
        omega, vecs = fem.fem_solve(mu=mu, beta_deg=BETA_DEG, n_modes=n_modes)
        return omega * fem.scale, vecs


def fem_free_dof_axial_mask() -> np.ndarray:
    n1 = n2 = fem.N_ELEM
    ndof = 3 * (n1 + n2 + 1)
    bc = {0, 1, 2, 3 * (n1 + n2), 3 * (n1 + n2) + 1, 3 * (n1 + n2) + 2}
    free = np.array(sorted(set(range(ndof)) - bc), dtype=int)
    return (free % 3) == 0


def fem_axial_fractions(vecs: np.ndarray) -> np.ndarray:
    axial_mask = fem_free_dof_axial_mask()
    fractions = []
    for idx in range(vecs.shape[1]):
        phi = vecs[:, idx]
        total = float(np.dot(phi, phi))
        axial = float(np.dot(phi[axial_mask], phi[axial_mask]))
        fractions.append(axial / total if total > 0 else np.nan)
    return np.asarray(fractions, dtype=float)


def mode_kind(axial_fraction: float) -> str:
    if axial_fraction > 0.9:
        return "axial"
    if axial_fraction < 0.1:
        return "bending"
    return "mixed"


def fem_mode_rows(params: BeamParams, mu: float, n_modes: int) -> list[dict[str, float | int | str]]:
    freqs, vecs = solve_fem_modes(params, mu=mu, n_modes=n_modes)
    fractions = fem_axial_fractions(vecs)

    rows = []
    for idx, (freq, fraction) in enumerate(zip(freqs, fractions), start=1):
        rows.append(
            {
                "mode_id": idx,
                "frequency_hz": float(freq),
                "axial_fraction": float(fraction),
                "mode_type": mode_kind(float(fraction)),
            }
        )
    return rows


def select_mode_rows(
    mode_rows: list[dict[str, float | int | str]],
    mode_type: str,
    count: int,
) -> list[dict[str, float | int | str]]:
    selected = [row for row in mode_rows if row["mode_type"] == mode_type]
    if len(selected) < count:
        raise RuntimeError(f"Needed {count} {mode_type} modes, got {len(selected)}")
    return selected[:count]


def first_axial_seed_from_mu0(
    params: BeamParams,
    mode_rows: list[dict[str, float | int | str]] | None = None,
    analytic_hz_override: float | None = None,
) -> dict[str, float | int | str]:
    rows = mode_rows if mode_rows is not None else fem_mode_rows(params, mu=0.0, n_modes=N_FEM_TYPED)
    axial_row = select_mode_rows(rows, mode_type="axial", count=1)[0]
    fem_mode_id = int(axial_row["mode_id"])
    analytic_hz = analytic_hz_override
    if analytic_hz is None:
        analytic_freqs = analytic_mu0_frequencies(params, fem_mode_id)
        analytic_hz = float(analytic_freqs[fem_mode_id - 1])
    return {
        "analytic_branch_id": "axial_1",
        "analytic_hz": float(analytic_hz),
        "fem_mode_id": fem_mode_id,
        "fem_hz": float(axial_row["frequency_hz"]),
        "fem_axial_fraction": float(axial_row["axial_fraction"]),
        "mode_type": "axial",
    }


def load_fem_frequency_table(csv_path: Path, n_modes: int) -> tuple[np.ndarray, np.ndarray]:
    data = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=float, encoding=None)
    required = ["mu"] + [f"f{i + 1}_Hz_beta0" for i in range(n_modes)]
    missing = [name for name in required if name not in data.dtype.names]
    if missing:
        raise KeyError(f"Missing columns in {csv_path}: {missing}")

    mu_values = np.asarray(data["mu"], dtype=float)
    if mu_values.ndim != 1 or len(mu_values) == 0:
        raise ValueError(f"Bad mu column in {csv_path}")
    if not np.all(np.diff(mu_values) >= 0):
        raise ValueError(f"mu column is not monotone in {csv_path}")

    fem_freqs = np.vstack([np.asarray(data[f"f{i + 1}_Hz_beta0"], dtype=float) for i in range(n_modes)])
    return mu_values, fem_freqs


def compute_baseline_fem_sweep(mu_values: np.ndarray, n_modes: int) -> np.ndarray:
    with fem_parameter_override(BASE_PARAMS):
        tracked = fem.track_modes(np.asarray(mu_values, dtype=float), beta_deg=BETA_DEG, n_track=n_modes)
        return tracked.T * fem.scale


def load_or_compute_fem_sweep(mu_values: np.ndarray | None, n_modes: int) -> tuple[np.ndarray, np.ndarray, str]:
    requested_mu = np.asarray(analysis_mu_grid(), dtype=float) if mu_values is None else np.asarray(mu_values, dtype=float)

    if FEM_CSV_PATH.exists():
        try:
            csv_mu, fem_freqs = load_fem_frequency_table(FEM_CSV_PATH, n_modes=n_modes)
            fem_mu0, _ = solve_fem_modes(BASE_PARAMS, mu=0.0, n_modes=n_modes)
            if abs(csv_mu[0]) > 1e-12:
                raise ValueError("existing FEM CSV does not start at mu=0")
            if len(csv_mu) != len(requested_mu) or not np.allclose(csv_mu, requested_mu, atol=1e-10, rtol=0.0):
                raise ValueError("existing FEM CSV uses a different mu grid than the current analysis policy")
            if not np.allclose(fem_freqs[:, 0], fem_mu0, rtol=1e-4, atol=1e-3):
                raise ValueError("existing FEM CSV does not match the current baseline FEM at mu=0")
            return csv_mu, fem_freqs, str(FEM_CSV_PATH.relative_to(REPO_ROOT))
        except Exception as exc:
            print(f"[info] Existing FEM CSV is unsuitable for beta=0 comparison: {exc}")

    fem_freqs = compute_baseline_fem_sweep(requested_mu, n_modes=n_modes)
    return requested_mu, fem_freqs, "recomputed in-memory from src/my_project/fem/python_fem.py"


def build_mu0_rows() -> list[dict[str, float | int | str]]:
    analytic = analytic_mu0_frequencies(BASE_PARAMS, N_BENDING)
    fem_freqs, _ = solve_fem_modes(BASE_PARAMS, mu=0.0, n_modes=N_BENDING)

    rows = []
    for idx, (ana, fem_freq) in enumerate(zip(analytic, fem_freqs), start=1):
        abs_err = abs(ana - fem_freq)
        rows.append(
            {
                "mode": idx,
                "identifier": f"sorted mode {idx}",
                "analytic_hz": ana,
                "fem_hz": fem_freq,
                "abs_error_hz": abs_err,
                "rel_error": relative_error(ana, fem_freq),
            }
        )
    return rows


def build_mode_type_rows() -> list[dict[str, float | int | str]]:
    analytic = analytic_mu0_frequencies(BASE_PARAMS, N_CLASSIFY)
    fem_freqs, vecs = solve_fem_modes(BASE_PARAMS, mu=0.0, n_modes=N_CLASSIFY)
    fractions = fem_axial_fractions(vecs)

    rows = []
    for idx, (ana, fem_freq, fraction) in enumerate(zip(analytic, fem_freqs, fractions), start=1):
        rows.append(
            {
                "mode": idx,
                "analytic_hz": ana,
                "fem_hz": fem_freq,
                "axial_fraction": float(fraction),
                "kind": mode_kind(float(fraction)),
            }
        )
    return rows


def collect_mu_type_aware_matches(mu_values: np.ndarray, analytic_axial_hz: float) -> dict[str, object]:
    analytic_bending = analytic_bending_mu_sweep(BASE_PARAMS, mu_values=mu_values, n_modes=N_BENDING)
    axial_seed = first_axial_seed_from_mu0(BASE_PARAMS, analytic_hz_override=analytic_axial_hz)
    analytic_axial = np.full(len(mu_values), float(analytic_axial_hz), dtype=float)

    fem_bending = np.full((N_BENDING, len(mu_values)), np.nan, dtype=float)
    fem_bending_mode_ids = np.full((N_BENDING, len(mu_values)), -1, dtype=int)
    fem_axial = np.full(len(mu_values), np.nan, dtype=float)
    fem_axial_mode_ids = np.full(len(mu_values), -1, dtype=int)
    rows: list[dict[str, float | int | str]] = []

    for col, mu in enumerate(mu_values):
        mode_rows = fem_mode_rows(BASE_PARAMS, mu=float(mu), n_modes=N_FEM_TYPED)
        bending_rows = select_mode_rows(mode_rows, mode_type="bending", count=N_BENDING)
        axial_row = select_mode_rows(mode_rows, mode_type="axial", count=1)[0]

        for branch_idx, fem_row in enumerate(bending_rows):
            analytic_hz = float(analytic_bending[branch_idx, col])
            fem_hz = float(fem_row["frequency_hz"])
            fem_bending[branch_idx, col] = fem_hz
            fem_bending_mode_ids[branch_idx, col] = int(fem_row["mode_id"])
            rows.append(
                {
                    "parameter_kind": "mu",
                    "mu": float(mu),
                    "radius_mm": BASE_PARAMS.r * 1e3,
                    "analytic_branch_id": f"bending_{branch_idx + 1}",
                    "fem_mode_id": int(fem_row["mode_id"]),
                    "mode_type": "bending",
                    "analytic_hz": analytic_hz,
                    "fem_hz": fem_hz,
                    "relative_error": relative_error(analytic_hz, fem_hz),
                    "fem_axial_fraction": float(fem_row["axial_fraction"]),
                }
            )

        fem_axial[col] = float(axial_row["frequency_hz"])
        fem_axial_mode_ids[col] = int(axial_row["mode_id"])
        rows.append(
            {
                "parameter_kind": "mu",
                "mu": float(mu),
                "radius_mm": BASE_PARAMS.r * 1e3,
                "analytic_branch_id": "axial_1",
                "fem_mode_id": int(axial_row["mode_id"]),
                "mode_type": "axial",
                "analytic_hz": float(analytic_axial[col]),
                "fem_hz": float(axial_row["frequency_hz"]),
                "relative_error": relative_error(float(analytic_axial[col]), float(axial_row["frequency_hz"])),
                "fem_axial_fraction": float(axial_row["axial_fraction"]),
            }
        )

    return {
        "analytic_bending": analytic_bending,
        "analytic_axial": analytic_axial,
        "fem_bending": fem_bending,
        "fem_bending_mode_ids": fem_bending_mode_ids,
        "fem_axial": fem_axial,
        "fem_axial_mode_ids": fem_axial_mode_ids,
        "rows": rows,
        "axial_seed": axial_seed,
    }


def collect_radius_low_matches(radius_values: tuple[float, ...], analytic_axial_hz: float) -> dict[str, object]:
    rows: list[dict[str, float | int | str]] = []
    low_rows: list[dict[str, float | int | str]] = []
    axial_rows: list[dict[str, float | int | str]] = []

    for radius in radius_values:
        params = build_params(radius)
        analytic_bending = analytic_mu0_frequencies(params, N_BENDING)
        mode_rows = fem_mode_rows(params, mu=0.0, n_modes=N_FEM_TYPED)
        bending_rows = select_mode_rows(mode_rows, mode_type="bending", count=N_BENDING)
        axial_seed = first_axial_seed_from_mu0(
            params,
            mode_rows=mode_rows,
            analytic_hz_override=analytic_axial_hz,
        )

        for branch_idx, fem_row in enumerate(bending_rows):
            analytic_hz = float(analytic_bending[branch_idx])
            fem_hz = float(fem_row["frequency_hz"])
            rows.append(
                {
                    "parameter_kind": "radius_low",
                    "mu": 0.0,
                    "radius_mm": radius * 1e3,
                    "analytic_branch_id": f"bending_{branch_idx + 1}",
                    "fem_mode_id": int(fem_row["mode_id"]),
                    "mode_type": "bending",
                    "analytic_hz": analytic_hz,
                    "fem_hz": fem_hz,
                    "relative_error": relative_error(analytic_hz, fem_hz),
                    "fem_axial_fraction": float(fem_row["axial_fraction"]),
                }
            )
            low_rows.append(
                {
                    "radius_mm": radius * 1e3,
                    "mode": branch_idx + 1,
                    "kind": "bending",
                    "analytic_hz": analytic_hz,
                    "fem_hz": fem_hz,
                }
            )

        rows.append(
            {
                "parameter_kind": "radius_low",
                "mu": 0.0,
                "radius_mm": radius * 1e3,
                "analytic_branch_id": "axial_1",
                "fem_mode_id": int(axial_seed["fem_mode_id"]),
                "mode_type": "axial",
                "analytic_hz": float(axial_seed["analytic_hz"]),
                "fem_hz": float(axial_seed["fem_hz"]),
                "relative_error": relative_error(float(axial_seed["analytic_hz"]), float(axial_seed["fem_hz"])),
                "fem_axial_fraction": float(axial_seed["fem_axial_fraction"]),
            }
        )
        axial_rows.append(
            {
                "radius_mm": radius * 1e3,
                "first_axial_sorted_mode": int(axial_seed["fem_mode_id"]),
                "analytic_hz": float(axial_seed["analytic_hz"]),
                "fem_hz": float(axial_seed["fem_hz"]),
                "axial_fraction": float(axial_seed["fem_axial_fraction"]),
            }
        )

    return {"rows": rows, "low_rows": low_rows, "axial_rows": axial_rows}


def collect_axial_radius_sweep(radius_values: np.ndarray, analytic_axial_hz: float) -> dict[str, object]:
    radii_mm = []
    analytic_axial = []
    fem_axial = []
    fem_mode_ids = []
    rows: list[dict[str, float | int | str]] = []

    for radius in radius_values:
        params = build_params(float(radius))
        mode_rows = fem_mode_rows(params, mu=0.0, n_modes=N_FEM_TYPED)
        axial_seed = first_axial_seed_from_mu0(
            params,
            mode_rows=mode_rows,
            analytic_hz_override=analytic_axial_hz,
        )
        radius_mm = float(radius) * 1e3

        radii_mm.append(radius_mm)
        analytic_axial.append(float(axial_seed["analytic_hz"]))
        fem_axial.append(float(axial_seed["fem_hz"]))
        fem_mode_ids.append(int(axial_seed["fem_mode_id"]))
        rows.append(
            {
                "parameter_kind": "radius_axial",
                "mu": 0.0,
                "radius_mm": radius_mm,
                "analytic_branch_id": "axial_1",
                "fem_mode_id": int(axial_seed["fem_mode_id"]),
                "mode_type": "axial",
                "analytic_hz": float(axial_seed["analytic_hz"]),
                "fem_hz": float(axial_seed["fem_hz"]),
                "relative_error": relative_error(float(axial_seed["analytic_hz"]), float(axial_seed["fem_hz"])),
                "fem_axial_fraction": float(axial_seed["fem_axial_fraction"]),
            }
        )

    return {
        "radius_mm": np.asarray(radii_mm, dtype=float),
        "analytic_axial": np.asarray(analytic_axial, dtype=float),
        "fem_axial": np.asarray(fem_axial, dtype=float),
        "fem_mode_ids": np.asarray(fem_mode_ids, dtype=int),
        "rows": rows,
    }


def write_csv_rows(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot_bending_vs_mu(mu_values: np.ndarray, analytic_bending: np.ndarray, fem_bending: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(11.0, 6.4))

    mode_handles = []
    for idx in range(analytic_bending.shape[0]):
        (line,) = ax.plot(mu_values, analytic_bending[idx], linewidth=2.0)
        color = line.get_color()
        ax.plot(
            mu_values,
            fem_bending[idx],
            linestyle="None",
            marker="o",
            markersize=3.6,
            color=color,
            markerfacecolor=color,
            markeredgecolor=color,
        )
        mode_handles.append(Line2D([0], [0], color=color, lw=2.0, label=f"bending_{idx + 1}"))

    style_handles = [
        Line2D([0], [0], color="black", lw=2.0, label="analytic"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=5.0, label="FEM"),
    ]

    ax.set_xlabel("mu")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title("beta = 0: first bending branches only\nsolid lines = analytic, markers = FEM, type-aware matching")
    ax.grid(True, alpha=0.3)
    ax.legend(handles=style_handles + mode_handles, ncols=2, fontsize=9)
    fig.tight_layout()
    fig.savefig(BENDING_PLOT_PATH, dpi=200, bbox_inches="tight")
    fig.savefig(LEGACY_BETA0_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_axial_vs_mu(mu_values: np.ndarray, analytic_axial: np.ndarray, fem_axial: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(10.0, 5.2))
    ax.plot(mu_values, analytic_axial, linewidth=2.0, label="analytic axial_1")
    ax.plot(mu_values, fem_axial, linestyle="None", marker="o", markersize=3.8, label="FEM axial_1")
    ax.set_xlabel("mu")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title("beta = 0: first axial branch vs mu")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(AXIAL_MU_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_axial_vs_radius(radius_mm: np.ndarray, analytic_axial: np.ndarray, fem_axial: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(10.0, 5.2))
    ax.plot(radius_mm, analytic_axial, linewidth=2.0, label="analytic axial_1")
    ax.plot(radius_mm, fem_axial, linestyle="None", marker="o", markersize=4.0, label="FEM axial_1")
    ax.set_xlabel("Radius (mm)")
    ax.set_ylabel("Frequency (Hz)")
    ax.set_title("beta = 0, mu = 0: first axial branch vs radius")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(AXIAL_RADIUS_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def print_mu0_table(rows: list[dict[str, float | int | str]]) -> None:
    print("\nmu = 0 comparison (matched by ascending frequency order):")
    print("mode  analytic_hz   fem_hz        abs_error_hz  rel_error")
    for row in rows:
        print(
            f"{int(row['mode']):>4d}  "
            f"{float(row['analytic_hz']):>11.6f}  "
            f"{float(row['fem_hz']):>11.6f}  "
            f"{float(row['abs_error_hz']):>13.6e}  "
            f"{float(row['rel_error']):>9.3e}"
        )


def print_radius_summary(
    low_rows: list[dict[str, float | int | str]],
    axial_rows: list[dict[str, float | int | str]],
) -> None:
    print("\nRadius sensitivity at mu = 0 (first 6 bending modes, type-aware):")
    print("radius_mm  mode  kind      analytic_hz   fem_hz")
    for row in low_rows:
        print(
            f"{float(row['radius_mm']):>9.3f}  "
            f"{int(row['mode']):>4d}  "
            f"{str(row['kind']):>7s}  "
            f"{float(row['analytic_hz']):>11.6f}  "
            f"{float(row['fem_hz']):>11.6f}"
        )

    print("\nFirst axial branch at mu = 0:")
    print("radius_mm  sorted_mode  analytic_hz   fem_hz       axial_fraction")
    for row in axial_rows:
        print(
            f"{float(row['radius_mm']):>9.3f}  "
            f"{int(row['first_axial_sorted_mode']):>11d}  "
            f"{float(row['analytic_hz']):>11.6f}  "
            f"{float(row['fem_hz']):>11.6f}  "
            f"{float(row['axial_fraction']):>15.6f}"
        )


def constant_id_ranges(parameter_values: np.ndarray, mode_ids: np.ndarray) -> list[dict[str, float | int]]:
    ranges = []
    start = 0
    for idx in range(1, len(mode_ids) + 1):
        if idx == len(mode_ids) or mode_ids[idx] != mode_ids[start]:
            ranges.append(
                {
                    "mode_id": int(mode_ids[start]),
                    "param_start": float(parameter_values[start]),
                    "param_end": float(parameter_values[idx - 1]),
                }
            )
            start = idx
    return ranges


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)

    mu_values, fem_low_from_csv, fem_source = load_or_compute_fem_sweep(mu_values=None, n_modes=N_BENDING)
    mu0_rows = build_mu0_rows()
    mode_type_rows = build_mode_type_rows()
    analytic_axial_hz = float(first_axial_seed_from_mu0(BASE_PARAMS)["analytic_hz"])
    mu_matches = collect_mu_type_aware_matches(mu_values, analytic_axial_hz=analytic_axial_hz)
    radius_low_matches = collect_radius_low_matches(RADIUS_LOW_VALUES, analytic_axial_hz=analytic_axial_hz)
    axial_radius_sweep = collect_axial_radius_sweep(RADIUS_AXIAL_SWEEP, analytic_axial_hz=analytic_axial_hz)

    write_csv_rows(MU0_TABLE_PATH, mu0_rows)
    write_csv_rows(MODE_TYPES_PATH, mode_type_rows)
    write_csv_rows(RADIUS_LOW_PATH, radius_low_matches["low_rows"])
    write_csv_rows(RADIUS_AXIAL_PATH, radius_low_matches["axial_rows"])
    write_csv_rows(TYPE_MATCH_MU_PATH, mu_matches["rows"])
    write_csv_rows(TYPE_MATCH_RADIUS_PATH, radius_low_matches["rows"])
    write_csv_rows(AXIAL_RADIUS_SWEEP_PATH, axial_radius_sweep["rows"])

    plot_bending_vs_mu(
        mu_values=mu_values,
        analytic_bending=mu_matches["analytic_bending"],
        fem_bending=mu_matches["fem_bending"],
    )
    plot_axial_vs_mu(
        mu_values=mu_values,
        analytic_axial=mu_matches["analytic_axial"],
        fem_axial=mu_matches["fem_axial"],
    )
    plot_axial_vs_radius(
        radius_mm=axial_radius_sweep["radius_mm"],
        analytic_axial=axial_radius_sweep["analytic_axial"],
        fem_axial=axial_radius_sweep["fem_axial"],
    )

    csv_vs_typed = np.abs(fem_low_from_csv - mu_matches["fem_bending"])
    bending_abs_err = np.abs(mu_matches["analytic_bending"] - mu_matches["fem_bending"])
    bending_rel_err = np.divide(
        bending_abs_err,
        np.abs(mu_matches["fem_bending"]),
        out=np.full_like(bending_abs_err, np.nan),
        where=np.abs(mu_matches["fem_bending"]) > 0,
    )
    axial_rel_err_mu = np.divide(
        np.abs(mu_matches["analytic_axial"] - mu_matches["fem_axial"]),
        np.abs(mu_matches["fem_axial"]),
        out=np.full_like(mu_matches["fem_axial"], np.nan),
        where=np.abs(mu_matches["fem_axial"]) > 0,
    )
    axial_rel_err_radius = np.divide(
        np.abs(axial_radius_sweep["analytic_axial"] - axial_radius_sweep["fem_axial"]),
        np.abs(axial_radius_sweep["fem_axial"]),
        out=np.full_like(axial_radius_sweep["fem_axial"], np.nan),
        where=np.abs(axial_radius_sweep["fem_axial"]) > 0,
    )

    mu_axial_ranges = constant_id_ranges(mu_values, mu_matches["fem_axial_mode_ids"])
    radius_axial_ranges = constant_id_ranges(axial_radius_sweep["radius_mm"], axial_radius_sweep["fem_mode_ids"])

    print(f"FEM source for beta=0 low-mode sweep: {fem_source}")
    print(f"mu grid: {mu_values[0]:.6f} .. {mu_values[-1]:.6f} ({len(mu_values)} points)")
    print(f"analysis mu base step: {ANALYSIS_MU_STEP:.4f}")
    print("local mu refinement windows: none")
    print(
        "type-aware matching: bending rows use the first six modes classified as bending; "
        "axial rows use the first mode classified as axial."
    )
    print(
        "analytic axial branch: seeded once from the base mu=0 determinant root at the first FEM-identified axial slot; "
        "for beta=0 and fixed L_total this branch is carried as the decoupled axial branch across mu and r."
    )
    print(f"saved plot: {BENDING_PLOT_PATH}")
    print(f"saved plot: {AXIAL_MU_PLOT_PATH}")
    print(f"saved plot: {AXIAL_RADIUS_PLOT_PATH}")
    print(f"saved mu=0 table: {MU0_TABLE_PATH}")
    print(f"saved mode-type table: {MODE_TYPES_PATH}")
    print(f"saved type-aware mu table: {TYPE_MATCH_MU_PATH}")
    print(f"saved type-aware radius table: {TYPE_MATCH_RADIUS_PATH}")
    print(f"saved axial radius sweep table: {AXIAL_RADIUS_SWEEP_PATH}")
    print(f"saved low-radius summary tables: {RADIUS_LOW_PATH}, {RADIUS_AXIAL_PATH}")
    print(
        f"existing FEM CSV vs type-aware FEM bending check: max_abs_diff={np.nanmax(csv_vs_typed):.6e} Hz, "
        f"max_rel_diff={np.nanmax(np.divide(csv_vs_typed, np.abs(fem_low_from_csv), out=np.full_like(csv_vs_typed, np.nan), where=np.abs(fem_low_from_csv) > 0)):.6e}"
    )
    print(
        f"bending mu-sweep error stats: max_abs_err={np.nanmax(bending_abs_err):.6f} Hz, "
        f"max_rel_err={np.nanmax(bending_rel_err):.6e}, mean_rel_err={np.nanmean(bending_rel_err):.6e}"
    )
    print(
        f"axial mu-sweep error stats: max_rel_err={np.nanmax(axial_rel_err_mu):.6e}, "
        f"mean_rel_err={np.nanmean(axial_rel_err_mu):.6e}"
    )
    print(
        f"axial radius-sweep error stats: max_rel_err={np.nanmax(axial_rel_err_radius):.6e}, "
        f"mean_rel_err={np.nanmean(axial_rel_err_radius):.6e}"
    )
    print("first axial FEM mode index over mu:", mu_axial_ranges)
    print("first axial FEM mode index over radius (mm):", radius_axial_ranges)
    print_mu0_table(mu0_rows)
    print_radius_summary(radius_low_matches["low_rows"], radius_low_matches["axial_rows"])


if __name__ == "__main__":
    main()
