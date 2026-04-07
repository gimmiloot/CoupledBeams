from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from formulas import BeamParams, det_clamped_coupled, lambdas_to_frequencies, segment_lengths
    from solvers import (
        find_first_n_roots as _shared_find_first_n_roots,
        find_roots_scan_bisect as _shared_find_roots_scan_bisect,
        fixed_fixed_lambdas,
        track_branches,
    )
else:
    from .formulas import BeamParams, det_clamped_coupled, lambdas_to_frequencies, segment_lengths
    from .solvers import (
        find_first_n_roots as _shared_find_first_n_roots,
        find_roots_scan_bisect as _shared_find_roots_scan_bisect,
        fixed_fixed_lambdas,
        track_branches,
    )


def find_roots_scan_bisect(
    beta: float,
    mu: float,
    eps: float,
    n_roots: int,
    Lmin: float,
    Lmax: float,
    scan_step: float,
    bisect_iters: int = 70,
) -> list[float]:
    return _shared_find_roots_scan_bisect(beta, mu, eps, n_roots, Lmin, Lmax, scan_step, bisect_iters)


def find_first_n_roots(
    beta: float,
    mu: float,
    eps: float,
    n_roots: int,
    Lmin: float = 0.2,
    Lmax0: float = 35.0,
    scan_step: float = 0.03,
    grow_factor: float = 1.35,
    max_tries: int = 7,
) -> np.ndarray:
    return _shared_find_first_n_roots(
        beta,
        mu,
        eps,
        n_roots,
        Lmin=Lmin,
        Lmax0=Lmax0,
        scan_step=scan_step,
        grow_factor=grow_factor,
        max_tries=max_tries,
    )


def build_beta_sweep_data(
    params: BeamParams,
    mu: float,
    beta_deg_grid: np.ndarray | None,
    n_modes_main: int,
    Lmax0: float,
    scan_step: float,
) -> tuple[np.ndarray, np.ndarray]:
    if beta_deg_grid is None:
        beta_deg_grid = np.arange(0, 91, 1)
    beta_deg_grid = np.asarray(beta_deg_grid, dtype=float)
    beta_rad_grid = np.deg2rad(beta_deg_grid)

    freqs_raw = np.full((n_modes_main, len(beta_rad_grid)), np.nan, dtype=float)
    for j, beta in enumerate(beta_rad_grid):
        roots = find_first_n_roots(
            beta,
            mu,
            params.eps,
            n_roots=n_modes_main,
            Lmin=0.2,
            Lmax0=Lmax0,
            scan_step=scan_step,
        )
        freqs_raw[:, j] = lambdas_to_frequencies(roots, params)

    freqs_sorted = np.sort(freqs_raw, axis=0)
    return beta_deg_grid, track_branches(freqs_sorted)


def build_long_rod_reference(params: BeamParams, mu: float, n_dashed_main: int) -> np.ndarray:
    L1, L2 = segment_lengths(params, mu)
    L_long = max(L1, L2)
    lam_d = fixed_fixed_lambdas(n_dashed_main)
    return (lam_d**2 / (L_long**2)) * params.f_scale


def plot_branches_vs_beta(
    beta_deg_grid: np.ndarray,
    freqs_tr: np.ndarray,
    mu: float,
    n_modes_to_plot: int,
    dashed_freqs: np.ndarray | None = None,
    dashed_count: int = 0,
) -> None:
    plt.figure(figsize=(10.5, 6.2))
    for i in range(n_modes_to_plot):
        plt.plot(beta_deg_grid, freqs_tr[i], label=f"branch {i + 1}")

    if dashed_freqs is not None:
        for i in range(min(dashed_count, len(dashed_freqs))):
            plt.hlines(
                dashed_freqs[i],
                beta_deg_grid[0],
                beta_deg_grid[-1],
                linestyles="--",
                linewidth=1.2 if n_modes_to_plot > 6 else 1.3,
                label="long rod (fixed-fixed)" if i == 0 else None,
            )

    plt.xlabel("beta (degrees)")
    plt.ylabel("Frequency (Hz)")
    title = f"Coupled (clamped): first {n_modes_to_plot} branches vs beta (0..90 deg), mu={mu:+.3f}"
    if dashed_freqs is not None:
        title += f"\nDashed: long rod fixed-fixed (first {dashed_count} modes)"
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_vs_beta_for_mu(
    params: BeamParams,
    mu: float,
    beta_deg_grid=None,
    n_modes_main: int = 10,
    n_modes_small: int = 6,
    n_dashed_main: int = 5,
    n_dashed_small: int = 3,
    show_dashed: bool = True,
    Lmax0: float = 35.0,
    scan_step: float = 0.03,
):
    beta_deg_grid, freqs_tr = build_beta_sweep_data(
        params=params,
        mu=mu,
        beta_deg_grid=beta_deg_grid,
        n_modes_main=n_modes_main,
        Lmax0=Lmax0,
        scan_step=scan_step,
    )

    dashed_freqs = build_long_rod_reference(params, mu, n_dashed_main) if show_dashed else None
    plot_branches_vs_beta(
        beta_deg_grid=beta_deg_grid,
        freqs_tr=freqs_tr,
        mu=mu,
        n_modes_to_plot=n_modes_main,
        dashed_freqs=dashed_freqs,
        dashed_count=n_dashed_main,
    )
    plot_branches_vs_beta(
        beta_deg_grid=beta_deg_grid,
        freqs_tr=freqs_tr,
        mu=mu,
        n_modes_to_plot=n_modes_small,
        dashed_freqs=dashed_freqs if show_dashed else None,
        dashed_count=n_dashed_small,
    )


def main() -> None:
    p = BeamParams(E=2.1e11, rho=7800.0, r=0.005, L_total=2.0)
    print("eps =", p.eps)

    for mu in [0]:
        plot_vs_beta_for_mu(
            params=p,
            mu=mu,
            beta_deg_grid=np.arange(0, 91, 1),
            show_dashed=True,
            n_dashed_main=5,
            n_dashed_small=3,
            scan_step=0.03,
            Lmax0=35.0,
        )


if __name__ == "__main__":
    main()
