from pathlib import Path
import sys

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from formulas import BeamParams, det_clamped_coupled
    from solvers import (
        find_first_n_roots as _shared_find_first_n_roots,
        find_roots_scan_bisect as _shared_find_roots_scan_bisect,
        track_branches as _shared_track_branches,
        tracked_lambdas_vs_mu,
    )
else:
    from .formulas import BeamParams, det_clamped_coupled
    from .solvers import (
        find_first_n_roots as _shared_find_first_n_roots,
        find_roots_scan_bisect as _shared_find_roots_scan_bisect,
        track_branches as _shared_track_branches,
        tracked_lambdas_vs_mu,
    )


def _find_roots_scan_bisect(
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
    Lmax0: float = 55.0,
    scan_step: float = 0.02,
    grow_factor: float = 1.35,
    max_tries: int = 8,
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


def track_branches(values_sorted: np.ndarray) -> np.ndarray:
    return _shared_track_branches(values_sorted, method="greedy")


def coupled_lambda_vs_mu(params: BeamParams, beta_deg: float, mu: np.ndarray, n_modes: int = 6) -> np.ndarray:
    return tracked_lambdas_vs_mu(
        params=params,
        beta_deg=beta_deg,
        mu_values=mu,
        n_modes=n_modes,
        Lmin=0.2,
        Lmax0=55.0,
        scan_step=0.02,
        grow_factor=1.35,
        max_tries=8,
        tracking_method="greedy",
    )


# ============================================================
# Single beam roots (Euler–Bernoulli)
# ============================================================
def roots_clamped_supported(n: int) -> np.ndarray:
    """Roots for clamped-supported (CS): tan(alpha) = tanh(alpha)."""
    mp.mp.dps = 60

    def f(a):
        return mp.tan(a) - mp.tanh(a)

    roots = []
    for k in range(1, n + 1):
        a0 = k * mp.pi - mp.mpf("0.2")
        roots.append(float(mp.findroot(f, a0)))
    return np.array(roots, dtype=float)


def single_lambda(alpha: np.ndarray, l: float, L: np.ndarray) -> np.ndarray:
    """Since \bar f = alpha^2 (l/L)^2 and \bar f = Lambda^2, we plot Lambda = alpha (l/L)."""
    return alpha[:, None] * l / L[None, :]


# ============================================================
# Combined plot
# ============================================================
def plot_combined(
    beta_deg_coupled: float = 15.0,
    mu_min: float = 0.0,
    mu_max: float = 0.9,
    mu_step: float = 0.01,
    n_coupled: int = 6,
    n_Lplus_CS: int = 6,
    n_Lminus_CS: int = 3,
    y_max: int = 13,
    save_path: str | None = None,
):
    """
    One combined plot:
      - first 6 coupled-rods frequencies, bright solid lines;
      - first 6 single-beam CS frequencies for L+ = l(1+mu), muted dashed lines;
      - first 3 single-beam CS frequencies for L- = l(1-mu), muted dashed lines.

    The vertical axis is Lambda (not Lambda^2).
    The angle beta_deg_coupled remains adjustable in the function call.
    """
    params = BeamParams(E=2.1e11, rho=7800.0, r=0.04, L_total=2.0)
    l = params.L_base
    mu = np.arange(mu_min, mu_max + 1e-12, mu_step)

    Lambda_c = coupled_lambda_vs_mu(params, beta_deg_coupled, mu, n_modes=n_coupled)

    a_cs = roots_clamped_supported(max(n_Lplus_CS, n_Lminus_CS))
    Lp = l * (1.0 + mu)
    Lm = l * (1.0 - mu)
    Lambda_Lp_CS = single_lambda(a_cs[:n_Lplus_CS], l, Lp)
    Lambda_Lm_CS = single_lambda(a_cs[:n_Lminus_CS], l, Lm)

    fig, ax = plt.subplots(figsize=(11.5, 6.6))

    for i in range(Lambda_c.shape[0]):
        ax.plot(mu, Lambda_c[i], lw=2.4, ls="-")

    bright_colors = [ax.lines[i].get_color() for i in range(n_coupled)]

    for i in range(n_Lplus_CS):
        c = to_rgba(bright_colors[i], alpha=0.48)
        ax.plot(mu, Lambda_Lp_CS[i], lw=1.8, ls="--", color=c)

    for i in range(n_Lminus_CS):
        c = to_rgba(bright_colors[i], alpha=0.48)
        ax.plot(mu, Lambda_Lm_CS[i], lw=1.8, ls=(0, (2.0, 2.4)), color=c)

    ax.set_xlabel("μ")
    ax.set_ylabel(r"Безразмерная частота $\,\Lambda $")
    if y_max is not None:
        ax.set_ylim(0, y_max)
    ax.grid(True, alpha=0.3)
    ax.set_title(
        f"Сопряжённые стержни (β={beta_deg_coupled}°) и одиночные стержни CS, "
        fr"график для $\,\Lambda$, r={params.r}"
    )

    handles = [
        Line2D([0], [0], color="black", lw=2.4, ls="-", label="сплошная — сопряжённые стержни, первые 6 частот"),
        Line2D(
            [0],
            [0],
            color=to_rgba("black", alpha=0.38),
            lw=1.8,
            ls="--",
            label=r"пунктир — одиночный стержень, CS (заделка–шарнир), $L^+=\ell(1+\mu)$, первые 6",
        ),
        Line2D(
            [0],
            [0],
            color=to_rgba("black", alpha=0.28),
            lw=1.8,
            ls=(0, (2.0, 2.4)),
            label=r"пунктир — одиночный стержень, CS (заделка–шарнир), $L^-=\ell(1-\mu)$, первые 3",
        ),
    ]
    ax.legend(handles=handles, fontsize=9, loc="upper left")

    fig.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=200, bbox_inches="tight")
    plt.show()


def main() -> None:
    plot_combined(beta_deg_coupled=15.0)


if __name__ == "__main__":
    main()