import argparse
import inspect
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from matplotlib.colors import to_rgba
from matplotlib.lines import Line2D


REPO_ROOT = Path(__file__).resolve().parents[3]
SRC_ROOT = REPO_ROOT / "src"
DEFAULT_OUTPUT_PATH = REPO_ROOT / "results" / "lambda_mu_fixed_beta.png"
DEFAULT_EPSILON = 0.0025
DEFAULT_L_TOTAL = 2.0
DEFAULT_E = 2.1e11
DEFAULT_RHO = 7800.0

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from formulas import BeamParams, det_clamped_coupled
    from solvers import (
        find_first_n_roots as _shared_find_first_n_roots,
        find_roots_scan_bisect as _shared_find_roots_scan_bisect,
        track_branches as _shared_track_branches,
    )
else:
    from .formulas import BeamParams, det_clamped_coupled
    from .solvers import (
        find_first_n_roots as _shared_find_first_n_roots,
        find_roots_scan_bisect as _shared_find_roots_scan_bisect,
        track_branches as _shared_track_branches,
    )

from scripts.lib.analytic_branch_tracking import (  # noqa: E402
    DEFAULT_MU_STEPS,
    DEFAULT_N_SOLVE,
    DEFAULT_N_TRACK,
    branch_id_from_base_sorted_index,
    dense_mu_values_for_targets,
    track_mu_sweep,
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


def coupled_lambda_vs_mu(
    params: BeamParams,
    beta_deg: float,
    mu: np.ndarray,
    n_modes: int = 6,
    *,
    allow_low_mac: bool = False,
) -> np.ndarray:
    n_track = max(DEFAULT_N_TRACK, int(n_modes))
    requested_mu = np.asarray(mu, dtype=float)
    tracking_mu = dense_mu_values_for_targets(
        requested_mu,
        mu_steps=max(DEFAULT_MU_STEPS, int(len(requested_mu))),
    )
    branch_ids = [branch_id_from_base_sorted_index(index) for index in range(1, int(n_modes) + 1)]
    result = track_mu_sweep(
        epsilon=params.eps,
        beta=beta_deg,
        mu_values=tracking_mu,
        n_track=n_track,
        n_solve=max(DEFAULT_N_SOLVE, n_track),
        shape_metric="full",
        allow_low_mac=allow_low_mac,
        required_branch_ids=branch_ids,
    )
    return result.lambda_grid(branch_ids, requested_mu, beta=float(beta_deg))


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


def beam_params_from_epsilon(
    epsilon: float,
    L_total: float = DEFAULT_L_TOTAL,
    E: float = DEFAULT_E,
    rho: float = DEFAULT_RHO,
) -> BeamParams:
    """Build project beam parameters from epsilon = sqrt(I/S) / l for a circular section."""
    l_base = 0.5 * float(L_total)
    radius = 2.0 * l_base * float(epsilon)
    return BeamParams(E=float(E), rho=float(rho), r=radius, L_total=float(L_total))


# ============================================================
# Combined plot
# ============================================================
def plot_combined(
    beta_deg_coupled: float = 45.0,
    epsilon: float = 0.01,
    L_total: float = DEFAULT_L_TOTAL,
    mu_min: float = 0.0,
    mu_max: float = 0.9,
    mu_step: float = 0.01,
    n_coupled: int = 10,
    n_Lplus_CS: int = 9,
    n_Lminus_CS: int = 5,
    n_dashed_lines: int | None = None,
    y_max: int = 20,
    save_path: str | Path | None = None,
    show: bool = True,
    allow_low_mac: bool = False,
) -> Path | None:
    """
    One combined plot:
      - first n_coupled coupled-rods frequencies, bright solid lines;
      - single-beam CS frequencies for L+ = l(1+mu), muted dashed lines;
      - single-beam CS frequencies for L- = l(1-mu), muted dashed lines.

    The vertical axis is Lambda (not Lambda^2).
    The angle, epsilon, number of coupled branches, and number of dashed
    reference curves remain adjustable in the function call or CLI.
    """
    if n_dashed_lines is not None:
        n_Lplus_CS = int(n_dashed_lines)
        n_Lminus_CS = int(n_dashed_lines)
    if n_coupled < 1:
        raise ValueError("n_coupled must be positive.")
    if n_Lplus_CS < 0 or n_Lminus_CS < 0:
        raise ValueError("Dashed reference-line counts must be non-negative.")

    params = beam_params_from_epsilon(epsilon=epsilon, L_total=L_total)
    l = params.L_base
    mu = np.arange(mu_min, mu_max + 1e-12, mu_step)

    Lambda_c = coupled_lambda_vs_mu(
        params,
        beta_deg_coupled,
        mu,
        n_modes=n_coupled,
        allow_low_mac=allow_low_mac,
    )

    n_reference = max(n_Lplus_CS, n_Lminus_CS)
    a_cs = roots_clamped_supported(n_reference) if n_reference else np.array([], dtype=float)
    Lp = l * (1.0 + mu)
    Lm = l * (1.0 - mu)
    Lambda_Lp_CS = single_lambda(a_cs[:n_Lplus_CS], l, Lp)
    Lambda_Lm_CS = single_lambda(a_cs[:n_Lminus_CS], l, Lm)

    fig, ax = plt.subplots(figsize=(11.5, 6.6))

    for i in range(Lambda_c.shape[0]):
        ax.plot(mu, Lambda_c[i], lw=2.4, ls="-")

    bright_colors = [ax.lines[i].get_color() for i in range(n_coupled)]

    for i in range(n_Lplus_CS):
        c = to_rgba(bright_colors[i % len(bright_colors)], alpha=0.48)
        ax.plot(mu, Lambda_Lp_CS[i], lw=1.8, ls="--", color=c)

    for i in range(n_Lminus_CS):
        c = to_rgba(bright_colors[i % len(bright_colors)], alpha=0.48)
        ax.plot(mu, Lambda_Lm_CS[i], lw=1.8, ls=(0, (2.0, 2.4)), color=c)

    ax.set_xlabel("μ")
    ax.set_ylabel(r"Безразмерная частота $\,\Lambda $")
    if y_max is not None:
        ax.set_ylim(0, y_max)
    ax.grid(True, alpha=0.3)
    ax.set_title(
        f"Сопряжённые стержни (β={beta_deg_coupled}°) и одиночные стержни CS, "
        fr"график для $\,\Lambda$, $\varepsilon={params.eps:.5g}$"
    )

    handles = [
        Line2D(
            [0],
            [0],
            color="black",
            lw=2.4,
            ls="-",
            label=f"сплошная — сопряжённые стержни, первые {n_coupled} частот",
        ),
        Line2D(
            [0],
            [0],
            color=to_rgba("black", alpha=0.38),
            lw=1.8,
            ls="--",
            label=rf"пунктир — одиночный стержень, CS (заделка–шарнир), $L^+=\ell(1+\mu)$, первые {n_Lplus_CS}",
        ),
        Line2D(
            [0],
            [0],
            color=to_rgba("black", alpha=0.28),
            lw=1.8,
            ls=(0, (2.0, 2.4)),
            label=rf"пунктир — одиночный стержень, CS (заделка–шарнир), $L^-=\ell(1-\mu)$, первые {n_Lminus_CS}",
        ),
    ]
    ax.legend(handles=handles, fontsize=9, loc="upper left")

    fig.tight_layout()

    saved_path: Path | None = None
    if save_path is not None:
        saved_path = Path(save_path)
        saved_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(saved_path, dpi=200, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return saved_path


def plot_combined_default(parameter_name: str) -> object:
    """Return the current default from plot_combined by parameter name."""
    return inspect.signature(plot_combined).parameters[parameter_name].default


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot Lambda(mu) at fixed beta with single-rod CS dashed reference families.",
    )
    parser.add_argument(
        "--beta",
        type=float,
        default=None,
        help=f"fixed coupling angle in degrees; omitted value uses plot_combined default {plot_combined_default('beta_deg_coupled')}",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=None,
        help=f"thickness parameter; omitted value uses plot_combined default {plot_combined_default('epsilon')}",
    )
    parser.add_argument(
        "--num-modes",
        type=int,
        default=None,
        help=f"number of coupled analytic branches to plot; omitted value uses plot_combined default {plot_combined_default('n_coupled')}",
    )
    parser.add_argument(
        "--num-dashed-lines",
        type=int,
        default=None,
        help="number of dashed CS roots per family; overrides the two family-specific counts",
    )
    parser.add_argument(
        "--num-lplus-dashed-lines",
        type=int,
        default=None,
        help=f"number of dashed CS roots for L+ = l(1+mu); omitted value uses plot_combined default {plot_combined_default('n_Lplus_CS')}",
    )
    parser.add_argument(
        "--num-lminus-dashed-lines",
        type=int,
        default=None,
        help=f"number of dashed CS roots for L- = l(1-mu); omitted value uses plot_combined default {plot_combined_default('n_Lminus_CS')}",
    )
    parser.add_argument(
        "--mu-min",
        type=float,
        default=None,
        help=f"minimum mu value; omitted value uses plot_combined default {plot_combined_default('mu_min')}",
    )
    parser.add_argument(
        "--mu-max",
        type=float,
        default=None,
        help=f"maximum mu value; omitted value uses plot_combined default {plot_combined_default('mu_max')}",
    )
    parser.add_argument(
        "--mu-step",
        type=float,
        default=None,
        help=f"mu grid step; omitted value uses plot_combined default {plot_combined_default('mu_step')}",
    )
    parser.add_argument(
        "--y-max",
        type=float,
        default=None,
        help=f"upper y-axis limit; use <= 0 for automatic; omitted value uses plot_combined default {plot_combined_default('y_max')}",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT_PATH, help="PNG output path")
    parser.add_argument("--show", action="store_true", help="also display the figure window after saving")
    parser.add_argument(
        "--allow-low-mac",
        action="store_true",
        help="Allow exploratory plots even if analytic branch assignment falls below the MAC warning threshold.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    plot_kwargs: dict[str, object] = {
        "save_path": args.output,
        "show": args.show,
        "allow_low_mac": args.allow_low_mac,
    }
    if args.beta is not None:
        plot_kwargs["beta_deg_coupled"] = args.beta
    if args.epsilon is not None:
        plot_kwargs["epsilon"] = args.epsilon
    if args.mu_min is not None:
        plot_kwargs["mu_min"] = args.mu_min
    if args.mu_max is not None:
        plot_kwargs["mu_max"] = args.mu_max
    if args.mu_step is not None:
        plot_kwargs["mu_step"] = args.mu_step
    if args.num_modes is not None:
        plot_kwargs["n_coupled"] = args.num_modes
    if args.num_lplus_dashed_lines is not None:
        plot_kwargs["n_Lplus_CS"] = args.num_lplus_dashed_lines
    if args.num_lminus_dashed_lines is not None:
        plot_kwargs["n_Lminus_CS"] = args.num_lminus_dashed_lines
    if args.num_dashed_lines is not None:
        plot_kwargs["n_dashed_lines"] = args.num_dashed_lines
    if args.y_max is not None:
        plot_kwargs["y_max"] = None if args.y_max <= 0 else args.y_max

    saved_path = plot_combined(**plot_kwargs)
    if saved_path is not None:
        print(f"saved figure: {saved_path}")
    print(f"beta={plot_kwargs.get('beta_deg_coupled', plot_combined_default('beta_deg_coupled'))} deg")
    print(f"epsilon={plot_kwargs.get('epsilon', plot_combined_default('epsilon'))}")
    print(f"num_modes={plot_kwargs.get('n_coupled', plot_combined_default('n_coupled'))}")
    if args.num_dashed_lines is not None:
        print(f"num_dashed_lines_per_family={args.num_dashed_lines}")
    else:
        print(f"num_lplus_dashed_lines={plot_kwargs.get('n_Lplus_CS', plot_combined_default('n_Lplus_CS'))}")
        print(f"num_lminus_dashed_lines={plot_kwargs.get('n_Lminus_CS', plot_combined_default('n_Lminus_CS'))}")


if __name__ == "__main__":
    main()
