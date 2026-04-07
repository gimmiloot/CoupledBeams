from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from formulas import BeamParams, det_clamped_coupled, frequency_scale, lambdas_to_frequencies
    from solvers import (
        find_close_pair_candidates,
        find_first_n_roots as _shared_find_first_n_roots,
        find_roots_scan_bisect as _shared_find_roots_scan_bisect,
        golden_minimize,
        refine_tracked_pair,
        root_by_min_abs_det,
        track_branches,
        tracked_lambdas_vs_mu,
    )
else:
    from .formulas import BeamParams, det_clamped_coupled, frequency_scale, lambdas_to_frequencies
    from .solvers import (
        find_close_pair_candidates,
        find_first_n_roots as _shared_find_first_n_roots,
        find_roots_scan_bisect as _shared_find_roots_scan_bisect,
        golden_minimize,
        refine_tracked_pair,
        root_by_min_abs_det,
        track_branches,
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


def integer_ratio_marks(mu_max: float, n_max: int = 200):
    mus = []
    labels = []
    for n in range(1, n_max + 1):
        mu = (n - 1) / (n + 1)
        if mu <= mu_max + 1e-12:
            mus.append(mu)
            labels.append(f"1:{n}")
        else:
            break
    return np.array(mus), labels


def build_mu_sweep_data(
    params: BeamParams,
    beta_deg: float,
    mu_min: float,
    mu_max: float,
    mu_step: float,
    n_modes: int,
    Lmax0: float,
    scan_step: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mu_grid = np.arange(mu_min, mu_max + 1e-12, mu_step)
    lambdas_tr = tracked_lambdas_vs_mu(
        params=params,
        beta_deg=beta_deg,
        mu_values=mu_grid,
        n_modes=n_modes,
        Lmin=0.2,
        Lmax0=Lmax0,
        scan_step=scan_step,
        grow_factor=1.35,
        max_tries=8,
        tracking_method="auto",
    )
    f_tr = lambdas_to_frequencies(lambdas_tr, params)
    return mu_grid, lambdas_tr, f_tr


def plot_tracked_frequencies_vs_mu(
    mu_grid: np.ndarray,
    f_tr: np.ndarray,
    beta_deg: float,
    show_ratio_lines: bool,
    label_ratio_upto: int,
) -> None:
    mu_marks, labels = integer_ratio_marks(mu_grid[-1])

    plt.figure(figsize=(11, 6.5))
    for i in range(f_tr.shape[0]):
        plt.plot(mu_grid, f_tr[i], label=f"branch {i + 1}")

    if show_ratio_lines:
        for mu_m in mu_marks:
            plt.axvline(mu_m, linestyle=":", linewidth=1.0, alpha=0.45)

        for k, mu_m in enumerate(mu_marks[:label_ratio_upto]):
            plt.text(mu_m, plt.ylim()[0], labels[k], rotation=90, va="bottom", ha="right", fontsize=8)
        if len(mu_marks) > label_ratio_upto:
            plt.text(mu_marks[-1], plt.ylim()[0], labels[-1], rotation=90, va="bottom", ha="right", fontsize=8)

    plt.xlabel("mu  (L1=1-mu, L2=1+mu, L1+L2=2)")
    plt.ylabel("Frequency (Hz)")
    plt.title(
        f"First {f_tr.shape[0]} frequencies vs mu at fixed beta={beta_deg} deg (clamped ends)\n"
        "Branch colors are tracked (stay the same even if branches change order)"
    )
    plt.grid(True, alpha=0.3)
    plt.legend(ncols=2, fontsize=9)
    plt.tight_layout()
    plt.show()


def plot_refined_pair_result(result: dict[str, object], beta_deg: float) -> None:
    i = int(result["i"])
    j = int(result["j"])
    mu_zoom = result["mu_zoom"]
    fi_ref = result["fi_ref"]
    fj_ref = result["fj_ref"]
    mu_min_gap = float(result["mu_min_gap"])

    plt.figure(figsize=(10, 5.2))
    plt.plot(mu_zoom, fi_ref, label=f"refined branch {i + 1}")
    plt.plot(mu_zoom, fj_ref, label=f"refined branch {j + 1}")
    plt.axvline(mu_min_gap, linestyle="--", alpha=0.6)
    plt.xlabel("mu")
    plt.ylabel("Frequency (Hz)")
    plt.title(f"Refined pair ({i + 1},{j + 1}) around mu* (beta={beta_deg} deg)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_gap_result(result: dict[str, object]) -> None:
    i = int(result["i"])
    j = int(result["j"])
    mu_zoom = result["mu_zoom"]
    gap_ref = result["gap_ref"]
    gap_tol = float(result["gap_tol"])
    mu_min_gap = float(result["mu_min_gap"])

    plt.figure(figsize=(10, 4.8))
    plt.plot(mu_zoom, gap_ref, label="gap = |fj-fi|")
    plt.axhline(gap_tol, linestyle="--", alpha=0.7, label="gap_tol")
    plt.axvline(mu_min_gap, linestyle="--", alpha=0.6)
    plt.xlabel("mu")
    plt.ylabel("Gap (Hz)")
    plt.title(f"Gap diagnostic for pair ({i + 1},{j + 1})")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_det_sanity_result(result: dict[str, object]) -> None:
    mu_zoom = result["mu_zoom"]
    detabs_i = result["detabs_i"]
    detabs_j = result["detabs_j"]
    i = int(result["i"])
    j = int(result["j"])

    plt.figure(figsize=(10, 4.4))
    plt.semilogy(mu_zoom, np.maximum(detabs_i, 1e-300), label=f"|det| at refined branch {i + 1}")
    plt.semilogy(mu_zoom, np.maximum(detabs_j, 1e-300), label=f"|det| at refined branch {j + 1}")
    plt.xlabel("mu")
    plt.ylabel("|det| (log scale)")
    plt.title("Sanity: how close we are to det=0 in refinement")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


def _coarse_gap_summary(
    params: BeamParams,
    mu_grid: np.ndarray,
    lambdas_tr: np.ndarray,
    i: int,
    j: int,
) -> tuple[int, float, float]:
    fi = lambdas_to_frequencies(lambdas_tr[i], params)
    fj = lambdas_to_frequencies(lambdas_tr[j], params)
    gap = np.abs(fj - fi)
    k0 = int(np.nanargmin(gap))
    return k0, float(mu_grid[k0]), float(gap[k0])


def diagnose_crossing_vs_avoided(
    params: BeamParams,
    beta_deg: float,
    mu_grid: np.ndarray,
    lambdas_tr: np.ndarray,
    branch2_idx: int = 1,
    branch3_idx: int = 2,
    mu_half_width: float = 0.04,
    mu_step: float = 2e-4,
    lambda_window: float = 0.35,
    tol_lambda: float = 1e-10,
):
    k0, mu_star, gap_coarse = _coarse_gap_summary(params, mu_grid, lambdas_tr, branch2_idx, branch3_idx)
    print(f"[diagnose] coarse mu* ~ {mu_star:.6f}, min gap ~ {gap_coarse:.6g} Hz")

    L2_prev = float(lambdas_tr[branch2_idx, k0])
    L3_prev = float(lambdas_tr[branch3_idx, k0])
    if not (np.isfinite(L2_prev) and np.isfinite(L3_prev) and L2_prev < L3_prev):
        raise RuntimeError("Bad initialization for refinement (NaN or wrong order).")

    result = refine_tracked_pair(
        params=params,
        beta_deg=beta_deg,
        mu_grid=mu_grid,
        lambdas_tr=lambdas_tr,
        i=branch2_idx,
        j=branch3_idx,
        mu_half_width=mu_half_width,
        mu_step=mu_step,
        lambda_window=lambda_window,
        tol_lambda=tol_lambda,
    )
    if result["verdict"] == "FAILED":
        raise RuntimeError("Bad initialization for refinement (NaN or wrong order).")

    print(f"[diagnose] refined mu_min_gap = {result['mu_min_gap']:.8f}")
    print(f"[diagnose] refined min gap    = {result['gap_min']:.10g} Hz")
    print(f"[diagnose] numerical gap_tol  = {result['gap_tol']:.10g} Hz")
    print("=>", result["verbose_verdict"])

    plot_refined_pair_result(result, beta_deg)
    plot_gap_result(result)
    plot_det_sanity_result(result)
    return result["verbose_verdict"]


def diagnose_pair_refined(
    params: BeamParams,
    beta_deg: float,
    mu_grid: np.ndarray,
    lambdas_tr: np.ndarray,
    i: int,
    j: int,
    mu_half_width: float,
    mu_step: float,
    lambda_window: float,
    tol_lambda: float,
    make_plots: bool,
):
    result = refine_tracked_pair(
        params=params,
        beta_deg=beta_deg,
        mu_grid=mu_grid,
        lambdas_tr=lambdas_tr,
        i=i,
        j=j,
        mu_half_width=mu_half_width,
        mu_step=mu_step,
        lambda_window=lambda_window,
        tol_lambda=tol_lambda,
    )

    if result["verdict"] == "FAILED":
        return "FAILED", {"mu_min_gap": np.nan, "gap_min": np.nan, "gap_tol": np.nan}

    if make_plots:
        plot_refined_pair_result(result, beta_deg)
        plot_gap_result(result)

    info = {
        "mu_min_gap": float(result["mu_min_gap"]),
        "gap_min": float(result["gap_min"]),
        "gap_tol": float(result["gap_tol"]),
    }
    return result["verdict"], info


def auto_crossing_diagnostics_close_pairs(
    params: BeamParams,
    beta_deg: float,
    mu_grid: np.ndarray,
    lambdas_tr: np.ndarray,
    gap_threshold_hz: float = 5.0,
    gap_rel_threshold: float | None = None,
    max_pairs: int = 6,
    mu_half_width: float = 0.04,
    mu_step: float = 2e-4,
    lambda_window: float = 0.35,
    tol_lambda: float = 1e-10,
    make_plots: bool = True,
):
    freqs_tr = lambdas_to_frequencies(lambdas_tr, params)
    candidates = find_close_pair_candidates(
        freqs_tr=freqs_tr,
        mu_grid=mu_grid,
        gap_threshold_hz=gap_threshold_hz,
        gap_rel_threshold=gap_rel_threshold,
        max_pairs=max_pairs,
    )

    if len(candidates) == 0:
        print("[auto-diagnose] No close pairs found under the thresholds.")
        return []

    print("\n=== AUTO DIAGNOSTICS (close pairs, coarse) ===")
    for c in candidates:
        print(f"pair ({c['i'] + 1},{c['j'] + 1}): mu*~{c['mu_star']:.6f}, coarse gap~{c['gap_coarse']:.6g} Hz")

    results = []
    for c in candidates:
        verdict, info = diagnose_pair_refined(
            params=params,
            beta_deg=beta_deg,
            mu_grid=mu_grid,
            lambdas_tr=lambdas_tr,
            i=c["i"],
            j=c["j"],
            mu_half_width=mu_half_width,
            mu_step=mu_step,
            lambda_window=lambda_window,
            tol_lambda=tol_lambda,
            make_plots=make_plots,
        )
        results.append({**c, **info, "verdict": verdict})

    print("\n=== AUTO DIAGNOSTICS (refined summary) ===")
    for r in results:
        print(
            f"pair ({r['i'] + 1},{r['j'] + 1}): "
            f"mu_min_gap={r['mu_min_gap']:.8f}, "
            f"min_gap={r['gap_min']:.6g} Hz, "
            f"gap_tol={r['gap_tol']:.3g} Hz -> {r['verdict']}"
        )

    return results


def plot_freqs_vs_mu(
    params: BeamParams,
    beta_deg: float = 5.0,
    mu_min: float = 0.0,
    mu_max: float = 0.9,
    mu_step: float = 0.01,
    n_modes: int = 6,
    Lmax0: float = 55.0,
    scan_step: float = 0.02,
    show_ratio_lines: bool = True,
    label_ratio_upto: int = 10,
    do_crossing_diagnose: bool = True,
):
    mu_grid, lambdas_tr, f_tr = build_mu_sweep_data(
        params=params,
        beta_deg=beta_deg,
        mu_min=mu_min,
        mu_max=mu_max,
        mu_step=mu_step,
        n_modes=n_modes,
        Lmax0=Lmax0,
        scan_step=scan_step,
    )

    plot_tracked_frequencies_vs_mu(
        mu_grid=mu_grid,
        f_tr=f_tr,
        beta_deg=beta_deg,
        show_ratio_lines=show_ratio_lines,
        label_ratio_upto=label_ratio_upto,
    )

    if do_crossing_diagnose:
        auto_crossing_diagnostics_close_pairs(
            params=params,
            beta_deg=beta_deg,
            mu_grid=mu_grid,
            lambdas_tr=lambdas_tr,
            gap_threshold_hz=5.0,
            gap_rel_threshold=None,
            max_pairs=6,
            make_plots=True,
        )


def main() -> None:
    p = BeamParams(
        E=2.1e11,
        rho=7800.0,
        r=0.005,
        L_total=2.0,
    )
    print("eps =", p.eps)

    plot_freqs_vs_mu(
        params=p,
        beta_deg=20.0,
        mu_min=0.0,
        mu_max=0.9,
        mu_step=0.01,
        n_modes=6,
        Lmax0=55.0,
        scan_step=0.02,
        show_ratio_lines=True,
        label_ratio_upto=10,
        do_crossing_diagnose=True,
    )


if __name__ == "__main__":
    main()
