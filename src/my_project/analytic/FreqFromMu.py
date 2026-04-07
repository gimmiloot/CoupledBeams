import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass

# ============================================================
# 1) Parameters (same as your "clamped / zadelka" setup)
# ============================================================
@dataclass
class BeamParams:
    E: float
    rho: float
    r: float
    L_total: float  # fixed: L1+L2

    @property
    def S(self): return np.pi * self.r**2

    @property
    def I(self): return np.pi * self.r**4 / 4.0

    @property
    def L_base(self):  # L=(L1+L2)/2
        return 0.5 * self.L_total

    @property
    def eps(self):     # for circle = r/L_total
        return np.sqrt(self.I / self.S) / self.L_base

    @property
    def f_scale(self):
        # f = (Lambda^2/L^2)*sqrt(EI/(rho S))/(2*pi)
        return np.sqrt(self.E * self.I / (self.rho * self.S)) / (2.0 * np.pi)


# ============================================================
# 2) Determinant det M (your clamped matrix "zadelka")
# ============================================================
def det_clamped_coupled(Lambda: float, beta: float, mu: float, eps: float) -> float:
    L1f, L2f = 1.0 - mu, 1.0 + mu
    x1, x2 = Lambda * L1f, Lambda * L2f

    Cd1 = np.cos(x1) - np.cosh(x1)
    Sd1 = np.sin(x1) - np.sinh(x1)
    Cs1 = np.cos(x1) + np.cosh(x1)
    Ss1 = np.sin(x1) + np.sinh(x1)

    Cd2 = np.cos(x2) - np.cosh(x2)
    Sd2 = np.sin(x2) - np.sinh(x2)
    Cs2 = np.cos(x2) + np.cosh(x2)
    Ss2 = np.sin(x2) + np.sinh(x2)

    th1 = eps * (Lambda**2) * L1f
    th2 = eps * (Lambda**2) * L2f

    cb, sb = np.cos(beta), np.sin(beta)

    M = np.array([
        [Cd1,               Sd1,               -Cd2*cb,             Sd2*cb,             0.0,            np.sin(th2)*sb],
        [0.0,               0.0,                Cd2*sb,            -Sd2*sb,        np.sin(th1),        np.sin(th2)*cb],
        [-Ss1,              Cd1,                -Ss2,               -Cd2,               0.0,            0.0],
        [-Cs1,             -Ss1,                 Cs2,               -Ss2,               0.0,            0.0],
        [-eps*Lambda*Sd1,    eps*Lambda*Cs1,    -eps*Lambda*Sd2*cb, -eps*Lambda*Cs2*cb, 0.0,           -np.cos(th2)*sb],
        [0.0,               0.0,                 eps*Lambda*Sd2*sb,  eps*Lambda*Cs2*sb, np.cos(th1),    -np.cos(th2)*cb]
    ], dtype=float)

    return float(np.linalg.det(M))


# ============================================================
# 3) Root finding: scan + bisection + auto-grow
# ============================================================
def _find_roots_scan_bisect(beta: float, mu: float, eps: float,
                            n_roots: int,
                            Lmin: float,
                            Lmax: float,
                            scan_step: float,
                            bisect_iters: int = 70) -> list[float]:
    grid = np.arange(Lmin, Lmax + scan_step, scan_step)
    vals = np.array([det_clamped_coupled(L, beta, mu, eps) for L in grid], dtype=float)

    roots = []
    def sgn(x): return 1 if x > 0 else (-1 if x < 0 else 0)

    for i in range(len(grid) - 1):
        a, b = grid[i], grid[i+1]
        fa, fb = vals[i], vals[i+1]
        sa, sb = sgn(fa), sgn(fb)

        if sa == 0:
            roots.append(a)
        elif sa * sb < 0:
            left, right = a, b
            fl = fa
            for _ in range(bisect_iters):
                mid = 0.5 * (left + right)
                fm = det_clamped_coupled(mid, beta, mu, eps)
                sm = sgn(fm)
                if sm == 0:
                    left = right = mid
                    break
                if sgn(fl) * sm < 0:
                    right = mid
                else:
                    left = mid
                    fl = fm
            roots.append(0.5 * (left + right))

        if len(roots) >= n_roots:
            break

    roots = sorted(set([float(r) for r in roots]))
    return roots[:n_roots]


def find_first_n_roots(beta: float, mu: float, eps: float,
                       n_roots: int,
                       Lmin: float = 0.2,
                       Lmax0: float = 55.0,
                       scan_step: float = 0.02,
                       grow_factor: float = 1.35,
                       max_tries: int = 8) -> np.ndarray:
    Lmax = Lmax0
    roots = []
    for _ in range(max_tries):
        roots = _find_roots_scan_bisect(beta, mu, eps, n_roots, Lmin, Lmax, scan_step)
        if len(roots) >= n_roots:
            break
        Lmax *= grow_factor

    out = np.full(n_roots, np.nan, dtype=float)
    out[:len(roots)] = roots
    return out


# ============================================================
# 4) Branch tracking along mu (keeps color even if order swaps)
# ============================================================
def track_branches(values_sorted: np.ndarray) -> np.ndarray:
    """
    values_sorted: (N_modes, N_mu), each column sorted by magnitude.
    Returns tracked branches (rows continuous in mu).
    """
    N, M = values_sorted.shape
    tracked = np.full_like(values_sorted, np.nan)
    tracked[:, 0] = values_sorted[:, 0]

    try:
        from scipy.optimize import linear_sum_assignment
        use_hungarian = True
    except Exception:
        use_hungarian = False

    for j in range(1, M):
        prev = tracked[:, j-1]
        cur  = values_sorted[:, j]

        prev_idx = np.where(np.isfinite(prev))[0]
        cur_idx  = np.where(np.isfinite(cur))[0]
        if len(prev_idx) == 0 or len(cur_idx) == 0:
            continue

        C = np.abs(prev[prev_idx][:, None] - cur[cur_idx][None, :])

        assign = {}
        if use_hungarian:
            r, c = linear_sum_assignment(C)
            for rr, cc in zip(r, c):
                assign[prev_idx[rr]] = cur_idx[cc]
        else:
            used_r, used_c = set(), set()
            while True:
                best = None; bestv = np.inf
                for rr in range(C.shape[0]):
                    if rr in used_r: continue
                    for cc in range(C.shape[1]):
                        if cc in used_c: continue
                        if C[rr, cc] < bestv:
                            bestv = C[rr, cc]; best = (rr, cc)
                if best is None:
                    break
                rr, cc = best
                used_r.add(rr); used_c.add(cc)
                assign[prev_idx[rr]] = cur_idx[cc]
                if len(used_r) == C.shape[0] or len(used_c) == C.shape[1]:
                    break

        for i in range(N):
            tracked[i, j] = cur[assign[i]] if i in assign else np.nan

    return tracked


# ============================================================
# 5) Integer ratio marks: L2/L1 = n => mu = (n-1)/(n+1)
# ============================================================
def integer_ratio_marks(mu_max: float, n_max: int = 200):
    mus = []
    labels = []
    for n in range(1, n_max+1):
        mu = (n - 1) / (n + 1)
        if mu <= mu_max + 1e-12:
            mus.append(mu)
            labels.append(f"1:{n}")
        else:
            break
    return np.array(mus), labels


# ============================================================
# 6) High-accuracy "crossing vs avoided crossing" refinement
#    (inserted into the same program)
# ============================================================
def golden_minimize(g, a, b, tol=1e-10, max_iter=200):
    phi = (1 + 5 ** 0.5) / 2
    invphi = 1 / phi
    c = b - (b - a) * invphi
    d = a + (b - a) * invphi
    gc = g(c); gd = g(d)
    for _ in range(max_iter):
        if abs(b - a) < tol:
            break
        if gc < gd:
            b, d, gd = d, c, gc
            c = b - (b - a) * invphi
            gc = g(c)
        else:
            a, c, gc = c, d, gd
            d = a + (b - a) * invphi
            gd = g(d)
    x = 0.5 * (a + b)
    return x, g(x), 0.5 * (b - a)

def root_by_min_abs_det(det_func_1d, a, b, tol_lambda=1e-10):
    g = lambda x: abs(det_func_1d(x))
    x_star, g_star, err = golden_minimize(g, a, b, tol=tol_lambda)
    return x_star, g_star, err

def diagnose_crossing_vs_avoided(params: BeamParams,
                                 beta_deg: float,
                                 mu_grid: np.ndarray,
                                 lambdas_tr: np.ndarray,
                                 branch2_idx: int = 1,   # 0-based: branch 2
                                 branch3_idx: int = 2,   # 0-based: branch 3
                                 mu_half_width: float = 0.04,
                                 mu_step: float = 2e-4,
                                 lambda_window: float = 0.35,
                                 tol_lambda: float = 1e-10):
    """
    Uses the already tracked branches (so colors are consistent) to:
      - find mu* where gap is minimal
      - refine two roots by minimizing |det| in separated intervals
      - decide crossing vs avoided crossing
    """
    beta = np.deg2rad(beta_deg)
    eps = params.eps
    L = params.L_base
    fscale = params.f_scale
    K = fscale / (L**2)

    # frequencies from tracked branches
    f2 = (lambdas_tr[branch2_idx]**2) * K
    f3 = (lambdas_tr[branch3_idx]**2) * K
    gap = np.abs(f3 - f2)

    k0 = int(np.nanargmin(gap))
    mu_star = float(mu_grid[k0])
    print(f"[diagnose] coarse mu* ~ {mu_star:.6f}, min gap ~ {gap[k0]:.6g} Hz")

    mu_a = max(mu_grid[0], mu_star - mu_half_width)
    mu_b = min(mu_grid[-1], mu_star + mu_half_width)
    mu_zoom = np.arange(mu_a, mu_b + 1e-12, mu_step)

    # initialize lambdas from tracked solution at mu_star
    L2_prev = float(lambdas_tr[branch2_idx, k0])
    L3_prev = float(lambdas_tr[branch3_idx, k0])
    if not (np.isfinite(L2_prev) and np.isfinite(L3_prev) and L2_prev < L3_prev):
        raise RuntimeError("Bad initialization for refinement (NaN or wrong order).")

    L2_ref, L3_ref = [], []
    err2, err3 = [], []
    detabs2, detabs3 = [], []

    for mu in mu_zoom:
        mid = 0.5 * (L2_prev + L3_prev)

        a2 = max(0.05, L2_prev - lambda_window)
        b2 = max(a2 + 1e-6, mid)
        a3 = min(mid, L3_prev + lambda_window - 1e-6)
        b3 = L3_prev + lambda_window

        det_here = lambda x: det_clamped_coupled(x, beta, mu, eps)

        L2, g2, e2 = root_by_min_abs_det(det_here, a2, b2, tol_lambda=tol_lambda)
        L3, g3, e3 = root_by_min_abs_det(det_here, a3, b3, tol_lambda=tol_lambda)

        # enforce ordering
        if L2 > L3:
            L2, L3 = L3, L2
            g2, g3 = g3, g2
            e2, e3 = e3, e2

        L2_ref.append(L2); L3_ref.append(L3)
        detabs2.append(g2); detabs3.append(g3)
        err2.append(e2); err3.append(e3)

        L2_prev, L3_prev = L2, L3

    mu_zoom = np.array(mu_zoom)
    L2_ref = np.array(L2_ref); L3_ref = np.array(L3_ref)
    f2r = (L2_ref**2) * K
    f3r = (L3_ref**2) * K
    gapr = np.abs(f3r - f2r)

    # error estimate from lambda-bracket half-width: df ~ 2*K*Lambda*dLambda
    errf2 = 2 * K * L2_ref * np.array(err2)
    errf3 = 2 * K * L3_ref * np.array(err3)

    kk = int(np.argmin(gapr))
    mu_min_gap = float(mu_zoom[kk])
    gap_min = float(gapr[kk])
    gap_tol = float(10.0 * (errf2[kk] + errf3[kk]))  # conservative

    print(f"[diagnose] refined mu_min_gap = {mu_min_gap:.8f}")
    print(f"[diagnose] refined min gap    = {gap_min:.10g} Hz")
    print(f"[diagnose] numerical gap_tol  = {gap_tol:.10g} Hz")
    if gap_min <= gap_tol:
        verdict = "CROSSING (within numerical accuracy)"
    else:
        verdict = "AVOIDED CROSSING (nonzero minimum gap)"
    print("=>", verdict)

    # plots
    plt.figure(figsize=(10, 5.2))
    plt.plot(mu_zoom, f2r, label="refined branch 2")
    plt.plot(mu_zoom, f3r, label="refined branch 3")
    plt.axvline(mu_min_gap, linestyle="--", alpha=0.6)
    plt.xlabel("mu")
    plt.ylabel("Frequency (Hz)")
    plt.title(f"Refined branches 2/3 around mu* (beta={beta_deg}°)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

    plt.figure(figsize=(10, 4.8))
    plt.plot(mu_zoom, gapr, label="gap = |f3-f2| (refined)")
    plt.axhline(gap_tol, linestyle="--", alpha=0.7, label="gap_tol (10× numerical)")
    plt.axvline(mu_min_gap, linestyle="--", alpha=0.6)
    plt.xlabel("mu")
    plt.ylabel("Gap (Hz)")
    plt.title("Crossing vs avoided crossing diagnostic")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # optional: check if det minima are actually small
    plt.figure(figsize=(10, 4.4))
    plt.semilogy(mu_zoom, np.maximum(np.array(detabs2), 1e-300), label="|det| at refined branch 2")
    plt.semilogy(mu_zoom, np.maximum(np.array(detabs3), 1e-300), label="|det| at refined branch 3")
    plt.xlabel("mu")
    plt.ylabel("|det| (log scale)")
    plt.title("Sanity: how close we are to det=0 in refinement")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return verdict

def auto_crossing_diagnostics_close_pairs(params: BeamParams,
                                          beta_deg: float,
                                          mu_grid: np.ndarray,
                                          lambdas_tr: np.ndarray,
                                          gap_threshold_hz: float = 5.0,
                                          gap_rel_threshold: float | None = None,  # e.g. 0.01 for 1%
                                          max_pairs: int = 6,
                                          mu_half_width: float = 0.04,
                                          mu_step: float = 2e-4,
                                          lambda_window: float = 0.35,
                                          tol_lambda: float = 1e-10,
                                          make_plots: bool = True):
    """
    Runs refined avoided/crossing diagnostics ONLY for pairs that get 'close enough'
    on the coarse grid.

    Close-enough criteria:
      - abs gap <= gap_threshold_hz
      - optionally also abs gap <= gap_rel_threshold * min(freq_i, freq_j)
    Then keeps up to max_pairs smallest-gap pairs.
    """
    beta = np.deg2rad(beta_deg)
    eps = params.eps
    L = params.L_base
    fscale = params.f_scale
    K = fscale / (L**2)

    n_modes = lambdas_tr.shape[0]
    freqs_tr = (lambdas_tr**2) * K  # Hz

    candidates = []
    for i in range(n_modes):
        for j in range(i+1, n_modes):
            fi = freqs_tr[i]
            fj = freqs_tr[j]
            gap = np.abs(fj - fi)
            if np.all(~np.isfinite(gap)):
                continue
            k0 = int(np.nanargmin(gap))
            gap_min = float(gap[k0])
            mu_star = float(mu_grid[k0])

            # optional relative criterion
            ok = gap_min <= gap_threshold_hz
            if gap_rel_threshold is not None:
                fref = float(np.nanmin([fi[k0], fj[k0]]))
                if np.isfinite(fref) and fref > 0:
                    ok = ok and (gap_min <= gap_rel_threshold * fref)

            if ok:
                candidates.append({
                    "i": i, "j": j,
                    "mu_star": mu_star,
                    "gap_coarse": gap_min,
                    "k0": k0
                })

    candidates = sorted(candidates, key=lambda d: d["gap_coarse"])[:max_pairs]

    if len(candidates) == 0:
        print("[auto-diagnose] No close pairs found under the thresholds.")
        return []

    print("\n=== AUTO DIAGNOSTICS (close pairs, coarse) ===")
    for c in candidates:
        print(f"pair ({c['i']+1},{c['j']+1}): mu*~{c['mu_star']:.6f}, coarse gap~{c['gap_coarse']:.6g} Hz")

    results = []
    for c in candidates:
        verdict, info = diagnose_pair_refined(
            params=params,
            beta_deg=beta_deg,
            mu_grid=mu_grid,
            lambdas_tr=lambdas_tr,
            i=c["i"], j=c["j"],
            mu_half_width=mu_half_width,
            mu_step=mu_step,
            lambda_window=lambda_window,
            tol_lambda=tol_lambda,
            make_plots=make_plots
        )
        results.append({**c, **info, "verdict": verdict})

    print("\n=== AUTO DIAGNOSTICS (refined summary) ===")
    for r in results:
        print(f"pair ({r['i']+1},{r['j']+1}): "
              f"mu_min_gap={r['mu_min_gap']:.8f}, "
              f"min_gap={r['gap_min']:.6g} Hz, "
              f"gap_tol={r['gap_tol']:.3g} Hz -> {r['verdict']}")

    return results


def diagnose_pair_refined(params: BeamParams,
                          beta_deg: float,
                          mu_grid: np.ndarray,
                          lambdas_tr: np.ndarray,
                          i: int, j: int,
                          mu_half_width: float,
                          mu_step: float,
                          lambda_window: float,
                          tol_lambda: float,
                          make_plots: bool):
    """
    Refined decision for a specific tracked pair (i,j).
    Same logic as your branch(2,3) diagnose, but generalized.
    """
    beta = np.deg2rad(beta_deg)
    eps = params.eps
    L = params.L_base
    fscale = params.f_scale
    K = fscale / (L**2)

    fi = (lambdas_tr[i]**2) * K
    fj = (lambdas_tr[j]**2) * K
    gap = np.abs(fj - fi)

    k0 = int(np.nanargmin(gap))
    mu_star = float(mu_grid[k0])

    # zoom interval
    mu_a = max(mu_grid[0], mu_star - mu_half_width)
    mu_b = min(mu_grid[-1], mu_star + mu_half_width)
    mu_zoom = np.arange(mu_a, mu_b + 1e-12, mu_step)

    # init lambdas from tracked at mu_star
    Li_prev = float(lambdas_tr[i, k0])
    Lj_prev = float(lambdas_tr[j, k0])
    if not (np.isfinite(Li_prev) and np.isfinite(Lj_prev)):
        return "FAILED", {"mu_min_gap": np.nan, "gap_min": np.nan, "gap_tol": np.nan}

    # enforce order in Lambda (important for interval split)
    if Li_prev > Lj_prev:
        Li_prev, Lj_prev = Lj_prev, Li_prev

    Li_ref, Lj_ref = [], []
    erri, errj = [], []

    for mu in mu_zoom:
        mid = 0.5 * (Li_prev + Lj_prev)

        ai, bi = max(0.05, Li_prev - lambda_window), max(0.05, mid)
        aj, bj = min(mid, Lj_prev + lambda_window - 1e-6), Lj_prev + lambda_window

        det_here = lambda x: det_clamped_coupled(x, beta, mu, eps)

        Li, gi, ei = root_by_min_abs_det(det_here, ai, bi, tol_lambda=tol_lambda)
        Lj, gj, ej = root_by_min_abs_det(det_here, aj, bj, tol_lambda=tol_lambda)

        if Li > Lj:
            Li, Lj = Lj, Li
            ei, ej = ej, ei

        Li_ref.append(Li); Lj_ref.append(Lj)
        erri.append(ei); errj.append(ej)
        Li_prev, Lj_prev = Li, Lj

    mu_zoom = np.array(mu_zoom)
    Li_ref = np.array(Li_ref); Lj_ref = np.array(Lj_ref)
    fi_ref = (Li_ref**2) * K
    fj_ref = (Lj_ref**2) * K
    gap_ref = np.abs(fj_ref - fi_ref)

    # df error estimate
    errfi = 2 * K * Li_ref * np.array(erri)
    errfj = 2 * K * Lj_ref * np.array(errj)

    kk = int(np.argmin(gap_ref))
    mu_min_gap = float(mu_zoom[kk])
    gap_min = float(gap_ref[kk])
    gap_tol = float(10.0 * (errfi[kk] + errfj[kk]))

    verdict = "CROSSING" if gap_min <= gap_tol else "AVOIDED"

    if make_plots:
        plt.figure(figsize=(10, 5.2))
        plt.plot(mu_zoom, fi_ref, label=f"refined branch {i+1}")
        plt.plot(mu_zoom, fj_ref, label=f"refined branch {j+1}")
        plt.axvline(mu_min_gap, linestyle="--", alpha=0.6)
        plt.xlabel("mu")
        plt.ylabel("Frequency (Hz)")
        plt.title(f"Refined pair ({i+1},{j+1}) around mu* (beta={beta_deg}°)")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.show()

        plt.figure(figsize=(10, 4.8))
        plt.plot(mu_zoom, gap_ref, label="gap = |fj-fi|")
        plt.axhline(gap_tol, linestyle="--", alpha=0.7, label="gap_tol")
        plt.axvline(mu_min_gap, linestyle="--", alpha=0.6)
        plt.xlabel("mu")
        plt.ylabel("Gap (Hz)")
        plt.title(f"Gap diagnostic for pair ({i+1},{j+1})")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        plt.show()

    info = {"mu_min_gap": mu_min_gap, "gap_min": gap_min, "gap_tol": gap_tol}
    return verdict, info

# ============================================================
# 7) Main plot: frequencies vs mu at fixed beta
#    + optional crossing/avoided diagnostic (branches 2/3)
# ============================================================
def plot_freqs_vs_mu(params: BeamParams,
                     beta_deg: float = 5.0,
                     mu_min: float = 0.0,
                     mu_max: float = 0.9,
                     mu_step: float = 0.01,
                     n_modes: int = 6,
                     Lmax0: float = 55.0,
                     scan_step: float = 0.02,
                     show_ratio_lines: bool = True,
                     label_ratio_upto: int = 10,
                     do_crossing_diagnose: bool = True):
    beta = np.deg2rad(beta_deg)
    eps = params.eps
    L = params.L_base
    fscale = params.f_scale
    K = fscale / (L**2)

    mu_grid = np.arange(mu_min, mu_max + 1e-12, mu_step)

    # compute raw lambdas for each mu
    lambdas_raw = np.full((n_modes, len(mu_grid)), np.nan)
    for j, mu in enumerate(mu_grid):
        roots = find_first_n_roots(beta=beta, mu=mu, eps=eps,
                                   n_roots=n_modes, Lmax0=Lmax0, scan_step=scan_step)
        lambdas_raw[:, j] = roots

    # sort each column then track -> consistent colors even if order swaps
    lambdas_sorted = np.sort(lambdas_raw, axis=0)
    lambdas_tr = track_branches(lambdas_sorted)

    # frequencies for tracked branches
    f_tr = (lambdas_tr**2) * K

    # ratio marks
    mu_marks, labels = integer_ratio_marks(mu_max)

    # plot main
    plt.figure(figsize=(11, 6.5))
    for i in range(n_modes):
        plt.plot(mu_grid, f_tr[i], label=f"branch {i+1}")

    if show_ratio_lines:
        for mu_m in mu_marks:
            plt.axvline(mu_m, linestyle=":", linewidth=1.0, alpha=0.45)

        # label a few to avoid clutter
        for k, mu_m in enumerate(mu_marks[:label_ratio_upto]):
            plt.text(mu_m, plt.ylim()[0], labels[k], rotation=90, va="bottom", ha="right", fontsize=8)
        if len(mu_marks) > label_ratio_upto:
            plt.text(mu_marks[-1], plt.ylim()[0], labels[-1], rotation=90, va="bottom", ha="right", fontsize=8)

    plt.xlabel("μ  (L1=1-μ, L2=1+μ, L1+L2=2)")
    plt.ylabel("Frequency (Hz)")
    plt.title(f"First {n_modes} frequencies vs μ at fixed β={beta_deg}° (clamped ends)\n"
              "Branch colors are tracked (stay the same even if branches change order)")
    plt.grid(True, alpha=0.3)
    plt.legend(ncols=2, fontsize=9)
    plt.tight_layout()
    plt.show()

    # optionally run crossing/avoided diagnosis for branches 2/3
    if do_crossing_diagnose:
        auto_crossing_diagnostics_close_pairs(
            params=params,
            beta_deg=beta_deg,
            mu_grid=mu_grid,
            lambdas_tr=lambdas_tr,
            gap_threshold_hz=5.0,  # считаем "близко", если меньше 5 Гц
            gap_rel_threshold=None,  # или поставь 0.01 для 1% одновременно
            max_pairs=6,  # не больше 6 пар на уточнение
            make_plots=True
        )


# ============================================================
# Example run (your same parameters)
# ============================================================
if __name__ == "__main__":
    p = BeamParams(
        E=2.1e11,
        rho=7800.0,
        r=0.005,     # eps=0.0025 when L_total=2
        L_total=2.0
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
        do_crossing_diagnose=True
    )