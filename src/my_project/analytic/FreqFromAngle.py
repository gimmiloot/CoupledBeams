import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
import mpmath as mp

# =========================================
# Parameters
# =========================================
@dataclass
class BeamParams:
    E: float
    rho: float
    r: float
    L_total: float  # fixed, L1+L2

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
        return np.sqrt(self.E * self.I / (self.rho * self.S)) / (2.0 * np.pi)


# =========================================
# Determinant (your clamped matrix)
# =========================================
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


# =========================================
# Root finding: scan + bisection + auto-grow
# =========================================
def find_roots_scan_bisect(beta: float, mu: float, eps: float,
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
        a, b = grid[i], grid[i + 1]
        fa, fb = vals[i], vals[i + 1]
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
                       Lmax0: float = 35.0,
                       scan_step: float = 0.03,
                       grow_factor: float = 1.35,
                       max_tries: int = 7) -> np.ndarray:
    Lmax = Lmax0
    roots = []
    for _ in range(max_tries):
        roots = find_roots_scan_bisect(beta, mu, eps, n_roots, Lmin, Lmax, scan_step)
        if len(roots) >= n_roots:
            break
        Lmax *= grow_factor

    out = np.full(n_roots, np.nan, dtype=float)
    out[:len(roots)] = roots
    return out


# =========================================
# Branch tracking
# =========================================
def track_branches(freqs_raw: np.ndarray) -> np.ndarray:
    N, M = freqs_raw.shape
    tracked = np.full_like(freqs_raw, np.nan)
    tracked[:, 0] = freqs_raw[:, 0]

    try:
        from scipy.optimize import linear_sum_assignment
        use_hungarian = True
    except Exception:
        use_hungarian = False

    for j in range(1, M):
        prev = tracked[:, j - 1]
        cur = freqs_raw[:, j]

        prev_idx = np.where(np.isfinite(prev))[0]
        cur_idx = np.where(np.isfinite(cur))[0]
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
                min_val = np.inf
                min_pair = None
                for rr in range(C.shape[0]):
                    if rr in used_r:
                        continue
                    for cc in range(C.shape[1]):
                        if cc in used_c:
                            continue
                        if C[rr, cc] < min_val:
                            min_val = C[rr, cc]
                            min_pair = (rr, cc)
                if min_pair is None:
                    break
                rr, cc = min_pair
                used_r.add(rr)
                used_c.add(cc)
                assign[prev_idx[rr]] = cur_idx[cc]
                if len(used_r) == C.shape[0] or len(used_c) == C.shape[1]:
                    break

        for i in range(N):
            tracked[i, j] = cur[assign[i]] if i in assign else np.nan

    return tracked


# =========================================
# Long rod reference (fixed-fixed), bending only
# =========================================
def fixed_fixed_lambdas(n: int) -> np.ndarray:
    mp.mp.dps = 60

    def f(x): return mp.cosh(x) * mp.cos(x) - 1

    roots = []
    for m in range(1, n + 1):
        x0 = (m + 0.5) * mp.pi
        roots.append(float(mp.findroot(f, x0)))
    return np.array(roots, dtype=float)


# =========================================
# MAIN plotting function with show_dashed switch
# =========================================
def plot_vs_beta_for_mu(params: BeamParams,
                        mu: float,
                        beta_deg_grid=None,
                        n_modes_main: int = 10,
                        n_modes_small: int = 6,
                        n_dashed_main: int = 5,
                        n_dashed_small: int = 3,
                        show_dashed: bool = True,   # <-- переключатель
                        Lmax0: float = 35.0,
                        scan_step: float = 0.03):
    if beta_deg_grid is None:
        beta_deg_grid = np.arange(0, 91, 1)  # 0..90
    beta_rad_grid = np.deg2rad(beta_deg_grid)

    eps = params.eps
    fscale = params.f_scale
    Lbase = params.L_base

    # lengths (fixed total)
    L1 = Lbase * (1.0 - mu)
    L2 = Lbase * (1.0 + mu)
    L_long = max(L1, L2)

    # dashed reference frequencies (compute only if needed)
    if show_dashed:
        lam_d = fixed_fixed_lambdas(n_dashed_main)
        f_long = (lam_d**2 / (L_long**2)) * fscale

    # compute raw coupled frequencies
    freqs_raw = np.full((n_modes_main, len(beta_rad_grid)), np.nan, dtype=float)
    for j, beta in enumerate(beta_rad_grid):
        roots = find_first_n_roots(beta, mu, eps,
                                   n_roots=n_modes_main,
                                   Lmin=0.2,
                                   Lmax0=Lmax0,
                                   scan_step=scan_step)
        freqs_raw[:, j] = (roots**2 / (Lbase**2)) * fscale

    # sort per beta, then track branches
    freqs_sorted = np.sort(freqs_raw, axis=0)
    freqs_tr = track_branches(freqs_sorted)

    # ---------- Plot 1: first 10 ----------
    plt.figure(figsize=(10.5, 6.2))
    for i in range(n_modes_main):
        plt.plot(beta_deg_grid, freqs_tr[i], label=f"branch {i+1}")

    if show_dashed:
        for i in range(n_dashed_main):
            plt.hlines(f_long[i], beta_deg_grid[0], beta_deg_grid[-1],
                       linestyles="--", linewidth=1.2,
                       label="long rod (fixed-fixed)" if i == 0 else None)

    plt.xlabel("β (degrees)")
    plt.ylabel("Frequency (Hz)")
    title = f"Coupled (clamped): first {n_modes_main} branches vs β (0..90°), mu={mu:+.3f}"
    if show_dashed:
        title += f"\nDashed: single longer rod fixed–fixed (first {n_dashed_main} modes)"
    plt.title(title)
    plt.grid(True, alpha=0.3)
    #plt.legend(ncols=2, fontsize=9)
    plt.tight_layout()
    plt.show()

    # ---------- Plot 2: first 4 ----------
    plt.figure(figsize=(10.5, 6.2))
    for i in range(n_modes_small):
        plt.plot(beta_deg_grid, freqs_tr[i], label=f"branch {i+1}")

    if show_dashed:
        for i in range(min(n_dashed_small, n_dashed_main)):
            plt.hlines(f_long[i], beta_deg_grid[0], beta_deg_grid[-1],
                       linestyles="--", linewidth=1.3,
                       label="long rod (fixed-fixed)" if i == 0 else None)

    plt.xlabel("β (degrees)")
    plt.ylabel("Frequency (Hz)")
    title = f"Coupled (clamped): first {n_modes_small} branches vs β (0..90°), mu={mu:+.3f}"
    if show_dashed:
        title += f"\nDashed: long rod fixed–fixed (first {n_dashed_small} modes)"
    plt.title(title)
    plt.grid(True, alpha=0.3)
    #plt.legend(fontsize=10)
    plt.tight_layout()
    plt.show()


# =========================================
# Example
# =========================================
if __name__ == "__main__":
    p = BeamParams(E=2.1e11, rho=7800.0, r=0.005, L_total=2.0)
    print("eps =", p.eps)  # 0.0025 expected

    for mu in [0]:
        plot_vs_beta_for_mu(
            params=p,
            mu=mu,
            beta_deg_grid=np.arange(0, 91, 1),
            show_dashed=True,   # <-- False для отключения пунктиров
            n_dashed_main=5,
            n_dashed_small=3,
            scan_step=0.03,
            Lmax0=35.0
        )