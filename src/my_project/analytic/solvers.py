import mpmath as mp
import numpy as np

try:
    from .formulas import det_clamped_coupled, frequency_scale
except ImportError:
    from formulas import det_clamped_coupled, frequency_scale


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
    grid = np.arange(Lmin, Lmax + scan_step, scan_step)
    vals = np.array([det_clamped_coupled(L, beta, mu, eps) for L in grid], dtype=float)

    roots = []

    def sgn(x: float) -> int:
        return 1 if x > 0 else (-1 if x < 0 else 0)

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

    roots = sorted(set(float(r) for r in roots))
    return roots[:n_roots]


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
    Lmax = Lmax0
    roots = []
    for _ in range(max_tries):
        roots = find_roots_scan_bisect(beta, mu, eps, n_roots, Lmin, Lmax, scan_step)
        if len(roots) >= n_roots:
            break
        Lmax *= grow_factor

    out = np.full(n_roots, np.nan, dtype=float)
    out[: len(roots)] = roots
    return out


def track_branches(values_sorted: np.ndarray) -> np.ndarray:
    N, M = values_sorted.shape
    tracked = np.full_like(values_sorted, np.nan)
    tracked[:, 0] = values_sorted[:, 0]

    try:
        from scipy.optimize import linear_sum_assignment

        use_hungarian = True
    except Exception:
        use_hungarian = False

    for j in range(1, M):
        prev = tracked[:, j - 1]
        cur = values_sorted[:, j]

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


def fixed_fixed_lambdas(n: int) -> np.ndarray:
    mp.mp.dps = 60

    def f(x):
        return mp.cosh(x) * mp.cos(x) - 1

    roots = []
    for m in range(1, n + 1):
        x0 = (m + 0.5) * mp.pi
        roots.append(float(mp.findroot(f, x0)))
    return np.array(roots, dtype=float)


def golden_minimize(g, a: float, b: float, tol: float = 1e-10, max_iter: int = 200):
    phi = (1 + 5**0.5) / 2
    invphi = 1 / phi
    c = b - (b - a) * invphi
    d = a + (b - a) * invphi
    gc = g(c)
    gd = g(d)
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


def root_by_min_abs_det(det_func_1d, a: float, b: float, tol_lambda: float = 1e-10):
    g = lambda x: abs(det_func_1d(x))
    x_star, g_star, err = golden_minimize(g, a, b, tol=tol_lambda)
    return x_star, g_star, err


def find_close_pair_candidates(
    freqs_tr: np.ndarray,
    mu_grid: np.ndarray,
    gap_threshold_hz: float = 5.0,
    gap_rel_threshold: float | None = None,
    max_pairs: int = 6,
) -> list[dict[str, float | int]]:
    n_modes = freqs_tr.shape[0]
    candidates = []
    for i in range(n_modes):
        for j in range(i + 1, n_modes):
            fi = freqs_tr[i]
            fj = freqs_tr[j]
            gap = np.abs(fj - fi)
            if np.all(~np.isfinite(gap)):
                continue
            k0 = int(np.nanargmin(gap))
            gap_min = float(gap[k0])
            mu_star = float(mu_grid[k0])

            ok = gap_min <= gap_threshold_hz
            if gap_rel_threshold is not None:
                fref = float(np.nanmin([fi[k0], fj[k0]]))
                if np.isfinite(fref) and fref > 0:
                    ok = ok and (gap_min <= gap_rel_threshold * fref)

            if ok:
                candidates.append(
                    {
                        "i": i,
                        "j": j,
                        "mu_star": mu_star,
                        "gap_coarse": gap_min,
                        "k0": k0,
                    }
                )

    return sorted(candidates, key=lambda d: d["gap_coarse"])[:max_pairs]


def refine_tracked_pair(
    params,
    beta_deg: float,
    mu_grid: np.ndarray,
    lambdas_tr: np.ndarray,
    i: int,
    j: int,
    mu_half_width: float,
    mu_step: float,
    lambda_window: float,
    tol_lambda: float,
) -> dict[str, object]:
    beta = np.deg2rad(beta_deg)
    eps = params.eps
    K = frequency_scale(params)

    fi = (lambdas_tr[i] ** 2) * K
    fj = (lambdas_tr[j] ** 2) * K
    gap = np.abs(fj - fi)

    empty = np.array([], dtype=float)
    if np.all(~np.isfinite(gap)):
        return {
            "i": i,
            "j": j,
            "mu_star": np.nan,
            "mu_zoom": empty,
            "li_ref": empty,
            "lj_ref": empty,
            "fi_ref": empty,
            "fj_ref": empty,
            "gap_ref": empty,
            "detabs_i": empty,
            "detabs_j": empty,
            "erri": empty,
            "errj": empty,
            "mu_min_gap": np.nan,
            "gap_min": np.nan,
            "gap_tol": np.nan,
            "verdict": "FAILED",
            "verbose_verdict": "FAILED",
        }

    k0 = int(np.nanargmin(gap))
    mu_star = float(mu_grid[k0])

    mu_a = max(mu_grid[0], mu_star - mu_half_width)
    mu_b = min(mu_grid[-1], mu_star + mu_half_width)
    mu_zoom = np.arange(mu_a, mu_b + 1e-12, mu_step)

    Li_prev = float(lambdas_tr[i, k0])
    Lj_prev = float(lambdas_tr[j, k0])
    if not (np.isfinite(Li_prev) and np.isfinite(Lj_prev)):
        return {
            "i": i,
            "j": j,
            "mu_star": mu_star,
            "mu_zoom": mu_zoom,
            "li_ref": empty,
            "lj_ref": empty,
            "fi_ref": empty,
            "fj_ref": empty,
            "gap_ref": empty,
            "detabs_i": empty,
            "detabs_j": empty,
            "erri": empty,
            "errj": empty,
            "mu_min_gap": np.nan,
            "gap_min": np.nan,
            "gap_tol": np.nan,
            "verdict": "FAILED",
            "verbose_verdict": "FAILED",
        }

    if Li_prev > Lj_prev:
        Li_prev, Lj_prev = Lj_prev, Li_prev

    Li_ref, Lj_ref = [], []
    detabs_i, detabs_j = [], []
    erri, errj = [], []

    for mu in mu_zoom:
        mid = 0.5 * (Li_prev + Lj_prev)

        ai = max(0.05, Li_prev - lambda_window)
        bi = max(ai + 1e-6, mid)
        aj = min(mid, Lj_prev + lambda_window - 1e-6)
        bj = Lj_prev + lambda_window

        det_here = lambda x: det_clamped_coupled(x, beta, mu, eps)

        Li, gi, ei = root_by_min_abs_det(det_here, ai, bi, tol_lambda=tol_lambda)
        Lj, gj, ej = root_by_min_abs_det(det_here, aj, bj, tol_lambda=tol_lambda)

        if Li > Lj:
            Li, Lj = Lj, Li
            gi, gj = gj, gi
            ei, ej = ej, ei

        Li_ref.append(Li)
        Lj_ref.append(Lj)
        detabs_i.append(gi)
        detabs_j.append(gj)
        erri.append(ei)
        errj.append(ej)
        Li_prev, Lj_prev = Li, Lj

    Li_ref = np.array(Li_ref)
    Lj_ref = np.array(Lj_ref)
    detabs_i = np.array(detabs_i)
    detabs_j = np.array(detabs_j)
    erri = np.array(erri)
    errj = np.array(errj)

    fi_ref = (Li_ref**2) * K
    fj_ref = (Lj_ref**2) * K
    gap_ref = np.abs(fj_ref - fi_ref)

    errfi = 2 * K * Li_ref * erri
    errfj = 2 * K * Lj_ref * errj

    kk = int(np.argmin(gap_ref))
    mu_min_gap = float(mu_zoom[kk])
    gap_min = float(gap_ref[kk])
    gap_tol = float(10.0 * (errfi[kk] + errfj[kk]))

    is_crossing = gap_min <= gap_tol
    verdict = "CROSSING" if is_crossing else "AVOIDED"
    verbose_verdict = (
        "CROSSING (within numerical accuracy)"
        if is_crossing
        else "AVOIDED CROSSING (nonzero minimum gap)"
    )

    return {
        "i": i,
        "j": j,
        "mu_star": mu_star,
        "mu_zoom": mu_zoom,
        "li_ref": Li_ref,
        "lj_ref": Lj_ref,
        "fi_ref": fi_ref,
        "fj_ref": fj_ref,
        "gap_ref": gap_ref,
        "detabs_i": detabs_i,
        "detabs_j": detabs_j,
        "erri": erri,
        "errj": errj,
        "mu_min_gap": mu_min_gap,
        "gap_min": gap_min,
        "gap_tol": gap_tol,
        "verdict": verdict,
        "verbose_verdict": verbose_verdict,
    }


__all__ = [
    "find_roots_scan_bisect",
    "find_first_n_roots",
    "track_branches",
    "fixed_fixed_lambdas",
    "golden_minimize",
    "root_by_min_abs_det",
    "find_close_pair_candidates",
    "refine_tracked_pair",
]
