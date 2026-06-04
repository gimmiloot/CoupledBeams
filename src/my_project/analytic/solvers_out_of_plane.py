from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

import numpy as np

from .formulas_out_of_plane import det_out_of_plane


@dataclass(frozen=True)
class OutOfPlaneRootSearchResult:
    """Root-search result with diagnostic warnings."""

    roots: list[float]
    warnings: list[str]
    lambda_min: float
    lambda_max: float
    scan_step: float


def _positive_float(name: str, value: float) -> float:
    value_f = float(value)
    if not isfinite(value_f) or value_f <= 0.0:
        raise ValueError(f"{name} must be positive and finite.")
    return value_f


def _positive_int(name: str, value: int) -> int:
    value_i = int(value)
    if value_i <= 0:
        raise ValueError(f"{name} must be positive.")
    return value_i


def _sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def _det_value(
    Lambda: float,
    beta: float,
    mu: float,
    epsilon: float,
    eta: float,
    poisson: float,
) -> float:
    return float(
        det_out_of_plane(
            Lambda=float(Lambda),
            beta=float(beta),
            mu=float(mu),
            epsilon=float(epsilon),
            eta=float(eta),
            poisson=float(poisson),
        )
    )


def find_roots_scan_bisect_out_of_plane(
    beta: float,
    mu: float,
    epsilon: float,
    eta: float,
    poisson: float,
    n_roots: int,
    lambda_min: float,
    lambda_max: float,
    scan_step: float,
    bisect_iters: int = 70,
) -> tuple[list[float], list[str]]:
    """Find determinant sign-change roots by scanning and bisection."""

    n_roots_i = _positive_int("n_roots", n_roots)
    lambda_min_f = _positive_float("lambda_min", lambda_min)
    lambda_max_f = _positive_float("lambda_max", lambda_max)
    scan_step_f = _positive_float("scan_step", scan_step)
    if lambda_max_f <= lambda_min_f:
        raise ValueError("lambda_max must be greater than lambda_min.")

    grid = np.arange(lambda_min_f, lambda_max_f + 0.5 * scan_step_f, scan_step_f, dtype=float)
    values = np.array(
        [_det_value(value, beta, mu, epsilon, eta, poisson) for value in grid],
        dtype=float,
    )
    warnings: list[str] = []
    nonfinite = int(np.count_nonzero(~np.isfinite(values)))
    if nonfinite:
        warnings.append(f"{nonfinite} non-finite determinant samples were skipped.")

    roots: list[float] = []
    for index in range(len(grid) - 1):
        left = float(grid[index])
        right = float(grid[index + 1])
        f_left = float(values[index])
        f_right = float(values[index + 1])
        if not (isfinite(f_left) and isfinite(f_right)):
            continue
        s_left = _sign(f_left)
        s_right = _sign(f_right)
        if s_left == 0:
            roots.append(left)
        elif s_left * s_right < 0:
            a = left
            b = right
            fa = f_left
            for _ in range(int(bisect_iters)):
                mid = 0.5 * (a + b)
                fm = _det_value(mid, beta, mu, epsilon, eta, poisson)
                s_mid = _sign(fm)
                if s_mid == 0:
                    a = b = mid
                    break
                if _sign(fa) * s_mid < 0:
                    b = mid
                else:
                    a = mid
                    fa = fm
            roots.append(0.5 * (a + b))

        if len(roots) >= n_roots_i:
            break

    deduped: list[float] = []
    for root in sorted(float(root) for root in roots):
        if root <= 0.0:
            continue
        if not deduped or abs(root - deduped[-1]) > max(1e-8, 0.25 * scan_step_f):
            deduped.append(root)
        if len(deduped) >= n_roots_i:
            break

    if len(deduped) < n_roots_i:
        warnings.append(
            f"found {len(deduped)} roots, fewer than requested {n_roots_i}, "
            f"on Lambda in [{lambda_min_f:g}, {lambda_max_f:g}]"
        )
    return deduped, warnings


def find_first_n_roots_out_of_plane_with_warnings(
    beta: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
    n_roots: int = 8,
    lambda_max: float = 25.0,
    lambda_min: float = 0.05,
    scan_step: float = 0.01,
    bisect_iters: int = 70,
) -> OutOfPlaneRootSearchResult:
    """Return sorted positive out-of-plane roots plus root-search warnings."""

    roots, warnings = find_roots_scan_bisect_out_of_plane(
        beta=beta,
        mu=mu,
        epsilon=epsilon,
        eta=eta,
        poisson=poisson,
        n_roots=n_roots,
        lambda_min=lambda_min,
        lambda_max=lambda_max,
        scan_step=scan_step,
        bisect_iters=bisect_iters,
    )
    return OutOfPlaneRootSearchResult(
        roots=roots,
        warnings=warnings,
        lambda_min=float(lambda_min),
        lambda_max=float(lambda_max),
        scan_step=float(scan_step),
    )


def find_first_n_roots_out_of_plane(
    beta: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
    n_roots: int = 8,
    lambda_max: float = 25.0,
) -> list[float]:
    """Return the first sorted positive roots of the out-of-plane determinant."""

    return find_first_n_roots_out_of_plane_with_warnings(
        beta=beta,
        mu=mu,
        epsilon=epsilon,
        eta=eta,
        poisson=poisson,
        n_roots=n_roots,
        lambda_max=lambda_max,
    ).roots


__all__ = [
    "OutOfPlaneRootSearchResult",
    "find_first_n_roots_out_of_plane",
    "find_first_n_roots_out_of_plane_with_warnings",
    "find_roots_scan_bisect_out_of_plane",
]
