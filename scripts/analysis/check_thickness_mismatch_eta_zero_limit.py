from __future__ import annotations

import csv
from pathlib import Path
import sys
from typing import Iterable, Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import det_clamped_coupled  # noqa: E402
from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    det_eta,
    find_first_n_roots_eta,
    thickness_mismatch_factors,
)
from my_project.analytic.solvers import find_first_n_roots  # noqa: E402


RESULTS_DIR = REPO_ROOT / "results"
ETA_ZERO_ROOTS_CSV = RESULTS_DIR / "thickness_mismatch_eta_zero_roots_check.csv"
SWAP_SYMMETRY_CSV = RESULTS_DIR / "thickness_mismatch_swap_symmetry_check.csv"

ROOT_FIELDNAMES = [
    "beta_deg",
    "mu",
    "epsilon",
    "eta",
    "root_index",
    "Lambda_old",
    "Lambda_new",
    "abs_diff",
    "rel_diff",
]
SWAP_FIELDNAMES = [
    "beta_deg",
    "mu",
    "eta",
    "epsilon",
    "root_index",
    "Lambda_mu_eta",
    "Lambda_minus_mu_minus_eta",
    "abs_diff",
    "rel_diff",
]

EPSILON = 0.0025
NUM_ROOTS = 8
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0


def relative_diff(left: float, right: float) -> float:
    denom = max(abs(float(left)), abs(float(right)), 1e-15)
    return abs(float(left) - float(right)) / denom


def write_csv(path: Path, rows: Iterable[dict[str, float | int]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def max_mass_conservation_error() -> float:
    errors: list[float] = []
    for mu in np.linspace(-0.9, 0.9, 19):
        for eta in np.linspace(-0.5, 0.5, 21):
            factors = thickness_mismatch_factors(float(mu), float(eta))
            errors.append(abs(factors.mass_factor - 2.0))
    return float(max(errors))


def determinant_eta_zero_errors() -> tuple[float, float]:
    abs_errors: list[float] = []
    rel_errors: list[float] = []
    for beta_deg in (0.0, 5.0, 15.0, 45.0):
        beta = float(np.deg2rad(beta_deg))
        for mu in (0.0, 0.3, 0.6):
            for epsilon in (0.0025, 0.01):
                for Lambda in (0.5, 1.7, 3.9, 6.2, 10.5, 15.0):
                    old = det_clamped_coupled(float(Lambda), beta, float(mu), float(epsilon))
                    new = det_eta(float(Lambda), beta, float(mu), float(epsilon), 0.0)
                    abs_errors.append(abs(new - old))
                    rel_errors.append(relative_diff(new, old))
    return float(max(abs_errors)), float(max(rel_errors))


def eta_zero_root_rows() -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for beta_deg in (5.0, 15.0, 45.0):
        beta = float(np.deg2rad(beta_deg))
        for mu in (0.0, 0.3, 0.6):
            old_roots = find_first_n_roots(
                beta,
                float(mu),
                EPSILON,
                NUM_ROOTS,
                Lmax0=ROOT_LMAX0,
                scan_step=ROOT_SCAN_STEP,
            )
            new_roots = find_first_n_roots_eta(
                beta,
                float(mu),
                EPSILON,
                0.0,
                NUM_ROOTS,
                Lmax0=ROOT_LMAX0,
                scan_step=ROOT_SCAN_STEP,
            )
            if np.any(~np.isfinite(old_roots)) or np.any(~np.isfinite(new_roots)):
                raise RuntimeError(f"Missing eta=0 roots for beta={beta_deg:g}, mu={mu:g}.")
            for idx, (old, new) in enumerate(zip(old_roots, new_roots, strict=True), start=1):
                rows.append(
                    {
                        "beta_deg": float(beta_deg),
                        "mu": float(mu),
                        "epsilon": EPSILON,
                        "eta": 0.0,
                        "root_index": int(idx),
                        "Lambda_old": float(old),
                        "Lambda_new": float(new),
                        "abs_diff": abs(float(new) - float(old)),
                        "rel_diff": relative_diff(float(new), float(old)),
                    }
                )
    return rows


def swap_symmetry_rows() -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for beta_deg in (5.0, 15.0, 45.0):
        beta = float(np.deg2rad(beta_deg))
        for mu in (0.0, 0.3, 0.6):
            for eta in (-0.15, -0.1, 0.1, 0.15):
                roots = find_first_n_roots_eta(
                    beta,
                    float(mu),
                    EPSILON,
                    float(eta),
                    NUM_ROOTS,
                    Lmax0=ROOT_LMAX0,
                    scan_step=ROOT_SCAN_STEP,
                )
                swapped = find_first_n_roots_eta(
                    beta,
                    -float(mu),
                    EPSILON,
                    -float(eta),
                    NUM_ROOTS,
                    Lmax0=ROOT_LMAX0,
                    scan_step=ROOT_SCAN_STEP,
                )
                if np.any(~np.isfinite(roots)) or np.any(~np.isfinite(swapped)):
                    raise RuntimeError(f"Missing swap roots for beta={beta_deg:g}, mu={mu:g}, eta={eta:g}.")
                for idx, (left, right) in enumerate(zip(roots, swapped, strict=True), start=1):
                    rows.append(
                        {
                            "beta_deg": float(beta_deg),
                            "mu": float(mu),
                            "eta": float(eta),
                            "epsilon": EPSILON,
                            "root_index": int(idx),
                            "Lambda_mu_eta": float(left),
                            "Lambda_minus_mu_minus_eta": float(right),
                            "abs_diff": abs(float(left) - float(right)),
                            "rel_diff": relative_diff(float(left), float(right)),
                        }
                    )
    return rows


def max_field(rows: Sequence[dict[str, float | int]], field: str) -> float:
    return float(max(float(row[field]) for row in rows))


def main() -> dict[str, float | Path]:
    mass_error = max_mass_conservation_error()
    det_abs_error, det_rel_error = determinant_eta_zero_errors()

    root_rows = eta_zero_root_rows()
    swap_rows = swap_symmetry_rows()
    write_csv(ETA_ZERO_ROOTS_CSV, root_rows, ROOT_FIELDNAMES)
    write_csv(SWAP_SYMMETRY_CSV, swap_rows, SWAP_FIELDNAMES)

    eta_zero_root_abs = max_field(root_rows, "abs_diff")
    eta_zero_root_rel = max_field(root_rows, "rel_diff")
    swap_abs = max_field(swap_rows, "abs_diff")
    swap_rel = max_field(swap_rows, "rel_diff")

    print(f"saved eta=0 roots check CSV: {ETA_ZERO_ROOTS_CSV}")
    print(f"saved swap symmetry check CSV: {SWAP_SYMMETRY_CSV}")
    print(f"max mass conservation error: {mass_error:.6e}")
    print(f"max eta=0 determinant abs diff: {det_abs_error:.6e}")
    print(f"max eta=0 determinant rel diff: {det_rel_error:.6e}")
    print(f"max eta=0 root abs diff: {eta_zero_root_abs:.6e}")
    print(f"max eta=0 root rel diff: {eta_zero_root_rel:.6e}")
    print(f"max swap root abs diff: {swap_abs:.6e}")
    print(f"max swap root rel diff: {swap_rel:.6e}")

    if mass_error > 5e-13:
        raise RuntimeError("Mass conservation check exceeded tolerance.")
    if det_abs_error > 1e-8 and det_rel_error > 1e-12:
        raise RuntimeError("eta=0 determinant comparison exceeded tolerance.")
    if eta_zero_root_abs > 1e-8:
        raise RuntimeError("eta=0 root comparison exceeded tolerance.")
    if swap_abs > 1e-6:
        raise RuntimeError("swap symmetry root comparison exceeded tolerance.")

    return {
        "eta_zero_roots_csv": ETA_ZERO_ROOTS_CSV,
        "swap_symmetry_csv": SWAP_SYMMETRY_CSV,
        "max_mass_error": mass_error,
        "max_eta_zero_determinant_abs_diff": det_abs_error,
        "max_eta_zero_determinant_rel_diff": det_rel_error,
        "max_eta_zero_root_abs_diff": eta_zero_root_abs,
        "max_eta_zero_root_rel_diff": eta_zero_root_rel,
        "max_swap_root_abs_diff": swap_abs,
        "max_swap_root_rel_diff": swap_rel,
    }


if __name__ == "__main__":
    main()
