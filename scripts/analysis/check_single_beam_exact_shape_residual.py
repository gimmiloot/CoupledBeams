"""Sanity-check FEM residuals on exact single-beam and axial-bar modes.

This diagnostic intentionally stays outside the coupled-rods analytic model.
It uses the element matrices from ``src/my_project/fem/python_fem.py`` and
checks whether exact continuous single-member shapes give convergent FEM
residuals when sampled on the nodes.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.fem import python_fem as fem  # noqa: E402


RESULTS_DIR = REPO_ROOT / "results"
N_ELEM_VALUES = (20, 40, 80, 160)
L_TOTAL = 1.0

CLAMPED_FREE_BENDING_ALPHA = 1.875104068711961
CLAMPED_FREE_AXIAL_ALPHA = 0.5 * np.pi

CSV_FIELDS = [
    "test_name",
    "n_elem",
    "boundary_condition",
    "exact_alpha",
    "exact_omega_sq",
    "rayleigh_omega_sq",
    "relative_rayleigh_diff",
    "relative_residual",
    "residual_inf",
    "notes",
]


def _relative_diff(value: float, reference: float) -> float:
    denom = max(abs(reference), np.finfo(float).eps)
    return float(abs(value - reference) / denom)


def _mac_euclidean(a: np.ndarray, b: np.ndarray) -> float:
    numerator = float(np.dot(a, b) ** 2)
    denominator = float(np.dot(a, a) * np.dot(b, b))
    if denominator <= 0.0:
        return float("nan")
    return numerator / denominator


def _residual_metrics(
    K: np.ndarray,
    M: np.ndarray,
    q: np.ndarray,
    omega_sq: float,
) -> tuple[float, float, float]:
    kq = K @ q
    mq = M @ q
    residual = kq - omega_sq * mq
    denom = np.linalg.norm(kq) + abs(omega_sq) * np.linalg.norm(mq)
    if denom <= 0.0:
        relative_residual = float("nan")
    else:
        relative_residual = float(np.linalg.norm(residual) / denom)
    residual_inf = float(np.max(np.abs(residual))) if residual.size else 0.0
    qmq = float(q @ mq)
    rayleigh_omega_sq = float((q @ kq) / qmq) if abs(qmq) > 0.0 else float("nan")
    return relative_residual, residual_inf, rayleigh_omega_sq


def _first_fem_mode_metrics(
    K: np.ndarray,
    M: np.ndarray,
    q_exact: np.ndarray,
    exact_omega_sq: float,
) -> tuple[float, float, float]:
    eigvals, eigvecs = eigh(K, M, subset_by_index=[0, 0])
    fem_omega_sq = float(eigvals[0])
    fem_vec = np.asarray(eigvecs[:, 0], dtype=float)
    mac = _mac_euclidean(q_exact, fem_vec)
    rel = _relative_diff(fem_omega_sq, exact_omega_sq)
    return fem_omega_sq, rel, mac


def _assemble_bending_matrices(n_elem: int, length: float) -> tuple[np.ndarray, np.ndarray]:
    n_nodes = n_elem + 1
    ndof = 2 * n_nodes
    K = np.zeros((ndof, ndof), dtype=float)
    M = np.zeros((ndof, ndof), dtype=float)
    le = length / n_elem
    bending_idx = np.array([1, 2, 4, 5], dtype=int)
    Ke = fem.elem_K(le)[np.ix_(bending_idx, bending_idx)]
    Me = fem.elem_M(le)[np.ix_(bending_idx, bending_idx)]

    for elem_index in range(n_elem):
        dofs = np.array(
            [
                2 * elem_index,
                2 * elem_index + 1,
                2 * (elem_index + 1),
                2 * (elem_index + 1) + 1,
            ],
            dtype=int,
        )
        K[np.ix_(dofs, dofs)] += Ke
        M[np.ix_(dofs, dofs)] += Me
    return K, M


def _assemble_axial_matrices(n_elem: int, length: float) -> tuple[np.ndarray, np.ndarray]:
    n_nodes = n_elem + 1
    K = np.zeros((n_nodes, n_nodes), dtype=float)
    M = np.zeros((n_nodes, n_nodes), dtype=float)
    le = length / n_elem
    axial_idx = np.array([0, 3], dtype=int)
    Ke = fem.elem_K(le)[np.ix_(axial_idx, axial_idx)]
    Me = fem.elem_M(le)[np.ix_(axial_idx, axial_idx)]

    for elem_index in range(n_elem):
        dofs = np.array([elem_index, elem_index + 1], dtype=int)
        K[np.ix_(dofs, dofs)] += Ke
        M[np.ix_(dofs, dofs)] += Me
    return K, M


def _clamped_free_bending_shape(
    x: np.ndarray,
    length: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    alpha = CLAMPED_FREE_BENDING_ALPHA
    k = alpha / length
    z = k * x
    sigma = (np.cosh(alpha) + np.cos(alpha)) / (np.sinh(alpha) + np.sin(alpha))
    w = np.cosh(z) - np.cos(z) - sigma * (np.sinh(z) - np.sin(z))
    theta = k * (np.sinh(z) + np.sin(z) - sigma * (np.cosh(z) - np.cos(z)))
    omega_sq = k**4
    return w, theta, omega_sq


def _clamped_free_axial_shape(x: np.ndarray, length: float) -> tuple[np.ndarray, float]:
    alpha = CLAMPED_FREE_AXIAL_ALPHA
    k = alpha / length
    u = np.sin(k * x)
    omega_sq = (float(fem.EA_nd) / float(fem.rhoA_nd)) * k**2
    return u, omega_sq


def _free_reduction(K: np.ndarray, M: np.ndarray, fixed_dofs: list[int]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    fixed = set(fixed_dofs)
    free = np.array([idx for idx in range(K.shape[0]) if idx not in fixed], dtype=int)
    return K[np.ix_(free, free)], M[np.ix_(free, free)], free


def _bending_row(n_elem: int) -> dict[str, object]:
    x = np.linspace(0.0, L_TOTAL, n_elem + 1)
    w, theta, exact_omega_sq = _clamped_free_bending_shape(x, L_TOTAL)
    q = np.empty(2 * (n_elem + 1), dtype=float)
    q[0::2] = w
    q[1::2] = theta

    K, M = _assemble_bending_matrices(n_elem, L_TOTAL)
    Kf, Mf, free = _free_reduction(K, M, fixed_dofs=[0, 1])
    qf = q[free]

    relative_residual, residual_inf, rayleigh_omega_sq = _residual_metrics(
        Kf, Mf, qf, exact_omega_sq
    )
    fem_omega_sq, fem_rel_diff, exact_fem_mac = _first_fem_mode_metrics(
        Kf, Mf, qf, exact_omega_sq
    )
    notes = (
        "bending-only subblock from elem_K/elem_M; "
        f"fem_omega_sq={fem_omega_sq:.16e}; "
        f"fem_rel_diff={fem_rel_diff:.16e}; "
        f"exact_fem_mac={exact_fem_mac:.16e}"
    )
    return {
        "test_name": "single_beam_exact_bending",
        "n_elem": n_elem,
        "boundary_condition": "clamped_free",
        "exact_alpha": f"{CLAMPED_FREE_BENDING_ALPHA:.16e}",
        "exact_omega_sq": f"{exact_omega_sq:.16e}",
        "rayleigh_omega_sq": f"{rayleigh_omega_sq:.16e}",
        "relative_rayleigh_diff": f"{_relative_diff(rayleigh_omega_sq, exact_omega_sq):.16e}",
        "relative_residual": f"{relative_residual:.16e}",
        "residual_inf": f"{residual_inf:.16e}",
        "notes": notes,
    }


def _axial_row(n_elem: int) -> dict[str, object]:
    x = np.linspace(0.0, L_TOTAL, n_elem + 1)
    u, exact_omega_sq = _clamped_free_axial_shape(x, L_TOTAL)

    K, M = _assemble_axial_matrices(n_elem, L_TOTAL)
    Kf, Mf, free = _free_reduction(K, M, fixed_dofs=[0])
    qf = u[free]

    relative_residual, residual_inf, rayleigh_omega_sq = _residual_metrics(
        Kf, Mf, qf, exact_omega_sq
    )
    fem_omega_sq, fem_rel_diff, exact_fem_mac = _first_fem_mode_metrics(
        Kf, Mf, qf, exact_omega_sq
    )
    notes = (
        "axial-only subblock from elem_K/elem_M; "
        f"EA_nd={float(fem.EA_nd):.16e}; "
        f"fem_omega_sq={fem_omega_sq:.16e}; "
        f"fem_rel_diff={fem_rel_diff:.16e}; "
        f"exact_fem_mac={exact_fem_mac:.16e}"
    )
    return {
        "test_name": "single_bar_exact_axial",
        "n_elem": n_elem,
        "boundary_condition": "clamped_free",
        "exact_alpha": f"{CLAMPED_FREE_AXIAL_ALPHA:.16e}",
        "exact_omega_sq": f"{exact_omega_sq:.16e}",
        "rayleigh_omega_sq": f"{rayleigh_omega_sq:.16e}",
        "relative_rayleigh_diff": f"{_relative_diff(rayleigh_omega_sq, exact_omega_sq):.16e}",
        "relative_residual": f"{relative_residual:.16e}",
        "residual_inf": f"{residual_inf:.16e}",
        "notes": notes,
    }


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def _print_table(title: str, rows: list[dict[str, object]]) -> None:
    print(title)
    print("n_elem  relative_residual       relative_rayleigh_diff")
    for row in rows:
        print(
            f"{int(row['n_elem']):6d}  "
            f"{float(row['relative_residual']):.8e}  "
            f"{float(row['relative_rayleigh_diff']):.8e}"
        )
    print()


def main() -> int:
    bending_rows = [_bending_row(n_elem) for n_elem in N_ELEM_VALUES]
    axial_rows = [_axial_row(n_elem) for n_elem in N_ELEM_VALUES]

    bending_path = RESULTS_DIR / "single_beam_exact_shape_residual_bending.csv"
    axial_path = RESULTS_DIR / "single_beam_exact_shape_residual_axial.csv"
    _write_csv(bending_path, bending_rows)
    _write_csv(axial_path, axial_rows)

    _print_table("Single Euler-Bernoulli cantilever bending residual:", bending_rows)
    _print_table("Single axial cantilever bar residual:", axial_rows)
    print(f"Wrote {bending_path}")
    print(f"Wrote {axial_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
