"""Audit right-arm FEM local/global transform variants.

This script intentionally does not modify the baseline FEM implementation.  It
reassembles diagnostic K/M matrices with several right-arm transform formulas
and checks whether the current determinant-null article embedding is compatible
with any of them.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy.linalg import eigh

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import assemble_clamped_coupled_matrix  # noqa: E402
from my_project.fem import python_fem as fem  # noqa: E402
from scripts.analysis.check_analytic_shape_in_fem_residual import (  # noqa: E402
    analytic_local_fields,
    build_analytic_global_q,
    build_params_for_case,
    euclidean_mac,
    expand_full_vector,
    free_dof_indices,
    mass_weighted_mac,
    residual_metrics,
    resolve_analytic_branch,
)
from scripts.compare_beta0_analytic_vs_fem import fem_parameter_override  # noqa: E402
from scripts.lib.analytic_coupled_rods_shapes import NEAR_ZERO_NORM, analytic_null_vector  # noqa: E402


DEFAULT_BRANCH_ID = "bending_desc_05"
DEFAULT_BETA = 15.0
DEFAULT_EPSILON = 0.0025
DEFAULT_MU_VALUES = (0.0, 0.2)
DEFAULT_L_TOTAL = 2.0
DEFAULT_NUM_ANALYTIC_ROOTS = 20
DEFAULT_N_FEM_MODES = 30
DEFAULT_AUDIT_CSV = REPO_ROOT / "results" / "fem_right_transform_variant_audit_beta15_bending_desc_05.csv"
DEFAULT_SANITY_CSV = REPO_ROOT / "results" / "fem_right_transform_variant_sanity.csv"
ARTICLE_RIGHT_THETA_SIGN = 1.0

AUDIT_FIELDNAMES = [
    "mu",
    "variant_name",
    "fem_lambda",
    "analytic_lambda",
    "rel_lambda_diff",
    "analytic_in_variant_fem_residual",
    "rayleigh_lambda",
    "global_dof_mac",
    "mass_weighted_mac",
    "q_joint_mac",
    "notes",
]

SANITY_FIELDNAMES = [
    "variant_name",
    "assembly_formula",
    "local_to_global_inferred",
    "pure_local_axial_global_ux",
    "pure_local_axial_global_uy",
    "pure_local_transverse_global_ux",
    "pure_local_transverse_global_uy",
    "matches_article_kinematics",
    "notes",
]


@dataclass(frozen=True)
class TransformVariant:
    name: str
    assembly_formula: str
    assembly_t_kind: str
    energy_formula: str
    notes: str


VARIANTS = (
    TransformVariant(
        name="A_baseline",
        assembly_formula="T=R(+beta); K_g=T K_l T.T; M_g=T M_l T.T",
        assembly_t_kind="R_plus",
        energy_formula="TKTt",
        notes="same right-arm transform formula as the corrected python_fem.py baseline",
    ),
    TransformVariant(
        name="B_local_to_global_energy",
        assembly_formula="T=R(+beta); K_g=T K_l T.T; M_g=T M_l T.T",
        assembly_t_kind="R_plus",
        energy_formula="TKTt",
        notes="energy form for q_global=T q_local with T=R(+beta)",
    ),
    TransformVariant(
        name="C_transpose_T_baseline_formula",
        assembly_formula="T=R(+beta).T; K_g=T.T K_l T; M_g=T.T M_l T",
        assembly_t_kind="R_plus_transpose",
        energy_formula="TtKT",
        notes="baseline formula with transposed right-arm T",
    ),
    TransformVariant(
        name="D_beta_negative_baseline",
        assembly_formula="T=R(-beta); K_g=T.T K_l T; M_g=T.T M_l T",
        assembly_t_kind="R_minus",
        energy_formula="TtKT",
        notes="baseline formula with negative beta in right-arm T",
    ),
)


def parse_mu_values(raw: str) -> tuple[float, ...]:
    values = tuple(float(part.strip()) for part in raw.split(",") if part.strip())
    if not values:
        raise argparse.ArgumentTypeError("expected at least one comma-separated mu value")
    return values


def resolve_output_path(path: str | Path) -> Path:
    value = Path(path)
    return value if value.is_absolute() else REPO_ROOT / value


def rotation3(beta_rad: float) -> np.ndarray:
    c, s = np.cos(beta_rad), np.sin(beta_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)


def rotation6(beta_rad: float) -> np.ndarray:
    r3 = rotation3(beta_rad)
    t = np.zeros((6, 6), dtype=float)
    t[:3, :3] = r3
    t[3:, 3:] = r3
    return t


def variant_t_matrix(variant: TransformVariant, beta_deg: float) -> np.ndarray:
    beta_rad = float(np.deg2rad(beta_deg))
    if variant.assembly_t_kind == "R_plus":
        return rotation6(beta_rad)
    if variant.assembly_t_kind == "R_plus_transpose":
        return rotation6(beta_rad).T
    if variant.assembly_t_kind == "R_minus":
        return rotation6(-beta_rad)
    raise ValueError(f"unknown variant T kind: {variant.assembly_t_kind}")


def transform_element_matrices(
    variant: TransformVariant,
    k_local: np.ndarray,
    m_local: np.ndarray,
    beta_deg: float,
) -> tuple[np.ndarray, np.ndarray]:
    t = variant_t_matrix(variant, beta_deg)
    if variant.energy_formula == "TtKT":
        return t.T @ k_local @ t, t.T @ m_local @ t
    if variant.energy_formula == "TKTt":
        return t @ k_local @ t.T, t @ m_local @ t.T
    raise ValueError(f"unknown variant energy formula: {variant.energy_formula}")


def variant_q_local_from_global_matrix(variant: TransformVariant, beta_deg: float) -> np.ndarray:
    t = variant_t_matrix(variant, beta_deg)
    if variant.energy_formula == "TtKT":
        return t
    if variant.energy_formula == "TKTt":
        return t.T
    raise ValueError(f"unknown variant energy formula: {variant.energy_formula}")


def variant_global_from_local_matrix(variant: TransformVariant, beta_deg: float) -> np.ndarray:
    q_local_from_global = variant_q_local_from_global_matrix(variant, beta_deg)
    return q_local_from_global.T


def assemble_variant_matrices(
    *,
    params,
    mu: float,
    beta_deg: float,
    variant: TransformVariant,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with fem_parameter_override(params):
        n = fem.N_ELEM
        ndof = 3 * (2 * n + 1)
        length_left = max(float(fem.ell) * (1.0 - float(mu)), 1.0e-3)
        length_right = max(float(fem.ell) * (1.0 + float(mu)), 1.0e-3)
        le_left = length_left / n
        le_right = length_right / n
        k_left = fem.elem_K(le_left)
        m_left = fem.elem_M(le_left)
        k_right, m_right = transform_element_matrices(
            variant,
            fem.elem_K(le_right),
            fem.elem_M(le_right),
            beta_deg,
        )

        k_full = np.zeros((ndof, ndof), dtype=float)
        m_full = np.zeros((ndof, ndof), dtype=float)

        def assemble(dof_map: list[int], ke: np.ndarray, me: np.ndarray) -> None:
            for row, dof_row in enumerate(dof_map):
                for col, dof_col in enumerate(dof_map):
                    k_full[dof_row, dof_col] += ke[row, col]
                    m_full[dof_row, dof_col] += me[row, col]

        for elem in range(n):
            dofs = [
                3 * elem,
                3 * elem + 1,
                3 * elem + 2,
                3 * (elem + 1),
                3 * (elem + 1) + 1,
                3 * (elem + 1) + 2,
            ]
            assemble(dofs, k_left, m_left)

        for elem in range(n):
            node0 = n + elem
            node1 = n + elem + 1
            dofs = [
                3 * node0,
                3 * node0 + 1,
                3 * node0 + 2,
                3 * node1,
                3 * node1 + 1,
                3 * node1 + 2,
            ]
            assemble(dofs, k_right, m_right)

    free = free_dof_indices()
    return (
        k_full[np.ix_(free, free)],
        m_full[np.ix_(free, free)],
        k_full,
        m_full,
        free,
    )


def solve_variant_modes(
    k_free: np.ndarray,
    m_free: np.ndarray,
    *,
    n_modes: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    eigvals, eigvecs = eigh(k_free, m_free, subset_by_index=[0, int(n_modes) - 1])
    eigvals = np.maximum(np.asarray(eigvals, dtype=float), 0.0)
    omega = np.sqrt(eigvals)
    lambdas = np.sqrt(omega)
    return lambdas, omega, np.asarray(eigvecs, dtype=float)


def rayleigh_lambda(k_free: np.ndarray, m_free: np.ndarray, q_free: np.ndarray) -> float:
    q = np.asarray(q_free, dtype=float)
    denominator = float(q @ (m_free @ q))
    numerator = float(q @ (k_free @ q))
    if denominator <= NEAR_ZERO_NORM or numerator < 0.0:
        return float("nan")
    return float((numerator / denominator) ** 0.25)


def relative_difference(value: float, reference: float) -> float:
    return abs(float(value) - float(reference)) / abs(float(reference)) if abs(float(reference)) > NEAR_ZERO_NORM else np.nan


def analytic_article_shape(
    *,
    branch_id: str,
    beta_deg: float,
    mu: float,
    epsilon: float,
    num_analytic_roots: int,
) -> tuple[float, np.ndarray, np.ndarray, dict[str, float], object]:
    tracking = resolve_analytic_branch(
        branch_id=branch_id,
        beta_deg=float(beta_deg),
        mu=float(mu),
        epsilon=float(epsilon),
        num_analytic_roots=int(num_analytic_roots),
    )
    point = tracking["point"]
    analytic_lambda = float(point.Lambda)
    matrix = assemble_clamped_coupled_matrix(
        analytic_lambda,
        float(np.deg2rad(beta_deg)),
        float(mu),
        float(epsilon),
    )
    coeff, _, _ = analytic_null_vector(matrix)
    fields = analytic_local_fields(
        analytic_lambda,
        mu=float(mu),
        epsilon=float(epsilon),
        coeff=coeff,
        s_norm=np.linspace(0.0, 1.0, fem.N_ELEM + 1),
        right_theta_sign=ARTICLE_RIGHT_THETA_SIGN,
    )
    q_article, embedding = build_analytic_global_q(fields, beta_deg=float(beta_deg))
    return analytic_lambda, coeff, q_article, embedding, point


def joint_vector(q_full: np.ndarray) -> np.ndarray:
    n = fem.N_ELEM
    return np.asarray(q_full, dtype=float)[3 * n : 3 * n + 3]


def audit_rows(
    *,
    branch_id: str,
    beta_deg: float,
    epsilon: float,
    mu_values: Sequence[float],
    num_analytic_roots: int,
    n_fem_modes: int,
) -> list[dict[str, float | str]]:
    params = build_params_for_case(epsilon=float(epsilon), l_total=DEFAULT_L_TOTAL)
    rows: list[dict[str, float | str]] = []
    for mu in mu_values:
        analytic_lambda, _, q_article, embedding, point = analytic_article_shape(
            branch_id=branch_id,
            beta_deg=float(beta_deg),
            mu=float(mu),
            epsilon=float(epsilon),
            num_analytic_roots=int(num_analytic_roots),
        )
        analytic_free_cache: dict[tuple[int, ...], np.ndarray] = {}

        for variant in VARIANTS:
            k_free, m_free, _, _, free = assemble_variant_matrices(
                params=params,
                mu=float(mu),
                beta_deg=float(beta_deg),
                variant=variant,
            )
            lambdas, _, vecs = solve_variant_modes(k_free, m_free, n_modes=int(n_fem_modes))
            selected_idx = int(np.argmin(np.abs(lambdas - analytic_lambda)))
            fem_lambda = float(lambdas[selected_idx])
            fem_full = expand_full_vector(vecs[:, selected_idx])
            free_key = tuple(int(index) for index in free)
            q_free = analytic_free_cache.setdefault(free_key, np.asarray(q_article, dtype=float)[free])
            fem_free = fem_full[free]
            residual = residual_metrics(k_free, m_free, q_article, free, omega_sq=analytic_lambda**4)
            rayleigh = rayleigh_lambda(k_free, m_free, q_free)
            q_joint_analytic = joint_vector(q_article)
            q_joint_fem = joint_vector(fem_full)
            notes = (
                f"selected_fem_mode_index={selected_idx + 1}; "
                f"analytic_current_sorted_index={int(point.current_sorted_index)}; "
                f"article_embedding=R(+beta) local_to_global, right_theta_sign=+1; "
                f"joint_translation_mismatch={float(embedding['joint_translation_embedding_mismatch']):.3e}; "
                f"joint_theta_mismatch={float(embedding['joint_theta_embedding_mismatch']):.3e}; "
                f"{variant.notes}"
            )
            rows.append(
                {
                    "mu": float(mu),
                    "variant_name": variant.name,
                    "fem_lambda": fem_lambda,
                    "analytic_lambda": analytic_lambda,
                    "rel_lambda_diff": relative_difference(fem_lambda, analytic_lambda),
                    "analytic_in_variant_fem_residual": float(residual["relative_residual"]),
                    "rayleigh_lambda": rayleigh,
                    "global_dof_mac": euclidean_mac(q_free, fem_free),
                    "mass_weighted_mac": mass_weighted_mac(q_free, fem_free, m_free),
                    "q_joint_mac": euclidean_mac(q_joint_analytic, q_joint_fem),
                    "notes": notes,
                }
            )
    return rows


def sanity_rows(beta_deg: float) -> list[dict[str, float | str]]:
    article_axial = np.array([np.cos(np.deg2rad(beta_deg)), np.sin(np.deg2rad(beta_deg))], dtype=float)
    article_transverse = np.array([-np.sin(np.deg2rad(beta_deg)), np.cos(np.deg2rad(beta_deg))], dtype=float)
    rows: list[dict[str, float | str]] = []
    for variant in VARIANTS:
        global_from_local = variant_global_from_local_matrix(variant, beta_deg)
        axial = global_from_local[:2, :2] @ np.array([1.0, 0.0], dtype=float)
        transverse = global_from_local[:2, :2] @ np.array([0.0, 1.0], dtype=float)
        matches_article = bool(
            np.allclose(axial, article_axial, rtol=1.0e-12, atol=1.0e-12)
            and np.allclose(transverse, article_transverse, rtol=1.0e-12, atol=1.0e-12)
        )
        rows.append(
            {
                "variant_name": variant.name,
                "assembly_formula": variant.assembly_formula,
                "local_to_global_inferred": (
                    "global_from_local=(q_local_from_global).T inferred from the quadratic energy form"
                ),
                "pure_local_axial_global_ux": float(axial[0]),
                "pure_local_axial_global_uy": float(axial[1]),
                "pure_local_transverse_global_ux": float(transverse[0]),
                "pure_local_transverse_global_uy": float(transverse[1]),
                "matches_article_kinematics": "yes" if matches_article else "no",
                "notes": variant.notes,
            }
        )
    return rows


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Diagnostic-only audit of right-arm FEM transform assembly variants."
    )
    parser.add_argument("--branch-id", default=DEFAULT_BRANCH_ID)
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--mu-values", type=parse_mu_values, default=DEFAULT_MU_VALUES)
    parser.add_argument("--num-analytic-roots", type=int, default=DEFAULT_NUM_ANALYTIC_ROOTS)
    parser.add_argument("--n-fem-modes", type=int, default=DEFAULT_N_FEM_MODES)
    parser.add_argument("--audit-csv", type=Path, default=DEFAULT_AUDIT_CSV)
    parser.add_argument("--sanity-csv", type=Path, default=DEFAULT_SANITY_CSV)
    return parser


def print_summary(rows: Sequence[dict[str, float | str]], sanity: Sequence[dict[str, float | str]]) -> None:
    print("FEM right-arm transform variant audit")
    print("Sanity local-to-global directions:")
    for row in sanity:
        print(
            "  "
            f"{row['variant_name']}: axial=({float(row['pure_local_axial_global_ux']):+.6f}, "
            f"{float(row['pure_local_axial_global_uy']):+.6f}), transverse=("
            f"{float(row['pure_local_transverse_global_ux']):+.6f}, "
            f"{float(row['pure_local_transverse_global_uy']):+.6f}), "
            f"article={row['matches_article_kinematics']}"
        )
    print("Best residual per mu:")
    for mu in sorted({float(row["mu"]) for row in rows}):
        case_rows = [row for row in rows if abs(float(row["mu"]) - mu) <= 1.0e-12]
        best = min(case_rows, key=lambda row: float(row["analytic_in_variant_fem_residual"]))
        baseline = next(row for row in case_rows if row["variant_name"] == "A_baseline")
        print(
            "  "
            f"mu={mu:g}: best={best['variant_name']} residual="
            f"{float(best['analytic_in_variant_fem_residual']):.6e}; baseline="
            f"{float(baseline['analytic_in_variant_fem_residual']):.6e}"
        )
    print("Rows:")
    for row in rows:
        print(
            "  "
            f"mu={float(row['mu']):g}, {row['variant_name']}: "
            f"fem_lambda={float(row['fem_lambda']):.9g}, "
            f"analytic_lambda={float(row['analytic_lambda']):.9g}, "
            f"rel_lambda_diff={float(row['rel_lambda_diff']):.3e}, "
            f"residual={float(row['analytic_in_variant_fem_residual']):.3e}, "
            f"global_MAC={float(row['global_dof_mac']):.3e}, "
            f"q_joint_MAC={float(row['q_joint_mac']):.3e}"
        )


def main(argv: Sequence[str] | None = None) -> None:
    args = build_arg_parser().parse_args(argv)
    rows = audit_rows(
        branch_id=str(args.branch_id),
        beta_deg=float(args.beta),
        epsilon=float(args.epsilon),
        mu_values=tuple(float(value) for value in args.mu_values),
        num_analytic_roots=int(args.num_analytic_roots),
        n_fem_modes=int(args.n_fem_modes),
    )
    sanity = sanity_rows(float(args.beta))
    audit_csv = resolve_output_path(args.audit_csv)
    sanity_csv = resolve_output_path(args.sanity_csv)
    write_csv(audit_csv, AUDIT_FIELDNAMES, rows)
    write_csv(sanity_csv, SANITY_FIELDNAMES, sanity)
    print_summary(rows, sanity)
    print(f"audit CSV: {audit_csv.relative_to(REPO_ROOT)}")
    print(f"sanity CSV: {sanity_csv.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
