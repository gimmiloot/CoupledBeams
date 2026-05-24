from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.linalg import eigh


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    thickness_mismatch_factors,
    thickness_to_length_ratios,
)
from my_project.fem import python_fem as baseline_fem  # noqa: E402
from scripts.lib.thickness_mismatch_diagnostic_helpers import (  # noqa: E402
    mu_grid,
    plot_with_validity_split,
    rods_label,
    roots_by_mu_eta,
    thickness_ratio_summary,
    track_descendants_from_mu0,
    tracking_warning_rows,
    valid_mask,
)


# =========================
# User-editable parameters
# =========================
BETA_DEG = 15.0
EPSILON = 0.0025
ETA = 0.5

MU_MIN = 0.0
MU_MAX = 0.9
MU_STEP = 0.005
MU_FEM_VALUES = tuple(np.round(np.arange(0.0, 0.9 + 0.5 * 0.05, 0.05), 10))

N_MODES_COMPARE = 6
N_SORTED_SCAN = 12
N_DESCENDANTS_TRACK = 8
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401

FEM_ELEMENTS_PER_ROD = 80
FEM_CONVERGENCE_ELEMENTS = (40, 80)

MAC_WARNING_THRESHOLD = 0.9
MAC_MARGIN_WARNING_THRESHOLD = 0.05
MAX_SORTED_POSITION_JUMP = 1
THICKNESS_RATIO_LIMIT = 0.1

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_PNG = OUTPUT_DIR / "thickness_mismatch_fem_comparison_beta15_eps0p0025_eta_p0p5.png"
OUTPUT_REPORT = OUTPUT_DIR / "thickness_mismatch_fem_comparison_beta15_eps0p0025_eta_p0p5_report.md"

MU_VALUES = mu_grid(MU_MIN, MU_MAX, MU_STEP)
NEAR_ZERO_NORM = 1e-14


@dataclass(frozen=True)
class Geometry:
    mu: float
    tau1: float
    tau2: float
    length1_nd: float
    length2_nd: float
    area1_factor: float
    area2_factor: float
    inertia1_factor: float
    inertia2_factor: float
    thickness_ratio_1: float
    thickness_ratio_2: float

    @property
    def validity_status(self) -> str:
        if (
            self.thickness_ratio_1 <= THICKNESS_RATIO_LIMIT + 1e-12
            and self.thickness_ratio_2 <= THICKNESS_RATIO_LIMIT + 1e-12
        ):
            return "valid"
        return "criterion_violated"


@dataclass(frozen=True)
class FemSolveResult:
    elements_per_rod: int
    omega_nd: np.ndarray
    lambda_values: np.ndarray


def geometry_for(mu: float) -> Geometry:
    factors = thickness_mismatch_factors(float(mu), ETA)
    ratio1, ratio2 = thickness_to_length_ratios(EPSILON, float(mu), ETA)
    return Geometry(
        mu=float(mu),
        tau1=float(factors.tau1),
        tau2=float(factors.tau2),
        length1_nd=1.0 - float(mu),
        length2_nd=1.0 + float(mu),
        area1_factor=float(factors.tau1**2),
        area2_factor=float(factors.tau2**2),
        inertia1_factor=float(factors.tau1**4),
        inertia2_factor=float(factors.tau2**4),
        thickness_ratio_1=float(ratio1),
        thickness_ratio_2=float(ratio2),
    )


def elem_k_param(length: float, *, ea_nd: float, ei_nd: float) -> np.ndarray:
    """Euler-Bernoulli frame stiffness in the existing FEM convention."""
    le = float(length)
    matrix = np.zeros((6, 6), dtype=float)

    axial = float(ea_nd) / le
    matrix[0, 0] += axial
    matrix[0, 3] -= axial
    matrix[3, 0] -= axial
    matrix[3, 3] += axial

    bending = float(ei_nd) / le**3
    block = bending * np.array(
        [
            [12.0, 6.0 * le, -12.0, 6.0 * le],
            [6.0 * le, 4.0 * le**2, -6.0 * le, 2.0 * le**2],
            [-12.0, -6.0 * le, 12.0, -6.0 * le],
            [6.0 * le, 2.0 * le**2, -6.0 * le, 4.0 * le**2],
        ],
        dtype=float,
    )
    for i, row in enumerate([1, 2, 4, 5]):
        for j, col in enumerate([1, 2, 4, 5]):
            matrix[row, col] += block[i, j]
    return matrix


def elem_m_param(length: float, *, rhoa_nd: float) -> np.ndarray:
    """Consistent axial + bending frame mass in the existing FEM convention."""
    le = float(length)
    matrix = np.zeros((6, 6), dtype=float)

    axial = float(rhoa_nd) * le / 6.0
    matrix[0, 0] += 2.0 * axial
    matrix[0, 3] += axial
    matrix[3, 0] += axial
    matrix[3, 3] += 2.0 * axial

    bending = float(rhoa_nd) * le / 420.0
    block = bending * np.array(
        [
            [156.0, 22.0 * le, 54.0, -13.0 * le],
            [22.0 * le, 4.0 * le**2, 13.0 * le, -3.0 * le**2],
            [54.0, 13.0 * le, 156.0, -22.0 * le],
            [-13.0 * le, -3.0 * le**2, -22.0 * le, 4.0 * le**2],
        ],
        dtype=float,
    )
    for i, row in enumerate([1, 2, 4, 5]):
        for j, col in enumerate([1, 2, 4, 5]):
            matrix[row, col] += block[i, j]
    return matrix


def assemble_mismatch_fem(mu: float, elements_per_rod: int) -> tuple[np.ndarray, np.ndarray]:
    geom = geometry_for(float(mu))
    n1 = n2 = int(elements_per_rod)
    le1 = geom.length1_nd / n1
    le2 = geom.length2_nd / n2
    ndof = 3 * (n1 + n2 + 1)
    stiffness = np.zeros((ndof, ndof), dtype=float)
    mass = np.zeros((ndof, ndof), dtype=float)

    def assemble(dofs: Sequence[int], local_k: np.ndarray, local_m: np.ndarray) -> None:
        for row_idx, row_dof in enumerate(dofs):
            for col_idx, col_dof in enumerate(dofs):
                stiffness[row_dof, col_dof] += local_k[row_idx, col_idx]
                mass[row_dof, col_dof] += local_m[row_idx, col_idx]

    # The nondimensional axial scale matches the existing diagnostic FEM:
    # EA/(EJ0/l^2) = tau_i^2 / epsilon^2, EI/(EJ0) = tau_i^4.
    k1 = elem_k_param(le1, ea_nd=geom.area1_factor / EPSILON**2, ei_nd=geom.inertia1_factor)
    m1 = elem_m_param(le1, rhoa_nd=geom.area1_factor)
    k2_local = elem_k_param(le2, ea_nd=geom.area2_factor / EPSILON**2, ei_nd=geom.inertia2_factor)
    m2_local = elem_m_param(le2, rhoa_nd=geom.area2_factor)
    transform = baseline_fem.rotation_matrix_6x6(np.deg2rad(BETA_DEG))
    k2 = transform @ k2_local @ transform.T
    m2 = transform @ m2_local @ transform.T

    for elem in range(n1):
        dofs = [3 * elem, 3 * elem + 1, 3 * elem + 2, 3 * (elem + 1), 3 * (elem + 1) + 1, 3 * (elem + 1) + 2]
        assemble(dofs, k1, m1)

    for elem in range(n2):
        base = n1 + elem
        dofs = [
            3 * base,
            3 * base + 1,
            3 * base + 2,
            3 * (base + 1),
            3 * (base + 1) + 1,
            3 * (base + 1) + 2,
        ]
        assemble(dofs, k2, m2)

    fixed = {0, 1, 2, 3 * (n1 + n2), 3 * (n1 + n2) + 1, 3 * (n1 + n2) + 2}
    free = np.array(sorted(set(range(ndof)) - fixed), dtype=int)
    return stiffness[np.ix_(free, free)], mass[np.ix_(free, free)]


def solve_mismatch_fem(mu: float, elements_per_rod: int, n_modes: int) -> FemSolveResult:
    stiffness, mass = assemble_mismatch_fem(float(mu), int(elements_per_rod))
    eigenvalues, _eigenvectors = eigh(stiffness, mass, subset_by_index=[0, int(n_modes) - 1])
    omega_nd = np.sqrt(np.maximum(eigenvalues, 0.0))
    lambda_values = np.sqrt(np.maximum(omega_nd, 0.0))
    return FemSolveResult(
        elements_per_rod=int(elements_per_rod),
        omega_nd=omega_nd,
        lambda_values=lambda_values,
    )


def analytic_tracking() -> tuple[np.ndarray, list[dict[str, float | int | str]]]:
    roots = roots_by_mu_eta(
        beta_rad=float(np.deg2rad(BETA_DEG)),
        epsilon=EPSILON,
        eta=ETA,
        mu_values=MU_VALUES,
        n_roots=N_SORTED_SCAN,
        root_lmax0=ROOT_LMAX0,
        root_scan_step=ROOT_SCAN_STEP,
    )
    tracking = track_descendants_from_mu0(
        beta_rad=float(np.deg2rad(BETA_DEG)),
        epsilon=EPSILON,
        eta=ETA,
        mu_values=MU_VALUES,
        roots_by_mu=roots,
        num_descendants=N_DESCENDANTS_TRACK,
        num_shape_samples=NUM_SHAPE_SAMPLES,
        mac_warning_threshold=MAC_WARNING_THRESHOLD,
        mac_margin_warning_threshold=MAC_MARGIN_WARNING_THRESHOLD,
        max_sorted_position_jump=MAX_SORTED_POSITION_JUMP,
    )
    return tracking.tracked, tracking.rows


def mu_index(mu: float) -> int:
    matches = np.where(np.isclose(MU_VALUES, float(mu), rtol=0.0, atol=1e-12))[0]
    if len(matches) != 1:
        raise RuntimeError(f"mu={float(mu):g} is not on the analytic grid.")
    return int(matches[0])


def relative_difference(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / abs(float(right)) if abs(float(right)) > NEAR_ZERO_NORM else np.nan


def compute_fem_by_mesh() -> dict[int, dict[float, FemSolveResult]]:
    fem: dict[int, dict[float, FemSolveResult]] = {}
    for mesh in FEM_CONVERGENCE_ELEMENTS:
        fem[int(mesh)] = {}
        for mu in MU_FEM_VALUES:
            fem[int(mesh)][float(mu)] = solve_mismatch_fem(float(mu), int(mesh), N_MODES_COMPARE)
    return fem


def comparison_rows(
    *,
    tracked: np.ndarray,
    fem_by_mesh: dict[int, dict[float, FemSolveResult]],
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    fine = fem_by_mesh[int(FEM_ELEMENTS_PER_ROD)]
    coarse = fem_by_mesh[int(FEM_CONVERGENCE_ELEMENTS[0])]
    for mu in MU_FEM_VALUES:
        col = mu_index(float(mu))
        fem_lambda = fine[float(mu)].lambda_values
        coarse_lambda = coarse[float(mu)].lambda_values
        for mode in range(1, N_MODES_COMPARE + 1):
            analytic = float(tracked[mode - 1, col])
            fem_same = float(fem_lambda[mode - 1])
            nearest_idx = int(np.argmin(np.abs(fem_lambda - analytic)))
            rows.append(
                {
                    "mu": float(mu),
                    "descendant_branch": int(mode),
                    "Lambda_analytic_descendant": analytic,
                    "fem_sorted_mode_same_index": int(mode),
                    "Lambda_fem_same_index": fem_same,
                    "abs_diff": abs(analytic - fem_same),
                    "rel_diff": relative_difference(analytic, fem_same),
                    "nearest_fem_sorted_mode": int(nearest_idx) + 1,
                    "nearest_fem_abs_diff": abs(analytic - float(fem_lambda[nearest_idx])),
                    "nearest_fem_rel_diff": relative_difference(analytic, float(fem_lambda[nearest_idx])),
                    "Lambda_fem_coarse": float(coarse_lambda[mode - 1]),
                    "mesh_abs_change": abs(float(fem_lambda[mode - 1]) - float(coarse_lambda[mode - 1])),
                    "mesh_rel_change": relative_difference(float(fem_lambda[mode - 1]), float(coarse_lambda[mode - 1])),
                    "validity_status": geometry_for(float(mu)).validity_status,
                }
            )
    return rows


def plot_comparison(*, tracked: np.ndarray, fem_by_mesh: dict[int, dict[float, FemSolveResult]]) -> None:
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, ax = plt.subplots(figsize=(9.8, 5.7))
    valid = valid_mask(epsilon=EPSILON, eta=ETA, mu_values=MU_VALUES, limit=THICKNESS_RATIO_LIMIT)

    for branch in range(1, N_MODES_COMPARE + 1):
        color = colors[(branch - 1) % len(colors)]
        plot_with_validity_split(
            ax,
            MU_VALUES,
            tracked[branch - 1],
            valid,
            color=color,
            linewidth=1.9,
            label=f"desc {branch}",
            zorder=3,
        )

    fine = fem_by_mesh[int(FEM_ELEMENTS_PER_ROD)]
    for mode in range(1, N_MODES_COMPARE + 1):
        color = colors[(mode - 1) % len(colors)]
        marker = "s" if mode in (5, 6) else "o"
        size = 42 if mode in (5, 6) else 28
        ax.scatter(
            MU_FEM_VALUES,
            [fine[float(mu)].lambda_values[mode - 1] for mu in MU_FEM_VALUES],
            color=color,
            marker=marker,
            s=size,
            zorder=5,
        )

    handles = [
        Line2D([0], [0], color="black", lw=1.9, label="analytic descendants"),
        Line2D([0], [0], color="black", marker="o", lw=0.0, label="FEM sorted modes"),
        Line2D([0], [0], color="0.25", lw=1.7, ls="-", label=r"solid: $2r_i/l_i \leq 0.1$"),
        Line2D([0], [0], color="0.25", lw=1.7, ls="--", label="dashed: criterion violated"),
    ]
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\Lambda$")
    ax.set_title(rf"Analytic descendant branches vs FEM, $\eta={ETA:g}$, $\beta={BETA_DEG:g}^\circ$")
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(handles=handles, frameon=False, fontsize=9, loc="best")
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=240, bbox_inches="tight")
    plt.close(fig)


def max_row(rows: Sequence[dict[str, float | int | str]], key: str) -> dict[str, float | int | str]:
    return max(rows, key=lambda row: float(row[key]) if np.isfinite(float(row[key])) else -np.inf)


def write_report(
    *,
    rows: list[dict[str, float | int | str]],
    tracking_rows: list[dict[str, float | int | str]],
) -> dict[str, object]:
    OUTPUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    abs_row = max_row(rows, "abs_diff")
    rel_row = max_row(rows, "rel_diff")
    mesh_abs_row = max_row(rows, "mesh_abs_change")
    mesh_rel_row = max_row(rows, "mesh_rel_change")
    nearest_mismatch = [row for row in rows if int(row["nearest_fem_sorted_mode"]) != int(row["descendant_branch"])]
    tracking_warnings = tracking_warning_rows(tracking_rows)
    summary = thickness_ratio_summary(
        epsilon=EPSILON,
        eta=ETA,
        mu_values=MU_FEM_VALUES,
        limit=THICKNESS_RATIO_LIMIT,
    )

    lines = [
        "# Thickness-Mismatch FEM Comparison: eta=0.5",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g}",
        f"- analytic mu grid: {MU_MIN:g}..{MU_MAX:g}, step {MU_STEP:g}",
        f"- FEM mu values: {', '.join(f'{float(mu):g}' for mu in MU_FEM_VALUES)}",
        f"- compared modes/descendants: first {N_MODES_COMPARE}",
        f"- FEM elements per rod: primary {FEM_ELEMENTS_PER_ROD}; convergence {FEM_CONVERGENCE_ELEMENTS}",
        "",
        "## FEM Formulation",
        "",
        "The script uses the existing planar Euler-Bernoulli frame-element",
        "convention: axial bar stiffness, Euler-Bernoulli bending stiffness,",
        "consistent axial and bending mass, clamped external ends, and shared",
        "rigid-joint DOFs. The right arm uses the production transform from",
        "`src/my_project/fem/python_fem.py`: `q_global = T @ q_local`,",
        "`K_global = T @ K_local @ T.T`, and `M_global = T @ M_local @ T.T`.",
        "",
        "The existing FEM module is not modified. This diagnostic only assembles",
        "arm-specific area and inertia factors for the two-radius eta model.",
        "",
        "## Lambda Normalization",
        "",
        "The nondimensional FEM eigenproblem uses the equal-radius reference",
        "section. Its angular-frequency parameter satisfies",
        "`omega_nd = Omega*l^2*sqrt(rho*S0/(E*J0)) = Lambda^2`, so the plotted",
        "FEM value is `Lambda_fem = sqrt(omega_nd)`. Equivalently,",
        "`Lambda_fem = sqrt(Omega)*l*(rho*S0/(E*J0))^(1/4)`.",
        "",
        "## Branch/FEM Matching",
        "",
        "Analytic curves are descendant branches seeded at `mu=0` and tracked by",
        "analytic shape MAC. FEM markers are sorted FEM eigenfrequencies. The",
        "main plotted comparison overlays descendant branch `k` with FEM sorted",
        "mode `k`; nearest FEM mode is recorded as a warning diagnostic, not as",
        "a branch renumbering rule.",
        "",
        "## Errors",
        "",
        f"- max absolute difference: {float(abs_row['abs_diff']):.6e} "
        f"at mu={float(abs_row['mu']):g}, descendant/mode={int(abs_row['descendant_branch'])}",
        f"- max relative difference: {float(rel_row['rel_diff']):.6e} "
        f"at mu={float(rel_row['mu']):g}, descendant/mode={int(rel_row['descendant_branch'])}",
        f"- max mesh absolute change: {float(mesh_abs_row['mesh_abs_change']):.6e} "
        f"at mu={float(mesh_abs_row['mu']):g}, mode={int(mesh_abs_row['descendant_branch'])}",
        f"- max mesh relative change: {float(mesh_rel_row['mesh_rel_change']):.6e} "
        f"at mu={float(mesh_rel_row['mu']):g}, mode={int(mesh_rel_row['descendant_branch'])}",
        "",
        "## Warnings",
        "",
        f"- tracking warning rows on analytic grid: {len(tracking_warnings)}",
        f"- nearest-FEM-mode mismatches on FEM grid: {len(nearest_mismatch)}",
    ]
    if nearest_mismatch:
        for row in nearest_mismatch[:12]:
            lines.append(
                f"  - mu={float(row['mu']):g}, descendant={int(row['descendant_branch'])}: "
                f"nearest FEM mode {int(row['nearest_fem_sorted_mode'])}, "
                f"same-index abs diff={float(row['abs_diff']):.6e}, "
                f"nearest abs diff={float(row['nearest_fem_abs_diff']):.6e}"
            )
        if len(nearest_mismatch) > 12:
            lines.append(f"  - plus {len(nearest_mismatch) - 12} additional rows.")
    if tracking_warnings:
        lines.append("- analytic tracking has warning rows; these are unresolved diagnostics, not branch renumbering:")
        for row in tracking_warnings[:12]:
            lines.append(
                f"  - mu={float(row['mu']):g}, desc={int(row['branch_index_from_mu0'])}, "
                f"accepted pos={int(row['mac_sorted_root_index'])}, "
                f"candidate pos={int(row['diagnostic_candidate_sorted_position'])}, "
                f"candidate MAC={float(row['diagnostic_candidate_mac_to_previous']):.6g}"
            )
        if len(tracking_warnings) > 12:
            lines.append(f"  - plus {len(tracking_warnings) - 12} additional rows.")

    lines.extend(["", "## Thin-Rod Criterion", ""])
    if summary.has_violations:
        for segment in summary.segments:
            lines.append(
                f"- WARNING eta={ETA:g}: criterion violated on mu={segment.start_mu:g}..{segment.end_mu:g}, "
                f"rod(s) {rods_label(segment.rods)}"
            )
        lines.append(
            f"- max ratio={summary.max_ratio:.6g} at mu={summary.max_ratio_mu:g}, rod {summary.max_ratio_rod}"
        )
    else:
        lines.append(
            f"- no violations of `2*r_i/l_i <= 0.1` on the FEM grid; "
            f"max ratio={summary.max_ratio:.6g} at mu={summary.max_ratio_mu:g}, rod {summary.max_ratio_rod}"
        )

    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
            "",
            "This is diagnostic-only. It does not modify article files, article",
            "figures, the old determinant, old solvers, or the existing FEM",
            "physical model.",
            "",
        ]
    )
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")
    return {
        "max_abs": float(abs_row["abs_diff"]),
        "max_abs_mu": float(abs_row["mu"]),
        "max_abs_mode": int(abs_row["descendant_branch"]),
        "max_rel": float(rel_row["rel_diff"]),
        "max_rel_mu": float(rel_row["mu"]),
        "max_rel_mode": int(rel_row["descendant_branch"]),
        "nearest_mismatch_count": len(nearest_mismatch),
        "tracking_warning_count": len(tracking_warnings),
        "max_thickness_ratio": summary.max_ratio,
        "thickness_violations": summary.has_violations,
    }


def main() -> dict[str, object]:
    tracked, tracking_rows = analytic_tracking()
    fem_by_mesh = compute_fem_by_mesh()
    rows = comparison_rows(tracked=tracked, fem_by_mesh=fem_by_mesh)
    plot_comparison(tracked=tracked, fem_by_mesh=fem_by_mesh)
    report = write_report(rows=rows, tracking_rows=tracking_rows)

    print(f"saved FEM comparison PNG: {OUTPUT_PNG}")
    print(f"saved FEM comparison report: {OUTPUT_REPORT}")
    print(
        f"max abs diff={report['max_abs']:.6e} at mu={report['max_abs_mu']:g}, "
        f"mode={report['max_abs_mode']}"
    )
    print(
        f"max rel diff={report['max_rel']:.6e} at mu={report['max_rel_mu']:g}, "
        f"mode={report['max_rel_mode']}"
    )
    print(f"nearest FEM mode mismatches: {report['nearest_mismatch_count']}")
    print(f"analytic tracking warning rows: {report['tracking_warning_count']}")
    if report["thickness_violations"]:
        print("WARNING: thickness-ratio criterion violations detected on FEM grid")
    else:
        print(f"no thickness-ratio violations on FEM grid; max={report['max_thickness_ratio']:.6g}")

    return {"png": OUTPUT_PNG, "report": OUTPUT_REPORT, **report}


if __name__ == "__main__":
    main()
