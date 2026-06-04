from __future__ import annotations

import csv
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_out_of_plane import (  # noqa: E402
    fixed_fixed_bending_roots_total_length_two,
    torsion_roots_uniform_eta0_beta0,
)
from my_project.analytic.out_of_plane_fem_1d import (  # noqa: E402
    assemble_out_of_plane_fem_1d_matrices,
    first_beta0_eta0_torsion_fem_root,
    out_of_plane_fem_1d_matrix_warnings,
    solve_out_of_plane_fem_1d_modes,
)
from my_project.analytic.solvers_out_of_plane import (  # noqa: E402
    find_first_n_roots_out_of_plane_with_warnings,
)


EPSILON = 0.0025
POISSON = 0.3
ETA_VALUES = (0.0, 0.1)
MU_VALUES = (0.0, 0.3)
BETA_DEG_VALUES = (0.0, 15.0, 45.0)
N_MODES = 8
MESHES = (4, 8, 16, 32)
LAMBDA_MAX = 25.0
ROOT_SCAN_STEP = 0.01

OUTPUT_DIR = REPO_ROOT / "results" / "out_of_plane_fem_validation"
COMPARISON_CSV = OUTPUT_DIR / "out_of_plane_analytic_vs_1d_fem_comparison.csv"
CONVERGENCE_CSV = OUTPUT_DIR / "out_of_plane_analytic_vs_1d_fem_convergence.csv"
REPORT_MD = OUTPUT_DIR / "out_of_plane_analytic_vs_1d_fem_report.md"


@dataclass(frozen=True)
class ValidationCase:
    beta_deg: float
    beta_rad: float
    mu: float
    eta: float
    epsilon: float
    poisson: float


def _fmt(value: float) -> str:
    if not isfinite(float(value)):
        return "nan"
    return f"{float(value):.16e}"


def _rel_error(estimate: float, reference: float) -> float:
    estimate_f = float(estimate)
    reference_f = float(reference)
    if not (isfinite(estimate_f) and isfinite(reference_f)) or abs(reference_f) <= 0.0:
        return float("nan")
    return abs(estimate_f - reference_f) / abs(reference_f)


def _warning_text(warnings: list[str]) -> str:
    return "; ".join(warnings) if warnings else "none"


def _cases() -> list[ValidationCase]:
    cases: list[ValidationCase] = []
    for eta in ETA_VALUES:
        for mu in MU_VALUES:
            for beta_deg in BETA_DEG_VALUES:
                cases.append(
                    ValidationCase(
                        beta_deg=float(beta_deg),
                        beta_rad=float(np.deg2rad(beta_deg)),
                        mu=float(mu),
                        eta=float(eta),
                        epsilon=EPSILON,
                        poisson=POISSON,
                    )
                )
    return cases


def _fem_modes_and_warnings(case: ValidationCase, mesh: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    warnings: list[str] = []
    try:
        stiffness, mass = assemble_out_of_plane_fem_1d_matrices(
            beta=case.beta_rad,
            mu=case.mu,
            epsilon=case.epsilon,
            eta=case.eta,
            poisson=case.poisson,
            n_elements_per_rod=mesh,
        )
        warnings.extend(out_of_plane_fem_1d_matrix_warnings(stiffness, mass))
        modes = solve_out_of_plane_fem_1d_modes(
            beta=case.beta_rad,
            mu=case.mu,
            epsilon=case.epsilon,
            eta=case.eta,
            poisson=case.poisson,
            n_elements_per_rod=mesh,
            n_modes=N_MODES,
        )
    except Exception as exc:  # pragma: no cover - retained for diagnostic CSV reporting.
        empty = np.array([], dtype=float)
        return empty, empty, empty, [f"FEM solve failed: {exc}"]
    return (
        np.asarray(modes.lambdas, dtype=float),
        np.asarray(modes.bending_energy_fraction, dtype=float),
        np.asarray(modes.torsion_energy_fraction, dtype=float),
        warnings,
    )


def _comparison_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in _cases():
        analytic_result = find_first_n_roots_out_of_plane_with_warnings(
            beta=case.beta_rad,
            mu=case.mu,
            epsilon=case.epsilon,
            eta=case.eta,
            poisson=case.poisson,
            n_roots=N_MODES,
            lambda_max=LAMBDA_MAX,
            scan_step=ROOT_SCAN_STEP,
        )
        analytic = np.asarray(analytic_result.roots, dtype=float)
        root_warnings = list(analytic_result.warnings)

        for mesh in MESHES:
            fem, bending_fraction, torsion_fraction, fem_warnings = _fem_modes_and_warnings(case, int(mesh))
            for mode_zero in range(N_MODES):
                analytic_value = float(analytic[mode_zero]) if mode_zero < len(analytic) else float("nan")
                fem_value = float(fem[mode_zero]) if mode_zero < len(fem) else float("nan")
                bending_value = (
                    float(bending_fraction[mode_zero]) if mode_zero < len(bending_fraction) else float("nan")
                )
                torsion_value = (
                    float(torsion_fraction[mode_zero]) if mode_zero < len(torsion_fraction) else float("nan")
                )
                notes: list[str] = []
                if mode_zero >= len(analytic):
                    notes.append("missing analytic root")
                if mode_zero >= len(fem):
                    notes.append("missing FEM root")
                if not notes:
                    notes.append("ok")
                rows.append(
                    {
                        "beta_deg": case.beta_deg,
                        "mu": case.mu,
                        "eta": case.eta,
                        "epsilon": case.epsilon,
                        "poisson": case.poisson,
                        "n_elements_per_rod": int(mesh),
                        "mode_index": mode_zero + 1,
                        "Lambda_analytic": analytic_value,
                        "Lambda_fem": fem_value,
                        "fem_bending_energy_fraction": bending_value,
                        "fem_torsion_energy_fraction": torsion_value,
                        "abs_error": abs(fem_value - analytic_value)
                        if isfinite(fem_value) and isfinite(analytic_value)
                        else float("nan"),
                        "rel_error": _rel_error(fem_value, analytic_value),
                        "notes": "; ".join(notes),
                        "root_solver_warnings": _warning_text(root_warnings),
                        "fem_warnings": _warning_text(fem_warnings),
                        "num_analytic_roots_found": len(analytic),
                        "num_fem_roots_found": len(fem),
                    }
                )
    return rows


def _write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out = dict(row)
            for key, value in list(out.items()):
                if isinstance(value, float):
                    out[key] = _fmt(value)
            writer.writerow(out)


def _convergence_rows(comparison_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_key: dict[tuple[float, float, float, int], dict[int, float]] = {}
    for row in comparison_rows:
        key = (float(row["beta_deg"]), float(row["mu"]), float(row["eta"]), int(row["mode_index"]))
        by_key.setdefault(key, {})[int(row["n_elements_per_rod"])] = float(row["rel_error"])

    rows: list[dict[str, object]] = []
    for (beta_deg, mu, eta, mode_index), errors in sorted(by_key.items()):
        values = [errors.get(int(mesh), float("nan")) for mesh in MESHES]
        finite_pairs = [
            (left, right)
            for left, right in zip(values[:-1], values[1:])
            if isfinite(float(left)) and isfinite(float(right))
        ]
        decrease_count = sum(1 for left, right in finite_pairs if right <= left)
        final_le_initial = (
            isfinite(values[0]) and isfinite(values[-1]) and values[-1] <= values[0]
        )
        strict_monotone = bool(finite_pairs) and all(right <= left for left, right in finite_pairs)
        suspicious = mode_index <= 4 and not final_le_initial
        rows.append(
            {
                "beta_deg": beta_deg,
                "mu": mu,
                "eta": eta,
                "epsilon": EPSILON,
                "poisson": POISSON,
                "mode_index": mode_index,
                "rel_error_mesh_4": values[0],
                "rel_error_mesh_8": values[1],
                "rel_error_mesh_16": values[2],
                "rel_error_mesh_32": values[3],
                "adjacent_error_decrease_count": decrease_count,
                "final_error_not_larger_than_mesh4": "yes" if final_le_initial else "no",
                "strict_monotone_decrease": "yes" if strict_monotone else "no",
                "notes": "suspicious low-mode non-convergence" if suspicious else "ok",
            }
        )
    return rows


def _unique_non_none(rows: list[dict[str, object]], key: str) -> list[str]:
    values = sorted({str(row[key]) for row in rows if str(row[key]) != "none"})
    return values


def _straight_beam_summary(comparison_rows: list[dict[str, object]]) -> dict[str, float]:
    exact_bending = fixed_fixed_bending_roots_total_length_two()
    selected = [
        row
        for row in comparison_rows
        if int(row["n_elements_per_rod"]) == max(MESHES)
        and abs(float(row["beta_deg"])) <= 1e-12
        and abs(float(row["eta"])) <= 1e-12
        and int(row["mode_index"]) <= len(exact_bending)
    ]
    analytic_errors = []
    fem_errors = []
    for row in selected:
        exact = float(exact_bending[int(row["mode_index"]) - 1])
        analytic_errors.append(_rel_error(float(row["Lambda_analytic"]), exact))
        fem_errors.append(_rel_error(float(row["Lambda_fem"]), exact))
    return {
        "analytic_max_rel_error_first3_bending": max(analytic_errors) if analytic_errors else float("nan"),
        "fem_max_rel_error_first3_bending": max(fem_errors) if fem_errors else float("nan"),
        "torsion_root_1": torsion_roots_uniform_eta0_beta0(1, EPSILON, POISSON),
    }


def _targeted_torsion_summary() -> dict[str, float]:
    analytic = torsion_roots_uniform_eta0_beta0(1, EPSILON, POISSON)
    fem = first_beta0_eta0_torsion_fem_root(
        epsilon=EPSILON,
        poisson=POISSON,
        n_elements_per_rod=max(MESHES),
    )
    return {
        "analytic_torsion_root_1": float(analytic),
        "fem_torsion_root_1": float(fem),
        "relative_error": _rel_error(fem, analytic),
        "n_elements_per_rod": float(max(MESHES)),
    }


def _write_report(comparison_rows: list[dict[str, object]], convergence_rows: list[dict[str, object]]) -> dict[str, object]:
    finest = [
        row
        for row in comparison_rows
        if int(row["n_elements_per_rod"]) == max(MESHES) and isfinite(float(row["rel_error"]))
    ]
    rel_errors = np.asarray([float(row["rel_error"]) for row in finest], dtype=float)
    max_rel = float(np.max(rel_errors)) if len(rel_errors) else float("nan")
    median_rel = float(np.median(rel_errors)) if len(rel_errors) else float("nan")
    low_mode_convergence = [
        row for row in convergence_rows if int(row["mode_index"]) <= 4
    ]
    low_mode_decreasing = sum(
        1 for row in low_mode_convergence if row["final_error_not_larger_than_mesh4"] == "yes"
    )
    suspicious = [row for row in convergence_rows if str(row["notes"]).startswith("suspicious")]
    straight = _straight_beam_summary(comparison_rows)
    torsion = _targeted_torsion_summary()
    root_warnings = _unique_non_none(comparison_rows, "root_solver_warnings")
    fem_warnings = _unique_non_none(comparison_rows, "fem_warnings")
    default_finest = [
        row
        for row in comparison_rows
        if int(row["n_elements_per_rod"]) == max(MESHES)
        and isfinite(float(row["fem_bending_energy_fraction"]))
        and isfinite(float(row["fem_torsion_energy_fraction"]))
    ]
    bending_dominated_count = sum(
        1
        for row in default_finest
        if float(row["fem_bending_energy_fraction"]) >= float(row["fem_torsion_energy_fraction"])
    )

    lines = [
        "# Out-of-Plane Analytic vs 1D FEM Validation",
        "",
        "This diagnostic compares the out-of-plane Euler--Bernoulli plus",
        "Saint-Venant torsion determinant with an independent 1D finite element",
        "discretization using the same continuum assumptions.",
        "",
        "It is not a 3D FEM validation. No Gmsh, CalculiX, solid elements, or",
        "article figures are used here.",
        "",
        "## Parameters",
        "",
        f"- epsilon: {EPSILON:g}",
        f"- poisson: {POISSON:g}",
        f"- eta values: {', '.join(f'{value:g}' for value in ETA_VALUES)}",
        f"- mu values: {', '.join(f'{value:g}' for value in MU_VALUES)}",
        f"- beta values: {', '.join(f'{value:g}' for value in BETA_DEG_VALUES)} deg",
        f"- modes compared: first {N_MODES} sorted Lambda values",
        f"- elements per rod: {', '.join(str(value) for value in MESHES)}",
        "",
        "## Scaling",
        "",
        "The FEM uses `l0=E=rho=S0=1`, `J0=epsilon^2`, `Jp0=2*J0`,",
        "and `G=E/(2*(1+poisson))`. Its generalized eigenproblem is",
        "`K q = omega^2 M q`, converted to the analytic scale by",
        "`Lambda = (omega^2 / epsilon^2)^0.25`.",
        "",
        "## Finest Mesh Summary",
        "",
        f"- finest mesh: {max(MESHES)} elements per rod",
        f"- max relative error: {max_rel:.6e}",
        f"- median relative error: {median_rel:.6e}",
        "",
        "## beta=0, eta=0 Straight-Beam Checks",
        "",
        "- first three bending roots are checked against",
        "  `cosh(2 Lambda) cos(2 Lambda) - 1 = 0`.",
        f"- analytic max relative error, first three bending roots: "
        f"{straight['analytic_max_rel_error_first3_bending']:.6e}",
        f"- 1D FEM max relative error, first three bending roots on finest mesh: "
        f"{straight['fem_max_rel_error_first3_bending']:.6e}",
        f"- first closed-form torsion root for these parameters: "
        f"{straight['torsion_root_1']:.12g}",
        "- for the default `epsilon=0.0025`, that torsion root lies above the",
        "  first eight sorted modes compared in the CSV.",
        "- those default first-eight FEM rows are bending-dominated by the",
        f"  stiffness-energy diagnostic in {bending_dominated_count}/{len(default_finest)}",
        "  finest-mesh rows.",
        "",
        "## Targeted Torsion Check",
        "",
        "- a separate beta=0, eta=0 torsion-block FEM check was run so the",
        "  torsional family is visible even though it lies above the default",
        "  first eight sorted modes.",
        f"- analytic first torsion root: {torsion['analytic_torsion_root_1']:.12g}",
        f"- FEM first torsion root, {int(torsion['n_elements_per_rod'])} elements per rod: "
        f"{torsion['fem_torsion_root_1']:.12g}",
        f"- relative error: {torsion['relative_error']:.6e}",
        "",
        "## FEM Mode-Character Diagnostic",
        "",
        "- comparison CSV columns `fem_bending_energy_fraction` and",
        "  `fem_torsion_energy_fraction` record each FEM mode's split by local",
        "  EB bending stiffness energy and Saint-Venant torsion stiffness energy.",
        "",
        "## Mesh Convergence",
        "",
        f"- low-mode rows with final error not larger than mesh 4: "
        f"{low_mode_decreasing}/{len(low_mode_convergence)}",
        f"- suspicious low-mode non-convergence flags: {len(suspicious)}",
        "- strict monotone decrease is recorded in the convergence CSV but is not",
        "  required for every mode.",
        "",
        "## Warnings",
        "",
        f"- analytic root-finder warnings: {len(root_warnings)}",
    ]
    lines.extend([f"  - {warning}" for warning in root_warnings] or ["  - none"])
    lines.append(f"- FEM warnings: {len(fem_warnings)}")
    lines.extend([f"  - {warning}" for warning in fem_warnings] or ["  - none"])
    lines.extend(
        [
            "",
            "## Future 3D FEM Note",
            "",
            "After the 1D FEM validation is stable, a second-level comparison with",
            "3D FEM may be attempted. That comparison must be treated separately",
            "because 3D FEM includes joint geometry, shear/rotary inertia, possible",
            "warping, Poisson effects, mesh sensitivity, and possible mixing between",
            "in-plane and out-of-plane modes if the model is not perfectly symmetric.",
            "",
            "## Outputs",
            "",
            f"- `{COMPARISON_CSV.relative_to(REPO_ROOT)}`",
            f"- `{CONVERGENCE_CSV.relative_to(REPO_ROOT)}`",
            f"- `{REPORT_MD.relative_to(REPO_ROOT)}`",
            "",
        ]
    )
    REPORT_MD.write_text("\n".join(lines), encoding="utf-8")
    return {
        "max_rel_error_finest": max_rel,
        "median_rel_error_finest": median_rel,
        "straight": straight,
        "torsion": torsion,
        "bending_dominated_finest_rows": bending_dominated_count,
        "finest_rows": len(default_finest),
        "root_warning_count": len(root_warnings),
        "fem_warning_count": len(fem_warnings),
        "suspicious_convergence_count": len(suspicious),
    }


def main() -> dict[str, object]:
    comparison = _comparison_rows()
    convergence = _convergence_rows(comparison)
    comparison_fields = [
        "beta_deg",
        "mu",
        "eta",
        "epsilon",
        "poisson",
        "n_elements_per_rod",
        "mode_index",
        "Lambda_analytic",
        "Lambda_fem",
        "fem_bending_energy_fraction",
        "fem_torsion_energy_fraction",
        "abs_error",
        "rel_error",
        "notes",
        "root_solver_warnings",
        "fem_warnings",
        "num_analytic_roots_found",
        "num_fem_roots_found",
    ]
    convergence_fields = [
        "beta_deg",
        "mu",
        "eta",
        "epsilon",
        "poisson",
        "mode_index",
        "rel_error_mesh_4",
        "rel_error_mesh_8",
        "rel_error_mesh_16",
        "rel_error_mesh_32",
        "adjacent_error_decrease_count",
        "final_error_not_larger_than_mesh4",
        "strict_monotone_decrease",
        "notes",
    ]
    _write_csv(COMPARISON_CSV, comparison, comparison_fields)
    _write_csv(CONVERGENCE_CSV, convergence, convergence_fields)
    summary = _write_report(comparison, convergence)

    print(f"saved comparison CSV: {COMPARISON_CSV}")
    print(f"saved convergence CSV: {CONVERGENCE_CSV}")
    print(f"saved report: {REPORT_MD}")
    print(f"max relative error on finest mesh: {summary['max_rel_error_finest']:.6e}")
    print(f"median relative error on finest mesh: {summary['median_rel_error_finest']:.6e}")
    print(
        "targeted torsion root: "
        f"analytic={summary['torsion']['analytic_torsion_root_1']:.12g}, "
        f"fem={summary['torsion']['fem_torsion_root_1']:.12g}, "
        f"rel_error={summary['torsion']['relative_error']:.6e}"
    )
    print(
        "bending-dominated first-eight finest rows: "
        f"{summary['bending_dominated_finest_rows']}/{summary['finest_rows']}"
    )
    print(f"analytic root-finder warnings: {summary['root_warning_count']}")
    print(f"FEM warnings: {summary['fem_warning_count']}")
    print(f"suspicious convergence flags: {summary['suspicious_convergence_count']}")
    return summary


if __name__ == "__main__":
    main()
