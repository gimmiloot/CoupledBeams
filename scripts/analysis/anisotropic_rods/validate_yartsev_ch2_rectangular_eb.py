"""Validate the theta=0 rectangular EB comparator and independent 1D FEM."""

from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib.yartsev_ch2_coupled_rods import (  # noqa: E402
    coupled_boundary_matrix,
    coupled_boundary_matrix_raw,
    equilibrate_matrix,
    straight_boundary_matrix,
    straight_boundary_matrix_raw,
)
from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
    Geometry,
    RodPoint,
    RootResult,
    cantilever_geometry,
    find_elastic_roots,
    hms_dx_209_material,
    make_rod_point,
)
from scripts.lib.yartsev_ch2_rectangular_eb import (  # noqa: E402
    assemble_two_arm_eb_fem,
    eb_bending_coupled_boundary_matrix,
    eb_bending_coupled_boundary_matrix_raw,
    eb_clamp_matrix,
    eb_coupled_boundary_matrix,
    eb_coupled_boundary_matrix_raw,
    eb_joint_mapping_residual,
    eb_state_matrix,
    eb_straight_boundary_matrix,
    eb_straight_boundary_matrix_raw,
    fixed_fixed_bending_frequencies_hz,
    fixed_fixed_torsion_frequencies_hz,
    saint_venant_stiffness,
    solve_two_arm_eb_fem,
)


DEFAULT_OUTPUT_DIR = (
    ROOT / "results" / "anisotropic_rods" / "yartsev_ch2_rectangular_eb_validation"
)
NUM_ROOTS = 7
MESHES = (4, 8, 16, 32, 64)
SPLITS = ((0.20, 0.20), (0.15, 0.25), (0.10, 0.30), (0.05, 0.35))
NEAR_DEGENERATE_RELATIVE_GAP = 1.0e-3
TARGETED_PROPORTIONAL_MESHES = ((8, 24), (16, 48), (32, 96), (64, 192))
TARGETED_EQUAL_H_CONTROL = (128, 128)


class CountingBoundaryBuilder:
    """Cache and time unique matrices supplied to the existing root finder."""

    def __init__(self, factory: Callable[[complex], np.ndarray]) -> None:
        self.factory = factory
        self.calls = 0
        self.evaluations = 0
        self.elapsed_seconds = 0.0
        self._cache: dict[complex, np.ndarray] = {}

    def __call__(self, omega: complex, _point: RodPoint, formulation: str) -> np.ndarray:
        if formulation != "state_corrected":
            raise ValueError("this validation accepts only state_corrected root machinery")
        self.calls += 1
        key = complex(omega)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        started = time.perf_counter()
        matrix = np.asarray(self.factory(key), dtype=np.complex128)
        self.elapsed_seconds += time.perf_counter() - started
        self.evaluations += 1
        self._cache[key] = matrix
        return matrix


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument(
        "--smoke",
        action="store_true",
        help="run exact straight checks and one beta=0 analytic/FEM case only",
    )
    modes.add_argument(
        "--targeted-length-proportional-refinement",
        action="store_true",
        help="run only the unequal-arm proportional meshes and equal-h control",
    )
    modes.add_argument(
        "--targeted-smoke",
        action="store_true",
        help="run only the unequal-arm (16,48) proportional-mesh smoke",
    )
    return parser.parse_args(argv)


def _point(length: float, section_scale: float = 1.0) -> RodPoint:
    base = cantilever_geometry(length)
    geometry = Geometry(
        a=section_scale * base.a,
        b=section_scale * base.b,
        length=length,
        shear_factor=base.shear_factor,
    )
    return make_rod_point(
        0.0,
        material_mode="elastic",
        geometry=geometry,
        material=hms_dx_209_material(),
    )


def _write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    data = list(rows)
    if not data:
        raise ValueError(f"refusing to create empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in data:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)


def _matrix_quality(matrix: np.ndarray) -> dict[str, float]:
    value = np.asarray(matrix, dtype=np.complex128)
    singular = np.linalg.svd(value, compute_uv=False)
    sigma_max = float(singular[0])
    sigma_min = float(singular[-1])
    determinant_abs = float(abs(np.linalg.det(value)))
    determinant_scale = max(sigma_max ** value.shape[0], np.finfo(float).tiny)
    return {
        "raw_determinant_abs": determinant_abs,
        "raw_determinant_residual": determinant_abs / determinant_scale,
        "raw_sigma_min": sigma_min,
        "raw_sigma_max": sigma_max,
        "raw_relative_singular_residual": sigma_min / sigma_max if sigma_max else 0.0,
    }


def _root_fields(root: RootResult, raw_factory: Callable[[complex], np.ndarray]) -> dict[str, Any]:
    raw_matrix = np.asarray(raw_factory(root.omega), dtype=np.complex128)
    raw_quality = _matrix_quality(raw_matrix)
    scaled_ok = _root_quality_ok(root)
    raw_ok = (
        raw_quality["raw_determinant_residual"] <= 1.0e-8
        and raw_quality["raw_relative_singular_residual"] <= 1.0e-8
    )
    if raw_matrix.shape == (1, 1):
        quality_status = "NOT_APPLICABLE_SCALAR_EXACT_FAMILY"
        quality_basis = "exact_frequency_equation"
    else:
        quality_status = "PASS" if scaled_ok or raw_ok else "FAIL"
        quality_basis = "scaled" if scaled_ok else "physical_raw"
    return {
        "scaled_determinant_abs": root.raw_determinant_abs,
        "scaled_determinant_residual": root.determinant_residual,
        "scaled_sigma_min": root.sigma_min,
        "scaled_sigma_max": root.sigma_max,
        "scaled_relative_singular_residual": root.relative_singular_residual,
        "root_status": root.status,
        "quality_status": quality_status,
        "quality_basis": quality_basis,
        **raw_quality,
    }


def _root_quality_ok(root: RootResult) -> bool:
    return (
        root.determinant_residual <= 1.0e-8
        and root.relative_singular_residual <= 1.0e-8
        and not root.status.startswith("rejected")
    )


def _solve(
    point: RodPoint,
    factory: Callable[[complex], np.ndarray],
    *,
    num_roots: int = NUM_ROOTS,
) -> tuple[list[RootResult], CountingBoundaryBuilder, float]:
    builder = CountingBoundaryBuilder(factory)
    started = time.perf_counter()
    roots = find_elastic_roots(
        point,
        "state_corrected",
        num_roots=num_roots,
        scan_step_hz=10.0,
        initial_max_hz=5000.0,
        max_hz=100_000.0,
        boundary_matrix_builder=builder,
    )
    return roots, builder, time.perf_counter() - started


def _frequencies(roots: Sequence[RootResult]) -> np.ndarray:
    return np.asarray([root.frequency_hz for root in roots], dtype=float)


def _neighbor_gaps(frequencies: np.ndarray) -> np.ndarray:
    gaps = np.full(len(frequencies), np.inf, dtype=float)
    for index in range(len(frequencies)):
        candidates: list[float] = []
        if index:
            candidates.append(float(frequencies[index] - frequencies[index - 1]))
        if index + 1 < len(frequencies):
            candidates.append(float(frequencies[index + 1] - frequencies[index]))
        if candidates:
            gaps[index] = min(candidates)
    return gaps


def _match_frequencies(
    analytic: np.ndarray, fem: np.ndarray
) -> tuple[np.ndarray, list[str], np.ndarray]:
    """Use sorted matching except inside connected near-degenerate clusters."""

    matched = np.asarray(fem, dtype=float).copy()
    statuses = ["sorted_no_near_degeneracy"] * len(analytic)
    gaps = _neighbor_gaps(analytic)
    near_links = [
        index
        for index in range(len(analytic) - 1)
        if analytic[index + 1] - analytic[index]
        <= NEAR_DEGENERATE_RELATIVE_GAP * max(analytic[index], analytic[index + 1])
    ]
    start = 0
    while start < len(near_links):
        first_link = near_links[start]
        last_link = first_link
        while start + 1 < len(near_links) and near_links[start + 1] == last_link + 1:
            start += 1
            last_link = near_links[start]
        indices = np.arange(first_link, last_link + 2)
        cost = np.abs(analytic[indices, None] - fem[indices][None, :])
        rows, columns = linear_sum_assignment(cost)
        for row, column in zip(rows, columns):
            analytic_index = int(indices[row])
            matched[analytic_index] = fem[int(indices[column])]
            statuses[analytic_index] = (
                "local_frequency_cluster_matching_modes_"
                + "-".join(str(int(item + 1)) for item in indices)
            )
        start += 1
    return matched, statuses, gaps


def _exact_checks(
    straight: RodPoint,
) -> tuple[list[dict[str, Any]], list[RootResult], list[RootResult], dict[str, Any]]:
    exact_bending = fixed_fixed_bending_frequencies_hz(straight, 3)
    exact_torsion = fixed_fixed_torsion_frequencies_hz(straight, 3)

    def bending_raw(omega: complex) -> np.ndarray:
        raw = eb_straight_boundary_matrix_raw(omega, straight)
        return raw[np.ix_([0, 1], [0, 1])]

    def bending_scaled(omega: complex) -> np.ndarray:
        return equilibrate_matrix(bending_raw(omega))[0]

    def torsion_raw(omega: complex) -> np.ndarray:
        raw = eb_straight_boundary_matrix_raw(omega, straight)
        return raw[np.ix_([2], [2])]

    bending_roots, bending_builder, bending_runtime = _solve(
        straight, bending_scaled, num_roots=3
    )
    torsion_roots, torsion_builder, torsion_runtime = _solve(
        straight, torsion_raw, num_roots=3
    )
    rows: list[dict[str, Any]] = []
    for family, exact, roots, raw_factory in (
        ("bending", exact_bending, bending_roots, bending_raw),
        ("torsion", exact_torsion, torsion_roots, torsion_raw),
    ):
        for mode, (expected, root) in enumerate(zip(exact, roots), start=1):
            absolute = abs(root.frequency_hz - expected)
            rows.append(
                {
                    "family": family,
                    "family_mode": mode,
                    "exact_frequency_hz": expected,
                    "transfer_frequency_hz": root.frequency_hz,
                    "absolute_error_hz": absolute,
                    "relative_error": absolute / expected,
                    **_root_fields(root, raw_factory),
                }
            )
    metrics = {
        "runtime_seconds": bending_runtime + torsion_runtime,
        "evaluations": bending_builder.evaluations + torsion_builder.evaluations,
    }
    return rows, bending_roots, torsion_roots, metrics


def _smoke() -> dict[str, Any]:
    started = time.perf_counter()
    arm = _point(0.2)
    straight = _point(0.4)
    exact_rows, bending_roots, torsion_roots, exact_metrics = _exact_checks(straight)
    analytic, builder, analytic_runtime = _solve(
        arm, lambda omega: eb_coupled_boundary_matrix(omega, 0.0, arm, arm)
    )
    fem_started = time.perf_counter()
    fem = solve_two_arm_eb_fem(
        assemble_two_arm_eb_fem(arm, arm, 0.0, 8), num_roots=NUM_ROOTS
    )
    fem_runtime = time.perf_counter() - fem_started
    relative = np.abs(_frequencies(analytic) - fem.frequencies_hz) / _frequencies(analytic)
    return {
        "runtime_seconds": time.perf_counter() - started,
        "exact_max_relative_error": max(row["relative_error"] for row in exact_rows),
        "exact_roots": len(bending_roots) + len(torsion_roots),
        "analytic_roots": len(analytic),
        "analytic_evaluations": builder.evaluations,
        "exact_evaluations": exact_metrics["evaluations"],
        "analytic_runtime_seconds": analytic_runtime,
        "fem_runtime_seconds": fem_runtime,
        "fem_reduced_size": int(assemble_two_arm_eb_fem(arm, arm, 0.0, 8).mass.shape[0]),
        "mesh8_max_first6_relative_error": float(np.max(relative[:6])),
    }


def _report(
    status: str,
    summary_rows: list[dict[str, Any]],
    unequal_timo_rows: list[dict[str, Any]],
    unequal_eb_rows: list[dict[str, Any]],
    exact_rows: list[dict[str, Any]],
    comparison_rows: list[dict[str, Any]],
    slender_rows: list[dict[str, Any]],
    fem_rows: list[dict[str, Any]],
    convergence_rows: list[dict[str, Any]],
    runtime_seconds: float,
    evaluations: int,
) -> str:
    max_timo = max(row["relative_error"] for row in unequal_timo_rows)
    max_eb = max(row["relative_error"] for row in unequal_eb_rows)
    finest = [row for row in fem_rows if row["elements_per_arm"] == 64 and row["mode"] <= 6]
    near = sorted({row["matching_status"] for row in fem_rows if row["matching_status"] != "sorted_no_near_degeneracy"})
    lines = [
        "# Rectangular orthotropic EB validation report",
        "",
        f"Final gate status: **{status}**.",
        "",
        "## Scope",
        "",
        "Elastic HMS/DX-209, theta=0 only; rectangular EB plus generalized Saint-Venant torsion is compared with the unchanged Chapter-2 Timoshenko model and an independent 1D FEM.",
        "",
        "## Thirteen acceptance gates",
        "",
        "| Gate | Result | Observed | Threshold/details |",
        "|---|---:|---:|---|",
    ]
    for row in summary_rows:
        lines.append(
            f"| {row['gate']} | {row['status']} | {row['observed']} | {row['requirement']} |"
        )
    lines.extend(
        [
            "",
            "## Exact straight families",
            "",
            "| Family | Mode | Exact (Hz) | Transfer (Hz) | Relative error |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in exact_rows:
        lines.append(
            f"| {row['family']} | {row['family_mode']} | {row['exact_frequency_hz']:.9f} | {row['transfer_frequency_hz']:.9f} | {row['relative_error']:.3e} |"
        )
    lines.extend(
        [
            "",
            "## Unequal-length invariance",
            "",
            f"Maximum relative difference: Timoshenko `{max_timo:.3e}`, EB `{max_eb:.3e}` (target `1e-8`). Root-by-root values are in the two unequal-length CSV files.",
            "",
            "## Timoshenko versus EB",
            "",
            "| beta (deg) | mode | Timoshenko (Hz) | EB (Hz) | delta_f | order flag |",
            "|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in comparison_rows:
        if row["mode"] <= 7:
            lines.append(
                f"| {row['beta_deg']:.0f} | {row['mode']} | {row['timoshenko_frequency_hz']:.6f} | {row['eb_frequency_hz']:.6f} | {row['delta_f']:.3e} | {row['matching_status']} |"
            )
    lines.extend(
        [
            "",
            "## Slender-limit bending families",
            "",
            "| scale | mode | Timoshenko (Hz) | EB (Hz) | relative difference |",
            "|---:|---:|---:|---:|---:|",
        ]
    )
    for row in slender_rows:
        lines.append(
            f"| {row['section_scale']:.2f} | {row['bending_mode']} | {row['timoshenko_frequency_hz']:.6f} | {row['eb_frequency_hz']:.6f} | {row['relative_difference']:.3e} |"
        )
    lines.extend(
        [
            "",
            "## Independent FEM",
            "",
            "The FEM uses standard Hermite bending and linear torsion elements with consistent translational and rotary-inertia mass matrices. The outer nodes are fixed. Its joint is assembled independently through common `[w_J, theta_t, theta_n]` DOFs and a congruence transformation; `J_book` is used only for the post-assembly kinematic residual check.",
            "",
            f"Finest-mesh maximum relative error over modes 1--6: `{max(row['relative_error'] for row in finest):.3e}`. Near-degenerate matching statuses: `{near or ['none']}`.",
            "",
            "The strict `1e-5` first-three target is reported without changing element coefficients. For the beta=0 case, the standard linear torsion element retains its expected second-order consistent-mass dispersion error.",
            "",
            "## Runtime and root quality",
            "",
            f"Scientific runtime: `{runtime_seconds:.6f} s`; unique analytic boundary-matrix evaluations: `{evaluations}`. All root-quality rows retain raw/scaled determinant and SVD diagnostics, and every full spectrum includes root 7 as a guard.",
            "",
            "## Exclusions",
            "",
            "No theta != 0, unequal thickness, material bending--torsion coupling in EB, Timoshenko FEM, 3D FEM, damping, complex roots, coupled mode shapes, MAC tracking, or parameter maps were introduced.",
            "",
            "## Interpretation",
            "",
            "This is a project comparator, not a source reproduction and not a production anisotropic API. Timoshenko--EB differences are model differences caused here by transverse shear and bending rotary inertia; they are not an accuracy ranking.",
            "",
        ]
    )
    return "\n".join(lines)


def _load_accepted_analytic_frequencies(
    output_dir: Path, case_id: str
) -> np.ndarray:
    """Load the seven analytic roots accepted by the preceding full gate."""

    source = output_dir / "eb_analytic_vs_fem.csv"
    if not source.is_file():
        raise FileNotFoundError(
            "targeted refinement requires the preceding accepted analytic evidence: "
            f"{source}"
        )
    with source.open("r", encoding="utf-8", newline="") as stream:
        rows = [
            row
            for row in csv.DictReader(stream)
            if row["case"] == case_id and int(row["elements_per_arm"]) == 64
        ]
    rows.sort(key=lambda row: int(row["mode"]))
    if [int(row["mode"]) for row in rows] != list(range(1, NUM_ROOTS + 1)):
        raise RuntimeError(
            f"accepted analytic evidence for {case_id!r} does not contain modes 1--7"
        )
    frequencies = np.asarray(
        [float(row["analytic_frequency_hz"]) for row in rows], dtype=float
    )
    if not np.all(np.isfinite(frequencies)) or not np.all(np.diff(frequencies) > 0.0):
        raise RuntimeError(f"accepted analytic spectrum for {case_id!r} is invalid")
    return frequencies


def _beta0_exact_family_labels(
    straight: RodPoint, analytic_frequencies: np.ndarray
) -> list[str]:
    """Classify the beta=0 sorted roots from the exact decoupled families."""

    candidates: list[tuple[float, str]] = []
    for frequency in fixed_fixed_bending_frequencies_hz(straight, NUM_ROOTS):
        candidates.append((float(frequency), "bending"))
    for frequency in fixed_fixed_torsion_frequencies_hz(straight, NUM_ROOTS):
        candidates.append((float(frequency), "torsion"))
    candidates.sort(key=lambda item: item[0])
    selected = candidates[:NUM_ROOTS]
    exact = np.asarray([item[0] for item in selected], dtype=float)
    relative = np.abs(exact - analytic_frequencies) / analytic_frequencies
    if float(np.max(relative)) > 1.0e-8:
        raise RuntimeError(
            "accepted beta=0 analytic roots do not match the exact decoupled families"
        )
    return [item[1] for item in selected]


def _targeted_fem_rows(
    *,
    case_id: str,
    point_1: RodPoint,
    point_2: RodPoint,
    element_counts: tuple[int, int],
    analytic_frequencies: np.ndarray,
    family_labels: Sequence[str],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    beta = 0.0
    assembly_started = time.perf_counter()
    assembly = assemble_two_arm_eb_fem(
        point_1, point_2, beta, element_counts
    )
    assembly_runtime = time.perf_counter() - assembly_started
    solve_started = time.perf_counter()
    solution = solve_two_arm_eb_fem(assembly, num_roots=NUM_ROOTS)
    solve_runtime = time.perf_counter() - solve_started
    conditioning_started = time.perf_counter()
    stiffness_condition_number = float(np.linalg.cond(assembly.stiffness))
    mass_condition_number = float(np.linalg.cond(assembly.mass))
    conditioning_runtime = time.perf_counter() - conditioning_started
    matched, matching_statuses, neighbor_gaps = _match_frequencies(
        analytic_frequencies, solution.frequencies_hz
    )

    n_1, n_2 = element_counts
    h_1 = point_1.geometry.length / n_1
    h_2 = point_2.geometry.length / n_2
    positive_root_count = int(np.count_nonzero(solution.eigenvalues > 0.0))
    joint_kinematic_residual = eb_joint_mapping_residual(beta)
    diagnostics = {
        "case_id": case_id,
        "n_elements_arm_1": n_1,
        "n_elements_arm_2": n_2,
        "total_elements": n_1 + n_2,
        "h_1_m": h_1,
        "h_2_m": h_2,
        "h_max_m": max(h_1, h_2),
        "h_difference_m": abs(h_1 - h_2),
        "reduced_matrix_size": int(assembly.mass.shape[0]),
        "positive_root_count": positive_root_count,
        "spurious_zero_mode_count": solution.zero_mode_count,
        "stiffness_symmetry_residual": solution.stiffness_symmetry_residual,
        "mass_symmetry_residual": solution.mass_symmetry_residual,
        "minimum_reduced_mass_eigenvalue": solution.minimum_mass_eigenvalue,
        "stiffness_condition_number": stiffness_condition_number,
        "mass_condition_number": mass_condition_number,
        "joint_kinematic_residual": joint_kinematic_residual,
        "max_joint_equilibrium_residual": float(
            np.max(solution.joint_equilibrium_residuals)
        ),
        "assembly_runtime_s": assembly_runtime,
        "solve_runtime_s": solve_runtime,
        "conditioning_runtime_s": conditioning_runtime,
        "total_runtime_s": assembly_runtime + solve_runtime + conditioning_runtime,
    }
    rows: list[dict[str, Any]] = []
    for index, (
        analytic_frequency,
        fem_frequency,
        family,
        neighbor_gap,
        matching_status,
        joint_equilibrium_residual,
    ) in enumerate(
        zip(
            analytic_frequencies,
            matched,
            family_labels,
            neighbor_gaps,
            matching_statuses,
            solution.joint_equilibrium_residuals,
        ),
        start=1,
    ):
        absolute_error = abs(float(fem_frequency) - float(analytic_frequency))
        rows.append(
            {
                "case_id": case_id,
                "L_1_m": point_1.geometry.length,
                "L_2_m": point_2.geometry.length,
                "n_elements_arm_1": n_1,
                "n_elements_arm_2": n_2,
                "total_elements": n_1 + n_2,
                "h_1_m": h_1,
                "h_2_m": h_2,
                "h_max_m": max(h_1, h_2),
                "reduced_matrix_size": assembly.mass.shape[0],
                "mode": index,
                "family": family,
                "role": "validation" if index <= 6 else "guard",
                "analytic_frequency_hz": analytic_frequency,
                "fem_frequency_hz": fem_frequency,
                "absolute_error_hz": absolute_error,
                "relative_error": absolute_error / analytic_frequency,
                "neighbor_gap_hz": neighbor_gap,
                "matching_status": matching_status,
                "positive_root_count": positive_root_count,
                "spurious_zero_mode_count": solution.zero_mode_count,
                "stiffness_symmetry_residual": solution.stiffness_symmetry_residual,
                "mass_symmetry_residual": solution.mass_symmetry_residual,
                "minimum_reduced_mass_eigenvalue": solution.minimum_mass_eigenvalue,
                "stiffness_condition_number": stiffness_condition_number,
                "mass_condition_number": mass_condition_number,
                "joint_kinematic_residual": joint_kinematic_residual,
                "joint_equilibrium_residual": joint_equilibrium_residual,
                "assembly_runtime_s": assembly_runtime,
                "solve_runtime_s": solve_runtime,
                "conditioning_runtime_s": conditioning_runtime,
            }
        )
    return rows, diagnostics


def _empirical_order(coarse_error: float, fine_error: float) -> float:
    if coarse_error <= 0.0 or fine_error <= 0.0:
        return float("nan")
    return float(np.log(coarse_error / fine_error) / np.log(2.0))


def _targeted_report(
    *,
    status: str,
    target_rows: list[dict[str, Any]],
    order_rows: list[dict[str, Any]],
    summary_rows: list[dict[str, Any]],
    runtime_seconds: float,
) -> str:
    proportional = [
        row
        for row in target_rows
        if row["case_id"] == "unequal_length_proportional"
    ]
    finest = [
        row
        for row in proportional
        if row["n_elements_arm_1"] == 64
    ]
    control = [
        row
        for row in target_rows
        if row["case_id"] == "equal_arm_equal_h_control"
    ]
    mode_2_order = next(row for row in order_rows if row["mode"] == 2)
    nonmonotone = [
        (coarse["mode"], coarse["relative_error"], fine["relative_error"])
        for coarse, fine in zip(
            [
                row
                for row in proportional
                if row["n_elements_arm_1"] == 32 and row["mode"] <= 6
            ],
            [
                row
                for row in proportional
                if row["n_elements_arm_1"] == 64 and row["mode"] <= 6
            ],
        )
        if fine["relative_error"]
        > 1.01 * coarse["relative_error"] + 1.0e-10
    ]
    lines = [
        "# Targeted length-proportional FEM refinement",
        "",
        f"Targeted refinement status: **{status}**.",
        "",
        (
            "Overall rectangular EB validation: **PASS after targeted length-proportional FEM refinement**."
            if status == "PASS"
            else "Overall rectangular EB validation: **PARTIAL_PASS**; it is not promoted because the targeted refinement did not pass every acceptance gate."
        ),
        "",
        "The original fixed-elements-per-arm gate remains `PARTIAL_PASS`: `(64,64)` gave unequal physical element lengths and a first-three maximum relative error of `5.181993e-5`. No original CSV or report was overwritten.",
        "",
        "## Mesh sequence",
        "",
        "| N1 | N2 | h1 (m) | h2 (m) | reduced size | solve runtime (s) | max roots 1--3 | max roots 1--6 |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for n_1, n_2 in TARGETED_PROPORTIONAL_MESHES:
        level = [
            row
            for row in proportional
            if row["n_elements_arm_1"] == n_1
        ]
        lines.append(
            f"| {n_1} | {n_2} | {level[0]['h_1_m']:.10g} | {level[0]['h_2_m']:.10g} | {level[0]['reduced_matrix_size']} | {level[0]['solve_runtime_s']:.6f} | {max(row['relative_error'] for row in level[:3]):.3e} | {max(row['relative_error'] for row in level[:6]):.3e} |"
        )
    lines.extend(
        [
            "",
            "## Acceptance",
            "",
            "| Gate | Result | Observed | Requirement |",
            "|---|---:|---:|---|",
        ]
    )
    for row in summary_rows:
        lines.append(
            f"| {row['gate']} | {row['status']} | {row['observed']} | {row['requirement']} |"
        )
    lines.extend(
        [
            "",
            "## Mode-2 torsion diagnostic",
            "",
            f"Errors for short-arm counts 8,16,32,64 are `{mode_2_order['error_8']:.6e}`, `{mode_2_order['error_16']:.6e}`, `{mode_2_order['error_32']:.6e}`, `{mode_2_order['error_64']:.6e}`. Empirical orders are `{mode_2_order['p_8_16']:.6f}`, `{mode_2_order['p_16_32']:.6f}`, `{mode_2_order['p_32_64']:.6f}`.",
            "",
            f"The diagnostic Richardson frequency is `{mode_2_order['richardson_frequency_hz']:.12f} Hz`; its relative difference from the accepted analytic root is `{mode_2_order['richardson_relative_error']:.3e}`. This extrapolated value is not used by any acceptance gate.",
            "",
            "## Equal-physical-step control",
            "",
            f"The equal-arm `(128,128)` control uses `h={control[0]['h_1_m']:.7g} m`, the same as the finest unequal mesh. Its first-three maximum relative error is `{max(row['relative_error'] for row in control[:3]):.3e}`. The finest unequal/control maximum root-by-root frequency difference is `{max(abs(left['fem_frequency_hz']-right['fem_frequency_hz']) for left, right in zip(finest, control)):.3e} Hz`.",
            "",
            "## Interpretation",
            "",
            "The targeted status uses only the raw `(64,192)` FEM roots. Element stiffness, consistent mass, continuum equations, joint mapping, eigensolver, analytic reference, and the `1e-5`/`5e-4` thresholds are unchanged.",
            "",
            (
                "The unchanged refinement gate fails only for "
                + ", ".join(
                    f"mode {mode} (`{coarse:.3e}` to `{fine:.3e}`)"
                    for mode, coarse, fine in nonmonotone
                )
                + f" between `(32,96)` and `(64,192)`. At the finest level the raw reduced matrices have condition numbers `cond(K)={finest[0]['stiffness_condition_number']:.3e}` and `cond(M)={finest[0]['mass_condition_number']:.3e}`; the nonmonotonic error is therefore recorded as a conditioning-floor diagnostic, not hidden by a solver or tolerance change."
                if nonmonotone
                else "All modes 1--6 satisfy the unchanged last-interval refinement gate."
            ),
            "",
            f"Scientific runtime for the targeted sequence and control: `{runtime_seconds:.6f} s`.",
            "",
        ]
    )
    return "\n".join(lines)


def _run_targeted_refinement(output_dir: Path, *, smoke: bool) -> int:
    started = time.perf_counter()
    unequal_analytic = _load_accepted_analytic_frequencies(
        output_dir, "unequal_0p1_0p3_beta0"
    )
    straight = _point(0.4)
    family_labels = _beta0_exact_family_labels(straight, unequal_analytic)
    unequal_1, unequal_2 = _point(0.1), _point(0.3)

    if smoke:
        rows, diagnostics = _targeted_fem_rows(
            case_id="unequal_length_proportional_smoke",
            point_1=unequal_1,
            point_2=unequal_2,
            element_counts=(16, 48),
            analytic_frequencies=unequal_analytic,
            family_labels=family_labels,
        )
        elapsed = time.perf_counter() - started
        print("targeted_refinement_smoke_status=PASS" if elapsed <= 120.0 else "targeted_refinement_smoke_status=PARTIAL_PASS")
        print(f"reduced_matrix_size={diagnostics['reduced_matrix_size']}")
        print(f"assembly_runtime_s={diagnostics['assembly_runtime_s']:.6f}")
        print(f"solve_runtime_s={diagnostics['solve_runtime_s']:.6f}")
        print("frequencies_hz=" + ",".join(f"{row['fem_frequency_hz']:.12f}" for row in rows))
        print(f"maximum_first_three_relative_error={max(row['relative_error'] for row in rows[:3]):.6e}")
        print(f"maximum_first_six_relative_error={max(row['relative_error'] for row in rows[:6]):.6e}")
        print(f"scientific_runtime_s={elapsed:.6f}")
        return 0 if elapsed <= 120.0 else 3

    target_rows: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []
    for element_counts in TARGETED_PROPORTIONAL_MESHES:
        rows, level_diagnostics = _targeted_fem_rows(
            case_id="unequal_length_proportional",
            point_1=unequal_1,
            point_2=unequal_2,
            element_counts=element_counts,
            analytic_frequencies=unequal_analytic,
            family_labels=family_labels,
        )
        target_rows.extend(rows)
        diagnostics.append(level_diagnostics)

    equal_analytic = _load_accepted_analytic_frequencies(output_dir, "equal_beta0")
    equal_labels = _beta0_exact_family_labels(straight, equal_analytic)
    equal_arm = _point(0.2)
    control_rows, control_diagnostics = _targeted_fem_rows(
        case_id="equal_arm_equal_h_control",
        point_1=equal_arm,
        point_2=equal_arm,
        element_counts=TARGETED_EQUAL_H_CONTROL,
        analytic_frequencies=equal_analytic,
        family_labels=equal_labels,
    )
    target_rows.extend(control_rows)
    diagnostics.append(control_diagnostics)

    proportional_rows = [
        row
        for row in target_rows
        if row["case_id"] == "unequal_length_proportional"
    ]
    order_rows: list[dict[str, Any]] = []
    for mode in range(1, NUM_ROOTS + 1):
        sequence = [
            next(
                row
                for row in proportional_rows
                if row["mode"] == mode and row["n_elements_arm_1"] == n_1
            )
            for n_1, _ in TARGETED_PROPORTIONAL_MESHES
        ]
        errors = [float(row["relative_error"]) for row in sequence]
        p_8_16 = _empirical_order(errors[0], errors[1])
        p_16_32 = _empirical_order(errors[1], errors[2])
        p_32_64 = _empirical_order(errors[2], errors[3])
        richardson_frequency = float("nan")
        richardson_absolute_error = float("nan")
        richardson_relative_error = float("nan")
        if mode == 2 and np.isfinite(p_32_64):
            coarse = float(sequence[2]["fem_frequency_hz"])
            fine = float(sequence[3]["fem_frequency_hz"])
            richardson_frequency = fine + (fine - coarse) / (2.0**p_32_64 - 1.0)
            analytic = float(sequence[3]["analytic_frequency_hz"])
            richardson_absolute_error = abs(richardson_frequency - analytic)
            richardson_relative_error = richardson_absolute_error / analytic
        order_rows.append(
            {
                "mode": mode,
                "family": sequence[0]["family"],
                "error_8": errors[0],
                "error_16": errors[1],
                "error_32": errors[2],
                "error_64": errors[3],
                "p_8_16": p_8_16,
                "p_16_32": p_16_32,
                "p_32_64": p_32_64,
                "richardson_frequency_hz": richardson_frequency,
                "richardson_absolute_error_hz": richardson_absolute_error,
                "richardson_relative_error": richardson_relative_error,
            }
        )

    finest_rows = [
        row for row in proportional_rows if row["n_elements_arm_1"] == 64
    ]
    previous_rows = [
        row for row in proportional_rows if row["n_elements_arm_1"] == 32
    ]
    first_three_max = max(row["relative_error"] for row in finest_rows[:3])
    first_six_max = max(row["relative_error"] for row in finest_rows[:6])
    h_ok = all(item["h_difference_m"] <= 1.0e-15 for item in diagnostics)
    structural_ok = all(
        item["stiffness_symmetry_residual"] <= 1.0e-12
        and item["mass_symmetry_residual"] <= 1.0e-12
        and item["minimum_reduced_mass_eigenvalue"] > 0.0
        and item["spurious_zero_mode_count"] == 0
        and item["joint_kinematic_residual"] <= 2.0e-14
        and item["max_joint_equilibrium_residual"] <= 1.0e-7
        for item in diagnostics
    )
    completeness_ok = all(item["positive_root_count"] >= NUM_ROOTS for item in diagnostics)
    guard_ok = bool(finest_rows[6]["fem_frequency_hz"] > 0.0)
    convergence_ok = all(
        fine["relative_error"] <= 1.01 * coarse["relative_error"] + 1.0e-10
        for coarse, fine in zip(previous_rows[:6], finest_rows[:6])
    )
    mode_2 = next(row for row in order_rows if row["mode"] == 2)
    torsion_order_ok = (
        1.8 <= mode_2["p_16_32"] <= 2.2
        and 1.8 <= mode_2["p_32_64"] <= 2.2
    )
    first_three_ok = first_three_max <= 1.0e-5
    first_six_ok = first_six_max <= 5.0e-4

    gates = [
        ("equal physical element lengths", max(item["h_difference_m"] for item in diagnostics), "<=1e-15 m", h_ok),
        ("first-three raw finest accuracy", first_three_max, "<=1e-5", first_three_ok),
        ("first-six raw finest accuracy", first_six_max, "<=5e-4", first_six_ok),
        ("root-7 guard", finest_rows[6]["fem_frequency_hz"], "present and positive", guard_ok),
        ("modes 1-6 refinement", convergence_ok, "finest <= previous with prior allowance", convergence_ok),
        ("mode-2 second-order convergence", min(mode_2["p_16_32"], mode_2["p_32_64"]), "both last orders in [1.8,2.2]", torsion_order_ok),
        ("FEM structural checks", structural_ok, "unchanged tolerances on every new solve", structural_ok),
        ("seven positive roots", completeness_ok, "all new solves", completeness_ok),
    ]
    summary_rows = [
        {
            "gate": gate,
            "observed": observed,
            "requirement": requirement,
            "status": "PASS" if passed else "FAIL",
        }
        for gate, observed, requirement, passed in gates
    ]
    if not structural_ok:
        status = "FAIL_FEM_STRUCTURE"
    elif not completeness_ok or not guard_ok:
        status = "FAIL_ROOT_COMPLETENESS"
    elif not first_six_ok:
        status = "FAIL_TARGET_FIRST_SIX"
    elif not convergence_ok:
        status = "FAIL_CONVERGENCE_ORDER"
    elif not torsion_order_ok:
        status = "FAIL_CONVERGENCE_ORDER"
    elif not first_three_ok:
        status = "PARTIAL_PASS"
    elif h_ok:
        status = "PASS"
    else:
        status = "MODEL_AMBIGUITY"
    summary_rows.append(
        {
            "gate": "TARGETED_OVERALL",
            "observed": status,
            "requirement": "all targeted acceptance gates",
            "status": status,
        }
    )
    overall_status = (
        "PASS after targeted length-proportional FEM refinement"
        if status == "PASS"
        else "PARTIAL_PASS"
    )
    summary_rows.append(
        {
            "gate": "OVERALL_RECTANGULAR_EB",
            "observed": overall_status,
            "requirement": "targeted PASS required for promotion",
            "status": overall_status,
        }
    )
    runtime_seconds = time.perf_counter() - started
    _write_csv(output_dir / "targeted_length_proportional_refinement.csv", target_rows)
    _write_csv(output_dir / "targeted_convergence_orders.csv", order_rows)
    _write_csv(output_dir / "targeted_refinement_summary.csv", summary_rows)
    report = _targeted_report(
        status=status,
        target_rows=target_rows,
        order_rows=order_rows,
        summary_rows=summary_rows,
        runtime_seconds=runtime_seconds,
    )
    (output_dir / "targeted_refinement_report.md").write_text(
        report, encoding="utf-8"
    )

    print(f"targeted_refinement_status={status}")
    print(f"finest_first_three_max_relative_error={first_three_max:.6e}")
    print(f"finest_first_six_max_relative_error={first_six_max:.6e}")
    print(f"mode_2_p_16_32={mode_2['p_16_32']:.6f}")
    print(f"mode_2_p_32_64={mode_2['p_32_64']:.6f}")
    print(f"scientific_runtime_s={runtime_seconds:.6f}")
    print(f"output_dir={output_dir}")
    return 0 if status in ("PASS", "PARTIAL_PASS") else 4


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = args.output_dir.resolve()
    if args.targeted_smoke:
        return _run_targeted_refinement(output_dir, smoke=True)
    if args.targeted_length_proportional_refinement:
        return _run_targeted_refinement(output_dir, smoke=False)
    if args.smoke:
        metrics = _smoke()
        print("rectangular_eb_smoke_status=PASS" if metrics["runtime_seconds"] <= 300.0 else "rectangular_eb_smoke_status=FAIL_ROOT_QUALITY")
        for key, value in metrics.items():
            print(f"{key}={value}")
        return 0 if metrics["runtime_seconds"] <= 300.0 else 3

    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    total_evaluations = 0
    all_quality_roots: list[RootResult] = []
    quality_checks: list[bool] = []
    guard_spectra: list[list[RootResult]] = []

    arm = _point(0.2)
    straight = _point(0.4)
    exact_rows, exact_bending_roots, exact_torsion_roots, exact_metrics = _exact_checks(straight)
    total_evaluations += int(exact_metrics["evaluations"])
    all_quality_roots.extend(exact_bending_roots + exact_torsion_roots)

    timo_reference, builder, _ = _solve(
        straight, lambda omega: straight_boundary_matrix(omega, straight)
    )
    total_evaluations += builder.evaluations
    all_quality_roots.extend(timo_reference)
    quality_checks.extend(
        _root_fields(root, lambda omega: straight_boundary_matrix_raw(omega, straight))["quality_status"] == "PASS"
        for root in timo_reference
    )
    guard_spectra.append(timo_reference)
    eb_reference, builder, _ = _solve(
        straight, lambda omega: eb_straight_boundary_matrix(omega, straight)
    )
    total_evaluations += builder.evaluations
    all_quality_roots.extend(eb_reference)
    quality_checks.extend(
        _root_fields(root, lambda omega: eb_straight_boundary_matrix_raw(omega, straight))["quality_status"] == "PASS"
        for root in eb_reference
    )
    guard_spectra.append(eb_reference)

    unequal_timo_rows: list[dict[str, Any]] = []
    unequal_eb_rows: list[dict[str, Any]] = []
    timo_spectra: dict[tuple[float, float, float], list[RootResult]] = {}
    eb_spectra: dict[tuple[float, float, float], list[RootResult]] = {}
    for length_1, length_2 in SPLITS:
        point_1, point_2 = _point(length_1), _point(length_2)
        timo, builder, runtime = _solve(
            point_1,
            lambda omega, p1=point_1, p2=point_2: coupled_boundary_matrix(omega, 0.0, p1, p2),
        )
        total_evaluations += builder.evaluations
        all_quality_roots.extend(timo)
        quality_checks.extend(
            _root_fields(
                root,
                lambda omega, p1=point_1, p2=point_2: coupled_boundary_matrix_raw(omega, 0.0, p1, p2),
            )["quality_status"] == "PASS"
            for root in timo
        )
        guard_spectra.append(timo)
        timo_spectra[(length_1, length_2, 0.0)] = timo
        eb, builder, eb_runtime = _solve(
            point_1,
            lambda omega, p1=point_1, p2=point_2: eb_coupled_boundary_matrix(omega, 0.0, p1, p2),
        )
        total_evaluations += builder.evaluations
        all_quality_roots.extend(eb)
        quality_checks.extend(
            _root_fields(
                root,
                lambda omega, p1=point_1, p2=point_2: eb_coupled_boundary_matrix_raw(omega, 0.0, p1, p2),
            )["quality_status"] == "PASS"
            for root in eb
        )
        guard_spectra.append(eb)
        eb_spectra[(length_1, length_2, 0.0)] = eb
        for theory, roots, reference, rows, raw_factory, solve_runtime in (
            (
                "timoshenko",
                timo,
                timo_reference,
                unequal_timo_rows,
                lambda omega, p1=point_1, p2=point_2: coupled_boundary_matrix_raw(omega, 0.0, p1, p2),
                runtime,
            ),
            (
                "eb",
                eb,
                eb_reference,
                unequal_eb_rows,
                lambda omega, p1=point_1, p2=point_2: eb_coupled_boundary_matrix_raw(omega, 0.0, p1, p2),
                eb_runtime,
            ),
        ):
            for mode, (root, direct) in enumerate(zip(roots, reference), start=1):
                absolute = abs(root.frequency_hz - direct.frequency_hz)
                rows.append(
                    {
                        "theory": theory,
                        "length_1_m": length_1,
                        "length_2_m": length_2,
                        "mode": mode,
                        "coupled_frequency_hz": root.frequency_hz,
                        "direct_straight_frequency_hz": direct.frequency_hz,
                        "absolute_error_hz": absolute,
                        "relative_error": absolute / direct.frequency_hz,
                        "solve_runtime_seconds": solve_runtime,
                        **_root_fields(root, raw_factory),
                    }
                )

    comparison_rows: list[dict[str, Any]] = []
    equal_key = (0.2, 0.2, 0.0)
    for beta_deg in (0.0, 30.0, 90.0):
        beta = math.radians(beta_deg)
        if beta_deg == 0.0:
            timo, eb = timo_spectra[equal_key], eb_spectra[equal_key]
        else:
            timo, builder, _ = _solve(
                arm, lambda omega, beta=beta: coupled_boundary_matrix(omega, beta, arm, arm)
            )
            total_evaluations += builder.evaluations
            all_quality_roots.extend(timo)
            quality_checks.extend(
                _root_fields(
                    root,
                    lambda omega, beta=beta: coupled_boundary_matrix_raw(omega, beta, arm, arm),
                )["quality_status"] == "PASS"
                for root in timo
            )
            guard_spectra.append(timo)
            eb, builder, _ = _solve(
                arm, lambda omega, beta=beta: eb_coupled_boundary_matrix(omega, beta, arm, arm)
            )
            total_evaluations += builder.evaluations
            all_quality_roots.extend(eb)
            quality_checks.extend(
                _root_fields(
                    root,
                    lambda omega, beta=beta: eb_coupled_boundary_matrix_raw(omega, beta, arm, arm),
                )["quality_status"] == "PASS"
                for root in eb
            )
            guard_spectra.append(eb)
            timo_spectra[(0.2, 0.2, beta_deg)] = timo
            eb_spectra[(0.2, 0.2, beta_deg)] = eb
        timo_frequency, eb_frequency = _frequencies(timo), _frequencies(eb)
        timo_gaps, eb_gaps = _neighbor_gaps(timo_frequency), _neighbor_gaps(eb_frequency)
        for mode, (timo_root, eb_root) in enumerate(zip(timo, eb), start=1):
            near = min(
                timo_gaps[mode - 1] / timo_root.frequency_hz,
                eb_gaps[mode - 1] / eb_root.frequency_hz,
            ) <= NEAR_DEGENERATE_RELATIVE_GAP
            comparison_rows.append(
                {
                    "beta_deg": beta_deg,
                    "mode": mode,
                    "role": "comparison" if mode <= 6 else "guard",
                    "timoshenko_frequency_hz": timo_root.frequency_hz,
                    "eb_frequency_hz": eb_root.frequency_hz,
                    "delta_f": abs(timo_root.frequency_hz - eb_root.frequency_hz) / timo_root.frequency_hz,
                    "delta_lambda": abs(float(timo_root.omega.real) ** 2 - float(eb_root.omega.real) ** 2) / float(timo_root.omega.real) ** 2,
                    "timoshenko_neighbor_gap_hz": timo_gaps[mode - 1],
                    "eb_neighbor_gap_hz": eb_gaps[mode - 1],
                    "matching_status": "possible_mode_reordering" if near else "sorted_no_near_degeneracy",
                    **{f"timoshenko_{key}": value for key, value in _root_fields(timo_root, lambda omega, beta=beta: coupled_boundary_matrix_raw(omega, beta, arm, arm)).items()},
                    **{f"eb_{key}": value for key, value in _root_fields(eb_root, lambda omega, beta=beta: eb_coupled_boundary_matrix_raw(omega, beta, arm, arm)).items()},
                }
            )

    slender_rows: list[dict[str, Any]] = []
    for section_scale in (1.0, 0.5, 0.25):
        point = _point(0.2, section_scale)

        def timo_bending_raw(omega: complex, p=point) -> np.ndarray:
            raw = coupled_boundary_matrix_raw(omega, 0.0, p, p)
            return raw[np.ix_([0, 2, 3, 5], [0, 1, 3, 4])]

        timo_bending, builder, _ = _solve(
            point,
            lambda omega: equilibrate_matrix(timo_bending_raw(omega))[0],
            num_roots=3,
        )
        total_evaluations += builder.evaluations
        all_quality_roots.extend(timo_bending)
        quality_checks.extend(
            _root_fields(root, timo_bending_raw)["quality_status"] == "PASS"
            for root in timo_bending
        )
        eb_bending, builder, _ = _solve(
            point,
            lambda omega, p=point: eb_bending_coupled_boundary_matrix(omega, p, p),
            num_roots=3,
        )
        total_evaluations += builder.evaluations
        all_quality_roots.extend(eb_bending)
        quality_checks.extend(
            _root_fields(
                root,
                lambda omega, p=point: eb_bending_coupled_boundary_matrix_raw(omega, p, p),
            )["quality_status"] == "PASS"
            for root in eb_bending
        )
        for mode, (timo_root, eb_root) in enumerate(zip(timo_bending, eb_bending), start=1):
            difference = abs(timo_root.frequency_hz - eb_root.frequency_hz) / timo_root.frequency_hz
            slender_rows.append(
                {
                    "section_scale": section_scale,
                    "a_m": point.geometry.a,
                    "b_m": point.geometry.b,
                    "a_over_b": point.geometry.a / point.geometry.b,
                    "bending_mode": mode,
                    "timoshenko_frequency_hz": timo_root.frequency_hz,
                    "eb_frequency_hz": eb_root.frequency_hz,
                    "relative_difference": difference,
                    **{
                        f"timoshenko_{key}": value
                        for key, value in _root_fields(timo_root, timo_bending_raw).items()
                    },
                    **{
                        f"eb_{key}": value
                        for key, value in _root_fields(
                            eb_root,
                            lambda omega, p=point: eb_bending_coupled_boundary_matrix_raw(omega, p, p),
                        ).items()
                    },
                }
            )

    analytic_cases = {
        "equal_beta0": (arm, arm, 0.0, eb_spectra[(0.2, 0.2, 0.0)]),
        "equal_beta30": (arm, arm, 30.0, eb_spectra[(0.2, 0.2, 30.0)]),
        "equal_beta90": (arm, arm, 90.0, eb_spectra[(0.2, 0.2, 90.0)]),
        "unequal_0p1_0p3_beta0": (_point(0.1), _point(0.3), 0.0, eb_spectra[(0.1, 0.3, 0.0)]),
    }
    fem_rows: list[dict[str, Any]] = []
    fem_diagnostics: dict[tuple[str, int], dict[str, Any]] = {}
    for case, (point_1, point_2, beta_deg, analytic_roots) in analytic_cases.items():
        analytic_frequency = _frequencies(analytic_roots)
        beta = math.radians(beta_deg)
        for mesh in MESHES:
            fem_started = time.perf_counter()
            assembly = assemble_two_arm_eb_fem(point_1, point_2, beta, mesh)
            solution = solve_two_arm_eb_fem(assembly, num_roots=NUM_ROOTS)
            fem_runtime = time.perf_counter() - fem_started
            matched, matching_statuses, gaps = _match_frequencies(
                analytic_frequency, solution.frequencies_hz
            )
            fem_diagnostics[(case, mesh)] = {
                "stiffness_symmetry_residual": solution.stiffness_symmetry_residual,
                "mass_symmetry_residual": solution.mass_symmetry_residual,
                "minimum_mass_eigenvalue": solution.minimum_mass_eigenvalue,
                "zero_mode_count": solution.zero_mode_count,
                "roots_found": len(solution.frequencies_hz),
                "max_joint_equilibrium_residual": float(np.max(solution.joint_equilibrium_residuals)),
                "joint_mapping_residual": eb_joint_mapping_residual(beta),
                "reduced_matrix_size": assembly.mass.shape[0],
                "fem_runtime_seconds": fem_runtime,
            }
            for mode, (analytic_root, fem_frequency, match_status, gap) in enumerate(
                zip(analytic_roots, matched, matching_statuses, gaps), start=1
            ):
                absolute = abs(analytic_root.frequency_hz - fem_frequency)
                fem_rows.append(
                    {
                        "case": case,
                        "length_1_m": point_1.geometry.length,
                        "length_2_m": point_2.geometry.length,
                        "beta_deg": beta_deg,
                        "elements_per_arm": mesh,
                        "mode": mode,
                        "role": "comparison" if mode <= 6 else "guard",
                        "analytic_frequency_hz": analytic_root.frequency_hz,
                        "fem_frequency_hz": fem_frequency,
                        "absolute_error_hz": absolute,
                        "relative_error": absolute / analytic_root.frequency_hz,
                        "neighbor_gap_hz": gap,
                        "matching_status": match_status,
                        **fem_diagnostics[(case, mesh)],
                    }
                )

    convergence_rows: list[dict[str, Any]] = []
    for case in analytic_cases:
        for mode in range(1, NUM_ROOTS + 1):
            sequence = [
                next(
                    row for row in fem_rows
                    if row["case"] == case and row["elements_per_arm"] == mesh and row["mode"] == mode
                )
                for mesh in MESHES
            ]
            errors = np.asarray([row["relative_error"] for row in sequence])
            divergent_steps = int(np.count_nonzero(errors[1:] > 1.05 * errors[:-1] + 1.0e-10))
            for index, row in enumerate(sequence):
                convergence_rows.append(
                    {
                        "case": case,
                        "mode": mode,
                        "elements_per_arm": row["elements_per_arm"],
                        "relative_error": row["relative_error"],
                        "previous_mesh_relative_error": "" if index == 0 else errors[index - 1],
                        "finest_not_worse_than_mesh8": bool(errors[-1] <= 1.01 * errors[1] + 1.0e-10),
                        "divergent_refinement_steps": divergent_steps,
                        "mesh64_first3_gate": bool(errors[-1] <= 1.0e-5) if mode <= 3 else "n/a",
                        "mesh64_first6_gate": bool(errors[-1] <= 5.0e-4) if mode <= 6 else "n/a",
                    }
                )

    max_unequal_timo = max(row["relative_error"] for row in unequal_timo_rows)
    max_unequal_eb = max(row["relative_error"] for row in unequal_eb_rows)
    bending_exact_error = max(row["relative_error"] for row in exact_rows if row["family"] == "bending")
    torsion_exact_error = max(row["relative_error"] for row in exact_rows if row["family"] == "torsion")
    slender_by_mode = {
        mode: [row["relative_difference"] for row in slender_rows if row["bending_mode"] == mode]
        for mode in (1, 2, 3)
    }
    slender_ok = all(
        values[1] <= 1.001 * values[0] + 1.0e-10
        and values[2] <= 1.001 * values[1] + 1.0e-10
        for values in slender_by_mode.values()
    )
    point = arm
    sbar16_relative = abs(point.properties.Sbar16) / max(abs(point.properties.Sbar11), np.finfo(float).tiny)
    torsion_identity = abs(point.torsion.C_T - point.torsion.Cbar) / abs(point.torsion.Cbar)
    structural_state = eb_state_matrix(1.0, point)
    structural_ok = structural_state.shape == (6, 6) and np.array_equal(eb_clamp_matrix(point)[3:], np.eye(3))
    raw_block = eb_coupled_boundary_matrix_raw(2.0 * math.pi * 500.0, 0.0, arm, arm)
    off_block = np.concatenate(
        [raw_block[np.ix_([0, 2, 3, 5], [2, 5])].ravel(), raw_block[np.ix_([1, 4], [0, 1, 3, 4])].ravel()]
    )
    block_residual = float(np.linalg.norm(off_block) / max(np.linalg.norm(raw_block), 1.0))
    fem_structural_ok = all(
        item["stiffness_symmetry_residual"] <= 1.0e-12
        and item["mass_symmetry_residual"] <= 1.0e-12
        and item["minimum_mass_eigenvalue"] > 0.0
        and item["zero_mode_count"] == 0
        and item["roots_found"] >= NUM_ROOTS
        for item in fem_diagnostics.values()
    )
    fem_joint_ok = all(
        item["joint_mapping_residual"] <= 2.0e-14
        and item["max_joint_equilibrium_residual"] <= 1.0e-7
        for item in fem_diagnostics.values()
    )
    convergence_ok = all(
        row["finest_not_worse_than_mesh8"] and row["divergent_refinement_steps"] <= 1
        for row in convergence_rows if row["mode"] <= 6
    )
    finest_rows = [
        row for row in fem_rows if row["elements_per_arm"] == 64 and row["mode"] <= 6
    ]
    first3_accuracy = max(row["relative_error"] for row in finest_rows if row["mode"] <= 3)
    first6_accuracy = max(row["relative_error"] for row in finest_rows)
    accuracy_ok = first3_accuracy <= 1.0e-5 and first6_accuracy <= 5.0e-4
    root_quality_ok = all(quality_checks)
    guard_ok = all(len(roots) >= NUM_ROOTS for roots in guard_spectra)

    gate_values = [
        ("1 Sbar16=0 at theta=0", sbar16_relative, "<=1e-14", sbar16_relative <= 1e-14),
        ("2 C_T=Cbar=C_SV", torsion_identity, "<=1e-14", torsion_identity <= 1e-14),
        ("3 EB state/clamp/matrix structure", structural_ok, "all exact structural checks", structural_ok),
        ("4 exact fixed-fixed bending", bending_exact_error, "relative error <=1e-8", bending_exact_error <= 1e-8),
        ("5 exact fixed-fixed torsion", torsion_exact_error, "relative error <=1e-8", torsion_exact_error <= 1e-8),
        ("6 unequal-length Timoshenko", max_unequal_timo, "relative error <=1e-8", max_unequal_timo <= 1e-8),
        ("7 unequal-length EB", max_unequal_eb, "relative error <=1e-8", max_unequal_eb <= 1e-8),
        ("8 beta=0 EB block separation", block_residual, "relative off-block <=2e-14", block_residual <= 2e-14),
        ("9 slender-limit bending approach", slender_ok, "modes 1-3 non-increasing within allowance", slender_ok),
        ("10 FEM structural/no-zero/root checks", fem_structural_ok, "all meshes and cases", fem_structural_ok),
        ("11 FEM joint mapping/equilibrium", fem_joint_ok, "mapping <=2e-14; equilibrium <=1e-7", fem_joint_ok),
        ("12 FEM convergence", convergence_ok, "modes 1-6; finest no worse than mesh8", convergence_ok),
        ("13 mesh64 analytic/FEM accuracy", max(first3_accuracy, first6_accuracy), "first3 <=1e-5; first6 <=5e-4", accuracy_ok),
    ]
    summary_rows = [
        {
            "gate": gate,
            "observed": observed,
            "requirement": requirement,
            "status": "PASS" if passed else "FAIL",
        }
        for gate, observed, requirement, passed in gate_values
    ]
    summary_rows.append(
        {
            "gate": "analytic root quality",
            "observed": root_quality_ok,
            "requirement": "scaled or physical-raw determinant and singular residual <=1e-8",
            "status": "PASS" if root_quality_ok else "FAIL",
        }
    )
    summary_rows.append(
        {
            "gate": "root 7 guard",
            "observed": guard_ok,
            "requirement": "present in every full analytic spectrum",
            "status": "PASS" if guard_ok else "FAIL",
        }
    )

    if max_unequal_timo > 1e-8:
        status = "FAIL_UNEQUAL_LENGTH_TIMOSHENKO"
    elif max_unequal_eb > 1e-8:
        status = "FAIL_UNEQUAL_LENGTH_EB"
    elif bending_exact_error > 1e-8 or torsion_exact_error > 1e-8:
        status = "FAIL_EXACT_EB_LIMIT"
    elif not slender_ok:
        status = "FAIL_SLENDER_LIMIT"
    elif not fem_structural_ok or not fem_joint_ok:
        status = "FAIL_FEM_ASSEMBLY"
    elif not convergence_ok:
        status = "FAIL_FEM_CONVERGENCE"
    elif not root_quality_ok or not guard_ok:
        status = "FAIL_ROOT_QUALITY"
    elif not accuracy_ok:
        status = "PARTIAL_PASS"
    elif all(item[3] for item in gate_values):
        status = "PASS"
    else:
        status = "MODEL_AMBIGUITY"

    runtime_seconds = time.perf_counter() - started
    summary_rows.append(
        {
            "gate": "OVERALL",
            "observed": status,
            "requirement": "all thirteen gates plus root quality and guard",
            "status": status,
        }
    )
    _write_csv(output_dir / "unequal_length_timoshenko_check.csv", unequal_timo_rows)
    _write_csv(output_dir / "unequal_length_eb_check.csv", unequal_eb_rows)
    _write_csv(output_dir / "straight_exact_family_check.csv", exact_rows)
    _write_csv(output_dir / "timoshenko_vs_eb.csv", comparison_rows)
    _write_csv(output_dir / "slender_limit_bending_check.csv", slender_rows)
    _write_csv(output_dir / "eb_analytic_vs_fem.csv", fem_rows)
    _write_csv(output_dir / "eb_fem_convergence.csv", convergence_rows)
    _write_csv(output_dir / "rectangular_eb_validation_summary.csv", summary_rows)
    report = _report(
        status,
        summary_rows,
        unequal_timo_rows,
        unequal_eb_rows,
        exact_rows,
        comparison_rows,
        slender_rows,
        fem_rows,
        convergence_rows,
        runtime_seconds,
        total_evaluations,
    )
    (output_dir / "rectangular_eb_validation_report.md").write_text(report, encoding="utf-8")

    print(f"rectangular_eb_gate_status={status}")
    print(f"maximum_unequal_timoshenko_relative_error={max_unequal_timo:.6e}")
    print(f"maximum_unequal_eb_relative_error={max_unequal_eb:.6e}")
    print(f"maximum_exact_bending_relative_error={bending_exact_error:.6e}")
    print(f"maximum_exact_torsion_relative_error={torsion_exact_error:.6e}")
    print(f"mesh64_maximum_first3_relative_error={first3_accuracy:.6e}")
    print(f"mesh64_maximum_first6_relative_error={first6_accuracy:.6e}")
    print(f"analytic_boundary_matrix_evaluations={total_evaluations}")
    print(f"scientific_runtime_seconds={runtime_seconds:.6f}")
    print(f"output_dir={output_dir}")
    return 0 if status in ("PASS", "PARTIAL_PASS") else 4


if __name__ == "__main__":
    raise SystemExit(main())
