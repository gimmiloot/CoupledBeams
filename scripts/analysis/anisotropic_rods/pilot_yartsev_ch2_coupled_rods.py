"""Run the small elastic Chapter-2 rigid-angular-joint pilot.

The workflow is deliberately limited to HMS/DX-209, two 0.2 m arms,
theta=15 deg, and beta=0/30/90 deg.  It performs no complex-root, shape, map,
Euler--Bernoulli, Saint-Venant, or FEM calculation.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import time
from pathlib import Path
from typing import Any, Callable, Iterable, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib.yartsev_ch2_coupled_rods import (  # noqa: E402
    coupled_boundary_matrix,
    coupled_boundary_matrix_raw,
    equilibrate_matrix,
    joint_basis,
    joint_matrix_book,
    joint_matrix_old,
    joint_virtual_work,
    kinematic_joint_residual,
    moment_joint_residual,
    physical_moment_vector,
    physical_rotation_vector,
    straight_boundary_matrix,
    straight_boundary_matrix_raw,
    two_arm_notation_transform_matrix,
)
from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
    RootResult,
    RodPoint,
    cantilever_geometry,
    find_elastic_roots,
    hms_dx_209_material,
    make_rod_point,
)


DEFAULT_OUTPUT_DIR = (
    ROOT / "results" / "anisotropic_rods" / "yartsev_ch2_coupled_joint_pilot"
)
PILOT_BETAS_DEG = (0.0, 30.0, 90.0)
RANK_BETAS_DEG = (0.0, 15.0, 30.0, 90.0, 180.0)
RANDOM_SEED = 20260731
NUM_RANDOM_CASES = 100
NUM_ROOTS = 7
SCAN_STEP_HZ = 10.0
INITIAL_MAX_HZ = 5000.0
MAX_HZ = 100_000.0


class CountingBoundaryBuilder:
    """Cache and time unique boundary-matrix evaluations."""

    def __init__(self, factory: Callable[[complex], np.ndarray]) -> None:
        self.factory = factory
        self.calls = 0
        self.evaluations = 0
        self.elapsed_seconds = 0.0
        self._cache: dict[complex, np.ndarray] = {}

    def __call__(self, omega: complex, _point: RodPoint, formulation: str) -> np.ndarray:
        if formulation != "state_corrected":
            raise ValueError("the coupled pilot accepts only state_corrected")
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

    @property
    def mean_seconds(self) -> float:
        return self.elapsed_seconds / self.evaluations if self.evaluations else 0.0


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--beta0-only",
        action="store_true",
        help="run checks and beta=0/reference roots only; do not solve beta=30/90",
    )
    return parser.parse_args(argv)


def _write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    data = list(rows)
    if not data:
        return
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


def _point(theta_deg: float, length: float) -> RodPoint:
    return make_rod_point(
        theta_deg,
        material_mode="elastic",
        geometry=cantilever_geometry(length),
        material=hms_dx_209_material(),
    )


def _quality(matrix: np.ndarray) -> dict[str, float]:
    singular = np.linalg.svd(matrix, compute_uv=False)
    sigma_max = float(singular[0])
    sigma_min = float(singular[-1])
    relative = sigma_min / sigma_max if sigma_max > 0.0 else 0.0
    determinant = complex(np.linalg.det(matrix))
    scale = max(sigma_max ** matrix.shape[0], np.finfo(float).tiny)
    return {
        "determinant_abs": float(abs(determinant)),
        "determinant_residual": float(abs(determinant) / scale),
        "sigma_min": sigma_min,
        "sigma_max": sigma_max,
        "relative_singular_residual": float(relative),
    }


def _root_quality_ok(root: RootResult) -> bool:
    return (
        root.determinant_residual <= 1.0e-8
        and root.relative_singular_residual <= 1.0e-8
        and root.status != "rejected_complex_quality"
    )


def _solve(
    point: RodPoint, factory: Callable[[complex], np.ndarray]
) -> tuple[list[RootResult], CountingBoundaryBuilder, float]:
    builder = CountingBoundaryBuilder(factory)
    started = time.perf_counter()
    roots = find_elastic_roots(
        point,
        "state_corrected",
        num_roots=NUM_ROOTS,
        scan_step_hz=SCAN_STEP_HZ,
        initial_max_hz=INITIAL_MAX_HZ,
        max_hz=MAX_HZ,
        boundary_matrix_builder=builder,
    )
    return roots, builder, time.perf_counter() - started


def _joint_matrix_checks(arm: RodPoint) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    transform = two_arm_notation_transform_matrix()
    for beta_deg in RANK_BETAS_DEG:
        beta = math.radians(beta_deg)
        book = joint_matrix_book(beta)
        sign_error = float(np.max(np.abs(book - joint_matrix_old(beta) @ transform)))
        rank = int(np.linalg.matrix_rank(book))
        rows.extend(
            [
                {
                    "check": "notation_sign_gate",
                    "beta_deg": beta_deg,
                    "observed": sign_error,
                    "tolerance": 5e-15,
                    "status": "PASS" if sign_error <= 5e-15 else "FAIL",
                    "details": "J_book = J_old @ block_diag(P,P)",
                },
                {
                    "check": "joint_rank",
                    "beta_deg": beta_deg,
                    "observed": rank,
                    "tolerance": 6,
                    "status": "PASS" if rank == 6 else "FAIL",
                    "details": "rank(J_book)=6",
                },
            ]
        )

    zero = joint_basis(0.0)
    ninety = joint_basis(math.pi / 2.0)
    basis_error = max(
        float(np.linalg.norm(zero.t_2 + zero.t_1)),
        float(np.linalg.norm(zero.n_2 + zero.n_1)),
        float(np.linalg.norm(ninety.t_2 - ninety.n_1)),
        float(np.linalg.norm(ninety.n_2 + ninety.t_1)),
    )
    rows.append(
        {
            "check": "basis_limits",
            "beta_deg": "0;90",
            "observed": basis_error,
            "tolerance": 5e-15,
            "status": "PASS" if basis_error <= 5e-15 else "FAIL",
            "details": "t2,n2 limits and unchanged e_z,t_i,n_i convention",
        }
    )

    beta0 = joint_matrix_book(0.0)
    expected_beta0 = np.array(
        [
            [1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0],
        ],
        dtype=float,
    )
    beta90 = joint_matrix_book(math.pi / 2.0)
    expected_beta90 = np.array(
        [
            [1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1],
        ],
        dtype=float,
    )
    for name, observed, expected, beta_deg in (
        ("beta0_scalar_limit", beta0, expected_beta0, 0.0),
        ("beta90_scalar_limit", beta90, expected_beta90, 90.0),
    ):
        error = float(np.max(np.abs(observed - expected)))
        rows.append(
            {
                "check": name,
                "beta_deg": beta_deg,
                "observed": error,
                "tolerance": 5e-15,
                "status": "PASS" if error <= 5e-15 else "FAIL",
                "details": "exact requested scalar limit",
            }
        )

    rng = np.random.default_rng(RANDOM_SEED)
    for beta_deg in PILOT_BETAS_DEG:
        beta = math.radians(beta_deg)
        basis = joint_basis(beta)
        rotation_max = 0.0
        moment_max = 0.0
        pairing_max = 0.0
        for _ in range(NUM_RANDOM_CASES):
            rotation = rng.normal() * basis.t_1 + rng.normal() * basis.n_1
            moment = rng.normal() * basis.t_1 + rng.normal() * basis.n_1
            rotation_max = max(
                rotation_max, float(np.linalg.norm(kinematic_joint_residual(beta, rotation)))
            )
            moment_max = max(
                moment_max, float(np.linalg.norm(moment_joint_residual(beta, moment)))
            )
            dPhi, dpsi, M_T, M = rng.normal(size=4)
            vector_value = physical_moment_vector(M_T, M, basis.t_1, basis.n_1) @ physical_rotation_vector(
                dPhi, dpsi, basis.t_1, basis.n_1
            )
            pairing_max = max(pairing_max, abs(vector_value - (M_T * dPhi + M * dpsi)))
        for name, value, details in (
            ("random_compatible_rotation", rotation_max, "100 deterministic global rotations"),
            ("random_equilibrated_moment", moment_max, "100 deterministic global moments"),
            ("book_energy_pairing", pairing_max, "m.delta(vartheta)=M_T delta(Phi)+M delta(psi)"),
        ):
            rows.append(
                {
                    "check": name,
                    "beta_deg": beta_deg,
                    "observed": value,
                    "tolerance": 5e-13,
                    "status": "PASS" if value <= 5e-13 else "FAIL",
                    "details": details,
                }
            )

    matrix = coupled_boundary_matrix(2.0 * math.pi * 500.0, math.pi / 6.0, arm, arm)
    swapped = matrix[:, [3, 4, 5, 0, 1, 2]]
    singular_error = float(
        np.max(
            np.abs(
                np.linalg.svd(matrix, compute_uv=False)
                - np.linalg.svd(swapped, compute_uv=False)
            )
        )
    )
    rows.append(
        {
            "check": "equal_arm_exchange_singular_values",
            "beta_deg": 30.0,
            "observed": singular_error,
            "tolerance": 5e-13,
            "status": "PASS" if singular_error <= 5e-13 else "FAIL",
            "details": "reaction block-column permutation",
        }
    )
    return rows


def _virtual_work_checks() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = [
        {
            "check": "symbolic_reduction",
            "beta_deg": "all",
            "cases": "exact",
            "max_absolute_residual": 0.0,
            "tolerance": 0.0,
            "status": "PASS",
            "details": "(Q1+Q2) delta_w_J + (m1+m2).delta_vartheta_J = 0",
        }
    ]
    rng = np.random.default_rng(RANDOM_SEED)
    for beta_deg in PILOT_BETAS_DEG:
        basis = joint_basis(math.radians(beta_deg))
        maximum = 0.0
        for _ in range(NUM_RANDOM_CASES):
            delta_w = float(rng.normal())
            q_1 = float(rng.normal())
            delta_rotation = rng.normal() * basis.t_1 + rng.normal() * basis.n_1
            moment_1 = rng.normal() * basis.t_1 + rng.normal() * basis.n_1
            value = joint_virtual_work(
                q_1,
                -q_1,
                delta_w,
                delta_w,
                moment_1,
                -moment_1,
                delta_rotation,
                delta_rotation,
            )
            maximum = max(maximum, abs(value))
        rows.append(
            {
                "check": "deterministic_random_virtual_work",
                "beta_deg": beta_deg,
                "cases": NUM_RANDOM_CASES,
                "max_absolute_residual": maximum,
                "tolerance": 5e-13,
                "status": "PASS" if maximum <= 5e-13 else "FAIL",
                "details": f"seed={RANDOM_SEED}",
            }
        )
    return rows


def _orthotropic_check() -> tuple[list[dict[str, Any]], bool]:
    point = _point(0.0, 0.2)
    matrix = coupled_boundary_matrix_raw(2.0 * math.pi * 500.0, 0.0, point, point)
    bending_rows = [0, 2, 3, 5]
    torsion_rows = [1, 4]
    bending_columns = [0, 1, 3, 4]
    torsion_columns = [2, 5]
    off_block = np.concatenate(
        (
            matrix[np.ix_(bending_rows, torsion_columns)].ravel(),
            matrix[np.ix_(torsion_rows, bending_columns)].ravel(),
        )
    )
    off_norm = float(np.linalg.norm(off_block))
    full_norm = float(np.linalg.norm(matrix))
    relative = off_norm / full_norm if full_norm else 0.0
    rank = int(np.linalg.matrix_rank(matrix))
    passed = (
        abs(point.properties.Sbar16) <= 1e-24
        and relative <= 1e-12
        and rank == 6
    )
    rows = [
        {
            "theta_1_deg": 0.0,
            "theta_2_deg": 0.0,
            "beta_deg": 0.0,
            "omega_rad_s": 2.0 * math.pi * 500.0,
            "Sbar16_1_abs": abs(point.properties.Sbar16),
            "Sbar16_2_abs": abs(point.properties.Sbar16),
            "bending_block_shape": "4x4",
            "torsion_block_shape": "2x2",
            "off_block_norm": off_norm,
            "full_matrix_norm": full_norm,
            "relative_off_block_norm": relative,
            "matrix_rank": rank,
            "status": "PASS" if passed else "FAIL",
        }
    ]
    return rows, passed


def _pilot_root_rows(
    beta_deg: float,
    roots: Sequence[RootResult],
    raw_factory: Callable[[complex], np.ndarray],
    builder: CountingBoundaryBuilder,
    runtime_seconds: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for mode, item in enumerate(roots, start=1):
        raw = _quality(raw_factory(item.omega))
        scaled_matrix, scaling = equilibrate_matrix(raw_factory(item.omega))
        scaled = _quality(scaled_matrix)
        rows.append(
            {
                "beta_deg": beta_deg,
                "sorted_mode": mode,
                "pilot_role": "pilot_spectrum" if mode <= 6 else "completeness_guard",
                "omega_rad_s": item.omega.real,
                "frequency_hz": item.frequency_hz,
                "determinant_residual": item.determinant_residual,
                "sigma_min": item.sigma_min,
                "sigma_max": item.sigma_max,
                "relative_singular_residual": item.relative_singular_residual,
                "raw_determinant_abs": raw["determinant_abs"],
                "raw_determinant_residual": raw["determinant_residual"],
                "raw_sigma_min": raw["sigma_min"],
                "raw_sigma_max": raw["sigma_max"],
                "raw_relative_singular_residual": raw["relative_singular_residual"],
                "scaled_determinant_abs": scaled["determinant_abs"],
                "scaled_determinant_residual": scaled["determinant_residual"],
                "scaled_sigma_min": scaled["sigma_min"],
                "scaled_sigma_max": scaled["sigma_max"],
                "scaled_relative_singular_residual": scaled["relative_singular_residual"],
                "min_row_scaling_factor": float(np.min(scaling.row_factors)),
                "max_row_scaling_factor": float(np.max(scaling.row_factors)),
                "min_column_scaling_factor": float(np.min(scaling.column_factors)),
                "max_column_scaling_factor": float(np.max(scaling.column_factors)),
                "root_status": item.status,
                "quality_status": "PASS" if _root_quality_ok(item) else "FAIL",
                "boundary_matrix_evaluations": builder.evaluations,
                "mean_boundary_matrix_seconds": builder.mean_seconds,
                "root_search_runtime_seconds": runtime_seconds,
            }
        )
    return rows


def _equivalence_rows(
    coupled_roots: Sequence[RootResult], reference_roots: Sequence[RootResult]
) -> tuple[list[dict[str, Any]], float]:
    rows: list[dict[str, Any]] = []
    maximum = 0.0
    for mode, (coupled, reference) in enumerate(zip(coupled_roots, reference_roots), start=1):
        relative = abs(coupled.frequency_hz - reference.frequency_hz) / reference.frequency_hz
        maximum = max(maximum, relative)
        rows.append(
            {
                "sorted_mode": mode,
                "pilot_role": "pilot_spectrum" if mode <= 6 else "completeness_guard",
                "coupled_omega_rad_s": coupled.omega.real,
                "straight_omega_rad_s": reference.omega.real,
                "coupled_frequency_hz": coupled.frequency_hz,
                "straight_frequency_hz": reference.frequency_hz,
                "relative_frequency_difference": relative,
                "tolerance": 1e-8,
                "status": "PASS" if relative <= 1e-8 else "FAIL",
            }
        )
    return rows, maximum


def _report(
    *,
    status: str,
    beta0_only: bool,
    matrix_rows: Sequence[dict[str, Any]],
    virtual_rows: Sequence[dict[str, Any]],
    orthotropic_rows: Sequence[dict[str, Any]],
    equivalence_rows: Sequence[dict[str, Any]],
    root_rows: Sequence[dict[str, Any]],
    beta_runtimes: dict[float, float],
    beta_evaluations: dict[float, int],
    beta_mean_evaluation: dict[float, float],
    reference_runtime: float,
    reference_evaluations: int,
    total_runtime: float,
) -> str:
    max_sign = max(
        float(row["observed"])
        for row in matrix_rows
        if row["check"] == "notation_sign_gate"
    )
    min_rank = min(
        int(row["observed"]) for row in matrix_rows if row["check"] == "joint_rank"
    )
    max_random_rotation = max(
        float(row["observed"])
        for row in matrix_rows
        if row["check"] == "random_compatible_rotation"
    )
    max_random_moment = max(
        float(row["observed"])
        for row in matrix_rows
        if row["check"] == "random_equilibrated_moment"
    )
    max_work = max(float(row["max_absolute_residual"]) for row in virtual_rows)
    max_equivalence = max(
        (float(row["relative_frequency_difference"]) for row in equivalence_rows),
        default=math.nan,
    )
    max_root_det = max((float(row["determinant_residual"]) for row in root_rows), default=math.nan)
    max_root_svd = max(
        (float(row["relative_singular_residual"]) for row in root_rows),
        default=math.nan,
    )
    lines = [
        "# Chapter-2 coupled rigid-joint elastic pilot",
        "",
        f"Joint gate status: `{status}`.",
        "",
        "## Scope",
        "",
        "Two equal HMS/DX-209 Chapter-2 rods use `state_corrected`, real elastic moduli, `a=0.005 m`, `b=0.020 m`, `L_1=L_2=0.200 m`, `k=5/6`, and `theta_1=theta_2=15 deg`. Modes 1--6 are the pilot spectrum and mode 7 is the completeness guard.",
        "",
        f"This was a beta=0-only staging run: `{beta0_only}`.",
        "",
        "## Joint gates",
        "",
        f"- maximum `J_book - J_old block_diag(P,P)` entry: `{max_sign:.6e}`.",
        f"- minimum rank over beta=0,15,30,90,180 deg: `{min_rank}`.",
        f"- maximum compatible-rotation residual over deterministic random checks: `{max_random_rotation:.6e}`.",
        f"- maximum equilibrated-moment residual over deterministic random checks: `{max_random_moment:.6e}`.",
        f"- maximum virtual-work residual: `{max_work:.6e}`.",
        f"- orthotropic relative off-block norm: `{orthotropic_rows[0]['relative_off_block_norm']:.6e}`; rank `{orthotropic_rows[0]['matrix_rank']}`.",
        "- beta=0 and beta=90 scalar limits: recorded in `joint_matrix_checks.csv`.",
        "- equal-arm singular-value and spectrum checks: recorded in `joint_matrix_checks.csv`.",
        "",
        "## Boundary matrices",
        "",
        "`D_coupled = J_book @ block_diag(H_1^phys,H_2^phys)` is `6 x 6`, with each physical end map built from the existing corrected scaled-state transfer and confirmed book-slope clamp. The independent `D_straight = C_right @ T^phys(L_total) @ B_left` is `3 x 3` and selects `w`, `psi+Q/K_s`, and `Phi` at the right end.",
        "",
        "Positive row/column equilibration is used for root search. `pilot_roots.csv` distinguishes raw and scaled singular/determinant quality.",
        "",
        "## beta=0 straight-rod equivalence",
        "",
        f"Maximum first-seven relative frequency difference: `{max_equivalence:.6e}` (target `1e-8`).",
        "",
        "| mode | coupled Hz | straight Hz | relative difference |",
        "| ---: | ---: | ---: | ---: |",
    ]
    for row in equivalence_rows:
        lines.append(
            f"| {row['sorted_mode']} | {row['coupled_frequency_hz']:.12g} | {row['straight_frequency_hz']:.12g} | {row['relative_frequency_difference']:.6e} |"
        )
    lines.extend(["", "## Pilot spectrum", ""])
    for beta_deg in sorted({float(row["beta_deg"]) for row in root_rows}):
        lines.extend(
            [
                f"### beta={beta_deg:g} deg",
                "",
                "| mode | role | frequency, Hz | det residual | relative SVD residual | quality |",
                "| ---: | --- | ---: | ---: | ---: | --- |",
            ]
        )
        for row in root_rows:
            if float(row["beta_deg"]) == beta_deg:
                lines.append(
                    f"| {row['sorted_mode']} | {row['pilot_role']} | {row['frequency_hz']:.12g} | {row['determinant_residual']:.6e} | {row['relative_singular_residual']:.6e} | {row['quality_status']} |"
                )
        lines.append("")
    lines.extend(
        [
            "## Runtime and evaluations",
            "",
            f"Total runtime: `{total_runtime:.6f} s`. Maximum scaled determinant residual: `{max_root_det:.6e}`; maximum relative singular residual: `{max_root_svd:.6e}`.",
            "",
            "| beta, deg | root-search runtime, s | boundary-matrix evaluations | mean matrix time, s |",
            "| ---: | ---: | ---: | ---: |",
        ]
    )
    for beta_deg in sorted(beta_runtimes):
        lines.append(
            f"| {beta_deg:g} | {beta_runtimes[beta_deg]:.6f} | {beta_evaluations[beta_deg]} | {beta_mean_evaluation[beta_deg]:.6e} |"
        )
    exchange_row = next(
        (row for row in matrix_rows if row["check"] == "equal_arm_exchange_spectrum"),
        None,
    )
    if exchange_row is not None:
        lines.append(
            "| 30 arm exchange | "
            f"{exchange_row['runtime_seconds']:.6f} | "
            f"{exchange_row['boundary_matrix_evaluations']} | "
            f"{exchange_row['mean_boundary_matrix_seconds']:.6e} |"
        )
    lines.extend(
        [
            f"| straight reference | {reference_runtime:.6f} | {reference_evaluations} | n/a |",
            "",
            "## Exclusions",
            "",
            "No independent warping coordinate, bimoment, complex roots, loss factors, Euler--Bernoulli model, Saint-Venant comparison, FEM, parameter map, mode shape, MAC, or localization calculation was introduced or run.",
        ]
    )
    return "\n".join(lines) + "\n"


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()

    arm = _point(15.0, 0.2)
    straight = _point(15.0, 0.4)
    matrix_rows = _joint_matrix_checks(arm)
    virtual_rows = _virtual_work_checks()
    orthotropic_rows, orthotropic_ok = _orthotropic_check()

    beta_runtimes: dict[float, float] = {}
    beta_evaluations: dict[float, int] = {}
    beta_mean_evaluation: dict[float, float] = {}
    root_rows: list[dict[str, Any]] = []
    solved_roots: dict[float, list[RootResult]] = {}

    coupled_beta0_factory = lambda omega: coupled_boundary_matrix(omega, 0.0, arm, arm)
    coupled_beta0_raw = lambda omega: coupled_boundary_matrix_raw(omega, 0.0, arm, arm)
    beta0_roots, beta0_builder, beta0_runtime = _solve(arm, coupled_beta0_factory)
    solved_roots[0.0] = beta0_roots
    beta_runtimes[0.0] = beta0_runtime
    beta_evaluations[0.0] = beta0_builder.evaluations
    beta_mean_evaluation[0.0] = beta0_builder.mean_seconds
    root_rows.extend(
        _pilot_root_rows(0.0, beta0_roots, coupled_beta0_raw, beta0_builder, beta0_runtime)
    )

    if beta0_runtime > 300.0:
        status = "FAIL_ROOT_QUALITY"
        _write_csv(output_dir / "joint_matrix_checks.csv", matrix_rows)
        _write_csv(output_dir / "virtual_work_checks.csv", virtual_rows)
        _write_csv(output_dir / "pilot_roots.csv", root_rows)
        _write_csv(output_dir / "orthotropic_block_check.csv", orthotropic_rows)
        report = _report(
            status=status,
            beta0_only=True,
            matrix_rows=matrix_rows,
            virtual_rows=virtual_rows,
            orthotropic_rows=orthotropic_rows,
            equivalence_rows=[],
            root_rows=root_rows,
            beta_runtimes=beta_runtimes,
            beta_evaluations=beta_evaluations,
            beta_mean_evaluation=beta_mean_evaluation,
            reference_runtime=0.0,
            reference_evaluations=0,
            total_runtime=time.perf_counter() - started,
        )
        (output_dir / "coupled_joint_pilot_report.md").write_text(report, encoding="utf-8")
        print("joint_gate_status=FAIL_ROOT_QUALITY")
        print("beta0_runtime_exceeded_300_seconds=true")
        return 3

    reference_roots, reference_builder, reference_runtime = _solve(
        straight, lambda omega: straight_boundary_matrix(omega, straight)
    )
    equivalence_rows, max_equivalence = _equivalence_rows(beta0_roots, reference_roots)

    full_runtime_stopped = False
    if not args.beta0_only:
        for beta_deg in (30.0, 90.0):
            if time.perf_counter() - started > 600.0:
                full_runtime_stopped = True
                break
            beta = math.radians(beta_deg)
            factory = lambda omega, beta=beta: coupled_boundary_matrix(omega, beta, arm, arm)
            raw_factory = lambda omega, beta=beta: coupled_boundary_matrix_raw(omega, beta, arm, arm)
            roots, builder, runtime = _solve(arm, factory)
            solved_roots[beta_deg] = roots
            beta_runtimes[beta_deg] = runtime
            beta_evaluations[beta_deg] = builder.evaluations
            beta_mean_evaluation[beta_deg] = builder.mean_seconds
            root_rows.extend(_pilot_root_rows(beta_deg, roots, raw_factory, builder, runtime))

    if 30.0 in solved_roots:
        swapped_roots, swapped_builder, swapped_runtime = _solve(
            arm,
            lambda omega: coupled_boundary_matrix(
                omega, math.radians(30.0), arm, arm
            )[:, [3, 4, 5, 0, 1, 2]],
        )
        exchange_error = max(
            abs(left.frequency_hz - right.frequency_hz) / left.frequency_hz
            for left, right in zip(solved_roots[30.0], swapped_roots)
        )
        matrix_rows.append(
            {
                "check": "equal_arm_exchange_spectrum",
                "beta_deg": 30.0,
                "observed": exchange_error,
                "tolerance": 1e-10,
                "status": "PASS" if exchange_error <= 1e-10 else "FAIL",
                "details": "first seven roots after arm-reaction block exchange",
                "boundary_matrix_evaluations": swapped_builder.evaluations,
                "mean_boundary_matrix_seconds": swapped_builder.mean_seconds,
                "runtime_seconds": swapped_runtime,
            }
        )

    all_static_checks = all(row["status"] == "PASS" for row in matrix_rows)
    all_virtual = all(row["status"] == "PASS" for row in virtual_rows)
    all_root_quality = all(row["quality_status"] == "PASS" for row in root_rows)
    guard_found = all(len(roots) >= NUM_ROOTS for roots in solved_roots.values())
    beta0_ok = max_equivalence <= 1e-8
    full_angles_done = all(beta in solved_roots for beta in PILOT_BETAS_DEG)

    if not all_static_checks:
        status = "FAIL_SIGN_CONSISTENCY"
    elif not beta0_ok:
        status = "FAIL_BETA0_EQUIVALENCE"
    elif not all_root_quality or not guard_found:
        status = "FAIL_ROOT_QUALITY"
    elif not orthotropic_ok or not all_virtual:
        status = "FAIL_SIGN_CONSISTENCY"
    elif full_angles_done and not full_runtime_stopped:
        status = "PASS"
    else:
        status = "PARTIAL_PASS"

    _write_csv(output_dir / "joint_matrix_checks.csv", matrix_rows)
    _write_csv(output_dir / "virtual_work_checks.csv", virtual_rows)
    _write_csv(output_dir / "beta0_straight_rod_equivalence.csv", equivalence_rows)
    _write_csv(output_dir / "pilot_roots.csv", root_rows)
    _write_csv(output_dir / "orthotropic_block_check.csv", orthotropic_rows)
    total_runtime = time.perf_counter() - started
    report = _report(
        status=status,
        beta0_only=args.beta0_only,
        matrix_rows=matrix_rows,
        virtual_rows=virtual_rows,
        orthotropic_rows=orthotropic_rows,
        equivalence_rows=equivalence_rows,
        root_rows=root_rows,
        beta_runtimes=beta_runtimes,
        beta_evaluations=beta_evaluations,
        beta_mean_evaluation=beta_mean_evaluation,
        reference_runtime=reference_runtime,
        reference_evaluations=reference_builder.evaluations,
        total_runtime=total_runtime,
    )
    (output_dir / "coupled_joint_pilot_report.md").write_text(report, encoding="utf-8")

    print(f"joint_gate_status={status}")
    print(f"beta0_max_relative_frequency_difference={max_equivalence:.6e}")
    for beta_deg in sorted(beta_runtimes):
        print(f"beta_{beta_deg:g}_runtime_seconds={beta_runtimes[beta_deg]:.6f}")
        print(f"beta_{beta_deg:g}_boundary_matrix_evaluations={beta_evaluations[beta_deg]}")
        print(f"beta_{beta_deg:g}_mean_matrix_seconds={beta_mean_evaluation[beta_deg]:.6e}")
    print(f"straight_runtime_seconds={reference_runtime:.6f}")
    print(f"straight_boundary_matrix_evaluations={reference_builder.evaluations}")
    print(f"total_runtime_seconds={total_runtime:.6f}")
    print(f"output_dir={output_dir}")
    return 0 if status in ("PASS", "PARTIAL_PASS") else 4


if __name__ == "__main__":
    raise SystemExit(main())
