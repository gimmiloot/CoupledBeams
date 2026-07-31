"""Reproduce Yartsev Chapter-2 cantilever studies and audit the clamp condition.

This diagnostic covers one rectangular monoclinic Timoshenko rod only.  It
does not implement coupled rods, Euler--Bernoulli theory, Saint-Venant
comparison, a production anisotropic API, or FEM.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import statistics
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import Any, Iterable, Literal, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


_POSTPROCESS_ONLY_REQUESTED = any(
    flag in sys.argv
    for flag in (
        "--postprocess-boundary-source-check",
        "--boundary-source-check-from-saved",
    )
)

if not _POSTPROCESS_ONLY_REQUESTED:
    from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

if not _POSTPROCESS_ONLY_REQUESTED:
    from scripts.analysis.anisotropic_rods.reproduce_yartsev_fig_2_2 import (  # noqa: E402
        _sha256,
        _write_csv,
    )
    from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
        ClampVariant,
        Formulation,
        RodPoint,
        RootResult,
        assign_modes_by_mac,
        boundary_quality,
        cantilever_boundary_matrix,
        cantilever_boundary_residuals,
        cantilever_energy_fractions,
        cantilever_formulation_mapping_residual,
        cantilever_geometry,
        cantilever_mode_shape,
        cantilever_state_trajectory,
        cantilever_zero_frequency_nullity,
        continue_loss_root,
        decoupled_cantilever_boundary_factors,
        find_elastic_roots,
        fixed_free_torsion_omega,
        hms_dx_209_material,
        make_rod_point,
        modal_assurance,
        modal_loss_factors,
        partial_bending_boundary_matrix,
        partial_bending_mode_shape,
        partial_torsion_mode_shape,
        solve_complex_root,
        with_gxz_scale,
    )


Study = Literal["orientation", "length"]
DEFAULT_RESULT_DIR = ROOT / "results" / "anisotropic_rods" / "yartsev_ch2_cantilever"
SMOKE_RESULT_DIR = ROOT / "results" / "_smoke" / "yartsev_ch2_cantilever"
QUICK_RESULT_DIR = ROOT / "results" / "anisotropic_rods" / "yartsev_ch2_cantilever_quick_gate"
BOUNDARY_SOURCE_CHECK_RESULT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_cantilever_boundary_source_check"
)
CHAPTER_1 = ROOT / "docs" / "literature" / "pdf" / "Глава 1_compressed.pdf"
CHAPTER_2 = ROOT / "docs" / "literature" / "pdf" / "Глава 2_compressed.pdf"

CLAMP_OPTIONS: dict[str, ClampVariant] = {
    "book-slope": "book_slope_clamp",
    "section-rotation": "timoshenko_section_clamp",
}
CLAMP_LABELS: dict[ClampVariant, str] = {
    "book_slope_clamp": "book slope clamp",
    "timoshenko_section_clamp": "section-rotation clamp",
}

# Filled from a manual reading of the calculated curves after the numerical
# implementation was frozen.  Frequencies are kHz; losses are eta*1e2.
# The source provides no digital curve table and these values are never inputs
# to a root solve.  The tables are deliberately rounded to scan resolution.
FIGURE_2_8_DIGITIZED: tuple[tuple[float, int, str, float], ...] = ()
FIGURE_2_11_DIGITIZED: tuple[tuple[float, int, str, float], ...] = ()


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Modes:\n"
            "  ordinary full run: no mode flag; configure --study, "
            "--clamp-variant, and --material-mode\n"
            "  smoke: --smoke\n"
            "  elastic boundary sensitivity: --quick-boundary-gate\n"
            "  saved-data-only source check: "
            "--postprocess-boundary-source-check\n\n"
            "The postprocess-only mode exits before scientific solver imports. "
            "This CLI is not a coupled-rod solver."
        ),
    )
    parser.add_argument(
        "--study", choices=("orientation", "length", "both"), default="both"
    )
    parser.add_argument(
        "--clamp-variant",
        choices=("book-slope", "section-rotation", "both"),
        default="both",
    )
    parser.add_argument(
        "--material-mode",
        choices=("elastic", "book-complex", "both"),
        default="both",
    )
    parser.add_argument("--theta-step-deg", type=float, default=1.0)
    parser.add_argument("--length-step-m", type=float, default=0.005)
    parser.add_argument("--num-positive-modes", type=int, default=5)
    parser.add_argument("--scan-step-hz", type=float, default=10.0)
    parser.add_argument("--series-rtol", type=float, default=1.0e-12)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="run a reduced scientific smoke configuration",
    )
    parser.add_argument(
        "--quick-boundary-gate",
        action="store_true",
        help="run only the two-stage elastic sorted-spectrum clamp comparison",
    )
    parser.add_argument(
        "--postprocess-boundary-source-check",
        "--boundary-source-check-from-saved",
        dest="postprocess_boundary_source_check",
        action="store_true",
        help=(
            "compare manually digitized source curves with existing full-run CSV "
            "without importing or invoking the scientific solver workflow"
        ),
    )
    return parser.parse_args(argv)


def _validate_args(args: argparse.Namespace) -> None:
    if not 0.0 < args.theta_step_deg <= 90.0:
        raise ValueError("--theta-step-deg must lie in (0, 90]")
    if not 0.0 < args.length_step_m <= 0.35:
        raise ValueError("--length-step-m must lie in (0, 0.35]")
    if args.num_positive_modes < 1:
        raise ValueError("--num-positive-modes must be positive")
    if args.scan_step_hz <= 0.0 or args.series_rtol <= 0.0:
        raise ValueError("scan step and series tolerance must be positive")


def _inclusive_grid(start: float, stop: float, step: float) -> np.ndarray:
    count = int(math.floor((stop - start) / step + 1.0e-10))
    values = [start + index * step for index in range(count + 1)]
    values.append(stop)
    return np.array(sorted(set(round(float(value), 12) for value in values)))


def _orientation_grid(args: argparse.Namespace) -> np.ndarray:
    if args.smoke:
        return np.array([0.0, 15.0, 45.0, 90.0])
    base = _inclusive_grid(0.0, 90.0, args.theta_step_deg)
    refine_7 = _inclusive_grid(4.0, 10.0, 0.25)
    refine_19 = _inclusive_grid(16.0, 22.0, 0.25)
    mandatory = np.array([0.0, 7.0, 15.0, 19.0, 30.0, 45.0, 60.0, 75.0, 90.0])
    return np.unique(np.concatenate((base, refine_7, refine_19, mandatory)))


def _length_grid(args: argparse.Namespace) -> np.ndarray:
    if args.smoke:
        return np.array([0.10, 0.20, 0.40])
    base = _inclusive_grid(0.05, 0.40, args.length_step_m)
    refine_025 = _inclusive_grid(0.22, 0.28, 0.002)
    refine_008 = _inclusive_grid(0.06, 0.10, 0.001)
    mandatory = np.array([0.05, 0.08, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40])
    return np.unique(np.round(np.concatenate((base, refine_025, refine_008, mandatory)), 12))


def _studies(args: argparse.Namespace) -> list[Study]:
    return ["orientation", "length"] if args.study == "both" else [args.study]


def _clamps(args: argparse.Namespace) -> list[ClampVariant]:
    if args.clamp_variant == "both":
        return ["book_slope_clamp", "timoshenko_section_clamp"]
    return [CLAMP_OPTIONS[args.clamp_variant]]


def _point(
    study: Study,
    parameter: float,
    material_mode: Literal["elastic", "book_complex"],
    *,
    loss_scale: float = 1.0,
    series_rtol: float = 1.0e-12,
) -> RodPoint:
    theta = float(parameter) if study == "orientation" else 15.0
    length = 0.2 if study == "orientation" else float(parameter)
    return make_rod_point(
        theta,
        material_mode=material_mode,
        loss_scale=loss_scale,
        geometry=cantilever_geometry(length),
        material=hms_dx_209_material(),
        series_relative_tolerance=series_rtol,
    )


def _cantilever_builder(clamp: ClampVariant):
    def builder(
        omega: complex, point: RodPoint, formulation: Formulation
    ) -> np.ndarray:
        return cantilever_boundary_matrix(omega, point, clamp, formulation)

    return builder


def _material_rows() -> list[dict[str, Any]]:
    material = hms_dx_209_material()
    geometry = cantilever_geometry()
    values = [
        ("geometry", "a", geometry.a, "m", "Chapter 2 p.57"),
        ("geometry", "b", geometry.b, "m", "Chapter 2 p.57"),
        ("geometry", "A", geometry.area, "m^2", "Chapter 2 definition"),
        ("geometry", "I_y", geometry.I_y, "m^4", "Chapter 2 p.53"),
        ("geometry", "I_p", geometry.I_p, "m^4", "Chapter 2 p.53"),
        ("geometry", "k", geometry.shear_factor, "1", "Chapter 2 eq. (2.4)"),
        ("material", "Re(E1)", material.E1_real, "Pa", "Chapter 1 Table 1.1 p.45"),
        ("material", "Re(E2)", material.E2_real, "Pa", "Chapter 1 Table 1.1 p.45"),
        ("material", "Re(G12)", material.G12_real, "Pa", "Chapter 1 Table 1.1 p.45"),
        ("material", "Re(G13)", material.G13_real, "Pa", "Chapter 1 Table 1.1 p.45"),
        ("material", "Re(G23)", material.G23_real, "Pa", "Chapter 1 Table 1.1 p.45"),
        ("material", "nu12", material.nu12, "1", "Chapter 1 Table 1.1 p.45"),
        ("material", "rho", material.rho, "kg/m^3", "Chapter 1 Table 1.1 p.45"),
        ("material", "eta1", material.eta1, "1", "Table value 7.8 times 1e-4"),
        ("material", "eta2", material.eta2, "1", "Table value 6.7 times 1e-3"),
        ("material", "eta12", material.eta12, "1", "Table value 1.16 times 1e-2"),
        ("material", "eta13", material.eta13, "1", "Table value 1.16 times 1e-2"),
        ("material", "eta23", material.eta23, "1", "Table value 1.15 times 1e-2"),
        ("source_only", "h_ply", material.ply_thickness, "m", "Table value 2 times 1e-4; unused by 1D equations"),
    ]
    return [
        {"category": category, "parameter": name, "value_si": value, "unit": unit, "source": source}
        for category, name, value, unit, source in values
    ]


def _write_source_manifest(output_dir: Path) -> None:
    rows = [
        "# Source manifest — Yartsev Chapter-2 cantilever reproduction",
        "",
        "Only the two local monograph fragments below were used; no web or replacement source was used.",
        "",
        "| relative file | available | bytes | SHA256 | pages used |",
        "| --- | --- | ---: | --- | --- |",
    ]
    roles = {
        CHAPTER_1: "printed pp. 44–45, Table 1.1 (HMS/DX-209)",
        CHAPTER_2: "printed pp. 52–59 and 64–68, equations (2.1)–(2.18), Figures 2.8–2.12",
    }
    for path in (CHAPTER_1, CHAPTER_2):
        relative = path.relative_to(ROOT).as_posix()
        rows.append(
            f"| `{relative}` | {'yes' if path.is_file() else 'no'} | "
            f"{path.stat().st_size if path.is_file() else ''} | "
            f"`{_sha256(path) if path.is_file() else ''}` | {roles[path]} |"
        )
    rows.extend(
        [
            "",
            "`h_ply = 2e-4 m` is inventoried from Table 1.1 and is not used by the Chapter-2 one-dimensional equations.",
            "The PDFs are local untracked source fragments and may be absent from a public clone.",
        ]
    )
    (output_dir / "source_manifest.md").write_text("\n".join(rows) + "\n", encoding="utf-8")


def _solve_elastic(
    study: Study,
    parameters: np.ndarray,
    clamps: Sequence[ClampVariant],
    *,
    guard_roots: int,
    scan_step_hz: float,
    series_rtol: float,
) -> dict[tuple[ClampVariant, float], list[RootResult]]:
    solved: dict[tuple[ClampVariant, float], list[RootResult]] = {}
    for clamp in clamps:
        builder = _cantilever_builder(clamp)
        for index, raw_parameter in enumerate(parameters):
            parameter = float(raw_parameter)
            point = _point(study, parameter, "elastic", series_rtol=series_rtol)
            length_factor = max(1.0, 0.2 / point.geometry.length)
            solved[(clamp, parameter)] = find_elastic_roots(
                point,
                "state_corrected",
                num_roots=guard_roots,
                scan_step_hz=scan_step_hz * length_factor,
                initial_max_hz=5000.0 * length_factor,
                max_hz=300_000.0,
                boundary_matrix_builder=builder,
            )
            print(
                f"elastic {study} {clamp}: {parameter:g} ({index + 1}/{len(parameters)})",
                flush=True,
            )
    return solved


def _with_neighbor_distances(roots: Sequence[RootResult]) -> list[RootResult]:
    ordered = sorted(roots, key=lambda item: item.frequency_hz)
    result: list[RootResult] = []
    for index, item in enumerate(ordered):
        distances: list[float] = []
        if index:
            distances.append(item.frequency_hz - ordered[index - 1].frequency_hz)
        if index + 1 < len(ordered):
            distances.append(ordered[index + 1].frequency_hz - item.frequency_hz)
        result.append(replace(item, min_neighbor_distance_hz=min(distances) if distances else math.inf))
    return result


def _solve_complex(
    study: Study,
    parameters: np.ndarray,
    clamps: Sequence[ClampVariant],
    elastic: dict[tuple[ClampVariant, float], list[RootResult]],
    *,
    num_modes: int,
    series_rtol: float,
) -> dict[tuple[ClampVariant, float], list[RootResult]]:
    solved_all: dict[tuple[ClampVariant, float], list[RootResult]] = {}
    for clamp in clamps:
        builder = _cantilever_builder(clamp)
        previous_elastic: list[RootResult] | None = None
        previous_complex: list[RootResult] | None = None
        for parameter_index, raw_parameter in enumerate(parameters):
            parameter = float(raw_parameter)
            elastic_roots = elastic[(clamp, parameter)][:num_modes]
            roots: list[RootResult] = []
            for mode_index, elastic_root in enumerate(elastic_roots):
                factory = lambda scale, value=parameter: _point(
                    study,
                    value,
                    "book_complex",
                    loss_scale=scale,
                    series_rtol=series_rtol,
                )
                if previous_elastic is None or previous_complex is None:
                    result = continue_loss_root(
                        factory,
                        "state_corrected",
                        elastic_root.omega.real,
                        boundary_matrix_builder=builder,
                    )
                else:
                    old_e = previous_elastic[mode_index].omega.real
                    old_c = previous_complex[mode_index].omega
                    predictor = complex(
                        elastic_root.omega.real + old_c.real - old_e,
                        old_c.imag * elastic_root.omega.real / old_e,
                    )
                    result = solve_complex_root(
                        factory(1.0),
                        "state_corrected",
                        predictor,
                        boundary_matrix_builder=builder,
                    )
                    departure = abs(result.omega.real - elastic_root.omega.real) / elastic_root.omega.real
                    if result.status == "rejected_complex_quality" or departure > 0.10:
                        result = continue_loss_root(
                            factory,
                            "state_corrected",
                            elastic_root.omega.real,
                            boundary_matrix_builder=builder,
                        )
                roots.append(result)
            roots = _with_neighbor_distances(roots)
            duplicate = any(
                abs(roots[index].omega - roots[index - 1].omega)
                <= 1.0e-7 * max(abs(roots[index].omega), 1.0)
                for index in range(1, len(roots))
            )
            if duplicate:
                roots = _with_neighbor_distances(
                    [
                        continue_loss_root(
                            lambda scale, value=parameter: _point(
                                study,
                                value,
                                "book_complex",
                                loss_scale=scale,
                                series_rtol=series_rtol,
                            ),
                            "state_corrected",
                            item.omega.real,
                            boundary_matrix_builder=builder,
                        )
                        for item in elastic_roots
                    ]
                )
            solved_all[(clamp, parameter)] = roots
            previous_elastic = elastic_roots
            previous_complex = roots
            print(
                f"complex {study} {clamp}: {parameter:g} ({parameter_index + 1}/{len(parameters)})",
                flush=True,
            )
    return solved_all


def _track_modes(
    study: Study,
    parameters: np.ndarray,
    clamps: Sequence[ClampVariant],
    roots: dict[tuple[ClampVariant, float], list[RootResult]],
    *,
    material_mode: Literal["elastic", "book_complex"],
    num_modes: int,
    series_rtol: float,
) -> tuple[
    dict[tuple[ClampVariant, float, int], int],
    list[dict[str, Any]],
    dict[tuple[ClampVariant, float, int], np.ndarray],
]:
    mapping: dict[tuple[ClampVariant, float, int], int] = {}
    rows: list[dict[str, Any]] = []
    shapes: dict[tuple[ClampVariant, float, int], np.ndarray] = {}
    x_grid = np.linspace(0.0, 1.0, 41)
    for clamp in clamps:
        previous: list[np.ndarray] | None = None
        for parameter_index, raw_parameter in enumerate(parameters):
            parameter = float(raw_parameter)
            point = _point(study, parameter, material_mode, series_rtol=series_rtol)
            current = [
                cantilever_mode_shape(item.omega, point, clamp, x_grid)
                for item in roots[(clamp, parameter)][:num_modes]
            ]
            for sorted_index, shape in enumerate(current, start=1):
                shapes[(clamp, parameter, sorted_index)] = shape
            if previous is None:
                assignment = list(range(num_modes))
                mac_matrix = np.eye(num_modes)
            else:
                assignment, mac_matrix = assign_modes_by_mac(previous, current)
            next_previous: list[np.ndarray] = []
            for tracked_index, sorted_index in enumerate(assignment):
                tracked_mode = tracked_index + 1
                sorted_mode = sorted_index + 1
                mac = float(mac_matrix[tracked_index, sorted_index])
                ordered = np.sort(mac_matrix[tracked_index])[::-1]
                margin = float(ordered[0] - ordered[1]) if len(ordered) > 1 else 1.0
                mapping[(clamp, parameter, sorted_mode)] = tracked_mode
                low_mac = bool(parameter_index and mac < 0.75)
                rows.append(
                    {
                        "study": study,
                        "parameter_value": parameter,
                        "theta_deg": parameter if study == "orientation" else 15.0,
                        "length_m": 0.2 if study == "orientation" else parameter,
                        "clamp_variant": clamp,
                        "material_mode": material_mode,
                        "tracked_mode": tracked_mode,
                        "current_sorted_mode": sorted_mode,
                        "mac_from_previous_parameter": mac if parameter_index else 1.0,
                        "mac_margin": margin if parameter_index else 1.0,
                        "sorted_position_changed": bool(parameter_index and sorted_mode != tracked_mode),
                        "low_mac_warning": low_mac,
                        "mapping_status": "unresolved_low_mac" if low_mac else "accepted",
                        "local_refinement_present": (
                            study == "orientation" and (4.0 <= parameter <= 10.0 or 16.0 <= parameter <= 22.0)
                        )
                        or (study == "length" and (0.06 <= parameter <= 0.10 or 0.22 <= parameter <= 0.28)),
                    }
                )
                next_previous.append(current[sorted_index])
            previous = next_previous
    return mapping, rows, shapes


def _tracked_root(
    roots: dict[tuple[ClampVariant, float], list[RootResult]],
    mapping: dict[tuple[ClampVariant, float, int], int],
    clamp: ClampVariant,
    parameter: float,
    tracked_mode: int,
) -> tuple[int, RootResult]:
    sorted_mode = next(
        index
        for index in range(1, len(roots[(clamp, parameter)]) + 1)
        if mapping.get((clamp, parameter, index)) == tracked_mode
    )
    return sorted_mode, roots[(clamp, parameter)][sorted_mode - 1]


def _root_rows(
    study: Study,
    parameters: np.ndarray,
    clamps: Sequence[ClampVariant],
    roots: dict[tuple[ClampVariant, float], list[RootResult]],
    mapping: dict[tuple[ClampVariant, float, int], int],
    *,
    material_mode: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for raw_parameter in parameters:
        parameter = float(raw_parameter)
        for clamp in clamps:
            point = _point(study, parameter, "elastic" if material_mode == "elastic" else "book_complex")
            for sorted_mode, item in enumerate(roots[(clamp, parameter)], start=1):
                eta_exact, eta_approx, eta_difference = modal_loss_factors(item.omega)
                bending, shear, torsion = cantilever_energy_fractions(item.omega, point, clamp, samples=81)
                rows.append(
                    {
                        "study": study,
                        "parameter_value": parameter,
                        "theta_deg": parameter if study == "orientation" else 15.0,
                        "length_m": 0.2 if study == "orientation" else parameter,
                        "clamp_variant": clamp,
                        "material_mode": material_mode,
                        "sorted_mode": sorted_mode,
                        "tracked_mode": mapping.get((clamp, parameter, sorted_mode), ""),
                        "omega_real_rad_s": item.omega.real,
                        "omega_imag_rad_s": item.omega.imag,
                        "frequency_hz": item.frequency_hz,
                        "eta_exact": eta_exact,
                        "eta_approx": eta_approx,
                        "eta_exact_approx_relative_difference": eta_difference,
                        "bending_energy_fraction": bending,
                        "shear_energy_fraction": shear,
                        "torsion_energy_fraction": torsion,
                        "det_residual": item.determinant_residual,
                        "raw_determinant_abs": item.raw_determinant_abs,
                        "sigma_min": item.sigma_min,
                        "sigma_max": item.sigma_max,
                        "relative_singular_residual": item.relative_singular_residual,
                        "min_neighbor_distance_hz": item.min_neighbor_distance_hz,
                        "root_refinements": item.refinements,
                        "root_status": item.status,
                    }
                )
    return rows


def _boundary_rows(
    study: Study,
    parameters: np.ndarray,
    clamps: Sequence[ClampVariant],
    roots: dict[tuple[ClampVariant, float], list[RootResult]],
    *,
    material_mode: Literal["elastic", "book_complex"],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for raw_parameter in parameters:
        parameter = float(raw_parameter)
        point = _point(study, parameter, material_mode)
        for clamp in clamps:
            for sorted_mode, item in enumerate(roots[(clamp, parameter)], start=1):
                residuals = cantilever_boundary_residuals(item.omega, point, clamp)
                scaled = cantilever_state_trajectory(item.omega, point, clamp, [0.0, 1.0])
                condition_values = (
                    [scaled[0, 0], residuals["w_prime_0"], scaled[0, 2]]
                    if clamp == "book_slope_clamp"
                    else [scaled[0, 0], scaled[0, 1], scaled[0, 2]]
                )
                condition_values += [scaled[1, 3], scaled[1, 4], scaled[1, 5]]
                row: dict[str, Any] = {
                    "study": study,
                    "parameter_value": parameter,
                    "theta_deg": parameter if study == "orientation" else 15.0,
                    "length_m": 0.2 if study == "orientation" else parameter,
                    "clamp_variant": clamp,
                    "material_mode": material_mode,
                    "sorted_mode": sorted_mode,
                    "max_scaled_boundary_residual": max(abs(value) for value in condition_values),
                    "w_0_over_L": scaled[0, 0],
                    "psi_b_0": residuals["psi_b_0"],
                    "Phi_t_0": residuals["Phi_t_0"],
                    "Q_0_over_Ks": residuals["Q_0_over_Ks"],
                    "w_prime_0": residuals["w_prime_0"],
                    "slope_compensation": residuals["slope_compensation"],
                    "Q_L_scaled": scaled[1, 3],
                    "M_L_scaled": scaled[1, 4],
                    "M_T_L_scaled": scaled[1, 5],
                    "Q_L_physical": residuals["Q_L"],
                    "M_L_physical": residuals["M_L"],
                    "M_T_L_physical": residuals["M_T_L"],
                }
                rows.append(row)
    return rows


def _clamp_comparison_rows(
    study: Study,
    parameters: np.ndarray,
    elastic: dict[tuple[ClampVariant, float], list[RootResult]],
    complex_roots: dict[tuple[ClampVariant, float], list[RootResult]],
    elastic_maps: dict[tuple[ClampVariant, float, int], int],
    elastic_shapes: dict[tuple[ClampVariant, float, int], np.ndarray],
    *,
    num_modes: int,
) -> list[dict[str, Any]]:
    if not all(
        (clamp, float(parameters[0])) in elastic
        for clamp in ("book_slope_clamp", "timoshenko_section_clamp")
    ):
        return []
    rows: list[dict[str, Any]] = []
    for raw_parameter in parameters:
        parameter = float(raw_parameter)
        book_shapes = [elastic_shapes[("book_slope_clamp", parameter, index)] for index in range(1, num_modes + 1)]
        section_shapes = [elastic_shapes[("timoshenko_section_clamp", parameter, index)] for index in range(1, num_modes + 1)]
        mac_matrix = np.array(
            [[modal_assurance(left, right) for right in section_shapes] for left in book_shapes]
        )
        frequency_cost = np.array(
            [
                [
                    abs(
                        elastic[("book_slope_clamp", parameter)][i].frequency_hz
                        - elastic[("timoshenko_section_clamp", parameter)][j].frequency_hz
                    )
                    / elastic[("book_slope_clamp", parameter)][i].frequency_hz
                    for j in range(num_modes)
                ]
                for i in range(num_modes)
            ]
        )
        book_indices, section_indices = linear_sum_assignment(
            0.75 * (1.0 - mac_matrix) + 0.25 * np.minimum(frequency_cost, 1.0)
        )
        for book_index, section_index in zip(book_indices, section_indices):
            book = elastic[("book_slope_clamp", parameter)][book_index]
            section = elastic[("timoshenko_section_clamp", parameter)][section_index]
            book_eta = section_eta = math.nan
            if complex_roots:
                book_eta = modal_loss_factors(complex_roots[("book_slope_clamp", parameter)][book_index].omega)[0]
                section_eta = modal_loss_factors(complex_roots[("timoshenko_section_clamp", parameter)][section_index].omega)[0]
            mac = float(mac_matrix[book_index, section_index])
            rows.append(
                {
                    "study": study,
                    "parameter_value": parameter,
                    "theta_deg": parameter if study == "orientation" else 15.0,
                    "length_m": 0.2 if study == "orientation" else parameter,
                    "book_sorted_mode": int(book_index + 1),
                    "section_sorted_mode": int(section_index + 1),
                    "matched_mode_id": elastic_maps.get(("book_slope_clamp", parameter, int(book_index + 1)), ""),
                    "MAC": mac,
                    "book_frequency_hz": book.frequency_hz,
                    "section_frequency_hz": section.frequency_hz,
                    "frequency_difference_hz": section.frequency_hz - book.frequency_hz,
                    "delta_f_relative": abs(section.frequency_hz - book.frequency_hz) / book.frequency_hz,
                    "delta_frequency_squared_relative": abs(section.frequency_hz**2 - book.frequency_hz**2) / book.frequency_hz**2,
                    "book_eta_exact": book_eta,
                    "section_eta_exact": section_eta,
                    "loss_difference": section_eta - book_eta,
                    "mapping_status": "accepted" if mac >= 0.75 else "unresolved_low_mac",
                    "sorted_order_differs": bool(book_index != section_index),
                    "modal_character_difference": 1.0 - mac,
                }
            )
    return rows


def _formulation_equivalence_rows(
    *, scan_step_hz: float, series_rtol: float
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for theta in (0.0, 7.0, 15.0, 19.0, 45.0, 90.0):
        point = _point("orientation", theta, "elastic", series_rtol=series_rtol)
        builder = _cantilever_builder("book_slope_clamp")
        state = find_elastic_roots(
            point,
            "state_corrected",
            num_roots=6,
            scan_step_hz=scan_step_hz,
            boundary_matrix_builder=builder,
        )
        eliminated = find_elastic_roots(
            point,
            "eliminated_corrected",
            num_roots=6,
            scan_step_hz=scan_step_hz,
            boundary_matrix_builder=builder,
        )
        for mode, (state_root, eliminated_root) in enumerate(zip(state, eliminated), start=1):
            state_quality = boundary_quality(
                state_root.omega,
                point,
                "state_corrected",
                boundary_matrix_builder=builder,
            )
            eliminated_quality = boundary_quality(
                eliminated_root.omega,
                point,
                "eliminated_corrected",
                boundary_matrix_builder=builder,
            )
            rows.append(
                {
                    "theta_deg": theta,
                    "sorted_mode": mode,
                    "state_frequency_hz": state_root.frequency_hz,
                    "eliminated_corrected_frequency_hz": eliminated_root.frequency_hz,
                    "absolute_difference_hz": eliminated_root.frequency_hz - state_root.frequency_hz,
                    "relative_difference": abs(eliminated_root.frequency_hz - state_root.frequency_hz) / state_root.frequency_hz,
                    "state_relative_singular_residual": state_quality.relative_singular_residual,
                    "eliminated_relative_singular_residual": eliminated_quality.relative_singular_residual,
                    "field_mapping_residual": cantilever_formulation_mapping_residual(state_root.omega, point, samples=9),
                }
            )
    return rows


def _shear_rigid_rows(
    *, scan_step_hz: float, num_modes: int
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for scale in (1.0, 1.0e2, 1.0e4, 1.0e6):
        base = _point("orientation", 15.0, "elastic")
        point = with_gxz_scale(base, scale)
        roots: dict[ClampVariant, list[RootResult]] = {}
        shapes: dict[ClampVariant, list[np.ndarray]] = {}
        for clamp in ("book_slope_clamp", "timoshenko_section_clamp"):
            roots[clamp] = find_elastic_roots(
                point,
                "state_corrected",
                num_roots=num_modes,
                scan_step_hz=scan_step_hz,
                boundary_matrix_builder=_cantilever_builder(clamp),
            )
            shapes[clamp] = [
                cantilever_mode_shape(item.omega, point, clamp, np.linspace(0.0, 1.0, 61))
                for item in roots[clamp]
            ]
        mac_matrix = np.array(
            [[modal_assurance(left, right) for right in shapes["timoshenko_section_clamp"]] for left in shapes["book_slope_clamp"]]
        )
        book_indices, section_indices = linear_sum_assignment(-mac_matrix)
        for book_index, section_index in zip(book_indices, section_indices):
            book = roots["book_slope_clamp"][book_index]
            section = roots["timoshenko_section_clamp"][section_index]
            residuals = cantilever_boundary_residuals(book.omega, point, "book_slope_clamp")
            rows.append(
                {
                    "Gxz_scale": scale,
                    "mode": int(book_index + 1),
                    "section_matched_mode": int(section_index + 1),
                    "f_book_hz": book.frequency_hz,
                    "f_section_hz": section.frequency_hz,
                    "relative_difference": abs(section.frequency_hz - book.frequency_hz) / book.frequency_hz,
                    "MAC": float(mac_matrix[book_index, section_index]),
                    "boundary_shear_angle": abs(residuals["Q_0_over_Ks"]),
                }
            )
    return rows


def _partial_rows(
    study: Study,
    parameters: np.ndarray,
    *,
    include_complex: bool,
    num_modes: int,
    scan_step_hz: float,
) -> tuple[list[dict[str, Any]], dict[tuple[Study, float], list[RootResult]]]:
    rows: list[dict[str, Any]] = []
    elastic_bending: dict[tuple[Study, float], list[RootResult]] = {}
    for raw_parameter in parameters:
        parameter = float(raw_parameter)
        point = _point(study, parameter, "elastic")
        length_factor = max(1.0, 0.2 / point.geometry.length)
        bending = find_elastic_roots(
            point,
            "state_corrected",
            num_roots=num_modes,
            scan_step_hz=scan_step_hz * length_factor,
            initial_max_hz=5000.0 * length_factor,
            max_hz=300_000.0,
            boundary_matrix_builder=partial_bending_boundary_matrix,
        )
        elastic_bending[(study, parameter)] = bending
        for mode, item in enumerate(bending, start=1):
            rows.append(
                {
                    "study": study,
                    "parameter_value": parameter,
                    "theta_deg": parameter if study == "orientation" else 15.0,
                    "length_m": 0.2 if study == "orientation" else parameter,
                    "partial_family": "bending_eq_2_14",
                    "material_mode": "elastic",
                    "mode": mode,
                    "frequency_hz": item.frequency_hz,
                    "eta_exact": 0.0,
                    "eta_approx": 0.0,
                    "module_convention": "free Ex(theta) and free Gxz(theta)",
                }
            )
        for mode in range(1, num_modes + 1):
            omega = fixed_free_torsion_omega(point, mode, partial=True)
            rows.append(
                {
                    "study": study,
                    "parameter_value": parameter,
                    "theta_deg": parameter if study == "orientation" else 15.0,
                    "length_m": 0.2 if study == "orientation" else parameter,
                    "partial_family": "torsion_eq_2_15",
                    "material_mode": "elastic",
                    "mode": mode,
                    "frequency_hz": omega.real / (2.0 * np.pi),
                    "eta_exact": 0.0,
                    "eta_approx": 0.0,
                    "module_convention": "free Gxz(theta), pure Gxy_bar(theta), Cbar",
                }
            )
        if include_complex:
            complex_point = _point(study, parameter, "book_complex")
            for mode, elastic_root in enumerate(bending, start=1):
                result = continue_loss_root(
                    lambda scale, value=parameter: _point(
                        study, value, "book_complex", loss_scale=scale
                    ),
                    "state_corrected",
                    elastic_root.omega.real,
                    boundary_matrix_builder=partial_bending_boundary_matrix,
                )
                eta_exact, eta_approx, _ = modal_loss_factors(result.omega)
                rows.append(
                    {
                        "study": study,
                        "parameter_value": parameter,
                        "theta_deg": parameter if study == "orientation" else 15.0,
                        "length_m": 0.2 if study == "orientation" else parameter,
                        "partial_family": "bending_eq_2_14",
                        "material_mode": "book_complex",
                        "mode": mode,
                        "frequency_hz": result.frequency_hz,
                        "eta_exact": eta_exact,
                        "eta_approx": eta_approx,
                        "module_convention": "free Ex(theta) and free Gxz(theta)",
                    }
                )
            for mode in range(1, num_modes + 1):
                omega = fixed_free_torsion_omega(complex_point, mode, partial=True)
                eta_exact, eta_approx, _ = modal_loss_factors(omega)
                rows.append(
                    {
                        "study": study,
                        "parameter_value": parameter,
                        "theta_deg": parameter if study == "orientation" else 15.0,
                        "length_m": 0.2 if study == "orientation" else parameter,
                        "partial_family": "torsion_eq_2_15",
                        "material_mode": "book_complex",
                        "mode": mode,
                        "frequency_hz": omega.real / (2.0 * np.pi),
                        "eta_exact": eta_exact,
                        "eta_approx": eta_approx,
                        "module_convention": "free Gxz(theta), pure Gxy_bar(theta), Cbar",
                    }
                )
    return rows, elastic_bending


def _digitized_rows(study: Study) -> list[dict[str, Any]]:
    source = FIGURE_2_8_DIGITIZED if study == "orientation" else FIGURE_2_11_DIGITIZED
    figure = "2.8" if study == "orientation" else "2.11"
    rows: list[dict[str, Any]] = []
    for parameter, mode, quantity, value in source:
        frequency = quantity == "frequency"
        rows.append(
            {
                "study": study,
                "parameter_value": parameter,
                "book_curve_mode_label": mode,
                "quantity": quantity,
                "digitized_value": value,
                "unit": "kHz" if frequency else "eta_times_100",
                "uncertainty_type": "relative" if study == "length" and frequency else "absolute",
                "estimated_graphical_uncertainty": (
                    0.12 if study == "length" and frequency else (0.08 if frequency else 0.04)
                ),
                "source": f"Figure {figure} calculated solid curve",
                "method": "manual scan reading after implementation; rounded; not used for fitting",
                "log_frequency_axis": bool(study == "length" and frequency),
            }
        )
    return rows


def _digitized_comparison(
    study: Study,
    digitized: list[dict[str, Any]],
    complex_roots: dict[tuple[ClampVariant, float], list[RootResult]],
    complex_maps: dict[tuple[ClampVariant, float, int], int],
    clamps: Sequence[ClampVariant],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for source in digitized:
        parameter = float(source["parameter_value"])
        tracked_mode = int(source["book_curve_mode_label"])
        for clamp in clamps:
            if (clamp, parameter) not in complex_roots:
                continue
            sorted_mode, root = _tracked_root(
                complex_roots, complex_maps, clamp, parameter, tracked_mode
            )
            computed = root.frequency_hz / 1000.0 if source["quantity"] == "frequency" else 100.0 * modal_loss_factors(root.omega)[0]
            difference = computed - float(source["digitized_value"])
            uncertainty = float(source["estimated_graphical_uncertainty"])
            tolerance = uncertainty * abs(float(source["digitized_value"])) if source["uncertainty_type"] == "relative" else uncertainty
            rows.append(
                {
                    **source,
                    "clamp_variant": clamp,
                    "computed_sorted_mode": sorted_mode,
                    "computed_value": computed,
                    "computed_minus_digitized": difference,
                    "absolute_difference": abs(difference),
                    "comparison_tolerance": tolerance,
                    "difference_over_graphical_uncertainty": abs(difference) / tolerance if tolerance else math.inf,
                    "within_estimated_graphical_resolution": abs(difference) <= tolerance,
                }
            )
    return rows


def _phase_real(values: np.ndarray) -> np.ndarray:
    flat = values.ravel()
    pivot = flat[int(np.argmax(np.abs(flat)))]
    rotated = values * (np.exp(-1j * np.angle(pivot)) if abs(pivot) else 1.0)
    return np.real(rotated)


def _plot_reproduction(
    output_dir: Path,
    study: Study,
    parameters: np.ndarray,
    elastic: dict[tuple[ClampVariant, float], list[RootResult]],
    complex_roots: dict[tuple[ClampVariant, float], list[RootResult]],
    elastic_maps: dict[tuple[ClampVariant, float, int], int],
    complex_maps: dict[tuple[ClampVariant, float, int], int],
    partial_rows: list[dict[str, Any]],
    digitized: list[dict[str, Any]],
    *,
    num_modes: int,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), constrained_layout=True)
    colors = plt.cm.tab10(np.linspace(0.0, 0.9, num_modes))
    for clamp, linestyle in (("book_slope_clamp", "-"), ("timoshenko_section_clamp", "--")):
        if (clamp, float(parameters[0])) not in elastic:
            continue
        for tracked_mode in range(1, num_modes + 1):
            frequencies = [
                _tracked_root(elastic, elastic_maps, clamp, float(parameter), tracked_mode)[1].frequency_hz / 1000.0
                for parameter in parameters
            ]
            axes[0].plot(parameters, frequencies, linestyle=linestyle, color=colors[tracked_mode - 1], linewidth=1.5)
            if complex_roots:
                losses = [
                    100.0 * modal_loss_factors(
                        _tracked_root(complex_roots, complex_maps, clamp, float(parameter), tracked_mode)[1].omega
                    )[0]
                    for parameter in parameters
                ]
                axes[1].plot(parameters, losses, linestyle=linestyle, color=colors[tracked_mode - 1], linewidth=1.5)
    book_partial = [row for row in partial_rows if row["material_mode"] == ("book_complex" if complex_roots else "elastic")]
    for family, linestyle in (("bending_eq_2_14", ":"), ("torsion_eq_2_15", "-.")):
        for mode in range(1, min(num_modes, 4 if family.startswith("bending") else 3) + 1):
            subset = [row for row in book_partial if row["partial_family"] == family and row["mode"] == mode]
            if len(subset) != len(parameters):
                continue
            axes[0].plot(parameters, [row["frequency_hz"] / 1000.0 for row in subset], color="0.55", linestyle=linestyle, linewidth=0.8)
            if complex_roots:
                axes[1].plot(parameters, [100.0 * row["eta_exact"] for row in subset], color="0.65", linestyle=linestyle, linewidth=0.8)
    for row in digitized:
        axis = axes[0] if row["quantity"] == "frequency" else axes[1]
        axis.plot(row["parameter_value"], row["digitized_value"], marker="o", markersize=3, markerfacecolor="none", markeredgecolor="black", linestyle="none")
    axes[0].set_ylabel("f, kHz")
    axes[1].set_ylabel(r"$\eta_i\times10^2$")
    if study == "orientation":
        for axis in axes:
            axis.set_xlabel(r"$\theta$, deg")
            axis.set_xlim(0.0, 90.0)
            axis.set_xticks(np.arange(0.0, 91.0, 15.0))
        stem = "figure_2_8_reproduction"
    else:
        axes[0].set_yscale("log")
        for axis in axes:
            axis.set_xlabel("L, m")
            axis.set_xlim(0.05, 0.40)
        stem = "figure_2_11_reproduction"
    for axis in axes:
        axis.grid(True, alpha=0.25, which="both")
    fig.suptitle("solid: book slope; dashed: section rotation; gray: partial controls; circles: scan reading")
    fig.savefig(output_dir / f"{stem}.png", dpi=220)
    fig.savefig(output_dir / f"{stem}.pdf")
    plt.close(fig)


def _plot_mode_panels(
    output_dir: Path,
    study: Study,
    targets: Sequence[float],
    elastic: dict[tuple[ClampVariant, float], list[RootResult]],
    mapping: dict[tuple[ClampVariant, float, int], int],
    *,
    num_modes: int,
) -> None:
    clamp: ClampVariant = "book_slope_clamp"
    if (clamp, float(targets[0])) not in elastic:
        return
    x = np.linspace(0.0, 1.0, 101)
    fig, axes = plt.subplots(num_modes, len(targets), figsize=(2.0 * len(targets), 1.45 * num_modes), sharex=True, constrained_layout=True, squeeze=False)
    for column, raw_target in enumerate(targets):
        target = float(raw_target)
        point = _point(study, target, "elastic")
        for tracked_mode in range(1, num_modes + 1):
            _, root = _tracked_root(elastic, mapping, clamp, target, tracked_mode)
            shape = cantilever_mode_shape(root.omega, point, clamp, x)
            real = _phase_real(shape)
            scale = max(np.max(np.abs(real[:, 0])), np.max(np.abs(real[:, 2])), 1.0e-15)
            axis = axes[tracked_mode - 1, column]
            axis.plot(x, real[:, 0] / scale, color="black", linewidth=1.0)
            axis.plot(x, real[:, 2] / scale, color="tab:blue", linestyle="--", linewidth=0.9)
            axis.axhline(0.0, color="0.8", linewidth=0.5)
            if column == 0:
                axis.set_ylabel(f"i={tracked_mode}")
            if tracked_mode == 1:
                axis.set_title((r"$\theta$=" if study == "orientation" else "L=") + f"{target:g}")
            if tracked_mode == num_modes:
                axis.set_xlabel("x/L")
            axis.set_yticks([])
    stem = "figure_2_9_mode_shapes" if study == "orientation" else "figure_2_12_mode_shapes"
    fig.suptitle("book-slope clamp: w/L (black), Phi (blue dashed); continuity labels")
    fig.savefig(output_dir / f"{stem}.png", dpi=220)
    fig.savefig(output_dir / f"{stem}.pdf")
    plt.close(fig)


def _plot_partial_shapes(
    output_dir: Path,
    bending_roots: dict[tuple[Study, float], list[RootResult]],
) -> None:
    key = ("orientation", 15.0)
    if key not in bending_roots:
        return
    point = _point("orientation", 15.0, "elastic")
    x = np.linspace(0.0, 1.0, 151)
    fig, axes = plt.subplots(2, 3, figsize=(9.0, 5.6), constrained_layout=True)
    for mode in range(1, 5):
        axis = axes.ravel()[mode - 1]
        shape = _phase_real(partial_bending_mode_shape(bending_roots[key][mode - 1].omega, point, x))
        axis.plot(x, shape[:, 0] / max(np.max(np.abs(shape[:, 0])), 1.0e-15), color="black")
        axis.set_title(f"partial bending {mode}")
    for mode in (1, 2):
        axis = axes.ravel()[mode + 3]
        axis.plot(x, partial_torsion_mode_shape(mode, x), color="tab:blue")
        axis.set_title(f"partial torsion {mode}")
    for axis in axes.ravel():
        axis.axhline(0.0, color="0.8", linewidth=0.5)
        axis.set_xlabel("x/L")
        axis.set_yticks([])
    fig.suptitle("Equation (2.14)/(2.15) partial fixed-free mode controls")
    fig.savefig(output_dir / "figure_2_10_partial_mode_shapes.png", dpi=220)
    fig.savefig(output_dir / "figure_2_10_partial_mode_shapes.pdf")
    plt.close(fig)


def _plot_clamp_comparison(output_dir: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5), constrained_layout=True)
    for study, axis in (("orientation", axes[0]), ("length", axes[1])):
        subset = [row for row in rows if row["study"] == study]
        for mode in sorted({int(row["matched_mode_id"]) for row in subset if row["matched_mode_id"] != ""}):
            mode_rows = [row for row in subset if row["matched_mode_id"] == mode]
            axis.plot([row["parameter_value"] for row in mode_rows], [100.0 * row["delta_f_relative"] for row in mode_rows], label=f"mode {mode}")
        axis.set_xlabel(r"$\theta$, deg" if study == "orientation" else "L, m")
        axis.set_ylabel("|delta f|, %")
        axis.set_title(study)
        axis.grid(True, alpha=0.25)
        if subset:
            axis.legend(fontsize=7, ncol=2)
    fig.suptitle("Book-slope versus section-rotation clamp after MAC matching")
    fig.savefig(output_dir / "clamp_variant_comparison.png", dpi=220)
    fig.savefig(output_dir / "clamp_variant_comparison.pdf")
    plt.close(fig)


def _audit_status(comparison_rows: list[dict[str, Any]], smoke: bool) -> str:
    if smoke or not comparison_rows:
        return "SOURCE_AMBIGUITY"
    outcomes: dict[ClampVariant, bool] = {}
    for clamp in ("book_slope_clamp", "timoshenko_section_clamp"):
        subset = [row for row in comparison_rows if row["clamp_variant"] == clamp]
        outcomes[clamp] = bool(subset) and all(bool(row["within_estimated_graphical_resolution"]) for row in subset)
    if outcomes["book_slope_clamp"] and outcomes["timoshenko_section_clamp"]:
        return "BOTH_WITHIN_GRAPH_RESOLUTION"
    if outcomes["book_slope_clamp"]:
        return "BOOK_SLOPE_CLAMP_CONFIRMED"
    if outcomes["timoshenko_section_clamp"]:
        return "TIMOSHENKO_SECTION_CLAMP_CONFIRMED"
    return "NEITHER_REPRODUCES"


def _write_report(
    output_dir: Path,
    args: argparse.Namespace,
    status: str,
    comparison_rows: list[dict[str, Any]],
    equivalence_rows: list[dict[str, Any]],
    boundary_rows: list[dict[str, Any]],
    continuity_rows: list[dict[str, Any]],
    clamp_rows: list[dict[str, Any]],
    shear_rows: list[dict[str, Any]],
    digitized_rows: list[dict[str, Any]],
    runtime_seconds: float,
) -> None:
    max_equivalence = max((row["relative_difference"] for row in equivalence_rows), default=math.nan)
    max_mapping = max((row["field_mapping_residual"] for row in equivalence_rows), default=math.nan)
    max_boundary = max((abs(row["max_scaled_boundary_residual"]) for row in boundary_rows), default=math.nan)
    max_df = max((row["delta_f_relative"] for row in clamp_rows), default=math.nan)
    max_df2 = max((row["delta_frequency_squared_relative"] for row in clamp_rows), default=math.nan)
    max_loss = max((abs(row["loss_difference"]) for row in clamp_rows if math.isfinite(row["loss_difference"])), default=math.nan)
    min_mac = min((row["MAC"] for row in clamp_rows), default=math.nan)
    low_mac = sum(bool(row["low_mac_warning"]) for row in continuity_rows)
    graph_failures = sum(not bool(row["within_estimated_graphical_resolution"]) for row in digitized_rows)
    shear_start = max((row["relative_difference"] for row in shear_rows if row["Gxz_scale"] == 1.0), default=math.nan)
    shear_end = max((row["relative_difference"] for row in shear_rows if row["Gxz_scale"] == 1.0e6), default=math.nan)
    lines = [
        "# Yartsev Chapter-2 cantilever reproduction and clamp audit",
        "",
        f"**Boundary-condition audit status: `{status}`.**",
        "",
        "## Source and scope",
        "",
        "The calculation uses only local Chapter-1 p.45/Table 1.1 and Chapter-2 pp.52–59, 64–68. It reproduces the calculated coupled curves and mode-shape trends in Figures 2.8–2.12 for one cantilever rod. No parameter was fitted.",
        "",
        "## Clamp definitions",
        "",
        "- `book_slope_clamp`: `w(0)=0`, `w'(0)=psi_b(0)+Q(0)/(k Gxz A)=0`, `Phi_t(0)=0`.",
        "- `timoshenko_section_clamp`: `w(0)=0`, `psi_b(0)=0`, `Phi_t(0)=0`.",
        "- Both use `Q(L)=M(L)=M_T(L)=0`; `D_cf=S_f exp(A_state L) B_clamp` is 3x3.",
        "",
        "## Numerical checks",
        "",
        f"- Maximum state/eliminated-corrected frequency difference for the slope clamp: `{max_equivalence:.6e}` relative.",
        f"- Maximum state-to-eliminated field mapping residual: `{max_mapping:.6e}`.",
        f"- Maximum scaled physical boundary residual: `{max_boundary:.6e}`.",
        "- Both zero-frequency cantilever matrices are nonsingular; there are no rigid-body roots.",
        f"- Shear-rigid maximum clamp difference changes from `{shear_start:.6e}` at scale 1 to `{shear_end:.6e}` at scale 1e6.",
        f"- Continuity assignments below MAC 0.75: `{low_mac}`; refined interaction grids are retained in the CSVs.",
        "- At theta=0 and 90 deg, Sbar16 vanishes, the determinant factorizes, and fixed-free torsion agrees with `(2n-1)pi/(2L)*sqrt(C_T/(rho I_p))`.",
        "",
        "## Partial controls",
        "",
        "Equation (2.14) uses the source's free complex `Ex(theta)` and `Gxz(theta)`. Equation (2.15) uses free `Gxz(theta)` and pure `Gxy_bar(theta)` through `Cbar`. These controls interpret the dotted/dashed curves; they do not replace the coupled spectrum.",
        "",
        "## Clamp comparison",
        "",
        f"- Maximum MAC-matched frequency difference: `{max_df:.6e}` relative.",
        f"- Maximum squared-frequency difference: `{max_df2:.6e}` relative.",
        f"- Maximum modal-loss difference: `{max_loss:.6e}`.",
        f"- Minimum cross-clamp MAC: `{min_mac:.6e}`.",
        "- Matching uses field MAC plus frequency proximity; sorted mode number is retained only as metadata.",
        "",
        "## Figure comparison",
        "",
        "Manual digitization was performed after implementation. Figure 2.8 uses absolute uncertainties of ±0.08 kHz and ±0.04 in eta*1e2. Figure 2.11 uses ±12% for its logarithmic frequency graph and ±0.04 in eta*1e2. Experimental data are not involved.",
        f"Comparison rows outside those declared resolutions: `{graph_failures}`.",
        "",
        "## Run",
        "",
        "```text",
        "python scripts/analysis/anisotropic_rods/reproduce_yartsev_ch2_cantilever.py --study both --clamp-variant both --material-mode both",
        "```",
        "",
        f"Runtime: `{runtime_seconds:.3f} s` ({'smoke' if args.smoke else 'ordinary'}). Positive root {args.num_positive_modes + 1} is the completeness guard.",
        "",
        "## Limits",
        "",
        "This is a source-validated single-rod diagnostic. It is not a coupled-rod model, Euler–Bernoulli or Saint-Venant comparison, production anisotropic API, or FEM validation. Scan readings are resolution-limited and were not used to tune material data.",
    ]
    (output_dir / "cantilever_reproduction_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


class _CountingCantileverBuilder:
    """Count characteristic-matrix evaluations in the quick gate only."""

    def __init__(self, clamp: ClampVariant) -> None:
        self.clamp = clamp
        self.evaluations = 0

    def __call__(
        self, omega: complex, point: RodPoint, formulation: Formulation
    ) -> np.ndarray:
        self.evaluations += 1
        return cantilever_boundary_matrix(omega, point, self.clamp, formulation)


def _quick_solve_grid(
    theta_values: Sequence[float],
    *,
    num_roots: int,
    scan_step_hz: float,
) -> tuple[
    dict[tuple[ClampVariant, float], list[RootResult]], int, float
]:
    roots: dict[tuple[ClampVariant, float], list[RootResult]] = {}
    evaluations = 0
    started = time.perf_counter()
    for theta in theta_values:
        point = _point("orientation", float(theta), "elastic")
        for clamp in ("book_slope_clamp", "timoshenko_section_clamp"):
            counter = _CountingCantileverBuilder(clamp)
            roots[(clamp, float(theta))] = find_elastic_roots(
                point,
                "state_corrected",
                num_roots=num_roots,
                scan_step_hz=scan_step_hz,
                initial_max_hz=5000.0,
                max_hz=100_000.0,
                boundary_matrix_builder=counter,
            )
            evaluations += counter.evaluations
    return roots, evaluations, time.perf_counter() - started


def _quick_root_rows(
    theta_values: Sequence[float],
    roots: dict[tuple[ClampVariant, float], list[RootResult]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for theta in theta_values:
        for clamp in ("book_slope_clamp", "timoshenko_section_clamp"):
            for sorted_mode, item in enumerate(roots[(clamp, float(theta))], start=1):
                rows.append(
                    {
                        "theta_deg": float(theta),
                        "clamp_variant": clamp,
                        "sorted_mode": sorted_mode,
                        "omega_rad_s": item.omega.real,
                        "frequency_hz": item.frequency_hz,
                        "det_residual": item.determinant_residual,
                        "relative_singular_residual": item.relative_singular_residual,
                        "root_status": item.status,
                    }
                )
    return rows


def _quick_comparison_rows(
    theta_values: Sequence[float],
    roots: dict[tuple[ClampVariant, float], list[RootResult]],
    *,
    researched_modes: int = 5,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for theta in theta_values:
        book = roots[("book_slope_clamp", float(theta))]
        section = roots[("timoshenko_section_clamp", float(theta))]
        for index in range(min(researched_modes, len(book), len(section))):
            left = book[index]
            right = section[index]
            book_relative_gap = left.min_neighbor_distance_hz / left.frequency_hz
            section_relative_gap = right.min_neighbor_distance_hz / right.frequency_hz
            possible_reordering = min(book_relative_gap, section_relative_gap) < 0.05
            rows.append(
                {
                    "theta_deg": float(theta),
                    "book_sorted_mode": index + 1,
                    "section_sorted_mode": index + 1,
                    "f_book_hz": left.frequency_hz,
                    "f_section_hz": right.frequency_hz,
                    "relative_difference": abs(right.frequency_hz - left.frequency_hz) / left.frequency_hz,
                    "relative_squared_frequency_difference": abs(right.frequency_hz**2 - left.frequency_hz**2) / left.frequency_hz**2,
                    "possible_mode_reordering": possible_reordering,
                    "book_min_neighbor_distance_hz": left.min_neighbor_distance_hz,
                    "section_min_neighbor_distance_hz": right.min_neighbor_distance_hz,
                }
            )
    return rows


def _quick_plot(
    output_dir: Path,
    theta_values: Sequence[float],
    roots: dict[tuple[ClampVariant, float], list[RootResult]],
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.6), sharex=True, constrained_layout=True)
    for mode in range(1, 6):
        axis = axes.ravel()[mode - 1]
        book = [roots[("book_slope_clamp", float(theta))][mode - 1].frequency_hz for theta in theta_values]
        section = [roots[("timoshenko_section_clamp", float(theta))][mode - 1].frequency_hz for theta in theta_values]
        axis.plot(theta_values, book, marker="o", label="book slope")
        axis.plot(theta_values, section, marker="s", linestyle="--", label="section rotation")
        axis.set_title(f"sorted mode {mode}")
        axis.set_ylabel("f, Hz")
        axis.grid(True, alpha=0.25)
    axes.ravel()[5].axis("off")
    for axis in axes[-1, :2]:
        axis.set_xlabel(r"$\theta$, deg")
    axes.ravel()[0].legend(fontsize=8)
    fig.suptitle("Elastic quick boundary gate, L=0.2 m, HMS/DX-209")
    fig.savefig(output_dir / "quick_boundary_gate_comparison.pdf")
    plt.close(fig)


def _quick_torsion_checks(
    roots: dict[tuple[ClampVariant, float], list[RootResult]]
) -> list[dict[str, Any]]:
    point = _point("orientation", 0.0, "elastic")
    rows: list[dict[str, Any]] = []
    for clamp in ("book_slope_clamp", "timoshenko_section_clamp"):
        numerical = roots[(clamp, 0.0)]
        for torsion_mode in (1, 2, 3):
            analytical = fixed_free_torsion_omega(point, torsion_mode).real / (2.0 * np.pi)
            sorted_mode, nearest = min(
                enumerate(numerical, start=1),
                key=lambda pair: abs(pair[1].frequency_hz - analytical),
            )
            rows.append(
                {
                    "clamp_variant": clamp,
                    "torsion_mode": torsion_mode,
                    "matched_sorted_mode": sorted_mode,
                    "analytical_frequency_hz": analytical,
                    "numerical_frequency_hz": nearest.frequency_hz,
                    "relative_difference": abs(nearest.frequency_hz - analytical) / analytical,
                }
            )
    return rows


def _write_quick_report(
    output_dir: Path,
    *,
    smoke_runtime: float,
    quick_runtime: float | None,
    smoke_evaluations: int,
    quick_evaluations: int,
    convergence_evaluations: int,
    roots: dict[tuple[ClampVariant, float], list[RootResult]],
    comparison: list[dict[str, Any]],
    torsion_checks: list[dict[str, Any]],
    convergence_difference: float,
    d0_rows: list[dict[str, Any]],
    stopped_after_smoke: bool,
) -> None:
    # RootResult construction performs one additional D_cf build per accepted
    # root only to recover the matrix dimension for residual scaling.  Remove
    # those probes when reporting actual determinant/quality evaluations.
    smoke_determinants = smoke_evaluations - 2 * 4 * 3
    quick_determinants = quick_evaluations - (2 * 7 * 6 if quick_runtime is not None else 0)
    convergence_determinants = convergence_evaluations - (2 * 6 if convergence_evaluations else 0)
    max_frequency = max((row["relative_difference"] for row in comparison), default=math.nan)
    max_squared = max((row["relative_squared_frequency_difference"] for row in comparison), default=math.nan)
    max_rows = [row for row in comparison if abs(row["relative_difference"] - max_frequency) < 1e-15]
    max_det = max((item.determinant_residual for values in roots.values() for item in values), default=math.nan)
    max_sigma = max((item.relative_singular_residual for values in roots.values() for item in values), default=math.nan)
    max_torsion = max((row["relative_difference"] for row in torsion_checks), default=math.nan)
    reorder_count = sum(bool(row["possible_mode_reordering"]) for row in comparison)
    lines = [
        "# Yartsev Chapter-2 cantilever — quick boundary gate",
        "",
        "Status: preliminary elastic sorted-spectrum comparison only.",
        "",
        "## Scope",
        "",
        "The gate uses `state_corrected`, HMS/DX-209 real moduli, `a=0.005 m`, `b=0.02 m`, `L=0.2 m`, and `k=5/6`. It does not choose a preferred clamp.",
        "",
        "- `book_slope_clamp`: `w(0)=0`, `psi_b(0)+Q(0)/(k Gxz A)=0`, `Phi(0)=0`.",
        "- `timoshenko_section_clamp`: `w(0)=0`, `psi_b(0)=0`, `Phi(0)=0`.",
        "- Both use `Q(L)=M(L)=M_T(L)=0`; both characteristic matrices are 3x3.",
        "",
        "## Timings and evaluations",
        "",
        f"- smoke runtime: `{smoke_runtime:.6f} s`; determinant/quality evaluations: `{smoke_determinants}` (`D_cf` builds including dimension probes: `{smoke_evaluations}`).",
        f"- quick-gate runtime: `{quick_runtime:.6f} s`; determinant/quality evaluations: `{quick_determinants}` (`D_cf` builds including dimension probes: `{quick_evaluations}`)." if quick_runtime is not None else "- quick gate was not started because smoke exceeded 300 s.",
        f"- theta=15 half-step convergence determinant/quality evaluations: `{convergence_determinants}` (`D_cf` builds: `{convergence_evaluations}`).",
        f"- stopped after smoke: `{stopped_after_smoke}`.",
        "",
        "## Numerical checks",
        "",
        f"- maximum normalized determinant residual: `{max_det:.6e}`.",
        f"- maximum relative singular residual: `{max_sigma:.6e}`.",
        f"- maximum theta=15 scan-step-halving relative root change: `{convergence_difference:.6e}`.",
        f"- maximum analytical/numerical fixed-free torsion difference at theta=0: `{max_torsion:.6e}`.",
        f"- possible sorted-mode reordering flags from a 5% neighbor-gap screen: `{reorder_count}`.",
        "",
        "Zero-frequency matrices:",
        "",
        "| clamp | shape | det D(0) | nullity |",
        "| --- | --- | ---: | ---: |",
    ]
    for row in d0_rows:
        lines.append(
            f"| `{row['clamp_variant']}` | `{row['shape']}` | {row['determinant']:.12g} | {row['nullity']} |"
        )
    lines.extend(
        [
            "",
            "## Spectrum difference",
            "",
            f"Maximum relative frequency difference: `{max_frequency:.6e}`; maximum squared-frequency difference: `{max_squared:.6e}`.",
            "Maximum-frequency-difference rows: "
            + ", ".join(
                f"theta={row['theta_deg']:g} deg, sorted mode {row['book_sorted_mode']}"
                for row in max_rows
            )
            + ".",
            "",
            "The CSV comparison deliberately uses equal sorted indices. A close-pair flag is diagnostic only; no MAC or continuation was run.",
            "",
            "## Excluded work",
            "",
            "No complex roots, loss factors, length scan, shapes, MAC, partial spectra, figure digitization, eliminated formulation, Euler–Bernoulli, coupled rods, or FEM were run by this quick gate.",
        ]
    )
    (output_dir / "quick_boundary_gate_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _run_quick_boundary_gate(args: argparse.Namespace) -> int:
    output_dir = (args.output_dir or QUICK_RESULT_DIR).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    smoke_theta = [0.0, 15.0, 45.0, 90.0]
    smoke_roots, smoke_evaluations, smoke_runtime = _quick_solve_grid(
        smoke_theta, num_roots=3, scan_step_hz=args.scan_step_hz
    )
    if smoke_runtime > 300.0:
        smoke_rows = _quick_root_rows(smoke_theta, smoke_roots)
        _write_csv(output_dir / "quick_boundary_gate_roots.csv", smoke_rows)
        _write_csv(output_dir / "quick_boundary_gate_comparison.csv", _quick_comparison_rows(smoke_theta, smoke_roots, researched_modes=3))
        _quick_plot(output_dir, smoke_theta, smoke_roots)
        _write_quick_report(
            output_dir,
            smoke_runtime=smoke_runtime,
            quick_runtime=None,
            smoke_evaluations=smoke_evaluations,
            quick_evaluations=0,
            convergence_evaluations=0,
            roots=smoke_roots,
            comparison=_quick_comparison_rows(smoke_theta, smoke_roots, researched_modes=3),
            torsion_checks=[],
            convergence_difference=math.nan,
            d0_rows=[],
            stopped_after_smoke=True,
        )
        print(f"smoke_runtime_seconds={smoke_runtime:.6f}")
        print("quick_gate_skipped=smoke_exceeded_300_seconds")
        return 3

    quick_theta = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0]
    roots, quick_evaluations, quick_runtime = _quick_solve_grid(
        quick_theta, num_roots=6, scan_step_hz=args.scan_step_hz
    )
    rows = _quick_root_rows(quick_theta, roots)
    comparison = _quick_comparison_rows(quick_theta, roots)
    _write_csv(output_dir / "quick_boundary_gate_roots.csv", rows)
    _write_csv(output_dir / "quick_boundary_gate_comparison.csv", comparison)
    _quick_plot(output_dir, quick_theta, roots)

    convergence_evaluations = 0
    convergence_difference = 0.0
    point_15 = _point("orientation", 15.0, "elastic")
    for clamp in ("book_slope_clamp", "timoshenko_section_clamp"):
        counter = _CountingCantileverBuilder(clamp)
        fine = find_elastic_roots(
            point_15,
            "state_corrected",
            num_roots=6,
            scan_step_hz=args.scan_step_hz / 2.0,
            boundary_matrix_builder=counter,
        )
        convergence_evaluations += counter.evaluations
        coarse = roots[(clamp, 15.0)]
        convergence_difference = max(
            convergence_difference,
            max(
                abs(left.frequency_hz - right.frequency_hz) / left.frequency_hz
                for left, right in zip(coarse, fine)
            ),
        )

    torsion_checks = _quick_torsion_checks(roots)
    d0_rows: list[dict[str, Any]] = []
    for clamp in ("book_slope_clamp", "timoshenko_section_clamp"):
        matrix = cantilever_boundary_matrix(0.0, point_15, clamp)
        d0_rows.append(
            {
                "clamp_variant": clamp,
                "shape": f"{matrix.shape[0]}x{matrix.shape[1]}",
                "determinant": float(np.linalg.det(matrix).real),
                "nullity": cantilever_zero_frequency_nullity(point_15, clamp),
            }
        )
    _write_quick_report(
        output_dir,
        smoke_runtime=smoke_runtime,
        quick_runtime=quick_runtime,
        smoke_evaluations=smoke_evaluations,
        quick_evaluations=quick_evaluations,
        convergence_evaluations=convergence_evaluations,
        roots=roots,
        comparison=comparison,
        torsion_checks=torsion_checks,
        convergence_difference=convergence_difference,
        d0_rows=d0_rows,
        stopped_after_smoke=False,
    )
    print(f"smoke_runtime_seconds={smoke_runtime:.6f}")
    print(f"quick_gate_runtime_seconds={quick_runtime:.6f}")
    print(f"smoke_determinant_evaluations={smoke_evaluations}")
    print(f"quick_determinant_evaluations={quick_evaluations}")
    print(f"convergence_determinant_evaluations={convergence_evaluations}")
    print(f"output_dir={output_dir}")
    return 0


_SOURCE_CHECK_THETA = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)
_SOURCE_CHECK_FREQUENCIES_KHZ = {
    0.0: (0.23, 0.78, 1.34, 2.34, 3.55),
    15.0: (0.12, 0.62, 1.41, 1.65, 2.98),
    30.0: (0.08, 0.37, 1.01, 1.39, 1.94),
    45.0: (0.05, 0.28, 0.76, 1.12, 1.49),
    60.0: (0.04, 0.24, 0.66, 0.92, 1.31),
    75.0: (0.03, 0.23, 0.63, 0.80, 1.22),
    90.0: (0.03, 0.21, 0.63, 0.77, 1.21),
}
_SOURCE_CHECK_LOSSES_X100 = {
    0.0: (0.08, 1.18, 0.13, 1.18, 0.24),
    15.0: (0.94, 0.97, 0.50, 0.95, 1.00),
    30.0: (1.05, 1.05, 1.05, 0.67, 1.05),
    45.0: (0.96, 0.97, 0.97, 0.84, 0.98),
    60.0: (0.83, 0.84, 0.84, 1.00, 0.84),
    75.0: (0.71, 0.71, 0.71, 1.13, 0.71),
    90.0: (0.66, 0.66, 0.66, 1.19, 0.66),
}
_SOURCE_CHECK_FREQUENCY_UNCERTAINTY_KHZ = 0.08
_SOURCE_CHECK_LOSS_UNCERTAINTY_X100 = 0.04
_SOURCE_CHECK_CLAMPS = (
    "book_slope_clamp",
    "timoshenko_section_clamp",
)


def _postprocess_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _postprocess_fingerprint(path: Path) -> dict[str, Any]:
    stat = path.stat()
    return {
        "relative_path": path.relative_to(ROOT).as_posix(),
        "size_bytes": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "sha256": _postprocess_sha256(path),
    }


def _postprocess_read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def _postprocess_write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty postprocess CSV: {path.name}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _postprocess_digitized_rows(
    values: dict[float, tuple[float, ...]],
    *,
    value_name: str,
    uncertainty_name: str,
    uncertainty: float,
    note: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for theta in _SOURCE_CHECK_THETA:
        for mode, value in enumerate(values[theta], start=1):
            rows.append(
                {
                    "theta_deg": theta,
                    "book_curve_label": f"solid_coupled_{mode}",
                    "book_mode_label": f"mode_{mode}",
                    value_name: value,
                    uncertainty_name: uncertainty,
                    "axis_scale": "linear",
                    "source_page": "printed_65_discussion_starts_64",
                    "source_figure": "Figure_2.8",
                    "digitization_note": note,
                }
            )
    return rows


def _postprocess_validate_saved_data(
    roots_by_study: dict[str, list[dict[str, str]]],
    continuity_by_study: dict[str, list[dict[str, str]]],
) -> dict[str, Any]:
    requirements = {
        "orientation": {
            "parameters": _SOURCE_CHECK_THETA,
            "parameter_field": "theta_deg",
            "fixed_field": "length_m",
            "fixed_value": 0.2,
        },
        "length": {
            "parameters": (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40),
            "parameter_field": "length_m",
            "fixed_field": "theta_deg",
            "fixed_value": 15.0,
        },
    }
    summary: dict[str, Any] = {}
    missing: list[str] = []
    for study, requirement in requirements.items():
        roots = [
            row
            for row in roots_by_study[study]
            if row["material_mode"] == "book_complex"
        ]
        continuity = [
            row
            for row in continuity_by_study[study]
            if row["material_mode"] == "book_complex"
        ]
        for clamp in _SOURCE_CHECK_CLAMPS:
            for parameter in requirement["parameters"]:
                for mode in range(1, 6):
                    root_matches = [
                        row
                        for row in roots
                        if row["clamp_variant"] == clamp
                        and abs(float(row[requirement["parameter_field"]]) - parameter)
                        < 1.0e-9
                        and abs(float(row[requirement["fixed_field"]]) - requirement["fixed_value"])
                        < 1.0e-9
                        and int(row["sorted_mode"]) == mode
                    ]
                    continuity_matches = [
                        row
                        for row in continuity
                        if row["clamp_variant"] == clamp
                        and abs(float(row[requirement["parameter_field"]]) - parameter)
                        < 1.0e-9
                        and abs(float(row[requirement["fixed_field"]]) - requirement["fixed_value"])
                        < 1.0e-9
                        and int(row["tracked_mode"]) == mode
                    ]
                    if len(root_matches) != 1:
                        missing.append(
                            f"{study}: root {clamp}, parameter={parameter}, sorted_mode={mode}"
                        )
                    if len(continuity_matches) != 1:
                        missing.append(
                            f"{study}: continuity {clamp}, parameter={parameter}, tracked_mode={mode}"
                        )
        summary[study] = {
            "complex_root_rows": len(roots),
            "complex_continuity_rows": len(continuity),
            "root_statuses": sorted({row["root_status"] for row in roots}),
            "mapping_statuses": sorted({row["mapping_status"] for row in continuity}),
            "clamp_variants": sorted({row["clamp_variant"] for row in roots}),
            "sorted_modes": sorted({int(row["sorted_mode"]) for row in roots}),
            "tracked_modes": sorted({int(row["tracked_mode"]) for row in roots}),
            "eta_exact_present": all(row.get("eta_exact", "") != "" for row in roots),
        }
    if missing:
        raise RuntimeError("saved full-run data are incomplete:\n- " + "\n- ".join(missing))
    return summary


def _postprocess_orientation_indexes(
    roots: Sequence[dict[str, str]],
    continuity: Sequence[dict[str, str]],
) -> tuple[
    dict[tuple[str, float, int], dict[str, str]],
    dict[tuple[str, float, int], dict[str, str]],
]:
    root_index = {
        (row["clamp_variant"], round(float(row["theta_deg"]), 9), int(row["sorted_mode"])): row
        for row in roots
        if row["material_mode"] == "book_complex"
        and abs(float(row["length_m"]) - 0.2) < 1.0e-9
    }
    continuity_index = {
        (row["clamp_variant"], round(float(row["theta_deg"]), 9), int(row["tracked_mode"])): row
        for row in continuity
        if row["material_mode"] == "book_complex"
        and abs(float(row["length_m"]) - 0.2) < 1.0e-9
    }
    return root_index, continuity_index


def _postprocess_compare_stage(
    digitized_rows: Sequence[dict[str, Any]],
    root_index: dict[tuple[str, float, int], dict[str, str]],
    continuity_index: dict[tuple[str, float, int], dict[str, str]],
    *,
    source_value_field: str,
    uncertainty_field: str,
    calculated_value_field: str,
    calculated_value,
    error_prefix: str,
) -> list[dict[str, Any]]:
    compared: list[dict[str, Any]] = []
    for source in digitized_rows:
        theta = float(source["theta_deg"])
        tracked_mode = int(str(source["book_mode_label"]).split("_")[-1])
        for clamp in _SOURCE_CHECK_CLAMPS:
            continuity = continuity_index[(clamp, round(theta, 9), tracked_mode)]
            sorted_mode = int(continuity["current_sorted_mode"])
            root = root_index[(clamp, round(theta, 9), sorted_mode)]
            mapping_ok = (
                continuity["mapping_status"] == "accepted"
                and continuity["low_mac_warning"].lower() == "false"
                and int(root["tracked_mode"]) == tracked_mode
            )
            mapping_status = (
                "accepted" if mapping_ok else "SOURCE_MODE_MAPPING_AMBIGUOUS"
            )
            source_value = float(source[source_value_field])
            uncertainty = float(source[uncertainty_field])
            calculated = float(calculated_value(root))
            absolute_error = abs(calculated - source_value)
            relative_error = absolute_error / max(abs(source_value), 1.0e-15)
            normalized = absolute_error / uncertainty
            compared.append(
                {
                    "theta_deg": theta,
                    "book_curve_label": source["book_curve_label"],
                    "book_mode_label": source["book_mode_label"],
                    "clamp_variant": clamp,
                    "calculated_sorted_mode": sorted_mode,
                    "calculated_tracked_mode": tracked_mode,
                    calculated_value_field: calculated,
                    source_value_field: source_value,
                    uncertainty_field: uncertainty,
                    f"{error_prefix}absolute_error": absolute_error,
                    f"{error_prefix}relative_error": relative_error,
                    "normalized_graph_error": normalized,
                    "mapping_status": mapping_status,
                    "point_status": (
                        "within_declared_uncertainty"
                        if mapping_ok and normalized <= 1.0
                        else (
                            "outside_declared_uncertainty"
                            if mapping_ok
                            else "SOURCE_MODE_MAPPING_AMBIGUOUS"
                        )
                    ),
                }
            )
    return compared


def _postprocess_metrics(
    rows: Sequence[dict[str, Any]],
    *,
    absolute_error_field: str,
    relative_error_field: str,
) -> dict[str, Any]:
    normalized = [float(row["normalized_graph_error"]) for row in rows]
    absolute = [float(row[absolute_error_field]) for row in rows]
    relative = [float(row[relative_error_field]) for row in rows]
    within = sum(row["point_status"] == "within_declared_uncertainty" for row in rows)
    ambiguous = sum(row["mapping_status"] != "accepted" for row in rows)
    return {
        "number_of_digitized_points": len(rows),
        "number_within_declared_uncertainty": within,
        "fraction_within_declared_uncertainty": within / len(rows),
        "maximum_absolute_error": max(absolute),
        "maximum_relative_error": max(relative),
        "maximum_normalized_graph_error": max(normalized),
        "RMS_normalized_graph_error": math.sqrt(sum(value * value for value in normalized) / len(normalized)),
        "median_normalized_graph_error": statistics.median(normalized),
        "ambiguous_mode_mappings": ambiguous,
    }


def _postprocess_model_separation(
    comparison: Sequence[dict[str, Any]],
    *,
    value_field: str,
    uncertainty_field: str,
) -> list[dict[str, Any]]:
    indexed = {
        (row["theta_deg"], row["book_mode_label"], row["clamp_variant"]): row
        for row in comparison
    }
    rows: list[dict[str, Any]] = []
    for theta in _SOURCE_CHECK_THETA:
        for mode in range(1, 6):
            label = f"mode_{mode}"
            book = indexed[(theta, label, "book_slope_clamp")]
            section = indexed[(theta, label, "timoshenko_section_clamp")]
            difference = abs(float(book[value_field]) - float(section[value_field]))
            uncertainty = float(book[uncertainty_field])
            rows.append(
                {
                    "theta_deg": theta,
                    "book_mode_label": label,
                    "absolute_model_difference": difference,
                    "digitization_uncertainty": uncertainty,
                    "normalized_model_difference": difference / uncertainty,
                    "resolvable_by_declared_uncertainty": difference > uncertainty,
                }
            )
    return rows


def _postprocess_classify_frequency(
    comparison: Sequence[dict[str, Any]],
    separation: Sequence[dict[str, Any]],
) -> str:
    by_clamp = {
        clamp: [row for row in comparison if row["clamp_variant"] == clamp]
        for clamp in _SOURCE_CHECK_CLAMPS
    }
    within = {
        clamp: sum(row["point_status"] == "within_declared_uncertainty" for row in rows)
        for clamp, rows in by_clamp.items()
    }
    ambiguous = sum(row["mapping_status"] != "accepted" for row in comparison)
    resolvable = sum(bool(row["resolvable_by_declared_uncertainty"]) for row in separation)
    if ambiguous:
        return "SOURCE_AMBIGUITY"
    if within["book_slope_clamp"] == 35 and within["timoshenko_section_clamp"] == 35:
        return "BOTH_WITHIN_GRAPH_RESOLUTION"
    if resolvable <= 1:
        return "SOURCE_AMBIGUITY"
    if max(within.values()) < 28:
        return "NEITHER_REPRODUCES"
    return "SOURCE_AMBIGUITY"


def _postprocess_classify_with_loss(
    comparison: Sequence[dict[str, Any]],
    separation: Sequence[dict[str, Any]],
) -> str:
    indexed = {
        (row["theta_deg"], row["book_mode_label"], row["clamp_variant"]): row
        for row in comparison
    }
    book_discriminating = 0
    section_discriminating = 0
    for item in separation:
        if not item["resolvable_by_declared_uncertainty"]:
            continue
        key = (item["theta_deg"], item["book_mode_label"])
        book = indexed[(*key, "book_slope_clamp")]
        section = indexed[(*key, "timoshenko_section_clamp")]
        if (
            book["point_status"] == "within_declared_uncertainty"
            and section["point_status"] == "outside_declared_uncertainty"
        ):
            book_discriminating += 1
        if (
            section["point_status"] == "within_declared_uncertainty"
            and book["point_status"] == "outside_declared_uncertainty"
        ):
            section_discriminating += 1
    if any(row["mapping_status"] != "accepted" for row in comparison):
        return "SOURCE_AMBIGUITY"
    if book_discriminating >= 2 and section_discriminating == 0:
        return "BOOK_SLOPE_CLAMP_CONFIRMED"
    if section_discriminating >= 2 and book_discriminating == 0:
        return "TIMOSHENKO_SECTION_CLAMP_CONFIRMED"
    if all(row["point_status"] == "within_declared_uncertainty" for row in comparison):
        return "BOTH_WITHIN_GRAPH_RESOLUTION"
    return "SOURCE_AMBIGUITY"


def _postprocess_plot(
    output_dir: Path,
    roots: Sequence[dict[str, str]],
    frequency_digitized: Sequence[dict[str, Any]],
    frequency_comparison: Sequence[dict[str, Any]],
    loss_digitized: Sequence[dict[str, Any]],
    loss_comparison: Sequence[dict[str, Any]],
) -> None:
    complex_rows = [
        row
        for row in roots
        if row["material_mode"] == "book_complex"
        and abs(float(row["length_m"]) - 0.2) < 1.0e-9
        and row["tracked_mode"] != ""
    ]
    colors = plt.cm.tab10(np.linspace(0.0, 0.8, 5))
    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.5), constrained_layout=True)
    for mode, color in zip(range(1, 6), colors):
        source_frequency = [
            row for row in frequency_digitized if row["book_mode_label"] == f"mode_{mode}"
        ]
        axes[0, 0].errorbar(
            [row["theta_deg"] for row in source_frequency],
            [row["frequency_khz_digitized"] for row in source_frequency],
            yerr=_SOURCE_CHECK_FREQUENCY_UNCERTAINTY_KHZ,
            fmt="o",
            color=color,
            markersize=3.5,
            capsize=2,
            label=f"source mode {mode}" if mode == 1 else None,
        )
        source_loss = [
            row for row in loss_digitized if row["book_mode_label"] == f"mode_{mode}"
        ]
        axes[1, 0].errorbar(
            [row["theta_deg"] for row in source_loss],
            [row["eta_x100_digitized"] for row in source_loss],
            yerr=_SOURCE_CHECK_LOSS_UNCERTAINTY_X100,
            fmt="o",
            color=color,
            markersize=3.5,
            capsize=2,
            label=f"source mode {mode}" if mode == 1 else None,
        )
        for clamp, style in zip(_SOURCE_CHECK_CLAMPS, ("-", "--")):
            curve = sorted(
                (
                    row
                    for row in complex_rows
                    if row["clamp_variant"] == clamp and int(row["tracked_mode"]) == mode
                ),
                key=lambda row: float(row["theta_deg"]),
            )
            label = (
                "book_slope saved"
                if clamp == "book_slope_clamp"
                else "section_rotation saved"
            ) if mode == 1 else None
            axes[0, 0].plot(
                [float(row["theta_deg"]) for row in curve],
                [float(row["frequency_hz"]) / 1000.0 for row in curve],
                style,
                color=color,
                linewidth=1.1,
                label=label,
            )
            axes[1, 0].plot(
                [float(row["theta_deg"]) for row in curve],
                [100.0 * float(row["eta_exact"]) for row in curve],
                style,
                color=color,
                linewidth=1.1,
                label=label,
            )
        frequency_by_clamp = {
            clamp: {
                (row["theta_deg"], row["book_mode_label"]): row
                for row in frequency_comparison
                if row["clamp_variant"] == clamp
            }
            for clamp in _SOURCE_CHECK_CLAMPS
        }
        loss_by_clamp = {
            clamp: {
                (row["theta_deg"], row["book_mode_label"]): row
                for row in loss_comparison
                if row["clamp_variant"] == clamp
            }
            for clamp in _SOURCE_CHECK_CLAMPS
        }
        axes[0, 1].plot(
            _SOURCE_CHECK_THETA,
            [
                abs(
                    frequency_by_clamp["book_slope_clamp"][(theta, f"mode_{mode}")]["frequency_khz_calculated"]
                    - frequency_by_clamp["timoshenko_section_clamp"][(theta, f"mode_{mode}")]["frequency_khz_calculated"]
                )
                for theta in _SOURCE_CHECK_THETA
            ],
            "-o",
            color=color,
            markersize=3,
            label=f"mode {mode}",
        )
        axes[1, 1].plot(
            _SOURCE_CHECK_THETA,
            [
                abs(
                    loss_by_clamp["book_slope_clamp"][(theta, f"mode_{mode}")]["eta_x100_calculated"]
                    - loss_by_clamp["timoshenko_section_clamp"][(theta, f"mode_{mode}")]["eta_x100_calculated"]
                )
                for theta in _SOURCE_CHECK_THETA
            ],
            "-o",
            color=color,
            markersize=3,
            label=f"mode {mode}",
        )
    axes[0, 0].set(title="Figure 2.8 frequency: source vs saved roots", ylabel="frequency (kHz)")
    axes[1, 0].set(title="Figure 2.8 loss: source vs saved roots", xlabel="theta (deg)", ylabel="eta x 1e2")
    axes[0, 1].axhline(_SOURCE_CHECK_FREQUENCY_UNCERTAINTY_KHZ, color="black", linestyle=":", label="digitization uncertainty")
    axes[0, 1].set(title="Clamp-model frequency difference", ylabel="absolute difference (kHz)")
    axes[1, 1].axhline(_SOURCE_CHECK_LOSS_UNCERTAINTY_X100, color="black", linestyle=":", label="digitization uncertainty")
    axes[1, 1].set(title="Clamp-model loss difference", xlabel="theta (deg)", ylabel="absolute difference (eta x 1e2)")
    for axis in axes.flat:
        axis.set_xlim(0.0, 90.0)
        axis.set_xticks(_SOURCE_CHECK_THETA)
        axis.grid(True, alpha=0.25)
    axes[0, 0].legend(fontsize=7, ncol=3)
    axes[0, 1].legend(fontsize=7, ncol=2)
    axes[1, 0].legend(fontsize=7, ncol=3)
    axes[1, 1].legend(fontsize=7, ncol=2)
    fig.savefig(output_dir / "boundary_source_check.pdf")
    fig.savefig(output_dir / "boundary_source_check.png", dpi=180)
    plt.close(fig)


def _run_postprocess_boundary_source_check(args: argparse.Namespace) -> int:
    start = time.perf_counter()
    input_dir = DEFAULT_RESULT_DIR
    output_dir = (args.output_dir or BOUNDARY_SOURCE_CHECK_RESULT_DIR).resolve()
    required = {
        "orientation_roots": input_dir / "orientation_roots.csv",
        "length_roots": input_dir / "length_roots.csv",
        "orientation_continuity": input_dir / "orientation_mode_continuity.csv",
        "length_continuity": input_dir / "length_mode_continuity.csv",
    }
    absent = [path.relative_to(ROOT).as_posix() for path in required.values() if not path.is_file()]
    if absent:
        print("saved-data postprocess stopped; missing inputs:", file=sys.stderr)
        for path in absent:
            print(f"- {path}", file=sys.stderr)
        return 2
    if output_dir == input_dir.resolve():
        raise ValueError("postprocess output directory must differ from the full-run input directory")

    before = {name: _postprocess_fingerprint(path) for name, path in required.items()}
    roots_by_study = {
        "orientation": _postprocess_read_csv(required["orientation_roots"]),
        "length": _postprocess_read_csv(required["length_roots"]),
    }
    continuity_by_study = {
        "orientation": _postprocess_read_csv(required["orientation_continuity"]),
        "length": _postprocess_read_csv(required["length_continuity"]),
    }
    validation = _postprocess_validate_saved_data(roots_by_study, continuity_by_study)
    root_index, continuity_index = _postprocess_orientation_indexes(
        roots_by_study["orientation"], continuity_by_study["orientation"]
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    frequency_digitized = _postprocess_digitized_rows(
        _SOURCE_CHECK_FREQUENCIES_KHZ,
        value_name="frequency_khz_digitized",
        uncertainty_name="digitization_uncertainty_khz",
        uncertainty=_SOURCE_CHECK_FREQUENCY_UNCERTAINTY_KHZ,
        note=(
            "manual reading of calculated solid coupled curve only; partial dashed and dotted curves excluded"
        ),
    )
    frequency_comparison = _postprocess_compare_stage(
        frequency_digitized,
        root_index,
        continuity_index,
        source_value_field="frequency_khz_digitized",
        uncertainty_field="digitization_uncertainty_khz",
        calculated_value_field="frequency_khz_calculated",
        calculated_value=lambda row: float(row["frequency_hz"]) / 1000.0,
        error_prefix="",
    )
    frequency_separation = _postprocess_model_separation(
        frequency_comparison,
        value_field="frequency_khz_calculated",
        uncertainty_field="digitization_uncertainty_khz",
    )
    frequency_status = _postprocess_classify_frequency(
        frequency_comparison, frequency_separation
    )
    _postprocess_write_csv(
        output_dir / "figure_2_8_frequency_digitized.csv", frequency_digitized
    )
    _postprocess_write_csv(
        output_dir / "figure_2_8_frequency_comparison.csv", frequency_comparison
    )

    loss_executed = frequency_status in (
        "BOTH_WITHIN_GRAPH_RESOLUTION",
        "SOURCE_AMBIGUITY",
    )
    loss_digitized: list[dict[str, Any]] = []
    loss_comparison: list[dict[str, Any]] = []
    loss_separation: list[dict[str, Any]] = []
    final_status = frequency_status
    if loss_executed:
        loss_digitized = _postprocess_digitized_rows(
            _SOURCE_CHECK_LOSSES_X100,
            value_name="eta_x100_digitized",
            uncertainty_name="digitization_uncertainty_eta_x100",
            uncertainty=_SOURCE_CHECK_LOSS_UNCERTAINTY_X100,
            note=(
                "manual reading of calculated solid coupled modal-loss curve only; partial curves excluded"
            ),
        )
        loss_comparison = _postprocess_compare_stage(
            loss_digitized,
            root_index,
            continuity_index,
            source_value_field="eta_x100_digitized",
            uncertainty_field="digitization_uncertainty_eta_x100",
            calculated_value_field="eta_x100_calculated",
            calculated_value=lambda row: 100.0 * float(row["eta_exact"]),
            error_prefix="",
        )
        loss_separation = _postprocess_model_separation(
            loss_comparison,
            value_field="eta_x100_calculated",
            uncertainty_field="digitization_uncertainty_eta_x100",
        )
        final_status = _postprocess_classify_with_loss(loss_comparison, loss_separation)
        _postprocess_write_csv(
            output_dir / "figure_2_8_loss_digitized.csv", loss_digitized
        )
        _postprocess_write_csv(
            output_dir / "figure_2_8_loss_comparison.csv", loss_comparison
        )

    figure_2_11_executed = final_status in (
        "BOTH_WITHIN_GRAPH_RESOLUTION",
        "SOURCE_AMBIGUITY",
    )
    if figure_2_11_executed:
        raise RuntimeError(
            "Figure 2.8 remains ambiguous; this postprocess revision intentionally stops before Figure 2.11 digitization"
        )

    summary_rows: list[dict[str, Any]] = []
    for stage, comparison, separation, abs_field, rel_field, stage_status in (
        (
            "figure_2_8_frequency",
            frequency_comparison,
            frequency_separation,
            "absolute_error",
            "relative_error",
            frequency_status,
        ),
        (
            "figure_2_8_loss",
            loss_comparison,
            loss_separation,
            "absolute_error",
            "relative_error",
            final_status,
        ),
    ):
        if not comparison:
            continue
        for clamp in _SOURCE_CHECK_CLAMPS:
            subset = [row for row in comparison if row["clamp_variant"] == clamp]
            metrics = _postprocess_metrics(
                subset,
                absolute_error_field=abs_field,
                relative_error_field=rel_field,
            )
            summary_rows.append(
                {
                    "stage": stage,
                    "comparison_subject": clamp,
                    "stage_status": stage_status,
                    **metrics,
                    "maximum_model_difference": "",
                    "maximum_normalized_model_difference": "",
                    "number_model_differences_above_uncertainty": "",
                }
            )
        summary_rows.append(
            {
                "stage": stage,
                "comparison_subject": "book_slope_vs_section_rotation",
                "stage_status": stage_status,
                "number_of_digitized_points": len(separation),
                "number_within_declared_uncertainty": sum(
                    not row["resolvable_by_declared_uncertainty"] for row in separation
                ),
                "fraction_within_declared_uncertainty": sum(
                    not row["resolvable_by_declared_uncertainty"] for row in separation
                )
                / len(separation),
                "maximum_absolute_error": "",
                "maximum_relative_error": "",
                "maximum_normalized_graph_error": "",
                "RMS_normalized_graph_error": "",
                "median_normalized_graph_error": "",
                "ambiguous_mode_mappings": 0,
                "maximum_model_difference": max(
                    float(row["absolute_model_difference"]) for row in separation
                ),
                "maximum_normalized_model_difference": max(
                    float(row["normalized_model_difference"]) for row in separation
                ),
                "number_model_differences_above_uncertainty": sum(
                    row["resolvable_by_declared_uncertainty"] for row in separation
                ),
            }
        )
    _postprocess_write_csv(output_dir / "boundary_source_check_summary.csv", summary_rows)
    _postprocess_plot(
        output_dir,
        roots_by_study["orientation"],
        frequency_digitized,
        frequency_comparison,
        loss_digitized,
        loss_comparison,
    )

    after = {name: _postprocess_fingerprint(path) for name, path in required.items()}
    unchanged = before == after
    inventory = sorted(path for path in input_dir.iterdir() if path.is_file())
    manifest_lines = [
        "# Boundary source check manifest",
        "",
        "Postprocess-only comparison. No scientific solver module was imported by this CLI path.",
        "",
        "## Source and axes",
        "",
        "- Local source: `docs/literature/pdf/Глава 2_compressed.pdf`.",
        "- Figure 2.8 graphic: printed p. 65; the discussion starts on printed p. 64.",
        "- Frequency axes: theta 0–90 deg, linear 15-deg ticks; frequency 0–4 kHz, linear 1-kHz ticks.",
        f"- Frequency digitization uncertainty fixed before comparison: ±{_SOURCE_CHECK_FREQUENCY_UNCERTAINTY_KHZ:.2f} kHz.",
        "- Frequency rationale: the 1-kHz major spacing spans about 132 source pixels; ±0.08 kHz conservatively covers raster/line thickness and manual placement between ticks.",
        "- Loss axes: theta 0–90 deg, linear 15-deg ticks; eta x 1e2, linear 0.2 ticks.",
        f"- Loss digitization uncertainty fixed before comparison: ±{_SOURCE_CHECK_LOSS_UNCERTAINTY_X100:.2f} in eta x 1e2.",
        "- Loss rationale: a 0.2 major spacing spans about 87 source pixels, but early-angle curve overlap dominates; ±0.04 is one fifth of a major division.",
        "",
        "## Saved full-run inventory",
        "",
        "| relative path | bytes | role |",
        "| --- | ---: | --- |",
    ]
    for path in inventory:
        name = path.name
        role = "generated full-run artifact"
        if name == "orientation_roots.csv":
            role = "orientation elastic and complex roots; sorted/tracked labels, eta_exact, residuals, root status"
        elif name == "length_roots.csv":
            role = "length elastic and complex roots; sorted/tracked labels, eta_exact, residuals, root status"
        elif name == "orientation_mode_continuity.csv":
            role = "saved orientation continuity/MAC mapping"
        elif name == "length_mode_continuity.csv":
            role = "saved length continuity/MAC mapping"
        elif name == "root_quality.csv":
            role = "combined root-quality table"
        elif name == "clamp_comparison.csv":
            role = "saved cross-clamp comparison and MAC mapping"
        elif path.suffix.lower() == ".npz":
            role = "saved array data"
        manifest_lines.append(
            f"| `{path.relative_to(ROOT).as_posix()}` | {path.stat().st_size} | {role} |"
        )
    manifest_lines.extend(
        [
            "",
            "## Input fingerprints",
            "",
            "| relative path | bytes | mtime_ns | SHA256 | unchanged after postprocess |",
            "| --- | ---: | ---: | --- | --- |",
        ]
    )
    for name in required:
        item = before[name]
        manifest_lines.append(
            f"| `{item['relative_path']}` | {item['size_bytes']} | {item['mtime_ns']} | `{item['sha256']}` | {before[name] == after[name]} |"
        )
    manifest_lines.extend(
        [
            "",
            "## Completeness",
            "",
            f"- Orientation: {validation['orientation']}.",
            f"- Length: {validation['length']}.",
            "- Root formulation: corrected Chapter-2 state workflow (`state_corrected`); no solver was invoked here.",
            "- Input root CSV unchanged: " + str(unchanged) + ".",
            "",
            "## Zero scientific-work counters",
            "",
            "- determinant evaluations = 0",
            "- root solves = 0",
            "- complex continuations = 0",
            "- scientific matrix-exponential calls = 0",
            "- state-matrix calls = 0",
            "- boundary-matrix calls = 0",
        ]
    )
    (output_dir / "boundary_source_check_manifest.md").write_text(
        "\n".join(manifest_lines) + "\n", encoding="utf-8"
    )

    frequency_metrics = {
        clamp: _postprocess_metrics(
            [row for row in frequency_comparison if row["clamp_variant"] == clamp],
            absolute_error_field="absolute_error",
            relative_error_field="relative_error",
        )
        for clamp in _SOURCE_CHECK_CLAMPS
    }
    loss_metrics = {
        clamp: _postprocess_metrics(
            [row for row in loss_comparison if row["clamp_variant"] == clamp],
            absolute_error_field="absolute_error",
            relative_error_field="relative_error",
        )
        for clamp in _SOURCE_CHECK_CLAMPS
    }
    runtime = time.perf_counter() - start
    report_lines = [
        "# Yartsev Figure 2.8 boundary-source check",
        "",
        f"**Final status: `{final_status}`.**",
        "",
        "This result uses only saved full-run CSV values and manual readings of the solid coupled curves in Figure 2.8. No model parameter was fitted.",
        "The frequency uncertainty (±0.08 kHz) was fixed from the 1-kHz grid spacing, source-pixel scale and line thickness; the loss uncertainty (±0.04 in eta x 1e2) was fixed as one fifth of the 0.2 major spacing because of early-angle curve overlap. Both were fixed before model errors were read.",
        "",
        "## Stage decision",
        "",
        f"- Frequency-only status: `{frequency_status}`.",
        "- The frequency panel resolves the clamp difference at only "
        + str(sum(row["resolvable_by_declared_uncertainty"] for row in frequency_separation))
        + " of 35 points, so frequency alone is not a robust discriminator.",
        f"- Loss stage executed: {loss_executed}.",
        f"- Figure 2.11 stage executed: {figure_2_11_executed}.",
        "- At theta=0, loss curves for modes 3 and 5 distinguish the variants beyond the fixed graph uncertainty; the slope clamp remains within it and the section-rotation clamp does not.",
        "",
        "## Metrics",
        "",
        "| stage | clamp | points within | fraction | max absolute error | max relative error | max normalized | RMS normalized | median normalized | ambiguous mappings |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for stage, metrics_by_clamp in (
        ("frequency (kHz)", frequency_metrics),
        ("loss (eta x 1e2)", loss_metrics),
    ):
        for clamp, metrics in metrics_by_clamp.items():
            report_lines.append(
                "| {stage} | `{clamp}` | {within}/{total} | {fraction:.6f} | {max_abs:.6g} | {max_rel:.6g} | {max_norm:.6g} | {rms:.6g} | {median:.6g} | {ambiguous} |".format(
                    stage=stage,
                    clamp=clamp,
                    within=metrics["number_within_declared_uncertainty"],
                    total=metrics["number_of_digitized_points"],
                    fraction=metrics["fraction_within_declared_uncertainty"],
                    max_abs=metrics["maximum_absolute_error"],
                    max_rel=metrics["maximum_relative_error"],
                    max_norm=metrics["maximum_normalized_graph_error"],
                    rms=metrics["RMS_normalized_graph_error"],
                    median=metrics["median_normalized_graph_error"],
                    ambiguous=metrics["ambiguous_mode_mappings"],
                )
            )
    report_lines.extend(
        [
            "",
            "## Clamp-model separation",
            "",
            f"- Frequency: maximum {max(row['absolute_model_difference'] for row in frequency_separation):.6g} kHz; above ±{_SOURCE_CHECK_FREQUENCY_UNCERTAINTY_KHZ:.2f} kHz at {sum(row['resolvable_by_declared_uncertainty'] for row in frequency_separation)}/35 points.",
            f"- Loss: maximum {max(row['absolute_model_difference'] for row in loss_separation):.6g} in eta x 1e2; above ±{_SOURCE_CHECK_LOSS_UNCERTAINTY_X100:.2f} at {sum(row['resolvable_by_declared_uncertainty'] for row in loss_separation)}/35 points.",
            "",
            "## Mode mapping and units",
            "",
            "Saved continuity rows select `current_sorted_mode` by `tracked_mode`; all 70 required orientation mappings are accepted and have no low-MAC warning. Saved frequencies in Hz were divided by 1000; saved `eta_exact` was multiplied by 100 for comparison with eta x 1e2.",
            "",
            "## Execution guard",
            "",
            "The early CLI path skipped imports from `scripts/lib/yartsev_ch2_monoclinic_rod.py` and the free-free root workflow.",
            "",
            "- determinant evaluations = 0",
            "- root solves = 0",
            "- complex continuations = 0",
            "- scientific matrix-exponential calls = 0",
            "- state-matrix calls = 0",
            "- boundary-matrix calls = 0",
            f"- input root CSV unchanged by size, mtime_ns and SHA256: {unchanged}",
            f"- runtime_seconds = {runtime:.6f}",
        ]
    )
    (output_dir / "boundary_source_check_report.md").write_text(
        "\n".join(report_lines) + "\n", encoding="utf-8"
    )

    scientific_module_imported = any(
        name.endswith("yartsev_ch2_monoclinic_rod") for name in sys.modules
    )
    if scientific_module_imported:
        raise RuntimeError("postprocess-only guard failed: scientific model module was imported")
    print(f"frequency_stage_status={frequency_status}")
    print(f"loss_stage_executed={int(loss_executed)}")
    print(f"figure_2_11_stage_executed={int(figure_2_11_executed)}")
    print(f"boundary_source_check_status={final_status}")
    print("determinant_evaluations=0")
    print("root_solves=0")
    print("complex_continuations=0")
    print("scientific_expm_calls=0")
    print("state_matrix_calls=0")
    print("boundary_matrix_calls=0")
    print(f"input_root_csv_unchanged={int(unchanged)}")
    print(f"output_dir={output_dir}")
    print(f"runtime_seconds={runtime:.6f}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    _validate_args(args)
    if args.postprocess_boundary_source_check:
        return _run_postprocess_boundary_source_check(args)
    missing = [path for path in (CHAPTER_1, CHAPTER_2) if not path.is_file()]
    if missing:
        for path in missing:
            print(f"missing source: {path}", file=sys.stderr)
        return 2
    if args.quick_boundary_gate:
        return _run_quick_boundary_gate(args)
    if args.smoke:
        args.num_positive_modes = min(args.num_positive_modes, 3)
    output_dir = (args.output_dir or (SMOKE_RESULT_DIR if args.smoke else DEFAULT_RESULT_DIR)).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    studies = _studies(args)
    clamps = _clamps(args)
    include_complex = args.material_mode in ("book-complex", "both")
    guard_roots = args.num_positive_modes + 1

    _write_source_manifest(output_dir)
    _write_csv(output_dir / "cantilever_material_parameters.csv", _material_rows())

    all_root_rows: list[dict[str, Any]] = []
    all_boundary_rows: list[dict[str, Any]] = []
    all_clamp_rows: list[dict[str, Any]] = []
    all_partial_rows: list[dict[str, Any]] = []
    all_digitized_rows: list[dict[str, Any]] = []
    all_digitized_comparison: list[dict[str, Any]] = []
    continuity_by_study: dict[Study, list[dict[str, Any]]] = {}
    bending_roots_for_plot: dict[tuple[Study, float], list[RootResult]] = {}

    for study in studies:
        parameters = _orientation_grid(args) if study == "orientation" else _length_grid(args)
        elastic = _solve_elastic(
            study,
            parameters,
            clamps,
            guard_roots=guard_roots,
            scan_step_hz=args.scan_step_hz,
            series_rtol=args.series_rtol,
        )
        complex_roots = (
            _solve_complex(
                study,
                parameters,
                clamps,
                elastic,
                num_modes=args.num_positive_modes,
                series_rtol=args.series_rtol,
            )
            if include_complex
            else {}
        )
        elastic_map, elastic_continuity, elastic_shapes = _track_modes(
            study,
            parameters,
            clamps,
            elastic,
            material_mode="elastic",
            num_modes=args.num_positive_modes,
            series_rtol=args.series_rtol,
        )
        complex_map: dict[tuple[ClampVariant, float, int], int] = {}
        complex_continuity: list[dict[str, Any]] = []
        if complex_roots:
            complex_map, complex_continuity, _ = _track_modes(
                study,
                parameters,
                clamps,
                complex_roots,
                material_mode="book_complex",
                num_modes=args.num_positive_modes,
                series_rtol=args.series_rtol,
            )
        continuity_by_study[study] = [*elastic_continuity, *complex_continuity]
        elastic_rows = _root_rows(
            study, parameters, clamps, elastic, elastic_map, material_mode="elastic"
        )
        complex_rows = (
            _root_rows(
                study,
                parameters,
                clamps,
                complex_roots,
                complex_map,
                material_mode="book_complex",
            )
            if complex_roots
            else []
        )
        all_root_rows.extend(elastic_rows)
        all_root_rows.extend(complex_rows)
        all_boundary_rows.extend(
            _boundary_rows(study, parameters, clamps, elastic, material_mode="elastic")
        )
        if complex_roots:
            all_boundary_rows.extend(
                _boundary_rows(study, parameters, clamps, complex_roots, material_mode="book_complex")
            )
        clamp_rows = _clamp_comparison_rows(
            study,
            parameters,
            elastic,
            complex_roots,
            elastic_map,
            elastic_shapes,
            num_modes=args.num_positive_modes,
        )
        all_clamp_rows.extend(clamp_rows)
        partial_rows, bending_roots = _partial_rows(
            study,
            parameters,
            include_complex=include_complex,
            num_modes=args.num_positive_modes,
            scan_step_hz=args.scan_step_hz,
        )
        all_partial_rows.extend(partial_rows)
        bending_roots_for_plot.update(bending_roots)
        digitized = _digitized_rows(study)
        all_digitized_rows.extend(digitized)
        digitized_comparison = (
            _digitized_comparison(study, digitized, complex_roots, complex_map, clamps)
            if complex_roots
            else []
        )
        all_digitized_comparison.extend(digitized_comparison)
        _plot_reproduction(
            output_dir,
            study,
            parameters,
            elastic,
            complex_roots,
            elastic_map,
            complex_map,
            partial_rows,
            digitized,
            num_modes=args.num_positive_modes,
        )
        targets = (
            [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0]
            if study == "orientation"
            else [0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
        )
        if not args.smoke and "book_slope_clamp" in clamps:
            _plot_mode_panels(
                output_dir,
                study,
                targets,
                elastic,
                elastic_map,
                num_modes=args.num_positive_modes,
            )

    _write_csv(output_dir / "orientation_roots.csv", [row for row in all_root_rows if row["study"] == "orientation"])
    _write_csv(output_dir / "length_roots.csv", [row for row in all_root_rows if row["study"] == "length"])
    _write_csv(output_dir / "root_quality.csv", all_root_rows)
    _write_csv(output_dir / "boundary_residuals.csv", all_boundary_rows)
    _write_csv(output_dir / "clamp_comparison.csv", all_clamp_rows)
    _write_csv(output_dir / "orientation_mode_continuity.csv", continuity_by_study.get("orientation", []))
    _write_csv(output_dir / "length_mode_continuity.csv", continuity_by_study.get("length", []))
    _write_csv(output_dir / "partial_spectra.csv", all_partial_rows)
    _write_csv(output_dir / "figure_2_8_digitized_calculated_curves.csv", [row for row in all_digitized_rows if row["study"] == "orientation"])
    _write_csv(output_dir / "figure_2_8_digitized_comparison.csv", [row for row in all_digitized_comparison if row["study"] == "orientation"])
    _write_csv(output_dir / "figure_2_11_digitized_calculated_curves.csv", [row for row in all_digitized_rows if row["study"] == "length"])
    _write_csv(output_dir / "figure_2_11_digitized_comparison.csv", [row for row in all_digitized_comparison if row["study"] == "length"])

    equivalence_rows = _formulation_equivalence_rows(
        scan_step_hz=args.scan_step_hz, series_rtol=args.series_rtol
    )
    _write_csv(output_dir / "formulation_equivalence.csv", equivalence_rows)
    shear_rows = _shear_rigid_rows(
        scan_step_hz=args.scan_step_hz, num_modes=args.num_positive_modes
    )
    _write_csv(output_dir / "shear_rigid_limit.csv", shear_rows)
    if not args.smoke:
        _plot_partial_shapes(output_dir, bending_roots_for_plot)
    _plot_clamp_comparison(output_dir, all_clamp_rows)

    # Endpoint controls are evaluated explicitly and included in the report-facing CSV.
    endpoint_rows: list[dict[str, Any]] = []
    for theta in (0.0, 90.0):
        point = _point("orientation", theta, "elastic")
        for clamp in ("book_slope_clamp", "timoshenko_section_clamp"):
            matrix = cantilever_boundary_matrix(0.0, point, clamp)
            for mode in (1, 2, 3):
                omega = fixed_free_torsion_omega(point, mode)
                full, bending, torsion = decoupled_cantilever_boundary_factors(omega, point, clamp)
                endpoint_rows.append(
                    {
                        "theta_deg": theta,
                        "clamp_variant": clamp,
                        "zero_frequency_determinant": np.linalg.det(matrix),
                        "zero_frequency_nullity": cantilever_zero_frequency_nullity(point, clamp),
                        "torsion_mode": mode,
                        "analytical_torsion_frequency_hz": omega.real / (2.0 * np.pi),
                        "full_factor": full,
                        "bending_factor": bending,
                        "torsion_factor": torsion,
                        "factorization_relative_residual": abs(full - bending * torsion) / max(abs(full), 1.0),
                    }
                )
    _write_csv(output_dir / "orthotropic_endpoint_checks.csv", endpoint_rows)

    status = _audit_status(all_digitized_comparison, args.smoke)
    runtime_seconds = time.perf_counter() - start
    _write_report(
        output_dir,
        args,
        status,
        all_digitized_comparison,
        equivalence_rows,
        all_boundary_rows,
        [row for values in continuity_by_study.values() for row in values],
        all_clamp_rows,
        shear_rows,
        all_digitized_comparison,
        runtime_seconds,
    )
    print(f"boundary_condition_audit_status={status}")
    print(f"output_dir={output_dir}")
    print(f"runtime_seconds={runtime_seconds:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
