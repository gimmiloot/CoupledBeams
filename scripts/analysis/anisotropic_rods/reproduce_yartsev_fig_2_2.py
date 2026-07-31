"""Reproduce Yartsev Chapter-2 Figure 2.2 for one free-free monoclinic rod.

This is a diagnostic/source-reproduction entry point.  It does not implement
coupled rods, Euler--Bernoulli theory, Saint-Venant torsion, or FEM.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
import time
from dataclasses import replace
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
    BookMaterial,
    Formulation,
    Geometry,
    RootResult,
    assign_modes_by_mac,
    boundary_quality,
    continue_loss_root,
    find_elastic_roots,
    formulation_mapping_residual,
    make_rod_point,
    modal_loss_factors,
    mode_shape,
    rigid_body_nullity,
    solve_complex_root,
)


DEFAULT_RESULT_DIR = ROOT / "results" / "anisotropic_rods" / "yartsev_ch2_free_free"
SMOKE_RESULT_DIR = ROOT / "results" / "_smoke" / "yartsev_ch2_free_free"
CHAPTER_1 = ROOT / "docs" / "literature" / "pdf" / "Глава 1_compressed.pdf"
CHAPTER_2 = ROOT / "docs" / "literature" / "pdf" / "Глава 2_compressed.pdf"

# Manual post-implementation reading of the *calculated solid curves* in
# Figure 2.2 (printed p.56).  Values are deliberately rounded to the graphical
# resolution and are never solver inputs.  Frequencies are in kHz; loss values
# are eta*1e2.  Experimental symbols were not digitized in this gate.
_DIGITIZED_FREQUENCY_KHZ = {
    1: (0.48, 0.42, 0.35, 0.31, 0.30, 0.30, 0.30),
    2: (1.30, 1.15, 0.95, 0.85, 0.81, 0.82, 0.83),
    3: (1.70, 1.85, 1.83, 1.63, 1.57, 1.58, 1.60),
    4: (2.45, 2.20, 2.10, 2.08, 1.90, 1.73, 1.65),
    5: (3.40, 3.50, 3.00, 2.65, 2.55, 2.55, 2.60),
    6: (3.95, 3.80, 4.20, 4.15, 3.80, 3.45, 3.30),
    7: (5.10, 5.05, 4.30, 3.85, 3.70, 3.75, 3.80),
}
_DIGITIZED_ETA_TIMES_100 = {
    1: (0.45, 0.95, 1.55, 1.80, 1.70, 1.50, 1.40),
    2: (0.52, 1.02, 1.58, 1.77, 1.72, 1.53, 1.40),
    3: (2.13, 1.72, 1.64, 1.80, 1.75, 1.55, 1.45),
    4: (0.60, 1.10, 1.40, 1.55, 1.82, 2.10, 2.25),
    5: (2.13, 1.30, 1.68, 1.83, 1.78, 1.60, 1.50),
    6: (0.70, 1.58, 1.40, 1.54, 1.83, 2.12, 2.25),
    7: (2.13, 1.38, 1.74, 1.87, 1.82, 1.67, 1.58),
}
_DIGITIZED_THETA = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--material-mode",
        choices=("elastic", "book-complex", "both"),
        default="both",
    )
    parser.add_argument(
        "--equation-variant",
        choices=("corrected", "printed", "both"),
        default="both",
    )
    parser.add_argument("--theta-step-deg", type=float, default=1.0)
    parser.add_argument("--num-positive-modes", type=int, default=7)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--scan-step-hz", type=float, default=10.0)
    parser.add_argument("--series-rtol", type=float, default=1.0e-12)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args(argv)


def _validate_args(args: argparse.Namespace) -> None:
    if args.theta_step_deg <= 0.0 or args.theta_step_deg > 90.0:
        raise ValueError("--theta-step-deg must lie in (0, 90]")
    if args.num_positive_modes < 1:
        raise ValueError("--num-positive-modes must be positive")
    if args.scan_step_hz <= 0.0:
        raise ValueError("--scan-step-hz must be positive")
    if args.series_rtol <= 0.0:
        raise ValueError("--series-rtol must be positive")


def _theta_grid(args: argparse.Namespace) -> np.ndarray:
    if args.smoke:
        return np.array([0.0, 45.0, 90.0])
    values = list(np.arange(0.0, 90.0, args.theta_step_deg))
    values.append(90.0)
    return np.array(sorted(set(round(float(value), 12) for value in values)))


def _write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    data = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not data:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    for row in data:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _complex_columns(prefix: str, value: complex) -> dict[str, float]:
    return {f"{prefix}_real": float(np.real(value)), f"{prefix}_imag": float(np.imag(value))}


def _material_rows() -> list[dict[str, Any]]:
    material = BookMaterial()
    geometry = Geometry()
    return [
        {"category": "geometry", "parameter": "a", "value_si": geometry.a, "unit": "m", "source": "Chapter 2 p.56"},
        {"category": "geometry", "parameter": "b", "value_si": geometry.b, "unit": "m", "source": "Chapter 2 p.56"},
        {"category": "geometry", "parameter": "L", "value_si": geometry.length, "unit": "m", "source": "Chapter 2 p.56"},
        {"category": "geometry", "parameter": "A", "value_si": geometry.area, "unit": "m^2", "source": "Chapter 2 definition"},
        {"category": "geometry", "parameter": "I_y", "value_si": geometry.I_y, "unit": "m^4", "source": "Chapter 2 p.53"},
        {"category": "geometry", "parameter": "I_p", "value_si": geometry.I_p, "unit": "m^4", "source": "Chapter 2 p.53"},
        {"category": "geometry", "parameter": "k", "value_si": geometry.shear_factor, "unit": "1", "source": "Chapter 2 eq. (2.4)"},
        {"category": "material", "parameter": "Re(E1)", "value_si": material.E1_real, "unit": "Pa", "source": "Chapter 1 Table 1.2 p.46"},
        {"category": "material", "parameter": "Re(E2)", "value_si": material.E2_real, "unit": "Pa", "source": "Chapter 1 Table 1.2 p.46"},
        {"category": "material", "parameter": "Re(G12)", "value_si": material.G12_real, "unit": "Pa", "source": "Chapter 1 Table 1.2 p.46"},
        {"category": "material", "parameter": "Re(G13)", "value_si": material.G13_real, "unit": "Pa", "source": "Chapter 1 Table 1.2 p.46"},
        {"category": "material", "parameter": "Re(G23)", "value_si": material.G23_real, "unit": "Pa", "source": "Chapter 1 Table 1.2 p.46"},
        {"category": "material", "parameter": "nu12", "value_si": material.nu12, "unit": "1", "source": "Chapter 1 Table 1.2 p.46"},
        {"category": "material", "parameter": "rho", "value_si": material.rho, "unit": "kg/m^3", "source": "Chapter 1 Table 1.2 p.46"},
        {"category": "material", "parameter": "eta1", "value_si": material.eta1, "unit": "1", "source": "Chapter 1 Table 1.2 p.46"},
        {"category": "material", "parameter": "eta2", "value_si": material.eta2, "unit": "1", "source": "Chapter 1 Table 1.2 p.46"},
        {"category": "material", "parameter": "eta12", "value_si": material.eta12, "unit": "1", "source": "Chapter 1 Table 1.2 p.46"},
        {"category": "material", "parameter": "eta13", "value_si": material.eta13, "unit": "1", "source": "Chapter 1 Table 1.2 p.46"},
        {"category": "material", "parameter": "eta23", "value_si": material.eta23, "unit": "1", "source": "Chapter 1 Table 1.2 p.46"},
        {"category": "source_only", "parameter": "h_ply", "value_si": material.ply_thickness, "unit": "m", "source": "Chapter 1 Table 1.2 p.46; not used by Chapter-2 1D model"},
    ]


def _write_source_manifest(output_dir: Path) -> None:
    sources = [CHAPTER_1, CHAPTER_2]
    rows = [
        "# Source manifest — Yartsev Chapter-2 single-rod reproduction",
        "",
        "Only the two registered local monograph fragments below were used. No web or replacement rod source was used.",
        "Formula-consistency conclusions are recorded in the tracked source and reproduction notes.",
        "",
        "| relative file | available | bytes | SHA256 | role |",
        "| --- | --- | ---: | --- | --- |",
    ]
    roles = {
        CHAPTER_1: "Chapter 1 equations (1.32), (1.34), (1.41), (1.42), (1.50)–(1.56), Table 1.2 (printed pp. 24–25, 30–31, 40–46)",
        CHAPTER_2: "Chapter 2 equations (2.1)–(2.18), specimen and Figure 2.2 (printed pp. 52–57)",
    }
    for source in sources:
        relative_source = source.relative_to(ROOT).as_posix()
        rows.append(
            f"| `{relative_source}` | {'yes' if source.is_file() else 'no'} | "
            f"{source.stat().st_size if source.is_file() else ''} | "
            f"`{_sha256(source) if source.is_file() else ''}` | {roles[source]} |"
        )
    rows += [
        "",
        "Table 1.2 gives `h_ply = 3.5e-4 m`; it is recorded but is not an input to the one-dimensional Chapter-2 equations.",
        "Rendered inspection pages remained under a Git-ignored local results directory and are not public/tracked artifacts.",
    ]
    (output_dir / "source_manifest.md").write_text("\n".join(rows) + "\n", encoding="utf-8")


def _digitized_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for mode, values in _DIGITIZED_FREQUENCY_KHZ.items():
        for theta, value in zip(_DIGITIZED_THETA, values):
            rows.append(
                {
                    "theta_deg": theta,
                    "book_curve_mode_label": mode,
                    "quantity": "frequency",
                    "digitized_value": value,
                    "unit": "kHz",
                    "estimated_graphical_uncertainty": 0.06,
                    "source": "Figure 2.2a calculated solid line, printed p.56",
                    "method": "manual grid reading after implementation; rounded; not used for fitting",
                }
            )
    for mode, values in _DIGITIZED_ETA_TIMES_100.items():
        for theta, value in zip(_DIGITIZED_THETA, values):
            rows.append(
                {
                    "theta_deg": theta,
                    "book_curve_mode_label": mode,
                    "quantity": "eta_times_100",
                    "digitized_value": value,
                    "unit": "1",
                    "estimated_graphical_uncertainty": 0.06,
                    "source": "Figure 2.2b-d calculated solid line, printed p.56",
                    "method": "manual grid reading after implementation; rounded; not used for fitting",
                }
            )
    return rows


def _assess_digitized_curves(
    theta_grid: np.ndarray,
    complex_roots: dict[tuple[Formulation, float], list[RootResult]],
    *,
    num_modes: int,
    smoke: bool,
    tracked_map: dict[tuple[str, Formulation, float, int], int],
) -> tuple[str, list[dict[str, Any]]]:
    available_theta = {float(value) for value in theta_grid}
    rows: list[dict[str, Any]] = []
    for source_row in _digitized_rows():
        theta = float(source_row["theta_deg"])
        mode = int(source_row["book_curve_mode_label"])
        key = ("state_corrected", theta)
        if theta not in available_theta or key not in complex_roots or mode > num_modes:
            continue
        sorted_mode = next(
            (
                candidate
                for candidate in range(1, num_modes + 1)
                if tracked_map.get(
                    ("book_complex", "state_corrected", theta, candidate)
                )
                == mode
            ),
            mode,
        )
        solved = complex_roots[key][sorted_mode - 1]
        if source_row["quantity"] == "frequency":
            computed = solved.frequency_hz / 1000.0
        else:
            computed = 100.0 * modal_loss_factors(solved.omega)[0]
        difference = computed - float(source_row["digitized_value"])
        uncertainty = float(source_row["estimated_graphical_uncertainty"])
        rows.append(
            {
                **source_row,
                "computed_sorted_mode": sorted_mode,
                "computed_tracked_mode": mode,
                "computed_value": computed,
                "computed_minus_digitized": difference,
                "absolute_difference": abs(difference),
                "difference_over_graphical_uncertainty": abs(difference) / uncertainty,
                "within_estimated_graphical_resolution": abs(difference) <= uncertainty,
            }
        )
    full_count = 2 * 7 * 7
    if smoke or len(rows) < full_count:
        return "PARTIAL_REPRODUCTION", rows
    if rows and all(bool(row["within_estimated_graphical_resolution"]) for row in rows):
        return "PASS_WITHIN_GRAPH_RESOLUTION", rows
    return "NOT_REPRODUCED", rows


def _transformed_rows(theta_grid: np.ndarray, series_rtol: float, include_complex: bool) -> list[dict[str, Any]]:
    modes = ["elastic"] + (["book_complex"] if include_complex else [])
    rows: list[dict[str, Any]] = []
    for material_mode in modes:
        for theta in theta_grid:
            point = make_rod_point(
                float(theta),
                material_mode=material_mode,  # type: ignore[arg-type]
                series_relative_tolerance=series_rtol,
            )
            strict = make_rod_point(
                float(theta),
                material_mode=material_mode,  # type: ignore[arg-type]
                series_relative_tolerance=min(series_rtol * 0.01, 1e-14),
            )
            props = point.properties
            torsion = point.torsion
            row: dict[str, Any] = {
                "theta_deg": float(theta),
                "material_mode": material_mode,
                "loss_scale": point.loss_scale,
                "torsion_terms_used": torsion.terms_used,
                "torsion_estimated_relative_tail": torsion.estimated_relative_tail,
                "strict_torsion_terms_used": strict.torsion.terms_used,
                "Cbar_strict_relative_difference": abs(torsion.Cbar - strict.torsion.Cbar) / abs(strict.torsion.Cbar),
                "C_T_strict_relative_difference": abs(torsion.C_T - strict.torsion.C_T) / abs(strict.torsion.C_T),
            }
            for name in (
                "Sbar11", "Sbar16", "Sbar66", "Sbar55", "Ex", "Gxy", "Gxz",
                "Gxy_bar", "mu_x_xy", "mu_xy_x",
            ):
                row.update(_complex_columns(name, getattr(props, name)))
            row.update(_complex_columns("Cbar", torsion.Cbar))
            row.update(_complex_columns("C_T", torsion.C_T))
            rows.append(row)
    return rows


def _selected_formulations(variant: str) -> tuple[list[Formulation], list[Formulation]]:
    elastic: list[Formulation] = []
    complex_variants: list[Formulation] = []
    if variant in ("corrected", "both"):
        elastic += ["state_corrected", "eliminated_corrected"]
        complex_variants.append("state_corrected")
    if variant in ("printed", "both"):
        elastic.append("eliminated_printed")
        complex_variants.append("eliminated_printed")
    return elastic, complex_variants


def _neighbor_distances(roots: list[RootResult]) -> list[RootResult]:
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


def _solve_complex_sweep(
    theta_grid: np.ndarray,
    elastic: dict[tuple[Formulation, float], list[RootResult]],
    formulations: Sequence[Formulation],
    *,
    num_modes: int,
    series_rtol: float,
) -> dict[tuple[Formulation, float], list[RootResult]]:
    complex_roots: dict[tuple[Formulation, float], list[RootResult]] = {}
    for formulation in formulations:
        previous_elastic: list[RootResult] | None = None
        previous_complex: list[RootResult] | None = None
        for theta_index, theta_value in enumerate(theta_grid):
            theta = float(theta_value)
            elastic_roots = elastic[(formulation, theta)][:num_modes]
            solved_roots: list[RootResult] = []
            for mode_index, elastic_root in enumerate(elastic_roots):
                factory = lambda scale, angle=theta: make_rod_point(
                    angle,
                    material_mode="book_complex",
                    loss_scale=scale,
                    series_relative_tolerance=series_rtol,
                )
                if theta_index == 0 or previous_complex is None or previous_elastic is None:
                    solved = continue_loss_root(
                        factory, formulation, elastic_root.omega.real
                    )
                else:
                    old_elastic = previous_elastic[mode_index].omega.real
                    old_complex = previous_complex[mode_index].omega
                    predictor = complex(
                        elastic_root.omega.real + (old_complex.real - old_elastic),
                        old_complex.imag * elastic_root.omega.real / old_elastic,
                    )
                    solved = solve_complex_root(factory(1.0), formulation, predictor)
                    relative_departure = abs(solved.omega.real - elastic_root.omega.real) / elastic_root.omega.real
                    if solved.status == "rejected_complex_quality" or relative_departure > 0.08:
                        solved = continue_loss_root(
                            factory, formulation, elastic_root.omega.real
                        )
                solved_roots.append(solved)
            solved_roots = _neighbor_distances(solved_roots)
            if any(
                abs(solved_roots[index].omega - solved_roots[index - 1].omega)
                <= 1e-7 * max(abs(solved_roots[index].omega), 1.0)
                for index in range(1, len(solved_roots))
            ):
                # Duplicate continuation targets are re-seeded independently
                # from their current elastic roots.
                solved_roots = [
                    continue_loss_root(
                        lambda scale, angle=theta: make_rod_point(
                            angle,
                            material_mode="book_complex",
                            loss_scale=scale,
                            series_relative_tolerance=series_rtol,
                        ),
                        formulation,
                        item.omega.real,
                    )
                    for item in elastic_roots
                ]
                solved_roots = _neighbor_distances(solved_roots)
            complex_roots[(formulation, theta)] = solved_roots
            previous_elastic = elastic_roots
            previous_complex = solved_roots
            print(
                f"complex {formulation}: theta={theta:g} deg "
                f"({theta_index + 1}/{len(theta_grid)})",
                flush=True,
            )
    return complex_roots


def _track_modes(
    theta_grid: np.ndarray,
    roots: dict[tuple[Formulation, float], list[RootResult]],
    formulations: Sequence[Formulation],
    *,
    material_mode: str,
    num_modes: int,
    series_rtol: float,
) -> tuple[dict[tuple[str, Formulation, float, int], int], list[dict[str, Any]], dict[str, np.ndarray]]:
    tracked_map: dict[tuple[str, Formulation, float, int], int] = {}
    diagnostics: list[dict[str, Any]] = []
    saved_shapes: dict[str, np.ndarray] = {}
    x_grid = np.linspace(0.0, 1.0, 41)
    for formulation in formulations:
        all_shapes = np.zeros(
            (len(theta_grid), num_modes, len(x_grid), 3), dtype=np.complex128
        )
        previous_tracked_shapes: list[np.ndarray] | None = None
        for theta_index, theta_value in enumerate(theta_grid):
            theta = float(theta_value)
            point = make_rod_point(
                theta,
                material_mode=material_mode,  # type: ignore[arg-type]
                series_relative_tolerance=series_rtol,
            )
            current_shapes = [
                mode_shape(item.omega, point, formulation, x_grid)
                for item in roots[(formulation, theta)][:num_modes]
            ]
            all_shapes[theta_index] = np.asarray(current_shapes)
            if previous_tracked_shapes is None:
                assignment = list(range(num_modes))
                mac_matrix = np.eye(num_modes)
            else:
                assignment, mac_matrix = assign_modes_by_mac(
                    previous_tracked_shapes, current_shapes
                )
            next_tracked: list[np.ndarray] = []
            for tracked_index, sorted_index in enumerate(assignment):
                tracked_mode = tracked_index + 1
                sorted_mode = sorted_index + 1
                tracked_map[(material_mode, formulation, theta, sorted_mode)] = tracked_mode
                mac = float(mac_matrix[tracked_index, sorted_index])
                row_values = np.sort(mac_matrix[tracked_index])[::-1]
                margin = float(row_values[0] - row_values[1]) if len(row_values) > 1 else 1.0
                diagnostics.append(
                    {
                        "theta_deg": theta,
                        "material_mode": material_mode,
                        "equation_variant": formulation,
                        "tracked_mode": tracked_mode,
                        "current_sorted_mode": sorted_mode,
                        "mac_from_previous_theta": mac if theta_index else 1.0,
                        "mac_margin": margin if theta_index else 1.0,
                        "sorted_mode_swap_warning": bool(theta_index and sorted_mode != tracked_mode),
                        "low_mac_warning": bool(theta_index and mac < 0.75),
                    }
                )
                next_tracked.append(current_shapes[sorted_index])
            previous_tracked_shapes = next_tracked
        key = f"{material_mode}__{formulation}"
        saved_shapes[f"{key}__x_over_length"] = x_grid
        saved_shapes[f"{key}__theta_deg"] = theta_grid
        saved_shapes[f"{key}__sorted_shapes"] = all_shapes
    return tracked_map, diagnostics, saved_shapes


def _root_rows(
    theta_grid: np.ndarray,
    roots: dict[tuple[Formulation, float], list[RootResult]],
    formulations: Sequence[Formulation],
    *,
    material_mode: str,
    tracked_map: dict[tuple[str, Formulation, float, int], int],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for theta_value in theta_grid:
        theta = float(theta_value)
        for formulation in formulations:
            for mode_index, item in enumerate(roots[(formulation, theta)], start=1):
                eta_exact, eta_approx, eta_difference = modal_loss_factors(item.omega)
                rows.append(
                    {
                        "theta_deg": theta,
                        "material_mode": material_mode,
                        "equation_variant": formulation,
                        "sorted_mode": mode_index,
                        "tracked_mode_if_available": tracked_map.get(
                            (material_mode, formulation, theta, mode_index), ""
                        ),
                        "omega_real_rad_s": item.omega.real,
                        "omega_imag_rad_s": item.omega.imag,
                        "frequency_hz": item.frequency_hz,
                        "eta_exact": eta_exact,
                        "eta_approx": eta_approx,
                        "eta_relative_difference": eta_difference,
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


def _comparison_rows(
    theta_grid: np.ndarray,
    elastic: dict[tuple[Formulation, float], list[RootResult]],
    complex_roots: dict[tuple[Formulation, float], list[RootResult]],
    *,
    num_modes: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    printed_rows: list[dict[str, Any]] = []
    equivalence_rows: list[dict[str, Any]] = []
    for theta_value in theta_grid:
        theta = float(theta_value)
        if ("state_corrected", theta) in elastic and ("eliminated_corrected", theta) in elastic:
            for mode, (state, eliminated) in enumerate(
                zip(elastic[("state_corrected", theta)], elastic[("eliminated_corrected", theta)]), start=1
            ):
                point = make_rod_point(theta)
                equivalence_rows.append(
                    {
                        "theta_deg": theta,
                        "sorted_mode": mode,
                        "state_frequency_hz": state.frequency_hz,
                        "eliminated_corrected_frequency_hz": eliminated.frequency_hz,
                        "absolute_difference_hz": eliminated.frequency_hz - state.frequency_hz,
                        "relative_difference": abs(eliminated.frequency_hz - state.frequency_hz) / state.frequency_hz,
                        "state_relative_singular_residual": state.relative_singular_residual,
                        "eliminated_relative_singular_residual": eliminated.relative_singular_residual,
                        "state_to_eliminated_mapping_residual": formulation_mapping_residual(state.omega, point, samples=7),
                    }
                )
        for material_mode, data in (("elastic", elastic), ("book_complex", complex_roots)):
            if ("state_corrected", theta) not in data or ("eliminated_printed", theta) not in data:
                continue
            for mode, (corrected, printed) in enumerate(
                zip(data[("state_corrected", theta)][:num_modes], data[("eliminated_printed", theta)][:num_modes]), start=1
            ):
                corrected_eta = modal_loss_factors(corrected.omega)[0]
                printed_eta = modal_loss_factors(printed.omega)[0]
                printed_rows.append(
                    {
                        "theta_deg": theta,
                        "material_mode": material_mode,
                        "sorted_mode": mode,
                        "corrected_frequency_hz": corrected.frequency_hz,
                        "printed_frequency_hz": printed.frequency_hz,
                        "printed_minus_corrected_hz": printed.frequency_hz - corrected.frequency_hz,
                        "relative_frequency_difference": abs(printed.frequency_hz - corrected.frequency_hz) / corrected.frequency_hz,
                        "corrected_eta_exact": corrected_eta,
                        "printed_eta_exact": printed_eta,
                        "printed_minus_corrected_eta": printed_eta - corrected_eta,
                    }
                )
    return printed_rows, equivalence_rows


def _plot_reproduction(
    output_dir: Path,
    theta_grid: np.ndarray,
    elastic: dict[tuple[Formulation, float], list[RootResult]],
    complex_roots: dict[tuple[Formulation, float], list[RootResult]],
    *,
    num_modes: int,
    elastic_track_map: dict[tuple[str, Formulation, float, int], int],
    complex_track_map: dict[tuple[str, Formulation, float, int], int],
) -> None:
    if not any(key[0] == "state_corrected" for key in elastic):
        return
    colors = plt.cm.tab10(np.linspace(0.0, 0.9, max(num_modes, 1)))
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2), constrained_layout=True)
    frequency_axis = axes[0, 0]
    corrected_complex = all(("state_corrected", float(theta)) in complex_roots for theta in theta_grid)
    def root_for_tracked_mode(
        data: dict[tuple[Formulation, float], list[RootResult]],
        tracked: dict[tuple[str, Formulation, float, int], int],
        material_mode: str,
        theta: float,
        tracked_mode: int,
    ) -> RootResult:
        sorted_mode = next(
            (
                candidate
                for candidate in range(1, num_modes + 1)
                if tracked.get(
                    (material_mode, "state_corrected", theta, candidate)
                )
                == tracked_mode
            ),
            tracked_mode,
        )
        return data[("state_corrected", theta)][sorted_mode - 1]

    for mode in range(num_modes):
        tracked_mode = mode + 1
        elastic_values = [
            root_for_tracked_mode(
                elastic,
                elastic_track_map,
                "elastic",
                float(theta),
                tracked_mode,
            ).frequency_hz
            / 1000.0
            for theta in theta_grid
        ]
        frequency_axis.plot(theta_grid, elastic_values, linestyle="--", linewidth=1.0, color=colors[mode], alpha=0.65)
        if corrected_complex:
            complex_values = [
                root_for_tracked_mode(
                    complex_roots,
                    complex_track_map,
                    "book_complex",
                    float(theta),
                    tracked_mode,
                ).frequency_hz
                / 1000.0
                for theta in theta_grid
            ]
            frequency_axis.plot(theta_grid, complex_values, linewidth=1.7, color=colors[mode], label=f"mode {mode + 1}")
    frequency_axis.set_title("a) First positive frequencies")
    frequency_axis.set_ylabel("f, kHz")
    frequency_axis.legend(ncol=2, fontsize=8)

    groups = ((1, 3, 5), (2, 4), (6, 7))
    loss_axes = (axes[0, 1], axes[1, 0], axes[1, 1])
    labels = ("b) modes 1, 3, 5", "c) modes 2, 4", "d) modes 6, 7")
    for axis, group, title in zip(loss_axes, groups, labels):
        axis.set_title(title)
        plotted = False
        if corrected_complex:
            for mode_number in group:
                if mode_number > num_modes:
                    continue
                values = [
                    100.0 * modal_loss_factors(
                        root_for_tracked_mode(
                            complex_roots,
                            complex_track_map,
                            "book_complex",
                            float(theta),
                            mode_number,
                        ).omega
                    )[0]
                    for theta in theta_grid
                ]
                axis.plot(theta_grid, values, linewidth=1.7, color=colors[mode_number - 1], label=f"mode {mode_number}")
                plotted = True
            if plotted:
                axis.legend(fontsize=8)
        axis.set_ylabel(r"$\eta_i\times10^2$")
    for axis in axes.ravel():
        axis.set_xlim(0.0, 90.0)
        axis.set_xticks(np.arange(0.0, 91.0, 15.0))
        axis.set_xlabel(r"$\theta$, deg")
        axis.grid(True, alpha=0.25)
    fig.suptitle("Yartsev Chapter 2 — single free-free monoclinic rod\nsource-group continuity labels; solid: book_complex corrected; dashed: elastic corrected")
    fig.savefig(output_dir / "figure_2_2_reproduction.png", dpi=220)
    fig.savefig(output_dir / "figure_2_2_reproduction.pdf")
    plt.close(fig)


def _plot_printed_comparison(
    output_dir: Path,
    theta_grid: np.ndarray,
    elastic: dict[tuple[Formulation, float], list[RootResult]],
    complex_roots: dict[tuple[Formulation, float], list[RootResult]],
    *,
    num_modes: int,
) -> None:
    required = all(
        (variant, float(theta)) in elastic
        for theta in theta_grid
        for variant in ("state_corrected", "eliminated_printed")
    )
    if not required:
        return
    colors = plt.cm.tab10(np.linspace(0.0, 0.9, max(num_modes, 1)))
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2), constrained_layout=True)
    for mode in range(num_modes):
        corrected = np.array([elastic[("state_corrected", float(theta))][mode].frequency_hz / 1000.0 for theta in theta_grid])
        printed = np.array([elastic[("eliminated_printed", float(theta))][mode].frequency_hz / 1000.0 for theta in theta_grid])
        axes[0, 0].plot(theta_grid, corrected, color=colors[mode], linewidth=1.5)
        axes[0, 0].plot(theta_grid, printed, color=colors[mode], linewidth=1.0, linestyle="--")
        axes[0, 1].plot(theta_grid, 1000.0 * (printed - corrected), color=colors[mode], linewidth=1.3, label=f"mode {mode + 1}")
    axes[0, 0].set_title("elastic spectra: corrected solid, printed dashed")
    axes[0, 0].set_ylabel("f, kHz")
    axes[0, 1].set_title("printed − corrected elastic frequency")
    axes[0, 1].set_ylabel("Δf, Hz")
    complex_required = all(
        (variant, float(theta)) in complex_roots
        for theta in theta_grid
        for variant in ("state_corrected", "eliminated_printed")
    )
    if complex_required:
        for mode in range(num_modes):
            corrected_eta = np.array([100.0 * modal_loss_factors(complex_roots[("state_corrected", float(theta))][mode].omega)[0] for theta in theta_grid])
            printed_eta = np.array([100.0 * modal_loss_factors(complex_roots[("eliminated_printed", float(theta))][mode].omega)[0] for theta in theta_grid])
            axes[1, 0].plot(theta_grid, corrected_eta, color=colors[mode], linewidth=1.5)
            axes[1, 0].plot(theta_grid, printed_eta, color=colors[mode], linewidth=1.0, linestyle="--")
            axes[1, 1].plot(theta_grid, printed_eta - corrected_eta, color=colors[mode], linewidth=1.3)
    axes[1, 0].set_title("book_complex losses: corrected solid, printed dashed")
    axes[1, 0].set_ylabel(r"$\eta_i\times10^2$")
    axes[1, 1].set_title("printed − corrected modal loss")
    axes[1, 1].set_ylabel(r"$\Delta\eta_i\times10^2$")
    axes[0, 1].legend(ncol=2, fontsize=8)
    for axis in axes.ravel():
        axis.set_xlim(0.0, 90.0)
        axis.set_xticks(np.arange(0.0, 91.0, 15.0))
        axis.set_xlabel(r"$\theta$, deg")
        axis.grid(True, alpha=0.25)
    fig.suptitle("Diagnostic effect of the signs printed after equation (2.16)")
    fig.savefig(output_dir / "printed_vs_corrected.pdf")
    fig.savefig(output_dir / "printed_vs_corrected.png", dpi=220)
    plt.close(fig)


def _report(
    output_dir: Path,
    args: argparse.Namespace,
    theta_grid: np.ndarray,
    transformed_rows: list[dict[str, Any]],
    elastic: dict[tuple[Formulation, float], list[RootResult]],
    complex_roots: dict[tuple[Formulation, float], list[RootResult]],
    printed_rows: list[dict[str, Any]],
    equivalence_rows: list[dict[str, Any]],
    diagnostics: list[dict[str, Any]],
    reproduction_status: str,
    digitized_comparison: list[dict[str, Any]],
    runtime_seconds: float,
) -> None:
    max_equivalence = max((float(row["relative_difference"]) for row in equivalence_rows), default=math.nan)
    max_mapping = max((float(row["state_to_eliminated_mapping_residual"]) for row in equivalence_rows), default=math.nan)
    max_series = max((float(row["C_T_strict_relative_difference"]) for row in transformed_rows), default=math.nan)
    max_printed = max((float(row["relative_frequency_difference"]) for row in printed_rows if row["material_mode"] == "elastic"), default=math.nan)
    max_exact_approx = max(
        (
            modal_loss_factors(item.omega)[2]
            for roots in complex_roots.values()
            for item in roots
        ),
        default=0.0,
    )
    rigid = rigid_body_nullity(make_rod_point(45.0))
    low_mac = sum(bool(row["low_mac_warning"]) for row in diagnostics)
    swaps = sum(bool(row["sorted_mode_swap_warning"]) for row in diagnostics)
    max_graph_ratio = max(
        (float(row["difference_over_graphical_uncertainty"]) for row in digitized_comparison),
        default=math.nan,
    )
    max_frequency_difference = max(
        (float(row["absolute_difference"]) for row in digitized_comparison if row["quantity"] == "frequency"),
        default=math.nan,
    )
    max_eta_difference = max(
        (float(row["absolute_difference"]) for row in digitized_comparison if row["quantity"] == "eta_times_100"),
        default=math.nan,
    )
    lines = [
        "# Yartsev Chapter-2 single free-free rod reproduction",
        "",
        f"**Reproduction status: `{reproduction_status}`.**",
        "",
        "The internally consistent Chapter-2 model was implemented and compared with a separate, rounded manual reading of the calculated solid curves in Figure 2.2. The book supplies no digital curve table, so this is a graph-resolution assessment rather than an exact percentage comparison. No material parameter was fitted. Root CSV files remain sorted by Re(f) at every theta; the reproduction plot additionally uses MAC-continuity labels to preserve the book's smooth curve grouping through the 6/7 exchanges.",
        "",
        "## Source and printed/corrected distinction",
        "",
        "- Chapter 1 printed pp. 24–25: complex moduli, equations (1.32), (1.34), `M*=Re(M*)[1+i eta]`.",
        "- Chapter 1 printed pp. 30–31: exact and small-loss modal definitions, equations (1.41), (1.42).",
        "- Chapter 1 printed pp. 40–46: compliance data/rotation (1.50)–(1.56), Table 1.2.",
        "- Chapter 2 printed pp. 52–55: equations (2.1)–(2.18); pp. 56–57: specimen and Figure 2.2.",
        "- Literally printed (2.1) omits `I_y` in the rotary-inertia term and is dimensionally inconsistent. It was recorded, not run as a physical variant.",
        "- `state_corrected` restores `rho I_y psi_tt`. `eliminated_corrected` uses positive `d0` and `f0`; `eliminated_printed` preserves their negative printed signs.",
        "",
        "## Numerical gates",
        "",
        f"- theta grid: {len(theta_grid)} points from {theta_grid[0]:g} to {theta_grid[-1]:g} deg.",
        f"- conceptual/numerical rigid-body nullity at omega=0: {rigid}; these roots were excluded from positive tones.",
        f"- maximum state vs eliminated_corrected relative frequency difference: `{max_equivalence:.6e}`.",
        f"- maximum corrected formulation field-mapping residual: `{max_mapping:.6e}`.",
        f"- maximum strict-series change in C_T: `{max_series:.6e}`.",
        f"- maximum printed-vs-corrected elastic sorted-frequency relative difference: `{max_printed:.6e}`.",
        f"- maximum eta_exact-vs-eta_approx relative difference: `{max_exact_approx:.6e}`.",
        f"- MAC warnings below 0.75: {low_mac}; sorted-position swap warnings: {swaps}.",
        f"- maximum frequency discrepancy from rounded Figure-2.2 reading: `{max_frequency_difference:.6g} kHz`.",
        f"- maximum eta*1e2 discrepancy from rounded Figure-2.2 reading: `{max_eta_difference:.6g}`.",
        f"- maximum discrepancy / declared graphical uncertainty: `{max_graph_ratio:.6g}`.",
        "",
        f"Positive elastic root {args.num_positive_modes + 1} is retained as the completeness guard; plots use the requested positive modes only.",
        "",
        "## Outputs and command",
        "",
        f"Run time: `{runtime_seconds:.3f} s` ({'smoke' if args.smoke else 'ordinary'} mode).",
        "",
        "```text",
        "python scripts/analysis/anisotropic_rods/reproduce_yartsev_fig_2_2.py",
        "```",
        "",
        "Generated root CSV files use sorted modes as the primary spectrum and distinguish them from MAC-based tracked-mode diagnostics. `figure_2_2_reproduction.*` uses the continuity labels for the book-style curve grouping, while `mode_shapes.npz` stores normalized `[w/L, psi_b, Phi_t]` or eliminated `[w/L, w', Phi]` fields on a common `x/L` grid.",
        "",
        "## Limits",
        "",
        "- Figure 2.2 values were manually read after implementation, rounded to the graphical resolution, and were not used for fitting.",
        "- The generated digitization CSV covers calculated solid curves only; experimental symbols were not digitized.",
        "- The declared ±0.06 kHz and ±0.06 in eta*1e2 are estimated reading uncertainties, not statistical confidence intervals.",
        "- This is a source-reproduction model, not a promoted general anisotropic API.",
        "- Coupled rods, Euler–Bernoulli, Saint-Venant substitution, and FEM are outside this run and were not implemented/run.",
    ]
    (output_dir / "single_rod_reproduction_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    _validate_args(args)
    required = [CHAPTER_1, CHAPTER_2]
    missing = [path for path in required if not path.is_file()]
    if missing:
        print("Missing required source files:", file=sys.stderr)
        for path in missing:
            print(f"- {path}", file=sys.stderr)
        return 2

    if args.smoke:
        args.num_positive_modes = min(args.num_positive_modes, 3)
    output_dir = (args.output_dir or (SMOKE_RESULT_DIR if args.smoke else DEFAULT_RESULT_DIR)).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    theta_grid = _theta_grid(args)
    start = time.perf_counter()

    _write_source_manifest(output_dir)
    _write_csv(output_dir / "material_parameters.csv", _material_rows())
    _write_csv(
        output_dir / "figure_2_2_digitized_calculated_curves.csv",
        _digitized_rows(),
    )
    include_complex = args.material_mode in ("book-complex", "both")
    transformed_rows = _transformed_rows(theta_grid, args.series_rtol, include_complex)
    _write_csv(output_dir / "transformed_properties_vs_theta.csv", transformed_rows)

    elastic_formulations, complex_formulations = _selected_formulations(args.equation_variant)
    if args.material_mode == "elastic":
        complex_formulations = []
    guard_roots = args.num_positive_modes + 1
    elastic: dict[tuple[Formulation, float], list[RootResult]] = {}
    for formulation in elastic_formulations:
        for theta_index, theta_value in enumerate(theta_grid):
            theta = float(theta_value)
            point = make_rod_point(
                theta,
                material_mode="elastic",
                series_relative_tolerance=args.series_rtol,
            )
            elastic[(formulation, theta)] = find_elastic_roots(
                point,
                formulation,
                num_roots=guard_roots,
                scan_step_hz=args.scan_step_hz,
            )
            print(
                f"elastic {formulation}: theta={theta:g} deg "
                f"({theta_index + 1}/{len(theta_grid)})",
                flush=True,
            )

    complex_roots: dict[tuple[Formulation, float], list[RootResult]] = {}
    if include_complex:
        complex_roots = _solve_complex_sweep(
            theta_grid,
            elastic,
            complex_formulations,
            num_modes=args.num_positive_modes,
            series_rtol=args.series_rtol,
        )

    elastic_track_map, elastic_diagnostics, elastic_shapes = _track_modes(
        theta_grid,
        elastic,
        elastic_formulations,
        material_mode="elastic",
        num_modes=args.num_positive_modes,
        series_rtol=args.series_rtol,
    )
    complex_track_map: dict[tuple[str, Formulation, float, int], int] = {}
    complex_diagnostics: list[dict[str, Any]] = []
    complex_shapes: dict[str, np.ndarray] = {}
    if complex_roots:
        complex_track_map, complex_diagnostics, complex_shapes = _track_modes(
            theta_grid,
            complex_roots,
            complex_formulations,
            material_mode="book_complex",
            num_modes=args.num_positive_modes,
            series_rtol=args.series_rtol,
        )

    elastic_rows = _root_rows(
        theta_grid,
        elastic,
        elastic_formulations,
        material_mode="elastic",
        tracked_map=elastic_track_map,
    )
    complex_rows = _root_rows(
        theta_grid,
        complex_roots,
        complex_formulations,
        material_mode="book_complex",
        tracked_map=complex_track_map,
    ) if complex_roots else []
    _write_csv(output_dir / "elastic_roots.csv", elastic_rows)
    _write_csv(output_dir / "complex_roots.csv", complex_rows)
    _write_csv(output_dir / "root_quality.csv", [*elastic_rows, *complex_rows])
    _write_csv(output_dir / "mode_continuity.csv", [*elastic_diagnostics, *complex_diagnostics])
    np.savez_compressed(output_dir / "mode_shapes.npz", **elastic_shapes, **complex_shapes)

    printed_rows, equivalence_rows = _comparison_rows(
        theta_grid,
        elastic,
        complex_roots,
        num_modes=args.num_positive_modes,
    )
    _write_csv(output_dir / "printed_vs_corrected.csv", printed_rows)
    _write_csv(output_dir / "formulation_equivalence.csv", equivalence_rows)
    reproduction_status, digitized_comparison = _assess_digitized_curves(
        theta_grid,
        complex_roots,
        num_modes=args.num_positive_modes,
        smoke=args.smoke,
        tracked_map=complex_track_map,
    )
    _write_csv(
        output_dir / "figure_2_2_digitized_comparison.csv",
        digitized_comparison,
    )
    _plot_reproduction(
        output_dir,
        theta_grid,
        elastic,
        complex_roots,
        num_modes=args.num_positive_modes,
        elastic_track_map=elastic_track_map,
        complex_track_map=complex_track_map,
    )
    _plot_printed_comparison(
        output_dir,
        theta_grid,
        elastic,
        complex_roots,
        num_modes=args.num_positive_modes,
    )
    runtime_seconds = time.perf_counter() - start
    _report(
        output_dir,
        args,
        theta_grid,
        transformed_rows,
        elastic,
        complex_roots,
        printed_rows,
        equivalence_rows,
        [*elastic_diagnostics, *complex_diagnostics],
        reproduction_status,
        digitized_comparison,
        runtime_seconds,
    )
    print(f"output_dir={output_dir}")
    print(f"runtime_seconds={runtime_seconds:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
