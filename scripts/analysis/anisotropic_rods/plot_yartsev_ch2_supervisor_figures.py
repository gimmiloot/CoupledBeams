"""Build supervisor figures for the Yartsev Chapter-2 rod workflow.

This presentation-only entry point is deliberately isolated from the parallel
isotropic circular-rod article.  Scientific matrices come only from the
Chapter-2 monoclinic-rod, coupled-rod, and rectangular-EB helpers.  Figure 1
reuses the canonical verified Figure-2.2 CSV evidence; Figures 2--4 contain
independently sorted positive spectra at every beta and do no mode tracking.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Callable, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib.yartsev_ch2_coupled_rods import (  # noqa: E402
    coupled_boundary_matrix,
    coupled_boundary_matrix_raw,
    equilibrate_matrix,
    joint_matrix_book,
)
from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
    BookMaterial,
    Geometry,
    RodPoint,
    RootResult,
    cantilever_clamp_matrix,
    find_elastic_roots,
    hms_dx_209_material,
    make_rod_point,
    physical_state_transfer_matrix,
)
from scripts.lib.yartsev_ch2_rectangular_eb import (  # noqa: E402
    eb_clamp_matrix,
    eb_coupled_boundary_matrix,
    eb_coupled_boundary_matrix_raw,
    eb_state_transfer_matrix,
)
from scripts.lib.yartsev_ch2_fast_beta_sweep import (  # noqa: E402
    FAST_SOLVER_VERSION,
    ExactTransferLRU,
    FastSolverValidationError,
    FastSweepResult,
    FastSweepSettings,
    PerformanceCounters,
    SearchInterval,
    load_family_checkpoint,
    maximum_relative_frequency_error,
    predict_sorted_frequencies,
    run_fast_beta_sweep,
    stable_fingerprint,
    write_family_checkpoint,
)


DEFAULT_OUTPUT_DIR = (
    ROOT / "results" / "anisotropic_rods" / "yartsev_ch2_supervisor_figures"
)
CANONICAL_FIGURE_1_DIR = (
    ROOT / "results" / "anisotropic_rods" / "yartsev_ch2_free_free"
)
CANONICAL_FIGURE_1_ROOTS = CANONICAL_FIGURE_1_DIR / "complex_roots.csv"
CANONICAL_FIGURE_1_DIGITIZED = (
    CANONICAL_FIGURE_1_DIR / "figure_2_2_digitized_calculated_curves.csv"
)
CANONICAL_FIGURE_1_COMPARISON = (
    CANONICAL_FIGURE_1_DIR / "figure_2_2_digitized_comparison.csv"
)
CANONICAL_FIGURE_1_REPORT = (
    CANONICAL_FIGURE_1_DIR / "single_rod_reproduction_report.md"
)

FIGURE_BASENAMES = {
    1: "figure_01_yartsev_fig_2_2_reproduction",
    2: "figure_02_clamp_comparison_lambda_vs_beta_book_material",
    3: "figure_03_timoshenko_vs_eb_book_slope_lambda_vs_beta",
    4: "figure_04_timoshenko_vs_eb_section_clamp_lambda_vs_beta",
    5: "figure_05_timoshenko_vs_eb_unequal_lengths_book_slope",
    6: "figure_06_timoshenko_vs_eb_unequal_lengths_and_thickness_book_slope",
    7: "figure_07_monoclinic_theta5_vs_orthotropic_eb_approximation",
    8: "figure_08_chapter2_theta15_vs_theta0",
}
DATA_FILENAMES = {number: f"figure_{number:02d}_data.csv" for number in range(1, 9)}

DEFAULT_BETA_STEP_DEG = 0.5
PLOTTED_ROOT_COUNT = 6
GUARD_ROOT_COUNT = 7
ROOT_DETERMINANT_TOLERANCE = 1.0e-8
ROOT_SINGULAR_TOLERANCE = 1.0e-8
SCAN_STEP_HZ = 10.0
BASE_POINT_SCAN_STEP_HZ = 0.5
COMPLETENESS_REFINEMENT_STEPS_HZ = (2.0, 0.5)
LOCAL_HINT_MAX_SCAN_STEP_HZ = 0.5
LOCAL_HINT_FALLBACK_SCAN_STEP_HZ = 5.0e-3
CLOSE_ROOT_HINT_GAP_HZ = 4.0 * SCAN_STEP_HZ
SPECTRAL_STEP_JUMP_TRIGGER = 3.0e-2
INITIAL_MAX_HZ = 5000.0
MAX_HZ = 100_000.0
FAST_GLOBAL_INITIAL_MAX_HZ = 2000.0
FORMULATION = "state_corrected"
MODE_COLORS = (
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#D55E00",
    "#CC79A7",
    "#56B4E9",
)
FIGURE_1_COLORS = (*MODE_COLORS, "#7A7A7A")
SOLID_LINEWIDTH = 1.6
DASHED_LINEWIDTH = 1.3
DASH_PATTERN = (5.0, 2.4)
COMPARISON_FIGURE_SIZE_IN = (7.2, 4.8)
FIGURE_1_SIZE_IN = (11.0, 8.2)
PNG_DPI = 300
FONT_FAMILY = "DejaVu Sans"
REFERENCE_A_M = 0.005
REFERENCE_B_M = 0.020
FAST_FAMILY_DIRNAME = "fast_family_checkpoints"
FAST_VALIDATION_FILENAME = "fast_solver_validation.csv"
FAST_BENCHMARK_FILENAME = "fast_solver_benchmark.json"
LEGACY_RECORDED_RUNTIME_S = 860.7748033000062
ARTICLE_WORKSPACES = (
    ROOT / "paper_dorofeev_style",
    ROOT / "paper_thickness_mismatch_timoshenko",
)


@dataclass(frozen=True)
class FigurePreset:
    """Exact rectangular Chapter-2 material and two-arm geometry."""

    figure_numbers: tuple[int, ...]
    material_name: str
    material_factory: Callable[[], BookMaterial]
    a_m: float
    b_m: float
    length_1_m: float
    length_2_m: float
    a_2_m: float | None = None
    b_2_m: float | None = None
    theta_1_deg: float = 0.0
    theta_2_deg: float = 0.0
    material_mode: str = "elastic"
    mu: float = 0.0
    shear_factor: float = 5.0 / 6.0


FIGURE_2_PRESET = FigurePreset(
    figure_numbers=(2,),
    material_name="T-53(VM)-78/PN-609-21M (BookMaterial)",
    material_factory=BookMaterial,
    a_m=9.76e-3,
    b_m=2.524e-2,
    length_1_m=0.295,
    length_2_m=0.295,
)
FIGURES_3_4_PRESET = FigurePreset(
    figure_numbers=(3, 4),
    material_name="HMS/DX-209",
    material_factory=hms_dx_209_material,
    a_m=0.005,
    b_m=0.020,
    length_1_m=0.400,
    length_2_m=0.400,
)
FIGURE_5_PRESET = FigurePreset(
    figure_numbers=(5,),
    material_name="HMS/DX-209",
    material_factory=hms_dx_209_material,
    a_m=0.005,
    b_m=0.020,
    length_1_m=0.300,
    length_2_m=0.500,
    mu=0.25,
)
FIGURE_6_PRESET = FigurePreset(
    figure_numbers=(6,),
    material_name="HMS/DX-209",
    material_factory=hms_dx_209_material,
    a_m=0.004,
    b_m=0.020,
    length_1_m=0.300,
    length_2_m=0.500,
    a_2_m=0.006,
    b_2_m=0.020,
    mu=0.25,
)
FIGURE_7_PRESET = FigurePreset(
    figure_numbers=(7,),
    material_name="HMS/DX-209",
    material_factory=hms_dx_209_material,
    a_m=0.005,
    b_m=0.020,
    length_1_m=0.400,
    length_2_m=0.400,
    theta_1_deg=5.0,
    theta_2_deg=5.0,
)
FIGURE_8_PRESET = FigurePreset(
    figure_numbers=(8,),
    material_name="HMS/DX-209",
    material_factory=hms_dx_209_material,
    a_m=0.005,
    b_m=0.020,
    length_1_m=0.400,
    length_2_m=0.400,
    theta_1_deg=15.0,
    theta_2_deg=15.0,
)
FIGURE_PRESETS = (
    FIGURE_2_PRESET,
    FIGURES_3_4_PRESET,
    FIGURE_5_PRESET,
    FIGURE_6_PRESET,
    FIGURE_7_PRESET,
    FIGURE_8_PRESET,
)


@dataclass
class SpectrumResult:
    roots: list[RootResult]
    quality: list[dict[str, Any]]
    scaled_matrix_evaluations: int
    raw_quality_matrix_evaluations: int
    runtime_seconds: float
    scan_step_hz: float = SCAN_STEP_HZ
    completeness_refinement_attempts: int = 0


@dataclass
class SweepResult:
    spectra: dict[float, SpectrumResult]
    runtime_seconds: float
    scaled_matrix_evaluations: int
    raw_quality_matrix_evaluations: int
    completeness_refined_beta_count: int
    completeness_refinement_attempts: int


@dataclass
class SectionReferenceResult:
    preset_label: str
    coupled: SpectrumResult
    straight: SpectrumResult
    maximum_relative_frequency_difference: float
    status: str


class _CountingBoundaryBuilder:
    """Cache and count unique equilibrated matrices used by the root finder."""

    def __init__(self, factory: Callable[[complex], NDArray[np.complex128]]) -> None:
        self.factory = factory
        self.cache: dict[complex, NDArray[np.complex128]] = {}
        self.evaluations = 0

    def __call__(
        self, omega: complex, _point: RodPoint, _formulation: str
    ) -> NDArray[np.complex128]:
        key = complex(omega)
        if key not in self.cache:
            self.cache[key] = np.asarray(self.factory(key), dtype=np.complex128)
            self.evaluations += 1
        return self.cache[key]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--figure",
        choices=("1", "2", "3", "4", "5", "6", "7", "8", "new", "all"),
        default="all",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    data_group = parser.add_mutually_exclusive_group()
    data_group.add_argument("--reuse-data", action="store_true")
    data_group.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--beta-step-deg", type=float, default=DEFAULT_BETA_STEP_DEG)
    parser.add_argument("--solver-mode", choices=("legacy", "fast"), default="fast")
    action_group = parser.add_mutually_exclusive_group()
    action_group.add_argument("--validate-fast-solver", action="store_true")
    action_group.add_argument("--benchmark-fast-solver", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--jobs", type=int, default=1)
    return parser.parse_args(argv)


def beta_grid(beta_step_deg: float = DEFAULT_BETA_STEP_DEG) -> NDArray[np.float64]:
    """Return an inclusive 0--90 degree grid without assuming step divisibility."""

    if not math.isfinite(beta_step_deg) or beta_step_deg <= 0.0 or beta_step_deg > 90.0:
        raise ValueError("--beta-step-deg must lie in (0, 90]")
    values = list(np.arange(0.0, 90.0, float(beta_step_deg)))
    values.append(90.0)
    return np.asarray(sorted({round(float(value), 12) for value in values}), dtype=float)


def rectangular_reference_section(a_m: float, b_m: float) -> tuple[float, float]:
    """Return the required rectangular ``A=ab`` and ``I_y=a^3 b/12``."""

    if a_m <= 0.0 or b_m <= 0.0:
        raise ValueError("rectangular dimensions must be positive")
    return float(a_m * b_m), float(a_m**3 * b_m / 12.0)


def lambda_from_omega(
    omega_rad_s: float | NDArray[np.float64],
    *,
    rho_kg_m3: float,
    a_m: float,
    b_m: float,
    length_1_m: float,
    length_2_m: float,
    elastic_ex_pa: float,
) -> float | NDArray[np.float64]:
    """Evaluate ``(rho A omega^2 l^4 / (E_x I_y))^(1/4)``."""

    omega = np.asarray(omega_rad_s, dtype=float)
    if np.any(omega < 0.0):
        raise ValueError("omega must be non-negative")
    area, inertia_y = rectangular_reference_section(a_m, b_m)
    reference_length = 0.5 * (length_1_m + length_2_m)
    value = (
        rho_kg_m3
        * area
        * omega**2
        * reference_length**4
        / (elastic_ex_pa * inertia_y)
    ) ** 0.25
    return float(value) if value.ndim == 0 else value


def lambda_from_omega_equivalent(
    omega_rad_s: float | NDArray[np.float64],
    *,
    rho_kg_m3: float,
    a_m: float,
    b_m: float,
    length_1_m: float,
    length_2_m: float,
    elastic_ex_pa: float,
) -> float | NDArray[np.float64]:
    """Evaluate the algebraically equivalent ``l (rho A/E_x I_y)^1/4 sqrt(omega)``."""

    omega = np.asarray(omega_rad_s, dtype=float)
    if np.any(omega < 0.0):
        raise ValueError("omega must be non-negative")
    area, inertia_y = rectangular_reference_section(a_m, b_m)
    reference_length = 0.5 * (length_1_m + length_2_m)
    value = (
        reference_length
        * (rho_kg_m3 * area / (elastic_ex_pa * inertia_y)) ** 0.25
        * np.sqrt(omega)
    )
    return float(value) if value.ndim == 0 else value


def _point(preset: FigurePreset, arm: int, *, total_length: bool = False) -> RodPoint:
    if arm not in (1, 2):
        raise ValueError("arm must be 1 or 2")
    length = (
        preset.length_1_m + preset.length_2_m
        if total_length
        else (preset.length_1_m if arm == 1 else preset.length_2_m)
    )
    theta = preset.theta_1_deg if arm == 1 else preset.theta_2_deg
    a_m = preset.a_m if arm == 1 or preset.a_2_m is None else preset.a_2_m
    b_m = preset.b_m if arm == 1 or preset.b_2_m is None else preset.b_2_m
    geometry = Geometry(
        a=a_m,
        b=b_m,
        length=length,
        shear_factor=preset.shear_factor,
    )
    return make_rod_point(
        theta,
        material_mode="elastic",
        geometry=geometry,
        material=preset.material_factory(),
    )


def _section_clamp_coupled_boundary_matrix_raw(
    omega: complex, beta_rad: float, point_1: RodPoint, point_2: RodPoint
) -> NDArray[np.complex128]:
    """Script-local ``J_book blockdiag(T_i B_i_section)`` builder."""

    block = np.zeros((12, 6), dtype=np.complex128)
    for offset_row, offset_column, point in ((0, 0, point_1), (6, 3, point_2)):
        section_clamp = cantilever_clamp_matrix(
            point, "timoshenko_section_clamp", scaled=False
        )
        block[offset_row : offset_row + 6, offset_column : offset_column + 3] = (
            physical_state_transfer_matrix(omega, point) @ section_clamp
        )
    return joint_matrix_book(beta_rad) @ block


def _section_clamp_coupled_boundary_matrix(
    omega: complex,
    beta_rad: float,
    point_1: RodPoint,
    point_2: RodPoint,
    *,
    scaled: bool = True,
) -> NDArray[np.complex128]:
    raw = _section_clamp_coupled_boundary_matrix_raw(
        omega, beta_rad, point_1, point_2
    )
    return equilibrate_matrix(raw)[0] if scaled else raw


def _straight_section_clamp_boundary_matrix_raw(
    omega: complex, point_total: RodPoint
) -> NDArray[np.complex128]:
    """Independent homogeneous straight reference without ``J_book``."""

    selector = np.zeros((3, 6), dtype=np.complex128)
    selector[:, :3] = np.eye(3, dtype=np.complex128)
    return (
        selector
        @ physical_state_transfer_matrix(omega, point_total)
        @ cantilever_clamp_matrix(
            point_total, "timoshenko_section_clamp", scaled=False
        )
    )


def _straight_section_clamp_boundary_matrix(
    omega: complex, point_total: RodPoint, *, scaled: bool = True
) -> NDArray[np.complex128]:
    raw = _straight_section_clamp_boundary_matrix_raw(omega, point_total)
    return equilibrate_matrix(raw)[0] if scaled else raw


def _matrix_quality(matrix: NDArray[np.complex128]) -> dict[str, float]:
    value = np.asarray(matrix, dtype=np.complex128)
    singular = np.linalg.svd(value, compute_uv=False)
    sigma_max = float(singular[0])
    sigma_min = float(singular[-1])
    determinant = complex(np.linalg.det(value))
    determinant_scale = max(sigma_max ** value.shape[0], np.finfo(float).tiny)
    return {
        "physical_raw_determinant_abs": float(abs(determinant)),
        "physical_raw_normalized_determinant_residual": float(
            abs(determinant) / determinant_scale
        ),
        "physical_raw_sigma_min": sigma_min,
        "physical_raw_sigma_max": sigma_max,
        "physical_raw_relative_singular_residual": (
            sigma_min / sigma_max if sigma_max else 0.0
        ),
    }


def _root_quality_fields(
    root: RootResult, raw_factory: Callable[[complex], NDArray[np.complex128]]
) -> dict[str, Any]:
    """Apply the unchanged scaled-or-normalized-physical-raw acceptance rule."""

    raw = _matrix_quality(raw_factory(root.omega))
    status_ok = not root.status.startswith("rejected")
    scaled_ok = (
        root.determinant_residual <= ROOT_DETERMINANT_TOLERANCE
        and root.relative_singular_residual <= ROOT_SINGULAR_TOLERANCE
        and status_ok
    )
    raw_ok = (
        raw["physical_raw_normalized_determinant_residual"]
        <= ROOT_DETERMINANT_TOLERANCE
        and raw["physical_raw_relative_singular_residual"]
        <= ROOT_SINGULAR_TOLERANCE
        and status_ok
    )
    basis = "scaled" if scaled_ok else "physical_raw" if raw_ok else "none"
    return {
        "scaled_determinant_residual": root.determinant_residual,
        "scaled_relative_singular_residual": root.relative_singular_residual,
        **raw,
        "root_status": root.status,
        "quality_status": "PASS" if scaled_ok or raw_ok else "FAIL",
        "quality_basis": basis,
        "accepted_determinant_residual": (
            root.determinant_residual
            if scaled_ok
            else raw["physical_raw_normalized_determinant_residual"]
        ),
        "accepted_relative_singular_residual": (
            root.relative_singular_residual
            if scaled_ok
            else raw["physical_raw_relative_singular_residual"]
        ),
    }


def _validate_spectrum(roots: Sequence[RootResult], quality: Sequence[Mapping[str, Any]]) -> None:
    if len(roots) != GUARD_ROOT_COUNT:
        raise RuntimeError(
            f"expected {GUARD_ROOT_COUNT} positive roots, found {len(roots)}"
        )
    frequencies = np.asarray([root.frequency_hz for root in roots], dtype=float)
    if not np.all(np.isfinite(frequencies)) or np.any(frequencies <= 0.0):
        raise RuntimeError("spectrum contains a non-finite or non-positive root")
    if np.any(np.diff(frequencies) <= 0.0):
        raise RuntimeError("spectrum is not strictly sorted or contains duplicates")
    relative_gaps = np.diff(frequencies) / np.maximum(
        frequencies[:-1], frequencies[1:]
    )
    if np.any(relative_gaps <= 2.0e-8):
        raise RuntimeError("duplicate/indistinguishable roots detected")
    failures = [index + 1 for index, item in enumerate(quality) if item["quality_status"] != "PASS"]
    if failures:
        raise RuntimeError(f"root-quality gate failed for modes {failures}")


def _solve_spectrum(
    point_for_contract: RodPoint,
    scaled_factory: Callable[[complex], NDArray[np.complex128]],
    raw_factory: Callable[[complex], NDArray[np.complex128]],
    *,
    scan_step_hz: float = SCAN_STEP_HZ,
) -> SpectrumResult:
    builder = _CountingBoundaryBuilder(scaled_factory)
    started = time.perf_counter()
    roots = find_elastic_roots(
        point_for_contract,
        FORMULATION,  # type: ignore[arg-type]
        num_roots=GUARD_ROOT_COUNT,
        scan_step_hz=scan_step_hz,
        initial_max_hz=INITIAL_MAX_HZ,
        max_hz=MAX_HZ,
        boundary_matrix_builder=builder,
    )
    quality = [_root_quality_fields(root, raw_factory) for root in roots]
    _validate_spectrum(roots, quality)
    return SpectrumResult(
        roots=list(roots),
        quality=quality,
        scaled_matrix_evaluations=builder.evaluations,
        raw_quality_matrix_evaluations=len(roots),
        runtime_seconds=time.perf_counter() - started,
        scan_step_hz=float(scan_step_hz),
    )


def _spectral_step_difference(
    previous: SpectrumResult, current: SpectrumResult
) -> tuple[float, int]:
    previous_frequency = np.asarray(
        [root.frequency_hz for root in previous.roots], dtype=float
    )
    current_frequency = np.asarray(
        [root.frequency_hz for root in current.roots], dtype=float
    )
    relative = np.abs(current_frequency - previous_frequency) / np.maximum(
        previous_frequency, 1.0
    )
    index = int(np.argmax(relative))
    return float(relative[index]), index + 1


def _retain_final_attempt_with_accounting(
    attempts: Sequence[SpectrumResult],
) -> SpectrumResult:
    final = attempts[-1]
    final.scaled_matrix_evaluations = sum(
        item.scaled_matrix_evaluations for item in attempts
    )
    final.raw_quality_matrix_evaluations = sum(
        item.raw_quality_matrix_evaluations for item in attempts
    )
    final.runtime_seconds = sum(item.runtime_seconds for item in attempts)
    final.completeness_refinement_attempts = sum(
        item.completeness_refinement_attempts for item in attempts
    ) + len(attempts) - 1
    return final


def _recover_spectrum_from_sorted_hints(
    point_for_contract: RodPoint,
    scaled_factory: Callable[[complex], NDArray[np.complex128]],
    raw_factory: Callable[[complex], NDArray[np.complex128]],
    previous_previous: SpectrumResult,
    previous: SpectrumResult,
    incomplete: SpectrumResult,
    *,
    hint_indices: Sequence[int] | None = None,
) -> SpectrumResult:
    """Use the existing root finder in narrow offset windows around prior roots.

    This is only a completeness repair after the ordinary full scan and two
    finer full scans disagree sharply with the preceding sorted spectrum.  The
    prior sorted values define numerical windows, not modal identities.  Each
    local solve still calls the canonical ``find_elastic_roots`` and the final
    current-beta array is independently sorted and quality-gated.
    """

    older = np.asarray(
        [root.frequency_hz for root in previous_previous.roots], dtype=float
    )
    prior = np.asarray([root.frequency_hz for root in previous.roots], dtype=float)
    current = np.asarray([root.frequency_hz for root in incomplete.roots], dtype=float)
    relative = np.abs(current - prior) / np.maximum(prior, 1.0)
    suspicious = (
        np.asarray(sorted(set(int(value) for value in hint_indices)), dtype=int)
        if hint_indices is not None
        else np.flatnonzero(relative > SPECTRAL_STEP_JUMP_TRIGGER)
    )
    if suspicious.size == 0:
        return incomplete

    predicted = prior + (prior - older)
    windows: list[tuple[float, float, float]] = []
    for index in suspicious:
        half_width = max(
            5.0,
            4.0 * abs(prior[index] - older[index]),
            2.0e-3 * abs(predicted[index]),
        )
        windows.append(
            (
                max(1.0e-3, float(predicted[index] - half_width)),
                float(predicted[index] + half_width),
                float(predicted[index]),
            )
        )
    windows.sort()
    merged: list[dict[str, Any]] = []
    for lower, upper, center in windows:
        if not merged or lower > float(merged[-1]["upper"]):
            merged.append({"lower": lower, "upper": upper, "centers": [center]})
        else:
            merged[-1]["upper"] = max(float(merged[-1]["upper"]), upper)
            merged[-1]["centers"].append(center)

    candidates: list[RootResult] = list(incomplete.roots)
    evaluations = 0
    minimum_local_step = LOCAL_HINT_MAX_SCAN_STEP_HZ
    started = time.perf_counter()
    for interval in merged:
        lower = float(interval["lower"])
        upper = float(interval["upper"])
        centers = sorted(float(value) for value in interval["centers"])
        expected = len(centers)
        width = upper - lower
        predicted_gaps = [
            right - left
            for left, right in zip(centers, centers[1:])
            if right - left > 1.0e-8 * max(abs(left), abs(right), 1.0)
        ]
        gap_step = (
            min(predicted_gaps) / 8.0
            if predicted_gaps
            else LOCAL_HINT_FALLBACK_SCAN_STEP_HZ
        )
        local_scan_step = max(
            1.0e-4,
            min(
                LOCAL_HINT_MAX_SCAN_STEP_HZ,
                width / max(20.0 * expected, 20.0),
                gap_step,
            ),
        )
        minimum_local_step = min(minimum_local_step, local_scan_step)
        offset = 2.0 * math.pi * lower
        local_builder = _CountingBoundaryBuilder(
            lambda local_omega, shift=offset: scaled_factory(local_omega + shift)
        )
        local_roots = find_elastic_roots(
            point_for_contract,
            FORMULATION,  # type: ignore[arg-type]
            num_roots=expected,
            scan_step_hz=local_scan_step,
            initial_max_hz=width,
            max_hz=width,
            boundary_matrix_builder=local_builder,
        )
        evaluations += local_builder.evaluations
        for local_root in local_roots:
            candidates.append(
                replace(
                    local_root,
                    omega=local_root.omega + offset,
                    frequency_hz=local_root.frequency_hz + lower,
                    status=f"{local_root.status}_sorted_hint_window",
                )
            )

    candidates.sort(key=lambda root: root.frequency_hz)
    unique: list[RootResult] = []
    for candidate in candidates:
        if any(
            abs(candidate.frequency_hz - existing.frequency_hz)
            <= 2.0e-8
            * max(abs(candidate.frequency_hz), abs(existing.frequency_hz), 1.0)
            for existing in unique
        ):
            continue
        unique.append(candidate)
    final_roots = unique[:GUARD_ROOT_COUNT]
    for index, root in enumerate(final_roots):
        neighbors: list[float] = []
        if index:
            neighbors.append(root.frequency_hz - final_roots[index - 1].frequency_hz)
        if index + 1 < len(final_roots):
            neighbors.append(final_roots[index + 1].frequency_hz - root.frequency_hz)
        final_roots[index] = replace(
            root,
            min_neighbor_distance_hz=min(neighbors) if neighbors else math.inf,
        )
    quality = [_root_quality_fields(root, raw_factory) for root in final_roots]
    _validate_spectrum(final_roots, quality)
    return SpectrumResult(
        roots=final_roots,
        quality=quality,
        scaled_matrix_evaluations=evaluations,
        raw_quality_matrix_evaluations=len(final_roots),
        runtime_seconds=time.perf_counter() - started,
        scan_step_hz=minimum_local_step,
        completeness_refinement_attempts=1,
    )


def _factories(
    model: str,
    beta_deg: float,
    point_1: RodPoint,
    point_2: RodPoint,
) -> tuple[
    Callable[[complex], NDArray[np.complex128]],
    Callable[[complex], NDArray[np.complex128]],
]:
    beta_rad = math.radians(beta_deg)
    if model == "book_slope_clamp":
        return (
            lambda omega: coupled_boundary_matrix(
                omega, beta_rad, point_1, point_2, scaled=True
            ),
            lambda omega: coupled_boundary_matrix_raw(
                omega, beta_rad, point_1, point_2
            ),
        )
    if model == "timoshenko_section_clamp":
        return (
            lambda omega: _section_clamp_coupled_boundary_matrix(
                omega, beta_rad, point_1, point_2, scaled=True
            ),
            lambda omega: _section_clamp_coupled_boundary_matrix_raw(
                omega, beta_rad, point_1, point_2
            ),
        )
    if model == "rectangular_eb_saint_venant":
        return (
            lambda omega: eb_coupled_boundary_matrix(
                omega, beta_rad, point_1, point_2, scaled=True
            ),
            lambda omega: eb_coupled_boundary_matrix_raw(
                omega, beta_rad, point_1, point_2
            ),
        )
    raise ValueError(f"unsupported Chapter-2 spectrum model: {model}")


def _cached_factories(
    model: str,
    beta_deg: float,
    point_1: RodPoint,
    point_2: RodPoint,
    transfer_cache: ExactTransferLRU,
    counters: PerformanceCounters,
) -> tuple[
    Callable[[complex], NDArray[np.complex128]],
    Callable[[complex], NDArray[np.complex128]],
]:
    """Build the unchanged physical matrix from cached exact arm transfers."""

    beta_rad = math.radians(beta_deg)
    if model in ("book_slope_clamp", "timoshenko_section_clamp"):
        clamp_variant = (
            "book_slope_clamp"
            if model == "book_slope_clamp"
            else "timoshenko_section_clamp"
        )
        clamp_1 = cantilever_clamp_matrix(point_1, clamp_variant, scaled=False)
        clamp_2 = cantilever_clamp_matrix(point_2, clamp_variant, scaled=False)

        def transfer(omega: complex, point: RodPoint) -> NDArray[np.complex128]:
            return transfer_cache.get(
                "chapter2_state_corrected_timoshenko",
                omega,
                point,
                lambda: physical_state_transfer_matrix(omega, point),
            )

    elif model == "rectangular_eb_saint_venant":
        clamp_1 = eb_clamp_matrix(point_1)
        clamp_2 = eb_clamp_matrix(point_2)

        def transfer(omega: complex, point: RodPoint) -> NDArray[np.complex128]:
            return transfer_cache.get(
                "rectangular_orthotropic_eb_saint_venant",
                omega,
                point,
                lambda: eb_state_transfer_matrix(omega, point),
            )

    else:
        raise ValueError(f"unsupported Chapter-2 fast spectrum model: {model}")

    joint = joint_matrix_book(beta_rad)

    def raw_factory(omega: complex) -> NDArray[np.complex128]:
        counters.boundary_matrix_assemblies += 1
        block = np.zeros((12, 6), dtype=np.complex128)
        block[:6, :3] = transfer(omega, point_1) @ clamp_1
        block[6:, 3:] = transfer(omega, point_2) @ clamp_2
        return joint @ block

    def scaled_factory(omega: complex) -> NDArray[np.complex128]:
        return equilibrate_matrix(raw_factory(omega))[0]

    return scaled_factory, raw_factory


def _finalize_fast_roots(
    roots: Sequence[RootResult],
    raw_factory: Callable[[complex], NDArray[np.complex128]],
    counters: PerformanceCounters,
    *,
    scan_step_hz: float,
) -> SpectrumResult:
    ordered = sorted(roots, key=lambda root: root.frequency_hz)
    updated: list[RootResult] = []
    for index, root in enumerate(ordered):
        neighbors: list[float] = []
        if index:
            neighbors.append(root.frequency_hz - ordered[index - 1].frequency_hz)
        if index + 1 < len(ordered):
            neighbors.append(ordered[index + 1].frequency_hz - root.frequency_hz)
        updated.append(
            replace(
                root,
                min_neighbor_distance_hz=min(neighbors) if neighbors else math.inf,
            )
        )
    quality = [_root_quality_fields(root, raw_factory) for root in updated]
    counters.raw_quality_evaluations += len(updated)
    _validate_spectrum(updated, quality)
    return SpectrumResult(
        roots=updated,
        quality=quality,
        scaled_matrix_evaluations=0,
        raw_quality_matrix_evaluations=len(updated),
        runtime_seconds=0.0,
        scan_step_hz=scan_step_hz,
    )


def _solve_cached_global_spectrum(
    point_for_contract: RodPoint,
    scaled_factory: Callable[[complex], NDArray[np.complex128]],
    raw_factory: Callable[[complex], NDArray[np.complex128]],
    counters: PerformanceCounters,
    *,
    scan_step_hz: float = BASE_POINT_SCAN_STEP_HZ,
) -> SpectrumResult:
    builder = _CountingBoundaryBuilder(scaled_factory)
    started = time.perf_counter()
    roots = find_elastic_roots(
        point_for_contract,
        FORMULATION,  # type: ignore[arg-type]
        num_roots=GUARD_ROOT_COUNT,
        scan_step_hz=scan_step_hz,
        initial_max_hz=FAST_GLOBAL_INITIAL_MAX_HZ,
        max_hz=MAX_HZ,
        boundary_matrix_builder=builder,
    )
    counters.scaled_quality_evaluations += builder.evaluations
    result = _finalize_fast_roots(
        roots,
        raw_factory,
        counters,
        scan_step_hz=scan_step_hz,
    )
    result.scaled_matrix_evaluations = builder.evaluations
    result.runtime_seconds = time.perf_counter() - started
    return result


def _scan_fast_interval(
    point_for_contract: RodPoint,
    scaled_factory: Callable[[complex], NDArray[np.complex128]],
    interval: SearchInterval,
    counters: PerformanceCounters,
) -> list[RootResult]:
    width = interval.upper_hz - interval.lower_hz
    if width <= 0.0:
        raise ValueError("fast local interval is empty")
    predicted_gaps = np.diff(np.asarray(interval.predicted_hz, dtype=float))
    positive_gaps = predicted_gaps[predicted_gaps > 0.0]
    cluster_step = (
        float(np.min(positive_gaps)) / 12.0
        if positive_gaps.size
        else LOCAL_HINT_MAX_SCAN_STEP_HZ
    )
    scan_step_hz = max(
        1.0e-4,
        min(
            LOCAL_HINT_MAX_SCAN_STEP_HZ,
            width / max(30.0 * interval.expected_count, 30.0),
            cluster_step,
        ),
    )
    offset = 2.0 * math.pi * interval.lower_hz
    builder = _CountingBoundaryBuilder(
        lambda local_omega: scaled_factory(local_omega + offset)
    )
    roots = find_elastic_roots(
        point_for_contract,
        FORMULATION,  # type: ignore[arg-type]
        num_roots=interval.expected_count,
        scan_step_hz=scan_step_hz,
        initial_max_hz=width,
        max_hz=width,
        boundary_matrix_builder=builder,
    )
    counters.scaled_quality_evaluations += builder.evaluations
    return [
        replace(
            root,
            omega=root.omega + offset,
            frequency_hz=root.frequency_hz + interval.lower_hz,
            status=f"{root.status}_fast_local",
        )
        for root in roots
    ]


def _run_fast_family(
    preset: FigurePreset,
    beta_values: NDArray[np.float64],
    model: str,
    *,
    transfer_cache: ExactTransferLRU | None = None,
    counters: PerformanceCounters | None = None,
) -> FastSweepResult[SpectrumResult]:
    """Compute one sequential family without adding any physical equations."""

    shared_counters = counters or PerformanceCounters()
    cache = transfer_cache or ExactTransferLRU(
        FastSweepSettings().cache_maxsize, counters=shared_counters
    )
    if cache.counters is not shared_counters:
        raise ValueError("transfer cache and fast family must share counters")
    before = shared_counters.copy()
    point_1 = _point(preset, 1)
    point_2 = _point(preset, 2)
    beta_index = {float(value): index for index, value in enumerate(beta_values)}
    factory_cache: dict[
        float,
        tuple[
            Callable[[complex], NDArray[np.complex128]],
            Callable[[complex], NDArray[np.complex128]],
        ],
    ] = {}

    def factories(beta_deg: float):
        if beta_deg not in factory_cache:
            factory_cache[beta_deg] = _cached_factories(
                model,
                beta_deg,
                point_1,
                point_2,
                cache,
                shared_counters,
            )
        return factory_cache[beta_deg]

    def global_search(beta_deg: float) -> SpectrumResult:
        scaled_factory, raw_factory = factories(beta_deg)
        return _solve_cached_global_spectrum(
            point_1, scaled_factory, raw_factory, shared_counters
        )

    def scan_interval(beta_deg: float, interval: SearchInterval):
        scaled_factory, _ = factories(beta_deg)
        return _scan_fast_interval(
            point_1, scaled_factory, interval, shared_counters
        )

    def finalize(beta_deg: float, roots: Sequence[RootResult]) -> SpectrumResult:
        scaled_factory, raw_factory = factories(beta_deg)
        reconciled = list(roots)
        if beta_deg not in FastSweepSettings().mandatory_anchors_deg:
            shared_counters.global_inventory_checks += 1
            coarse = _solve_cached_global_spectrum(
                point_1,
                scaled_factory,
                raw_factory,
                shared_counters,
                scan_step_hz=SCAN_STEP_HZ,
            )
            candidates = list(coarse.roots) + reconciled
            candidates.sort(key=lambda root: root.frequency_hz)
            reconciled = []
            for candidate in candidates:
                if any(
                    abs(candidate.frequency_hz - existing.frequency_hz)
                    <= 2.0e-8
                    * max(
                        abs(candidate.frequency_hz),
                        abs(existing.frequency_hz),
                        1.0,
                    )
                    for existing in reconciled
                ):
                    continue
                reconciled.append(candidate)
            reconciled = reconciled[:GUARD_ROOT_COUNT]
        return _finalize_fast_roots(
            reconciled,
            raw_factory,
            shared_counters,
            scan_step_hz=LOCAL_HINT_MAX_SCAN_STEP_HZ,
        )

    def fallback_search(
        beta_deg: float,
        previous: SpectrumResult,
        older: SpectrumResult | None,
    ) -> SpectrumResult:
        scaled_factory, raw_factory = factories(beta_deg)
        global_result = _solve_cached_global_spectrum(
            point_1, scaled_factory, raw_factory, shared_counters
        )
        previous_frequencies = np.asarray(
            [root.frequency_hz for root in previous.roots], dtype=float
        )
        older_spectrum = older if older is not None else previous
        older_frequencies = np.asarray(
            [root.frequency_hz for root in older_spectrum.roots], dtype=float
        )
        current_index = beta_index[beta_deg]
        previous_beta_deg = float(beta_values[current_index - 1])
        predicted = (
            predict_sorted_frequencies(
                beta_deg,
                previous_beta_deg,
                previous_frequencies,
                older_beta_deg=float(beta_values[current_index - 2]),
                older_frequencies_hz=older_frequencies,
            )
            if current_index >= 2
            else previous_frequencies.copy()
        )
        relative_previous_gaps = np.diff(previous_frequencies) / np.maximum(
            previous_frequencies[:-1], previous_frequencies[1:]
        )
        predicted_denominator = np.maximum(
            np.abs(predicted[:-1]), np.abs(predicted[1:])
        )
        relative_predicted_gaps = np.diff(predicted) / np.maximum(
            predicted_denominator, 1.0
        )
        suspect_pairs = np.flatnonzero(
            np.logical_or(
                relative_previous_gaps < 5.0e-3,
                relative_predicted_gaps < 5.0e-3,
            )
        )
        hint_indices = sorted(
            {
                int(index)
                for pair in suspect_pairs
                for index in (pair, pair + 1)
            }
        )
        if not hint_indices:
            hint_indices = list(range(GUARD_ROOT_COUNT))
        recovered = _recover_spectrum_from_sorted_hints(
            point_1,
            scaled_factory,
            raw_factory,
            older_spectrum,
            previous,
            global_result,
            hint_indices=hint_indices,
        )
        shared_counters.scaled_quality_evaluations += (
            recovered.scaled_matrix_evaluations
        )
        shared_counters.raw_quality_evaluations += (
            recovered.raw_quality_matrix_evaluations
        )
        return recovered

    result = run_fast_beta_sweep(
        beta_values,
        global_search=global_search,
        scan_interval=scan_interval,
        finalize_local=finalize,
        spectrum_frequencies=lambda spectrum: [
            root.frequency_hz for root in spectrum.roots
        ],
        root_frequency=lambda root: root.frequency_hz,
        fallback_search=fallback_search,
        settings=FastSweepSettings(root_count=GUARD_ROOT_COUNT),
        counters=shared_counters,
    )
    result.counters = shared_counters.difference(before)
    return result


def _fast_family_fingerprint(
    preset: FigurePreset,
    beta_values: NDArray[np.float64],
    model: str,
    *,
    lambda_normalization: Mapping[str, Any],
) -> str:
    return stable_fingerprint(
        {
            "solver_version": FAST_SOLVER_VERSION,
            "settings": FastSweepSettings(),
            "model": model,
            "material": preset.material_factory(),
            "point_1": _point(preset, 1),
            "point_2": _point(preset, 2),
            "theta_1_deg": preset.theta_1_deg,
            "theta_2_deg": preset.theta_2_deg,
            "mu": preset.mu,
            "beta_grid_deg": beta_values.tolist(),
            "root_count": GUARD_ROOT_COUNT,
            "lambda_normalization": dict(lambda_normalization),
        }
    )


def _fast_family_checkpoint_rows(
    result: FastSweepResult[SpectrumResult], family_id: str
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for beta_deg, spectrum in sorted(result.spectra.items()):
        for mode, (root, quality) in enumerate(
            zip(spectrum.roots, spectrum.quality), start=1
        ):
            rows.append(
                {
                    "family_id": family_id,
                    "beta_deg": beta_deg,
                    "mode": mode,
                    "omega_real": float(root.omega.real),
                    "omega_imag": float(root.omega.imag),
                    "frequency_hz": root.frequency_hz,
                    "determinant_residual": root.determinant_residual,
                    "raw_determinant_abs": root.raw_determinant_abs,
                    "sigma_min": root.sigma_min,
                    "sigma_max": root.sigma_max,
                    "relative_singular_residual": root.relative_singular_residual,
                    "min_neighbor_distance_hz": root.min_neighbor_distance_hz,
                    "refinements": root.refinements,
                    "root_status": root.status,
                    "quality_json": json.dumps(
                        quality, ensure_ascii=False, sort_keys=True
                    ),
                    "data_origin": result.data_origins[beta_deg],
                    "fallback_reason": result.fallback_reasons.get(beta_deg, ""),
                    "anchor_relative_error": result.anchor_relative_errors.get(
                        beta_deg, ""
                    ),
                }
            )
    return rows


def _fast_family_from_checkpoint(
    rows: Sequence[Mapping[str, str]], manifest: Mapping[str, Any]
) -> FastSweepResult[SpectrumResult]:
    grouped: dict[float, list[Mapping[str, str]]] = {}
    for row in rows:
        grouped.setdefault(float(row["beta_deg"]), []).append(row)
    spectra: dict[float, SpectrumResult] = {}
    origins: dict[float, str] = {}
    fallbacks: dict[float, str] = {}
    anchors: dict[float, float] = {}
    for beta_deg, group in sorted(grouped.items()):
        ordered = sorted(group, key=lambda row: int(row["mode"]))
        roots: list[RootResult] = []
        quality: list[dict[str, Any]] = []
        for row in ordered:
            roots.append(
                RootResult(
                    omega=complex(
                        float(row["omega_real"]), float(row["omega_imag"])
                    ),
                    frequency_hz=float(row["frequency_hz"]),
                    determinant_residual=float(row["determinant_residual"]),
                    raw_determinant_abs=float(row["raw_determinant_abs"]),
                    sigma_min=float(row["sigma_min"]),
                    sigma_max=float(row["sigma_max"]),
                    relative_singular_residual=float(
                        row["relative_singular_residual"]
                    ),
                    min_neighbor_distance_hz=float(
                        row["min_neighbor_distance_hz"]
                    ),
                    refinements=int(row["refinements"]),
                    status=row["root_status"],
                )
            )
            quality.append(json.loads(row["quality_json"]))
        _validate_spectrum(roots, quality)
        spectra[beta_deg] = SpectrumResult(
            roots=roots,
            quality=quality,
            scaled_matrix_evaluations=0,
            raw_quality_matrix_evaluations=len(roots),
            runtime_seconds=0.0,
            scan_step_hz=LOCAL_HINT_MAX_SCAN_STEP_HZ,
        )
        origins[beta_deg] = ordered[0]["data_origin"]
        if ordered[0].get("fallback_reason"):
            fallbacks[beta_deg] = ordered[0]["fallback_reason"]
        if ordered[0].get("anchor_relative_error") not in (None, ""):
            anchors[beta_deg] = float(ordered[0]["anchor_relative_error"])
    counter_values = {
        name: manifest.get("performance_counters", {}).get(name, 0)
        for name in PerformanceCounters.__dataclass_fields__
    }
    counters = PerformanceCounters(**counter_values)
    return FastSweepResult(
        spectra=spectra,
        data_origins=origins,
        fallback_reasons=fallbacks,
        anchor_relative_errors=anchors,
        counters=counters,
        runtime_seconds=float(manifest.get("runtime_seconds", 0.0)),
        status="PASS",
    )


def _run_or_resume_fast_family(
    preset: FigurePreset,
    beta_values: NDArray[np.float64],
    model: str,
    *,
    family_id: str,
    output_dir: Path,
    resume: bool,
    transfer_cache: ExactTransferLRU,
    counters: PerformanceCounters,
    lambda_normalization: Mapping[str, Any],
) -> tuple[FastSweepResult[SpectrumResult], bool]:
    checkpoint_dir = output_dir / FAST_FAMILY_DIRNAME
    fingerprint = _fast_family_fingerprint(
        preset,
        beta_values,
        model,
        lambda_normalization=lambda_normalization,
    )
    expected_rows = len(beta_values) * GUARD_ROOT_COUNT
    if resume:
        loaded = load_family_checkpoint(
            checkpoint_dir,
            family_id,
            expected_fingerprint=fingerprint,
            expected_row_count=expected_rows,
        )
        if loaded is not None:
            rows, manifest = loaded
            return _fast_family_from_checkpoint(rows, manifest), True
    result = _run_fast_family(
        preset,
        beta_values,
        model,
        transfer_cache=transfer_cache,
        counters=counters,
    )
    rows = _fast_family_checkpoint_rows(result, family_id)
    write_family_checkpoint(
        checkpoint_dir,
        family_id,
        rows,
        fingerprint=fingerprint,
        metadata={
            "solver_version": FAST_SOLVER_VERSION,
            "model": model,
            "material": preset.material_name,
            "runtime_seconds": result.runtime_seconds,
            "performance_counters": result.counters.to_dict(),
            "fallback_reasons": {
                str(beta): reason
                for beta, reason in sorted(result.fallback_reasons.items())
            },
            "anchor_relative_errors": {
                str(beta): value
                for beta, value in sorted(result.anchor_relative_errors.items())
            },
            "lambda_normalization": dict(lambda_normalization),
        },
    )
    return result, False


def _fast_result_as_sweep(
    result: FastSweepResult[SpectrumResult],
) -> SweepResult:
    """Expose fast results to the frozen Figure-2--4 CSV row builders."""

    return SweepResult(
        spectra=result.spectra,
        runtime_seconds=result.runtime_seconds,
        scaled_matrix_evaluations=result.counters.scaled_quality_evaluations,
        raw_quality_matrix_evaluations=result.counters.raw_quality_evaluations,
        completeness_refined_beta_count=0,
        completeness_refinement_attempts=0,
    )


def _selected_sweep(
    preset: FigurePreset,
    beta_values: NDArray[np.float64],
    model: str,
    *,
    solver_mode: str,
    beta_zero_result: SpectrumResult | None = None,
    transfer_cache: ExactTransferLRU | None = None,
    counters: PerformanceCounters | None = None,
) -> SweepResult:
    """Select the preserved legacy path or the validated fast coordinator."""

    if solver_mode == "legacy":
        return _sweep(
            preset,
            beta_values,
            model,
            beta_zero_result=beta_zero_result,
        )
    if solver_mode != "fast":
        raise ValueError(f"unsupported solver mode: {solver_mode}")
    return _fast_result_as_sweep(
        _run_fast_family(
            preset,
            beta_values,
            model,
            transfer_cache=transfer_cache,
            counters=counters,
        )
    )


def _sweep(
    preset: FigurePreset,
    beta_values: NDArray[np.float64],
    model: str,
    *,
    beta_zero_result: SpectrumResult | None = None,
) -> SweepResult:
    point_1 = _point(preset, 1)
    point_2 = _point(preset, 2)
    spectra: dict[float, SpectrumResult] = {}
    refined_beta_count = 0
    refinement_attempt_count = 0
    started = time.perf_counter()
    for index, beta_value in enumerate(beta_values):
        beta_deg = float(beta_value)
        if beta_deg == 0.0 and beta_zero_result is not None:
            result = beta_zero_result
        else:
            scaled_factory, raw_factory = _factories(
                model, beta_deg, point_1, point_2
            )
            result = _solve_spectrum(
                point_1,
                scaled_factory,
                raw_factory,
                scan_step_hz=(BASE_POINT_SCAN_STEP_HZ if beta_deg == 0.0 else SCAN_STEP_HZ),
            )
            if spectra:
                previous = spectra[float(beta_values[index - 1])]
                previous_frequency = np.asarray(
                    [root.frequency_hz for root in previous.roots], dtype=float
                )
                close_pairs = np.flatnonzero(
                    np.diff(previous_frequency) <= CLOSE_ROOT_HINT_GAP_HZ
                )
                close_indices = sorted(
                    {
                        int(mode_index)
                        for pair_index in close_pairs
                        for mode_index in (pair_index, pair_index + 1)
                    }
                )
                if close_indices:
                    previous_previous = (
                        spectra[float(beta_values[index - 2])]
                        if index >= 2
                        else previous
                    )
                    verified = _recover_spectrum_from_sorted_hints(
                        point_1,
                        scaled_factory,
                        raw_factory,
                        previous_previous,
                        previous,
                        result,
                        hint_indices=close_indices,
                    )
                    result = _retain_final_attempt_with_accounting(
                        [result, verified]
                    )
                jump, jump_mode = _spectral_step_difference(previous, result)
                if jump > SPECTRAL_STEP_JUMP_TRIGGER:
                    print(
                        f"{model}: completeness refinement at beta={beta_deg:g} deg "
                        f"(initial jump={jump:.3e}, mode={jump_mode})",
                        flush=True,
                    )
                    attempts = [result]
                    for refined_step in COMPLETENESS_REFINEMENT_STEPS_HZ:
                        attempts.append(
                            _solve_spectrum(
                                point_1,
                                scaled_factory,
                                raw_factory,
                                scan_step_hz=refined_step,
                            )
                        )
                        result = attempts[-1]
                        jump, jump_mode = _spectral_step_difference(previous, result)
                        print(
                            f"{model}: scan_step={refined_step:g} Hz -> "
                            f"jump={jump:.3e}, mode={jump_mode}",
                            flush=True,
                        )
                        if jump <= SPECTRAL_STEP_JUMP_TRIGGER:
                            break
                    if jump > SPECTRAL_STEP_JUMP_TRIGGER and len(spectra) >= 2:
                        previous_previous = spectra[float(beta_values[index - 2])]
                        attempts.append(
                            _recover_spectrum_from_sorted_hints(
                                point_1,
                                scaled_factory,
                                raw_factory,
                                previous_previous,
                                previous,
                                result,
                            )
                        )
                        result = attempts[-1]
                        jump, jump_mode = _spectral_step_difference(previous, result)
                        print(
                            f"{model}: sorted-hint windows -> jump={jump:.3e}, "
                            f"mode={jump_mode}",
                            flush=True,
                        )
                    result = _retain_final_attempt_with_accounting(attempts)
                    if jump > SPECTRAL_STEP_JUMP_TRIGGER:
                        raise RuntimeError(
                            "suspected missing root after completeness refinement: "
                            f"model={model}, beta={beta_deg:g} deg, mode={jump_mode}, "
                            f"relative sorted-spectrum step={jump:.6e}, "
                            "previous_hz="
                            f"{[root.frequency_hz for root in previous.roots]}, "
                            "current_hz="
                            f"{[root.frequency_hz for root in result.roots]}"
                        )
            if result.completeness_refinement_attempts:
                refined_beta_count += 1
                refinement_attempt_count += result.completeness_refinement_attempts
        spectra[beta_deg] = result
        if index == 0 or index + 1 == len(beta_values) or index % max(1, len(beta_values) // 10) == 0:
            print(
                f"{model}: beta={beta_deg:g} deg ({index + 1}/{len(beta_values)})",
                flush=True,
            )
    return SweepResult(
        spectra=spectra,
        runtime_seconds=time.perf_counter() - started,
        scaled_matrix_evaluations=sum(
            item.scaled_matrix_evaluations for item in spectra.values()
        ),
        raw_quality_matrix_evaluations=sum(
            item.raw_quality_matrix_evaluations for item in spectra.values()
        ),
        completeness_refined_beta_count=refined_beta_count,
        completeness_refinement_attempts=refinement_attempt_count,
    )


def section_clamp_straight_reference(
    preset: FigurePreset = FIGURES_3_4_PRESET,
) -> SectionReferenceResult:
    """Compare seven beta-zero coupled roots with an independent 2L rod."""

    point_1 = _point(preset, 1)
    point_2 = _point(preset, 2)
    point_total = _point(preset, 1, total_length=True)
    coupled_scaled, coupled_raw = _factories(
        "timoshenko_section_clamp", 0.0, point_1, point_2
    )
    coupled = _solve_spectrum(
        point_1,
        coupled_scaled,
        coupled_raw,
        scan_step_hz=BASE_POINT_SCAN_STEP_HZ,
    )
    straight = _solve_spectrum(
        point_total,
        lambda omega: _straight_section_clamp_boundary_matrix(
            omega, point_total, scaled=True
        ),
        lambda omega: _straight_section_clamp_boundary_matrix_raw(
            omega, point_total
        ),
        scan_step_hz=BASE_POINT_SCAN_STEP_HZ,
    )
    left = np.asarray([root.frequency_hz for root in coupled.roots])
    right = np.asarray([root.frequency_hz for root in straight.roots])
    maximum = float(np.max(np.abs(left - right) / np.maximum(np.abs(right), 1.0)))
    status = "PASS" if maximum <= 1.0e-8 else "FAIL_SECTION_CLAMP_REFERENCE"
    return SectionReferenceResult(
        preset_label=preset.material_name,
        coupled=coupled,
        straight=straight,
        maximum_relative_frequency_difference=maximum,
        status=status,
    )


def _lambda_for_root(root: RootResult, preset: FigurePreset, point: RodPoint) -> float:
    elastic_ex = float(np.real(point.properties.Ex))
    return float(
        lambda_from_omega(
            float(root.omega.real),
            rho_kg_m3=point.material.rho,
            a_m=preset.a_m,
            b_m=preset.b_m,
            length_1_m=preset.length_1_m,
            length_2_m=preset.length_2_m,
            elastic_ex_pa=elastic_ex,
        )
    )


def _quality_columns(prefix: str, item: Mapping[str, Any]) -> dict[str, Any]:
    return {
        f"{prefix}_quality_status": item["quality_status"],
        f"{prefix}_quality_basis": item["quality_basis"],
        f"{prefix}_root_status": item["root_status"],
        f"{prefix}_accepted_determinant_residual": item[
            "accepted_determinant_residual"
        ],
        f"{prefix}_accepted_relative_singular_residual": item[
            "accepted_relative_singular_residual"
        ],
        f"{prefix}_scaled_determinant_residual": item[
            "scaled_determinant_residual"
        ],
        f"{prefix}_scaled_relative_singular_residual": item[
            "scaled_relative_singular_residual"
        ],
        f"{prefix}_physical_raw_normalized_determinant_residual": item[
            "physical_raw_normalized_determinant_residual"
        ],
        f"{prefix}_physical_raw_relative_singular_residual": item[
            "physical_raw_relative_singular_residual"
        ],
    }


def _figure_2_rows(
    beta_values: NDArray[np.float64],
    book: SweepResult,
    section: SweepResult,
) -> list[dict[str, Any]]:
    point = _point(FIGURE_2_PRESET, 1)
    rows: list[dict[str, Any]] = []
    for beta_value in beta_values:
        beta_deg = float(beta_value)
        book_result = book.spectra[beta_deg]
        section_result = section.spectra[beta_deg]
        for mode in range(1, GUARD_ROOT_COUNT + 1):
            left = book_result.roots[mode - 1]
            right = section_result.roots[mode - 1]
            left_lambda = _lambda_for_root(left, FIGURE_2_PRESET, point)
            right_lambda = _lambda_for_root(right, FIGURE_2_PRESET, point)
            rows.append(
                {
                    "beta_deg": beta_deg,
                    "mode": mode,
                    "root_role": "plotted" if mode <= PLOTTED_ROOT_COUNT else "guard",
                    "book_slope_frequency_hz": left.frequency_hz,
                    "book_slope_lambda": left_lambda,
                    "book_slope_scan_step_hz": book_result.scan_step_hz,
                    "book_slope_completeness_refinement_attempts": book_result.completeness_refinement_attempts,
                    **_quality_columns("book_slope", book_result.quality[mode - 1]),
                    "section_clamp_frequency_hz": right.frequency_hz,
                    "section_clamp_lambda": right_lambda,
                    "section_clamp_scan_step_hz": section_result.scan_step_hz,
                    "section_clamp_completeness_refinement_attempts": section_result.completeness_refinement_attempts,
                    **_quality_columns("section_clamp", section_result.quality[mode - 1]),
                    "relative_clamp_difference": abs(right_lambda - left_lambda)
                    / left_lambda,
                    "relative_clamp_frequency_difference": abs(
                        right.frequency_hz - left.frequency_hz
                    )
                    / left.frequency_hz,
                }
            )
    return rows


def _figure_3_4_rows(
    figure: int,
    beta_values: NDArray[np.float64],
    timoshenko: SweepResult,
    eb: SweepResult,
    clamp_type: str,
) -> list[dict[str, Any]]:
    point = _point(FIGURES_3_4_PRESET, 1)
    rows: list[dict[str, Any]] = []
    for beta_value in beta_values:
        beta_deg = float(beta_value)
        timo_result = timoshenko.spectra[beta_deg]
        eb_result = eb.spectra[beta_deg]
        for mode in range(1, GUARD_ROOT_COUNT + 1):
            timo_root = timo_result.roots[mode - 1]
            eb_root = eb_result.roots[mode - 1]
            timo_lambda = _lambda_for_root(timo_root, FIGURES_3_4_PRESET, point)
            eb_lambda = _lambda_for_root(eb_root, FIGURES_3_4_PRESET, point)
            rows.append(
                {
                    "figure": figure,
                    "clamp_type": clamp_type,
                    "beta_deg": beta_deg,
                    "mode": mode,
                    "root_role": "plotted" if mode <= PLOTTED_ROOT_COUNT else "guard",
                    "timoshenko_frequency_hz": timo_root.frequency_hz,
                    "timoshenko_lambda": timo_lambda,
                    "timoshenko_scan_step_hz": timo_result.scan_step_hz,
                    "timoshenko_completeness_refinement_attempts": timo_result.completeness_refinement_attempts,
                    **_quality_columns("timoshenko", timo_result.quality[mode - 1]),
                    "eb_frequency_hz": eb_root.frequency_hz,
                    "eb_lambda": eb_lambda,
                    "eb_scan_step_hz": eb_result.scan_step_hz,
                    "eb_completeness_refinement_attempts": eb_result.completeness_refinement_attempts,
                    **_quality_columns("eb", eb_result.quality[mode - 1]),
                    "relative_theory_difference": abs(eb_lambda - timo_lambda)
                    / timo_lambda,
                    "relative_theory_frequency_difference": abs(
                        eb_root.frequency_hz - timo_root.frequency_hz
                    )
                    / timo_root.frequency_hz,
                }
            )
    return rows


def _fixed_reference_lambda_for_root(
    root: RootResult, preset: FigurePreset
) -> float:
    reference_point = _point(FIGURES_3_4_PRESET, 1)
    return float(
        lambda_from_omega(
            float(root.omega.real),
            rho_kg_m3=preset.material_factory().rho,
            a_m=REFERENCE_A_M,
            b_m=REFERENCE_B_M,
            length_1_m=preset.length_1_m,
            length_2_m=preset.length_2_m,
            elastic_ex_pa=float(reference_point.properties.Ex.real),
        )
    )


def _saved_quality(row: Mapping[str, Any], prefix: str) -> dict[str, Any]:
    return {
        "quality_status": row[f"{prefix}_quality_status"],
        "quality_basis": row[f"{prefix}_quality_basis"],
        "root_status": row[f"{prefix}_root_status"],
        "accepted_determinant_residual": float(
            row[f"{prefix}_accepted_determinant_residual"]
        ),
        "accepted_relative_singular_residual": float(
            row[f"{prefix}_accepted_relative_singular_residual"]
        ),
        "scaled_determinant_residual": float(
            row[f"{prefix}_scaled_determinant_residual"]
        ),
        "scaled_relative_singular_residual": float(
            row[f"{prefix}_scaled_relative_singular_residual"]
        ),
        "physical_raw_normalized_determinant_residual": float(
            row[f"{prefix}_physical_raw_normalized_determinant_residual"]
        ),
        "physical_raw_relative_singular_residual": float(
            row[f"{prefix}_physical_raw_relative_singular_residual"]
        ),
    }


def _extended_figure_rows(
    figure: int,
    preset: FigurePreset,
    beta_values: NDArray[np.float64],
    left: FastSweepResult[SpectrumResult],
    *,
    comparison_type: str,
    left_model: str,
    left_theta_deg: float,
    right_model: str,
    right_theta_deg: float,
    right_fast: FastSweepResult[SpectrumResult] | None = None,
    reused_figure_3_rows: Sequence[Mapping[str, Any]] | None = None,
    reused_prefix: str | None = None,
) -> list[dict[str, Any]]:
    if (right_fast is None) == (reused_figure_3_rows is None):
        raise ValueError("exactly one right-side source must be supplied")
    reused_map = (
        {
            (float(row["beta_deg"]), int(row["mode"])): row
            for row in reused_figure_3_rows or []
        }
        if reused_figure_3_rows is not None
        else {}
    )
    rows: list[dict[str, Any]] = []
    for beta_value in beta_values:
        beta_deg = float(beta_value)
        left_spectrum = left.spectra[beta_deg]
        for mode in range(1, GUARD_ROOT_COUNT + 1):
            left_root = left_spectrum.roots[mode - 1]
            left_lambda = _fixed_reference_lambda_for_root(left_root, preset)
            if right_fast is not None:
                right_spectrum = right_fast.spectra[beta_deg]
                right_root = right_spectrum.roots[mode - 1]
                right_frequency = right_root.frequency_hz
                right_lambda = _fixed_reference_lambda_for_root(right_root, preset)
                right_quality = right_spectrum.quality[mode - 1]
                right_origin = right_fast.data_origins[beta_deg]
            else:
                source = reused_map[(beta_deg, mode)]
                if reused_prefix is None:
                    raise ValueError("reused source prefix is required")
                right_frequency = float(source[f"{reused_prefix}_frequency_hz"])
                right_lambda = float(source[f"{reused_prefix}_lambda"])
                right_quality = _saved_quality(source, reused_prefix)
                right_origin = "reused_figure_03"
            left_origin = left.data_origins[beta_deg]
            data_origin = (
                "global_fallback"
                if left_origin == "global_fallback"
                or right_origin == "global_fallback"
                else "reused_figure_03"
                if right_origin == "reused_figure_03"
                else "new_fast_solve"
            )
            absolute = abs(left_lambda - right_lambda)
            rows.append(
                {
                    "figure": figure,
                    "comparison_type": comparison_type,
                    "beta_deg": beta_deg,
                    "mode": mode,
                    "root_role": (
                        "plotted" if mode <= PLOTTED_ROOT_COUNT else "guard"
                    ),
                    "left_model": left_model,
                    "left_theta_deg": left_theta_deg,
                    "left_frequency_hz": left_root.frequency_hz,
                    "left_lambda": left_lambda,
                    **_quality_columns(
                        "left", left_spectrum.quality[mode - 1]
                    ),
                    "right_model": right_model,
                    "right_theta_deg": right_theta_deg,
                    "right_frequency_hz": right_frequency,
                    "right_lambda": right_lambda,
                    **_quality_columns("right", right_quality),
                    "absolute_lambda_difference": absolute,
                    "relative_lambda_difference": absolute
                    / max(abs(right_lambda), np.finfo(float).tiny),
                    "data_origin": data_origin,
                    "left_data_origin": left_origin,
                    "right_data_origin": right_origin,
                }
            )
    return rows


def _extended_rows_as_numbers(
    rows: Sequence[Mapping[str, Any]], figure: int
) -> list[dict[str, Any]]:
    converted: list[dict[str, Any]] = []
    float_fields = {
        "beta_deg",
        "left_theta_deg",
        "left_frequency_hz",
        "left_lambda",
        "right_theta_deg",
        "right_frequency_hz",
        "right_lambda",
        "absolute_lambda_difference",
        "relative_lambda_difference",
    }
    for row in rows:
        item = dict(row)
        item["figure"] = int(item["figure"])
        item["mode"] = int(item["mode"])
        if item["figure"] != figure:
            raise ValueError(f"saved row does not belong to Figure {figure}")
        for field in float_fields:
            item[field] = float(item[field])
        converted.append(item)
    return converted


def _validate_saved_extended_rows(
    rows: Sequence[Mapping[str, Any]],
    beta_values: NDArray[np.float64],
    figure: int,
) -> None:
    if len(rows) != len(beta_values) * GUARD_ROOT_COUNT:
        raise RuntimeError(f"Figure {figure} does not contain 181x7 complete rows")
    for beta_deg in beta_values:
        selected = sorted(
            (row for row in rows if float(row["beta_deg"]) == float(beta_deg)),
            key=lambda row: int(row["mode"]),
        )
        if [int(row["mode"]) for row in selected] != list(
            range(1, GUARD_ROOT_COUNT + 1)
        ):
            raise RuntimeError(f"Figure {figure} has missing/duplicate modes")
        for prefix in ("left", "right"):
            frequencies = np.asarray(
                [float(row[f"{prefix}_frequency_hz"]) for row in selected]
            )
            if np.any(~np.isfinite(frequencies)) or np.any(frequencies <= 0.0):
                raise RuntimeError(f"Figure {figure} contains invalid roots")
            if np.any(np.diff(frequencies) <= 0.0):
                raise RuntimeError(f"Figure {figure} contains unsorted roots")
            if any(row[f"{prefix}_quality_status"] != "PASS" for row in selected):
                raise RuntimeError(f"Figure {figure} contains rejected roots")


def _minimum_relative_neighbor_gap(
    rows: Sequence[Mapping[str, Any]], frequency_keys: Sequence[str]
) -> dict[str, Any]:
    minimum = math.inf
    location: dict[str, Any] = {}
    grouped: dict[float, list[Mapping[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(float(row["beta_deg"]), []).append(row)
    for beta_deg, group in grouped.items():
        ordered = sorted(group, key=lambda row: int(row["mode"]))
        for key in frequency_keys:
            frequencies = np.asarray([float(row[key]) for row in ordered])
            gaps = np.diff(frequencies) / np.maximum(
                frequencies[:-1], frequencies[1:]
            )
            index = int(np.argmin(gaps))
            if float(gaps[index]) < minimum:
                minimum = float(gaps[index])
                location = {
                    "minimum_relative_neighbor_gap": minimum,
                    "beta_deg": beta_deg,
                    "pair": [index + 1, index + 2],
                    "model_side": "left" if key.startswith("left") else "right",
                }
    return location


def _extended_figure_performance(
    figure: int, metadata: Mapping[str, Any]
) -> dict[str, Any]:
    family_ids = {
        5: ("figure05_timoshenko", "figure05_eb"),
        6: ("figure06_timoshenko", "figure06_eb"),
        7: ("figure07_timoshenko_theta5",),
        8: ("figure08_timoshenko_theta15",),
    }[figure]
    families = metadata.get("families", {})
    selected = [families[name] for name in family_ids if name in families]
    counter_names = tuple(PerformanceCounters.__dataclass_fields__)
    counters: dict[str, int | float] = {}
    for name in counter_names:
        counters[name] = sum(
            float(family.get("performance_counters", {}).get(name, 0))
            for family in selected
        )
        if name not in ("family_runtime_s", "total_scientific_runtime_s"):
            counters[name] = int(counters[name])
    requests = int(counters.get("transfer_cache_hits", 0)) + int(
        counters.get("transfer_cache_misses", 0)
    )
    counters["transfer_cache_requests"] = requests
    counters["transfer_cache_hit_rate"] = (
        int(counters.get("transfer_cache_hits", 0)) / requests if requests else 0.0
    )
    fallback_reasons: dict[str, str] = {}
    for family_id in family_ids:
        family = families.get(family_id, {})
        for beta_deg, reason in family.get("fallback_reasons", {}).items():
            fallback_reasons[f"{family_id}@{beta_deg}"] = str(reason)
    return {
        "scientific_runtime_s": sum(
            float(family.get("runtime_seconds", 0.0)) for family in selected
        ),
        "global_fallback_count": sum(
            int(family.get("global_fallback_count", 0)) for family in selected
        ),
        "fallback_reasons": fallback_reasons,
        "performance_counters": counters,
    }


def _write_csv(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    data = [dict(row) for row in rows]
    if not data:
        raise ValueError(f"refusing to write empty CSV: {path}")
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


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _ensure_canonical_figure_1_data(*, reuse_data: bool) -> None:
    required = (
        CANONICAL_FIGURE_1_ROOTS,
        CANONICAL_FIGURE_1_DIGITIZED,
        CANONICAL_FIGURE_1_COMPARISON,
        CANONICAL_FIGURE_1_REPORT,
    )
    missing = [path for path in required if not path.is_file()]
    if not missing:
        report = CANONICAL_FIGURE_1_REPORT.read_text(encoding="utf-8")
        if "PASS_WITHIN_GRAPH_RESOLUTION" not in report:
            raise RuntimeError("canonical Figure-2.2 evidence is not verified PASS data")
        return
    if reuse_data:
        raise FileNotFoundError(
            "--reuse-data forbids restoring missing canonical Figure-2.2 data: "
            + ", ".join(str(path) for path in missing)
        )
    from scripts.analysis.anisotropic_rods.reproduce_yartsev_fig_2_2 import (
        main as reproduce_figure_2_2,
    )

    status = reproduce_figure_2_2([])
    if status != 0:
        raise RuntimeError(f"canonical Figure-2.2 workflow failed with exit code {status}")
    _ensure_canonical_figure_1_data(reuse_data=True)


def _build_figure_1_rows(*, reuse_data: bool) -> list[dict[str, Any]]:
    _ensure_canonical_figure_1_data(reuse_data=reuse_data)
    roots = _read_csv(CANONICAL_FIGURE_1_ROOTS)
    source_rows = _read_csv(CANONICAL_FIGURE_1_DIGITIZED)
    comparison_rows = _read_csv(CANONICAL_FIGURE_1_COMPARISON)
    if len(source_rows) != 98 or len(comparison_rows) != 98:
        raise RuntimeError("canonical Figure-2.2 digitized evidence must contain 98 rows")
    if not all(row["within_estimated_graphical_resolution"] == "True" for row in comparison_rows):
        raise RuntimeError("canonical Figure-2.2 comparison is not within graph resolution")
    digitized: dict[tuple[float, int, str], dict[str, str]] = {
        (float(row["theta_deg"]), int(row["book_curve_mode_label"]), row["quantity"]): row
        for row in source_rows
    }
    selected = [
        row
        for row in roots
        if row["material_mode"] == "book_complex"
        and row["equation_variant"] == "state_corrected"
        and 1 <= int(row["tracked_mode_if_available"]) <= 7
    ]
    if len(selected) != 91 * 7:
        raise RuntimeError("canonical Figure-2.2 corrected complex root grid is incomplete")
    result: list[dict[str, Any]] = []
    for row in selected:
        theta = float(row["theta_deg"])
        mode = int(row["tracked_mode_if_available"])
        frequency_source = digitized.get((theta, mode, "frequency"))
        loss_source = digitized.get((theta, mode, "eta_times_100"))
        result.append(
            {
                "theta_deg": theta,
                "book_curve_mode_label": mode,
                "computed_sorted_mode": int(row["sorted_mode"]),
                "calculated_frequency_khz": float(row["frequency_hz"]) / 1000.0,
                "calculated_modal_loss_factor_times_100": 100.0 * float(row["eta_exact"]),
                "digitized_frequency_khz": (
                    float(frequency_source["digitized_value"])
                    if frequency_source is not None
                    else ""
                ),
                "digitized_frequency_uncertainty_khz": (
                    float(frequency_source["estimated_graphical_uncertainty"])
                    if frequency_source is not None
                    else ""
                ),
                "digitized_modal_loss_factor_times_100": (
                    float(loss_source["digitized_value"])
                    if loss_source is not None
                    else ""
                ),
                "digitized_loss_uncertainty_times_100": (
                    float(loss_source["estimated_graphical_uncertainty"])
                    if loss_source is not None
                    else ""
                ),
                "root_status": row["root_status"],
                "data_source": "canonical verified Figure-2.2 reproduction CSV",
            }
        )
    return sorted(
        result, key=lambda item: (float(item["theta_deg"]), int(item["book_curve_mode_label"]))
    )


def _float_rows(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for source in rows:
        row = dict(source)
        for field in fields:
            if field in row and row[field] != "":
                row[field] = float(row[field])
        if "mode" in row:
            row["mode"] = int(float(row["mode"]))
        if "book_curve_mode_label" in row:
            row["book_curve_mode_label"] = int(float(row["book_curve_mode_label"]))
        result.append(row)
    return result


def _comparison_rows_as_numbers(
    rows: Sequence[Mapping[str, Any]], figure: int
) -> list[dict[str, Any]]:
    fields = [
        "beta_deg",
        "book_slope_frequency_hz",
        "book_slope_lambda",
        "section_clamp_frequency_hz",
        "section_clamp_lambda",
        "relative_clamp_difference",
        "timoshenko_frequency_hz",
        "timoshenko_lambda",
        "eb_frequency_hz",
        "eb_lambda",
        "relative_theory_difference",
    ]
    converted = _float_rows(rows, fields)
    expected = "relative_clamp_difference" if figure == 2 else "relative_theory_difference"
    if not converted or expected not in converted[0]:
        raise RuntimeError(f"Figure {figure} CSV has an incompatible schema")
    return converted


def _validate_saved_comparison_rows(
    rows: Sequence[Mapping[str, Any]], beta_values: NDArray[np.float64], figure: int
) -> None:
    for beta in beta_values:
        subset = [row for row in rows if float(row["beta_deg"]) == float(beta)]
        modes = sorted(int(row["mode"]) for row in subset)
        if modes != list(range(1, GUARD_ROOT_COUNT + 1)):
            raise RuntimeError(
                f"Figure {figure}: incomplete/duplicate saved roots at beta={beta:g}"
            )
        prefixes = (
            ("book_slope", "section_clamp")
            if figure == 2
            else ("timoshenko", "eb")
        )
        for row in subset:
            for prefix in prefixes:
                if row.get(f"{prefix}_quality_status") != "PASS":
                    raise RuntimeError(
                        f"Figure {figure}: saved root quality failed at beta={beta:g}"
                    )


def _style_context() -> dict[str, Any]:
    return {
        "font.family": FONT_FAMILY,
        "font.size": 10.0,
        "axes.labelsize": 10.0,
        "xtick.labelsize": 9.0,
        "ytick.labelsize": 9.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    }


def create_figure_1(rows: Sequence[Mapping[str, Any]]) -> plt.Figure:
    """Create the four-panel verified Figure-2.2 reproduction without legends."""

    numeric = _float_rows(
        rows,
        (
            "theta_deg",
            "calculated_frequency_khz",
            "calculated_modal_loss_factor_times_100",
            "digitized_frequency_khz",
            "digitized_modal_loss_factor_times_100",
        ),
    )
    with plt.rc_context(_style_context()):
        fig, axes = plt.subplots(
            2, 2, figsize=FIGURE_1_SIZE_IN, constrained_layout=True
        )
        frequency_axis = axes[0, 0]
        groups = ((1, 3, 5), (2, 4), (6, 7))
        loss_axes = (axes[0, 1], axes[1, 0], axes[1, 1])
        panel_labels = ("a)", "b) 1, 3, 5", "c) 2, 4", "d) 6, 7")
        for mode in range(1, 8):
            selected = sorted(
                [row for row in numeric if int(row["book_curve_mode_label"]) == mode],
                key=lambda row: float(row["theta_deg"]),
            )
            theta = np.asarray([float(row["theta_deg"]) for row in selected])
            frequency = np.asarray(
                [float(row["calculated_frequency_khz"]) for row in selected]
            )
            frequency_axis.plot(
                theta, frequency, color=FIGURE_1_COLORS[mode - 1], linewidth=1.6
            )
            source = [row for row in selected if row.get("digitized_frequency_khz", "") != ""]
            frequency_axis.plot(
                [float(row["theta_deg"]) for row in source],
                [float(row["digitized_frequency_khz"]) for row in source],
                linestyle="none",
                marker="o",
                markersize=3.2,
                markerfacecolor="white",
                markeredgewidth=0.8,
                color=FIGURE_1_COLORS[mode - 1],
            )
        for axis, group in zip(loss_axes, groups):
            for mode in group:
                selected = sorted(
                    [row for row in numeric if int(row["book_curve_mode_label"]) == mode],
                    key=lambda row: float(row["theta_deg"]),
                )
                axis.plot(
                    [float(row["theta_deg"]) for row in selected],
                    [float(row["calculated_modal_loss_factor_times_100"]) for row in selected],
                    color=FIGURE_1_COLORS[mode - 1],
                    linewidth=1.6,
                )
                source = [
                    row
                    for row in selected
                    if row.get("digitized_modal_loss_factor_times_100", "") != ""
                ]
                axis.plot(
                    [float(row["theta_deg"]) for row in source],
                    [float(row["digitized_modal_loss_factor_times_100"]) for row in source],
                    linestyle="none",
                    marker="o",
                    markersize=3.2,
                    markerfacecolor="white",
                    markeredgewidth=0.8,
                    color=FIGURE_1_COLORS[mode - 1],
                )
            axis.set_ylabel(r"$\eta_i\times 10^2$")
        frequency_axis.set_ylabel(r"$f$, kHz")
        for axis, label in zip(axes.ravel(), panel_labels):
            axis.text(0.02, 0.96, label, transform=axis.transAxes, ha="left", va="top")
            axis.set_xlim(0.0, 90.0)
            axis.set_xticks(np.arange(0.0, 91.0, 15.0))
            axis.set_xlabel(r"$\theta$, $^\circ$")
            axis.grid(True, color="#D9D9D9", linewidth=0.5)
        return fig


def comparison_ylim(rows_3: Sequence[Mapping[str, Any]], rows_4: Sequence[Mapping[str, Any]]) -> tuple[float, float]:
    """Return the shared Figure-3/4 y-limits from the union of plotted data."""

    values: list[float] = []
    for row in (*rows_3, *rows_4):
        if int(row["mode"]) <= PLOTTED_ROOT_COUNT:
            values.extend((float(row["timoshenko_lambda"]), float(row["eb_lambda"])))
    if not values or not np.all(np.isfinite(values)):
        raise RuntimeError("cannot determine shared y-limits from finite comparison data")
    low = min(values)
    high = max(values)
    margin = 0.04 * max(high - low, 0.05 * max(abs(high), 1.0))
    return low - margin, high + margin


def _own_ylim(rows: Sequence[Mapping[str, Any]], left_key: str, right_key: str) -> tuple[float, float]:
    values = [
        float(row[key])
        for row in rows
        if int(row["mode"]) <= PLOTTED_ROOT_COUNT
        for key in (left_key, right_key)
    ]
    low, high = min(values), max(values)
    margin = 0.04 * max(high - low, 0.05 * max(abs(high), 1.0))
    return low - margin, high + margin


def create_comparison_figure(
    rows: Sequence[Mapping[str, Any]],
    *,
    solid_key: str,
    dashed_key: str,
    ylim: tuple[float, float],
) -> plt.Figure:
    """Create one 12-line Lambda(beta) axes without a legend."""

    with plt.rc_context(_style_context()):
        fig, axis = plt.subplots(figsize=COMPARISON_FIGURE_SIZE_IN, constrained_layout=True)
        for mode in range(1, PLOTTED_ROOT_COUNT + 1):
            selected = sorted(
                [row for row in rows if int(row["mode"]) == mode],
                key=lambda row: float(row["beta_deg"]),
            )
            beta = [float(row["beta_deg"]) for row in selected]
            color = MODE_COLORS[mode - 1]
            axis.plot(
                beta,
                [float(row[solid_key]) for row in selected],
                color=color,
                linestyle="-",
                linewidth=SOLID_LINEWIDTH,
            )
            axis.plot(
                beta,
                [float(row[dashed_key]) for row in selected],
                color=color,
                linestyle=(0.0, DASH_PATTERN),
                linewidth=DASHED_LINEWIDTH,
            )
        axis.set_xlim(0.0, 90.0)
        axis.set_ylim(*ylim)
        axis.set_xticks(np.arange(0.0, 91.0, 15.0))
        axis.set_xlabel(r"$\beta$, $^\circ$")
        axis.set_ylabel(r"$\Lambda$")
        axis.grid(True, color="#D9D9D9", linewidth=0.5)
        return fig


def _assert_no_legends(fig: plt.Figure) -> None:
    if fig.legends or any(axis.get_legend() is not None for axis in fig.axes):
        raise RuntimeError("legend object detected; supervisor figures forbid legends")


def _save_figure(fig: plt.Figure, output_dir: Path, basename: str) -> list[Path]:
    _assert_no_legends(fig)
    pdf = output_dir / f"{basename}.pdf"
    png = output_dir / f"{basename}.png"
    fig.savefig(
        pdf,
        format="pdf",
        metadata={"Creator": "CoupledBeams Chapter-2 supervisor workflow", "CreationDate": None},
    )
    fig.savefig(
        png,
        format="png",
        dpi=PNG_DPI,
        metadata={"Software": "CoupledBeams Chapter-2 supervisor workflow"},
    )
    plt.close(fig)
    return [pdf, png]


def _relative_diagnostics(
    rows: Sequence[Mapping[str, Any]], difference_key: str
) -> dict[str, Any]:
    plotted = [row for row in rows if int(row["mode"]) <= PLOTTED_ROOT_COUNT]
    maximum_row = max(plotted, key=lambda row: float(row[difference_key]))
    per_mode: dict[str, dict[str, float]] = {}
    for mode in range(1, PLOTTED_ROOT_COUNT + 1):
        row = max(
            (item for item in plotted if int(item["mode"]) == mode),
            key=lambda item: float(item[difference_key]),
        )
        per_mode[str(mode)] = {
            "maximum_relative_difference": float(row[difference_key]),
            "beta_deg": float(row["beta_deg"]),
        }
    return {
        "maximum_relative_difference": float(maximum_row[difference_key]),
        "beta_deg": float(maximum_row["beta_deg"]),
        "mode": int(maximum_row["mode"]),
        "per_mode": per_mode,
    }


def _quality_summary(rows: Sequence[Mapping[str, Any]], prefixes: Sequence[str]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for prefix in prefixes:
        determinant = [float(row[f"{prefix}_accepted_determinant_residual"]) for row in rows]
        singular = [float(row[f"{prefix}_accepted_relative_singular_residual"]) for row in rows]
        summary[prefix] = {
            "accepted_root_count": len(rows),
            "maximum_accepted_determinant_residual": max(determinant),
            "maximum_accepted_relative_singular_residual": max(singular),
            "quality_basis_counts": {
                basis: sum(row[f"{prefix}_quality_basis"] == basis for row in rows)
                for basis in ("scaled", "physical_raw")
            },
        }
    return summary


def _normalization_manifest(
    preset: FigurePreset, *, fixed_reference: bool
) -> dict[str, Any]:
    return {
        "formula": "(rho*A0*omega^2*l^4/(E_x0*I_y0))^(1/4)",
        "a0_m": REFERENCE_A_M if fixed_reference else preset.a_m,
        "b0_m": REFERENCE_B_M if fixed_reference else preset.b_m,
        "l_m": 0.5 * (preset.length_1_m + preset.length_2_m),
        "E_x0": "elastic E_x at theta=0",
        "fixed_reference": bool(fixed_reference),
    }


def _oracle_family_specs() -> list[dict[str, Any]]:
    return [
        {
            "family_id": "oracle_figure02_book_slope",
            "preset": FIGURE_2_PRESET,
            "model": "book_slope_clamp",
            "csv": DEFAULT_OUTPUT_DIR / DATA_FILENAMES[2],
            "frequency_key": "book_slope_frequency_hz",
            "lambda_key": "book_slope_lambda",
            "quality_prefix": "book_slope",
        },
        {
            "family_id": "oracle_figure02_section_clamp",
            "preset": FIGURE_2_PRESET,
            "model": "timoshenko_section_clamp",
            "csv": DEFAULT_OUTPUT_DIR / DATA_FILENAMES[2],
            "frequency_key": "section_clamp_frequency_hz",
            "lambda_key": "section_clamp_lambda",
            "quality_prefix": "section_clamp",
        },
        {
            "family_id": "oracle_figure03_timoshenko_book_slope",
            "preset": FIGURES_3_4_PRESET,
            "model": "book_slope_clamp",
            "csv": DEFAULT_OUTPUT_DIR / DATA_FILENAMES[3],
            "frequency_key": "timoshenko_frequency_hz",
            "lambda_key": "timoshenko_lambda",
            "quality_prefix": "timoshenko",
        },
        {
            "family_id": "oracle_figure04_timoshenko_section_clamp",
            "preset": FIGURES_3_4_PRESET,
            "model": "timoshenko_section_clamp",
            "csv": DEFAULT_OUTPUT_DIR / DATA_FILENAMES[4],
            "frequency_key": "timoshenko_frequency_hz",
            "lambda_key": "timoshenko_lambda",
            "quality_prefix": "timoshenko",
        },
        {
            "family_id": "oracle_figure03_shared_eb",
            "preset": FIGURES_3_4_PRESET,
            "model": "rectangular_eb_saint_venant",
            "csv": DEFAULT_OUTPUT_DIR / DATA_FILENAMES[3],
            "frequency_key": "eb_frequency_hz",
            "lambda_key": "eb_lambda",
            "quality_prefix": "eb",
        },
    ]


def _legacy_recorded_boundary_and_expm_counts(
    old_manifest: Mapping[str, Any],
) -> tuple[int, int]:
    counts = old_manifest.get("matrix_evaluation_counts", {})
    boundary_count = 0
    for family in counts.values():
        if not isinstance(family, Mapping):
            continue
        for key, value in family.items():
            if key.startswith("reference_"):
                continue
            if "scaled" in key or "raw_quality" in key:
                boundary_count += int(value)
    # All five validated oracle families are two-arm assemblies.  Before the
    # exact transfer cache, each boundary assembly evaluated both arm expm
    # calls even when the arms were identical.
    return boundary_count, 2 * boundary_count


def _validate_fast_solver_against_oracle(
    output_dir: Path,
    beta_values: NDArray[np.float64],
    *,
    jobs: int,
) -> dict[str, Any]:
    if jobs != 1:
        raise ValueError(
            "sequential validation is mandatory; --jobs must equal 1"
        )
    for spec in _oracle_family_specs():
        if not Path(spec["csv"]).is_file():
            raise FileNotFoundError(f"missing oracle CSV: {spec['csv']}")
    old_manifest_path = output_dir / "plot_manifest.json"
    if not old_manifest_path.is_file():
        raise FileNotFoundError(f"missing oracle manifest: {old_manifest_path}")
    old_manifest = json.loads(old_manifest_path.read_text(encoding="utf-8"))
    aggregate = PerformanceCounters()
    transfer_cache = ExactTransferLRU(
        FastSweepSettings().cache_maxsize, counters=aggregate
    )
    validation_rows: list[dict[str, Any]] = []
    family_reports: dict[str, Any] = {}
    started = time.perf_counter()
    for spec in _oracle_family_specs():
        family_id = str(spec["family_id"])
        preset = spec["preset"]
        print(f"fast oracle validation: {family_id}", flush=True)
        result, _ = _run_or_resume_fast_family(
            preset,
            beta_values,
            str(spec["model"]),
            family_id=family_id,
            output_dir=output_dir,
            resume=False,
            transfer_cache=transfer_cache,
            counters=aggregate,
            lambda_normalization=_normalization_manifest(
                preset, fixed_reference=False
            ),
        )
        oracle_rows = _read_csv(Path(spec["csv"]))
        oracle = {
            (float(row["beta_deg"]), int(row["mode"])): row
            for row in oracle_rows
        }
        point = _point(preset, 1)
        family_frequency_max = 0.0
        family_lambda_max = 0.0
        quality_matches = True
        for beta_deg, spectrum in sorted(result.spectra.items()):
            for mode, (root, quality) in enumerate(
                zip(spectrum.roots, spectrum.quality), start=1
            ):
                key = (beta_deg, mode)
                if key not in oracle:
                    raise FastSolverValidationError(
                        f"FAST_SOLVER_VALIDATION_FAIL: oracle row missing {family_id} {key}"
                    )
                oracle_row = oracle[key]
                oracle_frequency = float(oracle_row[str(spec["frequency_key"])])
                fast_lambda = _lambda_for_root(root, preset, point)
                oracle_lambda = float(oracle_row[str(spec["lambda_key"])])
                frequency_error = abs(root.frequency_hz - oracle_frequency) / max(
                    abs(oracle_frequency), 1.0
                )
                lambda_error = abs(fast_lambda - oracle_lambda) / max(
                    abs(oracle_lambda), 1.0
                )
                same_quality = (
                    quality["quality_status"]
                    == oracle_row[f"{spec['quality_prefix']}_quality_status"]
                    == "PASS"
                )
                family_frequency_max = max(family_frequency_max, frequency_error)
                family_lambda_max = max(family_lambda_max, lambda_error)
                quality_matches = quality_matches and same_quality
                validation_rows.append(
                    {
                        "family_id": family_id,
                        "model": spec["model"],
                        "beta_deg": beta_deg,
                        "mode": mode,
                        "fast_frequency_hz": root.frequency_hz,
                        "oracle_frequency_hz": oracle_frequency,
                        "relative_frequency_error": frequency_error,
                        "fast_lambda": fast_lambda,
                        "oracle_lambda": oracle_lambda,
                        "relative_lambda_error": lambda_error,
                        "fast_quality_status": quality["quality_status"],
                        "oracle_quality_status": oracle_row[
                            f"{spec['quality_prefix']}_quality_status"
                        ],
                        "quality_acceptance_identical": same_quality,
                        "data_origin": result.data_origins[beta_deg],
                        "fallback_reason": result.fallback_reasons.get(beta_deg, ""),
                        "status": (
                            "PASS"
                            if frequency_error <= 1.0e-8
                            and lambda_error <= 1.0e-8
                            and same_quality
                            else "FAIL"
                        ),
                    }
                )
        family_reports[family_id] = {
            "maximum_relative_frequency_error": family_frequency_max,
            "maximum_relative_lambda_error": family_lambda_max,
            "quality_acceptance_identical": quality_matches,
            "anchor_maximum_relative_error": max(
                result.anchor_relative_errors.values(), default=0.0
            ),
            "fallback_count": len(result.fallback_reasons),
            "fallback_reasons": {
                str(beta): reason
                for beta, reason in sorted(result.fallback_reasons.items())
            },
            "runtime_seconds": result.runtime_seconds,
            "performance_counters": result.counters.to_dict(),
        }
    fast_runtime = time.perf_counter() - started
    _write_csv(output_dir / FAST_VALIDATION_FILENAME, validation_rows)
    maximum_frequency_error = max(
        float(row["relative_frequency_error"]) for row in validation_rows
    )
    maximum_lambda_error = max(
        float(row["relative_lambda_error"]) for row in validation_rows
    )
    all_rows_pass = all(row["status"] == "PASS" for row in validation_rows)
    speedup = LEGACY_RECORDED_RUNTIME_S / fast_runtime
    performance_pass = fast_runtime <= 180.0 or speedup >= 5.0
    correctness_pass = (
        all_rows_pass
        and len(validation_rows)
        == len(_oracle_family_specs()) * len(beta_values) * GUARD_ROOT_COUNT
        and maximum_frequency_error <= 1.0e-8
        and maximum_lambda_error <= 1.0e-8
    )
    if not correctness_pass:
        status = "FAST_SOLVER_FAIL"
    elif not performance_pass:
        status = "FAST_SOLVER_CORRECT_BUT_PERFORMANCE_TARGET_NOT_MET"
    else:
        status = "FAST_SOLVER_PASS"
    legacy_boundary, legacy_expm = _legacy_recorded_boundary_and_expm_counts(
        old_manifest
    )
    benchmark = {
        "status": status,
        "solver_version": FAST_SOLVER_VERSION,
        "jobs": jobs,
        "beta_point_count": len(beta_values),
        "root_count": GUARD_ROOT_COUNT,
        "validated_family_count": len(_oracle_family_specs()),
        "validated_row_count": len(validation_rows),
        "oracle_frequency_relative_tolerance": 1.0e-8,
        "oracle_lambda_relative_tolerance": 1.0e-8,
        "maximum_relative_frequency_error": maximum_frequency_error,
        "maximum_relative_lambda_error": maximum_lambda_error,
        "legacy_recorded_runtime_s": LEGACY_RECORDED_RUNTIME_S,
        "fast_sequential_runtime_s": fast_runtime,
        "speedup": speedup,
        "performance_target": "runtime <= 180 s or speedup >= 5",
        "performance_target_pass": performance_pass,
        "legacy_recorded_boundary_matrix_evaluations": legacy_boundary,
        "legacy_transfer_expm_evaluations_for_five_oracle_families": legacy_expm,
        "legacy_expm_count_basis": "two uncached equal-arm transfers per coupled boundary assembly",
        "fast_performance_counters": aggregate.to_dict(),
        "global_anchor_count": aggregate.global_anchor_scans,
        "global_fallback_count": aggregate.global_fallback_scans,
        "fallback_rate": aggregate.global_fallback_scans
        / max(len(_oracle_family_specs()) * len(beta_values), 1),
        "family_reports": family_reports,
        "oracle_csv_paths": sorted(
            {str(Path(spec["csv"]).relative_to(ROOT)) for spec in _oracle_family_specs()}
        ),
        "validation_csv": str(
            (output_dir / FAST_VALIDATION_FILENAME).relative_to(ROOT)
        ),
    }
    (output_dir / FAST_BENCHMARK_FILENAME).write_text(
        json.dumps(benchmark, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return benchmark


def _require_fast_solver_pass(output_dir: Path) -> dict[str, Any]:
    path = output_dir / FAST_BENCHMARK_FILENAME
    if not path.is_file():
        raise RuntimeError(
            "Figures 5-8 require completed --validate-fast-solver evidence"
        )
    benchmark = json.loads(path.read_text(encoding="utf-8"))
    if benchmark.get("status") != "FAST_SOLVER_PASS":
        raise RuntimeError(
            "Figures 5-8 are blocked because fast solver status is not PASS"
        )
    if float(benchmark.get("maximum_relative_frequency_error", math.inf)) > 1.0e-8:
        raise RuntimeError("saved fast solver oracle error exceeds 1e-8")
    return benchmark


def _compute_extended_figure_data(
    output_dir: Path,
    beta_values: NDArray[np.float64],
    figures: set[int],
    *,
    resume: bool,
) -> tuple[dict[int, list[dict[str, Any]]], dict[str, Any], float]:
    _require_fast_solver_pass(output_dir)
    figure_3_path = output_dir / DATA_FILENAMES[3]
    if not figure_3_path.is_file():
        raise FileNotFoundError(
            "Figures 7-8 require the validated Figure-3 CSV for exact reuse"
        )
    rows_3 = _comparison_rows_as_numbers(_read_csv(figure_3_path), 3)
    aggregate = PerformanceCounters()
    transfer_cache = ExactTransferLRU(
        FastSweepSettings().cache_maxsize, counters=aggregate
    )
    rows_by_figure: dict[int, list[dict[str, Any]]] = {}
    family_metadata: dict[str, Any] = {}
    current_scientific_runtime = 0.0
    normalization = _normalization_manifest(
        FIGURES_3_4_PRESET, fixed_reference=True
    )

    def run_family(
        family_id: str, preset: FigurePreset, model: str
    ) -> FastSweepResult[SpectrumResult]:
        nonlocal current_scientific_runtime
        result, reused = _run_or_resume_fast_family(
            preset,
            beta_values,
            model,
            family_id=family_id,
            output_dir=output_dir,
            resume=resume,
            transfer_cache=transfer_cache,
            counters=aggregate,
            lambda_normalization=normalization,
        )
        if not reused:
            current_scientific_runtime += result.runtime_seconds
        family_metadata[family_id] = {
            "model": model,
            "preset": _preset_manifest(preset),
            "resumed_from_checkpoint": reused,
            "runtime_seconds": result.runtime_seconds,
            "performance_counters": result.counters.to_dict(),
            "global_fallback_count": len(result.fallback_reasons),
            "fallback_reasons": {
                str(beta): reason
                for beta, reason in sorted(result.fallback_reasons.items())
            },
            "maximum_anchor_relative_error": max(
                result.anchor_relative_errors.values(), default=0.0
            ),
        }
        return result

    if 5 in figures:
        left = run_family(
            "figure05_timoshenko", FIGURE_5_PRESET, "book_slope_clamp"
        )
        right = run_family(
            "figure05_eb", FIGURE_5_PRESET, "rectangular_eb_saint_venant"
        )
        rows_by_figure[5] = _extended_figure_rows(
            5,
            FIGURE_5_PRESET,
            beta_values,
            left,
            comparison_type="timoshenko_vs_eb_unequal_lengths",
            left_model="Chapter2_Timoshenko_state_corrected_generalized_torsion",
            left_theta_deg=0.0,
            right_model="rectangular_orthotropic_EB_Saint_Venant",
            right_theta_deg=0.0,
            right_fast=right,
        )
    if 6 in figures:
        left = run_family(
            "figure06_timoshenko", FIGURE_6_PRESET, "book_slope_clamp"
        )
        right = run_family(
            "figure06_eb", FIGURE_6_PRESET, "rectangular_eb_saint_venant"
        )
        rows_by_figure[6] = _extended_figure_rows(
            6,
            FIGURE_6_PRESET,
            beta_values,
            left,
            comparison_type="timoshenko_vs_eb_unequal_lengths_and_thickness",
            left_model="Chapter2_Timoshenko_state_corrected_generalized_torsion",
            left_theta_deg=0.0,
            right_model="rectangular_orthotropic_EB_Saint_Venant",
            right_theta_deg=0.0,
            right_fast=right,
        )
    if 7 in figures:
        left = run_family(
            "figure07_timoshenko_theta5", FIGURE_7_PRESET, "book_slope_clamp"
        )
        rows_by_figure[7] = _extended_figure_rows(
            7,
            FIGURE_7_PRESET,
            beta_values,
            left,
            comparison_type="diagnostic_orthotropic_theta0_approximation_for_weak_monoclinic_theta5",
            left_model="Chapter2_monoclinic_Timoshenko_theta5",
            left_theta_deg=5.0,
            right_model="rectangular_orthotropic_EB_theta0_Saint_Venant_approximation",
            right_theta_deg=0.0,
            reused_figure_3_rows=rows_3,
            reused_prefix="eb",
        )
    if 8 in figures:
        left = run_family(
            "figure08_timoshenko_theta15", FIGURE_8_PRESET, "book_slope_clamp"
        )
        rows_by_figure[8] = _extended_figure_rows(
            8,
            FIGURE_8_PRESET,
            beta_values,
            left,
            comparison_type="chapter2_rotated_material_axes_and_Sbar16_coupling",
            left_model="Chapter2_theta15",
            left_theta_deg=15.0,
            right_model="Chapter2_theta0",
            right_theta_deg=0.0,
            reused_figure_3_rows=rows_3,
            reused_prefix="timoshenko",
        )

    source_eb = np.asarray(
        [float(row["eb_frequency_hz"]) for row in rows_3], dtype=float
    )
    source_timo = np.asarray(
        [float(row["timoshenko_frequency_hz"]) for row in rows_3], dtype=float
    )
    reuse_checks = {
        "figure_7_eb_equals_figure_3_eb": (
            bool(
                np.array_equal(
                    np.asarray(
                        [row["right_frequency_hz"] for row in rows_by_figure[7]],
                        dtype=float,
                    ),
                    source_eb,
                )
            )
            if 7 in rows_by_figure
            else None
        ),
        "figure_8_theta0_timoshenko_equals_figure_3_timoshenko": (
            bool(
                np.array_equal(
                    np.asarray(
                        [row["right_frequency_hz"] for row in rows_by_figure[8]],
                        dtype=float,
                    ),
                    source_timo,
                )
            )
            if 8 in rows_by_figure
            else None
        ),
    }
    if any(value is False for value in reuse_checks.values()):
        raise RuntimeError("Figure 7/8 reused baseline arrays are not exactly identical")
    metadata = {
        "solver_version": FAST_SOLVER_VERSION,
        "fixed_reference_normalization": normalization,
        "families": family_metadata,
        "reuse_checks": reuse_checks,
        "current_run_performance_counters": aggregate.to_dict(),
    }
    return rows_by_figure, metadata, current_scientific_runtime


def eb_arrays_are_identical(
    rows_3: Sequence[Mapping[str, Any]], rows_4: Sequence[Mapping[str, Any]]
) -> bool:
    left = np.asarray(
        [
            (float(row["beta_deg"]), int(row["mode"]), float(row["eb_frequency_hz"]), float(row["eb_lambda"]))
            for row in rows_3
        ],
        dtype=float,
    )
    right = np.asarray(
        [
            (float(row["beta_deg"]), int(row["mode"]), float(row["eb_frequency_hz"]), float(row["eb_lambda"]))
            for row in rows_4
        ],
        dtype=float,
    )
    return bool(np.array_equal(left, right))


def _safe_output_dir(path: Path) -> Path:
    resolved = path.resolve()
    for workspace in ARTICLE_WORKSPACES:
        article = workspace.resolve()
        try:
            resolved.relative_to(article)
        except ValueError:
            continue
        raise ValueError(f"output directory must not be inside article workspace: {resolved}")
    return resolved


def _git(*args: str) -> str:
    completed = subprocess.run(
        ["git", *args],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return completed.stdout.strip()


def _git_context() -> dict[str, Any]:
    return {
        "cwd": str(ROOT),
        "git_root": _git("rev-parse", "--show-toplevel"),
        "branch": _git("branch", "--show-current"),
        "head": _git("rev-parse", "HEAD"),
        "origin_main": _git("rev-parse", "origin/main"),
        "status_short": _git("status", "--short").splitlines(),
    }


def _reference_manifest(reference: SectionReferenceResult) -> dict[str, Any]:
    return {
        "material_and_geometry": reference.preset_label,
        "status": reference.status,
        "maximum_relative_frequency_difference": reference.maximum_relative_frequency_difference,
        "criterion": 1.0e-8,
        "compared_root_count": GUARD_ROOT_COUNT,
        "coupled_frequencies_hz": [root.frequency_hz for root in reference.coupled.roots],
        "straight_frequencies_hz": [root.frequency_hz for root in reference.straight.roots],
        "scaled_matrix_evaluations": (
            reference.coupled.scaled_matrix_evaluations
            + reference.straight.scaled_matrix_evaluations
        ),
        "raw_quality_matrix_evaluations": (
            reference.coupled.raw_quality_matrix_evaluations
            + reference.straight.raw_quality_matrix_evaluations
        ),
        "runtime_seconds": reference.coupled.runtime_seconds + reference.straight.runtime_seconds,
        "independent_reference_uses_J_book": False,
    }


def _preset_manifest(preset: FigurePreset) -> dict[str, Any]:
    material = preset.material_factory()
    point = _point(preset, 1)
    arm_1_area, arm_1_inertia_y = rectangular_reference_section(
        preset.a_m, preset.b_m
    )
    a_2_m = preset.a_m if preset.a_2_m is None else preset.a_2_m
    b_2_m = preset.b_m if preset.b_2_m is None else preset.b_2_m
    arm_2_area, arm_2_inertia_y = rectangular_reference_section(a_2_m, b_2_m)
    return {
        "material_name": preset.material_name,
        "material_factory": preset.material_factory.__name__,
        "material_mode": preset.material_mode,
        "rho_kg_m3": material.rho,
        "elastic_Ex_theta0_pa": float(np.real(point.properties.Ex)),
        "theta_1_deg": preset.theta_1_deg,
        "theta_2_deg": preset.theta_2_deg,
        "a_1_m": preset.a_m,
        "a_2_m": a_2_m,
        "b_1_m": preset.b_m,
        "b_2_m": b_2_m,
        "L_1_m": preset.length_1_m,
        "L_2_m": preset.length_2_m,
        "mu": preset.mu,
        "arm_1_area_m2": arm_1_area,
        "arm_1_I_y_m4": arm_1_inertia_y,
        "arm_2_area_m2": arm_2_area,
        "arm_2_I_y_m4": arm_2_inertia_y,
    }


def _figure_1_source_manifest() -> dict[str, Any]:
    return {
        "workflow": "scripts/analysis/anisotropic_rods/reproduce_yartsev_fig_2_2.py",
        "status": "PASS_WITHIN_GRAPH_RESOLUTION",
        "material": "BookMaterial() / T-53(VM)-78/PN-609-21M",
        "material_mode": "book_complex",
        "formulation": "state_corrected",
        "roots_csv": str(CANONICAL_FIGURE_1_ROOTS.relative_to(ROOT)),
        "roots_csv_sha256": _sha256(CANONICAL_FIGURE_1_ROOTS),
        "digitized_csv": str(CANONICAL_FIGURE_1_DIGITIZED.relative_to(ROOT)),
        "digitized_csv_sha256": _sha256(CANONICAL_FIGURE_1_DIGITIZED),
        "comparison_csv": str(CANONICAL_FIGURE_1_COMPARISON.relative_to(ROOT)),
        "comparison_csv_sha256": _sha256(CANONICAL_FIGURE_1_COMPARISON),
        "new_material_fitting": False,
    }


def _report_text(manifest: Mapping[str, Any]) -> str:
    diagnostics = manifest["diagnostics"]
    fig2 = diagnostics.get("figure_2", {})
    fig3 = diagnostics.get("figure_3", {})
    fig4 = diagnostics.get("figure_4", {})
    references = manifest.get("section_clamp_reference", {})
    benchmark = manifest.get("fast_solver_benchmark", {})
    extended = manifest.get("extended_fast_solver", {})
    quality = manifest.get("root_quality", {})
    lines = [
        "# Рисунки для отчёта научному руководителю: глава 2 монографии Ярцева",
        "",
        f"**Supervisor figure workflow: {manifest['workflow_status']}**",
        "",
        f"**Fast beta-sweep optimization: {manifest.get('fast_solver_status', 'FAIL')}**",
        "",
        f"**Extended supervisor figures 5–8: {manifest.get('extended_figures_status', 'FAIL')}**",
        "",
        "Этот отчёт относится только к анизотропным/ортотропным прямоугольным стержням главы 2. Параллельная статья о круглых изотропных стержнях, её scripts и результаты не использовались и не изменялись.",
        "",
        "Для Figures 2–8 показаны первые шесть sorted spectral positions при каждом фиксированном угле. The plotted curves are sorted spectral positions 1–6 at every beta, not tracked modal descendants. Седьмой положительный корень является completeness guard.",
        "",
        "## Fast solver и oracle-validation",
        "",
        "Быстрый diagnostic-only solver использует предыдущую точку или линейный predictor по двум предыдущим точкам только для построения локальных окон. Близкие корни с относительным gap менее `1e-3` объединяются в connected cluster. Полные global inventories обязательны при beta `0, 15, 30, 45, 60, 75, 90` градусов и автоматически используются при любом отказе local path. После fallback predictor продолжается от принятого global spectrum. Legacy global solver сохранён без изменения как oracle и fallback path.",
        "",
        "Matrix exponential кэшируется bounded LRU только по точному IEEE-754 значению omega, типу модели и полному immutable fingerprint материала/геометрии. Округления omega и приближённого determinant cache нет.",
        "",
    ]
    if benchmark:
        counters = benchmark.get("fast_performance_counters", {})
        lines += [
            f"Полная проверка охватила `{benchmark.get('validated_row_count', 0)}` oracle roots. Максимальная относительная ошибка частоты: `{float(benchmark.get('maximum_relative_frequency_error', math.nan)):.12e}`; максимальная относительная ошибка Lambda: `{float(benchmark.get('maximum_relative_lambda_error', math.nan)):.12e}`.",
            f"Recorded legacy runtime: `{float(benchmark.get('legacy_recorded_runtime_s', 0.0)):.6f} s`; fast sequential runtime: `{float(benchmark.get('fast_sequential_runtime_s', 0.0)):.6f} s`; speedup: `{float(benchmark.get('speedup', 0.0)):.6f}`.",
            f"Global anchors: `{benchmark.get('global_anchor_count', 0)}`; global fallbacks: `{benchmark.get('global_fallback_count', 0)}`; fallback rate: `{float(benchmark.get('fallback_rate', 0.0)):.6f}`.",
            f"Legacy estimated transfer expm count: `{benchmark.get('legacy_transfer_expm_evaluations_for_five_oracle_families', 0)}`; fast expm count: `{counters.get('transfer_expm_evaluations', 0)}`; cache hit rate: `{float(counters.get('transfer_cache_hit_rate', 0.0)):.6f}`.",
            "",
        ]
    lines += [
        "## Нормировка",
        "",
        "Безразмерная частота для Figures 2–4:",
        "",
        "```text",
        "Lambda = (rho*A*omega^2*l^4/(E_x*I_y))^(1/4)",
        "l = (L1+L2)/2, A=a*b, I_y=a^3*b/12, omega=2*pi*f",
        "```",
        "",
        "Для Figures 5–8 используется единый fixed rectangular reference: `a0=5 mm`, `b0=20 mm`, `A0=a0*b0`, `I_y0=a0^3*b0/12`, `E_x0=E_x(theta=0)` и `l=(L1+L2)/2`. Эта нормировка одинакова для обеих кривых каждого рисунка и не зависит от фактической толщины плеча Figure 6 или угла theta Figures 7–8.",
        "",
        "Старый thickness-mismatch parameter с именем eta в presets отсутствует; обозначение модального коэффициента потерь на Figure 1 относится только к книжному Figure 2.2.",
        "",
        "## Готовые подписи",
        "",
        "**Рисунок 1.** Воспроизведение расчётных зависимостей рисунка 2.2 монографии для исходного материала T-53(VM)-78/PN-609-21M и книжной геометрии. Линии — расчёт по реализованной исправленной постановке `state_corrected` с книжными комплексными свойствами; маркеры — округлённые значения рассчитанных кривых, считанные с исходного рисунка. Экспериментальные маркеры монографии не были оцифрованы и здесь не показаны.",
        "",
        "**Рисунок 2.** Влияние постановки внешней заделки в новой Chapter-2 модели для книжного материала и геометрии. Сплошные линии — `book_slope_clamp`; пунктирные линии — `timoshenko_section_clamp`. Показаны первые шесть sorted frequencies. Кривые соответствуют первым шести sorted spectral positions, а не отслеженным модальным ветвям.",
        "",
        "**Рисунок 3.** Сравнение моделей для прямоугольных стержней HMS/DX-209 при книжной заделке. Сплошные линии — модель моноклинного прямоугольного стержня Тимошенко с обобщённым кручением. Пунктирные линии — прямоугольная ортотропная модель Эйлера–Бернулли с кручением Сен-Венана. На внешних концах используется `book_slope_clamp`. Кривые соответствуют первым шести sorted spectral positions, а не отслеженным модальным ветвям.",
        "",
        "**Рисунок 4.** Сравнение моделей для прямоугольных стержней HMS/DX-209 при полной заделке сечения. Сплошные линии — модель моноклинного прямоугольного стержня Тимошенко с обобщённым кручением. Пунктирные линии — прямоугольная ортотропная модель Эйлера–Бернулли с кручением Сен-Венана. На внешних концах используется `timoshenko_section_clamp`. Кривые соответствуют первым шести sorted spectral positions, а не отслеженным модальным ветвям.",
        "",
        "**Рисунок 5.** Сравнение моделей для прямоугольных стержней HMS/DX-209 с различными длинами плеч: `L1=0.3` м, `L2=0.5` м, `a1=a2=5` мм. Сплошные линии — модель моноклинного прямоугольного стержня Тимошенко с обобщённым кручением; пунктирные линии — прямоугольная ортотропная модель Эйлера–Бернулли с кручением Сен-Венана. Используется `book_slope_clamp`. Показаны первые шесть sorted spectral positions при каждом фиксированном beta; это не отслеженные модальные ветви.",
        "",
        "**Рисунок 6.** Сравнение моделей для прямоугольных стержней HMS/DX-209 с различными длинами и толщинами плеч: `L1=0.3` м, `L2=0.5` м, `a1=4` мм, `a2=6` мм. Геометрия задана непосредственно без массо-сохраняющей параметризации. Сплошные линии — модель Тимошенко с обобщённым кручением; пунктирные линии — модель Эйлера–Бернулли с кручением Сен-Венана. Показаны первые шесть sorted spectral positions при каждом фиксированном beta; это не отслеженные модальные ветви.",
        "",
        "**Рисунок 7.** Сопоставление полной моноклинной модели главы 2 при `theta1=theta2=5°` с ортотропным приближением Эйлера–Бернулли при `theta=0°`. Сплошные линии учитывают поворот эффективных свойств и материальную изгибно-крутильную связь; пунктирные линии соответствуют приближению, в котором эта связь полностью отсутствует. Рисунок является диагностикой применимости ортотропного приближения, а не чистым сравнением только сдвиговой и вращательно-инерционной поправок. Показаны первые шесть sorted spectral positions при каждом фиксированном beta; это не отслеженные модальные ветви.",
        "",
        "**Рисунок 8.** Влияние поворота материальных осей в полной модели главы 2. Сплошные линии — `theta1=theta2=15°`; пунктирные линии — `theta1=theta2=0°`. Обе группы кривых рассчитаны по одной модели Тимошенко с обобщённым кручением и `book_slope_clamp`. Показаны первые шесть sorted spectral positions при каждом фиксированном beta; это не отслеженные модальные ветви.",
        "",
        "Во всех восьми рисунках legends отсутствуют. Один цвет соответствует одной sorted spectral position; Figures 2–8 используют один детерминированный шестицветный цикл.",
        "",
        "## Численные diagnostics",
        "",
    ]
    if fig2:
        lines += [
            f"- Figure 2: maximum relative Lambda clamp difference `{fig2['maximum_relative_difference']:.12e}` at beta=`{fig2['beta_deg']:g}` deg, mode `{fig2['mode']}`.",
            f"- Figure 2 section-clamp straight-reference residual: `{references['figure_2_book_material']['maximum_relative_frequency_difference']:.12e}`.",
            f"- Figure 2 per-mode maxima: `{json.dumps(fig2['per_mode'], ensure_ascii=False, sort_keys=True)}`.",
        ]
    if fig3:
        lines += [
            f"- Figure 3: maximum relative Lambda Timoshenko–EB difference `{fig3['maximum_relative_difference']:.12e}` at beta=`{fig3['beta_deg']:g}` deg, mode `{fig3['mode']}`.",
            f"- Figure 3 per-mode maxima: `{json.dumps(fig3['per_mode'], ensure_ascii=False, sort_keys=True)}`.",
        ]
    if fig4:
        lines += [
            f"- Figure 4: maximum relative Lambda Timoshenko–EB difference `{fig4['maximum_relative_difference']:.12e}` at beta=`{fig4['beta_deg']:g}` deg, mode `{fig4['mode']}`.",
            f"- Figure 4 section-clamp straight-reference residual: `{references['figures_3_4_hms_dx_209']['maximum_relative_frequency_difference']:.12e}`.",
            f"- Figure 4 per-mode maxima: `{json.dumps(fig4['per_mode'], ensure_ascii=False, sort_keys=True)}`.",
            f"- EB arrays used in Figures 3 and 4 are exactly identical: `{manifest['eb_arrays_figures_3_4_exactly_identical']}`.",
        ]
    for figure in range(5, 9):
        item = diagnostics.get(f"figure_{figure}")
        if not item:
            continue
        gap = item["minimum_gap"]
        counters = item.get("performance_counters", {})
        q = quality.get(f"figure_{figure}", {})
        lines += [
            f"- Figure {figure}: maximum relative Lambda difference `{item['maximum_relative_difference']:.12e}` at beta=`{item['beta_deg']:g}` deg, mode `{item['mode']}`; per-mode maxima: `{json.dumps(item['per_mode'], ensure_ascii=False, sort_keys=True)}`.",
            f"- Figure {figure}: minimum relative neighbor gap `{gap['minimum_relative_neighbor_gap']:.12e}` at beta=`{gap['beta_deg']:g}` deg, pair `{gap['pair']}`, side `{gap['model_side']}`.",
            f"- Figure {figure}: global fallbacks `{item.get('global_fallback_count', 0)}`; scientific runtime `{float(item.get('scientific_runtime_s', 0.0)):.6f} s`; transfer expm `{counters.get('transfer_expm_evaluations', 0)}`; cache hit rate `{float(counters.get('transfer_cache_hit_rate', 0.0)):.6f}`.",
            f"- Figure {figure}: root-quality maxima `{json.dumps(q, ensure_ascii=False, sort_keys=True)}`.",
        ]
        if item.get("fallback_reasons"):
            lines.append(
                f"- Figure {figure}: fallback reasons `{json.dumps(item['fallback_reasons'], ensure_ascii=False, sort_keys=True)}`."
            )
        if item.get("interpretation"):
            lines.append(f"- Figure {figure}: {item['interpretation']}.")
    reuse = extended.get("reuse_checks", {})
    if reuse:
        lines += [
            f"- Figure 7 reused EB array is exactly identical to Figure-3 EB: `{reuse.get('figure_7_eb_equals_figure_3_eb')}`.",
            f"- Figure 8 reused theta=0 Timoshenko array is exactly identical to Figure-3 Timoshenko: `{reuse.get('figure_8_theta0_timoshenko_equals_figure_3_timoshenko')}`.",
        ]
    lines += [
        "",
        "Различия между теориями и заделками являются diagnostics без acceptance threshold и не ранжируют физическую точность моделей.",
        "",
        "## Параметры и ограничения",
        "",
        f"- beta grid: `{manifest['beta_grid']['point_count']}` points, 0…90 deg, step `{manifest['beta_grid']['requested_step_deg']}` deg.",
        "- Figure 2: `BookMaterial()`, `a=9.76 mm`, `b=25.24 mm`, `L1=L2=0.295 m`, `theta1=theta2=0`, elastic, `mu=0`.",
        "- Figures 3–4: `hms_dx_209_material()`, `a=5 mm`, `b=20 mm`, `L1=L2=0.4 m`, `theta1=theta2=0`, elastic, `mu=0`.",
        "- Figure 5: `hms_dx_209_material()`, `L1=0.3 m`, `L2=0.5 m`, `a1=a2=5 mm`, `b1=b2=20 mm`, `theta1=theta2=0`, elastic, `mu=0.25`, `book_slope_clamp`.",
        "- Figure 6: `hms_dx_209_material()`, `L1=0.3 m`, `L2=0.5 m`, `a1=4 mm`, `a2=6 mm`, `b1=b2=20 mm`, `theta1=theta2=0`, elastic, `mu=0.25`, direct geometry, `book_slope_clamp`.",
        "- Figure 7: full Chapter-2 Timoshenko at `theta1=theta2=5°` versus the exactly reused Figure-3 orthotropic EB `theta=0°` approximation; all geometric dimensions are `a=5 mm`, `b=20 mm`, `L1=L2=0.4 m`.",
        "- Figure 8: Chapter-2 Timoshenko at `theta1=theta2=15°` versus the exactly reused Figure-3 Chapter-2 Timoshenko array at `theta1=theta2=0°`; the same `5×20 mm`, `L1=L2=0.4 m` geometry is used.",
        "- Figure 1 uses the canonical verified Figure-2.2 source-reproduction CSV and remains in the book coordinates `(theta, f, modal loss factor)`; it is not converted to Lambda(beta).",
        "- Root quality passes when the scaled branch or normalized physical-raw branch has determinant and relative singular residuals no greater than `1e-8`, and the root status is not rejected.",
        "- No interpolation is used for missing or rejected roots.",
        "- No MAC, shape tracking, nearest-frequency continuation, mode shapes, FEM, or 3D calculation is part of this workflow.",
        "- Outputs are written only below `results/anisotropic_rods/yartsev_ch2_supervisor_figures/`, outside every article workspace.",
        "",
        "## Воспроизводимость",
        "",
        "```text",
        "python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py",
        "python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py --reuse-data",
        "python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py --validate-fast-solver --jobs 1",
        "python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py --benchmark-fast-solver --jobs 1",
        "python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py --figure new --solver-mode fast --resume --jobs 1",
        "python scripts/analysis/anisotropic_rods/plot_yartsev_ch2_supervisor_figures.py --figure all --reuse-data",
        "```",
        "",
        f"Scientific/root runtime recorded by the producing run: `{manifest['runtimes_seconds'].get('scientific_total', 0.0):.6f} s`; plotting/report runtime for this run: `{manifest['runtimes_seconds'].get('plotting_and_reporting', 0.0):.6f} s`.",
    ]
    return "\n".join(lines) + "\n"


def _write_failure(output_dir: Path, status: str, details: str) -> None:
    _write_csv(
        output_dir / "failure_diagnostics.csv",
        [{"workflow_status": status, "details": details}],
    )


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = _safe_output_dir(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.jobs < 1:
        raise ValueError("--jobs must be positive")
    if args.jobs != 1:
        raise ValueError(
            "family-level multiprocessing is not enabled because the sequential "
            "path must be validated and benchmarked first"
        )
    if args.figure == "all":
        selected = set(range(1, 9))
    elif args.figure == "new":
        selected = {5, 6, 7, 8}
    else:
        selected = {int(args.figure)}
    beta_values = beta_grid(args.beta_step_deg)
    old_manifest_path = output_dir / "plot_manifest.json"
    old_manifest = (
        json.loads(old_manifest_path.read_text(encoding="utf-8"))
        if old_manifest_path.is_file()
        else {}
    )
    compute_metadata: dict[str, Any] = {}
    references: dict[str, Any] = dict(old_manifest.get("section_clamp_reference", {}))
    scientific_runtime = 0.0
    scientific_counts: dict[str, Any] = dict(old_manifest.get("matrix_evaluation_counts", {}))

    if args.validate_fast_solver:
        try:
            benchmark = _validate_fast_solver_against_oracle(
                output_dir, beta_values, jobs=args.jobs
            )
        except Exception as error:
            _write_failure(output_dir, "FAST_SOLVER_VALIDATION_FAIL", str(error))
            print("Fast beta-sweep optimization: FAIL", file=sys.stderr)
            print(str(error), file=sys.stderr)
            return 1
        print(json.dumps(benchmark, ensure_ascii=False, indent=2, sort_keys=True))
        if benchmark["status"] == "FAST_SOLVER_PASS":
            stale_failure_path = output_dir / "failure_diagnostics.csv"
            if stale_failure_path.is_file():
                stale_failure_path.unlink()
            print("Fast beta-sweep optimization: PASS")
            return 0
        if benchmark["status"] == "FAST_SOLVER_CORRECT_BUT_PERFORMANCE_TARGET_NOT_MET":
            print(
                "Fast beta-sweep optimization: "
                "CORRECT_BUT_PERFORMANCE_TARGET_NOT_MET"
            )
            return 1
        print("Fast beta-sweep optimization: FAIL")
        return 1

    if args.benchmark_fast_solver:
        benchmark_path = output_dir / FAST_BENCHMARK_FILENAME
        if not benchmark_path.is_file():
            print(
                "Fast beta-sweep optimization: FAIL\n"
                "run --validate-fast-solver first",
                file=sys.stderr,
            )
            return 1
        benchmark = json.loads(benchmark_path.read_text(encoding="utf-8"))
        print(json.dumps(benchmark, ensure_ascii=False, indent=2, sort_keys=True))
        if benchmark.get("status") != "FAST_SOLVER_PASS":
            print("Fast beta-sweep optimization: FAIL", file=sys.stderr)
            return 1
        print("Fast beta-sweep optimization: PASS")
        return 0

    try:
        figure_1_path = output_dir / DATA_FILENAMES[1]
        if 1 in selected:
            if args.force_recompute or not figure_1_path.is_file():
                rows_1 = _build_figure_1_rows(reuse_data=args.reuse_data)
                _write_csv(figure_1_path, rows_1)
            else:
                rows_1 = _read_csv(figure_1_path)
        else:
            rows_1 = _read_csv(figure_1_path) if figure_1_path.is_file() else []

        figure_2_path = output_dir / DATA_FILENAMES[2]
        need_compute_2 = 2 in selected and (args.force_recompute or not figure_2_path.is_file())
        if need_compute_2:
            if args.reuse_data:
                raise FileNotFoundError(f"--reuse-data requires {figure_2_path}")
            if args.solver_mode == "fast":
                _require_fast_solver_pass(output_dir)
            reference_2 = section_clamp_straight_reference(FIGURE_2_PRESET)
            if reference_2.status != "PASS":
                raise RuntimeError(
                    "FAIL_SECTION_CLAMP_REFERENCE: Figure 2 residual "
                    f"{reference_2.maximum_relative_frequency_difference:.6e}"
                )
            references["figure_2_book_material"] = _reference_manifest(reference_2)
            counters_2 = PerformanceCounters()
            cache_2 = ExactTransferLRU(
                FastSweepSettings().cache_maxsize, counters=counters_2
            )
            book_2 = _selected_sweep(
                FIGURE_2_PRESET,
                beta_values,
                "book_slope_clamp",
                solver_mode=args.solver_mode,
                transfer_cache=cache_2,
                counters=counters_2,
            )
            section_2 = _selected_sweep(
                FIGURE_2_PRESET,
                beta_values,
                "timoshenko_section_clamp",
                solver_mode=args.solver_mode,
                beta_zero_result=reference_2.coupled,
                transfer_cache=cache_2,
                counters=counters_2,
            )
            rows_2 = _figure_2_rows(beta_values, book_2, section_2)
            _write_csv(figure_2_path, rows_2)
            scientific_runtime += (
                book_2.runtime_seconds
                + section_2.runtime_seconds
                + reference_2.coupled.runtime_seconds
                + reference_2.straight.runtime_seconds
            )
            scientific_counts["figure_2"] = {
                "book_slope_scaled": book_2.scaled_matrix_evaluations,
                "book_slope_raw_quality": book_2.raw_quality_matrix_evaluations,
                "section_clamp_scaled_including_reused_beta0": section_2.scaled_matrix_evaluations,
                "section_clamp_raw_quality_including_reused_beta0": section_2.raw_quality_matrix_evaluations,
                "book_slope_completeness_refined_beta_count": book_2.completeness_refined_beta_count,
                "book_slope_completeness_refinement_attempts": book_2.completeness_refinement_attempts,
                "section_clamp_completeness_refined_beta_count": section_2.completeness_refined_beta_count,
                "section_clamp_completeness_refinement_attempts": section_2.completeness_refinement_attempts,
                "reference_scaled": references["figure_2_book_material"]["scaled_matrix_evaluations"],
                "reference_raw_quality": references["figure_2_book_material"]["raw_quality_matrix_evaluations"],
                "solver_mode": args.solver_mode,
                "fast_performance_counters": (
                    counters_2.to_dict() if args.solver_mode == "fast" else {}
                ),
            }
        elif figure_2_path.is_file():
            rows_2 = _comparison_rows_as_numbers(_read_csv(figure_2_path), 2)
        else:
            rows_2 = []

        figure_3_path = output_dir / DATA_FILENAMES[3]
        figure_4_path = output_dir / DATA_FILENAMES[4]
        comparison_selected = bool(selected & {3, 4})
        need_compute_34 = comparison_selected and (
            args.force_recompute
            or not figure_3_path.is_file()
            or not figure_4_path.is_file()
        )
        if need_compute_34:
            if args.reuse_data:
                raise FileNotFoundError(
                    f"--reuse-data requires both {figure_3_path} and {figure_4_path}"
                )
            if args.solver_mode == "fast":
                _require_fast_solver_pass(output_dir)
            reference_34 = section_clamp_straight_reference(FIGURES_3_4_PRESET)
            if reference_34.status != "PASS":
                raise RuntimeError(
                    "FAIL_SECTION_CLAMP_REFERENCE: Figures 3-4 residual "
                    f"{reference_34.maximum_relative_frequency_difference:.6e}"
                )
            references["figures_3_4_hms_dx_209"] = _reference_manifest(reference_34)
            counters_34 = PerformanceCounters()
            cache_34 = ExactTransferLRU(
                FastSweepSettings().cache_maxsize, counters=counters_34
            )
            book_34 = _selected_sweep(
                FIGURES_3_4_PRESET,
                beta_values,
                "book_slope_clamp",
                solver_mode=args.solver_mode,
                transfer_cache=cache_34,
                counters=counters_34,
            )
            section_34 = _selected_sweep(
                FIGURES_3_4_PRESET,
                beta_values,
                "timoshenko_section_clamp",
                solver_mode=args.solver_mode,
                beta_zero_result=reference_34.coupled,
                transfer_cache=cache_34,
                counters=counters_34,
            )
            eb_34 = _selected_sweep(
                FIGURES_3_4_PRESET,
                beta_values,
                "rectangular_eb_saint_venant",
                solver_mode=args.solver_mode,
                transfer_cache=cache_34,
                counters=counters_34,
            )
            rows_3 = _figure_3_4_rows(
                3, beta_values, book_34, eb_34, "book_slope_clamp"
            )
            rows_4 = _figure_3_4_rows(
                4,
                beta_values,
                section_34,
                eb_34,
                "timoshenko_section_clamp",
            )
            if not eb_arrays_are_identical(rows_3, rows_4):
                raise RuntimeError("EB arrays for Figures 3 and 4 are not exactly identical")
            _write_csv(figure_3_path, rows_3)
            _write_csv(figure_4_path, rows_4)
            scientific_runtime += (
                book_34.runtime_seconds
                + section_34.runtime_seconds
                + eb_34.runtime_seconds
                + reference_34.coupled.runtime_seconds
                + reference_34.straight.runtime_seconds
            )
            scientific_counts["figures_3_4"] = {
                "book_slope_timoshenko_scaled": book_34.scaled_matrix_evaluations,
                "book_slope_timoshenko_raw_quality": book_34.raw_quality_matrix_evaluations,
                "section_clamp_timoshenko_scaled_including_reused_beta0": section_34.scaled_matrix_evaluations,
                "section_clamp_timoshenko_raw_quality_including_reused_beta0": section_34.raw_quality_matrix_evaluations,
                "shared_eb_scaled": eb_34.scaled_matrix_evaluations,
                "shared_eb_raw_quality": eb_34.raw_quality_matrix_evaluations,
                "book_slope_timoshenko_completeness_refined_beta_count": book_34.completeness_refined_beta_count,
                "book_slope_timoshenko_completeness_refinement_attempts": book_34.completeness_refinement_attempts,
                "section_clamp_timoshenko_completeness_refined_beta_count": section_34.completeness_refined_beta_count,
                "section_clamp_timoshenko_completeness_refinement_attempts": section_34.completeness_refinement_attempts,
                "shared_eb_completeness_refined_beta_count": eb_34.completeness_refined_beta_count,
                "shared_eb_completeness_refinement_attempts": eb_34.completeness_refinement_attempts,
                "reference_scaled": references["figures_3_4_hms_dx_209"]["scaled_matrix_evaluations"],
                "reference_raw_quality": references["figures_3_4_hms_dx_209"]["raw_quality_matrix_evaluations"],
                "solver_mode": args.solver_mode,
                "fast_performance_counters": (
                    counters_34.to_dict() if args.solver_mode == "fast" else {}
                ),
            }
        elif figure_3_path.is_file() and figure_4_path.is_file():
            rows_3 = _comparison_rows_as_numbers(_read_csv(figure_3_path), 3)
            rows_4 = _comparison_rows_as_numbers(_read_csv(figure_4_path), 4)
        else:
            rows_3, rows_4 = [], []

        if 2 in selected:
            rows_2 = _comparison_rows_as_numbers(rows_2, 2)
            _validate_saved_comparison_rows(rows_2, beta_values, 2)
            if "figure_2_book_material" not in references:
                raise RuntimeError("saved Figure 2 data lack section-clamp reference evidence")
        if comparison_selected:
            rows_3 = _comparison_rows_as_numbers(rows_3, 3)
            rows_4 = _comparison_rows_as_numbers(rows_4, 4)
            _validate_saved_comparison_rows(rows_3, beta_values, 3)
            _validate_saved_comparison_rows(rows_4, beta_values, 4)
            if not eb_arrays_are_identical(rows_3, rows_4):
                raise RuntimeError("saved EB arrays for Figures 3 and 4 differ")
            if "figures_3_4_hms_dx_209" not in references:
                raise RuntimeError("saved Figure 4 data lack section-clamp reference evidence")

        extended_selected = selected & {5, 6, 7, 8}
        extended_metadata: dict[str, Any] = dict(
            old_manifest.get("extended_fast_solver", {})
        )
        extended_rows: dict[int, list[dict[str, Any]]] = {}
        if extended_selected:
            missing = {
                figure
                for figure in extended_selected
                if not (output_dir / DATA_FILENAMES[figure]).is_file()
            }
            need_compute_extended = bool(
                args.force_recompute or missing
            )
            if need_compute_extended:
                if args.reuse_data:
                    missing_paths = [
                        str(output_dir / DATA_FILENAMES[figure])
                        for figure in sorted(missing)
                    ]
                    raise FileNotFoundError(
                        "--reuse-data requires saved extended CSV: "
                        + ", ".join(missing_paths)
                    )
                if args.solver_mode != "fast":
                    raise RuntimeError(
                        "Figures 5-8 may be computed only with --solver-mode fast"
                    )
                to_compute = set(extended_selected) if args.force_recompute else missing
                built, metadata, extended_runtime = _compute_extended_figure_data(
                    output_dir,
                    beta_values,
                    to_compute,
                    resume=args.resume,
                )
                for figure, rows in built.items():
                    _write_csv(output_dir / DATA_FILENAMES[figure], rows)
                extended_metadata = metadata
                scientific_runtime += extended_runtime
                scientific_counts["extended_fast_solver"] = metadata[
                    "current_run_performance_counters"
                ]
            for figure in sorted(extended_selected):
                path = output_dir / DATA_FILENAMES[figure]
                if not path.is_file():
                    raise FileNotFoundError(path)
                rows = _extended_rows_as_numbers(_read_csv(path), figure)
                _validate_saved_extended_rows(rows, beta_values, figure)
                extended_rows[figure] = rows

        plotting_started = time.perf_counter()
        output_paths: list[Path] = []
        y_limits: dict[str, list[float]] = {}
        if 1 in selected:
            figure_1 = create_figure_1(rows_1)
            output_paths.extend(
                _save_figure(figure_1, output_dir, FIGURE_BASENAMES[1])
            )
        if 2 in selected:
            ylim_2 = _own_ylim(rows_2, "book_slope_lambda", "section_clamp_lambda")
            figure_2 = create_comparison_figure(
                rows_2,
                solid_key="book_slope_lambda",
                dashed_key="section_clamp_lambda",
                ylim=ylim_2,
            )
            output_paths.extend(
                _save_figure(figure_2, output_dir, FIGURE_BASENAMES[2])
            )
            y_limits["figure_2"] = list(ylim_2)
        if comparison_selected:
            shared_ylim = comparison_ylim(rows_3, rows_4)
            y_limits["figures_3_4_shared"] = list(shared_ylim)
            if 3 in selected:
                figure_3 = create_comparison_figure(
                    rows_3,
                    solid_key="timoshenko_lambda",
                    dashed_key="eb_lambda",
                    ylim=shared_ylim,
                )
                output_paths.extend(
                    _save_figure(figure_3, output_dir, FIGURE_BASENAMES[3])
                )
            if 4 in selected:
                figure_4 = create_comparison_figure(
                    rows_4,
                    solid_key="timoshenko_lambda",
                    dashed_key="eb_lambda",
                    ylim=shared_ylim,
                )
                output_paths.extend(
                    _save_figure(figure_4, output_dir, FIGURE_BASENAMES[4])
                )
        if extended_selected & {5, 6}:
            pair_rows = [
                row
                for figure in (5, 6)
                for row in extended_rows.get(figure, [])
            ]
            shared_56 = _own_ylim(pair_rows, "left_lambda", "right_lambda")
            y_limits["figures_5_6_shared"] = list(shared_56)
            for figure in sorted(extended_selected & {5, 6}):
                extended_figure = create_comparison_figure(
                    extended_rows[figure],
                    solid_key="left_lambda",
                    dashed_key="right_lambda",
                    ylim=shared_56,
                )
                output_paths.extend(
                    _save_figure(
                        extended_figure, output_dir, FIGURE_BASENAMES[figure]
                    )
                )
        if extended_selected & {7, 8}:
            pair_rows = [
                row
                for figure in (7, 8)
                for row in extended_rows.get(figure, [])
            ]
            shared_78 = _own_ylim(pair_rows, "left_lambda", "right_lambda")
            y_limits["figures_7_8_shared"] = list(shared_78)
            for figure in sorted(extended_selected & {7, 8}):
                extended_figure = create_comparison_figure(
                    extended_rows[figure],
                    solid_key="left_lambda",
                    dashed_key="right_lambda",
                    ylim=shared_78,
                )
                output_paths.extend(
                    _save_figure(
                        extended_figure, output_dir, FIGURE_BASENAMES[figure]
                    )
                )
        plotting_runtime = time.perf_counter() - plotting_started

        diagnostics: dict[str, Any] = {}
        quality: dict[str, Any] = {}
        if rows_2:
            diagnostics["figure_2"] = _relative_diagnostics(
                rows_2, "relative_clamp_difference"
            )
            quality["figure_2"] = _quality_summary(
                rows_2, ("book_slope", "section_clamp")
            )
        if rows_3:
            diagnostics["figure_3"] = _relative_diagnostics(
                rows_3, "relative_theory_difference"
            )
            quality["figure_3"] = _quality_summary(rows_3, ("timoshenko", "eb"))
        if rows_4:
            diagnostics["figure_4"] = _relative_diagnostics(
                rows_4, "relative_theory_difference"
            )
            quality["figure_4"] = _quality_summary(rows_4, ("timoshenko", "eb"))
        for figure, rows in sorted(extended_rows.items()):
            diagnostics[f"figure_{figure}"] = {
                **_relative_diagnostics(rows, "relative_lambda_difference"),
                "minimum_gap": _minimum_relative_neighbor_gap(
                    rows, ("left_frequency_hz", "right_frequency_hz")
                ),
                **_extended_figure_performance(figure, extended_metadata),
                "data_origin_counts": {
                    origin: sum(row["data_origin"] == origin for row in rows)
                    for origin in (
                        "new_fast_solve",
                        "reused_figure_03",
                        "global_fallback",
                    )
                },
            }
            if figure == 7:
                diagnostics[f"figure_{figure}"]["interpretation"] = (
                    "maximum difference is not interpreted as a pure "
                    "EB-vs-Timoshenko truncation error"
                )
            if figure == 8:
                diagnostics[f"figure_{figure}"]["interpretation"] = (
                    "difference measures the combined effect of rotated "
                    "effective properties and Sbar16 coupling inside the "
                    "Chapter-2 model"
                )
            quality[f"figure_{figure}"] = _quality_summary(
                rows, ("left", "right")
            )

        figures_manifest = []
        for number in sorted(selected):
            figures_manifest.append(
                {
                    "number": number,
                    "basename": FIGURE_BASENAMES[number],
                    "pdf": str((output_dir / f"{FIGURE_BASENAMES[number]}.pdf").relative_to(ROOT)),
                    "png": str((output_dir / f"{FIGURE_BASENAMES[number]}.png").relative_to(ROOT)),
                    "data_csv": str((output_dir / DATA_FILENAMES[number]).relative_to(ROOT)),
                    "figure_size_inches": list(
                        FIGURE_1_SIZE_IN if number == 1 else COMPARISON_FIGURE_SIZE_IN
                    ),
                    "png_dpi": PNG_DPI,
                    "legend_present": False,
                }
            )
        manifest = {
            "workflow": "Yartsev Chapter-2 rectangular anisotropic supervisor figures",
            "workflow_status": "PASS",
            "fast_solver_status": (
                "PASS"
                if (output_dir / FAST_BENCHMARK_FILENAME).is_file()
                and json.loads(
                    (output_dir / FAST_BENCHMARK_FILENAME).read_text(
                        encoding="utf-8"
                    )
                ).get("status") == "FAST_SOLVER_PASS"
                else "FAIL"
            ),
            "extended_figures_status": (
                "PASS"
                if all(
                    (output_dir / DATA_FILENAMES[number]).is_file()
                    and (output_dir / f"{FIGURE_BASENAMES[number]}.pdf").is_file()
                    and (output_dir / f"{FIGURE_BASENAMES[number]}.png").is_file()
                    for number in range(5, 9)
                )
                and extended_metadata.get("reuse_checks", {}).get(
                    "figure_7_eb_equals_figure_3_eb"
                )
                is True
                and extended_metadata.get("reuse_checks", {}).get(
                    "figure_8_theta0_timoshenko_equals_figure_3_timoshenko"
                )
                is True
                else "FAIL"
            ),
            "research_line_separation": {
                "chapter_2_rectangular_anisotropic_only": True,
                "circular_isotropic_article_used": False,
                "circular_section_used": False,
                "isotropic_steel_defaults_used": False,
                "article_workspace_output": False,
            },
            "git_context": _git_context(),
            "figures": figures_manifest,
            "parameters": {
                "figure_1": _figure_1_source_manifest() if rows_1 else {},
                "figure_2": _preset_manifest(FIGURE_2_PRESET),
                "figures_3_4": _preset_manifest(FIGURES_3_4_PRESET),
                "figure_5": _preset_manifest(FIGURE_5_PRESET),
                "figure_6": _preset_manifest(FIGURE_6_PRESET),
                "figure_7": _preset_manifest(FIGURE_7_PRESET),
                "figure_8": _preset_manifest(FIGURE_8_PRESET),
            },
            "lambda_definition": {
                "formula": "(rho*A*omega^2*l^4/(E_x*I_y))^(1/4)",
                "equivalent_formula": "l*(rho*A/(E_x*I_y))^(1/4)*sqrt(omega)",
                "reference_length": "l=(L_1+L_2)/2",
                "reference_area": "A=a*b",
                "reference_second_moment": "I_y=a^3*b/12",
                "elastic_modulus": "E_x at theta=0",
                "omega": "2*pi*f",
                "figures_5_8_fixed_reference": {
                    "a0_m": REFERENCE_A_M,
                    "b0_m": REFERENCE_B_M,
                    "A0_m2": REFERENCE_A_M * REFERENCE_B_M,
                    "I_y0_m4": REFERENCE_A_M**3 * REFERENCE_B_M / 12.0,
                    "formula": "(rho*A0*omega^2*l^4/(E_x0*I_y0))^(1/4)",
                    "E_x0": "elastic E_x at theta=0",
                },
            },
            "beta_grid": {
                "start_deg": float(beta_values[0]),
                "stop_deg": float(beta_values[-1]),
                "requested_step_deg": float(args.beta_step_deg),
                "point_count": len(beta_values),
                "values_deg": beta_values.tolist(),
            },
            "root_contract": {
                "plotted_root_count": PLOTTED_ROOT_COUNT,
                "guard_root_count": GUARD_ROOT_COUNT,
                "curve_identity": "sorted spectral positions at every beta; not tracked modal descendants",
                "interpolation_of_failures": False,
                "ordinary_scan_step_hz": SCAN_STEP_HZ,
                "beta_zero_and_straight_reference_scan_step_hz": BASE_POINT_SCAN_STEP_HZ,
                "completeness_refinement_scan_steps_hz": list(COMPLETENESS_REFINEMENT_STEPS_HZ),
                "close_root_hint_gap_hz": CLOSE_ROOT_HINT_GAP_HZ,
                "local_hint_max_scan_step_hz": LOCAL_HINT_MAX_SCAN_STEP_HZ,
                "local_hint_fallback_scan_step_hz": LOCAL_HINT_FALLBACK_SCAN_STEP_HZ,
                "adjacent_sorted_spectrum_jump_trigger": SPECTRAL_STEP_JUMP_TRIGGER,
                "refinement_role": "numerical missing-root detection only; final arrays remain independently sorted at each beta",
            },
            "root_quality_thresholds": {
                "normalized_determinant_residual": ROOT_DETERMINANT_TOLERANCE,
                "relative_singular_residual": ROOT_SINGULAR_TOLERANCE,
                "acceptance": "scaled PASS or normalized physical-raw PASS; root status must not start with rejected",
            },
            "line_style_meaning": {
                "figure_2_solid": "book_slope_clamp Chapter-2 Timoshenko/generalized torsion",
                "figure_2_dashed": "timoshenko_section_clamp Chapter-2 Timoshenko/generalized torsion",
                "figures_3_4_solid": "state_corrected Timoshenko with generalized rectangular torsion",
                "figures_3_4_dashed": "rectangular orthotropic Euler-Bernoulli with Saint-Venant torsion",
                "figures_5_6_solid": "Chapter-2 state_corrected Timoshenko with generalized torsion",
                "figures_5_6_dashed": "rectangular orthotropic Euler-Bernoulli with Saint-Venant torsion",
                "figure_7_solid": "full Chapter-2 monoclinic Timoshenko theta=5 deg",
                "figure_7_dashed": "reused Figure-3 rectangular orthotropic EB theta=0 approximation",
                "figure_8_solid": "Chapter-2 Timoshenko theta=15 deg",
                "figure_8_dashed": "reused Figure-3 Chapter-2 Timoshenko theta=0 deg",
                "solid_linewidth": SOLID_LINEWIDTH,
                "dashed_linewidth": DASHED_LINEWIDTH,
                "dashed_pattern": list(DASH_PATTERN),
                "mode_colors": list(MODE_COLORS),
            },
            "legends_present": False,
            "y_limits": y_limits,
            "diagnostics": diagnostics,
            "root_quality": quality,
            "section_clamp_reference": references,
            "eb_arrays_figures_3_4_exactly_identical": (
                eb_arrays_are_identical(rows_3, rows_4) if rows_3 and rows_4 else None
            ),
            "matrix_evaluation_counts": scientific_counts,
            "fast_solver_benchmark": (
                json.loads(
                    (output_dir / FAST_BENCHMARK_FILENAME).read_text(
                        encoding="utf-8"
                    )
                )
                if (output_dir / FAST_BENCHMARK_FILENAME).is_file()
                else {}
            ),
            "extended_fast_solver": extended_metadata,
            "runtimes_seconds": {
                "scientific_total": scientific_runtime
                if scientific_runtime
                else float(old_manifest.get("runtimes_seconds", {}).get("scientific_total", 0.0)),
                "plotting_and_reporting": plotting_runtime,
            },
            "output_directory": str(output_dir.relative_to(ROOT)),
            "output_paths": [str(path.relative_to(ROOT)) for path in output_paths]
            + [
                str((output_dir / DATA_FILENAMES[number]).relative_to(ROOT))
                for number in sorted(selected)
            ]
            + [
                str((output_dir / "plot_manifest.json").relative_to(ROOT)),
                str((output_dir / "report.md").relative_to(ROOT)),
            ],
            "execution_mode": (
                "reuse_data" if args.reuse_data else "force_recompute" if args.force_recompute else "auto"
            ),
            "solver_mode": args.solver_mode,
            "jobs": args.jobs,
        }
        (output_dir / "plot_manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        (output_dir / "report.md").write_text(
            _report_text(manifest), encoding="utf-8"
        )
        stale_failure_path = output_dir / "failure_diagnostics.csv"
        if stale_failure_path.is_file():
            stale_failure_path.unlink()
    except Exception as error:
        status = (
            "FAIL_SECTION_CLAMP_REFERENCE"
            if "FAIL_SECTION_CLAMP_REFERENCE" in str(error)
            else "FAIL"
        )
        _write_failure(output_dir, status, str(error))
        print(f"Supervisor figure workflow: {status}", file=sys.stderr)
        print(str(error), file=sys.stderr)
        return 1

    print(f"output_dir={output_dir}")
    print(f"figures={','.join(str(value) for value in sorted(selected))}")
    print(f"beta_points={len(beta_values)}")
    print("Supervisor figure workflow: PASS")
    if extended_selected:
        print(
            "Extended supervisor figures 5–8: "
            f"{manifest['extended_figures_status']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
