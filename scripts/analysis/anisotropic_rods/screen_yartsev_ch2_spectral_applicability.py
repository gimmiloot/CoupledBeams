"""AP-0 small-grid spectral-applicability screening at theta=2 degrees.

This diagnostic-only entry point is isolated from the circular-isotropic
article workflows.  It reuses the validated Chapter-2 supervisor builders,
fast sorted-spectrum coordinator, quality gates, and fixed rectangular
normalization.  Importing this module performs no calculation or plotting.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis.anisotropic_rods import (  # noqa: E402
    plot_yartsev_ch2_supervisor_figures as supervisor,
)
from scripts.lib.yartsev_ch2_fast_beta_sweep import (  # noqa: E402
    FAST_SOLVER_VERSION,
    ExactTransferLRU,
    FastSweepSettings,
    PerformanceCounters,
    load_family_checkpoint,
    run_fast_beta_sweep,
    write_family_checkpoint,
)
from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
    Geometry,
    hms_dx_209_material,
)


DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_spectral_applicability_screening"
    / "theta2_small_grid"
)
DEFAULT_SMOKE_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_spectral_applicability_screening"
    / "theta2_smoke"
)
AP1_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_spectral_applicability_screening"
    / "theta5_small_grid"
)
AP1_SMOKE_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_spectral_applicability_screening"
    / "theta5_smoke"
)
AP2_OUTPUT_DIRS = {
    theta: (
        ROOT
        / "results"
        / "anisotropic_rods"
        / "yartsev_ch2_spectral_applicability_screening"
        / f"theta{theta}_small_grid"
    )
    for theta in (3, 4)
}
AP2_SMOKE_OUTPUT_DIRS = {
    theta: (
        ROOT
        / "results"
        / "anisotropic_rods"
        / "yartsev_ch2_spectral_applicability_screening"
        / f"theta{theta}_smoke"
    )
    for theta in (3, 4)
}
SAMPLED_THETA_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_spectral_applicability_screening"
    / "theta2_theta3_theta4_theta5_comparison"
)
AP0_OUTPUT_DIR = DEFAULT_OUTPUT_DIR
SUPERVISOR_OUTPUT_DIR = supervisor.DEFAULT_OUTPUT_DIR

MU_VALUES = (0.0, 0.25, 0.5)
TAU_VALUES = (-0.2, 0.0, 0.2)
MODEL_IDS = ("T2", "T0", "EB0")
MODEL_THETA_DEG = {"T2": 2.0, "T0": 0.0, "EB0": 0.0}
AP1_MODEL_IDS = ("T5", "T0", "EB0")
AP1_MODEL_THETA_DEG = {"T5": 5.0, "T0": 0.0, "EB0": 0.0}
AP2_MODEL_IDS = {
    3: ("T3", "T0", "EB0"),
    4: ("T4", "T0", "EB0"),
}
AP2_MODEL_THETA_DEG = {"T3": 3.0, "T4": 4.0}
MODEL_SOLVER_PATH = {
    "T2": "book_slope_clamp",
    "T3": "book_slope_clamp",
    "T4": "book_slope_clamp",
    "T5": "book_slope_clamp",
    "T0": "book_slope_clamp",
    "EB0": "rectangular_eb_saint_venant",
}
CLAMP = "book_slope_clamp"

KAPPA = 4.0
REFERENCE_LENGTH_M = 0.4
REFERENCE_A_M = 0.005
REFERENCE_B_M = 0.020
REFERENCE_AREA_M2 = REFERENCE_A_M * REFERENCE_B_M
REFERENCE_IY_M4 = REFERENCE_A_M**3 * REFERENCE_B_M / 12.0
REFERENCE_EX_PA = 191.0e9
MATERIAL_DENSITY_KG_M3 = 1580.0
REFERENCE_VOLUME_M3 = 2.0 * REFERENCE_AREA_M2 * REFERENCE_LENGTH_M
REFERENCE_MASS_KG = MATERIAL_DENSITY_KG_M3 * REFERENCE_VOLUME_M3

FULL_BETA_VALUES_DEG = tuple(float(value) for value in range(0, 91, 5))
SMOKE_BETA_VALUES_DEG = (0.0, 45.0, 90.0)
GUARD_ROOT_COUNT = 7
COMPARED_ROOT_COUNT = 6
APPLICABILITY_THRESHOLD = 0.10
VOLUME_RELATIVE_TOLERANCE = 1.0e-12
SECTION_RATIO_TOLERANCE = 1.0e-14
ROOT_DETERMINANT_TOLERANCE = supervisor.ROOT_DETERMINANT_TOLERANCE
ROOT_SINGULAR_TOLERANCE = supervisor.ROOT_SINGULAR_TOLERANCE

GEOMETRY_FILENAME = "geometry_manifest.csv"
SPECTRA_FILENAME = "screening_spectra.csv"
POINT_METRICS_FILENAME = "screening_point_metrics.csv"
GEOMETRY_SUMMARY_FILENAME = "screening_geometry_summary.csv"
REUSE_AUDIT_FILENAME = "screening_reuse_audit.csv"
CANDIDATES_FILENAME = "candidate_cases.json"
SUMMARY_FILENAME = "screening_summary.json"
REPORT_FILENAME = "screening_report.md"
FIGURE_S1_BASENAME = "screening_delta_simpl_theta2_heatmap"
FIGURE_S2_BASENAME = "screening_gap_opening_theta2_heatmap"
FIGURE_S3_BASENAME = "screening_baseline_gap_vs_orientation_effect"
PNG_DPI = 300
FIGURE_SIZE_IN = (6.4, 4.8)

AP0_GEOMETRY_PATH = AP0_OUTPUT_DIR / GEOMETRY_FILENAME
AP0_SPECTRA_PATH = AP0_OUTPUT_DIR / SPECTRA_FILENAME
AP0_POINT_METRICS_PATH = AP0_OUTPUT_DIR / POINT_METRICS_FILENAME
AP0_GEOMETRY_SUMMARY_PATH = AP0_OUTPUT_DIR / GEOMETRY_SUMMARY_FILENAME
AP0_SUMMARY_PATH = AP0_OUTPUT_DIR / SUMMARY_FILENAME
FIGURE_07_PATH = SUPERVISOR_OUTPUT_DIR / supervisor.DATA_FILENAMES[7]

PAIRWISE_FILENAME = "screening_pairwise_gap_metrics.csv"
POINT_COMPARISON_FILENAME = "theta2_theta5_point_comparison.csv"
GEOMETRY_COMPARISON_FILENAME = "theta2_theta5_geometry_comparison.csv"
AP1_FIGURE_S1_BASENAME = "screening_delta_simpl_theta5_heatmap"
AP1_FIGURE_S2_BASENAME = "screening_delta_simpl_theta2_vs_theta5"
AP1_FIGURE_S3_BASENAME = "screening_pair_gap_vs_change_theta5"

SAMPLED_POINT_FILENAME = "sampled_theta_point_comparison.csv"
SAMPLED_GEOMETRY_FILENAME = "sampled_theta_geometry_comparison.csv"
SAMPLED_COUNTS_FILENAME = "sampled_theta_classification_counts.csv"
SAMPLED_PAIRWISE_FILENAME = "sampled_theta_pairwise_gap_metrics.csv"
SAMPLED_REUSE_AUDIT_FILENAME = "sampled_theta_reuse_audit.csv"
SAMPLED_SUMMARY_FILENAME = "sampled_theta_summary.json"
SAMPLED_REPORT_FILENAME = "sampled_theta_report.md"
AP2_FIGURE_S3_BASENAME = "screening_delta_simpl_sampled_theta_family_curves"
AP2_FIGURE_S4_BASENAME = "screening_sampled_theta_exceedance_counts"
SAMPLED_THETAS_DEG = (2.0, 3.0, 4.0, 5.0)
SAMPLED_NONDECREASING_ABS_TOLERANCE = 1.0e-12

AP0_PRESERVATION_FILENAMES = (
    GEOMETRY_FILENAME,
    SPECTRA_FILENAME,
    POINT_METRICS_FILENAME,
    GEOMETRY_SUMMARY_FILENAME,
    REUSE_AUDIT_FILENAME,
    CANDIDATES_FILENAME,
    SUMMARY_FILENAME,
    REPORT_FILENAME,
    f"{FIGURE_S1_BASENAME}.pdf",
    f"{FIGURE_S1_BASENAME}.png",
    f"{FIGURE_S2_BASENAME}.pdf",
    f"{FIGURE_S2_BASENAME}.png",
    f"{FIGURE_S3_BASENAME}.pdf",
    f"{FIGURE_S3_BASENAME}.png",
)

AP1_PRESERVATION_FILENAMES = (
    GEOMETRY_FILENAME,
    SPECTRA_FILENAME,
    POINT_METRICS_FILENAME,
    GEOMETRY_SUMMARY_FILENAME,
    PAIRWISE_FILENAME,
    POINT_COMPARISON_FILENAME,
    GEOMETRY_COMPARISON_FILENAME,
    REUSE_AUDIT_FILENAME,
    CANDIDATES_FILENAME,
    SUMMARY_FILENAME,
    REPORT_FILENAME,
    f"{AP1_FIGURE_S1_BASENAME}.pdf",
    f"{AP1_FIGURE_S1_BASENAME}.png",
    f"{AP1_FIGURE_S2_BASENAME}.pdf",
    f"{AP1_FIGURE_S2_BASENAME}.png",
    f"{AP1_FIGURE_S3_BASENAME}.pdf",
    f"{AP1_FIGURE_S3_BASENAME}.png",
)

SOURCE_PATHS = (
    ROOT / "scripts" / "lib" / "yartsev_ch2_monoclinic_rod.py",
    ROOT / "scripts" / "lib" / "yartsev_ch2_coupled_rods.py",
    ROOT / "scripts" / "lib" / "yartsev_ch2_rectangular_eb.py",
    ROOT / "scripts" / "lib" / "yartsev_ch2_fast_beta_sweep.py",
    ROOT
    / "scripts"
    / "analysis"
    / "anisotropic_rods"
    / "plot_yartsev_ch2_supervisor_figures.py",
    Path(__file__).resolve(),
)


@dataclass(frozen=True)
class ScreeningConfiguration:
    target_theta_deg: float
    stage_id: str
    primary_model_id: str
    model_ids: tuple[str, str, str]
    output_dir: Path
    smoke_output_dir: Path

    @property
    def is_ap1(self) -> bool:
        return self.stage_id == "AP-1"

    @property
    def is_ap2(self) -> bool:
        return self.stage_id == "AP-2"


def screening_configuration(target_theta_deg: float) -> ScreeningConfiguration:
    theta = float(target_theta_deg)
    if theta == 2.0:
        return ScreeningConfiguration(
            target_theta_deg=2.0,
            stage_id="AP-0",
            primary_model_id="T2",
            model_ids=MODEL_IDS,
            output_dir=DEFAULT_OUTPUT_DIR,
            smoke_output_dir=DEFAULT_SMOKE_OUTPUT_DIR,
        )
    if theta == 5.0:
        return ScreeningConfiguration(
            target_theta_deg=5.0,
            stage_id="AP-1",
            primary_model_id="T5",
            model_ids=AP1_MODEL_IDS,
            output_dir=AP1_OUTPUT_DIR,
            smoke_output_dir=AP1_SMOKE_OUTPUT_DIR,
        )
    if theta in (3.0, 4.0):
        theta_int = int(theta)
        return ScreeningConfiguration(
            target_theta_deg=theta,
            stage_id="AP-2",
            primary_model_id=f"T{theta_int}",
            model_ids=AP2_MODEL_IDS[theta_int],
            output_dir=AP2_OUTPUT_DIRS[theta_int],
            smoke_output_dir=AP2_SMOKE_OUTPUT_DIRS[theta_int],
        )
    raise ValueError("screening target theta must be exactly 2, 3, 4, or 5 degrees")


def model_theta_deg(model_id: str) -> float:
    if model_id in MODEL_THETA_DEG:
        return MODEL_THETA_DEG[model_id]
    if model_id in AP1_MODEL_THETA_DEG:
        return AP1_MODEL_THETA_DEG[model_id]
    if model_id in AP2_MODEL_THETA_DEG:
        return AP2_MODEL_THETA_DEG[model_id]
    raise ValueError(f"unsupported screening model: {model_id}")


@dataclass(frozen=True)
class GeometryCase:
    geometry_id: str
    mu: float
    tau: float
    kappa: float
    L1_m: float
    L2_m: float
    a1_m: float
    a2_m: float
    b1_m: float
    b2_m: float
    A1_m2: float
    A2_m2: float
    Iy1_m4: float
    Iy2_m4: float
    Ip1_m4: float
    Ip2_m4: float
    volume_m3: float
    volume_relative_residual: float
    mass_kg: float
    a2_over_a1: float
    L1_over_a1: float
    L2_over_a2: float
    L1_over_b1: float
    L2_over_b2: float


@dataclass(frozen=True)
class ReuseSpec:
    mu: float
    tau: float
    model_id: str
    figure_number: int
    prefix: str
    data_origin: str

    @property
    def source_path(self) -> Path:
        return SUPERVISOR_OUTPUT_DIR / supervisor.DATA_FILENAMES[self.figure_number]


REUSE_SPECS = (
    ReuseSpec(0.0, 0.0, "T2", 10, "left", "reused_figure_10"),
    ReuseSpec(0.0, 0.0, "T0", 3, "timoshenko", "reused_figure_03"),
    ReuseSpec(0.0, 0.0, "EB0", 3, "eb", "reused_figure_03"),
    ReuseSpec(0.25, 0.0, "T0", 5, "left", "reused_figure_05"),
    ReuseSpec(0.25, 0.0, "EB0", 5, "right", "reused_figure_05"),
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--target-theta-deg",
        type=float,
        choices=(2.0, 3.0, 4.0, 5.0),
        default=2.0,
    )
    parser.add_argument(
        "--combine-sampled-theta",
        action="store_true",
        help="build the solver-free accepted theta=2,3,4,5 comparison",
    )
    data_group = parser.add_mutually_exclusive_group()
    data_group.add_argument("--reuse-data", action="store_true")
    data_group.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args(argv)


def beta_grid(*, smoke: bool = False) -> NDArray[np.float64]:
    return np.asarray(
        SMOKE_BETA_VALUES_DEG if smoke else FULL_BETA_VALUES_DEG, dtype=float
    )


def _number_token(value: float) -> str:
    sign = "m" if value < 0.0 else "p"
    magnitude = f"{abs(float(value)):.2f}".replace(".", "p")
    return f"{sign}{magnitude}"


def geometry_id(mu: float, tau: float) -> str:
    return f"mu_{_number_token(mu)}_tau_{_number_token(tau)}"


def make_geometry(mu: float, tau: float) -> GeometryCase:
    if mu not in MU_VALUES or tau not in TAU_VALUES:
        raise ValueError("geometry lies outside the fixed AP-0 grid")
    L1 = REFERENCE_LENGTH_M * (1.0 - mu)
    L2 = REFERENCE_LENGTH_M * (1.0 + mu)
    scale = REFERENCE_A_M / math.sqrt(1.0 + tau**2 + 2.0 * mu * tau)
    a1 = scale * (1.0 - tau)
    a2 = scale * (1.0 + tau)
    b1 = KAPPA * a1
    b2 = KAPPA * a2
    point_geometry_1 = Geometry(a=a1, b=b1, length=L1)
    point_geometry_2 = Geometry(a=a2, b=b2, length=L2)
    volume = point_geometry_1.area * L1 + point_geometry_2.area * L2
    relative_residual = abs(volume - REFERENCE_VOLUME_M3) / REFERENCE_VOLUME_M3
    values = (L1, L2, a1, a2, b1, b2)
    if any(not math.isfinite(value) or value <= 0.0 for value in values):
        raise RuntimeError("AP-0 geometry contains a non-positive dimension")
    if relative_residual > VOLUME_RELATIVE_TOLERANCE:
        raise RuntimeError("AP-0 volume-preservation gate failed")
    if abs(b1 / a1 - KAPPA) > SECTION_RATIO_TOLERANCE or abs(
        b2 / a2 - KAPPA
    ) > SECTION_RATIO_TOLERANCE:
        raise RuntimeError("AP-0 similar-section ratio gate failed")
    if not math.isclose(L1 + L2, 0.8, rel_tol=0.0, abs_tol=2.0e-16):
        raise RuntimeError("AP-0 total-length gate failed")
    return GeometryCase(
        geometry_id=geometry_id(mu, tau),
        mu=float(mu),
        tau=float(tau),
        kappa=KAPPA,
        L1_m=L1,
        L2_m=L2,
        a1_m=a1,
        a2_m=a2,
        b1_m=b1,
        b2_m=b2,
        A1_m2=point_geometry_1.area,
        A2_m2=point_geometry_2.area,
        Iy1_m4=point_geometry_1.I_y,
        Iy2_m4=point_geometry_2.I_y,
        Ip1_m4=point_geometry_1.I_p,
        Ip2_m4=point_geometry_2.I_p,
        volume_m3=volume,
        volume_relative_residual=relative_residual,
        mass_kg=MATERIAL_DENSITY_KG_M3 * volume,
        a2_over_a1=a2 / a1,
        L1_over_a1=L1 / a1,
        L2_over_a2=L2 / a2,
        L1_over_b1=L1 / b1,
        L2_over_b2=L2 / b2,
    )


def screening_geometries(*, smoke: bool = False) -> list[GeometryCase]:
    cases = [make_geometry(mu, tau) for mu in MU_VALUES for tau in TAU_VALUES]
    return cases[:3] if smoke else cases


def validate_geometry_manifest(
    cases: Sequence[GeometryCase], *, smoke: bool
) -> None:
    canonical = screening_geometries(smoke=smoke)
    if len(cases) != len(canonical):
        raise RuntimeError("AP-0 geometry manifest has the wrong row count")
    expected_by_id = {case.geometry_id: case for case in canonical}
    if {case.geometry_id for case in cases} != set(expected_by_id):
        raise RuntimeError("AP-0 geometry manifest has missing/extra cases")
    for case in cases:
        expected = expected_by_id[case.geometry_id]
        for field in GeometryCase.__dataclass_fields__:
            if field == "geometry_id":
                continue
            if not math.isclose(
                float(getattr(case, field)),
                float(getattr(expected, field)),
                rel_tol=2.0e-15,
                abs_tol=2.0e-16,
            ):
                raise RuntimeError(
                    f"{case.geometry_id}: geometry manifest field {field} changed"
                )
        if case.volume_relative_residual > VOLUME_RELATIVE_TOLERANCE:
            raise RuntimeError(f"{case.geometry_id}: volume gate failed")


def model_preset(case: GeometryCase, model_id: str) -> supervisor.FigurePreset:
    if model_id not in (*MODEL_IDS, "T3", "T4", "T5"):
        raise ValueError(f"unsupported AP-0 model: {model_id}")
    theta = model_theta_deg(model_id)
    return supervisor.FigurePreset(
        figure_numbers=(),
        material_name="HMS/DX-209",
        material_factory=hms_dx_209_material,
        a_m=case.a1_m,
        b_m=case.b1_m,
        length_1_m=case.L1_m,
        length_2_m=case.L2_m,
        a_2_m=case.a2_m,
        b_2_m=case.b2_m,
        theta_1_deg=theta,
        theta_2_deg=theta,
        material_mode="elastic",
        mu=case.mu,
    )


def fixed_lambda_from_frequency(
    frequency_hz: float | NDArray[np.float64],
) -> float | NDArray[np.float64]:
    frequency = np.asarray(frequency_hz, dtype=float)
    value = supervisor.lambda_from_omega(
        2.0 * np.pi * frequency,
        rho_kg_m3=MATERIAL_DENSITY_KG_M3,
        a_m=REFERENCE_A_M,
        b_m=REFERENCE_B_M,
        length_1_m=REFERENCE_LENGTH_M,
        length_2_m=REFERENCE_LENGTH_M,
        elastic_ex_pa=REFERENCE_EX_PA,
    )
    return value


def normalization_contract(
    case: GeometryCase,
    model_id: str,
    *,
    source_hashes: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    contract = {
        "formula": "(rho*A0*omega^2*l^4/(E_x0*I_y0))^(1/4)",
        "rho_kg_m3": MATERIAL_DENSITY_KG_M3,
        "a0_m": REFERENCE_A_M,
        "b0_m": REFERENCE_B_M,
        "A0_m2": REFERENCE_AREA_M2,
        "I_y0_m4": REFERENCE_IY_M4,
        "l_m": REFERENCE_LENGTH_M,
        "E_x0_pa": REFERENCE_EX_PA,
        "fixed_for_all_geometries_and_models": True,
        "screening_fingerprint_context": {
            "model_id": model_id,
            "theta_deg": model_theta_deg(model_id),
            "mu": case.mu,
            "tau": case.tau,
            "clamp": CLAMP,
            "root_count": GUARD_ROOT_COUNT,
            "quality_determinant_tolerance": ROOT_DETERMINANT_TOLERANCE,
            "quality_singular_tolerance": ROOT_SINGULAR_TOLERANCE,
        },
    }
    if source_hashes is not None:
        contract["screening_fingerprint_context"]["source_sha256"] = dict(
            source_hashes
        )
    return contract


def normalized_neighbor_gaps(values: Sequence[float]) -> NDArray[np.float64]:
    array = np.asarray(values, dtype=float)
    if array.shape != (GUARD_ROOT_COUNT,):
        raise ValueError("gap calculation requires seven roots")
    if np.any(~np.isfinite(array)) or np.any(array <= 0.0) or np.any(np.diff(array) <= 0.0):
        raise ValueError("gap calculation requires finite strictly sorted roots")
    return 2.0 * np.diff(array) / (array[1:] + array[:-1])


def point_metrics_from_lambdas(
    t2: Sequence[float], t0: Sequence[float], eb0: Sequence[float]
) -> dict[str, Any]:
    T2 = np.asarray(t2, dtype=float)
    T0 = np.asarray(t0, dtype=float)
    EB0 = np.asarray(eb0, dtype=float)
    for values in (T2, T0, EB0):
        if values.shape != (GUARD_ROOT_COUNT,):
            raise ValueError("point metrics require seven roots for every model")
    delta_beam = np.abs(T0[:COMPARED_ROOT_COUNT] - EB0[:COMPARED_ROOT_COUNT]) / T0[:COMPARED_ROOT_COUNT]
    delta_orient = np.abs(T2[:COMPARED_ROOT_COUNT] - T0[:COMPARED_ROOT_COUNT]) / T2[:COMPARED_ROOT_COUNT]
    delta_simpl = np.abs(T2[:COMPARED_ROOT_COUNT] - EB0[:COMPARED_ROOT_COUNT]) / T2[:COMPARED_ROOT_COUNT]
    gaps_0 = normalized_neighbor_gaps(T0)
    gaps_2 = normalized_neighbor_gaps(T2)
    baseline_pair_index = int(np.argmin(gaps_0))
    opening = gaps_2 - gaps_0
    open_index = int(np.argmax(opening))

    def maximum_with_mode(values: NDArray[np.float64]) -> tuple[float, int]:
        index = int(np.argmax(values))
        return float(values[index]), index + 1

    Delta_beam, beam_mode = maximum_with_mode(delta_beam)
    Delta_orient, orient_mode = maximum_with_mode(delta_orient)
    Delta_simpl, simpl_mode = maximum_with_mode(delta_simpl)
    return {
        "Delta_beam": Delta_beam,
        "Delta_beam_mode": beam_mode,
        "Delta_orient": Delta_orient,
        "Delta_orient_mode": orient_mode,
        "Delta_simpl": Delta_simpl,
        "Delta_simpl_mode": simpl_mode,
        "g_min_0": float(gaps_0[baseline_pair_index]),
        "baseline_min_gap_pair": baseline_pair_index + 1,
        "g_min_2": float(np.min(gaps_2)),
        "G_nearest": float(opening[baseline_pair_index]),
        "G_open": float(opening[open_index]),
        "G_open_pair": open_index + 1,
        "G_close": float(np.min(opening)),
        "classification_at_point": (
            "WITHIN_10_PERCENT_ON_SCREENING_GRID"
            if Delta_simpl <= APPLICABILITY_THRESHOLD
            else "EXCEEDS_10_PERCENT_ON_SCREENING_GRID"
        ),
    }


def _write_csv(path: Path, rows: Iterable[Mapping[str, Any]]) -> None:
    supervisor._write_csv(path, rows)


def _read_csv(path: Path) -> list[dict[str, str]]:
    return supervisor._read_csv(path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def supervisor_data_hashes() -> dict[str, str]:
    return {
        supervisor.DATA_FILENAMES[number]: _sha256(
            SUPERVISOR_OUTPUT_DIR / supervisor.DATA_FILENAMES[number]
        )
        for number in range(1, 13)
    }


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


def git_context() -> dict[str, Any]:
    return {
        "cwd": str(ROOT),
        "git_root": _git("rev-parse", "--show-toplevel"),
        "branch": _git("branch", "--show-current"),
        "head": _git("rev-parse", "HEAD"),
        "origin_main": _git("rev-parse", "origin/main"),
        "status_short": _git("status", "--short").splitlines(),
    }


def _reuse_spec(case: GeometryCase, model_id: str) -> ReuseSpec | None:
    return next(
        (
            spec
            for spec in REUSE_SPECS
            if spec.mu == case.mu
            and spec.tau == case.tau
            and spec.model_id == model_id
        ),
        None,
    )


def _quality_from_source(
    row: Mapping[str, Any], prefix: str
) -> dict[str, Any]:
    return supervisor._saved_quality(row, prefix)


def _spectrum_row(
    case: GeometryCase,
    *,
    beta_deg: float,
    model_id: str,
    mode: int,
    frequency_hz: float,
    lambda_ref: float,
    quality: Mapping[str, Any],
    data_origin: str,
    fallback_reason: str = "",
) -> dict[str, Any]:
    return {
        "geometry_id": case.geometry_id,
        "mu": case.mu,
        "tau": case.tau,
        "beta_deg": beta_deg,
        "model_id": model_id,
        "theta_deg": model_theta_deg(model_id),
        "mode": mode,
        "root_role": "plotted" if mode <= COMPARED_ROOT_COUNT else "guard",
        "frequency_hz": frequency_hz,
        "lambda_ref": lambda_ref,
        "quality_status": quality["quality_status"],
        "quality_basis": quality["quality_basis"],
        "root_status": quality["root_status"],
        "accepted_determinant_residual": quality[
            "accepted_determinant_residual"
        ],
        "accepted_relative_singular_residual": quality[
            "accepted_relative_singular_residual"
        ],
        "data_origin": data_origin,
        "global_fallback_reason": fallback_reason,
    }


def load_reused_family(
    case: GeometryCase,
    model_id: str,
    beta_values: NDArray[np.float64],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    spec = _reuse_spec(case, model_id)
    if spec is None:
        raise ValueError("requested AP-0 family has no approved reuse source")
    source_rows = _read_csv(spec.source_path)
    selected_beta = {float(value) for value in beta_values}
    selected = sorted(
        (
            row
            for row in source_rows
            if float(row["beta_deg"]) in selected_beta
        ),
        key=lambda row: (float(row["beta_deg"]), int(row["mode"])),
    )
    expected_count = len(beta_values) * GUARD_ROOT_COUNT
    if len(selected) != expected_count:
        raise RuntimeError(
            f"reuse source Figure {spec.figure_number} has incomplete beta/root rows"
        )
    rows: list[dict[str, Any]] = []
    source_array: list[tuple[float, int, float, float]] = []
    output_array: list[tuple[float, int, float, float]] = []
    maximum_normalization_residual = 0.0
    for source in selected:
        beta_deg = float(source["beta_deg"])
        mode = int(source["mode"])
        frequency_hz = float(source[f"{spec.prefix}_frequency_hz"])
        lambda_ref = float(source[f"{spec.prefix}_lambda"])
        quality = _quality_from_source(source, spec.prefix)
        if quality["quality_status"] != "PASS":
            raise RuntimeError("approved reuse source contains a rejected root")
        calculated_lambda = float(fixed_lambda_from_frequency(frequency_hz))
        normalization_residual = abs(calculated_lambda - lambda_ref) / max(
            abs(lambda_ref), np.finfo(float).tiny
        )
        maximum_normalization_residual = max(
            maximum_normalization_residual, normalization_residual
        )
        row = _spectrum_row(
            case,
            beta_deg=beta_deg,
            model_id=model_id,
            mode=mode,
            frequency_hz=frequency_hz,
            lambda_ref=lambda_ref,
            quality=quality,
            data_origin=spec.data_origin,
        )
        rows.append(row)
        source_array.append((beta_deg, mode, frequency_hz, lambda_ref))
        output_array.append(
            (
                float(row["beta_deg"]),
                int(row["mode"]),
                float(row["frequency_hz"]),
                float(row["lambda_ref"]),
            )
        )
    exact = bool(
        np.array_equal(np.asarray(source_array), np.asarray(output_array))
    )
    if not exact or maximum_normalization_residual > 5.0e-15:
        raise RuntimeError("approved reuse source normalization/equality gate failed")
    audit = {
        "geometry_id": case.geometry_id,
        "mu": case.mu,
        "tau": case.tau,
        "model_id": model_id,
        "source_figure": spec.figure_number,
        "source_csv": str(spec.source_path.relative_to(ROOT)),
        "source_csv_sha256": _sha256(spec.source_path),
        "data_origin": spec.data_origin,
        "expected_row_count": expected_count,
        "reused_row_count": len(rows),
        "beta_frequency_lambda_exactly_equal": exact,
        "maximum_fixed_normalization_relative_residual": maximum_normalization_residual,
        "geometry_fingerprint": supervisor.stable_fingerprint(case),
        "status": "PASS",
    }
    return rows, audit


def _family_id(
    case: GeometryCase,
    model_id: str,
    *,
    smoke: bool,
    stage_tag: str = "ap0",
) -> str:
    prefix = f"{stage_tag}_smoke" if smoke else stage_tag
    return f"{prefix}_{case.geometry_id}_{model_id.lower()}"


def _run_ap0_fast_family(
    preset: supervisor.FigurePreset,
    beta_values: NDArray[np.float64],
    model: str,
    *,
    transfer_cache: ExactTransferLRU,
    counters: PerformanceCounters,
):
    """Run the validated coordinator with an unconditional global fallback.

    All matrix construction, root finding, local windows, cluster handling,
    and quality finalization remain the existing supervisor implementations.
    AP-0 only supplies the documented fallback policy: if a local/predictor
    path fails, retain the already quality-checked global inventory directly.
    """

    before = counters.copy()
    point_1 = supervisor._point(preset, 1)
    point_2 = supervisor._point(preset, 2)
    factory_cache: dict[float, tuple[Any, Any]] = {}
    global_cache: dict[float, Any] = {}

    def factories(beta_deg: float):
        if beta_deg not in factory_cache:
            factory_cache[beta_deg] = supervisor._cached_factories(
                model,
                beta_deg,
                point_1,
                point_2,
                transfer_cache,
                counters,
            )
        return factory_cache[beta_deg]

    def global_search(beta_deg: float):
        if beta_deg not in global_cache:
            scaled_factory, raw_factory = factories(beta_deg)
            global_cache[beta_deg] = supervisor._solve_cached_global_spectrum(
                point_1, scaled_factory, raw_factory, counters
            )
        return global_cache[beta_deg]

    def scan_interval(beta_deg: float, interval):
        scaled_factory, _ = factories(beta_deg)
        return supervisor._scan_fast_interval(
            point_1, scaled_factory, interval, counters
        )

    def finalize(beta_deg: float, roots):
        _, raw_factory = factories(beta_deg)
        reconciled = list(roots)
        if beta_deg not in FastSweepSettings().mandatory_anchors_deg:
            counters.global_inventory_checks += 1
            candidates = list(global_search(beta_deg).roots) + reconciled
            candidates.sort(key=lambda root: root.frequency_hz)
            reconciled = []
            for candidate in candidates:
                duplicate = any(
                    abs(candidate.frequency_hz - existing.frequency_hz)
                    <= FastSweepSettings().duplicate_relative_tolerance
                    * max(
                        abs(candidate.frequency_hz),
                        abs(existing.frequency_hz),
                        1.0,
                    )
                    for existing in reconciled
                )
                if not duplicate:
                    reconciled.append(candidate)
            reconciled = reconciled[:GUARD_ROOT_COUNT]
        return supervisor._finalize_fast_roots(
            reconciled,
            raw_factory,
            counters,
            scan_step_hz=supervisor.LOCAL_HINT_MAX_SCAN_STEP_HZ,
        )

    result = run_fast_beta_sweep(
        beta_values,
        global_search=global_search,
        scan_interval=scan_interval,
        finalize_local=finalize,
        spectrum_frequencies=lambda spectrum: [
            root.frequency_hz for root in spectrum.roots
        ],
        root_frequency=lambda root: root.frequency_hz,
        fallback_search=lambda beta_deg, _previous, _older: global_search(
            beta_deg
        ),
        settings=FastSweepSettings(root_count=GUARD_ROOT_COUNT),
        counters=counters,
    )
    result.counters = counters.difference(before)
    return result


def _run_or_resume_ap0_family(
    preset: supervisor.FigurePreset,
    beta_values: NDArray[np.float64],
    model: str,
    *,
    family_id: str,
    output_dir: Path,
    resume: bool,
    transfer_cache: ExactTransferLRU,
    counters: PerformanceCounters,
    lambda_normalization: Mapping[str, Any],
):
    checkpoint_dir = output_dir / supervisor.FAST_FAMILY_DIRNAME
    fingerprint = supervisor._fast_family_fingerprint(
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
            return supervisor._fast_family_from_checkpoint(rows, manifest), True
    result = _run_ap0_fast_family(
        preset,
        beta_values,
        model,
        transfer_cache=transfer_cache,
        counters=counters,
    )
    write_family_checkpoint(
        checkpoint_dir,
        family_id,
        supervisor._fast_family_checkpoint_rows(result, family_id),
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
            "ap0_fallback_policy": "quality_checked_global_inventory",
        },
    )
    return result, False


def solve_new_family(
    case: GeometryCase,
    model_id: str,
    beta_values: NDArray[np.float64],
    *,
    output_dir: Path,
    resume: bool,
    transfer_cache: ExactTransferLRU,
    counters: PerformanceCounters,
    smoke: bool,
    stage_tag: str = "ap0",
    source_hashes: Mapping[str, str] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any], float]:
    preset = model_preset(case, model_id)
    solver_path = MODEL_SOLVER_PATH[model_id]
    normalization = normalization_contract(
        case, model_id, source_hashes=source_hashes
    )
    family_id = _family_id(
        case, model_id, smoke=smoke, stage_tag=stage_tag
    )
    started = time.perf_counter()
    result, resumed = _run_or_resume_ap0_family(
        preset,
        beta_values,
        solver_path,
        family_id=family_id,
        output_dir=output_dir,
        resume=resume,
        transfer_cache=transfer_cache,
        counters=counters,
        lambda_normalization=normalization,
    )
    current_runtime = 0.0 if resumed else time.perf_counter() - started
    rows: list[dict[str, Any]] = []
    for beta_value in beta_values:
        beta_deg = float(beta_value)
        spectrum = result.spectra[beta_deg]
        for mode, (root, quality) in enumerate(
            zip(spectrum.roots, spectrum.quality), start=1
        ):
            origin = result.data_origins[beta_deg]
            rows.append(
                _spectrum_row(
                    case,
                    beta_deg=beta_deg,
                    model_id=model_id,
                    mode=mode,
                    frequency_hz=root.frequency_hz,
                    lambda_ref=float(
                        fixed_lambda_from_frequency(root.frequency_hz)
                    ),
                    quality=quality,
                    data_origin=(
                        "global_fallback"
                        if origin == "global_fallback"
                        else "new_fast_solve"
                    ),
                    fallback_reason=result.fallback_reasons.get(beta_deg, ""),
                )
            )
    fingerprint = supervisor._fast_family_fingerprint(
        preset,
        beta_values,
        solver_path,
        lambda_normalization=normalization,
    )
    metadata = {
        "family_id": family_id,
        "geometry_id": case.geometry_id,
        "model_id": model_id,
        "theta_deg": model_theta_deg(model_id),
        "solver_path": solver_path,
        "fingerprint": fingerprint,
        "resumed_from_checkpoint": resumed,
        "recorded_runtime_seconds": result.runtime_seconds,
        "current_runtime_seconds": current_runtime,
        "performance_counters": result.counters.to_dict(),
        "global_fallback_count": len(result.fallback_reasons),
        "fallback_reasons": {
            str(beta): reason
            for beta, reason in sorted(result.fallback_reasons.items())
        },
        "maximum_anchor_relative_error": max(
            result.anchor_relative_errors.values(), default=0.0
        ),
        "normalization": normalization,
    }
    return rows, metadata, current_runtime


def validate_spectra(
    rows: Sequence[Mapping[str, Any]],
    geometries: Sequence[GeometryCase],
    beta_values: NDArray[np.float64],
    *,
    model_ids: Sequence[str] = MODEL_IDS,
) -> None:
    expected_count = (
        len(geometries) * len(beta_values) * len(model_ids) * GUARD_ROOT_COUNT
    )
    if len(rows) != expected_count:
        raise RuntimeError(
            f"screening spectra row count {len(rows)} != {expected_count}"
        )
    by_key: dict[tuple[str, float, str], list[Mapping[str, Any]]] = {}
    for row in rows:
        key = (
            str(row["geometry_id"]),
            float(row["beta_deg"]),
            str(row["model_id"]),
        )
        by_key.setdefault(key, []).append(row)
    expected_groups = len(geometries) * len(beta_values) * len(model_ids)
    if len(by_key) != expected_groups:
        raise RuntimeError("screening spectra contain missing/duplicate groups")
    for key, group in by_key.items():
        ordered = sorted(group, key=lambda row: int(row["mode"]))
        if [int(row["mode"]) for row in ordered] != list(
            range(1, GUARD_ROOT_COUNT + 1)
        ):
            raise RuntimeError(f"{key}: incomplete guard-root inventory")
        frequencies = np.asarray(
            [float(row["frequency_hz"]) for row in ordered], dtype=float
        )
        lambdas = np.asarray(
            [float(row["lambda_ref"]) for row in ordered], dtype=float
        )
        if (
            np.any(~np.isfinite(frequencies))
            or np.any(frequencies <= 0.0)
            or np.any(np.diff(frequencies) <= 0.0)
            or np.any(~np.isfinite(lambdas))
            or np.any(lambdas <= 0.0)
            or np.any(np.diff(lambdas) <= 0.0)
        ):
            raise RuntimeError(f"{key}: invalid or duplicate sorted roots")
        for row in ordered:
            if row["quality_status"] != "PASS" or str(
                row["root_status"]
            ).startswith("rejected"):
                raise RuntimeError(f"{key}: rejected root")
            if (
                float(row["accepted_determinant_residual"])
                > ROOT_DETERMINANT_TOLERANCE
                or float(row["accepted_relative_singular_residual"])
                > ROOT_SINGULAR_TOLERANCE
            ):
                raise RuntimeError(f"{key}: root-quality threshold exceeded")


def spectra_as_numbers(
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    floats = (
        "mu",
        "tau",
        "beta_deg",
        "theta_deg",
        "frequency_hz",
        "lambda_ref",
        "accepted_determinant_residual",
        "accepted_relative_singular_residual",
    )
    result: list[dict[str, Any]] = []
    for source in rows:
        row = dict(source)
        row["mode"] = int(row["mode"])
        for field in floats:
            row[field] = float(row[field])
        result.append(row)
    return result


def reuse_audit_as_numbers(
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    converted: list[dict[str, Any]] = []
    for source in rows:
        row = dict(source)
        for field in ("mu", "tau", "maximum_fixed_normalization_relative_residual"):
            if field in row and row[field] != "":
                row[field] = float(row[field])
        for field in ("source_figure", "expected_row_count", "reused_row_count"):
            if field in row and row[field] != "":
                row[field] = int(row[field])
        if "beta_frequency_lambda_exactly_equal" in row:
            exact = row["beta_frequency_lambda_exactly_equal"]
            row["beta_frequency_lambda_exactly_equal"] = (
                exact
                if isinstance(exact, bool)
                else str(exact).strip().lower() == "true"
            )
        converted.append(row)
    return converted


def build_point_metrics(
    spectra: Sequence[Mapping[str, Any]],
    geometries: Sequence[GeometryCase],
    beta_values: NDArray[np.float64],
) -> list[dict[str, Any]]:
    lookup: dict[tuple[str, float, str], list[Mapping[str, Any]]] = {}
    for row in spectra:
        lookup.setdefault(
            (
                str(row["geometry_id"]),
                float(row["beta_deg"]),
                str(row["model_id"]),
            ),
            [],
        ).append(row)
    rows: list[dict[str, Any]] = []
    for case in geometries:
        for beta_value in beta_values:
            beta_deg = float(beta_value)
            arrays: dict[str, NDArray[np.float64]] = {}
            for model_id in MODEL_IDS:
                selected = sorted(
                    lookup[(case.geometry_id, beta_deg, model_id)],
                    key=lambda row: int(row["mode"]),
                )
                arrays[model_id] = np.asarray(
                    [float(row["lambda_ref"]) for row in selected], dtype=float
                )
            metrics = point_metrics_from_lambdas(
                arrays["T2"], arrays["T0"], arrays["EB0"]
            )
            rows.append(
                {
                    "geometry_id": case.geometry_id,
                    "mu": case.mu,
                    "tau": case.tau,
                    "beta_deg": beta_deg,
                    **metrics,
                }
            )
    return rows


def point_metrics_as_numbers(
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    float_fields = (
        "mu",
        "tau",
        "beta_deg",
        "Delta_beam",
        "Delta_orient",
        "Delta_simpl",
        "g_min_0",
        "g_min_2",
        "G_nearest",
        "G_open",
        "G_close",
    )
    int_fields = (
        "Delta_beam_mode",
        "Delta_orient_mode",
        "Delta_simpl_mode",
        "baseline_min_gap_pair",
        "G_open_pair",
    )
    converted: list[dict[str, Any]] = []
    for source in rows:
        row = dict(source)
        for field in float_fields:
            row[field] = float(row[field])
        for field in int_fields:
            row[field] = int(row[field])
        converted.append(row)
    return converted


def _extreme_row(
    rows: Sequence[Mapping[str, Any]], key: str, *, maximum: bool
) -> Mapping[str, Any]:
    function = max if maximum else min
    return function(rows, key=lambda row: float(row[key]))


def build_geometry_summary(
    point_rows: Sequence[Mapping[str, Any]], geometries: Sequence[GeometryCase]
) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    for case in geometries:
        selected = [
            row for row in point_rows if row["geometry_id"] == case.geometry_id
        ]
        beam = _extreme_row(selected, "Delta_beam", maximum=True)
        orient = _extreme_row(selected, "Delta_orient", maximum=True)
        simpl = _extreme_row(selected, "Delta_simpl", maximum=True)
        baseline_gap = _extreme_row(selected, "g_min_0", maximum=False)
        nearest = _extreme_row(selected, "G_nearest", maximum=True)
        opening = _extreme_row(selected, "G_open", maximum=True)
        Delta_simpl_max = float(simpl["Delta_simpl"])
        summaries.append(
            {
                "geometry_id": case.geometry_id,
                "mu": case.mu,
                "tau": case.tau,
                "Delta_beam_max": float(beam["Delta_beam"]),
                "Delta_beam_max_beta": float(beam["beta_deg"]),
                "Delta_beam_max_mode": int(beam["Delta_beam_mode"]),
                "Delta_orient_max": float(orient["Delta_orient"]),
                "Delta_orient_max_beta": float(orient["beta_deg"]),
                "Delta_orient_max_mode": int(orient["Delta_orient_mode"]),
                "Delta_simpl_max": Delta_simpl_max,
                "Delta_simpl_max_beta": float(simpl["beta_deg"]),
                "Delta_simpl_max_mode": int(simpl["Delta_simpl_mode"]),
                "minimum_baseline_gap": float(baseline_gap["g_min_0"]),
                "minimum_baseline_gap_beta": float(baseline_gap["beta_deg"]),
                "minimum_baseline_gap_pair": int(
                    baseline_gap["baseline_min_gap_pair"]
                ),
                "maximum_G_nearest": float(nearest["G_nearest"]),
                "maximum_G_nearest_beta": float(nearest["beta_deg"]),
                "maximum_G_open": float(opening["G_open"]),
                "maximum_G_open_beta": float(opening["beta_deg"]),
                "maximum_G_open_pair": int(opening["G_open_pair"]),
                "screening_classification": (
                    "WITHIN_10_PERCENT_ON_SCREENING_GRID"
                    if Delta_simpl_max <= APPLICABILITY_THRESHOLD
                    else "EXCEEDS_10_PERCENT_ON_SCREENING_GRID"
                ),
            }
        )
    return summaries


def geometry_summary_as_numbers(
    rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    float_fields = (
        "mu",
        "tau",
        "Delta_beam_max",
        "Delta_beam_max_beta",
        "Delta_orient_max",
        "Delta_orient_max_beta",
        "Delta_simpl_max",
        "Delta_simpl_max_beta",
        "minimum_baseline_gap",
        "minimum_baseline_gap_beta",
        "maximum_G_nearest",
        "maximum_G_nearest_beta",
        "maximum_G_open",
        "maximum_G_open_beta",
    )
    int_fields = (
        "Delta_beam_max_mode",
        "Delta_orient_max_mode",
        "Delta_simpl_max_mode",
        "minimum_baseline_gap_pair",
        "maximum_G_open_pair",
    )
    converted: list[dict[str, Any]] = []
    for source in rows:
        row = dict(source)
        for field in float_fields:
            row[field] = float(row[field])
        for field in int_fields:
            row[field] = int(row[field])
        converted.append(row)
    return converted


def select_candidates(
    summaries: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    stable = _extreme_row(summaries, "Delta_simpl_max", maximum=False)
    sensitive = _extreme_row(summaries, "Delta_simpl_max", maximum=True)
    closest = min(
        summaries,
        key=lambda row: abs(
            float(row["Delta_simpl_max"]) - APPLICABILITY_THRESHOLD
        ),
    )
    smallest_gap = _extreme_row(
        summaries, "minimum_baseline_gap", maximum=False
    )
    largest_opening = _extreme_row(summaries, "maximum_G_open", maximum=True)

    def compact(row: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "geometry_id": row["geometry_id"],
            "mu": float(row["mu"]),
            "tau": float(row["tau"]),
            "Delta_simpl_max": float(row["Delta_simpl_max"]),
            "minimum_baseline_gap": float(row["minimum_baseline_gap"]),
            "maximum_G_open": float(row["maximum_G_open"]),
        }

    closest_distance = abs(
        float(closest["Delta_simpl_max"]) - APPLICABILITY_THRESHOLD
    )
    return {
        "stable_candidate": compact(stable),
        "sensitive_candidate": compact(sensitive),
        "borderline_candidate": (
            compact(closest)
            if closest_distance <= 0.03
            else "NO_CLOSE_BORDERLINE_CASE_FOUND"
        ),
        "borderline_distance_to_0p10": closest_distance,
        "smallest_gap_candidate": compact(smallest_gap),
        "largest_gap_opening_candidate": compact(largest_opening),
        "next_stage_started": False,
    }


def _rankdata(values: Sequence[float]) -> NDArray[np.float64]:
    array = np.asarray(values, dtype=float)
    order = np.argsort(array, kind="mergesort")
    ranks = np.empty(len(array), dtype=float)
    start = 0
    while start < len(array):
        stop = start + 1
        while stop < len(array) and array[order[stop]] == array[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * (start + stop - 1) + 1.0
        start = stop
    return ranks


def descriptive_gap_relation(
    point_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    gap = np.asarray([float(row["g_min_0"]) for row in point_rows])
    orient = np.asarray([float(row["Delta_orient"]) for row in point_rows])
    opening = np.asarray([float(row["G_open"]) for row in point_rows])

    def spearman(left: NDArray[np.float64], right: NDArray[np.float64]) -> float:
        ranked_left = _rankdata(left)
        ranked_right = _rankdata(right)
        return float(np.corrcoef(ranked_left, ranked_right)[0, 1])

    quartile = float(np.quantile(gap, 0.25))
    low = gap <= quartile
    return {
        "point_count": len(point_rows),
        "lowest_gap_quartile_upper_bound": quartile,
        "spearman_g_min_0_vs_Delta_orient": spearman(gap, orient),
        "spearman_g_min_0_vs_G_open": spearman(gap, opening),
        "lowest_gap_quartile_median_Delta_orient": float(np.median(orient[low])),
        "remaining_points_median_Delta_orient": float(np.median(orient[~low])),
        "lowest_gap_quartile_median_G_open": float(np.median(opening[low])),
        "remaining_points_median_G_open": float(np.median(opening[~low])),
        "causality_claimed": False,
        "regression_law_fitted": False,
    }


def _save_figure(figure: plt.Figure, output_dir: Path, basename: str) -> list[str]:
    pdf = output_dir / f"{basename}.pdf"
    png = output_dir / f"{basename}.png"
    figure.savefig(
        pdf,
        format="pdf",
        metadata={
            "Creator": "CoupledBeams AP-0 screening",
            "CreationDate": None,
        },
    )
    figure.savefig(
        png,
        format="png",
        dpi=PNG_DPI,
        metadata={"Software": "CoupledBeams AP-0 screening"},
    )
    plt.close(figure)
    return [str(pdf.relative_to(ROOT)), str(png.relative_to(ROOT))]


def _heatmap_matrix(
    summaries: Sequence[Mapping[str, Any]], key: str
) -> NDArray[np.float64]:
    mu_values = sorted({float(row["mu"]) for row in summaries})
    tau_values = sorted({float(row["tau"]) for row in summaries})
    lookup = {
        (float(row["mu"]), float(row["tau"])): float(row[key])
        for row in summaries
    }
    return np.asarray(
        [[lookup.get((mu, tau), np.nan) for tau in tau_values] for mu in mu_values],
        dtype=float,
    )


def create_heatmap(
    summaries: Sequence[Mapping[str, Any]],
    *,
    key: str,
    colorbar_label: str,
    threshold: float | None = None,
) -> plt.Figure:
    matrix = _heatmap_matrix(summaries, key)
    mu_values = sorted({float(row["mu"]) for row in summaries})
    tau_values = sorted({float(row["tau"]) for row in summaries})
    with plt.rc_context(supervisor._style_context()):
        figure, axis = plt.subplots(figsize=FIGURE_SIZE_IN, constrained_layout=True)
        upper = max(float(np.nanmax(matrix)), threshold or 0.0, 1.0e-12)
        image = axis.imshow(
            matrix,
            origin="lower",
            aspect="auto",
            cmap="viridis",
            vmin=min(0.0, float(np.min(matrix))),
            vmax=upper,
        )
        axis.set_xticks(
            range(len(tau_values)), [f"{value:g}" for value in tau_values]
        )
        axis.set_yticks(
            range(len(mu_values)), [f"{value:g}" for value in mu_values]
        )
        axis.set_xlabel(r"$\tau$")
        axis.set_ylabel(r"$\mu$")
        for row_index in range(matrix.shape[0]):
            for column_index in range(matrix.shape[1]):
                value = matrix[row_index, column_index]
                if not np.isfinite(value):
                    continue
                axis.text(
                    column_index,
                    row_index,
                    f"{value:.3f}",
                    ha="center",
                    va="center",
                    color="white" if value > 0.55 * upper else "black",
                )
        colorbar = figure.colorbar(image, ax=axis)
        colorbar.set_label(colorbar_label)
        if threshold is not None:
            colorbar.ax.axhline(threshold, color="red", linewidth=1.2)
            axis.text(
                0.99,
                1.015,
                "criterion: 0.10",
                color="red",
                ha="right",
                va="bottom",
                transform=axis.transAxes,
            )
        return figure


def create_gap_scatter(
    point_rows: Sequence[Mapping[str, Any]],
) -> plt.Figure:
    with plt.rc_context(supervisor._style_context()):
        figure, axis = plt.subplots(figsize=FIGURE_SIZE_IN, constrained_layout=True)
        x = np.asarray([float(row["g_min_0"]) for row in point_rows])
        y = np.asarray([float(row["Delta_orient"]) for row in point_rows])
        axis.scatter(x, y, s=22.0, color="#0072B2", alpha=0.75)
        axis.set_xlabel(r"$g_{\min,0}$")
        axis.set_ylabel(r"$\Delta_{\mathrm{orient}}$")
        axis.grid(True, color="#D9D9D9", linewidth=0.5)
        sensitive = sorted(
            point_rows,
            key=lambda row: float(row["Delta_orient"]),
            reverse=True,
        )[:3]
        for index, row in enumerate(sensitive):
            axis.annotate(
                (
                    f"mu={float(row['mu']):g}, tau={float(row['tau']):g}, "
                    f"beta={float(row['beta_deg']):g}"
                ),
                (float(row["g_min_0"]), float(row["Delta_orient"])),
                xytext=(5, -10 - 14 * index),
                textcoords="offset points",
                fontsize=7.0,
            )
        return figure


def create_diagnostic_figures(
    summaries: Sequence[Mapping[str, Any]],
    point_rows: Sequence[Mapping[str, Any]],
    output_dir: Path,
) -> list[str]:
    paths: list[str] = []
    paths.extend(
        _save_figure(
            create_heatmap(
                summaries,
                key="Delta_simpl_max",
                colorbar_label=r"$\Delta_{\mathrm{simpl,max}}$",
                threshold=APPLICABILITY_THRESHOLD,
            ),
            output_dir,
            FIGURE_S1_BASENAME,
        )
    )
    paths.extend(
        _save_figure(
            create_heatmap(
                summaries,
                key="maximum_G_open",
                colorbar_label=r"$\max_\beta G_{\mathrm{open}}$",
            ),
            output_dir,
            FIGURE_S2_BASENAME,
        )
    )
    paths.extend(
        _save_figure(
            create_gap_scatter(point_rows), output_dir, FIGURE_S3_BASENAME
        )
    )
    return paths


def _quality_summary(spectra: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    determinant = np.asarray(
        [float(row["accepted_determinant_residual"]) for row in spectra]
    )
    singular = np.asarray(
        [float(row["accepted_relative_singular_residual"]) for row in spectra]
    )
    return {
        "accepted_root_count": len(spectra),
        "maximum_accepted_determinant_residual": float(np.max(determinant)),
        "maximum_accepted_relative_singular_residual": float(np.max(singular)),
        "quality_basis_counts": {
            basis: sum(row["quality_basis"] == basis for row in spectra)
            for basis in ("scaled", "physical_raw")
        },
        "rejected_root_count": sum(
            row["quality_status"] != "PASS" for row in spectra
        ),
    }


def _aggregate_family_counters(
    families: Mapping[str, Mapping[str, Any]]
) -> dict[str, int | float]:
    counter_names = tuple(PerformanceCounters.__dataclass_fields__)
    result: dict[str, int | float] = {}
    for name in counter_names:
        value = sum(
            float(family.get("performance_counters", {}).get(name, 0.0))
            for family in families.values()
        )
        result[name] = (
            value
            if name in ("family_runtime_s", "total_scientific_runtime_s")
            else int(value)
        )
    hits = int(result.get("transfer_cache_hits", 0))
    misses = int(result.get("transfer_cache_misses", 0))
    result["transfer_cache_requests"] = hits + misses
    result["transfer_cache_hit_rate"] = hits / (hits + misses) if hits + misses else 0.0
    return result


def _source_hashes() -> dict[str, str]:
    return {str(path.relative_to(ROOT)): _sha256(path) for path in SOURCE_PATHS}


def _geometry_rows_as_numbers(
    rows: Sequence[Mapping[str, Any]],
) -> list[GeometryCase]:
    field_names = tuple(GeometryCase.__dataclass_fields__)
    result: list[GeometryCase] = []
    for row in rows:
        values: dict[str, Any] = {"geometry_id": row["geometry_id"]}
        for field in field_names:
            if field != "geometry_id":
                values[field] = float(row[field])
        result.append(GeometryCase(**values))
    return result


def _expected_output_paths(output_dir: Path) -> list[Path]:
    return [
        output_dir / GEOMETRY_FILENAME,
        output_dir / SPECTRA_FILENAME,
        output_dir / POINT_METRICS_FILENAME,
        output_dir / GEOMETRY_SUMMARY_FILENAME,
        output_dir / REUSE_AUDIT_FILENAME,
        output_dir / CANDIDATES_FILENAME,
        output_dir / SUMMARY_FILENAME,
        output_dir / REPORT_FILENAME,
        output_dir / f"{FIGURE_S1_BASENAME}.pdf",
        output_dir / f"{FIGURE_S1_BASENAME}.png",
        output_dir / f"{FIGURE_S2_BASENAME}.pdf",
        output_dir / f"{FIGURE_S2_BASENAME}.png",
        output_dir / f"{FIGURE_S3_BASENAME}.pdf",
        output_dir / f"{FIGURE_S3_BASENAME}.png",
    ]


def determine_workflow_status(
    *,
    smoke: bool,
    geometries: Sequence[GeometryCase],
    beta_values: NDArray[np.float64],
    spectra: Sequence[Mapping[str, Any]],
    point_rows: Sequence[Mapping[str, Any]],
    summaries: Sequence[Mapping[str, Any]],
    reuse_audit: Sequence[Mapping[str, Any]],
    supervisor_hashes_unchanged: bool,
    output_dir: Path,
) -> str:
    expected_geometry_count = 3 if smoke else 9
    expected_spectra = (
        expected_geometry_count
        * len(beta_values)
        * len(MODEL_IDS)
        * GUARD_ROOT_COUNT
    )
    expected_points = expected_geometry_count * len(beta_values)
    outputs_exist = all(path.is_file() for path in _expected_output_paths(output_dir))
    gates = (
        len(geometries) == expected_geometry_count,
        len(beta_values) == (3 if smoke else 19),
        len(spectra) == expected_spectra,
        len(point_rows) == expected_points,
        len(summaries) == expected_geometry_count,
        all(row.get("status") == "PASS" for row in reuse_audit),
        len(reuse_audit) == (3 if smoke else 5),
        all(row["quality_status"] == "PASS" for row in spectra),
        supervisor_hashes_unchanged,
        outputs_exist,
    )
    return "PASS" if all(gates) else "FAIL"


def screening_report(summary: Mapping[str, Any]) -> str:
    geometry_rows = summary["geometry_table"]
    scientific_rows = summary["geometry_results"]
    candidates = summary["candidates"]
    relation = summary["descriptive_gap_relation"]
    lines = [
        "# AP-0 — small-grid spectral-applicability screening at theta=2 deg",
        "",
        f"**AP-0 small-grid screening workflow: {summary['workflow_status']}**",
        "",
        f"Number of geometries within 10%: **{summary['scientific_result']['within_10_percent_count']} / {summary['scientific_result']['geometry_count']}**",
        "",
        f"Number of geometries exceeding 10%: **{summary['scientific_result']['exceeding_10_percent_count']} / {summary['scientific_result']['geometry_count']}**",
        "",
        "Этот screening относится только к Chapter-2 прямоугольным анизотропным стержням HMS/DX-209. T2 является одномерным reference для сопоставления моделей, а не абсолютной точной моделью реальной трёхмерной конструкции.",
        "",
        "Показатели относятся только к finite grid `mu={0,0.25,0.5}`, `tau={-0.2,0,0.2}`, `beta=0,5,...,90 deg`. Точная boundary по theta и значения между beta-grid points не определялись.",
        "",
        "Все сравнения выполнены между independently sorted spectral positions. MAC, mode shapes, energy classification и modal-descendant tracking не применялись.",
        "",
        "## Модели и нормировка",
        "",
        "- `T2`: Chapter-2 `state_corrected` Timoshenko, generalized rectangular torsion, rotated properties и `Sbar16` coupling при `theta1=theta2=2 deg`.",
        "- `T0`: ортотропный предел той же Chapter-2 модели при `theta=0`.",
        "- `EB0`: validated rectangular orthotropic Euler–Bernoulli + Saint-Venant comparator при `theta=0`.",
        "- Для всех моделей используется `book_slope_clamp` и неизменённый `J_book(beta)`.",
        "",
        "```text",
        "Lambda_ref = (rho*A0*omega^2*l^4/(E_x0*I_y0))^(1/4)",
        "a0=0.005 m, b0=0.020 m, l=0.4 m, E_x0=191 GPa, rho=1580 kg/m^3",
        "```",
        "",
        "Физическая T2 model использует фактические rotated properties при 2 deg; reference normalization остаётся theta-zero и geometry-independent.",
        "",
        "## Geometry manifest",
        "",
        "| geometry | mu | tau | L1 | L2 | a1, mm | a2, mm | b1, mm | b2, mm | volume residual | mass, kg |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in geometry_rows:
        lines.append(
            f"| {row['geometry_id']} | {float(row['mu']):g} | {float(row['tau']):g} | {float(row['L1_m']):.6f} | {float(row['L2_m']):.6f} | {1e3*float(row['a1_m']):.6f} | {1e3*float(row['a2_m']):.6f} | {1e3*float(row['b1_m']):.6f} | {1e3*float(row['b2_m']):.6f} | {float(row['volume_relative_residual']):.3e} | {float(row['mass_kg']):.6f} |"
        )
    lines += [
        "",
        f"## Результаты {len(scientific_rows)} геометрий",
        "",
        "| geometry | Delta_beam,max | Delta_orient,max | Delta_simpl,max | beta/mode simpl | min baseline gap | beta/pair gap | max G_open | beta/pair open | classification |",
        "|---|---:|---:|---:|---|---:|---|---:|---|---|",
    ]
    for row in scientific_rows:
        lines.append(
            f"| {row['geometry_id']} | {float(row['Delta_beam_max']):.6f} | {float(row['Delta_orient_max']):.6f} | {float(row['Delta_simpl_max']):.6f} | {float(row['Delta_simpl_max_beta']):g}/{int(row['Delta_simpl_max_mode'])} | {float(row['minimum_baseline_gap']):.6f} | {float(row['minimum_baseline_gap_beta']):g}/{int(row['minimum_baseline_gap_pair'])}-{int(row['minimum_baseline_gap_pair'])+1} | {float(row['maximum_G_open']):.6f} | {float(row['maximum_G_open_beta']):g}/{int(row['maximum_G_open_pair'])}-{int(row['maximum_G_open_pair'])+1} | {row['screening_classification']} |"
        )
    lines += [
        "",
        "Превышение 10% является scientific screening result, а не computational failure.",
        "",
        "## Исходная плотность спектра и раскрытие gaps",
        "",
        f"Spearman `g_min_0` versus `Delta_orient`: `{float(relation['spearman_g_min_0_vs_Delta_orient']):.6f}`; versus `G_open`: `{float(relation['spearman_g_min_0_vs_G_open']):.6f}`.",
        f"Для lowest-gap quartile median `Delta_orient={float(relation['lowest_gap_quartile_median_Delta_orient']):.6f}` и `G_open={float(relation['lowest_gap_quartile_median_G_open']):.6f}`; для остальных points соответственно `{float(relation['remaining_points_median_Delta_orient']):.6f}` и `{float(relation['remaining_points_median_G_open']):.6f}`.",
        "Это только описательная связь на 171 screening points; regression law не подгонялся и причинность не утверждается.",
        "",
        "## Candidates следующего gate",
        "",
        "```json",
        json.dumps(candidates, ensure_ascii=False, indent=2, sort_keys=True),
        "```",
        "",
        "Candidates только записаны. Дополнительные theta, beta refinement, MAC, shapes или energy analysis не запускались.",
        "",
        "## Root quality, reuse и performance",
        "",
        f"Accepted roots: `{summary['root_quality']['accepted_root_count']}`; maximum determinant residual `{summary['root_quality']['maximum_accepted_determinant_residual']:.12e}`; maximum singular residual `{summary['root_quality']['maximum_accepted_relative_singular_residual']:.12e}`; rejected roots: `0`.",
        f"New spectral families: `{summary['families']['new_family_count']}`; reused families: `{summary['families']['reused_family_count']}`; resumed checkpoints: `{summary['families']['resumed_family_count']}`.",
        f"Global anchors: `{summary['performance_counters'].get('global_anchor_scans', 0)}`; inventories: `{summary['performance_counters'].get('global_inventory_checks', 0)}`; fallbacks: `{summary['performance_counters'].get('global_fallback_scans', 0)}`; transfer expm: `{summary['performance_counters'].get('transfer_expm_evaluations', 0)}`.",
        f"Recorded scientific family runtime: `{summary['runtimes_seconds']['recorded_family_total']:.6f} s`; current command scientific runtime: `{summary['runtimes_seconds']['scientific_current']:.6f} s`; current command wall runtime: `{summary['runtimes_seconds']['workflow_wall_current']:.6f} s`.",
        f"Approved reused arrays exactly equal their source beta/frequency/Lambda values: `{summary['reuse_audit_pass']}`. Figure 6 was not used.",
        "",
        "## Ограничения",
        "",
        "- Точная допустимая boundary по theta не искалась.",
        "- Theta refinement и local beta refinement не выполнялись.",
        "- MAC, mode shapes, energy analysis и modal descendants не использовались.",
        "- FEM и 3D FEM не запускались.",
        "- Circular-rod workflows, isotropic steel defaults и старый parameter name не использовались.",
        "- Supervisor Figures 1–12 и их report/manifest не перезаписывались.",
        "",
        "## Воспроизводимость",
        "",
        "```text",
        "python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --smoke",
        "python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --resume",
        "python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --reuse-data",
        "```",
    ]
    return "\n".join(lines) + "\n"


def _safe_screening_output(path: Path) -> Path:
    resolved = supervisor._safe_output_dir(path)
    supervisor_resolved = SUPERVISOR_OUTPUT_DIR.resolve()
    try:
        resolved.relative_to(supervisor_resolved)
    except ValueError:
        return resolved
    raise ValueError("AP-0 output must not be inside the supervisor output directory")


def _json_write(path: Path, value: Mapping[str, Any]) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _build_summary(
    *,
    smoke: bool,
    execution_mode: str,
    geometries: Sequence[GeometryCase],
    beta_values: NDArray[np.float64],
    spectra: Sequence[Mapping[str, Any]],
    point_rows: Sequence[Mapping[str, Any]],
    summaries: Sequence[Mapping[str, Any]],
    reuse_audit: Sequence[Mapping[str, Any]],
    candidates: Mapping[str, Any],
    family_metadata: Mapping[str, Mapping[str, Any]],
    previous_summary: Mapping[str, Any],
    supervisor_hashes_before: Mapping[str, str],
    supervisor_hashes_after: Mapping[str, str],
    scientific_runtime: float,
    wall_runtime: float,
    figure_paths: Sequence[str],
    output_dir: Path,
) -> dict[str, Any]:
    within = sum(
        row["screening_classification"]
        == "WITHIN_10_PERCENT_ON_SCREENING_GRID"
        for row in summaries
    )
    if execution_mode != "reuse_data":
        families_block = {
            "new_family_count": len(family_metadata),
            "reused_family_count": len(reuse_audit),
            "resumed_family_count": sum(
                bool(item.get("resumed_from_checkpoint"))
                for item in family_metadata.values()
            ),
            "metadata": dict(family_metadata),
        }
    else:
        families_block = dict(previous_summary.get("families", {}))
        families_block.setdefault("new_family_count", 0)
        families_block.setdefault("reused_family_count", len(reuse_audit))
        families_block.setdefault("resumed_family_count", 0)
        families_block.setdefault("metadata", {})
    counters = (
        _aggregate_family_counters(family_metadata)
        if family_metadata
        else dict(previous_summary.get("performance_counters", {}))
    )
    quality = _quality_summary(spectra)
    volume_max = max(case.volume_relative_residual for case in geometries)
    return {
        "workflow": "AP-0 small-grid spectral-applicability screening at theta=2 deg",
        "workflow_status": "PENDING",
        "smoke": smoke,
        "smoke_is_scientific_baseline": False,
        "execution_mode": execution_mode,
        "research_line_separation": {
            "chapter_2_rectangular_anisotropic_only": True,
            "circular_isotropic_article_used": False,
            "isotropic_steel_defaults_used": False,
            "fem_used": False,
            "three_dimensional_fem_used": False,
            "complex_roots_used": False,
            "damping_used": False,
            "mode_tracking_used": False,
            "energy_analysis_used": False,
        },
        "git_context": git_context(),
        "models": {
            "T2": "Chapter-2 state_corrected Timoshenko/generalized torsion theta=2 deg",
            "T0": "Chapter-2 state_corrected orthotropic Timoshenko/generalized torsion theta=0",
            "EB0": "rectangular orthotropic Euler-Bernoulli/Saint-Venant theta=0",
            "clamp": CLAMP,
            "joint": "unchanged J_book(beta)",
            "material": "HMS/DX-209 elastic",
        },
        "grids": {
            "mu": [case.mu for case in geometries]
            if smoke
            else list(MU_VALUES),
            "tau": [case.tau for case in geometries]
            if smoke
            else list(TAU_VALUES),
            "geometry_count": len(geometries),
            "beta_deg": beta_values.tolist(),
            "beta_point_count": len(beta_values),
        },
        "geometry_contract": {
            "kappa": KAPPA,
            "L1": "l*(1-mu)",
            "L2": "l*(1+mu)",
            "thickness_parameter_definition": "(a2-a1)/(a1+a2)",
            "scale": "a0/sqrt(1+tau^2+2*mu*tau)",
            "reference_volume_m3": REFERENCE_VOLUME_M3,
            "reference_mass_kg": REFERENCE_MASS_KG,
            "maximum_volume_relative_residual": volume_max,
            "volume_relative_tolerance": VOLUME_RELATIVE_TOLERANCE,
            "hard_slenderness_threshold": None,
        },
        "normalization": normalization_contract(geometries[0], "T2")
        | {"screening_fingerprint_context": "family-specific in checkpoints"},
        "root_contract": {
            "guard_root_count": GUARD_ROOT_COUNT,
            "compared_root_count": COMPARED_ROOT_COUNT,
            "sorted_spectral_positions_not_descendants": True,
            "interpolation": False,
            "accepted_determinant_tolerance": ROOT_DETERMINANT_TOLERANCE,
            "accepted_singular_tolerance": ROOT_SINGULAR_TOLERANCE,
        },
        "metrics": {
            "Delta_beam": "abs(T0-EB0)/T0",
            "Delta_orient": "abs(T2-T0)/T2",
            "Delta_simpl": "abs(T2-EB0)/T2",
            "criterion": APPLICABILITY_THRESHOLD,
            "gap": "2*(Lambda[k+1]-Lambda[k])/(Lambda[k+1]+Lambda[k])",
        },
        "scientific_result": {
            "geometry_count": len(summaries),
            "within_10_percent_count": within,
            "exceeding_10_percent_count": len(summaries) - within,
            "scope": (
                "NOT_A_SCIENTIFIC_BASELINE"
                if smoke
                else "FINITE_AP0_SCREENING_GRID_ONLY"
            ),
        },
        "geometry_table": [asdict(case) for case in geometries],
        "geometry_results": [dict(row) for row in summaries],
        "candidates": dict(candidates),
        "descriptive_gap_relation": descriptive_gap_relation(point_rows),
        "root_quality": quality,
        "families": families_block,
        "performance_counters": counters,
        "reuse_audit_pass": all(
            row.get("status") == "PASS" for row in reuse_audit
        ),
        "reuse_audit": [dict(row) for row in reuse_audit],
        "figure_6_used": False,
        "source_sha256": _source_hashes(),
        "supervisor_figure_data_sha256_before": dict(supervisor_hashes_before),
        "supervisor_figure_data_sha256_after": dict(supervisor_hashes_after),
        "supervisor_figure_data_unchanged": (
            dict(supervisor_hashes_before) == dict(supervisor_hashes_after)
        ),
        "runtimes_seconds": {
            "scientific_current": scientific_runtime,
            "workflow_wall_current": wall_runtime,
            "recorded_family_total": sum(
                float(item.get("recorded_runtime_seconds", 0.0))
                for item in family_metadata.values()
            )
            if family_metadata
            else float(
                previous_summary.get("runtimes_seconds", {}).get(
                    "recorded_family_total", 0.0
                )
            ),
        },
        "row_counts": {
            GEOMETRY_FILENAME: len(geometries),
            SPECTRA_FILENAME: len(spectra),
            POINT_METRICS_FILENAME: len(point_rows),
            GEOMETRY_SUMMARY_FILENAME: len(summaries),
            REUSE_AUDIT_FILENAME: len(reuse_audit),
        },
        "output_directory": str(output_dir.relative_to(ROOT)),
        "diagnostic_figure_paths": list(figure_paths),
        "limitations": {
            "exact_theta_boundary_searched": False,
            "theta_refinement_run": False,
            "beta_refinement_run": False,
            "MAC_run": False,
            "mode_shapes_run": False,
            "energy_analysis_run": False,
            "FEM_run": False,
            "three_dimensional_FEM_run": False,
            "next_gate_started": False,
        },
    }


def preservation_hashes(directory: Path, names: Sequence[str]) -> dict[str, str]:
    return {name: _sha256(directory / name) for name in names}


def supervisor_preservation_hashes() -> dict[str, str]:
    names = [supervisor.DATA_FILENAMES[number] for number in range(1, 13)]
    names.extend(("plot_manifest.json", "report.md"))
    return preservation_hashes(SUPERVISOR_OUTPUT_DIR, names)


def ap1_source_hashes() -> dict[str, str]:
    paths = (*SOURCE_PATHS, AP0_GEOMETRY_PATH, AP0_SPECTRA_PATH,
             AP0_POINT_METRICS_PATH, AP0_GEOMETRY_SUMMARY_PATH,
             AP0_SUMMARY_PATH, FIGURE_07_PATH)
    return {
        str(path.relative_to(ROOT)): _sha256(path)
        for path in paths
        if path.is_file()
    }


def load_ap0_geometries(*, smoke: bool) -> list[GeometryCase]:
    if not AP0_GEOMETRY_PATH.is_file():
        raise FileNotFoundError("AP-1 requires the canonical AP-0 geometry manifest")
    cases = _geometry_rows_as_numbers(_read_csv(AP0_GEOMETRY_PATH))
    validate_geometry_manifest(cases, smoke=False)
    canonical = screening_geometries(smoke=False)
    for source, expected in zip(cases, canonical):
        if source.geometry_id != expected.geometry_id:
            raise RuntimeError("AP-0 geometry IDs differ from the canonical grid")
        for field in GeometryCase.__dataclass_fields__:
            if field == "geometry_id":
                continue
            if not math.isclose(
                float(getattr(source, field)),
                float(getattr(expected, field)),
                rel_tol=1.0e-14,
                abs_tol=2.0e-16,
            ):
                raise RuntimeError(
                    f"AP-0 geometry field mismatch: {source.geometry_id}/{field}"
                )
    return cases[:3] if smoke else cases


def _ap1_reuse_audit_row(
    case: GeometryCase,
    model_id: str,
    *,
    source_path: Path,
    beta_count: int,
    frequency_equal: bool,
    lambda_equal: bool,
    quality_equal: bool,
) -> dict[str, Any]:
    fingerprint = supervisor.stable_fingerprint(case)
    return {
        "geometry_id": case.geometry_id,
        "model_id": model_id,
        "theta_deg": model_theta_deg(model_id),
        "source_file": str(source_path.relative_to(ROOT)),
        "source_sha256": _sha256(source_path),
        "source_geometry_fingerprint": fingerprint,
        "target_geometry_fingerprint": fingerprint,
        "beta_count": beta_count,
        "root_count": GUARD_ROOT_COUNT,
        "frequency_array_equal": frequency_equal,
        "lambda_array_equal": lambda_equal,
        "quality_status_equal": quality_equal,
        "reuse_status": (
            "PASS"
            if frequency_equal and lambda_equal and quality_equal
            else "FAIL"
        ),
    }


def load_ap1_ap0_family(
    case: GeometryCase,
    model_id: str,
    beta_values: NDArray[np.float64],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if model_id not in ("T0", "EB0"):
        raise ValueError("AP-1 AP-0 reuse is limited to T0 and EB0")
    source_rows = spectra_as_numbers(_read_csv(AP0_SPECTRA_PATH))
    beta_set = {float(value) for value in beta_values}
    selected = sorted(
        (
            row
            for row in source_rows
            if row["geometry_id"] == case.geometry_id
            and row["model_id"] == model_id
            and float(row["beta_deg"]) in beta_set
        ),
        key=lambda row: (float(row["beta_deg"]), int(row["mode"])),
    )
    expected = len(beta_values) * GUARD_ROOT_COUNT
    if len(selected) != expected:
        raise RuntimeError(f"{case.geometry_id}/{model_id}: AP-0 reuse is incomplete")
    source_file = str(AP0_SPECTRA_PATH.relative_to(ROOT))
    source_hash = _sha256(AP0_SPECTRA_PATH)
    rows: list[dict[str, Any]] = []
    for source in selected:
        row = dict(source)
        row["data_origin"] = f"reused_ap0_{model_id}"
        row["source_file"] = source_file
        row["source_sha256"] = source_hash
        rows.append(row)
    fields = (
        "beta_deg",
        "mode",
        "frequency_hz",
        "lambda_ref",
        "quality_status",
        "quality_basis",
        "root_status",
        "accepted_determinant_residual",
        "accepted_relative_singular_residual",
    )
    source_array = np.asarray([[row[field] for field in fields] for row in selected], dtype=object)
    target_array = np.asarray([[row[field] for field in fields] for row in rows], dtype=object)
    frequency_equal = np.array_equal(
        np.asarray([row["frequency_hz"] for row in selected]),
        np.asarray([row["frequency_hz"] for row in rows]),
    )
    lambda_equal = np.array_equal(
        np.asarray([row["lambda_ref"] for row in selected]),
        np.asarray([row["lambda_ref"] for row in rows]),
    )
    quality_equal = np.array_equal(source_array, target_array)
    audit = _ap1_reuse_audit_row(
        case,
        model_id,
        source_path=AP0_SPECTRA_PATH,
        beta_count=len(beta_values),
        frequency_equal=frequency_equal,
        lambda_equal=lambda_equal,
        quality_equal=quality_equal,
    )
    if audit["reuse_status"] != "PASS":
        raise RuntimeError(f"{case.geometry_id}/{model_id}: exact AP-0 reuse failed")
    return rows, audit


def load_ap1_figure7_t5_family(
    case: GeometryCase,
    beta_values: NDArray[np.float64],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if (case.mu, case.tau) != (0.0, 0.0):
        raise ValueError("Figure 7 T5 reuse is only valid for mu=0,tau=0")
    source_rows = _read_csv(FIGURE_07_PATH)
    beta_set = {float(value) for value in beta_values}
    selected = sorted(
        (row for row in source_rows if float(row["beta_deg"]) in beta_set),
        key=lambda row: (float(row["beta_deg"]), int(row["mode"])),
    )
    expected = len(beta_values) * GUARD_ROOT_COUNT
    if len(selected) != expected:
        raise RuntimeError("Figure 7 T5 reuse source is incomplete")
    rows: list[dict[str, Any]] = []
    source_frequency: list[float] = []
    source_lambda: list[float] = []
    for source in selected:
        frequency = float(source["left_frequency_hz"])
        lambda_ref = float(source["left_lambda"])
        quality = supervisor._saved_quality(source, "left")
        row = _spectrum_row(
            case,
            beta_deg=float(source["beta_deg"]),
            model_id="T5",
            mode=int(source["mode"]),
            frequency_hz=frequency,
            lambda_ref=lambda_ref,
            quality=quality,
            data_origin="reused_figure_07_T5",
        )
        row["source_file"] = str(FIGURE_07_PATH.relative_to(ROOT))
        row["source_sha256"] = _sha256(FIGURE_07_PATH)
        rows.append(row)
        source_frequency.append(frequency)
        source_lambda.append(lambda_ref)
    frequency_equal = np.array_equal(
        np.asarray(source_frequency),
        np.asarray([row["frequency_hz"] for row in rows]),
    )
    lambda_equal = np.array_equal(
        np.asarray(source_lambda),
        np.asarray([row["lambda_ref"] for row in rows]),
    )
    quality_equal = all(row["quality_status"] == "PASS" for row in rows)
    normalization_residual = max(
        abs(float(fixed_lambda_from_frequency(row["frequency_hz"])) - row["lambda_ref"])
        / max(abs(row["lambda_ref"]), np.finfo(float).tiny)
        for row in rows
    )
    if normalization_residual > 5.0e-15:
        raise RuntimeError("Figure 7 uses an incompatible normalization")
    audit = _ap1_reuse_audit_row(
        case,
        "T5",
        source_path=FIGURE_07_PATH,
        beta_count=len(beta_values),
        frequency_equal=frequency_equal,
        lambda_equal=lambda_equal,
        quality_equal=quality_equal,
    )
    audit["maximum_fixed_normalization_relative_residual"] = normalization_residual
    if audit["reuse_status"] != "PASS":
        raise RuntimeError("Figure 7 T5 exact reuse failed")
    return rows, audit


@dataclass(frozen=True)
class SupervisorBaselineSource:
    theta_deg: float
    figure_number: int
    figure_identifier: str
    source_path: Path


def locate_ap2_supervisor_baseline(theta_deg: float) -> SupervisorBaselineSource:
    """Locate a small-theta source by verified manifest content, not filename."""

    theta = float(theta_deg)
    if theta not in (3.0, 4.0):
        raise ValueError("AP-2 supervisor reuse is limited to theta=3 or 4 deg")
    manifest_path = SUPERVISOR_OUTPUT_DIR / "plot_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    matches: list[SupervisorBaselineSource] = []
    figures = {
        int(item["number"]): item for item in manifest.get("figures", [])
    }
    for name, parameters in manifest.get("parameters", {}).items():
        if not name.startswith("figure_") or not isinstance(parameters, Mapping):
            continue
        if float(parameters.get("theta_1_deg", -1.0)) != theta or float(
            parameters.get("theta_2_deg", -1.0)
        ) != theta:
            continue
        checks = (
            parameters.get("material_name") == "HMS/DX-209",
            parameters.get("material_factory") == "hms_dx_209_material",
            parameters.get("material_mode") == "elastic",
            float(parameters.get("L_1_m", -1.0)) == REFERENCE_LENGTH_M,
            float(parameters.get("L_2_m", -1.0)) == REFERENCE_LENGTH_M,
            float(parameters.get("a_1_m", -1.0)) == REFERENCE_A_M,
            float(parameters.get("a_2_m", -1.0)) == REFERENCE_A_M,
            float(parameters.get("b_1_m", -1.0)) == REFERENCE_B_M,
            float(parameters.get("b_2_m", -1.0)) == REFERENCE_B_M,
            float(parameters.get("mu", -1.0)) == 0.0,
            float(parameters.get("rho_kg_m3", -1.0))
            == MATERIAL_DENSITY_KG_M3,
            float(parameters.get("elastic_Ex_theta0_pa", -1.0))
            == REFERENCE_EX_PA,
        )
        if not all(checks):
            continue
        figure_number = int(name.removeprefix("figure_"))
        figure = figures.get(figure_number)
        if figure is None:
            continue
        source_path = ROOT / str(figure["data_csv"])
        matches.append(
            SupervisorBaselineSource(
                theta_deg=theta,
                figure_number=figure_number,
                figure_identifier=str(figure["basename"]),
                source_path=source_path,
            )
        )
    if len(matches) != 1:
        raise RuntimeError(
            f"theta={theta:g}: expected one identity-matched supervisor baseline, "
            f"found {len(matches)}"
        )
    source = matches[0]
    if not source.source_path.is_file():
        raise FileNotFoundError(source.source_path)
    if supervisor.SMALL_THETA_CLAMP != CLAMP:
        raise RuntimeError("supervisor small-theta clamp is incompatible")
    preset = supervisor.SMALL_THETA_PRESETS.get(source.figure_number)
    if preset is None or float(preset.theta_1_deg) != theta or float(
        preset.theta_2_deg
    ) != theta:
        raise RuntimeError("supervisor preset does not match manifest theta")
    root_contract = manifest.get("root_contract", {})
    if int(root_contract.get("guard_root_count", -1)) != GUARD_ROOT_COUNT:
        raise RuntimeError("supervisor guard-root contract is incompatible")
    lambda_contract = manifest.get("lambda_definition", {}).get(
        "figures_5_12_fixed_reference", {}
    )
    if (
        float(lambda_contract.get("a0_m", -1.0)) != REFERENCE_A_M
        or float(lambda_contract.get("b0_m", -1.0)) != REFERENCE_B_M
        or float(lambda_contract.get("A0_m2", -1.0)) != REFERENCE_AREA_M2
        or float(lambda_contract.get("I_y0_m4", -1.0)) != REFERENCE_IY_M4
    ):
        raise RuntimeError("supervisor fixed Lambda normalization is incompatible")
    return source


def load_ap2_supervisor_family(
    case: GeometryCase,
    model_id: str,
    beta_values: NDArray[np.float64],
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if (case.mu, case.tau) != (0.0, 0.0):
        raise ValueError("AP-2 supervisor baseline reuse is only valid at mu=0,tau=0")
    theta = model_theta_deg(model_id)
    if model_id not in ("T3", "T4"):
        raise ValueError("AP-2 supervisor baseline model must be T3 or T4")
    source = locate_ap2_supervisor_baseline(theta)
    source_rows = _read_csv(source.source_path)
    beta_set = {float(value) for value in beta_values}
    selected = sorted(
        (row for row in source_rows if float(row["beta_deg"]) in beta_set),
        key=lambda row: (float(row["beta_deg"]), int(row["mode"])),
    )
    expected = len(beta_values) * GUARD_ROOT_COUNT
    if len(selected) != expected:
        raise RuntimeError(
            f"{source.figure_identifier}: supervisor source has {len(selected)} "
            f"selected rows, expected {expected}"
        )
    rows: list[dict[str, Any]] = []
    source_frequency: list[float] = []
    source_lambda: list[float] = []
    for source_row in selected:
        expected_model = (
            "Chapter2_monoclinic_Timoshenko_state_corrected_"
            f"generalized_torsion_theta{int(theta)}"
        )
        if (
            int(source_row["figure"]) != source.figure_number
            or source_row["left_model"] != expected_model
            or float(source_row["left_theta_deg"]) != theta
            or float(source_row["right_theta_deg"]) != 0.0
            or source_row["comparison_type"]
            != f"diagnostic_orthotropic_theta0_EB_baseline_for_monoclinic_theta{int(theta)}"
        ):
            raise RuntimeError(
                f"{source.figure_identifier}: scientific model identity mismatch"
            )
        frequency = float(source_row["left_frequency_hz"])
        lambda_ref = float(source_row["left_lambda"])
        quality = supervisor._saved_quality(source_row, "left")
        row = _spectrum_row(
            case,
            beta_deg=float(source_row["beta_deg"]),
            model_id=model_id,
            mode=int(source_row["mode"]),
            frequency_hz=frequency,
            lambda_ref=lambda_ref,
            quality=quality,
            data_origin=f"reused_supervisor_figure_{source.figure_number}_{model_id}",
        )
        row["source_file"] = str(source.source_path.relative_to(ROOT))
        row["source_sha256"] = _sha256(source.source_path)
        rows.append(row)
        source_frequency.append(frequency)
        source_lambda.append(lambda_ref)
    frequency_equal = np.array_equal(
        np.asarray(source_frequency),
        np.asarray([row["frequency_hz"] for row in rows]),
    )
    lambda_equal = np.array_equal(
        np.asarray(source_lambda),
        np.asarray([row["lambda_ref"] for row in rows]),
    )
    quality_equal = all(row["quality_status"] == "PASS" for row in rows)
    normalization_residual = max(
        abs(float(fixed_lambda_from_frequency(row["frequency_hz"])) - row["lambda_ref"])
        / max(abs(row["lambda_ref"]), np.finfo(float).tiny)
        for row in rows
    )
    if normalization_residual > 5.0e-15:
        raise RuntimeError(
            f"{source.figure_identifier}: incompatible fixed Lambda normalization"
        )
    audit = _ap1_reuse_audit_row(
        case,
        model_id,
        source_path=source.source_path,
        beta_count=len(beta_values),
        frequency_equal=frequency_equal,
        lambda_equal=lambda_equal,
        quality_equal=quality_equal,
    )
    audit.update(
        {
            "source_figure": source.figure_number,
            "source_figure_identifier": source.figure_identifier,
            "maximum_fixed_normalization_relative_residual": normalization_residual,
            "scientific_identity_status": "PASS",
        }
    )
    if audit["reuse_status"] != "PASS":
        raise RuntimeError(f"{source.figure_identifier}: exact-array reuse failed")
    return rows, audit


def rejected_ap2_supervisor_reuse_audit(
    case: GeometryCase, model_id: str, reason: str, beta_count: int
) -> dict[str, Any]:
    return {
        "geometry_id": case.geometry_id,
        "model_id": model_id,
        "theta_deg": model_theta_deg(model_id),
        "source_file": "",
        "source_sha256": "",
        "source_geometry_fingerprint": "",
        "target_geometry_fingerprint": supervisor.stable_fingerprint(case),
        "beta_count": beta_count,
        "root_count": GUARD_ROOT_COUNT,
        "frequency_array_equal": False,
        "lambda_array_equal": False,
        "quality_status_equal": False,
        "reuse_status": "SUPERVISOR_BASELINE_REUSE_REJECTED",
        "reuse_rejection_reason": reason,
    }


def attempt_ap2_supervisor_reuse(
    case: GeometryCase,
    model_id: str,
    beta_values: NDArray[np.float64],
) -> tuple[list[dict[str, Any]] | None, dict[str, Any]]:
    try:
        rows, audit = load_ap2_supervisor_family(case, model_id, beta_values)
    except (FileNotFoundError, RuntimeError, ValueError) as error:
        return None, rejected_ap2_supervisor_reuse_audit(
            case, model_id, str(error), len(beta_values)
        )
    return rows, audit


def ap2_point_metrics_from_lambdas(
    target: Sequence[float],
    t0: Sequence[float],
    eb0: Sequence[float],
    theta_deg: float,
) -> dict[str, Any]:
    theta = int(theta_deg)
    if theta not in (3, 4):
        raise ValueError("AP-2 metrics require theta=3 or theta=4")
    base = point_metrics_from_lambdas(target, t0, eb0)
    return {
        "Delta_beam": base["Delta_beam"],
        "Delta_beam_mode": base["Delta_beam_mode"],
        f"Delta_orient_{theta}": base["Delta_orient"],
        f"Delta_orient_{theta}_mode": base["Delta_orient_mode"],
        f"Delta_simpl_{theta}": base["Delta_simpl"],
        f"Delta_simpl_{theta}_mode": base["Delta_simpl_mode"],
        "g_min_0": base["g_min_0"],
        "baseline_min_gap_pair": base["baseline_min_gap_pair"],
        f"g_min_{theta}": base["g_min_2"],
        f"G_nearest_{theta}": base["G_nearest"],
        f"G_open_{theta}": base["G_open"],
        f"G_open_{theta}_pair": base["G_open_pair"],
        f"G_close_{theta}": base["G_close"],
        f"classification_at_point_theta{theta}": base[
            "classification_at_point"
        ],
    }


def build_ap2_metrics(
    spectra: Sequence[Mapping[str, Any]],
    geometries: Sequence[GeometryCase],
    beta_values: NDArray[np.float64],
    theta_deg: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    theta = int(theta_deg)
    model_id = f"T{theta}"
    target_lookup = _spectrum_lookup(spectra)
    ap0_spectra = spectra_as_numbers(_read_csv(AP0_SPECTRA_PATH))
    ap0_lookup = _spectrum_lookup(ap0_spectra)
    ap0_points = point_metrics_as_numbers(_read_csv(AP0_POINT_METRICS_PATH))
    ap0_point_lookup = {
        (str(row["geometry_id"]), float(row["beta_deg"])): row
        for row in ap0_points
    }
    point_rows: list[dict[str, Any]] = []
    pair_rows: list[dict[str, Any]] = []
    for case in geometries:
        for beta_value in beta_values:
            beta = float(beta_value)
            arrays: dict[str, NDArray[np.float64]] = {}
            for current_model in AP2_MODEL_IDS[theta]:
                selected = sorted(
                    target_lookup[(case.geometry_id, beta, current_model)],
                    key=lambda row: int(row["mode"]),
                )
                arrays[current_model] = np.asarray(
                    [float(row["lambda_ref"]) for row in selected], dtype=float
                )
            selected_t2 = sorted(
                ap0_lookup[(case.geometry_id, beta, "T2")],
                key=lambda row: int(row["mode"]),
            )
            t2 = np.asarray(
                [float(row["lambda_ref"]) for row in selected_t2], dtype=float
            )
            metrics = ap2_point_metrics_from_lambdas(
                arrays[model_id], arrays["T0"], arrays["EB0"], theta
            )
            ap0_point = ap0_point_lookup[(case.geometry_id, beta)]
            if (
                metrics["Delta_beam"] != ap0_point["Delta_beam"]
                or metrics["Delta_beam_mode"] != ap0_point["Delta_beam_mode"]
            ):
                raise RuntimeError(
                    f"{case.geometry_id}/beta={beta:g}: Delta_beam changed from AP-0"
                )
            point_rows.append(
                {
                    "geometry_id": case.geometry_id,
                    "mu": case.mu,
                    "tau": case.tau,
                    "beta_deg": beta,
                    **metrics,
                }
            )
            gaps_0 = normalized_neighbor_gaps(arrays["T0"])
            gaps_2 = normalized_neighbor_gaps(t2)
            gaps_target = normalized_neighbor_gaps(arrays[model_id])
            for lower in range(COMPARED_ROOT_COUNT):
                delta_2 = float(gaps_2[lower] - gaps_0[lower])
                delta_target = float(gaps_target[lower] - gaps_0[lower])
                pair_rows.append(
                    {
                        "geometry_id": case.geometry_id,
                        "mu": case.mu,
                        "tau": case.tau,
                        "beta_deg": beta,
                        "pair_lower_mode": lower + 1,
                        "pair_upper_mode": lower + 2,
                        "g_pair_T0": float(gaps_0[lower]),
                        "g_pair_T2": float(gaps_2[lower]),
                        f"g_pair_T{theta}": float(gaps_target[lower]),
                        "delta_g_theta2": delta_2,
                        f"delta_g_theta{theta}": delta_target,
                        "abs_delta_g_theta2": abs(delta_2),
                        f"abs_delta_g_theta{theta}": abs(delta_target),
                        f"change_g_2_to_{theta}": float(
                            gaps_target[lower] - gaps_2[lower]
                        ),
                        "pair_opened_at_theta2": delta_2 > 0.0,
                        f"pair_opened_at_theta{theta}": delta_target > 0.0,
                        f"Delta_orient_{theta}_at_same_point": metrics[
                            f"Delta_orient_{theta}"
                        ],
                        f"Delta_simpl_{theta}_at_same_point": metrics[
                            f"Delta_simpl_{theta}"
                        ],
                    }
                )
    return point_rows, pair_rows


def build_ap2_geometry_summary(
    point_rows: Sequence[Mapping[str, Any]],
    geometries: Sequence[GeometryCase],
    theta_deg: float,
) -> list[dict[str, Any]]:
    theta = int(theta_deg)
    ap0_rows = geometry_summary_as_numbers(_read_csv(AP0_GEOMETRY_SUMMARY_PATH))
    ap0_lookup = {str(row["geometry_id"]): row for row in ap0_rows}
    summaries: list[dict[str, Any]] = []
    for case in geometries:
        selected = [
            row for row in point_rows if row["geometry_id"] == case.geometry_id
        ]
        beam = _extreme_row(selected, "Delta_beam", maximum=True)
        orient = _extreme_row(selected, f"Delta_orient_{theta}", maximum=True)
        simpl = _extreme_row(selected, f"Delta_simpl_{theta}", maximum=True)
        baseline = _extreme_row(selected, "g_min_0", maximum=False)
        nearest = _extreme_row(selected, f"G_nearest_{theta}", maximum=True)
        opening = _extreme_row(selected, f"G_open_{theta}", maximum=True)
        ap0 = ap0_lookup[case.geometry_id]
        if float(beam["Delta_beam"]) != float(ap0["Delta_beam_max"]):
            raise RuntimeError(f"{case.geometry_id}: Delta_beam maximum changed")
        classification = (
            "WITHIN_10_PERCENT_ON_SCREENING_GRID"
            if float(simpl[f"Delta_simpl_{theta}"])
            <= APPLICABILITY_THRESHOLD
            else "EXCEEDS_10_PERCENT_ON_SCREENING_GRID"
        )
        summaries.append(
            {
                "geometry_id": case.geometry_id,
                "mu": case.mu,
                "tau": case.tau,
                "Delta_beam_max": float(beam["Delta_beam"]),
                f"Delta_orient_{theta}_max": float(
                    orient[f"Delta_orient_{theta}"]
                ),
                f"Delta_orient_{theta}_max_beta": float(orient["beta_deg"]),
                f"Delta_orient_{theta}_max_mode": int(
                    orient[f"Delta_orient_{theta}_mode"]
                ),
                f"Delta_simpl_{theta}_max": float(
                    simpl[f"Delta_simpl_{theta}"]
                ),
                f"Delta_simpl_{theta}_max_beta": float(simpl["beta_deg"]),
                f"Delta_simpl_{theta}_max_mode": int(
                    simpl[f"Delta_simpl_{theta}_mode"]
                ),
                "minimum_baseline_gap": float(baseline["g_min_0"]),
                "minimum_baseline_gap_beta": float(baseline["beta_deg"]),
                "minimum_baseline_gap_pair": int(
                    baseline["baseline_min_gap_pair"]
                ),
                f"maximum_G_nearest_{theta}": float(
                    nearest[f"G_nearest_{theta}"]
                ),
                f"maximum_G_nearest_{theta}_beta": float(nearest["beta_deg"]),
                f"maximum_G_open_{theta}": float(opening[f"G_open_{theta}"]),
                f"maximum_G_open_{theta}_beta": float(opening["beta_deg"]),
                f"maximum_G_open_{theta}_pair": int(
                    opening[f"G_open_{theta}_pair"]
                ),
                f"screening_classification_theta{theta}": classification,
                "Delta_simpl_2_max": float(ap0["Delta_simpl_max"]),
                f"Delta_simpl_increment_{theta}_minus_2": float(
                    simpl[f"Delta_simpl_{theta}"]
                )
                - float(ap0["Delta_simpl_max"]),
                "classification_theta2": str(ap0["screening_classification"]),
            }
        )
    return summaries


def select_ap2_candidates(
    summaries: Sequence[Mapping[str, Any]], theta_deg: float
) -> dict[str, Any]:
    theta = int(theta_deg)
    key = f"Delta_simpl_{theta}_max"
    stable = _extreme_row(summaries, key, maximum=False)
    sensitive = _extreme_row(summaries, key, maximum=True)
    nearest = min(
        summaries,
        key=lambda row: abs(float(row[key]) - APPLICABILITY_THRESHOLD),
    )

    def compact(row: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "geometry_id": row["geometry_id"],
            "mu": float(row["mu"]),
            "tau": float(row["tau"]),
            key: float(row[key]),
        }

    distance = abs(float(nearest[key]) - APPLICABILITY_THRESHOLD)
    return {
        f"stable_theta{theta}_candidate": compact(stable),
        f"sensitive_theta{theta}_candidate": compact(sensitive),
        "nearest_to_threshold_candidate": compact(nearest),
        "nearest_to_threshold_distance": distance,
        "borderline_status": (
            "CLOSE_BORDERLINE_CASE_FOUND"
            if distance <= 0.03
            else "NO_CLOSE_BORDERLINE_CASE_FOUND"
        ),
        "next_gate_started": False,
    }


def ap2_pairwise_association(
    pair_rows: Sequence[Mapping[str, Any]], theta_deg: float
) -> dict[str, Any]:
    theta = int(theta_deg)
    baseline = np.asarray([float(row["g_pair_T0"]) for row in pair_rows])
    change = np.asarray(
        [float(row[f"abs_delta_g_theta{theta}"]) for row in pair_rows]
    )
    opened = np.asarray(
        [bool(row[f"pair_opened_at_theta{theta}"]) for row in pair_rows]
    )

    def spearman(left: NDArray[np.float64], right: NDArray[np.float64]) -> float:
        return float(np.corrcoef(_rankdata(left), _rankdata(right))[0, 1])

    quartile = float(np.quantile(baseline, 0.25))
    low = baseline <= quartile
    per_pair: dict[str, float] = {}
    for lower in range(1, 7):
        selected = np.asarray(
            [int(row["pair_lower_mode"]) == lower for row in pair_rows]
        )
        per_pair[f"{lower}-{lower + 1}"] = spearman(
            baseline[selected], change[selected]
        )
    return {
        "observation_count": len(pair_rows),
        f"pooled_spearman_g_pair_T0_vs_abs_delta_g_theta{theta}": spearman(
            baseline, change
        ),
        "per_pair_spearman": per_pair,
        "lowest_gap_quartile_upper_bound": quartile,
        "lowest_gap_quartile_count": int(np.count_nonzero(low)),
        f"lowest_gap_quartile_median_abs_delta_g_theta{theta}": float(
            np.median(change[low])
        ),
        f"remaining_median_abs_delta_g_theta{theta}": float(
            np.median(change[~low])
        ),
        "lowest_gap_quartile_opened_pair_fraction": float(
            np.mean(opened[low])
        ),
        "remaining_opened_pair_fraction": float(np.mean(opened[~low])),
        "inferential_p_values_computed": False,
        "independence_claimed": False,
        "regression_law_fitted": False,
        "causality_claimed": False,
    }


def ap1_point_metrics_from_lambdas(
    t5: Sequence[float], t0: Sequence[float], eb0: Sequence[float]
) -> dict[str, Any]:
    base = point_metrics_from_lambdas(t5, t0, eb0)
    return {
        "Delta_beam": base["Delta_beam"],
        "Delta_beam_mode": base["Delta_beam_mode"],
        "Delta_orient_5": base["Delta_orient"],
        "Delta_orient_5_mode": base["Delta_orient_mode"],
        "Delta_simpl_5": base["Delta_simpl"],
        "Delta_simpl_5_mode": base["Delta_simpl_mode"],
        "g_min_0": base["g_min_0"],
        "baseline_min_gap_pair": base["baseline_min_gap_pair"],
        "g_min_5": base["g_min_2"],
        "G_nearest_5": base["G_nearest"],
        "G_open_5": base["G_open"],
        "G_open_5_pair": base["G_open_pair"],
        "G_close_5": base["G_close"],
        "classification_at_point_theta5": base["classification_at_point"],
    }


def _spectrum_lookup(
    rows: Sequence[Mapping[str, Any]],
) -> dict[tuple[str, float, str], list[Mapping[str, Any]]]:
    lookup: dict[tuple[str, float, str], list[Mapping[str, Any]]] = {}
    for row in rows:
        key = (
            str(row["geometry_id"]),
            float(row["beta_deg"]),
            str(row["model_id"]),
        )
        lookup.setdefault(key, []).append(row)
    return lookup


def build_ap1_metrics(
    spectra: Sequence[Mapping[str, Any]],
    geometries: Sequence[GeometryCase],
    beta_values: NDArray[np.float64],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    ap1_lookup = _spectrum_lookup(spectra)
    ap0_spectra = spectra_as_numbers(_read_csv(AP0_SPECTRA_PATH))
    ap0_lookup = _spectrum_lookup(ap0_spectra)
    ap0_points = point_metrics_as_numbers(_read_csv(AP0_POINT_METRICS_PATH))
    ap0_point_lookup = {
        (str(row["geometry_id"]), float(row["beta_deg"])): row
        for row in ap0_points
    }
    point_rows: list[dict[str, Any]] = []
    pair_rows: list[dict[str, Any]] = []
    comparison_rows: list[dict[str, Any]] = []
    for case in geometries:
        for beta_value in beta_values:
            beta_deg = float(beta_value)
            arrays: dict[str, NDArray[np.float64]] = {}
            for model_id in AP1_MODEL_IDS:
                selected = sorted(
                    ap1_lookup[(case.geometry_id, beta_deg, model_id)],
                    key=lambda row: int(row["mode"]),
                )
                arrays[model_id] = np.asarray(
                    [float(row["lambda_ref"]) for row in selected], dtype=float
                )
            t2_selected = sorted(
                ap0_lookup[(case.geometry_id, beta_deg, "T2")],
                key=lambda row: int(row["mode"]),
            )
            t2 = np.asarray(
                [float(row["lambda_ref"]) for row in t2_selected], dtype=float
            )
            metrics = ap1_point_metrics_from_lambdas(
                arrays["T5"], arrays["T0"], arrays["EB0"]
            )
            ap0_point = ap0_point_lookup[(case.geometry_id, beta_deg)]
            if (
                metrics["Delta_beam"] != ap0_point["Delta_beam"]
                or metrics["Delta_beam_mode"] != ap0_point["Delta_beam_mode"]
            ):
                raise RuntimeError(
                    f"{case.geometry_id}/beta={beta_deg:g}: Delta_beam changed from AP-0"
                )
            point_row = {
                "geometry_id": case.geometry_id,
                "mu": case.mu,
                "tau": case.tau,
                "beta_deg": beta_deg,
                **metrics,
            }
            point_rows.append(point_row)
            comparison_rows.append(
                {
                    "geometry_id": case.geometry_id,
                    "mu": case.mu,
                    "tau": case.tau,
                    "beta_deg": beta_deg,
                    "Delta_beam_theta2": ap0_point["Delta_beam"],
                    "Delta_beam_theta5": metrics["Delta_beam"],
                    "Delta_beam_array_equal": True,
                    "Delta_orient_2": ap0_point["Delta_orient"],
                    "Delta_orient_5": metrics["Delta_orient_5"],
                    "Delta_orient_increment_5_minus_2": (
                        metrics["Delta_orient_5"] - ap0_point["Delta_orient"]
                    ),
                    "Delta_simpl_2": ap0_point["Delta_simpl"],
                    "Delta_simpl_5": metrics["Delta_simpl_5"],
                    "Delta_simpl_increment_5_minus_2": (
                        metrics["Delta_simpl_5"] - ap0_point["Delta_simpl"]
                    ),
                    "classification_at_point_theta2": ap0_point[
                        "classification_at_point"
                    ],
                    "classification_at_point_theta5": metrics[
                        "classification_at_point_theta5"
                    ],
                }
            )
            gaps_0 = normalized_neighbor_gaps(arrays["T0"])
            gaps_2 = normalized_neighbor_gaps(t2)
            gaps_5 = normalized_neighbor_gaps(arrays["T5"])
            for lower in range(COMPARED_ROOT_COUNT):
                delta_2 = float(gaps_2[lower] - gaps_0[lower])
                delta_5 = float(gaps_5[lower] - gaps_0[lower])
                pair_rows.append(
                    {
                        "geometry_id": case.geometry_id,
                        "mu": case.mu,
                        "tau": case.tau,
                        "beta_deg": beta_deg,
                        "pair_lower_mode": lower + 1,
                        "pair_upper_mode": lower + 2,
                        "g_pair_T0": float(gaps_0[lower]),
                        "g_pair_T2": float(gaps_2[lower]),
                        "g_pair_T5": float(gaps_5[lower]),
                        "delta_g_theta2": delta_2,
                        "delta_g_theta5": delta_5,
                        "abs_delta_g_theta2": abs(delta_2),
                        "abs_delta_g_theta5": abs(delta_5),
                        "additional_change_theta2_to_theta5": float(
                            gaps_5[lower] - gaps_2[lower]
                        ),
                        "pair_opened_at_theta2": delta_2 > 0.0,
                        "pair_opened_at_theta5": delta_5 > 0.0,
                        "Delta_orient_5_at_same_point": metrics[
                            "Delta_orient_5"
                        ],
                        "Delta_simpl_5_at_same_point": metrics[
                            "Delta_simpl_5"
                        ],
                    }
                )
    return point_rows, pair_rows, comparison_rows


def _classification_transition(theta2: str, theta5: str) -> str:
    if (
        theta2 == "WITHIN_10_PERCENT_ON_SCREENING_GRID"
        and theta5 == "WITHIN_10_PERCENT_ON_SCREENING_GRID"
    ):
        return "REMAINS_WITHIN_AT_2_AND_5"
    if (
        theta2 == "WITHIN_10_PERCENT_ON_SCREENING_GRID"
        and theta5 == "EXCEEDS_10_PERCENT_ON_SCREENING_GRID"
    ):
        return "WITHIN_AT_2_EXCEEDS_AT_5"
    return f"{theta2}_TO_{theta5}"


def build_ap1_geometry_summary(
    point_rows: Sequence[Mapping[str, Any]],
    geometries: Sequence[GeometryCase],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    ap0_rows = geometry_summary_as_numbers(
        _read_csv(AP0_GEOMETRY_SUMMARY_PATH)
    )
    ap0_lookup = {str(row["geometry_id"]): row for row in ap0_rows}
    summaries: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []
    for case in geometries:
        selected = [
            row for row in point_rows if row["geometry_id"] == case.geometry_id
        ]
        beam = _extreme_row(selected, "Delta_beam", maximum=True)
        orient = _extreme_row(selected, "Delta_orient_5", maximum=True)
        simpl = _extreme_row(selected, "Delta_simpl_5", maximum=True)
        baseline = _extreme_row(selected, "g_min_0", maximum=False)
        nearest = _extreme_row(selected, "G_nearest_5", maximum=True)
        opening = _extreme_row(selected, "G_open_5", maximum=True)
        ap0 = ap0_lookup[case.geometry_id]
        if float(beam["Delta_beam"]) != float(ap0["Delta_beam_max"]):
            raise RuntimeError(f"{case.geometry_id}: Delta_beam maximum changed")
        classification_5 = (
            "WITHIN_10_PERCENT_ON_SCREENING_GRID"
            if float(simpl["Delta_simpl_5"]) <= APPLICABILITY_THRESHOLD
            else "EXCEEDS_10_PERCENT_ON_SCREENING_GRID"
        )
        classification_2 = str(ap0["screening_classification"])
        transition = _classification_transition(classification_2, classification_5)
        summary = {
            "geometry_id": case.geometry_id,
            "mu": case.mu,
            "tau": case.tau,
            "Delta_beam_max": float(beam["Delta_beam"]),
            "Delta_orient_5_max": float(orient["Delta_orient_5"]),
            "Delta_orient_5_max_beta": float(orient["beta_deg"]),
            "Delta_orient_5_max_mode": int(orient["Delta_orient_5_mode"]),
            "Delta_simpl_5_max": float(simpl["Delta_simpl_5"]),
            "Delta_simpl_5_max_beta": float(simpl["beta_deg"]),
            "Delta_simpl_5_max_mode": int(simpl["Delta_simpl_5_mode"]),
            "minimum_baseline_gap": float(baseline["g_min_0"]),
            "minimum_baseline_gap_beta": float(baseline["beta_deg"]),
            "minimum_baseline_gap_pair": int(
                baseline["baseline_min_gap_pair"]
            ),
            "maximum_G_nearest_5": float(nearest["G_nearest_5"]),
            "maximum_G_nearest_5_beta": float(nearest["beta_deg"]),
            "maximum_G_open_5": float(opening["G_open_5"]),
            "maximum_G_open_5_beta": float(opening["beta_deg"]),
            "maximum_G_open_5_pair": int(opening["G_open_5_pair"]),
            "screening_classification_theta5": classification_5,
            "Delta_simpl_2_max": float(ap0["Delta_simpl_max"]),
            "Delta_simpl_increment_5_minus_2": float(
                simpl["Delta_simpl_5"]
            )
            - float(ap0["Delta_simpl_max"]),
            "classification_theta2": classification_2,
            "classification_transition": transition,
        }
        summaries.append(summary)
        comparisons.append(
            {
                "geometry_id": case.geometry_id,
                "mu": case.mu,
                "tau": case.tau,
                "Delta_orient_2_max": float(ap0["Delta_orient_max"]),
                "Delta_orient_2_max_beta": float(ap0["Delta_orient_max_beta"]),
                "Delta_orient_2_max_mode": int(ap0["Delta_orient_max_mode"]),
                "Delta_orient_5_max": summary["Delta_orient_5_max"],
                "Delta_orient_5_max_beta": summary["Delta_orient_5_max_beta"],
                "Delta_orient_5_max_mode": summary["Delta_orient_5_max_mode"],
                "Delta_orient_increment_5_minus_2": (
                    summary["Delta_orient_5_max"]
                    - float(ap0["Delta_orient_max"])
                ),
                "Delta_simpl_2_max": float(ap0["Delta_simpl_max"]),
                "Delta_simpl_2_max_beta": float(ap0["Delta_simpl_max_beta"]),
                "Delta_simpl_2_max_mode": int(ap0["Delta_simpl_max_mode"]),
                "Delta_simpl_5_max": summary["Delta_simpl_5_max"],
                "Delta_simpl_5_max_beta": summary["Delta_simpl_5_max_beta"],
                "Delta_simpl_5_max_mode": summary["Delta_simpl_5_max_mode"],
                "Delta_simpl_increment_5_minus_2": summary[
                    "Delta_simpl_increment_5_minus_2"
                ],
                "classification_theta2": classification_2,
                "classification_theta5": classification_5,
                "classification_transition": transition,
            }
        )
    return summaries, comparisons


def ap1_pairwise_association(
    pair_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    baseline = np.asarray([float(row["g_pair_T0"]) for row in pair_rows])
    change = np.asarray(
        [float(row["abs_delta_g_theta5"]) for row in pair_rows]
    )

    def spearman(left: NDArray[np.float64], right: NDArray[np.float64]) -> float:
        return float(np.corrcoef(_rankdata(left), _rankdata(right))[0, 1])

    quartile = float(np.quantile(baseline, 0.25))
    low = baseline <= quartile
    opened = np.asarray(
        [bool(row["pair_opened_at_theta5"]) for row in pair_rows], dtype=bool
    )
    per_pair: dict[str, float] = {}
    for lower in range(1, 7):
        indices = np.asarray(
            [int(row["pair_lower_mode"]) == lower for row in pair_rows]
        )
        per_pair[f"{lower}-{lower + 1}"] = spearman(
            baseline[indices], change[indices]
        )
    return {
        "observation_count": len(pair_rows),
        "pooled_spearman_g_pair_T0_vs_abs_delta_g_theta5": spearman(
            baseline, change
        ),
        "per_pair_spearman": per_pair,
        "lowest_gap_quartile_upper_bound": quartile,
        "lowest_gap_quartile_count": int(np.count_nonzero(low)),
        "lowest_gap_quartile_median_abs_delta_g_theta5": float(
            np.median(change[low])
        ),
        "remaining_median_abs_delta_g_theta5": float(np.median(change[~low])),
        "lowest_gap_quartile_opened_pair_fraction": float(np.mean(opened[low])),
        "remaining_opened_pair_fraction": float(np.mean(opened[~low])),
        "inferential_p_values_computed": False,
        "independence_claimed": False,
        "regression_law_fitted": False,
        "causality_claimed": False,
    }


def select_ap1_candidates(
    summaries: Sequence[Mapping[str, Any]],
    pair_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    stable = min(summaries, key=lambda row: float(row["Delta_simpl_5_max"]))
    sensitive = max(summaries, key=lambda row: float(row["Delta_simpl_5_max"]))
    nearest = min(
        summaries,
        key=lambda row: abs(float(row["Delta_simpl_5_max"]) - 0.10),
    )
    nearest_distance = abs(float(nearest["Delta_simpl_5_max"]) - 0.10)
    smallest_gap = min(pair_rows, key=lambda row: float(row["g_pair_T0"]))
    largest_change = max(
        pair_rows, key=lambda row: float(row["abs_delta_g_theta5"])
    )

    def geometry(row: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "geometry_id": row["geometry_id"],
            "mu": float(row["mu"]),
            "tau": float(row["tau"]),
            "Delta_simpl_5_max": float(row["Delta_simpl_5_max"]),
        }

    def pair(row: Mapping[str, Any]) -> dict[str, Any]:
        return {
            "geometry_id": row["geometry_id"],
            "mu": float(row["mu"]),
            "tau": float(row["tau"]),
            "beta_deg": float(row["beta_deg"]),
            "pair": f"{int(row['pair_lower_mode'])}-{int(row['pair_upper_mode'])}",
            "g_pair_T0": float(row["g_pair_T0"]),
            "abs_delta_g_theta5": float(row["abs_delta_g_theta5"]),
        }

    return {
        "stable_theta5_candidate": geometry(stable),
        "sensitive_theta5_candidate": geometry(sensitive),
        "nearest_to_threshold_candidate": (
            geometry(nearest)
            if nearest_distance <= 0.03
            else "NO_CLOSE_BORDERLINE_CASE_FOUND"
        ),
        "nearest_to_threshold_distance": nearest_distance,
        "transition_candidates": [
            geometry(row)
            for row in summaries
            if row["classification_transition"]
            == "WITHIN_AT_2_EXCEEDS_AT_5"
        ],
        "smallest_pairwise_gap_candidate": pair(smallest_gap),
        "largest_pairwise_change_candidate": pair(largest_change),
        "next_stage_started": False,
    }


def _coerce_csv_rows(
    rows: Sequence[Mapping[str, Any]],
    *,
    float_fields: Sequence[str] = (),
    int_fields: Sequence[str] = (),
    bool_fields: Sequence[str] = (),
) -> list[dict[str, Any]]:
    converted: list[dict[str, Any]] = []
    for source in rows:
        row = dict(source)
        for field in float_fields:
            row[field] = float(row[field])
        for field in int_fields:
            row[field] = int(row[field])
        for field in bool_fields:
            value = row[field]
            row[field] = (
                value
                if isinstance(value, bool)
                else str(value).strip().lower() == "true"
            )
        converted.append(row)
    return converted


def create_ap1_figures(
    summaries: Sequence[Mapping[str, Any]],
    comparisons: Sequence[Mapping[str, Any]],
    pair_rows: Sequence[Mapping[str, Any]],
    output_dir: Path,
) -> list[str]:
    paths: list[str] = []
    paths.extend(
        _save_figure(
            create_heatmap(
                summaries,
                key="Delta_simpl_5_max",
                colorbar_label=r"$\Delta_{\mathrm{simpl},5,\max}$",
                threshold=APPLICABILITY_THRESHOLD,
            ),
            output_dir,
            AP1_FIGURE_S1_BASENAME,
        )
    )
    with plt.rc_context(supervisor._style_context()):
        figure, axis = plt.subplots(
            figsize=(7.2, 4.8), constrained_layout=True
        )
        x = np.arange(len(comparisons), dtype=float)
        theta2 = np.asarray(
            [float(row["Delta_simpl_2_max"]) for row in comparisons]
        )
        theta5 = np.asarray(
            [float(row["Delta_simpl_5_max"]) for row in comparisons]
        )
        for index in range(len(comparisons)):
            axis.plot(
                [x[index], x[index]],
                [theta2[index], theta5[index]],
                color="#777777",
                linewidth=0.9,
                zorder=1,
            )
        axis.scatter(x, theta2, color="#0072B2", s=32, label="theta=2 deg")
        axis.scatter(x, theta5, color="#D55E00", s=32, label="theta=5 deg")
        axis.axhline(0.10, color="#CC0000", linewidth=1.1, linestyle="--")
        axis.set_xticks(
            x,
            [
                f"{float(row['mu']):g}/{float(row['tau']):g}"
                for row in comparisons
            ],
            rotation=45,
            ha="right",
        )
        axis.set_xlabel(r"$\mu/\tau$")
        axis.set_ylabel(r"$\Delta_{\mathrm{simpl},\max}$")
        axis.grid(True, axis="y", color="#D9D9D9", linewidth=0.5)
        axis.legend(frameon=False, loc="upper left")
        paths.extend(_save_figure(figure, output_dir, AP1_FIGURE_S2_BASENAME))
    with plt.rc_context(supervisor._style_context()):
        figure, axis = plt.subplots(
            figsize=(7.2, 4.8), constrained_layout=True
        )
        colors = supervisor.MODE_COLORS
        for lower in range(1, 7):
            selected = [
                row
                for row in pair_rows
                if int(row["pair_lower_mode"]) == lower
            ]
            axis.scatter(
                [float(row["g_pair_T0"]) for row in selected],
                [float(row["abs_delta_g_theta5"]) for row in selected],
                s=17,
                alpha=0.65,
                color=colors[lower - 1],
                label=f"{lower}-{lower + 1}",
            )
        extremes = sorted(
            pair_rows,
            key=lambda row: float(row["abs_delta_g_theta5"]),
            reverse=True,
        )[:3]
        for index, row in enumerate(extremes):
            axis.annotate(
                (
                    f"mu={float(row['mu']):g}, tau={float(row['tau']):g}, "
                    f"beta={float(row['beta_deg']):g}, "
                    f"{int(row['pair_lower_mode'])}-{int(row['pair_upper_mode'])}"
                ),
                (
                    float(row["g_pair_T0"]),
                    float(row["abs_delta_g_theta5"]),
                ),
                xytext=(5, -8 - 13 * index),
                textcoords="offset points",
                fontsize=7,
            )
        axis.set_xlabel(r"$g_{k,0}$")
        axis.set_ylabel(r"$|g_{k,5}-g_{k,0}|$")
        axis.grid(True, color="#D9D9D9", linewidth=0.5)
        axis.legend(
            title="sorted pair",
            frameon=False,
            ncol=2,
            fontsize=8,
            title_fontsize=8,
        )
        paths.extend(_save_figure(figure, output_dir, AP1_FIGURE_S3_BASENAME))
    return paths


def _ap1_expected_output_paths(output_dir: Path) -> list[Path]:
    names = (
        GEOMETRY_FILENAME,
        SPECTRA_FILENAME,
        POINT_METRICS_FILENAME,
        GEOMETRY_SUMMARY_FILENAME,
        PAIRWISE_FILENAME,
        POINT_COMPARISON_FILENAME,
        GEOMETRY_COMPARISON_FILENAME,
        REUSE_AUDIT_FILENAME,
        CANDIDATES_FILENAME,
        SUMMARY_FILENAME,
        REPORT_FILENAME,
        f"{AP1_FIGURE_S1_BASENAME}.pdf",
        f"{AP1_FIGURE_S1_BASENAME}.png",
        f"{AP1_FIGURE_S2_BASENAME}.pdf",
        f"{AP1_FIGURE_S2_BASENAME}.png",
        f"{AP1_FIGURE_S3_BASENAME}.pdf",
        f"{AP1_FIGURE_S3_BASENAME}.png",
    )
    return [output_dir / name for name in names]


def determine_ap1_status(
    *,
    smoke: bool,
    geometries: Sequence[GeometryCase],
    beta_values: NDArray[np.float64],
    spectra: Sequence[Mapping[str, Any]],
    point_rows: Sequence[Mapping[str, Any]],
    pair_rows: Sequence[Mapping[str, Any]],
    point_comparisons: Sequence[Mapping[str, Any]],
    summaries: Sequence[Mapping[str, Any]],
    geometry_comparisons: Sequence[Mapping[str, Any]],
    reuse_audit: Sequence[Mapping[str, Any]],
    ap0_unchanged: bool,
    supervisor_unchanged: bool,
    output_dir: Path,
) -> str:
    geometry_count = 3 if smoke else 9
    beta_count = 3 if smoke else 19
    reuse_count = 7 if smoke else 19
    if not FIGURE_07_PATH.is_file():
        reuse_count -= 1
    gates = (
        len(geometries) == geometry_count,
        len(beta_values) == beta_count,
        len(spectra) == geometry_count * beta_count * 3 * 7,
        len(point_rows) == geometry_count * beta_count,
        len(pair_rows) == geometry_count * beta_count * 6,
        len(point_comparisons) == geometry_count * beta_count,
        len(summaries) == geometry_count,
        len(geometry_comparisons) == geometry_count,
        len(reuse_audit) == reuse_count,
        all(row["reuse_status"] == "PASS" for row in reuse_audit),
        all(row["quality_status"] == "PASS" for row in spectra),
        all(bool(row["Delta_beam_array_equal"]) for row in point_comparisons),
        ap0_unchanged,
        supervisor_unchanged,
        all(path.is_file() for path in _ap1_expected_output_paths(output_dir)),
    )
    return "PASS" if all(gates) else "FAIL"


def build_ap1_summary(
    *,
    smoke: bool,
    execution_mode: str,
    geometries: Sequence[GeometryCase],
    beta_values: NDArray[np.float64],
    spectra: Sequence[Mapping[str, Any]],
    point_rows: Sequence[Mapping[str, Any]],
    pair_rows: Sequence[Mapping[str, Any]],
    point_comparisons: Sequence[Mapping[str, Any]],
    summaries: Sequence[Mapping[str, Any]],
    geometry_comparisons: Sequence[Mapping[str, Any]],
    reuse_audit: Sequence[Mapping[str, Any]],
    candidates: Mapping[str, Any],
    association: Mapping[str, Any],
    family_metadata: Mapping[str, Mapping[str, Any]],
    previous_summary: Mapping[str, Any],
    ap0_before: Mapping[str, str],
    ap0_after: Mapping[str, str],
    supervisor_before: Mapping[str, str],
    supervisor_after: Mapping[str, str],
    scientific_runtime: float,
    wall_runtime: float,
    figure_paths: Sequence[str],
    output_dir: Path,
) -> dict[str, Any]:
    if execution_mode == "reuse_data":
        families = dict(previous_summary.get("families", {}))
        counters = dict(previous_summary.get("performance_counters", {}))
        recorded_runtime = float(
            previous_summary.get("runtimes_seconds", {}).get(
                "recorded_family_total", 0.0
            )
        )
    else:
        families = {
            "new_family_count": len(family_metadata),
            "reused_family_count": len(reuse_audit),
            "resumed_family_count": sum(
                bool(item.get("resumed_from_checkpoint"))
                for item in family_metadata.values()
            ),
            "metadata": dict(family_metadata),
        }
        counters = _aggregate_family_counters(family_metadata)
        recorded_runtime = sum(
            float(item.get("recorded_runtime_seconds", 0.0))
            for item in family_metadata.values()
        )
    families.setdefault("new_family_count", 0)
    families.setdefault("reused_family_count", len(reuse_audit))
    families.setdefault("resumed_family_count", 0)
    families.setdefault("metadata", {})
    within = sum(
        row["screening_classification_theta5"]
        == "WITHIN_10_PERCENT_ON_SCREENING_GRID"
        for row in summaries
    )
    transitions = sum(
        row["classification_transition"] == "WITHIN_AT_2_EXCEEDS_AT_5"
        for row in summaries
    )
    return {
        "workflow": "AP-1 same-grid spectral-applicability screening at theta=5 deg",
        "workflow_status": "PENDING",
        "execution_mode": execution_mode,
        "smoke": smoke,
        "smoke_is_scientific_baseline": False,
        "git_context": git_context(),
        "configuration": {
            "target_theta_deg": 5.0,
            "models": list(AP1_MODEL_IDS),
            "material": "HMS/DX-209 elastic",
            "clamp": CLAMP,
            "joint": "unchanged J_book(beta)",
            "geometry_source": str(AP0_GEOMETRY_PATH.relative_to(ROOT)),
            "geometry_source_sha256": _sha256(AP0_GEOMETRY_PATH),
            "geometry_manifest_content_identical_to_ap0": (
                smoke
                or _sha256(output_dir / GEOMETRY_FILENAME)
                == _sha256(AP0_GEOMETRY_PATH)
            ),
            "fixed_lambda_normalization": normalization_contract(
                geometries[0], "T5"
            ),
        },
        "grid": {
            "mu": list(MU_VALUES),
            "tau": list(TAU_VALUES),
            "geometry_count": len(geometries),
            "beta_deg": beta_values.tolist(),
            "beta_count": len(beta_values),
            "root_count": GUARD_ROOT_COUNT,
            "compared_root_count": COMPARED_ROOT_COUNT,
        },
        "scientific_result": {
            "within_10_percent_at_theta5": within,
            "exceeding_10_percent_at_theta5": len(summaries) - within,
            "classification_transitions_within2_exceeds5": transitions,
            "geometry_count": len(summaries),
            "scope": "NOT_A_SCIENTIFIC_BASELINE" if smoke else "FINITE_AP1_GRID_ONLY",
        },
        "geometry_table": [asdict(case) for case in geometries],
        "geometry_results": [dict(row) for row in summaries],
        "theta2_theta5_geometry_comparison": [
            dict(row) for row in geometry_comparisons
        ],
        "pairwise_descriptive_association": dict(association),
        "candidates": dict(candidates),
        "root_quality": _quality_summary(spectra),
        "families": families,
        "expected_family_accounting": {
            "new_T5": 2 if smoke else 8,
            "reused_T5": 1,
            "reused_T0": len(geometries),
            "reused_EB0": len(geometries),
            "new_root_count": (2 if smoke else 8) * len(beta_values) * 7,
        },
        "performance_counters": counters,
        "reuse_audit": [dict(row) for row in reuse_audit],
        "reuse_audit_pass": all(
            row["reuse_status"] == "PASS" for row in reuse_audit
        ),
        "Delta_beam_point_array_equal_to_ap0": all(
            bool(row["Delta_beam_array_equal"]) for row in point_comparisons
        ),
        "source_sha256": ap1_source_hashes(),
        "ap0_sha256_before": dict(ap0_before),
        "ap0_sha256_after": dict(ap0_after),
        "ap0_outputs_unchanged": dict(ap0_before) == dict(ap0_after),
        "supervisor_sha256_before": dict(supervisor_before),
        "supervisor_sha256_after": dict(supervisor_after),
        "supervisor_outputs_unchanged": (
            dict(supervisor_before) == dict(supervisor_after)
        ),
        "runtimes_seconds": {
            "scientific_current": scientific_runtime,
            "recorded_family_total": recorded_runtime,
            "workflow_wall_current": wall_runtime,
        },
        "row_counts": {
            GEOMETRY_FILENAME: len(geometries),
            SPECTRA_FILENAME: len(spectra),
            POINT_METRICS_FILENAME: len(point_rows),
            GEOMETRY_SUMMARY_FILENAME: len(summaries),
            PAIRWISE_FILENAME: len(pair_rows),
            POINT_COMPARISON_FILENAME: len(point_comparisons),
            GEOMETRY_COMPARISON_FILENAME: len(geometry_comparisons),
            REUSE_AUDIT_FILENAME: len(reuse_audit),
        },
        "output_directory": str(output_dir.relative_to(ROOT)),
        "diagnostic_figure_paths": list(figure_paths),
        "research_line_separation": {
            "chapter_2_rectangular_anisotropic_only": True,
            "circular_rod_workflow_used": False,
            "isotropic_steel_defaults_used": False,
            "FEM_used": False,
            "three_dimensional_FEM_used": False,
            "damping_used": False,
            "complex_roots_used": False,
        },
        "limitations": {
            "exact_theta_boundary_searched": False,
            "theta3_or_theta4_full_grid_run": False,
            "theta_refinement_run": False,
            "beta_refinement_run": False,
            "MAC_run": False,
            "mode_shapes_run": False,
            "energy_analysis_run": False,
            "FEM_run": False,
            "next_gate_started": False,
        },
    }


def ap1_report(summary: Mapping[str, Any]) -> str:
    science = summary["scientific_result"]
    association = summary["pairwise_descriptive_association"]
    lines = [
        "# AP-1 — same-grid spectral-applicability screening at theta=5 deg",
        "",
        f"**AP-1 theta=5 same-grid screening workflow: {summary['workflow_status']}**",
        "",
        f"Number of geometries within 10% at theta=5: **{science['within_10_percent_at_theta5']} / {science['geometry_count']}**",
        "",
        f"Number of geometries exceeding 10% at theta=5: **{science['exceeding_10_percent_at_theta5']} / {science['geometry_count']}**",
        "",
        f"Number of transitions WITHIN_AT_2_EXCEEDS_AT_5: **{science['classification_transitions_within2_exceeds5']} / {science['geometry_count']}**",
        "",
        "AP-1 uses exactly the AP-0 volume-preserving geometry grid and beta=0,5,...,90 deg. Results are independently sorted spectral positions, not modal descendants.",
        "",
        "`T5` is Chapter-2 state-corrected Timoshenko/generalized-torsion at theta=5 deg; `T0` and rectangular EB/Saint-Venant `EB0` remain theta-zero references. All use HMS/DX-209 elastic, book_slope_clamp and unchanged J_book(beta).",
        "",
        "The fixed AP-0 normalization is unchanged; T5 physical equations retain actual rotated material properties at 5 deg.",
        "",
        "## Geometry-level theta=2 versus theta=5 results",
        "",
        "| geometry | Delta_beam | Delta_orient,2 | Delta_orient,5 | increment | Delta_simpl,2 | Delta_simpl,5 | increment | beta/mode at 5 | theta2 | theta5 | transition |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|",
    ]
    comparisons = {
        row["geometry_id"]: row
        for row in summary["theta2_theta5_geometry_comparison"]
    }
    for row in summary["geometry_results"]:
        comparison = comparisons[row["geometry_id"]]
        lines.append(
            f"| {row['geometry_id']} | {float(row['Delta_beam_max']):.6f} | {float(comparison['Delta_orient_2_max']):.6f} | {float(row['Delta_orient_5_max']):.6f} | {float(comparison['Delta_orient_increment_5_minus_2']):+.6f} | {float(row['Delta_simpl_2_max']):.6f} | {float(row['Delta_simpl_5_max']):.6f} | {float(row['Delta_simpl_increment_5_minus_2']):+.6f} | {float(row['Delta_simpl_5_max_beta']):g}/{int(row['Delta_simpl_5_max_mode'])} | {row['classification_theta2']} | {row['screening_classification_theta5']} | {row['classification_transition']} |"
        )
    lines += [
        "",
        "A threshold exceedance is a scientific result, not a computational failure. The result is limited to the finite grid; no values between beta points or theta=2 and 5 were inferred.",
        "",
        "## Pairwise same-gap descriptive analysis",
        "",
        f"Pooled Spearman `g_pair_T0` versus `abs_delta_g_theta5`: `{float(association['pooled_spearman_g_pair_T0_vs_abs_delta_g_theta5']):.6f}` on `{association['observation_count']}` non-independent screening observations.",
        "",
        "Per-pair Spearman: "
        + ", ".join(
            f"{pair}: {float(value):.6f}"
            for pair, value in association["per_pair_spearman"].items()
        )
        + ".",
        "",
        f"Lowest baseline-gap quartile median abs change: `{float(association['lowest_gap_quartile_median_abs_delta_g_theta5']):.6f}`; remaining observations: `{float(association['remaining_median_abs_delta_g_theta5']):.6f}`.",
        "",
        f"Opened-pair fraction in the lowest-gap quartile: `{float(association['lowest_gap_quartile_opened_pair_fraction']):.6f}`; outside it: `{float(association['remaining_opened_pair_fraction']):.6f}`.",
        "",
        "This is a descriptive pairwise tendency on the finite screening grid. No p-values, independence claim, regression law or causal interpretation is used.",
        "",
        "## Candidates",
        "",
        "```json",
        json.dumps(summary["candidates"], ensure_ascii=False, indent=2, sort_keys=True),
        "```",
        "",
        "Candidates are recorded only; no refinement or next-stage calculation was started.",
        "",
        "## Reuse, quality and performance",
        "",
        f"Reused families: `{summary['families']['reused_family_count']}`; new T5 families: `{summary['families']['new_family_count']}`; resumed: `{summary['families']['resumed_family_count']}`.",
        f"Accepted roots: `{summary['root_quality']['accepted_root_count']}`; maximum determinant residual `{summary['root_quality']['maximum_accepted_determinant_residual']:.12e}`; maximum singular residual `{summary['root_quality']['maximum_accepted_relative_singular_residual']:.12e}`; rejected roots `0`.",
        f"Global anchors `{summary['performance_counters'].get('global_anchor_scans', 0)}`, inventories `{summary['performance_counters'].get('global_inventory_checks', 0)}`, fallbacks `{summary['performance_counters'].get('global_fallback_scans', 0)}`.",
        f"Recorded scientific family runtime `{summary['runtimes_seconds']['recorded_family_total']:.6f} s`; current command wall runtime `{summary['runtimes_seconds']['workflow_wall_current']:.6f} s`.",
        f"Exact reuse audit: `{summary['reuse_audit_pass']}`; Delta_beam point array equals AP-0: `{summary['Delta_beam_point_array_equal_to_ap0']}`.",
        f"AP-0 outputs unchanged: `{summary['ap0_outputs_unchanged']}`; supervisor outputs unchanged: `{summary['supervisor_outputs_unchanged']}`.",
        "",
        "## Limitations",
        "",
        "No exact theta boundary, theta=3/4 full-grid run, theta/beta refinement, MAC, modal descendants, mode shapes, energy analysis, damping, complex roots, FEM or 3D FEM was performed. Circular-rod workflows and article workspaces were not used.",
        "",
        "## Commands",
        "",
        "```text",
        "python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg 5 --smoke",
        "python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg 5 --resume",
        "python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg 5 --reuse-data",
        "```",
    ]
    return "\n".join(lines) + "\n"


def _ap1_output_dir(
    args: argparse.Namespace, config: ScreeningConfiguration
) -> Path:
    uses_default = args.output_dir.resolve() == DEFAULT_OUTPUT_DIR.resolve()
    requested = config.output_dir if uses_default else args.output_dir
    if args.smoke:
        requested = config.smoke_output_dir if uses_default else requested / "smoke"
    return _safe_screening_output(requested)


def run_ap1(
    args: argparse.Namespace, config: ScreeningConfiguration
) -> int:
    output_dir = _ap1_output_dir(args, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    ap0_before = preservation_hashes(
        AP0_OUTPUT_DIR, AP0_PRESERVATION_FILENAMES
    )
    supervisor_before = supervisor_preservation_hashes()
    previous_summary_path = output_dir / SUMMARY_FILENAME
    previous_summary = (
        json.loads(previous_summary_path.read_text(encoding="utf-8"))
        if previous_summary_path.is_file()
        else {}
    )
    family_metadata: dict[str, dict[str, Any]] = {}
    scientific_runtime = 0.0
    try:
        beta_values = beta_grid(smoke=args.smoke)
        if args.reuse_data:
            required = (
                GEOMETRY_FILENAME,
                SPECTRA_FILENAME,
                POINT_METRICS_FILENAME,
                GEOMETRY_SUMMARY_FILENAME,
                PAIRWISE_FILENAME,
                POINT_COMPARISON_FILENAME,
                GEOMETRY_COMPARISON_FILENAME,
                REUSE_AUDIT_FILENAME,
                CANDIDATES_FILENAME,
            )
            missing = [
                str(output_dir / name)
                for name in required
                if not (output_dir / name).is_file()
            ]
            if missing:
                raise FileNotFoundError(
                    "--reuse-data requires saved AP-1 data: " + ", ".join(missing)
                )
            geometries = _geometry_rows_as_numbers(
                _read_csv(output_dir / GEOMETRY_FILENAME)
            )
            spectra = spectra_as_numbers(_read_csv(output_dir / SPECTRA_FILENAME))
            point_rows = _coerce_csv_rows(
                _read_csv(output_dir / POINT_METRICS_FILENAME),
                float_fields=(
                    "mu", "tau", "beta_deg", "Delta_beam",
                    "Delta_orient_5", "Delta_simpl_5", "g_min_0",
                    "g_min_5", "G_nearest_5", "G_open_5", "G_close_5",
                ),
                int_fields=(
                    "Delta_beam_mode", "Delta_orient_5_mode",
                    "Delta_simpl_5_mode", "baseline_min_gap_pair",
                    "G_open_5_pair",
                ),
            )
            summaries = _coerce_csv_rows(
                _read_csv(output_dir / GEOMETRY_SUMMARY_FILENAME),
                float_fields=(
                    "mu", "tau", "Delta_beam_max", "Delta_orient_5_max",
                    "Delta_orient_5_max_beta", "Delta_simpl_5_max",
                    "Delta_simpl_5_max_beta", "minimum_baseline_gap",
                    "minimum_baseline_gap_beta", "maximum_G_nearest_5",
                    "maximum_G_nearest_5_beta", "maximum_G_open_5",
                    "maximum_G_open_5_beta", "Delta_simpl_2_max",
                    "Delta_simpl_increment_5_minus_2",
                ),
                int_fields=(
                    "Delta_orient_5_max_mode", "Delta_simpl_5_max_mode",
                    "minimum_baseline_gap_pair", "maximum_G_open_5_pair",
                ),
            )
            pair_rows = _coerce_csv_rows(
                _read_csv(output_dir / PAIRWISE_FILENAME),
                float_fields=(
                    "mu", "tau", "beta_deg", "g_pair_T0", "g_pair_T2",
                    "g_pair_T5", "delta_g_theta2", "delta_g_theta5",
                    "abs_delta_g_theta2", "abs_delta_g_theta5",
                    "additional_change_theta2_to_theta5",
                    "Delta_orient_5_at_same_point",
                    "Delta_simpl_5_at_same_point",
                ),
                int_fields=("pair_lower_mode", "pair_upper_mode"),
                bool_fields=("pair_opened_at_theta2", "pair_opened_at_theta5"),
            )
            point_comparisons = _coerce_csv_rows(
                _read_csv(output_dir / POINT_COMPARISON_FILENAME),
                float_fields=(
                    "mu", "tau", "beta_deg", "Delta_beam_theta2",
                    "Delta_beam_theta5", "Delta_orient_2", "Delta_orient_5",
                    "Delta_orient_increment_5_minus_2", "Delta_simpl_2",
                    "Delta_simpl_5", "Delta_simpl_increment_5_minus_2",
                ),
                bool_fields=("Delta_beam_array_equal",),
            )
            geometry_comparisons = _coerce_csv_rows(
                _read_csv(output_dir / GEOMETRY_COMPARISON_FILENAME),
                float_fields=(
                    "mu", "tau", "Delta_orient_2_max",
                    "Delta_orient_2_max_beta", "Delta_orient_5_max",
                    "Delta_orient_5_max_beta",
                    "Delta_orient_increment_5_minus_2", "Delta_simpl_2_max",
                    "Delta_simpl_2_max_beta", "Delta_simpl_5_max",
                    "Delta_simpl_5_max_beta",
                    "Delta_simpl_increment_5_minus_2",
                ),
                int_fields=(
                    "Delta_orient_2_max_mode", "Delta_orient_5_max_mode",
                    "Delta_simpl_2_max_mode", "Delta_simpl_5_max_mode",
                ),
            )
            reuse_audit = _coerce_csv_rows(
                _read_csv(output_dir / REUSE_AUDIT_FILENAME),
                float_fields=("theta_deg",),
                int_fields=("beta_count", "root_count"),
                bool_fields=(
                    "frequency_array_equal", "lambda_array_equal",
                    "quality_status_equal",
                ),
            )
            candidates = json.loads(
                (output_dir / CANDIDATES_FILENAME).read_text(encoding="utf-8")
            )
            association = ap1_pairwise_association(pair_rows)
            execution_mode = "reuse_data"
        else:
            supervisor._require_fast_solver_pass(SUPERVISOR_OUTPUT_DIR)
            geometries = load_ap0_geometries(smoke=args.smoke)
            if args.smoke:
                _write_csv(
                    output_dir / GEOMETRY_FILENAME,
                    [asdict(case) for case in geometries],
                )
            else:
                shutil.copyfile(AP0_GEOMETRY_PATH, output_dir / GEOMETRY_FILENAME)
            aggregate = PerformanceCounters()
            transfer_cache = ExactTransferLRU(
                FastSweepSettings().cache_maxsize, counters=aggregate
            )
            spectra = []
            reuse_audit = []
            fingerprint_sources = ap1_source_hashes()
            for case in geometries:
                for model_id in AP1_MODEL_IDS:
                    if model_id in ("T0", "EB0"):
                        rows, audit = load_ap1_ap0_family(
                            case, model_id, beta_values
                        )
                        spectra.extend(rows)
                        reuse_audit.append(audit)
                        continue
                    if case.mu == 0.0 and case.tau == 0.0 and FIGURE_07_PATH.is_file():
                        rows, audit = load_ap1_figure7_t5_family(case, beta_values)
                        spectra.extend(rows)
                        reuse_audit.append(audit)
                        continue
                    rows, metadata, runtime = solve_new_family(
                        case,
                        "T5",
                        beta_values,
                        output_dir=output_dir,
                        resume=args.resume and not args.force_recompute,
                        transfer_cache=transfer_cache,
                        counters=aggregate,
                        smoke=args.smoke,
                        stage_tag="ap1",
                        source_hashes=fingerprint_sources,
                    )
                    for row in rows:
                        if case.mu == 0.0 and case.tau == 0.0:
                            if row["data_origin"] != "global_fallback":
                                row["data_origin"] = "recomputed_missing_source"
                        row["source_file"] = ""
                        row["source_sha256"] = ""
                    spectra.extend(rows)
                    family_metadata[metadata["family_id"]] = metadata
                    scientific_runtime += runtime
            spectra.sort(
                key=lambda row: (
                    str(row["geometry_id"]),
                    float(row["beta_deg"]),
                    AP1_MODEL_IDS.index(str(row["model_id"])),
                    int(row["mode"]),
                )
            )
            _write_csv(output_dir / SPECTRA_FILENAME, spectra)
            _write_csv(output_dir / REUSE_AUDIT_FILENAME, reuse_audit)
            point_rows, pair_rows, point_comparisons = build_ap1_metrics(
                spectra, geometries, beta_values
            )
            summaries, geometry_comparisons = build_ap1_geometry_summary(
                point_rows, geometries
            )
            candidates = select_ap1_candidates(summaries, pair_rows)
            association = ap1_pairwise_association(pair_rows)
            _write_csv(output_dir / POINT_METRICS_FILENAME, point_rows)
            _write_csv(output_dir / GEOMETRY_SUMMARY_FILENAME, summaries)
            _write_csv(output_dir / PAIRWISE_FILENAME, pair_rows)
            _write_csv(output_dir / POINT_COMPARISON_FILENAME, point_comparisons)
            _write_csv(
                output_dir / GEOMETRY_COMPARISON_FILENAME,
                geometry_comparisons,
            )
            _json_write(output_dir / CANDIDATES_FILENAME, candidates)
            execution_mode = (
                "force_recompute" if args.force_recompute else "scientific"
            )

        validate_geometry_manifest(geometries, smoke=args.smoke)
        validate_spectra(
            spectra, geometries, beta_values, model_ids=AP1_MODEL_IDS
        )
        figure_paths = create_ap1_figures(
            summaries, geometry_comparisons, pair_rows, output_dir
        )
        ap0_after = preservation_hashes(
            AP0_OUTPUT_DIR, AP0_PRESERVATION_FILENAMES
        )
        supervisor_after = supervisor_preservation_hashes()
        summary = build_ap1_summary(
            smoke=args.smoke,
            execution_mode=execution_mode,
            geometries=geometries,
            beta_values=beta_values,
            spectra=spectra,
            point_rows=point_rows,
            pair_rows=pair_rows,
            point_comparisons=point_comparisons,
            summaries=summaries,
            geometry_comparisons=geometry_comparisons,
            reuse_audit=reuse_audit,
            candidates=candidates,
            association=association,
            family_metadata=family_metadata,
            previous_summary=previous_summary,
            ap0_before=ap0_before,
            ap0_after=ap0_after,
            supervisor_before=supervisor_before,
            supervisor_after=supervisor_after,
            scientific_runtime=scientific_runtime,
            wall_runtime=time.perf_counter() - started,
            figure_paths=figure_paths,
            output_dir=output_dir,
        )
        _json_write(output_dir / SUMMARY_FILENAME, summary)
        (output_dir / REPORT_FILENAME).write_text(
            ap1_report(summary), encoding="utf-8"
        )
        status = determine_ap1_status(
            smoke=args.smoke,
            geometries=geometries,
            beta_values=beta_values,
            spectra=spectra,
            point_rows=point_rows,
            pair_rows=pair_rows,
            point_comparisons=point_comparisons,
            summaries=summaries,
            geometry_comparisons=geometry_comparisons,
            reuse_audit=reuse_audit,
            ap0_unchanged=ap0_before == ap0_after,
            supervisor_unchanged=supervisor_before == supervisor_after,
            output_dir=output_dir,
        )
        summary["workflow_status"] = status
        summary["runtimes_seconds"]["workflow_wall_current"] = (
            time.perf_counter() - started
        )
        _json_write(output_dir / SUMMARY_FILENAME, summary)
        (output_dir / REPORT_FILENAME).write_text(
            ap1_report(summary), encoding="utf-8"
        )
    except Exception as error:
        failure = {
            "workflow": "AP-1 theta=5 same-grid screening",
            "workflow_status": "FAIL",
            "smoke": bool(args.smoke),
            "error": str(error),
            "git_context": git_context(),
        }
        _json_write(output_dir / SUMMARY_FILENAME, failure)
        print("AP-1 theta=5 same-grid screening workflow: FAIL", file=sys.stderr)
        print(str(error), file=sys.stderr)
        return 1
    science = summary["scientific_result"]
    label = "AP-1 theta=5 smoke workflow" if args.smoke else "AP-1 theta=5 same-grid screening workflow"
    print(f"output_dir={output_dir}")
    print(f"geometry_count={len(geometries)}")
    print(f"beta_point_count={len(beta_values)}")
    print(f"spectrum_row_count={len(spectra)}")
    print(f"pairwise_row_count={len(pair_rows)}")
    print(f"{label}: {summary['workflow_status']}")
    if not args.smoke:
        print(
            "Number of geometries within 10% at theta=5: "
            f"{science['within_10_percent_at_theta5']} / 9"
        )
        print(
            "Number of geometries exceeding 10% at theta=5: "
            f"{science['exceeding_10_percent_at_theta5']} / 9"
        )
        print(
            "Number of transitions WITHIN_AT_2_EXCEEDS_AT_5: "
            f"{science['classification_transitions_within2_exceeds5']} / 9"
        )
    return 0 if summary["workflow_status"] == "PASS" else 1


def accepted_artifact_hashes() -> dict[str, dict[str, str]]:
    return {
        "ap0": preservation_hashes(
            AP0_OUTPUT_DIR, AP0_PRESERVATION_FILENAMES
        ),
        "ap1": preservation_hashes(
            AP1_OUTPUT_DIR, AP1_PRESERVATION_FILENAMES
        ),
        "supervisor": supervisor_preservation_hashes(),
    }


def ap2_source_hashes(theta_deg: float) -> dict[str, str]:
    paths = [
        *SOURCE_PATHS,
        AP0_GEOMETRY_PATH,
        AP0_SPECTRA_PATH,
        AP0_POINT_METRICS_PATH,
        AP0_GEOMETRY_SUMMARY_PATH,
        AP0_SUMMARY_PATH,
        AP1_OUTPUT_DIR / SPECTRA_FILENAME,
        AP1_OUTPUT_DIR / POINT_METRICS_FILENAME,
        AP1_OUTPUT_DIR / GEOMETRY_SUMMARY_FILENAME,
        AP1_OUTPUT_DIR / SUMMARY_FILENAME,
        SUPERVISOR_OUTPUT_DIR / "plot_manifest.json",
    ]
    try:
        paths.append(locate_ap2_supervisor_baseline(theta_deg).source_path)
    except (FileNotFoundError, RuntimeError, ValueError):
        pass
    return {
        str(path.relative_to(ROOT)): _sha256(path)
        for path in paths
        if path.is_file()
    }


def _ap2_expected_output_paths(output_dir: Path, theta_deg: float) -> list[Path]:
    theta = int(theta_deg)
    names = (
        GEOMETRY_FILENAME,
        SPECTRA_FILENAME,
        POINT_METRICS_FILENAME,
        GEOMETRY_SUMMARY_FILENAME,
        PAIRWISE_FILENAME,
        REUSE_AUDIT_FILENAME,
        CANDIDATES_FILENAME,
        SUMMARY_FILENAME,
        REPORT_FILENAME,
        f"screening_delta_simpl_theta{theta}_heatmap.pdf",
        f"screening_delta_simpl_theta{theta}_heatmap.png",
    )
    return [output_dir / name for name in names]


def create_ap2_figure(
    summaries: Sequence[Mapping[str, Any]],
    theta_deg: float,
    output_dir: Path,
) -> list[str]:
    theta = int(theta_deg)
    figure = create_heatmap(
        summaries,
        key=f"Delta_simpl_{theta}_max",
        colorbar_label=rf"$\Delta_{{\mathrm{{simpl}},{theta},\max}}$",
        threshold=APPLICABILITY_THRESHOLD,
    )
    return _save_figure(
        figure, output_dir, f"screening_delta_simpl_theta{theta}_heatmap"
    )


def _ap2_fallback_diagnostics(
    spectra: Sequence[Mapping[str, Any]], theta_deg: float
) -> list[dict[str, Any]]:
    theta = int(theta_deg)
    grouped: dict[tuple[str, float], list[Mapping[str, Any]]] = {}
    for row in spectra:
        if row["model_id"] == f"T{theta}" and row["data_origin"] == "global_fallback":
            grouped.setdefault(
                (str(row["geometry_id"]), float(row["beta_deg"])), []
            ).append(row)
    result: list[dict[str, Any]] = []
    case_lookup = {case.geometry_id: case for case in screening_geometries()}
    for (geometry, beta), rows in sorted(grouped.items()):
        ordered = sorted(rows, key=lambda row: int(row["mode"]))
        case = case_lookup[geometry]
        result.append(
            {
                "theta_deg": float(theta),
                "geometry_id": geometry,
                "mu": case.mu,
                "tau": case.tau,
                "beta_deg": beta,
                "predictor_failure_reason": str(
                    ordered[0].get("global_fallback_reason", "")
                ),
                "accepted_global_inventory_frequency_hz": [
                    float(row["frequency_hz"]) for row in ordered
                ],
                "accepted_determinant_residuals": [
                    float(row["accepted_determinant_residual"])
                    for row in ordered
                ],
                "accepted_relative_singular_residuals": [
                    float(row["accepted_relative_singular_residual"])
                    for row in ordered
                ],
                "quality_statuses": [str(row["quality_status"]) for row in ordered],
            }
        )
    return result


def _classification_decomposition(
    point_rows: Sequence[Mapping[str, Any]], theta_deg: float
) -> dict[str, int]:
    theta = int(theta_deg)
    orient_key = f"Delta_orient_{theta}"
    simpl_key = f"Delta_simpl_{theta}"
    a = sum(float(row[orient_key]) > APPLICABILITY_THRESHOLD for row in point_rows)
    b = sum(
        float(row[orient_key]) <= APPLICABILITY_THRESHOLD
        and float(row[simpl_key]) > APPLICABILITY_THRESHOLD
        for row in point_rows
    )
    c = sum(float(row[simpl_key]) <= APPLICABILITY_THRESHOLD for row in point_rows)
    if a + b + c != len(point_rows):
        raise RuntimeError("AP-2 classification decomposition is incomplete")
    return {
        "A_orientation_effect_exceeds": int(a),
        "B_orientation_within_but_simplification_exceeds": int(b),
        "C_simplification_within": int(c),
        "A_plus_B": int(a + b),
        "total": int(a + b + c),
    }


def determine_ap2_status(
    *,
    smoke: bool,
    theta_deg: float,
    geometries: Sequence[GeometryCase],
    beta_values: NDArray[np.float64],
    spectra: Sequence[Mapping[str, Any]],
    point_rows: Sequence[Mapping[str, Any]],
    pair_rows: Sequence[Mapping[str, Any]],
    summaries: Sequence[Mapping[str, Any]],
    reuse_audit: Sequence[Mapping[str, Any]],
    preservation_unchanged: bool,
    output_dir: Path,
) -> str:
    geometry_count = 3 if smoke else 9
    beta_count = 3 if smoke else 19
    reuse_count = 7 if smoke else 19
    acceptable_audits = all(
        row["reuse_status"]
        in ("PASS", "SUPERVISOR_BASELINE_REUSE_REJECTED")
        for row in reuse_audit
    )
    gates = (
        int(theta_deg) in (3, 4),
        len(geometries) == geometry_count,
        len(beta_values) == beta_count,
        len(spectra) == geometry_count * beta_count * 3 * GUARD_ROOT_COUNT,
        len(point_rows) == geometry_count * beta_count,
        len(pair_rows) == geometry_count * beta_count * COMPARED_ROOT_COUNT,
        len(summaries) == geometry_count,
        len(reuse_audit) == reuse_count,
        acceptable_audits,
        all(row["quality_status"] == "PASS" for row in spectra),
        preservation_unchanged,
        all(path.is_file() for path in _ap2_expected_output_paths(output_dir, theta_deg)),
    )
    return "PASS" if all(gates) else "FAIL"


def build_ap2_summary(
    *,
    smoke: bool,
    theta_deg: float,
    execution_mode: str,
    geometries: Sequence[GeometryCase],
    beta_values: NDArray[np.float64],
    spectra: Sequence[Mapping[str, Any]],
    point_rows: Sequence[Mapping[str, Any]],
    pair_rows: Sequence[Mapping[str, Any]],
    summaries: Sequence[Mapping[str, Any]],
    reuse_audit: Sequence[Mapping[str, Any]],
    candidates: Mapping[str, Any],
    association: Mapping[str, Any],
    family_metadata: Mapping[str, Mapping[str, Any]],
    previous_summary: Mapping[str, Any],
    preservation_before: Mapping[str, Mapping[str, str]],
    preservation_after: Mapping[str, Mapping[str, str]],
    scientific_runtime: float,
    wall_runtime: float,
    figure_paths: Sequence[str],
    output_dir: Path,
) -> dict[str, Any]:
    theta = int(theta_deg)
    if execution_mode == "reuse_data":
        families = dict(previous_summary.get("families", {}))
        counters = dict(previous_summary.get("performance_counters", {}))
        recorded_runtime = float(
            previous_summary.get("runtimes_seconds", {}).get(
                "recorded_family_total", 0.0
            )
        )
    else:
        families = {
            "new_family_count": len(family_metadata),
            "reused_family_count": sum(
                row["reuse_status"] == "PASS" for row in reuse_audit
            ),
            "rejected_reuse_recomputed_family_count": sum(
                row["reuse_status"] == "SUPERVISOR_BASELINE_REUSE_REJECTED"
                for row in reuse_audit
            ),
            "resumed_family_count": sum(
                bool(item.get("resumed_from_checkpoint"))
                for item in family_metadata.values()
            ),
            "metadata": dict(family_metadata),
        }
        counters = _aggregate_family_counters(family_metadata)
        recorded_runtime = sum(
            float(item.get("recorded_runtime_seconds", 0.0))
            for item in family_metadata.values()
        )
    families.setdefault("new_family_count", 0)
    families.setdefault("reused_family_count", len(reuse_audit))
    families.setdefault("rejected_reuse_recomputed_family_count", 0)
    families.setdefault("resumed_family_count", 0)
    families.setdefault("metadata", {})
    simpl_key = f"Delta_simpl_{theta}"
    class_key = f"screening_classification_theta{theta}"
    point_within = sum(
        float(row[simpl_key]) <= APPLICABILITY_THRESHOLD for row in point_rows
    )
    family_within = sum(
        row[class_key] == "WITHIN_10_PERCENT_ON_SCREENING_GRID"
        for row in summaries
    )
    decomposition = _classification_decomposition(point_rows, theta)
    quality = _quality_summary(spectra)
    fallbacks = _ap2_fallback_diagnostics(spectra, theta)
    geometry_hash = _sha256(output_dir / GEOMETRY_FILENAME)
    return {
        "workflow": f"AP-2 same-grid spectral-applicability screening at theta={theta} deg",
        "workflow_status": "PENDING",
        "smoke": smoke,
        "smoke_is_scientific_result": False,
        "execution_mode": execution_mode,
        "output_directory": str(output_dir.relative_to(ROOT)),
        "git_context": git_context(),
        "configuration": {
            "target_theta_deg": float(theta),
            "models": list(AP2_MODEL_IDS[theta]),
            "material": "HMS/DX-209 elastic",
            "clamp": CLAMP,
            "joint": "unchanged J_book(beta)",
            "geometry_source": str(AP0_GEOMETRY_PATH.relative_to(ROOT)),
            "geometry_source_sha256": _sha256(AP0_GEOMETRY_PATH),
            "geometry_manifest_sha256": geometry_hash,
            "geometry_manifest_content_identical_to_ap0": (
                smoke or geometry_hash == _sha256(AP0_GEOMETRY_PATH)
            ),
            "fixed_lambda_normalization": normalization_contract(
                geometries[0], f"T{theta}"
            ),
        },
        "grid": {
            "mu": list(MU_VALUES),
            "tau": list(TAU_VALUES),
            "beta_deg": [float(value) for value in beta_values],
            "geometry_count": len(geometries),
            "beta_count": len(beta_values),
            "root_count": GUARD_ROOT_COUNT,
            "compared_root_count": COMPARED_ROOT_COUNT,
        },
        "row_counts": {
            GEOMETRY_FILENAME: len(geometries),
            SPECTRA_FILENAME: len(spectra),
            POINT_METRICS_FILENAME: len(point_rows),
            GEOMETRY_SUMMARY_FILENAME: len(summaries),
            PAIRWISE_FILENAME: len(pair_rows),
            REUSE_AUDIT_FILENAME: len(reuse_audit),
        },
        "families": families,
        "root_quality": quality,
        "performance_counters": counters,
        "new_ap2_fallbacks": fallbacks,
        "new_ap2_fallback_count": len(fallbacks),
        "historical_reused_fallbacks_counted_as_new": False,
        "runtimes_seconds": {
            "recorded_family_total": recorded_runtime,
            "scientific_current": scientific_runtime,
            "workflow_wall_current": wall_runtime,
        },
        "scientific_result": {
            "scope": "FINITE_AP2_GRID_ONLY",
            f"pointwise_within_theta{theta}": point_within,
            f"pointwise_exceeding_theta{theta}": len(point_rows) - point_within,
            f"families_uniformly_within_theta{theta}": family_within,
            f"families_exceeding_theta{theta}": len(summaries) - family_within,
            "classification_decomposition": decomposition,
        },
        "geometry_table": [dict(row) for row in summaries],
        "pairwise_descriptive_association": dict(association),
        "candidates": dict(candidates),
        "reuse_audit": [dict(row) for row in reuse_audit],
        "reuse_audit_pass_or_recomputed": all(
            row["reuse_status"]
            in ("PASS", "SUPERVISOR_BASELINE_REUSE_REJECTED")
            for row in reuse_audit
        ),
        "preservation_sha256_before": dict(preservation_before),
        "preservation_sha256_after": dict(preservation_after),
        "ap0_outputs_unchanged": preservation_before["ap0"]
        == preservation_after["ap0"],
        "ap1_outputs_unchanged": preservation_before["ap1"]
        == preservation_after["ap1"],
        "supervisor_outputs_unchanged": preservation_before["supervisor"]
        == preservation_after["supervisor"],
        "source_sha256": ap2_source_hashes(theta),
        "diagnostic_figure_paths": list(figure_paths),
        "limitations": {
            "exact_theta_boundary_searched": False,
            "theta_interpolation_performed": False,
            "beta_refinement_performed": False,
            "modal_descendants_tracked": False,
            "mac_or_shapes_or_energy_performed": False,
            "fem_or_3d_fem_performed": False,
        },
    }


def ap2_report(summary: Mapping[str, Any]) -> str:
    theta = int(summary["configuration"]["target_theta_deg"])
    science = summary["scientific_result"]
    quality = summary["root_quality"]
    families = summary["families"]
    association = summary["pairwise_descriptive_association"]
    lines = [
        f"# AP-2 same-grid screening at theta={theta} deg",
        "",
        f"**AP-2 theta={theta} computational workflow: {summary['workflow_status']}**",
        "",
        "Theta is the material-axis orientation relative to each local rod axis; beta is the geometric joint angle. T"
        f"{theta} is a more complete one-dimensional Chapter-2 reference model, not an exact 3D model.",
        "",
        "## Computational audit",
        "",
        f"Rows: spectra `{summary['row_counts'][SPECTRA_FILENAME]}`, points `{summary['row_counts'][POINT_METRICS_FILENAME]}`, pairs `{summary['row_counts'][PAIRWISE_FILENAME]}`.",
        f"Families new/reused/resumed: `{families['new_family_count']}/{families['reused_family_count']}/{families['resumed_family_count']}`.",
        f"Accepted/rejected roots: `{quality['accepted_root_count']}/{quality['rejected_root_count']}`; determinant maximum `{quality['maximum_accepted_determinant_residual']:.6e}`; singular maximum `{quality['maximum_accepted_relative_singular_residual']:.6e}`.",
        f"Quality bases: `{quality['quality_basis_counts']}`.",
        f"New AP-2 fallbacks: `{summary['new_ap2_fallback_count']}`; recorded family runtime `{summary['runtimes_seconds']['recorded_family_total']:.6f} s`; current scientific runtime `{summary['runtimes_seconds']['scientific_current']:.6f} s`.",
        f"Anchors/inventories/local roots/transfer exponentials/cache hit rate: `{summary['performance_counters'].get('global_anchor_scans', 0)}/{summary['performance_counters'].get('global_inventory_checks', 0)}/{summary['performance_counters'].get('local_root_solves', 0)}/{summary['performance_counters'].get('transfer_expm_evaluations', 0)}/{summary['performance_counters'].get('transfer_cache_hit_rate', 0.0):.6f}`.",
        f"AP-0/AP-1/supervisor unchanged: `{summary['ap0_outputs_unchanged']}/{summary['ap1_outputs_unchanged']}/{summary['supervisor_outputs_unchanged']}`.",
        "",
        "## Scientific result",
        "",
        f"Pointwise within/exceeding: `{science[f'pointwise_within_theta{theta}']}/{science[f'pointwise_exceeding_theta{theta}']}` out of `{summary['row_counts'][POINT_METRICS_FILENAME]}`.",
        f"Families uniformly within/exceeding: `{science[f'families_uniformly_within_theta{theta}']}/{science[f'families_exceeding_theta{theta}']}` out of `{summary['row_counts'][GEOMETRY_SUMMARY_FILENAME]}`.",
        f"Classification decomposition A/B/C: `{science['classification_decomposition']['A_orientation_effect_exceeds']}/{science['classification_decomposition']['B_orientation_within_but_simplification_exceeds']}/{science['classification_decomposition']['C_simplification_within']}`.",
        "",
        "The decomposition is a classification decomposition, not an additive error budget. Curves and metrics use independently sorted spectral positions, not tracked physical modal descendants.",
        "",
        "| mu | tau | Delta orient max | Delta simpl max | beta/k at simpl max | family status |",
        "|---:|---:|---:|---:|---:|---|",
    ]
    for row in summary["geometry_table"]:
        lines.append(
            f"| {float(row['mu']):.2f} | {float(row['tau']):+.2f} | {float(row[f'Delta_orient_{theta}_max']):.6f} | {float(row[f'Delta_simpl_{theta}_max']):.6f} | {float(row[f'Delta_simpl_{theta}_max_beta']):g}/{int(row[f'Delta_simpl_{theta}_max_mode'])} | {row[f'screening_classification_theta{theta}']} |"
        )
    lines.extend(
        [
            "",
            "## New AP-2 fallback diagnostics",
            "",
            "Historical fallback metadata in reused AP-0/supervisor arrays is not counted here.",
            "",
            "| mu | tau | beta | predictor/local failure reason |",
            "|---:|---:|---:|---|",
        ]
    )
    for row in summary["new_ap2_fallbacks"]:
        lines.append(
            f"| {float(row['mu']):.2f} | {float(row['tau']):+.2f} | {float(row['beta_deg']):g} | {row['predictor_failure_reason']} |"
        )
    lines.extend(
        [
        "",
        "## Pairwise same-pair diagnostic",
        "",
        f"Pooled Spearman: `{association[f'pooled_spearman_g_pair_T0_vs_abs_delta_g_theta{theta}']:.9f}` on `{association['observation_count']}` finite-grid observations. No p-values, independence claim, regression law, or causal claim is made.",
        f"Per-pair Spearman: `{association['per_pair_spearman']}`.",
        f"Lowest-gap/remaining median absolute changes: `{association[f'lowest_gap_quartile_median_abs_delta_g_theta{theta}']:.9f}/{association[f'remaining_median_abs_delta_g_theta{theta}']:.9f}`; opened fractions: `{association['lowest_gap_quartile_opened_pair_fraction']:.9f}/{association['remaining_opened_pair_fraction']:.9f}`.",
        "",
        "## Limitations",
        "",
        "No exact theta boundary, theta interpolation, beta refinement, MAC, mode shapes, energy fractions, FEM, 3D FEM, damping, or circular-rod workflow was run.",
        "",
        "## Commands",
        "",
        "```text",
        f"python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg {theta} --smoke",
        f"python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg {theta} --resume",
        f"python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --target-theta-deg {theta} --reuse-data",
        "```",
        ]
    )
    return "\n".join(lines) + "\n"


def run_ap2(args: argparse.Namespace, config: ScreeningConfiguration) -> int:
    theta = int(config.target_theta_deg)
    output_dir = _ap1_output_dir(args, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    preservation_before = accepted_artifact_hashes()
    previous_summary_path = output_dir / SUMMARY_FILENAME
    previous_summary = (
        json.loads(previous_summary_path.read_text(encoding="utf-8"))
        if previous_summary_path.is_file()
        else {}
    )
    family_metadata: dict[str, dict[str, Any]] = {}
    scientific_runtime = 0.0
    try:
        beta_values = beta_grid(smoke=args.smoke)
        if args.reuse_data:
            required = (
                GEOMETRY_FILENAME,
                SPECTRA_FILENAME,
                POINT_METRICS_FILENAME,
                GEOMETRY_SUMMARY_FILENAME,
                PAIRWISE_FILENAME,
                REUSE_AUDIT_FILENAME,
                CANDIDATES_FILENAME,
            )
            missing = [
                str(output_dir / name)
                for name in required
                if not (output_dir / name).is_file()
            ]
            if missing:
                raise FileNotFoundError(
                    "--reuse-data requires saved AP-2 data: " + ", ".join(missing)
                )
            geometries = _geometry_rows_as_numbers(
                _read_csv(output_dir / GEOMETRY_FILENAME)
            )
            spectra = spectra_as_numbers(_read_csv(output_dir / SPECTRA_FILENAME))
            point_rows = _coerce_csv_rows(
                _read_csv(output_dir / POINT_METRICS_FILENAME),
                float_fields=(
                    "mu", "tau", "beta_deg", "Delta_beam",
                    f"Delta_orient_{theta}", f"Delta_simpl_{theta}",
                    "g_min_0", f"g_min_{theta}", f"G_nearest_{theta}",
                    f"G_open_{theta}", f"G_close_{theta}",
                ),
                int_fields=(
                    "Delta_beam_mode", f"Delta_orient_{theta}_mode",
                    f"Delta_simpl_{theta}_mode", "baseline_min_gap_pair",
                    f"G_open_{theta}_pair",
                ),
            )
            summaries = _coerce_csv_rows(
                _read_csv(output_dir / GEOMETRY_SUMMARY_FILENAME),
                float_fields=(
                    "mu", "tau", "Delta_beam_max",
                    f"Delta_orient_{theta}_max",
                    f"Delta_orient_{theta}_max_beta",
                    f"Delta_simpl_{theta}_max",
                    f"Delta_simpl_{theta}_max_beta", "minimum_baseline_gap",
                    "minimum_baseline_gap_beta", f"maximum_G_nearest_{theta}",
                    f"maximum_G_nearest_{theta}_beta", f"maximum_G_open_{theta}",
                    f"maximum_G_open_{theta}_beta", "Delta_simpl_2_max",
                    f"Delta_simpl_increment_{theta}_minus_2",
                ),
                int_fields=(
                    f"Delta_orient_{theta}_max_mode",
                    f"Delta_simpl_{theta}_max_mode", "minimum_baseline_gap_pair",
                    f"maximum_G_open_{theta}_pair",
                ),
            )
            pair_rows = _coerce_csv_rows(
                _read_csv(output_dir / PAIRWISE_FILENAME),
                float_fields=(
                    "mu", "tau", "beta_deg", "g_pair_T0", "g_pair_T2",
                    f"g_pair_T{theta}", "delta_g_theta2",
                    f"delta_g_theta{theta}", "abs_delta_g_theta2",
                    f"abs_delta_g_theta{theta}", f"change_g_2_to_{theta}",
                    f"Delta_orient_{theta}_at_same_point",
                    f"Delta_simpl_{theta}_at_same_point",
                ),
                int_fields=("pair_lower_mode", "pair_upper_mode"),
                bool_fields=(
                    "pair_opened_at_theta2", f"pair_opened_at_theta{theta}"
                ),
            )
            reuse_audit = _coerce_csv_rows(
                _read_csv(output_dir / REUSE_AUDIT_FILENAME),
                float_fields=("theta_deg",),
                int_fields=("beta_count", "root_count"),
                bool_fields=(
                    "frequency_array_equal", "lambda_array_equal",
                    "quality_status_equal",
                ),
            )
            candidates = json.loads(
                (output_dir / CANDIDATES_FILENAME).read_text(encoding="utf-8")
            )
            association = ap2_pairwise_association(pair_rows, theta)
            execution_mode = "reuse_data"
        else:
            supervisor._require_fast_solver_pass(SUPERVISOR_OUTPUT_DIR)
            geometries = load_ap0_geometries(smoke=args.smoke)
            if args.smoke:
                _write_csv(
                    output_dir / GEOMETRY_FILENAME,
                    [asdict(case) for case in geometries],
                )
            else:
                shutil.copyfile(AP0_GEOMETRY_PATH, output_dir / GEOMETRY_FILENAME)
            aggregate = PerformanceCounters()
            transfer_cache = ExactTransferLRU(
                FastSweepSettings().cache_maxsize, counters=aggregate
            )
            spectra: list[dict[str, Any]] = []
            reuse_audit: list[dict[str, Any]] = []
            source_hashes = ap2_source_hashes(theta)
            for case in geometries:
                for model_id in config.model_ids:
                    if model_id in ("T0", "EB0"):
                        rows, audit = load_ap1_ap0_family(case, model_id, beta_values)
                        spectra.extend(rows)
                        reuse_audit.append(audit)
                        continue
                    reuse_rejection_reason = ""
                    if case.mu == 0.0 and case.tau == 0.0:
                        rows, audit = attempt_ap2_supervisor_reuse(
                            case, model_id, beta_values
                        )
                        reuse_audit.append(audit)
                        if rows is not None:
                            spectra.extend(rows)
                            continue
                        reuse_rejection_reason = str(
                            audit.get("reuse_rejection_reason", "")
                        )
                    rows, metadata, runtime = solve_new_family(
                        case,
                        model_id,
                        beta_values,
                        output_dir=output_dir,
                        resume=args.resume and not args.force_recompute,
                        transfer_cache=transfer_cache,
                        counters=aggregate,
                        smoke=args.smoke,
                        stage_tag=f"ap2_theta{theta}",
                        source_hashes=source_hashes,
                    )
                    for row in rows:
                        if reuse_rejection_reason and row["data_origin"] != "global_fallback":
                            row["data_origin"] = "recomputed_missing_or_rejected_supervisor_source"
                        row["source_file"] = ""
                        row["source_sha256"] = ""
                    spectra.extend(rows)
                    family_metadata[metadata["family_id"]] = metadata
                    scientific_runtime += runtime
            spectra.sort(
                key=lambda row: (
                    str(row["geometry_id"]),
                    float(row["beta_deg"]),
                    config.model_ids.index(str(row["model_id"])),
                    int(row["mode"]),
                )
            )
            _write_csv(output_dir / SPECTRA_FILENAME, spectra)
            _write_csv(output_dir / REUSE_AUDIT_FILENAME, reuse_audit)
            point_rows, pair_rows = build_ap2_metrics(
                spectra, geometries, beta_values, theta
            )
            summaries = build_ap2_geometry_summary(point_rows, geometries, theta)
            candidates = select_ap2_candidates(summaries, theta)
            association = ap2_pairwise_association(pair_rows, theta)
            _write_csv(output_dir / POINT_METRICS_FILENAME, point_rows)
            _write_csv(output_dir / GEOMETRY_SUMMARY_FILENAME, summaries)
            _write_csv(output_dir / PAIRWISE_FILENAME, pair_rows)
            _json_write(output_dir / CANDIDATES_FILENAME, candidates)
            execution_mode = (
                "force_recompute" if args.force_recompute else "scientific"
            )

        validate_geometry_manifest(geometries, smoke=args.smoke)
        validate_spectra(
            spectra, geometries, beta_values, model_ids=config.model_ids
        )
        figure_paths = create_ap2_figure(summaries, theta, output_dir)
        preservation_after = accepted_artifact_hashes()
        summary = build_ap2_summary(
            smoke=args.smoke,
            theta_deg=theta,
            execution_mode=execution_mode,
            geometries=geometries,
            beta_values=beta_values,
            spectra=spectra,
            point_rows=point_rows,
            pair_rows=pair_rows,
            summaries=summaries,
            reuse_audit=reuse_audit,
            candidates=candidates,
            association=association,
            family_metadata=family_metadata,
            previous_summary=previous_summary,
            preservation_before=preservation_before,
            preservation_after=preservation_after,
            scientific_runtime=scientific_runtime,
            wall_runtime=time.perf_counter() - started,
            figure_paths=figure_paths,
            output_dir=output_dir,
        )
        _json_write(output_dir / SUMMARY_FILENAME, summary)
        (output_dir / REPORT_FILENAME).write_text(
            ap2_report(summary), encoding="utf-8"
        )
        status = determine_ap2_status(
            smoke=args.smoke,
            theta_deg=theta,
            geometries=geometries,
            beta_values=beta_values,
            spectra=spectra,
            point_rows=point_rows,
            pair_rows=pair_rows,
            summaries=summaries,
            reuse_audit=reuse_audit,
            preservation_unchanged=preservation_before == preservation_after,
            output_dir=output_dir,
        )
        summary["workflow_status"] = status
        summary["runtimes_seconds"]["workflow_wall_current"] = (
            time.perf_counter() - started
        )
        _json_write(output_dir / SUMMARY_FILENAME, summary)
        (output_dir / REPORT_FILENAME).write_text(
            ap2_report(summary), encoding="utf-8"
        )
    except Exception as error:
        failure = {
            "workflow": f"AP-2 theta={theta} same-grid screening",
            "workflow_status": "FAIL",
            "smoke": bool(args.smoke),
            "error": str(error),
            "git_context": git_context(),
        }
        _json_write(output_dir / SUMMARY_FILENAME, failure)
        print(f"AP-2 theta={theta} computational workflow: FAIL", file=sys.stderr)
        print(str(error), file=sys.stderr)
        return 1
    print(f"output_dir={output_dir}")
    print(f"geometry_count={len(geometries)}")
    print(f"beta_point_count={len(beta_values)}")
    print(f"spectrum_row_count={len(spectra)}")
    print(f"pairwise_row_count={len(pair_rows)}")
    print(f"AP-2 theta={theta} computational workflow: {summary['workflow_status']}")
    if not args.smoke:
        science = summary["scientific_result"]
        print(
            f"Pointwise within/exceeding at theta={theta}: "
            f"{science[f'pointwise_within_theta{theta}']}/"
            f"{science[f'pointwise_exceeding_theta{theta}']} out of 171"
        )
        print(
            f"Families within/exceeding at theta={theta}: "
            f"{science[f'families_uniformly_within_theta{theta}']}/"
            f"{science[f'families_exceeding_theta{theta}']} out of 9"
        )
    return 0 if summary["workflow_status"] == "PASS" else 1


def _sampled_status(values: Sequence[float]) -> tuple[list[str], str, float | str, bool, bool]:
    statuses = [
        "W" if float(value) <= APPLICABILITY_THRESHOLD else "E"
        for value in values
    ]
    first_index = next(
        (index for index, status in enumerate(statuses) if status == "E"), None
    )
    first_theta: float | str = (
        "" if first_index is None else SAMPLED_THETAS_DEG[first_index]
    )
    reentry = bool(
        first_index is not None and "W" in statuses[first_index + 1 :]
    )
    nondecreasing = all(
        float(right) + SAMPLED_NONDECREASING_ABS_TOLERANCE >= float(left)
        for left, right in zip(values, values[1:])
    )
    if first_index is None:
        label = "WITHIN_AT_ALL_SAMPLED_THETAS"
    else:
        label = f"FIRST_EXCEEDS_AT_{int(SAMPLED_THETAS_DEG[first_index])}_DEG"
    return statuses, label, first_theta, reentry, nondecreasing


def _sampled_source_rows() -> dict[int, dict[str, list[dict[str, str]]]]:
    directories = {
        2: AP0_OUTPUT_DIR,
        3: AP2_OUTPUT_DIRS[3],
        4: AP2_OUTPUT_DIRS[4],
        5: AP1_OUTPUT_DIR,
    }
    result: dict[int, dict[str, list[dict[str, str]]]] = {}
    for theta, directory in directories.items():
        required = (
            SPECTRA_FILENAME,
            POINT_METRICS_FILENAME,
            GEOMETRY_SUMMARY_FILENAME,
        )
        missing = [name for name in required if not (directory / name).is_file()]
        if missing:
            raise FileNotFoundError(
                f"theta={theta}: accepted source files missing: {', '.join(missing)}"
            )
        result[theta] = {
            "spectra": _read_csv(directory / SPECTRA_FILENAME),
            "points": _read_csv(directory / POINT_METRICS_FILENAME),
            "summaries": _read_csv(directory / GEOMETRY_SUMMARY_FILENAME),
        }
    return result


def build_sampled_point_comparison(
    sources: Mapping[int, Mapping[str, Sequence[Mapping[str, Any]]]],
) -> list[dict[str, Any]]:
    lookup = {
        theta: {
            (str(row["geometry_id"]), float(row["beta_deg"])): row
            for row in source["points"]
        }
        for theta, source in sources.items()
    }
    rows: list[dict[str, Any]] = []
    for case in screening_geometries():
        for beta in beta_grid():
            key = (case.geometry_id, float(beta))
            source_rows = {theta: lookup[theta][key] for theta in (2, 3, 4, 5)}
            beam_values = [float(source_rows[theta]["Delta_beam"]) for theta in (2, 3, 4, 5)]
            if not all(value == beam_values[0] for value in beam_values[1:]):
                raise RuntimeError(
                    f"{case.geometry_id}/beta={beta:g}: Delta_beam is not exact across sampled theta"
                )
            orient = {
                2: float(source_rows[2]["Delta_orient"]),
                3: float(source_rows[3]["Delta_orient_3"]),
                4: float(source_rows[4]["Delta_orient_4"]),
                5: float(source_rows[5]["Delta_orient_5"]),
            }
            simpl = {
                2: float(source_rows[2]["Delta_simpl"]),
                3: float(source_rows[3]["Delta_simpl_3"]),
                4: float(source_rows[4]["Delta_simpl_4"]),
                5: float(source_rows[5]["Delta_simpl_5"]),
            }
            orient_modes = {
                2: int(source_rows[2]["Delta_orient_mode"]),
                3: int(source_rows[3]["Delta_orient_3_mode"]),
                4: int(source_rows[4]["Delta_orient_4_mode"]),
                5: int(source_rows[5]["Delta_orient_5_mode"]),
            }
            simpl_modes = {
                2: int(source_rows[2]["Delta_simpl_mode"]),
                3: int(source_rows[3]["Delta_simpl_3_mode"]),
                4: int(source_rows[4]["Delta_simpl_4_mode"]),
                5: int(source_rows[5]["Delta_simpl_5_mode"]),
            }
            statuses, label, first, reentry, simpl_nondecreasing = _sampled_status(
                [simpl[theta] for theta in (2, 3, 4, 5)]
            )
            _, _, _, _, orient_nondecreasing = _sampled_status(
                [orient[theta] for theta in (2, 3, 4, 5)]
            )
            row: dict[str, Any] = {
                "geometry_id": case.geometry_id,
                "mu": case.mu,
                "tau": case.tau,
                "beta_deg": float(beta),
                "Delta_beam": beam_values[0],
            }
            for theta in (2, 3, 4, 5):
                row[f"Delta_orient_theta{theta}"] = orient[theta]
                row[f"Delta_simpl_theta{theta}"] = simpl[theta]
                row[f"argmax_k_orient_theta{theta}"] = orient_modes[theta]
                row[f"argmax_k_simpl_theta{theta}"] = simpl_modes[theta]
                row[f"status_theta{theta}"] = (
                    "WITHIN_POINT" if statuses[theta - 2] == "W" else "EXCEEDS_POINT"
                )
            for left, right in ((2, 3), (3, 4), (4, 5)):
                row[f"increment_orient_{left}_to_{right}"] = orient[right] - orient[left]
                row[f"increment_simpl_{left}_to_{right}"] = simpl[right] - simpl[left]
            row.update(
                {
                    "sampled_theta_status_sequence": "-".join(statuses),
                    "sampled_theta_status_label": label,
                    "first_exceeding_sampled_theta_deg": first,
                    "has_reentry_after_exceedance": reentry,
                    "delta_orient_sampled_nondecreasing": orient_nondecreasing,
                    "delta_simpl_sampled_nondecreasing": simpl_nondecreasing,
                }
            )
            rows.append(row)
    return rows


def build_sampled_geometry_comparison(
    sources: Mapping[int, Mapping[str, Sequence[Mapping[str, Any]]]],
) -> list[dict[str, Any]]:
    lookup = {
        theta: {
            str(row["geometry_id"]): row for row in source["summaries"]
        }
        for theta, source in sources.items()
    }
    rows: list[dict[str, Any]] = []
    for case in screening_geometries():
        source_rows = {theta: lookup[theta][case.geometry_id] for theta in (2, 3, 4, 5)}
        simpl = {
            2: float(source_rows[2]["Delta_simpl_max"]),
            3: float(source_rows[3]["Delta_simpl_3_max"]),
            4: float(source_rows[4]["Delta_simpl_4_max"]),
            5: float(source_rows[5]["Delta_simpl_5_max"]),
        }
        orient = {
            2: float(source_rows[2]["Delta_orient_max"]),
            3: float(source_rows[3]["Delta_orient_3_max"]),
            4: float(source_rows[4]["Delta_orient_4_max"]),
            5: float(source_rows[5]["Delta_orient_5_max"]),
        }
        simpl_beta = {
            2: float(source_rows[2]["Delta_simpl_max_beta"]),
            3: float(source_rows[3]["Delta_simpl_3_max_beta"]),
            4: float(source_rows[4]["Delta_simpl_4_max_beta"]),
            5: float(source_rows[5]["Delta_simpl_5_max_beta"]),
        }
        simpl_mode = {
            2: int(source_rows[2]["Delta_simpl_max_mode"]),
            3: int(source_rows[3]["Delta_simpl_3_max_mode"]),
            4: int(source_rows[4]["Delta_simpl_4_max_mode"]),
            5: int(source_rows[5]["Delta_simpl_5_max_mode"]),
        }
        orient_beta = {
            2: float(source_rows[2]["Delta_orient_max_beta"]),
            3: float(source_rows[3]["Delta_orient_3_max_beta"]),
            4: float(source_rows[4]["Delta_orient_4_max_beta"]),
            5: float(source_rows[5]["Delta_orient_5_max_beta"]),
        }
        orient_mode = {
            2: int(source_rows[2]["Delta_orient_max_mode"]),
            3: int(source_rows[3]["Delta_orient_3_max_mode"]),
            4: int(source_rows[4]["Delta_orient_4_max_mode"]),
            5: int(source_rows[5]["Delta_orient_5_max_mode"]),
        }
        statuses, _, first, reentry, nondecreasing = _sampled_status(
            [simpl[theta] for theta in (2, 3, 4, 5)]
        )
        row: dict[str, Any] = {
            "geometry_id": case.geometry_id,
            "mu": case.mu,
            "tau": case.tau,
        }
        for theta in (2, 3, 4, 5):
            row[f"Delta_simpl_max_theta{theta}"] = simpl[theta]
            row[f"Delta_orient_max_theta{theta}"] = orient[theta]
            row[f"beta_at_simpl_max_theta{theta}"] = simpl_beta[theta]
            row[f"k_at_simpl_max_theta{theta}"] = simpl_mode[theta]
            row[f"beta_at_orient_max_theta{theta}"] = orient_beta[theta]
            row[f"k_at_orient_max_theta{theta}"] = orient_mode[theta]
            row[f"family_status_theta{theta}"] = (
                "WITHIN_FAMILY" if statuses[theta - 2] == "W" else "EXCEEDS_FAMILY"
            )
        for left, right in ((2, 3), (3, 4), (4, 5)):
            row[f"increment_simpl_max_{left}_to_{right}"] = simpl[right] - simpl[left]
        row.update(
            {
                "sampled_theta_family_status_sequence": "-".join(statuses),
                "first_exceeding_sampled_theta_deg": first,
                "has_reentry_after_exceedance": reentry,
                "sampled_nondecreasing": nondecreasing,
            }
        )
        rows.append(row)
    return rows


def build_sampled_classification_counts(
    point_rows: Sequence[Mapping[str, Any]],
    geometry_rows: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    counts: list[dict[str, Any]] = []
    for theta in (2, 3, 4, 5):
        point_within = sum(
            row[f"status_theta{theta}"] == "WITHIN_POINT" for row in point_rows
        )
        family_within = sum(
            row[f"family_status_theta{theta}"] == "WITHIN_FAMILY"
            for row in geometry_rows
        )
        counts.append(
            {
                "theta_deg": float(theta),
                "pointwise_within_count": point_within,
                "pointwise_exceeding_count": len(point_rows) - point_within,
                "pointwise_total": len(point_rows),
                "pointwise_within_fraction": point_within / len(point_rows),
                "pointwise_exceeding_fraction": (len(point_rows) - point_within)
                / len(point_rows),
                "uniform_families_within_count": family_within,
                "families_exceeding_count": len(geometry_rows) - family_within,
                "family_total": len(geometry_rows),
            }
        )
    if (
        counts[0]["pointwise_within_count"] != 171
        or counts[0]["pointwise_exceeding_count"] != 0
        or counts[0]["uniform_families_within_count"] != 9
        or counts[3]["pointwise_within_count"] != 66
        or counts[3]["pointwise_exceeding_count"] != 105
        or counts[3]["uniform_families_within_count"] != 0
    ):
        raise RuntimeError("accepted AP-0/AP-1 classification counts changed")
    return counts


def build_sampled_pairwise_metrics(
    sources: Mapping[int, Mapping[str, Sequence[Mapping[str, Any]]]],
) -> list[dict[str, Any]]:
    spectra = {
        theta: _spectrum_lookup(spectra_as_numbers(source["spectra"]))
        for theta, source in sources.items()
    }
    rows: list[dict[str, Any]] = []
    for case in screening_geometries():
        for beta in beta_grid():
            arrays: dict[int | str, NDArray[np.float64]] = {}
            for theta, model_id in ((2, "T2"), (3, "T3"), (4, "T4"), (5, "T5")):
                selected = sorted(
                    spectra[theta][(case.geometry_id, float(beta), model_id)],
                    key=lambda row: int(row["mode"]),
                )
                arrays[theta] = np.asarray(
                    [float(row["lambda_ref"]) for row in selected], dtype=float
                )
            selected_t0 = sorted(
                spectra[2][(case.geometry_id, float(beta), "T0")],
                key=lambda row: int(row["mode"]),
            )
            arrays["T0"] = np.asarray(
                [float(row["lambda_ref"]) for row in selected_t0], dtype=float
            )
            gaps = {key: normalized_neighbor_gaps(value) for key, value in arrays.items()}
            for lower in range(COMPARED_ROOT_COUNT):
                delta = {
                    theta: float(gaps[theta][lower] - gaps["T0"][lower])
                    for theta in (2, 3, 4, 5)
                }
                signs = {
                    theta: int(np.sign(delta[theta])) for theta in (2, 3, 4, 5)
                }
                row: dict[str, Any] = {
                    "geometry_id": case.geometry_id,
                    "mu": case.mu,
                    "tau": case.tau,
                    "beta_deg": float(beta),
                    "pair_lower_mode": lower + 1,
                    "pair_upper_mode": lower + 2,
                    "g_pair_T0": float(gaps["T0"][lower]),
                }
                for theta in (2, 3, 4, 5):
                    row[f"g_pair_T{theta}"] = float(gaps[theta][lower])
                    row[f"delta_g_theta{theta}"] = delta[theta]
                    row[f"abs_delta_g_theta{theta}"] = abs(delta[theta])
                    row[f"sign_delta_g_theta{theta}"] = signs[theta]
                    row[f"pair_opened_at_theta{theta}"] = delta[theta] > 0.0
                for left, right in ((2, 3), (3, 4), (4, 5)):
                    change = float(gaps[right][lower] - gaps[left][lower])
                    row[f"change_g_{left}_to_{right}"] = change
                    row[f"abs_change_g_{left}_to_{right}"] = abs(change)
                    row[f"sign_reversal_{left}_to_{right}"] = (
                        signs[left] * signs[right] < 0
                    )
                rows.append(row)
    return rows


def sampled_pairwise_summary(
    pair_rows: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    baseline = np.asarray([float(row["g_pair_T0"]) for row in pair_rows])
    quartile = float(np.quantile(baseline, 0.25))
    low = baseline <= quartile

    def spearman(left: NDArray[np.float64], right: NDArray[np.float64]) -> float:
        return float(np.corrcoef(_rankdata(left), _rankdata(right))[0, 1])

    by_theta: dict[str, Any] = {}
    for theta in (2, 3, 4, 5):
        changes = np.asarray(
            [float(row[f"abs_delta_g_theta{theta}"]) for row in pair_rows]
        )
        opened = np.asarray(
            [bool(row[f"pair_opened_at_theta{theta}"]) for row in pair_rows]
        )
        per_pair: dict[str, float] = {}
        for lower in range(1, 7):
            selected = np.asarray(
                [int(row["pair_lower_mode"]) == lower for row in pair_rows]
            )
            per_pair[f"{lower}-{lower + 1}"] = spearman(
                baseline[selected], changes[selected]
            )
        by_theta[str(theta)] = {
            "pooled_spearman": spearman(baseline, changes),
            "per_pair_spearman": per_pair,
            "pooled_median_abs_delta_g": float(np.median(changes)),
            "lowest_gap_quartile_median_abs_delta_g": float(
                np.median(changes[low])
            ),
            "remaining_median_abs_delta_g": float(np.median(changes[~low])),
            "lowest_gap_quartile_opened_pair_fraction": float(
                np.mean(opened[low])
            ),
            "remaining_opened_pair_fraction": float(np.mean(opened[~low])),
        }
    increases = {}
    reversals = {}
    for left, right in ((2, 3), (3, 4), (4, 5)):
        increases[f"fraction_abs_delta_g_theta{right}_greater_theta{left}"] = float(
            np.mean(
                [
                    float(row[f"abs_delta_g_theta{right}"])
                    > float(row[f"abs_delta_g_theta{left}"])
                    for row in pair_rows
                ]
            )
        )
        reversals[f"sign_reversal_count_{left}_to_{right}"] = sum(
            bool(row[f"sign_reversal_{left}_to_{right}"]) for row in pair_rows
        )
    return {
        "scope": "descriptive pairwise tendency on the finite screening grid",
        "observation_count": len(pair_rows),
        "lowest_gap_quartile_upper_bound": quartile,
        "lowest_gap_quartile_count": int(np.count_nonzero(low)),
        "by_theta": by_theta,
        "interval_abs_change_increase_fractions": increases,
        "sign_reversal_counts": reversals,
        "inferential_p_values_computed": False,
        "independence_claimed": False,
        "regression_law_fitted": False,
        "causality_claimed": False,
    }


def _sampled_reuse_audit() -> list[dict[str, Any]]:
    directories = {
        2: AP0_OUTPUT_DIR,
        3: AP2_OUTPUT_DIRS[3],
        4: AP2_OUTPUT_DIRS[4],
        5: AP1_OUTPUT_DIR,
    }
    rows: list[dict[str, Any]] = []
    for theta, directory in directories.items():
        spectra = spectra_as_numbers(_read_csv(directory / SPECTRA_FILENAME))
        expected_models = MODEL_IDS if theta == 2 else (
            AP1_MODEL_IDS if theta == 5 else AP2_MODEL_IDS[theta]
        )
        validate_spectra(
            spectra,
            screening_geometries(),
            beta_grid(),
            model_ids=expected_models,
        )
        rows.append(
            {
                "theta_deg": float(theta),
                "source_directory": str(directory.relative_to(ROOT)),
                "spectra_source_file": str(
                    (directory / SPECTRA_FILENAME).relative_to(ROOT)
                ),
                "spectra_source_sha256": _sha256(directory / SPECTRA_FILENAME),
                "point_source_sha256": _sha256(directory / POINT_METRICS_FILENAME),
                "geometry_source_sha256": _sha256(
                    directory / GEOMETRY_SUMMARY_FILENAME
                ),
                "spectrum_row_count": len(spectra),
                "root_count": GUARD_ROOT_COUNT,
                "quality_validation_status": "PASS",
                "scientific_solver_called": False,
                "reuse_status": "PASS",
            }
        )
    return rows


def create_sampled_figures(
    geometry_rows: Sequence[Mapping[str, Any]],
    counts: Sequence[Mapping[str, Any]],
    output_dir: Path,
) -> list[str]:
    paths: list[str] = []
    with plt.rc_context(supervisor._style_context()):
        figure, axis = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
        colors = plt.get_cmap("tab10").colors
        for index, row in enumerate(geometry_rows):
            values = [float(row[f"Delta_simpl_max_theta{theta}"]) for theta in (2, 3, 4, 5)]
            axis.plot(
                (2, 3, 4, 5), values, marker="o", markersize=3.5,
                linewidth=1.2, color=colors[index],
                label=f"mu={float(row['mu']):g}, tau={float(row['tau']):g}",
            )
        axis.axhline(APPLICABILITY_THRESHOLD, color="#CC0000", linestyle="--", linewidth=1.1)
        axis.set_xticks((2, 3, 4, 5))
        axis.set_xlabel(r"sampled $\theta$, degrees")
        axis.set_ylabel(r"$\Delta_{\mathrm{simpl},\max}$")
        axis.grid(True, color="#D9D9D9", linewidth=0.5)
        axis.legend(frameon=False, ncol=3, fontsize=7)
        paths.extend(_save_figure(figure, output_dir, AP2_FIGURE_S3_BASENAME))
    with plt.rc_context(supervisor._style_context()):
        figure, axes = plt.subplots(1, 2, figsize=(8.0, 4.8), constrained_layout=True)
        theta = [int(float(row["theta_deg"])) for row in counts]
        point = [int(row["pointwise_exceeding_count"]) for row in counts]
        family = [int(row["families_exceeding_count"]) for row in counts]
        axes[0].bar(theta, point, color="#0072B2")
        axes[0].set_ylabel("exceeding configurations / 171")
        axes[0].set_ylim(0, 171)
        axes[1].bar(theta, family, color="#D55E00")
        axes[1].set_ylabel("exceeding geometry families / 9")
        axes[1].set_ylim(0, 9)
        for axis in axes:
            axis.set_xticks((2, 3, 4, 5))
            axis.set_xlabel(r"sampled $\theta$, degrees")
            axis.grid(True, axis="y", color="#D9D9D9", linewidth=0.5)
        paths.extend(_save_figure(figure, output_dir, AP2_FIGURE_S4_BASENAME))
    return paths


def _sampled_expected_output_paths(output_dir: Path) -> list[Path]:
    names = (
        SAMPLED_POINT_FILENAME,
        SAMPLED_GEOMETRY_FILENAME,
        SAMPLED_COUNTS_FILENAME,
        SAMPLED_PAIRWISE_FILENAME,
        SAMPLED_REUSE_AUDIT_FILENAME,
        SAMPLED_SUMMARY_FILENAME,
        SAMPLED_REPORT_FILENAME,
        f"{AP2_FIGURE_S3_BASENAME}.pdf",
        f"{AP2_FIGURE_S3_BASENAME}.png",
        f"{AP2_FIGURE_S4_BASENAME}.pdf",
        f"{AP2_FIGURE_S4_BASENAME}.png",
    )
    return [output_dir / name for name in names]


def sampled_theta_report(summary: Mapping[str, Any]) -> str:
    counts = summary["classification_counts"]
    dynamics = summary["sampled_dynamics"]
    lines = [
        "# AP-2 sampled-theta comparison: 2, 3, 4, 5 degrees",
        "",
        f"**AP-2 computational workflow: {summary['workflow_status']}**",
        "",
        "Theta is the sampled material-axis orientation; beta is the geometric joint angle. The connected figure points visualize sampled values only and do not establish a continuous threshold angle.",
        "",
        "## Classification counts",
        "",
        "| theta | pointwise within/171 | pointwise exceeds/171 | families within/9 | families exceeds/9 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in counts:
        lines.append(
            f"| {float(row['theta_deg']):g} | {row['pointwise_within_count']} | {row['pointwise_exceeding_count']} | {row['uniform_families_within_count']} | {row['families_exceeding_count']} |"
        )
    lines.extend(
        [
            "",
            "## Sampled dynamics",
            "",
            f"First pointwise exceedance counts at theta=3/4/5: `{dynamics['first_exceeding_point_counts']['3']}/{dynamics['first_exceeding_point_counts']['4']}/{dynamics['first_exceeding_point_counts']['5']}`; within at all sampled theta: `{dynamics['points_within_all_sampled']}`.",
            f"Pointwise sampled nondecreasing/nonmonotone: `{dynamics['pointwise_sampled_nondecreasing']}/{dynamics['pointwise_sampled_nonmonotone']}`; re-entry after exceedance: `{dynamics['pointwise_reentry_count']}`.",
            f"Family sampled nondecreasing/nonmonotone: `{dynamics['family_sampled_nondecreasing']}/{dynamics['family_sampled_nonmonotone']}`; family re-entry: `{dynamics['family_reentry_count']}`.",
            f"Pointwise maximizing-position changes: `{dynamics['pointwise_maximizing_sorted_position_changes']}`; family maximizing-position changes: `{dynamics['family_maximizing_sorted_position_changes']}`; family beta-at-maximum changes: `{dynamics['family_beta_at_max_changes']}`.",
            f"Pointwise increment statistics: `{dynamics['pointwise_increment_statistics']}`.",
            f"Family-maximum increment statistics: `{dynamics['family_maximum_increment_statistics']}`.",
            "",
            "The first exceeding sampled theta is not a critical angle. No interpolation or monotonicity claim between sampled orientations is made. Equal sorted indices do not establish physical modal-descendant identity.",
            "",
            "## Geometry-family sampled dynamics",
            "",
            "| mu | tau | Delta simpl max at 2/3/4/5 deg | status sequence | first exceeding sampled theta | beta/k at 2/3/4/5 deg |",
            "|---:|---:|---|---|---:|---|",
        ]
    )
    for row in summary["geometry_table"]:
        delta = "/".join(
            f"{float(row[f'Delta_simpl_max_theta{theta}']):.6f}"
            for theta in (2, 3, 4, 5)
        )
        beta_mode = ", ".join(
            f"{float(row[f'beta_at_simpl_max_theta{theta}']):g}/{int(row[f'k_at_simpl_max_theta{theta}'])}"
            for theta in (2, 3, 4, 5)
        )
        first = row["first_exceeding_sampled_theta_deg"] or "none"
        lines.append(
            f"| {float(row['mu']):.2f} | {float(row['tau']):+.2f} | {delta} | {row['sampled_theta_family_status_sequence']} | {first} | {beta_mode} |"
        )
    lines.extend(
        [
            "",
            "## Classification decomposition",
            "",
            f"A/B/C at theta 3, 4, 5: `{summary['classification_decomposition']}`. This is a classification decomposition, not an additive error budget.",
            "",
            "## Pairwise gaps",
            "",
            f"Full descriptive summary: `{summary['pairwise_descriptive_summary']}`.",
            "The pairwise statistics describe the same adjacent sorted pair on the finite grid. No p-values, observation-independence claim, regression law, or causal claim is made.",
            "",
            "## Preservation",
            "",
            f"AP-0/AP-1/supervisor unchanged: `{summary['ap0_outputs_unchanged']}/{summary['ap1_outputs_unchanged']}/{summary['supervisor_outputs_unchanged']}`.",
            "",
            "## Commands",
            "",
            "```text",
            "python scripts/analysis/anisotropic_rods/screen_yartsev_ch2_spectral_applicability.py --combine-sampled-theta",
            "```",
        ]
    )
    return "\n".join(lines) + "\n"


def run_sampled_theta_comparison(args: argparse.Namespace) -> int:
    if args.smoke or args.force_recompute or args.resume:
        raise ValueError(
            "--combine-sampled-theta is solver-free and does not accept smoke, resume, or force-recompute"
        )
    uses_default = args.output_dir.resolve() == DEFAULT_OUTPUT_DIR.resolve()
    output_dir = _safe_screening_output(
        SAMPLED_THETA_OUTPUT_DIR if uses_default else args.output_dir
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    preservation_before = accepted_artifact_hashes()
    try:
        sources = _sampled_source_rows()
        point_rows = build_sampled_point_comparison(sources)
        geometry_rows = build_sampled_geometry_comparison(sources)
        count_rows = build_sampled_classification_counts(point_rows, geometry_rows)
        pair_rows = build_sampled_pairwise_metrics(sources)
        reuse_audit = _sampled_reuse_audit()
        _write_csv(output_dir / SAMPLED_POINT_FILENAME, point_rows)
        _write_csv(output_dir / SAMPLED_GEOMETRY_FILENAME, geometry_rows)
        _write_csv(output_dir / SAMPLED_COUNTS_FILENAME, count_rows)
        _write_csv(output_dir / SAMPLED_PAIRWISE_FILENAME, pair_rows)
        _write_csv(output_dir / SAMPLED_REUSE_AUDIT_FILENAME, reuse_audit)
        figure_paths = create_sampled_figures(geometry_rows, count_rows, output_dir)
        pairwise = sampled_pairwise_summary(pair_rows)
        first_counts = {
            str(theta): sum(
                row["first_exceeding_sampled_theta_deg"] == float(theta)
                for row in point_rows
            )
            for theta in (3, 4, 5)
        }
        dynamics = {
            "first_exceeding_point_counts": first_counts,
            "points_within_all_sampled": sum(
                row["first_exceeding_sampled_theta_deg"] == ""
                for row in point_rows
            ),
            "pointwise_sampled_nondecreasing": sum(
                bool(row["delta_simpl_sampled_nondecreasing"])
                for row in point_rows
            ),
            "pointwise_sampled_nonmonotone": sum(
                not bool(row["delta_simpl_sampled_nondecreasing"])
                for row in point_rows
            ),
            "pointwise_reentry_count": sum(
                bool(row["has_reentry_after_exceedance"]) for row in point_rows
            ),
            "family_sampled_nondecreasing": sum(
                bool(row["sampled_nondecreasing"]) for row in geometry_rows
            ),
            "family_sampled_nonmonotone": sum(
                not bool(row["sampled_nondecreasing"]) for row in geometry_rows
            ),
            "family_reentry_count": sum(
                bool(row["has_reentry_after_exceedance"])
                for row in geometry_rows
            ),
            "pointwise_maximizing_sorted_position_changes": sum(
                int(row[f"argmax_k_simpl_theta{left}"])
                != int(row[f"argmax_k_simpl_theta{right}"])
                for row in point_rows
                for left, right in ((2, 3), (3, 4), (4, 5))
            ),
            "family_maximizing_sorted_position_changes": sum(
                int(row[f"k_at_simpl_max_theta{left}"])
                != int(row[f"k_at_simpl_max_theta{right}"])
                for row in geometry_rows
                for left, right in ((2, 3), (3, 4), (4, 5))
            ),
            "family_beta_at_max_changes": sum(
                float(row[f"beta_at_simpl_max_theta{left}"])
                != float(row[f"beta_at_simpl_max_theta{right}"])
                for row in geometry_rows
                for left, right in ((2, 3), (3, 4), (4, 5))
            ),
            "pointwise_increment_statistics": {
                f"{left}_to_{right}": {
                    "minimum": float(
                        np.min(
                            [
                                float(row[f"increment_simpl_{left}_to_{right}"])
                                for row in point_rows
                            ]
                        )
                    ),
                    "median": float(
                        np.median(
                            [
                                float(row[f"increment_simpl_{left}_to_{right}"])
                                for row in point_rows
                            ]
                        )
                    ),
                    "maximum": float(
                        np.max(
                            [
                                float(row[f"increment_simpl_{left}_to_{right}"])
                                for row in point_rows
                            ]
                        )
                    ),
                }
                for left, right in ((2, 3), (3, 4), (4, 5))
            },
            "family_maximum_increment_statistics": {
                f"{left}_to_{right}": {
                    "minimum": float(
                        np.min(
                            [
                                float(row[f"increment_simpl_max_{left}_to_{right}"])
                                for row in geometry_rows
                            ]
                        )
                    ),
                    "median": float(
                        np.median(
                            [
                                float(row[f"increment_simpl_max_{left}_to_{right}"])
                                for row in geometry_rows
                            ]
                        )
                    ),
                    "maximum": float(
                        np.max(
                            [
                                float(row[f"increment_simpl_max_{left}_to_{right}"])
                                for row in geometry_rows
                            ]
                        )
                    ),
                }
                for left, right in ((2, 3), (3, 4), (4, 5))
            },
        }
        decomposition = {
            str(theta): {
                "A_orientation_effect_exceeds": sum(
                    float(row[f"Delta_orient_theta{theta}"])
                    > APPLICABILITY_THRESHOLD
                    for row in point_rows
                ),
                "B_orientation_within_but_simplification_exceeds": sum(
                    float(row[f"Delta_orient_theta{theta}"])
                    <= APPLICABILITY_THRESHOLD
                    and float(row[f"Delta_simpl_theta{theta}"])
                    > APPLICABILITY_THRESHOLD
                    for row in point_rows
                ),
                "C_simplification_within": sum(
                    float(row[f"Delta_simpl_theta{theta}"])
                    <= APPLICABILITY_THRESHOLD
                    for row in point_rows
                ),
            }
            for theta in (3, 4, 5)
        }
        preservation_after = accepted_artifact_hashes()
        summary: dict[str, Any] = {
            "workflow": "AP-2 intermediate-angle sampled-theta comparison",
            "workflow_status": "PENDING",
            "execution_mode": "reuse_only_postprocessing",
            "scientific_solver_called": False,
            "output_directory": str(output_dir.relative_to(ROOT)),
            "git_context": git_context(),
            "row_counts": {
                SAMPLED_POINT_FILENAME: len(point_rows),
                SAMPLED_GEOMETRY_FILENAME: len(geometry_rows),
                SAMPLED_COUNTS_FILENAME: len(count_rows),
                SAMPLED_PAIRWISE_FILENAME: len(pair_rows),
                SAMPLED_REUSE_AUDIT_FILENAME: len(reuse_audit),
            },
            "classification_counts": count_rows,
            "geometry_table": geometry_rows,
            "sampled_dynamics": dynamics,
            "classification_decomposition": decomposition,
            "pairwise_descriptive_summary": pairwise,
            "reuse_audit": reuse_audit,
            "reuse_audit_pass": all(row["reuse_status"] == "PASS" for row in reuse_audit),
            "sampled_nondecreasing_absolute_tolerance": SAMPLED_NONDECREASING_ABS_TOLERANCE,
            "figure_paths": figure_paths,
            "runtimes_seconds": {
                "scientific_current": 0.0,
                "workflow_wall_current": time.perf_counter() - started,
            },
            "preservation_sha256_before": preservation_before,
            "preservation_sha256_after": preservation_after,
            "ap0_outputs_unchanged": preservation_before["ap0"] == preservation_after["ap0"],
            "ap1_outputs_unchanged": preservation_before["ap1"] == preservation_after["ap1"],
            "supervisor_outputs_unchanged": preservation_before["supervisor"] == preservation_after["supervisor"],
            "limitations": {
                "first_exceeding_sampled_theta_is_critical_angle": False,
                "continuous_monotonicity_established": False,
                "interpolation_performed": False,
                "additional_theta_sweep_performed": False,
                "modal_descendants_tracked": False,
            },
        }
        _json_write(output_dir / SAMPLED_SUMMARY_FILENAME, summary)
        (output_dir / SAMPLED_REPORT_FILENAME).write_text(
            sampled_theta_report(summary), encoding="utf-8"
        )
        gates = (
            len(point_rows) == 171,
            len(geometry_rows) == 9,
            len(count_rows) == 4,
            len(pair_rows) == 1026,
            len(reuse_audit) == 4,
            summary["reuse_audit_pass"],
            summary["ap0_outputs_unchanged"],
            summary["ap1_outputs_unchanged"],
            summary["supervisor_outputs_unchanged"],
            all(path.is_file() for path in _sampled_expected_output_paths(output_dir)),
        )
        summary["workflow_status"] = "PASS" if all(gates) else "FAIL"
        summary["runtimes_seconds"]["workflow_wall_current"] = time.perf_counter() - started
        _json_write(output_dir / SAMPLED_SUMMARY_FILENAME, summary)
        (output_dir / SAMPLED_REPORT_FILENAME).write_text(
            sampled_theta_report(summary), encoding="utf-8"
        )
    except Exception as error:
        failure = {
            "workflow": "AP-2 intermediate-angle sampled-theta comparison",
            "workflow_status": "FAIL",
            "error": str(error),
            "scientific_solver_called": False,
            "git_context": git_context(),
        }
        _json_write(output_dir / SAMPLED_SUMMARY_FILENAME, failure)
        print("AP-2 combined sampled-theta comparison: FAIL", file=sys.stderr)
        print(str(error), file=sys.stderr)
        return 1
    print(f"output_dir={output_dir}")
    print(f"sampled_point_row_count={len(point_rows)}")
    print(f"sampled_geometry_row_count={len(geometry_rows)}")
    print(f"sampled_pairwise_row_count={len(pair_rows)}")
    print(f"AP-2 combined sampled-theta comparison: {summary['workflow_status']}")
    return 0 if summary["workflow_status"] == "PASS" else 1


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.combine_sampled_theta:
        return run_sampled_theta_comparison(args)
    configuration = screening_configuration(args.target_theta_deg)
    if configuration.is_ap1:
        return run_ap1(args, configuration)
    if configuration.is_ap2:
        return run_ap2(args, configuration)
    requested_output = args.output_dir
    if args.smoke:
        requested_output = (
            DEFAULT_SMOKE_OUTPUT_DIR
            if args.output_dir.resolve() == DEFAULT_OUTPUT_DIR.resolve()
            else args.output_dir / "smoke"
        )
    output_dir = _safe_screening_output(requested_output)
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    supervisor_hashes_before = supervisor_data_hashes()
    previous_summary_path = output_dir / SUMMARY_FILENAME
    previous_summary = (
        json.loads(previous_summary_path.read_text(encoding="utf-8"))
        if previous_summary_path.is_file()
        else {}
    )
    scientific_runtime = 0.0
    family_metadata: dict[str, dict[str, Any]] = {}
    try:
        if args.reuse_data:
            required = (
                output_dir / GEOMETRY_FILENAME,
                output_dir / SPECTRA_FILENAME,
                output_dir / REUSE_AUDIT_FILENAME,
            )
            missing = [str(path) for path in required if not path.is_file()]
            if missing:
                raise FileNotFoundError(
                    "--reuse-data requires saved AP-0 CSV: " + ", ".join(missing)
                )
            geometries = _geometry_rows_as_numbers(
                _read_csv(output_dir / GEOMETRY_FILENAME)
            )
            beta_values = beta_grid(smoke=args.smoke)
            spectra = spectra_as_numbers(_read_csv(output_dir / SPECTRA_FILENAME))
            reuse_audit: list[dict[str, Any]] = [
                dict(row) for row in _read_csv(output_dir / REUSE_AUDIT_FILENAME)
            ]
            reuse_audit = reuse_audit_as_numbers(reuse_audit)
            execution_mode = "reuse_data"
        else:
            supervisor._require_fast_solver_pass(SUPERVISOR_OUTPUT_DIR)
            geometries = screening_geometries(smoke=args.smoke)
            beta_values = beta_grid(smoke=args.smoke)
            _write_csv(
                output_dir / GEOMETRY_FILENAME,
                [asdict(case) for case in geometries],
            )
            aggregate = PerformanceCounters()
            transfer_cache = ExactTransferLRU(
                FastSweepSettings().cache_maxsize, counters=aggregate
            )
            spectra = []
            reuse_audit = []
            for case in geometries:
                for model_id in MODEL_IDS:
                    reuse = _reuse_spec(case, model_id)
                    if reuse is not None:
                        reused_rows, audit = load_reused_family(
                            case, model_id, beta_values
                        )
                        spectra.extend(reused_rows)
                        reuse_audit.append(audit)
                        continue
                    rows, metadata, runtime = solve_new_family(
                        case,
                        model_id,
                        beta_values,
                        output_dir=output_dir,
                        resume=args.resume and not args.force_recompute,
                        transfer_cache=transfer_cache,
                        counters=aggregate,
                        smoke=args.smoke,
                    )
                    spectra.extend(rows)
                    family_metadata[metadata["family_id"]] = metadata
                    scientific_runtime += runtime
            spectra.sort(
                key=lambda row: (
                    str(row["geometry_id"]),
                    float(row["beta_deg"]),
                    MODEL_IDS.index(str(row["model_id"])),
                    int(row["mode"]),
                )
            )
            _write_csv(output_dir / SPECTRA_FILENAME, spectra)
            _write_csv(output_dir / REUSE_AUDIT_FILENAME, reuse_audit)
            execution_mode = (
                "force_recompute" if args.force_recompute else "scientific"
            )

        validate_geometry_manifest(geometries, smoke=args.smoke)
        validate_spectra(spectra, geometries, beta_values)
        point_rows = build_point_metrics(spectra, geometries, beta_values)
        summaries = build_geometry_summary(point_rows, geometries)
        candidates = select_candidates(summaries)
        _write_csv(output_dir / POINT_METRICS_FILENAME, point_rows)
        _write_csv(output_dir / GEOMETRY_SUMMARY_FILENAME, summaries)
        _json_write(output_dir / CANDIDATES_FILENAME, candidates)
        figure_paths = create_diagnostic_figures(
            summaries, point_rows, output_dir
        )
        supervisor_hashes_after = supervisor_data_hashes()
        wall_runtime = time.perf_counter() - started
        summary = _build_summary(
            smoke=args.smoke,
            execution_mode=execution_mode,
            geometries=geometries,
            beta_values=beta_values,
            spectra=spectra,
            point_rows=point_rows,
            summaries=summaries,
            reuse_audit=reuse_audit,
            candidates=candidates,
            family_metadata=family_metadata,
            previous_summary=previous_summary,
            supervisor_hashes_before=supervisor_hashes_before,
            supervisor_hashes_after=supervisor_hashes_after,
            scientific_runtime=scientific_runtime,
            wall_runtime=wall_runtime,
            figure_paths=figure_paths,
            output_dir=output_dir,
        )
        # Write provisional text so the final output-existence gate can include
        # the report and summary themselves, then replace PENDING atomically at
        # the file level with the final status.
        _json_write(output_dir / SUMMARY_FILENAME, summary)
        (output_dir / REPORT_FILENAME).write_text(
            screening_report(summary), encoding="utf-8"
        )
        status = determine_workflow_status(
            smoke=args.smoke,
            geometries=geometries,
            beta_values=beta_values,
            spectra=spectra,
            point_rows=point_rows,
            summaries=summaries,
            reuse_audit=reuse_audit,
            supervisor_hashes_unchanged=(
                supervisor_hashes_before == supervisor_hashes_after
            ),
            output_dir=output_dir,
        )
        summary["workflow_status"] = status
        summary["runtimes_seconds"]["workflow_wall_current"] = (
            time.perf_counter() - started
        )
        _json_write(output_dir / SUMMARY_FILENAME, summary)
        (output_dir / REPORT_FILENAME).write_text(
            screening_report(summary), encoding="utf-8"
        )
    except Exception as error:
        failure = {
            "workflow": "AP-0 small-grid spectral-applicability screening",
            "workflow_status": "FAIL",
            "smoke": bool(args.smoke),
            "error": str(error),
            "git_context": git_context(),
        }
        _json_write(output_dir / SUMMARY_FILENAME, failure)
        print("AP-0 small-grid screening workflow: FAIL", file=sys.stderr)
        print(str(error), file=sys.stderr)
        return 1

    label = "AP-0 smoke workflow" if args.smoke else "AP-0 small-grid screening workflow"
    print(f"output_dir={output_dir}")
    print(f"geometry_count={len(geometries)}")
    print(f"beta_point_count={len(beta_values)}")
    print(f"spectrum_row_count={len(spectra)}")
    print(f"{label}: {summary['workflow_status']}")
    if not args.smoke:
        print(
            "Number of geometries within 10%: "
            f"{summary['scientific_result']['within_10_percent_count']} / 9"
        )
        print(
            "Number of geometries exceeding 10%: "
            f"{summary['scientific_result']['exceeding_10_percent_count']} / 9"
        )
    return 0 if summary["workflow_status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
