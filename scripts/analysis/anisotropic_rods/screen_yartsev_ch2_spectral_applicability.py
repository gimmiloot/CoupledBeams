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
SUPERVISOR_OUTPUT_DIR = supervisor.DEFAULT_OUTPUT_DIR

MU_VALUES = (0.0, 0.25, 0.5)
TAU_VALUES = (-0.2, 0.0, 0.2)
MODEL_IDS = ("T2", "T0", "EB0")
MODEL_THETA_DEG = {"T2": 2.0, "T0": 0.0, "EB0": 0.0}
MODEL_SOLVER_PATH = {
    "T2": "book_slope_clamp",
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
    if model_id not in MODEL_IDS:
        raise ValueError(f"unsupported AP-0 model: {model_id}")
    theta = MODEL_THETA_DEG[model_id]
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


def normalization_contract(case: GeometryCase, model_id: str) -> dict[str, Any]:
    return {
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
            "theta_deg": MODEL_THETA_DEG[model_id],
            "mu": case.mu,
            "tau": case.tau,
            "clamp": CLAMP,
            "root_count": GUARD_ROOT_COUNT,
            "quality_determinant_tolerance": ROOT_DETERMINANT_TOLERANCE,
            "quality_singular_tolerance": ROOT_SINGULAR_TOLERANCE,
        },
    }


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
        "theta_deg": MODEL_THETA_DEG[model_id],
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


def _family_id(case: GeometryCase, model_id: str, *, smoke: bool) -> str:
    prefix = "ap0_smoke" if smoke else "ap0"
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
) -> tuple[list[dict[str, Any]], dict[str, Any], float]:
    preset = model_preset(case, model_id)
    solver_path = MODEL_SOLVER_PATH[model_id]
    normalization = normalization_contract(case, model_id)
    family_id = _family_id(case, model_id, smoke=smoke)
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
        "theta_deg": MODEL_THETA_DEG[model_id],
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
) -> None:
    expected_count = (
        len(geometries) * len(beta_values) * len(MODEL_IDS) * GUARD_ROOT_COUNT
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
    expected_groups = len(geometries) * len(beta_values) * len(MODEL_IDS)
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


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
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
