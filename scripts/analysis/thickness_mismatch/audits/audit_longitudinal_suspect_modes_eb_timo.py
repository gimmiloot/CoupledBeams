from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    assemble_clamped_coupled_matrix_eta,
)
from scripts.analysis.thickness_mismatch.maps import (  # noqa: E402
    plot_eb_vs_timoshenko_lambda_beta_cases as beta_workflow,
)
from scripts.lib import in_plane_shape_geometry as DISPLAY  # noqa: E402
from scripts.lib import variable_length_timoshenko as TIMO  # noqa: E402
from scripts.lib.analytic_coupled_rods_shapes import analytic_null_vector  # noqa: E402


MODEL_EB = beta_workflow.MODEL_EB
MODEL_TIMO = beta_workflow.MODEL_TIMO
MODELS = (MODEL_EB, MODEL_TIMO)

DEFAULT_BETA_DEG = 45.0
DEFAULT_ETA = 0.0
DEFAULT_SUSPECTS = ("0.03:5", "0.05:4")
DEFAULT_MU_VALUES = (0.0, 0.35, 0.7)
DEFAULT_OUTPUT_DIR = Path("results") / "eb_vs_timoshenko_longitudinal_suspect_modes"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "eb_vs_timoshenko_longitudinal_suspect_modes"
DEFAULT_LAMBDA_CSV_DIRS = (
    Path("results") / "eb_vs_timoshenko_lambda_mu_beta45_eta0_eps_scan",
    Path("results") / "eb_vs_timoshenko_lambda_mu_cases",
)
AUTO_GRID_STEP = 0.01
DEFAULT_N_POINTS = 801
SMOKE_N_POINTS = 201
DEFAULT_N_ROOTS = 6
MODE_SCALE = 0.16
ENERGY_SUM_TOL = 5.0e-4
JOINT_KINEMATIC_TOL = 1.0e-6
JOINT_FORCE_TOL = 1.0e-6

ENERGY_FIELDS = [
    "epsilon",
    "beta_deg",
    "eta",
    "mu",
    "model",
    "sorted_index",
    "Lambda",
    "U_axial",
    "U_bending",
    "U_shear",
    "U_total",
    "axial_energy_fraction",
    "bending_energy_fraction",
    "shear_energy_fraction",
    "bending_shear_fraction",
    "axial_displacement_fraction",
    "transverse_displacement_fraction",
    "classification",
    "notes",
]

SHAPE_SUMMARY_FIELDS = [
    "epsilon",
    "mu",
    "model",
    "sorted_index",
    "Lambda",
    "figure_file",
    "sign_normalization",
    "warnings",
    "notes",
]

TIMO_JOINT_CONTINUITY_FIELDS = [
    "epsilon",
    "mu",
    "beta_deg",
    "eta",
    "model",
    "sorted_index",
    "Lambda",
    "gap_w",
    "gap_u",
    "gap_psi",
    "gap_M",
    "gap_Q",
    "gap_N",
    "max_abs_kinematic_gap",
    "max_abs_force_gap",
    "pass_kinematic",
    "pass_force",
    "notes",
]


@dataclass(frozen=True)
class SuspectSpec:
    epsilon: float
    sorted_index: int


@dataclass(frozen=True)
class RootPair:
    eb: float
    timo: float
    source: str
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class AutoPoint:
    mu: float
    reason: str
    source: str


@dataclass(frozen=True)
class PlotComponents:
    u_left: np.ndarray
    w_left: np.ndarray
    u_right: np.ndarray
    w_right: np.ndarray


@dataclass(frozen=True)
class ModeDiagnostic:
    epsilon: float
    beta_deg: float
    eta: float
    mu: float
    model: str
    sorted_index: int
    Lambda: float
    plot_components: PlotComponents
    energy_row: dict[str, object]
    shape_row: dict[str, object]
    warnings: tuple[str, ...]
    figure_file: Path | None = None
    joint_row: dict[str, object] | None = None


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def fmt(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        value_f = float(value)
        if not isfinite(value_f):
            return "nan"
        return f"{value_f:.16e}"
    return value


def write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})
    return path


def token(value: float, *, min_decimals: int = 2, max_decimals: int = 5) -> str:
    text = f"{float(value):.{max_decimals}f}".rstrip("0").rstrip(".")
    if "." not in text:
        if min_decimals > 0:
            text += "." + "0" * min_decimals
    else:
        decimals = len(text.split(".", 1)[1])
        if decimals < min_decimals:
            text += "0" * (min_decimals - decimals)
    return text.replace("-", "m").replace(".", "p")


def eps_token(value: float) -> str:
    return token(float(value), min_decimals=2, max_decimals=5)


def mu_token(value: float) -> str:
    return token(float(value), min_decimals=2, max_decimals=5)


def number_text(value: object, digits: int = 6) -> str:
    try:
        value_f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not isfinite(value_f):
        return "nan"
    return f"{value_f:.{digits}g}"


def parse_suspect(value: str) -> SuspectSpec:
    parts = str(value).split(":", 1)
    if len(parts) != 2:
        raise argparse.ArgumentTypeError("suspect must look like epsilon:sorted_index")
    try:
        epsilon = float(parts[0])
        sorted_index = int(parts[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid suspect specification: {value!r}") from exc
    if epsilon <= 0.0 or sorted_index <= 0:
        raise argparse.ArgumentTypeError("epsilon and sorted_index must be positive")
    return SuspectSpec(epsilon=epsilon, sorted_index=sorted_index)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Diagnostic-only audit of suspicious sorted Lambda(mu) EB/Timoshenko "
            "branches using analytic mode shapes and energy fractions."
        ),
    )
    parser.add_argument("--beta-deg", type=float, default=DEFAULT_BETA_DEG)
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--suspect", type=parse_suspect, nargs="+", default=[parse_suspect(v) for v in DEFAULT_SUSPECTS])
    parser.add_argument("--mu-values", type=float, nargs="+", default=list(DEFAULT_MU_VALUES))
    parser.add_argument("--include-auto-mu-points", dest="include_auto_mu_points", action="store_true", default=True)
    parser.add_argument("--no-include-auto-mu-points", dest="include_auto_mu_points", action="store_false")
    parser.add_argument("--include-control-modes", dest="include_control_modes", action="store_true", default=True)
    parser.add_argument("--no-include-control-modes", dest="include_control_modes", action="store_false")
    parser.add_argument("--n-points", type=int, default=DEFAULT_N_POINTS)
    parser.add_argument("--n-roots", type=int, default=DEFAULT_N_ROOTS)
    parser.add_argument("--lambda-csv-dir", type=Path, action="append", default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--compile-tex", dest="compile_tex", action="store_true", default=True)
    parser.add_argument("--no-compile-tex", dest="compile_tex", action="store_false")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    if bool(args.smoke):
        args.suspect = [SuspectSpec(epsilon=0.05, sorted_index=4)]
        args.mu_values = [0.35]
        args.include_auto_mu_points = False
        args.include_control_modes = False
        args.n_points = SMOKE_N_POINTS
        args.output_dir = SMOKE_OUTPUT_DIR
        args.compile_tex = False

    args.output_dir = repo_path(Path(args.output_dir))
    user_dirs = tuple(repo_path(Path(path)) for path in (args.lambda_csv_dir or ()))
    args.lambda_csv_dirs = tuple(repo_path(path) for path in DEFAULT_LAMBDA_CSV_DIRS) + user_dirs
    args.suspect = tuple(args.suspect)
    args.mu_values = tuple(float(value) for value in args.mu_values)
    validate_args(args)
    return args


def validate_args(args: argparse.Namespace) -> None:
    if not (-1.0 < float(args.eta) < 1.0):
        raise ValueError("--eta must lie inside (-1, 1).")
    if int(args.n_points) < 51:
        raise ValueError("--n-points must be at least 51.")
    if int(args.n_roots) < max(spec.sorted_index for spec in args.suspect):
        raise ValueError("--n-roots must include every suspect sorted index.")
    if any(not (0.0 <= float(mu) <= 0.7) for mu in args.mu_values):
        raise ValueError("--mu-values must lie in [0, 0.7] for this audit.")
    for spec in args.suspect:
        TIMO.tau_factors(float(args.mu_values[0]), float(args.eta))
        if spec.epsilon <= 0.0:
            raise ValueError("epsilon values must be positive.")


def regular_grid(start: float, end: float, step: float) -> np.ndarray:
    values = np.arange(float(start), float(end) + 0.5 * float(step), float(step), dtype=float)
    values = values[values <= float(end) + 1.0e-12]
    if values.size == 0 or not np.isclose(values[0], float(start), rtol=0.0, atol=1.0e-12):
        values = np.insert(values, 0, float(start))
    if not np.isclose(values[-1], float(end), rtol=0.0, atol=1.0e-12):
        values = np.append(values, float(end))
    return np.unique(np.round(values, 12))


def lambda_csv_path(csv_dir: Path, beta_deg: float, eta: float, epsilon: float) -> Path:
    return Path(csv_dir) / (
        f"lambda_mu_eb_vs_timo_beta{token(float(beta_deg), min_decimals=0, max_decimals=5)}_"
        f"eta{token(float(eta), min_decimals=0, max_decimals=5)}_eps{eps_token(float(epsilon))}.csv"
    )


def read_lambda_rows(
    csv_dirs: Sequence[Path],
    *,
    beta_deg: float,
    eta: float,
    epsilon: float,
) -> tuple[list[dict[str, str]], Path | None]:
    for csv_dir in csv_dirs:
        path = lambda_csv_path(Path(csv_dir), beta_deg, eta, epsilon)
        if not path.exists():
            continue
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle)), path
    return [], None


def row_float(row: dict[str, str], *names: str) -> float:
    for name in names:
        if name not in row:
            continue
        try:
            value = float(row[name])
        except (TypeError, ValueError):
            continue
        if isfinite(value):
            return value
    return float("nan")


def rows_at_mu(rows: Sequence[dict[str, str]], mu: float, *, tol: float = 5.0e-10) -> list[dict[str, str]]:
    out = []
    for row in rows:
        try:
            row_mu = float(row["mu"])
        except (KeyError, ValueError):
            continue
        if abs(row_mu - float(mu)) <= tol:
            out.append(row)
    return out


def lambda_from_rows(rows: Sequence[dict[str, str]], *, mu: float, sorted_index: int) -> RootPair | None:
    selected = [
        row
        for row in rows_at_mu(rows, float(mu))
        if int(float(row.get("sorted_index", "nan"))) == int(sorted_index)
    ]
    if not selected:
        return None
    row = selected[0]
    eb = row_float(row, "Lambda_EB_plot", "Lambda_EB_raw")
    timo = row_float(row, "Lambda_Timoshenko_plot", "Lambda_Timoshenko_raw")
    if not (isfinite(eb) and isfinite(timo)):
        return None
    warnings = tuple(
        item
        for item in (
            row.get("root_warning_EB", ""),
            row.get("root_warning_Timoshenko", ""),
            row.get("notes_EB", ""),
            row.get("notes_Timoshenko", ""),
        )
        if item and item.lower() != "ok"
    )
    return RootPair(eb=eb, timo=timo, source="lambda_csv", warnings=warnings)


def solve_roots_at_point(
    *,
    beta_deg: float,
    eta: float,
    epsilon: float,
    mu: float,
    n_roots: int,
) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
    case = beta_workflow.CaseSpec(mu=float(mu), eta=float(eta), epsilon=float(epsilon))
    warnings: list[str] = []

    eb_result = beta_workflow.solve_model(case, float(beta_deg), int(n_roots), MODEL_EB)
    if eb_result.root_count_found < int(n_roots):
        eb_result = beta_workflow.retry_missing_roots(eb_result, case, float(beta_deg), int(n_roots), MODEL_EB)
    warnings.extend(f"EB: {item}" for item in eb_result.warnings)

    finite_eb = [float(root) for root in eb_result.roots if isfinite(float(root))]
    upper_hint = max(finite_eb) if finite_eb else None
    timo_result = beta_workflow.solve_model(
        case,
        float(beta_deg),
        int(n_roots),
        MODEL_TIMO,
        upper_hint=upper_hint,
    )
    if timo_result.root_count_found < int(n_roots):
        timo_result = beta_workflow.retry_missing_roots(
            timo_result,
            case,
            float(beta_deg),
            int(n_roots),
            MODEL_TIMO,
            upper_hint=upper_hint,
        )
    warnings.extend(f"Timoshenko: {item}" for item in timo_result.warnings)
    return np.asarray(eb_result.roots, dtype=float), np.asarray(timo_result.roots, dtype=float), tuple(warnings)


class RootProvider:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.csv_rows: dict[float, tuple[list[dict[str, str]], Path | None]] = {}
        self.root_cache: dict[tuple[float, float], tuple[np.ndarray, np.ndarray, tuple[str, ...]]] = {}

    def rows_for_epsilon(self, epsilon: float) -> tuple[list[dict[str, str]], Path | None]:
        key = round(float(epsilon), 12)
        if key not in self.csv_rows:
            self.csv_rows[key] = read_lambda_rows(
                self.args.lambda_csv_dirs,
                beta_deg=float(self.args.beta_deg),
                eta=float(self.args.eta),
                epsilon=float(epsilon),
            )
        return self.csv_rows[key]

    def roots(self, epsilon: float, mu: float) -> tuple[np.ndarray, np.ndarray, str, tuple[str, ...]]:
        cache_key = (round(float(epsilon), 12), round(float(mu), 12))
        if cache_key in self.root_cache:
            eb, timo, warnings = self.root_cache[cache_key]
            return eb.copy(), timo.copy(), "solver_cache", warnings

        rows, path = self.rows_for_epsilon(float(epsilon))
        if rows:
            by_index: dict[int, tuple[float, float]] = {}
            for row in rows_at_mu(rows, float(mu)):
                try:
                    index = int(float(row["sorted_index"]))
                except (KeyError, ValueError):
                    continue
                eb = row_float(row, "Lambda_EB_plot", "Lambda_EB_raw")
                timo = row_float(row, "Lambda_Timoshenko_plot", "Lambda_Timoshenko_raw")
                if isfinite(eb) and isfinite(timo):
                    by_index[index] = (eb, timo)
            if len(by_index) >= int(self.args.n_roots):
                eb_roots = np.array([by_index[index][0] for index in range(1, int(self.args.n_roots) + 1)], dtype=float)
                timo_roots = np.array([by_index[index][1] for index in range(1, int(self.args.n_roots) + 1)], dtype=float)
                warnings = (f"read sorted roots from {rel(path)}",) if path is not None else ("read sorted roots from CSV",)
                self.root_cache[cache_key] = (eb_roots.copy(), timo_roots.copy(), warnings)
                return eb_roots, timo_roots, "lambda_csv", warnings

        eb_roots, timo_roots, warnings = solve_roots_at_point(
            beta_deg=float(self.args.beta_deg),
            eta=float(self.args.eta),
            epsilon=float(epsilon),
            mu=float(mu),
            n_roots=int(self.args.n_roots),
        )
        self.root_cache[cache_key] = (eb_roots.copy(), timo_roots.copy(), warnings)
        return eb_roots, timo_roots, "solver", warnings

    def root_pair(self, spec: SuspectSpec, mu: float) -> RootPair:
        rows, path = self.rows_for_epsilon(spec.epsilon)
        pair = lambda_from_rows(rows, mu=float(mu), sorted_index=spec.sorted_index) if rows else None
        if pair is not None:
            source = f"lambda_csv:{rel(path)}" if path is not None else "lambda_csv"
            return RootPair(eb=pair.eb, timo=pair.timo, source=source, warnings=pair.warnings)
        eb_roots, timo_roots, source, warnings = self.roots(spec.epsilon, float(mu))
        index = spec.sorted_index - 1
        if index >= len(eb_roots) or index >= len(timo_roots):
            raise RuntimeError(f"Requested sorted index {spec.sorted_index} is outside solved root arrays.")
        eb = float(eb_roots[index])
        timo = float(timo_roots[index])
        if not (isfinite(eb) and isfinite(timo)):
            raise RuntimeError(
                f"Non-finite root for epsilon={spec.epsilon:g}, mu={float(mu):g}, sorted={spec.sorted_index}."
            )
        return RootPair(eb=eb, timo=timo, source=source, warnings=warnings)


def auto_points_for_spec(provider: RootProvider, spec: SuspectSpec) -> list[AutoPoint]:
    rows, path = provider.rows_for_epsilon(spec.epsilon)
    if rows:
        return auto_points_from_rows(rows, spec, path)
    return auto_points_from_solver(provider, spec)


def auto_points_from_rows(rows: Sequence[dict[str, str]], spec: SuspectSpec, path: Path | None) -> list[AutoPoint]:
    source = f"lambda_csv:{rel(path)}" if path is not None else "lambda_csv"
    branch_rows = [
        row
        for row in rows
        if int(float(row.get("sorted_index", "nan"))) == int(spec.sorted_index)
    ]
    rel_candidates: list[tuple[float, float]] = []
    for row in branch_rows:
        rel_diff = row_float(row, "rel_diff_Timo_vs_EB", "rel_diff_abs_Timoshenko_vs_EB")
        mu = row_float(row, "mu")
        if isfinite(rel_diff) and isfinite(mu):
            rel_candidates.append((abs(rel_diff), mu))
    auto: list[AutoPoint] = []
    if rel_candidates:
        rel_diff, mu = min(rel_candidates, key=lambda item: item[0])
        auto.append(AutoPoint(mu=float(mu), reason=f"minimal_EB_Timoshenko_relative_difference={rel_diff:.6g}", source=source))

    values_by_mu_model: dict[float, dict[str, dict[int, float]]] = {}
    for row in rows:
        mu = row_float(row, "mu")
        if not isfinite(mu):
            continue
        index = int(float(row.get("sorted_index", "nan")))
        entry = values_by_mu_model.setdefault(float(mu), {MODEL_EB: {}, MODEL_TIMO: {}})
        eb = row_float(row, "Lambda_EB_plot", "Lambda_EB_raw")
        timo = row_float(row, "Lambda_Timoshenko_plot", "Lambda_Timoshenko_raw")
        if isfinite(eb):
            entry[MODEL_EB][index] = eb
        if isfinite(timo):
            entry[MODEL_TIMO][index] = timo
    gap_candidates: list[tuple[float, float, str]] = []
    for mu, model_values in values_by_mu_model.items():
        for model, roots_by_index in model_values.items():
            current = roots_by_index.get(spec.sorted_index)
            if current is None:
                continue
            gaps = [
                abs(float(current) - float(roots_by_index[neighbor]))
                for neighbor in (spec.sorted_index - 1, spec.sorted_index + 1)
                if neighbor in roots_by_index
            ]
            if gaps:
                gap_candidates.append((min(gaps), float(mu), model))
    if gap_candidates:
        gap, mu, model = min(gap_candidates, key=lambda item: item[0])
        auto.append(AutoPoint(mu=float(mu), reason=f"closest_neighbor_gap={gap:.6g} in {model}", source=source))
    return auto


def auto_points_from_solver(provider: RootProvider, spec: SuspectSpec) -> list[AutoPoint]:
    rel_candidates: list[tuple[float, float]] = []
    gap_candidates: list[tuple[float, float, str]] = []
    for mu in regular_grid(0.0, 0.7, AUTO_GRID_STEP):
        eb_roots, timo_roots, source, _warnings = provider.roots(spec.epsilon, float(mu))
        index = spec.sorted_index - 1
        if index < len(eb_roots) and index < len(timo_roots):
            eb = float(eb_roots[index])
            timo = float(timo_roots[index])
            if isfinite(eb) and isfinite(timo) and abs(eb) > 1.0e-30:
                rel_candidates.append((abs(timo - eb) / abs(eb), float(mu)))
        for model, roots in ((MODEL_EB, eb_roots), (MODEL_TIMO, timo_roots)):
            if index >= len(roots) or not isfinite(float(roots[index])):
                continue
            gaps = []
            for neighbor in (index - 1, index + 1):
                if 0 <= neighbor < len(roots) and isfinite(float(roots[neighbor])):
                    gaps.append(abs(float(roots[index]) - float(roots[neighbor])))
            if gaps:
                gap_candidates.append((min(gaps), float(mu), model))
    auto: list[AutoPoint] = []
    if rel_candidates:
        rel_diff, mu = min(rel_candidates, key=lambda item: item[0])
        auto.append(AutoPoint(mu=float(mu), reason=f"minimal_EB_Timoshenko_relative_difference={rel_diff:.6g}", source="solver_auto_grid"))
    if gap_candidates:
        gap, mu, model = min(gap_candidates, key=lambda item: item[0])
        auto.append(AutoPoint(mu=float(mu), reason=f"closest_neighbor_gap={gap:.6g} in {model}", source="solver_auto_grid"))
    return auto


def selected_mu_values(args: argparse.Namespace, provider: RootProvider, spec: SuspectSpec) -> tuple[list[float], list[str]]:
    values = [float(mu) for mu in args.mu_values]
    notes: list[str] = [f"requested_mu={number_text(mu)}" for mu in values]
    if bool(args.include_auto_mu_points):
        for point in auto_points_for_spec(provider, spec):
            if any(abs(float(point.mu) - existing) <= 5.0e-10 for existing in values):
                notes.append(f"auto_duplicate_skipped:{point.reason} at mu={number_text(point.mu)}")
                continue
            values.append(float(point.mu))
            notes.append(f"auto_added:{point.reason} at mu={number_text(point.mu)} from {point.source}")
    return sorted(values), notes


def section_pair(epsilon: float, mu: float, eta: float) -> tuple[TIMO.Section, TIMO.Section]:
    factors = TIMO.tau_factors(float(mu), float(eta))
    return (
        TIMO.section_from_epsilon_tau(float(epsilon), factors.tau1),
        TIMO.section_from_epsilon_tau(float(epsilon), factors.tau2),
    )


def trapz(values: np.ndarray, coordinates: np.ndarray) -> float:
    integrate = getattr(np, "trapezoid", None)
    if integrate is not None:
        return float(integrate(np.asarray(values, dtype=float), np.asarray(coordinates, dtype=float)))
    y = np.asarray(values, dtype=float)
    x = np.asarray(coordinates, dtype=float)
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if abs(float(denominator)) > 1.0e-30 else float("nan")


def endpoint_value(rod: dict[str, object], name: str, index: int) -> float:
    return float(np.asarray(rod[name], dtype=float)[index])


def timo_joint_continuity_row(
    *,
    epsilon: float,
    beta_deg: float,
    eta: float,
    mu: float,
    sorted_index: int,
    Lambda: float,
    fields: dict[str, object],
) -> dict[str, object]:
    rod1 = fields["rod1"]
    rod2 = fields["rod2"]
    beta_rad = np.deg2rad(float(beta_deg))
    cb = float(np.cos(beta_rad))
    sb = float(np.sin(beta_rad))

    # Timoshenko coupling rows are written at rod1 x=+l1 and rod2 x=-l2.
    i1 = -1
    i2 = -1
    w1 = endpoint_value(rod1, "w", i1)
    u1 = endpoint_value(rod1, "u", i1)
    psi1 = endpoint_value(rod1, "psi", i1)
    m1 = endpoint_value(rod1, "M", i1)
    q1 = endpoint_value(rod1, "Q", i1)
    n1 = endpoint_value(rod1, "N", i1)

    w2 = endpoint_value(rod2, "w", i2)
    u2 = endpoint_value(rod2, "u", i2)
    psi2 = endpoint_value(rod2, "psi", i2)
    m2 = endpoint_value(rod2, "M", i2)
    q2 = endpoint_value(rod2, "Q", i2)
    n2 = endpoint_value(rod2, "N", i2)

    gap_w = w1 - cb * w2 + sb * u2
    gap_u = u1 - sb * w2 - cb * u2
    gap_psi = psi1 - psi2
    gap_m = m1 - m2
    gap_q = -q1 + cb * q2 - sb * n2
    gap_n = n1 - sb * q2 - cb * n2

    kinematic_scale = max(1.0, *(abs(value) for value in (w1, u1, psi1, w2, u2, psi2)))
    force_scale = max(1.0, *(abs(value) for value in (m1, q1, n1, m2, q2, n2)))
    max_abs_kinematic_gap = max(abs(gap_w), abs(gap_u), abs(gap_psi)) / kinematic_scale
    max_abs_force_gap = max(abs(gap_m), abs(gap_q), abs(gap_n)) / force_scale

    return {
        "epsilon": float(epsilon),
        "mu": float(mu),
        "beta_deg": float(beta_deg),
        "eta": float(eta),
        "model": MODEL_TIMO,
        "sorted_index": int(sorted_index),
        "Lambda": float(Lambda),
        "gap_w": float(gap_w),
        "gap_u": float(gap_u),
        "gap_psi": float(gap_psi),
        "gap_M": float(gap_m),
        "gap_Q": float(gap_q),
        "gap_N": float(gap_n),
        "max_abs_kinematic_gap": float(max_abs_kinematic_gap),
        "max_abs_force_gap": float(max_abs_force_gap),
        "pass_kinematic": bool(max_abs_kinematic_gap <= JOINT_KINEMATIC_TOL),
        "pass_force": bool(max_abs_force_gap <= JOINT_FORCE_TOL),
        "notes": (
            "rod1 endpoint x=+l1, rod2 endpoint x=-l2; "
            "gaps normalized by max(1, endpoint amplitudes)"
        ),
    }


def sign_factor_from_components(components: PlotComponents) -> tuple[float, str, str]:
    labels = ("u_left", "w_left", "u_right", "w_right")
    arrays = [getattr(components, label) for label in labels]
    flat = np.concatenate([np.asarray(array, dtype=float) for array in arrays])
    if flat.size == 0 or not np.any(np.isfinite(flat)):
        return 1.0, "largest absolute displacement component positive", "nonfinite_or_empty_displacement"
    index = int(np.nanargmax(np.abs(flat)))
    value = float(flat[index])
    if abs(value) <= 1.0e-14:
        return 1.0, "largest absolute displacement component positive", "zero_sign_reference"
    return (-1.0 if value < 0.0 else 1.0), "largest absolute displacement component positive", ""


def scale_components(components: PlotComponents, factor: float) -> PlotComponents:
    factor_f = float(factor)
    return PlotComponents(
        u_left=factor_f * np.asarray(components.u_left, dtype=float),
        w_left=factor_f * np.asarray(components.w_left, dtype=float),
        u_right=factor_f * np.asarray(components.u_right, dtype=float),
        w_right=factor_f * np.asarray(components.w_right, dtype=float),
    )


def max_component_amplitudes(components: PlotComponents) -> tuple[float, float]:
    max_u = max(float(np.max(np.abs(components.u_left))), float(np.max(np.abs(components.u_right))))
    max_w = max(float(np.max(np.abs(components.w_left))), float(np.max(np.abs(components.w_right))))
    return max_u, max_w


def normalized_plot_components(components: PlotComponents) -> PlotComponents:
    max_u, max_w = max_component_amplitudes(components)
    scale = max(max_u, max_w, 1.0e-14)
    return scale_components(components, 1.0 / scale)


def classify_mode(model: str, axial_fraction: float, bending_fraction: float, shear_fraction: float) -> str:
    bending_shear = float(bending_fraction) + float(shear_fraction)
    if isfinite(float(axial_fraction)) and float(axial_fraction) >= 0.70:
        return "longitudinal_dominated"
    if model == MODEL_TIMO:
        if isfinite(bending_shear) and bending_shear >= 0.70:
            return "bending_dominated"
    elif isfinite(float(bending_fraction)) and float(bending_fraction) >= 0.70:
        return "bending_dominated"
    return "mixed"


def energy_row(
    *,
    epsilon: float,
    beta_deg: float,
    eta: float,
    mu: float,
    model: str,
    sorted_index: int,
    Lambda: float,
    U_axial: float,
    U_bending: float,
    U_shear: float,
    axial_displacement_fraction: float,
    transverse_displacement_fraction: float,
    notes: str,
) -> dict[str, object]:
    U_total = float(U_axial) + float(U_bending) + float(U_shear)
    axial_fraction = safe_ratio(U_axial, U_total)
    bending_fraction = safe_ratio(U_bending, U_total)
    shear_fraction = safe_ratio(U_shear, U_total)
    bending_shear = bending_fraction + shear_fraction if isfinite(bending_fraction) and isfinite(shear_fraction) else float("nan")
    return {
        "epsilon": float(epsilon),
        "beta_deg": float(beta_deg),
        "eta": float(eta),
        "mu": float(mu),
        "model": model,
        "sorted_index": int(sorted_index),
        "Lambda": float(Lambda),
        "U_axial": float(U_axial),
        "U_bending": float(U_bending),
        "U_shear": float(U_shear),
        "U_total": float(U_total),
        "axial_energy_fraction": float(axial_fraction),
        "bending_energy_fraction": float(bending_fraction),
        "shear_energy_fraction": float(shear_fraction),
        "bending_shear_fraction": float(bending_shear),
        "axial_displacement_fraction": float(axial_displacement_fraction),
        "transverse_displacement_fraction": float(transverse_displacement_fraction),
        "classification": classify_mode(model, axial_fraction, bending_fraction, shear_fraction),
        "notes": notes,
    }


def eb_mode_diagnostic(
    *,
    epsilon: float,
    beta_deg: float,
    eta: float,
    mu: float,
    sorted_index: int,
    Lambda: float,
    n_points: int,
    root_source: str,
    root_warnings: Sequence[str],
) -> ModeDiagnostic:
    factors = TIMO.tau_factors(float(mu), float(eta))
    section1, section2 = section_pair(float(epsilon), float(mu), float(eta))
    l1, l2 = TIMO.segment_lengths(float(mu))
    matrix = assemble_clamped_coupled_matrix_eta(
        float(Lambda),
        float(np.deg2rad(float(beta_deg))),
        float(mu),
        float(epsilon),
        float(eta),
    )
    coeff, smallest, ratio = analytic_null_vector(matrix)
    residual = np.asarray(matrix, dtype=float) @ np.asarray(coeff, dtype=float)
    A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]

    x1 = np.linspace(0.0, l1, int(n_points), dtype=float)
    x2 = np.linspace(0.0, l2, int(n_points), dtype=float)
    z1 = float(Lambda) * x1 / np.sqrt(factors.tau1)
    th1 = float(epsilon) * float(Lambda) ** 2 * x1
    z2 = -float(Lambda) * x2 / np.sqrt(factors.tau2)
    th2 = -float(epsilon) * float(Lambda) ** 2 * x2

    w1 = A1 * (np.cos(z1) - np.cosh(z1)) + B1 * (np.sin(z1) - np.sinh(z1))
    u1 = P1 * np.sin(th1)
    u1_prime = P1 * float(epsilon) * float(Lambda) ** 2 * np.cos(th1)
    w1_second = (float(Lambda) / np.sqrt(factors.tau1)) ** 2 * (
        A1 * (-np.cos(z1) - np.cosh(z1)) + B1 * (-np.sin(z1) - np.sinh(z1))
    )

    w2_ext = A2 * (np.cos(z2) - np.cosh(z2)) + B2 * (np.sin(z2) - np.sinh(z2))
    u2_ext = P2 * np.sin(th2)
    u2_prime = -P2 * float(epsilon) * float(Lambda) ** 2 * np.cos(th2)
    w2_second = (float(Lambda) / np.sqrt(factors.tau2)) ** 2 * (
        A2 * (-np.cos(z2) - np.cosh(z2)) + B2 * (-np.sin(z2) - np.sinh(z2))
    )

    plot_components = PlotComponents(
        u_left=u1,
        w_left=w1,
        u_right=u2_ext[::-1],
        w_right=w2_ext[::-1],
    )
    sign_factor, sign_normalization, sign_warning = sign_factor_from_components(plot_components)
    plot_components = scale_components(plot_components, sign_factor)

    U_axial = 0.5 * TIMO.E * section1.area * trapz(u1_prime**2, x1)
    U_axial += 0.5 * TIMO.E * section2.area * trapz(u2_prime**2, x2)
    U_bending = 0.5 * section1.bending_stiffness * trapz(w1_second**2, x1)
    U_bending += 0.5 * section2.bending_stiffness * trapz(w2_second**2, x2)
    disp_u = trapz(u1**2, x1) + trapz(u2_ext**2, x2)
    disp_w = trapz(w1**2, x1) + trapz(w2_ext**2, x2)

    warnings = list(root_warnings)
    if sign_warning:
        warnings.append(sign_warning)
    if float(np.max(np.abs(residual))) > 1.0e-7:
        warnings.append("EB_matrix_residual_gt_1e-7")
    notes = (
        f"source={root_source}; sign_factor={sign_factor:+.0f}; "
        f"smallest_singular_value={smallest:.6g}; singular_value_ratio={ratio:.6g}; "
        "sorted frequency, no descendant tracking"
    )
    row = energy_row(
        epsilon=epsilon,
        beta_deg=beta_deg,
        eta=eta,
        mu=mu,
        model=MODEL_EB,
        sorted_index=sorted_index,
        Lambda=Lambda,
        U_axial=max(float(U_axial), 0.0),
        U_bending=max(float(U_bending), 0.0),
        U_shear=0.0,
        axial_displacement_fraction=safe_ratio(disp_u, disp_u + disp_w),
        transverse_displacement_fraction=safe_ratio(disp_w, disp_u + disp_w),
        notes=notes,
    )
    shape_row = {
        "epsilon": float(epsilon),
        "mu": float(mu),
        "model": MODEL_EB,
        "sorted_index": int(sorted_index),
        "Lambda": float(Lambda),
        "figure_file": "",
        "sign_normalization": sign_normalization,
        "warnings": "; ".join(warnings) if warnings else "none",
        "notes": notes,
    }
    return ModeDiagnostic(
        epsilon=float(epsilon),
        beta_deg=float(beta_deg),
        eta=float(eta),
        mu=float(mu),
        model=MODEL_EB,
        sorted_index=int(sorted_index),
        Lambda=float(Lambda),
        plot_components=plot_components,
        energy_row=row,
        shape_row=shape_row,
        warnings=tuple(warnings),
    )


def timo_mode_diagnostic(
    *,
    epsilon: float,
    beta_deg: float,
    eta: float,
    mu: float,
    sorted_index: int,
    Lambda: float,
    n_points: int,
    root_source: str,
    root_warnings: Sequence[str],
) -> ModeDiagnostic:
    mode = TIMO.timo_mode_coefficients(
        float(Lambda),
        float(beta_deg),
        float(mu),
        float(epsilon),
        float(eta),
    )
    fields = TIMO.timo_mode_fields(
        float(Lambda),
        float(beta_deg),
        float(mu),
        float(epsilon),
        float(eta),
        coeff=mode.coeff,
        n_points=int(n_points),
    )
    rod1 = fields["rod1"]
    rod2 = fields["rod2"]
    section1 = fields["section1"]
    section2 = fields["section2"]
    joint_row = timo_joint_continuity_row(
        epsilon=epsilon,
        beta_deg=beta_deg,
        eta=eta,
        mu=mu,
        sorted_index=sorted_index,
        Lambda=Lambda,
        fields=fields,
    )
    plot_components = PlotComponents(
        u_left=np.asarray(rod1["u"], dtype=float),
        w_left=np.asarray(rod1["w"], dtype=float),
        u_right=np.asarray(rod2["u"], dtype=float)[::-1],
        w_right=np.asarray(rod2["w"], dtype=float)[::-1],
    )
    sign_factor, sign_normalization, sign_warning = sign_factor_from_components(plot_components)
    plot_components = scale_components(plot_components, sign_factor)

    x1 = np.abs(np.asarray(rod1["x"], dtype=float))
    x2 = np.abs(np.asarray(rod2["x"], dtype=float))
    gamma1 = np.asarray(rod1["w_prime"], dtype=float) - np.asarray(rod1["psi"], dtype=float)
    gamma2 = np.asarray(rod2["w_prime"], dtype=float) - np.asarray(rod2["psi"], dtype=float)
    U_axial = 0.5 * TIMO.E * section1.area * trapz(np.asarray(rod1["u_prime"], dtype=float) ** 2, x1)
    U_axial += 0.5 * TIMO.E * section2.area * trapz(np.asarray(rod2["u_prime"], dtype=float) ** 2, x2)
    U_bending = 0.5 * section1.bending_stiffness * trapz(np.asarray(rod1["psi_prime"], dtype=float) ** 2, x1)
    U_bending += 0.5 * section2.bending_stiffness * trapz(np.asarray(rod2["psi_prime"], dtype=float) ** 2, x2)
    U_shear = 0.5 * section1.shear_stiffness * trapz(gamma1**2, x1)
    U_shear += 0.5 * section2.shear_stiffness * trapz(gamma2**2, x2)
    disp_u = trapz(np.asarray(rod1["u"], dtype=float) ** 2, x1) + trapz(np.asarray(rod2["u"], dtype=float) ** 2, x2)
    disp_w = trapz(np.asarray(rod1["w"], dtype=float) ** 2, x1) + trapz(np.asarray(rod2["w"], dtype=float) ** 2, x2)

    matrix, matrix_warnings = TIMO.timo_coupling_matrix(
        float(Lambda),
        float(beta_deg),
        float(mu),
        float(epsilon),
        float(eta),
    )
    residual = TIMO.row_normalized(matrix) @ np.asarray(mode.coeff, dtype=float)
    warnings = list(root_warnings)
    warnings.extend(matrix_warnings)
    warnings.extend(mode.warnings)
    warnings.extend(fields.get("warnings", ()))
    if sign_warning:
        warnings.append(sign_warning)
    if float(np.max(np.abs(residual))) > 1.0e-7:
        warnings.append("Timoshenko_row_normalized_matrix_residual_gt_1e-7")
    if not bool(joint_row["pass_kinematic"]):
        warnings.append("Timoshenko_joint_kinematic_gap_gt_1e-6")
    if not bool(joint_row["pass_force"]):
        warnings.append("Timoshenko_joint_force_gap_gt_1e-6")
    warnings = list(dict.fromkeys(warnings))
    notes = (
        f"source={root_source}; sign_factor={sign_factor:+.0f}; "
        f"smallest_singular_value={mode.smallest_singular_value:.6g}; "
        f"singular_value_ratio={mode.singular_value_ratio:.6g}; "
        f"joint_kinematic_gap={float(joint_row['max_abs_kinematic_gap']):.6g}; "
        "right-rod plot transform uses the reflected display basis; "
        "sorted frequency, no descendant tracking"
    )
    row = energy_row(
        epsilon=epsilon,
        beta_deg=beta_deg,
        eta=eta,
        mu=mu,
        model=MODEL_TIMO,
        sorted_index=sorted_index,
        Lambda=Lambda,
        U_axial=max(float(U_axial), 0.0),
        U_bending=max(float(U_bending), 0.0),
        U_shear=max(float(U_shear), 0.0),
        axial_displacement_fraction=safe_ratio(disp_u, disp_u + disp_w),
        transverse_displacement_fraction=safe_ratio(disp_w, disp_u + disp_w),
        notes=notes,
    )
    shape_row = {
        "epsilon": float(epsilon),
        "mu": float(mu),
        "model": MODEL_TIMO,
        "sorted_index": int(sorted_index),
        "Lambda": float(Lambda),
        "figure_file": "",
        "sign_normalization": sign_normalization,
        "warnings": "; ".join(warnings) if warnings else "none",
        "notes": notes,
    }
    return ModeDiagnostic(
        epsilon=float(epsilon),
        beta_deg=float(beta_deg),
        eta=float(eta),
        mu=float(mu),
        model=MODEL_TIMO,
        sorted_index=int(sorted_index),
        Lambda=float(Lambda),
        plot_components=plot_components,
        energy_row=row,
        shape_row=shape_row,
        warnings=tuple(warnings),
        joint_row=joint_row,
    )


def build_suspect_diagnostics(
    args: argparse.Namespace,
    provider: RootProvider,
    selected_by_spec: dict[SuspectSpec, list[float]],
) -> list[ModeDiagnostic]:
    diagnostics: list[ModeDiagnostic] = []
    for spec in selected_by_spec:
        for mu in selected_by_spec[spec]:
            pair = provider.root_pair(spec, float(mu))
            diagnostics.append(
                eb_mode_diagnostic(
                    epsilon=spec.epsilon,
                    beta_deg=float(args.beta_deg),
                    eta=float(args.eta),
                    mu=float(mu),
                    sorted_index=spec.sorted_index,
                    Lambda=pair.eb,
                    n_points=int(args.n_points),
                    root_source=pair.source,
                    root_warnings=pair.warnings,
                )
            )
            diagnostics.append(
                timo_mode_diagnostic(
                    epsilon=spec.epsilon,
                    beta_deg=float(args.beta_deg),
                    eta=float(args.eta),
                    mu=float(mu),
                    sorted_index=spec.sorted_index,
                    Lambda=pair.timo,
                    n_points=int(args.n_points),
                    root_source=pair.source,
                    root_warnings=pair.warnings,
                )
            )
    return diagnostics


def build_control_rows(args: argparse.Namespace, provider: RootProvider) -> list[dict[str, object]]:
    if not bool(args.include_control_modes):
        return []
    rows: list[dict[str, object]] = []
    for spec in args.suspect:
        indices = (spec.sorted_index - 1, spec.sorted_index, spec.sorted_index + 1)
        eb_roots, timo_roots, source, warnings = provider.roots(spec.epsilon, 0.35)
        for sorted_index in indices:
            if sorted_index <= 0:
                continue
            idx = sorted_index - 1
            if idx >= len(eb_roots) or idx >= len(timo_roots):
                continue
            if isfinite(float(eb_roots[idx])):
                rows.append(
                    eb_mode_diagnostic(
                        epsilon=spec.epsilon,
                        beta_deg=float(args.beta_deg),
                        eta=float(args.eta),
                        mu=0.35,
                        sorted_index=sorted_index,
                        Lambda=float(eb_roots[idx]),
                        n_points=int(args.n_points),
                        root_source=source,
                        root_warnings=warnings,
                    ).energy_row
                )
            if isfinite(float(timo_roots[idx])):
                rows.append(
                    timo_mode_diagnostic(
                        epsilon=spec.epsilon,
                        beta_deg=float(args.beta_deg),
                        eta=float(args.eta),
                        mu=0.35,
                        sorted_index=sorted_index,
                        Lambda=float(timo_roots[idx]),
                        n_points=int(args.n_points),
                        root_source=source,
                        root_warnings=warnings,
                    ).energy_row
                )
    return rows


def base_coordinates(mu: float, beta_deg: float, n_points: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    l1, l2 = TIMO.segment_lengths(float(mu))
    x1_grid = np.linspace(0.0, l1, int(n_points), dtype=float)
    x2_grid = np.linspace(-l2, 0.0, int(n_points), dtype=float)
    zeros1 = np.zeros_like(x1_grid)
    zeros2 = np.zeros_like(x2_grid)
    rod1_curve = DISPLAY.rod1_local_fields_to_display(x1_grid, zeros1, zeros1, scale=0.0)
    rod2_curve = DISPLAY.rod2_local_fields_to_display(
        x2_grid,
        zeros2,
        zeros2,
        l2=l2,
        x_joint=l1,
        beta_deg=float(beta_deg),
        scale=0.0,
    )
    return rod1_curve.x_base, rod1_curve.y_base, rod2_curve.x_base, rod2_curve.y_base


def local_in_plane_fields_to_global_geometry(
    *,
    model: str,
    rod: int,
    mu: float,
    beta_deg: float,
    u: np.ndarray,
    w: np.ndarray,
    deformation_scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    u_arr = np.asarray(u, dtype=float)
    w_arr = np.asarray(w, dtype=float)
    l1, l2 = TIMO.segment_lengths(float(mu))
    rod1_mapper = (
        DISPLAY.eb_rod1_local_fields_to_display
        if model == MODEL_EB
        else DISPLAY.rod1_local_fields_to_display
    )
    rod2_mapper = (
        DISPLAY.eb_rod2_local_fields_to_display
        if model == MODEL_EB
        else DISPLAY.rod2_local_fields_to_display
    )

    if int(rod) == 1:
        curve = rod1_mapper(
            np.linspace(0.0, l1, len(u_arr), dtype=float),
            u_arr,
            w_arr,
            scale=float(deformation_scale),
        )
    elif int(rod) == 2:
        curve = rod2_mapper(
            np.linspace(-l2, 0.0, len(u_arr), dtype=float),
            u_arr,
            w_arr,
            l2=l2,
            x_joint=l1,
            beta_deg=float(beta_deg),
            scale=float(deformation_scale),
        )
    else:
        raise ValueError("rod must be 1 or 2")
    return curve.x_base, curve.y_base, curve.x_deformed, curve.y_deformed


def deformed_coordinates(
    diagnostic: ModeDiagnostic,
    *,
    scale: float = MODE_SCALE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    components = normalized_plot_components(diagnostic.plot_components)
    l1, l2 = TIMO.segment_lengths(diagnostic.mu)
    length_scale = max(l1, l2, 1.0)
    amp = float(scale) * length_scale
    _x_left, _y_left, x_left_def, y_left_def = local_in_plane_fields_to_global_geometry(
        model=diagnostic.model,
        rod=1,
        mu=diagnostic.mu,
        beta_deg=diagnostic.beta_deg,
        u=components.u_left,
        w=components.w_left,
        deformation_scale=amp,
    )
    _x_right, _y_right, x_right_def, y_right_def = local_in_plane_fields_to_global_geometry(
        model=diagnostic.model,
        rod=2,
        mu=diagnostic.mu,
        beta_deg=diagnostic.beta_deg,
        u=components.u_right,
        w=components.w_right,
        deformation_scale=amp,
    )
    return x_left_def, y_left_def, x_right_def, y_right_def


def axis_limits(diagnostics: Sequence[ModeDiagnostic]) -> tuple[tuple[float, float], tuple[float, float]]:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for diagnostic in diagnostics:
        n_points = len(diagnostic.plot_components.u_left)
        x_left, y_left, x_right, y_right = base_coordinates(diagnostic.mu, diagnostic.beta_deg, n_points)
        xld, yld, xrd, yrd = deformed_coordinates(diagnostic)
        xs.extend([x_left, x_right, xld, xrd])
        ys.extend([y_left, y_right, yld, yrd])
    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    x_span = max(float(np.max(x_all) - np.min(x_all)), 1.0)
    y_span = max(float(np.max(y_all) - np.min(y_all)), 0.5)
    return (
        (float(np.min(x_all) - 0.07 * x_span), float(np.max(x_all) + 0.07 * x_span)),
        (float(np.min(y_all) - 0.18 * y_span), float(np.max(y_all) + 0.18 * y_span)),
    )


def draw_shape(
    ax: plt.Axes,
    diagnostic: ModeDiagnostic,
    *,
    limits: tuple[tuple[float, float], tuple[float, float]],
    labels: bool,
) -> None:
    n_points = len(diagnostic.plot_components.u_left)
    x_left, y_left, x_right, y_right = base_coordinates(diagnostic.mu, diagnostic.beta_deg, n_points)
    xld, yld, xrd, yrd = deformed_coordinates(diagnostic)
    color = "#006BA4" if diagnostic.model == MODEL_EB else "#C85200"
    ax.plot(x_left, y_left, color="0.70", linestyle="--", linewidth=1.0, label="undeformed" if labels else None)
    ax.plot(x_right, y_right, color="0.70", linestyle="--", linewidth=1.0)
    ax.plot(xld, yld, color=color, linewidth=2.0, label=diagnostic.model if labels else None)
    ax.plot(xrd, yrd, color=color, linewidth=2.0)
    ax.scatter([x_left[-1]], [y_left[-1]], color="black", s=18, zorder=5, label="joint" if labels else None)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.grid(True, color="0.88", linewidth=0.6)


def shape_title(diagnostic: ModeDiagnostic) -> str:
    row = diagnostic.energy_row
    return (
        f"{diagnostic.model}, eps={diagnostic.epsilon:g}, beta={diagnostic.beta_deg:g} deg, "
        f"eta={diagnostic.eta:g}, mu={diagnostic.mu:g}, sorted {diagnostic.sorted_index}\n"
        f"Lambda={diagnostic.Lambda:.7g}, axial energy fraction={float(row['axial_energy_fraction']):.3f}"
    )


def individual_png_path(output_dir: Path, diagnostic: ModeDiagnostic) -> Path:
    model_token = "EB" if diagnostic.model == MODEL_EB else "Timoshenko"
    return output_dir / (
        f"shape_eps{eps_token(diagnostic.epsilon)}_mu{mu_token(diagnostic.mu)}_"
        f"sorted{diagnostic.sorted_index}_{model_token}.png"
    )


def plot_individual(diagnostic: ModeDiagnostic, output_dir: Path) -> Path:
    output = individual_png_path(output_dir, diagnostic)
    output.parent.mkdir(parents=True, exist_ok=True)
    limits = axis_limits([diagnostic])
    fig, (ax_shape, ax_bar) = plt.subplots(
        1,
        2,
        figsize=(8.6, 4.5),
        gridspec_kw={"width_ratios": [3.2, 1.0]},
    )
    draw_shape(ax_shape, diagnostic, limits=limits, labels=True)
    ax_shape.set_xlabel("x")
    ax_shape.set_ylabel("y")
    ax_shape.legend(loc="best", fontsize=8, frameon=False)
    ax_shape.set_title(shape_title(diagnostic), fontsize=9.5)

    max_u, max_w = max_component_amplitudes(diagnostic.plot_components)
    denom = max(max_u, max_w, 1.0e-14)
    ax_bar.barh([0, 1], [max_u / denom, max_w / denom], color=["#4C78A8", "#F58518"])
    ax_bar.set_yticks([0, 1], labels=["max |u|", "max |w|"])
    ax_bar.set_xlim(0.0, 1.05)
    ax_bar.grid(True, axis="x", color="0.88", linewidth=0.6)
    ax_bar.tick_params(labelsize=8)
    ax_bar.set_title("displacement scale", fontsize=9)
    fig.tight_layout()
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_grid(output_dir: Path, spec: SuspectSpec, diagnostics: Sequence[ModeDiagnostic]) -> Path:
    case_diags = [
        diagnostic
        for diagnostic in diagnostics
        if abs(diagnostic.epsilon - spec.epsilon) <= 1.0e-12 and diagnostic.sorted_index == spec.sorted_index
    ]
    mu_values = sorted({round(float(diagnostic.mu), 12) for diagnostic in case_diags})
    by_key = {(round(float(diagnostic.mu), 12), diagnostic.model): diagnostic for diagnostic in case_diags}
    n_rows = len(mu_values)
    n_cols = 2
    output = output_dir / f"suspect_shapes_eps{eps_token(spec.epsilon)}_sorted{spec.sorted_index}_grid.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    limits = axis_limits(case_diags)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(8.5, 2.35 * n_rows + 0.7), squeeze=False)
    for row_idx, mu in enumerate(mu_values):
        for col_idx, model in enumerate(MODELS):
            ax = axes[row_idx, col_idx]
            diagnostic = by_key[(mu, model)]
            draw_shape(ax, diagnostic, limits=limits, labels=False)
            energy = diagnostic.energy_row
            ax.set_title(
                (
                    f"{model}, mu={float(mu):g}, Lambda={diagnostic.Lambda:.5g}\n"
                    f"axial={float(energy['axial_energy_fraction']):.3f}, "
                    f"class={energy['classification']}"
                ),
                fontsize=8.5,
            )
            ax.tick_params(labelsize=7)
            if row_idx == n_rows - 1:
                ax.set_xlabel("x", fontsize=8)
            if col_idx == 0:
                ax.set_ylabel("y", fontsize=8)
    fig.suptitle(
        f"Sorted suspect mode shapes: eps={spec.epsilon:g}, beta=45 deg, eta=0, sorted {spec.sorted_index}",
        fontsize=11,
        y=0.992,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.955), h_pad=0.8, w_pad=0.55)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def attach_shape_figures(
    output_dir: Path,
    specs: Sequence[SuspectSpec],
    diagnostics: Sequence[ModeDiagnostic],
) -> tuple[list[ModeDiagnostic], list[Path]]:
    updated: list[ModeDiagnostic] = []
    paths: list[Path] = []
    for diagnostic in diagnostics:
        path = plot_individual(diagnostic, output_dir)
        paths.append(path)
        shape_row = dict(diagnostic.shape_row)
        shape_row["figure_file"] = rel(path)
        updated.append(
            ModeDiagnostic(
                epsilon=diagnostic.epsilon,
                beta_deg=diagnostic.beta_deg,
                eta=diagnostic.eta,
                mu=diagnostic.mu,
                model=diagnostic.model,
                sorted_index=diagnostic.sorted_index,
                Lambda=diagnostic.Lambda,
                plot_components=diagnostic.plot_components,
                energy_row=diagnostic.energy_row,
                shape_row=shape_row,
                warnings=diagnostic.warnings,
                figure_file=path,
                joint_row=diagnostic.joint_row,
            )
        )
    for spec in specs:
        paths.append(plot_grid(output_dir, spec, updated))
    return updated, paths


def validate_energy_rows(rows: Sequence[dict[str, object]]) -> tuple[bool, list[str]]:
    warnings: list[str] = []
    ok = True
    for row in rows:
        model = str(row["model"])
        total = float(row["U_total"])
        if not isfinite(total) or total <= 0.0:
            ok = False
            warnings.append(
                f"nonpositive energy total for eps={row['epsilon']}, mu={row['mu']}, model={model}, sorted={row['sorted_index']}"
            )
            continue
        fractions = float(row["axial_energy_fraction"]) + float(row["bending_energy_fraction"]) + float(row["shear_energy_fraction"])
        if abs(fractions - 1.0) > ENERGY_SUM_TOL:
            ok = False
            warnings.append(
                f"energy fractions sum to {fractions:.8g} for eps={row['epsilon']}, mu={row['mu']}, model={model}, sorted={row['sorted_index']}"
            )
    return ok, warnings


def joint_rows_from_diagnostics(diagnostics: Sequence[ModeDiagnostic]) -> list[dict[str, object]]:
    return [dict(diagnostic.joint_row) for diagnostic in diagnostics if diagnostic.joint_row is not None]


def validate_joint_rows(rows: Sequence[dict[str, object]]) -> tuple[bool, list[str]]:
    warnings: list[str] = []
    ok = True
    for row in rows:
        kin_gap = float(row["max_abs_kinematic_gap"])
        force_gap = float(row["max_abs_force_gap"])
        if not bool(row["pass_kinematic"]):
            ok = False
            warnings.append(
                (
                    "Timoshenko joint kinematic gap exceeds tolerance: "
                    f"eps={float(row['epsilon']):g}, mu={float(row['mu']):g}, "
                    f"sorted={int(row['sorted_index'])}, gap={kin_gap:.6g}"
                )
            )
        if not bool(row["pass_force"]):
            warnings.append(
                (
                    "Timoshenko joint force gap exceeds diagnostic tolerance: "
                    f"eps={float(row['epsilon']):g}, mu={float(row['mu']):g}, "
                    f"sorted={int(row['sorted_index'])}, gap={force_gap:.6g}"
                )
            )
    return ok, warnings


def conclusion_for_spec(rows: Sequence[dict[str, object]], spec: SuspectSpec) -> str:
    selected = [
        row
        for row in rows
        if abs(float(row["epsilon"]) - spec.epsilon) <= 1.0e-12
        and int(row["sorted_index"]) == spec.sorted_index
    ]
    axial_values = [float(row["axial_energy_fraction"]) for row in selected if isfinite(float(row["axial_energy_fraction"]))]
    if axial_values and min(axial_values) >= 0.70:
        return "longitudinal_dominated"
    if axial_values and max(axial_values) >= 0.40:
        return "mixed_longitudinal_influenced"
    return "not_longitudinal_dominated"


def tex_escape(value: object) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    return text


def frequency_table_rows(rows: Sequence[dict[str, object]]) -> list[str]:
    grouped: dict[tuple[float, float, int], dict[str, float]] = {}
    for row in rows:
        key = (float(row["epsilon"]), float(row["mu"]), int(row["sorted_index"]))
        entry = grouped.setdefault(key, {})
        entry[str(row["model"])] = float(row["Lambda"])
    lines = [
        r"\begin{tabular}{rrrrr}",
        r"\hline",
        r"$\epsilon$ & $\mu$ & sorted & $\Lambda_{\mathrm{EB}}$ & $\Lambda_{\mathrm{Timo}}$ \\",
        r"\hline",
    ]
    for (epsilon, mu, sorted_index), entry in sorted(grouped.items()):
        eb = entry.get(MODEL_EB, float("nan"))
        timo = entry.get(MODEL_TIMO, float("nan"))
        lines.append(
            f"{epsilon:.5g} & {mu:.5g} & {sorted_index:d} & {eb:.7g} & {timo:.7g} \\\\"
        )
    lines.extend([r"\hline", r"\end{tabular}"])
    return lines


def energy_table_rows(rows: Sequence[dict[str, object]]) -> list[str]:
    lines = [
        r"\begin{tabular}{rrrrrrrrl}",
        r"\hline",
        (
            r"$\epsilon$ & $\mu$ & model & sorted & "
            r"$f_a$ & $f_b$ & $f_s$ & $f_{b+s}$ & class \\"
        ),
        r"\hline",
    ]
    for row in sorted(rows, key=lambda item: (float(item["epsilon"]), float(item["mu"]), str(item["model"]))):
        model = "EB" if row["model"] == MODEL_EB else "Timo"
        lines.append(
            (
                f"{float(row['epsilon']):.5g} & {float(row['mu']):.5g} & {model} & "
                f"{int(row['sorted_index'])} & {float(row['axial_energy_fraction']):.3f} & "
                f"{float(row['bending_energy_fraction']):.3f} & {float(row['shear_energy_fraction']):.3f} & "
                f"{float(row['bending_shear_fraction']):.3f} & {tex_escape(row['classification'])} \\\\"
            )
        )
    lines.extend([r"\hline", r"\end{tabular}"])
    return lines


def control_table_rows(rows: Sequence[dict[str, object]]) -> list[str]:
    lines = [
        r"\begin{tabular}{rrrrrl}",
        r"\hline",
        r"$\epsilon$ & model & sorted & $\Lambda$ & $f_a$ & class \\",
        r"\hline",
    ]
    for row in sorted(rows, key=lambda item: (float(item["epsilon"]), str(item["model"]), int(item["sorted_index"]))):
        model = "EB" if row["model"] == MODEL_EB else "Timo"
        lines.append(
            (
                f"{float(row['epsilon']):.5g} & {model} & {int(row['sorted_index'])} & "
                f"{float(row['Lambda']):.7g} & {float(row['axial_energy_fraction']):.3f} & "
                f"{tex_escape(row['classification'])} \\\\"
            )
        )
    lines.extend([r"\hline", r"\end{tabular}"])
    return lines


def joint_audit_table_rows(rows: Sequence[dict[str, object]]) -> list[str]:
    lines = [
        r"\begin{tabular}{rrrrrrl}",
        r"\hline",
        (
            r"$\epsilon$ & $\mu$ & sorted & $\Lambda_{\mathrm{Timo}}$ & "
            r"$\max |g_{\mathrm{kin}}|$ & $\max |g_{\mathrm{force}}|$ & pass \\"
        ),
        r"\hline",
    ]
    for row in sorted(rows, key=lambda item: (float(item["epsilon"]), float(item["mu"]), int(item["sorted_index"]))):
        passed = "yes" if bool(row["pass_kinematic"]) else "no"
        lines.append(
            (
                f"{float(row['epsilon']):.5g} & {float(row['mu']):.5g} & {int(row['sorted_index'])} & "
                f"{float(row['Lambda']):.7g} & {float(row['max_abs_kinematic_gap']):.3e} & "
                f"{float(row['max_abs_force_gap']):.3e} & {passed} \\\\"
            )
        )
    lines.extend([r"\hline", r"\end{tabular}"])
    return lines


def write_tex_report(
    output_dir: Path,
    specs: Sequence[SuspectSpec],
    selected_by_spec: dict[SuspectSpec, list[float]],
    selection_notes: Sequence[str],
    energy_rows: Sequence[dict[str, object]],
    control_rows: Sequence[dict[str, object]],
    joint_rows: Sequence[dict[str, object]],
    grid_paths: Sequence[Path],
    validation_warnings: Sequence[str],
) -> Path:
    path = output_dir / "longitudinal_suspect_modes_report.tex"
    lines: list[str] = [
        r"\documentclass[11pt,a4paper]{article}",
        r"\usepackage[T2A]{fontenc}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[russian]{babel}",
        r"\usepackage{geometry}",
        r"\usepackage{graphicx}",
        r"\geometry{margin=20mm}",
        r"\begin{document}",
        r"\title{Диагностическая проверка подозрительных частот EB/Timoshenko}",
        r"\author{}",
        r"\date{}",
        r"\maketitle",
        r"\section*{Цель}",
        (
            "Проверяется гипотеза, что подозрительные sorted-частоты на графиках "
            r"$\Lambda(\mu)$ являются продольно-доминированными или существенно "
            "продольно-смешанными."
        ),
        r"\section*{Параметры}",
        (
            r"$\beta=45^\circ$, $\eta=0$, $\mu\in[0,0.7]$. "
            "Проверены случаи: "
            + ", ".join(f"$\\epsilon={spec.epsilon:g}$, sorted {spec.sorted_index}" for spec in specs)
            + "."
        ),
        "Выбор частот выполнен по sorted index в каждой точке, без descendant tracking.",
        r"\section*{Почему продольная ветвь может давать малую разницу}",
        (
            "Модель Тимошенко меняет изгибно-сдвиговую часть колебаний: учитываются "
            "сдвиговая деформация и поворот сечения. Продольное уравнение в текущих "
            "EB и Timoshenko реализациях остается одним и тем же. Поэтому если форма "
            "моды определяется главным образом продольной деформацией, поправка "
            "Тимошенко к частоте может оказаться малой даже для высокой sorted-частоты."
        ),
        r"\section*{Исправление визуализации Timoshenko и контроль стыка}",
        (
            "В ходе повторного diagnostic-only аудита обнаружена ошибка именно в построении "
            "Timoshenko-форм: правый стержень отображался в глобальные координаты как обычная "
            "tangent/normal пара. Для коэффициентов стержня 2, используемых в матрице сопряжения "
            r"при \(x_2=-l_2\), нужно применять преобразование "
            r"The determinant relation $u_1=c u_2+s w_2$, $w_1=-s u_2+c w_2$ must not be used directly as Cartesian display components. "
            r"The corrected display mapping is $dX_1=u_1$, $dY_1=-w_1$, $dX_2=c u_2+s w_2$, $dY_2=s u_2-c w_2$. "
            "Старые Timoshenko PNG с разрывом в узле не используются как физическое "
            "свидетельство; новые графики принимаются только после проверки совместимости стыка."
        ),
        (
            "This display-only correction invalidates all older Timoshenko global centerline and vector-field PNGs from this workflow for physical interpretation. "
            "The local fields, frequencies, null vectors, joint residuals, and energy fractions remain valid; formulas, determinants, root solvers, and k' are unchanged."
        ),
        (
            "Энергетические доли пересчитаны теми же выражениями, что и раньше. "
            "Исправление относится к визуализации и явному аудиту узловой совместимости, "
            "а не к аналитическим формулам, определителю, root solver или коэффициенту сдвига."
        ),
        r"\section*{Выбранные точки по \(\mu\)}",
    ]
    for spec in specs:
        values = ", ".join(number_text(mu) for mu in selected_by_spec[spec])
        lines.append(f"$\\epsilon={spec.epsilon:g}$, sorted {spec.sorted_index}: $\\mu={values}$.")
    if selection_notes:
        lines.append(r"\begin{itemize}")
        for note in selection_notes:
            lines.append(r"\item " + tex_escape(note))
        lines.append(r"\end{itemize}")

    lines.extend([r"\section*{Частоты}", *frequency_table_rows(energy_rows)])
    if joint_rows:
        lines.extend(
            [
                r"\section*{Аудит совместимости стыка Timoshenko}",
                (
                    "В таблице приведены нормированные максимальные зазоры для кинематических "
                    r"условий \(w_1-cw_2+s u_2=0\), \(u_1-s w_2-c u_2=0\), "
                    r"\(\psi_1-\psi_2=0\) и для силовых условий. "
                    r"Критерий для кинематики: $\max |g_{\mathrm{kin}}|<10^{-6}$."
                ),
                *joint_audit_table_rows(joint_rows),
                "Все показанные Timoshenko-формы прошли кинематическую проверку стыка.",
            ]
        )
    lines.extend([r"\section*{Энергетические доли}", *energy_table_rows(energy_rows)])
    if control_rows:
        lines.extend(
            [
                r"\section*{Контрольные соседние sorted-моды при \(\mu=0.35\)}",
                *control_table_rows(control_rows),
            ]
        )

    lines.append(r"\section*{Формы колебаний}")
    for grid_path in grid_paths:
        lines.extend(
            [
                r"\begin{figure}[htbp]",
                r"\centering",
                rf"\includegraphics[width=0.98\linewidth]{{{Path(grid_path).name}}}",
                rf"\caption{{{tex_escape(Path(grid_path).name)}}}",
                r"\end{figure}",
            ]
        )

    lines.append(r"\section*{Вывод}")
    for spec in specs:
        conclusion = conclusion_for_spec(energy_rows, spec)
        selected = [
            row
            for row in energy_rows
            if abs(float(row["epsilon"]) - spec.epsilon) <= 1.0e-12
            and int(row["sorted_index"]) == spec.sorted_index
        ]
        axial_values = [float(row["axial_energy_fraction"]) for row in selected]
        axial_range = f"{min(axial_values):.3f}--{max(axial_values):.3f}" if axial_values else "nan"
        if conclusion == "longitudinal_dominated":
            text = "диагностика поддерживает классификацию longitudinal-dominated"
        elif conclusion == "mixed_longitudinal_influenced":
            text = "мода выглядит смешанной, но с заметным продольным вкладом"
        else:
            text = "диагностика не поддерживает продольное доминирование"
        lines.append(
            (
                f"Для $\\epsilon={spec.epsilon:g}$, sorted {spec.sorted_index}: {text}; "
                f"диапазон axial energy fraction по проанализированным EB/Timoshenko точкам: {axial_range}."
            )
        )
    lines.append(
        "Эти выводы являются диагностическими: классификация основана на энергетических долях, "
        "а не только на визуальном виде формы."
    )
    if validation_warnings:
        lines.append(r"\section*{Предупреждения}")
        lines.append(r"\begin{itemize}")
        for warning in validation_warnings:
            lines.append(r"\item " + tex_escape(warning))
        lines.append(r"\end{itemize}")
    lines.extend(
        [
            r"\section*{Явные ограничения}",
            r"\begin{itemize}",
            r"\item FEM, 3D FEM, Gmsh и CalculiX не запускались.",
            r"\item Аналитические формулы, determinant/root solver и коэффициент сдвига \(k'\) не изменялись.",
            r"\item Использовались sorted frequencies; descendant tracking не использовался.",
            r"\item Article workspace, main.tex, article figures, old determinant, old solvers и baseline results не затрагивались.",
            r"\end{itemize}",
            r"\end{document}",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def compile_tex(tex_path: Path) -> tuple[Path | None, str]:
    if shutil.which("pdflatex") is None:
        return None, "pdflatex not found"
    log_path = tex_path.with_suffix(".pdflatex.log")
    combined_stdout = []
    for _ in range(2):
        completed = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tex_path.name],
            cwd=tex_path.parent,
            check=False,
            capture_output=True,
            text=True,
            timeout=180,
        )
        combined_stdout.append(completed.stdout)
        combined_stdout.append(completed.stderr)
        if completed.returncode != 0:
            log_path.write_text("\n".join(combined_stdout), encoding="utf-8", errors="replace")
            return None, f"pdflatex failed; see {rel(log_path)}"
    log_path.write_text("\n".join(combined_stdout), encoding="utf-8", errors="replace")
    pdf_path = tex_path.with_suffix(".pdf")
    if pdf_path.exists():
        return pdf_path, f"compiled; log {rel(log_path)}"
    return None, f"pdflatex finished but {rel(pdf_path)} was not created"


def copy_grid_paths_for_tex(grid_paths: Sequence[Path], output_dir: Path) -> None:
    for path in grid_paths:
        source = Path(path)
        target = output_dir / source.name
        if source.resolve() != target.resolve() and source.exists():
            shutil.copyfile(source, target)


def validate_tex_figure_refs(tex_path: Path, grid_paths: Sequence[Path]) -> list[str]:
    text = tex_path.read_text(encoding="utf-8")
    warnings = []
    for path in grid_paths:
        name = Path(path).name
        escaped_name = tex_escape(name)
        if name not in text and escaped_name not in text:
            warnings.append(f"TeX report does not reference {name}")
        if not Path(path).exists():
            warnings.append(f"referenced PNG is missing: {rel(Path(path))}")
    return warnings


def print_summary(
    output_paths: Sequence[Path],
    specs: Sequence[SuspectSpec],
    selected_by_spec: dict[SuspectSpec, list[float]],
    energy_rows: Sequence[dict[str, object]],
    control_rows: Sequence[dict[str, object]],
    joint_rows: Sequence[dict[str, object]],
    pdf_path: Path | None,
    compile_status: str,
) -> None:
    print("longitudinal suspect mode audit outputs:")
    for path in output_paths:
        print(f"  {rel(path)}")
    for spec in specs:
        conclusion = conclusion_for_spec(energy_rows, spec)
        print(
            f"case epsilon={spec.epsilon:g}, sorted={spec.sorted_index}: "
            f"mu points={', '.join(number_text(mu) for mu in selected_by_spec[spec])}; "
            f"conclusion={conclusion}"
        )
        for model in MODELS:
            values = [
                float(row["axial_energy_fraction"])
                for row in energy_rows
                if abs(float(row["epsilon"]) - spec.epsilon) <= 1.0e-12
                and int(row["sorted_index"]) == spec.sorted_index
                and row["model"] == model
            ]
            if values:
                print(f"  {model} axial energy fraction range: {min(values):.6g}..{max(values):.6g}")
    if joint_rows:
        max_kin_gap = max(float(row["max_abs_kinematic_gap"]) for row in joint_rows)
        max_force_gap = max(float(row["max_abs_force_gap"]) for row in joint_rows)
        kin_pass_count = sum(1 for row in joint_rows if bool(row["pass_kinematic"]))
        print(
            f"Timoshenko joint continuity: {kin_pass_count}/{len(joint_rows)} kinematic pass; "
            f"max normalized kinematic gap={max_kin_gap:.6g}; max normalized force gap={max_force_gap:.6g}"
        )
    if control_rows:
        print("control modes at mu=0.35:")
        for row in control_rows:
            print(
                f"  eps={float(row['epsilon']):g}, {row['model']}, sorted={int(row['sorted_index'])}, "
                f"axial={float(row['axial_energy_fraction']):.6g}, class={row['classification']}"
            )
    print(f"TeX compilation: {compile_status}")
    if pdf_path is not None:
        print(f"PDF: {rel(pdf_path)}")


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    provider = RootProvider(args)

    selected_by_spec: dict[SuspectSpec, list[float]] = {}
    selection_notes: list[str] = []
    for spec in args.suspect:
        values, notes = selected_mu_values(args, provider, spec)
        selected_by_spec[spec] = values
        selection_notes.extend(f"epsilon={spec.epsilon:g}, sorted={spec.sorted_index}: {note}" for note in notes)

    diagnostics = build_suspect_diagnostics(args, provider, selected_by_spec)
    joint_rows = joint_rows_from_diagnostics(diagnostics)
    joint_csv = write_csv(
        args.output_dir / "timoshenko_joint_continuity_audit.csv",
        joint_rows,
        TIMO_JOINT_CONTINUITY_FIELDS,
    )
    joint_ok, joint_warnings = validate_joint_rows(joint_rows)
    if not joint_ok:
        raise RuntimeError(
            "Timoshenko plotted shapes failed joint kinematic continuity; "
            f"see {rel(joint_csv)}. " + "; ".join(joint_warnings)
        )
    diagnostics, figure_paths = attach_shape_figures(args.output_dir, args.suspect, diagnostics)
    grid_paths = [
        args.output_dir / f"suspect_shapes_eps{eps_token(spec.epsilon)}_sorted{spec.sorted_index}_grid.png"
        for spec in args.suspect
    ]

    energy_rows = [diagnostic.energy_row for diagnostic in diagnostics]
    shape_rows = [diagnostic.shape_row for diagnostic in diagnostics]
    control_rows = build_control_rows(args, provider)
    energy_ok, energy_warnings = validate_energy_rows([*energy_rows, *control_rows])
    warning_rows = sorted({warning for diagnostic in diagnostics for warning in diagnostic.warnings})
    validation_warnings = list(energy_warnings)
    validation_warnings.extend(joint_warnings)
    validation_warnings.extend(f"mode reconstruction warning: {warning}" for warning in warning_rows)

    energy_csv = write_csv(args.output_dir / "suspect_mode_energy_fractions.csv", energy_rows, ENERGY_FIELDS)
    control_csv = write_csv(args.output_dir / "control_mode_energy_fractions.csv", control_rows, ENERGY_FIELDS)
    shape_csv = write_csv(args.output_dir / "suspect_mode_shape_summary.csv", shape_rows, SHAPE_SUMMARY_FIELDS)
    tex_path = write_tex_report(
        args.output_dir,
        args.suspect,
        selected_by_spec,
        selection_notes,
        energy_rows,
        control_rows,
        joint_rows,
        grid_paths,
        validation_warnings,
    )
    copy_grid_paths_for_tex(grid_paths, args.output_dir)
    tex_ref_warnings = validate_tex_figure_refs(tex_path, grid_paths)
    if tex_ref_warnings:
        validation_warnings.extend(tex_ref_warnings)
        tex_path = write_tex_report(
            args.output_dir,
            args.suspect,
            selected_by_spec,
            selection_notes,
            energy_rows,
            control_rows,
            joint_rows,
            grid_paths,
            validation_warnings,
        )

    pdf_path: Path | None = None
    compile_status = "skipped by --no-compile-tex"
    if bool(args.compile_tex):
        pdf_path, compile_status = compile_tex(tex_path)

    output_paths = [energy_csv, control_csv, shape_csv, joint_csv, *figure_paths, tex_path]
    if pdf_path is not None:
        output_paths.append(pdf_path)
    print_summary(output_paths, args.suspect, selected_by_spec, energy_rows, control_rows, joint_rows, pdf_path, compile_status)
    return {
        "energy_csv": energy_csv,
        "control_csv": control_csv,
        "shape_csv": shape_csv,
        "joint_csv": joint_csv,
        "figure_paths": figure_paths,
        "tex_path": tex_path,
        "pdf_path": pdf_path,
        "energy_rows": energy_rows,
        "control_rows": control_rows,
        "joint_rows": joint_rows,
        "selected_by_spec": selected_by_spec,
        "compile_status": compile_status,
    }


if __name__ == "__main__":
    main()
