from __future__ import annotations

import argparse
import csv
import importlib.util
import math
from pathlib import Path
import sys
from types import ModuleType
from typing import Sequence


REPO_ROOT = Path(__file__).resolve().parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import find_first_n_roots_eta  # noqa: E402
from my_project.analytic.solvers_out_of_plane import find_first_n_roots_out_of_plane_with_warnings  # noqa: E402


DEFAULT_EPSILON = 0.0025
DEFAULT_POISSON = 0.3
DEFAULT_TOTAL_LENGTH = 2.0
DEFAULT_N_ANALYTIC_ROOTS_PER_SUBSYSTEM = 12
DEFAULT_N_FEM_MODES = 40
DEFAULT_ROOT_SCAN_STEP = 0.02
DEFAULT_LAMBDA_MAX = 35.0
DEFAULT_OUTPUT_DIR = Path("results") / "full_spectrum_analytic_vs_3d_fem" / "smoke_straight_uniform"
WORKFLOW = "straight_fixed_fixed_cylinder_from_single_rod_helpers"


def _load_module(module_name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _repo_output_dir(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def _fmt(value: object) -> object:
    if isinstance(value, float):
        if not math.isfinite(value):
            return "nan"
        return f"{value:.16e}"
    return value


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row.get(key, "")) for key in fieldnames})


def build_analytic_rows(
    *,
    epsilon: float,
    poisson: float,
    n_roots_per_subsystem: int,
    lambda_max: float,
    root_scan_step: float,
) -> tuple[list[dict[str, object]], tuple[str, ...]]:
    compare_module = _load_module(
        "compare_full_spectrum_analytic_vs_3d_fem_for_smoke",
        REPO_ROOT / "scripts" / "analysis" / "thickness_mismatch" / "audits" / "compare_full_spectrum_analytic_vs_3d_fem.py",
    )
    in_plane = find_first_n_roots_eta(
        beta=0.0,
        mu=0.0,
        epsilon=float(epsilon),
        eta=0.0,
        n_roots=int(n_roots_per_subsystem),
        Lmax0=float(lambda_max),
        scan_step=float(root_scan_step),
    )
    out = find_first_n_roots_out_of_plane_with_warnings(
        beta=0.0,
        mu=0.0,
        epsilon=float(epsilon),
        eta=0.0,
        poisson=float(poisson),
        n_roots=int(n_roots_per_subsystem),
        lambda_max=float(lambda_max),
        scan_step=float(root_scan_step),
    )
    union = compare_module.build_analytic_full_spectrum_union(in_plane, out.roots)
    rows: list[dict[str, object]] = []
    for row in union:
        rows.append(
            {
                "case_id": "smoke_straight_uniform",
                "beta_deg": 0.0,
                "mu": 0.0,
                "eta": 0.0,
                "epsilon": float(epsilon),
                "poisson": float(poisson),
                "analytic_full_index": int(row["analytic_full_index"]),
                "analytic_subsystem": row["analytic_subsystem"],
                "analytic_subsystem_index": int(row["analytic_subsystem_index"]),
                "Lambda": float(row["Lambda"]),
                "notes": "ok",
            }
        )
    warnings = tuple(f"out_of_plane: {item}" for item in out.warnings)
    return rows, warnings


def compare_module() -> ModuleType:
    return _load_module(
        "compare_full_spectrum_analytic_vs_3d_fem_for_smoke_matches",
        REPO_ROOT / "scripts" / "analysis" / "thickness_mismatch" / "audits" / "compare_full_spectrum_analytic_vs_3d_fem.py",
    )


def classify_displacements(mode_shape: dict[int, tuple[float, float, float]]) -> tuple[str, float, float]:
    total = 0.0
    in_plane = 0.0
    out_of_plane = 0.0
    for ux, uy, uz in mode_shape.values():
        total += ux * ux + uy * uy + uz * uz
        in_plane += ux * ux + uy * uy
        out_of_plane += uz * uz
    if total <= 0.0:
        return "mixed_or_unclassified", float("nan"), float("nan")
    in_fraction = in_plane / total
    out_fraction = out_of_plane / total
    if in_fraction >= 0.8:
        return "mostly_in_plane", in_fraction, out_fraction
    if out_fraction >= 0.8:
        return "mostly_out_of_plane", in_fraction, out_fraction
    return "mixed_or_unclassified", in_fraction, out_fraction


def build_raw_fem_rows(
    *,
    single: ModuleType,
    omegas: Sequence[float],
    parse_result: object | None,
    epsilon: float,
    status: str,
    warning: str,
) -> list[dict[str, object]]:
    shape_by_mode = getattr(parse_result, "mode_shapes", {}) if parse_result is not None else {}
    rows: list[dict[str, object]] = []
    if not omegas:
        return [
            {
                "case_id": "smoke_straight_uniform",
                "fem_mode_index": "",
                "Omega_3d_fem": float("nan"),
                "Lambda_3d_fem": float("nan"),
                "fem_mode_character": "classification_unavailable",
                "fem_in_plane_fraction": float("nan"),
                "fem_out_of_plane_fraction": float("nan"),
                "status": status,
                "warning": warning,
                "notes": "no parsed FEM frequencies",
            }
        ]
    for index, omega in enumerate(omegas, start=1):
        shape = shape_by_mode.get(index, {})
        if shape:
            character, in_fraction, out_fraction = classify_displacements(shape)
            notes = "classification from CalculiX FRD displacement vectors"
        else:
            character, in_fraction, out_fraction = "classification_unavailable", float("nan"), float("nan")
            notes = "FRD displacement vectors unavailable for this mode"
        rows.append(
            {
                "case_id": "smoke_straight_uniform",
                "fem_mode_index": index,
                "Omega_3d_fem": float(omega),
                "Lambda_3d_fem": single.lambda_from_omega(float(omega), float(epsilon)),
                "fem_mode_character": character,
                "fem_in_plane_fraction": in_fraction,
                "fem_out_of_plane_fraction": out_fraction,
                "status": status,
                "warning": warning,
                "notes": notes,
            }
        )
    return rows


def build_sorted_match_rows(
    analytic_rows: Sequence[dict[str, object]],
    fem_rows: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    computed_fem = [
        row
        for row in fem_rows
        if isinstance(row.get("fem_mode_index"), int) and math.isfinite(float(row.get("Lambda_3d_fem", float("nan"))))
    ]
    rows: list[dict[str, object]] = []
    for analytic in analytic_rows:
        full_index = int(analytic["analytic_full_index"])
        fem = computed_fem[full_index - 1] if full_index <= len(computed_fem) else None
        lambda_analytic = float(analytic["Lambda"])
        lambda_fem = float("nan") if fem is None else float(fem["Lambda_3d_fem"])
        abs_error = abs(lambda_fem - lambda_analytic) if math.isfinite(lambda_fem) else float("nan")
        rel_error = abs_error / abs(lambda_analytic) if math.isfinite(abs_error) and lambda_analytic != 0.0 else float("nan")
        rows.append(
            {
                "case_id": "smoke_straight_uniform",
                "analytic_full_index": full_index,
                "analytic_subsystem": analytic["analytic_subsystem"],
                "analytic_subsystem_index": analytic["analytic_subsystem_index"],
                "Lambda_analytic": lambda_analytic,
                "fem_mode_index": "" if fem is None else fem["fem_mode_index"],
                "Lambda_3d_fem": lambda_fem,
                "abs_error": abs_error,
                "rel_error": rel_error,
                "fem_mode_character": "classification_unavailable" if fem is None else fem["fem_mode_character"],
                "fem_in_plane_fraction": float("nan") if fem is None else fem["fem_in_plane_fraction"],
                "fem_out_of_plane_fraction": float("nan") if fem is None else fem["fem_out_of_plane_fraction"],
                "match_method": "sorted_index" if fem is not None else "fem_mode_missing",
                "warning": "" if fem is not None else "No FEM mode at this sorted index.",
                "notes": "frequency-order smoke comparison only; do not use as branch identity",
            }
        )
    return rows


def _first_values(rows: Sequence[dict[str, object]], key: str, limit: int = 10) -> list[float]:
    values: list[float] = []
    for row in rows:
        try:
            value = float(row[key])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(value):
            values.append(value)
        if len(values) >= limit:
            break
    return values


def _error_summary(rows: Sequence[dict[str, object]]) -> tuple[float, float]:
    values = sorted(
        float(row["rel_error"])
        for row in rows
        if math.isfinite(float(row.get("rel_error", float("nan"))))
    )
    if not values:
        return float("nan"), float("nan")
    midpoint = len(values) // 2
    if len(values) % 2:
        median = values[midpoint]
    else:
        median = 0.5 * (values[midpoint - 1] + values[midpoint])
    return max(values), median


def _adjacent_pair_summary(values: Sequence[float], *, rel_tol: float = 2.0e-3, limit: int = 10) -> str:
    pairs: list[str] = []
    limited = list(values[:limit])
    for index in range(len(limited) - 1):
        left = limited[index]
        right = limited[index + 1]
        mean = 0.5 * (left + right)
        if mean > 0.0 and abs(right - left) / mean <= rel_tol:
            pairs.append(f"{index + 1}-{index + 2}")
    return ", ".join(pairs) if pairs else "none in first 10"


def _label_match_summary(rows: Sequence[dict[str, object]]) -> str:
    total = 0
    compatible = 0
    unavailable = 0
    mixed = 0
    for row in rows:
        character = str(row.get("fem_mode_character", ""))
        subsystem = str(row.get("analytic_subsystem", ""))
        if character == "classification_unavailable":
            unavailable += 1
            continue
        if character == "mixed_or_unclassified":
            mixed += 1
            continue
        if character not in {"mostly_in_plane", "mostly_out_of_plane"}:
            continue
        total += 1
        if (subsystem == "in_plane" and character == "mostly_in_plane") or (
            subsystem == "out_of_plane" and character == "mostly_out_of_plane"
        ):
            compatible += 1
    if total == 0:
        return f"unavailable or mixed only (mixed={mixed}, unavailable={unavailable})"
    return f"{compatible}/{total} sorted-index matches with definite 3D character; mixed={mixed}, unavailable={unavailable}"


def _shift_summary(analytic_values: Sequence[float], fem_values: Sequence[float]) -> str:
    if not analytic_values or not fem_values:
        return "unavailable"
    ratio = fem_values[0] / analytic_values[0] if analytic_values[0] != 0.0 else float("nan")
    if math.isfinite(ratio) and ratio > 1.5:
        return f"first FEM Lambda is {ratio:.4g} times the first analytic Lambda; clear upward shift"
    if math.isfinite(ratio) and ratio < 0.67:
        return f"first FEM Lambda is {ratio:.4g} times the first analytic Lambda; clear downward shift"
    return f"first FEM/analytic Lambda ratio is {ratio:.4g}"


def write_report(
    *,
    path: Path,
    args: argparse.Namespace,
    analytic_rows: Sequence[dict[str, object]],
    fem_rows: Sequence[dict[str, object]],
    match_rows: Sequence[dict[str, object]],
    nearest_rows: Sequence[dict[str, object]],
    messages: Sequence[str],
    analytic_warnings: Sequence[str],
    generated_files: Sequence[Path],
    actual_3d_run: bool,
    workflow_choice: str,
    radius: float,
    mesh_size: float,
    gmsh_exe: str | None,
    ccx_exe: str | None,
    parse_message: str,
) -> None:
    first_analytic = _first_values(analytic_rows, "Lambda", limit=10)
    first_fem = _first_values(fem_rows, "Lambda_3d_fem", limit=10)
    max_error, median_error = _error_summary(match_rows)
    classes = sorted({str(row.get("fem_mode_character", "")) for row in fem_rows if row.get("fem_mode_character")})
    shift_summary = _shift_summary(first_analytic, first_fem)
    analytic_pairs = _adjacent_pair_summary(first_analytic)
    fem_pairs = _adjacent_pair_summary(first_fem)
    label_summary = _label_match_summary(match_rows)
    nearest_errors = sorted(
        float(row["rel_error"])
        for row in nearest_rows
        if math.isfinite(float(row.get("rel_error", float("nan"))))
    )
    nearest_median = (
        float("nan")
        if not nearest_errors
        else nearest_errors[len(nearest_errors) // 2]
        if len(nearest_errors) % 2
        else 0.5 * (nearest_errors[len(nearest_errors) // 2 - 1] + nearest_errors[len(nearest_errors) // 2])
    )
    lines = [
        "# Full-Spectrum 3D FEM Smoke Case",
        "",
        "Diagnostic-only straight-uniform smoke comparison.",
        "",
        "## Case",
        "",
        f"- epsilon: {args.epsilon:g}",
        f"- poisson: {args.poisson:g}",
        "- beta: 0 deg",
        "- mu: 0",
        "- eta: 0",
        f"- total length: {args.total_length:g}",
        f"- circular radius: {radius:.16g}",
        f"- radius relation: `r0 = 2*epsilon`",
        f"- material: `E = 1`, `rho = 1`, `nu = {args.poisson:g}`",
        f"- requested FEM modes: {args.n_fem_modes}",
        f"- analytic roots per subsystem: {args.n_analytic_roots_per_subsystem}",
        "",
        "## Workflow",
        "",
        f"- selected workflow: `{workflow_choice}`",
        "- reason: the fused two-rod beta=0 helper places the two cylinders on top of each other,",
        "  while this smoke case is the total-length-2 straight uniform limit.",
        f"- Gmsh executable: `{gmsh_exe}`" if gmsh_exe else "- Gmsh executable: not found",
        f"- CalculiX executable: `{ccx_exe}`" if ccx_exe else "- CalculiX executable: not found",
        f"- actual 3D FEM run: {'yes' if actual_3d_run else 'no'}",
        f"- smoke mesh size: {mesh_size:.16g}",
        "- mesh note: coarse smoke mesh, not a convergence mesh.",
        "",
        "## Roots",
        "",
        f"- first 10 analytic full-spectrum roots: {', '.join(f'{value:.8g}' for value in first_analytic)}",
        f"- first 10 3D FEM roots: {', '.join(f'{value:.8g}' for value in first_fem) if first_fem else 'unavailable'}",
        f"- max relative error: {max_error:.8g}" if math.isfinite(max_error) else "- max relative error: unavailable",
        f"- median relative error: {median_error:.8g}" if math.isfinite(median_error) else "- median relative error: unavailable",
        f"- obvious systematic shift: {shift_summary}",
        f"- analytic close adjacent pairs in first 10: {analytic_pairs}",
        f"- FEM close adjacent pairs in first 10: {fem_pairs}",
        f"- analytic-label/3D-character sorted-index compatibility: {label_summary}",
        f"- nearest-frequency max relative error: {max(nearest_errors):.8g}" if nearest_errors else "- nearest-frequency max relative error: unavailable",
        f"- nearest-frequency median relative error: {nearest_median:.8g}" if math.isfinite(nearest_median) else "- nearest-frequency median relative error: unavailable",
        "",
        "## Mode Classification",
        "",
        f"- classification available: {'yes' if any('mostly_' in item or item == 'mixed_or_unclassified' for item in classes) else 'no'}",
        f"- classes present: {', '.join(classes) if classes else 'none'}",
        f"- parser status: {parse_message}",
        "",
        "## Messages",
        "",
    ]
    lines.extend([f"- {message}" for message in messages] or ["- none"])
    lines.extend(["", "## Analytic Warnings", ""])
    lines.extend([f"- {warning}" for warning in analytic_warnings] or ["- none"])
    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- This is one coarse, diagnostic 3D smoke case only.",
            "- Sorted-index matching is reported for orientation and is not branch identity.",
            "- No determinant, old solver, production FEM workflow, or article figure was modified.",
            "- Do not use these numbers as mesh-converged validation.",
            "",
            "## Generated Files",
            "",
        ]
    )
    for item in [*generated_files, path]:
        lines.append(f"- `{item.relative_to(REPO_ROOT)}`")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run one straight-uniform full-spectrum 3D FEM smoke case.")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--poisson", type=float, default=DEFAULT_POISSON)
    parser.add_argument("--total-length", type=float, default=DEFAULT_TOTAL_LENGTH)
    parser.add_argument("--n-analytic-roots-per-subsystem", type=int, default=DEFAULT_N_ANALYTIC_ROOTS_PER_SUBSYSTEM)
    parser.add_argument("--n-fem-modes", type=int, default=DEFAULT_N_FEM_MODES)
    parser.add_argument("--lambda-max", type=float, default=DEFAULT_LAMBDA_MAX)
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    parser.add_argument("--mesh-size-factor", type=float, default=1.0, help="Smoke mesh size as a multiple of r0.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-3d-fem", action="store_true")
    parser.add_argument("--ccx-timeout-seconds", type=int, default=240)
    args = parser.parse_args(argv)
    if args.epsilon <= 0.0:
        raise ValueError("epsilon must be positive.")
    if args.total_length <= 0.0:
        raise ValueError("total length must be positive.")
    if args.mesh_size_factor <= 0.0:
        raise ValueError("mesh size factor must be positive.")
    args.output_dir = _repo_output_dir(args.output_dir)
    return args


def run(args: argparse.Namespace) -> dict[str, object]:
    output_dir = Path(args.output_dir)
    analytic_rows, analytic_warnings = build_analytic_rows(
        epsilon=args.epsilon,
        poisson=args.poisson,
        n_roots_per_subsystem=args.n_analytic_roots_per_subsystem,
        lambda_max=args.lambda_max,
        root_scan_step=args.root_scan_step,
    )

    single = _load_module(
        "solid_fem_single_rod_fixed_fixed_smoke",
        REPO_ROOT / "scripts" / "analysis" / "solid_fem_single_rod_fixed_fixed.py",
    )
    single.L = float(args.total_length)
    single.NU = float(args.poisson)
    single.SOLID_MODES_REQUESTED = int(args.n_fem_modes)
    single.CCX_TIMEOUT_SECONDS = int(args.ccx_timeout_seconds)
    radius = float(args.epsilon) * float(args.total_length)
    mesh_size = float(args.mesh_size_factor) * radius
    single.mesh_size_for_epsilon = lambda _epsilon: mesh_size

    case_dir = output_dir / "solid_fem_case"
    stem = "straight_uniform_L2_eps0p0025"
    paths = single.CasePaths(
        case_dir=case_dir,
        geo=case_dir / f"{stem}.geo",
        msh=case_dir / f"{stem}_mesh.msh",
        gmsh_inp=case_dir / f"{stem}_mesh.inp",
        ccx_mesh_inp=case_dir / f"{stem}_ccx_mesh.inp",
        ccx_modal_inp=case_dir / f"{stem}_ccx_modal.inp",
        ccx_dat=case_dir / f"{stem}_ccx_modal.dat",
        ccx_frd=case_dir / f"{stem}_ccx_modal.frd",
        ccx_stdout=case_dir / f"{stem}_ccx_modal.stdout.txt",
        ccx_stderr=case_dir / f"{stem}_ccx_modal.stderr.txt",
    )

    audit = single.audit_tools()
    messages: list[str] = []
    actual_3d_run = False
    parse_result: object | None = None
    parse_message = "mode-shape parsing not attempted"
    omegas: tuple[float, ...] = ()
    fem_status = "not_run"
    fem_warning = "run without --run-3d-fem"

    if args.run_3d_fem:
        single.write_gmsh_geo(float(args.epsilon), paths)
        messages.append(f"wrote Gmsh geometry: {paths.geo.relative_to(REPO_ROOT)}")
        if not audit.gmsh_exe:
            fem_warning = "Gmsh executable not found"
            messages.append(fem_warning)
        elif not audit.calculix_ccx:
            fem_warning = "CalculiX executable not found"
            messages.append(fem_warning)
        else:
            mesh_ok, mesh_message, gmsh_messages = single.generate_mesh_with_gmsh_cli(paths, audit.gmsh_exe)
            messages.append(mesh_message)
            messages.extend(gmsh_messages)
            if mesh_ok and paths.gmsh_inp.exists():
                ccx_input = single.write_calculix_inputs(paths, float(args.epsilon))
                messages.append(
                    "wrote CalculiX input with "
                    f"{ccx_input.left_fixed_node_count} left fixed nodes and "
                    f"{ccx_input.right_fixed_node_count} right fixed nodes"
                )
                result = single.run_calculix(paths, audit.calculix_ccx, float(args.epsilon))
                actual_3d_run = True
                fem_status = "computed" if result.success and result.parsed_omegas else "run_failed_or_unparsed"
                fem_warning = "" if result.success and result.parsed_omegas else result.message
                messages.append(result.message)
                omegas = tuple(result.parsed_omegas)
                if result.success and result.parsed_omegas:
                    _metrics, parse_result = single.compute_mode_metrics_for_case(
                        epsilon=float(args.epsilon),
                        paths=paths,
                        omegas=result.parsed_omegas,
                    )
                    parse_message = getattr(parse_result, "message", "mode-shape parse result unavailable")
                    messages.append(parse_message)
            else:
                fem_warning = mesh_message
    else:
        messages.append("3D FEM not run because --run-3d-fem was not supplied.")

    fem_rows = build_raw_fem_rows(
        single=single,
        omegas=omegas,
        parse_result=parse_result,
        epsilon=float(args.epsilon),
        status=fem_status,
        warning=fem_warning,
    )
    match_rows = build_sorted_match_rows(analytic_rows, fem_rows)
    nearest_rows = compare_module().nearest_frequency_match_rows(analytic_rows, fem_rows)

    analytic_path = output_dir / "analytic_union.csv"
    fem_path = output_dir / "fem_3d_raw_modes.csv"
    match_path = output_dir / "analytic_vs_3d_fem_matched.csv"
    nearest_path = output_dir / "nearest_match_analytic_vs_3d_fem.csv"
    report_path = output_dir / "report.md"
    write_csv(
        analytic_path,
        analytic_rows,
        [
            "case_id",
            "beta_deg",
            "mu",
            "eta",
            "epsilon",
            "poisson",
            "analytic_full_index",
            "analytic_subsystem",
            "analytic_subsystem_index",
            "Lambda",
            "notes",
        ],
    )
    write_csv(
        fem_path,
        fem_rows,
        [
            "case_id",
            "fem_mode_index",
            "Omega_3d_fem",
            "Lambda_3d_fem",
            "fem_mode_character",
            "fem_in_plane_fraction",
            "fem_out_of_plane_fraction",
            "status",
            "warning",
            "notes",
        ],
    )
    write_csv(
        match_path,
        match_rows,
        [
            "case_id",
            "analytic_full_index",
            "analytic_subsystem",
            "analytic_subsystem_index",
            "Lambda_analytic",
            "fem_mode_index",
            "Lambda_3d_fem",
            "abs_error",
            "rel_error",
            "fem_mode_character",
            "fem_in_plane_fraction",
            "fem_out_of_plane_fraction",
            "match_method",
            "warning",
            "notes",
        ],
    )
    write_csv(
        nearest_path,
        nearest_rows,
        [
            "case_id",
            "fem_mode_index",
            "Lambda_3d_fem",
            "nearest_analytic_full_index",
            "nearest_analytic_subsystem",
            "nearest_analytic_subsystem_index",
            "Lambda_analytic_nearest",
            "abs_error",
            "rel_error",
            "ambiguity_count_within_tolerance",
            "warning",
            "notes",
        ],
    )
    generated_files = [analytic_path, fem_path, match_path, nearest_path]
    write_report(
        path=report_path,
        args=args,
        analytic_rows=analytic_rows,
        fem_rows=fem_rows,
        match_rows=match_rows,
        nearest_rows=nearest_rows,
        messages=messages,
        analytic_warnings=analytic_warnings,
        generated_files=generated_files,
        actual_3d_run=actual_3d_run,
        workflow_choice=WORKFLOW,
        radius=radius,
        mesh_size=mesh_size,
        gmsh_exe=audit.gmsh_exe,
        ccx_exe=audit.calculix_ccx,
        parse_message=parse_message,
    )
    print(f"saved analytic union CSV: {analytic_path}")
    print(f"saved raw FEM CSV: {fem_path}")
    print(f"saved matched comparison CSV: {match_path}")
    print(f"saved nearest-frequency CSV: {nearest_path}")
    print(f"saved report: {report_path}")
    print(f"actual 3D FEM run: {'yes' if actual_3d_run else 'no'}")
    print(f"workflow: {WORKFLOW}")
    print(f"radius: {radius:.16g}")
    print(f"mesh size: {mesh_size:.16g}")
    return {
        "analytic_path": analytic_path,
        "fem_path": fem_path,
        "match_path": match_path,
        "report_path": report_path,
        "actual_3d_run": actual_3d_run,
    }


def main(argv: list[str] | None = None) -> int:
    run(parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
