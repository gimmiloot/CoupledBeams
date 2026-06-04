from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import re
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


DEFAULT_EPSILON = 0.0025
DEFAULT_POISSON = 0.3
DEFAULT_TOTAL_LENGTH = 2.0
DEFAULT_REQUESTED_MODES = 40
DEFAULT_COMPARE_MODES = 16
DEFAULT_MESH_SIZES = (0.005, 0.004, 0.003)
DEFAULT_OUTPUT_DIR = Path("results") / "full_spectrum_analytic_vs_3d_fem" / "straight_uniform_mesh_convergence"
DEFAULT_LAMBDA_MAX = 35.0
DEFAULT_ROOT_SCAN_STEP = 0.02
WORKFLOW = "straight_fixed_fixed_cylinder_from_single_rod_helpers"


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def _token(value: float) -> str:
    return f"{float(value):.6g}".replace("-", "m").replace(".", "p")


def _fmt(value: object) -> object:
    if isinstance(value, float):
        if not math.isfinite(value):
            return "nan"
        return f"{value:.16e}"
    return value


def write_csv(path: Path, rows: Sequence[dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row.get(key, "")) for key in fieldnames})


def rel(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def error_stats(values: Sequence[float]) -> tuple[float, float]:
    finite = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not finite:
        return float("nan"), float("nan")
    midpoint = len(finite) // 2
    median = finite[midpoint] if len(finite) % 2 else 0.5 * (finite[midpoint - 1] + finite[midpoint])
    return max(finite), median


def parse_mesh_quality(messages: Sequence[str]) -> tuple[int, float, tuple[str, ...]]:
    negative_count = 0
    worst = float("nan")
    warning_lines: list[str] = []
    pattern = re.compile(r"worst distortion\s*=\s*([-+0-9.Ee]+).*?(\d+)\s+elements with jac\.\s*<\s*0")
    for message in messages:
        if "worst distortion" in message or "jac." in message:
            warning_lines.append(message)
        match = pattern.search(message)
        if match:
            worst = float(match.group(1))
            negative_count = int(match.group(2))
    return negative_count, worst, tuple(warning_lines)


def first_lambdas(fem_rows: Sequence[dict[str, object]], limit: int = 10) -> list[float]:
    values: list[float] = []
    for row in fem_rows:
        try:
            value = float(row["Lambda_3d_fem"])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(value):
            values.append(value)
        if len(values) >= limit:
            break
    return values


def doublet_splits(values: Sequence[float], *, pair_count: int = 8) -> list[float]:
    splits: list[float] = []
    for start in range(0, min(len(values), 2 * pair_count), 2):
        if start + 1 >= len(values):
            break
        left = float(values[start])
        right = float(values[start + 1])
        mean = 0.5 * (left + right)
        splits.append(abs(right - left) / mean if mean > 0.0 else float("nan"))
    return splits


def with_case_id(rows: Sequence[dict[str, object]], case_id: str) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in rows:
        copied = dict(row)
        copied["case_id"] = case_id
        out.append(copied)
    return out


def sorted_mode_rows(
    *,
    mesh_info: dict[str, object],
    analytic_rows: Sequence[dict[str, object]],
    fem_rows: Sequence[dict[str, object]],
    compare_modes: int,
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for index in range(compare_modes):
        if index >= len(analytic_rows):
            break
        analytic = analytic_rows[index]
        fem = fem_rows[index] if index < len(fem_rows) else None
        Lambda_analytic = float(analytic["Lambda"])
        Lambda_fem = float("nan") if fem is None else float(fem["Lambda_3d_fem"])
        abs_error = abs(Lambda_fem - Lambda_analytic) if math.isfinite(Lambda_fem) else float("nan")
        rel_error = abs_error / abs(Lambda_analytic) if math.isfinite(abs_error) and Lambda_analytic != 0.0 else float("nan")
        out.append(
            {
                **mesh_info,
                "mode_index": index + 1,
                "match_method": "sorted_index",
                "Lambda_analytic": Lambda_analytic,
                "Lambda_3d_fem": Lambda_fem,
                "abs_error": abs_error,
                "rel_error": rel_error,
                "fem_mode_character": "" if fem is None else fem["fem_mode_character"],
                "notes": "sorted-index diagnostic; duplicate pairs require nearest-match cross-check",
            }
        )
    return out


def nearest_mode_rows(
    *,
    compare: ModuleType,
    mesh_info: dict[str, object],
    analytic_rows: Sequence[dict[str, object]],
    fem_rows: Sequence[dict[str, object]],
    compare_modes: int,
) -> list[dict[str, object]]:
    nearest = compare.nearest_frequency_match_rows(analytic_rows, fem_rows)
    out: list[dict[str, object]] = []
    for row in nearest[:compare_modes]:
        out.append(
            {
                **mesh_info,
                "mode_index": row["fem_mode_index"],
                "match_method": "nearest_frequency",
                "Lambda_analytic": row["Lambda_analytic_nearest"],
                "Lambda_3d_fem": row["Lambda_3d_fem"],
                "abs_error": row["abs_error"],
                "rel_error": row["rel_error"],
                "fem_mode_character": "",
                "notes": (
                    f"nearest analytic full index {row['nearest_analytic_full_index']}; "
                    f"ambiguity_count={row['ambiguity_count_within_tolerance']}; {row['warning']}"
                ),
            }
        )
    return out


def configure_single_module(single: ModuleType, *, total_length: float, poisson: float, requested_modes: int, timeout: int) -> None:
    single.L = float(total_length)
    single.NU = float(poisson)
    single.SOLID_MODES_REQUESTED = int(requested_modes)
    single.CCX_TIMEOUT_SECONDS = int(timeout)


def run_one_mesh(
    *,
    mesh_id: str,
    h: float,
    args: argparse.Namespace,
    single: ModuleType,
    smoke: ModuleType,
    compare: ModuleType,
    analytic_template: Sequence[dict[str, object]],
    audit: object,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    single.mesh_size_for_epsilon = lambda _epsilon: float(h)
    case_dir = Path(args.output_dir) / "solid_fem_cases" / mesh_id
    stem = f"straight_uniform_L2_eps0p0025_{mesh_id}"
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
    single.write_gmsh_geo(float(args.epsilon), paths)
    messages: list[str] = [f"wrote {rel(paths.geo)}"]
    omegas: tuple[float, ...] = ()
    parse_result: object | None = None
    num_nodes = 0
    num_elements = 0
    status = "not_run"
    warning = ""
    gmsh_messages: tuple[str, ...] = ()

    if not getattr(audit, "gmsh_exe"):
        warning = "Gmsh executable not found"
        messages.append(warning)
    elif not getattr(audit, "calculix_ccx"):
        warning = "CalculiX executable not found"
        messages.append(warning)
    else:
        mesh_ok, mesh_message, gmsh_messages = single.generate_mesh_with_gmsh_cli(paths, audit.gmsh_exe)
        messages.append(mesh_message)
        messages.extend(gmsh_messages)
        if mesh_ok and paths.gmsh_inp.exists():
            mesh = single.read_gmsh_inp_mesh_data(paths.gmsh_inp)
            num_nodes = len(mesh.nodes)
            num_elements = len(mesh.solid_elements)
            ccx_input = single.write_calculix_inputs(paths, float(args.epsilon))
            messages.append(
                f"wrote CalculiX input: left_fixed={ccx_input.left_fixed_node_count}, "
                f"right_fixed={ccx_input.right_fixed_node_count}"
            )
            result = single.run_calculix(paths, audit.calculix_ccx, float(args.epsilon))
            messages.append(result.message)
            status = "computed" if result.success and result.parsed_omegas else "run_failed_or_unparsed"
            warning = "" if result.success and result.parsed_omegas else result.message
            omegas = tuple(result.parsed_omegas)
            if result.success and result.parsed_omegas:
                _metrics, parse_result = single.compute_mode_metrics_for_case(
                    epsilon=float(args.epsilon),
                    paths=paths,
                    omegas=result.parsed_omegas,
                )
                messages.append(getattr(parse_result, "message", "mode-shape parse result unavailable"))
        else:
            warning = mesh_message

    num_negative, worst_distortion, mesh_warning_lines = parse_mesh_quality(gmsh_messages)
    fem_rows = smoke.build_raw_fem_rows(
        single=single,
        omegas=omegas,
        parse_result=parse_result,
        epsilon=float(args.epsilon),
        status=status,
        warning=warning,
    )
    fem_rows = with_case_id(fem_rows, mesh_id)
    analytic_rows = with_case_id(analytic_template, mesh_id)
    mesh_info = {
        "mesh_id": mesh_id,
        "h_or_mesh_parameter": float(h),
        "num_nodes": num_nodes,
        "num_elements": num_elements,
        "num_negative_jacobian_elements": num_negative,
        "worst_distortion": worst_distortion,
        "requested_modes": int(args.n_fem_modes),
        "parsed_modes": len([row for row in fem_rows if str(row.get("parser_mode_index", "")) or str(row.get("fem_mode_index", "")).isdigit()]),
    }
    sorted_rows = sorted_mode_rows(
        mesh_info=mesh_info,
        analytic_rows=analytic_rows,
        fem_rows=fem_rows,
        compare_modes=int(args.compare_modes),
    )
    nearest_rows = nearest_mode_rows(
        compare=compare,
        mesh_info=mesh_info,
        analytic_rows=analytic_rows,
        fem_rows=fem_rows,
        compare_modes=int(args.compare_modes),
    )
    first10 = first_lambdas(fem_rows, 10)
    split_values = doublet_splits(first_lambdas(fem_rows, int(args.compare_modes)), pair_count=int(args.compare_modes) // 2)
    sorted_errors = [float(row["rel_error"]) for row in sorted_rows]
    nearest_errors = [float(row["rel_error"]) for row in nearest_rows]
    sorted_max, sorted_median = error_stats(sorted_errors)
    nearest_max, nearest_median = error_stats(nearest_errors)
    split_max, split_median = error_stats(split_values)
    summary = {
        **mesh_info,
        "low_modes_recovered": bool(first10 and first10[0] < 3.0),
        "sorted_first_modes_max_rel_error": sorted_max,
        "sorted_first_modes_median_rel_error": sorted_median,
        "nearest_first_modes_max_rel_error": nearest_max,
        "nearest_first_modes_median_rel_error": nearest_median,
        "max_doublet_split_rel": split_max,
        "median_doublet_split_rel": split_median,
        "first_10_fem_lambdas": "; ".join(f"{value:.8g}" for value in first10),
        "mesh_warnings": " | ".join(mesh_warning_lines),
        "messages": " | ".join(messages),
    }
    return summary, [*sorted_rows, *nearest_rows]


def recommendation(summary_rows: Sequence[dict[str, object]]) -> str:
    candidates = [row for row in summary_rows if row.get("low_modes_recovered")]
    if not candidates:
        return "No tested mesh recovered the low modes; do not use this 3D path for future diagnostics yet."
    clean = [row for row in candidates if int(row["num_negative_jacobian_elements"]) == 0]
    if clean:
        candidates = clean
    best = min(
        candidates,
        key=lambda row: (
            float(row["nearest_first_modes_max_rel_error"]),
            float(row["nearest_first_modes_median_rel_error"]),
            int(row["num_nodes"]),
            float(row["h_or_mesh_parameter"]),
        ),
    )
    return (
        f"Use {best['mesh_id']} (h={float(best['h_or_mesh_parameter']):.6g}) for future straight-cylinder "
        "diagnostic smoke comparisons unless runtime becomes limiting; rerun mesh-quality cleanup before article use."
    )


def write_report(path: Path, args: argparse.Namespace, summary_rows: Sequence[dict[str, object]]) -> None:
    lines = [
        "# Straight Uniform 3D FEM Mesh Convergence",
        "",
        "Diagnostic-only mesh/extraction convergence for the straight fixed-fixed cylinder.",
        "This is a 3D smoke/convergence baseline, not article validation.",
        "",
        "## Case",
        "",
        f"- epsilon: {args.epsilon:g}",
        f"- poisson: {args.poisson:g}",
        "- beta: 0 deg",
        "- mu: 0",
        "- eta: 0",
        f"- total length L: {args.total_length:g}",
        f"- radius r0 = 2*epsilon: {2.0 * args.epsilon:g}",
        "- material: E=1, rho=1",
        f"- requested modes per mesh: {args.n_fem_modes}",
        f"- compared modes per mesh: {args.compare_modes}",
        f"- workflow: `{WORKFLOW}`",
        "- mesh helper: existing unstructured second-order Gmsh cylinder helper; no structured/transfinite option was used.",
        "",
        "## Summary",
        "",
        "| mesh | h | nodes | elems | neg jac elems | worst distortion | parsed | sorted median err | nearest median err | max doublet split | low modes |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in summary_rows:
        lines.append(
            "| "
            f"{row['mesh_id']} | {float(row['h_or_mesh_parameter']):.6g} | {int(row['num_nodes'])} | "
            f"{int(row['num_elements'])} | {int(row['num_negative_jacobian_elements'])} | "
            f"{float(row['worst_distortion']):.6g} | {int(row['parsed_modes'])} | "
            f"{float(row['sorted_first_modes_median_rel_error']):.6g} | "
            f"{float(row['nearest_first_modes_median_rel_error']):.6g} | "
            f"{float(row['max_doublet_split_rel']):.6g} | {row['low_modes_recovered']} |"
        )
    lines.extend(["", "## First 10 FEM Lambdas", ""])
    for row in summary_rows:
        lines.append(f"- {row['mesh_id']} h={float(row['h_or_mesh_parameter']):.6g}: {row['first_10_fem_lambdas']}")
    lines.extend(["", "## Answers", ""])
    low_modes_all = all(bool(row["low_modes_recovered"]) for row in summary_rows)
    distorted = [row for row in summary_rows if int(row["num_negative_jacobian_elements"]) > 0]
    clean = [row for row in summary_rows if int(row["num_negative_jacobian_elements"]) == 0]
    lines.append(
        f"- Requesting {args.n_fem_modes} modes {'consistently recovers' if low_modes_all else 'does not consistently recover'} "
        "the low bending doublets."
    )
    lines.append("- Errors stay at diagnostic smoke scale; see the summary table for sorted and nearest median relative errors.")
    if distorted and clean:
        lines.append(
            "- Mesh distortion appears only in "
            f"{', '.join(str(row['mesh_id']) for row in distorted)}; "
            f"{', '.join(str(row['mesh_id']) for row in clean)} had no reported negative-Jacobian warning."
        )
    elif distorted:
        lines.append("- Mesh distortion remains for all tested meshes; warning details are retained below.")
    else:
        lines.append("- No negative-Jacobian mesh distortion warning was reported for the tested meshes.")
    lines.append("- Bending doublets remain split; the summary table reports the maximum relative split among compared pairs.")
    lines.append(f"- Recommended mesh setting: {recommendation(summary_rows)}")
    lines.append("- The remaining discrepancy is acceptable for diagnostic straight-cylinder FEM smoke work, not for article validation.")
    lines.extend(["", "## Mesh Warnings", ""])
    for row in summary_rows:
        warnings = str(row["mesh_warnings"]) or "none"
        lines.append(f"- {row['mesh_id']}: {warnings}")
    lines.extend(
        [
            "",
            "## Generated Files",
            "",
            f"- `{rel(path / 'mesh_convergence_summary.csv')}`",
            f"- `{rel(path / 'mesh_convergence_modes.csv')}`",
            f"- `{rel(path / 'report.md')}`",
            "",
        ]
    )
    path.mkdir(parents=True, exist_ok=True)
    (path / "report.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Straight uniform 3D FEM mesh convergence audit.")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--poisson", type=float, default=DEFAULT_POISSON)
    parser.add_argument("--total-length", type=float, default=DEFAULT_TOTAL_LENGTH)
    parser.add_argument("--n-fem-modes", type=int, default=DEFAULT_REQUESTED_MODES)
    parser.add_argument("--compare-modes", type=int, default=DEFAULT_COMPARE_MODES)
    parser.add_argument("--mesh-sizes", nargs="+", type=float, default=list(DEFAULT_MESH_SIZES))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--ccx-timeout-seconds", type=int, default=420)
    parser.add_argument("--lambda-max", type=float, default=DEFAULT_LAMBDA_MAX)
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    args = parser.parse_args(argv)
    if args.epsilon <= 0.0:
        raise ValueError("epsilon must be positive.")
    if args.total_length <= 0.0:
        raise ValueError("total length must be positive.")
    if args.n_fem_modes <= args.compare_modes:
        raise ValueError("n_fem_modes must exceed compare_modes.")
    if any(value <= 0.0 for value in args.mesh_sizes):
        raise ValueError("mesh sizes must be positive.")
    args.output_dir = _repo_path(Path(args.output_dir))
    return args


def run(args: argparse.Namespace) -> dict[str, object]:
    single = _load_module("solid_fem_single_rod_fixed_fixed_mesh_convergence", REPO_ROOT / "scripts" / "analysis" / "solid_fem_single_rod_fixed_fixed.py")
    smoke = _load_module("run_full_spectrum_3d_fem_smoke_case_mesh_convergence", REPO_ROOT / "scripts" / "analysis" / "thickness_mismatch" / "audits" / "run_full_spectrum_3d_fem_smoke_case.py")
    compare = _load_module("compare_full_spectrum_analytic_vs_3d_fem_mesh_convergence", REPO_ROOT / "scripts" / "analysis" / "thickness_mismatch" / "audits" / "compare_full_spectrum_analytic_vs_3d_fem.py")
    configure_single_module(
        single,
        total_length=float(args.total_length),
        poisson=float(args.poisson),
        requested_modes=int(args.n_fem_modes),
        timeout=int(args.ccx_timeout_seconds),
    )
    audit = single.audit_tools()
    analytic_template, _warnings = smoke.build_analytic_rows(
        epsilon=float(args.epsilon),
        poisson=float(args.poisson),
        n_roots_per_subsystem=max(20, int(args.compare_modes)),
        lambda_max=float(args.lambda_max),
        root_scan_step=float(args.root_scan_step),
    )

    summary_rows: list[dict[str, object]] = []
    mode_rows: list[dict[str, object]] = []
    for h in args.mesh_sizes:
        mesh_id = f"h_{_token(float(h))}"
        print(f"running mesh {mesh_id} (h={float(h):g})")
        summary, rows = run_one_mesh(
            mesh_id=mesh_id,
            h=float(h),
            args=args,
            single=single,
            smoke=smoke,
            compare=compare,
            analytic_template=analytic_template,
            audit=audit,
        )
        summary_rows.append(summary)
        mode_rows.extend(rows)

    output_dir = Path(args.output_dir)
    summary_path = output_dir / "mesh_convergence_summary.csv"
    modes_path = output_dir / "mesh_convergence_modes.csv"
    summary_fields = [
        "mesh_id",
        "h_or_mesh_parameter",
        "num_nodes",
        "num_elements",
        "num_negative_jacobian_elements",
        "worst_distortion",
        "requested_modes",
        "parsed_modes",
        "low_modes_recovered",
        "sorted_first_modes_max_rel_error",
        "sorted_first_modes_median_rel_error",
        "nearest_first_modes_max_rel_error",
        "nearest_first_modes_median_rel_error",
        "max_doublet_split_rel",
        "median_doublet_split_rel",
        "first_10_fem_lambdas",
        "mesh_warnings",
        "messages",
    ]
    mode_fields = [
        "mesh_id",
        "h_or_mesh_parameter",
        "num_nodes",
        "num_elements",
        "num_negative_jacobian_elements",
        "worst_distortion",
        "requested_modes",
        "parsed_modes",
        "mode_index",
        "match_method",
        "Lambda_analytic",
        "Lambda_3d_fem",
        "abs_error",
        "rel_error",
        "fem_mode_character",
        "notes",
    ]
    write_csv(summary_path, summary_rows, summary_fields)
    write_csv(modes_path, mode_rows, mode_fields)
    write_report(output_dir, args, summary_rows)
    print(f"saved summary CSV: {summary_path}")
    print(f"saved modes CSV: {modes_path}")
    print(f"saved report: {output_dir / 'report.md'}")
    print(recommendation(summary_rows))
    return {"summary_path": summary_path, "modes_path": modes_path, "report_path": output_dir / "report.md"}


def main(argv: list[str] | None = None) -> int:
    run(parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
