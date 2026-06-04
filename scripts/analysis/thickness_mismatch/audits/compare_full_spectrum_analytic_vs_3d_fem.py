from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import sys
from typing import Sequence

import numpy as np


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
DEFAULT_N_ANALYTIC_ROOTS_PER_SUBSYSTEM = 12
DEFAULT_N_FEM_MODES = 40
DEFAULT_OUTPUT_DIR = Path("results") / "full_spectrum_analytic_vs_3d_fem"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "full_spectrum_analytic_vs_3d_fem"
ROOT_SCAN_STEP = 0.02
ROOT_LMAX = 35.0
NEAREST_AMBIGUITY_REL_TOL = 1.0e-3


@dataclass(frozen=True)
class SpectrumCase:
    case_id: str
    beta_deg: float
    mu: float
    eta: float


@dataclass(frozen=True)
class Args:
    epsilon: float
    poisson: float
    cases_arg: str
    cases: tuple[SpectrumCase, ...]
    n_analytic_roots_per_subsystem: int
    n_fem_modes: int
    output_dir: Path
    smoke: bool
    skip_3d_fem: bool
    reuse_existing_fem_results: bool
    force: bool


DEFAULT_CASES = (
    SpectrumCase("case_1_beta0_mu0_eta0", beta_deg=0.0, mu=0.0, eta=0.0),
    SpectrumCase("case_2_beta15_mu0p3_eta0", beta_deg=15.0, mu=0.3, eta=0.0),
    SpectrumCase("case_3_beta15_mu0p3_eta0p1", beta_deg=15.0, mu=0.3, eta=0.1),
    SpectrumCase("case_4_beta45_mu0p6_eta0p5", beta_deg=45.0, mu=0.6, eta=0.5),
)
SMOKE_CASES = (SpectrumCase("smoke_beta0_mu0_eta0", beta_deg=0.0, mu=0.0, eta=0.0),)


def _repo_output_dir(path: Path) -> Path:
    path_obj = Path(path)
    return path_obj if path_obj.is_absolute() else REPO_ROOT / path_obj


def _fmt(value: float) -> str:
    value_f = float(value)
    if not isfinite(value_f):
        return "nan"
    return f"{value_f:.16e}"


def _read_case_csv(path: Path) -> tuple[SpectrumCase, ...]:
    cases: list[SpectrumCase] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for index, row in enumerate(reader, start=1):
            cases.append(
                SpectrumCase(
                    case_id=row.get("case_id") or f"case_{index}",
                    beta_deg=float(row["beta_deg"]),
                    mu=float(row["mu"]),
                    eta=float(row["eta"]),
                )
            )
    if not cases:
        raise ValueError(f"No cases found in {path}")
    return tuple(cases)


def parse_cases(cases_arg: str, *, smoke: bool) -> tuple[SpectrumCase, ...]:
    if smoke:
        return SMOKE_CASES
    if cases_arg == "default":
        return DEFAULT_CASES
    if cases_arg == "smoke":
        return SMOKE_CASES
    path = Path(cases_arg)
    if not path.is_absolute():
        path = REPO_ROOT / path
    return _read_case_csv(path)


def parse_args(argv: list[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(description="Experimental analytic full-spectrum vs 3D FEM audit scaffold.")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--poisson", type=float, default=DEFAULT_POISSON)
    parser.add_argument("--cases", default="default")
    parser.add_argument("--n-analytic-roots-per-subsystem", type=int, default=DEFAULT_N_ANALYTIC_ROOTS_PER_SUBSYSTEM)
    parser.add_argument("--n-fem-modes", type=int, default=DEFAULT_N_FEM_MODES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--skip-3d-fem", action="store_true")
    parser.add_argument("--reuse-existing-fem-results", action="store_true")
    parser.add_argument("--force", action="store_true")
    ns = parser.parse_args(argv)

    output_dir = _repo_output_dir(Path(ns.output_dir))
    n_analytic = int(ns.n_analytic_roots_per_subsystem)
    n_fem_modes = int(ns.n_fem_modes)
    cases_arg = str(ns.cases)
    if ns.smoke:
        output_dir = _repo_output_dir(SMOKE_OUTPUT_DIR)
        n_analytic = 6
        n_fem_modes = 12
        cases_arg = "smoke"
    args = Args(
        epsilon=float(ns.epsilon),
        poisson=float(ns.poisson),
        cases_arg=cases_arg,
        cases=parse_cases(cases_arg, smoke=bool(ns.smoke)),
        n_analytic_roots_per_subsystem=n_analytic,
        n_fem_modes=n_fem_modes,
        output_dir=output_dir,
        smoke=bool(ns.smoke),
        skip_3d_fem=bool(ns.skip_3d_fem),
        reuse_existing_fem_results=bool(ns.reuse_existing_fem_results),
        force=bool(ns.force),
    )
    validate_args(args)
    return args


def validate_args(args: Args) -> None:
    if not isfinite(args.epsilon) or args.epsilon <= 0.0:
        raise ValueError("epsilon must be positive and finite.")
    if not isfinite(args.poisson) or args.poisson <= -1.0:
        raise ValueError("poisson must be finite and greater than -1.")
    if args.n_analytic_roots_per_subsystem <= 0:
        raise ValueError("n_analytic_roots_per_subsystem must be positive.")
    if args.n_fem_modes <= 0:
        raise ValueError("n_fem_modes must be positive.")
    for case in args.cases:
        if not (-1.0 < case.mu < 1.0):
            raise ValueError("all mu values must lie inside (-1, 1).")
        if not isfinite(case.beta_deg):
            raise ValueError("all beta values must be finite.")


def build_analytic_full_spectrum_union(
    in_plane_roots: Sequence[float],
    out_of_plane_roots: Sequence[float],
) -> list[dict[str, float | int | str]]:
    """Return a sorted union of in-plane and out-of-plane roots.

    Close or equal values are intentionally not deduplicated because the full
    ideal spectrum can contain distinct modes at the same Lambda.
    """

    rows: list[dict[str, float | int | str]] = []
    for subsystem, roots in (("in_plane", in_plane_roots), ("out_of_plane", out_of_plane_roots)):
        for index, value in enumerate(roots, start=1):
            value_f = float(value)
            if isfinite(value_f):
                rows.append(
                    {
                        "analytic_subsystem": subsystem,
                        "analytic_subsystem_index": index,
                        "Lambda": value_f,
                    }
                )
    rows.sort(key=lambda row: (float(row["Lambda"]), str(row["analytic_subsystem"]), int(row["analytic_subsystem_index"])))
    for full_index, row in enumerate(rows, start=1):
        row["analytic_full_index"] = full_index
    return rows


def classify_3d_displacement_mode(
    displacements: np.ndarray,
    threshold: float = 0.8,
) -> tuple[str, float, float]:
    """Classify a 3D mode from displacement components.

    The model plane is the xy-plane and out-of-plane is z.
    """

    values = np.asarray(displacements, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError("displacements must have shape (n_nodes, 3).")
    total = float(np.sum(values**2))
    if total <= 0.0:
        return "mixed_or_unclassified", float("nan"), float("nan")
    in_plane = float(np.sum(values[:, :2] ** 2) / total)
    out_of_plane = float(np.sum(values[:, 2] ** 2) / total)
    if in_plane >= float(threshold):
        return "mostly_in_plane", in_plane, out_of_plane
    if out_of_plane >= float(threshold):
        return "mostly_out_of_plane", in_plane, out_of_plane
    return "mixed_or_unclassified", in_plane, out_of_plane


def existing_3d_workflow_summary() -> list[str]:
    candidates = [
        "scripts/analysis/solid_fem_coupled_equal_rods_point_joint.py",
        "scripts/analysis/solid_fem_coupled_equal_rods.py",
        "scripts/analysis/solid_fem_single_rod_fixed_fixed.py",
        "scripts/analysis/audit_3d_fem_point_joint_mu0_eps0p01.py",
        "scripts/analysis/audit_3d_fem_joint_coupling_models.py",
    ]
    found = [path for path in candidates if (REPO_ROOT / path).exists()]
    return found


def analytic_union_rows(args: Args) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in args.cases:
        beta_rad = float(np.deg2rad(case.beta_deg))
        in_plane = find_first_n_roots_eta(
            beta=beta_rad,
            mu=case.mu,
            epsilon=args.epsilon,
            eta=case.eta,
            n_roots=args.n_analytic_roots_per_subsystem,
            Lmax0=ROOT_LMAX,
            scan_step=ROOT_SCAN_STEP,
        )
        out = find_first_n_roots_out_of_plane_with_warnings(
            beta=beta_rad,
            mu=case.mu,
            epsilon=args.epsilon,
            eta=case.eta,
            poisson=args.poisson,
            n_roots=args.n_analytic_roots_per_subsystem,
            lambda_max=ROOT_LMAX,
            scan_step=ROOT_SCAN_STEP,
        )
        union = build_analytic_full_spectrum_union(in_plane, out.roots)
        for row in union:
            rows.append(
                {
                    "case_id": case.case_id,
                    "beta_deg": case.beta_deg,
                    "mu": case.mu,
                    "eta": case.eta,
                    "epsilon": args.epsilon,
                    "poisson": args.poisson,
                    "analytic_full_index": int(row["analytic_full_index"]),
                    "analytic_subsystem": row["analytic_subsystem"],
                    "analytic_subsystem_index": int(row["analytic_subsystem_index"]),
                    "Lambda": float(row["Lambda"]),
                    "notes": "out-of-plane root warning" if out.warnings and row["analytic_subsystem"] == "out_of_plane" else "ok",
                }
            )
    return rows


def raw_fem_rows(args: Args, status: str, warning: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in args.cases:
        rows.append(
            {
                "case_id": case.case_id,
                "beta_deg": case.beta_deg,
                "mu": case.mu,
                "eta": case.eta,
                "epsilon": args.epsilon,
                "poisson": args.poisson,
                "fem_mode_index": "",
                "Lambda_3d_fem": float("nan"),
                "fem_mode_character": "classification_unavailable",
                "fem_in_plane_fraction": float("nan"),
                "fem_out_of_plane_fraction": float("nan"),
                "status": status,
                "warning": warning,
                "notes": "3D FEM not run in this diagnostic invocation",
            }
        )
    return rows


def comparison_rows(analytic_rows: list[dict[str, object]], fem_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    fem_by_case: dict[str, list[dict[str, object]]] = {}
    for row in fem_rows:
        if isfinite(float(row["Lambda_3d_fem"])):
            fem_by_case.setdefault(str(row["case_id"]), []).append(row)

    rows: list[dict[str, object]] = []
    for analytic in analytic_rows:
        case_fem = sorted(fem_by_case.get(str(analytic["case_id"]), []), key=lambda row: float(row["Lambda_3d_fem"]))
        fem = None
        if case_fem:
            fem = min(case_fem, key=lambda row: abs(float(row["Lambda_3d_fem"]) - float(analytic["Lambda"])))
        if fem is None:
            rows.append(
                {
                    "case_id": analytic["case_id"],
                    "beta_deg": analytic["beta_deg"],
                    "mu": analytic["mu"],
                    "eta": analytic["eta"],
                    "epsilon": analytic["epsilon"],
                    "poisson": analytic["poisson"],
                    "analytic_full_index": analytic["analytic_full_index"],
                    "analytic_subsystem": analytic["analytic_subsystem"],
                    "analytic_subsystem_index": analytic["analytic_subsystem_index"],
                    "Lambda_analytic": analytic["Lambda"],
                    "fem_mode_index": "",
                    "Lambda_3d_fem": float("nan"),
                    "abs_error": float("nan"),
                    "rel_error": float("nan"),
                    "fem_mode_character": "classification_unavailable",
                    "fem_in_plane_fraction": float("nan"),
                    "fem_out_of_plane_fraction": float("nan"),
                    "match_method": "not_run",
                    "warning": "3D FEM skipped or unavailable",
                    "notes": "analytic-only row",
                }
            )
        else:
            abs_error = abs(float(fem["Lambda_3d_fem"]) - float(analytic["Lambda"]))
            rel_error = abs_error / abs(float(analytic["Lambda"]))
            rows.append(
                {
                    "case_id": analytic["case_id"],
                    "beta_deg": analytic["beta_deg"],
                    "mu": analytic["mu"],
                    "eta": analytic["eta"],
                    "epsilon": analytic["epsilon"],
                    "poisson": analytic["poisson"],
                    "analytic_full_index": analytic["analytic_full_index"],
                    "analytic_subsystem": analytic["analytic_subsystem"],
                    "analytic_subsystem_index": analytic["analytic_subsystem_index"],
                    "Lambda_analytic": analytic["Lambda"],
                    "fem_mode_index": fem["fem_mode_index"],
                    "Lambda_3d_fem": fem["Lambda_3d_fem"],
                    "abs_error": abs_error,
                    "rel_error": rel_error,
                    "fem_mode_character": fem["fem_mode_character"],
                    "fem_in_plane_fraction": fem["fem_in_plane_fraction"],
                    "fem_out_of_plane_fraction": fem["fem_out_of_plane_fraction"],
                    "match_method": "nearest_frequency",
                    "warning": "",
                    "notes": "exploratory frequency-only match",
                }
            )
    return rows


def nearest_frequency_match_rows(
    analytic_rows: Sequence[dict[str, object]],
    fem_rows: Sequence[dict[str, object]],
    *,
    ambiguity_rel_tol: float = NEAREST_AMBIGUITY_REL_TOL,
) -> list[dict[str, object]]:
    """Greedy one-to-one nearest-frequency diagnostic matches.

    Analytic duplicate roots are kept as separate rows. Each FEM row can match
    one analytic row and each analytic row is used at most once within a case.
    """

    analytic_by_case: dict[str, list[dict[str, object]]] = {}
    for row in analytic_rows:
        try:
            Lambda = float(row["Lambda"])
        except (KeyError, TypeError, ValueError):
            continue
        if isfinite(Lambda):
            analytic_by_case.setdefault(str(row["case_id"]), []).append(row)
    for rows in analytic_by_case.values():
        rows.sort(key=lambda row: (float(row["Lambda"]), int(row["analytic_full_index"])))

    fem_by_case: dict[str, list[dict[str, object]]] = {}
    for row in fem_rows:
        try:
            Lambda = float(row["Lambda_3d_fem"])
        except (KeyError, TypeError, ValueError):
            continue
        if isfinite(Lambda):
            fem_by_case.setdefault(str(row["case_id"]), []).append(row)
    for rows in fem_by_case.values():
        rows.sort(key=lambda row: int(row["fem_mode_index"]) if str(row["fem_mode_index"]).isdigit() else 10**9)

    output: list[dict[str, object]] = []
    for case_id, case_fem_rows in sorted(fem_by_case.items()):
        available = list(analytic_by_case.get(case_id, []))
        all_analytic = list(available)
        for fem in case_fem_rows:
            Lambda_fem = float(fem["Lambda_3d_fem"])
            if not available:
                output.append(
                    {
                        "case_id": case_id,
                        "fem_mode_index": fem["fem_mode_index"],
                        "Lambda_3d_fem": Lambda_fem,
                        "nearest_analytic_full_index": "",
                        "nearest_analytic_subsystem": "",
                        "nearest_analytic_subsystem_index": "",
                        "Lambda_analytic_nearest": float("nan"),
                        "abs_error": float("nan"),
                        "rel_error": float("nan"),
                        "ambiguity_count_within_tolerance": 0,
                        "warning": "no unused analytic root available",
                        "notes": "nearest-frequency diagnostic; not branch identity",
                    }
                )
                continue
            nearest = min(available, key=lambda row: abs(float(row["Lambda"]) - Lambda_fem))
            available.remove(nearest)
            Lambda_analytic = float(nearest["Lambda"])
            abs_error = abs(Lambda_fem - Lambda_analytic)
            rel_error = abs_error / abs(Lambda_analytic) if Lambda_analytic != 0.0 else float("nan")
            ambiguity_count = 0
            for analytic in all_analytic:
                candidate = float(analytic["Lambda"])
                scale = max(abs(candidate), 1.0e-12)
                if abs(Lambda_fem - candidate) / scale <= float(ambiguity_rel_tol):
                    ambiguity_count += 1
            warning = "ambiguous duplicate/clustered analytic roots" if ambiguity_count > 1 else ""
            output.append(
                {
                    "case_id": case_id,
                    "fem_mode_index": fem["fem_mode_index"],
                    "Lambda_3d_fem": Lambda_fem,
                    "nearest_analytic_full_index": nearest["analytic_full_index"],
                    "nearest_analytic_subsystem": nearest["analytic_subsystem"],
                    "nearest_analytic_subsystem_index": nearest["analytic_subsystem_index"],
                    "Lambda_analytic_nearest": Lambda_analytic,
                    "abs_error": abs_error,
                    "rel_error": rel_error,
                    "ambiguity_count_within_tolerance": ambiguity_count,
                    "warning": warning,
                    "notes": "nearest-frequency diagnostic; not branch identity",
                }
            )
    return output


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out = dict(row)
            for key, value in list(out.items()):
                if isinstance(value, float):
                    out[key] = _fmt(value)
            writer.writerow(out)


def write_report(
    args: Args,
    analytic_rows: list[dict[str, object]],
    fem_rows: list[dict[str, object]],
    comparison: list[dict[str, object]],
    output_paths: list[Path],
    workflow_summary: list[str],
    fem_status: str,
    nearest_rows: list[dict[str, object]] | None = None,
) -> Path:
    report_path = args.output_dir / "full_spectrum_analytic_vs_3d_fem_report.md"
    cases = ", ".join(case.case_id for case in args.cases)
    lines = [
        "# Full-Spectrum Analytic vs 3D FEM Diagnostic",
        "",
        "This is an exploratory second-level comparison scaffold. The analytic",
        "full EB-level spectrum is the sorted union of in-plane EB roots and",
        "out-of-plane EB + Saint-Venant torsion roots.",
        "",
        "3D FEM includes effects outside the EB determinants, including finite",
        "joint geometry, shear/rotary inertia, possible warping, Poisson effects,",
        "mesh sensitivity, and mode mixing.",
        "",
        "## Parameters",
        "",
        f"- epsilon: {args.epsilon:g}",
        f"- poisson: {args.poisson:g}",
        f"- cases: {cases}",
        f"- analytic roots per subsystem: {args.n_analytic_roots_per_subsystem}",
        f"- requested FEM modes: {args.n_fem_modes}",
        "",
        "## Existing 3D FEM Workflows Found",
        "",
    ]
    lines.extend([f"- `{path}`" for path in workflow_summary] or ["- none found"])
    lines.extend(
        [
            "",
            "## 3D FEM Status",
            "",
            f"- actual 3D FEM run: {fem_status}",
            "- no production FEM/Gmsh/CalculiX workflow was modified.",
            "- no old FEM outputs were overwritten.",
            "",
            "The existing workflows are substantial diagnostic scripts, not a small",
            "arbitrary `(beta, mu, eta)` API for full-spectrum comparison. This",
            "wrapper therefore starts with an analytic union and an explicit",
            "skip/reuse boundary. Actual 3D runs should be launched only after",
            "choosing the appropriate existing point-joint or solid workflow and",
            "confirming runtime, mesh settings, radius scaling, and output paths.",
            "",
            "## Analytic Union",
            "",
            f"- analytic union rows written: {len(analytic_rows)}",
            f"- comparison rows written: {len(comparison)}",
            f"- nearest-frequency match rows written: {0 if nearest_rows is None else len(nearest_rows)}",
            "",
            "## Mode Classification",
            "",
            "- raw 3D mode classification is unavailable because 3D FEM was not run",
            "  in this invocation.",
            "- helper `classify_3d_displacement_mode` is available for parsed",
            "  displacement vectors when a reusable 3D result path is selected.",
            "",
            "## Limitations",
            "",
            "- frequency-only nearest matching is not used unless FEM rows exist.",
            "- no claims are made about determinant errors, crossings, or branch",
            "  identity.",
            "- this is not a replacement for the validated 1D EB+torsion FEM check.",
            "",
            "## Generated Files",
            "",
        ]
    )
    for path in [*output_paths, report_path]:
        lines.append(f"- `{path.relative_to(REPO_ROOT)}`")
    lines.append("")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path


def run(args: Args) -> dict[str, object]:
    workflow_summary = existing_3d_workflow_summary()
    actual_3d_run = False
    fem_status = "skipped by --skip-3d-fem" if args.skip_3d_fem else "not run by scaffold"
    if args.reuse_existing_fem_results:
        fem_status = "reuse requested, but no standardized reusable full-spectrum result parser is implemented yet"
    if not args.skip_3d_fem and not args.reuse_existing_fem_results:
        fem_status = "not run; use --skip-3d-fem for analytic-only smoke or extend wrapper after selecting workflow"

    analytic = analytic_union_rows(args)
    fem = raw_fem_rows(args, status=fem_status, warning="actual 3D FEM not run")
    comparison = comparison_rows(analytic, fem)
    nearest = nearest_frequency_match_rows(analytic, fem)

    output_paths: list[Path] = []
    union_path = args.output_dir / "full_spectrum_analytic_union.csv"
    raw_fem_path = args.output_dir / "full_spectrum_3d_fem_raw_modes.csv"
    comparison_path = args.output_dir / "full_spectrum_analytic_vs_3d_fem_comparison.csv"
    nearest_path = args.output_dir / "nearest_match_analytic_vs_3d_fem.csv"
    write_csv(
        union_path,
        analytic,
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
        raw_fem_path,
        fem,
        [
            "case_id",
            "beta_deg",
            "mu",
            "eta",
            "epsilon",
            "poisson",
            "fem_mode_index",
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
        comparison_path,
        comparison,
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
        nearest,
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
    output_paths.extend([union_path, raw_fem_path, comparison_path, nearest_path])
    report_path = write_report(args, analytic, fem, comparison, output_paths, workflow_summary, fem_status, nearest)
    output_paths.append(report_path)

    print(f"saved analytic union CSV: {union_path}")
    print(f"saved raw FEM CSV: {raw_fem_path}")
    print(f"saved comparison CSV: {comparison_path}")
    print(f"saved nearest-frequency CSV: {nearest_path}")
    print(f"saved report: {report_path}")
    print(f"existing 3D workflows found: {len(workflow_summary)}")
    print(f"actual 3D FEM run: {'yes' if actual_3d_run else 'no'}")
    print(f"3D FEM status: {fem_status}")
    return {
        "output_paths": output_paths,
        "workflow_summary": workflow_summary,
        "actual_3d_run": actual_3d_run,
        "fem_status": fem_status,
        "analytic_rows": len(analytic),
        "comparison_rows": len(comparison),
    }


def main(argv: list[str] | None = None) -> dict[str, object]:
    return run(parse_args(argv))


if __name__ == "__main__":
    main()
