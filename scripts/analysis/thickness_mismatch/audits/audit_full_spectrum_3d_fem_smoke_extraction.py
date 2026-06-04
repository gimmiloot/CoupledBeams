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


DEFAULT_SMOKE_DIR = Path("results") / "full_spectrum_analytic_vs_3d_fem" / "smoke_straight_uniform"
DEFAULT_EPSILON = 0.0025
DEFAULT_TOTAL_LENGTH = 2.0
FIRST_CLAMPED_CLAMPED_ALPHA = 4.730040744862704


def _load_module(module_name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


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


def read_csv(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def find_single(path: Path, pattern: str) -> Path | None:
    matches = sorted(path.glob(pattern))
    return matches[0] if matches else None


def frequency_card_lines(input_path: Path) -> list[str]:
    if not input_path.exists():
        return []
    lines = input_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    for index, line in enumerate(lines):
        if line.strip().upper().startswith("*FREQUENCY"):
            out = [line.strip()]
            for follow in lines[index + 1 :]:
                if follow.strip().startswith("*"):
                    break
                if follow.strip():
                    out.append(follow.strip())
            return out
    return []


def requested_mode_count(card_lines: Sequence[str]) -> int | None:
    for line in card_lines[1:]:
        first = line.split(",", 1)[0].strip()
        try:
            return int(first)
        except ValueError:
            continue
    return None


def has_frequency_shift_or_bounds(card_lines: Sequence[str]) -> bool:
    text = " ".join(card_lines).upper()
    markers = ("LOWER", "UPPER", "SHIFT", "STORAGE", "SOLVER=", "PERTURBATION")
    if any(marker in text for marker in markers):
        return True
    for line in card_lines[1:]:
        if len([part for part in line.split(",") if part.strip()]) > 1:
            return True
    return False


def parse_calculix_frequency_table(dat_path: Path) -> list[dict[str, float | int]]:
    if not dat_path.exists():
        return []
    rows: list[dict[str, float | int]] = []
    in_eigen_table = False
    saw_frequency_header = False
    for line in dat_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        upper = line.upper()
        if "E I G E N V A L U E" in upper and "O U T P U T" in upper:
            in_eigen_table = True
            saw_frequency_header = False
            continue
        if in_eigen_table and "FREQUENCY" in upper and "MODE" in upper:
            saw_frequency_header = True
            continue
        if in_eigen_table and "P A R T I C I P A T I O N" in upper:
            break
        if not in_eigen_table or not saw_frequency_header:
            continue
        parts = line.split()
        if len(parts) < 5 or not parts[0].isdigit():
            continue
        try:
            mode = int(parts[0])
            eigenvalue = float(parts[1].replace("D", "E").replace("d", "E"))
            angular = float(parts[2].replace("D", "E").replace("d", "E"))
            cyclic = float(parts[3].replace("D", "E").replace("d", "E"))
        except ValueError:
            continue
        if mode > 0 and eigenvalue > 0.0 and angular > 0.0:
            rows.append(
                {
                    "raw_mode_number": mode,
                    "eigenvalue": eigenvalue,
                    "angular_frequency": angular,
                    "cyclic_frequency": cyclic,
                }
            )
    return rows


def parsed_mode_map(fem_rows: Sequence[dict[str, object]]) -> dict[int, dict[str, object]]:
    out: dict[int, dict[str, object]] = {}
    for row in fem_rows:
        try:
            mode = int(str(row["fem_mode_index"]))
        except (KeyError, TypeError, ValueError):
            continue
        out[mode] = row
    return out


def sorted_match_mode_numbers(match_rows: Sequence[dict[str, object]]) -> set[int]:
    out: set[int] = set()
    for row in match_rows:
        try:
            out.add(int(str(row["fem_mode_index"])))
        except (KeyError, TypeError, ValueError):
            continue
    return out


def raw_frequency_audit_rows(
    raw_rows: Sequence[dict[str, float | int]],
    fem_rows: Sequence[dict[str, object]],
    match_rows: Sequence[dict[str, object]],
    *,
    epsilon: float,
) -> list[dict[str, object]]:
    parsed = parsed_mode_map(fem_rows)
    matched_modes = sorted_match_mode_numbers(match_rows)
    rows: list[dict[str, object]] = []
    for raw in raw_rows:
        mode = int(raw["raw_mode_number"])
        angular = float(raw["angular_frequency"])
        converted = math.sqrt(angular / float(epsilon))
        parsed_row = parsed.get(mode)
        parser_included = parsed_row is not None
        reason = "" if parser_included else "not present in fem_3d_raw_modes.csv"
        if parser_included:
            try:
                parsed_lambda = float(parsed_row["Lambda_3d_fem"])
            except (KeyError, TypeError, ValueError):
                parsed_lambda = float("nan")
            if math.isfinite(parsed_lambda) and abs(parsed_lambda - converted) > 1.0e-5 * max(1.0, abs(converted)):
                reason = "parser Lambda differs from raw table conversion"
        notes = [
            f"eigenvalue={float(raw['eigenvalue']):.8g}",
            f"cyclic_frequency={float(raw['cyclic_frequency']):.8g}",
            f"included_in_sorted_match={'yes' if mode in matched_modes else 'no'}",
        ]
        rows.append(
            {
                "raw_mode_number": mode,
                "raw_frequency_value": angular,
                "raw_frequency_units_if_known": "angular_frequency_rad_per_time",
                "converted_Lambda": converted,
                "parser_included": "yes" if parser_included else "no",
                "parser_mode_index": mode if parser_included else "",
                "reason_if_excluded": reason,
                "notes": "; ".join(notes),
            }
        )
    return rows


def expected_first_bending(*, epsilon: float, total_length: float) -> tuple[float, float, float]:
    Lambda = FIRST_CLAMPED_CLAMPED_ALPHA / float(total_length)
    omega = float(epsilon) * Lambda**2
    cyclic = omega / (2.0 * math.pi)
    return Lambda, omega, cyclic


def first_values(rows: Sequence[dict[str, object]], key: str, limit: int = 10) -> list[float]:
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


def low_mode_status(raw_lambdas: Sequence[float], expected_roots: Sequence[float], *, rel_tol: float = 0.05) -> str:
    missing: list[float] = []
    for expected in expected_roots[:6]:
        if not any(abs(raw - expected) / expected <= rel_tol for raw in raw_lambdas):
            missing.append(expected)
    if missing:
        return "missing low EB-like roots: " + ", ".join(f"{value:.7g}" for value in missing)
    return "low EB-like roots present within tolerance"


def mesh_warning_lines(smoke_dir: Path) -> list[str]:
    report_path = smoke_dir / "report.md"
    if not report_path.exists():
        return []
    lines = []
    for line in report_path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if "worst distortion" in line or "jac." in line:
            lines.append(line.lstrip("- ").strip())
    return lines


def write_report(
    *,
    path: Path,
    smoke_dir: Path,
    input_path: Path | None,
    dat_path: Path | None,
    frd_path: Path | None,
    frequency_card: Sequence[str],
    requested_modes: int | None,
    raw_rows: Sequence[dict[str, float | int]],
    audit_rows: Sequence[dict[str, object]],
    analytic_rows: Sequence[dict[str, object]],
    fem_rows: Sequence[dict[str, object]],
    nearest_rows: Sequence[dict[str, object]],
    epsilon: float,
    total_length: float,
    mesh_warnings: Sequence[str],
) -> None:
    expected_Lambda, expected_omega, expected_cyclic = expected_first_bending(
        epsilon=epsilon,
        total_length=total_length,
    )
    raw_lambdas = first_values(audit_rows, "converted_Lambda", limit=10)
    parsed_lambdas = first_values(fem_rows, "Lambda_3d_fem", limit=10)
    analytic_lambdas = first_values(analytic_rows, "Lambda", limit=12)
    low_status = low_mode_status(raw_lambdas, analytic_lambdas)
    parser_included_count = sum(1 for row in audit_rows if row["parser_included"] == "yes")
    nearest_errors = sorted(
        float(row["rel_error"])
        for row in nearest_rows
        if math.isfinite(float(row.get("rel_error", float("nan"))))
    )
    nearest_summary = (
        "unavailable"
        if not nearest_errors
        else f"max={max(nearest_errors):.6g}, median={nearest_errors[len(nearest_errors) // 2]:.6g}"
    )
    lines = [
        "# 3D FEM Extraction Audit",
        "",
        "Diagnostic-only audit of the straight-uniform smoke FEM extraction path.",
        "",
        "## Files",
        "",
        f"- smoke directory: `{smoke_dir.relative_to(REPO_ROOT)}`",
        f"- CalculiX input: `{input_path.relative_to(REPO_ROOT)}`" if input_path else "- CalculiX input: missing",
        f"- CalculiX dat: `{dat_path.relative_to(REPO_ROOT)}`" if dat_path else "- CalculiX dat: missing",
        f"- CalculiX frd: `{frd_path.relative_to(REPO_ROOT)}`" if frd_path else "- CalculiX frd: missing",
        "",
        "## Frequency Request",
        "",
        "```text",
        *(frequency_card or ["not found"]),
        "```",
        "",
        f"- eigenmodes requested: {requested_modes if requested_modes is not None else 'unknown'}",
        f"- lower frequency bound / shift / interval present: {has_frequency_shift_or_bounds(frequency_card)}",
        "- rigid-body modes requested or filtered: no explicit request/filter found; both end faces are fixed.",
        "",
        "## Unit Conversion",
        "",
        f"- expected first bending Lambda: {expected_Lambda:.16g}",
        f"- expected first bending angular frequency omega = epsilon*Lambda^2: {expected_omega:.16g}",
        f"- expected first bending cyclic frequency f = omega/(2*pi): {expected_cyclic:.16g}",
        "- CalculiX table header reports the third numeric column as `FREQUENCY REAL PART (RAD/TIME)`",
        "  and the fourth numeric column as cyclic frequency.",
        "- Parser conversion used here: `Lambda = sqrt(omega/epsilon)` with omega from the",
        "  RAD/TIME column. This is equivalent to `sqrt(sqrt(eigenvalue)/epsilon)`.",
        "",
        "## Raw Versus Parsed",
        "",
        f"- raw CalculiX modes listed: {len(raw_rows)}",
        f"- parser-included raw modes: {parser_included_count}",
        f"- parser skipped low modes: {'no' if parser_included_count == len(raw_rows) else 'yes'}",
        f"- raw output low-mode status: {low_status}",
        f"- first 10 raw FEM Lambdas: {', '.join(f'{value:.8g}' for value in raw_lambdas) if raw_lambdas else 'unavailable'}",
        f"- first 10 parsed FEM Lambdas: {', '.join(f'{value:.8g}' for value in parsed_lambdas) if parsed_lambdas else 'unavailable'}",
        "",
        "## Nearest-Frequency Diagnostic",
        "",
        f"- nearest-frequency relative errors: {nearest_summary}",
        "- nearest matching is diagnostic only and does not replace sorted-index output or branch identity.",
        "",
        "## Mesh Quality",
        "",
    ]
    lines.extend([f"- {warning}" for warning in mesh_warnings] or ["- no mesh warning line found in smoke report"])
    lines.append("- TODO: after extraction is fixed, rerun with a cleaner mesh or convergence check.")
    lines.extend(["", "## Conclusion", ""])
    if "low EB-like roots present" in low_status:
        lines.extend(
            [
                "The raw CalculiX frequency table contains the low EB-like roots and the parser",
                "includes the raw modes it sees. This supports a request/extraction issue in",
                "the earlier low-mode-missing run rather than a parser, unit-conversion, or",
                "analytic-reference issue.",
            ]
        )
    elif parser_included_count == len(raw_rows):
        lines.extend(
            [
                "The raw CalculiX frequency table starts at the high converted Lambda values",
                "reported by `fem_3d_raw_modes.csv`. The parser includes the raw modes it sees;",
                "it is not dropping the missing low EB-like roots. The issue is upstream of the",
                "parser, most likely in the solver request/extraction behavior for this run,",
                "with geometry/boundary-condition/solid-model scaling/mesh still to be checked",
                "if an explicit larger mode request does not recover the low modes.",
            ]
        )
    else:
        lines.extend(
            [
                "The raw CalculiX table and parsed FEM CSV disagree. This points to a parser or",
                "post-processing issue that should be fixed before interpreting FEM errors.",
            ]
        )
    lines.extend(
        [
            "",
            "## Generated Files",
            "",
            f"- `{(smoke_dir / 'fem_3d_raw_frequency_audit.csv').relative_to(REPO_ROOT)}`",
            f"- `{(smoke_dir / 'nearest_match_analytic_vs_3d_fem.csv').relative_to(REPO_ROOT)}`",
            f"- `{path.relative_to(REPO_ROOT)}`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit raw CalculiX frequency extraction for the 3D smoke case.")
    parser.add_argument("--smoke-dir", type=Path, default=DEFAULT_SMOKE_DIR)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--total-length", type=float, default=DEFAULT_TOTAL_LENGTH)
    return parser.parse_args(argv)


def run(args: argparse.Namespace) -> dict[str, object]:
    smoke_dir = _repo_path(args.smoke_dir)
    solid_dir = smoke_dir / "solid_fem_case"
    input_path = find_single(solid_dir, "*_ccx_modal.inp")
    dat_path = find_single(solid_dir, "*_ccx_modal.dat")
    frd_path = find_single(solid_dir, "*_ccx_modal.frd")
    frequency_card = frequency_card_lines(input_path) if input_path else []
    raw_rows = parse_calculix_frequency_table(dat_path) if dat_path else []

    analytic_rows = read_csv(smoke_dir / "analytic_union.csv")
    fem_rows = read_csv(smoke_dir / "fem_3d_raw_modes.csv")
    match_rows = read_csv(smoke_dir / "analytic_vs_3d_fem_matched.csv")
    audit_rows = raw_frequency_audit_rows(raw_rows, fem_rows, match_rows, epsilon=float(args.epsilon))

    raw_audit_path = smoke_dir / "fem_3d_raw_frequency_audit.csv"
    write_csv(
        raw_audit_path,
        audit_rows,
        [
            "raw_mode_number",
            "raw_frequency_value",
            "raw_frequency_units_if_known",
            "converted_Lambda",
            "parser_included",
            "parser_mode_index",
            "reason_if_excluded",
            "notes",
        ],
    )

    compare = _load_module(
        "compare_full_spectrum_analytic_vs_3d_fem_for_extraction_audit",
        REPO_ROOT / "scripts" / "analysis" / "thickness_mismatch" / "audits" / "compare_full_spectrum_analytic_vs_3d_fem.py",
    )
    nearest_rows = compare.nearest_frequency_match_rows(analytic_rows, fem_rows)
    nearest_path = smoke_dir / "nearest_match_analytic_vs_3d_fem.csv"
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

    report_path = smoke_dir / "3d_fem_extraction_audit.md"
    write_report(
        path=report_path,
        smoke_dir=smoke_dir,
        input_path=input_path,
        dat_path=dat_path,
        frd_path=frd_path,
        frequency_card=frequency_card,
        requested_modes=requested_mode_count(frequency_card),
        raw_rows=raw_rows,
        audit_rows=audit_rows,
        analytic_rows=analytic_rows,
        fem_rows=fem_rows,
        nearest_rows=nearest_rows,
        epsilon=float(args.epsilon),
        total_length=float(args.total_length),
        mesh_warnings=mesh_warning_lines(smoke_dir),
    )
    print(f"saved raw frequency audit CSV: {raw_audit_path}")
    print(f"saved nearest-frequency CSV: {nearest_path}")
    print(f"saved extraction audit report: {report_path}")
    print(f"raw modes listed: {len(raw_rows)}")
    print(f"parser-included modes: {sum(1 for row in audit_rows if row['parser_included'] == 'yes')}")
    return {"raw_audit_path": raw_audit_path, "nearest_path": nearest_path, "report_path": report_path}


def main(argv: list[str] | None = None) -> int:
    run(parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
