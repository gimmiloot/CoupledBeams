from __future__ import annotations

import argparse
import csv
import importlib.util
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import ModuleType
from typing import Sequence


REPO_ROOT = Path(__file__).resolve().parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))


DEFAULT_OUTPUT_DIR = Path("results") / "_smoke" / "3d_fem_environment_check"
DEFAULT_EPSILON = 0.02
DEFAULT_POISSON = 0.3
DEFAULT_TOTAL_LENGTH = 2.0
DEFAULT_MESH_SIZE = 0.04
DEFAULT_REQUESTED_MODES = 40
DEFAULT_TIMEOUT_SECONDS = 300
DEFAULT_GMSH_CANDIDATES = (
    Path(r"D:\PHD\gmsh-4.15.2-Windows64\gmsh.exe"),
    Path(r"D:\PHD\gmsh\gmsh.exe"),
)
DEFAULT_CCX_CANDIDATES = (
    Path(r"D:\PHD\calculix\ccx.exe"),
    Path(r"D:\PHD\calculix\ccx_static.exe"),
    Path(r"D:\PHD\calculix\calculix_2.22_4win\ccx_static.exe"),
)
CALCULIX_RELEASE_ZIPS = (
    Path(r"D:\PHD\CalculiX-Windows-master\releases\CalculiX-2.23.0-win-x64.zip"),
    Path(r"D:\PHD\CalculiX-Windows-master\releases\CalculiX-2.22.0-win-x64.zip"),
    Path(r"D:\PHD\CalculiX-Windows-master\releases\CalculiX-2.21.0-win-x64.zip"),
    Path(r"D:\PHD\CalculiX-Windows-master\releases\CalculiX-2.20.0-win-x64.zip"),
)
WORKFLOW = "straight_fixed_fixed_cylinder_from_single_rod_helpers"


EXPECTED_ALPHA_ROOTS = (
    4.730040744862704,
    7.853204624095838,
    10.99560783800167,
    14.137165491257,
    17.278759532088237,
    20.42035225104125,
)


class ExecutableResolution:
    def __init__(
        self,
        *,
        path: Path | None,
        source: str,
        candidates_checked: Sequence[Path],
        archive_candidates: Sequence[Path] = (),
    ) -> None:
        self.path = path
        self.source = source
        self.candidates_checked = tuple(candidates_checked)
        self.archive_candidates = tuple(archive_candidates)


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


def resolve_existing(path_text: str | None) -> Path | None:
    if not path_text:
        return None
    path = Path(path_text.strip().strip('"'))
    if path.is_file():
        return path
    return None


def find_under(root: Path, names: Sequence[str], *, limit: int = 16) -> list[Path]:
    if not root.exists():
        return []
    wanted = {name.lower() for name in names}
    out: list[Path] = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [name for name in dirnames if name not in {"$RECYCLE.BIN", "System Volume Information"}]
        for filename in filenames:
            if filename.lower() in wanted:
                out.append(Path(dirpath) / filename)
                if len(out) >= limit:
                    return out
    return out


def release_zips_with_ccx() -> tuple[Path, ...]:
    return tuple(path for path in CALCULIX_RELEASE_ZIPS if path.exists())


def resolve_executable(
    *,
    explicit: str | None,
    env_name: str,
    path_names: Sequence[str],
    default_candidates: Sequence[Path],
    search_root: Path | None,
    archive_candidates: Sequence[Path] = (),
) -> ExecutableResolution:
    checked: list[Path] = []
    explicit_path = resolve_existing(explicit)
    if explicit_path is not None:
        return ExecutableResolution(path=explicit_path, source="explicit CLI path", candidates_checked=[explicit_path])
    if explicit:
        checked.append(Path(explicit.strip().strip('"')))

    env_path = resolve_existing(os.environ.get(env_name))
    if env_path is not None:
        return ExecutableResolution(path=env_path, source=f"environment variable {env_name}", candidates_checked=[env_path])
    if os.environ.get(env_name):
        checked.append(Path(os.environ[env_name].strip().strip('"')))

    for name in path_names:
        found = shutil.which(name)
        if found:
            return ExecutableResolution(path=Path(found), source=f"PATH lookup for {name}", candidates_checked=checked)

    for candidate in default_candidates:
        checked.append(candidate)
        if candidate.is_file():
            return ExecutableResolution(path=candidate, source="known D:\\PHD candidate", candidates_checked=checked)

    if search_root is not None:
        found_under_root = find_under(search_root, path_names)
        checked.extend(found_under_root)
        if found_under_root:
            return ExecutableResolution(path=found_under_root[0], source=f"recursive search under {search_root}", candidates_checked=checked)

    return ExecutableResolution(
        path=None,
        source="not found",
        candidates_checked=checked,
        archive_candidates=archive_candidates,
    )


def run_probe(command: Sequence[str], *, cwd: Path, timeout: int) -> dict[str, object]:
    try:
        completed = subprocess.run(
            list(command),
            cwd=str(cwd),
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=timeout,
        )
        return {
            "command": " ".join(str(part) for part in command),
            "returncode": completed.returncode,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "error": "",
        }
    except Exception as exc:
        return {
            "command": " ".join(str(part) for part in command),
            "returncode": "",
            "stdout": "",
            "stderr": "",
            "error": str(exc),
        }


def first_probe_line(probe: dict[str, object]) -> str:
    text = "\n".join(str(probe.get(key, "")) for key in ("stdout", "stderr", "error"))
    for line in text.splitlines():
        stripped = line.strip()
        if stripped:
            return stripped
    return "no output"


def collect_warning_lines(text: str) -> tuple[str, ...]:
    out: list[str] = []
    for line in text.splitlines():
        lower = line.lower()
        if "warning" in lower or "error" in lower or "jac." in lower or "distortion" in lower:
            out.append(line.strip())
    return tuple(out)


def run_gmsh(
    *,
    gmsh_exe: Path,
    paths: object,
    output_dir: Path,
    timeout: int,
) -> dict[str, object]:
    stdout_path = output_dir / "gmsh_stdout.txt"
    stderr_path = output_dir / "gmsh_stderr.txt"
    commands = (
        (str(gmsh_exe), str(paths.geo), "-3", "-format", "msh4", "-o", str(paths.msh)),
        (str(gmsh_exe), str(paths.geo), "-3", "-format", "inp", "-o", str(paths.gmsh_inp)),
    )
    stdout_parts: list[str] = []
    stderr_parts: list[str] = []
    returncodes: list[int | str] = []
    for command in commands:
        result = run_probe(command, cwd=paths.case_dir, timeout=timeout)
        stdout_parts.extend([f"$ {' '.join(command)}", str(result["stdout"])])
        stderr_parts.extend([f"$ {' '.join(command)}", str(result["stderr"]), str(result["error"])])
        returncodes.append(result["returncode"])
        if result["returncode"] != 0:
            break
    stdout_text = "\n".join(stdout_parts)
    stderr_text = "\n".join(stderr_parts)
    stdout_path.write_text(stdout_text, encoding="utf-8")
    stderr_path.write_text(stderr_text, encoding="utf-8")
    ok = bool(returncodes) and all(code == 0 for code in returncodes)
    return {
        "ok": ok,
        "returncodes": tuple(returncodes),
        "stdout_path": stdout_path,
        "stderr_path": stderr_path,
        "warnings": collect_warning_lines(stdout_text + "\n" + stderr_text),
    }


def run_ccx(
    *,
    ccx_exe: Path,
    paths: object,
    output_dir: Path,
    timeout: int,
) -> dict[str, object]:
    stdout_path = output_dir / "ccx_stdout.txt"
    stderr_path = output_dir / "ccx_stderr.txt"
    command = (str(ccx_exe), paths.ccx_modal_inp.stem)
    result = run_probe(command, cwd=paths.case_dir, timeout=timeout)
    stdout_text = str(result["stdout"])
    stderr_text = "\n".join(part for part in (str(result["stderr"]), str(result["error"])) if part)
    stdout_path.write_text(stdout_text, encoding="utf-8")
    stderr_path.write_text(stderr_text, encoding="utf-8")
    paths.ccx_stdout.write_text(stdout_text, encoding="utf-8")
    paths.ccx_stderr.write_text(stderr_text, encoding="utf-8")
    return {
        "ok": result["returncode"] == 0,
        "returncode": result["returncode"],
        "stdout_path": stdout_path,
        "stderr_path": stderr_path,
        "warnings": collect_warning_lines(stdout_text + "\n" + stderr_text),
    }


def parse_raw_calculix_modes(
    *,
    dat_path: Path,
    stdout_path: Path,
    epsilon: float,
    fallback_parser: object,
) -> tuple[list[dict[str, object]], str]:
    extraction = _load_module(
        "audit_full_spectrum_3d_fem_smoke_extraction_for_environment_check",
        REPO_ROOT / "scripts" / "analysis" / "thickness_mismatch" / "audits" / "audit_full_spectrum_3d_fem_smoke_extraction.py",
    )
    table_rows = extraction.parse_calculix_frequency_table(dat_path)
    if table_rows:
        rows = []
        for row in table_rows:
            omega = float(row["angular_frequency"])
            rows.append(
                {
                    "raw_mode_number": int(row["raw_mode_number"]),
                    "omega_rad_per_time": omega,
                    "Lambda": math.sqrt(omega / float(epsilon)),
                    "parser_included": "yes",
                    "notes": "parsed from CalculiX .dat FREQUENCY REAL PART (RAD/TIME) column",
                }
            )
        return rows, "CalculiX .dat frequency table RAD/TIME column"

    omegas, parse_source = fallback_parser(dat_path, stdout_path)
    rows = []
    for index, omega in enumerate(omegas, start=1):
        rows.append(
            {
                "raw_mode_number": index,
                "omega_rad_per_time": float(omega),
                "Lambda": math.sqrt(float(omega) / float(epsilon)),
                "parser_included": "fallback",
                "notes": parse_source,
            }
        )
    return rows, f"fallback helper parser: {parse_source}"


def expected_doublet_lambdas(limit: int) -> list[float]:
    values: list[float] = []
    for alpha in EXPECTED_ALPHA_ROOTS:
        values.extend([float(alpha) / 2.0, float(alpha) / 2.0])
    return values[:limit]


def comparison_rows(raw_rows: Sequence[dict[str, object]], *, expected_count: int) -> list[dict[str, object]]:
    available = [
        (int(row["raw_mode_number"]), float(row["Lambda"]))
        for row in raw_rows
        if math.isfinite(float(row.get("Lambda", float("nan"))))
    ]
    rows: list[dict[str, object]] = []
    for index, expected in enumerate(expected_doublet_lambdas(expected_count), start=1):
        if not available:
            rows.append(
                {
                    "expected_mode_index": index,
                    "Lambda_expected_EB": expected,
                    "nearest_fem_mode_index": "",
                    "Lambda_3D_FEM": float("nan"),
                    "abs_error": float("nan"),
                    "rel_error": float("nan"),
                    "notes": "no parsed FEM mode available",
                }
            )
            continue
        mode, Lambda = min(available, key=lambda item: abs(item[1] - expected))
        available.remove((mode, Lambda))
        abs_error = abs(Lambda - expected)
        rows.append(
            {
                "expected_mode_index": index,
                "Lambda_expected_EB": expected,
                "nearest_fem_mode_index": mode,
                "Lambda_3D_FEM": Lambda,
                "abs_error": abs_error,
                "rel_error": abs_error / abs(expected) if expected else float("nan"),
                "notes": "greedy nearest-frequency smoke comparison; not branch identity",
            }
        )
    return rows


def write_report(
    *,
    path: Path,
    args: argparse.Namespace,
    gmsh_resolution: ExecutableResolution,
    ccx_resolution: ExecutableResolution,
    gmsh_probe: dict[str, object],
    ccx_probe: dict[str, object],
    gmsh_run: dict[str, object] | None,
    ccx_run: dict[str, object] | None,
    mesh_info: dict[str, object],
    ccx_input_info: dict[str, object],
    raw_rows: Sequence[dict[str, object]],
    compare_rows: Sequence[dict[str, object]],
    parse_source: str,
    generated_files: Sequence[Path],
) -> None:
    first10 = [float(row["Lambda"]) for row in raw_rows[:10]]
    low_bands = {
        "2.3-2.5": any(2.3 <= float(row["Lambda"]) <= 2.5 for row in raw_rows),
        "3.8-4.1": any(3.8 <= float(row["Lambda"]) <= 4.1 for row in raw_rows),
        "5.3-5.8": any(5.3 <= float(row["Lambda"]) <= 5.8 for row in raw_rows),
    }
    ready = (
        gmsh_resolution.path is not None
        and ccx_resolution.path is not None
        and gmsh_run is not None
        and bool(gmsh_run["ok"])
        and ccx_run is not None
        and bool(ccx_run["ok"])
        and len(raw_rows) >= 20
        and all(low_bands.values())
    )
    blocker = "none"
    if gmsh_resolution.path is None:
        blocker = "Gmsh executable was not found."
    elif ccx_resolution.path is None:
        blocker = "CalculiX ccx executable was not found."
    elif gmsh_run is None or not bool(gmsh_run["ok"]):
        blocker = "Gmsh mesh generation did not complete successfully."
    elif ccx_run is None or not bool(ccx_run["ok"]):
        blocker = "CalculiX modal run did not complete successfully."
    elif len(raw_rows) < 20:
        blocker = f"Parsed only {len(raw_rows)} modes; expected at least 20 for this smoke check."
    elif not all(low_bands.values()):
        blocker = "Parsed modes did not include all low bending Lambda bands 2.3-2.5, 3.8-4.1, and 5.3-5.8."

    lines = [
        "# 3D FEM Environment Check",
        "",
        "Diagnostic-only Gmsh + CalculiX smoke check for the straight fixed-fixed cylinder.",
        "",
        "## Parameters",
        "",
        "- beta: 0 deg",
        "- mu: 0",
        "- eta: 0",
        f"- epsilon: {args.epsilon:g}",
        f"- poisson: {args.poisson:g}",
        "- E: 1",
        "- rho: 1",
        f"- total length L: {args.total_length:g}",
        f"- radius r0 = 2*epsilon: {2.0 * args.epsilon:.16g}",
        f"- mesh size h: {args.mesh_size:.16g}",
        f"- requested modes: {args.n_fem_modes}",
        f"- workflow: `{WORKFLOW}`",
        "",
        "## Executables",
        "",
        f"- gmsh path: `{gmsh_resolution.path}`" if gmsh_resolution.path else "- gmsh path: not found",
        f"- gmsh resolution source: {gmsh_resolution.source}",
        f"- gmsh version/probe command: `{gmsh_probe['command']}`",
        f"- gmsh version/probe return code: {gmsh_probe['returncode']}",
        f"- gmsh version/probe first line: {first_probe_line(gmsh_probe)}",
        f"- ccx path: `{ccx_resolution.path}`" if ccx_resolution.path else "- ccx path: not found",
        f"- ccx resolution source: {ccx_resolution.source}",
        f"- ccx version/probe command: `{ccx_probe['command']}`" if ccx_probe else "- ccx version/probe command: not run",
        f"- ccx version/probe return code: {ccx_probe['returncode']}" if ccx_probe else "- ccx version/probe return code: not run",
        f"- ccx version/probe first line: {first_probe_line(ccx_probe)}" if ccx_probe else "- ccx version/probe first line: not run",
        "",
    ]
    if ccx_resolution.archive_candidates and ccx_resolution.path is None:
        lines.extend(["CalculiX release archives found but not unpacked as `ccx.exe`:", ""])
        lines.extend(f"- `{item}`" for item in ccx_resolution.archive_candidates)
        lines.extend(
            [
                "",
                "Unpack one release and either add its `bin` folder to PATH or pass",
                "`--ccx-exe path\\to\\ccx.exe`.",
                "",
            ]
        )

    lines.extend(
        [
            "## Run Status",
            "",
            f"- Gmsh completed: {bool(gmsh_run and gmsh_run['ok'])}",
            f"- Gmsh return codes: {gmsh_run['returncodes'] if gmsh_run else 'not run'}",
            f"- CalculiX completed: {bool(ccx_run and ccx_run['ok'])}",
            f"- CalculiX return code: {ccx_run['returncode'] if ccx_run else 'not run'}",
            f"- parser source: {parse_source}",
            f"- environment ready for straight-cylinder 3D FEM smoke work: {ready}",
            f"- blocker: {blocker}",
            "",
            "## Mesh",
            "",
            f"- nodes: {mesh_info.get('nodes', 0)}",
            f"- solid elements: {mesh_info.get('elements', 0)}",
            f"- inferred length: {mesh_info.get('inferred_length', 'unavailable')}",
            f"- inferred radius: {mesh_info.get('inferred_radius', 'unavailable')}",
            f"- LEFT_FIXED nodes: {ccx_input_info.get('left_fixed_nodes', 0)}",
            f"- RIGHT_FIXED nodes: {ccx_input_info.get('right_fixed_nodes', 0)}",
            "",
            "Mesh/Gmsh warnings:",
            "",
        ]
    )
    gmsh_warnings = tuple(gmsh_run.get("warnings", ())) if gmsh_run else ()
    lines.extend([f"- {warning}" for warning in gmsh_warnings] or ["- none"])
    lines.extend(
        [
            "",
            "CalculiX warnings/errors:",
            "",
        ]
    )
    ccx_warnings = tuple(ccx_run.get("warnings", ())) if ccx_run else ()
    lines.extend([f"- {warning}" for warning in ccx_warnings] or ["- none"])
    lines.extend(
        [
            "",
            "## Parsed Frequencies",
            "",
            f"- parsed modes: {len(raw_rows)}",
            f"- first 10 Lambdas: {', '.join(f'{value:.8g}' for value in first10) if first10 else 'unavailable'}",
            f"- low band 2.3-2.5 present: {low_bands['2.3-2.5']}",
            f"- low band 3.8-4.1 present: {low_bands['3.8-4.1']}",
            f"- low band 5.3-5.8 present: {low_bands['5.3-5.8']}",
            "",
            "## Expected Bending Doublet Comparison",
            "",
            "| expected mode | EB Lambda | FEM mode | FEM Lambda | rel error |",
            "| ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in compare_rows:
        fem_mode = row["nearest_fem_mode_index"] if row["nearest_fem_mode_index"] != "" else "-"
        fem_lambda = float(row["Lambda_3D_FEM"])
        fem_text = f"{fem_lambda:.8g}" if math.isfinite(fem_lambda) else "-"
        rel = float(row["rel_error"])
        rel_text = f"{rel:.6g}" if math.isfinite(rel) else "-"
        lines.append(
            f"| {int(row['expected_mode_index'])} | {float(row['Lambda_expected_EB']):.8g} | "
            f"{fem_mode} | {fem_text} | {rel_text} |"
        )
    lines.extend(
        [
            "",
            "## Generated Files",
            "",
        ]
    )
    for item in [*generated_files, path]:
        try:
            lines.append(f"- `{item.relative_to(REPO_ROOT)}`")
        except ValueError:
            lines.append(f"- `{item}`")
    lines.extend(
        [
            "",
            "## Boundaries",
            "",
            "- No beta scan was run.",
            "- No thickness, eta, or mu scan was run.",
            "- No coupled-angle 3D FEM case was run.",
            "- No analytic formulas, determinants, old solvers, article workspace, article figures, or baseline results were modified.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Diagnostic-only Gmsh + CalculiX 3D FEM environment smoke check.")
    parser.add_argument("--gmsh-exe", default=None)
    parser.add_argument("--ccx-exe", default=None)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--poisson", type=float, default=DEFAULT_POISSON)
    parser.add_argument("--total-length", type=float, default=DEFAULT_TOTAL_LENGTH)
    parser.add_argument("--mesh-size", type=float, default=DEFAULT_MESH_SIZE)
    parser.add_argument("--n-fem-modes", type=int, default=DEFAULT_REQUESTED_MODES)
    parser.add_argument("--timeout-seconds", type=int, default=DEFAULT_TIMEOUT_SECONDS)
    parser.add_argument("--compare-expected-modes", type=int, default=6)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--also-check-eps0p0025", action="store_true")
    args = parser.parse_args(argv)
    if args.epsilon <= 0.0 or args.mesh_size <= 0.0 or args.total_length <= 0.0:
        raise ValueError("epsilon, mesh size, and total length must be positive.")
    if args.n_fem_modes < 40:
        raise ValueError("this environment smoke check requires at least 40 requested FEM modes.")
    args.output_dir = _repo_path(Path(args.output_dir))
    return args


def run_primary_case(args: argparse.Namespace) -> dict[str, object]:
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    single = _load_module(
        "solid_fem_single_rod_fixed_fixed_for_environment_check",
        REPO_ROOT / "scripts" / "analysis" / "solid_fem_single_rod_fixed_fixed.py",
    )
    single.L = float(args.total_length)
    single.NU = float(args.poisson)
    single.SOLID_MODES_REQUESTED = int(args.n_fem_modes)
    single.CCX_TIMEOUT_SECONDS = int(args.timeout_seconds)
    single.mesh_size_for_epsilon = lambda _epsilon: float(args.mesh_size)

    gmsh_resolution = resolve_executable(
        explicit=args.gmsh_exe,
        env_name="GMSH_EXE",
        path_names=("gmsh.exe", "gmsh"),
        default_candidates=DEFAULT_GMSH_CANDIDATES,
        search_root=Path(r"D:\PHD"),
    )
    ccx_resolution = resolve_executable(
        explicit=args.ccx_exe,
        env_name="CCX_EXE",
        path_names=("ccx.exe", "ccx_MT.exe", "ccx_static.exe", "ccx_dynamic.exe", "ccx"),
        default_candidates=DEFAULT_CCX_CANDIDATES,
        search_root=Path(r"D:\PHD"),
        archive_candidates=release_zips_with_ccx(),
    )
    gmsh_probe = (
        run_probe((str(gmsh_resolution.path), "--version"), cwd=output_dir, timeout=30)
        if gmsh_resolution.path
        else {"command": "gmsh --version", "returncode": "not run", "stdout": "", "stderr": "", "error": "gmsh not found"}
    )
    ccx_probe = (
        run_probe((str(ccx_resolution.path), "-v"), cwd=output_dir, timeout=30)
        if ccx_resolution.path
        else {"command": "ccx -v", "returncode": "not run", "stdout": "", "stderr": "", "error": "ccx not found"}
    )

    case_dir = output_dir / "solid_fem_case"
    stem = f"straight_uniform_L2_eps{_token(float(args.epsilon))}"
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
    generated_files: list[Path] = [paths.geo]

    gmsh_run = None
    ccx_run = None
    mesh_info: dict[str, object] = {"nodes": 0, "elements": 0, "inferred_length": "unavailable", "inferred_radius": "unavailable"}
    ccx_input_info: dict[str, object] = {"left_fixed_nodes": 0, "right_fixed_nodes": 0}
    raw_rows: list[dict[str, object]] = []
    compare_rows: list[dict[str, object]] = comparison_rows(raw_rows, expected_count=int(args.compare_expected_modes))
    parse_source = "not attempted"

    if gmsh_resolution.path is not None:
        gmsh_run = run_gmsh(gmsh_exe=gmsh_resolution.path, paths=paths, output_dir=output_dir, timeout=int(args.timeout_seconds))
        generated_files.extend([output_dir / "gmsh_stdout.txt", output_dir / "gmsh_stderr.txt"])
        if gmsh_run["ok"] and paths.gmsh_inp.exists():
            mesh = single.read_gmsh_inp_mesh_data(paths.gmsh_inp)
            length = mesh.bbox_max[0] - mesh.bbox_min[0] if mesh.bbox_min and mesh.bbox_max else float("nan")
            mesh_info = {
                "nodes": len(mesh.nodes),
                "elements": len(mesh.solid_elements),
                "inferred_length": length,
                "inferred_radius": mesh.radial_max,
            }
            ccx_input = single.write_calculix_inputs(paths, float(args.epsilon))
            ccx_input_info = {
                "left_fixed_nodes": ccx_input.left_fixed_node_count,
                "right_fixed_nodes": ccx_input.right_fixed_node_count,
            }
            generated_files.extend([paths.gmsh_inp, paths.ccx_mesh_inp, paths.ccx_modal_inp])

            if ccx_resolution.path is not None:
                ccx_run = run_ccx(ccx_exe=ccx_resolution.path, paths=paths, output_dir=output_dir, timeout=int(args.timeout_seconds))
                generated_files.extend([output_dir / "ccx_stdout.txt", output_dir / "ccx_stderr.txt"])
                if ccx_run["ok"]:
                    raw_rows, parse_source = parse_raw_calculix_modes(
                        dat_path=paths.ccx_dat,
                        stdout_path=paths.ccx_stdout,
                        epsilon=float(args.epsilon),
                        fallback_parser=single.parse_calculix_eigen_omegas,
                    )
                    compare_rows = comparison_rows(raw_rows, expected_count=int(args.compare_expected_modes))
                    generated_files.extend([paths.ccx_dat, paths.ccx_frd])

    placeholder_logs = {
        output_dir / "gmsh_stdout.txt": "Gmsh was not run because gmsh.exe was not resolved.\n",
        output_dir / "gmsh_stderr.txt": "Gmsh was not run because gmsh.exe was not resolved.\n",
        output_dir / "ccx_stdout.txt": "CalculiX was not run because ccx.exe was not resolved or mesh generation failed.\n",
        output_dir / "ccx_stderr.txt": "CalculiX was not run because ccx.exe was not resolved or mesh generation failed.\n",
    }
    for log_path, text in placeholder_logs.items():
        if not log_path.exists():
            log_path.write_text(text, encoding="utf-8")
            generated_files.append(log_path)

    raw_path = output_dir / "fem_3d_raw_modes.csv"
    comparison_path = output_dir / "fem_3d_environment_check_comparison.csv"
    report_path = output_dir / "3d_fem_environment_check_report.md"
    write_csv(
        raw_path,
        raw_rows,
        ["raw_mode_number", "omega_rad_per_time", "Lambda", "parser_included", "notes"],
    )
    write_csv(
        comparison_path,
        compare_rows,
        ["expected_mode_index", "Lambda_expected_EB", "nearest_fem_mode_index", "Lambda_3D_FEM", "abs_error", "rel_error", "notes"],
    )
    generated_files.extend([raw_path, comparison_path])
    write_report(
        path=report_path,
        args=args,
        gmsh_resolution=gmsh_resolution,
        ccx_resolution=ccx_resolution,
        gmsh_probe=gmsh_probe,
        ccx_probe=ccx_probe,
        gmsh_run=gmsh_run,
        ccx_run=ccx_run,
        mesh_info=mesh_info,
        ccx_input_info=ccx_input_info,
        raw_rows=raw_rows,
        compare_rows=compare_rows,
        parse_source=parse_source,
        generated_files=generated_files,
    )

    print(f"gmsh: {gmsh_resolution.path if gmsh_resolution.path else 'not found'}")
    print(f"ccx: {ccx_resolution.path if ccx_resolution.path else 'not found'}")
    print(f"wrote report: {report_path}")
    print(f"wrote raw modes CSV: {raw_path}")
    print(f"wrote comparison CSV: {comparison_path}")
    print(f"parsed modes: {len(raw_rows)}")
    return {
        "report_path": report_path,
        "raw_path": raw_path,
        "comparison_path": comparison_path,
        "gmsh_path": gmsh_resolution.path,
        "ccx_path": ccx_resolution.path,
        "parsed_modes": len(raw_rows),
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    run_primary_case(args)
    if args.also_check_eps0p0025:
        extra = parse_args(
            [
                "--epsilon",
                "0.0025",
                "--mesh-size",
                "0.004",
                "--n-fem-modes",
                str(args.n_fem_modes),
                "--timeout-seconds",
                str(args.timeout_seconds),
                "--output-dir",
                str(Path(args.output_dir) / "eps0p0025"),
                *(["--gmsh-exe", args.gmsh_exe] if args.gmsh_exe else []),
                *(["--ccx-exe", args.ccx_exe] if args.ccx_exe else []),
            ]
        )
        run_primary_case(extra)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
