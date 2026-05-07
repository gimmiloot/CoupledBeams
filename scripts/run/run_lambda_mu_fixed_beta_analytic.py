from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Sequence

import numpy as np


# ============================================================
# USER PARAMETERS
# Change these values for ordinary runs.
# CLI arguments, if provided, override these defaults.
# ============================================================

DEFAULT_BETA = 15.0
DEFAULT_EPSILON = 0.0025
DEFAULT_NUM_MODES = 7
DEFAULT_NUM_DASHED_LINES = None  # None means use DEFAULT_NUM_MODES.
DEFAULT_MU_MIN = 0.0
DEFAULT_MU_MAX = 0.9
DEFAULT_MU_STEP = 0.01
DEFAULT_Y_MAX = 20.0
DEFAULT_OUTPUT = None
DEFAULT_CSV_OUTPUT = None
DEFAULT_SHOW = False
DEFAULT_ALLOW_LOW_MAC = False


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.FreqMuNet import plot_combined  # noqa: E402
from scripts.lib.analytic_branch_tracking import (  # noqa: E402
    DEFAULT_MU_STEPS,
    DEFAULT_N_SOLVE,
    DEFAULT_N_TRACK,
    branch_id_from_base_sorted_index,
    dense_mu_values_for_targets,
    track_mu_sweep,
)


def filename_number_token(value: float) -> str:
    return f"{float(value):g}".replace("-", "m").replace("+", "").replace(".", "p")


def default_output_path(*, beta: float, epsilon: float, num_modes: int, num_dashed_lines: int) -> Path:
    return REPO_ROOT / "results" / (
        f"lambda_mu_beta{filename_number_token(beta)}_eps{filename_number_token(epsilon)}"
        f"_modes{int(num_modes)}_dash{int(num_dashed_lines)}.png"
    )


def default_csv_path(*, beta: float, epsilon: float, num_modes: int, num_dashed_lines: int) -> Path:
    return REPO_ROOT / "results" / (
        f"lambda_mu_beta{filename_number_token(beta)}_eps{filename_number_token(epsilon)}"
        f"_modes{int(num_modes)}_dash{int(num_dashed_lines)}.csv"
    )


def resolve_repo_path(value: str | Path | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Plot analytic Lambda(mu) at fixed beta for an arbitrary number of canonical "
            "tracked coupled-rods branches, using the project FreqMuNet visual style."
        ),
    )
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA, help="Fixed coupling angle in degrees.")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON, help="Slenderness parameter epsilon.")
    parser.add_argument("--num-modes", type=int, default=DEFAULT_NUM_MODES, help="Number of coupled analytic modes to draw.")
    parser.add_argument(
        "--num-dashed-lines",
        type=int,
        default=DEFAULT_NUM_DASHED_LINES,
        help="Number of dashed CS reference lines per family. Default: same as --num-modes.",
    )
    parser.add_argument("--mu-min", type=float, default=DEFAULT_MU_MIN)
    parser.add_argument("--mu-max", type=float, default=DEFAULT_MU_MAX)
    parser.add_argument("--mu-step", type=float, default=DEFAULT_MU_STEP)
    parser.add_argument("--y-max", type=float, default=DEFAULT_Y_MAX, help="Use <= 0 for automatic y limit.")
    parser.add_argument("--output", default=DEFAULT_OUTPUT, help="PNG output path. Default is deterministic under results/.")
    parser.add_argument("--csv-output", default=DEFAULT_CSV_OUTPUT, help="CSV output path. Default matches the PNG stem.")
    parser.add_argument("--show", action="store_true", default=DEFAULT_SHOW, help="Also display the plot window.")
    parser.add_argument(
        "--allow-low-mac",
        action="store_true",
        default=DEFAULT_ALLOW_LOW_MAC,
        help="Allow exploratory output even if canonical tracking falls below the MAC warning threshold.",
    )
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    if args.epsilon <= 0.0:
        parser.error("--epsilon must be positive.")
    if args.num_modes <= 0:
        parser.error("--num-modes must be positive.")
    if args.num_dashed_lines is not None and args.num_dashed_lines < 0:
        parser.error("--num-dashed-lines must be non-negative.")
    if args.mu_step <= 0.0:
        parser.error("--mu-step must be positive.")
    if args.mu_max < args.mu_min:
        parser.error("--mu-max must be greater than or equal to --mu-min.")
    if args.mu_min < -1e-12:
        parser.error("--mu-min must be non-negative for canonical beta-then-mu analytic tracking.")
    return args


def mu_grid(*, mu_min: float, mu_max: float, mu_step: float) -> np.ndarray:
    values = np.arange(float(mu_min), float(mu_max) + 0.5 * float(mu_step), float(mu_step), dtype=float)
    if values.size == 0 or abs(float(values[-1]) - float(mu_max)) > 1e-10:
        values = np.append(values, float(mu_max))
    values[0] = float(mu_min)
    values[-1] = float(mu_max)
    return np.unique(np.round(values, 12))


def collect_csv_rows(
    *,
    beta: float,
    epsilon: float,
    num_modes: int,
    mu_values: np.ndarray,
    allow_low_mac: bool = False,
) -> list[dict[str, float | int | str]]:
    n_track = max(DEFAULT_N_TRACK, int(num_modes))
    tracking_mu_values = dense_mu_values_for_targets(
        mu_values,
        mu_steps=max(DEFAULT_MU_STEPS, int(len(mu_values))),
    )
    branch_ids = [branch_id_from_base_sorted_index(index) for index in range(1, int(num_modes) + 1)]
    result = track_mu_sweep(
        epsilon=float(epsilon),
        beta=float(beta),
        mu_values=tracking_mu_values,
        n_track=n_track,
        n_solve=max(DEFAULT_N_SOLVE, n_track),
        shape_metric="full",
        allow_low_mac=allow_low_mac,
        required_branch_ids=branch_ids,
    )

    rows: list[dict[str, float | int | str]] = []
    for mode_rank in range(1, int(num_modes) + 1):
        branch_id = branch_ids[mode_rank - 1]
        for mu in mu_values:
            point = result.point_at(branch_id, beta=float(beta), mu=float(mu))
            rows.append(
                {
                    "beta": float(beta),
                    "epsilon": float(epsilon),
                    "mu": float(mu),
                    "mode_rank": int(mode_rank),
                    "branch_id": branch_id,
                    "base_sorted_index": int(mode_rank),
                    "current_sorted_index": int(point.current_sorted_index),
                    "Lambda": float(point.Lambda),
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["beta", "epsilon", "mu", "mode_rank", "branch_id", "base_sorted_index", "current_sorted_index", "Lambda"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> dict[str, Path | int | float]:
    args = parse_args(argv)
    num_dashed_lines = int(args.num_modes if args.num_dashed_lines is None else args.num_dashed_lines)
    output_path = resolve_repo_path(args.output) or default_output_path(
        beta=args.beta,
        epsilon=args.epsilon,
        num_modes=args.num_modes,
        num_dashed_lines=num_dashed_lines,
    )
    csv_path = resolve_repo_path(args.csv_output) or default_csv_path(
        beta=args.beta,
        epsilon=args.epsilon,
        num_modes=args.num_modes,
        num_dashed_lines=num_dashed_lines,
    )

    mu_values = mu_grid(mu_min=args.mu_min, mu_max=args.mu_max, mu_step=args.mu_step)
    plot_combined(
        beta_deg_coupled=float(args.beta),
        epsilon=float(args.epsilon),
        mu_min=float(args.mu_min),
        mu_max=float(args.mu_max),
        mu_step=float(args.mu_step),
        n_coupled=int(args.num_modes),
        n_dashed_lines=num_dashed_lines,
        y_max=None if float(args.y_max) <= 0.0 else float(args.y_max),
        save_path=output_path,
        show=bool(args.show),
        allow_low_mac=bool(args.allow_low_mac),
    )
    rows = collect_csv_rows(
        beta=float(args.beta),
        epsilon=float(args.epsilon),
        num_modes=int(args.num_modes),
        mu_values=mu_values,
        allow_low_mac=bool(args.allow_low_mac),
    )
    write_csv(csv_path, rows)

    print(f"saved PNG: {output_path}")
    print(f"saved CSV: {csv_path}")
    print(f"beta: {float(args.beta):g} deg")
    print(f"epsilon: {float(args.epsilon):g}")
    print(f"num_modes: {int(args.num_modes)}")
    print(f"num_dashed_lines_per_family: {num_dashed_lines}")
    print(f"mu points: {len(mu_values)}")
    print(f"allow_low_mac: {bool(args.allow_low_mac)}")
    print("tracking debug CSV: not saved")
    return {
        "output_png": output_path,
        "output_csv": csv_path,
        "beta": float(args.beta),
        "epsilon": float(args.epsilon),
        "num_modes": int(args.num_modes),
        "num_dashed_lines": int(num_dashed_lines),
    }


if __name__ == "__main__":
    main()
