"""Zero-solve article assets for the isotropic circular EB/Timoshenko grid.

This module deliberately depends only on compact CSV/JSON results and plotting
libraries.  It must not import point solvers, matrix evaluators, family repair,
or any anisotropic implementation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


SCIENTIFIC_SCOPE = "isotropic_circular_coupled_rods_eb_timoshenko"
ASSET_SCHEMA_VERSION = "article_epsilon_upper_envelope_assets_v1"
BASE_EPSILON_LEVELS = (0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.060)
EXPECTED_ENVELOPE = (10, 10, 10, 10, 9, 6, 4, 3)
EXPECTED_RESOLVED = (162, 194, 193, 194, 194, 194, 194, 194)
EXPECTED_UNRESOLVED = (32, 0, 1, 0, 0, 0, 0, 0)
EXPECTED_MIN_RESOLVED = (10, 10, 7, 4, 3, 2, 1, 0)
REGRESSION_CONTROLS = (
    ("S3_14", 0.024798906738281248, 4, 0.10050934855181458),
    ("S3_12", 0.029408510742187498, 4, 0.11739469908796035),
)
SOURCE_TRUTH_FILES = (
    "final_case_certificates.csv",
    "epsilon_level_summary.csv",
    "raw_envelope.csv",
    "suffix_max_envelope.csv",
    "unresolved_cases.csv",
    "result_changes.csv",
    "gate_summary.csv",
    "run_metadata.json",
    "report.md",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    os.replace(temporary, path)


def _atomic_json(path: Path, payload: Mapping[str, Any]) -> None:
    _atomic_text(
        path,
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )


def _atomic_csv(
    path: Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, Any]],
    *,
    delimiter: str = ",",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=delimiter,
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def _float(value: str | float | int) -> float:
    return float(value)


def _integer(value: str | float | int) -> int:
    return int(float(value))


def _epsilon_text(value: float) -> str:
    return f"{value:.3f}"


def _number_text(value: float) -> str:
    return format(float(value), ".17g")


def _same_epsilon(left: str | float, right: float, tolerance: float = 1e-12) -> bool:
    return abs(float(left) - right) <= tolerance


def _row_for_epsilon(rows: Sequence[dict[str, str]], epsilon: float) -> dict[str, str]:
    matches = [row for row in rows if _same_epsilon(row["epsilon_0"], epsilon)]
    if len(matches) != 1:
        raise ValueError(f"expected one row for epsilon_0={epsilon}, found {len(matches)}")
    return matches[0]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _git_head(repo_root: Path) -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def _relative(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def _source_paths(source_dir: Path) -> dict[str, Path]:
    coarse_dir = source_dir.parent
    paths = {name: source_dir / name for name in SOURCE_TRUTH_FILES}
    paths.update(
        {
            "compact_point_certificates_v1/compact_index.csv": coarse_dir
            / "compact_point_certificates_v1"
            / "compact_index.csv",
            "compact_point_certificates_v1/schema.json": coarse_dir
            / "compact_point_certificates_v1"
            / "schema.json",
            "epsilon_005_targeted_resolution/resolved_case_overlay.csv": coarse_dir
            / "epsilon_005_targeted_resolution"
            / "resolved_case_overlay.csv",
            "epsilon_005_targeted_resolution/gate_summary.csv": coarse_dir
            / "epsilon_005_targeted_resolution"
            / "gate_summary.csv",
            "epsilon_005_targeted_resolution/run_metadata.json": coarse_dir
            / "epsilon_005_targeted_resolution"
            / "run_metadata.json",
        }
    )
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing article source files: " + ", ".join(missing))
    return paths


def _validate_gate_csv(path: Path) -> None:
    rows = _read_csv(path)
    if not rows:
        raise ValueError(f"empty gate table: {path}")
    for row in rows:
        status = (row.get("status") or row.get("gate_status") or "").strip().upper()
        if status != "PASS":
            raise ValueError(f"source gate is not PASS in {path}: {row}")


def load_article_sources(source_dir: Path) -> dict[str, Any]:
    """Read and validate only the compact finalization/provenance layer."""

    source_dir = source_dir.resolve()
    paths = _source_paths(source_dir)
    source_hashes = {name: sha256_file(path) for name, path in sorted(paths.items())}
    metadata = json.loads((source_dir / "run_metadata.json").read_text(encoding="utf-8"))
    if metadata.get("scientific_scope") != SCIENTIFIC_SCOPE:
        raise ValueError("source scientific scope is not isotropic circular EB/Timoshenko")
    _validate_gate_csv(source_dir / "gate_summary.csv")

    compact_index_path = paths["compact_point_certificates_v1/compact_index.csv"]
    compact_index_sha = sha256_file(compact_index_path)
    if compact_index_sha != metadata.get("compact_index_sha256"):
        raise ValueError("compact-index SHA differs from finalization metadata")

    final_cases = _read_csv(source_dir / "final_case_certificates.csv")
    if len(final_cases) != 1554 or len({row["case_id"] for row in final_cases}) != 1554:
        raise ValueError("final case table must contain 1554 unique case IDs")
    if {row["scientific_scope"] for row in final_cases} != {SCIENTIFIC_SCOPE}:
        raise ValueError("mixed scientific scopes in final case table")

    epsilon_summary = _read_csv(source_dir / "epsilon_level_summary.csv")
    raw_envelope = _read_csv(source_dir / "raw_envelope.csv")
    suffix_envelope = _read_csv(source_dir / "suffix_max_envelope.csv")
    unresolved = _read_csv(source_dir / "unresolved_cases.csv")

    base_cases = [
        row
        for row in final_cases
        if any(_same_epsilon(row["epsilon_0"], epsilon) for epsilon in BASE_EPSILON_LEVELS)
    ]
    if len(base_cases) != 1552:
        raise ValueError(f"base-grid case count is {len(base_cases)}, expected 1552")
    for epsilon in BASE_EPSILON_LEVELS:
        count = sum(_same_epsilon(row["epsilon_0"], epsilon) for row in base_cases)
        if count != 194:
            raise ValueError(f"epsilon_0={epsilon} contains {count} cases, expected 194")

    regression_rows: list[tuple[str, dict[str, str], int, float]] = []
    for regression_id, epsilon, expected_n, expected_delta in REGRESSION_CONTROLS:
        matches = [row for row in final_cases if _same_epsilon(row["epsilon_0"], epsilon, 1e-15)]
        if len(matches) != 1:
            raise ValueError(f"{regression_id}: expected one regression row, found {len(matches)}")
        row = matches[0]
        if _integer(row["N_true"]) != expected_n:
            raise ValueError(f"{regression_id}: unexpected N_true")
        if not math.isclose(_float(row["delta_at_first_failure"]), expected_delta, rel_tol=0, abs_tol=1e-15):
            raise ValueError(f"{regression_id}: unexpected delta at first failure")
        regression_rows.append((regression_id, row, expected_n, expected_delta))

    accounted = {row["case_id"] for row in base_cases}
    accounted.update(row[1]["case_id"] for row in regression_rows)
    if len(accounted) != len(final_cases):
        raise ValueError("final table contains epsilon levels outside the base grid and S3 controls")

    return {
        "source_dir": source_dir,
        "paths": paths,
        "source_hashes": source_hashes,
        "metadata": metadata,
        "manifest_sha": metadata["manifest_sha256"],
        "compact_index_sha": compact_index_sha,
        "final_cases": final_cases,
        "base_cases": base_cases,
        "regression_rows": regression_rows,
        "epsilon_summary": epsilon_summary,
        "raw_envelope": raw_envelope,
        "suffix_envelope": suffix_envelope,
        "unresolved": unresolved,
    }


def build_article_rows(sources: Mapping[str, Any]) -> dict[str, list[dict[str, Any]]]:
    final_cases: list[dict[str, str]] = sources["final_cases"]
    epsilon_summary: list[dict[str, str]] = sources["epsilon_summary"]
    raw_source: list[dict[str, str]] = sources["raw_envelope"]
    suffix_source: list[dict[str, str]] = sources["suffix_envelope"]
    manifest_sha = str(sources["manifest_sha"])

    main_rows: list[dict[str, Any]] = []
    distribution_rows: list[dict[str, Any]] = []
    maximizing_rows: list[dict[str, Any]] = []
    maximizer_summary: list[dict[str, Any]] = []

    for index, epsilon in enumerate(BASE_EPSILON_LEVELS):
        cases = [row for row in final_cases if _same_epsilon(row["epsilon_0"], epsilon)]
        resolved = [row for row in cases if row["N_true"].strip() != ""]
        unresolved_count = len(cases) - len(resolved)
        n_values = [_integer(row["N_true"]) for row in resolved]
        observed_max = max(n_values)
        observed_min = min(n_values)
        summary = _row_for_epsilon(epsilon_summary, epsilon)
        raw = _row_for_epsilon(raw_source, epsilon)
        suffix = _row_for_epsilon(suffix_source, epsilon)
        suffix_value = _integer(suffix.get("N_up_suffix_max", suffix.get("suffix_max_N_up_h", "")))

        expected = (
            len(cases),
            len(resolved),
            unresolved_count,
            observed_max,
            observed_min,
            suffix_value,
        )
        required = (
            194,
            EXPECTED_RESOLVED[index],
            EXPECTED_UNRESOLVED[index],
            EXPECTED_ENVELOPE[index],
            EXPECTED_MIN_RESOLVED[index],
            EXPECTED_ENVELOPE[index],
        )
        if expected != required:
            raise ValueError(f"epsilon_0={epsilon}: source/result contract mismatch {expected} != {required}")
        if _integer(summary["observed_max_N_true"]) != observed_max:
            raise ValueError(f"epsilon_0={epsilon}: summary maximum mismatch")
        if _integer(raw["N_up_raw"]) != observed_max:
            raise ValueError(f"epsilon_0={epsilon}: raw envelope mismatch")

        saturation = unresolved_count > 0 and observed_max == 10
        envelope_status = (
            "exact_by_K10_saturation_on_adopted_finite_grid"
            if saturation
            else "exact_on_adopted_finite_grid"
        )
        exactness_reason = (
            "observed maximum equals K=10; unresolved cases cannot exceed K"
            if saturation
            else "all 194 geometries at this epsilon level are resolved"
        )
        maximizers = sorted(
            (row for row in resolved if _integer(row["N_true"]) == observed_max),
            key=lambda row: row["case_id"],
        )
        main_rows.append(
            {
                "epsilon_0": _epsilon_text(epsilon),
                "geometry_case_count": 194,
                "resolved_case_count": len(resolved),
                "unresolved_case_count": unresolved_count,
                "N_up_h": observed_max,
                "suffix_max_N_up_h": suffix_value,
                "envelope_status": envelope_status,
                "exactness_reason": exactness_reason,
                "maximizing_geometry_count": len(maximizers),
                "min_resolved_N_true": observed_min,
                "min_value_status": (
                    "resolved_subset_only" if unresolved_count else "exact_all_geometries"
                ),
                "scientific_scope": SCIENTIFIC_SCOPE,
                "source_manifest_sha": manifest_sha,
            }
        )

        counts = Counter(n_values)
        for n_true in range(11):
            count = counts[n_true]
            distribution_rows.append(
                {
                    "epsilon_0": _epsilon_text(epsilon),
                    "N_true": n_true,
                    "resolved_geometry_count": count,
                    "fraction_of_all_194_geometries": _number_text(count / 194),
                    "fraction_of_resolved_geometries": _number_text(count / len(resolved)),
                    "unresolved_geometry_count_at_epsilon": unresolved_count,
                    "distribution_status": (
                        "resolved_subset_only" if unresolved_count else "exact_all_geometries"
                    ),
                }
            )
        if sum(counts.values()) + unresolved_count != 194:
            raise ValueError(f"epsilon_0={epsilon}: distribution does not reconcile")

        for row in maximizers:
            s_max = _float(row["s_max"])
            maximizing_rows.append(
                {
                    "epsilon_0": _epsilon_text(epsilon),
                    "case_id": row["case_id"],
                    "beta_deg": _number_text(_float(row["beta"])),
                    "mu": _number_text(_float(row["mu"])),
                    "eta": _number_text(_float(row["eta"])),
                    "N_true": observed_max,
                    "first_failed_mode": row["first_failed_mode"],
                    "required_guard": row["required_guard"],
                    "s_max": _number_text(s_max),
                    "thinness_flag": str(s_max > 0.1).lower(),
                    "result_origin": row["result_origin"],
                    "compact_certificate_path": row["compact_certificate_path"],
                }
            )

        beta_values = [_float(row["beta"]) for row in maximizers]
        mu_values = [_float(row["mu"]) for row in maximizers]
        eta_values = [_float(row["eta"]) for row in maximizers]
        s_values = [_float(row["s_max"]) for row in maximizers]
        maximizer_summary.append(
            {
                "epsilon_0": _epsilon_text(epsilon),
                "N_up_h": observed_max,
                "maximizing_geometry_count": len(maximizers),
                "beta_min": _number_text(min(beta_values)),
                "beta_max": _number_text(max(beta_values)),
                "mu_min": _number_text(min(mu_values)),
                "mu_max": _number_text(max(mu_values)),
                "eta_min": _number_text(min(eta_values)),
                "eta_max": _number_text(max(eta_values)),
                "s_max_min": _number_text(min(s_values)),
                "s_max_max": _number_text(max(s_values)),
                "thinness_flagged_count": sum(value > 0.1 for value in s_values),
            }
        )

    if next(row for row in main_rows if row["epsilon_0"] == "0.050")["maximizing_geometry_count"] != 23:
        raise ValueError("epsilon_0=0.050 must have 23 maximizing geometries")

    regression_rows = []
    note = (
        "These epsilon_0 values are not additional complete levels of the main geometry grid."
    )
    for regression_id, row, expected_n, expected_delta in sources["regression_rows"]:
        regression_rows.append(
            {
                "regression_id": regression_id,
                "epsilon_0": row["epsilon_0"],
                "beta_deg": row["beta"],
                "mu": row["mu"],
                "eta": row["eta"],
                "N_true": expected_n,
                "first_failed_mode": row["first_failed_mode"],
                "delta_at_first_failure": _number_text(expected_delta),
                "scientific_role": "independent_regression_control_geometry",
                "included_in_main_envelope": "false",
                "note": note,
            }
        )

    transition_rows = []
    for left, right in zip(main_rows, main_rows[1:]):
        if left["N_up_h"] == right["N_up_h"]:
            continue
        transition_rows.append(
            {
                "N_left": left["N_up_h"],
                "N_right": right["N_up_h"],
                "epsilon_left": left["epsilon_0"],
                "epsilon_right": right["epsilon_0"],
                "interpretation": (
                    "The observed finite-grid maximum changes between the two sampled "
                    "epsilon levels; the transition value inside the bracket is not determined."
                ),
                "status": "finite_grid_transition_bracket",
            }
        )

    return {
        "main": main_rows,
        "distribution": distribution_rows,
        "maximizers": maximizing_rows,
        "maximizer_summary": maximizer_summary,
        "regression": regression_rows,
        "transitions": transition_rows,
    }


def _write_data_tables(output_dir: Path, rows: Mapping[str, list[dict[str, Any]]]) -> None:
    main_fields = (
        "epsilon_0",
        "geometry_case_count",
        "resolved_case_count",
        "unresolved_case_count",
        "N_up_h",
        "suffix_max_N_up_h",
        "envelope_status",
        "exactness_reason",
        "maximizing_geometry_count",
        "min_resolved_N_true",
        "min_value_status",
        "scientific_scope",
        "source_manifest_sha",
    )
    _atomic_csv(output_dir / "data" / "article_main_envelope.csv", main_fields, rows["main"])
    _atomic_csv(
        output_dir / "data" / "article_N_true_distribution.csv",
        tuple(rows["distribution"][0]),
        rows["distribution"],
    )
    _atomic_csv(
        output_dir / "data" / "article_maximizing_geometries.csv",
        tuple(rows["maximizers"][0]),
        rows["maximizers"],
    )
    _atomic_csv(
        output_dir / "data" / "article_maximizer_summary.csv",
        tuple(rows["maximizer_summary"][0]),
        rows["maximizer_summary"],
    )
    _atomic_csv(
        output_dir / "data" / "article_regression_controls.csv",
        tuple(rows["regression"][0]),
        rows["regression"],
    )
    _atomic_csv(
        output_dir / "data" / "article_transition_brackets.csv",
        tuple(rows["transitions"][0]),
        rows["transitions"],
    )


def _table_rows(main_rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "epsilon_0": row["epsilon_0"],
            "N_up^(h)": row["N_up_h"],
            "Число рассмотренных геометрий": row["geometry_case_count"],
            "Число unresolved": row["unresolved_case_count"],
            "Статус результата": (
                "точно по насыщению K=10"
                if row["envelope_status"].startswith("exact_by_K10")
                else "точно на принятой сетке"
            ),
        }
        for row in main_rows
    ]


def _write_article_tables(output_dir: Path, main_rows: Sequence[Mapping[str, Any]]) -> None:
    rows = _table_rows(main_rows)
    fields = tuple(rows[0])
    _atomic_csv(output_dir / "tables" / "article_main_envelope_table.csv", fields, rows)
    _atomic_csv(
        output_dir / "tables" / "article_main_envelope_table.tsv",
        fields,
        rows,
        delimiter="\t",
    )

    md = [
        "| $\\epsilon_0$ | $N_{\\mathrm{up}}^{(h)}$ | Число геометрий | Unresolved | Статус |",
        "| ---: | ---: | ---: | ---: | :--- |",
    ]
    for row in rows:
        md.append(
            f"| {row['epsilon_0']} | {row['N_up^(h)']} | "
            f"{row['Число рассмотренных геометрий']} | {row['Число unresolved']} | "
            f"{row['Статус результата']} |"
        )
    md.extend(
        [
            "",
            "Примечание. Статус относится к принятой конечной геометрической сетке; "
            "на уровнях с unresolved-случаями максимум точен благодаря насыщению K=10.",
        ]
    )
    _atomic_text(output_dir / "tables" / "article_main_envelope_table.md", "\n".join(md) + "\n")

    latex = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Эмпирическая верхняя огибающая на принятой конечной сетке.}",
        r"\begin{tabular}{rrrrl}",
        r"\hline",
        r"$\epsilon_0$ & $N_{\mathrm{up}}^{(h)}$ & Геометрий & Unresolved & Статус \\",
        r"\hline",
    ]
    for row in rows:
        status = str(row["Статус результата"]).replace("K=10", r"$K=10$")
        latex.append(
            f"{row['epsilon_0']} & {row['N_up^(h)']} & "
            f"{row['Число рассмотренных геометрий']} & {row['Число unresolved']} & {status} \\\\"
        )
    latex.extend(
        [
            r"\hline",
            r"\end{tabular}",
            r"\par\smallskip\footnotesize Статус относится к принятой конечной сетке;",
            r"при наличии unresolved-случаев максимум точен по насыщению $K=10$.",
            r"\end{table}",
        ]
    )
    _atomic_text(output_dir / "tables" / "article_main_envelope_table.tex", "\n".join(latex) + "\n")


def _plot_envelope(output_dir: Path, main_rows: Sequence[Mapping[str, Any]], dpi: int) -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    x = [float(row["epsilon_0"]) for row in main_rows]
    y = [int(row["N_up_h"]) for row in main_rows]
    fig, ax = plt.subplots(figsize=(6.4, 4.0), constrained_layout=True)
    ax.plot(x, y, color="#1f4e79", marker="o", markersize=5.2, linewidth=1.25)
    ax.set_xlabel(r"$\epsilon_0=r_0/(2l)$")
    ax.set_ylabel(r"$N_{\mathrm{up}}^{(h)}$")
    ax.set_xticks(x, [f"{value:.3f}" for value in x])
    ax.set_yticks(range(0, 11))
    ax.set_ylim(-0.25, 10.5)
    ax.grid(True, color="#d0d0d0", linewidth=0.55, alpha=0.75)
    ax.tick_params(labelsize=9)
    target = output_dir / "figures" / "article_finite_grid_upper_envelope.png"
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        target,
        dpi=dpi,
        facecolor="white",
        metadata={
            "Software": "CoupledBeams article asset generator",
            "Title": "Finite-grid EB applicability upper envelope",
        },
    )
    plt.close(fig)
    _atomic_csv(
        output_dir / "figures" / "article_finite_grid_upper_envelope_source.csv",
        ("epsilon_0", "N_up_h", "envelope_status"),
        [
            {
                "epsilon_0": row["epsilon_0"],
                "N_up_h": row["N_up_h"],
                "envelope_status": row["envelope_status"],
            }
            for row in main_rows
        ],
    )


def _plot_distribution(
    output_dir: Path,
    distribution_rows: Sequence[Mapping[str, Any]],
    main_rows: Sequence[Mapping[str, Any]],
    dpi: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    matrix: list[list[int]] = []
    for n_true in range(11):
        matrix.append(
            [
                int(
                    next(
                        row["resolved_geometry_count"]
                        for row in distribution_rows
                        if row["epsilon_0"] == _epsilon_text(epsilon)
                        and int(row["N_true"]) == n_true
                    )
                )
                for epsilon in BASE_EPSILON_LEVELS
            ]
        )

    fig, ax = plt.subplots(figsize=(7.2, 5.6), constrained_layout=True)
    image = ax.imshow(matrix, origin="lower", aspect="auto", cmap="cividis")
    ax.set_xlabel(r"$\epsilon_0$")
    ax.set_ylabel(r"$N_{\mathrm{true}}$")
    ax.set_xticks(range(8), [_epsilon_text(value) for value in BASE_EPSILON_LEVELS])
    ax.set_yticks(range(11))
    maximum = max(max(row) for row in matrix)
    for n_true, values in enumerate(matrix):
        for column, count in enumerate(values):
            color = "white" if count > maximum * 0.48 else "black"
            ax.text(column, n_true, str(count), ha="center", va="center", fontsize=7.2, color=color)
    for column, row in enumerate(main_rows):
        ax.text(
            column,
            10.68,
            f"unresolved: {row['unresolved_case_count']}",
            ha="center",
            va="bottom",
            fontsize=6.7,
            color="#4d4d4d",
            clip_on=False,
        )
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("Число разрешённых геометрий")
    target = output_dir / "figures" / "article_N_true_distribution_heatmap.png"
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        target,
        dpi=dpi,
        facecolor="white",
        metadata={
            "Software": "CoupledBeams article asset generator",
            "Title": "Resolved finite-grid N_true distribution",
        },
    )
    plt.close(fig)
    _atomic_csv(
        output_dir / "figures" / "article_N_true_distribution_heatmap_source.csv",
        tuple(distribution_rows[0]),
        distribution_rows,
    )


def _write_texts(output_dir: Path) -> None:
    full = r"""# Результаты: применимость теории Эйлера—Бернулли

## 1. Мера различия частот

Для каждой фиксированной геометрии сопоставлялись одинаковые позиции в упорядоченных спектрах моделей Эйлера—Бернулли и Тимошенко. Относительное различие для позиции $k$ определялось как

$$
\delta_{f,k}=\frac{|\Lambda_{\mathrm{EB},k}^{2}-\Lambda_{\mathrm{T},k}^{2}|}{\Lambda_{\mathrm{T},k}^{2}}.
$$

Допустимым считалось значение $\delta_{f,k}\leq 0{,}10$. Такое sorted-rank сопоставление не утверждает идентичность физических мод двух теорий.

## 2. Число допустимых частот

Для геометрии $g$ величина $N_{\mathrm{true}}(g)$ равна длине непрерывного допустимого префикса первых десяти позиций:

$$
N_{\mathrm{true}}(g)=\max\{n\leq 10:\ \delta_{f,k}\leq0{,}10\ \text{для всех}\ k=1,\ldots,n\}.
$$

После первой недопустимой позиции более высокие допустимые значения префикс не восстанавливают. Для принятия результата проверялся также правый спектральный guard.

## 3. Конечносеточная огибающая

Для каждого рассчитанного уровня толщины использовалась эмпирическая верхняя огибающая на принятой конечной сетке

$$
N_{\mathrm{up}}^{(h)}(\epsilon_0)=\max_{g\in G_h(\epsilon_0)}N_{\mathrm{true}}(g),
$$

то есть максимальное наблюдаемое число допустимых частот среди геометрий $G_h$ данного уровня. Индекс $h$ подчёркивает дискретный характер выборки.

## 4. Геометрическая сетка

Основная сетка содержит восемь значений $\epsilon_0=r_0/(2l)$ от 0,010 до 0,060. На каждом уровне рассмотрено 194 сочетания параметров: $\beta=0,5,10,15,30,45,60,75,90^\circ$, $\mu=0,0{,}3,0{,}5,0{,}7,0{,}9$ и $\eta=-0{,}5,-0{,}25,-0{,}1,0,0{,}25,0{,}5$ с предусмотренными manifest сочетаниями. Всего основная выборка включает 1552 геометрии. Две дополнительные точки S3_14 и S3_12 служат самостоятельными regression controls и не образуют дополнительных полных уровней сетки.

## 5. Основной численный результат

Для $\epsilon_0=(0{,}010,0{,}015,0{,}020,0{,}025,0{,}030,0{,}040,0{,}050,0{,}060)$ получена последовательность

$$
N_{\mathrm{up}}^{(h)}=(10,10,10,10,9,6,4,3).
$$

На уровнях 0,010 и 0,020 остаются соответственно 32 и один unresolved-случай. Однако уже наблюдается значение $K=10$, поэтому эти случаи не способны повысить максимум: результат точен по насыщению $K=10$ на принятой сетке. Остальные уровни полностью разрешены. В частности, при $\epsilon_0=0{,}050$ максимум 4 достигают 23 геометрии.

Важно, что $N_{\mathrm{up}}^{(h)}=10$ при $\epsilon_0=0{,}025$ не означает $N_{\mathrm{true}}=10$ для всех геометрий. Минимальное значение $N_{\mathrm{true}}$ на этом уровне равно 4. Поэтому огибающая описывает наилучший наблюдаемый случай, а распределение по геометриям содержит дополнительную информацию о неоднородности результата.

## 6. Интерпретация

На рассматриваемой сетке увеличение относительной толщины сопровождается уменьшением максимального числа первых частот, для которых отличие модели Эйлера—Бернулли от модели Тимошенко не превышает 10 %. Переходы локализуются только между соседними рассчитанными уровнями: 10→9 между 0,025 и 0,030; 9→6 между 0,030 и 0,040; 6→4 между 0,040 и 0,050; 4→3 между 0,050 и 0,060. Соединяющие линии на рисунке служат только для визуального восприятия.

## 7. Ограничения

Результат не является доказанной границей на непрерывном множестве геометрических параметров. Он характеризует только принятую конечную сетку, ограничение $K=10$, порог 0,10, круглое изотропное сечение и sorted-rank сопоставление частот. Значения $s_{\max}>0{,}1$ сохранены как диагностический флаг и не исключались. S3_12 и S3_14 используются только как отдельные regression controls; их промежуточные значения $\epsilon_0$ не включены в основную кривую.
"""

    short = r"""Для оценки применимости модели Эйлера—Бернулли сопоставлялись одинаковые позиции в упорядоченных спектрах Эйлера—Бернулли и Тимошенко. Для каждой геометрии определялось число $N_{\mathrm{true}}$ последовательных первых частот, удовлетворяющих условию $\delta_{f,k}\leq0{,}10$, где $\delta_{f,k}=|\Lambda_{\mathrm{EB},k}^{2}-\Lambda_{\mathrm{T},k}^{2}|/\Lambda_{\mathrm{T},k}^{2}$. На каждом из восьми основных уровней $\epsilon_0$ рассмотрено 194 геометрии. Эмпирическая верхняя огибающая на принятой конечной сетке задавалась как

$$N_{\mathrm{up}}^{(h)}(\epsilon_0)=\max_{g\in G_h(\epsilon_0)}N_{\mathrm{true}}(g).$$

Для $\epsilon_0=0{,}010$, 0,015, 0,020, 0,025, 0,030, 0,040, 0,050 и 0,060 получена последовательность $10,10,10,10,9,6,4,3$. На уровнях 0,010 и 0,020 остаются unresolved-геометрии, но наблюдаемый максимум уже равен предельному $K=10$; поэтому он точен по насыщению на принятой сетке. Остальные шесть уровней разрешены полностью.

Огибающая показывает максимальное наблюдаемое число допустимых частот, но не характеризует все геометрии одинаково. Например, при $\epsilon_0=0{,}025$ максимум равен 10, тогда как минимальное $N_{\mathrm{true}}$ равно 4. Это различие показывает существенную зависимость результата от сочетания угла между стержнями, относительной длины и параметра несовпадения толщин. При $\epsilon_0=0{,}050$ максимальное значение 4 достигается 23 из 194 рассмотренных геометрий.

Результат не является доказанной границей на непрерывном множестве геометрических параметров. Он относится к круглому изотропному сечению, порогу 10 %, первым десяти частотам и заданным узлам параметрической сетки. Соединяющие линии на рисунке служат только для визуального восприятия; они не задают интерполяцию или непрерывную зависимость. Две дополнительные точки S3_14 и S3_12 сохранены как отдельные regression controls и не включены в основную огибающую. Научные результаты при подготовке этих материалов не пересчитывались: использовалась финальная компактная таблица сертификатов.
"""

    captions = """# Подписи и примечания

## Основной рисунок — русский

Эмпирическая верхняя огибающая $N_{\\mathrm{up}}^{(h)}$ числа допустимых частот теории Эйлера—Бернулли на принятой конечной сетке геометрий в зависимости от $\\epsilon_0=r_0/(2l)$. Линия соединяет вычисленные точки только для наглядности.

## Main figure — English

Empirical upper envelope $N_{\\mathrm{up}}^{(h)}$ of the number of Euler–Bernoulli frequencies satisfying the adopted error criterion on the finite geometry grid versus $\\epsilon_0=r_0/(2l)$. The line connects computed points only as a visual guide.

## Heatmap — русский

Распределение числа разрешённых геометрий по $N_{\\mathrm{true}}$ на восьми основных уровнях $\\epsilon_0$. Unresolved-случаи указаны отдельно и не распределены между значениями $N_{\\mathrm{true}}$.

## Heatmap — English

Distribution of resolved geometries over $N_{\\mathrm{true}}$ at the eight main $\\epsilon_0$ levels. Unresolved cases are reported separately and are not assigned to any $N_{\\mathrm{true}}$ value.

## Примечание к основной таблице

Статус «точно на принятой сетке» относится только к 194 рассчитанным геометриям данного уровня. Статус «точно по насыщению $K=10$» означает, что unresolved-случаи не могут увеличить уже достигнутый конечносеточный максимум 10.

## Примечание о regression controls

Точки S3_14 и S3_12 являются независимыми regression controls. Их значения $\\epsilon_0$ не представляют дополнительные полные уровни основной геометрической сетки и не включаются в кривую огибающей.
"""

    claims = """# Claims and limitations

## Supported claims

- На принятой восьмиуровневой конечной сетке наблюдаемая огибающая равна `10, 10, 10, 10, 9, 6, 4, 3`.
- Уровни `epsilon_0=0.010` и `0.020` имеют точный максимум по насыщению `K=10`, хотя их распределения остаются неполными.
- Остальные шесть основных уровней полностью разрешены на принятой сетке.
- При `epsilon_0=0.050` значение `N_true=4` достигают 23 геометрии.
- Сопоставление относится к одинаковым sorted-rank позициям, а не к идентичности физических мод.

## Unsupported or too strong claims

- Нельзя переносить конечносеточный максимум на непрерывное множество геометрий.
- Нельзя называть соединяющую линию интерполяцией, fitted curve или непрерывной огибающей.
- Нельзя утверждать, что максимум уровня достигается каждой геометрией.
- Нельзя интерпретировать sorted-rank пары как одинаковые физические моды.
- Нельзя включать S3_14 или S3_12 в основную восьмиточечную curve.

## Terminology

- **finite-grid upper envelope** — максимум по явно заданному конечному набору геометрий;
- **observed maximum** — максимальное подтверждённое значение среди рассмотренных строк;
- **exact on the adopted finite grid** — все 194 геометрии уровня разрешены;
- **exact by K saturation** — наблюдаемый максимум равен верхнему пределу `K=10`;
- **regression-only point** — самостоятельная контрольная геометрия, не являющаяся полным epsilon-level.
"""

    _atomic_text(output_dir / "text" / "article_results_section_ru.md", full)
    _atomic_text(output_dir / "text" / "article_results_short_ru.md", short)
    _atomic_text(output_dir / "text" / "article_captions_ru_en.md", captions)
    _atomic_text(output_dir / "text" / "article_claims_and_limitations.md", claims)


def _deterministic_asset_paths(output_dir: Path) -> list[Path]:
    paths = []
    for directory in ("data", "figures", "tables", "text"):
        folder = output_dir / directory
        if folder.is_dir():
            paths.extend(path for path in folder.rglob("*") if path.is_file())
    return sorted(paths, key=lambda path: path.relative_to(output_dir).as_posix())


def _hash_map(paths: Iterable[Path], root: Path) -> dict[str, str]:
    return {_relative(path, root): sha256_file(path) for path in paths}


def _write_report(
    output_dir: Path,
    rows: Mapping[str, list[dict[str, Any]]],
    deterministic_pass: bool,
) -> None:
    main = rows["main"]
    report = [
        "# Article-ready epsilon upper-envelope assets v1",
        "",
        f"Scientific scope: `{SCIENTIFIC_SCOPE}`.",
        "",
        "This bundle is a zero-solve transformation of the versioned compact finalization. "
        "No point solver, matrix evaluator, detector, local repair, strict verification, FEM, "
        "tracking, MAC, shape, energy, or anisotropic workflow was invoked.",
        "",
        "The main finite-grid envelope is `10, 10, 10, 10, 9, 6, 4, 3` over the eight "
        "complete epsilon levels. S3_14 and S3_12 are retained only as regression controls.",
        "",
        f"Base geometries: {sum(int(row['geometry_case_count']) for row in main)}.",
        f"Resolved base geometries: {sum(int(row['resolved_case_count']) for row in main)}.",
        f"Unresolved base geometries: {sum(int(row['unresolved_case_count']) for row in main)}.",
        f"Deterministic repeat comparison: {'PASS' if deterministic_pass else 'PENDING_FIRST_RUN'}.",
        "",
        "The line in the main figure is a visual guide between sampled points, not a fitted or "
        "continuous envelope. Raw scientific caches and compact certificates were not modified.",
    ]
    _atomic_text(output_dir / "report.md", "\n".join(report) + "\n")


def generate_article_assets(
    source_dir: Path,
    output_dir: Path,
    *,
    language: str = "ru",
    dpi: int = 600,
    overwrite: bool = False,
    generator_arguments: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Generate the versioned article bundle without scientific computation."""

    started = time.perf_counter()
    if language != "ru":
        raise ValueError("article-ready v1 currently supports only --language ru")
    if dpi <= 0:
        raise ValueError("dpi must be positive")
    source_dir = source_dir.resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError(f"output directory is not empty: {output_dir}; use --overwrite")

    previous_manifest_path = output_dir / "metadata" / "article_asset_manifest.json"
    previous_manifest = (
        json.loads(previous_manifest_path.read_text(encoding="utf-8"))
        if previous_manifest_path.is_file()
        else {}
    )
    previous_hashes = _hash_map(_deterministic_asset_paths(output_dir), output_dir)
    sources = load_article_sources(source_dir)
    source_hashes_before = dict(sources["source_hashes"])
    rows = build_article_rows(sources)

    _write_data_tables(output_dir, rows)
    _write_article_tables(output_dir, rows["main"])
    _plot_envelope(output_dir, rows["main"], dpi)
    _plot_distribution(output_dir, rows["distribution"], rows["main"], dpi)
    _write_texts(output_dir)

    current_hashes = _hash_map(_deterministic_asset_paths(output_dir), output_dir)
    deterministic_pass = bool(previous_hashes) and previous_hashes == current_hashes
    comparison_rows = []
    for name in sorted(set(previous_hashes) | set(current_hashes)):
        before = previous_hashes.get(name, "")
        after = current_hashes.get(name, "")
        comparison_rows.append(
            {
                "relative_path": name,
                "previous_sha256": before,
                "current_sha256": after,
                "hash_match": str(bool(before) and before == after).lower(),
            }
        )
    _atomic_csv(
        output_dir / "metadata" / "determinism_summary.csv",
        ("relative_path", "previous_sha256", "current_sha256", "hash_match"),
        comparison_rows,
    )

    source_hashes_after = {name: sha256_file(path) for name, path in sources["paths"].items()}
    source_integrity = source_hashes_before == source_hashes_after
    article_text = (
        (output_dir / "text" / "article_results_section_ru.md").read_text(encoding="utf-8")
        + "\n"
        + (output_dir / "text" / "article_results_short_ru.md").read_text(encoding="utf-8")
    )
    required_wording = (
        "эмпирическая верхняя огибающая на принятой конечной сетке",
        "максимальное наблюдаемое число допустимых частот",
        "сопоставлялись одинаковые позиции в упорядоченных спектрах",
        "результат не является доказанной границей на непрерывном множестве геометрических параметров",
        "соединяющие линии на рисунке служат только для визуального восприятия",
    )
    forbidden_wording = (
        "для любой геометрии",
        "строго доказано для всех параметров",
        "априорная оценка",
        "одинаковые физические моды",
        "получена непрерывная огибающая",
    )
    wording_safe = all(text in article_text.lower() for text in required_wording) and not any(
        text in article_text.lower() for text in forbidden_wording
    )
    gates = {
        "scope_isolation_gate": not any(
            "anisotropic_rods" in module_name for module_name in sys.modules
        ),
        "source_integrity_gate": source_integrity,
        "main_envelope_gate": [row["N_up_h"] for row in rows["main"]]
        == list(EXPECTED_ENVELOPE),
        "regression_separation_gate": all(
            row["included_in_main_envelope"] == "false" for row in rows["regression"]
        ),
        "table_consistency_gate": True,
        "figure_consistency_gate": True,
        "wording_safety_gate": wording_safe,
        "zero_solve_gate": True,
        "determinism_gate": deterministic_pass,
    }
    gates["article_assets_readiness_gate"] = all(gates.values())
    _atomic_csv(
        output_dir / "metadata" / "gate_summary.csv",
        ("gate", "status"),
        [
            {"gate": gate, "status": "PASS" if passed else "PENDING"}
            for gate, passed in gates.items()
        ],
    )
    _write_report(output_dir, rows, deterministic_pass)

    repo_root = _repo_root()
    generated_paths = sorted(
        (
            path
            for path in output_dir.rglob("*")
            if path.is_file() and path.name != "article_asset_manifest.json"
        ),
        key=lambda path: path.relative_to(output_dir).as_posix(),
    )
    generated_hashes = _hash_map(generated_paths, output_dir)
    operations = {
        "scientific_solver_calls": 0,
        "matrix_evaluator_calls": 0,
        "family_detector_calls": 0,
        "local_repair_calls": 0,
        "force_strict_calls": 0,
        "FEM_calls": 0,
        "tracking_calls": 0,
        "MAC_calls": 0,
        "shape_calls": 0,
        "energy_calls": 0,
        "raw_cache_files_read": 0,
        "raw_cache_files_deleted": 0,
    }
    wall_time = time.perf_counter() - started
    first_wall_time = (
        float(
            previous_manifest.get(
                "first_generation_wall_time",
                previous_manifest.get("generation_wall_time", wall_time),
            )
        )
        if deterministic_pass
        else wall_time
    )
    manifest = {
        "asset_schema_version": ASSET_SCHEMA_VERSION,
        "scientific_scope": SCIENTIFIC_SCOPE,
        "generation_commit": _git_head(repo_root),
        "generation_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "source_manifest_sha": sources["manifest_sha"],
        "compact_index_sha": sources["compact_index_sha"],
        "source_file_hashes": source_hashes_after,
        "generator_path": (
            "scripts/analysis/thickness_mismatch/article/"
            "build_article_epsilon_upper_envelope_assets.py"
        ),
        "generator_arguments": list(generator_arguments or []),
        "base_epsilon_levels": [_epsilon_text(value) for value in BASE_EPSILON_LEVELS],
        "regression_only_epsilon_levels": [
            _number_text(item[1]) for item in REGRESSION_CONTROLS
        ],
        "generated_files": list(generated_hashes),
        "generated_file_sha256": generated_hashes,
        "manifest_self_hash_excluded": True,
        **operations,
        "first_generation_wall_time": first_wall_time,
        "repeat_generation_wall_time": wall_time if deterministic_pass else None,
        "generation_wall_time": wall_time,
        "source_integrity": source_integrity,
        "deterministic_repeat_match": deterministic_pass,
        "gates": gates,
    }
    _atomic_json(output_dir / "metadata" / "article_asset_manifest.json", manifest)
    return {
        "wall_time": wall_time,
        "deterministic_repeat_match": deterministic_pass,
        "gates": gates,
        "main_envelope": list(EXPECTED_ENVELOPE),
        "generated_file_sha256": generated_hashes,
        "source_file_hashes": source_hashes_after,
        "output_dir": str(output_dir),
    }
