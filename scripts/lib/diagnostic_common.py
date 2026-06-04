from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


def number_token(value: float) -> str:
    """Return a compact filename-safe decimal token."""
    return f"{float(value):.10g}".replace("-", "m").replace(".", "p")


def number_text(value: float | int | str | None, *, nonfinite_text: str = "") -> str:
    """Return a compact CSV/report value, with configurable non-finite text."""
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    value_f = float(value)
    if not np.isfinite(value_f):
        return str(nonfinite_text)
    return f"{value_f:.12g}"


def inclusive_grid(start: float, stop: float, step: float, *, step_name: str = "grid step") -> np.ndarray:
    start_f = float(start)
    stop_f = float(stop)
    step_f = float(step)
    if step_f <= 0.0:
        raise ValueError(f"{step_name} must be positive.")
    values = np.arange(start_f, stop_f + 0.5 * step_f, step_f, dtype=float)
    if values.size == 0:
        values = np.array([start_f, stop_f], dtype=float)
    if abs(float(values[0]) - start_f) > 1.0e-12:
        values = np.insert(values, 0, start_f)
    if abs(float(values[-1]) - stop_f) > 1.0e-12:
        values = np.append(values, stop_f)
    values[0] = start_f
    values[-1] = stop_f
    return np.unique(np.round(values, 12))


def ensure_output_dir(path: str | Path) -> Path:
    output_dir = Path(path)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def finite_or_nan(value: object) -> float:
    try:
        value_f = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return value_f if np.isfinite(value_f) else float("nan")


def write_dict_rows_csv(
    path: str | Path,
    rows: Iterable[Mapping[str, object]],
    fieldnames: Sequence[str],
) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return output


def parse_float_list(values: Sequence[str | float | int] | None) -> list[float]:
    if values is None:
        return []
    return [float(value) for value in values]


def finite_range(values: object) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return float("nan"), float("nan")
    return float(np.min(finite)), float(np.max(finite))
