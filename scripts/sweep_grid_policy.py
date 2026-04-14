from __future__ import annotations

import numpy as np


PRESENTATION_BETA_STEP = 2.0
PRESENTATION_MU_STEP = 0.03
ANALYSIS_BETA_STEP = 0.5
ANALYSIS_MU_STEP = 0.0125
LOCAL_BETA_REFINEMENT_STEP = 0.1
LOCAL_MU_REFINEMENT_STEP = 0.0025


def parameter_grid(start: float, stop: float, step: float) -> np.ndarray:
    start_f = float(start)
    stop_f = float(stop)
    step_f = float(step)
    if step_f <= 0.0:
        raise ValueError("Grid step must be positive.")

    values = np.arange(start_f, stop_f + 0.5 * step_f, step_f, dtype=float)
    if len(values) == 0 or abs(values[-1] - stop_f) > 1e-10:
        values = np.append(values, stop_f)
    values[0] = start_f
    values[-1] = stop_f
    return np.unique(np.round(values, 10))


def refine_parameter_grid(
    base_grid: np.ndarray,
    local_windows: tuple[tuple[float, float], ...] = (),
    local_step: float | None = None,
    anchor_points: tuple[float, ...] = (),
) -> np.ndarray:
    chunks = [np.asarray(base_grid, dtype=float)]

    if local_step is not None:
        for left, right in local_windows:
            chunks.append(parameter_grid(left, right, local_step))

    if anchor_points:
        chunks.append(np.asarray(anchor_points, dtype=float))

    return np.unique(np.round(np.concatenate(chunks), 10))


def presentation_beta_grid(
    start: float = 0.0,
    stop: float = 90.0,
    local_windows: tuple[tuple[float, float], ...] = (),
    anchor_points: tuple[float, ...] = (),
) -> np.ndarray:
    base = parameter_grid(start, stop, PRESENTATION_BETA_STEP)
    return refine_parameter_grid(
        base_grid=base,
        local_windows=local_windows,
        local_step=LOCAL_BETA_REFINEMENT_STEP if local_windows else None,
        anchor_points=anchor_points,
    )


def presentation_mu_grid(
    start: float = 0.0,
    stop: float = 0.9,
    local_windows: tuple[tuple[float, float], ...] = (),
    anchor_points: tuple[float, ...] = (),
) -> np.ndarray:
    base = parameter_grid(start, stop, PRESENTATION_MU_STEP)
    return refine_parameter_grid(
        base_grid=base,
        local_windows=local_windows,
        local_step=LOCAL_MU_REFINEMENT_STEP if local_windows else None,
        anchor_points=anchor_points,
    )


def analysis_beta_grid(
    start: float = 0.0,
    stop: float = 15.0,
    local_windows: tuple[tuple[float, float], ...] = (),
    anchor_points: tuple[float, ...] = (),
) -> np.ndarray:
    base = parameter_grid(start, stop, ANALYSIS_BETA_STEP)
    return refine_parameter_grid(
        base_grid=base,
        local_windows=local_windows,
        local_step=LOCAL_BETA_REFINEMENT_STEP if local_windows else None,
        anchor_points=anchor_points,
    )


def analysis_mu_grid(
    start: float = 0.0,
    stop: float = 0.9,
    local_windows: tuple[tuple[float, float], ...] = (),
    anchor_points: tuple[float, ...] = (),
) -> np.ndarray:
    base = parameter_grid(start, stop, ANALYSIS_MU_STEP)
    return refine_parameter_grid(
        base_grid=base,
        local_windows=local_windows,
        local_step=LOCAL_MU_REFINEMENT_STEP if local_windows else None,
        anchor_points=anchor_points,
    )


def nominal_step(values: np.ndarray) -> float:
    diffs = np.diff(np.asarray(values, dtype=float))
    positive = diffs[diffs > 1e-12]
    if len(positive) == 0:
        return 0.0
    return float(np.min(positive))
