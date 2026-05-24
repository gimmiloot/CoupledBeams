from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np

from my_project.analytic.FreqMuNet import roots_clamped_supported, single_lambda
from my_project.analytic.formulas_thickness_mismatch import (
    find_first_n_roots_eta,
    thickness_mismatch_factors,
    thickness_to_length_ratios,
    thin_rod_validity,
)
from scripts.lib.thickness_mismatch_mac_tracking import (
    ShapeMacTrackingResult,
    track_mu_branches_shape_mac,
)


VALID_TOL = 1e-12


@dataclass(frozen=True)
class ThicknessViolationSegment:
    start_mu: float
    end_mu: float
    rods: tuple[int, ...]


@dataclass(frozen=True)
class ThicknessRatioSummary:
    eta: float
    segments: tuple[ThicknessViolationSegment, ...]
    max_ratio: float
    max_ratio_rod: int
    max_ratio_mu: float
    max_ratio_1: float
    max_ratio_2: float

    @property
    def has_violations(self) -> bool:
        return bool(self.segments)


def number_token(value: float) -> str:
    return f"{float(value):.10g}".replace("-", "m").replace(".", "p")


def mu_grid(mu_min: float, mu_max: float, mu_step: float) -> np.ndarray:
    values = np.arange(float(mu_min), float(mu_max) + 0.5 * float(mu_step), float(mu_step), dtype=float)
    if values.size == 0 or abs(float(values[-1]) - float(mu_max)) > 1e-10:
        values = np.append(values, float(mu_max))
    values[0] = float(mu_min)
    values[-1] = float(mu_max)
    return np.unique(np.round(values, 12))


def roots_by_mu_eta(
    *,
    beta_rad: float,
    epsilon: float,
    eta: float,
    mu_values: Sequence[float],
    n_roots: int,
    root_lmax0: float,
    root_scan_step: float,
) -> dict[float, np.ndarray]:
    roots: dict[float, np.ndarray] = {}
    for mu in np.asarray(mu_values, dtype=float):
        values = find_first_n_roots_eta(
            float(beta_rad),
            float(mu),
            float(epsilon),
            float(eta),
            int(n_roots),
            Lmax0=float(root_lmax0),
            scan_step=float(root_scan_step),
        )
        if np.any(~np.isfinite(values)):
            raise RuntimeError(f"Missing roots for mu={float(mu):g}, eta={float(eta):g}.")
        roots[float(mu)] = values
    return roots


def track_descendants_from_mu0(
    *,
    beta_rad: float,
    epsilon: float,
    eta: float,
    mu_values: Sequence[float],
    roots_by_mu: Mapping[float, np.ndarray],
    num_descendants: int,
    num_shape_samples: int,
    mac_warning_threshold: float,
    mac_margin_warning_threshold: float,
    max_sorted_position_jump: int,
) -> ShapeMacTrackingResult:
    mu_values_array = np.asarray(mu_values, dtype=float)
    if not np.isclose(float(mu_values_array[0]), 0.0, rtol=0.0, atol=1e-12):
        raise ValueError("Fixed-eta Lambda(mu) descendants must be seeded at mu=0.")
    return track_mu_branches_shape_mac(
        beta_rad=float(beta_rad),
        epsilon=float(epsilon),
        eta=float(eta),
        mu_values=mu_values_array,
        roots_by_mu=roots_by_mu,
        num_tracked_branches=int(num_descendants),
        num_shape_samples=int(num_shape_samples),
        mac_warning_threshold=float(mac_warning_threshold),
        mac_margin_warning_threshold=float(mac_margin_warning_threshold),
        max_sorted_position_jump=int(max_sorted_position_jump),
        max_mu_step_for_confidence=float(mu_values_array[1] - mu_values_array[0])
        if len(mu_values_array) > 1
        else None,
    )


def thickness_ratio_values(
    *,
    epsilon: float,
    eta: float,
    mu_values: Sequence[float],
) -> np.ndarray:
    return np.array(
        [thickness_to_length_ratios(float(epsilon), float(mu), float(eta)) for mu in mu_values],
        dtype=float,
    )


def valid_mask(
    *,
    epsilon: float,
    eta: float,
    mu_values: Sequence[float],
    limit: float,
) -> np.ndarray:
    return np.array(
        [thin_rod_validity(float(epsilon), float(mu), float(eta), limit=float(limit)) for mu in mu_values],
        dtype=bool,
    )


def violating_rods(
    *,
    epsilon: float,
    eta: float,
    mu: float,
    limit: float,
) -> tuple[int, ...]:
    ratio1, ratio2 = thickness_to_length_ratios(float(epsilon), float(mu), float(eta))
    rods: list[int] = []
    if ratio1 > float(limit) + VALID_TOL:
        rods.append(1)
    if ratio2 > float(limit) + VALID_TOL:
        rods.append(2)
    return tuple(rods)


def thickness_violation_segments(
    *,
    epsilon: float,
    eta: float,
    mu_values: Sequence[float],
    limit: float,
) -> tuple[ThicknessViolationSegment, ...]:
    mu_array = np.asarray(mu_values, dtype=float)
    segments: list[ThicknessViolationSegment] = []
    start_idx: int | None = None
    current_rods: tuple[int, ...] = ()

    for idx, mu in enumerate(mu_array):
        rods = violating_rods(epsilon=float(epsilon), eta=float(eta), mu=float(mu), limit=float(limit))
        if rods:
            if start_idx is None:
                start_idx = idx
                current_rods = rods
            elif rods != current_rods:
                segments.append(
                    ThicknessViolationSegment(
                        start_mu=float(mu_array[start_idx]),
                        end_mu=float(mu_array[idx - 1]),
                        rods=current_rods,
                    )
                )
                start_idx = idx
                current_rods = rods
        elif start_idx is not None:
            segments.append(
                ThicknessViolationSegment(
                    start_mu=float(mu_array[start_idx]),
                    end_mu=float(mu_array[idx - 1]),
                    rods=current_rods,
                )
            )
            start_idx = None
            current_rods = ()

    if start_idx is not None:
        segments.append(
            ThicknessViolationSegment(
                start_mu=float(mu_array[start_idx]),
                end_mu=float(mu_array[-1]),
                rods=current_rods,
            )
        )
    return tuple(segments)


def thickness_ratio_summary(
    *,
    epsilon: float,
    eta: float,
    mu_values: Sequence[float],
    limit: float,
) -> ThicknessRatioSummary:
    mu_array = np.asarray(mu_values, dtype=float)
    ratios = thickness_ratio_values(epsilon=float(epsilon), eta=float(eta), mu_values=mu_array)
    flat_idx = int(np.argmax(ratios))
    mu_idx, rod_idx = np.unravel_index(flat_idx, ratios.shape)
    return ThicknessRatioSummary(
        eta=float(eta),
        segments=thickness_violation_segments(
            epsilon=float(epsilon),
            eta=float(eta),
            mu_values=mu_array,
            limit=float(limit),
        ),
        max_ratio=float(ratios[mu_idx, rod_idx]),
        max_ratio_rod=int(rod_idx) + 1,
        max_ratio_mu=float(mu_array[mu_idx]),
        max_ratio_1=float(np.max(ratios[:, 0])),
        max_ratio_2=float(np.max(ratios[:, 1])),
    )


def rods_label(rods: tuple[int, ...]) -> str:
    return "both" if rods == (1, 2) else ", ".join(str(rod) for rod in rods)


def contiguous_true_runs(mask: Sequence[bool]) -> list[tuple[int, int]]:
    mask_array = np.asarray(mask, dtype=bool)
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for idx, value in enumerate(mask_array):
        if bool(value) and start is None:
            start = idx
        elif not bool(value) and start is not None:
            runs.append((start, idx - 1))
            start = None
    if start is not None:
        runs.append((start, len(mask_array) - 1))
    return runs


def plot_with_validity_split(
    ax: plt.Axes,
    x_values: Sequence[float],
    y_values: Sequence[float],
    valid: Sequence[bool],
    *,
    color: str,
    linewidth: float,
    label: str | None = None,
    alpha: float = 1.0,
    zorder: float | None = None,
) -> None:
    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    valid_array = np.asarray(valid, dtype=bool)

    valid_runs = contiguous_true_runs(valid_array)
    for run_idx, (start, end) in enumerate(valid_runs):
        local_label = label if run_idx == 0 else None
        if end > start:
            ax.plot(
                x[start : end + 1],
                y[start : end + 1],
                color=color,
                lw=float(linewidth),
                ls="-",
                alpha=float(alpha),
                label=local_label,
                zorder=zorder,
            )

    for start, end in contiguous_true_runs(~valid_array):
        plot_start = max(0, start - 1)
        plot_end = min(len(x) - 1, end + 1)
        if plot_end > plot_start:
            ax.plot(
                x[plot_start : plot_end + 1],
                y[plot_start : plot_end + 1],
                color=color,
                lw=float(linewidth),
                ls="--",
                alpha=float(alpha),
                label=None,
                zorder=zorder,
            )


def tracking_warning_rows(rows: Sequence[dict[str, float | int | str]]) -> list[dict[str, float | int | str]]:
    return [
        row
        for row in rows
        if str(row.get("requires_refined_check", "no")) == "yes"
        or str(row.get("frequency_mac_disagreement", "no")) == "yes"
        or str(row.get("suspicious_assignment", "no")) == "yes"
    ]


def tracking_warning_summary(rows: Sequence[dict[str, float | int | str]]) -> dict[str, int]:
    return {
        "low_mac": sum(1 for row in rows if str(row.get("low_mac", "no")) == "yes"),
        "low_margin": sum(1 for row in rows if str(row.get("low_margin", "no")) == "yes"),
        "unresolved": sum(1 for row in rows if str(row.get("unresolved_assignment", "no")) == "yes"),
        "suspicious": sum(1 for row in rows if str(row.get("suspicious_assignment", "no")) == "yes"),
        "requires_refined": sum(1 for row in rows if str(row.get("requires_refined_check", "no")) == "yes"),
        "frequency_mac_disagreement": sum(
            1 for row in rows if str(row.get("frequency_mac_disagreement", "no")) == "yes"
        ),
    }


def isolated_clamped_supported_references(
    *,
    epsilon: float,
    eta: float,
    mu_values: Sequence[float],
    n_curves: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return isolated-rod CS Lambda references for rods 1 and 2.

    The convention matches the existing Lambda(mu) dashed references:
    clamped-supported single rods with roots tan(alpha)=tanh(alpha).  The
    thickness-mismatch common Lambda scaling adds sqrt(tau_i) to the
    equal-radius references, so Lambda_i = alpha*sqrt(tau_i)/(1 +/- mu).
    """
    del epsilon
    mu_array = np.asarray(mu_values, dtype=float)
    alphas = roots_clamped_supported(int(n_curves))
    rod1 = np.full((int(n_curves), len(mu_array)), np.nan, dtype=float)
    rod2 = np.full_like(rod1, np.nan)
    for col, mu in enumerate(mu_array):
        factors = thickness_mismatch_factors(float(mu), float(eta))
        rod1[:, col] = alphas * np.sqrt(factors.tau1) / (1.0 - float(mu))
        rod2[:, col] = alphas * np.sqrt(factors.tau2) / (1.0 + float(mu))
    if abs(float(eta)) < 1e-15:
        equal_rod1 = single_lambda(alphas, 1.0, 1.0 - mu_array)
        equal_rod2 = single_lambda(alphas, 1.0, 1.0 + mu_array)
        if not (np.allclose(rod1, equal_rod1) and np.allclose(rod2, equal_rod2)):
            raise RuntimeError("Eta=0 isolated references must reduce to alpha/(1 +/- mu).")
    return alphas, rod1, rod2


__all__ = [
    "ThicknessRatioSummary",
    "ThicknessViolationSegment",
    "contiguous_true_runs",
    "isolated_clamped_supported_references",
    "mu_grid",
    "number_token",
    "plot_with_validity_split",
    "rods_label",
    "roots_by_mu_eta",
    "thickness_ratio_summary",
    "thickness_ratio_values",
    "thickness_violation_segments",
    "track_descendants_from_mu0",
    "tracking_warning_rows",
    "tracking_warning_summary",
    "valid_mask",
    "violating_rods",
]
