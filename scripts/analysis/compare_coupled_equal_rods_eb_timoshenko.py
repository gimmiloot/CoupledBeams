from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Callable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq


# =========================
# User-editable parameters
# =========================
L_SEGMENT = 1.0
E = 1.0
RHO = 1.0
NU = 0.3
EPSILON_VALUES = [0.0025, 0.005, 0.01, 0.015, 0.02, 0.025, 0.0375, 0.05]
N_MODES = 6
THICKNESS_RATIO_LIMIT = 0.1
BETA_DEG = 15.0
MU = 0.0


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import assemble_clamped_coupled_matrix  # noqa: E402
from my_project.analytic.solvers import find_first_n_roots  # noqa: E402


OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_CSV = OUTPUT_DIR / "coupled_equal_rods_beta15_eb_timoshenko_frequencies.csv"
OUTPUT_PNG = OUTPUT_DIR / "coupled_equal_rods_beta15_eb_timoshenko_lambda_vs_epsilon.png"
OUTPUT_REPORT = OUTPUT_DIR / "coupled_equal_rods_beta15_eb_timoshenko_report.md"
OUTPUT_SENSITIVITY_CSV = OUTPUT_DIR / "coupled_equal_rods_beta15_timoshenko_kappa_sensitivity.csv"

ROOT_SCAN_START = 0.2
ROOT_SCAN_STEP = 0.002
ROOT_DEDUP_TOL = 1.0e-7
SINGULAR_RESIDUAL_TOL = 1.0e-7
CUTOFF_WARNING_RATIO = 0.8
KAPPA_SENSITIVITY_EPSILONS = [0.025, 0.05]


@dataclass(frozen=True)
class Section:
    radius: float
    area: float
    inertia: float
    shear_modulus: float
    kappa: float

    @property
    def shear_stiffness(self) -> float:
        return self.kappa * self.shear_modulus * self.area

    @property
    def bending_stiffness(self) -> float:
        return E * self.inertia

    @property
    def mass_per_length(self) -> float:
        return RHO * self.area

    @property
    def rotary_inertia_per_length(self) -> float:
        return RHO * self.inertia


@dataclass(frozen=True)
class TimoshenkoBasis:
    a: float
    b: float
    h: float
    q: float
    warnings: tuple[str, ...]


def circular_shear_coefficient(nu: float) -> float:
    """Circular-rod shear coefficient used by Diaz-de-Anda et al."""
    return (6.0 + 12.0 * nu + 6.0 * nu**2) / (7.0 + 12.0 * nu + 4.0 * nu**2)


def section_from_epsilon(epsilon: float, kappa: float | None = None) -> Section:
    radius = 2.0 * float(epsilon) * L_SEGMENT
    area = np.pi * radius**2
    inertia = np.pi * radius**4 / 4.0
    shear_modulus = E / (2.0 * (1.0 + NU))
    return Section(
        radius=radius,
        area=area,
        inertia=inertia,
        shear_modulus=shear_modulus,
        kappa=circular_shear_coefficient(NU) if kappa is None else float(kappa),
    )


def project_omega(Lambda: float, epsilon: float) -> float:
    return float(epsilon) * float(Lambda) ** 2


def omega_cutoff(section: Section) -> float:
    return float(np.sqrt(section.shear_stiffness / (RHO * section.inertia)))


def lambda_cutoff(epsilon: float, section: Section) -> float:
    return float(np.sqrt(omega_cutoff(section) / float(epsilon)))


def diameter_to_segment_length(epsilon: float) -> float:
    return 4.0 * float(epsilon)


def thin_rod_valid(epsilon: float) -> bool:
    return diameter_to_segment_length(epsilon) <= THICKNESS_RATIO_LIMIT + 1.0e-15


def row_normalized(matrix: np.ndarray) -> np.ndarray:
    out = np.array(matrix, dtype=float, copy=True)
    norms = np.linalg.norm(out, axis=1)
    for index, norm in enumerate(norms):
        if norm > 0.0 and np.isfinite(norm):
            out[index, :] /= norm
    return out


def normalized_det(matrix: np.ndarray) -> float:
    return float(np.linalg.det(row_normalized(matrix)))


def singular_residual(matrix: np.ndarray) -> float:
    scaled = row_normalized(matrix)
    singular_values = np.linalg.svd(scaled, compute_uv=False)
    largest = float(np.max(singular_values))
    if largest == 0.0 or not np.isfinite(largest):
        return float("nan")
    return float(np.min(singular_values) / largest)


def timo_basis(Lambda: float, epsilon: float, section: Section) -> TimoshenkoBasis:
    omega = project_omega(float(Lambda), float(epsilon))
    k_stiff = section.shear_stiffness
    b_stiff = section.bending_stiffness
    m = section.mass_per_length
    j = section.rotary_inertia_per_length

    c2 = omega**2 * (b_stiff * m / k_stiff + j)
    c0 = j * m * omega**4 / k_stiff - m * omega**2
    roots = np.roots([b_stiff, c2, c0])
    warnings: list[str] = []

    scale = max(1.0, float(np.max(np.abs(roots))))
    if np.max(np.abs(np.imag(roots))) > 1.0e-8 * scale:
        warnings.append(
            "Timoshenko lambda^2 roots have non-negligible imaginary parts; "
            f"Lambda={float(Lambda):.8g}, roots={roots}."
        )
    z_values = np.real(roots)
    positives = [float(z) for z in z_values if z > 1.0e-12]
    negatives = [float(z) for z in z_values if z < -1.0e-12]

    if len(positives) != 1 or len(negatives) != 1:
        warnings.append(
            "Expected one positive and one negative lambda^2 root below cut-off; "
            f"Lambda={float(Lambda):.8g}, roots={z_values}."
        )
        z_positive = float(np.max(z_values))
        z_negative = float(np.min(z_values))
        if z_positive <= 0.0 or z_negative >= 0.0:
            raise ValueError(warnings[-1])
    else:
        z_positive = positives[0]
        z_negative = negatives[0]

    a = float(np.sqrt(z_positive))
    b = float(np.sqrt(-z_negative))
    if a <= 0.0 or b <= 0.0:
        raise ValueError(f"Invalid Timoshenko basis parameters a={a:g}, b={b:g}.")

    h = a + m * omega**2 / (k_stiff * a)
    q = b - m * omega**2 / (k_stiff * b)
    if abs(q) <= 1.0e-12 * max(1.0, abs(b)):
        warnings.append(
            "Timoshenko trigonometric coupling q is close to zero; "
            f"Lambda={float(Lambda):.8g}, q={q:.8g}."
        )
    return TimoshenkoBasis(a=a, b=b, h=float(h), q=float(q), warnings=tuple(warnings))


def bending_endpoint_columns(x: float, basis: TimoshenkoBasis) -> dict[str, np.ndarray]:
    a = basis.a
    b = basis.b
    h = basis.h
    q = basis.q
    ax = a * float(x)
    bx = b * float(x)

    cosh_ax = np.cosh(ax)
    sinh_ax = np.sinh(ax)
    cos_bx = np.cos(bx)
    sin_bx = np.sin(bx)

    w = np.array(
        [
            cosh_ax - cos_bx,
            sinh_ax - (h / q) * sin_bx,
        ],
        dtype=float,
    )
    psi = np.array(
        [
            h * sinh_ax + q * sin_bx,
            h * (cosh_ax - cos_bx),
        ],
        dtype=float,
    )
    w_prime = np.array(
        [
            a * sinh_ax + b * sin_bx,
            a * cosh_ax - (h / q) * b * cos_bx,
        ],
        dtype=float,
    )
    psi_prime = np.array(
        [
            h * a * cosh_ax + q * b * cos_bx,
            h * (a * sinh_ax + b * sin_bx),
        ],
        dtype=float,
    )
    return {"w": w, "psi": psi, "w_prime": w_prime, "psi_prime": psi_prime}


def endpoint_columns(
    x: float,
    theta: float,
    basis: TimoshenkoBasis,
    section: Section,
) -> dict[str, np.ndarray]:
    bend = bending_endpoint_columns(float(x), basis)
    w = np.zeros(3, dtype=float)
    psi = np.zeros(3, dtype=float)
    w_prime = np.zeros(3, dtype=float)
    psi_prime = np.zeros(3, dtype=float)
    u = np.zeros(3, dtype=float)
    u_prime = np.zeros(3, dtype=float)

    w[:2] = bend["w"]
    psi[:2] = bend["psi"]
    w_prime[:2] = bend["w_prime"]
    psi_prime[:2] = bend["psi_prime"]
    u[2] = np.sin(theta * float(x))
    u_prime[2] = theta * np.cos(theta * float(x))

    q_force = section.shear_stiffness * (w_prime - psi)
    moment = section.bending_stiffness * psi_prime
    normal = E * section.area * u_prime
    return {"w": w, "u": u, "psi": psi, "Q": q_force, "M": moment, "N": normal}


def timo_coupling_matrix(
    Lambda: float,
    epsilon: float,
    *,
    kappa: float | None = None,
) -> tuple[np.ndarray, tuple[str, ...]]:
    section = section_from_epsilon(float(epsilon), kappa=kappa)
    basis = timo_basis(float(Lambda), float(epsilon), section)
    theta = project_omega(float(Lambda), float(epsilon))
    rod1 = endpoint_columns(1.0, theta, basis, section)
    rod2 = endpoint_columns(-1.0, theta, basis, section)
    cb = float(np.cos(np.deg2rad(BETA_DEG)))
    sb = float(np.sin(np.deg2rad(BETA_DEG)))

    matrix = np.zeros((6, 6), dtype=float)

    # Parity-adjusted displacement block matching assemble_clamped_coupled_matrix.
    # The local theory writes the joint with x2 = -1; the implemented EB
    # determinant has already absorbed odd-function signs in the rod-2 columns.
    # This block is the convention that recovers the EB determinant as
    # epsilon -> 0.
    # w1(1) = w2(-1)*cos(beta) - u2(-1)*sin(beta)
    matrix[0, 0:3] = rod1["w"]
    matrix[0, 3:6] = -cb * rod2["w"] + sb * rod2["u"]

    # u1(1) = w2(-1)*sin(beta) + u2(-1)*cos(beta)
    matrix[1, 0:3] = rod1["u"]
    matrix[1, 3:6] = -sb * rod2["w"] - cb * rod2["u"]

    # Timoshenko rotation continuity replaces the Euler-Bernoulli slope condition.
    matrix[2, 0:3] = rod1["psi"]
    matrix[2, 3:6] = -rod2["psi"]

    # M1(1) = M2(-1)
    matrix[3, 0:3] = rod1["M"]
    matrix[3, 3:6] = -rod2["M"]

    # -Q1(1) = -Q2(-1)*cos(beta) + N2(-1)*sin(beta)
    matrix[4, 0:3] = -rod1["Q"]
    matrix[4, 3:6] = cb * rod2["Q"] - sb * rod2["N"]

    # N1(1) = Q2(-1)*sin(beta) + N2(-1)*cos(beta)
    matrix[5, 0:3] = rod1["N"]
    matrix[5, 3:6] = -sb * rod2["Q"] - cb * rod2["N"]

    return matrix, basis.warnings


def timo_det_for_scan(Lambda: float, epsilon: float, *, kappa: float | None = None) -> float:
    try:
        matrix, _ = timo_coupling_matrix(float(Lambda), float(epsilon), kappa=kappa)
    except (FloatingPointError, ValueError, OverflowError):
        return float("nan")
    if not np.all(np.isfinite(matrix)):
        return float("nan")
    return normalized_det(matrix)


def find_roots_by_sign_scan(
    func: Callable[[float], float],
    n_roots: int,
    *,
    start: float,
    upper: float,
    scan_step: float,
    grow_factor: float = 1.35,
    max_tries: int = 8,
) -> tuple[np.ndarray, list[str]]:
    warnings: list[str] = []
    roots: list[float] = []
    current_upper = float(upper)

    for _ in range(max_tries):
        roots.clear()
        grid = np.arange(float(start), current_upper + scan_step, scan_step)
        values = np.array([func(float(x)) for x in grid], dtype=float)

        for left, right, f_left, f_right in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
            if not (np.isfinite(f_left) and np.isfinite(f_right)):
                continue
            candidate: float | None = None
            if f_left == 0.0:
                candidate = float(left)
            elif f_left * f_right < 0.0:
                try:
                    candidate = float(
                        brentq(
                            func,
                            float(left),
                            float(right),
                            xtol=1.0e-12,
                            rtol=1.0e-12,
                            maxiter=120,
                        )
                    )
                except ValueError as exc:
                    warnings.append(
                        f"Skipped bracket [{float(left):.8g}, {float(right):.8g}]: {exc}"
                    )
                    continue

            if candidate is None:
                continue
            if roots and abs(candidate - roots[-1]) <= ROOT_DEDUP_TOL:
                continue
            roots.append(candidate)
            if len(roots) >= n_roots:
                return np.array(roots[:n_roots], dtype=float), warnings

        current_upper *= grow_factor

    warnings.append(
        f"Found only {len(roots)} roots below Lambda={current_upper:.8g}; "
        "increase scan range or reduce scan step."
    )
    out = np.full(int(n_roots), np.nan, dtype=float)
    out[: len(roots)] = roots[:n_roots]
    return out, warnings


def eb_roots(epsilon: float) -> np.ndarray:
    return find_first_n_roots(
        beta=float(np.deg2rad(BETA_DEG)),
        mu=MU,
        eps=float(epsilon),
        n_roots=N_MODES,
        Lmin=ROOT_SCAN_START,
        Lmax0=35.0,
        scan_step=0.005,
        grow_factor=1.35,
        max_tries=8,
    )


def timo_roots(
    epsilon: float,
    eb_values: np.ndarray,
    *,
    kappa: float | None = None,
) -> tuple[np.ndarray, list[str]]:
    upper = max(float(np.nanmax(eb_values)) * 1.25, 12.0)
    return find_roots_by_sign_scan(
        lambda value: timo_det_for_scan(value, float(epsilon), kappa=kappa),
        N_MODES,
        start=ROOT_SCAN_START,
        upper=upper,
        scan_step=ROOT_SCAN_STEP,
    )


def eb_residual(Lambda: float, epsilon: float) -> float:
    matrix = assemble_clamped_coupled_matrix(
        float(Lambda),
        float(np.deg2rad(BETA_DEG)),
        MU,
        float(epsilon),
    )
    return singular_residual(matrix)


def timo_residual_and_warnings(
    Lambda: float,
    epsilon: float,
    *,
    kappa: float | None = None,
) -> tuple[float, tuple[str, ...]]:
    if not np.isfinite(Lambda):
        return float("nan"), ("Timoshenko root was not found.",)
    matrix, warnings = timo_coupling_matrix(float(Lambda), float(epsilon), kappa=kappa)
    return singular_residual(matrix), warnings


def format_float(value: float) -> str:
    if not np.isfinite(value):
        return "nan"
    return f"{float(value):.12g}"


def compute_rows() -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    warnings: list[str] = []

    for epsilon in EPSILON_VALUES:
        eb_values = eb_roots(float(epsilon))
        timo_values, root_warnings = timo_roots(float(epsilon), eb_values)
        warnings.extend(f"epsilon={float(epsilon):g}: {warning}" for warning in root_warnings)
        section = section_from_epsilon(float(epsilon))
        omega_c = omega_cutoff(section)
        lambda_c = lambda_cutoff(float(epsilon), section)
        d_over_l = diameter_to_segment_length(float(epsilon))
        valid = thin_rod_valid(float(epsilon))

        for mode_index in range(N_MODES):
            lambda_eb = float(eb_values[mode_index])
            lambda_timo = float(timo_values[mode_index])
            omega_timo = project_omega(lambda_timo, float(epsilon)) if np.isfinite(lambda_timo) else float("nan")
            omega_ratio = omega_timo / omega_c if np.isfinite(omega_timo) else float("nan")
            lambda_ratio = lambda_timo / lambda_c if np.isfinite(lambda_timo) else float("nan")
            below_cutoff = bool(np.isfinite(omega_ratio) and omega_ratio < 1.0)
            residual_eb = eb_residual(lambda_eb, float(epsilon)) if np.isfinite(lambda_eb) else float("nan")
            residual_timo, basis_warnings = timo_residual_and_warnings(lambda_timo, float(epsilon))
            warnings.extend(f"epsilon={float(epsilon):g}, mode={mode_index + 1}: {item}" for item in basis_warnings)
            abs_diff = lambda_timo - lambda_eb if np.isfinite(lambda_timo) and np.isfinite(lambda_eb) else float("nan")
            rel_diff = abs_diff / lambda_eb if np.isfinite(abs_diff) and lambda_eb != 0.0 else float("nan")

            notes: list[str] = []
            if not valid:
                notes.append("outside_thin_rod_applicability")
            if np.isfinite(residual_timo) and residual_timo > SINGULAR_RESIDUAL_TOL:
                notes.append("large_timoshenko_residual")
            if not np.isfinite(lambda_timo):
                notes.append("missing_timoshenko_root")
            if basis_warnings:
                notes.append("basis_warning")
            if np.isfinite(omega_ratio) and omega_ratio >= CUTOFF_WARNING_RATIO:
                notes.append("near_cutoff")
            if not below_cutoff:
                notes.append("not_below_cutoff")

            rows.append(
                {
                    "epsilon": float(epsilon),
                    "diameter_to_segment_length": d_over_l,
                    "mode": mode_index + 1,
                    "Lambda_EB": lambda_eb,
                    "Lambda_Timoshenko": lambda_timo,
                    "Omega_Timoshenko": omega_timo,
                    "Omega_cutoff": omega_c,
                    "Omega_over_cutoff": omega_ratio,
                    "Lambda_cutoff": lambda_c,
                    "Lambda_over_cutoff": lambda_ratio,
                    "below_cutoff": below_cutoff,
                    "abs_diff": abs_diff,
                    "rel_diff": rel_diff,
                    "thin_rod_valid": valid,
                    "det_residual_EB": residual_eb,
                    "det_residual_Timoshenko": residual_timo,
                    "notes": ";".join(notes),
                }
            )

    return rows, sorted(set(warnings))


def write_csv(rows: Sequence[dict[str, object]]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "epsilon",
        "diameter_to_segment_length",
        "mode",
        "Lambda_EB",
        "Lambda_Timoshenko",
        "Omega_Timoshenko",
        "Omega_cutoff",
        "Omega_over_cutoff",
        "Lambda_cutoff",
        "Lambda_over_cutoff",
        "below_cutoff",
        "abs_diff",
        "rel_diff",
        "thin_rod_valid",
        "det_residual_EB",
        "det_residual_Timoshenko",
        "notes",
    ]
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def kappa_sensitivity_options() -> list[tuple[str, float]]:
    return [
        ("baseline_circular", circular_shear_coefficient(NU)),
        ("kappa_1p0", 1.0),
    ]


def compute_kappa_sensitivity_rows() -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    warnings: list[str] = []

    for epsilon in KAPPA_SENSITIVITY_EPSILONS:
        eb_values = eb_roots(float(epsilon))
        baseline_kappa = circular_shear_coefficient(NU)
        baseline_roots, baseline_warnings = timo_roots(
            float(epsilon),
            eb_values,
            kappa=baseline_kappa,
        )
        warnings.extend(
            f"sensitivity epsilon={float(epsilon):g}, baseline_circular: {warning}"
            for warning in baseline_warnings
        )

        for kappa_label, kappa_value in kappa_sensitivity_options():
            if kappa_label == "baseline_circular":
                roots = baseline_roots
                root_warnings: list[str] = []
            else:
                roots, root_warnings = timo_roots(float(epsilon), eb_values, kappa=kappa_value)
            warnings.extend(
                f"sensitivity epsilon={float(epsilon):g}, {kappa_label}: {warning}"
                for warning in root_warnings
            )

            for mode_index in range(N_MODES):
                lambda_eb = float(eb_values[mode_index])
                lambda_timo = float(roots[mode_index])
                lambda_baseline = float(baseline_roots[mode_index])
                delta = (
                    lambda_timo - lambda_baseline
                    if np.isfinite(lambda_timo) and np.isfinite(lambda_baseline)
                    else float("nan")
                )
                rel_delta = delta / lambda_baseline if np.isfinite(delta) and lambda_baseline != 0.0 else float("nan")
                baseline_diff = (
                    lambda_baseline - lambda_eb
                    if np.isfinite(lambda_baseline) and np.isfinite(lambda_eb)
                    else float("nan")
                )
                ratio_to_eb_timo = (
                    abs(delta) / abs(baseline_diff)
                    if np.isfinite(delta) and np.isfinite(baseline_diff) and baseline_diff != 0.0
                    else float("nan")
                )

                rows.append(
                    {
                        "epsilon": float(epsilon),
                        "mode": mode_index + 1,
                        "kappa_label": kappa_label,
                        "kappa": float(kappa_value),
                        "Lambda_EB": lambda_eb,
                        "Lambda_Timoshenko": lambda_timo,
                        "Lambda_baseline_kappa": lambda_baseline,
                        "delta_Lambda_vs_baseline": delta,
                        "rel_delta_Lambda_vs_baseline": rel_delta,
                        "baseline_Timoshenko_minus_EB": baseline_diff,
                        "abs_sensitivity_over_abs_EB_Timoshenko_diff": ratio_to_eb_timo,
                    }
                )

    return rows, sorted(set(warnings))


def write_kappa_sensitivity_csv(rows: Sequence[dict[str, object]]) -> None:
    fieldnames = [
        "epsilon",
        "mode",
        "kappa_label",
        "kappa",
        "Lambda_EB",
        "Lambda_Timoshenko",
        "Lambda_baseline_kappa",
        "delta_Lambda_vs_baseline",
        "rel_delta_Lambda_vs_baseline",
        "baseline_Timoshenko_minus_EB",
        "abs_sensitivity_over_abs_EB_Timoshenko_diff",
    ]
    with OUTPUT_SENSITIVITY_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_rows(rows: Sequence[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    mode_colors = plt.cm.viridis(np.linspace(0.08, 0.92, N_MODES))

    x_values = np.array([float(row["epsilon"]) for row in rows], dtype=float)
    x_min = float(np.nanmin(x_values))
    x_max = float(np.nanmax(x_values))
    validity_boundary = THICKNESS_RATIO_LIMIT / 4.0
    if x_max > validity_boundary:
        ax.axvspan(validity_boundary, x_max, color="0.9", alpha=0.75, label="outside 4*epsilon <= 0.1")
    ax.axvline(validity_boundary, color="0.45", linestyle=":", linewidth=1.1)

    for mode in range(1, N_MODES + 1):
        mode_rows = [row for row in rows if int(row["mode"]) == mode]
        x = np.array([float(row["epsilon"]) for row in mode_rows], dtype=float)
        y_eb = np.array([float(row["Lambda_EB"]) for row in mode_rows], dtype=float)
        y_timo = np.array([float(row["Lambda_Timoshenko"]) for row in mode_rows], dtype=float)
        valid = np.array([bool(row["thin_rod_valid"]) for row in mode_rows], dtype=bool)
        color = mode_colors[mode - 1]

        ax.plot(
            x,
            y_eb,
            color=color,
            linestyle="--",
            linewidth=1.15,
            alpha=0.55,
            label=f"EB sorted mode {mode}" if mode == 1 else None,
        )
        ax.plot(x[valid], y_timo[valid], marker="o", color=color, linewidth=1.8, label=f"Timo mode {mode}")
        if np.any(~valid):
            ax.plot(
                x[~valid],
                y_timo[~valid],
                marker="o",
                markerfacecolor="white",
                color=color,
                linestyle=":",
                linewidth=1.8,
            )

    ax.set_xlabel("epsilon = r/(2 L_segment)")
    ax.set_ylabel("Lambda")
    ax.set_title("Coupled equal rods, beta=15 deg: sorted EB vs Timoshenko frequencies")
    ax.grid(True, alpha=0.3)
    ax.legend(ncols=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=180, bbox_inches="tight")
    plt.close(fig)


def max_diff_rows(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for epsilon in EPSILON_VALUES:
        subset = [row for row in rows if abs(float(row["epsilon"]) - float(epsilon)) <= 1.0e-15]
        finite = [row for row in subset if np.isfinite(float(row["rel_diff"]))]
        if not finite:
            continue
        selected = max(finite, key=lambda row: abs(float(row["rel_diff"])))
        out.append(selected)
    return out


def markdown_numeric_summary(rows: Sequence[dict[str, object]]) -> str:
    lines = [
        "| epsilon | 4*epsilon | valid | max abs rel diff | min cut-off margin | mode | Lambda_EB | Lambda_Timoshenko |",
        "| ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in max_diff_rows(rows):
        subset = [item for item in rows if abs(float(item["epsilon"]) - float(row["epsilon"])) <= 1.0e-15]
        min_margin = min(1.0 - float(item["Omega_over_cutoff"]) for item in subset)
        lines.append(
            "| "
            f"{float(row['epsilon']):.5g} | "
            f"{float(row['diameter_to_segment_length']):.5g} | "
            f"{bool(row['thin_rod_valid'])} | "
            f"{abs(float(row['rel_diff'])):.6g} | "
            f"{min_margin:.6g} | "
            f"{int(row['mode'])} | "
            f"{float(row['Lambda_EB']):.8g} | "
            f"{float(row['Lambda_Timoshenko']):.8g} |"
        )
    return "\n".join(lines)


def cutoff_report_lines(rows: Sequence[dict[str, object]]) -> list[str]:
    finite_rows = [
        row
        for row in rows
        if np.isfinite(float(row["Omega_over_cutoff"])) and np.isfinite(float(row["Omega_cutoff"]))
    ]
    if not finite_rows:
        return ["## Cut-Off Frequency Check", "", "No finite cut-off ratios were computed."]

    closest = max(finite_rows, key=lambda row: float(row["Omega_over_cutoff"]))
    max_ratio = float(closest["Omega_over_cutoff"])
    min_margin = 1.0 - max_ratio
    near_rows = [row for row in finite_rows if float(row["Omega_over_cutoff"]) >= CUTOFF_WARNING_RATIO]
    at_or_above_rows = [row for row in finite_rows if not bool(row["below_cutoff"])]

    lines = [
        "## Cut-Off Frequency Check",
        "",
        "For each Timoshenko root this diagnostic computes",
        "",
        "```text",
        "Omega_c = sqrt(kappa*G*A/(rho*I))",
        "```",
        "",
        "With the coupled equal-rod normalization `Omega = epsilon*Lambda^2`, "
        "`r = 2*epsilon*L_segment`, and `L_segment = 1`, this becomes",
        "",
        "```text",
        "Omega_c = sqrt(kappa*G/rho)/epsilon",
        "Lambda_c = (kappa/(2*(1 + nu)))**0.25/epsilon",
        "```",
        "",
        "The CSV records `Omega_Timoshenko`, `Omega_cutoff`, `Omega_over_cutoff`, "
        "`Lambda_cutoff`, `Lambda_over_cutoff`, and `below_cutoff`.",
        "",
        "Minimum cut-off margin `1 - Omega/Omega_c`:",
        "",
        "```text",
        f"{min_margin:.6e}",
        "```",
        "",
        "Closest row to cut-off:",
        "",
        "```text",
        (
            f"epsilon={float(closest['epsilon']):g}, mode={int(closest['mode'])}, "
            f"Omega/Omega_c={max_ratio:.6e}, Lambda/Lambda_c={float(closest['Lambda_over_cutoff']):.6e}"
        ),
        "```",
        "",
    ]

    if not at_or_above_rows:
        lines.append("All computed Timoshenko roots lie below the cut-off frequency.")
    else:
        lines.append(
            "ERROR: at least one computed Timoshenko root is at or above the cut-off frequency."
        )
    if near_rows:
        lines.append(
            f"WARNING: {len(near_rows)} computed roots have Omega/Omega_c >= {CUTOFF_WARNING_RATIO:g}."
        )
    else:
        lines.append(f"No computed root has Omega/Omega_c >= {CUTOFF_WARNING_RATIO:g}.")

    return lines


def kappa_sensitivity_report_lines(rows: Sequence[dict[str, object]]) -> list[str]:
    comparison_rows = [
        row
        for row in rows
        if row["kappa_label"] != "baseline_circular"
        and np.isfinite(float(row["delta_Lambda_vs_baseline"]))
    ]
    if not comparison_rows:
        return [
            "## Kappa Sensitivity",
            "",
            "No finite non-baseline kappa sensitivity rows were computed.",
        ]

    max_delta_row = max(comparison_rows, key=lambda row: abs(float(row["delta_Lambda_vs_baseline"])))
    max_rel_row = max(comparison_rows, key=lambda row: abs(float(row["rel_delta_Lambda_vs_baseline"])))
    finite_baseline_diffs = [
        abs(float(row["baseline_Timoshenko_minus_EB"]))
        for row in comparison_rows
        if np.isfinite(float(row["baseline_Timoshenko_minus_EB"]))
    ]
    max_baseline_diff = max(finite_baseline_diffs) if finite_baseline_diffs else float("nan")
    max_delta = abs(float(max_delta_row["delta_Lambda_vs_baseline"]))
    factor = max_delta / max_baseline_diff if np.isfinite(max_baseline_diff) and max_baseline_diff != 0.0 else float("nan")
    relation = "smaller than" if np.isfinite(factor) and factor < 1.0 else "larger than or comparable to"

    lines = [
        "## Kappa Sensitivity",
        "",
        "A small diagnostic sensitivity check was run for `epsilon = 0.025, 0.05` "
        "and the first six sorted modes.",
        "",
        "Kappa options:",
        "",
        f"- baseline circular: `{circular_shear_coefficient(NU):.12g}`",
        "- unit shear coefficient: `1.0`",
        "",
        "A Cowper-like circular option is not included here because no Cowper "
        "circular formula has been promoted to the local literature notes as a "
        "project diagnostic choice.",
        "",
        f"Output CSV: `{OUTPUT_SENSITIVITY_CSV.relative_to(REPO_ROOT)}`",
        "",
        "Largest absolute Lambda change from changing kappa:",
        "",
        "```text",
        (
            f"{max_delta:.6e} at epsilon={float(max_delta_row['epsilon']):g}, "
            f"mode={int(max_delta_row['mode'])}, option={max_delta_row['kappa_label']}"
        ),
        "```",
        "",
        "Largest relative Lambda change from changing kappa:",
        "",
        "```text",
        (
            f"{abs(float(max_rel_row['rel_delta_Lambda_vs_baseline'])):.6e} "
            f"at epsilon={float(max_rel_row['epsilon']):g}, mode={int(max_rel_row['mode'])}"
        ),
        "```",
        "",
        "For the same epsilon/mode set, the maximum baseline `|Lambda_Timoshenko - Lambda_EB|` is",
        "",
        "```text",
        f"{max_baseline_diff:.6e}",
        "```",
        "",
        (
            f"The kappa sensitivity is {relation} the EB-vs-Timoshenko difference "
            f"for this diagnostic subset (ratio `{factor:.6g}`)."
        ),
        "",
        "`kappa` choice is a modeling choice, not a unique material constant.",
    ]
    return lines


def report_warning_lines(warnings: Sequence[str]) -> list[str]:
    if not warnings:
        return ["No Timoshenko basis or root-search warnings were recorded."]
    lines = ["Numerical/root-search warnings:"]
    for warning in list(warnings)[:20]:
        lines.append(f"- {warning}")
    if len(warnings) > 20:
        lines.append(f"- ... {len(warnings) - 20} more warnings")
    return lines


def write_report(
    rows: Sequence[dict[str, object]],
    warnings: Sequence[str],
    sensitivity_rows: Sequence[dict[str, object]],
    sensitivity_warnings: Sequence[str],
) -> None:
    kappa = circular_shear_coefficient(NU)
    valid_rows = [row for row in rows if bool(row["thin_rod_valid"]) and np.isfinite(float(row["rel_diff"]))]
    finite_rows = [row for row in rows if np.isfinite(float(row["rel_diff"]))]
    max_valid = max(valid_rows, key=lambda row: abs(float(row["rel_diff"]))) if valid_rows else None
    max_full = max(finite_rows, key=lambda row: abs(float(row["rel_diff"]))) if finite_rows else None
    small_rows = [row for row in rows if abs(float(row["epsilon"]) - min(EPSILON_VALUES)) <= 1.0e-15]
    small_max = max((abs(float(row["rel_diff"])) for row in small_rows if np.isfinite(float(row["rel_diff"]))), default=float("nan"))
    invalid_eps = [epsilon for epsilon in EPSILON_VALUES if not thin_rod_valid(float(epsilon))]

    lines: list[str] = []
    lines.append("# Coupled Equal Rods EB/Timoshenko Diagnostic")
    lines.append("")
    lines.append("Diagnostic only. This is a sorted frequency comparison, not descendant branch tracking.")
    lines.append("No determinant, formulas.py, old solver, FEM physical model, article file, or article figure was changed.")
    lines.append("")
    lines.append("## Configuration")
    lines.append("")
    lines.append(f"- beta: {BETA_DEG:g} deg")
    lines.append(f"- mu: {MU:g}")
    lines.append(f"- eta: 0")
    lines.append(f"- tau1 = tau2 = 1")
    lines.append(f"- L_segment: {L_SEGMENT:g}")
    lines.append(f"- epsilon values: {', '.join(f'{float(value):g}' for value in EPSILON_VALUES)}")
    lines.append("- project normalization for equal rods: Omega = epsilon*Lambda^2")
    lines.append("- diameter-to-segment-length ratio: 2r/L_segment = 4*epsilon")
    lines.append(f"- thin-rod criterion: 4*epsilon <= {THICKNESS_RATIO_LIMIT:g}, so epsilon <= {THICKNESS_RATIO_LIMIT / 4.0:g}")
    lines.append("")
    lines.append("## Euler-Bernoulli Baseline")
    lines.append("")
    lines.append(
        "The EB roots are the first sorted roots of "
        "`src/my_project/analytic/formulas.py::det_clamped_coupled(Lambda, beta, mu, eps)` "
        "found through the existing `find_first_n_roots` helper."
    )
    lines.append("The old determinant and helper implementation were used as source-of-truth baseline inputs and were not edited.")
    lines.append("")
    lines.append("## Timoshenko Trial Model")
    lines.append("")
    lines.append("For each rod, the longitudinal field is `u_i(x) = P_i sin(theta*x)` with `theta = epsilon*Lambda^2`.")
    lines.append("For bending, `K = kappa*G*A`, `B = E*I`, `m = rho*A`, `j = rho*I`, and `Omega = epsilon*Lambda^2`.")
    lines.append("")
    lines.append("The Timoshenko spatial parameters come from")
    lines.append("")
    lines.append("```text")
    lines.append("B*lambda^4 + Omega^2*(B*m/K + j)*lambda^2 + (j*m*Omega^4/K - m*Omega^2) = 0")
    lines.append("```")
    lines.append("")
    lines.append("Below cut-off the expected roots are one positive and one negative value of `lambda^2`, giving `a` and `b`.")
    lines.append("The rotation-to-deflection ratios are `h = a + m*Omega^2/(K*a)` and `q = b - m*Omega^2/(K*b)`.")
    lines.append("")
    lines.append("The clamped-end bending basis is")
    lines.append("")
    lines.append("```text")
    lines.append("W_A(x)   = cosh(a*x) - cos(b*x)")
    lines.append("W_B(x)   = sinh(a*x) - (h/q)*sin(b*x)")
    lines.append("PSI_A(x) = h*sinh(a*x) + q*sin(b*x)")
    lines.append("PSI_B(x) = h*(cosh(a*x) - cos(b*x))")
    lines.append("```")
    lines.append("")
    lines.append("## Coupling Conditions")
    lines.append("")
    lines.append(
        "The first two displacement compatibility conditions are unchanged relative to the existing determinant's "
        "parity-adjusted `x2 = -1` column convention:"
    )
    lines.append("")
    lines.append("```text")
    lines.append("w1(1) = w2(-1)*cos(beta) - u2(-1)*sin(beta)")
    lines.append("u1(1) = w2(-1)*sin(beta) + u2(-1)*cos(beta)")
    lines.append("```")
    lines.append("")
    lines.append(
        "This sign form is intentionally tied to `assemble_clamped_coupled_matrix`; using the literal pre-substitution "
        "kinematic lines with the explicit `sin(theta*x2)` at `x2=-1` does not recover the current EB determinant."
    )
    lines.append("")
    lines.append("The Euler-Bernoulli slope condition is replaced by independent Timoshenko rotation continuity:")
    lines.append("")
    lines.append("```text")
    lines.append("psi1(1) = psi2(-1)")
    lines.append("M1(1) = M2(-1)")
    lines.append("-Q1(1) = -Q2(-1)*cos(beta) + N2(-1)*sin(beta)")
    lines.append("N1(1) = Q2(-1)*sin(beta) + N2(-1)*cos(beta)")
    lines.append("```")
    lines.append("")
    lines.append("Here `N_i = E*A*u_i'`, `M_i = E*I*psi_i'`, and `Q_i = kappa*G*A*(w_i' - psi_i)`. For root scanning the six rows are normalized by positive row norms; this changes scaling only, not the zero set.")
    lines.append("")
    lines.append("## Flexural And Axial Model Scope")
    lines.append("")
    lines.append("The current diagnostic Timoshenko extension modifies only the flexural part:")
    lines.append("- bending moment: `M = E I psi'`")
    lines.append("- shear force: `Q = kappa G A (w' - psi)`")
    lines.append("")
    lines.append("The axial part remains classical:")
    lines.append("- `N = E A u'`")
    lines.append("")
    lines.append(
        "This is acceptable for the current low flexural-frequency diagnostics, but higher-frequency "
        "longitudinal refinements may require Mindlin-Herrmann / related rod models."
    )
    lines.append("")
    lines.append("## Shear Coefficient")
    lines.append("")
    lines.append("The circular-rod shear coefficient is used as the working baseline from the circular Timoshenko-rod sources:")
    lines.append("")
    lines.append("```text")
    lines.append("kappa = (6 + 12*nu + 6*nu^2)/(7 + 12*nu + 4*nu^2)")
    lines.append("```")
    lines.append("")
    lines.append(f"For `nu = {NU:g}`, `kappa = {kappa:.12g}`. This is a modeling choice for the diagnostic circular rods, not a claim of uniqueness.")
    lines.append("")
    lines.extend(cutoff_report_lines(rows))
    lines.append("")
    lines.extend(kappa_sensitivity_report_lines(sensitivity_rows))
    if sensitivity_warnings:
        lines.append("")
        lines.append("Kappa sensitivity numerical/root-search warnings:")
        for warning in list(sensitivity_warnings)[:20]:
            lines.append(f"- {warning}")
        if len(sensitivity_warnings) > 20:
            lines.append(f"- ... {len(sensitivity_warnings) - 20} more warnings")
    lines.append("")
    lines.append("## Numeric Summary")
    lines.append("")
    lines.append(markdown_numeric_summary(rows))
    lines.append("")
    if max_valid is not None:
        lines.append(
            "Maximum relative difference in the valid region "
            f"(`4*epsilon <= 0.1`) is {abs(float(max_valid['rel_diff'])):.6g} "
            f"at epsilon={float(max_valid['epsilon']):g}, mode {int(max_valid['mode'])}."
        )
    if max_full is not None:
        lines.append(
            "Maximum relative difference in the full tested range is "
            f"{abs(float(max_full['rel_diff'])):.6g} "
            f"at epsilon={float(max_full['epsilon']):g}, mode {int(max_full['mode'])}."
        )
    lines.append(
        f"At the smallest tested epsilon={min(EPSILON_VALUES):g}, the maximum sorted-mode relative difference is {small_max:.6g}, "
        "so the Timoshenko roots tend toward the Euler-Bernoulli baseline in the thin limit for this diagnostic."
    )
    lines.append("")
    lines.append("## Applicability Warning")
    lines.append("")
    if invalid_eps:
        lines.append(
            "The following epsilon values are outside the current thin-rod applicability criterion "
            f"`4*epsilon <= 0.1`: {', '.join(f'{float(value):g}' for value in invalid_eps)}."
        )
    else:
        lines.append("All tested epsilon values satisfy the current thin-rod applicability criterion.")
    lines.append("")
    lines.extend(report_warning_lines(warnings))
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    lines.append(f"- CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`")
    lines.append(f"- kappa sensitivity CSV: `{OUTPUT_SENSITIVITY_CSV.relative_to(REPO_ROOT)}`")
    lines.append(f"- PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`")
    lines.append(f"- report: `{OUTPUT_REPORT.relative_to(REPO_ROOT)}`")
    lines.append("")
    OUTPUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    rows, warnings = compute_rows()
    sensitivity_rows, sensitivity_warnings = compute_kappa_sensitivity_rows()
    write_csv(rows)
    write_kappa_sensitivity_csv(sensitivity_rows)
    plot_rows(rows)
    write_report(rows, warnings, sensitivity_rows, sensitivity_warnings)

    print(f"saved CSV: {OUTPUT_CSV}")
    print(f"saved kappa sensitivity CSV: {OUTPUT_SENSITIVITY_CSV}")
    print(f"saved plot: {OUTPUT_PNG}")
    print(f"saved report: {OUTPUT_REPORT}")
    print(markdown_numeric_summary(rows))
    if warnings:
        print("warnings:")
        for warning in warnings[:10]:
            print(f"- {warning}")
        if len(warnings) > 10:
            print(f"- ... {len(warnings) - 10} more warnings")


if __name__ == "__main__":
    main()
