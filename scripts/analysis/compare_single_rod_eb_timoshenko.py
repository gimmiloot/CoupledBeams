from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq


# =========================
# User-editable parameters
# =========================
L = 1.0
NU = 0.3
E = 1.0
RHO = 1.0
EPSILON_VALUES = [0.0025, 0.01, 0.025, 0.05, 0.075, 0.1]
N_MODES = 6
THICKNESS_RATIO_LIMIT = 0.1


REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_CSV = OUTPUT_DIR / "single_rod_fixed_fixed_eb_timoshenko_frequencies.csv"
OUTPUT_PNG = OUTPUT_DIR / "single_rod_fixed_fixed_eb_timoshenko_lambda_vs_epsilon.png"
OUTPUT_REPORT = OUTPUT_DIR / "single_rod_fixed_fixed_eb_timoshenko_report.md"

ROOT_SCAN_POINTS = 25000
ROOT_SCAN_START = 1.0e-12
ROOT_DEDUP_TOL = 1.0e-8
CUTOFF_WARNING_RATIO = 0.8


@dataclass(frozen=True)
class Section:
    radius: float
    area: float
    inertia: float
    shear_modulus: float
    kappa: float


def circular_shear_coefficient(nu: float) -> float:
    """Diaz-de-Anda et al. circular-rod coefficient."""
    return (6.0 + 12.0 * nu + 6.0 * nu**2) / (7.0 + 12.0 * nu + 4.0 * nu**2)


def section_from_epsilon(epsilon: float) -> Section:
    radius = epsilon * L
    area = np.pi * radius**2
    inertia = np.pi * radius**4 / 4.0
    shear_modulus = E / (2.0 * (1.0 + NU))
    return Section(
        radius=radius,
        area=area,
        inertia=inertia,
        shear_modulus=shear_modulus,
        kappa=circular_shear_coefficient(NU),
    )


def fixed_fixed_eb_roots(n_modes: int) -> np.ndarray:
    """Positive roots of cos(alpha) cosh(alpha) = 1."""
    roots: list[float] = []

    def characteristic(alpha: float) -> float:
        return float(np.cos(alpha) * np.cosh(alpha) - 1.0)

    for mode in range(1, n_modes + 1):
        center = (mode + 0.5) * np.pi
        left = center - 0.25 * np.pi
        right = center + 0.25 * np.pi
        roots.append(brentq(characteristic, left, right, xtol=1.0e-14, rtol=1.0e-14))
    return np.array(roots, dtype=float)


def omega_eb(alpha: float, section: Section) -> float:
    return alpha**2 / L**2 * np.sqrt(E * section.inertia / (RHO * section.area))


def lambda_from_omega(omega: float, section: Section) -> float:
    return np.sqrt(omega) * (L / 2.0) * (RHO * section.area / (E * section.inertia)) ** 0.25


def omega_cutoff(section: Section) -> float:
    return float(np.sqrt(section.kappa * section.shear_modulus * section.area / (RHO * section.inertia)))


def timoshenko_state_matrix(omega: float, section: Section) -> np.ndarray:
    """State matrix for y = [w, psi, V, M].

    Convention:
        V = kappa*G*A*(w' - psi)
        M = E*I*psi'

    Harmonic motion with angular frequency omega gives:
        w'   = psi + V/(kappa*G*A)
        psi' = M/(E*I)
        V'   = -rho*A*omega**2*w
        M'   = -V - rho*I*omega**2*psi
    """
    kga = section.kappa * section.shear_modulus * section.area
    ei = E * section.inertia
    return np.array(
        [
            [0.0, 1.0, 1.0 / kga, 0.0],
            [0.0, 0.0, 0.0, 1.0 / ei],
            [-RHO * section.area * omega**2, 0.0, 0.0, 0.0],
            [0.0, -RHO * section.inertia * omega**2, -1.0, 0.0],
        ],
        dtype=float,
    )


def timoshenko_boundary_determinant(omega: float, section: Section) -> float:
    transfer = expm(timoshenko_state_matrix(omega, section) * L)
    block = transfer[:2, 2:4]
    scale = max(float(np.linalg.norm(block[0, :]) * np.linalg.norm(block[1, :])), 1.0e-300)
    return float(np.linalg.det(block) / scale)


def find_timoshenko_roots(section: Section, eb_omegas: np.ndarray, n_modes: int) -> np.ndarray:
    """Find the first clamped-clamped Timoshenko roots by sign-change scanning."""
    roots: list[float] = []

    for factor in (1.1, 1.5, 2.0, 3.0):
        roots.clear()
        previous_root = -np.inf
        upper = max(float(eb_omegas[n_modes - 1] * factor), 1.0)
        grid = np.linspace(ROOT_SCAN_START, upper, ROOT_SCAN_POINTS)
        x_prev = float(grid[0])
        y_prev = timoshenko_boundary_determinant(x_prev, section)

        for x in grid[1:]:
            x = float(x)
            y = timoshenko_boundary_determinant(x, section)
            if not np.isfinite(y_prev) or not np.isfinite(y):
                x_prev, y_prev = x, y
                continue

            if y_prev == 0.0:
                candidate = x_prev
            elif y_prev * y < 0.0:
                candidate = brentq(
                    lambda omega: timoshenko_boundary_determinant(omega, section),
                    x_prev,
                    x,
                    xtol=1.0e-12,
                    rtol=1.0e-12,
                    maxiter=100,
                )
            else:
                x_prev, y_prev = x, y
                continue

            if candidate - previous_root > ROOT_DEDUP_TOL:
                roots.append(float(candidate))
                previous_root = float(candidate)
            if len(roots) >= n_modes:
                return np.array(roots[:n_modes], dtype=float)

            x_prev, y_prev = x, y

    raise RuntimeError(
        f"Found only {len(roots)} Timoshenko roots below {upper:g}; "
        f"increase ROOT_SCAN_POINTS or the scan factor."
    )


def compute_rows() -> list[dict[str, float | int | bool]]:
    alpha_roots = fixed_fixed_eb_roots(N_MODES)
    rows: list[dict[str, float | int | bool]] = []

    for epsilon in EPSILON_VALUES:
        section = section_from_epsilon(float(epsilon))
        eb_omegas = np.array([omega_eb(alpha, section) for alpha in alpha_roots], dtype=float)
        eb_lambdas = np.array([lambda_from_omega(omega, section) for omega in eb_omegas], dtype=float)
        timo_omegas = find_timoshenko_roots(section, eb_omegas, N_MODES)
        timo_lambdas = np.array([lambda_from_omega(omega, section) for omega in timo_omegas], dtype=float)
        omega_c = omega_cutoff(section)
        lambda_c = lambda_from_omega(omega_c, section)
        diameter_to_length = 2.0 * float(epsilon)
        thin_rod_valid = diameter_to_length <= THICKNESS_RATIO_LIMIT

        for mode_index in range(N_MODES):
            omega_ratio = float(timo_omegas[mode_index] / omega_c)
            lambda_ratio = float(timo_lambdas[mode_index] / lambda_c)
            rows.append(
                {
                    "epsilon": float(epsilon),
                    "diameter_to_length": diameter_to_length,
                    "mode": mode_index + 1,
                    "alpha_EB": float(alpha_roots[mode_index]),
                    "Omega_EB": float(eb_omegas[mode_index]),
                    "Lambda_EB": float(eb_lambdas[mode_index]),
                    "Omega_Timoshenko": float(timo_omegas[mode_index]),
                    "Omega_cutoff": omega_c,
                    "Omega_over_cutoff": omega_ratio,
                    "Lambda_Timoshenko": float(timo_lambdas[mode_index]),
                    "Lambda_cutoff": lambda_c,
                    "Lambda_over_cutoff": lambda_ratio,
                    "below_cutoff": omega_ratio < 1.0,
                    "rel_diff_Lambda": float((timo_lambdas[mode_index] - eb_lambdas[mode_index]) / eb_lambdas[mode_index]),
                    "thin_rod_valid": thin_rod_valid,
                }
            )

    return rows


def write_csv(rows: list[dict[str, float | int | bool]]) -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "epsilon",
        "diameter_to_length",
        "mode",
        "alpha_EB",
        "Omega_EB",
        "Lambda_EB",
        "Omega_Timoshenko",
        "Omega_cutoff",
        "Omega_over_cutoff",
        "Lambda_Timoshenko",
        "Lambda_cutoff",
        "Lambda_over_cutoff",
        "below_cutoff",
        "rel_diff_Lambda",
        "thin_rod_valid",
    ]
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_rows(rows: list[dict[str, float | int | bool]]) -> None:
    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    mode_colors = plt.cm.viridis(np.linspace(0.08, 0.92, N_MODES))

    x_min = min(float(row["diameter_to_length"]) for row in rows)
    x_max = max(float(row["diameter_to_length"]) for row in rows)
    if x_max > THICKNESS_RATIO_LIMIT:
        ax.axvspan(THICKNESS_RATIO_LIMIT, x_max, color="0.9", alpha=0.75, label="outside 2r/L <= 0.1")
    ax.axvline(THICKNESS_RATIO_LIMIT, color="0.45", linestyle=":", linewidth=1.1)

    for mode in range(1, N_MODES + 1):
        mode_rows = [row for row in rows if int(row["mode"]) == mode]
        x = np.array([float(row["diameter_to_length"]) for row in mode_rows])
        y_timo = np.array([float(row["Lambda_Timoshenko"]) for row in mode_rows])
        y_eb = np.array([float(row["Lambda_EB"]) for row in mode_rows])
        valid = x <= THICKNESS_RATIO_LIMIT
        color = mode_colors[mode - 1]

        ax.hlines(
            y_eb[0],
            x_min,
            x_max,
            colors=[color],
            linestyles="--",
            linewidth=1.0,
            alpha=0.55,
        )
        ax.plot(x[valid], y_timo[valid], marker="o", color=color, linewidth=1.7, label=f"mode {mode}")
        if np.any(~valid):
            ax.plot(
                x[~valid],
                y_timo[~valid],
                marker="o",
                markerfacecolor="white",
                color=color,
                linestyle=":",
                linewidth=1.7,
            )

    ax.set_xlabel(r"diameter-to-length ratio $2r/L$")
    ax.set_ylabel(r"project frequency parameter $\Lambda$")
    ax.set_title("Single fixed-fixed circular rod: EB vs Timoshenko")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=200)
    plt.close(fig)


def format_summary_table(rows: list[dict[str, float | int | bool]]) -> str:
    lines = [
        "| epsilon | 2r/L | max rel diff Lambda | min cut-off margin | valid? |",
        "|---:|---:|---:|---:|:---:|",
    ]
    for epsilon in EPSILON_VALUES:
        eps_rows = [row for row in rows if np.isclose(float(row["epsilon"]), float(epsilon))]
        max_abs = max(abs(float(row["rel_diff_Lambda"])) for row in eps_rows)
        min_margin = min(1.0 - float(row["Omega_over_cutoff"]) for row in eps_rows)
        valid = "yes" if bool(eps_rows[0]["thin_rod_valid"]) else "no"
        lines.append(f"| {epsilon:.4g} | {2.0 * epsilon:.4g} | {max_abs:.6e} | {min_margin:.6e} | {valid} |")
    return "\n".join(lines)


def cutoff_report_text(rows: list[dict[str, float | int | bool]]) -> str:
    finite_rows = [
        row
        for row in rows
        if np.isfinite(float(row["Omega_over_cutoff"])) and np.isfinite(float(row["Omega_cutoff"]))
    ]
    if not finite_rows:
        return "## Cut-Off Frequency Check\n\nNo finite cut-off ratios were computed."

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
        "For this single circular rod with `r = epsilon*L`, this is also",
        "",
        "```text",
        "Omega_c = 2*sqrt(kappa*G/rho)/(epsilon*L)",
        "Lambda_c = (kappa/(2*(1 + nu)))**0.25/epsilon",
        "```",
        "",
        "The CSV records `Omega_cutoff`, `Omega_over_cutoff`, `Lambda_cutoff`, "
        "`Lambda_over_cutoff`, and `below_cutoff` for each Timoshenko root.",
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

    return "\n".join(lines)


def write_report(rows: list[dict[str, float | int | bool]]) -> None:
    alpha_roots = fixed_fixed_eb_roots(N_MODES)
    kappa = circular_shear_coefficient(NU)
    eb_identity_errors = [abs(float(row["Lambda_EB"]) - float(row["alpha_EB"]) / 2.0) for row in rows]
    valid_rows = [row for row in rows if bool(row["thin_rod_valid"])]
    invalid_rows = [row for row in rows if not bool(row["thin_rod_valid"])]
    max_valid = max(abs(float(row["rel_diff_Lambda"])) for row in valid_rows)
    max_full = max(abs(float(row["rel_diff_Lambda"])) for row in rows)
    first_eps_rows = [row for row in rows if np.isclose(float(row["epsilon"]), float(EPSILON_VALUES[0]))]
    max_first_eps = max(abs(float(row["rel_diff_Lambda"])) for row in first_eps_rows)
    invalid_epsilons = sorted({float(row["epsilon"]) for row in invalid_rows})
    invalid_epsilon_text = ", ".join(f"{epsilon:g}" for epsilon in invalid_epsilons) if invalid_epsilons else "none"

    report = f"""# Single Fixed-Fixed Rod: Euler-Bernoulli vs Timoshenko

## Purpose

Diagnostic-only check for one straight circular rod of length `L = 2l` with
fixed-fixed boundary conditions. This is a preparation step before any
Timoshenko correction is considered for coupled beams.

No baseline determinant, old solver, article file, article figure, or FEM model
is used or modified by this calculation.

## Parameters

- `L = {L:g}`
- `E = {E:g}`
- `RHO = {RHO:g}`
- `NU = {NU:g}`
- `G = E/(2*(1 + NU)) = {E / (2.0 * (1.0 + NU)):.12g}`
- `epsilon = r/L`
- `diameter_to_length = 2*epsilon`
- tested epsilons: `{EPSILON_VALUES}`
- thin-rod validity criterion: `2r/L <= {THICKNESS_RATIO_LIMIT:g}`

## Shear Coefficient

For the circular section this diagnostic uses

```text
kappa = (6 + 12*NU + 6*NU**2)/(7 + 12*NU + 4*NU**2)
```

For `NU = {NU:g}`, this gives

```text
kappa = {kappa:.12g}
```

This coefficient is adopted from Diaz-de-Anda et al. for circular
aluminum Timoshenko rods. It is a working diagnostic choice, not a claim that
the shear coefficient is unique.

## Euler-Bernoulli Model

Fixed-fixed roots `alpha_n` solve

```text
cos(alpha_n)*cosh(alpha_n) = 1
```

The dimensional angular frequency and project-normalized frequency are

```text
Omega_EB = alpha_n**2/L**2 * sqrt(E*I/(RHO*A))
Lambda_EB = sqrt(Omega_EB)*(L/2)*(RHO*A/(E*I))**0.25
```

For `L = 2l`, the identity `Lambda_EB = alpha_n/2` should hold.
The maximum numerical error in this identity over the CSV rows is
`{max(eb_identity_errors):.3e}`.

EB roots used:

```text
{", ".join(f"{alpha:.12g}" for alpha in alpha_roots)}
```

## Timoshenko State-Space Model

The state is

```text
y = [w, psi, V, M]
```

with convention

```text
V = kappa*G*A*(w' - psi)
M = E*I*psi'
```

For harmonic frequency `Omega`, the first-order system is

```text
w'   = psi + V/(kappa*G*A)
psi' = M/(E*I)
V'   = -RHO*A*Omega**2*w
M'   = -V - RHO*I*Omega**2*psi
```

Boundary conditions:

```text
w(0) = 0, psi(0) = 0
w(L) = 0, psi(L) = 0
```

For each trial `Omega`, the transfer matrix `T = expm(B(Omega)*L)` is formed.
The determinant of the `2 x 2` block mapping initial `[V(0), M(0)]` to final
`[w(L), psi(L)]` is set to zero.

{cutoff_report_text(rows)}

## Flexural And Axial Model Scope

The current diagnostic Timoshenko extension modifies only the flexural part:
- bending moment: `M = E I psi'`
- shear force: `Q = kappa G A (w' - psi)`

The axial part remains classical:
- `N = E A u'`

This is acceptable for the current low flexural-frequency diagnostics, but
higher-frequency longitudinal refinements may require Mindlin-Herrmann /
related rod models.

## Numeric Summary

{format_summary_table(rows)}

Maximum absolute relative Lambda difference in the valid region
`2r/L <= {THICKNESS_RATIO_LIMIT:g}`:

```text
{max_valid:.6e}
```

Maximum absolute relative Lambda difference over the full tested range:

```text
{max_full:.6e}
```

At the smallest tested epsilon, `{EPSILON_VALUES[0]:g}`, the maximum absolute
relative Lambda difference is `{max_first_eps:.6e}`. This confirms that the
Timoshenko roots tend to the Euler-Bernoulli roots in the low-thickness limit
under the sign convention used here.

## Thin-Rod Warning

The rows with `epsilon > 0.05` have `2r/L > 0.1` and are outside the current
thin-rod applicability criterion. In this run those epsilons are:

```text
{invalid_epsilon_text}
```

Interpret those rows as stress tests of the diagnostic model, not as accepted
thin-rod predictions.
"""
    OUTPUT_REPORT.write_text(report, encoding="utf-8")


def main() -> None:
    rows = compute_rows()
    write_csv(rows)
    plot_rows(rows)
    write_report(rows)
    print(f"Wrote {OUTPUT_CSV}")
    print(f"Wrote {OUTPUT_PNG}")
    print(f"Wrote {OUTPUT_REPORT}")


if __name__ == "__main__":
    main()
