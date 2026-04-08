from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.FreqFromMu import BeamParams  # noqa: E402
from my_project.analytic.formulas import lambdas_to_frequencies  # noqa: E402
from my_project.analytic.solvers import tracked_lambdas_vs_mu  # noqa: E402


MODEL_LABEL = "\u0441\u043f\u043b\u043e\u0448\u043d\u0430\u044f \u2014 \u043c\u043e\u0434\u0435\u043b\u044c"
FEM_LABEL = "\u043a\u0440\u0443\u0436\u043a\u0438 \u2014 FEM"
MODE_LABEL_TEMPLATE = "\u043c\u043e\u0434\u0430 {}"
X_LABEL = "\u03bc"
Y_LABEL = "\u0427\u0430\u0441\u0442\u043e\u0442\u0430, \u0413\u0446"
TITLE_LINE_1 = "\u0421\u0440\u0430\u0432\u043d\u0435\u043d\u0438\u0435 \u0447\u0430\u0441\u0442\u043e\u0442: \u043c\u043e\u0434\u0435\u043b\u044c \u0441\u043e\u043f\u0440\u044f\u0436\u0451\u043d\u043d\u044b\u0445 \u0441\u0442\u0435\u0440\u0436\u043d\u0435\u0439 \u0438 FEM"


def load_fem_frequency_table(csv_path: Path, beta_deg: float, n_modes: int) -> tuple[np.ndarray, np.ndarray]:
    beta_suffix = f"beta{int(beta_deg)}" if float(beta_deg).is_integer() else f"beta{beta_deg}"
    data = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=float, encoding=None)

    mu_values = np.asarray(data["mu"], dtype=float)
    freq_columns = [f"f{i + 1}_Hz_{beta_suffix}" for i in range(n_modes)]
    missing = [name for name in freq_columns if name not in data.dtype.names]
    if missing:
        raise KeyError(f"Missing FEM columns in {csv_path}: {missing}")

    fem_freqs = np.vstack([np.asarray(data[name], dtype=float) for name in freq_columns])
    return mu_values, fem_freqs


def compute_analytic_frequencies(
    params: BeamParams,
    beta_deg: float,
    mu_values: np.ndarray,
    n_modes: int,
) -> np.ndarray:
    lambdas_tr = tracked_lambdas_vs_mu(
        params=params,
        beta_deg=beta_deg,
        mu_values=mu_values,
        n_modes=n_modes,
        Lmin=0.2,
        Lmax0=55.0,
        scan_step=0.02,
        grow_factor=1.35,
        max_tries=8,
        tracking_method="auto",
    )
    return lambdas_to_frequencies(lambdas_tr, params)


def plot_comparison(
    mu_values: np.ndarray,
    analytic_freqs: np.ndarray,
    fem_freqs: np.ndarray,
    beta_deg: float,
    params: BeamParams,
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(11.2, 6.6))

    mode_handles = []
    for i in range(analytic_freqs.shape[0]):
        (line,) = ax.plot(mu_values, analytic_freqs[i], linewidth=2.0)
        color = line.get_color()
        ax.plot(
            mu_values,
            fem_freqs[i],
            linestyle="None",
            marker="o",
            color=color,
            markerfacecolor=color,
            markeredgecolor=color,
            markersize=3.8,
        )
        mode_handles.append(
            Line2D([0], [0], color=color, lw=2.0, label=MODE_LABEL_TEMPLATE.format(i + 1))
        )

    style_handles = [
        Line2D([0], [0], color="black", lw=2.0, label=MODEL_LABEL),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=5, label=FEM_LABEL),
    ]

    ax.set_xlabel(X_LABEL)
    ax.set_ylabel(Y_LABEL)
    ax.set_title(
        TITLE_LINE_1 + "\n" + f"\u03b2={beta_deg:.1f}\u00b0, L_total={params.L_total:.1f}, r={params.r:.3f}"
    )
    ax.grid(True)
    ax.legend(handles=style_handles + mode_handles, fontsize=9, ncols=2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.show()


def main() -> None:
    beta_deg = 15.0
    n_modes = 6
    params = BeamParams(E=2.1e11, rho=7800.0, r=0.005, L_total=2.0)

    csv_path = REPO_ROOT / "results" / "fem_spectrum.csv"
    output_path = REPO_ROOT / "results" / "freq_mu_comparison_beta15.png"

    mu_values, fem_freqs = load_fem_frequency_table(csv_path, beta_deg=beta_deg, n_modes=n_modes)
    analytic_freqs = compute_analytic_frequencies(params, beta_deg=beta_deg, mu_values=mu_values, n_modes=n_modes)

    plot_comparison(
        mu_values=mu_values,
        analytic_freqs=analytic_freqs,
        fem_freqs=fem_freqs,
        beta_deg=beta_deg,
        params=params,
        output_path=output_path,
    )

    print(f"FEM CSV: {csv_path}")
    print(f"Saved plot: {output_path}")
    print(f"Compared modes: {n_modes}")
    print("Analytic curves: tracked frequencies from the coupled-beam model at beta=15 deg.")
    print("FEM markers: tracked FEM frequencies from columns f*_Hz_beta15 in fem_spectrum.csv.")


if __name__ == "__main__":
    main()