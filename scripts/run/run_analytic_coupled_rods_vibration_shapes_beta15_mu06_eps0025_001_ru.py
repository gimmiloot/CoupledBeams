from __future__ import annotations

from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.analysis.plot_desc05_full_shapes_beta15_eps_sweep import main as run_desc05_shape_sweep  # noqa: E402


# ============================================================
# USER PARAMETERS
# Change these values for ordinary runs.
# ============================================================

BRANCH_NUMBER = 5
BETA = 45.0
MU = 0.6

# Each tuple is (epsilon, legend root label). The root label is only used in
# the legend; the branch identity is still tracked as descendant BRANCH_NUMBER.
EPSILON_ROOT_CASES = (
    (0.0025, 6),
    (0.01, 5),
)

TITLE_PREFIX = "Формы колебаний: потомок 5-й изгибной ветви"
OUTPUT = "results/vibration_shapes_descendant5_beta15_mu0p6_eps0p0025_0p01.png"
L_TOTAL = 2.0
MODE_SCALE = 0.12
DPI = 240
SHAPE_METRIC = "full"
HIDE_GEOMETRY_IN_LEGEND = True
ALLOW_LOW_MAC = False

BETA_STEPS = 200
MU_STEPS = 260
MAX_REFINEMENT_DEPTH = 8
MIN_BETA_STEP = 1e-3
MIN_MU_STEP = 1e-4


def build_args() -> list[str]:
    epsilons = [epsilon for epsilon, _root_label in EPSILON_ROOT_CASES]
    root_labels = [root_label for _epsilon, root_label in EPSILON_ROOT_CASES]

    args = [
        "--branch-number",
        f"{int(BRANCH_NUMBER)}",
        "--beta",
        f"{float(BETA):g}",
        "--mus",
        f"{float(MU):g}",
        "--epsilons",
        *(f"{float(epsilon):g}" for epsilon in epsilons),
        "--case-root-labels",
        *(f"{int(root_label)}" for root_label in root_labels),
        "--title-prefix",
        TITLE_PREFIX,
        "--output",
        OUTPUT,
        "--l-total",
        f"{float(L_TOTAL):g}",
        "--mode-scale",
        f"{float(MODE_SCALE):g}",
        "--dpi",
        f"{int(DPI)}",
        "--shape-metric",
        SHAPE_METRIC,
        "--beta-steps",
        f"{int(BETA_STEPS)}",
        "--mu-steps",
        f"{int(MU_STEPS)}",
        "--max-refinement-depth",
        f"{int(MAX_REFINEMENT_DEPTH)}",
        "--min-beta-step",
        f"{float(MIN_BETA_STEP):g}",
        "--min-mu-step",
        f"{float(MIN_MU_STEP):g}",
    ]
    if HIDE_GEOMETRY_IN_LEGEND:
        args.append("--hide-geometry-legend")
    if ALLOW_LOW_MAC:
        args.append("--allow-low-mac")
    return args


def main() -> list[dict[str, float | int | str]]:
    return run_desc05_shape_sweep(build_args())


if __name__ == "__main__":
    main()
