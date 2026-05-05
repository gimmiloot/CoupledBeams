from __future__ import annotations

from pathlib import Path
import sys
from typing import Sequence


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.plot_mu_sweep_radius_fixed_four_betas_analytic import main as target_main  # noqa: E402


def main(argv: Sequence[str] | None = None) -> None:
    target_main(list(argv) if argv is not None else None)


if __name__ == "__main__":
    main(sys.argv[1:])
