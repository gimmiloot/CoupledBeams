from __future__ import annotations

from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.plot_tracked_bending_descendant_shapes_ru import main  # noqa: E402


if __name__ == "__main__":
    main(
        [
            "--branch-id",
            "bending_desc_04",
            "--target-mus",
            "0.0",
            "0.1",
            "0.2",
            "--output",
            "results/flat_mu_bending_shapes_beta15_bending_desc_04_mu0_0p1_0p2_ru.png",
            "--axis-label-style",
            "xy",
        ]
    )
