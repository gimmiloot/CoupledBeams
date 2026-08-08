"""Build the versioned zero-solve epsilon upper-envelope article bundle."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.lib.article_epsilon_upper_envelope_assets import generate_article_assets


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate article-ready plots, tables, and text from compact isotropic-circular "
            "EB/Timoshenko results without scientific solves."
        )
    )
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--language", choices=("ru",), default="ru")
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    recorded_arguments = list(sys.argv[1:] if argv is None else argv)
    result = generate_article_assets(
        args.source_dir,
        args.output_dir,
        language=args.language,
        dpi=args.dpi,
        overwrite=args.overwrite,
        generator_arguments=recorded_arguments,
    )
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
