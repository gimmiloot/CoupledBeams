"""Fixed-volume individual-section RL_FULL/Timoshenko study.

Script Proliferation Control: the new two-section physical/output contract needs
individual pole catalogs and weighted full/leading field checks. This stable
entry point reuses the old Run and numerical helpers without changing the old
CLI fingerprint or invalidating its existing immutable research bundles.
"""
import argparse
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['gates','preflight','all'],default='all')
    parser.add_argument('--plot-only',type=Path)
    parser.add_argument('--verify-cache',type=Path)
    args=parser.parse_args()
    from scripts.lib.two_springs_thickness_asymmetry import main as run
    run(args)


if __name__=='__main__':main()
