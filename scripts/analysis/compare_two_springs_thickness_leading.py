"""Compare saved full spectra to the unchanged individual-section leading model.

Separate read-full/solve-leading contract preserves immutable existing workflows.
No physical assembly or full-spectrum implementation is introduced here.
"""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['sources','regression','all'],default='all')
    parser.add_argument('--verify-cache',type=Path)
    parser.add_argument('--plot-only',type=Path)
    args=parser.parse_args()
    from scripts.lib.two_springs_thickness_leading_comparison import main as run
    run(args)


if __name__=='__main__':main()
