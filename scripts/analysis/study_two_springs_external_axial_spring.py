"""Stable diagnostic entry point for the external local axial support."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.lib.two_springs_external_axial_spring import main

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('preflight','run','plot_only','finalize'),default='run')
    parser.add_argument('--folder',type=Path)
    main(parser.parse_args())
