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
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument('--high-frequency-ranges', action='store_true',
                        help='Extend archived axial spectra to 90 degrees; sampled 90/95 percent ranges')
    modes.add_argument('--small-stiffness-transition', action='store_true',
                       help='Bounded free-to-weak-axial-support experiment with archived endpoints')
    modes.add_argument('--transition-mechanism', action='store_true',
                       help='Three intermediate axial stiffnesses, extrema, sign brackets and full-model curvature')
    args=parser.parse_args()
    if args.transition_mechanism:
        from scripts.lib.two_springs_external_axial_mechanism import main as mechanism_main
        mechanism_main(args)
    elif args.small_stiffness_transition:
        from scripts.lib.two_springs_external_axial_transition import main as transition_main
        transition_main(args)
    elif args.high_frequency_ranges:
        from scripts.lib.two_springs_external_axial_high_ranges import main as ranges_main
        ranges_main(args)
    else:
        main(args)
