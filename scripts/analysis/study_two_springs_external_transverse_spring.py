"""Bounded transverse-support studies with retained or free outer rotation."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=('preflight','run','diagnose','plot_only','finalize'),default='run')
    parser.add_argument('--folder',type=Path)
    parser.add_argument('--leading-resume',action='store_true',help='Leading-only continuation with analytic curvature and read-only full spectra')
    parser.add_argument('--rotational-replacement',action='store_true',help='Free outer rotation: k_end=0, bounded full and leading study')
    parser.add_argument('--high-frequency-ranges',action='store_true',help='Read existing spectra and complete discrete 2/5-percent high ranges to 90 degrees')
    parser.add_argument('--high-frequency-ranges-90-95',action='store_true',help='Reuse completed 20-construction map; refine only new 10-percent boundaries')
    args=parser.parse_args()
    if sum([args.leading_resume,args.rotational_replacement,args.high_frequency_ranges,args.high_frequency_ranges_90_95])>1:
        parser.error('Select only one study mode')
    if args.high_frequency_ranges or args.high_frequency_ranges_90_95:
        from scripts.lib.two_springs_high_frequency_angle_ranges import main
    elif args.rotational_replacement:
        from scripts.lib.two_springs_external_rotational_replacement import main
    elif args.leading_resume:
        from scripts.lib.two_springs_external_transverse_leading_resume import main
    else:
        from scripts.lib.two_springs_external_transverse_spring import main
    main(args)
