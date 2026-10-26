"""Bounded diagnosis of the saved transverse-support Ritz strong residual."""
from pathlib import Path
import argparse
import json
import re
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from scripts.lib import two_springs_transverse_leading_diagnostic as diagnostic

def supplement(folder):
    """Local quadrature of the saved residual only; no eigenvalue calculation."""
    if not folder.resolve().is_relative_to(diagnostic.BASE.resolve()):
        raise ValueError('Supplement is restricted to the new diagnostic family')
    np=diagnostic.np
    solution=json.loads((folder/'original_failed_point.json').read_text())['solution']
    c=np.array(solution['coefficients'])[:,0];z=solution['zeta'][0];rows=[]
    x,w=np.polynomial.legendre.leggauss(64)
    for arm,L in enumerate((.99,1.01)):
        for lo,hi in zip([0,.05,.1,.9,.95],[.05,.1,.9,.95,1]):
            t=lo+(x+1)*(hi-lo)/2;f=diagnostic.fields(24,t,L,arm)
            W=f[0]@c[arm*24:(arm+1)*24]+(c[-1] if arm else 0.)
            r=f[4]@c[arm*24:(arm+1)*24]-z*W
            rows.append(dict(arm=arm+1,t0=lo,t1=hi,numerator_squared=float(L*(hi-lo)/2*(w@(r*r)))))
    io=diagnostic.api();io.csvwrite(folder/'localization.csv',rows)
    manifest=json.loads((folder/'manifest.json').read_text())
    manifest['artifacts']['localization.csv']=io.sha(folder/'localization.csv')
    io.write(folder/'manifest.json',manifest)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',type=Path,help='Verify saved artifacts only; no spectral solutions')
    parser.add_argument('--replica',help='Independent reproduction label; creates a new folder without overwriting evidence')
    parser.add_argument('--supplement',type=Path,help='Re-evaluate spatial contributions from saved coefficients only')
    args=parser.parse_args()
    if args.verify:
        print(diagnostic.verify_source(args.verify))
    elif args.supplement:
        supplement(args.supplement)
    else:
        if args.replica:
            if not re.fullmatch(r'[A-Za-z0-9_-]+',args.replica):parser.error('Use letters, digits, underscore or hyphen')
            diagnostic.BASE=diagnostic.BASE/'replicas'/args.replica
        supplement(diagnostic.run())
