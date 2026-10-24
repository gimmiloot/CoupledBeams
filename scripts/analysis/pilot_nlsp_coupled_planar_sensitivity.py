"""Scoped coupled-rod tangent/full nonlinear pilot; never launches 3D FEM."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
def main():
    if '--six-vs-seven-out-stability' in sys.argv[1:]:
        from scripts.lib.nlsp_six_vs_seven_out_stability import main as run
        return run([a for a in sys.argv[1:] if a != '--six-vs-seven-out-stability'])
    if '--six-vs-seven-planar-orbit' in sys.argv[1:]:
        from scripts.lib.nlsp_six_vs_seven_planar_orbit import main as run
        return run([a for a in sys.argv[1:] if a != '--six-vs-seven-planar-orbit'])
    if '--six-vs-seven-linear' in sys.argv[1:]:
        from scripts.lib.nlsp_six_vs_seven_linear_comparison import main as run
        return run([a for a in sys.argv[1:] if a != '--six-vs-seven-linear'])
    from scripts.lib.nlsp_coupled_sensitivity_pilot import main as run
    return run()

if __name__=='__main__':main()
