"""Scoped coupled-rod tangent/full nonlinear pilot; never launches 3D FEM."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.lib.nlsp_coupled_sensitivity_pilot import main
if __name__=='__main__':main()
