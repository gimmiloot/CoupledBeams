"""Free outer rotation, reusing the verified transverse-support assembly.

Full RL/Timoshenko equations already admit k_end=0 and are used unchanged.
The leading constructor historically requires positive rotational parameters;
reuse its basis/mass/curvature assembly, then omit Ke BEFORE any eigensolve.
No penalty, subtraction of nearly equal matrices, or extra constraint is used.
"""
import numpy as np
from scripts.lib.isotropic_two_springs_external_transverse_spring import Leading


def leading_model(k_joint, k_w, *, mu=.01, order=24, quadrature=64):
    model = Leading(k_joint, 1., mu=mu, order=order,
                    quadrature=quadrature, k_w=k_w, analytic_curvature=True)
    # The positive constructor value only permits creation of the existing
    # space. There is no solution with that value in this research workflow.
    model.k_end = 0.
    model.Ke = np.zeros_like(model.Kb)
    model.K0 = model.Kb + model.Kj
    if k_w != 'RIGID':
        model.K0 = model.K0 + model.Kw
    model.matrix_checks['symmetry'] = max(
        float(np.max(abs(a-a.T))) for a in (model.K0, model.M, model.P))
    return model
