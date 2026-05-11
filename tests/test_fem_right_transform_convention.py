import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import assemble_clamped_coupled_matrix
from my_project.fem import python_fem as fem
from scripts.analysis.check_analytic_shape_in_fem_residual import (
    analytic_local_fields,
    assemble_fem_matrices,
    build_analytic_global_q,
    build_params_for_case,
    residual_metrics,
    resolve_analytic_branch,
)
from scripts.lib.analytic_coupled_rods_shapes import analytic_null_vector


class FemRightTransformConventionTest(unittest.TestCase):
    def test_rotation_matrix_maps_local_to_article_global(self) -> None:
        beta = np.deg2rad(15.0)
        transform = fem.rotation_matrix_6x6(beta)[:2, :2]
        c, s = np.cos(beta), np.sin(beta)

        axial_global = transform @ np.array([1.0, 0.0])
        transverse_global = transform @ np.array([0.0, 1.0])

        np.testing.assert_allclose(axial_global, np.array([c, s]), rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(transverse_global, np.array([-s, c]), rtol=1e-14, atol=1e-14)

    def test_article_determinant_shape_has_small_corrected_fem_residual(self) -> None:
        beta_deg = 15.0
        epsilon = 0.0025
        params = build_params_for_case(epsilon=epsilon, l_total=2.0)
        for mu in (0.0, 0.2):
            with self.subTest(mu=mu):
                tracking = resolve_analytic_branch(
                    branch_id="bending_desc_05",
                    beta_deg=beta_deg,
                    mu=mu,
                    epsilon=epsilon,
                    num_analytic_roots=20,
                )
                Lambda = float(tracking["point"].Lambda)
                matrix = assemble_clamped_coupled_matrix(Lambda, np.deg2rad(beta_deg), mu, epsilon)
                coeff, _, _ = analytic_null_vector(matrix)
                fields = analytic_local_fields(
                    Lambda,
                    mu=mu,
                    epsilon=epsilon,
                    coeff=coeff,
                    s_norm=np.linspace(0.0, 1.0, fem.N_ELEM + 1),
                    right_theta_sign=1.0,
                )
                q_analytic, _ = build_analytic_global_q(fields, beta_deg=beta_deg)
                k_free, m_free, _, _, free = assemble_fem_matrices(params, mu=mu, beta_deg=beta_deg)
                residual = residual_metrics(k_free, m_free, q_analytic, free, omega_sq=Lambda**4)

                self.assertLess(float(residual["relative_residual"]), 1e-5)


if __name__ == "__main__":
    unittest.main()
