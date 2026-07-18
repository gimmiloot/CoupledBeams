import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.analysis.thickness_mismatch.audits import audit_timoshenko_shape_construction as shape_audit
from scripts.lib import in_plane_shape_geometry as geometry


def mac(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    denominator = float(np.dot(a, a) * np.dot(b, b))
    return float(np.dot(a, b) ** 2 / denominator)


def display_vectors(result: shape_audit.ModeResult) -> tuple[np.ndarray, np.ndarray]:
    if result.model == shape_audit.MODEL_EB:
        dx1, dy1 = geometry.eb_rod1_local_displacement_to_display(result.rod1.u, result.rod1.w)
        dx2, dy2 = geometry.eb_rod2_local_displacement_to_display(
            result.rod2.u,
            result.rod2.w,
            beta_deg=result.beta_deg,
        )
    else:
        dx1, dy1 = geometry.rod1_local_displacement_to_display(result.rod1.u, result.rod1.w)
        dx2, dy2 = geometry.rod2_local_displacement_to_display(
            result.rod2.u,
            result.rod2.w,
            beta_deg=result.beta_deg,
        )
    rod2 = np.concatenate([dx2, dy2])
    return np.concatenate([dx1, dy1, rod2]), rod2


class InPlaneShapeGeometryTest(unittest.TestCase):
    def test_rod2_display_basis_and_pure_components(self) -> None:
        for beta_deg in (0.0, 15.0, 45.0, 90.0):
            with self.subTest(beta_deg=beta_deg):
                basis = geometry.rod2_display_basis(beta_deg)
                self.assertAlmostEqual(float(np.dot(basis.tangent, basis.normal)), 0.0, places=14)
                self.assertAlmostEqual(float(np.linalg.norm(basis.tangent)), 1.0, places=14)
                self.assertAlmostEqual(float(np.linalg.norm(basis.normal)), 1.0, places=14)

                axial = np.array(geometry.rod2_local_displacement_to_display(
                    np.array([1.0]),
                    np.array([0.0]),
                    beta_deg=beta_deg,
                ))[:, 0]
                transverse = np.array(geometry.rod2_local_displacement_to_display(
                    np.array([0.0]),
                    np.array([1.0]),
                    beta_deg=beta_deg,
                ))[:, 0]
                np.testing.assert_allclose(axial, basis.tangent, rtol=0.0, atol=1.0e-14)
                np.testing.assert_allclose(transverse, basis.normal, rtol=0.0, atol=1.0e-14)
                self.assertAlmostEqual(float(np.dot(transverse, basis.tangent)), 0.0, places=14)

    def test_beta45_transverse_is_not_tangent(self) -> None:
        basis = geometry.rod2_display_basis(45.0)
        self.assertFalse(np.allclose(basis.normal, basis.tangent, rtol=0.0, atol=1.0e-14))
        self.assertAlmostEqual(float(np.dot(basis.normal, basis.tangent)), 0.0, places=14)

    def test_signed_rod2_grid_runs_from_joint_to_clamp(self) -> None:
        l2 = 1.4
        x_joint = 0.6
        x2_grid = np.linspace(-l2, 0.0, 11)
        curve = geometry.rod2_local_fields_to_display(
            x2_grid,
            np.zeros_like(x2_grid),
            np.zeros_like(x2_grid),
            l2=l2,
            x_joint=x_joint,
            beta_deg=45.0,
            scale=1.0,
        )
        self.assertAlmostEqual(float(curve.x_base[0]), x_joint, places=14)
        self.assertAlmostEqual(float(curve.y_base[0]), 0.0, places=14)
        expected_end = np.array([x_joint, 0.0]) + l2 * geometry.rod2_display_basis(45.0).tangent
        np.testing.assert_allclose(
            np.array([curve.x_base[-1], curve.y_base[-1]]),
            expected_end,
            rtol=0.0,
            atol=1.0e-14,
        )

    def test_thin_limit_eb_timoshenko_display_mac(self) -> None:
        beta_deg = 45.0
        epsilon = 0.0025
        eta = 0.0
        provider = shape_audit.RootProvider(beta_deg, eta, 6)
        for mu in (0.0, 0.2, 0.4, 0.6):
            eb_roots, timo_roots, warnings = provider.roots(epsilon, mu)
            for sorted_index in (4, 5, 6):
                with self.subTest(mu=mu, sorted_index=sorted_index):
                    eb = shape_audit.eb_mode_result(
                        epsilon=epsilon,
                        beta_deg=beta_deg,
                        eta=eta,
                        mu=mu,
                        sorted_index=sorted_index,
                        Lambda=float(eb_roots[sorted_index - 1]),
                        n_points=401,
                        case_role="display_geometry_regression",
                        root_source="sorted_global_solver",
                        root_warnings=warnings,
                    )
                    timo = shape_audit.timo_mode_result(
                        epsilon=epsilon,
                        beta_deg=beta_deg,
                        eta=eta,
                        mu=mu,
                        sorted_index=sorted_index,
                        Lambda=float(timo_roots[sorted_index - 1]),
                        n_points=401,
                        case_role="display_geometry_regression",
                        root_source="sorted_global_solver",
                        root_warnings=warnings,
                    )
                    eb_full, eb_rod2 = display_vectors(eb)
                    timo_full, timo_rod2 = display_vectors(timo)
                    self.assertGreater(mac(eb_full, timo_full), 0.98)
                    self.assertGreater(mac(eb_rod2, timo_rod2), 0.98)


if __name__ == "__main__":
    unittest.main()
