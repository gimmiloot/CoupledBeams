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

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_eb_validity_fixed_epsilon_geometry_scan as scan,
)
from scripts.lib import variable_length_timoshenko as TIMO  # noqa: E402


def synthetic_result(
    *,
    x: np.ndarray,
    u: np.ndarray,
    w: np.ndarray,
    u_prime: np.ndarray,
    w_prime: np.ndarray,
    w_second: np.ndarray,
    w_third: np.ndarray,
    epsilon: float = 0.02,
    Lambda: float = np.pi,
    scale: float = 1.0,
) -> scan.ModeResult:
    zeros = np.zeros_like(x)
    rod1 = scan.RodFields(
        x=x,
        u=scale * u,
        w=scale * w,
        u_prime=scale * u_prime,
        w_prime=scale * w_prime,
        w_second=scale * w_second,
        w_third=scale * w_third,
        psi=None,
        psi_prime=None,
        gamma=None,
    )
    rod2 = scan.RodFields(
        x=-x,
        u=zeros,
        w=zeros,
        u_prime=zeros,
        w_prime=zeros,
        w_second=zeros,
        w_third=zeros,
        psi=None,
        psi_prime=None,
        gamma=None,
    )
    section = TIMO.section_from_epsilon_tau(epsilon, 1.0)
    U_axial = 0.5 * TIMO.E * section.area * scan.trapz(rod1.u_prime**2, np.abs(rod1.x))
    U_bending = 0.5 * section.bending_stiffness * scan.trapz(np.asarray(rod1.w_second) ** 2, np.abs(rod1.x))
    return scan.ModeResult(
        model=scan.MODEL_EB,
        epsilon=epsilon,
        beta_deg=0.0,
        mu=0.0,
        eta=0.0,
        sorted_index=1,
        Lambda=Lambda,
        coeff=np.ones(6),
        rod1=rod1,
        rod2=rod2,
        energy=scan.energy_from_terms(scan.MODEL_EB, U_axial, U_bending, 0.0),
        warnings=(),
    )


class FixedEpsilonScanTest(unittest.TestCase):
    def test_theta_max_eb_uses_eb_lambda(self) -> None:
        theta = scan.theta_max_eb(lambda_eb=2.0, epsilon=0.02, mu=0.0, eta=0.0)
        expected = scan.theta_factor() * (0.02 * 2.0) ** 2
        not_timo = scan.theta_factor() * (0.02 * 5.0) ** 2
        self.assertAlmostEqual(theta, expected)
        self.assertNotAlmostEqual(theta, not_timo)

    def test_chi_eff_uses_eb_bending_energy_weights(self) -> None:
        value = scan.chi_eff_from_eb_bending_energies(chi1=1.0, chi2=3.0, bending1=3.0, bending2=1.0)
        self.assertAlmostEqual(value, np.sqrt((3.0 * 1.0**2 + 1.0 * 3.0**2) / 4.0))

    def test_pi_eb_is_invariant_under_mode_scaling(self) -> None:
        x = np.linspace(0.0, 1.0, 801)
        lam = np.pi
        base = synthetic_result(
            x=x,
            u=np.zeros_like(x),
            w=np.sin(lam * x),
            u_prime=np.zeros_like(x),
            w_prime=lam * np.cos(lam * x),
            w_second=-(lam**2) * np.sin(lam * x),
            w_third=-(lam**3) * np.cos(lam * x),
            Lambda=lam,
            scale=1.0,
        )
        scaled = synthetic_result(
            x=x,
            u=np.zeros_like(x),
            w=np.sin(lam * x),
            u_prime=np.zeros_like(x),
            w_prime=lam * np.cos(lam * x),
            w_second=-(lam**2) * np.sin(lam * x),
            w_third=-(lam**3) * np.cos(lam * x),
            Lambda=lam,
            scale=7.0,
        )
        for key in ("Pi_shear_EB", "Pi_rotary_EB", "Pi_EB"):
            self.assertAlmostEqual(scan.pi_eb_metrics(base)[key], scan.pi_eb_metrics(scaled)[key], delta=1.0e-12)

    def test_pure_axial_synthetic_mode_has_near_zero_pi_eb(self) -> None:
        x = np.linspace(0.0, 1.0, 401)
        result = synthetic_result(
            x=x,
            u=np.sin(np.pi * x),
            w=np.zeros_like(x),
            u_prime=np.pi * np.cos(np.pi * x),
            w_prime=np.zeros_like(x),
            w_second=np.zeros_like(x),
            w_third=np.zeros_like(x),
        )
        self.assertAlmostEqual(scan.pi_eb_metrics(result)["Pi_EB"], 0.0, delta=1.0e-14)

    def test_sinusoidal_bending_pi_matches_theta_scaling(self) -> None:
        x = np.linspace(0.0, 1.0, 2001)
        lam = np.pi
        result = synthetic_result(
            x=x,
            u=np.zeros_like(x),
            w=np.sin(lam * x),
            u_prime=np.zeros_like(x),
            w_prime=lam * np.cos(lam * x),
            w_second=-(lam**2) * np.sin(lam * x),
            w_third=-(lam**3) * np.cos(lam * x),
            Lambda=lam,
        )
        self.assertAlmostEqual(scan.pi_eb_metrics(result)["Pi_EB"], scan.theta_max_eb(lam, 0.02, 0.0, 0.0), delta=1.0e-5)

    def test_delta_f_uses_squared_lambda(self) -> None:
        self.assertAlmostEqual(scan.frequency_metrics(2.0, 1.0)["delta_f"], 3.0)

    def test_zero_and_one_percent_false_safe_threshold_logic(self) -> None:
        rows = [
            {"p": 1.0, "delta_f": 0.02},
            {"p": 2.0, "delta_f": 0.08},
            {"p": 3.0, "delta_f": 0.12},
            {"p": 4.0, "delta_f": 0.14},
        ]
        self.assertEqual(scan.conservative_threshold(rows, "p", 0.0), 2.0)
        self.assertEqual(scan.conservative_threshold(rows, "p", 0.01), 2.0)

    def test_mass_weighted_mac_is_sign_and_normalization_invariant(self) -> None:
        x = np.linspace(0.0, 1.0, 101)
        base = synthetic_result(
            x=x,
            u=np.sin(np.pi * x),
            w=np.cos(np.pi * x),
            u_prime=np.pi * np.cos(np.pi * x),
            w_prime=-np.pi * np.sin(np.pi * x),
            w_second=-(np.pi**2) * np.cos(np.pi * x),
            w_third=(np.pi**3) * np.sin(np.pi * x),
        )
        flipped = synthetic_result(
            x=x,
            u=np.sin(np.pi * x),
            w=np.cos(np.pi * x),
            u_prime=np.pi * np.cos(np.pi * x),
            w_prime=-np.pi * np.sin(np.pi * x),
            w_second=-(np.pi**2) * np.cos(np.pi * x),
            w_third=(np.pi**3) * np.sin(np.pi * x),
            scale=-4.0,
        )
        self.assertAlmostEqual(scan.mac_value(scan.weighted_uw_vector(base), scan.weighted_uw_vector(flipped)), 1.0)

    def test_subspace_mac_for_rotated_basis(self) -> None:
        e1 = np.array([1.0, 0.0, 0.0])
        e2 = np.array([0.0, 1.0, 0.0])
        angle = 0.37
        r1 = np.cos(angle) * e1 + np.sin(angle) * e2
        r2 = -np.sin(angle) * e1 + np.cos(angle) * e2
        self.assertAlmostEqual(scan.subspace_mac([e1, e2], [r1, r2]), 1.0)

    def test_candidate_boundary_status(self) -> None:
        status = scan.matching_status(
            mac_uw=0.99,
            mac_margin=0.2,
            candidate_boundary=True,
            root_warning=False,
            cluster_member=False,
            subspace_value=float("nan"),
        )
        self.assertEqual(status, "candidate_boundary")

    def test_smoke_scan_produces_expected_rows_and_finite_predictors(self) -> None:
        out_dir = ROOT / "results" / "_smoke" / "eb_validity_fixed_epsilon_geometry_scan_unittest"
        result = scan.main(
            [
                "--smoke",
                "--beta-values",
                "45",
                "--mu-values",
                "0",
                "--eta-values",
                "0",
                "--n-reported-modes",
                "2",
                "--n-candidate-roots",
                "4",
                "--n-shape-points",
                "101",
                "--output-dir",
                str(out_dir),
            ]
        )
        mode_rows = result["mode_rows"]
        point_rows = result["point_rows"]
        self.assertEqual(len(point_rows), 1)
        self.assertEqual(sum(1 for row in mode_rows if row["comparison_type"] == "sorted_index"), 2)
        self.assertEqual(sum(1 for row in mode_rows if row["comparison_type"] == "homologous_mode"), 2)
        self.assertEqual(point_rows[0]["root_warning_count"], 0)
        for row in mode_rows:
            if row["comparison_type"] == "homologous_mode":
                self.assertTrue(np.isfinite(float(row["Theta_max_EB"])))
                self.assertTrue(np.isfinite(float(row["Pi_EB"])))


if __name__ == "__main__":
    unittest.main()
