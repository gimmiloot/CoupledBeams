import dataclasses
import inspect
import json
import math
import os
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

from my_project.analytic.common.thickness_scaling import (  # noqa: E402
    mass_preserving_denominator_squared,
    mass_preserving_tau_values_from_denominator,
)
from my_project.analytic import formulas_thickness_mismatch as EB  # noqa: E402
from scripts.lib import variable_length_timoshenko as TIMO  # noqa: E402


MU_VALUES = (-0.9, -0.7, -0.3, 0.0, 0.3, 0.7, 0.9)
ETA_VALUES = (-0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9)
SNAPSHOT_ROOT = os.environ.get("COUPLEDBEAMS_THICKNESS_SNAPSHOT_DIR")


class CommonThicknessScalingTest(unittest.TestCase):
    def test_canonical_algebra_and_mass_conservation(self) -> None:
        for mu in MU_VALUES:
            for eta in ETA_VALUES:
                with self.subTest(mu=mu, eta=eta):
                    denominator_squared = mass_preserving_denominator_squared(mu, eta)
                    self.assertEqual(1.0 + 2.0 * mu * eta + eta**2, denominator_squared)
                    denominator = math.sqrt(denominator_squared)
                    tau1, tau2 = mass_preserving_tau_values_from_denominator(eta, denominator)
                    self.assertEqual((1.0 - eta) / denominator, tau1)
                    self.assertEqual((1.0 + eta) / denominator, tau2)
                    mass_factor = (1.0 - mu) * tau1**2 + (1.0 + mu) * tau2**2
                    self.assertAlmostEqual(2.0, mass_factor, delta=4.0e-15)

    def test_eta_zero_limit(self) -> None:
        for mu in MU_VALUES:
            with self.subTest(mu=mu):
                denominator = math.sqrt(mass_preserving_denominator_squared(mu, 0.0))
                self.assertEqual((1.0, 1.0), mass_preserving_tau_values_from_denominator(0.0, denominator))

    def test_swap_symmetry(self) -> None:
        for mu in MU_VALUES:
            for eta in ETA_VALUES:
                with self.subTest(mu=mu, eta=eta):
                    denominator = math.sqrt(mass_preserving_denominator_squared(mu, eta))
                    swapped_denominator = math.sqrt(mass_preserving_denominator_squared(-mu, -eta))
                    tau1, tau2 = mass_preserving_tau_values_from_denominator(eta, denominator)
                    swapped_tau1, swapped_tau2 = mass_preserving_tau_values_from_denominator(
                        -eta, swapped_denominator
                    )
                    self.assertEqual(tau1, swapped_tau2)
                    self.assertEqual(tau2, swapped_tau1)

    def test_public_wrappers_are_exactly_equivalent_on_full_grid(self) -> None:
        for mu in MU_VALUES:
            for eta in ETA_VALUES:
                with self.subTest(mu=mu, eta=eta):
                    eb = EB.thickness_mismatch_factors(mu, eta)
                    timo = TIMO.tau_factors(mu, eta)
                    self.assertEqual(eb.denom, timo.denom)
                    self.assertEqual(eb.tau1, timo.tau1)
                    self.assertEqual(eb.tau2, timo.tau2)
                    self.assertEqual(eb.mass_factor, timo.mass_factor)

    def test_public_api_contract(self) -> None:
        eb = EB.thickness_mismatch_factors(0.3, 0.1)
        timo = TIMO.tau_factors(0.3, 0.1)
        self.assertIs(type(eb), EB.ThicknessMismatchFactors)
        self.assertIs(type(timo), TIMO.TauFactors)
        expected_fields = ("mu", "eta", "denom", "tau1", "tau2")
        self.assertEqual(expected_fields, tuple(field.name for field in dataclasses.fields(type(eb))))
        self.assertEqual(expected_fields, tuple(field.name for field in dataclasses.fields(type(timo))))
        self.assertTrue(EB.ThicknessMismatchFactors.__dataclass_params__.frozen)
        self.assertTrue(TIMO.TauFactors.__dataclass_params__.frozen)
        self.assertEqual(
            "ThicknessMismatchFactors(mu=0.3, eta=0.1, denom=1.03440804327886, "
            "tau1=0.8700628401410972, tau2=1.06341013795023)",
            repr(eb),
        )
        self.assertEqual(
            "TauFactors(mu=0.3, eta=0.1, denom=1.03440804327886, "
            "tau1=0.8700628401410972, tau2=1.06341013795023)",
            repr(timo),
        )
        self.assertEqual(
            "(mu: 'float', eta: 'float') -> 'ThicknessMismatchFactors'",
            str(inspect.signature(EB.thickness_mismatch_factors)),
        )
        self.assertEqual(
            "(mu: 'float', eta: 'float') -> 'TauFactors'",
            str(inspect.signature(TIMO.tau_factors)),
        )
        self.assertEqual(
            [
                "ThicknessMismatchFactors",
                "assemble_clamped_coupled_matrix_eta",
                "assemble_clamped_coupled_matrix_eta_stable",
                "det_eta",
                "find_first_n_roots_eta",
                "find_roots_scan_bisect_eta",
                "local_epsilons",
                "thickness_mismatch_factors",
                "thickness_to_length_ratios",
                "thin_rod_validity",
            ],
            EB.__all__,
        )
        self.assertFalse(hasattr(TIMO, "__all__"))
        self.assertFalse(hasattr(TIMO, "mass_preserving_denominator_squared"))
        self.assertFalse(hasattr(TIMO, "mass_preserving_tau_values_from_denominator"))

    def test_public_exception_contract_and_validation_order(self) -> None:
        invalid_values = (-1.0, 1.0, -1.1, 1.1, float("nan"), float("inf"), float("-inf"))
        for function in (EB.thickness_mismatch_factors, TIMO.tau_factors):
            for value in invalid_values:
                with self.subTest(function=function.__name__, parameter="mu", value=value):
                    with self.assertRaisesRegex(
                        ValueError,
                        r"^mu must lie inside \(-1, 1\) for positive segment lengths\.$",
                    ):
                        function(value, 0.0)
                with self.subTest(function=function.__name__, parameter="eta", value=value):
                    with self.assertRaisesRegex(
                        ValueError,
                        r"^eta must lie inside \(-1, 1\) for positive radii\.$",
                    ):
                        function(0.0, value)
            with self.subTest(function=function.__name__, validation_order=True):
                with self.assertRaisesRegex(
                    ValueError,
                    r"^mu must lie inside \(-1, 1\) for positive segment lengths\.$",
                ):
                    function(-1.0, -1.0)

    def test_downstream_geometry_helpers(self) -> None:
        for epsilon, mu, eta in (
            (0.0025, 0.0, 0.0),
            (0.0025, 0.3, 0.1),
            (0.01, -0.3, -0.1),
            (0.02, 0.7, 0.5),
            (0.05, -0.7, -0.5),
        ):
            with self.subTest(epsilon=epsilon, mu=mu, eta=eta):
                factors = EB.thickness_mismatch_factors(mu, eta)
                self.assertEqual(
                    (epsilon * factors.tau1, epsilon * factors.tau2),
                    EB.local_epsilons(epsilon, mu, eta),
                )
                expected_ratios = (
                    4.0 * epsilon * factors.tau1 / (1.0 - mu),
                    4.0 * epsilon * factors.tau2 / (1.0 + mu),
                )
                self.assertEqual(expected_ratios, EB.thickness_to_length_ratios(epsilon, mu, eta))
                self.assertEqual(
                    max(expected_ratios) <= 0.1,
                    EB.thin_rod_validity(epsilon, mu, eta),
                )
                for tau in (factors.tau1, factors.tau2):
                    section = TIMO.section_from_epsilon_tau(epsilon, tau)
                    radius = 2.0 * epsilon * TIMO.L_SEGMENT * tau
                    self.assertEqual(radius, section.radius)
                    self.assertEqual(math.pi * radius**2, section.area)
                    self.assertEqual(math.pi * radius**4 / 4.0, section.inertia)

    @unittest.skipUnless(SNAPSHOT_ROOT, "external pre-change snapshot path was not supplied")
    def test_public_factors_match_prechange_hex_snapshot(self) -> None:
        snapshot_path = Path(SNAPSHOT_ROOT) / "prechange_factor_outputs.json"
        payload = json.loads(snapshot_path.read_text(encoding="utf-8"))
        functions = {
            (EB.__name__, EB.thickness_mismatch_factors.__name__): EB.thickness_mismatch_factors,
            (TIMO.__name__, TIMO.tau_factors.__name__): TIMO.tau_factors,
        }
        for record in payload["factor_cases"]:
            function = functions[(record["module"], record["function"])]
            result = function(record["mu"], record["eta"])
            with self.subTest(module=record["module"], mu=record["mu"], eta=record["eta"]):
                self.assertEqual(record["result_type_module"], type(result).__module__)
                self.assertEqual(record["result_type_name"], type(result).__name__)
                self.assertEqual(record["repr"], repr(result))
                self.assertEqual(record["hex"]["denom"], float(result.denom).hex())
                self.assertEqual(record["hex"]["tau1"], float(result.tau1).hex())
                self.assertEqual(record["hex"]["tau2"], float(result.tau2).hex())
                self.assertEqual(record["hex"]["mass_factor"], float(result.mass_factor).hex())

    @unittest.skipUnless(SNAPSHOT_ROOT, "external pre-change snapshot path was not supplied")
    def test_fixed_lambda_matrices_match_prechange_snapshot_exactly(self) -> None:
        snapshot_path = Path(SNAPSHOT_ROOT) / "prechange_matrix_outputs.npz"
        with np.load(snapshot_path, allow_pickle=False) as snapshot:
            metadata = json.loads(str(snapshot["metadata_json"].item()))
            for record in metadata:
                if record["model"] == "EB":
                    matrix = EB.assemble_clamped_coupled_matrix_eta(
                        record["Lambda"],
                        record["beta_argument"],
                        record["mu"],
                        record["epsilon"],
                        record["eta"],
                    )
                else:
                    matrix, warnings = TIMO.timo_coupling_matrix(
                        record["Lambda"],
                        record["beta_argument"],
                        record["mu"],
                        record["epsilon"],
                        record["eta"],
                    )
                    self.assertEqual(tuple(record["warnings"]), warnings)
                expected = snapshot[record["key"]]
                with self.subTest(model=record["model"], name=record["name"]):
                    self.assertEqual(tuple(record["shape"]), matrix.shape)
                    self.assertEqual(record["dtype"], str(matrix.dtype))
                    self.assertTrue(np.array_equal(expected, matrix))


if __name__ == "__main__":
    unittest.main()
