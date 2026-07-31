import math
import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
    BookMaterial,
    Geometry,
    assign_modes_by_mac,
    boundary_quality,
    cantilever_boundary_matrix,
    cantilever_boundary_residuals,
    cantilever_clamp_matrix,
    cantilever_energy_fractions,
    cantilever_formulation_mapping_residual,
    cantilever_geometry,
    cantilever_mode_shape,
    cantilever_state_trajectory,
    cantilever_zero_frequency_nullity,
    continue_loss_root,
    decoupled_boundary_factors,
    decoupled_cantilever_boundary_factors,
    eliminated_coefficients,
    find_elastic_roots,
    fixed_free_torsion_omega,
    formulation_mapping_residual,
    free_free_torsion_omega,
    generalized_torsional_stiffness,
    hms_dx_209_material,
    make_rod_point,
    modal_loss_factors,
    partial_bending_boundary_matrix,
    rigid_body_nullity,
    rotate_material,
    solve_complex_root,
    with_gxz_scale,
)


class YartsevChapter2SingleRodTest(unittest.TestCase):
    @staticmethod
    def _cantilever_point(theta: float, length: float = 0.2):
        return make_rod_point(
            theta,
            geometry=cantilever_geometry(length),
            material=hms_dx_209_material(),
        )

    @staticmethod
    def _cantilever_builder(clamp):
        return lambda omega, point, formulation: cantilever_boundary_matrix(
            omega, point, clamp, formulation
        )

    def test_geometry_values_and_dimensions(self) -> None:
        geometry = Geometry()
        self.assertAlmostEqual(geometry.area, 2.463424e-4, delta=1e-16)
        self.assertAlmostEqual(geometry.I_y, 1.955498816853333e-9, delta=1e-21)
        self.assertAlmostEqual(geometry.I_p, 1.503335699370666e-8, delta=1e-20)
        self.assertGreater(geometry.area, 0.0)  # m^2
        self.assertGreater(geometry.I_y, 0.0)  # m^4
        self.assertGreater(geometry.I_p, 0.0)  # m^4

    def test_material_rotation_reciprocity_and_positive_definiteness(self) -> None:
        material = BookMaterial()
        endpoint_0 = rotate_material(0.0)
        endpoint_90 = rotate_material(90.0)
        self.assertAlmostEqual(abs(endpoint_0.Sbar16), 0.0, delta=1e-24)
        self.assertAlmostEqual(abs(endpoint_90.Sbar16), 0.0, delta=1e-24)

        for theta in (0.0, 15.0, 37.0, 60.0, 90.0):
            with self.subTest(theta=theta):
                properties = rotate_material(theta)
                self.assertAlmostEqual(
                    properties.mu_x_xy,
                    properties.Gxy * properties.Sbar16,
                    delta=1e-15,
                )
                self.assertAlmostEqual(
                    properties.mu_xy_x,
                    properties.Ex * properties.Sbar16,
                    delta=1e-15,
                )
                in_plane = np.array(
                    [
                        [properties.Sbar11.real, properties.Sbar16.real],
                        [properties.Sbar16.real, properties.Sbar66.real],
                    ]
                )
                self.assertTrue(np.all(np.linalg.eigvalsh(in_plane) > 0.0))
                self.assertGreater(properties.Sbar55.real, 0.0)

        natural = np.array(
            [
                [1.0 / material.E1_real, -material.nu12 / material.E1_real, 0.0],
                [-material.nu12 / material.E1_real, 1.0 / material.E2_real, 0.0],
                [0.0, 0.0, 1.0 / material.G12_real],
            ]
        )
        self.assertTrue(np.all(np.linalg.eigvalsh(natural) > 0.0))

    def test_torsional_series_convergence_and_admissibility(self) -> None:
        properties = rotate_material(45.0)
        normal = generalized_torsional_stiffness(
            properties, relative_tolerance=1e-12
        )
        strict = generalized_torsional_stiffness(
            properties, relative_tolerance=1e-14
        )
        self.assertGreater(normal.terms_used, 1)
        self.assertGreater(strict.terms_used, normal.terms_used)
        self.assertLess(
            abs(normal.Cbar - strict.Cbar) / abs(strict.Cbar), 2e-12
        )
        self.assertLess(abs(normal.C_T - strict.C_T) / abs(strict.C_T), 2e-12)
        self.assertTrue(np.isfinite(normal.Cbar))
        self.assertTrue(np.isfinite(normal.C_T))
        self.assertGreater(normal.Cbar.real, 0.0)  # N m^2
        self.assertGreater(normal.C_T.real, 0.0)  # N m^2
        self.assertLessEqual(normal.C_T.real, normal.Cbar.real)

    def test_elastic_limit_is_real_and_has_zero_modal_loss(self) -> None:
        point = make_rod_point(33.0, material_mode="elastic")
        roots = find_elastic_roots(
            point, "state_corrected", num_roots=3, scan_step_hz=12.0
        )
        for item in roots:
            self.assertAlmostEqual(item.omega.imag, 0.0, delta=1e-12)
            exact, approximate, difference = modal_loss_factors(item.omega)
            self.assertEqual(exact, 0.0)
            self.assertEqual(approximate, 0.0)
            self.assertEqual(difference, 0.0)
            quality = boundary_quality(item.omega, point, "state_corrected")
            self.assertLess(abs(quality.determinant.imag), 1e-10)

    def test_book_complex_convention_and_physical_decay_sign(self) -> None:
        material = BookMaterial()
        E1, E2, G12, G13, G23 = material.moduli("book_complex")
        self.assertEqual(E1, material.E1_real * (1.0 + 1j * material.eta1))
        self.assertEqual(E2, material.E2_real * (1.0 + 1j * material.eta2))
        self.assertEqual(G12, material.G12_real * (1.0 + 1j * material.eta12))
        self.assertEqual(G13, material.G13_real * (1.0 + 1j * material.eta13))
        self.assertEqual(G23, material.G23_real * (1.0 + 1j * material.eta23))

        elastic = make_rod_point(45.0)
        elastic_root = find_elastic_roots(
            elastic, "state_corrected", num_roots=1
        )[0]
        result = continue_loss_root(
            lambda scale: make_rod_point(
                45.0, material_mode="book_complex", loss_scale=scale
            ),
            "state_corrected",
            elastic_root.omega.real,
        )
        self.assertGreater(result.omega.real, 0.0)
        self.assertGreater(result.omega.imag, 0.0)
        eta_exact, eta_approx, relative_difference = modal_loss_factors(result.omega)
        self.assertGreater(eta_exact, 0.0)
        self.assertGreater(eta_approx, 0.0)
        self.assertLess(relative_difference, 2e-4)

    def test_orthotropic_decoupling_factorization_and_torsion_roots(self) -> None:
        for theta in (0.0, 90.0):
            point = make_rod_point(theta)
            with self.subTest(theta=theta):
                self.assertLess(abs(point.properties.Sbar16), 1e-24)
                omega = 2.0 * math.pi * 1234.0
                full, bending, torsion = decoupled_boundary_factors(omega, point)
                self.assertLess(
                    abs(full - bending * torsion) / max(abs(full), 1.0), 2e-13
                )
                for mode_number in (1, 2, 3):
                    torsion_omega = free_free_torsion_omega(point, mode_number)
                    full_factor, _, torsion_factor = decoupled_boundary_factors(
                        torsion_omega, point
                    )
                    self.assertLess(abs(torsion_factor), 2e-12)
                    _, bending_factor, _ = decoupled_boundary_factors(
                        torsion_omega, point
                    )
                    self.assertLess(
                        abs(full_factor) / max(abs(bending_factor), 1.0), 2e-12
                    )
                    wave_number = torsion_omega * np.sqrt(
                        point.material.rho * point.geometry.I_p
                        / point.torsion.C_T
                    )
                    self.assertAlmostEqual(
                        wave_number.real * point.geometry.length,
                        mode_number * math.pi,
                        delta=2e-12,
                    )

    def test_state_and_corrected_eliminated_formulations_agree(self) -> None:
        for theta in (0.0, 37.0, 90.0):
            point = make_rod_point(theta)
            state_roots = find_elastic_roots(
                point, "state_corrected", num_roots=8, scan_step_hz=10.0
            )
            eliminated_roots = find_elastic_roots(
                point, "eliminated_corrected", num_roots=8, scan_step_hz=10.0
            )
            state_frequencies = np.array([item.frequency_hz for item in state_roots])
            eliminated_frequencies = np.array(
                [item.frequency_hz for item in eliminated_roots]
            )
            with self.subTest(theta=theta):
                np.testing.assert_allclose(
                    state_frequencies,
                    eliminated_frequencies,
                    rtol=2e-9,
                    atol=2e-6,
                )
                self.assertLess(
                    formulation_mapping_residual(
                        state_roots[3].omega, point, samples=9
                    ),
                    2e-8,
                )
                for state_item, eliminated_item in zip(
                    state_roots, eliminated_roots
                ):
                    self.assertLess(state_item.relative_singular_residual, 1e-8)
                    self.assertLess(eliminated_item.relative_singular_residual, 1e-8)

    def test_printed_signs_remain_a_distinct_diagnostic(self) -> None:
        point = make_rod_point(45.0)
        corrected = eliminated_coefficients(2.0 * math.pi * 1000.0, point, printed=False)
        printed = eliminated_coefficients(2.0 * math.pi * 1000.0, point, printed=True)
        self.assertEqual(printed["d0"], -corrected["d0"])
        self.assertEqual(printed["f0"], -corrected["f0"])
        for name in ("a0", "b0", "c0", "e0", "ag", "bg", "cg"):
            self.assertEqual(printed[name], corrected[name])

    def test_three_rigid_modes_are_not_positive_tones(self) -> None:
        point = make_rod_point(45.0)
        self.assertEqual(rigid_body_nullity(point), 3)
        roots = find_elastic_roots(
            point, "state_corrected", num_roots=8, scan_step_hz=10.0
        )
        self.assertEqual(len(roots), 8)
        self.assertTrue(all(item.omega.real > 0.0 for item in roots))
        self.assertGreater(roots[0].frequency_hz, 100.0)

    def test_root_scan_and_series_convergence(self) -> None:
        point = make_rod_point(45.0)
        coarse = find_elastic_roots(
            point, "state_corrected", num_roots=8, scan_step_hz=10.0
        )
        fine = find_elastic_roots(
            point, "state_corrected", num_roots=8, scan_step_hz=5.0
        )
        np.testing.assert_allclose(
            [item.frequency_hz for item in coarse],
            [item.frequency_hz for item in fine],
            rtol=2e-10,
            atol=2e-6,
        )

        strict_point = make_rod_point(
            45.0, series_relative_tolerance=1e-14
        )
        strict = find_elastic_roots(
            strict_point, "state_corrected", num_roots=8, scan_step_hz=10.0
        )
        np.testing.assert_allclose(
            [item.frequency_hz for item in coarse],
            [item.frequency_hz for item in strict],
            rtol=2e-10,
            atol=2e-6,
        )

    def test_free_free_representative_root_regression(self) -> None:
        roots = find_elastic_roots(
            make_rod_point(45.0),
            "state_corrected",
            num_roots=3,
            scan_step_hz=10.0,
        )
        np.testing.assert_allclose(
            [item.frequency_hz for item in roots],
            [309.321329443630, 844.113322677715, 1629.785138186661],
            rtol=2e-11,
            atol=2e-7,
        )

    def test_cantilever_geometry_and_hms_dx_209_table_values(self) -> None:
        geometry = cantilever_geometry()
        material = hms_dx_209_material()
        self.assertEqual((geometry.a, geometry.b, geometry.length), (0.005, 0.02, 0.2))
        self.assertAlmostEqual(geometry.area, 1.0e-4, delta=1e-18)
        self.assertAlmostEqual(geometry.I_y, 2.0833333333333336e-10, delta=1e-23)
        self.assertAlmostEqual(geometry.I_p, 3.541666666666667e-9, delta=1e-22)
        self.assertEqual(
            (
                material.E1_real,
                material.E2_real,
                material.G12_real,
                material.G13_real,
                material.G23_real,
                material.nu12,
                material.rho,
            ),
            (191e9, 5e9, 3e9, 3e9, 2.5e9, 0.279, 1580.0),
        )
        self.assertEqual(
            (material.eta1, material.eta2, material.eta12, material.eta13, material.eta23),
            (7.8e-4, 6.7e-3, 1.16e-2, 1.16e-2, 1.15e-2),
        )
        self.assertEqual(material.ply_thickness, 2e-4)

    def test_cantilever_boundary_matrices_and_zero_frequency(self) -> None:
        point = self._cantilever_point(15.0)
        shear_stiffness = (
            point.geometry.shear_factor
            * point.properties.Gxz
            * point.geometry.area
        )
        expected_book = np.array(
            [
                [0.0, 0.0, 0.0],
                [-1.0 / shear_stiffness, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ]
        )
        expected_section = expected_book.copy()
        expected_section[1, 0] = 0.0
        np.testing.assert_allclose(
            cantilever_clamp_matrix(point, "book_slope_clamp"), expected_book
        )
        np.testing.assert_allclose(
            cantilever_clamp_matrix(point, "timoshenko_section_clamp"),
            expected_section,
        )
        for clamp in ("book_slope_clamp", "timoshenko_section_clamp"):
            matrix = cantilever_boundary_matrix(0.0, point, clamp)
            self.assertEqual(matrix.shape, (3, 3))
            self.assertNotEqual(np.linalg.det(matrix), 0.0)
            self.assertEqual(cantilever_zero_frequency_nullity(point, clamp), 0)

    def test_cantilever_boundary_residuals(self) -> None:
        point = self._cantilever_point(15.0)
        for clamp in ("book_slope_clamp", "timoshenko_section_clamp"):
            roots = find_elastic_roots(
                point,
                "state_corrected",
                num_roots=6,
                boundary_matrix_builder=self._cantilever_builder(clamp),
            )
            for item in roots:
                residuals = cantilever_boundary_residuals(item.omega, point, clamp)
                scaled = cantilever_state_trajectory(item.omega, point, clamp, [0.0, 1.0])
                self.assertLess(abs(scaled[0, 0]), 1e-12)
                self.assertLess(abs(scaled[0, 2]), 1e-12)
                self.assertLess(max(abs(scaled[1, 3:])), 2e-8)
                if clamp == "book_slope_clamp":
                    self.assertLess(abs(residuals["slope_compensation"]), 2e-12)
                else:
                    self.assertLess(abs(residuals["psi_b_0"]), 2e-12)

    def test_cantilever_state_eliminated_equivalence(self) -> None:
        builder = self._cantilever_builder("book_slope_clamp")
        for theta in (0.0, 7.0, 15.0, 19.0, 45.0, 90.0):
            point = self._cantilever_point(theta)
            state = find_elastic_roots(
                point,
                "state_corrected",
                num_roots=6,
                boundary_matrix_builder=builder,
            )
            eliminated = find_elastic_roots(
                point,
                "eliminated_corrected",
                num_roots=6,
                boundary_matrix_builder=builder,
            )
            with self.subTest(theta=theta):
                np.testing.assert_allclose(
                    [item.frequency_hz for item in state],
                    [item.frequency_hz for item in eliminated],
                    rtol=3e-9,
                    atol=3e-6,
                )
                for state_item, eliminated_item in zip(state, eliminated):
                    self.assertLess(state_item.relative_singular_residual, 2e-8)
                    self.assertLess(eliminated_item.relative_singular_residual, 2e-8)
                    self.assertLess(
                        cantilever_formulation_mapping_residual(
                            state_item.omega, point, samples=7
                        ),
                        3e-8,
                    )

    def test_cantilever_orthotropic_decoupling_torsion_and_energy(self) -> None:
        for theta in (0.0, 90.0):
            point = self._cantilever_point(theta)
            self.assertLess(abs(point.properties.Sbar16), 1e-24)
            for clamp in ("book_slope_clamp", "timoshenko_section_clamp"):
                for mode in (1, 2, 3):
                    omega = fixed_free_torsion_omega(point, mode)
                    full, bending, torsion = decoupled_cantilever_boundary_factors(
                        omega, point, clamp
                    )
                    self.assertLess(abs(full - bending * torsion) / max(abs(full), 1.0), 2e-13)
                    self.assertLess(abs(torsion), 3e-12)
                    wave_number = omega * np.sqrt(
                        point.material.rho * point.geometry.I_p / point.torsion.C_T
                    )
                    self.assertAlmostEqual(
                        wave_number.real * point.geometry.length,
                        (2 * mode - 1) * math.pi / 2.0,
                        delta=3e-12,
                    )
                bending, shear, torsion = cantilever_energy_fractions(
                    fixed_free_torsion_omega(point, 1), point, clamp
                )
                self.assertLess(bending + shear, 2e-20)
                self.assertAlmostEqual(torsion, 1.0, delta=2e-13)

    def test_cantilever_shear_rigid_limit(self) -> None:
        differences = []
        for scale in (1.0e2, 1.0e4, 1.0e6):
            point = with_gxz_scale(self._cantilever_point(15.0), scale)
            book = find_elastic_roots(
                point,
                "state_corrected",
                num_roots=3,
                boundary_matrix_builder=self._cantilever_builder("book_slope_clamp"),
            )
            section = find_elastic_roots(
                point,
                "state_corrected",
                num_roots=3,
                boundary_matrix_builder=self._cantilever_builder("timoshenko_section_clamp"),
            )
            book_shapes = [
                cantilever_mode_shape(item.omega, point, "book_slope_clamp", np.linspace(0.0, 1.0, 41))
                for item in book
            ]
            section_shapes = [
                cantilever_mode_shape(item.omega, point, "timoshenko_section_clamp", np.linspace(0.0, 1.0, 41))
                for item in section
            ]
            assignment, _ = assign_modes_by_mac(book_shapes, section_shapes)
            differences.append(
                max(
                    abs(section[assignment[index]].frequency_hz - book[index].frequency_hz)
                    / book[index].frequency_hz
                    for index in range(3)
                )
            )
        self.assertGreater(differences[0], differences[1])
        self.assertGreater(differences[1], differences[2])
        self.assertLess(differences[-1], 2e-8)

    def test_cantilever_root_and_complex_limits(self) -> None:
        point = self._cantilever_point(15.0)
        builder = self._cantilever_builder("book_slope_clamp")
        coarse = find_elastic_roots(
            point,
            "state_corrected",
            num_roots=6,
            scan_step_hz=10.0,
            boundary_matrix_builder=builder,
        )
        fine = find_elastic_roots(
            point,
            "state_corrected",
            num_roots=6,
            scan_step_hz=5.0,
            boundary_matrix_builder=builder,
        )
        strict_point = make_rod_point(
            15.0,
            geometry=cantilever_geometry(),
            material=hms_dx_209_material(),
            series_relative_tolerance=1e-14,
        )
        strict = find_elastic_roots(
            strict_point,
            "state_corrected",
            num_roots=6,
            boundary_matrix_builder=builder,
        )
        np.testing.assert_allclose(
            [item.frequency_hz for item in coarse],
            [item.frequency_hz for item in fine],
            rtol=3e-10,
            atol=3e-6,
        )
        np.testing.assert_allclose(
            [item.frequency_hz for item in coarse],
            [item.frequency_hz for item in strict],
            rtol=3e-10,
            atol=3e-6,
        )
        zero_loss = make_rod_point(
            15.0,
            material_mode="book_complex",
            loss_scale=0.0,
            geometry=cantilever_geometry(),
            material=hms_dx_209_material(),
        )
        zero_root = solve_complex_root(
            zero_loss,
            "state_corrected",
            coarse[0].omega,
            boundary_matrix_builder=builder,
        )
        self.assertAlmostEqual(zero_root.omega.imag, 0.0, delta=2e-8)
        damped = continue_loss_root(
            lambda scale: make_rod_point(
                15.0,
                material_mode="book_complex",
                loss_scale=scale,
                geometry=cantilever_geometry(),
                material=hms_dx_209_material(),
            ),
            "state_corrected",
            coarse[0].omega.real,
            boundary_matrix_builder=builder,
        )
        self.assertGreater(damped.omega.imag, 0.0)
        exact, approximate, difference = modal_loss_factors(damped.omega)
        self.assertGreater(exact, 0.0)
        self.assertGreater(approximate, 0.0)
        self.assertLess(difference, 2e-4)

    def test_partial_bending_matrix_and_mac_assignment(self) -> None:
        point = self._cantilever_point(15.0)
        self.assertEqual(
            partial_bending_boundary_matrix(2.0 * math.pi * 1000.0, point).shape,
            (2, 2),
        )
        basis = [
            np.array([[1.0], [0.0]], dtype=complex),
            np.array([[0.0], [1.0]], dtype=complex),
        ]
        assignment, mac = assign_modes_by_mac(basis, basis[::-1])
        self.assertEqual(assignment, [1, 0])
        np.testing.assert_allclose(mac, [[0.0, 1.0], [1.0, 0.0]])


if __name__ == "__main__":
    unittest.main()
