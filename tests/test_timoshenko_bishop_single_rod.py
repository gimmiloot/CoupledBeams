"""Exact kinematic audit; a successful test is not acceptance of a hybrid theory."""
from fractions import Fraction as F
import json

import numpy as np
import pytest

from scripts.analysis import audit_timoshenko_bishop_single_rod as a
from scripts.lib import bishop_longitudinal as bishop
from scripts.lib import isotropic_rectangular_timoshenko_coupled_beams as timo


def jet(name, dx=0, dt=0):
    return name, dx, dt


def coefficient(matrix, left, right):
    return matrix.get((left, right), F(0))


def matrices(nu=F(3, 10), full=False, zc=F(0)):
    b, h = F(1, 5), F(1, 20)
    return tuple(a.integrated(m, b, h, zc=zc)
                 for m in a.candidate_hessians(F(1), F(1), nu, b, h, full))


@pytest.mark.parametrize("b,h", [(F(1, 5), F(1, 20)), (F(7, 3), F(2, 7)), (F(1), F(1))])
def test_exact_rectangle_moments_and_parallel_axis_negative_control(b, h):
    assert a.moment(0, 0, b, h) == b*h
    for i, j in ((1, 0), (0, 1), (1, 1)):
        assert a.moment(i, j, b, h) == 0
    assert a.moment(2, 0, b, h) == h*b**3/12  # I_z, NOT I_y
    assert a.moment(0, 2, b, h) == b*h**3/12
    yc, zc = F(2, 9), F(-3, 7)
    assert a.moment(1, 0, b, h, yc, zc) == b*h*yc
    assert a.moment(0, 1, b, h, yc, zc) == b*h*zc
    assert a.moment(1, 1, b, h, yc, zc) == b*h*yc*zc
    assert a.moment(0, 2, b, h, yc, zc) == b*h**3/12+b*h*zc**2


def test_strains_differentiated_from_displacements_with_project_signs():
    nu = F(3, 10)
    actual = a.strains(a.displacement(nu, F(1, 5), F(1, 20)))
    assert actual == (
        {jet("u", 1): {(0, 0): F(1)}, jet("psi", 1): {(0, 1): F(-1)}},
        {jet("u", 1): {(0, 0): -nu}},
        {jet("u", 1): {(0, 0): -nu}}, {},
        {jet("w", 1): {(0, 0): F(1)}, jet("psi"): {(0, 0): F(-1)}, jet("u", 2): {(0, 1): -nu}},
        {jet("u", 2): {(1, 0): -nu}},
    )


@pytest.mark.parametrize("nu", [F(0), F(3, 10), F(-1, 5), F(9, 20)])
@pytest.mark.parametrize("full", [False, True])
def test_all_centered_mixed_coefficients_are_exactly_zero(nu, full):
    # All matrix pairs inspected, not just a preselected list of expected zeros.
    for matrix in matrices(nu, full):
        assert all(value == 0 for key, value in matrix.items() if a.axial_bending(key))
        assert all(value == matrix.get((right, left), F(0))
                   for (left, right), value in matrix.items())


def test_every_minimal_mixed_term_before_centroid_substitution():
    nu = F(3, 10)
    raw_m, raw_k = a.candidate_hessians(F(1), F(1), nu, F(1, 5), F(1, 20))
    cross = lambda matrix: {key: p for key, p in matrix.items() if a.axial_bending(key) and key[0][0] == "u"}
    assert cross(raw_m) == {
        (jet("u", dt=1), jet("psi", dt=1)): {(0, 1): F(-1)},
        (jet("u", 1, 1), jet("w", dt=1)): {(0, 1): -nu},
    }
    G = 1/(2*(1+nu))
    assert cross(raw_k) == {
        (jet("u", 1), jet("psi", 1)): {(0, 1): F(-1)},
        (jet("u", 2), jet("w", 1)): {(0, 1): -nu*G},
        (jet("u", 2), jet("psi")): {(0, 1): nu*G},
    }
    mass, stiffness = matrices(zc=F(1, 100))
    Qz = F(1, 5)*F(1, 20)*F(1, 100)
    assert coefficient(mass, jet("u", dt=1), jet("psi", dt=1)) == -Qz != 0
    assert coefficient(mass, jet("u", 1, 1), jet("w", dt=1)) == -nu*Qz != 0
    assert coefficient(stiffness, jet("u", 1), jet("psi", 1)) == -Qz != 0
    assert coefficient(stiffness, jet("u", 2), jet("w", 1)) == -nu*G*Qz != 0
    assert coefficient(stiffness, jet("u", 2), jet("psi")) == nu*G*Qz != 0


def test_raw_normal_stress_and_self_term_expose_failed_timoshenko_limit():
    E, nu = F(1), F(3, 10)
    C = a.constitutive(E, nu)
    strains = a.strains(a.displacement(nu, F(1, 5), F(1, 20)))
    for row in (1, 2):
        stress = a.add(*(a.scale(strains[i].get(jet("u", 1), {}), C[row][i]) for i in range(6)))
        assert stress == {}  # free lateral normal stress for pure axial
        stress = a.add(*(a.scale(strains[i].get(jet("psi", 1), {}), C[row][i]) for i in range(6)))
        assert stress == {(0, 1): -F(15, 26)}  # -lambda_L*z, not zero
    _, stiffness = matrices()
    Iy, area = F(1, 5)*F(1, 20)**3/12, F(1, 100)
    assert coefficient(stiffness, jet("psi", 1), jet("psi", 1))/(E*Iy) == F(35, 26)
    assert coefficient(stiffness, jet("w", 1), jet("w", 1))/(F(5, 6)*C[4][4]*area) == F(6, 5)
    # At nu=0 the normal mismatch disappears; the shear-correction mismatch does not.
    _, zero = matrices(nu=F(0))
    assert coefficient(zero, jet("psi", 1), jet("psi", 1)) == E*Iy


def test_full_poisson_compatibility_requires_extra_bending_terms():
    nu, b, h = F(3, 10), F(1, 5), F(1, 20)
    fields = a.displacement(nu, b, h, True)
    strain = a.strains(fields)
    for transverse in (strain[1], strain[2]):
        assert transverse == {key: a.scale(p, -nu) for key, p in strain[0].items()}
    assert strain[3] == {}  # no gamma_yz: both transverse warp components required
    assert a.integrate(fields[2][jet("psi", 1)], b, h) == 0  # w is still centroidal
    mass, stiffness = matrices(full=True)
    K4 = b*h*(b**4+5*b*b*h*h+h**4)/720
    assert K4 > 0
    assert coefficient(mass, jet("psi", 1, 1), jet("psi", 1, 1)) == nu**2*K4 > 0
    assert coefficient(stiffness, jet("psi", 2), jet("psi", 2)) == nu**2*K4/(2*(1+nu)) > 0
    assert coefficient(stiffness, jet("psi", 1), jet("psi", 1)) == b*h**3/12


def test_full_poisson_additional_cross_densities_are_derived_not_discarded():
    nu, b, h = F(3, 10), F(1, 5), F(1, 20)
    mass, stiffness = a.candidate_hessians(F(1), F(1), nu, b, h, True)
    c = (h*h-b*b)/12
    # y^2*z + z*r = (y^2*z + z^3 - c*z)/2, an odd cubic.
    odd = {(2, 1): F(1, 2), (0, 3): F(1, 2), (0, 1): -c/2}
    assert mass[(jet("u", 1, 1), jet("psi", 1, 1))] == a.scale(odd, -nu**2)
    assert stiffness[(jet("u", 2), jet("psi", 2))] == a.scale(odd, -nu**2/(2*(1+nu)))
    assert a.integrate(odd, b, h) == 0
    assert a.integrate(odd, b, h, zc=F(1, 100)) != 0


def test_material_factors_are_retained_in_general_energy():
    args = (F(3, 10), F(1, 5), F(1, 20))
    unit = a.candidate_hessians(F(1), F(1), *args)
    scaled = a.candidate_hessians(F(7), F(11), *args)
    for first, second, factor in zip(unit, scaled, (11, 7)):
        assert second == {key: a.scale(p, factor) for key, p in first.items()}


def test_shear_profile_parity_is_required_beyond_kappa_scalar():
    b, h = F(1, 5), F(1, 20)
    # An even xz profile is orthogonal to Bishop's z; an odd xy profile in
    # both y,z is orthogonal to Bishop's y after integration over z.
    even_xz, odd_xy = {(0, 0): F(1), (0, 2): -4/h**2}, {(1, 1): F(1)}
    assert a.integrate(a.multiply(even_xz, {(0, 1): F(1)}), b, h) == 0
    assert a.integrate(a.multiply(odd_xy, {(1, 0): F(1)}), b, h) == 0
    # Same self norm for +/- odd contamination, different cross coefficient:
    # a scalar shear correction cannot by itself specify the cross term.
    profiles = [a.add(even_xz, {(0, 1): F(sign)/h}) for sign in (-1, 1)]
    norms = [a.integrate(a.multiply(p, p), b, h) for p in profiles]
    crosses = [a.integrate(a.multiply(p, {(0, 1): F(1)}), b, h) for p in profiles]
    assert norms[0] == norms[1]
    assert crosses[0] == -crosses[1] != 0


def test_pure_axial_coefficients_match_existing_bishop_and_pure_bending_inertia():
    p = a.PARAMETERS
    b, h, nu = p["b"], p["h"], p["nu"]
    area, Iy, Iz = b*h, b*h**3/12, h*b**3/12
    mass, stiffness = matrices()
    J, H = nu**2*(Iy+Iz), nu**2*(Iy+Iz)/(2*(1+nu))
    s = bishop.Segment(1., float(area), float(area), float(H), float(J))
    for actual, expected in [(coefficient(mass, jet("u", dt=1), jet("u", dt=1)), s.m),
                             (coefficient(mass, jet("u", 1, 1), jet("u", 1, 1)), s.J),
                             (coefficient(stiffness, jet("u", 1), jet("u", 1)), s.EA),
                             (coefficient(stiffness, jet("u", 2), jet("u", 2)), s.H)]:
        assert float(actual) == expected
    section = timo.rectangular_section(E=1., nu=.3, rho=1., width=.2, thickness=.05, K=5/6)
    assert float(coefficient(mass, jet("w", dt=1), jet("w", dt=1))) == pytest.approx(section.rhoA, rel=1e-14)
    assert float(coefficient(mass, jet("psi", dt=1), jet("psi", dt=1))) == pytest.approx(section.rhoI, rel=1e-14)
    assert float(Iy) == pytest.approx(section.I_y, rel=1e-14)
    assert timo.LEGACY_SHEAR_CONVENTION == "Q=KGA*(dw_dx-psi)"


def physical_operator(mass, stiffness):
    """Negative Euler derivative of T-V; exact constant-coefficient differential entries."""
    result = {}
    for matrix, sign in ((mass, -1), (stiffness, 1)):
        for ((left, dx, dt), (right, ex, et)), value in matrix.items():
            key = left, right, dx+ex, dt+et
            result[key] = result.get(key, F(0)) + sign*(-1)**(dx+dt)*value
    return {key: value for key, value in result.items() if value}


def test_variation_of_raw_general_energy_has_bishop_block_but_wrong_bending_self():
    mass, stiffness = matrices()
    operator = physical_operator(mass, stiffness)
    p = a.PARAMETERS
    area = p["b"]*p["h"]
    Iy = p["b"]*p["h"]**3/12
    Ip = Iy+p["h"]*p["b"]**3/12
    G = 1/(2*(1+p["nu"]))
    assert operator == {
        ("u", "u", 0, 2): area, ("u", "u", 2, 2): -p["nu"]**2*Ip,
        ("u", "u", 2, 0): -area, ("u", "u", 4, 0): p["nu"]**2*G*Ip,
        ("w", "w", 0, 2): area, ("w", "w", 2, 0): -G*area,
        ("w", "psi", 1, 0): G*area,
        ("psi", "psi", 0, 2): Iy, ("psi", "psi", 2, 0): -F(35, 26)*Iy,
        ("psi", "psi", 0, 0): G*area, ("psi", "w", 1, 0): -G*area,
    }


@pytest.mark.parametrize("model", ["elementary", "rayleigh_love", "bishop"])
def test_existing_axial_order_reductions_without_a_combined_spectrum(model):
    p = a.PARAMETERS
    Ip = p["b"]*p["h"]*(p["b"]**2+p["h"]**2)/12
    H = float(p["nu"]**2*Ip/(2*(1+p["nu"]))) if model == "bishop" else 0.
    J = float(p["nu"]**2*Ip) if model != "elementary" else 0.
    s = bishop.Segment(1., .01, .01, H, J)
    omega, x = .7, [.1, .6]
    # Arbitrary-frequency differential check, not a new eigenvalue calculation.
    U, U2, U4 = [bishop.basis(s, omega, x, d) for d in (0, 2, 4)]
    terms = [H*U4, (J*omega**2-s.EA)*U2, -s.m*omega**2*U]
    assert np.max(abs(sum(terms)))/max(np.max(abs(t)) for t in terms) < 1e-12
    assert U.shape[1] == (4 if model == "bishop" else 2)
    if model != "bishop":
        assert s.H == 0.0  # exact reduced order, never epsilon regularization


def test_report_does_not_promote_conditional_zeros_to_accepted_combined_model():
    result = a.audit()
    assert result["centered_mixed_coefficients_zero"]
    assert result["scientific_status"] == "COMBINED_KINEMATICS_NOT_UNIQUELY_DEFINED"
    assert result["spectrum_and_factorization"] == "NOT_RUN_KINEMATICS_HARD_GATE"
    assert result["raw_to_timoshenko_bending_rigidity_ratio"] == "35/26"
    assert result["raw_to_timoshenko_shear_rigidity_ratio"] == "6/5"
    json.dumps(result, allow_nan=False)
