"""Bounded literature checks; no chat-derived golden frequencies or ranking assertions."""
from dataclasses import replace
import json

import numpy as np
import pytest
from scipy.optimize import brentq

from scripts.analysis import reproduce_bishop_literature as cli
from scripts.lib import bishop_longitudinal as b


@pytest.fixture(scope="module")
def fixture():
    return json.loads(cli.FIXTURE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def marais(fixture):
    p = fixture["sources"]["marais"]["parameters_si"]
    segments = [b.circular_segment(L, p["E"], p["rho"], p["nu"], r)
                for L, r in zip(p["lengths"], p["radii"])]
    roots, search = b.bounded_roots(lambda f: b.characteristic(segments, f, ("C", "F")),
                                  [1, 31000], 620, 5, fixture["numerical_contract"])
    assert search["status"] == "PASS"
    return segments, roots, search


def test_source_hashes_and_printed_precision(fixture):
    _, checked = cli.source_check()
    assert set(checked) == {"marais", "popov"}
    p = fixture["sources"]["popov"]
    assert len(p["ratio_strings"]) == 30
    assert all(len(v.split(".")[1]) == 3 for v in p["ratio_strings"])
    assert p["ratio_strings"][10] == "10.998"  # checked page image, never smoothed
    assert p["parameters_printed"]["f30_exp_hz_separately_printed"] == "38557.99"


@pytest.mark.parametrize("model", ["wave", "rayleigh", "bishop"])
def test_uniform_formula_vs_boundary_and_split(model):
    s = b.speed_segment(2.006, 5206.822, .34, .0248, model)
    ends = ("UP", "UP") if s.H else ("U", "U")
    split = [replace(s, L=s.L*.4), replace(s, L=s.L*.6)]
    for exact in b.explicit_frequencies(s, 5):
        root = brentq(lambda f: b.characteristic([s], f, ends), exact*.99, exact*1.01)
        other = brentq(lambda f: b.characteristic(split, f, ends), exact*.99, exact*1.01)
        assert root == pytest.approx(exact, rel=1e-10)
        assert other == pytest.approx(root, rel=1e-10)


@pytest.mark.parametrize("model", ["wave", "rayleigh", "bishop"])
def test_rigid_mode_and_second_order_boundaries(model):
    s = b.speed_segment(2.006, 5206.822, .34, .0248, model)
    coeff = np.array([1., 0., 0., 0.] if s.H else [1., 0.])
    matrix, _ = b.boundary_matrix([s], 0., ("F", "F"), balanced=False)
    np.testing.assert_array_equal(matrix@coeff, 0.)
    np.testing.assert_array_equal(b.basis(s, 0, [0, s.L/2, s.L], 1)@coeff, 0.)
    if not s.H:
        with pytest.raises(ValueError, match="fourth order"):
            b.boundary_matrix([s], 1000., ("C", "F"))
        if s.J:
            with pytest.raises(ValueError, match="cutoff"):
                b.basis(s, np.sqrt(s.EA/s.J)*1.01, [0.])


@pytest.mark.parametrize("model", ["wave", "rayleigh"])
def test_reduced_free_positive_frequencies_use_their_own_traction(model):
    s = b.speed_segment(2.006, 5206.822, .34, .0248, model)
    for exact in b.explicit_frequencies(s, 30)[[0, 29]]:
        root = brentq(lambda f: b.characteristic([s], f, ("F", "F")), exact*.999, exact*1.001)
        assert root == pytest.approx(exact, rel=1e-10)
        coefficients, _ = b.mode_coefficients([s], root, ("F", "F"))
        amplitude = max(abs(b.basis(s, 2*np.pi*root, [0., s.L])@coefficients[0]))
        for x in [0., s.L]:
            state = b.state_basis(s, 2*np.pi*root, x)@coefficients[0]
            assert abs(state[2])/(s.EA/s.L*amplitude) < 1e-8


def test_gamma_sign_through_boundary_work():
    s = b.circular_segment(.2, 70e9, 2700, .33, .11)
    omega = 2*np.pi*4000
    coefficients = np.array([.7, -.2, .3, -.1])
    z, w = b.quadrature(100)
    x, w = (z+1)*s.L/2, w*s.L/2
    u = [b.basis(s, omega, x, d)@coefficients for d in range(3)]
    # Arbitrary polynomial virtual displacement: checks boundary signs independently.
    v, dv, ddv = 1+x+x*x, 1+2*x, 2*np.ones_like(x)
    volume = np.dot(w, s.EA*u[1]*dv+s.H*u[2]*ddv-omega**2*(s.m*u[0]*v+s.J*u[1]*dv))
    edge = []
    for xx in [0., s.L]:
        state = b.state_basis(s, omega, xx)@coefficients
        edge.append(state[2]*(1+xx+xx*xx)+state[3]*(1+2*xx))
        du = b.basis(s, omega, [xx], 1)[0]@coefficients
        d3 = b.basis(s, omega, [xx], 3)[0]@coefficients
        gamma = (-s.EA+s.J*omega**2)*du+s.H*d3
        assert gamma == pytest.approx(-state[2], rel=1e-12)
    assert volume == pytest.approx(edge[1]-edge[0], rel=1e-11)


def test_reconstructed_wave_numbers_satisfy_source_equation():
    s = b.circular_segment(.2, 70e9, 2700, .33, .11)
    for frequency in [1., 4000., 29000., 70000.]:
        omega = 2*np.pi*frequency
        a, beta = b.wave_numbers(s, omega)
        for root in (a, -a, 1j*beta, -1j*beta):
            terms = [s.H*root**4, (s.J*omega**2-s.EA)*root**2, -s.m*omega**2]
            assert abs(sum(terms))/sum(abs(t) for t in terms) < 1e-12


def test_dimensions_frequency_and_normalized_material_mapping():
    s = b.circular_segment(.2, 70e9, 2700, .33, .11)
    scaled = b.circular_segment(.2*7, 70e9, 2700, .33, .11*7)
    np.testing.assert_allclose(b.explicit_frequencies(scaled, 5)*7, b.explicit_frequencies(s, 5), rtol=1e-13)
    normalized = b.speed_segment(.2, np.sqrt(70e9/2700), .33, .22)
    np.testing.assert_allclose([s.EA/s.m, s.H/s.m, s.J/s.m], [normalized.EA, normalized.H, normalized.J], rtol=1e-14)
    wave = b.Segment(.2, s.EA, s.m)
    assert b.explicit_frequencies(wave, 1)[0] == pytest.approx(np.sqrt(s.EA/s.m)/(2*s.L))


def test_long_rod_basis_does_not_overflow():
    s = b.speed_segment(2.006, 5206.822, .34, .0248)
    a, _ = b.wave_numbers(s, 2*np.pi*1300)
    assert a*s.L > 700  # unscaled cosh overflows here
    for d in range(5):
        assert np.isfinite(b.basis(s, 2*np.pi*1300, [0., s.L/2, s.L], d)).all()


def test_printed_popov_five_and_fifteen_direct_substitution():
    for nu in [.31, .34]:
        n = np.arange(1, 31)
        eta = np.pi**2*n**2*nu**2*.0248**2/(8*2.006**2)
        rayleigh = n*5206.822/(2*2.006)/np.sqrt(1+eta)
        bishop_up = rayleigh*np.sqrt(1+eta/(2*(1+nu)))
        for model, formula in [("rayleigh", rayleigh), ("bishop", bishop_up)]:
            s = b.speed_segment(2.006, 5206.822, nu, .0248, model)
            np.testing.assert_allclose(b.explicit_frequencies(s, 30), formula, rtol=1e-14)


def test_marais_energy_orthogonality_conditions_and_norms(marais, fixture):
    segments, roots, _ = marais
    v, profiles = b.verify_modes(segments, roots, ("C", "F"), quadrature_order=128)
    assert cli.numerical_status(v, fixture["numerical_contract"]) == "PASS"
    for n in range(1, 6):
        maximum = max(abs(row["Y"]) for row in profiles if row["mode"] == n)
        assert .9999 < maximum <= 1+1e-12
    # Independent quadrature order must retain the energy and Gram result.
    other, _ = b.verify_modes(segments, roots, ("C", "F"), quadrature_order=64)
    np.testing.assert_allclose(other["mass_gram"], v["mass_gram"], atol=1e-10)


def test_marais_bounded_completeness(marais):
    segments, roots, search = marais
    count = b.marais_argument_count(segments, horizontal_samples=2048)
    assert count["count"] == len(roots)
    assert count["branch_points_outside_contour"]
    assert count["max_phase_step_rad"] < np.pi/2
    assert len(search["attempts"]) == 1
    assert all(item["evaluations"] <= 102 for item in search["attempts"][0]["brackets"])


def test_independent_state_exponential_one_addressed_root(marais, fixture):
    _, roots, _ = marais
    run = b.independent_marais_transfer(fixture["sources"]["marais"]["parameters_si"], [float(roots[0])], 50)
    assert run["roots"][0]["status"] == "PASS"
    assert run["roots"][0]["relative_to_double"] < 1e-9


@pytest.mark.parametrize("end", ["F", "C"])
def test_popov_printed_equation13_vs_boundary_matrix(end):
    s = b.speed_segment(2.006, 5206.822, .34, .0248)
    for estimate in b.explicit_frequencies(s, 30)[[0, 14, 29]]:
        bounds = (estimate*.999, estimate*1.004)
        first = brentq(lambda f: b.characteristic([s], f, (end, end)), *bounds)
        second = brentq(lambda f: b.popov_characteristic(s, f, end), *bounds)
        assert first == pytest.approx(second, rel=1e-10)


def test_rayleigh_fit_uses_rounded_table_without_ranking_gate(fixture):
    source = fixture["sources"]["popov"]
    p = source["parameters_si"]
    experimental = np.array(source["ratio_strings"], float)*p["f1_exp"]
    fit = cli.fit_rayleigh(p, experimental)
    assert fit["status"] == "PASS"
    assert .31 <= fit["nu"] <= .36
    for nu in [.31, .36, .337]:
        s = b.speed_segment(p["L"], p["c_printed"], nu, p["d"], "rayleigh")
        assert fit["mae_hz"] <= np.mean(abs(b.explicit_frequencies(s, 30)-experimental))+1e-8


def test_print_match_uses_significant_digits_not_solver_tolerance():
    near = cli.print_comparison(12711., "12710")
    assert near["print_step_hz"] == 10
    assert near["nearest_status"] == "PRINT_MATCH"
    assert cli.print_comparison(12716., "12710")["nearest_status"] == "PRINT_MISMATCH"
    assert cli.print_comparison(12716., "12710")["truncation_status"] == "PRINT_MATCH"


def test_bounded_failure_is_reported_not_silently_accepted(fixture):
    roots, report = b.bounded_roots(lambda f: f-3, [1, 5], 10, 2, fixture["numerical_contract"])
    assert report["status"] == "UNRESOLVED"
    assert len(report["attempts"]) == 2 and len(roots) == 1
    assert report["attempts"][1]["intervals"] == 20


def test_cache_rejects_stale_identity_and_modified_data(tmp_path):
    artifact = tmp_path/"frequencies.csv"
    artifact.write_text("f\n1\n", encoding="utf-8")
    cli.write_json(tmp_path/"manifest.json", {"fingerprint": "checked", "artifacts": {artifact.name: cli.sha(artifact)}})
    cli.validate_bundle(tmp_path, "checked")
    with pytest.raises(ValueError, match="Stale"):
        cli.validate_bundle(tmp_path, "changed")
    artifact.write_text("f\n2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="artifact changed"):
        cli.validate_bundle(tmp_path, "checked")


def test_plot_only_never_calls_solver(tmp_path, monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("plot-only called solver")
    for name in ("bounded_roots", "boundary_matrix", "verify_modes"):
        monkeypatch.setattr(b, name, forbidden)
    rows = [{"mode": n, "x_m": x, "Y": x} for n in range(1, 6) for x in [0., .3]]
    cli.write_csv(tmp_path/"profiles.csv", rows)
    cli.render("marais", tmp_path)
    assert (tmp_path/"marais_fig2.png").stat().st_size > 1000
