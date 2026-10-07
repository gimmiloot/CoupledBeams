"""Lightweight planar action/RHS controls; full 5*T1 runs belong to the CLI."""
from pathlib import Path
import copy
import hashlib
import json
import math
import sys

import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib import weakly_nonlinear_planar_dynamics as planar
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section


ROOT = Path(__file__).resolve().parents[1]
PLANAR_INDICES = (0, 1, 5, 6)


@pytest.fixture(scope="module")
def audited_model():
    return rod.derive_polynomials()


@pytest.fixture(scope="module")
def coefficients():
    # C_T is not present in the planar restriction; retain the accepted G20
    # value rather than introduce a new torsional prescription.
    from scripts.lib import yartsev_ch2_monoclinic_rod as book
    G = 1/2.6
    material = book.BookMaterial(E1_real=1., E2_real=1., G12_real=G,
        G13_real=G, G23_real=G, nu12=.3, rho=1., eta1=0., eta2=0.,
        eta12=0., eta13=0., eta23=0.)
    point = book.make_rod_point(0., material=material,
        geometry=book.Geometry(a=.20, b=.05, length=1., shear_factor=5/6))
    return rod.RodCoefficients.rectangular(1., 1., .3, .20, .05,
                                          float(point.torsion.C_T.real))


@pytest.fixture(scope="module")
def old_reference():
    """Read the preserved analytic reference; no old research recomputation."""
    bundle = ROOT/"results/mindlin_herrmann_timoshenko_single_rod/342ce44bff81c36f"
    if not bundle.is_dir():
        pytest.skip("Preserved local single-rod reference bundle unavailable")
    manifest = json.loads((bundle/"manifest.json").read_text(encoding="utf-8"))
    for name, digest in manifest["artifact_hashes"].items():
        assert hashlib.sha256((bundle/name).read_bytes()).hexdigest() == digest
    data = json.loads((bundle/"result.json").read_text(encoding="utf-8"))
    model = mh.project_jang_reduced_rectangular(rectangular_section(E=1., nu=.3,
        rho=1., width=.20, thickness=.05, K=5/6))
    omega = data["timoshenko"]["roots"][0]["omega"]
    coef = mh.finite_mode(model, 1., omega, "timoshenko", 200)["coefficients"]
    stationary = brentq(lambda x: (mh.finite_state_basis(model, 1., omega, x,
        "timoshenko", 1)@coef)[0], .4, .6, xtol=1e-14)
    peak = (mh.finite_state_basis(model, 1., omega, stationary, "timoshenko")@coef)[0]
    coef = coef/peak  # One common multiplier for the complete eigenpair.
    def pair(points):
        return mh.finite_state_basis(model, 1., omega, points, "timoshenko")@coef
    def initial_fields(points, amplitude=1.):
        points = np.atleast_1d(points)
        values = pair(points)
        return np.column_stack((points*0, amplitude*values[:, 0],
                                amplitude*values[:, 1], points*0))
    return {"data": data, "omega": omega, "T1": 2*math.pi/omega,
            "pair": pair, "initial_fields": initial_fields, "max_location": stationary}


@pytest.fixture(scope="module")
def discretization(coefficients, audited_model):
    return planar.PlanarGalerkin(coefficients, p=6, length=1., model=audited_model)


def physical_state(discretization):
    """All four independent fields/velocities, deterministic small values."""
    def qfield(x):
        bubble = x*(1-x)
        return np.column_stack((2e-4*bubble*(2*x-1), 3e-3*bubble,
                                4e-3*bubble*(2*x-1), 3e-4*bubble))
    def vfield(x):
        bubble = x*(1-x)
        return np.column_stack((3e-4*bubble*(1+.3*(2*x-1)),
                                2e-3*bubble*(1+.4*(2*x-1)),
                                3e-3*bubble*(.5+2*x-1),
                                -2e-4*bubble*(1+.2*(2*x-1))))
    return discretization.project(qfield), discretization.project(vfield)


def physical_derivative(discretization, q, order=0):
    basis = discretization.basis_at(discretization.x, derivative=order)
    n = discretization.n
    return np.column_stack([basis[field]@q[i*n:(i+1)*n]
        for i, field in enumerate(discretization.fields)])


def embedded_jet(discretization, q, velocity, acceleration):
    planar_values = (physical_derivative(discretization, q),
        physical_derivative(discretization, q, 1),
        physical_derivative(discretization, velocity),
        physical_derivative(discretization, q, 2),
        physical_derivative(discretization, velocity, 1),
        physical_derivative(discretization, acceleration))
    full = []
    for values in planar_values:
        embedded = np.zeros((discretization.nq, 7))
        embedded[:, PLANAR_INDICES] = values
        full.append(embedded)
    return [rod.FieldJet(*(part[j] for part in full)) for j in range(discretization.nq)]


@pytest.mark.parametrize("p", (6, 16, 24, 32))
def test_basis_four_independent_fields_essential_BC_without_slope_constraints(coefficients, audited_model, p):
    d = planar.PlanarGalerkin(coefficients, p=p, model=audited_model)
    assert d.fields == ("u", "w", "theta", "c")
    assert d.n == p-1 and d.ndof == 4*(p-1)
    endpoints = d.basis_at(np.array([0., 1.]))
    slopes = d.basis_at(np.array([0., 1.]), derivative=1)
    for field in d.fields:
        np.testing.assert_allclose(endpoints[field], 0., rtol=0, atol=2e-11)
        assert np.linalg.norm(slopes[field]) > 0
    # No c=-nu*u_s or theta=w_s algebraic constraint reduces the dimension.
    for i in range(4):
        q = np.zeros(d.ndof); q[i*d.n] = 1e-7
        values = d.reconstruct(q, np.array([.27, .61]))
        assert np.linalg.norm(values[:, i]) > 0
        assert np.count_nonzero(values[:, [j for j in range(4) if j != i]]) == 0


def test_action_quadrature_exactness_and_independent_quadrature_increase(discretization, coefficients, audited_model):
    d = discretization
    assert d.nq >= 2*d.p+1  # Quartic action degree <=4*p; Gaussian degree2*nq-1.
    q, velocity = physical_state(d)
    jets = embedded_jet(d, q, velocity, np.zeros(d.ndof))
    density = [rod.polynomial_evaluate(jet, coefficients, audited_model) for jet in jets]
    reference = d.weights@np.array([value["T"]+value["V"] for value in density])
    assert d.energy(q, velocity) == pytest.approx(reference, rel=2e-12, abs=2e-20)
    higher = planar.PlanarGalerkin(coefficients, p=d.p, nq=d.nq+8, model=audited_model)
    higher_q = higher.project(lambda x: d.reconstruct(q, x))
    higher_v = higher.project(lambda x: d.reconstruct(velocity, x))
    assert higher.energy(higher_q, higher_v) == pytest.approx(reference, rel=2e-12, abs=2e-20)


def test_action_assembly_matches_independent_continuum_weak_projection(discretization, coefficients, audited_model):
    d = discretization
    q, velocity = physical_state(d)
    acceleration = d.project(lambda x: np.column_stack((x*(1-x)*.001,
        x*(1-x)*-.002, x*(1-x)*.003, x*(1-x)*.004)))
    jets = embedded_jet(d, q, velocity, acceleration)
    residual = np.array([rod.polynomial_evaluate(jet, coefficients, audited_model)["residual"][list(PLANAR_INDICES)] for jet in jets])
    basis = d.basis_at(d.x)
    reference = np.concatenate([basis[field].T@(d.weights*residual[:, i]) for i, field in enumerate(d.fields)])
    actual = d.mass_matrix(q)@acceleration+d.inertial_terms(q, velocity)+d.potential(q)["gradient"]
    np.testing.assert_allclose(actual, reference, rtol=2e-11, atol=2e-16)
    np.testing.assert_allclose(d.weak_residual(q, velocity, acceleration), reference, rtol=2e-11, atol=2e-16)


def test_true_variable_mass_and_both_coordinate_inertia_terms(discretization, coefficients):
    d = discretization
    q, velocity = physical_state(d)
    mass0, mass = d.mass_matrix(np.zeros(d.ndof)), d.mass_matrix(q)
    np.testing.assert_allclose(mass, mass.T, rtol=0, atol=2e-15)
    assert np.linalg.eigvalsh(mass0).min() > 0
    assert np.linalg.norm(mass-mass0) > 0
    values, velocities = d.reconstruct(q, d.x), d.reconstruct(velocity, d.x)
    basis = d.basis_at(d.x)
    c, ct, theta_t = values[:, 3], velocities[:, 3], velocities[:, 2]
    expected = np.zeros(d.ndof)
    expected[2*d.n:3*d.n] = basis["theta"].T@(d.weights*2*coefficients.jp*(1+c)*ct*theta_t)
    expected[3*d.n:] = basis["c"].T@(d.weights*-coefficients.jp*(1+c)*theta_t**2)
    np.testing.assert_allclose(d.inertial_terms(q, velocity), expected, rtol=2e-12, atol=2e-20)
    force = expected+d.potential(q)["gradient"]
    np.testing.assert_allclose(d.acceleration(q, velocity), np.linalg.solve(mass, -force), rtol=2e-12, atol=2e-14)


def test_energy_RHS_identity_using_audited_density_derivatives(discretization, coefficients, audited_model):
    d = discretization
    q, velocity = physical_state(d)
    acceleration = d.acceleration(q, velocity)
    derivative = (audited_model.T4+audited_model.V4).total_derivative("t")
    jets = embedded_jet(d, q, velocity, acceleration)
    rate = d.weights@np.array([derivative.evaluate(jet.values()|coefficients.values()) for jet in jets])
    # Use the unsigned work contributions. Spatial parity may cancel each
    # signed dot product before the energy identity is compared; that would
    # make a near-zero denominator an unsuitable numerical scale.
    scale = np.sum(np.abs(velocity*d.potential(q)["gradient"]))+np.sum(np.abs(velocity*(d.mass_matrix(q)@acceleration)))
    assert abs(rate)/max(scale, 1e-30) < 2e-12


def test_analytic_RHS_Jacobian_against_central_spot_check(discretization):
    d = discretization
    q, velocity = physical_state(d)
    state = np.r_[q, velocity]
    analytic = d.jacobian(0., state)
    columns = (0, d.n+1, 2*d.n, 3*d.n+2, d.ndof+2*d.n+1, d.ndof+3*d.n)
    for column in columns:
        h = 1e-7
        plus, minus = state.copy(), state.copy()
        plus[column] += h; minus[column] -= h
        observed = (d.rhs(0., plus)-d.rhs(0., minus))/(2*h)
        error = np.linalg.norm(observed-analytic[:, column])/max(np.linalg.norm(analytic[:, column]), 1.)
        assert error < 3e-7


def test_zero_state_and_pure_axial_RHS_are_linear(discretization):
    d = discretization
    zero = np.zeros(2*d.ndof)
    np.testing.assert_array_equal(d.rhs(0., zero), zero)
    assert d.energy(zero[:d.ndof], zero[d.ndof:]) == 0
    q, velocity = physical_state(d)
    q[d.n:3*d.n] = 0; velocity[d.n:3*d.n] = 0
    state = np.r_[q, velocity]
    np.testing.assert_allclose(d.rhs(0., state), d.jacobian(0., zero)@state, rtol=2e-12, atol=2e-15)
    np.testing.assert_array_equal(d.acceleration(q, velocity)[d.n:3*d.n], np.zeros(2*d.n))


@pytest.mark.parametrize("p,relative_tolerance", ((16, 4e-5), (24, 3e-7), (32, 3e-9)))
def test_family_specific_linear_frequency_convergence(coefficients, audited_model, old_reference, p, relative_tolerance):
    d = planar.PlanarGalerkin(coefficients, p=p, model=audited_model)
    for block, family in (("mh", "mh"), ("timoshenko", "timoshenko")):
        roots = np.array([r["omega"] for r in old_reference["data"][family]["roots"][:3]])
        obtained = d.linear_eigenpairs(block=block)["omega"][:3]
        tol = relative_tolerance if block == "mh" else 3e-9
        np.testing.assert_allclose(obtained, roots, rtol=tol, atol=0)


@pytest.mark.parametrize("p", (16, 24, 32))
def test_same_continuous_initial_pair_common_normalization_and_projection(coefficients, audited_model, old_reference, p):
    reference = old_reference
    grid = np.linspace(0., 1., 1001)
    pair = reference["pair"](grid)
    assert reference["max_location"] == pytest.approx(.5, abs=1e-12)
    assert np.max(np.abs(pair[:, 0])) == pytest.approx(1., abs=2e-14)
    # A single original eigenvector scaling fixes theta; no separate angle norm.
    d = planar.PlanarGalerkin(coefficients, p=p, model=audited_model)
    projected = d.project(reference["initial_fields"])
    assert np.all(projected[:d.n] == 0) and np.all(projected[3*d.n:] == 0)
    values = d.reconstruct(projected, grid)
    target = reference["initial_fields"](grid)
    for field in (1, 2):
        relative = np.linalg.norm(values[:, field]-target[:, field])/np.linalg.norm(target[:, field])
        assert relative < 2e-10
    assert reference["T1"] == pytest.approx(2*math.pi/reference["omega"], rel=1e-15)


def test_full_semidiscrete_linear_reference_exact_time_against_short_integrator(discretization):
    d = discretization
    q, velocity = physical_state(d)
    zero = np.zeros(2*d.ndof)
    matrix = d.jacobian(0., zero)
    times = np.linspace(0., .004, 9)
    exact = d.linear_reference(q, velocity, times)
    solved = solve_ivp(lambda t, y: matrix@y, (0., times[-1]), np.r_[q, velocity],
        method="Radau", jac=matrix, t_eval=times, rtol=1e-10, atol=1e-15, max_step=.0005)
    assert solved.success
    np.testing.assert_allclose(solved.y[:d.ndof].T, exact["q"], rtol=2e-8, atol=2e-14)
    np.testing.assert_allclose(solved.y[d.ndof:].T, exact["velocity"], rtol=2e-8, atol=2e-12)


def test_pure_axial_short_time_control_matches_full_linear_reference(discretization):
    d = discretization
    q, velocity = physical_state(d)
    q[d.n:3*d.n] = 0; velocity[d.n:3*d.n] = 0
    times = np.linspace(0., .004, 7)
    exact = d.linear_reference(q, velocity, times)
    solved = solve_ivp(d.rhs, (0., times[-1]), np.r_[q, velocity],
        method="Radau", jac=d.jacobian, t_eval=times,
        rtol=1e-10, atol=1e-15, max_step=.0005)
    assert solved.success
    np.testing.assert_allclose(solved.y[:d.ndof].T, exact["q"], rtol=2e-8, atol=2e-14)
    np.testing.assert_allclose(solved.y[d.ndof:].T, exact["velocity"], rtol=2e-8, atol=2e-12)
    np.testing.assert_array_equal(solved.y[d.n:3*d.n], np.zeros((2*d.n, len(times))))


def test_four_field_embedding_has_zero_out_of_plane_residuals(discretization, coefficients, audited_model):
    d = discretization
    q, velocity = physical_state(d)
    acceleration = d.acceleration(q, velocity)
    for jet in embedded_jet(d, q, velocity, acceleration)[::4]:
        residual = rod.polynomial_evaluate(jet, coefficients, audited_model)["residual"]
        np.testing.assert_array_equal(residual[[2, 3, 4]], np.zeros(3))
    # No transverse perturbation or stability experiment is performed here.


def test_tiny_planar_nonlinear_smoke_keeps_all_four_fields(discretization):
    d = discretization
    q = d.project(lambda x: np.column_stack((x*0, 1e-5*x*(1-x),
        2e-5*x*(1-x)*(2*x-1), x*0)))
    initial = np.r_[q, np.zeros(d.ndof)]
    solved = solve_ivp(d.rhs, (0., .004), initial, method="Radau", jac=d.jacobian,
        t_eval=np.linspace(0., .004, 5), rtol=1e-9, atol=1e-15, max_step=.0005)
    assert solved.success and np.all(np.isfinite(solved.y))
    assert solved.y.shape == (2*d.ndof, 5)
    # u/c are retained; the nonlinear equations are permitted to excite them.
    assert np.linalg.norm(solved.y[:d.n, -1]) > 0
    assert np.linalg.norm(solved.y[3*d.n:d.ndof, -1]) > 0
    energy = np.array([d.energy(y[:d.ndof], y[d.ndof:]) for y in solved.y.T])
    assert np.max(np.abs(energy/energy[0]-1)) < 1e-6


def test_potential_remains_quartic_action_instead_of_full_trigonometric_model(discretization, coefficients, audited_model):
    d = discretization
    q = d.project(lambda x: np.column_stack((x*0, .008*4*x*(1-x),
                                             .09*4*x*(1-x), x*0)))
    jets = embedded_jet(d, q, np.zeros(d.ndof), np.zeros(d.ndof))
    quartic = d.weights@np.array([rod.polynomial_evaluate(jet, coefficients, audited_model)["V"] for jet in jets])
    full = d.weights@np.array([rod.full_evaluate(jet, coefficients)["V"] for jet in jets])
    assert d.potential(q)["V"] == pytest.approx(quartic, rel=2e-12, abs=2e-20)
    assert not math.isclose(quartic, full, rel_tol=1e-8, abs_tol=1e-18)


def test_locked_pilot_config_has_bounded_grid_and_no_physical_condensation():
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as pilot
    config = json.loads(pilot.CONFIG.read_text(encoding="utf-8"))
    assert config["fields"] == ["u", "w", "theta", "c"]
    assert config["spatial"]["degrees"] == [16, 24, 32]
    assert config["spatial"]["allowed_extra_degree"] == 48
    assert config["amplitude_over_h"] == [.05, .025]
    assert config["periods"] == 5
    assert config["budget"]["status"] == "FIXED_AFTER_SMOKE"
    assert config["budget"]["total_wall_seconds"] == 2400
    assert config["budget"]["per_case_wall_seconds"] == 480
    assert config["integrator"]["method"] == "Radau"
    assert config["integrator"]["blas_threads"] == 1
    for flag in ("static_initial_correction", "nonlinear_modal_truncation", "phase_alignment",
                 "period_fit", "energy_classification", "out_of_plane_stability", "parameter_maps"):
        assert config["semantics"][flag] is False


@pytest.mark.parametrize("changed_field", ("amplitude", "spatial_degree", "quadrature", "time_step", "tolerance", "model"))
def test_cache_identity_changes_with_declared_inputs(tmp_path, changed_field):
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as pilot
    config = json.loads(pilot.CONFIG.read_text(encoding="utf-8"))
    path = tmp_path/"config.json"
    pilot.write_json(path, config)
    baseline, _ = pilot.identity(path)
    changed = copy.deepcopy(config)
    if changed_field == "amplitude":
        changed["amplitude_over_h"][0] = .049
    elif changed_field == "spatial_degree":
        changed["spatial"]["degrees"][0] = 18
    elif changed_field == "quadrature":
        changed["spatial"]["quadrature"] = "3*p+7 Gauss points; identity test only"
    elif changed_field == "time_step":
        changed["time_levels"]["tight"]["max_step_cutoff_period_fraction"] /= 2
    elif changed_field == "tolerance":
        changed["time_levels"]["tight"]["rtol"] /= 2
    else:
        changed["model"] += "-identity-test"
    pilot.write_json(path, changed)
    altered, _ = pilot.identity(path)
    assert altered != baseline


def test_cache_identity_includes_actual_model_file_hash(tmp_path, monkeypatch):
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as pilot
    config = json.loads(pilot.CONFIG.read_text(encoding="utf-8"))
    path = tmp_path/"config.json"
    pilot.write_json(path, config)
    baseline, _ = pilot.identity(path)
    original = pilot.sha
    monkeypatch.setattr(pilot, "sha", lambda p: "f"*64 if Path(p).name == "weakly_nonlinear_spatial_rod.py" else original(p))
    altered, _ = pilot.identity(path)
    assert altered != baseline


def synthetic_cached_bundle(tmp_path):
    """Clearly synthetic cache IO fixture, never an accepted trajectory."""
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as pilot
    key, identity = pilot.identity()
    bundle = tmp_path/key
    bundle.mkdir()
    config = copy.deepcopy(identity["config"])
    name = pilot.case_name(config["spatial"]["degrees"][-1], .05, "tight")
    folder = bundle/"cases"/name
    folder.mkdir(parents=True)
    times = np.linspace(0., .1, 5)
    observations = np.zeros((len(times), 3, 4))
    observations[:, 1, 1] = .0025*np.cos(.317474290788*times)
    np.savez_compressed(folder/"trajectory.npz", time=times, observations=observations,
                        energy_drift=np.zeros(len(times)))
    summary = {"statuses": {"UNIT_TEST_CACHE": "SYNTHETIC"}, "config": config,
               "cases": {name: {"amplitude": .0025}},
               "initial_eigenpair": {"omega": .317474290788, "T1": 19.7911625902},
               "synthetic_test_fixture": True}
    pilot.write_json(bundle/"summary.json", summary)
    pilot.write_json(bundle/"manifest.json", {"identity": identity,
        "artifact_hashes": {"summary.json": pilot.sha(bundle/"summary.json"),
                            f"cases/{name}/trajectory.npz": pilot.sha(folder/"trajectory.npz")}})
    return bundle, summary


@pytest.mark.parametrize("action", ("compute", "report-only", "plot-only"))
def test_cached_CLI_actions_do_zero_integrations_roots_and_derivations(tmp_path, monkeypatch, capsys, action):
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as pilot
    bundle, summary = synthetic_cached_bundle(tmp_path)
    def forbidden(*args, **kwargs):
        raise AssertionError("Cached/report/plot path evaluated dynamics, roots or symbolic derivation")
    monkeypatch.setattr(pilot, "integrate_case", forbidden)
    monkeypatch.setattr(pilot, "run_compute", forbidden)
    monkeypatch.setattr(rod, "derive_polynomials", forbidden)
    monkeypatch.setattr(mh, "finite_roots", forbidden)
    monkeypatch.setattr(mh, "finite_mode", forbidden)
    if action == "compute":
        args = ["planar-pilot", "--compute", "--output-dir", str(tmp_path)]
    else:
        args = ["planar-pilot", "--"+action, str(bundle)]
    monkeypatch.setattr(sys, "argv", args)
    pilot.main()
    reported = json.loads(capsys.readouterr().out.strip())
    assert reported["statuses"] == summary["statuses"]
    assert reported["time_integrations"] == reported["root_solves"] == reported["symbolic_derivations"] == 0
    if action == "plot-only":
        assert len(reported["figures"]) == 3
        for path in reported["figures"]:
            assert Path(path).is_file() and Path(path).with_suffix(".png").is_file()


def test_cache_rejects_corrupted_artifacts_and_changed_identity(tmp_path):
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as pilot
    bundle, _ = synthetic_cached_bundle(tmp_path)
    manifest = json.loads((bundle/"manifest.json").read_text(encoding="utf-8"))
    expected = copy.deepcopy(manifest["identity"])
    expected["version"] += "-different"
    with pytest.raises(ValueError, match="identity"):
        pilot.validate_bundle(bundle, expected)
    (bundle/"summary.json").write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="artifact hash"):
        pilot.validate_bundle(bundle)


def test_runtime_controls_preserve_existing_audit_and_reference_artifacts(discretization):
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as pilot
    config = json.loads(pilot.CONFIG.read_text(encoding="utf-8"))
    protected = [ROOT/"scripts/lib/weakly_nonlinear_spatial_rod.py",
                 ROOT/"docs/theory/weakly_nonlinear_spatial_rod_expansion_generated.md",
                 ROOT/"scripts/lib/mindlin_herrmann_longitudinal.py"]
    for key in ("audit_bundle", "linear_reference_bundle"):
        bundle = ROOT/config[key]
        protected.extend((bundle/"manifest.json", bundle/"result.json"))
    before = {path: pilot.sha(path) for path in protected}
    q, velocity = physical_state(discretization)
    discretization.rhs(0., np.r_[q, velocity])
    pilot.identity()
    assert {path: pilot.sha(path) for path in protected} == before


def read_saved_snapshot_energy(case, config):
    """Independent read-only audit of completed saved trajectories.

    The raw Legendre basis and its exact Gram matrix are rebuilt here without
    PlanarGalerkin, energy(), mass_matrix(), or a time integration. The stored
    audited seven-field T4/V4 polynomials provide the continuum energy density.
    """
    from numpy.polynomial.legendre import Legendre, leggauss, legvander
    from scipy.linalg import solve_triangular
    case = Path(case)
    meta = json.loads((case/"case.json").read_text(encoding="utf-8"))
    accepted = json.loads((ROOT/config["audit_bundle"]/"result.json").read_text(encoding="utf-8"))
    coefficients = accepted["coefficients"]
    T = rod.Polynomial.deserialize(accepted["polynomials"]["T4"])
    V = rod.Polynomial.deserialize(accepted["polynomials"]["V4"])
    trajectory = case/"trajectory.npz"
    assert hashlib.sha256(trajectory.read_bytes()).hexdigest() == meta["artifact_hashes"]["trajectory.npz"]
    p, n, length = meta["p"], meta["p"]-1, config["material_geometry"]["L"]
    nodes, weights = leggauss(2*p+1)
    weights *= length/2
    vander = legvander(nodes, p)
    B = vander[:, :n]-vander[:, 2:n+2]
    D = np.column_stack([(Legendre.basis(k)-Legendre.basis(k+2)).deriv()(nodes)
                         for k in range(n)])*2/length
    # Orthogonality of the Legendre polynomials gives exact basis integrals.
    gram = np.zeros((n, n))
    for k in range(n):
        gram[k, k] = length/(2*k+1)+length/(2*k+5)
        if k+2 < n:
            gram[k, k+2] = gram[k+2, k] = -length/(2*k+5)
    with np.load(trajectory, allow_pickle=False) as data:
        indices = data["snapshot_indices"]
        rawq = data["raw_snapshot_coefficients"]
        white_velocity = data["velocity"][indices]
        masses = (coefficients["m"], coefficients["m"], coefficients["jp"], coefficients["jp"])
        raw_velocity = np.column_stack([
            white_velocity[:, j*n:(j+1)*n]@solve_triangular(
                np.linalg.cholesky(mass*gram).T, np.eye(n), lower=False).T
            for j, mass in enumerate(masses)])
        energies = data["energy"][indices]
        snapshots, points = data["snapshots"], data["snapshot_points"]
        S = legvander(2*points/length-1, p)
        S = S[:, :n]-S[:, 2:n+2]
        independent = []
        max_endpoint = max_snapshot_difference = 0.
        for row, (q, velocity) in enumerate(zip(rawq, raw_velocity)):
            physical_q = np.column_stack([B@q[j*n:(j+1)*n] for j in range(4)])
            physical_qs = np.column_stack([D@q[j*n:(j+1)*n] for j in range(4)])
            physical_velocity = np.column_stack([B@velocity[j*n:(j+1)*n] for j in range(4)])
            physical_snapshot = np.column_stack([S@q[j*n:(j+1)*n] for j in range(4)])
            max_endpoint = max(max_endpoint, float(np.max(np.abs(physical_snapshot[[0, -1]]))))
            max_snapshot_difference = max(max_snapshot_difference,
                float(np.max(np.abs(physical_snapshot-snapshots[row]))))
            integral = 0.
            for node, weight in enumerate(weights):
                values = {name: 0. for name in rod.SYMBOL_ORDER}
                values.update(coefficients)
                for j, name in enumerate(("u", "w", "theta", "c")):
                    values[name] = float(physical_q[node, j])
                    values[name+"_s"] = float(physical_qs[node, j])
                    values[name+"_t"] = float(physical_velocity[node, j])
                integral += weight*(T.evaluate(values)+V.evaluate(values))
            independent.append(integral)
        times = data["time"]
        assert np.all(np.isfinite(times)) and np.all(np.diff(times) > 0)
        assert len(rawq) == len(indices) == len(energies)
        assert np.all(energies > 0)
        np.testing.assert_array_equal(snapshots[0, :, [0, 3]], 0.)
        return {"energy_relative_max": float(np.max(np.abs(np.array(independent)/energies-1))),
                "endpoint_essential_max": max_endpoint,
                "raw_snapshot_physical_difference_max": max_snapshot_difference,
                "snapshot_count": len(indices), "time_end": float(times[-1]),
                "generated_u_max": float(np.max(np.abs(snapshots[1:, :, 0]))),
                "generated_c_max": float(np.max(np.abs(snapshots[1:, :, 3])))}


@pytest.mark.parametrize("case_name", ("p16_Aoverh0p05_tight", "p32_Aoverh0p05_tight", "p32_Aoverh0p025_tight"))
def test_saved_snapshot_energy_independent_source_polynomial_reader(case_name):
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as pilot
    config = json.loads(pilot.CONFIG.read_text(encoding="utf-8"))
    case = ROOT/"results/weakly_nonlinear_planar_time_pilot/c97287772bc461ef/cases"/case_name
    if not (case/"case.json").is_file() or not (case/"trajectory.npz").is_file():
        pytest.skip("Requested local trajectory is unavailable or still integrating")
    meta = json.loads((case/"case.json").read_text(encoding="utf-8"))
    if meta["time_end"] < meta["target_time_end"]*(1-32*np.finfo(float).eps):
        pytest.skip("Saved trajectory has not reached the requested final time")
    result = read_saved_snapshot_energy(case, config)
    assert result["snapshot_count"] >= 2
    assert result["energy_relative_max"] < 2e-12
    assert result["endpoint_essential_max"] < 2e-12
    assert result["raw_snapshot_physical_difference_max"] < 2e-12
