"""Fixed-case physics and saved-result checks; never repeat an eigenvalue solve."""
import numpy as np
import pytest

from scripts.analysis import compare_longitudinal_kelvin_voigt as run
from scripts.lib import contraction_joint_sensitivity as old
from scripts.lib import longitudinal_kelvin_voigt_comparison as study


@pytest.fixture(scope="module")
def inputs():
    config = run.read(run.CONFIG)
    model, previous, protection = run.check_source(config)
    return config, model, previous, protection


@pytest.fixture(scope="module")
def completed(inputs):
    config, _, _, protection = inputs
    fingerprint, identity = run.make_identity(config, protection)
    directory = run.OUTPUT/fingerprint
    manifest = run.read(directory/"manifest.json")
    assert manifest["identity"] == identity
    assert all(run.sha(directory/name) == digest for name, digest in manifest["artifact_hashes"].items())
    return run.read(directory/"result.json")


TARGETS = [(variant, k) for variant in run.VARIANTS[:2] for k in (1, 6)]


@pytest.mark.parametrize("variant,k", TARGETS)
def test_complex_operator_harmonic_recovery_and_natural_force(inputs, variant, k):
    _, model, previous, _ = inputs
    item = run.records(previous)[variant]
    case, omega = old.Case(**item["case"]), item["roots"][k-1]["omega"]
    np.testing.assert_allclose(study.complex_operator(model, case, 1j*omega),
                               old.state_operator(model, omega, case), rtol=1e-14, atol=0.)
    seg = study.reduced.segment(model, .5, variant)
    p = -.001+1j*omega
    operator = study.complex_operator(model, case, p)
    assert np.isclose(operator[0, 3]*(seg.EA+seg.J*p*p), 1., rtol=1e-14)
    assert operator[3, 0] == seg.m*p*p
    assert operator[5, 4] == -1  # Jang/Tim M'=r*p^2*theta-Q, not the EB/RLB sign.


@pytest.mark.parametrize("variant,k", TARGETS)
def test_analytic_boundary_derivative_in_both_complex_directions(inputs, variant, k):
    config, model, previous, _ = inputs
    item = run.records(previous)[variant]
    provider = study.ReducedComplexBasis(model, [.5, .5], item["roots"][k-1]["omega"],
        old.Case(**item["case"]), previous["k_theta"], previous["t_ref"])
    assert provider.matrix(1j*provider.omega*provider.tref, 0.).shape == (12, 12)
    assert study.provider_preflight(provider, previous["c_theta"], config["policy"])["status"] == "PASS"


@pytest.mark.parametrize("variant", run.VARIANTS)
def test_mass_is_theory_specific_without_new_shape_solve(inputs, variant):
    config, model, previous, _ = inputs
    item = run.records(previous)[variant]
    case, root = old.Case(**item["case"]), item["roots"][5]
    d = study.saved_elastic_diagnostics(model, [.5, .5], case, root, previous["k_theta"], 200)
    y, dx, weight, c = d["values"], d["gradients"], d["weights"], model.coefficients
    u, w, theta = (0, 2, 3) if case.is_mh else (0, 1, 2)
    common = np.sum(weight*(c["m"]*(abs(y[:, :, u])**2+abs(y[:, :, w])**2)+c["r"]*abs(y[:, :, theta])**2))
    extra = (np.sum(weight*c["j"]*abs(y[:, :, 1])**2) if case.is_mh else
             np.sum(weight*model.section.nu**2*c["j"]*abs(dx[:, :, 0])**2) if variant == "rayleigh_love_planar" else 0.)
    assert abs(common+extra-1) < config["policy"]["mass_normalization_tol"]
    assert abs(d["mass"]-common-extra) < 1e-14
    assert np.isclose(d["contraction_or_gradient_mass"], extra, rtol=1e-13, atol=1e-16)
    if variant != "elementary":
        assert extra > 3e-5  # Omitting this inertia cannot pass the independent mass calculation.


@pytest.mark.parametrize("variant,k", TARGETS)
def test_saved_complex_solution_full_acceptance(completed, variant, k):
    value = next(v for v in completed["complex_modes"] if v["variant"] == variant and v["k"] == k)
    assert value["status"] == "PASS" and all(value["gates"].values())
    p = complex(value["p"]["real"], value["p"]["imag"])
    d = value["diagnostics"]
    delta = complex(d["delta_theta"]["real"], d["delta_theta"]["imag"])
    damping = completed["c_theta"]*abs(delta)**2
    residual = p*p*d["mass"]+p*damping+d["elastic_energy"]
    assert abs(residual)/(abs(p*p)*d["mass"]+abs(p)*damping+d["elastic_energy"]) < 5e-8
    assert abs(-p.real/(damping/(2*d["mass"]))-1) < 5e-8
    assert value["elastic_frequency_relative_error"] < 2e-8
    assert value["mass_MAC_with_elastic"] > .98


def test_activity_projection_and_exact_solve_budget(completed):
    assert completed["status"] == "COMPLETE"
    assert completed["solve_ledger"] == dict(elastic_spectrum_searches=0, elastic_shape_solves=0,
        reused_elastic_modes=18, reused_MH_complex_roots=2, elastic_provider_correctors=4, new_damped_correctors=4)
    for variant in run.VARIANTS:
        rows = [r for r in completed["stage_a"]["table"] if r["variant"] == variant]
        assert [r["k"] for r in rows if r["active"]] == [1, 3, 6]
        assert all(r["alpha_pred_vs_MH_percent"] is None for r in rows if not r["active"])
    for value in completed["stage_a"]["correspondence"].values():
        assert value["status"] == "PASS"
        assert list(np.argmax(value["matrix"], axis=1)) == list(range(6))


def test_cache_blocks_every_spectral_solver(monkeypatch, completed, capsys):
    def forbidden(*args, **kwargs):
        pytest.fail("A validated bundle must not launch eigenvalue or shape solves")
    monkeypatch.setattr(study, "complex_target", forbidden)
    monkeypatch.setattr(old, "prefix_roots", forbidden)
    monkeypatch.setattr(old, "elastic_mode", forbidden)
    monkeypatch.setattr(old, "complex_mode", forbidden)
    run.main(["--compute"])
    assert '"new_root_correctors_this_invocation": 0' in capsys.readouterr().out


def test_source_artifact_mismatch_never_falls_back_to_solving(monkeypatch, inputs):
    config = inputs[0]
    original = run.sha
    def wrong(path):
        return "invalid" if path.name == "profiles.csv" else original(path)
    monkeypatch.setattr(run, "sha", wrong)
    with pytest.raises(ValueError, match="Source artifact mismatch"):
        run.check_source(config)


def test_complex_stage_requires_stage_a(inputs, tmp_path):
    config, model, previous, _ = inputs
    with pytest.raises(ArithmeticError, match="Stage A hard gate"):
        run.stage_b(model, previous, config, {"stage_a": {"status": "NUMERICAL_UNRESOLVED"}}, {}, tmp_path)


def test_no_new_mh_provider(inputs):
    with pytest.raises(ValueError, match="never solved here"):
        study.complex_operator(inputs[1], old.Case(45., "FREE"), 1j)
