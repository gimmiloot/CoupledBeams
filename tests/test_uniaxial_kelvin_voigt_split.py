"""Linear split and cached prefix checks; no spectral searches in this suite."""
import copy
from fractions import Fraction
import numpy as np
import pytest

from scripts.analysis import check_uniaxial_kelvin_voigt_split as run
from scripts.lib import uniaxial_kelvin_voigt_split as s


@pytest.fixture(scope="module")
def data():
    config = run.read(run.CONFIG)
    spatial, kv, model, protection, _ = run.sources(config)
    fp, identity = run.identity(config, protection)
    directory = run.OUTPUT/fp
    manifest = run.read(directory/"manifest.json")
    assert manifest["identity"] == identity
    assert all(run.sha(directory/n) == h for n, h in manifest["artifact_hashes"].items())
    return config, spatial, kv, model, run.read(directory/"result.json")


def test_exact_quadratic_action_and_all_bulk_cross_blocks(data):
    c, operator, audit = s.quadratic_action(data[1])
    assert all(value == 0 for value in audit["mixed_blocks"].values())
    for p in (1j, -.001+2j):
        a = operator(p)
        assert np.count_nonzero(a[np.ix_(s.PY, s.OY)]) == 0
        assert np.count_nonzero(a[np.ix_(s.OY, s.PY)]) == 0
    assert c.CT == data[1]["coefficients"]["CT"]


@pytest.mark.parametrize("beta", [0., 45., 90.])
def test_rank_trace_arm_exchange_and_virtual_work(data, beta):
    k, c = data[2]["k_theta"], data[2]["c_theta"]
    for stiffness in (k, k+(.001+2j)*c):
        j = s.joint_operator(beta, stiffness)
        assert j.shape == (14, 28)
        assert np.linalg.matrix_rank(j/np.linalg.norm(j, axis=1)[:, None]) == 14
    assert np.linalg.matrix_rank(s.trace_map(beta)) == 9
    saved = next(x for x in data[4]["structure"]["geometries"] if x["beta"] == beta)
    assert saved["joint_virtual_work_absolute"] < 2e-12
    assert saved["power_balance_absolute"] < 2e-12
    assert saved["arm_swap_max_error"] == 0.


@pytest.mark.parametrize("beta", [0., 45., 90.])
def test_physical_rotation_moment_duality(beta):
    z = np.array([.3, -.7, .2])
    covector = np.array([-.4, .8, .9])
    for transform in s.spatial.frame_transforms(beta):
        t, n, k = transform[:3, :3].T
        rotation = z[0]*t-z[1]*n+z[2]*k
        moment = covector[0]*t-covector[1]*n+covector[2]*k
        np.testing.assert_allclose(rotation, transform[3:6, 3:6]@z, atol=1e-15)
        np.testing.assert_allclose(moment, transform[3:6, 3:6]@covector, atol=1e-15)
        assert abs(moment@rotation-covector@z) < 1e-15
        assert np.allclose(np.cross(t, n), k, atol=1e-15)


@pytest.mark.parametrize("beta", [0., 45., 90.])
def test_exact_out_joint_independence_and_same_rigid_operator(beta):
    rigid = s.out_joint(beta, rigid=True)
    for impedance in (0., 1., 2.+3j):  # algebraic coefficients, no parameter spectrum
        np.testing.assert_array_equal(s.out_joint(beta, impedance), rigid)
    np.testing.assert_array_equal(rigid, s.book_joint.joint_matrix_book(np.deg2rad(beta)))
    for value in (s.joint_operator(beta, 0.), s.joint_operator(beta, 2.+3j)):
        assert np.count_nonzero(value[np.ix_(s.PR, s.OC)]) == 0
        assert np.count_nonzero(value[np.ix_(s.OR, s.PC)]) == 0


@pytest.mark.parametrize("beta", [0., 90.])
def test_exact_integer_endpoint_ranks(beta):
    transforms = [np.rint(t) for t in s.spatial.frame_transforms(beta)]
    assert s._rank_fraction(s.joint_operator(beta, 2., transforms).astype(int)) == 14


def test_negative_rotation_coupling_is_detected(data):
    altered = copy.deepcopy(data[1])
    potential = s.rod.Polynomial.deserialize(altered["polynomials"]["V4"])
    potential += Fraction(1, 4)*s.rod.Polynomial.symbol("theta")*s.rod.Polynomial.symbol("Phi")
    altered["polynomials"]["V4"] = potential.serialize()
    with pytest.raises(ArithmeticError, match="quadratic action"):
        s.quadratic_action(altered)
    j = s.joint_operator(45., 1.)
    j[10, 3] += .25
    assert np.max(abs(j[np.ix_(s.PR, s.OC)])) == .25
    out = s.out_joint(45.)
    out[4, 2] += .25
    assert np.max(abs(out-s.out_joint(45., rigid=True))) == .25


def test_section_clamp_is_not_book_slope_clamp_and_CT_is_not_GIp(data):
    c = data[1]["coefficients"]
    # State (v,psi,Phi,Qv,Mpsi,MT): admissible section clamp with shear reaction.
    y = np.array([0., 0., 0., 1., 0., 0.])
    assert np.max(abs(y[:3])) == 0.
    assert y[1]+y[3]/c["S"] != 0.
    g_ip = (c["jp"]+c["jb"])/(2*(1+.3))
    assert abs(c["CT"]/g_ip-1) > .5


@pytest.mark.parametrize("k", range(1, 8))
def test_cached_physical_out_mode_for_both_joint_assemblies(data, k):
    root = data[4]["spectral"]["modes"][k-1]
    d = root["diagnostics"]
    assert root["omega"] > 0 and root["guard"] == (k == 7)
    assert d["rigid_kv_boundary_max"] == 0.
    assert d["alpha_exact"] == d["delta_theta_exact"] == 0.
    assert d["complex_energy_relative"] < 5e-8
    assert d["clamp_scaled"] < 1e-9 and max(d["joint_scaled"]) < 1e-9
    assert d["ode_scaled"] < 1e-9 and abs(d["mass"]-1) < 5e-7
    # Residual cancellation also needs the same arithmetic, not equality of
    # separately rounded real/complex dot products near a numerical zero.
    endpoint = np.array(d["end_states"], dtype=complex)[:, 1].ravel()
    impedance = data[2]["k_theta"]+1j*root["omega"]*data[2]["c_theta"]
    rigid = s.out_joint(45., rigid=True).astype(complex)@endpoint
    kv = s.out_joint(45., impedance).astype(complex)@endpoint
    np.testing.assert_array_equal(kv, rigid)


def test_complete_prefix_count_and_single_actual_search(data):
    result = data[4]
    spectrum = result["spectral"]
    assert result["status"] == "COMPLETE" and result["source_preservation"] == "PASS"
    assert spectrum["catalog"]["status"] == "PASS"
    assert len(spectrum["catalog"]["tim_roots"]) == spectrum["catalog"]["tim_count_bound"]["upper_count"] == 2
    assert spectrum["search"]["lower_count"] == 0
    assert spectrum["search"]["upper_count"] >= 7
    assert len(spectrum["search"]["roots"]) == 7
    assert spectrum["mass_gram_error"] < 5e-7
    original = result["prefix_reuse"]["original_solve_ledger"]
    assert original["global_prefix_searches"] == 1 and original["global_out_roots"] == 7
    assert result["solve_ledger"]["global_prefix_searches"] == 0
    assert result["solve_ledger"]["reused_out_roots"] == 7


def test_inplane_regression_is_matrix_only(data):
    for geo in data[4]["structure"]["geometries"]:
        for check in geo["spectral_matrices"]:
            assert check["parallel_old_matrix_max"] == 0.
            assert check["boundary_mixed_max"] == 0.
            assert check["out_rigid_matrix_max"] == 0.
            assert check["out_parameter_derivative_max"] == 0.
    assert data[4]["solve_ledger"]["inplane_root_solves"] == 0
    assert data[4]["solve_ledger"]["complex_correctors"] == 0


def test_cache_cannot_repeat_the_spectrum(data, monkeypatch, capsys):
    def forbidden(*args, **kwargs):
        pytest.fail("Validated cache attempted a spectral search")
    monkeypatch.setattr(s.OutOfPlane, "pole_catalog", forbidden)
    monkeypatch.setattr(s.OutOfPlane, "prefix", forbidden)
    monkeypatch.setattr(s.OutOfPlane, "mode", forbidden)
    run.main(["--compute"])
    assert '"new_spectral_solves_this_invocation": 0' in capsys.readouterr().out
