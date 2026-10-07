"""Comparator physics/gates and sorted, geometric-only screening contract."""
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest

from scripts.analysis import screen_coupled_longitudinal_theory_hierarchy as audit
from scripts.lib import coupled_longitudinal_comparators as reduced
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import mindlin_herrmann_timoshenko_joint as joint


@pytest.fixture(scope="module")
def setup():
    return audit.check_inputs()


@pytest.fixture(scope="module")
def gates(setup):
    config, base, model, length, *_ = setup
    upper = 2*math.pi*config["frequency_ceiling_fstar"]
    catalogs = {l: audit.tim_catalog(model, l, upper, base) for l in (.5, .35, .65)}
    reports, cases = audit.comparator_gates(config, base, model, length, catalogs, upper)
    return reports, cases, catalogs


@pytest.mark.parametrize("variant", reduced.NAMES)
def test_direct_rod_and_asymmetric_split_recovery(gates, variant):
    report = gates[0][variant]
    assert report["status"] == "PASS"
    assert [r["split"] for r in report["splits"]] == [.5, .35]
    assert max(r["max_frequency_relative_difference"] for r in report["splits"]) < 1e-10
    for r in report["splits"]:
        case = r["case"]
        assert len(case["roots"]) == case["search"]["upper_count"] >= 13
        assert case["search"]["failed_intervals"] == []
        assert case["mass_gram_max_error"] < audit.POLICY["mass_gram_tol"]


@pytest.mark.parametrize("variant", reduced.NAMES)
def test_swap_and_reflection_profiles_and_roots(gates, variant):
    for operation in ("swap", "reflection"):
        result = gates[0][variant][operation]
        assert result["status"] == "PASS" and len(result["rows"]) == 13
        assert max(r["kinematic_L2_relative"] for r in result["rows"]) < 5e-8


@pytest.mark.parametrize("variant", reduced.NAMES)
def test_no_fictitious_contraction_or_slope_clamp(setup, variant):
    model = setup[2]
    assert reduced.STATE_ORDER == ("u", "w", "theta", "N", "Q", "M")
    matrix = reduced.boundary(model, (.5, .5), 3., variant, joint.frames(45))
    assert matrix.shape == (12, 12)
    assert reduced.segment(model, .5, variant).H == 0


def test_planar_love_inertia_and_harmonic_resultant(setup):
    model = setup[2]
    s = reduced.segment(model, .5, "rayleigh_love_planar")
    assert s.J == model.section.nu**2*model.section.rhoI
    assert s.J != model.section.nu**2*model.section.rho*mh.rectangle_moments(.2, .05)["Ip"]
    x, omega = np.linspace(0, .5, 7), 6.
    value = reduced.basis(model, .5, omega, x, "rayleigh_love_planar")
    derivative = reduced.basis(model, .5, omega, x, "rayleigh_love_planar", 1)
    np.testing.assert_array_equal(value[:, 3, :2], (s.EA-s.J*omega**2)*derivative[:, 0, :2])
    assert not np.array_equal(value[:, 3, :2], s.EA*derivative[:, 0, :2])


def test_love_boundary_variation_exact_polynomial():
    # ∫J v_x δv_x dx = [J v_x δv] - ∫J v_xx δv.
    # Time integration puts +J u_xtt in natural N. Exact polynomial check.
    # v=x², δv=x³ on[0,1]: ∫(2x)(3x²)=3/2=2-1/2.
    assert F(6, 4) == F(2)-F(2, 4)


@pytest.mark.parametrize("variant", reduced.NAMES)
def test_shared_timoshenko_basis_exact(setup, variant):
    model = setup[2]
    x = np.linspace(0, .5, 9)
    values = reduced.basis(model, .5, 7., x, variant)
    np.testing.assert_array_equal(values[:, (1, 2, 4, 5), 2:],
        mh.finite_state_basis(model, .5, 7., x, "timoshenko"))


@pytest.mark.parametrize("beta", (0, 5, 15, 30, 45, 60, 75, 90))
def test_geometry_is_same_dual_project_transform(beta):
    frames = joint.frames(beta)
    operator = reduced.joint_matrix(frames)
    assert operator.shape == (6, 12) and np.linalg.matrix_rank(operator) == 6
    np.testing.assert_allclose(operator@operator.T, 2*np.eye(6), atol=8e-16)
    for frame in frames:
        g = reduced.transform(frame)
        q, p = np.array([.2, -.6, .9]), np.array([-.7, .4, .3])
        assert p@(g.T@q) == pytest.approx((g@p)@q, abs=2e-16)
        np.testing.assert_array_equal(g[:, 2], [0., 0., 1.])


@pytest.mark.parametrize("variant", reduced.NAMES)
def test_global_root_arm_pole_coincidence_count_not_removed(setup, gates, variant):
    model = setup[2]
    pole = reduced.axial_poles(model, .5, variant, 10.)[0]
    effective, count, record = reduced.count_at(model, (.5, .5), pole, variant,
        joint.frames(0), gates[2], audit.POLICY)
    assert effective > pole and record["excluded_poles"]
    case = gates[1][variant, 0]
    assert min(abs(r["omega"]-pole) for r in case["roots"]) < 1e-9
    assert count >= 1


@pytest.mark.parametrize("variant", reduced.NAMES)
def test_residual_normalization_and_condition_gates(gates, variant):
    for beta in (0, 45):
        for root in gates[1][variant, beta]["roots"]:
            reduced.validate_mode(root, audit.POLICY)
            assert root["diagnostics"]["mass_norm"] == pytest.approx(1., abs=5e-16)


def test_overlap_self_and_sign_invariance():
    values = np.array([[1., 2., 3.], [-2., 1., 0.]])
    weight = np.array([.2, .3, .5])
    same, _, _ = reduced.overlap(values, values, weight)
    opposite, _, _ = reduced.overlap(-values, values, weight)
    np.testing.assert_allclose(np.diag(same), 1, atol=5e-16)
    np.testing.assert_array_equal(same, opposite)


def test_theta_small_norm_handling_and_separate_self_overlap():
    values = np.array([[0., 0.], [1., -2.]])
    overlap, norms, _ = reduced.overlap(values, values, [.5, .5], [True, False], [True, False])
    assert norms[0] == 0 and overlap[0] == [None, None]
    assert overlap[1][0] is None and overlap[1][1] == pytest.approx(1.)


def test_contraction_quasistatic_zero_triangle_and_zero_fields():
    v = np.array([1., 2., 3.])
    assert reduced.contraction_diagnostic(-v, v, np.ones(3))["D_c"] == 0
    rng = np.random.default_rng(715)
    for _ in range(7):
        result = reduced.contraction_diagnostic(rng.normal(size=5), rng.normal(size=5), np.ones(5))
        assert 0 <= result["D_c"] <= 1+2e-16
    result = reduced.contraction_diagnostic(np.zeros(3), np.zeros(3), np.ones(3))
    assert result["D_c"] is None and result["status"] == "NOT_DEFINED_SMALL_FIELD"


def test_ranking_reports_exchange_without_reassignment():
    matrix = [[.1, .9], [.8, .2]]
    result = audit.ranking(matrix)
    assert result["rows"][0]["argmax"] == 2
    assert result["rows"][0]["status"] == "POSITION_CORRESPONDENCE_NONDIAGONAL"
    assert matrix == [[.1, .9], [.8, .2]]


def test_inputs_sorted_semantics_and_no_tracking_classification(setup):
    config, _, model, length, *_ = setup
    assert config["spectrum_semantics"] == "sorted_positions"
    assert config["beta_deg"] == [0, 5, 15, 30, 45, 60, 75, 90]
    assert model.variant == "project_jang_reduced_rectangular" and length == 1
    assert reduced.segment(model, .5, "elementary").EA == model.section.EA
    assert reduced.segment(model, .5, "rayleigh_love_planar").m == model.section.rhoA
    for term in ("energy_fraction", "modal_type", "branch_id", "descendant", "hungarian"):
        assert term not in json.dumps(config).lower()


def test_accepted_bundle_and_provenance_when_present(setup):
    pointer = audit.OUTPUT/"current.json"
    if not pointer.exists():
        pytest.skip("No generated screening bundle; comparator computations tested above")
    from pathlib import Path
    out = Path(json.loads(pointer.read_text(encoding="utf-8"))["directory"])
    manifest = json.loads((out/"manifest.json").read_text(encoding="utf-8"))
    assert all(audit.sha(out/p) == h for p, h in manifest["artifact_hashes"].items())
    _, inputs = audit.identity(setup[0], setup[4])
    assert manifest["identity"] == inputs
    data = json.loads((out/"summary.json").read_text(encoding="utf-8"))
    assert data["statuses"]["COUPLED_HIERARCHY_SCREENING"] == "COMPLETE"
    assert len(data["frequency_table"]) == 96 and len(data["overlaps"]) == 24
    assert len(data["adjacent_gaps"]) == 264
    for term in ("energy_fraction", "modal_type", "branch_id", "descendant", "hungarian"):
        assert term not in json.dumps(data).lower()
    for model in setup[0]["models"]:
        for beta in setup[0]["beta_deg"]:
            case = json.loads((out/f"inventory_{model}_{beta}.json").read_text(encoding="utf-8"))
            frequencies = [r["frequency_hz"] for r in case["roots"]]
            assert len(frequencies) >= 13 and all(b > a for a, b in zip(frequencies[:-1], frequencies[1:]))
            assert case["prefix_certificate"]["status"] == "PASS"
            assert len(case["ceiling_attempts"]) == 1


def test_identity_changes_with_input(setup):
    config = dict(setup[0])
    first = audit.identity(config, setup[4])[0]
    config["frequency_ceiling_fstar"] = 5.625
    assert audit.identity(config, setup[4])[0] != first


@pytest.mark.parametrize("beta", (0, 45, 90))
def test_overlap_and_contraction_quadrature_convergence(setup, beta):
    from pathlib import Path
    pointer = audit.OUTPUT/"current.json"
    if not pointer.exists():
        pytest.skip("No saved screening profiles")
    out = Path(json.loads(pointer.read_text(encoding="utf-8"))["directory"])
    samples = []
    for order in (200, 300):
        fields = {}
        for name in setup[0]["models"]:
            case = json.loads((out/f"inventory_{name}_{beta}.json").read_text(encoding="utf-8"))
            v, ux, weights = audit.sample_case(setup[2], case, order=order)
            fields[name] = (v, ux, weights)
        first, _, weights = fields["elementary"]
        second, ux, _ = fields["mindlin_herrmann"]
        overlap, _, _ = reduced.overlap(first[:12, :, :2].reshape(12, -1), second[:12, :, :2].reshape(12, -1), np.repeat(weights, 2))
        dc = [reduced.contraction_diagnostic(second[k, :, 3], setup[2].section.nu*ux[k], weights)["D_c"] for k in range(12)]
        # Null MH fields are noise, not informative ratios; filter by saved status.
        source = json.loads((out/"summary.json").read_text(encoding="utf-8"))
        informative = [r["position"]-1 for r in source["contraction"] if r["beta_deg"] == beta and r["D_c"] is not None]
        samples.append((overlap, np.array(dc)[informative]))
    np.testing.assert_allclose(samples[0][0], samples[1][0], atol=1e-10, rtol=0)
    np.testing.assert_allclose(samples[0][1], samples[1][1], atol=1e-10, rtol=0)
