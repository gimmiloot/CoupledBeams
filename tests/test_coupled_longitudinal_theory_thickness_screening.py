"""Exact coefficient scaling and the bounded, independently sorted thickness gate."""
from fractions import Fraction as F
import json
from pathlib import Path

import numpy as np
import pytest

from scripts.analysis import screen_coupled_longitudinal_theory_thickness as audit
from scripts.analysis import screen_coupled_longitudinal_theory_hierarchy as hierarchy
from scripts.lib import coupled_longitudinal_comparators as reduced
from scripts.lib import mindlin_herrmann_longitudinal as mh


@pytest.fixture(scope="module")
def setup():
    return audit.check_inputs()


@pytest.fixture(scope="module")
def saved(setup):
    pointer = audit.OUTPUT/"current.json"
    if not pointer.exists():
        pytest.skip("Run the documented thickness CLI to create local spectral evidence")
    out = Path(json.loads(pointer.read_text(encoding="utf-8"))["directory"])
    manifest = json.loads((out/"manifest.json").read_text(encoding="utf-8"))
    assert all(audit.sha(out/name) == digest for name, digest in manifest["artifact_hashes"].items())
    _, identity = audit.identity(*setup[:3])
    assert manifest["identity"] == identity
    return out, json.loads((out/"summary.json").read_text(encoding="utf-8"))


@pytest.mark.parametrize("sh", (F(1), F(5, 4), F(3, 2), F(7, 4), F(2)))
def test_exact_rectangle_and_quadratic_scale(sh):
    b, h0, arm = F(1, 5), F(1, 20), F(1, 2)
    h = h0*sh
    moments = mh.rectangle_moments(b, h)
    assert moments["A"] == b*h
    assert moments["Iy"] == b*h**3/12
    assert moments["Iy"]/moments["A"] == h**2/12
    q = moments["Iy"]/moments["A"]/arm**2
    assert q/(h0**2/12/arm**2) == sh**2


@pytest.mark.parametrize("h", (F(1, 20), F(1, 16), F(3, 40), F(7, 80), F(1, 10)))
def test_width_cancellation_and_normalized_coefficients_exact(h):
    E, rho, nu, kappa, b = F(1), F(1), F(3, 10), F(5, 6), F(1, 5)
    ratios = audit.coefficient_ratios(E, rho, nu, kappa, b, h)
    assert ratios == audit.coefficient_ratios(E, rho, nu, kappa, b*3, h)
    ell2 = h**2/12
    assert ratios["EA_over_m"] == E/rho
    assert ratios["J_over_m"] == nu**2*ell2
    assert ratios["j_over_m"] == ratios["r_over_m"] == ell2
    assert ratios["H_over_C"] == kappa*(1-nu)*ell2/2
    assert ratios["B_over_m"] == E/rho*ell2
    assert ratios["S_over_r"] == kappa*E/(2*(1+nu)*rho*ell2)
    assert ratios["static_MH_layer_squared"] == kappa*h**2/(24*(1+nu))


@pytest.mark.parametrize("h", (.05, .0625, .075, .0875, .1))
def test_same_bending_inputs_and_basis_all_theories(setup, h):
    model = audit.make_model(setup[0], h)
    x = np.linspace(0, .5, 5)
    exact = mh.finite_state_basis(model, .5, 8., x, "timoshenko")
    for name in reduced.NAMES:
        np.testing.assert_array_equal(reduced.basis(model, .5, 8., x, name)[:, (1, 2, 4, 5), 2:], exact)
    assert model.section.K == 5/6 == model.mh_shear_factor
    assert model.mh_inertia_factor == model.tim_rotary_factor == 1
    assert reduced.segment(model, .5, "rayleigh_love_planar").J == model.section.nu**2*model.section.rhoI
    assert reduced.segment(model, .5, "elementary").H == 0


def test_mass_grows_without_inter_geometry_preservation(setup):
    baseline = audit.make_model(setup[0], .05)
    thick = audit.make_model(setup[0], .1)
    assert thick.section.rhoA/baseline.section.rhoA == 2
    assert thick.section.rhoI/baseline.section.rhoI == 8


def test_scaling_audit_and_domain_metadata(setup):
    record = audit.scaling_audit(setup[0])
    assert record["status"] == "PASS"
    for r in record["records"]:
        assert r["bending_inputs_by_model"]["elementary"] == r["bending_inputs_by_model"]["mindlin_herrmann"]
        assert r["bending_inputs_by_model"]["rayleigh_love_planar"] == r["bending_inputs_by_model"]["mindlin_herrmann"]
    last = record["records"][-1]
    assert last["L_arm_over_h"] == 5 and last["b_over_h"] == 2
    assert last["Tim_cutoff_fstar"] < setup[0]["frequency_ceiling_fstar"]


@pytest.mark.parametrize("sh", (1., 1.5, 2.))
@pytest.mark.parametrize("name", ("elementary", "rayleigh_love_planar", "mindlin_herrmann"))
def test_fresh_direct_beta0_profiles_and_joint_transparency(setup, saved, sh, name):
    out, _ = saved
    h = setup[0]["h"][setup[0]["s_h"].index(sh)]
    model = audit.make_model(setup[0], h)
    case = json.loads((out/f"inventory_{sh:g}_{name}_0.json").read_text(encoding="utf-8"))
    gate = audit.direct_spot(model, case)
    assert gate["status"] == "PASS" and len(gate["rows"]) == 13
    assert max(r["kinematic_profile_L2_relative"] for r in gate["rows"]) < 5e-8
    assert max(r["relative_frequency_difference"] for r in gate["rows"]) < 1e-10


def test_every_thickness_direct_checks_saved(saved):
    out, _ = saved
    checks = json.loads((out/"direct_beta0_checks.json").read_text(encoding="utf-8"))
    assert len(checks) == 15
    assert all(c["status"] == "PASS" and len(c["rows"]) == 13 for c in checks)


def test_all45_prefixes_and_quality_gates(saved):
    out, data = saved
    files = list(out.glob("inventory_*.json"))
    assert len(files) == 45
    for file in files:
        case = json.loads(file.read_text(encoding="utf-8"))
        assert len(case["roots"]) >= 13
        assert case["prefix_certificate"]["status"] == case["search"]["status"] == "PASS"
        assert case["search"]["lower_count"] == 0
        assert case["search"]["failed_intervals"] == []
        frequencies = [r["frequency_hz"] for r in case["roots"]]
        assert all(b > a for a, b in zip(frequencies[:-1], frequencies[1:]))
        for r in case["roots"]:
            reduced.validate_mode(r, audit.POLICY)
        assert len(case["ceiling_attempts"]) == 1
    assert data["statuses"]["COUPLED_THICKNESS_SCREENING"] == "COMPLETE"


def test_baseline_hierarchy_roots_unchanged(setup, saved):
    out, _ = saved
    for beta in setup[0]["beta_deg"]:
        for name in setup[0]["models"]:
            baseline = json.loads((setup[3]/f"inventory_{name}_{beta}.json").read_text(encoding="utf-8"))
            fresh = json.loads((out/f"inventory_1_{name}_{beta}.json").read_text(encoding="utf-8"))
            assert fresh["roots"] == baseline["roots"]


def test_geometric_overlap_self_sign_and_dc_zero():
    v = np.array([[1., 2., -.5], [.2, -.7, .9]])
    weights = np.array([.5, .2, .3])
    first, _, _ = reduced.overlap(v, v, weights)
    second, _, _ = reduced.overlap(v, -v, weights)
    np.testing.assert_allclose(np.diag(first), 1, atol=5e-16)
    np.testing.assert_array_equal(first, second)
    assert reduced.contraction_diagnostic(-v[0], v[0], weights)["D_c"] == 0


def test_fixed_case_semantics_and_no_classification(saved):
    _, data = saved
    assert len(data["frequency_table"]) == 180
    assert len(data["overlaps"]) == 45
    assert len(data["adjacent_gaps"]) == 495
    assert len(data["summary"]) == 15
    for term in ("energy_fraction", "modal_type", "branch_id", "descendant", "hungarian", "tracking_assignment"):
        assert term not in json.dumps(data).lower()
    for row in data["frequency_table"]:
        assert row["abs_delta_E"] == abs(row["delta_E_signed"])
        assert row["abs_delta_RL"] == abs(row["delta_RL_signed"])
        assert row["D_c"] is None or 0 <= row["D_c"] <= 1+1e-14


@pytest.mark.parametrize("beta", (0, 45, 90))
def test_thickest_shape_and_contraction_quadrature_convergence(setup, saved, beta):
    out, data = saved
    model = audit.make_model(setup[0], .1)
    samples = []
    for order in (200, 300):
        fields = {}
        for name in setup[0]["models"]:
            case = json.loads((out/f"inventory_2_{name}_{beta}.json").read_text(encoding="utf-8"))
            fields[name] = hierarchy.sample_case(model, case, order=order)
        first, _, weights = fields["elementary"]
        second, ux, _ = fields["mindlin_herrmann"]
        overlap, _, _ = reduced.overlap(first[:12, :, :2].reshape(12, -1), second[:12, :, :2].reshape(12, -1), np.repeat(weights, 2))
        dc = [reduced.contraction_diagnostic(second[k, :, 3], .3*ux[k], weights)["D_c"] for k in range(12)]
        informative = [r["position"]-1 for r in data["contraction"] if r["s_h"] == 2 and r["beta_deg"] == beta and r["D_c"] is not None]
        samples.append((overlap, np.array(dc)[informative]))
    np.testing.assert_allclose(samples[0][0], samples[1][0], atol=1e-10, rtol=0)
    np.testing.assert_allclose(samples[0][1], samples[1][1], atol=1e-10, rtol=0)


def test_cache_identity_changes_with_config(setup):
    first = audit.identity(*setup[:3])[0]
    config = dict(setup[0])
    config["frequency_ceiling_fstar"] = 5.625
    assert audit.identity(config, setup[1], setup[2])[0] != first
