"""Fixed-reference normalization and production paired-arm map gates."""
from fractions import Fraction as F
import json
import math
from pathlib import Path

import numpy as np
import pytest
from scripts.analysis import verify_mindlin_herrmann_timoshenko_lambda_beta_large_checks as audit
from scripts.lib import mindlin_herrmann_timoshenko_joint as joint


@pytest.fixture(scope="module")
def setup():
    return audit.check_inputs()


@pytest.fixture(scope="module")
def saved(setup):
    pointer = audit.OUTPUT/"current.json"
    if not pointer.exists():
        pytest.skip("Run the map command for generated local evidence")
    out = Path(json.loads(pointer.read_text(encoding="utf-8"))["directory"])
    manifest = json.loads((out/"manifest.json").read_text(encoding="utf-8"))
    assert manifest["identity"] == audit.identity(setup[0],setup[2])[1]
    assert all(audit.sha(out/n)==h for n,h in manifest["artifact_hashes"].items())
    return out,json.loads((out/"result.json").read_text(encoding="utf-8"))


def test_canonical_Lambda_exact_and_not_old_fstar(setup):
    A,I,l = F(1,100),F(1,5)*F(1,20)**3/12,F(1,2)
    assert A*l**4/I == 300
    assert audit.lambda_factor(setup[0]) == pytest.approx(300,rel=3e-15)
    f = np.array([.05,.5,1.,2.])
    result = audit.lambda_from_frequency(f,setup[0])
    np.testing.assert_allclose(result**4,300*(2*math.pi*f)**2,rtol=1e-15)
    assert not np.array_equal(result,f)
    np.testing.assert_allclose(result**2,2*math.pi*f*math.sqrt(300),rtol=5e-16)


@pytest.mark.parametrize("mu",(0,.25,.5,-.25,-.5))
def test_length_total_mass_and_geometry(setup,mu):
    models,lengths,description = audit.geometry(setup[0],"length_asymmetry",mu)
    assert lengths == (.5*(1-mu),.5*(1+mu))
    assert sum(lengths)==1
    assert description["total_mass"] == pytest.approx(.01,abs=3e-18)
    assert all(m.section.K==5/6 and m.mh_shear_factor==5/6 and m.mh_inertia_factor==1 for m in models)


@pytest.mark.parametrize("contrast",(0,.2,.4,-.2,-.4))
def test_sections_mass_and_preset_unchanged(setup,contrast):
    models,lengths,description = audit.geometry(setup[0],"thickness_contrast",contrast)
    assert lengths == (.5,.5)
    assert description["heights"] == (.05*(1-contrast),.05*(1+contrast))
    assert min(description["heights"])>0
    assert description["total_mass"] == pytest.approx(.01,abs=3e-18)
    assert all(m.variant=="project_jang_reduced_rectangular" and m.section.K==5/6 for m in models)


@pytest.mark.parametrize("beta",(0,22.5,45,67.5,90))
def test_composition_is_exact_existing_boundary_when_models_equal(setup,beta):
    models,lengths,_ = audit.geometry(setup[0],"length_asymmetry",.25)
    np.testing.assert_array_equal(audit.boundary(models,lengths,7.,beta),joint.frame_boundary_matrix(models[0],lengths,7.,beta_deg=beta))


@pytest.mark.parametrize("beta",(0,22.5,45,67.5,90))
def test_arm_swap_is_bisector_reflection_with_signed_theta(beta):
    frames = joint.frames(beta)
    angle = math.radians(beta)
    reflection = np.array([[-math.cos(angle),-math.sin(angle)],[-math.sin(angle),math.cos(angle)]])
    for first,second in zip(frames,frames[::-1]):
        np.testing.assert_allclose(reflection@first.t,second.t,atol=3e-16)
        np.testing.assert_allclose(-reflection@first.n,second.n,atol=3e-16)
    np.testing.assert_array_equal(joint.MIRROR_STATE,[1,1,-1,-1,1,1,-1,-1])


def test_baseline_frequency_and_Lambda_regressions(saved):
    result = saved[1]
    assert [r["beta_deg"] for r in result["baseline_regression"]]==[0,5,15,30,45,60,75,90]
    assert max(v for r in result["baseline_regression"] for v in r["frequency_relative_differences"])<1e-10
    assert max(v for r in result["baseline_regression"] for v in r["Lambda_relative_differences"])<1e-10


@pytest.mark.parametrize("mu",(0,.25,.5))
def test_length_beta0_direct_collapse_and_full_profiles(saved,mu):
    check = next(c for c in saved[1]["reference_checks"] if c["family"]=="length_asymmetry" and c["geometry_parameter"]==mu)
    assert check["comparison"]["status"]=="PASS"
    assert check["comparison"]["direct_homogeneous_max_frequency_relative_difference"]<1e-10
    assert max(r["kinematic_L2_relative"] for r in check["comparison"]["rows"])<5e-8
    assert max(r["resultant_L2_relative"] for r in check["comparison"]["rows"])<5e-8


@pytest.mark.parametrize("contrast",(.2,.4))
def test_independent_straight_stepped_match(saved,contrast):
    check = next(c for c in saved[1]["reference_checks"] if c["family"]=="thickness_contrast" and c["geometry_parameter"]==contrast)
    assert "no angle transformation" in check["reference"]["method"]
    assert max(r["frequency_relative_difference"] for r in check["comparison"]["rows"])<5e-10
    assert max(r["resultant_L2_relative"] for r in check["comparison"]["rows"])<5e-8


@pytest.mark.parametrize("family,parameter",(("length_asymmetry",.25),("length_asymmetry",.5),("thickness_contrast",.2),("thickness_contrast",.4)))
def test_sparse_swap_frequencies_shapes_and_resultants(saved,family,parameter):
    checks = [c for c in saved[1]["swap_checks"] if c["family"]==family and c["geometry_parameter"]==parameter]
    assert [c["beta_deg"] for c in checks]==[0,22.5,45,67.5,90]
    assert all(c["comparison"]["status"]=="PASS" for c in checks)
    assert max(r["frequency_relative_difference"] for c in checks for r in c["comparison"]["rows"])<5e-10
    assert max(r["resultant_L2_relative"] for c in checks for r in c["comparison"]["rows"])<5e-8


def test_all_map_points_and_unequal_section_quality(saved):
    result = saved[1]
    assert len(result["series"])==5
    for series in result["series"].values():
        assert len(series["cases"])==37
        for case in series["cases"]:
            assert len(case["roots"])==13
            assert case["search"]["prefix_certificate"]=="PASS" and case["search"]["guard_right_count"]==13
            assert case["mass_gram_max_error"]<5e-7
            assert case["search"]["failed_intervals"]==[]
            frequencies = [r["omega"] for r in case["roots"]]
            assert all(b>a for a,b in zip(frequencies[:-1],frequencies[1:]))
            for r in case["roots"]:
                d = r["diagnostics"]
                assert max(d["joint_residual_scaled"].values())<1e-9
                assert d["clamp_scaled_residual"]<1e-9 and d["equation_scaled_residual"]<1e-9
                assert d["mass_norm"]==pytest.approx(1.,abs=5e-16)
                assert d["nonzero_singular_condition"]<1e8
    assert result["performance"]["local_continuations"]==200
    assert result["performance"]["full_scans"]==5
    assert result["performance"]["fallback_full_scans"]==0


def test_curvature_flags_have_independent_same_point_review(saved):
    flags = [f for s in saved[1]["series"].values() for f in s["anomaly_flags"]]
    for flag in flags:
        assert flag["status"]=="REVIEWED_INDEPENDENT_QR_COUNT_SVD_RESIDUAL_PASS"
        assert flag["independent_same_point_review"]["relative_frequency_difference"]<1e-10


def test_no_tracking_or_energy_classification_metadata(saved):
    text = json.dumps(saved[1]).lower()
    for forbidden in ("energy_fraction","modal_type","branch_id","descendant_id","hungarian","tracking_assignment"):
        assert forbidden not in text


def test_normalization_change_invalidates_cache_identity(setup):
    first = audit.identity(setup[0],setup[2])[0]
    config = dict(setup[0]);config["reference"]={**config["reference"],"h_ref":.06}
    assert audit.identity(config,setup[2])[0]!=first


def test_plot_only_has_zero_root_calls(saved,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError("Plot-only called a solver")
    monkeypatch.setattr(audit,"compute",forbidden)
    monkeypatch.setattr(audit,"solve_case",forbidden)
    monkeypatch.setattr(audit,"localized_roots",forbidden)
    monkeypatch.setattr(audit,"check_inputs",forbidden)
    assert audit.main(["--plot-only",str(saved[0])])==0
