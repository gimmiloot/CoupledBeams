"""Saved-data FEM-3C checks; no real FEM/ODE/equilibrium jobs."""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pytest

from scripts.lib import nlsp_fem3c_diagnostics as diagnostic


@pytest.fixture
def policy():
    return {"T1_fraction": .25, "points": 201,
        "primary_interpolation": "linear", "diagnostic_interpolation": "PCHIP",
        "baseline_model_discrepancy": 2.824717e-7,
        "temporal_ratio_limit": .25, "spatial_ratio_limit": .25,
        "interpolation_ratio_to_effect_limit": .25,
        "interpolation_ratio_to_baseline_limit": .25,
        "phase_amplitude_fitting": False}


def case(root, evolution=2e-6, *, nonlinear=False, horizon=1.):
    root.mkdir(parents=True)
    time = np.linspace(.01, horizon, 101)
    x = np.linspace(0., 1., 41)
    profile = np.sin(np.pi*x)
    fields = np.zeros((len(time), len(x), 7))
    fields[..., 1] = (0.005 - .001*time[:, None])*profile
    initial = np.zeros((len(x), 7))
    initial[:, 1] = .005*profile
    if nonlinear:
        fields[..., 1] += (-9e-6 + evolution*time[:, None])*profile
        initial[:, 1] -= 9e-6*profile
    velocity = np.zeros((len(time), len(x), 3))
    velocity[..., 1] = -.001*profile
    np.savez(root/"section_history.npz", time=time, x=x, fields=fields,
        translation_velocities=velocity)
    np.savez(root/"initial_sections.npz", x=x, fields=initial,
        source_static_end_time=1., not_a_native_dynamic_zero_frame=True)


def test_preregistered_comparison_grid(policy):
    grid = diagnostic.comparison_grid(policy, 4.)
    assert len(grid) == 201
    assert grid[0] == 0 and grid[-1] == 1.
    for key, value in (("points", 202), ("T1_fraction", .5),
            ("baseline_model_discrepancy", 3e-7), ("temporal_ratio_limit", .3),
            ("spatial_ratio_limit", .3), ("phase_amplitude_fitting", True)):
        changed = dict(policy, **{key: value})
        with pytest.raises(ValueError, match="Preregistered"):
            diagnostic.comparison_grid(changed, 4.)


def test_seven_field_mapping_and_planar_qualification():
    active = np.arange(24., dtype=float).reshape(2, 3, 4)
    full = diagnostic.seven_field_one_d(active)
    np.testing.assert_array_equal(full[..., [0, 1, 5, 6]], active)
    np.testing.assert_array_equal(full[..., [2, 3, 4]], 0.)
    assert diagnostic.CANONICAL_FIELDS[-1] == "c_eff_diagnostic"
    with pytest.raises(ValueError): diagnostic.seven_field_one_d(np.zeros((4, 7)))


def test_static_zero_and_interpolated_native_distinction(tmp_path):
    case(tmp_path/"case", nonlinear=True)
    data = diagnostic.transfer_case(tmp_path/"case", np.linspace(0, 1, 201))
    assert data["metadata"]["zero_source"] == "confirmed_STATIC_preload_not_native_DYNAMIC_frame"
    assert data["metadata"]["positive_comparison_values"] == "postprocessing_interpolated"
    np.testing.assert_array_equal(data["fields"][0], data["initial_fields"])
    np.testing.assert_array_equal(data["translation_velocities"][0], 0.)
    assert data["native_time"][0] > 0
    with pytest.raises(ValueError, match="coverage"):
        diagnostic.transfer_case(tmp_path/"case", np.linspace(0, 1.01, 201))
    with pytest.raises(ValueError, match="Only preregistered"):
        diagnostic.transfer_case(tmp_path/"case", np.linspace(0, 1, 201), "cubic")


def test_invalid_initial_provenance_is_not_a_native_dynamic_frame(tmp_path):
    root = tmp_path/"case"
    case(root)
    data = diagnostic.load_arrays(root/"initial_sections.npz")
    data["not_a_native_dynamic_zero_frame"] = False
    np.savez(root/"initial_sections.npz", **data)
    with pytest.raises(ValueError, match="invalid"):
        diagnostic.transfer_case(root, np.linspace(0, 1, 201))


def test_evolution_subtracts_static_correction_preserving_signs(tmp_path):
    case(tmp_path/"L")
    case(tmp_path/"N", nonlinear=True)
    time = np.linspace(0, 1, 201)
    pair = diagnostic.correction_pair(diagnostic.transfer_case(tmp_path/"L", time),
        diagnostic.transfer_case(tmp_path/"N", time))
    assert pair["initial_correction"][20, 1] == pytest.approx(-9e-6)
    assert pair["correction"][-1, 20, 1] == pytest.approx(-7e-6)
    assert pair["evolution"][-1, 20, 1] == pytest.approx(2e-6)
    np.testing.assert_array_equal(pair["evolution"][0], 0.)


def test_interpolation_comparability_and_zero_effect(policy):
    assert diagnostic.interpolation_comparability({"old": 1e-10}, 1e-8, 2e-8, policy)["status"] == "PASS"
    result = diagnostic.interpolation_comparability({"old": 4e-9}, 1e-8, 2e-8, policy)
    assert result["status"] == "INTERPOLATION_UNRESOLVED"
    assert result["diagnostic_comparability_only"]
    assert not result["continuum_error_bound"]
    assert diagnostic.interpolation_comparability({"old": 0.}, 0., 1e-8, policy)["status"] == "PASS"
    assert diagnostic.interpolation_comparability({"old": 1e-30}, 0., 1e-8, policy)["status"] == "INTERPOLATION_UNRESOLVED"


def test_complete_robustness_same_fixed_scales_and_source_preservation(tmp_path, policy, monkeypatch):
    roots = {"existing_medium": tmp_path/"old", "medium_refined_time": tmp_path/"medium",
             "fine_refined_time": tmp_path/"fine"}
    for root, evolution in zip(roots.values(), (2e-6, 2.01e-6, 2.03e-6)):
        for kind in ("linear", "nonlinear"):
            case(root/kind, evolution=evolution, nonlinear=kind == "nonlinear")
    original = {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in roots["existing_medium"].rglob("*") if p.is_file()}
    calls = []
    def reference(time):
        calls.append(time.copy())
        L = diagnostic.transfer_case(roots["existing_medium"]/"linear", time)
        N = diagnostic.transfer_case(roots["existing_medium"]/"nonlinear", time)
        return {"times": time.copy(), "x": L["x"],
            "linear_fields": L["fields"][..., [0, 1, 5, 6]],
            "nonlinear_fields": N["fields"][..., [0, 1, 5, 6]],
            "initial_linear_fields": L["initial_fields"][..., [0, 1, 5, 6]],
            "initial_nonlinear_fields": N["initial_fields"][..., [0, 1, 5, 6]]}
    output = tmp_path/"output"
    output.mkdir()
    result = diagnostic.robustness(output, {"comparison": policy, "T1": 4.},
        case_roots=roots, parent_bundle=roots["existing_medium"], one_d_evaluator=reference)
    assert len(calls) == 1
    assert result["full_period_numerical_robustness_gate"]
    assert result["refinement"]["temporal"]["evolution_w_max"] == pytest.approx(1e-8)
    assert result["refinement"]["spatial"]["evolution_w_max"] == pytest.approx(2e-8)
    assert result["refinement"]["temporal"]["ratio_to_fixed_baseline_model_discrepancy"] == pytest.approx(1e-8/2.824717e-7)
    for quantity in ("linear", "nonlinear", "correction", "evolution"):
        scales = [result["resolutions"][label]["model_comparisons"][quantity]["w"]["fixed_all_resolutions_characteristic_scale"] for label in roots]
        assert scales[0] == scales[1] == scales[2]
    assert result["scientific_calls"] == {"CCX": 0, "Gmsh": 0, "nonlinear_ODE": 0,
        "eigenanalysis": 0, "static_equilibrium": 0}
    assert {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in original} == original
    assert all((output/name).is_file() for name in ("robustness_comparison.json", "robustness_comparison.npz", "robustness_summary.csv"))


def test_robustness_partial_does_not_relax_threshold(tmp_path, policy):
    result = diagnostic.interpolation_comparability({"medium": 2e-8}, 1e-7, 1e-9, policy)
    assert result["status"] == "INTERPOLATION_UNRESOLVED"
    assert result["limit_from_effect"] == 2.5e-10


def test_full_period_saved_fields_and_native_csv_are_distinguished(tmp_path):
    parent = tmp_path/"old"
    output = tmp_path/"new"
    for kind in ("linear", "nonlinear"):
        case(parent/"cases"/kind, nonlinear=kind == "nonlinear", horizon=1.)
        case(output/"cases"/"full_period_medium"/kind, nonlinear=kind == "nonlinear", horizon=4.)
    item = {"config": {"omega1": 2*np.pi/4.}, "validation_config": {
        "parent_completed": {"bundle": str(parent)},
        "full_period_comparison": {"points": 401, "primary_interpolation": "linear",
            "diagnostic_interpolation": "PCHIP", "snapshot_T1_fractions": [0., .25, .5, .75, 1.],
            "phase_amplitude_fitting": False},
        "seven_field_observations": {"u": .25, "w": .5, "v": .5, "Phi": .25, "psi": .25,
            "theta": .25, "c": .25, "inactive_one_d_fields": ["v", "Phi", "psi"],
            "three_d_c": "effective_contraction_proxy_not_identical_generalized_coordinate"}}}
    def reference(time):
        root = output/"cases"/"full_period_medium"
        L, N = (diagnostic.transfer_case(root/kind, time) for kind in ("linear", "nonlinear"))
        return {"times": time.copy(), "x": L["x"],
            "linear_fields": L["fields"][..., [0, 1, 5, 6]],
            "nonlinear_fields": N["fields"][..., [0, 1, 5, 6]],
            "linear_velocities": np.zeros((len(time), len(L["x"]), 4)),
            "nonlinear_velocities": np.zeros((len(time), len(L["x"]), 4)),
            "initial_linear_fields": L["initial_fields"][..., [0, 1, 5, 6]],
            "initial_nonlinear_fields": N["initial_fields"][..., [0, 1, 5, 6]]}
    result = diagnostic.full_period(output, item, {"p48_spatial_control": {"status": "PARTIAL"}}, one_d_evaluator=reference)
    assert result["physical_interval"] == [0., 4.]
    assert result["comparison_points"] == 401
    assert result["quarter_period_robustness_not_extended_to_full_period"]
    assert not result["nonlinear_periodic_orbit_assumed_or_found"]
    assert not result["general_seven_field_nonlinear_validation"]
    assert result["one_d_spatial_status"] == "PARTIAL"
    assert result["canonical_one_d_field_order"] == ("u", "w", "v", "Phi", "psi", "theta", "c")
    assert result["recovered_three_d_field_order"][-1] == "c_eff_diagnostic"
    interpolation = result["full_period_interpolation_uncertainty"]
    assert interpolation["quarter_period_robustness_not_extended"]
    assert interpolation["status"] == "OBSERVED_DIAGNOSTIC_NOT_FULL_PERIOD_CONVERGENCE_CERTIFICATION"
    assert not interpolation["continuum_error_bound_claimed"]
    for quantity in ("correction", "evolution"):
        row = interpolation["w"][quantity]
        assert row["absolute_max"] == result["interpolation_diagnostic_differences"][quantity]["w"]["absolute_max"]
        assert row["fixed_full_period_characteristic_scale"] == result["model_comparisons"][quantity]["w"]["fixed_all_resolutions_characteristic_scale"]
        assert not row["scientific_acceptance_threshold_assigned"]
    for field in ("v", "Phi", "psi"):
        assert result["inactive_one_d_fields"][field]["one_d_identically_zero_by_planar_subspace"]
        assert not result["inactive_one_d_fields"][field]["out_of_plane_stability_verified"]
    assert result["observations"]["c_eff_diagnostic"]["material_x"] == .25
    assert result["energy"]["status"] == "PARTIAL"
    assert result["energy"]["independent_internal_StVK_energy"] == "NOT_RUN"
    raw = (output/"full_period_native_observations.csv").read_text(encoding="utf8")
    transferred = (output/"full_period_observations.csv").read_text(encoding="utf8")
    assert "actual_native_DYNAMIC_sample" in raw
    assert "confirmed_STATIC_preload_not_native_DYNAMIC_frame" in raw
    assert "postprocessing_linear_interpolation" in transferred
    assert "actual_native_DYNAMIC_sample" not in transferred


def test_full_period_requires_frozen_sampling_and_observations():
    c = {"full_period_comparison": {"points": 800}}
    with pytest.raises(ValueError, match="Preregistered"):
        diagnostic._full_grid(c, 1.)


def test_independent_kinetic_saved_increment_matching(tmp_path):
    data = np.array([1e-8, 2e-8])
    np.savez(tmp_path/"independent_kinetic_energy.npz", time=[.1, .2], kinetic_energy=data)
    np.savez(tmp_path/"section_history.npz", time=[.1, .2], increments=[1, 3])
    diagnostic.write_json(tmp_path/"energy.json", {"records": [
        {"step": 2, "increment": 1, "kinetic_energy": 1.00001e-8},
        {"step": 2, "increment": 2, "kinetic_energy": 1.5e-8},
        {"step": 2, "increment": 3, "kinetic_energy": 2e-8}]})
    result = diagnostic._saved_kinetic_comparison(tmp_path)
    assert result["maximum_absolute_difference"] == pytest.approx(1e-13)
    assert result["maximum_relative_difference_fixed_K_scale"] == pytest.approx(5e-6)
    assert not result["native_internal_energy_reference_corrected"]
