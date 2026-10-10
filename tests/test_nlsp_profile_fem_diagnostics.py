"""Synthetic saved-field audit controls; no real FEM or ODE jobs."""
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from scripts.lib import nlsp_profile_fem_diagnostics as audit


def tetrahedron():
    corners = np.array(((0., 0., 0.), (1., 0., 0.), (0., .1, 0.), (0., 0., .2)))
    mids = np.array([(corners[i]+corners[j])/2 for i, j in audit.fem1.TET10_EDGES])
    xyz = np.vstack((corners, mids))
    mesh = SimpleNamespace(nodes={i+1: p for i, p in enumerate(xyz)},
                           solid_elements={1: tuple(range(1, 11))})
    return mesh, xyz


def volume_samples(count=41, F=None, thickness_quadratic=0., weights_unequal=False):
    rows = []
    for k in range(count):
        for d in (-.4, -.15, .15, .4):
            for eta in (-.04, -.01, .03):
                for zeta in (-.08, .02, .08):
                    rows.append(((k+.5+d)/count, eta, zeta))
    xyz = np.array(rows)
    F = np.eye(3) if F is None else np.asarray(F)
    displacement = xyz@(F-np.eye(3)).T
    displacement[:, 1] += thickness_quadratic*xyz[:, 1]**2
    grad = np.broadcast_to(F-np.eye(3), (len(xyz), 3, 3)).copy()
    grad[:, 1, 1] += 2*thickness_quadratic*xyz[:, 1]
    green = .5*(grad+grad.swapaxes(-1, -2)+grad.swapaxes(-1, -2)@grad)
    linear = .5*(grad+grad.swapaxes(-1, -2))
    volume = np.full(len(xyz), .02/len(xyz))
    if weights_unequal:
        volume *= 1.+xyz[:, 1]*8
        volume *= .02/volume.sum()
    deformation = np.eye(3)+grad
    transverse = deformation[:, :, 1:3]
    eig, vec = np.linalg.eigh(transverse.swapaxes(-1, -2)@transverse)
    stretch = (vec*np.sqrt(eig)[:, None, :])@vec.swapaxes(-1, -2)
    prepared = {"length": 1., "thickness": .1, "width": .2, "rho": 1.,
                "reference_volume": .02}
    native = {"global_points": xyz*audit.LOCAL_SIGNS,
        "global_displacement": displacement*audit.LOCAL_SIGNS,
        "local_points": xyz, "local_displacement": displacement,
        "mass_weights": volume, "volume_weights": volume,
        "gradient_local": grad, "linear_local": linear, "green_local": green,
        "point_polar_thickness": stretch[:, 0, 0]-1., "point_polar_width": stretch[:, 1, 1]-1.}
    return prepared, native


@pytest.mark.parametrize("rho", (1., 2.7))
def test_positive_reference_volume_mass_weighting(rho):
    mesh, xyz = tetrahedron()
    prepared = audit.prepare_mesh(mesh, rho)
    assert prepared["reference_volume"] == pytest.approx(.02/6)
    native = audit.native_quadrature_fields(prepared, xyz*0.)
    np.testing.assert_allclose(native["volume_weights"]*rho, native["mass_weights"])
    assert np.all(native["volume_weights"] > 0)


def test_global_to_local_axes_exact_affine_gradient():
    mesh, xyz = tetrahedron()
    local_gradient = np.array(((.002, .013, -.007), (.001, -.003, .005), (.002, .001, .004)))
    global_gradient = local_gradient*audit.LOCAL_SIGNS[:, None]*audit.LOCAL_SIGNS[None, :]
    native = audit.native_quadrature_fields(audit.prepare_mesh(mesh), xyz@global_gradient.T)
    expected = np.broadcast_to(local_gradient, native["gradient_local"].shape)
    np.testing.assert_allclose(native["gradient_local"], expected, atol=1e-14)
    green = .5*(local_gradient+local_gradient.T+local_gradient.T@local_gradient)
    np.testing.assert_allclose(native["green_local"], np.broadcast_to(green, native["green_local"].shape), atol=1e-14)
    assert native["same_gradient_green_polar_relation_max_error"] < 1e-14


@pytest.mark.parametrize("angle", (.02, -.23))
def test_rigid_rotation_zero_finite_strain_not_false_small_strain(angle):
    mesh, xyz = tetrahedron()
    R = Rotation.from_rotvec((0., 0., angle)).as_matrix()
    global_R = R*audit.LOCAL_SIGNS[:, None]*audit.LOCAL_SIGNS[None, :]
    native = audit.native_quadrature_fields(audit.prepare_mesh(mesh), xyz@(global_R-np.eye(3)).T)
    assert np.max(np.abs(native["green_local"])) < 1e-14
    assert np.max(np.abs(native["point_polar_thickness"])) < 1e-14
    np.testing.assert_allclose(native["gradient_local"][..., 1, 1], np.cos(angle)-1., atol=1e-14)
    assert native["minimum_det_F"] == pytest.approx(1., abs=1e-14)


def test_exact_C3D10_quadratic_displacement_gradient_and_bin_moment():
    mesh, xyz = tetrahedron()
    # P2 is exactly representable in the existing synthetic C3D10 element.
    # Local eta=-global Y; its squared displacement maps back to global -Y.
    local = xyz*audit.LOCAL_SIGNS
    local_U = np.zeros_like(xyz)
    local_U[:, 1] = .03*local[:, 1]**2
    native = audit.native_quadrature_fields(audit.prepare_mesh(mesh), local_U*audit.LOCAL_SIGNS)
    expected = .06*native["local_points"][..., 1]
    np.testing.assert_allclose(native["gradient_local"][..., 1, 1], expected, atol=1e-14)
    weight = native["volume_weights"]
    average = np.average(native["gradient_local"][..., 1, 1], weights=weight)
    moment = np.average(native["local_points"][..., 1], weights=weight)
    assert average == pytest.approx(.06*moment, abs=1e-14)
    assert abs(average) > 1e-5  # asymmetric sampling changes the mean despite exact P2.


@pytest.mark.parametrize("count", (21, 41, 81))
def test_uniform_transverse_stretch_exact_raw_and_native(count):
    F = np.diag((1.001, .98, 1.03))
    prepared, native = volume_samples(count, F)
    out = audit.recover_level(prepared, native, count)
    np.testing.assert_allclose(out["c_eff"], -.02, atol=1e-13)
    np.testing.assert_allclose(out["c_small"], -.02, atol=1e-13)
    np.testing.assert_allclose(out["width_effective"], .03, atol=1e-13)
    np.testing.assert_allclose(out["native_columns"]["thickness_small"], -.02, atol=1e-13)
    np.testing.assert_allclose(out["native_columns"]["thickness_green"], -.02+.5*.02**2, atol=1e-13)
    assert np.all(out["fit_rank"] == 8)
    assert out["summary"]["fit_green_polar_relation_error_max"] < 1e-13
    assert out["historical_profile"]["additional_derivative_constraints"] is False
    assert out["interpolation"]["raw_knot_reproduction_max_error"] < 1e-15
    assert out["roughness"]["c_eff"]["raw_extrema_count"] == 0


def test_historical_recovery_exactly_reused():
    prepared, native = volume_samples(F=np.array(((1., -.1, .02), (.01, .99, .03), (.02, .01, 1.01))))
    actual = audit.recover_level(prepared, native)
    original = audit.fem2.fem2_recover_reference_samples(native["global_points"], native["global_displacement"],
        native["mass_weights"], 1., .1, .2)
    np.testing.assert_array_equal(actual["historical_profile"]["fields"], original["fields"])
    np.testing.assert_array_equal(actual["historical_profile"]["x"], original["x"])


def test_uniform_rigid_section_rotation_recovery_no_polar_contraction():
    angle = .19
    R = Rotation.from_rotvec((0., 0., angle)).as_matrix()
    prepared, native = volume_samples(F=R)
    out = audit.recover_level(prepared, native)
    np.testing.assert_allclose(out["c_eff"], 0., atol=1e-13)
    np.testing.assert_allclose(out["c_small"], np.cos(angle)-1., atol=1e-13)
    np.testing.assert_allclose(out["raw_fields"][:, 5], angle, atol=1e-13)
    assert np.max(abs(out["native_columns"]["thickness_green"])) < 1e-13


def test_volume_weighted_independent_strain_mean_and_geometry_moments():
    prepared, native = volume_samples(thickness_quadratic=.002, weights_unequal=True)
    out = audit.recover_level(prepared, native)
    np.testing.assert_allclose(out["native_columns"]["thickness_small"],
        .004*out["geometric_centroid_local"][:, 1], atol=1e-14)
    assert abs(out["geometric_centroid_local"][0, 1]) > 1e-4
    assert out["summary"]["native_transverse_heterogeneity_std_max"] > 0
    np.testing.assert_allclose(out["fit_volume"].sum(), .02, atol=1e-14)
    np.testing.assert_allclose(out["slab_volume_to_nominal_ratio"], 1., atol=1e-14)
    assert "not exact clipped" in audit.DEFINITIONS["slab_quadrature_qualification"]


def test_quadratic_omitted_from_fit_leaks_via_weighted_third_moment():
    coefficient = .002
    prepared, native = volume_samples(thickness_quadratic=coefficient, weights_unequal=True)
    out = audit.recover_level(prepared, native)
    mean = out["geometric_centroid_local"][:, 1]
    variance = out["geometric_covariance_local"][:, 1, 1]
    third = out["geometric_third_central_moment_local"][:, 1]
    predicted = coefficient*(2*mean+third/variance)
    np.testing.assert_allclose(out["c_small"], predicted, atol=1e-13)
    assert np.max(abs(out["c_small"]-2*coefficient*mean)) > 1e-6
    assert out["summary"]["c_small_minus_native_gradient_mean_max"] > 1e-6


def test_separate_quadratic_probe_exact_P2_removes_false_affine_slope():
    prepared, native = volume_samples(thickness_quadratic=.002, weights_unequal=True)
    historical = audit.recover_level(prepared, native)
    probe = audit.quadratic_transverse_fit_probe(prepared, native)
    assert probe["status"] == "COMPLETE_DIAGNOSTIC_PROBE"
    assert not probe["replaces_historical_recovery"]
    assert not probe["added_dynamic_constraints"]
    assert np.all(probe["fit_rank"] == 11)
    assert np.max(abs(probe["c_eff"])) < 1e-13
    assert np.max(abs(probe["c_small"])) < 1e-13
    assert np.max(abs(probe["fit_transverse_gradient_residual_RMS"])) < 1e-13
    assert np.max(abs(historical["c_small"])) > 1e-6
    np.testing.assert_allclose(historical["c_small"]-probe["c_small"],
        probe["omitted_quadratic_alias_c_small"], atol=1e-13)
    assert probe["nested_fit_alias_coefficient_identity_max_error"] < 1e-14


def test_quadratic_probe_preserves_affine_rotation_and_contraction():
    R = Rotation.from_rotvec((0., 0., .11)).as_matrix()
    prepared, native = volume_samples(F=R@np.diag((1.002, .98, 1.03)))
    primary = audit.recover_level(prepared, native)
    probe = audit.quadratic_transverse_fit_probe(prepared, native)
    np.testing.assert_allclose(probe["raw_fields"], primary["raw_fields"], atol=1e-13)
    assert not probe["smoothed"]


def test_cubic_interpolation_preserves_raw_teeth_and_additional_extrema():
    x = np.array((.1, .3, .5, .7, .9))
    y = np.array((0., 1., 0., 1., 0.))*1e-5
    out = audit.interpolation_diagnostics(x, y, 1.)
    assert out["raw_roughness"]["raw_extrema_count"] == 3
    assert out["interior_extrema_count"] >= 3
    assert out["raw_knot_reproduction_max_error"] < 1e-18
    assert out["maximum_neighbor_range_overshoot"] > 0
    assert out["cubic_minus_linear_max"] > 0
    assert out["smoothed"] is False
    assert out["interpolated_values_are_native_samples"] is False
    for row in out["stationary_extrema"]:
        assert 0 < row["x"] < 1.


def test_coordinate_scaled_roughness_unequal_spacing():
    for x in (np.linspace(.1, .9, 41), np.linspace(.1, .9, 81), np.linspace(0., 1., 81)**1.4):
        out = audit.profile_roughness(x, 2*x*x+3*x-1.)
        assert out["second_derivative_RMS"] == pytest.approx(4., abs=2e-10)
        assert out["total_variation"] == pytest.approx(abs(2*x[-1]**2+3*x[-1]-1.-(2*x[0]**2+3*x[0]-1.)))


def test_41_81_comparison_same_coordinates_no_strain_endpoint_constraint():
    a = audit.recover_level(*volume_samples(41, np.diag((1., .99, 1.02))), 41)
    b = audit.recover_level(*volume_samples(81, np.diag((1., .99, 1.02))), 81)
    out = audit.compare_levels(a, b)
    assert out["different_sample_centers"]
    assert not out["alignment_used"]
    assert out["native_strain_comparison_span"][0] > 0
    assert out["native_strain_comparison_span"][1] < 1.
    assert out["native_measures"]["thickness_small"]["absolute_max"] < 1e-14
    assert out["metrics"]["c_eff"]["linear"]["absolute_max"] > 0  # zero FACE interpolation differs from interior constant.


def test_complete_finite_U_required():
    mesh, xyz = tetrahedron()
    prepared = audit.prepare_mesh(mesh)
    with pytest.raises(ValueError, match="Complete"):
        audit.native_quadrature_fields(prepared, xyz[:-1])
    xyz[0, 1] = np.nan
    with pytest.raises(ValueError, match="Complete"):
        audit.native_quadrature_fields(prepared, xyz)


def test_no_scientific_execution_routes_or_history_reads(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("A saved-field diagnostic must not solve physics")
    monkeypatch.setattr(audit.fem2, "fem2_static_newton", forbidden)
    monkeypatch.setattr(audit.fem2, "run_static_case", forbidden, raising=False)
    monkeypatch.setattr(audit.fem1, "run_case", forbidden, raising=False)
    mesh, xyz = tetrahedron()
    out = audit.native_quadrature_fields(audit.prepare_mesh(mesh), xyz*0)
    assert out["minimum_det_F"] == pytest.approx(1.)
    prepared, native = volume_samples()
    recovered = audit.recover_level(prepared, native)
    assert recovered["summary"]["minimum_fit_rank"] == 8


def test_outside_span_and_invalid_window_rejected():
    with pytest.raises(ValueError, match="without extrapolation"):
        audit.interpolation_diagnostics(np.linspace(.1, .9, 4), np.ones(4), 1., np.array((0., 1.1)))
    with pytest.raises(ValueError, match="21/41/81"):
        audit.analyze_state({}, np.zeros((3, 3)), (41, 96))
