"""Synthetic contracts only: no native jobs, new meshes or dynamic solves."""
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from scripts.lib import nlsp_spatial_fem_protocol as protocol


def tetrahedron():
    corners = np.array(((0., 0., 0.), (1., 0., 0.), (0., .1, 0.), (0., 0., .2)))
    mids = np.array([(corners[i]+corners[j])/2 for i, j in protocol.fem1.TET10_EDGES])
    xyz = np.vstack((corners, mids))
    return SimpleNamespace(nodes={i+1: p for i, p in enumerate(xyz)},
        solid_elements={1: tuple(range(1, 11))}), xyz


def synthetic_samples(rotation):
    xyz = np.array([((k+.5+d)/5, eta, zeta) for k in range(5)
        for d in (-.4, -.15, .15, .4) for eta in (-.04, -.01, .03)
        for zeta in (-.08, .02, .08)])
    displacement = xyz@(rotation-np.eye(3)).T
    return protocol.static.fem2_recover_reference_samples(xyz*protocol.LOCAL_SIGNS,
        displacement*protocol.LOCAL_SIGNS, np.full(len(xyz), .02/len(xyz)), 1., .1, .2, 5)


def deck(tmp_path, nonlinear=False):
    mesh, _ = tetrahedron()
    directory = tmp_path/"mesh"; directory.mkdir(exist_ok=True)
    include = directory/"solid_mesh.inp"; include.write_text("** immutable synthetic mesh fixture\n")
    audit = {"status": "PASS", "fixed_left_ids": [1, 3, 4], "fixed_right_ids": [2]}
    path = tmp_path/("nonlinear.inp" if nonlinear else "linear.inp")
    science, gate = protocol.write_spatial_input(path, directory, mesh, audit,
        material={"E": 1., "rho": 1., "nu": .3, "kappa": 5/6},
        g_n=.0011379800853485065, g_k=.001422475106685633,
        omega1=.6054167303477958, static_settings=dict(protocol.static.STATIC_CONTROL_DEFAULTS),
        nonlinear=nonlinear)
    return path, science, gate


def test_local_global_dead_load_resultant_contract():
    contract = protocol.load_contract(2., 3.)
    np.testing.assert_array_equal(contract["global_acceleration"], [0., -2., -3.])
    np.testing.assert_allclose(contract["total_force_global"], [0., -.04, -.06])
    np.testing.assert_allclose(contract["line_load_local"], [0., .04, .06])
    assert contract["distributed_torque"] == 0 and not contract["follower_load"]
    assert np.linalg.norm(contract["global_direction"]) == pytest.approx(1.)


@pytest.mark.parametrize("g_n,g_k", [(np.nan, 1.), (1., np.inf), (0., 0.)])
def test_invalid_load_stops_before_generation(g_n, g_k):
    with pytest.raises(ValueError):
        protocol.load_contract(g_n, g_k)


@pytest.mark.parametrize("nonlinear", [False, True])
def test_safe_generator_reuses_unchanged_material_clamps_and_release(tmp_path, nonlinear):
    path, science, gate = deck(tmp_path, nonlinear)
    text = path.read_text(); static, dynamic = text.split("*END STEP", 1)
    assert gate["status"] == "PASS" and gate["numeric_width_max"] <= 20
    assert "ELKE" not in static if not nonlinear else "ELSE,ELKE" in static
    assert "ELSE,ELKE" in dynamic and "*DLOAD, OP=NEW\nSOLID,GRAV,0." in dynamic
    assert "*DYNAMIC, ALPHA=0" in dynamic and "AMPLITUDE=STEP" in dynamic
    assert "LEFT_FIXED,1,3,0" in static and "RIGHT_FIXED,1,3,0" in static
    assert "*MPC" not in text and "*SPRING" not in text and "*DAMPING" not in text
    assert "FREQUENCY=2" in dynamic and "*EL FILE, GLOBAL=YES, FREQUENCY=2\nS,E,ENER" in dynamic
    assert science["horizon_T1"] == .25 and science["dynamic"] == protocol.DEFAULT_DYNAMIC
    assert not gate["historical_preload_equality_required"]


def test_linear_nonlinear_physical_decks_differ_only_routing_and_static_energy(tmp_path):
    first, _, _ = deck(tmp_path, False); second, _, _ = deck(tmp_path, True)
    def normalize(text):
        return text.replace(", NLGEOM=NO", "").replace(", NLGEOM", "").replace("ELSE,ELKE", "ELSE")
    assert normalize(first.read_text()) == normalize(second.read_text())


@pytest.mark.parametrize("replacement", ["SOLID,GRAV,1.,0.,-1.,0.", "SOLID,GRAV,0.,0.,0.,0."])
def test_nonzero_or_invalid_release_rejected(tmp_path, replacement):
    path, science, _ = deck(tmp_path)
    text = path.read_text().replace("SOLID,GRAV,0.,0.,-1.,0.", replacement)
    with pytest.raises(ValueError):
        protocol.input_contract(text, science, nonlinear=False)


def test_overwidth_numeric_field_not_whole_line_rejected(tmp_path):
    path, science, _ = deck(tmp_path)
    text = path.read_text().replace("*DENSITY\n1.", "*DENSITY\n1.000000000000000000000")
    with pytest.raises(ValueError, match="20-character"):
        protocol.input_contract(text, science, nonlinear=False)
    assert protocol._numeric_cards("*CONTROLS, PARAMETERS=FIELD\n1.e-8,1.e-8,,,1.e-8,,1.e-8,1.e-8\n") < 20


@pytest.mark.parametrize("old,new", [("*DENSITY\n1.", "*DENSITY\n2."),
    ("RIGHT_FIXED,1,3,0", "RIGHT_FIXED,1,2,0"),
    ("*DYNAMIC, ALPHA=0", "*BOUNDARY\nALL_NODES,3,3,0\n*DYNAMIC, ALPHA=0")])
def test_frozen_material_and_no_artificial_constraints(tmp_path, old, new):
    path, science, _ = deck(tmp_path)
    with pytest.raises(ValueError):
        protocol.input_contract(path.read_text().replace(old, new), science, nonlinear=False)


def test_refined_time_and_bounded_resource_scope():
    assert protocol.DEFAULT_DYNAMIC["initial_T1_fraction"] == 1/8000
    assert protocol.DEFAULT_DYNAMIC["maximum_T1_fraction"] == 1/4000
    assert protocol.DEFAULT_DYNAMIC["output_frequency"] == 2
    assert protocol.RESOURCE_POLICY["maximum_production_jobs"] == 4
    assert protocol.RESOURCE_POLICY["job_timeouts_seconds"] == {"medium": 1800, "fine": 4200}
    assert not protocol.RESOURCE_POLICY["automatic_retry"]


def test_quadrature_vector_bodyload_exact_consistent_resultant():
    mesh, xyz = tetrahedron()
    ids, nodes, load, volume = protocol.consistent_bodyloads(mesh, 2., [0., -2., -3.])
    np.testing.assert_array_equal(nodes, xyz)
    np.testing.assert_array_equal(ids, np.arange(1, 11))
    np.testing.assert_allclose(load.sum(axis=0), 2.*volume*np.array([0., -2., -3.]), atol=1e-17)
    weights = np.r_[np.full(4, -volume/20), np.full(6, volume/5)]
    np.testing.assert_allclose(load, 2.*weights[:, None]*np.array([0., -2., -3.]), atol=1e-17)
    old_ids, old_nodes, old_load, old_volume = protocol.static.consistent_gravity_loads(mesh, 2., 2.)
    np.testing.assert_allclose(load[:, 1], old_load[:, 1], atol=1e-17)
    assert volume == pytest.approx(old_volume)


def test_independent_support_recovery_detects_force_and_moment_imbalance():
    mesh, xyz = tetrahedron()
    audit = {"fixed_left_ids": list(range(1, 6)), "fixed_right_ids": list(range(6, 11))}
    # All fixture DOFs fixed: RF=body+support=0; support is independently -body.
    RF = {"LEFT_FIXED": np.zeros((5, 3)), "RIGHT_FIXED": np.zeros((5, 3))}
    args = dict(rho=1., acceleration_global=[0., -2., -3.], nonlinear=False)
    result = protocol.support_equilibrium(mesh, audit, np.zeros_like(xyz), RF, **args)
    assert result["status"] == "PASS"
    assert result["force_imbalance_relative"] < 1e-14
    assert result["moment_imbalance_relative"] < 1e-14
    RF["RIGHT_FIXED"][0, 0] = 1e-4
    assert protocol.support_equilibrium(mesh, audit, np.zeros_like(xyz), RF, **args)["status"] == "FAIL"


def test_equilibrium_gate_cannot_be_loosened():
    mesh, xyz = tetrahedron()
    with pytest.raises(ValueError, match="unchanged"):
        protocol.support_equilibrium(mesh, {}, np.zeros_like(xyz), {}, rho=1.,
            acceleration_global=[0., -1., 0.], nonlinear=False, equilibrium_relative=1e-3)


@pytest.mark.parametrize("a", [(0., 0., .05), (.04, -.03, .06)])
def test_spatial_canonical_rotation_and_planar_restriction(a):
    R = Rotation.from_rotvec(a).as_matrix()
    original = synthetic_samples(R)
    fixed_original = np.asarray(original["fields"]).copy()
    recovered = protocol.canonical_section_profile(original)
    np.testing.assert_array_equal(original["fields"], fixed_original)
    np.testing.assert_allclose(recovered["raw_rotation_matrices"], np.broadcast_to(R, (5, 3, 3)), atol=2e-15)
    np.testing.assert_allclose(recovered["fields"][1:-1, 3:6], np.broadcast_to(np.array(a)*[1., -1., 1.], (5, 3)), atol=2e-15)
    np.testing.assert_allclose(recovered["fields"][:, 6], 0., atol=2e-15)
    np.testing.assert_array_equal(recovered["fields"][[0, -1]], np.zeros((2, 7)))
    if a[0] == 0:
        np.testing.assert_allclose(recovered["fields"][:, 5], original["fields"][:, 5], atol=1e-15)
    else:
        assert np.max(abs(recovered["fields"][:, 5]-recovered["historical_projected_theta"])) > 1e-4


def test_rotation_recovery_preserves_uniform_stretch_and_axes():
    R = Rotation.from_rotvec([.02, -.01, .03]).as_matrix()
    F = R@np.diag([1., .98, 1.03])
    profile = synthetic_samples(F)
    recovered = protocol.canonical_section_profile(profile)
    np.testing.assert_allclose(recovered["fields"][1:-1, 6], -.02, atol=2e-15)
    np.testing.assert_allclose([row["effective_width_strain"] for row in recovered["section_rows"]], .03, atol=2e-15)


def test_native_completion_failure_retains_original_outputs(tmp_path):
    raw = tmp_path/"motion.stdout.txt"; raw.write_text("*ERROR actual native failure\n")
    with pytest.raises(ValueError, match="completion"):
        protocol.recover_saved_outputs(tmp_path, None, {}, {}, nonlinear=False)
    assert raw.read_text() == "*ERROR actual native failure\n"
    assert not (tmp_path/"section_history.npz").exists()


def test_no_native_execution_or_ODE_entrypoint():
    import inspect
    source = inspect.getsource(protocol)
    assert "subprocess" not in source and "Popen" not in source
    assert "solve_ivp" not in source and "Radau(" not in source
    assert "run_job(" not in source and "static_newton(" not in source
    assert "historical_preload_equality_required\": False" in source


def mock_finished_output(tmp_path, monkeypatch, *, prefix=False):
    """Feed synthetic native-reader blocks; no input or solver is executed."""
    mesh, xyz = tetrahedron(); ids = np.arange(1, 11)
    audit = {"fixed_left_ids": list(range(1, 6)), "fixed_right_ids": list(range(6, 11))}
    science = {"material": {"rho": 1.}, "omega1": .6054167303477958,
        "horizon_T1": .25, "dynamic": protocol.DEFAULT_DYNAMIC,
        "load": protocol.load_contract(.001, .00125)}
    end = protocol.pilot.dynamic_settings(science)["duration"]
    (tmp_path/"motion.stdout.txt").write_text("JOB FINISHED\n")
    records = [{"step": 1, "increment": 1, "step_time": 1., "total_time": 1.},
        {"step": 2, "increment": 1, "step_time": end-.01 if prefix else end,
         "total_time": 1.+end, "time_rounding_bounds": {"step_time": 5e-7}}]
    monkeypatch.setattr(protocol.io, "read_transient_sta", lambda path:
        {"accepted_increments": records, "reported_cutbacks": 0})
    dat = []
    for step in (1, 2):
        dat.append({"step": step, "increment": 1, "set": "ALL_NODES", "name": "DISP", "values": np.zeros_like(xyz)})
    for name in ("LEFT_FIXED", "RIGHT_FIXED"):
        dat.append({"step": 1, "increment": 1, "set": name, "name": "FORC", "values": np.zeros((5, 3))})
    monkeypatch.setattr(protocol.io, "iter_transient_dat", lambda *args, **kwargs: iter(dat))
    frames = [{"step": step, "increment": 1, "fixed_displacement_max": 0.,
        "fixed_velocity_max": 0., "dynamic_time": end if step == 2 else None,
        "fields": {"DISP": np.zeros_like(xyz), "VELO": np.zeros_like(xyz),
                   "STRESS": np.zeros((10, 6)), "TOSTRAIN": np.zeros((10, 6))}}
        for step in (1, 2)]
    monkeypatch.setattr(protocol.io, "iter_transient_frd", lambda *args, **kwargs: iter(frames))
    def profile(amplitude):
        fields = np.zeros((5, 7)); fields[1:-1, 1:3] = amplitude
        return {"x": np.linspace(0., 1., 5), "fields": fields,
            "raw_rotation_matrices": np.broadcast_to(np.eye(3), (3, 3, 3)).copy(),
            "historical_projected_theta": np.zeros(5)}
    recovery_order = iter((.5, 1.))
    monkeypatch.setattr(protocol, "recover_spatial_sections", lambda *args, **kwargs: profile(next(recovery_order)))
    monkeypatch.setattr(protocol.static, "fem2_recover_reference_samples", lambda *args, **kwargs: profile(0.))
    monkeypatch.setattr(protocol.io, "parse_transient_dat_energies", lambda *args, **kwargs:
        {"records": [{"step": 1, "total_time": 1., "internal_energy": 1.},
            {"step": 2, "mechanical_energy": 2.}]})
    monkeypatch.setattr(protocol.io, "parse_transient_stdout_energies", lambda *args, **kwargs:
        {"records": [{"step": 2, "external_work": 0., "damping_work": 0., "initial_step_energy": 2.}]})
    return mesh, audit, science


def test_finished_output_recovery_keeps_static_anchor_actual_frames_and_raw_energy(tmp_path, monkeypatch):
    mesh, audit, science = mock_finished_output(tmp_path, monkeypatch)
    result = protocol.recover_saved_outputs(tmp_path, mesh, audit, science,
                                          nonlinear=False, section_count=5)
    assert result["status"] == "PASS" and result["dynamic_output_frames"] == 1
    assert result["energy_status"] == "PARTIAL"
    assert result["native_energy_static_to_dynamic_reference_jump_relative"] == 1.
    assert result["preload_equilibrium"]["status"] == "PASS"
    assert not result["historical_preload_equality_required"]
    with np.load(tmp_path/"initial_sections.npz") as saved:
        assert saved["not_a_native_dynamic_zero_frame"]
    with np.load(tmp_path/"section_history.npz") as saved:
        assert np.all(saved["time"] > 0.)
        assert saved["time"][-1] == protocol.pilot.dynamic_settings(science)["duration"]
        np.testing.assert_array_equal(saved["raw_rotation_matrices"], np.broadcast_to(np.eye(3), (1, 3, 3, 3)))


def test_incomplete_accepted_prefix_is_not_padded_or_retried(tmp_path, monkeypatch):
    mesh, audit, science = mock_finished_output(tmp_path, monkeypatch, prefix=True)
    with pytest.raises(ValueError, match="prefix"):
        protocol.recover_saved_outputs(tmp_path, mesh, audit, science,
                                      nonlinear=False, section_count=5)
    assert (tmp_path/"motion.stdout.txt").read_text() == "JOB FINISHED\n"
    assert not (tmp_path/"section_history.npz").exists()
