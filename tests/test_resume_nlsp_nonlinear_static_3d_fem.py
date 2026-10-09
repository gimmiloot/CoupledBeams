"""FEM-2R contract tests: synthetic I/O and saved evidence, no scientific jobs."""
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from scripts.analysis import verify_nlsp_nonlinear_static_3d_fem as old

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT / "results/nlsp_nonlinear_static_3d_fem/6714bd9f2778e6d7"


def _mesh():
    corners = np.array(((0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.)))
    edges = ((0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3))
    xyz = np.vstack((corners, [(corners[a] + corners[b]) / 2 for a, b in edges]))
    return SimpleNamespace(nodes={i + 1: tuple(v) for i, v in enumerate(xyz)},
                           solid_elements={1: tuple(range(1, 11))})


def _audit():
    return {"status": "PASS", "fixed_left_ids": [1, 3, 4], "fixed_right_ids": [2]}


def _frd_block(counter, name, values, time=1., increment=10):
    header = list(" " * 75)
    header[0:7] = "  100CL"
    header[7:12] = f"{100 + increment:5d}"
    header[12:24] = f"{time:12.9f}"
    header[24:36] = f"{len(values):12d}"
    result = f"    1PSTEP{counter:26d}{increment:12d}{1:12d}\n" + "".join(header) + "\n"
    result += f" -4  {name:<12} {values.shape[1]:d}    1\n"
    for node, value in enumerate(values, 1):
        result += f" -1{node:10d}" + "".join(f"{float(x):12.5E}" for x in value) + "\n"
    return result + " -3\n"


def _static_outputs(tmp_path, nonlinear=False, time=1.):
    """Parser fixture with analytically known consistent load and moments."""
    mesh = _mesh(); audit = _audit(); rho = 1.; g = .0014224751066856333
    ids, xyz, body, volume = old.consistent_gravity_loads(mesh, rho, g)
    U = np.zeros((len(ids), 3))
    reaction = np.zeros_like(body)
    magnitude = rho * g * volume
    # tetra centroid=(.25,.25,.25); support positions give exact force/moments.
    reaction[0, 1] = .5 * magnitude
    reaction[1, 1] = .25 * magnitude
    reaction[3, 1] = .25 * magnitude
    RF = body + reaction
    stem = tmp_path / "static"
    stem.with_suffix(".frd").write_text("".join(
        _frd_block(i, name, values, time) for i, (name, values) in enumerate((
            ("DISP", U), ("FORC", RF), ("STRESS", np.zeros((10, 6))),
            ("TOSTRAIN", np.zeros((10, 6)))), 1)), encoding="utf8")
    dat = f" displacements (vx,vy,vz) for set ALL_NODES and time {time:.12E}\n"
    dat += "".join(f"{n} " + " ".join(f"{v:.12E}" for v in U[n-1]) + "\n" for n in ids)
    for label, nodes in (("LEFT_FIXED", audit["fixed_left_ids"]),
                         ("RIGHT_FIXED", audit["fixed_right_ids"])):
        dat += f" forces (fx,fy,fz) for set {label} and time {time:.12E}\n"
        dat += "".join(f"{n} " + " ".join(f"{v:.12E}" for v in RF[n-1]) + "\n" for n in nodes)
    stem.with_suffix(".dat").write_text(dat, encoding="utf8")
    stem.with_suffix(".sta").write_text(f" STEP INC ATT ITRS TOT TIME STEP TIME INC TIME\n1 10 1 3 {time} {time} .1\n", encoding="utf8")
    stem.with_suffix(".stdout.txt").write_text(
        ("Nonlinear geometric effects are taken into account\n" if nonlinear else "") + "JOB FINISHED\n", encoding="utf8")
    stem.with_suffix(".stderr.txt").write_text("", encoding="utf8")
    return stem, mesh, audit, rho, g, body, reaction


def test_historical_failed_input_and_preview_are_distinct():
    bad = (PARENT / "cases/medium/linear/static.inp").read_text(encoding="utf8")
    fixed = (PARENT / "remediation_preview/medium_linear_corrected.inp").read_text(encoding="utf8")
    bad_fields = bad.split("*STATIC\n", 1)[1].splitlines()[0].split(",")
    fixed_fields = fixed.split("*STATIC\n", 1)[1].splitlines()[0].split(",")
    assert max(len(x.strip()) for x in bad_fields) > 20
    assert max(len(x.strip()) for x in fixed_fields) <= 20
    assert old.read_json(PARENT / "summary.json")["cases"]["medium"]["linear"]["job"]["returncode"] == 201


@pytest.mark.parametrize("value", [0., -1., 1., 1e-6, -1e-8, 1e-200, -1e-200, 1e200, -1e200, .0014224751066856333])
def test_existing_bounded_numeric_serializer_preserves_native_width_and_value(value):
    token = old.fem1.single.ccx_float(value)
    assert len(token) <= 20
    assert np.isfinite(float(token))
    assert float(token) == pytest.approx(value, rel=1e-12, abs=0.)


def test_all_relevant_static_cards_have_individually_bounded_fields(tmp_path):
    path = tmp_path / "linear.inp"
    old.write_static_input(path, tmp_path / "mesh.inp", _mesh(), _audit(),
                           {"E": 1., "rho": 1., "nu": .3}, .0014224751066856333, False)
    lines = path.read_text().splitlines()
    for card in ("*STATIC", "*CONTROLS, PARAMETERS=FIELD", "*DLOAD", "*ELASTIC", "*DENSITY"):
        payload = lines[lines.index(card) + 1]
        for token in payload.split(","):
            token = token.strip()
            if token and token not in ("SOLID", "GRAV"):
                assert len(token) <= 20, (card, token)
                assert np.isfinite(float(token)), (card, token)
    fields = lines[lines.index("*STATIC") + 1].split(",")
    np.testing.assert_allclose([float(s) for s in fields], [.1, 1., 1e-6, .1], rtol=1e-12)


def test_linear_nonlinear_decks_share_all_physical_inputs(tmp_path):
    arguments = (tmp_path / "mesh.inp", _mesh(), _audit(), {"E": 1., "rho": 1., "nu": .3}, .0014224751066856333)
    old.write_static_input(tmp_path / "lin.inp", *arguments, False)
    old.write_static_input(tmp_path / "nl.inp", *arguments, True)
    linear = (tmp_path / "lin.inp").read_text()
    nonlinear = (tmp_path / "nl.inp").read_text()
    assert nonlinear.replace(", NLGEOM", "") == linear
    for forbidden in ("*MPC", "*RIGID", "*SPRING", "*CONTACT", "*DYNAMIC", "*FREQUENCY", "*CLOAD"):
        assert forbidden not in nonlinear
    assert "LEFT_FIXED,1,3,0" in linear and "RIGHT_FIXED,1,3,0" in linear


@pytest.mark.parametrize("nonlinear", [False, True])
def test_static_output_reactions_remove_consistent_load_independently(tmp_path, nonlinear):
    stem, mesh, audit, rho, g, body, expected = _static_outputs(tmp_path, nonlinear)
    diagnostic, arrays = old.parse_static_outputs(stem, mesh, audit, rho, g, nonlinear)
    assert diagnostic["status"] == "PASS", diagnostic["failures"]
    assert diagnostic["final_dat"]["time"] == diagnostic["final_frd"]["time"] == 1.
    assert diagnostic["actual_NLGEOM"] is nonlinear
    assert diagnostic["support_force_balance_relative"] < 1e-12
    assert diagnostic["support_moment_balance_relative"] < 1e-12
    np.testing.assert_allclose(diagnostic["left"]["reaction_force_global"], expected[[0, 2, 3]].sum(axis=0), atol=1e-16)
    np.testing.assert_allclose(diagnostic["right"]["reaction_force_global"], expected[[1]].sum(axis=0), atol=1e-16)
    assert np.linalg.norm(diagnostic["left"]["raw_RF_global"]) != pytest.approx(np.linalg.norm(diagnostic["left"]["reaction_force_global"]), rel=1e-8, abs=0.)
    assert diagnostic["maximum_fixed_displacement"] == 0.
    assert arrays["U"].shape == (10, 3)


def test_incomplete_final_load_is_not_accepted(tmp_path):
    stem, mesh, audit, rho, g, _, _ = _static_outputs(tmp_path, True, .9)
    diagnostic, _ = old.parse_static_outputs(stem, mesh, audit, rho, g, True)
    assert diagnostic["status"] == "FAIL"
    assert "FINAL_LOAD_FACTOR_NOT_ONE" in diagnostic["failures"]
    assert "NO_COMPLETE_NL_LOAD_INCREMENT_HISTORY" in diagnostic["failures"]


def test_actual_nlgeom_routing_required(tmp_path):
    stem, mesh, audit, rho, g, _, _ = _static_outputs(tmp_path, False)
    diagnostic, _ = old.parse_static_outputs(stem, mesh, audit, rho, g, True)
    assert "NLGEOM_ACTUAL_ROUTING_MISMATCH" in diagnostic["failures"]


def test_job_finished_and_warning_gates_are_not_exit_code_only(tmp_path):
    stem, mesh, audit, rho, g, _, _ = _static_outputs(tmp_path)
    stem.with_suffix(".stdout.txt").write_text("*WARNING fixture problem\n")
    diagnostic, _ = old.parse_static_outputs(stem, mesh, audit, rho, g, False)
    assert diagnostic["status"] == "FAIL"
    assert "SOLVER_WARNING_OR_ERROR" in diagnostic["failures"]
    assert "NO_REAL_JOB_FINISHED_MARKER" in diagnostic["failures"]


def test_missing_nodal_displacement_cannot_be_filled_with_previous_value(tmp_path):
    stem, mesh, audit, _, _, _, _ = _static_outputs(tmp_path)
    path = stem.with_suffix(".dat")
    path.write_text(path.read_text().replace("10 0.000000000000E+00 0.000000000000E+00 0.000000000000E+00\n", "", 1))
    with pytest.raises(ValueError, match="Incomplete"):
        old.read_static_dat(path, list(mesh.nodes), audit["fixed_left_ids"], audit["fixed_right_ids"])


def test_correction_arithmetic_keeps_sign_and_common_scale():
    x = np.linspace(0, 1, 801)
    fields = np.zeros((801, 7)); fields[:, 1] = .08 * x*x*(1-x)**2
    nl = fields.copy(); nl[:, 1] *= .998
    f3d = fields*.97; n3d = f3d.copy(); n3d[:, 1] *= .997
    profiles = lambda value: {"x": x, "fields": value}
    result = old.compare_static_profiles(profiles(fields), profiles(nl), profiles(f3d), profiles(n3d), correction_scales={n: 1e-5 for n in old.FEM2_STATIC_FIELD_ORDER})
    assert result["metrics"]["w"]["correction_signed_midpoint_1D"] < 0.
    assert result["metrics"]["w"]["correction_signed_midpoint_3D"] < 0.
    assert not result["amplitude_or_space_alignment"]
    assert not result["frequency_mesh_criterion_applied"]


def test_nonzero_signal_does_not_hide_larger_numerical_uncertainty():
    result = old.fem2_signal_resolution(8e-6, 1e-5, 1e-7, 2e-8)
    assert result["status"] == "UNRESOLVED_SIGNAL"
    assert not result["continuum_error_bound_claimed"]


def test_saved_p64_profiles_are_loaded_without_new_static_preparation(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("No repeated 1D static preparation")
    monkeypatch.setattr(old, "build_static_preflight", forbidden)
    monkeypatch.setattr(old, "fem2_static_newton", forbidden)
    before = old.sha(PARENT / "one_d_p64.npz")
    linear = old.load_one_d_static_profile(PARENT, 64, False)
    nonlinear = old.load_one_d_static_profile(PARENT, 64, True)
    assert before == old.sha(PARENT / "one_d_p64.npz")
    assert linear["fields"].shape == nonlinear["fields"].shape == (1001, 7)
    assert linear["fields"][500, 1] == pytest.approx(.005, rel=1e-12)
    assert nonlinear["fields"][500, 1] == pytest.approx(.00499112992069312, rel=1e-12)
    assert nonlinear["fields"][500, 1] - linear["fields"][500, 1] == pytest.approx(-8.8700793068e-6, rel=1e-10)
    assert np.max(abs(linear["fields"][:, [2, 3, 4]])) == 0.


# Continuation contracts; all references are historical read-only artifacts.
from copy import deepcopy

RESUME_CONFIG = ROOT / "data/input/nlsp_nonlinear_static_3d_fem_resume.json"


def _resume():
    from scripts.analysis import resume_nlsp_nonlinear_static_3d_fem as module
    return module


@pytest.fixture(scope="module")
def resume_config():
    return old.read_json(RESUME_CONFIG)


def test_explicit_continuation_authorization_is_separate_from_historical_failure(resume_config):
    authorization = resume_config["authorization"]
    assert authorization["id"] == "explicit_user_FEM2R_2026_10_09"
    assert authorization["maximum_new_ccx_jobs"] == 6
    assert authorization["case_order"] == [
        "medium_linear", "medium_nonlinear", "fine_linear", "fine_nonlinear",
        "refined_linear", "refined_nonlinear"]
    assert authorization["automatic_retry_after_new_solver_failure"] is False
    assert resume_config["parent_failed_bundle"] == str(PARENT.relative_to(ROOT)).replace("\\", "/")
    assert old.sha(PARENT / "manifest.json") == resume_config["parent_manifest_sha256"]


def test_resume_geometry_material_and_load_are_exactly_frozen(resume_config):
    assert resume_config["geometry"] == {"L": 1, "b": .2, "h": .1}
    assert resume_config["material"] == {"E": 1, "rho": 1, "nu": .3, "kappa": 5/6}
    load = resume_config["frozen_load"]
    assert load["g"] == .0014224751066856333
    assert load["q"] == load["F_total"] == 2.844950213371267e-5
    assert load["q"] == pytest.approx(resume_config["material"]["rho"]*.2*.1*load["g"], rel=3*np.finfo(float).eps, abs=0.)


def test_continuation_cannot_authorize_other_physics_or_hidden_jobs(resume_config):
    assert resume_config["semantics"] == {
        "new_static_jobs_only": True, "new_meshes": False, "new_modal_jobs": False,
        "new_dynamics": False, "new_one_d_solves": False, "model_fitting": False}
    assert resume_config["one_d_policy"] == "reuse_parent_p48_p64_without_equilibrium_solves"
    assert resume_config["threads"] == 1
    assert resume_config["job_timeout_seconds"] <= 1200
    assert resume_config["job_memory_limit_bytes"] <= 4*1024**3
    assert resume_config["numerical_budget_seconds"] <= 3600


def test_all_parent_artifact_hashes_are_unchanged(resume_config):
    manifest = old.read_json(PARENT / "manifest.json")
    assert old.sha(PARENT / "manifest.json") == resume_config["parent_manifest_sha256"]
    for name, digest in manifest["artifact_hashes"].items():
        assert old.sha(PARENT / name) == digest, name
    assert old.sha(old.__file__) == resume_config["corrected_generator_sha256"]


def _reference_samples(count=41):
    local = np.array([( (k+.5+d)/count, y, z)
                      for k in range(count) for d in (-.4, -.15, .15, .4)
                      for y in (-.04, .04) for z in (-.08, .08)])
    return local, local*np.array((1., -1., -1.)), np.full(len(local), .02/len(local))


def test_reference_sections_use_original_coordinates_and_transverse_sign():
    local, xyz, weights = _reference_samples()
    global_U = np.zeros_like(xyz)
    global_U[:, 1] = -.005 * local[:, 0]**2
    result = old.fem2_recover_reference_samples(xyz, global_U, weights, 1., .1, .2,
                                               enforce_clamped_faces=False)
    np.testing.assert_allclose(result["fields"][:, 1], .005*result["x"]**2, atol=1e-12)
    assert result["reference_geometry_used_for_sections"]
    assert result["additional_derivative_constraints"] is False
    assert result["recovery_policy"] == old.FEM2_RECOVERY_POLICY


def test_effective_contraction_is_a_section_diagnostic_not_a_cartesian_dof():
    local, xyz, weights = _reference_samples()
    global_U = np.zeros_like(xyz)
    global_U[:, 1] = .002 * local[:, 1]
    result = old.fem2_recover_reference_samples(xyz, global_U, weights, 1., .1, .2,
                                               enforce_clamped_faces=False)
    np.testing.assert_allclose(result["fields"][:, 6], -.002, atol=1e-12)
    assert result["contraction_status"].endswith("NOT_MH_DOF")


def test_same_polar_rotation_policy_is_available_for_each_static_case():
    local, xyz, weights = _reference_samples()
    angle = .012
    rotation = np.array(((np.cos(angle), -np.sin(angle), 0.),
                         (np.sin(angle), np.cos(angle), 0.), (0., 0., 1.)))
    global_U = (local@(rotation-np.eye(3)).T)*np.array((1., -1., -1.))
    result = old.fem2_recover_reference_samples(xyz, global_U, weights, 1., .1, .2,
                                               enforce_clamped_faces=False)
    np.testing.assert_allclose(result["fields"][:, 5], angle, atol=1e-12)
    assert "both_linear_and_nonlinear" in result["theta_policy"]


def test_historical_strict_float64_qualification_is_not_reclassified():
    summary = old.read_json(PARENT / "summary.json")
    cases = summary["preflight"]["cases"]
    assert [r["p"] for r in cases] == [48, 64]
    statuses = []
    for row in cases:
        for kind in ("linear", "nonlinear"):
            record = row[kind]
            assert record["strong_action_strict_threshold"] == 2e-12
            statuses.append(record["strong_action_strict_status"])
        assert row["nonlinear_status"] == "PASS"
        assert row["nonlinear"]["relative_residual"] <= 1e-10
    assert "PARTIAL" in statuses


def test_resume_config_validates_without_weakening_frozen_contract(resume_config):
    assert _resume().validate_resume_config(resume_config) == resume_config


@pytest.mark.parametrize("path,value", [
    (("geometry", "h"), .12),
    (("material", "nu"), .31),
    (("frozen_load", "g"), .0015),
    (("frozen_load", "q"), 3e-5),
    (("authorization", "maximum_new_ccx_jobs"), 7),
    (("authorization", "automatic_retry_after_new_solver_failure"), True),
    (("threads",), 2),
    (("job_timeout_seconds",), 1201),
    (("job_memory_limit_bytes",), 4*1024**3+1),
    (("numerical_budget_seconds",), 3601),
    (("semantics", "new_meshes"), True),
    (("semantics", "new_one_d_solves"), True),
    (("semantics", "new_dynamics"), True),
])
def test_resume_rejects_unrequested_physics_budget_or_retry_changes(resume_config, path, value):
    altered = deepcopy(resume_config)
    target = altered
    for part in path[:-1]:
        target = target[part]
    target[path[-1]] = value
    with pytest.raises(ValueError):
        _resume().validate_resume_config(altered)


def test_resume_rejects_changed_case_order(resume_config):
    altered = deepcopy(resume_config)
    altered["authorization"]["case_order"][0:2] = ["medium_nonlinear", "medium_linear"]
    with pytest.raises(ValueError):
        _resume().validate_resume_config(altered)


def test_corrupted_parent_manifest_is_not_reconstructed(resume_config, monkeypatch):
    module = _resume()
    def forbidden(*args, **kwargs):
        raise AssertionError("Damaged historical artifacts must not be recreated")
    monkeypatch.setattr(old, "build_static_preflight", forbidden)
    monkeypatch.setattr(old, "fem2_static_newton", forbidden)
    monkeypatch.setattr(old.fem1, "run_job", forbidden)
    altered = deepcopy(resume_config)
    altered["parent_manifest_sha256"] = "0"*64
    with pytest.raises(ValueError):
        module.load_resume_parent(altered)


@pytest.fixture(scope="module")
def parent_context(resume_config):
    return _resume().load_resume_parent(resume_config)


def _fresh_summary(config, context):
    _, historical, science = context
    module = _resume()
    return {
        "preflight": historical["preflight"], "science_config": science,
        "authorization": config["authorization"], "cases": {}, "attempt_ledger": [],
        "job_calls": {"ccx": 0, "gmsh": 0, "modal": 0, "nonlinear_ODE": 0, "one_d_static": 0},
        "runtime": {"numerical_seconds": 0., "budget_seconds": 3600},
        "statuses": {"NLSP_FEM2_"+n: "NOT_RUN" for n in old.FEM2_STATUS_NAMES},
        "resume_statuses": {"NLSP_FEM2R_"+n: "NOT_RUN" for n in module.STATUS_NAMES}}


def test_parent_loader_reuses_complete_sources_without_scientific_calls(resume_config, parent_context, monkeypatch):
    module = _resume()
    def forbidden(*args, **kwargs):
        raise AssertionError("Source loading must not launch any scientific calculation")
    monkeypatch.setattr(old, "build_static_preflight", forbidden)
    monkeypatch.setattr(old, "fem2_static_newton", forbidden)
    monkeypatch.setattr(old.fem1, "run_job", forbidden)
    monkeypatch.setattr(old.fem1.single, "generate_mesh_with_gmsh_cli", forbidden)
    parent, historical, science = module.load_resume_parent(resume_config)
    assert parent == PARENT
    assert historical["cases"]["medium"]["linear"]["status"] == "FAIL"
    assert science["mesh_levels"] == ["medium", "fine", "refined"]
    assert historical["preflight"]["cases"][1]["p"] == 64


def test_input_gate_recreates_corrected_cards_without_solver_calls(mesh_tmp, resume_config, parent_context, monkeypatch):
    tmp_path = mesh_tmp
    module = _resume()
    def forbidden(*args, **kwargs):
        raise AssertionError("Input validation must not solve or mesh")
    monkeypatch.setattr(old.fem1, "run_job", forbidden)
    monkeypatch.setattr(old.fem1.single, "generate_mesh_with_gmsh_cli", forbidden)
    monkeypatch.setattr(old, "build_static_preflight", forbidden)
    result = module.input_gate(resume_config, tmp_path, parent_context[1])
    assert result["status"] == "PASS"
    assert result["new_solver_calls"] == 0
    assert result["preview_agrees_except_include"]
    assert result["linear_NL_only_NLGEOM_difference"]
    assert [r["case"] for r in result["rows"]] == list(module.ORDER)
    assert max(result["old_failed_STATIC_widths"]) == 22
    for row in result["rows"]:
        assert row["maximum_numeric_token_width"] <= 20
        assert abs(row["relative_g_rounding"]) <= 1e-12
        assert set(row["cards"]) == {"STATIC", "CONTROLS", "ELASTIC", "DENSITY", "DLOAD"}


@pytest.mark.parametrize("bad", ["9.9999999999999995E-07", "NaN", "INF", "-INF"])
def test_numeric_card_gate_rejects_invalid_native_field(tmp_path, bad):
    module = _resume()
    p = tmp_path / "input.inp"
    old.write_static_input(p, tmp_path/"mesh.inp", _mesh(), _audit(),
                          {"E": 1., "rho": 1., "nu": .3}, .0014224751066856333, False)
    text = p.read_text()
    data = text.split("*STATIC\n", 1)[1].splitlines()[0]
    with pytest.raises(ValueError):
        module.numeric_cards(text.replace(data, f".1,1.,{bad},.1"))


def test_new_medium_linear_failure_stops_every_later_job_and_preserves_ledger(mesh_tmp, resume_config, parent_context, monkeypatch):
    tmp_path = mesh_tmp
    module = _resume()
    summary = _fresh_summary(resume_config, parent_context)
    calls = []
    def failed_job(command, cwd, timeout, memory, prefix, env):
        calls.append((command, timeout, memory, env))
        Path(cwd, "static.stdout.txt").write_text("*ERROR synthetic input failure\n")
        Path(cwd, "static.stderr.txt").write_text("")
        return SimpleNamespace(returncode=201), {"failure": None, "returncode": 201}
    monkeypatch.setattr(old.fem1, "run_job", failed_job)
    result = module.run_resume(resume_config, tmp_path, summary, parent_context[0])
    assert len(calls) == result["job_calls"]["ccx"] == 1
    assert set(result["cases"]) == {"medium"}
    assert set(result["cases"]["medium"]) == {"linear"}
    assert result["cases"]["medium"]["linear"]["status"] == "FAILED_SOLVER"
    assert result["attempt_ledger"][0]["case"] == "medium_linear"
    assert result["attempt_ledger"][0]["status"] == "FAILED_SOLVER"
    assert result["resume_statuses"]["NLSP_FEM2R_MEDIUM_LINEAR"] == "FAIL"
    assert result["resume_statuses"]["NLSP_FEM2R_MEDIUM_NONLINEAR"] == "NOT_RUN"
    assert result["resume_statuses"]["NLSP_FEM2R_FINE_PAIR"] == "NOT_RUN"
    assert result["resume_statuses"]["NLSP_FEM2R_REFINED_PAIR"] == "NOT_RUN"
    command, timeout, memory, env = calls[0]
    assert command == [parent_context[2]["ccx_exe"], "static"]
    assert timeout <= 1200 and memory <= 4*1024**3
    assert all(env[k] == "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBER_OF_CPUS"))
    module.run_resume(resume_config, tmp_path, result, parent_context[0])
    assert len(calls) == 1 and len(result["attempt_ledger"]) == 1
    assert not (tmp_path/"cases/medium/nonlinear").exists()
    assert not (tmp_path/"cases/fine").exists() and not (tmp_path/"cases/refined").exists()


def test_medium_nonlinear_failure_stops_finer_pairs(mesh_tmp, resume_config, parent_context, monkeypatch):
    tmp_path = mesh_tmp
    module = _resume(); summary = _fresh_summary(resume_config, parent_context)
    summary["cases"] = {"medium": {"linear": {"status": "PASS"}}}
    summary["attempt_ledger"] = [{"case": "medium_linear", "ordinal": 1, "status": "PASS"}]
    calls = []
    def failed_job(command, cwd, timeout, memory, prefix, env):
        calls.append(command)
        Path(cwd, "static.stdout.txt").write_text("*ERROR synthetic nonlinear failure\n")
        Path(cwd, "static.stderr.txt").write_text("")
        return SimpleNamespace(returncode=201), {"failure": None, "returncode": 201}
    monkeypatch.setattr(old.fem1, "run_job", failed_job)
    result = module.run_resume(resume_config, tmp_path, summary, parent_context[0])
    assert len(calls) == 1
    assert result["attempt_ledger"][-1]["case"] == "medium_nonlinear"
    assert result["cases"]["medium"]["nonlinear"]["status"] == "FAILED_SOLVER"
    assert set(result["cases"]) == {"medium"}
    assert result["resume_statuses"]["NLSP_FEM2R_FINE_PAIR"] == "NOT_RUN"
    assert result["resume_statuses"]["NLSP_FEM2R_REFINED_PAIR"] == "NOT_RUN"


def test_scoped_parent_profile_callback_restores_old_helper(parent_context):
    module = _resume(); original = old.load_one_d_static_profile
    with module.parent_profile_loader(parent_context[0]):
        linear = old.load_one_d_static_profile(Path("no-local-copy-needed"), 64, False)
        nonlinear = old.load_one_d_static_profile(Path("no-local-copy-needed"), 64, True)
        assert linear["fields"][500, 1] == pytest.approx(.005, rel=1e-12)
        assert nonlinear["fields"][500, 1] < linear["fields"][500, 1]
    assert old.load_one_d_static_profile is original



@pytest.fixture
def mesh_tmp(tmp_path):
    """Relative source includes require test decks on the source mesh volume."""
    if tmp_path.drive.lower() == ROOT.drive.lower():
        yield tmp_path
        return
    import tempfile
    directory = ROOT / "Temp"
    directory.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="fem2r_contract_", dir=directory) as owned:
        resolved = Path(owned).resolve()
        assert resolved.is_relative_to(ROOT.resolve())
        yield resolved


@pytest.fixture
def stopped_resume_cache(tmp_path, resume_config, parent_context, monkeypatch):
    module = _resume()
    monkeypatch.setattr(module, "OUTPUT", tmp_path)
    bundle = tmp_path / "historical_code_hash"
    bundle.mkdir()
    _, item = module.identity(RESUME_CONFIG)
    item = deepcopy(item)
    item["code_sha256"] = "0"*64
    summary = _fresh_summary(resume_config, parent_context)
    summary["cases"] = {"medium": {"linear": {"status": "FAILED_SOLVER", "failure": "synthetic fixture"}}}
    summary["job_calls"]["ccx"] = 1
    summary["attempt_ledger"] = [{"case": "medium_linear", "ordinal": 1, "status": "FAILED_SOLVER"}]
    module.write_json(bundle/"provenance.json", item)
    module.save_progress(bundle, summary)
    module.finalize_manifest(bundle, item)
    return bundle, item, summary


def test_identity_pins_parent_load_sources_binary_code_and_environment(resume_config):
    module = _resume()
    first, item = module.identity(RESUME_CONFIG)
    second, repeated = module.identity(RESUME_CONFIG)
    assert first == second and item == repeated
    assert item["parent_manifest_sha256"] == resume_config["parent_manifest_sha256"]
    assert item["generator_sha256"] == resume_config["corrected_generator_sha256"]
    assert item["authorization"] == resume_config["authorization"]
    assert item["source_1D_sha256"] == {f"one_d_p{p}.npz": old.sha(PARENT/f"one_d_p{p}.npz") for p in (48, 64)}
    assert item["ccx_sha256"] == old.sha(item["config"]["ccx_exe"])
    assert item["code_sha256"] == old.sha(module.__file__)
    assert set(item["dependencies"]) == {"numpy", "scipy", "matplotlib"}
    assert item["sources"] == item["config"]["sources"]
    assert item["continuation_config"]["frozen_load"] == resume_config["frozen_load"]


def test_authorization_cache_replays_failure_despite_different_code_hash(stopped_resume_cache, resume_config):
    bundle, item, summary = stopped_resume_cache
    assert item["code_sha256"] != old.sha(_resume().__file__)
    found, result = _resume().find_authorized_attempt(resume_config, bundle.parent)
    assert found == bundle and result == summary
    assert result["attempt_ledger"][0]["status"] == "FAILED_SOLVER"


@pytest.mark.parametrize("mode", ["--check-source", "--run-fem", "--report-only", "--plot-only"])
def test_cached_routes_have_zero_solver_or_preparation_calls(stopped_resume_cache, resume_config, monkeypatch, mode):
    module = _resume()
    bundle, _, expected = stopped_resume_cache
    def forbidden(*args, **kwargs):
        raise AssertionError("Cached route must not execute or prepare a new case")
    monkeypatch.setattr(module, "run_resume", forbidden)
    monkeypatch.setattr(module, "input_gate", forbidden)
    monkeypatch.setattr(module, "identity", forbidden)
    monkeypatch.setattr(old, "build_static_preflight", forbidden)
    monkeypatch.setattr(old, "fem2_static_newton", forbidden)
    monkeypatch.setattr(old.fem1, "run_job", forbidden)
    monkeypatch.setattr(old.fem1.single, "generate_mesh_with_gmsh_cli", forbidden)
    monkeypatch.setattr(old.fem2_rod, "derive_polynomials", forbidden)
    arguments = [mode, str(bundle)] if mode in ("--report-only", "--plot-only") else [mode, "--output-dir", str(bundle.parent)]
    before = old.sha(bundle/"manifest.json")
    result = module.main(arguments)
    assert result == expected
    assert result["job_calls"]["ccx"] == 1
    assert before == old.sha(bundle/"manifest.json")


def test_resume_cache_detects_own_artifact_corruption(stopped_resume_cache):
    bundle, _, _ = stopped_resume_cache
    (bundle/"summary.json").write_text("{}\n")
    with pytest.raises(ValueError, match="artifact"):
        _resume().validate_cache(bundle)


def test_unmanifested_attempt_cannot_be_retried(tmp_path, resume_config, monkeypatch):
    module = _resume()
    monkeypatch.setattr(module, "OUTPUT", tmp_path)
    bundle = tmp_path/"interrupted"
    bundle.mkdir()
    module.write_json(bundle/"provenance.json", {
        "authorization": resume_config["authorization"],
        "parent_manifest_sha256": resume_config["parent_manifest_sha256"],
        "continuation_config": resume_config})
    with pytest.raises(RuntimeError, match="no new solver attempt"):
        module.find_authorized_attempt(resume_config, tmp_path)


def test_complete_execution_cannot_hide_partial_signal_resolution(resume_config, parent_context):
    summary = _fresh_summary(resume_config, parent_context)
    summary["cases"] = {level: {kind: {"status": "PASS"} for kind in ("linear", "nonlinear")}
                        for level in ("medium", "fine", "refined")}
    summary["completed_levels"] = ["medium", "fine", "refined"]
    summary["statuses"].update(NLSP_FEM2_NONLINEAR_CORRECTION_MESH_CHECK="PARTIAL",
                               NLSP_FEM2_1D_3D_COMPARISON="PARTIAL")
    _resume().update_resume_statuses(summary)
    assert summary["overall"] == "FEM2_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS"
    assert summary["resume_statuses"]["NLSP_FEM2R_NONLINEAR_SIGNAL_RESOLUTION"] == "PARTIAL"
    assert summary["resume_statuses"]["NLSP_FEM2R_1D_3D_COMPARISON"] == "PARTIAL"


def test_finished_output_can_be_reparsed_without_solver_retry(tmp_path, resume_config, parent_context, monkeypatch):
    module = _resume()
    summary = _fresh_summary(resume_config, parent_context)
    case = tmp_path/"cases/medium/linear"
    case.mkdir(parents=True)
    _, mesh, audit, _, _, _, _ = _static_outputs(case, False)
    summary["cases"] = {"medium": {"linear": {"status": "RECOVERY_PENDING", "solver_finished": True}}}
    def forbidden(*args, **kwargs):
        raise AssertionError("Read-only recovery must not rerun CCX")
    monkeypatch.setattr(old.fem1, "run_job", forbidden)
    monkeypatch.setattr(old, "recover_static_sections", lambda *a, **k: {
        "x": [0., 1.], "fields": np.zeros((2, 7)).tolist(), "recovery_policy": old.FEM2_RECOVERY_POLICY})
    assert module.recover_case(resume_config, tmp_path, summary, "medium", "linear", mesh, audit)
    assert summary["cases"]["medium"]["linear"]["status"] == "PASS"
    assert summary["job_calls"]["ccx"] == 0
    assert (case/"static_nodal_results.npz").exists()
    assert (case/"static_diagnostics.json").exists()


def test_parser_exception_remains_pending_and_does_not_fabricate_equilibrium(tmp_path, resume_config, parent_context, monkeypatch):
    module = _resume(); summary = _fresh_summary(resume_config, parent_context)
    case = tmp_path/"cases/medium/linear"; case.mkdir(parents=True)
    summary["cases"] = {"medium": {"linear": {"status": "RECOVERY_PENDING", "solver_finished": True}}}
    def unsupported(*args, **kwargs):
        raise ValueError("Synthetic unsupported output record")
    monkeypatch.setattr(module, "parse_static_outputs", unsupported)
    assert not module.recover_case(resume_config, tmp_path, summary, "medium", "linear", _mesh(), _audit())
    record = summary["cases"]["medium"]["linear"]
    assert record["status"] == "RECOVERY_PENDING"
    assert record["read_only_reparse_permitted"]
    assert summary["resume_statuses"]["NLSP_FEM2R_MEDIUM_LINEAR"] == "PARTIAL"
    assert summary["job_calls"]["ccx"] == 0


def test_measured_equilibrium_failure_is_not_reclassified_as_parser_pending(tmp_path, resume_config, parent_context, monkeypatch):
    module = _resume(); summary = _fresh_summary(resume_config, parent_context)
    case = tmp_path/"cases/medium/linear"; case.mkdir(parents=True)
    summary["cases"] = {"medium": {"linear": {"status": "RECOVERY_PENDING", "solver_finished": True}}}
    monkeypatch.setattr(module, "parse_static_outputs", lambda *a, **k: (
        {"status": "FAIL", "failures": ["SUPPORT_FORCE_BALANCE"]}, {"U": np.zeros((10, 3))}))
    assert not module.recover_case(resume_config, tmp_path, summary, "medium", "linear", _mesh(), _audit())
    record = summary["cases"]["medium"]["linear"]
    assert record["status"] == "FAIL"
    assert not record.get("read_only_reparse_permitted", False)
    assert summary["resume_statuses"]["NLSP_FEM2R_EQUILIBRIUM"] == "FAIL"
    assert summary["job_calls"]["ccx"] == 0


@pytest.fixture(scope="module")
def actual_resume_evidence(resume_config):
    module = _resume()
    candidates = []
    if module.OUTPUT.is_dir():
        for directory in module.OUTPUT.iterdir():
            if not directory.is_dir() or not (directory/"manifest.json").exists():
                continue
            manifest = module.read_json(directory/"manifest.json")
            if manifest["identity"].get("continuation_config") == resume_config:
                candidates.append(directory)
    if not candidates:
        pytest.skip("No actual FEM-2R artifact present; tests must not recreate it")
    assert len(candidates) == 1, "One explicit authorization must have only one continuation bundle"
    bundle = candidates[0]
    return bundle, module.validate_cache(bundle)


def test_actual_new_attempt_is_separate_from_old_input_failure(actual_resume_evidence):
    bundle, summary = actual_resume_evidence
    assert bundle != PARENT
    assert summary["parent_reference"]["path"] == str(PARENT.relative_to(ROOT)).replace("\\", "/")
    assert summary["job_calls"]["one_d_static"] == summary["job_calls"]["gmsh"] == summary["job_calls"]["modal"] == summary["job_calls"]["nonlinear_ODE"] == 0
    assert old.read_json(PARENT/"summary.json")["cases"]["medium"]["linear"]["job"]["returncode"] == 201


def test_actual_corrected_medium_linear_has_full_load_and_complete_fields(actual_resume_evidence):
    bundle, summary = actual_resume_evidence
    record = summary["cases"].get("medium", {}).get("linear", {})
    if record.get("status") != "PASS":
        pytest.skip("Actual medium linear not accepted; no extra job in test")
    diag = record["diagnostics"]
    assert record["solver_finished"] and record["job"]["returncode"] == 0
    assert diag["status"] == "PASS" and diag["failures"] == []
    assert diag["final_dat"]["time"] == diag["final_frd"]["time"] == 1.
    assert diag["actual_NLGEOM"] is False
    assert diag["maximum_fixed_displacement"] == 0.
    with np.load(bundle/"cases/medium/linear/static_nodal_results.npz", allow_pickle=False) as arrays:
        assert arrays["U"].shape == (record["mesh_audit"]["nodes"], 3)
        assert np.all(np.isfinite(arrays["U"]))
        assert np.min(arrays["U"][:, 1]) < 0.


def test_actual_attempt_count_order_and_resource_budget(actual_resume_evidence):
    _, summary = actual_resume_evidence
    attempts = summary["attempt_ledger"]
    assert [row["case"] for row in attempts] == list(_resume().ORDER[:len(attempts)])
    assert [row["ordinal"] for row in attempts] == list(range(1, len(attempts)+1))
    assert 0 < summary["job_calls"]["ccx"] <= 6
    assert summary["runtime"]["numerical_seconds"] <= 3600
    for rows in summary["cases"].values():
        for record in rows.values():
            if "job" in record:
                assert record["job"]["seconds"] <= 1200
                assert record["job"]["peak_working_set_bytes"] <= 4*1024**3


def test_actual_all_completed_solutions_keep_source_mesh_and_frozen_load(actual_resume_evidence, parent_context):
    _, summary = actual_resume_evidence
    sources = old.load_fem2_sources(parent_context[2])
    for level, rows in summary["cases"].items():
        source, mesh, audit = old.fem2_source_mesh(parent_context[2], level, sources)
        for kind, record in rows.items():
            assert record["source_mesh_include_sha256"] == old.sha(source/"solid_mesh.inp")
            assert record["mesh_audit"]["solid_element_types"] == ["C3D10"]
            assert record["mesh_audit"]["nodes"] == len(mesh.nodes)
            if record.get("status") == "PASS":
                assert record["diagnostics"]["equilibrium_relative_gate"] == 1e-5
                assert record["diagnostics"]["support_force_balance_relative"] <= 1e-5
                assert record["diagnostics"]["support_moment_balance_relative"] <= 1e-5
                np.testing.assert_allclose(record["diagnostics"]["applied_force_global"], [0., -2.844950213371267e-5, 0.], rtol=1e-12, atol=1e-18)


def test_actual_nonlinear_load_history_and_routing(actual_resume_evidence):
    _, summary = actual_resume_evidence
    nonlinear = [rows["nonlinear"] for rows in summary["cases"].values()
                 if rows.get("nonlinear", {}).get("status") == "PASS"]
    if not nonlinear:
        pytest.skip("No actual nonlinear equilibrium yet; no automatic solve")
    for record in nonlinear:
        diag = record["diagnostics"]
        assert diag["actual_NLGEOM"]
        assert diag["load_history"]["status"] == "PARSED"
        assert diag["load_history"]["last_time"] == 1.
        assert diag["load_history"]["accepted_increments"][-1]["total_time"] == 1.
        assert diag["solver_warnings"] == []


def test_actual_section_profiles_use_material_geometry_and_polar_policy(actual_resume_evidence):
    bundle, summary = actual_resume_evidence
    count = 0
    for level, rows in summary["cases"].items():
        for kind, record in rows.items():
            if record.get("status") != "PASS":
                continue
            count += 1
            profile = old.read_json(bundle/"cases"/level/kind/"recovered_sections.json")
            assert profile["reference_geometry_used_for_sections"]
            assert not profile["additional_derivative_constraints"]
            assert "both_linear_and_nonlinear" in profile["theta_policy"]
            assert profile["contraction_status"].endswith("NOT_MH_DOF")
            assert np.all(np.isfinite(np.asarray(profile["fields"])))
            assert profile["recovery_sensitivity"]["alternate_section_count"] == 81
            assert len(profile["recovery_sensitivity"]["alternate_x"]) == 83  # plus the two clamped face rows
    assert count >= 1


def test_actual_comparison_arithmetic_preserves_signed_corrections(actual_resume_evidence):
    _, summary = actual_resume_evidence
    rows = summary.get("comparison_rows", [])
    if not rows:
        pytest.skip("No completed actual pair to compare")
    for record in rows:
        assert record["FEM_NL_midspan"] - record["FEM_linear_midspan"] == pytest.approx(record["FEM_delta_w_midspan"], abs=1e-16)
        assert record["one_D_NL_midspan"] - record["one_D_linear_midspan"] == pytest.approx(record["one_D_delta_w_midspan"], abs=1e-16)
        assert record["FEM_nonlinear_effect"] == pytest.approx(record["FEM_delta_w_midspan"]/record["FEM_linear_midspan"], rel=1e-12)
        assert record["one_D_nonlinear_effect"] == pytest.approx(record["one_D_delta_w_midspan"]/record["one_D_linear_midspan"], rel=1e-12)


def test_actual_finished_series_retains_small_signal_uncertainty_qualification(actual_resume_evidence):
    _, summary = actual_resume_evidence
    if summary.get("completed_levels") != ["medium", "fine", "refined"]:
        pytest.skip("Series is an actual prefix; tests must not complete it")
    assert summary["job_calls"]["ccx"] == 6
    assert len(summary["mesh_changes"]) == 2
    assert [(r["from"], r["to"]) for r in summary["mesh_changes"]] == [("medium", "fine"), ("fine", "refined")]
    resolution = summary["signal_resolution"]
    assert not resolution["continuum_error_bound_claimed"]
    assert resolution["recovery_difference"]["absolute_max"] >= 0.
    assert summary["overall"] == "FEM2_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS"
    status = summary["resume_statuses"]["NLSP_FEM2R_NONLINEAR_SIGNAL_RESOLUTION"]
    assert status in ("PASS", "PARTIAL")


