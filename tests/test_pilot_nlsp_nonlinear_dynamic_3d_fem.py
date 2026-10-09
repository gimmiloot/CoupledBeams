"""Scoped FEM-3A orchestration tests; all native/ODE jobs are mocked or forbidden."""
import copy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from scripts.analysis import pilot_nlsp_nonlinear_dynamic_3d_fem as cli

ROOT = Path(__file__).resolve().parents[1]
PREFLIGHT = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_pilot/a69310e3bb30bab7"


@pytest.fixture
def config():
    return cli.read_json(cli.CONFIG)


@pytest.fixture
def decks():
    return {"linear": (PREFLIGHT / "remediation_preview" / "linear_corrected_NOT_RUN.inp").read_text(encoding="utf8"),
        "nonlinear": (PREFLIGHT / "input_gate" / "nonlinear.inp").read_text(encoding="utf8")}


def _forbid(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("report/preflight must perform zero scientific calls")
    monkeypatch.setattr(cli, "run_case", forbidden)
    monkeypatch.setattr(cli, "prepare", forbidden)
    monkeypatch.setattr(cli.one, "integrate_nonlinear_reference", forbidden)
    monkeypatch.setattr(cli.one, "exact_linear_reference", forbidden)
    monkeypatch.setattr(cli.one.dynamics.PlanarGalerkin, "linear_eigenpairs", forbidden)
    monkeypatch.setattr(cli.one.rod, "derive_polynomials", forbidden)
    monkeypatch.setattr(cli.base, "fem2_static_newton", forbidden)
    monkeypatch.setattr(cli.base.fem1, "run_job", forbidden)
    monkeypatch.setattr(cli.base.fem1.single, "generate_mesh_with_gmsh_python", forbidden)
    monkeypatch.setattr(cli.base.fem1.single, "generate_mesh_with_gmsh_cli", forbidden)


def test_exact_frozen_geometry_material_load_and_horizon(config):
    assert cli.validate_config(config) is config
    assert config["geometry"] == {"L": 1., "b": .2, "h": .1}
    assert config["g"] == .0014224751066856333
    assert config["q"] == 2.844950213371267e-5
    settings = cli.dynamic_settings(config)
    assert settings["T1"] == pytest.approx(10.37828159055014, rel=1e-15)
    assert settings["duration"] == .05 * settings["T1"]
    assert settings["initial_increment"] == settings["T1"] / 4000
    assert settings["maximum_increment"] == settings["T1"] / 2000
    assert settings["minimum_increment"] == settings["initial_increment"] * 1e-4


@pytest.mark.parametrize("key,value", [("g", .001), ("q", 1e-5), ("horizon_T1", .1),
    ("mesh_level", "fine"), ("omega1", .3174742907880648), ("threads", 2),
    ("job_timeout_seconds", 1201), ("execution_mode", "STRICT_ADMITTED"), ("admitted", True)])
def test_scope_changes_rejected(config, key, value):
    config[key] = value
    with pytest.raises(ValueError):
        cli.validate_config(config)


def test_authorization_and_no_retry_contract(config):
    assert config["authorization"]["maximum_production_CCX_jobs"] == 2
    assert config["authorization"]["maximum_nonlinear_1D_integrations"] == 1
    assert config["authorization"]["automatic_retry"] is False
    config["authorization"]["automatic_retry"] = True
    with pytest.raises(ValueError, match="authorization"):
        cli.validate_config(config)


def test_linear_nonlinear_decks_only_route_nlgeom(config, decks):
    def normalize(text):
        first,rest=text.split("*END STEP",1)
        return (first.replace("ELSE,ELKE","ELSE")+"*END STEP"+rest).replace(", NLGEOM=NO", "").replace(", NLGEOM", "")
    assert normalize(decks["linear"]) == normalize(decks["nonlinear"])
    for kind, text in decks.items():
        assert text.count("*STEP,") == 2
        assert text.count("*END STEP") == 2
        assert text.count("*STATIC") == 1
        assert text.count("*DYNAMIC, ALPHA=0") == 1
        assert cli.input_contract(text, config)["status"] == "PASS"


def test_initial_velocity_and_static_state_transfer_contract(decks):
    for text in decks.values():
        assert "*INITIAL CONDITIONS, TYPE=VELOCITY\nALL_NODES,1,0.\nALL_NODES,2,0.\nALL_NODES,3,0." in text
        assert text.index("*STATIC") < text.index("*DYNAMIC")
        assert text.count("*INCLUDE") == 1
        assert "*INITIAL CONDITIONS, TYPE=DISPLACEMENT" not in text
        assert "*RESTART" not in text
        for forbidden in ("*MPC", "*SPRING", "*CONTACT", "*DAMPING", "*CLOAD", "*FREQUENCY", "*MODAL DYNAMIC"):
            assert forbidden not in text.upper()


def test_preload_same_physics_and_output_extension_only(config, decks):
    old = cli.read_json(ROOT / config["source_resume"]["bundle"] / "frozen_config.json")
    science = old["science_config"] if "science_config" in old else cli.read_json(ROOT / config["source_resume"]["bundle"] / "summary.json")["science_config"]
    assert science["geometry"] == config["geometry"]
    assert science["material"] == config["material"]
    for text in decks.values():
        preload = text[:text.index("** Free motion")]
        load_line = next(line for line in preload.splitlines() if line.startswith("SOLID,GRAV,"))
        parts = [value.strip() for value in load_line.split(",")]
        assert float(parts[2]) == pytest.approx(config["g"], rel=1e-12)
        assert list(map(float, parts[3:])) == [0., -1., 0.]
        assert "*BOUNDARY" in preload
        assert "LEFT_FIXED,1,3" in preload and "RIGHT_FIXED,1,3" in preload
        assert "*EL PRINT, ELSET=SOLID, TOTALS=ONLY, FREQUENCY=1\nELSE" in preload
        assert ("ELKE" in preload) == ("NLGEOM" in next(row for row in preload.splitlines() if row.startswith("*STEP")))
        assert "S,E,ENER" in preload


def test_numeric_fields_all_twenty_characters(config, decks):
    for text in decks.values():
        cli.resume.numeric_cards(text)
        gate = cli.input_contract(text, config)
        assert all(len(token) <= 20 for token in gate["dynamic_native_numeric_fields"])
        assert all(np.isfinite(float(token)) for token in gate["dynamic_native_numeric_fields"])
        assert gate["maximum_numeric_width"] <= 20
    bad = decks["linear"].replace("*DYNAMIC, ALPHA=0\n", "*DYNAMIC, ALPHA=0\n" + "0.0000000000000000000000," ,1)
    with pytest.raises(ValueError):
        cli.input_contract(bad, config)


@pytest.mark.parametrize("replacement", ["*DLOAD, OP=MOD", "*DLOAD"])
def test_explicit_old_bodyload_removal_required(config, decks, replacement):
    with pytest.raises(ValueError, match="release/alpha"):
        cli.input_contract(decks["linear"].replace("*DLOAD, OP=NEW", replacement), config)


def test_step_release_zero_grav_alpha_no_damping(config, decks):
    for text in decks.values():
        dynamic = text[text.index("** Free motion"):]
        assert "AMPLITUDE=STEP" in dynamic
        assert "*DLOAD, OP=NEW\nSOLID,GRAV,0.,0.,-1.,0." in dynamic
        assert "*DYNAMIC, ALPHA=0" in dynamic
        gate = cli.input_contract(text, config)
        assert gate["external_GRAV_after_release"] == 0
        assert gate["alpha"] == 0
    with pytest.raises(ValueError, match="release/alpha"):
        cli.input_contract(decks["linear"].replace("ALPHA=0", "ALPHA=-0.05"), config)


def test_parent_hashes_preserved_by_source_audit(config):
    before = {key: cli.sha(ROOT / config[key]["bundle"] / "manifest.json")
        for key in ("source_resume", "source_static", "source_fem1", "source_action")}
    old, source, mesh, audit = cli.verify_sources(config)
    assert len(mesh.nodes) == 5649 and len(mesh.solid_elements) == 3120
    assert audit["status"] == "PASS"
    assert source.name == "medium"
    assert old["completed_levels"] == ["medium", "fine", "refined"]
    after = {key: cli.sha(ROOT / config[key]["bundle"] / "manifest.json") for key in before}
    assert before == after
    assert all(before[key] == config[key]["manifest_sha256"] for key in before)


def test_source_manifest_corruption_stops_before_any_job(config, monkeypatch):
    _forbid(monkeypatch)
    config["source_resume"]["manifest_sha256"] = "corrupt"
    with pytest.raises(ValueError, match="Source manifest mismatch"):
        cli.verify_sources(config)


def test_saved_initial_q_and_release_preflight(config):
    pre = cli.read_json(PREFLIGHT / "one_d_preflight.json")
    with np.load(ROOT / config["source_static"]["bundle"] / "one_d_p64.npz") as data:
        for kind in ("linear", "nonlinear"):
            row = pre["initial_states"][kind]
            np.testing.assert_array_equal(row["q0"], data["q_" + kind])
            assert np.count_nonzero(row["v0"]) == 0
            assert row["midspan_w_acceleration"] < 0
            assert row["mass_cholesky_positive"]
            assert row["released_action_relative_residual"] < 2e-12
    assert pre["strict_float64_strong_weak"] == "PARTIAL"
    assert pre["admitted"] is False
    assert pre["slope_constraints"] is False


@pytest.mark.parametrize("flag", ["--preflight", "--report-only", "--plot-only"])
def test_cached_preflight_report_plot_zero_scientific_calls(flag, monkeypatch):
    _forbid(monkeypatch)
    args = [flag] if flag == "--preflight" else [flag, str(PREFLIGHT)]
    result = cli.main(args)
    assert result["execution_mode"] == "EXPLORATORY_NOT_CERTIFIED"
    assert result["admitted"] is False
    # A later real pilot may update this bundle; cached calls still do zero work.
    assert result["overall"] in ("NOT_RUN", "PARTIAL", "BLOCKED_BY_SOLVER", "PILOT_COMPLETE_WITH_QUALIFICATIONS")


def test_preflight_manifest_verifies_without_solves(monkeypatch):
    _forbid(monkeypatch)
    summary = cli.validate_cache(PREFLIGHT)
    manifest = cli.read_json(PREFLIGHT / "manifest.json")
    assert manifest["identity"]["config"]["horizon_T1"] == .05
    assert "input_gate/linear.inp" in manifest["artifact_hashes"]
    assert "one_d_preflight.json" in manifest["artifact_hashes"]
    assert summary["strict_float64_qualification"] == "PARTIAL"


def test_actual_timestamp_matching_without_interpolation():
    all_times = np.array([0., .001, .003, .01])
    np.testing.assert_array_equal(cli._matches(all_times, np.array([.001, .01])), [1, 3])
    with pytest.raises(ValueError, match="Actual-time"):
        cli._matches(all_times, np.array([.0015]))
    assert cli._rounded_endpoint(.01 + 1e-9, .01) == .01
    assert cli._rounded_endpoint(.01 + 1e-5, .01) != .01


def test_nonlinear_correction_common_scale_and_coordinates():
    x = np.linspace(0., 1., 5)
    linear = np.array([[0., .1, .2, .1, 0.], [0., .09, .18, .09, 0.]])
    nonlinear = .99 * linear
    delta = nonlinear - linear
    compared = cli.field_difference(delta, delta * 1.1, x, scale=.0022)
    assert compared["absolute_max"] == pytest.approx(.0002)
    assert compared["relative_max"] == pytest.approx(.0002 / .0022)
    assert compared["phase_amplitude_alignment"] is False
    assert compared["sampled_maxima_only"] is True
    assert compared["x_at_max"] == .5


def test_first_native_failure_stops_and_no_retry(config, tmp_path, monkeypatch):
    summary = {"cases": {}, "attempts": [], "job_calls": {"CCX_production": 0},
        "numerical_seconds": 0., "overall": "NOT_RUN", "statuses": {}}
    calls = []
    monkeypatch.setattr(cli, "verify_sources", lambda c:
        ({"science_config": {"ccx_exe": "synthetic_ccx"}}, tmp_path, None, None))
    def input_only(path, *a):
        Path(path).write_text("synthetic-only deck", encoding="utf8")
        return {"status": "PASS"}
    monkeypatch.setattr(cli, "write_input", input_only)
    monkeypatch.setattr(cli, "finalize", lambda *a: None)
    monkeypatch.setattr(cli, "recover_case", lambda *a: pytest.fail("no output recovery after solver fail"))
    def failed(command, directory, *a):
        calls.append(command)
        (directory / "motion.stdout.txt").write_text("*ERROR reading fixture input")
        return SimpleNamespace(returncode=1), {"failure": "synthetic_native_failure"}
    monkeypatch.setattr(cli.base.fem1, "run_job", failed)
    assert cli.run_case(tmp_path, config, {}, summary, "linear") is False
    assert summary["overall"] == "BLOCKED_BY_SOLVER"
    assert summary["cases"]["linear"]["status"] == "FAIL"
    assert summary["attempts"][0]["status"] == "FAIL"
    assert summary["job_calls"]["CCX_production"] == 1
    assert cli.run_case(tmp_path, config, {}, summary, "linear") is False
    assert len(calls) == 1
    assert "nonlinear" not in summary["cases"]


def test_medium_linear_gate_prevents_nonlinear_after_failure(config, tmp_path, monkeypatch):
    summary = {"cases": {}, "attempts": [], "overall": "NOT_RUN", "statuses": {}, "job_calls": {}, "numerical_seconds": 0.}
    monkeypatch.setattr(cli, "existing_attempt", lambda c: (tmp_path, {}, summary))
    calls = []
    monkeypatch.setattr(cli, "run_case", lambda *a: calls.append(a[-1]) or False)
    monkeypatch.setattr(cli, "finalize", lambda *a: None)
    monkeypatch.setattr(cli, "finish_references_and_comparison", lambda *a: pytest.fail("no 1D after medium failure"))
    cli.main(["--run-pilot"])
    assert calls == ["linear"]


def test_blocked_cached_attempt_cannot_retry(config, tmp_path, monkeypatch):
    _forbid(monkeypatch)
    summary = {"overall": "BLOCKED_BY_SOLVER"}
    monkeypatch.setattr(cli, "existing_attempt", lambda c: (tmp_path, {}, summary))
    assert cli.main(["--run-pilot"]) is summary



def test_native_linear_static_elke_defect_is_guarded(config):
    failed=(PREFLIGHT/'cases/linear/motion.inp').read_text(encoding='utf8')
    with pytest.raises(ValueError,match='ELKE is unsafe in linear STATIC'):
        cli.input_contract(failed,config)
    assert cli.read_json(PREFLIGHT/'native_failure_diagnosis.json')['status']=='LOCALIZED_NATIVE_LINEAR_STATIC_ELKE_FREED_VELOCITY_READ'


def test_corrected_preview_retains_dynamic_energy_and_frozen_physics(config,decks):
    first,dynamic=decks['linear'].split('*END STEP',1)
    assert 'ELKE' not in first
    assert 'S,E,ENER' in first and '\nELSE\n' in first
    assert 'ELSE,ELKE' in dynamic
    original=(PREFLIGHT/'cases/linear/motion.inp').read_text(encoding='utf8')
    strip_include=lambda t:'\n'.join(row for row in t.splitlines() if not row.startswith('*INCLUDE'))
    expected=original.replace('ELSE,ELKE\n*END STEP','ELSE\n*END STEP',1)
    assert strip_include(expected)==strip_include(decks['linear'])
    assert cli.read_json(PREFLIGHT/'remediation_preview/input_gate.json')['execution_status']=='NOT_RUN'


def test_actual_failed_attempt_has_no_fabricated_prefix():
    s=cli.read_json(PREFLIGHT/'summary.json')
    assert s['overall']=='BLOCKED_BY_SOLVER'
    assert s['cases']['linear']['status']=='FAIL'
    assert s['statuses']['NLSP_FEM3A_3D_LINEAR_DYNAMIC']=='FAIL'
    assert s['job_calls']['CCX_production']==1
    assert s['job_calls']['1D_nonlinear_ODE']==0
    assert 'nonlinear' not in s['cases']
    assert (PREFLIGHT/'cases/linear/motion.dat').stat().st_size==0
    assert (PREFLIGHT/'cases/linear/motion.sta').stat().st_size==0
    assert not (PREFLIGHT/'one_d_nonlinear.npz').exists()
    assert not (PREFLIGHT/'figures').exists()


def test_actual_execution_code_and_failed_input_are_preserved():
    m=cli.read_json(PREFLIGHT/'provenance.json')
    name='scripts/analysis/pilot_nlsp_nonlinear_dynamic_3d_fem.py'
    assert cli.sha(PREFLIGHT/'execution_code/pilot_nlsp_nonlinear_dynamic_3d_fem.py')==m['helper_sha256'][name]
    assert cli.sha(cli.__file__)!=m['helper_sha256'][name]
    assert cli.sha(PREFLIGHT/'cases/linear/motion.inp')=='a38fb741b7638103bedee3893280266d32053cdae1b50d66804dd99369650998'


@pytest.mark.parametrize('args',[['--preflight'],['--run-pilot'],['--report-only',str(PREFLIGHT)],['--plot-only',str(PREFLIGHT)]])
def test_actual_failed_attempt_replay_survives_code_change_without_retry(args,monkeypatch):
    _forbid(monkeypatch)
    assert cli.main(args)['overall']=='BLOCKED_BY_SOLVER'
