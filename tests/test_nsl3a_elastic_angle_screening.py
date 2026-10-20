"""Bounded ANGLE diagnostics; no new eigenvalue calculation in tests."""
import numpy as np
import pytest
from scripts.analysis.nonsymmetric_laminated_beams import screen_elastic_angles as run


@pytest.fixture(scope="module")
def context():
    config = run.read(run.CONFIG)
    return config, run.context(config)


@pytest.mark.parametrize("angle", [15, 30, 45, 60, 75])
def test_reflection_is_structural_not_assigned_to_same_local(context, angle):
    config, (old, _, props, ref, _, _) = context
    for case in config["cases"]:
        model = run.make_model(old, props, ref, case, angle)
        result = run.structural(model, case, config["policy"])
        assert result["status"] == "PASS"
        assert result["symmetry_claim"] == (case == "MIRROR")
        assert result["protected_subspace_jump"] == 0


def test_saved_beta30_predictors_and_physical_overlap(context):
    config, (old, data, props, ref, _, kv) = context
    for state in kv["states"]:
        case = state["case"]; model = run.make_model(old, props, ref, case, 30)
        saved = data["spectral"][case+"_SPRING"]["modes"][state["elastic_position"]-1]
        mode = model.evaluate(saved["omega"], np.array(saved["coefficients"]), 160)
        protected = case == "MIRROR" and state["elastic_position"] in (4, 6)
        if protected:
            reflection = run.parity(mode, model.properties)
            assert reflection["eta"] == -1 and reflection["mass_error"] < 1e-7
        values = run.indicators(saved["omega"], mode["mass"], mode["delta_psi"],
                                ref, config["illustrative_d_theta"], protected)
        assert values["alpha_pred"] == pytest.approx(state["alpha_pred"], rel=1e-10, abs=0)
        assert values["predictor_algebra_relative"] < 1e-12
        assert run.global_pair(mode, mode, model.beta, model.beta, model.properties) == pytest.approx(mode["mass"], rel=1e-13)


def test_global_mac_transforms_local_vectors_before_comparison(context):
    _, (old, _, props, ref, _, _) = context
    ma = run.make_model(old, props, ref, "SAME_LOCAL", 15)
    mb = run.make_model(old, props, ref, "SAME_LOCAL", 75)
    fields = []
    for model in (ma, mb):
        arms = []
        for i, prop in enumerate(model.properties, 1):
            global_field = np.array([[.2, 1., -.3], [.7, -.1, .4]])
            local = global_field@run.joint.global_kinematic_map(model.beta, i)
            arms.append(dict(states=np.column_stack([local, np.zeros((2,3))]),
                             weights=np.array([.5,.5]), x=np.array([.2,.8])))
        fields.append(dict(arms=arms))
    same = run.global_pair(fields[0], fields[0], ma.beta, ma.beta, ma.properties)
    transformed = run.global_pair(fields[0], fields[1], ma.beta, mb.beta, ma.properties)
    wrong = run.global_pair(fields[0], fields[1], ma.beta, ma.beta, ma.properties)
    assert transformed == pytest.approx(same, rel=1e-14)
    assert abs(wrong/same-1) > .05


def test_small_same_local_participation_is_not_zeroed(context):
    config, (_, _, _, ref, _, _) = context
    active = run.indicators(.7, 1., 1e-12, ref, .001, False)
    protected = run.indicators(.7, 1., 1e-12, ref, .001, True)
    assert active["P"] > 0 and active["alpha_pred"] > 0
    assert protected["P"] == 0 and protected["P_numeric"] > 0
    assert protected["alpha_pred"] == 0


def test_ambiguous_mac_never_becomes_branch_identity(context):
    config, _ = context
    values = np.eye(6)
    values[:2,:2] = [[.93,.91],[.91,.93]]
    assigned = run.assignments(values, config["policy"])
    assert all(r["status"] == "UNRESOLVED" for r in assigned[:2])
    assert all(r["status"] == "RESOLVED" for r in assigned[2:])


def test_completed_cache_never_calls_root_or_recovery(monkeypatch, capsys):
    def forbidden(*args, **kwargs):
        raise AssertionError("Completed cache must not solve or recover modes")
    monkeypatch.setattr(run.joint.CoupledBeams, "prefix", forbidden)
    monkeypatch.setattr(run.recovery, "recover_mode", forbidden)
    monkeypatch.setattr(run.sys, "argv", [str(run.__file__), "--compute"])
    run.main()
    assert '"new_eigenvalue_solves": 0' in capsys.readouterr().out


def test_saved_prefixes_guard_and_exact_beta30_coefficients(context):
    _, (_, source, _, _, _, _) = context
    out = run.OUTPUT/run.read(run.OUTPUT/"current.json")["fingerprint"]
    result = run.validate_cache(out)
    assert result["status"] == "PASS"
    assert len(result["table"]) == 60 and len(result["guard"]) == 10
    assert all(row["position"] <= 6 for row in result["table"])
    for point in result["points"]:
        assert len(point["modes"]) == 7
        assert [r["counts"] for r in point["count"]["roots"]] == [[j,j+1] for j in range(7)]
        if point["beta_deg"] == 30:
            saved = source["spectral"][point["case"]+"_SPRING"]["modes"]
            assert point["new_roots"] == 0
            for new, old in zip(point["modes"], saved):
                assert new["omega"] == old["omega"]
                assert new["coefficients"] == old["coefficients"]
