"""Read-only presentation contracts; no physical solvers or spectral calls."""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

SPEC = importlib.util.spec_from_file_location(
    "reflection_publication", Path(__file__).resolve().parents[1] /
    "scripts/analysis/plot_rl_full_rr_reflection_publication.py")
plot = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(plot)


def arm(multiplier=1.):
    y = np.zeros((2, 6))
    y[:, 1] = multiplier*np.array([1., 2.])
    return dict(x=np.array([.25, .75]), weights=np.array([.5, .5]),
                y=y, dy=np.zeros_like(y), ends=np.zeros_like(y))


def test_manifest_rejects_changed_numerical_data(tmp_path):
    folder = tmp_path / "archive"
    folder.mkdir()
    source = tmp_path / "fields.npz"
    source.write_bytes(b"original")
    (folder / "manifest.json").write_text(json.dumps(dict(
        status="COMPLETE", source_sha256={"fields.npz": plot.sha(source)}, artifacts={})))
    plot.verify_manifest(folder, root=tmp_path)
    source.write_bytes(b"changed")
    with pytest.raises(ValueError, match="SHA256"):
        plot.verify_manifest(folder, root=tmp_path)


def test_first_class_selection_does_not_use_full_sorted_rank():
    rows = [dict(symmetry_class=c, beta_deg="6.2", class_rank="1", Lambda=v,
                 original_sorted_position=p, prefix_complete_below="True", classification_status=s)
            for c, v, p, s in [("+", "3.937", "2", "PURE_PLUS"), ("-", "3.924", "1", "PURE_MINUS")]]
    selected = plot.first_classes(rows)
    assert selected["+"][0]["original_sorted_position"] == "2"
    assert selected["-"][0]["Lambda"] == "3.924"
    with pytest.raises(ValueError, match="duplicate"):
        plot.first_classes(rows + [rows[0]])
    rows[0]["classification_status"] = "UNRESOLVED_MIXED"
    with pytest.raises(ValueError, match="Unconfirmed"):
        plot.first_classes(rows)


def test_mass_retains_rotation_and_longitudinal_derivative():
    a = arm()
    a["y"][:] = 0
    a["y"][:, 2] = 2
    a["dy"][:, 0] = 3
    contract = dict(l=1., epsilon=.005, rho=1., nu=.3)
    inertia = (np.sqrt(12)*.005)**4/12
    expected = 2*(inertia*4 + .3**2*2*inertia*9)
    assert plot.mass_inner([a, a], [a, a], contract) == pytest.approx(expected)
    assert expected > 0


def test_tracking_embedding_uses_full_inertia_and_global_translations():
    a, b = arm(), arm(2)
    for v in (a, b):
        v["y"][:, 0] = [.5, 1.]
        v["y"][:, 2] = [3., 4.]
        v["dy"][:, 0] = [5., 6.]
    contract = dict(l=1., epsilon=.005, rho=1., nu=.3)
    tracker = plot.tracking_helpers()
    for beta in (0., 6.13, 90.):
        x = plot.mass_features([a, a], beta, contract)
        z = plot.mass_features([b, b], beta, contract)
        inner = sum(float(np.dot(x[k], z[k])) for k in x)
        assert inner == pytest.approx(plot.mass_inner([a, a], [b, b], contract), rel=3e-15)
        assert np.linalg.norm(tracker.shape_vector_from_components(x, shape_metric="full")) == pytest.approx(1.)


def tracking_state(vector, position, frequency, beta):
    t = plot.tracking_helpers()
    components = dict(u_left=np.array(vector, dtype=float), w_left=np.zeros(2),
                      u_right=np.zeros(2), w_right=np.zeros(2))
    return t.AnalyticModeState(epsilon=.005, beta=beta, mu=0., current_sorted_index=position,
        Lambda=frequency, coeff=np.empty(0), components=components,
        shape_vector=t.shape_vector_from_components(components, shape_metric="full"),
        smallest_singular_value=np.nan, singular_value_ratio=np.nan)


def test_continuation_keeps_identity_when_sorted_indices_swap_and_sign_flips():
    t = plot.tracking_helpers()
    states = [tracking_state([1., 0.], 1, 3., 0.), tracking_state([0., 1.], 2, 4., 0.)]
    previous = [t.make_branch_point(state=s, branch_id=f"branch_{i}", base_sorted_index=i,
        mac_to_previous=1., relative_lambda_jump=0., sorted_lambdas=[3., 4.], step_type="base",
        step_index=0, mac_warning_threshold=.8) for i, s in enumerate(states, 1)]
    candidates = [tracking_state([0., -1.], 1, 3.9, .1), tracking_state([-1., 0.], 2, 4.1, .1)]
    result, margins = plot.continue_class(previous, candidates, [3.9, 4.1], 1)
    assert [p.branch_id for p in result] == ["branch_1", "branch_2"]
    assert [p.current_sorted_index for p in result] == [2, 1]
    assert [p.Lambda for p in result] == [4.1, 3.9]
    assert margins == [1., 1.]
    assert [p.mac_to_previous for p in result] == [1., 1.]


def test_ambiguous_continuation_is_rejected():
    t = plot.tracking_helpers()
    s = tracking_state([1., 0.], 1, 3., 0.)
    previous = [t.make_branch_point(state=s, branch_id="branch_1", base_sorted_index=1,
        mac_to_previous=1., relative_lambda_jump=0., sorted_lambdas=[3.], step_type="base",
        step_index=0, mac_warning_threshold=.8)]
    candidates = [tracking_state([1., 1.], 1, 3., .1), tracking_state([1., -1.], 2, 3.1, .1)]
    with pytest.raises(ValueError, match="Ambiguous|depends"):
        plot.continue_class(previous, candidates, [3., 3.1], 1)


@pytest.mark.skipif(not (plot.MAPS / "_audit_reflection_classification_20261107/class_spectrum.csv").exists(),
                    reason="Optional private accepted archive is not distributed")
def test_saved_branch_regression_no_spectral_calls(monkeypatch):
    t = plot.tracking_helpers()
    def forbidden(*args, **kwargs):
        raise AssertionError("Spectral/basis evaluation is forbidden in this plot-only regression")
    for name in ("solve_sorted_lambdas", "solve_sorted_lambdas_with_warning", "solve_mode_states",
                 "mode_state_from_lambda", "find_first_n_roots", "assemble_clamped_coupled_matrix",
                 "reconstruct_analytic_components", "analytic_null_vector", "initialize_base_points",
                 "svd_refined_candidate", "augment_mode_states_with_svd_candidates", "track_path"):
        monkeypatch.setattr(t, name, forbidden)
    folder = plot.MAPS / "b775d0980d7bad77"
    _, verified = plot.sorted_spectrum(folder, "square_beta_mu0_tau0_refined")
    rows = plot.table(plot.MAPS / "_audit_reflection_classification_20261107/class_spectrum.csv")
    source = [r for r in plot.table(folder / "spectrum.csv") if r["map"] == "square_beta_mu0_tau0_refined"]
    contract = plot.read(plot.MAPS / "_audit_reflection_intersection_20261107/manifest.json")["physical_contract"]
    curves, checks = plot.archived_branches(rows, contract, verified["verified_sha256"], source)
    assert all(len(rs) == 197 for rs in curves.values())
    assert checks["minimum_mass_overlap"] > .9993
    assert checks["class_rank_changes"] == 0
    assert checks["minimum_assignment_margin"] > .96
    for branch, rs in curves.items():
        assert all(r["base_sorted_index"] == branch for r in rs)
        assert {r["symmetry_class"] for r in rs} == ({"+"} if branch % 2 else {"-"})
    pair = {b: {float(r["beta_deg"]): r for r in curves[b]} for b in (1, 2)}
    assert [pair[b][6.1]["original_sorted_position"] for b in (1, 2)] == ["1", "2"]
    assert [pair[b][6.2]["original_sorted_position"] for b in (1, 2)] == ["2", "1"]
    assert float(pair[1][6.1]["Lambda"]) < float(pair[2][6.1]["Lambda"])
    assert float(pair[1][6.2]["Lambda"]) > float(pair[2][6.2]["Lambda"])
    local = plot.first_classes(plot.table(plot.MAPS / "_audit_reflection_intersection_20261107/block_roots.csv"), local=True)
    solution = plot.read(plot.MAPS / "_audit_reflection_intersection_20261107/crossing_solution.json")
    augmented = plot.insert_branch_local_roots(curves, local, solution)
    assert len(augmented[1]) == len(augmented[2]) == 201
    for b in (1, 2):
        crossing = next(r for r in augmented[b] if float(r["beta_deg"]) == solution["beta_deg"])
        assert crossing["Lambda"] == solution["Lambda"]
        assert crossing["branch_id"] == curves[b][0]["branch_id"]
    assert all(len(augmented[b]) == 197 for b in range(3, 9))


def test_single_scale_preserves_relative_mode_amplitudes_and_input():
    modes = {"+": [arm(), arm()], "-": [arm(3), arm(3)]}
    original = modes["+"][0]["y"].copy()
    scenes, info = plot.presentation_scenes(modes, 6.13, 1., .28)
    scale = info["common_displacement_multiplier"]
    maxima = {c: max(np.linalg.norm(scale*d, axis=1).max() for _, d in scene)
              for c, scene in scenes.items()}
    assert maxima["-"]/maxima["+"] == pytest.approx(3.)
    assert maxima["-"] == pytest.approx(.28)
    np.testing.assert_array_equal(modes["+"][0]["y"], original)


def test_geometry_inward_coordinates_and_endpoint_continuity():
    beta = 30.
    a = arm()
    scenes = plot.coordinates([a, a], beta, 1.)
    np.testing.assert_allclose(scenes[0][0][0], [-1., 0.])
    np.testing.assert_allclose(scenes[1][0][0], [np.sqrt(3)/2, .5])
    for x, d in scenes:
        np.testing.assert_array_equal(x[-1], [0., 0.])
        np.testing.assert_array_equal(d[0], [0., 0.])
        assert len(x) == 4  # two saved Gauss nodes and two saved endpoints


def test_agreed_sign_choice_is_independent_of_archived_global_sign():
    modes = {"+": [arm(), arm()], "-": [arm(2), arm(2)]}
    flipped = {c: [{k: (-v if k in ("y", "dy", "ends") else v).copy()
                    for k, v in a.items()} for a in arms] for c, arms in modes.items()}
    before, b = plot.presentation_scenes(modes, 6.13, 1., .28)
    after, a = plot.presentation_scenes(flipped, 6.13, 1., .28)
    assert a["common_displacement_multiplier"] == b["common_displacement_multiplier"]
    for c in modes:
        for (_, d), (_, e) in zip(before[c], after[c]):
            np.testing.assert_array_equal(d, e)


def guarded_point():
    roots = [dict(sorted_position=i, Lambda=float(i), omega=i*i/200, Omega=float(i*i),
                  role="plot" if i <= 8 else "guard") for i in range(1, 10)]
    record = dict(contract=dict(
        model=dict(epsilon=.005, nu=.3, aspect=1., E=1., rho=1., l=1.,
                   ks=5/6, variant="RL_FULL", mu=0., tau=0.), K_plot=8, K_guard=9,
        case=dict(beta=30., mu=0., tau=0., k_joint="RIGID", k_end="RIGID"),
        spectrum_semantics="sorted_positions", coordinates="outer clamp to joint",
        state_order=["u", "w", "psi", "N", "Q", "M"], time_reference=200.),
        beta_deg=30., mu=0., tau=0., status="PASS", roots=roots,
        search=dict(status="PASS", lower_count=0, requested_guard=9, computed_ranks=9,
                    guard_right_count=9, unresolved_intervals=[], guard_right_omega=.5),
        quality=dict(status="PASS", gates=[dict(boundary=True, mass=True) for _ in roots]))
    rows = [dict(rank=str(r["sorted_position"]), Lambda=str(r["Lambda"]), omega=str(r["omega"]),
                 Omega=str(r["Omega"]), status="PASS") for r in roots]
    return record, rows


@pytest.mark.parametrize("fault", ["missing_guard", "bad_count", "wrong_model", "changed_csv", "failed_quality"])
def test_sorted_reader_refuses_incomplete_or_incompatible_saved_spectrum(fault):
    record, rows = guarded_point()
    plot.validate_guard_record(record, rows, 30.)
    if fault == "missing_guard":
        record["roots"].pop()
    elif fault == "bad_count":
        record["search"]["guard_right_count"] = 8
    elif fault == "wrong_model":
        record["contract"]["model"]["variant"] = "RL_PLANAR_LEGACY"
    elif fault == "changed_csv":
        rows[0]["Lambda"] = "1.001"
    else:
        record["quality"]["gates"][3]["mass"] = False
    with pytest.raises(ValueError):
        plot.validate_guard_record(record, rows, 30.)


def test_saved_crossing_adds_only_two_known_roots_and_keeps_original_values():
    curves = {k: [dict(beta_deg="6.1", Lambda=str(k+2.1)), dict(beta_deg="6.2", Lambda=str(k+2.2))]
              for k in range(1, 9)}
    solution = dict(first_class_roots_certified=True, beta_deg=6.13, Lambda=3.924, omega=.077)
    result = plot.insert_saved_crossing(curves, solution)
    for k in range(1, 9):
        original = [r for r in result[k] if r["plot_source"] == "spectrum.csv"]
        assert [(r["beta_deg"], r["Lambda"]) for r in original] == [(r["beta_deg"], r["Lambda"]) for r in curves[k]]
        assert len(result[k]) == (3 if k <= 2 else 2)
    assert result[1][1]["Lambda"] == result[2][1]["Lambda"] == solution["Lambda"]
    solution["first_class_roots_certified"] = False
    with pytest.raises(ValueError):
        plot.insert_saved_crossing(curves, solution)


def test_branch_descriptions_require_observed_frequency_behaviour():
    curves = {"+": [dict(Lambda=2.), dict(Lambda=4.)],
              "-": [dict(Lambda=3.9242), dict(Lambda=3.9236)]}
    plot.branch_behaviour(curves)
    curves["+"] = list(reversed(curves["+"]))
    with pytest.raises(ValueError):
        plot.branch_behaviour(curves)
