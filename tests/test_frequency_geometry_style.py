"""Geometry and saved-data presentation only; no spectra or tracking runs."""
import csv
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from scripts.lib import frequency_geometry_style as g


@pytest.fixture(autouse=True)
def no_spectral_work(monkeypatch):
    from scripts.analysis import map_rl_full_rr as maps
    from scripts.lib import rl_full_rr_sections_prefix as prefix
    def forbidden(*args, **kwargs):
        raise AssertionError("Spectral call forbidden in geometry/presentation tests")
    for obj, name in [(maps, "solve"), (prefix, "counted_prefix"), (prefix, "catalogues"),
                      (prefix.p, "ritz"), (prefix.p, "boundary"), (prefix.p, "count_at"),
                      (prefix.p, "mode"), (prefix.p.base, "basis"),
                      (prefix.p.base, "state_matrix"), (prefix.p.base, "joint_matrix")]:
        monkeypatch.setattr(obj, name, forbidden)
    yield
    plt.close("all")


@pytest.mark.parametrize("r1,r2,expected", [
    (11, 12, "solid"), (9, 12, "dashed"), (12, 9, "dashed"),
    (9, 8, "dashed"), (10, 12, "solid"), (12, 10, "solid"), (10, 10, "solid"),
    (10-1e-12, 12, "dashed")])
def test_individual_arms_and_inclusive_equality(r1, r2, expected):
    node = g.classify(r1, r2)
    assert node["line_style"] == expected
    assert node["slenderness_criterion_pass"] == (expected == "solid")
    assert node["chi"] == min(r1, r2)
    assert node["slenderness_threshold"] == 10


@pytest.mark.parametrize("mu,epsilon,style,failed_arm", [
    (0., .005, "solid", None), (.9, .005, "dashed", 1),
    (-.9, .005, "dashed", 2), (0., .05, "dashed", "both")])
def test_canonical_geometries(mu, epsilon, style, failed_arm):
    node = g.geometry_node(mu=mu, epsilon=epsilon)
    assert node["line_style"] == style
    failures = [i for i in (1, 2) if node[f"L{i}_over_h{i}"] < 10]
    assert failures == ([] if failed_arm is None else [1, 2] if failed_arm == "both" else [failed_arm])
    # Uniform length scaling must not change classification.
    scaled = g.geometry_node(mu=mu, epsilon=epsilon, l=17)
    assert scaled["chi"] == pytest.approx(node["chi"])


def test_canonical_exact_equality():
    node = g.geometry_node(epsilon=1/(10*math.sqrt(12)))
    assert node["chi"] == 10
    assert node["line_style"] == "solid"


@pytest.mark.parametrize("epsilon,tau,expected", [
    (.005, -.5, "dashed"), (.005, 0, "dashed"), (.005, .5, "solid"),
    (.01, -.5, "dashed"), (.01, 0, "dashed"), (.01, .5, "dashed")])
def test_plan_geometries_and_reflection(epsilon, tau, expected):
    node = g.geometry_node(mu=.9, tau=tau, epsilon=epsilon, beta_deg=15)
    mirror = g.geometry_node(mu=-.9, tau=-tau, epsilon=epsilon, beta_deg=90)
    assert node["line_style"] == mirror["line_style"] == expected
    assert node["L1_over_h1"] == pytest.approx(mirror["L2_over_h2"])


@pytest.mark.parametrize("reverse", [False, True])
def test_geometric_boundary_between_nodes_no_new_observation(reverse):
    x, y = [.8, .9], [2., 4.]
    if reverse:
        x.reverse(); y.reverse()
    original = (x.copy(), y.copy())
    geometry = g.FrequencyGeometry(x, axis="MU", fixed={"tau": 0}, material={"epsilon": .005})
    boundary = 1-10*math.sqrt(12)*.005
    assert geometry.boundaries == pytest.approx([boundary], abs=2e-15)
    assert len(geometry.nodes) == 2
    fig, ax = plt.subplots()
    lines = g.plot_frequency_curve(ax, geometry, y, color="#123456", marker="o", label="B1")
    assert [a.get_linestyle() for a in lines] == (["--", "-"] if reverse else ["-", "--"])
    assert {a.get_color() for a in lines} == {"#123456"}
    expected_y = 2+(boundary-.8)*20
    for line in lines:
        vertices = list(zip(line.get_xdata(), line.get_ydata()))
        assert any(v == pytest.approx((boundary, expected_y)) for v in vertices)
        assert all(vertices[i][0] in x for i in line.get_markevery())
    assert original == (x, y)
    assert g.add_geometry_caption(fig) == g.CAPTION


def test_boundary_at_saved_node_belongs_to_solid():
    epsilon = 1/(20*math.sqrt(12))
    geometry = g.FrequencyGeometry([0., .5, .75], axis="mu", fixed={}, material={"epsilon": epsilon})
    assert geometry.nodes[1]["chi"] == 10
    assert geometry.nodes[1]["line_style"] == "solid"
    assert geometry.boundaries == [.5]
    assert len(geometry.nodes) == 3
    _, ax = plt.subplots()
    lines = g.plot_frequency_curve(ax, geometry, [1, 2, 3], color="black", marker="o")
    marked = [(line.get_linestyle(), line.get_xdata()[i]) for line in lines for i in line.get_markevery()]
    assert ("-", .5) in marked and ("--", .5) not in marked


def test_tau_two_transitions_even_when_both_endpoints_fail():
    geometry = g.FrequencyGeometry([-.8, .8], axis="tau", fixed={"mu": 0}, material={"epsilon": .025})
    assert all(n["line_style"] == "dashed" for n in geometry.nodes)
    assert len(geometry.boundaries) == 2
    assert [s["line_style"] for s in geometry.segments] == ["dashed", "solid", "dashed"]
    for value in geometry.boundaries:
        assert geometry.at(value)["chi"] == pytest.approx(10, abs=1e-12)
    reverse = g.FrequencyGeometry([.8, -.8], axis="tau", fixed={"mu": 0}, material={"epsilon": .025})
    assert reverse.boundaries == geometry.boundaries


def test_epsilon_sweep_and_beta_independence():
    geometry = g.FrequencyGeometry([.005, .01], axis="EPSILON", fixed={"mu": .8}, material={})
    assert geometry.boundaries == pytest.approx([.2/(10*math.sqrt(12))])
    assert [n["line_style"] for n in geometry.nodes] == ["solid", "dashed"]
    angle = g.FrequencyGeometry([90, 0], axis="BETA", fixed={"mu": .9, "tau": .5}, material={"epsilon": .005})
    assert angle.boundaries == []
    assert angle.nodes[0]["chi"] == angle.nodes[1]["chi"]


def test_all_branches_share_classification_and_caption_is_conditional():
    geometry = g.FrequencyGeometry([.8, .9], axis="mu", fixed={}, material={"epsilon": .005})
    fig, ax = plt.subplots()
    for rank in range(10):
        lines = g.plot_frequency_curve(ax, geometry, [rank+1, rank+2], color=f"C{rank}")
        assert [line.get_linestyle() for line in lines] == ["-", "--"]
        assert len({line.get_color() for line in lines}) == 1
    assert g.add_geometry_caption(fig)
    fig2, axes = plt.subplots(1, 2)
    for ax, mu in zip(axes, (0, .9)):
        plan = g.FrequencyGeometry([0, 90], axis="beta", fixed={"mu": mu}, material={"epsilon": .005})
        g.plot_frequency_curve(ax, plan, [1, 2], color="black")
    assert g.add_geometry_caption(fig2)  # solid panel + dashed panel
    fig3, ax = plt.subplots()
    g.plot_frequency_curve(ax, plan, [1, 2], color="black")
    assert not g.add_geometry_caption(fig3)  # entirely dashed
    ax.set_xlim(0, 90)
    solid = g.FrequencyGeometry([0, 90], axis="beta", fixed={}, material={})
    fig4, ax = plt.subplots()
    g.plot_frequency_curve(ax, solid, [1, 2], color="black")
    assert not g.add_geometry_caption(fig4)


def test_diagnostics_preserve_nodes_and_partial_style(tmp_path):
    geometry = g.FrequencyGeometry([.8, .9], axis="mu", fixed={}, material={})
    report = g.write_geometry_diagnostics(tmp_path, {"example": geometry})
    rows = list(csv.DictReader((tmp_path/"geometry_nodes.csv").open(encoding="utf-8")))
    assert [float(r["x"]) for r in rows] == [.8, .9]
    assert report["series"]["example"]["nodes"] == geometry.nodes
    assert not report["series"]["example"]["boundary_vertices_are_spectral_samples"]
    assert json.loads((tmp_path/"geometry_diagnostics.json").read_text(encoding="utf-8")) == report
    with pytest.raises(FileExistsError):
        g.write_geometry_diagnostics(tmp_path, {"example": geometry})
    _, ax = plt.subplots()
    partial, = g.plot_partial_frequency(ax, [.8, .9], [1, 2], label="CC, arm 1")
    assert partial.get_linestyle() == ":"


def test_nan_gap_is_never_filled_by_geometric_split():
    geometry = g.FrequencyGeometry([.8, .85, .9], axis="mu", fixed={}, material={})
    for _, _, y in geometry.drawing_runs([1, np.nan, 3]):
        assert not any(math.isfinite(a) and math.isfinite(b) for a, b in zip(y, y[1:]))


@pytest.mark.parametrize("mode,grid,fixed", [
    ("MU", [.8, .9], {"beta_deg": 15, "tau": 0}),
    ("TAU", [-.5, .5], {"beta_deg": 90, "mu": .9}),
    ("EPSILON", [.005, .01], {"beta_deg": 15, "mu": .8, "tau": 0})])
def test_common_renderer_uses_geometry_and_writes_manifest(tmp_path, monkeypatch, mode, grid, fixed):
    from scripts.analysis import map_rl_full_rr as maps
    folder = tmp_path/"source"; folder.mkdir()
    rows = [dict(map="fixture", rank=k, x=x, Lambda=k+x, status="PASS") for x in grid for k in range(1, 10)]
    with (folder/"spectrum.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
    maps.publish(folder/"manifest.json", dict(artifacts={"spectrum.csv": maps.sha(folder/"spectrum.csv")},
        identity=dict(material=dict(epsilon=.005, aspect=1.),
                      maps=[dict(name="fixture", mode=mode, fixed=fixed, grid=grid)])))
    before = {p.name: maps.sha(p) for p in folder.iterdir()}
    monkeypatch.setattr(maps, "OUTPUT", tmp_path/"derived")
    output = maps.render(folder, dict(layout="split", dpi=35))
    report = maps.read(output/"manifest.json")
    nodes = report["geometry_diagnostics"]["series"]["fixture"]["nodes"]
    assert [n["x"] for n in nodes] == grid
    assert {n["line_style"] for n in nodes} == {"solid", "dashed"}
    assert report["geometry_captions"]["fixture_split"]
    assert before == {p.name: maps.sha(p) for p in folder.iterdir()}


def test_publication_draw_replaces_style_for_all_branches():
    from scripts.analysis import plot_rl_full_rr_publication as pub
    data = dict(spec=dict(name="fixture", mode="BETA", fixed=dict(mu=.9, tau=0)),
                material=dict(epsilon=.005, aspect=1.), x=[0, 90], y=[[k for k in range(1,10)]]*2,
                window=[0,90])
    fig, ax = plt.subplots()
    report = pub.draw(ax, data, range(1,9))
    assert len(ax.lines) == 8
    assert all(line.get_linestyle() == "--" for line in ax.lines)
    assert all(node["line_style"] == "dashed" for node in report["geometry"]["nodes"])


def test_reflection_presenter_uses_shared_rule_without_field_processing(tmp_path):
    from scripts.analysis import plot_rl_full_rr_reflection_publication as pub
    rows = [{"beta_deg": x, "Lambda": 3+x/100} for x in (0., 90.)]
    curves = {"+": rows, "-": rows}
    axis = np.array([[0., 0.], [1., 0.]])
    scenes = {s: [(axis, np.zeros_like(axis))] for s in ("+", "-")}
    report = pub.draw(tmp_path, curves, curves, dict(beta_deg=6.13, Lambda=3.924),
        scenes, dict(common_displacement_multiplier=1), 35, 9,
        contract=dict(epsilon=.005, aspect=1., mu=.9, tau=0.))
    assert all(n["line_style"] == "dashed" for s in report["geometry_diagnostics"]["series"].values() for n in s["nodes"])
    assert report["geometry_captions"]["class_frequencies"] == ""


def test_continuation_presenter_reuses_saved_contract_without_assignment(tmp_path, monkeypatch):
    from scripts.analysis import postprocess_rl_full_rr_continuation as post
    def forbidden(*args, **kwargs):
        raise AssertionError("Tracking is not part of geometric presentation")
    monkeypatch.setattr(post, "continue_series", forbidden)
    specs = [dict(id=f"s{i}", axis="mu") for i in range(4)]
    points = [dict(beta=15., mu=mu, json=dict(contract=dict(model=dict(epsilon=.005, aspect=1., mu=.8, tau=0.)))) for mu in (.8, .9)]
    rows = [dict(series=s["id"], branch_id=f"B{k}", x=p["mu"], Lambda=k+p["mu"], status="CONFIRMED_FROM_SEED")
            for s in specs for k in (1, 2) for p in points]
    post.draw(tmp_path, dict(series=specs, dpi=35), rows, [], {s["id"]: points for s in specs})
    report = json.loads((tmp_path/"geometry_presentation.json").read_text(encoding="utf-8"))
    assert all(report["geometry_captions"].values())
    assert [n["line_style"] for n in report["geometry_diagnostics"]["series"]["s0"]["nodes"]] == ["solid", "dashed"]
