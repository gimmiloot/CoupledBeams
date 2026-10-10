"""Bounded saved-profile orchestration contracts; no scientific solves."""
from __future__ import annotations

import copy
import inspect
from pathlib import Path

import numpy as np
import pytest

from scripts.lib import nlsp_spatial_profile_audit as audit

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/nlsp_nonlinear_dynamic_validation/c6256269eb8143ef"


def _config():
    return copy.deepcopy(audit.read(audit.CONFIG))


def _no_science(monkeypatch):
    from scripts.lib import nlsp_fem3a_1d_reference as one
    from scripts.analysis import resume_nlsp_nonlinear_dynamic_3d_fem as resume

    def forbidden(*args, **kwargs):
        pytest.fail("A saved-profile test attempted a physical or scientific solve")

    for target, name in ((audit.fem.fem1, "run_job"),
                         (audit.fem.fem1.single, "generate_mesh_with_gmsh_cli"),
                         (audit.fem.fem1.single, "generate_mesh_with_gmsh_python"),
                         (audit.fem.fem2, "fem2_static_newton"),
                         (one.runner, "integrate_case"),
                         (one.dynamics.PlanarGalerkin, "linear_eigenpairs"),
                         (one.rod, "derive_polynomials"),
                         (resume.base, "run_case")):
        monkeypatch.setattr(target, name, forbidden)
    return forbidden


def _cache(tmp_path, monkeypatch):
    root = tmp_path / "fixture_root"
    root.mkdir()
    monkeypatch.setattr(audit, "ROOT", root)
    parent = root / "parent"
    parent.mkdir()
    audit.write(parent / "manifest.json", {"artifact_hashes": {}})
    (parent / "selected_frame.bin").write_bytes(b"immutable saved nodal field fixture")
    (parent / "dense.bin").write_bytes(b"immutable accepted dense polynomial fixture")
    (root / "frozen_model.py").write_bytes(b"frozen physics fixture")
    identity = {"config": {"source_bundles": {"FEM3C": {
        "path": "parent", "manifest_sha256": audit.sha(parent / "manifest.json")}}},
        "source_artifacts": {"parent/selected_frame.bin": audit.sha(parent / "selected_frame.bin")},
        "additional_postprocessing_sources": {"parent/dense.bin": audit.sha(parent / "dense.bin")},
        "frozen_helpers": {"frozen_model.py": audit.sha(root / "frozen_model.py")}}
    bundle = root / "diagnostic"
    bundle.mkdir()
    audit.write(bundle / "provenance.json", identity)
    audit.write(bundle / "config.json", identity["config"])
    np.savez(bundle / "derived.npz", values=np.arange(4.))
    summary = {"status": "DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS", "completed": True,
        "scientific_calls": 0, "states": {}, "one_d_c_spatial_status": "PARTIAL"}
    audit.save(bundle, identity, summary)
    return root, bundle, identity, summary


def _time_fixture(tmp_path, monkeypatch):
    root = tmp_path / "times_root"
    root.mkdir()
    monkeypatch.setattr(audit, "ROOT", root)
    source = root / "saved"
    provenance = {"config": {"omega1": 2 * np.pi}}  # T1=1 for this fixture.
    for stage, times in (("full_period_medium", [.001, .249, .501, .749, 1.]),
                         ("medium_refined_time", [.001, .124, .25]),
                         ("fine_refined_time", [.001, .124, .25])):
        path = source / "cases" / stage / "nonlinear"
        path.mkdir(parents=True)
        np.savez(path / "section_history.npz", time=np.array(times),
            increments=np.arange(1, len(times)+1) * 2)
    return source, provenance


def test_default_policy_is_postprocessing_only_and_frozen():
    config = _config()
    assert audit.validate_config(config) is config
    assert config["new_physical_solves"] is False
    assert config["new_meshes"] is False
    assert config["smoothing"] is False
    assert config["phase_amplitude_fitting"] is False
    assert config["historical_statuses_unchanged"] is True
    assert config["primary_section_count"] == 41
    assert config["maximum_main_figures"] == 4
    assert config["source_bundles"] == audit.SOURCE_BUNDLES


@pytest.mark.parametrize("key,value", [
    ("section_counts", [41, 81, 161]), ("primary_section_count", 81),
    ("snapshot_T1_fractions", [0., .25, .5, .75, 2.]),
    ("spatial_comparison_T1_fractions", [0., .25, .5]),
    ("maximum_main_figures", 5), ("one_d_boundary_partition_ell_c", 2.),
    ("FEM_spatial_grid_points", 1601), ("one_d_spatial_grid_points", 1001),
    ("smoothing", True), ("phase_amplitude_fitting", True),
    ("new_physical_solves", True), ("new_meshes", True),
    ("historical_statuses_unchanged", False), ("authorization", "old_FEM3C_permission"),
    ("native_time_selection", "interpolate native nodal histories"),
    ("original_figure_time_policy", "replace old plots"),
    ("independent_strain_quadrature", "single centroid"),
    ("quadratic_fit_probe", "replace historical recovery"),
])
def test_changed_authorization_or_numerical_scope_is_rejected(key, value):
    config = _config()
    config[key] = value
    with pytest.raises(ValueError, match="bounded"):
        audit.validate_config(config)


@pytest.mark.parametrize("change", ["bundle", "manifest", "extra_source"])
def test_parent_source_identity_cannot_silently_change(change):
    config = _config()
    if change == "bundle":
        config["source_bundles"]["FEM3C"]["path"] = "other_geometry"
    elif change == "manifest":
        config["source_bundles"]["FEM3C"]["manifest_sha256"] = "0" * 64
    else:
        config["source_bundles"]["FEM3D"] = {"path": "unrequested", "manifest_sha256": "0" * 64}
    with pytest.raises(ValueError, match="bounded"):
        audit.validate_config(config)


def test_explicit_historical_manifest_hashes_remain_unchanged():
    for evidence in audit.SOURCE_BUNDLES.values():
        assert audit.sha(ROOT / evidence["path"] / "manifest.json") == evidence["manifest_sha256"]


def test_checked_artifact_records_hash_and_rejects_tampering(tmp_path, monkeypatch):
    monkeypatch.setattr(audit, "ROOT", tmp_path)
    source = tmp_path / "source"
    source.mkdir()
    path = source / "native.bin"
    path.write_bytes(b"original nodal displacement evidence")
    digest = audit.sha(path)
    manifest = {"artifact_hashes": {"native.bin": digest}}
    registry = {}
    assert audit.checked(path, manifest, source, registry) == path
    assert registry == {"source/native.bin": digest}
    path.write_bytes(b"tampered")
    with pytest.raises(ValueError, match="corrupt"):
        audit.checked(path, manifest, source, registry)


def test_checked_artifact_cannot_escape_source_bundle(tmp_path, monkeypatch):
    monkeypatch.setattr(audit, "ROOT", tmp_path)
    source = tmp_path / "source"; source.mkdir()
    other = tmp_path / "elsewhere.bin"; other.write_bytes(b"wrong source")
    with pytest.raises(ValueError):
        audit.checked(other, {"artifact_hashes": {}}, source, {})


def test_saved_time_selection_preserves_actual_coordinates_and_static_datum(tmp_path, monkeypatch):
    source, provenance = _time_fixture(tmp_path, monkeypatch)
    period, states = audit.selected_states(source, provenance)
    assert period == 1.
    assert len(states) == 11
    full_quarter = next(s for s in states if s["stage"] == "full_period_medium" and s["requested_tau"] == .25)
    assert full_quarter["actual_time"] == .249
    assert full_quarter["requested_time"] == .25
    assert full_quarter["time_offset"] == pytest.approx(-.001)
    assert full_quarter["exact_requested_time_available"] is False
    assert full_quarter["time_policy"] == "nearest_actual_native_frame_no_time_interpolation"
    static = [s for s in states if s["requested_tau"] == 0]
    assert all(s["step"] == 1 and s["actual_time"] == 0 and s["history_index"] is None for s in static)
    assert all(s["time_policy"] == "actual_STATIC_preload" for s in static)
    fine = [s for s in states if s["mesh"] == "fine"]
    assert [s["requested_tau"] for s in fine] == [0., .125, .25]
    assert not any(s["actual_time"] > .25 for s in fine)


@pytest.mark.parametrize("bad_times,bad_increments", [
    ([], []), ([.001, np.nan, .25], [2, 4, 6]),
    ([.001, .124, .124], [2, 4, 6]), ([.001, .124, .24], [2, 4, 6]),
    ([.001, .124, .25], [2, 4.5, 6]), ([.001, .124, .25], [2, 4, 4]),
    ([.001, .124, .25], [2, 4]),
])
def test_invalid_or_incomplete_saved_time_coverage_is_not_fabricated(tmp_path, monkeypatch, bad_times, bad_increments):
    source, provenance = _time_fixture(tmp_path, monkeypatch)
    np.savez(source / "cases/fine_refined_time/nonlinear/section_history.npz",
        time=np.asarray(bad_times), increments=np.asarray(bad_increments))
    with pytest.raises(ValueError, match="saved-time coverage"):
        audit.selected_states(source, provenance)


def test_medium_fine_time_pairing_requires_actual_equality(tmp_path, monkeypatch):
    source, provenance = _time_fixture(tmp_path, monkeypatch)
    np.savez(source / "cases/fine_refined_time/nonlinear/section_history.npz",
        time=np.array([.001, .125, .25]), increments=np.array([2, 4, 6]))
    with pytest.raises(ValueError, match="time pairing"):
        audit.selected_states(source, provenance)


def test_actual_saved_fine_coverage_stops_at_quarter_period():
    provenance = audit.read(SOURCE / "provenance.json")
    T, states = audit.selected_states(SOURCE, provenance)
    assert len(states) == 11
    assert sum(s["mesh"] == "fine" for s in states) == 3
    assert not any(s["mesh"] == "fine" and s["actual_time"] > .25*T for s in states)
    for stage in ("medium_refined_time", "fine_refined_time"):
        with np.load(SOURCE / "cases" / stage / "nonlinear/section_history.npz") as saved:
            assert saved["time"][-1] == .25*T
    for state in states:
        assert (ROOT / state["frame"]).is_file()


def test_reordered_saved_nodes_stop_before_recovery(tmp_path, monkeypatch):
    _no_science(monkeypatch)
    root, bundle = tmp_path / "root", tmp_path / "bundle"
    root.mkdir(); bundle.mkdir()
    monkeypatch.setattr(audit, "ROOT", root)
    code = root / "scripts/lib"
    code.mkdir(parents=True)
    (code / "nlsp_profile_1d_diagnostics.py").write_bytes(b"pure helper fixture")
    np.savez(root / "field.npz", node_ids=np.array([2, 1]), U=np.zeros((2, 3)))
    monkeypatch.setattr(audit.fem.fem1.single, "read_gmsh_inp_mesh_data", lambda path: {})
    monkeypatch.setattr(audit.fem, "prepare_mesh", lambda mesh, rho: {"node_ids": np.array([1, 2])})
    monkeypatch.setattr(audit.fem, "analyze_state", lambda *args, **kwargs: pytest.fail("Reordered nodes reached recovery"))
    item = {"source_meshes": {"medium": {"include": "mesh.inp"}},
        "selected_states": [{"name": "one", "frame": "field.npz", "mesh": "medium"}]}
    with pytest.raises(ValueError, match="field ordering"):
        audit.compute(bundle, item, {"states": {}})


@pytest.mark.parametrize("registry,name", [
    ("source_artifacts", "parent/selected_frame.bin"),
    ("additional_postprocessing_sources", "parent/dense.bin"),
    ("frozen_helpers", "frozen_model.py"),
])
def test_cached_report_rechecks_sources_even_if_parent_manifest_is_unchanged(tmp_path, monkeypatch, registry, name):
    root, bundle, identity, summary = _cache(tmp_path, monkeypatch)
    parent_sha = audit.sha(root / "parent/manifest.json")
    assert audit.validate_cache(bundle) == summary
    (root / name).write_bytes(b"tampered after parent hashing")
    assert audit.sha(root / "parent/manifest.json") == parent_sha
    with pytest.raises(ValueError, match="source artifact changed"):
        audit.validate_cache(bundle)


def test_cache_rejects_changed_historical_manifest(tmp_path, monkeypatch):
    root, bundle, identity, summary = _cache(tmp_path, monkeypatch)
    (root / "parent/manifest.json").write_bytes(b"changed")
    with pytest.raises(ValueError, match="Historical manifest changed"):
        audit.validate_cache(bundle)


def test_cache_rejects_derived_array_tampering(tmp_path, monkeypatch):
    root, bundle, identity, summary = _cache(tmp_path, monkeypatch)
    np.savez(bundle / "derived.npz", values=np.zeros(4))
    with pytest.raises(ValueError, match="cache artifact mismatch"):
        audit.validate_cache(bundle)


@pytest.mark.parametrize("mode", ["compute", "report-only", "plot-only"])
def test_matching_cached_modes_have_zero_scientific_calls(tmp_path, monkeypatch, mode):
    root, bundle, identity, summary = _cache(tmp_path, monkeypatch)
    forbidden = _no_science(monkeypatch)
    monkeypatch.setattr(audit.fem, "analyze_state", forbidden)
    monkeypatch.setattr(audit.fem, "prepare_mesh", forbidden)
    monkeypatch.setattr(audit, "prepare", lambda config: (bundle, identity, audit.validate_cache(bundle)))
    if hasattr(audit, "complete_saved_evidence"):
        # The production callback rechecks saved evidence/code phases. This
        # synthetic fixture replaces only that read-only check, while every
        # scientific execution route remains forbidden above.
        monkeypatch.setattr(audit, "complete_saved_evidence", lambda bundle, item, summary: summary)
    figures = []

    def mock_plot(path):
        assert Path(path) == bundle
        figures.extend(f"figure_{i}" for i in range(4))
        return figures

    monkeypatch.setattr(audit, "plot", mock_plot)
    before_source = audit.sha(root / "parent/selected_frame.bin")
    args = ["--compute"] if mode == "compute" else ["--"+mode, str(bundle)]
    result = audit.main(args)
    assert result["scientific_calls"] == 0 and result["completed"] is True
    assert len(figures) == (0 if mode == "report-only" else 4)
    assert audit.sha(root / "parent/selected_frame.bin") == before_source


def test_completed_compute_returns_without_recovery_or_new_states(monkeypatch):
    forbidden = _no_science(monkeypatch)
    monkeypatch.setattr(audit.fem, "prepare_mesh", forbidden)
    monkeypatch.setattr(audit.fem, "analyze_state", forbidden)
    monkeypatch.setattr(audit.fem.fem1.single, "read_gmsh_inp_mesh_data", forbidden)
    summary = {"completed": True, "status": "DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS"}
    assert audit.compute(Path("unused_cache"), {}, summary) is summary


@pytest.mark.parametrize("flag,target", [
    ("--profile-audit", "scripts.lib.nlsp_spatial_profile_audit"),
    ("--validation", "scripts.lib.nlsp_fem3c_validation"),
    ("--long-horizon", "scripts.lib.nlsp_fem3b_continuation"),
])
def test_existing_dispatcher_routes_old_and_new_flags_without_scientific_calls(monkeypatch, flag, target):
    import importlib
    from scripts.analysis import resume_nlsp_nonlinear_dynamic_3d_fem as dispatcher
    _no_science(monkeypatch)
    calls = []
    module = importlib.import_module(target)
    sentinel = {"saved": True}
    monkeypatch.setattr(module, "main", lambda arguments: calls.append(arguments) or sentinel)
    assert dispatcher.main([flag, "--report-only", "saved_bundle"]) is sentinel
    assert calls == [["--report-only", "saved_bundle"]]


def test_existing_default_dispatcher_cache_path_is_preserved(monkeypatch):
    from scripts.analysis import resume_nlsp_nonlinear_dynamic_3d_fem as dispatcher
    _no_science(monkeypatch)
    summary = {"resume_statuses": {"historical": "PARTIAL"}, "overall": "PILOT_COMPLETE_WITH_QUALIFICATIONS"}
    monkeypatch.setattr(dispatcher, "validate_cache", lambda bundle: summary)
    assert dispatcher.main(["--report-only", "old_bundle"]) is summary


def test_manifest_contains_explicit_source_and_derived_hashes(tmp_path, monkeypatch):
    root, bundle, identity, summary = _cache(tmp_path, monkeypatch)
    manifest = audit.read(bundle / "manifest.json")
    assert manifest["scientific_calls"] == 0
    assert manifest["source_manifest_hashes"] == identity["config"]["source_bundles"]
    for name in ("provenance.json", "config.json", "summary.json", "derived.npz"):
        assert manifest["artifact_hashes"][name] == audit.sha(bundle / name)


def test_json_rejects_nonfinite_values_instead_of_silent_nan(tmp_path):
    with pytest.raises(ValueError):
        audit.write(tmp_path / "bad.json", {"field": np.nan})


def test_orchestration_source_does_not_invoke_physical_solvers():
    source = inspect.getsource(audit)
    for forbidden in ("run_job(", "run_case(", "integrate_case(", "solve_ivp(",
                      "fem2_static_newton(", "derive_polynomials(", "linear_eigenpairs(",
                      "generate_mesh_with_gmsh"):
        assert forbidden not in source


def test_checked_accepts_legacy_serialized_action_and_current_manifest_keys(tmp_path, monkeypatch):
    monkeypatch.setattr(audit, "ROOT", tmp_path)
    source = tmp_path / "serialized_action"
    source.mkdir()
    artifact = source / "result.json"
    artifact.write_bytes(b'{"polynomials":"saved immutable action fixture"}')
    digest = audit.sha(artifact)
    for key in ("artifacts", "artifact_hashes"):
        registry = {}
        assert audit.checked(artifact, {key: {"result.json": digest}}, source, registry) == artifact
        assert registry == {"serialized_action/result.json": digest}
        with pytest.raises(ValueError, match="corrupt"):
            audit.checked(artifact, {key: {"result.json": "0" * 64}}, source, {})


def test_interrupted_compute_preserves_original_code_phase_before_copy(tmp_path, monkeypatch):
    _no_science(monkeypatch)
    root = tmp_path / "fixture_root"
    root.mkdir()
    monkeypatch.setattr(audit, "ROOT", root)
    bundle = tmp_path / "interrupted_audit"
    phase = bundle / "execution_code"
    phase.mkdir(parents=True)
    snapshot = phase / Path(audit.__file__).name
    original = b"original executed postprocessing source must remain immutable"
    snapshot.write_bytes(original)
    copies = []

    def unexpected_copy(*args, **kwargs):
        copies.append(args)
        pytest.fail("Changed interrupted code phase was copied before preservation gate")

    monkeypatch.setattr(audit.shutil, "copyfile", unexpected_copy)
    with pytest.raises(ValueError, match="Interrupted diagnostic code phase changed"):
        audit.compute(bundle, {"source_meshes": {}, "selected_states": []}, {"states": {}})
    assert copies == []
    assert snapshot.read_bytes() == original
    assert list(phase.iterdir()) == [snapshot]


def test_saved_historical_time_policy_is_distinct_from_native_samples(monkeypatch):
    _no_science(monkeypatch)
    bundle = ROOT / "results/nlsp_spatial_profile_audit/c19f5a82203c0260"
    if not (bundle / "historical_time_policy.json").exists():
        pytest.skip("Saved bounded profile audit is not present in this checkout")
    evidence = audit.read(bundle / "historical_time_policy.json")
    with np.load(bundle / "historical_time_policy.npz", allow_pickle=False) as native, np.load(
            bundle / "historical_figure_profiles.npz", allow_pickle=False) as original:
        assert len(evidence["rows"]) == 5
        for i, row in enumerate(evidence["rows"]):
            fraction = row["right_fraction"]
            assert 0 <= fraction <= 1
            replay = ((1-fraction)*native["left_native_fields"][i]
                + fraction*native["right_native_fields"][i])
            assert np.max(abs(replay-original["three_d_nonlinear_fields"][i])) < 1e-13
            difference = np.max(abs(original["three_d_nonlinear_fields"][i]
                - native["nearest_actual_fields"][i]), axis=0)
            np.testing.assert_allclose(difference, row["historical_minus_nearest_actual_field_max"], rtol=0, atol=1e-20)
            if row["requested_tau"] in (0., 1.):
                assert row["requested_time"] == row["nearest_actual_time"]
                assert difference[6] == 0
                assert np.max(difference) < 1e-20  # tiny inactive-field interpolation roundoff
            else:
                assert row["requested_time"] != row["nearest_actual_time"]
                assert difference[6] > 0
