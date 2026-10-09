"""FEM-1R bounded continuation checks; no meshing, solver jobs, or 1D solves."""
import copy
import math
from pathlib import Path

import numpy as np
import pytest

from scripts.analysis import verify_nlsp_linear_rectangular_3d_fem_refinement as workflow


@pytest.fixture(scope='module')
def config():
    return workflow.fem.read_json(workflow.CONFIG)


@pytest.fixture(scope='module')
def parent(config):
    path = workflow.ROOT / config['parent_bundle']
    if not (path / 'manifest.json').is_file():
        pytest.skip('Immutable FEM-1 bundle is absent; tests never recreate it')
    before = workflow.fem.sha(path / 'manifest.json')
    summary, pre, profiles = workflow.load_parent_bundle(path, config['parent_manifest_sha256'])
    assert workflow.fem.sha(path / 'manifest.json') == before
    return path, summary, pre, profiles


def forbidden_numerics(*args, **kwargs):
    raise AssertionError('Continuation test attempted a new scientific computation')


def block_numerics(monkeypatch):
    for name in ('run_refinement',):
        monkeypatch.setattr(workflow, name, forbidden_numerics)
    for name in ('run_job', 'run_fem', 'build_preflight', 'write_preflight'):
        monkeypatch.setattr(workflow.fem, name, forbidden_numerics)
    monkeypatch.setattr(workflow.fem.single, 'generate_mesh_with_gmsh_cli', forbidden_numerics)
    monkeypatch.setattr(workflow.fem.single, 'run_calculix', forbidden_numerics)
    monkeypatch.setattr(workflow.fem.mh, 'finite_roots', forbidden_numerics)
    monkeypatch.setattr(np.linalg, 'eig', forbidden_numerics)
    monkeypatch.setattr(np.linalg, 'eigh', forbidden_numerics)


def test_authorized_continuation_config(config):
    assert workflow.validate_config(copy.deepcopy(config)) == config
    assert config['geometry'] == {'L': 1., 'b': .2, 'h': .1}
    assert config['material'] == {'E': 1., 'rho': 1., 'nu': .3, 'kappa': 5 / 6}
    assert config['mesh_level'] == {'name': 'refined', 'target_size': .02, 'through_h_nominal': 5}
    assert config['element'] == 'C3D10'
    assert config['requested_eigenpairs'] == 24
    assert config['frozen_omega_window'] == 3.4651360027859885
    assert config['numerical_mesh_convergence_relative'] == .001
    assert config['threads'] == 1
    assert config['job_timeout_seconds'] <= 900
    assert config['job_memory_limit_bytes'] <= 4 * 1024**3
    assert config['numerical_budget_seconds'] <= 1200
    assert config['semantics'] == {'linear_only': True, 'model_fitting': False,
        'new_1d_solutions': False, 'new_nonlinear_jobs': False,
        'automatic_mesh_or_eigenpair_extension': False}
    assert 'mesh_levels' not in config
    assert 'one_extension_eigenpairs' not in config


@pytest.mark.parametrize('path,value', [
    (('geometry', 'L'), 2.), (('geometry', 'b'), .1), (('geometry', 'h'), .12),
    (('material', 'E'), 2.), (('material', 'rho'), 2.), (('material', 'nu'), .2),
    (('material', 'kappa'), .8), (('element',), 'C3D4'),
    (('mesh_level', 'target_size'), .025), (('mesh_level', 'name'), 'fine'),
    (('mesh_level', 'through_h_nominal'), 4), (('requested_eigenpairs',), 36),
    (('frozen_omega_window',), 4.), (('numerical_mesh_convergence_relative',), .002),
    (('threads',), 2), (('numerical_budget_seconds',), 1201),
    (('job_timeout_seconds',), 901), (('job_memory_limit_bytes',), 4 * 1024**3 + 1),
    (('semantics', 'linear_only'), False), (('semantics', 'model_fitting'), True),
    (('semantics', 'new_1d_solutions'), True), (('semantics', 'new_nonlinear_jobs'), True),
    (('semantics', 'automatic_mesh_or_eigenpair_extension'), True),
])
def test_unauthorized_scope_changes_are_rejected(config, path, value):
    changed = copy.deepcopy(config)
    target = changed
    for name in path[:-1]:
        target = target[name]
    target[path[-1]] = value
    with pytest.raises(ValueError):
        workflow.validate_config(changed)


def test_parent_manifest_hash_mismatch_blocks_loading(parent):
    path, _, _, _ = parent
    with pytest.raises(ValueError):
        workflow.load_parent_bundle(path, '0' * 64)


def test_parent_source_and_old_configuration_are_immutable(parent):
    path, summary, pre, profiles = parent
    manifest = workflow.fem.read_json(path / 'manifest.json')
    frozen = workflow.fem.read_json(path / 'frozen_config.json')
    assert frozen == summary['config'] == manifest['identity']['config']
    assert workflow.fem.sha(workflow.fem.__file__) == manifest['identity']['code_sha256']
    assert workflow.fem.sha(workflow.fem.CONFIG) == manifest['identity']['config_sha256']
    assert [level['name'] for level in frozen['mesh_levels']] == ['coarse', 'medium', 'fine']
    assert len(frozen['mesh_levels']) == 3
    assert summary['statuses']['NLSP_FEM1_3D_MESH_CONVERGENCE'] == 'PARTIAL'
    assert summary['extensions_used'] == 0
    assert summary['job_calls'] == {'gmsh': 6, 'ccx': 3}
    assert pre['first_axial_sorted_index'] == 8
    assert len(pre['merged_spectrum']) == len(profiles) == 8
    for relative, digest in frozen['model_hashes'].items():
        assert workflow.fem.sha(workflow.ROOT / relative) == digest


def test_parent_profiles_load_exactly_without_analytic_solutions(parent, config, monkeypatch):
    path, summary, pre, expected = parent
    block_numerics(monkeypatch)
    loaded_summary, loaded_pre, loaded_profiles = workflow.load_parent_bundle(
        path, config['parent_manifest_sha256'])
    assert loaded_summary == summary
    assert loaded_pre == pre
    assert set(loaded_profiles) == set(expected)
    for key in expected:
        np.testing.assert_array_equal(loaded_profiles[key]['x'], expected[key]['x'])
        np.testing.assert_array_equal(loaded_profiles[key]['q'], expected[key]['q'])
        assert loaded_profiles[key]['q'].shape == (401, 7)


def test_geometry_generator_retains_one_monolithic_rectangle_and_algorithm(config):
    text = workflow.fem.rectangular_geo(**config['geometry'], target_size=.02)
    assert 'Box(1)' in text
    assert 'Physical Volume("SOLID",1)' in text
    assert 'Physical Surface("FIXED_LEFT",2)' in text
    assert 'Physical Surface("FIXED_RIGHT",3)' in text
    assert 'Mesh.ElementOrder = 2' in text
    assert 'Cylinder(' not in text
    assert 'RIGID BODY' not in text
    assert 'NLGEOM' not in text
    assert 'DYNAMIC' not in text
    assert 'STATIC' not in text


def test_shape_assignment_uses_mass_shapes_and_keeps_duplicate_ambiguity():
    accepted = workflow.fem.nlsp_shape_assignment([[.01, .99], [.98, .02]])
    assert [row['column'] for row in accepted] == [1, 0]
    assert all(row['status'] == 'MATCHED' for row in accepted)
    duplicate = workflow.fem.nlsp_shape_assignment([[.95, .01], [.94, .02]])
    assert all(row['conflict'] for row in duplicate)
    assert all(row['status'] == 'DUPLICATE_INDEPENDENT_MATCH' for row in duplicate)
    ambiguous = workflow.fem.nlsp_shape_assignment([[.51, .49], [.49, .51]])
    assert all(row['ambiguous'] for row in ambiguous)
    assert all(row['status'] != 'MATCHED' for row in ambiguous)


def test_shared_frequency_difference_retains_signed_3d_denominator():
    result = workflow.fem.nlsp_relative_frequency_difference(3., 4.)
    assert result['signed_relative_difference'] == -.25
    assert result['absolute_relative_difference'] == .25


def test_shared_complete_nodal_vector_guard():
    vectors = {1: {3: (1., 0., 0.), 5: (0., 1., 0.)}}
    valid = workflow.fem.nlsp_validate_nodal_modes([3, 5], vectors, [1])
    np.testing.assert_array_equal(valid[1], [[1., 0., 0.], [0., 1., 0.]])
    with pytest.raises(ValueError, match='missing_nodes'):
        workflow.fem.nlsp_validate_nodal_modes([3, 5, 9], vectors, [1])
    with pytest.raises(ValueError, match='Missing FEM'):
        workflow.fem.nlsp_validate_nodal_modes([3, 5], vectors, [1, 2])


def synthetic_four_mesh_evidence(frequencies=(1100., 1050., 1001., 1000.)):
    pre = {'merged_spectrum': []}
    source = {'meshes': {name: {'modal': {'matches': []}} for name in workflow.OLD_LEVELS}}
    case = {'modal': {'matches': []}}
    cross = {'rows': [], 'all_consistent': True}
    families = [('inplane_bending', 1), ('outplane_bending', 1), ('torsion', 1),
                ('inplane_bending', 2), ('outplane_bending', 2), ('inplane_bending', 3),
                ('torsion', 2), ('axial_mh', 1)]
    for index, (family, local) in enumerate(families, 1):
        reference = {'sorted_index': index, 'family': family, 'local_mode': local, 'omega': 900. * index}
        pre['merged_spectrum'].append(reference)
        for level, frequency in zip(workflow.OLD_LEVELS, frequencies[:3]):
            source['meshes'][level]['modal']['matches'].append({'omega_3d': frequency * index})
        case['modal']['matches'].append({'status': 'MATCHED', 'omega_3d': frequencies[3] * index,
            'fem_mode': index, 'mac': .99, 'margin': .98,
            'section_residual_fraction': .01, 'axial_warp_fraction': .005})
        cross['rows'].append({'status': 'MATCHED', 'agrees_with_independent_1d_assignment': True,
                             'mac': .999, 'margin': .998})
    return pre, source, case, cross


def test_four_mesh_formula_uses_refined_denominator_and_keeps_all_levels():
    rows = workflow.compare_four_meshes(*synthetic_four_mesh_evidence())
    assert len(rows) == 8
    row = rows[0]
    assert row['fine_refined_relative'] == .001
    assert row['fine_refined_absolute'] == 1.
    assert row['omega_coarse'] == 1100.
    assert row['omega_medium'] == 1050.
    assert row['omega_fine'] == 1001.
    assert row['omega_refined'] == 1000.
    assert row['signed_relative_difference'] == -.1
    assert row['absolute_relative_difference'] == .1
    assert row['mesh_status'] == 'MESH_ACCEPTED_AT_PRESET_TOLERANCE'
    assert row['frequency_monotone'] is True
    assert row['absolute_changes_decrease'] is True
    assert row['regularity'] == 'REGULAR_OBSERVED_TREND'
    previous = abs(900. - 1001.) / 1001.
    assert row['model_difference_change_percentage_points'] == pytest.approx(100 * (.1 - previous))
    assert all(row['mesh_status'] == 'MESH_ACCEPTED_AT_PRESET_TOLERANCE' for row in rows)


def test_mesh_criterion_is_not_loosened_above_point_one_percent():
    rows = workflow.compare_four_meshes(*synthetic_four_mesh_evidence((1100., 1050., 1001., 999.99)))
    assert all(row['fine_refined_relative'] > .001 for row in rows)
    assert all(row['mesh_status'] == 'MESH_UNRESOLVED' for row in rows)


def test_small_last_change_with_deteriorating_successive_change_is_visible():
    rows = workflow.compare_four_meshes(*synthetic_four_mesh_evidence((1003., 1000.6, 1000.5, 1000.)))
    assert all(row['fine_refined_relative'] < .001 for row in rows)
    assert all(row['fine_refined_relative'] > row['medium_fine_relative'] for row in rows)
    assert all(row['mesh_status'] == 'MESH_UNRESOLVED' for row in rows)
    assert all(row['absolute_changes_decrease'] is False for row in rows)


def test_nonmonotone_frequencies_are_not_hidden_as_regular_convergence():
    rows = workflow.compare_four_meshes(*synthetic_four_mesh_evidence((1003., 1001., 1002., 1001.5)))
    assert all(row['frequency_monotone'] is False for row in rows)
    assert all(row['regularity'] == 'IRREGULAR_OBSERVED_TREND' for row in rows)


@pytest.mark.parametrize('failure', ('new_match', 'cross_match', 'cross_identity'))
def test_unresolved_shape_cannot_create_a_frequency_comparison(failure):
    pre, source, case, cross = synthetic_four_mesh_evidence()
    if failure == 'new_match':
        case['modal']['matches'][0]['status'] = 'AMBIGUOUS_SHAPE_OR_SUBSPACE'
    elif failure == 'cross_match':
        cross['rows'][0]['status'] = 'LOW_MAC'
    else:
        cross['rows'][0]['agrees_with_independent_1d_assignment'] = False
    rows = workflow.compare_four_meshes(pre, source, case, cross)
    assert rows[0]['mesh_status'] == 'MATCH_UNRESOLVED'
    assert 'omega_refined' not in rows[0]
    assert len(rows) == 8
    assert all(row['mesh_status'] == 'MESH_ACCEPTED_AT_PRESET_TOLERANCE' for row in rows[1:])


def test_changed_parent_artifact_is_rejected_without_scientific_recovery(tmp_path, parent, monkeypatch):
    _, source, _, _ = parent
    fake = tmp_path / 'historical'
    fake.mkdir()
    workflow.write_json(fake / 'summary.json', source)
    digest = workflow.sha(fake / 'summary.json')
    workflow.write_json(fake / 'manifest.json', {'identity': {}, 'artifact_hashes': {'summary.json': digest}})
    manifest_digest = workflow.sha(fake / 'manifest.json')
    workflow.write_json(fake / 'summary.json', {'changed': 'invalid historical data'})
    block_numerics(monkeypatch)
    with pytest.raises(ValueError, match='artifact hash mismatch'):
        workflow.load_parent_bundle(fake, manifest_digest)


def make_synthetic_cache(path, config, *, item=None):
    path.mkdir(parents=True)
    item = {'schema': 'synthetic-cache-unit-test'} if item is None else item
    summary = {'config': config, 'comparisons': [], 'statuses': {
        'NLSP_FEM1R_MODAL_EXECUTION': 'NOT_RUN'},
        'scope': 'Synthetic cache test; no actual frequencies or jobs'}
    workflow.write_json(path / 'summary.json', summary)
    workflow.write_json(path / 'manifest.json', {'identity': item,
        'artifact_hashes': {'summary.json': workflow.sha(path / 'summary.json')}})
    return item, summary


def test_continuation_cache_rejects_identity_and_corrupt_artifact(tmp_path, config, monkeypatch):
    path = tmp_path / 'cache'
    item, expected = make_synthetic_cache(path, config)
    monkeypatch.setattr(workflow, 'load_parent_bundle', lambda *a, **kw: (None, None, None))
    assert workflow.validate_cache(path, item) == expected
    with pytest.raises(ValueError, match='cache identity mismatch'):
        workflow.validate_cache(path, {'changed': True})
    workflow.write_json(path / 'summary.json', {'corrupted': True})
    with pytest.raises(ValueError, match='artifact hash mismatch'):
        workflow.validate_cache(path, item)


@pytest.mark.parametrize('mode', ('--report-only', '--plot-only'))
def test_report_and_plot_cache_replay_have_zero_numerical_calls(tmp_path, config, monkeypatch, mode):
    path = tmp_path / 'cache'
    _, expected = make_synthetic_cache(path, config)
    monkeypatch.setattr(workflow, 'load_parent_bundle', lambda *a, **kw: (None, None, None))
    monkeypatch.setattr(workflow, 'identity', forbidden_numerics)
    block_numerics(monkeypatch)
    assert workflow.main([mode, str(path)]) == expected


def test_matching_compute_cache_makes_zero_mesh_modal_or_analytic_calls(tmp_path, config, monkeypatch):
    output = tmp_path / 'cache_root'
    key, item = 'synthetic-fingerprint', {'schema': 'synthetic-cache-unit-test'}
    _, expected = make_synthetic_cache(output / key, config, item=item)
    monkeypatch.setattr(workflow, 'identity', lambda *a, **kw: (key, item))
    monkeypatch.setattr(workflow, 'load_parent_bundle', lambda *a, **kw: (None, None, None))
    block_numerics(monkeypatch)
    assert workflow.main(['--run-fem', '--output-dir', str(output)]) == expected


def test_unmanifested_partial_attempt_does_not_start_a_retry(tmp_path, config, monkeypatch):
    output = tmp_path / 'cache_root'
    key, item = 'synthetic-fingerprint', {'schema': 'synthetic-cache-unit-test'}
    target = output / key
    target.mkdir(parents=True)
    (target / 'partial.txt').write_text('Interrupted synthetic attempt')
    monkeypatch.setattr(workflow, 'identity', lambda *a, **kw: (key, item))
    monkeypatch.setattr(workflow, 'load_parent_bundle', lambda *a, **kw: (None, None, None))
    block_numerics(monkeypatch)
    with pytest.raises(RuntimeError, match='no automatic retry'):
        workflow.main(['--run-fem', '--output-dir', str(output)])


def test_failed_generation_retains_diagnostics_and_never_runs_calculix(tmp_path, config, parent, monkeypatch):
    _, source, pre, profiles = parent
    calls = []
    def failed_mesh(*args, **kwargs):
        calls.append('mesh')
        return False, 'Synthetic generation failure', []
    monkeypatch.setattr(workflow.fem.single, 'generate_mesh_with_gmsh_cli', failed_mesh)
    monkeypatch.setattr(workflow.fem.single, 'run_calculix', forbidden_numerics)
    monkeypatch.setattr(workflow.fem, 'run_job', forbidden_numerics)
    monkeypatch.setattr(workflow.fem, 'build_preflight', forbidden_numerics)
    single = workflow.fem.single
    names = ('L', 'E', 'RHO', 'NU', 'SOLID_MODES_REQUESTED')
    original = {name: getattr(single, name) for name in names}
    original_subprocess = single.subprocess.run
    result = workflow.run_refinement(copy.deepcopy(config), tmp_path, source, pre, profiles)
    assert calls == ['mesh']
    assert result['job_calls'] == {'gmsh': 0, 'ccx': 0}
    assert result['refined_case']['status'] == 'FAIL'
    assert result['refined_case']['failure'] == 'Synthetic generation failure'
    assert result['statuses']['NLSP_FEM1R_NEW_MESH_QUALITY'] == 'FAIL'
    assert result['statuses']['NLSP_FEM1R_MODAL_EXECUTION'] == 'NOT_RUN'
    assert result['comparisons'] == []
    assert (tmp_path / 'summary.json').is_file()
    assert (tmp_path / 'meshes/refined/case.json').is_file()
    assert {name: getattr(single, name) for name in names} == original
    assert single.subprocess.run is original_subprocess


def test_failed_mesh_quality_blocks_the_only_modal_call(tmp_path, config, parent, monkeypatch):
    _, source, pre, profiles = parent
    monkeypatch.setattr(workflow.fem.single, 'generate_mesh_with_gmsh_cli', lambda *a, **kw: (True, 'Synthetic only', []))
    monkeypatch.setattr(workflow.fem, 'audit_rectangular_mesh', lambda *a, **kw: ({'status': 'FAIL'}, None))
    monkeypatch.setattr(workflow.fem.single, 'write_calculix_template', forbidden_numerics)
    monkeypatch.setattr(workflow.fem.single, 'run_calculix', forbidden_numerics)
    monkeypatch.setattr(workflow.fem, 'run_job', forbidden_numerics)
    result = workflow.run_refinement(copy.deepcopy(config), tmp_path, source, pre, profiles)
    assert result['refined_case']['status'] == 'FAIL'
    assert result['statuses']['NLSP_FEM1R_NEW_MESH_QUALITY'] == 'FAIL'
    assert result['statuses']['NLSP_FEM1R_MODAL_EXECUTION'] == 'NOT_RUN'
    assert result['job_calls']['ccx'] == 0
    assert result['comparisons'] == []


def test_cache_identity_explicitly_covers_parent_helpers_binaries_and_policy(config, parent):
    key, item = workflow.identity()
    assert len(key) == 16
    assert item['config'] == config
    assert item['parent_bundle'] == config['parent_bundle']
    assert item['parent_manifest_sha256'] == config['parent_manifest_sha256']
    assert item['executables'] == {name: workflow.sha(config[name]) for name in ('gmsh_exe', 'ccx_exe')}
    assert item['code_sha256'] == workflow.sha(workflow.__file__)
    assert item['helper_sha256']['scripts/analysis/verify_nlsp_linear_rectangular_3d_fem.py'] == workflow.sha(workflow.fem.__file__)
    assert item['python']
    assert set(item['dependencies']) == {'numpy', 'scipy', 'matplotlib'}
    assert item['config']['mesh_level']['target_size'] == .02
    assert item['config']['requested_eigenpairs'] == 24
    assert item['config']['numerical_mesh_convergence_relative'] == .001
    assert item['config']['matching'] == parent[1]['config']['matching']


def test_cache_identity_changes_if_config_bytes_change(tmp_path, config):
    first = tmp_path / 'first.json'
    second = tmp_path / 'second.json'
    workflow.write_json(first, config)
    # Same physical/matching policy, visibly different provenance input bytes.
    changed = copy.deepcopy(config)
    changed['synthetic_cache_test_provenance'] = 'unit test only; no computation'
    workflow.write_json(second, changed)
    key1, item1 = workflow.identity(first)
    key2, item2 = workflow.identity(second)
    assert key1 != key2
    assert item1['config_sha256'] != item2['config_sha256']
    assert item1['parent_manifest_sha256'] == item2['parent_manifest_sha256']
    assert item1['executables'] == item2['executables']


def test_guard_blocks_a_second_modal_job_even_from_a_bad_helper(tmp_path, config, parent, monkeypatch):
    import subprocess
    _, source, pre, profiles = parent
    monkeypatch.setattr(workflow.fem.single, 'generate_mesh_with_gmsh_cli', lambda *a, **kw: (True, 'Synthetic only', []))
    monkeypatch.setattr(workflow.fem, 'audit_rectangular_mesh', lambda *a, **kw: ({'status': 'PASS'}, None))
    monkeypatch.setattr(workflow.fem.single, 'write_calculix_template', lambda *a, **kw: None)
    calls = []
    def fake_job(command, *args, **kwargs):
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, '', ''), {'failure': None,
            'returncode': 0, 'seconds': .001, 'peak_working_set_bytes': 0}
    monkeypatch.setattr(workflow.fem, 'run_job', fake_job)
    def bad_modal_helper(paths, executable, *args):
        workflow.fem.single.subprocess.run([executable], cwd=paths.case_dir)
        workflow.fem.single.subprocess.run([executable], cwd=paths.case_dir)
        raise AssertionError('Second fake modal job should have been blocked')
    monkeypatch.setattr(workflow.fem.single, 'run_calculix', bad_modal_helper)
    result = workflow.run_refinement(copy.deepcopy(config), tmp_path, source, pre, profiles)
    assert len(calls) == 1
    assert result['job_calls'] == {'gmsh': 0, 'ccx': 1}
    assert result['refined_case']['status'] == 'FAIL'
    assert result['refined_case']['failure'] == 'Only one new CalculiX job authorized'
    assert result['statuses']['NLSP_FEM1R_MODAL_EXECUTION'] == 'FAIL'
    assert result['statuses']['NLSP_FEM1R_MODE_IDENTIFICATION'] == 'NOT_RUN'


@pytest.fixture(scope='module')
def actual_refinement(config):
    available = []
    for path in workflow.OUTPUT.glob('*'):
        if not (path / 'manifest.json').is_file():
            continue
        manifest = workflow.read_json(path / 'manifest.json')
        item = manifest['identity']
        if (item.get('code_sha256') == workflow.sha(workflow.__file__)
                and item.get('config') == config):
            available.append(path)
    if not available:
        pytest.skip('Saved FEM-1R result is absent; tests never generate it')
    assert len(available) == 1, 'More than one matching bounded continuation bundle'
    path = available[0]
    summary = workflow.validate_cache(path)
    return path, summary


def test_actual_exactly_one_level_one_modal_job_and_no_extensions(actual_refinement, config, parent):
    path, summary = actual_refinement
    assert summary['config'] == config
    assert summary['job_calls'] == {'gmsh': 2, 'ccx': 1}
    assert summary['parent_reference']['one_D_recomputed'] is False
    assert summary['parent_reference']['old_FEM_jobs_repeated'] is False
    assert summary['parent_reference']['manifest_sha256'] == config['parent_manifest_sha256']
    assert sorted(p.name for p in (path / 'meshes').iterdir()) == ['refined']
    assert summary['refined_case']['target_size'] == .02
    assert summary['refined_case']['status'] == 'PASS'
    assert len(summary['refined_case']['jobs']) == 3
    assert all(job['returncode'] == 0 and job['failure'] is None for job in summary['refined_case']['jobs'])
    assert summary['runtime']['numerical_seconds'] <= summary['runtime']['budget_seconds'] == 1200
    assert workflow.sha(parent[0] / 'manifest.json') == config['parent_manifest_sha256']
    assert workflow.fem.validate_cache(parent[0]) == parent[1]


def test_actual_refined_mesh_quality_geometry_and_resolution(actual_refinement):
    _, summary = actual_refinement
    audit = summary['refined_case']['mesh_audit']
    assert audit['status'] == 'PASS'
    assert audit['nodes'] == 20752
    assert audit['c3d10_elements'] == 12687
    assert audit['solid_element_types'] == ['C3D10']
    assert audit['bbox_matches'] is True
    assert audit['volume'] == pytest.approx(.02, rel=1e-12)
    assert audit['mass'] == pytest.approx(.02, rel=1e-12)
    assert audit['volume_face_connected_components'] == 1
    assert audit['nonconforming_quadratic_face_count'] == 0
    assert audit['negative_or_zero_jacobian_elements'] == 0
    assert audit['minimum_quadratic_jacobian'] > 0
    assert audit['quality_min_median_max'][0] > 0
    assert audit['straight_midpoint_geometry'] is True
    assert audit['rigid_body_kinematic_constraint_check'] == 'PASS'
    assert audit['fixed_left_count'] > 0 and audit['fixed_right_count'] > 0
    assert audit['gmsh_warnings_errors'] == []
    for direction, expected in [('thickness_h', 5), ('width_b', 10)]:
        edges = audit['actual_resolution'][direction]['box_edge_resolution']['edges']
        assert len(edges) == 4
        assert all(edge['linear_segment_count'] == expected for edge in edges)


def test_actual_linear_input_has_two_complete_faces_and_no_internal_joint(actual_refinement):
    path, _ = actual_refinement
    case = path / 'meshes/refined'
    text = (case / 'modal.inp').read_text(encoding='utf-8').upper()
    solid = (case / 'solid_mesh.inp').read_text(encoding='utf-8').upper()
    assert 'LEFT_FIXED, 1, 3, 0' in text
    assert 'RIGHT_FIXED, 1, 3, 0' in text
    assert '*ELASTIC' in text and '*DENSITY' in text
    assert '*FREQUENCY\n24\n' in text
    assert '*NODE FILE\nU\n' in text
    assert 'TYPE=C3D10' in solid
    for prohibited in ('NLGEOM', '*STATIC', '*DYNAMIC', '*MPC', '*RIGID BODY', '*SPRING', '*DAMPING'):
        assert prohibited not in text
        assert prohibited not in solid
    expected = workflow.fem.rectangular_geo(L=1., b=.2, h=.1, target_size=.02)
    assert (case / 'rod.geo').read_text(encoding='utf-8') == expected + '\nGeneral.NumThreads=1;\nMesh.RandomSeed=1;\n'


def test_actual_24_positive_frequencies_full_vectors_units_and_clamps(actual_refinement):
    path, summary = actual_refinement
    case = summary['refined_case']
    modal, audit = case['modal'], case['mesh_audit']
    assert modal['eigenpairs'] == modal['parsed_eigenvectors'] == 24
    assert modal['full_frequency_window_covered'] is True
    assert modal['maximum_clamp_relative_residual'] == 0.
    assert modal['frequency_unit'] == 'rad per normalized time'
    assert modal['maximum_printed_frequency_unit_error'] < modal['unit_gate_relative'] == 1e-5
    rows = modal['raw_frequencies']
    assert len(rows) == 24
    assert len({row['raw_mode_number'] for row in rows}) == 24
    assert np.all(np.diff([row['angular_frequency'] for row in rows]) > 0)
    assert rows[-1]['angular_frequency'] > summary['config']['frozen_omega_window']
    with np.load(path / 'meshes/refined/modal_vectors.npz', allow_pickle=False) as vectors:
        ids = vectors['node_ids']
        assert len(ids) == audit['nodes']
        assert len(np.unique(ids)) == len(ids)
        assert vectors['nodes'].shape == (len(ids), 3)
        fixed = np.isin(ids, audit['fixed_left_ids'] + audit['fixed_right_ids'])
        assert np.count_nonzero(fixed) == audit['fixed_left_count'] + audit['fixed_right_count']
        vector_keys = [key for key in vectors.files if key.startswith('U_')]
        assert len(vector_keys) == 24
        for row in rows:
            assert row['eigenvalue'] > 0
            assert abs(row['angular_frequency']**2 / row['eigenvalue'] - 1) < 1e-5
            assert abs(2 * math.pi * row['cyclic_frequency'] / row['angular_frequency'] - 1) < 1e-5
            value = vectors['U_' + str(row['raw_mode_number'])]
            assert value.shape == (len(ids), 3)
            assert np.all(np.isfinite(value))
            assert np.max(np.abs(value)) > 0
            assert np.max(np.abs(value[fixed])) == 0


def test_actual_eight_unique_shape_matches_keep_family_and_full_mass_evidence(actual_refinement):
    path, summary = actual_refinement
    modal = summary['refined_case']['modal']
    assert modal['all_eight_identified'] is True
    assert modal['matching_uses_frequencies'] is False
    assert modal['additional_modes_in_window'] == []
    matches = modal['matches']
    assert len(matches) == 8
    assert len({row['fem_mode'] for row in matches}) == 8
    assert {row['family'] for row in matches} == {'axial_mh', 'inplane_bending', 'outplane_bending', 'torsion'}
    policy = summary['config']['matching']
    for row in matches:
        assert row['status'] == 'MATCHED'
        assert row['fem_family'] == row['family']
        assert row['mac'] >= policy['minimum_mac']
        assert row['margin'] >= policy['minimum_margin']
        assert row['conflict'] is False
        assert row['low_mac'] is False
        assert row['ambiguous'] is False
        assert row['fem_mass_norm'] > 0
        assert row['one_D_lift_mass_norm'] > 0
        assert 0 <= row['section_residual_fraction'] <= 1
        assert 0 <= row['axial_warp_fraction'] <= 1
    with np.load(path / 'meshes/refined/modal_vectors.npz', allow_pickle=False) as vectors:
        assignment = workflow.fem.nlsp_shape_assignment(vectors['MAC'], policy['minimum_mac'], policy['minimum_margin'])
    raw_ids = [row['raw_mode_number'] for row in modal['raw_frequencies']]
    assert [raw_ids[row['column']] for row in assignment] == [row['fem_mode'] for row in matches]
    assert all(row['status'] == 'MATCHED' for row in assignment)


def test_actual_cross_mesh_matching_reproduces_shape_only_consistency(actual_refinement, parent):
    path, summary = actual_refinement
    cross = workflow.read_json(path / 'fine_refined_shape_correspondence.json')
    recalculated = workflow.cross_mesh_matching(parent[0], path / 'meshes/refined', parent[2], summary['config']['matching'])
    assert recalculated == cross
    assert cross['all_consistent'] is True
    assert len(cross['rows']) == 8
    assert len({row['refined_mode'] for row in cross['rows']}) == 8
    assert all(row['status'] == 'MATCHED' and row['agrees_with_independent_1d_assignment'] for row in cross['rows'])
    assert all(row['mac'] >= summary['config']['matching']['minimum_mac'] for row in cross['rows'])


def test_actual_cross_mesh_identification_does_not_use_frequency_values(tmp_path, actual_refinement, parent):
    path, summary = actual_refinement
    original = workflow.read_json(path / 'meshes/refined/modal_analysis.json')
    altered = copy.deepcopy(original)
    for row in altered['raw_frequencies']:
        row['angular_frequency'] *= 1000.
        row['cyclic_frequency'] *= 1000.
        row['eigenvalue'] *= 1e6
    for row in altered['matches']:
        row['omega_3d'] *= 1000.
    workflow.write_json(tmp_path / 'modal_analysis.json', altered)
    with np.load(path / 'meshes/refined/modal_vectors.npz', allow_pickle=False) as data:
        profiles_only = {key: data[key].copy() for key in data.files if key.startswith(('profile_x_', 'profile_q_'))}
    np.savez_compressed(tmp_path / 'modal_vectors.npz', **profiles_only)
    actual = workflow.cross_mesh_matching(parent[0], tmp_path, parent[2], summary['config']['matching'])
    expected = workflow.read_json(path / 'fine_refined_shape_correspondence.json')
    assert actual == expected


def test_actual_four_mesh_rows_use_frozen_1d_and_preserve_signed_differences(actual_refinement, parent):
    path, summary = actual_refinement
    cross = workflow.read_json(path / 'fine_refined_shape_correspondence.json')
    rows = summary['comparisons']
    assert workflow.compare_four_meshes(parent[2], parent[1], summary['refined_case'], cross) == rows
    assert len(rows) == 8
    threshold = summary['config']['numerical_mesh_convergence_relative']
    for index, row in enumerate(rows):
        assert row['omega_1d'] == parent[2]['merged_spectrum'][index]['omega']
        for level in workflow.OLD_LEVELS:
            assert row['omega_' + level] == parent[1]['meshes'][level]['modal']['matches'][index]['omega_3d']
        fine, refined = row['omega_fine'], row['omega_refined']
        assert row['fine_refined_relative'] == pytest.approx(abs(refined - fine) / refined, rel=1e-14)
        assert row['signed_relative_difference'] == pytest.approx((row['omega_1d'] - refined) / refined, rel=1e-14)
        assert row['absolute_relative_difference'] == abs(row['signed_relative_difference'])
        assert row['absolute_changes_decrease'] is True
        assert row['frequency_monotone'] is True
        assert row['regularity'] == 'REGULAR_OBSERVED_TREND'
        accepted = row['fine_refined_relative'] <= threshold and row['fine_refined_relative'] <= row['medium_fine_relative']
        expected = 'MESH_ACCEPTED_AT_PRESET_TOLERANCE' if accepted else 'MESH_UNRESOLVED'
        assert row['mesh_status'] == expected
    all_accepted = all(row['mesh_status'] == 'MESH_ACCEPTED_AT_PRESET_TOLERANCE' for row in rows)
    assert summary['statuses']['NLSP_FEM1R_MESH_CONVERGENCE'] == ('PASS' if all_accepted else 'PARTIAL')
    assert summary['statuses']['NLSP_FEM1R_ALL_FAMILY_COMPARISON'] == ('PASS' if all_accepted else 'PARTIAL')


@pytest.mark.parametrize('mode', ('--report-only', '--plot-only'))
def test_actual_readonly_replay_has_zero_scientific_calls(actual_refinement, monkeypatch, mode):
    path, expected = actual_refinement
    monkeypatch.setattr(workflow, 'identity', forbidden_numerics)
    block_numerics(monkeypatch)
    from matplotlib.figure import Figure
    saved = []
    monkeypatch.setattr(Figure, 'savefig', lambda self, target, *args, **kwargs: saved.append(Path(target).name))
    assert workflow.main([mode, str(path)]) == expected
    if mode == '--plot-only':
        assert sorted(saved) == sorted(name + '.' + extension for name in
            ('four_mesh_convergence', 'refined_model_difference_and_mesh_change') for extension in ('pdf', 'png'))
    else:
        assert saved == []
