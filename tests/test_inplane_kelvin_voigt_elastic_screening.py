"""Sparse elastic screening; never run the spectral grid or a damped solve."""
import json
from pathlib import Path

import numpy as np
import pytest

from scripts.analysis.laminated_beams import screen_inplane_kelvin_voigt_elastic as flow
from scripts.lib import inplane_kelvin_voigt as kv


@pytest.fixture(scope='module')
def saved():
    if not (flow.OUTPUT/'diagnostics.json').exists():
        pytest.skip('local scientific screening artifacts are not installed')
    return flow.read_json(flow.OUTPUT/'diagnostics.json')


@pytest.mark.parametrize('amplitude', [.031, -17., 2+3j])
def test_participation_invariant_to_amplitude_and_phase(amplitude):
    a = flow.participation(4-3j, 2., 7.)
    b = flow.participation(amplitude*(4-3j), 2., abs(amplitude)**2*7.)
    assert b == pytest.approx(a, rel=2e-15)
    assert a['s_joint'] == pytest.approx(flow.REF.D * 25/28)
    assert a['a_slope_pred'] == pytest.approx(flow.REF.m * 25/14)


def test_exact_reflection_relation_implies_zero_joint_rotation():
    first = np.arange(30.).reshape(5, 6)
    second = -first @ flow.modes.REFLECTION
    assert np.array_equal(first[:, 2], second[:, 2])
    assert first[-1, 2]-second[-1, 2] == 0
    assert flow.classify(-1, 0., identical_arms=True, root_confirmed=True) == 'EXACT_INACTIVE_BY_SYMMETRY'


def test_small_active_participation_is_not_structural_inactivity():
    assert flow.participation(1e-20, 1., 1.)['s_joint'] > 0
    assert flow.classify(1, 0., identical_arms=True, root_confirmed=True) == 'ACTIVE'
    assert flow.classify(-1, 0., identical_arms=False, root_confirmed=True) == 'ACTIVE'
    with pytest.raises(ValueError, match='unconfirmed'):
        flow.classify(-1, .01, identical_arms=True, root_confirmed=True)


def test_fixed_configuration_and_no_damped_parameter():
    p, c = flow.configuration()
    assert c['beta_deg'] == [0., 5., 45., 75.]
    assert c['model'] == 'EB' and c['d_theta'] == c['mu'] == 0
    assert c['k_theta'] == flow.REF.D != p.D
    ctx = flow.ElasticMatrices(5., p)
    assert all(a.rotational_mass == 0 for a in ctx.pair)
    with pytest.raises(ValueError, match='real elastic'):
        ctx.assembly(.1+.01j)
    with pytest.raises(ValueError, match='angular refinement'):
        flow.ElasticMatrices(5.1, p)


def test_beta5_participation_and_read_only_k12_control(saved):
    a, b = saved['groups']['5.0']['rows'][:2]
    assert a['activity_status'] == 'ACTIVE' and a['abs_Delta_psi'] > 1
    assert min(a[k] for k in ('s_joint', 'a_slope_pred', 'zeta_slope_pred')) > 0
    assert b['activity_status'] == 'EXACT_INACTIVE_BY_SYMMETRY'
    assert b['abs_Delta_psi'] < 1e-8 and b['s_joint'] < 1e-20
    control = flow.k12_control([a, b])
    assert control['ACTIVE']['accepted']
    assert control['ACTIVE']['relative_difference'] < 1e-3
    assert control['INACTIVE']['accepted'] and control['INACTIVE']['max_abs_a'] < 1e-8
    assert control['INACTIVE']['max_abs_Delta'] < 1e-8


def test_first_six_sorted_positions_not_old_bending_branches(saved):
    for beta, group in saved['groups'].items():
        rows = group['rows']
        assert len(rows) == 6 and [r['sorted_mode'] for r in rows] == [f'sorted_{i:02d}' for i in range(1, 7)]
        assert all(a['Omega'] < b['Omega'] for a, b in zip(rows, rows[1:]))
        expected = [2, 4, 6] if float(beta) <= 5 else [1, 4, 6]
        assert [i for i, r in enumerate(rows, 1) if r['activity_status'] == 'EXACT_INACTIVE_BY_SYMMETRY'] == expected
        assert group['guard']['Omega'] > rows[-1]['Omega']
    # The sixth beta=0 root is axial, not the seventh (bending) source root.
    r = saved['groups']['0.0']['rows'][5]
    with np.load(flow.OUTPUT/'screening_shapes.npz') as shapes:
        y = shapes[r['shape_key']+'__states']
    assert np.linalg.norm(y[:, :, 1:3]) < 1e-9*np.linalg.norm(y[:, :, 0])


def test_forms_real_global_mass_and_gates(saved):
    p = flow.mechanics.section()[1]
    pair = flow.mechanics.arms('EB', 0., p)
    _, weights = flow.modes.quadrature()
    with np.load(flow.OUTPUT/'screening_shapes.npz') as shapes:
        for group in saved['groups'].values():
            for r in group['rows']:
                y = shapes[r['shape_key']+'__states']
                assert not np.iscomplexobj(y)
                v = flow.mechanics.physical_vector(y, pair, weights)
                assert np.vdot(v, v) == pytest.approx(1., abs=1e-12)
                assert r['root_residual'] <= 1e-9 and r['sigma_ratio'] <= 1e-9
                assert r['physical_residual'] <= 1e-9 and r['detected_nullity'] == 1
                assert all(abs(q) <= 1e-10 for q in group['physical_conditions'][r['shape_key']][:2])
                scaled = flow.participation(3j*r['Delta_psi'], r['omega'], 9*r['mass_M'])
                assert scaled['s_joint'] == pytest.approx(r['s_joint'], rel=1e-14, abs=1e-40)


def test_reused_frequencies_and_shapes_are_exact_copies(saved):
    old = flow.read_json(flow.SOURCE/'diagnostics.json')['points']
    count = 0
    with np.load(flow.SOURCE/'shapes.npz') as source, np.load(flow.OUTPUT/'screening_shapes.npz') as result:
        for g in saved['groups'].values():
            for r in g['rows']:
                if r['reuse_status'] != 'REUSED_ROOT_AND_FORM':
                    continue
                key = r['source_provenance'].split('#')[1]
                original = next(q for p in old.values() for q in p['roots'] if q['shape_key'] == key)
                assert (r['omega'], r['Omega'], r['Lambda']) == (original['omega'], original['Omega'], original['Lambda'])
                np.testing.assert_array_equal(result[r['shape_key']+'__states'], source[key+'__states'])
                count += 1
    assert count == 14


def test_missing_only_no_computational_calls_or_mutation(saved, monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail('completed screening must perform no computation')
    for name in ('configuration', 'search_missing', 'diagnose_shape', 'ElasticMatrices'):
        monkeypatch.setattr(flow, name, forbidden)
    monkeypatch.setattr(kv, 'correct', forbidden)
    monkeypatch.setattr(kv.Provider, 'matrices', forbidden)
    before = {p.name: flow.sha(p) for p in flow.OUTPUT.iterdir() if p.is_file()}
    result = flow.compute()
    assert all(result[k] == 0 for k in ('new_root_groups', 'new_elastic_solver_calls', 'form_recoveries', 'matrix_builds', 'new_positive_d_complex_roots'))
    assert before == {p.name: flow.sha(p) for p in flow.OUTPUT.iterdir() if p.is_file()}
    assert saved['protected_sources'] == flow.protected_hashes()


def test_positive_d_complex_solver_is_not_in_execution_graph(saved):
    # No KV provider/corrector is imported by this real EB workflow.
    import ast
    source = Path(flow.__file__).read_text()
    tree = ast.parse(source)
    imports = [ast.unparse(node) for node in ast.walk(tree) if isinstance(node, (ast.Import, ast.ImportFrom))]
    assert not any('inplane_kelvin_voigt ' in name or 'pilot_inplane_kelvin_voigt' in name for name in imports)
    assert saved['new_positive_d_complex_roots'] == 0
    manifest = flow.read_json(flow.OUTPUT/'run_manifest.json')
    assert manifest['counters']['new_positive_d_complex_roots'] == 0
    assert all(r['d_theta'] == 0 for g in saved['groups'].values() for r in g['rows'])


def test_source_directory_and_descendants_are_read_only():
    for path in (flow.SOURCE, flow.K12, flow.SOURCE/'nested', flow.K12/'nested'):
        with pytest.raises(ValueError, match='read-only'):
            flow.compute(path)
