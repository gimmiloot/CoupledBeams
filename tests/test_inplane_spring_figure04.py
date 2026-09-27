"""Focused saved-data/figure checks; never solve or recover an eigenmode."""
import sys

import numpy as np
import pytest

from scripts.lib import inplane_spring_figure04 as plot


@pytest.fixture(scope='module')
def saved():
    if not all(p.exists() for p in plot.source_paths()):
        pytest.skip('local K15/K22/K23 artifacts required; no automatic reconstruction')
    return plot.load_saved()


def test_exact_integer_angles_and_saved_matching(saved):
    assert plot.BETAS == (45, 75)
    assert set(saved['data']) == {(b, t) for b in (45, 75) for t in ('EB', 'RLB')}
    for beta in plot.BETAS:
        m = saved['matching'][beta]
        assert m['eb_sorted_mode'] == 'sorted_05' and m['rlb_sorted_mode'] == 'rlb_sorted_05'
        assert m['eta'] == '1' and m['match_status'] == 'CONFIRMED'
        assert float(m['MAC_common']) > .99
        for theory in plot.THEORIES:
            r = saved['data'][(beta, theory)]['row']
            assert float(r['beta_deg']) == beta and float(r['d_theta']) == 0


def test_saved_mass_and_no_independent_shape_rescaling(saved):
    for theory in plot.THEORIES:
        with np.load(plot.BASE/plot.FOLDERS[theory]/plot.SHAPES[theory], allow_pickle=False) as archive:
            for beta in plot.BETAS:
                item = saved['data'][(beta, theory)]
                raw = archive[item['row']['shape_key']+'__states']
                assert item['sign'] in (-1., 1.)
                np.testing.assert_array_equal(item['states'], item['sign']*raw)
                assert abs(plot.saved_mass(item['states'], saved['configurations'][theory], theory)-1) < 1e-10
                if theory == 'RLB':
                    # A read-only check on saved arrays: rotary mass must not be omitted.
                    assert plot.saved_mass(item['states'], saved['configurations'][theory], 'EB') < .999


def test_csv_values_are_exact_source_copies(saved):
    rows = plot.figure_rows(saved)
    assert len(rows) == 4
    for row in rows:
        beta, theory = row['beta_deg'], row['theory']
        source = saved['data'][(beta, theory)]['row']
        for target, original in [('Omega', 'Omega'), ('Lambda', 'Lambda'), ('psi_joint', 'psi1_joint'),
                                 ('Delta_psi', 'Delta_psi'), ('G', 'zeta_slope_pred')]:
            assert row[target] == source[original]
        assert row['zeta_at_reference_d'] == saved['context'][beta]['zeta_'+theory]
        assert row['reference_d_theta'] == .001 and row['elastic_d_theta'] == 0
        assert row['Delta_psi_display'] == pytest.approx(2*row['psi_joint_display'])


def test_geometry_overlay_and_rotation_are_unaltered_fields(saved, monkeypatch):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    def no_derivative(*args, **kwargs):
        pytest.fail('plot saved EB psi=-w prime and independent RLB psi; no finite differences')
    monkeypatch.setattr(np, 'gradient', no_derivative)
    fig = plot.draw(saved)
    for col, beta in enumerate(plot.BETAS):
        top, bottom = fig.axes[col], fig.axes[2+col]
        assert top.get_title() == r'$\beta='+str(beta)+r'^{\circ}$'
        geometry = {line.get_label(): line for line in top.lines}
        rotations = {line.get_label(): line for line in bottom.lines}
        for theory in plot.THEORIES:
            states = saved['data'][(beta, theory)]['states']
            ref, disp = plot.global_centrelines(beta, 0, states)
            for arm in range(2):
                actual = np.array(geometry[f'{theory}_{arm}'].get_data()).T
                np.testing.assert_array_equal(actual, ref[arm]+saved['display_scale']*disp[arm])
            np.testing.assert_array_equal(rotations[theory].get_ydata(), states[0, :, 2])
            np.testing.assert_array_equal(rotations[theory].get_xdata(), np.linspace(0, 1, 129))
        assert top.get_aspect() == 1
    assert [t.get_text() for t in fig.legends[0].get_texts()] == ['EB', 'RLB']
    plt.close(fig)


def test_plot_only_generation_and_source_preservation(saved, tmp_path, monkeypatch):
    from scripts.analysis.laminated_beams import check_inplane_spring_robustness as run
    def forbidden(*args, **kwargs):
        pytest.fail('no roots, matrices, forms or tracking in plot-only')
    for obj, name in ((run, 'Run'), (run, 'preflight'), (run, 'track_window'), (run.mechanics, 'recover'),
                      (run.eb, 'state_matrix'), (run.eb, 'transfer_matrix'),
                      (run.rlb, 'state_matrix'), (run.rlb, 'transfer_matrix')):
        monkeypatch.setattr(obj, name, forbidden)
    before = {p: plot.sha(p) for p in plot.source_paths()}
    original = plot.render
    monkeypatch.setattr(plot, 'render', lambda: original(output=tmp_path))
    monkeypatch.setattr(sys, 'argv', ['check_inplane_spring_robustness.py', 'plot-only', '--figure04'])
    run.main()
    assert before == {p: plot.sha(p) for p in plot.source_paths()}
    assert plot.csv_rows(tmp_path/'figure04_data.csv') == [
        {k: str(v) for k, v in row.items()} for row in plot.figure_rows(saved)]
    manifest = plot.read_json(tmp_path/'figure_manifest.json')
    assert manifest['reused_elastic_shapes'] == 4
    assert all(manifest[k] == 0 for k in ('new_roots', 'new_beta', 'new_d', 'matrix_calls',
                                         'shape_recoveries', 'tracking_calls', 'interpolation', 'solver_changes'))
    from PIL import Image
    with Image.open(tmp_path/'figure04_eb_rlb_shapes.png') as image:
        assert image.size == (3000, 2160) and abs(image.info['dpi'][0]-300) < .1
    pdf = (tmp_path/'figure04_eb_rlb_shapes.pdf').read_bytes()
    assert pdf.startswith(b'%PDF-') and b'/Subtype /Image' not in pdf


def test_figure_preset_cannot_trigger_compute(monkeypatch):
    from scripts.analysis.laminated_beams import check_inplane_spring_robustness as run
    monkeypatch.setattr(run, 'Run', lambda: pytest.fail('must reject before computation'))
    monkeypatch.setattr(sys, 'argv', ['check_inplane_spring_robustness.py', 'compute', '--figure04'])
    with pytest.raises(SystemExit) as error:
        run.main()
    assert error.value.code == 2
