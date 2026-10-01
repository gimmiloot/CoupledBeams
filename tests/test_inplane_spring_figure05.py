"""Figure 5 source transformations and plot-only checks, no eigenproblems."""
import csv
import json
import math
import sys

import pytest

from scripts.lib import inplane_spring_figure05 as plot


@pytest.fixture(scope='module')
def rows():
    if not plot.SOURCE.exists():
        pytest.skip('local K23 CSV required; never compute missing inputs')
    return plot.load_rows()


def raw_rows():
    with plot.SOURCE.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def write_source(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def test_exact_cases_complex_fields_and_percent_formulas(rows):
    assert [(r['case_id'], r['beta_deg'], r['mode']) for r in rows] == list(plot.CASES)
    sources = {r['case_id']: r for r in raw_rows()}
    for r in rows:
        source = sources[r['case_id']]
        assert r['d_theta'] == '0.001'
        for key in ('Omega_d_EB', 'Omega_d_RLB', 'zeta_EB', 'zeta_RLB'):
            assert r[key] == source[key]
        assert r['delta_Omega_d_percent'] == pytest.approx(100*float(source['damped_delta_Omega']), abs=1e-12)
        assert r['delta_zeta_percent'] == pytest.approx(100*float(source['complex_delta_zeta']), abs=1e-12)
    assert rows[-1]['delta_Omega_d_percent'] < 0 < rows[-1]['delta_zeta_percent']


def test_elastic_predictors_and_R0_cannot_enter_figure(rows, tmp_path):
    raw = raw_rows()
    for row in raw:
        for key in ('Omega0_EB', 'Omega0_RLB', 'elastic_delta_Omega', 'G_EB', 'G_RLB', 'elastic_delta_G'):
            row[key] = 'NOT_A_COMPLEX_RESULT'
        if row['case_id'] == 'R0':
            row['d_theta'] = row['Omega_d_EB'] = row['zeta_RLB'] = 'NOT_SELECTED'
    file = tmp_path/'source.csv'; write_source(file, raw)
    selected = plot.load_rows(file)
    for original, transformed in zip(rows, selected, strict=True):
        assert {k: v for k, v in original.items() if k != 'source'} == {
            k: v for k, v in transformed.items() if k != 'source'}


def test_wrong_reference_d_is_rejected(rows, tmp_path):
    raw = raw_rows()
    next(r for r in raw if r['case_id'] == 'R3')['d_theta'] = '.005'
    file = tmp_path/'source.csv'; write_source(file, raw)
    with pytest.raises(AssertionError):
        plot.load_rows(file)


def test_panels_share_case_positions_and_keep_signed_values(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 11}):
        fig = plot.draw(rows)
        left, right = fig.axes
        assert left.get_shared_y_axes().joined(left, right)
        assert left.get_xlim() != right.get_xlim()
        for ax, key in zip(fig.axes, ('delta_Omega_d_percent', 'delta_zeta_percent'), strict=True):
            assert len(ax.patches) == 4
            for index, (bar, row) in enumerate(zip(ax.patches, rows, strict=True)):
                assert bar.get_width() == row[key]
                assert bar.get_x() == 0
                assert bar.get_y()+bar.get_height()/2 == pytest.approx(index)
            assert list(ax.lines[0].get_xdata()) == [0, 0]
        assert [t.get_text() for t in left.get_yticklabels()] == [
            r'$0^{\circ}\,/\,05$', r'$45^{\circ}\,/\,05$', r'$75^{\circ}\,/\,05$', r'$45^{\circ}\,/\,03$']
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for ax in fig.axes:
            for label in [*ax.texts, ax.xaxis.label, *ax.get_yticklabels()]:
                if not label.get_visible(): continue
                box = label.get_window_extent(renderer)
                assert box.x0 >= 0 and box.x1 <= fig.bbox.width
                assert box.y0 >= 0 and box.y1 <= fig.bbox.height
        plt.close(fig)


def test_existing_entry_generates_files_without_science(rows, tmp_path, monkeypatch):
    from scripts.analysis.laminated_beams import check_inplane_spring_robustness as run
    def forbidden(*args, **kwargs):
        pytest.fail('Figure 5 must only transform saved complex CSV rows')
    for obj, name in ((run, 'Run'), (run, 'preflight'), (run, 'track_window'), (run.mechanics, 'recover'),
                      (run.eb, 'state_matrix'), (run.eb, 'transfer_matrix'),
                      (run.rlb, 'state_matrix'), (run.rlb, 'transfer_matrix')):
        monkeypatch.setattr(obj, name, forbidden)
    before = plot.sha(plot.SOURCE)
    original = plot.render
    monkeypatch.setattr(plot, 'render', lambda: original(output=tmp_path))
    monkeypatch.setattr(sys, 'argv', ['check_inplane_spring_robustness.py', 'plot-only', '--figure05', '--figure05-view', 'corrections'])
    run.main()
    assert plot.sha(plot.SOURCE) == before
    with (tmp_path/'figure05_data.csv').open(encoding='utf-8', newline='') as stream:
        assert list(csv.DictReader(stream)) == [{k: str(v) for k, v in r.items()} for r in rows]
    manifest = json.loads((tmp_path/'figure05_manifest.json').read_text(encoding='utf-8'))
    assert manifest['reused_comparison_rows'] == 4
    assert all(manifest[k] == 0 for k in ('new_roots', 'new_forms', 'new_beta', 'new_d', 'matrix_calls',
                                         'shape_recoveries', 'interpolation', 'solver_changes', 'high_precision'))
    from PIL import Image
    with Image.open(tmp_path/(plot.STEM+'.png')) as image:
        assert image.size == (3000, 1350) and abs(image.info['dpi'][0]-300) < .1
    pdf = (tmp_path/(plot.STEM+'.pdf')).read_bytes()
    assert pdf.startswith(b'%PDF-') and b'/Subtype /Image' not in pdf


def test_figure_preset_rejects_scientific_command(monkeypatch):
    from scripts.analysis.laminated_beams import check_inplane_spring_robustness as run
    monkeypatch.setattr(run, 'Run', lambda: pytest.fail('must reject before calculation'))
    monkeypatch.setattr(sys, 'argv', ['check_inplane_spring_robustness.py', 'compute', '--figure05'])
    with pytest.raises(SystemExit) as error:
        run.main()
    assert error.value.code == 2


@pytest.fixture(scope='module')
def envelope_rows():
    if not plot.SOURCE.exists():
        pytest.skip('K23 data required; never regenerate scientific inputs')
    return plot.load_envelope_rows()


def test_envelope_scope_and_dimensional_alpha(envelope_rows):
    assert [(r['case_id'], r['beta_deg'], r['theory']) for r in envelope_rows] == [
        ('R2', 45, 'EB'), ('R2', 45, 'RLB'), ('R3', 75, 'EB'), ('R3', 75, 'RLB')]
    raw = {r['case_id']: r for r in raw_rows()}
    stored = json.loads((plot.SOURCE.parent/'diagnostics.json').read_text(encoding='utf-8'))['all_states']
    for row in envelope_rows:
        assert row['mode'] == '05' and row['d_theta'] == '0.001'
        source = raw[row['case_id']]; theory = row['theory']
        assert row['Omega_d'] == source['Omega_d_'+theory] and row['zeta'] == source['zeta_'+theory]
        accepted, = [r for r in stored if r['case_id'] == row['case_id'] and r['theory'] == theory]
        for name in ('alpha', 'a', 'omega_d'):
            assert row[name] == pytest.approx(accepted[name], rel=1e-13)
        assert row['alpha']*row['t_ref'] == pytest.approx(row['a'], rel=1e-14)


def test_exact_envelope_and_half_decay(envelope_rows):
    for row in envelope_rows:
        zeta = row['zeta']; n = row['N_half']
        assert plot.envelope(0, zeta) == 1
        assert plot.envelope(n, zeta) == pytest.approx(.5, abs=1e-15)
        assert plot.envelope(2*n, zeta) == pytest.approx(.25, abs=1e-15)
        time = 2*math.pi*n/row['omega_d']
        assert plot.envelope(n, zeta) == pytest.approx(math.exp(-row['alpha']*time), abs=1e-15)
    # A purely mathematical check makes accidental small-zeta substitution visible.
    assert plot.envelope(1, .6) == pytest.approx(math.exp(-1.5*math.pi))
    assert plot.envelope(1, .6) != pytest.approx(math.exp(-1.2*math.pi))


def test_envelope_panels_markers_and_styles(envelope_rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from scripts.lib.inplane_spring_figure04 import COLORS, STYLES
    assert plot.ENVELOPE_COLORS == COLORS and plot.ENVELOPE_STYLES == STYLES
    with plt.rc_context({'font.family': 'DejaVu Sans', 'font.size': 11}):
        fig = plot.draw_envelopes(envelope_rows)
        for ax, beta in zip(fig.axes, (45, 75), strict=True):
            assert ax.get_xlim() == (0, plot.PERIOD_LIMITS[beta])
            curves = {line.get_label(): line for line in ax.lines if line.get_label() in ('EB', 'RLB')}
            assert set(curves) == {'EB', 'RLB'}
            for row in (r for r in envelope_rows if r['beta_deg'] == beta):
                n, amplitude = curves[row['theory']].get_data()
                assert len(n) == len(amplitude) == 501
                assert list(amplitude) == [plot.envelope(float(x), row['zeta']) for x in n]
                assert row['N_half'] < ax.get_xlim()[1]
                assert any(len(line.get_xdata()) == 1 and float(line.get_xdata()[0]) == row['N_half']
                           and float(line.get_ydata()[0]) == .5 for line in ax.lines)
        fig.canvas.draw(); renderer = fig.canvas.get_renderer()
        for ax in fig.axes:
            lo, hi = ax.get_xlim()
            visible_ticks = [label for x, label in zip(ax.get_xticks(), ax.get_xticklabels(), strict=True)
                             if lo <= x <= hi]
            for label in [*ax.texts, ax.title, ax.xaxis.label, ax.yaxis.label, *visible_ticks]:
                if not label.get_visible() or not label.get_text(): continue
                box = label.get_window_extent(renderer)
                assert 0 <= box.x0 <= box.x1 <= fig.bbox.width
                assert 0 <= box.y0 <= box.y1 <= fig.bbox.height
        plt.close(fig)


def test_default_figure05_envelopes_are_plot_only(envelope_rows, tmp_path, monkeypatch):
    from scripts.analysis.laminated_beams import check_inplane_spring_robustness as run
    def forbidden(*args, **kwargs):
        pytest.fail('no scientific work or historical bar render in envelope mode')
    for obj, name in ((run, 'Run'), (run, 'preflight'), (run, 'track_window'), (run.mechanics, 'recover'),
                      (run.eb, 'state_matrix'), (run.eb, 'transfer_matrix'),
                      (run.rlb, 'state_matrix'), (run.rlb, 'transfer_matrix'), (plot, 'render')):
        monkeypatch.setattr(obj, name, forbidden)
    protected = list(plot.OUTPUT.glob('*'))+[plot.SOURCE, plot.SOURCE.parent/'run_manifest.json']
    hashes = {p: plot.sha(p) for p in protected if p.is_file()}
    original = plot.render_envelopes
    monkeypatch.setattr(plot, 'render_envelopes', lambda: original(output=tmp_path))
    monkeypatch.setattr(sys, 'argv', ['check_inplane_spring_robustness.py', 'plot-only', '--figure05'])
    run.main()
    assert all(plot.sha(p) == h for p, h in hashes.items())
    with (tmp_path/'figure05_data.csv').open(encoding='utf-8', newline='') as stream:
        assert list(csv.DictReader(stream)) == [{k: str(v) for k, v in row.items()} for row in envelope_rows]
    manifest = json.loads((tmp_path/'figure05_manifest.json').read_text(encoding='utf-8'))
    assert manifest['reused_comparison_rows'] == 2 and manifest['displayed_modal_states'] == 4
    assert all(manifest[k] == 0 for k in ('new_roots', 'new_forms', 'new_beta', 'new_d', 'new_tracking',
                                         'new_solver_work', 'matrix_calls', 'shape_recoveries',
                                         'interpolation', 'solver_changes', 'high_precision'))
    from PIL import Image
    with Image.open(tmp_path/(plot.ENVELOPE_STEM+'.png')) as image:
        assert image.size == (3000, 1350) and abs(image.info['dpi'][0]-300) < .1
    pdf = (tmp_path/(plot.ENVELOPE_STEM+'.pdf')).read_bytes()
    assert pdf.startswith(b'%PDF-') and b'/Subtype /Image' not in pdf


def test_unadorned_plot_only_keeps_existing_route(monkeypatch):
    from scripts.analysis.laminated_beams import check_inplane_spring_robustness as run
    calls = []
    monkeypatch.setattr(run, 'render', lambda: calls.append('old plot-only'))
    monkeypatch.setattr(run, 'Run', lambda: pytest.fail('no scientific work'))
    monkeypatch.setattr(sys, 'argv', ['check_inplane_spring_robustness.py', 'plot-only'])
    run.main()
    assert calls == ['old plot-only']
