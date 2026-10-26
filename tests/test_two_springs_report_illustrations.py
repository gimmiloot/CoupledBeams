"""Archive-backed safeguards for report normalization and zero-solve rendering."""
import copy
import inspect

import numpy as np
import pytest

from scripts.analysis import prepare_two_springs_report_illustrations as report


@pytest.fixture
def saved():
    if not (report.BASE/report.SOURCES['growth'][0]/'manifest.json').exists():
        pytest.skip('Local scientific archives are not included in Git')
    return report.SavedData()


def test_report_archive_controls_and_render_without_scientific_calls(saved, tmp_path, monkeypatch):
    from scripts.lib import isotropic_two_springs_leading as leading
    def forbidden(*args, **kwargs):
        raise AssertionError('Scientific computation called during plot_only')
    for name, function in inspect.getmembers(report.archive.physics, inspect.isfunction):
        if function.__module__ == report.archive.physics.__name__:
            monkeypatch.setattr(report.archive.physics, name, forbidden)
    for name in ('mechanics_fields', 'variation_analysis', 'leading_postprocess', 'brentq'):
        monkeypatch.setattr(report.archive, name, forbidden)
    monkeypatch.setattr(leading.Ritz, '__init__', forbidden)
    monkeypatch.setattr(leading.Ritz, 'solve', forbidden)
    monkeypatch.setattr(np.linalg, 'svd', forbidden)
    tables, checks, inputs = report.prepare()
    assert checks['status'] == 'PASS'
    assert checks['overlap_prefix_comparisons'] == 234
    assert checks['maximum_overlap_relative_difference'] == 0
    assert [r['beta_pl_deg'] for r in checks['figure1']] == [6, 12]
    assert checks['figure4']['displacement_reduction'] == pytest.approx(8.135178352495185)
    report.render(tables, tmp_path)
    from PIL import Image
    for name in report.FIGURES:
        assert (tmp_path/(name+'.pdf')).read_bytes().startswith(b'%PDF')
        with Image.open(tmp_path/(name+'.png')) as image:
            assert min(image.info['dpi']) >= 300
            assert min(image.size) >= 1500
    inputs.unchanged()


def test_report_rejects_wrong_leading_error_denominator(saved, monkeypatch):
    original = saved.read
    def changed(name, filename):
        rows = copy.deepcopy(original(name, filename))
        if (name, filename) == ('leading', 'comparison.csv'):
            r = next(r for r in rows if r['ID'] == report.PAIR and float(r['xi']) > 20)
            # Plausible-looking but wrong normalization must not be plotted.
            r['E'] = str(1-float(r['F_lim'])/float(r['F_full']))
        return rows
    monkeypatch.setattr(saved, 'read', changed)
    with pytest.raises(ValueError, match='denominator'):
        report.leading_data(saved)


def test_report_checks_tail_outside_visible_interval(saved, monkeypatch):
    original = report.archive.stiffness_slopes
    def changed(beta, values):
        slope, names = original(beta, values)
        slope[-1] = 1  # Hidden beta=90 point; the plotted interval ends at 20.
        return slope, names
    monkeypatch.setattr(report.archive, 'stiffness_slopes', changed)
    with pytest.raises(ValueError, match='identity|derivative'):
        report.thickness_data(saved)
