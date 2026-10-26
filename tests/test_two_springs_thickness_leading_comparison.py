"""Normalization, signed effects, source-only access and saved-data guards."""
from types import SimpleNamespace
import json
import math
import numpy as np
import pytest
from scripts.lib import two_springs_thickness_leading_comparison as study


def test_frequency_error_denominator_and_quartic_normalization():
    r=study.frequency_comparison(3.,2.,81.,16.)
    assert r['F_full']==r['F_lim']==2.25
    assert r['e']==r['delta_percent']==0
    r=study.frequency_comparison(3.,2.,64.,16.)
    assert r['e']==pytest.approx(.125)
    assert r['delta_percent']==pytest.approx(12.5)


def test_signed_thickness_effect_and_percentage_points():
    r=study.sign_effect(dict(F_full=2.,F_lim=2.1),dict(F_full=2.5,F_lim=2.5))
    assert r['D_full']==pytest.approx(-.2)
    assert r['D_lim']==pytest.approx(-.16)
    assert r['Delta_D_pp']==pytest.approx(4.)
    assert study.sign_effect(dict(F_full=1,F_lim=1),dict(F_full=1,F_lim=1))['Delta_D_pp']==0


def test_full_calculation_guard_restores_previous_functions():
    original=study.physics.prefix
    with study.forbid_full_calculations():
        for module in (study.physics,study.physics.base):
            for name in ('prefix','mode','ritz'):
                with pytest.raises(RuntimeError,match='forbidden'):getattr(module,name)()
        with pytest.raises(RuntimeError,match='forbidden'):study.api().Run.solve(None)
    assert study.physics.prefix is original


@pytest.mark.parametrize('operation',['plot_only','verify_cache'])
def test_saved_modes_do_not_initialize_leading_calculation(operation,monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('Unexpected matrix or spectral calculation')
    monkeypatch.setattr(study,'Experiment',forbidden)
    monkeypatch.setattr(study.physics.Leading,'solve',forbidden)
    events=[]
    monkeypatch.setattr(study,'verify_bundle',lambda folder,figures=True:events.append(('verify',figures)))
    monkeypatch.setattr(study,'render',lambda folder:events.append(('render',)))
    args=SimpleNamespace(plot_only=None,verify_cache=None)
    setattr(args,operation,study.BASE/'87cd1ba35fa1b99e')
    study.main(args)
    assert events==([('verify',False),('render',),('verify',True)] if operation=='plot_only' else [('verify',True)])


def test_saved_actual_nodes_and_control_coverage():
    folder=study.BASE/'87cd1ba35fa1b99e'
    if not (folder/'summary.json').exists():pytest.skip('Local comparison arrays absent')
    rows=study.read_csv(folder/'comparison.csv');source=study.read_csv(folder/'reused_full_points.csv')
    assert len(rows)==len(source)==232
    keys=lambda r:(r['pair'],float(r['mu']),float(r['tau']),float(r['beta_deg']))
    assert {keys(r) for r in rows}=={keys(r) for r in source}
    for r in rows:
        assert float(r['xi'])==pytest.approx(math.radians(float(r['beta_deg']))/.005,abs=1e-14)
        assert 0<=float(r['xi'])<=30
        assert float(r['delta_percent'])==pytest.approx(100*abs(float(r['F_full'])/float(r['F_lim'])-1),abs=1e-13)
    controls=study.read_csv(folder/'convergence.csv');controlled={keys(r) for r in controls}
    for stat in study.read_csv(folder/'error_statistics.csv'):
        assert (stat['pair'],float(stat['mu']),float(stat['tau']),float(stat['max_beta_deg'])) in controlled
    assert all(r['status']=='PASS' for r in study.read_csv(folder/'leading_quality.csv'))
    assert json.loads((folder/'summary.json').read_text())['full_model_solves']==0
