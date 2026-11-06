"""Map orchestration and saved-data contract; no new spectral benchmarks."""
import copy
import csv
import json
from pathlib import Path
import numpy as np
import pytest
from scripts.analysis import map_rl_full_rr as m
from scripts.lib import rl_full_rr_sections_prefix as e


def contract(mu=.01, tau=0., beta=5.):
    c, policy = m.configuration(m.CONFIG)
    return e.point_contract(e.p.Model(**c['material'], mu=mu, tau=tau),
        e.p.Case(beta, mu, 'RIGID', 'RIGID', tau), 8, 9, policy)


def test_reuse_contract_numeric_equivalence_and_separation():
    assert m.physical_key(contract(beta=5)) == m.physical_key(contract(beta=5.))
    keys = [m.physical_key(contract()), m.physical_key(contract(mu=.02)),
            m.physical_key(contract(tau=.05)), m.physical_key(contract(beta=30))]
    assert len(set(keys)) == 4
    changed = contract(); changed['policy']['root_rtol'] *= 2
    assert m.physical_key(changed) != keys[0]
    assert m.digest(dict(mode='fast_plot', key=keys[0])) != m.digest(dict(mode='certified_audit', key=keys[0]))


def test_accepted_reuse_all_thirteen_demo_entries():
    c, _ = m.configuration(m.CONFIG)
    known, _, _ = m.accepted_sources(c)
    counts = []
    for spec in c['maps']:
        counts.append(0)
        for x in spec['grid']:
            q = dict(spec['fixed'], **{m.AXES[spec['mode']]:x})
            counts[-1] += m.physical_key(contract(q['mu'],q['tau'],q['beta_deg'])) in known
    assert counts == [6,2,2,1,2]


def test_geometry_and_lambda_common_reference():
    c, _ = m.configuration(m.CONFIG)
    for mu,tau in [(0.,0.),(.1,-.1),(.1,.1)]:
        model = e.p.Model(**c['material'],mu=mu,tau=tau)
        case = e.p.Case(5,mu,'RIGID','RIGID',tau)
        g = m.geometry(model,case)
        assert g['volume_relative_error'] < 1e-14
        assert model.time == pytest.approx(200)
        for a in g['arms']:
            assert a['b']/a['h'] == 4
            assert a['J'] == pytest.approx(17*.3**2*a['rhoI'])
        w = .5
        lam = model.Lambda(w)
        assert lam == pytest.approx(10)
        assert lam**4 == pytest.approx(model.section.rhoA*w*w/model.section.EI)


def test_validation_no_unordered_or_epsilon_sweep(tmp_path):
    c = m.read(m.CONFIG)
    c['maps'][0]['grid'] = [0,5,4]
    path = tmp_path/'config.json'; path.write_text(json.dumps(c))
    with pytest.raises(ValueError): m.configuration(path)
    c = m.read(m.CONFIG); c['material']['epsilon'] = .006
    path.write_text(json.dumps(c))
    with pytest.raises(ValueError): m.configuration(path)


def test_close_distinct_and_multiplicity_agreement():
    a = dict(roots=[dict(omega=1.,multiplicity=1),dict(omega=1.+1e-8,multiplicity=1)])
    b = copy.deepcopy(a)
    assert m.agreement(a,b,1e-8)
    b['roots'][1]['multiplicity'] = 2
    assert not m.agreement(a,b,1e-8)
    b = copy.deepcopy(a); b['roots'][0]['omega'] *= 1.25
    assert m.jump(a,b) == pytest.approx(.25)


def saved_fixture(tmp_path):
    folder = tmp_path/'bundle'; folder.mkdir()
    rows=[]
    for x,status in [(0,'PASS'),(1,'FAILED'),(2,'PASS')]:
        for k in range(1,10):
            rows.append(dict(map='example',rank=k,x=x,Lambda=k+x if status=='PASS' else '',status=status))
    with (folder/'spectrum.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    m.publish(folder/'manifest.json',dict(artifacts={'spectrum.csv':m.sha(folder/'spectrum.csv')},
        identity=dict(maps=[dict(name='example',mode='BETA',fixed=dict(mu=.01,tau=0))])))
    return folder


def test_plot_tables_only_zero_solves_no_guard_no_gap_fill(tmp_path,monkeypatch):
    folder = saved_fixture(tmp_path)
    before={p.name:m.sha(p) for p in folder.iterdir()}
    def forbidden(*a,**k): raise AssertionError('Physical call in saved-data mode')
    for obj,name in [(m,'solve'),(e,'counted_prefix'),(e,'catalogues'),(e.p,'ritz'),
                     (e.p,'boundary'),(e.p,'count_at'),(e.p,'mode'),(e.p.base,'basis')]:
        monkeypatch.setattr(obj,name,forbidden)
    monkeypatch.setattr(m,'OUTPUT',tmp_path/'output')
    assert m.render(folder,{},True)['root_calculations']==0
    from matplotlib.axes import Axes
    original=Axes.plot; calls=[]
    def capture(self,x,y,*a,**k):
        calls.append((x,y,k)); return original(self,x,y,*a,**k)
    monkeypatch.setattr(Axes,'plot',capture)
    out=m.render(folder,dict(layout='all',dpi=60))
    assert len(calls)==8
    assert all(np.isnan(y[1]) for _,y,_ in calls)
    assert m.read(out/'manifest.json')['root_calculations']==0
    assert len(list(out.glob('*.pdf')))==1
    assert before=={p.name:m.sha(p) for p in folder.iterdir()}


def test_tampered_table_refused(tmp_path):
    folder=saved_fixture(tmp_path)
    with (folder/'spectrum.csv').open('a') as f:f.write('tampered')
    with pytest.raises(ValueError,match='checksum'):m.saved_table(folder)


def test_publication_does_not_overwrite(tmp_path):
    p=tmp_path/'record.json'; m.publish(p,dict(status='PASS'))
    before=m.sha(p)
    with pytest.raises(FileExistsError):m.publish(p,dict(status='FAILED'))
    assert m.sha(p)==before


def test_numerical_failure_is_sealed_and_resumed_without_solve(tmp_path,monkeypatch):
    def fail(*a,**k):
        raise e.PrefixFailure('COUNT_OUTSIDE_CATALOG',requested=2.,domain=1.)
    monkeypatch.setattr(e,'catalogues',fail)
    path,rec,calculated=m.solve(contract(),tmp_path,'test-source',1.6)
    assert calculated and rec['status']=='COUNT_OUTSIDE_CATALOG'
    assert rec['failure']==dict(requested=2.,domain=1.)
    before=m.sha(path)
    def forbidden(*a,**k):raise AssertionError('Recomputed a saved failure')
    monkeypatch.setattr(e,'catalogues',forbidden)
    _,again,calculated=m.solve(contract(),tmp_path,'test-source',1.6)
    assert not calculated and again==rec and m.sha(path)==before


def test_v1_cache_migration_requires_unchanged_scientific_contract():
    old=m.read(m.ROOT/'results/rl_full_rr_maps/58bfb20de7bbd76e/manifest.json')['identity']
    current=copy.deepcopy(old);current['engine_id']='changed-orchestration-only'
    current['sources']['scripts/analysis/map_rl_full_rr.py']=m.sha(m.__file__)
    assert m.cache_compatible(old,current)
    current['policy']['energy_tol']*=2
    assert not m.cache_compatible(old,current)
