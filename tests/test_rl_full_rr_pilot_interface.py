"""K6/K8 and two-thickness interface regressions; no physical root searches."""
import copy
import csv
import json
import math
from pathlib import Path

import numpy as np
import pytest

from scripts.analysis import map_rl_full_rr as m
from scripts.lib import rl_full_rr_sections_prefix as e
from scripts.lib.rl_full_rr_aspect import make_model
from scripts.lib.rl_full_rr_square_checks import check_coefficients


@pytest.fixture
def policy():
    return m.configuration(m.CONFIG)[1]


def contract(policy, *, k=6, epsilon=.005, mu=.1, tau=0, beta=15):
    model=make_model(aspect=1.,epsilon=epsilon,mu=mu,tau=tau)
    return e.point_contract(model,e.p.Case(beta,mu,'RIGID','RIGID',tau),k,k+1,policy)


@pytest.mark.parametrize('k,guard,accepted',[(6,7,True),(8,9,True),(10,11,False),(6,9,False),(8,7,False),(6.,7,False)])
def test_only_authorized_prefix_pairs(k,guard,accepted):
    if accepted:assert m.prefix_ranks(dict(K_plot=k,K_guard=guard))==(k,guard)
    else:
        with pytest.raises(ValueError):m.prefix_ranks(dict(K_plot=k,K_guard=guard))
    assert m.prefix_ranks({})==(8,9)


def test_reuse_direction_and_full_physics(policy):
    upper=contract(policy,k=8);lower=contract(policy)
    assert m.reusable_prefix(upper,lower)
    assert not m.reusable_prefix(lower,upper)
    for key,value in [('epsilon',.01),('mu',.2),('tau',.1),('beta',30)]:
        assert not m.reusable_prefix(upper,contract(policy,**{key:value}))
    for parent,key,value in [('model','ks',.8),('model','nu',.2),('model','rho',2.),
                              ('model','E',2.),('model','variant','ELEMENTARY_LEGACY'),
                              ('case','k_joint',0.)]:
        other=copy.deepcopy(lower);other[parent][key]=value
        assert not m.reusable_prefix(upper,other)
    for key,value in [('normalization','omega=Lambda'),('time_reference',100.),('coordinates','reversed')]:
        other=copy.deepcopy(lower);other[key]=value
        assert not m.reusable_prefix(upper,other)
    other=copy.deepcopy(lower);other['arms'][0]['J']*=2
    assert not m.reusable_prefix(upper,other)


@pytest.mark.parametrize('epsilon',[.005,.01])
@pytest.mark.parametrize('tau',[-.5,0,.1,.5])
def test_real_coefficients_volume_and_regular_domain(epsilon,tau):
    model=make_model(aspect=1.,epsilon=epsilon,mu=.3,tau=tau)
    case=e.p.Case(15.,.3,'RIGID','RIGID',tau)
    checked=check_coefficients(model,case)
    assert checked['status']=='PASS'
    assert sum(a.section.area*L for a,L in zip(model.arms,case.lengths(model)))==pytest.approx(24*epsilon**2)
    for arm in model.arms:
        s=arm.section
        assert s.width==s.thickness
        assert s.area==pytest.approx(s.thickness**2)
        assert s.second_moment==pytest.approx(s.thickness**4/12)
        assert arm.J==pytest.approx(2*.3**2*s.rhoI)
    assert model.time==pytest.approx(1/epsilon)
    assert model.omega(7.)==pytest.approx(epsilon*49)


def test_epsilon_changes_stiffness_inertia_and_basis_domain():
    a,b=[make_model(aspect=1.,epsilon=x,mu=.3,tau=.1) for x in (.005,.01)]
    for x,y in zip(a.arms,b.arms):
        for field,ratio in [('EA',4),('KGA',4),('rhoA',4),('EI',16),('rhoI',16)]:
            assert getattr(y.section,field)/getattr(x.section,field)==pytest.approx(ratio)
        assert y.J/x.J==pytest.approx(16)
        assert math.sqrt(y.section.EA/y.J)/math.sqrt(x.section.EA/x.J)==pytest.approx(.5)
        assert math.sqrt(y.section.KGA/y.section.rhoI)/math.sqrt(x.section.KGA/x.section.rhoI)==pytest.approx(.5)


def archived_point():
    folder=m.ROOT/'results/rl_full_rr_maps/b775d0980d7bad77'
    manifest=m.verify_manifest(folder)
    item=next(p for p in manifest['points'] if p['beta_deg']==15 and p['mu']==.1 and p['tau']==0)
    path=m.ROOT/item['record'];assert m.sha(path)==item['record_sha256']
    return path,m.verify_record(path)


def test_saved_K8_supplies_seven_complete_normalized_forms(tmp_path,monkeypatch):
    path,old=archived_point();target=dict(old['contract'],K_plot=6,K_guard=7)
    hashes={p:m.sha(p) for p in (path,path.parent.parent/old['shape_file'])}
    def forbidden(*a,**kw):raise AssertionError('Spectral work during reuse')
    for obj,name in [(m,'solve'),(e,'catalogues'),(e,'counted_prefix'),(e.p,'boundary'),(e.p,'ritz'),(e.p.base,'basis')]:
        monkeypatch.setattr(obj,name,forbidden)
    derived,record=m.prefix_view(path,target,tmp_path,'test')
    assert len(record['roots'])==7
    assert [r['Lambda'] for r in record['roots']]==[r['Lambda'] for r in old['roots'][:7]]
    assert [r['role'] for r in record['roots']]==['plot']*6+['guard']
    assert record['quality']['status']=='PASS'
    assert record['search']['guard_right_count']==7
    assert record['search']['new_count_queries']==0
    assert record['reuse']['source_record_sha256']==hashes[path]
    with np.load(tmp_path/record['shape_file']) as z, np.load(path.parent.parent/old['shape_file']) as original:
        assert len(z.files)==70
        for key in z.files:np.testing.assert_array_equal(z[key],original[key])
    assert hashes=={p:m.sha(p) for p in hashes}
    assert m.prefix_view(path,target,tmp_path,'test')[0]==derived


def test_failure_records_use_requested_prefix_without_retry(tmp_path,monkeypatch,policy):
    seen=[]
    def fail(*a,**kw):seen.append(True);raise e.PrefixFailure('COUNT_OUTSIDE_CATALOG',domain=.1)
    monkeypatch.setattr(e,'catalogues',fail)
    args=(contract(policy,epsilon=.01),tmp_path,'test',1.6)
    path,rec,new=m.solve(*args)
    assert new and rec['status']=='COUNT_OUTSIDE_CATALOG'
    assert rec['contract']['K_guard']==7
    assert not m.solve(*args)[2] and len(seen)==1


def test_cache_different_epsilon_rejected_even_same_engine(policy):
    old=dict(engine_id='same',versions={},policy=policy,initial_upper_omega=1.6,material=dict(epsilon=.005),sources={})
    new=dict(old,material=dict(epsilon=.01))
    assert not m.cache_compatible(old,new)


@pytest.mark.parametrize('separation',[0.,1e-7])
def test_K6_guard_cluster_kept_complete_without_merging_close_roots(policy,separation):
    roots=np.array([1.,2.,3.,4.,5.,6.,7.,7.+separation,10.])
    events,search=e.counted_prefix(lambda w:np.diag(np.r_[roots-w,1.]),
        lambda w:(w,int(sum(roots<w)),{}),.1,11.,7,policy)
    assert search['computed_ranks']==search['guard_right_count']==8
    assert events[-1]['last_rank']==8
    assert events[-1]['multiplicity']==(2 if separation==0 else 1)
    if separation:
        assert events[-1]['omega']>events[-2]['omega']


def test_K6_plot_hides_guard_and_never_solves(tmp_path,monkeypatch):
    folder=tmp_path/'source';folder.mkdir()
    rows=[dict(map='pilot',rank=k,x=x,Lambda=k+x,status='PASS') for x in (.8,.9) for k in range(1,8)]
    with (folder/'spectrum.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    m.publish(folder/'manifest.json',dict(identity=dict(K_plot=6,K_guard=7,material=dict(epsilon=.01,aspect=1),
        maps=[dict(name='pilot',mode='MU',fixed=dict(beta_deg=45,tau=0),grid=[.8,.9])]),artifacts={'spectrum.csv':m.sha(folder/'spectrum.csv')}))
    monkeypatch.setattr(m,'OUTPUT',tmp_path/'render')
    def forbidden(*a,**kw):raise AssertionError('Solver called by plot-only')
    monkeypatch.setattr(m,'solve',forbidden);monkeypatch.setattr(e,'catalogues',forbidden)
    from matplotlib.axes import Axes
    original=Axes.plot;calls=[]
    def plot(self,x,y,*args,**kwargs):calls.append(list(y));return original(self,x,y,*args,**kwargs)
    monkeypatch.setattr(Axes,'plot',plot)
    m.render(folder,dict(layout='all',dpi=30))
    assert len(calls)==6 and max(max(v) for v in calls)<7
