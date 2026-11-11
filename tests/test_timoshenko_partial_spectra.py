"""Boundary adaptation and regression at the explicitly requested square lengths."""
import json
from pathlib import Path
import numpy as np
import pytest
from scripts.lib import timoshenko_partial_spectra as part
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib.rl_full_rr_aspect import make_model

ROOT=Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def setup():
    model=make_model(epsilon=.005,aspect=1.,E=1.,rho=1.,l=1.,nu=.3,ks=5/6,variant='RL_FULL').arms[0].carrier
    policy=json.loads((ROOT/'data/input/rl_full_rr_geometry_c.json').read_text())['policy']
    inv=json.loads((ROOT/'private/dissertation/_history/20261110_section142_dual_sources/postprocess/cc_catalog_inventory.json').read_text())
    def catalog(L):
        item=next(i for i in inv if abs(i['length']-L)<1e-12)
        point=json.loads((ROOT/item['record']).read_text())
        return next(c for c in point['catalogs'].values() if c['catalog_id']==item['catalog_id'])
    return model,policy,catalog


def test_requested_ode_exact_coefficients(setup):
    model,_,_=setup;p=model.coefficients;w=.2
    expected=np.array([[0,1,1/p['S'],0],[0,0,0,1/p['B']],[-p['m']*w**2,0,0,0],[0,-p['r']*w**2,-1,0]])
    np.testing.assert_array_equal(mh.harmonic_state_matrix(model,w,'timoshenko'),expected)


def test_private_adapter_leaves_cc_function_unchanged(setup):
    model,_,_=setup
    assert mh.finite_boundary_matrix is not part.cs_boundary_matrix
    assert part._cs.finite_boundary_matrix is part.cs_boundary_matrix
    cc=mh.finite_boundary_matrix(model,1.,.2,'timoshenko')
    np.testing.assert_array_equal(part.boundary_matrix(model,1.,.2,'CC'),cc)
    assert not np.allclose(cc,part.cs_boundary_matrix(model,1.,.2))


def test_cs_rejects_other_physics(setup):
    with pytest.raises(ValueError):part.cs_boundary_matrix(setup[0],1.,.2,'mh')


@pytest.mark.parametrize('length',[.1,.3,.65,1.,1.35,1.9])
def test_count_normalization_and_independent_boundaries(setup,length):
    model,policy,catalog=setup;cc=catalog(length)
    roots,cert=part.solve_cs(model,length,cc,.605,policy)
    modes,quality=part.verify_modes(model,length,roots,'CS',policy)
    assert len(roots)==cert['upper_count'] and quality['status']=='PASS'
    for i,(root,mode) in enumerate(zip(roots,modes)):
        assert root['omega']<=cc['poles'][i]
        assert abs(mode['ends'][1,3])<1e-9*np.max(abs(mode['values'][:,3]))
        # Rotation is fixed; the Timoshenko transverse slope is not fixed.
        slope=mh.finite_state_basis(model,length,root['omega'],0.,'timoshenko',1)@mode['coefficients']
        assert abs(slope[0])>1e-8
    # Existing CC roots checked without invoking any CC root finder.
    original=[r for r in cc['roots'] if r['omega']<=.605]
    _,qcc=part.verify_modes(model,length,original,'CC',policy)
    assert qcc['status']=='PASS'


def test_missing_count_is_rejected(setup,monkeypatch):
    model,policy,catalog=setup
    monkeypatch.setattr(part._cs,'finite_roots',lambda *a,**k:([],{}))
    with pytest.raises(ArithmeticError,match='not saturated'):
        part.solve_cs(model,1.,catalog(1.),.605,policy)


def test_short_length_needs_no_root_search(setup,monkeypatch):
    def forbidden(*a,**k):raise AssertionError('unexpected search')
    monkeypatch.setattr(part._cs,'finite_roots',forbidden)
    roots,cert=part.solve_cs(setup[0],.1,setup[2](.1),.605,setup[1])
    assert not roots and cert['search_calls']==0


@pytest.fixture(scope='module')
def completed():
    path=ROOT/'results/rl_full_rr_maps/_section142_partials_20261110'
    if not (path/'manifest.json').exists():pytest.skip('local completed partial bundle absent')
    return path


def test_completed_band_and_equal_lengths(completed):
    from scripts.analysis import map_timoshenko_partial_frequencies as run
    data=run.rows(completed/'partial_frequencies.csv')
    assert sum(r['bc']=='CS' and r['role']=='plot' for r in data)==1003
    assert sum(r['bc']=='CC' and r['role']=='plot' for r in data)==917
    for r in data:
        assert float(r['omega'])==pytest.approx(.005*float(r['Lambda'])**2,rel=5e-15)
        assert (r['role']=='plot')==(float(r['Lambda'])<=11)
    rows=[r for r in run.rows(completed/'arm_frequencies.csv') if float(r['mu'])==0]
    first={(r['bc'],r['rank']):r['omega'] for r in rows if r['arm']=='1'}
    second={(r['bc'],r['rank']):r['omega'] for r in rows if r['arm']=='2'}
    assert first==second


def test_plot_only_forbids_spectral_calls_and_reproduces(completed,tmp_path,monkeypatch):
    from scripts.analysis import map_timoshenko_partial_frequencies as run
    def forbidden(*a,**k):raise AssertionError('solver/basis called by plot-only')
    for owner,name in [(part,'solve_cs'),(part._cs,'finite_roots'),(mh,'finite_roots'),(mh,'finite_state_basis')]:
        monkeypatch.setattr(owner,name,forbidden)
    run.render(completed,tmp_path/'redraw')
    for name in ('frequency_A_CS','frequency_B_CS_CC'):
        for extension in ('pdf','png'):
            filename=name+'.'+extension
            assert run.sha(tmp_path/'redraw'/filename)==run.sha(completed/'figures'/filename)


def test_plot_rejects_modified_input_table(completed,tmp_path):
    import shutil
    from scripts.analysis import map_timoshenko_partial_frequencies as run
    for name in ('manifest.json','contract.json','partial_frequencies.csv','coupled_frequencies.csv','arm_frequencies.csv'):
        shutil.copy2(completed/name,tmp_path/name)
    with (tmp_path/'partial_frequencies.csv').open('a') as f:f.write('tampered\n')
    with pytest.raises(AssertionError):run.render(tmp_path,tmp_path/'bad_plot')


def test_other_thickness_cannot_reuse_this_scope(tmp_path,monkeypatch):
    from scripts.analysis import map_timoshenko_partial_frequencies as run
    def forbidden(*a,**k):raise AssertionError('search must not start for another physical scope')
    monkeypatch.setattr(part,'solve_cs',forbidden)
    config=run.read(run.DEFAULT);config['material']['epsilon']=.01
    path=tmp_path/'other.json';run.write(path,config)
    with pytest.raises(AssertionError):run.compute(path,tmp_path/'other_bundle')
