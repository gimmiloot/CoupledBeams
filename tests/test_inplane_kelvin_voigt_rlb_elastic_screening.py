"""Sparse screening checks without repeated root searches or damped solves."""
import ast
import csv
import inspect
import numpy as np
import pytest
from scripts.analysis.laminated_beams import screen_inplane_kelvin_voigt_rlb_elastic as run


@pytest.fixture(scope='module')
def saved():
    path=run.OUTPUT/'diagnostics.json'
    if not path.exists(): pytest.skip('local screening artifacts absent; no automatic regeneration')
    return run.eb.read_json(path)


@pytest.fixture(scope='module')
def rows(): return run.eb.read_csv(run.OUTPUT/'rlb_elastic_screening.csv')


@pytest.fixture(scope='module')
def shapes():
    with np.load(run.OUTPUT/'rlb_elastic_shapes.npz',allow_pickle=False) as archive:
        return dict(archive)


def test_exact_sparse_scope_sorted_and_guard(saved,rows):
    assert run.BETAS==(0.,5.,45.,75.) and len(rows)==24
    for beta in run.BETAS:
        group=[r for r in rows if float(r['beta_deg'])==beta]
        assert [r['rlb_sorted_mode'] for r in group]==[f'rlb_sorted_{i:02d}' for i in range(1,7)]
        frequencies=[float(r['Omega']) for r in group]
        assert 0<frequencies[0] and all(a<b for a,b in zip(frequencies,frequencies[1:]))
        assert saved['groups'][str(beta)]['guard']['status']=='CONFIRMED'
        assert saved['groups'][str(beta)]['guard']['Omega']>frequencies[-1]
    assert all(r['model']=='RLB' and float(r['d_theta'])==0 and float(r['kappa_theta'])==1 for r in rows)
    assert saved['new_root_groups']==2 and saved['new_states']==10 and saved['reused_states']==14


def test_production_reduced_real_matrices_only(monkeypatch):
    p,_=run.eb.configuration(); context=run.ElasticContext(5.,p)
    def forbidden(*a,**k): pytest.fail('no full search or complex Newton')
    monkeypatch.setattr(context.full,'matrices',forbidden)
    monkeypatch.setattr(run.production,'solve_mode',forbidden)
    monkeypatch.setattr(run.kv,'correct',forbidden)
    assert context.config.d_theta==0 and run.production.route(context.config)=='SYMMETRY_REDUCED'
    for eta in (1,-1):
        B=context.matrix(.1,eta)
        assert B.shape==(3,3) and not np.iscomplexobj(B)
    with pytest.raises(ValueError): run.ElasticContext(90.,p)
    with pytest.raises(ValueError): context.matrix(.1+.01j,1)
    assert context.calls.full_B==0 and context.calls.half_B==2
    print('\nD21_TEST_MATRIX_CALLS',context.calls.snapshot())


def test_real_reduction_cannot_hide_imaginary_component():
    np.testing.assert_array_equal(run.real_array(np.array([1.+0j])),[1.])
    with pytest.raises(ValueError): run.real_array(np.array([1.+1e-30j]))


def test_structural_inactive_fields_and_predictors(rows,shapes):
    inactive=[r for r in rows if r['eta']=='-1']
    assert len(inactive)==12
    for r in inactive:
        y=shapes[r['shape_key']+'__states']
        np.testing.assert_array_equal(y[1],-y[0]@run.sd.F)
        assert y[0,-1,2]==y[1,-1,2]
        assert r['activity_status']=='EXACT_INACTIVE_BY_SYMMETRY'
        assert all(float(r[k])==0 for k in ('Delta_psi','P_joint','s_joint','G','a_slope_pred','zeta_slope_pred'))
    assert run.eb.classify(1,0.,identical_arms=True,root_confirmed=True)=='ACTIVE'


def test_mass_rotation_global_normalization_and_gates(saved,rows,shapes):
    p,_=run.eb.configuration(); arms=(run.kv.Arm.reduced('RLB',p),)*2
    _,w=run.kv.quadrature(); rotary=[]
    for r in rows:
        y=shapes[r['shape_key']+'__states']; assert not np.iscomplexobj(y)
        vector=run.kv.mass_vector(y,arms,w); mass=float(np.vdot(vector,vector))
        assert mass==pytest.approx(1.,abs=1e-12) and mass==pytest.approx(float(r['mass_M']),abs=1e-15)
        rotational=sum(a.L*a.J*np.dot(w,v[:,2]**2) for a,v in zip(arms,y));rotary.append(rotational)
        translational=sum(a.L*a.m*np.dot(w,v[:,0]**2+v[:,1]**2) for a,v in zip(arms,y))
        assert mass==pytest.approx(translational+rotational,abs=1e-12)
        assert float(r['physical_residual'])<=run.kv.CRITERIA['physical_residual']
        assert float(r['energy_residual'])<=run.kv.CRITERIA['energy_residual']
        assert float(r['next_sigma_ratio'])>=run.kv.CRITERIA['simple_sigma_separation']
        assert float(r['G'])==float(r['zeta_slope_pred'])
    assert max(rotary)>1e-3  # omission of J changes mass measurably


@pytest.mark.parametrize('factor',[.031,-17.,2+3j])
def test_scale_invariance(rows,shapes,factor):
    r=next(r for r in rows if r['eta']=='1');y=shapes[r['shape_key']+'__states']
    p,_=run.eb.configuration();arms=(run.kv.Arm.reduced('RLB',p),)*2;_,w=run.kv.quadrature()
    v=run.kv.mass_vector(y*factor,arms,w);mass=float(np.vdot(v,v).real)
    delta=(y[0,-1,2]-y[1,-1,2])*factor
    assert abs(delta)**2/mass==pytest.approx(float(r['P_joint']),rel=1e-13)
    values=run.eb.participation(delta,float(r['omega']),mass)
    for k,v in values.items(): assert v==pytest.approx(float(r[k]),rel=1e-13)


def test_reuse_preserves_source_roots_and_forms(rows,shapes):
    with np.load(run.eb.SOURCE/'shapes.npz',allow_pickle=False) as source:
        pool=run.eb.read_csv(run.eb.SOURCE/'verified_roots.csv')
        for row in rows:
            if row['reuse_status']!='REUSED_ROOT_AND_FORM': continue
            key=row['source_provenance'].split('#')[-1]
            old=next(r for r in pool if r['shape_key']==key)
            assert row['Omega']==old['Omega']
            for name in ('states','reactions'):
                # Reaction storage is flattened; physical values unchanged.
                np.testing.assert_array_equal(shapes[row['shape_key']+'__'+name].ravel(),source[key+'__'+name].ravel())


def test_match_by_fields_not_indices_or_frequency():
    left=[dict(common=np.eye(4)[i],symmetry_class=1 if i<2 else -1) for i in range(4)]
    right=[left[2],left[1],left[3],left[0]]
    result=run.match_pairs(left,right)
    assert [r[0] for r in result]==[3,1,0,2]
    assert all(r[3]=='CONFIRMED' for r in result)
    # Perfect similarity in the wrong class cannot override symmetry.
    other=[dict(common=left[0]['common'],symmetry_class=-1),dict(common=left[2]['common'],symmetry_class=1)]
    result=run.match_pairs([left[0],left[2]],other)
    assert [r[0] for r in result]==[1,0] and all(r[3]=='MATCH_AMBIGUOUS' for r in result)


def test_ambiguous_assignment_not_forced():
    left=[dict(common=x,symmetry_class=1) for x in np.eye(2)]
    right=[dict(common=x,symmetry_class=1) for x in np.array([[1,1],[1,-1]])/np.sqrt(2)]
    assert all(r[3]=='MATCH_AMBIGUOUS' for r in run.match_pairs(left,right))


def test_ambiguous_comparison_suppresses_scientific_differences(saved,shapes,monkeypatch):
    rows=[r for g in saved['groups'].values() for r in g['rows']]
    original=run.match_pairs
    monkeypatch.setattr(run,'match_pairs',lambda a,b:[(j,m,g,'MATCH_AMBIGUOUS') for j,m,g,s in original(a,b)])
    result=run.compare(rows,shapes)
    assert len(result)==24 and all(r['delta_G'] is r['delta_Omega'] is None for r in result)


def test_saved_matching_one_to_one_read_only_eb(saved,rows):
    compared=run.eb.read_csv(run.OUTPUT/'eb_rlb_matched_comparison.csv')
    references=run.eb.read_csv(run.eb.OUTPUT/'elastic_screening.csv')
    assert len(compared)==24 and saved['matching_ambiguous']==0
    for beta in run.BETAS:
        group=[r for r in compared if float(r['beta_deg'])==beta]
        assert len({r['eb_sorted_mode'] for r in group})==len({r['rlb_sorted_mode'] for r in group})==6
    for c in compared:
        e=next(r for r in references if r['beta_deg']==c['beta_deg'] and r['sorted_mode']==c['eb_sorted_mode'])
        r=next(r for r in rows if r['beta_deg']==c['beta_deg'] and r['rlb_sorted_mode']==c['rlb_sorted_mode'])
        assert e['symmetry_eta']==c['eta']==r['eta']
        assert e['Omega']==c['Omega_EB'] and r['Omega']==c['Omega_RLB']
        assert float(c['MAC_common'])>=run.CRITERIA['comparison_MAC']
        assert float(c['competing_margin'])>=run.CRITERIA['comparison_margin']
        assert float(c['delta_Omega'])==pytest.approx((float(r['Omega'])-float(e['Omega']))/float(e['Omega']))
        if c['eta']=='-1': assert c['delta_G']==c['G_ratio']=='' and c['comparison_status']=='INACTIVE_BOTH_THEORIES'
        else: assert float(c['delta_G'])==pytest.approx(float(c['G_RLB'])/float(c['G_EB'])-1)


def test_k12_finite_d_read_only_check(saved):
    c=saved['K12_control'];assert c['ACTIVE']['status']==c['INACTIVE']['status']=='PASS'
    assert c['ACTIVE']['relative_a_difference']<=1e-3 and c['ACTIVE']['relative_G_difference']<=1e-3
    assert c['INACTIVE']['predicted_G']==c['INACTIVE']['predicted_a']==0
    assert c['ACTIVE']['source_key']=='RLB_ACTIVE_d1'


def test_scope_and_immutable_physics(saved):
    for key in ('new_positive_d_complex_roots','complex_newton_calls','full_root_searches',
                'cross_beta_tracking','new_d','beta_outside_selection','asymmetric_damping','FEM'):
        assert saved[key]==0
    assert saved['calls']['full_B']==2 and saved['calls']['full_B_z']==0
    assert all(run.eb.sha(run.ROOT/p)==h for p,h in saved['protected_sources'].items())
    tree=ast.parse(inspect.getsource(run))
    called={ast.unparse(n.func) for n in ast.walk(tree) if isinstance(n,ast.Call)}
    assert not ({'production.solve_mode','kv.correct','eb.compute','arch.compute','mechanics.recover'} & called)
    assert all(g['calls']['B']+g['calls']['B_z']<=run.CRITERIA['max_point_matrices'] for g in saved['groups'].values())


def test_missing_only_no_scientific_calls_or_mutation(saved,monkeypatch):
    def forbidden(*a,**k): pytest.fail('missing-only must not compute')
    for obj,name in ((run,'ElasticContext'),(run,'search_missing'),(run,'shape_at'),(run,'compare'),
                     (run.production,'solve_mode'),(run.kv,'correct'),(run.eb,'compute')):
        monkeypatch.setattr(obj,name,forbidden)
    before={p.name:run.eb.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
    result=run.compute()
    assert result['missing_only'] and all(v==0 for k,v in result.items() if k!='missing_only')
    assert before=={p.name:run.eb.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
