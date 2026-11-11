"""Local block identities, archive consistency and crossing regression.

No full spectrum, new pole catalogues, Ritz or other geometries are evaluated.
The archive-specific tests require the accepted private result bundle.
"""
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.linalg import block_diag
from scripts.analysis import audit_rl_full_rr_reflection_intersection as d


@pytest.fixture(scope='module')
def audited():
    return d.stage_a()


@pytest.fixture
def engine(audited):
    return d.Blocks(audited[2][(6.1,'+')])


def test_all_archive_records_and_independent_sorting(audited):
    rows,checks,first,_,_=audited
    assert len(rows)==394 and len(first)==394
    assert checks['forms']==1773 and checks['spectral_or_basis_calls']==0
    assert checks['nine_root_unions_match']
    assert all(s['distinct_float_values']==197 for s in checks['summary'].values())
    assert checks['summary']['-']['status']=='BRANCH_VARIATION_CONFIRMED'


@pytest.mark.parametrize('beta',[0.,6.15,90.])
def test_restrictions_against_production_joint(engine,beta):
    # The physical J and explicit equations independently determine C+ and C-.
    f=np.diag(d.sections.REFLECTION_STATE);eye=np.eye(6)
    transform=np.block([[eye,eye],[f,-f]])
    left=d.row_transform(beta)@d.p.joint_matrix(engine.model,engine.case(beta))@transform
    np.testing.assert_allclose(left,block_diag(*d.restrictions(beta)),atol=3e-16,rtol=2e-15)
    assert np.linalg.det(d.row_transform(beta))==pytest.approx(1/64,rel=2e-15)


@pytest.mark.parametrize('beta,omega',[(0.,.06),(6.15,.077),(90.,.10)])
def test_clamp_elimination_and_six_to_three_block_equivalence(engine,beta,omega):
    with d.local_only() as calls:
        pair,s0,tl=engine.pair(beta,omega)
        pp=np.diag([1.,-1.,-1.]);eye=np.eye(3)
        vc=np.block([[eye,eye],[pp,-pp]])
        full_after_clamps=d.p.joint_matrix(engine.model,engine.case(beta))@block_diag(tl,tl)
        raw=d.row_transform(beta)@full_after_clamps@vc
        target=block_diag(*(c@tl for c in d.restrictions(beta)))
        np.testing.assert_allclose(raw,target,atol=np.linalg.norm(target)*4e-16,rtol=4e-15)
        assert np.linalg.det(vc)==pytest.approx(-8.)
        assert s0.shape==(6,3) and np.linalg.matrix_rank(s0)==3
        assert engine.regularity[-1]['balanced_backward_error']<64*np.finfo(float).eps
        assert engine.regularity[-1]['Z0_balanced_condition']<2
        assert not calls


def test_frequency_normalization_and_domain_guard(engine):
    for lam in (2.,3.924,7.):
        assert d.frequency(.005*lam**2,engine.reference)==pytest.approx(lam,rel=3e-16)
    with pytest.raises(ValueError,match='outside proved'):
        engine.transfer(engine.domain*1.001)


def test_local_solution_stability_and_first_root_certification(engine,audited):
    with d.local_only() as calls:
        rows,solution=d.local_roots(engine,audited[2])
    assert not calls and len(rows)==12
    assert solution['beta_deg']==pytest.approx(6.1300435048849,abs=1e-10)
    assert solution['Lambda']==pytest.approx(3.924236598658322,abs=1e-12)
    assert max(abs(v['beta_deg']-solution['beta_deg']) for v in solution['nested'])<1e-9
    assert [v['nullity'] for v in solution['block_checks']]==[1,1]
    assert [v['count'] for v in solution['count_certificates'][-2:]]==[0,2]


def test_restored_full_states_mass_and_orthogonality(engine,audited):
    solution=d.read(d.OUTPUT/'crossing_solution.json')
    b,w=solution['beta_deg'],solution['omega'];modes=[]
    with d.local_only() as calls:
        for parity in (1,-1):
            mode,_,_=d.restore(engine,b,w,parity)
            check=d.form_checks(engine,b,w,mode,parity)
            assert check['reflection_defect']==0 and check['mass']==pytest.approx(1.,abs=5e-15)
            for arm in mode['arms']:
                h=d.base.state_matrix(engine.arm,w)
                # Independent residual evaluation checks the complete six-field state.
                assert np.linalg.norm(arm['dy']-arm['y']@h.T)/np.linalg.norm(arm['dy'])<1e-12
            reflected=d.sections.reflect_mode(mode)
            assert d.saved.mass(engine.model,reflected)==pytest.approx(1.,abs=5e-15)
            modes.append(mode)
        assert d.saved.mass_inner(engine.model,*modes)==0
        assert not calls


def test_no_false_multiplicity_for_nearby_distinct_class_roots(engine):
    rows=d.table(d.OUTPUT/'block_roots.csv')
    with d.local_only():
        for beta in (6.1,6.2):
            pair=[r for r in rows if float(r['beta_deg'])==beta]
            assert float(pair[0]['omega'])!=float(pair[1]['omega'])
            for row in pair:
                index=0 if row['symmetry_class']=='+' else 1
                checks=engine.checks(beta,float(row['omega']))
                assert checks[index]['nullity']==1
                assert checks[1-index]['nullity']==0
        # A close frequency alone never authorizes reconstruction of an exact kernel.
        with pytest.raises(ValueError,match='one-dimensional'):
            d.restore(engine,6.15,.07701,1)


def test_full_solver_is_explicitly_forbidden(engine):
    with d.local_only() as calls:
        with pytest.raises(AssertionError,match='prohibited'):
            d.p.boundary(engine.model,engine.case(6.15),.077)
    assert len(calls)==1
