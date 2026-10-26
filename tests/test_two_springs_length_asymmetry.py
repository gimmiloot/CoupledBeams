"""Addressed checks for the bounded length study, without spectral solves."""
from fractions import Fraction
import json
from types import SimpleNamespace

import numpy as np
import pytest

from scripts.lib import two_springs_length_asymmetry as study
from scripts.lib import isotropic_two_springs_rl_timoshenko as physics


def grid():
    cfg = json.loads(study.CONFIG.read_text())
    return cfg['map']['parameter_grid']+cfg['tail_grid']


def test_negative_mu_preserves_arm_supports_and_coordinates():
    model = physics.Model()
    for pair in ((.1,10),(10,.1)):
        pos,neg = [physics.Case(6,mu,*pair) for mu in (.1,-.1)]
        assert neg.lengths(model)==(1.1,.9)
        assert pos.lengths(model)==(.9,1.1)
        for arm in (0,1):
            np.testing.assert_array_equal(physics.external_matrix(model,pos,arm),
                                          physics.external_matrix(model,neg,arm))
        np.testing.assert_array_equal(physics.joint_matrix(model,pos),physics.joint_matrix(model,neg))
        np.testing.assert_array_equal(physics.external_matrix(model,neg,0),np.eye(6)[:3])
        row=physics.external_matrix(model,neg,1)[2]
        assert row[2]==pair[1]*model.section.EI/model.l and row[5]==-1


def test_axial_minimum_fixed_sum_is_independent_of_mu():
    # Exact rational verification of the constrained quadratic minimum.
    for mu in map(Fraction,('-0.1','-0.01','0','0.01','0.1')):
        L1,L2=1-mu,1+mu
        EA,defect=Fraction(7),Fraction(3,10)
        d1,d2=defect*L1/(L1+L2),defect*L2/(L1+L2)
        assert d1+d2==defect
        assert EA*d1/L1==EA*d2/L2
        assert EA*(d1*d1/L1+d2*d2/L2)/2==EA*defect**2/4


def test_nonuniform_join_stencils_and_zero_symmetry():
    beta=np.array(grid())
    slope,names=study.length_slopes(beta,beta**2)
    np.testing.assert_allclose(slope,2*beta,atol=1e-11,rtol=0)
    assert names[0]=='reflection_symmetry'
    i=list(beta).index(20)
    assert 'nonuniform' in names[i]
    # A large unused intermediate ordinate cannot affect the archived beta=3 join.
    y=beta**2
    y[list(beta).index(2.5)]+=100
    changed,_=study.length_slopes(beta,y)
    assert changed[list(beta).index(3)]==6


def test_late_rebound_disqualifies_early_crossing():
    beta=np.array(grid())
    q=np.where(beta<6,.2,.1)
    q[list(beta).index(80)]=.3
    pick,_=study.api().refined_pick(beta,q,.15,30)
    # A late rebound leaves less than the required 30-degree qualifying tail.
    assert pick['beta_pl_deg'] is None


def test_archive_identities_signed_ratios_and_leading_acceptance():
    folder=study.BASE/'length_asymmetry/ddf22cb06e8563e2'
    if not (folder/'numerical_summary.json').exists():
        pytest.skip('Local research evidence is not distributed with Git')
    rows=study.read_csv(folder/'spectrum.csv')
    by={(r['pair'],float(r['mu']),float(r['beta_deg'])):r for r in rows}
    assert len(by)==510
    for r in rows:
        base=by[r['pair'],float(r['mu']),0.]
        assert float(r['omega1'])<float(r['omega2'])<float(r['omega3'])
        assert float(r['F_omega'])==pytest.approx((float(r['Lambda1'])/float(base['Lambda1']))**2,rel=1e-13)
    for r in study.read_csv(folder/'signed_comparison.csv'):
        m,b=float(r['abs_mu']),float(r['beta_deg'])
        pos,neg=by[r['pair'],m,b],by[r['pair'],-m,b]
        assert float(r['D_Lambda'])==pytest.approx(float(pos['Lambda1'])/float(neg['Lambda1'])-1,abs=1e-15)
        assert float(r['D_F'])==pytest.approx(float(pos['F_omega'])/float(neg['F_omega'])-1,abs=1e-15)
    quality=study.read_csv(folder/'leading_all_quality.csv')
    assert len(quality)==50
    assert all(r['status']=='PASS' for r in quality if r['role']=='production')
    assert all(r['status']=='DIAGNOSTIC_NOT_ACCEPTED' for r in quality if r['role']=='convergence_probe')


def test_plot_only_dispatch_never_initializes_experiment(monkeypatch):
    calls=[]
    def forbidden():raise AssertionError('Experiment/solver reached in plot-only')
    monkeypatch.setattr(study,'Experiment',forbidden)
    monkeypatch.setattr(study,'verify_bundle',lambda p,figures=True:calls.append(('verify',figures)))
    monkeypatch.setattr(study,'render',lambda p:calls.append(('render',None)))
    study.main(SimpleNamespace(verify_cache=None,plot_only=study.BASE/'length_asymmetry/test',
                               reuse_points_from=None,plot_set='original'))
    assert calls==[('verify',False),('render',None),('verify',True)]
