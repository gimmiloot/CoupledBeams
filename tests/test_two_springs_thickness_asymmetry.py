"""Physical and algebraic checks for individual sections; no broad sweeps."""
import json
from types import SimpleNamespace
import numpy as np
import pytest
from scripts.lib import isotropic_two_springs_thickness_asymmetry as p
from scripts.lib import two_springs_thickness_asymmetry as study


def test_exact_polynomial_identities_and_grid_symmetries():
    assert study.polynomial_identity()['status']=='PASS'
    for mu in (-.1,0,.1):
        for tau in (-.1,0,.1):
            m=p.Model(mu=mu,tau=tau)
            h,a,i=p.ratios(mu,tau)
            assert a@np.array([1-mu,1+mu])==pytest.approx(2,abs=2e-15)
            assert m.S_ax==pytest.approx(p.Model(mu=-mu,tau=-tau).S_ax,abs=2e-15)
            assert m.S_ax==pytest.approx(p.Model(mu=mu,tau=-tau).S_ax,abs=2e-15)
            assert i[1]/i[0]==pytest.approx(((1+tau)/(1-tau))**4,rel=1e-14)
            assert (h[1]-h[0])/(h[1]+h[0])==pytest.approx(tau,abs=1e-15)


def test_reference_springs_frequency_and_unchanged_parent():
    for mu in (-.1,0,.1):
        for tau in (-.1,0,.1):
            m=p.Model(mu=mu,tau=tau);c=p.Case(6,mu,.1,10,tau);base=p.base.Model()
            assert m.time==base.time
            np.testing.assert_array_equal(p.trace_map(m,c)[1],p.base.trace_map(base,c)[1])
            np.testing.assert_array_equal(p.external_matrix(m,c,1),p.base.external_matrix(base,c,1))
            for a in m.arms:assert a.J==pytest.approx(17*a.nu**2*a.section.rhoI)
            if tau==0:
                np.testing.assert_array_equal(p.boundary(m,c,.08),p.base.boundary(base,c,.08))
    assert p.base.boundary is not p._engine.boundary
    assert p.base.count_at is not p._engine.count_at


def test_individual_mass_and_energy_against_polynomial_integrals():
    m=p.Model(mu=.1,tau=-.1);case=p.Case(6,.1,.1,10,-.1)
    nodes,weights=np.polynomial.legendre.leggauss(8);arms=[];expected_mass=[];expected_elastic=[]
    constants=[(2.,3.,5.),(-4.,2.,-3.)]
    for arm,L,(a,b,c) in zip(m.arms,case.lengths(m),constants):
        x=(nodes+1)*L/2;y=np.zeros((len(x),6));dy=np.zeros_like(y);s=arm.section
        y[:,0]=a*x;y[:,1]=b*x;y[:,2]=c;dy[:,0]=a;dy[:,1]=b
        ends=np.array([[0,0,c,0,0,0],[a*L,b*L,c,0,0,0]])
        arms.append(dict(x=x,weights=weights*L/2,y=y,dy=dy,ends=ends))
        expected_mass.append([s.rhoA*a*a*L**3/3,s.rhoA*b*b*L**3/3,s.rhoI*c*c*L,arm.J*a*a*L])
        expected_elastic.append([s.EA*a*a*L,0,s.KGA*(b-c)**2*L])
    mass,elastic,springs=p.energies(m,case,arms)
    np.testing.assert_allclose(mass,expected_mass,rtol=2e-14,atol=1e-18)
    np.testing.assert_allclose(elastic,expected_elastic,rtol=2e-14,atol=1e-18)
    np.testing.assert_allclose(springs,[.1*m.section.EI*8**2,10*m.section.EI*3**2],rtol=2e-14)


def test_section_specific_pole_keys_do_not_alias_equal_lengths():
    m=p.Model(mu=0,tau=.1)
    assert p.pole_key(m.arms[0],1)!=p.pole_key(m.arms[1],1)
    assert p.pole_key(p.Model(tau=0).arms[0],1)==p.pole_key(p.base.Model(),1)


def test_leading_weighted_matrices_and_first_tau_variation():
    mu=.1;h=1e-5
    zero=p.Leading(.1,10,mu=mu)
    old=p.leading.Ritz(.1,10,mu=mu,order=20)
    np.testing.assert_array_equal(zero.K0,old.K0)
    np.testing.assert_array_equal(zero.M,old.M)
    plus=p.Leading(.1,10,mu=mu,tau=h);minus=p.Leading(.1,10,mu=mu,tau=-h)
    da=np.array([-2*(1+mu),2*(1-mu)])
    np.testing.assert_allclose((plus.K0-minus.K0)/(2*h),sum(2*a*K for a,K in zip(da,zero.arm_K)),rtol=1e-8,atol=1e-8)
    np.testing.assert_allclose((plus.M-minus.M)/(2*h),sum(a*M for a,M in zip(da,zero.arm_M)),rtol=1e-8,atol=1e-8)
    assert plus.S_ax==pytest.approx(minus.S_ax,abs=1e-15)
    # Energy minimization: half the added quadratic stiffness is potential energy.
    xi,y=3.,2.
    L=np.array(plus.lengths);A=plus.area
    extension=-xi*y*L/A/plus.S_ax
    assert extension.sum()==pytest.approx(-xi*y)
    assert .5*np.sum(A*extension**2/L)==pytest.approx(xi**2*y**2/(2*plus.S_ax))


def test_saved_spectrum_normalization_signed_comparison_and_quality():
    folder=study.BASE/'12ff3907cf7cd6da'
    if not (folder/'numerical_summary.json').exists():pytest.skip('Local research arrays not present')
    rows=study.read_csv(folder/'spectrum.csv');assert len(rows)==918
    by={(r['pair'],float(r['mu']),float(r['tau']),float(r['beta_deg'])):r for r in rows}
    for r in rows:
        base=by[r['pair'],float(r['mu']),float(r['tau']),0.]
        assert float(r['F_omega'])==pytest.approx((float(r['Lambda1'])/float(base['Lambda1']))**2,rel=1e-13)
        assert float(r['omega1'])<float(r['omega2'])<float(r['omega3'])
    for r in study.read_csv(folder/'signed_comparison.csv'):
        plus,minus=[by[r['pair'],float(r['mu']),tau,float(r['beta_deg'])] for tau in (.1,-.1)]
        assert float(r['D_F'])==pytest.approx(float(plus['F_omega'])/float(minus['F_omega'])-1,abs=1e-15)
    assert all(r['status']=='PASS' for r in study.read_csv(folder/'leading_quality.csv'))


@pytest.mark.parametrize('operation',['plot_only','verify_cache'])
def test_saved_data_modes_never_construct_a_calculation(operation,monkeypatch):
    def forbidden():
        raise AssertionError('Saved-data operation initialized a spectral experiment')
    events=[]
    monkeypatch.setattr(study,'Experiment',forbidden)
    monkeypatch.setattr(study,'verify_bundle',lambda folder,figures=True:events.append(('verify',figures)))
    monkeypatch.setattr(study,'render',lambda folder:events.append(('render',)))
    args=SimpleNamespace(plot_only=None,verify_cache=None)
    setattr(args,operation,study.BASE/'12ff3907cf7cd6da')
    study.main(args)
    assert events==([('verify',False),('render',),('verify',True)] if operation=='plot_only' else [('verify',True)])
