"""Boundary work, exact free rotation, coercivity, and accepted saved evidence."""
from fractions import Fraction
import json
import math
import numpy as np
import pytest
from scipy.linalg import null_space, cholesky
from scripts.lib import isotropic_two_springs_external_transverse_spring as p
from scripts.lib.isotropic_two_springs_external_rotational_replacement import leading_model


@pytest.mark.parametrize('kw',['RIGID',1.,0.])
def test_full_boundary_virtual_work_and_free_rotation(kw):
    m=p.Model();c=p.Case(6,.01,.1,0.,kw)
    op=np.zeros((12,24));op[:3,:6]=p.external_matrix(m,c,0);op[3:6,12:18]=p.external_matrix(m,c,1)
    op[6:,6:12]=p.joint_matrix(m,c)[:,:6];op[6:,18:]=p.joint_matrix(m,c)[:,6:]
    kernel=null_space(op/np.linalg.norm(op,axis=1)[:,None]);trace,s=p.trace_map(m,c)
    rng=np.random.default_rng(311);state=(kernel@rng.normal(size=kernel.shape[1])).reshape(4,6)
    v=rng.normal(size=trace.shape[1]);delta=(trace@v).reshape(4,3)
    coord=np.linalg.lstsq(trace,state[:,:3].flatten(),rcond=None)[0]
    work=sum(sign*a[3:]@b for sign,a,b in zip((-1,1,-1,1),state,delta))+coord@s@v
    assert abs(work)<1e-12
    assert np.count_nonzero(trace[8])>0  # psi2(0) is admissible
    assert np.count_nonzero(trace[6])==0  # u2(0) is fixed
    assert np.count_nonzero(trace[7])==(0 if kw=='RIGID' else 1)
    assert np.array_equal(p.external_matrix(m,c,1)[2],[0,0,0,0,0,-1])


@pytest.mark.parametrize('kj',[.1,10.])
@pytest.mark.parametrize('kw',['RIGID',1.,0.])
def test_leading_exact_same_space_mass_and_omitted_rotation(kj,kw):
    new=leading_model(kj,kw);old=p.Leading(kj,10.,order=24,k_w=kw,analytic_curvature=True)
    for name in ('Kb','Kj','M','P','Z','constraint'):assert np.array_equal(getattr(new,name),getattr(old,name))
    assert new.k_end==0 and np.count_nonzero(new.Ke)==0
    assert np.array_equal(new.K0,new.Kb+new.Kj+(new.Kw if kw!='RIGID' else 0))
    assert np.linalg.norm(new.outer[1][1]@new.Z)>0
    assert np.linalg.norm(new.outer[1][0]@new.Z)>0 if kw!='RIGID' else np.count_nonzero(new.outer[1][0])==0
    cholesky(new.Z.T@new.K0@new.Z);cholesky(new.Z.T@new.M@new.Z)


def test_axial_energy_exact_and_coefficient_unchanged():
    L1,L2=Fraction(99,100),Fraction(101,100);q=Fraction(7,3)
    force=-q/(L1+L2);d1,d2=force*L1,force*L2
    assert d1+d2==-q
    assert (d1*d1/L1+d2*d2/L2)/2==q*q/4


def test_saved_map_acceptance_and_own_normalizations():
    from scripts.lib import two_springs_external_rotational_replacement as w
    exp=w.Experiment();path=exp.folder/'spectrum.csv'
    if not path.exists():pytest.skip('Main map not yet accepted')
    rows=w.read_csv(path);assert len(rows)==510
    zero={(r['pair'],r['k_w']):float(r['Lambda1']) for r in rows if float(r['beta_deg'])==0}
    for r in rows:
        c=exp.case(r['pair'],w.stiffness(r['k_w']),float(r['beta_deg']));rec=json.loads((w.ROOT/r['point_file']).read_text());exp.validate(rec,c)
        assert float(r['F_omega'])==pytest.approx((float(r['Lambda1'])/zero[r['pair'],r['k_w']])**2,rel=2e-14)
    for r in w.read_csv(exp.folder/'leading_comparison.csv'):
        assert float(r['xi'])==pytest.approx(math.radians(float(r['beta_deg']))/.005)
        assert 0<=float(r['xi'])<=30
        assert float(r['delta_percent'])==pytest.approx(100*abs(float(r['F_full'])/float(r['F_lim'])-1))
    for r in w.read_csv(exp.folder/'mechanics.csv'):
        assert float(r['mass'])==pytest.approx(1.,abs=1e-11)
        assert float(r['end_rotation_fraction'])==0.
    for r in w.read_csv(exp.folder/'leading_curvature.csv'):
        assert float(r['C'])==pytest.approx(float(r['node0_squared'])/(4*float(r['zeta0'])),rel=1e-14)
