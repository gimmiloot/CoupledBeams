"""Addressed physical, leading-space and archive regression checks."""
from fractions import Fraction
import json
import math
import numpy as np
import pytest
from scripts.lib import isotropic_two_springs_external_transverse_spring as p
from scripts.lib import two_springs_external_transverse_spring as w


def test_rigid_assembly_and_no_global_mutation():
    m=p.Model();c=p.Case(6,.01,.1,10)
    assert np.array_equal(p.boundary(m,c,m.omega(3)),p.base.boundary(m,c,m.omega(3)))
    assert p.base.external_matrix is not p.external_matrix
    for a,b in zip(p.trace_map(m,c),p.base.trace_map(m,c)):assert np.array_equal(a,b)
    model=p.Leading(.1,10,order=20);old=p.leading.Ritz(.1,10,order=20)
    for field in ('M','K0','P','constraint'):assert np.array_equal(getattr(model,field),getattr(old,field))


@pytest.mark.parametrize('k',[0.,.1,10.])
def test_external_sign_dimension_and_exact_axial_support(k):
    m=p.Model(l=2.);c=p.Case(6,.01,.1,10,k);L=c.lengths(m)[1];a=2.
    state=np.array([0,a,0,0,-m.section.KGA*a/L,0])
    value=p.external_matrix(m,c,1)[1]@state
    assert value==pytest.approx((m.section.KGA/L+k*m.section.EI/m.l**3)*a)
    assert np.array_equal(p.external_matrix(m,c,1)[0],np.eye(6)[0])
    trace,s=p.trace_map(m,c)
    assert np.array_equal(trace[:,-1],np.eye(12)[:,7]) and np.count_nonzero(trace[6])==0
    assert s[-1,-1]==k*m.section.EI/m.l**3
    assert np.array_equal(p.joint_matrix(m,c),p.base.joint_matrix(m,c))


def test_leading_free_constant_has_consistent_mass_and_no_essential_constraint():
    model=p.Leading(.1,10,order=20,k_w=.1)
    assert model.M.shape==(41,41) and model.Z.shape==(41,40)
    assert model.M[-1,-1]==pytest.approx(1.01)
    assert np.linalg.norm(model.M[:-1,-1])>0
    assert model.outer[1][0][-1]==1 and model.constraint[0,-1]==1
    assert model.K0[-1,-1]==pytest.approx(.1)
    assert model.P[-1,-1]==1
    assert np.linalg.eigvalsh(model.Z.T@model.M@model.Z).min()>0


def test_leading_natural_force_sign_by_integration():
    # A cubic W and a variation with only delta W(0)=1 nonzero at the ends.
    L=1.01;a,b,c,d=.7,-.2,.3,.4;k=.1
    nodes,weights=np.polynomial.legendre.leggauss(8);x=(nodes+1)*L/2
    W2=2*c+6*d*x;delta2=(-6+12*x/L)/L**2
    variation=weights@(W2*delta2)*L/2+k*a
    assert variation==pytest.approx(6*d+k*a,rel=1e-13)


def test_unchanged_axial_minimum_exact_rational_arithmetic():
    L1,L2=Fraction(99,100),Fraction(101,100);q=Fraction(7,3)
    force=-q/(L1+L2);d1,d2=force*L1,force*L2
    assert d1+d2==-q
    assert (d1*d1/L1+d2*d2/L2)/2==q*q/4


def test_saved_spectra_normalization_and_all_physical_gates():
    exp=w.Experiment();path=exp.folder/'spectrum.csv'
    if not path.exists():pytest.skip('Bounded study not computed')
    rows=w.read_csv(path)
    if len(rows)!=510:pytest.skip('Bounded study is still running')
    zeros={(r['pair'],r['k_w']):float(r['Lambda1']) for r in rows if float(r['beta_deg'])==0}
    for r in rows:
        pair,k,b=r['pair'],w.stiffness(r['k_w']),float(r['beta_deg']);rec=json.loads((w.ROOT/r['point_file']).read_text())
        case=p.base.Case(b,.01,*exp.cfg['cases'][pair]) if k=='RIGID' else exp.case(pair,k,b)
        exp.validate(rec,case)
        assert float(r['F_omega'])==pytest.approx((float(r['Lambda1'])/zeros[pair,r['k_w']])**2,rel=2e-14)
    assert sum(r['origin']=='archive' for r in rows)==102


def test_saved_leading_coefficient_and_actual_coordinates():
    exp=w.Experiment();path=exp.folder/'leading_comparison.csv'
    if not path.exists():pytest.skip('Leading study not computed')
    for r in w.read_csv(path):
        assert float(r['xi'])==pytest.approx(math.radians(float(r['beta_deg']))/.005)
        assert float(r['beta_deg']) in exp.cfg['map']['parameter_grid'] and float(r['xi'])<=30
        assert float(r['delta_percent'])==pytest.approx(100*abs(float(r['F_full'])/float(r['F_lim'])-1))
    for r in w.read_csv(exp.folder/'leading_curvature.csv'):
        assert float(r['C'])==pytest.approx(float(r['node0'])**2/(4*float(r['zeta0'])),rel=1e-14)


def test_numerical_stop_is_retained_without_restarting_solves():
    exp=w.Experiment();path=exp.folder/'leading_gate.json'
    if not path.exists():pytest.skip('No saved leading diagnosis')
    gate=json.loads(path.read_text())
    assert gate['status']=='STOP' and gate['strong_ode']>gate['tolerance']
    assert not gate['coefficient_family_accepted'] and not gate['full_leading_comparison_completed']
    assert gate['identity']['pair']=='J0E4' and gate['identity']['k_w']==1 and gate['identity']['order']==24
    assert gate['strong_per_position'][0]==pytest.approx(gate['strong_ode'],rel=1e-14)


def test_saved_data_guard_blocks_full_and_leading_calculation():
    old=p.Leading.solve
    with w.prohibit_spectral_calls():
        with pytest.raises(RuntimeError,match='forbidden'):p.prefix(None)
        with pytest.raises(RuntimeError,match='forbidden'):p.Leading.solve(None,0)
        with pytest.raises(RuntimeError,match='forbidden'):w.api().Run.solve(None,None)
    assert p.Leading.solve is old
