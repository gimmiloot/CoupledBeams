"""Addressed checks of the new outer support, normalization and reuse."""
import json
import math
from pathlib import Path
import numpy as np
import pytest
from scripts.lib import isotropic_two_springs_external_axial_spring as p
from scripts.lib import two_springs_external_axial_spring as workflow


def test_rigid_preserves_canonical_assembly_and_private_hooks():
    m=p.Model();case=p.Case(6,.01,.1,10)
    assert np.array_equal(p.boundary(m,case,m.omega(3)),p.base.boundary(m,case,m.omega(3)))
    assert p.base.external_matrix is not p.external_matrix
    assert p._engine.external_matrix is p.external_matrix
    for a,b in zip(p.trace_map(m,case),p.base.trace_map(m,case)):assert np.array_equal(a,b)


@pytest.mark.parametrize('k',[0.,.01,1.])
def test_local_outer_dof_and_independent_endpoint_energy_sign(k):
    m=p.Model();case=p.Case(11,.01,.1,10,k);L=case.lengths(m)[1]
    trace,stiffness=p.trace_map(m,case)
    assert np.array_equal(trace[:,-1],np.eye(12)[:,6])
    assert stiffness[-1,-1]==k*m.section.EA/m.l
    a=2.;omega=.03;effective=m.section.EA-m.J*omega**2
    state=np.array([a,0,0,-effective*a/L,0,0])
    endpoint=p.external_matrix(m,case,1)[0]@state
    assert endpoint==pytest.approx((effective/L+k*m.section.EA/m.l)*a,rel=1e-14)
    assert np.array_equal(p.joint_matrix(m,case),p.base.joint_matrix(m,case))


def test_leading_compliance_limits_and_exact_energy_identity():
    assert len(workflow.compliance_check())==3
    assert p.leading_scale(0)==0 and p.leading_scale('RIGID')==1
    for k in [.01,.1,1.]:assert p.leading_scale(k)**2/2==pytest.approx(1/(2+1/k))
    with pytest.raises(ValueError):p.Case(0,chi_u=-.1)


def test_saved_physical_results_and_full_mass():
    exp=workflow.Experiment();folder=exp.folder
    if not (folder/'numerical_summary.json').exists():pytest.skip('Run the bounded numerical study first')
    data=workflow.read_csv(folder/'spectrum.csv')
    assert len(data)==510 and sum(r['origin']=='archive' for r in data)==102
    for row in data:
        pair,k,b=row['pair'],workflow.chi(row['chi_u']),float(row['beta_deg'])
        rec=json.loads((workflow.ROOT/row['point_file']).read_text())
        case=p.base.Case(b,.01,*exp.cfg['cases'][pair]) if k=='RIGID' else exp.case(pair,k,b)
        exp.validate(rec,case)
        zero=next(r for r in data if (r['pair'],r['chi_u'],float(r['beta_deg']))==(pair,row['chi_u'],0.))
        assert float(row['F_omega'])==pytest.approx((float(row['Lambda1'])/float(zero['Lambda1']))**2,rel=2e-14)
    # Independent integration of every inertial term, including Rayleigh--Love.
    row=next(r for r in data if r['pair']=='J0E4' and workflow.chi(r['chi_u'])==.01 and float(r['beta_deg'])==6)
    with np.load((workflow.ROOT/row['point_file']).with_suffix('.npz')) as z:
        arms=[{f:z[f'mode1_arm{i}_{f}'] for f in ('weights','y','dy','ends')} for i in (1,2)]
        m=exp.run.model;mass=0.
        for arm in arms:
            y,d,wt=arm['y'],arm['dy'],arm['weights']
            mass+=wt@(m.section.rhoA*(y[:,0]**2+y[:,1]**2)+m.section.rhoI*y[:,2]**2+m.J*d[:,0]**2)
        assert mass==pytest.approx(1.,abs=1e-12)
        components=p.energy_components(m,exp.case('J0E4',.01,6.),arms)
        assert components[-1]==pytest.approx(.01*m.section.EA/m.l*arms[1]['ends'][0,0]**2)
        assert components.sum()==pytest.approx(float(row['omega1'])**2,rel=5e-8)


def test_location_comparison_uses_actual_common_nodes():
    exp=workflow.Experiment();file=exp.folder/'location_comparison.csv'
    if not file.exists():pytest.skip('Study not yet calculated')
    for r in workflow.read_csv(file):
        pair,k,b=r['pair'],workflow.chi(r['chi']),float(r['beta_deg'])
        assert (pair,k,b) in exp.joint and b in exp.cfg['map']['parameter_grid']
        assert float(r['D_location'])==pytest.approx(float(r['F_end'])/float(r['F_joint'])-1,abs=1e-15)
    for r in workflow.read_csv(exp.folder/'leading_comparison.csv'):
        assert float(r['xi'])==pytest.approx(math.radians(float(r['beta_deg']))/.005)
        assert float(r['delta_percent'])==pytest.approx(100*abs(float(r['F_full'])/float(r['F_lim'])-1))
        if workflow.chi(r['chi_u'])==0:assert float(r['F_lim'])==1.


def test_plot_only_explicitly_forbids_spectral_work():
    old=p.prefix
    with workflow.plots_without_calculation():
        with pytest.raises(RuntimeError,match='forbidden'):p.prefix(None)
        with pytest.raises(RuntimeError,match='forbidden'):workflow.Ritz.solve(None,0)
        with pytest.raises(RuntimeError,match='forbidden'):workflow.api().Run.solve(None,None)
    assert p.prefix is old
