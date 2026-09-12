"""Bounded mechanics/orchestration checks; no spectrum search or parameter map."""
import json
import numpy as np
import pytest
from scipy.linalg import expm

from scripts.lib import inplane_spring_modes as m
from scripts.analysis.laminated_beams import check_inplane_spring_robustness as flow


def test_native_heterogeneous_reduction():
    section,p=m.section()
    assert np.linalg.norm(section.B)==0 and section.I1==0
    assert p.A==pytest.approx(.011)
    assert p.D==pytest.approx(1.1*.20*.05**3/12*1.3)
    assert p.S==pytest.approx((5/6)*(.20*.05)/2.6)
    assert p.m==pytest.approx(.01) and p.J==pytest.approx(.20*.05**3/12)
    assert max(p.axial_reduction.relative_difference,p.bending_reduction.relative_difference,
               p.shear_reduction_before_K.relative_difference)<1e-12


@pytest.mark.parametrize('model',['EB','RLB'])
def test_actual_lengths_in_transfer_and_shapes(model):
    p=m.section()[1];pair=m.arms(model,.01,p)
    assert [a.L for a in pair]==[.99,1.01]
    r=np.array([.002,1e-6,-2e-6]);omega=.3
    computed=[]
    for arm in pair:
        states,_,errors=m.states_along_arm(omega,arm,r,check=True)
        assert max(errors)<1e-9
        scale=arm.scale();h=arm.matrix(omega)*scale[None,:]/scale[:,None]
        expected=scale*(expm(h*arm.L)@np.r_[np.zeros(3),r/scale[3:]])
        np.testing.assert_allclose(states[-1]/scale,expected/scale,atol=1e-9,rtol=1e-9)
        computed.append(states[-1])
    assert not np.allclose(computed[0],computed[1],rtol=1e-4,atol=1e-10)


@pytest.mark.parametrize('model',['EB','RLB'])
def test_mass_contains_each_length_and_only_rlb_rotation(model):
    p=m.section()[1];pair=m.arms(model,.01,p)
    xi,w=m.modes.quadrature()
    y=np.zeros((2,len(xi),6));y[0,:,0]=xi;y[1,:,0]=3*xi;y[:,:,2]=2
    vector=m.physical_vector(y,pair,w)
    expected=p.m*(.99+9*1.01)/3+(8*p.J if model=='RLB' else 0)
    assert np.vdot(vector,vector).real==pytest.approx(expected)
    normalized=y/np.sqrt(expected)
    np.testing.assert_allclose(normalized[1,:,0],3*normalized[0,:,0])
    assert np.linalg.norm(m.physical_vector(normalized,pair,w))==pytest.approx(1)


def test_rlb_never_uses_eb_shape_function(monkeypatch):
    monkeypatch.setattr(m.modes,'arm_states',lambda *a,**k:pytest.fail('EB reconstruction called for RLB'))
    p=m.section()[1];arm=m.Arm('RLB',p,1.01)
    states,calls,_=m.states_along_arm(.2,arm,[1e-3,1e-6,1e-6])
    assert states.shape==(129,6) and calls=={'analytic':0,'expm':1}


def test_no_projection_or_class_constraint_for_asymmetry():
    p=m.section()[1];pair=m.arms('RLB',.01,p)
    with pytest.raises(ValueError,match='projection forbidden'):
        m.recover(None,.2,pair,.3,m.eb.Joint('SPRING',0),parity=1)
    left=[dict(common=np.array([1,0]),symmetry_class=1,branch_id='mode_01'),
          dict(common=np.array([0,1]),symmetry_class=-1,branch_id='mode_02')]
    right=[dict(common=np.array([0,1]),symmetry_class=None),dict(common=np.array([1,0]),symmetry_class=None)]
    indices,mac,margin=m.match(left,right)
    assert indices.tolist()==[1,0]
    assert min(mac)==1 and min(margin)==1
    with pytest.raises(ValueError,match='broken symmetry'):
        m.match(left,right,restrict_symmetry=True)


@pytest.mark.parametrize('model,Omega',[('EB',18.07914233661534),('RLB',17.851330709400546)])
def test_mixed_full_matrix_mode_is_not_projected_or_rejected(model,Omega):
    # Regression frequencies from this bounded check, not an independent
    # reference solver. No root search is executed by this test.
    p=m.section()[1];pair=m.arms(model,.01,p)
    joint=m.eb.Joint('SPRING',m.REFERENCE.D);beta=np.pi/6
    assembled=m.assembly(Omega/m.FREQUENCY_SCALE,pair,beta,joint)
    result=m.recover(assembled,Omega/m.FREQUENCY_SCALE,pair,beta,joint)
    assert not result['failures'] and result['symmetry_class'] is None
    assert .4<result['reflection_weight_plus']<.6
    assert result['M']==pytest.approx(1) and result['s']>0
    assert max(result['physical_residuals'])<1e-9
    assert not np.allclose(result['states'][1],result['states'][0]@m.modes.REFLECTION)
    assert not np.allclose(result['states'][1],-result['states'][0]@m.modes.REFLECTION)


def test_hinge_matrix_and_independent_sensitivity_at_zero():
    matrix=m.eb.joint_matrix(.4,m.eb.Joint('SPRING',0))
    np.testing.assert_array_equal(matrix[2],np.eye(12)[5])
    np.testing.assert_array_equal(matrix[5],np.eye(12)[5]+np.eye(12)[11])
    delta,omega,mass=3.,2.,5.
    s=m.REFERENCE.D*delta**2/(omega**2*mass)
    assert s>0
    assert m.REFERENCE.D*(7*delta)**2/(omega**2*(49*mass))==pytest.approx(s)


def test_common_metric_excludes_rigidities_lengths_and_rotational_mass():
    _,w=m.modes.quadrature();rng=np.random.default_rng(31)
    y=rng.normal(size=(2,129,6));z=y.copy();z[:,:,2:]*=1e8
    np.testing.assert_array_equal(m.common_vector(y,w),m.common_vector(z,w))
    assert abs(np.vdot(m.common_vector(2j*y,w),m.common_vector(y,w)))==pytest.approx(1)


def test_seed_assignment_can_retain_seventh_position_and_multiplicity():
    vectors=np.eye(8)
    left=[dict(common=vectors[j],symmetry_class=None) for j in (0,1,2,3,4,6)]
    right=[dict(common=v,current_sorted_position=j+1,Omega=4. if j in (3,4) else j+1,
                symmetry_class=None) for j,v in enumerate(vectors)]
    indices,_,_=m.match(left,right)
    assert indices.tolist()==[0,1,2,3,4,6]
    assert right[3]['Omega']==right[4]['Omega'] and len(indices)==6
    with pytest.raises(ValueError,match='distinct'):
        m.modes.orthonormal_subspace(vectors[:2],[1.,1.00001])


@pytest.mark.parametrize('bad',[float('nan'),float('inf'),1.,-1.])
def test_invalid_mu(bad):
    with pytest.raises(ValueError):m.arms('RLB',bad,m.section()[1])


def test_curve_keeps_rejected_node_as_gap_without_interpolation():
    def row(beta,frequency,status='CONFIRMED'):
        return dict(model='RLB',mu=.01,kappa=1,branch_id='mode_01',beta_deg=beta,
                    Lambda=frequency,mapping_status=status,root_status='CONFIRMED')
    rows=[row(20,4),row(21,4.1,'MAPPING_AMBIGUOUS'),row(22,4.2)]
    x,y=flow.curve(rows,'RLB',.01,1,'mode_01')
    assert x==[20,21,22] and np.isnan(y[1]) and y[0]==4 and y[2]==4.2


def test_plot_only_zero_computational_calls_and_immutable_sources(tmp_path,monkeypatch):
    def forbidden(*args,**kwargs):pytest.fail('plot_only made a computational call')
    monkeypatch.setattr(flow,'OUTPUT',tmp_path)
    for obj,name in [(flow,'Run'),(flow.PointMatrices,'matrix'),(m,'recover'),(m,'match'),
                     (flow.pilot.roots,'_scan_candidates')]:monkeypatch.setattr(obj,name,forbidden)
    rows=[]
    for model in ('EB','RLB'):
        for window in ('A','B'):
            for mu in (0.,.01):
                for kappa in ((0.,1.,100.) if window=='A' else (1.,)):
                    for branch in (('mode_01',) if window=='A' else ('mode_01','mode_02')):
                        for beta in (8.,9.):rows.append(dict(model=model,window=window,mu=mu,kappa=kappa,
                            branch_id=branch,beta_deg=beta,Lambda=3+beta/20,mapping_status='CONFIRMED',root_status='CONFIRMED'))
    flow.write_csv(tmp_path/'local_modes.csv',rows)
    before=(tmp_path/'local_modes.csv').read_bytes()
    flow.render()
    assert (tmp_path/'local_modes.csv').read_bytes()==before
    result=json.loads((tmp_path/'render_manifest.json').read_text())
    assert [result[k] for k in ('solver_calls','matrix_calls','shape_calls','tracking_calls')]==[0,0,0,0]
