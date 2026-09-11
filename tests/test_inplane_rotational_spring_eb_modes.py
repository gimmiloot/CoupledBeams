"""Targeted physical-shape/algebra/orchestration tests; no angular sweep."""
import csv
import json

import numpy as np
import pytest
from scipy.linalg import expm

from scripts.lib import inplane_rotational_spring_eb_modes as m
from scripts.analysis.laminated_beams import track_inplane_rotational_spring_eb as workflow

ARM=workflow.ARM


@pytest.mark.parametrize('omega',[0.,.07,.5,1.5])
def test_analytic_shape_equals_full_transfer_and_physical_signs(omega):
    xi=np.array([0.,.17,.53,1.]);r=np.array([.02,.000003,-.000002])
    shape=m.arm_states(omega,ARM,r,xi)
    scale=m.eb.state_scale(ARM);H=m.eb.state_matrix(omega,ARM)*scale[None,:]/scale[:,None]
    expected=np.array([scale*(expm(H*x)@np.r_[np.zeros(3),r/scale[3:]]) for x in xi])
    np.testing.assert_allclose(shape/scale,expected/scale,rtol=2e-9,atol=2e-9)
    assert np.max(abs(shape[0,:3]))==0
    np.testing.assert_allclose(shape[0,3:],r,rtol=1e-14)


def test_whole_structure_mass_and_relative_arm_amplitudes():
    xi,w=m.quadrature(65)
    states=np.zeros((2,65,6));states[0,:,0]=xi;states[1,:,0]=3*xi
    normalized,reactions,vector,mass=m.normalize(states,np.ones((2,3)),ARM,w)
    assert mass==pytest.approx(10*ARM.m*ARM.L/3)
    assert np.vdot(vector,vector)==pytest.approx(1)
    np.testing.assert_allclose(normalized[1,:,0],3*normalized[0,:,0])
    assert reactions[0,0]==reactions[1,0]


def test_mass_excludes_rotations_and_forces():
    _,w=m.quadrature(65);states=np.ones((2,65,6));other=states.copy();other[:,:,2:]*=1e8
    np.testing.assert_array_equal(m.mass_vector(states,ARM,w),m.mass_vector(other,ARM,w))


def test_mac_amplitude_phase_and_global_assignment():
    left=np.eye(3,dtype=complex);right=np.array([2j*left[2],-7*left[0],.01*left[1]])
    columns,mac,margin=m.assign(left,right)
    assert columns.tolist()==[1,2,0] and len(set(columns))==3
    np.testing.assert_allclose(mac[np.arange(3),columns],1)
    np.testing.assert_allclose(margin,1)


@pytest.mark.parametrize('complex_basis',[False,True])
def test_repeated_subspace_basis_invariance_and_no_distinct_root_mixing(complex_basis):
    rng=np.random.default_rng(42);matrix=rng.normal(size=(20,3))
    if complex_basis:matrix=matrix+1j*rng.normal(size=(20,3))
    v=np.linalg.qr(matrix)[0].T
    transform=np.array([[2.,1,0],[0,3,1],[1,0,2]],dtype=complex if complex_basis else float)
    if complex_basis:transform[0,1]+=2j
    a=m.orthonormal_subspace(v,[4,4,4]);b=m.orthonormal_subspace(transform@v,[4,4,4])
    np.testing.assert_allclose(m.principal_correlations(a,b),1,atol=1e-14)
    with pytest.raises(ValueError,match='distinct'):
        m.orthonormal_subspace(v,[4,4.001,4])
    with pytest.raises(ValueError,match='distinct'):
        m.orthonormal_subspace(v,[4,4+1e-13,4])


@pytest.mark.parametrize('beta',[0.,.1,24.5,47.,90.])
def test_reflection_and_full_joint_row_space(beta):
    k=100*ARM.D/ARM.L;result=m.symmetry_equivalence(np.deg2rad(beta),k,ARM)
    assert result['row_rank']==6 and result['row_error']<1e-12
    assert result['reflection_commutator']==0
    values=np.arange(2*7*6).reshape(2,7,6)
    np.testing.assert_array_equal(m.reflect(m.reflect(values)),values)
    for parity in (1,-1):
        rng=np.random.default_rng(11);endpoint=rng.normal(size=(6,3))
        J=m.eb.joint_matrix(np.deg2rad(beta),m.eb.Joint('SPRING',k))
        full=J@np.vstack([endpoint,parity*m.REFLECTION@endpoint])
        reduced=m.class_conditions(np.deg2rad(beta),k,parity)@endpoint
        assert np.linalg.matrix_rank(full)==np.linalg.matrix_rank(reduced)==3


def test_known_crossing_is_not_nearest_frequency_relabeling():
    left=np.eye(2);right=np.eye(2)[::-1]
    # Stable eigenvectors, but sorted eigenvalues exchanged.
    columns,_,_=m.assign(left,right,[1,-1],[-1,1])
    values=np.array([1.,2.]);assigned=values[columns]
    assert assigned.tolist()==[2.,1.]
    assert sorted(assigned)==sorted(values)


def test_veering_requires_small_steps_instead_of_character_jump():
    def basis(theta):return np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
    a=basis(0);b=basis(np.pi/4)
    _,mac,margin=m.assign(a,b)
    assert np.max(mac)<workflow.CRITERIA['mac'] and np.max(margin)<workflow.CRITERIA['margin']
    previous=a
    for theta in np.linspace(0,np.pi/2,33)[1:]:
        current=basis(theta);columns,mac,_=m.assign(previous,current)
        assert columns.tolist()==[0,1] and min(mac[np.arange(2),columns])>.99
        previous=current


def test_column_equilibration_and_physical_reaction_roundtrip():
    omega=.7;beta=.31;joint=m.eb.Joint('SPRING',ARM.D/ARM.L)
    assembly=m.eb.boundary_assembly(omega,ARM,ARM,beta,joint,ARM)
    eq=m.eb.positively_equilibrate_matrix(assembly.dimensionless)
    z=np.array([.2,-.3,.4,.6,-.2,.7]);hat=eq.column_factors*z
    physical=assembly.reaction_scales*hat
    ends=np.concatenate([m.arm_states(omega,ARM,r,[1.])[0] for r in physical.reshape(2,3)])
    np.testing.assert_allclose(ends,assembly.endpoint_map@hat,rtol=1e-10,atol=1e-10)
    np.testing.assert_allclose(assembly.physical@physical,m.eb.scalar_joint_residuals(ends,beta,joint),rtol=1e-10,atol=1e-10)
    np.testing.assert_allclose(eq.scaled_matrix@z,eq.row_factors*(assembly.dimensionless@hat),rtol=1e-12,atol=1e-12)


def test_full_boundary_and_two_blocks_have_same_transformed_equations():
    beta=.62;k=ARM.D/ARM.L
    context=workflow.PointMatrices(1,np.rad2deg(beta));a=context.assembly(.7)
    eye=np.eye(6);Q=np.block([[eye,eye],[m.REFLECTION,-m.REFLECTION]])/np.sqrt(2)
    f=np.diag([1.,-1.,-1.]);R=np.block([[np.eye(3),np.eye(3)],[f,-f]])/np.sqrt(2)
    np.testing.assert_allclose(a.endpoint_map@R,Q@a.endpoint_map,rtol=1e-14,atol=1e-14)
    # Each column in a class spans precisely the original full state image.
    for index,parity in enumerate((1,-1)):
        selected=a.endpoint_map@R[:,3*index:3*(index+1)]
        np.testing.assert_allclose(selected[6:],parity*m.REFLECTION@selected[:6],rtol=1e-14,atol=1e-14)
        np.testing.assert_allclose(m.class_conditions(beta,k,parity)@selected[:6]*np.sqrt(2),
            m.class_conditions(beta,k,parity)@a.endpoint_map[:6,:3],rtol=1e-14,atol=1e-14)


def test_two_recoveries_per_event_not_per_recursive_midpoint(monkeypatch):
    monkeypatch.setattr(workflow,'PLOT_KAPPAS',[1]);monkeypatch.setattr(workflow,'save',lambda *a,**k:None)
    monkeypatch.setattr(workflow,'write_csv',lambda *a,**k:None)
    def roots(beta):
        return [dict(kappa=1,beta_deg=beta,current_sorted_position=j+1,Omega=j+1.,omega=j+1.,Lambda=np.sqrt(j+1),
            root_status='CONFIRMED',symmetry_class=1,source='SYNTHETIC',shape_key=f'b{beta}_r{j}',cluster_id='') for j in range(6)]
    state=dict(points={workflow.point_id(1,b):dict(kappa=1,beta_deg=b,roots=roots(b)) for b in (0.,1.)},
        added_points=[],tracking_attempts=[])
    def always_ambiguous(previous,point,shapes):return point['roots'],np.ones(6)*.5,np.zeros(6),list(range(6))
    monkeypatch.setattr(workflow,'match',always_ambiguous)
    def insert(state,shapes,k,b,trigger):
        state['added_points'].append(workflow.point_id(k,b));state['points'][workflow.point_id(k,b)]=dict(kappa=k,beta_deg=b,roots=roots(b))
    monkeypatch.setattr(workflow,'search_point',insert)
    rows=workflow.track(state,{})
    assert len(state['tracking_attempts'])==len(state['added_points'])==2
    assert len(set(r['event_id'] for r in state['tracking_attempts']))==1
    assert all(r['Omega'] is None for r in rows if r['beta_deg']>0)


def test_plot_gaps_and_case_specific_angles_and_fixed_colors():
    rows=[dict(kappa=1,beta_deg=b,branch_id='mode_01',Lambda=v,tracking_status=s)
          for b,v,s in [(0,2,'SEED_CONFIRMED'),(1,None,'TRACKING_AMBIGUOUS'),(2,3,'TRACKED')]]
    rows.append(dict(kappa=100,beta_deg=.5,branch_id='mode_01',Lambda=9,tracking_status='TRACKED'))
    x,y=workflow.curve(rows,1,'mode_01')
    assert x==[0,1,2] and np.isnan(y[1]) and y[0]==2 and y[2]==3
    assert len(workflow.BRANCH_COLORS)==6 and workflow.BRANCH_COLORS['mode_01']=='#0072B2'


def test_plot_only_zero_compute_calls(tmp_path,monkeypatch):
    import matplotlib.figure
    monkeypatch.setattr(workflow,'OUTPUT',tmp_path)
    def forbidden(*a,**kw):raise AssertionError('compute called from plot-only')
    for name in ('PointMatrices','reconstruct','search_point','track','load'):monkeypatch.setattr(workflow,name,forbidden)
    monkeypatch.setattr(matplotlib.figure.Figure,'savefig',lambda *a,**kw:None)
    rows=[dict(kappa=k,beta_deg=b,branch_id=f'mode_{j:02d}',Lambda=j+b/100,tracking_status='TRACKED')
          for k in (1,100) for b in (0,90) for j in range(1,7)]
    workflow.write_csv(tmp_path/'tracked_branches.csv',rows)
    before=(tmp_path/'tracked_branches.csv').read_bytes();workflow.render()
    assert before==(tmp_path/'tracked_branches.csv').read_bytes()
    result=json.loads((tmp_path/'run_manifest.json').read_text())['render']
    assert result['matrix_calls']==result['root_calls']==result['shape_calls']==result['tracking_calls']==0


def test_saved_doublet_import_preserves_every_position_and_class():
    rows=[dict(Omega=w,multiplicity=n) for w,n in [(1.,1),(2.,2),(2.,2),(2.0001,1)]]
    assert m.saved_root_classes(rows)==[None,1,-1,None]
    assert [r['Omega'] for r in rows]==[1.,2.,2.,2.0001]
    with pytest.raises(ValueError,match='incomplete'):m.saved_root_classes([dict(Omega=2.,multiplicity=2)])


@pytest.mark.parametrize('Omega',[15.370552823847948,49.31932972572842,94.73879177000134,114.03153713097169])
def test_hinge_endpoint_two_independent_classes_full_matrix_and_zero_moments(Omega):
    # Saved full-precision benchmark roots; no root search in this test.
    context=workflow.PointMatrices(0,90);omega=Omega/workflow.FS;a=context.assembly(omega)
    recovered=[m.recover(a,omega,np.pi/2,m.eb.Joint('SPRING',0.),ARM,parity=eta) for eta in (1,-1)]
    for eta,r in zip((1,-1),recovered):
        assert r['detected_nullity']==2 and not r['failures'] and r['symmetry_class']==eta
        ends=r['states'][:,-1,:]
        units=np.array([ARM.L,ARM.L,1.,ARM.D/ARM.L**2,ARM.D/ARM.L**2,ARM.D/ARM.L])
        amplitude=np.max(abs(ends/units))
        assert max(abs(ends[:,5]))/(amplitude*ARM.D/ARM.L)<1e-9
    vectors=np.array([r['vector'] for r in recovered])
    np.testing.assert_allclose(vectors.conj()@vectors.T,np.eye(2),atol=1e-12)
    # At beta=90 the blocks are sign-conjugate, not literally identical.
    plus,minus=(context.block(omega,eta) for eta in (1,-1))
    np.testing.assert_allclose(minus,np.diag([1,1,-1])@plus@np.diag([1,-1,-1]),rtol=1e-12,atol=1e-10)
    assert m.symmetry_equivalence(np.pi/2,0.,ARM)['row_rank']==6


def test_unique_shape_keys_at_one_frequency(monkeypatch):
    def fake(context,frequency,parity):
        return dict(states=np.ones((2,129,6))*parity,reactions=np.ones((2,3))*parity,
                    vector=np.array([1.,parity]),symmetry_class=parity,detected_nullity=2,failures=[])
    monkeypatch.setattr(workflow,'reconstruct',fake)
    state={};shapes={};roots=[]
    for j,eta in enumerate((1,-1),1):
        root=dict(kappa=0,beta_deg=90,Omega=3.,current_sorted_position=j,multiplicity=2)
        roots.append(workflow.attach_shape(state,shapes,None,root,eta))
    assert len(shapes)==6 and len({r['shape_key'] for r in roots})==2
    assert not np.array_equal(shapes[roots[0]['shape_key']+'__vector'],shapes[roots[1]['shape_key']+'__vector'])


def test_seed_bijection_is_from_mass_mac_not_names_or_frequencies():
    reference=[dict(branch_id=f'mode_{i+1:02d}',symmetry_class=1 if i%2==0 else -1,vector=np.eye(6)[i]) for i in range(6)]
    permutation=[2,3,4,5,0,1]
    candidates=[dict(reference[i],branch_id=f'mode_{j+1:02d}') for j,i in enumerate(permutation)]
    mapping=workflow.seed_assignment(reference,candidates)
    assert [r['source_branch_id'] for r in mapping]==['mode_05','mode_06','mode_01','mode_02','mode_03','mode_04']
    assert all(r['status']=='CONFIRMED' and r['MAC']==1 for r in mapping)


def test_frozen_mapping_exact_source_frequencies_own_grids_and_gaps():
    mapping=[dict(comparison_mode_id='comparison_mode_01',kappa=k,source_branch_id='mode_03',status='CONFIRMED',source=f'case{k}') for k in (0,1)]
    def row(k,b,status='TRACKED'):
        return dict(kappa=k,beta_deg=b,branch_id='mode_03',omega='0.12345678901234567',Omega='12.345678901234567',Lambda='3.513641828',
                    current_sorted_position='5',symmetry_class='1',root_status='CONFIRMED',tracking_status=status,shape_key=f'{k}_{b}')
    source={0:[row(0,0),row(0,.05),row(0,.1,'TRACKING_AMBIGUOUS'),row(0,.2)],1:[row(1,0),row(1,.1),row(1,.2)]}
    before=json.dumps(source,sort_keys=True);rows=workflow.mapped_rows(mapping,source)
    assert json.dumps(source,sort_keys=True)==before
    assert all(r['current_sorted_position']=='5' and r['Omega']=='12.345678901234567' for r in rows)
    x,y=workflow.comparison_curve(rows,0,'comparison_mode_01');assert x==[0,.05,.1,.2] and np.isnan(y[2])
    x,y=workflow.comparison_curve(rows,1,'comparison_mode_01');assert x==[0,.1,.2] and not np.isnan(y).any()
    assert [s['linestyle'] for s in workflow.COMPARISON_STYLES.values()]==['-','--','-.']
    assert len({s['color'] for s in workflow.COMPARISON_STYLES.values()})==3


def test_comparison_plot_only_zero_compute_and_mapping_calls(tmp_path,monkeypatch):
    import matplotlib.figure
    monkeypatch.setattr(workflow,'COMPARISON_OUTPUT',tmp_path)
    def forbidden(*a,**kw):raise AssertionError('numerical operation in plot-only')
    for name in ('PointMatrices','reconstruct','search_point','track','load','comparison_inputs','seed_assignment','prepare_seed_mapping'):
        monkeypatch.setattr(workflow,name,forbidden)
    monkeypatch.setattr(matplotlib.figure.Figure,'savefig',lambda *a,**kw:None)
    rows=[dict(comparison_mode_id=f'comparison_mode_{j:02d}',kappa=k,beta_deg=b,Lambda=j+b/100,
        mapping_status='CONFIRMED',root_status='CONFIRMED',tracking_status='TRACKED') for j in range(1,7) for k in (0,1,100) for b in (0,90)]
    workflow.write_csv(tmp_path/'comparison_branches.csv',rows)
    workflow.write_csv(tmp_path/'seed_mode_mapping.csv',[dict(id=j) for j in range(18)])
    (tmp_path/'comparison_manifest.json').write_text('{}')
    before={p:p.read_bytes() for p in tmp_path.glob('*.csv')};workflow.render_comparison()
    assert all(p.read_bytes()==v for p,v in before.items())
    r=json.loads((tmp_path/'comparison_manifest.json').read_text())['render']
    assert r['matrix_calls']==r['root_calls']==r['shape_calls']==r['tracking_calls']==r['seed_matching_calls']==0


def test_protected_old_files_are_only_read_and_changes_detected(tmp_path,monkeypatch):
    oldfile=tmp_path/'old.csv';oldfile.write_text('original frequencies and branches')
    snapshot={str(oldfile):workflow.sha(oldfile)}
    (tmp_path/'run_manifest.json').write_text(json.dumps({'protected_source_files_sha256':snapshot}))
    monkeypatch.setattr(workflow,'HINGE_OUTPUT',tmp_path)
    workflow.assert_protected_sources();assert oldfile.read_text()=='original frequencies and branches'
    oldfile.write_text('changed')
    with pytest.raises(ValueError,match='Protected source changed'):workflow.assert_protected_sources()


def test_seed_guard_qualification_retains_original_warning():
    roots=[dict(Omega=float(j),root_status='CONFIRMED') for j in range(1,8)]
    warning=dict(interval=[6.8,7.1],reason='BOUNDARY_MINIMUM')
    p=dict(roots=roots,status='CONFIRMED',guard_warnings=[warning],common_class_upper=7.4)
    result=workflow.seed_guard_qualification(p)
    assert result['ROOT']=='CONFIRMED' and result['GUARD']=='QUALIFIED_DETECTOR_WARNING'
    assert result['original_helper_status']=='CONFIRMED' and result['source_warning_retained']==[warning]
    assert p['guard_warnings']==[warning]


def test_local_seed_parameter_and_combined_matrix_budget(monkeypatch):
    context=workflow.PointMatrices(10,0)
    assert context.point['k_theta']==10*ARM.D/ARM.L and context.point['beta_rad']==0
    monkeypatch.setitem(workflow.CRITERIA,'max_B_per_point',2)
    context.assembly(.1);context.block(.1,1)
    with pytest.raises(workflow.old.CostLimit):context.block(.1,-1)
    for k in (-1,np.nan,np.inf):
        with pytest.raises(ValueError):workflow.PointMatrices(k,0)


def test_import_is_missing_only_including_doublet_keys(tmp_path,monkeypatch):
    # The existing completed entries are never reconstructed, even when two
    # rows at a frequency have different cached symmetry forms.
    source=dict(contract=workflow.old.contract(),points={})
    point='k0_b90';source['points'][point]=dict(status='COMPLETED',rows=[dict(Omega=2.),dict(Omega=2.)])
    (tmp_path/'diagnostics.json').write_text(json.dumps(source));(tmp_path/'spectrum_roots.csv').write_text('point_id,Omega\n')
    monkeypatch.setattr(workflow.old,'OUTPUT',tmp_path);monkeypatch.setattr(workflow.old,'grid_tenths',lambda:[900])
    # contract grid changes with the mock, but only physical fields are checked.
    monkeypatch.setattr(workflow,'PLOT_KAPPAS',[0]);monkeypatch.setattr(workflow,'save',lambda *a:None)
    def forbidden(*a,**kw):raise AssertionError('completed group recalculated')
    monkeypatch.setattr(workflow,'PointMatrices',forbidden)
    roots=[dict(shape_key='class_plus'),dict(shape_key='class_minus')]
    state=dict(points={point:dict(roots=roots)});shapes={'class_plus__vector':np.array([1,0]),'class_minus__vector':np.array([0,1])}
    workflow.import_saved(state,shapes)
    assert state['points'][point]['roots']==roots
    assert state['source_data_audit']['source_frequency_rows']==2
