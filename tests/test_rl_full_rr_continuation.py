"""Saved-field continuation: no spectral routines in these tests."""
import numpy as np
import pytest
import copy

from scripts.analysis import postprocess_rl_full_rr_continuation as p


def synthetic(mu=0., beta=0.):
    xi=np.array([.2,.5,.8]);w=np.array([.3,.4,.3])
    sections=[];arms=[]
    for length in (1-mu,1+mu):
        sections.append(dict(length=length,rhoA=2.,rhoI=.4,J=.3))
        y=np.zeros((2,3,6));dy=np.zeros_like(y)
        y[0,:,0]=[1,2,3];y[0,:,2]=[.5,1,1.5];dy[0,:,0]=2/length
        y[1,:,1]=[1,3,2]
        arms.append(dict(x=length*xi,weights=length*w,y=y,dy=dy))
    return dict(beta=beta,mu=mu,sections=sections,arms=arms)


def test_full_inertia_contains_rotation_and_axial_gradient():
    point=synthetic();f=p.features(point)
    expected=0.
    for a,s in zip(point['arms'],point['sections']):
        expected+=np.dot(a['weights'],2*a['y'][0,:,0]**2+.4*a['y'][0,:,2]**2+.3*a['dy'][0,:,0]**2)
    assert f[0]@f[0]==pytest.approx(expected)
    assert f[0]@f[1]==pytest.approx(0)


def test_common_coordinate_derivative_transport_and_positive_metric():
    a,b=synthetic(0.),synthetic(.4)
    mid=p.metric_between(a,b)
    # Identical normalized-coordinate fields include L*U_x, not U_x alone.
    np.testing.assert_allclose(p.features(a,mid),p.features(b,mid))
    assert np.all(np.linalg.norm(p.features(a,mid),axis=1)>0)


def test_assignment_follows_shapes_through_frequency_order_change():
    a=np.eye(3);b=a[[1,0,2]]
    with p.saved_only():
        d=p.assignment(a,b,[1.,2.,3.],[1.1,1.9,3.1])
    np.testing.assert_array_equal(d.assignment,[1,0,2])
    assert np.allclose(d.mac[np.arange(3),d.assignment],1)


def test_sign_change_does_not_change_assignment():
    a=np.eye(3);b=a[[2,0,1]]
    with p.saved_only():
        d=p.assignment(a,b,[1,2,3],[1,2,3])
        e=p.assignment(-a,b*np.array([-1,1,-1])[:,None],[1,2,3],[1,2,3])
    np.testing.assert_array_equal(d.assignment,e.assignment)
    np.testing.assert_allclose(d.mac,e.mac)


def test_rotated_pair_is_valid_subspace_but_not_individual_assignment():
    a=np.eye(2);b=np.array([[1.,1.],[-1.,1.]])/np.sqrt(2)
    with p.saved_only():d=p.assignment(a,b,[1,1.00001],[1,1.00001])
    np.testing.assert_allclose(p.principal_cosines(a,b),1)
    for i,k in enumerate(d.assignment):
        assert not p.accept_transition(d.mac[i,k],0,True,0)
    # No degeneracy claim is made: two distinct frequencies stay distinct.
    assert d.sorted_lambdas==(1.,1.00001)


def test_prefix_loss_cannot_be_replaced_by_nearby_frequency():
    a=np.eye(3)[[2]];b=np.eye(3)[:2]
    with p.saved_only():d=p.assignment(a,b,[3.],[2.9,3.01])
    assert not p.accept_transition(d.mac[0,d.assignment[0]],0,True,0)


def test_endpoint_instability_or_large_jump_is_not_accepted():
    assert not p.accept_transition(.999,.5,False,0)
    assert not p.accept_transition(.999,.5,True,2)
    assert p.accept_transition(.999,.5,True,0)


def test_solver_entry_is_blocked_before_any_spectral_work():
    with p.saved_only(),pytest.raises(RuntimeError,match='Forbidden'):
        p.tracker.solve_sorted_lambdas(beta=0.,mu=0.,epsilon=.005,n_roots=1)


def test_midpoint_and_endpoint_assignments_for_small_length_change():
    a,b=synthetic(.1),synthetic(.105)
    found=[]
    with p.saved_only():
        for t in (0,.5,1):
            metric=p.metric_between(a,b,t)
            found.append(p.assignment(p.features(a,metric),p.features(b,metric),[1,2],[1.01,2.01]).assignment)
    for assignment in found:np.testing.assert_array_equal(assignment,[0,1])


def test_saved_result_regression_preserves_ambiguity_and_seed_classes():
    folder=p.ROOT/p.read(p.DEFAULT)['output']
    if not (folder/'manifest.json').exists():
        pytest.skip('Local numerical archives are not distributed with Git')
    rows=p.table(folder/'branch_records.csv');quality=p.read(folder/'quality_checks.json')
    assert quality['beta_mu001']['minimum_confirmed_seed_pair_gap']['x']==6.1
    assert quality['beta_mu01']['minimum_confirmed_seed_pair_gap']['x']==5.4
    for r in rows:
        if r['series']=='mu_beta30' and r['branch_id'] in ('B5','B6') and float(r['x'])>=.31:
            assert r['status']=='UNRESOLVED_FROM_SEED'
            assert r['current_sorted_position']==''
        if r['series']=='mu_beta30' and r['step']=='0' and r['branch_id']=='B2':
            assert r['seed_class']=='-' and r['current_sorted_position']=='1'
    assert all(v['minimum_confirmed_seed_pair_gap']['complete_interval'] for v in quality.values())
    manifest=p.read(folder/'manifest.json')
    assert manifest['spectral_calls']==manifest['basis_calls']==manifest['ritz_calls']==0
    for name,h in manifest['outputs'].items():assert p.sha(folder/name)==h


def new_archive_rows():
    path=p.MAPS/'2edc45fece76e252/manifest.json'
    if not path.exists():pytest.skip('Local spectral archives not distributed')
    return p.read(path)['points']


@pytest.mark.parametrize('name,mu',[('B0_01',0.),('B03_01',.3)])
def test_unequal_square_sections_k6_load_full_guard_and_canonical_coefficients(name,mu):
    row=next(r for r in new_archive_rows() if r['map']==name)
    with p.saved_only():point=p.Archive().point(row)
    assert len(point['frequencies'])==point['K_guard']==7
    assert point['mu']==mu and point['tau']==.1
    assert point['sections'][0]['EA']!=point['sections'][1]['EA']
    np.testing.assert_allclose(p.features(point)@p.features(point).T,np.eye(7),atol=5e-7)


@pytest.mark.parametrize('field',['J','rhoI','EI','length'])
def test_formal_metadata_match_cannot_hide_bad_section_coefficients(monkeypatch,field):
    row=new_archive_rows()[0];path=p.ROOT/row['record'];original=p.read
    bad=copy.deepcopy(original(path));bad['contract']['arms'][0][field]*=1.01
    monkeypatch.setattr(p,'read',lambda q:bad if p.Path(q)==path else original(q))
    with p.saved_only(),pytest.raises(ValueError,match='canonical physical contract'):
        p.Archive().point(row)


def test_full_metric_matches_existing_global_inertia_for_nonzero_tau():
    rows=new_archive_rows();archive=p.Archive()
    with p.saved_only():
        a,b=[archive.point(next(r for r in rows if r['map']=='B03_01' and r['x']==x)) for x in (10.,10.5)]
    def mode(point,k):return dict(arms=[{f:arm[f][k] if f in ('y','dy','ends') else arm[f] for f in ('y','dy','ends','x','weights')} for arm in point['arms']])
    model=p.make_model(**a['json']['contract']['model'])
    value=p.sections_engine.global_inner(model,mode(a,1),mode(b,2),a['beta'],b['beta'])
    assert p.features(a)[1]@p.features(b)[2]==pytest.approx(value,rel=1e-13,abs=1e-14)


@pytest.mark.parametrize('name',['basis','state_matrix'])
def test_geometry_allowlist_never_permits_physical_basis_or_matrix(name):
    with p.saved_only(),pytest.raises(RuntimeError,match='Forbidden physical evaluation'):
        getattr(p.sections_engine.p.base,name)(None,None,None,None) if name=='basis' else getattr(p.sections_engine.p.base,name)(None,None)


def test_legacy_k8_saved_archive_remains_readable():
    path=p.MAPS/'b775d0980d7bad77/manifest.json'
    if not path.exists():pytest.skip('Local spectral archives not distributed')
    row=next(r for r in p.read(path)['points'] if r['map']=='square_beta_mu01_tau0_refined' and r['x']==15.)
    with p.saved_only():point=p.Archive().point(row)
    assert point['K_guard']==len(point['frequencies'])==9 and point['tau']==0


def test_near_pair_current_grid_does_not_silently_confirm_low_overlap():
    rows=new_archive_rows();archive=p.Archive()
    with p.saved_only():
        a,b=[archive.point(next(r for r in rows if r['map']=='B03_01' and r['x']==x)) for x in (10.,10.5)]
        diag=p.assignment(p.features(a)[:6],p.features(b),a['frequencies'][:6],b['frequencies'])
    for i in (1,2):
        k=diag.assignment[i];overlap=diag.mac[i,k];margin=overlap-np.max(np.delete(diag.mac[i],k))
        assert .947<overlap<.95 and not p.accept_transition(overlap,margin,True,abs(k-i))
    assert p.principal_cosines(p.features(a)[[1,2]],p.features(b)[[1,2]]).min()>.99


def test_asymmetric_saved_regression_does_not_recover_lost_seed_identity():
    folder=p.MAPS/'_asymmetric_beta_publication_20261109_final'
    if not (folder/'manifest.json').exists():
        pytest.skip('Local derived results not distributed')
    rows=p.table(folder/'branch_records.csv');transitions=p.table(folder/'transitions.csv')
    assert len(rows)==2*261*6 and len(transitions)==2*260*6
    failed=[t for t in transitions if t['local_status']!='ACCEPTED']
    assert {(t['branch_id'],float(t['x_from']),float(t['x_to'])) for t in failed}=={('B2',10.,10.5),('B3',10.,10.5)}
    assert all(t['reverse_stable']=='True' for t in transitions)
    for r in rows:
        lost=r['series']=='beta_mu03_tau01' and r['branch_id'] in ('B2','B3') and float(r['x'])>=10.5
        assert (r['status']=='UNRESOLVED_FROM_SEED')==lost
        assert (r['current_sorted_position']=='')==lost
    manifest=p.read(folder/'manifest.json')
    assert manifest['spectral_calls']==manifest['basis_calls']==manifest['ritz_calls']==0
    for name,h in manifest['outputs'].items():assert p.sha(folder/name)==h


def test_shape_panels_follow_confirmed_descendant_instead_of_old_sorted_slot():
    rows=[dict(series='s',beta_deg=10.5,seed_position=2,status='CONFIRMED_FROM_SEED',current_sorted_position=3,branch_id='B2'),
          dict(series='s',beta_deg=10.5,seed_position=3,status='UNRESOLVED_FROM_SEED',current_sorted_position=None,branch_id='B3')]
    k,record=p.displayed_mode(rows,'s',10.5,2)
    assert k==2 and record['branch_id']=='B2'
    k,record=p.displayed_mode(rows,'s',10.5,3)
    assert k==2 and record is None


def test_refined_saved_pair_keeps_order_without_forcing_coarse_assignment():
    path=p.MAPS/'73b65e9a33d121c7/manifest.json'
    if not path.exists():pytest.skip('Local refined spectra not distributed')
    cfg=p.read(p.ROOT/'private/dissertation/research/chapter1/chapter1_asymmetric_beta_publication_refined.json')
    archive=p.Archive();spec=cfg['series'][1]
    with p.saved_only():
        points=[archive.point(r) for r in p.read(path)['points'] if 10<=r['x']<=10.5]
        rows,tr,_,_=p.continue_series(spec,points,[(f'B{k}','',0) for k in range(1,7)],cfg)
        coarse=p.assignment(p.features(points[0])[:6],p.features(points[-1]),points[0]['frequencies'][:6],points[-1]['frequencies'])
    assert len(points)==11
    assert all(t['local_status']=='ACCEPTED' and t['reverse_stable'] for t in tr)
    assert .958<min(t['overlap'] for t in tr)<.959
    assert all(r['current_sorted_position']==r['seed_position'] for r in rows)
    np.testing.assert_array_equal(coarse.assignment[1:3],[2,1])


def test_refined_delivery_preserves_identity_and_monochrome_mode_axes():
    folder=p.MAPS/'_asymmetric_beta_refined_20261109'
    if not (folder/'manifest.json').exists():pytest.skip('Local refined postprocessing not distributed')
    rows=p.table(folder/'branch_records.csv')
    assert len(rows)==(261+270)*6
    assert all(r['status']=='CONFIRMED_FROM_SEED' for r in rows)
    assert all(r['current_sorted_position']==r['seed_position'] for r in rows)
    steps=p.read(folder/'step_refinement.json')
    assert steps['forward_reverse_all_nodes_agree'] and steps['reverse_rejected_transitions']==0
    modes=p.table(folder/'plotted_mode_panels.csv')
    assert len(modes)==12 and {r['deformed_axis_color'] for r in modes}=={'.15'}
    assert {r['deformed_axis_linestyle'] for r in modes}=={'solid'}
    assert all(r['canonical_branch'] for r in modes)


def redraw_fixture():
    point=dict(beta=0.,mu=.1,frequencies=np.arange(1.,10.),record='point.json',npz='shape.npz')
    rows=[dict(series='s',step='0',seed_position=str(k),current_sorted_position=str(k),
               status='CONFIRMED_FROM_SEED',beta_deg='0',mu='.1',Lambda=str(float(k)),
               sign=str(float((-1)**k)),source_record='point.json',source_npz='shape.npz') for k in range(1,9)]
    # Saved identities need not equal sorted positions: retain an accepted permutation.
    rows[1]['current_sorted_position'],rows[1]['Lambda']='3','3.0'
    rows[2]['current_sorted_position'],rows[2]['Lambda']='2','2.0'
    return point,rows


def test_selected_redraw_keeps_six_of_eight_saved_identities_and_phases():
    point,rows=redraw_fixture()
    selected=p.confirmed_redraw_rows(rows,dict(id='s',branches=6),[point])
    assert len(selected)==6 and len(rows)==8 and len(point['frequencies'])==9
    assert [r['current_sorted_position'] for r in selected]==[1,3,2,4,5,6]
    assert [r['sign'] for r in selected]==[-1.,1.,-1.,1.,-1.,1.]
    assert rows[1]['current_sorted_position']=='3'  # Inputs stay unmodified.


@pytest.mark.parametrize('field,value',[
    ('status','UNRESOLVED_FROM_SEED'),('source_record','different.json'),
    ('source_npz','different.npz'),('Lambda','3.01'),('current_sorted_position','7')])
def test_selected_redraw_rejects_unconfirmed_or_mismatched_saved_data(field,value):
    point,rows=redraw_fixture();rows[1][field]=value
    with pytest.raises(ValueError):p.confirmed_redraw_rows(rows,dict(id='s',branches=6),[point])


def test_selected_redraw_rejects_missing_and_duplicate_records():
    point,rows=redraw_fixture()
    for incomplete in (rows[1:],rows+[rows[0]]):
        with pytest.raises(ValueError):p.confirmed_redraw_rows(incomplete,dict(id='s',branches=6),[point])
