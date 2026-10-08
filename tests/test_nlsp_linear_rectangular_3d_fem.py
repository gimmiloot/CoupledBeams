"""FEM-1 contract tests; no executable mesh/solid-solver jobs."""
import math
from fractions import Fraction
from pathlib import Path
from itertools import permutations

import numpy as np
import pytest

from scripts.analysis import verify_nlsp_linear_rectangular_3d_fem as workflow
from scripts.analysis import solid_fem_single_rod_fixed_fixed as historical_solid


@pytest.fixture(scope='module')
def preflight():
    # One cheap existing analytical preflight; no Gmsh, CCX or 3D eigensolve.
    return workflow.build_preflight()


def test_preflight_geometry_material_and_orientation(preflight):
    result,_=preflight
    assert result['geometry']=={'L':1.,'b':.2,'h':.1}
    assert result['material']=={'E':1.,'rho':1.,'nu':.3,'kappa':5/6}
    b,h=Fraction(1,5),Fraction(1,10)
    expected={'A0':b*h,'I_parallel':b*h**3/12,'I_perp':h*b**3/12,
              'Ip':b*h*(b*b+h*h)/12}
    for key,value in expected.items():
        assert result['section'][key]==pytest.approx(float(value),rel=1e-14)
    assert result['section']['I_perp']/result['section']['I_parallel']==pytest.approx(4.)


def test_original_section_material_coefficients(preflight):
    result,_=preflight
    area=Fraction(1,50)
    ip=Fraction(1,60000)
    ib=Fraction(1,15000)
    G=Fraction(5,13)
    expected={'m':area,'jp':ip,'jb':ib,'C':area/Fraction(91,100),
              'H':Fraction(5,6)*G*ip,'S':Fraction(5,6)*G*area,'Bp':ip,'Bb':ib}
    for key,value in expected.items():
        assert result['coefficients'][key]==pytest.approx(float(value),rel=1e-14)
    assert result['coefficient_mapping_error']<1e-14
    assert result['outplane_yartsev_state_mapping_max_error']<1e-14


def test_generalized_torsion_has_source_and_is_not_polar_substitution(preflight):
    result,_=preflight
    provenance=result['C_T_provenance']
    assert 'generalized_torsional_stiffness' in provenance['source']
    assert provenance['book_geometry']=={'a':.2,'b':.1,'length':1.}
    assert provenance['Sbar16']==0.
    assert provenance['C_T']==provenance['Cbar']
    assert provenance['C_T']==pytest.approx(1.759089824002232e-5,rel=2e-11)
    assert provenance['estimated_relative_tail']<=1e-12
    assert provenance['terms_used']>0
    assert provenance['substituted_G_Ip'] is False
    assert provenance['independent_dynamic_warping_field'] is False
    gip=1/2.6*result['section']['Ip']
    assert abs(provenance['C_T']/gip-1)>.4


def test_seven_essential_end_values_are_not_slope_constraints(preflight):
    result,profiles=preflight
    bc=result['boundary_conditions']
    assert bc['essential_fields']==['u','w','v','Phi','psi','theta','c']
    assert bc['endpoints']==[0.,1.]
    assert bc['rotation_clamp']=='section_rotation'
    assert bc['source_book_slope_clamp_used'] is False
    for profile in profiles.values():
        q=profile['q']
        scale=np.max(np.abs(q),axis=0)
        active=scale>0
        assert np.max(np.abs(q[[0,-1]][:,active])/scale[active])<1e-9


def test_full_merged_spectrum_has_all_four_families_and_local_identity(preflight):
    result,_=preflight
    rows=result['merged_spectrum']
    expected=[('inplane_bending',1,.6054167303477958),
        ('outplane_bending',1,1.038923588954992),('torsion',1,1.443392697870566),
        ('inplane_bending',2,1.551535482991772),('outplane_bending',2,2.3781017311578183),
        ('inplane_bending',3,2.8042758204627085),('torsion',2,2.886785395741132),
        ('axial_mh',1,3.150123638896353)]
    assert len(rows)==len(expected)
    assert [row['sorted_index'] for row in rows]==list(range(1,9))
    for row,(family,local,omega) in zip(rows,expected):
        assert (row['family'],row['local_mode'])==(family,local)
        assert row['omega']==pytest.approx(omega,rel=1e-10)
    assert len({(row['family'],row['local_mode']) for row in rows})==8
    assert result['first_axial_sorted_index']==8
    assert result['geometry_gate']=='PASS'


def test_completeness_is_saturated_independent_bound(preflight):
    result,_=preflight
    completeness=result['completeness']
    assert completeness['status']=='PASS'
    assert completeness['total_found']==completeness['total_upper_count']==8
    for certificate in completeness['certificates'].values():
        assert certificate['saturated'] is True
        assert certificate['found_distinct_modes']==certificate['upper_count']
    window=result['frequency_window']
    assert window['counts_by_family']=={'axial_mh':1,'inplane_bending':3,'outplane_bending':2,'torsion':2}
    assert window['guard_multiplier']==1.1
    assert window['omega_max']==pytest.approx(1.1*window['omega_required_max'])
    assert window['initial_requested_3d_eigenpairs']>=24
    assert result['contraction_dominated_window']=='OUTSIDE_BOUNDED_WINDOW'
    assert result['contraction_optical_cutoff_omega']>window['omega_max']


def test_acoustic_mh_is_independent_contraction_coordinate(preflight):
    result,profiles=preflight
    row=next(row for row in result['merged_spectrum'] if row['family']=='axial_mh')
    assert row['branch_character']=='acoustic_axial'
    assert row['kinetic_primary_fraction']>.99
    assert 0<row['kinetic_secondary_fraction']<.01
    q=profiles['axial_mh:1']['q']
    assert np.max(np.abs(q[:,0]))>0
    assert np.max(np.abs(q[:,6]))>0
    assert np.max(np.abs(q[:,[1,2,3,4,5]]))==0
    # Pure symmetry makes the first acoustic u even and contraction c odd.
    assert np.max(np.abs(q[:,0]-q[::-1,0]))<1e-10
    assert np.max(np.abs(q[:,6]+q[::-1,6]))<1e-9


def test_profiles_keep_fields_and_analytic_boundary_energy_evidence(preflight):
    result,profiles=preflight
    allowed={'axial_mh':{0,6},'inplane_bending':{1,5},'outplane_bending':{2,4},'torsion':{3}}
    for row in result['merged_spectrum']:
        profile=profiles[f"{row['family']}:{row['local_mode']}"]
        x,q=profile['x'],profile['q']
        assert x.shape==(401,)
        assert q.shape==(401,7)
        assert np.all(np.isfinite(q))
        assert np.all(np.diff(x)>0)
        assert x[0]==0 and x[-1]==1
        forbidden=sorted(set(range(7))-allowed[row['family']])
        assert np.max(np.abs(q[:,forbidden]))==0
        if row['family']!='torsion':
            diagnostic=row['diagnostics']
            for name in ('boundary_scaled_residual','equation_scaled_residual','energy_relative_error',
                         'analytic_singular_ratio','independent_state_shooting_singular_ratio'):
                assert diagnostic[name]<1e-8


def test_normalized_omega_and_cycles_per_time_are_separate(preflight):
    result,_=preflight
    for row in result['merged_spectrum']:
        assert row['frequency_cycles_per_time']*2*math.pi==pytest.approx(row['omega'],rel=1e-14)
    window=result['frequency_window']
    assert window['frequency_max_cycles_per_time']*2*math.pi==pytest.approx(window['omega_max'],rel=1e-14)


@pytest.mark.parametrize('geometry',[{'L':1.,'b':.2,'h':.11},{'L':1.,'b':.4,'h':.1}])
def test_no_free_geometry_tuning(geometry):
    with pytest.raises(ValueError,match='authorized'):
        workflow.build_preflight({'geometry':geometry})




def synthetic_c3d10_box(path, *, invert_first=False, solid_type='C3D10', shift=0.):
    L,h,b=1.,.1,.2
    nodes,elements,lookup={}, {}, {}
    def node_id(point):
        key=tuple(np.round(point,12))
        if key not in lookup:
            number=len(nodes)+1
            lookup[key]=number
            nodes[number]=tuple(point)
        return lookup[key]
    for permutation in permutations(range(3)):
        vertices=[np.zeros(3)]
        point=np.zeros(3)
        for axis in permutation:
            point=point.copy()
            point[axis]=1
            vertices.append(point)
        xyz=np.asarray(vertices)*[L,h,b]+[shift,-h/2,-b/2]
        if np.linalg.det((xyz[1:]-xyz[:1]).T)<0:
            xyz[[1,2]]=xyz[[2,1]]
        if invert_first and not elements:
            xyz[[1,2]]=xyz[[2,1]]
        connectivity=[node_id(point) for point in xyz]
        connectivity += [node_id((xyz[i]+xyz[j])/2) for i,j in workflow.TET10_EDGES]
        elements[len(elements)+1]=connectivity
    lines=['*NODE']+[f'{number}, '+', '.join(str(value) for value in point) for number,point in nodes.items()]
    lines += [f'*ELEMENT,TYPE={solid_type},ELSET=SOLID']
    lines += [f'{number}, '+', '.join(str(value) for value in connectivity) for number,connectivity in elements.items()]
    path.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return path


def test_positive_tet10_quadrature_is_exact_through_degree_five():
    bary,weights=workflow.tet10_quadrature()
    assert np.all(weights>0)
    assert weights.sum()==pytest.approx(1/6,abs=1e-14)
    for a in range(6):
        for b in range(6-a):
            for c in range(6-a-b):
                exact=math.factorial(a)*math.factorial(b)*math.factorial(c)/math.factorial(a+b+c+3)
                actual=np.sum(weights*bary[:,1]**a*bary[:,2]**b*bary[:,3]**c)
                assert actual==pytest.approx(exact,abs=1e-14)
    values,gradient=workflow.tet10_shape(bary)
    np.testing.assert_allclose(values.sum(axis=1),1.,atol=1e-14)
    np.testing.assert_allclose(gradient.sum(axis=1),0.,atol=1e-14)


def test_synthetic_box_mesh_has_real_fixed_faces_and_volume_metadata(tmp_path):
    path=synthetic_c3d10_box(tmp_path/'box.inp')
    metadata,mesh=workflow.audit_rectangular_mesh(path,1.,.2,.1,1.,historical_solid.read_gmsh_inp_mesh_data)
    assert metadata['status']=='PASS'
    assert metadata['nodes']==27
    assert metadata['c3d10_elements']==6
    assert metadata['solid_element_types']==['C3D10']
    assert metadata['volume_face_connected_components']==1
    assert metadata['fixed_left_count']==metadata['fixed_right_count']==9
    assert set(metadata['fixed_left_ids']).isdisjoint(metadata['fixed_right_ids'])
    assert all(mesh.nodes[number][0]==0. for number in metadata['fixed_left_ids'])
    assert all(mesh.nodes[number][0]==1. for number in metadata['fixed_right_ids'])
    assert metadata['volume']==pytest.approx(.02,abs=1e-14)
    assert metadata['mass']==pytest.approx(.02,abs=1e-14)
    assert metadata['rigid_body_kinematic_constraint_check']=='PASS'
    assert metadata['negative_or_zero_jacobian_elements']==0
    assert metadata['minimum_quadratic_jacobian']>0
    assert metadata['straight_midpoint_geometry'] is True
    assert metadata['quality_min_median_max'][0]>0
    for direction in ('thickness_h','width_b'):
        edges=metadata['actual_resolution'][direction]['box_edge_resolution']['edges']
        assert all(edge['linear_segment_count']==1 for edge in edges)
        assert all(edge['corner_vertex_count']==2 for edge in edges)


def test_quadrature_mass_integrates_quadratic_displacement_product(tmp_path):
    path=synthetic_c3d10_box(tmp_path/'mass_box.inp')
    _,mesh=workflow.audit_rectangular_mesh(path,1.,.2,.1,1.,historical_solid.read_gmsh_inp_mesh_data)
    samples=workflow.quadrature_arrays(mesh)
    assert samples['weights'].sum()==pytest.approx(.02,abs=1e-14)
    displacement=np.asarray([mesh.nodes[int(number)] for number in samples['node_ids']])**2
    field=np.einsum('qi,eic->eqc',samples['N'],displacement[samples['conn']])
    norm=np.sum(samples['weights']*np.sum(field**2,axis=2))
    L,h,b=1.,.1,.2
    exact=h*b*L**5/5+L*b*h**5/80+L*h*b**5/80
    assert norm==pytest.approx(exact,abs=1e-14)


def test_inverted_c3d10_is_blocked_before_eigenanalysis(tmp_path):
    path=synthetic_c3d10_box(tmp_path/'inverted.inp',invert_first=True)
    metadata,mesh=workflow.audit_rectangular_mesh(path,1.,.2,.1,1.,historical_solid.read_gmsh_inp_mesh_data)
    assert metadata['status']=='FAIL'
    assert metadata['negative_or_zero_jacobian_elements']==1
    assert any('Jacobian' in failure for failure in metadata['failures'])
    with pytest.raises(ValueError,match='Nonpositive'):
        workflow.quadrature_arrays(mesh)


def test_other_solid_type_cannot_be_relabelled_c3d10(tmp_path):
    path=synthetic_c3d10_box(tmp_path/'wrong_type.inp',solid_type='C3D4')
    with pytest.raises(ValueError,match='only C3D10'):
        workflow.audit_rectangular_mesh(path,1.,.2,.1,1.,historical_solid.read_gmsh_inp_mesh_data)


def test_missing_fixed_face_is_explicit_failure(tmp_path):
    path=synthetic_c3d10_box(tmp_path/'shifted.inp',shift=.1)
    metadata,_=workflow.audit_rectangular_mesh(path,1.,.2,.1,1.,historical_solid.read_gmsh_inp_mesh_data)
    assert metadata['status']=='FAIL'
    assert metadata['fixed_left_count']==0
    assert metadata['bbox_matches'] is False
    assert metadata['rigid_body_kinematic_constraint_check']=='FAIL'


def test_geo_contract_is_monolithic_box_with_two_physical_faces():
    text=workflow.rectangular_geo(1.,.2,.1,.05)
    assert 'Box(1)' in text
    assert 'Physical Volume("SOLID",1)' in text
    assert 'Physical Surface("FIXED_LEFT",2)' in text
    assert 'Physical Surface("FIXED_RIGHT",3)' in text
    assert 'Mesh.ElementOrder = 2' in text
    assert 'Cylinder(' not in text
    assert 'RIGID BODY' not in text

def test_quadratic_shared_faces_must_share_the_same_midside_nodes(tmp_path):
    path=synthetic_c3d10_box(tmp_path/'nonconforming.inp')
    original=historical_solid.read_gmsh_inp_mesh_data(path)
    nodes=dict(original.nodes)
    elements={number:list(connectivity) for number,connectivity in original.solid_elements.items()}
    first=elements[min(elements)]
    old_mid=first[4]
    new_mid=max(nodes)+1
    nodes[new_mid]=nodes[old_mid]
    first[4]=new_mid
    lines=['*NODE']+[f'{number}, '+', '.join(str(value) for value in point) for number,point in nodes.items()]
    lines += ['*ELEMENT,TYPE=C3D10,ELSET=SOLID']
    lines += [f'{number}, '+', '.join(str(value) for value in connectivity) for number,connectivity in elements.items()]
    path.write_text('\n'.join(lines)+'\n',encoding='utf-8')
    metadata,_=workflow.audit_rectangular_mesh(path,1.,.2,.1,1.,historical_solid.read_gmsh_inp_mesh_data)
    assert metadata['status']=='FAIL'
    assert metadata['nonconforming_quadratic_face_count']>0
    assert metadata['negative_or_zero_jacobian_elements']==0


def test_gmsh_warnings_are_retained_and_critical_errors_block_mesh(tmp_path):
    path=synthetic_c3d10_box(tmp_path/'warnings.inp')
    warning='Warning : example numerical warning retained for review'
    metadata,_=workflow.audit_rectangular_mesh(path,1.,.2,.1,1.,historical_solid.read_gmsh_inp_mesh_data,warning)
    assert metadata['status']=='PASS'
    assert metadata['gmsh_warnings_errors']==[warning]
    error='Error : invalid element encountered'
    metadata,_=workflow.audit_rectangular_mesh(path,1.,.2,.1,1.,historical_solid.read_gmsh_inp_mesh_data,error)
    assert metadata['status']=='FAIL'
    assert metadata['gmsh_warnings_errors']==[error]


def sample_box(nbin=41):
    x=np.array([(k+t)/nbin for k in range(nbin) for t in (.2,.5,.8)])
    y=np.array([-.04,0.,.04]); z=np.array([-.08,0.,.08])
    xyz=np.array(np.meshgrid(x,y,z,indexing='ij')).reshape(3,-1).T
    return xyz,np.ones(len(xyz))*.02/len(xyz)


def test_coordinate_signs_and_contraction_not_cartesian_dof():
    x=np.array([0.,1.]);q=np.tile(np.arange(1.,8.),(2,1))
    actual=workflow.nlsp_lift_fields(x,q,np.array([[.5,.02,.04]]))
    np.testing.assert_allclose(actual,[[1.32,-2.02,-2.92]],rtol=0,atol=1e-14)
    purec=np.zeros((2,7));purec[:,6]=.1
    actual=workflow.nlsp_lift_fields(x,purec,np.array([[.5,.02,.04]]))
    np.testing.assert_allclose(actual,[[0.,.002,0.]],atol=1e-15)


def test_rigid_section_fit_recovers_rotation_signs_and_x_variation():
    xyz,w=sample_box()
    x=np.array([0.,1.]);q=np.array([[1.,2.,3.,4.,5.,6.,0.],[2.,3.,4.,5.,6.,7.,0.]])
    u=workflow.nlsp_lift_fields(x,q,xyz)
    projected=workflow.nlsp_project_sections(xyz,u,w,1.,.02,.02*.1**2/12,.02*.2**2/12)
    target=np.array([np.interp(projected['x'][1:-1],x,q[:,k]) for k in range(7)]).T
    np.testing.assert_allclose(projected['fields'][1:-1],target,atol=1e-11)
    assert projected['residual_fraction']<1e-25


def test_contraction_is_explicit_diagnostic_and_axial_structure_remains():
    xyz,w=sample_box()
    x=np.array([0.,1.]);q=np.zeros((2,7));q[:,0]=1.;q[:,6]=-.3
    u=workflow.nlsp_lift_fields(x,q,xyz)
    p=workflow.nlsp_project_sections(xyz,u,w,1.,.02,.02*.1**2/12,.02*.2**2/12)
    np.testing.assert_allclose(p['fields'][1:-1,6],-.3,atol=1e-12)
    assert p['family']=='axial_mh'
    assert p['contraction_status']=='DIAGNOSTIC_EFFECTIVE_THICKNESS_STRAIN_NOT_MH_DOF'
    assert p['transverse_residual_fraction']>0


def test_local_section_mode_is_not_forced_into_four_families():
    xyz,w=sample_box();local=workflow.nlsp_local_vectors(xyz)
    u=np.zeros_like(xyz);u[:,1]=local[:,1]**2-np.mean(local[:,1]**2)
    u=workflow.nlsp_local_vectors(u)
    p=workflow.nlsp_project_sections(xyz,u,w,1.,.02,.02*.1**2/12,.02*.2**2/12)
    assert p['family']=='cross_section_or_local'
    assert p['residual_fraction']>.99


def test_torsion_and_warping_are_distinct_diagnostics():
    xyz,w=sample_box();local=workflow.nlsp_local_vectors(xyz)
    x=np.array([0.,1.]);q=np.zeros((2,7));q[:,3]=1.
    u=workflow.nlsp_lift_fields(x,q,xyz);u[:,0]+=3*local[:,1]*local[:,2]
    p=workflow.nlsp_project_sections(xyz,u,w,1.,.02,.02*.1**2/12,.02*.2**2/12)
    assert p['family']=='torsion'
    assert 0<p['axial_warp_fraction']<.35
    np.testing.assert_allclose(p['fields'][1:-1,3],1.,atol=1e-12)


def test_mac_uses_complete_vector_and_mass_weights():
    a=np.array([[1.,0.,0.],[0.,2.,0.]]);b=np.array([[1.,0.,0.],[0.,-2.,0.]])
    assert workflow.nlsp_weighted_mac(a,-2*a,np.array([1.,3.]))==pytest.approx(1.)
    assert workflow.nlsp_weighted_mac(a,b,np.array([1.,3.]))==pytest.approx(121/169)
    assert workflow.nlsp_weighted_mac(a,b,np.array([1.,1.]))==pytest.approx(9/25)
    with pytest.raises(ValueError):workflow.nlsp_weighted_mac(a,b,np.array([1.,-1.]))


def test_one_to_one_shape_assignment_and_visible_duplicate_conflicts():
    assigned=workflow.nlsp_shape_assignment([[.01,.99],[.98,.02]])
    assert [r['column'] for r in assigned]==[1,0]
    assert all(r['status']=='MATCHED' for r in assigned)
    duplicate=workflow.nlsp_shape_assignment([[.95,.01],[.94,.02]])
    assert all(r['conflict'] and r['status']=='DUPLICATE_INDEPENDENT_MATCH' for r in duplicate)
    assert len({r['column'] for r in duplicate})==2


def test_ambiguous_subspaces_have_no_asserted_individual_identity():
    assigned=workflow.nlsp_shape_assignment([[.51,.49],[.49,.51]])
    assert all(r['ambiguous'] and r['status'] in ('LOW_MAC','AMBIGUOUS_SHAPE_OR_SUBSPACE') for r in assigned)
    a=np.array([[[1.,0.,0.],[0.,0.,0.]],[[0.,0.,0.],[1.,0.,0.]]])
    b=np.array([(a[0]+a[1])/np.sqrt(2),(a[0]-a[1])/np.sqrt(2)])
    p=workflow.nlsp_subspace_overlap(a,b,np.ones(2))
    np.testing.assert_allclose(p['principal_mac'],[1.,1.],atol=1e-14)
    assert p['status']=='DIAGNOSTIC_SUBSPACE_ONLY'


def test_missing_and_duplicate_fem_modes_cannot_be_silent():
    full={1:{3:(1.,0.,0.),5:(0.,1.,0.)}}
    actual=workflow.nlsp_validate_nodal_modes([3,5],full,[1])
    np.testing.assert_equal(actual[1],[[1.,0.,0.],[0.,1.,0.]])
    with pytest.raises(ValueError,match='Missing FEM'):workflow.nlsp_validate_nodal_modes([3,5],full,[1,2])
    with pytest.raises(ValueError,match='Duplicate expected'):workflow.nlsp_validate_nodal_modes([3,5],full,[1,1])
    with pytest.raises(ValueError,match='missing_nodes'):workflow.nlsp_validate_nodal_modes([3,5,9],full,[1])
    assert workflow.nlsp_shape_assignment(np.empty((2,0)))[0]['status']=='MISSING_FEM_MODE'


def test_model_difference_preserves_3d_denominator():
    result=workflow.nlsp_relative_frequency_difference(3.,4.)
    assert result['signed_relative_difference']==-.25
    assert result['absolute_relative_difference']==.25



def test_fixed_width_frd_keeps_adjacent_node_id_and_positive_ux(tmp_path):
    path=tmp_path/'sample.frd'
    def row(k,u):return ' -1'+f'{k:10d}'+''.join(f'{v:12.5E}' for v in u)
    text='    1PMODE   1\n -4  DISP  4 1\n'+row(74,(.626111,-4.26211,-3.77924))+'\n'+row(1,(0.,0.,0.))+'\n -3\n'
    path.write_text(text)
    modes=workflow.nlsp_read_frd_modes(path)
    assert sorted(modes[1])==[1,74]
    np.testing.assert_allclose(modes[1][74],[.626111,-4.26211,-3.77924])
    path.write_text(text.replace(' -3\n',row(74,(0.,0.,0.))+'\n -3\n'))
    with pytest.raises(ValueError,match='Duplicate node'):workflow.nlsp_read_frd_modes(path)
    path.write_text(text+text)
    with pytest.raises(ValueError,match='Duplicate FEM mode'):workflow.nlsp_read_frd_modes(path)


def test_frd_three_digit_exponent_overflow_has_all_cartesian_components(tmp_path):
    path=tmp_path/'three_digit.frd'
    raw=' -1        746.26111E-001-4.26211E+000-3.77924E+000'
    path.write_text('    1PMODE   1\n -4  DISP  4 1\n'+raw+'\n -3\n',encoding='utf-8')
    modes=workflow.nlsp_read_frd_modes(path)
    assert sorted(modes[1])==[74]
    np.testing.assert_allclose(modes[1][74],[.626111,-4.26211,-3.77924])


def test_dat_frequency_columns_preserve_eigenvalue_omega_and_cycles_per_time(tmp_path):
    path=tmp_path/'units.dat'
    omega=2.
    cycles=omega/(2*math.pi)
    path.write_text(' E I G E N V A L U E   O U T P U T\n MODE FREQUENCY\n'+
                    f' 1 4.000000000D+00 2.000000000D+00 {cycles:.12E} 0.0\n'+
                    ' P A R T I C I P A T I O N\n',encoding='utf-8')
    rows=workflow.parse_calculix_frequency_table(path)
    assert len(rows)==1
    assert rows[0]['raw_mode_number']==1
    assert rows[0]['eigenvalue']==pytest.approx(omega**2)
    assert rows[0]['angular_frequency']==pytest.approx(omega)
    assert rows[0]['cyclic_frequency']*2*math.pi==pytest.approx(omega,rel=1e-11)


def test_cache_identity_covers_config_model_and_executable_contents(tmp_path,monkeypatch):
    import json
    monkeypatch.setattr(workflow,'ROOT',tmp_path)
    gmsh=tmp_path/'gmsh.bin'
    ccx=tmp_path/'ccx.bin'
    model=tmp_path/'model.py'
    gmsh.write_bytes(b'gmsh synthetic identity input')
    ccx.write_bytes(b'ccx synthetic identity input')
    model.write_text('frozen synthetic source',encoding='utf-8')
    config={'schema':'identity-test-only','geometry':{'L':1.,'b':.2,'h':.1},
            'gmsh_exe':str(gmsh),'ccx_exe':str(ccx),'model_hashes':{'model.py':'synthetic-source-contract'}}
    path=tmp_path/'config.json'
    path.write_text(json.dumps(config),encoding='utf-8')
    baseline,item=workflow.identity(path)
    repeated,repeated_item=workflow.identity(path)
    assert repeated==baseline and repeated_item==item
    model.write_text('changed synthetic source',encoding='utf-8')
    changed_source,_=workflow.identity(path)
    assert changed_source!=baseline
    model.write_text('frozen synthetic source',encoding='utf-8')
    gmsh.write_bytes(b'changed synthetic executable contents')
    changed_exe,_=workflow.identity(path)
    assert changed_exe!=baseline
    gmsh.write_bytes(b'gmsh synthetic identity input')
    config['geometry']['h']=.12
    path.write_text(json.dumps(config),encoding='utf-8')
    changed_config,_=workflow.identity(path)
    assert changed_config!=baseline


def make_synthetic_cached_bundle(path):
    path.mkdir()
    summary={'statuses':{'NLSP_FEM1_3D_MODAL_EXECUTION':'NOT_RUN'},'preflight':{},'comparisons':[],
             'scope':'Synthetic unit-test cache, no scientific frequencies'}
    workflow.write_json(path/'summary.json',summary)
    item={'schema':'synthetic-cache-unit-test'}
    manifest={'identity':item,'artifact_hashes':{'summary.json':workflow.sha(path/'summary.json')}}
    workflow.write_json(path/'manifest.json',manifest)
    return item,summary


def test_cache_rejects_identity_mismatch_and_changed_artifact(tmp_path):
    bundle=tmp_path/'cache'
    item,summary=make_synthetic_cached_bundle(bundle)
    assert workflow.validate_cache(bundle,item)==summary
    with pytest.raises(ValueError,match='Cache identity mismatch'):
        workflow.validate_cache(bundle,{'schema':'other-identity'})
    workflow.write_json(bundle/'summary.json',{'altered':'synthetic data'})
    with pytest.raises(ValueError,match='Cache artifact hash mismatch'):
        workflow.validate_cache(bundle,item)


@pytest.mark.parametrize('mode',('--report-only','--plot-only'))
def test_cached_reports_and_plots_make_zero_solver_and_preflight_calls(tmp_path,monkeypatch,mode):
    bundle=tmp_path/'replay'
    _,expected=make_synthetic_cached_bundle(bundle)
    def forbidden(*args,**kwargs):
        raise AssertionError('Read-only cache replay attempted new numerical computation')
    for name in ('run_job','run_fem','build_preflight','write_preflight','identity'):
        monkeypatch.setattr(workflow,name,forbidden)
    monkeypatch.setattr(workflow.single,'generate_mesh_with_gmsh_cli',forbidden)
    monkeypatch.setattr(workflow.single,'run_calculix',forbidden)
    monkeypatch.setattr(workflow.mh,'finite_roots',forbidden)
    monkeypatch.setattr(np.linalg,'eig',forbidden)
    monkeypatch.setattr(np.linalg,'eigh',forbidden)
    assert workflow.main([mode,str(bundle)])==expected


def test_cached_preflight_recovers_exact_saved_seven_field_profiles(tmp_path,preflight):
    result,profiles=preflight
    arrays={key.replace(':','_')+'_'+name:value for key,profile in profiles.items()
            for name,value in profile.items() if isinstance(value,np.ndarray)}
    np.savez_compressed(tmp_path/'one_D_profiles.npz',**arrays)
    recovered=workflow.load_profiles(tmp_path,result)
    assert set(recovered)==set(profiles)
    for key in profiles:
        np.testing.assert_array_equal(recovered[key]['x'],profiles[key]['x'])
        np.testing.assert_array_equal(recovered[key]['q'],profiles[key]['q'])


@pytest.fixture(scope='module')
def actual_fem_evidence():
    # Optional local artifact checks. Missing ignored outputs never trigger FEM.
    bundle=Path(__file__).resolve().parents[1]/'results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad'
    if not (bundle/'summary.json').exists():
        pytest.skip('Bounded actual FEM-1 evidence is not present in this checkout')
    summary=workflow.validate_cache(bundle)
    return bundle,summary


def test_actual_pre_fem_geometry_and_frozen_config_are_preserved(actual_fem_evidence):
    bundle,summary=actual_fem_evidence
    config=workflow.read_json(bundle/'frozen_config.json')
    preflight=workflow.read_json(bundle/'preflight.json')
    manifest=workflow.read_json(bundle/'manifest.json')
    assert config==summary['config']==manifest['identity']['config']
    assert preflight==summary['preflight']
    assert config['geometry']==preflight['geometry']=={'L':1.,'b':.2,'h':.1}
    assert config['material']==preflight['material']=={'E':1.,'rho':1.,'nu':.3,'kappa':5/6}
    assert preflight['configuration_frozen_before_fem'] is True
    assert config['geometry_selection']['backup_used'] is False
    assert preflight['first_axial_sorted_index']==config['geometry_selection']['first_axial_sorted_index']==8
    assert config['frozen_omega_window']==preflight['frequency_window']['omega_max']
    assert config['semantics']['linear_only'] is True
    assert config['semantics']['new_nonlinear_jobs'] is False
    for relative,digest in config['model_hashes'].items():
        assert workflow.sha(workflow.ROOT/relative)==digest


@pytest.mark.parametrize('level,through_h,through_b,nodes,elements',
                         [('coarse',2,4,2092,1054),('medium',3,6,5649,3120),('fine',4,8,11553,6670)])
def test_actual_three_meshes_have_quality_mass_and_real_edge_resolution(
        actual_fem_evidence,level,through_h,through_b,nodes,elements):
    _,summary=actual_fem_evidence
    case=summary['meshes'][level]
    audit=case['mesh_audit']
    assert case['status']==audit['status']=='PASS'
    assert audit['nodes']==nodes
    assert audit['c3d10_elements']==elements
    assert audit['solid_element_types']==['C3D10']
    assert audit['bbox_matches'] is True
    assert audit['volume']==pytest.approx(.02,rel=1e-12)
    assert audit['mass']==pytest.approx(.02,rel=1e-12)
    assert audit['volume_face_connected_components']==1
    assert audit['nonconforming_quadratic_face_count']==0
    assert audit['negative_or_zero_jacobian_elements']==0
    assert audit['minimum_quadratic_jacobian']>0
    assert audit['rigid_body_kinematic_constraint_check']=='PASS'
    assert audit['fixed_left_count']>0 and audit['fixed_right_count']>0
    for direction,expected in [('thickness_h',through_h),('width_b',through_b)]:
        edges=audit['actual_resolution'][direction]['box_edge_resolution']['edges']
        assert len(edges)==4
        assert all(edge['linear_segment_count']==expected for edge in edges)
    assert len(case['jobs'])==3
    assert all(job['returncode']==0 and job['failure'] is None for job in case['jobs'])


@pytest.mark.parametrize('level',('coarse','medium','fine'))
def test_actual_modal_results_have_24_complete_vectors_and_units(actual_fem_evidence,level):
    bundle,summary=actual_fem_evidence
    case=summary['meshes'][level]
    modal=case['modal']
    assert modal['eigenpairs']==modal['parsed_eigenvectors']==24
    assert modal['maximum_clamp_relative_residual']==0.
    assert modal['maximum_printed_frequency_unit_error']<modal['unit_gate_relative']==1e-5
    assert modal['frequency_unit']=='rad per normalized time'
    assert modal['full_frequency_window_covered'] is True
    assert modal['legacy_FRD_reader_used'] is False
    rows=modal['raw_frequencies']
    assert len(rows)==24
    assert [int(row['raw_mode_number']) for row in rows]==list(range(1,25))
    assert np.all(np.diff([row['angular_frequency'] for row in rows])>0)
    for row in rows:
        assert row['eigenvalue']>0
        assert abs(row['angular_frequency']**2/row['eigenvalue']-1)<1e-5
        assert abs(2*math.pi*row['cyclic_frequency']/row['angular_frequency']-1)<1e-5
    with np.load(bundle/'meshes'/level/'modal_vectors.npz',allow_pickle=False) as vectors:
        ids=vectors['node_ids']
        assert len(ids)==case['mesh_audit']['nodes']
        assert len(np.unique(ids))==len(ids)
        assert vectors['nodes'].shape==(len(ids),3)
        fixed=np.isin(ids,case['mesh_audit']['fixed_left_ids']+case['mesh_audit']['fixed_right_ids'])
        for row in rows:
            displacement=vectors['U_'+str(row['raw_mode_number'])]
            assert displacement.shape==(len(ids),3)
            assert np.all(np.isfinite(displacement))
            assert np.max(np.abs(displacement))>0
            assert np.max(np.abs(displacement[fixed]))==0


def test_actual_shape_matches_cover_four_families_without_conflicts(actual_fem_evidence):
    _,summary=actual_fem_evidence
    assert summary['extensions_used']==0
    assert summary['job_calls']=={'gmsh':6,'ccx':3}
    for case in summary['meshes'].values():
        modal=case['modal']
        assert modal['all_eight_identified'] is True
        assert modal['additional_modes_in_window']==[]
        assert modal['matching_uses_frequencies'] is False
        matches=modal['matches']
        assert len(matches)==8
        assert len({row['fem_mode'] for row in matches})==8
        assert set(row['family'] for row in matches)=={'axial_mh','inplane_bending','outplane_bending','torsion'}
        for row in matches:
            assert row['status']=='MATCHED'
            assert row['fem_family']==row['family']
            assert row['conflict'] is False
            assert row['low_mac'] is False
            assert row['ambiguous'] is False
            assert row['mac']>=summary['config']['matching']['minimum_mac']
            assert row['margin']>=summary['config']['matching']['minimum_margin']


def test_actual_mesh_unresolved_modes_are_retained_as_partial(actual_fem_evidence):
    _,summary=actual_fem_evidence
    assert summary['config']['numerical_mesh_convergence_relative']==.001
    rows=summary['comparisons']
    unresolved={(row['family'],row['local_mode']) for row in rows if row['mesh_status']=='MESH_UNRESOLVED'}
    resolved={(row['family'],row['local_mode']) for row in rows if row['mesh_status']=='PASS'}
    assert unresolved=={('inplane_bending',1),('inplane_bending',2),('inplane_bending',3),('torsion',1),('torsion',2)}
    assert resolved=={('outplane_bending',1),('outplane_bending',2),('axial_mh',1)}
    for row in rows:
        recalculated=abs(row['omega_fine']-row['omega_medium'])/row['omega_fine']
        assert row['medium_fine_relative']==pytest.approx(recalculated,rel=1e-14)
        assert row['medium_fine_relative']<row['coarse_medium_relative']
        if (row['family'],row['local_mode']) in unresolved:
            assert row['medium_fine_relative']>.001
        else:
            assert row['medium_fine_relative']<=.001
    statuses=summary['statuses']
    assert statuses['NLSP_FEM1_3D_MESH_QUALITY']=='PASS'
    assert statuses['NLSP_FEM1_3D_MODAL_EXECUTION']=='PASS'
    assert statuses['NLSP_FEM1_MODE_IDENTIFICATION']=='PASS'
    assert statuses['NLSP_FEM1_3D_MESH_CONVERGENCE']=='PARTIAL'
    assert statuses['NLSP_FEM1_ALL_FAMILY_COMPARISON']=='PARTIAL'


@pytest.mark.parametrize('mode',('--report-only','--plot-only'))
def test_actual_cached_replay_is_readonly_and_has_zero_new_solver_calls(
        actual_fem_evidence,monkeypatch,mode):
    bundle,expected=actual_fem_evidence
    def forbidden(*args,**kwargs):
        raise AssertionError('Actual artifact replay attempted new numerical computation')
    for name in ('run_job','run_fem','build_preflight','write_preflight','identity'):
        monkeypatch.setattr(workflow,name,forbidden)
    monkeypatch.setattr(workflow.single,'generate_mesh_with_gmsh_cli',forbidden)
    monkeypatch.setattr(workflow.single,'run_calculix',forbidden)
    monkeypatch.setattr(workflow.mh,'finite_roots',forbidden)
    monkeypatch.setattr(np.linalg,'eig',forbidden)
    monkeypatch.setattr(np.linalg,'eigh',forbidden)
    # Exercise all actual-data plotting branches without rewriting scientific figures.
    from matplotlib.figure import Figure
    saved=[]
    monkeypatch.setattr(Figure,'savefig',lambda self,path,*args,**kwargs:saved.append(Path(path).name))
    assert (bundle/'figures').is_dir()
    assert workflow.main([mode,str(bundle)])==expected
    if mode=='--plot-only':
        assert sorted(saved)==sorted(name+'.'+extension for name in
            ('combined_linear_spectrum','frequency_difference_and_mesh_convergence','representative_all_family_shapes')
            for extension in ('pdf','png'))
    else:
        assert saved==[]
