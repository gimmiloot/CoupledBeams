"""Synthetic 1D FEM-2 tests; no 3D jobs, ODE or eigenfrequency roots."""
from pathlib import Path
from types import SimpleNamespace
import json
import numpy as np
import pytest
from scripts.analysis import verify_nlsp_nonlinear_static_3d_fem as sut
from scripts.lib import weakly_nonlinear_planar_dynamics as dyn
from scripts.lib import weakly_nonlinear_spatial_rod as rod


@pytest.fixture(scope='module')
def static_disc():
    root=Path(__file__).resolve().parents[1]
    # In final tests this resolves repository root. The model is DESERIALIZED,
    # not rederived, and no mass/dynamics inversion is used to define statics.
    if not (root/'results').exists(): root=Path.cwd()
    saved=json.loads((root/'results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295/result.json').read_text(encoding='utf8'))
    pol=saved['polynomials']
    model=SimpleNamespace(T4=rod.Polynomial.deserialize(pol['T4']),V4=rod.Polynomial.deserialize(pol['V4']),
        residual_a=tuple(rod.Polynomial.deserialize(p) for p in pol['residuals_A']),
        symbols={n:rod.Polynomial.symbol(n) for n in rod.SYMBOL_ORDER})
    coef=json.loads((root/'results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad/preflight.json').read_text(encoding='utf8'))['coefficients']
    return dyn.PlanarGalerkin(rod.RodCoefficients(**coef),p=10,nq=21,model=model)


def test_static_line_load_physical_mapping(static_disc):
    q=3e-5
    force=sut.fem2_line_load(static_disc,q)
    assert np.count_nonzero(force[:static_disc.slices['w'].start])==0
    assert np.count_nonzero(force[static_disc.slices['w'].stop:])==0
    vector=np.zeros(static_disc.ndof)
    # An arbitrary admissible w has the expected external work.
    vector[static_disc.slices['w']]=np.arange(static_disc.n)*1e-7
    assert force@vector==pytest.approx(q*static_disc.weights@static_disc.reconstruct(vector)[:,1],rel=2e-15)


def test_static_uniform_tim_is_independent_of_mh(static_disc):
    x=np.linspace(0,1,101);p=static_disc.coefficients;q=3e-5
    fields,first=sut.fem2_tim_uniform(x,q,1,p.Bp,p.S)
    assert np.max(abs(fields[:,[0,3]]))==0
    assert np.max(abs(fields[[0,-1]]))<1e-16
    assert fields[50,1]==pytest.approx(q*(1/(384*p.Bp)+1/(8*p.S)),rel=3e-15)
    assert first[0,1]==pytest.approx(q/(2*p.S),rel=3e-15)
    assert first[-1,1]==pytest.approx(-q/(2*p.S),rel=3e-15)
    assert first[0,1]!=0 # No spurious slope clamp.


def test_static_linear_coefficients_match_exact_tim(static_disc):
    q=3e-5
    a=np.linalg.solve(static_disc.linear_stiffness,sut.fem2_line_load(static_disc,q))
    x=np.linspace(0,1,101)
    expected,_=sut.fem2_tim_uniform(x,q,1,static_disc.coefficients.Bp,static_disc.coefficients.S)
    np.testing.assert_allclose(static_disc.reconstruct(a,x),expected,rtol=2e-12,atol=2e-16)


def test_static_quartic_gradient_directional_identity(static_disc):
    rng=np.random.default_rng(20261009)
    q=rng.normal(size=static_disc.ndof)*1e-7
    d=rng.normal(size=static_disc.ndof);d/=np.linalg.norm(d);h=1e-8
    base=static_disc.potential(q)
    fd=(static_disc.potential(q+h*d,gradient=False)['V']-static_disc.potential(q-h*d,gradient=False)['V'])/(2*h)
    assert fd==pytest.approx(float(base['gradient']@d),rel=1e-7,abs=1e-12)


def test_static_analytic_hessian_consistent(static_disc):
    rng=np.random.default_rng(10091)
    q=rng.normal(size=static_disc.ndof)*1e-7
    d=rng.normal(size=static_disc.ndof);d/=np.linalg.norm(d);h=1e-8
    value=static_disc.potential(q,hessian=True)
    fd=(static_disc.potential(q+h*d)['gradient']-static_disc.potential(q-h*d)['gradient'])/(2*h)
    np.testing.assert_allclose(fd,value['hessian']@d,rtol=2e-8,atol=1e-9)
    symmetry=np.linalg.norm(value['hessian']-value['hessian'].T)/np.linalg.norm(value['hessian'])
    assert symmetry<=64*np.finfo(float).eps


def test_static_unloaded_newton_without_rhs(static_disc,monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('No dynamic inversion or ODE in statics')
    monkeypatch.setattr(static_disc,'rhs',forbidden)
    monkeypatch.setattr(static_disc,'acceleration',forbidden)
    monkeypatch.setattr(static_disc,'_solve_mass',forbidden)
    result=sut.fem2_static_newton(static_disc,np.zeros(static_disc.ndof),{'load_steps':2})
    assert result['status']=='PASS'
    assert result['load_factor']==1
    assert np.max(abs(result['coordinate']))==0
    assert all(r['relative_residual']==0 for r in result['history'])


def test_static_flux_reactions_analytic_sign_and_balance(static_disc):
    p=static_disc.coefficients;q=3e-5
    a=np.linalg.solve(static_disc.linear_stiffness,sut.fem2_line_load(static_disc,q))
    result,profile=sut.fem2_static_diagnostics(static_disc,a,q,{'h':.1},linear=True)
    r=np.asarray(result['support_reactions_local_u_w_M_Rc'])
    np.testing.assert_allclose(r[:,1],[-q/2,-q/2],rtol=2e-12,atol=1e-17)
    np.testing.assert_allclose(r[:,2],[-q/12,+q/12],rtol=2e-12,atol=1e-17)
    assert result['force_balance_relative']<1e-10
    assert result['moment_balance_relative']<1e-10
    assert result['endpoint_BC_absolute_max']==0
    assert result['slope_constraints'] is False
    assert 'endpoint partial' in result['reaction_source']


def test_static_load_selected_before_fem_and_keeps_primary(static_disc):
    c={'geometry':{'L':1.,'b':.2,'h':.1},'material':{'E':1.,'rho':1.,'nu':.3,'kappa':5/6},'load_policy':{'primary_w_over_h':.05,'backup_w_over_h':.03,'bending_surface_strain_ceiling':.01}}
    r=sut.fem2_load_selection(c,static_disc.coefficients)
    assert r['status']=='PASS'
    assert r['g']==pytest.approx(.0014224751066856333,rel=2e-15)
    assert r['q']==pytest.approx(c['material']['rho']*.2*.1*r['g'],rel=2e-15)
    assert r['F_total']==r['q']
    assert len(r['candidate_history'])==1
    assert r['target_w_over_h']==.05
    assert r['follower_load'] is False
    assert r['body_direction_global']==[0.,-1.,0.]


def test_static_load_only_predeclared_backup(static_disc):
    c={'geometry':{'L':1.,'b':.2,'h':.1},'material':{'rho':1.},'load_policy':{'primary_w_over_h':.05,'backup_w_over_h':.03,'bending_surface_strain_ceiling':.005}}
    r=sut.fem2_load_selection(c,static_disc.coefficients)
    assert r['target_w_over_h']==.03
    assert len(r['candidate_history'])==2
    c['load_policy']['bending_surface_strain_ceiling']=.001
    assert sut.fem2_load_selection(c,static_disc.coefficients)['status']=='FAIL'


def test_static_all_four_independent_coordinates(static_disc):
    assert static_disc.ndof==4*(static_disc.p-1)
    assert tuple(static_disc.slices)==('u','w','theta','c')
    x=np.array([0.,1.]);r=np.random.default_rng(99).normal(size=static_disc.ndof)
    assert np.max(abs(static_disc.reconstruct(r,x)))==0
    assert np.max(abs(static_disc.reconstruct(r,x,1)))>0


def test_static_profile_difference_no_alignment():
    x=np.linspace(0,1,101)
    profiles=np.column_stack((x*(1-x),2*x*(1-x),3*x*(1-x),4*x*(1-x)))
    rows=sut.fem2_static_comparison(profiles,-profiles,x,1.)
    assert all(row['relative_max']==pytest.approx(2.) for row in rows)
    assert all(row['relative_L2']==pytest.approx(2.) for row in rows)
    assert [row['field'] for row in rows]==list(sut.FEM2_FIELDS)


def test_static_strict_work_scale_not_redefined(static_disc):
    p=static_disc.coefficients;q=3e-5
    a=np.linalg.solve(static_disc.linear_stiffness,sut.fem2_line_load(static_disc,q))
    result,_=sut.fem2_static_diagnostics(static_disc,a,q,{'h':.1},linear=True)
    assert result['strong_action_strict_threshold']==2e-12
    assert result['strong_action_uncancelled_work_scale']>0
    assert result['dynamics_rhs_calls']==0
    assert result['dynamics_jacobian_calls']==0
    assert result['linear_eigendecompositions']==0

"""Synthetic static I/O tests; no executable FEM jobs."""
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest


def _io():
    from scripts.analysis import verify_nlsp_nonlinear_static_3d_fem as module
    return module


def _io_tet_mesh():
    corners=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    edges=((0,1),(1,2),(0,2),(0,3),(1,3),(2,3))
    xyz=np.concatenate((corners,np.array([(corners[i]+corners[j])/2 for i,j in edges])))
    return SimpleNamespace(nodes={i+1:tuple(x) for i,x in enumerate(xyz)},solid_elements={1:tuple(range(1,11))})


def _frd_header(counter,inc,step,time,nodes):
    line=list(' '*75)
    line[0:7]='  100CL'
    line[7:12]=f'{100+inc:5d}'
    line[12:24]=f'{time:12.9f}'
    line[24:36]=f'{nodes:12d}'
    return f'    1PSTEP{counter:26d}{inc:12d}{step:12d}\n'+''.join(line)+'\n'


def _frd_block(counter,inc,time,name,values):
    text=_frd_header(counter,inc,1,time,len(values))
    text+=f' -4  {name:<12} {values.shape[1]:d}    1\n'
    for i,row in enumerate(values,1):
        text+=f' -1{i:10d}'+''.join(f'{float(v):12.5E}' for v in row)+'\n'
    return text+' -3\n'


def test_static_deck_pair_only_nlgeom_difference(tmp_path):
    m=_io();mesh=_io_tet_mesh();audit={'status':'PASS','fixed_left_ids':[1,3,4],'fixed_right_ids':[2]}
    args=(tmp_path/'mesh.inp',mesh,audit,{'E':1.,'rho':1.,'nu':.3},.0002)
    m.write_static_input(tmp_path/'linear.inp',*args,False)
    m.write_static_input(tmp_path/'nonlinear.inp',*args,True)
    lin=(tmp_path/'linear.inp').read_text();nl=(tmp_path/'nonlinear.inp').read_text()
    assert nl.replace(', NLGEOM','')==lin
    assert 'SOLID,GRAV,2.000000000000E-04,0,-1,0' in lin
    assert '*STATIC' in lin and '*DYNAMIC' not in lin and '*FREQUENCY' not in lin
    assert 'LEFT_FIXED,1,3,0' in lin and 'RIGHT_FIXED,1,3,0' in lin
    assert '*MPC' not in lin and '*CLOAD' not in lin
    assert 'ALL_NODES' in lin and 'S,E' in lin
    assert 'PARAMETERS=FIELD' in lin


def test_static_deck_requires_source_mesh_gate(tmp_path):
    with pytest.raises(ValueError,match='mesh'):
        _io().write_static_input(tmp_path/'bad.inp',tmp_path/'mesh.inp',_io_tet_mesh(),
            {'status':'FAIL'}, {'E':1.,'rho':1.,'nu':.3},.0002,True)


def test_consistent_bodyload_exact_simplex_moments():
    ids,xyz,load,volume=_io().consistent_gravity_loads(_io_tet_mesh(),2.,3.)
    assert volume==pytest.approx(1/6)
    np.testing.assert_allclose(load[:4,1],np.full(4,.05))
    np.testing.assert_allclose(load[4:,1],np.full(6,-.2))
    np.testing.assert_allclose(np.sum(load,axis=0),[0.,-1.,0.],atol=2e-16)
    np.testing.assert_allclose(np.sum(np.cross(xyz,load),axis=0),[.25,0.,-.25],atol=2e-16)


def test_bodyload_requires_affine_geometry():
    mesh=_io_tet_mesh();mesh.nodes[5]=(0.51,0.,0.)
    with pytest.raises(ValueError,match='Nonaffine'):_io().consistent_gravity_loads(mesh,1.,1.)


def test_static_frd_final_time_increment_all_fields(tmp_path):
    values=np.array([[.002,-.003,.004],[.005,.006,-.007]])
    text=_frd_block(1,1,.1,'DISP',values*0.)
    for counter,(name,value) in enumerate((('DISP',values),('FORC',values),('STRESS',np.c_[values,values]),('TOSTRAIN',np.c_[values,values])),2):
        text+=_frd_block(counter,10,1.,name,value)
    path=tmp_path/'final.frd';path.write_text(text)
    fields,meta=_io().read_static_frd(path,[1,2])
    np.testing.assert_allclose(fields['DISP'],values)
    assert meta['step']==1 and meta['increment']==10 and meta['time']==1.
    assert fields['STRESS'].shape==(2,6)


def test_static_frd_keeps_positive_payload_separate_from_node():
    node,value=_io()._static_frd_record(f' -1{9999:10d}'+'1.12345E-003-1.23456E-0021.65432E+000',3)
    assert node==9999
    np.testing.assert_allclose(value,[.00112345,-.0123456,1.65432])


def test_static_frd_rejects_modal(tmp_path):
    path=tmp_path/'modal.frd';path.write_text('    1PMODE           1\n')
    with pytest.raises(ValueError,match='Modal'):_io().read_static_frd(path,[1])


def test_static_frd_rejects_missing_final_field(tmp_path):
    path=tmp_path/'bad.frd';path.write_text(_frd_block(1,1,1.,'DISP',np.ones((2,3))))
    with pytest.raises(ValueError,match='Missing final'):_io().read_static_frd(path,[1,2])


def test_static_dat_final_u_and_support_rf(tmp_path):
    text=' displacements (vx,vy,vz) for set ALL_NODES and time 0.1000000E+00\n1 0.0 0.0 0.0\n2 0.0 0.0 0.0\n'
    text+=' displacements (vx,vy,vz) for set ALL_NODES and time 0.1000000E+01\n1 0.01 -0.02 0.03\n2 0.04 0.05 0.06\n'
    text+=' forces (fx,fy,fz) for set LEFT_FIXED and time 0.1000000E+01\n1 0.1 0.2 0.3\n total force for set LEFT_FIXED\n 0.1 0.2 0.3\n'
    text+=' forces (fx,fy,fz) for set RIGHT_FIXED and time 0.1000000E+01\n2 -0.1 0.2 -0.3\n'
    path=tmp_path/'static.dat';path.write_text(text)
    data,meta=_io().read_static_dat(path,[1,2],[1],[2])
    assert meta['time']==1.
    np.testing.assert_allclose(data['ALL_NODES'][0],[.01,-.02,.03])
    np.testing.assert_allclose(data['LEFT_FIXED'][0],[.1,.2,.3])


def test_static_dat_rejects_missing_node(tmp_path):
    path=tmp_path/'static.dat';path.write_text(' displacements (vx,vy,vz) for set ALL_NODES and time 1.0\n1 0.0 0.0 0.0\n')
    with pytest.raises(ValueError,match='Incomplete'):_io().read_static_dat(path,[1,2],[1],[2])


def test_static_sta_actual_increment_times(tmp_path):
    path=tmp_path/'static.sta';path.write_text(' STEP INC ATT ITRS TOT TIME STEP TIME INC TIME\n1 1 1 3 .1 .1 .1\n1 2 1 4 .2 .2 .1\n1 10 1 3 1.0 1.0 .1\n')
    data=_io().read_static_sta(path)
    assert data['last_time']==1.
    assert len(data['accepted_increments'])==3
    assert data['accepted_increments'][-1]['increment']==10


def test_static_reaction_contract_not_blind_rf_sum():
    facts=_io().static_documentation_evidence()
    assert 'RF_support - integral_reference' in facts['reaction_recovery']
    assert facts['RF']['pages']==[558,562]
    assert 'StVK' in facts['constitutive_comparison']


"""Synthetic no-job tests for the scoped FEM-2 recovery fragment."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from scripts.analysis import verify_nlsp_nonlinear_static_3d_fem as rec


def _samples(count=41):
    rows = []
    for k in range(count):
        for d in (-.4, -.15, .15, .4):
            for y in (-.04, .04):
                for z in (-.08, .08):
                    rows.append(((k+.5+d)/count, y, z))
    local = np.array(rows)
    return local, local*np.array((1., -1., -1.)), np.ones(len(rows))*.02/len(rows)


def _recover_affine(F, T=None):
    local, xyz, wt = _samples()
    translate = np.zeros(3) if T is None else np.asarray(T)
    disp = (local@(F-np.eye(3)).T+translate)*np.array((1., -1., -1.))
    return rec.fem2_recover_reference_samples(xyz, disp, wt, 1., .1, .2,
                                             enforce_clamped_faces=False)


def _rec_tet_mesh():
    corners = np.array(((0., 0., 0.), (1., 0., 0.), (0., 1., 0.), (0., 0., 1.)))
    mids = np.array([(corners[i]+corners[j])/2 for i, j in rec.fem1.TET10_EDGES])
    xyz = np.vstack((corners, mids))
    return SimpleNamespace(nodes={i+1: row for i, row in enumerate(xyz)},
                           solid_elements={1: tuple(range(1, 11))}), xyz


def test_zero_and_reference_geometry():
    out = _recover_affine(np.eye(3))
    assert np.max(np.abs(out['fields'])) < 1e-12
    assert out['reference_geometry_used_for_sections']
    assert not out['additional_derivative_constraints']
    assert out['recovery_policy'] == rec.FEM2_RECOVERY_POLICY


def test_translation_axis_displacement_and_sign():
    T = np.array((.001, .005, -.002))
    out = _recover_affine(np.eye(3), T)
    np.testing.assert_allclose(out['fields'][:, :3], np.tile(T, (41, 1)), atol=1e-12)


@pytest.mark.parametrize('angle', [.001, .18, -.27])
def test_finite_inplane_rotation_no_small_angle_contamination(angle):
    R = Rotation.from_rotvec((0., 0., angle)).as_matrix()
    out = _recover_affine(R)
    np.testing.assert_allclose(out['fields'][:, 5], angle, atol=1e-12)
    np.testing.assert_allclose(out['fields'][:, 6], 0., atol=1e-12)
    np.testing.assert_allclose(out['small_rotation_theta'], np.sin(angle), atol=1e-12)
    np.testing.assert_allclose(out['fields'][:, 0], (np.cos(angle)-1)*out['x'], atol=1e-12)
    np.testing.assert_allclose(out['fields'][:, 1], np.sin(angle)*out['x'], atol=1e-12)


def test_finite_torsion_and_rotation_signs():
    R = Rotation.from_rotvec((.12, 0., 0.)).as_matrix()
    out = _recover_affine(R)
    np.testing.assert_allclose(out['fields'][:, 3], .12, atol=1e-12)
    np.testing.assert_allclose(out['fields'][:, 4:6], 0., atol=1e-12)


def test_outplane_psi_sign():
    R = Rotation.from_rotvec((0., -.12, 0.)).as_matrix()
    out = _recover_affine(R)
    np.testing.assert_allclose(out['fields'][:, 4], .12, atol=1e-12)


def test_affine_stretch_separate_from_rotation():
    R = Rotation.from_rotvec((0., 0., .12)).as_matrix()
    out = _recover_affine(R@np.diag((1.001, .98, 1.02)))
    np.testing.assert_allclose(out['fields'][:, 5], .12, atol=1e-12)
    np.testing.assert_allclose(out['fields'][:, 6], -.02, atol=1e-12)
    assert max(abs(r['effective_width_strain']-.02) for r in out['section_rows']) < 1e-12
    assert out['contraction_status'].endswith('NOT_MH_DOF')


def test_linear_rotation_physical_orientation_policy_explicit():
    F = np.eye(3)
    F[0, 1] = -.2
    out = _recover_affine(F)
    np.testing.assert_allclose(out['fields'][:, 5], np.arctan(.2), atol=1e-12)
    np.testing.assert_allclose(out['small_rotation_theta'], .2, atol=1e-12)
    assert 'both_linear_and_nonlinear' in out['theta_policy']


def test_xcubic_axis_variation_not_mislabeled_rotation():
    local, xyz, wt = _samples()
    x = local[:, 0]
    displacement = np.column_stack((.001*x*x, .005*x**3, np.zeros(len(x))))
    out = rec.fem2_recover_reference_samples(xyz, displacement*np.array((1., -1., -1.)), wt,
                                            1., .1, .2, enforce_clamped_faces=False)
    np.testing.assert_allclose(out['fields'][:, 0], .001*out['x']**2, atol=1e-12)
    np.testing.assert_allclose(out['fields'][:, 1], .005*out['x']**3, atol=1e-12)
    np.testing.assert_allclose(out['fields'][:, 3:], 0., atol=1e-12)


def test_polar_frame_right_handed_and_orthogonal():
    out = rec.fem2_polar_section_orientation(np.array(((.1, 0.), (.99, .1), (.03, 1.))))
    np.testing.assert_allclose(out['matrix'].T@out['matrix'], np.eye(3), atol=1e-12)
    assert np.linalg.det(out['matrix']) > 0.


def test_polar_degenerate_section_rejected():
    with pytest.raises(ValueError, match='Degenerate'):
        rec.fem2_polar_section_orientation(np.zeros((3, 2)))


def test_reference_positive_weight_required():
    local, xyz, wt = _samples()
    wt[0] = 0.
    with pytest.raises(ValueError, match='positive'):
        rec.fem2_recover_reference_samples(xyz, xyz*0, wt, 1., .1, .2)


def test_independent_FE_gradient_affine_strain_measures():
    mesh, xyz = _rec_tet_mesh()
    gradient = np.array(((.001, .01, 0.), (.002, -.003, 0.), (0., 0., .004)))
    out = rec.fem2_fe_strain_diagnostics(mesh, xyz@gradient.T)
    linear = .5*(gradient+gradient.T)
    green = linear+.5*gradient.T@gradient
    assert out['linear_strain_max_abs'] == pytest.approx(np.max(np.abs(linear)), abs=1e-14)
    assert out['green_lagrange_max_abs'] == pytest.approx(np.max(np.abs(green)), abs=1e-14)
    assert not out['CCX_stress_strain_measures_replaced']


def test_rigid_finite_rotation_zero_green_lagrange():
    mesh, xyz = _rec_tet_mesh()
    R = Rotation.from_rotvec((.08, -.10, .13)).as_matrix()
    out = rec.fem2_fe_strain_diagnostics(mesh, xyz@(R-np.eye(3)).T)
    assert out['green_lagrange_max_abs'] < 1e-14
    assert out['linear_strain_max_abs'] > 1e-3
    assert out['minimum_det_deformation_gradient'] == pytest.approx(1., abs=1e-13)


def test_difference_metrics_no_alignment_and_fixed_scale():
    x = np.linspace(0, 1, 801)
    out = rec.fem2_curve_difference(x, x, .99*x, fixed_scale=.1, characteristic_scale=2.)
    assert out['absolute_max'] == pytest.approx(.01)
    assert out['absolute_L2'] == pytest.approx(.01/np.sqrt(3), rel=1e-6)
    assert out['relative_max'] == pytest.approx(.005)
    assert out['fixed_scale_max'] == pytest.approx(.1)
    assert not out['alignment_used']


def test_zero_characteristic_scale_not_divided():
    x = np.linspace(0, 1, 9)
    out = rec.fem2_curve_difference(x, x*0, x*0, characteristic_scale=0.)
    assert out['relative_max'] is None


def test_signal_resolution_keeps_measured_uncertainty_visible():
    assert rec.fem2_signal_resolution(1e-5, 2e-5)['status'] == 'UNRESOLVED_SIGNAL'
    out = rec.fem2_signal_resolution(1e-5, 2e-6, 4e-6)
    assert out['status'] == 'SIGNAL_EXCEEDS_OBSERVED_UNCERTAINTY'
    assert out['signal_to_largest_observed_uncertainty'] == pytest.approx(2.5)
    assert not out['continuum_error_bound_claimed']


def test_profile_policy_same_comparison_and_corrections():
    x = np.linspace(0, 1, 801)
    lin = np.zeros((801, 7))
    lin[:, 1] = .005*x*x*(1-x)**2
    nl = lin.copy()
    nl[:, 1] -= .00001*x*x*(1-x)**2
    p, q = {'x': x, 'fields': lin}, {'x': x, 'fields': nl}
    out = rec.compare_static_profiles(p, q, p, q, correction_scales={k: 1. for k in rec.FEM2_STATIC_FIELD_ORDER})
    assert out['metrics']['w']['nonlinear_correction_1D_vs_3D']['absolute_max'] == 0.
    assert out['metrics']['w']['correction_signed_midpoint_3D'] < 0.
    assert not out['amplitude_or_space_alignment']
    assert not out['frequency_mesh_criterion_applied']


# Orchestration and actual stopped-at-input evidence; no real job runs.
from scripts.analysis import verify_nlsp_nonlinear_static_3d_fem as workflow

def test_ccx_numeric_fields_fit_twenty_character_parser(tmp_path):
    mesh=_io_tet_mesh();audit={'status':'PASS','fixed_left_ids':[1,3,4],'fixed_right_ids':[2]}
    p=tmp_path/'fixed.inp'
    workflow.write_static_input(p,tmp_path/'mesh.inp',mesh,audit,{'E':1.,'rho':1.,'nu':.3},.0014224751066856333,False)
    lines=p.read_text().splitlines()
    values=lines[lines.index('*STATIC')+1].split(',')
    assert all(len(v.strip())<=20 for v in values)
    assert [float(v) for v in values]==pytest.approx([.1,1.,1e-6,.1],rel=1e-12)
    grav=lines[lines.index('*DLOAD')+1].split(',')
    assert len(grav[2])<=20
    assert float(grav[2])==pytest.approx(.0014224751066856333,rel=1e-12)
    # The rejected original .17g exponent is demonstrably truncated mid-exponent.
    old=f'{1e-6:.17g}'
    assert len(old)>20
    with pytest.raises(ValueError):float(old[:20])


@pytest.fixture(scope='module')
def stopped_static_evidence():
    b=workflow.FEM2_OUTPUT/'6714bd9f2778e6d7'
    if not (b/'manifest.json').exists():pytest.skip('Stopped FEM2 evidence absent; never recreate')
    return b,workflow.validate_fem2_cache(b)


def test_failed_medium_is_not_a_physical_solution(stopped_static_evidence):
    b,s=stopped_static_evidence
    case=s['cases']['medium']['linear']
    assert case['status']=='FAIL'
    assert case['job']['returncode']==201
    assert s['job_calls']=={'ccx':1,'gmsh':0,'modal':0,'nonlinear_ODE':0}
    assert s['statuses']['NLSP_FEM2_3D_LINEAR_STATIC']=='FAIL'
    assert s['statuses']['NLSP_FEM2_3D_NONLINEAR_STATIC']=='NOT_RUN'
    log=(b/'cases/medium/linear/static.stdout.txt').read_text()
    assert '*ERROR reading *STATIC' in log
    assert 'JOB FINISHED' not in log
    assert not (b/'cases/medium/nonlinear').exists()
    assert not (b/'cases/fine').exists() and not (b/'cases/refined').exists()


@pytest.mark.parametrize('mode',['--preflight','--run-fem'])
def test_code_fix_cannot_silently_repeat_failed_attempt(stopped_static_evidence,monkeypatch,mode):
    _,expected=stopped_static_evidence
    def blocked(*args,**kwargs):raise AssertionError('Unexpected new computation')
    monkeypatch.setattr(workflow,'build_static_preflight',blocked)
    monkeypatch.setattr(workflow,'run_fem2_cases',blocked)
    monkeypatch.setattr(workflow,'fem2_identity',blocked)
    monkeypatch.setattr(workflow.fem1,'run_job',blocked)
    assert workflow.main([mode])==expected


@pytest.mark.parametrize('mode',['--report-only','--plot-only'])
def test_stopped_report_and_plot_have_zero_solver_calls(stopped_static_evidence,monkeypatch,mode):
    b,expected=stopped_static_evidence
    def blocked(*args,**kwargs):raise AssertionError('Unexpected numerical call')
    monkeypatch.setattr(workflow,'build_static_preflight',blocked)
    monkeypatch.setattr(workflow,'run_fem2_cases',blocked)
    monkeypatch.setattr(workflow.fem1,'run_job',blocked)
    monkeypatch.setattr(workflow.fem1.single,'generate_mesh_with_gmsh_cli',blocked)
    monkeypatch.setattr(workflow.fem2_rod,'derive_polynomials',blocked)
    # Test plotting branches without rewriting the artifact figures.
    from matplotlib.figure import Figure
    saved=[];monkeypatch.setattr(Figure,'savefig',lambda self,path,*a,**kw:saved.append(str(path)))
    assert workflow.main([mode,str(b)])==expected
    assert len(saved)==(4 if mode=='--plot-only' else 0)


def test_early_failure_is_fail_and_unattempted_stages_not_run(tmp_path):
    c=workflow.read_json(workflow.FEM2_CONFIG)
    s={'cases':{'medium':{'linear':{'status':'FAIL'}}},
       'statuses':{'NLSP_FEM2_'+n:'NOT_RUN' for n in workflow.FEM2_STATUS_NAMES}}
    result=workflow.finalize_fem2_summary(c,tmp_path,s)
    assert result['statuses']['NLSP_FEM2_3D_LINEAR_STATIC']=='FAIL'
    assert result['statuses']['NLSP_FEM2_3D_NONLINEAR_STATIC']=='NOT_RUN'
    assert result['statuses']['NLSP_FEM2_STATIC_EQUILIBRIUM']=='NOT_RUN'
    assert result['overall']=='PARTIAL'


def test_sources_and_frozen_physics_are_verified_read_only(stopped_static_evidence):
    b,s=stopped_static_evidence;c=workflow.read_json(workflow.FEM2_CONFIG)
    before={k:workflow.sha(workflow.ROOT/v['bundle']/'manifest.json') for k,v in c['sources'].items()}
    sources=workflow.load_fem2_sources(c)
    assert set(sources)=={'fem1','fem1r'}
    after={k:workflow.sha(workflow.ROOT/v['bundle']/'manifest.json') for k,v in c['sources'].items()}
    assert before==after=={k:v['manifest_sha256'] for k,v in c['sources'].items()}
    m=workflow.read_json(b/'manifest.json')
    for p,d in m['identity']['helper_sha256'].items():
        assert workflow.sha(workflow.ROOT/p)==d


def test_actual_static_preflight_is_same_load_and_all_four_independent_fields(stopped_static_evidence):
    _,s=stopped_static_evidence;pre=s['preflight']
    assert pre['status']=='PASS'
    assert pre['load']['target_w_over_h']==.05
    assert pre['load']['linear_bending_surface_strain']<.01
    assert pre['load']['q']==pytest.approx(pre['load']['g']*.2*.1,rel=4*np.finfo(float).eps,abs=0)
    assert pre['load']['F_total']==pre['load']['q']
    assert [r['p'] for r in pre['cases']]==[48,64]
    for row in pre['cases']:
        assert row['linear_status']==row['nonlinear_status']=='PASS'
        assert row['fields']==['u','w','theta','c']
        assert row['linear']['endpoint_BC_absolute_max']==row['nonlinear']['endpoint_BC_absolute_max']==0
        assert row['nonlinear']['tangent_positive']
        assert row['nonlinear']['relative_residual']<1e-10
        assert row['delta_w_midspan']<0
    assert pre['spatial_comparison']['delta_w']['relative_max']<1e-10


def test_no_production_model_change_is_hidden_in_static_fix(stopped_static_evidence):
    b,s=stopped_static_evidence
    m=workflow.read_json(b/'manifest.json')
    execution=b/'execution_code'/Path(workflow.__file__).name
    assert workflow.sha(execution)==m['identity']['code_sha256']
    assert workflow.sha(execution)!=workflow.sha(workflow.__file__)
    assert s['cases']['medium']['linear']['source_mesh_include_sha256']==workflow.sha(
        workflow.ROOT/s['cases']['medium']['linear']['source_mesh']/'solid_mesh.inp')

