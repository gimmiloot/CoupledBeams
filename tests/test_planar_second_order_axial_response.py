"""Leading axial-response controls with zero numerical time integrations.

All finite-dimensional MH coordinates are retained. The small p8/p16
validation eigendecompositions are separate from the main one per p.
"""
from __future__ import annotations
import ast
import json
from pathlib import Path
import numpy as np
import pytest
from scipy.linalg import expm
from scripts.lib import planar_second_order_axial_response as leading
from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib import weakly_nonlinear_planar_dynamics as planar
ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / 'data/input/weakly_nonlinear_planar_time_pilot.json'

@pytest.fixture(scope='module')
def source_config():
    return json.loads(CONFIG.read_text(encoding='utf8'))

@pytest.fixture(scope='module')
def symbolic_model():
    return rod.derive_polynomials()

@pytest.fixture(scope='module')
def coefficients(source_config):
    path = ROOT / source_config['audit_bundle'] / 'result.json'
    return rod.RodCoefficients(**json.loads(path.read_text(encoding='utf8'))['coefficients'])

@pytest.fixture(scope='module')
def background(source_config):
    return leading.background_from_pilot(source_config, root=ROOT)

@pytest.fixture(scope='module', params=(8,16))
def forced(request, coefficients, background, symbolic_model):
    return leading.SecondOrderAxial(coefficients,request.param,background,model=symbolic_model)

@pytest.mark.parametrize('omega,drive',((0.,0.),(1.7,0.),(1.7,1.7),(1.7,1.7+1e-13),(2.3,.8)))
def test_exact_time_kernel_satisfies_initial_data_and_scalar_ode(omega,drive):
    times=np.array([0.,1e-10,.13,.9,3.])
    value=leading.response_kernel(omega,drive,times,derivative=0)
    speed=leading.response_kernel(omega,drive,times,derivative=1)
    acceleration=leading.response_kernel(omega,drive,times,derivative=2)
    assert value[0]==speed[0]==0.
    assert acceleration[0]==1.
    np.testing.assert_allclose(acceleration+omega**2*value,np.cos(drive*times),rtol=2e-13,atol=2e-13)
    assert np.all(np.isfinite(value)) and np.all(np.isfinite(speed))
    np.testing.assert_allclose(value[1]/times[1]**2,.5,rtol=2e-13)
    np.testing.assert_allclose(speed[1]/times[1],1.,rtol=2e-13)

def test_exact_resonance_kernel_and_derivatives():
    omega=1.9;times=np.array([0.,1e-8,.2,2.,17.])
    expected=times*np.sin(omega*times)/(2*omega)
    speed=(np.sin(omega*times)+omega*times*np.cos(omega*times))/(2*omega)
    acceleration=np.cos(omega*times)-omega*times*np.sin(omega*times)/2
    for order,reference in enumerate((expected,speed,acceleration)):
        np.testing.assert_allclose(leading.response_kernel(omega,omega,times,derivative=order),reference,rtol=3e-13,atol=3e-13)

def test_kernel_near_resonance_and_small_time_against_high_precision():
    import mpmath as mp
    with mp.workdps(70):
        omega,drive=2.,2.+1e-12;times=np.array([0.,1e-12,.5,4.,20.]);w,d=mp.mpf(omega),mp.mpf(drive)
        references=[[],[],[]]
        for t_float in times:
            t=mp.mpf(t_float);denominator=w*w-d*d
            references[0].append(float((mp.cos(d*t)-mp.cos(w*t))/denominator))
            references[1].append(float((w*mp.sin(w*t)-d*mp.sin(d*t))/denominator))
            references[2].append(float((w*w*mp.cos(w*t)-d*d*mp.cos(d*t))/denominator))
    for order,reference in enumerate(references):
        actual=leading.response_kernel(omega,drive,times,derivative=order)
        np.testing.assert_allclose(actual,reference,rtol=3e-13,atol=1e-25 if order==0 else 3e-13)

def test_constant_forcing_kernel_does_not_subtract_nearly_equal_cosines():
    omega=3.;times=np.array([0.,1e-14,1e-9,.3,2.]);expected=2*np.sin(omega*times/2)**2/omega**2
    np.testing.assert_allclose(leading.response_kernel(omega,0.,times),expected,rtol=3e-13,atol=0)
    np.testing.assert_allclose(leading.response_kernel(omega,0.,times,derivative=1),np.sin(omega*times)/omega,rtol=3e-13,atol=0)

def test_no_test_or_helper_executes_an_ode_integrator():
    for path in (Path(__file__),Path(leading.__file__)):
        tree=ast.parse(path.read_text(encoding='utf8'))
        for node in ast.walk(tree):
            if isinstance(node,ast.Call):
                name=node.func.id if isinstance(node.func,ast.Name) else(node.func.attr if isinstance(node.func,ast.Attribute) else None)
                assert name not in {'solve_ivp','Radau','integrate_case'}

def test_exact_fraction_extraction_and_bending_second_order_zero(symbolic_model):
    extracted=leading.derive_second_order(symbolic_model);s=symbolic_model.symbols
    assert all(extracted['checks'].values())
    assert extracted['linear_u']==s['m']*s['u_tt']-s['C']*s['u_ss']-s['nu']*s['C']*s['c_s']
    assert extracted['linear_c']==s['jp']*s['c_tt']-s['H']*s['c_ss']+s['C']*(s['c']+s['nu']*s['u_s'])
    expectedZ=(s['C']-s['S'])*s['theta']*s['w_s']+(s['S']-s['C']/2)*s['theta']**2
    assert extracted['forcing_u']==expectedZ.total_derivative('s')
    assert extracted['forcing_c']==-s['nu']*s['C']*(s['theta']*s['w_s']-s['theta']**2/2)+s['jp']*s['theta_t']**2
    assert all(p==rod.Polynomial()for p in extracted['bending_second_order'])
    for polynomial in(extracted['u_equation'],extracted['c_equation']):
        assert all(isinstance(value,__import__('fractions').Fraction)for value in polynomial.terms.values())
    assert leading.derive_second_order(symbolic_model) is extracted

def test_background_uses_original_common_shape_frequency_and_dimensionless_amplitude(background,source_config):
    reference=json.loads((ROOT/source_config['linear_reference_bundle']/'result.json').read_text(encoding='utf8'))
    assert background.omega==reference['timoshenko']['roots'][0]['omega']
    assert background.h0==source_config['material_geometry']['h']
    np.testing.assert_allclose(background.evaluate([.5])[0,0],background.h0,rtol=2e-13)
    np.testing.assert_allclose(background.evaluate([0.,1.]),0.,rtol=0,atol=2e-13)
    grid=np.linspace(0,1,101);physical=background.evaluate(grid)
    big=.05*physical;small=.025*physical
    np.testing.assert_array_equal(big,2*small)
    assert background.T1==2*np.pi/background.omega
    assert 'epsilon_a=A/h0' in background.as_dict()['amplitude_convention']

@pytest.mark.parametrize('fraction',(0.,.125,.25,.5))
def test_forcing_constant_twice_frequency_decomposition_from_protected_residual(
        fraction,background,coefficients,symbolic_model):
    points=np.array([0.,.13,.5,.82,1.]);field=background.evaluate(points);first=background.evaluate(points,1);second=background.evaluate(points,2)
    t=fraction*background.T1;cos=np.cos(background.omega*t);sin=np.sin(background.omega*t)
    derived=leading.derive_second_order(symbolic_model);p=coefficients
    W,Theta=field.T;Ws,Ts=first.T;Wss=second[:,0]
    Zs=(p.C-p.S)*(Ts*Ws+Theta*Wss)+(2*p.S-p.C)*Theta*Ts
    D=Theta*Ws-Theta**2/2
    Fu0=Zs/2;Fu2=Zs/2
    Fc0=(-p.nu*p.C*D+p.jp*background.omega**2*Theta**2)/2
    Fc2=(-p.nu*p.C*D-p.jp*background.omega**2*Theta**2)/2
    for index in range(len(points)):
        q=np.zeros(7);qs=np.zeros(7);qss=np.zeros(7);qt=np.zeros(7)
        q[[1,5]]=field[index]*cos;qs[[1,5]]=first[index]*cos;qss[[1,5]]=second[index]*cos;qt[[1,5]]=-background.omega*field[index]*sin
        jet=rod.FieldJet(q,qs,qt,qss,np.zeros(7),np.zeros(7));values=jet.values()|p.values()
        expected_u=Fu0[index]+Fu2[index]*np.cos(2*background.omega*t)
        expected_c=Fc0[index]+Fc2[index]*np.cos(2*background.omega*t)
        np.testing.assert_allclose(derived['forcing_u'].evaluate(values),expected_u,rtol=2e-12,atol=2e-15)
        np.testing.assert_allclose(derived['forcing_c'].evaluate(values),expected_c,rtol=2e-12,atol=2e-15)
    if fraction==.25:
        np.testing.assert_allclose(Fc0-Fc2,p.jp*background.omega**2*Theta**2,rtol=2e-12,atol=2e-15)
        assert np.max(Fc0-Fc2)>0

def test_weak_axial_source_minus_sign_and_strong_assembly_agree(forced):
    nq=4*forced.p+9;nodes,weights=np.polynomial.legendre.leggauss(nq);points=(nodes+1)*forced.length/2;weights*=forced.length/2
    field=forced.background.evaluate(points);first=forced.background.evaluate(points,1)
    p=forced.coefficients;Theta=field[:,1];Ws=first[:,0]
    Z=(p.C-p.S)*Theta*Ws+(p.S-p.C/2)*Theta**2
    expected=-forced.basis_at(points,1)['u'].T@(weights*Z/2)
    weak=forced.assemble_forcing(nq);strong=forced.assemble_forcing(nq,strong=True)
    for source in weak:np.testing.assert_allclose(source[:forced.n],expected,rtol=2e-12,atol=2e-15)
    for a,b in zip(weak,strong):np.testing.assert_allclose(a,b,rtol=2e-12,atol=2e-15)

def test_mass_stiffness_match_frozen_mh_block_and_retain_all_coordinates(forced):
    n=forced.n;assert forced.ndof==2*(forced.p-1)
    assert forced.vectors.shape==(forced.ndof,forced.ndof)
    assert len(forced.omega)==len(forced.b0)==len(forced.b2)==forced.ndof
    assert forced.counters()['modal_reduction'] is False
    assert forced.counters()['time_integrations']==0
    assert forced.eigen_decompositions==1
    ids=np.r_[forced.disc._indices['u'],forced.disc._indices['c']]
    np.testing.assert_allclose(forced.M,forced.disc.M0[np.ix_(ids,ids)],rtol=2e-12,atol=2e-12)
    np.testing.assert_allclose(forced.K,forced.disc.K[np.ix_(ids,ids)],rtol=2e-12,atol=2e-12)
    np.testing.assert_allclose(forced.M,forced.M.T,rtol=0,atol=2e-12)
    np.testing.assert_allclose(forced.K,forced.K.T,rtol=0,atol=2e-12)
    np.linalg.cholesky(forced.M);np.linalg.cholesky(forced.K)
    np.testing.assert_allclose(forced.vectors.T@forced.M@forced.vectors,np.eye(forced.ndof),rtol=2e-12,atol=2e-12)
    residual=forced.K@forced.vectors-forced.M@forced.vectors*forced.eigenvalues
    scale=np.linalg.norm(forced.K@forced.vectors)+np.linalg.norm(forced.M@forced.vectors*forced.eigenvalues)
    assert np.linalg.norm(residual)/scale<2e-12

def test_exact_time_zero_initial_conditions_and_initial_acceleration(forced):
    np.testing.assert_array_equal(forced.evaluate([0.]),0.)
    np.testing.assert_array_equal(forced.evaluate([0.],1),0.)
    np.testing.assert_allclose(forced.M@forced.evaluate([0.],2)[0],forced.f0+forced.f2,rtol=2e-12,atol=2e-15)
    physical=forced.reconstruct(forced.evaluate([.23])[0],[0.,1.])
    np.testing.assert_array_equal(physical,0.)

def test_exact_time_against_independent_augmented_matrix_exponential(forced):
    n=forced.ndof;matrix=np.zeros((2*n+3,2*n+3));matrix[:n,n:2*n]=np.eye(n)
    matrix[n:2*n,:n]=-np.linalg.solve(forced.M,forced.K)
    matrix[n:2*n,2*n]=np.linalg.solve(forced.M,forced.f0)
    matrix[n:2*n,2*n+1]=np.linalg.solve(forced.M,forced.f2)
    matrix[2*n+1,2*n+2]=-forced.driving_omega;matrix[2*n+2,2*n+1]=forced.driving_omega
    initial=np.r_[np.zeros(2*n),1.,1.,0.]
    times=np.array([0.,.001,.13,.73,1.9])
    states=np.asarray([expm(matrix*t)@initial for t in times])
    for derivative,expected in((0,states[:,:n]),(1,states[:,n:2*n]),(2,(states@matrix.T)[:,n:2*n])):
        actual=forced.evaluate(times,derivative)
        scale=max(np.linalg.norm(expected),1e-30)
        assert np.linalg.norm(actual-expected)/scale<2e-11
        np.testing.assert_allclose(actual,expected,rtol=2e-11,atol=2e-13*scale)

def test_forced_energy_has_power_identity_not_conservation_requirement(forced):
    times=np.array([0.,.01,.3,1.]);q=forced.evaluate(times);v=forced.evaluate(times,1);a=forced.evaluate(times,2);f=forced.forcing(times)
    mass=a@forced.M.T;stiffness=q@forced.K.T;residual=mass+stiffness-f
    scale=np.linalg.norm(mass,axis=1)+np.linalg.norm(stiffness,axis=1)+np.linalg.norm(f,axis=1)
    assert np.max(np.linalg.norm(residual,axis=1)/np.maximum(scale,1e-30))<2e-12
    rate=np.einsum('ti,ti->t',v,mass+stiffness);power=np.einsum('ti,ti->t',v,f)
    power_scale=np.maximum(np.linalg.norm(v,axis=1)*scale,1e-30)
    assert np.max(abs(rate-power)/power_scale)<2e-12
    energy=(np.einsum('ti,ij,tj->t',v,forced.M,v)+np.einsum('ti,ij,tj->t',q,forced.K,q))/2
    assert energy[0]==0 and energy[-1]>0

def test_analytic_forcing_quadrature_converges_without_another_eigendecomposition(forced):
    count=forced.eigen_decompositions
    report=forced.forcing_checks((2*forced.p+1,3*forced.p+7,4*forced.p+9))
    assert report['quadrature_is_exact_for_matrices_not_analytic_background'] is True
    for row in report['rows']:
        assert max(row['strong_weak_relative'].values())<2e-12
        if row['relative_changes']:assert max(row['relative_changes'].values())<2e-12
    assert forced.eigen_decompositions==count

def test_amplitude_rescaling_reuses_one_leading_solution_without_any_solve(forced,monkeypatch):
    count=forced.eigen_decompositions;times=np.array([0.,.013,.41]);coordinates=forced.evaluate(times)
    def forbidden(*args,**kwargs):raise AssertionError('Amplitude rescaling attempted another solve')
    monkeypatch.setattr(leading,'eigh',forbidden)
    big=.05**2*coordinates;small=.025**2*coordinates
    np.testing.assert_array_equal(big,4*small)
    assert forced.eigen_decompositions==count

def test_each_leading_equation_and_forcing_monomial_has_expected_units(symbolic_model):
    extracted=leading.derive_second_order(symbolic_model)
    units={'m':(1,-1,0),'jp':(1,1,0),'C':(1,1,-2),'S':(1,1,-2),'H':(1,3,-2),'nu':(0,0,0)}
    shifts={'':(0,0,0),'_s':(0,-1,0),'_t':(0,0,-1),'_ss':(0,-2,0),'_st':(0,-1,-1),'_tt':(0,0,-2)}
    for name,base in(('u',(0,1,0)),('w',(0,1,0)),('theta',(0,0,0)),('c',(0,0,0))):
        for suffix,shift in shifts.items():units[name+suffix]=tuple(a+b for a,b in zip(base,shift))
    for key,target in(('u_equation',(1,0,-2)),('forcing_u',(1,0,-2)),('c_equation',(1,1,-2)),('forcing_c',(1,1,-2))):
        for monomial in extracted[key].terms:
            actual=tuple(sum(units[rod.SYMBOL_ORDER[index]][axis]for index in monomial)for axis in range(3))
            assert actual==target,(key,monomial,actual)

def test_second_order_keeps_original_continuous_initial_acceleration_mismatch(background,coefficients):
    first=background.evaluate([0.,background.length],1)
    trace=(coefficients.C-coefficients.S)/coefficients.m*first[:,1]*first[:,0]
    assert trace[0]>0 and trace[1]<0
    physical=.05**2*trace
    path=ROOT/'results/weakly_nonlinear_planar_recovery/054874a4a4c9c9ff/initial_compatibility.json'
    if path.is_file():
        old=json.loads(path.read_text(encoding='utf8'))
        expected=np.array([row['axial_A2_expression']for row in old['records']if row['A']==.0025])
        np.testing.assert_allclose(physical,expected,rtol=2e-12,atol=2e-16)
    # Halving epsilon rescales this unchanged trace; no endpoint source correction.
    np.testing.assert_array_equal(physical,4*(.025**2*trace))

@pytest.mark.parametrize('change',('amplitude','sampling'))
def test_cache_identity_invalidates_amplitude_and_sampling(tmp_path,monkeypatch,change):
    from scripts.analysis import verify_planar_second_order_axial_response as cli
    config=json.loads(cli.CONFIG.read_text(encoding='utf8'));path=tmp_path/'config.json'
    path.write_text(json.dumps(config),encoding='utf8')
    real_sha=cli.sha
    monkeypatch.setattr(cli,'sha',lambda item:real_sha(item)if Path(item).resolve()==path.resolve()else 'stable-code-or-history-hash')
    before,before_identity=cli.identity(path)
    if change=='amplitude':config['amplitude_over_h'][0]=.04
    else:config['sampling']['samples_per_fastest_retained_period']+=1
    path.write_text(json.dumps(config),encoding='utf8')
    after,after_identity=cli.identity(path)
    assert before!=after
    assert before_identity['code_hashes']==after_identity['code_hashes']
    assert before_identity['historical_manifests']==after_identity['historical_manifests']
    assert before_identity['config']!=after_identity['config']

def test_cached_artifacts_require_their_own_hashes_and_matching_identity(tmp_path):
    from scripts.analysis import verify_planar_second_order_axial_response as cli
    summary={'statuses':{'NLSP_SECOND_ORDER_AXIAL_DIAGNOSTIC':'COMPLETE'},'counters':{'new_ODE_integrations':0}}
    cli.write_json(tmp_path/'summary.json',summary)
    identity={'fixture':'exact-time'}
    cli.write_json(tmp_path/'manifest.json',{'identity':identity,'artifact_hashes':{'summary.json':cli.sha(tmp_path/'summary.json')}})
    assert cli.validate_cache(tmp_path,identity)==summary
    with pytest.raises(ValueError,match='identity mismatch'):cli.validate_cache(tmp_path,{'fixture':'different'})
    (tmp_path/'summary.json').write_text('{}',encoding='utf8')
    with pytest.raises(ValueError,match='artifact hash mismatch'):cli.validate_cache(tmp_path)

@pytest.mark.parametrize('action',('compute','report-only','plot-only'))
def test_cached_entrypoints_have_zero_ode_eigen_derivation_and_analytic_work(tmp_path,monkeypatch,capsys,action):
    from scripts.analysis import verify_planar_second_order_axial_response as cli
    import scipy.linalg
    import scipy.integrate
    output=tmp_path/'results';bundle=output/'fixture';bundle.mkdir(parents=True)
    identity={'synthetic':'cached-exact-time'}
    fields={name:{'relative_L2':.003,'relative_max':.004}for name in cli.COMPONENTS}
    summary={'statuses':{'NLSP_SECOND_ORDER_AXIAL_DIAGNOSTIC':'COMPLETE'},'background':{'T1':2.},
             'convergence':{'p24_p32':{'pair':[24,32],'fields':fields}}}
    cli.write_json(bundle/'summary.json',summary)
    (bundle/'models').mkdir()
    times=np.linspace(0,10,21);values=np.stack([np.sin(times),np.cos(times)],axis=-1)[:,None,:]*.001
    cli.save_npz(bundle/'models/p24_short_observations.npz',time=times,q=values,velocity=values)
    cli.write_json(bundle/'manifest.json',{'identity':identity,'artifact_hashes':
        {str(path.relative_to(bundle)):cli.sha(path)for path in bundle.rglob('*')if path.is_file()}})
    def forbidden(*args,**kwargs):raise AssertionError('Cached entrypoint attempted new numerical work')
    monkeypatch.setattr(cli,'run_compute',forbidden)
    monkeypatch.setattr(leading,'SecondOrderAxial',forbidden)
    monkeypatch.setattr(leading,'derive_second_order',forbidden)
    monkeypatch.setattr(leading,'background_from_pilot',forbidden)
    monkeypatch.setattr(leading,'eigh',forbidden)
    monkeypatch.setattr(rod,'derive_polynomials',forbidden)
    monkeypatch.setattr(scipy.linalg,'eigh',forbidden)
    monkeypatch.setattr(scipy.linalg,'expm',forbidden)
    monkeypatch.setattr(scipy.integrate,'solve_ivp',forbidden)
    monkeypatch.setattr(scipy.integrate,'Radau',forbidden)
    monkeypatch.setattr(cli,'identity',lambda *args:('fixture',identity))
    argv=['--compute','--output-dir',str(output)]if action=='compute'else['--'+action,str(bundle)]
    assert cli.main(argv)==0
    result=json.loads(capsys.readouterr().out)
    for name in('new_time_integrations','eigendecompositions','analytic_evaluations','derivations'):assert result[name]==0
    if action=='plot-only':
        assert (bundle/'figures/leading_response.pdf').is_file()
        assert (bundle/'figures/spatial_convergence.png').is_file()


def test_report_and_plot_dispatch_precede_any_compute_dispatch():
    from scripts.analysis import verify_planar_second_order_axial_response as cli
    tree=ast.parse(Path(cli.__file__).read_text(encoding='utf8'))
    main=next(node for node in tree.body if isinstance(node,ast.FunctionDef)and node.name=='main')
    gate=next(node for node in main.body if isinstance(node,ast.If)and 'args.report_only' in ast.unparse(node.test))
    assert any(isinstance(node,ast.Return)for node in gate.body)
    assert all(not(isinstance(node,ast.Call)and isinstance(node.func,ast.Name)and node.func.id=='run_compute')for node in ast.walk(gate))
    calls=[node for node in ast.walk(main)if isinstance(node,ast.Call)and isinstance(node.func,ast.Name)and node.func.id=='run_compute']
    assert calls and gate.lineno<min(node.lineno for node in calls)


def _saved_response_arrays(forced):
    return {'M':forced.M.copy(),'K':forced.K.copy(),'f0':forced.f0.copy(),'f2':forced.f2.copy(),
            'omega':forced.omega.copy(),'vectors':forced.vectors.copy(),'b0':forced.b0.copy(),'b2':forced.b2.copy(),
            'transform_u':forced.transforms[0].copy(),'transform_c':forced.transforms[1].copy()}


def test_saved_response_roundtrip_retains_every_coordinate_without_an_eigensolve(forced,tmp_path,monkeypatch):
    path=tmp_path/'saved.npz';np.savez_compressed(path,**_saved_response_arrays(forced))
    def forbidden(*args,**kwargs):raise AssertionError('Restore invoked an eigendecomposition')
    monkeypatch.setattr(leading,'eigh',forbidden)
    monkeypatch.setattr(planar,'eigh',forbidden)
    restored=leading.SecondOrderAxial.from_saved(forced.coefficients,forced.p,forced.background,path,model=forced.model)
    assert restored.eigen_decompositions==0
    assert restored.counters()['time_integrations']==0
    assert restored.ndof==forced.ndof==2*(forced.p-1)
    assert restored.vectors.shape==(forced.ndof,forced.ndof)
    for name in('M','K','f0','f2','omega','vectors','b0','b2'):
        np.testing.assert_array_equal(getattr(restored,name),getattr(forced,name))
    times=np.array([0.,.001,.13,1.97])
    for order in(0,1,2):
        np.testing.assert_array_equal(restored.evaluate(times,order),forced.evaluate(times,order))
    np.testing.assert_array_equal(restored.reconstruct_series(restored.evaluate(times),[0.,.25,.5,1.]),
                                  forced.reconstruct_series(forced.evaluate(times),[0.,.25,.5,1.]))
    assert restored.restore_checks['eigendecompositions']==0


@pytest.mark.parametrize('corruption',('K','f0','b0'))
def test_saved_response_rejects_corrupted_operator_or_forcing_coefficients(forced,tmp_path,monkeypatch,corruption):
    arrays=_saved_response_arrays(forced)
    index=np.unravel_index(np.argmax(abs(arrays[corruption])),arrays[corruption].shape)
    arrays[corruption][index]*=1.01
    path=tmp_path/'corrupt.npz';np.savez_compressed(path,**arrays)
    def forbidden(*args,**kwargs):raise AssertionError('Invalid cache attempted a new eigensolve')
    monkeypatch.setattr(leading,'eigh',forbidden)
    with pytest.raises(ValueError,match='differs from declared model/background'):
        leading.SecondOrderAxial.from_saved(forced.coefficients,forced.p,forced.background,path,model=forced.model)


def test_saved_response_rejects_different_declared_physical_coefficient(forced,tmp_path,monkeypatch):
    path=tmp_path/'saved.npz';np.savez_compressed(path,**_saved_response_arrays(forced))
    values=forced.coefficients.values();values['C']*=1.001
    changed=rod.RodCoefficients(**values)
    def forbidden(*args,**kwargs):raise AssertionError('Mismatched cache attempted a new eigensolve')
    monkeypatch.setattr(leading,'eigh',forbidden)
    with pytest.raises(ValueError,match='differs from declared model/background'):
        leading.SecondOrderAxial.from_saved(changed,forced.p,forced.background,path,model=forced.model)


def _conditional_coverage_fixture(cli,output):
    base={'version':'synthetic-same-model','config':{'degrees':[16,24,32,48,64],
          'amplitude_convention':'epsilon_a=A/h0'},'code_hashes':{'helper':'accepted-code'}}
    child=output/'conditional-child';child.mkdir(parents=True)
    summary={'statuses':{'NLSP_SECOND_ORDER_AXIAL_DIAGNOSTIC':'COMPLETE'},
             'models':{str(p):{}for p in(16,24,32,48,64,96)}}
    cli.write_json(child/'summary.json',summary)
    cli.write_json(child/'manifest.json',{'identity':{**base,'conditional_p96_parent':{'manifest_sha256':'preserved-parent'}},
        'artifact_hashes':{'summary.json':cli.sha(child/'summary.json')}})
    cli.write_json(output/'current.json',{'bundle':str(child),'fingerprint':'conditional-child'})
    return base,child


def test_current_conditional_coverage_reuses_primary_grid_without_new_numerical_work(tmp_path,monkeypatch,capsys):
    from scripts.analysis import verify_planar_second_order_axial_response as cli
    output=tmp_path/'results';base,child=_conditional_coverage_fixture(cli,output)
    monkeypatch.setattr(cli,'identity',lambda *args:('base-not-yet-on-disk',base))
    def forbidden(*args,**kwargs):raise AssertionError('Conditional coverage recomputed an existing primary grid')
    monkeypatch.setattr(cli,'run_compute',forbidden)
    monkeypatch.setattr(cli,'run_extend_p96',forbidden)
    monkeypatch.setattr(leading,'SecondOrderAxial',forbidden)
    monkeypatch.setattr(leading,'derive_second_order',forbidden)
    monkeypatch.setattr(leading,'eigh',forbidden)
    monkeypatch.setattr(rod,'derive_polynomials',forbidden)
    assert cli.main(['--compute','--output-dir',str(output)])==0
    result=json.loads(capsys.readouterr().out)
    assert result['cache']=='reused_conditional_coverage'
    assert Path(result['bundle'])==child
    for name in('new_time_integrations','eigendecompositions','analytic_evaluations','derivations'):assert result[name]==0
    assert not(output/'base-not-yet-on-disk').exists()


def test_normalization_or_code_change_invalidates_conditional_coverage(tmp_path,monkeypatch):
    from scripts.analysis import verify_planar_second_order_axial_response as cli
    import copy
    output=tmp_path/'results';base,_=_conditional_coverage_fixture(cli,output)
    config_path=tmp_path/'config.json';cli.write_json(config_path,base['config'])
    calls=[]
    def sentinel(config,bundle):
        calls.append(str(bundle))
        raise RuntimeError('CACHE_INVALIDATED_NO_ACTUAL_COMPUTE')
    monkeypatch.setattr(cli,'run_compute',sentinel)
    for changed_property in('normalization','code'):
        expected=copy.deepcopy(base)
        if changed_property=='normalization':expected['config']['amplitude_convention']='different-reference-scale'
        else:expected['code_hashes']['helper']='changed-code'
        monkeypatch.setattr(cli,'identity',lambda *args,expected=expected:('invalidated-'+changed_property,expected))
        with pytest.raises(RuntimeError,match='CACHE_INVALIDATED_NO_ACTUAL_COMPUTE'):
            cli.main(['--compute','--config',str(config_path),'--output-dir',str(output)])
    assert len(calls)==2


@pytest.mark.parametrize("reference_kind", ("audit_bundle", "linear_reference_bundle"))
def test_identity_invalidates_changed_reference_result(tmp_path, monkeypatch, reference_kind):
    from scripts.analysis import verify_planar_second_order_axial_response as cli
    config = json.loads(cli.CONFIG.read_text(encoding="utf8"))
    path = tmp_path/"config.json"
    path.write_text(json.dumps(config), encoding="utf8")
    pilot = cli.read_json(cli.ROOT/config["pilot_config"])
    target = (cli.ROOT/pilot[reference_kind]/"result.json").resolve()
    real_sha, changed = cli.sha, [False]
    def source_hash(item):
        item = Path(item).resolve()
        if item == path.resolve():
            return real_sha(item)
        return "changed-reference-result" if changed[0] and item == target else "stable-source-or-code"
    monkeypatch.setattr(cli, "sha", source_hash)
    before, old = cli.identity(path)
    changed[0] = True
    after, new = cli.identity(path)
    assert before != after
    assert old["code_hashes"] == new["code_hashes"]
    assert old["reference_source_hashes"][reference_kind]["result"] != new["reference_source_hashes"][reference_kind]["result"]
