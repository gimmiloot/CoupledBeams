"""Targeted source/sign/scaling tests; no repeated literature spectrum."""
import numpy as np
import pytest
from scripts.lib import inplane_kelvin_voigt_literature_benchmarks as lit
from scripts.lib import inplane_kelvin_voigt as kv


def test_failla_time_and_dimensional_mapping():
    L,m,D=2.,3.,5.
    time=lit.reference_time(L,m,D)
    omega=21.9+.263j
    z=lit.failla_to_project(omega)
    assert z==-0.263+21.9j
    assert time==pytest.approx(4*np.sqrt(3/5))
    assert np.exp(1j*(omega/time)*.12)==pytest.approx(np.exp(z*.12/time))
    cu,cr=.7,.8
    assert cu*L/np.sqrt(m*D)==pytest.approx(cu*L**3/(D*time))
    assert cr/(L*np.sqrt(m*D))==pytest.approx(cr*L/(D*time))


@pytest.mark.parametrize('ku,gu,kr,gr',[(3.,.2,None,0.),(0.,0.,4.,.1),(3.,.2,4.,.1)])
def test_failla_interfaces(ku,gu,kr,gr):
    z=-.2+3j;minus=np.array([.1+.2j,.4,-.3j,2.])
    plus=lit.failla_jump(z,ku,gu,kr,gr)@minus
    expected=minus.copy()
    if kr is not None:expected[1]-=minus[2]/(kr+gr*z)
    expected[3]+=(ku+gu*z)*minus[0]
    np.testing.assert_allclose(plus,expected)
    L,R=lit.failla_interface(z,ku,gu,kr,gr)
    np.testing.assert_allclose(L@minus+R@plus,0,atol=1e-14)


def test_inactive_bare_mode_is_not_assigned_damping():
    beta=4*np.pi;z=1j*beta**2
    for xi in (.25,.5,.75):
        state=np.array([np.sin(beta*xi),beta*np.cos(beta*xi),
                        beta**2*np.sin(beta*xi),beta**3*np.cos(beta*xi)],complex)
        np.testing.assert_allclose(lit.failla_jump(z)@state,state,atol=1e-11)
    # A positive damping coefficient remains present in the operator.
    Lz,Rz=lit.failla_interface(z,derivative=True)
    assert Lz[1,1]==-.1 and Rz[1,1]==.1 and Lz[3,0]==-.1


def test_hong_geometry_shape_factor_and_independent_moduli():
    h=lit.Hong()
    assert h.A==pytest.approx(.000625)
    assert h.I==pytest.approx(3.255208333333334e-8)
    assert h.K==pytest.approx(13/15.3)
    assert h.G==80e9 and h.G!=h.E/(2*(1+h.nu))


@pytest.mark.parametrize('s',[-.3+20j,1.+13j,-2.+800j])
def test_hong_source_matrix_and_project_mapping(s):
    h=lit.Hong()
    expected=np.array([[0,1,-1/(h.K*h.A*h.G),0],
       [0,0,0,1/(h.E*h.I)],[-h.rho*h.A*s*s,0,0,0],
       [0,h.rho*h.I*s*s,1,0]],complex)
    np.testing.assert_allclose(lit.hong_state(s),expected,rtol=kv.CRITERIA['H_rtol'],atol=0.)
    arm=kv.Arm('RLB',h.E*h.A,h.D,h.m,h.L,1/h.S,h.J)
    H=kv.state_matrix(s,arm)[np.ix_([1,2,4,5],[1,2,4,5])]
    S=np.diag([1,-1,-1,-1])
    np.testing.assert_array_equal(lit.hong_state(s),S@H@S)
    assert lit.hong_state(s)[2,0].imag!=0


def test_failla_project_mapping():
    z=-.1+12j
    arm=kv.Arm('EB',1,1,1,1,0,0)
    H=kv.state_matrix(z,arm)[np.ix_([1,2,4,5],[1,2,4,5])]
    P=np.array([[1,0,0,0],[0,-1,0,0],[0,0,0,1],[0,0,1,0]])
    np.testing.assert_array_equal(lit.failla_state(z),P@H@P.T)


@pytest.mark.parametrize('case',['failla','hong_hh','hong_ff','hong_damped'])
def test_boundary_scalar_equations_derivative_and_conjugacy(case):
    beam=lit.Beam(case);z=-.15+24j
    B,Bz=beam.matrices(z,derivative=True)
    step=1e-5
    fd=(beam.matrices(z+step)[0]-beam.matrices(z-step)[0])/(2*step)
    np.testing.assert_allclose(Bz,fd,rtol=1e-7,atol=1e-8)
    np.testing.assert_allclose(beam.matrices(z.conjugate())[0],B.conj(),rtol=1e-12,atol=1e-12)
    rng=np.random.default_rng(316)
    a=rng.normal(size=16)+1j*rng.normal(size=16)
    y,norm_a=beam.recover(z,a)
    scale=np.max(abs(y/beam.scale))
    np.testing.assert_allclose(beam.physical_conditions(z,y),B@norm_a/scale,atol=1e-12,rtol=1e-10)


def test_hong_boundaries_and_support_stiffness():
    s=-3+12j;h=lit.Hong();z=s*h.time
    hh=lit.Beam('hong_hh');ff=lit.Beam('hong_ff');support=lit.Beam('hong_damped')
    L,R=hh.boundary(z)
    np.testing.assert_array_equal(L,np.eye(4)[[0,3]])
    np.testing.assert_array_equal(R,L)
    np.testing.assert_array_equal(ff.boundary(z)[0],np.eye(4)[[2,3]])
    L,R=support.boundary(z)
    K=(2e6+20*s)*h.L**3/h.D
    assert L[0,0]==pytest.approx(K/1000)
    assert R[0,0]==pytest.approx(-K/1000)
    assert L[1,3]==R[1,3]==1


@pytest.mark.parametrize('printed,half',[('21.9037',5e-5),('238.023',5e-4),
    ('3.55778498e+002',5e-7),('-6.6651e-002',5e-7),('0.0120',5e-5)])
def test_published_rounding(printed,half):
    v=float(printed)
    assert lit.rounding(printed,v)['rounding_tolerance']==pytest.approx(half)
    assert lit.rounding(printed,v+.4*half)['rounding_pass']
    assert not lit.rounding(printed,v+3*half)['rounding_pass']
    assert not lit.rounding(printed,np.nan)['rounding_pass']


def test_literal_missing_zero_and_complex_dtype():
    assert lit.rounding(None,1e-11)['rounding_pass']
    assert not lit.rounding(None,1e-4)['rounding_pass']
    assert np.iscomplexobj(lit.Beam('failla').matrices(2j)[0])


def test_d13_inactive_exact_value_and_printed_ratios():
    assert (4*np.pi)**2==pytest.approx(157.91367041742973,abs=1e-13)
    assert format((4*np.pi)**2,'.4f')=='157.9137'
    expected=['0.0120','0.0388','0.0797','0.0','0.0607']
    rows=[lit.printed_zeta(*target) for target in lit.FAILLA]
    assert [r['ordinary_rounded_ratio'] for r in rows]==expected
    assert [r['displayed_consistent'] for r in rows]==[True,True,False,True,True]
    assert rows[2]['implied_ratio']==pytest.approx(6.541/(81.8244**2+6.541**2)**.5)


@pytest.mark.parametrize('case,roots',[
    ('hong_hh',[355.7784981198391,1418.824486875438,3176.508755760719,
                5608.549959938625,8688.066141281099]),
    ('hong_ff',[805.5329396844996,2211.427382729008,4309.664856724562,
                7068.984777032399,10460.217657249696])])
def test_table2_analytical_footnote_values_in_matrix(case,roots):
    # Independent scalar-formula values retained from the K13 check, no new roots.
    beam=lit.Beam(case)
    for omega in roots:
        singular=np.linalg.svd(beam.matrices(1j*omega*beam.time)[0],compute_uv=False)
        assert singular[-1]/singular[0]<kv.CRITERIA['sigma_ratio']


@pytest.mark.parametrize('failure',[None,'equation','mapping','parameters','criteria'])
def test_d13_gate_is_equation_based_and_preserves_print_fail(failure):
    args=dict(equation_status='PASS_EQUATION_LEVEL',mapping_unchanged=True,
              parameters_unchanged=True,criteria_unchanged=True)
    if failure=='equation':args['equation_status']='FAIL_EQUATION_LEVEL'
    elif failure:args[failure+'_unchanged']=False
    assert lit.d13_gate(**args)==(failure is None)


def test_second_pass_hong_targets_mapping_and_separate_statuses():
    expected=[(-.066651,334.44),(-2.7327,1107.9),(-12.133,1927.1),
              (-20.106,2954.2),(-20.135,4711.1)]
    assert [(float(r),float(i)) for r,i in lit.HONG_DAMPED]==expected
    d=dict(status='CONVERGED',solver_residual=1e-16,sigma_ratio=1e-16,
           physical_residual=1e-14,conjugate_residual=1e-16,steps=2,
           last_delta_z=1e-14,left_Bz_right=.1)
    s=complex(*expected[0])
    row=lit.hong_second_pass_comparison(s,d,1)
    assert row['computed_real']==s.real and row['computed_imag']==s.imag
    assert row['printed_rounding_status']=='PRINT_MATCH'
    assert row['equation_solver_status']=='SOLVER_PASS'
    row=lit.hong_second_pass_comparison(s+7e-7,d,1)
    assert row['printed_rounding_status']=='PRINT_MISMATCH'
    assert row['equation_solver_status']=='SOLVER_PASS'
    assert row['last_printed_place_only']  # description, NOT PRINT_MATCH


def test_d13_audit_strict_formula_tolerance_and_csv_consistency():
    failla=[dict(mode=n,printed_value_real=p,printed_value_imag=q or '',
                 printed_value_ratio=r,computed_real=float(p),inactive_confirmed='True')
            for n,(p,q,r) in enumerate(lit.FAILLA,1)]
    hong=[];checks=[]
    for case,targets in [('hong_hh',lit.HONG_HH),('hong_ff',lit.HONG_FF)]:
        for n,target in enumerate(targets,1):
            w=float(target)
            hong.append(dict(case=case,mode=n,computed_imag=w,printed_value_imag=target,
                             status='LITERATURE_MISMATCH',rounding_pass='False'))
            checks.append(dict(case=case,mode=n,matrix_omega=w,formula_omega=w))
    d=dict(local_formula_check=dict(checks=checks))
    a=lit.source_precision_audit(failla,hong,d)
    assert a['hong_printed_status']=='FAIL' and a['hong_equation_status']=='PASS_EQUATION_LEVEL'
    checks[0]['formula_omega']+=1e-8
    assert lit.source_precision_audit(failla,hong,d)['hong_equation_status']=='FAIL_EQUATION_LEVEL'
    checks[0]['matrix_omega']+=1
    with pytest.raises(ValueError,match='CSV'):
        lit.source_precision_audit(failla,hong,d)


def test_second_pass_reuse_no_computation_and_source_mutation(monkeypatch,tmp_path):
    import json
    from scripts.analysis.laminated_beams import benchmark_inplane_kelvin_voigt_literature as run
    monkeypatch.setattr(run,'ROOT',tmp_path)
    runner=tmp_path/'runner.py';runner.write_text('frozen runner')
    helper=tmp_path/'scripts/lib/inplane_kelvin_voigt_literature_benchmarks.py'
    helper.parent.mkdir(parents=True);helper.write_text('frozen helper')
    monkeypatch.setattr(run,'__file__',str(runner))
    protected=tmp_path/'first.csv';protected.write_text('first pass unchanged')
    output=tmp_path/'computed.csv';output.write_text('second pass unchanged')
    versions={str(p.relative_to(tmp_path)):run.sha(p) for p in (runner,helper)}
    manifest=tmp_path/'second_pass_manifest.json'
    manifest.write_text(json.dumps(dict(source_versions=versions,finished=True,
        protected_sources={'first.csv':run.sha(protected)},output_hashes={'computed.csv':run.sha(output)})))
    def forbidden(*args,**kwargs):pytest.fail('completed reuse called computation')
    monkeypatch.setattr(run,'preflight',forbidden)
    monkeypatch.setattr(lit,'solve',forbidden)
    monkeypatch.setattr(lit.Beam,'matrices',forbidden)
    monkeypatch.setattr(lit.Beam,'recover',forbidden)
    monkeypatch.setattr(lit,'source_precision_audit',forbidden)
    before={p.name:p.read_bytes() for p in tmp_path.iterdir() if p.is_file()}
    run.second_pass(tmp_path)
    assert before=={p.name:p.read_bytes() for p in tmp_path.iterdir() if p.is_file()}
    protected.write_text('mutated')
    with pytest.raises(ValueError,match='provenance changed'):run.second_pass(tmp_path)
