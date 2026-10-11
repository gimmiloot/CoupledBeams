"""Modal algebra and saved source checks; no physical ODE/FEM in tests."""
from types import SimpleNamespace
import numpy as np
import pytest
from scripts.lib import nlsp_resonant_spatial_coupling_screening as s
from scripts.lib import nlsp_coupled_sensitivity_analysis as a

def test_scope_and_zero_ODE_budget():
    config=s.p.read(s.CONFIG);assert s.validate_scope(config)==config
    assert config['new_ODE_integrations']==0 and config['numerical_budget_seconds']==900
    for key,value in [('degree',16),('amplitude',.01),('mu_points',[.1,.2,.3]),('automatic_retry',True)]:
        with pytest.raises(ValueError):s.validate_scope({**config,key:value})

def test_modal_energy_unit_duality_and_conjugate_branches():
    omega=np.array([.2,.7]);M=np.diag([2.,3.]);phi=np.diag(1/np.sqrt(np.diag(M)));K=M@np.diag(omega**2)
    R,L,E=s.state_modes(K,M,omega,phi);A=np.block([[np.zeros((2,2)),np.eye(2)],[-np.linalg.solve(M,K),np.zeros((2,2))]])
    assert np.allclose(L@R,np.eye(4),atol=2e-15)
    assert np.array_equal(R[:,2:],R[:,:2].conj())
    assert np.allclose(A@R,R*np.r_[1j*omega,-1j*omega],atol=2e-15)

def test_harmonic_extraction_convention():
    phase=np.arange(32)*2*np.pi/32;values=(3+2*np.cos(phase)+4*np.sin(2*phase))[:,None,None]
    H=s.harmonics(values,phase,[-2,-1,0,1,2])[:,0,0]
    assert np.allclose(H,[2j,1,3,1,-2j],atol=3e-15)

def test_even_odd_amplitude_and_inverse_mass_higher_order_removal():
    amplitudes=np.array([.005,.0025,.00125]);orders=[-2,-1,0,1,2];powers=np.array([2,1,2,1,2])
    leading=np.arange(1,6)[:,None,None];higher=(np.arange(1,6)*2)[:,None,None]
    H=amplitudes[:,None,None,None]**powers[None,:,None,None]*leading+amplitudes[:,None,None,None]**(powers[None,:,None,None]+2)*higher
    first,second,p=s.amplitude_orders(H,amplitudes,orders)
    assert np.allclose(first,leading,atol=2e-15) and np.allclose(second,leading,atol=2e-15)

def test_negative_frequency_coefficient_is_not_simple_complex_conjugation():
    omega=np.array([.2]);R,L,E=s.state_modes(np.diag(omega**2),np.eye(1),omega,np.eye(1));phase=np.arange(32)*2*np.pi/32
    values=np.array([L@np.array([[0.,0.],[-np.cos(t),0.]])@R for t in phase]);H=s.harmonics(values,phase,[-1,0,1]);swap=[1,0]
    assert np.allclose(H[0],H[2][np.ix_(swap,swap)].conj(),atol=2e-15)
    assert not np.allclose(H[0],H[2].conj())

def test_difference_exchange_vs_sum_direct_growth_algebra():
    config=s.p.read(s.CONFIG);k=2;H=np.zeros((5,4,4),complex);orders=config['harmonics']
    H[orders.index(1),1,0]=.1j;H[orders.index(-1),0,1]=.1j
    H[orders.index(1),0,2]=-.1j;H[orders.index(-1),2,0]=.1j
    rows=s.pair_rows({'omega':np.array([.5,1.5]),'Omega':1.,'mu':.2},H,np.zeros_like(H.real),config)
    diff=next(r for r in rows if r['kind']=='difference' and r['harmonic']==1)
    self0=next(r for r in rows if r['kind']=='self' and r['harmonic']==1 and r['i_sorted']==1)
    assert diff['product']['real']<0 and not diff['direct_term_candidate'] and diff['classification']=='difference_exchange'
    assert self0['product']['real']>0 and self0['direct_term_candidate']

@pytest.mark.parametrize('mode',['--compute','--report-only','--plot-only'])
def test_cached_replay_no_physical_calls(monkeypatch,tmp_path,mode):
    summary={'completed':True,'hard_stop':False,'overall':'DONE','stages':{},'calls':{},'numerical_seconds':1.}
    monkeypatch.setattr(s,'cache',lambda *x:summary);monkeypatch.setattr(s,'prepare',lambda *x:(tmp_path,{},summary))
    monkeypatch.setattr(s.p,'sources',lambda *x:pytest.fail('source reconstruction'))
    monkeypatch.setattr(s,'modal_problem',lambda *x:pytest.fail('new eigensolve'))
    monkeypatch.setattr(s,'projected_at',lambda *x:pytest.fail('new Jacobian'))
    monkeypatch.setattr(s.p,'integrate_attempt',lambda *x,**kw:pytest.fail('new nonlinear ODE'))
    monkeypatch.setattr(a,'tangent_attempt',lambda *x,**kw:pytest.fail('new tangent ODE'))
    plot_calls=[]
    monkeypatch.setattr(s,'plot_saved',lambda b:plot_calls.append(b))
    monkeypatch.setattr(s.p,'read',lambda *x:{})
    monkeypatch.setattr(s.p,'finalize',lambda *x:None)
    s.main([mode] if mode=='--compute' else [mode,str(tmp_path)])
    assert len(plot_calls)==(1 if mode=='--plot-only' else 0)

def test_source_hashes_and_frozen_kernels():
    config=s.p.read(s.CONFIG)
    for key in ('modal_parent','longer_parent'):assert s.p.sha(s.p.ROOT/config[key]/'manifest.json')==config[key+'_manifest_sha256']
    parent=s.p.read(s.p.ROOT/config['modal_parent']/'provenance.json')
    assert all(s.p.sha(s.p.ROOT/f)==parent['code_hashes'][f] for f in s.r.KERNELS)

@pytest.fixture(scope='module')
def problems():
    config=s.p.read(s.CONFIG);source=s.p.read(s.p.ROOT/config['modal_parent']/'provenance.json');model=s.p.sources(source['physical_source_config'])[0]
    calls={'linear_eigenproblems':0,'pointwise_Jacobian':0}
    return {mu:s.modal_problem(model,mu,config,calls) for mu in (.2,0.)}

def test_actual_modal_projection_includes_variable_mass_and_full_jacobian(problems,monkeypatch):
    monkeypatch.setattr(s.p,'integrate_attempt',lambda *x,**kw:pytest.fail('nonlinear ODE'))
    monkeypatch.setattr(a,'tangent_attempt',lambda *x,**kw:pytest.fail('tangent ODE'))
    problem=problems[.2];d=problem['d'];q=problem['q'];x=np.r_[q*.7,-q*problem['Omega']*.3];counter={'pointwise_Jacobian':0}
    actual=s.projected_at(problem,x,counter);J=d.jacobian(0,x)[np.ix_(d.state_out,d.state_out)]
    assert np.allclose(actual,problem['L']@(J-problem['A0'])@problem['R'],atol=1e-14)
    assert problem['duality']<2e-10 and problem['eigen_residual']<2e-8
    assert s.p.polynomial_w_max(d,q)[0]==pytest.approx(.005,rel=2e-14)
    assert len(problem['guard_omega'])==1

def test_real_equal_arm_mirror_selection_rule(problems):
    config=s.p.read(s.CONFIG);counter={'pointwise_Jacobian':0};problem=problems[0.]
    H,_=s.harmonic_sample(problem,1.,32,config['harmonics'],counter,1e30)
    result=s.mirror_selection(problem,H,config,counter)
    assert result['status']=='PASS' and result['planar_parity']==pytest.approx(-1,abs=1e-12)
    assert result['coefficients']['1']['forbidden_pairs_max_coefficient']<1e-11
    assert result['coefficients']['1']['allowed_pairs_max_coefficient']>1e-4

def test_harmonic_phase_resolution_and_source_geometry(problems):
    problem=problems[.2];config=s.p.read(s.CONFIG);counter={'pointwise_Jacobian':0}
    full,values=s.harmonic_sample(problem,1.,32,config['harmonics'],counter,1e30)
    coarse=s.harmonics(values[::2],np.arange(16)*2*np.pi/16,config['harmonics'])
    assert np.max(abs(full-coarse))<2e-10
    assert problem['d'].lengths==(.8,1.2) and all(arm.coefficients.Bp/arm.coefficients.Bb==pytest.approx(4) for arm in problem['d'].arms)

def test_unresolved_mode_matching_not_forced(problems):
    mac=s.physical_MAC(problems[.2],problems[.2]);assert np.allclose(np.diag(mac),1,atol=2e-15)
    config=s.p.read(s.CONFIG)
    assert config['correspondence_MAC_gate']==.8
