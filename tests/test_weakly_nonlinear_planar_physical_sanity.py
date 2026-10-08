"""Limited physical-sanity tests; no extra ODE, BVP or eigensolves.

Formal bulk stretching and amplitude-sign checks are separate from numerical
certification and from full 3D or experimental validation.
"""
from __future__ import annotations
import ast
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from scripts.lib import planar_prepared_initial_state as prep
from scripts.lib import planar_second_order_axial_response as leading
from scripts.lib import weakly_nonlinear_planar_dynamics as planar
from scripts.lib import weakly_nonlinear_spatial_rod as rod

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"results/planar_prepared_one_T1/795dcb14d3cd3a55"
PREPARED=ROOT/"results/planar_prepared_initial_state/5ea8d41faf8ede54"
AUDIT=ROOT/"results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295"
SIGNS=np.array([1.,-1.,-1.,1.])


def read(path):
    return json.loads(Path(path).read_text(encoding="utf8"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture(autouse=True)
def no_extra_integrations_bvp_or_eigensystems(monkeypatch):
    import scipy.integrate
    import scipy.linalg
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as runner
    def forbidden(*args,**kwargs):
        raise AssertionError("Physical-sanity tests may not start ODE/BVP/eigensolver runs")
    for module,names in ((scipy.integrate,("solve_ivp","Radau")),
                         (scipy.linalg,("eigh","eig")),
                         (np.linalg,("eigh","eig","eigvals")),
                         (leading,("eigh",)),(planar,("eigh",)),
                         (runner,("integrate_case",))):
        for name in names:
            monkeypatch.setattr(module,name,forbidden)


@pytest.fixture(scope="module")
def frozen_state():
    return prep.load_frozen_prepared_state(PREPARED)


@pytest.fixture(scope="module")
def frozen_action():
    source=read(AUDIT/"result.json")
    manifest=read(AUDIT/"manifest.json")
    hashes=manifest.get("artifact_hashes",manifest.get("artifacts"))
    assert sha(AUDIT/"result.json")==hashes["result.json"]
    pol=source["polynomials"]
    return SimpleNamespace(T4=rod.Polynomial.deserialize(pol["T4"]),
        V4=rod.Polynomial.deserialize(pol["V4"]),
        residual_a=tuple(rod.Polynomial.deserialize(x) for x in pol["residuals_A"]),
        flux_a=tuple(rod.Polynomial.deserialize(x) for x in pol["fluxes"]),
        symbols={key:rod.Polynomial.symbol(key) for key in rod.SYMBOL_ORDER})


def test_original_physical_helpers_and_coefficients_remain_frozen():
    manifest=read(SOURCE/"manifest.json")
    hashes=manifest["identity"]["code_hashes"]
    for path in ("scripts/lib/weakly_nonlinear_spatial_rod.py",
                 "scripts/lib/weakly_nonlinear_planar_dynamics.py",
                 "scripts/lib/planar_prepared_initial_state.py"):
        assert sha(ROOT/path)==hashes[path]
    old=read(SOURCE/"summary.json")
    assert old["state_admitted_flag"] is False
    assert old["source_strict_status"]=="PARTIAL"
    assert old["execution_mode"]=="EXPLORATORY_NOT_CERTIFIED"


@pytest.mark.parametrize("derivative",(0,1,2))
def test_half_amplitude_evaluator_has_distinct_field_powers_and_fixed_theta3(frozen_state,derivative):
    state,_,_=frozen_state
    x=np.linspace(0.,state.length,17);large=.05;small=.025
    big=state.evaluate(x,large,derivative,require_admitted=False)
    half=state.evaluate(x,small,derivative,require_admitted=False)
    np.testing.assert_array_equal(half[:,[0,3]],big[:,[0,3]]/4)
    np.testing.assert_array_equal(half[:,1],big[:,1]/2)
    expected=big[:,2]/2-(3/8)*large**3*state.correction.evaluate(x,derivative)
    np.testing.assert_allclose(half[:,2],expected,rtol=2e-14,atol=2e-16)
    assert np.linalg.norm(half[:,0]-big[:,0]/2)>0.
    assert np.linalg.norm(half[:,3]-big[:,3]/2)>0.
    assert state.admitted is False


@pytest.mark.parametrize("epsilon",(.05,.025))
@pytest.mark.parametrize("derivative",(0,1,2))
def test_prepared_evaluator_obeys_amplitude_sign_symmetry_not_spatial_tracking(frozen_state,epsilon,derivative):
    state,_,_=frozen_state
    x=np.linspace(0.,state.length,17)
    positive=state.evaluate(x,epsilon,derivative,require_admitted=False)
    negative=state.evaluate(x,-epsilon,derivative,require_admitted=False)
    np.testing.assert_array_equal(negative,positive*SIGNS)
    np.testing.assert_array_equal(state.initial_velocities(x),0.)


def test_quartic_action_and_cubic_residuals_are_exactly_covariant_under_bending_sign(frozen_action):
    model=frozen_action
    mapping={field+suffix:(-model.symbols[field+suffix])
             for field in ("w","theta")
             for suffix in ("","_s","_ss","_t","_st","_tt")}
    # All other planar fields and every out-of-plane jet are held unchanged.
    plane={field+suffix:0 for field in ("v","Phi","psi")
           for suffix in ("","_s","_ss","_t","_st","_tt")}
    for expression in (model.T4,model.V4):
        restricted=expression.substitute(plane)
        assert restricted.substitute(mapping)==restricted
    for index,sign in zip((0,1,5,6),SIGNS.astype(int)):
        restricted=model.residual_a[index].substitute(plane)
        assert restricted.substitute(mapping)==sign*restricted


def test_canonical_translation_fluxes_come_from_same_quartic_energy(frozen_action):
    plane={field+suffix:0 for field in ("v","Phi","psi")
           for suffix in ("","_s","_ss","_t","_st","_tt")}
    energy=frozen_action.V4.substitute(plane)
    for field,index in (("u",0),("w",1),("theta",5),("c",6)):
        assert frozen_action.flux_a[index].substitute(plane)==energy.derivative(field+"_s")


def test_linear_reaction_flux_convention_and_c_R_independence(frozen_action,frozen_state):
    _,coeff,_=frozen_state
    p=frozen_action.symbols
    plane={field+suffix:0 for field in ("v","Phi","psi")
           for suffix in ("","_s","_ss","_t","_st","_tt")}
    expected=(p["C"]*(p["u_s"]+p["nu"]*p["c"]),
              p["S"]*(p["w_s"]-p["theta"]),p["Bp"]*p["theta_s"],p["H"]*p["c_s"])
    for index,e in zip((0,1,5,6),expected):
        assert frozen_action.flux_a[index].substitute(plane).truncate(1)==e
    assert coeff.H>0. and coeff.jp>0.  # independent contraction retained
    # Outward boundary work has sigma=-1,+1, not equality of the two forces.
    right,left=np.array([.2,-.1]),np.array([-.3,.4])
    traction=np.stack((-left,right))
    np.testing.assert_array_equal(traction.sum(axis=0),right-left)
    assert np.linalg.norm(right-left)>0.


@pytest.fixture(scope="module")
def parity_disc(frozen_state,frozen_action):
    _,coeff,_=frozen_state
    return planar.PlanarGalerkin(coeff,8,model=frozen_action)


def test_semidiscrete_rhs_has_same_amplitude_parity_with_variable_mass(parity_disc):
    d=parity_disc
    rng=np.random.default_rng(809)
    q=rng.normal(size=d.ndof)*1e-8
    velocity=rng.normal(size=d.ndof)*1e-9
    sign=np.repeat(SIGNS,d.p-1)
    state=np.r_[q,velocity];reflected=np.r_[sign*q,sign*velocity]
    original=d.rhs(0.,state)
    other=d.rhs(0.,reflected)
    expected=np.r_[sign,sign]*original
    scale=max(np.max(abs(expected)),1e-30)
    np.testing.assert_allclose(other,expected,rtol=2e-12,atol=2e-12*scale)
    np.testing.assert_allclose(d.mass_matrix(sign*q),d.mass_matrix(q),rtol=2e-13,atol=2e-13)
    assert d.mass_matrix(q).shape==(4*(d.p-1),)*2


def test_formal_local_contraction_elimination_recovers_EA_without_changing_clamps():
    C,nu,gamma=Fraction(13,7),Fraction(3,10),Fraction(2,9)
    c=-nu*gamma
    energy=C*(gamma**2+2*nu*gamma*c+c**2)/2
    assert C*(c+nu*gamma)==0
    assert energy==C*(1-nu**2)*gamma**2/2
    # Bulk stationarity usually violates the resolved production c-clamp.
    assert c!=0  # this is not a uniform finite-rod limit under c(0)=c(L)=0


def _poly_product(a,b):
    out=[Fraction(0)]*(len(a)+len(b)-1)
    for i,aa in enumerate(a):
        for j,bb in enumerate(b):out[i+j]+=aa*bb
    return out


def _poly_derivative(a):
    return [i*aa for i,aa in enumerate(a) if i]


def _poly_integral_zero_one(a):
    return sum((aa/Fraction(i+1) for i,aa in enumerate(a)),Fraction(0))


def test_classical_midplane_stretching_coefficient_sign_scaling_and_variation_exact_fraction():
    # One exact fixed-fixed displacement polynomial, unrelated to solver modes.
    shape=[Fraction(0),Fraction(0),Fraction(1),Fraction(-2),Fraction(1)]
    test=[Fraction(0),Fraction(1),Fraction(-1)]
    slope=_poly_derivative(shape);test_slope=_poly_derivative(test)
    J=_poly_integral_zero_one(_poly_product(slope,slope))
    EA,L=Fraction(7,3),Fraction(1)
    gamma=J/(2*L);N=EA*gamma
    energy=EA*L*gamma**2/2
    assert energy==EA*J**2/(8*L) and energy>0
    # Quasistatic u_s=gamma-w_s^2/2 has integral zero at the fixed ends.
    assert L*gamma-J/2==0
    assert EA*L*(4*J/(2*L))**2/2==16*energy  # A doubled gives A^4
    assert _poly_integral_zero_one(_poly_product([-x for x in slope],[-x for x in slope]))==J
    deltaJ=2*_poly_integral_zero_one(_poly_product(slope,test_slope))
    variation=EA*J*deltaJ/(4*L)
    assert variation==N*_poly_integral_zero_one(_poly_product(slope,test_slope))
    curvature=_poly_derivative(slope)
    assert variation==-N*_poly_integral_zero_one(_poly_product(curvature,test))


def test_tests_never_call_integrators_or_eigen_decomposition():
    tree=ast.parse(Path(__file__).read_text(encoding="utf8"))
    forbidden={"integrate_case","solve_ivp","Radau","eigh","eig","eigvals","linear_eigenpairs"}
    for node in ast.walk(tree):
        if isinstance(node,ast.Call):
            name=node.func.id if isinstance(node.func,ast.Name) else (
                node.func.attr if isinstance(node.func,ast.Attribute) else None)
            assert name not in forbidden



def test_quartic_angular_defect_is_retained_and_through_cubic_balance_is_not_overclaimed(frozen_action):
    plane={field+suffix:0 for field in ("v","Phi","psi")
           for suffix in ("","_s","_ss","_t","_st","_tt")}
    V=frozen_action.V4.substitute(plane)
    p=frozen_action.symbols
    Fx,Fw=V.derivative("u_s"),V.derivative("w_s")
    defect=V.derivative("theta")+(1+p["u_s"])*Fw-p["w_s"]*Fx
    expected=(p["theta"]**3*p["u_s"]*(-p["C"]/2+2*p["S"]/3)
              +p["nu"]*p["C"]*p["c"]*p["w_s"]*p["theta"]**2/2
              +2*(p["C"]-p["S"])*p["u_s"]*p["w_s"]*p["theta"]**2)
    assert defect.truncate(3)==rod.Polynomial(0)
    assert defect.homogeneous(4)==expected
    assert expected!=rod.Polynomial(0)
    # Every term has weight five for this prepared ?/?? ordering.
    for monomial in expected.terms:
        degree=sum(2 if rod.SYMBOL_ORDER[index].split("_")[0] in ("u","c") else 1
                   for index in monomial if index<len(rod.FIELD_ORDER)*len(rod.JET_ORDER))
        assert degree==5
    # Generic degree-four defect is not discarded or forced to zero.


def test_cubic_material_force_rotation_requires_truncation_to_energy_order(frozen_action):
    p=frozen_action.symbols
    g1=p["u_s"]+p["w_s"]*p["theta"]-p["theta"]**2/2-p["u_s"]*p["theta"]**2/2
    g2=p["w_s"]-p["theta"]-p["u_s"]*p["theta"]-p["w_s"]*p["theta"]**2/2+p["theta"]**3/6
    N=p["C"]*(g1+p["nu"]*p["c"]);Q=p["S"]*g2
    cosine=1-p["theta"]**2/2;sine=p["theta"]-p["theta"]**3/6
    plane={field+suffix:0 for field in ("v","Phi","psi")
           for suffix in ("","_s","_ss","_t","_st","_tt")}
    Fx=frozen_action.flux_a[0].substitute(plane)
    Fw=frozen_action.flux_a[1].substitute(plane)
    assert (N*cosine-Q*sine).truncate(3)==Fx
    assert (N*sine+Q*cosine).truncate(3)==Fw
    assert N*cosine-Q*sine!=Fx  # a blind product includes discarded orders
    assert N*sine+Q*cosine!=Fw



@pytest.fixture(scope="module")
def cli():
    from scripts.analysis import check_weakly_nonlinear_planar_physics
    return check_weakly_nonlinear_planar_physics


@pytest.fixture(scope="module")
def physical_config(cli):
    return read(cli.CONFIG)


@pytest.fixture(scope="module")
def inputs(cli,physical_config):
    return cli.load_inputs(physical_config)


def test_retained_measures_are_cubic_and_match_frozen_action_not_quartic_plot_expression(cli,frozen_action):
    p=frozen_action.symbols
    g1,g2=cli.cubic_measures(p["u_s"],p["w_s"],p["theta"])
    expected1=p["u_s"]+p["theta"]*p["w_s"]-p["theta"]**2/2-p["u_s"]*p["theta"]**2/2
    expected2=p["w_s"]-p["theta"]-p["u_s"]*p["theta"]-p["w_s"]*p["theta"]**2/2+p["theta"]**3/6
    assert g1==expected1 and g2==expected2
    assert g1==g1.truncate(3) and g2==g2.truncate(3)
    assert g1.homogeneous(4)==g2.homogeneous(4)==rod.Polynomial(0)


def test_compiled_canonical_force_moment_and_scalar_flux_match_independent_action(
        cli,frozen_action,frozen_state):
    _,coeff,_=frozen_state
    rng=np.random.default_rng(47)
    fields=rng.normal(size=(5,4))*1e-3
    gradient=rng.normal(size=(5,4))*1e-3
    actual=cli.canonical_flux(frozen_action,coeff,fields,gradient)
    V=planar._restrict_to_plane(frozen_action.V4)
    p=frozen_action.symbols
    Fx,Fw=V.derivative("u_s"),V.derivative("w_s")
    expressions=(Fx,Fw,V.derivative("theta_s"),V.derivative("c_s"),
                 V.derivative("theta"),V.derivative("c"),
                 V.derivative("theta")+(1+p["u_s"])*Fw-p["w_s"]*Fx)
    expected=[]
    for row,grad in zip(fields,gradient):
        values={key:0. for key in rod.SYMBOL_ORDER}|coeff.values()
        values.update(dict(zip(("u","w","theta","c"),row)))
        values.update(dict(zip(("u_s","w_s","theta_s","c_s"),grad)))
        expected.append([expression.evaluate(values) for expression in expressions])
    np.testing.assert_allclose(actual,expected,rtol=2e-12,atol=2e-12*np.max(abs(actual)))
    assert actual.shape==(5,7)
    np.testing.assert_allclose(actual[:,2],coeff.Bp*gradient[:,2],rtol=2e-15,atol=0)
    np.testing.assert_allclose(actual[:,3],coeff.H*gradient[:,3],rtol=2e-15,atol=0)


def test_local_to_global_reaction_map_preserves_virtual_work_and_c_is_not_cartesian():
    local_force=np.array([.7,-.4])
    delta_local=np.array([-.3,.2])
    global_force=local_force*np.array([1.,-1.])
    delta_global=delta_local*np.array([1.,-1.])
    assert local_force@delta_local==global_force@delta_global
    M,Rc,delta_theta,delta_c=.2,.3,-.1,.4
    assert M*delta_theta==(-M)*(-delta_theta)
    assert Rc*delta_c==.12  # scalar work conjugate to c, not added to xy force


def test_limited_classical_report_keeps_formal_bulk_and_nonzero_rotation_defect(cli,frozen_action):
    proof=cli.limited_math(frozen_action)
    assert proof["status"]=="PASS" and all(proof["checks"].values())
    assert proof["classical"]["bulk_limit_only"] is True
    assert proof["classical"]["finite_c_clamp_elimination"] is False
    assert "8L" in proof["classical"]["Vstretch"]
    assert "K4" in proof["angular_balance"]
    K=rod.Polynomial.deserialize(proof["polynomials"]["K"])
    assert K!=rod.Polynomial(0) and K.truncate(3)==rod.Polynomial(0)


def test_loaded_leading_profiles_are_saved_exactly_and_no_bvp_or_eigenpair_is_rebuilt(inputs):
    assert inputs["summary"]["state_admitted_flag"] is False
    assert inputs["state"].admitted is False
    assert inputs["provenance"]["regenerated_profiles"] is False
    assert inputs["provenance"]["regenerated_Theta3"] is False
    np.testing.assert_array_equal(inputs["stat"]+inputs["harm"],inputs["state"].profiles.coefficients)
    with np.load(PREPARED/"coordinates/p96.npz",allow_pickle=False) as data:
        np.testing.assert_array_equal(inputs["stat"],data["legendre_stat"])
        np.testing.assert_array_equal(inputs["harm"],data["legendre_harm"])


@pytest.mark.parametrize("velocity",(False,True))
def test_asymptotic_evaluator_uses_exact_timestamps_and_fixed_amplitude_orders(cli,inputs,velocity):
    x=np.array([0.,.25,.5,.75,1.])
    omega=inputs["state"].background.omega
    times=np.array([0.,.031,.271,.903])*inputs["state"].background.T1
    epsilon=.025
    result=cli.asymptotic_fields(inputs,x,times,epsilon,velocity)
    stat=prep.LegendreProfiles(inputs["stat"]).evaluate(x)
    harm=prep.LegendreProfiles(inputs["harm"]).evaluate(x)
    pair=inputs["state"].background.evaluate(x)
    first=-omega*np.sin(omega*times) if velocity else np.cos(omega*times)
    second=-2*omega*np.sin(2*omega*times) if velocity else np.cos(2*omega*times)
    uc=second[:,None,None]*harm[None,:,:]
    if not velocity:uc=uc+stat[None,:,:]
    np.testing.assert_allclose(result[:,:,[0,3]],epsilon**2*uc,rtol=2e-15,atol=2e-15*np.max(abs(result)))
    np.testing.assert_allclose(result[:,:,1:3],epsilon*first[:,None,None]*pair[None,:,:],rtol=2e-15,atol=2e-15*np.max(abs(result)))
    # No interpolation of cached sampled analytic histories is used.


def test_initial_theta3_offset_from_linear_reference_is_known_initial_data(cli,inputs):
    x=np.linspace(0.,1.,17);epsilon=.05
    leading=cli.asymptotic_fields(inputs,x,[0.],epsilon)[0]
    prepared=inputs["state"].evaluate(x,epsilon,require_admitted=False)
    np.testing.assert_allclose(prepared[:,:2],leading[:,:2],rtol=2e-12,atol=2e-16)
    np.testing.assert_allclose(prepared[:,3],leading[:,3],rtol=2e-12,atol=2e-16)
    np.testing.assert_allclose(prepared[:,2]-leading[:,2],
                               epsilon**3*inputs["state"].correction.evaluate(x),rtol=2e-12,atol=2e-16)


def test_physical_comparison_scales_preserve_dimensions_and_velocity_factors(cli):
    e,h,L,omega=.025,.05,2.,.317
    expected_q=np.array([e**2*h,e*h,e*h/L,e**2])
    expected_v=expected_q*np.array([2*omega,omega,omega,2*omega])
    np.testing.assert_array_equal(cli.amplitude_scales(e,h,L,omega),np.r_[expected_q,expected_v])


def test_physical_sanity_configuration_is_one_half_run_with_frozen_science(physical_config):
    c=physical_config
    assert (c["p"],c["large_amplitude"],c["small_amplitude"],c["time_level"],c["periods"])==(64,.05,.025,"tight",1.)
    assert c["projection_policy"]==prep.CONSTRAINED_PROJECTION and c["projection_dps"]==70
    assert c["execution_mode"]=="EXPLORATORY_NOT_CERTIFIED"
    assert c["budget"]["max_new_integrations"]==1
    assert c["budget"]["numerical_wall_seconds"]==600
    assert c["semantics"]["physical_validation"] is False
    assert c["semantics"]["negative_amplitude_ODE"] is False
    assert c["semantics"]["classical_limit_bulk_only"] is True
    assert c["semantics"]["momentum_includes_retained_action_rotation_defect"] is True


@pytest.mark.parametrize("key,value",(("p",48),("small_amplitude",.0125),("time_level","allowed_extra"),
                                     ("periods",5.),("projection_policy",prep.UNCONSTRAINED_PROJECTION)))
def test_unrequested_case_is_rejected_before_historical_data_loading(cli,physical_config,monkeypatch,key,value):
    c=json.loads(json.dumps(physical_config));c[key]=value
    def forbidden(*args,**kwargs):
        raise AssertionError("Unrequested case must stop before source loading")
    monkeypatch.setattr(cli.saved,"validate_cache",forbidden)
    with pytest.raises(ValueError):
        cli.load_inputs(c)


@pytest.mark.parametrize("changed",("amplitude","projection","source","time"))
def test_physical_cache_identity_includes_sources_and_comparison_policy(cli,physical_config,tmp_path,changed):
    c=json.loads(json.dumps(physical_config));path=tmp_path/"config.json"
    cli.write_json(path,c);before,item=cli.identity(path)
    assert item["source_manifests"]["one_T1_source"]==sha(SOURCE/"manifest.json")
    if changed=="amplitude":c["small_amplitude"]=.03
    elif changed=="projection":c["projection_dps"]=80
    elif changed=="time":c["periods"]=.5
    else:
        folder=tmp_path/"source";folder.mkdir()
        original=read(SOURCE/"manifest.json");original["test_changed_source"]=True
        cli.write_json(folder/"manifest.json",original);c["one_T1_source"]=str(folder)
    cli.write_json(path,c);after,_=cli.identity(path)
    assert after!=before


@pytest.mark.parametrize("action",("compute","report-only","plot-only"))
def test_cached_physical_modes_evaluate_no_ode_asymptotics_reactions_or_projection(
        cli,tmp_path,monkeypatch,capsys,action):
    item={"synthetic_cache":"physical-sanity-zero-work"};bundle=tmp_path/"out"/"fixture"
    bundle.mkdir(parents=True)
    summary={"schema":"nlsp-planar-physical-sanity-v1","synthetic_fixture":True,
             "statuses":{"NLSP_PLANAR_PHYSICAL_SANITY_CHECKS":"DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS"},
             "new_ODE_integrations":1,"execution_mode":"EXPLORATORY_NOT_CERTIFIED","state_admitted_flag":False}
    cli.write_json(bundle/"summary.json",summary)
    cli.write_json(bundle/"manifest.json",cli.saved.manifest_for(bundle,item))
    def forbidden(*args,**kwargs):
        raise AssertionError("Cached physical-sanity mode may only read")
    for name in ("run_compute","load_inputs","history_analysis","reaction_checks","rhs_parity","asymptotic_fields"):
        monkeypatch.setattr(cli,name,forbidden)
    monkeypatch.setattr(prep,"stable_initial_projection",forbidden)
    plots=[]
    monkeypatch.setattr(cli,"plot_only",lambda path:plots.append(Path(path)))
    monkeypatch.setattr(cli,"identity",lambda *args:("fixture",item))
    args=["--compute","--output-dir",str(bundle.parent)] if action=="compute" else ["--"+action,str(bundle)]
    returned=cli.main(args);printed=json.loads(capsys.readouterr().out)
    assert returned==summary
    assert printed["new_ODE_BVP_eigen_symbolic_calls"]==0
    assert plots==([bundle] if action=="plot-only" else [])
    assert returned["new_ODE_integrations"]==1  # inherited saved count, not current work



@pytest.fixture(scope="module")
def actual_result(cli):
    # One authorized physical trajectory; fingerprint may change only with an
    # explicitly documented post-execution metadata/provenance revision.
    candidates=[path.parents[2] for path in cli.OUTPUT.glob("*/cases/p64_half_tight/state.npy")
                if (path.parents[2]/"manifest.json").is_file()]
    if not candidates:
        pytest.skip("Physical half-amplitude result absent; tests never reproduce it")
    bundle=max(candidates,key=lambda path:(path/"summary.json").stat().st_mtime)
    summary=cli.saved.validate_cache(bundle)
    history=cli.saved.one_T1_history(bundle,"p64_half_tight")
    return bundle,summary,history


def test_actual_physical_stage_reuses_big_run_and_keeps_certification_qualified(actual_result,inputs):
    _,summary,_=actual_result
    assert summary["new_ODE_integrations"]==1
    assert summary["new_BVP_eigensolves"]==0
    assert summary["runtime"]["new_ODE_integrations"]==1
    assert summary["runtime"]["numerical_wall_seconds"]<=600.
    assert summary["execution_mode"]=="EXPLORATORY_NOT_CERTIFIED"
    assert summary["state_admitted_flag"] is False
    assert summary["source_strict_qualification"]=="PARTIAL"
    assert summary["statuses"]["NLSP_PLANAR_PHYSICAL_SANITY_CHECKS"]=="DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS"
    assert summary["statuses"]["NLSP_SANITY_REACTIONS_AND_MOMENTUM"]=="PARTIAL"
    assert summary["statuses"]["NLSP_SANITY_AMPLITUDE_SCALING"]=="PARTIAL"
    assert "no arbitrary new physical tolerance" in summary["status_interpretation"]
    assert inputs["summary"]["source_strict_status"]=="PARTIAL"


def test_actual_half_run_has_same_T1_grid_and_recomputed_old_atol_policy(actual_result,inputs):
    _,summary,hist=actual_result
    old=inputs["summary"]["cases"]["p64_tight"];case=summary["half_case"]
    source_time=np.load(inputs["source"]/"cases/p64_tight/time.npy",mmap_mode="r")
    np.testing.assert_array_equal(hist["time"],source_time)
    assert hist["time"][-1]==case["time_end"]==case["target_time_end"]==inputs["state"].background.T1
    assert case["status"]=="PASS" and case["p"]==64 and case["time_level"]=="tight"
    assert case["amplitude_over_h"]==.025
    assert case["amplitude"]==.025*inputs["pilot"]["material_geometry"]["h"]
    assert case["rtol"]==old["rtol"] and case["max_step"]==old["max_step"]
    np.testing.assert_array_equal(case["atol"],np.asarray(old["atol"])/2)
    np.testing.assert_array_equal(case["atol_coordinate_scales"],np.asarray(old["atol_coordinate_scales"])/2)
    assert case["velocity_scale_multiplier"]==old["velocity_scale_multiplier"]
    assert hist["q"].shape==hist["velocity"].shape==(len(source_time),4*63)


def test_actual_half_initial_state_matches_saved_projection_and_is_not_half_all_coordinates(actual_result,inputs):
    bundle,summary,hist=actual_result
    with np.load(bundle/"half_initial.npz",allow_pickle=False) as initial:
        np.testing.assert_array_equal(hist["q"][0],initial["q"])
        np.testing.assert_array_equal(hist["velocity"][0],initial["velocity"])
    large=cli_read_large_initial(inputs)
    n=63
    np.testing.assert_allclose(hist["q"][0,:n],large[:n]/4,rtol=2e-12,atol=2e-12*np.max(abs(large[:n])))
    np.testing.assert_allclose(hist["q"][0,3*n:],large[3*n:]/4,rtol=2e-12,atol=2e-12*np.max(abs(large[3*n:])))
    np.testing.assert_allclose(hist["q"][0,n:2*n],large[n:2*n]/2,rtol=2e-12,atol=2e-12*np.max(abs(large[n:2*n])))
    assert np.linalg.norm(hist["q"][0]-large/2)>0.
    np.testing.assert_array_equal(hist["velocity"][0],0.)
    assert summary["initial_projection"]["pass"] is True
    provenance=read(bundle/"initial_provenance.json")
    assert provenance["Theta3_regenerated"] is False
    assert provenance["state_admitted"] is False
    assert provenance["projection_policy"]==prep.CONSTRAINED_PROJECTION


def cli_read_large_initial(inputs):
    source=np.load(inputs["source"]/"cases/p64_tight/state.npy",mmap_mode="r")
    return source[0,:4*63]


@pytest.mark.parametrize("label",("large","half"))
def test_actual_asymptotic_comparison_uses_real_timestamps_and_absolute_fixed_scales(
        actual_result,inputs,label):
    bundle,summary,hist=actual_result
    analysis=summary[label+"_analysis"]
    assert analysis["epsilon_a"]==(.05 if label=="large" else .025)
    assert analysis["actual_time_end"]==inputs["state"].background.T1
    assert len(analysis["deviations"])==8
    assert "higher orders and discretization" in analysis["qualification"]
    with np.load(bundle/(label+"_asymptotic.npz"),allow_pickle=False) as data:
        np.testing.assert_array_equal(data["time"],hist["time"])
        assert data["difference_L2"].shape==data["difference_max"].shape==(len(hist["time"]),8)
        names=[part+"_"+f for part in ("q","velocity") for f in ("u","w","theta","c")]
        for k,key in enumerate(names):
            row=analysis["deviations"][key]
            assert row["absolute_L2"]==np.max(data["difference_L2"][:,k])
            assert row["absolute_max"]==np.max(data["difference_max"][:,k])
            assert row["fixed_physical_scale"]==data["physical_scales"][k]
            assert row["fixed_scaled_L2"]==row["absolute_L2"]/row["fixed_physical_scale"]
            assert row["fixed_scaled_max"]==row["absolute_max"]/row["fixed_physical_scale"]
        assert np.all(np.isfinite(data["difference_L2"]))


def test_actual_normalized_amplitude_comparison_has_no_pointwise_zero_ratios_or_power_fit(actual_result):
    _,summary,_=actual_result
    report=summary["amplitude_comparison"]
    assert report["same_physical_time"] is True
    assert report["pointwise_ratios_near_zero"] is False
    assert report["power_law_fit"] is False
    assert len(report["rows"])==8
    for key,row in report["rows"].items():
        field=key.split("_")[-1];power=2 if field in ("u","c") else 1
        assert row["amplitude_power"]==power
        assert row["expected_characteristic_ratio"]==2.**power
        assert row["observed_max_ratio"]==row["large_max"]/row["small_max"]
        assert np.isfinite(row["observed_max_ratio"]) and row["small_max"]>0.
        # Post-execution unit conversion may divide twice; allow only
        # arithmetic representation error, not a physical acceptance tolerance.
        np.testing.assert_allclose(row["fixed_scaled_normalized_max"],
            row["normalized_difference_max"]/row["fixed_normalized_scale"],
            rtol=4*np.finfo(float).eps,atol=0.)


@pytest.mark.parametrize("label",("large","half"))
def test_actual_strain_scales_are_reported_with_locations_without_new_applicability_threshold(
        actual_result,label):
    _,summary,_=actual_result
    rows=summary[label+"_analysis"]["strains"]
    assert set(rows)=={"u_s","w_s","theta","c","L_c_s","L_theta_s","Gamma1","Gamma2","surface_bending_strain"}
    for row in rows.values():
        assert np.isfinite(row["max_abs"]) and row["max_abs"]>=0.
        assert 0.<=row["time_tau"]<=1.
        assert 0.<=row["s_over_L"]<=1.
        assert "physical_threshold" not in row
    assert rows["surface_bending_strain"]["max_abs"]==.025*rows["L_theta_s"]["max_abs"]


def test_actual_energy_ratio_is_observed_not_forced_and_safety_mass_pass(actual_result,inputs):
    _,summary,hist=actual_result
    case=summary["half_case"];old=inputs["summary"]["cases"]["p64_tight"]
    assert case["energy_and_mass_pass"] and case["safety_pass"]
    assert case["mass_lower_bound_min"]>=inputs["pilot"]["safety"]["min_relative_mass_eigenvalue"]
    assert case["relative_energy_drift_max"]<=inputs["pilot"]["gates"]["energy_relative_drift"]
    assert summary["energy_ratio_large_over_half"]==old["initial_energy"]/case["initial_energy"]
    assert "energy_matching" not in summary
    assert len(hist["time"])==case["actual_valid_rows"]


@pytest.mark.parametrize("label",("large","half"))
def test_actual_reactions_momentum_and_retained_defect_remain_distinct_quantities(actual_result,label):
    _,summary,_=actual_result
    report=summary[label+"_reactions"]
    assert report["outward_signs"]==[-1,1]
    assert "finite-p strong residual" in report["qualification"]
    assert "K4" in report["qualification"]
    for row in report["snapshots"]:
        local=np.asarray(row["endpoint_on_rod_local_Fu_Fw_M_Rc"])
        global_force=np.asarray(row["physical_global_force_xy"])
        np.testing.assert_array_equal(global_force,local[:,:2]*np.array([1.,-1.]))
        np.testing.assert_array_equal(row["physical_global_couple_z"],-local[:,2])
        np.testing.assert_array_equal(np.asarray(row["linear_momentum_derivative"])-row["endpoint_force_sum"],
                                      row["linear_balance_difference"])
        assert row["angular_momentum_derivative"]-row["endpoint_torque_about_left"]==row["angular_balance_difference"]
        assert row["weak_reaction_balance_is_bookkeeping_not_independent_validation"] is True
        assert "reaction_fit_to_balance" not in row
        assert 0.<=row["tau"]<=1.
        assert len(row["strong_residual_L2"])==4
    assert any(np.linalg.norm(row["endpoint_force_sum"])>0. for row in report["snapshots"])


def test_actual_amplitude_symmetry_is_rhs_checked_without_negative_extra_run(actual_result):
    _,summary,_=actual_result
    parity=summary["reflection"]
    assert parity["status"]=="PASS" and parity["negative_amplitude_integrations"]==0
    assert all(row["relative"]<=2e-12 for row in parity["rows"])



def test_strict_execution_label_is_rejected_before_loading_unqualified_source(cli,physical_config,monkeypatch):
    c=json.loads(json.dumps(physical_config));c["execution_mode"]="STRICT_ADMITTED"
    def forbidden(*args,**kwargs):
        raise AssertionError("Uncertified strict label must stop before source loading")
    monkeypatch.setattr(cli.saved,"validate_cache",forbidden)
    with pytest.raises(ValueError,match="Strict admission is not certified"):
        cli.load_inputs(c)


def test_audit_dependency_is_frozen_in_identity_and_rejected_on_input_hash_change(
        cli,physical_config,monkeypatch):
    _,identity=cli.identity(cli.CONFIG)
    assert identity["audit_result_sha256"]==sha(AUDIT/"result.json")
    original=cli.sha
    target=AUDIT/"result.json"
    monkeypatch.setattr(cli,"sha",lambda path:
                        "0"*64 if Path(path).resolve()==target.resolve() else original(path))
    with pytest.raises(ValueError,match="Frozen action archive hash mismatch"):
        cli.load_inputs(physical_config)


@pytest.mark.parametrize("label",("large","half"))
def test_saved_asymptotic_spatial_snapshots_keep_actual_time_and_difference_orders(actual_result,inputs,cli,label):
    bundle,_,hist=actual_result
    epsilon=.05 if label=="large" else .025
    with np.load(bundle/(label+"_comparison_snapshots.npz"),allow_pickle=False) as data:
        assert data["fields"].shape==data["leading_fields"].shape==data["difference"].shape==(6,501,4)
        assert np.all(np.isin(data["time"],hist["time"]))
        np.testing.assert_array_equal(data["difference"],data["fields"]-data["leading_fields"])
        np.testing.assert_array_equal(data["velocity_difference"],data["velocities"]-data["leading_velocities"])
        expected=cli.asymptotic_fields(inputs,data["s"],data["time"],epsilon)
        np.testing.assert_array_equal(data["leading_fields"],expected)
        source=(inputs["source"]/"cases/p64_tight/snapshots.npz" if label=="large"
                else bundle/"cases/p64_half_tight/snapshots.npz")
        with np.load(source,allow_pickle=False) as original:
            np.testing.assert_array_equal(data["fields"],original["fields"])
            np.testing.assert_array_equal(data["time"],original["time"])
        assert data["Gamma1_cubic"].shape==data["Gamma2_cubic"].shape==(6,501)


def test_velocity_normalized_scales_have_frequency_units_and_source_uncertainty_is_not_removed(actual_result,inputs):
    _,summary,_=actual_result
    h,L,omega=inputs["state"].background.h0,inputs["state"].length,inputs["state"].background.omega
    q_scales=np.array([h,h,h/L,1.])
    v_scales=q_scales*np.array([2*omega,omega,omega,2*omega])
    for k,field in enumerate(("u","w","theta","c")):
        assert summary["amplitude_comparison"]["rows"]["q_"+field]["fixed_normalized_scale"]==q_scales[k]
        assert summary["amplitude_comparison"]["rows"]["velocity_"+field]["fixed_normalized_scale"]==v_scales[k]
    assert summary["small_amplitude_independent_p_time_check"] is False
    for key,value in summary["source_numerical_uncertainty"].items():
        assert value==inputs["summary"][key]
    assert summary["post_execution_revision"]["ODE_reexecuted"] is False


@pytest.mark.parametrize("action",("compute","report-only"))
def test_actual_physical_cache_reuses_final_identity_with_zero_new_work(actual_result,cli,monkeypatch,capsys,action):
    bundle,summary,_=actual_result
    def forbidden(*args,**kwargs):
        raise AssertionError("Final physical cache cannot rerun any stage")
    for name in ("run_compute","load_inputs","history_analysis","reaction_checks","rhs_parity","asymptotic_fields"):
        monkeypatch.setattr(cli,name,forbidden)
    item=read(bundle/"manifest.json")["identity"]
    monkeypatch.setattr(cli,"identity",lambda *args:(bundle.name,item))
    args=["--compute","--output-dir",str(bundle.parent)] if action=="compute" else ["--report-only",str(bundle)]
    returned=cli.main(args);printed=json.loads(capsys.readouterr().out)
    assert returned==summary and printed["new_ODE_BVP_eigen_symbolic_calls"]==0

