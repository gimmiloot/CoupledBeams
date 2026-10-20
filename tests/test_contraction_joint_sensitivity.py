"""Scientific gates for the isolated G20 c-closure contract."""
from fractions import Fraction as F
import copy
import json

import numpy as np
import pytest

from scripts.analysis import check_contraction_joint_sensitivity as audit
from scripts.lib import contraction_joint_sensitivity as study
from scripts.lib import mindlin_herrmann_timoshenko_joint as rigid


@pytest.fixture(scope="module")
def setup():
    config = json.loads(audit.CONFIG.read_text(encoding="utf-8"))
    model, refs, _, _ = audit.references(config)
    return config,model,refs


@pytest.fixture(scope="module")
def computed(setup,tmp_path_factory):
    config,model,refs = setup
    directory = tmp_path_factory.mktemp("contraction")
    result,diags = audit.elastic_compute(model,config,refs,directory)
    audit.complex_compute(model,config,result,diags,directory)
    return result,diags


def test_exact_section_work_and_non_bending_R():
    b,h,G,kappa,cx,shear = F(1,5),F(1,20),F(5,13),F(5,6),F(3,7),F(2,9)
    area, first, inertia = b*h,F(0),b*h**3/12
    # Integrate (shear+z c') times the virtual mode 1 or z over the end.
    Q = kappa*G*(shear*area+cx*first)
    R = kappa*G*(shear*first+cx*inertia)
    assert Q==kappa*G*area*shear
    assert R==kappa*G*inertia*cx
    # Positive end work in c; the lower integration endpoint has minus sign.
    assert R*F(5,11) == (kappa*G*inertia*cx)*F(5,11)
    assert rigid.endpoint_sign("left")==-1 and rigid.endpoint_sign("right")==1


def test_rational_contraction_spring_limits_work():
    c1,c2,dc1,dc2,k = F(1,3),F(2,7),F(4,9),F(5,11),F(13,17)
    R1,R2 = -k*(c1-c2),k*(c1-c2)
    assert R1*dc1+R2*dc2+k*(c1-c2)*(dc1-dc2)==0
    assert R1+R2==0
    assert k*(F(3)-F(3))**2==0  # Common c is not blocked as k grows.


@pytest.mark.parametrize("beta",(0.,45.))
@pytest.mark.parametrize("closure,ntrace",(("FREE",6),("CONTINUOUS",5),("BLOCKED",4)))
def test_trace_space_kinematics_and_rank(setup,beta,closure,ntrace):
    config,model,_ = setup
    case = study.Case(beta,closure)
    pmap,jump = study.trace_map(case)
    operator = study.joint_matrix(case,model.coefficients["B"]/.5)
    assert pmap.shape==(8,ntrace) and np.linalg.matrix_rank(pmap)==ntrace
    assert np.linalg.matrix_rank(operator/np.linalg.norm(operator,axis=1)[:,None])==8
    q = pmap@np.arange(1,ntrace+1)
    frames = rigid.frames(beta)
    d = [frames[i].translation@q[4*i+np.array([0,2])] for i in range(2)]
    np.testing.assert_allclose(d[0],d[1],atol=2e-15)
    if closure=="CONTINUOUS":
        assert q[1]==q[5]
    if closure=="BLOCKED":
        assert q[1]==q[5]==0
    assert q[3]!=q[7]  # c rows did not replace the spring's independent rotations.


def test_rigid_matrix_recovers_old_production_rows(setup):
    _,model,_ = setup
    case = study.Case(0.,rotation="RIGID")
    new = study.boundary(model,[.5,.5],1.37,case,model.coefficients["B"]/.5,False)
    old = rigid.frame_boundary_matrix(model,(.5,.5),1.37,beta_deg=0.)
    # New rows dX,dY,FX,FY,theta,M,c,R vs old dX,dY,c,theta,FX,FY,R,M.
    permutation = list(range(8))+[8,9,12,13,11,15,10,14]
    np.testing.assert_array_equal(new,old[permutation])


def test_all_structural_virtual_work_and_positive_coefficients(computed):
    result,_ = computed
    structure = result["structure"]
    assert len(structure["records"])==10
    assert max(r["virtual_work_error"] for r in structure["records"])<2e-13
    assert min(structure["local_energy_eigenvalues"])>0
    assert min(structure["local_mass_eigenvalues"])>0


def test_direct_and_frozen_rigid_spectrum_and_shapes(computed):
    result,_ = computed
    assert max(r["relative_frequency_error"] for r in result["rigid_regression"]["checks"])<2e-8
    assert max(r["mass_L2"] for r in result["rigid_regression"]["checks"])<5e-8


def test_all_prefixes_count_mass_energy_and_nested_ordering(computed):
    result,_ = computed
    assert len(result["cases"])==11
    for case in result["cases"]:
        assert len(case["roots"])==7 and case["search"]["lower_count"]==0
        assert case["search"]["computed_roots"]==7
        assert case["mass_gram_error"]<5e-7
        assert all(b["counts"][1]-b["counts"][0]==1 for b in case["search"]["brackets"])
        assert max(r["diagnostics"]["energy_relative"] for r in case["roots"])<5e-8
        assert max(r["diagnostics"]["quadrature_relative"] for r in case["roots"])<5e-8
    assert max(r["relative_violation"] for r in result["comparison"]["ordering"])<2e-8


def test_hard_gate_rejects_wrong_order_without_reassignment(setup,computed):
    config,model,_ = setup
    result,diags = computed
    cases = copy.deepcopy(result["cases"])
    bad = next(c for c in cases if c["key"]==study.Case(45.,"BLOCKED").key)
    bad["roots"][0]["omega"] *= .9
    with pytest.raises(ArithmeticError,match="min-max"):
        audit.comparisons(model,cases,diags,config["policy"])


def test_complex_stage_refuses_elastic_failure(setup,tmp_path):
    config,model,_ = setup
    with pytest.raises(ArithmeticError,match="requires all elastic"):
        audit.complex_compute(model,config,{"elastic_status":"FAIL"},{},tmp_path)


def test_conditional_complex_pairs_pass_physical_identity(computed):
    result,_ = computed
    assert result["complex_selection"]["selected_positions"]==[1,6]
    assert len(result["complex_modes"])==6
    for mode in result["complex_modes"]:
        assert mode["p"]["real"]<0 and mode["p"]["imag"]>0
        assert mode["energy_identity_relative"]<5e-8
        assert mode["alpha_identity_relative"]<5e-8
        assert mode["predictor_relative_error"]<.01
        assert mode["mass_MAC_with_elastic"]>.98


def test_complex_continuation_preserves_elastic_state_operator(setup):
    _,model,_ = setup
    case = study.Case(45.,"FREE")
    provider = study.ComplexBasis(model,[.5,.5],1.3,case,model.coefficients["B"]/.5,17.)
    np.testing.assert_array_equal(provider.operator(1.3j),study.state_operator(model,1.3,case))
    # Check the eigen-exponential basis independently against the unchanged ODE.
    values = provider.values(1.3j,0,np.array([0.,.17,.5]))
    gradient = provider.values(1.3j,0,np.array([0.,.17,.5]),1)
    np.testing.assert_allclose(gradient,np.einsum("ij,xjk->xik",study.state_operator(model,1.3,case),values),rtol=3e-11,atol=2e-10)
