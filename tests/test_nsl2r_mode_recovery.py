"""Saved-root NSL-2R regressions; no structural eigenvalue/root searches."""
import numpy as np
import pytest
from scripts.analysis.nonsymmetric_laminated_beams import diagnose_mode_recovery as diagnostic
from scripts.lib import nonsymmetric_laminated_beam as single


@pytest.fixture(scope="module")
def saved():
    return diagnostic.source(diagnostic.read(diagnostic.CONFIG))


@pytest.mark.parametrize("kind,position",[("COUPLED",7),("SINGLE",6),("SINGLE",7)])
def test_problem_form_passes_original_physical_tolerance_without_root_change(saved,kind,position):
    config,data,props,models=saved
    row=data["spectral"]["STRAIGHT_RIGID"]["modes"][position-1] if kind=="COUPLED" else data["straight_single"]["count"]["roots"][position-1]
    omega=row["omega"];model=models[kind]
    c,d=diagnostic.recover_mode(model,omega,config["policy"])
    assert all(diagnostic.physical_gates(d,config["policy"]).values())
    assert d["boundary_scaled"]<1e-9
    # All physical checks are made by the untouched original evaluator.
    independent=model.evaluate(omega,c,config["policy"]["mode_quadrature"])
    assert independent["boundary_scaled"]==d["boundary_scaled"]
    assert abs(independent["mass"]-1)<1e-12
    assert row["omega"]==omega


def test_column_change_preserves_boundary_equation(saved):
    config,data,props,models=saved
    model=models["SINGLE"];omega=data["straight_single"]["count"]["roots"][6]["omega"]
    b=model.boundary(omega);scales=np.linalg.norm(b,axis=0)
    v=np.array([.1,-.2,.3,.4,-.5,.6])
    assert np.allclose((b/scales)@v,b@(v/scales),rtol=1e-14,atol=1e-14)
    assert abs(np.linalg.det(np.diag(1/scales)))>0
    diag,_,_=diagnostic.null_diagnostics(diagnostic.raw_boundary(model,omega))
    assert diag["column_balanced_nonzero_condition"]<diag["nonzero_condition"]/100


def test_bounded_comparison_basis_satisfies_same_ode(saved):
    config,data,props,models=saved
    omega=data["straight_single"]["count"]["roots"][6]["omega"]
    beam=diagnostic.BoundedBeam(props["HHLL"],2.)
    y=beam.basis(omega,[0.,.37,1.,2.])
    derivative=beam.basis(omega,[0.,.37,1.,2.],True)
    rhs=single.control_state_matrix(1j*omega,beam.a).real@y
    error=np.max(abs(derivative-rhs))/np.max(abs(rhs))
    assert error<1e-12
