"""Targeted checks without any full-model spectral computation."""
from fractions import Fraction as F
import json
import numpy as np
import pytest
from scripts.lib.isotropic_two_springs_external_transverse_spring import Leading
from scripts.lib import two_springs_transverse_leading_diagnostic as d

RESULT=d.BASE/'bcc1660d87bd57a3'

def test_exact_curvature_gram_and_same_kinematic_space():
    old=Leading(.1,10.,order=24,k_w=1.)
    fixed=Leading(.1,10.,order=24,k_w=1.,analytic_curvature=True)
    for key in ('M','Kj','Ke','Kw','P','constraint','Z'):
        assert np.array_equal(getattr(old,key),getattr(fixed,key))
    assert np.array_equal(fixed.Kb,np.diag(np.diag(fixed.Kb)))
    for arm,L in enumerate(fixed.lengths):
        rows=d.basis(24,arm)
        for i,ra in enumerate(rows):
            for j,rb in enumerate(rows):
                exact=d.dot_integral(d.deriv(d.deriv(ra)),d.deriv(d.deriv(rb)))
                assert fixed.Kb[arm*24+i,arm*24+j]==pytest.approx(float(exact)/L**3,rel=3e-16,abs=0)
    assert fixed.outer[1][0][-1]==1
    assert fixed.Z.shape==(49,48)
    assert np.linalg.eigvalsh(fixed.Z.T@fixed.M@fixed.Z).min()>0

def test_polynomials_and_independent_derivatives():
    result=d.polynomial_checks()
    assert result['monomial_relative_max']<1e-50
    assert result['existing_fields_relative_max']<1e-10
    assert result['constant_derivatives']==[1,0,0,0,0]

def test_green_identity_and_endpoint_signs_exact():
    # W=t^5+2t^3+3t^2, Z=1-3t^2+2t^3; only Z(0)=1 is nonzero.
    W=[F(0),F(0),F(3),F(2),F(0),F(1)]
    Z=[F(1),F(0),F(-3),F(2)]
    def derivative(a):return [i*a[i] for i in range(1,len(a))]
    def product_integral(a,b):return sum(x*y/F(i+j+1) for i,x in enumerate(a) for j,y in enumerate(b))
    W2=derivative(derivative(W));W3=derivative(W2);W4=derivative(W3)
    Z2=derivative(derivative(Z))
    assert product_integral(W2,Z2)-product_integral(W4,Z)==W3[0]
    # Variation with only Z'(0)=1: Z=t-2t^2+t^3 gives -W''(0).
    Z=[F(0),F(1),F(-2),F(1)]
    assert product_integral(W2,derivative(derivative(Z)))-product_integral(W4,Z)==-W2[0]

def test_saved_failure_and_default_regression():
    original=json.loads((RESULT/'original_failed_point.json').read_text())
    reference=json.loads((RESULT/'derivative_reference.json').read_text())
    regression=json.loads((RESULT/'default_regression.json').read_text())
    assert original['solution']['diagnostics']['strong_ode']>1e-9
    assert reference[-1]['strong']>1e-9
    assert reference[-1]['independent_derivatives_max_abs']<1e-65
    assert regression['zeta_max_relative']==0 and regression['coefficient_max_abs']==0

def test_saved_corrected_physical_gates_and_control():
    decision=json.loads((RESULT/'decision.json').read_text())
    assert decision['status']=='DIAGNOSTIC_PASS' and decision['tolerance']==1e-9
    for p in (16,20,22,24,26,28):
        data=json.loads((RESULT/f'corrected_order_{p}.json').read_text())
        assert data['diagnostics']['strong_ode']<1e-9
        assert data['diagnostics']['natural_boundary']<1e-9
        assert data['diagnostics']['mass_positive']
    for p in (20,24):
        data=json.loads((RESULT/f'control_order_{p}.json').read_text())
        assert data['xi']==10 and data['diagnostics']['strong_ode']<1e-9
    assert decision['control_frequency_order_change']<1e-8
    assert json.loads((RESULT/'execution.json').read_text())['new_full_model_solves']==0

def test_original_archive_is_immutable():
    source=d.ROOT/json.loads(d.CONFIG.read_text())['source']
    assert d.verify_source(source)['artifacts_checked']==977
