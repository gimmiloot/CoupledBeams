"""HS-LA2 regression of accepted finite-space identities; no spectrum solves."""
import csv
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
import numpy as np
import pytest
from scripts.lib import two_springs_horizontal_spectral_structure as audit
from scripts.lib.isotropic_two_springs_horizontal_leading import Leading

CFG=json.loads(audit.CONFIG.read_text())
FOLDER=audit.BASE/audit.identity(CFG)[0]


def rows(name):
    with open(FOLDER/name,encoding='utf8') as f:return list(csv.DictReader(f))


def point(p,tau):
    name=hashlib.sha256(repr((p,tau)).encode()).hexdigest()[:16]
    r=json.loads((FOLDER/'points'/f'{name}.json').read_text())
    with np.load(FOLDER/'points'/f'{name}.npz') as a:c=a['coefficients'].copy()
    return r,c


def supplemental_metrics():
    records=json.loads((FOLDER/'point_records.json').read_text())
    models={p:Leading(order=p) for p in CFG['orders']}
    shape_error=0.;minimum_increment=math.inf
    for tau in {r['tau'] for r in records}:
        a,ca=point(20,tau);b,cb=point(24,tau)
        # Compare fields on a common quadrature, not coefficient vectors with
        # different dimensions. A sign choice only; no branch tracking.
        va=[];vb=[];weights=[]
        for arm,L in enumerate(models[24].lengths):
            va.append(models[20].field_matrices(arm,models[24].nodes)[0]@ca)
            vb.append(models[24].values[arm][0]@cb)
            weights.extend(models[24].weights*L)
        va=np.vstack(va);vb=np.vstack(vb);w=np.array(weights)
        signs=np.sign(np.sum(w[:,None]*va*vb,axis=0))
        shape_error=max(shape_error,float(max(np.sqrt(np.sum(w[:,None]*(va-vb*signs)**2,axis=0)))))
    for p in models:
        r=sorted((r for r in records if r['order']==p),key=lambda r:math.inf if r['tau']=='RIGID' else r['tau'])
        z=np.array([a['zeta'] for a in r])
        minimum_increment=min(minimum_increment,float(np.min(np.diff(z,axis=0))))
    return dict(maximum_mass_L2_shape_difference=shape_error,minimum_adjacent_zeta_increment=minimum_increment,
        all_mass_positive=all(r['diagnostics']['mass_positive'] for r in records),
        all_stiffness_positive=all(r['diagnostics']['stiffness_positive'] for r in records),
        full_model_solves=0,additional_eigenvalue_solves=0)


def test_archive_provenance_and_unchanged_hs_l1():
    m=json.loads((FOLDER/'manifest.json').read_text())
    assert m['status']=='PASS' and m['full_model_solves']==0
    assert m['new_low_prefixes']==80 and m['reused_low_prefixes']==10
    assert m['full_finite_basis_solves']==8
    audit.Audit(CFG).verify()


@pytest.mark.parametrize('p',[20,24])
def test_all_N_pairs_in_authoritative_arithmetic(p):
    key=hashlib.sha256(repr((p,0.)).encode()).hexdigest()[:16]
    with np.load(FOLDER/f'full_basis_{key}.npz') as b:
        K=b['K'].copy();M=b['M'].copy();d=b['d'].copy()
    r=json.loads((FOLDER/f'full_basis_{key}_80digits.json').read_text())
    with mp.workdps(80):
        c=mp.matrix(r['coefficients']);z=[mp.mpf(a) for a in r['zeta']]
        km=mp.matrix(K.tolist());mm=mp.matrix(M.tolist())
        assert len(z)==2*p
        assert max(abs(v) for v in c.T*mm*c-mp.eye(len(z)))<mp.mpf('1e-60')
        for j in range(len(z)):
            lhs=km*c[:,j];rhs=z[j]*mm*c[:,j]
            assert mp.norm(lhs-rhs)/(mp.norm(lhs)+mp.norm(rhs))<mp.mpf('1e-60')
        projections=mp.matrix(d.tolist()).T*c
        assert max(abs(projections[j]-mp.mpf(r['d_n'][j])) for j in range(len(z)))<mp.mpf('1e-60')


def test_rank_one_determinant_identity_off_poles():
    m=Leading(order=20);Z=m.Z;K=Z.T@m.K0@Z;M=Z.T@m.M@Z;d=Z.T@m.D
    assert np.linalg.matrix_rank(Z.T@m.P@Z)==1
    assert np.max(abs(Z.T@m.P@Z-np.outer(d,d)))<1e-14
    for z,tau in [(-1.,1.),(20.,24.),(300.,6.135887525019995)]:
        A=K-z*M;B=A+tau*np.outer(d,d)
        sa,la=np.linalg.slogdet(A);sb,lb=np.linalg.slogdet(B)
        ratio=sb/sa*math.exp(lb-la)
        assert ratio==pytest.approx(1+tau*(d@np.linalg.solve(A,d)),rel=1e-11,abs=1e-11)
    for r in rows('secular.csv'):
        assert abs(float(r['equation_all_modes']))<1e-8
        assert float(r['frequency_relative_error'])<1e-7
        assert r['status']=='ACTIVE_REGULAR'


def test_exact_constraint_reaction_and_natural_forces():
    m=Leading(order=24);r,c=point(24,'RIGID')
    cinf,force=audit.reaction(m,dict(zeta=r['zeta'],coefficients=c))
    assert force['r']==pytest.approx(1.284214328402106,rel=1e-11)
    assert force['kkt_residual']<1e-9
    assert abs(m.D@cinf)<1e-9
    assert abs(m.outer[1][0]@cinf)>.01  # D=0 does not mean W2(0)=0.
    assert force['r_outer']==pytest.approx(force['r'],rel=1e-8)
    assert force['r_joint']==pytest.approx(force['r'],rel=1e-8)
    for r in json.loads((FOLDER/'reaction.json').read_text()):
        assert r['r_squared']*r['Gprime_all_modes']==pytest.approx(1,rel=1e-8)


def test_energy_sensitivity_two_steps_and_scaling():
    for r in rows('sensitivity.csv'):
        z,t,D=float(r['zeta']),float(r['tau']),float(r['D_abs'])
        assert float(r['P_D'])==pytest.approx(t*D*D/z,abs=1e-14)
        assert 0<=float(r['P_D'])<=1
    for r in rows('finite_differences.csv'):assert float(r['relative_error'])<5e-5
    for r in rows('curvature.csv'):
        assert r['status']=='PASS' and float(r['curvature_all_modes'])<0
    for tau in [2.9124184246057174,6.135887525019995]:
        beta5=math.sqrt(tau*(1/20000+2*CFG['epsilon']**2))
        beta6=math.sqrt(tau*2*CFG['epsilon']**2)
        assert beta5/beta6==pytest.approx(math.sqrt(2),rel=1e-14)


def test_saturation_error_is_reported_not_fitted_away():
    r=[r for r in rows('threshold_asymptotics.csv') if r['order']=='24']
    assert float(r[0]['zeta_prediction_relative_error'])==pytest.approx(.15528303243012,rel=1e-10)
    assert float(r[1]['omega_prediction_relative_error'])==pytest.approx(.04188161754456,rel=1e-10)
    tail=[r for r in rows('saturation.csv') if r['order']=='24']
    assert np.all(np.diff([abs(float(r['zeta_limit_relative_error'])) for r in tail])<0)


def test_shapes_positive_forms_and_neighbour_audit():
    s=supplemental_metrics()
    assert s['maximum_mass_L2_shape_difference']<CFG['policy']['shape_convergence_tol']
    assert s['minimum_adjacent_zeta_increment']>0
    assert s['all_mass_positive'] and s['all_stiffness_positive']


def test_plot_only_does_not_solve(tmp_path,monkeypatch):
    for name in ['saturation.csv','sensitivity.csv']:(tmp_path/name).write_bytes((FOLDER/name).read_bytes())
    def forbidden(*a,**k):raise AssertionError('plot_only must not calculate spectra')
    monkeypatch.setattr(Leading,'__post_init__',forbidden)
    monkeypatch.setattr(audit,'all_pairs',forbidden)
    audit.plot_only(tmp_path)
    from PIL import Image
    with Image.open(next((tmp_path/'figures').glob('*.png'))) as im:assert min(im.info['dpi'])>=300


# HS-LA3 adds postprocessing checks; the original HS-LA2 tests are unchanged.
from scripts.lib import two_springs_horizontal_saturation as second


def second_archive():
    cfg=json.loads(second.CONFIG.read_text())
    return second.BASE/second.identity(cfg)[0],cfg


def test_second_coefficient_indefinite_linear_benchmark():
    # Exact two-dimensional algebra, unrelated to a new beam parameter case.
    # G=1/(1-z)+4/(4-z), zero z=8/5, a=36/125, b=-324/3125.
    result,_,_=second.coefficient_from_linear(np.diag([1.,4.]),np.eye(2),np.array([1.,2.]),8/5)
    assert result['a']==pytest.approx(36/125,rel=1e-14)
    assert result['b']==pytest.approx(-324/3125,rel=1e-14)
    assert max(result['residual_u'],result['residual_v'])<1e-14


@pytest.mark.parametrize('p',[20,24])
def test_second_coefficient_authoritative_sum_and_saved_linear_vectors(p):
    folder,cfg=second_archive();data=json.loads((folder/'coefficients.json').read_text())
    result=next(r for r in data if r['order']==p)
    key=hashlib.sha256(repr((p,0.)).encode()).hexdigest()[:16]
    authoritative=json.loads((FOLDER/f'full_basis_{key}_80digits.json').read_text())
    modal=second.coefficient_from_modes(authoritative,result['zeta_infty'])
    with np.load(FOLDER/f'full_basis_{key}.npz') as source:
        K=source['K'];M=source['M'];d=source['d'];A=K-result['zeta_infty']*M
    with np.load(folder/'linear_vectors.npz') as arrays:u=arrays[f'p{p}_u'];v=arrays[f'p{p}_v']
    assert np.linalg.norm(A@u-d)/(np.linalg.norm(A@u)+np.linalg.norm(d))<1e-9
    assert np.linalg.norm(A@v-M@u)/(np.linalg.norm(A@v)+np.linalg.norm(M@u))<1e-9
    b=float(u@M@v/(u@M@u)**3)
    assert b==pytest.approx(modal['b'],rel=1e-10)
    assert b==pytest.approx(-.639617945455777,abs=1e-14)
    assert -.643974<b<=-.638619
    assert result['analytic_lower']<b<=result['analytic_upper']
    assert modal['N']==2*p and modal['negative_directions']==1


def test_second_coefficient_large_tau_inverse_and_angle_scaling():
    folder,cfg=second_archive()
    with open(folder/'threshold_predictions.csv') as f:r=list(csv.DictReader(f))
    coefficients=json.loads((folder/'coefficients.json').read_text())
    for item in r:
        c=next(c for c in coefficients if c['order']==int(item['order']))
        a=c['a_reaction'];b=c['linear']['b'];t=float(item['tau_two']);delta=float(item['delta'])
        assert a/t+b/t**2==pytest.approx(delta,rel=1e-14)
        assert t>float(item['tau_other']) and t>-2*b/a
        assert abs(float(item['tau_two_relative_error']))<.002
        assert abs(float(item['beta_two_relative_error']))<.001
        assert abs(float(item['tau_two_relative_error']))<abs(float(item['tau_one_relative_error']))
    for p in [20,24]:
        for q in [.9,.95]:
            a,b=[next(v for v in r if int(v['order'])==p and float(v['q'])==q and v['id']==id) for id in ['H5','H6']]
            assert float(a['beta_two_deg'])/float(b['beta_two_deg'])==pytest.approx(math.sqrt(2),rel=1e-14)


def test_second_coefficient_comparison_and_zero_new_spectra():
    folder,cfg=second_archive();v=json.loads((folder/'validation.json').read_text())
    assert v['new_full_model_eigensolves']==v['new_leading_eigensolves']==0
    assert not v['fit'] and v['new_linear_system_solves']==4
    assert v['coefficient_order_relative_error']<cfg['tolerances']['order_convergence_relative']
    with open(folder/'asymptotic_comparison.csv') as f:r=list(csv.DictReader(f))
    assert len(r)==8
    for item in r:
        assert abs(float(item['delta_two_relative_error']))<abs(float(item['delta_one_relative_error']))
        assert abs(float(item['omega_two_relative_error']))<abs(float(item['omega_one_relative_error']))


def test_second_coefficient_cached_run_cannot_solve(monkeypatch):
    def forbidden(*a,**k):raise AssertionError('cached HS-LA3 must only verify')
    monkeypatch.setattr(second,'coefficient_from_linear',forbidden)
    monkeypatch.setattr(second,'coefficient_from_modes',forbidden)
    monkeypatch.setattr(Leading,'solve',forbidden)
    monkeypatch.setattr(audit,'all_pairs',forbidden)
    assert second.run()==second_archive()[0]
