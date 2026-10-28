"""HS-LA4: archived physical states, independent counts and exact constraints."""
import csv
import hashlib
import json
import math
import numpy as np
import pytest
from scipy.linalg import null_space,solve
from scripts.lib import two_springs_horizontal_joint_stiffness as study

CFG=json.loads(study.CONFIG.read_text());FOLDER=study.BASE/study.identity(CFG)[0]


def rows(name):
    with open(FOLDER/name,encoding='utf8') as f:return list(csv.DictReader(f))


def supplemental_checks():
    records=json.loads((FOLDER/'point_records.json').read_text())
    coefficients=json.loads((FOLDER/'coefficients.json').read_text())
    finite=[r for r in records if r['order']==24 and r['role']=='primary' and r['k_joint']!='RIGID' and r['tau']!='RIGID']
    concavity=0.
    for x in finite:
        for y in finite:
            # The tangent plane of a simple first eigenvalue is a global
            # upper bound by the Rayleigh principle. No new solve is needed.
            mode=x['per_mode'][0]
            plane=x['zeta'][0]+mode['J_squared']*(y['k_joint']-x['k_joint'])+mode['D_squared']*(y['tau']-x['tau'])
            concavity=max(concavity,(y['zeta'][0]-plane)/max(1,abs(plane)))
    primary=[r for r in records if r['role']=='primary']
    minimum_increment=math.inf;frequency_convergence=0.
    for p in CFG['orders']:
        for t in [0,'RIGID']:
            spectra=[next(r['zeta'] for r in primary if r['order']==p and r['k_joint']==k and r['tau']==t) for k in CFG['k_joint']]
            minimum_increment=min(minimum_increment,float(np.diff(spectra,axis=0).min()))
        for k in CFG['k_joint']:
            states=sorted((r for r in primary if r['order']==p and r['k_joint']==k),key=lambda r: math.inf if r['tau']=='RIGID' else r['tau'])
            minimum_increment=min(minimum_increment,float(np.diff([r['zeta'] for r in states],axis=0).min()))
    for k in CFG['k_joint']:
        states=[]
        for p in CFG['orders']:
            states.append(sorted((r for r in primary if r['order']==p and r['k_joint']==k),key=lambda r: math.inf if r['tau']=='RIGID' else r['tau']))
        for x,y in zip(*states):
            a=np.array(x['zeta']);b=np.array(y['zeta']);positive=a>0
            frequency_convergence=max(frequency_convergence,float(np.max(abs(np.sqrt(a[positive]/b[positive])-1))))
    coefficients_d=[];derivatives=[];secular=[]
    for c in coefficients:
        if c['order']!=24 or c['k_joint']=='RIGID':continue
        k=c['k_joint'];m=study.JointLeading(k);key=hashlib.sha256(repr((k,24)).encode()).hexdigest()[:16]
        with np.load(FOLDER/f'linear_{key}.npz') as data:
            u=data['u'];v=data['v'];K=data['K0'];M=data['M'];d=data['d']
        inf=next(r for r in records if r['order']==24 and r['k_joint']==k and r['tau']=='RIGID')
        j=m.Z.T@m.jump;gp=c['linear']['Gprime'];gpp=c['linear']['Gsecond'];mu_prime=inf['per_mode'][0]['J_squared']
        g_zeta_k=-2*float(j@u)*float(j@v)
        aprime=-(g_zeta_k+gpp*mu_prime)/gp**2
        coefficients_d.append(dict(k_joint=k,a_prime=aprime,mu_prime=mu_prime,
            mu_prime_resolvent=float((j@u)**2/gp),mu_prime_identity_relative_error=abs(float((j@u)**2/gp)/mu_prime-1)))
        for q in CFG['levels']:
            row=next(r for r in rows('sensitivities.csv') if r['order']=='24' and float(r['k_joint'])==k and float(r['q'])==q)
            z=float(row['target_zeta']);tau=float(row['tau']);B=m.Z.T@m.Kb@m.Z-z*M
            X=solve(B,np.column_stack((j,d)),assume_a='sym');gJJ=float(j@X[:,0]);gDD=float(d@X[:,1]);gJD=float(j@X[:,1])
            determinant=(1+k*gJJ)*(1+tau*gDD)-k*tau*gJD**2
            secular.append(dict(k_joint=k,q=q,rank_two_secular_residual=determinant))
        if k not in CFG['finite_difference']['joints']:continue
        for h in (CFG['finite_difference']['zero_steps'] if k==0 else [k*f for f in CFG['finite_difference']['relative_steps']]):
            def archived_a(kk):
                r=next(r for r in records if r['order']==24 and r['k_joint']==kk and r['tau']=='RIGID')
                return r['per_mode'][0]['r_D']**2
            fd=(-3*c['a']+4*archived_a(h)-archived_a(2*h))/(2*h) if k==0 else (archived_a(k+h)-archived_a(k-h))/(2*h)
            derivatives.append(dict(k_joint=k,step=h,a_prime=aprime,difference=fd,relative_error=abs(fd/aprime-1)))
    return dict(status='PASS',new_eigensolves=0,concavity_tangent_max_violation=concavity,
        neighbour_minimum_zeta_increment=minimum_increment,prefix_frequency_convergence=frequency_convergence,
        a_derivatives=coefficients_d,a_difference_checks=derivatives,rank_two_checks=secular)


def test_exact_mechanism_mass_and_one_dimensional_kernel():
    m=study.JointLeading(0);v=np.zeros(len(m.D));v[m.order]=1.01;v[-1]=-1.01
    assert np.linalg.norm(m.K0@v)==0
    assert abs(m.constraint@v).max()<1e-14
    assert float(v@m.M@v)==pytest.approx(1.01**3/3,rel=1e-13)
    assert float(m.jump@v)==pytest.approx(-1)
    assert float(m.D@v)==pytest.approx(-1.01)
    assert len(m.Z.T@m.Kb@m.Z)-np.linalg.matrix_rank(m.Z.T@m.Kb@m.Z)==1
    r=next(r for r in json.loads((FOLDER/'point_records.json').read_text()) if r['k_joint']==0 and r['order']==24 and r['tau']==0)
    assert r['zeta'][0]==0 and r['diagnostics']['nullity']==1
    assert r['per_mode'][0]['strong']==r['per_mode'][0]['energy']==0


def test_exact_rigid_constraints_commute_and_cantilever_limit():
    m=study.JointLeading('RIGID')
    A=null_space(np.vstack((m.base_constraint,m.jump,m.D)))
    B=null_space(np.vstack((m.base_constraint,m.D,m.jump)))
    assert A.shape[1]==46
    assert np.max(abs(A@A.T-B@B.T))<1e-13
    r=next(r for r in rows('summary.csv') if r['order']=='24' and r['k_joint']=='RIGID')
    assert float(r['zeta0'])==pytest.approx(CFG['cantilever_zeta'],rel=1e-12)


def test_prefix_counts_and_independent_ode():
    checks=json.loads((FOLDER/'independent_checks.json').read_text())
    assert len(checks)==57
    for point in checks:
        for mode in point['checks']:
            if mode['exact_zero']:assert mode['zero_nullity']==mode['count_above']==1
            else:
                assert mode['count_lower']['count']==mode['position']-1
                assert mode['count_upper']['count']==mode['position']
                if 'frequency_error' in mode:assert mode['frequency_error']<1e-7


def test_coefficients_reaction_and_preserved_k01_reference():
    coefficients=json.loads((FOLDER/'coefficients.json').read_text())
    for r in coefficients:
        assert r['a']==pytest.approx(r['reaction']**2,rel=1e-14)
        assert r['a']==pytest.approx(r['linear']['a'],rel=1e-8)
        assert r['b']==pytest.approx(r['linear']['Gsecond']/(2*r['linear']['Gprime']**3),rel=1e-14)
        assert max(abs(r['linear']['G']),r['linear']['residual_u'],r['linear']['residual_v'])<1e-9
        if r['k_joint']==.1:
            assert r['a']==pytest.approx(1.649206441273272,rel=1e-12)
            assert r['b']==pytest.approx(-.639617945455777,rel=1e-12)


def test_threshold_sensitivities_energies_and_limits():
    for r in rows('sensitivities.csv'):
        assert abs(float(r['level_residual']))<1e-9
        assert 0<float(r['P_D'])<1
        assert float(r['beta_H5_deg'])/float(r['beta_H6_deg'])==pytest.approx(math.sqrt(2),rel=1e-14)
        if r['k_joint']=='RIGID':
            assert r['d_tau_d_k']==r['d_log_beta_d_log_k']==''
            assert float(r['P_J'])==float(r['P_J_infty'])==0
        else:
            derivative=(float(r['q'])**2*float(r['J_infty_squared'])-float(r['J_squared']))/float(r['D_squared'])
            assert derivative==pytest.approx(float(r['d_tau_d_k']),rel=1e-13)
            assert derivative>0
            if float(r['k_joint'])>0:
                assert float(r['P_J_infty'])>float(r['P_J'])
                assert float(r['d_log_beta_d_log_k'])==pytest.approx(float(r['log_sensitivity_energy']),rel=1e-10)
            else:assert r['d_log_beta_d_log_k']==''
    for r in rows('finite_differences.csv'):assert float(r['relative_error'])<CFG['policy']['sensitivity_relative_tol']


def test_additional_analytic_identities_from_archived_states():
    r=supplemental_checks()
    assert r['concavity_tangent_max_violation']<1e-9
    assert r['neighbour_minimum_zeta_increment']>=-1e-9
    assert r['prefix_frequency_convergence']<CFG['policy']['frequency_convergence_tol']
    assert max(x['relative_error'] for x in r['a_difference_checks'])<5e-5
    assert max(abs(x['rank_two_secular_residual']) for x in r['rank_two_checks'])<1e-9
    assert max(x['mu_prime_identity_relative_error'] for x in r['a_derivatives'])<1e-9


def test_frozen_policy_all_physical_records_and_saved_forms():
    records=json.loads((FOLDER/'point_records.json').read_text())
    assert len(records)==81
    for r in records:
        d=r['diagnostics']
        for name in ['strong_ode','natural_boundary','essential','matrix_residual','kkt']:assert d[name]<1e-9
        assert d['mass_gram']<5e-7 and d['energy']<5e-8
        for mode in r['per_mode']:
            if mode['P_J'] is not None:assert mode['P_J']+mode['P_D']+mode['P_bending']==pytest.approx(1,abs=5e-8)
    assert len(rows('first_forms.csv'))==6*3*2*80


def test_cached_workflow_does_not_solve(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('cache must not solve spectra')
    monkeypatch.setattr(study.JointLeading,'solve',forbidden)
    assert study.run()==FOLDER
