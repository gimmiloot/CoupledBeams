"""Targeted linear NLSP risks; no new eigensolves, nonlinear action or ODE."""
import shutil

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from scripts.lib import nlsp_six_vs_seven_linear_comparison as n


@pytest.fixture(scope='module')
def pair():
    return n.LinearCoupled(12,6),n.LinearCoupled(12,7)


@pytest.fixture(scope='module')
def saved():
    candidates=list(n.OUTPUT.glob('*/manifest.json'))
    if not candidates:
        pytest.skip('Optional completed local linear comparison bundle required')
    bundle=candidates[-1].parent
    return bundle,n.io.cache(bundle)


def test_explicit_energy_difference_and_mass_difference(pair):
    """Independent field integrals of the exact V7-V6 and T7-T6 identities."""
    d6,d7=pair; E=n.field_embedding(d6,d7)
    rng=np.random.default_rng(731)
    q7=rng.normal(size=d7.ndof); q6=E.T@q7
    extra_v=extra_t=0.
    xi,wg=leggauss(39)
    for i,(L,c) in enumerate(zip(d7.lengths,d7.coefficients)):
        f=d7.fields_at(q7,i,(xi+1)*L/2)[:,:,0]
        s=d7.fields_at(q7,i,(xi+1)*L/2,1)[:,:,0]
        extra_v+=wg*L/2@(c.C*(f[:,6]+c.nu*s[:,0])**2+c.H*s[:,6]**2)
        extra_t+=wg*L/2@(c.jp*f[:,6]**2)
    assert q7@d7.Kraw@q7-q6@d6.Kraw@q6==pytest.approx(extra_v,rel=2e-12)
    assert q7@d7.Mraw@q7-q6@d6.Mraw@q6==pytest.approx(extra_t,rel=2e-12)


def test_EA_is_not_C_and_no_axial_slope_clamp(pair):
    d6,d7=pair
    a=n.matrix_checks(d6,n.io.read(n.CONFIG))
    assert a['pure_u_stiffness_coefficient']==pytest.approx(.005,rel=2e-12)
    q=np.zeros(d6.ndof); q[d6.labels.index('arm1:u:0')]=1.
    assert d6.fields_at(q,0,np.array([0.]))[0,0,0]==0.
    assert abs(d6.fields_at(q,0,np.array([0.]),1)[0,0,0])>1.
    E=n.field_embedding(d6,d7)
    assert (E@q)@d7.Kraw@(E@q)>q@d6.Kraw@q


def test_nu_zero_common_block_and_decoupling_without_spectrum(pair):
    coeff=tuple(n.RodCoefficients.rectangular(1.,1.,0.,.05,.10,c.CT*1.3)
                for c in pair[0].coefficients)
    a,b=n.LinearCoupled(12,6,coeff),n.LinearCoupled(12,7,coeff)
    E=n.field_embedding(a,b)
    cids=[i for i,l in enumerate(b.labels) if ':c:' in l or l=='joint:c']
    for name in ('Kraw','Mraw'):
        np.testing.assert_array_equal(getattr(a,name),E.T@getattr(b,name)@E)
        np.testing.assert_array_equal(E.T@getattr(b,name)[:,cids],0.)


def test_MAC_uses_common_mass_only_and_is_scale_invariant():
    a=np.array([[[[1.,0.],[0.,1.],[1.,1.]]]])
    b=a*np.array([-7.,.3])
    weights=np.array([[[2.,2.,.1]]])
    mac=n.common_mac(a,b,weights)[0]
    np.testing.assert_allclose(np.diag(mac),1.,atol=1e-15)
    assert mac[0,1]==pytest.approx(.1**2/2.1**2)


def test_saved_dimensions_out_identity_recovery_and_checks(saved):
    bundle,s=saved
    for p in (12,14,16):
        with np.load(bundle/f'modes_6_p{p}.npz') as a,np.load(bundle/f'modes_7_p{p}.npz') as b:
            assert len(a['plane'])==3*(2*p-1)
            assert len(b['plane'])==4*(2*p-1)
            for key in ('Kraw','Mraw'):
                np.testing.assert_array_equal(a[key][np.ix_(a['out'],a['out'])],b[key][np.ix_(b['out'],b['out'])])
            for z in (a,b):
                np.testing.assert_allclose(z['raw_vectors'],z['W']@z['vectors'],atol=1e-14)
                np.testing.assert_allclose(z['W'].T@z['Mraw']@z['W'],z['M0'],atol=1e-14)
                nf=len(z['fields'])
                local=np.array([P@z['raw_vectors'][:,:6] for P in z['raw_maps']]).reshape(2,nf,p,6)
                np.testing.assert_allclose(z['transforms'][0]@local[0,:,-1],z['transforms'][1]@local[1,:,-1],atol=1e-13)
        for nf in ('6','7'):
            c=s['checks'][str(p)][nf]
            assert c['independent_energy_relative']<2e-12
            assert c['K_minimum_Cholesky_pivot']>0 and c['M_minimum_Cholesky_pivot']>0
            assert max(c['modal']['eigen_residual'])<2e-8
            assert max(c['modal']['boundary_relative'])<2e-10


def test_saved_matching_and_bounded_refinement(saved):
    _,s=saved
    assert s['calls']==dict(linear_eigenproblems=4,seven_field_eigenproblems=1,HB=0,Hill=0,Floquet=0,ODE=0,FEM=0)
    assert s['p16_reason']
    assert s['resonances'][1]['difference_over_indicator']<10
    assert s['resonances'][2]['difference_over_indicator']>25
    for r in s['matches']:
        assert r['status']=='UNAMBIGUOUS'
        assert r['MAC']-r['competitor_MAC']>.99
        assert r['delta_omega']<0 and 0<r['P_u_7']<r['P_u_6']<1
    assert s['matches'][-2]['omega_6']==pytest.approx(.4348073138328,abs=2e-11)
    assert s['matches'][-1]['omega_6']==pytest.approx(1.3062932741739,abs=2e-11)
    assert s['resonances'][-1]['delta_Delta13']==pytest.approx(4.19503825e-6,abs=2e-11)


def test_completed_replay_and_reporting_never_assemble_or_solve(saved,monkeypatch,tmp_path):
    bundle,_=saved
    def forbidden(*args,**kwargs):
        raise AssertionError('A completed replay must not assemble or solve')
    monkeypatch.setattr(n.LinearCoupled,'__init__',forbidden)
    monkeypatch.setattr(n.LinearCoupled,'eigenpairs',forbidden)
    monkeypatch.setattr(n,'compare_degree',forbidden)
    found,_,cached=n.run()
    assert found==bundle and cached
    copy=tmp_path/'bundle'; shutil.copytree(bundle,copy)
    n.render_saved(copy)
    n.io.cache(copy)


def test_tampered_result_is_rejected(saved,tmp_path):
    bundle,_=saved; copy=tmp_path/'tampered'; shutil.copytree(bundle,copy)
    with (copy/'spectrum.csv').open('a',encoding='utf8') as stream:
        stream.write('changed\n')
    with pytest.raises(ValueError,match='artifact changed'):
        n.io.cache(copy)
