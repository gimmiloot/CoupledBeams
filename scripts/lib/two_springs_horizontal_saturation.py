"""HS-LA3: postprocess frozen HS-LA2 matrices; no eigensolver or fit.

Separate small helper preserves source-hashed HS-LA2 workflow/cache identity.
The new contract is two indefinite linear solves plus an archived modal sum,
not another spectral runner. Reuses HS-LA2 serialization and hash helpers.
Entry: python -c "from scripts.lib.two_springs_horizontal_saturation import run; run()"
"""
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import mpmath as mp
import numpy as np
import scipy
from scipy.linalg import solve
from scripts.lib.two_springs_horizontal_spectral_structure import ROOT,sha,write,csvwrite,require

CONFIG=ROOT/'data/input/isotropic_two_springs_horizontal_leading_second_coefficient.json'
BASE=ROOT/'results/isotropic_two_springs_rl_timoshenko/horizontal_leading_second_coefficient'


def identity(cfg):
    record=dict(config=cfg,source_sha256=sha(__file__),config_sha256=sha(CONFIG),
        initial_sha256=sha(ROOT/cfg['initial_workspace']),
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__))
    return hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()[:16],record


def coefficient_from_linear(K,M,d,zeta_infty):
    """Only K0, M, d and the constrained eigenvalue enter this calculation.

    LAPACK symmetric indefinite solve (pivoting); never Cholesky on A.
    No finite-tau eigenvalue is an input.
    """
    A=K-zeta_infty*M
    u=solve(A,d,assume_a='sym');Mu=M@u
    v=solve(A,Mu,assume_a='sym')
    gp=float(u@Mu);half_gpp=float(u@M@v)
    residual=lambda x,rhs:float(np.linalg.norm(A@x-rhs)/(np.linalg.norm(A@x)+np.linalg.norm(rhs)))
    return dict(G=float(d@u),Gprime=gp,Gsecond=2*half_gpp,a=1/gp,b=half_gpp/gp**3,
        residual_u=residual(u,d),residual_v=residual(v,Mu),condition_A_2=float(np.linalg.cond(A)),
        symmetry_A=float(np.max(abs(A-A.T)))),u,v


def coefficient_from_modes(authoritative,zeta_infty,digits=80):
    """Use decimal strings of ALL authoritative pairs, never rounded vectors."""
    with mp.workdps(digits):
        z=mp.mpf(float(zeta_infty))
        poles=[mp.mpf(v) for v in authoritative['zeta']]
        dn=[mp.mpf(v) for v in authoritative['d_n']]
        require(len(dn)==len(poles),'Incomplete modal data')
        distance=min(abs(v-z) for v in poles);require(distance>0,'Original pole')
        G=mp.fsum(d*d/(p-z) for p,d in zip(poles,dn))
        gp=mp.fsum(d*d/(p-z)**2 for p,d in zip(poles,dn))
        gpp=2*mp.fsum(d*d/(p-z)**3 for p,d in zip(poles,dn))
        values=dict(G=G,Gprime=gp,Gsecond=gpp,a=1/gp,b=gpp/(2*gp**3))
        return dict(**{k:float(v) for k,v in values.items()},decimal={k:mp.nstr(v,digits) for k,v in values.items()},
            N=len(poles),negative_directions=sum(p<z for p in poles),pole_distance=float(distance))


def threshold_predictions(a,b,zeta_infty,q):
    delta=zeta_infty*(1-q*q);disc=a*a+4*b*delta
    require(delta>0 and disc>0,'Unresolved asymptotic inversion')
    # Plus sign is the large-tau inverse; the other root is diagnostic only.
    return dict(delta=delta,discriminant=disc,tau_one=a/delta,
        tau_two=(a+math.sqrt(disc))/(2*delta),tau_other=(a-math.sqrt(disc))/(2*delta))


def verify_sources(cfg):
    src=ROOT/cfg['source_archive'];manifest=json.loads((src/'manifest.json').read_text())
    require(sha(src/'manifest.json')==cfg['source_manifest_sha256'],'HS-LA2 manifest changed')
    for name,digest in manifest['artifacts'].items():require(sha(src/name)==digest,'HS-LA2 artifact '+name)
    for name,digest in cfg['frozen_files'].items():require(sha(ROOT/name)==digest,'Frozen source '+name)
    state=json.loads((ROOT/cfg['initial_workspace']).read_text())
    for name,digest in state['source_archive_files'].items():require(sha(ROOT/name)==digest,'Historical archive '+name)
    require(sha(ROOT/cfg['attachment']['copy_path'])==cfg['attachment']['sha256'],'Analytical source')
    # HS-L1 control values are read only after validating its own artifacts.
    l1=ROOT/cfg['HS_L1_archive'];m1=json.loads((l1/'manifest.json').read_text())
    for name,digest in m1['artifacts'].items():require(sha(l1/name)==digest,'HS-L1 artifact '+name)
    return src,l1


def run():
    cfg=json.loads(CONFIG.read_text());identifier,record=identity(cfg);folder=BASE/identifier
    src,l1=verify_sources(cfg)
    if (folder/'manifest.json').exists():
        m=json.loads((folder/'manifest.json').read_text());require(m['identity']==record,'Cache identity')
        for name,digest in m['artifacts'].items():require(sha(folder/name)==digest,'Artifact '+name)
        print('HS-LA3 verified',folder);return folder
    folder.mkdir(parents=True,exist_ok=True);write(folder/'contract.json',record)
    coefficients=[];vectors={};tol=cfg['tolerances']
    try:
        reactions=json.loads((src/'reaction.json').read_text())
        limits=json.loads((src/'gate_A.json').read_text())['initial']
        for p in cfg['orders']:
            key=hashlib.sha256(repr((p,0.)).encode()).hexdigest()[:16]
            with np.load(src/f'full_basis_{key}.npz') as data:
                K=data['K'].copy();M=data['M'].copy();d=data['d'].copy()
                # Rounded high-mode coefficients are deliberately never loaded.
            modes=json.loads((src/f'full_basis_{key}_80digits.json').read_text())
            r=next(v for v in reactions if v['order']==p);z=r['zeta_infty'];a=r['r_squared']
            direct,u,v=coefficient_from_linear(K,M,d,z)
            spectral=coefficient_from_modes(modes,z,cfg['arithmetic']['digits'])
            initial=next(v for v in limits if v['order']==p)
            s43=sum(1/(x-z) for x in initial['zeta0'][:4])-sum(1/(x-z) for x in initial['zeta_infty'][1:3])
            upper=a*a*s43;lower=upper-a*a/(initial['zeta0'][3]-z)
            error=abs(direct['b']/spectral['b']-1)
            row=dict(order=p,zeta_infty=z,r_infty=r['r'],a_reaction=a,linear=direct,modal=spectral,
                coefficient_methods_relative_error=error,Gprime_reaction_relative_error=abs(direct['Gprime']*a-1),
                analytic_lower=lower,analytic_upper=upper,requested_interval_pass=cfg['bounds']['requested_lower']<direct['b']<=cfg['bounds']['requested_upper'])
            coefficients.append(row);write(folder/'coefficients.json',coefficients)
            vectors[f'p{p}_u']=u;vectors[f'p{p}_v']=v
            require(spectral['N']==2*p and spectral['negative_directions']==1,'Indefinite system inertia')
            require(max(direct['residual_u'],direct['residual_v'],direct['symmetry_A'])<tol['linear_residual'],'Linear systems')
            require(max(abs(direct['G']),abs(spectral['G']))<tol['G_absolute'],'Zero of G')
            require(row['Gprime_reaction_relative_error']<tol['reaction_identity_relative'],'Reaction normalization')
            require(error<tol['coefficient_methods_relative'],'Coefficient methods')
            require(lower<direct['b']<=upper and row['requested_interval_pass'],'Analytical bound mismatch: retain result and diagnose')
        np.savez_compressed(folder/'linear_vectors.npz',**vectors)
        order_error=abs(coefficients[0]['linear']['b']/coefficients[1]['linear']['b']-1)
        require(order_error<tol['order_convergence_relative'],'Order convergence')
        # Only now read finite-tau data: a and b above cannot be fitted to them.
        with open(src/'saturation.csv',encoding='utf8') as f:tail=list(csv.DictReader(f))
        thresholds=json.loads((l1/'leading_thresholds.json').read_text())
        with open(l1/'thresholds.csv',encoding='utf8') as f:angles=list(csv.DictReader(f))
        comparisons=[];predictions=[]
        for c in coefficients:
            p=c['order'];z=c['zeta_infty'];a=c['a_reaction'];b=c['linear']['b']
            taus=[thresholds['0.9']['tau'],thresholds['0.95']['tau'],12.,24.]
            for t in taus:
                source=next(v for v in tail if int(v['order'])==p and float(v['tau'])==t)
                exact=float(source['zeta']);delta=z-exact;one=a/t;two=one+b/t**2
                require(z-one>0 and z-two>0,'Frequency approximation not positive')
                comparisons.append(dict(order=p,tau=t,zeta_archive=exact,delta_archive=delta,delta_one=one,delta_two=two,
                    delta_one_relative_error=one/delta-1,delta_two_relative_error=two/delta-1,
                    omega_one_relative_error=math.sqrt((z-one)/exact)-1,omega_two_relative_error=math.sqrt((z-two)/exact)-1,
                    residual_one=delta-one,residual_two=delta-two))
            for q in cfg['levels']:
                actual=thresholds[str(q)]['tau'];prediction=threshold_predictions(a,b,z,q)
                for construction,factor in [('H5',4),('H6',2)]:
                    ref=next(v for v in angles if v['id']==construction and float(v['q'])==q)
                    beta=float(ref['beta_lead_deg']);scale=180/math.pi*cfg['epsilon']*math.sqrt(factor)
                    require(abs(scale*math.sqrt(actual)/beta-1)<1e-12,'Archived angle normalization')
                    first=scale*math.sqrt(prediction['tau_one']);second=scale*math.sqrt(prediction['tau_two'])
                    predictions.append(dict(order=p,id=construction,q=q,tau_archive=actual,**prediction,
                        tau_one_relative_error=prediction['tau_one']/actual-1,tau_two_relative_error=prediction['tau_two']/actual-1,
                        beta_archive_deg=beta,beta_one_deg=first,beta_two_deg=second,
                        beta_one_error_deg=first-beta,beta_two_error_deg=second-beta,
                        beta_one_relative_error=first/beta-1,beta_two_relative_error=second/beta-1))
        csvwrite(folder/'asymptotic_comparison.csv',comparisons);csvwrite(folder/'threshold_predictions.csv',predictions)
        write(folder/'validation.json',dict(status='PASS',coefficient_order_relative_error=order_error,
            new_linear_system_solves=4,new_leading_eigensolves=0,new_full_model_eigensolves=0,
            reused_authoritative_spectra=2,reused_constrained_first_pairs=2,reused_finite_tau_prefixes=8,
            fit=False))
        write(folder/'manifest.json',dict(status='PASS',identity=record,artifacts={p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file()}))
    except Exception as e:
        write(folder/'failure.json',dict(status='FAIL',error=repr(e)));raise
    print('HS-LA3 PASS',folder);return folder
