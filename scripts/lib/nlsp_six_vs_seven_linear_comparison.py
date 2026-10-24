"""Bounded NLSP linear comparison; no nonlinear action construction or evolution.

Six fields retain EA, physical nu and ordinary first amplitude order for u.
All coefficient selection is in the labelled, unwhitened physical space.
The old symmetric generalized eigensolver and immutable bundle utilities are
reused. Historical seven-field eigenpairs are verified before interpretation.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss
from scipy.linalg import solve_triangular
from scipy.optimize import linear_sum_assignment
import scipy

from scripts.lib import nlsp_coupled_sensitivity_pilot as io
from scripts.lib.weakly_nonlinear_coupled_spatial_dynamics import (
    CoupledSpatialGalerkin, EndpointArm, rectangular_coefficients,
)
from scripts.lib.weakly_nonlinear_spatial_rod import RodCoefficients, FIELD_ORDER
from scripts.lib.reddy_inplane_geometry import reddy_inplane_geometry, k

ROOT = io.ROOT
CONFIG = ROOT / 'data/input/nlsp_six_vs_seven_linear_comparison.json'
OUTPUT = ROOT / 'results/nlsp_six_vs_seven_linear_comparison'
VERSION = 'nlsp-six-seven-linear-v1'
KERNELS = [f'scripts/lib/{name}.py' for name in (
    'weakly_nonlinear_spatial_rod', 'weakly_nonlinear_planar_dynamics',
    'weakly_nonlinear_spatial_dynamics', 'weakly_nonlinear_coupled_spatial_dynamics',
    'reddy_inplane_geometry', 'yartsev_ch2_monoclinic_rod')]
CODE = KERNELS + [
    'scripts/lib/nlsp_six_vs_seven_linear_comparison.py',
    'scripts/lib/nlsp_coupled_sensitivity_pilot.py',
    'scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py',
    'data/input/nlsp_six_vs_seven_linear_comparison.json',
]


def relative(a, b):
    return float(np.linalg.norm(a-b) / max(np.linalg.norm(a), np.linalg.norm(b), 1e-300))


def require(value, gate, label):
    if not np.isfinite(value) or value > gate:
        raise ArithmeticError(f'{label}: {value:.6g} > {gate:.6g}')
    return float(value)


class LinearCoupled:
    """Endpoint-lifted physical coordinates; same order as the original assembly.

    Local order: FIELD_ORDER[:nf]. Global order: all arm-1 interiors,
    all arm-2 interiors, then shared d_X,d_Y,d_Z,alpha_X,alpha_Y,alpha_Z[,c].
    EndpointArm._raw_basis is reused without initializing a nonlinear arm.
    W maps mass-scaled coordinates to these raw physical coefficients.
    """
    eigenpairs = CoupledSpatialGalerkin.eigenpairs

    def __init__(self, p, nf, coefficients=None):
        if nf not in (6, 7) or p not in (12, 14, 16):
            raise ValueError('Only six/seven fields and p12/p14/conditional p16')
        self.p, self.nf = p, nf
        self.lengths = (.9, 1.1)
        self.coefficients = coefficients or tuple(rectangular_coefficients(.05, .10, L)[0]
                                                   for L in self.lengths)
        self.fields = FIELD_ORDER[:nf]
        self.ndof = nf*(2*p-1)
        self.joint = np.arange(2*nf*(p-1), self.ndof)
        self.labels = [f'arm{i+1}:{f}:{j}' for i in range(2) for f in self.fields
                       for j in range(p-1)] + ['joint:'+f for f in
                       ('dX','dY','dZ','aX','aY','aZ','c')[:nf]]
        geometry = reddy_inplane_geometry(15.)
        self.bases = np.array([np.column_stack((a.t, a.n, k)) for a in
                               (geometry.arm1, geometry.arm2)])
        self.maps, self.transforms, self.arm_K, self.arm_M = [], [], [], []
        self.Kraw = np.zeros((self.ndof, self.ndof))
        self.Mraw = np.zeros_like(self.Kraw)
        groups = {'plane': [], 'out': []}
        for i, (L, c, B) in enumerate(zip(self.lengths, self.coefficients, self.bases)):
            T = np.zeros((nf, nf)); T[:3,:3] = B
            T[3:6,3:6] = B @ np.diag([1., -1., 1.])
            if nf == 7:
                T[6,6] = 1.
            P = np.zeros((nf*p, self.ndof))
            for f in range(nf):
                ids = np.arange(i*nf*(p-1)+f*(p-1), i*nf*(p-1)+(f+1)*(p-1))
                P[f*p:f*p+p-1, ids] = np.eye(p-1)
                P[(f+1)*p-1, self.joint] = T.T[f]
                groups['out' if f in (2,3,4) else 'plane'].extend(ids)
            self.maps.append(P); self.transforms.append(T)
            xi, wg = leggauss(2*p+3)
            x, weights = (xi+1)*L/2, wg*L/2
            val, der = self.basis(i, x), self.basis(i, x, 1)
            F = [np.zeros((len(x), nf*p)) for _ in range(nf)]
            D = [a.copy() for a in F]
            for f in range(nf):
                F[f][:,f*p:(f+1)*p] = val
                D[f][:,f*p:(f+1)*p] = der
            # Assemble V6 directly with EA; add the explicit square only for V7.
            strains = [(c.C*(1-c.nu**2), D[0]), (c.S, D[1]-F[5]),
                       (c.S, D[2]-F[4]), (c.Bp, D[5]), (c.Bb, D[4]), (c.CT, D[3])]
            if nf == 7:
                strains += [(c.C, F[6]+c.nu*D[0]), (c.H, D[6])]
            Ka = sum(coef*A.T@(weights[:,None]*A) for coef,A in strains)
            masses = (c.m,c.m,c.m,c.jp+c.jb,c.jb,c.jp,c.jp)[:nf]
            Ma = sum(coef*A.T@(weights[:,None]*A) for coef,A in zip(masses,F))
            self.arm_K.append(Ka); self.arm_M.append(Ma)
            self.Kraw += P.T@Ka@P; self.Mraw += P.T@Ma@P
        groups['plane'].extend(self.joint[[0,1,5]+([6] if nf == 7 else [])])
        groups['out'].extend(self.joint[[2,3,4]])
        self.plane = np.array(sorted(groups['plane']))
        self.out = np.array(sorted(groups['out']))
        self.W = np.zeros_like(self.Mraw)
        for ids in (self.plane, self.out):
            lower = np.linalg.cholesky(self.Mraw[np.ix_(ids,ids)])
            self.W[np.ix_(ids,ids)] = solve_triangular(lower.T, np.eye(len(ids)), lower=False)
        self.K = self.W.T@self.Kraw@self.W
        self.M0 = self.W.T@self.Mraw@self.W

    def basis(self, arm, points, derivative=0):
        template = SimpleNamespace(interior_degree=self.p, length=self.lengths[arm])
        return EndpointArm._raw_basis(template, points, derivative)

    def fields_at(self, raw, arm, points, derivative=0):
        """Return (point, field, mode) for one or several raw vectors."""
        local = (self.maps[arm] @ np.atleast_2d(raw).reshape(self.ndof, -1)).reshape(
            self.nf, self.p, -1)
        return np.einsum('xp,fpm->xfm', self.basis(arm, points, derivative), local)

    def common_profiles(self, raw, nq=41):
        xi, wg = leggauss(nq)
        profiles, weights = [], []
        for i,(L,c) in enumerate(zip(self.lengths,self.coefficients)):
            profiles.append(self.fields_at(raw,i,(xi+1)*L/2)[:,[0,1,5],:])
            weights.append(wg[:,None]*L/2*np.array([c.m,c.m,c.jp]))
        return np.array(profiles), np.array(weights)


def field_embedding(small, large):
    """Label-based embedding BEFORE whitening; no c constraint on u_s."""
    return np.eye(large.ndof)[:, [large.labels.index(label) for label in small.labels]]


def independent_energy(d, raw, nq):
    """Expanded original densities evaluated on fields, independent of strain matrices."""
    xi,wg = leggauss(nq); potential=kinetic=0.
    for i,(L,c) in enumerate(zip(d.lengths,d.coefficients)):
        x=(xi+1)*L/2; w=wg*L/2
        f=d.fields_at(raw,i,x)[:,:,0]; s=d.fields_at(raw,i,x,1)[:,:,0]
        axial = c.C*(1-c.nu**2)*s[:,0]**2 if d.nf==6 else c.C*(
            s[:,0]**2+2*c.nu*s[:,0]*f[:,6]+f[:,6]**2)+c.H*s[:,6]**2
        V=axial+c.S*((s[:,1]-f[:,5])**2+(s[:,2]-f[:,4])**2)+c.Bp*s[:,5]**2+c.Bb*s[:,4]**2+c.CT*s[:,3]**2
        T=c.m*np.sum(f[:,:3]**2,axis=1)+(c.jp+c.jb)*f[:,3]**2+c.jb*f[:,4]**2+c.jp*f[:,5]**2
        if d.nf==7:
            T+=c.jp*f[:,6]**2
        potential+=w@V; kinetic+=w@T
    return float(potential),float(kinetic)


def matrix_checks(d, cfg):
    result={}
    for name,A in [('K',d.Kraw),('M',d.Mraw)]:
        result[name+'_symmetry']=require(relative(A,A.T),cfg['symmetry_relative_gate'],name+' symmetry')
        factor=np.linalg.cholesky(A)
        result[name+'_minimum_Cholesky_pivot']=float(np.min(np.diag(factor)))
    errors=[]
    rng=np.random.default_rng(20261021)
    for _ in range(8):
        raw=rng.normal(size=d.ndof)
        V,T=independent_energy(d,raw,2*d.p+7)
        errors.extend([abs(V-raw@d.Kraw@raw)/V,abs(T-raw@d.Mraw@raw)/T])
    result['independent_energy_relative']=require(max(errors),cfg['energy_relative_gate'],'energy')
    result['plane_size']=len(d.plane); result['out_size']=len(d.out)
    assert len(d.plane)==(d.nf-3)*(2*d.p-1) and len(d.out)==3*(2*d.p-1)
    result['PL_OUT_K_cross']=require(float(np.max(abs(d.Kraw[np.ix_(d.plane,d.out)]))),0.,'reflection K')
    result['PL_OUT_M_cross']=require(float(np.max(abs(d.Mraw[np.ix_(d.plane,d.out)]))),0.,'reflection M')
    # Pure interior u: no endpoint mixing; detects an erroneous retained C.
    q=np.zeros(d.ndof); q[d.labels.index('arm1:u:0')]=1.
    xi,wg=leggauss(2*d.p+7); L=d.lengths[0]
    us=d.fields_at(q,0,(xi+1)*L/2,1)[:,0,0]
    result['pure_u_stiffness_coefficient']=float(q@d.Kraw@q/(wg*L/2@us**2))
    expected=.005 if d.nf==6 else d.coefficients[0].C
    require(abs(result['pure_u_stiffness_coefficient']-expected)/expected,cfg['energy_relative_gate'],'EA or C')
    return result


def common_mac(a,b,weights):
    A=a.reshape(-1,a.shape[-1]); B=b.reshape(-1,b.shape[-1]); w=weights.ravel()
    cross=A.T@(w[:,None]*B)
    na=np.sum(w[:,None]*A*A,axis=0); nb=np.sum(w[:,None]*B*B,axis=0)
    return cross**2/(na[:,None]*nb[None,:]),cross,na,nb


def modal_checks(d,omega,vectors,cfg):
    raw=d.W@vectors
    residual=[]; boundary=[]; joint_weak=[]
    for j in range(min(6,len(omega))):
        q=raw[:,j]; lam=omega[j]**2
        r=d.Kraw@q-lam*d.Mraw@q
        denominator=np.linalg.norm(d.Kraw@q)+lam*np.linalg.norm(d.Mraw@q)
        residual.append(np.linalg.norm(r)/denominator)
        joint_weak.append(np.linalg.norm(r[d.joint])/denominator)
        ends=[d.fields_at(q,i,np.array([0.,L]))[:,:,0] for i,L in enumerate(d.lengths)]
        scale=max(np.max(abs(e)) for e in ends)
        defects=[np.max(abs(e[0])) for e in ends]
        defects.append(np.max(abs(d.transforms[0]@ends[0][1]-d.transforms[1]@ends[1][1])))
        boundary.append(max(defects)/max(scale,1e-300))
    require(max(residual),cfg['eigen_residual_gate'],'raw eigen residual')
    require(max(boundary),cfg['boundary_relative_gate'],'essential BC and joint')
    gram=relative(raw.T@d.Mraw@raw,np.eye(len(omega)))
    require(gram,cfg['orthogonality_gate'],'mass orthogonality')
    profiles,weights=d.common_profiles(raw[:,:6],cfg['common_quadrature_points'])
    common=np.sum(weights[:,:,:,None]*profiles**2,axis=(0,1,2))
    pu=np.sum(weights[:,:,0,None]*profiles[:,:,0,:]**2,axis=(0,1))/common
    return raw,profiles,weights,pu,dict(eigen_residual=residual,boundary_relative=boundary,
        discrete_joint_balance_relative=joint_weak,mass_orthogonality=gram,
        contraction_mass_fraction=(1-common).tolist() if d.nf==7 else [0.]*6)


def source_inputs(cfg):
    artifacts={}; manifests={}
    for degree,folder in cfg['sources'].items():
        path=ROOT/folder; manifest=io.read(path/'manifest.json')
        manifests[degree]={'path':folder,'sha256':io.sha(path/'manifest.json')}
        hashes=manifest.get('artifact_hashes',manifest.get('artifacts',{}))
        for name in ('linear_reference.npz','config.json','provenance.json'):
            digest=io.sha(path/name)
            if hashes.get(name)!=digest:
                raise ValueError('Historical source hash mismatch: '+str(path/name))
            artifacts[(path/name).relative_to(ROOT).as_posix()]=digest
        old=io.read(path/'config.json')
        for key in ('lengths','beta_deg','mu','sections','material'):
            if old[key]!=cfg[key]:
                raise ValueError('Historical parameter mismatch: '+key)
        if old['degree']!=int(degree):
            raise ValueError('Historical degree mismatch')
        old_hashes=io.read(path/'provenance.json')['code_hashes']
        for name in KERNELS:
            if old_hashes[name]!=io.sha(ROOT/name):
                raise ValueError('Historical physical implementation changed: '+name)
    return artifacts,manifests


def seven_modes(d,cfg):
    if str(d.p) not in cfg['sources']:
        omega,v=d.eigenpairs(); full=np.zeros((d.ndof,len(omega))); full[d.plane]=v
        return omega,full,{'origin':'new linear eigensolve (conditional p16)','new_eigensolves':1}
    path=ROOT/cfg['sources'][str(d.p)]/'linear_reference.npz'
    with np.load(path) as z:
        require(relative(d.K,z['K']),cfg['matrix_relative_gate'],'saved whitened K')
        require(relative(d.M0,z['M0']),cfg['matrix_relative_gate'],'saved whitened M')
        assert np.array_equal(d.plane,z['plane']) and np.array_equal(d.out,z['out'])
        omega=z['planar_omega' if d.p==12 else 'omega'].copy()
        vectors=z['planar_Phi' if d.p==12 else 'Phi'].copy()
        result=dict(origin=str(path.relative_to(ROOT)),new_eigensolves=0,
                    K_relative=relative(d.K,z['K']),M_relative=relative(d.M0,z['M0']))
    references={12:(.434820730491217,1.306329394911),14:(.434820604535,1.306328966049)}
    result['rounded_reference_absolute']=float(np.max(abs(omega[[1,3]]-references[d.p])))
    require(result['rounded_reference_absolute'],cfg['reference_frequency_absolute_gate'],'seven reference')
    return omega,vectors,result


def save_modes(bundle,d,omega,vectors,raw,profiles,weights):
    grid=np.linspace(0.,1.,501)
    np.savez_compressed(bundle/f'modes_{d.nf}_p{d.p}.npz',
        omega=omega,vectors=vectors,raw_vectors=raw,W=d.W,raw_maps=np.array(d.maps),
        Kraw=d.Kraw,Mraw=d.Mraw,K=d.K,M0=d.M0,plane=d.plane,out=d.out,
        labels=np.array(d.labels),fields=np.array(d.fields),lengths=d.lengths,
        bases=d.bases,transforms=np.array(d.transforms),p=d.p,nq=2*d.p+3,
        common_profiles=profiles,common_weights=weights,grid_fraction=grid,
        grid_profiles=np.array([d.fields_at(raw[:,:6],i,grid*L) for i,L in enumerate(d.lengths)]),
        coefficients=np.array([list(asdict(c).values()) for c in d.coefficients]),
        coefficient_names=np.array(list(asdict(d.coefficients[0]))))


def compare_degree(bundle,degree,cfg):
    started=time.perf_counter(); d6,d7=LinearCoupled(degree,6),LinearCoupled(degree,7)
    checks={'6':matrix_checks(d6,cfg),'7':matrix_checks(d7,cfg)}
    for name in ('Kraw','Mraw'):
        A=getattr(d6,name)[np.ix_(d6.out,d6.out)]; B=getattr(d7,name)[np.ix_(d7.out,d7.out)]
        checks['OUT_'+name+'_relative']=require(relative(A,B),cfg['matrix_relative_gate'],'OUT '+name)
        checks['OUT_'+name+'_maximum_absolute']=float(np.max(abs(A-B)))
    # nu=0 changes G,S,H,CT consistently; no spectral solve for this control.
    coeff0=tuple(RodCoefficients.rectangular(1.,1.,0.,.05,.10,c.CT*1.3) for c in d6.coefficients)
    z6,z7=LinearCoupled(degree,6,coeff0),LinearCoupled(degree,7,coeff0)
    E=field_embedding(z6,z7); cids=[i for i,l in enumerate(z7.labels) if ':c:' in l or l=='joint:c']
    for name in ('Kraw','Mraw'):
        A=getattr(z6,name); B=getattr(z7,name)
        checks['nu0_'+name+'_common']=require(relative(A,E.T@B@E),cfg['matrix_relative_gate'],'nu0 common')
        checks['nu0_'+name+'_cross']=require(float(np.max(abs(E.T@B[:,cids]))),0.,'nu0 decoupling')
    w7,v7,checks['seven_source']=seven_modes(d7,cfg)
    w6,v=d6.eigenpairs(); v6=np.zeros((d6.ndof,len(w6))); v6[d6.plane]=v
    r6,f6,weights,pu6,checks['6']['modal']=modal_checks(d6,w6,v6,cfg)
    r7,f7,_,pu7,checks['7']['modal']=modal_checks(d7,w7,v7,cfg)
    mac,cross,_,_=common_mac(f7,f6,weights)
    rows,cols=linear_sum_assignment(-mac); assignment=dict(zip(rows,cols))
    matches=[]
    for j in (1,3):
        target=assignment[j]; competitors=[i for i in range(6) if i!=target]
        other=max(competitors,key=lambda i:mac[j,i])
        gap=min(abs(w6[target]-np.delete(w6,target)))/w6[target]
        gap7=min(abs(w7[j]-np.delete(w7,j)))/w7[j]
        accepted=(mac[j,target]>=cfg['mode_identity_MAC_gate'] and
                  mac[j,target]-mac[j,other]>=cfg['mode_identity_margin_gate'] and
                  min(gap,gap7)>cfg['cluster_relative_gap'])
        matches.append(dict(p=degree,source_target=f'PL{j+1}',position_7=j+1,position_6=target+1,
            omega_7=w7[j],omega_6=w6[target],delta_omega=w6[target]-w7[j],
            relative_delta_omega=(w6[target]-w7[j])/w7[j],MAC=mac[j,target],
            competitor_position_6=other+1,competitor_MAC=mac[j,other],
            status='UNAMBIGUOUS' if accepted else 'UNRESOLVED',
            nearest_relative_gap_6=gap,nearest_relative_gap_7=gap7,
            P_u_7=pu7[j],P_u_6=pu6[target],
            c_mass_fraction_7=checks['7']['modal']['contraction_mass_fraction'][j]))
    np.savez_compressed(bundle/f'matching_p{degree}.npz',MAC=mac,cross=cross,
                        assignment=np.array([assignment[j] for j in range(6)]))
    spectrum=[]
    for d,w,vectors,raw,profiles,pu in ((d6,w6,v6,r6,f6,pu6),(d7,w7,v7,r7,f7,pu7)):
        save_modes(bundle,d,w,vectors,raw,profiles,weights)
        for j in range(6):
            spectrum.append(dict(model=d.nf,p=degree,family='PL',sorted_position=j+1,
                omega=w[j],eigen_residual=checks[str(d.nf)]['modal']['eigen_residual'][j],P_u=pu[j]))
    checks['elapsed_seconds']=time.perf_counter()-started
    return matches,spectrum,checks


def resonances(matches):
    result=[]
    for degree in sorted({r['p'] for r in matches}):
        a,b=[r for r in matches if r['p']==degree]
        row=dict(p=degree,position_6_low=a['position_6'],position_6_high=b['position_6'],
                 position_7_low=2,position_7_high=4,confirmed=a['status']==b['status']=='UNAMBIGUOUS')
        for nf in (6,7):
            row[f'Delta13_{nf}']=b[f'omega_{nf}']-3*a[f'omega_{nf}'] if row['confirmed'] else None
            row[f'r13_{nf}']=b[f'omega_{nf}']/(3*a[f'omega_{nf}'])-1 if row['confirmed'] else None
            for name in ('eta2','eta4','etaDelta'):
                row[f'{name}_{nf}']=None
        row['delta_Delta13']=row['Delta13_6']-row['Delta13_7'] if row['confirmed'] else None
        row['eta_diff']=None; row['difference_over_indicator']=None
        if result and row['confirmed'] and result[-1]['confirmed']:
            old=[r for r in matches if r['p']==result[-1]['p']]
            for nf in (6,7):
                for label,new,prev in [('eta2',a,old[0]),('eta4',b,old[1])]:
                    row[f'{label}_{nf}']=abs(new[f'omega_{nf}']-prev[f'omega_{nf}'])
                row[f'etaDelta_{nf}']=row[f'eta4_{nf}']+3*row[f'eta2_{nf}']
            row['eta_diff']=row['etaDelta_6']+row['etaDelta_7']
            row['difference_over_indicator']=abs(row['delta_Delta13'])/max(row['eta_diff'],1e-300)
        result.append(row)
    return result


def write_csv(path,rows):
    with path.open('w',newline='',encoding='utf8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def render_saved(bundle):
    """Read saved profiles/tables only; never assemble or solve."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    s=io.cache(bundle); high=max(r['p'] for r in s['matches'])
    matches=[r for r in s['matches'] if r['p']==high]
    fig,axes=plt.subplots(3,2,figsize=(11,8),sharex=True)
    colors=('#1764a0','#d56724')
    arrays={nf:dict(np.load(bundle/f'modes_{nf}_p{high}.npz')) for nf in (6,7)}
    for col,row in enumerate(matches):
        for nf in (7,6):
            z=arrays[nf]; j=row[f'position_{nf}']-1
            q=z['grid_profiles'][:,:,:,j].copy()
            # Exact extrema of the represented w polynomial, not component-wise normalization.
            maximum=0.
            for i in range(2):
                local=(z['raw_maps'][i]@z['raw_vectors'][:,j]).reshape(nf,high)[1]
                co=np.zeros(high+1); co[:high-1]+=local[:-1]; co[2:]-=local[:-1]
                co[:2]+=local[-1]/2
                poly=Legendre(co); roots=poly.deriv().roots()
                roots=roots[np.isreal(roots)].real; pts=np.r_[-1.,1.,roots[abs(roots)<=1]]
                maximum=max(maximum,float(np.max(abs(poly(pts)))))
            q/=maximum
            if nf==6:
                with np.load(bundle/f'matching_p{high}.npz') as match:
                    q*=1 if match['cross'][row['position_7']-1,j]>=0 else -1
            for i,L in enumerate(z['lengths']):
                for a,f in enumerate((0,1,5)):
                    axes[a,col].plot(z['grid_fraction']*L,q[i,:,f],
                        color=colors[i],linestyle='-' if nf==7 else '--',
                        label=f'{nf} fields, arm {i+1}')
        axes[0,col].set_title(f"PL{row['position_7']} (7) / PL{row['position_6']} (6), p{high}")
        for a,label in enumerate(('u / max|w|','w / max|w|','theta / max|w|')):
            axes[a,col].set_ylabel(label); axes[a,col].grid(alpha=.25)
        axes[2,col].set_xlabel('s_i: outer clamp -> joint (model length)')
    axes[0,0].legend(fontsize=8,ncol=2)
    fig.suptitle('Common normalization: max over both arms |w| = 1; mass-overlap sign alignment',fontsize=11)
    fig.tight_layout(); fig.savefig(bundle/'target_shapes.png',dpi=200); fig.savefig(bundle/'target_shapes.pdf'); plt.close(fig)
    lines=['# Saved linear comparison', '', 'No eigensolves in report/plot generation.', '',
           '|p|target|position 6|omega 7|omega 6|delta omega|MAC|Pu 7|Pu 6|',
           '|---:|---|---:|---:|---:|---:|---:|---:|---:|']
    for r in s['matches']:
        lines.append(f"|{r['p']}|{r['source_target']}|{r['position_6']}|{r['omega_7']:.12g}|{r['omega_6']:.12g}|{r['delta_omega']:.8g}|{r['MAC']:.12g}|{r['P_u_7']:.9g}|{r['P_u_6']:.9g}|")
    lines+=['','|p|Delta 7|Delta 6|difference|eta diff|','|---:|---:|---:|---:|---:|']
    for r in s['resonances']:
        lines.append(f"|{r['p']}|{r['Delta13_7']}|{r['Delta13_6']}|{r['delta_Delta13']}|{r['eta_diff']}|")
    (bundle/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    io.finalize(bundle,io.read(bundle/'provenance.json'),s)
    return s


def run():
    cfg=io.read(CONFIG); artifacts,manifests=source_inputs(cfg)
    inputs=dict(version=VERSION,config=cfg,code_hashes={p:io.sha(ROOT/p) for p in CODE},
                source_artifacts=artifacts,source_manifests=manifests,
                numpy=np.__version__,scipy=scipy.__version__)
    fingerprint=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()[:16]
    bundle=OUTPUT/fingerprint
    if (bundle/'manifest.json').exists():
        summary=io.cache(bundle)
        if not summary['completed']:
            raise ValueError('Preserve incomplete attempt; no automatic restart')
        return bundle,summary,True
    if bundle.exists():
        raise ValueError('Preserve unmanifested attempt')
    bundle.mkdir(parents=True)
    def git(*args):
        return subprocess.check_output(['git',*args],cwd=ROOT,text=True,encoding='utf8').strip()
    provenance=dict(inputs,schema=VERSION,HEAD=git('rev-parse','HEAD'),
        cwd=str(ROOT),git_root=git('rev-parse','--show-toplevel'),branch=git('branch','--show-current'),
        status_short_at_launch=git('status','--short'),
        initial_session={'HEAD':'c6f6e4ead444a86e19c2d5d24c9e820cbc95f722','branch':'main','status_short':''},
        started_UTC=datetime.now(timezone.utc).isoformat(),interpreter=sys.executable,
        units='previous NLSP model length/time, angular omega; Lref=1',
        solver='unchanged CoupledSpatialGalerkin.eigenpairs: scipy.linalg.eigh(K,M), full PL block',
        constraints='outer values zero; B1 U1=B2 U2, B1 a1=B2 a2; seven only: c(0)=0,c1=c2',
        basis='p-1 Shen P_n-P_(n+2) plus s/L; nq=2p+3; raw=W@scaled',
        torsion=[rectangular_coefficients(.05,.10,L)[1] for L in (.9,1.1)],
        torsion_settings={'series_relative_tolerance':1e-12,'series_max_terms':100000},
        coefficients=[asdict(rectangular_coefficients(.05,.10,L)[0]) for L in (.9,1.1)])
    io.write(bundle/'provenance.json',provenance); io.write(bundle/'config.json',cfg)
    for name in CODE:
        dest=bundle/'execution_code'/name; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(ROOT/name,dest)
    started=time.perf_counter(); matches=[]; spectrum=[]; checks={}
    degrees=list(cfg['degrees']); p16_reason=None
    for degree in degrees:
        m,sp,ch=compare_degree(bundle,degree,cfg)
        matches+=m; spectrum+=sp; checks[str(degree)]=ch
        if degree==14:
            last=resonances(matches)[-1]
            if not last['confirmed']:
                p16_reason='Target correspondence unresolved at p12/p14.'
            elif abs(last['delta_Delta13'])<=cfg['refinement_indicator_factor']*last['eta_diff']:
                p16_reason='Model detuning change <= 10 times observed p12/p14 eta_diff; resolve sign/significant digits.'
            if p16_reason:
                io.write(bundle/'p16_authorized_trigger.json',dict(reason=p16_reason,p12_p14=last))
                degrees.append(16)
    summary=dict(overall='LINEAR_COMPARISON_COMPLETE',completed=True,matches=matches,
        resonances=resonances(matches),checks=checks,p16_reason=p16_reason,
        calls={'linear_eigenproblems':len(degrees)+sum(c['seven_source']['new_eigensolves'] for c in checks.values()),
               'seven_field_eigenproblems':sum(c['seven_source']['new_eigensolves'] for c in checks.values()),
               'HB':0,'Hill':0,'Floquet':0,'ODE':0,'FEM':0},numerical_seconds=time.perf_counter()-started)
    # Inter-resolution correspondence on the same physical quadrature.
    summary['resolution_MAC']={}
    for low,high in zip(degrees[:-1],degrees[1:]):
        for nf in (6,7):
            with np.load(bundle/f'modes_{nf}_p{low}.npz') as a, np.load(bundle/f'modes_{nf}_p{high}.npz') as b:
                mac=common_mac(a['common_profiles'],b['common_profiles'],a['common_weights'])[0]
                summary['resolution_MAC'][f'{nf}:{low}:{high}']=mac.tolist()
                for target in ('PL2','PL4'):
                    j=next(r[f'position_{nf}']-1 for r in matches if r['p']==low and r['source_target']==target)
                    h=next(r[f'position_{nf}']-1 for r in matches if r['p']==high and r['source_target']==target)
                    require(1-mac[j,h],1-cfg['mode_identity_MAC_gate'],'p-refinement mode identity')
    write_csv(bundle/'spectrum.csv',spectrum); write_csv(bundle/'mode_matches.csv',matches)
    write_csv(bundle/'resonance_comparison.csv',summary['resonances'])
    io.write(bundle/'verification.json',checks); io.finalize(bundle,provenance,summary)
    render_saved(bundle)
    return bundle,summary,False


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--compute',action='store_true')
    mode.add_argument('--report-only',type=Path)
    mode.add_argument('--plot-only',type=Path)
    args=parser.parse_args(argv)
    if args.compute:
        bundle,summary,cached=run()
    else:
        bundle=(args.report_only or args.plot_only).resolve()
        summary=render_saved(bundle); cached=True
    print(json.dumps(io.plain(dict(bundle=str(bundle),overall=summary['overall'],
        reused_completed=cached,new_eigensolves_this_invocation=0 if cached else summary['calls']['linear_eigenproblems'],
        matches=summary['matches'],resonances=summary['resonances'],p16_reason=summary['p16_reason'])),indent=2))
