"""Four bounded OUT Hill decompositions on saved six-field p14 HB orbits.

No new physical law, linear spectrum, HB solve, ODE, PL spectrum or initial-
condition projection. All matching is in common physical (eta, velocity).
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np
from scipy.linalg import block_diag, cho_factor, cho_solve, eig

from scripts.lib import nlsp_six_vs_seven_planar_orbit as orbit
from scripts.lib import nlsp_instability_feasibility as variation

io=orbit.io; hb=orbit.hb; linear=orbit.linear; ROOT=orbit.ROOT
CONFIG=ROOT/'data/input/nlsp_six_vs_seven_out_stability.json'
OUTPUT=ROOT/'results/nlsp_six_vs_seven_out_stability'
VERSION='six-seven-out-stability-v1'
CODE=list(dict.fromkeys(orbit.CODE+[
    'scripts/lib/nlsp_six_vs_seven_out_stability.py',
    'scripts/lib/nlsp_instability_feasibility.py',
    'data/input/nlsp_six_vs_seven_out_stability.json']))


def cdict(z):
    return dict(real=float(np.real(z)),imag=float(np.imag(z)))


def uncomplex(row):
    return complex(row['real'],row['imag'])


def check_time(deadline):
    if time.perf_counter()>deadline:raise TimeoutError('Bounded OUT stage exhausted 1800s')


def sources(cfg):
    selected={}; manifests={}
    names={
        'orbit_parent':['config.json','provenance.json','HB11.npz','HB13.npz','HB11.json','HB13.json','linear_6_p14.npz','recovery.npz','six_field_action.json'],
        'linear_parent':['config.json','modes_7_p14.npz'],
        'seven_parent':['config.json','linear_reference.npz','OUT_Fourier.npz','OUT_k81_N9.npz','OUT_k81_N9.json','OUT_checks.json','physical_Floquet_and_displacement.npz','projection_checks.json','HB13.json','HB13.npz']}
    for key,files in names.items():
        folder=ROOT/cfg[key]; manifest=io.read(folder/'manifest.json')
        hashes=manifest.get('artifact_hashes',manifest.get('artifacts',{}))
        for name in files:
            digest=io.sha(folder/name)
            if hashes.get(name)!=digest:raise ValueError('Source artifact changed: '+str(folder/name))
            selected[(folder/name).relative_to(ROOT).as_posix()]=digest
        manifests[key]=dict(path=cfg[key],sha256=io.sha(folder/'manifest.json'))
        other=io.read(folder/'config.json')
        for par in ('lengths','sections','beta_deg','material'):
            if other[par]!=cfg[par]:raise ValueError('Source parameter mismatch: '+par)
    # Frozen action and physical code provenance are checked through the orbit cache.
    io.cache(ROOT/cfg['orbit_parent'])
    parent=io.read(ROOT/cfg['orbit_parent']/'provenance.json')
    for name in linear.KERNELS+['scripts/lib/nlsp_six_vs_seven_planar_orbit.py']:
        if io.sha(ROOT/name)!=parent['code_hashes'][name]:raise ValueError('Frozen physical source changed: '+name)
    return selected,manifests


def out_transfer(source,target):
    """Mass Galerkin transfer of v,Phi,psi before any modal truncation."""
    cross=np.zeros((81,81)); gram=np.zeros((81,81))
    xi,w=np.polynomial.legendre.leggauss(41)
    for a,b,Ps,Pt in zip(source.arms,target.arms,source.maps,target.maps):
        x=(xi+1)*a.length/2; weights=w*a.length/2
        for f,m in (('v',b.coefficients.m),('Phi',b.coefficients.jp+b.coefficients.jb),('psi',b.coefficients.jb)):
            A=a.basis_at(x)[f]@Ps[a.slices[f]][:,source.out]
            B=b.basis_at(x)[f]@Pt[b.slices[f]][:,target.out]
            cross+=m*B.T@(weights[:,None]*A);gram+=m*B.T@(weights[:,None]*B)
    return np.linalg.solve(gram,cross)


def load_problem(cfg):
    base=ROOT/cfg['orbit_parent']; ref=ROOT/cfg['seven_parent']
    old=orbit.load_action(io.read(base/'config.json'))
    d=orbit.load_operator(base)
    seven=orbit.SavedActionSystem(orbit.load_npz(ROOT/cfg['linear_parent']/'modes_7_p14.npz'),old,7)
    transfer=out_transfer(seven,d); oldlin=orbit.load_npz(ref/'linear_reference.npz')
    K=d.K[np.ix_(d.out,d.out)];M=d.M0[np.ix_(d.out,d.out)]
    K7=seven.K[np.ix_(seven.out,seven.out)];M7=seven.M0[np.ix_(seven.out,seven.out)]
    checks=dict(K_transfer=linear.relative(transfer.T@K@transfer,K7),M_transfer=linear.relative(transfer.T@M@transfer,M7))
    checks['seven_source_K']=linear.relative(seven.K,oldlin['K'])
    checks['seven_source_M']=linear.relative(seven.M0,oldlin['M0'])
    Phi=transfer@oldlin['out_Phi'];omega=oldlin['out_omega']
    checks['modal_mass']=linear.relative(Phi.T@M@Phi,np.eye(81))
    checks['modal_stiffness']=linear.relative(Phi.T@K@Phi,np.diag(omega**2))
    physical=[]
    for a,b,Ps,Pt in zip(seven.arms,d.arms,seven.maps,d.maps):
        x=np.linspace(0,a.length,37)
        for f in ('v','Phi','psi'):
            A=a.basis_at(x)[f]@Ps[a.slices[f]][:,seven.out]
            B=b.basis_at(x)[f]@Pt[b.slices[f]][:,d.out]@transfer
            physical.append(linear.relative(A,B))
    checks['physical_profiles']=max(physical)
    for name,value in checks.items():linear.require(value,cfg['operator_relative_gate'],'common OUT '+name)
    E=block_diag(K,M); state_transfer=block_diag(transfer,transfer)
    oldphysical=orbit.load_npz(ref/'physical_Floquet_and_displacement.npz')
    reference=state_transfer@oldphysical['OUT_unit'];reference/=np.sqrt(reference@E@reference)
    # Saved full model is used with effective EA coefficients already in d.arms.
    action=variation.QuadraticAction(d,old.T4,old.V4)
    sixmodel,algebra=orbit.six_action(old)
    _,_,second=variation.second_variation(sixmodel)
    R,L=hb.energy_modal_maps(d,d.out,Phi,omega)
    return dict(d=d,seven=seven,action=action,reference=reference,E=E,transfer=transfer,
        state_transfer=state_transfer,Phi=Phi,omega=omega,R=R,L=L,checks=checks,
        algebra=algebra,second_variation=second,
        orbits={H:orbit.load_npz(base/f'HB{H}.npz') for H in (11,13)})


def canonical_blocks(blocks):
    """L_vq is B (historical key C); it is NOT the damping-like tangent C."""
    M,B,D,Mdot,Bdot=(blocks[k] for k in ('M','C','Q','Mdot','Cdot'))
    n=len(M);I=np.eye(n);Z=np.zeros_like(M)
    factor=cho_factor(M,lower=True,check_finite=False)
    solve=lambda x:cho_solve(factor,x,check_finite=False)
    C=Mdot+B-B.T;K=Bdot-D
    Av=np.block([[Z,I],[-solve(K),-solve(C)]])
    MB=solve(B);MI=solve(I)
    Ac=np.block([[-MB,MI],[D-B.T@MB,MB.T]])
    S=np.block([[I,Z],[B,M]]);Sdot=np.block([[Z,Z],[Bdot,Mdot]])
    J=np.block([[Z,I],[-I,Z]])
    return dict(M=M,B=B,D=D,C=C,K=K,Av=Av,Ac=Ac,S=S,Sdot=Sdot,J=J,
                mass_pivot=float(np.min(np.diag(factor[0]))))


def orbit_states(z,phases):
    return hb.hb_phases(z['U'],float(z['Omega']),z['Phi'],np.asarray(phases))


def canonical_fourier(problem,H,count,order,cfg,deadline):
    """Adapt the historical canonical_fourier to full supplied-acceleration checks."""
    d=problem['d'];action=problem['action'];R,L=problem['R'],problem['L'];z=problem['orbits'][H]
    phases=np.arange(count)*2*np.pi/count;q,v,a=orbit_states(z,phases)
    values=[];pivots=[];diagnostics=[];sampled={key:[] for key in ('M','B','D','C','K','S','Sdot')}
    scale0=d.M0[np.ix_(d.out,d.out)]
    for i,(qi,vi,ai) in enumerate(zip(q,v,a)):
        check_time(deadline)
        o=canonical_blocks(action.blocks(qi,vi,ai));values.append(L@o['Ac']@R);pivots.append(o['mass_pivot'])
        rec=dict(Hamiltonian=linear.relative(o['Ac'].T@o['J'],-o['J']@o['Ac']),
            symmetry_M=linear.relative(o['M'],o['M'].T),symmetry_D=linear.relative(o['D'],o['D'].T))
        if i in (0,7,31):
            Q,V,Ma=d.tangents(qi,vi,ai);ids=np.ix_(d.out,d.out)
            rec.update(massform_M=linear.relative(o['M'],Ma[ids]),massform_C=linear.relative(o['C'],V[ids]),massform_K=linear.relative(o['K'],Q[ids]),
                canonical_transform=linear.relative(np.linalg.solve(o['S'].T,(o['Sdot']+o['S']@o['Av']).T).T,o['Ac']),
                forbidden_PL_OUT=float(max(np.max(abs(X[np.ix_(d.out,d.plane)])) for X in (Q,V,Ma))))
            flow=np.linalg.solve(Ma,-d.inertial_terms(qi,vi)-d.potential(qi)['gradient'])
            residual=d.residual(qi,vi,ai);delta=ai-flow
            force_scale=np.linalg.norm(Ma@ai)+np.linalg.norm(d.potential(qi)['gradient'])
            rec['base_acceleration_identity']=float(np.linalg.norm(Ma@delta-residual)/max(force_scale,1e-30))
            Qflow,Vflow,Mflow=d.tangents(qi,vi,flow)
            flow_action=canonical_blocks(action.blocks(qi,vi,flow))
            rec['same_acceleration_flow_action']=linear.relative(flow_action['K'],Qflow[ids])
            rec['base_residual_norm']=float(np.linalg.norm(residual))
            rec['HB_minus_flow_acceleration_norm']=float(np.linalg.norm(delta))
            rec['HB_minus_flow_K_relative']=linear.relative(Q,Qflow)
            # Independent directional M_q*(a_HB-a_flow). No off-orbit equality is imposed.
            direction=np.zeros(d.ndof);direction[d.out]=np.sin(np.arange(81)+.2)
            h=2e-5;term=((d.mass_matrix(qi+h*direction)-d.mass_matrix(qi-h*direction))/(2*h))@delta
            rec['residual_acceleration_tangent_correction']=float(np.linalg.norm((Q-Qflow)@direction-term)/max(np.linalg.norm(Q@direction),1e-30))
        diagnostics.append(rec)
        for key in sampled:sampled[key].append(o[key])
    values=np.array(values);orders=np.arange(-order,order+1)
    F=np.einsum('mt,tij->mij',np.exp(-1j*orders[:,None]*phases)/count,values,optimize=True)
    fft=np.fft.fft(values,axis=0)/count;fftk=np.fft.fftfreq(count)*count
    tail=float(np.linalg.norm(fft[abs(fftk)>order])/np.linalg.norm(fft))
    independent=np.arange(cfg['independent_phases'])*2*np.pi/cfg['independent_phases']+cfg['independent_phase_offset']
    iq,iv,ia=orbit_states(z,independent);errors=[];imag=[];indpiv=[]
    for ph,qi,vi,ai in zip(independent,iq,iv,ia):
        check_time(deadline);o=canonical_blocks(action.blocks(qi,vi,ai));actual=L@o['Ac']@R
        rebuilt=np.einsum('m,mij->ij',np.exp(1j*orders*ph),F)
        errors.append(linear.relative(rebuilt,actual));imag.append(float(np.linalg.norm(rebuilt.imag)/np.linalg.norm(actual)));indpiv.append(o['mass_pivot'])
    pq,pv,pa=orbit_states(z,[0,2*np.pi]);end=[canonical_blocks(action.blocks(qi,vi,ai)) for qi,vi,ai in zip(pq,pv,pa)]
    periodic=linear.relative(end[0]['Ac'],end[1]['Ac'])
    maxima={k:max(r[k] for r in diagnostics if k in r) for k in diagnostics[0]}
    for k in diagnostics[7]:maxima[k]=max(r[k] for r in diagnostics if k in r)
    non_gates={'base_residual_norm','HB_minus_flow_acceleration_norm','HB_minus_flow_K_relative'}
    operator_pass=all(v<=cfg['operator_relative_gate'] for k,v in maxima.items() if k not in non_gates)
    duality=linear.relative(L@R,np.eye(162));conjugacy=linear.relative(F,F[::-1].conj())
    fourier_pass=max(max(errors),tail,max(imag),conjugacy,periodic)<=cfg['Fourier_relative_gate']
    diag=dict(HB=H,count=count,order=order,operator_checks=maxima,operator_pass=bool(operator_pass and duality<=cfg['operator_relative_gate']),
        Fourier_pass=bool(fourier_pass),independent_reconstruction=max(errors),discarded_sampled_Fourier_tail=tail,
        reconstruction_imaginary=max(imag),conjugacy=conjugacy,periodicity=periodic,modal_duality=duality,
        minimum_mass_Cholesky_pivot=min(pivots+indpiv),mass_change_from_rest=max(linear.relative(x,scale0) for x in sampled['M']),
        independent_phases=len(independent),basis_residual_is_not_zero=True)
    arrays=dict(F=F,orders=orders,R=R,L=L,phases=phases,Ac_energy=values,Omega=z['Omega'],
        independent_phases=independent,independent_errors=errors,**{k:np.array(v) for k,v in sampled.items()})
    return diag,arrays


def real_subspace(y,E,gate=1e-8):
    """Use a real line only when the weighted real/imaginary Gram supports it."""
    B=np.c_[y.real,y.imag];Gram=B.T@E@B
    vals,vecs=np.linalg.eigh(Gram);ratio=max(float(vals[0]),0)/max(float(vals[-1]),1e-300)
    chosen=np.flatnonzero(vals>max(vals[-1]*gate,1e-300))
    return B@vecs[:,chosen]/np.sqrt(vals[chosen]),ratio


def subspace_match(y,reference,E,gate):
    U,ratio=real_subspace(y,E,gate);ref=reference/np.sqrt(reference@E@reference)
    projection=U.T@E@ref;mac=float(np.sum(projection**2))
    aligned=U@projection/np.sqrt(max(mac,1e-300))
    return min(1.,mac),ratio,U.shape[1],aligned


def floquet_values(coeff,lam,Omega,phases):
    N=(len(coeff)-1)//2;n=np.arange(-N,N+1)
    periodic=np.exp(1j*np.asarray(phases)[:,None]*n)@coeff
    derivative=np.exp(1j*np.asarray(phases)[:,None]*n)@((lam+1j*n*Omega)[:,None]*coeff)
    return periodic,derivative


def direct_phase_operators(problem,H,cfg):
    phases=np.arange(cfg['Floquet_phases'])*2*np.pi/cfg['Floquet_phases']+cfg['Floquet_phase_offset']
    q,v,a=orbit_states(problem['orbits'][H],phases)
    ops=[canonical_blocks(problem['action'].blocks(qi,vi,ai)) for qi,vi,ai in zip(q,v,a)]
    return phases,ops


def spectral_task(problem,arrays,H,k,N,cfg,deadline):
    """One eig decomposition, using the frozen Hill assembly/convolution.

    Historical hill_spectrum clips sigma and discards alternative vectors;
    this adapter retains signed spectra, matched and maximum qualified modes.
    """
    check_time(deadline);start=time.perf_counter();F,orders=arrays['F'],arrays['orders'];Omega=float(arrays['Omega']);T=2*np.pi/Omega
    matrix,coeff=hb.build_hill(F,orders,Omega,N,k)
    eigen,vectors=eig(matrix,check_finite=False,overwrite_a=False)
    eig_seconds=time.perf_counter()-start;check_time(deadline)
    folded=hb.fold_exponents(eigen,Omega);states=vectors.reshape(2*N+1,2*k,-1)
    edge,aliases=hb.alias_edge_weights(eigen,folded,states,Omega);found=np.flatnonzero(edge<=cfg['Hill_edge_weight_gate'])
    ids=np.r_[np.arange(k),np.arange(k)+81];R=problem['R'][:,ids]
    q0,v0,a0=orbit_states(problem['orbits'][H],[0]);zero=canonical_blocks(problem['action'].blocks(q0[0],v0[0],a0[0]))
    physical=np.linalg.solve(zero['S'],R@states[:,:,found].sum(axis=0))
    matches=[subspace_match(y,problem['reference'],problem['E'],cfg['real_Gram_ratio_gate']) for y in physical.T]
    if not len(found):raise ArithmeticError('No interior Fourier replica; no mode assignment')
    # Similarity ranks physical subspaces; frequency proximity is not an acceptance gate.
    target=int(found[np.argmax([m[0] for m in matches])])
    phases,ops=direct_phase_operators(problem,H,cfg);details={};vectors_to_save={}
    matrixnorm=np.linalg.norm(matrix)
    def assess(index):
        if index in details:return details[index]
        check_time(deadline);lam=eigen[index];z=vectors[:,index];cu=z.reshape(2*N+1,2*k)
        hz=hb.hill_action(coeff,orders,Omega,N,z);dense=matrix@z
        local=float(np.linalg.norm(hz-lam*z)/max(np.linalg.norm(hz)+abs(lam)*np.linalg.norm(z),1e-300))
        conv=linear.relative(hz,dense)
        # Full physical dimension, direct action at independent phases, all omitted modes included.
        vals,dots=floquet_values(cu,lam,Omega,phases);pv=vals@R.T;pd=dots@R.T
        full_errors=[float(np.linalg.norm(di-o['Ac']@vi)/max(np.linalg.norm(di)+np.linalg.norm(o['Ac']@vi),1e-30)) for vi,di,o in zip(pv,pd,ops)]
        errors=[float(np.linalg.norm(di-R@(problem['L'][ids]@o['Ac']@vi))/max(np.linalg.norm(di)+np.linalg.norm(R@(problem['L'][ids]@o['Ac']@vi)),1e-30)) for vi,di,o in zip(pv,pd,ops)] if k<81 else full_errors
        y=np.linalg.solve(zero['S'],R@cu.sum(axis=0));mac,rank,dim,aligned=subspace_match(y,problem['reference'],problem['E'],cfg['real_Gram_ratio_gate'])
        quartet=max(float(np.min(abs(eigen-t))) for t in (lam.conjugate(),-lam,-lam.conjugate()))
        rho=np.exp(lam*T)
        qualified=local<=cfg['Hill_residual_gate'] and edge[index]<=cfg['Hill_edge_weight_gate'] and max(errors)<=cfg['Floquet_pointwise_relative_gate'] and quartet<=cfg['spectral_control_absolute_gate'] and conv<=cfg['operator_relative_gate']
        rec=dict(index=int(index),raw_lambda=cdict(lam),folded_lambda=cdict(folded[index]),alias=int(aliases[index]),sigma=float(lam.real),sigma_T=float(lam.real*T),rho=cdict(rho),
            edge_weight=float(edge[index]),local_residual=local,global_residual=float(np.linalg.norm(hz-lam*z)/matrixnorm),convolution_dense_relative=conv,
            quartet_defect=quartet,pointwise_residual=max(errors),full81_pointwise_residual=max(full_errors),physical_MAC=mac,real_Gram_ratio=rank,real_subspace_dimension=dim,
            qualified=bool(qualified),negative_real_multiplier=bool(rho.real<0 and abs(rho.imag)/abs(rho)<=cfg['real_multiplier_relative_gate']))
        details[index]=rec;vectors_to_save[index]=dict(vector=z,physical_complex=y,physical_aligned=aligned,pointwise_errors=errors)
        return rec
    matched=assess(target)
    maximum=None
    for index in found[np.argsort(eigen[found].real)[::-1]]:
        rec=assess(int(index))
        if rec['qualified']:maximum=rec;break
    reciprocal=int(found[np.argmin(abs(eigen[found]+eigen[target].conjugate()))]);stable=assess(reciprocal)
    competitor=[]
    for index,match in zip(found,matches):
        # Quotient by multiplier equivalence, so conjugate boundary copies are not competitors.
        if abs(np.exp(eigen[index]*T)-np.exp(eigen[target]*T))>cfg['replica_multiplier_tolerance']:
            competitor.append((match[0],int(index)))
    other=max(competitor,default=(0.,-1))
    identity=matched['physical_MAC']>=cfg['mode_identity_MAC_gate'] and matched['physical_MAC']-other[0]>=cfg['mode_identity_margin_gate']
    pair=abs(eigen[reciprocal]+eigen[target].conjugate())
    result=dict(HB=H,k=k,N=N,Omega=Omega,period=T,dimension=len(matrix),dense_MiB=matrix.nbytes/2**20,
        eigen_seconds=eig_seconds,seconds=time.perf_counter()-start,interior_replica_count=len(found),target=matched,
        maximum_qualified=maximum,maximum_is_target_mode=bool(maximum is not None and abs(uncomplex(maximum['rho'])-uncomplex(matched['rho']))<=cfg['replica_multiplier_tolerance']),
        reciprocal=stable,reciprocal_pair_defect=float(pair),reciprocal_multiplier_product_defect=float(abs(np.exp((eigen[target]+eigen[reciprocal])*T)-1)),
        competing_MAC=float(other[0]),competing_index=other[1],identity_established=bool(identity),
        assessed_candidates=list(details.values()),signed_real_parts_preserved=True)
    saved=dict(eigenvalues=eigen,folded_eigenvalues=folded,aliases=aliases,edge_weights=edge,accepted=edge<=cfg['Hill_edge_weight_gate'],
        candidate_indices=found,candidate_physical_MAC=np.array([m[0] for m in matches]),
        target_eigenvector=vectors_to_save[target]['vector'],reciprocal_eigenvector=vectors_to_save[reciprocal]['vector'],
        target_lambda=eigen[target],reciprocal_lambda=eigen[reciprocal],R_reduced=R,S0=zero['S'],metric=problem['E'],
        target_physical_complex=vectors_to_save[target]['physical_complex'],target_physical_aligned=vectors_to_save[target]['physical_aligned'],
        reciprocal_physical_complex=vectors_to_save[reciprocal]['physical_complex'],
        pointwise_phases=phases,target_pointwise_errors=vectors_to_save[target]['pointwise_errors'],reciprocal_pointwise_errors=vectors_to_save[reciprocal]['pointwise_errors'])
    if maximum is not None:saved['maximum_eigenvector']=vectors_to_save[maximum['index']]['vector']
    return result,saved


def physical_composition(problem,y):
    d=problem['d'];q=np.zeros(d.ndof);v=np.zeros(d.ndof);q[d.out]=y[:81];v[d.out]=y[81:]
    parts=orbit.CoupledSpatialGalerkin.reference_energy_parts(d,q,v)
    M=d.M0[np.ix_(d.out,d.out)];Phi=problem['Phi'];w=problem['omega']
    Q=Phi.T@M@y[:81];V=Phi.T@M@y[81:];weights=(w*Q)**2+V**2;weights/=weights.sum()
    profiles=np.array([a.reconstruct(P@q,points=np.linspace(0,a.length,201))[:,:6] for a,P in zip(d.arms,d.maps)])
    speeds=np.array([a.reconstruct(P@v,points=np.linspace(0,a.length,201))[:,:6] for a,P in zip(d.arms,d.maps)])
    return dict(bending_fraction=float(parts[:,0].sum()/parts.sum()),torsion_fraction=float(parts[:,1].sum()/parts.sum()),
        leading_OUT=[dict(position=int(j+1),fraction=float(weights[j])) for j in np.argsort(weights)[::-1][:6]]),dict(profiles=profiles,speeds=speeds,modal_weights=weights,reference_parts=parts)


def comparison(problem,records,solutions,operator_checks,cfg):
    ref=ROOT/cfg['seven_parent'];old=io.read(ref/'OUT_checks.json');oldrow=io.read(ref/'OUT_k81_N9.json')
    last=records[2];target=last['target'];sigma=target['sigma'];period=last['period']
    changes=dict(Hill_N=abs(sigma-records[1]['target']['sigma']),modal_72_81=abs(records[1]['target']['sigma']-records[0]['target']['sigma']),
        HB11_HB13=abs(sigma-records[3]['target']['sigma']))
    eta=max(*changes.values(),cfg['spectral_numerical_floor']);eta7=old['uncertainty_indicator']
    oldsigma=uncomplex(oldrow['dominant_raw_lambda']).real;T7=2*np.pi/io.read(ref/'HB13.json')['Omega'];rho7=np.exp(uncomplex(oldrow['dominant_raw_lambda'])*T7)
    controls=[];y=solutions[2]['target_physical_aligned'];E=problem['E']
    for sol in solutions:
        other=sol['target_physical_aligned'];controls.append(float((y@E@other)**2/((y@E@y)*(other@E@other))))
    identity=all(r['identity_established'] for r in records) and min(controls)>=cfg['mode_identity_MAC_gate']
    operators=all(r['operator_pass'] and r['Fourier_pass'] and r['minimum_mass_Cholesky_pivot']>0 for r in operator_checks.values())
    spectral=all(r['target']['qualified'] and r['reciprocal']['qualified'] and r['reciprocal_pair_defect']<=cfg['spectral_control_absolute_gate'] for r in records)
    resolved=operators and spectral and identity and max(changes.values())<=cfg['spectral_control_absolute_gate']
    positive=resolved and sigma>cfg['growth_uncertainty_factor']*eta
    comp6,data6=physical_composition(problem,y);comp7,data7=physical_composition(problem,problem['reference'])
    rho=uncomplex(target['rho']);growth6=sigma*period;growth7=oldsigma*T7
    rows=[dict(model=7,sigma=oldsigma,period=T7,sigma_T=growth7,rho_real=rho7.real,rho_imag=rho7.imag,eta=eta7,
               efolding_time=1/oldsigma,efolding_periods=1/growth7,amplification_5T=np.exp(5*growth7),**comp7),
          dict(model=6,sigma=sigma,period=period,sigma_T=growth6,rho_real=rho.real,rho_imag=rho.imag,eta=eta,
               efolding_time=1/sigma if positive else None,efolding_periods=1/growth6 if positive else None,
               amplification_5T=np.exp(5*growth6),**comp6)]
    result=dict(rows=rows,operator_status='PASS' if operators else 'FAIL',mode_qualification='PASS' if resolved else 'UNRESOLVED',
        growth_status='RESOLVED_POSITIVE' if positive else ('NO_RESOLVED_POSITIVE_TARGET_GROWTH' if resolved else 'UNRESOLVED'),
        correspondence='ESTABLISHED' if identity else 'UNRESOLVED',resolution_changes=changes,control_physical_MAC=controls,
        eta6=eta,eta7=eta7,eta_diff=eta+eta7,eta7_definition='max(saved last72/all81 modal difference, N7/N9, HB11/HB13 k24 difference, 1e-8); p-error excluded',
        sigma_difference=sigma-oldsigma,sigma_relative_difference=(sigma-oldsigma)/oldsigma,
        sigma_T_difference=growth6-growth7,sigma_T_relative_difference=(growth6-growth7)/growth7,
        rho_difference=cdict(rho-rho7),rho_relative_difference=cdict((rho-rho7)/rho7),
        difference_to_observed_indicator=abs(sigma-oldsigma)/(eta+eta7),growth_to_indicator=sigma/eta,
        physical_MAC=target['physical_MAC'],negative_real_multiplier=target['negative_real_multiplier'],
        mode_composition_metric='phase-zero unit E0, common physical velocity and displacement; quadratic reference energy only',
        no_initial_projection=True,no_ODE=True,p_convergence='NOT_TESTED')
    data=dict(unit6=y,unit7=problem['reference'],metric=E,
        **{k+'_6':v for k,v in data6.items()},**{k+'_7':v for k,v in data7.items()})
    # Full complex coefficients remain available even when a real line is qualified.
    N=9;lam=complex(solutions[2]['target_lambda']);z=solutions[2]['target_eigenvector'].reshape(19,162);R=problem['R']
    phases=np.linspace(0,2*np.pi,65);vals,_=floquet_values(z,lam,float(problem['orbits'][13]['Omega']),phases)
    q,v,a=orbit_states(problem['orbits'][13],phases);physical=[]
    for ph,val,qi,vi,ai in zip(phases,vals,q,v,a):
        S=canonical_blocks(problem['action'].blocks(qi,vi,ai))['S']
        physical.append(np.linalg.solve(S,R@val)*np.exp(1j*lam.imag*ph/float(problem['orbits'][13]['Omega'])))
    physical=np.array(physical)
    result['growth_removed_antiperiodic_relative_defect']=linear.relative(physical[-1],-physical[0])
    data.update(phases=phases,growth_removed_complex_physical=physical)
    return result,data


def render_saved(bundle):
    """No operator assembly or eigenproblem in reports or plots."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    b=Path(bundle);s=io.cache(b)
    if 'comparison' not in s:return s
    c=s['comparison'];z=orbit.load_npz(b/'comparison.npz')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    rows=c['rows'];axes[0].errorbar([7,6],[r['sigma'] for r in rows],yerr=[r['eta'] for r in rows],fmt='o',capsize=5)
    axes[0].set_xticks([7,6],['7 fields','6 fields']);axes[0].set_ylabel('sigma (inverse model time)')
    periods=np.linspace(0,5,151)
    for r,style in zip(rows,('-','--')):axes[1].plot(periods,np.exp(r['sigma_T']*periods),style,label=f"{r['model']} fields")
    axes[1].set_xlabel('own orbit periods');axes[1].set_ylabel('single Floquet component amplification');axes[1].legend()
    for ax in axes:ax.grid(alpha=.25)
    fig.suptitle('Spectral prediction; observed indicators, no time integration')
    fig.tight_layout();fig.savefig(b/'growth_comparison.png',dpi=180);fig.savefig(b/'growth_comparison.pdf');plt.close(fig)
    fig,axes=plt.subplots(3,2,figsize=(10,8),sharex=True)
    for arm in range(2):
        for row,(index,label) in enumerate(((2,'v'),(3,'Phi'),(4,'psi'))):
            for model,style in ((7,'-'),(6,'--')):axes[row,arm].plot(np.linspace(0,1,201),z[f'profiles_{model}'][arm,:,index],style,label=f'{model} fields')
            axes[row,arm].set_ylabel(label);axes[row,arm].grid(alpha=.25)
        axes[0,arm].set_title(f'Arm {arm+1}');axes[2,arm].set_xlabel('s / L (clamp to joint)')
    axes[0,0].legend();fig.suptitle('Matched physical mode at tau=0; common unit E0, sign aligned')
    fig.tight_layout();fig.savefig(b/'physical_mode_comparison.png',dpi=180);fig.savefig(b/'physical_mode_comparison.pdf');plt.close(fig)
    lines=['# Six-field OUT stability','',s['overall'],'','|model|sigma|sigma T|rho real|eta|','|---:|---:|---:|---:|---:|']
    for r in rows:lines.append('|'+ '|'.join(str(r[k]) for k in ('model','sigma','sigma_T','rho_real','eta'))+'|')
    (b/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf8');io.finalize(b,io.read(b/'provenance.json'),s)
    return s


def run():
    cfg=io.read(CONFIG)
    if cfg['tasks']!=[[13,72,7],[13,81,7],[13,81,9],[11,81,9]] or cfg['degree']!=14 or cfg['A']!=.005 or cfg['numerical_budget_seconds']>1800:
        raise ValueError('Only the four explicitly authorized OUT tasks are supported')
    selected,manifests=sources(cfg)
    inputs=dict(schema=VERSION,config=cfg,source_artifacts=selected,source_manifests=manifests,code_hashes={p:io.sha(ROOT/p) for p in CODE})
    fp=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()[:16];b=OUTPUT/fp
    if (b/'manifest.json').exists():return b,io.cache(b),True
    if b.exists() or (OUTPUT.exists() and list(OUTPUT.glob('*/provenance.json'))):raise ValueError('Preserve prior attempt; no automatic repetition of spectral tasks')
    b.mkdir(parents=True);initial=io.read(ROOT/'results/_smoke/nlsp_six_vs_seven_out_stability/initial_state.json')
    prov=dict(inputs,HEAD=initial['HEAD'],initial_state=initial,interpreter=sys.executable,started_UTC=datetime.now(timezone.utc).isoformat())
    io.write(b/'config.json',cfg);io.write(b/'provenance.json',prov)
    for name in CODE:
        dest=b/'execution_code'/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/name,dest)
    s=dict(overall='PREPARED',completed=False,calls=dict(Hill=0,**cfg['forbidden_calls']),operator_checks={},tasks=[])
    io.write(b/'cost_forecast.json',dict(total_seconds=cfg['forecast_seconds'],budget_seconds=1800,
        task_dimensions=[2160,2430,3078,3078],largest_dense_MiB=3078**2*16/2**20,memory_reserve_MiB=cfg['memory_reserve_MiB'],
        basis='saved seven-field p14 N7 18.36s, N9 37.40s; no benchmark eigensolve'))
    preflight=io.read(ROOT/'results/_smoke/nlsp_six_vs_seven_out_stability/preflight.json')
    io.write(b/'preflight.json',preflight)
    start=time.perf_counter();deadline=start+cfg['numerical_budget_seconds']-preflight['seconds']
    def save():
        s['numerical_seconds']=time.perf_counter()-start+preflight['seconds'];io.finalize(b,prov,s)
    save()
    try:
        problem=load_problem(cfg);s['basis_checks']=problem['checks']
        io.write(b/'second_variation_checks.json',problem['second_variation'])
        np.savez_compressed(b/'recovery.npz',transfer_7_to_6=problem['transfer'],state_transfer=problem['state_transfer'],
            Phi=problem['Phi'],omega=problem['omega'],R=problem['R'],L=problem['L'],metric=problem['E'],reference7=problem['reference'],
            W=problem['d'].W,raw_maps=problem['d'].raw_maps,local_maps=np.array(problem['d'].maps),
            plane=problem['d'].plane,out=problem['d'].out,labels=problem['d'].saved['labels'],fields=problem['d'].saved['fields'],lengths=problem['d'].lengths,bases=problem['d'].bases)
        fourier={}
        # Both orbital operators must pass before ANY spectral solve.
        for H in (13,11):
            attempts=[]
            for count,order in (cfg['Fourier_initial'],cfg['Fourier_single_refinement']):
                diag,arrays=canonical_fourier(problem,H,count,order,cfg,deadline);attempts.append(diag)
                io.write(b/f'HB{H}_Fourier_{count}_{order}.json',diag)
                np.savez_compressed(b/f'HB{H}_Fourier_{count}_{order}.npz',**arrays)
                s['operator_checks'][str(H)]=diag;save()
                if not diag['operator_pass']:raise ArithmeticError('OUT variational/canonical operator gate failed')
                if diag['Fourier_pass']:break
                if count==cfg['Fourier_initial'][0]:io.write(b/f'HB{H}_Fourier_refinement_reason.json',dict(reason='Initial Fourier representation gate failed',initial=diag))
            if not diag['Fourier_pass']:raise ArithmeticError('Fourier representation unresolved at bounded refinement')
            fourier[H]=arrays
            print(json.dumps(dict(stage='operator',HB=H,count=count,order=order,reconstruction=diag['independent_reconstruction'])),flush=True)
        solutions=[]
        for H,k,N in cfg['tasks']:
            check_time(deadline)
            if deadline-time.perf_counter()<75:raise TimeoutError('Insufficient budget reserve before next eig decomposition')
            s['calls']['Hill']+=1;save()
            record,arrays=spectral_task(problem,fourier[H],H,k,N,cfg,deadline)
            name=f'HB{H}_OUT{k}_N{N}';io.write(b/(name+'.json'),record);np.savez_compressed(b/(name+'.npz'),**arrays)
            s['tasks'].append(record);solutions.append(arrays);save()
            print(json.dumps(dict(stage=name,sigma=record['target']['sigma'],MAC=record['target']['physical_MAC'],qualified=record['target']['qualified'],seconds=record['seconds'])),flush=True)
        c,data=comparison(problem,s['tasks'],solutions,s['operator_checks'],cfg);s['comparison']=c
        io.write(b/'comparison.json',c);np.savez_compressed(b/'comparison.npz',**data)
        linear.write_csv(b/'comparison.csv',[{k:v for k,v in r.items() if k!='leading_OUT'} for r in c['rows']])
        linear.write_csv(b/'controls.csv',[dict(HB=r['HB'],k=r['k'],N=r['N'],**{k:v for k,v in r['target'].items() if not isinstance(v,dict)}) for r in s['tasks']])
        s['overall']='SIX_FIELD_OUT_GROWTH_REPRODUCED' if c['growth_status']=='RESOLVED_POSITIVE' else 'SIX_FIELD_OUT_BOUNDED_RESULT_'+c['mode_qualification']
        s['completed']=True;save();render_saved(b)
    except (ArithmeticError,ValueError,TimeoutError,np.linalg.LinAlgError) as exc:
        s.update(overall='BOUNDED_OUT_ATTEMPT_STOPPED',completed=True,failure=repr(exc));save();raise
    return b,s,False


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);m=parser.add_mutually_exclusive_group(required=True)
    m.add_argument('--compute',action='store_true');m.add_argument('--report-only',type=Path);m.add_argument('--plot-only',type=Path)
    args=parser.parse_args(argv)
    if args.compute:b,s,cached=run()
    else:b=(args.report_only or args.plot_only).resolve();s=render_saved(b);cached=True
    print(json.dumps(io.plain(dict(bundle=str(b),overall=s['overall'],reused=cached,new_Hill_solves=0 if cached else s['calls']['Hill'],comparison=s.get('comparison'))),indent=2))
