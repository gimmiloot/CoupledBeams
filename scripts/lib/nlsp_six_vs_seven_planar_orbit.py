"""Six-field quartic NLSP action and one full p14 reversible planar orbit.

The unchanged seven-field *local evaluator* is pulled back through c=0,
with its C coefficient changed to EA at every amplitude degree. The actual
discrete space has 162 coordinates, no c coordinate and no singular dummy mass.
No linear eigensolve, IVP, Hill or Floquet calculation is provided here.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np
from scipy.spatial.transform import Rotation

from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib import nlsp_coupled_sensitivity_pilot as io
from scripts.lib import nlsp_six_vs_seven_linear_comparison as linear
from scripts.lib import nlsp_planar_1to3_internal_resonance as hb
from scripts.lib.weakly_nonlinear_coupled_spatial_dynamics import EndpointArm, CoupledSpatialGalerkin

ROOT=io.ROOT
CONFIG=ROOT/'data/input/nlsp_six_vs_seven_planar_orbit.json'
OUTPUT=ROOT/'results/nlsp_six_vs_seven_planar_orbit'
MODEL='nlsp-six-field-EA-J0-quartic-action-v1'
VERSION='six-seven-planar-orbit-v1'
CODE=linear.KERNELS+[
    'scripts/lib/nlsp_six_vs_seven_planar_orbit.py',
    'scripts/lib/nlsp_planar_1to3_internal_resonance.py',
    'scripts/lib/nlsp_coupled_sensitivity_pilot.py',
    'scripts/lib/nlsp_six_vs_seven_linear_comparison.py',
    'scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py',
    'data/input/nlsp_six_vs_seven_planar_orbit.json']


def load_npz(path):
    with np.load(path) as z:
        return {k:z[k].copy() for k in z.files}


def restrict(poly,fields):
    ids={i for i,name in enumerate(rod.SYMBOL_ORDER[:42]) if name.split('_')[0] in fields}
    return rod.Polynomial({key:value for key,value in poly.terms.items() if not ids.intersection(key)})


def load_action(cfg):
    polynomials=io.read(ROOT/cfg['action_parent']/'result.json')['polynomials']
    return SimpleNamespace(T4=rod.Polynomial.deserialize(polynomials['T4']),
        V4=rod.Polynomial.deserialize(polynomials['V4']),
        residual_a=tuple(rod.Polynomial.deserialize(p) for p in polynomials['residuals_A']))


def six_action(old):
    """Exact rational audit: independent exponential series vs saved correction."""
    P=rod.Polynomial; p={name:P.symbol(name) for name in rod.SYMBOL_ORDER}
    a=(p['Phi'],-p['psi'],p['theta']); hat=rod._hat(a)
    R=rod._identity(); power=rod._identity()
    for degree,divisor in ((1,1),(2,2),(3,6)):
        power=rod._matmul(power,hat)
        R=rod._matrix_add(R,rod._matrix_scale(Fraction(1,divisor),power))
    RT=rod._transpose(R); base=(1+p['u_s'],p['w_s'],p['v_s'])
    gamma=rod._truncate_vector(rod._add(rod._matvec(RT,base),(-P(1),P(),P())),3)
    measures=[]
    for axis in ('s','t'):
        derivative=tuple(tuple(x.total_derivative(axis) for x in row) for row in R)
        measures.append(rod._truncate_vector(rod._vee(rod._matmul(RT,derivative)),3))
    chi,omega=measures; K=p['C']*(1-p['nu']**2)
    T=(p['m']*sum((p[f+'_t']**2 for f in ('u','w','v')),P())+
       sum((j*x*x for j,x in zip((p['jp']+p['jb'],p['jb'],p['jp']),omega)),P()))/2
    V=(sum((d*x*x for d,x in zip((K,p['S'],p['S']),gamma)),P())+
       sum((d*x*x for d,x in zip((p['CT'],p['Bb'],p['Bp']),chi)),P()))/2
    oldT,oldV=restrict(old.T4,{'c'}),restrict(old.V4,{'c'})
    correction=p['C']*p['nu']**2*gamma[0]**2/2
    if T!=oldT or V!=oldV-correction:
        raise ArithmeticError('Six-field action disagrees with independent matrix exponential')
    # Independent user-specified planar expression; K is still C*(1-nu^2).
    us,ws,th=p['u_s'],p['w_s'],p['theta']; S=p['S']
    expectedT=(p['m']*(p['u_t']**2+p['w_t']**2)+p['jp']*p['theta_t']**2)/2
    expectedV=(K*us**2+S*(ws-th)**2+p['Bp']*p['theta_s']**2)/2
    expectedV+=(K-S)*us*ws*th-(K-2*S)*us*th**2/2
    expectedV+=(K-S)*th**2*(ws**2-us**2)/2-(3*K-4*S)*ws*th**3/6+(3*K-4*S)*th**4/24
    if restrict(T,{'v','Phi','psi'})!=expectedT or restrict(V,{'v','Phi','psi'})!=expectedV:
        raise ArithmeticError('Independent planar energy formula mismatch')
    if T.substitute({'nu':0})!=oldT.substitute({'nu':0}) or V.substitute({'nu':0})!=oldV.substitute({'nu':0}):
        raise ArithmeticError('nu=0 action restriction mismatch')
    L=T-V
    residual=tuple(L.derivative(f+'_t').total_derivative('t')+
        L.derivative(f+'_s').total_derivative('s')-L.derivative(f) for f in rod.FIELD_ORDER[:6])
    if any(actual!=restrict(reference,{'c'}).substitute({'C':K})
           for actual,reference in zip(residual,old.residual_a[:6])):
        raise ArithmeticError('Restriction and variation disagree for the six physical fields')
    return SimpleNamespace(T4=T,V4=V,residual_a=residual),dict(
        exact_action_identity=True,exact_planar_identity=True,nu0_restricted_action_identity=True,
        correction_terms_by_degree={str(d):len(correction.homogeneous(d).terms) for d in (2,3,4)},
        residual_max_degree=max(rod._degree(k) for r in residual for k in r.terms),
        formal_proof='exact rational polynomial comparison; Lean unavailable on PATH')


class SavedActionSystem:
    """Variational pullback using saved physical maps and mass coordinate change."""
    local=CoupledSpatialGalerkin.local
    potential=CoupledSpatialGalerkin.potential
    inertial_terms=CoupledSpatialGalerkin.inertial_terms
    kinetic=CoupledSpatialGalerkin.kinetic
    energy=CoupledSpatialGalerkin.energy

    def __init__(self,saved,model,nf):
        self.saved=saved; self.nf=nf; self.p=int(saved['p'])
        if self.p!=14 or nf not in (6,7):
            raise ValueError('Only saved six/seven-field p14 allowed')
        self.W=saved['W']; self.ndof=len(self.W)
        self.plane=saved['plane']; self.out=saved['out']; self.lengths=tuple(saved['lengths'])
        self.bases=saved['bases']; self.state_plane=np.r_[self.plane,self.ndof+self.plane]
        self.state_out=np.r_[self.out,self.ndof+self.out]
        self.arms=[]; self.maps=[]; self.raw_maps=saved['raw_maps']
        self.coefficients=[]
        for i,L in enumerate(self.lengths):
            c=rod.RodCoefficients(**dict(zip(saved['coefficient_names'],saved['coefficients'][i])))
            self.coefficients.append(c)
            effective=replace(c,C=c.C*(1-c.nu**2)) if nf==6 else c
            self.arms.append(EndpointArm(effective,self.p,L,model))
            # The only internal embedding: named raw fields, BEFORE whitening.
            P=np.zeros((7*self.p,self.ndof)); P[:nf*self.p]=saved['raw_maps'][i]@self.W
            self.maps.append(P)
        self.M0=self.mass_matrix(np.zeros(self.ndof))
        self.K=self.potential(np.zeros(self.ndof),hessian=True)['hessian']
        linear.require(linear.relative(self.M0,saved['M0']),2e-10,'saved resting M')
        linear.require(linear.relative(self.K,saved['K']),2e-10,'saved linear K')
        if nf==6 and (self.ndof,len(self.plane),len(self.out))!=(162,81,81):
            raise ArithmeticError('Six-field p14 dimension mismatch')

    def mass_matrix(self,q):
        return sum(P.T@a.mass_matrix(qi)@P for a,P,qi in zip(self.arms,self.maps,self.local(q)))

    def residual(self,q,velocity,acceleration):
        return self.mass_matrix(q)@acceleration+self.inertial_terms(q,velocity)+self.potential(q)['gradient']

    def tangents(self,q,velocity,acceleration):
        """R_q, R_v, R_a at the supplied triple, including M_q*acceleration."""
        return hb.full_massform_tangents(self,q,velocity,acceleration)

    def energy_rate(self,q,velocity,acceleration):
        return sum(a.energy_rate(qi,vi,ai) for a,qi,vi,ai in
                   zip(self.arms,self.local(q),self.local(velocity),self.local(acceleration)))

    def joint_diagnostics(self,q,velocity=None):
        ends=[a.reconstruct(qi,np.array([a.length]))[0] for a,qi in zip(self.arms,self.local(q))]
        positions=[B@e[:3] for B,e in zip(self.bases,ends)]
        orientations=[B@Rotation.from_rotvec(e[3:6]*[1,-1,1]).as_matrix()@B.T for B,e in zip(self.bases,ends)]
        return dict(translation_difference=float(np.max(abs(positions[0]-positions[1]))),
            orientation_difference=float(np.max(abs(orientations[0]-orientations[1]))),
            contraction_difference=None if self.nf==6 else float(abs(ends[0][6]-ends[1][6])),
            contraction_status='NOT_APPLICABLE: no c coordinate or c constraint' if self.nf==6 else 'seven-field c continuity')

    def profiles(self,coordinates,points):
        return np.array([a.reconstruct_series(coordinates@P.T,points=points*L)[:,:,:self.nf]
                         for a,P,L in zip(self.arms,self.maps,self.lengths)])


def common_transfer(source,target):
    """Physical Galerkin transfer for exactly u,w,theta; no coefficient deletion."""
    cross=np.zeros((target.ndof,source.ndof)); gram=np.zeros((target.ndof,target.ndof))
    for a,b,Ps,Pt in zip(source.arms,target.arms,source.maps,target.maps):
        xi,w=np.polynomial.legendre.leggauss(41); x=(xi+1)*a.length/2; w=w*a.length/2
        for f,m in (('u',b.coefficients.m),('w',b.coefficients.m),('theta',b.coefficients.jp)):
            A=a.basis_at(x)[f]@Ps[a.slices[f]]; B=b.basis_at(x)[f]@Pt[b.slices[f]]
            cross+=m*B.T@(w[:,None]*A); gram+=m*B.T@(w[:,None]*B)
    transfer=np.zeros_like(cross); ids=target.plane
    transfer[ids]=np.linalg.solve(gram[np.ix_(ids,ids)],cross[ids])
    return transfer


def operator_checks(d,model6,old,cfg):
    out=dict(dimensions=[d.ndof,len(d.plane),len(d.out)],nq=[a.nq for a in d.arms],
        maximum_integrand_degree=4*d.p,gauss_exact_degree=2*d.arms[0].nq-1,
        K_saved_relative=linear.relative(d.K,d.saved['K']),M_saved_relative=linear.relative(d.M0,d.saved['M0']))
    rng=np.random.default_rng(20261021); records=[]
    for planar,amplitude in ((True,2e-3),(False,1e-3),(False,4e-3)):
        def state():
            raw=rng.normal(size=d.ndof)
            if planar:
                raw[d.out]=0
            # Bound local physical components while retaining all coordinate directions.
            q=np.linalg.solve(d.W,raw)
            scale=max(np.max(abs(a.reconstruct(P@q))) for a,P in zip(d.arms,d.maps))
            return q*amplitude/scale
        q,v,a=state(),state(),state(); dq,dv,da=state(),state(),state()
        F=d.potential(q,hessian=True); M=d.mass_matrix(q); Q,V,Ma=d.tangents(q,v,a)
        h=cfg['directional_step']
        vp=d.potential(q+h*dq,gradient=False)['V']; vm=d.potential(q-h*dq,gradient=False)['V']
        fdF=(d.potential(q+h*dq)['gradient']-d.potential(q-h*dq)['gradient'])/(2*h)
        fdR=(d.residual(q+h*dq,v+h*dv,a+h*da)-d.residual(q-h*dq,v-h*dv,a-h*da))/(2*h)
        expected=Q@dq+V@dv+Ma@da
        # Arbitrary acceleration is used throughout; no equilibrium substitution.
        fdE=(d.energy(q+h*v,v+h*a)-d.energy(q-h*v,v-h*a))/(2*h)
        rate=d.energy_rate(q,v,a); work=v@d.residual(q,v,a)
        rec=dict(planar=planar,physical_state_scale=amplitude,
            energy_direction=abs((vp-vm)/(2*h)-F['gradient']@dq)/max(abs(F['gradient']@dq),1e-30),
            force_direction=linear.relative(fdF,F['hessian']@dq),
            residual_direction=linear.relative(fdR,expected),
            kinetic_energy_rate_FD=abs(fdE-rate)/max(abs(rate),1e-30),
            energy_power_identity=abs(rate-work)/max(abs(rate),abs(work),1e-30),
            M_symmetry=linear.relative(M,M.T),V_hessian_symmetry=linear.relative(F['hessian'],F['hessian'].T),
            minimum_mass_Cholesky_pivot=float(np.min(np.diag(np.linalg.cholesky(M)))))
        # Verify derivatives of the tangents themselves (third variation of action).
        # First directional tangent is the reliable FD check; second differences
        # at a separate larger h avoid cancellation of the small nonlinear part.
        hh=.02
        Qp,Vp,Mp=d.tangents(q+hh*dq,v+hh*dv,a+hh*da)
        Qm,Vm,Mm=d.tangents(q-hh*dq,v-hh*dv,a-hh*da)
        second=((Qp-Qm)@dq+(Vp-Vm)@dv+(Mp-Mm)@da)/(2*hh)
        fd2=(d.residual(q+hh*dq,v+hh*dv,a+hh*da)-2*d.residual(q,v,a)+d.residual(q-hh*dq,v-hh*dv,a-hh*da))/hh**2
        rec['tangent_direction_second_variation']=linear.relative(second,fd2)
        if planar:
            rec['planar_mass_constant']=linear.relative(M[np.ix_(d.plane,d.plane)],d.M0[np.ix_(d.plane,d.plane)])
            rec['planar_inertia']=float(np.max(abs(d.inertial_terms(q,v))))
            rec['OUT_residual']=float(np.max(abs(d.residual(q,v,a)[d.out])))
            rec['cross_tangents']=max(float(np.max(abs(A[np.ix_(d.plane,d.out)]))) for A in (Q,V,M))
            for key in ('planar_mass_constant','planar_inertia','OUT_residual','cross_tangents'):
                linear.require(rec[key],cfg['operator_relative_gate'],key)
        else:
            rec['spatial_inertia_norm']=float(np.linalg.norm(d.inertial_terms(q,v)))
            rec['spatial_mass_change']=linear.relative(M,d.M0)
        for key in ('energy_direction','force_direction','residual_direction','kinetic_energy_rate_FD','tangent_direction_second_variation'):
            linear.require(rec[key],cfg['directional_derivative_relative_gate'],key)
        for key in ('M_symmetry','V_hessian_symmetry'):
            linear.require(rec[key],cfg['symmetry_relative_gate'],key)
        linear.require(rec['energy_power_identity'],cfg['operator_relative_gate'],'energy power')
        records.append(rec)
    # Numerical production action vs independently generated six-field polynomial.
    errors=[]
    for arm,c in zip(d.arms,d.coefficients):
        values={name:float(rng.normal()*2e-3) for name in rod.SYMBOL_ORDER[:42]}
        values.update(c.values()); production=dict(values); production['C']=c.C*(1-c.nu**2)
        for name in rod.SYMBOL_ORDER[:42]:
            if name.split('_')[0]=='c':production[name]=0.; values[name]=0.
        for actual,reference in ((old.T4,model6.T4),(old.V4,model6.V4)):
            errors.append(abs(actual.evaluate(production)-reference.evaluate(values))/max(abs(reference.evaluate(values)),1e-30))
    out['production_vs_independent_action']=max(errors); out['directional_checks']=records
    linear.require(max(errors),cfg['operator_relative_gate'],'production six action')
    return out


class SixPlanarHB(hb.FullPlanarHB):
    """Reuse unchanged equations/Newton; replace seven-field validation only."""
    def validate(self,U,Omega,A):
        cfg=self.cfg; d=self.d; H=len(U)-1
        phases=np.arange(cfg['independent_phase_count'])*2*np.pi/cfg['independent_phase_count']+cfg['independent_phase_offset']
        q,v,acc,jerk=hb.orbit_with_jerk(U,Omega,self.Phi,phases)
        res=[]; energies=[]; pivots=[]; joints=[]; analytic=[]; derivative_scales=[]
        for qi,vi,ai,ji in zip(q,v,acc,jerk):
            self.check_time(); M=d.mass_matrix(qi); F=d.potential(qi,hessian=True)
            res.append(M@ai+d.inertial_terms(qi,vi)+F['gradient']); energies.append(d.energy(qi,vi))
            pivots.append(np.min(np.diag(np.linalg.cholesky(M))))
            jt=d.joint_diagnostics(qi,vi)
            joints.append(max(jt['translation_difference'],jt['orientation_difference']))
            analytic.append(M@ji+F['hessian']@vi)
            derivative_scales.append(np.linalg.norm(M@ji)+np.linalg.norm(F['hessian']@vi))
        res=np.array(res); analytic=np.array(analytic); energies=np.array(energies)
        orders=np.arange(-3*H,3*H+1)
        coefficients=np.exp(-1j*orders[:,None]*phases)@res/len(phases)
        reconstructed=np.exp(1j*phases[:,None]*orders)@coefficients
        independent_dot=np.exp(1j*phases[:,None]*orders)@(1j*orders[:,None]*Omega*coefficients)
        scale=Omega**2*np.sqrt(U[1,1]**2+81*U[3,3]**2)
        amp,_,peak=hb.physical_amplitude_gradient(d,self.Phi,U)
        residual_max=float(np.max(np.linalg.norm(res,axis=1)))
        # The floor is explicit and common to the physical identity checks.
        denom=np.maximum(derivative_scales,cfg['zero_scale_floor'])
        identity=float(np.max(np.linalg.norm(independent_dot-analytic,axis=1)/denom))
        phase=float(np.angle(complex(U[3,3]))-3*np.angle(complex(U[1,1])))
        phase=float(np.angle(np.exp(1j*phase)))
        row=dict(A=A,Omega=Omega,period=2*np.pi/Omega,actual_A=amp,peak=peak,
            Q21=float(U[1,1]),Q43=float(U[3,3]),ratio=float(U[3,3]/U[1,1]),phase_Psi=phase,
            full_residual_absolute=residual_max,full_residual_relative=residual_max/scale,
            full_residual_RMS_relative=float(np.linalg.norm(res)/np.sqrt(len(phases))/scale),
            residual_scale=scale,energy_relative_variation=float(np.ptp(energies)/abs(np.mean(energies))),
            minimum_mass_Cholesky_pivot=float(min(pivots)),joint_compatibility=float(max(joints)),
            contraction_difference=None,contraction_status='NOT_APPLICABLE: six-field action has no c',
            residual_time_derivative_max=float(np.max(np.linalg.norm(independent_dot,axis=1))),
            residual_time_derivative_RMS=float(np.linalg.norm(independent_dot)/np.sqrt(len(phases))),
            time_shift_identity_relative_defect=identity,
            residual_Fourier_reconstruction_absolute=float(np.max(abs(reconstructed-res))),
            independent_phase_count=len(phases),residual_harmonic_limit=3*H,
            harmonic_common_mass_norm=np.linalg.norm(U,axis=1).tolist(),
            harmonic_reference_position_energy_norm=np.linalg.norm(U*self.omega,axis=1).tolist(),
            last_two_harmonic_energy_norm_ratio=float(np.linalg.norm(U[-2:]*self.omega)/np.linalg.norm(U*self.omega)),
            initial_velocity_max=float(np.max(abs(hb.hb_phases(U,Omega,self.Phi,np.array([0.]))[1]))))
        linear.require(identity,cfg['operator_relative_gate'],'time shift residual derivative identity')
        return row,dict(U=U,Omega=Omega,Phi=self.Phi,phases=phases,q=q,v=v,acceleration=acc,jerk=jerk,
            massform_residual=res,energy=energies,residual_orders=orders,residual_coefficients=coefficients,
            residual_time_derivative=independent_dot.real,time_shift_identity_defect=independent_dot.real-analytic)


def acceptance(row,cfg):
    gates=dict(projected=row['Newton']['scaled_projected_residual']<=cfg['HB_projected_relative_target'],
        full_residual=row['full_residual_relative']<=cfg['HB_full_relative_gate'],
        energy=row['energy_relative_variation']<=cfg['energy_relative_gate'],
        amplitude=abs(row['actual_A']-cfg['A'])<=cfg['amplitude_absolute_gate'],
        joint=row['joint_compatibility']<=cfg['joint_absolute_gate'],
        positive_mass=row['minimum_mass_Cholesky_pivot']>0,
        phase_zero=abs(row['phase_Psi'])<1e-12 and row['Q21']>0 and row['Q43']>0,
        time_shift_identity=row['time_shift_identity_relative_defect']<=cfg['operator_relative_gate'])
    return dict(gates=gates,passed=all(gates.values()))


def verified_sources(cfg):
    selected={}; manifests={}
    files={
        'linear_parent':['config.json','provenance.json','modes_6_p14.npz','modes_7_p14.npz','matching_p14.npz'],
        'orbit_parent':['config.json','provenance.json','linear_reference.npz','HB11.npz','HB13.npz','HB11.json','HB13.json'],
        'action_parent':['result.json']}
    for key,names in files.items():
        folder=ROOT/cfg[key]; manifest=io.read(folder/'manifest.json')
        manifests[key]=dict(path=cfg[key],sha256=io.sha(folder/'manifest.json'))
        hashes=manifest.get('artifact_hashes',manifest.get('artifacts',{}))
        for name in names:
            digest=io.sha(folder/name)
            if hashes.get(name)!=digest:
                raise ValueError('Source artifact mismatch: '+str(folder/name))
            selected[(folder/name).relative_to(ROOT).as_posix()]=digest
    for key in ('linear_parent','orbit_parent'):
        old=io.read(ROOT/cfg[key]/'config.json')
        for name in ('lengths','sections','beta_deg','material'):
            if old[name]!=cfg[name]:
                raise ValueError('Source parameter mismatch: '+name)
        old_hashes=io.read(ROOT/cfg[key]/'provenance.json')['code_hashes']
        for name in linear.KERNELS:
            if old_hashes[name]!=io.sha(ROOT/name):
                raise ValueError('Historical physical kernel changed: '+name)
    orbitcfg=io.read(ROOT/cfg['orbit_parent']/'config.json')
    if orbitcfg['degree']!=14 or orbitcfg['A']!=cfg['A'] or orbitcfg['phase_Psi']!=0:
        raise ValueError('The source must be the saved p14 A=.005 Psi=0 orbit')
    return selected,manifests


def build_problem(cfg):
    """Load all complete linear bases; never call a spectral routine."""
    old=load_action(cfg); model6,algebra=six_action(old)
    linear_folder=ROOT/cfg['linear_parent']; orbit_folder=ROOT/cfg['orbit_parent']
    data6=load_npz(linear_folder/'modes_6_p14.npz'); data7=load_npz(linear_folder/'modes_7_p14.npz')
    six=SavedActionSystem(data6,old,6); seven=SavedActionSystem(data7,old,7)
    source=load_npz(orbit_folder/'linear_reference.npz'); U7=load_npz(orbit_folder/'HB13.npz')['U']
    Omega7=io.read(orbit_folder/'HB13.json')['Omega']
    if abs(Omega7-.435115458303789)>2e-12:
        raise ValueError('Saved orbit frequency mismatch')
    linear.require(linear.relative(source['K'],seven.K),cfg['matrix_relative_gate'],'source orbit K basis')
    linear.require(linear.relative(source['M0'],seven.M0),cfg['matrix_relative_gate'],'source orbit M basis')
    transfer=common_transfer(seven,six)
    Phi6=data6['vectors'].copy(); Phi7=source['Phi']; omega6=data6['omega']
    overlaps=Phi6.T@six.M0@transfer@Phi7
    signs=np.ones(Phi6.shape[1]); signs[:6]=np.sign(np.diag(overlaps[:6,:6]))
    Phi6*=signs
    Ucommon=U7@Phi7.T@transfer.T@six.M0@Phi6
    errors=[]
    xi,_=np.polynomial.legendre.leggauss(41); points=(xi+1)/2
    physical7=seven.profiles(U7@Phi7.T,points)[:,:,:,[0,1,5]]
    physical6=six.profiles(Ucommon@Phi6.T,points)[:,:,:,[0,1,5]]
    errors.append(linear.relative(physical6,physical7))
    linear.require(max(errors),cfg['matrix_relative_gate'],'physical seed transfer')
    return dict(six=six,seven=seven,old_action=old,model6=model6,algebra=algebra,
        Phi6=Phi6,Phi7=Phi7,omega6=omega6,omega7=source['omega'],U7=U7,Omega7=Omega7,
        Ucommon=Ucommon,transfer=transfer,signs=signs,transfer_profile_relative=max(errors))


def compare_saved_orbits(problem,U6,Omega6,lowU,lowOmega,cfg):
    d=problem['six']; d7=problem['seven']; Phi=problem['Phi6']; Phi7=problem['Phi7']
    U7=problem['U7']; common=problem['Ucommon']; O7=problem['Omega7']
    phases=np.arange(128)*2*np.pi/128+cfg['independent_phase_offset']
    q6,v6,a6=hb.hb_phases(U6,Omega6,Phi,phases)
    q7,v7,a7=hb.hb_phases(U7,O7,Phi7,phases)
    xi,wg=np.polynomial.legendre.leggauss(41); points=(xi+1)/2
    f6=d.profiles(q6,points)[:,:,:,[0,1,5]]
    f7=d7.profiles(q7,points)[:,:,:,[0,1,5]]
    squares=[]; references=[]
    for i,(L,c) in enumerate(zip(d.lengths,d.coefficients)):
        weights=(wg*L/2)[:,None]*[c.m,c.m,c.jp]
        squares.append(np.sum(weights[None,:,:]*(f6[i]-f7[i])**2,axis=(0,1)))
        references.append(np.sum(weights[None,:,:]*f7[i]**2,axis=(0,1)))
    squares=np.sum(squares,axis=0); references=np.sum(references,axis=0)
    rows=[]
    for model,U,O,ownomega,disc,P in ((7,common,O7,problem['omega7'],d7,Phi7),(6,U6,Omega6,problem['omega6'],d,Phi)):
        original=U7 if model==7 else U
        amp,_,peak=hb.physical_amplitude_gradient(disc,P,original)
        phase=float(np.angle(np.exp(1j*(np.angle(complex(U[3,3]))-3*np.angle(complex(U[1,1]))))))
        rows.append(dict(model=model,Omega=O,period=2*np.pi/O,
            nonlinear_shift=O-ownomega[1],actual_A=amp,peak_arm=peak['arm'],peak_s=peak['s'],
            Q21_common=U[1,1],Q43_common=U[3,3],ratio_common=U[3,3]/U[1,1],phase_Psi_common=phase,
            constant_common_mass_norm=np.linalg.norm(U[0]),
            harmonic2_common_mass_norm=np.linalg.norm(U[2]),harmonic3_common_mass_norm=np.linalg.norm(U[3])))
    Ulo=np.pad(lowU,((0,len(U6)-len(lowU)),(0,0))); factors=np.r_[1.,np.full(len(U6)-1,.5)]
    convergence=dict(frequency_absolute_difference=abs(Omega6-lowOmega),
        common_mass_waveform_relative=float(np.sqrt(np.sum(factors[:,None]*(U6-Ulo)**2)/np.sum(factors[:,None]*U6**2))),
        Q21_difference=float(U6[1,1]-lowU[1,1]),Q43_difference=float(U6[3,3]-lowU[3,3]),
        ratio_difference=float(U6[3,3]/U6[1,1]-lowU[3,3]/lowU[1,1]),
        same_phase_branch=bool(U6[1,1]>0 and U6[3,3]>0 and lowU[1,1]>0 and lowU[3,3]>0))
    data=dict(rows=rows,frequency_difference=Omega6-O7,frequency_relative_difference=(Omega6-O7)/O7,
        common_mass_phase_RMS_relative=float(np.sqrt(sum(squares)/sum(references))),
        component_phase_RMS_relative=dict(zip(('u','w','theta'),np.sqrt(squares/references).tolist())),
        common_modal_basis=cfg['common_modal_basis'],harmonic_convergence=convergence)
    arrays=dict(phases=phases,material_points_fraction=points,common_profiles_6=f6,common_profiles_7=f7,
        U6_common=U6,U7_common=common,omega_reference=problem['omega6'],Phi_common=Phi,
        harmonic_common_mass_norms=np.array([np.linalg.norm(common,axis=1),np.linalg.norm(U6,axis=1)]),
        harmonic_reference_energy_norms=np.array([np.linalg.norm(common*problem['omega6'],axis=1),np.linalg.norm(U6*problem['omega6'],axis=1)]))
    return data,arrays


def reconstruct_orbit(bundle,times,harmonics=13):
    """Saved six-field q(t), qdot(t), qddot(t); no assembly or solve."""
    z=load_npz(Path(bundle)/f'HB{harmonics}.npz')
    return hb.hb_phases(z['U'],float(z['Omega']),z['Phi'],np.atleast_1d(times)*float(z['Omega']))


def load_operator(bundle):
    """Full spatial residual/tangent interface for a separately authorized next stage."""
    bundle=Path(bundle); io.cache(bundle)
    cfg=io.read(bundle/'config.json')
    return SavedActionSystem(load_npz(bundle/'linear_6_p14.npz'),load_action(cfg),6)


def render_saved(bundle):
    """Two plots and a concise local report from saved arrays, with no physics calls."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    bundle=Path(bundle); s=io.cache(bundle)
    if not (bundle/'comparison.npz').exists():
        return s
    z=load_npz(bundle/'comparison.npz'); x=z['material_points_fraction']; phases=z['phases']
    # Physical fields at a common phase, shown by their midpoint waveforms over a full period.
    fig,axes=plt.subplots(3,2,figsize=(11,8),sharex=True)
    for arm in range(2):
        for f,label in enumerate(('u','w','theta')):
            for model,style in ((7,'-'),(6,'--')):
                values=z[f'common_profiles_{model}'][arm]
                mid=np.array([np.interp(.5,x,row[:,f]) for row in values])
                axes[f,arm].plot(phases/(2*np.pi),mid,style,label=f'{model} fields')
            axes[f,arm].set_ylabel(label+' (model units)'); axes[f,arm].grid(alpha=.25)
        axes[0,arm].set_title(f'Arm {arm+1}, s/L=0.5')
        axes[2,arm].set_xlabel('common phase tau / (2 pi)')
    axes[0,0].legend(); fig.suptitle('Physical planar orbit; each model at its own period, A=0.005')
    fig.tight_layout(); fig.savefig(bundle/'physical_orbits.png',dpi=180); fig.savefig(bundle/'physical_orbits.pdf'); plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for index,(model,style) in enumerate(((7,'o-'),(6,'s--'))):
        U=z[f'U{model}_common']; norms=np.linalg.norm(U,axis=1)
        axes[0].semilogy(np.arange(len(U)),np.maximum(norms,1e-20),style,label=f'{model} fields')
        axes[1].plot([1,3],[U[1,1],U[3,3]],style,label=f'{model} fields')
    axes[0].set_ylabel('common mass norm of cosine coefficient'); axes[0].set_xlabel('harmonic n (0 = mean)')
    axes[1].set_xticks([1,3],['PL2 at n=1','PL4 at n=3']); axes[1].set_ylabel('signed common modal coefficient')
    for ax in axes:ax.grid(alpha=.25); ax.legend()
    fig.suptitle('Composition in the same signed six-field linear PL basis')
    fig.tight_layout(); fig.savefig(bundle/'harmonic_composition.png',dpi=180); fig.savefig(bundle/'harmonic_composition.pdf'); plt.close(fig)
    lines=['# Six-field nonlinear operator and planar orbit','',
        f"Operator: {s['operator_status']}; orbit: {s['orbit_status']}; stability: NOT_EVALUATED.",'',
        '|model|Omega|period|nonlinear shift|Q21 common|Q43 common|ratio|',
        '|---:|---:|---:|---:|---:|---:|---:|']
    for r in s['comparison']['rows']:
        lines.append('|'+ '|'.join(str(r[k]) for k in ('model','Omega','period','nonlinear_shift','Q21_common','Q43_common','ratio_common'))+'|')
    (bundle/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    io.finalize(bundle,io.read(bundle/'provenance.json'),s)
    return s


def run():
    cfg=io.read(CONFIG)
    if (cfg['model']!=MODEL or cfg['degree']!=14 or cfg['A']!=.005 or cfg['HB_harmonics']!=[11,13]
        or cfg['max_Newton_steps']>5 or cfg['numerical_budget_seconds']>1800 or cfg['linear_eigensolves'] or cfg['ODE_Hill_Floquet_FEM']):
        raise ValueError('Only the explicitly bounded p14 A=.005 HB11/HB13 task is supported')
    selected,manifests=verified_sources(cfg)
    inputs=dict(schema=VERSION,model=MODEL,config=cfg,source_artifacts=selected,source_manifests=manifests,
                code_hashes={f:io.sha(ROOT/f) for f in CODE})
    fingerprint=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()[:16]
    bundle=OUTPUT/fingerprint
    if (bundle/'manifest.json').exists():
        return bundle,io.cache(bundle),True
    if bundle.exists():
        raise ValueError('Preserve unfinished unmanifested attempt')
    # A code change does not itself authorize a second orbit attempt.
    previous=list(OUTPUT.glob('*/provenance.json')) if OUTPUT.exists() else []
    if previous:
        raise ValueError('An earlier attempt exists; inspect its numerical state before an authorized retry')
    bundle.mkdir(parents=True)
    provenance=dict(inputs,HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        status_short=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        initial_state=io.read(ROOT/'results/_smoke/nlsp_six_vs_seven_planar_orbit/initial_state.json'),
        interpreter=sys.executable,started_UTC=datetime.now(timezone.utc).isoformat(),
        conventions='saved NLSP clamp-to-joint axes, q6=(u,w,v,Phi,psi,theta), a=(Phi,-psi,theta)',
        independent_c='absent; embedded zero only inside seven-field local evaluator')
    io.write(bundle/'config.json',cfg); io.write(bundle/'provenance.json',provenance)
    for name in CODE:
        dest=bundle/'execution_code'/name; dest.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(ROOT/name,dest)
    s=dict(overall='PREPARED',completed=False,operator_status='NOT_CHECKED',orbit_status='NOT_RUN',
        spatial_stability='NOT_EVALUATED',calls=dict(linear_eigensolves=0,HB_solves=0,HB13_retries=0,ODE=0,Hill=0,Floquet=0,FEM=0))
    io.finalize(bundle,provenance,s)
    start=time.perf_counter()
    try:
        problem=build_problem(cfg); d=problem['six']; Phi=problem['Phi6']; omega=problem['omega6']
        operator=operator_checks(d,problem['model6'],problem['old_action'],cfg)
        operator.update(problem['algebra']); operator['seed_transfer_relative']=problem['transfer_profile_relative']
        io.write(bundle/'operator_checks.json',operator)
        io.write(bundle/'six_field_action.json',dict(model=MODEL,fields=list(rod.FIELD_ORDER[:6]),
            coefficient_contract='C=EA/(1-nu^2); K=C*(1-nu^2) is explicit in serialized polynomials',
            T4=problem['model6'].T4.serialize(),V4=problem['model6'].V4.serialize(),
            residuals=[r.serialize() for r in problem['model6'].residual_a]))
        shutil.copyfile(ROOT/cfg['linear_parent']/'modes_6_p14.npz',bundle/'linear_6_p14.npz')
        np.savez_compressed(bundle/'recovery.npz',W=d.W,raw_maps=d.raw_maps,local_evaluator_maps=np.array(d.maps),
            plane=d.plane,out=d.out,state_plane=d.state_plane,state_out=d.state_out,fields=np.array(rod.FIELD_ORDER[:6]),
            labels=d.saved['labels'],Phi=Phi,omega=omega,transfer_7_to_6=problem['transfer'],
            signs=problem['signs'],lengths=d.lengths,bases=d.bases,p=14,nq=31)
        s['operator_status']='PASS'; s['operator_seconds']=time.perf_counter()-start
        io.finalize(bundle,provenance,s)
        orbit_start=time.perf_counter(); deadline=orbit_start+cfg['numerical_budget_seconds']
        solver=SixPlanarHB(d,Phi,omega,cfg,deadline)
        seed=problem['Ucommon'][:12].copy(); Oseed=problem['Omega7']
        initial,arrays=solver.validate(seed,Oseed,cfg['A'])
        initial['projected_residual_unscaled_max']=float(np.max(abs(solver.equations(seed,Oseed,cfg['A']))))
        io.write(bundle/'transferred_seed.json',initial); np.savez_compressed(bundle/'transferred_seed.npz',**arrays)
        U,O=seed,Oseed; rows=[]; solutions=[]
        for H in cfg['HB_harmonics']:
            t0=time.perf_counter(); s['calls']['HB_solves']+=1; io.finalize(bundle,provenance,s)
            U,O,newton=solver.solve(U,O,cfg['A'],H)
            row,arrays=solver.validate(U,O,cfg['A']); row['Newton']=newton
            row['acceptance']=acceptance(row,cfg); row['seconds']=time.perf_counter()-t0
            rows.append(row); solutions.append((U.copy(),O)); s['HB']=rows
            np.savez_compressed(bundle/f'HB{H}.npz',**arrays); io.write(bundle/f'HB{H}.json',row)
            io.finalize(bundle,provenance,s)
        comparison,arrays=compare_saved_orbits(problem,U,O,*solutions[0],cfg)
        # Seven-field residual/energy are reused with their original verification grid.
        oldrow=io.read(ROOT/cfg['orbit_parent']/'HB13.json')
        comparison['rows'][0].update(full_residual_relative=oldrow['full_residual_relative'],
            energy_relative_variation=oldrow['energy_relative_variation'],verification_phases=oldrow['independent_phase_count'])
        comparison['rows'][1].update(full_residual_relative=rows[-1]['full_residual_relative'],
            energy_relative_variation=rows[-1]['energy_relative_variation'],verification_phases=128)
        s['comparison']=comparison; np.savez_compressed(bundle/'comparison.npz',**arrays)
        linear.write_csv(bundle/'orbit_comparison.csv',comparison['rows']); io.write(bundle/'comparison.json',comparison)
        conv=comparison['harmonic_convergence']
        temporal=conv['same_phase_branch'] and conv['frequency_absolute_difference']<=cfg['frequency_control_gate']
        s['temporal_convergence']=temporal
        passed=all(r['acceptance']['passed'] for r in rows) and temporal
        s['orbit_status']='PASS' if passed else 'PARTIAL: see independent acceptance gates'
        s['overall']='SIX_FIELD_OPERATOR_AND_PLANAR_ORBIT_COMPLETE' if passed else 'SIX_FIELD_OPERATOR_READY_ORBIT_PARTIAL'
        s['completed']=True; s['orbit_seconds']=time.perf_counter()-orbit_start
        s['operator_checks']=operator; s['phase_evaluations']=solver.calls
        io.finalize(bundle,provenance,s); render_saved(bundle)
    except (ArithmeticError,ValueError,TimeoutError,np.linalg.LinAlgError) as exc:
        s.update(overall='BOUNDED_ATTEMPT_STOPPED',completed=True,failure=repr(exc),elapsed_seconds=time.perf_counter()-start)
        io.finalize(bundle,provenance,s)
        raise
    return bundle,s,False


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--compute',action='store_true')
    modes.add_argument('--report-only',type=Path); modes.add_argument('--plot-only',type=Path)
    args=parser.parse_args(argv)
    if args.compute:
        bundle,s,cached=run()
    else:
        bundle=(args.report_only or args.plot_only).resolve(); s=render_saved(bundle); cached=True
    print(json.dumps(io.plain(dict(bundle=str(bundle),overall=s['overall'],operator=s['operator_status'],
        orbit=s['orbit_status'],spatial_stability=s['spatial_stability'],reused=cached,
        new_HB_solves=0 if cached else s['calls']['HB_solves'],comparison=s.get('comparison'))),indent=2))
