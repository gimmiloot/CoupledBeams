"""Bounded horizontal-support verification, with four sequential hard gates.

Script Proliferation Control: a new physical admissible space and zero-inclusive
contract cannot be run by the historical positive-only map CLI. This is the
single reusable verification entry point; it creates no parameter maps.
"""
from pathlib import Path
import argparse
from dataclasses import asdict, replace
import hashlib
import json
import math
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import numpy as np
import scipy
from fractions import Fraction
from scipy.optimize import brentq
from scripts.lib import isotropic_two_springs_horizontal_spring as h
from scripts.analysis.verify_isotropic_two_springs_rl_timoshenko import sha,write,csvwrite

CONFIG = ROOT/'data/input/isotropic_two_springs_horizontal_spring.json'
OUTPUT = ROOT/'results/isotropic_two_springs_rl_timoshenko/horizontal_spring_verification'


def identity(cfg):
    files = [CONFIG,Path(__file__),Path(h.__file__),Path(h.base.__file__),
             Path(h.base.axial.__file__),Path(h.base.tim.__file__),
             Path(h.base.geometry.__file__),Path(h.base.old_joint.__file__),
             ROOT/'scripts/lib/reddy_symmetric_coupled_beams_ritz.py',
             ROOT/'scripts/lib/reddy_inplane_geometry.py',
             ROOT/'scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py',
             ROOT/'scripts/analysis/verify_isotropic_two_springs_rl_timoshenko.py',
             ROOT/cfg['pole_source']]
    record = {'config':cfg,'files':{p.relative_to(ROOT).as_posix():sha(p) for p in files},
              'versions':{'python':platform.python_version(),'numpy':np.__version__,
                          'scipy':scipy.__version__}}
    key = hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()[:16]
    return key,record


def mathematical_gate():
    # Direct Taylor/variational proof is in the canonical note. Exact-rational
    # matrix checks supplement it; neither SymPy nor Lean is required/claimed.
    c,s=Fraction(3,5),Fraction(4,5)
    matrix=np.outer([c,s],[c,s])
    force=np.array([[-c,-s],[-s,c]],dtype=object)@np.array([c,s])
    assert matrix[0,0]*matrix[1,1]-matrix[0,1]*matrix[1,0]==0
    assert matrix.trace()==1 and force.tolist()==[-1,0]
    # Static first-order solution terminates at cubic order.
    model=h.Model()
    assert np.max(abs(np.linalg.matrix_power(h.state_matrix(model,0),4)))==0
    return {'status':'PASS','quadratic_spring_energy_twice':'k*X**2',
            'local_spring_hessian':'k*[[c*c,c*s],[c*s,s*s]]',
            'boundary_variation':['-N+k*d*c','-Q+k*d*s','-M'],
            'physical_internal_force':['-k*d',0],'spring_rank_positive_k':1,
            'zero_kernel_proof':'all rod strain energies zero plus first clamp; only second-arm rigid rotation remains',
            'mathematical_mismatch':None,
            'formal_system':'direct variational/Taylor derivation and exact rational checks; SymPy unavailable in the project interpreter; Lean not used'}


class Run:
    def __init__(self,cfg,folder):
        self.cfg,self.folder,self.p = cfg,folder,cfg['policy']
        self.m = h.Model(**cfg['material'])
        if cfg['low_frequency_basis_Lambda'] != h.LOW_FREQUENCY_BASIS_LAMBDA:
            raise ValueError('Low-frequency basis switch differs from contract')
        self.cases = [h.Case(b,cfg['mu'],j,k) for b,j,k in cfg['controls']]
        source=ROOT/cfg['pole_source']
        old=json.loads((source.parent/'contract.json').read_text(encoding='utf8'))
        oldcfg=old.get('config',old)
        if oldcfg['material'] != cfg['material'] or oldcfg['mu'] != cfg['mu']:
            raise ValueError('Arm-pole material/geometry mismatch')
        self.catalogs = {float(k):v for k,v in json.loads(source.read_text()).items()}
        if set(self.catalogs) != set(self.cases[0].lengths(self.m)):
            raise ValueError('Arm-pole lengths mismatch')
        self.solved={}
        self.ritz_solves=0

    def name(self,case):
        return hashlib.sha256(json.dumps(asdict(case),sort_keys=True).encode()).hexdigest()[:16]

    def kernel_gate(self):
        rows=[]
        for case in self.cases:
            static=h.static_kernel(self.m,case,self.p['kernel_rank_relative_tol'])
            if static['status']!='PASS':
                raise ArithmeticError(f'Static kernel: {case} {static}')
            ritz={str(order):h.ritz(self.m,case,order,rank_tol=self.p['kernel_rank_relative_tol'],kernel_only=True)
                  for order in self.p['ritz_orders']}
            for rd in ritz.values():
                if rd['kernel_residual']>self.p['kernel_residual_tol']:
                    raise ArithmeticError('Ritz zero residual failed')
            record={'case':asdict(case),'static':static,'ritz':ritz}
            if static['expected']:
                mode=h.mode(self.m,case,0)
                L2=case.lengths(self.m)[1]
                J2=self.m.section.rhoA*L2**3/3+self.m.section.rhoI*L2
                err=abs(mode['diagnostics']['raw_mass']/J2-1)
                if err>self.p['zero_mass_tol']:
                    raise ArithmeticError('Zero-mode inertia mismatch')
                record.update(J2=J2,mass_error=err,zero_diagnostics=mode['diagnostics'])
            rows.append(record)
        result={'status':'PASS','cases':rows,'new_static_checks':len(rows),'new_ritz_kernel_checks':2*len(rows)}
        write(self.folder/'gate2_kernel.json',result)
        return result

    def solve(self,case,label,ritz=False):
        key=self.name(case)
        if key in self.solved:
            return self.solved[key]
        record={'case':asdict(case),'label':label,'status':'UNRESOLVED'}
        write(self.folder/'points'/f'{key}.json',record)
        roots,count=h.prefix(self.m,case,self.catalogs,self.p)
        omega=np.array([r['omega'] for r in roots])
        modes=[h.mode(self.m,case,w,self.p['quadrature_order']) for w in omega]
        diags=[a['diagnostics'] for a in modes]
        gram=np.array([[h.inner(self.m,a,b) for a in modes] for b in modes])
        error=float(np.max(abs(gram-np.eye(3))))
        record.update(roots=roots,count=count,omega=omega.tolist(),Lambda=self.m.Lambda(omega).tolist(),
                      diagnostics=diags,mass_gram=gram.tolist(),mass_gram_error=error)
        write(self.folder/'points'/f'{key}.json',record)
        for d in diags:
            metrics={'boundary_tol':max(np.max(d['joint_scaled']),np.max(d['external_scaled'])),
                     'ode_tol':d['ode_error'],'energy_tol':d['energy_error'],
                     'singular_tol':d['singular_ratio'],'nonzero_condition_max':d['nonzero_condition']}
            for name,value in metrics.items():
                if value>self.p[name]:
                    raise ArithmeticError(f'{case}: {name}={value} > {self.p[name]}')
        if error>self.p['mass_gram_tol']:
            raise ArithmeticError(f'Mass orthogonality failed: {error}')
        archive={}
        for pos,mo in enumerate(modes,1):
            for arm,a in enumerate(mo['arms'],1):
                for name,val in a.items():
                    archive[f'mode{pos}_arm{arm}_{name}']=val
        if ritz:
            checks=[]
            for order in self.p['ritz_orders']:
                rw,co,rd=h.ritz(self.m,case,order,rank_tol=self.p['kernel_rank_relative_tol'])
                self.ritz_solves+=1
                positive=omega>0
                err=float(np.max(abs(rw[positive]/omega[positive]-1)))
                overlaps=[abs(h.inner(self.m,mo,h.ritz_fields(self.m,case,order,co[:,i]))) for i,mo in enumerate(modes)]
                row={'order':order,'omega':rw.tolist(),'frequency_error':err,'mass_overlaps':overlaps,**rd}
                checks.append(row)
                record['ritz']=checks
                write(self.folder/'points'/f'{key}.json',record)
                if err>self.p['ritz_frequency_tol'] or max(abs(np.array(overlaps)-1))>self.p['shape_overlap_tol']:
                    raise ArithmeticError(f'Independent Ritz comparison failed: {row}')
                if rd['mass_gram_error']>self.p['mass_gram_tol'] or rd['energy_error']>self.p['energy_tol']:
                    raise ArithmeticError('Ritz mass/energy check failed')
                if rd['constraint_residual']>self.p['boundary_tol'] or rd['matrix_residual']>self.p['singular_tol']:
                    raise ArithmeticError('Ritz constraint/matrix residual failed')
                archive[f'ritz{order}_coefficients']=co
            record['ritz_order_frequency_error']=float(max(abs(np.array(checks[0]['omega'])[positive]/np.array(checks[1]['omega'])[positive]-1)))
        record['status']='PASS'
        write(self.folder/'points'/f'{key}.json',record)
        np.savez_compressed(self.folder/'points'/f'{key}.npz',**archive)
        self.solved[key]=(record,modes)
        print(label,asdict(case),'PASS',flush=True)
        return record,modes

    def control_gate(self):
        records=[self.solve(case,'control',ritz=True)[0] for case in self.cases]
        # No modal tracking: compare sorted positions at unchanged geometry.
        monotonic=[]
        for parameter,other in [('kappa_horizontal','k_joint'),('k_joint','kappa_horizontal')]:
            for beta in (0,5,90):
                groups={}
                for r in records:
                    if r['case']['beta']==beta:
                        groups.setdefault(str(r['case'][other]),[]).append(r)
                for key,rows in groups.items():
                    rows.sort(key=lambda r:float('inf') if r['case'][parameter]==h.RIGID else r['case'][parameter])
                    for left,right in zip(rows[:-1],rows[1:]):
                        delta=np.array(right['omega'])-np.array(left['omega'])
                        lower=-self.p['ritz_frequency_tol']*np.maximum(right['omega'],1e-10)
                        if np.any(delta<lower):
                            raise ArithmeticError('Sorted monotonicity failed')
                        monotonic.append({'parameter':parameter,'beta':beta,'other':key,
                                          'from':left['case'][parameter],'to':right['case'][parameter],
                                          'omega_difference':delta.tolist()})
        straight=[r for r in records if r['case']['beta']==0 and r['case']['k_joint']==.1]
        invariance=float(np.max(abs(np.array([r['omega'] for r in straight])/straight[0]['omega']-1)))
        if invariance>self.p['ritz_frequency_tol']:
            raise ArithmeticError('Straight bending independence failed')
        # Reflection tests a genuinely different angle, no sign-copying of spectra.
        reflection=self.solve(h.Case(-5,self.cfg['mu'],.1,1),'reflection',ritz=True)[0]
        plus=self.solve(h.Case(5,self.cfg['mu'],.1,1),'control')[0]
        even=float(max(abs(np.array(reflection['omega'])/plus['omega']-1)))
        if even>self.p['ritz_frequency_tol']:
            raise ArithmeticError('Reflection symmetry failed')
        # Independent SINGLE-arm cantilever characteristic, not the node matrix.
        cant=[r for r in records if r['case']['beta']==0 and r['case']['k_joint']==h.RIGID][0]
        def characteristic(w):
            b=h.base.tim.finite_state_basis(self.m.carrier,2*self.m.l,w,[0.,2*self.m.l],'timoshenko')
            a=np.vstack((b[0,:2],b[1,2:]))
            return np.linalg.det(h.base.balance(a)[0])
        single=[brentq(characteristic,.98*w,1.02*w,xtol=self.p['root_xtol'],rtol=self.p['root_rtol']) for w in cant['omega']]
        ce=float(max(abs(np.array(single)/cant['omega']-1)))
        if ce>self.p['ritz_frequency_tol']:
            raise ArithmeticError('Single cantilever limit failed')
        # Inertia convention immediately on either side of EACH arm pole.
        pole_audits=[]
        test=h.Case(5,self.cfg['mu'],.1,1)
        upper=float(self.m.omega(self.p['Lambda_upper']))
        for pole in h.poles(self.m,test,self.catalogs):
            if pole>=upper:
                continue
            eps=4*self.p['pole_exclusion']*max(1,pole)
            left=h.count_at(self.m,test,pole-eps,self.catalogs,self.p)[2]
            right=h.count_at(self.m,test,pole+eps,self.catalogs,self.p)[2]
            if left['count']!=right['count'] or right['J0']-left['J0']!=1:
                raise ArithmeticError('Arm-pole inertia cancellation unresolved')
            pole_audits.append({'pole':pole,'left':left,'right':right})
        result={'status':'PASS','controls':len(records),'reflection_frequency_error':even,
                'straight_bending_independence_error':invariance,'cantilever_omega':single,
                'cantilever_frequency_error':ce,'monotonicity':monotonic,'pole_audits':pole_audits}
        write(self.folder/'gate3_controls.json',result)
        return result

    def asymptotic_gate(self):
        soft=[]
        L1,L2=self.cases[0].lengths(self.m)
        J2=self.m.section.rhoA*L2**3/3+self.m.section.rhoI*L2
        for item in self.cfg['soft_local_checks']:
            kh=item['kappa_horizontal']
            k=None if kh==h.RIGID else kh*self.m.section.EI/self.m.l**3
            keff=self.m.section.EA/(L1+L2) if k is None else 1/(1/k+(L1+L2)/self.m.section.EA)
            predicted=keff*L2**2/J2
            errors=[]
            for angle in item['angles']:
                r,_=self.solve(h.Case(angle,self.cfg['mu'],0,kh),'soft_local',ritz=(angle==item['angles'][-1]))
                measured=r['omega'][0]**2/math.radians(angle)**2
                error=measured/predicted-1
                errors.append(error)
                soft.append({'kappa_horizontal':kh,'beta':angle,'coefficient_predicted':predicted,
                             'coefficient_observed':measured,'signed_relative_error':error})
            if abs(errors[-1])>self.p['asymptotic_relative_tol']:
                raise ArithmeticError('Soft-hinge coefficient unresolved')
        finite=[]
        cfg=self.cfg['finite_joint_local_checks']
        free_case=h.Case(0,self.cfg['mu'],cfg['kappa_joint'],0)
        r0,mo0=self.solve(free_case,'control')
        w0=r0['omega'][0]
        pred0=h.curvature_prediction(self.m,free_case,w0,mo0[0])
        for kh in cfg['kappa_horizontal']:
            case=replace(free_case,kappa_horizontal=kh)
            r,mos=self.solve(case,'control')
            direct=h.curvature_prediction(self.m,case,w0,mos[0])
            ca=pred0['axial_dynamic_compliance']
            k=None if kh==h.RIGID else kh*self.m.section.EI/self.m.l**3
            factor=1/ca if k is None else k/(1+k*ca)
            predicted=pred0['C2_direct']+pred0['gamma0']**2*factor
            cross=abs(direct['C2_direct']-predicted)/max(abs(predicted),1e-30)
            if cross>self.p['ritz_frequency_tol']:
                raise ArithmeticError(f'Perturbation formula mismatch: {kh} {cross}')
            errors=[]
            for angle in cfg['angles']:
                rr,_=self.solve(replace(case,beta=angle),'finite_local')
                measured=(rr['omega'][0]**2-r['omega'][0]**2)/math.radians(angle)**2
                error=measured/predicted-1
                errors.append(error)
                finite.append({'kappa_horizontal':kh,'beta':angle,'C2_predicted':predicted,
                               'C2_direct':direct['C2_direct'],'C2_observed':measured,
                               'signed_relative_error':error,'formula_crosscheck':cross})
            if abs(errors[-1])>self.p['asymptotic_relative_tol']:
                raise ArithmeticError('Finite-joint curvature unresolved')
        result={'status':'PASS','soft_hinge':soft,'finite_joint':finite,'unloaded_prediction':pred0}
        write(self.folder/'gate4_asymptotics.json',result)
        csvwrite(self.folder/'soft_hinge.csv',soft)
        csvwrite(self.folder/'finite_joint_curvature.csv',finite)
        return result

    def numerical_recovery_record(self):
        """Reconstruct the two diagnosed numerical failures, not new physics."""
        case=h.Case(.0625,self.cfg['mu'],0,1)
        record,modes=self.solved[self.name(case)]
        a,_=h.assembled(self.m,case,0)
        oldscale=h.congruence(self.m,case)
        oldvals=np.linalg.eigvalsh((a+a.T)/2*oldscale[:,None]*oldscale)
        static=h.static_kernel(self.m,case,self.p['kernel_rank_relative_tol'])
        # Saved failed root, never substituted for an accepted result.
        failed_omega=9.399404750661146e-06
        original=h.basis
        try:
            h.basis=h.base.basis
            oldmode=h.mode(self.m,case,failed_omega)
        finally:
            h.basis=original
        newmode=modes[0]
        diag={'case':asdict(case),'original_archives':['fe9526b6428e4031','5c2f6d2583135787'],
              'old_static_eigenvalues':oldvals.tolist(),
              'old_static_relative_minimum':float(min(abs(oldvals))/max(abs(oldvals))),
              'new_static':static,'old_frequency':failed_omega,'new_frequency':record['omega'][0],
              'old_mode_diagnostics':oldmode['diagnostics'],'new_mode_diagnostics':newmode['diagnostics'],
              'absolute_mass_overlap':abs(h.inner(self.m,oldmode,newmode)),
              'unchanged_tolerances':True,'mechanical_equations_changed':False}
        write(self.folder/'numerical_recovery.json',diag)
        return diag


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['kernel','controls','all'],default='all')
    parser.add_argument('--verify-cache',type=Path)
    args=parser.parse_args()
    cfg=json.loads(CONFIG.read_text(encoding='utf8'))
    key,ident=identity(cfg)
    if args.verify_cache:
        manifest=json.loads((args.verify_cache/'manifest.json').read_text(encoding='utf8'))
        if manifest['identity']!=ident or any(sha(args.verify_cache/f)!=v for f,v in manifest['artifacts'].items()):
            raise ValueError('Archive identity/hash mismatch')
        print('Cache verified; no eigenvalue calculations');return
    folder=OUTPUT/key
    folder.mkdir(parents=True,exist_ok=True)
    if (folder/'manifest.json').exists():
        raise FileExistsError('Completed immutable archive exists; use --verify-cache')
    write(folder/'contract.json',ident)
    gates={}
    run=None
    try:
        gates['gate1']=mathematical_gate()
        write(folder/'gate1_mathematics.json',gates['gate1'])
        run=Run(cfg,folder)
        gates['gate2']=run.kernel_gate()
        print('Gate 1/2 PASS',flush=True)
        if args.stage!='kernel':
            gates['gate3']=run.control_gate()
            print('Gate 3 PASS',flush=True)
        if args.stage=='all':
            gates['gate4']=run.asymptotic_gate()
            run.numerical_recovery_record()
            print('Gate 4 PASS',flush=True)
    except Exception as exc:
        write(folder/'failure.json',{'status':'PARTIAL_PASS' if 'gate2' in gates else 'FAIL',
              'error':repr(exc),'completed_gates':list(gates),
              'accepted_spectral_points':len(run.solved) if run else 0})
        print(folder,flush=True)
        raise
    rows=[]
    for r,_ in run.solved.values():
        rows.append({**r['case'],'label':r['label'],'zero_count':r['count']['zero_gate']['expected'],
                     **{f'omega{i+1}':w for i,w in enumerate(r['omega'])},
                     **{f'Lambda{i+1}':w for i,w in enumerate(r['Lambda'])},
                     'mass_gram_error':r['mass_gram_error']})
    csvwrite(folder/'spectra.csv',rows)
    status='PASS' if args.stage=='all' else 'PARTIAL_PASS'
    manifest={'status':status,'identity':ident,'git_head':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
              'initial_git_status':cfg['initial_git_status'],'new_spectral_points':len(rows),
              'reused_spectral_points':0,'reused_arm_pole_catalogs':2,'ritz_solves':run.ritz_solves,
              'stage':args.stage,'artifacts':{p.relative_to(folder).as_posix():sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}}
    write(folder/'manifest.json',manifest)
    print(folder,flush=True)


if __name__=='__main__':
    main()
