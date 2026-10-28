"""HS-P1: fixed zero-inclusive pilot, reusing the verified horizontal solver.

Script Proliferation Control: this module is the stable pilot entry point
(`python -m ...`), not a new spectral solver. Its map/reuse/checkpoint/plot-only
contract differs from the frozen four-gate verification CLI. That CLI and its
source-hashed archive remain byte-identical; no second analysis script is added.
"""
from pathlib import Path
from dataclasses import asdict
import argparse
import hashlib
import json
import math
import platform
import time

import numpy as np
import scipy
from scripts.analysis import verify_two_springs_horizontal_support as verified
from scripts.lib import isotropic_two_springs_horizontal_spring as physics

ROOT=Path(__file__).resolve().parents[2]
CONFIG=ROOT/'data/input/isotropic_two_springs_horizontal_pilot.json'
OUTPUT=ROOT/'results/isotropic_two_springs_rl_timoshenko/horizontal_spring_pilot'
sha,write,csvwrite=verified.sha,verified.write,verified.csvwrite


def require(ok,message):
    if not ok:
        raise ArithmeticError(message)


def key(case):
    return (float(case.beta),case.k_joint if case.k_joint==physics.RIGID else float(case.k_joint),case.kappa_horizontal if
            case.kappa_horizontal==physics.RIGID else float(case.kappa_horizontal))


def identity(cfg):
    record={'config':cfg,'workflow_sha256':sha(__file__),'config_sha256':sha(CONFIG),
            'initial_state_sha256':sha(ROOT/cfg['initial_workspace']),
            'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__}}
    return hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()[:16],record


def validate_record(record,case,policy):
    require(record['status']=='PASS','Unaccepted source point')
    require(record['case']==asdict(case),'Physical case mismatch')
    require(record['count']['status']=='PASS' and record['count']['physical_count']==3,'Incomplete physical prefix')
    require(record['count']['zero_gate']['status']=='PASS','Unresolved zero gate')
    omega=np.asarray(record['omega'])
    require(len(omega)==3 and np.all(np.diff(omega)>0),'Unsorted/duplicate physical spectrum')
    require(np.count_nonzero(omega==0)==physics.expected_nullity(case),'Wrong exact-zero multiplicity')
    require([r['physical_position'] for r in record['roots']]==[1,2,3],'Physical numbering mismatch')
    require(record['mass_gram_error']<=policy['mass_gram_tol'],'Mass orthogonality')
    for d in record['diagnostics']:
        require(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))<=policy['boundary_tol'],'Boundary residual')
        for field,tol in [('ode_error','ode_tol'),('energy_error','energy_tol'),
                          ('singular_ratio','singular_tol'),('nonzero_condition','nonzero_condition_max')]:
            require(d[field]<=policy[tol],field)


def load_modes(point_file):
    with np.load(point_file.with_suffix('.npz')) as z:
        modes=[]
        for pos in range(1,4):
            arms=[{name:z[f'mode{pos}_arm{arm}_{name}'].copy() for name in ('x','weights','y','dy','ends')}
                  for arm in (1,2)]
            modes.append({'arms':arms})
    return modes


def energy_parts(model,case,mode,omega):
    """Twice potential energies and full mass, no modal-descendant assignment."""
    s=model.section
    strain=[];masses=[]
    for a in mode['arms']:
        y,dy,wt=a['y'],a['dy'],a['weights']
        strain.append([float(wt@(s.EA*dy[:,0]**2)),float(wt@(s.EI*dy[:,2]**2)),
                       float(wt@(s.KGA*(dy[:,1]-y[:,2])**2))])
        masses.append(float(wt@(s.rhoA*(y[:,0]**2+y[:,1]**2)+s.rhoI*y[:,2]**2+model.J*dy[:,0]**2)))
    jump=mode['arms'][0]['ends'][1,2]-mode['arms'][1]['ends'][1,2]
    c,ss=physics.cs(case)
    d=c*mode['arms'][1]['ends'][0,0]+ss*mode['arms'][1]['ends'][0,1]
    ej=0. if case.k_joint==physics.RIGID else case.k_joint*s.EI/model.l*jump**2
    eh=0. if case.kappa_horizontal==physics.RIGID else case.kappa_horizontal*s.EI/model.l**3*d**2
    array=np.array(strain)
    parts=dict(zip(('axial','bending','shear'),array.sum(axis=0)))
    parts.update(joint=float(ej),horizontal=float(eh))
    total=float(sum(parts.values()));mass=sum(masses)
    a=mode['arms'][1];L=case.lengths(model)[1]
    x,wt,y,dy=a['x'],a['weights'],a['y'],a['dy']
    relu=y[:,0]-a['ends'][1,0];relw=y[:,1]-a['ends'][1,1]
    J2=s.rhoA*L**3/3+s.rhoI*L
    rotation=float(wt@(s.rhoA*(x-L)*relw+s.rhoI*y[:,2])/J2)
    numerator=float(wt@(s.rhoA*(relu**2+(relw-rotation*(x-L))**2)+s.rhoI*(y[:,2]-rotation)**2+model.J*dy[:,0]**2))
    denominator=float(wt@(s.rhoA*(relu**2+relw**2)+s.rhoI*y[:,2]**2+model.J*dy[:,0]**2))
    out={'mass':mass,'energy_twice':total,'energy_error':abs(total/(omega**2*mass)-1) if omega else total,
         'arm1_mass_fraction':masses[0]/mass,'arm2_mass_fraction':masses[1]/mass,
         'arm1_strain_fraction':float(array[0].sum()/total) if total else None,
         'arm2_strain_fraction':float(array[1].sum()/total) if total else None,
         'rigid_rotation_residual':math.sqrt(max(0,numerator/denominator)) if denominator else None,
         'best_relative_rotation':rotation}
    for name,value in parts.items():
        out[name+'_energy_twice']=float(value)
        out[name+'_fraction']=float(value/total) if total else None
    return out


class Pilot:
    def __init__(self,cfg):
        self.cfg=cfg; self.ident_key,self.ident=identity(cfg)
        self.folder=OUTPUT/self.ident_key
        self.source=ROOT/cfg['source_archive']
        self.model=physics.Model(**cfg['material'])
        self.p=cfg['policy']; self.rows=[]; self.records={}
        self.reference={}; self.locators=[]; self.audits=[]; self.failures=[]
        self.new=self.reused=self.resumed=self.ritz_count=0

    def case(self,construction,beta):
        return physics.Case(beta,self.cfg['mu'],construction['k_joint'],construction['kappa_horizontal'])

    def audit_source(self):
        require(sha(self.source/'manifest.json')==self.cfg['source_manifest_sha256'],'Source manifest changed')
        manifest=json.loads((self.source/'manifest.json').read_text())
        for name,digest in manifest['artifacts'].items():
            require(sha(self.source/name)==digest,'Changed source artifact '+name)
        for name,digest in self.cfg['frozen_files'].items():
            require(sha(ROOT/name)==digest,'Changed verified source '+name)
        sourcecfg=manifest['identity']['config']
        require(sourcecfg['material']==self.cfg['material'] and sourcecfg['mu']==self.cfg['mu'],'Source geometry/material mismatch')
        require(sourcecfg['policy']==self.p,'Acceptance tolerances changed')
        for path in (self.source/'points').glob('*.json'):
            record=json.loads(path.read_text())
            case=physics.Case(**record['case']);validate_record(record,case,self.p)
            self.reference[key(case)]=(record,path)
        inventory=[]
        for construction in self.cfg['constructions']:
            angles=[b for b in self.cfg['angles_deg'] if key(self.case(construction,b)) in self.reference]
            inventory.append({'id':construction['id'],'available_angles':angles,
                              'missing_angles':[b for b in self.cfg['angles_deg'] if b not in angles]})
        write(self.folder/'source_audit.json',{'status':'PASS','inventory':inventory,
              'source_artifacts_checked':len(manifest['artifacts']),
              'source_status':'new-model gates PASS; previous overall PARTIAL_PASS due to unrelated historical CLI hash',
              'encoding_note':'Some pre-existing supplementary README/memory text contains question-mark replacements; numeric JSON/NPZ and canonical equations are valid; old text left unchanged.'})
        return sourcecfg

    def add_row(self,construction,beta,record,path,origin):
        case=self.case(construction,beta);validate_record(record,case,self.p)
        omega=np.array(record['omega']); lam=self.model.Lambda(omega)
        diags=record['diagnostics']
        row={'id':construction['id'],'kappa_j':case.k_joint,'kappa_s':case.kappa_horizontal,'beta_deg':beta,
             **{f'omega{i+1}':float(w) for i,w in enumerate(omega)},
             **{f'Omega{i+1}':float(v) for i,v in enumerate(lam**2)},
             **{f'Lambda{i+1}':float(v) for i,v in enumerate(lam)},
             'exact_zero_first':bool(omega[0]==0),'zero_count':int(sum(omega==0)),
             'F_omega':None,'D0_percent':None,'status':'PASS','origin':origin,
             'point_file':path.relative_to(ROOT).as_posix(),
             'boundary_max':max(max(np.max(d['joint_scaled']),np.max(d['external_scaled'])) for d in diags),
             'ode_max':max(d['ode_error'] for d in diags),'energy_max':max(d['energy_error'] for d in diags),
             'singular_max':max(d['singular_ratio'] for d in diags),'mass_gram_error':record['mass_gram_error'],
             'g12':float((omega[1]-omega[0])/omega[0]) if omega[0] else None}
        self.rows.append(row);self.records[(construction['id'],beta)]=record
        return row

    def neighbour_check(self,construction,previous,current):
        if previous is None:return
        jumps=[]
        for position in (1,2):
            a,b=previous[f'omega{position}'],current[f'omega{position}']
            jumps.append(None if a==0 else abs(b-a)/max(a,b))
        flags=[j for j in jumps if j is not None and j>self.cfg['frequency_map']['neighbour_jump_relative']]
        event={'id':construction['id'],'left_beta':previous['beta_deg'],'right_beta':current['beta_deg'],
               'relative_changes':jumps,'triggered':bool(flags),'status':'PASS'}
        if flags:
            case=self.case(construction,current['beta_deg'])
            record=self.records[(construction['id'],current['beta_deg'])]
            saved=[r for r in record.get('ritz',[]) if r['order']==20]
            if saved:
                rw=np.array(saved[0]['omega']); rd=saved[0];source='accepted archive Ritz p20'
            else:
                rw,_,rd=physics.ritz(self.model,case,20,rank_tol=self.p['kernel_rank_relative_tol'])
                self.ritz_count+=1;source='triggered independent Ritz p20 at existing angle'
            exact=np.array(record['omega']);positive=exact>0
            err=float(max(abs(rw[positive]/exact[positive]-1)))
            event.update(ritz_source=source,ritz_frequency_error=err,ritz_diagnostics=rd)
            require(err<=self.p['ritz_frequency_tol'],'Neighbour jump fails independent Ritz check')
            require(rd['nullity']==physics.expected_nullity(case),'Jump audit zero multiplicity')
            require(rd['mass_gram_error']<=self.p['mass_gram_tol'] and rd['energy_error']<=self.p['energy_tol'],'Jump audit Ritz energy/mass')
        self.audits.append(event)
        write(self.folder/'neighbour_audit.json',self.audits)

    def run(self):
        self.folder.mkdir(parents=True,exist_ok=True)
        if (self.folder/'manifest.json').exists():
            self.verify();print('REUSED completed HS-P1; zero new solves',self.folder);return
        write(self.folder/'contract.json',self.ident)
        basecfg=self.audit_source()
        engine=verified.Run(basecfg,self.folder/'computed')
        start=time.perf_counter()
        for construction in self.cfg['constructions']:
            previous=None
            for beta in self.cfg['angles_deg']:
                case=self.case(construction,beta)
                try:
                    match=self.reference.get(key(case))
                    target=engine.folder/'points'/f'{engine.name(case)}.json'
                    if match:
                        record,path=match;origin='verification_archive';self.reused+=1
                    elif target.exists() and json.loads(target.read_text()).get('status')=='PASS' and target.with_suffix('.npz').exists():
                        record=json.loads(target.read_text());path=target;origin='new_HS_P1';self.resumed+=1
                    else:
                        upper=self.p['Lambda_upper']
                        if previous is not None:
                            candidate=min(upper,float(self.model.Lambda(previous['omega3']*self.cfg['frequency_map']['guard_predictor_factor'])))
                            _,n,diag=physics.count_at(self.model,case,float(self.model.omega(candidate)),engine.catalogs,self.p)
                            used=candidate if n>=3 else upper
                            self.locators.append({'id':construction['id'],'beta':beta,'candidate_Lambda':candidate,
                                                  'candidate_count':n,'used_Lambda':used,'fallback':n<3,'count_diagnostic':diag})
                            upper=used
                        engine.p={**self.p,'Lambda_upper':upper}
                        record,_=engine.solve(case,'HS-P1',ritz=False)
                        record['locator_Lambda_upper']=upper
                        write(target,record)
                        path=target;origin='new_HS_P1';self.new+=1
                    row=self.add_row(construction,beta,record,path,origin)
                    self.neighbour_check(construction,previous,row)
                    previous=row
                    write(self.folder/'checkpoint.json',{'accepted':len(self.rows),'new_this_invocation':self.new,'archive_reused':self.reused,'resumed':self.resumed})
                except Exception as exc:
                    failure={'id':construction['id'],'beta':beta,'error':repr(exc),'series_stopped':True}
                    self.failures.append(failure);write(self.folder/'failures.json',self.failures)
                    print('SERIES STOP',failure,flush=True);break
        self.normalize()
        checks=self.monotonicity()
        csvwrite(self.folder/'spectra.csv',self.rows)
        csvwrite(self.folder/'summary.csv',self.summary())
        write(self.folder/'locator_audit.json',self.locators)
        write(self.folder/'monotonicity.json',checks)
        write(self.folder/'failures.json',self.failures)
        status='PASS' if len(self.rows)==105 and not self.failures and checks['status']=='PASS' else 'PARTIAL_PASS'
        manifest={'status':status,'identity':self.ident,'main_points':len(self.rows),
                  'new_distinct_points':sum(r['origin']=='new_HS_P1' for r in self.rows),
                  'archive_reused_points':self.reused,'new_this_invocation':self.new,'resumed_points':self.resumed,
                  'new_roots_this_invocation':sum(3-r['zero_count'] for r in self.rows if r['origin']=='new_HS_P1') if not self.resumed else None,
                  'triggered_ritz_solves':self.ritz_count,'added_angular_nodes':[],
                  'spectrum_seconds':time.perf_counter()-start,'frequency_map':self.cfg['frequency_map'],
                  'git_HEAD':json.loads((ROOT/self.cfg['initial_workspace']).read_text())['HEAD'],
                  'artifacts':{p.relative_to(self.folder).as_posix():sha(p) for p in self.folder.rglob('*') if p.is_file() and p.name!='manifest.json'}}
        write(self.folder/'manifest.json',manifest)
        print(status,self.folder,flush=True)

    def normalize(self):
        zero={r['id']:r['omega1'] for r in self.rows if r['beta_deg']==0}
        for r in self.rows:
            if zero.get(r['id'],0)>0:
                r['F_omega']=r['omega1']/zero[r['id']]
                r['D0_percent']=100*(r['F_omega']-1)

    def summary(self):
        result=[]
        for construction in self.cfg['constructions']:
            rows=sorted([r for r in self.rows if r['id']==construction['id']],key=lambda r:r['beta_deg'])
            if not rows:continue
            index=max(range(len(rows)),key=lambda i:rows[i]['omega1'])
            maximum=rows[index]
            result.append({'id':construction['id'],'kappa_j':construction['k_joint'],'kappa_s':construction['kappa_horizontal'],
                           'omega0':rows[0]['omega1'] if rows[0]['beta_deg']==0 else None,
                           'omega_max_grid':maximum['omega1'],'beta_max_grid':maximum['beta_deg'],
                           'max_left_node':rows[index-1]['beta_deg'] if index else None,
                           'max_right_node':rows[index+1]['beta_deg'] if index+1<len(rows) else None,
                           'omega90':rows[-1]['omega1'] if rows[-1]['beta_deg']==90 else None,
                           'Gmax_percent':maximum['D0_percent'],'nodes':len(rows),
                           'min_D0_percent':min(r['D0_percent'] for r in rows) if rows[0]['omega1']>0 else None,
                           'decreasing_steps':sum(b['omega1']<a['omega1']*(1-self.p['ritz_frequency_tol']) for a,b in zip(rows[:-1],rows[1:]))})
        return result

    def monotonicity(self):
        checks=[]
        for beta in self.cfg['angles_deg']:
            for axis,ids in [('external',['H2','H3','H4','H5','H6']),('joint',['H1','H5','H7'])]:
                rr={r['id']:r for r in self.rows if r['beta_deg']==beta}
                for low,high in zip(ids[:-1],ids[1:]):
                    if low not in rr or high not in rr:continue
                    delta=[rr[high][f'omega{i}']-rr[low][f'omega{i}'] for i in (1,2,3)]
                    ok=all(d>=-self.p['ritz_frequency_tol']*max(rr[high][f'omega{i}'],1e-10) for i,d in enumerate(delta,1))
                    checks.append({'beta':beta,'axis':axis,'low':low,'high':high,'differences':delta,'status':'PASS' if ok else 'FAIL'})
        return {'status':'PASS' if all(r['status']=='PASS' for r in checks) else 'FAIL','checks':checks}

    def verify(self):
        manifest=json.loads((self.folder/'manifest.json').read_text())
        require(manifest['identity']==self.ident,'Pilot source/config identity changed')
        for name,digest in manifest['artifacts'].items():
            require(sha(self.folder/name)==digest,'Changed pilot artifact '+name)
        return manifest

    def analyze(self):
        self.verify()
        import csv
        with (self.folder/'spectra.csv').open(encoding='utf8',newline='') as f:raw=list(csv.DictReader(f))
        rows=[{**r,**{k:float(r[k]) for k in ['beta_deg','omega1','omega2','omega3']}} for r in raw]
        plan=json.loads((self.folder/'analysis_plan.json').read_text(encoding='utf8'))
        start=time.perf_counter();energies=[];profiles=[]
        for selection in plan['shapes']:
            row=next(r for r in rows if r['id']==selection['id'] and r['beta_deg']==selection['beta_deg'])
            source=ROOT/row['point_file'];record=json.loads(source.read_text())
            case=physics.Case(**record['case']);mode=load_modes(source)[0]
            parts=energy_parts(self.model,case,mode,row['omega1'])
            require(abs(parts['mass']-1)<=self.p['mass_gram_tol'],'Selected form mass')
            require(parts['energy_error']<=self.p['energy_tol'],'Selected form energy')
            energies.append({'id':row['id'],'beta_deg':row['beta_deg'],'omega1':row['omega1'],'selection_reason':selection['reason'],**parts})
            for i,arm in enumerate(mode['arms'],1):
                for j,x in enumerate(arm['x']):
                    profiles.append({'id':row['id'],'beta_deg':row['beta_deg'],'arm':i,'x':x,
                                     **dict(zip(['u','w','psi','N','Q','M'],arm['y'][j]))})
        csvwrite(self.folder/'energy_distribution.csv',energies)
        csvwrite(self.folder/'selected_first_shapes.csv',profiles)
        source=json.loads((self.source/'gate4_asymptotics.json').read_text())['unloaded_prediction']
        C0,gamma,ca=source['C2_direct'],source['gamma0'],source['axial_dynamic_compliance']
        L2=self.model.l*(1+self.cfg['mu']);J2=self.model.section.rhoA*L2**3/3+self.model.section.rhoI*L2
        asym=[];coefficients=[]
        for construction in self.cfg['constructions']:
            name=construction['id']
            if name not in ['H1','H3','H4','H5','H6']:continue
            zero=next(r for r in rows if r['id']==name and r['beta_deg']==0)['omega1']
            kh=construction['kappa_horizontal']
            k=None if kh==physics.RIGID else kh*self.model.section.EI/self.model.l**3
            if name=='H1':
                coefficient=L2**2/J2/(1/k+2*self.model.l/self.model.section.EA)
            else:
                coefficient=C0+gamma**2*(1/ca if k is None else k/(1+k*ca))
            coefficients.append({'id':name,'omega0':zero,'coefficient_beta_rad_squared':coefficient,'source':'verified analytical coefficient; no fitting'})
            for angle in self.cfg['small_angle_nodes_deg']:
                row=next(r for r in rows if r['id']==name and r['beta_deg']==angle)
                b=math.radians(angle);actual=row['omega1']**2;approx=zero**2+coefficient*b*b
                resolution=4*max(row['omega1'],zero)*(self.p['root_xtol']+self.p['root_rtol']*max(row['omega1'],zero))
                asym.append({'id':name,'beta_deg':angle,'omega2_full':actual,'omega2_approx':approx,
                             'omega_full':row['omega1'],'omega_approx':math.sqrt(approx),
                             'signed_relative_error_lambda':actual/approx-1,
                             'absolute_error_lambda':abs(actual-approx),
                             'increment_resolution_indicator':resolution,
                             'increment_resolved_by_root_tolerance':abs(actual-zero**2)>20*resolution})
        csvwrite(self.folder/'small_angle_comparison.csv',asym)
        csvwrite(self.folder/'small_angle_coefficients.csv',coefficients)
        comparison=[]
        for beta in self.cfg['angles_deg']:
            a=next(r for r in rows if r['id']=='H5' and r['beta_deg']==beta)
            b=next(r for r in rows if r['id']=='H6' and r['beta_deg']==beta)
            comparison.append({'beta_deg':beta,'H5_omega1':a['omega1'],'H6_omega1':b['omega1'],'H5_below_H6_percent':100*(1-a['omega1']/b['omega1'])})
        csvwrite(self.folder/'H5_H6_comparison.csv',comparison)
        write(self.folder/'analysis_manifest.json',{'new_eigenvalue_solves':0,'seconds':time.perf_counter()-start,
              'source_manifest_sha256':sha(self.folder/'manifest.json'),'shape_count':len(energies),
              'artifacts':{name:sha(self.folder/name) for name in ['analysis_plan.json','energy_distribution.csv','selected_first_shapes.csv','small_angle_comparison.csv','small_angle_coefficients.csv','H5_H6_comparison.csv']}})
        print('Analysis complete; no new eigenvalue solves',self.folder)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage',choices=['run','analyze','plot_only','verify'],default='run')
    args=parser.parse_args()
    cfg=json.loads(CONFIG.read_text());pilot=Pilot(cfg)
    if args.stage=='run':pilot.run()
    elif args.stage=='analyze':pilot.analyze()
    elif args.stage=='verify':pilot.verify();print('PASS: archive verified, zero new solves',pilot.folder)
    else:
        from scripts.lib.two_springs_horizontal_pilot_figures import render
        render(pilot.folder)


if __name__=='__main__':main()
