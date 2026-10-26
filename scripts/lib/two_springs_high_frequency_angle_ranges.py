"""Discrete high-frequency ranges, reusing full spectra and the existing solver.

No physical assembly, branch tracking, leading solve or universal runner.
Refinement adds actual points only; all ranges refer to the resulting grid.
"""
import hashlib
import json
import math
from pathlib import Path
import shutil
import time
import numpy as np
import scipy
from scripts.lib import two_springs_external_transverse_spring as previous
from scripts.lib import two_springs_external_transverse_leading_resume as saved
from scripts.lib import isotropic_two_springs_external_transverse_spring as physics

ROOT=previous.ROOT
BASE=previous.BASE.parent/'high_frequency_angle_ranges'
CONFIG=ROOT/'data/input/isotropic_two_springs_high_frequency_angle_ranges.json'
REVISED_BASE=BASE.with_name('high_frequency_angle_ranges_90_95')
REVISED_CONFIG=CONFIG.with_name('isotropic_two_springs_high_frequency_angle_ranges_90_95.json')
api=previous.api
require=previous.require
read_csv=previous.read_csv
stiffness=previous.stiffness
FIGURES=('01_full_angle_frequency','02_frequency_deficit','03_high_range_boundaries','04_maximum_and_gain')


def high_ranges(beta,omega,level):
    """All contiguous runs of qualifying GRID nodes, including terminal status."""
    b=np.asarray(beta,float);w=np.asarray(omega,float)
    require(len(b)==len(w) and len(b)>=2 and b[0]==0 and b[-1]==90 and np.all(np.diff(b)>0),'Incomplete or unsorted 0--90 grid')
    require(np.all(np.isfinite(w)) and np.all(w>0) and 0<=level<100,'Invalid frequencies or threshold')
    R=100*(1-w/w.max());mask=w>=w.max()*(1-level/100);intervals=[]
    starts=np.flatnonzero(mask & np.r_[True,~mask[:-1]])
    ends=np.flatnonzero(mask & np.r_[~mask[1:],True])
    for i,j in zip(starts,ends):
        intervals.append(dict(start_deg=float(b[i]),end_deg=float(b[j]),nodes=int(j-i+1),
            entry_outside_deg=float(b[i-1]) if i else None,entry_resolution_deg=float(b[i]-b[i-1]) if i else None,
            exit_outside_deg=float(b[j+1]) if j+1<len(b) else None,exit_resolution_deg=float(b[j+1]-b[j]) if j+1<len(b) else None,
            singleton=bool(i==j),terminal=bool(j==len(b)-1)))
    terminal=intervals[-1] if mask[-1] else None
    status='NO_TERMINAL_HIGH_RANGE' if terminal is None else 'TERMINAL_POINT_ONLY' if terminal['singleton'] else 'TERMINAL_HIGH_RANGE'
    return dict(omega_max=float(w.max()),beta_max=float(b[np.argmax(w)]),G_max_percent=float(100*(w.max()/w[0]-1)),R=R,intervals=intervals,status=status,
                beta_high=terminal['start_deg'] if terminal else None,resolution=terminal['entry_resolution_deg'] if terminal else None,
                R90=float(R[-1]),first_high=float(b[starts[0]]),terminal_nodes=terminal['nodes'] if terminal else 0)


class Experiment(previous.Experiment):
    def __init__(self,revised=False):
        self.revised=revised
        self.a=api();self.cfg=json.loads((REVISED_CONFIG if revised else CONFIG).read_text());self.inputs={};self.sources={};self.manifests={};self.started=time.perf_counter();self.ritz_new=0;self.leading_new=0
        for name,digest in self.cfg['frozen_files'].items():self.check(ROOT/name,digest)
        for name,r in self.cfg['sources'].items():
            folder=ROOT/r['folder'];self.check(folder/'manifest.json',r['manifest_sha256']);saved.diagnostic.verify_source(folder)
            self.sources[name]=folder;self.manifests[name]=json.loads((folder/'manifest.json').read_text())
            require(self.manifests[name].get('full_model_status',self.manifests[name]['status'])=='PASS','Unaccepted source')
            cfg=json.loads(self.source(name,'contract.json').read_text())['config']
            require(cfg['material']==self.cfg['material'] and cfg.get('variant','RL_FULL')=='RL_FULL','Source material or axial theory mismatch')
            if name!='length':require(cfg['mu']==self.cfg['mu'],'Source geometry mismatch')
            if name=='prior':
                require(cfg['cases']==self.cfg['cases'] and cfg['k_w_values']==self.cfg['k_w_values'],'Source supports differ')
                require(cfg['map']['spectrum_semantics']=='sorted_positions' and cfg['map']['K_guard']==3,'Source spectral semantics')
        identity=dict(config=self.cfg,physical_sha256=self.a.sha(physics.__file__),versions=dict(numpy=np.__version__,scipy=scipy.__version__))
        self.key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16];self.folder=(REVISED_BASE if revised else BASE)/self.key;self.folder.mkdir(parents=True,exist_ok=True)
        target=self.folder/'contract.json'
        if target.exists():require(json.loads(target.read_text())==identity,'Immutable contract mismatch')
        else:self.a.write(target,identity)
        self.run=self.a.Run(self.cfg,self.folder,engine=physics);target=self.folder/'pole_catalogs.json'
        if revised:
            # Same material, lengths, physical solver and fixed-endpoint local
            # problems: the verified catalogs apply without a new pole search.
            for filename in ('pole_catalogs.json','triggered_counts.json'):
                source=self.source('prior',filename)
                if not (self.folder/filename).exists():shutil.copyfile(source,self.folder/filename)
        if target.exists():self.run.catalogs={float(k):v for k,v in json.loads(target.read_text()).items()}
        self.archive={};self.overlap=[];self.checked_files={}
        sources=[('prior','spectrum.csv')] if revised else [('transverse','spectrum.csv'),('replacement','spectrum.csv'),('length','complete_spectrum.csv'),('growth','spectrum.csv')]
        for source,fn in sources:
            for row in read_csv(self.source(source,fn)):
                pair=row.get('pair',row.get('ID'))
                if source in ('length','growth') and (pair not in ('J0E4','J4E0') or float(row.get('mu',.01))!=.01):continue
                self.ingest(ROOT/row['point_file'],row['point_sha256'],source,row)
        # Accepted finite-difference control points and repeated exact-support
        # preflights are legitimate existing spectra, not interpolated data.
        for source in (() if revised else ('transverse','replacement')):
            for fn,digest in self.manifests[source]['artifacts'].items():
                if fn.startswith('points/') and fn.endswith('.json'):self.ingest(self.sources[source]/fn,digest,source)
        inventory=[]
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                have=sorted(b for p,kw,b in self.archive if (p,kw)==(pair,k))
                if revised:require(have and have[-1]==90,'Source does not cover 90 degrees')
                require(have and have[0]==0 and all((pair,k,b) in self.archive for b in self.cfg['map']['parameter_grid'] if b<=20),'Incomplete accepted 0--20 data')
                missing=[b for b in self.cfg['map']['parameter_grid'] if b not in have]
                inventory.append(dict(pair=pair,k_w=k,available_count=len(have),available_beta_deg=have,missing_initial_beta_deg=missing,archive_max_beta=max(have)))
        self.a.write(self.folder/'archive_inventory.json',inventory);self.csv('archive_overlaps',self.overlap)
        if revised:require(len(self.archive)==self.cfg['expected_source_points'],'Source point inventory differs')
        self.allowed={};path=self.folder/'refinement_plan.json'
        self.plan=json.loads(path.read_text()) if path.exists() else []
        for r in self.plan:self.allowed.setdefault((r['pair'],r['k_w']),set()).update(r['new_beta_deg'])

    def ingest(self,path,digest,source,row=None):
        self.check(path,digest)
        if path not in self.checked_files:
            rec=json.loads(path.read_text());self.checked_files[path]=rec
        else:rec=self.checked_files[path]
        c=rec['case'];kw=c.get('k_w','RIGID');pair=next((p for p,v in self.cfg['cases'].items() if v==[c['k_joint'],c['k_end']]),None)
        if pair is None or c['mu']!=self.cfg['mu'] or kw not in self.cfg['k_w_values'] or c['variant']!='RL_FULL' or not 0<=c['beta']<=90:return
        self.check(path.with_suffix('.npz'),rec['shape_sha256'])
        cls=physics.Case if 'k_w' in c else physics.base.Case
        case=cls(**{n:v for n,v in c.items() if n!='variant'});w=self.validate(rec,case)
        if row is not None:
            require(all(abs(float(row[f'omega{i}'])/v-1)<1e-13 for i,v in enumerate(w,1)),'CSV/root mismatch')
        key=pair,kw,float(c['beta'])
        if key in self.archive:
            old=self.archive[key];err=float(max(abs(w/np.array([r['omega'] for r in old[0]['roots']])-1)))
            require(err<=self.cfg['policy']['same_physics_frequency_tol'],'Overlapping archives disagree')
            self.overlap.append(dict(pair=pair,k_w=kw,beta_deg=c['beta'],source=source,relative_frequency_error=err,same_file=path==old[1],status='PASS'))
        else:self.archive[key]=(rec,path,source)

    def point(self,pair,k,beta,previous_w=None):
        key=pair,k,float(beta)
        if key in self.archive:return self.archive[key]
        require((not self.revised and beta in self.cfg['map']['parameter_grid']) or beta in self.allowed.get((pair,k),set()),'Point lacks an approved grid/refinement reason')
        case=self.case(pair,k,beta);rec=self.run.solve(case,neighbour=previous_w);self.validate(rec,case)
        self.a.write(self.folder/'pole_catalogs.json',self.run.catalogs)
        return rec,previous.previous.lengths.point_path(self.folder,case),'new_study'

    def preflight(self):
        if self.revised:return self.revised_preflight()
        path=self.folder/'preflight.json'
        if path.exists():require(json.loads(path.read_text())['status']=='PASS','Previous preflight failure');return
        rows=[];ritz=[]
        for pair,k in self.cfg['preflight_cases']:
            for b in self.cfg['preflight_angles']:
                rec,file,source=self.point(pair,k,b);roots=rec['roots'];samples=rec['count']['samples']
                rows.append(dict(pair=pair,k_w=k,beta_deg=b,source=source,lower_count=samples[0]['count'],upper_count=samples[1]['count'],Lambda1=roots[0]['Lambda'],Lambda3=roots[2]['Lambda'],status='PASS'))
                if b in self.cfg['ritz_angles']:
                    for p in self.cfg['policy']['ritz_orders']:
                        w,coef,d=physics.ritz(self.run.model,self.case(pair,k,b),p);self.ritz_new+=1;err=float(max(abs(w/np.array([r['omega'] for r in roots])-1)))
                        ritz.append(dict(pair=pair,k_w=k,beta_deg=b,order=p,relative_frequency_error=err,**d));self.csv('preflight_ritz',ritz)
                        require(err<=self.cfg['policy']['ritz_frequency_tol'] and d['mass_gram_error']<=self.cfg['policy']['mass_gram_tol'] and d['energy_error']<=self.cfg['policy']['energy_tol'] and d['constraint_residual']<=self.cfg['policy']['boundary_tol'],'Large-angle independent Ritz gate')
        self.csv('preflight_points',rows)
        self.a.write(path,dict(status='PASS',points=len(rows),ritz_solves=len(ritz),max_ritz_error=max(r['relative_frequency_error'] for r in ritz),min_Lambda1=min(r['Lambda1'] for r in rows),max_Lambda3=max(r['Lambda3'] for r in rows)))

    def sweep(self):
        rows=[];quality=[];neighbours=[];audits_path=self.folder/'triggered_counts.json'
        audits=json.loads(audits_path.read_text()) if audits_path.exists() else {}
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                grid=sorted(set(self.cfg['map']['parameter_grid'])|{b for p,kw,b in self.archive if (p,kw)==(pair,k)}|self.allowed.get((pair,k),set()));prev=None
                for b in grid:
                    rec,file,source=self.point(pair,k,b,prev);w=np.array([r['omega'] for r in rec['roots']]);jump=float(max(abs(w/prev-1))) if prev is not None else 0.
                    triggered=jump>self.cfg['policy']['neighbour_relative_jump'];gap=w[1]/w[0]-1
                    if triggered or gap<self.cfg['refinement']['small_gap12']:
                        key=f'{pair}|{k}|{b}'
                        if key not in audits:
                            case=self.case(pair,k,b);cats=self.run.cats(self.run.model,case)
                            tests=[physics.count_at(self.run.model,case,x,cats,self.cfg['policy']) for x in (.9*w[0],(w[0]+w[1])/2,(w[1]+w[2])/2)]
                            require([r[1] for r in tests]==[0,1,2],'Triggered count failure');audits[key]=[r[2] for r in tests];self.a.write(audits_path,audits)
                    neighbours.append(dict(pair=pair,k_w=k,beta_deg=b,jump=jump,gap12=gap,triggered=triggered,status='PASS'));prev=w
                    rows.append(dict(pair=pair,variant=pair.split('_')[0],k_joint=self.cfg['cases'][pair][0],k_end=self.cfg['cases'][pair][1],k_w=k,beta_deg=b,**{f'{n}{i}':r[n] for i,r in enumerate(rec['roots'],1) for n in ('omega','Lambda')},gap12=gap,origin=source,point_file=file.relative_to(ROOT).as_posix(),point_sha256=self.a.sha(file),shape_sha256=rec['shape_sha256']))
                    for i,r in enumerate(rec['roots'],1):
                        d=r['diagnostics'];quality.append(dict(pair=pair,k_w=k,beta_deg=b,position=i,physical=float(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))),ode=d['ode_error'],energy=d['energy_error'],singular=d['singular_ratio'],mass_gram=rec['mass_gram_error'],condition=d['nonzero_condition'],count_lower=r['counts'][0],count_upper=r['counts'][1],low_count=rec['count']['samples'][0]['count'],upper_count=rec['count']['samples'][1]['count'],status='PASS'))
                for name,data in [('spectrum',rows),('quality',quality),('neighbours',neighbours)]:self.csv(name,data)
                print(f'{pair}, kw={k}: {len(grid)} checked nodes',flush=True)
        return rows

    def analyze(self,rows,prefix=''):
        summary=[];intervals=[];deficits=[]
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                rr=[r for r in rows if (r['pair'],r['k_w'])==(pair,k)];b=np.array([r['beta_deg'] for r in rr]);w=np.array([r['omega1'] for r in rr]);lam=np.array([r['Lambda1'] for r in rr])
                require(max(abs(w/w.max()-(lam/lam.max())**2))<2e-14,'Circular-frequency normalization')
                result=high_ranges(b,w,self.cfg['levels_percent'][0]);index=int(np.argmax(w))
                row=dict(pair=pair,variant=pair.split('_')[0],k_joint=self.cfg['cases'][pair][0],k_end=self.cfg['cases'][pair][1],k_w=k,nodes=len(rr),Lambda0=lam[0],Lambda_max=lam[index],omega0=w[0],omega_max=w[index],beta_max=b[index],G_max_percent=result['G_max_percent'],R90=result['R90'],min_gap12=min(r['gap12'] for r in rr))
                for level in self.cfg['levels_percent']:
                    result=high_ranges(b,w,level);tag=int(level)
                    row.update({f'beta_high_{tag}':result['beta_high'],f'status_{tag}':result['status'],f'resolution_{tag}':result['resolution'],f'first_high_{tag}':result['first_high'],f'terminal_nodes_{tag}':result['terminal_nodes']})
                    intervals.extend(dict(pair=pair,k_w=k,delta_percent=level,**r) for r in result['intervals'])
                summary.append(row)
                deficits.extend(dict(pair=pair,k_w=k,beta_deg=float(x),omega1=float(v),Lambda1=float(l),R_percent=float(R),omega_max=w[index],**{f'qualifies_{int(d)}':bool(v>=w[index]*(1-d/100)) for d in self.cfg['levels_percent']}) for x,v,l,R in zip(b,w,lam,result['R']))
        self.csv(prefix+'summary',summary);self.csv(prefix+'high_intervals',intervals);self.csv(prefix+'deficit',deficits)
        return summary

    def refine(self,rows,round_number):
        cfg=self.cfg['refinement'];new=[]
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                rr=[r for r in rows if (r['pair'],r['k_w'])==(pair,k)];b=np.array([r['beta_deg'] for r in rr]);w=np.array([r['omega1'] for r in rr]);reasons=[]
                for level in cfg.get('levels_percent',self.cfg['levels_percent']):
                    mask=w>=w.max()*(1-level/100)
                    for i in np.flatnonzero(mask[1:]!=mask[:-1]):reasons.append((b[i],b[i+1],cfg['crossing_step_deg'],f'crossing_{level:g}_percent'))
                if not self.revised:
                    i=int(np.argmax(w));lo=b[max(0,i-1)];hi=b[min(len(b)-1,i+1)]
                    reasons.append((lo,hi,cfg['maximum_step_deg'],'sampled_maximum'))
                    for i in range(1,len(b)-1):
                        left=w[i]-w[i-1];right=w[i+1]-w[i]
                        if left*right<0 and min(abs(left),abs(right))/w[i]>cfg['extremum_frequency_rtol']:
                            reasons.append((b[i-1],b[i+1],cfg['extremum_step_deg'],'resolved_local_extremum'))
                    gaps=np.array([r['gap12'] for r in rr]);i=int(np.argmin(gaps))
                    if gaps[i]<cfg['small_gap12']:reasons.append((b[max(i-1,0)],b[min(i+1,len(b)-1)],cfg['crossing_step_deg'],'small_first_gap'))
                    for i in range(len(b)-1):
                        if max(abs(rr[i+1][f'omega{j}']/rr[i][f'omega{j}']-1) for j in (1,2))>self.cfg['policy']['neighbour_relative_jump']:
                            reasons.append((b[i],b[i+1],cfg['extremum_step_deg'],'neighbour_jump'))
                for lo,hi,h,reason in reasons:
                    if hi-lo<=h+1e-10:continue
                    add=[round(float(x),10) for x in np.arange(math.ceil((lo+1e-9)/h)*h,hi-1e-9,h) if min(abs(b-x))>1e-8]
                    if not add:continue
                    r=dict(round=round_number,pair=pair,k_w=k,reason=reason,left_deg=float(lo),right_deg=float(hi),target_step_deg=h,new_beta_deg=add)
                    self.plan.append(r);new.append(r);self.allowed.setdefault((pair,k),set()).update(add)
        self.a.write(self.folder/'refinement_plan.json',self.plan)
        return new

    def run_study(self):
        self.preflight();rows=self.sweep();self.analyze(rows,'initial_') if not (self.folder/'initial_summary.csv').exists() else None
        for round_number in range(1,self.cfg['refinement']['max_rounds']+1):
            additions=self.refine(rows,round_number)
            if not additions:break
            rows=self.sweep()
        else:require(not self.refine(rows,round_number+1),'Refinement budget reached with unresolved brackets')
        summary=self.analyze(rows)
        for r in summary:
            for level in map(int,self.cfg['levels_percent']):require(r[f'resolution_{level}'] is None or r[f'resolution_{level}']<=.25+1e-10,'Terminal boundary unresolved')
        if self.revised:self.revised_comparison(rows,summary)
        self.a.write(self.folder/'calculation_complete.json',dict(status='PASS',constructions=len(summary),points=len(rows),reused_points=sum(r['origin']!='new_study' for r in rows),new_points=sum(r['origin']=='new_study' for r in rows),refinement_rounds=max((r['round'] for r in self.plan),default=0),max_delta_levels=self.cfg['levels_percent']))

    def finish(self,stage,status='PASS'):
        super().finish(stage,status)
        path=self.folder/'executions.json';runs=json.loads(path.read_text());runs[-1]['workflow_sha256']=self.a.sha(__file__);self.a.write(path,runs)

    def revised_preflight(self):
        """Read-only reproduction of the old maxima and both historical levels."""
        path=self.folder/'preflight.json'
        if path.exists():require(json.loads(path.read_text())['status']=='PASS','Previous regression failed');return
        with saved.forbid_full_calculations(leading_too=True):rows=self.sweep()
        old=read_csv(self.source('prior','summary.csv'));regression=[]
        require(len(old)==20,'Historical summary size')
        for r in old:
            rr=[v for v in rows if (v['pair'],v['k_w'])==(r['pair'],stiffness(r['k_w']))]
            b=[v['beta_deg'] for v in rr];w=[v['omega1'] for v in rr]
            for d in (2,5):
                q=high_ranges(b,w,d)
                require(q['omega_max']==float(r['omega_max']) and q['beta_max']==float(r['beta_max']),'Historical maximum differs')
                require(q['status']==r[f'status_{d}'] and q['beta_high']==(float(r[f'beta_high_{d}']) if r[f'beta_high_{d}'] else None),'Historical boundary differs')
                regression.append(dict(pair=r['pair'],k_w=r['k_w'],delta_percent=d,beta_high=q['beta_high'],omega_max=q['omega_max'],maximum_error=0.,boundary_error_deg=0.,status='PASS'))
        summary=self.analyze(rows,'initial_');self.check_nested(summary)
        self.csv('historical_regression',regression)
        self.a.write(path,dict(status='PASS',archive_points=len(rows),historical_boundaries=len(regression),new_full_prefixes=self.run.new_points,new_roots=self.run.new_roots))

    @staticmethod
    def check_nested(summary):
        for r in summary:
            if r['beta_high_5'] is not None:
                require(r['beta_high_10'] is not None and r['beta_high_10']<=r['beta_high_5'],'Non-nested high ranges')

    def revised_comparison(self,rows,summary):
        self.check_nested(summary);old=read_csv(self.source('prior','summary.csv'));comparisons=[];boundaries=[];historical=[]
        transitions=read_csv(self.source('length','transitions.csv'))
        for r in summary:
            prior=next(x for x in old if (x['pair'],stiffness(x['k_w']))==(r['pair'],r['k_w']))
            rr=[v for v in rows if (v['pair'],v['k_w'])==(r['pair'],r['k_w'])];b=[v['beta_deg'] for v in rr];w=[v['omega1'] for v in rr]
            old2=high_ranges(b,w,2)
            change=r['omega_max']/float(prior['omega_max'])-1
            comparisons.append(dict(pair=r['pair'],k_w=r['k_w'],omega_max_old=float(prior['omega_max']),omega_max_new=r['omega_max'],maximum_relative_change=change,
                beta_high_5_old=float(prior['beta_high_5']),beta_high_5_new=r['beta_high_5'],beta_high_2_old=float(prior['beta_high_2']),beta_high_2_new=old2['beta_high'],
                reason='same maximum and original nodes retained; added 10-percent bracket nodes' if change==0 else 'new actual spectral point raised the sampled maximum'))
            for d in (10,5):
                q=high_ranges(b,w,d);tail=next((i for i in q['intervals'] if i['terminal']),None)
                boundaries.append(dict(pair=r['pair'],k_w=r['k_w'],delta_percent=d,status=q['status'],beta_high=q['beta_high'],entry_outside_deg=tail['entry_outside_deg'] if tail else None,
                    resolution_deg=q['resolution'],starts_at_zero=q['beta_high']==0,checked_tail_nodes=q['terminal_nodes'],tail_beta_deg=[v for v in b if q['beta_high'] is not None and v>=q['beta_high']],max_tail_deficit_percent=max([R for v,R in zip(b,q['R']) if q['beta_high'] is not None and v>=q['beta_high']],default=None)))
            if r['pair'] in ('I_J01','I_J10') and r['k_w']=='RIGID':
                pair='J0E4' if r['pair']=='I_J01' else 'J4E0'
                t=next(t for t in transitions if t['pair']==pair and float(t['mu'])==self.cfg['mu'])
                require(t['status']=='PASS' and t['tail_to_90']=='True','Unverified historical derivative criterion')
                at_pl=next(v['omega1'] for v in rr if v['beta_deg']==float(t['beta_pl_deg']))
                historical.append(dict(pair=r['pair'],k_joint=r['k_joint'],k_end=r['k_end'],beta_pl_deg=float(t['beta_pl_deg']),q_threshold=.15,R_at_beta_pl_percent=100*(1-at_pl/r['omega_max']),beta_high_10=r['beta_high_10'],beta_high_5=r['beta_high_5'],beta_high_2_historical=float(prior['beta_high_2']),source=self.source('length','transitions.csv').relative_to(ROOT).as_posix()))
        self.csv('historical_comparison',comparisons);self.csv('derivative_criterion_comparison',historical);self.a.write(self.folder/'boundary_evidence.json',boundaries)


def postprocess_saved(folder):
    a=api();cfg=json.loads((folder/'contract.json').read_text())['config'];rows=read_csv(folder/'spectrum.csv');summary=read_csv(folder/'summary.csv')
    levels=list(map(int,cfg['levels_percent']))
    initial={(r['pair'],stiffness(r['k_w'])):r for r in read_csv(folder/'initial_summary.csv')};effects=[];changes=[];checks=[]
    by={(r['pair'],stiffness(r['k_w']),float(r['beta_deg'])):r for r in rows}
    for r in summary:
        old=initial[r['pair'],stiffness(r['k_w'])]
        effects.append(dict(pair=r['pair'],k_w=r['k_w'],initial_nodes=old['nodes'],final_nodes=r['nodes'],initial_beta_max=old['beta_max'],final_beta_max=r['beta_max'],relative_maximum_change=float(r['omega_max'])/float(old['omega_max'])-1,
            **{f'{prefix}_beta_high_{d}':data[f'beta_high_{d}'] for d in levels for prefix,data in [('initial',old),('final',r)]}))
    for j in ('J01','J10'):
        for k in cfg['k_w_values']:
            old=next(r for r in summary if r['pair']=='I_'+j and stiffness(r['k_w'])==k);new=next(r for r in summary if r['pair']=='II_'+j and stiffness(r['k_w'])==k)
            changes.append(dict(k_joint=old['k_joint'],k_w=k,omega_max_I=old['omega_max'],omega_max_II=new['omega_max'],omega_max_change_percent=100*(float(new['omega_max'])/float(old['omega_max'])-1),G_max_I=old['G_max_percent'],G_max_II=new['G_max_percent'],
                **{f'beta_high_{d}_change_deg':float(new[f'beta_high_{d}'])-float(old[f'beta_high_{d}']) if old[f'beta_high_{d}'] and new[f'beta_high_{d}'] else None for d in levels}))
            for pair1,pair2,kw1,kw2,kind in [('II_'+j,'I_'+j,k,k,'release_rotation')]:
                grid={b for p,kw,b in by if (p,kw)==(pair1,kw1)} & {b for p,kw,b in by if (p,kw)==(pair2,kw2)}
                for b in sorted(grid):
                    diff=max(float(by[pair1,kw1,b][f'omega{i}'])/float(by[pair2,kw2,b][f'omega{i}'])-1 for i in (1,2,3))
                    require(diff<=cfg['policy']['same_physics_frequency_tol'],'Absolute frequency increased after release')
                    checks.append(dict(kind=kind,pair=pair1,k_w_lower=k,k_w_upper=k,beta_deg=b,max_violation=diff,status='PASS'))
    for pair in cfg['cases']:
        for lo,hi in zip([0.,.1,1.,10.],[.1,1.,10.,'RIGID']):
            grid={b for p,kw,b in by if (p,kw)==(pair,lo)} & {b for p,kw,b in by if (p,kw)==(pair,hi)}
            for b in sorted(grid):
                diff=max(float(by[pair,lo,b][f'omega{i}'])/float(by[pair,hi,b][f'omega{i}'])-1 for i in (1,2,3))
                require(diff<=cfg['policy']['same_physics_frequency_tol'],'Absolute stiffness monotonicity')
                checks.append(dict(kind='increase_k_w',pair=pair,k_w_lower=lo,k_w_upper=hi,beta_deg=b,max_violation=diff,status='PASS'))
    a.csvwrite(folder/'refinement_effects.csv',effects);a.csvwrite(folder/'rotation_comparison.csv',changes);a.csvwrite(folder/'monotonicity_checks.csv',checks)


def render(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    cfg=json.loads((folder/'contract.json').read_text())['config'];full=read_csv(folder/'spectrum.csv');deficit=read_csv(folder/'deficit.csv');summary=read_csv(folder/'summary.csv');intervals=read_csv(folder/'high_intervals.csv')
    revised=10 in cfg['levels_percent'];levels=list(map(int,cfg['levels_percent']));reuse=[]
    if revised:
        comparisons=read_csv(folder/'historical_comparison.csv')
        if all(float(r['maximum_relative_change'])==0 for r in comparisons):reuse=[0,3]
    plt.rcParams.update({'font.size':10,'axes.labelsize':11,'legend.fontsize':8,'pdf.fonttype':42,'font.family':'DejaVu Sans'})
    colors=['#222222','#0072B2','#D55E00','#009E73','#8B4A9C'];styles=['-','--','-.',':',(0,(5,1,1,1))]
    comma=FuncFormatter(lambda x,pos:f'{x:g}'.replace('.',','));kwlabels=['Точная опора',r'$\kappa_w=10$',r'$\kappa_w=1$',r'$\kappa_w=0{,}1$',r'$\kappa_w=0$']
    def group(rows,pair,k):return [r for r in rows if r['pair']==pair and stiffness(r['k_w'])==k]
    def v(rows,key):return [float(r[key]) for r in rows]
    def num(x):return f'{float(x):g}'.replace('.',r'{,}')
    def title(pair):
        j,e=cfg['cases'][pair];return r'$\kappa_{\mathrm{j}}='+num(j)+r',\quad\kappa_{\mathrm{e}}='+num(e)+'$'
    def save(fig,axes,index):
        for ax in np.asarray(axes).flat:ax.grid(alpha=.2);ax.yaxis.set_major_formatter(comma)
        if index in reuse:
            for suffix in ('.pdf','.png'):shutil.copyfile(ROOT/cfg['sources']['prior']['folder']/(FIGURES[index]+suffix),folder/(FIGURES[index]+suffix))
        else:fig.savefig(folder/(FIGURES[index]+'.pdf'));fig.savefig(folder/(FIGURES[index]+'.png'),dpi=320)
        plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(7.8,6.0),sharex=True,sharey=True,layout='constrained')
    for ax,pair,letter in zip(axes.flat,cfg['cases'],'абвг'):
        for i,k in enumerate(cfg['k_w_values']):
            rr=group(full,pair,k);ax.plot(v(rr,'beta_deg'),v(rr,'Lambda1'),color=colors[i],ls=styles[i],lw=1.4,label=kwlabels[i])
        ax.set_title(letter+') '+title(pair),fontsize=10);ax.set_xlim(0,90);ax.xaxis.set_major_formatter(comma)
    for ax in axes[:,0]:ax.set_ylabel(r'$\Lambda_1$')
    for ax in axes[1]:ax.set_xlabel(r'$\beta$, град.')
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False);save(fig,axes,0)
    # Select observed early and late attainment, plus the matching rotation
    # comparison; no a priori label of a construction as delayed growth.
    selection_level=10 if revised else 2
    earliest=min(summary,key=lambda r:float(r[f'beta_high_{selection_level}']) if r[f'beta_high_{selection_level}'] else math.inf)
    latest=max(summary,key=lambda r:float(r[f'beta_high_{selection_level}']) if r[f'beta_high_{selection_level}'] else -math.inf)
    pair='I_'+earliest['pair'].split('_',1)[1] if earliest['variant']=='II' else 'II_'+earliest['pair'].split('_',1)[1]
    counterpart=next(r for r in summary if r['pair']==pair and r['k_w']==earliest['k_w'])
    finite=max([r for r in summary if r['k_w']!='RIGID'],key=lambda r:float(r[f'beta_high_{selection_level}']) if r[f'beta_high_{selection_level}'] else -math.inf)
    selected=[earliest,counterpart,finite,latest]
    api().csvwrite(folder/'figure_selection.csv',[dict(panel=i+1,pair=r['pair'],k_w=r['k_w'],reason=reason) for i,(r,reason) in enumerate(zip(selected,[f'earliest_{100-selection_level}_percent','same_transverse_support_other_rotation',f'latest_finite_support_{100-selection_level}_percent',f'latest_overall_{100-selection_level}_percent']))])
    fig,axes=plt.subplots(2,2,figsize=(7.8,6.1),sharex=True,sharey=True,layout='constrained')
    for ax,r,letter in zip(axes.flat,selected,'абвг'):
        pair,k=r['pair'],stiffness(r['k_w']);rr=group(deficit,pair,k)
        ax.plot(v(rr,'beta_deg'),v(rr,'R_percent'),color='#222222',lw=1.4)
        for level,color in [(max(levels),'#D55E00'),(min(levels),'#0072B2')]:
            ax.axhline(level,color=color,ls='--',lw=.9)
            for interval in intervals:
                if interval['pair']==pair and stiffness(interval['k_w'])==k and float(interval['delta_percent'])==level:
                    ax.axvspan(float(interval['start_deg']),float(interval['end_deg']),color=color,alpha=.06)
            if r[f'beta_high_{level}']:
                x=float(r[f'beta_high_{level}']);point=next(q for q in rr if float(q['beta_deg'])==x)
                ax.plot([x],[float(point['R_percent'])],marker='*' if revised and x==0 else 'o' if level==min(levels) else 's',mfc='white',mec=color,ms=8 if revised and x==0 else 4,clip_on=False)
        ranges='\n'.join(f'{100-level if revised else level}%: [{r[f"beta_high_{level}"].replace(".",",")}; 90°]' if r[f'beta_high_{level}'] else f'{level}%: нет хвостового диапазона' for level in sorted(levels,reverse=True))
        ax.text(.46,.72,ranges,transform=ax.transAxes,fontsize=9,va='top')
        kw='точная опора' if k=='RIGID' else r'$\kappa_w='+num(k)+'$'
        ax.set_title(letter+') '+title(pair)+', '+kw,fontsize=9);ax.set_xlim(0,12 if revised else 20);ax.set_ylim(-.3,20 if revised else 10);ax.xaxis.set_major_formatter(comma)
    for ax in axes[:,0]:ax.set_ylabel(r'$R$, %')
    for ax in axes[1]:ax.set_xlabel(r'$\beta$, град.')
    save(fig,axes,1)
    order=[0.,.1,1.,10.,'RIGID'];x=np.arange(5);ticks=['0','0,1','1','10','Точная\nопора']
    fig,axes=plt.subplots(1,2,figsize=(7.8,3.9),sharey=True,layout='constrained')
    for ax,j in zip(axes,['J01','J10']):
        for variant,color,marker in [('I','#0072B2','o'),('II','#D55E00','s')]:
            rr={stiffness(r['k_w']):r for r in summary if r['pair']==variant+'_'+j}
            for level,ls in zip(sorted(levels),['-','--']):
                y=[float(rr[k][f'beta_high_{level}']) if rr[k][f'beta_high_{level}'] else np.nan for k in order]
                style=dict(ms=6 if variant=='I' else 3.5,mfc='white' if variant=='I' else color,lw=2.5 if variant=='I' else 1.1,zorder=3 if variant=='I' else 4) if revised else dict(ms=4)
                ax.plot(x,y,color=color,marker=marker,ls=ls,label=f'{variant}, {100-level if revised else level}%',**style)
                if revised:
                    for i,value in enumerate(y):
                        if value==0:ax.plot(i,0,marker='*',ms=10,mfc='white',mec='black',ls='None',zorder=5,label='Начало при 0°')
                for i,k in enumerate(order):
                    if rr[k][f'status_{level}']=='NO_TERMINAL_HIGH_RANGE':ax.plot(i,1.03,marker='x',color=color,transform=ax.get_xaxis_transform(),clip_on=False)
        ax.set_xticks(x,ticks);ax.set_xlabel(r'$\kappa_w$');ax.set_title(r'$\kappa_{\mathrm{j}}='+num(cfg['cases']['I_'+j][0])+'$',fontsize=10)
    axes[0].set_ylabel(r'$\beta_{\mathrm{high}}^{\mathrm{grid}}$, град.');fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=3 if revised else 4,frameon=False);save(fig,axes,2)
    fig,axes=plt.subplots(2,2,figsize=(7.8,5.8),sharex=True,sharey='row',layout='constrained')
    for col,j in enumerate(['J01','J10']):
        for variant,color,marker,ls in [('I','#0072B2','o','-'),('II','#D55E00','s','--')]:
            rr={stiffness(r['k_w']):r for r in summary if r['pair']==variant+'_'+j}
            for row,key in enumerate(['omega_max','G_max_percent']):axes[row,col].plot(x,[float(rr[k][key]) for k in order],color=color,marker=marker,ms=4,ls=ls,label='Поворотная пружина' if variant=='I' else 'Свободный поворот')
        axes[0,col].set_title(r'$\kappa_{\mathrm{j}}='+num(cfg['cases']['I_'+j][0])+'$',fontsize=10);axes[1,col].set_xticks(x,ticks);axes[1,col].set_xlabel(r'$\kappa_w$')
    axes[0,0].set_ylabel(r'$\omega_{\max}^{\mathrm{grid}}$');axes[1,0].set_ylabel(r'$G_{\max}^{\mathrm{grid}}$, %')
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=2,frameon=False);save(fig,axes,3)
    if revised:
        api().write(folder/'figure_provenance.json',[dict(figure=i+1,reused=i in reuse,source=(Path(cfg['sources']['prior']['folder'])/(name+'.pdf')).as_posix() if i in reuse else 'current CSV',data_source='prior/spectrum.csv' if i==0 and i in reuse else 'summary.csv' if i in (2,3) else 'deficit.csv',reason='same physical frequencies and sampled maxima; new points only refine existing crossings' if i in reuse else 'new primary levels 90 and 95 percent') for i,name in enumerate(FIGURES)])


def finalize(folder):
    a=api();cfg=json.loads((folder/'contract.json').read_text())['config'];complete=json.loads((folder/'calculation_complete.json').read_text())
    require(complete['status']=='PASS' and not (folder/'failure.json').exists(),'Unaccepted study')
    for name in FIGURES:
        for suffix in ('.pdf','.png'):require((folder/(name+suffix)).is_file(),'Missing figure')
    for name,digest in json.loads((folder/'source_hashes.json').read_text()).items():require(a.sha(ROOT/name)==digest,'Protected source changed '+name)
    a.write(folder/'source_integrity.json',{name:saved.diagnostic.verify_source(ROOT/r['folder']) for name,r in cfg['sources'].items()})
    runs=json.loads((folder/'executions.json').read_text());q=read_csv(folder/'quality.csv');summary=read_csv(folder/'summary.csv');overlap=read_csv(folder/'archive_overlaps.csv') if (folder/'archive_overlaps.csv').exists() else []
    decision=dict(**complete,new_full_prefixes=sum(r['new_full_prefixes'] for r in runs),new_full_roots=sum(r['new_full_roots'] for r in runs),independent_ritz_solves=sum(r['new_ritz_solves'] for r in runs),new_leading_solves=0,quality_maxima={n:max(float(r[n]) for r in q) for n in ('physical','ode','energy','singular','mass_gram','condition')},overlap_count=len(overlap),max_overlap_error=max((float(r['relative_frequency_error']) for r in overlap),default=0.),no_terminal_cases={str(d):sum(r[f'status_{d}']=='NO_TERMINAL_HIGH_RANGE' for r in summary) for d in map(int,cfg['levels_percent'])})
    if 10 in cfg['levels_percent']:
        rows=read_csv(folder/'spectrum.csv');keys={(r['pair'],r['k_w'],r['beta_deg']) for r in rows if r['origin']=='new_study'}
        newq=[r for r in q if (r['pair'],r['k_w'],r['beta_deg']) in keys];comparison=read_csv(folder/'historical_comparison.csv')
        require(all(r['reason']=='crossing_10_percent' for r in json.loads((folder/'refinement_plan.json').read_text())),'Refinement exceeded the revised scope')
        require(all(not r['beta_high_5'] or r['beta_high_10'] and float(r['beta_high_10'])<=float(r['beta_high_5']) for r in summary),'Range nesting failed')
        decision.update(new_quality_maxima={n:max((float(r[n]) for r in newq),default=0.) for n in ('physical','ode','energy','singular','mass_gram','condition')},
            zero_start_cases=[dict(pair=r['pair'],k_w=r['k_w']) for r in summary if r['beta_high_10'] and float(r['beta_high_10'])==0],
            maxima_changed=sum(float(r['maximum_relative_change'])!=0 for r in comparison),historical_5_changed=sum(r['beta_high_5_old']!=r['beta_high_5_new'] for r in comparison),
            historical_2_changed=sum(r['beta_high_2_old']!=r['beta_high_2_new'] for r in comparison),new_pole_roots=sum(r['new_pole_roots'] for r in runs),range_nesting='PASS')
    require(decision['new_full_prefixes']==decision['new_points'],'New point accounting mismatch');a.write(folder/'decision.json',decision)
    a.write(folder/'manifest.json',dict(status='PASS',artifacts={p.relative_to(folder).as_posix():a.sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    return saved.diagnostic.verify_source(folder)


def main(args):
    revised=getattr(args,'high_frequency_ranges_90_95',False)
    if args.stage in ('plot_only','finalize'):
        base=REVISED_BASE if revised else BASE
        require(args.folder is not None and args.folder.resolve().is_relative_to(base.resolve()),'Explicit new archive required')
        with saved.forbid_full_calculations(leading_too=True):
            if args.stage=='plot_only':
                started=time.perf_counter();postprocess_saved(args.folder);render(args.folder)
                api().write(args.folder/'rendering.json',dict(seconds=time.perf_counter()-started,new_spectral_solves=0,dpi=320))
            else:print(finalize(args.folder))
        return
    require(args.stage in ('preflight','run'),'Supported stages: preflight/run/plot_only/finalize')
    exp=Experiment(revised=revised);print(exp.folder,flush=True)
    try:
        require(not (exp.folder/'failure.json').exists(),'Saved failure requires addressed diagnosis')
        if args.stage=='preflight':exp.preflight()
        else:exp.run_study()
        exp.finish(args.stage)
    except Exception as exc:
        exp.a.write(exp.folder/'failure.json',dict(status='PARTIAL',stage=args.stage,error=repr(exc)));exp.finish(args.stage,'PARTIAL');raise
