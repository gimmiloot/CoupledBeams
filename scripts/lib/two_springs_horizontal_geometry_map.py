"""HS-B geometry contract around the accepted HS-A transaction workflow.

No physical assembly, counter or root solver is reimplemented. Geometry-specific
arm signatures and independent stopping/resume scopes cannot be put into the
frozen HS-A identity (its Run constructor explicitly requires the old geometry).
"""
from pathlib import Path
import argparse
import copy
import csv
import hashlib
import json
import math
import platform
import time

import numpy as np
import scipy
from scipy.optimize import brentq
from scripts.lib import two_springs_horizontal_angle_map as a
from scripts.analysis.screen_coupled_longitudinal_theory_hierarchy import tim_catalog
from scripts.lib.isotropic_two_springs_horizontal_leading import Leading, independent_roots
from scripts.lib.two_springs_horizontal_right_angle import RightAngle, independent_checks

ROOT=a.ROOT
CONFIG=ROOT/'data/input/isotropic_two_springs_horizontal_geometry_map.json'
OUTPUT=ROOT/'results/isotropic_two_springs_rl_timoshenko/horizontal_spring_geometry_map'
read,atomic,sha,require,csvwrite=a.read,a.atomic,a.sha,a.require,a.csvwrite


def variants(cfg):
    out=[]
    def add(regime,group,suffix='',epsilon=.005,mu=.01,aspect=4):
        chi=regime['chi_s']
        out.append(dict(id=f'{group}_{regime["id"]}{suffix}',regime=regime['id'],group=group,
                        epsilon=epsilon,mu=mu,aspect=aspect,k_joint=regime['k_joint'],chi_s=chi,
                        kappa_horizontal=a.h.RIGID if chi==a.h.RIGID else chi/epsilon**2))
    for r in cfg['regimes']:
        add(r,'B0')
        for e in cfg['variations']['B1_epsilon']:add(r,'B1',f'_e{e:g}',epsilon=e)
        for m in cfg['variations']['B2_mu']:add(r,'B2',f'_m{m:g}',mu=m)
        for s in cfg['variations']['B3_aspect']:add(r,'B3',f'_a{s:g}',aspect=s)
        v=cfg['variations']['B4']
        if r['id'] in v['regimes']:add(r,'B4',epsilon=v['epsilon'],mu=v['mu'],aspect=v['aspect'])
    require(len(out)==34 and len({v['id'] for v in out})==34,'Variant contract')
    return out


def geometry(v,cfg):
    model=a.h.Model(**{**cfg['material'],'epsilon':v['epsilon'],'aspect':v['aspect']})
    s=model.section;h=math.sqrt(12)*v['epsilon']*model.l;b=v['aspect']*h
    low=model.l*(1-abs(v['mu']))
    return model,dict(**v,h=h,b=b,L1=model.l*(1-v['mu']),L2=model.l*(1+v['mu']),
        h_over_Lmin=h/low,b_over_Lmin=b/low,A=b*h,I_y=b*h**3/12,
        EA=s.EA,EI=s.EI,KGA=s.KGA,rhoA=s.rhoA,rhoI=s.rhoI,J_RL=model.J,
        k_joint_dimensional=v['k_joint']*s.EI/model.l,
        k_horizontal_dimensional=None if v['chi_s']==a.h.RIGID else v['chi_s']*s.EA/model.l,
        limiter=max(b,h)/low,status='PASS' if max(b,h)/low<=cfg['geometry_limiter'] else 'GEOMETRY_STOP')


def arm_signature(model,L):
    s=model.section
    return dict(length=L,EA=s.EA,EI=s.EI,KGA=s.KGA,rhoA=s.rhoA,rhoI=s.rhoI,J_RL=model.J)


class GeometryStudy(a.Study):
    def __init__(self,parent,v):
        self.parent=parent;self.v=v
        self.cfg=copy.deepcopy(parent.cfg);self.cfg['mu']=v['mu']
        self.cfg['material'].update(epsilon=v['epsilon'],aspect=v['aspect'])
        self.p=self.cfg['policy'];self.folder=parent.folder/'variants'/v['id']
        self.folder.mkdir(parents=True,exist_ok=True)
        self.model=geometry(v,self.cfg)[0];self.constructions=[v];self.reference={}
        if v['group']=='B0':
            oldid=a.construction_id(v['k_joint'],v['kappa_horizontal'])
            for e in parent.baseline_entries.values():
                if e['id']==oldid:
                    rec=read(ROOT/e['path']);a.pilot.validate_record(rec,self.case(v,e['beta_deg']),self.p)
                    self.reference[a.pilot.key(self.case(v,e['beta_deg']))]=(rec,ROOT/e['path'])
        # Dependency injection replaces only the old constructor's catalog lookup.
        self.engine=a.verified.Run.__new__(a.verified.Run)
        self.engine.cfg=self.cfg;self.engine.folder=self.folder/'staging';self.engine.p=self.p
        self.engine.m=self.model;self.engine.solved={};self.engine.ritz_solves=0
        self.engine.cases=[];self.engine.catalogs={}
        for attr,default in [('entries',{}),('events',[]),('failures',[]),('independent',{})]:
            path=self.folder/(attr+'.json');setattr(self,attr,read(path) if path.exists() else default)
        self.new_this_run=0

    def catalogs(self):
        if self.engine.catalogs:return
        records={}
        for L in self.case(self.v,0).lengths(self.model):
            signature=arm_signature(self.model,L)
            key=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()[:20]
            path=self.parent.folder/'catalogs'/(key+'.json')
            if path.exists():
                record=read(path);require(record['signature']==signature,'Pole physical identity')
            else:
                start=time.perf_counter()
                cat=tim_catalog(self.model.carrier,L,self.model.omega(self.p['Lambda_upper']),
                    {'policy':{'scan_intervals':self.p['pole_scan_intervals'],'root_xtol':self.p['root_xtol'],'root_rtol':self.p['root_rtol']}})
                record=dict(signature=signature,catalog=cat,seconds=time.perf_counter()-start,status='PASS')
                atomic(path,record)
            require(record['status']==record['catalog']['status']=='PASS','Pole completeness')
            self.engine.catalogs[L]=record['catalog'];records[str(L)]=dict(signature=signature,file=path.relative_to(ROOT).as_posix(),sha256=sha(path))
        atomic(self.folder/'catalog_identity.json',records)

    def point(self,c,beta,reason,previous=None):
        if a.pilot.key(self.case(c,beta)) not in self.reference:self.catalogs()
        r=super().point(c,beta,reason,previous)
        e=self.entries[self.engine.name(self.case(c,beta))]
        if e['origin']=='new_HS_A':e['origin']='new_HS_B';self.flush()
        return r

    def stop(self,stage,error,status='UNRESOLVED'):
        self.failures.append(dict(id=self.v['id'],stage=stage,error=str(error),status=status,series_stopped=True));self.flush()

    def investigate(self):
        c=self.v;pairs=self.series(c)
        b=[e['beta_deg'] for e,r in pairs];w=[r['omega'][0] for e,r in pairs]
        internal=[i for i in range(1,len(w)-1) if w[i]>w[i-1] and w[i]>w[i+1]]
        if internal:
            self.extrema(c)
            self.stop('angular extrema','Interior maximum detected and diagnosed; further series stopped','STOP_INTERIOR_MAXIMUM')
            return
        for target in [w[-1],max(w)]:
            for q in self.cfg['analysis']['q']:
                if len(a.crossing_brackets(b,w,q*target))>1:
                    atomic(self.folder/'multiple_crossings.json',dict(q=q,target=target,brackets=a.crossing_brackets(b,w,q*target)))
                    self.stop('thresholds','Multiple bracketed crossings','STOP_MULTIPLE_CROSSINGS');return
        self.neighbour_audit(c)
        if c['group']=='B0':
            old=read(self.parent.baseline/'thresholds'/(a.construction_id(c['k_joint'],c['kappa_horizontal'])+'.json'))
            rows=[{**x,'id':c['id'],'source':'accepted_HS_A'} for x in old]
        else:rows=self.thresholds(c)
        require(all(r['status']!='UNRESOLVED_CROSSING' for r in rows),'Threshold signs unresolved')
        require(all(len(r['crossings'])<=1 for r in rows),'Multiple refined crossings')
        self.neighbour_audit(c)
        atomic(self.folder/'thresholds.json',rows)


class Study:
    def __init__(self):
        self.cfg=read(CONFIG);self.variants=variants(self.cfg)
        self.identity=dict(config=self.cfg,config_sha256=sha(CONFIG),initial_sha256=sha(ROOT/self.cfg['initial_workspace']),
            versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
        self.run_id=hashlib.sha256(json.dumps(self.identity,sort_keys=True).encode()).hexdigest()[:16]
        self.folder=OUTPUT/self.run_id;self.folder.mkdir(parents=True,exist_ok=True)
        source=self.folder/'workflow_versions'/(sha(__file__)+'.txt');source.parent.mkdir(exist_ok=True)
        if not source.exists():source.write_bytes(Path(__file__).read_bytes())
        self.baseline=ROOT/self.cfg['baseline_archive']
        self.baseline_entries=read(self.baseline/'entries.json')
        self.children={v['id']:GeometryStudy(self,v) for v in self.variants}

    def prepare(self):
        for path,digest in self.cfg['sources'].items():require(sha(ROOT/path)==digest,'Frozen source '+path)
        require(sha(self.baseline/'manifest.json')==self.cfg['baseline_manifest_sha256'],'HS-A manifest')
        m=read(self.baseline/'manifest.json');require(m['status']=='PASS','HS-A status')
        require(m['identity']['config']['material']==self.cfg['material'],'Baseline material')
        require(m['identity']['config']['mu']==self.cfg['baseline_mu'],'Baseline lengths')
        require(m['identity']['config']['policy']==self.cfg['policy'],'Physical tolerances')
        for name,digest in m['artifacts'].items():require(sha(self.baseline/name)==digest,'HS-A artifact '+name)
        for e in self.baseline_entries.values():
            require(sha(ROOT/e['path'])==e['sha256'] and sha(ROOT/e['npz_path'])==e['npz_sha256'],'HS-A accepted point')
        grids={};rows=[]
        for v in self.variants:
            child=self.children[v['id']];_,row=geometry(v,self.cfg);rows.append(row)
            raw=[(n['beta_deg'],o) for n in a.base_nodes(v['kappa_horizontal'],child.cfg) for o in n['origins']]
            raw += [(key[0],'accepted_HS_A') for key in child.reference]
            grids[v['id']]=a.merge_nodes(raw,self.cfg['grid']['duplicate_absolute_deg'])
        if (self.folder/'grid.json').exists():require(read(self.folder/'grid.json')==grids,'Frozen geometry grids')
        atomic(self.folder/'grid.json',grids);atomic(self.folder/'contract.json',self.identity)
        csvwrite(self.folder/'geometry.csv',rows)
        atomic(self.folder/'audit.json',dict(status='PASS',baseline_artifacts=len(m['artifacts']),baseline_points=len(self.baseline_entries),
            eligible_baseline_points=sum(len(c.reference) for c in self.children.values()),geometry_limit_pass=all(r['status']=='PASS' for r in rows),
            max_limiter=max(r['limiter'] for r in rows),model='unchanged HS-P1 full physics and HS-A point transactions',
            aspect_coefficients=dict(unchanged=['EA/EI=1/(epsilon*l)^2','KGA/EI','rhoA/EI','rhoI/EI'],
                changed='J_RL/(rhoA*l^2)=nu^2*epsilon^2*(1+aspect^2)')))
        for r in rows:
            if r['status']!='PASS':self.children[r['id']].stop('geometry limiter',r['limiter'],'GEOMETRY_STOP')
        print(self.folder,'frozen',sum(map(len,grids.values())),'nodes',flush=True)
        return grids

    def preflight(self):
        grid=read(self.folder/'grid.json');path=self.folder/'preflight.json'
        if path.exists():return read(path)
        rows=[];start=time.perf_counter()
        for name,node in self.cfg['preflight']:
            s=self.children[name]
            if s.failures:continue
            beta=next(n['beta_deg'] for n in grid[name] if node in n['origins']) if isinstance(node,str) else node
            t=time.perf_counter()
            try:
                r=s.point(s.v,beta,'preflight');check=s.independent_check(s.v,beta,'geometry preflight',orders=(16,20))
                rows.append(dict(id=name,beta=beta,omega=r['omega'],independent=check,seconds=time.perf_counter()-t,status='PASS'))
            except Exception as exc:
                s.stop('preflight',repr(exc));rows.append(dict(id=name,beta=beta,error=repr(exc),status='UNRESOLVED'))
            print('preflight',name,rows[-1]['status'],flush=True)
        result=dict(status='PASS' if all(r['status']=='PASS' for r in rows) else 'PARTIAL',checks=rows,seconds=time.perf_counter()-start)
        atomic(path,result);return result

    def run_map(self):
        grid=self.prepare();self.preflight()
        for s in self.children.values():
            if s.failures:continue
            try:
                previous=None
                for n in grid[s.v['id']]:previous=s.point(s.v,n['beta_deg'],'frozen_grid',previous)
                print(s.v['id'],'grid',len(grid[s.v['id']]),'PASS',flush=True)
                if not (s.folder/'thresholds.json').exists():s.investigate()
                s.export()
                print(s.v['id'],'analysis','STOP' if s.failures else 'PASS',flush=True)
            except Exception as exc:s.stop('map or analysis',repr(exc));print('SERIES STOP',s.v['id'],repr(exc),flush=True)
        self.export()

    def export(self):
        spectra=[];summary=[];thresholds=[];failures=[]
        for s in self.children.values():
            rr,ss=s.export();meta={k:s.v[k] for k in ['regime','group','epsilon','mu','aspect']}
            spectra += [{**r,**meta} for r in rr]
            summary += [{**r,**meta,'status':'STOPPED' if s.failures else ('PASS' if len(rr)>=len(read(self.folder/'grid.json')[s.v['id']]) else 'INCOMPLETE')} for r in ss]
            path=s.folder/'thresholds.json'
            if path.exists():thresholds += [{**r,**meta} for r in read(path)]
            failures += s.failures
        csvwrite(self.folder/'spectra.csv',spectra);csvwrite(self.folder/'summary.csv',summary)
        atomic(self.folder/'thresholds.json',thresholds);atomic(self.folder/'failures.json',failures)
        csvwrite(self.folder/'thresholds.csv',[{**{k:v for k,v in r.items() if k not in ['crossings','grid_intervals']},
            'crossings_json':json.dumps(r['crossings']),'grid_intervals_json':json.dumps(r['grid_intervals'])} for r in thresholds])
        return spectra,summary,thresholds


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','preflight','map','compare','contrasts','finalize','verify','plot_only'])
    args=parser.parse_args();s=Study()
    if args.action=='prepare':s.prepare()
    elif args.action=='preflight':s.prepare();print(s.preflight()['status'])
    elif args.action=='map':s.run_map()
    elif args.action=='compare':compare(s)
    elif args.action=='contrasts':contrasts(s)
    elif args.action=='finalize':finalize(s)
    elif args.action=='verify':
        m=read(s.folder/'manifest.json')
        for name,digest in m['artifacts'].items():require(sha(s.folder/name)==digest,'Artifact '+name)
        print(m['status'],s.folder)
    else:
        from scripts.lib.two_springs_horizontal_geometry_figures import render
        render(s.folder)


def checked_leading(model,tau,policy):
    result=model.solve(tau);d=result['diagnostics']
    for key,tol in [('natural_boundary','boundary_tol'),('strong_ode','strong_tol'),('essential','essential_tol'),
                    ('matrix_residual','matrix_tol'),('energy','energy_tol'),('mass_gram','mass_tol'),('symmetry','matrix_tol')]:
        require(d[key]<=policy[tol],f'Leading {model.k_joint},{model.mu},{tau}: {key}={d[key]}')
    return result


def comparison_case(folder,k,mu,cfg):
    """One dimensionless pair; epsilon/aspect do not enter the two limit models."""
    from scripts.lib.two_springs_horizontal_joint_stiffness import endpoint_count
    name=f'j{k:g}_mu{mu:g}';path=folder/(name+'.json')
    if path.exists():return read(path)
    pol=cfg['leading_policy'];rp=cfg['leading_root'];models={p:Leading(k,mu,p) for p in rp['orders']}
    lows={};outs={};records=[];targets=[]
    for p,model in models.items():
        for tau in [0.,a.h.RIGID]:
            sol=checked_leading(model,tau,pol);lows[p,tau]=sol
            checks=independent_roots(sol['zeta'],tau,k,mu) if p==24 else []
            require(all(abs(r['relative_frequency_error'])<=pol['independent_frequency_tol'] for r in checks),'Leading ODE limit')
            records.append(dict(order=p,tau=tau,zeta=sol['zeta'].tolist(),diagnostics=sol['diagnostics'],independent=checks))
        outer=RightAngle(k,mu,p);sol=outer.solve();d=sol['diagnostics'];outs[p]=sol
        for key,tol in [('natural','boundary_tol'),('strong','strong_tol'),('essential','essential_tol'),
                        ('matrix_residual','matrix_tol'),('energy','energy_tol'),('mass_gram','mass_tol'),('symmetry','matrix_tol')]:
            require(d[key]<=pol[tol],f'Outer {k},{mu}: {key}')
        checks=independent_checks(sol['zeta'],k,mu,pol) if p==24 else None
        records.append(dict(order=p,tau='OUTER',zeta=sol['zeta'].tolist(),diagnostics=d,independent=checks))
    convergence=[]
    for label,aa,bb in [(t,lows[20,t],lows[24,t]) for t in [0.,a.h.RIGID]]+[('OUTER',outs[20],outs[24])]:
        err=float(max(abs(np.sqrt(aa['zeta']/bb['zeta'])-1)))
        require(err<=pol['frequency_convergence_tol'],'Limit order convergence')
        convergence.append(dict(tau=label,relative_frequency_error=err))
    for kind in ['leading','combined']:
        for q in cfg['analysis']['q']:
            roots=[];checks=[]
            for p,model in models.items():
                zmax=lows[p,a.h.RIGID]['zeta'][0]
                target=q*q*(zmax if kind=='leading' else outs[p]['zeta'][0])
                if target>=zmax:
                    roots.append(None);continue
                require(target>lows[p,0.]['zeta'][0],'Target already high at zero')
                evaluations=[]
                def fun(tau):
                    sol=checked_leading(model,tau,pol)
                    evaluations.append(dict(tau=float(tau),zeta=float(sol['zeta'][0]),diagnostics=sol['diagnostics']))
                    return sol['zeta'][0]/target-1
                hi=1.
                while fun(hi)<0:
                    hi*=2
                    require(hi<=rp['maximum_tau'],'Leading target not bracketed')
                root=float(brentq(fun,0.,hi,xtol=rp['xtol'],rtol=rp['rtol']))
                sol=checked_leading(model,root,pol);res=float(sol['zeta'][0]/target-1)
                require(abs(res)<=rp['relative_equation_tol'],'Leading root residual')
                independent=independent_roots(sol['zeta'],root,k,mu)
                require(max(abs(r['relative_frequency_error']) for r in independent)<=pol['independent_frequency_tol'],'Leading root ODE')
                count=[]
                for i,z in enumerate(sol['zeta']):
                    low=endpoint_count(z*(1-pol['count_bracket_relative']),k,root,mu,pol['inertia_relative_zero_tol'],pol['pole_relative_guard'])
                    high=endpoint_count(z*(1+pol['count_bracket_relative']),k,root,mu,pol['inertia_relative_zero_tol'],pol['pole_relative_guard'])
                    require(low['count']==i and high['count']==i+1,'Leading first-three count')
                    count.append(dict(lower=low,upper=high))
                checks.append(dict(order=p,tau=root,target=float(target),residual=res,diagnostics=sol['diagnostics'],
                    independent=independent,count=count,evaluations=evaluations));roots.append(root)
            if roots==[None,None]:status='TARGET_AT_OR_ABOVE_LEADING_SATURATION';error=None
            else:
                require(all(t is not None for t in roots),'Order-dependent target existence')
                error=abs(math.sqrt(roots[0]/roots[1])-1)
                require(error<=pol['frequency_convergence_tol'],'Threshold order convergence');status='PASS'
            targets.append(dict(kind=kind,q=q,status=status,tau=roots[-1],order_angle_relative_error=error,checks=checks))
    result=dict(k_joint=k,mu=mu,zeta_zero=float(lows[24,0.]['zeta'][0]),zeta_infty=float(lows[24,a.h.RIGID]['zeta'][0]),
        zeta_outer=float(outs[24]['zeta'][0]),records=records,convergence=convergence,targets=targets,status='PASS')
    atomic(path,result);return result


def compare(s):
    folder=s.folder/'leading_comparison';folder.mkdir(exist_ok=True)
    source=ROOT/'scripts/lib/two_springs_horizontal_joint_stiffness.py'
    contract=dict(policy=s.cfg['leading_policy'],root=s.cfg['leading_root'],additional_counter_source=source.relative_to(ROOT).as_posix(),
        additional_counter_sha256=sha(source),no_new_asymptotics=True,
        qualification='numerical error at each requested point; no universal xi validity boundary',new_full_solves=0)
    if (folder/'contract.json').exists():require(read(folder/'contract.json')==contract,'Comparison contract changed')
    atomic(folder/'contract.json',contract)
    _,summary,thresholds=s.export();models={};errors=[];rows=[]
    for v in s.variants:
        if s.children[v['id']].failures:continue
        key=v['k_joint'],v['mu']
        if key not in models:
            try:models[key]=comparison_case(folder,*key,s.cfg)
            except Exception as exc:
                models[key]=None;errors.append(dict(k_joint=key[0],mu=key[1],status='UNRESOLVED',error=repr(exc)))
                atomic(folder/'failures.json',errors)
                print('LEADING STOP',key,repr(exc),flush=True)
        model=models[key]
        if model is None:continue
        child=s.children[v['id']]
        sm=next(r for r in summary if r['id']==v['id'])
        scale=math.sqrt(2*v['epsilon']**2+(0 if v['chi_s']==a.h.RIGID else 1/v['kappa_horizontal']))
        for t in [t for t in thresholds if t['id']==v['id']]:
            full=t['terminal_beta_root_deg'];xi=None if full is None else math.radians(full)/v['epsilon']
            values={}
            for kind in ['leading','combined']:
                rr=next(x for x in model['targets'] if x['kind']==kind and x['q']==t['q'])
                beta=math.degrees(scale*math.sqrt(rr['tau'])) if rr['tau'] is not None else None
                values.update({kind+'_tau':rr['tau'],kind+'_beta_deg':beta,kind+'_status':rr['status'],
                    kind+'_in_checked_angle_domain':None if beta is None else 0<=beta<=90,
                    kind+'_error_percent':None if beta is None or full is None else 100*(beta/full-1),
                    kind+'_prediction_xi':None if beta is None else math.radians(beta)/v['epsilon']})
            # Pointwise check of the existing inner curve at the FULL root, not calibration.
            spectral=None
            if full is not None:
                tau=(math.radians(full)/scale)**2
                inner=checked_leading(Leading(v['k_joint'],v['mu'],24),tau,s.cfg['leading_policy'])
                fullomega=t['q']*t['omega_target']
                spectral=100*(v['epsilon']*math.sqrt(inner['zeta'][0])/fullomega-1)
            rows.append(dict(id=v['id'],regime=v['regime'],group=v['group'],epsilon=v['epsilon'],mu=v['mu'],aspect=v['aspect'],
                q=t['q'],normalization=t['normalization'],full_beta_deg=full,full_xi=xi,**values,
                inner_frequency_error_at_full_root_percent=spectral,outer_frequency_error_percent=100*(v['epsilon']*math.sqrt(model['zeta_outer'])/sm['omega90']-1),
                qualification='FORMAL_EXTENSION_POINTWISE_CHECKED_NO_UNIFORM_ERROR_BOUND',
                aspect_correction_in_limits=False))
        print('comparison',v['id'],'PASS',flush=True)
    atomic(folder/'failures.json',errors);csvwrite(s.folder/'predictions.csv',rows)


def contrasts(s):
    rows,summary,thresholds=s.export();out=[];diagnostics=[]
    for v in s.variants:
        if v['group']=='B0':continue
        base={r['beta_deg']:r for r in rows if r['id']=='B0_'+v['regime']}
        for r in [r for r in rows if r['id']==v['id'] and r['beta_deg'] in base]:
            ref=base[r['beta_deg']]
            out.append(dict(id=v['id'],regime=v['regime'],group=v['group'],beta_deg=r['beta_deg'],
                relative_omega1=r['omega1']/ref['omega1']-1,relative_Lambda1=r['Lambda1']/ref['Lambda1']-1,
                reference_id=ref['id']))
    # Observed aspect differences are exceptionally small: compare independent
    # Ritz orders at the largest positive shared-node difference for each regime.
    # These are the same accepted angles, no full spectral calls or extra grids.
    for regime in ['A1','A2','A3','A4']:
        selected=max((r for r in out if r['regime']==regime and r['group']=='B3'),key=lambda r:r['relative_omega1'])
        beta=selected['beta_deg'];checks=[]
        for name in ['B0_'+regime,'B3_'+regime+'_a2']:
            child=s.children[name]
            check=child.independent_check(child.v,beta,'resolve very small aspect effect',orders=(16,20));checks.append(check)
        differences=[checks[1]['checks'][i]['omega'][0]/checks[0]['checks'][i]['omega'][0]-1 for i in [0,1]]
        child=s.children['B0_'+regime];entry=next(e for e in child.entries.values() if e['beta_deg']==beta)
        with np.load(ROOT/entry['npz_path']) as data:
            fraction=child.model.J*sum(float(data[f'mode1_arm{i}_weights']@data[f'mode1_arm{i}_dy'][:,0]**2) for i in [1,2])
        # Exact Rayleigh upper bound after dividing common aspect factor:
        # K2=K4, M2=M4-(12/17)M_RL4; no fitted or new asymptotic coefficient.
        mass_loss=12/17*fraction
        bound=mass_loss/(math.sqrt(1-mass_loss)*(1+math.sqrt(1-mass_loss)))
        diagnostics.append(dict(regime=regime,beta_deg=beta,full_relative_difference=selected['relative_omega1'],
            ritz_relative_differences=differences,ritz_order_spread=abs(differences[1]-differences[0]),
            RL_mass_fraction=fraction,rayleigh_relative_frequency_upper_bound=bound,
            interpretation='Exact normalized mass ordering; tiny finite differences reported with independent order spread',status='PASS'))
    csvwrite(s.folder/'geometry_contrasts.csv',out);atomic(s.folder/'aspect_diagnostics.json',diagnostics)
    print('Aspect diagnostics',[(r['regime'],r['full_relative_difference'],r['ritz_relative_differences']) for r in diagnostics])


def finalize(s):
    rows,summary,thresholds=s.export();grid=read(s.folder/'grid.json')
    limits=[];neighbours=[];entries=[];independent=[]
    for child in s.children.values():
        entries+=list(child.entries.values());independent+=list(child.independent.values())
        if child.failures:continue
        done={e['beta_deg'] for e in child.entries.values()}
        require(all(n['beta_deg'] in done for n in grid[child.v['id']]),'Frozen grid incomplete')
        for e,r in child.series(child.v):
            require(sha(ROOT/e['path'])==e['sha256'] and sha(ROOT/e['npz_path'])==e['npz_sha256'],'Accepted point hash')
            a.pilot.validate_record(r,child.case(child.v,e['beta_deg']),child.p)
            require(r['count']['zero_gate']['expected']==0,'Unexpected zero mode')
            require(np.max(abs(np.array(r['omega'])/(child.v['epsilon']*np.array(r['Lambda'])**2)-1))<1e-13,'Frequency normalization')
        for beta in [0.,90.]:
            case=child.case(child.v,beta);ext=a.h.external_matrix(child.model,case,1)
            trace,_=a.h.trace_map(child.model,case);co,si=a.h.cs(case)
            if case.kappa_horizontal==a.h.RIGID:
                err=float(max(abs(co*trace[6]+si*trace[7])));require(err<=child.p['boundary_tol'] and trace.shape[1]==6,'Exact horizontal constraint')
            else:
                err=0.;require(np.linalg.matrix_rank(ext[:2,:2])==1 and trace.shape[1]==7,'Rank-one finite support')
            require(np.array_equal(ext[2],np.eye(6)[5]),'External moment-free condition')
            limits.append(dict(id=child.v['id'],beta=beta,trace_dimension=trace.shape[1],constraint_error=err,status='PASS'))
        for (ea,aa),(eb,bb) in zip(child.series(child.v)[:-1],child.series(child.v)[1:]):
            jump=float(max(abs(np.array(aa['omega'][:2])-bb['omega'][:2])/np.maximum(aa['omega'][:2],bb['omega'][:2])))
            neighbours.append(dict(id=child.v['id'],left=ea['beta_deg'],right=eb['beta_deg'],jump=jump,triggered=jump>s.cfg['frequency_map']['neighbour_jump_relative']))
    require(all(c['status']=='PASS' for c in independent),'Independent verification unresolved')
    for t in thresholds:
        require(t['status']!='UNRESOLVED_CROSSING','Unresolved root')
        for r in t['crossings']:
            require(r['status']=='PASS' and r['signs'][0]*r['signs'][1]<0 and abs(r['equation_residual'])<=s.cfg['analysis']['threshold_residual_tol'],'Threshold acceptance')
    initial=read(ROOT/s.cfg['initial_workspace']);preserved=0
    for path,digest in initial['protected_archives'].items():
        require(sha(ROOT/path)==digest,'Historical archive modified '+path);preserved+=1
    for path,digest in s.cfg['sources'].items():require(sha(ROOT/path)==digest,'Frozen source changed '+path)
    atomic(s.folder/'limit_checks.json',limits);atomic(s.folder/'neighbour_audit.json',neighbours)
    failures=read(s.folder/'failures.json');leadingfail=read(s.folder/'leading_comparison/failures.json')
    maxres={key:max(r[key] for r in rows) for key in ['boundary_max','ode_max','energy_max','mass_gram_error','singular_max']}
    result=dict(status='PASS' if not failures and not leadingfail and len(thresholds)==136 else 'PARTIAL',identity=s.identity,
        HEAD=initial['HEAD'],initial_dirty_files=len(initial['files']),protected_archive_files=preserved,
        variants=len(s.variants),completed_variants=sum(r['status']=='PASS' for r in summary),
        frozen_grid_points=sum(map(len,grid.values())),accepted_points=len(entries),
        new_full_points=sum(e['origin']=='new_HS_B' for e in entries),reused_full_points=sum(e['origin']=='accepted_archive' for e in entries),
        independent_ritz_solves=sum(len(c['checks']) for c in independent),
        new_full_seconds=sum(e['seconds'] for e in entries if e['origin']=='new_HS_B'),
        pole_catalogs=len(list((s.folder/'catalogs').glob('*.json'))),neighbour_pairs=len(neighbours),neighbour_triggers=sum(r['triggered'] for r in neighbours),
        minimum_g12=min(r['g12'] for r in rows),max_residuals=maxres,
        max_independent_frequency_error=max(c['frequency_error'] for x in independent for c in x['checks']),
        threshold_rows=len(thresholds),physical_thresholds=len(thresholds)//2,failures=failures,leading_failures=leadingfail)
    result['artifacts']={p.relative_to(s.folder).as_posix():sha(p) for p in s.folder.rglob('*') if p.is_file()
        and 'staging' not in p.relative_to(s.folder).parts and 'figures' not in p.relative_to(s.folder).parts
        and p.name not in ['manifest.json','README.md','figures_manifest.json','completion_checks.json','documentation_manifest.json']}
    atomic(s.folder/'manifest.json',result);print(result['status'],s.folder,len(entries),'points',flush=True)


if __name__=='__main__':main()
