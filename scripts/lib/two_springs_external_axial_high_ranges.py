"""Axial-support preset: immutable archived spectra plus bounded angle extension.

Reuses the canonical physical solver, count/field gates and high_ranges.
No physical assembly or spectral infrastructure is introduced here.
"""
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
import scipy
from scripts.lib import two_springs_external_axial_spring as previous
from scripts.lib import two_springs_high_frequency_angle_ranges as ranges
from scripts.lib import isotropic_two_springs_external_transverse_spring as transverse

ROOT=previous.ROOT
BASE=previous.BASE.with_name('external_axial_high_frequency_ranges')
CONFIG=ROOT/'data/input/isotropic_two_springs_external_axial_high_frequency_ranges.json'
physics=previous.physics
api=previous.api
require=previous.require
read_csv=previous.read_csv
chi=previous.chi
FIGURES=ranges.FIGURES


def archive_case(c):
    """Only an exact transverse support is physically identical to axial RIGID."""
    require(c.get('k_w','RIGID')=='RIGID', 'Transverse-compliant spectrum is not axial-support data')
    fields={k:v for k,v in c.items() if k!='variant'}
    cls=physics.Case if 'chi_u' in c else transverse.Case if 'k_w' in c else physics.base.Case
    return cls(**fields)


class Experiment(previous.Experiment):
    def __init__(self):
        self.a=api();self.cfg=json.loads(CONFIG.read_text());self.inputs={};self.sources={};self.manifests={}
        self.started=time.perf_counter();self.ritz_new=0;self.leading_new=0
        for name,digest in self.cfg['frozen_files'].items():self.check(ROOT/name,digest)
        for name,r in self.cfg['sources'].items():
            folder=ROOT/r['folder'];self.check(folder/'manifest.json',r['manifest_sha256'])
            ranges.saved.diagnostic.verify_source(folder)
            self.sources[name]=folder;self.manifests[name]=json.loads((folder/'manifest.json').read_text())
            require(self.manifests[name]['status']=='PASS','Source status')
            c=json.loads(self.source(name,'contract.json').read_text())['config']
            require(c['material']==self.cfg['material'] and c['mu']==self.cfg['mu'] and c['variant']=='RL_FULL','Source geometry/material')
            require(c['map']['spectrum_semantics']=='sorted_positions' and c['map']['K_guard']==3,'Source semantics')
            if name=='axial':require(c['cases']==self.cfg['cases'] and c['chi_values']==self.cfg['chi_values'],'Axial supports')
        identity=dict(config=self.cfg,physical_sha256=self.a.sha(physics.__file__),versions=dict(numpy=np.__version__,scipy=scipy.__version__))
        self.key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16]
        self.folder=BASE/self.key;self.folder.mkdir(parents=True,exist_ok=True)
        target=self.folder/'contract.json'
        if target.exists():require(json.loads(target.read_text())==identity,'Immutable contract')
        else:self.a.write(target,identity)
        self.run=self.a.Run(self.cfg,self.folder,engine=physics)
        # Same two individual rods and same fixed-displacement dynamic stiffness.
        # Only the assembled trace space depends on chi; the source uses it already.
        self.run.catalogs={float(k):v for k,v in json.loads(self.source('axial','pole_catalogs.json').read_text()).items()}
        self.a.write(self.folder/'pole_catalogs.json',self.run.catalogs)
        self.archive={};self.overlap=[]
        for source in ('axial','ranges'):
            for row in read_csv(self.source(source,'spectrum.csv')):
                if source=='ranges' and (row['pair'] not in ('I_J01','I_J10') or row['k_w']!='RIGID'):continue
                self.ingest(ROOT/row['point_file'],row['point_sha256'],source,row)
        # Existing exact-support regressions are checks of coincident points.
        for name,digest in self.manifests['axial']['artifacts'].items():
            if name.startswith('points/') and name.endswith('.json'):
                self.ingest(self.sources['axial']/name,digest,'axial')
        inventory=[]
        for p in self.cfg['cases']:
            for k in self.cfg['chi_values']:
                have=sorted(b for pp,kk,b in self.archive if (pp,kk)==(p,k))
                require(all(b in have for b in self.cfg['map']['parameter_grid'] if b<=20),'Missing archived main grid')
                if k=='RIGID':require(have[-1]==90,'Rigid reference not complete')
                inventory.append(dict(pair=p,chi_u=k,available_count=len(have),available_beta_deg=have,
                                      missing_initial_beta_deg=[b for b in self.cfg['map']['parameter_grid'] if b not in have]))
        self.a.write(self.folder/'archive_inventory.json',inventory);self.csv('archive_overlaps',self.overlap)
        self.plan=json.loads((self.folder/'refinement_plan.json').read_text()) if (self.folder/'refinement_plan.json').exists() else []
        self.allowed={}
        for r in self.plan:self.allowed.setdefault((r['pair'],r['chi_u']),set()).update(r['new_beta_deg'])

    def ingest(self,path,digest,source,row=None):
        self.check(path,digest);rec=json.loads(path.read_text());c=rec['case'];case=archive_case(c)
        pair=next((p for p,v in self.cfg['cases'].items() if v==[c['k_joint'],c['k_end']]),None)
        require(pair is not None and c['mu']==self.cfg['mu'],'Unexpected source construction')
        self.check(path.with_suffix('.npz'),rec['shape_sha256']);w=self.validate(rec,case)
        if row is not None:require(max(abs(w/np.array([float(row[f'omega{i}']) for i in (1,2,3)])-1))<1e-13,'CSV/root mismatch')
        key=pair,c.get('chi_u','RIGID'),float(c['beta'])
        if key in self.archive:
            old=self.archive[key];err=float(max(abs(w/np.array([r['omega'] for r in old[0]['roots']])-1)))
            require(err<=self.cfg['policy']['same_physics_frequency_tol'],'Archive overlap differs')
            self.overlap.append(dict(pair=pair,chi_u=key[1],beta_deg=key[2],source=source,relative_frequency_error=err,status='PASS'))
        else:self.archive[key]=(rec,path,source)

    def point(self,pair,k,beta,previous_w=None):
        if (pair,k,beta) in self.archive:return self.archive[pair,k,beta]
        require(k!='RIGID','Do not recalculate the accepted rigid map')
        require(beta in self.cfg['map']['parameter_grid'] or beta in self.allowed.get((pair,k),set()),'Unplanned angle')
        case=self.case(pair,k,beta);rec=self.run.solve(case,neighbour=previous_w);self.validate(rec,case)
        return rec,previous.lengths.point_path(self.folder,case),'new_study'

    def case(self,pair,k,beta):
        # Canonical JSON spelling also matters to the existing cache hash.
        return physics.Case(float(beta),float(self.cfg['mu']),*self.cfg['cases'][pair],chi(k))

    def preflight(self):
        path=self.folder/'preflight.json'
        if path.exists():require(json.loads(path.read_text())['status']=='PASS','Saved preflight failure');return
        rows=[];ritz=[];identity=[]
        for p in self.cfg['cases']:
            for b in self.cfg['preflight_angles']:
                c=self.case(p,'RIGID',b);base=physics.base.Case(b,c.mu,c.k_joint,c.k_end)
                require(np.array_equal(physics.external_matrix(self.run.model,c,1),physics.base.external_matrix(self.run.model,base,1)),'Rigid boundary identity')
                require(all(np.array_equal(a,b) for a,b in zip(physics.trace_map(self.run.model,c),physics.base.trace_map(self.run.model,base))),'Rigid trace identity')
                identity.append(dict(pair=p,beta_deg=b,status='PASS'))
            for k in ['RIGID']+self.cfg['preflight_chi']:
                for b in self.cfg['preflight_angles']:
                    rec,_,source=self.point(p,k,b);roots=rec['roots'];w=np.array([r['omega'] for r in roots])
                    rows.append(dict(pair=p,chi_u=k,beta_deg=b,source=source,Lambda1=roots[0]['Lambda'],Lambda3=roots[2]['Lambda'],lower_count=rec['count']['samples'][0]['count'],upper_count=rec['count']['samples'][1]['count'],status='PASS'))
                    if [k,b] in self.cfg['ritz_cases']:
                        for order in self.cfg['policy']['ritz_orders']:
                            wr,_,d=physics.ritz(self.run.model,self.case(p,k,b),order);self.ritz_new+=1
                            err=float(max(abs(wr/w-1)));ritz.append(dict(pair=p,chi_u=k,beta_deg=b,order=order,relative_frequency_error=err,**d))
                            self.csv('preflight_ritz',ritz)
                            require(err<=self.cfg['policy']['ritz_frequency_tol'] and d['positive_K_M'] and d['mass_gram_error']<=self.cfg['policy']['mass_gram_tol'] and d['energy_error']<=self.cfg['policy']['energy_tol'] and d['constraint_residual']<=self.cfg['policy']['boundary_tol'],'Independent Ritz gate')
                    self.csv('preflight_points',rows)
        self.csv('rigid_identity',identity)
        self.a.write(path,dict(status='PASS',points=len(rows),new_points=self.run.new_points,ritz_solves=len(ritz),max_ritz_error=max(r['relative_frequency_error'] for r in ritz),zero_mechanism=False,zero_evidence='Positive Ritz K/M at chi=0,beta=90; zero lower spectral counts; unchanged zero-energy argument'))

    def sweep(self):
        rows=[];quality=[];neighbours=[];ap=self.folder/'triggered_counts.json'
        audits=json.loads(ap.read_text()) if ap.exists() else {}
        for p in self.cfg['cases']:
            for k in self.cfg['chi_values']:
                grid={b for pp,kk,b in self.archive if (pp,kk)==(p,k)}|self.allowed.get((p,k),set())
                if k!='RIGID':grid|=set(self.cfg['map']['parameter_grid'])
                prev=None
                for b in sorted(grid):
                    rec,file,source=self.point(p,k,b,prev);w=np.array([r['omega'] for r in rec['roots']]);jump=float(max(abs(w/prev-1))) if prev is not None else 0.;gap=w[1]/w[0]-1
                    trigger=jump>self.cfg['policy']['neighbour_relative_jump'] or gap<self.cfg['refinement']['small_gap12']
                    if trigger:
                        key=f'{p}|{k}|{b}'
                        if key not in audits:
                            c=self.case(p,k,b);cats=self.run.cats(self.run.model,c)
                            tests=[physics.count_at(self.run.model,c,x,cats,self.cfg['policy']) for x in (.9*w[0],(w[0]+w[1])/2,(w[1]+w[2])/2)]
                            require([r[1] for r in tests]==[0,1,2],'Triggered count');audits[key]=[r[2] for r in tests];self.a.write(ap,audits)
                    neighbours.append(dict(pair=p,chi_u=k,beta_deg=b,jump=jump,gap12=gap,triggered=trigger,status='PASS'));prev=w
                    rows.append(dict(pair=p,k_joint=self.cfg['cases'][p][0],k_end=self.cfg['cases'][p][1],chi_u=k,beta_deg=b,**{f'{n}{i}':r[n] for i,r in enumerate(rec['roots'],1) for n in ('omega','Lambda')},gap12=gap,origin=source,point_file=file.relative_to(ROOT).as_posix(),point_sha256=self.a.sha(file),shape_sha256=rec['shape_sha256']))
                    for i,r in enumerate(rec['roots'],1):
                        d=r['diagnostics'];quality.append(dict(pair=p,chi_u=k,beta_deg=b,position=i,origin=source,physical=float(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))),ode=d['ode_error'],energy=d['energy_error'],singular=d['singular_ratio'],mass_gram=rec['mass_gram_error'],condition=d['nonzero_condition'],count_lower=r['counts'][0],count_upper=r['counts'][1],low_count=rec['count']['samples'][0]['count'],upper_count=rec['count']['samples'][1]['count'],status='PASS'))
                for name,data in [('spectrum',rows),('quality',quality),('neighbours',neighbours)]:self.csv(name,data)
                print(f'{p}, chi_u={k}: {len(grid)} checked nodes',flush=True)
        return rows

    def analyze(self,rows,prefix=''):
        return analyze(self.folder,self.cfg,rows,prefix)

    def refine(self,rows,round_number):
        # Reuse the established bounded refinement logic on a temporary column
        # view; k_w here is only its legacy column key, never a physical parameter.
        from types import SimpleNamespace
        cfg={**self.cfg,'k_w_values':self.cfg['chi_values']}
        oldplan=[{**{n:v for n,v in r.items() if n!='chi_u'},'k_w':r['chi_u']} for r in self.plan]
        proxy=SimpleNamespace(cfg=cfg,revised=False,plan=oldplan,allowed=self.allowed,a=self.a,folder=self.folder)
        view=[{**r,'k_w':r['chi_u']} for r in rows if r['chi_u']!='RIGID']
        proxy.cfg={**cfg,'k_w_values':[k for k in cfg['k_w_values'] if k!='RIGID']}
        added=ranges.Experiment.refine(proxy,view,round_number)
        self.plan=[{**{n:v for n,v in r.items() if n!='k_w'},'chi_u':r['k_w']} for r in proxy.plan]
        self.a.write(self.folder/'refinement_plan.json',self.plan)
        return added

    def run_study(self):
        self.preflight();rows=self.sweep()
        if not (self.folder/'initial_summary.csv').exists():self.analyze(rows,'initial_')
        for n in range(1,self.cfg['refinement']['max_rounds']+1):
            if not self.refine(rows,n):break
            rows=self.sweep()
        else:require(not self.refine(rows,n+1),'Unresolved refinement')
        summary=self.analyze(rows)
        for r in summary:
            for d in self.cfg['levels_percent']:require(r[f'resolution_{d}'] is None or r[f'resolution_{d}']<=.25+1e-10,'Unresolved terminal boundary')
        self.comparisons(rows,summary)
        mechanics(self.folder,self.cfg,rows,summary)
        self.a.write(self.folder/'calculation_complete.json',dict(status='PASS',constructions=len(summary),points=len(rows),reused_points=sum(r['origin']!='new_study' for r in rows),new_points=sum(r['origin']=='new_study' for r in rows),refinement_rounds=max((r['round'] for r in self.plan),default=0)))

    def comparisons(self,rows,summary):
        transverse_rows=[r for r in read_csv(self.source('ranges','summary.csv')) if r['variant']=='I']
        self.csv('transverse_comparison',transverse_rows)
        reference=[];monotonic=[]
        for p,p_old in [('J0E4','I_J01'),('J4E0','I_J10')]:
            old=next(r for r in transverse_rows if r['pair']==p_old and r['k_w']=='RIGID')
            new=next(r for r in summary if r['pair']==p and r['chi_u']=='RIGID')
            for n in ('omega0','omega_max','beta_max','G_max_percent','beta_high_10','beta_high_5'):
                require(float(old[n])==new[n],'Rigid archival range changed')
            reference.append(dict(pair=p,reference_pair=p_old,omega_max=new['omega_max'],beta_high_10=new['beta_high_10'],beta_high_5=new['beta_high_5'],status='PASS'))
            by={(r['chi_u'],r['beta_deg']):r for r in rows if r['pair']==p}
            common=set.intersection(*[{b for kk,b in by if kk==k} for k in self.cfg['chi_values']])
            order=[0.,.01,.1,1.,'RIGID']
            for b in sorted(common):
                for lo,hi in zip(order,order[1:]):
                    change=[by[hi,b][f'omega{i}']/by[lo,b][f'omega{i}']-1 for i in (1,2,3)]
                    require(min(change)>=-self.cfg['policy']['same_physics_frequency_tol'],'Stiffness monotonicity violated')
                    monotonic.append(dict(pair=p,beta_deg=b,chi_lower=lo,chi_upper=hi,minimum_relative_increment=min(change),status='PASS'))
        self.csv('rigid_range_regression',reference);self.csv('stiffness_monotonicity',monotonic)

    def finish(self,stage,status='PASS'):
        super().finish(stage,status)
        p=self.folder/'executions.json';runs=json.loads(p.read_text());runs[-1]['workflow_sha256']=self.a.sha(__file__);self.a.write(p,runs)


def analyze(folder,cfg,rows,prefix=''):
    summary=[];intervals=[];deficits=[]
    for p in cfg['cases']:
        for k in cfg['chi_values']:
            rr=sorted([r for r in rows if r['pair']==p and chi(r['chi_u'])==k],key=lambda r:float(r['beta_deg']))
            b=np.array([float(r['beta_deg']) for r in rr]);w=np.array([float(r['omega1']) for r in rr]);lam=np.array([float(r['Lambda1']) for r in rr]);i=int(np.argmax(w))
            require(max(abs(w-.005*lam**2)/w)<2e-14,'Frequency normalization')
            result=ranges.high_ranges(b,w,10)
            row=dict(pair=p,k_joint=cfg['cases'][p][0],k_end=cfg['cases'][p][1],chi_u=k,nodes=len(rr),Lambda0=lam[0],Lambda_max=lam[i],omega0=w[0],omega_max=w[i],beta_max=b[i],G_max_percent=result['G_max_percent'],R90=result['R90'],min_gap12=min(float(r['gap12']) for r in rr),omega_min=float(min(w)),beta_min=float(b[np.argmin(w)]))
            for d in cfg['levels_percent']:
                result=ranges.high_ranges(b,w,d)
                row.update({f'beta_high_{d}':result['beta_high'],f'status_{d}':result['status'],f'resolution_{d}':result['resolution'],f'first_high_{d}':result['first_high'],f'terminal_nodes_{d}':result['terminal_nodes']})
                intervals.extend(dict(pair=p,chi_u=k,delta_percent=d,**v) for v in result['intervals'])
            if row['beta_high_10'] is not None and row['beta_high_5'] is not None:require(row['beta_high_10']<=row['beta_high_5'],'Nested thresholds')
            summary.append(row)
            deficits.extend(dict(pair=p,chi_u=k,beta_deg=float(x),omega1=float(y),Lambda1=float(z),R_percent=float(R),qualifies_10=bool(y>=.9*w[i]),qualifies_5=bool(y>=.95*w[i])) for x,y,z,R in zip(b,w,lam,result['R']))
    for name,data in [('summary',summary),('high_intervals',intervals),('deficit',deficits)]:api().csvwrite(folder/(prefix+name+'.csv'),data)
    return summary


def finalize(folder):
    folder=folder.resolve()
    a=api();cfg=json.loads((folder/'contract.json').read_text())['config'];done=json.loads((folder/'calculation_complete.json').read_text())
    require(done['status']=='PASS' and not (folder/'failure.json').exists(),'Unaccepted study')
    for name in FIGURES:
        for ext in ('.pdf','.png'):require((folder/(name+ext)).exists(),'Missing figure')
    for name,digest in json.loads((folder/'source_hashes.json').read_text()).items():require(a.sha(ROOT/name)==digest,'Changed source '+name)
    a.write(folder/'source_integrity.json',{n:ranges.saved.diagnostic.verify_source(ROOT/r['folder']) for n,r in cfg['sources'].items()})
    runs=json.loads((folder/'executions.json').read_text());q=read_csv(folder/'quality.csv');s=read_csv(folder/'summary.csv')
    newq=[r for r in q if r['origin']=='new_study'];ov=read_csv(folder/'archive_overlaps.csv')
    decision=dict(**done,new_full_prefixes=sum(r['new_full_prefixes'] for r in runs),new_full_roots=sum(r['new_full_roots'] for r in runs),new_pole_roots=sum(r['new_pole_roots'] for r in runs),independent_ritz_solves=sum(r['new_ritz_solves'] for r in runs),new_leading_solves=0,
                  quality_maxima={n:max(float(r[n]) for r in q) for n in ('physical','ode','energy','singular','mass_gram','condition')},new_quality_maxima={n:max(float(r[n]) for r in newq) for n in ('physical','ode','energy','singular','mass_gram','condition')},overlaps=len(ov),max_overlap_error=max(float(r['relative_frequency_error']) for r in ov),no_terminal_cases={str(d):sum(r[f'status_{d}']=='NO_TERMINAL_HIGH_RANGE' for r in s) for d in cfg['levels_percent']})
    # Retain the six initial chi=0 integer-key preflight records. They represent
    # repeated checks of the same physical points, not additional grid nodes.
    canonical={json.dumps({**json.loads((ROOT/r['point_file']).read_text())['case'],'chi_u':chi(r['chi_u'])},sort_keys=True):r for r in read_csv(folder/'spectrum.csv') if r['origin']=='new_study'}
    duplicates=[]
    for path in (folder/'points').glob('*.json'):
        rec=json.loads(path.read_text());c=rec['case'];key=json.dumps({**c,'chi_u':chi(c['chi_u'])},sort_keys=True)
        require(key in canonical,'Unaccounted calculated point');row=canonical[key]
        if path.resolve()!=(ROOT/row['point_file']).resolve():
            err=max(abs(r['omega']/float(row[f'omega{i}'])-1) for i,r in enumerate(rec['roots'],1))
            require(err<=cfg['policy']['same_physics_frequency_tol'],'Repeated preflight differs')
            duplicates.append(dict(point_file=path.relative_to(ROOT).as_posix(),map_point_file=row['point_file'],relative_frequency_error=err,reason='Initial preflight integer 0 versus map float 0.0 cache key',status='PASS'))
    a.csvwrite(folder/'repeated_preflight_checks.csv',duplicates)
    decision['repeated_preflight_prefixes']=len(duplicates)
    decision['workflow_sha256']=a.sha(__file__)
    decision['entry_point_sha256']=a.sha(ROOT/'scripts/analysis/study_two_springs_external_axial_spring.py')
    require(decision['new_full_prefixes']==decision['new_points']+len(duplicates),'Point accounting')
    a.write(folder/'decision.json',decision)
    a.write(folder/'manifest.json',dict(status='PASS',artifacts={p.relative_to(folder).as_posix():a.sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    return ranges.saved.diagnostic.verify_source(folder)


def main(args):
    if args.stage in ('plot_only','finalize'):
        require(args.folder is not None and args.folder.resolve().is_relative_to(BASE.resolve()),'Explicit new archive required')
        with previous.plots_without_calculation(),ranges.saved.forbid_full_calculations(leading_too=True):
            if args.stage=='plot_only':
                cfg=json.loads((args.folder/'contract.json').read_text())['config'];rows=read_csv(args.folder/'spectrum.csv');summary=analyze(args.folder,cfg,rows);mechanics(args.folder,cfg,rows,summary);render(args.folder)
            else:print(finalize(args.folder))
        return
    exp=Experiment();print(exp.folder,flush=True)
    try:
        require(not (exp.folder/'failure.json').exists(),'Saved failure needs addressed diagnosis')
        if args.stage=='preflight':exp.preflight()
        else:exp.run_study()
        exp.finish(args.stage)
    except Exception as exc:
        exp.a.write(exp.folder/'failure.json',dict(status='PARTIAL',stage=args.stage,error=repr(exc)));exp.finish(args.stage,'PARTIAL');raise


def mechanics(folder,cfg,rows,summary):
    """Read accepted mass-normalized fields; no new eigenpair or field solve."""
    m=physics.Model(**cfg['material']);out=[]
    for s in summary:
        for b in sorted({0.,20.,float(s['beta_max']),90.}):
            r=next(r for r in rows if r['pair']==s['pair'] and chi(r['chi_u'])==s['chi_u'] and float(r['beta_deg'])==b)
            path=ROOT/r['point_file'];require(api().sha(path.with_suffix('.npz'))==r['shape_sha256'],'Shape provenance')
            with np.load(path.with_suffix('.npz')) as a:
                arms=[{f:a[f'mode1_arm{i}_{f}'] for f in ('x','weights','y','dy','ends')} for i in (1,2)]
                c=physics.Case(b,cfg['mu'],*cfg['cases'][s['pair']],chi(s['chi_u']))
                en=physics.energy_components(m,c,arms);mass=physics.inner(m,dict(arms=arms),dict(arms=arms));om=float(r['omega1'])
                require(abs(mass-1)<=cfg['policy']['mass_gram_tol'] and abs(en.sum()/(mass*om**2)-1)<=cfg['policy']['energy_tol'],'Read-only energy/normalization check')
                out.append(dict(pair=r['pair'],chi_u=r['chi_u'],beta_deg=b,omega1=om,mass=mass,energy_error=abs(en.sum()/(mass*om**2)-1),u20=float(arms[1]['ends'][0,0]),axial_fraction=en[0]/en.sum(),spring_fraction=en[-1]/en.sum()))
    api().csvwrite(folder/'mechanics.csv',out)


def render(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    cfg=json.loads((folder/'contract.json').read_text())['config'];full=read_csv(folder/'spectrum.csv');deficit=read_csv(folder/'deficit.csv');summary=read_csv(folder/'summary.csv')
    plt.rcParams.update({'font.size':10,'axes.labelsize':11,'legend.fontsize':9,'pdf.fonttype':42,'font.family':'DejaVu Sans'})
    colors=['#222222','#0072B2','#D55E00','#009E73','#8B4A9C'];styles=['-','--','-.',':',(0,(5,1,1,1))]
    labels=['Точная продольная опора',r'$\chi_u=1$',r'$\chi_u=0{,}1$',r'$\chi_u=0{,}01$',r'$\chi_u=0$']
    comma=FuncFormatter(lambda x,pos:f'{x:g}'.replace('.',','))
    def group(rows,p,k):return [r for r in rows if r['pair']==p and chi(r['chi_u'])==k]
    def val(rows,n):return [float(r[n]) for r in rows]
    def num(x):return f'{x:g}'.replace('.',r'{,}')
    def title(p):
        j,e=cfg['cases'][p];return r'$\kappa_{\mathrm{j}}='+num(j)+r',\quad\kappa_{\mathrm{e}}='+num(e)+'$'
    def save(fig,axes,i):
        for ax in np.asarray(axes).flat:ax.grid(alpha=.2);ax.yaxis.set_major_formatter(comma)
        fig.savefig(folder/(FIGURES[i]+'.pdf'));fig.savefig(folder/(FIGURES[i]+'.png'),dpi=320);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(8.5,3.9),sharex=True,sharey=True,layout='constrained')
    for ax,p,letter in zip(axes,cfg['cases'],'аб'):
        for i,k in enumerate(cfg['chi_values']):
            rr=group(full,p,k);ax.plot(val(rr,'beta_deg'),val(rr,'Lambda1'),color=colors[i],ls=styles[i],lw=1.6,label=labels[i])
        ax.set_title(letter+') '+title(p),fontsize=10);ax.set_xlim(0,90);ax.set_xlabel(r'$\beta$, град.');ax.xaxis.set_major_formatter(comma)
    axes[0].set_ylabel(r'$\Lambda_1$');fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False);save(fig,axes,0)
    fig,axes=plt.subplots(1,2,figsize=(8.5,4.0),sharex=True,sharey=True,layout='constrained')
    for ax,p,letter in zip(axes,cfg['cases'],'аб'):
        for i,k in enumerate(cfg['chi_values']):
            rr=group(deficit,p,k);s=group(summary,p,k)[0]
            ax.plot(val(rr,'beta_deg'),val(rr,'R_percent'),color=colors[i],ls=styles[i],lw=1.5,label=labels[i])
            for d,marker in [(10,'s'),(5,'o')]:
                if s[f'beta_high_{d}']:
                    b=float(s[f'beta_high_{d}']);r=next(r for r in rr if float(r['beta_deg'])==b)
                    ax.plot(b,float(r['R_percent']),marker=marker,ms=4,mfc='white',mec=colors[i])
        for d in (5,10):
            ax.axhline(d,color='#555555',lw=.8,ls='--');ax.text(90.5,d,f'{d}%',va='center',fontsize=8)
        ax.set_title(letter+') '+title(p),fontsize=10);ax.set_xlim(0,90);ax.set_ylim(-1,77);ax.set_xlabel(r'$\beta$, град.');ax.xaxis.set_major_formatter(comma)
    axes[0].set_ylabel(r'$R$, %');fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False);save(fig,axes,1)
    order=[0.,.01,.1,1.,'RIGID'];x=np.arange(5);ticks=['0','0,01','0,1','1','Точная\nопора']
    fig,axes=plt.subplots(1,2,figsize=(8.5,4.0),sharey=True,layout='constrained')
    for ax,p,letter in zip(axes,cfg['cases'],'аб'):
        ss={chi(r['chi_u']):r for r in summary if r['pair']==p}
        for d,color,marker,ls in [(10,'#D55E00','s','--'),(5,'#0072B2','o','-')]:
            y=[float(ss[k][f'beta_high_{d}']) if ss[k][f'beta_high_{d}'] else np.nan for k in order]
            ax.plot(x,y,color=color,marker=marker,ls=ls,lw=1.5,label=f'Не менее {100-d}%')
        for i,k in enumerate(order):
            if any(ss[k][f'status_{d}']=='NO_TERMINAL_HIGH_RANGE' for d in (10,5)):
                ax.plot(i,.88,marker='x',color='#8B4A9C',ms=8,transform=ax.get_xaxis_transform());ax.text(i,.95,'Нет хвостового\nдиапазона',transform=ax.get_xaxis_transform(),ha='left',va='center',fontsize=8,color='#8B4A9C')
        ax.set_title(letter+') '+title(p),fontsize=10);ax.set_xticks(x,ticks);ax.set_xlabel(r'$\chi_u$');ax.set_ylim(0,65)
    axes[0].set_ylabel(r'$\beta_{\mathrm{high}}^{\mathrm{grid}}$, град.');fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=2,frameon=False);save(fig,axes,2)
    fig,axes=plt.subplots(2,2,figsize=(8.5,5.5),sharex=True,sharey='row',layout='constrained')
    for col,p in enumerate(cfg['cases']):
        ss={chi(r['chi_u']):r for r in summary if r['pair']==p}
        for row,n in enumerate(('omega_max','G_max_percent')):
            axes[row,col].plot(x,[float(ss[k][n]) for k in order],color='#0072B2',marker='o',lw=1.5)
        axes[0,col].set_title(title(p),fontsize=10);axes[1,col].set_xticks(x,ticks);axes[1,col].set_xlabel(r'$\chi_u$')
    axes[0,0].set_ylabel(r'$\omega_{\max}^{\mathrm{grid}}$');axes[1,0].set_ylabel(r'$G_{\max}^{\mathrm{grid}}$, %');save(fig,axes,3)
    api().write(folder/'figure_provenance.json',[dict(figure=name,data_source=source,dpi=320,curves='segments between verified nodes only; no frequency interpolation') for name,source in zip(FIGURES,['spectrum.csv','deficit.csv','summary.csv','summary.csv'])])
    api().write(folder/'rendering.json',dict(status='PASS',new_spectral_solves=0,dpi=320))
