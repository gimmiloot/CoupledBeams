"""Bounded small-chi preset; inherited full solver, physics, gates and ranges.

The archive endpoints are read-only; no new physical or leading assembly.
"""
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
import scipy
from scripts.lib import two_springs_external_axial_high_ranges as previous

ROOT=previous.ROOT
BASE=previous.BASE.with_name('external_axial_small_stiffness_transition')
CONFIG=ROOT/'data/input/isotropic_two_springs_external_axial_small_stiffness_transition.json'
physics=previous.physics
api=previous.api
require=previous.require
read_csv=previous.read_csv
chi=previous.chi
FIGURES=('01_first_frequency','02_signed_growth','03_high_range_boundaries','04_maximum_and_gain')


class Experiment(previous.Experiment):
    def __init__(self):
        self.a=api();self.cfg=json.loads(CONFIG.read_text());self.inputs={};self.sources={};self.manifests={}
        self.started=time.perf_counter();self.ritz_new=0;self.leading_new=0
        for name,digest in self.cfg['frozen_files'].items():self.check(ROOT/name,digest)
        for name,r in self.cfg['sources'].items():
            folder=ROOT/r['folder'];self.check(folder/'manifest.json',r['manifest_sha256']);previous.ranges.saved.diagnostic.verify_source(folder)
            self.sources[name]=folder;self.manifests[name]=json.loads((folder/'manifest.json').read_text())
            c=json.loads(self.source(name,'contract.json').read_text())['config']
            require(self.manifests[name]['status']=='PASS' and c['material']==self.cfg['material'] and c['mu']==self.cfg['mu'] and c['cases']==self.cfg['cases'] and c['variant']==self.cfg['variant'],'Source physical identity')
            require(c['map']['spectrum_semantics']=='sorted_positions' and c['map']['K_guard']==3,'Source semantics')
        identity=dict(config=self.cfg,physical_sha256=self.a.sha(physics.__file__),versions=dict(numpy=np.__version__,scipy=scipy.__version__))
        self.key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16];self.folder=BASE/self.key;self.folder.mkdir(parents=True,exist_ok=True)
        target=self.folder/'contract.json'
        if target.exists():require(json.loads(target.read_text())==identity,'Immutable contract')
        else:self.a.write(target,identity)
        self.run=self.a.Run(self.cfg,self.folder,engine=physics)
        self.run.catalogs={float(k):v for k,v in json.loads(self.source('ranges','pole_catalogs.json').read_text()).items()}
        self.a.write(self.folder/'pole_catalogs.json',self.run.catalogs)
        self.archive={};self.overlap=[]
        for name in ('ranges','axial'):
            for row in read_csv(self.source(name,'spectrum.csv')):
                if chi(row['chi_u']) in self.cfg['archived_chi_values']:self.ingest(ROOT/row['point_file'],row['point_sha256'],name,row)
        require(len(self.archive)==426,'Unexpected endpoint inventory')
        inventory=[]
        for p in self.cfg['cases']:
            for k in self.cfg['chi_values']:
                have=sorted(b for pp,kk,b in self.archive if (pp,kk)==(p,k))
                if k in self.cfg['archived_chi_values']:require(have[0]==0 and have[-1]==90,'Endpoint coverage')
                inventory.append(dict(pair=p,chi_u=k,available_beta_deg=have,missing_main_beta_deg=[b for b in self.cfg['map']['parameter_grid'] if b not in have]))
        self.a.write(self.folder/'archive_inventory.json',inventory);self.csv('archive_overlaps',self.overlap)
        self.plan=json.loads((self.folder/'refinement_plan.json').read_text()) if (self.folder/'refinement_plan.json').exists() else []
        self.allowed={}
        for r in self.plan:self.allowed.setdefault((r['pair'],r['chi_u']),set()).update(r['new_beta_deg'])
        for p in self.cfg['cases']:
            for k in self.cfg['new_chi_values']:self.allowed.setdefault((p,k),set()).update(self.cfg['preflight_angles'])
        scale=[dict(chi_u=k,kappa_u_b=k/self.cfg['material']['epsilon']**2) for k in self.cfg['chi_values']]
        require(np.allclose([r['kappa_u_b'] for r in scale],self.cfg['stiffness_scale']['expected_kappa_b'],rtol=1e-14,atol=0),'Stiffness scale')
        self.csv('stiffness_scale',scale)

    def point(self,pair,k,beta,previous_w=None):
        if (pair,k,beta) in self.archive:return self.archive[pair,k,beta]
        require(k in self.cfg['new_chi_values'],'Archived or unrequested stiffness cannot be recomputed')
        return super().point(pair,k,beta,previous_w)

    def preflight(self):
        target=self.folder/'preflight.json'
        if target.exists():require(json.loads(target.read_text())['status']=='PASS','Saved preflight failure');return
        rows=[];ritz=[];reg=[]
        for p in self.cfg['cases']:
            ref=self.archive[p,0.,0.][0]['roots'][0]['omega']
            for k in self.cfg['new_chi_values']:
                for b in self.cfg['preflight_angles']:
                    rec,file,_=self.point(p,k,b);w=np.array([r['omega'] for r in rec['roots']])
                    rows.append(dict(pair=p,chi_u=k,beta_deg=b,Lambda1=rec['roots'][0]['Lambda'],Lambda3=rec['roots'][2]['Lambda'],omega1=w[0],relative_to_archived_zero=w[0]/ref-1,gap12=w[1]/w[0]-1,lower_count=rec['count']['samples'][0]['count'],upper_count=rec['count']['samples'][1]['count'],point_file=file.relative_to(ROOT).as_posix(),status='PASS'))
                    if b==0:
                        err=abs(w[0]/ref-1);reg.append(dict(pair=p,chi_u=k,relative_omega_error=err,Lambda1=rec['roots'][0]['Lambda'],status='PASS'));require(err<=self.cfg['policy']['reference_frequency_tol'],'Zero-angle regression failed')
                    if [k,b] in self.cfg['ritz_cases']:
                        for order in self.cfg['policy']['ritz_orders']:
                            wr,_,d=physics.ritz(self.run.model,self.case(p,k,b),order);self.ritz_new+=1;err=float(max(abs(wr/w-1)))
                            ritz.append(dict(pair=p,chi_u=k,beta_deg=b,order=order,relative_frequency_error=err,**d));self.csv('preflight_ritz',ritz)
                            require(err<=self.cfg['policy']['ritz_frequency_tol'] and d['positive_K_M'] and d['energy_error']<=self.cfg['policy']['energy_tol'] and d['mass_gram_error']<=self.cfg['policy']['mass_gram_tol'] and d['constraint_residual']<=self.cfg['policy']['boundary_tol'],'Small-stiffness Ritz gate')
                    self.csv('preflight_points',rows);self.csv('zero_angle_regression',reg)
        self.a.write(target,dict(status='PASS',points=len(rows),new_points=self.run.new_points,ritz_solves=len(ritz),max_ritz_error=max(r['relative_frequency_error'] for r in ritz),zero_angle_relative_error=max(r['relative_omega_error'] for r in reg),min_Lambda1=min(r['Lambda1'] for r in rows),max_Lambda3=max(r['Lambda3'] for r in rows),min_gap12=min(r['gap12'] for r in rows),zero_mechanism=False))

    def refine(self,rows,round_number):
        cfg=self.cfg
        try:
            self.cfg={**cfg,'chi_values':cfg['new_chi_values']}
            return super().refine(rows,round_number)
        finally:self.cfg=cfg

    def comparisons(self,rows,summary):
        old=read_csv(self.source('ranges','summary.csv'));reg=[];mon=[];limit=[];transition=[]
        by={(r['pair'],r['chi_u'],r['beta_deg']):r for r in rows}
        for s in summary:
            if s['chi_u'] in self.cfg['archived_chi_values']:
                o=next(r for r in old if r['pair']==s['pair'] and chi(r['chi_u'])==s['chi_u'])
                for key in ('omega0','omega_max','beta_max','G_max_percent','beta_high_10','beta_high_5'):
                    require((float(o[key]) if o[key] else None)==s[key],'Archived summary changed')
                reg.append(dict(pair=s['pair'],chi_u=s['chi_u'],status='PASS'))
        for p in self.cfg['cases']:
            common=set.intersection(*[{b for pp,kk,b in by if (pp,kk)==(p,k)} for k in self.cfg['chi_values']])
            for b in sorted(common):
                free=by[p,0.,b];w0=float(free['omega1'])
                for lo,hi in zip(self.cfg['chi_values'],self.cfg['chi_values'][1:]):
                    delta=[by[p,hi,b][f'omega{i}']/by[p,lo,b][f'omega{i}']-1 for i in (1,2,3)]
                    require(min(delta)>=-self.cfg['policy']['same_physics_frequency_tol'],'Stiffness monotonicity')
                    mon.append(dict(pair=p,beta_deg=b,chi_lower=lo,chi_upper=hi,min_relative_increment=min(delta),status='PASS'))
                for k in self.cfg['chi_values'][1:]:limit.append(dict(pair=p,beta_deg=b,chi_u=k,omega_free=w0,omega1=by[p,k,b]['omega1'],relative_to_free=by[p,k,b]['omega1']/w0-1))
            for b in self.cfg['transition_angles']:
                vals=[(k,by[p,k,b]['omega1']/by[p,k,0.]['omega1']-1) for k in self.cfg['chi_values']]
                for (lo,dl),(hi,dh) in zip(vals,vals[1:]):
                    if dl<0<dh:transition.append(dict(pair=p,beta_deg=b,chi_lower=lo,chi_upper=hi,kappa_b_lower=lo/self.cfg['material']['epsilon']**2,kappa_b_upper=hi/self.cfg['material']['epsilon']**2,D0_lower=dl,D0_upper=dh,meaning='sign bracket at this fixed angle; not a universal critical stiffness'))
        for n,data in [('archive_range_regression',reg),('stiffness_monotonicity',mon),('free_limit_main_grid',limit),('transition_brackets',transition)]:self.csv(n,data)

    def run_study(self):
        self.preflight();rows=self.sweep()
        if not (self.folder/'initial_summary.csv').exists():self.analyze(rows,'initial_')
        for n in range(1,self.cfg['refinement']['max_rounds']+1):
            if not self.refine(rows,n):break
            rows=self.sweep()
        else:require(not self.refine(rows,n+1),'Unresolved refinement')
        summary=self.analyze(rows);self.comparisons(rows,summary);details(self.folder,self.cfg,rows,summary)
        self.a.write(self.folder/'calculation_complete.json',dict(status='PASS',constructions=len(summary),points=len(rows),reused_points=sum(r['origin']!='new_study' for r in rows),new_points=sum(r['origin']=='new_study' for r in rows),refinement_rounds=max((r['round'] for r in self.plan),default=0)))

    def finish(self,stage,status='PASS'):
        super().finish(stage,status)
        path=self.folder/'executions.json';runs=json.loads(path.read_text());runs[-1]['workflow_sha256']=self.a.sha(__file__);self.a.write(path,runs)

    def addressed(self):
        """Six fixed-angle sign brackets and one smaller chi, never a new map."""
        target=self.folder/'addressed_complete.json'
        if target.exists():require(json.loads(target.read_text())['status']=='PASS','Saved addressed failure');return
        rows=read_csv(self.folder/'spectrum.csv');by={(r['pair'],chi(r['chi_u']),float(r['beta_deg'])):r for r in rows}
        brackets=read_csv(self.folder/'transition_brackets.csv')
        initial=read_csv(self.folder/'free_limit_main_grid.csv')
        worst={p:max(float(r['relative_to_free']) for r in initial if r['pair']==p and float(r['chi_u'])==min(self.cfg['new_chi_values'])) for p in self.cfg['cases']}
        plan=dict(defined_before_addressed_solves=True,main_spectrum_sha256=self.a.sha(self.folder/'spectrum.csv'),
                  sign_brackets=brackets,log10_bisections_per_bracket=3,angles=self.cfg['transition_angles'],
                  smaller_chi=1e-6,reason='The smallest main stiffness still differs from the free frequency by the measured amounts; one decade smaller tests contraction.',smallest_chi_max_relative_to_free=worst,full_additional_maps=False)
        self.a.write(self.folder/'addressed_plan.json',plan)
        added={};requests=[];free=[];final=[];history=[];quality=[]
        def point(p,k,b,reason):
            key=p,float(k),float(b)
            if key in by:return float(by[key]['omega1'])
            if key in added:return added[key]['omega1']
            require(b in [0.]+self.cfg['transition_angles'] and k>0,'Addressed scope')
            requests.append(dict(pair=p,chi_u=k,beta_deg=b,reason=reason));self.a.write(self.folder/'addressed_requests.json',requests)
            c=self.case(p,k,b);rec=self.run.solve(c);ww=self.validate(rec,c);path=previous.previous.lengths.point_path(self.folder,c)
            row=dict(pair=p,chi_u=k,beta_deg=b,reason=reason,**{f'{n}{i}':r[n] for i,r in enumerate(rec['roots'],1) for n in ('omega','Lambda')},point_file=path.relative_to(ROOT).as_posix(),point_sha256=self.a.sha(path),shape_sha256=rec['shape_sha256'])
            if b==0:require(abs(ww[0]/float(by[p,0.,0.]['omega1'])-1)<=self.cfg['policy']['reference_frequency_tol'],'Addressed zero-angle regression')
            added[key]=row
            for i,r in enumerate(rec['roots'],1):
                d=r['diagnostics'];quality.append(dict(pair=p,chi_u=k,beta_deg=b,position=i,physical=float(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))),ode=d['ode_error'],energy=d['energy_error'],singular=d['singular_ratio'],mass_gram=rec['mass_gram_error'],condition=d['nonzero_condition'],count_lower=r['counts'][0],count_upper=r['counts'][1],status='PASS'))
            self.csv('addressed_spectrum',list(added.values()));self.csv('addressed_quality',quality)
            return ww[0]
        for r in brackets:
            p=r['pair'];b=float(r['beta_deg']);lo=float(r['chi_lower']);hi=float(r['chi_upper']);dl=float(r['D0_lower']);dh=float(r['D0_upper'])
            for step in range(3):
                mid=math.sqrt(lo*hi);d=point(p,mid,b,'logarithmic sign bracket')/point(p,mid,0.,'own zero-angle normalization')-1
                require(abs(d)>self.cfg['policy']['same_physics_frequency_tol'],'Sign unresolved at numerical accuracy')
                history.append(dict(pair=p,beta_deg=b,step=step+1,chi_lower=lo,chi_upper=hi,chi_trial=mid,D0_trial=d))
                if d>0:hi,dh=mid,d
                else:lo,dl=mid,d
            final.append(dict(pair=p,beta_deg=b,chi_lower=lo,chi_upper=hi,kappa_b_lower=lo/self.cfg['material']['epsilon']**2,kappa_b_upper=hi/self.cfg['material']['epsilon']**2,D0_lower=dl,D0_upper=dh,log10_width=math.log10(hi/lo)))
        for p in self.cfg['cases']:
            for b in self.cfg['transition_angles']:
                r=by[p,0.,b];path=ROOT/r['point_file'];require(self.a.sha(path.with_suffix('.npz'))==r['shape_sha256'],'Free shape provenance')
                with np.load(path.with_suffix('.npz')) as data:u20=float(data['mode1_arm2_ends'][0,0])
                w0=float(r['omega1']);old=float(by[p,1e-5,b]['omega1']);new=point(p,1e-6,b,'one smaller stiffness for the free-limit check')
                upper=math.sqrt(w0*w0+1e-6*self.run.model.section.EA/self.run.model.l*u20*u20)
                require(new>=w0*(1-self.cfg['policy']['same_physics_frequency_tol']) and new<=upper*(1+self.cfg['policy']['same_physics_frequency_tol']) and abs(new/w0-1)<abs(old/w0-1),'Free-limit contraction/variational upper bound')
                free.append(dict(pair=p,beta_deg=b,omega_free=w0,omega_1e_5=old,omega_1e_6=new,relative_1e_5=old/w0-1,relative_1e_6=new/w0-1,contraction=(new/w0-1)/(old/w0-1),free_form_Rayleigh_upper=upper,status='PASS'))
        for n,v in [('transition_refined_brackets',final),('transition_refinement_history',history),('free_limit_addressed',free)]:self.csv(n,v)
        self.a.write(target,dict(status='PASS',unique_prefixes=len(added),sign_brackets=len(final),free_limit_points=len(free),full_additional_maps=0))


def details(folder,cfg,rows,summary):
    growth=[];extrema=[];shape=[]
    for s in summary:
        rr=sorted([r for r in rows if r['pair']==s['pair'] and chi(r['chi_u'])==s['chi_u']],key=lambda r:float(r['beta_deg']))
        b=np.array([float(r['beta_deg']) for r in rr]);w=np.array([float(r['omega1']) for r in rr]);d=w/w[0]-1
        growth.extend(dict(pair=s['pair'],chi_u=s['chi_u'],beta_deg=x,omega1=y,F=f+1,D0=f,D0_percent=100*f) for x,y,f in zip(b,w,d))
        resolved=[]
        for i in range(1,len(b)-1):
            dl=w[i]-w[i-1];dr=w[i+1]-w[i]
            if dl*dr<0 and min(abs(dl),abs(dr))/w[i]>cfg['refinement']['extremum_frequency_rtol']:
                kind='maximum' if dl>0 else 'minimum';resolved.append(kind)
                extrema.append(dict(pair=s['pair'],chi_u=s['chi_u'],kind=kind,beta_deg=b[i],omega1=w[i],left_deg=b[i-1],right_deg=b[i+1]))
        scale=cfg['policy']['same_physics_frequency_tol'];diff=np.diff(w)/w[:-1]
        shape.append(dict(pair=s['pair'],chi_u=s['chi_u'],min_D0_percent=100*min(d),max_D0_percent=100*max(d),D90_percent=100*d[-1],sampled_nondecreasing=bool(min(diff)>=-scale),sampled_nonincreasing=bool(max(diff)<=scale),resolved_interior_extrema=len(resolved),first_resolved_positive_angle=next((float(x) for x,y in zip(b,d) if y>scale),None)))
    for n,data in [('normalized_growth',growth),('interior_extrema',extrema),('curve_character',shape)]:api().csvwrite(folder/(n+'.csv'),data)
    previous.mechanics(folder,cfg,rows,summary)


def main(args):
    if args.stage in ('plot_only','finalize'):
        require(args.folder is not None and args.folder.resolve().is_relative_to(BASE.resolve()),'Explicit new archive required')
        with previous.previous.plots_without_calculation(),previous.ranges.saved.forbid_full_calculations(leading_too=True):
            if args.stage=='plot_only':render(args.folder)
            else:print(finalize(args.folder))
        return
    exp=Experiment();print(exp.folder,flush=True)
    try:
        require(not (exp.folder/'failure.json').exists(),'Saved failure needs diagnosis')
        if args.stage=='preflight':exp.preflight()
        else:exp.run_study();exp.addressed()
        exp.finish(args.stage)
    except Exception as exc:
        exp.a.write(exp.folder/'failure.json',dict(status='PARTIAL',stage=args.stage,error=repr(exc)));exp.finish(args.stage,'PARTIAL');raise


def finalize(folder):
    folder=folder.resolve();a=api();cfg=json.loads((folder/'contract.json').read_text())['config']
    done=json.loads((folder/'calculation_complete.json').read_text());addressed=json.loads((folder/'addressed_complete.json').read_text())
    require(done['status']==addressed['status']=='PASS' and not (folder/'failure.json').exists(),'Incomplete study')
    for n,digest in json.loads((folder/'source_hashes.json').read_text()).items():require(a.sha(ROOT/n)==digest,'Protected source changed')
    a.write(folder/'source_integrity.json',{n:previous.ranges.saved.diagnostic.verify_source(ROOT/r['folder']) for n,r in cfg['sources'].items()})
    for n in FIGURES:
        for ext in ('.pdf','.png'):require((folder/(n+ext)).exists(),'Missing figure')
    runs=json.loads((folder/'executions.json').read_text());q=read_csv(folder/'quality.csv');qa=read_csv(folder/'addressed_quality.csv');s=read_csv(folder/'summary.csv')
    prefixes=sum(r['new_full_prefixes'] for r in runs)
    require(prefixes==done['new_points']+addressed['unique_prefixes']==len(list((folder/'points').glob('*.json'))),'Unique point accounting')
    allq=q+qa
    addressed_rows=read_csv(folder/'addressed_spectrum.csv');full=read_csv(folder/'spectrum.csv')+addressed_rows;mon=[];growth=[]
    for p in cfg['cases']:
        for b in cfg['transition_angles']:
            rr=sorted([r for r in full if r['pair']==p and float(r['beta_deg'])==b],key=lambda r:float(r['chi_u']))
            for lo,hi in zip(rr,rr[1:]):
                delta=min(float(hi[f'omega{i}'])/float(lo[f'omega{i}'])-1 for i in (1,2,3))
                require(delta>=-cfg['policy']['same_physics_frequency_tol'],'Addressed spectral monotonicity')
                mon.append(dict(pair=p,beta_deg=b,chi_lower=lo['chi_u'],chi_upper=hi['chi_u'],min_relative_increment=delta,status='PASS'))
        zero={float(r['chi_u']):float(r['omega1']) for r in addressed_rows if r['pair']==p and float(r['beta_deg'])==0}
        for r in addressed_rows:
            k=float(r['chi_u'])
            if r['pair']==p and k in zero and float(r['beta_deg'])>0:growth.append(dict(pair=p,chi_u=k,beta_deg=float(r['beta_deg']),omega1=float(r['omega1']),omega0=zero[k],D0_percent=100*(float(r['omega1'])/zero[k]-1)))
    a.csvwrite(folder/'addressed_monotonicity.csv',mon);a.csvwrite(folder/'addressed_growth.csv',growth)
    for r in read_csv(folder/'high_intervals.csv'):
        for key in ('entry_resolution_deg','exit_resolution_deg'):require(not r[key] or float(r[key])<=.25+1e-10,'Unresolved threshold bracket')
    d=dict(**done,addressed_prefixes=addressed['unique_prefixes'],new_full_prefixes=prefixes,new_full_roots=sum(r['new_full_roots'] for r in runs),new_pole_roots=sum(r['new_pole_roots'] for r in runs),independent_ritz_solves=sum(r['new_ritz_solves'] for r in runs),new_leading_solves=0,quality_maxima={n:max(float(r[n]) for r in allq) for n in ('physical','ode','energy','singular','mass_gram','condition')},no_terminal_cases={str(v):sum(r[f'status_{v}']=='NO_TERMINAL_HIGH_RANGE' for r in s) for v in cfg['levels_percent']},monotonicity_rows=len(read_csv(folder/'stiffness_monotonicity.csv')),workflow_sha256=a.sha(__file__),entry_point_sha256=a.sha(ROOT/'scripts/analysis/study_two_springs_external_axial_spring.py'))
    a.write(folder/'decision.json',d)
    a.write(folder/'manifest.json',dict(status='PASS',artifacts={p.relative_to(folder).as_posix():a.sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    return previous.ranges.saved.diagnostic.verify_source(folder)


def render(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    from matplotlib.lines import Line2D
    started=time.perf_counter();cfg=json.loads((folder/'contract.json').read_text())['config'];rows=read_csv(folder/'spectrum.csv');growth=read_csv(folder/'normalized_growth.csv');s=read_csv(folder/'summary.csv');intervals=read_csv(folder/'high_intervals.csv')
    plt.rcParams.update({'font.size':10,'axes.labelsize':11,'legend.fontsize':9,'font.family':'DejaVu Sans','pdf.fonttype':42})
    colors=['#222222','#0072B2','#56B4E9','#D55E00','#009E73','#8B4A9C'];styles=['-', '--', '-.', ':', (0,(5,1,1,1)), '-']
    labels=[r'$\chi_u=0$',r'$10^{-5}$',r'$2{,}5\cdot10^{-5}$',r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$']
    comma=FuncFormatter(lambda x,pos:f'{x:g}'.replace('.',','))
    def group(data,p,k):return [r for r in data if r['pair']==p and chi(r['chi_u'])==k]
    def title(p):
        j,e=cfg['cases'][p];num=lambda v:f'{v:g}'.replace('.',r'{,}')
        return r'$\kappa_{\mathrm{j}}='+num(j)+r',\quad\kappa_{\mathrm{e}}='+num(e)+'$'
    def save(fig,axes,n):
        for ax in np.asarray(axes).flat:ax.grid(alpha=.2);ax.yaxis.set_major_formatter(comma)
        fig.savefig(folder/(FIGURES[n]+'.pdf'));fig.savefig(folder/(FIGURES[n]+'.png'),dpi=320);plt.close(fig)
    for n,data,ykey,ylabel in [(0,rows,'Lambda1',r'$\Lambda_1$'),(1,growth,'D0_percent',r'$100(F-1)$, %')]:
        fig,axes=plt.subplots(1,2,figsize=(8.5,4.1),sharex=True,sharey=True,layout='constrained')
        for ax,p,letter in zip(axes,cfg['cases'],'аб'):
            for i,k in enumerate(cfg['chi_values']):
                rr=group(data,p,k);ax.plot([float(r['beta_deg']) for r in rr],[float(r[ykey]) for r in rr],color=colors[i],ls=styles[i],lw=1.5,label=labels[i])
            if n==1:ax.axhline(0,color='gray',lw=.7)
            ax.set_title(letter+') '+title(p),fontsize=10);ax.set_xlim(0,90);ax.set_xlabel(r'$\beta$, град.')
        axes[0].set_ylabel(ylabel);fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False);save(fig,axes,n)
    fig,axes=plt.subplots(1,2,figsize=(8.5,4.8),sharex=True,sharey=True,layout='constrained')
    for ax,p,letter in zip(axes,cfg['cases'],'аб'):
        for i,k in enumerate(cfg['chi_values']):
            for d,color,off in [(10,'#D55E00',-.11),(5,'#0072B2',.11)]:
                for r in group(intervals,p,k):
                    if float(r['delta_percent'])!=d:continue
                    a,b=float(r['start_deg']),float(r['end_deg']);terminal=r['terminal']=='True'
                    ax.plot([a,b],[i+off,i+off],color=color,lw=2,ls='-' if terminal else '--',marker='|' if terminal else 'o',ms=4)
        ax.axhline(.5,color='gray',lw=.8,ls=':');ax.set_yticks(range(6),labels);ax.set_xlim(-1,91);ax.set_ylim(5.5,-.5);ax.set_xlabel(r'$\beta$, град.');ax.set_title(letter+') '+title(p),fontsize=10);ax.grid(axis='x',alpha=.2)
    handles=[Line2D([],[],color='#D55E00',lw=2,label='Не менее 90%'),Line2D([],[],color='#0072B2',lw=2,label='Не менее 95%'),Line2D([],[],color='gray',ls='--',label='Начальный участок без хвоста')]
    fig.legend(handles=handles,loc='outside lower center',ncol=2,frameon=False)
    # Keep categorical chi labels; the common save helper formats numeric y only.
    fig.savefig(folder/(FIGURES[2]+'.pdf'));fig.savefig(folder/(FIGURES[2]+'.png'),dpi=320);plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(9,5.7),sharey='row',gridspec_kw={'width_ratios':[.6,2.5,.6,2.5]},layout='constrained')
    for j,p in enumerate(cfg['cases']):
        ss={chi(r['chi_u']):r for r in s if r['pair']==p};positive=cfg['chi_values'][1:]
        for row,key in enumerate(('omega_max','G_max_percent')):
            zero=axes[row,2*j];ax=axes[row,2*j+1]
            zero.plot([0],[float(ss[0.][key])],marker='o',color='#222222');zero.set_xlim(-.4,.4);zero.set_xticks([0],['0']);zero.spines['right'].set_visible(False)
            ax.semilogx(positive,[float(ss[k][key]) for k in positive],color='#0072B2',marker='o',lw=1.5);ax.set_xlim(8e-6,.013);ax.set_xticks([1e-5,1e-4,1e-3,.01],[r'$10^{-5}$',r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$']);ax.spines['left'].set_visible(False)
        axes[0,2*j+1].set_title(title(p),fontsize=10);axes[1,2*j].set_xlabel(r'$\chi_u=0$');axes[1,2*j+1].set_xlabel(r'$\chi_u>0$')
    axes[0,0].set_ylabel(r'$\omega_{\max}^{\mathrm{grid}}$');axes[1,0].set_ylabel(r'$G_{\max}^{\mathrm{grid}}$, %');save(fig,axes,3)
    api().write(folder/'figure_provenance.json',[dict(figure=n,data_source=src,dpi=320,semantics='verified discrete nodes only') for n,src in zip(FIGURES,['spectrum.csv','normalized_growth.csv','high_intervals.csv','summary.csv'])])
    api().write(folder/'rendering.json',dict(seconds=time.perf_counter()-started,new_spectral_solves=0,dpi=320))
