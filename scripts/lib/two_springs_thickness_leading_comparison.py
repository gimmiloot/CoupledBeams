"""Addressed comparison; no new mechanical or spectral solver.

Full frequencies are immutable archive observations. Only physics.Leading.solve
is allowed to calculate spectra; all three sorted Ritz positions are checked.
"""
from contextlib import contextmanager
import hashlib
import json
import math
from pathlib import Path
import time
from types import SimpleNamespace
import numpy as np
import scipy
from scripts.lib import two_springs_thickness_asymmetry as previous
from scripts.lib import isotropic_two_springs_thickness_asymmetry as physics

ROOT=previous.ROOT
CONFIG=ROOT/'data/input/isotropic_two_springs_thickness_asymmetry_leading_comparison.json'
BASE=previous.lengths.BASE/'thickness_asymmetry_leading_comparison'
require=previous.require
read_csv=previous.read_csv
api=previous.api
FIGURES=('01_normalized_frequency','02_model_discrepancy','03_thickness_sign_effect')
QUALITY=('symmetry','eigen_residual','energy','mass_gram','constraint','natural_boundary','strong_ode','essential')


def atomic(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp');api().write(temporary,data);temporary.replace(path)


@contextmanager
def forbid_full_calculations():
    def forbidden(*args,**kwargs):raise RuntimeError('Full-model calculation is forbidden in this comparison')
    targets=[(m,n) for m in (physics,physics.base) for n in ('prefix','mode','ritz')]+[(api().Run,'solve')]
    saved=[getattr(m,n) for m,n in targets]
    try:
        for m,n in targets:setattr(m,n,forbidden)
        yield
    finally:
        for (m,n),f in zip(targets,saved):setattr(m,n,f)


def frequency_comparison(full_lambda,full_zero,leading_zeta,leading_zero):
    F_full=(full_lambda/full_zero)**2;F_lim=math.sqrt(leading_zeta/leading_zero)
    e=F_full/F_lim-1
    return dict(F_full=F_full,F_lim=F_lim,e=e,delta_percent=100*abs(e))


def sign_effect(plus,minus):
    full=plus['F_full']/minus['F_full']-1;lim=plus['F_lim']/minus['F_lim']-1
    return dict(D_full=full,D_lim=lim,Delta_D=lim-full,Delta_D_pp=100*(lim-full))


class Experiment:
    def __init__(self):
        self.cfg=json.loads(CONFIG.read_text());self.hashes={};self.sources={};self.manifests={}
        self.new=[];self.reused=0;self.started=time.perf_counter();self.render_seconds=0.
        for path,digest in self.cfg['frozen_files'].items():self.check(ROOT/path,digest)
        for name,source in self.cfg['sources'].items():
            path=ROOT/source['folder'];self.check(path/'manifest.json',source['manifest_sha256'])
            manifest=json.loads((path/'manifest.json').read_text());require(manifest['status']=='PASS','Source status')
            self.sources[name]=path;self.manifests[name]=manifest
        # Reuse the established complete-bundle and per-root verification.
        previous.verify_bundle(self.sources['thickness'])
        source_cfg=json.loads(self.source('thickness','contract.json').read_text())['config']
        require(source_cfg['material']==self.cfg['material'] and source_cfg['cases']==self.cfg['cases'],'Material or springs differ')
        require(source_cfg['variant']=='RL_FULL' and source_cfg['map']['spectrum_semantics']=='sorted_positions','Wrong full model')
        self.full=[];self.archive_checks=[]
        validator=SimpleNamespace(run=SimpleNamespace(policy=source_cfg['policy'],model=physics.base.Model(**self.cfg['material'])))
        for r in read_csv(self.source('thickness','spectrum.csv')):
            pair,mu,tau,beta=r['pair'],float(r['mu']),float(r['tau']),float(r['beta_deg'])
            xi=math.radians(beta)/self.cfg['material']['epsilon']
            if pair not in self.cfg['cases'] or mu not in self.cfg['mu_values'] or tau not in self.cfg['tau_values'] or not 0<=xi<=self.cfg['xi_max']:continue
            path=ROOT/r['point_file'];self.check(path,r['point_sha256']);self.check(path.with_suffix('.npz'),r['shape_sha256'])
            rec=json.loads(path.read_text());case=physics.Case(beta,mu,*self.cfg['cases'][pair],tau)
            previous.lengths.Experiment.validate(validator,rec,case)
            require(all(abs(float(r[f'omega{k}'])/p['omega']-1)<1e-14 and abs(float(r[f'Lambda{k}'])/p['Lambda']-1)<1e-14 for k,p in enumerate(rec['roots'],1)),'CSV/root mismatch')
            self.full.append(dict(pair=pair,mu=mu,tau=tau,beta_deg=beta,xi=xi,
                **{f'{name}{k}':float(r[f'{name}{k}']) for name in ('omega','Lambda') for k in (1,2,3)},
                g12_full=float(r['gap12_relative']),point_file=r['point_file'],point_sha256=r['point_sha256'],shape_sha256=r['shape_sha256']))
            self.archive_checks.append(dict(pair=pair,mu=mu,tau=tau,beta_deg=beta,xi=xi,
                physical=max(max(np.max(p['diagnostics']['joint_scaled']),np.max(p['diagnostics']['external_scaled'])) for p in rec['roots']),
                ode=max(p['diagnostics']['ode_error'] for p in rec['roots']),energy=max(p['diagnostics']['energy_error'] for p in rec['roots']),
                count_status=rec['count']['status'],status='PASS'))
        groups={(p,m,t) for p in self.cfg['cases'] for m in self.cfg['mu_values'] for t in self.cfg['tau_values']}
        require(set((r['pair'],r['mu'],r['tau']) for r in self.full)==groups,'Missing construction')
        grids=[sorted(r['beta_deg'] for r in self.full if (r['pair'],r['mu'],r['tau'])==g) for g in sorted(groups)]
        require(all(g==grids[0] and g[0]==0 and len(g)==len(set(g)) for g in grids),'Mismatched or duplicate archive nodes')
        self.grid=grids[0];self.groups=sorted(groups)
        identity=dict(config=self.cfg,beta_grid=self.grid,xi_grid=[math.radians(b)/self.cfg['material']['epsilon'] for b in self.grid],
                      versions=dict(numpy=np.__version__,scipy=scipy.__version__))
        self.key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16]
        self.folder=BASE/self.key;self.folder.mkdir(parents=True,exist_ok=True)
        path=self.folder/'contract.json'
        if path.exists():require(json.loads(path.read_text())==identity,'Contract changed')
        else:atomic(path,identity)
        self.csv('reused_full_points',self.full);self.csv('archive_quality',self.archive_checks)
        self.csv('angle_nodes',[dict(beta_deg=b,xi=x) for b,x in zip(self.grid,identity['xi_grid'])])
        atomic(self.folder/'source_gate.json',dict(status='PASS',constructions=len(groups),points=len(self.full),nodes_per_construction=len(self.grid),
            beta_max_deg=max(self.grid),xi_max=max(identity['xi_grid']),excluded_invalid_points=[],full_model_solves=0))

    def check(self,path,expected):
        actual=api().sha(path);require(actual==expected,'Changed source '+str(path))
        self.hashes[path.relative_to(ROOT).as_posix()]=actual
        return path

    def source(self,family,name):return self.check(self.sources[family]/name,self.manifests[family]['artifacts'][name])
    def csv(self,name,rows):api().csvwrite(self.folder/(name+'.csv'),rows)

    def solve(self,pair,mu,tau,xi,order=None,quadrature=None):
        cfg=self.cfg['leading'];order=order or cfg['primary_order'];quadrature=quadrature or cfg['quadrature']
        ident=dict(pair=pair,mu=float(mu),tau=float(tau),xi=float(xi),order=order,quadrature=quadrature)
        key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16];path=self.folder/'points'/(key+'.json')
        if path.exists():
            saved=json.loads(path.read_text());require(saved['identity']==ident and saved['status']=='PASS','Invalid cached point');r=saved['solution'];self.reused+=1
        else:
            model=physics.Leading(*self.cfg['cases'][pair],mu=mu,tau=tau,order=order,quadrature=quadrature)
            started=time.perf_counter();r=model.solve(xi,count=3)
            self.new.append(dict(**ident,seconds=time.perf_counter()-started))
        d=r['diagnostics'];z=np.asarray(r['zeta']);c=np.asarray(r['coefficients'])
        okay=all(d[k]<=cfg['quality'][k+'_tol'] for k in QUALITY)
        okay=okay and len(z)==3 and min(z)>0 and min(np.diff(z))>0 and d['mass_cholesky'] and d['constraint_rank']==d['P_rank']==1
        saved=dict(identity=ident,solution=r,status='PASS' if okay else 'FAIL')
        if not okay:atomic(self.folder/'failed_points'/(key+'.json'),saved)
        require(okay,'Leading physical gate failed: '+str(ident)+' '+str(d))
        if not path.exists():atomic(path,saved)
        for name in ('zeta','omega_lim','Lambda','coefficients','node_value','s_xi'):r[name]=np.asarray(r[name])
        return r

    def regression(self):
        for family,name in [('thickness','leading_curvature.csv'),('length','leading_curvature.csv'),('equal','spectrum.csv'),('equal','curvature.csv')]:self.source(family,name)
        file=self.folder/'regression.json'
        if file.exists():require(json.loads(file.read_text())['status']=='PASS','Unaccepted regression');return
        rows=[];old={(r['pair'],float(r['mu']),float(r['tau'])):r for r in read_csv(self.source('thickness','leading_curvature.csv'))}
        for pair,mu,tau in self.groups:
            r=self.solve(pair,mu,tau,0.);a=old[pair,mu,tau]
            dz=abs(r['zeta'][0]/float(a['zeta0'])-1);dc=abs(r['C_at_zero']/float(a['C'])-1)
            require(dz<1e-10 and dc<1e-10,'Heterogeneous zero regression')
            # Matrices remain those of the frozen physical module; verify formula normalization explicitly.
            m=physics.Leading(*self.cfg['cases'][pair],mu=mu,tau=tau,order=24)
            v=r['coefficients'][:,0];mass=float(v@m.M@v);C=float(v@m.P@v)/(2*m.S_ax*r['zeta'][0]*mass)
            require(abs(C/r['C_at_zero']-1)<1e-12 and abs(mass-1)<1e-12,'C or mass normalization')
            rows.append(dict(family='thickness',pair=pair,mu=mu,tau=tau,xi=0,zeta_relative_error=dz,C_relative_error=dc,status='PASS'))
        old={(r['pair'],float(r['mu'])):r for r in read_csv(self.source('length','leading_curvature.csv'))}
        legacy=read_csv(self.source('equal','spectrum.csv'));curvature=read_csv(self.source('equal','curvature.csv'))
        for pair in self.cfg['cases']:
            for mu in self.cfg['leading']['homogeneous_regression_mu']:
                r=self.solve(pair,mu,0.,0.);a=old[pair,mu]
                dz=abs(r['zeta'][0]/float(a['zeta0'])-1);dc=abs(r['C_at_zero']/float(a['C'])-1)
                require(max(dz,dc)<1e-10,'Homogeneous length regression')
                rows.append(dict(family='length',pair=pair,mu=mu,tau=0,xi=0,zeta_relative_error=dz,C_relative_error=dc,status='PASS'))
            mu=self.cfg['leading']['legacy_regression_mu']
            a=physics.Leading(*self.cfg['cases'][pair],mu=mu,tau=0.,order=24);b=physics.leading.Ritz(*self.cfg['cases'][pair],mu=mu,order=24)
            require(a.S_ax==2 and all(np.array_equal(getattr(a,k),getattr(b,k)) for k in ('K0','M','P','constraint','Z')),'tau=0 matrix limit')
            for xi in self.cfg['leading']['legacy_regression_xi']:
                r=self.solve(pair,mu,0.,xi);source=next(v for v in legacy if v['ID']==pair and float(v['xi'])==xi)
                dz=max(abs(r['zeta'][k-1]/float(source[f'zeta{k}'])-1) for k in (1,2,3));dc=0.
                if xi==0:
                    C=next(float(v['C']) for v in curvature if v['ID']==pair and int(v['order'])==20)
                    dc=abs(r['C_at_zero']/C-1)
                require(max(dz,dc)<1e-10,'Legacy leading regression')
                rows.append(dict(family='equal',pair=pair,mu=mu,tau=0,xi=xi,zeta_relative_error=dz,C_relative_error=dc,status='PASS'))
        self.csv('regression',rows);atomic(file,dict(status='PASS',points=len(rows),max_zeta_relative_error=max(r['zeta_relative_error'] for r in rows),max_C_relative_error=max(r['C_relative_error'] for r in rows)))
        print('Regression PASS:',len(rows),'leading problems',flush=True)

    def calculate(self):
        self.regression();rows=[];neighbours=[];controls=[];local=[]
        cfg=self.cfg['leading'];quality=cfg['quality']
        for pair,mu,tau in self.groups:
            source=sorted([r for r in self.full if (r['pair'],r['mu'],r['tau'])==(pair,mu,tau)],key=lambda r:r['xi'])
            zero=self.solve(pair,mu,tau,0.);previous_w=None;series=[];solutions={}
            for f in source:
                xi=f['xi'];r=self.solve(pair,mu,tau,xi);solutions[xi]=r
                w=r['omega_lim'];jump=0 if previous_w is None else max(abs(w/previous_w-1))
                require(previous_w is None or min(w/previous_w-1)>=-quality['monotonicity_rtol'],'Leading sorted monotonicity')
                previous_w=w
                row=dict(pair=pair,mu=mu,tau=tau,beta_deg=f['beta_deg'],xi=xi,
                    **{f'zeta{k}':float(r['zeta'][k-1]) for k in (1,2,3)},g12_full=f['g12_full'],g12_lim=float(w[1]/w[0]-1),
                    **frequency_comparison(f['Lambda1'],source[0]['Lambda1'],r['zeta'][0],zero['zeta'][0]))
                require(abs(row['F_full']/(f['omega1']/source[0]['omega1'])-1)<2e-14,'Full frequency conversion')
                series.append(row);neighbours.append(dict(pair=pair,mu=mu,tau=tau,xi=xi,jump=float(jump),triggered=jump>cfg['neighbour_jump'],status='PASS'))
            maximum=max(series,key=lambda r:r['delta_percent']);mingap=min(series,key=lambda r:min(r['g12_full'],r['g12_lim']))
            for row in series:
                xi=row['xi'];r=solutions[xi];reasons=[]
                if row['beta_deg'] in cfg['fixed_control_beta_deg']:reasons.append('fixed_control')
                if row is maximum:reasons.append('max_discrepancy')
                if row is mingap:reasons.append('min_gap')
                if min(row['g12_full'],row['g12_lim'])<cfg['small_gap_trigger']:reasons.append('small_gap')
                if any(r['diagnostics'][k]>cfg['residual_control_fraction']*quality[k+'_tol'] for k in QUALITY):reasons.append('residual_margin')
                if next(v['jump'] for v in neighbours if (v['pair'],v['mu'],v['tau'],v['xi'])==(pair,mu,tau,xi))>cfg['neighbour_jump']:reasons.append('neighbour_jump')
                if not reasons:continue
                lower=self.solve(pair,mu,tau,xi,order=cfg['comparison_order'])
                err=max(abs(r['omega_lim']/lower['omega_lim']-1))
                require(err<quality['frequency_convergence_rtol'],'Order convergence')
                qerr=0.;qchecked=('max_discrepancy' in reasons or 'min_gap' in reasons or row['beta_deg'] in (0,max(self.grid)))
                if qchecked:
                    q=self.solve(pair,mu,tau,xi,quadrature=cfg['quadrature_check'])
                    qerr=max(abs(q['omega_lim']/r['omega_lim']-1));require(qerr<quality['frequency_convergence_rtol'],'Quadrature convergence')
                controls.append(dict(pair=pair,mu=mu,tau=tau,beta_deg=row['beta_deg'],xi=xi,reasons=';'.join(reasons),
                    frequency_20_24=float(err),quadrature_checked=qchecked,frequency_64_96=float(qerr),status='PASS'))
            # Reuse actual, already checked local-xi source solutions; no fit over finite angles.
            for path in sorted((self.sources['thickness']/'leading_points').glob('*.json')):
                saved=json.loads(path.read_text());ident=saved['identity']
                if (ident['pair'],ident['mu'],ident['tau'],ident['order'],ident['quadrature'])!=(pair,mu,tau,24,64) or ident['xi'] not in cfg['local_C_xi']:continue
                self.source('thickness',path.relative_to(self.sources['thickness']).as_posix());r=saved['solution']
                require(all(r['diagnostics'][k]<=quality[k+'_tol'] for k in QUALITY),'Source local C quality')
                estimate=math.log(math.sqrt(r['zeta'][0]/zero['zeta'][0]))/ident['xi']**2
                local.append(dict(pair=pair,mu=mu,tau=tau,xi=ident['xi'],C=zero['C_at_zero'],logF_over_xi2=estimate,
                    relative_difference=estimate/zero['C_at_zero']-1,origin='checked_thickness_archive'))
            first=series[1];local.append(dict(pair=pair,mu=mu,tau=tau,xi=first['xi'],C=zero['C_at_zero'],
                logF_over_xi2=math.log(first['F_lim'])/first['xi']**2,relative_difference=math.log(first['F_lim'])/first['xi']**2/zero['C_at_zero']-1,origin='first_actual_comparison_node'))
            rows+=series;self.csv('comparison',rows);print(f'Compared {pair}, mu={mu:g}, tau={tau:g}: {len(series)} nodes',flush=True)
        self.csv('convergence',controls);self.csv('neighbours',neighbours);self.csv('curvature_check',local)
        statistics=[];samples=[]
        for pair,mu,tau in self.groups:
            rr=[r for r in rows if (r['pair'],r['mu'],r['tau'])==(pair,mu,tau)];worst=max(rr,key=lambda r:r['delta_percent'])
            signs=[np.sign(r['e']) for r in rr if r['xi']>0 and abs(r['e'])>1e-10]
            statistics.append(dict(pair=pair,mu=mu,tau=tau,nodes=len(rr),max_delta_percent=worst['delta_percent'],median_delta_percent=float(np.median([r['delta_percent'] for r in rr])),
                max_beta_deg=worst['beta_deg'],max_xi=worst['xi'],e_min=min(r['e'] for r in rr),e_max=max(r['e'] for r in rr),sign_changes=sum(a!=b for a,b in zip(signs,signs[1:])),
                delta_decreases=sum(b['delta_percent']<a['delta_percent']-1e-8 for a,b in zip(rr,rr[1:])),min_gap_full=min(r['g12_full'] for r in rr),min_gap_lim=min(r['g12_lim'] for r in rr)))
            samples += [r for r in rr if r['beta_deg'] in cfg['characteristic_beta_deg']]
        by={(r['pair'],r['mu'],r['tau'],r['xi']):r for r in rows};effects=[];effect_stats=[]
        for pair in self.cfg['cases']:
            for mu in self.cfg['mu_values']:
                rr=[]
                for beta in self.grid:
                    xi=math.radians(beta)/self.cfg['material']['epsilon'];e=sign_effect(by[pair,mu,.1,xi],by[pair,mu,-.1,xi])
                    rr.append(dict(pair=pair,mu=mu,beta_deg=beta,xi=xi,**e))
                resolved=[r for r in rr if r['xi']>0 and min(abs(r['D_full']),abs(r['D_lim']))>self.cfg['interpretation']['D_zero_sign_resolution']]
                worst=max(rr,key=lambda r:abs(r['Delta_D_pp']));end=rr[-1]
                effect_stats.append(dict(pair=pair,mu=mu,nonzero_nodes=len(rr)-1,resolved_nodes=len(resolved),sign_disagreements=sum(np.sign(r['D_full'])!=np.sign(r['D_lim']) for r in resolved),
                    max_abs_Delta_D_pp=abs(worst['Delta_D_pp']),max_beta_deg=worst['beta_deg'],max_xi=worst['xi'],end_D_full_percent=100*end['D_full'],end_D_lim_percent=100*end['D_lim'],end_Delta_D_pp=end['Delta_D_pp']))
                effects+=rr
        self.csv('error_statistics',statistics);self.csv('characteristic_nodes',samples);self.csv('thickness_sign_comparison',effects);self.csv('thickness_sign_statistics',effect_stats)
        all_quality=[]
        for path in sorted((self.folder/'points').glob('*.json')):
            r=json.loads(path.read_text());all_quality.append(dict(**r['identity'],**{k:r['solution']['diagnostics'][k] for k in QUALITY},gap12_relative=r['solution']['diagnostics']['gap12_relative'],status=r['status']))
        self.csv('leading_quality',all_quality)
        atomic(self.folder/'summary.json',dict(status='PASS',full_model_solves=0,full_field_reconstructions=0,main_points=len(rows),
            unique_leading_problems=len(all_quality),order_controls=len(controls),quadrature_controls=sum(r['quadrature_checked'] for r in controls),
            quality_max={k:max(r[k] for r in all_quality) for k in QUALITY},max_frequency_20_24=max(r['frequency_20_24'] for r in controls),
            max_frequency_64_96=max(r['frequency_64_96'] for r in controls),error_statistics=statistics,sign_statistics=effect_stats))

    def finish(self,stage,status):
        source_file=self.folder/'source_hashes.json'
        old=json.loads(source_file.read_text()) if source_file.exists() else {}
        old.update(self.hashes);self.hashes=old
        for name,digest in self.hashes.items():require(api().sha(ROOT/name)==digest,'Source changed during run '+name)
        atomic(self.folder/'source_hashes.json',self.hashes)
        file=self.folder/'executions.json';records=json.loads(file.read_text()) if file.exists() else []
        records.append(dict(stage=stage,status=status,full_model_solves=0,full_field_reconstructions=0,new_leading_problems=len(self.new),new_solutions=self.new,
                            cached_calls=self.reused,seconds=time.perf_counter()-self.started,render_seconds=self.render_seconds))
        atomic(file,records)


def verify_bundle(folder,figures=True):
    m=json.loads((folder/'manifest.json').read_text());require(m['status']=='PASS','Unaccepted bundle')
    for name,digest in m['artifacts'].items():
        if not figures and Path(name).suffix in ('.pdf','.png'):continue
        require(api().sha(folder/name)==digest,'Changed result '+name)
    for name,digest in {**m['sources'],**m['workflow']}.items():require(api().sha(ROOT/name)==digest,'Changed source '+name)
    return m


def finalize(folder):
    failure=folder/'failure.json'
    if failure.exists():
        event=json.loads(failure.read_text())
        if event['status']=='STOPPED':
            require('Unknown symbol' in event['error'] and 'mathrm' in event['error'],'Unresolved failure')
            event.update(status='RESOLVED_RENDER_LABEL',resolution='Mathtext roman subscripts require braces; numerical results unaffected.',physical_gate_failure=False)
            atomic(failure,event)
    paths=[CONFIG,Path(__file__),ROOT/'scripts/analysis/compare_two_springs_thickness_leading.py']
    atomic(folder/'manifest.json',dict(status='PASS',schema='thickness-leading-comparison-results-v1',
        sources=json.loads((folder/'source_hashes.json').read_text()),workflow={p.relative_to(ROOT).as_posix():api().sha(p) for p in paths},
        artifacts={p.relative_to(folder).as_posix():api().sha(p) for p in sorted(folder.rglob('*')) if p.is_file() and p.name!='manifest.json'}))


def render(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,'legend.fontsize':9,
        'pdf.fonttype':42,'lines.linewidth':1.4,'axes.linewidth':.7})
    fmt=FuncFormatter(lambda x,pos:f'{x:g}'.replace('.',','))
    pairs=['J0E4','J4E0'];mus=[-.1,.1];taus=[-.1,.1]
    kappas=[r'$(\kappa_{\mathrm{j}};\kappa_{\mathrm{e}})=(0{,}1;10)$',r'$(\kappa_{\mathrm{j}};\kappa_{\mathrm{e}})=(10;0{,}1)$']
    def label(symbol,x):return '$'+symbol+'='+f'{x:g}'.replace('.','{,}')+'$'
    def arr(rows,key):return np.array([float(r[key]) for r in rows])
    def style(ax):
        ax.grid(True,color='.88',lw=.5);ax.set_axisbelow(True);ax.set_xlim(0,30);ax.set_xticks([0,10,20,30])
        ax.xaxis.set_major_formatter(fmt);ax.yaxis.set_major_formatter(fmt)
    def save(fig,name):
        fig.savefig(folder/(name+'.pdf'),metadata={'CreationDate':None,'ModDate':None})
        fig.savefig(folder/(name+'.png'),dpi=320);plt.close(fig)
    rows=read_csv(folder/'comparison.csv')
    fig,axs=plt.subplots(2,2,figsize=(7.2,6.7),sharex=True,sharey=True)
    for i,mu in enumerate(mus):
        for j,tau in enumerate(taus):
            ax=axs[i,j]
            for k,pair in enumerate(pairs):
                rr=[r for r in rows if r['pair']==pair and float(r['mu'])==mu and float(r['tau'])==tau]
                ax.plot(arr(rr,'xi'),arr(rr,'F_lim'),color=['#155b86','#a65e13'][k],ls=['-','--'][k],label='Ведущая, '+kappas[k])
                ax.plot(arr(rr,'xi'),arr(rr,'F_full'),ls='none',marker=['o','s'][k],markersize=3.2,mfc='white',color=['#155b86','#a65e13'][k],label='Полная, '+kappas[k])
            style(ax);ax.text(.03,.91,'абвг'[2*i+j]+')  '+label(r'\mu',mu)+', '+label(r'\tau',tau),transform=ax.transAxes,fontsize=9.5)
            if j==0:ax.set_ylabel(r'$F(\xi;\mu,\tau)$')
            if i==1:ax.set_xlabel(r'$\xi=\beta_{\mathrm{rad}}/\varepsilon_0$')
    fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',ncol=2,frameon=False,fontsize=8.7)
    fig.subplots_adjust(left=.10,right=.98,bottom=.09,top=.87,hspace=.12,wspace=.08)
    save(fig,FIGURES[0])
    fig,axs=plt.subplots(1,2,figsize=(7.2,4.3),sharey=True)
    for k,pair in enumerate(pairs):
        ax=axs[k]
        for i,mu in enumerate(mus):
            for j,tau in enumerate(taus):
                rr=[r for r in rows if r['pair']==pair and float(r['mu'])==mu and float(r['tau'])==tau]
                ax.plot(arr(rr,'xi'),arr(rr,'delta_percent'),color=['#222222','#155b86'][i],ls=['-','--'][j],marker=['o','s','^','d'][2*i+j],markersize=3.1,mfc='white',markevery=3,label=label(r'\mu',mu)+', '+label(r'\tau',tau))
        style(ax);ax.set_ylim(bottom=0);ax.set_xlabel(r'$\xi$');ax.text(.01,1.04,'аб'[k]+')',transform=ax.transAxes)
        ax.text(.5,1.04,kappas[k],ha='center',transform=ax.transAxes,fontsize=9.5)
    axs[0].set_ylabel(r'$\delta(\xi)$, %')
    fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',ncol=2,frameon=False)
    fig.subplots_adjust(left=.11,right=.98,bottom=.15,top=.77,wspace=.1)
    save(fig,FIGURES[1])
    rows=read_csv(folder/'thickness_sign_comparison.csv')
    fig,axs=plt.subplots(2,2,figsize=(7.2,6.3),sharex=True,sharey=True)
    for i,mu in enumerate(mus):
        for j,pair in enumerate(pairs):
            ax=axs[i,j];rr=[r for r in rows if r['pair']==pair and float(r['mu'])==mu]
            ax.plot(arr(rr,'xi'),100*arr(rr,'D_lim'),color='#222222',label='Ведущая модель')
            ax.plot(arr(rr,'xi'),100*arr(rr,'D_full'),ls='none',marker='o',mfc='white',color='#b46918',markersize=3.7,label='Полная модель')
            style(ax);ax.axhline(0,color='.6',lw=.6)
            ax.text(.04,.90,'абвг'[2*i+j]+')  '+label(r'\mu',mu),transform=ax.transAxes)
            if i==0:ax.text(.5,1.04,kappas[j],ha='center',transform=ax.transAxes,fontsize=9.5)
            if i==1:ax.set_xlabel(r'$\xi$')
            if j==0:ax.set_ylabel(r'$100D_F^{(\tau)}(\xi)$, %')
    fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',ncol=2,frameon=False)
    fig.subplots_adjust(left=.12,right=.98,bottom=.1,top=.88,hspace=.1,wspace=.08)
    save(fig,FIGURES[2])


def write_readme(folder):
    text='''# Проверка ведущей задачи для разных толщин

[Научный отчёт](../../../../docs/theory/isotropic_two_springs_thickness_asymmetry_leading_comparison.md).
Восемь конструкций: μ=±0,1, τ=±0,1; J0E4=(κj;κe)=(0,1;10),
J4E0=(10;0,1). ε0=0,005; b_i/h_i=4; ν=0,3; k_s=5/6; l=E=ρ=1.
Суммарный объём и размерные пружины фиксированы через A0,I0,l.
Внешняя поворотная пружина принадлежит стержню 2; узловые связи точные.

## Источники и нормировки

Полная модель: `thickness_asymmetry/12ff3907cf7cd6da`, таблица `spectrum.csv`
и соответствующие JSON/NPZ. Использованы 232 точки: 29 на конструкцию,
β=0:0,5:3 и 3,25:0,25:8,5 градусов. Исключений по качеству нет.
ξ=β·π/(180ε0), последний узел ξ=29,6705972839. Точки при ξ=30 нет.
Новые полные спектры и восстановление полных форм запрещены программным
ограничением. Поля NPZ проверены по хешам, но не восстанавливались.

При ξ=0 ведущие ζ0 и C воспроизведены по тому же архиву. Однородный предел
проверен по `length_asymmetry/ddf22cb06e8563e2` при μ=±0,1 и по
`leading_validation/f1c670fe62102196` при μ=0,01, ξ=0;10.
Все пути относительно `results/isotropic_two_springs_rl_timoshenko/`.
Прежние физические модули, отчёты, сценарии и архивы проверяются хешами.

F_full=(Lambda1(beta)/Lambda1(0))², F_lim=sqrt(zeta1(xi)/zeta1(0)).
Обе нормировки выполняются внутри одной конструкции.
e=F_full/F_lim−1, delta=100|e|.
D_full/lim=F_full/lim(τ=+0,1)/F_full/lim(τ=−0,1)−1.
Delta_D=D_lim−D_full; Delta_D_pp=100Delta_D — процентные пункты.
Значения при нуле не учитываются как проверка знака эффекта толщин.

## Численные проверки

Неизменённая `physics.Leading`: интегрированный базис Лежандра, основной
порядок 24, квадратура 64. Проверки повышения 20→24 выполнены на 59 узлах,
квадратуры 64→96 — на 22. Выбор включает нуль, промежуточные узлы,
последний узел, максимумы расхождения, минимальные и малые зазоры,
срабатывания контроля невязок и изменения частоты между соседними точками.
Проверены три положительные отсортированные позиции, точные связи,
симметрия, положительность массы, энергия, ортогональность, сильные уравнения
и естественные условия. Допуски не менялись ради согласия моделей.

321 новое решение ведущей задачи: 232 основных, 8 дополнительных однородных
регрессий, 59 проверок порядка и 22 проверки квадратуры. Восемь неоднородных
регрессий входят в основные 232. Локальные проверки C используют 16 ранее
сохранённых решений при ξ=0,01;0,02. Новых решений полной модели — 0.
Политика `frequency-map-v1`, `fast_plot`, `sorted_positions`, K_plot=2,
K_guard=3; рисунки показывают первую позицию. Полнота матричного префикса
Ритца не объявляется отдельной сертификацией всего континуального спектра.

## Таблицы

CSV: UTF-8, разделитель запятая, десятичная точка; углы beta_deg в градусах,
ξ и остальные частотные отношения безразмерны. На рисунках десятичные запятые.

| Файл | Содержание |
| --- | --- |
| `contract.json` | Полный контракт, фактическая сетка и версии библиотек |
| `reused_full_points.csv`, `archive_quality.csv` | Частоты полной модели, пути/хеши JSON/NPZ и сохранённые проверки |
| `angle_nodes.csv` | Точное соответствие β и ξ без интерполяции |
| `points/*.json`, `leading_quality.csv` | Три собственных значения, коэффициенты ведущих форм и все численные невязки |
| `regression.csv/json` | Воспроизведение прежних ζ0, C и однородного предела |
| `comparison.csv` | F_full, F_lim, e, delta_percent, зазоры обеих моделей |
| `error_statistics.csv`, `characteristic_nodes.csv` | Сводка восьми конструкций и значения при β=0;0,5;3;6;8,5° |
| `thickness_sign_comparison.csv`, `thickness_sign_statistics.csv` | Подписанное влияние ±τ, разность в процентных пунктах, проверка знака вне нуля |
| `convergence.csv`, `neighbours.csv` | Порядки/квадратуры, основания адресных проверок и соседние точки |
| `curvature_check.csv` | Сохранённый C и ln(F_lim)/ξ²; малые архивные ξ отделены от первого фактического узла |
| `summary.json`, `executions.json` | Итоговые показатели, количество новых решений и время расчёта/рисования |
| `source_hashes.json`, `manifest.json` | Неизменность источников, кода и результатов |

## Рисунки и предлагаемые подписи

1. `01_normalized_frequency.pdf/png`. **Сопоставление нормированных первых
   частот полной и ведущей моделей.** Четыре поля соответствуют сочетаниям
   μ=±0,1, τ=±0,1. На каждом поле показаны две пары поворотных жёсткостей.
   Линии соединяют решения ведущей задачи, маркеры показывают фактически
   рассчитанные архивные точки полной модели.
2. `02_model_discrepancy.pdf/png`. **Расхождение нормированных частот.**
   δ=100|F_full/F_lim−1|; слева (0,1;10), справа (10;0,1).
   Максимумы относятся к конечной проверенной сетке ξ≤29,670597.
3. `03_thickness_sign_effect.pdf/png`. **Влияние знака асимметрии толщин при
   одинаковой продольной податливости.** Показано 100D_F; строки соответствуют
   μ=−0,1 и +0,1, столбцы — двум парам пружин. Линии — ведущая модель,
   маркеры — полная. Положительный знак означает больший нормированный рост
   при более толстом втором стержне. Численные различия приведены в таблицах.

PDF векторные, PNG 320 dpi. Частоты не интерполируются и не сглаживаются.
Четвёртый рисунок не требуется: разности Delta_D представлены в CSV и отчёте.

## Воспроизведение

Из корня проекта, Python с NumPy, SciPy, Matplotlib:

```powershell
python scripts/analysis/compare_two_springs_thickness_leading.py --stage sources
python scripts/analysis/compare_two_springs_thickness_leading.py --stage regression
python scripts/analysis/compare_two_springs_thickness_leading.py
python scripts/analysis/compare_two_springs_thickness_leading.py --verify-cache BUNDLE
python scripts/analysis/compare_two_springs_thickness_leading.py --plot-only BUNDLE
python -m pytest tests/test_two_springs_thickness_leading_comparison.py -q
```

BUNDLE=`results/isotropic_two_springs_rl_timoshenko/thickness_asymmetry_leading_comparison/KEY`.
Использованный Python: `D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe`.
Повторный основной запуск читает готовые ведущие решения. Режимы чтения
не собирают матриц и не решают собственных задач. Отсутствующие полные точки
не пересчитываются. Отдельный сценарий нужен для контракта «читать полную /
решать ведущую» и сохраняет хеши прежних сценариев завершённых исследований.
Общая асимптотическая сходимость, другие ε0, высокие моды, характерные углы
и другая модель опоры в данной работе не исследуются.
'''
    (folder/'README.md').write_text(text.replace('KEY',folder.name),encoding='utf8')


def main(args):
    if args.verify_cache or args.plot_only:
        require(not(args.verify_cache and args.plot_only),'Choose one saved-data mode')
        folder=(args.verify_cache or args.plot_only).resolve();require(folder.parent==BASE,'Wrong result family')
        verify_bundle(folder,not bool(args.plot_only))
        if args.plot_only:render(folder);verify_bundle(folder)
        print('Saved data verified; 0 full and 0 leading solves');return
    exp=None
    with forbid_full_calculations():
        try:
            exp=Experiment();print(exp.folder,flush=True)
            if args.stage in ('regression','all'):exp.regression()
            if args.stage=='all':
                exp.calculate();write_readme(exp.folder)
                started=time.perf_counter();render(exp.folder);exp.render_seconds=time.perf_counter()-started
            exp.finish(args.stage,'PASS')
            if args.stage=='all':finalize(exp.folder)
        except Exception as error:
            if exp is not None:
                atomic(exp.folder/'failure.json',dict(status='STOPPED',error=repr(error)));exp.finish(args.stage,'STOPPED')
            raise
