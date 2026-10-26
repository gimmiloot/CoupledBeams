"""Bounded free-rotation study; shared full solver, validators and archive I/O."""
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
import scipy
from scripts.lib import two_springs_external_transverse_spring as previous
from scripts.lib import two_springs_external_transverse_leading_resume as resumed
from scripts.lib import isotropic_two_springs_external_transverse_spring as physics
from scripts.lib.isotropic_two_springs_external_rotational_replacement import leading_model

ROOT=previous.ROOT
BASE=previous.BASE.parent/'external_rotational_spring_replacement'
CONFIG=ROOT/'data/input/isotropic_two_springs_external_rotational_replacement.json'
api=previous.api
require=previous.require
read_csv=previous.read_csv
stiffness=previous.stiffness
FIGURES=('01_first_frequency','02_normalized_growth','03_leading_comparison','04_initial_growth_coefficient')


class Experiment(previous.Experiment):
    def __init__(self):
        self.a=api();self.cfg=json.loads(CONFIG.read_text());self.inputs={};self.manifests={};self.sources={}
        self.started=time.perf_counter();self.leading_new=0;self.ritz_new=0
        for name,digest in self.cfg['frozen_files'].items():self.check(ROOT/name,digest)
        for name,s in self.cfg['sources'].items():
            folder=ROOT/s['folder'];self.check(folder/'manifest.json',s['manifest_sha256'])
            resumed.diagnostic.verify_source(folder)
            self.sources[name]=folder;self.manifests[name]=json.loads((folder/'manifest.json').read_text())
        require(self.manifests['full']['full_model_status']=='PASS' and self.manifests['leading']['status']=='PASS','Source acceptance')
        old=json.loads(self.source('full','contract.json').read_text())['config']
        require(old['material']==self.cfg['material'] and old['mu']==self.cfg['mu'] and old['variant']=='RL_FULL','Source physics')
        require(old['map']==self.cfg['map'],'Grid or semantics differs')
        for pair,ref in self.cfg['reference_pairs'].items():
            require(old['cases'][ref]==[self.cfg['cases'][pair][0],self.cfg['reference_k_end'][pair]],'Reference rotational stiffness')
        identity=dict(config=self.cfg,specialization_sha256=self.a.sha(ROOT/'scripts/lib/isotropic_two_springs_external_rotational_replacement.py'),versions=dict(numpy=np.__version__,scipy=scipy.__version__))
        self.key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16]
        self.folder=BASE/self.key;self.folder.mkdir(parents=True,exist_ok=True)
        target=self.folder/'contract.json'
        if target.exists():require(json.loads(target.read_text())==identity,'Immutable contract differs')
        else:self.a.write(target,identity)
        self.run=self.a.Run(self.cfg,self.folder,engine=physics)
        target=self.folder/'pole_catalogs.json'
        if target.exists():self.run.catalogs={float(k):v for k,v in json.loads(target.read_text()).items()}
        self.old={}
        for r in read_csv(self.source('full','spectrum.csv')):
            pair=next(k for k,v in self.cfg['reference_pairs'].items() if v==r['pair']);kw=stiffness(r['k_w']);b=float(r['beta_deg'])
            c=physics.Case(b,self.cfg['mu'],self.cfg['cases'][pair][0],self.cfg['reference_k_end'][pair],kw)
            if kw=='RIGID':c=physics.base.Case(b,c.mu,c.k_joint,c.k_end)
            self.archived(r,c)
            self.old[pair,kw,b]={**r,'pair':pair,'k_w':kw,'beta_deg':b,**{n:float(r[n]) for n in r if n.startswith(('omega','Lambda','F_','gap'))}}
        require(len(self.old)==510,'Incomplete previous map')

    def point(self,pair,k,beta,previous_w=None,fresh=False):
        require(pair in self.cfg['cases'] and k in self.cfg['k_w_values'] and beta in self.cfg['map']['parameter_grid'],'Unrequested point')
        case=self.case(pair,k,beta);rec=self.run.solve(case,neighbour=previous_w);self.validate(rec,case)
        self.a.write(self.folder/'pole_catalogs.json',self.run.catalogs)
        return rec,previous.previous.lengths.point_path(self.folder,case),'new_study'

    def preflight(self):
        target=self.folder/'preflight.json'
        if target.exists():require(json.loads(target.read_text())['status']=='PASS','Failed preflight');return
        ref=self.cfg['hinge_reference'];path=self.check(ROOT/ref['point_file'],ref['point_sha256']);old=json.loads(path.read_text())
        self.check(path.with_suffix('.npz'),old['shape_sha256'])
        c=physics.Case(30.,self.cfg['mu'],1.,0.,'RIGID');rec=self.run.solve(c);self.validate(rec,c)
        error=max(abs(a['omega']/b['omega']-1) for a,b in zip(rec['roots'],old['roots']))
        require(error<=self.cfg['policy']['same_physics_frequency_tol'],'Canonical hinge regression')
        self.a.write(self.folder/'hinge_regression.json',dict(status='PASS',relative_frequency_error=error,source=ref,case=rec['case'],Lambda=[r['Lambda'] for r in rec['roots']]))
        points={};rows=[];ritz=[];m=self.run.model
        for pair in self.cfg['cases']:
            for b in self.cfg['preflight_angles']:
                for k in self.cfg['k_w_values']:
                    c=self.case(pair,k,b);ext=physics.external_matrix(m,c,1)
                    require(np.array_equal(ext[0],np.eye(6)[0]) and np.array_equal(ext[2],[0,0,0,0,0,-1]),'Outer axial or moment condition')
                    if k=='RIGID':require(np.array_equal(ext[1],np.eye(6)[1]),'Rigid transverse condition')
                    else:require(np.array_equal(ext[1],[0,k*m.section.EI/m.l**3,0,0,-1,0]),'Outer transverse sign')
                    require(np.array_equal(physics.joint_matrix(m,c),physics.base.joint_matrix(m,c)),'Joint changed')
                    rec,_,_=self.point(pair,k,b);points[pair,k,b]=rec
                    require(rec['count']['samples'][0]['count']==0,'Incomplete lower spectrum')
                    rows.append(dict(pair=pair,k_w=k,beta_deg=b,lower_Lambda=self.cfg['policy']['Lambda_lower'],lower_count=0,upper_count=rec['count']['samples'][1]['count'],Lambda1=rec['roots'][0]['Lambda'],Lambda3=rec['roots'][2]['Lambda'],status='PASS'))
            print('Full preflight: '+pair,flush=True)
            for k,b in self.cfg['ritz_cases']:
                expected=np.array([r['omega'] for r in points[pair,k,b]['roots']])
                for p in self.cfg['policy']['ritz_orders']:
                    w,coef,d=physics.ritz(m,self.case(pair,k,b),p);self.ritz_new+=1;err=float(max(abs(w/expected-1)))
                    ritz.append(dict(pair=pair,k_w=k,beta_deg=b,order=p,relative_frequency_error=err,**d));self.csv('ritz_checks',ritz)
                    require(err<=self.cfg['policy']['ritz_frequency_tol'] and d['mass_gram_error']<=self.cfg['policy']['mass_gram_tol'] and d['energy_error']<=self.cfg['policy']['energy_tol'] and d['constraint_residual']<=self.cfg['policy']['boundary_tol'],'Independent full Ritz gate')
        self.csv('preflight_points',rows);self.csv('preflight_monotonicity',self.monotonic(points,self.cfg['preflight_angles']))
        self.a.write(target,dict(status='PASS',new_grid_points=len(rows),hinge_regression_error=error,independent_ritz_solves=len(ritz),max_ritz_error=max(r['relative_frequency_error'] for r in ritz),zero_mechanism=False,min_Lambda1=min(r['Lambda1'] for r in rows),lower_count_all_zero=True))

    def comparisons(self,rows):
        comparisons=[];interactions=[];by={(r['pair'],r['k_w'],r['beta_deg']):r for r in rows}
        for r in rows:
            pair,k,b=r['pair'],r['k_w'],r['beta_deg'];B=self.old[pair,k,b];A=self.old[pair,'RIGID',b];C=by[pair,'RIGID',b]
            ratios=[r[f'omega{i}']/B[f'omega{i}']-1 for i in (1,2,3)]
            require(max(ratios)<=self.cfg['policy']['same_physics_frequency_tol'],'Removal of rotation spring increased an absolute frequency')
            comparisons.append(dict(pair=pair,k_w=k,beta_deg=b,**{f'Lambda_{n}':v['Lambda1'] for n,v in zip('ABCD',[A,B,C,r])},**{f'F_{n}':v['F_omega'] for n,v in zip('ABCD',[A,B,C,r])},D_Lambda_CA=C['Lambda1']/A['Lambda1']-1,D_Lambda_DB=r['Lambda1']/B['Lambda1']-1,D_Lambda_DC=r['Lambda1']/C['Lambda1']-1,D_F_CA=C['F_omega']/A['F_omega']-1,D_F_DB=r['F_omega']/B['F_omega']-1,D_F_DC=r['F_omega']/C['F_omega']-1,**{f'omega{i}_removal_relative_change':v for i,v in enumerate(ratios,1)}))
            if b in (0,3,10,20) and k!='RIGID':interactions.append(dict(pair=pair,k_w=k,beta_deg=b,I_zeta=r['Lambda1']**4-C['Lambda1']**4-B['Lambda1']**4+A['Lambda1']**4))
        self.csv('four_supports',comparisons);self.csv('interaction',interactions)
        summary,mechanics=self.analyze(rows)
        for r in summary:
            pair,k=r['pair'],r['k_w'];old=self.old
            r.update(gain_10_20_percent=100*(r['F20']/r['F10']-1),old_gain_0_3_percent=100*(old[pair,k,3]['F_omega']-1),Lambda0_relative_removal=r['Lambda0']/old[pair,k,0]['Lambda1']-1)
        self.csv('summary',summary)
        return summary

    def finish(self,stage,status='PASS'):
        super().finish(stage,status)
        path=self.folder/'executions.json';runs=json.loads(path.read_text())
        runs[-1]['workflow_sha256']=self.a.sha(__file__)
        self.a.write(path,runs)


class LeadingExperiment(resumed.Experiment):
    """Reuse accepted leading point transactions/checks with the new K0."""
    def __init__(self,full):
        self.parent=full;self.folder=full.folder;self.a=full.a;self.new=0
        self.cfg={**full.cfg['leading'],'quality':full.cfg['leading']['quality'],'quadrature':full.cfg['leading']['quadrature']}

    def leading(self,pair,k,xi,p=24):
        return leading_model(self.parent.cfg['cases'][pair][0],k,mu=self.parent.cfg['mu'],order=p,quadrature=self.cfg['quadrature'])

    def preflight(self):
        path=self.folder/'leading_preflight.json'
        if path.exists():require(json.loads(path.read_text())['status']=='PASS','Leading preflight failed');return
        rows=[]
        for pair in self.parent.cfg['cases']:
            for k in self.parent.cfg['k_w_values']:
                xs=self.cfg['preflight_xi'] if k in self.cfg['preflight_k_w'] else [0.]
                for xi in xs:rows.append(self.convergence(pair,k,xi,'preflight'))
            print('Leading preflight: '+pair,flush=True)
        self.csv('leading_preflight_convergence',rows)
        self.a.write(path,dict(status='PASS',compared_points=len(rows),max_frequency_error=max(r['relative_frequency_error'] for r in rows)))

    def run(self,full):
        self.preflight();cfg=self.parent.cfg;zeros={};curvature=[];local=[];rows=[]
        old={(r['pair'],stiffness(r['k_w'])):r for r in read_csv(self.parent.source('leading','leading_curvature.csv'))}
        convergence=read_csv(self.folder/'leading_preflight_convergence.csv')
        for pair in cfg['cases']:
            for k in cfg['k_w_values']:
                z=self.point(pair,k,0.);zeros[pair,k]=z
                r=old[cfg['reference_pairs'][pair],k];oldC=float(r['C']);oldz=float(r['zeta0']);oldy=float(r['node0_squared'])
                C=z['C_at_zero'];y2=z['node_value'][0]**2
                curvature.append(dict(pair=pair,k_w=k,zeta0=z['zeta'][0],Lambda0=z['Lambda'][0],node0_squared=y2,C=C,old_C=oldC,old_zeta0=oldz,old_node0_squared=oldy,C_new_to_old=C/oldC,node_square_new_to_old=y2/oldy,zeta_new_to_old=z['zeta'][0]/oldz,gap12=z['diagnostics']['gap12_relative']))
                for xi in self.cfg['local_C_xi']:
                    s=self.point(pair,k,xi);estimate=math.log(math.sqrt(s['zeta'][0]/z['zeta'][0]))/xi**2;err=abs(estimate/C-1)
                    require(err<=self.cfg['local_C_rtol'],'Local C check')
                    local.append(dict(pair=pair,k_w=k,xi=xi,C=C,logF_over_xi2=estimate,relative_error=err))
        self.csv('leading_curvature',curvature);self.csv('local_C_checks',local)
        for r in full:
            pair,k,b=r['pair'],r['k_w'],r['beta_deg'];xi=math.radians(b)/cfg['material']['epsilon']
            if xi>self.cfg['xi_max']:continue
            s=self.point(pair,k,xi);z=zeros[pair,k];F=math.sqrt(s['zeta'][0]/z['zeta'][0]);e=r['F_omega']/F-1
            rows.append(dict(pair=pair,k_w=k,beta_deg=b,xi=xi,**{f'zeta{i}':v for i,v in enumerate(s['zeta'],1)},F_full=r['F_omega'],F_lim=F,e=e,delta_percent=100*abs(e),gap12_full=r['gap12_relative'],gap12_lim=s['diagnostics']['gap12_relative'],point_file=r['point_file'],point_sha256=r['point_sha256']))
            if b in self.cfg['check_beta_deg'] and b!=0:convergence.append(self.convergence(pair,k,xi,'representative_archive_node'))
        stats=[];growth=[]
        for pair in cfg['cases']:
            for k in cfg['k_w_values']:
                rr=[r for r in rows if (r['pair'],r['k_w'])==(pair,k)];worst=max(rr,key=lambda r:r['delta_percent']);last=rr[-1]
                convergence.append(self.convergence(pair,k,worst['xi'],'maximum_model_difference'))
                stats.append(dict(pair=pair,k_w=k,nodes=len(rr),max_delta_percent=worst['delta_percent'],median_delta_percent=float(np.median([r['delta_percent'] for r in rr])),worst_beta_deg=worst['beta_deg'],worst_xi=worst['xi'],e_last=last['e'],min_gap12_full=min(r['gap12_full'] for r in rr),min_gap12_lim=min(r['gap12_lim'] for r in rr)))
                r=next(r for r in rr if r['beta_deg']==3.)
                growth.append(dict(pair=pair,k_w=k,gain_full_percent=100*(r['F_full']-1),gain_lim_percent=100*(r['F_lim']-1)))
        by={(r['pair'],r['k_w'],r['xi']):r for r in rows};grid=sorted({r['xi'] for r in rows});monotonic=[]
        for pair in cfg['cases']:
            for xi in grid:
                for lo,hi in zip([0.,.1,1.,10.],[.1,1.,10.,'RIGID']):
                    diff=[math.sqrt(by[pair,hi,xi][f'zeta{i}']/by[pair,lo,xi][f'zeta{i}'])-1 for i in (1,2,3)]
                    require(min(diff)>=-self.cfg['quality']['monotonicity_rtol'],'Leading kw monotonicity')
                    monotonic.append(dict(pair=pair,parameter='k_w',fixed=xi,lower=lo,upper=hi,min_increment=min(diff)))
            for k in cfg['k_w_values']:
                for lo,hi in zip(grid,grid[1:]):
                    diff=[math.sqrt(by[pair,k,hi][f'zeta{i}']/by[pair,k,lo][f'zeta{i}'])-1 for i in (1,2,3)]
                    require(min(diff)>=-self.cfg['quality']['monotonicity_rtol'],'Leading xi monotonicity')
                    monotonic.append(dict(pair=pair,parameter='xi',fixed=k,lower=lo,upper=hi,min_increment=min(diff)))
        for name,data in [('leading_comparison',rows),('leading_summary',stats),('leading_convergence',convergence),('growth_at_3deg',growth),('leading_monotonicity',monotonic)]:self.csv(name,data)
        points=[json.loads(p.read_text()) for p in (self.folder/'leading_points').glob('*.json')]
        quality=[dict(**r['identity'],**{n:r['solution']['diagnostics'][n] for n in resumed.QUALITY},status=r['status']) for r in points]
        require(all(r['status']=='PASS' for r in quality),'Unaccepted leading point')
        self.csv('leading_quality',quality)
        self.a.write(self.folder/'leading_numerical_summary.json',dict(status='PASS',accepted_solves=len(points),main_solutions=len(rows),quality_maxima={n:max(r[n] for r in quality) for n in resumed.QUALITY},max_order_frequency_error=max(float(r['relative_frequency_error']) for r in convergence),max_C_order_error=max(float(r['relative_C_error']) for r in convergence if r['relative_C_error'] not in ('',None)),max_local_C_error=max(r['relative_error'] for r in local),monotonicity_checks=len(monotonic)))


def render(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    cfg=json.loads((folder/'contract.json').read_text())['config']
    full=read_csv(folder/'spectrum.csv');four=read_csv(folder/'four_supports.csv');leading=read_csv(folder/'leading_comparison.csv');C=read_csv(folder/'leading_curvature.csv')
    plt.rcParams.update({'font.size':10,'axes.labelsize':11,'legend.fontsize':8,'pdf.fonttype':42,'font.family':'DejaVu Sans'})
    colors=['#222222','#0072B2','#D55E00','#009E73','#8B4A9C'];styles=['-','--','-.',':',(0,(5,1,1,1))];markers=['o','s','^','v','D']
    comma=FuncFormatter(lambda x,pos:f'{x:g}'.replace('.',','))
    labels=['Точная опора',r'$\kappa_w=10$',r'$\kappa_w=1$',r'$\kappa_w=0{,}1$',r'$\kappa_w=0$']
    def group(rows,pair,k):return [r for r in rows if r['pair']==pair and stiffness(r['k_w'])==k]
    def values(rows,key):return [float(r[key]) for r in rows]
    def annotation(ax,pair):
        ax.text(.04,.95,r'$\kappa_{\mathrm{j}}='+f'{cfg["cases"][pair][0]:g}'.replace('.',r'{,}')+'$',transform=ax.transAxes,va='top')
    def save(fig,axes,index):
        for ax in np.asarray(axes).flat:
            ax.grid(alpha=.2);ax.xaxis.set_major_formatter(comma);ax.yaxis.set_major_formatter(comma)
        fig.savefig(folder/(FIGURES[index]+'.pdf'));fig.savefig(folder/(FIGURES[index]+'.png'),dpi=320);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(7.6,3.6),sharey=True,layout='constrained')
    for ax,pair in zip(axes,cfg['cases']):
        for j,k in enumerate(cfg['k_w_values']):
            rr=group(full,pair,k);ax.plot(values(rr,'beta_deg'),values(rr,'Lambda1'),color=colors[j],ls=styles[j],lw=1.4,label=labels[j])
        annotation(ax,pair);ax.set_xlabel(r'$\beta$, град.');ax.set_xlim(0,20)
    axes[0].set_ylabel(r'$\Lambda_1$');fig.legend(*axes[0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False)
    save(fig,axes,0)
    fig,axes=plt.subplots(2,3,figsize=(8.0,5.7),sharex=True,sharey='row',layout='constrained')
    for i,pair in enumerate(cfg['cases']):
        for j,k in enumerate(['RIGID',1.,0.]):
            ax=axes[i,j];rr=group(four,pair,k);x=values(rr,'beta_deg')
            ax.plot(x,values(rr,'F_D'),color='#0072B2',lw=1.6,label='Свободный поворот')
            ax.plot(x,values(rr,'F_B'),color='#333333',ls='--',lw=1.4,label='Поворотная пружина')
            annotation(ax,pair)
            xy=(.46,.16) if k=='RIGID' else (.04,.81)
            ax.text(*xy,labels[[0,2,4][j]],transform=ax.transAxes,va='top',fontsize=9)
            ax.set_xlim(0,20)
            if i==1:ax.set_xlabel(r'$\beta$, град.')
        axes[i,0].set_ylabel(r'$F=\omega_1(\beta)/\omega_1(0)$')
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=2,frameon=False);save(fig,axes,1)
    fig,axes=plt.subplots(2,2,figsize=(7.6,6.0),sharex=True,sharey='row',layout='constrained')
    for j,pair in enumerate(cfg['cases']):
        for i,k in enumerate(cfg['k_w_values']):
            rr=group(leading,pair,k);x=values(rr,'xi')
            axes[0,j].plot(x,values(rr,'F_lim'),color=colors[i],ls=styles[i],label=labels[i],lw=1.3)
            axes[0,j].plot(x,values(rr,'F_full'),ls='none',marker=markers[i],mfc='none',mec=colors[i],ms=3.0,mew=.65)
            axes[1,j].plot(x,values(rr,'delta_percent'),color=colors[i],ls=styles[i],lw=1.3)
        annotation(axes[0,j],pair);axes[1,j].set_xlabel(r'$\xi=\beta_{\mathrm{rad}}/\varepsilon$')
        axes[0,j].set_xlim(0,30)
    axes[0,0].set_ylabel(r'$F_{\mathrm{full}},\ F_{\mathrm{lim}}$');axes[1,0].set_ylabel(r'$\delta$, %')
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=3,frameon=False);save(fig,axes,2)
    fig,axes=plt.subplots(2,3,figsize=(7.6,5.5),sharey='row',gridspec_kw={'width_ratios':[1,3,1]},layout='constrained')
    for i,pair in enumerate(cfg['cases']):
        rr={stiffness(r['k_w']):r for r in C if r['pair']==pair}
        for key,color,style,marker,label in [('C','#0072B2','-','o','Свободный поворот'),('old_C','#333333','--','s','Поворотная пружина')]:
            axes[i,1].semilogx([.1,1,10],[float(rr[k][key]) for k in [.1,1,10]],color=color,ls=style,marker=marker,label=label)
            for col,k in [(0,0.),(2,'RIGID')]:axes[i,col].plot([0],[float(rr[k][key])],color=color,marker=marker)
        axes[i,0].set_ylabel(r'$C$');annotation(axes[i,1],pair);axes[i,1].set_xlabel(r'$\kappa_w$')
    # Save separately so categorical endpoints are never formatted as numbers.
    for i in range(2):
        axes[i,0].set_xticks([0],['0']);axes[i,1].set_xticks([.1,1,10],['0,1','1','10']);axes[i,2].set_xticks([0],['Точная\nопора'])
    for ax in axes.flat:ax.grid(alpha=.2);ax.yaxis.set_major_formatter(comma)
    fig.legend(*axes[0,1].get_legend_handles_labels(),loc='outside lower center',ncol=2,frameon=False)
    fig.savefig(folder/(FIGURES[3]+'.pdf'));fig.savefig(folder/(FIGURES[3]+'.png'),dpi=320);plt.close(fig)


def finalize(folder):
    a=api();cfg=json.loads((folder/'contract.json').read_text())['config']
    for n in ('preflight','leading_preflight','leading_numerical_summary'):require(json.loads((folder/(n+'.json')).read_text())['status']=='PASS','Unaccepted stage '+n)
    require(not (folder/'failure.json').exists(),'Unresolved failure')
    for name in FIGURES:
        for suffix in ('.pdf','.png'):require((folder/(name+suffix)).is_file(),'Missing figure')
    hashes=json.loads((folder/'source_hashes.json').read_text())
    for name,digest in hashes.items():require(a.sha(ROOT/name)==digest,'Protected source changed '+name)
    integrity={name:resumed.diagnostic.verify_source(ROOT/r['folder']) for name,r in cfg['sources'].items()}
    a.write(folder/'source_integrity.json',integrity)
    q=read_csv(folder/'quality.csv');runs=json.loads((folder/'executions.json').read_text());lead=json.loads((folder/'leading_numerical_summary.json').read_text())
    decision=dict(status='PASS',main_full_points=len(read_csv(folder/'spectrum.csv')),main_leading_points=lead['main_solutions'],new_full_prefixes=sum(r['new_full_prefixes'] for r in runs),new_full_roots=sum(r['new_full_roots'] for r in runs),independent_full_ritz_solves=sum(r['new_ritz_solves'] for r in runs),new_leading_solutions=sum(r['new_leading_solves'] for r in runs),full_quality_maxima={n:max(float(r[n]) for r in q) for n in ('physical','ode','energy','singular','mass_gram','condition')},leading=lead)
    require(decision['new_leading_solutions']==lead['accepted_solves'],'Leading solve accounting')
    a.write(folder/'decision.json',decision)
    a.write(folder/'manifest.json',dict(status='PASS',artifacts={p.relative_to(folder).as_posix():a.sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))
    return resumed.diagnostic.verify_source(folder)


def main(args):
    if args.stage in ('plot_only','finalize'):
        require(args.folder is not None and args.folder.resolve().is_relative_to(BASE.resolve()),'Explicit new study folder required')
        with resumed.forbid_full_calculations(leading_too=True):
            if args.stage=='plot_only':
                started=time.perf_counter();render(args.folder)
                api().write(args.folder/'rendering.json',dict(seconds=time.perf_counter()-started,new_spectral_solves=0,dpi=320))
            else:print(finalize(args.folder))
        return
    require(args.stage in ('preflight','run'),'Replacement study supports preflight/run/plot_only/finalize')
    exp=Experiment();print(exp.folder,flush=True)
    try:
        require(not (exp.folder/'failure.json').exists(),'Saved failure requires addressed diagnosis')
        exp.preflight()
        if args.stage=='run':
            rows,quality=exp.sweep();exp.comparisons(rows)
            leading=LeadingExperiment(exp)
            try:leading.run(rows)
            finally:exp.leading_new+=leading.new
        exp.finish(args.stage)
    except Exception as error:
        exp.a.write(exp.folder/'failure.json',dict(status='PARTIAL',stage=args.stage,error=repr(error)))
        exp.finish(args.stage,'PARTIAL');raise
