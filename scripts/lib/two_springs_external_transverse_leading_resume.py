"""Leading-only continuation of the stopped transverse-support study.

The full model is read-only. Reuses canonical Leading, archive validators,
quartic frequency conversion and writers; no new physical assembly or solver.
"""
from contextlib import contextmanager
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import numpy as np
import scipy
from scripts.lib import two_springs_external_transverse_spring as previous
from scripts.lib import isotropic_two_springs_external_transverse_spring as physics
from scripts.lib import two_springs_transverse_leading_diagnostic as diagnostic
from scripts.lib.two_springs_thickness_leading_comparison import frequency_comparison
from scripts.lib.isotropic_two_springs_leading import fields

ROOT=previous.ROOT
BASE=previous.previous.lengths.BASE/'external_transverse_leading_resume'
CONFIG=ROOT/'data/input/isotropic_two_springs_external_transverse_leading_resume.json'
api=previous.api
require=previous.require
read_csv=previous.read_csv
stiffness=previous.stiffness
QUALITY=previous.QUALITY
FIGURES=('03_leading_comparison','04_initial_growth_coefficient')

@contextmanager
def forbid_full_calculations(leading_too=False):
    def forbidden(*args,**kwargs):raise RuntimeError('Spectral calculation forbidden in this operation')
    targets=[(m,n) for m in (physics,physics.base) for n in ('prefix','mode','ritz')]+[(api().Run,'solve')]
    if leading_too:targets.append((physics.Leading,'solve'))
    saved=[getattr(m,n) for m,n in targets]
    try:
        for m,n in targets:setattr(m,n,forbidden)
        yield
    finally:
        for (m,n),f in zip(targets,saved):setattr(m,n,f)

def detailed_checks(model,solution,xi):
    """Per-position strong/trace evidence and finite Ritz prefix inertia."""
    coefficients=np.asarray(solution['coefficients']);zeta=np.asarray(solution['zeta']);rows=[]
    finite=model.k_w!='RIGID';p=model.order;t=xi*xi/2
    for j,z in enumerate(zeta):
        c=coefficients[:,j];ends=np.array([[v@c for v in a] for a in model.ends]);outer=np.array([[v@c for v in a] for a in model.outer])
        residual=dict(end_rotation=-outer[1,2]+model.k_end*outer[1,1],joint_moment_1=ends[0,2]+model.k_joint*(model.jump@c),
                      joint_moment_2=ends[1,2]-model.k_joint*(model.jump@c),joint_force=ends[0,3]-ends[1,3]+t*(model.y@c))
        scale=max(1.,np.max(abs(ends)),np.max(abs(outer)),abs(t*(model.y@c)))
        if finite:
            residual['outer_transverse_force']=outer[1,3]+model.k_w*outer[1,0]
            scale=max(scale,abs(model.k_w*outer[1,0]))
        arms=[]
        for i,L in enumerate(model.lengths):
            f=fields(p,model.nodes,L,i);v=c[i*p:(i+1)*p]
            W=f[0]@v+(c[-1] if finite and i==1 else 0.);D4=f[4]@v
            num=float(L*(model.weights@((D4-z*W)**2)));den=float(L*(model.weights@(D4**2+(z*W)**2)))
            arms.append(dict(arm=i+1,numerator_squared=num,denominator_squared=den))
        rows.append(dict(sorted_position=j+1,strong=math.sqrt(sum(a['numerator_squared'] for a in arms)/sum(a['denominator_squared'] for a in arms)),
                         natural=max(map(abs,residual.values()))/scale,natural_raw=residual,natural_scale=scale,arms=arms))
    K=model.K0+t*model.P;Z=model.Z;Kr=Z.T@K@Z;Mr=Z.T@model.M@Z;scale=np.sqrt(np.diag(Kr));counts=[]
    for shift,expected in zip([zeta[0]/2,(zeta[0]+zeta[1])/2,(zeta[1]+zeta[2])/2,zeta[2]*(1+1e-5)],[0,1,2,3]):
        B=(Kr-shift*Mr)/scale[:,None]/scale[None,:];eig=np.linalg.eigvalsh((B+B.T)/2)
        margin=float(min(abs(eig))/max(abs(eig)));count=int(sum(eig<0))
        require(count==expected and margin>100*np.finfo(float).eps,'Ritz prefix inertia ambiguous or incomplete')
        counts.append(dict(zeta_shift=float(shift),negative_count=count,relative_sign_margin=margin))
    return dict(positions=rows,inertia=counts,scope='Complete first three positions of the stated finite Ritz space; physical checks and order comparison remain separate.')

class Experiment:
    def __init__(self):
        self.a=api();self.cfg=json.loads(CONFIG.read_text());self.sources={};self.hashes={};self.new=0;self.started=time.perf_counter()
        for path,digest in self.cfg['frozen_files'].items():self.check(ROOT/path,digest)
        for name,record in self.cfg['sources'].items():
            folder=ROOT/record['folder'];self.check(folder/'manifest.json',record['manifest_sha256']);self.sources[name]=folder
            diagnostic.verify_source(folder)
        full=self.sources['full'];manifest=json.loads((full/'manifest.json').read_text());diag=json.loads((self.sources['diagnostic']/'decision.json').read_text())
        require(manifest['full_model_status']=='PASS' and manifest['leading_status']=='STOP' and diag['status']=='DIAGNOSTIC_PASS','Wrong source status')
        source_cfg=json.loads((full/'contract.json').read_text())['config'];self.source_cfg=source_cfg
        require(source_cfg['material']==self.cfg['material'] and source_cfg['cases']==self.cfg['cases'] and source_cfg['mu']==self.cfg['mu'],'Physical contract mismatch')
        require(source_cfg['variant']=='RL_FULL' and source_cfg['map']['spectrum_semantics']=='sorted_positions','Full spectrum semantics')
        require(self.cfg['quality']==source_cfg['leading']['quality'] and self.cfg['analytic_curvature'],'Tolerance or assembly changed')
        self.model=physics.Model(**self.cfg['material']);validator=SimpleNamespace(run=SimpleNamespace(model=self.model,policy=source_cfg['policy']))
        self.full=[];self.archive_checks=[];self.all_full=read_csv(full/'spectrum.csv')
        for r in self.all_full:
            pair,k,b=r['pair'],stiffness(r['k_w']),float(r['beta_deg']);xi=math.radians(b)/self.cfg['material']['epsilon']
            if xi>self.cfg['xi_max']:continue
            require(pair in self.cfg['cases'] and k in self.cfg['k_w_values'],'Unrequested source parameters')
            path=ROOT/r['point_file'];self.check(path,r['point_sha256']);self.check(path.with_suffix('.npz'),r['shape_sha256']);record=json.loads(path.read_text())
            case=physics.base.Case(b,self.cfg['mu'],*self.cfg['cases'][pair]) if k=='RIGID' else physics.Case(b,self.cfg['mu'],*self.cfg['cases'][pair],k)
            previous.previous.Experiment.validate(validator,record,case)
            require(all(abs(float(r[f'{name}{j}'])/v[name]-1)<1e-14 for j,v in enumerate(record['roots'],1) for name in ('omega','Lambda')),'CSV/root mismatch')
            self.full.append(dict(pair=pair,k_w=k,beta_deg=b,xi=xi,**{f'{name}{j}':float(r[f'{name}{j}']) for name in ('omega','Lambda') for j in (1,2,3)},
                                 F_omega=float(r['F_omega']),gap12_full=float(r['gap12_relative']),point_file=r['point_file'],point_sha256=r['point_sha256'],shape_sha256=r['shape_sha256']))
            self.archive_checks.append(dict(pair=pair,k_w=k,beta_deg=b,status='PASS',physical=max(max(np.max(v['diagnostics']['joint_scaled']),np.max(v['diagnostics']['external_scaled'])) for v in record['roots']),
                strong=max(v['diagnostics']['ode_error'] for v in record['roots']),energy=max(v['diagnostics']['energy_error'] for v in record['roots']),mass_gram=record['mass_gram_error']))
        grids=[]
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                rows=[r for r in self.full if (r['pair'],r['k_w'])==(pair,k)];grid=[r['beta_deg'] for r in rows]
                require(grid==sorted(set(grid)) and grid[0]==0 and grid[-1]<=math.degrees(.005*30),'Invalid archive nodes');grids.append(grid)
                for r in rows:require(abs(r['F_omega']-(r['Lambda1']/rows[0]['Lambda1'])**2)<1e-13,'Full normalization')
        require(all(g==grids[0] for g in grids),'Different archived grids')
        identity=dict(config=self.cfg,physical_sha256=self.a.sha(physics.__file__),versions=dict(numpy=np.__version__,scipy=scipy.__version__))
        self.key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16];self.folder=BASE/self.key
        self.folder.mkdir(parents=True,exist_ok=True);path=self.folder/'contract.json'
        if path.exists():require(json.loads(path.read_text())==identity,'Immutable input contract')
        else:self.a.write(path,identity)
        self.csv('reused_full_points',self.full);self.csv('archive_checks',self.archive_checks)
        self.csv('beta_xi',[dict(beta_deg=b,beta_rad=math.radians(b),xi=math.radians(b)/.005) for b in grids[0]])

    def check(self,path,digest):
        require(self.a.sha(path)==digest,'Source changed: '+str(path));self.hashes[path.relative_to(ROOT).as_posix()]=digest

    def csv(self,name,rows):self.a.csvwrite(self.folder/(name+'.csv'),rows)

    def leading(self,pair,k,xi,p=24):return physics.Leading(*self.cfg['cases'][pair],mu=self.cfg['mu'],order=p,quadrature=self.cfg['quadrature'],k_w=k,analytic_curvature=True)

    def point(self,pair,k,xi,p=24):
        identity=dict(pair=pair,k_w=k,xi=float(xi),order=p,quadrature=self.cfg['quadrature'],analytic_curvature=True)
        key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16];path=self.folder/'leading_points'/(key+'.json')
        if path.exists():
            record=json.loads(path.read_text());require(record['identity']==identity and record['status']=='PASS','Unaccepted saved leading point')
        else:
            model=self.leading(pair,k,xi,p);solution=model.solve(xi);self.new+=1;d=solution['diagnostics']
            record=dict(identity=identity,solution=solution,status='FAIL',workflow_sha256=self.a.sha(__file__))
            self.a.write(path,record)
            detail=detailed_checks(model,solution,xi);record['checks']=detail
            require(d['mass_cholesky'] and d['constraint_rank']==1 and d['P_rank']==1 and min(solution['zeta'])>0 and min(np.diff(solution['zeta']))>0,'Leading matrices or sorted prefix')
            for name in QUALITY:require(d[name]<=self.cfg['quality'][name+'_tol'],'Leading gate '+name+': '+str(identity)+' value='+str(d[name]))
            require(all(v['strong']<=self.cfg['quality']['strong_ode_tol'] and v['natural']<=self.cfg['quality']['natural_boundary_tol'] for v in detail['positions']),'Detailed leading gate')
            record['status']='PASS';self.a.write(path,record)
        return record['solution']

    def convergence(self,pair,k,xi,reason):
        hi=self.point(pair,k,xi,24);lo=self.point(pair,k,xi,20)
        error=float(max(abs(np.sqrt(np.array(lo['zeta'])/hi['zeta'])-1)))
        errorC=abs(lo['C_at_zero']/hi['C_at_zero']-1) if xi==0 else None
        require(error<=self.cfg['quality']['frequency_convergence_rtol'] and (errorC is None or errorC<=self.cfg['quality']['C_convergence_rtol']),'Order convergence')
        return dict(pair=pair,k_w=k,xi=xi,reason=reason,relative_frequency_error=error,relative_C_error=errorC,status='PASS')

    def preflight(self):
        path=self.folder/'preflight.json'
        if path.exists():require(json.loads(path.read_text())['status']=='PASS','Preflight failed');return
        oldmodel=physics.Leading(.1,10.,mu=.01,order=24,quadrature=64,k_w=1.,analytic_curvature=False)
        old=oldmodel.solve(0.);self.new+=1;archived=json.loads((self.sources['diagnostic']/'original_failed_point.json').read_text())['solution']
        self.a.write(self.folder/'old_failure_reproduction.json',dict(solution=old,status='EXPECTED_FAIL',identity=dict(pair='J0E4',k_w=1.,xi=0.,order=24,analytic_curvature=False)))
        require(old['diagnostics']['strong_ode']>1e-9 and abs(old['diagnostics']['strong_ode']/archived['diagnostics']['strong_ode']-1)<1e-10,'Old failure not reproduced')
        require(np.array_equal(old['coefficients'],archived['coefficients']),'Old eigenpair changed')
        corrected=self.point('J0E4',1.,0.);control=self.point('J0E4',1.,10.)
        for xi,sol,name in [(0.,corrected,'corrected_order_24.json'),(10.,control,'control_order_24.json')]:
            reference=json.loads((self.sources['diagnostic']/name).read_text())['solution']
            require(max(abs(np.sqrt(np.array(sol['zeta'])/reference['zeta'])-1))<1e-12,'Corrected diagnostic frequency regression')
            require(abs(sol['diagnostics']['strong_ode']/reference['diagnostics']['strong_ode']-1)<1e-6,'Corrected diagnostic strong regression')
        fixedmodel=self.leading('J0E4',1.,0.)
        for name in ('M','Kj','Ke','Kw','P','constraint','Z'):require(np.array_equal(getattr(oldmodel,name),getattr(fixedmodel,name)),'Mechanics altered: '+name)
        require(max(abs(np.sqrt(np.array(old['zeta'])/corrected['zeta'])-1))<=self.cfg['quality']['frequency_convergence_rtol'],'Old/new frequency discrepancy')
        convergence=[]
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                convergence.append(self.convergence(pair,k,0.,'all zero cases'))
                if k in self.cfg['preflight_k_w']:
                    for xi in self.cfg['preflight_xi']:convergence.append(self.convergence(pair,k,xi,'preflight nonzero'))
                self.csv('preflight_convergence',convergence)
            print('Leading preflight accepted: '+pair,flush=True)
        # Existing rigid coefficient regressions from the accepted axial study.
        refs=read_csv(ROOT/'results/isotropic_two_springs_rl_timoshenko/external_axial_spring/64b26f99cb3c7aba/leading_regression.csv')
        regression=[]
        for pair in self.cfg['cases']:
            ref=next(r for r in refs if r['pair']==pair);sol=self.point(pair,'RIGID',0.)
            error=max(abs(sol['zeta'][0]/float(ref['zeta0'])-1),abs(sol['C_at_zero']/float(ref['C_rigid'])-1))
            require(error<1e-10,'Canonical rigid regression');regression.append(dict(pair=pair,error=error,status='PASS'))
        for k,target in [('RIGID',.03945143),(10.,.03808395)]:require(abs(self.point('J0E4',k,0.)['C_at_zero']-target)<=.5e-8,'Rounded C control')
        self.csv('rigid_regression',regression)
        self.a.write(path,dict(status='PASS',cases=10,order_comparisons=len(convergence),old_failure=old['diagnostics']['strong_ode'],corrected=corrected['diagnostics']['strong_ode'],control_xi10=control['diagnostics']['strong_ode']))

    def run(self):
        self.preflight();curvature=[];local=[];convergence=read_csv(self.folder/'preflight_convergence.csv');rows=[]
        zeros={}
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                zero=self.point(pair,k,0.);zeros[pair,k]=zero
                model=self.leading(pair,k,0.);c=np.asarray(zero['coefficients'])[:,0];mass=float(c@model.M@c);y=float(model.y@c);C=y*y/(4*zero['zeta'][0]*mass)
                require(abs(C/zero['C_at_zero']-1)<1e-12,'Mass normalization of C')
                curvature.append(dict(pair=pair,k_w=k,zeta0=zero['zeta'][0],Lambda0=zero['Lambda'][0],omega_lim0=zero['omega_lim'][0],node0=y,node0_abs=abs(y),mass=mass,node0_squared=y*y,C=C,gap12=zero['diagnostics']['gap12_relative']))
                for xi in self.cfg['local_C_xi']:
                    sol=self.point(pair,k,xi);estimate=math.log(math.sqrt(sol['zeta'][0]/zero['zeta'][0]))/xi**2;error=abs(estimate/C-1)
                    local.append(dict(pair=pair,k_w=k,xi=xi,C=C,logF_over_xi2=estimate,relative_error=error))
                    require(error<=self.cfg['local_C_rtol'],'Local C check')
        self.csv('leading_curvature',curvature);self.csv('local_curvature_checks',local)
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                for r in [v for v in self.full if (v['pair'],v['k_w'])==(pair,k)]:
                    sol=self.point(pair,k,r['xi']);zero=zeros[pair,k];fullzero=next(v for v in self.full if (v['pair'],v['k_w'],v['beta_deg'])==(pair,k,0.))
                    compare=frequency_comparison(r['Lambda1'],fullzero['Lambda1'],sol['zeta'][0],zero['zeta'][0])
                    rows.append(dict(pair=pair,k_w=k,beta_deg=r['beta_deg'],xi=r['xi'],**{f'zeta{i}':z for i,z in enumerate(sol['zeta'],1)},**compare,gap12_lim=sol['diagnostics']['gap12_relative'],gap12_full=r['gap12_full']))
                    if r['beta_deg'] in self.cfg['convergence_beta_deg']:convergence.append(self.convergence(pair,k,r['xi'],'archive endpoint/intermediate'))
                self.csv('comparison',rows);self.csv('convergence',convergence);print('Leading archive nodes accepted: '+pair+' '+str(k),flush=True)
        stats=[];growth=[]
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                group=[r for r in rows if (r['pair'],r['k_w'])==(pair,k)];worst=max(group,key=lambda r:r['delta_percent']);last=group[-1]
                convergence.append(self.convergence(pair,k,worst['xi'],'maximum model discrepancy'))
                nonzero=[r['e'] for r in group if r['xi']>0]
                stats.append(dict(pair=pair,k_w=k,nodes=len(group),max_delta_percent=worst['delta_percent'],median_delta_percent=float(np.median([r['delta_percent'] for r in group])),
                    max_beta_deg=worst['beta_deg'],max_xi=worst['xi'],last_e=last['e'],last_delta_percent=last['delta_percent'],min_gap12_lim=min(r['gap12_lim'] for r in group),min_gap12_full=min(r['gap12_full'] for r in group),e_sign_change=min(nonzero)<0<max(nonzero)))
                at3=next(r for r in group if r['beta_deg']==3.)
                growth.append(dict(pair=pair,k_w=k,beta_deg=3.,xi=at3['xi'],gain_full_percent=100*(at3['F_full']-1),gain_lim_percent=100*(at3['F_lim']-1)))
        self.csv('convergence',convergence);self.csv('error_statistics',stats);self.csv('growth_at_3deg',growth)
        monotonic=[];by={(r['pair'],r['k_w'],r['xi']):r for r in rows};grid=sorted({r['xi'] for r in rows})
        for pair in self.cfg['cases']:
            for xi in grid:
                for lo,hi in zip([0.,.1,1.,10.],[.1,1.,10.,'RIGID']):
                    diff=[math.sqrt(by[pair,hi,xi][f'zeta{i}']/by[pair,lo,xi][f'zeta{i}'])-1 for i in (1,2,3)]
                    require(min(diff)>=-self.cfg['quality']['monotonicity_rtol'],'Absolute stiffness monotonicity')
                    monotonic.append(dict(pair=pair,parameter='k_w',fixed=xi,lower=lo,upper=hi,min_relative_frequency_increment=min(diff),status='PASS'))
            for k in self.cfg['k_w_values']:
                for lo,hi in zip(grid,grid[1:]):
                    diff=[math.sqrt(by[pair,k,hi][f'zeta{i}']/by[pair,k,lo][f'zeta{i}'])-1 for i in (1,2,3)]
                    require(min(diff)>=-self.cfg['quality']['monotonicity_rtol'],'Absolute xi monotonicity')
                    monotonic.append(dict(pair=pair,parameter='xi',fixed=k,lower=lo,upper=hi,min_relative_frequency_increment=min(diff),status='PASS'))
        self.csv('monotonicity',monotonic)
        points=[json.loads(p.read_text()) for p in sorted((self.folder/'leading_points').glob('*.json'))]
        quality=[dict(**r['identity'],**{key:r['solution']['diagnostics'][key] for key in QUALITY},gap12=r['solution']['diagnostics']['gap12_relative'],status=r['status']) for r in points]
        self.csv('leading_quality',quality)
        self.a.write(self.folder/'numerical_summary.json',dict(status='PASS',constructions=10,comparison_nodes=len(rows),nodes_per_construction=len(grid),xi_max=max(grid),beta_max=max(r['beta_deg'] for r in rows),
            accepted_leading_solves=len(points),old_failure_regression_solves=1,full_model_solves=0,quality_maxima={key:max(r[key] for r in quality) for key in QUALITY},
            order_frequency_error_max=max(float(r['relative_frequency_error']) for r in convergence),C_order_error_max=max(float(r['relative_C_error']) for r in convergence if r['relative_C_error'] not in (None,'')),
            monotonicity_checks=len(monotonic),max_delta_percent=max(r['max_delta_percent'] for r in stats),max_local_C_error=max(r['relative_error'] for r in local)))
        return rows

    def finish(self,stage,status):
        for name,digest in self.hashes.items():require(self.a.sha(ROOT/name)==digest,'Protected source changed')
        self.a.write(self.folder/'source_hashes.json',self.hashes)
        self.a.write(self.folder/'source_integrity.json',{name:diagnostic.verify_source(path) for name,path in self.sources.items()})
        path=self.folder/'executions.json';history=json.loads(path.read_text()) if path.exists() else []
        history.append(dict(stage=stage,status=status,new_leading_solves=self.new,new_full_model_solves=0,seconds=time.perf_counter()-self.started,workflow_sha256=self.a.sha(__file__)))
        self.a.write(path,history)

def render(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    rows=read_csv(folder/'comparison.csv');curves=read_csv(folder/'leading_curvature.csv');cfg=json.loads((folder/'contract.json').read_text())['config']
    plt.rcParams.update({'font.size':10,'axes.labelsize':11,'legend.fontsize':8,'pdf.fonttype':42,'font.family':'DejaVu Sans'})
    comma=FuncFormatter(lambda x,pos:f'{x:g}'.replace('.',','))
    colors=['#222222','#0072B2','#D55E00','#009E73','#8B4A9C'];styles=['-','--','-.',':',(0,(5,1,1,1))];markers=['o','s','^','v','D']
    labels=['Точная опора',r'$\kappa_w=10$',r'$\kappa_w=1$',r'$\kappa_w=0{,}1$',r'$\kappa_w=0$']
    fig,axes=plt.subplots(2,2,figsize=(8.1,6.2),sharex=True,layout='constrained')
    for column,(pair,springs) in enumerate(cfg['cases'].items()):
        for j,k in enumerate(cfg['k_w_values']):
            group=[r for r in rows if r['pair']==pair and stiffness(r['k_w'])==k];x=[float(r['xi']) for r in group]
            axes[0,column].plot(x,[float(r['F_lim']) for r in group],color=colors[j],ls=styles[j],label=labels[j],lw=1.35)
            axes[0,column].plot(x,[float(r['F_full']) for r in group],ls='none',marker=markers[j],ms=3.1,mfc='none',mec=colors[j],mew=.65)
            axes[1,column].plot(x,[float(r['delta_percent']) for r in group],color=colors[j],ls=styles[j],lw=1.35)
        text=r'$(\kappa_{\mathrm{j}},\kappa_{\mathrm{e}})=('+f'{springs[0]:g};{springs[1]:g}'.replace('.',r'{,}')+r')$'
        axes[0,column].text(.04,.96,text,transform=axes[0,column].transAxes,va='top',fontsize=10)
        axes[1,column].set_xlabel(r'$\xi=\beta_{\mathrm{rad}}/\varepsilon$')
    axes[0,0].set_ylabel(r'$F_{\mathrm{full}},\ F_{\mathrm{lim}}$');axes[1,0].set_ylabel(r'$\delta$, %')
    ymin=0.;ymax=max(float(r['delta_percent']) for r in rows)*1.08
    for i,a in enumerate(axes.flat):
        a.grid(alpha=.2);a.xaxis.set_major_formatter(comma);a.yaxis.set_major_formatter(comma);a.set_xlim(0,30)
        if i>=2:a.set_ylim(ymin,ymax)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='outside lower center',ncol=3,frameon=False)
    fig.savefig(folder/(FIGURES[0]+'.pdf'));fig.savefig(folder/(FIGURES[0]+'.png'),dpi=320);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(7.4,3.7),sharey=True,gridspec_kw={'width_ratios':[1,3,1]},layout='constrained')
    for j,(pair,springs) in enumerate(cfg['cases'].items()):
        data={stiffness(r['k_w']):float(r['C']) for r in curves if r['pair']==pair};label=r'$('+f'{springs[0]:g};{springs[1]:g}'.replace('.',r'{,}')+r')$'
        for ax,key in [(axes[0],0.),(axes[2],'RIGID')]:ax.plot([0],[data[key]],marker=markers[j],color=colors[j+1],ms=6)
        axes[1].semilogx([.1,1,10],[data[k] for k in [.1,1,10]],styles[j],marker=markers[j],color=colors[j+1],label=label)
    axes[0].set_xticks([0],['0']);axes[2].set_xticks([0],['Точная\nопора']);axes[1].set_xticks([.1,1,10],['0,1','1','10'])
    axes[0].set_ylabel(r'$C$');axes[1].set_xlabel(r'$\kappa_w$');axes[1].legend(title=r'$(\kappa_{\mathrm{j}},\kappa_{\mathrm{e}})$')
    for a in axes:a.grid(alpha=.2);a.yaxis.set_major_formatter(comma)
    fig.savefig(folder/(FIGURES[1]+'.pdf'));fig.savefig(folder/(FIGURES[1]+'.png'),dpi=320);plt.close(fig)

def verify(folder):
    return diagnostic.verify_source(folder)

def postprocess_saved(folder):
    """Mechanical ratios and archived 0--3 degree controls; no spectral calls."""
    a=api();cfg=json.loads((folder/'contract.json').read_text())['config'];source=ROOT/cfg['sources']['full']['folder']
    old={(r['pair'],stiffness(r['k_w'])):float(r['gain_0_3_percent']) for r in read_csv(source/'summary.csv')}
    growth=read_csv(folder/'growth_at_3deg.csv');curves=read_csv(folder/'leading_curvature.csv');mechanics=[]
    for pair in cfg['cases']:
        data={stiffness(r['k_w']):r for r in curves if r['pair']==pair};rigid=data['RIGID'];free=data[0.]
        mechanics.append(dict(pair=pair,node_square_free_to_rigid=float(free['node0_squared'])/float(rigid['node0_squared']),
            zeta_free_to_rigid=float(free['zeta0'])/float(rigid['zeta0']),C_free_to_rigid=float(free['C'])/float(rigid['C']),
            C_min_k_w=min(data,key=lambda k:float(data[k]['C'])),growth_min_k_w=min([r for r in growth if r['pair']==pair],key=lambda r:float(r['gain_lim_percent']))['k_w']))
    checks=[dict(pair=r['pair'],k_w=r['k_w'],absolute_percent_difference=abs(float(r['gain_full_percent'])-old[r['pair'],stiffness(r['k_w'])])) for r in growth]
    require(max(r['absolute_percent_difference'] for r in checks)<1e-10,'Archived full gain mismatch')
    a.csvwrite(folder/'mechanical_ratios.csv',mechanics);a.csvwrite(folder/'full_gain_regression.csv',checks)
    source_cfg=json.loads((source/'contract.json').read_text())['config'];ref_info=source_cfg['sources']['rigid'];ref=ROOT/ref_info['folder']
    require(a.sha(ref/'manifest.json')==ref_info['manifest_sha256'],'Rigid regression source changed')
    manifest=json.loads((ref/'manifest.json').read_text());name='leading_regression.csv'
    require(a.sha(ref/name)==manifest['artifacts'][name],'Rigid coefficient source changed')
    a.write(folder/'regression_source.json',dict(file=(ref/name).relative_to(ROOT).as_posix(),sha256=a.sha(ref/name),manifest_sha256=a.sha(ref/'manifest.json'),status='PASS'))

def finalize(folder):
    a=api();summary=json.loads((folder/'numerical_summary.json').read_text());require(summary['status']=='PASS','Cannot finalize incomplete study')
    for name in FIGURES:
        for suffix in ('.pdf','.png'):require((folder/(name+suffix)).is_file(),'Missing figure')
    a.write(folder/'manifest.json',dict(status='PASS',full_model_recalculated=False,artifacts={p.relative_to(folder).as_posix():a.sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}))

def main(args):
    if args.stage in ('plot_only','finalize'):
        require(args.folder is not None,'Explicit new archive folder required')
        require(args.folder.resolve().is_relative_to(BASE.resolve()),'Only the new leading-resume family may be written')
        with forbid_full_calculations(leading_too=True):
            if args.stage=='plot_only':postprocess_saved(args.folder);render(args.folder)
            else:finalize(args.folder);print(verify(args.folder))
        return
    require(args.stage in ('preflight','run'),'Use preflight/run/plot_only/finalize for leading resume')
    with forbid_full_calculations():
        exp=Experiment();status='PASS'
        print('ARCHIVE',exp.folder,flush=True)
        try:
            if (exp.folder/'failure.json').exists():raise RuntimeError('Saved failure requires addressed diagnosis before continuation')
            if args.stage=='preflight':exp.preflight()
            else:exp.run();postprocess_saved(exp.folder)
        except Exception as exc:
            status='PARTIAL';exp.a.write(exp.folder/'failure.json',dict(status=status,stage=args.stage,error=str(exc)));raise
        finally:exp.finish(args.stage,status)
