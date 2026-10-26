"""Bounded transverse-support study, reusing canonical and axial-study checks.

New physical trace space and leading contract justify a separate diagnostic
entry point. No prior solver, CLI or archive is modified.
"""
import hashlib
from contextlib import contextmanager
import json
import math
from pathlib import Path
import time
import numpy as np
import scipy
from scripts.lib import two_springs_external_axial_spring as previous
from scripts.lib import isotropic_two_springs_external_transverse_spring as physics

ROOT=previous.ROOT
BASE=previous.lengths.BASE/'external_transverse_spring'
CONFIG=ROOT/'data/input/isotropic_two_springs_external_transverse_spring.json'
api=previous.api
require=previous.require
read_csv=previous.read_csv
stiffness=previous.chi
QUALITY=previous.QUALITY
FIGURES=('01_first_frequency','02_normalized_growth','03_leading_comparison','04_initial_growth_coefficient')
PARTIAL='PARTIAL_FULL_PASS_LEADING_STOP'


class Experiment(previous.Experiment):
    """Reuse source hashes, CSV I/O and per-root validation, not axial physics."""
    def __init__(self):
        self.a=api();self.cfg=json.loads(CONFIG.read_text());self.inputs={};self.manifests={};self.sources={}
        self.started=time.perf_counter();self.leading_new=0;self.ritz_new=0
        for name,digest in self.cfg['frozen_files'].items():self.check(ROOT/name,digest)
        for name,s in self.cfg['sources'].items():
            folder=ROOT/s['folder'];self.check(folder/'manifest.json',s['manifest_sha256'])
            manifest=json.loads((folder/'manifest.json').read_text());require(manifest['status']=='PASS','Source status')
            self.sources[name]=folder;self.manifests[name]=manifest
        identity=dict(config=self.cfg,physical_sha256=self.a.sha(physics.__file__),versions=dict(numpy=np.__version__,scipy=scipy.__version__))
        self.key=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16]
        self.folder=BASE/self.key;self.folder.mkdir(parents=True,exist_ok=True)
        path=self.folder/'contract.json'
        if path.exists():require(json.loads(path.read_text())==identity,'Immutable contract differs')
        else:self.a.write(path,identity)
        self.run=self.a.Run(self.cfg,self.folder,engine=physics)
        path=self.folder/'pole_catalogs.json'
        if path.exists():self.run.catalogs={float(k):v for k,v in json.loads(path.read_text()).items()}
        cfg=json.loads(self.source('rigid','contract.json').read_text())['config']
        require(cfg['material']==self.cfg['material'] and cfg['cases']==self.cfg['cases'] and cfg['mu']==self.cfg['mu'] and cfg['variant']=='RL_FULL','Rigid source physical contract')
        self.rigid={(r['pair'],float(r['beta_deg'])):r for r in read_csv(self.source('rigid','spectrum.csv')) if r['chi_u']=='RIGID'}
        self.growth={(r['ID'],float(r['beta_deg'])):r for r in read_csv(self.source('growth','spectrum.csv')) if r['ID'] in self.cfg['cases']}
        require(all((p,b) in self.rigid for p in self.cfg['cases'] for b in self.cfg['map']['parameter_grid']),'Incomplete rigid grid')

    def case(self,pair,k,beta):return physics.Case(float(beta),self.cfg['mu'],*self.cfg['cases'][pair],k)

    def point(self,pair,k,beta,previous_w=None,fresh=False):
        case=self.case(pair,k,beta)
        if k=='RIGID' and not fresh:
            return self.archived(self.rigid[pair,beta],physics.base.Case(beta,self.cfg['mu'],*self.cfg['cases'][pair]))
        extra={round(b+sign*n*h,10) for b in self.cfg['variation']['beta_deg'] for h in self.cfg['variation']['steps_deg'] for sign in (-1,1) for n in (1,2)}
        require(k in self.cfg['k_w_values'] and (beta in self.cfg['map']['parameter_grid'] or k in self.cfg['variation']['k_w_values'] and beta in extra),'Unrequested point')
        rec=self.run.solve(case,neighbour=previous_w);self.validate(rec,case)
        self.a.write(self.folder/'pole_catalogs.json',self.run.catalogs)
        return rec,previous.lengths.point_path(self.folder,case),'new_study'

    def monotonic(self,points,angles):
        rows=[];order=[0.,.1,1.,10.,'RIGID']
        for pair in self.cfg['cases']:
            for b in angles:
                for lo,hi in zip(order,order[1:]):
                    a=np.array([r['omega'] for r in points[pair,lo,b]['roots']]);z=np.array([r['omega'] for r in points[pair,hi,b]['roots']]);diff=z/a-1
                    require(min(diff)>=-self.cfg['policy']['same_physics_frequency_tol'],'Nonmonotone absolute frequency')
                    rows.append(dict(pair=pair,beta_deg=b,k_w_lower=lo,k_w_upper=hi,**{f'omega{k}_relative_increment':float(v) for k,v in enumerate(diff,1)},status='PASS'))
        return rows

    def preflight(self):
        path=self.folder/'preflight.json'
        if path.exists():
            r=json.loads(path.read_text());require(r['status']=='PASS','Previous preflight failure');return r
        points={};regression=[];ritz=[];structural=[];m=self.run.model
        for pair in self.cfg['cases']:
            for b in self.cfg['preflight_angles']:
                old,_,_=self.point(pair,'RIGID',b);fresh,_,_=self.point(pair,'RIGID',b,fresh=True)
                growth,_,_=self.archived(self.growth[pair,b],physics.base.Case(b,self.cfg['mu'],*self.cfg['cases'][pair]))
                error=max(abs(a['omega']/z['omega']-1) for a,z in zip(fresh['roots'],old['roots']))
                overlap=max(abs(a['omega']/z['omega']-1) for a,z in zip(growth['roots'],old['roots']))
                require(max(error,overlap)<=self.cfg['policy']['same_physics_frequency_tol'],'Rigid/archive regression')
                c=self.case(pair,'RIGID',b);omega=old['roots'][0]['omega']
                require(np.array_equal(physics.boundary(m,c,omega),physics.base.boundary(m,c,omega)),'Rigid matrix changed')
                regression.append(dict(pair=pair,beta_deg=b,relative_frequency_error=error,archive_overlap_error=overlap,status='PASS'));points[pair,'RIGID',b]=old
                for k in self.cfg['k_w_values'][1:]:
                    c=self.case(pair,k,b);p,s=physics.trace_map(m,c);ext=physics.external_matrix(m,c,1)
                    require(np.array_equal(ext[0],np.eye(6)[0]) and np.array_equal(ext[1],[0,k*m.section.EI/m.l**3,0,0,-1,0]),'External support sign or axial condition')
                    require(np.array_equal(physics.joint_matrix(m,c),physics.base.joint_matrix(m,c)),'Joint changed')
                    require(np.array_equal(p[:,-1],np.eye(12)[:,7]) and np.linalg.matrix_rank(p)==6 and s[-1,-1]==k*m.section.EI/m.l**3,'Transverse trace degree of freedom')
                    require(physics.boundary(m,c,m.omega(.4)).shape==(12,12) and np.linalg.matrix_rank(physics.boundary(m,c,m.omega(.4)))==12,'Boundary condition dimension/rank')
                    rec,_,_=self.point(pair,k,b);points[pair,k,b]=rec;samples=rec['count']['samples']
                    require(samples[0]['count']==0 and samples[1]['count']>=3,'Search bounds do not enclose prefix')
                    structural.append(dict(pair=pair,k_w=k,beta_deg=b,low_count=samples[0]['count'],upper_count=samples[1]['count'],Lambda1=rec['roots'][0]['Lambda'],Lambda3=rec['roots'][2]['Lambda'],status='PASS'))
            print('Physical preflight: '+pair,flush=True)
            for k,b in self.cfg['ritz_cases']:
                expected=np.array([r['omega'] for r in points[pair,k,b]['roots']])
                for order in self.cfg['policy']['ritz_orders']:
                    w,coef,d=physics.ritz(m,self.case(pair,k,b),order);self.ritz_new+=1;error=float(max(abs(w/expected-1)))
                    ritz.append(dict(pair=pair,k_w=k,beta_deg=b,order=order,relative_frequency_error=error,**d));self.csv('ritz_checks',ritz)
                    require(error<=self.cfg['policy']['ritz_frequency_tol'] and d['mass_gram_error']<=self.cfg['policy']['mass_gram_tol'] and d['energy_error']<=self.cfg['policy']['energy_tol'] and d['constraint_residual']<=self.cfg['policy']['boundary_tol'],'Direct energy Ritz comparison')
        self.csv('regression',regression);self.csv('preflight_points',structural);self.csv('preflight_monotonicity',self.monotonic(points,self.cfg['preflight_angles']))
        result=dict(status='PASS',rigid_regressions=len(regression),finite_points=len(structural),ritz_solves=len(ritz),max_rigid_error=max(r['relative_frequency_error'] for r in regression),max_ritz_error=max(r['relative_frequency_error'] for r in ritz),min_Lambda1=min(r['Lambda1'] for r in structural),max_Lambda3=max(r['Lambda3'] for r in structural),zero_mechanism=False,
            range_evidence='Both bound counts checked for every preliminary case; zero-energy argument and positive definite Ritz K/M',pole_evidence='Fresh fixed-endpoint single-arm DtN pole catalogs; external support enters only admissible traces and spring matrix')
        self.a.write(path,result);return result

    def finish(self,stage,status='PASS'):
        for name,digest in self.inputs.items():require(self.a.sha(ROOT/name)==digest,'Protected source changed: '+name)
        path=self.folder/'source_hashes.json';old=json.loads(path.read_text()) if path.exists() else {};old.update(self.inputs);self.a.write(path,old)
        path=self.folder/'executions.json';runs=json.loads(path.read_text()) if path.exists() else []
        runs.append(dict(stage=stage,status=status,seconds=time.perf_counter()-self.started,new_full_prefixes=self.run.new_points,new_full_roots=self.run.new_roots,new_pole_roots=self.run.pole_roots,new_ritz_solves=self.ritz_new,new_leading_solves=self.leading_new,workflow_sha256=self.a.sha(__file__)))
        self.a.write(path,runs)

    def sweep(self):
        rows=[];quality=[];neighbours=[];points={};path=self.folder/'triggered_counts.json'
        audits=json.loads(path.read_text()) if path.exists() else {}
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                prev=None
                for b in self.cfg['map']['parameter_grid']:
                    rec,file,origin=self.point(pair,k,b,prev);points[pair,k,b]=rec
                    w=np.array([r['omega'] for r in rec['roots']]);jump=float(max(abs(w/prev-1))) if prev is not None else 0.
                    triggered=jump>self.cfg['policy']['neighbour_relative_jump']
                    if triggered:
                        key=f'{pair}|{k}|{b}'
                        if key not in audits:
                            case=self.case(pair,k,b);cats=self.run.cats(self.run.model,case)
                            test=[physics.count_at(self.run.model,case,x,cats,self.cfg['policy']) for x in (.9*w[0],(w[0]+w[1])/2,(w[1]+w[2])/2)]
                            require([r[1] for r in test]==[0,1,2],'Triggered prefix count')
                            audits[key]=[r[2] for r in test];self.a.write(path,audits)
                    neighbours.append(dict(pair=pair,k_w=k,beta_deg=b,jump=jump,triggered=triggered,status='PASS'));prev=w
                    if b==0:zero=w[0]
                    rows.append(dict(pair=pair,k_w=k,beta_deg=b,**{f'{name}{i}':r[name] for i,r in enumerate(rec['roots'],1) for name in ('omega','Lambda')},F_omega=w[0]/zero,gap12_relative=w[1]/w[0]-1,origin=origin,point_file=file.relative_to(ROOT).as_posix(),point_sha256=self.a.sha(file),shape_sha256=rec['shape_sha256']))
                    for i,r in enumerate(rec['roots'],1):
                        d=r['diagnostics'];quality.append(dict(pair=pair,k_w=k,beta_deg=b,position=i,physical=float(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))),ode=d['ode_error'],energy=d['energy_error'],singular=d['singular_ratio'],mass_gram=rec['mass_gram_error'],condition=d['nonzero_condition'],count_lower=r['counts'][0],count_upper=r['counts'][1],low_count=rec['count']['samples'][0]['count'],upper_count=rec['count']['samples'][1]['count'],status='PASS'))
                for name,data in [('spectrum',rows),('quality',quality),('neighbours',neighbours)]:self.csv(name,data)
                print(f'Map: {pair}, k_w={k}: 51 accepted nodes',flush=True)
        self.csv('monotonicity',self.monotonic(points,self.cfg['map']['parameter_grid']))
        return rows,quality

    def analyze(self,rows):
        summary=[];mechanics=[]
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                rr=[r for r in rows if (r['pair'],r['k_w'])==(pair,k)];at={r['beta_deg']:r for r in rr}
                summary.append(dict(pair=pair,k_w=k,Lambda0=at[0]['Lambda1'],Lambda3=at[3]['Lambda1'],Lambda10=at[10]['Lambda1'],Lambda20=at[20]['Lambda1'],gain_0_3_percent=100*(at[3]['F_omega']-1),gain_3_10_percent=100*(at[10]['F_omega']/at[3]['F_omega']-1),F10=at[10]['F_omega'],F20=at[20]['F_omega'],min_gap12=min(r['gap12_relative'] for r in rr),gap_min_beta=min(rr,key=lambda r:r['gap12_relative'])['beta_deg'],S0_per_radian=math.log(at[5]['F_omega'])/math.radians(5),beta_pl_status='not_established_no_tail'))
        m=self.run.model
        for r in rows:
            if r['beta_deg'] not in self.cfg['mechanics_angles']:continue
            file=ROOT/r['point_file'];self.check(file.with_suffix('.npz'),r['shape_sha256'])
            with np.load(file.with_suffix('.npz')) as data:arms=[{f:data[f'mode1_arm{i}_{f}'].copy() for f in ('x','weights','y','dy','ends')} for i in (1,2)]
            mass=physics.inner(m,dict(arms=arms),dict(arms=arms));ec=physics.energy_components(m,self.case(r['pair'],r['k_w'],r['beta_deg']),arms)
            require(abs(mass-1)<1e-11 and abs(ec.sum()/(mass*r['omega1']**2)-1)<self.cfg['policy']['energy_tol'],'Mechanical normalization/energy')
            u,w,_,N,Q,_=arms[0]['ends'][1];w20,Q20=arms[1]['ends'][0,[1,4]]
            arm_bend=[float(a['weights']@(m.section.EI*a['dy'][:,2]**2)/ec.sum()) for a in arms]
            mechanics.append(dict(pair=r['pair'],k_w=r['k_w'],beta_deg=r['beta_deg'],mass=mass,V=float(ec.sum()/2),d_joint=math.hypot(u,w),w20=w20,Q20=Q20,u1=u,w1=w,N1=N,Q1=Q,W_N=N*w,W_Q=-Q*u,s_per_radian=(N*w-Q*u)/ec.sum(),arm1_bend_fraction=arm_bend[0],arm2_bend_fraction=arm_bend[1],**{name:float(e/ec.sum()) for name,e in zip(('axial_fraction','bend_fraction','shear_fraction','joint_rotation_fraction','end_rotation_fraction','transverse_spring_fraction'),ec)}))
        self.csv('summary',summary);self.csv('mechanics',mechanics)
        return summary,mechanics

    def variation_checks(self,mechanics):
        cfg=self.cfg['variation'];rows=[];points={};m=self.run.model
        for r in mechanics:
            pair,k,b=r['pair'],r['k_w'],r['beta_deg']
            if k not in cfg['k_w_values'] or b not in cfg['beta_deg']:continue
            for h in cfg['steps_deg']:
                values=[]
                for offset in (-2,-1,1,2):
                    beta=round(b+offset*h,10);rec,path,origin=self.point(pair,k,beta)
                    values.append(math.log(rec['roots'][0]['omega']))
                    points[pair,k,beta]=dict(pair=pair,k_w=k,beta_deg=beta,omega1=rec['roots'][0]['omega'],point_file=path.relative_to(ROOT).as_posix(),point_sha256=self.a.sha(path),shape_sha256=rec['shape_sha256'])
                fd=(values[0]-8*values[1]+8*values[2]-values[3])/(12*math.radians(h));error=abs(fd-r['s_per_radian'])
                rows.append(dict(pair=pair,k_w=k,beta_deg=b,step_deg=h,s_variational_per_radian=r['s_per_radian'],s_fd_per_radian=fd,absolute_error=error,relative_error=error/max(abs(r['s_per_radian']),1e-30)))
                self.csv('variation_checks',rows);self.csv('variation_points',list(points.values()))
            require(rows[-1]['absolute_error']<=cfg['derivative_atol']+cfg['derivative_rtol']*abs(r['s_per_radian']),'Variational slope/finite difference disagreement')
            require(rows[-1]['absolute_error']<=rows[-2]['absolute_error']+cfg['derivative_atol'],'No finite-difference stabilization')
        return rows

    def leading_solution(self,pair,k,xi,order):
        ident=dict(pair=pair,k_w=k,xi=xi,order=order,quadrature=self.cfg['leading']['quadrature'])
        key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16];path=self.folder/'leading_points'/(key+'.json')
        if path.exists():
            rec=json.loads(path.read_text());require(rec['identity']==ident and rec['status']=='PASS','Leading cache identity');solution=rec['solution']
        else:
            model=physics.Leading(*self.cfg['cases'][pair],mu=self.cfg['mu'],order=order,quadrature=ident['quadrature'],k_w=k)
            solution=model.solve(xi);self.leading_new+=1;d=solution['diagnostics']
            passed=all(d[n]<=self.cfg['leading']['quality'][n+'_tol'] for n in QUALITY)
            self.a.write(path,dict(identity=ident,solution=solution,status='PASS' if passed else 'FAIL'))
        d=solution['diagnostics'];require(all(d[n]<=self.cfg['leading']['quality'][n+'_tol'] for n in QUALITY),'Leading physical gate: '+str(ident)+' '+str(d))
        require(d['mass_cholesky'] and d['constraint_rank']==1 and d['P_rank']==1 and min(solution['zeta'])>0 and min(np.diff(solution['zeta']))>0,'Leading matrix/sorted prefix')
        return solution

    def leading_study(self,full):
        cfg=self.cfg['leading'];refs=read_csv(self.source('rigid','leading_regression.csv'))
        regression=[];zeros={};curvature=[];local=[];convergence=[];comparison=[];monotonic=[]
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                lo=self.leading_solution(pair,k,0.,cfg['order']);hi=self.leading_solution(pair,k,0.,cfg['check_order']);zeros[pair,k]=lo
                errC=abs(lo['C_at_zero']/hi['C_at_zero']-1);require(errC<=cfg['quality']['C_convergence_rtol'],'C order convergence')
                curvature.append(dict(pair=pair,k_w=k,zeta0=lo['zeta'][0],Lambda0=lo['Lambda'][0],node0=lo['node_value'][0],C=lo['C_at_zero'],gap12=lo['diagnostics']['gap12_relative'],relative_C_order_error=errC))
                if k=='RIGID':
                    ref=next(r for r in refs if r['pair']==pair);error=max(abs(lo['zeta'][0]/float(ref['zeta0'])-1),abs(lo['C_at_zero']/float(ref['C_rigid'])-1));require(error<1e-10,'Rigid leading regression')
                    regression.append(dict(pair=pair,relative_zeta_C_error=error,status='PASS'))
                for xi in cfg['local_C_xi']:
                    sol=self.leading_solution(pair,k,xi,cfg['check_order']);estimate=math.log(math.sqrt(sol['zeta'][0]/hi['zeta'][0]))/xi**2
                    error=abs(estimate/hi['C_at_zero']-1);require(error<=cfg['local_C_rtol'],'Local curvature coefficient check')
                    local.append(dict(pair=pair,k_w=k,xi=xi,C=hi['C_at_zero'],logF_over_xi2=estimate,relative_error=error))
        self.csv('leading_regression',regression);self.csv('leading_curvature',curvature);self.csv('local_curvature_checks',local)
        for r in full:
            pair,k,b=r['pair'],r['k_w'],r['beta_deg'];xi=math.radians(b)/self.cfg['material']['epsilon']
            if xi>cfg['xi_max']:continue
            sol=self.leading_solution(pair,k,xi,cfg['order']);zero=zeros[pair,k]
            if b in cfg['check_beta_deg']:
                high=self.leading_solution(pair,k,xi,cfg['check_order']);error=float(max(abs(np.sqrt(np.array(sol['zeta'])/high['zeta'])-1)))
                require(error<=cfg['quality']['frequency_convergence_rtol'],'Leading order convergence')
                convergence.append(dict(pair=pair,k_w=k,beta_deg=b,xi=xi,relative_frequency_error=error,status='PASS'))
            Fl=math.sqrt(sol['zeta'][0]/zero['zeta'][0]);e=r['F_omega']/Fl-1
            comparison.append(dict(pair=pair,k_w=k,beta_deg=b,xi=xi,**{f'zeta{i}':float(z) for i,z in enumerate(sol['zeta'],1)},F_full=r['F_omega'],F_lim=Fl,e=e,delta_percent=100*abs(e),gap12_lim=sol['diagnostics']['gap12_relative']))
        by={(r['pair'],r['k_w'],r['beta_deg']):r for r in comparison};grid=sorted(set(r['beta_deg'] for r in comparison))
        for pair in self.cfg['cases']:
            for b in grid:
                for lo,hi in zip([0.,.1,1.,10.],[.1,1.,10.,'RIGID']):
                    diff=[math.sqrt(by[pair,hi,b][f'zeta{i}']/by[pair,lo,b][f'zeta{i}'])-1 for i in (1,2,3)]
                    require(min(diff)>=-cfg['quality']['monotonicity_rtol'],'Leading stiffness monotonicity')
                    monotonic.append(dict(pair=pair,beta_deg=b,k_w_lower=lo,k_w_upper=hi,**{f'omega{i}_relative_increment':v for i,v in enumerate(diff,1)},status='PASS'))
        quality=[]
        for path in sorted((self.folder/'leading_points').glob('*.json')):
            r=json.loads(path.read_text());quality.append(dict(**r['identity'],**{k:r['solution']['diagnostics'][k] for k in QUALITY},status=r['status']))
        stats=[]
        for pair in self.cfg['cases']:
            for k in self.cfg['k_w_values']:
                rr=[r for r in comparison if (r['pair'],r['k_w'])==(pair,k)];worst=max(rr,key=lambda r:r['delta_percent'])
                stats.append(dict(pair=pair,k_w=k,nodes=len(rr),max_delta_percent=worst['delta_percent'],median_delta_percent=float(np.median([r['delta_percent'] for r in rr])),worst_beta_deg=worst['beta_deg'],worst_xi=worst['xi'],e_min=min(r['e'] for r in rr),e_max=max(r['e'] for r in rr)))
        for name,data in [('leading_comparison',comparison),('leading_convergence',convergence),('leading_monotonicity',monotonic),('leading_quality',quality),('leading_summary',stats)]:self.csv(name,data)
        return stats,quality


def main(args):
    if args.stage in ('plot_only','finalize'):
        require(args.folder is not None,'Specify the result folder')
        folder=args.folder.resolve();verify_saved(folder)
        if args.stage=='plot_only':
            with prohibit_spectral_calls():render_full(folder)
        else:finalize_partial(folder)
        return
    exp=Experiment();print('Result: '+str(exp.folder),flush=True)
    if args.stage=='diagnose':
        with prohibit_spectral_calls():diagnose_saved(exp)
        exp.finish('diagnose',PARTIAL);numerical_manifest(exp.folder);bundle_readme(exp.folder)
        return
    require(not (exp.folder/'diagnostic_failure.json').exists(),'Recorded numerical STOP: no automatic restart; use diagnose or plot_only')
    try:
        exp.preflight()
        if args.stage=='preflight':exp.finish(args.stage);return
        rows,quality=exp.sweep();summary,mechanics=exp.analyze(rows);variation=exp.variation_checks(mechanics)
        stats,leading_quality=exp.leading_study(rows)
        exp.a.write(exp.folder/'numerical_summary.json',dict(status='PASS',map_points=len(rows),archived_map_points=sum(r['origin']=='archive' for r in rows),new_map_points=sum(r['origin']=='new_study' for r in rows),physical_maxima={k:max(r[k] for r in quality) for k in ('physical','ode','energy','singular','mass_gram','condition')},leading_quality_maxima={k:max(r[k] for r in leading_quality) for k in QUALITY},leading_solutions=len(leading_quality),tail_points=0,branch_tracking=False))
        exp.finish(args.stage);previous.numerical_manifest(exp.folder)
    except Exception as exc:
        api().write(exp.folder/'diagnostic_failure.json',dict(stage=args.stage,error=repr(exc),status='STOP'))
        exp.finish(args.stage,'FAIL');raise


@contextmanager
def prohibit_spectral_calls():
    def forbidden(*args,**kwargs):raise RuntimeError('Spectral calls forbidden after STOP and in saved-data modes')
    targets=[(physics,n) for n in ('prefix','count_at','mode','ritz')]+[(api().Run,'solve'),(physics.Leading,'solve')]
    old=[getattr(m,n) for m,n in targets]
    try:
        for m,n in targets:setattr(m,n,forbidden)
        yield
    finally:
        for (m,n),f in zip(targets,old):setattr(m,n,f)


def diagnose_saved(exp):
    """Read accepted full data and the failed leading point, without any solve."""
    folder=exp.folder;cfg=exp.cfg
    require((folder/'diagnostic_failure.json').exists(),'No recorded STOP to diagnose')
    rows=read_csv(folder/'spectrum.csv');require(len(rows)==510,'Full map incomplete')
    quality=[];paths={ROOT/r['point_file'] for r in rows}|set((folder/'points').glob('*.json'))
    for path in sorted(paths):
        rec=json.loads(path.read_text());ident={k:v for k,v in rec['case'].items() if k!='variant'}
        case=physics.Case(**ident) if 'k_w' in ident else physics.base.Case(**ident)
        exp.validate(rec,case);exp.check(path.with_suffix('.npz'),rec['shape_sha256'])
        for i,r in enumerate(rec['roots'],1):
            d=r['diagnostics'];quality.append(dict(point_file=path.relative_to(ROOT).as_posix(),k_w=ident.get('k_w','RIGID'),beta_deg=ident['beta'],position=i,
                physical=float(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))),ode=d['ode_error'],energy=d['energy_error'],singular=d['singular_ratio'],mass_gram=rec['mass_gram_error'],condition=d['nonzero_condition'],status='PASS'))
    exp.csv('all_full_quality',quality)
    fullgate=dict(status='PASS',map_points=510,archived_map_points=102,new_map_points=408,regression_points=10,variation_extra_points=len(read_csv(folder/'variation_points.csv')),
        physical_maxima={k:max(r[k] for r in quality) for k in ('physical','ode','energy','singular','mass_gram','condition')},
        monotonicity_scalar_checks=3*len(read_csv(folder/'monotonicity.csv')),zero_mechanism=False,tail_points=0,characteristic_angles=False)
    exp.a.write(folder/'full_model_gate.json',fullgate)
    points=[(p,json.loads(p.read_text())) for p in sorted((folder/'leading_points').glob('*.json'))]
    failed=[(p,r) for p,r in points if r['status']=='FAIL'];require(len(failed)==1,'Unexpected failed leading set')
    path,bad=failed[0];ident=bad['identity'];s=bad['solution'];low=next(r for _,r in points if r['identity']=={**ident,'order':cfg['leading']['order']})
    c=np.array(s['coefficients']);nodes,weights=np.polynomial.legendre.leggauss(ident['quadrature']);nodes=(nodes+1)/2;weights=weights/2
    per_mode=[]
    for j,zeta in enumerate(s['zeta']):
        num=den=0.
        for arm,L in enumerate((1-cfg['mu'],1+cfg['mu'])):
            f=physics.leading.fields(ident['order'],nodes,L,arm);ci=c[arm*ident['order']:(arm+1)*ident['order'],j]
            W=f[0]@ci+(c[-1,j] if arm==1 else 0.);D4=f[4]@ci
            num+=weights@((D4-zeta*W)**2)*L;den+=weights@(D4**2+(zeta*W)**2)*L
        per_mode.append(float(np.sqrt(num/den)))
    leadinggate=dict(status='STOP',failed_point=path.relative_to(folder).as_posix(),identity=ident,strong_ode=s['diagnostics']['strong_ode'],tolerance=cfg['leading']['quality']['strong_ode_tol'],strong_per_position=per_mode,
        max_relative_frequency_order_change=float(max(abs(np.sqrt(np.array(s['zeta'])/low['solution']['zeta'])-1))),
        max_natural_boundary=s['diagnostics']['natural_boundary'],mass_gram=s['diagnostics']['mass_gram'],
        actual_leading_solves=len(points),individual_passes=len(points)-1,coefficient_family_accepted=False,full_leading_comparison_completed=False,
        diagnosis='Fourth-derivative residual exceeds frozen tolerance in position 1; frequency agreement alone does not override it. Roundoff/conditioning is a possibility, not an established cause. No new spectral calculations after failure.')
    exp.a.write(folder/'leading_gate.json',leadinggate)
    exp.csv('leading_diagnostic_points',[dict(**r['identity'],point_file=p.relative_to(folder).as_posix(),status=r['status'],C_at_zero=r['solution']['C_at_zero'],**{k:r['solution']['diagnostics'][k] for k in QUALITY}) for p,r in points])
    # Preserve individually checked local coefficients with an explicit scope.
    local=[]
    for k in ('RIGID',10.):
        zero=next(r['solution'] for _,r in points if r['identity']['k_w']==k and r['identity']['xi']==0 and r['identity']['order']==24)
        for _,r in points:
            i=r['identity']
            if i['k_w']!=k or i['xi']==0:continue
            estimate=math.log(math.sqrt(r['solution']['zeta'][0]/zero['zeta'][0]))/i['xi']**2
            local.append(dict(pair=i['pair'],k_w=k,xi=i['xi'],C=zero['C_at_zero'],logF_over_xi2=estimate,relative_error=abs(estimate/zero['C_at_zero']-1),scope='individual_diagnostic_only_family_incomplete'))
    exp.csv('partial_curvature_diagnostics',local)
    exp.a.write(folder/'study_status.json',dict(status=PARTIAL,full_model='PASS',leading_model='STOP',accepted_figures=FIGURES[:2],missing_figures=FIGURES[2:],new_solves_after_stop=0))


def numerical_manifest(folder):
    excluded={'numerical_manifest.json','manifest.json','executions.json','source_hashes.json','delivery_checks.json','test_results.json'}
    files={p.relative_to(folder).as_posix():api().sha(p) for p in folder.rglob('*') if p.is_file() and p.suffix in ('.json','.csv','.npz') and p.name not in excluded}
    api().write(folder/'numerical_manifest.json',dict(status='INTEGRITY_SNAPSHOT',study_status=PARTIAL,artifacts=files))


def verify_saved(folder):
    m=json.loads((folder/'numerical_manifest.json').read_text());require(m['study_status']==PARTIAL,'Unexpected study scope')
    for name,digest in m['artifacts'].items():require(api().sha(folder/name)==digest,'Saved artifact changed: '+name)
    for name,digest in json.loads((folder/'source_hashes.json').read_text()).items():require(api().sha(ROOT/name)==digest,'Protected source changed: '+name)
    require(json.loads((folder/'full_model_gate.json').read_text())['status']=='PASS','Unaccepted full model')
    require(json.loads((folder/'leading_gate.json').read_text())['status']=='STOP','Leading failure must remain visible')


def render_full(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,'legend.fontsize':9,
        'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':320})
    data=read_csv(folder/'spectrum.csv');ks=['RIGID',10.,1.,.1,0.]
    colors=['#222222','#176b8b','#9a5100','#8157a1','#007b52'];styles=['-','--','-.',':',(0,(5,1,1,1))]
    labels=['Точная опора',r'$\kappa_w=10$',r'$\kappa_w=1$',r'$\kappa_w=0{,}1$',r'$\kappa_w=0$']
    comma=FuncFormatter(lambda x,pos:f'{x:g}'.replace('.',','))
    for field,ylabel,name in [('Lambda1',r'$\Lambda_1$',FIGURES[0]),('F_omega',r'$F=\omega_1(\beta)/\omega_1(0)$',FIGURES[1])]:
        fig,axes=plt.subplots(1,2,figsize=(7.2,3.55),sharey=True,layout='constrained')
        for ax,letter,pair,pairlabel in zip(axes,'аб',['J0E4','J4E0'],[r'$(0{,}1;10)$',r'$(10;0{,}1)$']):
            for k,c,s,label in zip(ks,colors,styles,labels):
                rr=[r for r in data if r['pair']==pair and stiffness(r['k_w'])==k]
                ax.plot([float(r['beta_deg']) for r in rr],[float(r[field]) for r in rr],color=c,linestyle=s,lw=1.6,label=label)
            ax.set_xlim(0,20);ax.set_xlabel(r'$\beta$, град.');ax.grid(alpha=.22,linewidth=.5)
            ax.xaxis.set_major_formatter(comma);ax.yaxis.set_major_formatter(comma)
            ax.text(.02,1.035,letter+r')  $(\kappa_{\mathrm{j}},\kappa_{\mathrm{e}})=$'+pairlabel,transform=ax.transAxes,fontsize=9)
        handles,legend_labels=axes[0].get_legend_handles_labels();fig.legend(handles,legend_labels,loc='outside lower center',frameon=False,ncol=5,fontsize=9)
        axes[0].set_ylabel(ylabel)
        fig.savefig(folder/(name+'.pdf'),bbox_inches='tight',metadata={'Title':name,'Creator':'CoupledBeams, accepted full-model observations'})
        fig.savefig(folder/(name+'.png'),dpi=320,bbox_inches='tight');plt.close(fig)


def finalize_partial(folder):
    verify_saved(folder)
    require(all((folder/(n+ext)).exists() for n in FIGURES[:2] for ext in ('.pdf','.png')),'Missing accepted full-model figures')
    require(not any((folder/(n+ext)).exists() for n in FIGURES[2:] for ext in ('.pdf','.png')),'Unqualified leading figure')
    report=ROOT/'docs/theory/isotropic_two_springs_external_transverse_spring.md'
    files={p.relative_to(folder).as_posix():api().sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}
    api().write(folder/'manifest.json',dict(status=PARTIAL,full_model_status='PASS',leading_status='STOP',artifacts=files,scientific_report=report.relative_to(ROOT).as_posix(),scientific_report_sha256=api().sha(report),physical_sha256=api().sha(physics.__file__),workflow_sha256=api().sha(__file__)))


def bundle_readme(folder):
    relative=folder.relative_to(ROOT).as_posix()
    (folder/'README.md').write_text(f'''# Поперечная податливость внешней опоры

**{PARTIAL}.** Полная модель проверена; ведущий этап остановлен
по сильной невязке. Общего статуса PASS у исследования нет.
[Научный отчёт](../../../../docs/theory/isotropic_two_springs_external_transverse_spring.md).

## Постановка и источники

RL_FULL/Тимошенко, две точные поступательные связи в узле, обе поворотные
пружины сохранены. u2(0)=0, -Q2(0)+kw*w2(0)=0; kw действует локально
поперёк второго стержня. kappa_w=kw*l^3/(EI), не kw*l/(EA).
mu=0,01; epsilon=0,005; b/h=4; nu=0,3; ks=5/6; l=E=rho=1.
J0E4 означает kappa_j=0,1, kappa_e=10; J4E0 — 10 и 0,1.
kappa_w=RIGID;10;1;0,1;0. RIGID задаётся точным ограничением.

Исходные архивы под `results/isotropic_two_springs_rl_timoshenko/`:

- `external_axial_spring/64b26f99cb3c7aba/spectrum.csv`: только 102 строки
  chi_u=RIGID, совпадающие с канонической моделью; исходные point JSON/NPZ
  указаны в `spectrum.csv`. Конечная продольная податливость не используется.
- `growth_stiffness_5x5/361e9a839ef5be30/spectrum.csv`: независимая проверка
  совпадения архивов при 0°, 3°, 6°, 10°, 20° обеих выбранных пар.
- `external_axial_spring/64b26f99cb3c7aba/leading_regression.csv`: исходные
  нулевые собственные пары ведущей задачи для точного закрепления.

Входной контракт и хеш нового физического модуля закреплены в `contract.json`;
идентификатор — первые 16 знаков SHA256 упорядоченного JSON контракта.
Все прежние источники неизменны; `source_hashes.json` фиксирует прочитанные файлы.

## Результаты и остановка

Принята карта 510 точек: 408 новых, 102 архивных. Дополнительно выполнены
10 новых жёстких регрессий, 48 точек для производной и 20 независимых
проверок Ритца полной задачи. Итого 466 новых полных префиксов по три частоты.

На десятом решении ведущей задачи остановлен контроль p=24 при J0E4,
kappa_w=1, xi=0: strong_ode=1,024892875e-9 при допуске 1e-9.
`diagnostic_failure.json` — первоначальная ошибка, `leading_gate.json` —
диагностика сохранённых коэффициентов. Девять отдельных решений прошли
свои проверки, но семейство C и сравнение F_full/F_lim не завершены.
Допуск, порядок и физические коэффициенты после отказа не менялись.
Новых спектральных решений после остановки нет.

## Воспроизведение

В текущем окружении Python: `D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe -X utf8`.
Зависимости: NumPy, SciPy, Matplotlib, pytest. Команды из корня проекта:

```powershell
python scripts/analysis/study_two_springs_external_transverse_spring.py --stage diagnose
python scripts/analysis/study_two_springs_external_transverse_spring.py --stage plot_only --folder {relative}
python -m pytest tests/test_two_springs_external_transverse_spring.py -q
python scripts/analysis/study_two_springs_external_transverse_spring.py --stage finalize --folder {relative}
```

Все перечисленные режимы работают с сохранёнными решениями. `diagnose` и
`plot_only` программно запрещают спектральные вызовы. Исходная численная
последовательность была `--stage preflight`, затем `--stage run`; журнал
`executions.json` сохраняет количество решений и момент отказа.
Повторный `run` при наличии записанного STOP не возобновляет вычисления.
Для продолжения необходим отдельный разбор численной невязки, не ослабление
допуска ради получения ожидаемого согласия. Старые данные удалять не нужно.

## Таблицы и нормировки

CSV используют десятичную точку. beta_deg — градусы, xi=beta_rad/epsilon.
Lambda^4=rho*A*omega^2*l^4/(EI), omega в единицах l=E=rho=1.

| Файл | Содержание |
| --- | --- |
| `spectrum.csv` | Первые три положительные отсортированные omega и Lambda; F_omega=omega1/собственная omega1(0), gap12_relative=(omega2-omega1)/omega1; происхождение и хеши |
| `summary.csv` | Lambda при 0,3,10,20°; прирост 0–3° относительно нуля и 3–10° относительно 3° в процентах; F10/F20; минимальный зазор; S0 на радиан; beta_pl не установлен |
| `mechanics.csv` | Первая форма с полной M=1; V, d_joint, w20,Q20; u1,w1,N1,Q1 в узле; W_N=N1*w1, W_Q=-Q1*u1; s_per_radian=(W_N+W_Q)/(2V); доли всех энергий, включая поперечную пружину |
| `variation_checks.csv` | Вариационная производная и центрированные разности четвёртого порядка при шагах 0,05° и 0,025°; абсолютное/относительное отличие |
| `variation_points.csv` | 48 фактических новых углов для разностных проверок; ссылки на полные префиксы |
| `quality.csv`, `all_full_quality.csv` | Невязки карты и всех использованных полных решений, соответственно |
| `preflight.json`, `regression.csv`, `ritz_checks.csv` | Жёсткий предел, диапазон поиска, отсутствие механизма, прямые Ritz-проверки |
| `monotonicity.csv` | Неубывание трёх абсолютных частот по kappa_w, включая RIGID; 1224 скалярных сравнения |
| `neighbours.csv`, `triggered_counts.json` | Пять срабатываний проверки соседей и разрешающие счётные проверки |
| `leading_diagnostic_points.csv` | Все десять ведущих проб с явными PASS/FAIL; это не принятая частотная карта |
| `partial_curvature_diagnostics.csv` | Только две локально проверенные конструкции J0E4; не полный набор C и не итоговая иллюстрация |
| `full_model_gate.json`, `leading_gate.json`, `study_status.json` | Раздельные статусы полной модели, ведущего контроля и всего задания |

`numerical_manifest.json` проверяет целостность и сохраняет частичный статус,
а не переименовывает отказ в PASS. `manifest.json` описывает итоговый комплект.

## Готовые рисунки и подписи

1. `01_first_frequency.pdf/png`: первая отсортированная частота Lambda1
   при поперечной податливости внешней опоры. Панели а, б — две пары поворотных
   жёсткостей; обе поворотные пружины и точное продольное закрепление сохранены.
2. `02_normalized_growth.pdf/png`: первая частота, нормированная на собственное
   значение при нулевом угле для каждого kappa_w. Монотонность абсолютных
   частот по жёсткости не переносится на эту нормированную зависимость.

Источник обоих рисунков — `spectrum.csv`; линии соединяют фактические узлы,
без сглаживания и интерполяции новых частот. PDF векторные, PNG 320 dpi.
Рисунки 3 и 4 не созданы: сравнение ведущей функции и полный набор C
не прошли предусмотренный этап контроля. Хвостов до 90°, beta_pl,
FEM, модального продолжения и замены внешней поворотной пружины нет.
''',encoding='utf8')
