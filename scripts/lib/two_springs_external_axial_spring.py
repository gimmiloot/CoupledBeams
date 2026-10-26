"""Bounded external-axial-support experiment; reuses existing count and Ritz.

The new physical boundary DOF requires a separate input/output contract.
This diagnostic workflow leaves the archived two-spring CLI untouched.
"""
from dataclasses import asdict
from fractions import Fraction
from contextlib import contextmanager
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
import scipy
from scripts.lib import two_springs_length_asymmetry as lengths
from scripts.lib import isotropic_two_springs_external_axial_spring as physics
from scripts.lib import isotropic_two_springs_longitudinal_slip as slip
from scripts.lib.isotropic_two_springs_leading import Ritz

ROOT=lengths.ROOT
CONFIG=ROOT/'data/input/isotropic_two_springs_external_axial_spring.json'
BASE=lengths.BASE/'external_axial_spring'
require=lengths.require
api=lengths.api
read_csv=lengths.read_csv
QUALITY=('symmetry','eigen_residual','energy','mass_gram','constraint','natural_boundary','strong_ode','essential')
FIGURES=('01_first_frequency','02_normalized_growth','03_spring_location','04_leading_comparison')


@contextmanager
def plots_without_calculation():
    def forbidden(*args,**kwargs):raise RuntimeError('Spectral calculation forbidden in plot_only')
    targets=[(physics,n) for n in ('prefix','mode','ritz','count_at')]+[(api().Run,'solve'),(Ritz,'solve')]
    old=[getattr(m,n) for m,n in targets]
    try:
        for m,n in targets:setattr(m,n,forbidden)
        yield
    finally:
        for (m,n),f in zip(targets,old):setattr(m,n,f)


def chi(value):return value if value=='RIGID' else float(value)


def compliance_check():
    """Exact rational constrained energy minimization, independent of matrices."""
    rows=[]
    a,b=Fraction(99,100),Fraction(101,100)
    for k in (Fraction(1),Fraction(1,10),Fraction(1,100)):
        S=a+b+1/k;force=-1/S;d1,d2,v=force*a,force*b,force/k
        energy=(d1*d1/a+d2*d2/b+k*v*v)/2
        slip_value=d1+d2+1
        joint_energy=(d1*d1/a+d2*d2/b+k*slip_value*slip_value)/2
        require(d1+d2+v==-1 and energy==joint_energy==1/(2*S),'Compliance derivation')
        rows.append(dict(chi_u=float(k),S_ax=float(S),normalized_energy=float(energy),exact_rational_identity=True))
    return rows


class Experiment:
    def validate(self,rec,case):
        w=lengths.Experiment.validate(self,rec,case)
        require(rec['count']['samples'][0]['count']==0,'Nonzero lower count')
        return w

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
        target=self.folder/'contract.json'
        if target.exists():require(json.loads(target.read_text())==identity,'Immutable contract differs')
        else:self.a.write(target,identity)
        self.run=self.a.Run(self.cfg,self.folder,engine=physics)
        target=self.folder/'pole_catalogs.json'
        if target.exists():self.run.catalogs={float(k):v for k,v in json.loads(target.read_text()).items()}
        self.rigid={(r['pair'],float(r['beta_deg'])):r for r in read_csv(self.source('rigid','complete_spectrum.csv')) if float(r['mu'])==self.cfg['mu']}
        self.joint={(r['pair'],chi(r['chi']),float(r['beta_deg'])):r for r in read_csv(self.source('joint','spectrum.csv')) if chi(r['chi']) in self.cfg['chi_values']}
        for name in ('rigid','joint'):
            c=json.loads(self.source(name,'contract.json').read_text())
            c=c.get('config',c.get('identity',{}).get('config'))
            require(c['material']==self.cfg['material'] and c['variant']=='RL_FULL','Source physical model')
        require(all((p,b) in self.rigid for p in self.cfg['cases'] for b in self.cfg['map']['parameter_grid']),'Incomplete rigid grid')

    def check(self,path,digest):
        actual=self.a.sha(path);require(actual==digest,'Changed source '+str(path))
        self.inputs[path.relative_to(ROOT).as_posix()]=actual
        return path

    def source(self,name,file):return self.check(self.sources[name]/file,self.manifests[name]['artifacts'][file])

    def csv(self,name,rows):self.a.csvwrite(self.folder/(name+'.csv'),rows)

    def case(self,pair,k,beta):return physics.Case(float(beta),self.cfg['mu'],*self.cfg['cases'][pair],k)

    def archived(self,row,case):
        path=ROOT/row['point_file'];self.check(path,row['point_sha256']);rec=json.loads(path.read_text())
        self.check(path.with_suffix('.npz'),rec['shape_sha256']);w=self.validate(rec,case)
        require(all(abs(float(row[f'omega{i}'])/v-1)<1e-13 for i,v in enumerate(w,1)),'CSV/root discrepancy')
        return rec,path,'archive'

    def point(self,pair,k,beta,previous=None,fresh=False):
        case=self.case(pair,k,beta)
        if k=='RIGID' and not fresh:
            return self.archived(self.rigid[pair,beta],physics.base.Case(beta,self.cfg['mu'],*self.cfg['cases'][pair]))
        require(k in self.cfg['chi_values'] and beta in self.cfg['map']['parameter_grid'],'Unrequested point')
        rec=self.run.solve(case,neighbour=previous);self.validate(rec,case)
        self.a.write(self.folder/'pole_catalogs.json',self.run.catalogs)
        return rec,lengths.point_path(self.folder,case),'new_study'

    def monotonic(self,points,angles):
        rows=[];ordered=[0.,.01,.1,1.,'RIGID']
        for pair in self.cfg['cases']:
            for b in angles:
                for lo,hi in zip(ordered,ordered[1:]):
                    a=np.array([r['omega'] for r in points[pair,lo,b]['roots']]);z=np.array([r['omega'] for r in points[pair,hi,b]['roots']])
                    diff=z/a-1;require(min(diff)>=-self.cfg['policy']['same_physics_frequency_tol'],'Stiffness monotonicity violated')
                    rows.append(dict(pair=pair,beta_deg=b,chi_lower=lo,chi_upper=hi,**{f'omega{k}_relative_increment':float(v) for k,v in enumerate(diff,1)},status='PASS'))
        return rows

    def preflight(self):
        file=self.folder/'preflight.json'
        if file.exists():
            r=json.loads(file.read_text());require(r['status']=='PASS','Previous preflight failure');return r
        points={};regression=[];ritz=[];structural=[];m=self.run.model
        self.csv('compliance_identity',compliance_check())
        for pair in self.cfg['cases']:
            for b in self.cfg['preflight_angles']:
                old,_,_=self.point(pair,'RIGID',b);fresh,_,_=self.point(pair,'RIGID',b,fresh=True)
                error=max(abs(a['omega']/z['omega']-1) for a,z in zip(fresh['roots'],old['roots']))
                require(error<=self.cfg['policy']['same_physics_frequency_tol'],'Rigid regression')
                c=self.case(pair,'RIGID',b);w=old['roots'][0]['omega']
                require(np.array_equal(physics.boundary(m,c,w),physics.base.boundary(m,c,w)),'Rigid boundary changed')
                regression.append(dict(pair=pair,beta_deg=b,relative_frequency_error=error,status='PASS'))
                points[pair,'RIGID',b]=old
                for k in self.cfg['chi_values'][1:]:
                    c=self.case(pair,k,b);p,s=physics.trace_map(m,c);joint=physics.joint_matrix(m,c)
                    require(np.array_equal(joint,physics.base.joint_matrix(m,c)),'Joint was altered')
                    require(p[6,-1]==1 and np.count_nonzero(p[:,-1])==1 and s[-1,-1]==k*m.section.EA/m.l,'Outer local degree of freedom')
                    row=physics.external_matrix(m,c,1)[0];require(np.array_equal(row,[k*m.section.EA/m.l,0,0,-1,0,0]),'Spring force sign')
                    rec,_,_=self.point(pair,k,b);points[pair,k,b]=rec
                    structural.append(dict(pair=pair,chi_u=k,beta_deg=b,low_count=rec['count']['samples'][0]['count'],status='PASS'))
            print('Physical preflight: '+pair,flush=True)
            for k,b in self.cfg['ritz_cases']:
                rec=points[pair,k,b];expected=np.array([r['omega'] for r in rec['roots']])
                for order in self.cfg['policy']['ritz_orders']:
                    w,coef,d=physics.ritz(m,self.case(pair,k,b),order);self.ritz_new+=1
                    error=float(max(abs(w/expected-1)))
                    ritz.append(dict(pair=pair,chi_u=k,beta_deg=b,order=order,relative_frequency_error=error,**d))
                    self.csv('ritz_checks',ritz)
                    require(error<=self.cfg['policy']['ritz_frequency_tol'] and d['mass_gram_error']<=self.cfg['policy']['mass_gram_tol'] and d['energy_error']<=self.cfg['policy']['energy_tol'] and d['constraint_residual']<=self.cfg['policy']['boundary_tol'],'Independent Ritz gate failed')
        self.csv('regression',regression);self.csv('preflight_points',structural)
        self.csv('preflight_monotonicity',self.monotonic(points,self.cfg['preflight_angles']))
        result=dict(status='PASS',rigid_regression_points=len(regression),finite_points=len(structural),ritz_solves=len(ritz),zero_mechanism=False,
            zero_mechanism_evidence='Positive stiffness and mass Cholesky in free-end Ritz; low spectral count zero; zero-energy argument in theory',
            pole_catalogs='New fixed-displacement single-arm catalogs, independent of outer structure constraints; same arms and DtN, changed admissible trace space',
            max_rigid_error=max(r['relative_frequency_error'] for r in regression),max_ritz_error=max(r['relative_frequency_error'] for r in ritz))
        self.a.write(file,result);return result

    def finish(self,stage,status='PASS'):
        for name,digest in self.inputs.items():require(self.a.sha(ROOT/name)==digest,'Protected source changed: '+name)
        file=self.folder/'source_hashes.json';old=json.loads(file.read_text()) if file.exists() else {};old.update(self.inputs);self.a.write(file,old)
        file=self.folder/'executions.json';runs=json.loads(file.read_text()) if file.exists() else []
        runs.append(dict(stage=stage,status=status,seconds=time.perf_counter()-self.started,new_full_prefixes=self.run.new_points,new_full_roots=self.run.new_roots,new_pole_roots=self.run.pole_roots,new_ritz_solves=self.ritz_new,new_leading_solves=self.leading_new,workflow_sha256=self.a.sha(__file__)))
        self.a.write(file,runs)

    def sweep(self):
        rows=[];quality=[];neighbours=[];points={};file=self.folder/'triggered_counts.json'
        audits=json.loads(file.read_text()) if file.exists() else {}
        for pair in self.cfg['cases']:
            for k in self.cfg['chi_values']:
                previous=None
                for b in self.cfg['map']['parameter_grid']:
                    rec,path,origin=self.point(pair,k,b,previous);points[pair,k,b]=rec
                    w=np.array([r['omega'] for r in rec['roots']]);jump=float(max(abs(w/previous-1))) if previous is not None else 0.
                    triggered=jump>self.cfg['policy']['neighbour_relative_jump']
                    if triggered:
                        key=f'{pair}|{k}|{b}'
                        if key not in audits:
                            c=self.case(pair,k,b);cats=self.run.cats(self.run.model,c)
                            tests=[physics.count_at(self.run.model,c,x,cats,self.cfg['policy']) for x in (.9*w[0],(w[0]+w[1])/2,(w[1]+w[2])/2)]
                            require([r[1] for r in tests]==[0,1,2],'Triggered count disagreement')
                            audits[key]=[r[2] for r in tests];self.a.write(file,audits)
                    neighbours.append(dict(pair=pair,chi_u=k,beta_deg=b,jump=jump,triggered=triggered,status='PASS'));previous=w
                    if b==0:zero=w[0]
                    rows.append(dict(pair=pair,chi_u=k,beta_deg=b,**{f'{name}{i}':r[name] for i,r in enumerate(rec['roots'],1) for name in ('omega','Lambda')},F_omega=w[0]/zero,gap12_relative=w[1]/w[0]-1,origin=origin,point_file=path.relative_to(ROOT).as_posix(),point_sha256=self.a.sha(path),shape_sha256=rec['shape_sha256']))
                    for i,r in enumerate(rec['roots'],1):
                        d=r['diagnostics'];quality.append(dict(pair=pair,chi_u=k,beta_deg=b,position=i,physical=float(max(np.max(d['joint_scaled']),np.max(d['external_scaled']))),ode=d['ode_error'],energy=d['energy_error'],singular=d['singular_ratio'],mass_gram=rec['mass_gram_error'],condition=d['nonzero_condition'],count_lower=r['counts'][0],count_upper=r['counts'][1],low_count=rec['count']['samples'][0]['count'],status='PASS'))
                self.csv('spectrum',rows);self.csv('quality',quality);self.csv('neighbours',neighbours)
                print(f'Map: {pair}, chi_u={k}: 51 accepted nodes',flush=True)
        self.csv('monotonicity',self.monotonic(points,self.cfg['map']['parameter_grid']))
        return rows,quality

    def analyze(self,rows):
        by={(r['pair'],r['chi_u'],r['beta_deg']):r for r in rows};summary=[];location=[];locstats=[];mechanics=[]
        for pair in self.cfg['cases']:
            for k in self.cfg['chi_values']:
                rr=[r for r in rows if (r['pair'],r['chi_u'])==(pair,k)];at={r['beta_deg']:r for r in rr}
                summary.append(dict(pair=pair,chi_u=k,Lambda0=at[0]['Lambda1'],Lambda3=at[3]['Lambda1'],Lambda10=at[10]['Lambda1'],Lambda20=at[20]['Lambda1'],gain_0_3_percent=100*(at[3]['F_omega']-1),gain_3_10_percent=100*(at[10]['F_omega']/at[3]['F_omega']-1),F10=at[10]['F_omega'],F20=at[20]['F_omega'],min_gap12=min(r['gap12_relative'] for r in rr),gap_min_beta=min(rr,key=lambda r:r['gap12_relative'])['beta_deg'],S0_per_radian=math.log(at[5]['F_omega'])/math.radians(5),beta_pl_status='not_established_no_tail'))
                zero=float(self.joint[pair,k,0]['omega1']);local=[]
                for b in sorted(b for p,ch,b in self.joint if (p,ch)==(pair,k) and (pair,k,b) in by):
                    row=self.joint[pair,k,b]
                    case=physics.base.Case(b,self.cfg['mu'],*self.cfg['cases'][pair]) if k=='RIGID' else slip.Case(b,self.cfg['mu'],*self.cfg['cases'][pair],k)
                    rec,path,_=self.archived(row,case);Fj=rec['roots'][0]['omega']/zero;Fe=by[pair,k,b]['F_omega']
                    local.append(dict(pair=pair,chi=k,beta_deg=b,F_end=Fe,F_joint=Fj,D_location=Fe/Fj-1,D_location_percent=100*(Fe/Fj-1),joint_point_file=path.relative_to(ROOT).as_posix(),joint_point_sha256=self.a.sha(path)))
                location+=local;worst=max(local,key=lambda r:abs(r['D_location']))
                locstats.append(dict(pair=pair,chi=k,common_nodes=len(local),min_D_percent=min(r['D_location_percent'] for r in local),max_D_percent=max(r['D_location_percent'] for r in local),max_abs_D_percent=abs(worst['D_location_percent']),worst_beta_deg=worst['beta_deg'],D3_percent=next(r['D_location_percent'] for r in local if r['beta_deg']==3),D10_percent=next(r['D_location_percent'] for r in local if r['beta_deg']==10),D20_percent=next(r['D_location_percent'] for r in local if r['beta_deg']==20)))
        m=self.run.model
        for r in rows:
            if r['beta_deg'] not in self.cfg['mechanics_angles']:continue
            path=ROOT/r['point_file'];self.check(path.with_suffix('.npz'),r['shape_sha256'])
            with np.load(path.with_suffix('.npz')) as data:arms=[{f:data[f'mode1_arm{i}_{f}'].copy() for f in ('x','weights','y','dy','ends')} for i in (1,2)]
            mass=physics.inner(m,dict(arms=arms),dict(arms=arms));ec=physics.energy_components(m,self.case(r['pair'],r['chi_u'],r['beta_deg']),arms)
            require(abs(mass-1)<1e-11 and abs(ec.sum()/(mass*r['omega1']**2)-1)<self.cfg['policy']['energy_tol'],'Mechanical energy check')
            u,w,_,N,Q,_=arms[0]['ends'][1];u20,N20=arms[1]['ends'][0,[0,3]]
            mechanics.append(dict(pair=r['pair'],chi_u=r['chi_u'],beta_deg=r['beta_deg'],mass=mass,V=float(ec.sum()/2),d_joint=math.hypot(u,w),u20=u20,N20=N20,**{name:float(e/ec.sum()) for name,e in zip(('axial_fraction','bend_fraction','shear_fraction','joint_rotation_fraction','end_rotation_fraction','axial_spring_fraction'),ec)}))
        for name,data in [('summary',summary),('location_comparison',location),('location_summary',locstats),('mechanics',mechanics)]:self.csv(name,data)
        return summary,locstats

    def leading_solution(self,pair,xi,order):
        ident=dict(pair=pair,xi_eff=xi,order=order,quadrature=self.cfg['leading']['quadrature'])
        key=hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16];path=self.folder/'leading_points'/(key+'.json')
        if path.exists():
            rec=json.loads(path.read_text());require(rec['identity']==ident and rec['status']=='PASS','Leading cache identity');solution=rec['solution']
        else:
            solver=Ritz(*self.cfg['cases'][pair],mu=self.cfg['mu'],order=order,quadrature=ident['quadrature']);solution=solver.solve(xi)
            self.leading_new+=1;d=solution['diagnostics']
            passed=all(d[k]<=self.cfg['leading']['quality'][k+'_tol'] for k in QUALITY)
            self.a.write(path,dict(identity=ident,solution=solution,status='PASS' if passed else 'FAIL'))
        d=solution['diagnostics'];require(all(d[k]<=self.cfg['leading']['quality'][k+'_tol'] for k in QUALITY),'Leading physical gate: '+str(ident)+' '+str(d))
        require(d['mass_cholesky'] and d['constraint_rank']==1 and d['P_rank']==1 and min(solution['zeta'])>0 and min(np.diff(solution['zeta']))>0,'Leading matrix/sorted prefix')
        return solution

    def leading_study(self,rows):
        cfg=self.cfg['leading'];by={(r['pair'],r['chi_u'],r['beta_deg']):r for r in rows}
        references=[r for r in read_csv(self.source('leading','curvature.csv')) if int(r['order'])==20]
        regression=[];comparison=[];convergence=[];curvature=[]
        for pair in self.cfg['cases']:
            zero=self.leading_solution(pair,0.,cfg['order']);old=next(r for r in references if r['ID']==pair)
            require(abs(zero['zeta'][0]/float(old['zeta0'])-1)<1e-10 and abs(zero['C_at_zero']/float(old['C'])-1)<1e-10,'Leading zero regression')
            regression.append(dict(pair=pair,zeta0=float(zero['zeta'][0]),C_rigid=zero['C_at_zero'],relative_zeta_error=zero['zeta'][0]/float(old['zeta0'])-1,relative_C_error=zero['C_at_zero']/float(old['C'])-1))
            for k in cfg['chi_values']:
                factor=physics.leading_scale(k);curvature.append(dict(pair=pair,chi_u=k,C=zero['C_at_zero']*factor**2))
                for b in cfg['beta_deg']:
                    xi=math.radians(b)/self.cfg['material']['epsilon'];require(xi<=30,'Outside addressed leading domain')
                    sol=self.leading_solution(pair,xi*factor,cfg['order'])
                    if b in cfg['check_beta_deg']:
                        high=self.leading_solution(pair,xi*factor,cfg['check_order']);error=float(max(abs(np.sqrt(np.array(sol['zeta'])/high['zeta'])-1)))
                        require(error<=cfg['quality']['frequency_convergence_rtol'],'Leading order convergence')
                        convergence.append(dict(pair=pair,chi_u=k,beta_deg=b,xi=xi,relative_frequency_error=error,status='PASS'))
                    Fe=by[pair,k,b]['F_omega'];Fl=math.sqrt(sol['zeta'][0]/zero['zeta'][0]);e=Fe/Fl-1
                    comparison.append(dict(pair=pair,chi_u=k,beta_deg=b,xi=xi,xi_eff=xi*factor,**{f'zeta{i}':float(z) for i,z in enumerate(sol['zeta'],1)},F_full=Fe,F_lim=Fl,e=e,delta_percent=100*abs(e),gap12_lim=sol['diagnostics']['gap12_relative']))
        quality=[]
        for path in sorted((self.folder/'leading_points').glob('*.json')):
            r=json.loads(path.read_text());quality.append(dict(**r['identity'],**{k:r['solution']['diagnostics'][k] for k in QUALITY},status=r['status']))
        stats=[]
        for pair in self.cfg['cases']:
            for k in cfg['chi_values']:
                rr=[r for r in comparison if (r['pair'],r['chi_u'])==(pair,k)];worst=max(rr,key=lambda r:r['delta_percent'])
                stats.append(dict(pair=pair,chi_u=k,nodes=len(rr),max_delta_percent=worst['delta_percent'],median_delta_percent=float(np.median([r['delta_percent'] for r in rr])),worst_beta_deg=worst['beta_deg'],worst_xi=worst['xi']))
        for name,data in [('leading_regression',regression),('leading_comparison',comparison),('leading_convergence',convergence),('leading_curvature',curvature),('leading_quality',quality),('leading_summary',stats)]:self.csv(name,data)
        return stats


def main(args):
    if args.stage in ('plot_only','finalize'):
        require(args.folder is not None,'Specify the accepted result folder')
        folder=args.folder.resolve();verify_numerical(folder)
        if args.stage=='plot_only':
            with plots_without_calculation():render(folder)
        else:finalize(folder)
        return
    exp=Experiment();print('Result: '+str(exp.folder),flush=True)
    try:
        exp.preflight()
        if args.stage=='preflight':exp.finish(args.stage);return
        rows,quality=exp.sweep();summary,location=exp.analyze(rows);leading=exp.leading_study(rows)
        exp.a.write(exp.folder/'numerical_summary.json',dict(status='PASS',map_points=len(rows),archived_map_points=sum(r['origin']=='archive' for r in rows),new_map_points=sum(r['origin']=='new_study' for r in rows),physical_maxima={k:max(r[k] for r in quality) for k in ('physical','ode','energy','singular','mass_gram','condition')},leading_solutions=len(list((exp.folder/'leading_points').glob('*.json'))),tail_points=0,branch_tracking=False))
        exp.finish(args.stage)
        numerical_manifest(exp.folder)
        bundle_readme(exp.folder)
    except Exception as exc:
        api().write(exp.folder/'diagnostic_failure.json',dict(stage=args.stage,error=repr(exc),status='STOP'))
        exp.finish(args.stage,'FAIL');raise


def numerical_manifest(folder):
    excluded={'numerical_manifest.json','manifest.json','executions.json','source_hashes.json','delivery_checks.json','test_results.json'}
    files={p.relative_to(folder).as_posix():api().sha(p) for p in folder.rglob('*') if p.is_file() and p.suffix in ('.json','.csv','.npz') and p.name not in excluded}
    api().write(folder/'numerical_manifest.json',dict(status='PASS',artifacts=files))


def verify_numerical(folder):
    manifest=json.loads((folder/'numerical_manifest.json').read_text());require(manifest['status']=='PASS','Unaccepted bundle')
    for name,digest in manifest['artifacts'].items():require(api().sha(folder/name)==digest,'Numerical artifact changed '+name)
    for name,digest in json.loads((folder/'source_hashes.json').read_text()).items():require(api().sha(ROOT/name)==digest,'Source changed '+name)


def render(folder):
    """Only table observations, no interpolation and no new physical fields."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,'legend.fontsize':8.5,
        'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':320})
    spectra=read_csv(folder/'spectrum.csv');location=read_csv(folder/'location_comparison.csv');leading=read_csv(folder/'leading_comparison.csv')
    pairs=['J0E4','J4E0'];ks=['RIGID',1.,.1,.01,0.]
    colors=['#222222','#176b8b','#9a5100','#8157a1','#007b52'];styles=['-','--','-.',':',(0,(5,1,1,1))]
    labels=['Точная опора',r'$\chi_u=1$',r'$\chi_u=0{,}1$',r'$\chi_u=0{,}01$',r'$\chi_u=0$']
    comma=FuncFormatter(lambda x,pos:f'{x:g}'.replace('.',','))
    def decorate(ax,letter,pair,xlabel=r'$\beta$, град.'):
        ax.xaxis.set_major_formatter(comma);ax.yaxis.set_major_formatter(comma)
        ax.grid(alpha=.22,linewidth=.5);ax.set_xlabel(xlabel)
        pairlabel=r'$(\kappa_{\mathrm{j}},\kappa_{\mathrm{e}})=(0{,}1;10)$' if pair=='J0E4' else r'$(\kappa_{\mathrm{j}},\kappa_{\mathrm{e}})=(10;0{,}1)$'
        ax.text(.02,1.035,letter+')  '+pairlabel,transform=ax.transAxes,fontsize=9)
    def save(fig,name):
        fig.savefig(folder/(name+'.pdf'),bbox_inches='tight',metadata={'Title':name,'Creator':'CoupledBeams, saved checked data'})
        fig.savefig(folder/(name+'.png'),dpi=320,bbox_inches='tight');plt.close(fig)
    for field,ylabel,name in [('Lambda1',r'$\Lambda_1$',FIGURES[0]),('F_omega',r'$F=\omega_1(\beta)/\omega_1(0)$',FIGURES[1])]:
        fig,axes=plt.subplots(1,2,figsize=(7.2,3.45),sharey=True,layout='constrained')
        for ax,letter,pair in zip(axes,'аб',pairs):
            for k,c,s,label in zip(ks,colors,styles,labels):
                rr=[r for r in spectra if r['pair']==pair and chi(r['chi_u'])==k]
                ax.plot([float(r['beta_deg']) for r in rr],[float(r[field]) for r in rr],color=c,linestyle=s,lw=1.6,label=label)
            ax.set_xlim(0,20);decorate(ax,letter,pair)
        handles,legend_labels=axes[0].get_legend_handles_labels()
        fig.legend(handles,legend_labels,loc='outside lower center',frameon=False,ncol=5,fontsize=9)
        axes[0].set_ylabel(ylabel);save(fig,name)
    fig,axes=plt.subplots(2,2,figsize=(7.2,5.8),sharex=True,layout='constrained')
    for row,k in enumerate([.01,0.]):
        for col,pair in enumerate(pairs):
            ax=axes[row,col];rr=[r for r in location if r['pair']==pair and chi(r['chi'])==k];x=[float(r['beta_deg']) for r in rr]
            ax.plot(x,[float(r['F_end']) for r in rr],color='#222222',lw=1.5,label='Пружина на внешней опоре')
            ax.plot(x,[float(r['F_joint']) for r in rr],color='#176b8b',ls='--',marker='o',markersize=3.5,markerfacecolor='white',lw=1.1,label='Пружина в узле')
            decorate(ax,'абвг'[2*row+col],pair);ax.set_xlim(0,20)
            ax.text(.04,.90 if k else .08,r'$\chi=0{,}01$' if k else r'$\chi=0$',transform=ax.transAxes)
            if row==0:ax.legend(loc='lower right',fontsize=7.7,frameon=False)
        axes[row,0].set_ylabel(r'$F=\omega_1(\beta)/\omega_1(0)$')
    save(fig,FIGURES[2])
    fig,axes=plt.subplots(2,2,figsize=(7.2,5.6),sharex=True,sharey='row',layout='constrained')
    for col,pair in enumerate(pairs):
        ax=axes[0,col]
        for k,marker,i in [('RIGID','o',0),(.01,'s',3),(0.,'^',4)]:
            rr=[r for r in leading if r['pair']==pair and chi(r['chi_u'])==k];x=[float(r['xi']) for r in rr]
            ax.plot(x,[float(r['F_lim']) for r in rr],color=colors[i],ls=styles[i],lw=1.5,label=labels[i])
            ax.plot(x,[float(r['F_full']) for r in rr],color=colors[i],ls='none',marker=marker,ms=4,mfc='white',mew=.9)
            axes[1,col].plot(x,[float(r['delta_percent']) for r in rr],color=colors[i],ls=styles[i],lw=1.3,marker=marker,ms=3,mfc='white')
        for row in (0,1):decorate(axes[row,col],'абвг'[row*2+col],pair,r'$\xi=\beta_{\mathrm{rad}}/\varepsilon$');axes[row,col].set_xlim(0,30)
        ax.legend(loc='upper left',frameon=False)
    axes[0,0].set_ylabel(r'$F=\omega_1(\beta)/\omega_1(0)$');axes[1,0].set_ylabel(r'$\delta$, %');save(fig,FIGURES[3])


def finalize(folder):
    verify_numerical(folder)
    require(all((folder/(n+ext)).exists() for n in FIGURES for ext in ('.pdf','.png')),'Missing figure')
    report=ROOT/'docs/theory/isotropic_two_springs_external_axial_spring.md'
    files={p.relative_to(folder).as_posix():api().sha(p) for p in folder.rglob('*') if p.is_file() and p.name!='manifest.json'}
    api().write(folder/'manifest.json',dict(status='PASS',contract_sha256=api().sha(folder/'contract.json'),artifacts=files,
        scientific_report=report.relative_to(ROOT).as_posix(),scientific_report_sha256=api().sha(report),
        physical_implementation_sha256=api().sha(physics.__file__),workflow_sha256=api().sha(__file__)))


def bundle_readme(folder):
    relative=folder.relative_to(ROOT).as_posix()
    (folder/'README.md').write_text(f'''# Продольная податливость внешней опоры

Принятая ограниченная серия `{folder.name}`, 2026-10-25.
[Научный отчёт](../../../../docs/theory/isotropic_two_springs_external_axial_spring.md).
Модель RL_FULL/Тимошенко, две точные поступательные связи узла и две
поворотные пружины. Новая продольная пружина действует вдоль локальной оси
второго стержня на внешнем конце. mu=0,01; epsilon=0,005; b/h=4;
nu=0,3; ks=5/6; l=E=rho=1. Пары J0E4=(0,1;10), J4E0=(10;0,1).
chi_u=RIGID;1;0,1;0,01;0. RIGID является точным ограничением.

## Источники и воспроизведение

- Точная опора: `length_asymmetry/ddf22cb06e8563e2/complete_spectrum.csv`,
  mu=0,01, 102 проверенные точки основной сетки. CSV ссылается на исходные
  point JSON/NPZ более ранних серий; эти пути и хеши сохранены в `spectrum.csv`.
- Пружина в узле: `longitudinal_slip/253c9655120041e1/spectrum.csv`;
  только совпадающие углы, 138 сравнений. Узловых расчётов нет.
- Ведущая регрессия: `leading_validation/f1c670fe62102196/curvature.csv`,
  нулевые собственные пары и C порядка 20.
- Новая внешняя модель: `points/*.json` и `.npz`; 408 точек карты и десять
  проверок RIGID. Ритц полной задачи: 20 адресных проверок.
- Ведущая задача: `leading_points/*.json`, 38 решений порядка 20 и десять
  проверок порядка 24. Повторяющиеся аргументы, включая нуль, переиспользованы.

Все пути выше относятся к `results/isotropic_two_springs_rl_timoshenko/`.
`contract.json` содержит входной контракт, закреплённые источники и SHA256
нового физического модуля; идентификатор — первые 16 знаков SHA256 этого
упорядоченного JSON. `source_hashes.json` фиксирует прочитанные источники.
`numerical_manifest.json` проверяет численные артефакты, `manifest.json` —
итоговый комплект. Прежние файлы не перезаписываются.

Из корня проекта, в Python с NumPy, SciPy, Matplotlib и pytest:

```powershell
python scripts/analysis/study_two_springs_external_axial_spring.py --stage preflight
python scripts/analysis/study_two_springs_external_axial_spring.py --stage run
python scripts/analysis/study_two_springs_external_axial_spring.py --stage plot_only --folder {relative}
python -m pytest tests/test_two_springs_external_axial_spring.py -q
python scripts/analysis/study_two_springs_external_axial_spring.py --stage finalize --folder {relative}
```

В текущем окружении использован `D:/python/Pycharm/pythonProject/.venv/Scripts/python.exe -X utf8`.
Повторный `run` читает принятый кэш; новые решения учитываются в `executions.json`.
`plot_only` запрещает вызовы спектральных функций и работает только с таблицами.
При нарушении физических допусков расчёт останавливается с диагностикой.
Ни одна стадия не редактирует научный отчёт автоматически.

## Таблицы, единицы и нормировки

CSV записаны с десятичной точкой; рисунки используют десятичные запятые.
`beta_deg` измеряется в градусах, `xi=beta_rad/epsilon`, `xi_eff` — аргумент
прежнего ведущего решателя, а не физический угол. Порядок частот отсортированный.

| Файл | Содержание |
| --- | --- |
| `spectrum.csv` | 510 точек; omega1–3 в единицах l=E=rho=1; Lambda1–3 безразмерны; Lambda^4=rho*A*omega^2*l^4/(EI); F_omega=omega1/собственная omega1(0); gap12_relative=(omega2-omega1)/omega1; источник и хеши |
| `summary.csv` | Абсолютные Lambda, прирост 0–3° относительно нуля и 3–10° относительно 3°, F10/F20, минимум зазора; S0_per_radian=ln(F(5°))/(5*pi/180); beta_pl не установлен |
| `location_comparison.csv` | F_end/F_joint с собственной нормировкой каждой модели; D_location=F_end/F_joint-1, D_location_percent=100*D; точный источник узловой точки |
| `location_summary.csv` | Подписанные пределы и максимумы модуля только на общей конечной сетке, без интерполяции |
| `mechanics.csv` | Первая форма с полной M=1; V=omega1^2/2, d_joint=sqrt(u1^2+w1^2), u20 и N20; доли шести положительных энергетических членов. Знак формы произволен |
| `leading_comparison.csv` | Действительные beta и xi; zeta1–3=Lambda_lim^4; F_lim=sqrt(zeta1/zeta1(0)); e=F_full/F_lim-1; delta_percent=100*abs(e) |
| `leading_summary.csv` | Максимум и медиана delta по десяти узлам каждой кривой |
| `leading_curvature.csv` | C в ln(F_lim)=C*xi^2+o(xi^2), из прежней нулевой формы и точного коэффициента податливости |
| `quality.csv`, `neighbours.csv`, `triggered_counts.json` | Физические невязки и полный массовый Gram; счёт 0→1→2→3; четыре адресные проверки соседей |
| `regression.csv`, `ritz_checks.csv`, `preflight*.json/csv` | Предварительные проверки, точный жёсткий предел, независимый прямой Ритц и отсутствие механизма |
| `monotonicity.csv` | Проверки неубывания трёх абсолютных частот по chi_u, включая RIGID; 408 строк по три позиции |
| `leading_quality.csv`, `leading_convergence.csv`, `leading_regression.csv` | Сильные/естественные невязки всех трёх позиций, положительная масса, точные ограничения, энергетическое равенство, порядки 20/24 и регрессия |
| `compliance_identity.csv` | Независимая точная рациональная проверка минимума энергии внешней и узловой моделей |

## Рисунки и предлагаемые подписи

Каждый рисунок сохранён в векторном PDF и PNG 320 dpi. Параметры общие,
как указано выше. Линии только соединяют фактические точки; сглаживания нет.

1. `01_first_frequency`: первая отсортированная частота Lambda1 при изменении
   угла и продольной жёсткости внешней опоры. Панели а, б — две пары поворотных
   жёсткостей. Источник: `spectrum.csv`, без нормировки на начальную частоту.
2. `02_normalized_growth`: относительное изменение первой частоты F при
   податливой внешней опоре. Каждая зависимость нормирована на собственное
   значение при beta=0. Источник: `spectrum.csv`.
3. `03_spring_location`: сравнение внешней и узловой продольных пружин при
   одинаковой chi=0,01 (а, б) и chi=0 (в, г). При chi=0 обе поворотные пружины
   и вторая поступательная связь узла сохраняются. Показаны только общие
   расчётные узлы. Источник: `location_comparison.csv`.
4. `04_leading_comparison`: нормированные частоты полной модели (маркеры)
   и ведущей задачи (линии), а, б; абсолютное относительное расхождение
   delta=100*abs(F_full/F_lim-1), в, г. Источник: `leading_comparison.csv`;
   xi<=29,670597. Равенство ведущих задач двух положений пружины следует
   из энергетического вывода, а приведённые расхождения являются наблюдениями
   на выбранной численной сетке при одной epsilon.

Хвосты до 90°, характерные углы, модальное продолжение, FEM и новые
геометрии не рассчитывались. Результат не доказывает сходимость собственных
пар и не устанавливает равенство полных динамических задач.
''',encoding='utf8')
