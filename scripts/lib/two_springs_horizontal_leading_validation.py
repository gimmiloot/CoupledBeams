"""HS-L1 bounded leading validation and four full-model threshold searches.

Script Proliferation Control: new leading-space/threshold contract, not a
pilot preset or a general map runner. Reuses integrated fields, frozen full
Run/prefix/count/physical checks and HS-P1 records. Old source identities stay
unchanged. Stable module entry: run, verify, plot_only. No global angular map.
"""
from pathlib import Path
from dataclasses import asdict
import argparse
import csv
import hashlib
import json
import math
import platform
import subprocess
import time
from fractions import Fraction as F
import numpy as np
import scipy
from scipy.optimize import brentq
from scripts.lib import isotropic_two_springs_horizontal_leading as lead
from scripts.lib.two_springs_horizontal_pilot import validate_record, key
from scripts.analysis import verify_two_springs_horizontal_support as full

ROOT=Path(__file__).resolve().parents[2]
CONFIG=ROOT/'data/input/isotropic_two_springs_horizontal_leading.json'
BASE=ROOT/'results/isotropic_two_springs_rl_timoshenko/horizontal_spring_leading'
sha,write,csvwrite=full.sha,full.write,full.csvwrite


def require(ok,message):
    if not ok:raise ArithmeticError(message)


def identity(cfg):
    pilot=ROOT/cfg['pilot_archive']
    record=dict(config=cfg,source_manifest=sha(pilot/'manifest.json'),
        source_analysis=sha(pilot/'analysis_manifest.json'),initial_state=sha(ROOT/cfg['initial_workspace']),
        files={p.relative_to(ROOT).as_posix():sha(p) for p in [CONFIG,Path(__file__),Path(lead.__file__)]},
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    return hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()[:16],record


def mathematical_check():
    # Independent exact rational minimizer: A=U1(L1), B=U2(0), C=U2(L2).
    # A+C=-xi*W2(L2); r=B+xi*W2(0); S=A+C-B=xi*D-r.
    for chi in [F(1,2),F(1),F(3,7)]:
        xi,D=F(7,3),F(-4,5);r=xi*D/(1+2*chi);S=xi*D-r
        assert S*S/2+chi*r*r==chi/(1+2*chi)*xi*xi*D*D
        assert -S+2*chi*r==0
    return dict(status='PASS',g='chi_s/(1+2*chi_s)',D='W2(0)-W2(lambda2)',
        saturation='D=0; W1third(L1)-W2third(L2)+W2third(0)=0',
        independent_condition_count_finite=8,independent_condition_count_saturated=8,
        formal_status='formal controlled scaling on bounded xi; no spectral-convergence theorem claimed',
        checking='direct variational derivation, exact rational elimination, independent 8-row ODE; Lean unavailable')


class Study:
    def __init__(self,cfg):
        self.cfg=cfg;self.identifier,self.ident=identity(cfg);self.folder=BASE/self.identifier
        self.folder.mkdir(parents=True,exist_ok=True);self.lp=cfg['leading'];self.eps=cfg['material']['epsilon']
        self.models={p:lead.Leading(cfg['k_joint'],cfg['mu'],p,self.lp['quadrature']) for p in self.lp['orders']}
        self.leading_records={};self.leading_solutions={};self.full_records={};self.used_archive={};self.comparisons={}
        self.convergence=[];self.full_new=0;self.recovered_full=set();self.recovered_leading=0;self.full_diagnostics=[];self.events=[]
        physical=json.loads((ROOT/'data/input/isotropic_two_springs_horizontal_spring.json').read_text())
        require(physical['policy']==cfg['policy_full'],'Full-model tolerances changed')
        self.runner=full.Run(physical,self.folder/'full');self.start=time.perf_counter()

    def audit(self):
        for name,digest in self.cfg['frozen_files'].items():require(sha(ROOT/name)==digest,'Frozen source '+name)
        pilot=ROOT/self.cfg['pilot_archive']
        for name in ['manifest.json','analysis_manifest.json']:
            for artifact,digest in json.loads((pilot/name).read_text())['artifacts'].items():
                require(sha(pilot/artifact)==digest,'Source artifact '+artifact)
        self.source_rows=list(csv.DictReader(open(pilot/'spectra.csv',encoding='utf8')))
        self.source={}
        for row in self.source_rows:
            if row['id'] not in ['H5','H6']:continue
            p=ROOT/row['point_file'];record=json.loads(p.read_text());case=full.h.Case(**record['case'])
            validate_record(record,case,self.cfg['policy_full']);self.source[key(case)]=(record,p)
        self.recovery={}
        if self.cfg.get('recovery_archive'):
            previous=ROOT/self.cfg['recovery_archive']
            for name,digest in self.cfg['recovery_sha256'].items():require(sha(previous/name)==digest,'Recovery artifact '+name)
            for p in (previous/'full/points').glob('*.json'):
                record=json.loads(p.read_text());case=full.h.Case(**record['case'])
                validate_record(record,case,self.cfg['policy_full']);self.recovery[key(case)]=(record,p)
            if not (self.folder/'manifest.json').exists():
                for p in (previous/'leading').glob('*.json'):
                    record=json.loads(p.read_text());tau=record['tau'];k=lead.RIGID if tau==lead.RIGID else float(tau).hex()
                    high=record['orders'][-1]
                    with np.load(p.with_suffix('.npz')) as z:coef=z['p24_coefficients'].copy()
                    self.leading_records[k]=record
                    self.leading_solutions[k]=dict(zeta=np.array(high['zeta']),coefficients=coef,diagnostics=high['diagnostics'],per_mode=high['per_mode'])
                    for source in [p,p.with_suffix('.npz')]:
                        target=self.folder/'leading'/source.name;target.parent.mkdir(exist_ok=True);target.write_bytes(source.read_bytes())
                    self.recovered_leading+=1
        write(self.folder/'source_audit.json',dict(status='PASS',pilot=self.cfg['pilot_archive'],available_points=len(self.source),HEAD=subprocess.run(['git','rev-parse','HEAD'],capture_output=True,text=True,check=True).stdout.strip()))

    def leading(self,tau,independent=False):
        k=lead.RIGID if tau==lead.RIGID else float(tau).hex()
        if k not in self.leading_records:
            solutions=[];records=[];arrays={}
            for p,m in self.models.items():
                s=m.solve(tau);d=s['diagnostics']
                row=dict(order=p,tau=tau,zeta=s['zeta'].tolist(),diagnostics=d,per_mode=s['per_mode'])
                records.append(row);solutions.append(s);arrays[f'p{p}_coefficients']=s['coefficients']
                for metric,tol in [('strong_ode','strong_tol'),('natural_boundary','boundary_tol'),('essential','essential_tol'),('matrix_residual','matrix_tol'),('energy','energy_tol'),('mass_gram','mass_tol'),('symmetry','matrix_tol')]:
                    if d[metric]>self.lp[tol]:
                        write(self.folder/'leading_failure.json',row);raise ArithmeticError(f'Leading gate {tau} p{p}: {metric}={d[metric]}')
                require(np.all(s['zeta']>0) and np.all(np.diff(s['zeta'])>0),'Leading ordering')
            low,high=solutions
            frequency_error=float(np.max(abs(np.sqrt(low['zeta']/high['zeta'])-1)))
            overlap=[]
            for i in range(3):
                total=0.
                for arm,L in enumerate(self.models[self.lp['orders'][0]].lengths):
                    m0,m1=list(self.models.values());a=m0.field_matrices(arm,m1.nodes)[0]@low['coefficients'][:,i];b=m1.values[arm][0]@high['coefficients'][:,i]
                    total+=float(m1.weights@(a*b)*L)
                overlap.append(abs(total))
            shape_error=float(np.max(abs(np.array(overlap)-1)))
            require(frequency_error<=self.lp['frequency_convergence_tol'],'p20/24 frequencies')
            require(shape_error<=self.lp['shape_convergence_tol'],'p20/24 shapes')
            record=dict(tau=tau,status='PASS',orders=records,frequency_convergence=frequency_error,shape_convergence=shape_error)
            self.leading_records[k]=record;self.leading_solutions[k]=high
            name=hashlib.sha256(k.encode()).hexdigest()[:16]
            write(self.folder/'leading'/f'{name}.json',record)
            np.savez_compressed(self.folder/'leading'/f'{name}.npz',**arrays)
        record=self.leading_records[k];s=self.leading_solutions[k]
        if independent and 'independent_ODE' not in record:
            checks=lead.independent_roots(s['zeta'],tau,self.cfg['k_joint'],self.cfg['mu'])
            require(max(abs(a['relative_frequency_error']) for a in checks)<=self.lp['ode_independent_frequency_tol'],'Independent leading ODE')
            record['independent_ODE']=checks
            write(self.folder/'leading'/f'{hashlib.sha256(k.encode()).hexdigest()[:16]}.json',record)
        return s

    def full_point(self,case_id,beta,label,ritz=False):
        construction=next(c for c in self.cfg['constructions'] if c['id']==case_id)
        case=full.h.Case(float(beta),self.cfg['mu'],self.cfg['k_joint'],construction['kappa_horizontal']);k=key(case)
        if k not in self.full_records:
            if k in self.source:
                rec,path=self.source[k];self.used_archive[k]=path;origin='HS-P1'
            elif k in self.recovery:
                rec,source=self.recovery[k];self.recovered_full.add(k);origin='HS-L1_recovered_after_table_error'
                path=self.folder/'full/points'/source.name;path.parent.mkdir(parents=True,exist_ok=True)
                for a,b in [(source,path),(source.with_suffix('.npz'),path.with_suffix('.npz'))]:b.write_bytes(a.read_bytes())
            else:
                rec,_=self.runner.solve(case,label,ritz=ritz);self.full_new+=1
                path=self.folder/'full/points'/f'{self.runner.name(case)}.json';origin='new'
            validate_record(rec,case,self.cfg['policy_full'])
            self.full_records[k]=(rec,path,origin)
        return self.full_records[k][0]

    def compare(self,case_id,beta,label):
        c=next(c for c in self.cfg['constructions'] if c['id']==case_id)
        fullrec=self.full_point(case_id,beta,label)
        xi=math.radians(beta)/self.eps;tau=c['g']*xi*xi;s=self.leading(tau)
        row=dict(id=case_id,beta_deg=float(beta),xi=xi,tau=tau,label=label,omega_full=fullrec['omega'][0],omega_lead=self.eps*math.sqrt(s['zeta'][0]))
        row['relative_frequency_error']=row['omega_lead']/row['omega_full']-1
        self.comparisons[(case_id,float(beta))]=row;return row

    def run(self):
        if (self.folder/'manifest.json').exists():self.verify();return
        write(self.folder/'contract.json',self.ident);self.audit();write(self.folder/'mathematical_gate.json',mathematical_check())
        try:
            for tau in self.lp['preflight_tau']:self.leading(tau,independent=True)
            saturated=self.leading(lead.RIGID);base=self.leading(0)
            omega_inf=self.eps*math.sqrt(saturated['zeta'][0]);omega0=self.eps*math.sqrt(base['zeta'][0])
            for tau in self.lp['saturation_tau']:self.leading(tau)
            finite=sorted((float(r['tau']),self.leading_solutions[k]['zeta']) for k,r in self.leading_records.items() if r['tau']!=lead.RIGID)
            require(all(np.all(b[1]>=a[1]*(1-1e-10)) for a,b in zip(finite,finite[1:])),'Tau monotonicity')
            require(all(np.all(z<=saturated['zeta']*(1+1e-10)) for t,z in finite),'Saturation upper bound')
            coeff=list(csv.DictReader(open(ROOT/self.cfg['pilot_archive']/'small_angle_coefficients.csv')))
            initial=[]
            for c in self.cfg['constructions']:
                old=self.full_point(c['id'],0,'initial');tail=self.full_point(c['id'],90,'reference level')
                oldC=float(next(a for a in coeff if a['id']==c['id'])['coefficient_beta_rad_squared'])
                C=c['g']*base['per_mode'][0]['D']**2
                initial.append(dict(id=c['id'],omega0_lead=omega0,omega0_full=old['omega'][0],omega0_relative_error=omega0/old['omega'][0]-1,omega_inf_lead=omega_inf,omega90_full=tail['omega'][0],saturation_relative_error=omega_inf/tail['omega'][0]-1,C2_lead=C,C2_full=oldC,C2_relative_error=C/oldC-1))
            require(initial[0]['C2_lead']/initial[1]['C2_lead']==.5,'Initial coefficient ratio')
            csvwrite(self.folder/'initial_and_limit.csv',initial)
            threshold_tau={};thresholds=[]
            for q in self.lp['threshold_levels']:
                fun=lambda t:math.sqrt(self.leading(t)['zeta'][0]/saturated['zeta'][0])-q
                lo,hi=self.lp['tau_initial_bracket'];require(fun(lo)<0<fun(hi),'Leading threshold bracket')
                t=brentq(fun,lo,hi,xtol=self.lp['tau_root_xtol'],rtol=self.lp['tau_root_rtol'])
                self.leading(t,independent=True);threshold_tau[q]=t
                write(self.folder/'leading_thresholds.json',{str(level):dict(tau=value,F_residual=float(math.sqrt(self.leading(value)['zeta'][0]/saturated['zeta'][0])-level)) for level,value in threshold_tau.items()})
            # Prior forecasts are now fixed, before any new full-model solves.
            for c in self.cfg['constructions']:
                case_id=c['id'];ref=self.full_point(case_id,90,'reference level')['omega'][0]
                available=sorted((float(r['beta_deg']),float(r['omega1'])) for r in self.source_rows if r['id']==case_id)
                for q,t in threshold_tau.items():
                    pred=math.degrees(self.eps*math.sqrt(t/c['g']))
                    crossings=[(a[0],b[0]) for a,b in zip(available,available[1:]) if a[1]<q*ref<=b[1]]
                    require(len(crossings)==1,'Full pilot threshold not uniquely bracketed')
                    a,b=crossings[0]
                    fun=lambda angle:self.full_point(case_id,angle,f'{case_id} q{q} local root')['omega'][0]/ref-q
                    root=brentq(fun,a,b,xtol=self.cfg['full_threshold']['angle_xtol_deg'],rtol=self.cfg['full_threshold']['angle_rtol'],maxiter=self.cfg['full_threshold']['max_iterations'])
                    h=self.cfg['full_threshold']['bracket_halfwidth_deg'];lower,upper=max(a,root-h),min(b,root+h)
                    fl,fu=fun(lower),fun(upper)
                    require(fl<0<fu,'Full threshold signs unresolved')
                    require(min(abs(fl),abs(fu))>2e-10,'Full threshold sign margin unresolved')
                    require(abs(fun(root))<=self.cfg['full_threshold']['frequency_equation_tol'],'Full threshold equation residual')
                    forecast=self.compare(case_id,pred,'apriori forecast')
                    actual=self.compare(case_id,root,'full threshold')
                    # Diagnostic level replacement is separate from the a priori prediction.
                    target=q*ref/self.eps
                    diagnostic=None
                    if target<math.sqrt(saturated['zeta'][0]):
                        f=lambda t:math.sqrt(self.leading(t)['zeta'][0])-target
                        td=brentq(f,0,100,xtol=self.lp['tau_root_xtol'],rtol=self.lp['tau_root_rtol'])
                        diagnostic=math.degrees(self.eps*math.sqrt(td/c['g']))
                    row=dict(id=case_id,q=q,tau_lead=t,beta_lead_deg=pred,beta_full_deg=root,beta_full_bracket_lower=lower,beta_full_bracket_upper=upper,bracket_F_lower=fl,bracket_F_upper=fu,pilot_bracket_lower=a,pilot_bracket_upper=b,xi_lead=math.radians(pred)/self.eps,xi_full=math.radians(root)/self.eps,angle_signed_error_deg=pred-root,angle_absolute_error_deg=abs(pred-root),angle_relative_error=pred/root-1,frequency_relative_error_at_forecast=forecast['relative_frequency_error'],omega_full_at_forecast=forecast['omega_full'],omega_lead_at_forecast=forecast['omega_lead'],beta_diagnostic_full_level_deg=diagnostic,equation_residual=fun(root))
                    thresholds.append(row);csvwrite(self.folder/'thresholds.csv',thresholds)
                    for factor in self.cfg['full_threshold']['near_threshold_factors']:self.compare(case_id,root*factor,'threshold neighbourhood')
            ratios=[]
            for q in threshold_tau:
                a,b=[next(r for r in thresholds if r['id']==name and r['q']==q) for name in ['H5','H6']]
                ratios.append(dict(q=q,leading_ratio=a['beta_lead_deg']/b['beta_lead_deg'],full_ratio=a['beta_full_deg']/b['beta_full_deg']))
            csvwrite(self.folder/'angle_ratios.csv',ratios)
            for row in self.source_rows:
                if row['id'] in ['H5','H6'] and float(row['beta_deg'])<=self.cfg['full_threshold']['pilot_compare_max_angle']:
                    self.compare(row['id'],float(row['beta_deg']),'archived pilot')
            csvwrite(self.folder/'comparison.csv',sorted(self.comparisons.values(),key=lambda r:(r['id'],r['beta_deg'])))
            records=[]
            for k,(r,path,origin) in self.full_records.items():
                records.append(dict(id='H6' if k[2]=='RIGID' else 'H5',beta_deg=k[0],origin=origin,path=path.relative_to(ROOT).as_posix(),omega_1=r['omega'][0],omega_2=r['omega'][1],omega_3=r['omega'][2],status=r['status'],boundary=max(max(np.max(d['joint_scaled']),np.max(d['external_scaled'])) for d in r['diagnostics']),ode=max(d['ode_error'] for d in r['diagnostics']),energy=max(d['energy_error'] for d in r['diagnostics']),mass_gram=r['mass_gram_error']))
            csvwrite(self.folder/'full_points.csv',records)
            neighbour=[]
            for case_id in ['H5','H6']:
                local=sorted((r for r in records if r['id']==case_id and r['beta_deg']<=5),key=lambda r:r['beta_deg'])
                for a,b in zip(local,local[1:]):
                    jump=max(abs(b[f'omega_{i}']-a[f'omega_{i}'])/max(b[f'omega_{i}'],a[f'omega_{i}']) for i in [1,2])
                    row=dict(id=case_id,left=a['beta_deg'],right=b['beta_deg'],jump=jump,status='PASS',triggered=jump>.35)
                    if jump>.35:
                        case=full.h.Case(b['beta_deg'],self.cfg['mu'],self.cfg['k_joint'],20000 if case_id=='H5' else 'RIGID')
                        w,_,d=full.h.ritz(self.runner.m,case,20,rank_tol=self.cfg['policy_full']['kernel_rank_relative_tol'])
                        error=float(np.max(abs(w/np.array([b[f'omega_{i}'] for i in [1,2,3]])-1)))
                        self.runner.ritz_solves+=1;row.update(ritz_order=20,ritz_frequency_error=error,ritz_diagnostics=d)
                        require(error<=self.cfg['policy_full']['ritz_frequency_tol'],'Neighbour Ritz frequency')
                        require(d['energy_error']<=self.cfg['policy_full']['energy_tol'],'Neighbour Ritz energy')
                        require(d['mass_gram_error']<=self.cfg['policy_full']['mass_gram_tol'],'Neighbour Ritz mass')
                    neighbour.append(row)
            write(self.folder/'neighbour_audit.json',neighbour)
            write(self.folder/'leading_records.json',list(self.leading_records.values()))
            write(self.folder/'failures.json',[])
            output=dict(status='PASS',identity=self.ident,leading_distinct_tau=len(self.leading_records),leading_Ritz_solves=2*len(self.leading_records),leading_recovered_points=self.recovered_leading,full_new_points=self.full_new+len(self.recovered_full),full_new_this_invocation=self.full_new,full_recovered_points=len(self.recovered_full),full_reused_points=len(self.used_archive),full_ritz_solves=self.runner.ritz_solves,completion_seconds=time.perf_counter()-self.start,artifacts={p.relative_to(self.folder).as_posix():sha(p) for p in self.folder.rglob('*') if p.is_file()})
            write(self.folder/'manifest.json',output)
        except Exception as error:
            write(self.folder/'failure.json',dict(status='PARTIAL_PASS',error=repr(error),full_new=self.full_new,leading_accepted=len(self.leading_records)));raise

    def verify(self):
        self.audit();m=json.loads((self.folder/'manifest.json').read_text())
        require(m['identity']==self.ident,'Identity mismatch')
        for file,digest in m['artifacts'].items():require(sha(self.folder/file)==digest,'Archive artifact '+file)


def plot_only(folder):
    # CSV-only renderer: no model construction, matrix assembly or root search.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=list(csv.DictReader(open(folder/'comparison.csv',encoding='utf8')))
    levels=list(csv.DictReader(open(folder/'thresholds.csv',encoding='utf8')))
    matplotlib.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'pdf.fonttype':42})
    output=folder/'figures';output.mkdir(exist_ok=True)
    for kind in ['frequency','error']:
        fig,axes=plt.subplots(1,2,figsize=(10,3.8),layout='constrained')
        for ax,name in zip(axes,['H5','H6']):
            data=sorted([r for r in rows if r['id']==name],key=lambda r:float(r['beta_deg']))
            x=[float(r['beta_deg']) for r in data]
            if kind=='frequency':
                ax.plot(x,[float(r['omega_lead']) for r in data],label='Ведущая модель',color='#0072B2')
                ax.plot(x,[float(r['omega_full']) for r in data],'o',ms=3,label='Полная модель',color='#D55E00')
                ax.set_ylabel('Первая круговая частота')
            else:
                ax.plot(x,[100*float(r['relative_frequency_error']) for r in data],'o-',ms=3,color='#0072B2')
                ax.axhline(0,color='black',lw=.7);ax.set_ylabel('Расхождение частоты, %')
            for r in levels:
                if r['id']==name:ax.axvline(float(r['beta_full_deg']),ls=':',color='gray',lw=1)
            ax.set_title(name);ax.set_xlabel('Угол сопряжения, град');ax.grid(alpha=.2)
        if kind=='frequency':axes[0].legend()
        stem='01_frequency' if kind=='frequency' else '02_error'
        fig.savefig(output/(stem+'.pdf'));fig.savefig(output/(stem+'.png'),dpi=320);plt.close(fig)
        (output/(stem+'.csv')).write_bytes((folder/'comparison.csv').read_bytes())
    write(folder/'figures_manifest.json',dict(calculation_mode='plot_only',new_eigenvalue_solves=0,dpi=320,artifacts={p.name:sha(p) for p in output.iterdir()}))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['run','verify','plot_only'],default='run');args=parser.parse_args()
    cfg=json.loads(CONFIG.read_text());identifier,_=identity(cfg);folder=BASE/identifier
    if args.stage=='plot_only':plot_only(folder)
    else:
        study=Study(cfg)
        if args.stage=='run':study.run()
        else:study.verify()
    print(args.stage,folder)


if __name__=='__main__':main()
