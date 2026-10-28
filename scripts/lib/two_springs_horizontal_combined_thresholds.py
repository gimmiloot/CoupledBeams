"""HS-LA6: close the threshold study using two accepted limit problems.

This bounded archival-comparison contract seals a priori forecasts before
parsing full-model references. It reuses frozen HS-L1 Ritz/ODE and HS-LA4
prefix checks; it neither imports a full-model runner nor adds physics.
Separate entry preserves all accepted source-hashed workflows.
"""
import csv
import hashlib
import json
import math
import platform
import numpy as np
import scipy
from scipy.optimize import brentq
from scripts.lib.isotropic_two_springs_horizontal_leading import Leading
from scripts.lib.two_springs_horizontal_joint_stiffness import check_prefix
from scripts.lib.two_springs_horizontal_spectral_structure import ROOT,sha,write,csvwrite,require

CONFIG=ROOT/'data/input/isotropic_two_springs_horizontal_combined_thresholds.json'
BASE=ROOT/'results/isotropic_two_springs_rl_timoshenko/horizontal_combined_threshold_prediction'
ARCHIVES=ROOT/'results/isotropic_two_springs_rl_timoshenko'


def identity(cfg):
    record=dict(config=cfg,source_sha256=sha(__file__),config_sha256=sha(CONFIG),
        initial_sha256=sha(ROOT/cfg['initial_workspace']),
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    return hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()[:16],record


def read_csv(path):
    with open(path,encoding='utf8') as stream:return list(csv.DictReader(stream))


class Study:
    def __init__(self,cfg):
        self.cfg=cfg;identifier,self.identity=identity(cfg);self.folder=BASE/identifier
        self.models={p:Leading(cfg['k_joint'],cfg['mu'],p,cfg['quadrature']) for p in cfg['orders']}
        self.records=[];self.solutions={};self.index={};self.new=0;self.reused=0;self.counts=[];self.events=[]

    def audit(self):
        cfg=self.cfg;state=json.loads((ROOT/cfg['initial_workspace']).read_text())
        for name,digest in cfg['frozen_files'].items():require(sha(ROOT/name)==digest,'Frozen source '+name)
        for name,digest in state['archive_files'].items():require(sha(ROOT/name)==digest,'Historical archive '+name)
        for name in cfg['archives'].values():
            manifest=json.loads((ARCHIVES/name/'manifest.json').read_text())
            require(manifest['status']=='PASS','Unaccepted source '+name)
        outer=ARCHIVES/cfg['archives']['outer']
        outer_contract=json.loads((outer/'contract.json').read_text())['config']
        require(outer_contract['mu']==cfg['mu'] and outer_contract['epsilon']==cfg['epsilon'],'Outer normalization')
        rows=json.loads((outer/'point_records.json').read_text())
        self.outer={p:next(r for r in rows if r['order']==p and r['k_joint']==cfg['k_joint']) for p in cfg['orders']}
        self.zout=self.outer[24]['zeta'][0]  # One authoritative target for both Ritz orders.
        require(abs(math.sqrt(self.outer[20]['zeta'][0]/self.zout)-1)<cfg['policy']['frequency_convergence_tol'],'Outer convergence')
        coef=ARCHIVES/cfg['archives']['coefficients']/'coefficients.json'
        self.coefficients={r['order']:r for r in json.loads(coef.read_text())}
        self.original_config=json.loads((ROOT/'data/input/isotropic_two_springs_horizontal_leading.json').read_text())
        require(self.original_config['mu']==cfg['mu'] and self.original_config['k_joint']==cfg['k_joint'],'Inner geometry')
        require(self.original_config['material']['epsilon']==cfg['epsilon'],'Inner units')
        require(self.original_config['constructions']==cfg['constructions'],'H5/H6 stiffness/scaling')
        for path in (ARCHIVES/cfg['archives']['leading_cache']/'points').glob('*.json'):
            r=json.loads(path.read_text())
            if r['tau']!='RIGID':self.index[r['order'],float(r['tau'])]=(path,r)
        self.events.append(dict(phase='leading_sources_only',full_thresholds_parsed=False))

    def leading(self,p,tau,independent=False):
        key=p,float(tau);cfg=self.cfg;pol=cfg['policy'];m=self.models[p]
        if key not in self.solutions:
            if key in self.index:
                path,old=self.index[key]
                with np.load(path.with_suffix('.npz')) as data:coefficients=data['coefficients'].copy()
                s=dict(zeta=np.array(old['zeta']),coefficients=coefficients,diagnostics=old['diagnostics'],per_mode=old['per_mode'])
                origin=str(path.relative_to(ROOT));self.reused+=1
            else:s=m.solve(tau);origin='new_HS_LA6';self.new+=1
            for field,tol in [('strong_ode','strong_tol'),('natural_boundary','boundary_tol'),('essential','essential_tol'),('matrix_residual','matrix_tol'),('energy','energy_tol'),('mass_gram','mass_tol'),('symmetry','matrix_tol')]:
                require(s['diagnostics'][field]<pol[tol],f'p{p},tau={tau}: {field}')
            require(s['zeta'][0]>0 and np.all(np.diff(s['zeta'])>0),'Positive sorted prefix')
            if p==24:self.counts.append(dict(order=p,tau=tau,checks=check_prefix(cfg['k_joint'],tau,s['zeta'],pol)))
            self.records.append(dict(order=p,tau=tau,origin=origin,zeta=s['zeta'].tolist(),diagnostics=s['diagnostics'],per_mode=s['per_mode']))
            self.solutions[key]=s
        s=self.solutions[key]
        if independent:
            checks=check_prefix(cfg['k_joint'],tau,s['zeta'],pol,refine=True)
            self.counts.append(dict(order=p,tau=tau,role='final independent ODE',checks=checks))
        return s

    def predictions(self):
        cfg=self.cfg;pol=cfg['policy'];rows=[];convergence=[]
        for p in cfg['orders']:
            c=self.coefficients[p];zinf=c['zeta_infty'];a=c['a_reaction'];b=c['modal']['b']
            for q in cfg['levels']:
                target=q*q*self.zout;delta=zinf-target
                require(delta>0,'Target at/above inner saturation: no finite threshold')
                low,high=cfg['bracket']
                fun=lambda t:self.leading(p,float(t))['zeta'][0]/target-1
                f_low=fun(low);f_high=fun(high)
                require(f_low<0<f_high,'Scalar root not bracketed')
                tau=brentq(fun,low,high,xtol=pol['tau_root_xtol'],rtol=pol['tau_root_rtol'],maxiter=pol['max_root_iterations'])
                s=self.leading(p,tau,independent=(p==24));residual=float(fun(tau));h=pol['tau_bracket_halfwidth']
                left,right=tau-h,tau+h;fl,fr=float(fun(left)),float(fun(right))
                require(fl<-pol['bracket_sign_min'] and fr>pol['bracket_sign_min'],'Threshold signs unresolved')
                require(abs(residual)<pol['threshold_equation_tol'],'Direct threshold equation')
                disc=a*a+4*b*delta;require(disc>0,'Two-term discriminant nonpositive')
                tau_two=(a+math.sqrt(disc))/(2*delta);tau_other=(a-math.sqrt(disc))/(2*delta)
                D=float(self.models[p].D@s['coefficients'][:,0]);PD=tau*D*D/s['zeta'][0]
                # Independent forward evaluation checks, never fits, the two-term forecast.
                two_s=self.leading(p,tau_two)
                rows.append(dict(order=p,q=q,target_zeta=target,zeta_infty=zinf,zeta_outer=self.zout,
                    omega_outer=cfg['epsilon']*math.sqrt(self.zout),delta=delta,a=a,b=b,
                    tau=tau,tau_lower=left,tau_upper=right,sign_lower=fl,sign_upper=fr,equation_residual=residual,
                    discriminant=disc,tau_two=tau_two,tau_other=tau_other,tau_one=a/delta,
                    tau_two_relative_error=tau_two/tau-1,angle_two_relative_error=math.sqrt(tau_two/tau)-1,
                    forward_frequency_error_two=math.sqrt(two_s['zeta'][0]/target)-1,
                    PD=PD,D_squared=D*D))
        for q in cfg['levels']:
            low,high=[next(r for r in rows if r['q']==q and r['order']==p) for p in cfg['orders']]
            # Same tau, not two fitted equal-level eigenvalues, for p-convergence.
            t=high['tau'];s0=self.leading(20,t);s1=self.leading(24,t)
            freq=float(max(abs(np.sqrt(s0['zeta']/s1['zeta'])-1)));va=[];vb=[];weights=[]
            for arm,L in enumerate(self.models[24].lengths):
                va.append(self.models[20].field_matrices(arm,self.models[24].nodes)[0]@s0['coefficients'])
                vb.append(self.models[24].values[arm][0]@s1['coefficients']);weights.extend(self.models[24].weights*L)
            va=np.vstack(va);vb=np.vstack(vb);w=np.array(weights);sign=np.sign(np.sum(w[:,None]*va*vb,axis=0))
            shape=float(max(np.sqrt(np.sum(w[:,None]*(va-vb*sign)**2,axis=0))))
            tau_error=abs(low['tau']/high['tau']-1)
            require(freq<pol['frequency_convergence_tol'] and tau_error<pol['frequency_convergence_tol'] and shape<pol['shape_convergence_tol'],'Order convergence')
            convergence.append(dict(q=q,frequency_relative_error=freq,tau_relative_error=tau_error,shape_mass_error=shape))
        angles=[]
        for c in cfg['constructions']:
            for r in rows:
                if r['order']!=24:continue
                beta=lambda t:math.degrees(cfg['epsilon']*math.sqrt(t/c['g']))
                angles.append(dict(id=c['id'],q=r['q'],beta_exact=beta(r['tau']),beta_two=beta(r['tau_two']),
                    beta_one=beta(r['tau_one']),beta_lower=beta(r['tau_lower']),beta_upper=beta(r['tau_upper']),
                    xi_exact=math.sqrt(r['tau']/c['g']),xi_two=math.sqrt(r['tau_two']/c['g'])))
        write(self.folder/'apriori_predictions.json',dict(status='PASS',inputs='HS-LA3 coefficients and HS-LA5 outer eigenvalue only',
            full_data_used=False,thresholds=rows,angles=angles,convergence=convergence))
        self.events.append(dict(phase='predictions_sealed',sha256=sha(self.folder/'apriori_predictions.json'),full_thresholds_parsed=False))
        write(self.folder/'access_order.json',self.events)
        return rows,angles

    def compare_archives(self,rows,angles):
        seal=self.events[-1]
        require(seal['phase']=='predictions_sealed' and seal['sha256']==sha(self.folder/'apriori_predictions.json'),'Forecasts must be frozen first')
        cfg=self.cfg;src=ARCHIVES/cfg['archives']['full_thresholds']
        refs=read_csv(src/'thresholds.csv')
        pilot=read_csv(ARCHIVES/cfg['archives']['pilot']/'spectra.csv')
        full_points=read_csv(src/'full_points.csv')
        self.events.append(dict(phase='full_reference_values_parsed',thresholds_sha256=sha(src/'thresholds.csv')))
        comparison=[];decomposition=[];frequency=[];audits=[]
        for forecast in angles:
            name,q=forecast['id'],forecast['q']
            ref=next(r for r in refs if r['id']==name and float(r['q'])==q)
            beta_full=float(ref['beta_full_deg']);diagnostic=float(ref['beta_diagnostic_full_level_deg'])
            require(float(ref['bracket_F_lower'])<0<float(ref['bracket_F_upper']),'Archived full bracket signs')
            require(abs(float(ref['equation_residual']))<=self.original_config['full_threshold']['frequency_equation_tol'],'Archived full equation')
            point=next(r for r in full_points if r['id']==name and float(r['beta_deg'])==beta_full)
            physical=json.loads((ROOT/point['path']).read_text())
            require(point['status']==physical['status']=='PASS' and physical['count']['status']=='PASS','Full prefix not accepted')
            maxrow=next(r for r in pilot if r['id']==name and float(r['beta_deg'])==90)
            require(maxrow['status']=='PASS','Full reference maximum status')
            wfull=float(maxrow['omega1']);wout=cfg['epsilon']*math.sqrt(self.zout)
            main=next(r for r in rows if r['order']==24 and r['q']==q)
            # Exact additive decomposition in angles, common full-angle denominator.
            curve=diagnostic-beta_full;level=forecast['beta_exact']-diagnostic
            asym=forecast['beta_two']-forecast['beta_exact']
            decomp=dict(id=name,q=q,beta_full_level_diagnostic=diagnostic,
                curve_deg=curve,level_deg=level,asymptotic_deg=asym,
                curve_percent=100*curve/beta_full,level_percent=100*level/beta_full,
                asymptotic_percent=100*asym/beta_full,
                total_exact_percent=100*(curve+level)/beta_full,total_two_percent=100*(curve+level+asym)/beta_full)
            result=dict(**forecast,beta_full=beta_full,beta_old=float(ref['beta_lead_deg']),
                full_bracket_lower=float(ref['beta_full_bracket_lower']),full_bracket_upper=float(ref['beta_full_bracket_upper']))
            for method in ['exact','two','old','one']:
                delta=result['beta_'+method]-beta_full
                result[method+'_signed_error_deg']=delta;result[method+'_absolute_error_deg']=abs(delta)
                result[method+'_relative_error_percent']=100*delta/beta_full
            require(abs(decomp['total_two_percent']-result['two_relative_error_percent'])<1e-12,'Error decomposition')
            comparison.append(result);decomposition.append(decomp)
            frequency.append(dict(id=name,q=q,omega_full_90=wfull,omega_outer=wout,
                level_frequency_error_percent=100*(wout/wfull-1),PD=main['PD'],
                level_angle_linear_estimate_percent=100*math.log(wout/wfull)/main['PD'],
                exact_level_angle_percent=decomp['level_percent']))
            audits.append(dict(id=name,q=q,full_point_path=point['path'],full_point_sha256=sha(ROOT/point['path']),
                full_equation_residual=float(ref['equation_residual']),bracket_width=result['full_bracket_upper']-result['full_bracket_lower'],
                diagnostic_angle_role='archived leading curve normalized by full level, used only after sealing'))
        csvwrite(self.folder/'comparison.csv',comparison);csvwrite(self.folder/'error_components.csv',decomposition)
        csvwrite(self.folder/'level_sensitivity.csv',frequency);write(self.folder/'full_reference_audit.json',audits)
        write(self.folder/'access_order.json',self.events)

    def run(self):
        self.audit()
        if (self.folder/'manifest.json').exists():
            manifest=json.loads((self.folder/'manifest.json').read_text());require(manifest['identity']==self.identity,'Identity')
            for path,digest in manifest['artifacts'].items():require(sha(self.folder/path)==digest,'Artifact '+path)
            print('HS-LA6 verified',self.folder);return self.folder
        self.folder.mkdir(parents=True,exist_ok=True);write(self.folder/'contract.json',self.identity)
        try:
            rows,angles=self.predictions()
            print('HS-LA6 forecasts sealed before full-threshold parsing',flush=True)
            self.compare_archives(rows,angles)
            write(self.folder/'leading_points.json',self.records);write(self.folder/'independent_checks.json',self.counts)
            arrays={f'p{p}_{float(t).hex()}':s['coefficients'] for (p,t),s in self.solutions.items()}
            np.savez_compressed(self.folder/'leading_coefficients.npz',**arrays)
            write(self.folder/'manifest.json',dict(status='PASS',identity=self.identity,new_leading_prefixes=self.new,
                reused_leading_prefixes=self.reused,new_full_model_solves=0,new_outer_solves=0,
                full_thresholds_reused=4,closed_form_predictions=4,stop=True,
                artifacts={p.relative_to(self.folder).as_posix():sha(p) for p in self.folder.rglob('*') if p.is_file()}))
            print('HS-LA6 PASS',self.folder);return self.folder
        except Exception as e:
            write(self.folder/'failure.json',dict(status='FAIL',error=repr(e),new_leading_prefixes=self.new,reused_leading_prefixes=self.reused,events=self.events))
            raise


def run():return Study(json.loads(CONFIG.read_text())).run()
