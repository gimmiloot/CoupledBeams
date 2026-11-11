"""HS-D1 bounded eigenvalue composition; no new physical assembly or full solves.

Reuses frozen HS-L1, HS-C, saturation coefficients and endpoint counters.
The two eigenvalues are combined, never their matrices or eigenvectors.
"""
import argparse
import hashlib
import json
import math
import platform
import time
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import brentq
from scripts.lib import isotropic_two_springs_horizontal_leading as inner
from scripts.lib import two_springs_horizontal_finite_angle as outer
from scripts.lib.two_springs_horizontal_joint_stiffness import endpoint_count
from scripts.lib.two_springs_horizontal_saturation import coefficient_from_linear
from scripts.lib.two_springs_horizontal_spectral_structure import reaction
from scripts.lib.two_springs_horizontal_angle_map import ROOT, read, atomic, sha, require, csvwrite, crossing_brackets

CONFIG=ROOT/'data/input/isotropic_two_springs_horizontal_composite_approximation.json'
BASE=ROOT/'results/isotropic_two_springs_rl_timoshenko/horizontal_composite_approximation'
RIGID='RIGID'


def g(chi):return .5 if chi==RIGID else chi/(1+2*chi)


def composition(zc,zl,zinf):return zc+zl-zinf


def identity(config=CONFIG):
    cfg=read(config)
    return dict(config=cfg,config_sha256=sha(config),initial_sha256=sha(ROOT/cfg['initial_workspace']),
        analytical_gate_sha256=sha(ROOT/cfg['analytical_gate']),versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))


def output_folder(config=CONFIG,base=BASE):return base/hashlib.sha256(json.dumps(identity(config),sort_keys=True).encode()).hexdigest()[:16]


def inner_count(z,k,tau,mu,policy,refine=False):
    rows=[]
    for i,value in enumerate(z):
        lo=value*(1-policy['count_bracket_relative']);hi=value*(1+policy['count_bracket_relative'])
        a=endpoint_count(lo,k,tau,mu,policy['inertia_relative_zero_tol'],policy['pole_relative_guard'])
        b=endpoint_count(hi,k,tau,mu,policy['inertia_relative_zero_tol'],policy['pole_relative_guard'])
        require(a['count']==i and b['count']==i+1,'HS-L1 prefix count')
        rows.append(dict(position=i+1,lower=a,upper=b))
    checks=inner.independent_roots(z,tau,k,mu) if refine else []
    require(all(abs(c['relative_frequency_error'])<=policy['independent_frequency_tol'] for c in checks),'HS-L1 independent roots')
    return dict(counts=rows,roots=checks)


class Study:
    def __init__(self,config=CONFIG,base=BASE):
        self.cfg=read(config);self.p=self.cfg['policy'];self.models={};self.solutions={}
        require(read(ROOT/self.cfg['analytical_gate'])['status']=='PASS','Analytical gate')
        self.identity=identity(config)
        self.id=hashlib.sha256(json.dumps(self.identity,sort_keys=True).encode()).hexdigest()[:16]
        self.folder=base/self.id;self.folder.mkdir(parents=True,exist_ok=True)
        self.points=self.folder/'points';self.points.mkdir(exist_ok=True)
        version=self.folder/'workflow_versions'/(sha(__file__)+'.txt');version.parent.mkdir(exist_ok=True)
        if not version.exists():version.write_bytes(Path(__file__).read_bytes())
        if (self.folder/'contract.json').exists():require(read(self.folder/'contract.json')==self.identity,'Identity changed')
        else:atomic(self.folder/'contract.json',self.identity)
        self.records={p.stem:read(p) for p in self.points.glob('*.json')}
        self.old_c=read(ROOT/self.cfg['outer_source']/'points.json');self.old_l={}
        for path in (ROOT/self.cfg['inner_source']/'leading').glob('*.json'):
            row=read(path);self.old_l[row['tau']]=(row,path)

    def audit(self):
        for path,digest in self.cfg['sources'].items():require(sha(ROOT/path)==digest,'Frozen source '+path)
        for name,digest in self.cfg['archives'].items():
            folder=ROOT/name;require(sha(folder/'manifest.json')==digest,'Archive manifest')
            m=read(folder/'manifest.json');require(m['status']=='PASS','Source status')
            for p,h in m['artifacts'].items():require(sha(folder/p)==h,'Source artifact '+p)
        require(sha(ROOT/self.cfg['full_source']/'spectra.csv')==self.cfg['full_source_csv_sha256'],'Full CSV')
        atomic(self.folder/'source_audit.json',dict(status='PASS',new_full_solves=0))

    def model(self,kind,k,mu,p,param,ks):
        key=(kind,k,mu,p,param if kind=='C' else None,ks)
        if key not in self.models:
            self.models[key]=inner.Leading(k,mu,p) if kind=='L1' else outer.FiniteAngle(k,ks,param,mu,p)
        return self.models[key]

    def validate(self,row,kind):
        d=row['diagnostics'];names=[('natural_boundary' if kind=='L1' else 'natural','boundary_tol'),
            ('strong_ode' if kind=='L1' else 'strong','strong_tol'),('essential','essential_tol'),('matrix_residual','matrix_tol'),
            ('symmetry','matrix_tol'),('mass_gram','mass_tol'),('energy','energy_tol')]
        for name,tol in names:require(d[name]<=self.p[tol],f'{kind}: {name}={d[name]}')
        require(d['mass_positive'] and d['stiffness_positive'],'Positive forms')
        require(all(v>0 for v in row['zeta']) and all(a<b for a,b in zip(row['zeta'][:-1],row['zeta'][1:])),'Positive sorted prefix')

    def point(self,kind,k,mu,param,ks=None,role='grid',independent=False):
        param=param if param==RIGID else float(param)
        key=hashlib.sha256(repr((kind,float(k),float(mu),param,ks)).encode()).hexdigest()[:16]
        if key in self.records:return self.records[key]
        rows=[];arrays={};origin='new';source=None
        try:
            if kind=='L1' and k==.1 and mu==.01 and param in self.old_l:
                old,path=self.old_l[param];rows=old['orders'];origin='reused_HS_L1';source=path.relative_to(ROOT).as_posix()
                with np.load(path.with_suffix('.npz')) as f:arrays={n:f[n].copy() for n in f.files}
                freq=old['frequency_convergence'];shape=old['shape_convergence']
            elif kind=='C':
                match=next(((name,r) for name,r in self.old_c.items() if r['k_joint']==k and r['mu']==mu and r['beta']==param and r['kappa_s']==ks),None)
                if match:
                    name,old=match;rows=old['orders'];origin='reused_HS_C';source=self.cfg['outer_source']+'/points.json'
                    with np.load(ROOT/self.cfg['outer_source']/(name+'.npz')) as f:arrays={'p24_coefficients':f['coefficients'].copy()}
                    freq=old['frequency_convergence'];shape=old['shape_convergence']
            if not rows:
                for p in self.cfg['orders']:
                    m=self.model(kind,k,mu,p,param,ks);s=m.solve(param) if kind=='L1' else m.solve()
                    rows.append(dict(order=p,zeta=s['zeta'].tolist(),diagnostics=s['diagnostics'],per_mode=s['per_mode']))
                    arrays[f'p{p}_coefficients']=s['coefficients']
                freq=float(max(abs(np.sqrt(np.array(rows[0]['zeta'])/rows[1]['zeta'])-1)))
                # Same-point p-convergence in the physical mass norm, no tracking.
                fields=[];m1=self.model(kind,k,mu,24,param,ks);b=m1 if kind=='L1' else m1.base
                weights=np.r_[b.weights*(1-mu),b.weights*(1+mu)]
                for p in self.cfg['orders']:
                    m=self.model(kind,k,mu,p,param,ks);bb=m if kind=='L1' else m.base;c=arrays[f'p{p}_coefficients']
                    values=[bb.field_matrices(i,b.nodes)[0]@c for i in [0,1]]
                    if kind=='C':
                        A=m.A@c;values[1]-=m.c*A;values.append((math.sqrt(1+mu)*m.s*A)[None,:])
                    fields.append(np.vstack(values))
                if kind=='C':weights=np.r_[weights,1.]
                signs=np.sign(np.sum(weights[:,None]*fields[0]*fields[1],axis=0))
                shape=float(max(np.sqrt(np.sum(weights[:,None]*(fields[0]-fields[1]*signs)**2,axis=0))))
            for row in rows:self.validate(row,kind)
            require(freq<=self.p['frequency_convergence_tol'],'p20/24 frequency')
            require(shape<=self.p['shape_convergence_tol'],'p20/24 shape')
            z=np.array(rows[-1]['zeta'])
            if kind=='L1':checks=inner_count(z,k,param,mu,self.p,independent)
            elif origin=='new':checks=outer.independent(z,k,ks,param,mu,self.p)
            else:checks=rows[-1]['independent']
            r=dict(kind=kind,k_joint=k,mu=mu,parameter=param,kappa_s=ks,role=role,origin=origin,source=source,
                source_sha256=sha(ROOT/source) if source else None,status='PASS',orders=rows,zeta=z.tolist(),
                frequency_convergence=freq,shape_convergence=shape,independent=checks,point_key=key)
            temp=self.points/(key+'.tmp.npz');np.savez_compressed(temp,**arrays);temp.replace(self.points/(key+'.npz'))
            atomic(self.points/(key+'.json'),r);self.records[key]=r;return r
        except Exception as exc:
            atomic(self.folder/'failure.json',dict(status='UNRESOLVED',kind=kind,k=k,mu=mu,parameter=param,ks=ks,role=role,partial_orders=rows,error=repr(exc)))
            raise

    def coefficients(self):
        path=self.folder/'coefficients.json'
        if path.exists():return read(path)
        out=[]
        for case in self.cfg['cases']:
            k=case['k_joint'];mu=case['mu'];eps=case['epsilon']
            inf=self.point('L1',k,mu,RIGID,role='exact common limit',independent=True)
            c0=self.point('C',k,mu,0.,RIGID,role='right limit 0+')
            c90=self.point('C',k,mu,90.,RIGID,role='positivity endpoint')
            zero=self.point('L1',k,mu,0.,role='positive lower bound')
            orders=[]
            for p in self.cfg['orders']:
                j=self.cfg['orders'].index(p);m=self.model('L1',k,mu,p,RIGID,None)
                with np.load(self.points/(inf['point_key']+'.npz')) as f:c=f[f'p{p}_coefficients']
                sol=dict(zeta=np.array(inf['orders'][j]['zeta']),coefficients=c)
                _,r=reaction(m,sol);Z=m.Z
                coeff,_,_=coefficient_from_linear(Z.T@m.K0@Z,Z.T@m.M@Z,Z.T@m.D,float(sol['zeta'][0]))
                require(max(coeff['residual_u'],coeff['residual_v'])<=self.p['matrix_tol'],'Saturation linear equations')
                require(abs(coeff['G'])<=self.p['matrix_tol'],'Saturation G=0')
                require(abs(coeff['a']/r['r_squared']-1)<=self.p['reaction_relative_tol'],'Reaction coefficient')
                require(r['kkt_residual']<=self.p['matrix_tol'],'Saturation KKT')
                zinf=float(sol['zeta'][0]);zc=c0['orders'][j]['zeta'][0];pm=c0['orders'][j]['per_mode'][0]
                C2=zc*pm['mixed_mass']/2
                require(abs(math.sqrt(zc/zinf)-1)<=self.p['independent_frequency_tol'],'Common limiting eigenvalue')
                lower=min(zc,c90['orders'][j]['zeta'][0])+zero['orders'][j]['zeta'][0]-zinf
                require(lower>0,'Global Rayleigh positivity bound')
                orders.append(dict(order=p,zeta_infty=zinf,zeta_C_zero_plus=zc,zeta_C90=c90['orders'][j]['zeta'][0],
                    zeta_L1_zero=zero['orders'][j]['zeta'][0],a=r['r_squared'],b=coeff['b'],C2=C2,A=pm['A'],
                    integral_V=-pm['mixed_mass']/(2*pm['A']),reaction=r,linear=coeff,global_composite_lower_bound=lower))
            for name in ['a','b','C2']:
                require(abs(orders[0][name]/orders[1][name]-1)<=self.p['coefficient_relative_tol'],'p20/24 '+name)
            tests=[];zinf=orders[-1]['zeta_C_zero_plus'];C2=orders[-1]['C2']
            for h in self.cfg['diagnostics']['C2_angles_deg']:
                v=self.point('C',k,mu,h,RIGID,role='C2 check')['zeta'][0]
                estimate=(v-zinf)/math.radians(h)**2
                tests.append(dict(beta_deg=h,C2_quotient=estimate,error_relative=estimate/C2-1))
            rich=(4*tests[-1]['C2_quotient']-tests[-2]['C2_quotient'])/3
            require(abs(rich/C2-1)<=self.cfg['diagnostics']['C2_relative_tol'],'C2 independent difference')
            transition=(orders[-1]['a']/(g(case['chi_s'])*abs(C2)))**.25*math.sqrt(eps) if C2 else None
            out.append(dict(id=case['id'],orders=orders,C2_checks=tests,C2_Richardson=rich,
                C2_Richardson_relative_error=rich/C2-1,beta_transition_rad=transition,
                beta_transition_deg=None if transition is None else math.degrees(transition),
                transition_role='not applied to control I' if case['role']=='control_I' else 'formal balance; not validity boundary'))
        atomic(path,out);return out

    def value(self,case,beta,role='grid'):
        k=case['k_joint'];mu=case['mu'];eps=case['epsilon'];beta=float(beta);rad=math.radians(beta)
        tau=g(case['chi_s'])*(rad/eps)**2;control=case['role']=='control_I';ks=case['kappa_s'] if control else RIGID
        l=self.point('L1',k,mu,tau,role=role,independent=role!='grid')
        c=self.point('C',k,mu,beta,ks,role=role)
        zinf=self.point('L1',k,mu,RIGID)['zeta'][0]
        zl=l['zeta'][0];zc=c['zeta'][0];z=zc if control else composition(zc,zl,zinf)
        require(z>0,'Nonpositive composite value')
        lowc=c['orders'][0]['zeta'][0];lowl=l['orders'][0]['zeta'][0]
        lowinf=self.point('L1',k,mu,RIGID)['orders'][0]['zeta'][0]
        low=lowc if control else composition(lowc,lowl,lowinf)
        require(abs(math.sqrt(low/z)-1)<=self.p['frequency_convergence_tol'],'Composite p20/24')
        pm=c['orders'][-1]['per_mode'][0]
        # B(X,X)=-2*A*integral(V), evaluated without dividing by cos(beta).
        if abs(math.cos(rad))>1e-12:B=pm['mixed_mass']/math.cos(rad)
        else:
            with np.load(self.points/(c['point_key']+'.npz')) as f:coef=f['p24_coefficients'][:,0]
            model=self.model('C',k,mu,24,beta,ks);B=float(coef@model.cross@coef)
        dzc=zc*math.sin(rad)*B
        if control:dz=dzc+2*ks*math.sin(rad)*math.cos(rad)*pm['V0']**2
        else:dz=dzc+2*g(case['chi_s'])*rad/eps**2*l['orders'][-1]['per_mode'][0]['D']**2
        return dict(id=case['id'],beta_deg=beta,tau=tau,xi=rad/eps,zeta_L1=zl,zeta_C=zc,zeta_infty=zinf,
            zeta_comp=None if control else z,zeta_primary=z,omega_L1=eps*math.sqrt(zl),omega_C=eps*math.sqrt(zc),
            omega_comp=None if control else eps*math.sqrt(z),omega_primary=eps*math.sqrt(z),
            slope_log_per_rad=dz/(2*z),slope_log_per_degree=dz/(2*z)*math.pi/180,
            p20_p24_frequency_error=abs(math.sqrt(low/z)-1),inner_key=l['point_key'],outer_key=c['point_key'],
            primary='C_I' if control else 'composite',zero_C_semantics='ordinary' if control else 'right_limit_0plus')

    def preflight(self):
        path=self.folder/'preflight.json'
        if path.exists():return
        start=time.perf_counter();self.audit();self.coefficients();rows=[]
        for case in self.cfg['cases']:
            for beta in [.1,90.]:rows.append(self.value(case,beta,'preflight high tau and small angle'))
        atomic(path,dict(status='PASS',elapsed_seconds=time.perf_counter()-start,values=rows,finite_tau_not_replaced_by_infinity=True))
        print('PREFLIGHT PASS',self.folder,flush=True)

    def predictions(self,curves):
        path=self.folder/'predictions.json'
        if path.exists():return read(path)
        results=[]
        for case in self.cfg['cases']:
            data=sorted([r for r in curves if r['id']==case['id']],key=lambda r:r['beta_deg'])
            angles=[r['beta_deg'] for r in data];freq=[r['omega_primary'] for r in data];level=freq[-1]
            for q in self.cfg['thresholds']['q']:
                target=q*level;crossings=[]
                def fun(b):return self.value(case,b,'model threshold')['omega_primary']/target-1
                for lo,hi,direction in crossing_brackets(angles,freq,target):
                    beta=brentq(fun,lo,hi,xtol=self.cfg['thresholds']['root_xtol_deg'],rtol=self.cfg['thresholds']['root_rtol'])
                    half=self.cfg['thresholds']['bracket_width_deg']/2;left=beta-half;right=beta+half
                    fl,fr=fun(left),fun(right);v=self.value(case,beta,'threshold slope');res=fun(beta)
                    require(fl*fr<0 and abs(res)<=self.cfg['thresholds']['residual_tol'],'Threshold bracket/residual')
                    differences=[]
                    for h in self.cfg['diagnostics']['slope_steps_deg']:
                        plus=self.value(case,beta+h,'slope check')['omega_primary'];minus=self.value(case,beta-h,'slope check')['omega_primary']
                        derivative=math.log(plus/minus)/(2*h);err=derivative/v['slope_log_per_degree']-1
                        differences.append(dict(h_deg=h,slope_log_per_degree=derivative,relative_error=err))
                    require(abs(differences[-1]['relative_error'])<=self.cfg['diagnostics']['slope_relative_tol'],'Threshold slope check')
                    crossings.append(dict(beta_deg=beta,direction=direction,bracket=[left,right],signs=[fl,fr],residual=res,
                        xi=v['xi'],slope_log_per_degree=v['slope_log_per_degree'],slope_log_per_rad=v['slope_log_per_rad'],
                        frequency_derivative_per_degree=target*v['slope_log_per_degree'],differences=differences))
                status='PASS' if crossings else 'ALL_SAMPLED_HIGH' if min(freq)>=target else 'NO_CROSSING_FOUND'
                results.append(dict(id=case['id'],q=q,omega_reference=level,primary=data[0]['primary'],status=status,crossings=crossings,
                    decreasing_steps=sum(b<a for a,b in zip(freq[:-1],freq[1:])),count_up=sum(r['direction']=='up' for r in crossings),
                    count_down=sum(r['direction']=='down' for r in crossings)))
                print('prediction',case['id'],q,status,[r['beta_deg'] for r in crossings],flush=True)
        atomic(path,results);return results

    def run(self):
        if (self.folder/'manifest.json').exists():self.verify();print('REUSED no new solutions');return
        self.preflight();start=time.perf_counter();curves=[]
        for case in self.cfg['cases']:
            for beta in self.cfg['angles_deg'][case['id']]:curves.append(self.value(case,beta))
            print('model curve',case['id'],'PASS',flush=True)
        csvwrite(self.folder/'model_curves.csv',curves)
        predictions=self.predictions(curves);prediction_hash=sha(self.folder/'predictions.json')
        source=ROOT/self.cfg['full_source'];full=outer.rows_csv(source/'spectra.csv')
        refs={(r['id'],float(r['beta_deg'])):r for r in full};comparison=[]
        for v in curves:
            f=refs[v['id'],v['beta_deg']];case=next(c for c in self.cfg['cases'] if c['id']==v['id'])
            require(all(float(f[a])==case[b] for a,b in [('epsilon','epsilon'),('mu','mu'),('aspect','aspect'),('kappa_j','k_joint')]),'Full physical identity')
            require(f['chi_s']==RIGID if case['chi_s']==RIGID else float(f['chi_s'])==case['chi_s'],'Full support identity')
            record=read(ROOT/f['point_file']);require(record['status']=='PASS' and record['count']['status']=='PASS','Full source point')
            omega=float(f['omega1']);require(omega==record['omega'][0],'Full CSV frequency')
            row=dict(**v,omega_full=omega,zeta_full=(omega/case['epsilon'])**2,full_point=f['point_file'],full_sha256=sha(ROOT/f['point_file']))
            for kind in ['L1','C','comp','primary']:
                w=v['omega_'+kind];row['error_'+kind+'_percent']=None if w is None else 100*(w/omega-1)
            if v['omega_comp'] is not None:
                errors=[abs(row['error_'+k+'_percent']) for k in ['L1','C','comp']]
                row['worse_than_both']=errors[2]>max(errors[:2]);row['better_than_both']=errors[2]<min(errors[:2])
            comparison.append(row)
        csvwrite(self.folder/'comparison.csv',comparison)
        summaries=[]
        for case in self.cfg['cases']:
            for label,lo,hi in [('0_5',0,5),('5_30',5,30),('30_90',30,90),('all',0,90)]:
                rr=[r for r in comparison if r['id']==case['id'] and lo<=r['beta_deg']<=hi]
                s=dict(id=case['id'],interval=label,n=len(rr))
                for kind in ['L1','C','comp','primary']:
                    valid=[r for r in rr if r['error_'+kind+'_percent'] is not None]
                    s['max_abs_'+kind+'_percent']=max(abs(r['error_'+kind+'_percent']) for r in valid) if valid else None
                    s['min_'+kind+'_percent']=min(r['error_'+kind+'_percent'] for r in valid) if valid else None
                    s['max_'+kind+'_percent']=max(r['error_'+kind+'_percent'] for r in valid) if valid else None
                s['worse_than_both_points']=sum(r.get('worse_than_both',False) for r in rr);summaries.append(s)
        csvwrite(self.folder/'interval_errors.csv',summaries)
        thresholds=read(source/'thresholds.json');threshold_rows=[]
        for prediction in predictions:
            case=next(c for c in self.cfg['cases'] if c['id']==prediction['id'])
            fullrow=next(r for r in thresholds if r['id']==case['id'] and r['normalization']=='omega90' and r['q']==prediction['q'])
            beta_full=fullrow['terminal_beta_root_deg'];value_full=self.value(case,beta_full,'archived full threshold comparison')
            flevel=fullrow['omega_target'];level_log=math.log(prediction['omega_reference']/flevel)
            curve_log=math.log(value_full['omega_primary']/(prediction['q']*flevel))
            slope=value_full['slope_log_per_degree']
            for crossing in prediction['crossings']:
                root=crossing['beta_deg'];threshold_rows.append(dict(id=case['id'],q=prediction['q'],status=prediction['status'],direction=crossing['direction'],
                    beta_model=root,beta_full=beta_full,error_deg=root-beta_full,error_percent=100*(root/beta_full-1),xi=crossing['xi'],
                    slope_log_per_degree=crossing['slope_log_per_degree'],full_angle_model_slope=slope,reference_frequency_error_percent=100*(math.exp(level_log)-1),
                    curve_frequency_error_at_full_percent=100*(math.exp(curve_log)-1),linearized_level_shift_deg=level_log/slope,
                    linearized_curve_shift_deg=-curve_log/slope,linearized_total_shift_deg=(level_log-curve_log)/slope,
                    slope_FD_relative_error=crossing['differences'][-1]['relative_error']))
        csvwrite(self.folder/'threshold_comparison.csv',threshold_rows)
        require(sha(self.folder/'predictions.json')==prediction_hash,'Predictions changed after full comparison')
        atomic(self.folder/'prediction_provenance.json',dict(predictions_sha256=prediction_hash,full_thresholds_sha256=sha(source/'thresholds.json'),
            forecast_before_threshold_table_read=True,prior_HSB_knowledge=True,fit_parameters=0))
        neighbours=[]
        for case in self.cfg['cases']:
            rr=[r for r in curves if r['id']==case['id']]
            for a,b in zip(rr[:-1],rr[1:]):
                jump=abs(b['omega_primary']-a['omega_primary'])/max(a['omega_primary'],b['omega_primary'])
                neighbours.append(dict(id=case['id'],left=a['beta_deg'],right=b['beta_deg'],jump=jump,
                    triggered=jump>self.cfg['frequency_map']['neighbour_jump_relative'],resolution='both component prefixes independently counted'))
        atomic(self.folder/'neighbour_audit.json',neighbours)
        manifest=dict(status='PASS',identity=self.identity,unique_leading_points=len(self.records),new_points=sum(r['origin']=='new' for r in self.records.values()),
            reused_HS_C=sum(r['origin']=='reused_HS_C' for r in self.records.values()),reused_HS_L1=sum(r['origin']=='reused_HS_L1' for r in self.records.values()),
            full_archived_grid_points=len(curves),full_thresholds=12,new_full_solves=0,production_seconds=time.perf_counter()-start,
            min_composite_zeta=min(r['zeta_comp'] for r in curves if r['zeta_comp'] is not None),
            unresolved_predictions=[r for r in predictions if r['status']!='PASS'])
        if manifest['unresolved_predictions']:manifest['status']='PARTIAL'
        analyze_saved(self.folder)
        manifest['artifacts']={p.relative_to(self.folder).as_posix():sha(p) for p in self.folder.rglob('*') if p.is_file() and p.name not in ['manifest.json','README.md','completion_checks.json','documentation_manifest.json','figures_manifest.json'] and 'figures' not in p.parts}
        atomic(self.folder/'manifest.json',manifest);print(manifest['status'],self.folder,flush=True)

    def verify(self):
        m=read(self.folder/'manifest.json')
        for path,digest in m['artifacts'].items():require(sha(self.folder/path)==digest,'Artifact '+path)
        print('VERIFIED',m['status'],self.folder,flush=True)


def analyze_saved(folder):
    folder=Path(folder);records=[read(p) for p in (folder/'points').glob('*.json')];coefficients=read(folder/'coefficients.json')
    maxima={};counts=0;root_errors=[];margins=[];tau=[];gaps=[]
    for r in records:
        if r['kind']=='L1' and r['parameter']!=RIGID:tau.append(r['parameter'])
        for row in r['orders']:
            d=row['diagnostics']
            for key in ['matrix_residual','mass_gram','energy','essential','symmetry']:
                maxima[key]=max(maxima.get(key,0),d[key])
            for target,source in [('natural','natural_boundary' if r['kind']=='L1' else 'natural'),('strong','strong_ode' if r['kind']=='L1' else 'strong')]:
                maxima[target]=max(maxima.get(target,0),d[source])
            z=row['zeta'];gaps.append((z[1]-z[0])/z[0])
        checks=r['independent']
        if r['kind']=='L1':
            counts+=len(checks['counts']);root_errors += [abs(x['relative_frequency_error']) for x in checks['roots']]
            for x in checks['counts']:margins += [x['lower']['inertia_margin'],x['upper']['inertia_margin']]
        else:
            counts+=len(checks);root_errors += [abs(x['frequency_error']) for x in checks]
            for x in checks:margins += [x['lower']['inertia_margin'],x['upper']['inertia_margin']]
    predictions=read(folder/'predictions.json');roots=[x for r in predictions for x in r['crossings']]
    atomic(folder/'numerical_checks.json',dict(status='PASS',max_residuals=maxima,max_p_frequency_error=max(r['frequency_convergence'] for r in records),
        max_p_shape_error=max(r['shape_convergence'] for r in records),counted_positions=counts,independent_refined_roots=len(root_errors),
        max_independent_frequency_error=max(root_errors),minimum_count_margin=min(margins),maximum_finite_tau=max(tau),
        minimum_relative_zeta_gap12=min(gaps),max_C2_Richardson_error=max(abs(c['C2_Richardson_relative_error']) for c in coefficients),
        max_threshold_residual=max(abs(x['residual']) for x in roots),max_slope_FD_small_step_error=max(abs(x['differences'][-1]['relative_error']) for x in roots),
        directions=[x['direction'] for x in roots],min_global_positive_bound=min(c['orders'][-1]['global_composite_lower_bound'] for c in coefficients),
        no_finite_tau_replaced_by_saturation=True,arithmetic='binary64, all frozen gates passed; no high-precision repair triggered'))
    # Only archived full points form these local secants; no full interpolation
    # is accepted as a spectrum or root. The secants are slope diagnostics.
    cfg=read(CONFIG);full=outer.rows_csv(ROOT/cfg['full_source']/'spectra.csv');thresholds=outer.rows_csv(folder/'threshold_comparison.csv');slopes=[]
    archived_thresholds=read(ROOT/cfg['full_source']/'thresholds.json')
    for row in thresholds:
        beta=float(row['beta_full']);rr=sorted([r for r in full if r['id']==row['id']],key=lambda r:float(r['beta_deg']))
        below=[r for r in rr if float(r['beta_deg'])<beta];above=[r for r in rr if float(r['beta_deg'])>beta]
        left,right=below[-1],above[0];a,b=float(left['beta_deg']),float(right['beta_deg'])
        nearest=math.log(float(right['omega1'])/float(left['omega1']))/(b-a)
        original=next(r for r in archived_thresholds if r['id']==row['id'] and r['q']==float(row['q']) and r['normalization']=='omega90')
        crossing=next(c for c in original['crossings'] if c['beta_deg']==beta);bl,br=crossing['bracket'];fl,fr=crossing['signs']
        secant=(math.log1p(fr)-math.log1p(fl))/(br-bl)
        slopes.append(dict(id=row['id'],q=row['q'],beta_full=beta,left=a,right=b,width_deg=b-a,
            nearest_secant_log_per_degree=nearest,archived_bracket_left=bl,archived_bracket_right=br,
            archived_bracket_secant_log_per_degree=secant,nearest_vs_bracket_relative_difference=nearest/secant-1,
            model_slope_at_full_degree=float(row['full_angle_model_slope']),
            role='original accepted full bracket is preferred; nearest-root iterates can be subtraction-sensitive; neither is an exact full derivative'))
    csvwrite(folder/'archived_slope_comparison.csv',slopes)
    print('saved-data analysis PASS',folder,flush=True)


def render(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder=Path(folder);rows=outer.rows_csv(folder/'comparison.csv');cfg=read(CONFIG);figdir=folder/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':True,'grid.alpha':.25,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    labels={'L1':'HS-L1','C':'HS-C II/III','comp':'Составное'};colors={'L1':'#0072B2','C':'#D55E00','comp':'#009E73'}
    titles={'B0_A1':'A1','B0_A3':'A3','B2_A3_m0.1':'A3, μ=0,1','B0_A4':'A4','B2_A4_m-0.1':'A4, μ=−0,1','B0_A2':'A2: контроль HS-C I'}
    for name,lo,hi,errors in [('01_frequency',0,90,False),('02_error_0_5',0,5,True),('03_error_5_30',5,30,True),('04_error_30_90',30,90,True)]:
        fig,axes=plt.subplots(2,3,figsize=(14,8))
        for ax,case in zip(axes.flat,cfg['cases']):
            rr=sorted([r for r in rows if r['id']==case['id'] and lo<=float(r['beta_deg'])<=hi],key=lambda r:float(r['beta_deg']))
            x=[float(r['beta_deg']) for r in rr]
            kinds=['L1','C'] if case['role']=='control_I' else ['L1','C','comp']
            for kind in kinds:
                y=[abs(float(r['error_'+kind+'_percent'])) if errors else float(r['omega_'+kind]) for r in rr]
                label='HS-C I' if kind=='C' and case['role']=='control_I' else labels[kind]
                ax.plot(x,y,'-' if kind=='comp' else '--' if kind=='L1' else ':',color=colors[kind],lw=1.6,label=label)
            if errors:ax.set_yscale('log')
            else:ax.plot(x,[float(r['omega_full']) for r in rr],'o',color='black',ms=2.5,markevery=max(1,len(rr)//25),label='Полная: архив')
            ax.set(title=titles[case['id']],xlim=(lo,hi),xlabel='Угол β, град',ylabel='Модуль ошибки частоты, %' if errors else 'Первая круговая частота ω₁')
            ax.legend(fontsize=8)
        fig.suptitle('Ошибки на архивных узлах; интервалы не являются границами применимости' if errors else 'HS-C II/III при β=0 означает предел β→0+; A2 не составляется',fontsize=11)
        fig.tight_layout();fig.savefig(figdir/(name+'.pdf'),bbox_inches='tight');fig.savefig(figdir/(name+'.png'),dpi=320,bbox_inches='tight');plt.close(fig)
    for name in ['comparison.csv','interval_errors.csv','threshold_comparison.csv']:(figdir/name).write_bytes((folder/name).read_bytes())
    atomic(folder/'figures_manifest.json',dict(calculation_mode='plot_only',new_spectral_solves=0,dpi=320,
        artifacts={p.name:sha(p) for p in figdir.iterdir() if p.is_file()}))
    print('plot_only',figdir,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['preflight','run','verify','analyze','plot_only','finalize','diagnose'])
    p.add_argument('--study',choices=['HS-D1','HS-D2'],default='HS-D1');args=p.parse_args()
    if args.study=='HS-D2':
        from scripts.lib.two_springs_horizontal_composite_geometry import dispatch
        dispatch(args.action);return
    if args.action=='diagnose':p.error('diagnose is a triggered HS-D2 action')
    if args.action=='plot_only':render(output_folder());return
    if args.action=='analyze':analyze_saved(output_folder());return
    s=Study()
    if args.action=='preflight':s.preflight()
    elif args.action=='run':s.run()
    elif args.action=='finalize':
        s.verify();m=read(s.folder/'manifest.json')
        for p in list((s.folder/'workflow_versions').glob('*.txt'))+[s.folder/'numerical_checks.json',s.folder/'archived_slope_comparison.csv']:
            m['artifacts'][p.relative_to(s.folder).as_posix()]=sha(p)
        atomic(s.folder/'manifest.json',m);print('sealed',s.folder)
    else:s.verify()


if __name__=='__main__':main()
