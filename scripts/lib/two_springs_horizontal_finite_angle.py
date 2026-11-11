"""HS-C fixed-angle thin-beam limits. See the canonical finite-angle note.

Energy Ritz uses F,V. Independent boundary/inertia checks use F,W=V-c*A.
The accepted RightAngle eigensolver and integrated fields are reused; no full
RL/Timoshenko spectral solver is called. beta=0,RIGID denotes 0+ ONLY.
"""
import argparse
import csv
import hashlib
import json
import math
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import expm, null_space, solve
from scipy.optimize import brentq
from scripts.lib.isotropic_two_springs_horizontal_leading import Leading
from scripts.lib.two_springs_horizontal_right_angle import RightAngle
from scripts.lib.two_springs_horizontal_angle_map import ROOT, read, atomic, sha, require, csvwrite, crossing_brackets

RIGID='RIGID'
CONFIG=ROOT/'data/input/isotropic_two_springs_horizontal_finite_angle_limit.json'
OUTPUT=ROOT/'results/isotropic_two_springs_rl_timoshenko/horizontal_finite_angle_limit'


def trig(beta):
    require(0<=beta<=90,'Finite-angle domain')
    return (0.,1.) if beta==90 else (1.,0.) if beta==0 else (math.cos(math.radians(beta)),math.sin(math.radians(beta)))


class FiniteAngle(RightAngle):
    def __init__(self,k_joint,kappa_s,beta,mu=.01,order=24,quadrature=80):
        require(k_joint>0,'HS-C positive joint stiffness contract')
        super().__init__(k_joint,mu,order,quadrature)
        self.kappa_s=kappa_s;self.beta=beta;self.c,self.s=trig(beta)
        self.M=self.M+self.c*self.cross
        if kappa_s!=RIGID:
            require(kappa_s>=0,'Nonnegative spring')
            self.C=self.base.ends[1][0][None,:]
            self.Z=null_space(self.C)
            d=self.base.outer[1][0]
            self.K=self.K+kappa_s*self.s**2*np.outer(d,d)

    def diagnose(self,z,coefficients):
        b=self.base;per=[]
        for i,ev in enumerate(z):
            x=coefficients[:,i];e=np.array([a@x for a in b.ends]);o=np.array([a@x for a in b.outer])
            A=float(self.A@x);J=float(b.jump@x);H=float(self.H@x)
            natural=[o[1,2],e[0,2]+self.k_joint*J,e[1,2]-self.k_joint*J,
                e[0,3]+ev*(b.lengths[1]*A-self.c*H)]
            spring=0.
            if self.kappa_s!=RIGID:
                spring=self.kappa_s*self.s**2
                natural.append(o[1,3]+spring*o[1,0])
            scale=max(1.,np.max(abs(e)),np.max(abs(o)),abs(ev*b.lengths[1]*A),abs(ev*self.c*H),abs(spring*o[1,0]),abs(self.k_joint*J))
            ni=di=0.;bending=[];mass=[]
            for arm,L in enumerate(b.lengths):
                values=[v@x for v in b.values[arm]]
                physical=values[0]-(self.c*A if arm else 0.)
                ni+=float(b.weights@(values[4]-ev*physical)**2)*L
                di+=float(b.weights@(values[4]**2+(ev*physical)**2))*L
                bending.append(float(b.weights@values[2]**2)*L)
                mass.append(float(b.weights@physical**2)*L)
            residual=self.Z.T@(self.K@x-ev*self.M@x)
            den=np.linalg.norm(self.Z.T@self.K@x)+ev*np.linalg.norm(self.Z.T@self.M@x)
            axial_mass=b.lengths[1]*self.s**2*A*A
            per.append(dict(natural=float(max(abs(np.array(natural)))/scale),natural_rows=np.array(natural).tolist(),strong=math.sqrt(ni/di),
                energy=abs(float(x@self.K@x)/ev-1),matrix_residual=float(np.linalg.norm(residual)/den),A=A,J=J,V0=float(o[1,0]),
                mass_F=mass[0],mass_W2=mass[1],mass_axial_translation=axial_mass,
                mass_square_identity_error=abs(mass[0]+mass[1]+axial_mass-float(x@self.M@x)),
                mixed_mass=-2*self.c*A*H,bending_F=bending[0],bending_V=bending[1],joint_energy=self.k_joint*J*J,
                spring_energy=spring*float(o[1,0])**2))
        d=dict(mass_gram=float(np.max(abs(coefficients.T@self.M@coefficients-np.eye(len(z))))),
            essential=float(np.max(abs(self.C@coefficients))),symmetry=float(max(np.max(abs(self.K-self.K.T)),np.max(abs(self.M-self.M.T)))))
        for key in ['natural','strong','energy','matrix_residual','mass_square_identity_error']:d[key]=max(r[key] for r in per)
        return d,per


def boundary_matrix(zeta,k,kappa_s,beta,mu):
    """Eight equations in the homogeneous F,W initial states; no Ritz data."""
    c,s=trig(beta);T=np.zeros((4,4));T[0,1]=T[1,2]=T[2,3]=1;T[3,0]=zeta
    o=[np.eye(8)[:4],np.eye(8)[4:]];e=[expm(T*L)@v for L,v in zip((1-mu,1+mu),o)]
    A=e[0][0];J=e[0][1]-e[1][1];d=o[1][0]+c*A
    rows=[o[0][0],o[0][1],e[1][0]+c*A,o[1][2],e[0][2]+k*J,e[1][2]-k*J,
        e[0][3]-c*e[1][3]+c*o[1][3]+zeta*(1+mu)*s*s*A,
        d if kappa_s==RIGID else o[1][3]+kappa_s*s*s*d]
    B=np.array(rows);return B/np.linalg.norm(B,axis=1)[:,None]


def endpoint_count(zeta,k,kappa_s,beta,mu,policy):
    """Model-specific extension of HS-LA5's EB endpoint Schur inertia.

Homogeneous W permits the SAME fixed-fixed beam pole catalogs. Add just the
translation mass L2*s^2*A^2 and impose W(L2)+c*A=0, plus one optional trace.
"""
    require(zeta>0,'Positive count probe');K=np.zeros((8,8));fixed=0;nearest=math.inf
    for arm,L in enumerate((1-mu,1+mu)):
        alpha=zeta**.25*L;n=1
        while (n+.5)*math.pi-.1<alpha+.2:
            pole=brentq(lambda x:math.cos(x)-1/math.cosh(x),(n+.5)*math.pi-.1,(n+.5)*math.pi+.1)
            nearest=min(nearest,abs(alpha-pole)/max(1,alpha));fixed+=int(pole<alpha);n+=1
        require(nearest>policy['pole_relative_guard'],'Auxiliary pole too close')
        T=np.zeros((4,4));T[0,1]=T[1,2]=T[2,3]=1;T[3,0]=zeta;T=expm(T*L)
        X=np.zeros((4,4));X[:2,:2]=np.eye(2);X[2:]=solve(T[:2,2:],np.column_stack((-T[:2,:2],np.eye(2))))
        Y=T@X;beam=np.vstack((X[3],-X[2],-Y[3],Y[2]));K[4*arm:4*arm+4,4*arm:4*arm+4]=(beam+beam.T)/2
    c,s=trig(beta);unit=np.eye(8);A=unit[2];J=unit[3]-unit[7];d=unit[4]+c*A
    K+=k*np.outer(J,J)-zeta*(1+mu)*s*s*np.outer(A,A)
    constraints=[unit[0],unit[1],unit[6]+c*A]
    if kappa_s==RIGID:constraints.append(d)
    else:K+=kappa_s*s*s*np.outer(d,d)
    Z=null_space(np.array(constraints));ev=np.linalg.eigvalsh(Z.T@K@Z)
    margin=float(min(abs(ev))/max(1,max(abs(ev))))
    require(margin>policy['inertia_relative_zero_tol'],'Endpoint inertia unresolved')
    return dict(count=int(fixed+sum(ev<0)),fixed_fixed_count=fixed,negative_inertia=int(sum(ev<0)),inertia_margin=margin,
        nearest_pole_distance=None if math.isinf(nearest) else nearest)


def independent(z,k,kappa_s,beta,mu,policy):
    require(np.linalg.matrix_rank(boundary_matrix(0,k,kappa_s,beta,mu))==8,'Unexpected zero mode')
    require(endpoint_count(1e-5,k,kappa_s,beta,mu,policy)['count']==0,'Missing low root')
    out=[]
    for i,value in enumerate(z):
        lo=value*(1-policy['count_bracket_relative']);hi=value*(1+policy['count_bracket_relative'])
        lower=endpoint_count(lo,k,kappa_s,beta,mu,policy);upper=endpoint_count(hi,k,kappa_s,beta,mu,policy)
        require(lower['count']==i and upper['count']==i+1,'First-three prefix count')
        root=brentq(lambda x:np.linalg.det(boundary_matrix(x,k,kappa_s,beta,mu)),lo,hi,xtol=1e-13,rtol=2e-14)
        err=abs(math.sqrt(root/value)-1);require(err<=policy['independent_frequency_tol'],'Independent ODE root')
        out.append(dict(position=i+1,zeta=float(root),frequency_error=err,lower=lower,upper=upper))
    return out


class Study:
    def __init__(self):
        self.cfg=read(CONFIG);self.p=self.cfg['policy'];self.source=ROOT/self.cfg['full_source']
        require(read(ROOT/self.cfg['analytical_gate'])['status']=='PASS','Analytical gate not passed')
        self.identity=dict(config=self.cfg,config_sha256=sha(CONFIG),initial_sha256=sha(ROOT/self.cfg['initial_workspace']),
            analytical_gate_sha256=sha(ROOT/self.cfg['analytical_gate']),versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
        self.run_id=hashlib.sha256(json.dumps(self.identity,sort_keys=True).encode()).hexdigest()[:16]
        self.folder=OUTPUT/self.run_id;self.folder.mkdir(parents=True,exist_ok=True)
        snapshot=self.folder/'workflow_versions'/(sha(__file__)+'.txt');snapshot.parent.mkdir(exist_ok=True)
        if not snapshot.exists():snapshot.write_bytes(Path(__file__).read_bytes())
        atomic(self.folder/'contract.json',self.identity)
        self.records=read(self.folder/'points.json') if (self.folder/'points.json').exists() else {}

    def point(self,k,ks,beta,mu,role):
        key=hashlib.sha256(repr((float(k),ks,float(beta),float(mu))).encode()).hexdigest()[:16]
        if key in self.records:return self.records[key]
        rows=[];solutions=[];models=[]
        try:
            for p in self.cfg['orders']:
                model=FiniteAngle(k,ks,beta,mu,p,self.cfg['quadrature']);sol=model.solve(3);d=sol['diagnostics']
                for name,tol in [('natural','boundary_tol'),('strong','strong_tol'),('essential','essential_tol'),('matrix_residual','matrix_tol'),
                                  ('symmetry','matrix_tol'),('mass_gram','mass_tol'),('energy','energy_tol'),('mass_square_identity_error','mass_tol')]:
                    require(d[name]<=self.p[tol],f'{k},{ks},{beta},{mu},p{p}: {name}={d[name]}')
                require(np.all(sol['zeta']>0) and np.all(np.diff(sol['zeta'])>0),'Positive ordered spectrum')
                checks=independent(sol['zeta'],k,ks,beta,mu,self.p) if p==24 else None
                rows.append(dict(order=p,zeta=sol['zeta'].tolist(),diagnostics=d,per_mode=sol['per_mode'],independent=checks))
                models.append(model);solutions.append(sol)
            error=float(max(abs(np.sqrt(solutions[0]['zeta']/solutions[1]['zeta'])-1)))
            require(error<=self.p['frequency_convergence_tol'],'p20/24 frequency convergence')
            # Compare the same sorted simple forms in the physical mass metric,
            # only between discretizations at one point, never across beta.
            c0=solutions[0]['coefficients'];c1=solutions[1]['coefficients'];m0,m1=models
            weights=np.r_[m1.base.weights*(1-mu),m1.base.weights*(1+mu),1.]
            fields=[]
            for m,coef in zip(models,[c0,c1]):
                A=m.A@coef;F=m.base.field_matrices(0,m1.base.nodes)[0]@coef
                W=m.base.field_matrices(1,m1.base.nodes)[0]@coef-m.c*A
                fields.append(np.vstack((F,W,math.sqrt(1+mu)*m.s*A)))
            signs=np.sign(np.sum(weights[:,None]*fields[0]*fields[1],axis=0))
            shape=float(max(np.sqrt(np.sum(weights[:,None]*(fields[0]-fields[1]*signs)**2,axis=0))))
            require(shape<=self.p['shape_convergence_tol'],'p20/24 mass-form convergence')
            result=dict(k_joint=k,kappa_s=ks,beta=beta,mu=mu,role=role,status='PASS',zero_semantics='0plus' if beta==0 and ks==RIGID else 'ordinary',
                zeta=rows[-1]['zeta'],orders=rows,frequency_convergence=error,shape_convergence=shape)
            np.savez_compressed(self.folder/(key+'.npz'),coefficients=solutions[1]['coefficients'])
            self.records[key]=result;atomic(self.folder/'points.json',self.records);return result
        except Exception as exc:
            atomic(self.folder/'failure.json',dict(status='UNRESOLVED',case=[k,ks,beta,mu],role=role,partial_orders=rows,error=repr(exc)));raise

    def audit(self):
        for path,digest in self.cfg['sources'].items():require(sha(ROOT/path)==digest,'Frozen source '+path)
        for folder,digest in self.cfg['archives'].items():
            path=ROOT/folder;require(sha(path/'manifest.json')==digest,'Source manifest')
            m=read(path/'manifest.json');require(m['status']=='PASS','Source not accepted')
            for name,h in m['artifacts'].items():require(sha(path/name)==h,'Source artifact '+name)
        atomic(self.folder/'source_audit.json',dict(status='PASS',new_full_solves=0,full_source=self.cfg['full_source']))

    def limits(self):
        path=self.folder/'limits.json'
        if path.exists():return
        rows=[]
        for mu in self.cfg['limits']['mu']:
            for k in self.cfg['limits']['k_joint']:
                finite=self.point(k,4.,0.,mu,'finite beta->0')
                saturated=self.point(k,RIGID,0.,mu,'rigid beta->0+')
                right=self.point(k,RIGID,90.,mu,'HS-LA5 exact limit')
                L=Leading(k,mu,24);O=RightAngle(k,mu,24)
                for label,new,old in [('finite0',finite,L.solve(0.)),('rigid0plus',saturated,L.solve(RIGID)),('rigid90',right,O.solve())]:
                    err=float(max(abs(np.sqrt(np.array(new['zeta'])/old['zeta'])-1)))
                    require(err<=self.p['independent_frequency_tol'],'Exact limit '+label)
                    rows.append(dict(mu=mu,k_joint=k,limit=label,frequency_error=err))
                r=FiniteAngle(k,RIGID,90,mu);require(np.array_equal(r.K,O.K) and np.array_equal(r.M,O.M) and np.array_equal(r.C,O.C),'HS-LA5 matrix identity')
        previous=None;ref=self.point(.1,RIGID,45.,.01,'exact spring limit')['zeta'];tail=[]
        for ks in self.cfg['limits']['stiffness_to_rigid']:
            row=self.point(.1,float(ks),45.,.01,'finite spring -> exact constraint');z=np.array(row['zeta'])
            require(np.all(z<=np.array(ref)*(1+1e-10)),'Penalty upper limit')
            if previous is not None:require(np.all(z>=previous*(1-1e-10)),'Penalty monotonicity')
            previous=z;tail.append(dict(kappa_s=ks,frequency_gap=(1-np.sqrt(z/ref)).tolist()))
        atomic(path,dict(status='PASS',exact_limits=rows,spring_to_rigid=tail,analytical_gate_sha256=sha(ROOT/self.cfg['analytical_gate'])))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['limits','run','verify','plot_only','finalize']);args=p.parse_args();s=Study()
    if args.action=='limits':s.audit();s.limits();print('LIMITS PASS',s.folder)
    elif args.action=='run':run(s)
    elif args.action=='verify':
        m=read(s.folder/'manifest.json')
        for name,digest in m['artifacts'].items():require(sha(s.folder/name)==digest,'Artifact '+name)
        print(m['status'],s.folder)
    elif args.action=='finalize':
        m=read(s.folder/'manifest.json')
        for path,digest in m['artifacts'].items():require(sha(s.folder/path)==digest,'Accepted numerical artifact changed')
        m['artifacts'].update({p.relative_to(s.folder).as_posix():sha(p) for p in (s.folder/'workflow_versions').glob('*.txt')})
        atomic(s.folder/'manifest.json',m);print('sealed',s.folder)
    else:render(s.folder)


def rows_csv(path):
    with Path(path).open(encoding='utf8',newline='') as f:return list(csv.DictReader(f))


def model_thresholds(study,case,ks,kind,eligible):
    cfg=study.cfg['threshold'];k=case['k_joint'];mu=case['mu']
    points=[study.point(k,ks,b,mu,'own-curve threshold classification') for b in study.cfg['angles_deg']]
    angles=[r['beta'] for r in points];freq=[math.sqrt(r['zeta'][0]) for r in points]
    own=freq[-1];outer=math.sqrt(study.point(k,RIGID,90.,mu,'independent HS-LA5 level')['zeta'][0])
    output=[]
    for ref,level in [('own90',own),('HS_LA5',outer)]:
        for q in cfg['q']:
            base=dict(id=case['id'],model=kind,q=q,reference=ref,omega_reference=case['epsilon']*level,
                      eligible_late_comparison=eligible,decreasing_sample_steps=sum(y<x for x,y in zip(freq[:-1],freq[1:])))
            if not eligible:
                output.append(dict(**base,status='SKIPPED_LATE_ACCURACY_GATE',beta_deg=None,crossings=[]));continue
            target=q*level;crossings=[]
            def fun(beta):return math.sqrt(study.point(k,ks,float(beta),mu,'apriori model threshold')['zeta'][0])/target-1
            for lo,hi,direction in crossing_brackets(angles,freq,target):
                root=float(brentq(fun,lo,hi,xtol=cfg['root_xtol_deg'],rtol=cfg['root_rtol']))
                half=cfg['bracket_width_deg']/2;left=max(0.,root-half);right=min(90.,root+half)
                fl,fr=fun(left),fun(right);res=fun(root)
                require(fl*fr<0 and abs(res)<=cfg['relative_equation_tol'],'Limit threshold signs/residual')
                crossings.append(dict(beta_deg=root,direction=direction,bracket=[left,right],signs=[fl,fr],residual=res))
            if not crossings and freq[0]>=target and all(w>=target for w in freq):
                status='ALREADY_HIGH_AT_ZERO_PLUS_NOT_FULL_THRESHOLD' if ks==RIGID else 'ALL_SAMPLED_ANGLES_HIGH'
                beta=None if ks==RIGID else 0.
            elif freq[-1]<target:status='NO_TERMINAL_HIGH_RANGE';beta=None
            elif len(crossings)==1 and crossings[0]['direction']=='up':status='PASS';beta=crossings[0]['beta_deg']
            elif crossings:status='MULTIPLE_CROSSINGS';beta=None
            else:status='UNRESOLVED';beta=None
            output.append(dict(**base,status=status,beta_deg=beta,crossings=crossings))
    return output


def run(s):
    if (s.folder/'manifest.json').exists():
        m=read(s.folder/'manifest.json')
        for path,digest in m['artifacts'].items():require(sha(s.folder/path)==digest,'Completed artifact '+path)
        print('REUSED',m['status'],s.folder,'no new spectral solutions');return
    s.audit();s.limits()
    full=rows_csv(s.source/'spectra.csv');comparisons=[];curves=[];eligibility=[];predictions=[]
    for case in s.cfg['cases']:
        family=[('III_rigid',RIGID)] if case['kappa_s']==RIGID else [('I_finite',case['kappa_s']),('II_rigid',RIGID)]
        low=Leading(case['k_joint'],case['mu'],24).solve(RIGID)['zeta'][0]
        for kind,ks in family:
            previous=None
            for beta in s.cfg['angles_deg']:
                r=s.point(case['k_joint'],ks,float(beta),case['mu'],'finite-angle curve')
                omega=case['epsilon']*np.sqrt(r['zeta']);pm=r['orders'][-1]['per_mode'][0]
                jump=None if previous is None else float(max(abs(omega[:2]-previous[:2])/np.maximum(omega[:2],previous[:2])))
                curves.append(dict(id=case['id'],model=kind,kappa_s=ks,beta=beta,omega1=float(omega[0]),omega2=float(omega[1]),omega3=float(omega[2]),
                    zero_semantics=r['zero_semantics'],inner_saturation_omega=case['epsilon']*math.sqrt(low),
                    neighbour_jump=jump,neighbour_trigger=jump is not None and jump>s.cfg['frequency_map']['neighbour_jump_relative'],
                    independent_at_every_point=True,**{x:pm[x] for x in ['mass_axial_translation','mass_W2','mass_F','mixed_mass','spring_energy']}))
                previous=omega
                if beta not in s.cfg['full_comparison_angles_deg']:continue
                f=next(row for row in full if row['id']==case['id'] and float(row['beta_deg'])==beta)
                require(float(f['epsilon'])==case['epsilon'] and float(f['mu'])==case['mu'] and float(f['aspect'])==case['aspect'],'Full geometry')
                expected=f['kappa_s'] if f['kappa_s']=='RIGID' else float(f['kappa_s'])
                require(expected==case['kappa_s'] and float(f['kappa_j'])==case['k_joint'],'Full spring normalization')
                record=read(ROOT/f['point_file']);require(record['status']=='PASS' and record['count']['status']=='PASS','Archive not accepted')
                w=float(f['omega1']);require(w==record['omega'][0],'CSV/point consistency')
                comparisons.append(dict(id=case['id'],model=kind,beta=beta,omega_limit=float(omega[0]),omega_full=w,
                    error_percent=100*(omega[0]/w-1),full_point=f['point_file'],full_sha256=sha(ROOT/f['point_file']),
                    inner_saturation_error_percent=100*(case['epsilon']*math.sqrt(low)/w-1)))
            late=[r for r in comparisons if r['id']==case['id'] and r['model']==kind and r['beta'] in s.cfg['threshold']['late_angles_deg']]
            err=float(max(abs(r['error_percent']) for r in late));ok=bool(err<=s.cfg['threshold']['late_frequency_acceptance_percent'])
            eligibility.append(dict(id=case['id'],model=kind,max_late_error_percent=err,eligible=ok))
            predictions+=model_thresholds(s,case,ks,kind,ok)
            print(case['id'],kind,'late_error',err,'eligible',ok,flush=True)
    csvwrite(s.folder/'curves.csv',curves);csvwrite(s.folder/'frequency_comparison.csv',comparisons)
    atomic(s.folder/'threshold_eligibility.json',eligibility)
    # Predictions are frozen BEFORE opening the historical threshold table.
    atomic(s.folder/'apriori_predictions.json',predictions)
    prediction_hash=sha(s.folder/'apriori_predictions.json')
    full_thresholds=read(s.source/'thresholds.json');old_predictions=rows_csv(s.source/'predictions.csv')
    threshold_comparison=[];frequency_at_full=[]
    for case in s.cfg['cases']:
        for q in s.cfg['threshold']['q']:
            fullrow=next(r for r in full_thresholds if r['id']==case['id'] and r['normalization']=='omega90' and r['q']==q)
            root=fullrow['terminal_beta_root_deg'];w=q*fullrow['omega_target']
            family=[('III_rigid',RIGID)] if case['kappa_s']==RIGID else [('I_finite',case['kappa_s']),('II_rigid',RIGID)]
            for kind,ks in family:
                result=s.point(case['k_joint'],ks,root,case['mu'],'at archived full threshold')
                wl=case['epsilon']*math.sqrt(result['zeta'][0])
                old=next(r for r in old_predictions if r['id']==case['id'] and r['normalization']=='omega90' and float(r['q'])==q)
                frequency_at_full.append(dict(id=case['id'],model=kind,q=q,full_beta_deg=root,omega_limit=wl,omega_full=w,
                    error_percent=100*(wl/w-1),HS_L1_error_percent=float(old['inner_frequency_error_at_full_root_percent'])))
                for prediction in [r for r in predictions if r['id']==case['id'] and r['model']==kind and r['q']==q]:
                    beta=prediction['beta_deg']
                    threshold_comparison.append(dict(id=case['id'],model=kind,q=q,reference=prediction['reference'],status=prediction['status'],
                        predicted_beta_deg=beta,full_beta_deg=root,error_deg=None if beta is None else beta-root,
                        error_percent=None if beta is None else 100*(beta/root-1)))
    require(sha(s.folder/'apriori_predictions.json')==prediction_hash,'Predictions changed after comparison')
    csvwrite(s.folder/'threshold_comparison.csv',threshold_comparison);csvwrite(s.folder/'frequency_at_full_thresholds.csv',frequency_at_full)
    atomic(s.folder/'prediction_provenance.json',dict(predictions_sha256=prediction_hash,full_thresholds_sha256=sha(s.source/'thresholds.json'),
        predictions_written_before_threshold_table_read=True,prior_HSB_threshold_knowledge=True,fit_parameters=0))
    initial=read(ROOT/s.cfg['initial_workspace'])
    for path,digest in initial['archives'].items():require(sha(ROOT/path)==digest,'Changed old archive '+path)
    diagnostics=[r['diagnostics'] for p in s.records.values() for r in p['orders']]
    maxres={key:max(d[key] for d in diagnostics) for key in ['natural','strong','essential','mass_gram','energy','matrix_residual','mass_square_identity_error']}
    independent_error=max(r['frequency_error'] for p in s.records.values() for r in p['orders'][-1]['independent'])
    m=dict(status='PASS',identity=s.identity,HEAD=initial['HEAD'],initial_git_status=initial['git_status'],new_full_solves=0,
        unique_limit_points=len(s.records),ritz_prefixes=2*len(s.records),independent_ode_roots=3*len(s.records),
        reused_full_common_points=len({r['full_point'] for r in comparisons}),reused_full_thresholds=18,protected_archive_files=len(initial['archives']),
        max_residuals=maxres,max_independent_frequency_error=independent_error,
        max_order_frequency_error=max(p['frequency_convergence'] for p in s.records.values()),
        max_order_shape_error=max(p['shape_convergence'] for p in s.records.values()),
        conditional_threshold_gate_failures=[r for r in eligibility if not r['eligible']],
        unresolved_predictions=[r for r in predictions if r['status']=='UNRESOLVED'])
    if m['unresolved_predictions']:m['status']='PARTIAL'
    m['artifacts']={p.relative_to(s.folder).as_posix():sha(p) for p in s.folder.rglob('*') if p.is_file() and p.name not in ['manifest.json','README.md','completion_checks.json','figures_manifest.json','documentation_manifest.json'] and 'figures' not in p.parts}
    atomic(s.folder/'manifest.json',m);print(m['status'],s.folder,len(s.records),'limit points',flush=True)


def render(folder):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder=Path(folder);curves=rows_csv(folder/'curves.csv');full=rows_csv(folder/'frequency_comparison.csv')
    figdir=folder/'figures';figdir.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':True,'grid.alpha':.25,'pdf.fonttype':42,
        'axes.spines.top':False,'axes.spines.right':False})
    labels={'I_finite':'I: конечная κₛ','II_rigid':'II: предел при χₛ>0','III_rigid':'III: точная связь'}
    colors={'I_finite':'#0072B2','II_rigid':'#D55E00','III_rigid':'#009E73'}
    for error,name in [(False,'01_frequency_comparison'),(True,'02_late_frequency_error')]:
        fig,axes=plt.subplots(2,2,figsize=(11,7.4))
        for ax,identifier in zip(axes.flat,['B0_A1','B0_A2','B0_A3','B0_A4']):
            kinds=list(dict.fromkeys(r['model'] for r in curves if r['id']==identifier))
            for kind in kinds:
                rr=sorted([r for r in (full if error else curves) if r['id']==identifier and r['model']==kind and (not error or float(r['beta'])>=30)],key=lambda r:float(r['beta']))
                ax.plot([float(r['beta']) for r in rr],[float(r['error_percent'] if error else r['omega1']) for r in rr],
                    '.-' if error else '-',color=colors[kind],label=labels[kind],lw=1.5)
            if not error:
                rr=[r for r in full if r['id']==identifier and r['model']==kinds[0]]
                ax.plot([float(r['beta']) for r in rr],[float(r['omega_full']) for r in rr],'o',color='black',ms=4,label='Полная модель: архив')
                sat=float(next(r['inner_saturation_omega'] for r in curves if r['id']==identifier))
                ax.axhline(sat,color='gray',ls=':',lw=1,label='Насыщение HS-L1')
            else:ax.axhline(0,color='black',lw=.7)
            ax.set(title=identifier.replace('B0_',''),xlim=(30,90) if error else (0,90),xlabel='Угол β, град',ylabel='Ошибка частоты, %' if error else 'Первая круговая частота ω₁')
            ax.legend(fontsize=8,loc='best')
        if not error:fig.suptitle('В точной связи левый предел означает β→0+, а не полную конфигурацию β=0',fontsize=10)
        fig.tight_layout();fig.savefig(figdir/(name+'.pdf'),bbox_inches='tight');fig.savefig(figdir/(name+'.png'),dpi=320,bbox_inches='tight');plt.close(fig)
    sources=['curves.csv','frequency_comparison.csv','threshold_comparison.csv','frequency_at_full_thresholds.csv']
    for name in sources:(figdir/name).write_bytes((folder/name).read_bytes())
    atomic(folder/'figures_manifest.json',dict(calculation_mode='plot_only',dpi=320,new_spectral_solves=0,
        sources={name:sha(folder/name) for name in sources},artifacts={p.name:sha(p) for p in figdir.iterdir() if p.is_file()}))
    print('plot_only',figdir)


if __name__=='__main__':main()
